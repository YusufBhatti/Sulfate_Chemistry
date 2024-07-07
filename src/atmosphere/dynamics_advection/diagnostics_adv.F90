! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate diagnostic quantities from the initial atmosphere dump
!

MODULE diagnostics_adv_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_ADV_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE Diagnostics_adv(                                             &
       row_length,rows,n_rows,                                          &
! primary wind fields:
       u,v,theta,q,qcl,qcf,qrain,qgraup,qcf2,                           &
       mv,mcl,mcf,mrain,mgraup,mcf2,cf,cfl,cff,                         &
! wind field increments after advection:
       R_u,R_v,r_w,                                                     &
       theta_star,q_star,qcl_star,qcf_star,                             &
       qrain_star,qgraup_star,qcf2_star,                                &
       m_star, mcl_star, mcf_star,                                      &
       mrain_star, mgraup_star, mcf2_star,                              &
       cf_star,cfl_star,cff_star,                                       &
       exner_theta_levels,                                              &
! w departure point information
       depart_lambda,depart_phi,depart_r,r_theta_levels,                &
       STASHwork)

USE ac_diagnostics_mod, ONLY: qcl_adv
USE acp_namel_mod,      ONLY: l_ac
USE atm_fields_bounds_mod, ONLY : udims, udims_s, vdims, vdims_s,       &
                                  tdims, tdims_s, tdims_l,              &
                                  wdims
USE missing_data_mod,     ONLY: rmdi
USE nlsizes_namelist_mod, ONLY: model_levels

USE adv_increments_mod, ONLY:                                                 &
    adv_inc_t, adv_inc_u, adv_inc_v,                                          &
    adv_inc_q, adv_inc_qcl, adv_inc_qcf, adv_inc_qrain, adv_inc_qgraup,       &
    adv_inc_qcf2, adv_inc_m_v, adv_inc_m_cl, adv_inc_m_cf, adv_inc_m_r,       &
    adv_inc_m_gr, adv_inc_m_cf2, adv_inc_cf, adv_inc_cfl, adv_inc_cff

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                              &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
!   Diagnostics_adv extracts diagnostics of (N+1) time-level estimates
!   of primary fields after advection has been called, to be processed
!   by STASH routines for UM section 12 (advection).
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays from stash_array_mod
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (  2,12) u wind               = u + R_u
!   (  3,12) v wind               = v + R_v
!   (  4,12) temperature          = theta_star/exner_theta_levels
!   ( 10,12) specific humidity    = q_star
!   (254,12) qcl                  = qcl_star
!   ( 12,12) qcf                  = qcf_star
!   (185,12) u increment          = delta(R_u) across advection
!   (186,12) v increment          = delta(R_v) across advection
!   (187,12) w increment          = delta(w  ) across advection
!   (181,12) T increment          = delta(theta)/exner
!   (182,12) q   increment        = delta(q)   across advection
!   (183,12) qcl increment        = delta(qcl) across advection
!   (184,12) qcf increment        = delta(qcf) across advection
!   (189,12) qrain increment      = delta(qrain) across advection
!   (190,12) qgraup increment     = delta(qgraup) across advection
!   (191,12) qcf2 increment       = delta(qcf2) across advection
!   (192,12) cf  increment        = delta(cf)   across advection
!   (193,12) cfl increment        = delta(cfl) across advection
!   (194,12) cff increment        = delta(cff) across advection
!   (195,12) mv  increment        = delta(mv)  across advection
!   (196,12) mcl increment        = delta(mcl) across advection
!   (197,12) mcf increment        = delta(mcf) across advection
!   (198,12) mrain increment      = delta(mrain) across advection
!   (199,12) mgruap increment     = delta(mgraup) across advection
!   (200,12) mcf2 increment       = delta(mcf2) across advection
!   (204,12) departure point (w)  = depart_lambda
!   (203,12) departure point (w)  = depart_phi
!   (205,12) model height diff    = depart_r - r_theta_levels
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from physics1.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::                                                        &
 row_length,rows                                                  &
                  ! horizontal dimensions
,n_rows
                  ! rows for last (N) row of pes

!   Array  arguments with intent(in):
REAL :: u      (udims_s%i_start:udims_s%i_end,        &
                udims_s%j_start:udims_s%j_end,        &
                udims_s%k_start:udims_s%k_end)
REAL :: v      (vdims_s%i_start:vdims_s%i_end,        &
                vdims_s%j_start:vdims_s%j_end,        &
                vdims_s%k_start:vdims_s%k_end)

REAL :: theta  (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)

REAL :: q      (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qcl    (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qcf    (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)

REAL :: qrain  (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qgraup (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: qcf2   (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)

REAL :: mv     (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)
REAL :: mcl    (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)
REAL :: mcf    (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)

REAL :: mrain  (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)
REAL :: mgraup (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)
REAL :: mcf2   (tdims_s%i_start:tdims_s%i_end,        &
                tdims_s%j_start:tdims_s%j_end,        &
                tdims_s%k_start:tdims_s%k_end)

REAL :: cf     (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: cfl    (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)
REAL :: cff    (tdims_l%i_start:tdims_l%i_end,        &
                tdims_l%j_start:tdims_l%j_end,        &
                tdims_l%k_start:tdims_l%k_end)

REAL, TARGET :: R_u    (udims_s%i_start:udims_s%i_end,        &
                        udims_s%j_start:udims_s%j_end,        &
                        udims_s%k_start:udims_s%k_end)
REAL, TARGET :: R_v    (vdims_s%i_start:vdims_s%i_end,        &
                        vdims_s%j_start:vdims_s%j_end,        &
                        vdims_s%k_start:vdims_s%k_end)
REAL, TARGET :: R_w    (wdims%i_start:wdims%i_end,            &
                        wdims%j_start:wdims%j_end,            &
                        wdims%k_start:wdims%k_end)

REAL :: theta_star  (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: q_star      (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: qcl_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: qcf_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: qrain_star  (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: qgraup_star (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: qcf2_star   (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: m_star      (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: mcl_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: mcf_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: mrain_star  (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: mgraup_star (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: mcf2_star   (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: cf_star     (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: cfl_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)
REAL :: cff_star    (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                     tdims%k_start:tdims%k_end)

REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  &
                           tdims_s%j_start:tdims_s%j_end,  &
                           tdims_s%k_start:tdims_s%k_end)

REAL :: depart_lambda(row_length, rows, model_levels),     &
        depart_phi(row_length, rows, model_levels),        &
        depart_r(row_length, rows, model_levels)

REAL :: r_theta_levels(tdims_l%i_start:tdims_l%i_end,      &
                       tdims_l%j_start:tdims_l%j_end,      &
                                     0:tdims_l%k_end)

REAL, ALLOCATABLE :: qcl_incr_store(:,:,:),                &
                     qcf_incr_store(:,:,:),                &
                     cfl_incr_store(:,:,:),                &
                     cff_incr_store(:,:,:)

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
REAL :: STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='DIAGNOSTICS_ADV')

INTEGER :: sect      ! STASH section for diagnostics
PARAMETER ( sect = 12 )

! Local scalars:
INTEGER ::                                                        &
 i,j,k,ji                                                         &
             !  loop indices
,im_index                                                         &
             !  internal model index for STASH arrays
,item                                                             &
             !  STASH item of diagnostic
,Errorstatus !  Error status

CHARACTER(LEN=errormessagelength) :: CMessage !  Error message

! Local dynamic arrays:
REAL :: work_1(tdims%i_start:tdims%i_end,        &
               tdims%j_start:tdims%j_end,        &
                           1:tdims%k_end)
REAL :: work_u(udims%i_start:udims%i_end,        &
               udims%j_start:udims%j_end,        &
               udims%k_start:udims%k_end)
REAL :: work_v(vdims%i_start:vdims%i_end,        &
               vdims%j_start:vdims%j_end,        &
               vdims%k_start:vdims%k_end)
REAL :: adv_inc_w                                &
              (wdims%i_start:wdims%i_end,        &
               wdims%j_start:wdims%j_end,        &
                           1:wdims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header

!
! 1. Initialisation
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
im_index    = 1
Errorstatus = 0
Cmessage    = ''

! 1.5. If calculating separate positive and negative increments for cloud
! fields, need to store the values into temporary arrays.

IF ( (sf(170,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(171,sect) .AND. Errorstatus == 0)) THEN
  ALLOCATE ( qcl_incr_store(tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                                        1:tdims%k_end) )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,qcl_incr_store,adv_inc_qcl)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcl_incr_store(i,j,k) = adv_inc_qcl(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

END IF ! if requesting separate positive or negative incr for qcl

IF ( (sf(172,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(173,sect) .AND. Errorstatus == 0)) THEN
  ALLOCATE ( qcf_incr_store(tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                                        1:tdims%k_end) )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,qcf_incr_store,adv_inc_qcf)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcf_incr_store(i,j,k) = adv_inc_qcf(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

END IF ! if requesting separate positive or negative incr for qcf

IF ( (sf(176,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(177,sect) .AND. Errorstatus == 0)) THEN
  ALLOCATE ( cfl_incr_store(tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                                        1:tdims%k_end) )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,cfl_incr_store,adv_inc_cfl)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cfl_incr_store(i,j,k) = adv_inc_cfl(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

END IF ! if requesting separate positive or negative incr for cfl

IF ( (sf(178,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(179,sect) .AND. Errorstatus == 0)) THEN
  ALLOCATE ( cff_incr_store(tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                                        1:tdims%k_end) )

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,cff_incr_store,adv_inc_cff)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cff_incr_store(i,j,k) = adv_inc_cff(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO
END IF ! if requesting separate positive or negative incr for cfl

!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! u wind estimate = u + physics1 increment
item = 2
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,work_u,u,R_u)
  DO k= udims%k_start, udims%k_end
    DO j= udims%j_start, udims%j_end
      DO i= udims%i_start, udims%i_end
        work_u(i,j,k) = u(i,j,k) + R_u(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_u,      &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! v wind estimate = v + physics1 increment
item = 3
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,work_v,v,R_v)
  DO k= vdims%k_start, vdims%k_end
    DO j= vdims%j_start, vdims%j_end
      DO i= vdims%i_start, vdims%i_end
        work_v(i,j,k) = v(i,j,k) + R_v(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_v,      &
        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)


! temperature estimate = theta / exner pressure
item = 4
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_1,theta_star,exner_theta_levels)
  DO k=1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        work_1(i,j,k)=theta_star(i,j,k)*exner_theta_levels(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_1,      &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! specific humidity
item = 10
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        q_star(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end,                         &
                             1:tdims%k_end),                      &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcl
item = 254
IF ((sf(item,sect) .OR. l_ac) .AND. Errorstatus == 0) THEN

  IF (.NOT. ALLOCATED(qcl_adv)) THEN
    ALLOCATE ( qcl_adv(row_length*rows,model_levels) )
    qcl_adv(1,1) = rmdi
  END IF
END IF

IF (sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        qcl_star(tdims%i_start:tdims%i_end,                       &
                 tdims%j_start:tdims%j_end,                       &
                               1:tdims%k_end),                    &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  ! Copy non-halo area of qcl_star into qcl_adv in module ac_diagnostics_mod

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i,j,k,ji)                                                       &
!$OMP SHARED(tdims,row_length,qcl_adv,qcl_star)
  DO k=1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        ji = (j-1)*row_length+i
        qcl_adv(ji,k) = qcl_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! sf(item,sect)

! qcf
item = 12
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        qcf_star(tdims%i_start:tdims%i_end,                       &
                 tdims%j_start:tdims%j_end,                       &
                               1:tdims%k_end),                    &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! u wind increment
item = 185
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(udims,adv_inc_u,R_u)
  DO k= udims%k_start, udims%k_end
    DO j= udims%j_start, udims%j_end
      DO i= udims%i_start, udims%i_end
        adv_inc_u(i,j,k) = R_u(i,j,k) - adv_inc_u(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_u,                                                &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! v wind increment
item = 186
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,adv_inc_v,R_v)
  DO k= vdims%k_start, vdims%k_end
    DO j= vdims%j_start, vdims%j_end
      DO i= vdims%i_start, vdims%i_end
        adv_inc_v(i,j,k) = R_v(i,j,k) - adv_inc_v(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_v,                                                &
        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! w wind increment
item = 187
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(wdims,adv_inc_w,R_w)
  DO k=             1, wdims%k_end
    DO j= wdims%j_start, wdims%j_end
      DO i= wdims%i_start, wdims%i_end
        adv_inc_w(i,j,k) = R_w(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_w,                                                &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! T increment
! theta_star now holds theta+dtheta whereas t_incr holds dtheta before
! advection
item = 181
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i,k,j)                                                          &
!$OMP SHARED(tdims,adv_inc_t,theta_star,theta,exner_theta_levels)
  DO k= 1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      !CDIR NOUNROLL
      DO i= tdims%i_start, tdims%i_end
        adv_inc_t(i,j,k) = (theta_star(i,j,k) -           &
                 (theta(i,j,k) + adv_inc_t(i,j,k)))       &
                                 *exner_theta_levels(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_t,                                                &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! q increment
item = 182
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_q,q_star,q)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_q(i,j,k) = q_star(i,j,k) -                &
                     (q(i,j,k) + adv_inc_q(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_q,                                                &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcl increment
item = 183
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qcl,qcl_star,qcl)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_qcl(i,j,k) = qcl_star(i,j,k) -            &
                     (qcl(i,j,k) + adv_inc_qcl(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcl wind increment: positive
item = 170
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qcl,qcl_star,qcl,qcl_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_qcl(i,j,k) = MAX( 0.0, qcl_star(i,j,k) -  &
          (qcl(i,j,k) + qcl_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcl wind increment: negative
item = 171
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_qcl,qcl_star,qcl,qcl_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_qcl(i,j,k) = MIN( 0.0, qcl_star(i,j,k) -  &
          (qcl(i,j,k) + qcl_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf increment
item = 184           ! qcf increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf,qcf_star,qcf)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_qcf(i,j,k) = qcf_star(i,j,k) -            &
                     (qcf(i,j,k) + adv_inc_qcf(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcf,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf increment: positive
item = 172
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf,qcf_star,qcf,qcf_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_qcf(i,j,k) = MAX( 0.0, qcf_star(i,j,k) -  &
          (qcf(i,j,k) + qcf_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcf,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf increment: negative
item = 173
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf,qcf_star,qcf,qcf_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_qcf(i,j,k) = MIN( 0.0, qcf_star(i,j,k) -  &
          (qcf(i,j,k) + qcf_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcf,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qrain increment
item = 189
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qrain,qrain_star,qrain)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_qrain(i,j,k) = qrain_star(i,j,k) -        &
                     (qrain(i,j,k) + adv_inc_qrain(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qrain,                                            &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qgraup increment
item = 190
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qgraup,qgraup_star,qgraup)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_qgraup(i,j,k) = qgraup_star(i,j,k) -      &
                   (qgraup(i,j,k) + adv_inc_qgraup(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qgraup,                                           &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf2 increment
item = 191
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_qcf2,qcf2_star,qcf2)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_qcf2(i,j,k) = qcf2_star(i,j,k) -      &
                   (qcf2(i,j,k) + adv_inc_qcf2(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_qcf2,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_v increment
item = 195
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_v,m_star,mv)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_v(i,j,k) = m_star(i,j,k) -              &
                     (mv(i,j,k) + adv_inc_m_v(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_v,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cl increment
item = 196
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cl,mcl_star,mcl)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_cl(i,j,k) = mcl_star(i,j,k) -            &
                     (mcl(i,j,k) + adv_inc_m_cl(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_cl,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cf increment
item = 197
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cf,mcf_star,mcf)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_cf(i,j,k) = mcf_star(i,j,k) -            &
                     (mcf(i,j,k) + adv_inc_m_cf(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_cf,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_r increment
item = 198
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_r,mrain_star,mrain)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_r(i,j,k) = mrain_star(i,j,k) -        &
                     (mrain(i,j,k) + adv_inc_m_r(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_r,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_gr increment
item = 199
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_gr,mgraup_star,mgraup)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_gr(i,j,k) = mgraup_star(i,j,k) -      &
                     (mgraup(i,j,k) + adv_inc_m_gr(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_gr,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cf2 increment
item = 200
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_m_cf2,mcf2_star,mcf2)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_m_cf2(i,j,k) = mcf2_star(i,j,k) -          &
                     (mcf2(i,j,k) + adv_inc_m_cf2(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_m_cf2,                                            &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cf increment
item = 192           ! cf increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cf,cf_star,cf)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_cf(i,j,k) = cf_star(i,j,k) -              &
                     (cf(i,j,k) + adv_inc_cf(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cf,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cfl increment
item = 193           ! cfl increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cfl,cfl_star,cfl)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_cfl(i,j,k) = cfl_star(i,j,k) -            &
                     (cfl(i,j,k) + adv_inc_cfl(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cfl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cfl increment: positive
item = 176
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cfl,cfl_star,cfl,cfl_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_cfl(i,j,k) = MAX( 0.0, cfl_star(i,j,k) -  &
          (cfl(i,j,k) + cfl_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cfl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cfl increment: positive
item = 177
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,i,j)                                                          &
!$OMP SHARED(tdims,adv_inc_cfl,cfl_star,cfl,cfl_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_cfl(i,j,k) = MIN( 0.0, cfl_star(i,j,k) -  &
          (cfl(i,j,k) + cfl_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cfl,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cff increment
item = 194           ! cff increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cff,cff_star,cff)
  DO k= 1, tdims%k_end
    !CDIR NOUNROLL
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        adv_inc_cff(i,j,k) = cff_star(i,j,k) -            &
                     (cff(i,j,k) + adv_inc_cff(i,j,k))
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cff,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cff increment: positive
item = 178
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cff,cff_star,cff,cff_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_cff(i,j,k) = MAX( 0.0, cff_star(i,j,k) -   &
          (cff(i,j,k) + cff_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cff,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! cff increment: negative
item = 179
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,adv_inc_cff,cff_star,cff,cff_incr_store)
  DO k=1, tdims%k_end
    !CDIR NOUNROLL
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        adv_inc_cff(i,j,k) = MIN( 0.0, cff_star(i,j,k) -   &
          (cff(i,j,k) + cff_incr_store(i,j,k)) )
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        adv_inc_cff,                                              &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! Departure point diagnostics for w
! (a) lambda
item = 204
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        depart_lambda,                                            &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)
! (b) phi
item = 205
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        depart_phi,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)
! (c) dr  difference from model height of departure point
item = 203
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_1,depart_r,r_theta_levels)
  DO k=  1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        work_1(i,j,k) = depart_r(i,j,k) - r_theta_levels(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        work_1,                                                   &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! 2.5. Deallocate extra arrays used if separate positive and negative
!      increments were requested.

IF ( (sf(178,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(179,sect) .AND. Errorstatus == 0)) THEN
  DEALLOCATE( cff_incr_store )
END IF ! sf(item,sect)

IF ( (sf(176,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(177,sect) .AND. Errorstatus == 0)) THEN
  DEALLOCATE( cfl_incr_store )
END IF ! sf(item,sect)

IF ( (sf(172,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(173,sect) .AND. Errorstatus == 0)) THEN
  DEALLOCATE( qcf_incr_store )
END IF ! sf(item,sect)

IF ( (sf(170,sect) .AND. Errorstatus == 0) .OR.                   &
     (sf(171,sect) .AND. Errorstatus == 0)) THEN
  DEALLOCATE( qcl_incr_store )
END IF ! sf(item,sect)

! 3. Error handling
!
IF (Errorstatus /= 0) THEN
  CALL Ereport(RoutineName,Errorstatus,Cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_adv

END MODULE diagnostics_adv_mod
