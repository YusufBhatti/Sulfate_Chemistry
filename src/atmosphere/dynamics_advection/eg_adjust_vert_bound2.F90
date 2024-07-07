! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:

MODULE eg_adjust_vert_bound_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_ADJUST_VERT_BOUND_MOD'

CONTAINS
SUBROUTINE eg_adjust_vert_bound(                                      &
            row_length, rows, model_levels,g_i_pe,                    &
            pnt_type, dep_row_len,                                    &
            dep_rows, off_i, off_j, off_k, offz,                      &
            etadot, etadot_np1, depart_xi1, depart_xi2,               &
            depart_eta )

USE dynamics_input_mod,    ONLY: l_fast_vert_adjust,                 &
                                  l_sl_bc_correction
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE timestep_mod,          ONLY: timestep
USE um_parvars,            ONLY: halo_i, halo_j,                     &
                                  datastart,at_extremity,            &
                                  gc_proc_row_group,nproc_x,         &
                                  gc_all_proc_group
USE um_parcore,            ONLY: mype, nproc

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE ereport_mod,           ONLY: ereport
USE Field_Types
USE horiz_grid_mod
USE atm_fields_bounds_mod

USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE umPrintMgr,            ONLY:                                     &
    umPrint,                                                          &
    umMessage,                                                        &
    printstatus,                                                      &
    prstatus_diag

USE nlsizes_namelist_mod,  ONLY: global_row_length    
USE halo_exchange, ONLY: swap_bounds
USE um_types, ONLY: integer32
USE mpl, ONLY: mpl_integer4, mpl_max
USE model_domain_mod, ONLY: model_type, mt_global


IMPLICIT NONE
!
! Description:
!
!   Adjust departure points crossing the top and bottom full level using
!   the method described by eqs (10.51) and (10.52) of ENDGame formulation.
!
! Method: ENDGame formulation version 1.03, section 10.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_ADJUST_VERT_BOUND'



INTEGER, INTENT(IN) :: row_length, rows, model_levels, pnt_type

! row_length, rows and indexing offset for u or v or w type point

INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,           &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) :: g_i_pe(1-halo_i:global_row_length+halo_i)

! Logical switches for advection

REAL, INTENT(IN) ::          etadot(wdims_s%i_start:wdims_s%i_end,    &
                                    wdims_s%j_start:wdims_s%j_end,    &
                                    wdims_s%k_start:wdims_s%k_end)
REAL, INTENT(IN) ::      etadot_np1(wdims_s%i_start:wdims_s%i_end,    &
                                    wdims_s%j_start:wdims_s%j_end,    &
                                    wdims_s%k_start:wdims_s%k_end)
! Departure point coordinates

REAL, INTENT(OUT) ::                                                  &
     depart_xi1(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

REAL, INTENT(OUT) ::                                                  &
     depart_xi2(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

REAL, INTENT(INOUT) ::                                                &
     depart_eta(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

! Local variables

INTEGER :: i, j, k
REAL :: deta_1,                                                        &
        nlmax(model_levels), nlmin(model_levels), tau_lk,              &
        numax(model_levels), numin(model_levels), tau_uk

REAL ::                                                               &
  etadot_d(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j)

REAL, TARGET :: etadot_tb(1-halo_i:row_length+halo_i,                 &
                          1-halo_j:rows+halo_j, 2)  ! for top/bottom levs
REAL, POINTER :: etadot_1(:,:)            ! point into etadot_tb
REAL, POINTER :: etadot_nm1(:,:)          ! point into etadot_tb

REAL    :: eta_low, eta_upp

INTEGER :: ierr, Nk

INTEGER (KIND=integer32) :: check_in(model_levels)
INTEGER (KIND=integer32) :: check(model_levels)

! End of header


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( pnt_type == fld_type_w ) THEN
  Nk = model_levels-1
ELSE
  Nk = model_levels
END IF

IF ( l_fast_vert_adjust .OR.                                         &
    (l_sl_bc_correction .AND. pnt_type /= fld_type_w) ) THEN
  SELECT CASE ( pnt_type )
  CASE ( fld_type_w )
    eta_low  = eta_theta_levels(0)
    eta_upp  = eta_theta_levels(model_levels)
  CASE ( fld_type_u, fld_type_v, fld_type_p )
    eta_low  = eta_rho_levels(1)
    eta_upp  = eta_rho_levels(model_levels)
  CASE DEFAULT
    ierr = 1
    CALL ereport("eg_adjust_vert_bound", ierr, "Invalid grid point type" )
  END SELECT

  DO k=1, Nk
    DO j=1-off_j, dep_rows-off_j
      DO i=1-off_i, dep_row_len-off_i
        depart_eta(i,j,k) = MAX(depart_eta(i,j,k), eta_low)
        depart_eta(i,j,k) = MIN(depart_eta(i,j,k), eta_upp)
      END DO
    END DO
  END DO

  GO TO 9999
END IF

! Copy etadot at lev 1 and ml-1 at single level extended arrays for use
! by eg_bi_linear_h()

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( pdims, etadot_tb, etadot, model_levels )                &
!$OMP& PRIVATE( i, j )
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    etadot_tb(i,j,1) = etadot(i,j,1)
    etadot_tb(i,j,2) = etadot(i,j,model_levels-1)
  END DO
END DO
!$OMP END PARALLEL DO

CALL swap_bounds( etadot_tb,row_length,rows,2,halo_i,               &
                                       halo_j,fld_type_p, swap_field_is_scalar)

etadot_1   => etadot_tb(:,:,1)
etadot_nm1 => etadot_tb(:,:,2)

! Adjust lower boundary.

deta_1  = eta_theta_levels(1) - eta_theta_levels(0)

check = 0
check_in = 0

DO k = 1, Nk
  ! bottom
  SELECT CASE ( pnt_type )

  CASE ( fld_type_w )
    nlmax(k) = eta_theta_levels(k)
    nlmin(k) = eta_theta_levels(1)
  CASE ( fld_type_u, fld_type_v, fld_type_p )
    nlmax(k) = MAX ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
    nlmin(k) = MIN ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
  CASE DEFAULT
    ierr = 1
    CALL ereport("eg_adjust_vert_bound", ierr,                         &
                 "Invalid grid point type" )
  END SELECT

  ! top
  SELECT CASE ( pnt_type )

  CASE ( fld_type_w )
    numax(k) = eta_theta_levels(model_levels-1)
    numin(k) = eta_theta_levels(k)
  CASE ( fld_type_u, fld_type_v, fld_type_p )
    numax(k) = MAX ( eta_rho_levels(k)  ,                          &
                  eta_theta_levels(model_levels-1) )
    numin(k) = MIN ( eta_rho_levels(k)  ,                          &
                  eta_theta_levels(model_levels-1) )
  CASE DEFAULT
    ierr = 1
    CALL ereport("eg_adjust_vert_bound", ierr,                        &
             "Invalid grid point type" )

  END SELECT

END DO

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( Nk, depart_eta, off_i,  dep_row_len, off_j,             &
!$OMP&         dep_rows, nlmin, numax, check_in )                      &
!$OMP& PRIVATE( k )
DO k = 1, Nk

  IF ( ANY(depart_eta(1-off_i: dep_row_len-off_i,                      &
                      1-off_j: dep_rows-off_j,k) < nlmin(k)) )         &
       check_in(k) = 1

  IF ( ANY(depart_eta(1-off_i: dep_row_len-off_i,                      &
                      1-off_j: dep_rows-off_j,k) > numax(k)) )         &
       check_in(k) = 2

END DO
!$OMP END PARALLEL DO

! If we want diagnostics, we'll need to do communications across all 
! processes. Otherwise, for global models just the row-group is sufficient,
! as communications within eg_bi_linear_h are only row-wise.
! Needless to say a reduced communication scope is cheaper!
! For LAMs no communication is required as all adjustments are PE local
! due to the lack of comms-on-demand.
! Note the printstatus here needs to match that in locate_hdps to ensure
! comms don't deadlock
IF (printstatus == prstatus_diag) THEN
  CALL mpl_allreduce(check_in, check, nk, mpl_integer4, mpl_max,         &
                     gc_all_proc_group, ierr)
ELSE
  IF (model_type == mt_global) THEN
    CALL mpl_allreduce(check_in, check, nk, mpl_integer4, mpl_max,         &
                       gc_proc_row_group, ierr)
  ELSE
    check(:) = check_in(:)
  END IF
END IF

IF ( mype == 0 .AND. printstatus == prstatus_diag) THEN
  WRITE(umMessage,*)  "***********************************************"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
  WRITE(umMessage,*)  "* vertical adjustment performed at levels     *"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
  WRITE(umMessage,*)  "* =========================================== *"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
  WRITE(umMessage,*)  "* level:  where: 1=adjust bottom 2=adjust top *"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
  DO k=1, model_levels-1
    IF (check(k)/=0) THEN
      WRITE(umMessage,'(A,2I4,A)') ' *  ',k,&
          check(k), "                                   *"
      CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
    END IF
  END DO
  WRITE(umMessage,*)  "* if no output above then no adjustments  *"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
  WRITE(umMessage,*)  "*******************************************"
  CALL umPrint(umMessage,src='eg_adjust_vert_bound2')
END IF

DO k=1, Nk

  ne0:IF (check(k)/=0) THEN

    chkk1:IF (check(k)==1) THEN

      deta_1  = eta_theta_levels(1) - eta_theta_levels(0)

      ! DEPENDS ON: eg_bi_linear_h
      CALL eg_bi_linear_h( etadot_1, depart_xi1(:,:,k:k),             &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1,                                 &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity,                  &
            gc_proc_row_group, etadot_d )

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( off_i, off_j, dep_row_len, dep_rows, depart_eta, nlmin, &
!$OMP&         nlmax, eta_theta_levels, timestep, deta_1, etadot_np1,  &
!$OMP&         etadot_d, k, Nk )                                       &
!$OMP& PRIVATE( i, j, tau_lk )
      DO j=1-off_j, dep_rows-off_j
        DO i=1-off_i, dep_row_len-off_i

          IF ( depart_eta(i,j,k) < nlmin(k) ) THEN

            tau_lk = (nlmax(k) -  eta_theta_levels(1)) /               &
              (nlmax(k)-MAX(eta_theta_levels(0),depart_eta(i,j,k) ))

            depart_eta(i,j,k) = eta_theta_levels(0) +                  &
                            (nlmin(k)-eta_theta_levels(0)) *           &
                            EXP(-timestep*(1.0 - tau_lk)/( deta_1 )    &
                            * ((1.0 - tau_lk)/2.0*etadot_np1(i,j,1) +  &
                            (1.0 + tau_lk)/2.0*etadot_d(i,j)))

            depart_eta(i,j,k) = MIN( depart_eta(i,j,k) ,               &
                                 eta_theta_levels(1) )

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF chkk1


    chkk2:IF ( check(k)==2 ) THEN

        ! Adjust upper boundary.

      deta_1 = eta_theta_levels(model_levels) -                     &
               eta_theta_levels(model_levels-1)

      ! DEPENDS ON: eg_bi_linear_h
      CALL eg_bi_linear_h( etadot_nm1, depart_xi1(:,:,k:k),           &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1,                                 &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity,                  &
            gc_proc_row_group, etadot_d )

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                      &
!$OMP& SHARED( off_i, off_j, dep_row_len, dep_rows,  depart_eta,       &
!$OMP&         numax, eta_theta_levels, numin, model_levels, timestep, &
!$OMP&         deta_1, etadot_np1, etadot_d, k, Nk )                   &
!$OMP& PRIVATE( i, j, tau_uk )
      DO j=1-off_j, dep_rows-off_j
        DO i=1-off_i, dep_row_len-off_i

          IF ( depart_eta(i,j,k) >  numax(k) ) THEN

            tau_uk = (eta_theta_levels(model_levels-1) - numin(k)) /   &
                 (MIN(eta_theta_levels(model_levels),                  &
                  depart_eta(i,j,k) )-numin(k))

            depart_eta(i,j,k) = eta_theta_levels(model_levels) -       &
              ( eta_theta_levels(model_levels) - numax(k) ) *          &
              EXP( timestep * (1.0-tau_uk)/( deta_1 ) *                &
              ((1.0-tau_uk)/2.0*etadot_np1(i,j,model_levels-1) +       &
                             (1.0+tau_uk)/2.0*etadot_d(i,j)))

            depart_eta(i,j,k) = MAX( depart_eta(i,j,k) ,               &
                                 eta_theta_levels(model_levels-1) )

          END IF
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF chkk2
  END IF ne0

END DO

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_adjust_vert_bound
END MODULE
