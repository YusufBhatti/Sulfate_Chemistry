! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate quantities from diffusion and divergence damping (sect 13)
!
! Subroutine Interface:
SUBROUTINE Diagnostics_dif(                                       &
 row_length, rows, n_rows, model_levels, bl_levels,               &
! primary fields:
       theta,q,                                                         &
! wind field increments after diffusion:
       R_u,R_v,R_w,                                                     &
! Current theta+dtheta and q+dq values
       theta_star,q_star,                                               &
       exner_theta_levels,                                              &
       STASHwork)

USE atm_fields_bounds_mod, ONLY: tdims, tdims_l, tdims_s,               &
                                 udims, udims_s,  vdims,  vdims_s, wdims
USE turb_diff_mod, ONLY: L_subfilter_horiz, L_subfilter_vert
USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, shear, rneutml_sq
USE submodel_mod, ONLY: atmos_im
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE stash_array_mod, ONLY:                                                    &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE errormessagelength_mod, ONLY: errormessagelength

USE diff_increments_mod, ONLY:                                                &
    diff_inc_t, diff_inc_u, diff_inc_v, diff_inc_w, diff_inc_q, w_local_mask

IMPLICIT NONE
!
! Description:
!   Diagnostics_fildif calculates diagnostics for divergence damping
!   and diffusion for output by STASH routines for UM section 13.
!
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays from stash_array_mod.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (181,13) T increment          = delta(theta)/exner
!   (182,13) q increment          = delta(q)   across diffusion
!   (185,13) u increment          = delta(R_u) across diffusion
!   (186,13) v increment          = delta(R_v) across diffusion
!   (187,13) w increment          = delta(R_w) across diffusion
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from all routines.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: diffusion and filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable  ! Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: row_length, rows                           &
                                        ! horizontal dimensions
,n_rows                                                           &
                  ! rows for last (N) row of pes
,model_levels                                                     &
                  ! vertical levels
,bl_levels        ! vertical levels with moisture

!   Array  arguments with intent(in):
REAL, INTENT(IN) :: theta  (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) :: q      (tdims_l%i_start:tdims_l%i_end,        &
                            tdims_l%j_start:tdims_l%j_end,        &
                            tdims_l%k_start:tdims_l%k_end)
! Wind field increments for u, v and w 
REAL, INTENT(IN) :: R_u    (udims_s%i_start:udims_s%i_end,        &
                            udims_s%j_start:udims_s%j_end,        &
                            udims_s%k_start:udims_s%k_end)
REAL, INTENT(IN) :: R_v    (vdims_s%i_start:vdims_s%i_end,        &
                            vdims_s%j_start:vdims_s%j_end,        &
                            vdims_s%k_start:vdims_s%k_end)
REAL, INTENT(IN) :: R_w    (wdims  %i_start:wdims  %i_end,        &
                            wdims  %j_start:wdims  %j_end,        &
                            wdims  %k_start:wdims  %k_end)
REAL, INTENT(IN) :: theta_star                                    &
                           (tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                            tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: q_star (tdims%i_start:tdims%i_end,            &
                            tdims%j_start:tdims%j_end,            &
                            tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: exner_theta_levels                            &
                           (tdims_s%i_start:tdims_s%i_end,        &
                            tdims_s%j_start:tdims_s%j_end,        &
                            tdims_s%k_start:tdims_s%k_end)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
REAL, INTENT(OUT) ::                                              &
 STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
CHARACTER (LEN=*), PARAMETER :: RoutineName='DIAGNOSTICS_DIF'

INTEGER, PARAMETER :: sect =13  ! STASH section for diagnostics

! Local scalars:
INTEGER  ::                                                       &
 i,j,k                                                            &
             !  loop indices
,im_index                                                         &
             !  internal model index for STASH arrays
,item        !  STASH item of diagnostic
INTEGER :: Errorstatus = 0  ! initial value for error code

CHARACTER (LEN=errormessagelength) :: CMessage !  Error message

REAL, ALLOCATABLE :: work_visc(:,:,:),rneutml(:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header

!
! 1. Initialisation
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
im_index    = 1
Cmessage    = ''
!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! T increment
item = 181           ! T increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k= 1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        diff_inc_t(i,j,k) = (theta_star(i,j,k) - diff_inc_t(i,j,k)) &
                                * exner_theta_levels(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        diff_inc_t,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! q  increment
item = 182           ! q increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k= 1, tdims%k_end
    DO j= tdims%j_start, tdims%j_end
      DO i= tdims%i_start, tdims%i_end
        diff_inc_q(i,j,k) = q_star(i,j,k) - diff_inc_q(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        diff_inc_q,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! u wind increment
item = 185           ! u increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k= udims%k_start, udims%k_end
    DO j= udims%j_start, udims%j_end
      DO i= udims%i_start, udims%i_end
        diff_inc_u(i,j,k) = R_u(i,j,k) - diff_inc_u(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        diff_inc_u,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! v wind increment
item = 186           ! v increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k= vdims%k_start, vdims%k_end
    DO j= vdims%j_start, vdims%j_end
      DO i= vdims%i_start, vdims%i_end
        diff_inc_v(i,j,k) = R_v(i,j,k) - diff_inc_v(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        diff_inc_v,                                               &
        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! w wind increment
item = 187           ! w increment
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k= 1, wdims%k_end
    DO j= wdims%j_start, wdims%j_end
      DO i= wdims%i_start, wdims%i_end
        diff_inc_w(i,j,k) = R_w(i,j,k) - diff_inc_w(i,j,k)
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        diff_inc_w,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

ALLOCATE (work_visc(tdims%i_start:tdims%i_end,     &
                    tdims%j_start:tdims%j_end,     &
                                1:tdims%k_end))

IF (.NOT. L_subfilter_vert .AND. .NOT. L_subfilter_horiz) THEN
  sf(190,sect)=.FALSE.
  sf(191,sect)=.FALSE.
  sf(192,sect)=.FALSE.
  sf(193,sect)=.FALSE.
  sf(194,sect)=.FALSE.
  sf(195,sect)=.FALSE.
  sf(196,sect)=.FALSE.
  sf(197,sect)=.FALSE.
ELSE IF (.NOT. L_subfilter_vert) THEN
  sf(196,sect)=.FALSE.
  sf(197,sect)=.FALSE.
END IF

item = 190           ! momentum viscosity coeff
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_visc(i,j,k) = visc_m(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        work_visc,                                                &
        row_length,rows,model_levels,0,0,0,0,                     &
        at_extremity,                                             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF

item = 191           ! scalar viscosity coeff
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_visc(i,j,k) = visc_h(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        work_visc,                                                &
        row_length,rows,model_levels,0,0,0,0,                     &
        at_extremity,                                             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF
item = 192 ! shear
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_visc(i,j,k) = shear(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        work_visc,                                                &
        row_length,rows,model_levels,0,0,0,0,                     &
        at_extremity,                                             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF

item = 193 ! mixing length
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ALLOCATE (rneutml(row_length,rows,model_levels))
  DO k=1,model_levels-1
    DO j=1,rows
      DO i=1,row_length
        rneutml(i,j,k) = SQRT(rneutml_sq(i,j,k))
      END DO
    END DO
  END DO
  k=model_levels
  DO j=1,rows
    DO i=1,row_length
      rneutml(i,j,k) = 0.0
    END DO
  END DO
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
        rneutml,                                                  &
        row_length,rows,model_levels,0,0,0,0,                     &
        at_extremity,                                             &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF

DEALLOCATE (work_visc)

! Counter for occurances of local q diffusion
item = 201           ! local q diffusion at a point
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(item,sect,im_index)),                &
        w_local_mask,                                             &
        row_length,rows,0,0,0,0,at_extremity,                     &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)
!
! 3. Error handling
!
IF (Errorstatus /= 0) THEN

  CALL Ereport(RoutineName,Errorstatus,Cmessage)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Diagnostics_dif
