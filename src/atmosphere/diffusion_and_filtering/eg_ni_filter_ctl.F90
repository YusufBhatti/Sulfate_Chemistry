! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_filter_Ctl
!
! Purpose:
!          Control polar filtering and diffusion
!
! Method:
!          Using 1-2-1 Shapiro filter to zap 2-grid waves
!          T. Davies.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
MODULE eg_NI_filter_Ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_NI_FILTER_CTL_MOD'

CONTAINS
SUBROUTINE eg_NI_filter_Ctl(                                      &
                         thetav, u, v, w, etadot, Exner,          &
                         exner_theta_levels,                      &
                         row_length, rows, n_rows, model_levels,  &
                         r_theta_levels, r_rho_levels,            &
                         r_at_u, r_at_v,                          &
                         max_filter_rows, u_sweeps, v_sweeps,     &
                         global_u_filter, global_v_filter,        &
                         u_begin, u_end, v_begin, v_end,          &
                         diff_coeff_phi,                          &
                         diff_coeff_u, diff_coeff_v,              &
                         first_constant_r_rho_level,              &
                         first_constant_r_rho_level_m1,           &
                         top_filt_start, top_filt_end,            &
                         up_diff, max_upd_levels,                 &
                         horizontal_level,                        &
                         off_x, off_y, halo_i, halo_j,            &
                         n_procy, at_extremity,                   &
                         L_pftheta, L_pfuv,                       &
                         L_pfw, L_pfexner, L_diff_exner,          &
                         L_diff_thermo, L_diff_wind, L_diff_w,    &
                         L_pofil_hadgem2,                         &
                         xi1_u, xi1_p, xi2_p, xi2_v,              &
                         pole_consts, gc_proc_row_group,          &
                         nproc, gc_proc_col_group,                &
                         global_row_length,                       &
                         csxi2_v, csxi2_p,                        &
                         delta_lambda, delta_phi,                 &
                         STASHwork )

USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE mpp_conf_mod,     ONLY: swap_field_is_vector
USE trignometric_mod, ONLY:  sin_theta_longitude,                &
                              cos_theta_longitude,                &
                         cos_theta_latitude, sec_theta_latitude,  &
                         cos_v_latitude, sec_v_latitude

USE swapable_field_mod, ONLY:                                    &
    swapable_field_pointer_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE eg_v_at_poles_mod
USE atm_fields_bounds_mod
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels
USE eg_pofil_vatp_mod
USE UM_ParParams
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE missing_data_mod, ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength
USE mpp_conf_mod,     ONLY: swap_field_is_scalar, swap_field_is_vector
USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.
LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

LOGICAL ::                                                        &
  L_pftheta                                                       &
                 ! switch for polar filter for theta
, L_pfw                                                           &
                 ! switch for polar filter for w
, L_pfuv                                                          &
                 ! switch for polar filter for horizontal winds
, L_pfexner                                                       &
                 ! switch for polar filter for Exner pressure
, L_diff_Exner                                                    &
                 ! switch for diffusion of Exner pressure
, L_diff_thermo                                                   &
                 ! switch for horiz. diffusion of theta
, L_diff_wind                                                     &
                 ! switch for horiz. diffusion of u,v
, L_diff_w       ! switch for horiz. diffusion of w

LOGICAL :: L_pofil_hadgem2  ! use hadgem2 polar filtering settings

INTEGER ::                                                        &
  max_filter_rows                                                 &
                     ! max array size for u_begin etc
, max_upd_levels                                                  &
                   ! max no. levels for upper diffusion
, global_u_filter                                                 &
                   ! number of filter sweeps; 0=diffusion only
, global_v_filter                                                 &
                   ! number of filter sweeps; 0=diffusion only
, row_length                                                      &
                   ! number of points on a row.
, rows                                                            &
                   ! number of rows.
, n_rows                                                          &
                   ! number of rows in a v field
, model_levels                                                    &
                   ! number of model levels.
, first_constant_r_rho_level                                      &
, first_constant_r_rho_level_m1                                   &
, off_x                                                           &
             ! Size of small halo in i
, off_y                                                           &
             ! Size of small halo in j.
, halo_i                                                          &
             ! Size of halo in i direction.
, halo_j
             ! Size of halo in j direction.


REAL ::                                                           &
     ! vertical co-ordinate arrays.
  r_theta_levels (1-halo_i:row_length+halo_i,                     &
                  1-halo_j:rows+halo_j, 0:model_levels)           &
, r_rho_levels (1-halo_i:row_length+halo_i,                       &
                1-halo_j:rows+halo_j, model_levels)               &
, r_at_u (1-halo_i:row_length+halo_i,                             &
          1-halo_j:rows+halo_j, model_levels)                     &
, r_at_v (1-halo_i:row_length+halo_i,                             &
          1-halo_j:n_rows+halo_j, model_levels)

REAL  ::  up_diff(max_upd_levels)
                            !upper-level diffusion coefficients

REAL  ::                                                          &
  diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y ),    &
                      !    diffusion coefficient for u/theta rows
  diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y ),  &
                      !    diffusion coefficient for v rows
  diff_coeff_phi_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y), &
                      !    NS diffusion coefficient for u/theta rows
  diff_coeff_phi_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)
                      !    NS diffusion coefficient for v rows

REAL  :: diff_coeff_phi   ! NS diffusion coeff

REAL :: delta_lambda, delta_phi

INTEGER ::                                                        &
  n_procy                                                         &
, horizontal_level                                                &
                               ! steep slope control
, top_filt_start                                                  &
                       ! start level for upper-level diffusion
, top_filt_end         ! END level for upper-level diffusion

INTEGER ::                                                        &
  u_sweeps(max_filter_rows)                                       &
                             ! sweeps for 1-2-1 filter
, v_sweeps(max_filter_rows)                                       &
                             ! sweeps for 1-2-1 filter
, u_begin(0:max_filter_rows)                                      &
                              ! row pointers for 1-2-1 filter
, u_end(0:max_filter_rows)                                        &
                              ! row pointers for 1-2-1 filter
, v_begin(0:max_filter_rows)                                      &
                              ! row pointers for 1-2-1 filter
, v_end(0:max_filter_rows)    ! row pointers for 1-2-1 filter

! Arguments with Intent IN/OUT.

REAL, TARGET ::                                                   &
  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     model_levels)                                                &
, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
      model_levels)                                               &
, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
      0:model_levels)                                             &
, etadot(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
      0:model_levels)                                             &
, thetav(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
      0:model_levels)

REAL ::                                                           &
  exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
         model_levels + 1)                                        &
, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
         0:model_levels )

! Varibles used for eg_v_at_poles
REAL :: pole_consts(4)

! grid locations
REAL :: xi1_p(1-halo_i:row_length+halo_i),                          &
      xi1_u(-halo_i:row_length-1+halo_i),                         &
      xi2_p(1-halo_j:rows+halo_j),                                &
      xi2_v(-halo_j:n_rows-1+halo_j)

INTEGER, INTENT(IN) :: gc_proc_row_group, global_row_length

REAL ::                                                            &
      csxi2_p(1-halo_j:rows+halo_j),                               &
      csxi2_v(-halo_j:n_rows-1+halo_j)

!   Array arguments with intent(out):
REAL, INTENT(OUT) ::                                              &
 STASHwork(*)   ! Output array holding diagnostic fields


! Local variables
LOGICAL :: l_u_diff  ! TRUE if diffusing u.

INTEGER ::                                                        &
  i, j, k                                                         &
                    ! loop counters
, j_start, j_stop                                                 &
                    ! Loop bounds
, active_levels  ! number of levels for upper-level diffusion

INTEGER :: i_field  ! counter for swapbounds

REAL, ALLOCATABLE ::                                              &
  r_theta_at_u(:,:,:),                                            &
  r_theta_at_v(:,:,:),                                            &
  r_uv_b(:,:,:)

REAL :: u_inc (row_length, rows, model_levels)
REAL :: v_inc (row_length, n_rows, model_levels)
REAL :: w_inc (row_length, rows, 0:model_levels)
REAL :: T_inc (row_length, rows, model_levels)
REAL :: exner_inc (row_length, rows, model_levels+1)
INTEGER :: sect
INTEGER :: item
INTEGER :: im_index      !  internal model index for STASH arrays
INTEGER :: Errorstatus = 0  ! initial value for error code
CHARACTER (LEN=errormessagelength) :: CMessage !  Error message

TYPE(swapable_field_pointer_type) :: fields_to_swap(4)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_NI_FILTER_CTL'

! z level test arrays and diagnostic
!       REAL ::                                                           &
!         diff_temp_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
!           model_levels),                                                &
!         diff_temp_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,     &
!           model_levels),                                                &
!         diff_temp_p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
!             model_levels),                                              &
!         diff_temp_w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
!             0:model_levels)
!       REAL :: val, flux, flux_out
!
!       LOGICAL :: L_eta_level_diff = .true.

! Reproducable sum varibles
INTEGER, INTENT(IN)    :: nproc
INTEGER, INTENT(IN)    :: gc_proc_col_group

INTEGER :: istat
REAL    :: temp
REAL    :: k_sum(row_length,rows)
REAL    :: row_sum(rows)

! ----------------------------------------------------------------------
! 0.1  Section for setting up diagnostic arrays
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
sect=13
im_index    = 1
Cmessage    = ''

IF (sf(0,sect)) THEN
  IF (sf(381,sect))t_inc(:,:,:)= thetav(1:row_length, 1:rows,  &
                                       1:model_levels)
  IF (sf(385,sect))u_inc(:,:,:)=u(1:row_length, 1:rows, :)
  IF (sf(386,sect))v_inc(:,:,:)=v(1:row_length, 1:n_rows, :)
  IF (sf(387,sect))w_inc(:,:,:)=w(1:row_length, 1:rows, :)
  IF (sf(388,sect))exner_inc(:,:,:)= exner(1:row_length, 1:rows, :)
END IF

! ----------------------------------------------------------------------
! 1.0  Section for polar filter
! ----------------------------------------------------------------------

ALLOCATE ( r_theta_at_u(1-off_x:row_length+off_x,                 &
                        1-off_y:rows+off_y, 0:model_levels) )
ALLOCATE ( r_theta_at_v(1-off_x:row_length+off_x,                 &
                        1-off_y:n_rows+off_y, 0:model_levels) )

DO k = 0, model_levels
  DO j = 1-off_y, rows+off_y
    DO i = 1-off_x, row_length+off_x
      ! END loop bound OK since r_theta_levels is halo_i > off_x
      r_theta_at_u(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
                                  r_theta_levels(i  ,j,k) )
    END DO
  END DO
END DO ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

DO k = 0, model_levels
  DO j = 1-off_y, n_rows+off_y
    ! END loop bound OK since r_theta_levels is halo_j > off_y
    DO i = 1-off_x, row_length+off_x
      r_theta_at_v(i,j,k) = .5 * ( r_theta_levels(i,j+1,k) +      &
                                   r_theta_levels(i,j  ,k) )
    END DO
  END DO
END DO ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

! Copy NS diffcoeff into NS diffcoeff arrays
DO j = 1-off_y,rows+off_y
  DO i = 1-off_x,row_length+off_x
    diff_coeff_phi_u(i,j) = diff_coeff_phi
  END DO
END DO
DO j = 1-off_y,n_rows+off_y
  DO i = 1-off_x,row_length+off_x
    diff_coeff_phi_v(i,j) = diff_coeff_phi
  END DO
END DO

! ----------------------------------------------------------------------
! 2.1  polar filter/diffusion theta/w  polar filter Exner
! ----------------------------------------------------------------------

! Swap_bounds are Done inside the polar filter sweeps
IF ( L_pftheta .OR. L_diff_thermo ) THEN

  l_u_diff = .FALSE.
  CALL eg_pofil_vatp(                                             &
                  thetav(1-off_x,1-off_y, 1),                     &
                  fld_type_p, 0, 0, 0, 1, 0, 0,                   &
                  model_levels, model_levels, model_levels,       &
                  first_constant_r_rho_level_m1,                  &
                  rows, n_rows, rows, row_length,                 &
                  r_theta_levels, r_rho_levels,                   &
                  r_theta_at_u, r_theta_at_v,                     &
                  off_x, off_y, off_x, off_y,                     &
                  halo_i, halo_j, halo_i, halo_j,                 &
                  off_x, off_y, off_x, off_y,                     &
                  sec_theta_latitude, cos_v_latitude,             &
                  n_procy, max_filter_rows, global_u_filter,      &
                  u_sweeps, u_begin, u_end,                       &
                  horizontal_level, diff_coeff_phi_u,             &
                  diff_coeff_u, L_diff_thermo,swap_field_is_scalar,&
                  L_pofil_hadgem2, csxi2_p, csxi2_v, 1,l_u_diff )

  DO i = 1-off_x,row_length+off_x
    DO j =  1-off_y,rows+off_y
      thetav(i,j,0) = thetav(i,j,1)
    END DO
  END DO
END IF !  L_pftheta .or. L_diff_thermo

IF ( L_pfw .OR. L_diff_w ) THEN
  l_u_diff = .FALSE.
  CALL eg_pofil_vatp(                                             &
                  etadot(1-off_x,1-off_y, 1), fld_type_p,         &
                  0, 0, 0, 1, 0, 0,                               &
                  model_levels, model_levels, model_levels-1,     &
                  first_constant_r_rho_level_m1,                  &
                  rows, n_rows, rows, row_length,                 &
                  r_theta_levels, r_rho_levels,                   &
                  r_theta_at_u, r_theta_at_v,                     &
                  off_x, off_y, off_x, off_y,                     &
                  halo_i, halo_j, halo_i, halo_j,                 &
                  off_x, off_y, off_x, off_y,                     &
                  sec_theta_latitude, cos_v_latitude,             &
                  n_procy, max_filter_rows, global_u_filter,      &
                  u_sweeps, u_begin, u_end,                       &
                  horizontal_level, diff_coeff_phi_u,             &
                  diff_coeff_u, L_diff_w, swap_field_is_scalar,   &
                  L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

END IF ! L_pfw .or. L_diff_w

IF ( L_pfexner ) THEN
  l_u_diff = .FALSE.
  CALL eg_pofil_vatp(                                             &
                  Exner, fld_type_p, 0, 0, 1, 0, 1, 1,            &
                  model_levels+1, model_levels, model_levels+1,   &
                  first_constant_r_rho_level,                     &
                  rows, n_rows, rows, row_length,                 &
                  r_rho_levels, r_theta_levels,                   &
                  r_at_u, r_at_v,                                 &
                  off_x, off_y, off_x, off_y,                     &
                  halo_i, halo_j, halo_i, halo_j,                 &
                  halo_i, halo_j, halo_i, halo_j,                 &
                  sec_theta_latitude, cos_v_latitude,             &
                  n_procy, max_filter_rows, global_u_filter,      &
                  u_sweeps, u_begin, u_end,                       &
                  horizontal_level, diff_coeff_phi_u,             &
                  diff_coeff_u, L_diff_Exner,swap_field_is_scalar,&
                  L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)
END IF !  L_pfexner

IF ( L_pfuv .OR. L_diff_wind ) THEN

  ALLOCATE ( r_uv_b(1-off_x:row_length+off_x,                     &
                    1-off_y:n_rows+off_y, model_levels) )

  !  r needed at centre of grid but no swap bound option for grid-type
  !      so fill required halos explicitly
  !      (fortunately have large halos to start)
  DO k = 1, model_levels
    DO j = 1-off_y, n_rows+off_y
      ! END loop bound OK since r_rho_levels is halo_j > off_y
      DO i = 1-off_x, row_length+off_x
        r_uv_b(i,j,k) = .5 * ( r_at_u(i,j,k) + r_at_u(i,j+1,k) )
      END DO
    END DO
  END DO ! k = 1, model_levels

  ! ----------------------------------------------------------------------
  ! 2.2  polar filter/diffusion u/v
  ! ----------------------------------------------------------------------
  l_u_diff = .TRUE.
  CALL eg_pofil_vatp(                                             &
                  u, fld_type_u, 1, 0, 1, 0, 1, 1,                &
                  model_levels, model_levels, model_levels,       &
                  first_constant_r_rho_level,                     &
                  rows, n_rows, rows, row_length,                 &
                  r_at_u, r_theta_at_u, r_rho_levels, r_uv_b,     &
                  off_x, off_y, off_x, off_y,                     &
                  halo_i, halo_j, off_x, off_y,                   &
                  halo_i, halo_j, off_x, off_y,                   &
                  sec_theta_latitude, cos_v_latitude,             &
                  n_procy, max_filter_rows, global_u_filter,      &
                  u_sweeps, u_begin, u_end,                       &
                  horizontal_level, diff_coeff_phi_u,             &
                  diff_coeff_u, L_diff_wind, swap_field_is_vector,&
                  L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

  l_u_diff = .FALSE.
  CALL eg_pofil_vatp(                                               &
                  v, fld_type_v, 0, 1, 1, 0, 1, 1,                  &
                  model_levels, model_levels, model_levels,         &
                  first_constant_r_rho_level,                       &
                  n_rows, rows, rows, row_length,                   &
                  r_at_v, r_theta_at_v, r_uv_b, r_rho_levels,       &
                  off_x, off_y, off_x, off_y,                       &
                  halo_i, halo_j, off_x, off_y,                     &
                  off_x, off_y, halo_i, halo_j,                     &
                  sec_v_latitude, cos_theta_latitude,               &
                  n_procy, max_filter_rows, global_v_filter,        &
                  v_sweeps, v_begin, v_end,                         &
                  horizontal_level, diff_coeff_phi_v,               &
                  diff_coeff_v, L_diff_wind, swap_field_is_vector,  &
                  L_pofil_hadgem2, csxi2_v, csxi2_p, 0, l_u_diff)

  DEALLOCATE ( r_uv_b )
  DEALLOCATE ( r_theta_at_u )
  DEALLOCATE ( r_theta_at_v )

  ! ----------------------------------------------------------------------
  ! Section 2.3  Set u wind at poles since v has changed
  !              only for global and IF upper level diffusion is inactive
  ! ----------------------------------------------------------------------
  IF (model_type == mt_global .AND.                               &
     (top_filt_start > model_levels .OR. top_filt_start==imdi)) THEN
    IF ( at_extremity(PSouth) ) THEN

      CALL eg_v_at_poles(u,v, 1.0, udims%j_start,vdims%j_start,&
                   udims_s,vdims_s)
    END IF

    IF ( at_extremity(PNorth) ) THEN

      CALL eg_v_at_poles(u,v,-1.0, udims%j_end, vdims%j_end,&
                   udims_s,vdims_s)
    END IF

  END IF


  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(u,row_length,rows,model_levels,off_x,off_y,     &
                   fld_type_u,swap_field_is_vector)

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(v,row_length,n_rows,model_levels,off_x,off_y,   &
                   fld_type_v,swap_field_is_vector)

END IF ! L_pfuv .or. L_diff_wind

! ----------------------------------------------------------------------
! 3.0  Section for upper-level diffusion
! ----------------------------------------------------------------------

IF ( top_filt_start < model_levels + 1 .AND. &
                                     top_filt_start > imdi) THEN

  active_levels = top_filt_end - top_filt_start + 1
  j_start = 1
  j_stop = rows
  l_u_diff = .FALSE.
  ! DEPENDS ON: eg_diffupper
  CALL eg_diffupper(                                              &
                 thetav(:,:,1:model_levels), 0,                   &
                 model_levels, active_levels,                     &
                 rows, n_rows, row_length,                        &
                 off_x, off_y, off_x, off_y,                      &
                 sec_theta_latitude, cos_v_latitude,              &
                 j_start, j_stop, top_filt_start,                 &
                 up_diff, max_upd_levels, csxi2_p, csxi2_v, 1,    &
                 l_u_diff )

  ! DEPENDS ON: eg_diffupper
  CALL eg_diffupper(                                              &
                 w(1-off_x, 1-off_y, 1), 0,                       &
                 model_levels - 1, active_levels - 1,             &
                 rows, n_rows, row_length,                        &
                 off_x, off_y, off_x, off_y,                      &
                 sec_theta_latitude, cos_v_latitude,              &
                 j_start, j_stop, top_filt_start,                 &
                 up_diff, max_upd_levels, csxi2_p, csxi2_v, 1,    &
                 l_u_diff)

  l_u_diff = .TRUE.
  ! DEPENDS ON: eg_diffupper
  CALL eg_diffupper(                                              &
                 u, 0,                                            &
                 model_levels, active_levels,                     &
                 rows, n_rows, row_length,                        &
                 off_x, off_y, off_x, off_y,                      &
                 sec_theta_latitude, cos_v_latitude,              &
                 j_start, j_stop, top_filt_start,                 &
                 up_diff, max_upd_levels, csxi2_p, csxi2_v, 1,    &
                 l_u_diff )

  IF (model_type == mt_global) THEN
    IF (at_extremity(PSouth)) j_start = 2
    IF (at_extremity(PNorth)) j_stop = n_rows - 1
  END IF

  l_u_diff = .FALSE.
  ! DEPENDS ON: eg_diffupper
  CALL eg_diffupper(                                              &
                 v, 1,                                            &
                 model_levels, active_levels,                     &
                 n_rows, rows, row_length,                        &
                 off_x, off_y, off_x, off_y,                      &
                 sec_v_latitude, cos_theta_latitude,              &
                 j_start, j_stop, top_filt_start,                 &
                 up_diff, max_upd_levels, csxi2_v, csxi2_p, 0,    &
                 l_u_diff )

  ! ----------------------------------------------------------------------
  ! Section 3.1  Set u wind at poles since v has changed
  ! ----------------------------------------------------------------------

  IF (model_type == mt_global) THEN

    IF ( at_extremity(PSouth) ) THEN

      CALL eg_v_at_poles(u,v, 1.0, udims%j_start, vdims%j_start,&
                   udims_s,vdims_s)
    END IF

    IF ( at_extremity(PNorth) ) THEN

      CALL eg_v_at_poles(u,v, -1.0, udims%j_end, vdims%j_end,&
                   udims_s,vdims_s)
    END IF

  END IF ! model_type == mt_global

  i_field = 0

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => thetav(:,:,:)
  fields_to_swap(i_field) % field_type = fld_type_p
  fields_to_swap(i_field) % levels     = model_levels+1
  fields_to_swap(i_field) % rows       = rows
  fields_to_swap(i_field) % vector     = .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => w(:,:,:)
  fields_to_swap(i_field) % field_type = fld_type_p
  fields_to_swap(i_field) % levels     = model_levels+1
  fields_to_swap(i_field) % rows       = rows
  fields_to_swap(i_field) % vector     = .FALSE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => u(:,:,:)
  fields_to_swap(i_field) % field_type = fld_type_u
  fields_to_swap(i_field) % levels     = model_levels
  fields_to_swap(i_field) % rows       = rows
  fields_to_swap(i_field) % vector     = .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field      => v(:,:,:)
  fields_to_swap(i_field) % field_type = fld_type_v
  fields_to_swap(i_field) % levels     = model_levels
  fields_to_swap(i_field) % rows       = rows
  fields_to_swap(i_field) % vector     = .TRUE.

  CALL swap_bounds_mv(fields_to_swap, i_field, row_length,        &
                      off_x, off_y)

END IF !  top_filt_start < model_levels + 1

! ******************  increments for STASH   **************************

IF (sf(0,sect)) THEN

  ! T increment
  item = 381
  IF (sf(item,sect)) THEN

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          T_inc(i,j,k) = (thetav(i,j,k) - T_inc(i,j,k))            &
                                  *exner_theta_levels(i,j,k)
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
        T_inc,                                                    &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  END IF ! sf(item,sect)

  ! u wind increment
  item=385
  IF ( sf(item,sect) ) THEN

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          u_inc(i,j,k) = u(i,j,k) - u_inc(i,j,k)
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
        u_inc,                                                    &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  END IF ! sf(item,sect)

  ! v wind increment
  item=386
  IF ( sf(item,sect) ) THEN

    DO k=1,model_levels
      DO j=1,n_rows
        DO i=1,row_length
          v_inc(i,j,k) = v(i,j,k) - v_inc(i,j,k)
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
        v_inc,                                                    &
        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  END IF ! sf(item,sect)

  ! w wind increment
  item=387
  IF ( sf(item,sect) ) THEN

    DO k=0,model_levels
      DO j=1,rows
        DO i=1,row_length
          w_inc(i,j,k) = w(i,j,k) - w_inc(i,j,k)
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
        w_inc,                                                    &
        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  END IF ! sf(item,sect)

  ! exner increment
  item=388
  IF ( sf(item,sect) ) THEN

    DO k=1,model_levels+1
      DO j=1,rows
        DO i=1,row_length
          exner_inc(i,j,k) = exner(i,j,k) - exner_inc(i,j,k)
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
        exner_inc,                                                &
        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

  END IF ! sf(item,sect)

END IF ! sf(0,sect)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_NI_filter_Ctl
END MODULE eg_NI_filter_Ctl_mod
