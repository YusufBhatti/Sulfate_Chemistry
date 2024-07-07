! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_NI_filter_incs_Ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_NI_FILTER_INCS_CTL_MOD'

CONTAINS

SUBROUTINE eg_NI_filter_incs_Ctl(                                 &
                         R_theta, R_u, R_v, R_w, R_rho,           &
                         S_theta, S_u, S_v, S_w,                  &
                         theta, u, v, w, rho, exner,              &
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
                         horizontal_level,                        &
                         off_x, off_y, halo_i, halo_j,            &
                         n_procy, at_extremity,                   &
                         L_polar_filter_incs, L_diff_incs,        &
                         L_pofil_hadgem2,                         &
                         l_eliminate_rho, outer_iter,             &
                         xi1_u, xi1_p, xi2_p,                     &
                         pole_consts, gc_proc_row_group,          &
                         global_row_length,                       &
                         csxi2_v, csxi2_p)

! Purpose:
!          CALL original polar filter code or new code
!          which also does horizontal diffusion as a filter
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE trignometric_mod, ONLY: cos_theta_latitude,                  &
                             sec_theta_latitude,                  &
                             cos_v_latitude,                      &
                             sec_v_latitude

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE eg_v_at_poles_mod
USE eg_pofil_vatp_mod
USE UM_ParParams

USE mpp_conf_mod,  ONLY: swap_field_is_vector, swap_field_is_scalar

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.
LOGICAL :: L_pofil_hadgem2  ! use hadgem2 polar filtering settings

LOGICAL ::                                                        &
  L_polar_filter_incs,                                            &
                       ! switch for polar filtering of increments
  L_diff_incs          ! switch for horiz. diffusion of increments

INTEGER ::                                                        &
  row_length,                                                     &
                   ! number of point on a row.
  rows,                                                           &
                   ! number of rows.
  n_rows,                                                         &
                   ! number of rows in a v field
  model_levels,                                                   &
                   ! number of model levels.
  first_constant_r_rho_level,                                     &
  first_constant_r_rho_level_m1,                                  &
  off_x,                                                          &
             ! Size of small halo in i
  off_y,                                                          &
             ! Size of small halo in j.
  halo_i,                                                         &
             ! Size of halo in i direction.
  halo_j,                                                         &
             ! Size of halo in j direction.
  l_datastart(3),                                                 &
             ! First gridpoints held by this processor
  me
             ! This processor



REAL ::                                                           &
!   diffusion coefficient for u/th rows
        diff_coeff_u(1-off_x:row_length+off_x,1-off_y:rows+off_y ),     &
!   diffusion coefficient for v rows
       diff_coeff_v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y ),    &
       diff_coeff_phi_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y),  &
       diff_coeff_phi_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

REAL :: diff_coeff_phi    ! NS diffusion coeff


INTEGER ::                                                        &
                               ! number of sweeps of filter to do
  max_filter_rows,                                                &
                     ! array size for u_begin etc
  global_u_filter,                                                &
                   ! sweep control for 1-2-1 filter
  global_v_filter,                                                &
                   ! sweep control for 1-2-1 filter
  u_begin(0:max_filter_rows),                                     &
                              ! row pointers for 1-2-1 filter
  u_end(0:max_filter_rows),                                       &
                              ! row pointers for 1-2-1 filter
  v_begin(0:max_filter_rows),                                     &
                              ! row pointers for 1-2-1 filter
  v_end(0:max_filter_rows),                                       &
                              ! row pointers for 1-2-1 filter
  u_sweeps(max_filter_rows),                                      &
                             ! sweeps for 1-2-1 filter
  v_sweeps(max_filter_rows),                                      &
                             ! sweeps for 1-2-1 filter
  horizontal_level
                         ! steep slope control


REAL ::                                                           &
     ! vertical co-ordinate arrays.
  r_theta_levels (1-halo_i:row_length+halo_i,                     &
                  1-halo_j:rows+halo_j, 0:model_levels),          &
  r_rho_levels (1-halo_i:row_length+halo_i,                       &
                1-halo_j:rows+halo_j, model_levels),              &
  r_at_u (1-halo_i:row_length+halo_i,                             &
          1-halo_j:rows+halo_j, model_levels),                    &
  r_at_v (1-halo_i:row_length+halo_i,                             &
          1-halo_j:n_rows+halo_j, model_levels)


INTEGER :: n_procy

LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Arguments with Intent IN/OUT.

! Indices have been shifted to allow field to be fixed sice
REAL, TARGET ::                                                   &
R_u(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels),        &
R_v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels),      &
R_w(row_length,rows,0:model_levels),                                  &
R_theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,0:model_levels),  &
R_rho(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

!     Fast Physics source terms
REAL , TARGET ::                                                 &
S_u    (1-off_x:row_length+off_x,1-off_y:  rows+off_y,model_levels), &
S_v    (1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels), &
S_w    (row_length,rows,0:model_levels),                             &
S_theta(1-off_x:row_length+off_x,1-off_y:  rows+off_y,0:model_levels)

REAL   ::                                                        &
w(1-off_x:row_length+off_x,1-off_y:rows+off_y,0:model_levels),        &
u(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels),          &
v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,model_levels),        &
theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,0:model_levels),    &
exner(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels+1),    &
rho(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels )



! Varibles used for eg_v_at_poles
REAL :: pole_consts(4)

! grid locations
REAL ::  xi1_p(1-halo_i:row_length+halo_i),                       &
         xi1_u(-halo_i:row_length-1+halo_i),                      &
         xi2_p(1-halo_j:rows+halo_j)

INTEGER, INTENT(IN) :: gc_proc_row_group, global_row_length

REAL ::                                                            &
      csxi2_p(1-halo_j:rows+halo_j),                               &
      csxi2_v(-halo_j:n_rows-1+halo_j)

LOGICAL, INTENT(IN) :: l_eliminate_rho

INTEGER, INTENT(IN) :: outer_iter

! !   Array arguments with intent(out):
!       REAL, INTENT(OUT) ::                                              &
!        STASHwork(*)   ! Output array holding diagnostic fields


! local  variables
INTEGER ::  i, j, k

LOGICAL :: l_u_diff   ! diffusing u field

! Local arrays
REAL ::                                                           &
  field( 1-off_x:row_length+off_x, 1-off_y:rows+off_y             &
                   , model_levels),                               &
  field_v( 1-off_x:row_length+off_x, 1-off_y:n_rows+off_y         &
                   , model_levels)

REAL, ALLOCATABLE ::                                              &
  r_theta_at_u(:,:,:),                                            &
  r_theta_at_v(:,:,:),                                            &
  r_uv_b(:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_NI_FILTER_INCS_CTL'


! ----------------------------------------------------------------------
! 0.1  Section for setting up diagnostic arrays
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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
! 2.0  Section for new polar filter
! ----------------------------------------------------------------------

IF ( L_polar_filter_incs .OR. L_diff_incs ) THEN

  IF ( outer_iter == 0 ) THEN ! diffuse timelevel n quantities
    !------------------------------------------------------------------------

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          field(i,j,k) = R_theta(i,j,k) - theta(i,j,k)
        END DO
      END DO
    END DO !  k = 1, model_levels


    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                             &
                     field, row_length, rows, model_levels,       &
                     off_x, off_y, fld_type_p, swap_field_is_scalar)
    
    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                           &
                    field(:,:,1:model_levels),                    &
                    fld_type_p, 0, 0, 0, 1, 0, 0,                 &
                    model_levels, model_levels, model_levels,     &
                    first_constant_r_rho_level_m1,                &
                    rows, n_rows, rows, row_length,               &
                    r_theta_levels, r_rho_levels,                 &
                    r_theta_at_u, r_theta_at_v,                   &
                    off_x, off_y, off_x, off_y,                   &
                    halo_i, halo_j, halo_i, halo_j,               &
                    off_x, off_y, off_x, off_y,                   &
                    sec_theta_latitude, cos_v_latitude,           &
                    n_procy, max_filter_rows, global_u_filter,    &
                    u_sweeps, u_begin, u_end,                     &
                    horizontal_level, diff_coeff_phi_u,           &
                    diff_coeff_u, L_diff_incs,swap_field_is_scalar,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          R_theta(i,j,k) = theta(i,j,k) + field(i,j,k)
        END DO
      END DO
    END DO


    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                               &
                       R_theta, row_length, rows, model_levels+1,   &
                       off_x, off_y, fld_type_p, swap_field_is_scalar)
    !------------------------------------------------------------------------

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          field(i,j,k) = R_rho(i,j,k) - rho(i,j,k)
        END DO
      END DO
    END DO !  k = 1, model_levels

    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                             &
                     field, row_length, rows, model_levels,       &
                     off_x, off_y, fld_type_p, swap_field_is_scalar)

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                           &
                    field, fld_type_p, 0, 0, 0, 1, 0, 0,          &
                    model_levels, model_levels, model_levels,     &
                    first_constant_r_rho_level_m1,                &
                    rows, n_rows, rows, row_length,               &
                    r_theta_levels, r_rho_levels,                 &
                    r_theta_at_u, r_theta_at_v,                   &
                    off_x, off_y, off_x, off_y,                   &
                    halo_i, halo_j, halo_i, halo_j,               &
                    off_x, off_y, off_x, off_y,                   &
                    sec_theta_latitude, cos_v_latitude,           &
                    n_procy, max_filter_rows, global_u_filter,    &
                    u_sweeps, u_begin, u_end,                     &
                    horizontal_level, diff_coeff_phi_u,           &
                    diff_coeff_u, L_diff_incs, swap_field_is_scalar,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          R_rho(i,j,k) = rho(i,j,k) + field(i,j,k)
        END DO
      END DO
    END DO
    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                             &
                       R_rho, row_length, rows, model_levels,     &
                       off_x, off_y, fld_type_p, swap_field_is_scalar)

    !------------------------------------------------------------------------

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          field(i,j,k) = R_w(i,j,k) - w(i,j,k)
        END DO
      END DO
    END DO  ! k = 1, model_levels - 1

    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                               &
                     field, row_length, rows, model_levels,         &
                     off_x, off_y, fld_type_p, swap_field_is_scalar)

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                             &
                    field(1-off_x,1-off_y, 1), fld_type_p,          &
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
                    diff_coeff_u, L_diff_incs, swap_field_is_scalar,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1,l_u_diff)

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          R_w(i,j,k) = w(i,j,k) + field(i,j,k)
        END DO
      END DO
    END DO  ! k = 1, model_levels - 1

    !------------------------------------------------------------------------

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

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          field(i,j,k) = R_u(i,j,k) - u(i,j,k)
        END DO
      END DO
    END DO
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(field,row_length,rows,model_levels,            &
                       off_x,off_y,fld_type_u,swap_field_is_vector)

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                             &
                    field, fld_type_u, 1, 0, 1, 0, 1, 1,            &
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
                    diff_coeff_u, L_diff_incs, swap_field_is_vector,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff )

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          R_u(i,j,k) = u(i,j,k) + field(i,j,k)
        END DO
      END DO
    END DO


    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(R_u,row_length,rows,model_levels,              &
                     off_x,off_y,fld_type_u,swap_field_is_vector)
    !------------------------------------------------------------------------
    DO k = 1, model_levels
      DO j = 1, n_rows
        DO i = 1, row_length
          field_v(i,j,k) = R_v(i,j,k) - v(i,j,k)
        END DO
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(field_v,row_length,n_rows,model_levels,        &
                     off_x,off_y,fld_type_v,swap_field_is_vector)

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                             &
                  field_v, fld_type_v, 0, 1, 1, 0, 1, 1,            &
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
                  diff_coeff_v, L_diff_incs, swap_field_is_vector,  &
                  L_pofil_hadgem2, csxi2_v, csxi2_p, 0,  l_u_diff)

    DO k = 1, model_levels
      DO j = 1, n_rows
        DO i = 1, row_length
          R_v(i,j,k) = v(i,j,k) + field_v(i,j,k)
        END DO
      END DO
    END DO
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(R_v,row_length,n_rows,model_levels,            &
                     off_x,off_y,fld_type_v,swap_field_is_vector)

    DEALLOCATE ( r_uv_b )
    !--------------------------------------------------------------------------
  ELSE ! outer_iter > 0

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                           &
                    S_theta(:,:,1:model_levels),                  &
                    fld_type_p, 0, 0, 0, 1, 0, 0,                 &
                    model_levels, model_levels, model_levels,     &
                    first_constant_r_rho_level_m1,                &
                    rows, n_rows, rows, row_length,               &
                    r_theta_levels, r_rho_levels,                 &
                    r_theta_at_u, r_theta_at_v,                   &
                    off_x, off_y, off_x, off_y,                   &
                    halo_i, halo_j, halo_i, halo_j,               &
                    off_x, off_y, off_x, off_y,                   &
                    sec_theta_latitude, cos_v_latitude,           &
                    n_procy, max_filter_rows, global_u_filter,    &
                    u_sweeps, u_begin, u_end,                     &
                    horizontal_level, diff_coeff_phi_u,           &
                    diff_coeff_u, L_diff_incs,swap_field_is_scalar,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

    !------------------------------------------------------------------------

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          field(i,j,k) = S_w(i,j,k)
        END DO
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL Swap_Bounds(                                               &
                     field, row_length, rows, model_levels,         &
                     off_x, off_y, fld_type_p, swap_field_is_scalar)

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                             &
                    field, fld_type_p,                              &
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
                    diff_coeff_u, L_diff_incs, swap_field_is_scalar,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff)

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          S_w(i,j,k) = field(i,j,k)
        END DO
      END DO
    END DO  ! k = 1, model_levels - 1

    ! S_w has no halo so doesn't need a swap bounds

    !------------------------------------------------------------------------

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

    l_u_diff  = .FALSE.
    CALL eg_pofil_vatp(                                             &
                    S_u, fld_type_u, 1, 0, 1, 0, 1, 1,              &
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
                    diff_coeff_u, L_diff_incs, swap_field_is_vector,&
                    L_pofil_hadgem2, csxi2_p, csxi2_v, 1, l_u_diff )

    !------------------------------------------------------------------------
    CALL eg_pofil_vatp(                                            &
                 S_v, fld_type_v, 0, 1, 1, 0, 1, 1,                &
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
                 diff_coeff_v, L_diff_incs, swap_field_is_vector,  &
                 L_pofil_hadgem2, csxi2_v, csxi2_p, 0, l_u_diff )

    DEALLOCATE ( r_uv_b )

  END IF ! outer_iter

END IF ! L_polar_filter_incs .or. L_diff_incs
!--------------------------------------------------------------------------

DEALLOCATE ( r_theta_at_u )
DEALLOCATE ( r_theta_at_v )

! ******************  increments for STASH   **************************

!       IF (sf(0,sect)) THEN
!
! ! T increment
!       item = 481
!         IF(sf(item,sect)) THEN
!
!           DO k=1,model_levels
!             DO j=1,rows
!               DO i=1,row_length
!                 T_inc(i,j,k) = (R_theta(i,j,k) - T_inc(i,j,k))       &
!                                         *0.5*(exner(i,j,k)+exner(i,j,k+1))
!               ENDDO  ! i
!             ENDDO  ! j
!           ENDDO  ! k
!
! ! DEPENDS ON: copydiag_3d
!           CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
!               T_inc,                                                    &
!               row_length,rows,model_levels,0,0,0,0,at_extremity,        &
!               stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
!               stash_levels,num_stash_levels+1,                          &
!               atmos_im,sect,item,                                       &
!               Errorstatus,cmessage)
!
!         ENDIF ! sf(item,sect)
!
! ! u wind increment
!         item=485
!         IF( sf(item,sect) ) THEN
!
!           DO k=1,model_levels
!             DO j=1,rows
!               DO i=1,row_length
!                 u_inc(i,j,k) = R_u(i,j,k) - u_inc(i,j,k)
!               ENDDO  ! i
!             ENDDO  ! j
!           ENDDO  ! k
!
! ! DEPENDS ON: copydiag_3d
!           CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
!               u_inc,                                                    &
!               row_length,rows,model_levels,0,0,0,0,at_extremity,        &
!               stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
!               stash_levels,num_stash_levels+1,                          &
!               atmos_im,sect,item,                                       &
!               Errorstatus,cmessage)
!
!         ENDIF ! sf(item,sect)
!
! ! v wind increment
!         item=486
!         IF( sf(item,sect) ) THEN
!
!           DO k=1,model_levels
!             DO j=1,n_rows
!               DO i=1,row_length
!                 v_inc(i,j,k) = R_v(i,j,k) - v_inc(i,j,k)
!               ENDDO  ! i
!             ENDDO  ! j
!           ENDDO  ! k
!
! ! DEPENDS ON: copydiag_3d
!           CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
!               v_inc,                                                    &
!               row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
!               stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
!               stash_levels,num_stash_levels+1,                          &
!               atmos_im,sect,item,                                       &
!               Errorstatus,cmessage)
!
!         ENDIF ! sf(item,sect)
!
! ! w wind increment
!         item=487
!         IF( sf(item,sect) ) THEN
!
!           DO k=1,model_levels
!             DO j=1,rows
!               DO i=1,row_length
!                 w_inc(i,j,k) = R_w(i,j,k) - w_inc(i,j,k)
!               ENDDO  ! i
!             ENDDO  ! j
!           ENDDO  ! k
!
! ! DEPENDS ON: copydiag_3d
!           CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
!               w_inc,                                                    &
!               row_length,rows,model_levels,0,0,0,0,at_extremity,        &
!               stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
!               stash_levels,num_stash_levels+1,                          &
!               atmos_im,sect,item,                                       &
!               Errorstatus,cmessage)
!
!         ENDIF ! sf(item,sect)
!
!       ENDIF ! sf(0,sect)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_NI_filter_incs_Ctl
END MODULE eg_NI_filter_incs_Ctl_mod
