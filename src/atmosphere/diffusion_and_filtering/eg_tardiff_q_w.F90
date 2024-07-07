! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_tardiff_q_w
SUBROUTINE eg_tardiff_q_w                                         &
                     (field, w_pred,                              &
                      r_theta_levels, r_rho_levels,               &
                      off_x, off_y, halo_i, halo_j,               &
                      rows, n_rows, row_length,                   &
                      model_levels, w_limit,                      &
                      factor, test_level, end_level,              &
                      csxi2_p,csxi2_v,                            &
                      n_procy, L_pofil_hadgem2,                   &
                      sweeps, j_sweeps_begin, j_sweeps_end,       &
                      max_filter_rows, horizontal_level,          &
                      first_constant_r_rho_level_m1)

! Purpose:
!          Calculates conservative horizontal diffusion increment to q
!          subject to w > w_limit
!
! Method:
!          Is described in ;
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE trignometric_mod, ONLY: sec_theta_latitude,                  &
                             cos_v_latitude
USE eg_pofil_vatp_mod
USE UM_ParParams
USE mpp_conf_mod,      ONLY: swap_field_is_scalar
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

LOGICAL, INTENT(IN) ::                                            &
  L_pofil_hadgem2 ! setting for hadgem2

INTEGER, INTENT(IN) ::                                            &
  row_length,                                                     &
                   ! number of point on a row.
  rows,                                                           &
                   ! number of rows.
  n_rows,                                                         &
                   ! number of v rows.
  model_levels,                                                   &
                   ! number of model levels.
  halo_i,                                                         &
                ! Size of halo in i direction.
  halo_j,                                                         &
                ! Size of halo in j direction.
  off_x,                                                          &
                ! Size of small halo in i
  off_y,                                                          &
                ! Size of small halo in j.
  n_procy,                                                        &
                ! Number of processors in latitude
  first_constant_r_rho_level_m1

INTEGER, INTENT(IN) ::                                            &
  test_level,                                                     &
  end_level,                                                      &
  horizontal_level   ! level at which steep slope test no
!                          ! longer operates

INTEGER, INTENT(IN) ::                                            &
                     ! number of sweeps of filter to do
  max_filter_rows,                                                &
                     ! array size for u_begin etc
  j_sweeps_begin(0:max_filter_rows),                              &
                     ! row pointers for 1-2-1 filter
  j_sweeps_end(0:max_filter_rows),                                &
  sweeps(max_filter_rows)


REAL, INTENT(IN) ::                                               &
  w_limit,                                                        &
           ! Vertical velocity test value
  factor   ! effective diffusion coefficient

REAL, INTENT(IN) ::                                               &
     ! vertical co-ordinate arrays.
  r_theta_levels (1-halo_i:row_length+halo_i,                     &
                  1-halo_j:rows+halo_j, 0:model_levels),          &
  r_rho_levels (1-halo_i:row_length+halo_i,                       &
                1-halo_j:rows+halo_j, model_levels)


REAL, INTENT(IN) ::                                               &
      csxi2_p(1-halo_j:rows+halo_j),                              &
      csxi2_v(-halo_j:n_rows-1+halo_j)

REAL, INTENT(IN) ::                                               &
  w_pred(1-off_x:row_length+off_x,1-off_y:rows+off_y,             &
         0:model_levels)
  ! Wpred is some predictor for w eg wnp1 or wFSP1


! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
REAL, INTENT(INOUT) ::                                            &
  field(1-off_x:row_length+off_x,1-off_y:rows+off_y,              &
        0:model_levels)!,                                         &
!         w_local_mask(row_length,rows)

! Local Variables.

      INTEGER :: i, j, k

LOGICAL :: l_tar_diff = .FALSE.
LOGICAL :: l_diff
LOGICAL :: l_u_diff  ! diffusing u field

INTEGER :: COUNT = 0

! Local arrays

REAL ::                                                           &
  diffc_i  (1-off_x:row_length+off_x, 1-off_y:rows+off_y),        &
  diffc_j  (1-off_x:row_length+off_x, 1-off_y:rows+off_y),        &
  w_col_max(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

INTEGER ::                                                        &
                     ! array size for u_begin etc
  j_q_sweeps_begin(0:max_filter_rows),                            &
                     ! row pointers for 1-2-1 filter
  j_q_sweeps_end(0:max_filter_rows),                              &
  q_sweeps(max_filter_rows)

REAL, ALLOCATABLE ::                                              &
  r_theta_at_u(:,:,:),                                            &
  r_theta_at_v(:,:,:)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TARDIFF_Q_W'



! No External Routines:

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

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

! ----------------------------------------------------------------------
! Section 1.   Calculate diffusion coefficients
! ----------------------------------------------------------------------

!  Initialise diffusion coefficients to zero
DO j = 1-off_y, rows+off_y
  DO i = 1-off_x, row_length+off_x
    diffc_i(i,j) = 0.0
    diffc_j(i,j) = 0.0
    w_col_max(i,j)   = 0.0
  END DO  ! i = 1-off_x, row_length+off_x
END DO    ! j = 1-off_y, rows+off_y

!  Check to see if w > w_print_limit at each point
DO k =  test_level, end_level
  DO j = 1-off_y, rows+off_y
    DO i = 1-off_x, row_length+off_x
      !  Find if vertical velocity above threshold at this point
      !  start at level 5 since delta_z small near surface
      IF ( w_col_max(i,j) < w_pred(i,j,k) ) THEN
        w_col_max(i,j) = w_pred(i,j,k)
      END IF ! w_col_max(i,j) < w_adv(i,j,k

    END DO  ! i = 1-off_x, row_length+off_x
  END DO    ! j = 1-off_y, rows+off_y
END DO  !   k =  test_level, end_level

DO j = 1-off_y, rows+off_y
  DO i = 1-off_x, row_length+off_x

    IF ( w_col_max(i,j) > w_limit) THEN
      diffc_i(i,j) = factor
      diffc_j(i,j) = factor
      l_tar_diff = .TRUE.
      COUNT = COUNT + 1
    END IF ! w_col_max(i,j) > w_limit

  END DO  ! i = 1-off_x, row_length+off_x
END DO    ! j = 1-off_y, rows+off_y

DO j = 1, rows+off_y
  DO i = 1, row_length+off_x

    IF ( w_col_max(i,j) > w_col_max(i-1,j)) THEN
      diffc_i(i-1,j) = diffc_i(i,j)
    END IF ! w_col_max(i,j) > w_col_max(i-1,j)

    IF ( w_col_max(i,j) > w_col_max(i,j-1)) THEN
      diffc_j(i,j-1) = diffc_j(i,j)
    END IF ! w_col_max(i,j) > w_col_max(i,j-1)

  END DO  !i = 1, row_length+off_x
END DO    !j = 1, rows+off_y

! ----------------------------------------------------------------------
! Section 2.1  Call diffusion
! ----------------------------------------------------------------------
!         IF ( l_tar_diff ) WRITE(6,*) 'Targeted diffusion on ',count,' columns'


j_q_sweeps_begin(0) = 1
j_q_sweeps_end(0)   = rows

q_sweeps(1)         = 1

j_q_sweeps_begin(1) =  -1
j_q_sweeps_end(1)   = -3
DO j = 2,max_filter_rows
  j_q_sweeps_begin(j) =  -1
  j_q_sweeps_end(j)   = -3
  q_sweeps(j) =  0
END DO

! diffusion on, field is not a vector and thus not u.
l_diff = .TRUE.
l_u_diff = .FALSE.
CALL eg_pofil_vatp(                                              &
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
                   n_procy, max_filter_rows, 0,                  &
                   q_sweeps, j_q_sweeps_begin, j_q_sweeps_end,   &
                   horizontal_level, diffc_j,                    &
                   diffc_i, l_diff, swap_field_is_scalar,        &
                   L_pofil_hadgem2, csxi2_p, csxi2_v, 1,l_u_diff )

!-------------------------------------------------------------------------
DEALLOCATE ( r_theta_at_u )
DEALLOCATE ( r_theta_at_v )

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_tardiff_q_w
