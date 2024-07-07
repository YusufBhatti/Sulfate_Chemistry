! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_turb_diff_v

SUBROUTINE eg_turb_diff_v(                                              &
                        field, fld_type,                                &
                        row_length, n_rows,                             &
                        model_levels,                                   &
                        off_i, off_j,                                   &
                        coeff_centre, coeff_th,                         &
                        delta_z, recip_r_squared_delz,                  &
                        j_start, j_stop,                                &
                        L_swap, field_star)

USE timestep_mod, ONLY: timestep
USE horiz_grid_mod, ONLY: xi1_p, xi1_u, xi2_p, xi2_v, csxi2_p, csxi2_v
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, vdims, vdims_s

! Purpose:
!          Turbulent diffusion of a field
!          Based on conservative diffusion operator.
!          v-at-poles version
!
! Method:
!          Is described in ;
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

LOGICAL, INTENT(IN) ::  L_swap   ! swap_bounds logical

INTEGER, INTENT(IN) ::                                                  &
  row_length                                                            &
                     ! number of point on a row.
, n_rows                                                                &
                     ! number of rows.
, model_levels                                                          &
                     ! number of model levels.
, off_i                                                                 &
                     ! Size of small halo for field_star
, off_j                                                                 &
                     ! Size of small halo for field_star
, j_start, j_stop                                                       &
                     ! j start/end points
, fld_type      ! field type.

REAL, INTENT(IN) ::                                                     &
  coeff_centre(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end, model_levels)             &
, coeff_th(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end, model_levels)                 &
, delta_z(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, model_levels)                  &
, recip_r_squared_delz(pdims_s%i_start:pdims%i_end,                     &
                       pdims_s%j_start:pdims%j_end, model_levels)

REAL, INTENT(IN) ::                                                     &
  field (vdims_s%i_start:vdims_s%i_end,                                 &
         vdims_s%j_start:vdims_s%j_end, model_levels )

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

REAL, INTENT(INOUT) ::                                                  &
  field_star (vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end, model_levels )

! Local Variables.

INTEGER ::                                                              &
  i, j, k     ! Loop indices

! Local arrays

REAL ::                                                                 &
  temp(pdims_s%i_start:pdims%i_end,                                     &
       pdims_s%j_start:pdims%j_end)                                     &
, diff_term(vdims%i_start:vdims%i_end,                                  &
            vdims%j_start:vdims%j_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TURB_DIFF_V'

! No External Routines:

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 2.0  Horizontal Diffusion
! ----------------------------------------------------------------------
!
! "conservative diffusion" diffuses dzQ rather than Q.  Use
! dz*(dQ/d_lambda) instead of 1/d_lambda * (dzQ) (and similarly for phi)
! to avoid instabilities ocurring in the presence of orography due to
! volume differences between neighbouring gridboxes.

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,temp,diff_term) &
!$OMP& SHARED(j_start,j_stop,field,delta_z,coeff_centre,csxi2_v,xi1_p,xi2_v, &
!$OMP& xi1_u,coeff_th,csxi2_p,field_star,timestep,recip_r_squared_delz, &
!$OMP& model_levels, xi2_p, vdims)
DO k = 1, model_levels

  ! ----------------------------------------------------------------------
  ! Section 2.1  Calculate lambda direction term.
  ! d/d_lambda( K_lambda/cos^2(phi) * dz * dQ/d_lambda )
  ! ----------------------------------------------------------------------
  DO j = j_start, j_stop
    DO i = vdims%i_start-1, vdims%i_end
! difference between v-point i -> i+1 is on u-point i
      temp(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *                   &
                                                 delta_z(i,j,k) *       &
                                                 coeff_centre(i,j,k) /  &
                                         ( csxi2_v(j) * csxi2_v(j) *    &
                                          (xi1_p(i+1) - xi1_p(i)) )
    END DO

! difference between u point i-1 -> i is on v-point i
    DO i = vdims%i_start, vdims%i_end
      diff_term(i,j) = ( temp(i,j) - temp(i-1,j) ) /                    &
                         ( xi1_u(i) - xi1_u(i-1) )
    END DO
  END DO

  ! ----------------------------------------------------------------------
  ! Section 2.2  Calculate phi direction term.
  ! 1/cos(phi) * d/d_phi( K_phi * cos(phi) * dz * dQ/d_phi )
  ! ----------------------------------------------------------------------

  DO j = j_start, j_stop+1
    DO i = vdims%i_start, vdims%i_end
! difference from v-point j-1 -> j is on p-point j
      temp(i,j) = ( field(i,j,k)  -  field(i,j-1,k) ) *                 &
                                                   delta_z(i,j,k) *     &
                                                 coeff_th(i,j,k) *      &
                                                 csxi2_p(j) /           &
                                                 (xi2_v(j) - xi2_v(j-1))
    END DO
  END DO

! difference from p-point j -> j+1 is on v-point j
  DO j = j_start, j_stop
    DO i = vdims%i_start, vdims%i_end
      diff_term(i,j) = diff_term(i,j) +                                 &
                       ( temp(i,j+1) - temp(i,j) ) /                    &
                       ( csxi2_v(j) * (xi2_p(j+1) - xi2_p(j)) )
    END DO
  END DO

  ! ----------------------------------------------------------------------
  ! Section 2.3   Calculate new variable.
  ! ----------------------------------------------------------------------

  DO j = j_start, j_stop
    DO i = vdims%i_start, vdims%i_end
      field_star(i,j,k) = field_star(i,j,k) + timestep *                &
                             recip_r_squared_delz(i,j,k) *              &
                                        diff_term(i,j)
    END DO
  END DO

END DO ! k = 1, model_levels
!$OMP END PARALLEL DO

! DEPENDS ON: swap_bounds
IF (off_i > 0 )CALL swap_bounds(                                        &
                         field_star, row_length, n_rows, model_levels,  &
                         off_i, off_j, fld_type, L_swap)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_turb_diff_v
