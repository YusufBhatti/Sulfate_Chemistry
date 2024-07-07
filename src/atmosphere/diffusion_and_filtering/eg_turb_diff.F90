! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_turb_diff

SUBROUTINE eg_turb_diff(                                                &
                        field, fld_type,                                &
                        row_length, rows,                               &
                        model_levels,                                   &
                        off_i, off_j,                                   &
                        coeff_u, coeff_v,                               &
                        delta_z, recip_r_squared_delz,                  &
                        L_swap, field_star)

USE timestep_mod, ONLY: timestep
USE atm_fields_bounds_mod, ONLY: pdims, udims, pdims_s, udims_s
USE horiz_grid_mod, ONLY: xi1_p, xi1_u, xi2_p, xi2_v, csxi2_p, csxi2_v

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
, rows                                                                  &
                     ! number of rows.
, model_levels                                                          &
                     ! number of model levels.
, off_i                                                                 &
                     ! Size of small halo for field_star
, off_j                                                                 &
                     ! Size of small halo for field_star
, fld_type      ! field type.

REAL, INTENT(IN) ::                                                     &
  coeff_u(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, model_levels)                  &
, coeff_v(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, model_levels)                  &
, delta_z(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, model_levels)                  &
, recip_r_squared_delz(pdims_s%i_start:pdims%i_end,                     &
                       pdims_s%j_start:pdims%j_end, model_levels)

REAL, INTENT(IN) ::                                                     &
  field (pdims_s%i_start:pdims_s%i_end,                                 &
         pdims_s%j_start:pdims_s%j_end, model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

REAL, INTENT(INOUT) ::                                                  &
  field_star (1-off_i:row_length+off_i,                                 &
              1-off_j:rows+off_j, model_levels )

! Local Variables.

INTEGER ::                                                              &
  i, j, k     ! Loop indices

! Local arrays

REAL ::                                                                 &
  temp(pdims_s%i_start:pdims%i_end,                                     &
       pdims_s%j_start:pdims%j_end)                                     &
, diff_term(pdims%i_start:pdims%i_end,                                  &
            pdims%j_start:pdims%j_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TURB_DIFF'

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
!$OMP& SHARED(pdims,field,delta_z,coeff_u,csxi2_p,xi1_p,xi2_p, &
!$OMP&   xi1_u,coeff_v,csxi2_v,field_star,timestep,recip_r_squared_delz, &
!$OMP& model_levels, xi2_v)
DO k = 1, model_levels

  ! ----------------------------------------------------------------------
  ! Section 2.1  Calculate lambda direction term.
  ! d/d_lambda( K_lambda/cos^2(phi) * dz * dQ/d_lambda )
  ! ----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start-1, pdims%i_end
! difference from p-point i -> i+1 is on u-point i
      temp(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *                   &
                                                 delta_z(i,j,k) *       &
                                                 coeff_u(i,j,k) /       &
                                         ( csxi2_p(j) * csxi2_p(j) *    &
                                          (xi1_p(i+1) - xi1_p(i)) )
    END DO

! difference from u-point i-1 -> i is on p-point i
    DO i = pdims%i_start, pdims%i_end
      diff_term(i,j) = ( temp(i,j) - temp(i-1,j) ) /                    &
                         ( xi1_u(i) - xi1_u(i-1) )
    END DO
  END DO

  ! ----------------------------------------------------------------------
  ! Section 2.2  Calculate phi direction term.
  ! 1/cos(phi) * d/d_phi( K_phi * cos(phi) * dz * dQ/d_phi )
  ! ----------------------------------------------------------------------

  DO j = pdims%j_start-1, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! difference from p-point j -> j+1 is on v-point j
      temp(i,j) = ( field(i,j+1,k)  -  field(i,j,k) ) *                 &
                                                   delta_z(i,j,k) *     &
                                                 coeff_v(i,j,k) *       &
                                                 csxi2_v(j) /           &
                                                 (xi2_p(j+1) - xi2_p(j))
    END DO
  END DO

! difference from v-point j-1 -> j is on p-point j
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      diff_term(i,j) = diff_term(i,j) +                                 &
                       ( temp(i,j) - temp(i,j-1) ) /                    &
                       ( csxi2_p(j) * (xi2_v(j) - xi2_v(j-1)) )
    END DO
  END DO

  ! ----------------------------------------------------------------------
  ! Section 2.3   Calculate new variable.
  ! ----------------------------------------------------------------------

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      field_star(i,j,k) = field_star(i,j,k) + timestep *                &
                             recip_r_squared_delz(i,j,k) *              &
                                        diff_term(i,j)
    END DO
  END DO

END DO ! k = 1, model_levels
!$OMP END PARALLEL DO

! DEPENDS ON: swap_bounds
IF (off_i > 0 )CALL swap_bounds(                                        &
                         field_star, row_length, rows, model_levels,    &
                         off_i, off_j, fld_type, L_swap)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_turb_diff
