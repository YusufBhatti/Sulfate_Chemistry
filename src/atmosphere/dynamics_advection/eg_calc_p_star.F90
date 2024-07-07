! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_calc_p_star_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CALC_P_STAR_MOD'

CONTAINS

SUBROUTINE eg_calc_p_star(                                            &
                 model_levels, row_length, rows,   p,                 &
                 thetav,  m_v, m_cl, m_cf, m_r, m_gr, m_cf2,          &
                 g,  p_star, w_surf, w_lid )

USE um_parvars,          ONLY: offx, offy, halo_i, halo_j
USE level_heights_mod,   ONLY: xi3_at_theta=>r_theta_levels,            &
                                xi3_at_rho=>r_rho_levels
USE metric_terms_mod,    ONLY:  h3_p_eta
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: cp
USE ref_pro_mod,         ONLY: p_ref_pro => exner_ref_pro
USE atm_fields_bounds_mod
USE Field_Types
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

IMPLICIT NONE
!
! Description:
!
!         Calculates surface pressure.
!
! Method: Place surface presure into balance so that
!         surface boundary condition holds
!         (ENDGame Documentation Section 7).
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CALC_P_STAR'


INTEGER, INTENT(IN) ::                                                &
  row_length,                                                         &
                   ! number of point on a row.
  rows,                                                               &
                   ! number of rows.
  model_levels

REAL, INTENT(IN) ::                                                   &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
   g(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels),         &
   w_surf(row_length,rows), w_lid(row_length,rows)

REAL :: p(1-offx:row_length+offx, 1-offy:rows+offy, model_levels+1)

REAL, INTENT(IN) ::      &
      m_v(1-offx:row_length+offx,1-offy:rows+offy, 0:model_levels),   &
      m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
      m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
      m_r(1-offx:row_length+offx,1-offy:rows+offy, 0:model_levels),   &
      m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
      m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Arguments with Intent OUT. ie: Output variables.

REAL, INTENT(INOUT) ::                                                &
  p_star (1-offx:row_length+offx, 1-offy:rows+offy)


! Local Variables.

INTEGER  :: i, j,k       ! Loop indices
REAL     :: a, b

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


k = model_levels
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    b = 1.0 + m_v(i,j,0) + m_cl(i,j,0) + m_cf(i,j,0)                &
            + m_r(i,j,0) + m_gr(i,j,0) + m_cf2(i,j,0)
    a = h3_p_eta(i,j,0)                                             &
             *(xi3_at_rho(i,j,1)-xi3_at_theta(i,j,0))               &
             *b/( cp*thetav(i,j,0) )

    p_star(i,j) = p(i,j,1) + a*( g(i,j,0) - w_surf(i,j) )


    b = 1.0 + m_v(i,j,k) + m_cl(i,j,k) + m_cf(i,j,k)                &
            + m_r(i,j,k) + m_gr(i,j,k) + m_cf2(i,j,k)

    a = -2.0*h3_p_eta(i,j,k)                                        &
             *(xi3_at_rho(i,j,k)-xi3_at_theta(i,j,k))               &
             *b/( cp*thetav(i,j,k) )

    p(i,j,k+1) = p(i,j,k) - a*( g(i,j,k) - w_lid(i,j) )

  END DO
END DO

! p_star is interpolated, so we do need the halo fill (apart from
! the fact that all fields with halos should have them filled or have no
! halos to start with!)

! DEPENDS ON: swap_bounds
CALL swap_bounds(p_star,row_length,rows,1,offx,offy,                &
                   fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_calc_p_star
END MODULE eg_calc_p_star_mod
