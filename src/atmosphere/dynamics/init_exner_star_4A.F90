! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE init_exner_star_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_EXNER_STAR_MOD'

CONTAINS

SUBROUTINE init_exner_star(p, thetav,  m_v, m_cl, m_cf, m_r, m_gr,       &
                           m_cf2, p_star)

USE um_parvars,          ONLY: offx, offy, halo_i, halo_j
USE level_heights_mod,   ONLY: xi3_at_theta=>r_theta_levels,            &
                                xi3_at_rho=>r_rho_levels
USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: g, cp
USE atm_fields_bounds_mod
USE mpp_conf_mod,        ONLY: swap_field_is_scalar
USE Field_Types

IMPLICIT NONE
!
! Description:
!
!         Calculates surface pressure.
!
! Method: Adaptation of ENDGame code to put
!         surface into hydro-static balance
!         (essentially a simplified version of
!          eg_calc_p_star.F90)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_EXNER_STAR'


REAL, INTENT(IN) ::                                                   &
  thetav(tdims_s%i_start:tdims_s%i_end,                               &
         tdims_s%j_start:tdims_s%j_end,                               &
         tdims_s%k_start:tdims_s%k_end)

REAL :: p(pdims_s%i_start:pdims_s%i_end,                                 &
       pdims_s%j_start:pdims_s%j_end,                                 &
       pdims_s%k_start:pdims_s%k_end+1)

REAL, INTENT(IN) ::                                                   &
       m_v(tdims_s%i_start:tdims_s%i_end,                             &
           tdims_s%j_start:tdims_s%j_end,                             &
           tdims_s%k_start:tdims_s%k_end),                            &
       m_cl(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_cf(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_r(tdims_s%i_start:tdims_s%i_end,                             &
           tdims_s%j_start:tdims_s%j_end,                             &
           tdims_s%k_start:tdims_s%k_end),                            &
       m_gr(tdims_s%i_start:tdims_s%i_end,                            &
            tdims_s%j_start:tdims_s%j_end,                            &
            tdims_s%k_start:tdims_s%k_end),                           &
       m_cf2(tdims_s%i_start:tdims_s%i_end,                           &
             tdims_s%j_start:tdims_s%j_end,                           &
             tdims_s%k_start:tdims_s%k_end)

! Arguments with Intent OUT. ie: Output variables.

REAL, INTENT(INOUT) ::                                                &
  p_star (pdims_s%i_start:pdims_s%i_end,                              &
          pdims_s%j_start:pdims_s%j_end)


! Local Variables.

INTEGER  :: i, j,k       ! Loop indices
INTEGER  :: row_length, rows
REAL     :: a, b

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    b = 1.0 + m_v(i,j,0) + m_cl(i,j,0) + m_cf(i,j,0)                &
            + m_r(i,j,0) + m_gr(i,j,0) + m_cf2(i,j,0)
    a = (xi3_at_rho(i,j,1)-xi3_at_theta(i,j,0))                     &
             *b/( cp*thetav(i,j,0) )

    p_star(i,j) = p(i,j,1) + a*g

  END DO
END DO

! p_star is interpolated, so we do need the halo fill (apart from
! the fact that all fields with halos should have them filled or have no
! halos to start with!)

row_length = pdims%i_len
rows       = pdims%j_len

! DEPENDS ON: swap_bounds
CALL swap_bounds(p_star,row_length,rows,1,offx,offy,                &
                   fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_exner_star
END MODULE init_exner_star_mod
