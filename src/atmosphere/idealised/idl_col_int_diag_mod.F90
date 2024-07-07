! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE idl_col_int_diag_mod

IMPLICIT NONE
  ! Description: Calculate column integral diagnostics from 
  !              forcing increments
  !
  ! Method: Column integrals done using information on rho levels as
  !         the model dry density is conserved using its integral on rho
  !         levels. If necessary fields on theta levels are interpolated in 
  !         the vertical to rho levels. 
  !         Choosing to calculate integrals due to du and dv on the
  !         U and V  horizontal grids as this will give more accurate values. 
  !         Interpolation to theta points would smooth the winds and reduce 
  !         the KE changes.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IDL_COL_INT_DIAG_MOD'

CONTAINS
SUBROUTINE idl_col_int_diag(l_mr_physics, u, v, dryrho, wetrho_r_sq_n )

USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod,     ONLY: tdims, tdims_s, pdims, pdims_s, &
                                     udims, udims_s, vdims, vdims_s, &
                                     tdims_l
USE field_types
USE level_heights_mod, ONLY: &
  r_theta_levels                ! Radii on theta levels (m)

USE timestep_mod,              ONLY: timestep
USE planet_constants_mod,      ONLY: r, cp
USE idealised_diag_mod,        ONLY: dt_inc_ideal_um, dq_inc_ideal_um, &
                                     du_inc_ideal_um, dv_inc_ideal_um, &
                                     dtheta_inc_ideal_um,              &
                                     dcolqdt_ideal_um, de_cvt_ideal_um, &
                                     de_u2_ideal_um, de_v2_ideal_um
USE eg_helmholtz_mod,          ONLY: ec_vol
USE horiz_grid_mod,            ONLY: intw_w2rho, cell_area_surface, intw_p2u,  &
                                     intw_p2v, xi1_p, xi2_p, xi1_u, xi2_v
USE metric_terms_mod,          ONLY: h1_p, h2_p, h3_p, h1_p_eta, h2_p_eta,     &
                                     h1_xi1_u, h2_xi2_v, deta_xi3_u, deta_xi3_v 
USE level_heights_mod,         ONLY: eta_theta_levels

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_mr_physics   ! .true. if mixing ratios 

REAL, INTENT (IN) ::                                 &
  u(udims_s%i_start:udims_s%i_end,                   & ! U wind (m/s)
    udims_s%j_start:udims_s%j_end,                   &
    udims_s%k_start:udims_s%k_end),                  &
  v(vdims_s%i_start:vdims_s%i_end,                   & ! V wind (m/s)
    vdims_s%j_start:vdims_s%j_end ,                  &
    vdims_s%k_start:vdims_s%k_end),                  &
  dryrho(pdims_s%i_start:pdims_s%i_end,              & ! dry rho (kg/m3)
         pdims_s%j_start:pdims_s%j_end,              &
         pdims_s%k_start:pdims_s%k_end),             &
  wetrho_r_sq_n(pdims_s%i_start:pdims_s%i_end,       & ! wet rho * r*r  (kg/m)
                pdims_s%j_start:pdims_s%j_end,       &
                pdims_s%k_start:pdims_s%k_end+1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IDL_COL_INT_DIAG'

INTEGER :: i, j, k



REAL    :: dqt_rho ! dqbydt on rho levels
REAL    :: dt_rho  ! dT inc on rho levels
REAL    :: cv      ! specific heat of dry air at constant volume
REAL    ::      &
  dx            & ! factor in x direction or EW
 ,dy            & ! factor in y direction or NS
 ,dz            & ! factor in z direction i.e. vertical
 ,tmp           & ! temporary variable
 ,du_ke         & ! Change in KE due to du
 ,dv_ke         & ! Change in KE due to dv
 ,u_vol         & ! volume of u gridbox
 ,v_vol         & ! volume of v gridbox
 ,rho_u         & ! dry density at u points
 ,rho_v         & ! dry density at v points
 ,u_surf_area   & ! Surface area of u gridbox
 ,v_surf_area     ! Surface area of v gridbox

!----------------------------------------------------------------------
! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Calculation of column q change per sec & change in cvT
! This calculation is using the zeroth level values. The assumption is
! these have sensible values from the above calculations.

cv = cp - r       ! value for dry air

DO k=1,tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dt_rho = intw_w2rho(k,1)*dt_inc_ideal_um(i,j,k)  +                     &
                 intw_w2rho(k,2)*dt_inc_ideal_um(i,j,k-1)
      ! Note choosing to use dry density as section 30 currently uses
      ! dry density in its integrals. 
      de_cvt_ideal_um(i,j) = de_cvt_ideal_um(i,j) +                          &
                               cv * dryrho(i,j,k) * dt_rho * ec_vol(i,j,k)
    END DO
  END DO
END DO

! Convert column integral quantities to be per unit area
DO j=pdims%j_start, pdims%j_end
  DO i=pdims%i_start, pdims%i_end
    de_cvt_ideal_um(i,j) = de_cvt_ideal_um(i,j) / cell_area_surface(i,j)
  END DO
END DO

! dq increment integrals depends on whether working in mixing ratios or
! specific quantities
IF (l_mr_physics) THEN    ! Mixing ratio increments
  DO k=1,tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dqt_rho = (intw_w2rho(k,1)*dq_inc_ideal_um(i,j,k)  +                 &
                   intw_w2rho(k,2)*dq_inc_ideal_um(i,j,k-1)) / timestep
        dcolqdt_ideal_um(i,j) = dcolqdt_ideal_um(i,j) +                      &
                                dryrho(i,j,k) * dqt_rho * ec_vol(i,j,k)
      END DO
    END DO
  END DO
  ! Convert column integral quantities to be per unit area
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      dcolqdt_ideal_um(i,j) = dcolqdt_ideal_um(i,j) / cell_area_surface(i,j)
    END DO
  END DO
ELSE                       ! Specific increments
  DO k=1,tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dqt_rho = (intw_w2rho(k,1)*dq_inc_ideal_um(i,j,k)  +                   &
                   intw_w2rho(k,2)*dq_inc_ideal_um(i,j,k-1)) / timestep
        ! Note cos(lat)*dlat*dphi term goes when going to per unit area
        dz =  r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
        dcolqdt_ideal_um(i,j) = dcolqdt_ideal_um(i,j) +                        &
                                wetrho_r_sq_n(i,j,k) * dqt_rho * dz
      END DO
    END DO
  END DO
  ! Convert column integral quantities to be per unit area by dividing
  ! by surface radius squared.
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      dcolqdt_ideal_um(i,j) = dcolqdt_ideal_um(i,j) /                   &
                      (r_theta_levels(i,j,0)*r_theta_levels(i,j,0))
    END DO
  END DO

END IF

! Change in KE due to dU i.e. u*du+0.5du*du
! Need density on u grid points in the horizontal
! Volume of gridbox about U point

DO k=1, udims%k_end
  dz = eta_theta_levels(k)-eta_theta_levels(k-1)
  DO j=udims%j_start, udims%j_end
    dy = xi2_v(j)-xi2_v(j-1)             
    DO i=udims%i_start, udims%i_end
      dx = xi1_p(i+1)-xi1_p(i)
      tmp = h1_xi1_u(i,j,k)*h2_p(i,j,k)*h3_p(i,j,k)*deta_xi3_u(i,j,k)
      u_vol = tmp*dx*dy*dz
      rho_u = intw_p2u(i,1)*dryrho(i,j,k)+intw_p2u(i,2)*dryrho(i+1,j,k)
      du_ke =( u(i,j,k)*du_inc_ideal_um(i,j,k) +                           &
                0.5*du_inc_ideal_um(i,j,k)*du_inc_ideal_um(i,j,k) )/timestep 
      de_u2_ideal_um(i,j) = de_u2_ideal_um(i,j) + rho_u *u_vol * du_ke
    END DO
  END DO
END DO
! Convert column integral quantities to be per unit area
DO j=udims%j_start, udims%j_end
  dy = xi2_v(j)-xi2_v(j-1)             
  DO i = udims%i_start, udims%i_end
    dx = xi1_p(i+1)-xi1_p(i)
    u_surf_area = h1_p_eta(i,j,0)*h2_p_eta(i,j,0)*dx*dy 
    de_u2_ideal_um(i,j) = de_u2_ideal_um(i,j) /u_surf_area 
  END DO
END DO

! Change in KE due to dV 
! Need density on V grid points in the horizontal
! Volume of gridbox about V point

DO k=1, vdims%k_end
  dz = eta_theta_levels(k)-eta_theta_levels(k-1)
  DO j=vdims%j_start, vdims%j_end
    dy = xi2_p(j+1)-xi2_p(j)             
    DO i=vdims%i_start, vdims%i_end
      dx = xi1_u(i)-xi1_u(i-1)
      tmp = h1_p(i,j,k)*h2_xi2_v(i,j,k)*h3_p(i,j,k)*deta_xi3_v(i,j,k)
      v_vol = tmp*dx*dy*dz
      rho_v = intw_p2v(j,1)*dryrho(i,j,k)+intw_p2v(j,2)*dryrho(i,j+1,k)
      dv_ke =( v(i,j,k)*dv_inc_ideal_um(i,j,k) +                           &
                0.5*dv_inc_ideal_um(i,j,k)*dv_inc_ideal_um(i,j,k) )/timestep 
      de_v2_ideal_um(i,j) = de_v2_ideal_um(i,j) + rho_v *v_vol * dv_ke
    END DO
  END DO
END DO
DO j=vdims%j_start, vdims%j_end
  dy = xi2_p(j+1)-xi2_p(j)             
  DO i=vdims%i_start, vdims%i_end
    dx = xi1_u(i+1)-xi1_u(i-1)
    v_surf_area = h1_p_eta(i,j,0)*h2_p_eta(i,j,0)*dx*dy 
    de_v2_ideal_um(i,j) = de_v2_ideal_um(i,j) /v_surf_area 
  END DO
END DO

!----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE idl_col_int_diag

END MODULE idl_col_int_diag_mod
