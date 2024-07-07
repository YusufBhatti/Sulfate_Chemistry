! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_resetcon_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SISL_RESETCON_MOD'

CONTAINS
SUBROUTINE eg_sisl_resetcon(                                          &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,      &
         offx, offy,  l_datastart, l_slice,l_test_tracer,             &
         l_shallow, l_const_grav,                                     &
         z_top_of_model, thetav, rho, p_star, exner,                  &
         m_v, m_r, m_gr, m_cl, m_cf, m_cf2,                           &
         f1_comp, f2_comp, f3_comp, ih )

USE atm_fields_bounds_mod
USE planet_constants_mod,  ONLY: r, p_zero, kappa
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE timestep_mod,        ONLY: timestep,timestep_number
USE level_heights_mod,   ONLY: eta_theta_levels,                     &
                                eta_rho_levels,                       &
                                xi3_at_theta=>r_theta_levels,         &
                                xi3_at_rho=>r_rho_levels,             &
                                xi3_at_u=>r_at_u,                     &
                                xi3_at_v=>r_at_v,                     &
                                xi3_at_u_w=>r_at_u_w,                 &
                                xi3_at_v_w=>r_at_v_w
USE eg_helmholtz_mod
USE eg_set_helmholtz_mod
USE eg_dry_static_adj_ref_pro_mod
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE Field_Types
USE umPrintMgr
USE eg_alpha_mod
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

IMPLICIT NONE
!
! Description: Resets the reference profile to last timestep state
!
!
! Method: Reduced implementation of eg_SISL_setcon
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

REAL, INTENT(IN) :: ih

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,        &
                       l_datastart(3), halo_i, halo_j, offx,          &
                       offy

! Loop index bounds for arrays defined on p, u, v points respectively


LOGICAL, INTENT(IN) :: l_slice, l_test_tracer, l_shallow, l_const_grav

REAL, INTENT(IN) ::  z_top_of_model

REAL :: t_surf_ref2d(1-offx:row_length+offx,                            &
                   1-offy:rows+offy)                                  &
,     p_surf_ref2d(1-offx:row_length+offx,                            &
                   1-offy:rows+offy)

REAL :: t_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)                                   &
,     p_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)


REAL, INTENT(IN)  ::                                                  &
  f1_comp (row_length,rows),                                          &
  f2_comp (row_length,rows),                                          &
  f3_comp (row_length,rows)

REAL, INTENT(INOUT) ::                                                &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(IN)    ::                                                &
  m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
  m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
  m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

!Output Arrays from this routine
REAL, INTENT(INOUT) ::                                                &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)


! Local variables

INTEGER :: i,j,k, option_bubble

REAL :: dtheta_dz1_ref(3)

! Tolerance for initial and reference surface pressures to be considered
! to differ:
REAL, PARAMETER :: p_tol=0.1

! Tolerance for initial and reference surface temperatures to be
! considered to differ:
REAL, PARAMETER :: t_tol=0.01


REAL, INTENT(IN) :: p_star(1-offx:row_length+offx, 1-offy:rows+offy)


REAL :: kp2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SISL_RESETCON'

INTEGER :: filter

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

reference_profile_changed = .TRUE.

kp2    = (1.0-kappa)/kappa

!-------------------------------------------------------------------------
! Initialise theta, rho, exner fields
!-------------------------------------------------------------------------

! If we are running with a different background state we need
! to recompute the reference profiles.
! Might as well always recompute!

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*)  'EG_SISL_Resetcon: calculate reference profile'
  CALL umPrint(umMessage,src='SISL_ReSetcon_4A')
END IF

dtheta_dz1_ref      = 0.0 ! Reset the theta gradient to zero

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( model_levels, thetav_ref_pro, thetav, exner_ref_pro,     &
!$OMP         exner, rho_ref_pro, rho )                                &
!$OMP PRIVATE( k )
DO k = 1, model_levels
  thetav_ref_pro(:,:,k) = thetav(:,:,k)
  exner_ref_pro(:,:,k)  = exner (:,:,k)
  rho_ref_pro(:,:,k)    = rho(:,:,k)
END DO
!$OMP END PARALLEL DO
k = model_levels+1
exner_ref_pro(:,:,0)  = p_star(:,:)
exner_ref_pro(:,:,k)  = exner(:,:,k)
thetav_ref_pro(:,:,0) = thetav(:,:,0)

CALL eg_dry_static_adj_ref_pro(thetav_ref_pro)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)        &
!$OMP SHARED(pdims,rho_ref_pro,p_zero,r,intw_w2rho,thetav_ref_pro,     &
!$OMP        exner_ref_pro,kappa)
DO k = pdims%k_start, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      !   recompute density to satisfy equation of state
      rho_ref_pro (i,j,k) = p_zero/(r*                            &
                      (intw_w2rho(k,1)*thetav_ref_pro(i,j,k)      &
                     + intw_w2rho(k,2)*thetav_ref_pro(i,j,k-1) )) &
                     *(exner_ref_pro(i,j,k))**((1.0-kappa)/kappa)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO


CALL swap_bounds(exner_ref_pro,row_length,rows,                   &
                        model_levels+2,                           &
                        offx,offy,fld_type_p,swap_field_is_scalar)
CALL swap_bounds(thetav_ref_pro,                                       &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)
CALL swap_bounds(rho_ref_pro,                                          &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)

!----------------------------------------------------------------------
! Compute Helmholtz coefficients
!----------------------------------------------------------------------

CALL eg_set_helmholtz (                                               &
       row_length, rows, n_rows, model_levels, halo_i, halo_j,        &
       offx, offy,  l_slice, l_test_tracer, l_shallow,                &
       f1_comp, f2_comp, f3_comp, ih )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sisl_resetcon
END MODULE eg_sisl_resetcon_mod
