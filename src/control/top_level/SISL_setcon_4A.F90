! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_setcon_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SISL_SETCON_MOD'

CONTAINS
SUBROUTINE eg_sisl_setcon(                                            &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,      &
         offx, offy,  l_datastart,                                    &
         l_shallow, l_const_grav,                                     &
         z_top_of_model,g_theta,g_rho,                                &
         pole_consts)

USE conversions_mod,       ONLY: pi
USE planet_constants_mod,  ONLY: g, planet_radius
USE init_vert_damp_mod,    ONLY: init_vert_damp
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

USE UM_ParVars,            ONLY: gc_proc_row_group, at_extremity

USE timestep_mod,          ONLY: timestep,                           &
                                 timestep_number

USE level_heights_mod,     ONLY: eta_theta_levels,                   &
                                 eta_rho_levels,                     &
                                 z_ref_theta, z_ref_rho,             &
                                 xi3_at_theta=>r_theta_levels,       &
                                 xi3_at_rho=>r_rho_levels,           &
                                 xi3_at_u=>r_at_u,                   &
                                 xi3_at_v=>r_at_v,                   &
                                 xi3_at_u_w=>r_at_u_w,               &
                                 xi3_at_v_w=>r_at_v_w
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE set_metric_terms_4A_mod
USE metric_terms_mod
USE ereport_mod, ONLY: ereport
USE field_types
USE UM_ParParams
USE umPrintMgr
USE atm_fields_bounds_mod
USE set_vert_interp_consts_mod
USE set_horiz_interp_consts_mod
USE dynamics_input_mod,   ONLY: damp_height, eg_vert_damp_coeff,     &
                                eg_vert_damp_profile, eta_s
USE dynamics_testing_mod, ONLY: problem_number, l_idealised_data
USE problem_mod,          ONLY: standard
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE lookup_table_mod

USE model_domain_mod, ONLY: model_type, mt_global
USE gravity_mod,      ONLY: g_ref_theta, g_ref_rho

USE science_fixes_mod, ONLY: l_eg_damp_height_lid

IMPLICIT NONE


!
! Description: computes the following constants for the SISL scheme:
!              metric terms, reference profiles, Helmholtz coefficients.
!
!
! Method:: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,    &
                       l_datastart(3), halo_i, halo_j, offx,      &
                       offy
!
LOGICAL, INTENT(IN) :: l_shallow, l_const_grav
!


REAL, INTENT(IN) ::  z_top_of_model

!Output Arrays from this routine

REAL, INTENT(OUT) :: pole_consts(4)

! Gravity arrays

REAL, INTENT(OUT) ::                                              &
  g_theta (1-offx:row_length+offx, 1-offy:rows+offy,              &
           0:model_levels)                                        &
, g_rho (1-offx:row_length+offx, 1-offy:rows+offy,                &
           model_levels)

! Local variables

INTEGER :: i,j,k, info

REAL :: dtheta_dz1_ref(3)

! Horizontal mesh spacing
REAL :: dxi1, dxi2

! Tolerance for initial and reference surface pressures to be considered
! to differ:
REAL, PARAMETER :: p_tol=0.1

! Tolerance for initial and reference surface temperatures to be
! considered to differ:
REAL, PARAMETER :: t_tol=0.01
REAL, PARAMETER :: delta = 1.0e-8
REAL :: dx, trig_chck, c, d, e, f, tmp_pole(0:row_length-1,4)

! Deep atmosphere
LOGICAL         :: l_deep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SISL_SETCON'

INTEGER :: ierr

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set 3D gravity acceleration arrays

IF (cartesian_grid .OR. l_shallow .OR. l_const_grav) THEN
  l_deep=.FALSE.
  IF (PrintStatus >= PrStatus_Normal) THEN
    ierr = -1
    CALL Ereport ('eg_SISL_setcon',ierr,                             &
                       ' Constant gravity enforced ')
  END IF
  k = 0
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      g_theta(i,j,k)=g
    END DO
  END DO

  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        g_theta(i,j,k)=g
        g_rho(i,j,k)  =g
      END DO
    END DO
  END DO
ELSE
  l_deep=.TRUE.
  k = 0
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      g_theta(i,j,k)=g*(planet_radius/xi3_at_theta(i,j,k))**2
    END DO
  END DO

  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        g_theta(i,j,k)=g*(planet_radius/xi3_at_theta(i,j,k))**2
        g_rho(i,j,k)  =g*(planet_radius/xi3_at_rho(i,j,k))**2
      END DO
    END DO
  END DO
END IF

CALL swap_bounds(g_theta,                                              &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(g_rho,                                                &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

IF (problem_number /= standard .AND. l_idealised_data) THEN
! Set 1D gravity acceleration reference arrays independent of orography
  g_ref_theta(0) = g
  IF ( l_deep ) THEN
    DO k=1, model_levels
      g_ref_theta(k) = g  * ( planet_radius /                          &
                            ( planet_radius + z_ref_theta(k) ) )**2
      g_ref_rho(k)   = g  * ( planet_radius /                          &
                            ( planet_radius + z_ref_rho(k) ) )**2
    END DO
  ELSE
    DO k=1, model_levels
      g_ref_theta(k) = g
      g_ref_rho(k)   = g
    END DO
  END IF
END IF

!-------------------------------------------------------------------------
! Compute metric terms
!-------------------------------------------------------------------------

CALL set_metric_terms_4A (l_shallow)

! Calculate constants used for v_at_poles (only done at the
! begining of the run so just use REPROD version

pole_consts = 0.0
tmp_pole    = 0.0
IF (.NOT. cartesian_grid .AND. model_type == mt_global) THEN
  IF ( at_extremity(pnorth) .OR. at_extremity(psouth) ) THEN
    DO i = udims%i_start, udims%i_end
      dx = ( xi1_p(i+1) - xi1_p(i) )/(2.0*pi)
      tmp_pole(i,1) = dx*COS(2.0*xi1_u(i))
      tmp_pole(i,2) = dx*SIN(2.0*xi1_u(i))
      tmp_pole(i,3) = dx*SIN(xi1_u(i))
      tmp_pole(i,4) = dx*COS(xi1_u(i))
    END DO
    CALL gcg_rvecsumr(row_length,row_length,1,4,tmp_pole,          &
                      gc_proc_row_group,info,pole_consts)
    c = pole_consts(1)
    d = pole_consts(2)
    e = pole_consts(3)
    f = pole_consts(4)
    trig_chck = SQRT((c + e**2 - f**2)**2 + (d-2.0*e*f)**2) + e**2 + f**2
    IF ( trig_chck > 1.0 - delta ) THEN
      ierr = 1
      CALL Ereport ('eg_SISL_setcon',ierr,' Problem at the poles ')
    END IF
  END IF
END IF

!-------------------------------------------------------------------------
! Set horizontal linear interpolation weights
!-------------------------------------------------------------------------

DO i=0, row_length
  intw_p2u(i,1) = (xi1_p(i+1)-xi1_u(i)) / (xi1_p(i+1)-xi1_p(i))
  intw_p2u(i,2) = 1.0 - intw_p2u(i,1)
END DO

DO j=0, n_rows
  intw_p2v(j,1) = (xi2_p(j+1)-xi2_v(j)) / (xi2_p(j+1)-xi2_p(j))
  intw_p2v(j,2) = 1.0 - intw_p2v(j,1)
END DO

DO i=1, row_length
  intw_u2p(i,1) = (xi1_u(i)-xi1_p(i)) / (xi1_u(i)-xi1_u(i-1))
  intw_u2p(i,2) = 1.0 -  intw_u2p(i,1)
END DO

DO j=1, rows
  intw_v2p(j,1) = (xi2_v(j)-xi2_p(j)) / (xi2_v(j)-xi2_v(j-1))
  intw_v2p(j,2) = 1.0 -intw_v2p(j,1)
END DO

!-------------------------------------------------------------------------
! Set vertical linear interpolation weights (eq C.10 ENDGame
! formulation vn 1.01)
!-------------------------------------------------------------------------
DO k=1, model_levels
  intw_w2rho(k,1) = ( eta_rho_levels(k)-eta_theta_levels(k-1) ) / &
                    ( eta_theta_levels(k)-eta_theta_levels(k-1) )
  intw_w2rho(k,2) = 1.0 - intw_w2rho(k,1)
END DO

DO k=1, model_levels-1
  intw_rho2w(k,1) = ( eta_theta_levels(k)-eta_rho_levels(k) ) /   &
                    ( eta_rho_levels(k+1)-eta_rho_levels(k) )
  intw_rho2w(k,2) = 1.0 - intw_rho2w(k,1)
END DO

intw_rho2w(model_levels,1) = 0.5
intw_rho2w(model_levels,2) = 0.5  ! as in New dynamics

!-------------------------------------------------------------------------
! Set area of grid cell at lower boundary
!-------------------------------------------------------------------------
DO j=pdims%j_start, pdims%j_end
  dxi2=xi2_v(j)-xi2_v(j-1)
  DO i=pdims%i_start, pdims%i_end
    dxi1=xi1_u(i)-xi1_u(i-1)
    cell_area_surface(i,j)=h1_p_eta(i,j,0)*h2_p_eta(i,j,0)*dxi1*dxi2
  END DO
END DO

!-------------------------------------------------------------------------
! Initialise vertical damping
!-------------------------------------------------------------------------
!
! If (temporary) logical l_eg_damp_height_lid is set then use a
! "reference" height of the top of the model. Otherwise, use the
! namelist input "damp_height"
!
IF (l_eg_damp_height_lid) THEN
  CALL init_vert_damp(eg_vert_damp_profile, eta_s,                  &
                    eg_vert_damp_coeff,z_top_of_model)
ELSE
  CALL init_vert_damp(eg_vert_damp_profile, eta_s,                  &
                    eg_vert_damp_coeff,damp_height)
END IF


! Set constants used by the interpolation routines

CALL set_vert_interp_consts(eta_rho_levels, eta_theta_levels, model_levels)
CALL set_horiz_interp_consts()
CALL set_look_up_table()

!-------------------------------------------------------------------------
! Compute Helmholtz coefficients
!-------------------------------------------------------------------------

CALL swap_bounds(h1_p,                                                 &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h2_p,                                                 &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h3_p,                                                 &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(phi_at_p,                                             &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h1_p_eta,                                             &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h2_p_eta,                                             &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h3_p_eta,                                             &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(phi_at_eta,                                           &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(h1_xi1_u,                                             &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u, swap_field_is_scalar)
CALL swap_bounds(h2_xi1_u,                                             &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u, swap_field_is_scalar)
CALL swap_bounds(h3_xi1_u,                                             &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u, swap_field_is_scalar)
CALL swap_bounds(phi_at_u,                                             &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u, swap_field_is_scalar)
CALL swap_bounds(h1_xi2_v,                                             &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
CALL swap_bounds(h2_xi2_v,                                             &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
CALL swap_bounds(h3_xi2_v,                                             &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
CALL swap_bounds(phi_at_v,                                             &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
CALL swap_bounds(deta_xi3,                                             &
                 pdims_s%i_len - 2*pdims_s%halo_i,                     &
                 pdims_s%j_len - 2*pdims_s%halo_j,                     &
                 pdims_s%k_len,                                        &
                 pdims_s%halo_i, pdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(deta_xi3_u,                                           &
                 udims_s%i_len - 2*udims_s%halo_i,                     &
                 udims_s%j_len - 2*udims_s%halo_j,                     &
                 udims_s%k_len,                                        &
                 udims_s%halo_i, udims_s%halo_j,                       &
                 fld_type_u, swap_field_is_scalar)
CALL swap_bounds(deta_xi3_v,                                           &
                 vdims_s%i_len - 2*vdims_s%halo_i,                     &
                 vdims_s%j_len - 2*vdims_s%halo_j,                     &
                 vdims_s%k_len,                                        &
                 vdims_s%halo_i, vdims_s%halo_j,                       &
                 fld_type_v, swap_field_is_scalar)
CALL swap_bounds(deta_xi3_theta,                                       &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(dxi1_xi3,                                             &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(dxi2_xi3,                                             &
                 tdims_s%i_len - 2*tdims_s%halo_i,                     &
                 tdims_s%j_len - 2*tdims_s%halo_j,                     &
                 tdims_s%k_len,                                        &
                 tdims_s%halo_i, tdims_s%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_setcon
END MODULE eg_sisl_setcon_mod
