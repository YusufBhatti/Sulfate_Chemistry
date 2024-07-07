! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_set_init_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_IDL_SET_INIT_MOD'

CONTAINS
SUBROUTINE eg_idl_set_init(                                              &
                      row_length, rows, n_rows, halo_i, halo_j,          &
                      offx, offy, model_levels,                          &
                      datastart, global_row_length, global_rows,         &
                      nproc, g_row_length, g_rows,                       &
                      u, v, w, u_adv, v_adv, w_adv,  rho, exner,         &
                      exner_surf,                                        &
                      m_v_np1, m_cl_np1, m_cf_np1,                       &
                      m_r_np1, m_gr_np1, m_cf2_np1, exner_surf_np1,      &
                      z_top_of_model,                                    &
                      l_RK_dps, l_dry, alpha_w, ih,                      &
                      etadot, psi_w_surf, psi_w_lid,                     &
                      thetav, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

USE parkind1,   ONLY: jpim, jprb       !DrHook
USE yomhook,    ONLY: lhook, dr_hook   !DrHook

USE timestep_mod,           ONLY: timestep_number
USE Field_Types, ONLY: fld_type_p, fld_type_u, fld_type_v

USE atm_fields_bounds_mod, ONLY : pdims, pdims_s, tdims_s, tdims,          &
                                  udims_l, udims_s, udims,                 &
                                  vdims_l, vdims_s, vdims,                 &
                                  wdims_l, wdims_s
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar, swap_field_is_vector

USE profiles_mod,       ONLY: num_theta_relax_times, num_theta_relax_heights,  &
                              num_mv_relax_times, num_mv_relax_heights,        &
                              num_uv_relax_times, num_uv_relax_heights,        &
                              theta_relax_timescale, mv_relax_timescale,       &
                              uv_relax_timescale,                              &
                              theta_relax_height, theta_relax_time,            &
                              theta_relax_data,                                &
                              mv_relax_height, mv_relax_time,                  &
                              mv_relax_data,                                   &
                              uv_relax_height, uv_relax_time,                  &
                              u_relax_data, v_relax_data,                      &
                              theta_relax, mv_relax, u_relax, v_relax,         &
                              num_theta_inc_times, num_theta_inc_heights,      &
                              theta_inc_field_type,                            &
                              num_mv_inc_times,    num_mv_inc_heights,         &
                              num_uv_inc_times,    num_uv_inc_heights,         &
                              theta_inc_time,      theta_inc_height,           &
                              mv_inc_time,         mv_inc_height,              &
                              uv_inc_time,         uv_inc_height,              &
                              theta_inc_data,      mv_inc_data,                &
                              u_inc_data,          v_inc_data,                 &
                              num_w_force_times,   num_w_force_heights,        &
                              w_force_time, w_force_height, w_force_data,      &
                              theta_inc, mv_inc, u_inc, v_inc, w_subs,         &
                              num_uv_geo_times, num_uv_geo_heights,            &
                              uv_geo_height, uv_geo_time, u_geo_data,          &
                              v_geo_data, u_geostrophic, v_geostrophic

USE copy_profile_mod,   ONLY: copy_profile
USE eg_set_adv_winds_mod

USE departure_pts_mod, ONLY: reset_dpt_pts

USE dynamics_testing_mod, ONLY: problem_number, l_idealised_data

USE problem_mod,         ONLY: idealised_problem

USE init_etadot_mod

USE missing_data_mod, ONLY: rmdi
USE init_psiw_mod
USE idealise_run_mod, ONLY: l_shallow

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

IMPLICIT NONE
!
! Description:
!
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
REAL ::                                                                 &
  etadot(wdims_s%i_start:wdims_s%i_end,                                 &
         wdims_s%j_start:wdims_s%j_end,                                 &
         wdims_s%k_start:wdims_s%k_end),                                &
  m_v   (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end),                                &
  m_cl  (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end),                                &
  m_cf  (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end),                                &
  m_r   (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end),                                &
  m_gr  (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end),                                &
  m_cf2 (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end,                                 &
         tdims_s%k_start:tdims_s%k_end)
REAL :: thetav(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)

INTEGER :: row_length, rows, n_rows, model_levels
INTEGER :: offx, offy, halo_i, halo_j

REAL ::                                                                 &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows)

INTEGER :: global_row_length, global_rows, datastart(3), nproc
INTEGER :: g_row_length(0:nproc-1), g_rows(0:nproc-1)

REAL ::  z_top_of_model

! Arrays
REAL ::                                                           &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  u_adv(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,         &
        model_levels),                                            &
  v_adv(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,       &
       model_levels),                                             &
  w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
        0:model_levels),                                          &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),  &
  exner_surf(1-offx:row_length+offx, 1-offy:rows+offy),           &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL :: alpha_w, ih

REAL ::                                                               &
  exner_surf_np1(1-offx:row_length+offx, 1-offy:rows+offy),           &
  m_v_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  m_cl_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_cf_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_r_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
  m_gr_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  m_cf2_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

LOGICAL :: l_RK_dps  ! Runge-Kutta departure point flag

! Local

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='EG_IDL_SET_INIT'

LOGICAL :: l_dry     ! dry thermodynamics flag

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up idealised forcing profiles, if required.
IF (problem_number == idealised_problem .AND. l_idealised_data) THEN
! - Newtonian relaxation profiles
  CALL copy_profile(num_theta_relax_times, num_theta_relax_heights,            &
                    theta_relax_time, theta_relax_height,                      &
                    theta_relax_data, theta_relax,                             &
                    timescale = theta_relax_timescale)

  CALL copy_profile(num_mv_relax_times, num_mv_relax_heights,                  &
                    mv_relax_time, mv_relax_height,                            &
                    mv_relax_data, mv_relax,                                   &
                    timescale = mv_relax_timescale)

  CALL copy_profile(num_uv_relax_times, num_uv_relax_heights,                  &
                    uv_relax_time, uv_relax_height,                            &
                    u_relax_data, u_relax,                                     &
                    timescale = uv_relax_timescale)

  CALL copy_profile(num_uv_relax_times, num_uv_relax_heights,                  &
                    uv_relax_time, uv_relax_height,                            &
                    v_relax_data, v_relax,                                     &
                    timescale = uv_relax_timescale)

! - Increment (tendency) profiles
  CALL copy_profile(num_theta_inc_times, num_theta_inc_heights,                &
                    theta_inc_time, theta_inc_height,                          &
                    theta_inc_data, theta_inc,                                 &
                    field_type = theta_inc_field_type)

  CALL copy_profile(num_mv_inc_times, num_mv_inc_heights,                      &
                    mv_inc_time, mv_inc_height,                                &
                    mv_inc_data, mv_inc)

  CALL copy_profile(num_uv_inc_times, num_uv_inc_heights,                      &
                    uv_inc_time, uv_inc_height,                                &
                    u_inc_data, u_inc)

  CALL copy_profile(num_uv_inc_times, num_uv_inc_heights,                      &
                    uv_inc_time, uv_inc_height,                                &
                    v_inc_data, v_inc)

  ! w profile for applying subsidence forcing
  CALL copy_profile(num_w_force_times, num_w_force_heights,                    &
                    w_force_time, w_force_height,                              &
                    w_force_data, w_subs)

  ! Geostropic forcing profiles for U and V
  CALL copy_profile(num_uv_geo_times, num_uv_geo_heights,                      &
                    uv_geo_time, uv_geo_height,                                &
                    u_geo_data, u_geostrophic)

  CALL copy_profile(num_uv_geo_times, num_uv_geo_heights,                      &
                    uv_geo_time, uv_geo_height,                                &
                    v_geo_data, v_geostrophic)

END IF

! the reconfiguration cannot do this due to the timestep dependency
! (and the shallow flag - but we have ignored that elsewhere)
! therefore the reconfiguration fills it with RMDI, which we pick up here
! and fill it if required.
IF (  psi_w_surf(pdims%i_start,pdims%j_end) == rmdi ) THEN

  CALL swap_bounds(u, &
                   udims_s%i_len - 2*udims_s%halo_i, &
                   udims_s%j_len - 2*udims_s%halo_j, &
                   udims_s%k_len, &
                   udims_s%halo_i, udims_s%halo_j,   &
                   fld_type_u,swap_field_is_vector)
  CALL swap_bounds(v, &
                   vdims_s%i_len - 2*vdims_s%halo_i, &
                   vdims_s%j_len - 2*vdims_s%halo_j, &
                   vdims_s%k_len, &
                   vdims_s%halo_i, vdims_s%halo_j,   &
                   fld_type_v,swap_field_is_vector)
  CALL swap_bounds(w, &
                   wdims_s%i_len - 2*wdims_s%halo_i, &
                   wdims_s%j_len - 2*wdims_s%halo_j, &
                   wdims_s%k_len, &
                   wdims_s%halo_i, wdims_s%halo_j,   &
                   fld_type_p, swap_field_is_scalar)

  CALL init_psiw (psi_w_surf, psi_w_lid, l_shallow)
END IF

CALL swap_bounds(u_adv,                                                &
                 udims_l%i_len - 2*udims_l%halo_i,                     &
                 udims_l%j_len - 2*udims_l%halo_j,                     &
                 udims_l%k_len,                                        &
                 udims_l%halo_i, udims_l%halo_j,                       &
                 fld_type_u,swap_field_is_vector)
CALL swap_bounds(v_adv,                                                &
                 vdims_l%i_len - 2*vdims_l%halo_i,                     &
                 vdims_l%j_len - 2*vdims_l%halo_j,                     &
                 vdims_l%k_len,                                        &
                 vdims_l%halo_i, vdims_l%halo_j,                       &
                 fld_type_v,swap_field_is_vector)
CALL swap_bounds(w_adv,                                                &
                 wdims_l%i_len - 2*wdims_l%halo_i,                     &
                 wdims_l%j_len - 2*wdims_l%halo_j,                     &
                 wdims_l%k_len,                                        &
                 wdims_l%halo_i, wdims_l%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)


! Departure points need to be initialised at the first run, i.e.
! when Runge Kutta departure points are used (Polar rows only, technically)
! because then the departure points are not in the chain dump and the
! departure point algorithm will not compute polar v departure points.
! This is not an issue with interpolated departure points, as they are
! derived from the w departure points. Since this is only called once at
! startup of the model and does not take long it does not matter.

! initialise u,v,w departure points

CALL reset_dpt_pts()

! Ensure unused moisture fields are zero if not in use

IF( .NOT. l_mcr_qrain  ) m_r(:,:,:)   = 0.0
IF( .NOT. l_mcr_qcf2   ) m_cf2(:,:,:) = 0.0
IF( .NOT. l_mcr_qgraup ) m_gr(:,:,:)  = 0.0


m_v_np1     = m_v
m_cl_np1    = m_cl
m_cf_np1    = m_cf
m_r_np1     = m_r
m_gr_np1    = m_gr
m_cf2_np1   = m_cf2
exner_surf_np1   = exner_surf

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_set_init
END MODULE eg_idl_set_init_mod
