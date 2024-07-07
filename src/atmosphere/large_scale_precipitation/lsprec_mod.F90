! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lsprec_mod

USE um_types, ONLY: real_lsprec

!Use statements for parameters
!General atmosphere modules
USE conversions_mod,      ONLY:                                               &
  zerodegc_orig => zerodegc,                                                  &
  pi_orig       => pi

USE water_constants_mod,  ONLY:                                               &
  lc_orig => lc,                                                              &
  lf_orig => lf

USE planet_constants_mod, ONLY:                                               &
  rv_orig => rv

USE pc2_constants_mod,    ONLY:                                               &
  wind_shear_factor_orig => wind_shear_factor

!Microphysics modules
USE mphys_ice_mod,        ONLY:                                               &
  t_scaling_orig => t_scaling,                                                &
  qcf0_orig      => qcf0,                                                     &
  t_agg_min_orig => t_agg_min,                                                &
  qcfmin_orig    => qcfmin,                                                   &
  thomo_orig     => thomo,                                                    &
  hm_t_min_orig  => hm_t_min,                                                 &
  hm_t_max_orig  => hm_t_max,                                                 &
  hm_decay_orig  => hm_decay,                                                 &
  hm_rqcl_orig   => hm_rqcl

USE lsp_dif_mod,          ONLY:                                               &
  cpwr_orig   => cpwr,                                                        &
  tw1_orig    => tw1,                                                         &
  tw2_orig    => tw2,                                                         &
  tw3_orig    => tw3,                                                         &
  tw4_orig    => tw4,                                                         &
  tw5_orig    => tw5,                                                         &
  tcor1_orig  => tcor1,                                                       &
  tcor2_orig  => tcor2

USE lsp_autoc_consts_mod, ONLY:                                               &
  auto_graup_qcf_thresh_orig  => auto_graup_qcf_thresh,                       &
  auto_graup_t_thresh_orig    => auto_graup_t_thresh,                         &
  auto_graup_coeff_orig       => auto_graup_coeff,                            &
  inhomog_rate_orig           => inhomog_rate,                                &
  inhomog_lim_orig            => inhomog_lim,                                 &
  r_thresh_orig               => r_thresh,                                    &
  r_auto_orig                 => r_auto,                                      &
  n_auto_orig                 => n_auto,                                      &
  power_droplet_auto_orig     => power_droplet_auto,                          &
  power_qcl_auto_orig         => power_qcl_auto,                              &
  power_rho_auto_orig         => power_rho_auto,                              &
  consts_auto_orig            => consts_auto,                                 &
  aut_pref_orig               => aut_pref,                                    &
  aut_qc_orig                 => aut_qc,                                      &
  aut_nc_orig                 => aut_nc

USE mphys_constants_mod,  ONLY:                                               &
  acc_pref_orig   => acc_pref,                                                &
  acc_qc_orig     => acc_qc,                                                  &
  acc_qr_orig     => acc_qr,                                                  &
  mprog_min_orig  => mprog_min,                                               &
  mprog_abs_orig  => mprog_abs,                                               &
  ntot_land_orig  => ntot_land,                                               &
  ntot_sea_orig   => ntot_sea,                                                &
  max_as_enh_orig => max_as_enh

USE mphys_ice_mod,        ONLY:                                               &
  m0_orig => m0

USE mphys_radar_mod,      ONLY:                                               &
  kliq_orig         => kliq,                                                  &
  kice_orig         => kice,                                                  &
  mm6m3_orig        => mm6m3,                                                 &
  ref_lim_orig      => ref_lim,                                               &
  ref_lim_lin_orig  => ref_lim_lin,                                           &
  mr_lim_orig       => mr_lim,                                                &
  cf_lim_orig       => cf_lim,                                                &
  nd_lim_orig       => nd_lim,                                                &
  ref_mom_orig      => ref_mom

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSPREC_MOD'

!Variables to be copied at lsprec precision, allowing full precision to be
!kept in their native modules as some are used all over the UM.

!Common magic numbers
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  zero = 0.0_real_lsprec,                                                     &
  half = 0.5_real_lsprec,                                                     &
  one  = 1.0_real_lsprec,                                                     &
  two  = 2.0_real_lsprec,                                                     &
  ten  = 10.0_real_lsprec

!Other PARAMETERS
!lsp_dif_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  cpwr  = REAL(cpwr_orig, real_lsprec),                                       &
  tw1   = REAL(tw1_orig, real_lsprec),                                        &
  tw2   = REAL(tw2_orig, real_lsprec),                                        &
  tw3   = REAL(tw3_orig, real_lsprec),                                        &
  tw4   = REAL(tw4_orig, real_lsprec),                                        &
  tw5   = REAL(tw5_orig, real_lsprec),                                        &
  tcor1 = REAL(tcor1_orig, real_lsprec),                                      &
  tcor2 = REAL(tcor2_orig, real_lsprec)

!lsp_autoc_consts_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  auto_graup_qcf_thresh = REAL(auto_graup_qcf_thresh_orig, real_lsprec),      &
  auto_graup_t_thresh   = REAL(auto_graup_t_thresh_orig, real_lsprec),        &
  auto_graup_coeff      = REAL(auto_graup_coeff_orig, real_lsprec),           &
  inhomog_rate          = REAL(inhomog_rate_orig, real_lsprec),               &
  inhomog_lim           = REAL(inhomog_lim_orig, real_lsprec),                &
  r_thresh              = REAL(r_thresh_orig, real_lsprec),                   &
  r_auto                = REAL(r_auto_orig, real_lsprec),                     &
  n_auto                = REAL(n_auto_orig, real_lsprec),                     &
  power_droplet_auto    = REAL(power_droplet_auto_orig, real_lsprec),         &
  power_qcl_auto        = REAL(power_qcl_auto_orig, real_lsprec),             &
  power_rho_auto        = REAL(power_rho_auto_orig, real_lsprec),             &
  consts_auto           = REAL(consts_auto_orig, real_lsprec),                &
  aut_pref              = REAL(aut_pref_orig, real_lsprec),                   &
  aut_qc                = REAL(aut_qc_orig, real_lsprec),                     &
  aut_nc                = REAL(aut_nc_orig, real_lsprec)

!mphys_constants_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  acc_pref    = REAL(acc_pref_orig, real_lsprec),                             &
  acc_qc      = REAL(acc_qc_orig, real_lsprec),                               &
  acc_qr      = REAL(acc_qr_orig, real_lsprec),                               &
  mprog_min   = REAL(mprog_min_orig, real_lsprec),                            &
  mprog_abs   = REAL(mprog_abs_orig, real_lsprec),                            &
  ntot_land   = REAL(ntot_land_orig, real_lsprec),                            &
  ntot_sea    = REAL(ntot_sea_orig, real_lsprec),                             &
  max_as_enh  = REAL(max_as_enh_orig, real_lsprec)

!conversions_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  zerodegc = REAL(zerodegc_orig, real_lsprec),                                &
  pi       = REAL(pi_orig, real_lsprec)

!mphys_ice_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  m0        = REAL(m0_orig, real_lsprec),                                     &
  t_scaling = REAL(t_scaling_orig, real_lsprec),                              &
  qcf0      = REAL(qcf0_orig, real_lsprec),                                   &
  t_agg_min = REAL(t_agg_min_orig, real_lsprec),                              &
  qcfmin    = REAL(qcfmin_orig, real_lsprec),                                 &
  thomo     = REAL(thomo_orig, real_lsprec),                                  &
  hm_t_min  = REAL(hm_t_min_orig, real_lsprec),                               &
  hm_t_max  = REAL(hm_t_max_orig, real_lsprec),                               &
  hm_decay  = REAL(hm_decay_orig, real_lsprec),                               &
  hm_rqcl   = REAL(hm_rqcl_orig, real_lsprec)

!water_constants_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  lc = REAL(lc_orig, real_lsprec),                                            &
  lf = REAL(lf_orig, real_lsprec)

!planet_constants_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  rv = REAL(rv_orig, real_lsprec)

!pc2_constants_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  wind_shear_factor = REAL(wind_shear_factor_orig, real_lsprec)

!mphys_radar_mod
REAL (KIND=real_lsprec), PARAMETER ::                                         &
  kliq        = REAL(kliq_orig, real_lsprec),                                 &
  kice        = REAL(kice_orig, real_lsprec),                                 &
  mm6m3       = REAL(mm6m3_orig, real_lsprec),                                &
  ref_lim     = REAL(ref_lim_orig, real_lsprec),                              &
  ref_lim_lin = REAL(ref_lim_lin_orig, real_lsprec),                          &
  mr_lim      = REAL(mr_lim_orig, real_lsprec),                               &
  cf_lim      = REAL(cf_lim_orig, real_lsprec),                               &
  nd_lim      = REAL(nd_lim_orig, real_lsprec),                               &
  ref_mom     = REAL(ref_mom_orig, real_lsprec)



!Non-parameter variables
!General atmosphere modules
REAL (KIND=real_lsprec) ::                                                    &
  !planet_constants_mod
  lcrcp, lfrcp, cp, r, repsilon, recip_epsilon, one_minus_epsilon, g,         &
  !timestep_mod
  timestep,                                                                   &
  !cloud_inputs_mod
  ice_width, cff_spread_rate,                                                 &
  !cderived_mod
  delta_lambda, delta_phi,                                                    &
  !fsd_parameters_mod
  f_cons(3), fsd_eff_lam, fsd_eff_phi,                                        &
  !rad_input_mod
  rad_mcica_sigma, two_d_fsd_factor,                                          &
  !stochastic_physics_run_mod
  ec_auto_rp

!Microphysics modules
REAL (KIND=real_lsprec) ::                                                    &
  !mphys_constants_mod
  timestep_mp, cx(200), constp(200), qclmin_rime, area_ratio_prefac,          &
  area_ratio_expn, rho_q_veloc, lam_evap_enh, ec_auto,                        &
  !lsp_dif_mod
  apb1, apb2, apb3, apb4, apb5, apb6,                                         &
  !mphys_inputs_mod
  ai, bi, tnuc, nscalesf

CONTAINS

SUBROUTINE lsprec_set_reals()

! Purpose:
!Module to manage conversions of all reals USEd below ls_ppnc to lsprec
!precision.

! Method:
! PARAMETERs are copied to PARAMETERs
! Normal REALS are copied to REALs

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: To follow

!All PARAMETERs are dealt with at the module level.
!All others are dealt with at the subroutine level.


!General atmosphere modules- not PARAMETER
USE planet_constants_mod, ONLY:                                               &
  lcrcp_orig              => lcrcp,                                           &
  lfrcp_orig              => lfrcp,                                           &
  cp_orig                 => cp,                                              &
  r_orig                  => r,                                               &
  repsilon_orig           => repsilon,                                        &
  recip_epsilon_orig      => recip_epsilon,                                   &
  one_minus_epsilon_orig  => one_minus_epsilon,                               &
  g_orig                  => g

USE timestep_mod,         ONLY:                                               &
  timestep_orig => timestep

USE cloud_inputs_mod,     ONLY:                                               &
  ice_width_orig        => ice_width,                                         &
  cff_spread_rate_orig  => cff_spread_rate

USE cderived_mod,         ONLY:                                               &
  delta_lambda_orig   => delta_lambda,                                        &
  delta_phi_orig      => delta_phi

USE fsd_parameters_mod,   ONLY:                                               &
  f_cons_orig         => f_cons,                                              &
  fsd_eff_lam_orig    => fsd_eff_lam,                                         &
  fsd_eff_phi_orig    => fsd_eff_phi

USE rad_input_mod,        ONLY:                                               &
  rad_mcica_sigma_orig  => rad_mcica_sigma,                                   &
  two_d_fsd_factor_orig => two_d_fsd_factor

USE stochastic_physics_run_mod, ONLY:                                         &
  ec_auto_rp_orig => ec_auto_rp

!Microphysics modules- not PARAMETER
USE mphys_constants_mod,  ONLY:                                               &
  timestep_mp_orig        => timestep_mp,                                     &
  cx_orig                 => cx,                                              &
  constp_orig             => constp,                                          &
  qclmin_rime_orig        => qclmin_rime,                                     &
  area_ratio_prefac_orig  => area_ratio_prefac,                               &
  area_ratio_expn_orig    => area_ratio_expn,                                 &
  rho_q_veloc_orig        => rho_q_veloc,                                     &
  lam_evap_enh_orig       => lam_evap_enh,                                    &
  ec_auto_orig            => ec_auto

USE lsp_dif_mod,          ONLY:                                               &
  apb1_orig => apb1,                                                          &
  apb2_orig => apb2,                                                          &
  apb3_orig => apb3,                                                          &
  apb4_orig => apb4,                                                          &
  apb5_orig => apb5,                                                          &
  apb6_orig => apb6

USE mphys_inputs_mod,     ONLY:                                               &
  ai_orig             => ai,                                                  &
  bi_orig             => bi,                                                  &
  tnuc_orig           => tnuc,                                                &
  nscalesf_orig       => nscalesf

IMPLICIT NONE

!===============================================================================
!Use in reals in lsprec precision, both microphysics related and general atmos
!USE lsprec_mod, ONLY:

!- logicals and integers
!===============================================================================

!In general, variables should be declared in the module so they can be used
!elsewhere

!End of header

!General atmosphere modules
  !planet_constants_mod
  lcrcp             = REAL(lcrcp_orig, real_lsprec)
  lfrcp             = REAL(lfrcp_orig, real_lsprec)
  cp                = REAL(cp_orig, real_lsprec)
  r                 = REAL(r_orig, real_lsprec)
  repsilon          = REAL(repsilon_orig, real_lsprec)
  recip_epsilon     = REAL(recip_epsilon_orig, real_lsprec)
  one_minus_epsilon = REAL(one_minus_epsilon_orig, real_lsprec)
  g                 = REAL(g_orig, real_lsprec)

  !timestep_mod
  timestep = REAL(timestep_orig, real_lsprec)

  !cloud_inputs_mod
  ice_width       = REAL(ice_width_orig, real_lsprec)
  cff_spread_rate = REAL(cff_spread_rate_orig, real_lsprec)

  !cderived_mod
  delta_lambda  = REAL(delta_lambda_orig, real_lsprec)
  delta_phi     = REAL(delta_phi_orig, real_lsprec)

  !fsd_parameters_mod
  f_cons(:)   = REAL(f_cons_orig(:), real_lsprec)
  fsd_eff_lam = REAL(fsd_eff_lam_orig, real_lsprec)
  fsd_eff_phi = REAL(fsd_eff_phi_orig, real_lsprec)

  !rad_input_mod
  rad_mcica_sigma   = REAL(rad_mcica_sigma_orig, real_lsprec)
  two_d_fsd_factor  = REAL(two_d_fsd_factor_orig, real_lsprec)

  !stochastic_physics_run_mod
  ec_auto_rp = REAL(ec_auto_rp_orig, real_lsprec)

!Microphysics modules
  !mphys_constants_mod
  timestep_mp       = REAL(timestep_mp_orig, real_lsprec)
  cx(:)             = REAL(cx_orig(:), real_lsprec)
  constp(:)         = REAL(constp_orig(:), real_lsprec)
  qclmin_rime       = REAL(qclmin_rime_orig, real_lsprec)
  area_ratio_prefac = REAL(area_ratio_prefac_orig, real_lsprec)
  area_ratio_expn   = REAL(area_ratio_expn_orig, real_lsprec)
  rho_q_veloc       = REAL(rho_q_veloc_orig, real_lsprec)
  lam_evap_enh      = REAL(lam_evap_enh_orig, real_lsprec)
  ec_auto           = REAL(ec_auto_orig, real_lsprec)

  !lsp_dif_mod
  apb1  = REAL(apb1_orig, real_lsprec)
  apb2  = REAL(apb2_orig, real_lsprec)
  apb3  = REAL(apb3_orig, real_lsprec)
  apb4  = REAL(apb4_orig, real_lsprec)
  apb5  = REAL(apb5_orig, real_lsprec)
  apb6  = REAL(apb6_orig, real_lsprec)

  !mphys_inputs_mod
  ai        = REAL(ai_orig, real_lsprec)
  bi        = REAL(bi_orig, real_lsprec)
  tnuc      = REAL(tnuc_orig, real_lsprec)
  nscalesf  = REAL(nscalesf_orig, real_lsprec)

END SUBROUTINE lsprec_set_reals
END MODULE lsprec_mod
