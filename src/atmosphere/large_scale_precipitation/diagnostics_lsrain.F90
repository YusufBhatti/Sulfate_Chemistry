! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large scale precip
MODULE diagnostics_lsrain_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_LSRAIN_MOD'

CONTAINS

SUBROUTINE diagnostics_lsrain(                                                &
                       lspice_dim1,lspice_dim2,lspice_dim3,                   &
                       p_layer_centres, deltaz, rhodz_dry,                    &
                       t, q, qcl, qcf, qrain, qgraup, qcf2,                   &
                       cf, cfl, cff,                                          &
                       t_n, q_n, qcl_n, qcf_n, qrain_n, qgraup_n,             &
                       qcf2_n, cf_n, cfl_n, cff_n,                            &
                       ls_rain, ls_snow, ls_graup,                            &
                       ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,          &
                       rnout_tracer,lscav_dust_all,lscav_tr,                  &
                       lscav_nh3,                                             &
                       rnout_soot, lscav_soot,                                &
                       rnout_bmass, lscav_bmass,                              &
                       rnout_ocff, lscav_ocff,                                &
                       rnout_nitrate, lscav_nitrate,                          &
                       dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_r, dbz_l,           &
                       psdep,psaut,psacw,psacr,                               &
                       psaci,psmlt,psmltevp,                                  &
                       praut,pracw,prevp,                                     &
                       pgaut,pgacw,pgacs,pgmlt,                               &
                       pifrw,pifrr,piprm,pidep,piacw,                         &
                       piacr,pimlt,pimltevp,                                  &
                       pifall,psfall,prfall,pgfall,                           &
                       plset,plevpset, n_drop_tpr, n_drop_3d,                 &
                       frac_agg, mphys_pts,                                   &
                       vm_cry, vm_agg, vtbranch_flag, vm_used,                &
                       sfwater,sfrain,sfsnow,                                 &
                       qcl_mpt, tau_d, inv_prt, disprate,                     &
                       inv_mt, si_avg, dcfl_mp, sigma2_s,                     &
                       cfl_incr_mp,  qcl_incr_mp_pc2,                         &
                       stashwork                                              &
  )

USE pws_diags_mod, ONLY: pws_precip_sym_ls_rain, pws_precip_sym_ls_snow
USE um_stashcode_mod, ONLY: stashcode_pws_sec, stashcode_pws_precip_sym
USE submodel_mod,     ONLY: atmos_im
USE stash_array_mod,  ONLY:                                                   &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE gen_phys_inputs_mod,   ONLY: l_mr_physics

USE conversions_mod,       ONLY: zerodegc

USE acp_namel_mod,         ONLY: l_ac
USE ac_diagnostics_mod,    ONLY: lsrr, lssr, tinc_ppn
USE missing_data_mod,      ONLY: rmdi

USE mphys_radar_mod,       ONLY: ref_lim

USE mphys_diags_mod,       ONLY: l_warn_gr_am, l_warn_gr_ra, l_warn_gr_ra3

USE mphys_inputs_mod,      ONLY: l_mcr_qgraup

USE dust_parameters_mod,   ONLY: l_dust, ndiv, tot_dust_dep_flux

USE timestep_mod,          ONLY: timestep, timestep_number

USE um_parvars,            ONLY: at_extremity

USE max_hail_size_mod,     ONLY: max_hail_size_2d, max_hail_size_3d

! Grid bounds module
USE atm_fields_bounds_mod, ONLY: tdims

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    qsat_wat_new     => qsat_wat,                       &
                    l_new_qsat_mphys !Currently defaults to FALSE

USE errormessagelength_mod, ONLY: errormessagelength

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

! Purpose:
!          Calculates diagnostics generated from large scale
!          precipitation (UM section 4).

! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the large scale
! precipitation routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing, except where indicated.
! NOTE: Although moisture field diagnostics are available from this
! section (and correspond to equivalent variables calculated at UM4.5
! and earlier), the most appropriate place to extract moisture
! diagnostics is after section 9, which is the end of the moist
! processes calculations during the timestep.

!  Diagnostics currently available: (in order calculated)

! STASH item (all section 4 )
! ------------------------------------------------------------
! 100 aggregate fraction (model levels)
! 101 flag for where microphysics is performed        (model levels)
! 102 mean fallspeed crystals
! 103 mean fallspeed aggregates
! 104 mean fallspeed used
! 105 fallspeed branch flag
! 181 temperature increment across ls precip routines (model levels)
! 182 humidity    increment across ls precip routines (model levels)
! 183 qcl         increment across ls precip routines (model levels)
! 184 qcf         increment across ls precip routines (model levels)
! 189 qrain       increment across ls precip routines (model levels)
! 190 qgraup      increment across ls precip routines (model levels)
! 191 qcf2        increment across ls precip routines (model levels)
! 192 cf          increment across ls precip routines (model levels)
! 193 cfl         increment across ls precip routines (model levels)
! 194 cff         increment across ls precip routines (model levels)
! 201 large scale rain amount (kg/m2 per timestep)    (surface)
! 202 large scale snow amount (kg/m2 per timestep)    (surface)
! 203 large scale rainfall rate (kg/m2/s)             (surface)
! 204 large scale snowfall rate (kg/m2/s)             (surface)
!   4 temperature           after ls precip           (model levels)
! 205 cloud water (qcl)     after ls precip           (model levels)
! 206 cloud ice (qcf)       after ls precip           (model levels)
! 207 relative humidity     (percent)                 (model levels)
! 209 large scale graupel amount (kg/m2 per timestep) (surface)
! 210 cloud droplet number / m3 (from autoconversion) (model levels)
! 211 cloud droplet number / m3 (from aerosols)       (model levels)
! 212 large scale graupel rate (kg/m2/s)              (surface)
! 213 Large-scale rainout of dissolved ammonium nitrate (kg[N]/m2/s)
! 214 Large-scale washout of dissolved ammonium nitrate (kg[N]/m2/s)
!  10 specific humidity (q) after ls precip           (model levels)
! 220 Large scale rainout of soot (kg/m2/s)           (surface)
! 221 Large scale washout of soot (kg/m2/s)           (surface)
! 223 snowfall rate (kg/m2/s)                         (model levels)
! 224 supercooled liquid water content                (model levels)
! 225 supercooled rainfall rate                       (model levels)
! 226 large scale graupel rate (kg/m2/s)              (model levels)
! 227 rain fraction                                   (model levels)
! 228 Large scale rainout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 229 Large scale washout of fossil-fuel organic carbon (kg/m2/s) (srf)
! 237 Large scale rainout of biomass smoke (kg/m2/s)  (surface)
! 238 Large scale washout of biomass smoke (kg/m2/s)  (surface)
! Microphysical process rate diagnostics all on wet model levels
! 240 Homogeneous nucleation rate (kg/kg/s)
! 241 Heterogeneous nucleation rate (kg/kg/s)
! 243 Ice deposition rate (kg/kg/s)
! 245 Snow deposition rate (kg/kg/s)
! 247 Ice collection rate of cloud liquid water (riming) (kg/kg/s)
! 248 Snow collection rate of cloud liquid water (riming) (kg/kg/s)
! 249 Ice collection rate of rain (capture) (kg/kg/s)
! 250 Snow collection rate of rain (capture) (kg/kg/s)
! 251 Evaporation rate of melting ice (kg/kg/s)
! 252 Evaporation rate of melting snow (kg/kg/s)
! 253 Melting rate for ice (kg/kg/s)
! 254 Melting rate for snow (kg/kg/s)
! 255 Snow aggregate autoconversion rate (kg/kg/s)
! 256 Snow collection rate of ice (capture) (kg/kg/s)
! 257 Rain autoconversion rate (kg/kg/s)
! 258 Rain collection rate of cloud liquid water (accretion) (kg/kg/s)
! 259 Evaporation rate of rain (kg/kg/s)
! 260 Graupel autoconversion rate  (kg/kg/s)
! 261 Graupel collection rate of cloud water (accretion) (kg/kg/s)
! 262 Graupel collection rate of snow (capture) (kg/kg/s)
! 263 Melting rate for graupel (kg/kg/s)
! 265 Ice crystal sedimentation rate (kg/kg/s)
! 266 Snow sedimentation rate (kg/kg/s)
! 267 Rain sedimentation rate (kg/kg/s)
! 268 Graupel sedimentation rate (kg/kg/s)
! 269 Droplet settling rate of liquid (kg/kg/s)
! 270 Evaporation rate for settling droplets (kg/kg/s)
! 271 Homogeneous freezing of rain (kg/kg/s)
! 275 Surface maximum hail size (mm)
! 276 Maximum hail size in column (mm)
! 277 Maximum hail size on model levels (mm)
! 294 qcl increment by mixed-phase scheme (kg/kg)
! 295 cfl diagnosed by mixed phase scheme
! 296 Turbulent decorrelation timescale
! 297 Inverse phase-relaxation timescale
! 298 Diagnosed turbulent dissipation rate
! 299 Inverse mixing timescale
! 300 Mean of subgrid supersaturation PDF
! 301 Variance of subgrid supersaturation PDF
! 302 Snowfall no graupel amount (kg/m2/ts) (surface)
! 303 qcl increment from mixed phase scheme (PC2)
! 304 Snowfall no graupel rate (kg/m2/s) (surface)
! 313 cfl increment from mixed phase scheme (PC2)
! 323 Snowfall no graupel rate (kg/m2/s) (model levels)
! 118 Total Reflectivity (dBZ)
! 113 Reflectivity due to graupel (dBZ)
! 114 Reflectivity due to ice aggregates (dBZ)
! 115 Reflectivity due to ice crystals (dBZ)
! 116 Reflectivity due to rain (dBZ)
! 117 Reflectivity due to liquid cloud (dBZ)
! 110 Surface reflectivity (dBZ)
! 111 Max Reflectivity in the column (dBZ)
! 112 Reflectivity at 1km AGL (dBZ)

!
! NOTE : The following PC2 diagnostics are part of section 4 but
!        are calculated elsewhere:
! In diagnostics_pc2checks:  141, 142, 143, 130, 131, 144, 132, 133
!                            152, 153, 136, 137, 154, 138, 139
! In pc2_turbulence_control: 281, 282, 283, 292, 293

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

USE ereport_mod, ONLY: ereport
USE umprintmgr,  ONLY: newline

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER, INTENT(IN) ::                                                        &
  lspice_dim1,                                                                &
                          ! Dimensions for 3D diagnostic arrays.
  lspice_dim2,                                                                &
                          !  These are set to 1 in order to save
  lspice_dim3         !  memory if the diagnostics are not used.

! Primary Arrays used in all models
REAL, INTENT(IN) ::                                                           &
  p_layer_centres( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               0 : tdims%k_end ),                             &
  deltaz(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  rhodz_dry(       tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  t(               tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  q(               tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  qcl(             tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  qcf(             tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  qrain(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  qgraup(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  qcf2(            tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  cf(              tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  cfl(             tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  cff(             tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &

! Time level n values for increment diagnostics
    t_n  (           tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    q_n  (           tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    qcl_n(           tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    qcf_n(           tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    qrain_n  (       tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    qgraup_n(        tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    qcf2_n(          tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    cf_n  (          tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    cfl_n  (         tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    cff_n  (         tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    ls_rain(         tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end ),                           &
    ls_snow(         tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end ),                           &
    ls_graup(        tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end ),                           &
    ls_snow3d(  lspice_dim1, lspice_dim2, lspice_dim3 ),                      &
    rainfrac3d( lspice_dim1, lspice_dim2, lspice_dim3 ),                      &
    ls_graup3d( lspice_dim1, lspice_dim2, lspice_dim3 ),                      &
    n_drop_tpr(      tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    n_drop_3d(       tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ),                           &
    frac_agg(        tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end )

LOGICAL, INTENT(IN) :: mphys_pts( tdims%i_start : tdims%i_end,                &
                                  tdims%j_start : tdims%j_end,                &
                                              1 : tdims%k_end )

! Microphysics Process Rate diagnostic arrays
REAL, INTENT(IN) ::                                                           &
  psdep(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psaut(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psacw(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psacr(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psaci(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psmlt(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psmltevp(        tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )

REAL, INTENT(IN) :: dbz_tot(tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end,                      &
                                        1 : tdims%k_end )
                    ! Total radar reflectivity (dBZ)

REAL, INTENT(IN) :: dbz_g(tdims%i_start : tdims%i_end,                        &
                          tdims%j_start : tdims%j_end,                        &
                                      1 : tdims%k_end )
                    ! Radar reflectivity from graupel (dBZ)

REAL, INTENT(IN) :: dbz_i(tdims%i_start : tdims%i_end,                        &
                          tdims%j_start : tdims%j_end,                        &
                                      1 : tdims%k_end )
                    ! Radar reflectivity from ice aggregates (dBZ)

REAL, INTENT(IN) :: dbz_i2(tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end,                       &
                                       1 : tdims%k_end )
                    ! Radar reflectivity from ice crystals (dBZ)

REAL, INTENT(IN) :: dbz_r(tdims%i_start : tdims%i_end,                        &
                          tdims%j_start : tdims%j_end,                        &
                                      1 : tdims%k_end )
                    ! Radar reflectivity from rain (dBZ)

REAL, INTENT(IN) :: dbz_l(tdims%i_start : tdims%i_end,                        &
                          tdims%j_start : tdims%j_end,                        &
                                      1 : tdims%k_end )
                    ! Radar reflectivity from liquid cloud (dBZ)

REAL, INTENT(INOUT) ::                                                        &
  praut(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pracw(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  prevp(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )
REAL, INTENT(INOUT) ::                                                        &
  pgaut(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pgacw(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pgacs(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pgmlt(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )
REAL, INTENT(INOUT) ::                                                        &
  pifrw(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pifrr(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  piprm(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pidep(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  piacw(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  piacr(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pimlt(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pimltevp(        tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )
REAL, INTENT(INOUT) ::                                                        &
  pifall(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  psfall(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  prfall(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  pgfall(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )
REAL, INTENT(INOUT) ::                                                        &
  plset(           tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  plevpset(        tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )

REAL, INTENT(INOUT) ::                                                        &
  vm_cry(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  vm_agg(          tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  vtbranch_flag(   tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ),                             &
  vm_used(         tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )

! Seeder feeder sub-grid orographic water and precipitation
REAL, INTENT(INOUT) :: sfwater( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                            1 : tdims%k_end )
REAL, INTENT(INOUT) :: sfrain( tdims%i_start : tdims%i_end,                   &
                               tdims%j_start : tdims%j_end,                   &
                                           1 : tdims%k_end )
REAL, INTENT(INOUT) :: sfsnow( tdims%i_start : tdims%i_end,                   &
                               tdims%j_start : tdims%j_end,                   &
                                           1 : tdims%k_end )

! Turbulent generation of mixed phase diagnostics
REAL, INTENT(INOUT) :: qcl_mpt(  tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: tau_d(    tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: inv_prt(  tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: disprate( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: inv_mt(   tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: si_avg(   tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: dcfl_mp(  tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: sigma2_s( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: cfl_incr_mp( tdims%i_start : tdims%i_end,              &
                                    tdims%j_start : tdims%j_end,              &
                                                1 : tdims%k_end )

REAL, INTENT(INOUT) :: qcl_incr_mp_pc2( tdims%i_start : tdims%i_end,          &
                                        tdims%j_start : tdims%j_end,          &
                                                    1 : tdims%k_end )

! Used as input and workspace
REAL, INTENT(INOUT) ::                                                        &
  ls_rain3d(     lspice_dim1, lspice_dim2, lspice_dim3 ),                     &
  rnout_tracer(  lspice_dim1, lspice_dim2 ),                                  &
  lscav_dust_all(tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end, ndiv ),                         &
                                           !scavenged mineral dust
  lscav_tr(      lspice_dim1, lspice_dim2 ),                                  &
  lscav_nh3(     lspice_dim1, lspice_dim2 ),                                  &
  rnout_soot(    lspice_dim1, lspice_dim2 ),                                  &
  lscav_soot(    lspice_dim1, lspice_dim2 ),                                  &
  rnout_bmass(   lspice_dim1, lspice_dim2 ),                                  &
  lscav_bmass(   lspice_dim1, lspice_dim2 ),                                  &
  rnout_ocff(    lspice_dim1, lspice_dim2 ),                                  &
  lscav_ocff(    lspice_dim1, lspice_dim2 ),                                  &
  rnout_nitrate( lspice_dim1, lspice_dim2 ),                                  &
  lscav_nitrate( lspice_dim1, lspice_dim2 )


! Diagnostic variables
REAL, INTENT(INOUT) ::                                                        &
 stashwork(*)     ! STASH workspace for section 4 (LS precip)

! Local variables
INTEGER ::                                                                    &
 i, j, k, ji,                                                                 &
     icode                ! Return code  = 0 Normal exit  >0 Error
                          !              < 0 Warning

INTEGER :: sect,item    ! STASH section, item no.s
PARAMETER (sect = 4) !  for microphysics - large scale rain

REAL :: work_3d( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end ) ! 3D work space

REAL :: work_2d( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end ) ! 2D work space

REAL :: dbz_surf( tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )
 ! Surface (model level 1) reflectivity

REAL :: dbz_1km(  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )
 ! Reflectivity at 1km AGL.

REAL :: dbz_max(  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end )
 ! Max Reflectivity in the column.

REAL :: z_agl( tdims%i_start : tdims%i_end,                                   &
               tdims%j_start : tdims%j_end )
 ! Height above ground level (metres)

LOGICAL :: dbz_def ( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end )
 ! Determines whether reflectivity needs definining or not

REAL, PARAMETER :: agl_1km = 1000.0 ! 1km AGL in metres

! Arrays for maximum hail size diagnostics:
! Surface maximum hail [mm]
REAL :: surface_hail_diam( tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end )

! Maximum graupel water content in column [kg m-3] - not output
! as a diagnostic, but used for diagnostic calculations.

REAL :: gwc_max( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end )

! Maximum hail in column [mm]
REAL :: max_col_hail_diam( tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end )

! Maximum hail diameter on model levels [mm]
REAL :: max_hail_diam( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end )

REAL :: rho_work(  tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end )

CHARACTER(LEN=errormessagelength):: cmessage

CHARACTER(LEN=*):: RoutineName
PARAMETER ( RoutineName='DIAGNOSTICS_LSRAIN')

INTEGER ::                                                                    &
  im_index        ! internal model index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0 ! Initialise error status
im_index = 1

! Set up a working dry density (rho_work) if diagnostics need it
IF ( sf(275, sect) .OR. sf(276, sect) .OR. sf(277, sect) ) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,rho_work,rhodz_dry,deltaz)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rho_work(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF ! sf diagnostics which require rho_work

!  Copy diagnostic information to STASHwork for STASH processing

!--------------------------------------------------------------------------
! Item 100: Aggregate Fraction
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(100,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(100,4,im_index)), frac_agg,                  &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,100,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,100,                                                        &
       icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(frac_agg)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 101: Flag for where microphysics is performed
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(101,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(101,4,im_index)), mphys_pts,                 &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,101,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,101,                                                        &
       icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(mphys_pts)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 102: Fallspeed crystal parameters
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(102,4)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(102,4,im_index)), vm_cry,                    &
                     tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,           &
                     at_extremity,                                            &
                     stlist(1,stindex(1,102,4,im_index)),len_stlist,          &
                     stash_levels,num_stash_levels+1,                         &
                     atmos_im,4,102,                                          &
                     icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(vm_cry)"
  END IF

END IF


!--------------------------------------------------------------------------
! Item 103: Fallspeed aggregate parameters
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(103,4)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(103,4,im_index)), vm_agg,                    &
                     tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,           &
                     at_extremity,                                            &
                     stlist(1,stindex(1,103,4,im_index)),len_stlist,          &
                     stash_levels,num_stash_levels+1,                         &
                     atmos_im,4,103,                                          &
                     icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(vm_agg)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 104: Fallspeed branch used
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(104,4)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(104,4,im_index)), vtbranch_flag,             &
                     tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,           &
                     at_extremity,                                            &
                     stlist(1,stindex(1,104,4,im_index)),len_stlist,          &
                     stash_levels,num_stash_levels+1,                         &
                     atmos_im,4,104,                                          &
                     icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(vtbranch_flag)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 105: Fallspeed used
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(105,4)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(105,4,im_index)), vm_used,                   &
                     tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,           &
                     at_extremity,                                            &
                     stlist(1,stindex(1,105,4,im_index)),len_stlist,          &
                     stash_levels,num_stash_levels+1,                         &
                     atmos_im,4,105,                                          &
                     icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(vm_used)"
  END IF

END IF


! increment diagnostics= modified - previous

item = 181  ! temperature increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,t,t_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        work_3d(i,j,k) = t(i,j,k) - t_n(i,j,k)

      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

END IF

IF (icode <= 0 .AND. (sf(item,sect) .OR. l_ac)) THEN
  IF (.NOT. ALLOCATED(tinc_ppn)) THEN
    ALLOCATE ( tinc_ppn(tdims%i_end*tdims%j_end,tdims%k_end) )
    tinc_ppn(1,1) = rmdi
  END IF
END IF

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0,                           &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 181)"//cmessage
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,ji)                                                       &
!$OMP SHARED(tdims,tinc_ppn,work_3d)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ji = (j-1)*tdims%i_end+i

        tinc_ppn(ji,k) = work_3d(i,j,k)

      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

END IF  !  sf(item,sect)

item = 182  ! humidity increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,q,q_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = q(i,j,k) - q_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 182)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 183  ! qcl increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcl,qcl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qcl(i,j,k) - qcl_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 183)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 170  ! qcl increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcl,qcl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MAX(0.0,qcl(i,j,k) - qcl_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 170)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 171  ! qcl increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcl,qcl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MIN(0.0,qcl(i,j,k) - qcl_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 171)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 184  ! qcf increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcf,qcf_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qcf(i,j,k) - qcf_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 184)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 172  ! qcf increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcf,qcf_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MAX(0.0,qcf(i,j,k) - qcf_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 172)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 173  ! qcf increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcf,qcf_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MIN(0.0,qcf(i,j,k) - qcf_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,           &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 173)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 189  ! qrain increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qrain,qrain_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qrain(i,j,k) - qrain_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 189)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 190  ! qgraup increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qgraup,qgraup_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qgraup(i,j,k) - qgraup_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
      work_3d,                                                                &
      tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                          &
      at_extremity,                                                           &
      stlist(1,stindex(1,item,sect,im_index)),len_stlist,                     &
      stash_levels,num_stash_levels+1,                                        &
      atmos_im,sect,item,                                                     &
      icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 190)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 191  ! qcf2 increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qcf2,qcf2_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qcf2(i,j,k) - qcf2_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 191)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 192  ! cf increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cf,cf_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = cf(i,j,k) - cf_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 192)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 193  ! cfl increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cfl,cfl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = cfl(i,j,k) - cfl_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 193)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 176  ! cfl increment:positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cfl,cfl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MAX(0.0,cfl(i,j,k) - cfl_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 176)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 177  ! cfl increment:negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cfl,cfl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MIN(0.0,cfl(i,j,k) - cfl_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 177)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 194  ! cff increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cff,cff_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = cff(i,j,k) - cff_n(i,j,k)
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 194)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 178  ! cff increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cff,cff_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MAX(0.0,cff(i,j,k) - cff_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 178)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 179  ! cff increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,cff,cff_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = MIN(0.0,cff(i,j,k) - cff_n(i,j,k))
      END DO ! i
    END DO   ! j
  END DO     ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 179)"//cmessage
  END IF

END IF  !  sf(item,sect)

! Item 201 Large scale rain

IF (icode <= 0 .AND. sf(201,4)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(201,4,im_index)),ls_rain,                        &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,201,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(ls_rain)"
  END IF
  ! Code to convert rate to amount for a given timestep

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i)                                                              &
!$OMP SHARED(tdims,stashwork,si,im_index,timestep)
  DO i = 1, tdims%i_end*tdims%j_end
    stashwork(si(201,4,im_index)+i-1)=                                        &
         stashwork(si(201,4,im_index)+i-1)* timestep
  END DO
!$OMP END PARALLEL DO

END IF


! Item 202 Large scale snow

IF (icode <= 0 .AND. sf(202,4)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(202,4,im_index)),ls_snow,                        &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,202,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(ls_snow)"
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i)                                                              &
!$OMP SHARED(tdims,stashwork,si,im_index,timestep)
  DO i = 1, tdims%i_end*tdims%j_end
    stashwork(si(202,4,im_index)+i-1)=                                        &
         stashwork(si(202,4,im_index)+i-1)* timestep
  END DO
!$OMP END PARALLEL DO

END IF


! Item 203 Large scale rain

! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_ls_rain(:,:) = ls_rain(:,:)
END IF

IF (icode <= 0 .AND. (sf(203,4) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(lsrr)) THEN
    ALLOCATE ( lsrr(tdims%i_end*tdims%j_end) )
    lsrr(1) = rmdi
  END IF
END IF

IF (sf(203,4)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(203,4,im_index)),ls_rain,                        &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,203,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(ls_rain)"
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ji)                                                         &
!$OMP SHARED(tdims,lsrr,ls_rain)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ji = (j-1)*tdims%i_end+i
      lsrr(ji) = ls_rain(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF


! Item 204 Large scale snow
! this is also used by stash 20014 precip_symbol
IF (icode <=0 .AND. (sf(stashcode_pws_precip_sym, stashcode_pws_sec))) THEN
  pws_precip_sym_ls_snow(:,:) = ls_snow(:,:)
END IF

IF (icode <= 0 .AND. (sf(204,4) .OR. l_ac)) THEN

  IF (.NOT. ALLOCATED(lssr)) THEN
    ALLOCATE ( lssr(tdims%i_end*tdims%j_end) )
    lssr(1) = rmdi
  END IF
END IF

IF (sf(204,4)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(204,4,im_index)),ls_snow,                        &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,204,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(ls_snow)"
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ji)                                                         &
!$OMP SHARED(tdims,lsrr,ls_snow,lssr)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ji = (j-1)*tdims%i_end+i
      lssr(ji) = ls_snow(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF



! Items 231-236 mineral dust scavenged by LS precip

IF (l_dust) THEN

  IF ( sf(440,3) ) THEN
    ! The first component of total dust deposition flux is
    ! calculated below, so allocate array in preparation
    ALLOCATE(tot_dust_dep_flux(tdims%i_len,tdims%j_len))
    tot_dust_dep_flux(:,:) = 0.0
  END IF

  DO k = 1, ndiv

    ! We require this calculation below if either this item is
    ! required or total dust deposition flux (3,440) is required
    IF ( (icode <= 0) .AND. (sf(440,3) .OR. sf(230+k,4)) ) THEN

      ! Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,lscav_dust_all,timestep,k)
      DO j = 1, lspice_dim2
        DO i = 1, lspice_dim1
          lscav_dust_all(i, j, k)=lscav_dust_all(i, j, k)/timestep
        END DO
      END DO
!$OMP END PARALLEL DO

      IF ( sf(230+k,4) ) THEN

        ! DEPENDS ON: copydiag
        CALL copydiag(stashwork(si(230+k,4,im_index)),                        &
          lscav_dust_all(tdims%i_start:tdims%i_end,                           &
                         tdims%j_start:tdims%j_end,k),                        &
          tdims%i_len,tdims%j_len,0,0,0,0, at_extremity,                      &
          atmos_im,4,230+k,                                                   &
          icode,cmessage)

        IF (icode  >   0) THEN
          cmessage=": ERROR IN COPYDIAG(LSCAV_DUST_ALL)"
        END IF

      END IF

      ! Add dust deposition from largescale precip to total
      ! dust dep flux
      IF ( sf(440,3) ) THEN

        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            tot_dust_dep_flux(i,j) = tot_dust_dep_flux(i,j) +                 &
                                     lscav_dust_all(i,j,k)
          END DO
        END DO

      END IF

    END IF

  END DO !NDIV

END IF !L_DUST


! Item 215 NH3 scavenged by LS precip

IF (icode <= 0 .AND. sf(215,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,lscav_nh3,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_nh3(i, j)=lscav_nh3(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(215,4,im_index)),lscav_nh3,                      &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,215,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="microphysics_ctl  : error in copydiag(lscav_nh3)"
  END IF

END IF

! Item 216 SO2 scavenged by LS precip

IF (icode <= 0 .AND. sf(216,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,lscav_tr,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_tr(i, j)=lscav_tr(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(216,4,im_index)),lscav_tr,                       &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,216,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(lscav_tr)"
  END IF

END IF

! Item 219 Dissolved SO4 aerosol scavenged by LS precip

IF (icode <= 0 .AND. sf(219,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,rnout_tracer,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      rnout_tracer(i, j)=rnout_tracer(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(219,4,im_index)),rnout_tracer,                   &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,219,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(rnout_tracer)"
  END IF

END IF

! Item 220 Soot scavenged by LS rainout

IF (icode <= 0 .AND. sf(220,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,rnout_soot,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      rnout_soot(i, j)=rnout_soot(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(220,4,im_index)),rnout_soot,                     &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,220,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(rnout_soot)"
  END IF

END IF

! Item 221 Soot scavenged by LS washout

IF (icode <= 0 .AND. sf(221,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,lscav_soot,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_soot(i, j)=lscav_soot(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(221,4,im_index)),lscav_soot,                     &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,221,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(lscav_soot)"
  END IF

END IF

! Item 228 Fossil-fuel organic carbon scavenged by LS rainout

IF (icode <= 0 .AND. sf(228,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,rnout_ocff,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      rnout_ocff(i, j)=rnout_ocff(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(228,4,im_index)),rnout_ocff,                     &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,228,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(rnout_ocff)"
  END IF

END IF

! Item 229 Fossil-fuel organic carbon scavenged by LS washout

IF (icode <= 0 .AND. sf(229,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim2,lspice_dim1,lscav_ocff,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_ocff(i, j)=lscav_ocff(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(229,4,im_index)),lscav_ocff,                     &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,229,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(lscav_ocff)"
  END IF

END IF

!  Copy T to STASHwork

IF (icode <= 0 .AND. sf(004,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(004,4,im_index)),t,                           &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0,                           &
       at_extremity,                                                          &
       stlist(1,stindex(1,004,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,004,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(T)"
  END IF

END IF



!  Copy Cloud water to STASHwork

IF (icode <= 0 .AND. sf(205,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(205,4,im_index)),qcl,                         &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,205,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,205,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(cloud water)"
  END IF

END IF

IF (icode <= 0 .AND. sf(206,4)) THEN

  !  Copy Cloud ice to STASHwork

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(206,4,im_index)),qcf,                         &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,206,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,206,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(cloud ice)"
  END IF

END IF

! -----------------------------
! 207 relative humidity wrt ice (T<0degC) and water (T>0degC) (mdl levs)

item = 207
IF (icode <= 0 .AND. sf(item,sect)) THEN

  IF (l_new_qsat_mphys) THEN
    IF (l_mr_physics) THEN
      CALL qsat_mix_new(work_3d,t,p_layer_centres,                            &
                        tdims%i_end,tdims%j_end,tdims%k_end)
    ELSE
      CALL qsat_new(work_3d,t,p_layer_centres,                                &
                        tdims%i_end,tdims%j_end,tdims%k_end)
    END IF
  ELSE
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(work_3d,t,p_layer_centres(1,1,1),                           &
              tdims%i_end*tdims%j_end*tdims%k_end, l_mr_physics )
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,q)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.0
        !  Supersaturation (>100%) can occur with mixed phase scheme but
        !  negative humidity is removed from the diagnostic:
        IF (work_3d(i,j,k)  <   0.0) THEN
          work_3d(i,j,k) = 0.0
        END IF

      END DO  ! i
    END DO    ! j
  END DO      ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,             &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 207)"//cmessage
  END IF

END IF   ! item 207

! -------------------------------
! 208 relative humidity wrt water (on model levels)

item = 208
IF (icode  <=  0 .AND. sf(item,sect)) THEN
       ! q saturation is put in work_3d array

  IF (l_new_qsat_mphys) THEN
    CALL qsat_wat_new(work_3d,t,p_layer_centres,                              &
                      tdims%i_end,tdims%j_end,tdims%k_end )
  ELSE
    ! DEPENDS ON: qsat_wat
    CALL qsat_wat(work_3d,t,p_layer_centres(1,1,1),                           &
              tdims%i_end*tdims%j_end*tdims%k_end )
  END IF

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,q)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = q(i,j,k)/work_3d(i,j,k)*100.0
             !  Supersaturation wrt water is limited to =< 100%
        IF (work_3d(i,j,k) > 100.0) THEN
          work_3d(i,j,k) = 100.0
        END IF
             !  Negative humidity also removed from the diagnostic
        IF (work_3d(i,j,k) < 0.0) THEN
          work_3d(i,j,k) = 0.0
        END IF
      END DO  ! i
    END DO    ! j
  END DO      ! k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work_3d,                                                               &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,             &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d (item 208, rhw)"//cmessage
  END IF

END IF   ! item 208


!--------------------------------------------------------------------------
! Item 209: Surface Graupel Amount
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(209,4)) THEN

  ! Calculate the amount falling out in a single timestep

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,work_2d,timestep,ls_graup)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = ls_graup(i,j) * timestep
    END DO
  END DO
!$END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(209,4,im_index)),work_2d,                        &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,209,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(STASH 4/209: ls_graup)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 210: 3D droplet number (autoconversion-derived)
!--------------------------------------------------------------------------

! Copy droplet number  to STASHwork

IF (icode <= 0 .AND. sf(210,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(210,4,im_index)), n_drop_3d,                 &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,210,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,210,                                                        &
       icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(n_drop_3d)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 211: 3D droplet number (aerosol-derived)
!--------------------------------------------------------------------------

! Copy droplet number  to STASHwork

IF (icode <= 0 .AND. sf(211,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(211,4,im_index)),n_drop_tpr,                 &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,211,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,211,                                                        &
       icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(n_drop_tpr)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 212: Surface Graupel Rate
!--------------------------------------------------------------------------

IF (icode <= 0 .AND. sf(212,4)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(212,4,im_index)),ls_graup,                       &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,4,212,                                                        &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(STASH 4/212: ls_graup)"
  END IF

END IF


! Item 213 Ammonium nitrate scavenged by LS rainout
IF (icode <= 0 .AND. sf(213,4)) THEN
  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim1,lspice_dim2,rnout_nitrate,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      rnout_nitrate(i, j)=rnout_nitrate(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(213,4,im_index)),rnout_nitrate,                  &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,213,                                                       &
        icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag(rnout_nitrate)"
  END IF
END IF

! Item 214 Ammonium nitrate scavenged by LS washout
IF (icode <= 0 .AND. sf(214,4)) THEN
  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim1,lspice_dim2,lscav_nitrate,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_nitrate(i, j)=lscav_nitrate(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(214,4,im_index)),lscav_nitrate,                  &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,214,                                                       &
        icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag(lscav_nitrate)"
  END IF

END IF



!  Copy q  to STASHwork

IF (icode <= 0 .AND. sf(010,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(010,4,im_index)),q,                           &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,010,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,010,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(q)"
  END IF

END IF


! Copy ls_rain3d  to STASHwork

IF (icode <= 0 .AND. sf(222,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(222,4,im_index)),ls_rain3d,                   &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,222,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,222,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(ls_rain3d)"
  END IF

END IF


! Copy ls_snoww3d  to STASHwork

IF (icode <= 0 .AND. sf(223,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(223,4,im_index)),ls_snow3d,                   &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,223,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,223,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(ls_snow3d)"
  END IF

END IF


! Need to produce diagnostic 225 before 224 in order to save memory.

IF (icode <= 0 .AND. sf(225,4)) THEN

  !  Supercooled 3D rain content. It is equal to
  !  the 3D rainrate at T < 0 and equal to 0 at T > 0
  !  Alter the array LS_RAIN3D directly

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,t,ls_rain3d)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (t(i,j,k)  >=  zerodegc) THEN
          ! Warm temperatures
          ls_rain3d(i,j,k)=0.0
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Copy supercooled rain (now in ls_rain3d)  to STASHwork

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(225,4,im_index)),ls_rain3d,                   &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,225,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,225,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(supercooled rain)"
  END IF

END IF

IF (icode <= 0 .AND. sf(224,4)) THEN

  !  Supercooled liquid water content at TIME LEVEL N. It is equal to
  !  the liquid water content at T < 0 and equal to 0 at T > 0
  !  Use LS_RAIN3D as the array in order to save memory

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,t,ls_rain3d,qcl_n)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (t(i,j,k)  <   zerodegc) THEN
          ! Supercooled temperatures
          ! Use time level n fields in this diagnostic
          ls_rain3d(i,j,k)=qcl_n(i,j,k)
        ELSE
          ! Warm temperatures
          ls_rain3d(i,j,k)=0.0
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Copy supercooled liquid (now in ls_rain3d)  to STASHwork

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(224,4,im_index)),ls_rain3d,                   &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,224,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,224,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(supercooled liq)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 226: 3D Graupel Fall Rate
!--------------------------------------------------------------------------

! Copy droplet number  to STASHwork

IF (icode <= 0 .AND. sf(226,4)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork(si(226,4,im_index)), ls_graup3d,                &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,226,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,226,                                                        &
       icode,cmessage )

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(n_drop_tpr)"
  END IF

END IF

IF (icode <= 0 .AND. sf(227,4)) THEN

  ! Copy rain fraction to stashwork

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(stashwork(si(227,4,im_index)),rainfrac3d,                  &
       tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0,                         &
       at_extremity,                                                          &
       stlist(1,stindex(1,227,4,im_index)),len_stlist,                        &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,4,227,                                                        &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=":error in copydiag_3d(rain fraction)"
  END IF

END IF

! Item 237 Biomass scavenged by LS rainout

IF (icode <= 0 .AND. sf(237,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim1,lspice_dim2,rnout_bmass,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      rnout_bmass(i, j)=rnout_bmass(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(237,4,im_index)),rnout_bmass,                    &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,237,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(rnout_bmass)"
  END IF

END IF

! Item 238 Biomass scavenged by LS washout

IF (icode <= 0 .AND. sf(238,4)) THEN

  !       Convert "per timestep" diagnostic to "per second":

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(lspice_dim1,lspice_dim2,lscav_bmass,timestep)
  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_bmass(i, j)=lscav_bmass(i, j)/timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(238,4,im_index)),lscav_bmass,                    &
        tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                        &
        atmos_im,4,238,                                                       &
        icode,cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(lscav_bmass)"
  END IF

END IF


      !---------------------------------------------------------------
      ! Homogeneous nucleation
      !---------------------------------------------------------------
item = 240

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pifrw,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Heterogeneous nucleation
      !---------------------------------------------------------------
item = 241

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    piprm,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Deposition of ice
      !---------------------------------------------------------------
item = 243

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pidep,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Deposition of snow aggregates
      !---------------------------------------------------------------
item = 245

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psdep,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF


      !---------------------------------------------------------------
      ! Ice collection of cloud liquid water (riming)
      !---------------------------------------------------------------
item = 247

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    piacw,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Snow collection of cloud liquid water (riming)
      !---------------------------------------------------------------
item = 248

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psacw,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Ice collection of rain (capture)
      !---------------------------------------------------------------
item = 249

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    piacr,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Snow collection of rain (capture)
      !---------------------------------------------------------------
item = 250

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psacr,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Evaporation of melting ice
      !---------------------------------------------------------------
item = 251

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pimltevp,                                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Evaporation of melting snow
      !---------------------------------------------------------------
item = 252

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psmltevp,                                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Melting ice
      !---------------------------------------------------------------
item = 253

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pimlt,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Melting snow
      !---------------------------------------------------------------
item = 254

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psmlt,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Snow aggregate autoconversion
      !---------------------------------------------------------------
item = 255

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psaut,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Snow collection of ice (capture)
      !---------------------------------------------------------------
item = 256

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psaci,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Rain autoconversion
      !---------------------------------------------------------------
item = 257

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    praut,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Rain collection of cloud liquid water (accretion)
      !---------------------------------------------------------------
item = 258

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pracw,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Evaporation of rain
      !---------------------------------------------------------------
item = 259

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    prevp,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Graupel autoconversion
      !---------------------------------------------------------------
item = 260

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pgaut,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Graupel collection of cloud liquid water (accretion)
      !---------------------------------------------------------------
item = 261

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pgacw,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Graupel collection of snow (capture)
      !---------------------------------------------------------------
item =262

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pgacs,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Melting graupel
      !---------------------------------------------------------------
item = 263

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pgmlt,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Ice crystal sedimentation
      !---------------------------------------------------------------
item = 265

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pifall,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Snow sedimentation
      !---------------------------------------------------------------
item = 266

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    psfall,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Rain sedimentation
      !---------------------------------------------------------------
item = 267

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    prfall,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Graupel sedimentation
      !---------------------------------------------------------------
item = 268

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pgfall,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Droplet settling of liquid
      !---------------------------------------------------------------
item = 269

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    plset,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Evaporated settled droplets
      !---------------------------------------------------------------
item = 270

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    plevpset,                                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,*)                                                         &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Homogeneous nucleation of rain
      !---------------------------------------------------------------
item = 271

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    pifrr,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Surface Maximum Hail Size (mm)
      !---------------------------------------------------------------
item = 275

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! In this case, work out the surface graupel water content (GWC) in
  ! kg m-3 by multiplying the air density by the graupel mass mixing ratio
  ! Pass this to max_hail_size_2d to generate the maximum hail size

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,work_2d,qgraup,rho_work)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = qgraup(i,j,1) * rho_work(i,j,1)
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL max_hail_size_2d(work_2d, surface_hail_diam)

! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), surface_hail_diam,         &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF

END IF

      !---------------------------------------------------------------
      ! Maximum Hail Size in column (mm)
      !---------------------------------------------------------------
item = 276

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! In this case, work out the maximum graupel water content (GWC) in
  ! the column (kg m-3) by multiplying the air density by the graupel
  ! mass mixing ratio. Pass this to max_hail_size_2d to generate the
  ! maximum hail size

  gwc_max(:,:) = 0.0
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = qgraup(i,j,k) * rho_work(i,j,k)
        gwc_max(i,j)   = MAX( gwc_max(i,j), work_3d(i,j,k) )
      END DO
    END DO
  END DO

  CALL max_hail_size_2d(gwc_max, max_col_hail_diam)

! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), max_col_hail_diam,         &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF

END IF

      !---------------------------------------------------------------
      ! Maximum Hail Size on model levels (mm)
      !---------------------------------------------------------------
item = 277

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! Re-use work_3d from above to save calculations
  IF (.NOT. sf(276, sect)) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,qgraup,rho_work)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work_3d(i,j,k) = qgraup(i,j,k) * rho_work(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF ! Not sf(276, 4)

  CALL max_hail_size_3d(work_3d, max_hail_diam)

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    max_hail_diam,                                                            &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF

END IF

      !---------------------------------------------------------------
      ! qcl increment by mixed phase scheme
      !---------------------------------------------------------------
item = 294

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    qcl_mpt,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! cfl diagnosed by mixed phase scheme
      !---------------------------------------------------------------
item = 295

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dcfl_mp,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Turbulent decorrelation timescale
      !---------------------------------------------------------------
item = 296

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    tau_d,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Inverse phase-relaxation timescale
      !---------------------------------------------------------------
item = 297

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    inv_prt,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Diagnosed turbulent dissipation rate
      !---------------------------------------------------------------
item = 298

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    disprate,                                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Inverse mixing timescale
      !---------------------------------------------------------------
item = 299

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    inv_mt,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Mean of subgrid supersaturation PDF (Turbulence)
      !---------------------------------------------------------------
item = 300

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    si_avg,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! Variance of subgrid supersaturation PDF (Turbulence)
      !---------------------------------------------------------------
item = 301

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    sigma2_s,                                                                 &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

!---------------------------------------------------------------
! 302: Surface Snowfall rate without graupel (2D)
!---------------------------------------------------------------
item = 302

IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,work_2d,ls_snow,ls_graup,timestep)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = (ls_snow(i,j) - ls_graup(i,j)) * timestep
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! qcl increment from mixed phase scheme (PC2)
      !---------------------------------------------------------------
item = 303

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    qcl_incr_mp_pc2,                                                          &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

!---------------------------------------------------------------
! 304: Surface Snowfall rate without graupel (2D)
!---------------------------------------------------------------
item = 304

IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,work_2d,ls_snow,ls_graup)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      work_2d(i,j) = ls_snow(i,j) - ls_graup(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), work_2d,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

      !---------------------------------------------------------------
      ! cfl increment from mixed phase scheme (PC2)
      !---------------------------------------------------------------
item = 313

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    cfl_incr_mp,                                                              &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

!---------------------------------------------------------------
! cfl increment from mixed phase scheme (PC2)
!---------------------------------------------------------------
item = 323

IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_3d,ls_snow3d,ls_graup3d)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_3d(i,j,k) = ls_snow3d(i,j,k) - ls_graup3d(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    work_3d,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF


!---------------------------------------------------------------
! Seeder Feeder sub-grid orographic water mixing ratio
!---------------------------------------------------------------
item = 400

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    sfwater,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

!---------------------------------------------------------------
! Seeder Feeder sub-grid orographic rain increase rate
!---------------------------------------------------------------
item = 401

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    sfrain,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF

!---------------------------------------------------------------
! Seeder Feeder sub-grid orographic snow increase rate
!---------------------------------------------------------------
item = 402

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    sfsnow,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF


      !---------------------------------------------------------------
      ! Total Reflectivity
      !---------------------------------------------------------------
item = 118

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_tot,                                                                  &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Reflectivity due to graupel
      !---------------------------------------------------------------
item = 113

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_g,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Reflectivity due to ice aggregates
      !---------------------------------------------------------------
item = 114

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_i,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Reflectivity due to ice crystals
      !---------------------------------------------------------------
item = 115

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_i2,                                                                   &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Reflectivity due to rain
      !---------------------------------------------------------------
item = 116

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_r,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF
      !---------------------------------------------------------------
      ! Reflectivity due to liquid cloud
      !---------------------------------------------------------------
item = 117

IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (                                                          &
    stashwork(si(item,sect,im_index)),                                        &
    dbz_l,                                                                    &
    tdims%i_end, tdims%j_end, tdims%k_end,0,0,0,0, at_extremity,              &
    stlist(1,stindex(1,item,sect,im_index)),                                  &
    len_stlist, stash_levels,num_stash_levels+1,                              &
    atmos_im,sect,item,                                                       &
    icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag_3d for diagnostic ',                                 &
      'section 4, item ',item
  END IF
END IF


      !---------------------------------------------------------------
      ! Surface (model level 1 reflectivity)
      !---------------------------------------------------------------
item = 110

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! First calculate dbz_surf: surface reflectivity

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(tdims,dbz_surf,dbz_tot)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dbz_surf(i,j) = dbz_tot(i,j,1)
    END DO ! i
  END DO ! j
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), dbz_surf,                  &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag for diagnostic ',                                    &
      'section 4, item ',item
  END IF
END IF ! icode <= 0

      !---------------------------------------------------------------
      ! Max reflectivity in the column
      !---------------------------------------------------------------
item = 111

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! Calculate Maximum Reflectivity

  dbz_max(:,:) = ref_lim

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (dbz_tot(i,j,k) > dbz_max(i,j)) THEN
          dbz_max(i,j) = dbz_tot(i,j,k)
        END IF
      END DO ! i
    END DO ! j
  END DO ! k

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), dbz_max,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag for diagnostic ',                                    &
      'section 4, item ',item
  END IF
END IF ! icode <= 0

      !---------------------------------------------------------------
      ! Reflectivity at 1km AGL
      !---------------------------------------------------------------
item = 112

IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! First calculate reflectivity at 1 km AGL using the layer thickness

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      z_agl(i,j)   = 0.0
      dbz_1km(i,j) = ref_lim
      dbz_def(i,j) = .TRUE.
    END DO
  END DO
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        z_agl(i,j) = z_agl(i,j) + deltaz(i,j,k)

        IF (dbz_def(i,j) .AND. z_agl(i,j) >= agl_1km ) THEN

          dbz_1km(i,j) = dbz_tot(i,j,k)
          dbz_def(i,j) = .FALSE. ! Prevent overwriting by higher levels

        END IF ! dbz_def

      END DO ! i
    END DO ! j
  END DO !k

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)), dbz_1km,                   &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                         &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode  >   0) THEN
    WRITE(cmessage,'(A,A,I3)')                                                &
      'ERROR in copydiag for diagnostic ',                                    &
      'section 4, item ',item
  END IF
END IF ! icode <= 0


! Single point exception handling
IF (icode /= 0) THEN
  CALL ereport(RoutineName,icode,cmessage)
END IF

! Warning messages where user includes graupel output at same time
! as snowfall output

IF (l_mcr_qgraup) THEN
  ! Don't bother warning if run doesn't include graupel
  ! as it may just confuse a user who doesn't have prognostic graupel
  ! turned on

  IF ( sf(209, 4) .AND. sf(202, 4) .AND. l_warn_gr_am ) THEN
    ! User has chosen to output both graupel amount and snow amount
    ! Add a warning to try and prevent them being added together
    ! in any post-run analysis.
    icode         = -1
    l_warn_gr_am  = .FALSE. ! Only call this warning once
    cmessage      =                                                           &
    'Model run outputs STASH 4/202 and 4/209 on same timestep. '//newline//   &
    'Please note that surface snowfall amount (STASH diagnostic'//newline//   &
    '4/202) already includes contributions from graupel        '//newline//   &
    '(STASH diagnostic 4/209). If you sum these diagnostics    '//newline//   &
    'together when analysing the model run, your answers will  '//newline//   &
    'be wrong!'

    CALL ereport(RoutineName,icode,cmessage)
  END IF ! sf 209/4

  IF ( sf(212, 4) .AND. sf(204, 4) .AND. l_warn_gr_ra ) THEN
    ! User has chosen to output both graupel rate and snow rate
    ! Add a warning to try and prevent them being added together
    ! in any post-run analysis.
    icode        = -1
    l_warn_gr_ra = .FALSE. ! Only call this warning once
    cmessage     =                                                            &
    'Model run outputs STASH 4/204 and 4/212 on same timestep. '//newline//   &
    'Please note that surface snowfall rate (STASH diagnostic  '//newline//   &
    '4/204) already includes contributions from graupel        '//newline//   &
    '(STASH diagnostic 4/212). If you sum these diagnostics    '//newline//   &
    'together when analysing the model run, your answers will  '//newline//   &
    'be wrong!'

    CALL ereport(RoutineName,icode,cmessage)
  END IF

  IF ( sf(226, 4) .AND. sf(223, 4) .AND. l_warn_gr_ra3 ) THEN
    ! User has chosen to output both graupel rate (3D) and snow rate (3D)
    ! Add a warning to try and prevent them being added together in any
    ! post-run analysis.
    icode          = -1
    l_warn_gr_ra3  = .FALSE. ! Only call this warning once
    cmessage       =                                                          &
    'Model run outputs STASH 4/223 and 4/226 on same timestep. '//newline//   &
    'Please note that surface snowfall rate (STASH diagnostic  '//newline//   &
    '4/223) already includes contributions from graupel        '//newline//   &
    '(STASH diagnostic 4/226). If you sum these diagnostics    '//newline//   &
    'together when analysing the model run, your answers will  '//newline//   &
    'be wrong!'

    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF ! l_mcr_qgraup

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_lsrain
END MODULE diagnostics_lsrain_mod
