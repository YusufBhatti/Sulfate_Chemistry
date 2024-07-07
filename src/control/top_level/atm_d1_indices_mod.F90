! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

MODULE atm_d1_indices_mod

! Purpose: Contains indices used in the d1 array. At present, the 
!          d1 array is a single array containing many fields, these 
!          indices are set to distinguish between the fields and enable 
!          one to use pointers to make them more human readable.

USE atm_fields_bounds_mod, ONLY: pdims, tdims, udims, vdims, wdims,        &
                                 o3dims2, tdims_s, wdims_s
USE nlsizes_namelist_mod,  ONLY: n_cca_lev, st_levels, sm_levels,          &
                                 model_levels, tpps_ozone_levels, tr_vars, &
                                 tr_ukca, tr_lbc_ukca, tr_lbc_vars
USE clmchfcg_scenario_mod, ONLY: nsulpat

USE parkind1,     ONLY: jprb, jpim
USE yomhook,      ONLY: lhook, dr_hook
USE ereport_mod,  ONLY: ereport
USE umPrintMgr,   ONLY: umMessage, newline

IMPLICIT NONE

INTEGER :: jsst                     ! SST in idealised tests
INTEGER :: jexnersurf               ! surface exner
INTEGER, ALLOCATABLE :: jdryrho (:) ! Dry Density
INTEGER, ALLOCATABLE :: jetadot (:) ! Etadot
INTEGER, ALLOCATABLE :: jthetav (:) ! Potential virtual dry temperature
INTEGER :: jpsiws                   ! W component of Psi at surface (ENDGame 
                                    ! dynamics variable)
INTEGER :: jpsiwl                   ! W component of Psi at lid of model 
                                    ! (ENDGame dynamics variable)
INTEGER, ALLOCATABLE :: jmv     (:) ! humidity mixing ratio
INTEGER, ALLOCATABLE :: jmcl    (:) ! humidity mixing ratio
INTEGER, ALLOCATABLE :: jmcf    (:) ! humidity mixing ratio
INTEGER, ALLOCATABLE :: jmcf2   (:) ! humidity mixing ratio
INTEGER, ALLOCATABLE :: jmrain  (:) ! humidity mixing ratio
INTEGER, ALLOCATABLE :: jmgraup (:) ! humidity mixing ratio


INTEGER, ALLOCATABLE :: ju      (:) ! u component of wind
INTEGER, ALLOCATABLE :: jv      (:) ! v component of wind
INTEGER, ALLOCATABLE :: jw      (:) ! w component of wind
INTEGER, ALLOCATABLE :: jrho    (:) ! Density
INTEGER, ALLOCATABLE :: jtheta  (:) ! Potential temperature
INTEGER, ALLOCATABLE :: jq      (:) ! Specific humidity
INTEGER, ALLOCATABLE :: jqcl    (:) ! qcl
INTEGER, ALLOCATABLE :: jqcf    (:) ! qcf
INTEGER, ALLOCATABLE :: jqcf2   (:) ! second ice
INTEGER, ALLOCATABLE :: jqrain  (:) ! rain
INTEGER, ALLOCATABLE :: jqgraup (:) ! graupel

! CASIM cloud and ice microphysics prognostics 
INTEGER, ALLOCATABLE :: jcloudnumber (:) ! cloud number
INTEGER, ALLOCATABLE :: jrainnumber  (:) ! rain number
INTEGER, ALLOCATABLE :: jrain3mom    (:) ! third moment for rain
INTEGER, ALLOCATABLE :: jicenumber   (:) ! ice number
INTEGER, ALLOCATABLE :: jsnownumber  (:) ! snow number
INTEGER, ALLOCATABLE :: jsnow3mom    (:) ! third moment for snow 
INTEGER, ALLOCATABLE :: jgraupnumber (:) ! graupel number
INTEGER, ALLOCATABLE :: jgraup3mom   (:) ! graupel number

! CASIM activated aerosol prognostics
INTEGER, ALLOCATABLE :: jactivesolliquid   (:)
! activated soluble aerosol in liquid
INTEGER, ALLOCATABLE :: jactivesolrain     (:) 
! activated soluble aerosol in rain
INTEGER, ALLOCATABLE :: jactiveinsolice    (:) 
! activated insoluble aerosol in ice      
INTEGER, ALLOCATABLE :: jactivesolice      (:) 
! activated soluble aerosol in ice      
INTEGER, ALLOCATABLE :: jactiveinsolliquid (:) 
! activated insoluble aerosol in liquid
INTEGER, ALLOCATABLE :: jactivesolnumber   (:) 
! activated soluble aerosol number
INTEGER, ALLOCATABLE :: jactiveinsolnumber (:) 
! activated insoluble aerosol number

INTEGER, ALLOCATABLE :: je_trb  (:) ! Turbulent Kinetic Energy (TKE)
INTEGER, ALLOCATABLE :: jtsq_trb(:) ! Self covariance of thetal'
INTEGER, ALLOCATABLE :: jqsq_trb(:) ! Self coveriance of qw'
INTEGER, ALLOCATABLE :: jcov_trb(:) ! Correlation between thetal' and qw'
INTEGER :: jzhpar_shcu              ! Height of mixed layer used to evaluate the
                                    ! non-grad buoy flux

!PV-tracers
INTEGER, ALLOCATABLE :: jdPV_rad(:)
INTEGER, ALLOCATABLE :: jdPV_sw(:)
INTEGER, ALLOCATABLE :: jdPV_lw(:)
INTEGER, ALLOCATABLE :: jdPV_mic(:)
INTEGER, ALLOCATABLE :: jdPV_gwd(:)
INTEGER, ALLOCATABLE :: jdPV_ph1(:)
INTEGER, ALLOCATABLE :: jdPV_conv(:)
INTEGER, ALLOCATABLE :: jdPV_bl(:)
INTEGER, ALLOCATABLE :: jdPV_stph(:)
INTEGER, ALLOCATABLE :: jdPV_cld(:)
INTEGER, ALLOCATABLE :: jdPV_iau(:)
INTEGER, ALLOCATABLE :: jdPV_nud(:)
INTEGER, ALLOCATABLE :: jdPV_tot(:)
INTEGER, ALLOCATABLE :: jdPV_adv(:)
INTEGER, ALLOCATABLE :: jdPV_sol(:)
INTEGER, ALLOCATABLE :: jdPV_mass(:)
INTEGER, ALLOCATABLE :: jadv_only_PV(:)

!Diabatic tracers
INTEGER, ALLOCATABLE :: jdtheta_0(:)
INTEGER, ALLOCATABLE :: jdtheta_bl(:)
INTEGER, ALLOCATABLE :: jdtheta_bl_mix(:)
INTEGER, ALLOCATABLE :: jdtheta_bl_LH(:)
INTEGER, ALLOCATABLE :: jdtheta_conv(:)
INTEGER, ALLOCATABLE :: jdtheta_mic(:)
INTEGER, ALLOCATABLE :: jdtheta_rad(:)
INTEGER, ALLOCATABLE :: jdtheta_SW(:)
INTEGER, ALLOCATABLE :: jdtheta_LW(:)
INTEGER, ALLOCATABLE :: jdtheta_slow(:)
INTEGER, ALLOCATABLE :: jdtheta_cld(:)

INTEGER, ALLOCATABLE :: jexner_rho_levels(:) ! Exner pressure on rho levels

! Convection prognostics
INTEGER, ALLOCATABLE :: jconv_prog_1     (:)
INTEGER, ALLOCATABLE :: jconv_prog_2     (:)
INTEGER, ALLOCATABLE :: jconv_prog_3     (:)
INTEGER, ALLOCATABLE :: jconv_prog_precip(:)

! Total precipitation rate at start of timestep
INTEGER :: jtotalppn

! Stochastic physics prognostics for BL perturbations
INTEGER :: jbl_pert_rand_fld ! random field for BL pert in stochastic physics
INTEGER :: jbl_pert_flag     ! flag for BL pert in stochastic physics

! 1.2: Data variables stored in secondary space.
INTEGER, ALLOCATABLE :: jp(:)  ! Pressure on rho le
INTEGER, ALLOCATABLE :: jp_theta_levels(:) ! Pressure on theta levels
INTEGER, ALLOCATABLE :: jexner_theta_levels(:) ! Exner pressure on theta levels

! 1.3: Cloud Fields
! n_cca_lev is set in dervsize
INTEGER, ALLOCATABLE :: jccw_rad  (:) ! CCW profile to radiation
INTEGER, ALLOCATABLE :: jcca      (:) ! Convective cloud amount
INTEGER, ALLOCATABLE :: jcca_dp   (:) ! Deep Convective cloud amount
INTEGER, ALLOCATABLE :: jcca_md   (:) ! Mid-level Convective cloud amount
INTEGER, ALLOCATABLE :: jcca_sh   (:) ! Shallow Convective cloud amount
INTEGER, ALLOCATABLE :: jcf_area  (:) ! Area Cloud Fraction
INTEGER, ALLOCATABLE :: jcf_bulk  (:) ! Bulk Cloud Fraction
INTEGER, ALLOCATABLE :: jcf_liquid(:) ! Liquid cloud fraction
INTEGER, ALLOCATABLE :: jcf_frozen(:) ! Frozen cloud fraction

! 1.4: Soil Ancillary fields
INTEGER, ALLOCATABLE :: j_deep_soil_temp(:) ! Deep soil temperature
INTEGER, ALLOCATABLE :: jsmcl           (:) ! Soil moisture content in layers
INTEGER, ALLOCATABLE :: jsthu           (:) ! Unfrozen soil moisture fraction
INTEGER, ALLOCATABLE :: jsthf           (:) ! Frozen soil moisture fraction

! 1.5: Radiation Increments
INTEGER, ALLOCATABLE :: jsw_incs(:) ! SW radiation increments
INTEGER, ALLOCATABLE :: jlw_incs(:) ! LW radiation increments
INTEGER :: jdirpar ! PAR radiation increment

! 1.6: Ozone
INTEGER, ALLOCATABLE :: jozone(:)     ! ozone
INTEGER, ALLOCATABLE :: jtppsozone(:) ! tropopause-based ozone

! 1.7: Tracer and aerosol fields
INTEGER, ALLOCATABLE :: jtracer (:,:)  ! Tracers
INTEGER, ALLOCATABLE :: jtr_ukca(:,:)  ! UKCA Tracers

INTEGER, ALLOCATABLE :: jmurk_source(:) ! multi-level murk source
INTEGER, ALLOCATABLE :: jmurk       (:) ! multi-level murk content

! Lightning potential and charge tracers
INTEGER, ALLOCATABLE :: jflash_pot (:) ! Flash potential

! Turbulent qcl production 
INTEGER, ALLOCATABLE :: jbl_w_var  (:) ! BL w-variance 

INTEGER, ALLOCATABLE :: jdust_div1 (:) ! dust mmr, division 1
INTEGER, ALLOCATABLE :: jdust_div2 (:) ! dust mmr, division 2
INTEGER, ALLOCATABLE :: jdust_div3 (:) ! dust mmr, division 3
INTEGER, ALLOCATABLE :: jdust_div4 (:) ! dust mmr, division 4
INTEGER, ALLOCATABLE :: jdust_div5 (:) ! dust mmr, division 5
INTEGER, ALLOCATABLE :: jdust_div6 (:) ! dust mmr, division 6

INTEGER, ALLOCATABLE :: jso2       (:) ! sulphur dioxide gas
INTEGER, ALLOCATABLE :: jdms       (:) ! dimethyl sulphide gas
INTEGER, ALLOCATABLE :: jso4_aitken(:) ! Aitken mode sulphate aer
INTEGER, ALLOCATABLE :: jso4_accu  (:) ! accumulation mode sulpha
INTEGER, ALLOCATABLE :: jso4_diss  (:) ! dissloved  sulphate aero
INTEGER, ALLOCATABLE :: jh2o2      (:) ! hydrogen peroxide mmr
INTEGER, ALLOCATABLE :: jnh3       (:) ! ammonia gas mmr

INTEGER, ALLOCATABLE :: jsoot_new  (:) ! fresh soot mmr
INTEGER, ALLOCATABLE :: jsoot_agd  (:) ! aged soot mmr
INTEGER, ALLOCATABLE :: jsoot_cld  (:) ! soot in cloud mmr

INTEGER, ALLOCATABLE :: jbmass_new (:) ! fresh biomass mmr
INTEGER, ALLOCATABLE :: jbmass_agd (:) ! aged biomass mmr
INTEGER, ALLOCATABLE :: jbmass_cld (:) ! cloud biomass mmr

INTEGER, ALLOCATABLE :: jocff_new  (:) ! fresh OCFF mmr
INTEGER, ALLOCATABLE :: jocff_agd  (:) ! aged OCFF mmr
INTEGER, ALLOCATABLE :: jocff_cld  (:) ! OCFF in cloud mmr

INTEGER, ALLOCATABLE :: jso2_natem (:) ! natural SO2 emissions
INTEGER, ALLOCATABLE :: joh        (:) ! hydroxyl radical ancilla
INTEGER, ALLOCATABLE :: jho2       (:) ! hydrogen dioxide ancilla
INTEGER, ALLOCATABLE :: jh2o2_limit(:) ! limiting H2O2 ancillary
INTEGER, ALLOCATABLE :: jo3_chem   (:) ! ozone for chemistry anci

INTEGER, ALLOCATABLE :: jco2       (:) ! 3D CO2 FIELD

INTEGER, ALLOCATABLE :: joh_ukca   (:) ! OH MMR from UKCA
INTEGER, ALLOCATABLE :: jho2_ukca  (:) ! HO2 MMR from UKCA
INTEGER, ALLOCATABLE :: jh2o2_ukca (:) ! H2O2 MMR from UKCA
INTEGER, ALLOCATABLE :: jo3_ukca   (:) ! O3 MMR from UKCA
INTEGER, ALLOCATABLE :: jhno3_ukca (:) ! HNO3 MMR from UKCA

INTEGER, ALLOCATABLE :: jozone_tracer(:) ! Prognostic O3 Tracer(Cariol)
INTEGER, ALLOCATABLE :: jo3_prod_loss(:) ! Cariol O3 Prod-Loss (P-L)
INTEGER, ALLOCATABLE :: jo3_p_l_vmr  (:) ! Cariol O3 P-L wrt VMR
INTEGER, ALLOCATABLE :: jo3_vmr      (:) ! Cariol O3 Vol Mix Ratio-VMR
INTEGER, ALLOCATABLE :: jo3_p_l_temp (:) ! Cariol O3 P-L wrt temp
INTEGER, ALLOCATABLE :: jo3_temp     (:) ! Cariol O3 temp
INTEGER, ALLOCATABLE :: jo3_p_l_colo3(:) ! Cariol O3 P-L wrt colO3
INTEGER, ALLOCATABLE :: jo3_colo3    (:) ! Cariol O3 column (colO3)

INTEGER, ALLOCATABLE :: jarclbiog_bg (:) ! Biogenic aerosol climatology
INTEGER, ALLOCATABLE :: jarclbiom_fr (:) ! Biomass burning (fresh) aerosol clim
INTEGER, ALLOCATABLE :: jarclbiom_ag (:) ! Biomass burning (aged) aerosol clim
INTEGER, ALLOCATABLE :: jarclbiom_ic (:) ! Biomass burning (in-cloud) 
                                         ! aerosol clim

INTEGER, ALLOCATABLE :: jarclblck_fr (:) ! Black carbon (fresh) aerosol clim
INTEGER, ALLOCATABLE :: jarclblck_ag (:) ! Black carbon (aged) aerosol clim
INTEGER, ALLOCATABLE :: jarclsslt_fi (:) ! Sea salt (film mode) aerosol clim
INTEGER, ALLOCATABLE :: jarclsslt_jt (:) ! Sea salt (jet mode) aerosol clim

INTEGER, ALLOCATABLE :: jarclsulp_ac (:) ! Sulphate (accumulation mode) 
                                         ! aero clim
INTEGER, ALLOCATABLE :: jarclsulp_ak (:) ! Sulphate (Aitken mode) aerosol clim
INTEGER, ALLOCATABLE :: jarclsulp_di (:) ! Sulphate (dissolved) aerosol clim

INTEGER, ALLOCATABLE :: jarcldust_b1 (:) ! Dust (bin 1) aerosol climatology
INTEGER, ALLOCATABLE :: jarcldust_b2 (:) ! Dust (bin 2) aerosol climatology
INTEGER, ALLOCATABLE :: jarcldust_b3 (:) ! Dust (bin 3) aerosol climatology
INTEGER, ALLOCATABLE :: jarcldust_b4 (:) ! Dust (bin 4) aerosol climatology
INTEGER, ALLOCATABLE :: jarcldust_b5 (:) ! Dust (bin 5) aerosol climatology
INTEGER, ALLOCATABLE :: jarcldust_b6 (:) ! Dust (bin 6) aerosol climatology

INTEGER, ALLOCATABLE :: jarclocff_fr (:) ! Org carbon fossil fuel (fresh) 
                                         ! aero clim
INTEGER, ALLOCATABLE :: jarclocff_ag (:) ! Org carbon fossil fuel (aged) 
                                         ! aero clim
INTEGER, ALLOCATABLE :: jarclocff_ic (:) ! Org carbon fossil fuel (in-cloud) 
                                         ! aero clim
INTEGER, ALLOCATABLE :: jarcldlta_dl (:) ! Delta aerosol climatology
INTEGER, ALLOCATABLE :: jnitr_acc    (:) ! Accumulation nitrate aerosol
INTEGER, ALLOCATABLE :: jnitr_diss   (:) ! Dissolved nitrate aerosol

! Section 54 GLOMAP_CLIM aerosol climatology

INTEGER, ALLOCATABLE :: jgc_nd_nuc_sol(:) ! GLOMAP_CLIM NUC (Sol) number density
INTEGER, ALLOCATABLE :: jgc_nuc_sol_su(:) ! GLOMAP_CLIM NUC (Sol) SO4
INTEGER, ALLOCATABLE :: jgc_nuc_sol_oc(:) ! GLOMAP_CLIM NUC (Sol) OC

INTEGER, ALLOCATABLE :: jgc_nd_ait_sol(:) ! GLOMAP_CLIM AIT (Sol) number density
INTEGER, ALLOCATABLE :: jgc_ait_sol_su(:) ! GLOMAP_CLIM AIT (Sol) SO4
INTEGER, ALLOCATABLE :: jgc_ait_sol_bc(:) ! GLOMAP_CLIM AIT (Sol) BC
INTEGER, ALLOCATABLE :: jgc_ait_sol_oc(:) ! GLOMAP_CLIM AIT (Sol) OC

INTEGER, ALLOCATABLE :: jgc_nd_acc_sol(:) ! GLOMAP_CLIM ACC (Sol) number density
INTEGER, ALLOCATABLE :: jgc_acc_sol_su(:) ! GLOMAP_CLIM ACC (Sol) SO4
INTEGER, ALLOCATABLE :: jgc_acc_sol_bc(:) ! GLOMAP_CLIM ACC (Sol) BC
INTEGER, ALLOCATABLE :: jgc_acc_sol_oc(:) ! GLOMAP_CLIM ACC (Sol) OC
INTEGER, ALLOCATABLE :: jgc_acc_sol_ss(:) ! GLOMAP_CLIM ACC (Sol) SS

INTEGER, ALLOCATABLE :: jgc_nd_cor_sol(:) ! GLOMAP_CLIM COR (Sol) number density
INTEGER, ALLOCATABLE :: jgc_cor_sol_su(:) ! GLOMAP_CLIM COR (Sol) SO4
INTEGER, ALLOCATABLE :: jgc_cor_sol_bc(:) ! GLOMAP_CLIM COR (Sol) BC
INTEGER, ALLOCATABLE :: jgc_cor_sol_oc(:) ! GLOMAP_CLIM COR (Sol) OC
INTEGER, ALLOCATABLE :: jgc_cor_sol_ss(:) ! GLOMAP_CLIM COR (Sol) SS

INTEGER, ALLOCATABLE :: jgc_nd_ait_ins(:) ! GLOMAP_CLIM AIT (Ins) number density
INTEGER, ALLOCATABLE :: jgc_ait_ins_bc(:) ! GLOMAP_CLIM AIT (Ins) BC
INTEGER, ALLOCATABLE :: jgc_ait_ins_oc(:) ! GLOMAP_CLIM AIT (Ins) OC

! 1.8: Multi-level user ancillary fields
INTEGER, ALLOCATABLE :: juser_mult1(:)
INTEGER, ALLOCATABLE :: juser_mult2(:)
INTEGER, ALLOCATABLE :: juser_mult3(:)
INTEGER, ALLOCATABLE :: juser_mult4(:)
INTEGER, ALLOCATABLE :: juser_mult5(:)
INTEGER, ALLOCATABLE :: juser_mult6(:)
INTEGER, ALLOCATABLE :: juser_mult7(:)
INTEGER, ALLOCATABLE :: juser_mult8(:)
INTEGER, ALLOCATABLE :: juser_mult9(:)
INTEGER, ALLOCATABLE :: juser_mult10(:)
INTEGER, ALLOCATABLE :: juser_mult11(:)
INTEGER, ALLOCATABLE :: juser_mult12(:)
INTEGER, ALLOCATABLE :: juser_mult13(:)
INTEGER, ALLOCATABLE :: juser_mult14(:)
INTEGER, ALLOCATABLE :: juser_mult15(:)
INTEGER, ALLOCATABLE :: juser_mult16(:)
INTEGER, ALLOCATABLE :: juser_mult17(:)
INTEGER, ALLOCATABLE :: juser_mult18(:)
INTEGER, ALLOCATABLE :: juser_mult19(:)
INTEGER, ALLOCATABLE :: juser_mult20(:)

! 1.9: Fields carried forward from previous version.
! Lateral Boundary Conditions
INTEGER :: jorog_lbc                       ! Orography LBC
INTEGER :: ju_lbc                          ! U LBC
INTEGER :: jv_lbc                          ! V LBC
INTEGER :: jw_lbc                          ! W LBC
INTEGER :: jrho_lbc                        ! RHO LBC
INTEGER :: jtheta_lbc                      ! Theta LBC
INTEGER :: jq_lbc                          ! Q LBC
INTEGER :: jqcl_lbc                        ! QCL LBC
INTEGER :: jqcf_lbc                        ! QCF LBC
INTEGER :: jqcf2_lbc                       ! 2nd Ice LBC
INTEGER :: jqrain_lbc                      ! Rain LBC
INTEGER :: jqgraup_lbc                     ! Graupel LBC
INTEGER :: jcf_bulk_lbc                    ! CF_BULK LBC
INTEGER :: jcf_liquid_lbc                  ! CF_LIQUID_LBC
INTEGER :: jcf_frozen_lbc                  ! CF_FROZEN_LBC
INTEGER :: jexner_lbc                      ! Exner LBC
INTEGER :: ju_adv_lbc                      ! U_ADV LBC
INTEGER :: jv_adv_lbc                      ! V_ADV LBC
INTEGER :: jw_adv_lbc                      ! W_ADV LBC
INTEGER :: jmurk_lbc                       ! Murk aerosol LBC
INTEGER, ALLOCATABLE :: jtracer_lbc (:)    ! Tracer LBCs
INTEGER, ALLOCATABLE :: jtr_ukca_lbc(:)    ! UKCA Tracer LBCs
INTEGER :: jdust_div1_lbc                  ! DUST_DIV1 LBC
INTEGER :: jdust_div2_lbc                  ! DUST_DIV2 LBC
INTEGER :: jdust_div3_lbc                  ! DUST_DIV3 LBC
INTEGER :: jdust_div4_lbc                  ! DUST_DIV4 LBC
INTEGER :: jdust_div5_lbc                  ! DUST_DIV5 LBC
INTEGER :: jdust_div6_lbc                  ! DUST_DIV6 LBC
INTEGER :: jso2_lbc                        ! SO2 LBC
INTEGER :: jdms_lbc                        ! DMS LBC
INTEGER :: jso4_aitken_lbc                 ! SO4_AITKEN LBC
INTEGER :: jso4_accu_lbc                   ! SO4_ACCU LBC
INTEGER :: jso4_diss_lbc                   ! SO4_DISS_LBC
INTEGER :: jnh3_lbc                        ! NH3 LBC
INTEGER :: jsoot_new_lbc                   ! SOOT_NEW LBC
INTEGER :: jsoot_agd_lbc                   ! SOOT_AGD LBC
INTEGER :: jsoot_cld_lbc                   ! SOOT_CLD LBC
INTEGER :: jbmass_new_lbc                  ! BMASS_NEW LBC
INTEGER :: jbmass_agd_lbc                  ! BMASS_AGD LBC
INTEGER :: jbmass_cld_lbc                  ! BMASS_CLD LBC
INTEGER :: jocff_new_lbc                   ! OCFF_NEW LBC
INTEGER :: jocff_agd_lbc                   ! OCFF_AGD LBC
INTEGER :: jocff_cld_lbc                   ! OCFF_CLD LBC
INTEGER :: jnitr_acc_lbc                   ! NITR_ACC_LBC
INTEGER :: jnitr_diss_lbc                  ! NITR_DISS_LBC
! Lateral Boundary Condition tendencies
INTEGER :: ju_lbc_tend                     ! U LBC  tendencies
INTEGER :: jv_lbc_tend                     ! V LBC tendencies
INTEGER :: jw_lbc_tend                     ! W LBC tendencies
INTEGER :: jrho_lbc_tend                   ! RHO LBC tendencies
INTEGER :: jtheta_lbc_tend                 ! Theta LBC tendencies
INTEGER :: jq_lbc_tend                     ! Q LBC tendencies
INTEGER :: jqcl_lbc_tend                   ! QCL LBC tendencies
INTEGER :: jqcf_lbc_tend                   ! QCF LBC tendencies
INTEGER :: jqcf2_lbc_tend                  ! 2nd Ice
INTEGER :: jqrain_lbc_tend                 ! Rain LBC tendencies
INTEGER :: jqgraup_lbc_tend                ! Graupel
INTEGER :: jcf_bulk_lbc_tend               ! CF_BULK LBC tend'cies
INTEGER :: jcf_liquid_lbc_tend             ! CF_LIQUID_LBC t'cies
INTEGER :: jcf_frozen_lbc_tend             ! CF_FROZEN_LBC t'cies
INTEGER :: jexner_lbc_tend                 ! Exner LBC tendencies
INTEGER :: ju_adv_lbc_tend                 ! U_ADV LBC tendencies
INTEGER :: jv_adv_lbc_tend                 ! V_ADV LBC tendencies
INTEGER :: jw_adv_lbc_tend                 ! W_ADV LBCtendencies
INTEGER :: jmurk_lbc_tend                  ! Murk aerosol LBC tend
INTEGER, ALLOCATABLE :: jtracer_lbc_tend(:) ! Tracer LBC tendencies
INTEGER, ALLOCATABLE :: jtr_ukca_lbc_tend(:)! UKCA Tracer LBC tend
INTEGER :: jdust_div1_lbc_tend             ! DUST_DIV1 LBC tend
INTEGER :: jdust_div2_lbc_tend             ! DUST_DIV2 LBC tend
INTEGER :: jdust_div3_lbc_tend             ! DUST_DIV3 LBC tend
INTEGER :: jdust_div4_lbc_tend             ! DUST_DIV4 LBC tend
INTEGER :: jdust_div5_lbc_tend             ! DUST_DIV5 LBC tend
INTEGER :: jdust_div6_lbc_tend             ! DUST_DIV6 LBC tend
INTEGER :: jso2_lbc_tend                   ! SO2 LBC tend
INTEGER :: jdms_lbc_tend                   ! DMS LBC tend
INTEGER :: jso4_aitken_lbc_tend            ! SO4_AITKEN LBC tend
INTEGER :: jso4_accu_lbc_tend              ! SO4_ACCU LBC tend
INTEGER :: jso4_diss_lbc_tend              ! SO4_DISS_LBC tend
INTEGER :: jnh3_lbc_tend                   ! NH3 LBC tend
INTEGER :: jsoot_new_lbc_tend              ! SOOT_NEW LBC tend
INTEGER :: jsoot_agd_lbc_tend              ! SOOT_AGD LBC tend
INTEGER :: jsoot_cld_lbc_tend              ! SOOT_CLD LBC tend
INTEGER :: jbmass_new_lbc_tend             ! BMASS_NEW LBC tend
INTEGER :: jbmass_agd_lbc_tend             ! BMASS_AGD LBC tend
INTEGER :: jbmass_cld_lbc_tend             ! BMASS_CLD LBC tend
INTEGER :: jocff_new_lbc_tend              ! OCFF_NEW LBC tend
INTEGER :: jocff_agd_lbc_tend              ! OCFF_AGD LBC tend
INTEGER :: jocff_cld_lbc_tend              ! OCFF_CLD LBC tend
INTEGER :: jnitr_acc_lbc_tend              ! NITR_ACC_LBC tend
INTEGER :: jnitr_diss_lbc_tend             ! NITR_DISS_LBC tend

! 2: Scalar Variables

! 2.1: Data variables stored in primary space.
INTEGER :: jtstar          ! Surface temperature
INTEGER :: jland           ! Land sea mask
INTEGER :: jtstar_anom     ! Surface temperature anomaly
!   2.15: Fields for coastal tiling
INTEGER :: jfrac_land      ! Land fraction in grid box
INTEGER :: jtstar_land     ! Land surface temperature
INTEGER :: jtstar_sea      ! Sea surface temperature
INTEGER :: jtstar_sice     ! Sea-ice surface temperature
INTEGER :: jtstar_sice_cat ! Sea-ice surface temperature on categories
INTEGER :: jsice_alb       ! Sea-ice albedo
INTEGER :: jland_alb       ! Mean land albedo

! 2.2: Data variables stored in secondary space.

INTEGER :: jpstar          ! Surface pressure

! 2.3: Cloud fields
INTEGER :: jlcbase         ! Lowest Convective cloud base
INTEGER :: jccb            ! Convective cloud base
INTEGER :: jcct            ! Convective cloud top

INTEGER :: jcclwp          ! Convective cloud liquid water path

INTEGER :: jdeepflag       ! Flag for history of deep convection
INTEGER :: jpastprecip     ! Past convective precipitation
INTEGER :: jpastconvht     ! Past convective height

! 2.4: Boundary layer fields

INTEGER :: jzh             ! Boundary layer depth
INTEGER :: jddmfx          ! Convective downdraught mass-flux at cloud-base

INTEGER :: jt1_sd          ! Standard deviation of turbulent fluctuations 
                           ! of layer 1 temperature
INTEGER :: jq1_sd          ! Standard deviation of turbulent fluctuations 
                           ! of layer 1 humidity
INTEGER :: jtscrndcl_tile  ! Decoupled screen-level temperature
INTEGER :: jtscrndcl_ssi
INTEGER :: jtstbtrans
INTEGER :: jntml           ! Number of model levels in the 
                           ! turbulently mixed layer
INTEGER :: jntdsc          ! Top level for turb mixing in any decoupled 
                           ! Sc layer
INTEGER :: jnbdsc          ! Bottom level for turb mixing in any decoupled 
                           ! Sc layer
INTEGER :: jcumulus        ! Boundary layer convection flag

! Convective cold-pool prognostics
INTEGER :: jux_ccp   ! x cmpt front speed vector sum
INTEGER :: juy_ccp   ! y cmpt front speed vector sum
INTEGER :: jum_ccp   ! front speed scalar sum
INTEGER :: jg_ccp    ! gridbox c.c.p. reduced gravity
INTEGER :: jh_ccp    ! gridbox c.c.p. depth
INTEGER :: jriso_ccp ! remain counter (isotropic)
INTEGER :: jrdir_ccp ! remain counter (directed)

! 2.4: Soil Ancillary fields

INTEGER :: jsat_soilw_suction ! Saturated soil water sucti
INTEGER :: jtherm_cap         ! Thermal capacity
INTEGER :: jtherm_cond        ! Thermal conductivity
INTEGER :: jvol_smc_crit      ! Vol smc at critical point
INTEGER :: jvol_smc_wilt      ! Vol smc at wilting point
INTEGER :: jvol_smc_sat       ! Vol smc at saturation
INTEGER :: jsat_soil_cond     ! Saturated soil conductivity
INTEGER :: jclapp_horn        ! Clapp-Hornberger B coefficient
INTEGER :: jz0m_soil          ! Bare soil roughness length (z0) for momentum

! 2.5: Other surface fields
INTEGER :: jcanopy_water      ! Canopy Water
INTEGER :: jz0                ! Roughness length; sea points on first timestep
INTEGER :: jgs                ! Gridbox mean canopy conductance

! 2.6: Orographic Ancillary fields

INTEGER :: jorog              ! Orographic height
INTEGER :: jorog_sd           ! Standard Deviation of orography
INTEGER :: jorog_sil          ! Silhouette area of orography
INTEGER :: jorog_ho2          ! Peak to trough height/(2*sqrt2)
INTEGER :: jorog_grad_x       ! Orographic gradient x
INTEGER :: jorog_grad_y       ! Orographic gradient y
INTEGER :: jorog_unfilt       ! Unfiltered orographic height
INTEGER :: jorog_grad_xx      ! Orographic gradient xx
INTEGER :: jorog_grad_xy      ! Orographic gradient xy
INTEGER :: jorog_grad_yy      ! Orographic gradient yy

! 2.7: Sea/Sea Ice

INTEGER :: ju_sea             ! Surface current (u component)
INTEGER :: jv_sea             ! Surface current (v component)
INTEGER :: ju_0_p             ! Surace  current (u) on p-grid
INTEGER :: jv_0_p             ! Surface current (v) on p-grid
INTEGER :: jice_fraction      ! Sea ice fraction
INTEGER :: jice_thickness     ! Sea ice depth
INTEGER :: jti                ! Sea ice temperature
INTEGER :: jice_fract_cat     ! Sea ice fraction on categories
INTEGER :: jice_thick_cat     ! Sea ice thickness on categories
INTEGER :: jpond_frac_cat     ! Sea ice meltpond fraction on categories 
INTEGER :: jpond_depth_cat    ! Sea ice meltpond depth on categories 
INTEGER :: jtfrz              ! Sea surface freezing temperature 
INTEGER :: jti_cat            ! Sea ice temperature on categories
INTEGER :: jice_k_cat         ! Sea ice effect conductivity on categories
INTEGER :: jchloro_sea        ! Sea near surface chlorophyll

! 2.8: Snow

INTEGER :: jsnodep            ! Snow depth on land
INTEGER :: jsnodep_sea        ! Snow depth on sea ice
INTEGER :: jsnodep_sea_cat    ! Snow depth on sea ice catagories
INTEGER :: jcatch_snow        ! Coniferous canopy snow capacity
INTEGER :: jsnow_grnd         ! Snow below canopy
INTEGER :: jsnsoot            ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

INTEGER :: jsoil_clay         ! soil clay fraction
INTEGER :: jsoil_silt         ! soil silt fraction
INTEGER :: jsoil_sand         ! soil sand fraction
INTEGER :: jdust_mrel1        ! soil rel mass, div 1
INTEGER :: jdust_mrel2        ! soil rel mass, div 2
INTEGER :: jdust_mrel3        ! soil rel mass, div 3
INTEGER :: jdust_mrel4        ! soil rel mass, div 4
INTEGER :: jdust_mrel5        ! soil rel mass, div 5
INTEGER :: jdust_mrel6        ! soil rel mass, div 6


INTEGER :: jso2_em            ! sulphur dioxide emission
INTEGER :: jdms_em            ! dimethyl sulphide emission
INTEGER :: jso2_hilem         ! high level SO2 emissions
INTEGER :: jnh3_em            ! ammonia gas surface emiss
INTEGER :: jsoot_em           ! fresh soot surface emissions
INTEGER :: jsoot_hilem        ! fresh soot high lev emissions
INTEGER :: jbmass_em          ! fresh bmass surface emissions
INTEGER :: jbmass_hilem       ! fresh bmass high lev emissions
INTEGER :: jocff_em           ! fresh OCFF surface emissions
INTEGER :: jocff_hilem        ! fresh OCFF high-level emissions
INTEGER :: jdms_conc          ! seawater dimethyl sulphide conc.
INTEGER :: jbmass_hilem_h1    ! min height for fresh bmass high lev emiss
INTEGER :: jbmass_hilem_h2    ! max height for fresh bmass high lev emiss

! 2.10: User ancillary fields
INTEGER :: juser_anc1
INTEGER :: juser_anc2
INTEGER :: juser_anc3
INTEGER :: juser_anc4
INTEGER :: juser_anc5
INTEGER :: juser_anc6
INTEGER :: juser_anc7
INTEGER :: juser_anc8
INTEGER :: juser_anc9
INTEGER :: juser_anc10
INTEGER :: juser_anc11
INTEGER :: juser_anc12
INTEGER :: juser_anc13
INTEGER :: juser_anc14
INTEGER :: juser_anc15
INTEGER :: juser_anc16
INTEGER :: juser_anc17
INTEGER :: juser_anc18
INTEGER :: juser_anc19
INTEGER :: juser_anc20

!   2.11: Store arrays for energy correction calculation
INTEGER :: jnet_flux          ! Net energy flux
INTEGER :: jnet_mflux         ! Net moisture flux

!   2.12: Tiled Vegetation and Triffid fields
INTEGER :: jfrac_typ          ! Fractions of surface type
INTEGER :: jfrac_con1         ! Fractions of surface type
INTEGER :: jfrac_con2         ! Fractions of surface type
INTEGER :: jfrac_con3         ! Fractions of surface type
INTEGER :: jfrac_con4         ! Fractions of surface type
INTEGER :: jfrac_con5         ! Fractions of surface type
INTEGER :: jfrac_con6         ! Fractions of surface type
INTEGER :: jfrac_con7         ! Fractions of surface type
INTEGER :: jfrac_con8         ! Fractions of surface type
INTEGER :: jfrac_con9         ! Fractions of surface type
INTEGER :: jlai_pft           ! LAI of plant functional types
INTEGER :: jcanht_pft         ! Canopy hght of plant func types
INTEGER :: jdisturb           ! Disturbed fraction of vegetation
INTEGER :: jdisturb_prev      ! Previous disturbed fraction of vegetation
INTEGER :: jpasture           ! Pasture fraction
INTEGER :: jpasture_prev      ! Previous pasture fraction
INTEGER :: jagr_crop          ! Crop fraction
INTEGER :: jagr_crop_prev     ! Previous Crop fraction
INTEGER :: jwoodprod_fast     ! Wood product pool (fast)                 
INTEGER :: jwoodprod_med      ! Wood product pool (med)                  
INTEGER :: jwoodprod_slow     ! Wood product pool (slow)                 
INTEGER :: jsoil_alb          ! Snow-free albedo of bare soil
INTEGER :: jobs_alb_sw        ! Observed snow-free SW albedo
INTEGER :: jobs_alb_vis       ! Observed snow-free VIS albedo
INTEGER :: jobs_alb_nir       ! Observed snow-free NIR albedo
INTEGER :: jsoil_carb         ! Soil carbon content
INTEGER :: jsoil_carb1        ! Soil carbon content DPM
INTEGER :: jsoil_carb2        ! Soil carbon content RPM
INTEGER :: jsoil_carb3        ! Soil carbon content BIO
INTEGER :: jsoil_carb4        ! Soil carbon content HUM
INTEGER :: jsoil_nitro1       ! Soil nitrogen content DPM
INTEGER :: jsoil_nitro2       ! Soil nitrogen content RPM
INTEGER :: jsoil_nitro3       ! Soil nitrogen content BIO
INTEGER :: jsoil_nitro4       ! Soil nitrogen content HUM
INTEGER :: jsoil_inorgnit     ! Soil inorganic nitrogen content
INTEGER :: j_n_deposition     ! Nitrogen deposition
INTEGER :: jnpp_pft_acc       ! Accumulated NPP on PFTs
INTEGER :: jg_lf_pft_acc      ! Accum. leaf turnover rate PFTs
INTEGER :: jg_phlf_pft_acc    ! Accumulated phenological leaf 
                              ! turnover rate on PFTs
INTEGER :: jrsp_w_pft_acc     ! Accum. wood respiration on PFTs
INTEGER :: jrsp_s_acc         ! Accumulated soil respiration
INTEGER :: jrsp_s_acc1        ! Accumulated soil respiration DPM
INTEGER :: jrsp_s_acc2        ! Accumulated soil respiration RPM
INTEGER :: jrsp_s_acc3        ! Accumulated soil respiration BIO
INTEGER :: jrsp_s_acc4        ! Accumulated soil respiration HUM
INTEGER :: jcan_water_tile    ! Canopy water content on tiles
INTEGER :: jcatch_tile        ! Canopy capacity on tiles
INTEGER :: jinfil_tile        ! Max infiltration rate on tiles
INTEGER :: jrgrain_tile       ! Snow grain size on tiles
INTEGER :: jsnodep_tile       ! Snow depth on tiles
INTEGER :: jtstar_tile        ! Surface temperature on tiles
INTEGER :: jtsurf_elev_surft  ! Temperature of elevated subsurface tiles (K)
INTEGER :: jz0_tile           ! Surface roughness on tiles
INTEGER :: jz0h_tile          ! Surface thermal roughness on tiles
INTEGER :: jdolr              ! TOA - surface upward LW at radiation timestep
INTEGER :: jlw_down           ! Surface downward LW at radiation timestep
INTEGER :: jsw_tile           ! Surface net SW on land tiles at 
                              ! radiation timestep
INTEGER :: jurbhgt            ! Building height
INTEGER :: jurbhwr            ! Urban H/W ratio
INTEGER :: jurbwrr            ! Width ratio
INTEGER :: jurbdisp           ! Displacement height
INTEGER :: jurbztm
INTEGER :: jurbalbwl          ! Wall albedo
INTEGER :: jurbalbrd          ! Road albedo
INTEGER :: jurbemisw          ! Wall emmissivity
INTEGER :: jurbemisr          ! Road emmissivity

!   2.14: Carbon cycle fields
INTEGER :: j_co2flux          ! Ocean CO2 flux (Kg CO2/m2/s1)
INTEGER :: j_co2_emits        ! Surface CO2 emissions (Kg CO2/m2/s1)
INTEGER :: j_triffid_co2      ! TRIFFID-derived CO2 fluxes for passing to the 
                              ! atmosphere in emissions-driven runs (kgC/m2/yr):
                              ! exudates + wood product pool flux + harvest flux

! Indices for ATMOSPHERE model constants. Scalars only.
INTEGER :: jetatheta
INTEGER :: jetarho
INTEGER :: jrhcrit
INTEGER :: jsoil_thickness
! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
INTEGER :: jzseak_theta       ! zsea(k) on theta levels
INTEGER :: jck_theta          ! C(k)    on theta levels
INTEGER :: jzseak_rho         ! zsea(k) on rho levels
INTEGER :: jck_rho            ! C(k)    on rho levels
! Addresses in Row and Col  dependent constants array.
INTEGER :: jlambda_input_p
INTEGER :: jlambda_input_u
INTEGER :: jphi_input_p
INTEGER :: jphi_input_v

!   2.16: Fields for large-scale hydrology scheme.
INTEGER :: jti_mean           ! Mean topographic index
INTEGER :: jti_sig            ! Standard dev. in topographic index
INTEGER :: jfexp              ! Exponential decay in soil saturated conductivity
INTEGER :: jgamma_int         ! Integrated gamma function
INTEGER :: jwater_table       ! Water table depth
INTEGER :: jfsfc_sat          ! Surface saturation fraction
INTEGER :: jf_wetland         ! Wetland fraction
INTEGER :: jsthzw             ! Soil moist fract. in deep-zw layer.
INTEGER :: ja_fsat            ! Fitting parameter for Fsat in LSH.
INTEGER :: jc_fsat            ! Fitting parameter for Fsat in LSH.
INTEGER :: ja_fwet            ! Fitting parameter for Fwet in LSH.
INTEGER :: jc_fwet            ! Fitting parameter for Fwet in LSH.

!   2.17: Fields for River routing.
INTEGER :: jriv_sequence      ! River sequence
INTEGER :: jriv_direction     ! River direction
INTEGER :: jriv_storage       ! River water storage
INTEGER :: jtot_surfroff      ! Accumulated surface runoff
INTEGER :: jtot_subroff       !     "       sub-surface runoff
INTEGER :: jriv_inlandatm     ! inland basin outflow
! Field for water conservation due to lake evaporation
INTEGER :: jacc_lake_evap     ! Accumulated lake evaporation

INTEGER :: jc_solar
INTEGER :: jc_blue
INTEGER :: jc_longwave
INTEGER :: jc_taux
INTEGER :: jc_tauy
INTEGER :: jc_w10
INTEGER :: jc_sensible
INTEGER :: jc_sublim
INTEGER :: jc_evap
INTEGER :: jc_fcondtopn
INTEGER :: jc_topmeltn
INTEGER :: jc_tstar_sicen
INTEGER :: jc_lsrain
INTEGER :: jc_lssnow
INTEGER :: jc_cvrain
INTEGER :: jc_cvsnow
INTEGER :: jc_riverout
INTEGER :: jc_mslp
INTEGER :: jc_calving
INTEGER :: jc_surf_co2
INTEGER :: jc_dust_dep

!   2.18: JULES variables
INTEGER :: jsnowdepth         ! Snow depth on ground on tiles (m)
INTEGER :: jrho_snow_grnd     ! Snowpack bulk density (kg/m3)
INTEGER :: jnsnow             ! Number of snow layers on ground on tiles
INTEGER :: jds                ! Snow layer thickness (m)
INTEGER :: jsice              ! Snow layer ice mass on tiles (Kg/m2)
INTEGER :: jsliq              ! Snow layer liquid mass on tiles (Kg/m2)
INTEGER :: jtsnowlayer        ! Snow layer temperature (K)
INTEGER :: jrho_snow          ! Snow layer densities (kg/m3)
INTEGER :: jrgrainl           ! Snow layer grain size on tiles (microns)
!  FLake lake scheme
INTEGER :: jlake_depth
INTEGER :: jlake_fetch
INTEGER :: jlake_t_mean
INTEGER :: jlake_t_mxl
INTEGER :: jlake_t_ice
INTEGER :: jlake_h_mxl
INTEGER :: jlake_h_ice
INTEGER :: jlake_shape
INTEGER :: jlake_g_dt

CHARACTER(LEN=*), PARAMETER, PRIVATE :: modulename = 'ATM_D1_INDICES_MOD'

CONTAINS

SUBROUTINE allocate_d1_indices()

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),   PARAMETER :: routinename = 'ALLOCATE_D1_INDICES'


IF (ALLOCATED(jdryrho)) DEALLOCATE(jdryrho)
IF (ALLOCATED(jetadot)) DEALLOCATE(jetadot)
IF (ALLOCATED(jthetav)) DEALLOCATE(jthetav)
IF (ALLOCATED(jmv))     DEALLOCATE(jmv)
IF (ALLOCATED(jmcl))    DEALLOCATE(jmcl)
IF (ALLOCATED(jmcf))    DEALLOCATE(jmcf)
IF (ALLOCATED(jmcf2))   DEALLOCATE(jmcf2)
IF (ALLOCATED(jmrain))  DEALLOCATE(jmrain)
IF (ALLOCATED(jmgraup)) DEALLOCATE(jmgraup)
IF (ALLOCATED(ju))      DEALLOCATE(ju)
IF (ALLOCATED(jv))      DEALLOCATE(jv)
IF (ALLOCATED(jw))      DEALLOCATE(jw)
IF (ALLOCATED(jrho))    DEALLOCATE(jrho)
IF (ALLOCATED(jtheta))  DEALLOCATE(jtheta)
IF (ALLOCATED(jq))      DEALLOCATE(jq)
IF (ALLOCATED(jqcl))    DEALLOCATE(jqcl)
IF (ALLOCATED(jqcf))    DEALLOCATE(jqcf)
IF (ALLOCATED(jqcf2))   DEALLOCATE(jqcf2)
IF (ALLOCATED(jqrain))  DEALLOCATE(jqrain)
IF (ALLOCATED(jqgraup)) DEALLOCATE(jqgraup)

IF (ALLOCATED(jcloudnumber)) DEALLOCATE(jcloudnumber)
IF (ALLOCATED(jrainnumber))  DEALLOCATE(jrainnumber)
IF (ALLOCATED(jrain3mom))    DEALLOCATE(jrain3mom)
IF (ALLOCATED(jicenumber))   DEALLOCATE(jicenumber)
IF (ALLOCATED(jsnownumber))  DEALLOCATE(jsnownumber)
IF (ALLOCATED(jsnow3mom))    DEALLOCATE(jsnow3mom)
IF (ALLOCATED(jgraupnumber)) DEALLOCATE(jgraupnumber)
IF (ALLOCATED(jgraup3mom))   DEALLOCATE(jgraup3mom)


IF (ALLOCATED(jactivesolliquid))   DEALLOCATE(jactivesolliquid)
IF (ALLOCATED(jactivesolrain))     DEALLOCATE(jactivesolrain)
IF (ALLOCATED(jactiveinsolice))    DEALLOCATE(jactiveinsolice)
IF (ALLOCATED(jactivesolice))      DEALLOCATE(jactivesolice)
IF (ALLOCATED(jactiveinsolliquid)) DEALLOCATE(jactiveinsolliquid)
IF (ALLOCATED(jactivesolnumber))   DEALLOCATE(jactivesolnumber)
IF (ALLOCATED(jactiveinsolnumber)) DEALLOCATE(jactiveinsolnumber)


IF (ALLOCATED(je_trb))            DEALLOCATE(je_trb)
IF (ALLOCATED(jtsq_trb))          DEALLOCATE(jtsq_trb)
IF (ALLOCATED(jqsq_trb))          DEALLOCATE(jqsq_trb)
IF (ALLOCATED(jcov_trb))          DEALLOCATE(jcov_trb)
IF (ALLOCATED(jexner_rho_levels)) DEALLOCATE(jexner_rho_levels)
IF (ALLOCATED(jconv_prog_1))      DEALLOCATE(jconv_prog_1)
IF (ALLOCATED(jconv_prog_2))      DEALLOCATE(jconv_prog_2)
IF (ALLOCATED(jconv_prog_3))      DEALLOCATE(jconv_prog_3)
IF (ALLOCATED(jconv_prog_precip)) DEALLOCATE(jconv_prog_precip)

IF (ALLOCATED(jdPV_rad))          DEALLOCATE(jdPV_rad)
IF (ALLOCATED(jdPV_sw))           DEALLOCATE(jdPV_sw)
IF (ALLOCATED(jdPV_lw))           DEALLOCATE(jdPV_lw)
IF (ALLOCATED(jdPV_mic))          DEALLOCATE(jdPV_mic)
IF (ALLOCATED(jdPV_gwd))          DEALLOCATE(jdPV_gwd)
IF (ALLOCATED(jdPV_ph1))          DEALLOCATE(jdPV_ph1)
IF (ALLOCATED(jdPV_conv))         DEALLOCATE(jdPV_conv)
IF (ALLOCATED(jdPV_bl))           DEALLOCATE(jdPV_bl)
IF (ALLOCATED(jdPV_stph))         DEALLOCATE(jdPV_stph)
IF (ALLOCATED(jdPV_cld))          DEALLOCATE(jdPV_cld)
IF (ALLOCATED(jdPV_iau))          DEALLOCATE(jdPV_iau)
IF (ALLOCATED(jdPV_nud))          DEALLOCATE(jdPV_nud)
IF (ALLOCATED(jdPV_tot))          DEALLOCATE(jdPV_tot)
IF (ALLOCATED(jdPV_adv))          DEALLOCATE(jdPV_adv)
IF (ALLOCATED(jdPV_sol))          DEALLOCATE(jdPV_sol)
IF (ALLOCATED(jdPV_mass))         DEALLOCATE(jdPV_mass)
IF (ALLOCATED(jadv_only_PV))      DEALLOCATE(jadv_only_PV)

IF (ALLOCATED(jdtheta_0))         DEALLOCATE(jdtheta_0)
IF (ALLOCATED(jdtheta_bl))        DEALLOCATE(jdtheta_bl)
IF (ALLOCATED(jdtheta_bl_mix))    DEALLOCATE(jdtheta_bl_mix)
IF (ALLOCATED(jdtheta_bl_LH))     DEALLOCATE(jdtheta_bl_LH)
IF (ALLOCATED(jdtheta_conv))      DEALLOCATE(jdtheta_conv)
IF (ALLOCATED(jdtheta_mic))       DEALLOCATE(jdtheta_mic)
IF (ALLOCATED(jdtheta_rad))       DEALLOCATE(jdtheta_rad)
IF (ALLOCATED(jdtheta_SW))        DEALLOCATE(jdtheta_SW)
IF (ALLOCATED(jdtheta_LW))        DEALLOCATE(jdtheta_LW)
IF (ALLOCATED(jdtheta_slow))      DEALLOCATE(jdtheta_slow)
IF (ALLOCATED(jdtheta_cld))       DEALLOCATE(jdtheta_cld)

IF (ALLOCATED(jp))                  DEALLOCATE(jp)
IF (ALLOCATED(jp_theta_levels))     DEALLOCATE(jp_theta_levels)
IF (ALLOCATED(jexner_theta_levels)) DEALLOCATE(jexner_theta_levels)
IF (ALLOCATED(jccw_rad))            DEALLOCATE(jccw_rad)
IF (ALLOCATED(jcca))                DEALLOCATE(jcca)
IF (ALLOCATED(jcca_dp))             DEALLOCATE(jcca_dp)
IF (ALLOCATED(jcca_md))             DEALLOCATE(jcca_md)
IF (ALLOCATED(jcca_sh))             DEALLOCATE(jcca_sh)

IF (ALLOCATED(jcf_area))         DEALLOCATE(jcf_area)
IF (ALLOCATED(jcf_bulk))         DEALLOCATE(jcf_bulk)
IF (ALLOCATED(jcf_liquid))       DEALLOCATE(jcf_liquid)
IF (ALLOCATED(jcf_frozen))       DEALLOCATE(jcf_frozen)
IF (ALLOCATED(j_deep_soil_temp)) DEALLOCATE(j_deep_soil_temp)
IF (ALLOCATED(jsmcl))            DEALLOCATE(jsmcl)
IF (ALLOCATED(jsthu))            DEALLOCATE(jsthu)
IF (ALLOCATED(jsthf))            DEALLOCATE(jsthf)

IF (ALLOCATED(jsw_incs))   DEALLOCATE(jsw_incs)
IF (ALLOCATED(jlw_incs))   DEALLOCATE(jlw_incs)
IF (ALLOCATED(jozone))     DEALLOCATE(jozone)
IF (ALLOCATED(jtppsozone)) DEALLOCATE(jtppsozone)

IF (ALLOCATED(jtracer))      DEALLOCATE(jtracer)
IF (ALLOCATED(jtr_ukca))     DEALLOCATE(jtr_ukca)
IF (ALLOCATED(jmurk_source)) DEALLOCATE(jmurk_source)
IF (ALLOCATED(jmurk))        DEALLOCATE(jmurk)
IF (ALLOCATED(jflash_pot))   DEALLOCATE(jflash_pot)
IF (ALLOCATED(jbl_w_var))    DEALLOCATE(jbl_w_var)

IF (ALLOCATED(jdust_div1)) DEALLOCATE(jdust_div1)
IF (ALLOCATED(jdust_div2)) DEALLOCATE(jdust_div2)
IF (ALLOCATED(jdust_div3)) DEALLOCATE(jdust_div3)
IF (ALLOCATED(jdust_div4)) DEALLOCATE(jdust_div4)
IF (ALLOCATED(jdust_div5)) DEALLOCATE(jdust_div5)
IF (ALLOCATED(jdust_div6)) DEALLOCATE(jdust_div6)


IF (ALLOCATED(jso2))        DEALLOCATE(jso2)
IF (ALLOCATED(jdms))        DEALLOCATE(jdms)
IF (ALLOCATED(jso4_aitken)) DEALLOCATE(jso4_aitken)
IF (ALLOCATED(jso4_accu))   DEALLOCATE(jso4_accu)
IF (ALLOCATED(jso4_diss))   DEALLOCATE(jso4_diss)
IF (ALLOCATED(jh2o2))       DEALLOCATE(jh2o2)
IF (ALLOCATED(jnh3))        DEALLOCATE(jnh3)

IF (ALLOCATED(jsoot_new))  DEALLOCATE(jsoot_new)
IF (ALLOCATED(jsoot_agd))  DEALLOCATE(jsoot_agd)
IF (ALLOCATED(jsoot_cld))  DEALLOCATE(jsoot_cld)
IF (ALLOCATED(jbmass_new)) DEALLOCATE(jbmass_new)
IF (ALLOCATED(jbmass_agd)) DEALLOCATE(jbmass_agd)
IF (ALLOCATED(jbmass_cld)) DEALLOCATE(jbmass_cld)
IF (ALLOCATED(jocff_new))  DEALLOCATE(jocff_new)
IF (ALLOCATED(jocff_agd))  DEALLOCATE(jocff_agd)
IF (ALLOCATED(jocff_cld))  DEALLOCATE(jocff_cld)



IF (ALLOCATED(jso2_natem))  DEALLOCATE(jso2_natem)
IF (ALLOCATED(joh))         DEALLOCATE(joh)
IF (ALLOCATED(jho2))        DEALLOCATE(jho2)
IF (ALLOCATED(jh2o2_limit)) DEALLOCATE(jh2o2_limit)
IF (ALLOCATED(jo3_chem))    DEALLOCATE(jo3_chem)
IF (ALLOCATED(jco2))        DEALLOCATE(jco2)


IF (ALLOCATED(joh_ukca))      DEALLOCATE(joh_ukca)
IF (ALLOCATED(jho2_ukca))     DEALLOCATE(jho2_ukca)
IF (ALLOCATED(jh2o2_ukca))    DEALLOCATE(jh2o2_ukca)
IF (ALLOCATED(jo3_ukca))      DEALLOCATE(jo3_ukca)
IF (ALLOCATED(jhno3_ukca))    DEALLOCATE(jhno3_ukca)
IF (ALLOCATED(jozone_tracer)) DEALLOCATE(jozone_tracer)
IF (ALLOCATED(jo3_prod_loss)) DEALLOCATE(jo3_prod_loss)
IF (ALLOCATED(jo3_p_l_vmr))   DEALLOCATE(jo3_p_l_vmr)
IF (ALLOCATED(jo3_vmr))       DEALLOCATE(jo3_vmr)
IF (ALLOCATED(jo3_p_l_temp))  DEALLOCATE(jo3_p_l_temp)
IF (ALLOCATED(jo3_temp))      DEALLOCATE(jo3_temp)
IF (ALLOCATED(jo3_p_l_colo3)) DEALLOCATE(jo3_p_l_colo3)
IF (ALLOCATED(jo3_colo3))     DEALLOCATE(jo3_colo3)


IF (ALLOCATED(jarclbiog_bg)) DEALLOCATE(jarclbiog_bg)
IF (ALLOCATED(jarclbiom_fr)) DEALLOCATE(jarclbiom_fr)
IF (ALLOCATED(jarclbiom_ag)) DEALLOCATE(jarclbiom_ag)
IF (ALLOCATED(jarclbiom_ic)) DEALLOCATE(jarclbiom_ic)
IF (ALLOCATED(jarclblck_fr)) DEALLOCATE(jarclblck_fr)
IF (ALLOCATED(jarclblck_ag)) DEALLOCATE(jarclblck_ag)

IF (ALLOCATED(jarclsslt_fi)) DEALLOCATE(jarclsslt_fi)
IF (ALLOCATED(jarclsslt_jt)) DEALLOCATE(jarclsslt_jt)

IF (ALLOCATED(jarclsulp_ac)) DEALLOCATE(jarclsulp_ac)
IF (ALLOCATED(jarclsulp_ak)) DEALLOCATE(jarclsulp_ak)
IF (ALLOCATED(jarclsulp_di)) DEALLOCATE(jarclsulp_di)

IF (ALLOCATED(jarcldust_b1)) DEALLOCATE(jarcldust_b1)
IF (ALLOCATED(jarcldust_b2)) DEALLOCATE(jarcldust_b2)
IF (ALLOCATED(jarcldust_b3)) DEALLOCATE(jarcldust_b3)
IF (ALLOCATED(jarcldust_b4)) DEALLOCATE(jarcldust_b4)
IF (ALLOCATED(jarcldust_b5)) DEALLOCATE(jarcldust_b5)
IF (ALLOCATED(jarcldust_b6)) DEALLOCATE(jarcldust_b6)

IF (ALLOCATED(jarclocff_fr)) DEALLOCATE(jarclocff_fr)
IF (ALLOCATED(jarclocff_ag)) DEALLOCATE(jarclocff_ag)
IF (ALLOCATED(jarclocff_ic)) DEALLOCATE(jarclocff_ic)

IF (ALLOCATED(jarcldlta_dl)) DEALLOCATE(jarcldlta_dl)
IF (ALLOCATED(jnitr_acc))    DEALLOCATE(jnitr_acc)
IF (ALLOCATED(jnitr_diss))   DEALLOCATE(jnitr_diss)


IF (ALLOCATED(jgc_nd_nuc_sol)) DEALLOCATE(jgc_nd_nuc_sol)
IF (ALLOCATED(jgc_nuc_sol_su)) DEALLOCATE(jgc_nuc_sol_su)
IF (ALLOCATED(jgc_nuc_sol_oc)) DEALLOCATE(jgc_nuc_sol_oc)
IF (ALLOCATED(jgc_nd_ait_sol)) DEALLOCATE(jgc_nd_ait_sol)
IF (ALLOCATED(jgc_ait_sol_su)) DEALLOCATE(jgc_ait_sol_su)
IF (ALLOCATED(jgc_ait_sol_bc)) DEALLOCATE(jgc_ait_sol_bc)
IF (ALLOCATED(jgc_ait_sol_oc)) DEALLOCATE(jgc_ait_sol_oc)

IF (ALLOCATED(jgc_nd_acc_sol)) DEALLOCATE(jgc_nd_acc_sol)
IF (ALLOCATED(jgc_acc_sol_su)) DEALLOCATE(jgc_acc_sol_su)
IF (ALLOCATED(jgc_acc_sol_bc)) DEALLOCATE(jgc_acc_sol_bc)
IF (ALLOCATED(jgc_acc_sol_oc)) DEALLOCATE(jgc_acc_sol_oc)
IF (ALLOCATED(jgc_acc_sol_ss)) DEALLOCATE(jgc_acc_sol_ss)

IF (ALLOCATED(jgc_nd_cor_sol)) DEALLOCATE(jgc_nd_cor_sol)
IF (ALLOCATED(jgc_cor_sol_su)) DEALLOCATE(jgc_cor_sol_su)
IF (ALLOCATED(jgc_cor_sol_bc)) DEALLOCATE(jgc_cor_sol_bc)
IF (ALLOCATED(jgc_cor_sol_oc)) DEALLOCATE(jgc_cor_sol_oc)
IF (ALLOCATED(jgc_cor_sol_ss)) DEALLOCATE(jgc_cor_sol_ss)

IF (ALLOCATED(jgc_nd_ait_ins)) DEALLOCATE(jgc_nd_ait_ins)
IF (ALLOCATED(jgc_ait_ins_bc)) DEALLOCATE(jgc_ait_ins_bc)
IF (ALLOCATED(jgc_ait_ins_oc)) DEALLOCATE(jgc_ait_ins_oc)

IF (ALLOCATED(juser_mult1))  DEALLOCATE(juser_mult1)
IF (ALLOCATED(juser_mult2))  DEALLOCATE(juser_mult2)
IF (ALLOCATED(juser_mult3))  DEALLOCATE(juser_mult3)
IF (ALLOCATED(juser_mult4))  DEALLOCATE(juser_mult4)
IF (ALLOCATED(juser_mult5))  DEALLOCATE(juser_mult5)
IF (ALLOCATED(juser_mult6))  DEALLOCATE(juser_mult6)
IF (ALLOCATED(juser_mult7))  DEALLOCATE(juser_mult7)
IF (ALLOCATED(juser_mult8))  DEALLOCATE(juser_mult8)
IF (ALLOCATED(juser_mult9))  DEALLOCATE(juser_mult9)
IF (ALLOCATED(juser_mult10)) DEALLOCATE(juser_mult10)
IF (ALLOCATED(juser_mult11)) DEALLOCATE(juser_mult11)
IF (ALLOCATED(juser_mult12)) DEALLOCATE(juser_mult12)
IF (ALLOCATED(juser_mult13)) DEALLOCATE(juser_mult13)
IF (ALLOCATED(juser_mult14)) DEALLOCATE(juser_mult14)
IF (ALLOCATED(juser_mult15)) DEALLOCATE(juser_mult15)
IF (ALLOCATED(juser_mult16)) DEALLOCATE(juser_mult16)
IF (ALLOCATED(juser_mult17)) DEALLOCATE(juser_mult17)
IF (ALLOCATED(juser_mult18)) DEALLOCATE(juser_mult18)
IF (ALLOCATED(juser_mult19)) DEALLOCATE(juser_mult19)
IF (ALLOCATED(juser_mult20)) DEALLOCATE(juser_mult20)

IF (ALLOCATED(jtracer_lbc))       DEALLOCATE(jtracer_lbc)
IF (ALLOCATED(jtr_ukca_lbc))      DEALLOCATE(jtr_ukca_lbc)
IF (ALLOCATED(jtracer_lbc_tend))  DEALLOCATE(jtracer_lbc_tend)
IF (ALLOCATED(jtr_ukca_lbc_tend)) DEALLOCATE(jtr_ukca_lbc_tend)

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

ALLOCATE(jdryrho (pdims%k_start:pdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdryrho]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jetadot (wdims%k_start:wdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jetadot]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jthetav (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jthetav]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmv (tdims%k_start:tdims%k_end)                                   &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmv]'//                       newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmcl (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmcl]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmcf (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmcf]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmcf2 (tdims%k_start:tdims%k_end)                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmcf2]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmrain (tdims%k_start:tdims%k_end)                                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmrain]'//                    newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmgraup (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmgraup]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(ju (udims%k_start:udims%k_end)                                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [ju]'//                        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jv (vdims%k_start:vdims%k_end)                                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jv]'//                        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jw (wdims_s%k_start:wdims_s%k_end)                                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jw]'//                        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jrho (pdims%k_start:pdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jrho]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtheta (tdims%k_start:tdims%k_end)                                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtheta]'//                    newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(jq (tdims%k_start:tdims%k_end)                                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jq]'//                        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqcl (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqcl]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqcf (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqcf]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqcf2 (tdims%k_start:tdims%k_end)                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqcf2]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqrain (tdims%k_start:tdims%k_end)                                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqrain]'//                    newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqgraup (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqgraup]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcloudnumber (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcloudnumber]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jrainnumber (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jrainnumber]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jrain3mom (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jrain3mom]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jicenumber (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jicenumber]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsnownumber (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsnownumber]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsnow3mom (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsnow3mom]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgraupnumber (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgraupnumber]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgraup3mom (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgraup3mom]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactivesolliquid (tdims%k_start:tdims%k_end)                      &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactivesolliquid]'//          newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactivesolrain (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactivesolrain]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactiveinsolice (tdims%k_start:tdims%k_end)                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactiveinsolice]'//           newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactivesolice (tdims%k_start:tdims%k_end)                         &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactivesolice]'//             newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactiveinsolliquid (tdims%k_start:tdims%k_end)                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactiveinsolliquid]'//        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactivesolnumber (tdims%k_start:tdims%k_end)                      &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactivesolnumber]'//          newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jactiveinsolnumber (tdims%k_start:tdims%k_end)                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jactiveinsolnumber]'//        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(je_trb (tdims%k_start:tdims%k_end)                                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [je_trb]'//                    newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtsq_trb (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtsq_trb]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jqsq_trb (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jqsq_trb]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcov_trb (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcov_trb]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_rad (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_rad]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_sw (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_sw]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_lw (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_lw]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_mic (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_mic]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_gwd (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_gwd]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_ph1 (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_ph1]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_conv (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_conv]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_bl (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_bl]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_stph (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_stph]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_cld (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_cdl]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_iau (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_iau]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_nud (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_nud]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_tot (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_tot]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_adv (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_adv]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_sol (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_sol]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdPV_mass (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdPV_mass]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jadv_only_PV (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jadv_only_PV]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_0 (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_0]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_bl (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_bl]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_bl_mix (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_bl_mix]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_bl_LH (tdims%k_start:tdims%k_end)                         &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_bl_LH]'//             newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_conv (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_conv]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_mic (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_mic]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_rad (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_rad]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_SW (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_SW]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_LW (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_LW]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_slow (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_slow]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdtheta_cld (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )

IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdtheta_cld]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jexner_rho_levels (pdims%k_start:pdims%k_end+1)                   &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jexner_rho_levels]'//         newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jconv_prog_1 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jconv_prog_1]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jconv_prog_2 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jconv_prog_2]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jconv_prog_3 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jconv_prog_3]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jconv_prog_precip (tdims%k_start:tdims%k_end)                     &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jconv_prog_precip]'//         newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(jp (pdims%k_start:pdims%k_end+1)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jp]'//                        newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jp_theta_levels (tdims%k_start:tdims%k_end)                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jp_theta_levels]'//           newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jexner_theta_levels (tdims%k_start:tdims%k_end)                   &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jexner_theta_levels]'//       newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(jccw_rad (tdims%k_end)                                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jccw_rad]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcca (n_cca_lev)                                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcca]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcca_dp (n_cca_lev)                                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcca_dp]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcca_md (n_cca_lev)                                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcca_md]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcca_sh (n_cca_lev)                                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcca_sh]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcf_area (tdims%k_end)                                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcf_area]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcf_bulk (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcf_bulk]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcf_liquid (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcf_liquid]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jcf_frozen (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jcf_frozen]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(j_deep_soil_temp (st_levels)                                      &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [j_deep_soil_temp]'//          newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsmcl (sm_levels)                                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsmcl]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsthu (sm_levels)                                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsthu]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsthf (sm_levels)                                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsthf]'//                     newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(jsw_incs (0:model_levels+1)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsw_incs]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jlw_incs (0:model_levels)                                         &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jlw_incs]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jozone (o3dims2%k_start:o3dims2%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jozone]'//                    newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtppsozone (tpps_ozone_levels)                                    &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtppsozone]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtracer (tdims_s%k_start:tdims_s%k_end,tr_vars+1)                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtracer]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtr_ukca (tdims_s%k_start:tdims_s%k_end,tr_ukca+1)                &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtr_ukca]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmurk_source (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmurk_source]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jmurk (tdims%k_start:tdims%k_end)                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jmurk]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jflash_pot (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jflash_pot]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jbl_w_var (tdims%k_end)                                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jbl_w_var]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div1 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div1]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div2 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div2]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div3 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div3]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div4 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div4]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div5 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div5]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdust_div6 (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdust_div6]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jso2 (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jso2]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jdms (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jdms]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jso4_aitken (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jso4_aitken]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jso4_accu (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jso4_accu]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jso4_diss (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jso4_diss]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jh2o2 (tdims%k_start:tdims%k_end)                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jh2o2]'//                     newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jnh3 (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jnh3]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsoot_new (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsoot_new]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsoot_agd (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsoot_agd]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jsoot_cld (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jsoot_cld]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jbmass_new (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jbmass_new]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jbmass_agd (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jbmass_agd]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jbmass_cld (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jbmass_cld]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jocff_new (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jocff_new]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jocff_agd (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jocff_agd]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jocff_cld (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jocff_cld]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jso2_natem (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jso2_natem]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(joh (tdims%k_start:tdims%k_end)                                   &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [joh]'//                       newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jho2 (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [joh2]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jh2o2_limit (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jh2o2_limit]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_chem (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_chem]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jco2 (tdims%k_start:tdims%k_end)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jco2]'//                      newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(joh_ukca (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [joh_ukca]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jho2_ukca (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jho2_ukca]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jh2o2_ukca (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jh2o2_ukca]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_ukca (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_ukca]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jhno3_ukca (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jhno3_ukca]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jozone_tracer (tdims%k_start:tdims%k_end)                         &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jozone_tracer]'//             newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_prod_loss (tdims%k_start:tdims%k_end)                         &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_prod_loss]'//             newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_p_l_vmr (tdims%k_start:tdims%k_end)                           &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_p_l_vmr]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_vmr (tdims%k_start:tdims%k_end)                               &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_vmr]'//                   newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_p_l_temp (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_p_l_temp]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_temp (tdims%k_start:tdims%k_end)                              &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_temp]'//                  newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_p_l_colo3 (tdims%k_start:tdims%k_end)                         &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_p_l_colo3]'//             newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jo3_colo3 (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jo3_colo3]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclbiog_bg (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclbiog_bg]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclbiom_fr (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclbiom_fr]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclbiom_ag (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclbiom_ag]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclbiom_ic (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclbiom_ic]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclblck_fr (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclblck_fr]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclblck_ag (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclblck_ag]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclsslt_fi (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclsslt_fi]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclsslt_jt (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclsslt_jt]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclsulp_ac (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclsulp_ac]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclsulp_ak (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclsulp_ak]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclsulp_di (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclsulp_di]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b1 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b1]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b2 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b2]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b3 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b3]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b4 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b4]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b5 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b5]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldust_b6 (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldust_b6]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclocff_fr (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclocff_fr]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclocff_ag (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclocff_ag]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarclocff_ic (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarclocff_ic]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jarcldlta_dl (tdims%k_start:tdims%k_end)                          &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jarcldlta_dl]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jnitr_acc (tdims%k_start:tdims%k_end)                             &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jnitr_acc]'//                 newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jnitr_diss (tdims%k_start:tdims%k_end)                            &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jnitr_diss]'//                newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nd_nuc_sol (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nd_nuc_sol]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nuc_sol_su (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nuc_sol_su]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nuc_sol_oc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nuc_sol_oc]'//            newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(jgc_nd_ait_sol (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nd_ait_sol]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_ait_sol_su (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_ait_sol_su]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_ait_sol_bc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_ait_sol_bc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_ait_sol_oc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_ait_sol_oc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nd_acc_sol (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nd_acc_sol]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_acc_sol_su (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_acc_sol_su]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_acc_sol_bc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_acc_sol_bc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_acc_sol_oc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_acc_sol_oc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_acc_sol_ss (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_acc_sol_ss]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nd_cor_sol (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nd_cor_sol]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_cor_sol_su (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_cor_sol_su]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_cor_sol_bc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_cor_sol_bc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_cor_sol_oc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_cor_sol_oc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_cor_sol_ss (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_cor_sol_ss]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_nd_ait_ins (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_nd_ait_ins]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_ait_ins_bc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_ait_ins_bc]'//            newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jgc_ait_ins_oc (tdims%k_start:tdims%k_end)                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jgc_ait_ins_oc]'//            newline//&
      TRIM(umMessage) )
END IF


ALLOCATE(juser_mult1 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult1]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult2 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult2]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult3 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult3]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult4 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult4]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult5 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult5]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult6 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult6]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult7 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult7]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult8 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult8]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult9 (model_levels)                                        &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult9]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult10 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult10]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult11 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult11]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult12 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult12]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult13 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult13]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult14 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult14]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult15 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult15]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult16 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult16]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult17 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult17]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult18 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult18]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult19 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult19]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(juser_mult20 (model_levels)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [juser_mult20]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtracer_lbc (tr_lbc_vars+1)                                       &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtracer_lbc]'//               newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtr_ukca_lbc (tr_lbc_ukca+1)                                      &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtr_ukca_lbc]'//              newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtracer_lbc_tend (tr_lbc_vars+1)                                  &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtracer_lbc_tend]'//          newline//&
      TRIM(umMessage) )
END IF

ALLOCATE(jtr_ukca_lbc_tend (tr_lbc_ukca+1)                                 &
         , STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode,                            &
      'Failure in allocating array [jtr_ukca_lbc_tend]'//         newline//&
      TRIM(umMessage) )
END IF

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END SUBROUTINE allocate_d1_indices

END MODULE atm_d1_indices_mod
