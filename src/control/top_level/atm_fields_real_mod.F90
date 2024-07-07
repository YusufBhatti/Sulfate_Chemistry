! DECLARE_ATM_FIELDS_REAL_MOD real arrays in atm_fields_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE atm_fields_real_mod

IMPLICIT NONE

!     TRACERS are a special case
      ! dummy array used to prevent crashing when there are no tracers
REAL, TARGET  :: dummy_field(1) = (/-1.0/)
REAL, TARGET  :: dummy_fld2(1,1) = -1.0
REAL, TARGET  :: dummy_fld3(1,1,1) = -1.0
REAL, TARGET  :: dummy_fld4(1,1,1,1) = -1.0
REAL, POINTER :: tracer(:,:,:,:)
REAL, POINTER :: tracer_ukca(:,:,:,:)
REAL, POINTER :: tracer_ukca_lbc(:,:,:)
REAL, POINTER :: tracer_ukca_lbc_tend(:,:,:)
REAL, POINTER :: tracer_lbc(:,:,:)
REAL, POINTER :: tracer_lbc_tend(:,:,:)
!     Add a specific pointer for ozone tracer and cariolle parameters
REAL, POINTER :: ozone_tracer(:,:,:)
REAL, POINTER :: o3_prod_loss(:,:)
REAL, POINTER :: o3_p_l_vmr(:,:)
REAL, POINTER :: o3_vmr(:,:)
REAL, POINTER :: o3_p_l_temp(:,:)
REAL, POINTER :: o3_temp(:,:)
REAL, POINTER :: o3_p_l_colo3(:,:)
REAL, POINTER :: o3_colo3(:,:)
! 1: Array Variables (dimensions are resolution dependent.)

! Oxidant concentrations from UKCA for use in HadGEM sulphur cycle
REAL, POINTER :: oh_ukca(:,:,:)
REAL, POINTER :: h2o2_ukca(:,:,:)
REAL, POINTER :: ho2_ukca(:,:,:)
REAL, POINTER :: o3_ukca(:,:,:)
REAL, POINTER :: hno3_ukca(:,:,:)

! 1.1: Data variables stored in primary space.

REAL, POINTER :: sst(:,:)

REAL, POINTER :: DryRho(:,:,:)
REAL, POINTER :: EtaDot(:,:,:)
REAL, POINTER :: ThetaV(:,:,:)
REAL, POINTER :: psi_w_surf(:,:)
REAL, POINTER :: psi_w_lid(:,:)
REAL, POINTER :: m_v(:,:,:)
REAL, POINTER :: m_cl(:,:,:)
REAL, POINTER :: m_cf(:,:,:)
REAL, POINTER :: m_cf2(:,:,:)
REAL, POINTER :: m_gr(:,:,:)
REAL, POINTER :: m_r(:,:,:)
REAL, POINTER :: exner_surf(:,:)
REAL, POINTER :: exner(:,:,:)

REAL, POINTER :: u(:,:,:)      ! u component of wind
REAL, POINTER :: v(:,:,:)      ! v component of wind
REAL, POINTER :: w(:,:,:)      ! w component of wind
REAL, POINTER :: wetrho_r_sq_n(:,:,:)    ! Density
REAL, POINTER :: rho(:,:,:)    ! Density
REAL, POINTER :: theta(:,:,:)  ! Potential temperature
REAL, POINTER :: q(:,:,:)        ! Specific humidity
REAL, POINTER :: qcl(:,:,:)      ! qcl
REAL, POINTER :: qcf(:,:,:)      ! qcf
REAL, POINTER :: qcf2(:,:,:)     ! second ice
REAL, POINTER :: qrain(:,:,:)    ! rain
REAL, POINTER :: qgraup(:,:,:)   ! graupel

! Atmospheric Electricity / Lightning scheme
REAL, POINTER :: flash_pot(:,:,:) ! lightning flash potential

! Boundary layer scheme w-variance for turbulent production
! of mixed phase clouds
REAL, POINTER :: bl_w_var(:,:,:) ! bl w-variance

! TKE based turbulence scheme
REAL, POINTER :: e_trb(:,:,:)     ! TKE
REAL, POINTER :: tsq_trb(:,:,:)   ! Self covariance of thetal'
REAL, POINTER :: qsq_trb(:,:,:)   ! Self coveriance of qw'
REAL, POINTER :: cov_trb(:,:,:)   ! Correlation between thetal' and qw'
REAL, POINTER :: zhpar_shcu(:,:)  ! Height of mixed layer used to
                                  ! evaluate the non-grad buoy flux

!PV-tracers
REAL, POINTER :: dPV_rad(:,:,:)
REAL, POINTER :: dPV_sw(:,:,:)
REAL, POINTER :: dPV_lw(:,:,:)
REAL, POINTER :: dPV_mic(:,:,:)
REAL, POINTER :: dPV_gwd(:,:,:)
REAL, POINTER :: dPV_ph1(:,:,:)
REAL, POINTER :: dPV_conv(:,:,:)
REAL, POINTER :: dPV_bl(:,:,:)
REAL, POINTER :: dPV_stph(:,:,:)
REAL, POINTER :: dPV_cld(:,:,:)
REAL, POINTER :: dPV_iau(:,:,:)
REAL, POINTER :: dPV_nud(:,:,:)
REAL, POINTER :: dPV_tot(:,:,:)
REAL, POINTER :: dPV_adv(:,:,:)
REAL, POINTER :: dPV_sol(:,:,:)
REAL, POINTER :: dPV_mass(:,:,:)
REAL, POINTER :: adv_only_PV(:,:,:)

!Diabatic tracers
REAL, POINTER :: dtheta_0(:,:,:)
REAL, POINTER :: dtheta_bl(:,:,:)
REAL, POINTER :: dtheta_bl_mix(:,:,:)
REAL, POINTER :: dtheta_bl_LH(:,:,:)
REAL, POINTER :: dtheta_conv(:,:,:)
REAL, POINTER :: dtheta_mic(:,:,:)
REAL, POINTER :: dtheta_rad(:,:,:)
REAL, POINTER :: dtheta_SW(:,:,:)
REAL, POINTER :: dtheta_LW(:,:,:)
REAL, POINTER :: dtheta_slow(:,:,:)
REAL, POINTER :: dtheta_cld(:,:,:)

! Exner pressure on rho levels
REAL, POINTER :: exner_rho_levels(:,:,:)

REAL, POINTER :: u_adv(:,:,:) ! Advective u component of wind
REAL, POINTER :: v_adv(:,:,:) ! Advective v component of wind
REAL, POINTER :: w_adv(:,:,:) ! Advective w component of wind
REAL, TARGET, ALLOCATABLE :: U_ADV_nodump(:)
REAL, TARGET, ALLOCATABLE :: V_ADV_nodump(:)
REAL, TARGET, ALLOCATABLE :: W_ADV_nodump(:)

! Convective prognostics
REAL, POINTER :: conv_prog_1(:,:,:)      ! Convection experimental prognostic 1
REAL, POINTER :: conv_prog_2(:,:,:)      ! Convection experimental prognostic 2
REAL, POINTER :: conv_prog_3(:,:,:)      ! Convection experimental prognostic 3
REAL, POINTER :: conv_prog_precip(:,:,:) ! Convection prognostic recent precip

! Total precipitation at start of timestep (for input to gw_ussp source calc.)
REAL, POINTER :: totalppn(:,:)

! Stochastic physics fields for BL perturbations
REAL, POINTER :: bl_pert_rand_fld(:,:)
REAL, POINTER :: bl_pert_flag(:,:)

! 1.2: Data variables stored in secondary space.
REAL, POINTER :: p(:,:,:)                  ! Pressure on rho levels

! Pressure on theta levels
REAL, POINTER :: p_theta_levels(:,:,:)

! Exner pressure on theta levels
REAL, POINTER :: exner_theta_levels(:,:,:)

! 1.3: Cloud Fields
REAL, POINTER :: cca(:,:,:)                ! Convective cloud amount
REAL, POINTER :: cf_area(:,:,:)            ! Area Cloud Fraction
REAL, POINTER :: cf_bulk(:,:,:)            ! Bulk Cloud Fraction
REAL, POINTER :: cf_liquid(:,:,:)          ! Liquid cloud fraction
REAL, POINTER :: cf_frozen(:,:,:)          ! Frozen cloud fraction

! 1.4: Soil Ancillary fields
REAL, POINTER :: deep_soil_temp(:,:)   ! Deep soil temperature

REAL, POINTER :: smcl(:,:)   ! Soil moisture content in layers
REAL, POINTER :: sthu(:,:)   ! Unfrozen soil moisture fraction
REAL, POINTER :: sthf(:,:)   ! Frozen soil moisture fraction

! 1.5: Radiation Increments
REAL, POINTER :: sw_incs(:,:,:)    ! SW radiation increments
REAL, POINTER :: lw_incs(:,:,:)    ! LW radiation increments
! PAR radiation increment
REAL, POINTER :: dirpar(:,:)

! 1.6: Ozone
REAL, POINTER :: o3(:)          ! Ozone
!  tropopause-based ozone
REAL, POINTER :: tppsozone(:)
! 1.7: Tracer and aerosol fields
REAL, POINTER :: murk_source(:,:,:)    ! multi-level murk source
REAL, POINTER :: murk(:,:,:)           ! multi-level murk concent
REAL, POINTER :: dust_div1(:,:,:)      ! dust mmr, division 1
REAL, POINTER :: dust_div2(:,:,:)      ! dust mmr, division 2
REAL, POINTER :: dust_div3(:,:,:)      ! dust mmr, division 3
REAL, POINTER :: dust_div4(:,:,:)      ! dust mmr, division 4
REAL, POINTER :: dust_div5(:,:,:)      ! dust mmr, division 5
REAL, POINTER :: dust_div6(:,:,:)      ! dust mmr, division 6
REAL, POINTER :: so2(:,:,:)            ! sulphur dioxide gas
REAL, POINTER :: dms(:,:,:)            ! dimethyl sulphide gas
REAL, POINTER :: so4_aitken(:,:,:)     ! Aitken mode sulphate aer
REAL, POINTER :: so4_accu(:,:,:)       ! accumulation mode sulpha
REAL, POINTER :: so4_diss(:,:,:)       ! dissloved  sulphate aero
REAL, POINTER :: h2o2(:,:,:)           ! hydrogen peroxide mmr
REAL, POINTER :: nh3(:,:,:)            ! ammonia gas mmr
REAL, POINTER :: soot_new(:,:,:)       ! fresh soot mmr
REAL, POINTER :: soot_agd(:,:,:)       ! aged soot mmr
REAL, POINTER :: soot_cld(:,:,:)       ! soot in cloud mmr
REAL, POINTER :: bmass_new(:,:,:)      ! fresh biomass mmr
REAL, POINTER :: bmass_agd(:,:,:)      ! aged biomass mmr
REAL, POINTER :: bmass_cld(:,:,:)      ! cloud biomass mmr
REAL, POINTER :: ocff_new(:,:,:)       ! fresh OCFF mmr
REAL, POINTER :: ocff_agd(:,:,:)       ! aged OCFF mmr
REAL, POINTER :: ocff_cld(:,:,:)       ! OCFF in cloud mmr
REAL, POINTER :: so2_natem(:,:,:)      ! natural SO2 emissions
REAL, POINTER :: oh(:,:,:)             ! hydroxyl radical ancilla
REAL, POINTER :: ho2(:,:,:)            ! hydrogen dioxide ancilla
REAL, POINTER :: h2o2_limit(:,:,:)     ! limiting H2O2 ancillary
REAL, POINTER :: o3_chem(:,:,:)        ! ozone for chemistry anci
REAL, POINTER :: nitr_acc(:,:,:)       ! accumulation ammonium nitrate
REAL, POINTER :: nitr_diss(:,:,:)      ! dissolved ammonium nitrate
REAL, POINTER :: co2(:,:,:)            ! 3D CO2 FIELD
! 1.7a: GLOMAP_CLIM aerosol section 54
REAL, POINTER :: gc_nd_nuc_sol(:,:,:)  ! GLOMAP_CLIM NUC (Sol) number
REAL, POINTER :: gc_nuc_sol_su(:,:,:)  ! GLOMAP_CLIM NUC (Sol) SO4
REAL, POINTER :: gc_nuc_sol_oc(:,:,:)  ! GLOMAP_CLIM NUC (Sol) OC
REAL, POINTER :: gc_nd_ait_sol(:,:,:)  ! GLOMAP_CLIM AIT (Sol) number
REAL, POINTER :: gc_ait_sol_su(:,:,:)  ! GLOMAP_CLIM AIT (Sol) SO4
REAL, POINTER :: gc_ait_sol_bc(:,:,:)  ! GLOMAP_CLIM AIT (Sol) BC
REAL, POINTER :: gc_ait_sol_oc(:,:,:)  ! GLOMAP_CLIM AIT (Sol) OC
REAL, POINTER :: gc_nd_acc_sol(:,:,:)  ! GLOMAP_CLIM ACC (Sol) number
REAL, POINTER :: gc_acc_sol_su(:,:,:)  ! GLOMAP_CLIM ACC (Sol) SO4
REAL, POINTER :: gc_acc_sol_bc(:,:,:)  ! GLOMAP_CLIM ACC (Sol) BC
REAL, POINTER :: gc_acc_sol_oc(:,:,:)  ! GLOMAP_CLIM ACC (Sol) OC
REAL, POINTER :: gc_acc_sol_ss(:,:,:)  ! GLOMAP_CLIM ACC (Sol) SS
REAL, POINTER :: gc_nd_cor_sol(:,:,:)  ! GLOMAP_CLIM COR (Sol) number
REAL, POINTER :: gc_cor_sol_su(:,:,:)  ! GLOMAP_CLIM COR (Sol) SO4
REAL, POINTER :: gc_cor_sol_bc(:,:,:)  ! GLOMAP_CLIM COR (Sol) BC
REAL, POINTER :: gc_cor_sol_oc(:,:,:)  ! GLOMAP_CLIM COR (Sol) OC
REAL, POINTER :: gc_cor_sol_ss(:,:,:)  ! GLOMAP_CLIM COR (Sol) SS
REAL, POINTER :: gc_nd_ait_ins(:,:,:)  ! GLOMAP_CLIM AIT (Ins) number
REAL, POINTER :: gc_ait_ins_bc(:,:,:)  ! GLOMAP_CLIM AIT (Ins) BC
REAL, POINTER :: gc_ait_ins_oc(:,:,:)  ! GLOMAP_CLIM AIT (Ins) OC
! 1.8: Multi-level user ancillary fields
REAL, POINTER :: user_mult1(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult2(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult3(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult4(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult5(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult6(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult7(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult8(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult9(:)     ! multi-level user ancilla
REAL, POINTER :: user_mult10(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult11(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult12(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult13(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult14(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult15(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult16(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult17(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult18(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult19(:)    ! multi-level user ancilla
REAL, POINTER :: user_mult20(:)    ! multi-level user ancilla

! 1.9: Fields carried forward from previous version.
! Lateral Boundary Conditions
REAL, POINTER :: orog_lbc(:)                         ! Orography LBC
REAL, POINTER :: u_lbc(:,:)                          ! U LBC
REAL, POINTER :: v_lbc(:,:)                          ! V LBC
REAL, POINTER :: w_lbc(:,:)                          ! W LBC
REAL, POINTER :: rho_lbc(:,:)                        ! RHO LBC
REAL, POINTER :: theta_lbc(:,:)                      ! Theta LBC
REAL, POINTER :: q_lbc(:,:)                          ! Q LBC
REAL, POINTER :: qcl_lbc(:,:)                        ! QCL LBC
REAL, POINTER :: qcf_lbc(:,:)                        ! QCF LBC
REAL, POINTER :: qcf2_lbc(:,:)                       ! 2nd Ice LBC
REAL, POINTER :: qrain_lbc(:,:)                      ! Rain LBC
REAL, POINTER :: qgraup_lbc(:,:)                     ! Graupel LBC
REAL, POINTER :: cf_bulk_lbc(:,:)                    ! CF_BULK LBC
REAL, POINTER :: cf_liquid_lbc(:,:)                  ! CF_LIQUID_LBC
REAL, POINTER :: cf_frozen_lbc(:,:)                  ! CF_FROZEN_LBC
REAL, POINTER :: exner_lbc(:,:)                      ! Exner LBC
REAL, POINTER :: u_adv_lbc(:,:)                      ! U_ADV LBC
REAL, POINTER :: v_adv_lbc(:,:)                      ! V_ADV LBC
REAL, POINTER :: w_adv_lbc(:,:)                      ! W_ADV LBC
REAL, POINTER :: murk_lbc(:,:)                       ! Murk aerosol LBC
REAL, POINTER :: dust_div1_lbc(:,:)      ! dust mmr, division 1 LBC
REAL, POINTER :: dust_div2_lbc(:,:)      ! dust mmr, division 2 LBC
REAL, POINTER :: dust_div3_lbc(:,:)      ! dust mmr, division 3 LBC
REAL, POINTER :: dust_div4_lbc(:,:)      ! dust mmr, division 4 LBC
REAL, POINTER :: dust_div5_lbc(:,:)      ! dust mmr, division 5 LBC
REAL, POINTER :: dust_div6_lbc(:,:)      ! dust mmr, division 6 LBC
REAL, POINTER :: so2_lbc(:,:)            ! sulphur dioxide gas LBC
REAL, POINTER :: dms_lbc(:,:)            ! dimethyl sulphide gas LBC
REAL, POINTER :: so4_aitken_lbc(:,:)     ! Aitken mode sulphate aerosol LBC
REAL, POINTER :: so4_accu_lbc(:,:)       ! accumulation mode sulphate LBC
REAL, POINTER :: so4_diss_lbc(:,:)       ! dissolved sulphate aero LBC
REAL, POINTER :: nh3_lbc(:,:)            ! ammonia gas LBC
REAL, POINTER :: soot_new_lbc(:,:)       ! fresh soot LBC
REAL, POINTER :: soot_agd_lbc(:,:)       ! aged soot LBC
REAL, POINTER :: soot_cld_lbc(:,:)       ! soot in cloud LBC
REAL, POINTER :: bmass_new_lbc(:,:)      ! fresh biomass LBC
REAL, POINTER :: bmass_agd_lbc(:,:)      ! aged biomass LBC
REAL, POINTER :: bmass_cld_lbc(:,:)      ! cloud biomass LBC
REAL, POINTER :: ocff_new_lbc(:,:)       ! fresh fossil fuel LBC
REAL, POINTER :: ocff_agd_lbc(:,:)       ! aged fossil fuel LBC
REAL, POINTER :: ocff_cld_lbc(:,:)       ! cloud fossil fuel LBC
REAL, POINTER :: nitr_acc_lbc(:,:)       ! accumulation ammonium nitrate LBC
REAL, POINTER :: nitr_diss_lbc(:,:)      ! dissolved ammonium nitrate LBC

! Lateral Boundary Condition tendencies
REAL, POINTER :: u_lbc_tend(:,:)                     ! U LBC  tendencies
REAL, POINTER :: v_lbc_tend(:,:)                     ! V LBC tendencies
REAL, POINTER :: w_lbc_tend(:,:)                     ! W LBC tendencies
REAL, POINTER :: rho_lbc_tend(:,:)                   ! RHO LBC tendencies
REAL, POINTER :: theta_lbc_tend(:,:)                 ! Theta LBC tendencies
REAL, POINTER :: q_lbc_tend(:,:)                     ! Q LBC tendencies
REAL, POINTER :: qcl_lbc_tend(:,:)                   ! QCL LBC tendencies
REAL, POINTER :: qcf_lbc_tend(:,:)                   ! QCF LBC tendencies
REAL, POINTER :: qcf2_lbc_tend(:,:)                  ! 2nd Ice
REAL, POINTER :: qrain_lbc_tend(:,:)                 ! Rain LBC tendencies
REAL, POINTER :: qgraup_lbc_tend(:,:)                ! Graupel
REAL, POINTER :: cf_bulk_lbc_tend(:,:)               ! CF_BULK LBC tend'cies
REAL, POINTER :: cf_liquid_lbc_tend(:,:)             ! CF_LIQUID_LBC t'cies
REAL, POINTER :: cf_frozen_lbc_tend(:,:)             ! CF_FROZEN_LBC t'cies
REAL, POINTER :: exner_lbc_tend(:,:)                 ! Exner LBC tendencies
REAL, POINTER :: u_adv_lbc_tend(:,:)                 ! U_ADV LBC tendencies
REAL, POINTER :: v_adv_lbc_tend(:,:)                 ! V_ADV LBC tendencies
REAL, POINTER :: w_adv_lbc_tend(:,:)                 ! W_ADV LBCtendencies
REAL, POINTER :: murk_lbc_tend(:,:)                  ! Murk aerosol LBC tend
REAL, POINTER :: dust_div1_lbc_tend(:,:)      ! dust mmr, division 1 LBC tend
REAL, POINTER :: dust_div2_lbc_tend(:,:)      ! dust mmr, division 2 LBC tend
REAL, POINTER :: dust_div3_lbc_tend(:,:)      ! dust mmr, division 3 LBC tend
REAL, POINTER :: dust_div4_lbc_tend(:,:)      ! dust mmr, division 4 LBC tend
REAL, POINTER :: dust_div5_lbc_tend(:,:)      ! dust mmr, division 5 LBC tend
REAL, POINTER :: dust_div6_lbc_tend(:,:)      ! dust mmr, division 6 LBC tend
REAL, POINTER :: so2_lbc_tend(:,:)            ! sulphur dioxide gas LBC tend
REAL, POINTER :: dms_lbc_tend(:,:)            ! dimethyl sulphide gas LBC tend
REAL, POINTER :: so4_aitken_lbc_tend(:,:)     ! Aitken mode sulphate aerosol LBC
REAL, POINTER :: so4_accu_lbc_tend(:,:)       ! accumulation mode sulphate LBC
REAL, POINTER :: so4_diss_lbc_tend(:,:)       ! dissolved sulphate aero LBC
REAL, POINTER :: nh3_lbc_tend(:,:)            ! ammonia gas LBC tend
REAL, POINTER :: soot_new_lbc_tend(:,:)       ! fresh soot LBC tend
REAL, POINTER :: soot_agd_lbc_tend(:,:)       ! aged soot LBC tend
REAL, POINTER :: soot_cld_lbc_tend(:,:)       ! soot in cloud LBC tend
REAL, POINTER :: bmass_new_lbc_tend(:,:)      ! fresh biomass LBC tend
REAL, POINTER :: bmass_agd_lbc_tend(:,:)      ! aged biomass LBC tend
REAL, POINTER :: bmass_cld_lbc_tend(:,:)      ! cloud biomass LBC tend
REAL, POINTER :: ocff_new_lbc_tend(:,:)       ! fresh fossil fuel LBC tend
REAL, POINTER :: ocff_agd_lbc_tend(:,:)       ! aged fossil fuel LBC tend
REAL, POINTER :: ocff_cld_lbc_tend(:,:)       ! cloud fossil fuel LBC tend
REAL, POINTER :: nitr_acc_lbc_tend(:,:)       ! accm'n ammonium nitrate LBC tend
REAL, POINTER :: nitr_diss_lbc_tend(:,:)      ! dissolved amm nitrate LBC tend

! 2: Scalar Variables

! 2.1: Data variables stored in primary space.
REAL, POINTER :: tstar(:,:)            ! Surface temperature
REAL, POINTER :: tstar_anom(:)         ! Surface temperature anomaly
!   2.15: Fields for coastal tiling
REAL, POINTER :: frac_land(:)          ! Land fraction in grid box
REAL, POINTER :: tstar_land(:,:)       ! Land surface temperature
REAL, POINTER :: tstar_sea(:,:)        ! Sea surface temperature
REAL, POINTER :: tstar_sice(:,:,:)     ! Sea-ice surface temperature
REAL, POINTER :: tstar_sice_cat(:,:,:) ! Sea-ice category surface temp

!   Set pointer for sea surface freezing temperature 
REAL, POINTER :: sstfrz(:,:)  

! Set pointers for sea-ice and land albedos
REAL, POINTER :: sice_alb(:,:)         ! Sea-ice albedo
REAL, POINTER :: land_alb(:,:)         ! Mean land albedo

! 2.2: Data variables stored in secondary space.

REAL, POINTER :: pstar(:,:)          ! Surface pressure

! 2.3: Cloud fields

REAL, POINTER :: cclwp(:,:)          ! Convective cloud liquid water path
REAL, POINTER :: deep_flag(:,:)      ! Deep convection flag
REAL, POINTER :: past_precip(:,:)    ! Past convective precipitation
REAL, POINTER :: past_conv_ht(:,:)   ! Past convective height

! 2.4: Boundary layer fields

REAL, POINTER :: zh(:,:)             ! Boundary layer depth

REAL, POINTER :: ddmfx(:,:)          ! Convective dowdraught
                                     ! mass-flux at cloud base

! Standard deviation of turbulent fluctuations of layer 1 temperature
REAL, POINTER :: t1_sd(:,:)

! Standard deviation of turbulent fluctuations of layer 1 humidity
REAL, POINTER :: q1_sd(:,:)

! Decoupled screen-level temperature
REAL, POINTER :: TScrnDcl_SSI(:,:)
REAL, POINTER :: TScrnDcl_TILE(:,:)
REAL, POINTER :: tStbTrans(:,:)

! Convective cold pools
REAL, POINTER :: ux_ccp  (:,:)
REAL, POINTER :: uy_ccp  (:,:)
REAL, POINTER :: um_ccp  (:,:)
REAL, POINTER :: g_ccp   (:,:)
REAL, POINTER :: h_ccp   (:,:)
REAL, POINTER :: riso_ccp(:,:)
REAL, POINTER :: rdir_ccp(:,:)

! 2.4: Soil Ancillary fields (ALL land_field, single level)

REAL, POINTER :: sat_soilw_suction(:) ! Saturated soil water suction
REAL, POINTER :: therm_cap    (:)     ! Thermal capacity
REAL, POINTER :: therm_cond   (:)     ! Thermal conductivity
REAL, POINTER :: vol_smc_crit (:)     ! Vol smc at critical point
REAL, POINTER :: vol_smc_wilt (:)     ! Vol smc at wilting point
REAL, POINTER :: vol_smc_sat  (:)     ! Vol smc at saturation
REAL, POINTER :: sat_soil_cond(:)     ! Saturated soil conductivity
REAL, POINTER :: clapp_horn   (:)     ! Clapp-Hornberger B coefficient
REAL, POINTER :: z0m_soil     (:)     ! Bare soil momentum roughness length

! 2.5: Other surface fields
REAL, POINTER :: canopy_water(:) ! Canopy Water                    (l_f, 1 lev)
REAL, POINTER :: gs(:)           ! Gridbox mean canopy conductance (l_f, 1 lev)
REAL, POINTER :: z0(:,:)         ! Roughness length, used for sea points

! 2.6: Orographic Ancillary fields

REAL, POINTER :: orography(:,:)     ! Orographic height
REAL, POINTER :: orog_sd(:)       ! Standard Deviation of orography
REAL, POINTER :: orog_sil(:)      ! Silhouette area of orography
REAL, POINTER :: orog_ho2(:)      ! Peak to trough height/(2*sqrt2)
REAL, POINTER :: orog_grad_x(:)
REAL, POINTER :: orog_grad_y(:)
REAL, POINTER :: orog_unfilt(:)
REAL, POINTER :: orog_grad_xx(:)  ! Orographic gradient xx
REAL, POINTER :: orog_grad_xy(:)  ! Orographic gradient xy
REAL, POINTER :: orog_grad_yy(:)  ! Orographic gradient yy

! 2.7: Sea/Sea Ice

REAL, POINTER :: u_sea(:,:)           ! Surface current (u component)
REAL, POINTER :: v_sea(:,:)           ! Surface current (v component)
REAL, POINTER :: u_0_p(:,:)           ! Surace  current (u) on p-grid
REAL, POINTER :: v_0_p(:,:)           ! Surface current (v) on p-grid
REAL, POINTER :: ice_fract_cat(:,:,:) ! Sea ice fraction on categories
REAL, POINTER :: ice_thick_cat(:,:,:) ! Sea ice thickness on categories
REAL, POINTER :: ti_cat(:,:,:)        ! Sea ice temperature on categories
REAL, POINTER :: ice_k_cat(:,:,:)     ! Sea ice effect cond on categories
REAL, POINTER :: chloro_sea(:,:)      ! Sea near surface sea chlorophyll
REAL, POINTER :: ice_fraction(:,:,:)  ! Sea ice fraction
REAL, POINTER :: ice_thickness(:,:,:) ! Sea ice depth
REAL, POINTER :: ti(:,:,:)            ! Sea ice temperature
REAL, POINTER :: pond_depth_cat(:,:,:)  ! Meltpond depth on categories 
REAL, POINTER :: pond_frac_cat(:,:,:)   ! Meltpond fraction on categories 

! 2.8: Snow

REAL, POINTER :: snodep(:,:)      ! Snow depth on land
REAL, POINTER :: snodep_sea(:,:,:)! Snow depth on sea ice (theta_pts_sea_only)
REAL, POINTER :: snodep_sea_cat(:,:,:) ! Snow depth on sea ice catagories
REAL, POINTER :: catch_snow(:,:)       ! Coniferous canopy snow capacity
REAL, POINTER :: snow_grnd(:,:)   ! Snow below canopy
REAL, POINTER :: snsoot(:,:)      ! Snow soot content

! 2.9: aerosol emission fields, including mineral dust parent soil props

REAL, POINTER :: soil_clay(:,:)                    ! soil clay fraction
REAL, POINTER :: soil_silt(:,:)                    ! soil silt fraction
REAL, POINTER :: soil_sand(:,:)                    ! soil sand fraction
REAL, POINTER :: dust_mrel1(:,:)                   ! soil rel mass, div 1
REAL, POINTER :: dust_mrel2(:,:)                   ! soil rel mass, div 2
REAL, POINTER :: dust_mrel3(:,:)                   ! soil rel mass, div 3
REAL, POINTER :: dust_mrel4(:,:)                   ! soil rel mass, div 4
REAL, POINTER :: dust_mrel5(:,:)                   ! soil rel mass, div 5
REAL, POINTER :: dust_mrel6(:,:)                   ! soil rel mass, div 6


REAL, POINTER :: so2_em(:,:)        ! sulphur dioxide emission
REAL, POINTER :: dms_em(:,:)        ! dimethyl sulphide emission
REAL, POINTER :: so2_hilem(:,:)     ! high level SO2 emissions
REAL, POINTER :: nh3_em(:,:)        ! ammonia gas surface emiss
REAL, POINTER :: soot_em(:,:)       ! fresh soot surface emissions
REAL, POINTER :: soot_hilem(:,:)    ! fresh soot high lev emissions
REAL, POINTER :: bmass_em(:,:)      ! fresh bmass surface emissions
REAL, POINTER :: bmass_hilem(:,:)   ! fresh bmass high lev emissions
REAL, POINTER :: ocff_em(:,:)       ! fresh OCFF surface emissions
REAL, POINTER :: ocff_hilem(:,:)    ! fresh OCFF high lev emissions
REAL, POINTER :: dms_conc(:,:)      ! seawater dimethyl sulphide conc.
! Minimum and maximum heights for injection of fresh bmass high lev emissions
REAL, POINTER :: bmass_hilem_h1 (:,:)
REAL, POINTER :: bmass_hilem_h2 (:,:)

! 2.10: User ancillary fields
REAL, POINTER :: user_anc1(:)         ! user ancillary field 1
REAL, POINTER :: user_anc2(:)         ! user ancillary field 2
REAL, POINTER :: user_anc3(:)         ! user ancillary field 3
REAL, POINTER :: user_anc4(:)         ! user ancillary field 4
REAL, POINTER :: user_anc5(:)         ! user ancillary field 5
REAL, POINTER :: user_anc6(:)         ! user ancillary field 6
REAL, POINTER :: user_anc7(:)         ! user ancillary field 7
REAL, POINTER :: user_anc8(:)         ! user ancillary field 8
REAL, POINTER :: user_anc9(:)         ! user ancillary field 9
REAL, POINTER :: user_anc10(:)        ! user ancillary field 10
REAL, POINTER :: user_anc11(:)        ! user ancillary field 11
REAL, POINTER :: user_anc12(:)        ! user ancillary field 12
REAL, POINTER :: user_anc13(:)        ! user ancillary field 13
REAL, POINTER :: user_anc14(:)        ! user ancillary field 14
REAL, POINTER :: user_anc15(:)        ! user ancillary field 15
REAL, POINTER :: user_anc16(:)        ! user ancillary field 16
REAL, POINTER :: user_anc17(:)        ! user ancillary field 17
REAL, POINTER :: user_anc18(:)        ! user ancillary field 18
REAL, POINTER :: user_anc19(:)        ! user ancillary field 19
REAL, POINTER :: user_anc20(:)        ! user ancillary field 20

!   2.11: Store arrays for energy correction calculation
REAL, POINTER :: net_flux(:,:)                   ! Net energy flux
REAL, POINTER :: net_mflux(:,:)                  ! Net moisture flux

!   2.12: Tiled Vegetation and Triffid fields
REAL, POINTER :: frac_typ(:,:)      ! Fractions of surface type
REAL, POINTER :: frac_con1(:)       ! Fractions of surface type
REAL, POINTER :: frac_con2(:)       ! Fractions of surface type
REAL, POINTER :: frac_con3(:)       ! Fractions of surface type
REAL, POINTER :: frac_con4(:)       ! Fractions of surface type
REAL, POINTER :: frac_con5(:)       ! Fractions of surface type
REAL, POINTER :: frac_con6(:)       ! Fractions of surface type
REAL, POINTER :: frac_con7(:)       ! Fractions of surface type
REAL, POINTER :: frac_con8(:)       ! Fractions of surface type
REAL, POINTER :: frac_con9(:)       ! Fractions of surface type

REAL, POINTER :: lai_pft(:,:)       ! LAI of plant functional types
REAL, POINTER :: canht_pft(:,:)     ! Canopy hght of plant func types
REAL, POINTER :: disturb_veg(:)     ! Disturbed fraction of vegetation (1-lev)
REAL, POINTER :: disturb_veg_prev(:) ! Previous Disturbed fraction of vegetation
REAL, POINTER :: pasture_frac_d1(:)      ! Pasture fraction of vegetation
REAL, POINTER :: pasture_frac_prev_d1(:) ! Previous pasture fraction
REAL, POINTER :: agr_crop_frac_d1(:)     ! Crop fraction of vegetation
REAL, POINTER :: agr_crop_frac_prev_d1(:)! Previous crop fraction
REAL, POINTER :: wood_prod_fast_d1(:)  ! Wood product pool (fast)  
REAL, POINTER :: wood_prod_med_d1(:)   ! Wood product pool (medium)
REAL, POINTER :: wood_prod_slow_d1(:)  ! Wood product pool (slow)  

REAL, POINTER :: soil_alb(:)        ! Snow-free albedo of bare soil    (1-lev)

REAL, POINTER :: obs_alb_sw(:)      ! Observed snow-free SW albedo
REAL, POINTER :: obs_alb_vis(:)     ! Observed snow-free VIS albedo
REAL, POINTER :: obs_alb_nir(:)     ! Observed snow-free NIR albedo

REAL, POINTER :: soil_nitro1(:)     ! Soil organic nitrogen content DPM
REAL, POINTER :: soil_nitro2(:)     ! Soil organic nitrogen content RPM
REAL, POINTER :: soil_nitro3(:)     ! Soil organic nitrogen content BIO
REAL, POINTER :: soil_nitro4(:)     ! Soil organic  nitrogen content HUM

REAL, POINTER :: soil_inorgnit(:)   ! Soil inorganic nitrogen content 
REAL, POINTER :: nitrogen_deposition_d1(:) ! Nitrogen deposition on land
 
REAL, POINTER :: soil_carb(:,:)     ! Soil carbon content
REAL, POINTER :: soil_carb1(:,:)    ! Soil carbon content DPM
REAL, POINTER :: soil_carb2(:,:)      ! Soil carbon content RPM
REAL, POINTER :: soil_carb3(:,:)      ! Soil carbon content BIO
REAL, POINTER :: soil_carb4(:,:)      ! Soil carbon content HUM

REAL, POINTER :: npp_pft_acc(:,:)     ! Accumulated NPP on PFTs
REAL, POINTER :: g_lf_pft_acc(:,:)    ! Accum. leaf turnover rate PFTs
REAL, POINTER :: g_phlf_pft_acc(:,:)  ! Accumulated phenological leaf
                                      !           turnover rate on PFTs
REAL, POINTER :: rsp_w_pft_acc(:,:)   ! Accum. wood respiration on PFTs
REAL, POINTER :: rsp_s_acc(:,:)       ! Accumulated soil respiration
REAL, POINTER :: rsp_s_acc1(:,:)      ! Accumulated soil respiration DPM
REAL, POINTER :: rsp_s_acc2(:,:)      ! Accumulated soil respiration RPM
REAL, POINTER :: rsp_s_acc3(:,:)      ! Accumulated soil respiration BIO
REAL, POINTER :: rsp_s_acc4(:,:)      ! Accumulated soil respiration HUM

REAL, POINTER :: can_water_tile(:,:)  ! Canopy water content on tiles
REAL, POINTER :: catch_tile(:,:)      ! Canopy capacity on tiles
REAL, POINTER :: infil_tile(:,:)      ! Max infiltration rate on tiles
REAL, POINTER :: rgrain_tile(:,:)     ! Snow grain size on tiles
REAL, POINTER :: snodep_tile(:,:)     ! Snow depth on tiles
REAL, POINTER :: tstar_tile(:,:)      ! Surface temperature on tiles
REAL, POINTER :: tsurf_elev_surft(:,:)! Temperature of elevated  
                                      ! subsurface tiles (K)
REAL, POINTER :: z0_tile(:,:)         ! Surface roughness on tiles
REAL, POINTER :: z0h_tile(:,:)        ! Snow-free thermal roughness
                                      ! on tiles
REAL, POINTER :: dolr_field(:,:)      ! TOA - surface upward LW at
                            ! radiation timestep
REAL, POINTER :: lw_down(:,:)         ! Surface downward LW at
                            ! radiation timestep
! (changed from SW_TILE b/c of naming conflicts)
REAL, POINTER :: sw_tile_rts(:,:)     ! Surface net SW on land tiles at
                            ! radiation timestep

! MORUSES - new urban two-tile scheme
REAL, POINTER :: hgt(:)      ! Building height (m)
REAL, POINTER :: hwr(:)      ! Height to width
REAL, POINTER :: wrr(:)      ! Width ratio
REAL, POINTER :: disp(:)     ! Displacement height (m)
REAL, POINTER :: ztm(:)      ! Eff roughness length of momentum (m)
REAL, POINTER :: albwl(:)    ! Wall albedo
REAL, POINTER :: albrd(:)    ! Road albedo
REAL, POINTER :: emisw(:)    ! Wall emmissivity
REAL, POINTER :: emisr(:)    ! Road emmissivity

!   2.14: Carbon cycle fields
REAL, POINTER :: triffid_co2_d1(:) 
                 ! TRIFFID-derived CO2 fluxes for passing to the 
                 ! atmosphere in emissions-driven runs (kgC/m2/yr):
                 ! exudates + wood product pool flux + harvest flux
REAL, POINTER :: co2flux(:,:)   ! Ocean CO2 flux (Kg CO2/m2/s1)
REAL, POINTER :: co2_emits(:,:) ! Surface CO2 emissions (Kg CO2/m2/s)

REAL, POINTER :: soil_thickness(:)

! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
! All 1d (in the vertical)
REAL, POINTER :: zseak_theta(:) ! zsea(k) on theta levels
REAL, POINTER :: Ck_theta(:)    ! C(k)    on theta levels
REAL, POINTER :: zseak_rho(:)   ! zsea(k) on rho levels
REAL, POINTER :: Ck_rho(:)      ! C(k)    on rho levels

!   2.16: Fields for large-scale hydrology scheme. All 1d (land_field).

REAL, POINTER :: ti_mean(:)          !Mean topographic index
REAL, POINTER :: ti_sig(:)           !Standard dev. in topographic index
REAL, POINTER :: fexp(:)             !Exponential decay in soil
!                                  ! saturated conductivity
REAL, POINTER :: gamma_int(:)        !Integrated gamma function
REAL, POINTER :: water_table(:)      !Water table depth
REAL, POINTER :: fsfc_sat(:)         !Surface saturation fraction
REAL, POINTER :: f_wetland(:)        !Wetland fraction

REAL, POINTER :: sthzw(:)
REAL, POINTER :: a_fsat(:)
REAL, POINTER :: c_fsat(:)
REAL, POINTER :: a_fwet(:)
REAL, POINTER :: c_fwet(:)

!   2.17: Fields for River routing.
REAL, POINTER :: riv_sequence(:)   ! River sequence
REAL, POINTER :: riv_direction(:)  ! River direction
REAL, POINTER :: riv_storage(:)    ! River water storage
REAL, POINTER :: tot_surfroff(:)   ! Accumulated surface runoff
REAL, POINTER :: tot_subroff(:)    !     "       sub-surface runoff
REAL, POINTER :: riv_inlandatm(:)       ! inland basin outflow

! Fields for water conservation correction due to lake evaporation:
REAL, POINTER :: acc_lake_evap(:)  ! Acc lake evaporation

! Fields for grid-to-grid river routing (river routing 2A)
REAL, POINTER :: riv_iarea(:)      ! Drainage area
REAL, POINTER :: riv_slope(:)      ! Grid-cell slope
REAL, POINTER :: riv_flowobs1(:)   ! Initial values of flow
REAL, POINTER :: riv_inext(:)      ! Flow direction (x)
REAL, POINTER :: riv_jnext(:)      ! Flow direction (y)
REAL, POINTER :: riv_land(:)       ! Land-type (land/river/sea)
REAL, POINTER :: riv_substore(:)   ! Subsurface storage
REAL, POINTER :: riv_surfstore(:)  ! Surface storage
REAL, POINTER :: riv_flowin(:)     ! Surface inflow
REAL, POINTER :: riv_bflowin(:)    ! Subsurface inflow

! Fields used when coupling using OASIS.
REAL, POINTER :: c_solar(:,:)       ! CPL solar radn
REAL, POINTER :: c_blue(:,:)        ! CPL blue radn
REAL, POINTER :: c_longwave(:,:)    ! CPL lw radn
REAL, POINTER :: c_taux(:,:)        ! CPL taux
REAL, POINTER :: c_tauy(:,:)        ! CPL tauy
REAL, POINTER :: c_w10(:,:)         ! CPL 10m wind
REAL, POINTER :: c_sensible(:,:)    ! CPL sensible ht flx
REAL, POINTER :: c_sublim(:,:,:)    ! CPL Sublimation rate
REAL, POINTER :: c_evap(:,:)        ! CPL Evaporation
REAL, POINTER :: c_fcondtopn(:,:,:) ! CPL Fcondtop
REAL, POINTER :: c_topmeltn(:,:,:)  ! CPL Topmelt
REAL, POINTER :: c_tstar_sicen(:,:,:)  ! CPL Sea ice surface temp
REAL, POINTER :: c_lsrain(:,:)
REAL, POINTER :: c_lssnow(:,:)
REAL, POINTER :: c_cvrain(:,:)
REAL, POINTER :: c_cvsnow(:,:)
REAL, POINTER :: c_riverout(:,:)
REAL, POINTER :: c_calving(:,:)
REAL, POINTER :: c_mslp(:,:)
REAL, POINTER :: c_surf_CO2(:,:)    ! CPL pp of CO2 
REAL, POINTER :: c_dust_dep(:,:)    ! CPL the total dust deposition

!   2.18: JULES variables
REAL, POINTER :: snowdepth(:,:)      ! Snow depth on ground on tiles (m)
REAL, POINTER :: rho_snow_grnd(:,:)  ! Snowpack bulk density (kg/m3)
REAL, POINTER :: nsnow(:,:)          ! Number of snow layers on ground on tiles
REAL, POINTER :: ds(:,:,:)           ! Snow layer thickness (m)
REAL, POINTER :: sice(:,:,:)         ! Snow layer ice mass on tiles (Kg/m2)
REAL, POINTER :: sliq(:,:,:)         ! Snow layer liquid mass on tiles (Kg/m2)
REAL, POINTER :: tsnowlayer(:,:,:)   ! Snow layer temperature (K)
REAL, POINTER :: rho_snow(:,:,:)     ! Snow layer densities (kg/m3)
REAL, POINTER :: rgrainl(:,:,:)      ! Snow layer grain size on tiles (microns)
!  FLake lake scheme (all 1d, land_field)
REAL, POINTER :: lake_depth(:)
REAL, POINTER :: lake_fetch(:)
REAL, POINTER :: lake_t_mean(:)
REAL, POINTER :: lake_t_mxl(:)
REAL, POINTER :: lake_t_ice(:)
REAL, POINTER :: lake_h_mxl(:)
REAL, POINTER :: lake_h_ice(:)
REAL, POINTER :: lake_shape(:)
REAL, POINTER :: lake_g_dt(:)

! Aerosol climatologies
REAL, POINTER :: arclbiog_bg(:,:,:) ! Biogenic aerosol climatology
REAL, POINTER :: arclbiom_fr(:,:,:) ! Biomass burning (fresh) aerosol clim
REAL, POINTER :: arclbiom_ag(:,:,:) ! Biomass burning (aged) aerosol clim
REAL, POINTER :: arclbiom_ic(:,:,:) ! Biomass burning (in-cloud) aerosol clim
REAL, POINTER :: arclblck_fr(:,:,:) ! Black carbon (fresh) aerosol clim
REAL, POINTER :: arclblck_ag(:,:,:) ! Black carbon (aged) aerosol clim
REAL, POINTER :: arclsslt_fi(:,:,:) ! Sea salt (film mode) aerosol clim
REAL, POINTER :: arclsslt_jt(:,:,:) ! Sea salt (jet mode) aerosol clim
REAL, POINTER :: arclsulp_ac(:,:,:) ! Sulphate (accumulation mode) aero clim
REAL, POINTER :: arclsulp_ak(:,:,:) ! Sulphate (Aitken mode) aerosol clim
REAL, POINTER :: arclsulp_di(:,:,:) ! Sulphate (dissolved) aerosol clim
REAL, POINTER :: arcldust_b1(:,:,:) ! Dust (bin 1) aerosol climatology
REAL, POINTER :: arcldust_b2(:,:,:) ! Dust (bin 2) aerosol climatology
REAL, POINTER :: arcldust_b3(:,:,:) ! Dust (bin 3) aerosol climatology
REAL, POINTER :: arcldust_b4(:,:,:) ! Dust (bin 4) aerosol climatology
REAL, POINTER :: arcldust_b5(:,:,:) ! Dust (bin 5) aerosol climatology
REAL, POINTER :: arcldust_b6(:,:,:) ! Dust (bin 6) aerosol climatology
REAL, POINTER :: arclocff_fr(:,:,:) ! Org carbon fossil fuel (fresh) aero clim
REAL, POINTER :: arclocff_ag(:,:,:) ! Org carbon fossil fuel (aged) aero clim
REAL, POINTER :: arclocff_ic(:,:,:) ! Org carb fossil fuel (in-cloud) aero clim
REAL, POINTER :: arcldlta_dl(:,:,:) ! Delta aerosol climatology

! Convective Cloud Fields
REAL, POINTER :: lcbase(:,:)
REAL, POINTER :: ccw_rad(:,:,:)
REAL, POINTER :: cca_dp(:,:,:)
REAL, POINTER :: cca_md(:,:,:)
REAL, POINTER :: cca_sh(:,:,:)

END MODULE atm_fields_real_mod


