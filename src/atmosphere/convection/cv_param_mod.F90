! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with convection.

MODULE cv_param_mod

  ! Description:
  !   Module containing parameters used by the convection code.
  !
  ! Method:
  !   Default values have been declared where appropriate.
  !
  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v7.4 programming standards.
  !
  ! Declarations:

USE conversions_mod, ONLY: rsec_per_day  ! Get seconds per day, to easily 
                                         ! express constants in units per day
                                         ! but store in units per second.

IMPLICIT NONE
SAVE

!===========================================================================
! The convection scheme could work in specific humidities or mixing ratios.
! The original UM code was designed to work in specific humidities, input
! in mixing ratios gets converted to specific humidities before calling
! shallow, deep or mid-level convection.
! The turbulence versions of the convection code have been written to work
! in mixing ratios.
! The new downdraught and evaporation code will work for either 
! specific humidity or mixing ratio as input, dependent on the setting of the 
! logical l_cv_mix_ratio.
! At present this has to be set to .false. as the updraught mass flux part of 
! convection scheme which calls the downdraught code is working in specific 
! humidity. 
!===========================================================================
   
LOGICAL :: l_cv_mix_ratio = .FALSE.  ! true if scheme to work in mixing ratios.

!===========================================================================
! Integer parameters to remove use of magic numbers
!===========================================================================

! Basis for 2d convective cloud amount
!   Total Condensed Water(TCW) original 4A
!   LES mb/wsc scalings (Grant and Lock, 2004) (shallow)
!   Surface rain from respective convective cloud type (Dp and Md-level)
INTEGER, PARAMETER :: cca2d_total_condensed_water = 0
INTEGER, PARAMETER :: cca2d_grant_lock            = 1 ! Shallow cnv (CCRad)
INTEGER, PARAMETER :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv (CCRad)

INTEGER, PARAMETER :: total_condensed_water = 0
INTEGER, PARAMETER :: grant_lock_over       = 1 ! Sh cu with overlap (CCRad)
INTEGER, PARAMETER :: grant_lock            = 2 ! Sh cu no overlap (CCRad)
INTEGER, PARAMETER :: srf_precip            = 1 ! mid/deep cnv (CCRad)

! Convective cloud decay
INTEGER, PARAMETER :: rad_decay_off           = 0
INTEGER, PARAMETER :: rad_decay_full_timestep = 1
INTEGER, PARAMETER :: rad_decay_conv_substep  = 2 ! CCRad only


! Convective cloud decay timescale
INTEGER, PARAMETER :: cld_life_constant = 0
INTEGER, PARAMETER :: cld_life_func_hgt = 1 ! CCRad only

! For all the options below except 3), the anvil base is at the freezing
! Level
INTEGER, PARAMETER :: anv_pressure      = 0
INTEGER, PARAMETER :: anv_height        = 1
INTEGER, PARAMETER :: anv_model_levels  = 2
INTEGER, PARAMETER :: anv_limited_pressure_depth = 3

! Energy correction for 6a convection
INTEGER, PARAMETER :: method_en_rho=1     ! Correct energy but not water on
                                          ! the rho grid
INTEGER, PARAMETER :: method_en_mx_rho=2  ! Correct water and then energy on
                                          ! the rho grid
INTEGER, PARAMETER :: method_en_qx_p=3    ! Correct water and then energy on
                                          ! the pressure level grid. Not
                                          ! recommended for normal use a
                                          ! useful check when developing
                                          ! scheme

! Shallow scheme cloud base closure options
INTEGER, PARAMETER ::  sh_wstar_closure = 0 ! standard closure on wstar
INTEGER, PARAMETER ::  sh_grey_closure  = 1 ! grey zone closure

! Parameters used to relate cca_2d of cld to its precipitation rate
REAL, PARAMETER    :: a_land = 0.3
REAL, PARAMETER    :: a_sea  = 0.3
REAL, PARAMETER    :: b_land = 0.025
REAL, PARAMETER    :: b_sea  = 0.025

! Application of convective cloud anvils
REAL, PARAMETER :: deep_dp = 50000.0  ! Min. depth for anvil criteria(Pa)

! Critical depth of cloud for the formation of
! convective precipitation over sea (m)
REAL, PARAMETER :: critdsea = 1.5e3

! critical depth of cloud for the formation of convective
! precipitation over land (m)
REAL, PARAMETER :: critdlnd = 4.0e3

! critical depth of a glaciated cloud for the formation of
! convective precipitation (m)
REAL, PARAMETER :: critdice = 1.0e3

! Parcel ascent in diagnosis, cloud water condensate
REAL, PARAMETER :: qlcrit = 1.0e-3   ! critical cloud water

! Timestep frequency for calling convection - hardwired to call every timestep
INTEGER,PARAMETER :: a_conv_step = 1

! Threshold magnitude of convective temperature tendency to determine whether 
! convection is "active" on a given model-level (used for the conv_prog_precip
! experimental prognostic field).
REAL,PARAMETER :: dthetadt_conv_active_threshold = 0.001 / rsec_per_day
                               ! Set to 0.001 K day-1

! Minimum bound applied to the convective precipitation rate when updating the
! conv_prog_precip experimental prognostic field / kg m-2 s-1
REAL,PARAMETER :: conv_prog_precip_min_threshold = 1.0E-5

!============================================================================
! Parcel ascent
!============================================================================

  ! coefficients used in calculation of entrainment rate
REAL, PARAMETER :: ae1     = 1.0        ! Not used
REAL, PARAMETER :: ae2     = 1.5

! Reference depth for variable entrainment and mixing detrainment for deep
! convection. It is the depth at which the variable entrainment
! and m.det will be the same as the original rates.
REAL, PARAMETER :: refdepth_dp = 8000.0

! Reference saturated humidity in kg/kg for variable entrainment.
REAL, PARAMETER :: refqsat = 0.02 

! minimum excess buoyancy to continue parcel ascent (K)

REAL,PARAMETER:: xsbmin = 0.2       ! Used in mid scheme

! initial excess potential temperature (k) and mixing ratio
! (kg/kg) for deep convection
REAL, PARAMETER :: thpixs_deep= 0.2
REAL, PARAMETER :: qpixs_deep =0.0

! initial excess potential temperature (k) and mixing ratio
! (kg/kg) for shallow convection
REAL, PARAMETER :: thpixs_shallow = 0.2
REAL, PARAMETER :: qpixs_shallow  = 0.0

! initial excess potential temperature (k) and mixing ratio
! (kg/kg) for mid-level convection
REAL, PARAMETER :: thpixs_mid= 0.2
REAL, PARAMETER :: qpixs_mid =0.0

! Minimum parcel buoyancy/layer thickness (K/Pa)
REAL, PARAMETER :: mparb = 1.0            ! Used in mid scheme

! Constants to determine initial convective mass flux from parcel buoyancy
! Deep convection
REAL, PARAMETER :: c_deep = 5.17e-4       ! No longer used
REAL, PARAMETER :: d_deep = 0.0           ! No longer used

! Shallow convection
REAL, PARAMETER :: c_shallow = 5.17e-4    ! No longer used
REAL, PARAMETER :: d_shallow = 0.0        ! No longer used

! Mid convection
REAL, PARAMETER :: c_mid = 5.17e-4
REAL, PARAMETER :: d_mid = 0.0

! limits on the initial convective parcel perturbations
REAL, PARAMETER :: max_diag_thpert  =  2.0
REAL, PARAMETER :: max_dp_thpert    =  2.0
REAL, PARAMETER :: min_dp_thpert    = -2.0
REAL, PARAMETER :: max_sh_thpert    =  2.0
REAL, PARAMETER :: min_sh_thpert    = -2.0
REAL, PARAMETER :: max_dp_qpert_fac =  0.2
REAL, PARAMETER :: max_sh_qpert_fac =  0.2

! mparfl = 1E-3 * minimum parcel buoyancy * mass flux parameter c
REAL, PARAMETER :: mparfl = 1.0e-3 * 1.0 * 3.33e-4

! Difference in potential temperature between levels above which the
! atmosphere is assumed to be too stable to convect (K)
REAL, PARAMETER :: delthst = 0.5

! Entrainment coefficient 

REAL, PARAMETER :: entcoef =3.0

!============================================================================
! Convective closure
!============================================================================

  ! Coefficient relating sub-cloud convective velocity scale to cumulus
  ! mass flux for shallow convection
REAL, PARAMETER :: c_mass=0.03

! Tuneable factor used in denominator of W_CAPE timescale calculation
REAL, PARAMETER :: wcape_fac = 3.0

! Parameter governing the speed of transition from parametrized to
! resolved convection, for use with cldbase_opt_sh = sh_grey_closure
REAL :: beta_cu = 0.15

! CAPE closure based on large-scale w
REAL, PARAMETER :: max_cape = 3600.0*4.0  ! 4 hours at present
REAL, PARAMETER :: a_cape = 3600.0*0.08  ! value from CASCADE fit
REAL, PARAMETER :: b_cape = -0.7        ! value from CASCADE fit

! Deep cloud closure based on wstar and large-scale w
! wup_cb   = a_cb * wstar                           (units m/s)
! sigma_cb = b_cb + c_cb * wls                      (fractional area)
! mf_cb    = rho_cb * sigma_cb * wup_cb             (Kg/m2/s)
! massflux_cb = g * rho_cb * sigma_cb * wup_cb         (Pa/s)
REAL, PARAMETER :: a_cb = 2.6          ! from fit to CRM simulations
REAL, PARAMETER :: b_cb = 0.02         !
REAL, PARAMETER :: c_cb = 0.2          !

!============================================================================
! Downdraught and evaporation below cloud base calculations
!============================================================================

  ! Coefficients used in calculation of downdraught entrainment
  ! rates
REAL, PARAMETER :: ddcoef1 = 1.8e6
REAL, PARAMETER :: ddcoef2 = 3.0

! Thickness level used in calculation of mixing detrainment for
! downdraught  (pa)
REAL, PARAMETER :: det_lyr = 10000.0

! exponents used in calculation of evaporation of liquid
REAL, PARAMETER :: p_lq1 = 0.52
REAL, PARAMETER :: p_lq2 = 0.67

! exponents used in calculation of evaporation of ice
REAL, PARAMETER :: p_ice1 = 0.55
REAL, PARAMETER :: p_ice2 = 0.76

! exponents and constants associated with density term in
! evaporation of liquid
REAL, PARAMETER :: rho_lqp1 = 0.26
REAL, PARAMETER :: rho_lqp2 = 0.59
REAL, PARAMETER :: rho_lqa  = 108.80
REAL, PARAMETER :: rho_lqb  = 830.73

! exponents and constants associated with density term in
! evaporation of ice
REAL, PARAMETER :: rho_icp1 = 0.28
REAL, PARAMETER :: rho_icp2 = 0.63
REAL, PARAMETER :: rho_icea = 1569.52
REAL, PARAMETER :: rho_iceb = 32069.02

! constants used in quadratic formula for evaporation of liquid
REAL, PARAMETER :: lq_a = 2.008e-9
REAL, PARAMETER :: lq_b = -1.385e-6
REAL, PARAMETER :: lq_c = 2.424e-4

! constants used in quadratic formula for evaporation of ice
REAL, PARAMETER :: ice_a = -5.2e-9
REAL, PARAMETER :: ice_b = 2.5332e-6
REAL, PARAMETER :: ice_c = -2.911e-4
REAL, PARAMETER :: ice_d = 1.7405e-5  ! value of A(T,p) at T=243.58
                                      ! See Gregory 1995 in UM DOC 27

! Scaling factor used to calculate area occupied by precip below cloud-base
! when no downdraught, used in evaporation and melting calculations
REAL, PARAMETER :: precip_area_fac = 1.0 ! precip area = CCA * precip_area_fac

! Scaling factor used to calculate area occupied by downdraught,
! used in evaporation calc
REAL, PARAMETER :: dd_area_fac = 0.5     ! dd area = CCA * dd_area_fac

! downdraught precipitation transfer efficiency factor
REAL, PARAMETER :: ddptef = 2.0

!============================================================================
! Convective momentum transport (CMT) calculations
!============================================================================

  ! Coefficient for "pressure term" relative to "shear term" as found in
  ! Kershaw & Gregory 1997

REAL, PARAMETER :: cpress_term = 0.7

! Deep turbulent CMT scheme
REAL, PARAMETER :: dp_cmt_gamma = 1.63
REAL, PARAMETER :: dp_cmt_delta = 2.0
REAL, PARAMETER :: top_press    = 15000.0




! Parameters specifying methods in calculation of w_eqn
INTEGER, PARAMETER :: SimpsonWiggert = 1
INTEGER, PARAMETER :: BuoyScale = 2

!------------------------------------------------------------------------------
! Current parameter settings for w-eqn options
! Currently hardwired for development
!------------------------------------------------------------------------------

REAL,    PARAMETER :: w2pi           = 1.0    ! Initial w^2 at cloud base
REAL,    PARAMETER :: gamma_in_w_eqn = 0.5    ! Virtual mass coefficient
REAL,    PARAMETER :: cumulus_r      = 100.0  ! Radius of cumulus tower ??
REAL,    PARAMETER :: drag_coeff     = 1.25   ! Aerodynamic drag coefficient
                                              ! for solid spheres
REAL,    PARAMETER :: k2_const       = 0.71   ! Entrainment coefficient
REAL,    PARAMETER :: fix_alpha      = 0.0001 ! Value for alpha if constant
REAL,    PARAMETER :: gamma_b        = 130.0  ! Scaling value for buoyancy
REAL,    PARAMETER :: NegBuoyMinW    = 0.1    ! Min. threshold of buoyancy

INTEGER, PARAMETER :: watload_opt    = 2      ! How to apply water loading
                                              ! for SimpsonWiggert W calc
                                              ! 1) SimpsonWiggert
                                              ! 2) As in UM
INTEGER, PARAMETER :: NegBuoyOpt     = 0      ! Option to decide what to do
                                              ! when encountering -ve buoyancy
INTEGER, PARAMETER :: alpha_opt      = 2      ! How to set alpha
                                              ! 1) As with SimpsonWiggert
                                              ! 2) Us UM entrainment values
                                              ! 3) Set as constant value
INTEGER, PARAMETER :: wCalcMethod    = 1      ! Sets method to calculate w
                                              ! 1) SimpsonWiggert
                                              ! 2) Buoyancy Scaling


!------------------------------------------------------------------------------
! END w-eqn options
!------------------------------------------------------------------------------

!============================================================================
! Diagnostics passed to the perturbation forecast (PF) model in VAR
! This parameter will not affect model evolution on its own but will 
! alter the evolution through the analysis cycle if l_fix_conv_diags_var=True
! because it affects the mass flux diagnostic that is passed to VAR.
! The chosen value works well but may not be optimal. It is possible that the 
! optimal value of this parameter will change slightly with vertical
! resolution.
!============================================================================

REAL, PARAMETER :: max_mf_fall = 0.3          ! The maximum fractional fall
                                              ! in mass flux between k and k+1
                                              ! before the mass flux is masked
                                              ! out

!============================================================================
! Convective cold pools
!============================================================================

! Scaling constant for convective cold-pool vertical kinetic energy scale
! in conv_diag
REAL, PARAMETER :: phi_ccp = 0.01

END MODULE cv_param_mod
