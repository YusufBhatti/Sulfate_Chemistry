! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with convection.

MODULE cv_run_mod

  ! Description:
  !   Module containing runtime logicals/options used by the convection code.
  !
  ! Method:
  !   All switches/options which are contained in the &Run_Convection
  !   namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate.
  !
  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.
  !
  ! Declarations:

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
SAVE

! To satisfy the requirements of ROSE, the defaults for all these
! switches/variables have been set to FALSE/MDI.
! The original defaults are noted in comments below.
! These values are also indicated in the ROSE help text.

INTEGER :: i_convection_vn = imdi
                      ! Switch to determine version of convection scheme
                      ! 5 => 5A scheme
                      ! 6 => 6A scheme
                      ! 10 => Lambert-Lewis scheme
INTEGER, PARAMETER :: i_convection_vn_5a = 5
INTEGER, PARAMETER :: i_convection_vn_6a = 6
INTEGER, PARAMETER :: i_cv_llcs = 10
INTEGER, PARAMETER :: ccp_off = 0

!===========================================================================
! Logical switches set from GUI
!===========================================================================

! Comments for TRUE status

LOGICAL :: l_fix_udfactor = .FALSE. ! Fix application of UD_FACTOR (lock)

LOGICAL :: l_cloud_deep   = .FALSE. ! Use depth criterion for convective anvils
                                    ! Original default TRUE
LOGICAL :: l_mom          = .FALSE. ! Use Convective Momentum Transport

LOGICAL :: l_mom_dd       = .FALSE. ! Use Convective Momentum Transport in
                                    ! new Downdraught scheme

LOGICAL :: l_eman_dd      = .FALSE. ! Use Emanuel downdraught scheme

LOGICAL :: l_new_dd       = .FALSE. ! Use new downdraught and evaporation 
                                    ! scheme 

LOGICAL :: l_use_dd       = .FALSE. ! Use the downdraught part of the new
                                    ! scheme. The alternative is to call the new
                                    ! downdraught and evaporation code just 
                                    ! allowing precipitation to evaporate but
                                    ! not allowing downdraughts to form.

LOGICAL :: l_rediagnosis  = .FALSE. ! If more than one convection step per model
                                    ! physics timestep then rediagnose cumulus
                                    ! points by recalling conv_diag before each
                                    ! convective sub-step.
                                    ! Original default TRUE

LOGICAL :: l_anvil        = .FALSE. ! Apply anvil scheme to 3D convective cloud.
                                    ! No effect unless l_3d_cca = .TRUE.

LOGICAL :: l_snow_rain    = .FALSE. ! If .TRUE. allows a mix of snow and rain
                                    ! between the freezing temperature
                                    ! (273.15K) and t_melt_snow

LOGICAL :: l_dcpl_cld4pc2 = .FALSE. ! Decouples section 0/5 cloud properties
                                    ! from each other, use with PC2. 
                                    ! May change results if used without PC2.

LOGICAL :: l_murk_conv    = .FALSE. ! Enables convective mixing of MURK as
                                    ! a tracer

LOGICAL :: l_safe_conv    = .FALSE. ! Safer convection, remove negative q
                                    ! before attempting convection, don't
                                    ! add increments for ascents with
                                    ! negative CAPE. (5A/6A only)

LOGICAL :: l_ccrad        = .FALSE. ! Main Logical, will include code
                                    ! connected with CCRad.
                                    ! (including bugfixes)

LOGICAL :: l_3d_cca       = .FALSE. ! Use 3D convective cloud amount

LOGICAL :: l_conv_hist    = .FALSE. ! True if 3 extra prognostics holding
                                    ! convective history information.

LOGICAL :: l_param_conv   = .FALSE. ! Run time switch for convection scheme

LOGICAL :: l_conv_prog_group_1  = .FALSE.
                                    ! Use convection experimental prognostics
                                    ! group 1 (currently 1 field, conv_prog_1)
LOGICAL :: l_conv_prog_group_2  = .FALSE.
                                    ! Use convection experimental prognostics
                                    ! group 2 (currently 1 field, conv_prog_2)
LOGICAL :: l_conv_prog_group_3  = .FALSE.
                                    ! Use convection experimental prognostics
                                    ! group 3 (currently 1 field, conv_prog_3)
LOGICAL :: l_conv_prog_precip  = .FALSE.
                                    ! Use convection 3D prognostic field of
                                    ! recent convective precipitation.
LOGICAL :: l_jules_flux = .FALSE.
                                    ! Use the real surface fluxes calculated
                                    ! by Jules in the diagnosis rather
                                    ! than the simple diagnosis estimate

!===========================================================================
! Logical switches not set from GUI
!===========================================================================

LOGICAL :: l_pc2_diag_sh  = .FALSE. ! If true uses  diagnostic convective
                                    ! shallow cloud in PC2 replacing
                                    ! prognotsic qcl etc

LOGICAL :: l_cv_conserve_check = .FALSE.
                                    ! Diagnostic conservation checks. Goes
                                    ! through energy correction code
                                    ! (5A/6A only)

LOGICAL :: l_cmt_heating  = .FALSE. ! Include the heating due to the
                                    ! dissipation of kinetic energy
                                    ! from the convective momentum transport

LOGICAL :: l_mr_conv = .FALSE.      ! Flag for convection running with
                                    ! mixing ratios rather than specifics.

LOGICAL :: l_mcr_conv = .FALSE.     ! Flag for convection adding increments
                                    ! to the optional microphysics prognostic
                                    ! hydrometeor fields qrain, qcf2, qgraup.

LOGICAL :: l_wvar_for_conv = .FALSE.! Flag to calculate boundary-layer
                                    ! turbulent vertical velocity variance
                                    ! for use by the convection scheme.

!===========================================================================
! Integer options set from GUI
!===========================================================================

! Convection integer options set from gui/namelist Convection Scheme

INTEGER :: n_conv_calls       = imdi ! Number of calls to convection
                                     ! per physics timestep
                                     ! Original default 1

INTEGER :: cld_life_opt       = imdi ! Convective cloud decay time
                                     ! Original default cld_life_constant (0)
INTEGER :: rad_cloud_decay_opt= imdi ! Convective cloud decay
                                     ! Original default rad_decay_off (0)
INTEGER :: anv_opt            = imdi ! Anvil cloud basis
                                     ! Original default anv_model_levels (2)

! Convection Scheme Options (5A)
! NOTE: These options were valid at the time of writing. They are used in
!       development code(5A) and so very likely to change.
!       Users should consult the Convection Group for current available
!       options.

INTEGER :: iconv_shallow  = imdi  ! Shallow (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: Turbulence scheme (non-precipitating)
                           !   3: Turbulence scheme (precipitating)

INTEGER :: iconv_congestus = imdi ! Congestus (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme untested for possible
                           !      future use. 

INTEGER :: iconv_mid       = imdi ! Mid-level (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: Future use

INTEGER :: iconv_deep      = imdi ! Deep (Original default 0)
                           !   0: no scheme,
                           !   1: Gregory-Rowntree scheme
                           !   2: turbulence scheme

INTEGER :: deep_cmt_opt    = imdi ! Deep CMT tunings (Original default 0)
                           !   0: turbulence based
                           !   1: Operational 70 level (modified turb based)
                           !   2: Gregory-Kershaw scheme
                           !   3: New turbulence based scheme.
                           !   4: Future use
                           !   6: Stabilized Gregory-Kershaw scheme

INTEGER :: mid_cmt_opt     = imdi ! Mid CMT scheme to be used (Orig default 0)
                           !   0: Gregory-Kershaw scheme
                           !   1: Diffusive scheme
                           !   2: Stabilized Gregory-Kershaw scheme

INTEGER :: icvdiag    = imdi ! Diagnosis calculation options (Orig default 1)
                           !   0: No longer allowed
                           !   1: 5A/6A scheme default
                           !   2: 0.55/z entrainment rates
                           !   3: 1/z entrainment
                           !   4: Diurnal cycle over land entrainment rates
                           !   5: Diurnal cycle over land entrainment rates,
                           !      Ocean undilute ascent.
                           ! 6A scheme only 
                           !   6: dilute, p/(p*)^2 entrainment rate dependent
                           !      on precipitation based 3d convective 
                           !      prognostic at the initiation level, as 
                           !      ent_opt_dp=6
                           !   7: dilute, p/(p*)^2 entrainment rate dependent
                           !      on precipitation based 3d convective 
                           !      prognostic at the current level, as 
                           !      ent_opt_dp=7
                           !   8: fixed p/(p*)^2 entrainment rate, as 
                           !      ent_opt_dp=0

INTEGER :: cvdiag_inv = imdi ! Inversion test in convective diagnosis? 5A & 6A
                           ! Original default 1
                           ! When doing an undilute parcel ascent
                           !   0: No inversion test
                           !   1: Original inversion test (Default 5A and 6A)
                           !   2: Future alternative inversion tests

INTEGER :: tv1_sd_opt = imdi ! Standard dev of level virtual temperature options
                           ! Original default 0
                           !   0: Assume BL depth of 300m (Default)
                           !   1: Use calculated BL depth
                           !   2: (1) plus stability dependence and coastal mean

INTEGER :: adapt    = imdi ! Adaptive detrainment/entrainment options
                           ! Original default 0
                           !   0: Original (Default)
                           !   1: Adaptive detrainment: mid and deep convection
                           !   2: Future/Experimental (En/Detrainment)
                           !   3: Adaptive detrainment: deep convection
                           !   4: Adaptive detrainment: shallow, mid and
                           !      deep convection
                           !   5: Smoothed adaptive detrainment:
                           !      mid and deep convection
                           !   6: Smoothed adaptive detrainment:
                           !      shallow, mid and deep convection
                           !   7: Improved smoothed adaptive detrainment:
                           !      mid and deep convection
                           !   8: Improved smoothed adaptive detrainment:
                           !      shallow, mid and deep convection

INTEGER :: fdet_opt = imdi ! Forced detrainment calculation options. 6A only.
                           !   0: calculate forced detrainment rate using
                           !      humidity continuity equation.
                           !   1: calculate forced detrainment rate using
                           !      theta continuity equation.
                           !   2: forced detrainment rate based on theta
                           !      equation with an improved treatment of
                           !      subsaturated conditions.
                           !   3: forced detrainment rate based on theta
                           !      equation with improved treatment of
                           !      subsaturated conditions and improved
                           !      numerical calculation of detrained parcel
                           !      properties.

INTEGER :: ent_opt_dp = imdi  ! Deep entrainment option (Orig default 0)
                           !   0: original Ap/(p*)^2
                           !   1: n/z dependence where n=ent_fac_dp
                           !   2: As Ap/(p*)^2 but multiplied by extra factor
                           !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))
                           !   3: factor * (A/p*)*((p/p*)^ent_dp_power)
                           !   Options 4 & 5 are experimental options that
                           !   depend on the diagnosed depth of convection.
                           !   Available for 6a only.
                           !   4: variable n/z style entrainment
                           !   5: variable p/(p*)^2  style entrainment
                           !   6: variable p/(p*)^2 entrainment dependent on 
                           !     precipitation based 3d convective prognostic 
                           !     at the initiation level.
                           !   7: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !      at the current level.
INTEGER :: ent_opt_md = imdi  ! mid entrainment option (Orig default 0)
                           !   0: original Ap/(p*)^2
                           !   1: n/z dependence where n=ent_fac_md
                           !   2: As Ap/(p*)^2 but multiplied by extra factor
                           !      f=1.+3.*(1.-(100000.-p(k))/(100000.-50000.))
                           !   3: factor * (A/p*)*((p/p*)^ent_md_power)
                           !   6: variable p/(p*)^2 entrainment dependent on 
                           !     precipitation based 3d convective prognostic 
                           !     at the initiation level.
                           !   7: variable p/(p*)^2 entrainment dependent on
                           !     precipitation based 3d convective prognostic
                           !      at the current level.
INTEGER :: mdet_opt_dp = imdi  ! Deep mixing detrainment option - 6a only
                           !   0: original entrainment/3
                           !   1: amdet_fac*entrainment*(1-rh)
                           !   2: dependent on diagnosed depth of convection.
                           !      Deeper/shallower convection has
                           !      lower/higher mixing detrainment.

! Deep cloud base closure options 5A & 6A
INTEGER :: cldbase_opt_dp = imdi
                           ! This option may have one of the following values:
                           !   0: RH based timescale
                           !   1: RH based timescale (timestep limited)
                           !   2: Fixed timescale
                           !   3: Fixed timescale reduced if w exceeds w max
                           !      CAPE closure
                           !   4: Area scaled CAPE closure
                           !   5: Experimental - not allowed
                           !   6: RH based timescale (timestep limited)
                           !      plus reduction if w exceeds w-max
                           !   7: CAPE timescale dependent on large-scale
                           !      vertical velocity
                           !   8: Closure based on boundary layer fluxes
                           !      and large-scale vertical velocity based CAPE
                           !      timescale
                           !   9: Closure based on boundary layer fluxes
                           !      and large-scale vertical velocity.

INTEGER :: cldbase_opt_md = imdi  ! Mid-level cloud base closure options 5A&6A
                           !   0: RH based timescale
                           !   1: RH based timescale (timestep limited)
                           !   2: Fixed timescale
                           !   3: Fixed timescale reduced if w exceeds w max
                           !      CAPE closure
                           !   4: Area scaled CAPE closure
                           !   5: Experimental - not allowed
                           !   6: RH based timescale (timestep limited)
                           !      plus reduction if w exceeds w-max
                           !   7: CAPE timescale dependent on large-scale
                           !      vertical velocity

INTEGER :: cldbase_opt_sh = imdi  ! Shallow cloud base closure options 5A&6A
                           !   0: Closure based on boundary layer fluxes
                           !   1: Grey zone combination of fixed timescale
                           !      CAPE closure and BL fluxes

INTEGER :: cape_bottom = imdi    ! Start level for w_max in column
                                 ! Original default IMDI

INTEGER :: cape_top    = imdi    ! End   level for w_max in column
                                 ! Original default IMDI

INTEGER :: sh_pert_opt  = imdi   ! Initial perturbation method for
                                 ! shallow cumulus (Orig default 0)
                                 ! 0 = Original code
                                 ! 1 = Revised  code

INTEGER :: limit_pert_opt = imdi ! Limits convective parcel perturbation
                                 ! to physically sensible values.
                                 ! Orig default 0
                                 ! 0 = original code - no limits
                                 ! 1 = apply limits to main ascent only
                                 ! 2 = apply limits to main ascent and in
                                 !     the convection diagnosis


INTEGER :: dd_opt         = imdi ! Downdraught scheme options
                                 ! Orig default 0
                                 ! 0 = Original code
                                 ! 1 = Revised  code

INTEGER :: termconv       = imdi ! Original default 0
                                 ! 0 for default
                                 ! 1 for modified termination condition

INTEGER :: bl_cnv_mix     = imdi ! Options for mixing convective increments in the BL
                                 ! Original default 0
                                 ! 0: original code
                                 ! 1: only mix the increments from the initial
                                 !    parcel perturbation

INTEGER :: cnv_wat_load_opt=imdi ! Options for including liquid and frozen water loading
                                 ! in the convective updraught buoyancy calculation
                                 ! Original default 0
                                 ! 0: Do not include water loading (default)
                                 ! 1: Include water loading

INTEGER :: cca2d_sh_opt = imdi   ! Method to evaluate cca2d (Shallow)
                                 ! Original default 0
INTEGER :: cca2d_md_opt = imdi   ! Method to evaluate cca2d (Mid-level)
                                 ! Original default 0
INTEGER :: cca2d_dp_opt = imdi   ! Method to evaluate cca2d (Deep)
                                 ! Original default 0

INTEGER :: ccw_for_precip_opt = imdi  ! Option controlling critical cloud water
                                      ! for the formation of precipitation
                                      ! Original default 0
                                      ! 0 - original code (dependent on
                                      !     fac_qsat & qlmin) Land sea
                                      !     dependence with dcrit.
                                      ! 1 - no Dcrit option
                                      ! 2 - Known as Manoj's first function
                                      ! 3 - Known as Manoj's congestus function
                                      ! 4 - no land sea diff dependence on
                                      !     fac_qsat & qlmin.

INTEGER :: plume_water_load = imdi ! Option for water loading in undilute
                                   ! parcel ascent
                                   ! Original default 0
                                   ! 0 - no water removal original undilute
                                   ! 1 - remove any water > 1g/kg
                                   ! 2 - remove any water > profile varying
                                   !     with qsat

INTEGER :: dil_plume_water_load  = imdi ! Option for water loading in dilute
                                        ! parcel ascent
                                        ! Original default 0
                                        ! 0 - no water removal
                                        ! 1 - remove any water > 1g/kg
                                        ! 2 - remove any water > profile varying
                                        !     with qsat

INTEGER :: cnv_cold_pools = imdi
                                ! Convective Cold Pool scheme options
                                ! 0 - scheme not used
                                ! 1 - CCP version 1

!===========================================================================
! Real values set from GUI
!===========================================================================

REAL :: cca_sh_knob = rmdi     ! Scales Shallow cloud fraction (CCRad)
                               ! Original default 1.0
REAL :: cca_md_knob = rmdi     ! Scales Mid     cloud fraction (CCRad)
                               ! Original default 1.0
REAL :: cca_dp_knob = rmdi     ! Scales Deep    cloud fraction (CCRad)
                               ! Original default 1.0
REAL :: ccw_sh_knob = rmdi     ! Scales Shallow cloud water (CCRad)
                               ! Original default 1.0
REAL :: ccw_md_knob = rmdi     ! Scales Mid     cloud water (CCRad)
                               ! Original default 1.0
REAL :: ccw_dp_knob = rmdi     ! Scales Deep    cloud water (CCRad)
                               ! Original default 1.0

REAL :: fixed_cld_life = rmdi ! Fixed convective cloud lifetime decay value
                              ! (seconds) (Original default 7200.0)

REAL :: cca_min = rmdi  ! Threshold value of convective cloud fraction
                        ! below which cca has neglible radiative impact and
                        ! is reset to zero (Original default 0.02)

REAL :: r_det   = rmdi  ! Parameter controlling adaptive detrainment -
                        ! Orig default 0.75 (operational)
                        ! HadGEM1a recommended 0.5

REAL :: cape_min     = rmdi ! Scale dependent min cape
                            ! Original default RMDI
REAL :: w_cape_limit = rmdi ! Test w for scale dependent cape timescale
                            ! Original default RMDI
REAL :: mparwtr      = rmdi ! Maximum value of the function that is used to calculate
                            ! the maximum convective cloud water/ice in a layer (kg/kg)
                            ! Original default 1.0e-3
REAL :: qlmin        = rmdi ! Minimum value of the function that is used to calculate
                            ! the maximum convective cloud water/ice in a layer (kg/kg)
                            ! Original default 2.0e-4
REAL :: fac_qsat     = rmdi ! Factor used to scale environmental qsat to give
                            ! the maximum convective cloud water/ice in a layer
                            ! Original default RMDI
REAL :: mid_cnv_pmin = rmdi ! The minimum pressure (max height) at which mid
                            ! level convection is allowed to initiate (Pa)
                            ! Original default 0.0
REAL :: amdet_fac    = rmdi ! Factor multiplying (1-rh) for adaptive mixing
                            ! detrainment rate.
                            ! Original default 1.0
! Scales with cca2d to determine convective cloud amount with anvil
REAL :: anvil_factor = rmdi ! x cca2d = max cca, cca is capped at 1.0
                            ! but this does not imply a cap to anvil_factor
                            ! Original default RMDI
REAL :: tower_factor = rmdi ! x cca2d = min cca
                            ! Original default RMDI
REAL :: ud_factor    = rmdi ! Updraught factor used in calculation of
                            ! convective water path Original default RMDI

REAL :: tice         = rmdi ! Phase change temperature in plume
                            ! Original default 273.15
REAL :: qstice       = rmdi ! Qsat at phase change temperature
                            ! (freeze/melting temperature)
                            ! Original default 3.5E-3

REAL :: t_melt_snow     = rmdi ! Temperature at which to melt all snow in
                               ! the downdraught.

! 5A & 6A code only
REAL :: ent_fac_dp      = rmdi ! Factor multiplying entrainment rate - deep
                               ! Original default 1.0
REAL :: ent_fac_md      = rmdi ! Factor multiplying entrainment rate - mid-level
                               ! Original default 1.0
REAL :: ent_dp_power    = rmdi ! Power n for (p/p*)^n for entrainment option 3
                               ! Original default 2.0
REAL :: ent_md_power    = rmdi ! Power n for (p/p*)^n for entrainment option 3
                               ! Original default 2.0
REAL :: cape_timescale  = rmdi ! Timescale for CAPE closure.
                               ! Original default RMDI
REAL :: cvdiag_sh_wtest = rmdi ! w for air above shallow convection must be
                               ! less than this. (Default value is 0.0)
                               ! Original default 0.0
! 6A code only
REAL :: eff_dcfl        = rmdi ! Factor that defines the efficiency by which
                               ! detrained liquid condensate creates liquid
                               ! cloud fraction
REAL :: eff_dcff        = rmdi ! Factor that defines the efficiency by which
                               ! detrained frozen condensate creates frozen
                               ! cloud fraction

! Convective prognostic options / variables
REAL :: tau_conv_prog_precip = rmdi
                               ! decay timescale for prognostic field of 
                               ! recent convective precipitation
! Factors used to calculate entrainment scaling from 3d prognostic field
! based on surface precipitation. Used under ent_opt_md,ent_opt_dp=6 or 7
REAL :: prog_ent_grad   = rmdi ! gradient term
REAL :: prog_ent_int    = rmdi ! intercept term
REAL :: prog_ent_max    = rmdi ! maximum scaling 
REAL :: prog_ent_min    = rmdi ! minimum scaling

!------------------------------------------------------------------------------
! Switches added for Warm TCS scheme
!------------------------------------------------------------------------------
! The following have been removed from the RUN_convection name list
!  - defaults are hard wired

INTEGER :: iwtcs_diag1 = 0  ! Options for WTCS diagnosis
INTEGER :: iwtcs_diag2 = 0  ! Options for WTCS diagnosis

!------------------------------------------------------------------------------
! Further changes

REAL, PARAMETER :: qmin_conv = 1.0e-8     ! Minimum allowed value of q after
   !  convection, also used for negative q check, and used in UKCA to ensure
   !  that stratospheric water vapour does not trigger messages from convection

!------------------------------------------------------------------------------
! Define namelist &Run_Convection read in from CNTLATM control file.
! Changes made to this list will affect both the Full UM and the SCM
!------------------------------------------------------------------------------

NAMELIST/Run_Convection/                                              &

i_convection_vn,                                                      &

! Logical switches
l_mom, l_mom_dd, l_fix_udfactor,                                      &
l_snow_rain,                                                          &
l_eman_dd, l_new_dd, l_use_dd, l_cloud_deep,                          &
          l_rediagnosis,  l_dcpl_cld4pc2,                             &
l_anvil,  l_murk_conv, l_safe_conv, l_cv_conserve_check,              &
l_ccrad,              l_3d_cca,               l_conv_hist,            &
l_param_conv,                                                         &
l_cmt_heating,                                                        &
l_conv_prog_group_1, l_conv_prog_group_2, l_conv_prog_group_3,        &
l_conv_prog_precip, l_jules_flux,                                     &

! General scheme options/variables
n_conv_calls,                                                         &
sh_pert_opt,                                                          &
dd_opt,               deep_cmt_opt,           mid_cmt_opt,            &
termconv,             adapt,                  fdet_opt,               &
r_det,                                                                &
tice,                 qstice,                                         &
ent_fac_dp,           ent_fac_md,                                     &
ent_opt_dp,           ent_opt_md,                                     &
ent_dp_power,         ent_md_power,                                   &
mdet_opt_dp,                                                          &
bl_cnv_mix,                                                           &
mid_cnv_pmin,         amdet_fac,              ccw_for_precip_opt,     &
cnv_wat_load_opt,     tv1_sd_opt,             limit_pert_opt,         &
cnv_cold_pools,                                                       &

! Convective diagnosis options
icvdiag,              plume_water_load,       dil_plume_water_load,   &
cvdiag_inv,           cvdiag_sh_wtest,                                &

! Closure & Cape related options/variables
cldbase_opt_dp,       cldbase_opt_md,         cldbase_opt_sh,         &
cape_bottom,          cape_top,               cape_timescale,         &
w_cape_limit,         cape_min,                                       &

! Convective cloud options/variables
cld_life_opt,         rad_cloud_decay_opt,    cca_min,                &
fixed_cld_life,       ud_factor,              mparwtr,                &
qlmin,                fac_qsat,               eff_dcfl,               &
eff_dcff,                                                             &

! CCRad options options/variables
cca2d_sh_opt,         cca_sh_knob,            ccw_sh_knob,            &
cca2d_md_opt,         cca_md_knob,            ccw_md_knob,            &
cca2d_dp_opt,         cca_dp_knob,            ccw_dp_knob,            &

! Anvil scheme options
anvil_factor,         tower_factor,           anv_opt,                &

! Downdraught control variables
t_melt_snow,                                                          &

! Convection type options
iconv_shallow,        iconv_mid,              iconv_deep,             &
iconv_congestus,                                                      &

! Convection prognostic options / variables
tau_conv_prog_precip,                                                 &

! Factors used to calculate entrainment scaling from 3d prognostic field
prog_ent_grad,        prog_ent_int,           prog_ent_max,           &
prog_ent_min


!------------------------------------------------------------------------------

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CV_RUN_MOD'

CONTAINS

SUBROUTINE check_run_convection()

USE chk_opts_mod, ONLY: chk_var, def_src
USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

INTEGER :: icode       ! error code

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_CONVECTION'
CHARACTER (LEN=errormessagelength) :: cmessage  ! used for ereport
CHARACTER (LEN=errormessagelength) :: mymessage ! Addtional error message

REAL(KIND=jprb)                    :: zhook_handle

!---------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------------
! required by chk_var routine
def_src = ModuleName//':'//RoutineName

! Which convection scheme 

CALL  chk_var(i_convection_vn,'i_convection_vn',                      &
              [i_convection_vn_5a,i_convection_vn_6a,                 &
               i_cv_llcs])

!---------------------------------------------------------------------------
! Check convective diagnosis options even when running without convection
!---------------------------------------------------------------------------

CALL chk_var(icvdiag,'icvdiag','[1:8]')
IF (i_convection_vn /= i_convection_vn_6a) THEN
  ! Check that not trying to use a 6A option
  CALL chk_var(icvdiag,'icvdiag','[1:5]')
END IF
CALL chk_var(cvdiag_inv,'cvdiag_inv','[0,1]')

CALL chk_var(cvdiag_sh_wtest,'cvdiag_sh_wtest','[-1.0:1.0]')

CALL chk_var(limit_pert_opt,'limit_pert_opt','[0,1,2]')
CALL chk_var(tv1_sd_opt,'tv1_sd_opt','[0,1,2]')

!---------------------------------------------------------------------------
! Only check the rest of the namelist if running with convection on
!---------------------------------------------------------------------------
IF (l_param_conv) THEN

  ! Recheck calling a convection scheme rather than just a convective
  ! diagnosis scheme  
  CALL  chk_var(i_convection_vn,'i_convection_vn',                      &
                   [i_convection_vn_5a,i_convection_vn_6a,i_cv_llcs])

  ! plume_water_loads are required for all convective schemes
  CALL chk_var(plume_water_load,'plume_water_load','[0,1,2]')
  CALL chk_var(dil_plume_water_load,'dil_plume_water_load','[0,1,2]')

  IF (i_convection_vn == i_convection_vn_5a .OR.                        &
      i_convection_vn == i_convection_vn_6a) THEN
  ! Sub-stepping - don't want to encourage a high value so now limiting 
  ! to 4
  mymessage ='Calling convection more than 4 times per timestep is not allowed'
  CALL chk_var(n_conv_calls,'n_conv_calls','[1:4]', cmessage=mymessage)

  ! Values only applying to 6A scheme
  IF (i_convection_vn == i_convection_vn_6a) THEN
    ! l_cmt_heating
    CALL chk_var(mdet_opt_dp,'mdet_opt_dp','[0:2]')
    CALL chk_var(fdet_opt,'fdet_opt','[0:3]')
    CALL chk_var(eff_dcff,'eff_dcff','[0.0:10.0]')
    CALL chk_var(eff_dcfl,'eff_dcfl','[0.0:10.0]')

    ! Convective memory 
    IF (l_conv_prog_precip) THEN
      ! Note meta data had no upper bound but checking routines require one
      ! so put in a day in seconds as not expecting to use values as large
      ! as this. 
      CALL chk_var(tau_conv_prog_precip,'tau_conv_prog_precip','[1.0:86400.0]')
      CALL chk_var(prog_ent_grad,'prog_ent_grad','[-3.0:0.0]')
      CALL chk_var(prog_ent_int,'prog_ent_int','[-5.0:5.0]')
      CALL chk_var(prog_ent_max,'prog_ent_max','[0.0:5.0]')
      CALL chk_var(prog_ent_min,'prog_ent_min','[0.0:5.0]')
      CALL chk_var(w_cape_limit,'w_cape_limit','[0.0:10000.0]')
    END IF

  END IF

  ! General convection panel
  ! Allowing setting of 1 to be used for testing mass flux congestus code
  CALL chk_var(iconv_congestus,'iconv_congestus','[0,1]')
  CALL chk_var(iconv_deep,'iconv_deep','[0:2]')
  CALL chk_var(iconv_mid,'iconv_mid','[0,1]')
  CALL chk_var(iconv_shallow,'iconv_shallow','[0:3]')

  CALL chk_var(bl_cnv_mix,'bl_cnv_mix','[0,1]')
  CALL chk_var(cnv_wat_load_opt,'cnv_wat_load_opt','[0,1]')
  CALL chk_var(mid_cnv_pmin,'mid_cnv_pmin','[0.0:25000.0]')

  CALL chk_var(qlmin,'qlmin','[0.0:1.0e-3]')
  CALL chk_var(qstice,'qstice','[0.0:0.5]')
  CALL chk_var(tice,'tice','[200.0:300.0]')

  ! Updraughts panel

  CALL chk_var(sh_pert_opt,'sh_pert_opt','[0,1]')
  CALL chk_var(termconv,'termconv','[0,1]')

  CALL chk_var(ccw_for_precip_opt,'ccw_for_precip_opt','[0:4]')
  CALL chk_var(mparwtr,'mparwtr','[1.0e-5:1.0e-2]')
  IF (ccw_for_precip_opt == 4) THEN
    CALL chk_var(fac_qsat,'fac_qsat','[0.0:5.0]')
  END IF
  CALL chk_var(ent_opt_dp,'ent_opt_dp','[0:9]')
  CALL chk_var(ent_fac_dp,'ent_fac_dp','[0.0:2.0]')
  IF (ent_opt_dp == 3) THEN
    CALL chk_var(ent_dp_power,'ent_dp_power','[0.0:10.0]')
  END IF 

  CALL chk_var(ent_opt_md,'ent_opt_md','[0,1,2,3,6,7]')
  CALL chk_var(ent_fac_md,'ent_fac_md','[0.0:2.0]')
  IF (ent_opt_md == 3) THEN
    CALL chk_var(ent_md_power,'ent_md_power','[0.0:10.0]')
  END IF 

  CALL chk_var(adapt,'adapt','[0:8]')
  CALL chk_var(r_det,'r_det','[0.0:1.0]')
  CALL chk_var(amdet_fac,'amdet_fac','[0.0:20.0]')

  ! Convective closure checking

  CALL chk_var(cldbase_opt_dp,'cldbase_opt_dp','[0,1,2,3,4,6,7,8,9]')
  CALL chk_var(cldbase_opt_md,'cldbase_opt_md','[0,1,2,3,4,6,7]')
  CALL chk_var(cldbase_opt_sh,'cldbase_opt_sh','[0,1]')
  CALL chk_var(cape_timescale,'cape_timescale','[1.0:9999999.0]')
  
  ! Variables only need checking for certain settings
  IF (cldbase_opt_dp == 3 .OR. cldbase_opt_dp == 4 .OR. &
      cldbase_opt_dp == 6 .OR. cldbase_opt_md == 3 .OR. &
      cldbase_opt_md == 4 .OR. cldbase_opt_md == 6) THEN
    CALL chk_var(cape_bottom,'cape_bottom','[1:10000]')
    CALL chk_var(cape_top,'cape_top','[2:10000]')
    ! Also check cape_top > cape_bottom or force error
    IF (cape_bottom >= cape_top) THEN
      cmessage='cape_top must be > cape_bottom'
      icode = 1
      CALL ereport(RoutineName,icode,cmessage)
    END IF   
    CALL chk_var(w_cape_limit,'w_cape_limit','[0.0:10000.0]')
  END IF
  IF (cldbase_opt_dp == 4 .OR. cldbase_opt_md == 4 ) THEN
    CALL chk_var(cape_min,'cape_min','[0.0:2000.0]')
  END IF

  ! Downdraughts
  CALL chk_var(dd_opt,'dd_opt','[0,1]')

  IF (l_snow_rain .OR. l_new_dd) THEN  
    CALL chk_var(t_melt_snow,'t_melt_snow','[273.15:300.0]')
  END IF

  ! CMT checking
  IF (l_mom) THEN
    CALL chk_var(deep_cmt_opt,'deep_cmt_opt','[0:6]')
    CALL chk_var(mid_cmt_opt,'mid_cmt_opt','[0:2]')

    ! Note l_mom_dd will only work with new DD scheme l_new_dd 
  END IF

  ! Cloud panel checks
  IF (l_3d_cca) THEN

    IF (l_anvil) THEN
      CALL chk_var(anv_opt,'anv_opt','[0,1,2]')
      CALL chk_var(anvil_factor,'anvil_factor','[0.0:3.0]')
      CALL chk_var(tower_factor,'tower_factor','[0.0:3.0]')
    END IF

    IF (l_ccrad) THEN
      CALL chk_var(cca2d_dp_opt,'cca2d_dp_opt','[0,1]')
      CALL chk_var(cca2d_md_opt,'cca2d_md_opt','[0,1]')
      CALL chk_var(cca2d_sh_opt,'cca2d_sh_opt','[0,1,2]')
      CALL chk_var(cca_dp_knob,'cca_dp_knob','[0.0:10.0]')
      CALL chk_var(cca_md_knob,'cca_md_knob','[0.0:10.0]')
      CALL chk_var(cca_sh_knob,'cca_sh_knob','[0.0:10.0]')
      CALL chk_var(ccw_dp_knob,'cca_dp_knob','[0.0:10.0]')
      CALL chk_var(ccw_md_knob,'cca_md_knob','[0.0:10.0]')
      CALL chk_var(ccw_sh_knob,'cca_sh_knob','[0.0:10.0]')
    END IF
  END IF

  CALL chk_var(rad_cloud_decay_opt,'rad_cloud_decay_opt','[0,1,2]')

  IF (rad_cloud_decay_opt==1 .OR. rad_cloud_decay_opt==2) THEN
    CALL chk_var(cld_life_opt,'cld_life_opt','[0,1]')
    CALL chk_var(cca_min,'cca_min','[0.001:0.4]')
    CALL chk_var(fixed_cld_life,'fixed_cld_life','[0.0:99999.0]')
  END IF

  ! Updraught factor only used for these switches
  IF (l_fix_udfactor .OR. (ccw_for_precip_opt == 0) .OR.     &
     (ccw_for_precip_opt == 1) ) THEN
    CALL chk_var(ud_factor,'ud_factor','[0.0:1.0]')
  END IF

  END IF ! 5a or 6a scheme checking
 
ELSE  ! Not planning to call a convection scheme so values should be
      ! trigger ignored which means logical must be false.
  ! Check no logicals set to true - at present prints warnings but should really
  ! be stopping the run.
  IF (l_mom) THEN
    icode = -1
    WRITE(cmessage,'(a29,a52)')' l_mom  SHOULD BE set .false.',        &
         ' as not calling convection - run will treat as false'
    CALL ereport(RoutineName,icode,cmessage)
  END IF
  IF (l_3d_cca) THEN
    icode = -2
    WRITE(cmessage,'(a32,a43)')' l_3d_cca  SHOULD BE set .false.',     &
         ' as not calling convection so wasting space'
    CALL ereport(RoutineName,icode,cmessage)
  END IF
  IF (l_ccrad) THEN
    icode = -3
    WRITE(cmessage,'(a31,a43)')' l_ccrad  SHOULD BE set .false.',      &
         ' as not calling convection so wasting space'
    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF

! Reset
def_src = ''

! Further checking
IF (L_ccrad) THEN

  IF (.NOT. l_3d_cca) THEN
    icode = 100
    CMessage    = '**ERROR**: CCRad is not yet available without'//           &
                            ' the anvil scheme (L_3D_CCA = .TRUE.)'

    CALL ereport(RoutineName, icode, CMessage)
  END IF

  IF (l_fix_udfactor) THEN
    icode = 101
    CMessage    = '**ERROR**: L_CCRad and l_fix_udfactor'//                   &
                            ' should not be both set to true.'

    CALL ereport(RoutineName, icode, CMessage)
  END IF

  IF (l_pc2_diag_sh) THEN
    icode = 102
    CMessage    = '**ERROR**: L_CCRad and l_pc2_diag_sh'//                    &
                            ' should not be both set to true.'

    CALL ereport(RoutineName, icode, CMessage)
  END IF


END IF
!

!---------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!---------------------------------------------------------------------------
RETURN
END SUBROUTINE check_run_convection


SUBROUTINE print_nlist_run_convection()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_CONVECTION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_convection', &
    src='cv_run_mod')

WRITE(lineBuffer,*)' i_convection_vn = ',i_convection_vn
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_mom = ',l_mom
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_mom_dd = ',l_mom_dd
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_fix_udfactor = ',l_fix_udfactor
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_snow_rain = ',l_snow_rain
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_eman_dd = ',l_eman_dd
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_new_dd = ',l_new_dd
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_use_dd = ',l_use_dd
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_cloud_deep = ',l_cloud_deep
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_rediagnosis = ',l_rediagnosis
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_dcpl_cld4pc2 = ',l_dcpl_cld4pc2
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_anvil = ',l_anvil
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_murk_conv = ',l_murk_conv
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_safe_conv = ',l_safe_conv
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_ccrad = ',l_ccrad
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_3d_cca = ',l_3d_cca
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_conv_hist = ',l_conv_hist
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_param_conv = ',l_param_conv
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_cv_conserve_check = ',l_cv_conserve_check
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_cmt_heating = ',l_cmt_heating
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_conv_prog_group_1 = ',l_conv_prog_group_1
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_conv_prog_group_2 = ',l_conv_prog_group_2
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_conv_prog_group_3 = ',l_conv_prog_group_3
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' l_conv_prog_precip = ',l_conv_prog_precip
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(A,L1)')' l_jules_flux = ',l_jules_flux
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' n_conv_calls = ',n_conv_calls
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' sh_pert_opt = ',sh_pert_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' dd_opt = ',dd_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' deep_cmt_opt = ',deep_cmt_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' mid_cmt_opt = ',mid_cmt_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' termconv = ',termconv
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' adapt = ',adapt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' fdet_opt = ',fdet_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' r_det = ',r_det
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' tice = ',tice
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' qstice = ',qstice
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_fac_dp = ',ent_fac_dp
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_fac_md = ',ent_fac_md
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_opt_dp = ',ent_opt_dp
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_opt_md = ',ent_opt_md
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_dp_power = ',ent_dp_power
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ent_md_power = ',ent_md_power
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' bl_cnv_mix = ',bl_cnv_mix
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' mid_cnv_pmin = ',mid_cnv_pmin
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' amdet_fac = ',amdet_fac
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ccw_for_precip_opt = ',ccw_for_precip_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cnv_wat_load_opt = ',cnv_wat_load_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' tv1_sd_opt = ',tv1_sd_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' limit_pert_opt = ',limit_pert_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' icvdiag = ',icvdiag
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' plume_water_load = ',plume_water_load
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' dil_plume_water_load = ',dil_plume_water_load
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cvdiag_inv = ',cvdiag_inv
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cvdiag_sh_wtest = ',cvdiag_sh_wtest
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cldbase_opt_dp = ',cldbase_opt_dp
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cldbase_opt_md = ',cldbase_opt_md
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cldbase_opt_sh = ',cldbase_opt_sh
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cape_bottom = ',cape_bottom
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cape_top = ',cape_top
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cape_timescale = ',cape_timescale
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' w_cape_limit = ',w_cape_limit
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cape_min = ',cape_min
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cld_life_opt = ',cld_life_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' rad_cloud_decay_opt = ',rad_cloud_decay_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca_min = ',cca_min
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' fixed_cld_life = ',fixed_cld_life
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ud_factor = ',ud_factor
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' mparwtr = ',mparwtr
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' qlmin = ',qlmin
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' fac_qsat = ',fac_qsat
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' eff_dcfl = ',eff_dcfl
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' eff_dcff = ',eff_dcff
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca2d_sh_opt = ',cca2d_sh_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca_sh_knob = ',cca_sh_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ccw_sh_knob = ',ccw_sh_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca2d_md_opt = ',cca2d_md_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca_md_knob = ',cca_md_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ccw_md_knob = ',ccw_md_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca2d_dp_opt = ',cca2d_dp_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' cca_dp_knob = ',cca_dp_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' ccw_dp_knob = ',ccw_dp_knob
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' anvil_factor = ',anvil_factor
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' tower_factor = ',tower_factor
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' anv_opt = ',anv_opt
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' t_melt_snow = ',t_melt_snow
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' iconv_shallow = ',iconv_shallow
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' iconv_mid = ',iconv_mid
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' iconv_deep = ',iconv_deep
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' iconv_congestus = ',iconv_congestus
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,*)' tau_conv_prog_precip = ',tau_conv_prog_precip
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(a17,f15.3)')' prog_ent_grad = ',prog_ent_grad
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(a17,f15.3)')' prog_ent_int =  ',prog_ent_int
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(a17,f15.3)')' prog_ent_max =  ',prog_ent_max
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(a17,f15.3)')' prog_ent_min =  ',prog_ent_min
CALL umPrint(lineBuffer,src='cv_run_mod')
WRITE(lineBuffer,'(A,I0)')' cnv_cold_pools =  ',cnv_cold_pools
CALL umPrint(lineBuffer,src='cv_run_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='cv_run_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_convection

#if !defined(LFRIC)
SUBROUTINE read_nml_run_convection(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_CONVECTION'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 37
INTEGER, PARAMETER :: n_real = 35
INTEGER, PARAMETER :: n_log = 24

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_convection_vn
  INTEGER :: n_conv_calls
  INTEGER :: sh_pert_opt
  INTEGER :: dd_opt
  INTEGER :: deep_cmt_opt
  INTEGER :: mid_cmt_opt
  INTEGER :: termconv
  INTEGER :: adapt
  INTEGER :: fdet_opt
  INTEGER :: ent_opt_dp
  INTEGER :: ent_opt_md
  INTEGER :: mdet_opt_dp
  INTEGER :: bl_cnv_mix
  INTEGER :: ccw_for_precip_opt
  INTEGER :: cnv_wat_load_opt
  INTEGER :: tv1_sd_opt
  INTEGER :: limit_pert_opt
  INTEGER :: icvdiag
  INTEGER :: plume_water_load
  INTEGER :: dil_plume_water_load
  INTEGER :: cvdiag_inv
  INTEGER :: cldbase_opt_dp
  INTEGER :: cldbase_opt_md
  INTEGER :: cldbase_opt_sh
  INTEGER :: cape_bottom
  INTEGER :: cape_top
  INTEGER :: cld_life_opt
  INTEGER :: rad_cloud_decay_opt
  INTEGER :: cca2d_sh_opt
  INTEGER :: cca2d_md_opt
  INTEGER :: cca2d_dp_opt
  INTEGER :: anv_opt
  INTEGER :: iconv_shallow
  INTEGER :: iconv_mid
  INTEGER :: iconv_deep
  INTEGER :: iconv_congestus
  INTEGER :: cnv_cold_pools
  REAL :: r_det
  REAL :: tice
  REAL :: qstice
  REAL :: ent_fac_dp
  REAL :: ent_fac_md
  REAL :: ent_dp_power
  REAL :: ent_md_power
  REAL :: mid_cnv_pmin
  REAL :: amdet_fac
  REAL :: cvdiag_sh_wtest
  REAL :: cape_timescale
  REAL :: w_cape_limit
  REAL :: cape_min
  REAL :: cca_min
  REAL :: fixed_cld_life
  REAL :: ud_factor
  REAL :: mparwtr
  REAL :: qlmin
  REAL :: fac_qsat
  REAL :: eff_dcfl
  REAL :: eff_dcff
  REAL :: cca_sh_knob
  REAL :: ccw_sh_knob
  REAL :: cca_md_knob
  REAL :: ccw_md_knob
  REAL :: cca_dp_knob
  REAL :: ccw_dp_knob
  REAL :: anvil_factor
  REAL :: tower_factor
  REAL :: t_melt_snow
  REAL :: tau_conv_prog_precip
  REAL :: prog_ent_grad
  REAL :: prog_ent_int
  REAL :: prog_ent_max
  REAL :: prog_ent_min
  LOGICAL :: l_mom
  LOGICAL :: l_mom_dd
  LOGICAL :: l_fix_udfactor
  LOGICAL :: l_snow_rain
  LOGICAL :: l_eman_dd
  LOGICAL :: l_new_dd
  LOGICAL :: l_use_dd
  LOGICAL :: l_cloud_deep
  LOGICAL :: l_rediagnosis
  LOGICAL :: l_dcpl_cld4pc2
  LOGICAL :: l_anvil
  LOGICAL :: l_murk_conv
  LOGICAL :: l_safe_conv
  LOGICAL :: l_cv_conserve_check
  LOGICAL :: l_ccrad
  LOGICAL :: l_3d_cca
  LOGICAL :: l_conv_hist
  LOGICAL :: l_param_conv
  LOGICAL :: l_cmt_heating
  LOGICAL :: l_conv_prog_group_1
  LOGICAL :: l_conv_prog_group_2
  LOGICAL :: l_conv_prog_group_3
  LOGICAL :: l_conv_prog_precip
  LOGICAL :: l_jules_flux
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Convection, IOSTAT=ErrorStatus,       &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Convection", iomessage)

  my_nml % i_convection_vn    = i_convection_vn
  my_nml % n_conv_calls       = n_conv_calls
  my_nml % sh_pert_opt        = sh_pert_opt
  my_nml % dd_opt             = dd_opt
  my_nml % deep_cmt_opt       = deep_cmt_opt
  my_nml % mid_cmt_opt        = mid_cmt_opt
  my_nml % termconv           = termconv
  my_nml % adapt              = adapt
  my_nml % fdet_opt           = fdet_opt
  my_nml % ent_opt_dp         = ent_opt_dp
  my_nml % ent_opt_md         = ent_opt_md
  my_nml % mdet_opt_dp        = mdet_opt_dp
  my_nml % bl_cnv_mix         = bl_cnv_mix
  my_nml % ccw_for_precip_opt = ccw_for_precip_opt
  my_nml % cnv_wat_load_opt   = cnv_wat_load_opt
  my_nml % tv1_sd_opt         = tv1_sd_opt
  my_nml % limit_pert_opt     = limit_pert_opt
  my_nml % icvdiag            = icvdiag
  my_nml % plume_water_load   = plume_water_load
  my_nml % dil_plume_water_load = dil_plume_water_load
  my_nml % cvdiag_inv         = cvdiag_inv
  my_nml % cldbase_opt_dp     = cldbase_opt_dp
  my_nml % cldbase_opt_md     = cldbase_opt_md
  my_nml % cldbase_opt_sh     = cldbase_opt_sh
  my_nml % cape_bottom        = cape_bottom
  my_nml % cape_top           = cape_top
  my_nml % cld_life_opt       = cld_life_opt
  my_nml % rad_cloud_decay_opt = rad_cloud_decay_opt
  my_nml % cca2d_sh_opt       = cca2d_sh_opt
  my_nml % cca2d_md_opt       = cca2d_md_opt
  my_nml % cca2d_dp_opt       = cca2d_dp_opt
  my_nml % anv_opt            = anv_opt
  my_nml % iconv_shallow      = iconv_shallow
  my_nml % iconv_mid          = iconv_mid
  my_nml % iconv_deep         = iconv_deep
  my_nml % iconv_congestus    = iconv_congestus
  my_nml % cnv_cold_pools     = cnv_cold_pools
  ! end of integers
  my_nml % r_det           = r_det
  my_nml % tice            = tice
  my_nml % qstice          = qstice
  my_nml % ent_fac_dp      = ent_fac_dp
  my_nml % ent_fac_md      = ent_fac_md
  my_nml % ent_dp_power    = ent_dp_power
  my_nml % ent_md_power    = ent_md_power
  my_nml % mid_cnv_pmin    = mid_cnv_pmin
  my_nml % amdet_fac       = amdet_fac
  my_nml % cvdiag_sh_wtest = cvdiag_sh_wtest
  my_nml % cape_timescale  = cape_timescale
  my_nml % w_cape_limit    = w_cape_limit
  my_nml % cape_min        = cape_min
  my_nml % cca_min         = cca_min
  my_nml % fixed_cld_life  = fixed_cld_life
  my_nml % ud_factor       = ud_factor
  my_nml % mparwtr         = mparwtr
  my_nml % qlmin           = qlmin
  my_nml % fac_qsat        = fac_qsat
  my_nml % eff_dcfl        = eff_dcfl
  my_nml % eff_dcff        = eff_dcff
  my_nml % cca_sh_knob     = cca_sh_knob
  my_nml % ccw_sh_knob     = ccw_sh_knob
  my_nml % cca_md_knob     = cca_md_knob
  my_nml % ccw_md_knob     = ccw_md_knob
  my_nml % cca_dp_knob     = cca_dp_knob
  my_nml % ccw_dp_knob     = ccw_dp_knob
  my_nml % anvil_factor    = anvil_factor
  my_nml % tower_factor    = tower_factor
  my_nml % t_melt_snow     = t_melt_snow
  my_nml % tau_conv_prog_precip = tau_conv_prog_precip
  my_nml % prog_ent_grad   = prog_ent_grad
  my_nml % prog_ent_int    = prog_ent_int
  my_nml % prog_ent_max    = prog_ent_max
  my_nml % prog_ent_min    = prog_ent_min
  ! end of reals
  my_nml % l_mom          = l_mom
  my_nml % l_mom_dd       = l_mom_dd
  my_nml % l_fix_udfactor = l_fix_udfactor
  my_nml % l_snow_rain    = l_snow_rain
  my_nml % l_eman_dd      = l_eman_dd
  my_nml % l_new_dd       = l_new_dd
  my_nml % l_use_dd       = l_use_dd
  my_nml % l_cloud_deep   = l_cloud_deep
  my_nml % l_rediagnosis  = l_rediagnosis
  my_nml % l_dcpl_cld4pc2 = l_dcpl_cld4pc2
  my_nml % l_anvil        = l_anvil
  my_nml % l_murk_conv    = l_murk_conv
  my_nml % l_safe_conv    = l_safe_conv
  my_nml % l_cv_conserve_check = l_cv_conserve_check
  my_nml % l_ccrad        = l_ccrad
  my_nml % l_3d_cca       = l_3d_cca
  my_nml % l_conv_hist    = l_conv_hist
  my_nml % l_param_conv   = l_param_conv
  my_nml % l_cmt_heating  = l_cmt_heating
  my_nml % l_conv_prog_group_1 = l_conv_prog_group_1
  my_nml % l_conv_prog_group_2 = l_conv_prog_group_2
  my_nml % l_conv_prog_group_3 = l_conv_prog_group_3
  my_nml % l_conv_prog_precip  = l_conv_prog_precip
  my_nml % l_jules_flux = l_jules_flux

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_convection_vn    = my_nml % i_convection_vn
  n_conv_calls       = my_nml % n_conv_calls
  sh_pert_opt        = my_nml % sh_pert_opt
  dd_opt             = my_nml % dd_opt
  deep_cmt_opt       = my_nml % deep_cmt_opt
  mid_cmt_opt        = my_nml % mid_cmt_opt
  termconv           = my_nml % termconv
  adapt              = my_nml % adapt
  fdet_opt           = my_nml % fdet_opt
  ent_opt_dp         = my_nml % ent_opt_dp
  ent_opt_md         = my_nml % ent_opt_md
  mdet_opt_dp        = my_nml % mdet_opt_dp
  bl_cnv_mix         = my_nml % bl_cnv_mix
  ccw_for_precip_opt = my_nml % ccw_for_precip_opt
  cnv_wat_load_opt   = my_nml % cnv_wat_load_opt
  tv1_sd_opt         = my_nml % tv1_sd_opt
  limit_pert_opt     = my_nml % limit_pert_opt
  icvdiag            = my_nml % icvdiag
  plume_water_load   = my_nml % plume_water_load
  dil_plume_water_load = my_nml % dil_plume_water_load
  cvdiag_inv         = my_nml % cvdiag_inv
  cldbase_opt_dp     = my_nml % cldbase_opt_dp
  cldbase_opt_md     = my_nml % cldbase_opt_md
  cldbase_opt_sh     = my_nml % cldbase_opt_sh
  cape_bottom        = my_nml % cape_bottom
  cape_top           = my_nml % cape_top
  cld_life_opt       = my_nml % cld_life_opt
  rad_cloud_decay_opt = my_nml % rad_cloud_decay_opt
  cca2d_sh_opt       = my_nml % cca2d_sh_opt
  cca2d_md_opt       = my_nml % cca2d_md_opt
  cca2d_dp_opt       = my_nml % cca2d_dp_opt
  anv_opt            = my_nml % anv_opt
  iconv_shallow      = my_nml % iconv_shallow
  iconv_mid          = my_nml % iconv_mid
  iconv_deep         = my_nml % iconv_deep
  iconv_congestus    = my_nml % iconv_congestus
  cnv_cold_pools     = my_nml % cnv_cold_pools
  ! end of integers
  r_det           = my_nml % r_det
  tice            = my_nml % tice
  qstice          = my_nml % qstice
  ent_fac_dp      = my_nml % ent_fac_dp
  ent_fac_md      = my_nml % ent_fac_md
  ent_dp_power    = my_nml % ent_dp_power
  ent_md_power    = my_nml % ent_md_power
  mid_cnv_pmin    = my_nml % mid_cnv_pmin
  amdet_fac       = my_nml % amdet_fac
  cvdiag_sh_wtest = my_nml % cvdiag_sh_wtest
  cape_timescale  = my_nml % cape_timescale
  w_cape_limit    = my_nml % w_cape_limit
  cape_min        = my_nml % cape_min
  cca_min         = my_nml % cca_min
  fixed_cld_life  = my_nml % fixed_cld_life
  ud_factor       = my_nml % ud_factor
  mparwtr         = my_nml % mparwtr
  qlmin           = my_nml % qlmin
  fac_qsat        = my_nml % fac_qsat
  eff_dcfl        = my_nml % eff_dcfl
  eff_dcff        = my_nml % eff_dcff
  cca_sh_knob     = my_nml % cca_sh_knob
  ccw_sh_knob     = my_nml % ccw_sh_knob
  cca_md_knob     = my_nml % cca_md_knob
  ccw_md_knob     = my_nml % ccw_md_knob
  cca_dp_knob     = my_nml % cca_dp_knob
  ccw_dp_knob     = my_nml % ccw_dp_knob
  anvil_factor    = my_nml % anvil_factor
  tower_factor    = my_nml % tower_factor
  t_melt_snow     = my_nml % t_melt_snow
  tau_conv_prog_precip = my_nml % tau_conv_prog_precip
  prog_ent_grad   = my_nml % prog_ent_grad
  prog_ent_int    = my_nml % prog_ent_int
  prog_ent_max    = my_nml % prog_ent_max
  prog_ent_min    = my_nml % prog_ent_min
  ! end of reals
  l_mom               = my_nml % l_mom
  l_mom_dd            = my_nml % l_mom_dd
  l_fix_udfactor      = my_nml % l_fix_udfactor
  l_snow_rain         = my_nml % l_snow_rain
  l_eman_dd           = my_nml % l_eman_dd
  l_new_dd            = my_nml % l_new_dd
  l_use_dd            = my_nml % l_use_dd
  l_cloud_deep        = my_nml % l_cloud_deep
  l_rediagnosis       = my_nml % l_rediagnosis
  l_dcpl_cld4pc2      = my_nml % l_dcpl_cld4pc2
  l_anvil             = my_nml % l_anvil
  l_murk_conv         = my_nml % l_murk_conv
  l_safe_conv         = my_nml % l_safe_conv
  l_cv_conserve_check = my_nml % l_cv_conserve_check
  l_ccrad             = my_nml % l_ccrad
  l_3d_cca            = my_nml % l_3d_cca
  l_conv_hist         = my_nml % l_conv_hist
  l_param_conv        = my_nml % l_param_conv
  l_cmt_heating       = my_nml % l_cmt_heating
  l_conv_prog_group_1 = my_nml % l_conv_prog_group_1
  l_conv_prog_group_2 = my_nml % l_conv_prog_group_2
  l_conv_prog_group_3 = my_nml % l_conv_prog_group_3
  l_conv_prog_precip  = my_nml % l_conv_prog_precip
  l_jules_flux        = my_nml % l_jules_flux

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_convection
#endif

END MODULE cv_run_mod
