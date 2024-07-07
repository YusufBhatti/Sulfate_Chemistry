! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE stochastic_physics_run_mod

!  Global data module for switches/options concerned with stochastic physics
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Stochastic Physics

  ! Description:
  !   Module containing runtime logicals/options used by the SKEB2
  !   and RP2 schemes.
  !   and the stochastic_physics_run_setup subroutine to control logic of
  !   selected options.

  ! Method:
  !   All switches/options which are contained in the &RUN_Stochastic
  !   sub-namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate, but mostly
  !   rmdi values are set to force users to specify values through Rose.

  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.

  !   Note, flags for other stochastic physics routines, which were
  !   embedded within other namelists, have now been moved here.

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength

! npft_max is defined in Jules code which the Utils do not depend on.
! fieldcalc uses this module, but not the arrays being declared as size
! npft_max. So rather than include a dependency on Jules in the Utils, the
! requirement to include npft_max is excluded from Utils builds.
#if !defined(UTILIO)
USE max_dimensions, ONLY: npft_max
USE jules_vegetation_mod, ONLY: l_triffid, l_phenol
#endif

IMPLICIT NONE
!======================================================================
! Logical switches set from RUN_Stochastic namelist via UI
!======================================================================

LOGICAL :: l_skeb2  = .FALSE.                                           &
           ! switch controls use of SKEB2 scheme
,           l_skeb2_psicdisp = .FALSE.                                  &
           ! TRUE = incl streamfunction modulation by convection
,           l_skeb2_psisdisp = .FALSE.                                  &
           ! TRUE = incl streamfunction modulation by SL-advection
,           l_skeb2_skeb1disp = .FALSE.                                 &
           ! TRUE = incl streamfunction modulation by SKEB1-type
,           l_skeb2_velpot = .FALSE.                                    &
           ! TRUE = Calc divergent-wind incr from VelPot forcing
,           l_rp2 = .FALSE.                                             &
           ! switch controls use of RP2 scheme
,           l_skebsmooth_adv = .FALSE.                                  &
           ! TRUE = Perform advanced smoothing of energy diss fields
,           l_skebprint = .FALSE.                                       &
           ! TRUE = Print global KE backscattered at each timestep
,           l_skeb2_conv_disp_mod=.FALSE.                               &
           ! TRUE= A horizontal resolution factor modulates the
           ! convective dissipation rate.
          
           
,           l_x_eq_sin_x= .FALSE.                                       &
           ! TRUE= Activates the approximation sin(x)=x for the
           ! computation of Lengrende Polynomials
,           l_rp2_cycle_out = .FALSE.                                   &
           ! TRUE = Write RP parameters to a file after 12 hours
,           l_rp2_cycle_in = .FALSE.                                    &
           ! TRUE = Initialise RP parameters from a specified file

!=========================================================================
! Logical switches for Stochastic Peturbation of Tendencies scheme (SPT)
!    (TRUE: Included, FALSE: Not included)
!=========================================================================

 ,           l_spt = .FALSE.                                            &
            ! Switch controls use of SPT scheme
 ,           l_spt_cfl=.FALSE.                                          &
             ! Switch to activate the CFL filtering to those
             ! points that breach CFL criteria
 ,           l_spt_rain = .FALSE.                                       &
             ! Switch to add LS rain (microphysics) tendency to SPT
 ,           l_spt_rad= .FALSE.                                         &
             ! Switch to add Radiation tendency to SPT
 ,           l_spt_gwd= .FALSE.                                         &
             ! Switch to add Gravity wave drag tendency to SPT
 ,           l_spt_conv= .FALSE.                                        &
             ! Switch to add convection tendency to SPT
 ,           l_spt_conv_mom= .FALSE.                                    &
             ! Switch to add tendency from convective momentum
             ! transport to SPT
 ,           l_spt_qcons = .FALSE.                                      &
             ! Conserve total water in the column
 ,           l_spt_mse = .FALSE.                                        
             ! Conserve moist static energy 
             
! Define numerical and convective scheme type, only one can be active
! at any time.
INTEGER :: type_cape = 4 ! Calc Conv Dissipation using CAPE
INTEGER :: type_mflx = 5 ! Calc Conv Dissipation using Mass-Flux
INTEGER :: skeb2_cdisp = imdi ! type of streamfn modulation by conv
INTEGER :: type_smag = 6 ! Calc num dissipation using Smagorinsky
INTEGER :: type_bihm = 7 ! Calc num dissipation using bi-harmonic
INTEGER :: skeb2_sdisp = imdi ! type of numerical dissipation scheme

!======================================================================
! Integer parameters to identify type of model run
!======================================================================
INTEGER, PARAMETER :: firsttimestep_true  = 1
         ! This is the first timestep in a new run (NRUN)
INTEGER, PARAMETER :: firsttimestep_false = 2
         ! This is after the first timestep
INTEGER, PARAMETER :: firsttimestep_crun  = -1
         ! This is the first timestep in a continuation run (CRUN)

!======================================================================
! Integers to indicate which stochastic physics information is
! output to the dump header
!======================================================================

INTEGER, PARAMETER :: stph_seed_present          = 1
INTEGER, PARAMETER :: stph_spt_data_present      = 2
INTEGER, PARAMETER :: stph_skeb2_data_present    = 4
INTEGER, PARAMETER :: stph_rp2_data_present      = 8

INTEGER :: stph_header_flag = stph_seed_present
           ! Flag to pass to ih_stochastic_flag indicating which
           ! stochastic physics information is output to the dump header.
           ! For each stochastic physics scheme in use, the parameter
           ! stph_<scheme>_data_present is added to stph_header_flag,
           ! resulting in a single integer value that can be used to determine
           ! which schemes have data stored in the dump header.
           ! It is initialised to stph_seed_present since the stochastic
           ! seed is stored if any stochastic physics scheme is in use.
           ! Data for the stochastic BL perturbation scheme is stored as
           ! a prognostic so is not included here.
INTEGER :: stph_spt_data_check = 0
INTEGER :: stph_skeb2_data_check = 0
INTEGER :: stph_rp2_data_check = 0

INTEGER, PARAMETER :: rp_max = 25  ! Max number of RP fields
                                   ! Used to set the size of the array
                                   ! for storing RP values in the dump

#if defined(UTILIO)
INTEGER, PARAMETER :: npft_max = 1
LOGICAL :: l_triffid = .FALSE.
LOGICAL :: l_phenol = .FALSE.
#endif

!======================================================================
! Integer options set from RUN_Stochastic namelist via UI
!======================================================================
!
!   stph_n1 and stph_n2 define the range of spherical harmonic orders over
!     which backscatter is applied (remember: N(N+1)/R*R is the
!     effective wavenumber squared). For ENDGame the suggested range
!     is [20;60] - horizontal wavelengths in the range 500 km to 2000 km.
!
!   skeb2_toplev defines the top level of the streamfunction modulation
!     field. This will change with vertical resolution. Typically it is
!     set to a level near 18km to avoid making changes above the stratosphere
!   skeb2_botlev defines the bottom level of the streamfunction
!     modulating field (we may want to avoid changes in the boundary
!     layer, certainly in level one)
!   nsmooth is the number of 1-2-1 spatial smoothing iterations (5)
!   rhcrit_ref_level is usually level 3

INTEGER :: stphseed = imdi       ! Control variable for random seed options:
           ! 0 => Use ensemble member and date/time of dump
           ! 1 => Use date/time from computer clock
           ! 2 => Use seed/coefficient values stored in dump file

INTEGER :: stph_n1=imdi          ! minimum wavenumber for backscatter
INTEGER :: stph_n2=imdi          ! maximum wavenumber for backscatter
INTEGER :: rp2_callfreq = imdi   ! RP2 calling frequency (in seconds)
INTEGER :: rp2_cycle_tm = imdi   ! time at which to write out RP2 parameters
                                 ! (in seconds) when l_rp2_cycle_out = TRUE

INTEGER :: i_rp_scheme = imdi    ! Switch to specify RP scheme option:
           ! 0 => RP2 scheme
INTEGER, PARAMETER :: i_rp2 = 0
           ! 1 => RP2b scheme (updated algorithm and additional parameters)
INTEGER, PARAMETER :: i_rp2b = 1

INTEGER :: skeb2_toplev=imdi     ! Top level of SKEB2 calculations
INTEGER :: skeb2_botlev=imdi     ! Bottom level of SKEB2 calculations
INTEGER :: nsmooth = imdi        ! Iteration count for spatial smoothing
INTEGER :: rhcrit_ref_level=imdi ! RHCrit reference level for RP2
INTEGER :: offx_stph = imdi      ! Halo size used in spatial smoothing
INTEGER :: offy_stph = imdi      ! Halo size used in spatial smoothing
                                 ! Both set to nsmooth in sthp_setup
INTEGER :: ran_max = imdi        ! number of independent AR1 processes in RP2
INTEGER :: ran_count = imdi      ! counter for random number array in RP
INTEGER :: stph_nens = imdi      ! ensemble member number
CHARACTER(LEN=8) :: ens_member   ! ensemble member (8-character string)

! ---- SPT integers ---
INTEGER :: spt_top_tap_lev= imdi ! Top level of spt calculations.
                                 ! Including tapering zone
INTEGER :: spt_toplev= imdi      ! Top level of spt calculations.
                                 ! Excluding tapering zone
INTEGER :: spt_bot_tap_lev= imdi ! Bottom level of spt calculations.
                                 ! Including tapering zone
INTEGER :: spt_botlev= imdi      ! Bottom level of spt calculations (min 2)
                                 ! Excluding tapering zone
INTEGER :: nsmooth_spt = imdi    ! Iteration count for spatial smoothing
                                 ! of SPT increments
INTEGER :: offx_spt, offy_spt    ! Halo size used in spatial smoothing
                                 ! for SPT

 !======================================================================
 ! Real values set from  RUN_Stochastic namelist via GUI
 !======================================================================
REAL :: tau_skeb = rmdi          ! decorrelation time for SKEB AR1 (in secs)
REAL :: tau_spt  = rmdi          ! decorrelation time for SPT AR1 (in secs)
REAL :: tot_backscat = rmdi      ! global-mean rate of energy backscatter
                                 ! in m**2 s**(-3) (usually 1.0e-4)
REAL :: br = rmdi                ! backscatter ratio (frac of assumed
                                 ! dissipated energy - usually 0.2)
REAL :: sdispfac = rmdi          ! Multiplication factor for numerical diss
                                 ! (determined empirically as 2.0)
REAL :: cdispfac = rmdi          ! Multiplication factor for convection diss
                                 ! (determined empirically as 1.0)
REAL :: kdispfac = rmdi          ! Multiplication factor for SKEB1 (KE) diss
                                 ! (determined empirically as 0.5)
REAL :: alphac = rmdi            ! Updraught proportion of gridbox
                                 ! (0.2% typical)

REAL :: rp2_decorr_ts = rmdi      ! RP2 de-correlation timescale (in seconds)

! RHCrit minimum and maximum default values (varies with resolution)
REAL :: rhcrit_max          = rmdi
REAL :: rhcrit_min          = rmdi

!------------------------------------

! Stochastic physics ice fallspeed multiplier m_ci,
! with minimum and maximum default values:
REAL :: m_ci                = rmdi
REAL :: m_ci_max            = rmdi
REAL :: m_ci_min            = rmdi
! Default, max and min value for rain particle size distribution, x1r
REAL :: x1r_rp              = rmdi ! suggested value = 0.22
REAL :: x1r_rp_max          = rmdi ! suggested value = 0.52
REAL :: x1r_rp_min          = rmdi ! suggested value = 0.07
! Default, max and min value for surface droplet number, ndrop_surf
REAL :: ndrop_surf_rp       = rmdi ! suggested value = 7.5e+07 
REAL :: ndrop_surf_rp_max   = rmdi ! suggested value = 10.0e+07
REAL :: ndrop_surf_rp_min   = rmdi ! suggested value = 2.0e+07
! Default, max and min value for autoconversion of cloud water to rain, ec_auto
REAL :: ec_auto_rp          = rmdi ! suggested value = 0.55
REAL :: ec_auto_rp_max      = rmdi ! suggested value = 0.6
REAL :: ec_auto_rp_min      = rmdi ! suggested value = 0.01

!------------------------------------

! ROSE - BOUNDARY LAYER stochastic items moved here from bl_option_mod
! Maximum and minimum values for the STPH_RP scheme
! Boundary Layer
! Max, mean and min value for the neutral mixing length
REAL :: par_mezcla_max= rmdi     ! suggested value = 0.5
REAL :: par_mezcla    = rmdi     ! suggested value = 0.15
REAL :: par_mezcla_min=rmdi      ! suggested value = 0.05
! Max,mean and min values for the flux profile parameter
REAL :: g0_rp_max = rmdi         ! suggested value = 20.0
REAL :: g0_rp  = rmdi            ! suggested value = 10.0
REAL :: g0_rp_min = rmdi         ! suggested value = 5.0
! Max and min values for the charnock parameter
REAL :: charnock_max= rmdi       ! suggested value = 0.026
REAL :: charnock_min= rmdi       ! suggested value = 0.01
! Max, mean and min values for the minimum mixing length
REAL :: lambda_min_rp_max= rmdi  ! suggested value = 100.0
REAL :: lambda_min_rp = rmdi     ! suggested value = 40.0
REAL :: lambda_min_rp_min= rmdi  ! suggested value = 20.0
! Max, mean and min values for the critical Ri
REAL :: ricrit_rp_max= rmdi      ! suggested value = 1.0
REAL :: ricrit_rp = rmdi         ! suggested value = 1.0
REAL :: ricrit_rp_min= rmdi      ! suggested value = 0.25
! Max, mean and min values for the entrainment parameter A1
REAL :: a_ent_1_rp_max = rmdi    ! suggested value = 0.4
REAL :: a_ent_1_rp  = rmdi       ! suggested value = 0.23
REAL :: a_ent_1_rp_min = rmdi    ! suggested value = 0.1
! Max and mean values for the entrainment parameter A1_shr
REAL :: a_ent_shr_rp_max = rmdi  ! suggested value = 5.0
REAL :: a_ent_shr_rp  = rmdi     ! suggested value = 1.6
! Max, mean and min values for the velocity scale parameter
REAL :: g1_rp_max= rmdi          ! suggested value = 1.5
REAL :: g1_rp = rmdi             ! suggested value = 0.85
REAL :: g1_rp_min= rmdi          ! suggested value = 0.5

! New parameter lam_meta to replace par_mezcla and lambda_min 
! (for use with i_rp_scheme == i_rp2b)
! Default, max and min values for lam_meta 
REAL :: lam_meta_rp        = rmdi ! suggested value = 1.0 
REAL :: lam_meta_rp_max    = rmdi ! suggested value = 3.0 
REAL :: lam_meta_rp_min    = rmdi ! suggested value = 0.2 

! Land-surface parameters for use with i_rp_scheme == i_rp2b
! 
! These are arrays corresponding to plant function types.
!
! Default, max and min values for JULES parameter dz0v_dh
REAL :: dz0v_dh_rp(npft_max)
REAL :: dz0v_dh_rp_max(npft_max)
REAL :: dz0v_dh_rp_min(npft_max)
DATA dz0v_dh_rp / npft_max * rmdi /
DATA dz0v_dh_rp_max / npft_max * rmdi /
DATA dz0v_dh_rp_min / npft_max * rmdi /
! Default, max and min values for JULES parameter z0hm_pft
REAL :: z0hm_pft_rp(npft_max)
REAL :: z0hm_pft_rp_max(npft_max)
REAL :: z0hm_pft_rp_min(npft_max)
DATA z0hm_pft_rp / npft_max * rmdi /
DATA z0hm_pft_rp_max / npft_max * rmdi /
DATA z0hm_pft_rp_min / npft_max * rmdi /
! LAI multiplier - used to scale lai
! Default, max and min values for lai_mult_rp
! Note lai_mult_rp default value is always 1.0
! so is initialised here and not read in from the namelist
REAL :: lai_mult_rp(npft_max)
REAL :: lai_mult_rp_max(npft_max)
REAL :: lai_mult_rp_min(npft_max)
DATA lai_mult_rp / npft_max * 1.0 /
DATA lai_mult_rp_max / npft_max * rmdi /
DATA lai_mult_rp_min / npft_max * rmdi /

!==============================================
! Stochastic forcing of theta in the BL options
!==============================================
! Switch to perturb theta: 0 = Off
!                          1 = with mag_pert_theta
!                          2 = using heat flux
!                          3 = theta and moisture
INTEGER :: i_pert_theta = imdi
INTEGER, PARAMETER :: pert_theta_mag = 1
INTEGER, PARAMETER :: pert_theta_star = 2
INTEGER, PARAMETER :: pert_theta_and_moist = 3

! Switch to perturb all points
! If false, perturb only cumulus points
LOGICAL :: l_pert_all_points = .FALSE.

! Switch to add vertical shape to perturbation profile if true
LOGICAL :: l_pert_shape = .FALSE.

! Switch to specify time variation of perturbations
! 0 = random sequence
! 1 = time correlated
INTEGER :: i_pert_theta_type = imdi
INTEGER, PARAMETER :: pert_theta_random_seq = 0                     
INTEGER, PARAMETER :: pert_theta_correl_seq = 1                  

! Maximum perturbation to theta (K)
REAL :: mag_pert_theta = rmdi

! Decorrelation timescale for theta perturbation (seconds)
REAL :: decorr_ts_pert_theta = rmdi

! Heights between which to pertub theta, converted to minlev_pert_theta 
! and maxlev_pert_theta in setcona
REAL :: zmin_pert_theta = rmdi
REAL :: zmax_pert_theta = rmdi

! Model levels between which to perturb theta, calculated in setcona
INTEGER :: minlev_pert_theta = imdi
INTEGER :: maxlev_pert_theta = imdi

! Number of points to apply the same perturbation over
INTEGER :: npts_pert_theta = imdi

  !======================================================================
  ! Arrays from Convection required by SKEB2
  !======================================================================
REAL, ALLOCATABLE ::                                                    &
      skeb2_up_flux(:, :, :)                                            &
           !  updraught mass flux
,     skeb2_dwn_flux(:, :, :)                                           &
           ! downdraught mass flux
,     skeb2_cape(:, :)                                                  &
           ! CAPE
,     Ymn(:,:,:)
           ! Matrix for the value of Legendre Polynomials


!========================================================================
! Standard deviation of the forcing pattern for the different
! parametrizations included in SPT.
!========================================================================

REAL ::                                                                 &

   rain_std= rmdi                                                       &
           ! Sigma for rain ECMWF 0.72
,  rad_std= rmdi                                                        &
           ! Sigma for radiation ECMWF 0.33
,  gwd_std= rmdi                                                        &
           ! Sigma for Gravity wave drag ECMWF 0.52
,  conv_std= rmdi                                                       &
           ! Sigma for convection ECMWF 0.62
,  sd_orog_thres = rmdi                                                 &
           ! Threshold of orographic sd to apply the capping
,  psif_orog_thres = rmdi                                                
           ! Threshold of psif above the one SPT pert are set to 0 over
           ! regions of high sd orog (specified by psif_orog_thres)
         
  !======================================================================
  ! Define Namelist RUN_Stochastic
  !======================================================================
NAMELIST/run_stochastic/                                                &
      l_skeb2, l_rp2, rp2_callfreq, i_rp_scheme, rp2_decorr_ts          &
,     l_rp2_cycle_out, l_rp2_cycle_in, rp2_cycle_tm                     &
,     stph_n1, stph_n2, br, tot_backscat, tau_skeb, alphac              &
,     l_skeb2_psicdisp, l_skeb2_psisdisp, l_skeb2_skeb1disp             &
,     sdispfac, cdispfac, kdispfac, skeb2_sdisp, skeb2_cdisp, nsmooth   &
,     skeb2_toplev, skeb2_botlev, l_skeb2_velpot, ran_max               &
,     rhcrit_ref_level, l_skebsmooth_adv, l_skebprint, stphseed         &
,     l_x_eq_sin_x, m_ci, m_ci_max, m_ci_min, rhcrit_max, rhcrit_min    &
,     x1r_rp_max,x1r_rp,x1r_rp_min                                      &
,     ndrop_surf_rp_max,ndrop_surf_rp,ndrop_surf_rp_min                 &
,     ec_auto_rp_max,ec_auto_rp,ec_auto_rp_min                          &
,     par_mezcla_max, par_mezcla, par_mezcla_min                        &
,     g0_rp_max, g0_rp, g0_rp_min, charnock_max, charnock_min           &
,     lambda_min_rp_max, lambda_min_rp, lambda_min_rp_min               &
,     lam_meta_rp, lam_meta_rp_min, lam_meta_rp_max                     &
,     ricrit_rp_max, ricrit_rp, ricrit_rp_min                           &
,     a_ent_1_rp_max, a_ent_1_rp, a_ent_1_rp_min                        &
,     a_ent_shr_rp, a_ent_shr_rp_max                                    &
,     g1_rp_max, g1_rp, g1_rp_min                                       &
,     lai_mult_rp_max, lai_mult_rp_min                                  &
,     dz0v_dh_rp_max, dz0v_dh_rp, dz0v_dh_rp_min                        &
,     z0hm_pft_rp_max, z0hm_pft_rp, z0hm_pft_rp_min                     &
,     l_skeb2_conv_disp_mod, l_pert_all_points, l_pert_shape            &
,     decorr_ts_pert_theta, i_pert_theta_type                           &
,     i_pert_theta, mag_pert_theta, zmin_pert_theta, zmax_pert_theta    &
,     npts_pert_theta                                                   &
! SPT namelist
,     l_spt, tau_spt, l_spt_rain, l_spt_rad, l_spt_gwd, l_spt_conv      &
,     l_spt_conv_mom, l_spt_cfl, spt_toplev, spt_botlev                 &
,     spt_bot_tap_lev, spt_top_tap_lev, nsmooth_spt                     &
,     rain_std, rad_std, gwd_std, conv_std,l_spt_qcons, l_spt_mse       &
,     sd_orog_thres, psif_orog_thres
 
  !======================================================================
  ! Array for spatial 1-2-1 smoothing
  !======================================================================
REAL, ALLOCATABLE, SAVE ::                                              &
      mask_pdamp(:, :)                                                  &
           !  Array of pattern damping coefficients for SKEB2
,     mask_smooth(:, :)                                                 &
           !  Array of smoothing coefficients for SKEB2
,     mask_smooth_spt(:, :)
           !  Array of smoothing coefficients for SPT
!======================================================================
! LOGICAL switches not set in namelist
!======================================================================


INTEGER  :: stphseed_unit = imdi
             ! Random seed file unit

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STOCHASTIC_PHYSICS_RUN_MOD'

CONTAINS

SUBROUTINE check_run_stochastic()

! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the run_stochastic namelist.

! Dr Hook Modules
USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY: umPrint, umMessage
IMPLICIT NONE

INTEGER :: icode    ! error code for ereport
! local temporary arrays
CHARACTER(LEN=errormessagelength)       :: cmessage      ! out error message
CHARACTER(LEN=*), PARAMETER  :: RoutineName='CHECK_RUN_STOCHASTIC'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! --------- Check that lai_mult_rp_min is not equal to zero in RP scheme
IF (l_rp2 .AND. i_rp_scheme == i_rp2b) THEN
  IF ( l_triffid .OR. l_phenol ) THEN 
      icode   = 1
      WRITE (cmessage,'(A)')                                            &
        '  **********************************  '//             newline//&
        'The RP2b scheme is on with either l_triffid'//        newline//&
        'or l_phenol set to true.  The RP2b scheme should'//   newline//&
        'not be used with either of these schemes.This is'//   newline//&
        'because LAI is a prognostic in the triffid and'//     newline//&
        'phenol schemes, and the RP2b scheme perturbs LAI'//   newline//&
        'with a method that assumes it is not a prognostic.'
      CALL ereport(routinename, icode, cmessage)
  END IF

  IF ( ANY (ABS(lai_mult_rp_min(:)) < TINY(0.0)) ) THEN 
      icode   = 1
      WRITE (cmessage,'(A)')                                            &
        '  **********************************  '//             newline//&
        'RP2b scheme is called with lai_mult_rp_min set'//     newline//&
        'to zero.  This may result in the model failing.'
      CALL ereport(routinename, icode, cmessage)
  END IF
END IF
! --------- Check that appropriate convective dissipation scheme is active
IF (l_skeb2_psicdisp) THEN
  IF (skeb2_cdisp /= type_mflx .AND. skeb2_cdisp /= type_cape) THEN
    WRITE(umMessage,'(A)')'  **********************************  '
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    WRITE(umMessage,'(A)')'SKEB Conv dissipation is on, but neither'
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    WRITE(umMessage,'(A)')'CAPE or mass-flux scheme is selected.'
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    icode   = 1
    WRITE (cmessage,'(A)')'SKEB2 CAPE and mass-flux convective ' &
          // 'dissipation are both off, but l_skeb2_psicdisp is true.'
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF
! --------- Check that appropriate numerical dissipation scheme is active
IF (l_skeb2_psisdisp) THEN
  IF (skeb2_sdisp /= type_smag .AND. skeb2_sdisp /= type_bihm) THEN
    WRITE(umMessage,'(A)')'  **********************************  '
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    WRITE(umMessage,'(A)')'SKEB Numerical dissipation is on, but neither'
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    WRITE(umMessage,'(A)')'Smagorinsky or Biharmonic scheme is selected.'
    CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    icode   = 1
    WRITE (cmessage,'(A)')'SKEB2 Biharmonic and Smagorinsky numerical ' &
          // 'dissipation are both off, but l_skeb2_psisdisp is true.'
    CALL ereport(routinename, icode, cmessage)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_stochastic

SUBROUTINE print_nlist_run_stochastic()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_STOCHASTIC'
CHARACTER(LEN=30) :: fmt_lsfc ! Format character for writing out
                              ! land-surface RP namelist values

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set write statement for land-surface random parameters based
! on the value of npft_max
WRITE(fmt_lsfc,'("(A,",I0, "ES12.4)")') npft_max

CALL umPrint('Contents of namelist run_stochastic',                     &
    src='stochastic_physics_run_mod')

WRITE(umMessage,'(A,L1)') 'l_skeb2 = ',l_skeb2
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_rp2 = ',l_rp2
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'rp2_callfreq = ',rp2_callfreq
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I10)') 'i_rp_scheme = ',i_rp_scheme
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'rp2_decorr_ts = ',rp2_decorr_ts
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_rp2_cycle_out = ',l_rp2_cycle_out
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_rp2_cycle_in = ',l_rp2_cycle_in
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I10)') 'rp2_cycle_tm = ',rp2_cycle_tm
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'stph_n1 = ',stph_n1
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'stph_n2 = ',stph_n2
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'br = ',br
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'tot_backscat = ',tot_backscat
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'tau_skeb = ',tau_skeb
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'tau_spt = ',tau_spt
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'alphac = ',alphac
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skeb2_psicdisp = ',l_skeb2_psicdisp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skeb2_psisdisp = ',l_skeb2_psisdisp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skeb2_skeb1disp = ',l_skeb2_skeb1disp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'sdispfac = ',sdispfac
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'cdispfac = ',cdispfac
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'kdispfac = ',kdispfac
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'skeb2_sdisp = ',skeb2_sdisp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'skeb2_cdisp = ',skeb2_cdisp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'nsmooth = ',nsmooth
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'skeb2_toplev = ',skeb2_toplev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'skeb2_botlev = ',skeb2_botlev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skeb2_velpot = ',l_skeb2_velpot
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'ran_max = ',ran_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'rhcrit_ref_level = ',rhcrit_ref_level
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skebsmooth_adv = ',l_skebsmooth_adv
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skebprint = ',l_skebprint
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_x_eq_sin_x = ',l_x_eq_sin_x
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'stphseed = ',stphseed
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'm_ci = ',m_ci
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'm_ci_max = ',m_ci_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'm_ci_min = ',m_ci_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'x1r_rp = ',x1r_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'x1r_rp_max = ',x1r_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'x1r_rp_min = ',x1r_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ndrop_surf_rp = ',ndrop_surf_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ndrop_surf_rp_max = ',ndrop_surf_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ndrop_surf_rp_min = ',ndrop_surf_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ec_auto_rp = ',ec_auto_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ec_auto_rp_max = ',ec_auto_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ec_auto_rp_min = ',ec_auto_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'rhcrit_max = ',rhcrit_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'rhcrit_min = ',rhcrit_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'par_mezcla_max = ',par_mezcla_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'par_mezcla = ',par_mezcla
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'par_mezcla_min = ',par_mezcla_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g0_rp_max = ',g0_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g0_rp = ',g0_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g0_rp_min = ',g0_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'charnock_max = ',charnock_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'charnock_min = ',charnock_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lambda_min_rp_max = ',lambda_min_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lambda_min_rp = ',lambda_min_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lambda_min_rp_min = ',lambda_min_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ricrit_rp_max = ',ricrit_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ricrit_rp = ',ricrit_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'ricrit_rp_min = ',ricrit_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'a_ent_1_rp_max = ',a_ent_1_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'a_ent_1_rp = ',a_ent_1_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'a_ent_1_rp_min = ',a_ent_1_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'a_ent_shr_rp_max = ',a_ent_shr_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'a_ent_shr_rp = ',a_ent_shr_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g1_rp_max = ',g1_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g1_rp = ',g1_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'g1_rp_min = ',g1_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lam_meta_rp_max = ',lam_meta_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lam_meta_rp = ',lam_meta_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'lam_meta_rp_min = ',lam_meta_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'lai_mult_rp_max = ',lai_mult_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'lai_mult_rp_min = ',lai_mult_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'dz0v_dh_rp_max = ',dz0v_dh_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'dz0v_dh_rp = ',dz0v_dh_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'dz0v_dh_rp_min = ',dz0v_dh_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'z0hm_pft_rp_max = ',z0hm_pft_rp_max
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'z0hm_pft_rp = ',z0hm_pft_rp
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,fmt_lsfc) 'z0hm_pft_rp_min = ',z0hm_pft_rp_min
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_skeb2_conv_disp_mod = ',l_skeb2_conv_disp_mod
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'i_pert_theta = ',i_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'i_pert_theta_type = ',i_pert_theta_type
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'decorr_ts_pert_theta = ',decorr_ts_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'mag_pert_theta = ',mag_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'zmin_pert_theta = ',zmin_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'zmax_pert_theta = ',zmax_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'npts_pert_theta = ',npts_pert_theta
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_pert_all_points = ',l_pert_all_points
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_pert_shape = ',l_pert_shape
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt = ',l_spt
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_rain = ',l_spt_rain
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_rad = ',l_spt_rad
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_gwd = ',l_spt_gwd
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_conv = ',l_spt_conv
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_conv_mom = ',l_spt_conv_mom
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_cfl = ',l_spt_cfl
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_qcons = ',l_spt_qcons
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,L1)') 'l_spt_mse = ',l_spt_mse
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'spt_toplev = ',spt_toplev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'spt_botlev = ',spt_botlev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'spt_bot_tap_lev = ',spt_bot_tap_lev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'spt_top_tap_lev = ',spt_top_tap_lev
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,I0)') 'nsmooth_spt  = ',nsmooth_spt
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'rain_std  = ',rain_std
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'rad_std  = ',rad_std
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'gwd_std  = ',gwd_std
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'conv_std  = ',conv_std
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'sd_orog_thres = ',sd_orog_thres
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
WRITE(umMessage,'(A,ES12.4)') 'psif_orog_thres = ',psif_orog_thres
CALL umPrint(umMessage,src='stochastic_physics_run_mod')
    
CALL umPrint('- - - - - - end of namelist - - - - - -',                 &
    src='stochastic_physics_run_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_stochastic

#if !defined(LFRIC)
SUBROUTINE read_nml_run_stochastic(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

USE physics_tendencies_mod, ONLY:                                       &
    l_retain_slow_tendencies,                                           &   
    l_retain_rad_tendencies, l_retain_mic_tendencies,                   &
    l_retain_gwd_tendencies, l_retain_conv_tendencies,                  &
    l_retain_conv_all_tendencies, l_retain_conv_mom_tendencies

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_STOCHASTIC'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 21
INTEGER, PARAMETER :: n_real = 58 + 8*npft_max
INTEGER, PARAMETER :: n_log = 23

TYPE my_namelist
  SEQUENCE
  INTEGER :: stph_n1
  INTEGER :: stph_n2
  INTEGER :: rp2_callfreq
  INTEGER :: rp2_cycle_tm
  INTEGER :: i_rp_scheme
  INTEGER :: skeb2_cdisp
  INTEGER :: skeb2_sdisp
  INTEGER :: nsmooth
  INTEGER :: skeb2_toplev
  INTEGER :: skeb2_botlev
  INTEGER :: ran_max
  INTEGER :: rhcrit_ref_level
  INTEGER :: stphseed
  INTEGER :: spt_top_tap_lev
  INTEGER :: spt_toplev
  INTEGER :: spt_bot_tap_lev
  INTEGER :: spt_botlev
  INTEGER :: nsmooth_spt
  INTEGER :: i_pert_theta
  INTEGER :: i_pert_theta_type
  INTEGER :: npts_pert_theta
  REAL :: br
  REAL :: tot_backscat
  REAL :: tau_skeb
  REAL :: tau_spt
  REAL :: alphac
  REAL :: sdispfac
  REAL :: cdispfac
  REAL :: kdispfac
  REAL :: rp2_decorr_ts
  REAL :: m_ci
  REAL :: m_ci_max
  REAL :: m_ci_min
  REAL :: x1r_rp
  REAL :: x1r_rp_max
  REAL :: x1r_rp_min
  REAL :: ndrop_surf_rp
  REAL :: ndrop_surf_rp_max
  REAL :: ndrop_surf_rp_min
  REAL :: ec_auto_rp
  REAL :: ec_auto_rp_max
  REAL :: ec_auto_rp_min
  REAL :: rhcrit_max
  REAL :: rhcrit_min
  REAL :: par_mezcla_max
  REAL :: par_mezcla
  REAL :: par_mezcla_min
  REAL :: g0_rp_max
  REAL :: g0_rp
  REAL :: g0_rp_min
  REAL :: charnock_max
  REAL :: charnock_min
  REAL :: lambda_min_rp_max
  REAL :: lambda_min_rp
  REAL :: lambda_min_rp_min
  REAL :: ricrit_rp_max
  REAL :: ricrit_rp
  REAL :: ricrit_rp_min
  REAL :: a_ent_1_rp_max
  REAL :: a_ent_1_rp
  REAL :: a_ent_1_rp_min
  REAL :: a_ent_shr_rp_max
  REAL :: a_ent_shr_rp
  REAL :: g1_rp_max
  REAL :: g1_rp
  REAL :: g1_rp_min
  REAL :: lam_meta_rp_max
  REAL :: lam_meta_rp
  REAL :: lam_meta_rp_min
  REAL :: lai_mult_rp_max(npft_max)
  REAL :: lai_mult_rp_min(npft_max)
  REAL :: dz0v_dh_rp_max(npft_max)
  REAL :: dz0v_dh_rp(npft_max)
  REAL :: dz0v_dh_rp_min(npft_max)
  REAL :: z0hm_pft_rp_max(npft_max)
  REAL :: z0hm_pft_rp(npft_max)
  REAL :: z0hm_pft_rp_min(npft_max)
  REAL :: rain_std
  REAL :: rad_std
  REAL :: gwd_std
  REAL :: conv_std
  REAL :: sd_orog_thres
  REAL :: psif_orog_thres
  REAL :: mag_pert_theta
  REAL :: zmin_pert_theta
  REAL :: zmax_pert_theta
  REAL :: decorr_ts_pert_theta
  LOGICAL :: l_skeb2
  LOGICAL :: l_pert_all_points
  LOGICAL :: l_pert_shape
  LOGICAL :: l_rp2
  LOGICAL :: l_rp2_cycle_out
  LOGICAL :: l_rp2_cycle_in
  LOGICAL :: l_x_eq_sin_x
  LOGICAL :: l_skeb2_psicdisp
  LOGICAL :: l_skeb2_psisdisp
  LOGICAL :: l_skeb2_skeb1disp
  LOGICAL :: l_skeb2_velpot
  LOGICAL :: l_skebsmooth_adv
  LOGICAL :: l_skebprint
  LOGICAL :: l_skeb2_conv_disp_mod
  LOGICAL :: l_spt
  LOGICAL :: l_spt_cfl
  LOGICAL :: l_spt_rain
  LOGICAL :: l_spt_rad
  LOGICAL :: l_spt_gwd
  LOGICAL :: l_spt_conv
  LOGICAL :: l_spt_conv_mom
  LOGICAL :: l_spt_qcons
  LOGICAL :: l_spt_mse
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,          &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Stochastic, IOSTAT=ErrorStatus,            &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Stochastic", iomessage)

  my_nml % stph_n1           = stph_n1
  my_nml % stph_n2           = stph_n2
  my_nml % rp2_callfreq      = rp2_callfreq
  my_nml % rp2_cycle_tm      = rp2_cycle_tm
  my_nml % i_rp_scheme       = i_rp_scheme
  my_nml % skeb2_cdisp       = skeb2_cdisp
  my_nml % skeb2_sdisp       = skeb2_sdisp
  my_nml % nsmooth           = nsmooth
  my_nml % skeb2_toplev      = skeb2_toplev
  my_nml % skeb2_botlev      = skeb2_botlev
  my_nml % ran_max           = ran_max
  my_nml % rhcrit_ref_level  = rhcrit_ref_level
  my_nml % stphseed          = stphseed
  my_nml % spt_top_tap_lev   = spt_top_tap_lev
  my_nml % spt_toplev        = spt_toplev
  my_nml % spt_bot_tap_lev   = spt_bot_tap_lev
  my_nml % spt_botlev        = spt_botlev
  my_nml % nsmooth_spt       = nsmooth_spt
  my_nml % i_pert_theta      = i_pert_theta
  my_nml % i_pert_theta_type = i_pert_theta_type
  my_nml % npts_pert_theta   = npts_pert_theta
  ! end of integers
  my_nml % br                = br
  my_nml % tot_backscat      = tot_backscat
  my_nml % tau_skeb          = tau_skeb
  my_nml % tau_spt           = tau_spt
  my_nml % alphac            = alphac
  my_nml % sdispfac          = sdispfac
  my_nml % cdispfac          = cdispfac
  my_nml % kdispfac          = kdispfac
  my_nml % rp2_decorr_ts     = rp2_decorr_ts
  my_nml % m_ci              = m_ci
  my_nml % m_ci_max          = m_ci_max
  my_nml % m_ci_min          = m_ci_min
  my_nml % x1r_rp            = x1r_rp
  my_nml % x1r_rp_max        = x1r_rp_max
  my_nml % x1r_rp_min        = x1r_rp_min
  my_nml % ndrop_surf_rp     = ndrop_surf_rp
  my_nml % ndrop_surf_rp_max = ndrop_surf_rp_max
  my_nml % ndrop_surf_rp_min = ndrop_surf_rp_min
  my_nml % ec_auto_rp        = ec_auto_rp
  my_nml % ec_auto_rp_max    = ec_auto_rp_max
  my_nml % ec_auto_rp_min    = ec_auto_rp_min
  my_nml % rhcrit_max        = rhcrit_max
  my_nml % rhcrit_min        = rhcrit_min
  my_nml % par_mezcla_max    = par_mezcla_max
  my_nml % par_mezcla        = par_mezcla
  my_nml % par_mezcla_min    = par_mezcla_min
  my_nml % g0_rp_max         = g0_rp_max
  my_nml % g0_rp             = g0_rp
  my_nml % g0_rp_min         = g0_rp_min
  my_nml % charnock_max      = charnock_max
  my_nml % charnock_min      = charnock_min
  my_nml % lambda_min_rp_max = lambda_min_rp_max
  my_nml % lambda_min_rp     = lambda_min_rp
  my_nml % lambda_min_rp_min = lambda_min_rp_min
  my_nml % ricrit_rp_max     = ricrit_rp_max
  my_nml % ricrit_rp         = ricrit_rp
  my_nml % ricrit_rp_min     = ricrit_rp_min
  my_nml % a_ent_1_rp_max    = a_ent_1_rp_max
  my_nml % a_ent_1_rp        = a_ent_1_rp
  my_nml % a_ent_1_rp_min    = a_ent_1_rp_min
  my_nml % a_ent_shr_rp_max  = a_ent_shr_rp_max
  my_nml % a_ent_shr_rp      = a_ent_shr_rp
  my_nml % g1_rp_max         = g1_rp_max
  my_nml % g1_rp             = g1_rp
  my_nml % g1_rp_min         = g1_rp_min
  my_nml % lam_meta_rp_max   = lam_meta_rp_max
  my_nml % lam_meta_rp       = lam_meta_rp
  my_nml % lam_meta_rp_min   = lam_meta_rp_min
  my_nml % lai_mult_rp_max   = lai_mult_rp_max
  my_nml % lai_mult_rp_min   = lai_mult_rp_min
  my_nml % dz0v_dh_rp_max    = dz0v_dh_rp_max
  my_nml % dz0v_dh_rp        = dz0v_dh_rp
  my_nml % dz0v_dh_rp_min    = dz0v_dh_rp_min
  my_nml % z0hm_pft_rp_max   = z0hm_pft_rp_max
  my_nml % z0hm_pft_rp       = z0hm_pft_rp
  my_nml % z0hm_pft_rp_min   = z0hm_pft_rp_min
  my_nml % rain_std          = rain_std
  my_nml % rad_std           = rad_std
  my_nml % gwd_std           = gwd_std
  my_nml % conv_std          = conv_std
  my_nml % sd_orog_thres     = sd_orog_thres
  my_nml % psif_orog_thres   = psif_orog_thres
  my_nml % mag_pert_theta    = mag_pert_theta
  my_nml % zmin_pert_theta   = zmin_pert_theta
  my_nml % zmax_pert_theta   = zmax_pert_theta
  my_nml % decorr_ts_pert_theta  = decorr_ts_pert_theta
  ! end of reals
  my_nml % l_skeb2           = l_skeb2
  my_nml % l_pert_all_points = l_pert_all_points
  my_nml % l_pert_shape      = l_pert_shape
  my_nml % l_rp2             = l_rp2
  my_nml % l_rp2_cycle_out   = l_rp2_cycle_out
  my_nml % l_rp2_cycle_in    = l_rp2_cycle_in
  my_nml % l_skeb2_psicdisp  = l_skeb2_psicdisp
  my_nml % l_skeb2_psisdisp  = l_skeb2_psisdisp
  my_nml % l_skeb2_skeb1disp = l_skeb2_skeb1disp
  my_nml % l_skeb2_velpot    = l_skeb2_velpot
  my_nml % l_skebsmooth_adv  = l_skebsmooth_adv
  my_nml % l_skebprint       = l_skebprint
  my_nml % l_x_eq_sin_x      = l_x_eq_sin_x
  my_nml % l_skeb2_conv_disp_mod = l_skeb2_conv_disp_mod
  my_nml % l_spt             = l_spt
  my_nml % l_spt_cfl         = l_spt_cfl
  my_nml % l_spt_rain        = l_spt_rain
  my_nml % l_spt_rad         = l_spt_rad
  my_nml % l_spt_gwd         = l_spt_gwd
  my_nml % l_spt_conv        = l_spt_conv
  my_nml % l_spt_conv_mom    = l_spt_conv_mom
  my_nml % l_spt_qcons       = l_spt_qcons
  my_nml % l_spt_mse         = l_spt_mse

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  stph_n1               = my_nml % stph_n1
  stph_n2               = my_nml % stph_n2
  rp2_callfreq          = my_nml % rp2_callfreq
  rp2_cycle_tm          = my_nml % rp2_cycle_tm
  i_rp_scheme           = my_nml % i_rp_scheme
  skeb2_cdisp           = my_nml % skeb2_cdisp
  skeb2_sdisp           = my_nml % skeb2_sdisp
  nsmooth               = my_nml % nsmooth
  skeb2_toplev          = my_nml % skeb2_toplev
  skeb2_botlev          = my_nml % skeb2_botlev
  ran_max               = my_nml % ran_max
  rhcrit_ref_level      = my_nml % rhcrit_ref_level
  stphseed              = my_nml % stphseed
  spt_top_tap_lev       = my_nml % spt_top_tap_lev
  spt_toplev            = my_nml % spt_toplev
  spt_bot_tap_lev       = my_nml % spt_bot_tap_lev
  spt_botlev            = my_nml % spt_botlev
  nsmooth_spt           = my_nml % nsmooth_spt
  i_pert_theta          = my_nml % i_pert_theta
  i_pert_theta_type     = my_nml % i_pert_theta_type
  npts_pert_theta       = my_nml % npts_pert_theta
  ! end of integers
  br                    = my_nml % br
  tot_backscat          = my_nml % tot_backscat
  tau_skeb              = my_nml % tau_skeb
  tau_spt               = my_nml % tau_spt
  alphac                = my_nml % alphac
  sdispfac              = my_nml % sdispfac
  cdispfac              = my_nml % cdispfac
  kdispfac              = my_nml % kdispfac
  rp2_decorr_ts         = my_nml % rp2_decorr_ts
  m_ci                  = my_nml % m_ci
  m_ci_max              = my_nml % m_ci_max
  m_ci_min              = my_nml % m_ci_min
  x1r_rp                = my_nml % x1r_rp
  x1r_rp_max            = my_nml % x1r_rp_max
  x1r_rp_min            = my_nml % x1r_rp_min
  ndrop_surf_rp         = my_nml % ndrop_surf_rp
  ndrop_surf_rp_max     = my_nml % ndrop_surf_rp_max
  ndrop_surf_rp_min     = my_nml % ndrop_surf_rp_min
  ec_auto_rp            = my_nml % ec_auto_rp
  ec_auto_rp_max        = my_nml % ec_auto_rp_max
  ec_auto_rp_min        = my_nml % ec_auto_rp_min
  rhcrit_max            = my_nml % rhcrit_max
  rhcrit_min            = my_nml % rhcrit_min
  par_mezcla_max        = my_nml % par_mezcla_max
  par_mezcla            = my_nml % par_mezcla
  par_mezcla_min        = my_nml % par_mezcla_min
  g0_rp_max             = my_nml % g0_rp_max
  g0_rp                 = my_nml % g0_rp
  g0_rp_min             = my_nml % g0_rp_min
  charnock_max          = my_nml % charnock_max
  charnock_min          = my_nml % charnock_min
  lambda_min_rp_max     = my_nml % lambda_min_rp_max
  lambda_min_rp         = my_nml % lambda_min_rp
  lambda_min_rp_min     = my_nml % lambda_min_rp_min
  ricrit_rp_max         = my_nml % ricrit_rp_max
  ricrit_rp             = my_nml % ricrit_rp
  ricrit_rp_min         = my_nml % ricrit_rp_min
  a_ent_1_rp_max        = my_nml % a_ent_1_rp_max
  a_ent_1_rp            = my_nml % a_ent_1_rp
  a_ent_1_rp_min        = my_nml % a_ent_1_rp_min
  a_ent_shr_rp_max      = my_nml % a_ent_shr_rp_max
  a_ent_shr_rp          = my_nml % a_ent_shr_rp
  g1_rp_max             = my_nml % g1_rp_max
  g1_rp                 = my_nml % g1_rp
  g1_rp_min             = my_nml % g1_rp_min
  lam_meta_rp_max       = my_nml % lam_meta_rp_max
  lam_meta_rp           = my_nml % lam_meta_rp
  lam_meta_rp_min       = my_nml % lam_meta_rp_min
  lai_mult_rp_max       = my_nml % lai_mult_rp_max
  lai_mult_rp_min       = my_nml % lai_mult_rp_min
  dz0v_dh_rp_max        = my_nml % dz0v_dh_rp_max
  dz0v_dh_rp            = my_nml % dz0v_dh_rp
  dz0v_dh_rp_min        = my_nml % dz0v_dh_rp_min
  z0hm_pft_rp_max       = my_nml % z0hm_pft_rp_max
  z0hm_pft_rp           = my_nml % z0hm_pft_rp
  z0hm_pft_rp_min       = my_nml % z0hm_pft_rp_min
  rain_std              = my_nml % rain_std
  rad_std               = my_nml % rad_std
  gwd_std               = my_nml % gwd_std
  conv_std              = my_nml % conv_std
  sd_orog_thres         = my_nml % sd_orog_thres   
  psif_orog_thres       = my_nml % psif_orog_thres 
  mag_pert_theta        = my_nml % mag_pert_theta
  zmin_pert_theta       = my_nml % zmin_pert_theta
  zmax_pert_theta       = my_nml % zmax_pert_theta
  decorr_ts_pert_theta  = my_nml % decorr_ts_pert_theta
  ! end of reals
  l_skeb2               = my_nml % l_skeb2
  l_pert_all_points     = my_nml % l_pert_all_points
  l_pert_shape          = my_nml % l_pert_shape
  l_rp2                 = my_nml % l_rp2
  l_rp2_cycle_out       = my_nml % l_rp2_cycle_out
  l_rp2_cycle_in        = my_nml % l_rp2_cycle_in
  l_skeb2_psicdisp      = my_nml % l_skeb2_psicdisp
  l_skeb2_psisdisp      = my_nml % l_skeb2_psisdisp
  l_skeb2_skeb1disp     = my_nml % l_skeb2_skeb1disp
  l_skeb2_velpot        = my_nml % l_skeb2_velpot
  l_skebsmooth_adv      = my_nml % l_skebsmooth_adv
  l_skebprint           = my_nml % l_skebprint
  l_x_eq_sin_x          = my_nml % l_x_eq_sin_x 
  l_skeb2_conv_disp_mod = my_nml % l_skeb2_conv_disp_mod
  l_spt                 = my_nml % l_spt
  l_spt_cfl             = my_nml % l_spt_cfl
  l_spt_rain            = my_nml % l_spt_rain
  l_spt_rad             = my_nml % l_spt_rad
  l_spt_gwd             = my_nml % l_spt_gwd
  l_spt_conv            = my_nml % l_spt_conv
  l_spt_conv_mom        = my_nml % l_spt_conv_mom
  l_spt_qcons           = my_nml % l_spt_qcons
  l_spt_mse             = my_nml % l_spt_mse 
            
END IF

CALL mpl_type_free(mpl_nml_type,icode)

! Set up retain tendencies flags
IF (l_spt)          l_retain_slow_tendencies     = .TRUE.
IF (l_spt)          l_retain_conv_all_tendencies = .TRUE.
IF (l_spt_rad)      l_retain_rad_tendencies      = .TRUE.
IF (l_spt_rain)     l_retain_mic_tendencies      = .TRUE.
IF (l_spt_gwd)      l_retain_gwd_tendencies      = .TRUE.
IF (l_spt_conv)     l_retain_conv_tendencies     = .TRUE.
IF (l_spt_conv_mom) l_retain_conv_mom_tendencies = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_stochastic
#endif

END MODULE stochastic_physics_run_mod
