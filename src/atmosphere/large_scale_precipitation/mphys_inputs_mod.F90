! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!  data module for switches/options concerned with the precipitation scheme.
  ! Description:
  !   Module containing runtime options/data used by the precipitation scheme

  ! Method:
  !   Switches and associated data values used by the precipitation scheme
  !   are defined here and assiged default values. These may be overridden
  !   by namelist input.

  !   A description of what each switch or number refers to is provided
  !   with the namelist

  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Large Scale Precipitation

  ! Code Description:
  !   Language: FORTRAN 90


MODULE mphys_inputs_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!===========================================================================
! INTEGER options set from RUN_PRECIP namelist
!===========================================================================

! Method for iterating microphysics scheme
INTEGER :: i_mcr_iter       = imdi
! Values to which this can be set
! No iterating the microphysics
INTEGER, PARAMETER :: i_mcr_iter_none   = 0
! Iterating on: user sets number of iterations per model timestep
INTEGER, PARAMETER :: i_mcr_iter_niters = 1
! Iterating on: user set microphysics "timestep" and number of
!               iterations is calculated
INTEGER, PARAMETER :: i_mcr_iter_tstep  = 2
!------------------------------------

! Number of iterations of microphysics scheme
!  Input parameter for i_mcr_iter=i_mcr_iter_niters
!  Calculated at run-time for i_mcr_iter=i_mcr_iter_tstep
INTEGER :: niters_mp        = imdi
!------------------------------------

! Requested (input) length of microphysics timestep (s)
!  Input parameter for i_mcr_iter=i_mcr_iter_tstep
INTEGER :: timestep_mp_in    = imdi

! Number of bins in the mixed phase calculation
INTEGER :: nbins_mp         = imdi
!------------------------------------

! Choice for each microphysical species for CASIM
INTEGER :: casim_moments_choice = imdi

!===========================================================================
! LOGICAL options set from RUN_PRECIP namelist
!===========================================================================

! Use improved warm rain microphysics scheme
LOGICAL :: l_warm_new       = .FALSE.
!--------------------------------------

! Use generic
! ice psd
LOGICAL :: l_psd            = .FALSE.
!------------------------------------

! Use global version
! (mid-lat selected if
! .false.)
LOGICAL :: l_psd_global     = .FALSE.
!------------------------------------

! Use murk aerosol to
! calculate the droplet
! number
LOGICAL :: l_autoconv_murk  = .FALSE.
!------------------------------------

! Enable tapering of cloud droplets
! towards surface
LOGICAL :: l_droplet_tpr    = .FALSE.
!------------------------------------

! Use the Clark et al (2008) aerosol
! scheme in MURK calculations
LOGICAL :: l_clark_aero     = .FALSE.

! New variant of taper curve
! with variable aerosol at surface
LOGICAL :: l_taper_new      = .FALSE.
!------------------------------------

! Use Abel & Shipway
! rain fall speeds
LOGICAL :: l_rainfall_as    = .FALSE.
!------------------------------------

! Allow snow-rain collisions to produce
! graupel
LOGICAL :: l_sr2graup       = .FALSE.
!------------------------------------

! Use Aerosol climatologies to generate drop number
LOGICAL :: l_mcr_arcl       = .FALSE.
!-----------------------------------

! Include prognostic rain
LOGICAL :: l_mcr_qrain      = .FALSE.
!-----------------------------------

! Include prognosic graupel
LOGICAL :: l_mcr_qgraup     = .FALSE.
!-----------------------------------

! Prognostic rain lbcs active
LOGICAL :: l_mcr_qrain_lbc  = .FALSE.
!-----------------------------------

! Prognostic graupel lbcs active
LOGICAL :: l_mcr_qgraup_lbc = .FALSE.
!-----------------------------------

! Turns precipitation code on/off
LOGICAL :: l_rain           = .FALSE.
!-----------------------------------

! Turns seeder-feeder code on/off
LOGICAL :: l_orograin       = .FALSE.
!-----------------------------------

! In seeder feeder code enhance riming also
LOGICAL :: l_orogrime       = .FALSE.
!-----------------------------------

! Account for blocking in seeder feeder
LOGICAL :: l_orograin_block = .FALSE.
!-----------------------------------

! Use same lwc FSD for autoconv and accretion as in cloud generator
LOGICAL :: l_fsd_generator = .FALSE.

! Use different fallspeed relations for
! crystals and aggregates with the generic psd
LOGICAL :: l_diff_icevt = .FALSE.
!-----------------------------------

! Use sulphate aerosol in microphysics
LOGICAL :: l_use_sulphate_autoconv = .FALSE.

! Use sea-salt aerosol in microphysics
LOGICAL :: l_use_seasalt_autoconv = .FALSE.

! Use biomass aerosol in microphysics
LOGICAL :: l_use_bmass_autoconv = .FALSE.

! Use fossil-fuel organic carbon in microphysics
LOGICAL :: l_use_ocff_autoconv = .FALSE.

! Use ammonium nitrate aerosol in microphysics
LOGICAL :: l_use_nitrate_autoconv = .FALSE.

! Use autoconversion de-biasing scheme in microphysics
LOGICAL :: l_auto_debias = .FALSE.

! Produce extra qcl by turbulent processes
LOGICAL :: l_subgrid_qcl_mp = .FALSE.

! Produce cfliquid by erosion
LOGICAL :: l_subgrid_cfl_mp_by_erosion = .FALSE.

! Apply a temperature limit to mixed phase calculations
LOGICAL :: l_mixed_phase_t_limit = .FALSE.

! Apply the shape-dependent riming parametrization
LOGICAL :: l_shape_rime = .FALSE.

! CASIM Microphysics- if switch is True, run with CASIM, otherwise
! use Wilson and Ballard Microphysics
LOGICAL :: l_casim = .FALSE.

!===========================================================================
! REAL values set from RUN_PRECIP namelist
!===========================================================================

! Rain particle size distribution
! values
REAL :: x1r                 = rmdi
REAL :: x2r                 = rmdi
!------------------------------------

! Ice mass-diameter relationship
! values
REAL :: ai                  = rmdi
REAL :: bi                  = rmdi
!------------------------------------

!-----------------------------------------------------------------------------
! Fallspeed parameters for crystals and aggregates to use with generic psd to
! us when l_diff_fallspeed = .true.
!-----------------------------------------------------------------------------
!-- Linear ice vt ------------------------------------------------------------
! Fallspeed parameters for crystals
REAL :: cic_input = rmdi !1042.18
REAL :: dic_input = rmdi !1.0
! Fallspeed parameters for aggregates
REAL :: ci_input = rmdi !14.2611
REAL :: di_input = rmdi !0.416351
!-----------------------------------------------------------------------------

! Droplet taper height

REAL :: z_peak_nd           = rmdi
!------------------------------------

! Droplet number at z_surf and below:

REAL :: ndrop_surf          = rmdi

! Height at which droplet number reaches ndrop_surf
REAL :: z_surf = rmdi

!------------------------------------
! Maximum droplet number assumed at model level 1:

REAL :: max_drop_surf = rmdi
!------------------------------------

! Cloud-rain correlation coefficient for inhomogeneity parametrization:
! i.e. when l_inhomog=.true.

REAL :: c_r_correl = rmdi
!------------------------------------

! Aerosol climatology scaling factor, to account for spatial/temporal
! inhomogeneity
REAL :: arcl_inhom_sc  = rmdi
!------------------------------------

! Axial ratios for aggregates and crystals
REAL :: ar  = rmdi
REAL :: arc = rmdi
!-----------------------------------

! Maximum Temp for ice nuclei nucleation (deg C)
! Typically minus 10 C.
REAL :: tnuc = rmdi
!-----------------------------------

! Temperature limit for mixed phase
REAL :: mp_t_limit   = rmdi

! Scaling value for grid spacing
REAL :: mp_dz_scal   = rmdi

! Limit for phase relaxation timescale
REAL :: mp_tau_d_lim = rmdi

! Constant used in
REAL :: mp_czero = rmdi

! Shape dependent riming rate parameters
!
REAL :: qclrime = rmdi
REAL :: a_ratio_fac = rmdi
REAL :: a_ratio_exp = rmdi

! Critical Froude number for blocking in Seeder Feeder scheme
REAL :: fcrit = rmdi
! Scaling parameter for sub-grid hill peak-to-trough height
REAL :: nsigmasf = rmdi
! Scaling parameter for horizontal wavelength of sub-grid sinusoidal
! ridge for use in vertical decay function
REAL :: nscalesf = rmdi


!----------------------------------------------------------------------

! Define the RUN_PRECIP namelist

NAMELIST/run_precip/                                                          &
       l_warm_new, l_psd, l_psd_global, l_autoconv_murk,                      &
       l_use_sulphate_autoconv, l_use_seasalt_autoconv, l_use_bmass_autoconv, &
       l_use_ocff_autoconv, l_use_nitrate_autoconv,                           &
       x1r, x2r, c_r_correl, ai, bi, ar, arc, tnuc,                           &
       cic_input,dic_input,ci_input,di_input,l_diff_icevt,                    &
       z_peak_nd, ndrop_surf, z_surf, l_droplet_tpr, l_clark_aero,            &
       l_taper_new, max_drop_surf, l_rainfall_as,                             &
       i_mcr_iter, niters_mp, timestep_mp_in,                                 &
       l_sr2graup, l_mcr_arcl, arcl_inhom_sc,                                 &
       l_mcr_qrain, l_mcr_qgraup, l_mcr_qrain_lbc, l_mcr_qgraup_lbc,          &
       l_rain, l_fsd_generator, l_subgrid_qcl_mp,                             &
       l_subgrid_cfl_mp_by_erosion, l_mixed_phase_t_limit, nbins_mp,          &
       mp_t_limit, mp_dz_scal, mp_tau_d_lim, mp_czero,                        &
       l_shape_rime, qclrime, a_ratio_fac, a_ratio_exp,                       &
       l_orograin, l_orogrime, l_orograin_block,                              &
       fcrit, nsigmasf, nscalesf, l_casim, casim_moments_choice

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================
! Do not use generic ice particle size distribution in calculations
LOGICAL, PARAMETER :: not_generic_size_dist = .FALSE.

! Second ice variable lbcs active
LOGICAL :: l_mcr_qcf2_lbc   = .FALSE.
!-----------------------------------

! Include second ice variable
LOGICAL :: l_mcr_qcf2       = .FALSE.
!-----------------------------------

! Use soot aerosol in microphysics (NOT set in namelist)
LOGICAL :: l_use_soot_autoconv = .FALSE.
!-----------------------------------

!==========================================================================
! CASIM Options for future inclusion in namelist
!==========================================================================

! These variables are intended to form part of the run_precip namelist
! in the near future and hence are set and initialised here to values which
! are known to work with CASIM at present.

! Choice for the aerosol modes option for CASIM
INTEGER, PARAMETER :: casim_aerosol_option = 0

! Choice for the aerosol coupling method in CASIM
INTEGER, PARAMETER :: casim_aerosol_couple_choice = 0

! Choice for the aerosol processing level in CASIM
INTEGER, PARAMETER :: casim_aerosol_process_level = 0

! CASIM aerosol processing has separate rain category
LOGICAL, PARAMETER :: l_separate_process_rain = .FALSE.

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

LOGICAL :: l_no_cf = .FALSE.
! A temporary logical for an advanced user to turn off the effect of cloud
! fraction in the radar reflectivity code. This is not intended to be added
! to the namelist (where it may confuse the user) and is just for testing
! and checking of the radar reflectivity code.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_INPUTS_MOD'

CONTAINS

SUBROUTINE print_nlist_run_precip()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_PRECIP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_precip',                               &
    src='mphys_inputs_mod')

WRITE(lineBuffer,'(A,L1)')' l_warm_new = ',l_warm_new
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_psd = ',l_psd
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_psd_global = ',l_psd_global
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_autoconv_murk = ',l_autoconv_murk
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_use_sulphate_autoconv = ',l_use_sulphate_autoconv
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_use_seasalt_autoconv = ',l_use_seasalt_autoconv
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_use_ocff_autoconv = ',l_use_ocff_autoconv
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_use_nitrate_autoconv = ',l_use_nitrate_autoconv
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_use_bmass_autoconv = ',l_use_bmass_autoconv
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' x1r = ',x1r
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' x2r = ',x2r
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' c_r_correl = ',c_r_correl
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' ai = ',ai
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' bi = ',bi
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' ar = ',ar
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' arc = ',arc
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' tnuc = ',tnuc
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' z_peak_nd = ',z_peak_nd
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' ndrop_surf = ',ndrop_surf
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' z_surf = ',z_surf
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_droplet_tpr = ',l_droplet_tpr
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_clark_aero = ',l_clark_aero
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_taper_new = ',l_taper_new
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' max_drop_surf = ',max_drop_surf
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_rainfall_as = ',l_rainfall_as
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,I0)')' i_mcr_iter = ',i_mcr_iter
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,I0)')' niters_mp = ',niters_mp
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,I0)')' timestep_mp_in = ',timestep_mp_in
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_sr2graup = ',l_sr2graup
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mcr_arcl = ',l_mcr_arcl
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' arcl_inhom_sc = ',arcl_inhom_sc
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mcr_qrain =', l_mcr_qrain
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mcr_qgraup =', l_mcr_qgraup
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mcr_qrain_lbc =', l_mcr_qrain_lbc
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mcr_qgraup_lbc =', l_mcr_qgraup_lbc
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_rain =',l_rain
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_fsd_generator =',l_fsd_generator
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_subgrid_qcl_mp =',l_subgrid_qcl_mp
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_subgrid_cfl_mp_by_erosion =',  &
                            l_subgrid_cfl_mp_by_erosion
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_mixed_phase_t_limit =',l_mixed_phase_t_limit
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,I0)')' nbins_mp =', nbins_mp
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' mp_t_limit =', mp_t_limit
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' mp_dz_scal =', mp_dz_scal
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' mp_tau_d_lim =', mp_tau_d_lim
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' mp_czero =', mp_czero
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_shape_rime =', l_shape_rime
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' qclrime =', qclrime
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' a_ratio_fac =', a_ratio_fac
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.4)')' a_ratio_exp =', a_ratio_exp
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_orograin =', l_orograin
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_orogrime =', l_orogrime
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_orograin_block =', l_orograin_block
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)')' fcrit =', fcrit
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)')' nsigmasf =', nsigmasf
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)')' nscalesf =', nscalesf
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_casim =', l_casim
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
WRITE(lineBuffer,'(A,I0)')' casim_moments_choice =', casim_moments_choice
CALL umPrint(lineBuffer,src='mphys_inputs_mod')
CALL umPrint('- - - - - - end of namelist - - - - - -',                       &
    src='mphys_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_precip

#if !defined(LFRIC)
SUBROUTINE read_nml_run_precip(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_PRECIP'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 5
INTEGER, PARAMETER :: n_real = 27
INTEGER, PARAMETER :: n_log = 30

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_mcr_iter
  INTEGER :: niters_mp
  INTEGER :: timestep_mp_in
  INTEGER :: nbins_mp
  INTEGER :: casim_moments_choice
  REAL :: x1r
  REAL :: x2r
  REAL :: c_r_correl
  REAL :: ai
  REAL :: bi
  REAL :: ar
  REAL :: arc
  REAL :: tnuc
  REAL :: cic_input
  REAL :: dic_input
  REAL :: ci_input
  REAL :: di_input
  REAL :: z_peak_nd
  REAL :: ndrop_surf
  REAL :: z_surf
  REAL :: max_drop_surf
  REAL :: arcl_inhom_sc
  REAL :: mp_t_limit
  REAL :: mp_dz_scal
  REAL :: mp_tau_d_lim
  REAL :: mp_czero
  REAL :: qclrime
  REAL :: a_ratio_fac
  REAL :: a_ratio_exp
  REAL :: fcrit
  REAL :: nsigmasf
  REAL :: nscalesf
  LOGICAL :: l_warm_new
  LOGICAL :: l_psd
  LOGICAL :: l_psd_global
  LOGICAL :: l_autoconv_murk
  LOGICAL :: l_use_sulphate_autoconv
  LOGICAL :: l_use_seasalt_autoconv
  LOGICAL :: l_use_bmass_autoconv
  LOGICAL :: l_use_ocff_autoconv
  LOGICAL :: l_use_nitrate_autoconv
  LOGICAL :: l_diff_icevt
  LOGICAL :: l_droplet_tpr
  LOGICAL :: l_clark_aero
  LOGICAL :: l_taper_new
  LOGICAL :: l_rainfall_as
  LOGICAL :: l_sr2graup
  LOGICAL :: l_mcr_arcl
  LOGICAL :: l_mcr_qrain
  LOGICAL :: l_mcr_qgraup
  LOGICAL :: l_mcr_qrain_lbc
  LOGICAL :: l_mcr_qgraup_lbc
  LOGICAL :: l_rain
  LOGICAL :: l_fsd_generator
  LOGICAL :: l_subgrid_qcl_mp
  LOGICAL :: l_subgrid_cfl_mp_by_erosion
  LOGICAL :: l_mixed_phase_t_limit
  LOGICAL :: l_shape_rime
  LOGICAL :: l_orograin
  LOGICAL :: l_orogrime
  LOGICAL :: l_orograin_block
  LOGICAL :: l_casim
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,                &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Precip, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Precip", iomessage)

  my_nml % i_mcr_iter    = i_mcr_iter
  my_nml % niters_mp     = niters_mp
  my_nml % timestep_mp_in= timestep_mp_in
  my_nml % nbins_mp = nbins_mp
  my_nml % casim_moments_choice = casim_moments_choice
  ! end of integers
  my_nml % x1r           = x1r
  my_nml % x2r           = x2r
  my_nml % c_r_correl    = c_r_correl
  my_nml % ai            = ai
  my_nml % bi            = bi
  my_nml % ar            = ar
  my_nml % arc           = arc
  my_nml % tnuc          = tnuc
  my_nml % cic_input     = cic_input
  my_nml % dic_input     = dic_input
  my_nml % ci_input      = ci_input
  my_nml % di_input      = di_input
  my_nml % z_peak_nd     = z_peak_nd
  my_nml % ndrop_surf    = ndrop_surf
  my_nml % z_surf        = z_surf
  my_nml % max_drop_surf = max_drop_surf
  my_nml % arcl_inhom_sc = arcl_inhom_sc
  my_nml % mp_t_limit    = mp_t_limit
  my_nml % mp_dz_scal    = mp_dz_scal
  my_nml % mp_tau_d_lim  = mp_tau_d_lim
  my_nml % mp_czero      = mp_czero
  my_nml % qclrime       = qclrime
  my_nml % a_ratio_fac   = a_ratio_fac
  my_nml % a_ratio_exp   = a_ratio_exp
  my_nml % fcrit         = fcrit
  my_nml % nsigmasf      = nsigmasf
  my_nml % nscalesf      = nscalesf
  ! end of reals
  my_nml % l_warm_new       = l_warm_new
  my_nml % l_psd            = l_psd
  my_nml % l_psd_global     = l_psd_global
  my_nml % l_autoconv_murk  = l_autoconv_murk
  my_nml % l_use_sulphate_autoconv = l_use_sulphate_autoconv
  my_nml % l_use_seasalt_autoconv  = l_use_seasalt_autoconv
  my_nml % l_use_bmass_autoconv    = l_use_bmass_autoconv
  my_nml % l_use_ocff_autoconv     = l_use_ocff_autoconv
  my_nml % l_use_nitrate_autoconv  = l_use_nitrate_autoconv
  my_nml % l_diff_icevt      = l_diff_icevt
  my_nml % l_droplet_tpr     = l_droplet_tpr
  my_nml % l_clark_aero      = l_clark_aero
  my_nml % l_taper_new       = l_taper_new
  my_nml % l_rainfall_as     = l_rainfall_as
  my_nml % l_sr2graup        = l_sr2graup
  my_nml % l_mcr_arcl        = l_mcr_arcl
  my_nml % l_mcr_qrain       = l_mcr_qrain
  my_nml % l_mcr_qgraup      = l_mcr_qgraup
  my_nml % l_mcr_qrain_lbc   = l_mcr_qrain_lbc
  my_nml % l_mcr_qgraup_lbc  = l_mcr_qgraup_lbc
  my_nml % l_rain            = l_rain
  my_nml % l_fsd_generator   = l_fsd_generator
  my_nml % l_subgrid_qcl_mp  = l_subgrid_qcl_mp
  my_nml % l_subgrid_cfl_mp_by_erosion = l_subgrid_cfl_mp_by_erosion
  my_nml % l_mixed_phase_t_limit = l_mixed_phase_t_limit
  my_nml % l_shape_rime = l_shape_rime
  my_nml % l_orograin = l_orograin
  my_nml % l_orogrime = l_orogrime
  my_nml % l_orograin_block = l_orograin_block
  my_nml % l_casim      = l_casim

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_mcr_iter = my_nml % i_mcr_iter
  niters_mp  = my_nml % niters_mp
  timestep_mp_in = my_nml % timestep_mp_in
  nbins_mp   = my_nml % nbins_mp
  casim_moments_choice = my_nml % casim_moments_choice
  ! end of integers
  x1r        = my_nml % x1r
  x2r        = my_nml % x2r
  c_r_correl = my_nml % c_r_correl
  ai         = my_nml % ai
  bi         = my_nml % bi
  ar         = my_nml % ar
  arc        = my_nml % arc
  tnuc       = my_nml % tnuc
  cic_input  = my_nml % cic_input
  dic_input  = my_nml % dic_input
  ci_input   = my_nml % ci_input
  di_input   = my_nml % di_input
  z_peak_nd  = my_nml % z_peak_nd
  ndrop_surf = my_nml % ndrop_surf
  z_surf     = my_nml % z_surf
  max_drop_surf = my_nml % max_drop_surf
  arcl_inhom_sc = my_nml % arcl_inhom_sc
  mp_t_limit    = my_nml % mp_t_limit
  mp_dz_scal    = my_nml % mp_dz_scal
  mp_tau_d_lim  = my_nml % mp_tau_d_lim
  mp_czero      = my_nml % mp_czero
  qclrime       = my_nml % qclrime
  a_ratio_fac   = my_nml % a_ratio_fac
  a_ratio_exp   = my_nml % a_ratio_exp
  fcrit         = my_nml % fcrit
  nsigmasf      = my_nml % nsigmasf
  nscalesf      = my_nml % nscalesf
  ! end of  reals
  l_warm_new       = my_nml % l_warm_new
  l_psd            = my_nml % l_psd
  l_psd_global     = my_nml % l_psd_global
  l_autoconv_murk  = my_nml % l_autoconv_murk
  l_use_sulphate_autoconv = my_nml % l_use_sulphate_autoconv
  l_use_seasalt_autoconv  = my_nml % l_use_seasalt_autoconv
  l_use_bmass_autoconv    = my_nml % l_use_bmass_autoconv
  l_use_ocff_autoconv     = my_nml % l_use_ocff_autoconv
  l_use_nitrate_autoconv  = my_nml % l_use_nitrate_autoconv
  l_diff_icevt      = my_nml % l_diff_icevt
  l_droplet_tpr     = my_nml % l_droplet_tpr
  l_clark_aero      = my_nml % l_clark_aero
  l_taper_new       = my_nml % l_taper_new
  l_rainfall_as     = my_nml % l_rainfall_as
  l_sr2graup        = my_nml % l_sr2graup
  l_mcr_arcl        = my_nml % l_mcr_arcl
  l_mcr_qrain       = my_nml % l_mcr_qrain
  l_mcr_qgraup      = my_nml % l_mcr_qgraup
  l_mcr_qrain_lbc   = my_nml % l_mcr_qrain_lbc
  l_mcr_qgraup_lbc  = my_nml % l_mcr_qgraup_lbc
  l_rain            = my_nml % l_rain
  l_fsd_generator   = my_nml % l_fsd_generator
  l_subgrid_qcl_mp  = my_nml % l_subgrid_qcl_mp
  l_subgrid_cfl_mp_by_erosion = my_nml % l_subgrid_cfl_mp_by_erosion
  l_mixed_phase_t_limit = my_nml % l_mixed_phase_t_limit
  l_shape_rime      = my_nml % l_shape_rime
  l_orograin        = my_nml % l_orograin
  l_orogrime        = my_nml % l_orogrime
  l_orograin_block  = my_nml % l_orograin_block
  l_casim           = my_nml % l_casim

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_precip
#endif

SUBROUTINE check_run_precip

USE umprintmgr,      ONLY: newline
USE ereport_mod,     ONLY: ereport
USE chk_opts_mod,    ONLY: chk_var, def_src
USE murk_inputs_mod, ONLY: l_murk

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_PRECIP'

CHARACTER(LEN=errormessagelength) :: comments
CHARACTER(LEN=100) :: ChkStr

INTEGER :: ErrorStatus

REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

ErrorStatus = 0

! Aerosol indirect effect
IF ( l_autoconv_murk .AND. l_mcr_arcl ) THEN
  comments='Cannot set both l_autoconv_murk and l_mcr_arcl to .true.'
  ErrorStatus=1
  CALL ereport(RoutineName,ErrorStatus,comments)
END IF

IF (l_rain) THEN
  ! Only bother checking the namelist if the precip scheme is in use.

  IF (l_casim) THEN
    CALL chk_var( casim_moments_choice, 'casim_moments_choice', '[0, 1, 8]')

    IF (l_psd .OR. l_psd_global) THEN

      ! This does not work with CASIM and can call users
      ! jobs to hang without crashing if used together
      ! Best is to Call Ereport to prevent this happening.
      ErrorStatus = 100

      comments = 'Generic ice PSD does not work with CASIM.     '// newline// &
      'Please open the Rose Gui and switch off l_psd, ensuring  '// newline// &
      'that l_psd_global is also set to false. This model run   '// newline// &
      'crash to save you the annoyance of it hanging and later  '// newline// &
      'timing out without any obvious error'

      CALL ereport(RoutineName, ErrorStatus, comments)

    END IF

  ELSE ! l_casim

    ! Check options for the Wilson and Ballard Microphysics

    WRITE(ChkStr, '(3(A,I1),A)')                                              &
         '[',i_mcr_iter_none,',',i_mcr_iter_niters,',',i_mcr_iter_tstep,']'

    ! First check i_mcr_iter and fall over if it is not set correctly
    CALL chk_var(i_mcr_iter, 'i_mcr_iter', ChkStr)

    SELECT CASE (i_mcr_iter)

    ! N.B. Case of i_mcr_iter_none has nothing to check

    CASE (i_mcr_iter_niters)
      ! Check number of iterations matches metadata
      CALL chk_var(niters_mp, 'niters_mp', '[1:999]')

    CASE (i_mcr_iter_tstep)
      ! Check microphysics timestep matches metadata
      CALL chk_var(timestep_mp_in, 'timestep_mp_in', '[1:3600]')

    END SELECT

    ! The options below will always need checking in the run
    CALL chk_var(x1r, 'x1r', '[>=0.0]')
    CALL chk_var(x2r, 'x2r', '[-10.0:10.0]')
    CALL chk_var(ai,  'ai', '[>=0.0]')
    CALL chk_var(bi,  'bi', '[>=0.0]')
    CALL chk_var(ar,  'ar', '[0.01:10.0]')
    CALL chk_var(tnuc, 'tnuc', '[-100.0:0.0]')

    IF (l_psd .AND. l_diff_icevt) THEN

      ! With both options set, need to check the input parameters to the
      ! split fall speed (cic_input, dic_input, ci_input and di_input)

      CALL chk_var(cic_input, 'cic_input', '[0.0:1.0E7]')
      CALL chk_var(dic_input, 'dic_input', '[0.0:1.0E7]')
      CALL chk_var(ci_input,  'ci_input', '[0.0:1.0E7]')
      CALL chk_var(di_input,  'di_input', '[0.0:1.0E7]')

    ELSE IF (.NOT. l_psd) THEN
      ! Run is not using the Generic Ice PSD. Must check all crystal options
      ! are correctly set.

      CALL chk_var(arc, 'arc', '[0.01:10.0]')

    END IF ! not l_psd

    IF (l_droplet_tpr) THEN

      CALL chk_var(z_peak_nd, 'z_peak_nd',   '[20.0:5000.0]')
      CALL chk_var(ndrop_surf, 'ndrop_surf', '[1.0E6:375.0E6]')
      CALL chk_var(z_surf, 'z_surf', '[>=0.0]')

      ! Check that z_surf <= z_peak_nd by writing to ChkStr:
      WRITE(ChkStr,'(A,F8.4,A)') '[<=',z_peak_nd,']'

      ! Then call chk_var with the above ChkStr
      CALL chk_var(z_surf, 'z_surf', ChkStr)

      IF (l_taper_new) THEN
        CALL chk_var(max_drop_surf, 'max_drop_surf', '[1.0E7:20.0E7]')
      END IF

    END IF ! l_droplet_tpr

    IF (l_warm_new) CALL chk_var(c_r_correl, 'c_r_correl', '[-1.0:1.0]')

    IF (l_mcr_arcl) CALL chk_var(arcl_inhom_sc, 'arcl_inhom_sc', '[0.01:10.0]')

    IF (l_subgrid_qcl_mp) THEN

      IF (l_mixed_phase_t_limit) THEN
        CALL chk_var(mp_t_limit, 'mp_t_limit','[150.0:300.0]')
      END IF

      CALL chk_var(mp_dz_scal, 'mp_dz_scal', '[0.01:10.0]')
      CALL chk_var(mp_tau_d_lim, 'mp_tau_d_lim','[1.0:1.0E5]')
      CALL chk_var(mp_czero, 'mp_czero','[0.1:1000.0]')
      CALL chk_var(nbins_mp, 'nbins_mp','[1:1000]')

    END IF ! l_subgrid_qcl_mp

    IF (l_shape_rime) THEN

      CALL chk_var(a_ratio_fac, 'a_ratio_fac','[0.0:1.0]')
      CALL chk_var(a_ratio_exp, 'a_ratio_exp','[-1.0:0.0]')
      CALL chk_var(qclrime,     'qclrime','[0.0:1.0E-2]')

    END IF

    IF (l_mcr_qcf2 .AND. l_psd) THEN

      ! These do not work together and have caused jobs
      ! to hang in the past when someone set these logicals
      ! to true. Thus it is probably best to throw an error
      ! right now.

      ErrorStatus = 100

      comments = 'Generic ice PSD does not work with second ice '// newline// &
      'crystal prognostic. Please open the Rose Gui and switch  '// newline// &
      'off l_psd.'

      CALL ereport(RoutineName, ErrorStatus, comments)

    END IF ! l_mcr_qcf2 / l_psd

    IF (l_autoconv_murk .AND. .NOT. l_murk) THEN

      ! Someone has set up MURK aerosol for autoconversion without the
      ! Murk aerosol on. This will not work.
      ! Note also that the l_murk variable must be set before this checking
      ! takes place.

      ErrorStatus = 100

      comments = 'l_autoconv_murk has been set without l_murk   '// newline// &
      'i.e. requesting Murk aerosol to do some autoconversion   '// newline// &
      'from cloud to rain when it is not available. Please open '// newline// &
      'the GUI and switch either l_autoconv_murk off or l_murk on.'

      CALL ereport(RoutineName, ErrorStatus, comments)

    END IF

    IF (l_orograin) THEN

      IF (l_orograin_block) THEN
        CALL chk_var(fcrit, 'fcrit','[>0.1]')
      END IF

      CALL chk_var(nsigmasf, 'nsigmasf', '[>=0.0]')
      CALL chk_var(nscalesf, 'nscalesf', '[0.01:6.0]')

    END IF ! l_orograin

  END IF ! l_casim


END IF ! l_rain

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_run_precip

END MODULE mphys_inputs_mod

