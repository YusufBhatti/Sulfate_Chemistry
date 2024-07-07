! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mphys_constants_mod

! Description:
! Holds constants required by the large-scale
! precipitation scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

SAVE

  ! Offset between ice types (crystals and aggregates)
INTEGER, PARAMETER :: ice_type_offset = 20

  ! Arrays cx and constp hold microphysics constants
  ! Set to missing data to avoid non-initialised data in
  ! the model
REAL :: cx(200)     = rmdi
REAL :: constp(200) = rmdi

  ! Saved value of m_ci; used to check for changes from
  ! stochastic physics and if we need to rerun code
  ! to calculate microphysics constants
REAL :: m_ci_sav = 1.0

  ! air density * specific humidity * bulk velocity
  ! required for Abel and Shipway (diagnostic scheme)
REAL :: rho_q_veloc = 3.0e-4 ! This is a typical value

  ! Logical to determine whether to calculate microphysics
  ! constants (only when constants are undefined or input
  ! from stochastic physics change, causing the constants
  ! to need recalculating).
LOGICAL :: l_calc_mp_const = .TRUE.

  ! logical to use inhomogeneity parametrization for autoconv
  ! and accretion as described in Boutle et al 2012 QJ
  ! hardwired to true, but left as an option for testing
LOGICAL, PARAMETER :: l_inhomog = .TRUE.

  ! Switch turns on process of generating ice particles (sometimes known as 
  ! 'Rime Splintering') according to the Hallett-Mossop process.
  ! Between temperatures of -3 and -8 degrees Celsius, collisions between
  ! ice crystals and water (riming) produce a number of small crystals (i.e. 
  ! frozen ice splinters) when water droplets freeze on impact.
  ! Switching on this option will produce more ice crystals in the UM.
  ! Reference: Hallett,J. and Mossop,S.C. 1974. Production of secondary ice 
  ! particles during the riming process. Nature.
  ! Hardwired to false, but left as an option for testing
LOGICAL, PARAMETER :: l_hallett_mossop = .FALSE.

  ! Depositional growth by ice crystal and aggregate category can deplete 
  ! the supersaturation. With this switch true the amount of supersaturation 
  ! each ice category can deplete is multiplied by fraction of ice mass in 
  ! each category, to stop supersaturation being depleted twice over.
  ! Hardwired to false, but left as an option for testing
LOGICAL, PARAMETER :: l_cry_agg_dep = .FALSE.

!-----------------------------------------------------------------------
! Tolerances for microphysics
!-----------------------------------------------------------------------
REAL, PARAMETER :: mprog_min = 1.0E-8
! Minimum value of a moist prognostic variable to be considered by the
! microphysics

REAL, PARAMETER :: mprog_abs = 0.01 * mprog_min
! Absolute value of a moist prognostic variable; below this threshold
! it should really be considered to be zero. Currently set to two orders
! of magnitude below mprog_min

!-----------------------------------------------------------------------
! Parameters for autoconversion scheme if l_warm_new=.true.
! Currently these are for Khairoutdinov & Kogan (2000, MWR)
!-----------------------------------------------------------------------
    ! Pre-factor
REAL, PARAMETER :: acc_pref = 67.0
  ! qc power
REAL, PARAMETER :: acc_qc  = 1.15
  ! qr power
REAL, PARAMETER :: acc_qr  = 1.15

!------------------------------------

! Mass diameter relationship for ice crystals:  m(D) = AI D^BI
! Recommended values for the generic particle size distribution are
! from Brown and Francis and are identical to AI/BI
! AIC = AI = 1.85E-2, BIC = BI = 1.9.  Copy done in MODULE lspcon_mod

!------------------------------------

! Best-Reynolds numbers for aggregates
REAL, PARAMETER :: lsp_ei      = 0.2072
REAL, PARAMETER :: lsp_fi      = 0.638

!------------------------------------

! Best-Reynolds numbers for ice crystals
REAL, PARAMETER :: lsp_eic     = 6.049000e-2
REAL, PARAMETER :: lsp_fic     = 8.3100e-1

!------------------------------------

! Ice intercept values
REAL            :: x1i         = 2.0e6
REAL            :: x1ic        = 4.0e7

!------------------------------------
! Rain particle size distribution
! values
REAL, PARAMETER :: x4r         = 0.00

!------------------------------------
!autoconversion efficiency
REAL :: ec_auto     = 0.55
!------------------------------------

!Droplet number over land
REAL, PARAMETER :: ntot_land   = 3.0e8
!------------------------------------

!Droplet number over sea
REAL, PARAMETER :: ntot_sea    = 1.0e8
!------------------------------------
! Maximum droplet number assumed throught the profile

REAL, PARAMETER :: max_drop    = 375.0e6
!------------------------------------

! MURK aerosol relationships
! Scaling concentration (m-3) in aerosol number concentration
REAL :: n0_murk

! Scaling mass (kg/kg) in aerosol number calculation from aerosol mass
REAL :: m0_murk

!-----------------------------------------------------------------
! Iterations of microphysics
!-----------------------------------------------------------------

! Height at which to start iterative melting

REAL, PARAMETER :: iter_z = 3380.0

! Corresponds to level 13 in 38 level data set
! (Change added for move to 70 levels)

! Maximum iterations for iterative melting

INTEGER, PARAMETER :: max_it_mlt = 5

! Timestep of each microphysics iteration (s) when iterating over columns
! This is calculated from the number of iterations over columns,
! which in turn can be calculated from a requested sub-step length
! Note, however, that this can differ from the requested sub-step
! length as it needs to divide exactly into the model timestep
REAL :: timestep_mp = rmdi

!-----------------------------------------------------------------
! Abel & Shipway (2007) Terms
!-----------------------------------------------------------------

! Specify the default Lambda (function of rain rate)
! This value is given assuming a rain rate of 10 mm per hour.


REAL :: lam_evap_enh = 343.9988

! Maximum possible enhancement of Abel & Shipway evaporation
! Largely to avoid issues with it having no upper bound
! and hence evaporating way too too much.

REAL, PARAMETER :: max_as_enh = 50.0

! Shape dependent riming rate parameters
!
REAL :: qclmin_rime = rmdi
REAL :: area_ratio_prefac = rmdi
REAL :: area_ratio_expn = rmdi

END MODULE mphys_constants_mod
