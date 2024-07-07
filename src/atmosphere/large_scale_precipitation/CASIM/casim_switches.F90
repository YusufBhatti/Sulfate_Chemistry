! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Cloud Aerosol Interacting Microphysics (CASIM)
! Module to hold the switches for CASIM- temporary until CASIM
! Microphysics is lodged on the trunk.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation (CASIM)

MODULE casim_switches

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

! Switches for the activated moments
LOGICAL :: l_mp_cloudnumber       = .FALSE.
LOGICAL :: l_mp_rainnumber        = .FALSE.
LOGICAL :: l_mp_rain3mom          = .FALSE.
LOGICAL :: l_mp_icenumber         = .FALSE.
LOGICAL :: l_mp_snownumber        = .FALSE.
LOGICAL :: l_mp_snow3mom          = .FALSE.
LOGICAL :: l_mp_graupnumber       = .FALSE.
LOGICAL :: l_mp_graup3mom         = .FALSE.

! Switches for the activated aerosols
LOGICAL :: l_mp_activesolliquid   = .FALSE.
LOGICAL :: l_mp_activesolrain     = .FALSE.
LOGICAL :: l_mp_activeinsolice    = .FALSE.
LOGICAL :: l_mp_activesolice      = .FALSE.
LOGICAL :: l_mp_activeinsolliquid = .FALSE.
LOGICAL :: l_mp_activesolnumber   = .FALSE.
LOGICAL :: l_mp_activeinsolnumber = .FALSE.

! Switches to control the aerosol coupling with CASIM
LOGICAL :: l_fix_aerosol          = .FALSE.
LOGICAL :: l_tracer_aerosol       = .FALSE.
LOGICAL :: l_ukca_aerosol         = .FALSE.
LOGICAL :: l_ukca_feeding_in      = .FALSE.
LOGICAL :: l_ukca_feeding_out     = .FALSE.

LOGICAL :: l_asci_progs           = .FALSE.

! Switch for a warm rain only simulation
LOGICAL :: l_casim_warm_only      = .FALSE.

! Number of extra CASIM prognostics
INTEGER :: n_casim_progs          = imdi

! CASIM moments option
INTEGER :: casim_moments_option   = imdi

! Initialise cloud and ice number the LBC with CASIM
LOGICAL :: l_init_lbc_number      = .FALSE.

! Switch for a diagnostic cloud scheme (under development)
LOGICAL :: l_cfrac_casim_diag_scheme = .FALSE.

! Switch for CASIM diagnostic call
LOGICAL :: l_casimmphys_diags = .FALSE.

! Switch for CASIM microphysics in the rim (should be removed)
LOGICAL :: l_micro_in_rim = .FALSE.

! CASIM moments - defined in casim_set_dependent_switches and
! for use in the interface
INTEGER :: cloud_mom
INTEGER :: rain_mom
INTEGER :: ice_mom
INTEGER :: snow_mom
INTEGER :: graup_mom

! Parameters for each moment type (single, double, triple)
INTEGER, PARAMETER :: no_moments    = 0
INTEGER, PARAMETER :: single_moment = 1
INTEGER, PARAMETER :: double_moment = 2
INTEGER, PARAMETER :: triple_moment = 3

! CASIM dimensions
INTEGER :: its   
INTEGER :: ite
INTEGER :: jts
INTEGER :: jte
INTEGER :: kts
INTEGER :: kte

! First and last dimensions over which CASIM is actually run. 
! If l_micro_in_rim = .TRUE., they are equivalent to tdims
! except kle is tdims%k_end - 2
! If l_micro_in_rim = .FALSE., they are the same as the
! rim dimensions below. 
INTEGER :: ils   
INTEGER :: ile
INTEGER :: jls
INTEGER :: jle
INTEGER :: kls
INTEGER :: kle

!----------------------------------------------------------
! Dimensions without the rim at the start and end. 
! So for a LAM these would be tdims%i_start + rimwidth and 
! tdims%i_end - rimwidth etc 

! Used for CASIM unless l_micro_in_rim = .TRUE.
INTEGER :: irs
INTEGER :: ire
INTEGER :: jrs
INTEGER :: jre
INTEGER :: krs
INTEGER :: kre

!--------------------------------------------------------
! Indexes for aerosol tracers used in the interface
!-------------------------------------------------------
! Soluble aerosol
INTEGER, PARAMETER :: i_AitkenSolMass    = 1
INTEGER, PARAMETER :: i_AitkenSolNumber  = 2
INTEGER, PARAMETER :: i_AccumSolMass     = 3
INTEGER, PARAMETER :: i_AccumSolNumber   = 4
INTEGER, PARAMETER :: i_CoarseSolMass    = 5
INTEGER, PARAMETER :: i_CoarseSolNumber  = 6

! Insoluble aerosol
INTEGER, PARAMETER :: i_CoarseDustMass   = 7
INTEGER, PARAMETER :: i_CoarseDustNumber = 8
INTEGER, PARAMETER :: i_AccumDustMass    = 9
INTEGER, PARAMETER :: i_AccumDustNumber  = 10

INTEGER, PARAMETER :: n_casim_tracers = 10

! Parameters for casim_aerosol_process_level
INTEGER, PARAMETER :: no_processing       = 0
INTEGER, PARAMETER :: passive_processing  = 1
INTEGER, PARAMETER :: full_processing     = 2
INTEGER, PARAMETER :: passive_ice_only    = 3 
INTEGER, PARAMETER :: passive_liquid_only = 4

END MODULE casim_switches
