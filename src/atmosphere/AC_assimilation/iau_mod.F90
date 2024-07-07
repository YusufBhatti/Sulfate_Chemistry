! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Parameters and variables for Incremental Analysis Update (IAU) scheme

MODULE IAU_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE filenamelength_mod, ONLY : &
    filenamelength

USE um_types
USE nlsizes_namelist_mod, ONLY: len_fixhd
USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-------------------------------------------------------------------------------
! [1]: Global constants.
!-------------------------------------------------------------------------------

INTEGER, PARAMETER :: MaxNum_IAU_incs = 9999 ! Max number of increment files

! Data for fields that may be read from the increment files:
INTEGER, PARAMETER :: IAU_NumFldCodes = 26
INTEGER, PARAMETER :: IAU_FldCodes(IAU_NumFldCodes) =         &
                        (/  2        ,  3        ,  4       , &
                            9        ,  10       ,  12      , &
                            20       ,  24       ,  90      , &
                            150      ,  253      ,  254     , &
                            255      ,  407      ,  431     , &
                            432      ,  433      ,  434     , &
                            435      ,  436      ,  480     , &
                            16207    ,  18001    ,  101     , &
                            34001    ,  34004 /)

CHARACTER(LEN=8), PARAMETER :: &
                      IAU_FldDescs(IAU_NumFldCodes) =          &
                        (/ 'u       ', 'v       ', 'theta   ', &
                           'soilm   ', 'q       ', 'qCF     ', &
                           'TSoil   ', 'sfctemp ', 'aerosol ', &
                           'w       ', 'rho     ', 'qCL     ', &
                           'exner   ', 'p       ', 'dust1   ', &
                           'dust2   ', 'dust3   ', 'dust4   ', &
                           'dust5   ', 'dust6   ', 'ozone   ', &
                           'qT      ', 'qT      ', 'SO2     ', &
                           'ukca_O3 ', 'ukca_NO2' /)

! TimeID choices, copied from VAR:
INTEGER, PARAMETER :: TimeID_Once     = 0 ! Instantaneous increment
INTEGER, PARAMETER :: TimeID_Constant = 3 ! Constant error tendency

! Call frequency codes:
INTEGER, PARAMETER :: CallFreq_Never              = 0
INTEGER, PARAMETER :: CallFreq_EveryCallIncsToAdd = 1
INTEGER, PARAMETER :: CallFreq_EndInc1Window      = 2
INTEGER, PARAMETER :: CallFreq_LastIAUCall        = 3

!-------------------------------------------------------------------------------
! [2]: Global type definitions.
!-------------------------------------------------------------------------------

! Structure to hold IAU increment data:
TYPE IAU_inc_type

  CHARACTER(filenamelength) :: FilePath

  LOGICAL :: FieldsToAdd    ! Any fields to be added to the model?

  LOGICAL :: TendencyIncs   ! Are increments (per second) tendencies?

  LOGICAL :: Contains_q_qCL_or_qCF ! Are there any increments to q, qCL or qCF?
  LOGICAL :: Contains_qT           ! Are there any increments to qT?
  LOGICAL :: Contains_sfctemp      ! Are there any increments to sfctemp?
  LOGICAL :: Contains_soilm        ! Are there any increments to soilm?

  INTEGER       :: InsertionStartTS ! Start timestep for insertion
  INTEGER       :: InsertionEndTS   ! End   timestep for insertion
  REAL, POINTER :: Weights(:)       ! Insertion weights

  INTEGER          :: FixHd(Len_FixHd)
  INTEGER          :: Len1Lookup
  INTEGER          :: Len2Lookup
  INTEGER, POINTER :: Lookup(:,:)

  LOGICAL :: StoreFields  ! Store the increment fields between IAU calls?
  LOGICAL :: StorePacked  ! Store them 32-bit-packed?
  LOGICAL :: FieldsStored ! Are the fields currently stored?

  INTEGER                    :: flds_len     ! Required length of storage array
  REAL,              POINTER :: flds     (:) ! Storage array for increments
  REAL(KIND=real32), POINTER :: flds32bit(:) ! 32-bit alternative

END TYPE IAU_inc_type

!-------------------------------------------------------------------------------
! [3]: Global variables.
!-------------------------------------------------------------------------------

TYPE (IAU_inc_type), ALLOCATABLE :: IAU_incs(:) ! IAU increments

INTEGER :: IAU_LocFldLens(IAU_NumFldCodes) ! Local lengths of fields that may
                                           ! be read from the increment files

REAL :: q_min ! Minimum value allowed for q after addition of q increments

INTEGER :: IAU_FirstCallTS = -1 ! Timestep number for first call to IAU
INTEGER :: IAU_LastCallTS  = -1 ! Timestep number for last  call to IAU

INTEGER :: ukca_O3_index  ! Index of O3  field in the tracer_ukca array
INTEGER :: ukca_NO2_index ! Index of NO2 field in the tracer_ukca array

!-------------------------------------------------------------------------------
! [3.1]: Namelist variables.
!-------------------------------------------------------------------------------

LOGICAL :: l_iau = .FALSE.  ! Activate IAU scheme?

INTEGER :: Num_IAU_incs = imdi ! Number of increment files

LOGICAL :: L_IAU_SpecifyInc1Filter = .FALSE. ! Specify filter details for first
                                            ! increment? Was TRUE

! Filter details for first increment:
INTEGER :: IAU_FilterType    = imdi   ! Filter type 1- uniform, 2- triangular,
                                      !             3- Lanczos, 4- Dolph
INTEGER :: IAU_StartMin      = imdi   ! Start minute of filtering period
                                      ! Was 0
INTEGER :: IAU_EndMin        = imdi   ! End   minute of filtering period
                                      ! Was 0
INTEGER :: IAU_ApexMin       = imdi   ! Apex minute for triangular filter
                                      ! Was 0
REAL    :: IAU_Cutoff_period = rmdi   ! Filter cutoff period in hours - was 0.0
REAL    :: IAU_SBE_period    = rmdi   ! Stop band edge period in hours - was 0.0

! Switches:
LOGICAL :: L_IAU_IncDiags       = .FALSE. ! Write increment diagnostics?
LOGICAL :: L_IAU_CalcExnerIncs  = .FALSE. ! Use p increments to calculate
                                          ! exner increments?
LOGICAL :: L_IAU_CalcThetaIncs  = .FALSE. ! Calculate theta increments using
                                          ! exner and q increments?
LOGICAL :: L_IAU_CalcRhoIncs    = .FALSE. ! Calculate rho increments using
                                          ! exner, theta and (possibly) q
                                          ! increments?
LOGICAL :: L_IAU_IncTStar       = .FALSE. ! Add level-one temperature
                                          ! increments to surface temperature
                                          ! and top-level soil temperature?
LOGICAL :: L_IAU_IncTStar_tile  = .FALSE. ! Add level-one temperature
                                          ! increments to surface temperature
                                          ! on land-surface tiles
                                          ! Was TRUE
LOGICAL :: L_IAU_IncTSurf_Snow  = .FALSE. ! Add level-one temperature
                                          ! increments to snow temperature
LOGICAL :: L_IAU_UseSfctempPerts= .FALSE. ! Include surface temperature
                                          ! perturbations and add to TStar
                                          ! and top-level soil temperature?
LOGICAL :: L_IAU_UseSoilmPerts  = .FALSE. ! Add soil moisture perturbations
                                          ! to soil-levels on land-surface?
LOGICAL :: L_IAU_UseTSoilPerts  = .FALSE. ! Add deep soil temperature 
                                          ! perturbations to soil-levels
                                          ! on land-surface?
LOGICAL :: L_IAU_ResetPoles     = .FALSE. ! Reset polar rows of relevant
                                          ! increment fields to their mean
                                          ! values? Was TRUE
LOGICAL :: L_IAU_RemoveSS       = .FALSE. ! Remove supersaturation wrt water?
LOGICAL :: L_IAU_DumpTS0State   = .FALSE. ! Write out model state
                                          ! immediately after timestep-zero
                                          ! call to IAU scheme
LOGICAL :: L_IAU_ApplyQTCorrections = .FALSE.     ! Apply corrections for
                                                  ! processing of qT incs.
LOGICAL :: L_IAU_Add_qT_prime_to_q = .FALSE. ! Replace q with qT_plus,
                                                  ! bypassing Var_DiagCloud.
LOGICAL :: L_IAU_IncrementIce   = .FALSE. ! Diagnose ice cloud increments as
                                          ! well as humidity and liquid cloud
                                          ! increments from qT increments?
LOGICAL :: L_IAU_ScaleCloud     = .FALSE. ! If diagnosing humidity and cloud
                                          ! increments from qT increments,
                                          ! scale increments to be within
                                          ! physical bounds?
LOGICAL :: L_IAU_LimitUpperThetaIncs = .FALSE. ! Apply limits to size of
                                               ! upper-level theta increments?
LOGICAL :: L_IAU_SetOzoneMin    = .FALSE. ! Reset ozone to oz_min if ozone was
                                          ! negative?
LOGICAL :: L_IAU_SMCChecks      = .FALSE. ! Apply limits and masks to SMC incs?
LOGICAL :: L_IAU_TSoilLandIceMask = .FALSE. ! Mask out TSoil incs at land-ice
                                            ! points?

! Parameters for use with L_IAU_LimitUpperThetaIncs:
REAL :: IAU_LimitUpperThetaIncs_pBound = rmdi ! Pressure boundary - was 200.0
REAL :: IAU_LimitUpperThetaIncs_maxInc = rmdi ! Maximum absolute value of
                                               ! increment after multiplication
                                               ! by IAU weight - was 100.0

! Parameters and switches for use with QLimits:
INTEGER :: IAU_QLimitsCallFreq  = imdi    ! Call frequency
                                          !  - was CallFreq_EndInc1Window
LOGICAL :: L_IAU_QLimitsDiags   = .FALSE. ! Write diagnostics?
LOGICAL :: L_IAU_RmNonTropQIncs = .FALSE. ! Remove non-trop q increments?
REAL    :: IAU_trop_min_RH      = rmdi    ! Lower limit to apply to trop RH
                                          ! Was 0.0
REAL    :: IAU_nonTrop_max_q    = rmdi    ! Upper limit to apply to non-trop q
                                          ! Was 3.0E-06
REAL    :: IAU_nonTrop_min_q    = rmdi    ! Lower limit to apply to non-trop q
                                          ! Was 1.0E-06
REAL    :: IAU_nonTrop_max_RH   = rmdi    ! Upper limit to apply to non-trop RH
                                          ! Was 0.1
REAL    :: IAU_trop_min_p       = rmdi    ! Minimum tropospheric pressure
                                          ! Was 1.0E+04
REAL    :: IAU_trop_max_PV      = rmdi    ! Maximum tropospheric ABS(PV)
                                          ! Was 2.5E-06
REAL    :: IAU_nonTrop_max_p    = rmdi    ! Maximum non-trop     pressure
                                          ! Was 4.0E+04

! Want to make this switch obsolete ASAP:
LOGICAL :: L_IAU_IgnoreTopLevThetaIncs = .FALSE.

! Parameters and switches for use within Var_DiagCloud itself.
! (Currently Var_DiagCloud is called only from IAU, so it makes sense to
! control them via the IAU namelist.)

REAL :: DiagCloud_Tol_q = rmdi  ! Was 0.0
  ! Tolerance for qcl==0 and qcf==0 tests in VarDiag_Cloud.
  ! (Tests suggest that a value of 1.0e-14 is more-or-less optimal.)

REAL :: DiagCloud_Tol_FM = rmdi  ! Was 0.0
  ! Tolerance for (1-FM)==0  tests in VarDiag_Cloud.
  ! (Tests suggest that a value of 1.0e-28 is more-or-less optimal.)

LOGICAL :: DiagCloud_ApplyCompLimits = .FALSE. ! Was TRUE
  ! If true, take measures to avoid machine precision issues in Var_DiagCloud.
  ! (Once we're confident that this option has no nasty side-effects, we should
  ! remove the switch and always use it.)

INTEGER :: DiagCloud_NumMaxLoops = imdi  ! Was 150
  ! Maximum number of loops for convergence in Var_DiagCloud.

REAL    :: DiagCloud_QN_CompRegimeLimit = rmdi ! Was 20.0
  ! When the magnitude of QN gets much beyond this parameter, the standard
  ! calculations for qCL_calc in Var_DiagCloud give deviations from its upper
  ! or lower bound that are highly sensitive to rounding. In the case where QN
  ! is negative, this makes the iterative algorithm used to calculate cloud
  ! quantities extremely sensitive to small changes in the input data, and when
  ! the IncrementIce option is used often leads to non-convergence. To get
  ! around this, when ApplyCompBounds is set to true, and the magnitude of QN
  ! exceeds the above limit, we set qCL_calc directly to the appropriate
  ! bounding value.
  !
  ! QN < -QN_CompRegimeLimit occurs at points very far from any saturation.
  ! QN >  QN_CompRegimeLimit occurs at points very far from any unsaturation.
  !
  ! The parameter value was determined experimentally by determining when the
  ! path through the code starts becoming affected by least-significant-bit
  ! changes to the LS states. The current value is appropriate for 64-bit
  ! reals.

!-------------------------------------------------------------------------------
! [4]: IAU namelist.
!-------------------------------------------------------------------------------

NAMELIST / IAU_nl /             &
l_iau,                          &
Num_IAU_incs,                   &
L_IAU_SpecifyInc1Filter,        &
IAU_FilterType,                 &
IAU_StartMin,                   &
IAU_EndMin,                     &
IAU_ApexMin,                    &
IAU_Cutoff_period,              &
IAU_SBE_period,                 &
L_IAU_IncDiags,                 &
L_IAU_CalcExnerIncs,            &
L_IAU_CalcThetaIncs,            &
L_IAU_CalcRhoIncs,              &
L_IAU_IncTStar,                 &
L_IAU_IncTStar_tile,            &
L_IAU_IncTSurf_Snow,            &
L_IAU_UseSoilmPerts,            &
L_IAU_UseTSoilPerts,            &
L_IAU_UseSfctempPerts,          &
L_IAU_ResetPoles,               &
L_IAU_RemoveSS,                 &
L_IAU_DumpTS0State,             &
L_IAU_ApplyQTCorrections,       &
L_IAU_Add_qT_prime_to_q,        &
L_IAU_IncrementIce,             &
L_IAU_ScaleCloud,               &
L_IAU_LimitUpperThetaIncs,      &
L_IAU_SetOzoneMin,              &
L_IAU_SMCChecks,                &
L_IAU_TSoilLandIceMask,         &
IAU_LimitUpperThetaIncs_pBound, &
IAU_LimitUpperThetaIncs_maxInc, &
IAU_QLimitsCallFreq,            &
L_IAU_QLimitsDiags,             &
L_IAU_RmNonTropQIncs,           &
IAU_trop_min_RH,                &
IAU_nonTrop_max_q,              &
IAU_nonTrop_min_q,              &
IAU_nonTrop_max_RH,             &
IAU_trop_min_p,                 &
IAU_trop_max_PV,                &
IAU_nonTrop_max_p,              &
L_IAU_IgnoreTopLevThetaIncs,    &
DiagCloud_Tol_q,                &
DiagCloud_Tol_FM,               &
DiagCloud_ApplyCompLimits,      &
DiagCloud_NumMaxLoops,          & ! must be greater than 50
DiagCloud_QN_CompRegimeLimit

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IAU_MOD'

CONTAINS

#if !defined(LFRIC)
SUBROUTINE read_nml_iau_nl(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER :: icode, errorstatus
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_IAU_NL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7
INTEGER, PARAMETER :: n_real = 14
INTEGER, PARAMETER :: n_log = 27

TYPE my_namelist
  SEQUENCE
  INTEGER :: Num_IAU_incs
  INTEGER :: IAU_FilterType
  INTEGER :: IAU_StartMin
  INTEGER :: IAU_EndMin
  INTEGER :: IAU_ApexMin
  INTEGER :: IAU_QLimitsCallFreq
  INTEGER :: DiagCloud_NumMaxLoops
  REAL ::  IAU_Cutoff_period
  REAL ::  IAU_SBE_period
  REAL ::  IAU_LimitUpperThetaIncs_pBound
  REAL ::  IAU_LimitUpperThetaIncs_maxInc
  REAL ::  IAU_trop_min_RH
  REAL ::  IAU_nonTrop_max_q
  REAL ::  IAU_nonTrop_min_q
  REAL ::  IAU_nonTrop_max_RH
  REAL ::  IAU_trop_min_p
  REAL ::  IAU_trop_max_PV
  REAL ::  IAU_nonTrop_max_p
  REAL ::  DiagCloud_Tol_q
  REAL ::  DiagCloud_Tol_FM
  REAL ::  DiagCloud_QN_CompRegimeLimit
  LOGICAL ::  l_iau
  LOGICAL ::  L_IAU_SpecifyInc1Filter
  LOGICAL ::  L_IAU_IncDiags
  LOGICAL ::  L_IAU_CalcExnerIncs
  LOGICAL ::  L_IAU_CalcThetaIncs
  LOGICAL ::  L_IAU_CalcRhoIncs
  LOGICAL ::  L_IAU_IncTStar
  LOGICAL ::  L_IAU_IncTStar_tile
  LOGICAL ::  L_IAU_IncTSurf_Snow
  LOGICAL ::  L_IAU_UseSoilmPerts
  LOGICAL ::  L_IAU_UseTSoilPerts
  LOGICAL ::  L_IAU_UseSfctempPerts
  LOGICAL ::  L_IAU_ResetPoles
  LOGICAL ::  L_IAU_RemoveSS
  LOGICAL ::  L_IAU_DumpTS0State
  LOGICAL ::  L_IAU_ApplyQTCorrections
  LOGICAL ::  L_IAU_Add_qT_prime_to_q
  LOGICAL ::  L_IAU_IncrementIce
  LOGICAL ::  L_IAU_ScaleCloud
  LOGICAL ::  L_IAU_LimitUpperThetaIncs
  LOGICAL ::  L_IAU_SetOzoneMin
  LOGICAL ::  L_IAU_SMCChecks
  LOGICAL ::  L_IAU_TSoilLandIceMask
  LOGICAL ::  L_IAU_QLimitsDiags
  LOGICAL ::  L_IAU_RmNonTropQIncs
  LOGICAL ::  L_IAU_IgnoreTopLevThetaIncs
  LOGICAL ::  DiagCloud_ApplyCompLimits
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,  &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=IAU_nl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist IAU_nl", iomessage)

  my_nml % Num_IAU_incs          = Num_IAU_incs
  my_nml % IAU_FilterType        = IAU_FilterType
  my_nml % IAU_StartMin          = IAU_StartMin
  my_nml % IAU_EndMin            = IAU_EndMin
  my_nml % IAU_ApexMin           = IAU_ApexMin
  my_nml % IAU_QLimitsCallFreq   = IAU_QLimitsCallFreq
  my_nml % DiagCloud_NumMaxLoops = DiagCloud_NumMaxLoops
  ! end of integers
  my_nml % IAU_Cutoff_period  = IAU_Cutoff_period
  my_nml % IAU_SBE_period     = IAU_SBE_period
  my_nml % IAU_LimitUpperThetaIncs_pBound = IAU_LimitUpperThetaIncs_pBound
  my_nml % IAU_LimitUpperThetaIncs_maxInc = IAU_LimitUpperThetaIncs_maxInc
  my_nml % IAU_trop_min_RH    = IAU_trop_min_RH
  my_nml % IAU_nonTrop_max_q  = IAU_nonTrop_max_q
  my_nml % IAU_nonTrop_min_q  = IAU_nonTrop_min_q
  my_nml % IAU_nonTrop_max_RH = IAU_nonTrop_max_RH
  my_nml % IAU_trop_min_p     = IAU_trop_min_p
  my_nml % IAU_trop_max_PV    = IAU_trop_max_PV
  my_nml % IAU_nonTrop_max_p  = IAU_nonTrop_max_p
  my_nml % DiagCloud_Tol_q    = DiagCloud_Tol_q
  my_nml % DiagCloud_Tol_FM   = DiagCloud_Tol_FM
  my_nml % DiagCloud_QN_CompRegimeLimit = DiagCloud_QN_CompRegimeLimit
  ! end of reals
  my_nml % l_iau                   = l_iau
  my_nml % L_IAU_SpecifyInc1Filter = L_IAU_SpecifyInc1Filter
  my_nml % L_IAU_IncDiags          = L_IAU_IncDiags
  my_nml % L_IAU_CalcExnerIncs     = L_IAU_CalcExnerIncs
  my_nml % L_IAU_CalcThetaIncs     = L_IAU_CalcThetaIncs
  my_nml % L_IAU_CalcRhoIncs       = L_IAU_CalcRhoIncs
  my_nml % L_IAU_IncTStar          = L_IAU_IncTStar
  my_nml % L_IAU_IncTStar_tile     = L_IAU_IncTStar_tile
  my_nml % L_IAU_IncTSurf_Snow     = L_IAU_IncTSurf_Snow
  my_nml % L_IAU_UseSoilmPerts     = L_IAU_UseSoilmPerts
  my_nml % L_IAU_UseTSoilPerts     = L_IAU_UseTSoilPerts
  my_nml % L_IAU_UseSfctempPerts   = L_IAU_UseSfctempPerts
  my_nml % L_IAU_ResetPoles        = L_IAU_ResetPoles
  my_nml % L_IAU_RemoveSS          = L_IAU_RemoveSS
  my_nml % L_IAU_DumpTS0State      = L_IAU_DumpTS0State
  my_nml % L_IAU_ApplyQTCorrections = L_IAU_ApplyQTCorrections
  my_nml % L_IAU_Add_qT_prime_to_q = L_IAU_Add_qT_prime_to_q
  my_nml % L_IAU_IncrementIce      = L_IAU_IncrementIce
  my_nml % L_IAU_ScaleCloud        = L_IAU_ScaleCloud
  my_nml % L_IAU_LimitUpperThetaIncs = L_IAU_LimitUpperThetaIncs
  my_nml % L_IAU_SetOzoneMin       = L_IAU_SetOzoneMin
  my_nml % L_IAU_SMCChecks         = L_IAU_SMCChecks
  my_nml % L_IAU_TSoilLandIceMask  = L_IAU_TSoilLandIceMask
  my_nml % L_IAU_QLimitsDiags      = L_IAU_QLimitsDiags
  my_nml % L_IAU_RmNonTropQIncs    = L_IAU_RmNonTropQIncs
  my_nml % L_IAU_IgnoreTopLevThetaIncs = L_IAU_IgnoreTopLevThetaIncs
  my_nml % DiagCloud_ApplyCompLimits = DiagCloud_ApplyCompLimits

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  Num_IAU_incs          = my_nml % Num_IAU_incs
  IAU_FilterType        = my_nml % IAU_FilterType
  IAU_StartMin          = my_nml % IAU_StartMin
  IAU_EndMin            = my_nml % IAU_EndMin
  IAU_ApexMin           = my_nml % IAU_ApexMin
  IAU_QLimitsCallFreq   = my_nml % IAU_QLimitsCallFreq
  DiagCloud_NumMaxLoops = my_nml % DiagCloud_NumMaxLoops
  ! end of integers
  IAU_Cutoff_period  = my_nml % IAU_Cutoff_period
  IAU_SBE_period     = my_nml % IAU_SBE_period
  IAU_LimitUpperThetaIncs_pBound = my_nml % IAU_LimitUpperThetaIncs_pBound
  IAU_LimitUpperThetaIncs_maxInc = my_nml % IAU_LimitUpperThetaIncs_maxInc
  IAU_trop_min_RH    = my_nml % IAU_trop_min_RH
  IAU_nonTrop_max_q  = my_nml % IAU_nonTrop_max_q
  IAU_nonTrop_min_q  = my_nml % IAU_nonTrop_min_q
  IAU_nonTrop_max_RH = my_nml % IAU_nonTrop_max_RH
  IAU_trop_min_p     = my_nml % IAU_trop_min_p
  IAU_trop_max_PV    = my_nml % IAU_trop_max_PV
  IAU_nonTrop_max_p  = my_nml % IAU_nonTrop_max_p
  DiagCloud_Tol_q    = my_nml % DiagCloud_Tol_q
  DiagCloud_Tol_FM   = my_nml % DiagCloud_Tol_FM
  DiagCloud_QN_CompRegimeLimit = my_nml % DiagCloud_QN_CompRegimeLimit
  ! end of reals
  l_iau                   = my_nml % l_iau
  L_IAU_SpecifyInc1Filter = my_nml % L_IAU_SpecifyInc1Filter
  L_IAU_IncDiags          = my_nml % L_IAU_IncDiags
  L_IAU_CalcExnerIncs     = my_nml % L_IAU_CalcExnerIncs
  L_IAU_CalcThetaIncs     = my_nml % L_IAU_CalcThetaIncs
  L_IAU_CalcRhoIncs       = my_nml % L_IAU_CalcRhoIncs
  L_IAU_IncTStar          = my_nml % L_IAU_IncTStar
  L_IAU_IncTStar_tile     = my_nml % L_IAU_IncTStar_tile
  L_IAU_IncTSurf_Snow     = my_nml % L_IAU_IncTSurf_Snow
  L_IAU_UseSoilmPerts     = my_nml % L_IAU_UseSoilmPerts
  L_IAU_UseTSoilPerts     = my_nml % L_IAU_UseTSoilPerts
  L_IAU_UseSfctempPerts   = my_nml % L_IAU_UseSfctempPerts
  L_IAU_ResetPoles        = my_nml % L_IAU_ResetPoles
  L_IAU_RemoveSS          = my_nml % L_IAU_RemoveSS
  L_IAU_DumpTS0State      = my_nml % L_IAU_DumpTS0State
  L_IAU_ApplyQTCorrections = my_nml % L_IAU_ApplyQTCorrections
  L_IAU_Add_qT_prime_to_q = my_nml % L_IAU_Add_qT_prime_to_q
  L_IAU_IncrementIce      = my_nml % L_IAU_IncrementIce
  L_IAU_ScaleCloud        = my_nml % L_IAU_ScaleCloud
  L_IAU_LimitUpperThetaIncs = my_nml % L_IAU_LimitUpperThetaIncs
  L_IAU_SetOzoneMin       = my_nml % L_IAU_SetOzoneMin
  L_IAU_SMCChecks         = my_nml % L_IAU_SMCChecks
  L_IAU_TSoilLandIceMask  = my_nml % L_IAU_TSoilLandIceMask
  L_IAU_QLimitsDiags      = my_nml % L_IAU_QLimitsDiags
  L_IAU_RmNonTropQIncs    = my_nml % L_IAU_RmNonTropQIncs
  L_IAU_IgnoreTopLevThetaIncs = my_nml % L_IAU_IgnoreTopLevThetaIncs
  DiagCloud_ApplyCompLimits = my_nml % DiagCloud_ApplyCompLimits

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_iau_nl
#endif

END MODULE IAU_mod
