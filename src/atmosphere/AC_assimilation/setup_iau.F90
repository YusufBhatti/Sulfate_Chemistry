! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  IAU set-up routine

MODULE setup_iau_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETUP_IAU_MOD'

CONTAINS

SUBROUTINE Setup_IAU

! Description:
!
!   Basically, set up a data structure for each of the IAU increment files
!   containing the IAU weights etc., and trap any setup errors.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE IAU_mod, ONLY:                 &
    MaxNum_IAU_incs,                &
    IAU_NumFldCodes,                &
    IAU_FldCodes,                   &
    IAU_FldDescs,                   &
    TimeID_Once,                    &
    TimeID_Constant,                &
    CallFreq_Never,                 &
    CallFreq_EveryCallIncsToAdd,    &
    CallFreq_EndInc1Window,         &
    CallFreq_LastIAUCall,           &
    IAU_incs,                       &
    IAU_LocFldLens,                 &
    q_min,                          &
    IAU_FirstCallTS,                &
    IAU_LastCallTS,                 &
    ukca_O3_index,                  &
    ukca_NO2_index,                 &
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
    L_IAU_UseSfctempPerts,          &
    L_IAU_UseSoilmPerts,            &
    L_IAU_UseTSoilPerts,            &
    L_IAU_ResetPoles,               &
    L_IAU_DumpTS0State,             &
    L_IAU_RemoveSS,                 &
    L_IAU_ApplyQTCorrections,       &
    L_IAU_Add_qT_prime_to_q,        &
    L_IAU_IncrementIce,             &
    L_IAU_ScaleCloud,               &
    L_IAU_LimitUpperThetaIncs,      &
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
    L_IAU_TSoilLandIceMask,         &
    L_IAU_SMCChecks

USE nlstcall_mod, ONLY: lcal360 
USE nlstgen_mod, ONLY: i_dump_output
USE yomhook, ONLY: lhook, dr_hook
USE io, ONLY:        &
    setpos, buffin
USE io_constants, ONLY: &
    ioOpenReadOnly, ioNameProvided, ioNoDelete
USE model_file, ONLY: model_file_open, model_file_close
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE Control_Max_Sizes
USE dust_parameters_mod, ONLY: l_twobin_dust
USE lookup_addresses
USE turb_diff_mod, ONLY: l_qpos, qlimit
USE nlcfiles_namelist_mod, ONLY: iau_inc_file => iau_inc
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE calc_tfiltwts_mod, ONLY: calc_tfiltwts
USE submodel_mod, ONLY: atmos_im
USE nlsizes_namelist_mod, ONLY: &
    len_fixhd, theta_field_size, u_field_size, v_field_size, model_levels

USE conversions_mod, ONLY: isec_per_day
USE errormessagelength_mod, ONLY: errormessagelength

USE model_time_mod, ONLY: &
    basis_time_days, basis_time_secs, secs_per_stepim

USE model_domain_mod, ONLY: model_type, mt_global

USE run_aerosol_mod, ONLY: l_sulpc_so2
USE ukca_nmspec_mod, ONLY: ukca_name2index
USE ukca_option_mod, ONLY: l_ukca
USE um_stashcode_mod, ONLY: stashcode_O3, stashcode_NO2, stashcode_SO2

USE packing_codes_mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE


! DEPENDS ON: read_flh
! DEPENDS ON: ioerror
! DEPENDS ON: time2sec

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SETUP_IAU'

! Local variables:

INTEGER :: i
INTEGER :: ICode
INTEGER :: TS_len_secs
INTEGER :: Code
INTEGER :: Level
INTEGER :: LenIO
INTEGER :: LenInc1Path
INTEGER :: LenInc1Name
INTEGER :: LenNumIncsStr
INTEGER :: LenIncNumStr
INTEGER :: FieldNum
INTEGER :: Lookup_len
INTEGER :: IncNum
INTEGER :: LocFldLen
INTEGER :: pack_code
INTEGER :: flds_len
INTEGER :: TS_num
INTEGER :: TimeId
INTEGER :: VAR_StartSec
INTEGER :: VAR_EndSec
INTEGER :: UM_StartSec
INTEGER :: UM_EndSec
INTEGER :: ElapsedDays
INTEGER :: ElapsedSecs
INTEGER :: DTSecs
INTEGER :: DTMins
INTEGER :: VTSecs
INTEGER :: VTMins
INTEGER :: NumWeights
INTEGER :: WeightNum
INTEGER :: FirstNonZeroWeightIndex
INTEGER :: LastNonZeroWeightIndex
INTEGER :: Sec
INTEGER :: IAUDustBinNum
INTEGER :: IAUNumDustBins
INTEGER :: IAU_unit

INTEGER :: VAR_StartSecSave(MaxNum_IAU_incs)
INTEGER :: VAR_EndSecSave  (MaxNum_IAU_incs)
INTEGER :: UM_StartSecSave (MaxNum_IAU_incs)
INTEGER :: UM_EndSecSave   (MaxNum_IAU_incs)

REAL :: a_io

REAL, ALLOCATABLE :: Weights(:)

LOGICAL :: contains_new_qT_code
LOGICAL :: contains_old_qT_code
LOGICAL :: AllPacked
LOGICAL :: Using_q_qCL_or_qCF
LOGICAL :: Using_qT
LOGICAL :: IAUHasDustBin(6)
LOGICAL :: Found_SO2 = .FALSE.
LOGICAL :: Found_ukca_O3 = .FALSE.
LOGICAL :: Found_ukca_NO2 = .FALSE.

CHARACTER(LEN=4)   :: NumIncsStr
CHARACTER(LEN=4)   :: IncNumStr
CHARACTER(LEN=9)   :: VAR_StartStr
CHARACTER(LEN=9)   :: VAR_EndStr
CHARACTER(LEN=9)   :: UM_StartStr
CHARACTER(LEN=9)   :: UM_EndStr
CHARACTER(LEN=9)   :: TimeStr
CHARACTER(LEN=9)   :: IAUNumDustBinsStr
CHARACTER(LEN=errormessagelength) :: CMessage = ' '
CHARACTER(LEN=errormessagelength) :: tempmessage = ' '

CHARACTER(LEN=60)  :: UsableFlds(MaxNum_IAU_incs)

CHARACTER(filenamelength) :: FilePath
CHARACTER(filenamelength) :: Inc1Path

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

!-------------------------------------------------------------------------------
! [1]: Check namelist variables.
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (Num_IAU_incs > MaxNum_IAU_incs) THEN
  ICode = 1
  WRITE (CMessage,'(A,I12,A,I4,A)')                                         &
    'Number of IAU increment files (', Num_IAU_incs, ') exceeds maximum (', &
    MaxNum_IAU_incs, ')'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

IF (IAU_QLimitsCallFreq /= CallFreq_Never              .AND. &
    IAU_QLimitsCallFreq /= CallFreq_EveryCallIncsToAdd .AND. &
    IAU_QLimitsCallFreq /= CallFreq_EndInc1Window      .AND. &
    IAU_QLimitsCallFreq /= CallFreq_LastIAUCall) THEN
  ICode = 1
  WRITE (CMessage,'(A,I6,A)') &
    'Value of IAU_QLimitsCallFreq (', IAU_QLimitsCallFreq, ') unsupported'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

IF (L_IAU_RmNonTropQIncs .AND. IAU_QLimitsCallFreq == CallFreq_Never) THEN
  ICode = 1
  CMessage = 'L_IAU_RmNonTropQIncs set, but not calling QLimits'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

IF (i_dump_output == 1 .AND. L_IAU_DumpTS0State) THEN
  ICode = -1
  CMessage = 'i_dump_output=1 (no dumping or meaning), so resetting '// &
             'L_IAU_DumpTS0State to .FALSE.'
  CALL EReport (RoutineName, ICode, CMessage)
  L_IAU_DumpTS0State = .FALSE.
END IF

!-------------------------------------------------------------------------------
! [2]: Set up q_min, depending on whether q_pos routine is being used.
!-------------------------------------------------------------------------------

IF (l_qpos) THEN
  q_min = qlimit
ELSE
  q_min = 1.0e-8
END IF

!-------------------------------------------------------------------------------
! [3]: Get local lengths of fields that may be read from the increment files.
!-------------------------------------------------------------------------------

DO i = 1, IAU_NumFldCodes

  Code = IAU_FldCodes(i)

  ! Default:
  IAU_LocFldLens(i) = theta_field_size

  ! Exceptions:
  IF (Code == 2) IAU_LocFldLens(i) = u_field_size ! u
  IF (Code == 3) IAU_LocFldLens(i) = v_field_size ! v

  IF (Code == 4 .AND. L_IAU_CalcThetaIncs)          &
    IAU_LocFldLens(i) = 0 ! Don't read theta

  IF (Code == 9 .AND. .NOT. L_IAU_UseSoilmPerts)    &
    IAU_LocFldLens(i) = 0 ! Don't read soilm perts

  IF (Code == 20 .AND. .NOT. L_IAU_UseTSoilPerts)   &
    IAU_LocFldLens(i) = 0 ! Don't read T Soil perts

  IF (Code == 24 .AND. .NOT. L_IAU_UseSfctempPerts) &
    IAU_LocFldLens(i) = 0 ! Don't read sfctemp perts

  IF (Code == 253 .AND. L_IAU_CalcRhoIncs)          &
    IAU_LocFldLens(i) = 0 ! Don't read rho

  IF (Code == 255 .AND. L_IAU_CalcExnerIncs)        &
    IAU_LocFldLens(i) = 0 ! Don't read exner

  IF (Code == 407 .AND. .NOT. L_IAU_CalcExnerIncs)   &
    IAU_LocFldLens(i) = 0 ! Don't read p

END DO

!-------------------------------------------------------------------------------
! [4]: Print IAU details common to all increments.
!-------------------------------------------------------------------------------

IF (PrintStatus >= PrStatus_Normal) THEN
  CALL umPrint( '',src='setup_iau')
  CALL umPrint( '<><><><><><><><><><><><><><><><><><><><><><><><><>',          &
      src='setup_iau')
  CALL umPrint( 'Start of Incremental Analysis Update (IAU) details',          &
      src='setup_iau')
  CALL umPrint( '',src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Write increment diagnostics? ',L_IAU_IncDiags
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Apply qT corrections?        ',                   &
      L_IAU_ApplyQTCorrections
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Replace q with q+qTprime?      ',                 &
      L_IAU_Add_qT_prime_To_q
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Calculate ice cloud incs?    ',L_IAU_IncrementIce
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Scale cloud incs?            ',L_IAU_ScaleCloud
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Calculate exner incs?        ',L_IAU_CalcExnerIncs
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Calculate theta incs?        ',L_IAU_CalcThetaIncs
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Calculate rho   incs?        ',L_IAU_CalcRhoIncs
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Ignore top-lev theta incs?   ',                   &
      L_IAU_IgnoreTopLevThetaIncs
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Remove supersaturation?      ',L_IAU_RemoveSS
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Increment TStar/TSoil?       ', L_IAU_IncTStar
  CALL umPrint(umMessage,src='setup_iau')
  WRITE(umMessage,'(A,L6)') 'Add TSoil Perturbations?     ', L_IAU_UseTSoilPerts
  CALL umPrint(umMessage,src='setup_iau')
  IF (L_IAU_UseTSoilPerts) THEN
    WRITE(umMessage,'(A,L6)') ' - Mask out TSoil incs at land-ice points?',   &
                              L_IAU_TSoilLandIceMask
    CALL umPrint(umMessage,src='setup_iau')
  END IF
  WRITE(umMessage,'(A,L6)') 'Add Soilm Perturbations?     ', L_IAU_UseSoilmPerts
  CALL umPrint(umMessage,src='setup_iau')
  IF (L_IAU_UseSoilmPerts) THEN
    WRITE(umMessage,'(A,L6)') ' - Apply limits and masks to SMC increments?', &
                              L_IAU_SMCChecks
    CALL umPrint(umMessage,src='setup_iau')
  END IF
  IF (model_type == mt_global) THEN
    WRITE(umMessage,'(A,L6)') 'Reset polar rows?            ', L_IAU_ResetPoles
    CALL umPrint(umMessage,src='setup_iau')
  END IF
  WRITE(umMessage,'(A,L6)') 'Limit upper-level theta incs?',                   &
      L_IAU_LimitUpperThetaIncs
  CALL umPrint(umMessage,src='setup_iau')
  IF (L_IAU_LimitUpperThetaIncs) THEN
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Pressure boundary: ', IAU_LimitUpperThetaIncs_pBound
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Max abs increment: ', IAU_LimitUpperThetaIncs_maxInc
    CALL umPrint(umMessage,src='setup_iau')
  END IF
  SELECT CASE (IAU_QLimitsCallFreq)
  CASE (CallFreq_Never)
    CALL umPrint( 'QLimits call frequency:       Never',src='setup_iau')
  CASE (CallFreq_EveryCallIncsToAdd)
    CALL umPrint( 'QLimits call frequency:       Every call on which '//     &
        'there are incs to add',   &
        src='setup_iau')
  CASE (CallFreq_EndInc1Window)
    CALL umPrint( 'QLimits call frequency:       End of Inc1 time window',   &
        src='setup_iau')
  CASE (CallFreq_LastIAUCall)
    CALL umPrint( 'QLimits call frequency:       Last IAU call',             &
        src='setup_iau')
  END SELECT
  IF (IAU_QLimitsCallFreq /= CallFreq_Never) THEN
    WRITE(umMessage,'(A,L6)')                                                  &
       '- Write diagnostics?                   ',  L_IAU_QLimitsDiags
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,L6)')                                                  &
       '- Remove non-trop q increments?        ',  L_IAU_RmNonTropQIncs
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Lower limit to apply to trop RH:     ', IAU_trop_min_RH
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Upper limit to apply to non-trop q:  ', IAU_nonTrop_max_q
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Lower limit to apply to non-trop q:  ', IAU_nonTrop_min_q
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Upper limit to apply to non-trop RH: ', IAU_nonTrop_max_RH
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Minimum tropospheric pressure (hPa): ', IAU_trop_min_p    * 100.0
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Maximum tropospheric ABS(PV) (PVU):  ', IAU_trop_max_PV   * 1.0e06
    CALL umPrint(umMessage,src='setup_iau')
    WRITE(umMessage,'(A,ES11.4)')                                              &
      ' - Maximum non-trop     pressure (hPa): ', IAU_nonTrop_max_p * 100.0
    CALL umPrint(umMessage,src='setup_iau')
  END IF
END IF ! (PrintStatus >= PrStatus_Normal)

!-------------------------------------------------------------------------------
! [5]: Loop over increment files, setting up a data structure for each one.
!-------------------------------------------------------------------------------

! Timestep length in seconds:
TS_len_secs = secs_per_stepim(atmos_im)

! Filepath of first increment:
Inc1Path = TRIM( iau_inc_file )

! If there is more than one increment file, the paths of the second onwards
! will be obtained by overwriting the end of the first's filename with the
! increment number. Thus we need to check that the filename is long enough to
! be overwritten with the largest increment number.
IF (Num_IAU_incs > 1) THEN
  LenInc1Path = LEN_TRIM(Inc1Path)
  LenInc1Name = LenInc1Path - INDEX(Inc1Path, "/", back=.TRUE.)

  WRITE (NumIncsStr,'(I4)') Num_IAU_incs
  LenNumIncsStr = LEN_TRIM(ADJUSTL(NumIncsStr))

  IF (LenNumIncsStr > LenInc1Name) THEN
    WRITE (CMessage,'(A,I4,A)')                                     &
      'Length of name of first IAU increment file (', LenInc1Name,  &
      ') too short to be overwritten by all increment numbers.'
    ICode = 1
    CALL EReport (RoutineName, ICode, CMessage)
  END IF
END IF

ALLOCATE (IAU_incs(Num_IAU_incs))

! Loop over increment files:
DO IncNum = 1, Num_IAU_incs

  ! Get increment number as a string:
  WRITE (IncNumStr,'(I4)') IncNum
  IncNumStr = ADJUSTL(IncNumStr)

  ! Initialise flags:
  IAU_incs(IncNum) % FieldsToAdd           = .FALSE.
  IAU_incs(IncNum) % Contains_q_qCL_or_qCF = .FALSE.
  IAU_incs(IncNum) % Contains_qT           = .FALSE.
  IAU_incs(IncNum) % FieldsStored          = .FALSE.

  !-----------------------------------------------------------------------------
  ! [5.1]: Read fixed header and lookup headers.
  !-----------------------------------------------------------------------------

  ! Construct filepath:
  FilePath = Inc1Path
  IF (IncNum > 1) THEN
    ! Overwrite end of path with increment number:
    LenIncNumStr = LEN_TRIM(IncNumStr)
    FilePath(LenInc1Path-LenIncNumStr+1:LenInc1Path) = TRIM(IncNumStr)
  END IF

  IAU_incs(IncNum) % FilePath = FilePath

  ! Open file:
  CALL assign_file_unit (FilePath, IAU_unit, handler="portio")
  CALL model_file_open (IAU_unit, FilePath, read_write=ioOpenReadOnly)

  ! Read in fixed-length header:
  CALL read_flh (IAU_unit, IAU_incs(IncNum) % FixHd, Len_FixHd, &
                 ICode, tempmessage)
  IF (ICode > 0) THEN
    CMessage(1:55)   = 'Error reading fixed header from increment no.'// &
                       TRIM(IncNumStr)
    CMessage(56:80) = 'Message from READ_FLD:'
    Cmessage(81:errormessagelength) = tempmessage(1:errormessagelength-80)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  IAU_incs(IncNum) % Len1Lookup = IAU_incs(IncNum) % FixHd(151)
  IAU_incs(IncNum) % Len2Lookup = IAU_incs(IncNum) % FixHd(152)

  ! Read in lookup tables:
  IF (IAU_incs(IncNum) % Len2Lookup /= 0) THEN

    ALLOCATE (IAU_incs(IncNum) % Lookup(IAU_incs(IncNum) % Len1Lookup,  &
                                        IAU_incs(IncNum) % Len2Lookup))

    CALL setpos (IAU_unit, IAU_incs(IncNum) % FixHd(150)-1, ICode)
    IF (ICode > 0) THEN
      CMessage = 'SETPOS error moving to start of lookup tables for '// &
                 'increment no.'// TRIM(IncNumStr)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    Lookup_len = IAU_incs(IncNum) % Len1Lookup &
               * IAU_incs(IncNum) % Len2Lookup

    CALL buffin ( IAU_unit,                       &
                  IAU_incs(IncNum) % Lookup, &
                  Lookup_len,                     &
                  LenIO,                          &
                  a_io )
    IF ( (a_io  /= -1.0      ) .OR.   &
         (LenIO /= Lookup_len) ) THEN
      CALL ioerror ('Buffer in of lookups from IAU increment', &
                    a_io, LenIO, Lookup_len)
      ICode    = NINT(a_io) + 1
      CMessage = 'Error reading lookup tables from increment no.'// &
                 TRIM(IncNumStr)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

  END IF

  ! Close file:
  CALL model_file_close (IAU_unit, FilePath, delete=ioNoDelete, error=ICode)
  CALL release_file_unit(IAU_unit, handler="portio")

  IF (ICode > 0) THEN
    CMessage(1:80) = 'Error closing IAU increment no.'// TRIM(IncNumStr) //':'
    CMessage(81:)  = TRIM(FilePath)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Cycle if there are no fields in the increment:
  IF (IAU_incs(IncNum) % Len2Lookup == 0) CYCLE

  !-----------------------------------------------------------------------------
  ! [5.2]: Get length of array required to hold data, and determine whether
  !        the data is all 32-bit packed, in which case the data will be
  !        stored in packed form.
  !-----------------------------------------------------------------------------

  AllPacked = .TRUE.
  flds_len  = 0

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup

    Code = IAU_incs(IncNum) % Lookup(item_code, FieldNum)

    ! Get local field length:
    LocFldLen = 0
    DO i = 1, IAU_NumFldCodes
      IF (IAU_FldCodes(i) == Code) LocFldLen = IAU_LocFldLens(i)
    END DO

    flds_len = flds_len + LocFldLen

    pack_code = MOD(IAU_incs(IncNum) % Lookup(lbpack, FieldNum), 10)
    IF (pack_code /= PC_Cray32_Packing) AllPacked = .FALSE.

  END DO

  ! Cycle if there are no relevant fields in the file:
  IF (flds_len == 0) THEN
    DEALLOCATE (IAU_incs(IncNum) % Lookup)
    CYCLE
  ELSE
    IAU_incs(IncNum) % FieldsToAdd = .TRUE.
  END IF

  IAU_incs(IncNum) % StorePacked = AllPacked
  IAU_incs(IncNum) % flds_len    = flds_len

  !-----------------------------------------------------------------------------
  ! [5.3]: See whether there are increments to qT, q, qCL, qCF, dust1-6, SO2,
  !        ukca_O3 or ukca_NO2, and do the necessary checks.
  !-----------------------------------------------------------------------------

  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code  = IAU_incs(IncNum) % Lookup(item_code, FieldNum)
    Level = IAU_incs(IncNum) % Lookup(lblev,     FieldNum)
    IF (Level <= model_levels) THEN
      IF (Code == 10 .OR. Code == 12 .OR. Code == 254) THEN
        IAU_incs(IncNum) % Contains_q_qCL_or_qCF = .TRUE.
      ELSE IF (Code == 16207 .OR. Code == 18001) THEN
        IAU_incs(IncNum) % Contains_qT           = .TRUE.
      END IF
    END IF
  END DO

  ! Abort if file contains qT fields with both the new and old STASH codes:
  contains_new_qT_code = .FALSE.
  contains_old_qT_code = .FALSE.
  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code  = IAU_incs(IncNum) % Lookup(item_code, FieldNum)
    Level = IAU_incs(IncNum) % Lookup(lblev,     FieldNum)
    IF (Level <= model_levels) THEN
      IF (Code == 16207) contains_new_qT_code = .TRUE.
      IF (Code == 18001) contains_old_qT_code = .TRUE.
    END IF
  END DO
  IF (contains_new_qT_code .AND. contains_old_qT_code) THEN
    ICode = 1
    CMessage = 'Old qT code (18001) being used with new one (16207) in '// &
               'increment no.'// TRIM(IncNumStr)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Abort if the file contains the wrong number of dust bins for the model
  ! being run (2 bin or 6 bin)
  IAUNumDustBins = 0
  DO i = 1, 6
    IAUHasDustBin(i) = .FALSE.
  END DO
  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code = IAU_incs(IncNum) % Lookup(item_code, FieldNum)
    ! All dust bins have codes in this range:
    IF (Code >= 431 .AND. Code <= 436) THEN
      IAUDustBinNum = Code-430
      IF (.NOT. IAUHasDustBin(IAUDustBinNum)) THEN
        IAUHasDustBin(IAUDustBinNum) = .TRUE.
        IAUNumDustBins = IAUNumDustBins + 1
      END IF
    END IF
  END DO
  WRITE (IAUNumDustBinsStr,'(I9)') IAUNumDustBins
  IF (IAUNumDustBins > 0) THEN
    IF (l_twobin_dust .AND. IAUNumDustBins /= 2) THEN
      ICode = 1
      CMessage = 'Number of dust bins in IAU is ' // IAUNumDustBinsStr  // &
                 ', but the dust scheme is being run with 2 size bins.'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
    IF (.NOT. l_twobin_dust .AND. IAUNumDustBins /= 6) THEN
      ICode = 1
      CMessage = 'Number of dust bins in IAU is ' // IAUNumDustBinsStr  // &
                 ', but the dust scheme is being run with 6 size bins.'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF
  END IF

  ! See whether the file contains SO2, ukca_O3 or ukca_NO2 increments.
  DO FieldNum = 1, IAU_incs(IncNum) % Len2Lookup
    Code  = IAU_incs(IncNum) % Lookup(item_code, FieldNum)
    Level = IAU_incs(IncNum) % Lookup(lblev, FieldNum)
    IF (Level >= 1 .AND. Level <= model_levels) THEN
      IF (Code == stashcode_SO2)   Found_SO2      = .TRUE.
      IF (Code == stashcode_O3) Found_ukca_O3  = .TRUE.
      IF (Code == stashcode_NO2) Found_ukca_NO2 = .TRUE.
    END IF
  END DO

  ! Abort if the file contains SO2 increments but there are no SO2 tracer
  ! fields to add them to.
  IF (Found_SO2 .AND. .NOT. l_sulpc_so2) THEN
    ICode = 1
    CMessage = 'l_sulpc_so2=.false., but SO2 (section 0, item 101) fields '// &
               'found in increment no.'// TRIM(IncNumStr)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Abort if the file contains ukca_O3 or ukca_NO2 increments but UKCA is not
  ! active.
  IF (Found_ukca_O3 .AND. .NOT. l_ukca) THEN
    ICode = 1
    CMessage = 'l_ukca=.false., but ukca_O3 (section 34, item 1) fields '// &
               'found in increment no.'// TRIM(IncNumStr)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF
  IF (Found_ukca_NO2 .AND. .NOT. l_ukca) THEN
    ICode = 1
    CMessage = 'l_ukca=.false., but ukca_NO2 (section 34, item 4) fields '// &
               'found in increment no.'// TRIM(IncNumStr)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! If necessary, get indices of O3 and NO2 fields in the tracer_ukca array.
  IF (Found_ukca_O3)  ukca_O3_index  = ukca_name2index('O3')
  IF (Found_ukca_NO2) ukca_NO2_index = ukca_name2index('NO2')

  !-----------------------------------------------------------------------------
  ! [5.4]: Get list of usable fields.
  !-----------------------------------------------------------------------------

  IF (PrintStatus >= PrStatus_Normal) THEN
    UsableFlds(IncNum) = ''
    DO i = 1, IAU_NumFldCodes
      IF (ANY(IAU_incs(IncNum) % Lookup(item_code,:) == IAU_FldCodes(i))) &
        UsableFlds(IncNum) = TRIM(UsableFlds(IncNum)) //','//             &
                             TRIM(IAU_FldDescs(i))
    END DO
    UsableFlds(IncNum) = UsableFlds(IncNum)(2:)
  END IF

  !-----------------------------------------------------------------------------
  ! [5.5]: Calculate filter details.
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [5.5.1]: Get increment data and validity times relative to basis time.
  !-----------------------------------------------------------------------------

  CALL time2sec ( IAU_incs(IncNum) % FixHd(21), & ! in
                  IAU_incs(IncNum) % FixHd(22), & ! in
                  IAU_incs(IncNum) % FixHd(23), & ! in
                  IAU_incs(IncNum) % FixHd(24), & ! in
                  IAU_incs(IncNum) % FixHd(25), & ! in
                  0,                            & ! in
                  basis_time_days,              & ! in
                  basis_time_secs,              & ! in
                  ElapsedDays,                  & ! out
                  ElapsedSecs,                  & ! out
                  lcal360 )                       ! in

  DTSecs = ElapsedSecs + ElapsedDays * isec_per_day
  DTMins = DTSecs / 60

  CALL time2sec ( IAU_incs(IncNum) % FixHd(28), & ! in
                  IAU_incs(IncNum) % FixHd(29), & ! in
                  IAU_incs(IncNum) % FixHd(30), & ! in
                  IAU_incs(IncNum) % FixHd(31), & ! in
                  IAU_incs(IncNum) % FixHd(32), & ! in
                  0,                            & ! in
                  basis_time_days,              & ! in
                  basis_time_secs,              & ! in
                  ElapsedDays,                  & ! out
                  ElapsedSecs,                  & ! out
                  lcal360 )                       ! in

  VTSecs = ElapsedSecs + ElapsedDays * isec_per_day
  VTMins = VTSecs / 60

  IF (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter) THEN

    !---------------------------------------------------------------------------
    ! [5.5.2]: Filter details are specified via the IAU namelist.
    !---------------------------------------------------------------------------

    IAU_incs(1) % TendencyIncs = .FALSE.

    IF ( (VTMins < IAU_StartMin) .OR. &
         (VTMins > IAU_EndMin) ) THEN
      ICode = 1
      CMessage = 'Validity time of first IAU increment outside IAU period'
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    !---------------------------------------------------------------------------
    ! [5.5.2.1]: Calculate weights covering the whole filter span.
    !---------------------------------------------------------------------------

    NumWeights = (IAU_EndMin - IAU_StartMin)*60 / TS_len_secs + 1

    ALLOCATE (Weights(NumWeights))

    CALL Calc_TFiltWts ( IAU_StartMin*60,    & ! in
                         IAU_EndMin*60,      & ! in
                         TS_len_secs,        & ! in
                         IAU_ApexMin*60,     & ! in
                         IAU_Cutoff_period,  & ! in
                         IAU_SBE_period,     & ! in
                         IAU_FilterType,     & ! in
                         NumWeights,         & ! in
                         Weights )             ! out

    !---------------------------------------------------------------------------
    ! [5.5.2.2]: Eliminate zero weights at the start and end of the span.
    !---------------------------------------------------------------------------

    FirstNonZeroWeightIndex = 0
    LastNonZeroWeightIndex  = 0

    DO i = 1, NumWeights
      IF (Weights(i) /= 0.0) THEN
        IF (FirstNonZeroWeightIndex == 0) FirstNonZeroWeightIndex = i
        LastNonZeroWeightIndex = i
      END IF
    END DO

    NumWeights = LastNonZeroWeightIndex - FirstNonZeroWeightIndex + 1

    ALLOCATE (IAU_incs(1) % Weights(NumWeights))

    IAU_incs(1) % Weights(1:NumWeights) = Weights(FirstNonZeroWeightIndex: &
                                                  LastNonZeroWeightIndex)

    DEALLOCATE (Weights)

    IAU_incs(1) % InsertionStartTS = IAU_StartMin * 60 / TS_len_secs &
                                   + FirstNonZeroWeightIndex - 1
    IAU_incs(1) % InsertionEndTS   = IAU_StartMin * 60 / TS_len_secs &
                                   + LastNonZeroWeightIndex  - 1

  ELSE

    !---------------------------------------------------------------------------
    ! [5.5.3]: Filter details to be deduced from information written to the
    !          fixed header by VAR.
    !
    ! The key parameter is the "TimeID", supplied as FixHd(10).
    !
    ! There are two main possibilities:
    !
    ! 1. TimeID = TimeID_Once
    !    Add the whole increment at the UM timestep nearest to the validity
    !    time supplied in the fixed header.
    !
    ! 2. TimeID /= TimeID_Once
    !    Treat the increments as (per second) tendencies for a time window
    !    starting at the validity time and ending at the data time supplied in
    !    the fixed header; the period over which the tendency will have been
    !    applied within VAR.
    !
    !    The start and end points of this window are adjusted so that they
    !    coincide with UM timesteps, and a correction factor calculated to
    !    account for the change in window length.
    !
    !    The shape of the weighting function to be used over the window is
    !    determined by the specific value of TimeID. Currently the only choice
    !    is uniform weights, but further choices may be added in the future.
    !---------------------------------------------------------------------------

    ! TimeID supplied by VAR:
    TimeID = IAU_incs(IncNum) % FixHd(10)

    IF (TimeID /= TimeID_Once .AND. &
        TimeID /= TimeID_Constant) THEN
      ICode = 1
      WRITE (CMessage,'(A,I6,A)') &
        'Unsupported TimeID (', TimeID, ') for increment no.'// TRIM(IncNumStr)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    VAR_StartSec = VTSecs
    VAR_EndSec   = DTSecs

    IF (VAR_EndSec < VAR_StartSec) THEN
      ICode = 1
      CMessage = 'VAR end time before start time for increment no.'// &
                 TRIM(IncNumStr)
      CALL EReport (RoutineName, ICode, CMessage)
    END IF

    ! Adjust start and end times to coincide with UM timesteps:
    UM_StartSec = TS_len_secs * (VAR_StartSec / TS_len_secs)
    UM_EndSec   = TS_len_secs * (VAR_EndSec   / TS_len_secs)
    IF ((VAR_StartSec - UM_StartSec) * 2 > TS_len_secs) THEN
      UM_StartSec = UM_StartSec + TS_len_secs
    END IF
    IF ((VAR_EndSec   - UM_EndSec)   * 2 > TS_len_secs) THEN
      UM_EndSec   = UM_EndSec   + TS_len_secs
    END IF

    ! Save start and end times for printout step:
    IF (PrintStatus >= PrStatus_Normal) THEN
      VAR_StartSecSave(IncNum) = VAR_StartSec
      VAR_EndSecSave  (IncNum) = VAR_EndSec
      UM_StartSecSave (IncNum) = UM_StartSec
      UM_EndSecSave   (IncNum) = UM_EndSec
    END IF

    IF (TimeID == TimeID_Once) THEN

      IAU_incs(IncNum) % TendencyIncs = .FALSE.

      IF (VAR_StartSec /= VAR_EndSec) THEN
        ICode = 1
        CMessage = 'TimeID=TimeID_Once, but VT and DT do not match for '// &
                   'increment no.'// TRIM(IncNumStr)
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      ALLOCATE (IAU_incs(IncNum) % Weights(1))

      IAU_incs(IncNum) % Weights(1)       = 1.0
      IAU_incs(IncNum) % InsertionStartTS = UM_StartSec / TS_len_secs
      IAU_incs(IncNum) % InsertionEndTS   = UM_EndSec   / TS_len_secs

    ELSE

      IAU_incs(IncNum) % TendencyIncs = .TRUE.

      IF (UM_StartSec == UM_EndSec) THEN
        ICode = 1
        CMessage = 'TimeID/=TimeID_Once, but zero period in UM for '// &
                   'increment no.'// TRIM(IncNumStr)
        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      ! Tendencies applied at the end of each timestep within the window:
      IAU_incs(IncNum) % InsertionStartTS = (UM_StartSec / TS_len_secs) + 1
      IAU_incs(IncNum) % InsertionEndTS   =  UM_EndSec   / TS_len_secs

      NumWeights = IAU_incs(IncNum) % InsertionEndTS   &
                 - IAU_incs(IncNum) % InsertionStartTS &
                 + 1

      ALLOCATE (IAU_incs(IncNum) % Weights(NumWeights))

      ! Calculate filter weights:
      SELECT CASE (TimeID)
      CASE (TimeID_Constant)
        IAU_incs(IncNum) % Weights(:) = 1.0
      END SELECT

      ! Scale weights to convert from (per second) tendencies to
      ! instantaneous increments. (This includes a correction factor to
      ! account for the change of window moving from VAR to the UM.)
      IAU_incs(IncNum) % Weights(:) = IAU_incs(IncNum) % Weights(:)   &
                                    * REAL(TS_len_secs)               &
                                    * REAL(VAR_EndSec - VAR_StartSec) &
                                    / REAL(UM_EndSec  - UM_StartSec)

    END IF ! (TimeID == TimeID_Once)

  END IF ! (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter)

  !-----------------------------------------------------------------------------
  ! [5.5.4]: Check that insertion period does not start before basis time.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % InsertionStartTS < 0) THEN
    ICode = 1
    CMessage = 'Insertion period starts before basis time basis time for '// &
               'increment no.'// TRIM(IncNumStr)
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  !-----------------------------------------------------------------------------
  ! [5.6]: Set StoreFields flag, indicating whether the increments should be
  !        stored in memory during the increment's insertion period.
  !-----------------------------------------------------------------------------

  IF (IAU_incs(IncNum) % InsertionEndTS   > &
      IAU_incs(IncNum) % InsertionStartTS) THEN
    IAU_incs(IncNum) % StoreFields = .TRUE.
  ELSE
    IAU_incs(IncNum) % StoreFields = .FALSE.
  END IF

  !-----------------------------------------------------------------------------
  ! [5.7]: Update timestep numbers for first and last IAU call.
  !-----------------------------------------------------------------------------

  IF (IAU_FirstCallTS == -1) THEN
    IAU_FirstCallTS = IAU_incs(IncNum) % InsertionStartTS
    IAU_LastCallTS  = IAU_incs(IncNum) % InsertionEndTS
  ELSE
    IAU_FirstCallTS = MIN(IAU_FirstCallTS, &
                          IAU_incs(IncNum) % InsertionStartTS)
    IAU_LastCallTS  = MAX(IAU_LastCallTS,  &
                          IAU_incs(IncNum) % InsertionEndTS)
  END IF

END DO ! (IncNum)

!-------------------------------------------------------------------------------
! [6]: Apply some restrictions to the way humidity and cloud are updated.
!
!      There are currently two methods for obtaining humidity and cloud
!      increments:
!        1. Diagnose them from qT increments in the IAU files.
!        2. Update them directly from the q, qCL and/or qCF increments in
!           the IAU files, and then apply further modifications if using the
!           PC2 cloud scheme.
!      For simplicity, we currently demand that only one of these methods is
!      used on an individual IAU timestep.
!-------------------------------------------------------------------------------

! Check that there are no increment files containing q, qCL or qCF if they
! also contain qT:
DO IncNum = 1, Num_IAU_incs

  IF (IAU_incs(IncNum) % Contains_q_qCL_or_qCF .AND. &
      IAU_incs(IncNum) % Contains_qT) THEN
    ICode = 1
    WRITE (CMessage(1:80),'(A,I2.2,A)') &
      'Increment no.', IncNum, ' contains q, qCL or qCF, but also qT.'
    CMessage(81:) = 'This combination is currently unsupported.'
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

END DO

! Check that q, qCL or qCF increments will not be used together with qT
! increments on the same IAU timestep:
DO TS_num = IAU_FirstCallTS, IAU_LastCallTS

  Using_q_qCL_or_qCF = .FALSE.
  Using_qT           = .FALSE.

  DO IncNum = 1, Num_IAU_incs
    IF (IAU_incs(IncNum) % FieldsToAdd) THEN
      IF (TS_num >= IAU_incs(IncNum) % InsertionStartTS .AND. &
          TS_num <= IAU_incs(IncNum) % InsertionEndTS) THEN
        IF (IAU_incs(IncNum) % Contains_q_qCL_or_qCF) &
          Using_q_qCL_or_qCF = .TRUE.
        IF (IAU_incs(IncNum) % Contains_qT)           &
          Using_qT           = .TRUE.
      END IF
    END IF
  END DO

  IF (Using_q_qCL_or_qCF .AND. Using_qT) THEN
    ICode = 1
    WRITE(CMessage(1:80),'(A,I6)')                                        &
      'qT incs being used together with q, qCL or qCF incs on timestep ', &
      TS_num
    CMessage(81:) = 'This combination is currently unsupported.'
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

END DO

! Print warning if asking to diagnose ice cloud increments from qT
! increments, but no qT increments available:
IF (L_IAU_IncrementIce .AND. .NOT. ANY(IAU_incs(:) % Contains_qT)) THEN
  ICode = -1
  CMessage = 'L_IAU_IncrementIce set, but no qT increments available'
  CALL EReport (RoutineName, ICode, CMessage)
END IF

!-------------------------------------------------------------------------------
! [7]: Print details for each increment.
!-------------------------------------------------------------------------------

IF (PrintStatus >= PrStatus_Normal) THEN

  CALL umPrint( '',src='setup_iau')
  CALL umPrint( '(Times given below are relative to the model basis time.)', &
      src='setup_iau')
  CALL umPrint( '',src='setup_iau')

  DO IncNum = 1, Num_IAU_incs

    ! Get increment number as a string:
    WRITE (IncNumStr,'(I4)') IncNum
    IncNumStr = ADJUSTL(IncNumStr)

    CALL umPrint( '** Details for increment no.'// TRIM(IncNumStr) //':', &
        src='Setup_IAU')
    IF (.NOT. IAU_incs(IncNum) % FieldsToAdd) THEN
      CALL umPrint( '   No usable fields',src='setup_iau')
      CYCLE
    ELSE
      CALL umPrint( '   Usable fields: '// TRIM(UsableFlds(IncNum)), &
          src='setup_iau')
    END IF

    IF (IncNum == 1 .AND. L_IAU_SpecifyInc1Filter) THEN

      IF (IAU_FilterType == 1) THEN
        CALL umPrint( ' Filter type: Uniform',src='setup_iau')
      ELSE IF (IAU_FilterType == 2) THEN
        CALL umPrint( ' Filter type: Triangular',src='setup_iau')
        WRITE(umMessage,'(A,I8)') ' Apex minute:           ', IAU_ApexMin
        CALL umPrint(umMessage,src='setup_iau')
      ELSE IF (IAU_FilterType == 3) THEN
        CALL umPrint( ' Filter type: Lanczos windowed',src='setup_iau')
        WRITE(umMessage,'(A,ES11.4,A)') ' Cut-off period: ', &
            IAU_Cutoff_period, ' hrs'
        CALL umPrint(umMessage,src='setup_iau')
      ELSE IF (IAU_FilterType == 4) THEN
        CALL umPrint( ' Filter type: Dolph',src='setup_iau')
        WRITE(umMessage,'(A,ES11.4,A)')                                        &
                             ' Stop band edge period: ', IAU_SBE_period, ' hrs'
        CALL umPrint(umMessage,src='setup_iau')
      END IF
      WRITE(umMessage,'(A,I8)') ' Insertion start min:   ', IAU_StartMin
      CALL umPrint(umMessage,src='setup_iau')
      WRITE(umMessage,'(A,I8)') ' Insertion end   min:   ', IAU_EndMin
      CALL umPrint(umMessage,src='setup_iau')

    ELSE

      VAR_StartSec = VAR_StartSecSave(IncNum)
      VAR_EndSec   = VAR_EndSecSave  (IncNum)
      UM_StartSec  = UM_StartSecSave (IncNum)
      UM_EndSec    = UM_EndSecSave   (IncNum)

      ! Time strings for start and end of time windows:
      WRITE (VAR_StartStr,'(I3,":",I2.2,":",I2.2)') &
        VAR_StartSec/3600, MOD(VAR_StartSec/60, 60), MOD(VAR_StartSec, 60)
      WRITE (VAR_EndStr,  '(I3,":",I2.2,":",I2.2)') &
        VAR_EndSec  /3600, MOD(VAR_EndSec  /60, 60), MOD(VAR_EndSec,   60)
      WRITE (UM_StartStr, '(I3,":",I2.2,":",I2.2)') &
        UM_StartSec /3600, MOD(UM_StartSec /60, 60), MOD(UM_StartSec,  60)
      WRITE (UM_EndStr,   '(I3,":",I2.2,":",I2.2)') &
        UM_EndSec   /3600, MOD(UM_EndSec   /60, 60), MOD(UM_EndSec,    60)

      WRITE(umMessage,'(A,I8)') '   TimeID:       ',                           &
          IAU_incs(IncNum) % FixHd(10)
      CALL umPrint(umMessage,src='setup_iau')
      WRITE(umMessage,'(A,L6)') '   Tendency?:    ',                           &
          IAU_incs(IncNum) % TendencyIncs
      CALL umPrint(umMessage,src='setup_iau')
      IF (IAU_incs(IncNum) % TendencyIncs) THEN
        CALL umPrint( '   VAR window (hr:mm:ss): '//                           &
                    VAR_StartStr// ' to '// VAR_EndStr,                        &
                    src='setup_iau')
        CALL umPrint( '   UM  window (hr:mm:ss): '//                           &
                    UM_StartStr //  ' to '// UM_EndStr,                        &
                    src='setup_iau')
      ELSE
        CALL umPrint( '   Increment VT   (hr:mm:ss): '// VAR_StartStr,         &
            src='setup_iau')
        CALL umPrint( '   Insertion time (hr:mm:ss): '// UM_StartStr,          &
            src='setup_iau')
      END IF

    END IF

    NumWeights = SIZE(IAU_incs(IncNum) % Weights)

    IF (NumWeights > 1) THEN
      CALL umPrint( '   IAU Weights:  hr:mm:ss  Weight',src='setup_iau')
      DO WeightNum = 1, NumWeights
        Sec = (IAU_incs(IncNum) % InsertionStartTS + WeightNum - 1) &
            * TS_len_secs
        WRITE (TimeStr,'(I3,":",I2.2,":",I2.2)')  &
          Sec/3600, MOD(Sec/60, 60), MOD(Sec, 60)
        WRITE(umMessage,'(A,A,ES18.10)')                                       &
          '                 ', TimeStr, IAU_incs(IncNum) % Weights(WeightNum)
        CALL umPrint(umMessage,src='setup_iau')
      END DO
    END IF

  END DO ! (IncNum)

  CALL umPrint( '',src='setup_iau')
  CALL umPrint( 'End of Incremental Analysis Update (IAU) details',            &
      src='setup_iau')
  CALL umPrint( '<><><><><><><><><><><><><><><><><><><><><><><><>',            &
      src='setup_iau')
  CALL umPrint( '',src='setup_iau')

END IF ! (PrintStatus >= PrStatus_Normal)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE Setup_IAU

END MODULE setup_iau_mod
