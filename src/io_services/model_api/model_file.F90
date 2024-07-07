! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A simple data module containing variables related to STASH buffering
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C84
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

MODULE model_file

USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim
USE ios, ONLY:                                                              &
    IOS_Init_PP_Lookup,                                                      &
    IOS_Set_Header,                                                          &
    IOS_Init_Header
USE IOS_Common, ONLY: minUnit, maxUnit                                   
USE IOS_Communicators, ONLY:                                                &
    model_rank, model_procs
USE IOS_Stash_Common, ONLY:                                                 &
    L_IO_Server,                                                             &
    isUsingAsyncStash
USE io, ONLY:                                                                &
    setpos, isRemote, returnFilename, readableWords, is_unit_open,           &
    isReadOnly, file_open, file_close, buffin
USE io_constants, ONLY: ioFileTypeUM
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore,  ONLY: mype
USE umPrintMgr
USE lookup_addresses
USE filenamelength_mod, ONLY: fileNameLength
USE file_manager, ONLY: get_file_by_unit, um_file_type, init_file_loop
USE MPPIO_file_utils, ONLY: file_touch, file_delete

USE io_configuration_mod, ONLY: io_external_control

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
PRIVATE
! Profiling
INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Locations in the pp file for interesting items, otherwise
! see f3.ps from the umdp
!
! File addresses
INTEGER, PARAMETER, PUBLIC   :: MF_Lookup_Address=150
INTEGER, PARAMETER, PUBLIC   :: MF_Data_Address=160
! Other Data
INTEGER, PARAMETER, PUBLIC   :: MF_Num_Preallocated_Headers=152
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CHARACTER (LEN=errormessagelength),PRIVATE  :: Cmessage ! General purpose string

! This value is part of the file format, and other utilities expect it
INTEGER, PARAMETER, PUBLIC   :: MF_Data_Missing     =-99

! Parts of the file
INTEGER, PARAMETER, PUBLIC   :: FixedHeader         = 1
INTEGER, PARAMETER, PUBLIC   :: FixedHeaderLen      = 256
INTEGER, PARAMETER, PUBLIC   :: IntConsts           = 2
INTEGER, PARAMETER, PUBLIC   :: RealConsts          = 3
INTEGER, PARAMETER, PUBLIC   :: RowConsts           = 4
INTEGER, PARAMETER, PUBLIC   :: LevConsts           = 5
INTEGER, PARAMETER, PUBLIC   :: ColConsts           = 6

! States of parts of the file
INTEGER, PARAMETER           :: stUnset             = 1
INTEGER, PARAMETER           :: stUpdFromModel      = 2
INTEGER, PARAMETER           :: stUpdFromIOS        = 3
INTEGER, PARAMETER           :: stPendingFromIOS    = 4

PUBLIC model_file_open
PUBLIC model_file_close
PUBLIC get_file_address
PUBLIC get_mf_information
PUBLIC InitHeader
PUBLIC attachHeader
PUBLIC SetHeader
PUBLIC loadHeader
PUBLIC storeHeader
PUBLIC InitLookups
PUBLIC attachLookups
PUBLIC SetLookups
PUBLIC StoreLookups
PUBLIC StoreAllLookups
PUBLIC SynchroniseAll
PUBLIC setRecordDiskLength
PUBLIC setRecordDataLength
PUBLIC setRecordDiskStart

! Semantics:
!
! Init           Initialise memory and reset, relays to IOS
! attach         Provide a pointer to the memory
! Set            Inform that data has changed, possibly providing data
!                relays to IOS
! Get            Recover data to user buffer
! load           load from disk
! store          save to disk

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MODEL_FILE'

CONTAINS


! Determine whether an operation has a relay to IO server action associated
LOGICAL FUNCTION shouldPropogate(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN)         :: UNIT
shouldPropogate=.FALSE.

IF (model_rank==0 .AND. isRemote(UNIT)) shouldPropogate=.TRUE.

END FUNCTION shouldPropogate


SUBROUTINE mf_printstate(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN)         :: UNIT
TYPE(um_file_type), POINTER :: um_file
    
NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")

!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I4,A)')'Managed unit ',UNIT,' *******************'
!$OMP END CRITICAL(internal_write)
CALL umPrint(umMessage,src='model_file')
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I8)')' lookup_address = ',                             &
    um_file % pp_meta % lookup_address
!$OMP END CRITICAL(internal_write)
CALL umPrint(umMessage,src='model_file')
IF (ASSOCIATED(um_file % pp_meta % fixed_header)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I8)')' fixed_header(size)',                          &
      SIZE(um_file % pp_meta % fixed_header)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='model_file')
ELSE
  CALL umPrint(' fixed_header not associated')
END IF
IF (ASSOCIATED(um_file % pp_meta % lookup_table)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I8)')' lookup_table(size)',                          &
      SIZE(um_file % pp_meta % lookup_table)
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='model_file')
ELSE
  CALL umPrint(' lookup_table not associated')
END IF

END SUBROUTINE mf_printstate

SUBROUTINE mf_ereport(r,code,m,UNIT)
USE io, ONLY: io_ereport
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: r
CHARACTER(LEN=*),INTENT(IN) :: m
INTEGER, INTENT(IN)         :: code ! code to pass through to ereport
INTEGER, INTENT(IN),                                                       &
    OPTIONAL                :: UNIT ! unit responsible for the error
INTEGER                     :: lcode
INTEGER                     :: i
REAL(KIND=jprb)             :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MF_EREPORT'
TYPE(um_file_type), POINTER :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
lcode=code ! fix ereport_mod/ereport API mismatch
WRITE(umMessage,'(A,A)')                                                           &
    '****************** Model_File Error Report ',                         &
    '***************************'
CALL umPrint(umMessage,src='model_file')
IF (PRESENT(UNIT)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(umMessage,'(A,I5)')'Unit Generating error=',UNIT
!$OMP END CRITICAL(internal_write)
  CALL umPrint(umMessage,src='model_file')
  !      IF (printstatus>=prstatus_diag.OR.code>0) THEN
  CALL mf_printstate(UNIT)

  NULLIFY(um_file)
  um_file => init_file_loop(handler="portio")
  DO WHILE (ASSOCIATED(um_file))
    IF (um_file % pp_meta % managed) THEN
      CALL mf_printstate(um_file % UNIT)
    END IF
    um_file => um_file % next
  END DO

  !      ELSE
  !        WRITE(6,'(A,A)')&
  !            '*** File states can be reported by ',&
  !            'setting diagnostic output levels **'
  !      END IF
END IF
CALL io_ereport(r,lcode,m,UNIT)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE mf_ereport


!----------------------------------------------------------------------
! Function: get_file_address
! Get the address of an object for a pp file
!----------------------------------------------------------------------
INTEGER FUNCTION get_file_address(UNIT,item)

IMPLICIT NONE
INTEGER,INTENT(IN) :: UNIT
INTEGER,INTENT(IN) :: item
CHARACTER (LEN=*), PARAMETER :: SubroutineName='MODEL_FILE:GET_FILE_ADDRESS'
TYPE(um_file_type), POINTER :: um_file

CALL checkUnit(UNIT,SubroutineName)

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")

IF (.NOT. ASSOCIATED(um_file % pp_meta % fixed_header)) THEN
  CALL mf_ereport(SubroutineName, 99,                                      &
      'Request for pp information before initalisation',UNIT)
END IF

SELECT CASE(item)
CASE (MF_Lookup_Address,MF_Data_Address)
  get_file_address = um_file % pp_meta % fixed_header(item)-1
CASE DEFAULT
  CALL mf_ereport(SubroutineName, 99,                                      &
      'Request for pp file address of an unknown type',UNIT)
END SELECT

END FUNCTION get_file_address


!----------------------------------------------------------------------
! Function: get_mf_information
! Provide information about a pp file
!----------------------------------------------------------------------
INTEGER FUNCTION get_mf_information(UNIT,item)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(IN) :: item
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:GET_MF_INFORMATION'
TYPE(um_file_type), POINTER :: um_file

CALL checkUnit(UNIT,SubroutineName)

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
IF (.NOT. ASSOCIATED(um_file % pp_meta % fixed_header)) THEN

  CALL mf_ereport(SubroutineName, 99,                                      &
      'Request for pp information before initalisation',UNIT)
END IF

IF ( item == MF_Num_Preallocated_Headers ) THEN
  get_mf_information = um_file % pp_meta % fixed_header(item)
ELSE

  CALL mf_ereport(SubroutineName, 99,                                      &
      'Request for pp information of an unknown type',UNIT)
END IF

END FUNCTION get_mf_information


!----------------------------------------------------------------------
! Subroutine: CheckUnit
! Sanity: ensure the unit is in the valid range
!----------------------------------------------------------------------
SUBROUTINE CheckUnit(UNIT,src)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
CHARACTER(LEN=*), INTENT(IN) :: src
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECKUNIT'
INTEGER                      :: errCode
CHARACTER (LEN=*), PARAMETER :: SubroutineName='MODEL_FILE:CHECKUNIT'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( UNIT < minUnit .OR. UNIT > maxUnit ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(Cmessage,'(A,i10,A,i3,A,i3,A)')'Unit ',UNIT,                       &
      ' supplied to a PP op was outside the allowed range [',              &
      minUnit,'-',maxUnit,']'
!$OMP END CRITICAL(internal_write)
  CALL mf_ereport(src//' via '//SubroutineName, 99,Cmessage,UNIT)
END IF

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
IF (.NOT. um_file % pp_meta % managed  ) THEN
  errCode=99
!$OMP CRITICAL(internal_write)
  WRITE(Cmessage,'(A,i4,A)')'Unit ',UNIT,                                  &
      ' supplied to model_file is not managed by model_file'
!$OMP END  CRITICAL(internal_write)
  CALL mf_ereport(src//' via '//SubroutineName, errCode,Cmessage,UNIT)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE CheckUnit

!----------------------------------------------------------------------
! Subroutine InitHeader
! Initialises fixed header structures
!----------------------------------------------------------------------
SUBROUTINE InitHeader(UNIT, header, headerLength)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, INTENT(IN)          :: header
INTEGER, INTENT(IN), OPTIONAL:: headerLength
INTEGER                      :: local_len
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITHEADER'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:INITHEADER'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL checkUnit(UNIT,SubroutineName)

IF (PRESENT(headerLength)) THEN
  local_len=headerLength
  IF (header==FixedHeader) THEN
    CALL umPrint('Cannot mandate the size of Fixed Header in UM files')
    local_len=fixedHeaderLen
  END IF
ELSE
  local_len=fixedHeaderLen
  IF (header/=FixedHeader) THEN
    CALL mf_ereport(SubroutineName,99,                                     &
        'Must specify the header length via the headLength argument')
  END IF
END IF

IF (shouldPropogate(UNIT)) CALL IOS_Init_Header(UNIT,local_len)

! Rank 0 in the atmos model needs to cause initialisation of the header
! on the remote server

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
ALLOCATE(um_file % pp_meta % fixed_header(local_len))
um_file % pp_meta % fixed_header(:) = MF_Data_Missing

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE InitHeader


!----------------------------------------------------------------------
! Subroutine SetHeader
! Sets the content of a fixed header from an input array
!----------------------------------------------------------------------
SUBROUTINE SetHeader(UNIT,header,data_in)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, INTENT(IN)          :: header
INTEGER, OPTIONAL            :: data_in(:)
INTEGER, POINTER             :: fixhd_internal(:)
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETHEADER'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:SETHEADER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
NULLIFY(fixhd_internal)
CALL checkUnit(UNIT,SubroutineName)
fixhd_internal=>attachHeader(UNIT,FixedHeader)

IF (PRESENT(data_in)) THEN

  IF ( SIZE(fixhd_internal) /= SIZE(data_in) ) THEN
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I3,A,I12,A,I12)')'ERROR: UNIT=',UNIT,              &
        ' LEN INTERNAL=',                                                  &
        SIZE(fixhd_internal),' LEN INPUT=',                                &
        SIZE(data_in)
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='model_file')

    CALL mf_ereport(SubroutineName, 99, 'bad parameters in setHeader',UNIT)
  END IF

  fixhd_internal(:)=data_in(:)

END IF

CALL IOS_Set_Header(UNIT,fixhd_internal)

NULLIFY(fixhd_internal)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE SetHeader

SUBROUTINE loadHeader(UNIT,header,data_out,icode)

IMPLICIT NONE
INTEGER, INTENT(IN)            :: UNIT
INTEGER, INTENT(IN)            :: header
INTEGER, OPTIONAL, INTENT(OUT) :: data_out(:)
INTEGER, OPTIONAL              :: icode
INTEGER, POINTER               :: fixhd_internal(:)
INTEGER                        :: mver_maj
INTEGER                        :: mver_min
REAL(KIND=jprb)                :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LOADHEADER'
CHARACTER (LEN=*), PARAMETER   :: SubroutineName=                          &
    'MODEL_FILE:LOADHEADER'
TYPE(um_file_type), POINTER    :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (PRESENT(icode)) icode=0

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
IF (.NOT. ASSOCIATED (um_file % pp_meta % fixed_header))                   &
    CALL initHeader(um_file % UNIT,header)

CALL setpos(um_file % UNIT,0)

IF (PRESENT(icode)) THEN
  icode=0
  IF (readableWords(um_file % UNIT) < &
      SIZE(um_file % pp_meta % fixed_header)) THEN
    icode=-1
!$OMP CRITICAL(internal_write)
    WRITE(umMessage,'(A,I3,A)')'loadHeader: The file on unit ',            &
        um_file % UNIT, ' does not have enough content to load the header'
!$OMP END CRITICAL(internal_write)
    CALL umPrint(umMessage,src='model_file')
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

CALL buffin(um_file % UNIT, um_file % pp_meta % fixed_header)

mver_maj = um_file % pp_meta % fixed_header(12) / 100
mver_min = um_file % pp_meta % fixed_header(12) - 100 * mver_maj
!$OMP CRITICAL(internal_write)
WRITE(umMessage,'(A,I0,A,I0)')'loadHeader: Model Version: ',               &
    mver_maj,'.',mver_min
!$OMP END CRITICAL(internal_write)
CALL umPrint(umMessage,src='model_file')

IF (PRESENT(data_out)) THEN
  data_out(:) = um_file % pp_meta % fixed_header(:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE loadHeader

!----------------------------------------------------------------------
! SUBROUTINE storeHeader
! Not yet implemented
!----------------------------------------------------------------------

SUBROUTINE storeHeader(UNIT,header)

IMPLICIT NONE
INTEGER, INTENT(IN)            :: UNIT
INTEGER, INTENT(IN)            :: header
REAL(KIND=jprb)                :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STOREHEADER'
CHARACTER (LEN=*), PARAMETER   :: SubroutineName=                          &
    'MODEL_FILE:STOREHEADER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL mf_ereport(SubroutineName,99,'Not Implemented',UNIT)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE storeHeader

!----------------------------------------------------------------------
! FUNCTION attachHeader
! Attaches a header to a given unit
!----------------------------------------------------------------------
FUNCTION attachHeader(UNIT,header) RESULT(hdr)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, INTENT(IN)          :: header
INTEGER, POINTER             :: hdr(:)
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATTACHHEADER'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:ATTACHHEADER'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
NULLIFY(hdr)
NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
CALL checkUnit(um_file % UNIT,SubroutineName)

IF (.NOT. ASSOCIATED(um_file % pp_meta % fixed_header))                     &
    CALL mf_ereport(SubroutineName,99,                                     &
    'Request for fixed header which is not allocated yet', um_file % UNIT)

hdr => um_file % pp_meta % fixed_header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION attachHeader


!----------------------------------------------------------------------
! SUBROUTINE InitLookups
! Initialises lookup table structures
!----------------------------------------------------------------------
SUBROUTINE InitLookups(ipplook, UNIT, len1,len2, address, step)

IMPLICIT NONE
INTEGER, POINTER                :: ipplook(:,:)
INTEGER, INTENT(IN)             :: UNIT
INTEGER, INTENT(IN)             :: step
INTEGER, INTENT(IN)             :: len1
INTEGER, INTENT(IN)             :: len2
INTEGER, INTENT(IN)             :: address
INTEGER                         :: ii, jj
INTEGER                         :: icode
REAL(KIND=jprb)                 :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITLOOKUPS'
INTEGER                         :: input_args(5)
INTEGER, PARAMETER              :: current_io_pe = 0
CHARACTER (LEN=*), PARAMETER    :: SubroutineName=                         &
    'MODEL_FILE:INITLOOKUPS'
TYPE(um_file_type), POINTER     :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(ipplook)
NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
CALL checkUnit(um_file % UNIT,SubroutineName)

! Rank 0 in the atmos model needs to cause initialisation of the lookup
! table on the remote server
IF ( shouldPropogate(um_file % UNIT) ) THEN
  input_args( 1 )=um_file % UNIT
  input_args( 2 )=len1
  input_args( 3 )=len2
  input_args( 4 )=address
  input_args( 5 )=step
  CALL IOS_Init_PP_Lookup(input_args)
END IF

IF ( step == 1 ) THEN

  IF (ASSOCIATED(um_file % pp_meta % lookup_table))                        &
      DEALLOCATE(um_file % pp_meta % lookup_table)

  IF ( mype == current_io_pe .OR. l_io_server ) THEN
    ALLOCATE(um_file % pp_meta % lookup_table(len1,len2))
  ELSE
    ! Save a bit of memory, we don't need lookups on
    ! general ranks - but we do want the table associated
    ALLOCATE(um_file % pp_meta % lookup_table(1,1))
  END IF

  ipplook => um_file % pp_meta % lookup_table
  ipplook(:,:) = MF_Data_Missing

ELSE IF ( step == 2 ) THEN

  IF ( mype == current_io_pe .OR. l_io_server )                            &
      um_file % pp_meta % lookup_address = address

ELSE   ! step /= 1 or 2
  Cmessage = 'InitLookups step should be 1 or 2'
  Icode = 10

  CALL mf_ereport(SubroutineName, Icode, Cmessage, um_file % UNIT)
END IF ! Test on step

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE InitLookups


!----------------------------------------------------------------------
! FUNCTION attachLookups
!
! This routine attaches a buffered lookup table to a pointer
! If we are using an IO server with A I O we send a corresponding
! instruction to the IO server too
!----------------------------------------------------------------------
FUNCTION attachLookups(UNIT) RESULT(lookup)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, POINTER             :: lookup(:,:)
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATTACHLOOKUPS'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:ATTACHLOOKUPS'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
CALL checkUnit(um_file % UNIT,SubroutineName)

NULLIFY (lookup)

IF (ASSOCIATED(um_file % pp_meta % lookup_table)) THEN
  lookup => um_file % pp_meta % lookup_table
ELSE
  CALL umPrint(                                                            &
      'DEBUG: attachLookups: doesnt exist for unit ', um_file % UNIT)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION attachLookups

!----------------------------------------------------------------------
! Subroutine SetRecordDiskLength
! Set the length of a field on the lookup table
! (meaning the disk length, the bytes written to disk including any
!  padding)
!----------------------------------------------------------------------
SUBROUTINE setRecordDiskLength(UNIT,record,SIZE)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(IN) :: record
INTEGER, INTENT(IN) :: SIZE
INTEGER, POINTER    :: lookup(:,:) => NULL()
lookup=>attachLookups(UNIT)
lookup(lbnrec,record )=SIZE
NULLIFY(lookup)
END SUBROUTINE setRecordDiskLength

!----------------------------------------------------------------------
! Subroutine SetRecordDataLength
! Set the length of a field on the lookup table
! (meaning the data length, the bytes written to disk excluding any
!  padding)
!----------------------------------------------------------------------
SUBROUTINE setRecordDataLength(UNIT,record,SIZE)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(IN) :: record
INTEGER, INTENT(IN) :: SIZE
INTEGER, POINTER    :: lookup(:,:) => NULL()
lookup=>attachLookups(UNIT)
lookup(lblrec,record )=SIZE
NULLIFY(lookup)

END SUBROUTINE setRecordDataLength

!----------------------------------------------------------------------
! Subroutine SetDiskRecordStart
! Set the start address of a field in the lookup table
! (meaning the data start address of a field)
!----------------------------------------------------------------------
SUBROUTINE setRecordDiskStart(UNIT,record,location)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER, INTENT(IN) :: record
INTEGER, INTENT(IN) :: Location
INTEGER, POINTER    :: lookup(:,:) => NULL()

lookup=>attachLookups(UNIT)
lookup(lbegin,record )=location
NULLIFY(lookup)

END SUBROUTINE setRecordDiskStart

!----------------------------------------------------------------------
! Subroutine Merge Lookup
! This routine merges a supplied lookup table into a stored one
! We overwrite all records with the supplied data except where
! the value is IDMI - in which case the original value is retained
!
! This is intended to be used by an IO server task to cope with
! situations where either the client or the server has data for a given
! field - e.g. data length, record length and location in the file, which
! may vary depending upon e.g. who did the compression.
!----------------------------------------------------------------------
SUBROUTINE SetLookups(UNIT,DATA)

IMPLICIT NONE
INTEGER, INTENT(IN)          :: UNIT
INTEGER, POINTER             :: DATA(:)
INTEGER, POINTER             :: ipplook_in(:,:)
INTEGER, POINTER             :: ipplook   (:,:)
INTEGER                      :: i1,i2,ii,icode
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETLOOKUPS'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:SETLOOKUPS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
NULLIFY(ipplook_in,ipplook)
! Fetch the internal lookup
ipplook=>attachLookups(UNIT)

! Check that the data input matches the dimention of the one we have
IF ( SIZE(DATA) /= SIZE(ipplook) ) THEN
!$OMP CRITICAL(internal_write)
  WRITE(cmessage,'(A,I12,A,I12)')'supplied data size(',SIZE(DATA),         &
      ') mismatch with internal size (',SIZE(ipplook),')'
!$OMP END CRITICAL(internal_write)
  icode=99

  CALL mf_ereport(SubroutineName, Icode, Cmessage,UNIT)
END IF

!Allocate a 2D buffer to expand the 1D input buffer
ALLOCATE( ipplook_in (SIZE(ipplook,1),SIZE(ipplook,2) ) )

!Expand to 2D to make f90 happy
ii=1
DO i2=1,SIZE(ipplook,2)
  DO i1=1,SIZE(ipplook,1)
    ipplook_in(i1,i2)=DATA(ii)
    ii=ii+1
  END DO
END DO

!Merge the records
DO i2=1,SIZE(ipplook,2)
  DO i1=1,SIZE(ipplook,1)
    IF ( ipplook(i1,i2) == MF_Data_Missing ) THEN
      ipplook(i1,i2)=ipplook_in(i1,i2)
    END IF
  END DO
END DO
DEALLOCATE( ipplook_in )
NULLIFY(ipplook_in,ipplook)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE SetLookups

!----------------------------------------------------------------------
! Subroutine: StoreLookups
! Write the lookup table from memory into the file at
! the correct address, having checked it is open, and a
! pp file
!----------------------------------------------------------------------
SUBROUTINE StoreLookups(UNIT)

USE io, ONLY:                                                             &
    setpos,                                                                &
    buffout,                                                               &
    is_unit_open
USE ios, ONLY:                                                            &
    IOS_Write_Lookup
USE stwork_aux
IMPLICIT NONE
INTEGER,INTENT(IN)           :: UNIT
INTEGER                      :: lookupAddress
INTEGER                      :: icode
INTEGER                      :: len_io
REAL                         :: IOSTAT
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STORELOOKUPS'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:STORELOOKUPS'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
CALL checkUnit(um_file % UNIT,SubroutineName)

IF (.NOT. L_IO_Server) THEN
  ! Only do the flushing if we actually are storing a lookup table.
  IF (isActivePPFile(um_file % UNIT)) THEN
    lookupAddress = um_file % pp_meta % lookup_address
    IF (isUsingAsyncStash()) THEN

      ! Make sure all client side stash data (ie fields) are flushed
      ! to the IO server
      CALL flushPendingStashForUnit(um_file % UNIT)

      ! And then then send the lookup metadata after it
      CALL IOS_Write_Lookup                                       &
          (um_file % UNIT, lookupAddress, um_file % pp_meta % lookup_table)

    ELSE
      CALL setpos(um_file % UNIT,lookupAddress,icode)
      CALL buffout(um_file % UNIT, um_file % pp_meta % lookup_table,       &
          SIZE (um_file % pp_meta % lookup_table ), len_io, IOSTAT)
    END IF
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE StoreLookups


!----------------------------------------------------------------------
! Subroutine: storeAllLookups
! cycle over all units that potentially have cached lookups
! and flush to disk
!----------------------------------------------------------------------
SUBROUTINE storeAllLookups()

IMPLICIT NONE
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STOREALLLOOKUPS'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:STOREALLLOOKUPS'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
NULLIFY(um_file)
um_file => init_file_loop(handler="portio")
DO WHILE (ASSOCIATED(um_file))
  IF (um_file % pp_meta % managed ) CALL storeLookups(um_file % UNIT)
  um_file => um_file % next
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE storeAllLookups


!----------------------------------------------------------------------
! Subroutine: SynchroniseAll
! cycle over all units that potentially have cached lookups
! and flush to disk
!----------------------------------------------------------------------
SUBROUTINE SynchroniseAll()

USE io, ONLY:                                                             &
    ioDiskSynchronise,                                                     &
    is_unit_open
IMPLICIT NONE
REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SYNCHRONISEALL'
CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
    'MODEL_FILE:SYNCHRONISEALL'
TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
NULLIFY(um_file)
um_file => init_file_loop(handler="portio")
DO WHILE (ASSOCIATED(um_file))
  IF (is_unit_open(um_file % UNIT)) THEN
    IF (ASSOCIATED(um_file % pp_meta % lookup_table)) THEN
      CALL ioDiskSynchronise(um_file % UNIT)
    END IF
  END IF
  um_file => um_file % next
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE SynchroniseAll

!----------------------------------------------------------------------
! Function: isActivePPFile
! convenience function return whether a unit has cached lookups
! associated with it
!----------------------------------------------------------------------
LOGICAL FUNCTION isActivePPFile(UNIT)

IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
TYPE(um_file_type), POINTER :: um_file

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")
isActivePPFile=.FALSE.
IF (ASSOCIATED(um_file % pp_meta % lookup_table))                          &
    isActivePPFile=.TRUE.

END FUNCTION isActivePPFile

SUBROUTINE model_file_open (                                                 &
    UNIT,env,env_len,read_write,name_in_environ,error,                       &
    ioLocality,ioPolicy,allowRemap,fileType)
USE IOS_Common, ONLY: IOS_local_ro_files

IMPLICIT NONE
INTEGER, INTENT(IN)               :: UNIT
INTEGER, INTENT(IN), OPTIONAL     :: env_len
INTEGER, INTENT(IN), OPTIONAL     :: read_write
INTEGER, INTENT(IN), OPTIONAL     :: name_in_environ
INTEGER, OPTIONAL                 :: ioPolicy
INTEGER, OPTIONAL                 :: ioLocality
LOGICAL, OPTIONAL                 :: allowRemap
INTEGER, OPTIONAL                 :: error
INTEGER, OPTIONAL                 :: fileType
INTEGER                           :: errCode
INTEGER                           :: lFileType
CHARACTER(LEN=*)                  :: env
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MODEL_FILE_OPEN'
CHARACTER (LEN=*), PARAMETER                                               &
    :: SubroutineName=                                                     &
    'MODEL_FILE:MODEL_FILE_OPEN'
CHARACTER(LEN=fileNameLength) :: filename
TYPE(um_file_type), POINTER   :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errCode=99

IF (PRESENT(fileType)) THEN
  lFileType=fileType
ELSE
  lFileType=ioFileTypeUM
END IF

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")

IF (is_unit_open(um_file % UNIT)) THEN
  CALL mf_ereport(SubroutineName, errcode ,                                &
                  'Unit already open', um_file % UNIT)
END IF

IF (um_file % pp_meta % managed ) THEN
  CALL mf_ereport(SubroutineName, errCode ,                                &
      'Unit already managed by model_file', um_file % UNIT)
END IF

CALL file_open (                                                           &
    um_file % UNIT,env,env_len,read_write,name_in_environ,error,           &
    ioLocality,ioPolicy,allowRemap,lFileType)

IF ((io_external_control) .AND. (.NOT. (isReadOnly(um_file % UNIT))) ) THEN
  filename=TRIM(returnFilename(um_file % UNIT))
  CALL completion_open(um_file % UNIT, filename)
END IF

um_file % pp_meta % managed = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE model_file_open


!----------------------------------------------------------------------
! Subroutine model_file_close
!----------------------------------------------------------------------
SUBROUTINE model_file_close(UNIT,env,env_len,name_in_environ,delete,error)

IMPLICIT NONE
INTEGER, INTENT(IN)              :: UNIT
CHARACTER(LEN=*), INTENT(IN)     :: env
INTEGER, INTENT(IN), OPTIONAL    :: env_len  ! see io.F90 for details
INTEGER, INTENT(IN), OPTIONAL    :: delete
INTEGER, INTENT(IN), OPTIONAL    :: name_in_environ
INTEGER, INTENT(OUT), OPTIONAL   :: error
INTEGER                          :: errCode
REAL(KIND=jprb)                  :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MODEL_FILE_CLOSE'
CHARACTER (LEN=*), PARAMETER     :: SubroutineName=                        &
    'MODEL_FILE:MODEL_FILE_CLOSE'

LOGICAL :: test_complete_close
CHARACTER(LEN=fileNameLength) :: filename
TYPE(um_file_type), POINTER   :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(um_file)
um_file => get_file_by_unit(UNIT, handler="portio")

IF (.NOT. um_file % pp_meta % managed ) THEN
  errCode=99
  CALL mf_ereport(SubroutineName, errCode ,                                &
      'Unit not managed by model_file', um_file % UNIT)
END IF

CALL storeLookups(um_file % UNIT)

um_file % pp_meta % managed = .FALSE.

IF ((io_external_control) .AND.  &
    (.NOT. (isReadOnly(um_file % UNIT))) &
   ) THEN
  test_complete_close=.TRUE.
  filename=TRIM(returnfilename(um_file % UNIT))
ELSE
  test_complete_close=.FALSE.
END IF

CALL file_close(um_file % UNIT,env,env_len,name_in_environ,delete,error)

IF (test_complete_close) THEN
  CALL completion_close(um_file % UNIT,filename)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE model_file_close

!----------------------------------------------------------------------
! SUBROUTINE completion_open
! Suite Control : Remove any postfixed control file when file is opened
!----------------------------------------------------------------------

SUBROUTINE completion_open(UNIT,filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
CHARACTER(LEN=fileNameLength), INTENT(IN) :: filename
CHARACTER(LEN=fileNameLength+5) :: postfixed_name

REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPLETION_OPEN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

postfixed_name=TRIM(filename)//'.done'
CALL file_delete(TRIM(postfixed_name),silent=.TRUE.)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE completion_open

!----------------------------------------------------------------------
! SUBROUTINE completion_close
! Suite Control : Generate a postfixed control file when a file is
! closed
!----------------------------------------------------------------------

SUBROUTINE completion_close(UNIT,filename)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
CHARACTER(LEN=fileNameLength), INTENT(IN) :: filename

CHARACTER(LEN=fileNameLength+5) :: postfixed_name

REAL(KIND=jprb)              :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPLETION_CLOSE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

postfixed_name=TRIM(filename)//'.done'

! File Op
CALL file_touch(TRIM(postfixed_name))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE completion_close

END MODULE model_file


