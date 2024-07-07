#if !defined(LFRIC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Description.
!
! An output message handler for writing text from the application.
! At application startup, all tasks have an open unit 6 and will write here,
! this can be modified via loaded options and API calls.
!
! External/Public:
!
!   umPrint(line,level,rank)       : Print the line (with options)
!
!   umPrintSetTarget(var,unit)     : Change the destination file/unit for output
!
!   umPrintLoadOptions()           : Load the options from the namelist
!
!   umPrintFlush(error,pe)         : Call flush on the output stream
!
!   umPrintFinalise()              : Any tidying up we need to do
!
!   umPrintSetLevel(level)         : Set the verbositty level of output
!
! Internal:
!
!   umPrintOpenStream()            : Open the destination file
!
!   umPrintError(source,message)   : Print message to error channel and abort
!
!   umPrintWarn(source,message)    : Print warning to error channel
!
!   umPrintExceptionHandler()      : To be called in case of fatal error
!                                  : recovery of buffered output.


MODULE umPrintMgr

USE um_types
USE filenamelength_mod, ONLY:                                                &
    filenamelength
USE UM_ParCore, ONLY:                                                        &
    mype,                                                                    &
    nproc_max
USE umFlush_mod, ONLY:                                                       &
    umFlush
USE yomhook, ONLY:                                                           &
    lhook,                                                                   &
    dr_hook
USE parkind1, ONLY:                                                          &
    jprb,                                                                    &
    jpim
USE errormessagelength_mod, ONLY:                                            &
    errormessagelength

IMPLICIT NONE


PRIVATE mype          ! These cause problems due to replication in
PRIVATE nproc_max     ! endgame, if they are not private

INTERFACE str
MODULE PROCEDURE                                                             &
   str_i
END INTERFACE

! ========================
! PARAMETERS
! ========================
INTEGER, PARAMETER            :: PrStatus_Min    = 1  ! Minimum output
INTEGER, PARAMETER            :: PrStatus_Normal = 2  ! Short info output
INTEGER, PARAMETER            :: PrStatus_Oper   = 3  ! Full info output
INTEGER, PARAMETER            :: PrStatus_Diag   = 4  ! Diagnostic output

!Shorter names for future use to avoid code clutter
INTEGER, PARAMETER            :: PrMin=PrStatus_Min
INTEGER, PARAMETER            :: PrNorm=PrStatus_Normal
INTEGER, PARAMETER            :: PrOper=PrStatus_Oper
INTEGER, PARAMETER            :: PrDiag=PrStatus_Diag

! Standard fortran units for default out/err
INTEGER, PARAMETER            :: oStreamUnitSystem=6
INTEGER, PARAMETER            :: eStreamUnitSystem=0

! Named Options for who writes
INTEGER, PARAMETER            :: outputAll  = 1     ! all pes
INTEGER, PARAMETER            :: outputZero = 2     ! pe 0 only
INTEGER, PARAMETER            :: outputZeroIOS = 3  ! pe 0 and head of IO server

! parameterisation of storage/buffering ammounts
INTEGER, PARAMETER :: maxLineLen=1024

! Unit to read the options file with
INTEGER, PARAMETER            :: prnt_namelist_unit=3

! Declare newline character
CHARACTER(LEN=1)              :: newline

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1


! ========================
! A buffer that clients
! can use for list
! directed writes
! ========================
CHARACTER (LEN=maxLineLen)    :: umMessage
!$OMP THREADPRIVATE (umMessage)


! ========================
! Runtime variables and
! initial defaults.
! ========================
INTEGER                       :: PrintStatus     = PrStatus_Normal
INTEGER                       :: printCount      = 0
INTEGER                       :: splitCount      = 0
INTEGER                       :: oStreamUnit     = 6
INTEGER                       :: eStreamUnit     = 0
LOGICAL                       :: printerActive   =.FALSE.
LOGICAL                       :: oStreamOpen     =.TRUE.
CHARACTER(LEN=filenamelength) :: oStreamFile
CHARACTER(LEN=filenamelength) :: oStreamVar      = 'STDOUT_FILE'

! ========================
! Options and controling
! namelist.
! ========================
LOGICAL                       :: prnt_src_pref    = .FALSE.
LOGICAL                       :: prnt_split_lines = .FALSE.
LOGICAL                       :: prnt_force_flush = .FALSE.
INTEGER                       :: prnt_writers     = outputAll
INTEGER                       :: prnt_paper_width = 80

NAMELIST /prnt_control/                                                      &
    prnt_src_pref,                                                           &
    prnt_split_lines,                                                        &
    prnt_force_flush,                                                        &
    prnt_paper_width,                                                        &
    prnt_writers

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UMPRINTMGR'

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintError(source,message)   : Print message to error channel and abort
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintError(source,message)
  ! We can't use ereport due to circular dependencies, so this is
  ! a mini error message handler that does as little as we can get away with!
USE um_abort_mod, ONLY: um_abort

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: message
CHARACTER(LEN=*), INTENT(IN) :: source
CHARACTER(LEN=*), PARAMETER  :: prefix="FATAL ERROR in umPrintMgr: "

!$OMP CRITICAL (ERROR_UNIT)
WRITE(eStreamUnit,'(A)')        prefix // source //' : '// message
!$OMP END CRITICAL (ERROR_UNIT)

CALL um_abort(1)

END SUBROUTINE umPrintError

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintWarn(source,message)    : Print warning to error channel
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintWarn(source,message)
  ! We can't use ereport due to circular dependencies, so this is
  ! a mini error message handler that does as little as we can get away with!
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: message
CHARACTER(LEN=*), INTENT(IN) :: source
CHARACTER(LEN=*), PARAMETER  :: prefix="Warning in umPrintMgr: "

!$OMP CRITICAL (ERROR_UNIT)
WRITE(eStreamUnit,'(A)')        prefix // source //' : '// message
!$OMP END CRITICAL (ERROR_UNIT)

END SUBROUTINE umPrintWarn

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintLoadOptions()           : load the options from the namelist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintLoadOptions()
USE FilenameLength_mod, ONLY:                                             &
    filenamelength
USE um_parcore, ONLY: mype, nproc
IMPLICIT NONE

INTEGER                        :: length     ! length of returned string
INTEGER                        :: error
CHARACTER (LEN=FileNameLength) :: NamelistFile = ' '
CHARACTER (LEN=errormessagelength) :: message = ' '
CHARACTER (LEN=errormessagelength) :: iomessage = ' '
REAL(KIND=jprb)                :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UMPRINTLOADOPTIONS'
INTEGER                        :: icode


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set newline character
newline = NEW_LINE('a')

! Call get_environment_variable directly to avoid circular dependencies
! from ereport_mod.
CALL GET_ENVIRONMENT_VARIABLE('IOSCNTL',NameListFile,length,error)

IF (error == 0 .AND. length /= 0) THEN
  IF (mype == 0) OPEN(Prnt_Namelist_Unit,FILE=NameListFile,ACTION='READ',  &
      FORM='FORMATTED', STATUS='OLD',IOSTAT=error, IOMSG=iomessage)

! broadcast error to all processes
  CALL gc_ibcast(101, 1, 0, nproc, icode, error)

  IF (error==0) THEN
    ! Read the I/O buffering control namelist too
    CALL read_nml_prnt_control(Prnt_Namelist_Unit, iomessage, error)
    IF (error /= 0 .AND. mype == 0) THEN
      WRITE(message,'(A,I6,A)')                                  newline //&
        'Error reading prnt_control namelist in '                        //&
         TRIM(NameListFile) //                                   newline //&
        'Read code was: ', error,                                newline //&
        'IoMsg: '//TRIM(iomessage)
      CALL umPrintError('umPrintLoadOptions',TRIM(message))
    END IF
    IF (mype == 0) CLOSE(Prnt_Namelist_Unit)
  ELSE
    IF (mype == 0) THEN
      WRITE(message,'(A,A,A)')'Failed to open IO control file (',          &
        TRIM(NameListFile),')'
      CALL umPrintWarn('umPrintLoadOptions',TRIM(message))
      CALL umPrintWarn('umPrintLoadOptions',TRIM(iomessage))
    END IF
  END IF
ELSE
  IF (mype == 0) THEN
    CALL umPrintWarn('umPrintLoadOptions',                                 &
      'Failed to get filename for IO control file from environment')
  END IF
END IF

! setup list of pe's which will output information
SELECT CASE (prnt_writers)
CASE (outputAll)
  printerActive=.TRUE.
CASE (outputZero)
  IF (mype==0) printerActive=.TRUE.
CASE (outputZeroIOS)
  IF (mype==0) printerActive=.TRUE.
! printeractive=.TRUE. is set on the head IOS PEs after all IO servers 
! have been initialised in IOS_init, not here
CASE DEFAULT
  CALL umPrintError('umPrintLoadOptions','Bad option for prnt_writers')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintLoadOptions

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintSetTarget(var,unit) : change the destination file/unit for output
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintSetTarget(var,UNIT)
IMPLICIT NONE

CHARACTER(LEN=*),OPTIONAL      :: var
INTEGER, OPTIONAL              :: UNIT
INTEGER                        :: length     ! length of returned string
INTEGER                        :: err
INTEGER                        :: fwidth
REAL(KIND=jprb)                :: zhook_handle
CHARACTER(LEN=filenamelength)  :: stdout_basename  ! base of filename
CHARACTER(LEN=20)              :: pe_descriptor
CHARACTER(LEN=*), PARAMETER    :: RoutineName='UMPRINTSETTARGET'
CHARACTER (LEN=maxLineLen)     :: locMessage
CHARACTER (LEN=errormessagelength)             :: cmessage

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(var)) THEN
  oStreamVar = var
END IF

! Call get_environment_variable directly to avoid circular dependencies
! from ereport_mod.
CALL GET_ENVIRONMENT_VARIABLE(oStreamVar,stdout_basename,length,err)

IF (err  /=  0 .OR. length == 0) THEN
  ! Environment variable for stdout has not been set, so abort.
  WRITE(cmessage, '(3A)') 'Failed to get value of ',TRIM(oStreamVar),      &
                          ' from environment'
  CALL umPrintError(RoutineName, cmessage)
END IF

! Discover the maximum number of characters needed to hold a pe number.
WRITE(pe_descriptor,'(I20)')nproc_max-1
pe_descriptor=ADJUSTL(pe_descriptor)
fwidth=LEN_TRIM(pe_descriptor)

! Construct an appropriate format descriptor for the file name
WRITE(pe_descriptor,'(A,I1,A,I1,A)')'(A,I',fwidth,                         &
    '.',fwidth,')'

! Construct the filename
WRITE(oStreamFile,pe_descriptor) TRIM(stdout_basename),mype
CALL umPrint(RoutineName//': Output set to '//TRIM(oStreamFile),           &
    level=PrStatus_Oper,src='umPrintMgr')

! Close any existing stream
IF (oStreamOpen) THEN
  WRITE(locMessage,'(A,I4)') RoutineName //                                &
      ': Closing unit ',oStreamUnit
  CALL umPrint(locMessage,level=PrStatus_Oper,src='umPrintMgr')
  CLOSE(oStreamUnit)
  oStreamOpen=.FALSE.
END IF

! Change the unit if requested
IF (PRESENT(UNIT)) THEN
  oStreamUnit=UNIT
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintSetTarget

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintOpenStream()            : Open the destination file
!
!   Note that this routine must be safe to be called from umPrint in a 
!   critical region. This means we can't do a umPrint from within as it can
!   lead to deadlock!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintOpenStream()
IMPLICIT NONE

REAL(KIND=jprb)                :: zhook_handle
CHARACTER(LEN=*), PARAMETER    :: RoutineName='UMPRINTOPENSTREAM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

OPEN(oStreamUnit,FILE=oStreamFile)
oStreamOpen=.TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintOpenStream

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintFlush(error,pe)         : call flush on the output stream
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintFlush(error,pe)
IMPLICIT NONE

INTEGER, INTENT(IN), OPTIONAL           :: pe
INTEGER, INTENT(OUT), OPTIONAL          :: error
LOGICAL                                 :: willFlush
INTEGER                                 :: Lerror
REAL(KIND=jprb)                         :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UMPRINTFLUSH'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

willFlush = .TRUE.
Lerror    = 0

IF (PRESENT(pe)) THEN
  IF (pe /= mype) willFlush=.FALSE.
END IF

IF (oStreamOpen .AND. willFlush) THEN
!$OMP CRITICAL (OUT_UNIT)
  CALL umFlush(oStreamUnit,Lerror)
!$OMP END CRITICAL (OUT_UNIT)
END IF

IF (PRESENT(error)) error = Lerror

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintFlush

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrint(line,level,rank)       : print the line (with options)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RECURSIVE SUBROUTINE umPrint(line,level,pe,src,model,                          &
                             UsrPrefix,HangIndent,stdErrorToo)
!$  USE omp_lib
IMPLICIT NONE

CHARACTER(LEN=*)            :: line
INTEGER, OPTIONAL           :: level
INTEGER, OPTIONAL           :: pe
CHARACTER(LEN=*), OPTIONAL  :: src
CHARACTER(LEN=*), OPTIONAL  :: model
CHARACTER(LEN=*), OPTIONAL  :: UsrPrefix       ! User prefix for each line
CHARACTER(LEN=*), OPTIONAL  :: HangIndent      ! Hanging indent for new lines
LOGICAL, OPTIONAL           :: stdErrorToo
INTEGER                     :: i
INTEGER                     :: BreakIndex
INTEGER                     :: thread
LOGICAL                     :: willPrint
LOGICAL                     :: first_print
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UMPRINT'
CHARACTER(LEN=maxLineLen)   :: finalBuffer
CHARACTER(LEN=maxLineLen)   :: LineToPrint


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Decide whether to print or not
willPrint=printerActive

IF (PRESENT(level)) THEN
  IF (printStatus < level) willPrint=.FALSE.
END IF

IF (PRESENT(pe)) THEN
  IF (pe /= mype) willPrint=.FALSE.
END IF

IF (.NOT. PRESENT(UsrPrefix)) THEN
  UsrPrefix = ''
END IF

IF (.NOT. PRESENT(HangIndent)) THEN
  HangIndent = ''
END IF

IF (willPrint) THEN
  ! Only want a single thread to do an open.
!$OMP CRITICAL (OUT_UNIT)
  IF (.NOT. oStreamOpen) THEN
    CALL umPrintOpenStream()
  END IF
!$OMP END CRITICAL (OUT_UNIT)

  ! Set BreakIndex to 1 so the while loop
  ! is performed at least once
  BreakIndex = 1
  first_print = .TRUE.
  DO WHILE (BreakIndex > 0)

    ! Check for linebreaks "newline" characters in the string to print
    !
    BreakIndex  = INDEX(line, newline)
    LineToPrint = ''

    IF (BreakIndex > 0) THEN
      LineToPrint(1:BreakIndex-1) = line(1:BreakIndex-1)
      line(1:BreakIndex) = ''
      line = ADJUSTL(line)
    ELSE
      LineToPrint = TRIM(line)
    END IF

    IF (first_print) THEN
      LineToPrint = UsrPrefix // ADJUSTL(LineToPrint)
    ELSE
      LineToPrint = UsrPrefix // HangIndent // ADJUSTL(LineToPrint)
    END IF

    IF ( PRESENT(src) .AND. prnt_src_pref ) THEN
      ! Prefix the output with src and recall dropping the argument
      finalBuffer=TRIM(ADJUSTL(src)) // ':' // TRIM(LineToPrint)
      CALL umPrint(TRIM(finalBuffer),                                      &
        model=model,                                                       &
        stdErrorToo=stdErrorToo)

    ELSE IF ( PRESENT(model) .AND. prnt_src_pref ) THEN
      ! Prefix the output with model and recall dropping the argument
      ! src has already been stripped above so no need to pass through.
      finalBuffer=TRIM(ADJUSTL(model)) // ':' // TRIM(LineToPrint)
      CALL umPrint(TRIM(finalBuffer),                                      &
        stdErrorToo=stdErrorToo)

    ELSE ! Nothing more to bolt to the string...

      IF (LEN_TRIM(LineToPrint) > prnt_paper_width) splitCount=splitCount+1

      IF (LEN_TRIM(LineToPrint) > prnt_paper_width .AND.                   &
         prnt_split_lines ) THEN

        CALL umPrint("NOTE: Long line has been split",                     &
          src='umPrintMgr',level=PrDiag)
        DO i=0,LEN_TRIM(LineToPrint)/prnt_paper_width
          CALL umPrint(LineToPrint(i*prnt_paper_width+1:                   &
            MIN( (i+1)*prnt_paper_width , LEN_TRIM(LineToPrint) )),        &
            stdErrorToo=stdErrorToo)
        END DO

      ELSE

!$OMP CRITICAL (INCR_COUNT)
        printCount=PrintCount+1
!$OMP END CRITICAL (INCR_COUNT)

        ! If we are in a parallel region, add the thread ID.

!$          IF ( omp_in_parallel() ) THEN
!$            thread=omp_get_thread_num()
!$            finalBuffer='['//TRIM(str_i(thread))//'] '//TRIM(LineToPrint)
!$          ELSE
        finalBuffer=TRIM(LineToPrint)
!$          END IF

        IF (PRESENT(stdErrorToo)) THEN
          IF (stdErrorToo) THEN
!$OMP CRITICAL (ERROR_UNIT)
            WRITE(eStreamUnit,'(A)')TRIM(finalBuffer)
!$OMP END CRITICAL (ERROR_UNIT)
          END IF
        END IF

!$OMP CRITICAL (OUT_UNIT)
        WRITE(oStreamUnit,'(A)')TRIM(finalBuffer)
!$OMP END CRITICAL (OUT_UNIT)

        ! flush the stream if options suggest to do so
        IF (prnt_force_flush .OR. PrintStatus >= PrDiag) CALL umPrintFlush()

      END IF
    END IF

    first_print = .FALSE.

  END DO ! While BreakIndex > 0
END IF ! willPrint

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrint

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintSetLevel(level)         : Set the verbositty level of output
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintSetLevel(level)
USE application_description_runtypes, ONLY: exe_type, exe_RCF
USE fort2c_portio_interfaces, ONLY: set_printstatus
IMPLICIT NONE

INTEGER, OPTIONAL             :: level
INTEGER, PARAMETER            :: maxPrStLen = 32
INTEGER                       :: length            ! length of returned string
INTEGER                       :: ErrorStatus       ! Error return code
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=maxPrStLen)     :: print_status    ! Content of env variable
CHARACTER(LEN=*), PARAMETER   :: RoutineName='UMPRINTSETLEVEL'
CHARACTER(LEN=maxPrStLen)     :: rcf_printstatus ! Content of env variable

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (PRESENT(level)) THEN
  PrintStatus = level
ELSE

  IF ( exe_type == exe_RCF ) THEN
    ! Call get_environment_variable directly to avoid circular dependencies
    ! from ereport_mod.
    CALL GET_ENVIRONMENT_VARIABLE( 'RCF_PRINTSTATUS',rcf_printstatus,length,  &
                                   ErrorStatus)
    IF (ErrorStatus /= 0 .OR. length == 0) THEN
      WRITE(umMessage,'(A,A,I4)')'Problem reading ',                      &
          '$RCF_PRINTSTATUS, output level will stay as ',printstatus
      CALL umPrintWarn(RoutineName,umMessage)
    ELSE   !  Env Var RCF_PRINTSTATUS has been set
      IF (rcf_printstatus(1:13) == 'PrStatus_Diag') THEN
        PrintStatus = PrStatus_Diag
      ELSE IF (rcf_printstatus(1:15) == 'PrStatus_Normal') THEN
        PrintStatus = PrStatus_Normal
      ELSE IF (rcf_printstatus(1:13) == 'PrStatus_Oper') THEN
        PrintStatus = PrStatus_Oper
      ELSE
        PrintStatus = PrStatus_Min
      END IF
    END IF

  ELSE  ! all applications except the RCF which has its own rcf_printstatus
    ! Call get_environment_variable directly to avoid circular dependencies
    ! from ereport_mod.
    CALL GET_ENVIRONMENT_VARIABLE('PRINT_STATUS',print_status,length,     &
                                  ErrorStatus)
    IF (ErrorStatus /= 0 .OR. length == 0) THEN
      WRITE(umMessage,'(A,A,I4)')'Problem reading ',                      &
          '$PRINT_STATUS, output level will stay as ',PrintStatus
      CALL umPrintWarn(RoutineName,umMessage)
    ELSE
      IF (print_status(1:13) == 'PrStatus_Diag') THEN
        PrintStatus = PrStatus_Diag
      ELSE IF (print_status(1:15) == 'PrStatus_Normal') THEN
        PrintStatus = PrStatus_Normal
      ELSE IF (print_status(1:13) == 'PrStatus_Oper') THEN
        PrintStatus = PrStatus_Oper
      ELSE
        PrintStatus = PrStatus_Min
      END IF
    END IF
  END IF
END IF
! Set PrintStatus in C code for control of I/O messages
CALL Set_Printstatus(PrintStatus)

WRITE(umMessage,'(A,A,I4)') RoutineName,                                  &
    ': PrintStatus initialised=',PrintStatus
CALL umPrint(umMessage,level=PrStatus_Oper,src='umPrintMgr')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintSetLevel

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   umPrintFinalise()              : Any tidying up we need to do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE umPrintFinalise()

IMPLICIT NONE

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UMPRINTFINALISE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,I10,A,I3,A)')'Note: ',splitCount,                      &
    ' output lines were too long for the paper width (',                   &
    prnt_paper_width,' columns)'
CALL umPrint(umMessage,src='umPrintMgr',level=PrOper)
WRITE(umMessage,'(A,I10,A)')'Congratulations: you printed things ',        &
    printCount,' times in this job. Goodbye :-)'
CALL umPrint(umMessage,src='umPrintMgr',level=PrDiag)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE umPrintFinalise

! Stringify an integer
CHARACTER(LEN=24) FUNCTION str_i(k)
INTEGER, INTENT(IN) :: k
!$OMP CRITICAL(internal_write)
WRITE (str_i, '(I20)') k
!$OMP END CRITICAL(internal_write)
str_i = ADJUSTL(str_i)
END FUNCTION str_i


SUBROUTINE umPrintExceptionHandler()
IMPLICIT NONE
!umPrintExceptionMode=.TRUE.
CALL umPrintWarn('umPrintExceptionHandler','Handler Invoked')
!IF (umPrintBuffering) THEN
!  CALL umPrintWarn('umPrintExceptionHandler','Dumping pending output')
!  CALL umPrintDumpBuffer()
!END IF
END SUBROUTINE umPrintExceptionHandler

SUBROUTINE read_nml_prnt_control(unit_in, iomessage, ErrorStatus)
! This routine uses a different approach than most to distribute namelists
! to ensure no MPL dependency. By using solely GCOM, the routines can be 
! used in serial programs. Rather than the more general case, where a type
! is defined, we instead just index an array. As only integer and logical 
! values are used this is OK.

USE um_parcore, ONLY: mype, nproc

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER, INTENT(OUT) :: ErrorStatus
CHARACTER(LEN=errormessagelength), INTENT(INOUT) :: iomessage
INTEGER :: icode

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_PRNT_CONTROL'

! set parameters for namelist reads
INTEGER, PARAMETER :: my_true            = 1
INTEGER, PARAMETER :: my_false           = 0
INTEGER, PARAMETER :: no_of_values       = 5
INTEGER, PARAMETER :: p_prnt_paper_width = 1
INTEGER, PARAMETER :: p_prnt_writers     = 2
INTEGER, PARAMETER :: p_prnt_src_pref    = 3
INTEGER, PARAMETER :: p_prnt_split_lines = 4
INTEGER, PARAMETER :: p_prnt_force_flush = 5

INTEGER :: my_nml(no_of_values)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=prnt_control, IOSTAT=ErrorStatus,          &
       IOMSG=iomessage)

  my_nml(p_prnt_paper_width) = prnt_paper_width
  my_nml(p_prnt_writers)     = prnt_writers
  IF (prnt_src_pref) THEN
    my_nml(p_prnt_src_pref) = my_true
  ELSE
    my_nml(p_prnt_src_pref) = my_false
  END IF

  IF (prnt_split_lines) THEN
    my_nml(p_prnt_split_lines) = my_true
  ELSE
    my_nml(p_prnt_split_lines) = my_false
  END IF

  IF (prnt_force_flush) THEN
    my_nml(p_prnt_force_flush) = my_true
  ELSE
    my_nml(p_prnt_force_flush) = my_false
  END IF

END IF

CALL gc_ibcast(102, no_of_values, 0, nproc, icode, my_nml)

IF (mype /= 0) THEN

  prnt_paper_width = my_nml(p_prnt_paper_width)
  prnt_writers     = my_nml(p_prnt_writers)

  IF (my_nml(p_prnt_src_pref) == my_true) THEN
    prnt_src_pref = .TRUE.
  ELSE
    prnt_src_pref = .FALSE.
  END IF

  IF (my_nml(p_prnt_split_lines) == my_true) THEN
    prnt_split_lines = .TRUE.
  ELSE
    prnt_split_lines = .FALSE.
  END IF

  IF (my_nml(p_prnt_force_flush) == my_true) THEN
    prnt_force_flush = .TRUE.
  ELSE
    prnt_force_flush = .FALSE.
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_prnt_control

END MODULE umPrintMgr
#endif
