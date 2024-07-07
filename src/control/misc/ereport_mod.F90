#if !defined(LFRIC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.

MODULE Ereport_Mod
USE um_types
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParCore
USE umPrintMgr
USE um_abort_mod, ONLY: um_abort
IMPLICIT NONE


INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! String Lengths
INTEGER, PARAMETER :: replinelen = 80
INTEGER, PARAMETER :: markerlen  = 3
INTEGER, PARAMETER :: repstrlen  = 7
INTEGER, PARAMETER :: minspaces  = 2 

! Banner strings
CHARACTER (LEN=replinelen), PARAMETER  :: astline = REPEAT("?", replinelen)
CHARACTER (LEN=markerlen),  PARAMETER  :: qmark   = REPEAT('?', markerlen)
CHARACTER (LEN=markerlen),  PARAMETER  :: exmark  = REPEAT('!', markerlen)
CHARACTER (LEN=markerlen),  PARAMETER  :: UsrPrefix  = '?'                  //&
                                          REPEAT(' ',markerlen-1)

CHARACTER (LEN=*), PARAMETER, PRIVATE  :: ModuleName = 'EREPORT_MOD'

! Create an interface to allow 32 or 64 bit errorstatuses to be used.
INTERFACE ereport
MODULE PROCEDURE ereport64, ereport32
END INTERFACE ereport

! Counter for number of warnings
INTEGER :: numwarnings  = 0

! Only make the interface available publicly
PRIVATE ereport64
PRIVATE ereport32

CONTAINS


!  Subroutine Ereport - handles errors/warnings
!
! Description:
!   The standard reconfiguration error and warning reporting routine.
!   Prints the relevant error and exits if required.
!
! Method:
!   ErrorStatus > 0 is a hard error
!   ErrorStatus < 0 is a warning
!   ErrorStatus = 0 is hunky-dory

SUBROUTINE Ereport64 ( RoutineName, ErrorStatus, Cmessage )

USE UM_Types, ONLY: &
  integer64
IMPLICIT NONE

! Arguments
CHARACTER (LEN=*), INTENT(IN)   :: RoutineName
CHARACTER (LEN=*), INTENT(IN)   :: Cmessage

INTEGER(KIND=Integer64) , INTENT(INOUT)        :: ErrorStatus
! We can't intent this because the UM passes constants.

! Local scalars
INTEGER                       :: UNIT,i
REAL(KIND=jprb)               :: zhook_handle

! Ensure that the write buffer used for sending cmessage to umPrint
! function is of adequate size
CHARACTER (LEN=LEN(cmessage)+repstrlen+10) :: localWriteBuffer

CHARACTER (LEN=*), PARAMETER  :: msg ='Job aborted from ereport.'
CHARACTER (LEN=*), PARAMETER  :: ThisRoutineName = 'EREPORT64'

! Temporary parameters to avoid causing problems with triple ? in the
! AIX Fortran preprocessor
LOGICAL                       :: stderr = .FALSE.
INTEGER                       :: repType,strlen,iter
INTEGER                       :: paddingl,paddingr
CHARACTER (LEN=14)            :: repString
CHARACTER (LEN=repstrlen)     :: repStringLow
CHARACTER (LEN=markerlen*2)   :: repBanner

IF (lhook) CALL dr_hook(ModuleName//':'//ThisRoutineName,zhook_in,zhook_handle)

repType=0

IF (ErrorStatus > 0) THEN
  ! Error
  repType      = 1
  repString    = 'ERROR'
  repStringLow = 'Error'
  repBanner    = qmark//exmark
  ! If ereport is reporting an error make sure pe output file is active
  printerActive = .TRUE.
  ! Ensure fatal error messages from all PEs go to stderr:
  stderr       = .TRUE.   
ELSE IF (ErrorStatus < 0) THEN
  ! Warning
  repType      = 2
  repString    = 'WARNING'
  repStringLow = 'Warning'
  repBanner    = qmark//qmark
END IF

repStringLow = ADJUSTR(repStringLow)

! ALL ereport messages from PE0 are duplicated to stderr:
IF (mype == 0) stderr = .TRUE.

IF (repType > 0) THEN ! some kind of report

  ! Process localWriteBuffer before any umPrint calls, in case umMessage was
  ! passed as the cmessage argument. Otherwise, if it was, we will over-write
  ! it later, thus losing the message.

  ! If changing the lengths of the strings prepended to cmessage please also
  ! update the LEN of localWriteBuffer (defined in the header of this routine)
  WRITE(localWriteBuffer,'(A)')  repStringLow // ' message: '// cmessage

  strlen = LEN(TRIM(repString))
  iter = (replinelen - strlen - (2*minspaces) ) / (4*markerlen) 
  paddingl = (replinelen - strlen - (iter*4*markerlen) ) / 2
  paddingr = (replinelen - strlen - (iter*4*markerlen) - paddingl)

  CALL umPrint(newline//astline, src='ereport_mod',stdErrorToo=stderr)
  CALL umPrint(                                                                &
    REPEAT(repBanner,iter)                                                   //&
    REPEAT(' ',paddingl) // TRIM(repString) // REPEAT(' ',paddingr)          //&
    REPEAT(repBanner,iter),                                                    &
    src='ereport_mod',stdErrorToo=stderr)

  WRITE(umMessage,'(A,I0)') repStringLow // ' code: ', ErrorStatus

  CALL umPrint(umMessage, src='ereport_mod', UsrPrefix=UsrPrefix,              &
               stdErrorToo=stderr)

  WRITE(umMessage,'(A)')  repStringLow // ' from routine: ' // RoutineName
  CALL umPrint(umMessage, src='ereport_mod',UsrPrefix=UsrPrefix,               &
               stdErrorToo=stderr)

  CALL umPrint(localWriteBuffer, src='ereport_mod',UsrPrefix=UsrPrefix,        &
       HangIndent=REPEAT(' ',strlen+1), stdErrorToo=stderr)

  WRITE(umMessage,'(A,I0)') repStringLow // ' from processor: ', mype
  CALL umPrint(umMessage, src='ereport_mod',UsrPrefix=UsrPrefix,               &
               stdErrorToo=stderr)

  WRITE(umMessage,'(A,I0)') repStringLow // ' number: ', numwarnings
  CALL umPrint(umMessage, src='ereport_mod',UsrPrefix=UsrPrefix,               &
               stdErrorToo=stderr)

  CALL umPrint(astline//newline,   src='ereport_mod',stdErrorToo=stderr)

  IF (repType==1) THEN
    ! Call registered callbacks at termination
    CALL um_abort(INT(ErrorStatus))
  END IF
  IF (repType==2) THEN
    ! Reset ErrorStatus, update warning count.
    numwarnings = numwarnings + 1
    ErrorStatus = 0
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//ThisRoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE Ereport64

! The 32 bit interface to the above routine
SUBROUTINE Ereport32 ( RoutineName, ErrorStatus32, Cmessage )

USE UM_Types, ONLY: &
  integer32,         &
  integer64
IMPLICIT NONE

! Arguments
CHARACTER (LEN=*), INTENT(IN)   :: RoutineName
CHARACTER (LEN=*), INTENT(IN)   :: Cmessage

INTEGER(KIND=Integer32)           :: ErrorStatus32

! Local variables
INTEGER(KIND=Integer64)           :: ErrorStatus64
CHARACTER (LEN=*), PARAMETER      :: ThisRoutineName = 'EREPORT32'

REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//ThisRoutineName,zhook_in,zhook_handle)

! Copy 32 bit errorstatus to 64 bit temporary
ErrorStatus64 = ErrorStatus32

! Call the ereport64 routine which does all the work.
CALL Ereport64(RoutineName, ErrorStatus64, Cmessage)

! Copy the 64 bit errorstatus (can be reset above) to the 32 bit version.
ErrorStatus32 = ErrorStatus64

IF (lhook) CALL dr_hook(ModuleName//':'//ThisRoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Ereport32


! Subroutine to display the number of warnings generated
SUBROUTINE Ereport_finalise( )

IMPLICIT NONE

INTEGER                       :: unset = -99
INTEGER                       :: UNIT,i
INTEGER                       :: out_units(2)
INTEGER                       :: strlen,iter
INTEGER                       :: paddingl,paddingr
CHARACTER (LEN=14)            :: repString
CHARACTER (LEN=8)             :: repStringLow
CHARACTER (LEN=6)             :: repBanner
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'EREPORT_FINALISE'
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (numwarnings > 0) THEN
  repString    = 'CAUTION'
  repStringLow = 'Caution'
  repBanner    = qmark//qmark

  repStringLow = ADJUSTR(repStringLow)

  strlen = LEN(TRIM(repString))
  iter = (replinelen - strlen - (2*minspaces) ) / (4*markerlen) 
  paddingl = (replinelen - strlen - (iter*4*markerlen) ) / 2
  paddingr = (replinelen - strlen - (iter*4*markerlen) - paddingl)

  CALL umPrint(newline//astline, src='ereport_mod')
  CALL umPrint(                                                             &
    REPEAT(repBanner,iter)                                                //&
    REPEAT(' ',paddingl) // TRIM(repString) // REPEAT(' ',paddingr)       //&
    REPEAT(repBanner,iter),                                                 &
    src='ereport_mod')

  WRITE(umMessage,'(A,I0,A)') repStringLow // ' ' //                        &
      'This run generated ',numwarnings,' warnings'
  CALL umPrint(umMessage,  UsrPrefix=UsrPrefix, src='ereport_mod')
  CALL umPrint(astline//newline, src='ereport_mod')
  CALL umPrintFlush()
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Ereport_finalise

! Note that code is PASS BY VALUE !
SUBROUTINE creport(routineName, code, message)                             &
    BIND(c,NAME="f_ereport")

USE, INTRINSIC :: ISO_C_BINDING
USE f_shum_string_conv_mod, ONLY: f_shum_c2f_string

IMPLICIT NONE

CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN)   :: RoutineName(*)
CHARACTER(KIND=C_CHAR,LEN=1), INTENT(IN)   :: message(*)
INTEGER(KIND=C_INT64_T), VALUE :: code
INTEGER(KIND=Integer64)  ::  icode

icode = code
CALL ereport64(f_shum_c2f_string(RoutineName), icode,                     &
               f_shum_c2f_string(message)                                )

END SUBROUTINE creport

END MODULE Ereport_Mod
#endif
