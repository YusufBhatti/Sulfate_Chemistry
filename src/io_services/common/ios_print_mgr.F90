! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

MODULE IOS_print_mgr

USE umPrintMgr, ONLY:                                                       &
    PrintStatus,                                                             &
    umPrintSetLevel,                                                         &
    umPrint,                                                                 &
    umPrintFlush,                                                            &
    PrMin,                                                                   &
    PrNorm,                                                                  &
    PrOper,                                                                  &
    PrDiag
USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!As well as literal meanings for print levels for
!IOS_Verbosity these are also codes for the intent
!argument of a 'config' action
!e.g. 1-5 sets the verbosity to that level,
!       6 switches on the c timing layer

  ! Minimum output, only essential messages
INTEGER,PARAMETER :: IOS_PrStatus_Min    = PrMin
! Normal informative messages + warnings
INTEGER,PARAMETER :: IOS_PrStatus_Normal = PrNorm
! Operational status, all informative messages
INTEGER,PARAMETER :: IOS_PrStatus_Oper   = PrOper
! All informative + extra diagnostic messages
INTEGER,PARAMETER :: IOS_PrStatus_Diag   = PrDiag
! All informative + diagnostic + debugging messages
INTEGER,PARAMETER :: IOS_PrStatus_Debug  = 5

INTEGER           :: IOS_Verbosity            = imdi

LOGICAL           :: IOS_print_start_time     = .FALSE.

INTEGER, PARAMETER, PRIVATE   :: max_line_len=1024

CHARACTER (LEN=max_line_len)  :: IOS_message
!$OMP THREADPRIVATE (IOS_message)

CONTAINS

  ! Set the printing level
SUBROUTINE IOS_PrintSetLevel(level,setParent,setFromParent)
IMPLICIT NONE
INTEGER, OPTIONAL          :: level
LOGICAL, OPTIONAL          :: setFromParent
LOGICAL, OPTIONAL          :: setParent

INTEGER                    :: Llevel
LOGICAL                    :: LsetFromParent
LOGICAL                    :: LsetParent

LLevel= imdi
IF (PRESENT(level)) LLevel=level
LsetFromParent=.FALSE.
IF (PRESENT(setFromParent)) LsetFromParent=setFromParent
LsetParent=.FALSE.
IF (PRESENT(setParent)) LsetParent=setParent

IF (LsetFromParent) THEN
  IOS_Verbosity=PrintStatus
ELSE IF (Llevel /= imdi) THEN
  IOS_Verbosity=Llevel
  IF (LsetParent) CALL umPrintSetLevel(Llevel)
ELSE
  CALL IOS_Print('IOS_PrintSetLevel: No level provided',src='ios_print_mgr')
END IF

END SUBROUTINE IOS_PrintSetLevel

SUBROUTINE IOS_print_flush(pe)
IMPLICIT NONE
INTEGER, INTENT(IN), OPTIONAL :: pe

CALL umPrintFlush(pe=pe)

END SUBROUTINE IOS_print_flush

SUBROUTINE IOS_print(line,src,level)
  ! Either send output to the umPrintMgr subsystem for coupled jobs
  ! or send it to stdout otherwise.

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: line
CHARACTER(LEN=*), OPTIONAL   :: src
INTEGER, OPTIONAL            :: level

CALL umPrint(TRIM(line),src=src,model='ios',level=level)

END SUBROUTINE IOS_print

END MODULE IOS_print_mgr
