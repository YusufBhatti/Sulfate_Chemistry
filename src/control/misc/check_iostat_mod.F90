! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE check_iostat_mod
! Description:
! This module is the interface to a utility subroutine
! that manages status of success/failure of
! FORTRAN namelist READs within the UM.
! Failures are handled in a standardised way across the UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, ONLY: newline

IMPLICIT NONE

CONTAINS

!--------------------------------------------

SUBROUTINE check_iostat(errorstatus, info, msg)

! if IOSTAT is nonzero from a READ then abort UM run.

IMPLICIT NONE

INTEGER,       INTENT(INOUT) :: errorstatus
CHARACTER(LEN=*), INTENT(IN) :: info
CHARACTER(LEN=*), INTENT(IN) :: msg

CHARACTER(LEN=*), PARAMETER  :: RoutineName='CHECK_IOSTAT'
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!Commented to prevent premature MPI initialisation by DrHook
!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Report fatal error (ABS(Errorstatus)) for all non-zero errors.
IF (ErrorStatus /= 0) THEN
  cmessage =                                                          newline//&
  'Error reading ' // TRIM(info) //                                   newline//&
  'IoMsg: '//TRIM(msg)//                                              newline//&
  'Please check input list against code.'
  errorstatus=ABS(ErrorStatus)
  CALL ereport (RoutineName, errorstatus, CMessage)
END IF

!Commented to prevent premature MPI initialisation by DrHook
!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_iostat

!--------------------------------------------

END MODULE check_iostat_mod
