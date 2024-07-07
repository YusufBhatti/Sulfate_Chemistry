! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_closeoutput_mod

USE umPrintMgr,  ONLY: umPrint, umMessage, newline
USE ereport_mod, ONLY: ereport

IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_CLOSEOUTPUT_MOD'

CONTAINS


SUBROUTINE stph_closeoutput()

! Closes the RPSEED output stream if it is open

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stochastic_physics_run_mod,  ONLY: stphseed_unit

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


! Local Variables
INTEGER :: icode ! Error return code

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'STPH_CLOSEOUTPUT'
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CLOSE(UNIT=stphseed_unit, IOSTAT=icode, IOMSG=iomessage)
IF (icode > 0) THEN
  WRITE(umMessage,'(A)')                                                       &
    "RPSEED OUT : Failed to close output file:"//                     newline//&
    "IoMsg: "//TRIM(iomessage)
  CALL ereport(routinename, icode, umMessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stph_closeoutput
END MODULE stph_closeoutput_mod
