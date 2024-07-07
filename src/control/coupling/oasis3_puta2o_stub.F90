#if !defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_puta2o()

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Description:
! Dummy version of oasis3_puta2o - should never be called.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!=====================================================================

CHARACTER(LEN=errormessagelength)   :: Message
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS3_PUTA2O'
INTEGER             :: ErrorStat        ! Return code:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!   0 = Normal exit
! +ve = Fatal Error
! -ve = Warning

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
ErrorStat   = 1
Message     = 'OASIS3-MCT Routines unavailable - see output.'

WRITE(umMessage,'(A)') '**ERROR**: oasis3_puta2o called but is unavailable.'
CALL umPrint(umMessage,src='oasis3_puta2o_stub')
WRITE(umMessage,'(A)') 'Check MCT cpp key is set'
CALL umPrint(umMessage,src='oasis3_puta2o_stub')


CALL ereport(RoutineName, ErrorStat, Message)
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis3_puta2o
#endif
