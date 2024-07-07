#if !defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: oasis_finalise
!
!  Purpose: Finalise MPI communications for OASIS3-MCT (dummy)
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Coupling

SUBROUTINE oasis_finalise

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ereport_mod, ONLY: ereport
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
!
!  Local variables
!
INTEGER :: ErrorStat

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS_FINALISE'
CHARACTER(LEN=errormessagelength)   :: Message

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

ErrorStat   = 1
Message     = 'OASIS3-MCT Routines unavailable - see output.'

WRITE(umMessage,'(A)') '**ERROR**: oasis_finalise called but is unavailable.'
CALL umPrint(umMessage,src='oasis_finalise_stub')
WRITE(umMessage,'(A)') 'Check MCT cpp key is set'
CALL umPrint(umMessage,src='oasis_finalise_stub')

CALL ereport(RoutineName, ErrorStat, Message)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE oasis_finalise
#endif
