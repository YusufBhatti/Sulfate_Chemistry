#if !defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_advance_date()

  ! Description: This subroutine is a stub routine for
  !              OASIS3_ADVANCE_DATE. Run time logic
  !              should prevent it from actually being called.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !=============================================================


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


CHARACTER(LEN=errormessagelength)   :: Message
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS3_ADVANCE_DATE'
INTEGER             :: ErrorStat        ! Return code:
!   0 = Normal exit
! +ve = Fatal Error
! -ve = Warning

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
ErrorStat   = 1
Message     = 'OASIS3-MCT Routines unavailable - see output.'

CALL umPrint( '**ERROR**: oasis3_advance_date unavailable.', &
    src='oasis3_advance_date_stub')
CALL umPrint( 'Check MCT cpp key is set',src='oasis3_advance_date_stub')

CALL ereport(RoutineName, ErrorStat, Message)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE oasis3_advance_date
#endif
