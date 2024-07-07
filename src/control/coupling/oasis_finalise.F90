#if defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: oasis_finalise
!
!  Purpose: Finalise MPI communications for OASIS3-MCT
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Coupling

SUBROUTINE oasis_finalise

USE mod_prism, ONLY: prism_terminate_proto

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE um_types

IMPLICIT NONE
!
!  Local variables
!
INTEGER(KIND=integer32) :: icode_OASIS ! 32-bit OASIS return code

! No calls to Dr Hook:
! This routine is called outside the top-level Dr Hook calipers in um_shell.

!----------------------------------------------------------------------
!
! Call prism routine to close OASIS3-MCT instead of GC_EXIT
!
!----------------------------------------------------------------------

CALL prism_terminate_proto(icode_OASIS)

WRITE(umMessage,'(A,I6)')'oasis_finalise: Called prism_terminate_proto ',  &
                          icode_OASIS
CALL umPrint(umMessage,src='oasis_finalise')

RETURN

END SUBROUTINE oasis_finalise
#endif
