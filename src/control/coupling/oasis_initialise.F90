#if defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: oasis_initialise
!
!  Purpose: Initialise MPI communications for OASIS3-MCT
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Coupling
!
!=====================================================================

SUBROUTINE oasis_initialise (mype,nproc_max,comm_in,thread_level_set)

USE oasis3_atmos_init_mod, ONLY: oasis3_atmos_init

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
!  Subroutine arguments
!
INTEGER :: mype
INTEGER :: nproc_max
INTEGER :: comm_in
INTEGER :: thread_level_set
!
!  Local variables
!
INTEGER :: icode       ! Work - Internal return code

CHARACTER(LEN=errormessagelength) :: cmessage ! Work - Internal error message


INTEGER :: dummy_comm     ! Dummy communicator for OASIS


! No calls to Dr Hook:
! This routine is called outside the top-level Dr Hook calipers in um_shell.

!----------------------------------------------------------------------
!
!  Call routine to initialise OASIS
!
!----------------------------------------------------------------------

      ! The key thing here is to get hold of the
      ! communicator defined for us by PRISM and then
      ! use that in GCOM rather than letting GCOM define
      ! its own MPI_COMM_WORLD.

comm_in=-999
CALL oasis3_atmos_init(comm_in,icode,cmessage)

    ! Check that MPI (or other) communication method
    ! is initialised. Discard the communicator returned
    ! from this call since we'll use an OASIS defined one.

CALL gc_init_intro_thread(dummy_comm, thread_level_set)

    ! Do all the initialisation with the correct communicator

CALL gc_init_final(mype,nproc_max,comm_in)

WRITE(umMessage,'(A,I7,I7,I8)') "oasis_initialise: GCOM for OASIS", &
                        mype, nproc_max, comm_in
CALL umPrint(umMessage,src='oasis_initialise')

    ! We need to find out which coupling fields to define, which of them
    ! are input and which output, which grids they're on, etc.
    ! We do this as soon as we possibly can because this infomation is
    ! potentially needed by IOS processes for synchronisation purposes and
    ! because if we do it here, we have the opportunity to employ the main
    ! NAMELIST file (since it's a namelist read) before the main UM code 
    ! needs to access it for the main UM control variables.

! DEPENDS ON: OASIS_read_translist
CALL OASIS_read_translist()

RETURN

END SUBROUTINE oasis_initialise
#endif
