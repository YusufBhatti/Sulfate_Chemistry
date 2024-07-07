MODULE OASIS3_split_comm_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OASIS3_SPLIT_COMM_MOD'

CONTAINS

SUBROUTINE OASIS3_split_comm(colour)
!
! Description: This routine defines the communicator used by
!              OASIS3-MCT when excluding certain processes from 
!              coupling (typically IO server procs.)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!==================================================================

USE mod_prism, ONLY: PRISM_create_couplcomm

USE mpl, ONLY: mpl_undefined
USE oasis_atm_data_mod, ONLY: comm_in, icpl
USE um_types
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

LOGICAL, INTENT(IN) :: colour
                         ! The "colour" of this process defined by IOS
                         ! (or potentailly any other sub-division of
                         ! processes states whether or not the task
                         ! participates in atmos processes

! Local variables...
INTEGER (KIND=integer32) :: ierror ! Return code for OASIS3-MCT calls

INTEGER (KIND=integer32) :: comm_in32 ! 32 bit version of comm_in
                                      ! - the communicator defined
                                      ! by prism/oasis_get_localcomm
INTEGER (KIND=integer32) :: couplcomm ! OASIS3-MCT communicator

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_SPLIT_COMM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! We need to create a special communicator for OASIS which excludes the
! IOS processes.
! Set default to "not involved in coupling". It's critical that this value
! matches what OASIS is expecting in terms of MPI_UNDEFINED.
icpl = mpl_undefined

! Set all atmos-only procs to be involved in coupling
IF ( colour ) icpl = 1

! We need a 32 bit version of the original global communicator
comm_in32 = comm_in

! Tell OASIS3-MCT to create its special communicator. In principle we could
! use an existing communicator and call PRISM_set_couplecomm
! instead, but nothing suitable is available in the UM
! at this point. The point being that we need a communicator
! in which some of the processors (the IOS ones) are undefined, but this
! and certain other oasis initialisation calls are collective.
CALL PRISM_create_couplcomm(icpl, comm_in32, couplcomm, ierror)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE OASIS3_split_comm

END MODULE OASIS3_split_comm_mod
