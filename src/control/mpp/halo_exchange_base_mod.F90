! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

!
! This module provides the base data for the methods used to updating
! boundaries (halos) on fields.

MODULE halo_exchange_base_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out = 1

! The MPL communicator used to exchange the halos.
INTEGER :: halos_comm


! The ID of the neighbour processes
INTEGER, PARAMETER :: all_neighbours = 8
INTEGER            :: he_neighbour(all_neighbours)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HALO_EXCHANGE_BASE_MOD'

CONTAINS

!==============================================================================
! Purpose:
! Initialise the data used for halo exchange.
! Must be called before any swap_bounds routines.
!==============================================================================
SUBROUTINE halo_exchange_base_init()

USE um_parvars, ONLY: full_neighbour
IMPLICIT NONE

INTEGER :: ierror
INTEGER :: global_comm

REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HALO_EXCHANGE_BASE_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

! Get the global communicator
CALL gc_get_communicator(global_comm,ierror)

!------------------------------------------------------------------
! Halos communicator

! Duplicate the global communicator so that the swap_bounds routines
! uses their own communicator.
CALL mpl_comm_dup(global_comm, halos_comm, ierror)


!------------------------------------------------------------------
! Neighbour processes

! The ID of the neighbour process is already available.
he_neighbour(1:8) = full_neighbour(1:8)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE halo_exchange_base_init


!==============================================================================
! Purpose:
! Finalise the data used for halo exchange. 
! Must be called after all swap_bounds routines have finished.
!==============================================================================
SUBROUTINE halo_exchange_base_finalise()

IMPLICIT NONE

REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HALO_EXCHANGE_BASE_FINALISE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE halo_exchange_base_finalise


END MODULE halo_exchange_base_mod
