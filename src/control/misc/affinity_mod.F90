! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains functionality to report affinity information.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C Code
!

!-------------------------------------------------------------------------------
! Full version of the module for Linux machines. Non-Linux version follows.
!-------------------------------------------------------------------------------

#if defined(__linux__)

MODULE affinity_mod
USE fort2c_affinity_interfaces_mod
USE umPrintMgr, ONLY: umPrint, umMessage
USE parkind1,   ONLY: jpim, jprb
USE yomhook,    ONLY: lhook, dr_hook
USE mpl
!$ USE omp_lib
IMPLICIT NONE
PRIVATE 

!-------------------------------------------------------------------------------
! Public routines
!-------------------------------------------------------------------------------

PUBLIC :: affinity_write

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER :: blank_string = ''

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AFFINITY_MOD'

CONTAINS

!-------------------------------------------------------------------------------
! Top-level routines - platform dependent
!-------------------------------------------------------------------------------

SUBROUTINE affinity_write(communicator, writer_rank, long_form)
  !-----------------------------------------------------------------------------
  ! SYNOPSIS
  !   call affinity_write(communicator, writer_rank, long_form)
  !
  ! DESCRIPTION
  !   Calls the appropriate routine to output the long- or short- form of
  !   affinity information, depending on logical input.
  !
  ! ARGUMENTS
  !   * communicator - The MPI communicator whose members will participate.
  !   * writer_rank  - The MPI rank that is to write the information to its PE
  !                    file.
  !   * long_form    - If true, calls the long-form routine. Default is
  !                    short-form.
  ! NOTES
  !   Broadcasts value of long_form from the writer rank.  Thereby handles the
  !   case where the print status is not the same on all PEs in the passed
  !   communicator.  Otherwise an MPI hang is likely, since different PEs could
  !   find themselves in different affinity routines.
  !-----------------------------------------------------------------------------
  IMPLICIT NONE

  !Arguments
  INTEGER, INTENT(IN) :: communicator
  INTEGER, INTENT(IN) :: writer_rank
  LOGICAL, INTENT(IN) :: long_form

  !Local variables
  REAL(KIND=jprb) :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AFFINITY_WRITE'
  INTEGER         :: ierr
  LOGICAL         :: bcast_long_form

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

  bcast_long_form = long_form
  CALL mpl_bcast(bcast_long_form, 1, mpl_logical,    &
                 writer_rank, communicator, ierr)

  IF(bcast_long_form) THEN
    CALL affinity_write_long(communicator, writer_rank)
  ELSE
    CALL affinity_write_short(communicator, writer_rank)
  END IF

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  RETURN
END SUBROUTINE affinity_write

!-------------------------------------------------------------------------------
! Short and long-form routines.
!-------------------------------------------------------------------------------

SUBROUTINE affinity_write_short(communicator, writer_rank)
  !-----------------------------------------------------------------------------
  ! SYNOPSIS
  !   call affinity_write_short(communicator, writer_rank)
  !
  ! DESCRIPTION
  !   Writes a few lines in the PE file for MPI rank writer_rank.  Gives:
  !     * The number of cores on a node of the machine
  !     * The maximum, minimum and mean number of cores that threads are able to
  !       migrate across, reduced over all MPI tasks. 
  !
  ! ARGUMENTS
  !   * communicator - The MPI communicator to reduce over
  !   * writer_rank  - The MPI rank that is to write the information to its PE
  !                    file.
  !-----------------------------------------------------------------------------
  IMPLICIT NONE

  !Arguments
  INTEGER, INTENT(IN) :: communicator
  INTEGER, INTENT(IN) :: writer_rank

  !Internal variables
  LOGICAL          :: on_writer
  INTEGER          :: max_available_cpus
  INTEGER          :: num_available_cpus

  INTEGER          :: ierr

  INTEGER          :: comm_size
  INTEGER          :: my_rank

  INTEGER          :: num_cpus_max
  INTEGER          :: num_cpus_min
  INTEGER          :: num_cpus_sum
  REAL             :: num_cpus_mean

  REAL(KIND=jprb)  :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AFFINITY_WRITE_SHORT'

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  CALL mpl_comm_size(communicator, comm_size, ierr)
  CALL mpl_comm_rank(communicator, my_rank,   ierr)

  on_writer = (my_rank == writer_rank)

  max_available_cpus = fort2c_affinity_max_available_cpus()
  num_available_cpus = fort2c_affinity_num_available_cpus()

  !Reductions and averaging
  CALL mpl_reduce(num_available_cpus, num_cpus_max, 1, mpl_integer,   &
                  mpl_max, writer_rank, communicator, ierr)
  CALL mpl_reduce(num_available_cpus, num_cpus_min, 1, mpl_integer,   &
                  mpl_min, writer_rank, communicator, ierr)
  CALL mpl_reduce(num_available_cpus, num_cpus_sum, 1, mpl_integer,   &
                  mpl_sum, writer_rank, communicator, ierr)

  num_cpus_mean =  REAL(num_cpus_sum) / REAL(comm_size)

  !Writer has the reduction results. So write them.
  IF(on_writer) THEN
    CALL umPrint(blank_string, src='affinity')
    WRITE(umMessage,'(a)') '--> AFFINITY REPORT (short form) <--'
    CALL umPrint(umMessage, src='affinity')
    CALL umPrint(blank_string, src='affinity')

    WRITE(umMessage,'(a,1x,i0)') 'Cores in a node:', max_available_cpus
    CALL umPrint(umMessage, src='affinity')
    CALL umPrint(blank_string, src='affinity')

    WRITE(umMessage,'(a)') 'Cores available to each thread across PEs:'
    CALL umPrint(umMessage, src='affinity')

    WRITE(umMessage,'(a,f12.2)') 'Mean:', num_cpus_mean
    CALL umPrint(umMessage,src='affinity')

    WRITE(umMessage,'(a,1x,i0)') 'Minimum:', num_cpus_min
    CALL umPrint(umMessage,src='affinity')

    WRITE(umMessage,'(a,1x,i0)') 'Maximum:', num_cpus_max
    CALL umPrint(umMessage,src='affinity')

    CALL umPrint(blank_string, src='affinity')
  END IF

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  
  RETURN
END SUBROUTINE affinity_write_short

SUBROUTINE affinity_write_long(communicator, writer_rank)
  !-----------------------------------------------------------------------------
  ! SYNOPSIS
  !   call affinity_write_long(communicator, writer_rank)
  !
  ! DESCRIPTION
  !   Writes long-form affinity information. Grid output, showing where
  !   individual threads are running, for each MPI task.  Also gives the number
  !   of cores that each MPI task's threads are free to migrate across.  A key
  !   to the output format is provided with the output itself.
  !
  ! ARGUMENTS
  !   * communicator - The MPI communicator, whose members will provide data.
  !   * writer_rank  - The MPI rank that is to write the information to its PE
  !                    file.
  !-----------------------------------------------------------------------------
  IMPLICIT NONE

  !Arguments
  INTEGER, INTENT(IN) :: communicator
  INTEGER, INTENT(IN) :: writer_rank

  !Internal variables
  LOGICAL          :: on_writer
  INTEGER          :: max_available_cpus
  INTEGER          :: num_available_cpus
  INTEGER          :: core_id
  INTEGER          :: thread_id
  INTEGER          :: buffer_size
  CHARACTER(LEN=1) :: hex_char
  CHARACTER(LEN=3) :: containment_string
  INTEGER          :: recv_status(mpl_status_size)
  INTEGER          :: send_status(mpl_status_size)
  INTEGER          :: send_request

  INTEGER          :: isubstr
  INTEGER          :: irank
  INTEGER          :: ierr

  INTEGER          :: comm_size
  INTEGER          :: my_rank

  ! The Cray compiler (CCE/8.3.4) doesn't correctly handle allocatable character
  ! strings: generates an internal compiler error when certain debugging flags 
  ! are used.  Workaround needed on the Cray for now.  Can remove once the
  ! compiler bug is fixed.

#if ! defined (CRAY_FORTRAN)
  CHARACTER(LEN=:), ALLOCATABLE :: mask
  CHARACTER(LEN=:), ALLOCATABLE :: send_buffer
  CHARACTER(LEN=:), ALLOCATABLE :: recv_buffer
#else
  CHARACTER(LEN=200) :: mask          = ' '
  CHARACTER(LEN=200) :: send_buffer   = ' '
  CHARACTER(LEN=200) :: recv_buffer   = ' '
#endif

  REAL(KIND=jprb) :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AFFINITY_WRITE_LONG'

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  CALL mpl_comm_size(communicator, comm_size, ierr)
  CALL mpl_comm_rank(communicator, my_rank,   ierr)

  on_writer = (my_rank == writer_rank)
  
  max_available_cpus = fort2c_affinity_max_available_cpus()
  buffer_size = max_available_cpus+50

#if ! defined (CRAY_FORTRAN)
  ALLOCATE(CHARACTER(LEN=max_available_cpus) :: mask)
  ALLOCATE(CHARACTER(LEN=buffer_size)        :: send_buffer)
  ALLOCATE(CHARACTER(LEN=buffer_size)        :: recv_buffer)
#endif
  
  DO isubstr = 1, max_available_cpus
    mask(isubstr:isubstr) = '.'
  END DO
  
!$OMP  PARALLEL DEFAULT(NONE)                            &
!$OMP& SHARED(my_rank, mask)                             &
!$OMP& PRIVATE(ierr, thread_id, core_id, hex_char, isubstr)
  
  thread_id=0
!$ thread_id = omp_get_thread_num()
  core_id = fort2c_affinity_running_on_core()

  !Threads 0-15 are hexadecimal. Above those, continue along the collating
  !sequence, starting from letter G.
  IF(thread_id>15) THEN
    WRITE(hex_char,'(A1)') ACHAR(ICHAR('G')+(thread_id-16))
  ELSE
    WRITE(hex_char,'(Z1)') thread_id
  END IF
  isubstr = core_id + 1

  !If more than one thread is running on the same core, show that with a hash
  !symbol. 
!$OMP CRITICAL
  IF( mask(isubstr:isubstr) == '.' ) THEN
    mask(isubstr:isubstr) = TRIM(hex_char)
  ELSE
    mask(isubstr:isubstr) = '#'
  END IF
!$OMP END CRITICAL

!$OMP END PARALLEL
  
  num_available_cpus = fort2c_affinity_num_available_cpus()
  
  !Everyone sends to writer
  WRITE(send_buffer,'(a,i8,1x,a,1x,a,1x,1a,1x,i0)')                    &
            'PE', my_rank, ':', TRIM(mask), ':', num_available_cpus
  CALL mpl_isend(send_buffer, buffer_size, mpl_character,              &
                 0, my_rank,                                           &
                 communicator, send_request, ierr)

  !Writer receives and writes to file
  IF(on_writer) THEN
    CALL umPrint(blank_string, src='affinity')
    WRITE(umMessage,'(a)') '--> AFFINITY REPORT (long form) <--'
    CALL umPrint(umMessage,src='affinity')
    CALL umPrint(blank_string, src='affinity')

  
    !Thread affinity map
    WRITE(umMessage,'(a)') 'Thread binding map, key:'
    CALL umPrint(umMessage,src='affinity')
    CALL umPrint(blank_string, src='affinity')

    WRITE(umMessage,'(a)')  'PE <RANK>'                    &
                       // ' : ...THREADS..ON..CORES... : ' &
                       // 'Num. cores available to threads'
    CALL umPrint(umMessage,src='affinity')
    WRITE(umMessage,'(a)')  '(# character denotes '        &
                       //   'multiple threads on same core)'
    CALL umPrint(umMessage,src='affinity')
    CALL umPrint(blank_string, src='affinity')

    WRITE(umMessage,'(14x,a)') 'Cores ---->'
    CALL umPrint(umMessage,src='affinity')
    
    DO irank = 0, comm_size-1 
      CALL mpl_recv(recv_buffer, buffer_size, mpl_character,  &
                    irank, irank,                             &
                    communicator, recv_status, ierr)
      CALL umPrint(TRIM(recv_buffer), src='affinity')
    END DO
  
    CALL umPrint(blank_string, src='affinity')
  END IF
  
  !Send completed?
  CALL mpl_wait(send_request, send_status, ierr)

#if ! defined (CRAY_FORTRAN)
  !Tidy up
  DEALLOCATE(recv_buffer)
  DEALLOCATE(send_buffer)
  DEALLOCATE(mask)
#endif

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  RETURN
END SUBROUTINE affinity_write_long

END MODULE affinity_mod

!-------------------------------------------------------------------------------
! Non-linux version
!-------------------------------------------------------------------------------
#else

MODULE affinity_mod
USE mpl
USE umPrintMgr, ONLY: umPrint, umMessage
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
!$ USE omp_lib
IMPLICIT NONE
PRIVATE 

PUBLIC :: affinity_write

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AFFINITY_MOD'

CONTAINS

SUBROUTINE affinity_write(communicator, writer_rank, long_form)
  !-----------------------------------------------------------------------------
  ! SYNOPSIS
  !   call affinity_write(communicator, writer_rank, long_form)
  !
  ! DESCRIPTION
  !   Dummy routine, used on platforms without current support for affinity
  !   reporting.  Reports unavailablity to the user.
  !
  ! ARGUMENTS
  !   * communicator - Dummy MPI communicator.
  !   * writer_rank  - Dummy MPI rank.
  !   * long_form    - Dummy logical.
  !-----------------------------------------------------------------------------
  IMPLICIT NONE

  !Arguments
  INTEGER, INTENT(IN) :: communicator
  INTEGER, INTENT(IN) :: writer_rank
  LOGICAL, INTENT(IN) :: long_form

  INTEGER :: my_rank
  INTEGER :: ierr

  CHARACTER(LEN=*), PARAMETER :: blank_string = ''

  REAL(KIND=jprb) :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='AFFINITY_WRITE'

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

  CALL mpl_comm_rank(communicator, my_rank, ierr)

  IF(my_rank == writer_rank) THEN

    WRITE(umMessage,'(a)') '--> AFFINITY REPORT <--'
    CALL umPrint(blank_string, src='affinity')
    CALL umPrint(umMessage,src='affinity')
    CALL umPrint(blank_string, src='affinity')
    
    WRITE(umMessage,'(a)') 'Unavailable for this platform.'
    CALL umPrint(umMessage,src='affinity')
    CALL umPrint(blank_string, src='affinity')
    
  END IF

  IF(lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

  RETURN
END SUBROUTINE affinity_write

END MODULE affinity_mod
#endif

