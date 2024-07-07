! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Provides convenience routines for async stash operations

MODULE stwork_aux

USE IOS_Stash
USE IOS_Stash_Common
USE IOS_Constants
USE IOS_Model_Geometry, ONLY:                                               &
    maxFieldDomain
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
CHARACTER (LEN=errormessagelength) :: stwork_aux_Message

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STWORK_AUX'

CONTAINS

SUBROUTINE completedBuffer(s)

IMPLICIT NONE
INTEGER, INTENT(IN) :: s
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPLETEDBUFFER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Release data buffers to the IO server
CALL ios_stash_dispatch( s )
! Send the control message to the IO Server.
IF (model_rank == 0)                                                       &
    CALL ios_stash_write_pp_data( s )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE completedBuffer


SUBROUTINE flushPendingStash()
IMPLICIT NONE
INTEGER         :: i
INTEGER         :: state
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLUSHPENDINGSTASH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (isUsingAsyncStash()) THEN
  DO i=1,IOS_AsyncNumSlots
    state=IOS_getSlotState(i)

    IF (slot(i)%bufferType==IOS_Action_StashWritePPData) THEN
      SELECT CASE (state)

      CASE (ios_queue_slot_partfilled)
        CALL completedBuffer(i)
      CASE (ios_queue_slot_initialized)
        IF (IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_Message,'(A,I3,A)')                                    &
              'stwork_aux: flush ',i,' initialised (resetting)'
!$OMP END CRITICAL(internal_write)
          CALL IOS_Print(IOS_Message,src='stwork_aux')
        END IF
        CALL resetBuffer(i)
      CASE (ios_queue_slot_dispatched)

      CASE (ios_queue_slot_unused)

      CASE (ios_queue_slot_ready)

      CASE DEFAULT
!$OMP CRITICAL(internal_write)
        WRITE(stwork_aux_message,'(A,I8)')                                 &
            'Slot :',i,' is in unknown state ',state
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( 'stwork_aux:flushPending',                       &
            -10, stwork_aux_Message )
      END SELECT
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE flushPendingStash


SUBROUTINE flushPendingStashForUnit(UNIT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER             :: i
INTEGER             :: state
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLUSHPENDINGSTASHFORUNIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (isUsingAsyncStash()) THEN
  DO i=1,IOS_AsyncNumSlots
    IF (UNIT==IOS_getSlotUnit(i)) THEN
      state=IOS_getSlotState(i)
      SELECT CASE (state)

      CASE (ios_queue_slot_partfilled)
        CALL completedBuffer(i)
      CASE (ios_queue_slot_initialized)
        IF (IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
!$OMP CRITICAL(internal_write)
          WRITE(IOS_Message,'(A,I3,A)')'stwork_aux: flush ',i,             &
              ' initialised (resetting)'
!$OMP END CRITICAL(internal_write)
          CALL IOS_Print(IOS_Message,src='stwork_aux')
        END IF
        CALL resetBuffer(i)
      CASE (ios_queue_slot_dispatched)

      CASE (ios_queue_slot_unused)

      CASE (ios_queue_slot_ready)

      CASE DEFAULT
!$OMP CRITICAL(internal_write)
        WRITE(stwork_aux_message,'(A,I8)')                                 &
            'Slot :',i,' is in unknown state ',state
!$OMP END CRITICAL(internal_write)
        CALL IOS_Ereport( 'stwork_aux:flushPending',                       &
            -10, stwork_aux_Message )
      END SELECT
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE flushPendingStashForUnit


FUNCTION getSlotForNextLevel(UNIT) RESULT (s)
IMPLICIT NONE
INTEGER, INTENT(IN) :: UNIT
INTEGER             :: s
INTEGER             :: current_levels_packed
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETSLOTFORNEXTLEVEL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

s=-1

IF (isUsingAsyncStash()) THEN

  s=IOS_getCurrentSlot(UNIT)
  IF (s<0) THEN

    CALL ios_stash_next_buffer( UNIT, s )
    IF (s==-1) THEN
      IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
        CALL IOS_Print( 'Warning: stwork_aux: No new slot returned',       &
            src='stwork_Aux')
        WRITE(IOS_Message,'(A,A)')'Warning: stwork_aux: ',                 &
            '  flushing everything - expensive '
        CALL IOS_print(IOS_message,src='stwork_aux')
      END IF

      ! next buffer failed, flush stuff and retry...
      CALL flushPendingStash()
      CALL ios_stash_next_buffer( UNIT, s )

      IF (s==-1) THEN
        WRITE(stwork_aux_message,'(A)')                                    &
            'Failure tring to locate a buffer (code fault)'
        CALL IOS_Ereport( 'stwork_aux:slotfornextlevel',                   &
            10, stwork_aux_Message )
      END IF

    END IF

    IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
!$OMP CRITICAL(internal_write)
      WRITE(IOS_Message,'(A,A,I3,A,I3)')'Info: stwork_aux: ',              &
          'No current slot for unit, ',UNIT,                               &
          ' new slot is ',s
!$OMP END CRITICAL(internal_write)
      CALL IOS_Print(IOS_Message,src='stwork_aux')
    END IF

  ELSE
    current_levels_packed=IOS_getLevsInPack(s)
    IF (current_levels_packed >= IOS_AsyncMaxFieldsInPack ) THEN
      ! This buffer is full need a new buffer

      ! Dispatch the old one
      CALL completedBuffer(s)

      !Get a new one
      CALL ios_stash_next_buffer( UNIT, s )

      IF (s==-1) THEN
        WRITE(stwork_aux_message,'(A)')                                    &
            'Failure trying to locate a buffer (code fault)'
        CALL IOS_Ereport( 'stwork_aux:slotfornextlevel',                   &
            10, stwork_aux_Message )
      END IF

      IF (IOS_Verbosity>=IOS_prstatus_diag) THEN
!$OMP CRITICAL(internal_write)
        WRITE(IOS_Message,'(A,A,I3,A,I3)')'Info: stwork_aux: ',            &
            'Buffer for unit, ',UNIT,                                      &
            ' was full, new slot is ',s
!$OMP END CRITICAL(internal_write)
        CALL IOS_Print(IOS_Message,src='stwork_aux')
      END IF

    END IF
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION getSlotForNextLevel

END MODULE stwork_aux
