! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: STACCUM  -------------------------------------------------
!
!    Purpose: Accumulates fields within STASH (temporal service routine)
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE staccum(fieldin,RESULT,fsize,masking)
!
USE missing_data_mod, ONLY: rmdi
USE stparam_mod, ONLY: s_nmdims, s_seams, s_lndms
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
INTEGER :: fsize             ! IN size of fieldin and result.
REAL :: fieldin(fsize)       ! IN  input field to be processed
REAL :: RESULT(fsize)        ! OUT where accum is done.
INTEGER :: masking           ! IN controls if masked or MDI treatment
!  ---------------------------------------------------------------------
!
!  Local variables
!
INTEGER :: i                 ! loop count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STACCUM'
! ----------------------------------------------------------------------
!  Loop over array size, masking to allow different treatment of MDI.
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!
!  BEST - If neither result or fieldin is rmdi accumulate in result,
!         If either is not rmdi, make use of value,
!         If both result and fieldin are rmdi, leave as rmdi,
!         as accumulating rmdi does not really make sense.
IF (masking == s_nmdims) THEN
  DO i=1,fsize
    IF ((RESULT(i) /= rmdi) .AND. (fieldin(i) /= rmdi)) THEN
      RESULT(i) = RESULT(i)+fieldin(i)
    ELSE IF ((fieldin(i) /= rmdi)) THEN
      RESULT(i) = fieldin(i)
    ELSE IF ((RESULT(i) /= rmdi)) THEN
      RESULT(i) = RESULT(i)
    END IF
  END DO

  !  If standard fixed sea or land point masking, only accumulate non-MDI
ELSE IF (masking==s_seams .OR. masking==s_lndms) THEN
  DO i=1,fsize
    IF ((RESULT(i) /= rmdi) .AND. (fieldin(i) /= rmdi)) THEN
      RESULT(i)=RESULT(i)+fieldin(i)
    ELSE
      RESULT(i)=rmdi
    END IF
  END DO

ELSE
  !  ORIGINAL CODE - accumulate result without checking for missing data
  !  but NOT very useful as accumulating MDI values does not make sense.
  DO i=1,fsize
    RESULT(i)=RESULT(i)+fieldin(i)
  END DO
END IF
!
!  ALTERNATIVE CODE - suggested but not used.
!  If neither result or fieldin is MDI accumulate in result,
!  If either is not MDI, make use of value,
!  If both result and fieldin are MDI and not masking, 
!  reset to zero, as accumulating MDI values does not make sense.
!      ELSE
!        DO i=1,fsize
!          IF ((result(i) /= rmdi).AND.(fieldin(i) /= rmdi)) THEN
!            result(i) = result(i)+fieldin(i)
!          ELSE IF ((fieldin(i) /= rmdi)) THEN
!            result(i) = fieldin(i)
!          ELSE IF ((result(i) /= rmdi)) THEN
!            result(i) = result(i)
!          ELSE
!            result(i) = 0.0
!          END IF
!        END DO
!      END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE staccum
