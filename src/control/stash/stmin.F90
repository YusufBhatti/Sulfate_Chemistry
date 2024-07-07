! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: STMIN ----------------------------------------------------
!
!    Purpose: Computes the point-by-point minimum in time of a field
!             by comparing the field at the current time with the
!             minimum so far (STASH TEMPORAL service routine)
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE stmin(fieldin,RESULT,fsize,masking)
!
USE missing_data_mod, ONLY: rmdi
USE stparam_mod, ONLY: s_lndms, s_seams, s_nmdims
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
INTEGER :: fsize            ! IN size of fieldin and result.
REAL :: fieldin(fsize)      ! IN input field
REAL :: RESULT(fsize)       ! OUT output field (minimum)
INTEGER :: masking          ! IN controls if masked or MDI treatment
! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER :: i ! loop count

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STMIN'
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!
!  1.1 loop over array size, if either result or fieldin is MDI set
!      result to MDI, else set result to minimum of fieldin and result
!
IF (masking==s_lndms .OR. masking==s_seams) THEN
  DO i=1,fsize
    IF ((RESULT(i) /= rmdi) .AND. (fieldin(i) /= rmdi)) THEN
      RESULT(i)=MIN(RESULT(i),fieldin(i))
    ELSE
      RESULT(i)=rmdi
    END IF
  END DO

ELSE IF (masking==s_nmdims) THEN
  !
  !  1.2 loop over array size, set result to minimum of fieldin and result
  !      while skipping missing data
  !
  DO i = 1,fsize
    IF ((RESULT(i) /= rmdi) .AND. (fieldin(i) /= rmdi)) THEN
      RESULT(i) = MIN(RESULT(i),fieldin(i))
    ELSE IF (fieldin(i) /= rmdi) THEN
      RESULT(i) = fieldin(i)
    END IF
  END DO

ELSE
  !
  !  1.3 loop over array size, set result to minimum of fieldin and result
  !      without checking for missing data
  !
  DO i=1,fsize
    RESULT(i)=MIN(RESULT(i),fieldin(i))
  END DO
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stmin
