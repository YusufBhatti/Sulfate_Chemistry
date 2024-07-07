! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_ZERO_LEVELS_LIST
!
! Purpose : To set up a list of levels at which a diagnostic is
!           required, including the zeroth level, using information
!           in the STASH list.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE set_zero_levels_list(levels, len_stlist, stlist, list, &
  stash_levels, len_stashlevels, icode, cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER ::                                                        &
  levels,                                                         &
               ! IN Number of levels in input data
  len_stlist,                                                     &
               ! IN
  stlist(len_stlist),                                             &
               ! IN STASH list
  len_stashlevels,                                                &
  stash_levels(len_stashlevels,*),                                &
               ! IN - list of levels required
  icode        ! OUT Return code =0 Normal exit
               !                 >1 Error message

LOGICAL ::                                                        &
  list(0:levels) ! OUT List of levels required.

CHARACTER(LEN=errormessagelength) :: cmessage ! Error message

!  Local variables
INTEGER :: k, kout

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_ZERO_LEVELS_LIST'


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  Initialise levels list to false
DO k=0, levels
  list(k)= .FALSE.
END DO

!  Check for method of levels selection
!  Levels list must be present.

IF (stlist(10) <  0) THEN

  ! Set logical array list to identify levels required.
  DO kout=2,stash_levels(1,-stlist(10))+1
    IF ((stash_levels(kout,-stlist(10)) >= 0) .AND. &
    (stash_levels(kout,-stlist(10)) <= levels)) THEN
      ! LEVEL IS IN THE RANGE OF LIST.
      list(stash_levels(kout,-stlist(10))) =.TRUE.
    ELSE
      ! LEVEL IS OUT OF THE RANGE OF LIST.
      cmessage=  ' SET_ZERO_LEVELS_LIST: level out of range'
      WRITE(umMessage,'(A)') ' SET_ZERO_LEVELS_LIST: level out of range'
      CALL umPrint(umMessage,src='set_zero_levels_list')
      WRITE(umMessage,'(A,I0)') ' level=',stash_levels(kout,-stlist(10))
      CALL umPrint(umMessage,src='set_zero_levels_list')
      WRITE(umMessage,'(A,2I0)') ' Section, Item =',stlist(2),stlist(1)
      CALL umPrint(umMessage,src='set_zero_levels_list')
      icode=2
    END IF
  END DO

ELSE IF (stlist(10) /= 100) THEN

  !  Set list of levels according to its definition in stlist :
  !l If stlist(10) positive, input on range of model levels starting
  !l at level stlist(10), finishing at level stlist(11)
  DO k=stlist(10),stlist(11)
    list(k)= .TRUE.
  END DO

ELSE

  !  Illegal control data
  icode=1
  cmessage='SET_ZERO_LEVELS_LIST: Illegal control data'
  WRITE(umMessage,'(A,2I0)') &
    'Illegal control data SET_ZERO_LEVELS_LIST,STLIST(10,11)=', &
    stlist(10) ,stlist(11)
  CALL umPrint(umMessage,src='set_zero_levels_list')
  WRITE(umMessage,'(A,2I0)') 'Section and item numbers ',stlist(2),stlist(1)
  CALL umPrint(umMessage,src='set_zero_levels_list')

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_zero_levels_list
