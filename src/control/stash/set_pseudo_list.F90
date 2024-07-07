! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_PSEUDO_LIST -----------------------------------------
!
! Purpose : To set up a list of pseudo levels at which a diagnostic
!           is required, using information in the STASH list.
!
! Copy of Subroutine SET_LEVELS_LIST (Deck SETLST1) taken and
! adapted for pseudo levels.
!
! Programming Standard : Unified Model Documentation paper number 3
!
! Documentation: U.M. Documentation paper number C4
!
! -----------------------------------------------------------------
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE set_pseudo_list                                        &
      (n_levels,len_stlist,stlist,pseudo_list,                    &
      stash_pseudo_levels,num_stash_pseudo,icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stparam_mod, ONLY: st_pseudo_in
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER ::                                                        &
   n_levels                                                       &
                       ! IN Number of possible pseudo levels
,  len_stlist                                                     &
                       ! IN Dimension of STLIST
,  stlist(len_stlist)                                             &
                       ! IN STASH list
,  num_stash_pseudo                                               &
                       ! IN Dimension for STASH_PSEUDO_LEVELS
,  stash_pseudo_levels(num_stash_pseudo+1,*)                      &
                                             ! IN Pseudo levels
,  icode               ! OUT Return code

LOGICAL ::                                                        &
   pseudo_list(n_levels) ! OUT List of pseudo levels required.

CHARACTER(LEN=errormessagelength) :: cmessage ! Error message

!  ---------------------------------------------------------------------

!  Local variables

INTEGER ::                                                        &
      jlev                                                        &
                 ! Loop counter over levels
,     level_no                                                    &
                 ! Level no in pseudo list
,     list_no    ! Pseudo level list number

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_PSEUDO_LIST'

!  ---------------------------------------------------------------------

!  Initialise pseudo levels list to false

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO jlev=1,n_levels
  pseudo_list(jlev)= .FALSE.
END DO

!  Get pseudo list number

list_no = stlist(st_pseudo_in)

!  Check that Pseudo list number is valid (should be GE 0)

IF (list_no <  0) THEN

  !       Illegal control data

  icode=1
  cmessage = 'SET_PSEUDO_LIST: Illegal control data'
  WRITE(umMessage,*) 'SET_PSEUDO_LIST: Illegal control data'
  CALL umPrint(umMessage,src='set_pseudo_list')
  WRITE(umMessage,*) 'ST_PSEUDO_IN         = ',st_pseudo_in
  CALL umPrint(umMessage,src='set_pseudo_list')
  WRITE(umMessage,*) 'STLIST(ST_PSEUDO_IN) = ',stlist(st_pseudo_in)
  CALL umPrint(umMessage,src='set_pseudo_list')
  WRITE(umMessage,*) 'Section and item numbers ',stlist(2),stlist(1)
  CALL umPrint(umMessage,src='set_pseudo_list')
  GO TO 999  !  Return

END IF

!  Set logical array list to identify pseudo levels required.

IF (list_no >  0) THEN

  DO jlev=2,stash_pseudo_levels(1,list_no)+1
    level_no = stash_pseudo_levels(jlev,list_no)
    IF (level_no >= 1 .AND. level_no <= n_levels) THEN

      !           Level is within range
      pseudo_list(level_no) =.TRUE.

    ELSE

      !           Level is out of range
      icode=2
      cmessage=  ' SET_PSEUDO_LIST : level out of range'
      WRITE(umMessage,*) ' SET_PSEUDO_LIST : level out of range'
      CALL umPrint(umMessage,src='set_pseudo_list')
      WRITE(umMessage,*) ' pseudo list no = ',list_no
      CALL umPrint(umMessage,src='set_pseudo_list')
      WRITE(umMessage,*) ' level = ',level_no
      CALL umPrint(umMessage,src='set_pseudo_list')
      WRITE(umMessage,*) ' Section, Item = ',stlist(2),stlist(1)
      CALL umPrint(umMessage,src='set_pseudo_list')
      GO TO 999   !  Return

    END IF
  END DO

END IF

999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_pseudo_list
