! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE AC_STASH -----------------------------------------------
!
!    Purpose : Stash processing for AC Scheme
!
!    Global and Limited area
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE ac_stash_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='AC_STASH_MOD'

CONTAINS

SUBROUTINE ac_stash (stash_item_no,level_no,                      &
                     group_no,n_groups,timestep_no,               &
                     stindex,stlist,len_stlist,si,sf,             &
                     stashwork,stash_levels,num_stash_levels,     &
                     stash_pseudo_levels,num_stash_pseudo,        &
                     field,len_fld,label,                         &
                     icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE stparam_mod, ONLY: st_input_bottom, st_pseudo_in
USE nlsizes_namelist_mod, ONLY: model_levels

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


INTEGER ::                                                        &
   stash_item_no                                                  &
                        ! Stash Item Number
,  level_no                                                       &
                        ! Level number
,  group_no                                                       &
                        ! Group number
,  n_groups                                                       &
                        ! No of groups
,  timestep_no                                                    &
                        ! Timestep number
,  len_fld                                                        &
                        ! Length of field to be stashed
,  len_stlist                                                     &
                        ! Dimension of STLIST
,  stindex(2,*)                                                   &
                        ! Start and no of items in STLIST
,  stlist(len_stlist,*)                                           &
                        ! Stash List of items to be output
,  si(*)                                                          &
                        ! Address of Item in STASHWORK
,  num_stash_levels                                               &
                        ! Number of levels lists
,  num_stash_pseudo                                               &
                        ! Number of pseudo lists
,  stash_levels(num_stash_levels+1,*)                             &
                                             ! Levels lists
,  stash_pseudo_levels(num_stash_pseudo+1,*) ! Pseudo lists

REAL ::                                                           &
   stashwork(*)                                                   &
                        ! Work array for stashed data
,  field(len_fld)       ! Field to be stashed

LOGICAL :: sf(*)           ! Stash Flags

CHARACTER(LEN=*) :: label     !  Label to indicate field being stashed

INTEGER :: icode           !  Return code
CHARACTER(LEN=errormessagelength) :: cmessage  !  Error message

!     Dynamic allocated arrays

LOGICAL :: levels_list(model_levels)  ! Expanded levels list
LOGICAL :: pseudo_list(n_groups)      ! Expanded pseudo list

!     Local variables

LOGICAL ::                                                        &
   l_single_lev                                                   &
,  l_levels_list                                                  &
,  l_pseudo_list                                                  &
,  lstash

INTEGER ::                                                        &
   j,jlev,jgrp                                                    &
                    !  Loop counters over levels/groups
,  j0                                                             &
                    !  Pointer in STASHWORK
,  ipos                                                           &
                    !  Position in STLIST for this STASH_ITEM_NO
,  n_levels_list                                                  &
                    !  No of levels in levels list
,  fld_no                                                         &
                    !  Field number in STASHWORK
,  grp_no           !  Group number in STASHWORK

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='AC_STASH'


!     Check that LEVEL_NO le N_LEVELS

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (level_no >  model_levels) THEN
  icode    = 1
  cmessage = ' AC_STASH : LEVEL_NO gt N_LEVELS ?'
  WRITE(umMessage,*) 'AC_STASH : LEVEL_NO must be LE to N_LEVELS'
  CALL umPrint(umMessage,src='ac_stash')
  WRITE(umMessage,*) 'Stash Item No = ',stash_item_no
  CALL umPrint(umMessage,src='ac_stash')
  WRITE(umMessage,*) 'LEVEL_NO = ',level_no,' MODEL_LEVELS = ',model_levels
  CALL umPrint(umMessage,src='ac_stash')
  GO TO 999   !  Return
END IF

!     Check that GROUP_NO le N_GROUPS
IF (group_no >  n_groups) THEN
  icode    = 1
  cmessage = ' AC_STASH : GROUP_NO gt N_GROUPS ?'
  WRITE(umMessage,*) 'AC_STASH : GROUP_NO must be LE to N_GROUPS'
  CALL umPrint(umMessage,src='ac_stash')
  WRITE(umMessage,*) 'Stash Item No = ',stash_item_no
  CALL umPrint(umMessage,src='ac_stash')
  WRITE(umMessage,*) 'GROUP_NO = ',group_no,' N_GROUPS = ',n_groups
  CALL umPrint(umMessage,src='ac_stash')
  GO TO 999   !  Return
END IF

DO jlev = 1,model_levels
  levels_list(jlev) = .FALSE.
END DO
DO jgrp = 1,n_groups
  pseudo_list(jgrp) = .FALSE.
END DO

!     Get position in STLIST
ipos = stindex(1,stash_item_no)

!     Determine if levels list used (Entry 10 in STLIST = negative)
l_levels_list = stlist(st_input_bottom,ipos) <  0

!     Determine if single level (Entry 10 in STLIST = 100)
l_single_lev  = stlist(st_input_bottom,ipos) == 100

!     Determine if pseudo list used (Entry 26 in STLIST = positive)
l_pseudo_list = stlist(st_pseudo_in,ipos) >  0

n_levels_list = 1

IF (l_levels_list) THEN

  n_levels_list = stash_levels(1,-stlist(st_input_bottom,ipos))

  !----   Get levels required for this field
  !----   (Sets up LEVELS_LIST from STASH_LEVELS)
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list (model_levels,len_stlist,stlist(1,ipos),       &
       levels_list,stash_levels,num_stash_levels+1,                   &
       icode,cmessage)
  IF (icode >  0) GO TO 999   !  Return

ELSE IF (.NOT. l_single_lev) THEN

  icode = 1
  cmessage = 'AC_STASH ; No levels list ?'
  GO TO 999   !  Return

ELSE
END IF

IF (l_pseudo_list) THEN

  !----   Get pseudo levels list required for this field
  !----   (Sets up PSEUDO_LIST from STASH_PSEUDO_LEVELS)
  ! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list (n_groups,len_stlist,stlist(1,ipos),       &
       pseudo_list,stash_pseudo_levels,num_stash_pseudo,          &
       icode,cmessage)
  IF (icode >  0) GO TO 999   !  Return

END IF

!     Determine if this field is to stashed
lstash = .FALSE.
IF ( l_single_lev .OR.                                            &
    (l_levels_list .AND. levels_list(level_no)) ) THEN
  IF (l_pseudo_list) THEN
    IF (pseudo_list(group_no)) lstash = .TRUE.
  ELSE
    lstash = .TRUE.
  END IF
END IF

IF (lstash) THEN

  !----   Determine position in STASHWORK for this field
  fld_no = 0
  IF (l_single_lev) THEN
    fld_no = 1
    !         Could have pseudo levels - look at later
  ELSE
    DO jlev=1,level_no
      IF (levels_list(jlev)) THEN
        fld_no = fld_no+1
      END IF
    END DO
  END IF

  grp_no = 0
  IF (l_pseudo_list) THEN
    DO jlev=1,group_no
      IF (pseudo_list(jlev)) THEN
        grp_no = grp_no+1
      END IF
    END DO
  ELSE
    grp_no = 1
  END IF

  !----   Check FLD_NO
  IF (fld_no == 0) THEN
    icode    = 1
    cmessage = ' AC_STASH : FLD_NO = 0 ?'
    WRITE(umMessage,*) 'AC_STASH : FLD_NO must be GT than 0'
    CALL umPrint(umMessage,src='ac_stash')
    WRITE(umMessage,*) 'Stash Item No = ',stash_item_no
    CALL umPrint(umMessage,src='ac_stash')
    WRITE(umMessage,*) 'Level,Group,FLD_NO = ',level_no,group_no,fld_no
    CALL umPrint(umMessage,src='ac_stash')
    GO TO 999   !  Return
  END IF

  !----   Check GRP_NO
  IF (grp_no == 0) THEN
    icode    = 1
    cmessage = ' AC_STASH : GRP_NO = 0 ?'
    WRITE(umMessage,*) 'AC_STASH : GRP_NO must be GT than 0'
    CALL umPrint(umMessage,src='ac_stash')
    WRITE(umMessage,*) 'Stash Item No = ',stash_item_no
    CALL umPrint(umMessage,src='ac_stash')
    WRITE(umMessage,*) 'Level,Group,FLD_NO = ',level_no,group_no,fld_no
    CALL umPrint(umMessage,src='ac_stash')
    GO TO 999   !  Return
  END IF

  !----   Set up pointer for this field in STASHWORK
  j0 = si(stash_item_no) - 1                                      &
     + (grp_no-1)*(n_levels_list*len_fld)                         &
     + (fld_no-1)*len_fld

  !----   Copy field into work space
  DO j=1,len_fld
    stashwork(j0 + j) = field (j)
  END DO


END IF

999  CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ac_stash
END MODULE ac_stash_mod
