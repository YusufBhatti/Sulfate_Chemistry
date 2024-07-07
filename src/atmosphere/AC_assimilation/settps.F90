! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE SETTPS--------------------------------------------------
!
!    Purpose : Sets up the list of AC Observation Types in the
!              order they are to be processed in the assimilation.
!              This routine is called each time the AC Observation
!              files are read in.
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!    Arguments----------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE settps_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETTPS_MOD'

CONTAINS

SUBROUTINE settps (icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx
USE umPrintMgr, ONLY:      &
    umPrint,str,            &
    umMessage
USE ac_control_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER ::    icode
CHARACTER(LEN=errormessagelength) :: cmessage

!    INTENT=OUT--------------------------------------------------------
!     ICODE        : Return Code
!     CMESSAGE     : Reason for failure
!    ------------------------------------------------------------------

!     The variables/arrays set up are :-

!     NACT  : Number of AC Observation Types to be processed.
!     LACT  : List of the Obs Types in the order to be processed.
!     N_GROUPS : Number of groups for processing in AC.
!     GROUP_NO : Group in which each type is to be processed.

!     The order of processing is controlled by the array AC_ORDER.

!     The array AC_OBS_TYPES in the ACP namelist is used to control
!     which observation types are to be processed.


!     Local variables

INTEGER :: this_type,j,jobt,jtype
INTEGER :: last_group,this_group
LOGICAL :: luse

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETTPS'


!     THIS_TYPE  : Obs type in AC_ORDER
!     J          : Loop counter for obs types in LACT
!     JOBT       : Loop counter for obs types in AC_ORDER
!     JTYPE      : Loop counter for obs types in AC_ORDER
!     THIS_GROUP : Indicator from AC_ORDER of grouping (current type)
!     LAST_GROUP : Indicator from AC_ORDER of grouping (previous type)

!     ------------------------------------------------------------------
!     1. Initialise arrays and variables set up by SETTPS
!     ------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nact  = 0
n_groups = 0
last_group =0
DO jobt=1,nobtypmx
  lact (jobt) = 0
  group_no(jobt) = 0
  group_first(jobt) = 0
  group_last(jobt) = 0
END DO
!     ----------------------------------------------------------------
!     2. Set up order of processing in LACT
!     -----------------------------------------------------------------
!     Loop over all AC Obs types known to AC Scheme

DO jtype=1,nobtypmx
  this_type  = MOD(def_ac_order(jtype),1000)
  this_group = (def_ac_order(jtype)-this_type)/1000

  !     This loop determines whether the observation type - THIS_TYPE -
  !     is to be used or not from the AC_OBS_TYPES array.

  IF (this_type >  0) THEN

    !       Use observation type if in namelist array AC_OBS_TYPES

    luse = .FALSE.
    DO jobt=1,nobtypmx
      IF (this_type == ac_obs_types(jobt)) THEN
        luse = .TRUE.
      END IF
    END DO

    IF (luse) THEN

      !         Set up to process this observation type
      nact = nact+1
      lact(nact) = this_type

      !         Group observation types ; Set up GROUP_NO and N_GROUPS

      IF (nact == 1 .OR. this_group /= last_group) THEN

        !           Start a new group.
        n_groups = n_groups+1
        group_index(n_groups) = this_group
        group_first(n_groups) = nact

      END IF
      last_group = this_group
      group_no(nact) = n_groups
      group_last(n_groups) = nact

      !         Find this type in MASTER_AC_TYPES ; Abort if not found.
      type_index(nact)=0
      DO jobt=1,nobtypmx
        IF (this_type == master_ac_types(jobt)) THEN
          type_index(nact)=jobt
        END IF
      END DO
      IF (type_index(nact) == 0) THEN
        icode = 1
        cmessage = 'SETTPS : Observation Type not in Master List ?'
        CALL umPrint(' Observation Type '//TRIM(str(this_type))//          &
            ' not in Master List ?',src='settps',pe=0)
        GO TO 999
      END IF

    END IF

  END IF
END DO   !   End of JTYPE loop.

!     -----------------------------------------------------------------
!     3. Print out list of AC Obs types to be processed
!     -----------------------------------------------------------------

IF (nact >  0 .AND. mype == 0) THEN

  CALL umPrint(' ',src='settps')
  CALL umPrint(' AC Obs Types to be processed this run',src='settps')
  CALL umPrint(' ',src='settps')
  WRITE(umMessage,'(A,(T12,13I5))') ' Type  No ',                        &
      (lact(j),j=1,nact)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,(T12,13I5))') ' Group No ',                        &
      (group_no(j),j=1,nact)
  CALL umPrint(umMessage,src='settps')
  !       WRITE(umMessage,'(A,15I5)') ' Position in Obs Type List    ',
  !    +  (TYPE_INDEX(J),J=1,NACT)
  !        CALL umPrint(umMessage,src='settps')

  CALL umPrint(' ',src='settps')
  WRITE(umMessage,'(A,15I5)') ' Group Number                ',           &
      (j,j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Group Index                 ',           &
      (group_index(j),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' First Type in Group         ',           &
      (group_first(j),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Last Type in Group          ',           &
      (group_last (j),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' No of iterations            ',           &
      (def_no_iterations(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Interval between Iterations ',           &
      (def_interval_iter(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Ratio of MG Rows to AG Rows ',           &
      (def_agres_rows(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Ratio of MG Pts  to AG Pts  ',           &
      (def_agres_pts(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' No of analysis levels       ',           &
      (def_no_anal_levs(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' No of weight levels         ',           &
      (def_no_wt_levs(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15I5)') ' Horizontal Analysis Mode    ',           &
      (def_mode_hanal(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  CALL umPrint(' ',src='settps')
  CALL umPrint(' Group Dep scaling FACTORS in FI  ',src='settps')
  WRITE(umMessage,'(A,15I11)')' Group No ',(j,j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  WRITE(umMessage,'(A,15E11.4)') '          ',                           &
      (def_fi_var_factor(group_index(j)),j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  CALL umPrint(' ',src='settps')
  CALL umPrint(' Nudging Coefficients ',src='settps')
  WRITE(umMessage,'(A,15I11)')' Group No ',(j,j=1,n_groups)
  CALL umPrint(umMessage,src='settps')
  IF (model_type == mt_global) THEN
    WRITE(umMessage,'(A,15E11.4)') ' NH       ',                         &
        (def_nudge_nh(group_index(j)),j=1,n_groups)
    CALL umPrint(umMessage,src='settps')
    WRITE(umMessage,'(A,15E11.4)') ' TR       ',                         &
        (def_nudge_tr(group_index(j)),j=1,n_groups)
    CALL umPrint(umMessage,src='settps')
    WRITE(umMessage,'(A,15E11.4)') ' SH       ',                         &
        (def_nudge_sh(group_index(j)),j=1,n_groups)
    CALL umPrint(umMessage,src='settps')
  ELSE
    WRITE(umMessage,'(A,15E11.4)') '          ',                         &
        (def_nudge_lam(group_index(j)),j=1,n_groups)
    CALL umPrint(umMessage,src='settps')
  END IF

ELSE IF (nact == 0 .AND. mype == 0) THEN

  CALL umPrint(' SETTPS : No observation types to process ?',src='settps')
  icode = 1
  cmessage = 'SETTPS : No obs types to process ?'
  GO TO 999

END IF

999   CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE settps
END MODULE settps_mod
