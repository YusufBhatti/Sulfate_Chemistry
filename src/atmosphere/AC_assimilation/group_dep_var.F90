! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE GROUP_DEP_VAR ------------------------------------------
!
!    Purpose : Process &ACP Namelist arrays which are Group Dependent.
!
!
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE group_dep_var_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GROUP_DEP_VAR_MOD'

CONTAINS

SUBROUTINE group_dep_var (ac_order,no_iterations,interval_iter,   &
                          n_anal_levs,n_wt_levs,                  &
                          nudge_nh,nudge_tr,nudge_sh,             &
                          nudge_lam,                              &
                          agres_rows,agres_pts,                   &
                          mode_hanal,fi_var_factor,               &
                          icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: imdi, rmdi 
USE ac_control_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


INTEGER :: ac_order     (nobtypmx)  !IN groupings
INTEGER :: no_iterations(nobtypmx)  !IN iterations per step
INTEGER :: interval_iter(nobtypmx)  !IN steps between AC
INTEGER :: n_anal_levs  (nobtypmx)  !IN analysis levels
INTEGER :: n_wt_levs    (nobtypmx)  !IN weights levels
INTEGER :: agres_rows   (nobtypmx)  !IN ratio modelrows/ anl rows
INTEGER :: agres_pts    (nobtypmx)  !IN ratio model pts/anl pts
INTEGER :: mode_hanal   (nobtypmx)  !IN mode of horizonatl analysis
REAL :: fi_var_factor(nobtypmx)  !IN group dep scaling in FI
REAL :: nudge_nh(nobtypmx)       !IN Nudging coeff NH
REAL :: nudge_tr(nobtypmx)       !IN Nudging coeff TR
REAL :: nudge_sh(nobtypmx)       !IN Nudging coeff SH
REAL :: nudge_lam(nobtypmx)      !IN Nudging coeff LAM
INTEGER :: icode                    !OUT error code and message
CHARACTER(LEN=errormessagelength) :: cmessage

! ---------------------------------------------------------------------
!  Local work space and variables
INTEGER :: jobt,jobt2,jg,n_obtyp,last_group,first_type,last_type
INTEGER :: this_type,this_type_def
INTEGER :: this_group,no_groups,ncount
INTEGER :: new_no_iters(nobtypmx),new_int_iter(nobtypmx)
INTEGER :: new_agres_rows(nobtypmx),new_agres_pts(nobtypmx)
INTEGER :: new_anal_levs(nobtypmx),new_wt_levs(nobtypmx)
INTEGER :: new_mode_hanal(nobtypmx)
LOGICAL :: l_new_groups,l_no_iters,l_int_iter
LOGICAL :: l_agres_rows,l_agres_pts,l_anal_levs,l_wt_levs
LOGICAL :: l_mode_hanal,l_fi_var_factor
REAL :: new_fi_var_factor(nobtypmx)
REAL :: new_nudge_nh(nobtypmx),new_nudge_tr(nobtypmx)
REAL :: new_nudge_sh(nobtypmx)
LOGICAL :: l_nudge_nh,l_nudge_tr,l_nudge_sh
REAL :: new_nudge_lam(nobtypmx)
LOGICAL :: l_nudge_lam
LOGICAL :: lknown,lfound
LOGICAL :: l_dummy  ! dummy logical

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GROUP_DEP_VAR'

! ---------------------------------------------------------------------

!     Check and count values in namelist AC_ORDER
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
n_obtyp=0
DO jobt=1,nobtypmx
  IF (ac_order(jobt) >  0) THEN
    n_obtyp = n_obtyp+1
  ELSE IF (ac_order(jobt) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Negative Value in AC_ORDER ?'
    GO TO 999
  ELSE
  END IF
END DO

l_new_groups = n_obtyp >  0
IF (l_new_groups) THEN

  !       Check validity of obs types in namelist AC_ORDER.
  DO jobt=1,n_obtyp
    IF (ac_order(jobt) >  0) THEN
      this_type = MOD(ac_order(jobt),1000)
      IF (this_type == 501) THEN
        this_type = 302
        this_group = ac_order(jobt)/1000
        ac_order(jobt) = this_group*1000 + this_type
        CALL umPrint('Type 501 in AC_ORDER changes to Type 302', &
            src='group_dep_var',pe=0)
      END IF
      IF (this_type == 502) THEN
        this_type = 305
        this_group = ac_order(jobt)/1000
        ac_order(jobt) = this_group*1000 + this_type
        CALL umPrint('Type 502 in AC_ORDER changes to Type 305', &
            src='group_dep_var',pe=0)
      END IF
      lknown  = .FALSE.
      DO jobt2 = 1,nobtypmx
        this_type_def = MOD(def_ac_order(jobt2),1000)
        IF (this_type == this_type_def) THEN
          lknown = .TRUE.
        END IF
      END DO
      IF (.NOT. lknown) THEN
        icode = 1
        cmessage =                                                &
        'GROUPDEP : Obs Type in AC_ORDER not known.'
        GO TO 999
      END IF
    END IF
  END DO

  !       Check that all obs types in AC_OBS_TYPES are in AC_ORDER.
  DO jobt=1,nobtypmx
    IF (ac_obs_types(jobt) >  0) THEN
      lfound = .FALSE.
      DO jobt2=1,n_obtyp
        this_type = MOD(ac_order(jobt2),1000)
        IF (this_type == ac_obs_types(jobt)) THEN
          lfound = .TRUE.
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage ='GROUPDEP : AC_OBS_TYPES / AC_ORDER mismatch'
        GO TO 999
      END IF
    END IF
  END DO

  !       Check validity and order of group numbers in AC_ORDER.
  last_group = 0
  DO jobt =1,nobtypmx
    IF (ac_order(jobt) >  0) THEN
      this_group = ac_order(jobt)/1000
      IF (this_group == 0) THEN
        icode = 1
        cmessage =                                                &
        'GROUPDEP : Obs Type in AC_ORDER with no group number.'
        GO TO 999
      END IF
      IF (this_group /= last_group .AND.                          &
          this_group /= last_group+1 ) THEN
        icode = 1
        cmessage =                                                &
        'GROUPDEP : Order of groups in AC_ORDER incorrect.'
        GO TO 999
      END IF
      last_group = this_group
    END IF
  END DO

ELSE   !  AC_ORDER not used ; Get no of obs types.

  n_obtyp = 0
  DO jobt=1,nobtypmx
    ac_order(jobt) = def_ac_order(jobt)
    IF (ac_order(jobt) >  0) THEN
      n_obtyp = n_obtyp+1
    END IF
  END DO

END IF

last_group = 0
no_groups  = 0
DO jobt=1, nobtypmx
  IF (ac_order(jobt) >  0) THEN
    this_group = ac_order(jobt)/1000
    IF (this_group /= last_group) THEN
      no_groups = no_groups+1
    END IF
    last_group = this_group
  END IF
END DO

IF (nprog == 1001 .AND. mype == 0) THEN

  WRITE(umMessage,'(10I7)') (ac_order(jobt),jobt=1,nobtypmx)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'N_OBTYP = ',n_obtyp
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'NO_ITERATIONS'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (no_iterations(jobt),jobt=1,nobtypmx)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'INTERVAL_ITER'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (interval_iter(jobt),jobt=1,nobtypmx)
  CALL umPrint(umMessage,src='group_dep_var')
END IF

!     Process NO_ITERATIONS
!     ---------------------
ncount = 0
DO jg=1,nobtypmx
  IF (no_iterations(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage =                                                  &
      'GROUPDEP : NO_ITERATIONS used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (no_iterations(jg) /= imdi) THEN
    icode = 1
    cmessage =                                                    &
    'GROUPDEP : Invalid value given in NO_ITERATIONS'
    GO TO 999
  ELSE
  END IF
END DO

l_no_iters = ncount >  0
IF (l_new_groups .AND. l_no_iters .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage =                                                      &
  'GROUPDEP : Wrong no of values given in NO_ITERATIONS'
  GO TO 999
END IF

!     Process INTERVAL_ITER
!     ---------------------
ncount = 0
DO jg=1,nobtypmx
  IF (interval_iter(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage =                                                  &
      'GROUPDEP : INTERVAL_ITER used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (interval_iter(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in INTERVAL_ITER'
    GO TO 999
  ELSE
  END IF
END DO

l_int_iter = ncount >  0
IF (l_new_groups .AND. l_int_iter .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in INTERVAL_ITER'
  GO TO 999
END IF

!     Process N_ANAL_LEVS
!     --------------------
ncount = 0
DO jg=1,nobtypmx
  IF (n_anal_levs(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage =                                                  &
      'GROUPDEP : N_ANAL_LEVS used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (n_anal_levs(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in N_ANAL_LEVS'
    GO TO 999
  ELSE
  END IF
END DO

l_anal_levs = ncount >  0
IF (l_new_groups .AND. l_anal_levs .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in N_ANAL_LEVS'
  GO TO 999
END IF

!     Process N_WT_LEVS
!     --------------------
ncount = 0
DO jg=1,nobtypmx
  IF (n_wt_levs(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage =                                                  &
      'GROUPDEP : N_WT_LEVS used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (n_wt_levs(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in N_WT_LEVS'
    GO TO 999
  ELSE
  END IF
END DO

l_wt_levs = ncount >  0
IF (l_new_groups .AND. l_wt_levs .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in N_WT_LEVS'
  GO TO 999
END IF

IF (model_type == mt_global) THEN
  !     Process NUDGE_NH
  !     ----------------
  ncount = 0
  DO jg=1,nobtypmx
    IF (nudge_nh(jg) >  0.0) THEN
      IF (jg >  no_groups) THEN
        icode = 1
        cmessage =                                                  &
        'GROUPDEP : NUDGE_NH used incorrectly.'
        GO TO 999
      END IF
      ncount = ncount+1
    ELSE IF (nudge_nh(jg) /= rmdi) THEN
      icode = 1
      cmessage = 'GROUPDEP : Invalid value in NUDGE_NH'
      GO TO 999
    ELSE
    END IF
  END DO

  l_nudge_nh = ncount >  0

  IF (l_new_groups .AND. l_nudge_nh .AND. ncount /= no_groups) THEN
    icode = 1
    cmessage = 'GROUPDEP : Wrong no of values in NUDGE_NH'
    GO TO 999
  END IF

  !     Process NUDGE_TR
  !     ----------------
  ncount = 0
  DO jg=1,nobtypmx
    IF (nudge_tr(jg) >  0.0) THEN
      IF (jg >  no_groups) THEN
        icode = 1
        cmessage = 'GROUPDEP : NUDGE_TR used incorrectly.'
        GO TO 999
      END IF
      ncount = ncount+1
    ELSE IF (nudge_tr(jg) /= rmdi) THEN
      icode = 1
      cmessage = 'GROUPDEP : Invalid value in NUDGE_TR'
      GO TO 999
    ELSE
    END IF
  END DO

  l_nudge_tr = ncount >  0

  IF (l_new_groups .AND. l_nudge_tr .AND. ncount /= no_groups) THEN
    icode = 1
    cmessage = 'GROUPDEP : Wrong no of values in NUDGE_TR'
    GO TO 999
  END IF

  !     Process NUDGE_SH
  !     ----------------
  ncount = 0
  DO jg=1,nobtypmx
    IF (nudge_sh(jg) >  0.0) THEN
      IF (jg >  no_groups) THEN
        icode = 1
        cmessage = 'GROUPDEP : NUDGE_SH used incorrectly.'
        GO TO 999
      END IF
      ncount = ncount+1
    ELSE IF (nudge_sh(jg) /= rmdi) THEN
      icode = 1
      cmessage = 'GROUPDEP : Invalid value in NUDGE_SH'
      GO TO 999
    ELSE
    END IF
  END DO

  l_nudge_sh = ncount >  0

  IF (l_new_groups .AND. l_nudge_sh .AND. ncount /= no_groups) THEN
    icode = 1
    cmessage = 'GROUPDEP : Wrong no of values in NUDGE_SH'
    GO TO 999
  END IF
ELSE
  !     Process NUDGE_LAM
  !     -----------------
  ncount = 0
  DO jg=1,nobtypmx
    IF (nudge_lam(jg) >  0.0) THEN
      IF (jg >  no_groups) THEN
        icode = 1
        cmessage = 'GROUPDEP : NUDGE_LAM used incorrectly.'
        GO TO 999
      END IF
      ncount = ncount+1
    ELSE IF (nudge_lam(jg) /= rmdi) THEN
      icode = 1
      cmessage = 'GROUPDEP : Invalid value in NUDGE_LAM'
      GO TO 999
    ELSE
    END IF
  END DO

  l_nudge_lam = ncount >  0

  IF (l_new_groups .AND. l_nudge_lam .AND. ncount /= no_groups) THEN
    icode = 1
    cmessage = 'GROUPDEP : Wrong no of values in NUDGE_LAM'
    GO TO 999
  END IF
END IF  ! if GLOBAL

!     Process AGRES_ROWS
!     ------------------
ncount = 0
DO jg=1,nobtypmx
  IF (agres_rows(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage = 'GROUPDEP : AGRES_ROWS used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (agres_rows(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in AGRES_ROWS'
    GO TO 999
  ELSE
  END IF
END DO

ncount = 0
DO jg=1,nobtypmx
  IF (agres_rows(jg) >  1) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,*) ' AGRES_ROWS /= 1 disallowed for MPP'
      CALL umPrint(umMessage,src='group_dep_var')
    END IF
  END IF
END DO
l_agres_rows = ncount >  0
IF (l_new_groups .AND. l_agres_rows .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in AGRES_ROWS'
  GO TO 999
END IF

!     Process AGRES_PTS
!     -----------------
ncount = 0
DO jg=1,nobtypmx
  IF (agres_pts(jg) >  0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage = 'GROUPDEP : AGRES_PTS used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (agres_pts(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in AGRES_PTS'
    GO TO 999
  ELSE
  END IF
END DO

ncount = 0
DO jg=1,nobtypmx
  IF (agres_pts(jg) >  1) THEN
    IF (mype == 0) THEN
      WRITE(umMessage,*) ' AGRES_PTS /= 1 disallowed for MPP'
      CALL umPrint(umMessage,src='group_dep_var')
    END IF
  END IF
END DO

l_agres_pts = ncount >  0

IF (l_new_groups .AND. l_agres_pts .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in AGRES_PTS'
  GO TO 999
END IF

!     Process MODE_HANAL
!     ------------------
ncount = 0
DO jg=1,nobtypmx
  IF (mode_hanal(jg) == 1 .OR. mode_hanal(jg) == 2) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage = 'GROUPDEP : MODE_HANAL used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (mode_hanal(jg) /= imdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in MODE_HANAL'
    GO TO 999
  ELSE
  END IF
END DO

l_mode_hanal = ncount >  0

IF (l_new_groups .AND. l_mode_hanal .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in MODE_HANAL'
  GO TO 999
END IF

!     Process FI_VAR_FACTOR
!     ---------------------
ncount = 0
DO jg=1,nobtypmx
  IF (fi_var_factor(jg) >= 0.0) THEN
    IF (jg >  no_groups) THEN
      icode = 1
      cmessage = 'GROUPDEP : FI_VAR_FACTOR used incorrectly.'
      GO TO 999
    END IF
    ncount = ncount+1
  ELSE IF (fi_var_factor(jg) /= rmdi) THEN
    icode = 1
    cmessage = 'GROUPDEP : Invalid value in FI_VAR_FACTOR'
    GO TO 999
  ELSE
  END IF
END DO

l_fi_var_factor = ncount >  0

IF (l_new_groups .AND. l_fi_var_factor                             &
    .AND. ncount /= no_groups) THEN
  icode = 1
  cmessage = 'GROUPDEP : Wrong no of values in FI_VAR_FACTOR'
  GO TO 999
END IF

IF (mype == 0) THEN
  IF (l_new_groups) CALL umPrint( 'AC_ORDER      used in namelist', &
      src='group_dep_var')
  IF (l_no_iters)   CALL umPrint( 'NO_ITERATIONS used in namelist', &
      src='group_dep_var')
  IF (l_int_iter)   CALL umPrint( 'INTERVAL_ITER used in namelist', &
      src='group_dep_var')
  IF (l_anal_levs)  CALL umPrint( 'N_ANAL_LEVS   used in namelist', &
      src='group_dep_var')
  IF (l_wt_levs)    CALL umPrint( 'N_WT_LEVS     used in namelist', &
      src='group_dep_var')
  IF (l_agres_rows) CALL umPrint( 'AGRES_ROWS    used in namelist', &
      src='group_dep_var')
  IF (l_agres_pts)  CALL umPrint( 'AGRES_PTS     used in namelist', &
      src='group_dep_var')
  IF (l_fi_var_factor)CALL umPrint('FI_VAR_FACTOR used in namelist', &
      src='group_dep_var')
  IF (model_type == mt_global) THEN
    IF (l_nudge_nh)   CALL umPrint( 'NUDGE_NH      used in namelist', &
        src='group_dep_var')
    IF (l_nudge_tr)   CALL umPrint( 'NUDGE_TR      used in namelist', &
    src='group_dep_var')
    IF (l_nudge_sh)   CALL umPrint( 'NUDGE_SH      used in namelist', &
        src='group_dep_var')
  ELSE
    IF (l_nudge_lam)  CALL umPrint( 'NUDGE_LAM     used in namelist', &
        src='group_dep_var')
  END IF
END IF

!     --------------------------------------------------------------
!     If AC_ORDER has been used to change the order or groups
!     but none of the group dependent arrays used, then this routine
!     attempts to derive new defaults for the new order or groups
!     from the existing defaults.
!     --------------------------------------------------------------
IF (l_new_groups) THEN  !  New order or groups.
  IF (model_type == mt_global) THEN
    l_dummy = .NOT. l_nudge_nh  .OR. .NOT. l_nudge_tr  .OR.        &
              .NOT. l_nudge_sh
  ELSE
    l_dummy = .NOT. l_nudge_lam
  END IF

  IF (.NOT. l_no_iters  .OR. .NOT. l_int_iter  .OR.                 &
      .NOT. l_anal_levs .OR. .NOT. l_wt_levs   .OR.                 &
      l_dummy .OR.                                                &
      .NOT. l_agres_rows .OR. .NOT. l_agres_pts .OR.                &
      .NOT. l_mode_hanal .OR. .NOT. l_fi_var_factor) THEN

    !         Set up new defaults in NEW_ arrays from existing
    !         defaults. New defaults are set up for all observation
    !         types in AC_ORDER from the corresponding default array.
    !         eg. NEW_NO_ITERS is set up from DEF_NO_ITERATIONS.
    !         New defaults are set up for all group dependent arrays
    !         here. Those not required are ignored later in this routine.

    IF (model_type == mt_global) THEN
      DO jobt=1,n_obtyp
        this_type = MOD(ac_order(jobt),1000)
        DO jobt2=1,nobtypmx
          this_type_def = MOD(def_ac_order(jobt2),1000)
          IF (this_type == this_type_def) THEN
            this_group = def_ac_order(jobt2)/1000
            new_no_iters(jobt)   = def_no_iterations(this_group)
            new_int_iter(jobt)   = def_interval_iter(this_group)
            new_anal_levs(jobt)  = def_no_anal_levs(this_group)
            new_wt_levs(jobt)    = def_no_wt_levs(this_group)
            new_nudge_nh(jobt)   = def_nudge_nh(this_group)
            new_nudge_tr(jobt)   = def_nudge_tr(this_group)
            new_nudge_sh(jobt)   = def_nudge_sh(this_group)
            new_agres_rows(jobt) = def_agres_rows(this_group)
            new_agres_pts(jobt)  = def_agres_pts(this_group)
            new_mode_hanal(jobt) = def_mode_hanal(this_group)
            new_fi_var_factor(jobt) = def_fi_var_factor(this_group)
          END IF
        END DO
      END DO
    ELSE
      DO jobt=1,n_obtyp
        this_type = MOD(ac_order(jobt),1000)
        DO jobt2=1,nobtypmx
          this_type_def = MOD(def_ac_order(jobt2),1000)
          IF (this_type == this_type_def) THEN
            this_group = def_ac_order(jobt2)/1000
            new_no_iters(jobt)   = def_no_iterations(this_group)
            new_int_iter(jobt)   = def_interval_iter(this_group)
            new_anal_levs(jobt)  = def_no_anal_levs(this_group)
            new_wt_levs(jobt)    = def_no_wt_levs(this_group)
            new_nudge_lam(jobt)  = def_nudge_lam(this_group)
            new_agres_rows(jobt) = def_agres_rows(this_group)
            new_agres_pts(jobt)  = def_agres_pts(this_group)
            new_mode_hanal(jobt) = def_mode_hanal(this_group)
            new_fi_var_factor(jobt) = def_fi_var_factor(this_group)
          END IF
        END DO
      END DO
    END IF  ! if GLOBAL

    IF (nprog == 1001 .AND. mype == 0) THEN
      WRITE(umMessage,*) 'NEW_NO_ITERS'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_no_iters(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_INT_ITER'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_int_iter(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_ANAL_LEVS'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_anal_levs(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_WT_LEVS'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_wt_levs(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      IF (model_type == mt_global) THEN
        WRITE(umMessage,*) 'NEW_NUDGE_NH'
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,'(10E10.3)') (new_nudge_nh(jobt),jobt=1,n_obtyp)
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,*) 'NEW_NUDGE_TR'
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,'(10E10.3)') (new_nudge_tr(jobt),jobt=1,n_obtyp)
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,*) 'NEW_NUDGE_SH'
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,'(10E10.3)') (new_nudge_sh(jobt),jobt=1,n_obtyp)
        CALL umPrint(umMessage,src='group_dep_var')
      ELSE
        WRITE(umMessage,*) 'NEW_NUDGE_LAM'
        CALL umPrint(umMessage,src='group_dep_var')
        WRITE(umMessage,'(10E10.3)') (new_nudge_lam(jobt),jobt=1,n_obtyp)
        CALL umPrint(umMessage,src='group_dep_var')
      END IF
      WRITE(umMessage,*) 'NEW_AGRES_ROWS'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_agres_rows(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_AGRES_PTS'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_agres_pts(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_MODE_HANAL'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10I7)') (new_mode_hanal(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,*) 'NEW_FI_VAR_FACTOR'
      CALL umPrint(umMessage,src='group_dep_var')
      WRITE(umMessage,'(10E10.3)')(new_fi_var_factor(jobt),jobt=1,n_obtyp)
      CALL umPrint(umMessage,src='group_dep_var')
    END IF

    !         Go through groups and check that the new default for each
    !         obs type in the group matches. If there is a mismatch, then
    !         the routine will abort and prompt the user to use the
    !         appropriate namelist array to specify new values for the
    !         new groups. If the new defaults for all the obs types in the
    !         group match, then this becomes the default for the new group.

    DO jg=1,no_groups

      first_type = 0
      last_type  = 0
      DO jobt=1,n_obtyp
        this_group = ac_order(jobt)/1000
        IF (this_group == jg) THEN
          IF (first_type == 0) THEN
            first_type = jobt
          END IF
          last_type = jobt
        END IF
      END DO

      IF (first_type == 0 .OR. last_type == 0) THEN
        icode = 1
        cmessage ='GROUPDEP : FIRST_TYPE=0 or LAST_TYPE=0 ?'
        GO TO 999
      END IF

      IF (.NOT. l_no_iters) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_no_iters(first_type) /= new_no_iters(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : NO_ITERATIONS must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_no_iterations(jg) = new_no_iters(first_type)
      END IF

      IF (.NOT. l_int_iter) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_int_iter(first_type) /= new_int_iter(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : INTERVAL_ITER must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_interval_iter(jg) = new_int_iter(first_type)
      END IF

      IF (.NOT. l_anal_levs) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_anal_levs(first_type) /= new_anal_levs(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : N_ANAL_LEVS must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_no_anal_levs(jg) = new_anal_levs(first_type)
      END IF

      IF (.NOT. l_wt_levs) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_wt_levs(first_type) /= new_wt_levs(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : N_WT_LEVS must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_no_wt_levs(jg) = new_wt_levs(first_type)
      END IF

      IF (.NOT. l_agres_rows) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_agres_rows(first_type)  /=                      &
                new_agres_rows(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : AGRES_ROWS must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_agres_rows(jg) = new_agres_rows(first_type)
      END IF

      IF (.NOT. l_agres_pts) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_agres_pts(first_type)  /=                       &
                new_agres_pts(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : AGRES_PTS must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_agres_pts(jg) = new_agres_pts(first_type)
      END IF

      IF (.NOT. l_mode_hanal) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_mode_hanal(first_type)  /=                      &
                new_mode_hanal(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : MODE_HANAL must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_mode_hanal(jg) = new_mode_hanal(first_type)
      END IF

      IF (.NOT. l_fi_var_factor) THEN
        IF (first_type <  last_type) THEN
          DO jobt=first_type+1,last_type
            IF (new_fi_var_factor(first_type)  /=                   &
                new_fi_var_factor(jobt)) THEN
              icode = 1
              cmessage =                                            &
              'GROUPDEP : FI_VAR_FACTOR must be used in namelist'
              GO TO 999
            END IF
          END DO
        END IF
        def_fi_var_factor(jg) = new_fi_var_factor(first_type)
      END IF

      IF (model_type == mt_global) THEN
        IF (.NOT. l_nudge_nh) THEN
          IF (first_type <  last_type) THEN
            DO jobt=first_type+1,last_type
              IF (new_nudge_nh(first_type) /= new_nudge_nh(jobt)) THEN
                icode = 1
                cmessage =                                            &
                'GROUPDEP : NUDGE_NH must be used in namelist'
                GO TO 999
              END IF
            END DO
          END IF
          def_nudge_nh(jg) = new_nudge_nh(first_type)
        END IF

        IF (.NOT. l_nudge_tr) THEN
          IF (first_type <  last_type) THEN
            DO jobt=first_type+1,last_type
              IF (new_nudge_tr(first_type) /= new_nudge_tr(jobt)) THEN
                icode = 1
                cmessage =                                            &
                'GROUPDEP : NUDGE_TR must be used in namelist'
                GO TO 999
              END IF
            END DO
          END IF
          def_nudge_tr(jg) = new_nudge_tr(first_type)
        END IF

        IF (.NOT. l_nudge_sh) THEN
          IF (first_type <  last_type) THEN
            DO jobt=first_type+1,last_type
              IF (new_nudge_sh(first_type) /= new_nudge_sh(jobt)) THEN
                icode = 1
                cmessage =                                            &
                'GROUPDEP : NUDGE_SH must be used in namelist'
                GO TO 999
              END IF
            END DO
          END IF
          def_nudge_sh(jg) = new_nudge_sh(first_type)
        END IF

      ELSE
        IF (.NOT. l_nudge_lam) THEN
          IF (first_type <  last_type) THEN
            DO jobt=first_type+1,last_type
              IF (new_nudge_lam(first_type)  /=                       &
                  new_nudge_lam(jobt)) THEN
                icode = 1
                cmessage =                                            &
                'GROUPDEP : NUDGE_LAM must be used in namelist'
                GO TO 999
              END IF
            END DO
          END IF
          def_nudge_lam(jg) = new_nudge_lam(first_type)
        END IF

      END IF  ! if GLOBAL
    END DO

  END IF
END IF

!     Overwrite existing defaults with new defaults.

IF (l_new_groups) THEN   !  New values for DEF_AC_ORDER
  DO jobt=1,nobtypmx
    def_ac_order(jobt) = ac_order(jobt)
  END DO
END IF

IF (l_no_iters) THEN   !  New values for DEF_NO_ITERATIONS
  DO jg=1,no_groups
    IF (no_iterations(jg) >  0) THEN
      def_no_iterations(jg) = no_iterations(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_no_iterations(jg) = imdi
    END DO
  END IF
END IF

IF (l_int_iter) THEN   !   New values for DEF_INTERVAL_ITER
  DO jg=1,no_groups
    IF (interval_iter(jg) >  0) THEN
      def_interval_iter(jg) = interval_iter(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_interval_iter(jg) = imdi
    END DO
  END IF
END IF

IF (l_anal_levs) THEN   !   New values for DEF_NO_ANAL_LEVS
  DO jg=1,no_groups
    IF (n_anal_levs(jg) >  0) THEN
      def_no_anal_levs(jg) = n_anal_levs(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_no_anal_levs(jg) = imdi
    END DO
  END IF
END IF

IF (l_wt_levs) THEN   !   New values for DEF_NO_WT_LEVS
  DO jg=1,no_groups
    IF (n_wt_levs(jg) >  0) THEN
      def_no_wt_levs(jg) = n_wt_levs(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_no_wt_levs(jg) = imdi
    END DO
  END IF
END IF

IF (l_agres_rows) THEN  !   New values for DEF_AGRES_ROWS
  DO jg=1,no_groups
    IF (agres_rows(jg) >  0) THEN
      def_agres_rows(jg) = agres_rows(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_agres_rows(jg) = imdi
    END DO
  END IF
END IF

IF (l_agres_pts) THEN  !   New values for DEF_AGRES_PTS
  DO jg=1,no_groups
    IF (agres_pts(jg) >  0) THEN
      def_agres_pts(jg) = agres_pts(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_agres_pts(jg) = imdi
    END DO
  END IF
END IF

IF (l_mode_hanal) THEN  !   New values for DEF_MODE_HANAL
  DO jg=1,no_groups
    IF (mode_hanal(jg) >  0) THEN
      def_mode_hanal(jg) = mode_hanal(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_mode_hanal(jg) = imdi
    END DO
  END IF
END IF

IF (l_fi_var_factor) THEN ! New values for DEF_FI_VAR_FACTOR
  DO jg=1,no_groups
    IF (fi_var_factor(jg) >  0.0) THEN
      def_fi_var_factor(jg) = fi_var_factor(jg)
    END IF
  END DO
  IF (no_groups <  nobtypmx) THEN
    DO jg=no_groups+1,nobtypmx
      def_fi_var_factor(jg) = rmdi
    END DO
  END IF
END IF

IF (model_type == mt_global) THEN
  IF (l_nudge_nh) THEN   !   New values for DEF_NUDGE_NH
    DO jg=1,no_groups
      IF (nudge_nh(jg) >  0.0) THEN
        def_nudge_nh(jg) = nudge_nh(jg)
      END IF
    END DO
    IF (no_groups <  nobtypmx) THEN
      DO jg=no_groups+1,nobtypmx
        def_nudge_nh(jg) = rmdi
      END DO
    END IF
  END IF

  IF (l_nudge_tr) THEN   !   New values for DEF_NUDGE_TR
    DO jg=1,no_groups
      IF (nudge_tr(jg) >  0.0) THEN
        def_nudge_tr(jg) = nudge_tr(jg)
      END IF
    END DO
    IF (no_groups <  nobtypmx) THEN
      DO jg=no_groups+1,nobtypmx
        def_nudge_tr(jg) = rmdi
      END DO
    END IF
  END IF

  IF (l_nudge_sh) THEN   !   New values for DEF_NUDGE_SH
    DO jg=1,no_groups
      IF (nudge_sh(jg) >  0.0) THEN
        def_nudge_sh(jg) = nudge_sh(jg)
      END IF
    END DO
    IF (no_groups <  nobtypmx) THEN
      DO jg=no_groups+1,nobtypmx
        def_nudge_sh(jg) = rmdi
      END DO
    END IF
  END IF
ELSE

  IF (l_nudge_lam) THEN   !   New values for DEF_NUDGE_LAM
    DO jg=1,no_groups
      IF (nudge_lam(jg) >  0.0) THEN
        def_nudge_lam(jg) = nudge_lam(jg)
      END IF
    END DO
    IF (no_groups <  nobtypmx) THEN
      DO jg=no_groups+1,nobtypmx
        def_nudge_lam(jg) = rmdi
      END DO
    END IF
  END IF
END IF  ! if GLOBAL

!     Check that no of wt levs is < or = of no of anal levs.
IF (l_anal_levs .OR. l_wt_levs) THEN
  DO jg = 1,no_groups
    IF (def_no_anal_levs(jg) <  def_no_wt_levs(jg)) THEN
      icode = 1
      cmessage = 'GROUPDEP: No of Anal Levs < No of Wt Levs ?'
      WRITE(umMessage,*) 'Group No ',jg,' No of anal/wt levs = ',        &
          def_no_anal_levs(jg),def_no_wt_levs(jg)
      CALL umPrint(umMessage,src='group_dep_var',pe=0)
      GO TO 999
    END IF
  END DO
END IF

IF (nprog == 1001 .AND. mype == 0) THEN

  WRITE(umMessage,*) 'DEF_AC_ORDER'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_ac_order(jobt),jobt=1,nobtypmx)
  CALL umPrint(umMessage,src='group_dep_var')

  WRITE(umMessage,*) 'GROUPNO'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (jg,jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'NO_ITERATIONS'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (no_iterations(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'INTERVAL_ITER'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (interval_iter(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_NO_ITERATIONS'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_no_iterations(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_INTERVAL_ITER'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_interval_iter(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_AGRES_ROWS'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_agres_rows(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_AGRES_PTS'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_agres_pts(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_MODE_HANAL'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10I7)') (def_mode_hanal(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,*) 'DEF_FI_VAR_FACTOR'
  CALL umPrint(umMessage,src='group_dep_var')
  WRITE(umMessage,'(10E10.3)') (def_fi_var_factor(jg),jg=1,no_groups)
  CALL umPrint(umMessage,src='group_dep_var')
  IF (model_type == mt_global) THEN
    WRITE(umMessage,*) 'DEF_NUDGE_NH'
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,'(10E10.3)') (def_nudge_nh(jg),jg=1,no_groups)
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,*) 'DEF_NUDGE_TR'
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,'(10E10.3)') (def_nudge_tr(jg),jg=1,no_groups)
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,*) 'DEF_NUDGE_SH'
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,'(10E10.3)') (def_nudge_sh(jg),jg=1,no_groups)
    CALL umPrint(umMessage,src='group_dep_var')
  ELSE
    WRITE(umMessage,*) 'DEF_NUDGE_LAM'
    CALL umPrint(umMessage,src='group_dep_var')
    WRITE(umMessage,'(10E10.3)') (def_nudge_lam(jg),jg=1,no_groups)
    CALL umPrint(umMessage,src='group_dep_var')
  END IF
END IF

999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE group_dep_var
END MODULE group_dep_var_mod
