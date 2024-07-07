! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE DEF_GROUP ----------------------------------------------
!
!    Purpose : Initialise Group Dependent arrays
!
!    For Global runs : Enable defs GLOBAL
!
!    Project Task : P3
!
!    Documentation:
!
!
!    ARGUMENTS :--------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE def_group_mod

USE missing_data_mod, ONLY: rmdi, imdi 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEF_GROUP_MOD'

CONTAINS

SUBROUTINE def_group (bl_levels,tr_levels,      &
                      icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global
USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

INTEGER ::                                                        &
   bl_levels                                                      &
                 ! IN - No of levels in boundary layer.
 , tr_levels                                                      &
                 ! IN - No of Tracer levels.
 , icode         ! OUT - Return Code
CHARACTER(LEN=errormessagelength) :: cmessage  ! OUT - Error Message




!    Workspace Usage ---------------------------------------------------
!     Local array.
INTEGER :: ac_groups (nobtypmx)  ! Stores default groups of AC Types.
!     Local variables.
INTEGER :: jobt  !  Loop counter.
INTEGER :: i     !  group identifier

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEF_GROUP'


!     Order of processing and grouping of Observation Types
!     =====================================================

!     Numbers = (Group No*1000) + Obs Type

!  Default groupings are:
!  1:pstar
!  2:upper level temperatures
!  3:surface temperatures
!  4:upper level winds
!  5:surface winds
!  6:upper level RH
!  7:surface RH
!  8:MOPS RH
!  9:Cloud histograms
!  10:MOPS precip rate/phase
!  11:Tracers (NB. when the assimilation is actually run, each tracer
!     must be in a separate group; they are only put together here for
!     convenience, to define common assimilation parameters.
!     Normally only a few would be used at once.)

DATA ac_groups /                                                  &
  1101,                                                           &
  2201,2203,2205,2206,2207,2208,2209,2211,                        &
  3202,3204,                                                      &
  4301,4303,4311,                                                 &
  5302,5304,5305,5306,                                            &
  6401,6403,6405,                                                 &
  7402,7404,                                                      &
  8406,                                                           &
  9407,                                                           &
  10506,                                                          &
  11601, 11602, 11603, 11604, 11605,                              &
  11606, 11607, 11608, 11609, 11610,                              &
  11611, 11612, 11613, 11614, 11615,                              &
  11616, 11617, 11618, 11619, 11620,                              &
  11621, 11622, 11623, 11624, 11625,                              &
  11626, 11627, 11628, 11629,                                     &
  12901, 126*0/
!     NB. the number of tracers (group 11) should correspond to
!     A_MAX_TRVARS (in free_tracers_inputs_mod), currently 29

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Copy above list into namelist ACP array DEF_AC_ORDER

DO jobt=1,nobtypmx
  def_ac_order(jobt) = ac_groups(jobt)
END DO

!     Group Dependent Variables
!     =========================

!     DEF_NO_ANAL_LEVS  No of analysis levels
!     DEF_NO_WT_LEVS    No of weight levels
!     DEF_NO_ITERATIONS No of iterations
!     DEF_INTERVAL_ITER Interval in timesteps between iterations
!     DEF_AGRES_ROWS Ratio of No of rows in Model Grid to Analysis Grid
!     DEF_AGRES_PTS  Ratio of No of pts in Model Grid to Analysis Grid
!     DEF_MODE_HANAL Mode of Horizontal Analysis
!     DEF_FI_VAR_FACTOR group dep scaling factor in FI
!     DEF_NUDGE_NH   Nudging Coefficients for NH
!     DEF_NUDGE_TR   Nudging Coefficients for TR
!     DEF_NUDGE_SH   Nudging Coefficients for SH
!     DEF_NUDGE_LAM  Nudging Coefficients for LAM

IF (model_type == mt_global) THEN

  i = 1                             ! pstar
  def_no_anal_levs(i)  = 1
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.8e-4
  END IF

  i = 2                             ! upper level temps
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.8e-4
  END IF

  i = 3                             ! surf temps
  def_no_anal_levs(i)  = MAX(bl_levels-2,3)
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.8e-4
  END IF

  i = 4                             !upper level winds
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35 * 3.0
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 6.0e-4
    def_nudge_tr(i) = 6.0e-4
    def_nudge_sh(i) = 6.0e-4
  ELSE
    def_nudge_nh(i) = 6.6e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.3e-4
  END IF

  i = 5                             ! surf winds (inc scat)
  def_no_anal_levs(i)  = bl_levels
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35 * 3.0
  IF (lac_uars) THEN
    def_nudge_nh(i) = 6.0e-4
    def_nudge_tr(i) = 6.0e-4
    def_nudge_sh(i) = 6.0e-4
  ELSE
    def_nudge_nh(i) = 6.3e-4
    def_nudge_tr(i) = 3.8e-4

    def_nudge_sh(i) = 4.0e-4
  END IF

  i = 6                             ! upper level RH
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.5e-4
  END IF

  i = 7                             ! surface RH
  def_no_anal_levs(i)  = MAX(bl_levels-2,3)
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.5e-4
  END IF

  i = 8                             ! MOPS RH
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.5e-4
  END IF

  i = 9                             ! Cloud histograms
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 2.0
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.5e-4
  END IF

  i = 10                  ! MOPS precip rate/phase
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = 1
  def_no_wt_levs(i)    = 1
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 1.0e6
    def_nudge_tr(i) = 1.0e6
    def_nudge_sh(i) = 1.0e6
  END IF

  i = 11                            ! tracers
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = tr_levels
  def_no_wt_levs(i)    = tr_levels
  IF (lac_uars) THEN
    def_nudge_nh(i) = 4.0e-4
    def_nudge_tr(i) = 4.0e-4
    def_nudge_sh(i) = 4.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.8e-4
  END IF

  i = 12                            ! surface LOG Visibility
  def_no_anal_levs(i)  = bl_levels
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  IF (lac_uars) THEN
    def_nudge_nh(i) = 3.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.0e-4
  ELSE
    def_nudge_nh(i) = 5.0e-4
    def_nudge_tr(i) = 3.0e-4
    def_nudge_sh(i) = 3.5e-4
  END IF

ELSE

  i = 1                             !pstar
  def_no_anal_levs(i)  = 1
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4

  i = 2                             !upper level temps
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4

  i = 3                             !surf temps
  def_no_anal_levs(i)  = MAX(bl_levels-2,3)
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4
  IF (lac_mes) THEN
    def_no_anal_levs(i)  = 6
  END IF

  i = 4                             !upper level winds
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35 * 3.0
  def_nudge_lam(i)   = 6.6e-4

  i = 5                            !surf winds (inc scat)
  def_no_anal_levs(i)  = bl_levels
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35 * 3.0
  def_nudge_lam(i)   = 6.3e-4
  IF (lac_mes) THEN
    def_no_anal_levs(i)  = 6
  END IF

  i = 6                            !upper level RH
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4

  i = 7                            !surface RH
  def_no_anal_levs(i)  = MAX(bl_levels-2,3)
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4
  IF (lac_mes) THEN
    def_no_anal_levs(i)  = 6
  END IF

  i = 8                            !MOPS RH
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 1.0e-3

  i = 9                             ! Cloud histograms
  def_no_anal_levs(i)  = model_levels
  def_no_wt_levs(i)    = model_levels
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 2.0
  def_nudge_lam(i)   = 5.0e-4

  i = 10                  ! MOPS precip rate/phase
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 2
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = 1
  def_no_wt_levs(i)    = 1
  def_nudge_lam(i)     = 1.0e6

  i = 11                            ! tracers
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_no_anal_levs(i)  = tr_levels
  def_no_wt_levs(i)    = tr_levels
  def_nudge_lam(i)  = 5.0e-4

  i = 12                           !surface LOG Visibility
  def_no_anal_levs(i)  = bl_levels
  def_no_wt_levs(i)    = 1
  def_no_iterations(i) = 1
  def_interval_iter(i) = 1
  def_agres_rows(i)    = 1
  def_agres_pts(i)     = 1
  def_mode_hanal(i)    = 1
  def_fi_var_factor(i) = 1.35
  def_nudge_lam(i)   = 5.0e-4
  IF (lac_mes) THEN
    def_no_anal_levs(i)  = 6
  END IF

END IF  ! if GLOBAL

DO jobt=i+1,nobtypmx
  def_no_anal_levs(jobt)  = imdi
  def_no_wt_levs(jobt)    = imdi
  def_no_iterations(jobt) = imdi
  def_interval_iter(jobt) = imdi
  def_agres_rows(jobt)    = imdi
  def_agres_pts(jobt)     = imdi
  def_mode_hanal(jobt)    = imdi
  def_fi_var_factor(jobt) = rmdi
  IF (model_type == mt_global) THEN
    def_nudge_nh(jobt) = rmdi
    def_nudge_tr(jobt) = rmdi
    def_nudge_sh(jobt) = rmdi

  ELSE
    def_nudge_lam(jobt) = rmdi
  END IF
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE def_group
END MODULE def_group_mod
