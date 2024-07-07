! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initial setting of interpolation logical switches.

MODULE Rcf_Set_Interp_Logicals_Mod

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

!  Subroutine Rcf_Set_Interp_Logicals
!
! Description:
!    Sets the h_int_active, v_int_active and v_int_active_soil
!    logical switches that control horizontal and vertical
!    interpolation. This is just an initial setting and may be
!    changed at a later time.
!
! Method:
!    h_int_active - Checks for domain and row/row-length changes
!    v_int_active - Checks for # of model/wet level changes and
!                   changes to eta values or height field constants
!    v_int_active_soil - Checks for changes in number and values of
!                        soil levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

PRIVATE :: rcf_check_horiz_change, lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_INTERP_LOGICALS_MOD'

CONTAINS

SUBROUTINE Rcf_Set_Interp_Logicals( Input_Grid, Output_Grid, Hdr_In)

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_UMhead_Mod, ONLY: &
    UM_Header_Type

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active,                  &
    h_int_active_u,                &
    h_int_active_v

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active,            &
    v_int_active_soil

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal,            &
    PrStatus_Oper

USE Rcf_HeadAddress_Mod, ONLY:                                &
    RC_LongSpacing,        RC_LatSpacing,         RC_FirstLat, &
    RC_FirstLong,          RC_PoleLat,            RC_PoleLong

USE Rcf_Generate_Heights_Mod, ONLY: &
    height_gen_original

USE UM_ParCore, ONLY: &
    mype

IMPLICIT NONE

! Arguments
TYPE( Grid_Type ), INTENT(IN)       :: Input_Grid
TYPE( Grid_Type ), INTENT(IN)       :: Output_Grid
TYPE( UM_Header_Type ), INTENT(IN)  :: Hdr_In

! Local variables
INTEGER                             :: i        ! Looper
CHARACTER (LEN=40)                  :: reason   ! why is interp.
                                                ! switched on.
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SET_INTERP_LOGICALS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!--------------------------------------------------------------------
! Horizontal Interpolation
!--------------------------------------------------------------------
! Default for all gridtypes.
h_int_active = .FALSE.
h_int_active_u = .FALSE.
h_int_active_v = .FALSE.

! First check P grid.
CALL rcf_check_horiz_change(input_grid, output_grid, "p grid", h_int_active)

! If P grid has changed we assume all gridtypes must change.
IF (h_int_active) THEN
  CALL umPrint( 'Horizontal interpolation is switched ON for u and v grid.', &
      src='rcf_set_interp_logicals_mod')
  IF (PrintStatus >= PrStatus_Oper) THEN
    CALL umPrint( 'Because of p grid being interpolated.', &
        src='rcf_set_interp_logicals_mod')
  END IF
  h_int_active_u = .TRUE.
  h_int_active_v = .TRUE.
ELSE
  ! If P grid has not changed we can check the other grids to make sure
  ! staggering has not changed.  ND LAM -> EG LAM makes this possible.
  CALL rcf_check_horiz_change(input_grid, output_grid, "u grid", h_int_active_u)
  CALL rcf_check_horiz_change(input_grid, output_grid, "v grid", h_int_active_v)
END IF

!---------------------------------------------------------------------
! Vertical interpolation
!---------------------------------------------------------------------
v_int_active = .FALSE.

! Turn on if number of levels varies
IF (Input_Grid % model_levels /= Output_Grid % model_levels) THEN

  v_int_active = .TRUE.
  reason = 'changed number of model levels'

ELSE   ! Turn on if eta levels vary
  DO i = 1, Hdr_In % Len1LevDepC
    IF ( ABS( Hdr_In % LevDepC( i, 1 ) -                            &
         Output_Grid % eta_theta_levels( i-1)) > EPSILON( 1.0 )) THEN

      v_int_active = .TRUE.
      reason = 'change in eta_theta levels'
    END IF
  END DO

  DO i = 1, Hdr_In % Len1LevDepC - 1
    IF ( ABS( Hdr_In % LevDepC( i, 2 ) -                            &
         Output_Grid % eta_rho_levels( i ) ) > EPSILON( 1.0 ) ) THEN

      v_int_active = .TRUE.
      reason = 'change in eta_rho levels'
    END IF
  END DO
END IF

! Turn on if height specifying constants vary
IF ( (Input_Grid  % first_constant_r_rho_level /=     &
      Output_Grid % first_constant_r_rho_level ) .OR. &
     (Input_Grid % z_top_of_model - Output_Grid % z_top_of_model ) > &
                                    EPSILON( 1.0 ) ) THEN
  v_int_active = .TRUE.
  reason = 'change in height defining constants'
END IF

! Turn on if number of boundary levels change
IF ( Input_Grid % Height_Gen_Method == height_gen_original .AND. &
     Input_Grid % BL_Levels /= Output_Grid % BL_Levels ) THEN
  v_int_active = .TRUE.
  reason = 'change in number of boundary levels'
END IF

! Turn on if change in way heights are calculated
IF ( Input_Grid  % Height_Gen_Method /=     &
     Output_Grid % Height_Gen_Method ) THEN
  v_int_active = .TRUE.
  reason = 'change in method of height generation'
END IF

IF (v_int_active) THEN
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) 'Vertical interpolation is switched ON'
    CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,*) 'Because of ', reason
      CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
    END IF
  END IF
ELSE
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) 'Vertical interpolation is switched OFF'
    CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
  END IF
END IF

!---------------------------------------------------------------------
! Vertical interpolation for soil levels
!---------------------------------------------------------------------
v_int_active_soil = .FALSE.
IF ( Input_Grid % sm_levels /= Output_Grid % sm_levels) THEN
  v_int_active_soil = .TRUE.
  reason = 'changed number of soil depths'
ELSE
  DO i = 1, Input_Grid % sm_levels
    IF ( ABS( Input_Grid % soil_depths(i) -                          &
              Output_Grid % soil_depths(i) ) > EPSILON(1.0) ) THEN

      v_int_active_soil = .TRUE.
      reason = 'changed soil depth'
    END IF
  END DO
END IF

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  IF (v_int_active_soil) THEN
    WRITE(umMessage,*) 'Vertical interpolation between soil depths is '//&
                 'switched ON'
    CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,*) 'Because of ', reason
      CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
    END IF
  ELSE
    WRITE(umMessage,*) 'Vertical interpolation between soil depths is '//&
                  'switched OFF'
    CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
  END IF
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='rcf_set_interp_logicals_mod')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_set_interp_logicals

!  Subroutine Rcf_check_horiz_change
!
! Description:
!   Given a input and output grid and a gridtype to check it will set a flag
!   corresponding to whether the grid has changed.
!
SUBROUTINE rcf_check_horiz_change(input_grid, output_grid, gridtype, flag)

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type
USE UM_ParCore, ONLY: &
    mype
USE umPrintMgr

IMPLICIT NONE
! Arguments
TYPE( Grid_Type ), INTENT(IN)       :: Input_Grid
TYPE( Grid_Type ), INTENT(IN)       :: Output_Grid
CHARACTER (LEN=*)                   :: gridtype
LOGICAL          , INTENT(OUT)      :: flag

REAL, POINTER :: lambda_in (:)
REAL, POINTER :: phi_in    (:)
REAL, POINTER :: lambda_out(:)
REAL, POINTER :: phi_out   (:)
INTEGER       :: row_length_in
INTEGER       :: rows_in
INTEGER       :: row_length_out
INTEGER       :: rows_out
INTEGER       :: i
CHARACTER (LEN=60)                  :: reason   ! why is interp.
                                                ! switched on.

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CHECK_HORIZ_CHANGE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

flag = .FALSE.

IF (gridtype == "p grid") THEN
  row_length_in  =  input_grid  % glob_p_row_length
  rows_in        =  input_grid  % glob_p_rows
  phi_in         => input_grid  % phi_p
  lambda_in      => input_grid  % lambda_p
  row_length_out =  output_grid % glob_p_row_length
  rows_out       =  output_grid % glob_p_rows
  phi_out        => output_grid % phi_p
  lambda_out     => output_grid % lambda_p
ELSE IF (gridtype == "u grid") THEN
  row_length_in  =  input_grid  % glob_u_row_length
  rows_in        =  input_grid  % glob_u_rows
  phi_in         => input_grid  % phi_p
  lambda_in      => input_grid  % lambda_u
  row_length_out =  output_grid % glob_u_row_length
  rows_out       =  output_grid % glob_u_rows
  phi_out        => output_grid % phi_p
  lambda_out     => output_grid % lambda_u
ELSE IF (gridtype == "v grid") THEN
  row_length_in  =  input_grid  % glob_v_row_length
  rows_in        =  input_grid  % glob_v_rows
  phi_in         => input_grid  % phi_v
  lambda_in      => input_grid  % lambda_p
  row_length_out =  output_grid % glob_v_row_length
  rows_out       =  output_grid % glob_v_rows
  phi_out        => output_grid % phi_v
  lambda_out     => output_grid % lambda_p
END IF


IF ( .NOT. Output_Grid % Global ) THEN
  IF ( ABS( phi_in(1) - phi_out(1) ) > EPSILON( 1.0 ) ) THEN
    flag = .TRUE.
    reason = 'differing first latitude for ' // gridtype
  END IF
END IF

IF ( ABS( lambda_in(1) - lambda_out(1) ) > EPSILON( 1.0 ) .OR. &
     ABS( input_grid % lambda_npole - output_grid % lambda_npole ) >       &
                                            EPSILON( 1.0 ) .OR. &
     ABS( input_grid % phi_npole    - output_grid % phi_npole )    >       &
                                            EPSILON( 1.0 ) ) THEN
  flag = .TRUE.
  reason = 'changed domain for ' // gridtype
END IF

! Check horizontal grid resolution.
IF ( rows_in       /= rows_out       .OR. &
     row_length_in /= row_length_out) THEN
  flag = .TRUE.
  reason = 'changed horizontal resolution for ' // gridtype
ELSE
  ! If the input and output grids have the same number of rows and points on rows
  ! we can compare the actual grid coordinates due to variable resolution is
  ! possible.
  DO i = 1, rows_in
    IF (ABS(phi_in(i) - phi_out(i)) > EPSILON(1.0)) THEN
      flag = .TRUE.
      reason = 'changed phi grid information for ' // gridtype
      EXIT
    END IF
  END DO

  DO i = 1, row_length_in
    IF (ABS(lambda_in(i) - lambda_out(i)) > EPSILON(1.0)) THEN
      flag = .TRUE.
      reason = 'changed lambda grid information for ' // gridtype
      EXIT
    END IF
  END DO
END IF

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  CALL umPrint('',src='rcf_set_interp_logicals_mod')
  IF (flag) THEN
    CALL umPrint('Horizontal interpolation is switched ON for ' // gridtype,   &
        src='rcf_set_interp_logicals_mod')
    IF (PrintStatus >= PrStatus_Oper) THEN
      CALL umPrint('Because of '// reason,src='rcf_set_interp_logicals_mod')
    END IF
  ELSE
    CALL umPrint('Horizontal interpolation is switched OFF for ' // gridtype)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_check_horiz_change
END MODULE Rcf_Set_Interp_Logicals_Mod
