! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Main control routine for the nudging

!  Loads variable, interpolates it onto model levels
!  Performs the nudging and calculates diagnostics

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_nudging_control_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,           &
               PRIVATE :: ModuleName = 'NUDGING_NUDGING_CONTROL_MOD'
CONTAINS

SUBROUTINE nudging_nudging_control(  &
  varname,                          &     ! Variable name
  grid_type,                        &     ! Grid type
  global_row_length,                &     ! (Global) row length
  global_rows,                      &     ! (Global) rows
  file_levels,                      &     ! File levels
  proc_row_length_min,              &     ! Min. column
  proc_row_length_max,              &     ! Max. column
  proc_rows_min,                    &     ! Min. row
  proc_rows_max,                    &     ! Max. row
  sin_theta_latitude_min,           &     ! sine of min latitude
  sin_theta_latitude_max,           &     ! sine of max latitude
  base_lambda,                      &     ! minimum longitude in model
  delta_lambda,                     &     ! longitude increment in model
  base_phi,                         &     ! minimum latitude in model
  delta_phi,                        &     ! latitude increment in model
  base_lambda_ana,                  &     ! minimum longitude in data
  delta_lambda_ana,                 &     ! longitude increment in data
  base_phi_ana,                     &     ! minimum latitude in data
  delta_phi_ana,                    &     ! latitude increment in data
  variable_file_data,               &     ! data on first timestep
  variable_file_data2,              &     ! data on second timestep
  logp_surf_data,                   &     ! pressure on first timestep
  logp_surf_data2,                  &     ! pressure on second timestep
  frac,                             &     ! Fraction between 1 & 2
  model_levels,                     &     ! Model levels
  model_pressure_model_levels,      &     ! Model pressure
  model_trop_pressure,              &     ! Tropopause Pressure
  model_variable_model_levels,      &     ! Variable nudged
  l_data_diag, diag_var_data,       &     ! Diag. 1 (analyses value)
  l_model_diag, diag_var_model,     &     ! Diag. 2 (model value)
  l_ntend_diag, diag_tend_nudging,  &     ! Diag. 3 (nudging tendency)
  l_mtend_diag, diag_tend_model,    &     ! Diag. 4 (model tendency)
  l_relax_diag, diag_relax)               ! Diag. 5 (relaxation paremeter)

USE timestep_mod, ONLY: timestep
USE nudging_control                       ! use standard nudging switches
USE umPrintMgr

USE nudging_calc_diags_mod, ONLY: nudging_calc_diags
USE nudging_calc_tropopause_level_mod, ONLY: nudging_calc_tropopause_level
USE nudging_call_relax_mod, ONLY: nudging_call_relax
USE nudging_interpolation_two_d_mod, ONLY: nudging_interpolation_two_d
USE nudging_nudging_three_d_mod, ONLY: nudging_nudging_three_d
USE nudging_varloader_pressure_lev_mod, ONLY: nudging_varloader_pressure_lev
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)  :: varname                 ! Variable name
INTEGER, INTENT(IN)       :: grid_type               ! Grid type
INTEGER, INTENT(IN)       :: global_row_length       ! (global) row length
INTEGER, INTENT(IN)       :: global_rows             ! (Global) rows
INTEGER, INTENT(IN)       :: file_levels             ! File levels
INTEGER, INTENT(IN)       :: proc_row_length_min     ! min column in grid
INTEGER, INTENT(IN)       :: proc_row_length_max     ! max column in grid
INTEGER, INTENT(IN)       :: proc_rows_min           ! min row in grid
INTEGER, INTENT(IN)       :: proc_rows_max           ! max row in grid
REAL, INTENT(IN)          :: sin_theta_latitude_min  ! sin min lat (T grid)
REAL, INTENT(IN)          :: sin_theta_latitude_max  ! sin max lat (T grid)
REAL, INTENT(IN)          :: frac                    ! Fraction between 1 & 2
INTEGER, INTENT(IN)       :: model_levels            ! number of model levels
REAL, INTENT(IN)  :: base_lambda             ! minimum longitude in the model
REAL, INTENT(IN)  :: delta_lambda            ! longitude increment in model
REAL, INTENT(IN)  :: base_phi                ! minimum latitude in model
REAL, INTENT(IN)  :: delta_phi               ! latitude increment in model
REAL, INTENT(IN)  :: base_lambda_ana         ! base longitude in the data
REAL, INTENT(IN)  :: delta_lambda_ana        ! longitude increment in the data
REAL, INTENT(IN)  :: base_phi_ana            ! base latitude in the data
REAL, INTENT(IN)  :: delta_phi_ana           ! latitude increment in the data

REAL, INTENT(IN)          :: variable_file_data                       &
  (global_row_length, global_rows, file_levels)

REAL, INTENT(IN)          :: variable_file_data2                      &
  (global_row_length, global_rows, file_levels)

! Analyses log surface pressure
REAL, INTENT(IN)        :: logp_surf_data                             &
  (global_row_length, global_rows)

! Analyses log surface pressure
REAL, INTENT(IN)        :: logp_surf_data2                            &
  (global_row_length, global_rows)

! Model pressure
REAL, INTENT(IN) :: model_pressure_model_levels(                      &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! tropopause pressure
REAL, INTENT(IN)  :: model_trop_pressure(                              &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max)

! Nudged variable
REAL, INTENT(INOUT) :: model_variable_model_levels(                   &
 proc_row_length_min:proc_row_length_max,                             &
 proc_rows_min:proc_rows_max,1:model_levels )

! Nudging Diagnostics (1-5)
LOGICAL, INTENT(IN) :: l_data_diag, l_model_diag, l_ntend_diag,       &
 l_mtend_diag, l_relax_diag         ! Switches for diagnostic output

REAL, INTENT(OUT) :: diag_var_data(                                   &
proc_row_length_min:proc_row_length_max,                              &
proc_rows_min:proc_rows_max,1:model_levels )

REAL, INTENT(INOUT) :: diag_var_model(                                &
proc_row_length_min:proc_row_length_max,                              &
proc_rows_min:proc_rows_max,1:model_levels )

REAL, INTENT(OUT) :: diag_tend_nudging(                               &
proc_row_length_min:proc_row_length_max,                              &
proc_rows_min:proc_rows_max,1:model_levels )

REAL, INTENT(OUT) :: diag_tend_model(                                 &
proc_row_length_min:proc_row_length_max,                              &
proc_rows_min:proc_rows_max,1:model_levels )

REAL, INTENT(OUT) :: diag_relax(                                      &
proc_row_length_min:proc_row_length_max,                              &
proc_rows_min:proc_rows_max,1:model_levels )

!--------------------
! Other variables (non IO)

! Analyses on model levels
REAL    :: data_variable_model_levels (                                &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

! Model variable before nudging
REAL    :: last_variable_model_levels (                                &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

INTEGER :: tropopause_level(                                           &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max)                     ! Relaxation Parameter

! Relaxation Parameter
REAL    :: relax_par (                                                 &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max, 1:model_levels )

! Analyses surface Pressure
REAL        :: data_logp_surf_level                                    &
  (global_row_length, global_rows)

! Analyses surface Pressure
REAL        :: data_pressure_surf_level                                &
  (global_row_length, global_rows)

! ******************************************************************
! End of Header

! If we are using hybrid pressure levels load surface pressure
IF (ndg_analysis_source == 0) THEN

  ! Interpolate log surface pressure onto timestep
  CALL nudging_interpolation_two_d(        &
    global_row_length,                     &  ! Row length
    global_rows,                           &  ! Rows
    logp_surf_data,                        &  ! Variable at start
    logp_surf_data2,                       &  ! Variable at end
    frac,                                  &  ! Frac. between start & end
    data_logp_surf_level,                  &  ! Interpolated variable
    0)                                        ! Linear Interpolation

  ! convert lnP to P
  data_pressure_surf_level=EXP(data_logp_surf_level)

END IF

!  Load data interpolated on to this model step
CALL nudging_varloader_pressure_lev( &
  varname,                           &    ! Variable name
  global_row_length,                 &    ! (global) Row length
  global_rows,                       &    ! (Global) rows
  file_levels,                       &    ! File levels
  proc_row_length_min,               &    ! Min. column in the PE
  proc_row_length_max,               &    ! Max. column in the PE
  proc_rows_min,                     &    ! Min. row in the PE
  proc_rows_max,                     &    ! Max. row in the PE
  base_lambda,                       &    ! minimum longitude in model
  delta_lambda,                      &    ! longitude increment in model
  base_phi,                          &    ! minimum latitude in model
  delta_phi,                         &    ! latitude increment in model
  base_lambda_ana,                   &    ! minimum longitude in data
  delta_lambda_ana,                  &    ! longitude increment in data
  base_phi_ana,                      &    ! minimum latitude in data
  delta_phi_ana,                     &    ! latitude increment in data
  variable_file_data,                &    ! Variable on first time step
  variable_file_data2,               &    ! Variable on second timestep
  frac,                              &    ! Fraction between  1 & 2
  model_levels,                      &    ! Model levels
  model_pressure_model_levels,       &    ! Model pressure
  data_pressure_surf_level,          &    ! Surface pressure (in analyses)
  data_variable_model_levels)             ! Variable on model levels

! Calculate tropopause level from tropopause pressure
CALL nudging_calc_tropopause_level( &
  proc_row_length_min,              &
  proc_row_length_max,              &
  proc_rows_min,                    &
  proc_rows_max,                    &
  model_levels,                     &
  model_pressure_model_levels,      &
  model_trop_pressure,              &
  tropopause_level)                  

! Load the relaxation parameter for variable
CALL nudging_call_relax(            &
  proc_row_length_min,              &
  proc_row_length_max,              &
  proc_rows_min,                    &
  proc_rows_max,                    &
  model_levels,                     &      ! Model levels
  varname,                          &      ! Variable name
  tropopause_level,                 &      ! Tropopause level
  relax_par)                               ! Relaxation parameter

! If in debug mode write out typical relaxation parameter
IF (PrintStatus > PrStatus_Oper) THEN
  WRITE(umMessage, '(4A,2E13.5)')                                              &
   ': NUDGING_NUDGING_CTL: Typical relaxation parameter',                      &
   ' for "', varname, '" is: ',                                                &
   relax_par(dbg_x,dbg_y,dbg_z), relax_par(dbg_x,dbg_y,dbg_z+5)
  CALL umPrint(umMessage,src='nudging_nudging_control')
END IF

IF (PrintStatus > PrStatus_Normal) THEN
  WRITE(umMessage,'(4A,E13.6)')                                                &
   ': NUDGING_NUDGING_CTL: Typical UNUDGED value',                             &
   ' for "', varname, '" is: ',                                                &
   model_variable_model_levels (dbg_x,dbg_y,dbg_z)
  CALL umPrint(umMessage,src='nudging_nudging_control')
END IF

! set the pre nudged model variable for diagnostic calculations
last_variable_model_levels = model_variable_model_levels

! Call function to nudge model towards data
CALL nudging_nudging_three_d(     &
  model_levels,                   &     ! Model levels
  proc_row_length_min,            &     ! Min. column in the PE
  proc_row_length_max,            &     ! Max. column in the PE
  proc_rows_min,                  &     ! Min. row in the PE
  proc_rows_max,                  &     ! Max. row in the PE
  model_variable_model_levels,    &     ! Model variable
  data_variable_model_levels,     &     ! analysis variable
  relax_par)                            ! size of timestep

IF (PrintStatus > PrStatus_Normal) THEN
  WRITE(umMessage,'(4A,E13.6)')                                                &
    ': NUDGING_NUDGING_CTL: Typical NUDGED value',                             &
    ' for "', varname, '" is: ',                                               &
    model_variable_model_levels (dbg_x,dbg_y,dbg_z)
  CALL umPrint(umMessage,src='nudging_nudging_control')
END IF

! Set nudging diagnostics
CALL nudging_calc_diags(       &
  model_levels,                &          ! Model levels
  proc_row_length_min,         &          ! Min. column in the PE
  proc_row_length_max,         &          ! Max. column in the PE
  proc_rows_min,               &          ! Min. row in the PE
  proc_rows_max,               &          ! Max. row in the PE
  last_variable_model_levels,  &          ! Model variable on last timestep
  model_variable_model_levels, &          ! Model variable
  data_variable_model_levels,  &          ! Analyses variable
  relax_par,                   &          ! Relaxation parameter (2D)
  l_data_diag, diag_var_data,      &          ! Diag 1 (variable in analyses)
  l_model_diag, diag_var_model,    &          ! Diag 2 (variable in model)
  l_ntend_diag, diag_tend_nudging, &          ! Diag 3 (nudging tend)
  l_mtend_diag, diag_tend_model,   &          ! Diag 4 (model tend)
  l_relax_diag, diag_relax)                   ! Diag 5 (relax par (3D))

RETURN
END SUBROUTINE nudging_nudging_control

END MODULE nudging_nudging_control_mod
