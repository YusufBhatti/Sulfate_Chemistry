! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_assign_vars_recon_horizontal_mod

IMPLICIT NONE

!  Subroutine rcf_assign_vars_recon_horizontal
!
! Description:
!   Having read the namelists, allocating space and assigning values
!   to the relavant bits of the Grid defined type
!
! Method:
!   Use the values read from the namelists directly from the modules
!   they're in and allocate space within the output_grid defined type
!   to store those values
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                   &
                        ModuleName='RCF_ASSIGN_VARS_RECON_HORIZONTAL_MOD'

CONTAINS

SUBROUTINE rcf_assign_vars_recon_horizontal ( )

USE Rcf_Grid_Type_Mod, ONLY:           Output_Grid

USE errormessagelength_mod, ONLY:      errormessagelength

USE model_domain_mod, ONLY: l_regular, l_cartesian, model_type, mt_global,   &
                            mt_bi_cyclic_lam, output_grid_stagger

USE rcf_readnl_horizont_mod, ONLY:     blend_weights,                     &
                                       orog_blend_width,                  &
                                       orog_blend_max,                    &
                                       orog_blend_weights

USE vrhoriz_grid_mod, ONLY: lambda_input_p, lambda_input_u,       &
    phi_input_p, phi_input_v, horizgrid, read_nml_horizgrid

USE Rcf_HeadAddress_Mod, ONLY: FH_GridStagger_C

USE lam_config_inputs_mod, ONLY: &
    polelata,  polelona,    &
    frstlona,  frstlata,    &
    delta_lon, delta_lat

USE nlsizes_namelist_mod, ONLY: &
    a_len2_rowdepc,              &
    a_len2_coldepc,              &
    var_grid

USE missing_data_mod, ONLY: rmdi

USE Ereport_Mod, ONLY: &
    Ereport

USE umprintmgr, ONLY: newline

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER                          :: errorStatus
INTEGER                          :: i
REAL                             :: small
! Dr Hook
INTEGER(KIND=jpim), PARAMETER    :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER    :: zhook_out = 1
REAL(KIND=jprb)                  :: zhook_handle

CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER       :: &
                             RoutineName = 'RCF_ASSIGN_VARS_RECON_HORIZONTAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Can't have an orography blending for a global model.
! Set width to zero.
IF (model_type == mt_global) orog_blend_width = 0

! Allocate and set the internal weights for orography blending
IF (orog_blend_width > 0) THEN

  ALLOCATE( blend_weights( orog_blend_width ), STAT=errorStatus )
  IF (errorStatus /= 0) THEN
     CALL ereport( modulename//":"//routinename, errorStatus              &
                   , 'Failure in allocating blend_weights')
  END IF

  ! blending zone weights come from namelist
  blend_weights( 1 : orog_blend_width ) = &
                 orog_blend_weights( 1 : orog_blend_width )
END IF

! Fill in gaps about Output Grid
Output_Grid % global = (model_type == mt_global)
Output_Grid % grid_stagger      = output_grid_stagger
! Set up grid u,v,p row_length and rows
Output_Grid % glob_u_row_length = Output_Grid % glob_p_row_length
Output_Grid % glob_v_row_length = Output_Grid % glob_p_row_length
Output_Grid % glob_u_rows       = Output_Grid % glob_p_rows
IF (output_grid_stagger == FH_GridStagger_C) THEN 
  Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows - 1
ELSE
  Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows + 1
END IF

! set up whether rotated or not and compensate for cyclic LAM
IF ( output_grid % global ) THEN  
  Output_Grid % rotated           = .FALSE.
ELSE    ! Atmos C grid LAM
  IF (l_cartesian) THEN  ! rotated pole make no sense
    Output_Grid % rotated           = .FALSE.
  ELSE
    IF (polelata /= 90 .OR. polelona /= 0 ) THEN
      Output_Grid % Rotated         = .TRUE.
    ELSE
      Output_Grid % Rotated         = .FALSE.
    END IF
  END IF
  IF (model_type == mt_bi_cyclic_lam) THEN
    Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows
  END IF
END IF

! set up field sizes, u,v,p
Output_Grid % glob_p_field = Output_Grid % glob_p_row_length *    &
                             Output_Grid % glob_p_rows
Output_Grid % glob_u_field = Output_Grid % glob_u_row_length *    &
                             Output_Grid % glob_u_rows
Output_Grid % glob_v_field = Output_Grid % glob_v_row_length *    &
                             Output_Grid % glob_v_rows
Output_Grid % glob_r_field = Output_Grid % glob_r_row_length *    &
                             Output_Grid % glob_r_rows

! Set some information which is not set via the namelists for global.
IF (output_grid % global) THEN
  delta_lon = 360.0/REAL( Output_Grid % glob_p_row_length )
  IF (output_grid_stagger == FH_GridStagger_C) THEN
    ! Read from lam_config but are RMDI for global jobs
    delta_lat = 180.0/REAL( Output_Grid % glob_p_rows - 1 )
  ELSE
    delta_lat = 180.0/REAL( Output_Grid % glob_p_rows )
  END IF
END IF

ALLOCATE(output_grid % lambda_p(output_grid % glob_p_row_length),        &
         STAT=errorStatus )
IF (errorStatus /= 0) THEN
   CALL ereport( modulename//":"//routinename, errorStatus              &
                 , 'Failure in allocating output_grid % lambda_p')
END IF
ALLOCATE(output_grid % lambda_u(output_grid % glob_u_row_length),        &
         STAT=errorStatus )
IF (errorStatus /= 0) THEN
   CALL ereport( modulename//":"//routinename, errorStatus              &
                 , 'Failure in allocating output_grid % lambda_u')
END IF
ALLOCATE(output_grid % phi_p   (output_grid % glob_p_rows), STAT=errorStatus )
IF (errorStatus /= 0) THEN
   CALL ereport( modulename//":"//routinename, errorStatus              &
                 , 'Failure in allocating output_grid % phi_p')
END IF
ALLOCATE(output_grid % phi_v   (output_grid % glob_v_rows), STAT=errorStatus )
IF (errorStatus /= 0) THEN
   CALL ereport( modulename//":"//routinename, errorStatus              &
                 , 'Failure in allocating output_grid % phi_v')
END IF

! ------------------------
! Horizontal Grid Namelist
! ------------------------
! Lets default the row and column dependent constants to have 0 entries
! for fixed resolution.
a_len2_rowdepc = 0
a_len2_coldepc = 0

IF (.NOT. l_regular) THEN
  ! Lets set the row dependent and column dependent sizes - originally set
  ! with the SIZES namelist but due to rationalisation we can just calculate it.
  ! Row dependent constants storing p and u
  a_len2_rowdepc = 2
  ! Column dependent constants storing p and v
  a_len2_coldepc = 2

  phi_input_v(Output_Grid % glob_v_rows + 1) = rmdi
  phi_input_p(Output_Grid % glob_p_rows + 1) = rmdi

  ! Check grid definition matches namelist sizes (similar check in CreateBC)
  IF (lambda_input_p(Output_Grid % glob_p_row_length+1) /= rmdi   .OR.   &
      lambda_input_p(Output_Grid % glob_p_row_length)   == rmdi   .OR.   &
      phi_input_p(Output_Grid % glob_p_rows+1)  /= rmdi           .OR.   &
      phi_input_p(Output_Grid % glob_p_rows)    == rmdi) THEN  
    WRITE(cmessage, '(A)') 'Dimensions of &HORIZGRID namelist'&
        // ' do not match those in &HORIZONT'
    errorstatus = 45
    CALL ereport(routinename, errorstatus, cmessage)
  END IF
  ! Set the following to BMDI for var grids model.
  delta_lon = rmdi
  delta_lat = rmdi
  frstlona  = rmdi
  frstlata  = rmdi
ELSE
  ! Endgame and New Dynamics has slightly different staggering so we need to
  ! make sure its setup correctly.
  IF (output_grid_stagger == FH_GridStagger_C) THEN
    DO i = 1, output_grid % glob_p_row_length
      lambda_input_p(i) = frstlona + (i-1)*delta_lon
    END DO
    DO i = 1, output_grid % glob_u_row_length
      lambda_input_u(i) = frstlona + (i-0.5)*delta_lon
    END DO
    DO i = 1, output_grid % glob_p_rows
      phi_input_p(i) = frstlata + (i-1)*delta_lat
    END DO
    DO i = 1, output_grid % glob_v_rows
      phi_input_v(i) = frstlata + (i-0.5)*delta_lat
    END DO
  ELSE
    DO i = 1, output_grid % glob_p_row_length
      lambda_input_p(i) = frstlona + (i-0.5)*delta_lon
    END DO
    DO i = 1, output_grid % glob_u_row_length
      lambda_input_u(i) = frstlona + (i-1)*delta_lon
    END DO
    DO i = 1, output_grid % glob_p_rows
      phi_input_p(i) = frstlata + (i-0.5)*delta_lat
    END DO
    DO i = 1, output_grid % glob_v_rows
      phi_input_v(i) = frstlata + (i-1)*delta_lat
    END DO
  END IF
END IF

! Copy grid information into reconfiguration versions rather than using the
! input variables.
output_grid % lambda_p(:) = lambda_input_p(1:output_grid % glob_p_row_length)
output_grid % lambda_u(:) = lambda_input_u(1:output_grid % glob_u_row_length)
output_grid % phi_p   (:) = phi_input_p   (1:output_grid % glob_p_rows)
output_grid % phi_v   (:) = phi_input_v   (1:output_grid % glob_v_rows)

! Prevent some unwise or inefficient LAM choices.
! When reconfiguring to an ENDGame grid the cases prevented below also avoid
! calling rcf_eg_poles and so polar rows will not be treated correctly.
! Doesn't test if variable resolution (rcf_eg_poles also assumes regular).
IF (l_regular) THEN
  IF (model_type == mt_bi_cyclic_lam) THEN
    small = EPSILON(1.0)
    ! Check if the Lat-Long bicyclic grid extends from pole to pole.
    ! If so, stop the job and use the suggested alternative domain:
    IF (ABS(delta_lat * Output_Grid % glob_p_rows - 180) < small) THEN

      IF (ABS(delta_lon * Output_Grid % glob_p_row_length - 360) < small) THEN
        WRITE(cmessage, '(A)')                                               &
          'Detected a regular bicyclic LAM that covers the whole globe.'     &
          // newline // 'Please use a global model type instead.'
        errorstatus = 90
        CALL ereport(routinename, errorstatus, cmessage)
      ELSE
        WRITE(cmessage, '(A)')                                               &
          'Detected a regular bicyclic LAM that extends from pole to pole.'  &
          // newline // 'Please use a rotated grid that avoids the poles '   &
          // 'to improve model performance.'
        errorstatus = 91
        CALL ereport(routinename, errorstatus, cmessage)
      END IF

    END IF  ! pole-to-pole
  END IF  ! bicyclic lam
END IF  ! l_regular

! Set pole coordinates
output_grid % phi_npole    = polelata
output_grid % lambda_npole = polelona

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE  rcf_assign_vars_recon_horizontal

END MODULE rcf_assign_vars_recon_horizontal_mod
