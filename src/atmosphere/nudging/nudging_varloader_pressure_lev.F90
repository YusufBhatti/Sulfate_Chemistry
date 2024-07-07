! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Load a variable from files, interpolate in time & height
!  and return on model levels & timestep

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_NUDGING_CONTROL.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_varloader_pressure_lev_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER,                 &
            PRIVATE :: ModuleName = 'NUDGING_VARLOADER_PRESSURE_LEV_MOD'
CONTAINS

SUBROUTINE nudging_varloader_pressure_lev( &
  varname,                                 & ! Variable name
  global_row_length,                       & ! Row length
  global_rows,                             & ! Rows
  file_levels,                             & ! Levels in file
  proc_row_length_min,                     & ! Min. column in the PE
  proc_row_length_max,                     & ! Max. column in the PE
  proc_rows_min,                           & ! Min. row in the PE
  proc_rows_max,                           & ! Max. row in the PE
  base_lambda,                             & ! minimum longitude in model
  delta_lambda,                            & ! longitude increment in model
  base_phi,                                & ! minimum latitude in model
  delta_phi,                               & ! latitude increment in model
  base_lambda_ana,                         & ! minimum longitude in data
  delta_lambda_ana,                        & ! longitude increment in data
  base_phi_ana,                            & ! minimum latitude in data
  delta_phi_ana,                           & ! latitude increment in data
  variable_file_data,                      & ! variable at first timestep
  variable_file_data2,                     & ! variable at second timestep
  frac,                                    & ! fraction between 1 & 2
  model_levels,                            & ! model levels
  pressure_model_levels,                   & ! model pressure
  pressure_surf_level,                     & ! surface pressure
  variable_model_levels)                     ! variable on model levels

USE nudging_control               ! Use standard nudging switches
USE nudging_ecmwf_60level_def     ! Use ECMWF model elvel information
USE nudging_jra_plevel_def        ! Use JRA pressure level
USE planet_constants_mod, ONLY: kappa, pref  ! Atmospheric constants

USE UM_ParVars
USE UM_ParParams, ONLY: pnorth

USE atm_fields_bounds_mod, ONLY: vdims
USE conversions_mod, ONLY: pi_over_180
USE trignometric_mod, ONLY: sin_v_latitude 

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE umPrintMgr
USE ereport_mod, ONLY: ereport
USE nudging_ana2mod_grid_mod, ONLY: nudging_ana2mod_grid
USE nudging_calc_ecmwf_60levs_mod, ONLY: nudging_calc_ecmwf_60levs
USE nudging_calc_ecmwf_91levs_mod, ONLY: nudging_calc_ecmwf_91levs
USE nudging_ecmwf_to_mod_mod, ONLY: nudging_ecmwf_to_mod
USE nudging_interpolation_three_d_mod, ONLY: nudging_interpolation_three_d
USE nudging_pres_to_mod_mod, ONLY: nudging_pres_to_mod
IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN)  :: varname              ! Variable name
INTEGER, INTENT(IN)       :: global_row_length    ! Row length
INTEGER, INTENT(IN)       :: global_rows ! Rows
INTEGER, INTENT(IN)       :: file_levels          ! Levels in file
INTEGER, INTENT(IN)       :: proc_row_length_min  ! Min. column in the PE
INTEGER, INTENT(IN)       :: proc_row_length_max  ! Max. column in the PE
INTEGER, INTENT(IN)       :: proc_rows_min        ! Min. row in the PE
INTEGER, INTENT(IN)       :: proc_rows_max        ! Max. row in the PE
INTEGER, INTENT(IN)       :: model_levels         ! Model levels
REAL, INTENT(IN)          :: base_lambda          ! MInimum longitude in model
REAL, INTENT(IN)          :: delta_lambda         ! longitude increment in model
REAL, INTENT(IN)          :: base_phi             ! Minimum latitude in model
REAL, INTENT(IN)          :: delta_phi            ! Latitude increment in model
REAL, INTENT(IN)          :: base_lambda_ana      ! Minimum longitude in data
REAL, INTENT(IN)          :: delta_lambda_ana     ! Longitude increment in data
REAL, INTENT(IN)          :: base_phi_ana         ! Minimum latitude in data
REAL, INTENT(IN)          :: delta_phi_ana        ! Latitude increment in data
REAL, INTENT(IN)          :: frac                 ! Frac between start & end

! Temporary Arrays used to transfer data from netcdf files to interpolation
REAL, INTENT(IN)          :: variable_file_data                        &
  (global_row_length, global_rows, file_levels)

REAL, INTENT(IN)          :: variable_file_data2                       &
  (global_row_length, global_rows, file_levels)

! Model pressure on model levels
REAL, INTENT(IN) :: pressure_model_levels(                             &
  proc_row_length_min:proc_row_length_max,                             &
  proc_rows_min:proc_rows_max, model_levels)

! Surface pressure from the data
REAL, INTENT(IN) :: pressure_surf_level                                &
 (global_row_length, global_rows)

! Output variable interpolated in time and on to model levels
REAL, INTENT(OUT) :: variable_model_levels(                            &
 proc_row_length_min:proc_row_length_max,                              &
 proc_rows_min:proc_rows_max,1:model_levels )

! Variable interpolated in time, but not in height
REAL :: variable_file_levels                                           &
  (global_row_length, global_rows, file_levels)

! Used in flipping analyses from top down to bottom up
! Legacy of sawitching from data obtained at the BADC to data from BDAN
REAL :: variable_file_levels_temp                                     &
  (global_row_length, global_rows, file_levels)

! Variable on ECMWF hybrid p-levels and model grid
REAL             :: variable_interp(                                  &
  proc_row_length_min:proc_row_length_max,                            &
  proc_rows_min:proc_rows_max, 1:file_levels)

! Pressure on ECMWF hybrid p-levels and anlysis grid
REAL             :: pressure_ecmwf_levels                              &
  (global_row_length, global_rows, file_levels)

! Pressure on ECMWF hybrid p-levels and model grid
REAL             :: pressure_interp(                                   &
  proc_row_length_min:proc_row_length_max,                             &
  proc_rows_min:proc_rows_max, 1:file_levels)

! Variables for handling v_at_poles
LOGICAL, SAVE    :: l_first = .TRUE.
REAL, ALLOCATABLE, SAVE :: v_lat(:) ! Latitude points
REAL             :: max_ana_lat, rfac, diff_lat 
                                    ! Max data lat, interp factor
REAL, SAVE       :: delta_lat       ! Delta_phi
INTEGER,SAVE     :: last_lat        ! Last model lat within analyses data 

INTEGER :: i, j, k, l, a, b         ! Loop variables
INTEGER :: dummy

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_VARLOADER_PRESSURE_LEV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! End of Header
! *****************************************************************************

! Initialise data to zero
variable_file_levels = 0

! Interpolate between the two input data files onto the model timestep
CALL nudging_interpolation_three_d(  &
  global_row_length,                 & ! Row length
  global_rows,                       & ! Rows
  file_levels,                       & ! File levels
  variable_file_data(:,:,:),         & ! Variable at start
  variable_file_data2(:,:,:),        & ! Variable at end
  frac,                              & ! Frac between start & end
  variable_file_levels,              & ! Interpolated variable
  0)                                   ! Linear interpolation

! If we know the data is in top down format flip it
! This is a legacy of BADC data which is bottom up
! Converted to BDAN data which is top down, so have to flip
IF (data_source == 0) THEN
  DO k=1, file_levels
    l = (file_levels+1) - k
    variable_file_levels_temp(:,:,l) = variable_file_levels(:,:,k)
  END DO
  DO k=1, file_levels
    variable_file_levels(:,:,k) = variable_file_levels_temp(:,:,k)
  END DO
END IF

! Depending on level structure of analyses source choose a method
! for converting from analyses levels to model levels
SELECT CASE(ndg_analysis_source)

  ! ECMWF hybrid pressure levels
CASE (0)

  ! Calculate levels on basis
  SELECT CASE(file_levels)

  CASE (60)

    ! Calculate the ECMWF hybrid pressure levels
    CALL nudging_calc_ecmwf_60levs(   &
      global_row_length,              & ! Length of model row
      global_rows,                    & ! Number of model rows
      file_levels,                    & ! Number of file levels
      pressure_surf_level,            & ! Surface pressure
      pressure_ecmwf_levels)            ! Return hybrid p levels

  CASE (91)

    ! Calculate the ECMWF hybrid pressure levels
    CALL nudging_calc_ecmwf_91levs(  &
      global_row_length,             & ! Length of model row
      global_rows,                   & ! Number of model rows
      file_levels,                   & ! Number of file levels
      pressure_surf_level,           & ! Surface pressure
      pressure_ecmwf_levels)           ! Return hybrid p levels

  CASE DEFAULT
    nmessage='Unknown number of levels in ECMWF Data'
    dummy = 999

    CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                  &
                 dummy,nmessage)
  END SELECT     ! file levels

  ! Convert temperature to potential temperature
  IF (varname == temp_name)                                             &
    variable_file_levels = variable_file_levels*                       &
                          (pref/pressure_ecmwf_levels)**kappa

  ! Convert from analysis to model grid
  CALL nudging_ana2mod_grid(       &
    global_row_length,             &  ! Global length of rows (in data)
    global_rows,                   &  ! Global numver of rows (in data)
    file_levels,                   &  ! Number of levels in data
    proc_row_length_min,           &  ! min column in model pe
    proc_row_length_max,           &  ! max column in model pe
    proc_rows_min,                 &  ! min row in model pe
    proc_rows_max,                 &  ! max row in model pe
    variable_file_levels,          &  ! Variable on data grid
    variable_interp,               &  ! Variable on model grid
    base_lambda,                   &  ! Min longitude in model
    delta_lambda,                  &  ! Longitude increment in model
    base_phi,                      &  ! Min latitude in model
    delta_phi,                     &  ! Latitude increment in model
    base_lambda_ana,               &  ! Min longitude in data
    delta_lambda_ana,              &  ! Longitude increment in data
    base_phi_ana,                  &  ! Min latitude in data
    delta_phi_ana)                    ! Latitude increment in data

  CALL nudging_ana2mod_grid(       &
    global_row_length,             &  ! Global length of rows (in data)
    global_rows,                   &  ! Global numver of rows (in data)
    file_levels,                   &  ! Number of levels in data
    proc_row_length_min,           &  ! min column in model pe
    proc_row_length_max,           &  ! max column in model pe
    proc_rows_min,                 &  ! min row in model pe
    proc_rows_max,                 &  ! max row in model pe
    pressure_ecmwf_levels,         &  ! Pressure on data grid
    pressure_interp,               &  ! Pressure on model grid
    base_lambda,                   &  ! Min longitude in model
    delta_lambda,                  &  ! Longitude increment in model
    base_phi,                      &  ! Min latitude in model
    delta_phi,                     &  ! Latitude increment in model
    base_lambda_ana,               &  ! Min longitude in data
    delta_lambda_ana,              &  ! Longitude increment in data
    base_phi_ana,                  &  ! Min latitude in data
    delta_phi_ana)                    ! Latitude increment in data

  ! Interpolate the analyses ftom ECMWF hybrid pressure levels
  ! to UM hybrid height levels
  CALL nudging_ecmwf_to_mod(       &
    global_row_length,             &  ! Row length
    global_rows,                   &  ! Rows
    model_levels,                  &  ! ECMWF hybrid levels
    file_levels,                   &  ! Model levels
    proc_row_length_min,           &  ! Min. column in the PE
    proc_row_length_max,           &  ! Max. column in the PE
    proc_rows_min,                 &  ! Min. row in the PE
    proc_rows_max,                 &  ! Max. row in the PE
    pressure_model_levels,         &  ! Pressure on model levels
    pressure_interp,               &  ! pressure on ECMWF hybrid levels
    variable_interp,               &  ! Variable on ecmwf levels
    variable_model_levels)            ! Variable on  model levels

  !variable_model_levels = pressure_interp

! For ENDGAME/ V_AT_POLES, last row at N pole (90.00N) is outside ECMWF
!  data (88.75N at N48), so extrapolate from previous two rows.
!  Very specific case, only for V. Also may need to be reviewed for 
!  higher resolution
 
  IF ( varname == v_name ) THEN
    IF (at_extremity(PNorth) ) THEN   
         
! Determine v_latitude points --one time
! ***Note***: assumes an uniform lat/long grid i.e. all longitude
!             points along a lat row have same latitude coordinates
 
      IF ( l_first ) THEN
        ALLOCATE(v_lat(proc_rows_min:proc_rows_max))
        v_lat(:)   = ASIN(sin_v_latitude(1,:))   
     
       ! Max v latitude point in analyses (radians)
        max_ana_lat = ana_base_phi_v + (ana_rows-2) * ana_delta_phi 
        max_ana_lat = max_ana_lat * pi_over_180
 
      ! Determine points lying beyond max analyses latitude
        last_lat = -1
        loop_lat: DO j = proc_rows_max, proc_rows_min, -1
          IF ( v_lat(j) <= max_ana_lat ) THEN
            last_lat = j          ! Reached max interp requirement 
            delta_lat = v_lat(last_lat) - v_lat(last_lat-1)
            EXIT loop_lat
          END IF
        END DO loop_lat

        l_first = .FALSE. 
      END IF           ! first
 
      IF ( last_lat > 0 ) THEN ! Not req if all rows in domain
 
! Extrapolate analysis V onto model V_at_pole
        DO j = last_lat+1, proc_rows_max
          diff_lat = v_lat(j) - v_lat(j-1) 
          DO k = 1, model_levels
           DO i = proc_row_length_min, proc_row_length_max
 
         ! Determine change_per_deg for last two 'within data' rows
             rfac = (variable_model_levels(i,last_lat,k)             &
                 - variable_model_levels(i,last_lat-1,k))/ delta_lat
 
         ! Apply factor to extrapolate
             variable_model_levels(i,j,k) =                          &
                  variable_model_levels(i,j-1,k) + (rfac * diff_lat)
 
           END DO     ! longitudes
          END DO      ! levels
         END DO     ! outlier rows
 
      END IF        ! interp required ( model in analyses domain)?
 
    END IF         ! 
  END IF           ! varname = V
 
! End special handling of v_at_poles

! ECMWF fixed pressure levels
CASE (1)

  IF (file_levels /= ecmwf_nplevs) THEN
    nmessage = 'Number of levels in file differ from standard'
    dummy = 999

    CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV', dummy, nmessage)
  END IF

  ! Convert temperature to potential temperature
  IF (varname == temp_name) THEN
    DO k=1, file_levels
      variable_file_levels(:,:,k) =   variable_file_levels(:,:,k)*     &
                                    (pref/data_pressure(k))**kappa
    END DO
  END IF

  ! Interpolate analyses onto UM grid

  ! Drop a line to let us know where we are
  IF (PrintStatus > PrStatus_Diag) THEN
    WRITE(umMessage, '(2A,E13.6)')                                      &
        ': NUDGING_VARLOADER_PRESSURE_LEV: Finished Interpolation,',        &
        'Typical "t" value equals ',variable_file_levels(dbg_x,dbg_y,dbg_z)
    CALL umPrint(umMessage,src='nudling_varloader_pressure_lev')
  END IF
  ! Interpolate analyses from ECMWF fixed p-levels to model levels
  CALL nudging_pres_to_mod(        &
    global_row_length,             & ! Row length
    global_rows,                   & ! Rows
    model_levels,                  & ! ECMWF fixed p levels
    file_levels,                   & ! model levels
    proc_row_length_min,           & ! Min. column in the PE
    proc_row_length_max,           & ! Max. column in the PE
    proc_rows_min,                 & ! Min. row in the PE
    proc_rows_max,                 & ! Max. row in the PE
    pressure_model_levels,         & ! Pressure on model levels
    data_pressure,                 & ! ECMWF fixed P levels
    variable_file_levels,          & ! Variable on ecmwf levels
    variable_model_levels)           ! Variable on model levels

  ! If variable is on UM hybrid height levels then trivial
CASE (2)

  ! Convert from analysis to model grid
  CALL nudging_ana2mod_grid(       &
   global_row_length,             &  ! Global length of rows (in data)
   global_rows,                   &  ! Global numver of rows (in data)
   file_levels,                   &  ! Number of levels in data
   proc_row_length_min,           &  ! min column in model pe
   proc_row_length_max,           &  ! max column in model pe
   proc_rows_min,                 &  ! min row in model pe
   proc_rows_max,                 &  ! max row in model pe
   variable_file_levels,          &  ! Variable on data grid
   variable_interp,               &  ! Variable on model grid
   base_lambda,                   &  ! Min longitude in model
   delta_lambda,                  &  ! Longitude increment in model
   base_phi,                      &  ! Min latitude in model
   delta_phi,                     &  ! Latitude increment in model
   base_lambda_ana,               &  ! Min longitude in data
   delta_lambda_ana,              &  ! Longitude increment in data
   base_phi_ana,                  &  ! Min latitude in data
   delta_phi_ana)                    ! Latitude increment in data

  variable_model_levels(:,:,:) = 0.0

  ! As on same levels no interpolation  required
  DO j= proc_rows_min, proc_rows_max
    DO i= proc_row_length_min, proc_row_length_max
      variable_model_levels(i,j, ndg_lev_bottom:ndg_lev_top) =          &
            variable_interp(i,j, ndg_lev_bottom:ndg_lev_top)
    END DO
  END DO

  ! Convert temperature to potential temperature
  IF (varname == temp_name) THEN
    DO k=ndg_lev_bottom, ndg_lev_top
      DO j = proc_rows_min, proc_rows_max
        DO i = proc_row_length_min, proc_row_length_max
          variable_model_levels(i,j,k) =                                   &
           variable_model_levels(i,j,k)*                                   &
           (pref/pressure_model_levels(i,j,k))**kappa
        END DO
      END DO
    END DO
  END IF

  ! Variable on JRA fixed Pressure levels
CASE (3)

  IF (file_levels /= jra_nplevs) THEN
    nmessage = 'Number of levels in file differ from standard'
    dummy = 999

    CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                     &
               dummy, nmessage)
  END IF

  ! Convert temperature to potential temperature
  IF (varname == temp_name) THEN
    DO k=1, file_levels
      variable_file_levels(:,:,k) =                                    &
       variable_file_levels(:,:,k)*                                    &
       (pref/jra_pressure(k))**kappa
    END DO
  END IF

  ! Interpolate analyses onto UM grid

  ! Drop a line to let us know where we are
  IF (PrintStatus > PrStatus_Diag) THEN
    WRITE(umMessage,'(2A,E13.6)')                                       &
        ': NUDGING_VARLOADER_PRESSURE_LEV: Finished Interpolation,',    &
        'Typical "t" value equals ',variable_file_levels(dbg_x,dbg_y,dbg_z)
    CALL umPrint(umMessage,src='nudling_varloader_pressure_lev')
  END IF
  ! Interpolate analyses from JRA fixed P levels to model levels
  CALL nudging_pres_to_mod(          &
    global_row_length,               & ! Row length
    global_rows,                     & ! Rows
    model_levels,                    & ! JRA levels
    file_levels,                     & ! model levels
    proc_row_length_min,             & ! Min. column in the PE
    proc_row_length_max,             & ! Max. column in the PE
    proc_rows_min,                   & ! Min. row in the PE
    proc_rows_max,                   & ! Max. row in the PE
    pressure_model_levels,           & ! Pressure on model levels
    jra_pressure,                    & ! JRA pressure levels
    variable_file_levels,            & ! Variable on JRA levels
    variable_model_levels)             ! Variable on model levels

CASE DEFAULT

  nmessage = 'Unkown Analysis data source'
  dummy = 999

  CALL ereport('NUDGING_VARLOADER_PRESSURE_LEV',                       &
               dummy, nmessage)

END SELECT         ! Pressure level scheme

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_varloader_pressure_lev

END MODULE nudging_varloader_pressure_lev_mod
