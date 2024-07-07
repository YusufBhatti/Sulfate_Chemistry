! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE fieldsfile_mod

USE datafile_mod,                 ONLY: datafile_type, new_dynamics, endgame
USE horizontal_lat_long_grid_mod, ONLY: horizontal_lat_long_grid_type
USE vertical_grid_mod,            ONLY: vertical_grid_type
USE missing_data_mod,             ONLY: imdi, rmdi
USE data_location_mod,            ONLY: data_location_type
USE field_mod,                    ONLY: field_type, ASSIGNMENT(=)
USE three_dimensional_grid_mod,   ONLY: three_dimensional_grid_type

USE io, ONLY: setpos, buffin, buffout, file_open, file_close
USE io_constants, ONLY: ioNameProvided, ioOpenReadWrite, ioOpenReadOnly

USE um_types, ONLY: real64, real32
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE conversions_mod,        ONLY: isec_per_day, isec_per_hour, isec_per_min
USE lbc_output_control_mod, ONLY: packing_unpacked_accuracy

USE Packing_Codes_Mod, ONLY: PC_WGDOS_Packing, PC_No_Packing,                  &
                             PC_Cray32_Packing, PC_BitMask_CompressType,       &
                             PC_LandMask_Compression

USE lookup_addresses, ONLY: lbpack

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Description:
!   A class representing fieldsfiles. This contains all the fieldsfile-specific
!   information, such as length of headers, and methods on how to read/write
!   the file format. See documentation UMDP F3 for further details on the FF
!   format.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FIELDSFILE_MOD'

TYPE, PUBLIC, EXTENDS(datafile_type) :: fieldsfile_type

  INTEGER, PUBLIC :: read_write_status = imdi

  TYPE(data_location_type), ALLOCATABLE :: data_locations(:)

  INTEGER :: vertical_grid_indicator
  INTEGER :: horizontal_grid_indicator
  INTEGER :: dataset_type

  ! Is the file on a variable resolution grid?
  LOGICAL :: variable_resolution

  INTEGER :: source_model_version
    
  TYPE(vertical_grid_type) :: file_rhcrit
  TYPE(vertical_grid_type) :: file_soil_thickness
  TYPE(vertical_grid_type) :: file_zsea_theta
  TYPE(vertical_grid_type) :: file_c_theta
  TYPE(vertical_grid_type) :: file_zsea_rho
  TYPE(vertical_grid_type) :: file_c_rho
  
  INTEGER :: time_1(7), time_2(7), time_3(7)

  INTEGER :: num_soil_levels = imdi
  INTEGER :: num_cloud_levels = imdi
  INTEGER :: num_tracer_levels = imdi
  INTEGER :: num_boundary_layer_levels = imdi
  INTEGER :: algorithm_to_generate_height_fields = imdi
  INTEGER :: num_land_points = imdi
  INTEGER :: num_soil_moisture_levels = imdi
 
  INTEGER :: num_2d_fields_in_file = 0
  
  INTEGER :: start_of_file = 0
  INTEGER :: integer_constants_start = 0
  INTEGER :: integer_constants_dimension = 0
  INTEGER :: real_constants_start = 0
  INTEGER :: real_constants_dimension = 0
  INTEGER :: level_dep_constants_start = 0
  INTEGER :: level_dep_constants_dimension1 = 0
  INTEGER :: level_dep_constants_dimension2 = 0
  INTEGER :: row_dep_constants_start = 0
  INTEGER :: row_dep_constants_dimension1 = 0
  INTEGER :: row_dep_constants_dimension2 = 0
  INTEGER :: column_dep_constants_start = 0
  INTEGER :: column_dep_constants_dimension1 = 0
  INTEGER :: column_dep_constants_dimension2 = 0
  INTEGER :: lookup_start = 0
  ! Set default for number of lookups, can overide from namelist
  INTEGER :: num_reserved_headers = 4096 

  INTEGER :: data_start = 0
  INTEGER :: data_length = 0

  INTEGER :: most_recent_data_position = imdi
  
  INTEGER :: fixed_header_length = 256
  INTEGER :: len_integer_constants = 46
  INTEGER :: len_real_constants = 38
  INTEGER :: len2_lev_dep_constants = 8
  INTEGER :: len2_row_dep_constants = 2
  INTEGER :: len2_col_dep_constants = 2
  INTEGER :: len_single_lookup = 64
  
  CONTAINS
  PROCEDURE, PASS :: read_header
  PROCEDURE, PASS :: read_field
  PROCEDURE, PASS :: write_header
  PROCEDURE, PASS :: write_field
  PROCEDURE, PASS :: open_file
  PROCEDURE, PASS :: close_file
  PROCEDURE, PASS :: add_field
  PROCEDURE, PASS :: sort_fields
  PROCEDURE, PASS :: field
  PROCEDURE, PASS :: clear_data_locations
  PROCEDURE, PASS :: allow_read_only
  PROCEDURE, PASS :: number_of_dust_bins_in_file
  PROCEDURE, PASS :: check_source_cyclic
  PROCEDURE, PASS :: process_times
  PROCEDURE, PASS :: set_header_sizes
  PROCEDURE, PASS :: number_of_unique_stash_items
  PROCEDURE, PASS :: number_of_tracers
  PROCEDURE, PASS :: copy_common_metadata
  PROCEDURE, PASS :: validate_grid
  PROCEDURE, PASS :: modify_packing_method
END TYPE fieldsfile_type

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE read_header(this)
USE field_mod, ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_utils_mod, ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: vertical_coord_hybrid_height,         &
                                     horiz_grid_type, vert_level_type,     &
                                     u_points, v_points, p_points,         &
                                     theta_levels, rho_levels,             &
                                     first_level_code,                     &
                                     last_level_code, halo_type_code
USE time_utils_mod, ONLY: set_calendar
USE umPrintMgr, ONLY: umPrint, umMessage, newline

IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this

! Variables containing the header read from disk
INTEGER              :: fixed_header(this%fixed_header_length)
INTEGER, ALLOCATABLE :: integer_constants(:)
REAL, ALLOCATABLE    :: real_constants(:)
REAL, ALLOCATABLE    :: level_dep_constants(:,:)
REAL, ALLOCATABLE    :: row_dep_constants(:,:)
REAL, ALLOCATABLE    :: column_dep_constants(:,:)
INTEGER, ALLOCATABLE :: lookup(:,:)
LOGICAL, ALLOCATABLE :: assigned_lookup(:)

TYPE(data_location_type), ALLOCATABLE :: local_data_locations(:)
TYPE(field_type), ALLOCATABLE :: local_field(:)

REAL :: real_part_of_lookup(19)
REAL :: startlat, startlong
INTEGER :: nlat, nlong
INTEGER :: run_ident
INTEGER :: expt_ident
INTEGER :: calendar
INTEGER :: source_model_version
REAL    :: height_at_top_theta_level
INTEGER :: first_rho_of_constant_height
INTEGER :: lookup_dimension1 = 64
INTEGER :: lookup_dimension2 = 4096

INTEGER, ALLOCATABLE :: local_levels(:)
INTEGER :: lookup_items_in_common(23)
INTEGER :: temp_hemisphere_indicator

CHARACTER(LEN=32) :: reason ! Explanation of LBREL difference
LOGICAL :: reset_seconds ! Flag to reset seconds for reading in LBREL=2 files
LOGICAL :: printed_lbrel_warning ! Flag to avoid spamming user with warnings
  
! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'READ_HEADER'

! Loop iterators/counters
INTEGER :: i, j, k, l, m, p
! These are used for the following:
!  i - initial loop over the lookup contents
!  j - counter for the number of 3D fields found
!  k - counter for the number of 2D fields making up the current 3D field 'j'
!  l - loop over lookup items which must be identical for a 2D field to be a
!         member of a 3D field.
!  m - loop over remaining lookup contents to search for other members of this
!         3D field.
!  p - loop over real header items to transfer from integer to real

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'ff:read_header', 5)
fixed_header = imdi
! The following lookup positions must match for fields to be considered the
! same except for level
lookup_items_in_common = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, &
                          16, 17, 18, 19, lbpack, 25, 26, 38, 42 ]

! Read in fixed header
CALL setpos(this%unit_num, this%start_of_file)
CALL buffin(this%unit_num, fixed_header, this%fixed_header_length)

this%vertical_grid_indicator   = fixed_header(3)
this%horizontal_grid_indicator = fixed_header(4)
CALL this%check_source_cyclic()
this%dataset_type         = fixed_header(5)
run_ident                 = fixed_header(6)
expt_ident                = fixed_header(7)
calendar                  = fixed_header(8)
CALL set_calendar(calendar)
this%grid_staggering      = fixed_header(9)
this%source_model_version = fixed_header(12)
this%time_1(1:7)          = fixed_header(21:27)
this%time_2(1:7)          = fixed_header(28:34)
this%time_3(1:7)          = fixed_header(35:41)

this%integer_constants_start         = fixed_header(100)
this%integer_constants_dimension     = fixed_header(101)
this%real_constants_start            = fixed_header(105)
this%real_constants_dimension        = fixed_header(106)
this%level_dep_constants_start       = fixed_header(110)
this%level_dep_constants_dimension1  = fixed_header(111)
this%level_dep_constants_dimension2  = fixed_header(112)
this%row_dep_constants_start         = fixed_header(115)
this%row_dep_constants_dimension1    = fixed_header(116)
this%row_dep_constants_dimension2    = fixed_header(117)
this%column_dep_constants_start      = fixed_header(120)
this%column_dep_constants_dimension1 = fixed_header(121)
this%column_dep_constants_dimension2 = fixed_header(122)
  
this%lookup_start = fixed_header(150)
lookup_dimension1 = fixed_header(151)
lookup_dimension2 = fixed_header(152)

this%data_start = fixed_header(160)
this%data_length = fixed_header(161)

! Read in integer and real constants
ALLOCATE(integer_constants(this%integer_constants_dimension))
integer_constants = imdi
CALL setpos(this%unit_num, this%integer_constants_start-1)
CALL buffin(this%unit_num, integer_constants, this%integer_constants_dimension)

ALLOCATE(real_constants(this%real_constants_dimension))
real_constants = rmdi
CALL setpos(this%unit_num, this%real_constants_start-1)
CALL buffin(this%unit_num, real_constants, this%real_constants_dimension)

this%num_model_levels                    = integer_constants(8)
this%num_soil_levels                     = integer_constants(10)
this%num_cloud_levels                    = integer_constants(11)
this%num_tracer_levels                   = integer_constants(12)
this%num_boundary_layer_levels           = integer_constants(13)
this%algorithm_to_generate_height_fields = integer_constants(17)
first_rho_of_constant_height             = integer_constants(24)
this%num_land_points                     = integer_constants(25)
this%num_soil_moisture_levels            = integer_constants(28)
this%pole_lat                            = real_constants(5)
this%pole_long                           = real_constants(6)
height_at_top_theta_level                = real_constants(16)

! Level-dependent constants and the corresponding vertical level information
ALLOCATE(level_dep_constants(this%level_dep_constants_dimension1,         &
                             this%level_dep_constants_dimension2))
level_dep_constants = rmdi
CALL setpos(this%unit_num, this%level_dep_constants_start-1)
CALL buffin(this%unit_num, level_dep_constants)

! Check that the number of model levels in the integer constants matches 
! the size of the level dependent constants
IF (this%num_model_levels + 1 /= this%level_dep_constants_dimension1) THEN
  icode = 10
  WRITE(cmessage, '(A,I8,A,I8)')   &
               "Inconsistency between level dependent constants dimension : ", &
               this%integer_constants_dimension,                               &
               " and model levels + 1 from fixed header : ",                   &
               this%num_model_levels + 1
  CALL ereport(routinename, icode, cmessage)
END IF

CALL this%file_theta_rho_levels%set_defined_vertical_grid(                     &
               integer_constants(8),                                           &
               num_model_levels              = integer_constants(8),           &
               height_at_top_theta_level     = height_at_top_theta_level,      &
               first_rho_of_constant_height  = first_rho_of_constant_height,   &
               list_eta_theta                =                                 &
                             level_dep_constants(1:this%num_model_levels+1,1), &
               list_eta_rho                  =                                 &
                               level_dep_constants(1:this%num_model_levels,2))

CALL this%file_rhcrit%set_defined_vertical_grid(                          &
                                        integer_constants(8)+1,           &
                                        list_z=level_dep_constants(:,3))

CALL this%file_soil_thickness%set_defined_vertical_grid(                  &
                                        this%num_soil_levels,             &
                                        list_z=level_dep_constants(:,4))

CALL this%file_zsea_theta%set_defined_vertical_grid(                      &
                                        integer_constants(8)+1,           &
                                        list_z=level_dep_constants(:,5))

CALL this%file_c_theta%set_defined_vertical_grid(                         &
                                        integer_constants(8)+1,           &
                                        list_z=level_dep_constants(:,6))

CALL this%file_zsea_rho%set_defined_vertical_grid(                        &
                                        integer_constants(8)+1,           &
                                        list_z=level_dep_constants(:,7))

CALL this%file_c_rho%set_defined_vertical_grid(                           &
                                        integer_constants(8)+1,           &
                                        list_z=level_dep_constants(:,8))


CALL this%p_grid%vert_grid%set_defined_vertical_grid(                         &
                integer_constants(8),                                         &
                num_model_levels             = integer_constants(8),          &
                height_at_top_theta_level    = height_at_top_theta_level,     &
                first_rho_of_constant_height = first_rho_of_constant_height,  &
                list_eta_theta               =                                &
                           level_dep_constants(1:this%num_model_levels+1,1),  &
                list_eta_rho                 =                                &
                               level_dep_constants(1:this%num_model_levels,2)) 

CALL this%u_grid%vert_grid%set_defined_vertical_grid(                         &
                    integer_constants(8),                                     &
                    num_model_levels             = integer_constants(8),      &
                    height_at_top_theta_level    = height_at_top_theta_level, &
                    first_rho_of_constant_height =                            &
                                           first_rho_of_constant_height,      &
                    list_eta_theta               =                            &
                           level_dep_constants(1:this%num_model_levels+1,1),  &
                    list_eta_rho                 =                            &
                               level_dep_constants(1:this%num_model_levels,2))

CALL this%v_grid%vert_grid%set_defined_vertical_grid(                         &
                    integer_constants(8),                                     &
                    num_model_levels             = integer_constants(8),      &
                    height_at_top_theta_level    = height_at_top_theta_level, &
                    first_rho_of_constant_height =                            &
                                                first_rho_of_constant_height, &
                    list_eta_theta               =                            &
                           level_dep_constants(1:this%num_model_levels+1,1),  &
                    list_eta_rho                 =                            &
                               level_dep_constants(1:this%num_model_levels,2))
  
! Set the hemisphere/domain indicator for the grid definitions contained
! in the file object The grid indicator in fixed header 4 is nearly the same
! as lookup 17 grid indicator BUT the lookup 17 value does not have 100 added
! for rotated grid
temp_hemisphere_indicator = fixed_header(4)
IF (temp_hemisphere_indicator >= 100) THEN
    temp_hemisphere_indicator = temp_hemisphere_indicator - 100
END IF
CALL this%p_grid%set_hemisphere_indicator(temp_hemisphere_indicator)
CALL this%u_grid%set_hemisphere_indicator(temp_hemisphere_indicator)
CALL this%v_grid%set_hemisphere_indicator(temp_hemisphere_indicator)
! Set the grid types
CALL this%p_grid%set_horiz_grid_type(p_points)
CALL this%u_grid%set_horiz_grid_type(u_points)
CALL this%v_grid%set_horiz_grid_type(v_points)

! Allocate row and column dependent constants if present (var-res grid)
IF (this%row_dep_constants_dimension1 > 0 .AND.                           &
    this%row_dep_constants_dimension2 > 0 .AND.                           &
    this%column_dep_constants_dimension1 > 0 .AND.                        &
    this%column_dep_constants_dimension2 > 0                     ) THEN  
  ALLOCATE(row_dep_constants(this%row_dep_constants_dimension1,           &
                               this%row_dep_constants_dimension2))
  row_dep_constants = rmdi
  CALL setpos(this%unit_num, this%row_dep_constants_start-1)
  CALL buffin(this%unit_num, row_dep_constants)

  ALLOCATE(column_dep_constants(this%column_dep_constants_dimension1,     &
                                this%column_dep_constants_dimension2))
  column_dep_constants = rmdi
  CALL setpos(this%unit_num, this%column_dep_constants_start-1)
  CALL buffin(this%unit_num, column_dep_constants)

  ! p-grid sizes    
  nlat  = integer_constants(7)
  nlong = integer_constants(6)

  CALL this%p_grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = nlong,                            &
                ny     = nlat,                             &
                list_x = column_dep_constants(1:nlong,1),  &
                list_y = row_dep_constants(1:nlat,1) )

  CALL this%u_grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = nlong,                            &
                ny     = nlat,                             &
                list_x = column_dep_constants(1:nlong,2),  &
                list_y = row_dep_constants(1:nlat,1) )

  ! v-grid size depends on whether we're Endgame and have an extra v-row 
  ! than we do in the p-grid, or New Dynamics which has one fewer.
  IF (this%grid_staggering == endgame) THEN
    CALL this%v_grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = nlong,                              &
                ny     = nlat + 1,                           &
                list_x = column_dep_constants(1:nlong,1),    &
                list_y = row_dep_constants(1:nlat+1,2) )
  ELSE
    CALL this%v_grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = nlong,                              &
                ny     = nlat - 1,                           &
                list_x = column_dep_constants(1:nlong,1),    &
                list_y = row_dep_constants(1:nlat-1,2) )

  END IF

  this%variable_resolution = .TRUE.
ELSE
  !     
  !     NEW DYNAMICS 4 x 3 grid                 ENDGAME 4 x 2 grid
  !
  !     Grid decompositions labels are based on the number of P points
  !  ________________________________     ________________________________
  !  |P   U   P   U   P   U   P   U |     |    V       V       V       V | 
  !  |                              |     |                              | 
  !  |V       V       V       V     |     |U   P   U   P   U   P   U   P | 
  !  |                              |     |                              | 
  !  |P   U   P   U   P   U   P   U |     |    V       V       V       V | 
  !  |                              |     |                              | 
  !  |V       V       V       V     |     |U   P   U   P   U   P   U   P | 
  !  |                              |     |                              | 
  !  |P   U   P   U   P   U   P   U |     |    V       V       V       V | 
  !  --------------------------------     -------------------------------- 
  !  
  !  i) Note that both these example grids cover the same area, and have the
  !  same grid start lat/long.  The start lat/long of the grid is the latitude
  !  of the first row and the longitude of the first column.  This is a P
  !  point for ND but lies between U and V points for EG.
  ! 
  !  ii) In order to cover the same size area the EG grid requires one less P
  !  row but has an extra V row.
    
  ! Calculate P grid start lat/long offsets from general grid start lat/long
  nlat  = integer_constants(7)
  nlong = integer_constants(6)
  IF (this%grid_staggering == endgame) THEN
    startlat  = real_constants(3) + 0.5 * real_constants(2)
    startlong = real_constants(4) + 0.5 * real_constants(1)
  ELSE 
    startlat  = real_constants(3)
    startlong = real_constants(4)
  END IF
  CALL this%p_grid%horiz_grid%set_regular_horizontal_grid( &
                nx     = nlong,                            &
                ny     = nlat,                             &
                startx = startlong,                        &
                starty = startlat,                         &
                dx     = real_constants(1),                &
                dy     = real_constants(2)  )

  ! Calculate U grid start lat/long offsets from general grid start lat/long
  IF (this%grid_staggering == endgame) THEN
    startlat  = real_constants(3) + 0.5 * real_constants(2)
    startlong = real_constants(4)
  ELSE
    startlat  = real_constants(3)
    startlong = real_constants(4) + 0.5 * real_constants(1)    
  END IF
  CALL this%u_grid%horiz_grid%set_regular_horizontal_grid( &
                nx     = nlong,                            &
                ny     = nlat,                             &
                startx = startlong,                        &
                starty = startlat,                         &
                dx     = real_constants(1),                &
                dy     = real_constants(2)  )

  ! Calculate V grid start lat/long offsets from general grid start lat/long
  IF (this%grid_staggering == endgame) THEN
    startlat  = real_constants(3)
    startlong = real_constants(4) + 0.5 * real_constants(1)
    nlat = integer_constants(7) + 1
  ELSE
    startlat  = real_constants(3) + 0.5 * real_constants(2)
    startlong = real_constants(4)
    nlat = integer_constants(7) - 1
  END IF
  CALL this%v_grid%horiz_grid%set_regular_horizontal_grid( &
                nx     = nlong,                            &
                ny     = nlat,                             &
                startx = startlong,                        &
                starty = startlat,                         &
                dx     = real_constants(1),                &
                dy     = real_constants(2)  )

  this%variable_resolution = .FALSE.
END IF

CALL this%p_grid%horiz_grid%set_pole(real_constants(5), real_constants(6))
CALL this%u_grid%horiz_grid%set_pole(real_constants(5), real_constants(6))
CALL this%v_grid%horiz_grid%set_pole(real_constants(5), real_constants(6))
IF ( real_constants(5) == 90.0 .AND. real_constants(6) == 0.0 ) THEN
  this%l_source_rotated = .FALSE.
ELSE
  this%l_source_rotated = .TRUE.
END IF

IF (lookup_dimension1 /= 64) THEN
  icode = 1
  WRITE(cmessage, '(A,I8,A)')                                             &
             'Only 64-word headers are supported, detected ',             &
             lookup_dimension1, ' words'
  CALL ereport(routinename, icode, cmessage)
END IF

! Read in the lookup table
CALL setpos(this%unit_num, this%lookup_start-1)
ALLOCATE(lookup(lookup_dimension1, lookup_dimension2))
ALLOCATE(assigned_lookup(lookup_dimension2))
ALLOCATE(local_field(lookup_dimension2))
ALLOCATE(local_data_locations(lookup_dimension2))
CALL buffin(this%unit_num, lookup)

printed_lbrel_warning = .FALSE.
assigned_lookup = .FALSE.
j = 0

DO i = 1, lookup_dimension2
  k = 1

  ! If this field is an additional level for a previously processed field,
  ! do nothing
  IF (assigned_lookup(i)) THEN
    CYCLE
  END IF

  ! Unused lookup entries seem to have -99. These are lookup items which have
  ! space allocated but aren't being used at the moment. (Typically a 
  ! fieldsfile allows 4096 2D fields. If only 1000 are present the remaining
  ! 3096 have -99 set in the lookup). Ignore them and go to the next lookup.
  ! It's probable that after the first -99 all remaining fields will be
  ! unused too, but I can't guarantee it.
  IF (lookup(1, i) < 1) THEN
    CYCLE
  END IF

  ! Sanity check, we don't support old file versions
  ! If we read a LBREL=2 file, we note that using reset_seconds and 
  ! modify the times in the lookup later accordingly (assuming that 
  ! fields 6 and 12 are zero), so we can write an LBREL=3 file.
  ! We also only print this warning once to avoid spamming the user.
  reset_seconds = .FALSE.
  IF ( (lookup(22, i)) /= 3) THEN
    IF ( .NOT. printed_lbrel_warning) THEN
      icode = -2
      IF (lookup(22, i) == 2) THEN
        reason = ' (day num in lookup)'
      ELSE
        reason = ' (meaning unknown)'
      END IF
      WRITE(cmessage,'(A,I1, A)') newline //                                  &
        'Release header version LBREL should be 3 (seconds in lookup),' //    &
        newline // 'but is currently ', lookup(22,i), TRIM(reason)
      CALL ereport(routinename, icode, cmessage)
      printed_lbrel_warning = .TRUE.
    END IF
    reset_seconds = .TRUE.
  END IF

  j = j + 1

  ! Allocate fixed-length temporary data structures
  ALLOCATE(local_data_locations(j)%data_start(lookup_dimension2))
  ALLOCATE(local_data_locations(j)%data_length(lookup_dimension2))
  ALLOCATE(local_data_locations(j)%disk_length(lookup_dimension2))
  ALLOCATE(local_data_locations(j)%lookup_position(lookup_dimension2))
  ALLOCATE(local_levels(lookup_dimension2))

  local_data_locations(j)%data_start = imdi
  local_data_locations(j)%data_length = imdi

  ! Copy the REAL part of the lookup into REAL
  DO p = 1, 19
    real_part_of_lookup(p) = TRANSFER(lookup(p+45, i), real_part_of_lookup(1))
  END DO
    
  ! We can now fill in the horizontal grid information for this field
  CALL local_field(j)%set_grid_code(lookup(16, i))
  CALL local_field(j)%set_hemisphere_indicator(lookup(17, i))
  CALL local_field(j)%set_halo_code(query_stashmaster(lookup(42, i),          &
                                    halo_type_code))
  CALL local_field(j)%set_horiz_grid_type(query_stashmaster(lookup(42, i),    &
                                          horiz_grid_type))

  IF (this%variable_resolution) THEN
    ! This assumes Arakawa-C grid staggering... 

    IF (local_field(j)%get_horiz_grid_type() == u_points) THEN
      CALL local_field(j)%grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = lookup(19,i),                                 &
                ny     = lookup(18,i),                                 &
                list_x = column_dep_constants(:,2),                    &
                list_y = row_dep_constants(:,1) )
    ELSE IF (local_field(j)%get_horiz_grid_type() == v_points) THEN
      CALL local_field(j)%grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = lookup(19,i),                                 &
                ny     = lookup(18,i),                                 &
                list_x = column_dep_constants(:,1),                    &
                list_y = row_dep_constants(:,2) )
    ELSE
      CALL local_field(j)%grid%horiz_grid%set_defined_horizontal_grid( &
                nx     = lookup(19,i),                                 &
                ny     = lookup(18,i),                                 &
                list_x = column_dep_constants(:,1),                    &
                list_y = row_dep_constants(:,1) )
    END IF

  ELSE
    CALL local_field(j)%grid%horiz_grid%set_regular_horizontal_grid( &
        nx     = lookup(19,i),                                      &
        ny     = lookup(18,i),                                      &
        startx = real_part_of_lookup(16) + real_part_of_lookup(17), &
        starty = real_part_of_lookup(14) + real_part_of_lookup(15), &
        dx     = real_part_of_lookup(17),                           &
        dy     = real_part_of_lookup(15) )
  END IF

  ! Finally, define the pole
  CALL local_field(j)%grid%horiz_grid%set_pole(real_part_of_lookup(11),   &
                                               real_part_of_lookup(12))


  ! Start filling the metadata
  local_field(j)%validity_time(1:6) = lookup(1:6, i)
  local_field(j)%data_time(1:6)     = lookup(7:12, i)
  IF (reset_seconds) THEN 
    ! If the input file is LBREL=2, set seconds to zero
    local_field(j)%validity_time(6) = 0
    local_field(j)%data_time(6) = 0
  END IF
  local_field(j)%calendar        = calendar       ! Copied from fixed header
  local_field(j)%time_indicator  = lookup(13, i)
  local_field(j)%fctime          = lookup(14, i)
  local_field(j)%packing_method  = lookup(lbpack, i)
  local_field(j)%processing_code = lookup(25, i)
  local_field(j)%user_reference  = lookup(28, i)
  local_field(j)%run_ident       = run_ident      ! Copied from fixed header
  local_field(j)%expt_ident      = expt_ident     ! Copied from fixed header
  ! Sets lookup 38, model identifier based on version of UM code
  CALL local_field(j)%update_um_vn() 
  local_field(j)%data_type        = lookup(39, i)
  local_field(j)%quantity_ident   = lookup(42, i)
  local_field(j)%datum_constant   = real_part_of_lookup(5)
  local_field(j)%packing_accuracy = real_part_of_lookup(6)
  local_field(j)%mdi              = real_part_of_lookup(18)

  ! Information about levels and data locations
  local_levels(k) = lookup(33, i)
  local_data_locations(j)%data_start(k)  = lookup(29, i)
  local_data_locations(j)%data_length(k) = lookup(15, i)
  local_data_locations(j)%disk_length(k) = lookup(30, i)
  local_data_locations(j)%lookup_position(k) = i

  ! Remap 9999 (surface) to level 0
  IF (local_levels(k) == 9999) local_levels(k) = 0
  IF (local_levels(k) == 8888) THEN
    local_levels(k) = lookup(43, i)
    local_field(j)%meaning_of_k_dimension = 1    ! Pseudo-levels or single 
  ELSE
    local_field(j)%meaning_of_k_dimension = 0    ! Real levels
  END IF

  ! Note we've now processed this lookup
  assigned_lookup(i) = .TRUE.

  ! Look through the rest of the lookup to find additional levels of 
  ! this field.
  next_lookup:    DO m = i+1, lookup_dimension2

          !ignore unused lookup
    IF (lookup(1, m) < 1) THEN
      CYCLE next_lookup
    END IF 
        
    ! Check lookups match
    DO l = 1, SIZE(lookup_items_in_common)
      IF (lookup(lookup_items_in_common(l), i) /=                             &
                                    lookup(lookup_items_in_common(l), m)) THEN
        CYCLE next_lookup
      END IF
    END DO
        
    ! At this point, we're looking at another level for this field
    k = k + 1

    ! Store the information about this level and carry on
    local_data_locations(j)%lookup_position(k) = m
    local_data_locations(j)%data_start(k)  = lookup(29, m)
    local_data_locations(j)%disk_length(k) = lookup(30, m)
    local_data_locations(j)%data_length(k) = lookup(15, m)
    local_levels(k) = lookup(33, m)

    ! Remap 9999 (surface) to level 0
    IF (local_levels(k) == 9999) local_levels(k) = 0
    ! Remap 8888 (single level) to pseudo-level
    IF (local_levels(k) == 8888) local_levels(k) = lookup(43, m)

    assigned_lookup(m) = .TRUE.

  END DO next_lookup
    
  ! We should now have found every level; sort them into ascending order
  CALL this%sort_fields(k, local_levels, local_data_locations(j))

  ! Now sort out the vertical structure
  ! To set up heights as well as levels we probably need a STASHmaster 
  ! lookup to see what levels map to which hybrid height
  ! But we need heights for interpolation, so we need to know if we're on
  ! rho or theta levels for each quantity.
  local_field(j)%grid%vert_grid%coordinate_type = lookup(26, i)
  local_field(j)%grid%vert_grid%level_type = query_stashmaster(lookup(42, i), &
                                                               vert_level_type)

  IF (local_field(j)%grid%vert_grid%coordinate_type ==                        &
                                            vertical_coord_hybrid_height) THEN
         
      ! Levels need to be specified properly by array index
    IF (local_field(j)%grid%vert_grid%level_type == rho_levels) THEN
      CALL local_field(j)%grid%vert_grid%set_defined_vertical_grid(           &
                    nz                           = k,                         &
                    num_model_levels             = this%num_model_levels,     &
                    levels                       = local_levels,              &
                    height_at_top_theta_level    = height_at_top_theta_level, &
                    first_rho_of_constant_height =                            &
                                                first_rho_of_constant_height, &
                    list_eta_theta               =                            &
                            level_dep_constants(1:this%num_model_levels+1,1), &
                    list_eta_rho                 =                            &
                               level_dep_constants(1:this%num_model_levels,2))
    ELSE IF (local_field(j)%grid%vert_grid%level_type == theta_levels) THEN
      CALL local_field(j)%grid%vert_grid%set_defined_vertical_grid(           &
                    nz                           = k,                         &
                    num_model_levels             = this%num_model_levels,     &
                    levels                       = local_levels,              &
                    height_at_top_theta_level    = height_at_top_theta_level, &
                    first_rho_of_constant_height =                            &
                                                first_rho_of_constant_height, &
                    list_eta_theta               =                            &
                            level_dep_constants(1:this%num_model_levels+1,1), &
                    list_eta_rho                 =                            &
                               level_dep_constants(1:this%num_model_levels,2))
    ELSE 
      CALL local_field(j)%grid%vert_grid%set_defined_vertical_grid(           &
                    k, levels=local_levels) ! add list_z here
    END IF         
  ELSE
    CALL local_field(j)%grid%vert_grid%set_defined_vertical_grid(             &
                    k, levels=local_levels)
  END IF
  ! Add first and last level information to the field vertical grid
  local_field(j)%grid%vert_grid%first_level_code = query_stashmaster(         &
       lookup(42, i), first_level_code)
  CALL local_field(j)%grid%vert_grid%calc_first_level(this%grid_staggering)
  local_field(j)%grid%vert_grid%last_level_code = query_stashmaster(          &
                                               lookup(42, i), last_level_code)
  CALL local_field(j)%grid%vert_grid%calc_last_level(this%grid_staggering)
  DEALLOCATE(local_levels)    

  this%num_2d_fields_in_file = i
END DO    


! We now have the data we need to create the field objects so create them.
this%num_fields = j
  
ALLOCATE(this%fields(this%num_fields))
ALLOCATE(this%data_locations(this%num_fields))
  
  
! Should this use move_alloc?
this%fields(1:this%num_fields) = local_field(1:this%num_fields)
this%data_locations(1:this%num_fields) = local_data_locations(1:this%num_fields)

! Tidy up memory, deallocating arrays which are no longer required
DO i = 1, this%num_fields
  DEALLOCATE(local_data_locations(i)%lookup_position)
  DEALLOCATE(local_data_locations(i)%disk_length)
  DEALLOCATE(local_data_locations(i)%data_length)
  DEALLOCATE(local_data_locations(i)%data_start)
END DO

DEALLOCATE(local_data_locations)
DEALLOCATE(local_field)
DEALLOCATE(assigned_lookup)
DEALLOCATE(lookup)

IF (ALLOCATED(column_dep_constants)) DEALLOCATE(column_dep_constants)
IF (ALLOCATED(row_dep_constants)) DEALLOCATE(row_dep_constants)
DEALLOCATE(level_dep_constants)
DEALLOCATE(real_constants)
DEALLOCATE(integer_constants)


WRITE(umMessage,'(A,I8,A,I8,A,A,A,I3,A)') '[INFO] read_header: Created ',   &
   this%num_fields, ' 3D fields from ',this%num_2d_fields_in_file,          &
   ' 2D fields in file ', TRIM(this%filename) ,' (',this%unit_num,')'
CALL umPrint(umMessage)

CALL timer( 'ff:read_header', 6)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_header

!-------------------------------------------------------------------------------

SUBROUTINE read_field(this, ifield)

USE stashmaster_constants_mod, ONLY: data_type_real, data_type_integer,     &
                                     data_type_logical
USE wgdos_packing_mod, ONLY: wgdos_expand_field

IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield
  
INTEGER :: i, j, k, m
INTEGER :: nx, ny, nz, num_words, datastart, packing_type, field_size
  
INTEGER, ALLOCATABLE :: idata(:)
INTEGER, ALLOCATABLE :: int_data(:)
! Work array for reading 32 bit data
REAL(KIND=real32), ALLOCATABLE :: work_array32(:)
INTEGER :: len_io
REAL :: err_io = 0.0
INTEGER :: icode = 0
LOGICAL :: pack_mask, land_mask ! Do we have a mask? Is it a land mask?
INTEGER :: idum ! Dummy integer
CHARACTER(LEN=errormessagelength) :: cmessage     ! Error message
INTEGER, PARAMETER :: model_code_atmos = 1
CHARACTER(LEN=*), PARAMETER :: routinename = 'READ_FIELD'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'ff:read_field', 5)

! Check the header has been read. If not, read it.
IF (.NOT. ALLOCATED(this%fields)) THEN
  icode = -10
  cmessage = 'Attempted to read field before reading header. ' //       &
             'Reading header first.'
  CALL ereport(routinename, icode, cmessage)
  CALL this%read_header()
END IF
  
packing_type = MOD(this%fields(ifield)%packing_method, 10)
nx = this%fields(ifield)%get_num_cols()
ny = this%fields(ifield)%get_num_rows()
nz = this%fields(ifield)%get_num_levels()

ALLOCATE(this%fields(ifield)%rdata(nx, ny, nz))

DO i = 1, nz
  num_words = this%data_locations(ifield)%data_length(i)
  IF (packing_type == 2) THEN    ! Cray 32-bit packing
    num_words = (num_words + 1) / 2
  END IF

  datastart = this%data_locations(ifield)%data_start(i)
  CALL setpos(this%unit_num, datastart)

  ! Check data type
  IF (this%fields(ifield)%data_type == data_type_real) THEN
    ! real data

        ! Check if this is packed data
    IF (this%fields(ifield)%packing_method > PC_No_Packing) THEN

      field_size = nx * ny

      pack_mask = MOD(this%fields(ifield)%packing_method / 10, 10) ==          &
                                                         PC_BitMask_CompressType
      land_mask = MOD(this%fields(ifield)%packing_method / 100, 10) ==         &
                                                         PC_LandMask_Compression

      IF (land_mask .AND. pack_mask) THEN
        ! We can't deal with land/sea-packed data at this level. This sort of
        ! packing breaks the assumption of independent fields.
        icode = 50
        cmessage = 'Unable to read land-packed data.'
        CALL ereport(routinename, icode, cmessage)
      END IF

      IF (packing_type == PC_Cray32_Packing) THEN    ! Cray 32-bit packing
        ALLOCATE(work_array32(2*num_words))
      ELSE IF (this%fields(ifield)%packing_method == PC_WGDOS_Packing) THEN
        ! Allocate array to contain unpacked data
        ALLOCATE(int_data(num_words))
      END IF

    END IF

    ! Read in data from disk
    CALL timer( 'ff:read_field_buffin', 5)
    IF (packing_type == PC_Cray32_Packing) THEN    ! Cray 32-bit packing
       CALL buffin(this%unit_num, work_array32, 2*num_words, len_io, err_io)
    ELSE IF (this%fields(ifield)%packing_method == PC_WGDOS_Packing) THEN
       CALL buffin(this%unit_num, int_data, num_words, len_io, err_io)
    ELSE
       CALL buffin(this%unit_num, this%fields(ifield)%rdata(:,:,i), num_words, &
                   len_io, err_io)
    END IF
    CALL timer( 'ff:read_field_buffin', 6)

        ! Check if this is packed data
    IF (this%fields(ifield)%packing_method > PC_No_Packing) THEN

      num_words = this%data_locations(ifield)%data_length(i) ! Reset this

      ! Check packing method
      IF (this%fields(ifield)%packing_method == PC_Cray32_Packing) THEN
        ! 32-bit packing for dumps

        ! Copy unpacked data back into master array
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(j,k)                                                            &
!$OMP& SHARED(this,work_array32,i,ifield,nx,ny)
        DO k = 1, ny
          DO j = 1, nx
            this%fields(ifield)%rdata(j,k,i) = work_array32((k-1)*nx+j)
          END DO
        END DO
!$OMP END PARALLEL DO

       DEALLOCATE(work_array32)

      ELSE IF (this%fields(ifield)%packing_method == PC_WGDOS_Packing) THEN
        ! WGDOS packing
        icode = 0
        CALL timer( 'ff:read_wgdos_expand_field', 5)

        CALL wgdos_expand_field(this%fields(ifield)%rdata(:,:,i), field_size,  &
                                int_data, num_words,                           &
                                idum, nx, ny, idum,                            &
                                this%fields(ifield)%mdi,                       &
                                this%fields(ifield)%quantity_ident, icode)

        CALL timer( 'ff:read_wgdos_expand_field', 6)
        IF (icode /= 0) THEN
          CALL ereport(routinename, icode, cmessage)
        END IF

        DEALLOCATE(int_data)

      ELSE
        ! An unsupported packing method
        icode = 10
        WRITE(cmessage, '(A,I8)')  'Error, this field has been packed' //      &
              ' using an unsupported method. LBPACK = ',                       &
              this%fields(ifield)%packing_method
        CALL ereport(routinename, icode, cmessage)
      END IF
    END IF    ! IF Packing

  ELSE IF (this%fields(ifield)%data_type == data_type_integer .OR.             &
           this%fields(ifield)%data_type == data_type_logical ) THEN
    ! Integer or logical data. For some reason logical data is integer...

    CALL buffin(this%unit_num, this%fields(ifield)%rdata(:,:,i), num_words,    &
                len_io, err_io)

    IF (this%fields(ifield)%packing_method > PC_No_Packing) THEN
      icode = 11
      cmessage = 'Cannot handle packed integer/logical data'
      CALL ereport(routinename, icode, cmessage)
    END IF

  END IF

END DO
CALL timer( 'ff:read_field', 6)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_field

!-------------------------------------------------------------------------------

SUBROUTINE write_header(this)
USE stashmaster_utils_mod, ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: meto8_code, pp_field_code
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER :: fixed_header(this%fixed_header_length)
INTEGER, ALLOCATABLE :: integer_constants(:)
REAL, ALLOCATABLE :: real_constants(:)
REAL, ALLOCATABLE :: level_dep_constants(:,:)
REAL, ALLOCATABLE :: row_dep_constants(:,:)
REAL, ALLOCATABLE :: column_dep_constants(:,:)
INTEGER, ALLOCATABLE :: lookup(:,:)
REAL :: real_part_of_lookup(19)
INTEGER :: start_of_data

INTEGER :: i, p, j, k
INTEGER :: file_position

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='WRITE_HEADER'
INTEGER :: icode
INTEGER :: nx, ny
REAL :: startx, starty, dx, dy, pole_lat, pole_long

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'ff:write_header', 5)
! Set the fixed, integer, real headers sizes and sizes of level,row and
! column dependent constants 
CALL this%set_header_sizes()
fixed_header = imdi
fixed_header(1) = 20
fixed_header(2) = 1
fixed_header(3) = this%vertical_grid_indicator
fixed_header(4) = this%horizontal_grid_indicator
fixed_header(5) = 3
IF (ALLOCATED(this%fields)) THEN
  fixed_header(6) = this%fields(1)%run_ident
  fixed_header(7) = this%fields(1)%expt_ident
  fixed_header(8) = this%fields(1)%calendar
END IF
fixed_header(9) = this%grid_staggering
fixed_header(12) = this%source_model_version
CALL this%process_times(fixed_header(21:26), fixed_header(28:33),             &
                        fixed_header(35:40))
fixed_header(27) = 0
fixed_header(34) = 0
fixed_header(41) = 0

! Integer constants dimensions
fixed_header(100) = this%integer_constants_start
fixed_header(101) = this%integer_constants_dimension
! Real constants dimensions
fixed_header(105) = this%real_constants_start
fixed_header(106) = this%real_constants_dimension
! Level dependent constants dimensions
fixed_header(110) = this%level_dep_constants_start
fixed_header(111) = this%level_dep_constants_dimension1
fixed_header(112) = this%level_dep_constants_dimension2
! Row dependent constants
fixed_header(115) = this%row_dep_constants_start
fixed_header(116) = this%row_dep_constants_dimension1
fixed_header(117) = this%row_dep_constants_dimension2
fixed_header(120) = this%column_dep_constants_start
fixed_header(121) = this%column_dep_constants_dimension1
fixed_header(122) = this%column_dep_constants_dimension2

fixed_header(150) = this%lookup_start
fixed_header(151) = this%len_single_lookup
fixed_header(152) = this%num_2d_fields_in_file
fixed_header(160) = this%data_start
fixed_header(161) = this%most_recent_data_position

IF (this%num_2d_fields_in_file > this%num_reserved_headers) THEN
  cmessage = 'Header too large for data, increase the number' //              &
             ' of reserved lookup headers'
  icode = 10
  CALL ereport(routinename, icode, cmessage)
END IF
  
CALL setpos(this%unit_num, this%start_of_file)
CALL buffout(this%unit_num, fixed_header, this%fixed_header_length, icode)

! Integer constants
ALLOCATE(integer_constants(this%len_integer_constants))
ALLOCATE(real_constants(this%len_real_constants))
integer_constants    = imdi
real_constants = rmdi
integer_constants(8) = this%file_theta_rho_levels%get_num_levels()
integer_constants(9) = this%file_theta_rho_levels%get_num_levels()
integer_constants(10) = this%num_soil_levels
integer_constants(11) = this%num_cloud_levels 
IF (this%number_of_tracers() > 0) THEN
  integer_constants(12) = this%num_tracer_levels
ELSE
  integer_constants(12) = 0
END IF
integer_constants(13) = this%num_boundary_layer_levels
integer_constants(15) = this%number_of_unique_stash_items()
integer_constants(17) = this%algorithm_to_generate_height_fields
integer_constants(24) = this%file_theta_rho_levels%first_rho_of_constant_height
integer_constants(25) = this%num_land_points 
integer_constants(28) = this%num_soil_moisture_levels
real_constants(16)    = this%file_theta_rho_levels%height_at_top_theta_level

integer_constants(6)  = this%p_grid%get_num_cols()
integer_constants(7)  = this%p_grid%get_num_rows()

IF (.NOT. this%variable_resolution) THEN    
  CALL this%p_grid%get_regular_horizontal_grid(nx, ny, startx, starty, dx, dy)
  real_constants(1) = dx
  real_constants(2) = dy
  IF (this%grid_staggering == new_dynamics) THEN
    real_constants(3) = starty
    real_constants(4) = startx
  ELSE IF (this%grid_staggering == endgame) THEN
    real_constants(3) = starty - 0.5 * dy ! Start lat/long is offset from
    real_constants(4) = startx - 0.5 * dx ! first P point for ENDGame
  END IF
END IF
real_constants(5) = this%pole_lat
real_constants(6) = this%pole_long

CALL setpos(this%unit_num, fixed_header(100) -1 )
CALL buffout(this%unit_num, integer_constants, this%len_integer_constants)
CALL setpos(this%unit_num, fixed_header(105) -1 )
CALL buffout(this%unit_num, real_constants, this%len_real_constants)
DEALLOCATE(real_constants)
DEALLOCATE(integer_constants)

ALLOCATE(level_dep_constants(this%file_theta_rho_levels%get_num_levels()+1,   &
                             this%len2_lev_dep_constants))
level_dep_constants = rmdi
level_dep_constants(1:this%file_theta_rho_levels%get_num_levels()+1,1) =      &
                                       this%file_theta_rho_levels%eta_theta(:)
level_dep_constants(1:this%file_theta_rho_levels%get_num_levels(),2) =        &
                                         this%file_theta_rho_levels%eta_rho(:)
level_dep_constants(1:this%file_rhcrit%get_num_levels(),3) =                  &
                                              this%file_rhcrit%level_values(:)
level_dep_constants(1:this%num_soil_levels,4) =                               &
                                      this%file_soil_thickness%level_values(:)
level_dep_constants(1:this%file_zsea_theta%get_num_levels(),5) =              &
                                          this%file_zsea_theta%level_values(:)
level_dep_constants(1:this%file_c_theta%get_num_levels(),6) =                 &
                                          this%file_c_theta%level_values(:)
level_dep_constants(1:this%file_zsea_rho%get_num_levels(),7) =                &
                                          this%file_zsea_rho%level_values(:)
level_dep_constants(1:this%file_c_rho%get_num_levels(),8) =                   &
                                          this%file_c_rho%level_values(:)
CALL setpos(this%unit_num, fixed_header(110) - 1)
CALL buffout(this%unit_num, level_dep_constants,                              &
             fixed_header(111) * fixed_header(112))
DEALLOCATE(level_dep_constants)

IF (this%variable_resolution) THEN
  ALLOCATE(row_dep_constants(fixed_header(116), fixed_header(117)))
  ALLOCATE(column_dep_constants(fixed_header(121), fixed_header(122)))
  row_dep_constants = rmdi
  column_dep_constants = rmdi

  row_dep_constants(1:this%p_grid%get_num_rows(),1) =                         &
                                                   this%p_grid%get_latitudes()
  row_dep_constants(1:this%v_grid%get_num_rows(),2) =                         &
                                                   this%v_grid%get_latitudes()
  column_dep_constants(1:this%p_grid%get_num_cols(),1) =                      &
                                                   this%p_grid%get_longitudes()
  column_dep_constants(1:this%u_grid%get_num_cols(),2) =                      &
                                                   this%u_grid%get_longitudes()

  CALL setpos(this%unit_num, fixed_header(115) -1 )
  CALL buffout(this%unit_num, row_dep_constants,                              &
                                          fixed_header(116)*fixed_header(117))
  CALL setpos(this%unit_num, fixed_header(120) -1 )
  CALL buffout(this%unit_num, column_dep_constants,                           &
                                          fixed_header(121)*fixed_header(122))
  DEALLOCATE(column_dep_constants)
  DEALLOCATE(row_dep_constants)
END IF

! Construct lookup header
ALLOCATE(lookup(this%len_single_lookup, this%num_reserved_headers))
lookup = -99

k = 1
DO i = 1, this%num_fields
  DO j = 1, this%fields(i)%get_num_levels()

    ! Constant part of lookup
    lookup(:,k) = imdi
    lookup(1:6,k) = this%fields(i)%validity_time(1:6)
    lookup(7:12,k) = this%fields(i)%data_time(1:6)
    lookup(13,k) = this%fields(i)%time_indicator
    lookup(14,k) = this%fields(i)%fctime
    lookup(16,k) = this%fields(i)%get_grid_code()
    lookup(17,k) = this%fields(i)%get_hemisphere_indicator()
    lookup(lbpack,k) = this%fields(i)%packing_method
    lookup(22,k) = this%fields(i)%calc_header_release()
    lookup(23,k) = query_stashmaster(this%fields(i)%quantity_ident,           &
                                     pp_field_code)
    lookup(25,k) = this%fields(i)%processing_code
    lookup(26,k) = this%fields(i)%grid%vert_grid%coordinate_type
    lookup(27,k) = 0
    lookup(28,k) = this%fields(i)%user_reference
    lookup(32,k) = query_stashmaster(this%fields(i)%quantity_ident,meto8_code)
    lookup(38,k) = this%fields(i)%source_model_identifier
    lookup(39,k) = this%fields(i)%data_type
    lookup(42,k) = this%fields(i)%quantity_ident

    lookup(43,k) = this%fields(i)%meaning_of_k_dimension
    lookup(45,k) = 1

    real_part_of_lookup = rmdi
    real_part_of_lookup(5) = this%fields(i)%datum_constant
    real_part_of_lookup(6) = this%fields(i)%packing_accuracy
    real_part_of_lookup(18) = this%fields(i)%mdi

    lookup(18,k) = this%fields(i)%get_num_rows()
    lookup(19,k) = this%fields(i)%get_num_cols()
    real_part_of_lookup(13) = 0
    IF (this%variable_resolution) THEN
      real_part_of_lookup(14) = rmdi
      real_part_of_lookup(15) = rmdi
      real_part_of_lookup(16) = rmdi
      real_part_of_lookup(17) = rmdi
    ELSE
      ! Set up grid ; remember start positions are zeroth, not first.
      CALL this%fields(i)%grid%horiz_grid%get_regular_horizontal_grid(nx, ny, &
                                                       startx, starty, dx, dy)
      real_part_of_lookup(14) = starty - dy
      real_part_of_lookup(15) = dy
      real_part_of_lookup(16) = startx - dx
      real_part_of_lookup(17) = dx
    END IF
    CALL this%fields(i)%get_pole(pole_lat, pole_long)
    real_part_of_lookup(11) = pole_lat
    real_part_of_lookup(12) = pole_long

    ! Copy the REAL part of the lookup into INTEGER
    DO p = 1, 19
      lookup(p+45,k) = TRANSFER(real_part_of_lookup(p), lookup(1,1))
    END DO

    ! Level dependent section

    ! Get level number and remap if appropriate
    lookup(33,k) = this%fields(i)%grid%vert_grid%levels(j)
    IF (lookup(33,k) == 0) lookup(33,k) = 9999
    IF (this%fields(i)%meaning_of_k_dimension == 1) THEN
      lookup(43,k) = lookup(33,k)
      lookup(33,k) = 8888
    END IF

    ! Disk locations - these should be set by write_fields
    lookup(29,k) = this%data_locations(i)%data_start(j)
    lookup(15,k) = this%data_locations(i)%data_length(j)
    lookup(30,k) = this%data_locations(i)%disk_length(j)
    lookup(40,k) = this%data_locations(i)%data_start(j)
    this%data_locations(i)%lookup_position(j) = k
    k = k + 1
  END DO
END DO

! Write lookup   
CALL setpos(this%unit_num, fixed_header(150) -1 )
CALL buffout(this%unit_num, lookup,                                           &
             this%len_single_lookup*this%num_reserved_headers)

DEALLOCATE(lookup)
  
CALL timer( 'ff:write_header', 6)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE write_header

!-------------------------------------------------------------------------------

SUBROUTINE write_field(this, ifield)
USE io_configuration_mod, ONLY: &
  io_field_padding
USE wgdos_packing_mod, ONLY: wgdos_compress_field
USE cray32_packing_mod, ONLY: ff_write_packed_field
USE stashmaster_utils_mod, ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: packing_acc
USE umPrintMgr, ONLY: newline
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield
INTEGER :: position_to_write, j, len_data, len_disk
INTEGER :: l,m,n
REAL :: stashmaster_packing_accuracy
  
INTEGER, ALLOCATABLE :: int_data(:)
INTEGER :: nx
INTEGER :: num_words, temp_packing, num_64_bit_words
  
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='WRITE_FIELD'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(this%fields)) THEN
  icode = 51
  WRITE(cmessage, '(A,I0,A)') 'Attempted to write field ', ifield,  &
                              ' but no fields exist.'
  CALL ereport(routinename, icode, cmessage)
END IF

IF (.NOT. ALLOCATED(this%fields(ifield)%rdata)) THEN
  icode = 52
  WRITE(cmessage, '(A,I0,A)') 'Attempted to write field ', ifield,  &
                              ' but field contains no data.'
  CALL ereport(routinename, icode, cmessage)
END IF

CALL timer( 'ff:write_field', 5)
IF (this%most_recent_data_position > 0) THEN 
  position_to_write = this%most_recent_data_position

  !The else if and else conditions are both used if no data has been written to
  !the file; therefore we clear any data locations - once we've written to a 
  !file we can no longer read from it.
ELSE IF (this%data_start > 0) THEN
  position_to_write = this%data_start - 1
  CALL this%clear_data_locations()
ELSE
  icode = 55
  cmessage = "Unable to write field to file. The value of data_start " //     &
             "in the file object is <= 0"
  CALL ereport(routinename, icode, cmessage)
END IF

stashmaster_packing_accuracy =  &
    query_stashmaster( this%fields(ifield)%quantity_ident, packing_acc)

IF (this%fields(ifield)%packing_method == PC_No_Packing) THEN
  ! Unpacked data has a packing accuracy of -99.0
  this%fields(ifield)%packing_accuracy = packing_unpacked_accuracy
ELSE 
  ! Warn if the packing accuracy doesn't match the STASHmaster and reset
  IF (this%fields(ifield)%packing_accuracy /= stashmaster_packing_accuracy) THEN
    WRITE(cmessage , '(A,I0,A,A,F16.8,A,F16.8,A)')                             &
       'Resetting packing accuracy of field (STASH code ',                     &
       this%fields(ifield)%quantity_ident , ')', newline //                    &
       ' from ',this%fields(ifield)%packing_accuracy,                          &
       ' to ', stashmaster_packing_accuracy, ' to match STASHmaster value'
    icode = -550
    CALL ereport(routinename, icode, cmessage)
    this%fields(ifield)%packing_accuracy = stashmaster_packing_accuracy
  END IF 
END IF

IF (INT(this%fields(ifield)%packing_accuracy) <=                               &
                                            INT(packing_unpacked_accuracy)) THEN
  ! Switch off packing if packing accuracy is -99
  temp_packing = this%fields(ifield)%packing_method / 10
  this%fields(ifield)%packing_method = temp_packing * 10
END IF


DO j = 1, this%fields(ifield)%get_num_levels()

  len_data = SIZE(this%fields(ifield)%rdata(:,:,j),1) *                    &
             SIZE(this%fields(ifield)%rdata(:,:,j),2)

  IF (MOD(this%fields(ifield)%packing_method, 10) /= PC_No_Packing) THEN
    ! Force WGDOS packing as it's the only type implemented
    temp_packing = this%fields(ifield)%packing_method / 10
    this%fields(ifield)%packing_method = temp_packing * 10 + PC_WGDOS_Packing
  END IF

  IF (MOD(this%fields(ifield)%packing_method, 10) == PC_No_Packing) THEN

    len_disk = len_data
    CALL setpos(this%unit_num, position_to_write)
    CALL buffout(this%unit_num, this%fields(ifield)%rdata(:,:,j), len_disk)

    this%data_locations(ifield)%data_start(j) = position_to_write
    this%data_locations(ifield)%data_length(j) = len_data
    this%data_locations(ifield)%disk_length(j) = len_disk
    position_to_write = position_to_write + len_disk
  ELSE IF (MOD(this%fields(ifield)%packing_method, 10) == PC_Cray32_Packing)   &
  THEN
      ! Cray 32-bit packing


    len_data = SIZE(this%fields(ifield)%rdata(:,:,j),1) *                    &
               SIZE(this%fields(ifield)%rdata(:,:,j),2)

    ! Number of 64-bit words required to store the 32-bit data
    len_disk = (len_data+1)/2

    ! round up to sector boundary
    len_disk=((len_disk+io_field_padding-1)/                                 &
                   io_field_padding)*io_field_padding

    CALL ff_write_packed_field(this%unit_num, this%fields(ifield)%rdata,     &
                    SIZE(this%fields(ifield)%rdata(:,:,j),1),                &
                    SIZE(this%fields(ifield)%rdata(:,:,j),2),                &
                    position_to_write, len_disk)

    this%data_locations(ifield)%data_start(j) = position_to_write
    this%data_locations(ifield)%data_length(j) = len_data
    this%data_locations(ifield)%disk_length(j) = len_disk
    position_to_write = position_to_write + len_disk

  ELSE 
    ! WGDOS packing
    IF (MOD(this%fields(ifield)%packing_method, 10) /= PC_WGDOS_Packing) THEN
      ! An unsupported packing method - default to WGDOS
      icode = -10
      WRITE(cmessage, '(A,I2,A)')  'Warning, this field needs to be packed'  &
           // ' using an unsupported method. LBPACK = ',                     &
           this%fields(ifield)%packing_method, ', attempting to WGDOS pack'  &
           // ' instead'
      CALL ereport(routinename, icode, cmessage)
      temp_packing = this%fields(ifield)%packing_method / 10
      this%fields(ifield)%packing_method = temp_packing * 10 + 1
    END IF

    nx = this%fields(ifield)%get_num_cols()

    ALLOCATE(int_data(len_data))
    icode = 0
    cmessage = 'No problems reported by wgdos_compress_field'
    CALL timer( 'ff:write_wgdos_compress_field', 5)
    CALL wgdos_compress_field(                                                &
                              this%fields(ifield)%rdata(:,:,j),               &
                              len_data,                                       &
                              int_data,                                       &
                              len_data,                                       &
                              nx, num_words,                                  &
                              INT(this%fields(ifield)%packing_accuracy),      &
                              this%fields(ifield)%mdi,                        &
                              this%fields(ifield)%quantity_ident, icode) 

    CALL timer( 'ff:write_wgdos_compress_field', 6)
    IF (icode /= 0) THEN
      CALL ereport(routinename, icode, cmessage)
    END IF

    num_64_bit_words = CEILING(REAL(num_words)/2.0)
    CALL setpos(this%unit_num, position_to_write)
    CALL buffout(this%unit_num, int_data(1:num_64_bit_words),               &
                  num_64_bit_words)
    this%data_locations(ifield)%data_start(j) = position_to_write
    this%data_locations(ifield)%data_length(j) = num_64_bit_words
    this%data_locations(ifield)%disk_length(j) = &
         ((((num_words+1)/2)+io_field_padding-1)/io_field_padding)*           &
                                                              io_field_padding
    position_to_write = position_to_write +                                   &
                                    this%data_locations(ifield)%disk_length(j)

    DEALLOCATE(int_data)

  END IF

END DO

this%most_recent_data_position = position_to_write 

CALL timer( 'ff:write_field', 6)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE write_field

!-------------------------------------------------------------------------------

SUBROUTINE allow_read_only(this, read_only)
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
LOGICAL, INTENT(IN) :: read_only
CHARACTER(LEN=*), PARAMETER :: routinename = 'ALLOW_READ_ONLY'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (read_only) THEN
  this%read_write_status = ioOpenReadOnly
ELSE
  this%read_write_status = ioOpenReadWrite
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE allow_read_only
!-------------------------------------------------------------------------------

SUBROUTINE open_file(this, filename)
  ! This overrides the "default" open_file in the datafile class so we can use
  ! the UM's IO system to open files "properly".

USE file_manager, ONLY: assign_file_unit
USE filenamelength_mod, ONLY: filenamelength
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
CHARACTER(LEN=filenamelength), INTENT(IN) :: filename

INTEGER :: icode,read_write_status
CHARACTER(LEN=*), PARAMETER :: routinename='OPEN_FILE'
CHARACTER(LEN=errormessagelength) :: cmessage

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

this%filename = filename
CALL assign_file_unit(this%filename, this%unit_num, handler="portio")
IF (this%read_write_status == imdi) THEN
  this%read_write_status = ioOpenReadWrite
END IF
CALL file_open(this%unit_num, this%filename, filenamelength,                  &
               read_write=this%read_write_status,                             &
               name_in_environ=ioNameProvided, error=icode)
IF (icode /= 0) THEN
  icode = ABS(icode)
  cmessage = 'Error attempting to open file '// this%filename
  CALL ereport(routinename, icode, cmessage)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE open_file

!-------------------------------------------------------------------------------

SUBROUTINE close_file(this)
  ! This overrides the "default" close_file in the datafile class so we can use
  ! the UM's IO system to close files "properly". 

USE filenamelength_mod, ONLY: filenamelength
USE file_manager, ONLY: release_file_unit
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'CLOSE_FILE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL file_close(this%unit_num, this%filename)
CALL release_file_unit(this%unit_num, handler="portio")
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE close_file


!-------------------------------------------------------------------------------

SUBROUTINE sort_fields(this, n, levels, data_addresses)

IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(INOUT) :: levels(n)
TYPE(data_location_type), INTENT(INOUT) :: data_addresses

LOGICAL :: swapped = .FALSE.
INTEGER :: i, num_left_to_sort, tmp, next_num_to_sort
CHARACTER(LEN=*), PARAMETER :: routinename = 'SORT_FIELDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_left_to_sort = n 
  
DO WHILE(num_left_to_sort > 0)
  next_num_to_sort = 0
  DO i = 2, num_left_to_sort - 1
    IF (levels(i-1) > levels(i)) THEN
      tmp = levels(i-1)
      levels(i-1) = levels(i)
      levels(i) = tmp
        
      tmp = data_addresses%lookup_position(i-1)
      data_addresses%lookup_position(i-1) = data_addresses%lookup_position(i)
      data_addresses%lookup_position(i) = tmp

      tmp = data_addresses%data_start(i-1)
      data_addresses%data_start(i-1) = data_addresses%data_start(i)
      data_addresses%data_start(i) = tmp


      tmp = data_addresses%disk_length(i-1)
      data_addresses%disk_length(i-1) = data_addresses%disk_length(i)
      data_addresses%disk_length(i) = tmp

      tmp = data_addresses%data_length(i-1)
      data_addresses%lookup_position(i-1) = data_addresses%data_length(i)
      data_addresses%data_length(i) = tmp
        
      next_num_to_sort = i
    END IF
  END DO
  num_left_to_sort = next_num_to_sort
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE sort_fields  

!-------------------------------------------------------------------------------

INTEGER FUNCTION add_field(this, new_field)
USE field_mod, ONLY: field_type, ASSIGNMENT(=)
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE
  
CLASS(fieldsfile_type), INTENT(INOUT) :: this
TYPE(field_type), INTENT(IN) :: new_field

TYPE(field_type), ALLOCATABLE :: tmp_fields(:)
TYPE(data_location_type), ALLOCATABLE :: tmp_data_locations(:)
INTEGER :: new_num_levels
INTEGER :: i
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename = 'ADD_FIELD'
CHARACTER(LEN=errormessagelength) :: cmessage
LOGICAL :: new_file  

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Abort if the new field contains no data
IF (.NOT. ALLOCATED(new_field%rdata) .AND.  &
    .NOT. ALLOCATED(new_field%lbc_rdata)) THEN
  icode= 79
  cmessage = 'Attempted to add field which contains no data to file'
  CALL ereport(routinename, icode, cmessage)
END IF

IF (ALLOCATED(this%fields)) THEN
  ALLOCATE(tmp_fields(this%num_fields))
  ALLOCATE(tmp_data_locations(this%num_fields))
  tmp_fields(1:this%num_fields) = this%fields(1:this%num_fields)
  tmp_data_locations(1:this%num_fields) = this%data_locations(1:this%num_fields)
  DEALLOCATE(this%fields)
  DEALLOCATE(this%data_locations)

  this%num_fields = this%num_fields + 1

  ALLOCATE(this%fields(this%num_fields))
  ALLOCATE(this%data_locations(this%num_fields))
  this%fields(1:this%num_fields-1) = tmp_fields(1:this%num_fields-1)  
  this%data_locations(1:this%num_fields-1) =                                  &
                                     tmp_data_locations(1:this%num_fields-1)  
  DEALLOCATE(tmp_data_locations)
  DEALLOCATE(tmp_fields)

  this%fields(this%num_fields) = new_field

  new_num_levels = new_field%get_num_levels()

  ALLOCATE(this%data_locations(this%num_fields)%lookup_position(new_num_levels))
  ALLOCATE(this%data_locations(this%num_fields)%data_start(new_num_levels))
  ALLOCATE(this%data_locations(this%num_fields)%data_length(new_num_levels))
  ALLOCATE(this%data_locations(this%num_fields)%disk_length(new_num_levels))

  ! These will be set on write (or reset, for existing fields)
  DO i = 1, new_num_levels
    this%data_locations(this%num_fields)%lookup_position(i) =                &
                                            this%num_2d_fields_in_file + 1
    this%data_locations(this%num_fields)%data_start(i) = imdi
    this%data_locations(this%num_fields)%data_length(i) = imdi
    this%data_locations(this%num_fields)%disk_length(i) = imdi 
    this%num_2d_fields_in_file = this%num_2d_fields_in_file + 1
  END DO
ELSE
  ! Empty file
  this%num_fields = 1
  ALLOCATE(this%fields(1))
  ALLOCATE(this%data_locations(1))
  ALLOCATE(this%data_locations(1)%lookup_position(new_field%get_num_levels()))
  ALLOCATE(this%data_locations(1)%data_start(new_field%get_num_levels()))
  ALLOCATE(this%data_locations(1)%data_length(new_field%get_num_levels()))
  ALLOCATE(this%data_locations(1)%disk_length(new_field%get_num_levels()))
  this%fields(1) = new_field  
  this%data_locations(1)%data_start = imdi 
  this%data_locations(1)%data_length = imdi 
  this%data_locations(1)%disk_length = imdi 
  DO i = 1, new_field%get_num_levels()
    this%data_locations(1)%lookup_position(i) = i
  END DO
  this%num_2d_fields_in_file = new_field%get_num_levels()
END IF
add_field = this%num_fields
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION add_field

!-------------------------------------------------------------------------------

TYPE(field_type) FUNCTION field(this, ifield)
IMPLICIT NONE
CLASS(fieldsfile_type) :: this
INTEGER :: ifield
CHARACTER(LEN=*), PARAMETER :: routinename = 'FIELD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
field = this%fields(ifield)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION field

!-------------------------------------------------------------------------------

SUBROUTINE clear_data_locations(this)
IMPLICIT NONE
CLASS(fieldsfile_type) :: this
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'CLEAR_DATA_LOCATIONS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (.NOT. ALLOCATED(this%fields)) THEN
  ! If the file doesn't contain any fields there are no data locations to clear.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
  
DO i = 1, this%num_fields
  this%data_locations(i)%disk_length = imdi
  this%data_locations(i)%data_length = imdi
  this%data_locations(i)%data_start = imdi
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE clear_data_locations

!-------------------------------------------------------------------------------

INTEGER FUNCTION number_of_dust_bins_in_file(this)
USE stashmaster_utils_mod, ONLY: number_of_dust_fields

IMPLICIT NONE
CLASS(fieldsfile_type) :: this

INTEGER :: stash_codes(this%num_fields)
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'NUMBER_OF_DUST_BINS_IN_FILE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, this%num_fields
  stash_codes(i) = this%fields(i)%quantity_ident
END DO

number_of_dust_bins_in_file = number_of_dust_fields(stash_codes)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION number_of_dust_bins_in_file

!-------------------------------------------------------------------------------

SUBROUTINE check_source_cyclic(this)
IMPLICIT NONE
CLASS(fieldsfile_type) :: this
INTEGER                :: icode
CHARACTER(LEN=errormessagelength)     :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'CHECK_SOURCE_CYCLIC'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Please UMDP F3, fixed header 4 for details of the horizontal grid
! indicators
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
SELECT CASE (this%horizontal_grid_indicator)
  ! Global or wrapping LAM, add 100 for rotated grids
CASE (0:2, 4, 100:102, 104)
  this%l_source_cyclic = .TRUE.
  ! LAM (No wrap)
CASE (3, 103)
  this%l_source_cyclic = .FALSE.
CASE DEFAULT
  icode = 45
  WRITE(cmessage, '(A,I5)')  'Invalid horizontal grid type indicator = ', &
       this%horizontal_grid_indicator
  CALL ereport(routinename, icode, cmessage)
END SELECT
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_source_cyclic
!-------------------------------------------------------------------------------

SUBROUTINE process_times(this, first_validity_time, last_validity_time,       &
                         validity_interval)
USE time_utils_mod, ONLY: timedate_to_seconds
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(IN) :: this
INTEGER :: first_validity_time(6), last_validity_time(6), validity_interval(6)
  
INTEGER :: i, first, last, first_secs, last_secs, i_secs, i1_secs, dsec
INTEGER :: nday, nhr, nmin, nsec
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS_TIMES'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
first_validity_time = 0
last_validity_time = 0
validity_interval = 0  
  
first = 1
last = 1

IF (.NOT. ALLOCATED(this%fields)) THEN
  ! If the file doesn't contain any fields we can't say anything about the 
  ! validity times until a field is in the file.
  first_validity_time = imdi
  last_validity_time = imdi
  validity_interval = imdi
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

first_secs = timedate_to_seconds(this%fields(first)%validity_time)
last_secs = timedate_to_seconds(this%fields(last)%validity_time)
i1_secs = first_secs ! Store first field time
dsec = 0.0

IF (this%num_fields > 1) THEN  
  DO i = 2, this%num_fields
    i_secs = timedate_to_seconds(this%fields(i)%validity_time)
    IF (i_secs < first_secs) THEN
      first_secs = i_secs
      first = i
    END IF

    IF (i_secs > last_secs) THEN
      last_secs = i_secs
      last = i
    END IF

    IF (dsec == 0.0 .AND. i_secs - first_secs > 0.0) THEN
      dsec = i_secs - first_secs
    END IF
  END DO
END IF

first_validity_time = this%fields(first)%validity_time
last_validity_time = this%fields(last)%validity_time

! Calculate interval
nday = dsec / isec_per_day
dsec = dsec - (nday * isec_per_day)
  
nhr = dsec / isec_per_hour
dsec = dsec - (nhr * isec_per_hour)
  
nmin = dsec / isec_per_min
dsec = dsec - (nmin*isec_per_min)
  
nsec = dsec
  
validity_interval(3) = nday
validity_interval(4) = nhr
validity_interval(5) = nmin
validity_interval(6) = nsec
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE process_times

!-------------------------------------------------------------------------------

SUBROUTINE set_header_sizes(this)
USE io_configuration_mod, ONLY: io_data_alignment
! Once a file contains the P, U and V grids objects it is possible to determine
! the size of the fixed header, integer header, real header, level, row and 
! column dependent constants and lookup headers
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER :: file_position = 0
INTEGER :: n_v_rows, n_p_cols, n_u_cols, n_p_rows
! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER  :: routinename = 'SET_HEADER_SIZES'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Get number rows and columns for grids
n_p_rows = this%p_grid%get_num_rows()
n_v_rows = this%v_grid%get_num_rows()
n_p_cols = this%p_grid%get_num_cols()
n_u_cols = this%u_grid%get_num_cols()

! Integer constants
file_position = this%fixed_header_length + 1
this%integer_constants_start     = file_position
this%integer_constants_dimension = this%len_integer_constants
file_position = file_position + this%integer_constants_dimension
! Real constants
this%real_constants_start        = file_position
this%real_constants_dimension    = this%len_real_constants
file_position = file_position + this%real_constants_dimension
! Level dependent constants
this%level_dep_constants_start  = file_position
this%level_dep_constants_dimension1 = this % file_theta_rho_levels %          &
                                                          num_model_levels + 1
this%level_dep_constants_dimension2 = this%len2_lev_dep_constants
file_position = file_position + (this%level_dep_constants_dimension1 *        &
                                          this%level_dep_constants_dimension2)

IF (this%variable_resolution) THEN
  ! Row dependent constants
  this%row_dep_constants_start = file_position
  this%row_dep_constants_dimension1 = MAX(n_v_rows, n_p_rows)
  this%row_dep_constants_dimension2 = this%len2_row_dep_constants
  file_position = file_position + (this%row_dep_constants_dimension1 *        &
                                            this%row_dep_constants_dimension2)
  ! Column dependent constants
  this%column_dep_constants_start = file_position
  this%column_dep_constants_dimension1 = MAX(n_p_cols, n_u_cols)
  this%column_dep_constants_dimension2 = this%len2_col_dep_constants
  file_position = file_position + (this%column_dep_constants_dimension1 *     &
                                         this%column_dep_constants_dimension2)
ELSE ! If not variable resolution then set to imdi
  this%row_dep_constants_start = imdi
  this%row_dep_constants_dimension1 = imdi
  this%row_dep_constants_dimension2 = imdi
  this%column_dep_constants_start = imdi
  this%column_dep_constants_dimension1 = imdi
  this%column_dep_constants_dimension2 = imdi
END IF
! Lookup headers
this%lookup_start = file_position
! Calculate the where the data should start.  Apply correction in order to
! round up to a disk sector boundary.
this%data_start = ((this%lookup_start + (this%len_single_lookup *             &
                  this%num_reserved_headers) -1 + io_data_alignment - 1) /    &
                   io_data_alignment)*io_data_alignment + 1
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_header_sizes

!-------------------------------------------------------------------------------

INTEGER FUNCTION number_of_unique_stash_items(this)
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(IN) :: this
INTEGER :: unique_stash_codes(9999) = imdi
INTEGER :: counter = 1
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'NUMBER_OF_UNIQUE_STASH_ITEMS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, SIZE(this%fields)
  IF (ANY(unique_stash_codes == this%fields(i)%quantity_ident)) THEN
    CYCLE
  ELSE
    unique_stash_codes(counter) = this%fields(i)%quantity_ident
    counter = counter + 1
  END IF
END DO
number_of_unique_stash_items = counter - 1
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION number_of_unique_stash_items

!-------------------------------------------------------------------------------

INTEGER FUNCTION number_of_tracers(this)
USE um_stashcode_mod, ONLY: stashcode_tracer_sec, stashcode_ukca_sec,        &
                            stashcode_lbc_free_tracer_1, stashcode_lbc_ukca_1
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(IN) :: this
LOGICAL :: ukca_tracers(999) = .FALSE.
LOGICAL :: free_tracers(999) = .FALSE.
INTEGER :: i, loop_quantity_ident
CHARACTER(LEN=*), PARAMETER :: routinename = 'NUMBER_OF_TRACERS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, this%num_fields
  loop_quantity_ident = this%fields(i)%quantity_ident
  IF (loop_quantity_ident / 1000 == stashcode_tracer_sec) THEN
    free_tracers(MOD(this%fields(i)%quantity_ident, 1000)) = .TRUE.
  ELSE IF (loop_quantity_ident / 1000 == stashcode_ukca_sec ) THEN
    ukca_tracers(MOD(this%fields(i)%quantity_ident, 1000)) = .TRUE.
  ELSE IF (loop_quantity_ident / 1000 == stashcode_lbc_free_tracer_1          &
                                                                  / 1000) THEN
    free_tracers(MOD(loop_quantity_ident, 1000)) = .TRUE.
  ELSE IF (loop_quantity_ident / 1000 == stashcode_lbc_ukca_1 / 1000 ) THEN
    ukca_tracers(MOD(this%fields(i)%quantity_ident, 1000)) = .TRUE.
  END IF
END DO

number_of_tracers = COUNT(ukca_tracers) + COUNT(free_tracers)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION number_of_tracers

!-------------------------------------------------------------------------------

SUBROUTINE copy_common_metadata(this, source_file)
! Copy common metadata from another fieldsfile
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
CLASS(fieldsfile_type), INTENT(IN)    :: source_file
CHARACTER(LEN=*), PARAMETER :: routinename = 'COPY_COMMON_METADATA'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%source_model_version                = source_file%source_model_version
this%grid_staggering                     = source_file%grid_staggering
this%vertical_grid_indicator             = source_file%vertical_grid_indicator
this%horizontal_grid_indicator           = source_file%horizontal_grid_indicator
this%file_theta_rho_levels               = source_file%file_theta_rho_levels
this%file_rhcrit                         = source_file%file_rhcrit
this%file_soil_thickness                 = source_file%file_soil_thickness
this%file_zsea_theta                     = source_file%file_zsea_theta
this%file_c_theta                        = source_file%file_c_theta
this%file_zsea_rho                       = source_file%file_zsea_rho
this%file_c_rho                          = source_file%file_c_rho
this%time_1                              = source_file%time_1
this%time_2                              = source_file%time_2
this%time_3                              = source_file%time_3
this%num_soil_levels                     = source_file%num_soil_levels
this%num_cloud_levels                    = source_file%num_cloud_levels
this%num_tracer_levels                   = source_file%num_tracer_levels
this%num_boundary_layer_levels           = source_file%num_boundary_layer_levels
this%algorithm_to_generate_height_fields =                                    &
                                 source_file%algorithm_to_generate_height_fields
this%num_land_points                     = source_file%num_land_points
this%num_soil_moisture_levels            = source_file%num_soil_moisture_levels
this%variable_resolution                 = source_file%variable_resolution
this%pole_lat                            = source_file%pole_lat
this%pole_long                           = source_file%pole_long
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE copy_common_metadata

!-------------------------------------------------------------------------------
SUBROUTINE validate_grid(this, field_num)
! Validate that the grid this field contains is consistent with the grid in
! the file object
USE stashmaster_constants_mod, ONLY: u_points, v_points, p_points
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: field_num
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='VALIDATE_GRID'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Ensure that the grid in the lookup (used later for the dimensions of the data)
! matches the grid in the integer constants (used for the interpolation)
  IF (this%fields(field_num)%get_horiz_grid_type() == u_points) THEN
    IF(this%fields(field_num)%grid%get_num_cols() /=                           &
                                               this%u_grid%get_num_cols()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_cols(),                         &
           ' columns but the integer constants',                               &
           ' indicate it should have ',this%u_grid%get_num_cols()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
    IF(this%fields(field_num)%grid%get_num_rows() /=                           &
                                               this%u_grid%get_num_rows()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_rows(),                         &
           ' rows but the integer constants',                                  &
           ' indicate it should have ',this%u_grid%get_num_rows()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
  ELSE IF (this%fields(field_num)%get_horiz_grid_type() == v_points) THEN
    IF(this%fields(field_num)%grid%get_num_cols() /=                           &
                                               this%v_grid%get_num_cols()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_cols(),                         &
           ' columns but the integer constants',                               &
           ' indicate it should have ',this%v_grid%get_num_cols()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
    IF(this%fields(field_num)%grid%get_num_rows() /=                           &
                                              this%v_grid%get_num_rows()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_rows(),                         &
           ' rows but the integer constants',                                  &
           ' indicate it should have ',this%v_grid%get_num_rows()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
  ELSE IF (this%fields(field_num)%get_horiz_grid_type() == p_points) THEN
    IF(this%fields(field_num)%grid%get_num_cols() /=                           &
                                               this%p_grid%get_num_cols()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_cols(),                         &
           ' columns but the integer constants',                               &
           ' indicate it should have ',this%p_grid%get_num_cols()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
    IF(this%fields(field_num)%grid%get_num_rows() /=                           &
                                               this%p_grid%get_num_rows()) THEN
      WRITE(cmessage,'(A,I0,A,I0,A,A,I0)') 'Field with STASH ',                &
           this%fields(field_num)%quantity_ident, ' has ',                     &
           this%fields(field_num)%grid%get_num_rows(),                         &
           ' rows but the integer constants',                                  &
           ' indicate it should have ',this%p_grid%get_num_rows()
      icode = 100
      CALL ereport(ModuleName//':'//routinename, icode, cmessage)
    END IF
  END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE validate_grid
!-------------------------------------------------------------------------------
SUBROUTINE modify_packing_method(this, field_num, new_packing_method)

USE stashmaster_utils_mod, ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: packing_acc
USE Packing_Codes_Mod, ONLY: PC_WGDOS_Packing, PC_No_Packing

IMPLICIT NONE
CLASS(fieldsfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: field_num
INTEGER, INTENT(IN) :: new_packing_method
REAL :: accuracy = rmdi
CHARACTER(LEN=*), PARAMETER :: routinename = 'MODIFY_PACKING_METHOD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (new_packing_method == PC_WGDOS_Packing) THEN
  accuracy =  &
      query_stashmaster( this%fields(field_num)%quantity_ident, packing_acc)
ELSE IF (new_packing_method == PC_No_Packing) THEN
  accuracy = -99.0
END IF

CALL this%fields(field_num)%modify_packing_method(new_packing_method, accuracy)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE  modify_packing_method
!-------------------------------------------------------------------------------
END MODULE fieldsfile_mod
