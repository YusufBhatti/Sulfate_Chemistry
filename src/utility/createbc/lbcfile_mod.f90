! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE lbcfile_mod

USE fieldsfile_mod,    ONLY: fieldsfile_type, new_dynamics, endgame
USE missing_data_mod,  ONLY: imdi, rmdi
USE data_location_mod, ONLY: data_location_type
USE field_mod,         ONLY: field_type, ASSIGNMENT(=)
USE io,                ONLY: setpos, buffout, buffin
USE ereport_mod,       ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE lbc_output_control_mod, ONLY: packing_unpacked_accuracy

USE lookup_addresses, ONLY: lbpack

USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing, PC_No_Packing

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Description:
!   A class representing LBC files. This contains all the LBC file-specific
!   information. See documentation UMDP F3 for further details on the LBC file
!   format. Quite a lot of the overhead is provided by fieldsfile_type, which
!   this class is a subclass of.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LBCFILE_MOD'

TYPE, PUBLIC, EXTENDS(fieldsfile_type) :: lbcfile_type
  PRIVATE  

  INTEGER :: rim_width
  INTEGER, ALLOCATABLE :: sorted_fields(:)
  
  CONTAINS
  PROCEDURE, PASS :: read_header
  PROCEDURE, PASS :: read_field
  PROCEDURE, PASS :: write_header
  PROCEDURE, PASS :: write_field
  PROCEDURE, PASS :: set_rim_width
END TYPE lbcfile_type

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE read_header(this)
USE field_mod, ONLY: field_type, ASSIGNMENT(=)
USE stashmaster_utils_mod, ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: vertical_coord_hybrid_height,         &
                                     horiz_grid_type, vert_level_type,     &
                                     u_points, v_points, theta_levels,     &
                                     rho_levels
USE umPrintMgr, ONLY: umPrint, umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
CLASS(lbcfile_type), INTENT(INOUT) :: this

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='READ_HEADER'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
cmessage = 'Reading LBC files is not implemented.'
icode =875
CALL ereport(ModuleName//':'//routinename, icode, cmessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_header

!-------------------------------------------------------------------------------

SUBROUTINE read_field(this, ifield)

USE stashmaster_constants_mod, ONLY: data_type_real, data_type_integer,     &
                                     data_type_logical

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
CLASS(lbcfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
CHARACTER(LEN=*), PARAMETER :: routinename='READ_FIELD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
cmessage = 'Reading LBC files is not implemented.'
icode =876
CALL ereport(ModuleName//':'//routinename, icode, cmessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_field

!-------------------------------------------------------------------------------

SUBROUTINE write_header(this)
USE stashmaster_utils_mod,     ONLY: query_stashmaster
USE stashmaster_constants_mod, ONLY: meto8_code, pp_field_code
IMPLICIT NONE
CLASS(lbcfile_type), INTENT(INOUT) :: this
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

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='WRITE_HEADER'
INTEGER :: icode
  
INTEGER :: nx, ny, n_p_rows, n_v_rows, n_p_cols, n_u_cols
INTEGER :: n_p_rows_halo, n_v_rows_halo, n_p_cols_halo, n_u_cols_halo
REAL :: startx, starty, dx, dy, pole_lat, pole_long
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL timer( 'lbc:write_header', 5)
! Set the fixed, integer, real headers sizes and sizes of level,row and column dependent
! constants 
CALL this%set_header_sizes()

fixed_header = imdi
fixed_header(1) = imdi
fixed_header(2) = 1
fixed_header(3) = this%vertical_grid_indicator
fixed_header(4) = this%horizontal_grid_indicator
fixed_header(5) = 5
IF (ALLOCATED(this%fields)) THEN
  fixed_header(6) = this%fields(1)%run_ident
  fixed_header(7) = this%fields(1)%expt_ident
  fixed_header(8) = this%fields(1)%calendar
END IF
fixed_header(9) = this%grid_staggering
fixed_header(10) = 1
fixed_header(12) = this%source_model_version
CALL this%process_times(fixed_header(21:26), fixed_header(28:33), fixed_header(35:40))
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
! Unlike fieldsfiles, this needs to be 3D fields as fields are multi-level in LBC files
fixed_header(152) = this%num_fields
fixed_header(160) = this%data_start

! This is set to zero for LBC files
fixed_header(161) = 0

IF (this%num_fields > this%num_reserved_headers) THEN
  cmessage = 'Header too large for data, increase the number of reserved lookup headers'
  icode = 10
  CALL ereport(ModuleName//':'//routinename, icode, cmessage)
END IF
  
CALL setpos(this%unit_num, this%start_of_file)
CALL buffout(this%unit_num, fixed_header, this%fixed_header_length, icode)

! Integer constants
ALLOCATE(integer_constants(this%len_integer_constants))
ALLOCATE(real_constants(this%len_real_constants))
integer_constants(:)    = imdi
real_constants(:) = rmdi
integer_constants(1) = 0
integer_constants(2) = imdi
integer_constants(3) = imdi
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
real_constants(16) = this%file_theta_rho_levels%height_at_top_theta_level

integer_constants(6) = this%p_grid%get_num_cols()
integer_constants(7) = this%p_grid%get_num_rows()

IF (.NOT. this%variable_resolution) THEN    
  CALL this%p_grid%horiz_grid%get_regular_horizontal_grid(nx, ny, startx, starty, dx, dy)
  real_constants(1) = dx
  real_constants(2) = dy
  IF (this%grid_staggering == new_dynamics) THEN
    real_constants(3) = starty
    real_constants(4) = startx
  ELSE IF (this%grid_staggering == endgame) THEN
    real_constants(3) = starty - 0.5 * dy ! Start lat/long is offset from first P
    real_constants(4) = startx - 0.5 * dx ! point for ENDGame
  END IF
END IF

real_constants(5) = this%pole_lat
real_constants(6) = this%pole_long

CALL setpos(this%unit_num, this%integer_constants_start -1 )
CALL buffout(this%unit_num, integer_constants, this%len_integer_constants)
CALL setpos(this%unit_num, this%real_constants_start -1 )
CALL buffout(this%unit_num, real_constants, this%len_real_constants)
DEALLOCATE(real_constants)
DEALLOCATE(integer_constants)

ALLOCATE(level_dep_constants(this%file_theta_rho_levels%get_num_levels()+1, this%len2_lev_dep_constants))
level_dep_constants = rmdi
level_dep_constants(1:this%file_theta_rho_levels%get_num_levels()+1,1) = this%file_theta_rho_levels%eta_theta(:)
level_dep_constants(1:this%file_theta_rho_levels%get_num_levels(),2) = this%file_theta_rho_levels%eta_rho(:)
CALL setpos(this%unit_num, this%level_dep_constants_start - 1)
CALL buffout(this%unit_num, level_dep_constants, this%level_dep_constants_dimension1 * this%level_dep_constants_dimension2)
DEALLOCATE(level_dep_constants)

IF (this%variable_resolution) THEN
  ALLOCATE(row_dep_constants(this%row_dep_constants_dimension1, this%row_dep_constants_dimension2))
  ALLOCATE(column_dep_constants(this%column_dep_constants_dimension1, this%column_dep_constants_dimension2))
  row_dep_constants = rmdi
  column_dep_constants = rmdi

  row_dep_constants(1:this%p_grid%get_num_rows(),1) = this%p_grid%get_latitudes()
  row_dep_constants(1:this%v_grid%get_num_rows(),2) = this%v_grid%get_latitudes()
  column_dep_constants(1:this%p_grid%get_num_cols(),1) = this%p_grid%get_longitudes()
  column_dep_constants(1:this%u_grid%get_num_cols(),2) = this%u_grid%get_longitudes()

  CALL setpos(this%unit_num, this%row_dep_constants_start -1 )
  CALL buffout(this%unit_num, row_dep_constants, this%row_dep_constants_dimension1*this%row_dep_constants_dimension2)
  CALL setpos(this%unit_num, this%column_dep_constants_start -1 )
  CALL buffout(this%unit_num, column_dep_constants, this%column_dep_constants_dimension1*this%column_dep_constants_dimension2)
  DEALLOCATE(column_dep_constants)
  DEALLOCATE(row_dep_constants)
END IF

! Construct lookup header
ALLOCATE(lookup(this%len_single_lookup, this%num_reserved_headers))
lookup = -99

k = 1
DO j = 1, this%num_fields
  ! If the fields have been sorted, write the header in that order
  IF (ALLOCATED(this%sorted_fields)) THEN
    i = this%sorted_fields(j)
  ELSE
    i = j
  END IF
 
  lookup(:,k) = 0
  lookup(1:6,k) = this%fields(i)%validity_time(1:6)
  lookup(7:12,k) = this%fields(i)%data_time(1:6)
  lookup(13,k) = this%fields(i)%time_indicator
  lookup(14,k) = this%fields(i)%fctime
  lookup(16,k) = this%fields(i)%get_grid_code()
  lookup(17,k) = 100 + this%fields(i)%get_num_levels()
  lookup(lbpack,k) = this%fields(i)%packing_method
  lookup(22,k) = this%fields(i)%calc_header_release()
  lookup(23,k) = query_stashmaster(this%fields(i)%quantity_ident,pp_field_code)
  lookup(25,k) = this%fields(i)%processing_code
  lookup(28,k) = 0
  lookup(32,k) = query_stashmaster(this%fields(i)%quantity_ident,meto8_code)
  lookup(38,k) = this%fields(i)%source_model_identifier
  lookup(39,k) = this%fields(i)%data_type
  lookup(42,k) = this%fields(i)%quantity_ident
  lookup(43,k) = this%fields(i)%meaning_of_k_dimension
  lookup(45,k) = 1
    
  lookup(41,k) = this%rim_width * 10000 + this%fields(i)%get_halo_ns() * 100 &
       + this%fields(i)%get_halo_ew()

  real_part_of_lookup = 0.0
  real_part_of_lookup(5) = this%fields(i)%datum_constant
  this%fields(i)%packing_accuracy = 0
  real_part_of_lookup(6) = this%fields(i)%packing_accuracy
  real_part_of_lookup(18) = this%fields(i)%mdi

  CALL this%fields(i)%get_pole(pole_lat, pole_long)
  real_part_of_lookup(11) = pole_lat
  real_part_of_lookup(12) = pole_long

  ! Lookup expects rows and row_length not including halos
  lookup(18,k) = this%fields(i)%get_num_rows() - this%fields(i)%get_halo_ns()*2
  lookup(19,k) = this%fields(i)%get_num_cols() - this%fields(i)%get_halo_ew()*2

  real_part_of_lookup(13) = 0 ! BGOR, not used set to zero

  IF (this%variable_resolution) THEN
    real_part_of_lookup(14) = rmdi
    real_part_of_lookup(15) = rmdi
    real_part_of_lookup(16) = rmdi
    real_part_of_lookup(17) = rmdi
  ELSE
    ! The LBC grids include any possible halos so need to take this into account
    ! when calculating zeroth positions for lookup
    CALL this%fields(i)%grid%horiz_grid%get_regular_horizontal_grid(nx, ny, startx, starty, dx, dy)
    real_part_of_lookup(14) = starty + ((this%fields(i)%get_halo_ns() - 1) * dy)
    real_part_of_lookup(15) = dy
    real_part_of_lookup(16) = startx + ((this%fields(i)%get_halo_ew()- 1) * dx)
    real_part_of_lookup(17) = dx
  END IF

  ! Copy the REAL part of the lookup into INTEGER
  DO p = 1, 19
    lookup(p+45,k) = TRANSFER(real_part_of_lookup(p), 0)
  END DO

  ! LBCs have 7777 as they're multilevel fields
  lookup(33,k) = 7777

  ! Disk locations - these should be set by write_fields
  lookup(29,k) = this%data_locations(i)%data_start(1)
  lookup(15,k) = this%data_locations(i)%data_length(1)
  lookup(30,k) = this%data_locations(i)%disk_length(1)
  IF (k == 1) THEN
    lookup(40,k) = 1
  ELSE
    lookup(40,k) = lookup(40,k-1) + lookup(15,k-1)
  END IF
  this%data_locations(i)%lookup_position(1) = k
  k = k + 1
END DO

! Write lookup   
CALL setpos(this%unit_num, this%lookup_start -1 )
CALL buffout(this%unit_num, lookup, this%len_single_lookup*this%num_reserved_headers)

DEALLOCATE(lookup)
  
CALL timer( 'lbc:write_header', 6)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE write_header

!-------------------------------------------------------------------------------

SUBROUTINE write_field(this, ifield)
USE io_configuration_mod, ONLY: &
    io_field_padding
USE cray32_packing_mod, ONLY: lbc_write_packed_field
IMPLICIT NONE
CLASS(lbcfile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield
INTEGER :: position_to_write, j, len_data, len_disk , k
INTEGER :: l,m,n,q
REAL, ALLOCATABLE :: work_array(:)
INTEGER :: nx, ny, idum, nz, ii
LOGICAL, PARAMETER :: compress = .TRUE.

INTEGER :: num_words, temp_packing, num_points
  
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='WRITE_FIELD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(this%fields)) THEN
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (.NOT. ALLOCATED(this%fields(ifield)%lbc_rdata)) THEN
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

CALL timer( 'lbc:write_field', 5)

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
  cmessage = "Unable to write field to file. The value of data_start in the file " // &
       "object is <= 0"
  CALL ereport(routinename, icode, cmessage)
END IF

nx = this%fields(ifield)%get_num_cols()
ny = this%fields(ifield)%get_num_rows()
nz = this%fields(ifield)%get_num_levels()
num_points = nx * ny

len_data = SIZE(this%fields(ifield)%lbc_rdata)

IF (MOD(this%fields(ifield)%packing_method, 10) == PC_No_Packing) THEN
  ! Unpacked data
  len_disk = len_data
  len_disk=((len_disk+io_field_padding-1)/                                 &
                   io_field_padding)*io_field_padding

  ALLOCATE(work_array(len_disk))
  work_array = rmdi
  m = 0
  DO l = 1, SIZE(this%fields(ifield)%lbc_rdata, 2)
    DO k = 1, SIZE(this%fields(ifield)%lbc_rdata, 1)
      m = m + 1
      work_array(m) = this%fields(ifield)%lbc_rdata(k,l)
    END DO
  END DO

  CALL setpos(this%unit_num, position_to_write)
  CALL buffout(this%unit_num, work_array, len_disk)
  this%data_locations(ifield)%data_start(1) = position_to_write
  this%data_locations(ifield)%data_length(1) = len_data
  this%data_locations(ifield)%disk_length(1) = len_disk
  position_to_write = position_to_write + len_disk
  DEALLOCATE(work_array)
ELSE
  ! CRAY 32-bit packing
  IF (MOD(this%fields(ifield)%packing_method, 10) /= PC_Cray32_Packing) THEN
    ! An unsupported packing method
    icode = -10
    WRITE(cmessage, '(A,I2,A)')  'Warning, this field needs to be packed'  &
         // ' using an unsupported method. LBPACK = ',                     &
         this%fields(ifield)%packing_method, ', attempting to CRAY32 pack' &
         // ' instead'
    CALL ereport(routinename, icode, cmessage)
    temp_packing = this%fields(ifield)%packing_method / 10
    this%fields(ifield)%packing_method = temp_packing * 10 + PC_Cray32_Packing
  END IF

  ! Number of 64-bit words required to store the 32-bit data
  len_disk = (len_data+1)/2

  ! round up to sector boundary
  len_disk=((len_disk+io_field_padding-1)/                                 &
                   io_field_padding)*io_field_padding

  CALL lbc_write_packed_field(this%unit_num, this%fields(ifield)%lbc_rdata,          &
                    len_data, position_to_write, len_disk)
  this%data_locations(ifield)%data_start(1) = position_to_write
  this%data_locations(ifield)%data_length(1) = len_data
  this%data_locations(ifield)%disk_length(1) = len_disk
  position_to_write = position_to_write + len_disk
END IF

this%most_recent_data_position = position_to_write 

CALL timer( 'lbc:write_field', 6)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE write_field

!-------------------------------------------------------------------------------

SUBROUTINE set_rim_width(this, rim_width)
IMPLICIT NONE
CLASS(lbcfile_type) :: this
INTEGER, INTENT(IN) :: rim_width
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_RIM_WIDTH'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%rim_width = rim_width

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_rim_width

!-------------------------------------------------------------------------------

END MODULE lbcfile_mod
