! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE datafile_mod

USE field_mod,                  ONLY: field_type, ASSIGNMENT(=)
USE file_mod,                   ONLY: file_type
USE three_dimensional_grid_mod, ONLY: three_dimensional_grid_type
USE vertical_grid_mod,          ONLY: vertical_grid_type
USE fieldsfile_constants_mod,   ONLY: endgame, new_dynamics
USE missing_data_mod,           ONLY: imdi, rmdi
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   A base class for datafiles. 
!     Subclasses of this class represent different types of files, e.g.
!     fieldsfiles, GRIB2, etc. Note the read/write file/header procedures must 
!     be defined by every subclass - this is enforced by using the DEFERRED
!     keyword.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DATAFILE_MOD'

TYPE, ABSTRACT, PUBLIC, EXTENDS(file_type) :: datafile_type

  INTEGER :: num_fields

  TYPE(field_type), ALLOCATABLE :: fields(:)

  LOGICAL :: l_source_rotated = .FALSE.   ! Input is on rotated grid
  INTEGER :: grid_staggering ! New Dynamics or ENDGame grid staggering

  ! The three potential grids contained in this file 
  TYPE(three_dimensional_grid_type) :: p_grid, v_grid, u_grid

  ! The main vertical grid associated with this file and number of model levels
  INTEGER :: num_model_levels = imdi
  TYPE(vertical_grid_type) :: file_theta_rho_levels

  ! The rotated pole locations for the fields in this file
  REAL    :: pole_lat = rmdi
  REAL    :: pole_long = rmdi

  ! Is the file grid cyclic? ie global or wrapping LAM
  LOGICAL :: l_source_cyclic = .FALSE.

  CONTAINS

  PROCEDURE(read_header_interface),     DEFERRED, PASS  :: read_header
  PROCEDURE(write_header_interface),    DEFERRED, PASS  :: write_header
  PROCEDURE(read_field_interface),      DEFERRED, PASS  :: read_field
  PROCEDURE(write_field_interface),     DEFERRED, PASS  :: write_field
  PROCEDURE(add_field_interface),       DEFERRED, PASS  :: add_field
  PROCEDURE(allow_read_only_interface), DEFERRED, PASS  :: allow_read_only
  PROCEDURE(set_header_sizes_interface),DEFERRED, PASS  :: set_header_sizes
  PROCEDURE, PASS :: read_all_fields
  PROCEDURE, PASS :: write_all_fields
  PROCEDURE, PASS :: unload_field

END TYPE datafile_type

ABSTRACT INTERFACE
SUBROUTINE read_header_interface(this)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
END SUBROUTINE read_header_interface

SUBROUTINE write_header_interface(this)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
END SUBROUTINE write_header_interface

SUBROUTINE read_field_interface(this, ifield)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield
END SUBROUTINE read_field_interface

SUBROUTINE write_field_interface(this, ifield)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: ifield
END SUBROUTINE write_field_interface

! This should return the field number the new field was added as.
INTEGER FUNCTION add_field_interface(this, new_field)
IMPORT :: datafile_type, field_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
TYPE(field_type), INTENT(IN) :: new_field
END FUNCTION add_field_interface
  
SUBROUTINE allow_read_only_interface(this, read_only)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
LOGICAL, INTENT(IN) :: read_only
END SUBROUTINE allow_read_only_interface

SUBROUTINE set_header_sizes_interface(this)
IMPORT :: datafile_type
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
END SUBROUTINE set_header_sizes_interface
END INTERFACE

CONTAINS

! The following may not be required for createbc but may be useful for other
! programs. 

!-------------------------------------------------------------------------------

SUBROUTINE read_all_fields(this)
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'READ_ALL_FIELDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, this%num_fields
  CALL this%read_field(i)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_all_fields

!-------------------------------------------------------------------------------

SUBROUTINE write_all_fields(this)
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'WRITE_ALL_FIELDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, this%num_fields
  CALL this%write_field(i)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE write_all_fields

!-------------------------------------------------------------------------------

SUBROUTINE unload_field(this, ifield)
IMPLICIT NONE
CLASS(datafile_type), INTENT(INOUT) :: this
INTEGER :: ifield
CHARACTER(LEN=*), PARAMETER :: routinename = 'UNLOAD_FIELD'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL this%fields(ifield)%unload_data()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE unload_field  

!-------------------------------------------------------------------------------

END MODULE datafile_mod
