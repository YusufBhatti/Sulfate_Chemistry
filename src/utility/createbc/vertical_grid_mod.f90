! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE vertical_grid_mod

USE missing_data_mod, ONLY: imdi

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   A class to describe a vertical one-dimensional grid. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'VERTICAL_GRID_MOD'

TYPE, PUBLIC :: vertical_grid_type
  PRIVATE
  REAL, ALLOCATABLE, PUBLIC :: level_values(:) ! A general array to contain values
                                               ! which vary with level number
  REAL, ALLOCATABLE, PUBLIC :: eta_theta(:)    ! Eta values for theta levels
  REAL, ALLOCATABLE, PUBLIC :: eta_rho(:)      ! Eta values for rho levels
  REAL, PUBLIC              :: height_at_top_theta_level
  INTEGER, ALLOCATABLE, PUBLIC :: levels(:)    ! Level numbers
  INTEGER, PUBLIC :: num_levels      ! Number of levels for which there is data, ie 1 for single level field
  INTEGER, PUBLIC :: num_model_levels    ! Number of model levels, could be different from num_levels
  INTEGER, PUBLIC  :: coordinate_type ! Vertical co-ordinate type lookup item 26
                                     ! i.e. height, pressure, freezing level
  INTEGER, PUBLIC :: level_type ! Type of levels, i.e. rho, theta, single level
  INTEGER, PUBLIC :: first_level_code ! STASHmaster code (LevelF)
  INTEGER, PUBLIC :: last_level_code  ! STASHmaster code (LevelL)
  INTEGER, PUBLIC :: first_level ! First level on which data is present
  INTEGER, PUBLIC :: last_level  ! Last level on which data is present
  INTEGER, PUBLIC :: first_rho_of_constant_height 
  CONTAINS
  PROCEDURE, PASS :: set_defined_vertical_grid
  PROCEDURE, PASS :: get_defined_vertical_grid
  PROCEDURE, PASS :: get_num_levels
  PROCEDURE, PASS :: set_num_levels
  PROCEDURE, PASS :: get_num_model_levels
  PROCEDURE, PASS :: set_num_model_levels
  PROCEDURE, PASS :: calc_first_level
  PROCEDURE, PASS :: calc_last_level
  PROCEDURE, PASS :: compare_grids
  PROCEDURE, PASS :: recalc_num_levels
  PROCEDURE, PASS :: check_num_levels
END TYPE vertical_grid_type  

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE set_defined_vertical_grid(this, nz, num_model_levels, list_z, levels, &
     first_rho_of_constant_height, height_at_top_theta_level,                    &
     list_eta_theta, list_eta_rho)
! This can be given a list of heights (real) or levels (integer) or both
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: nz
INTEGER, INTENT(IN), OPTIONAL :: num_model_levels
REAL, INTENT(IN), OPTIONAL :: list_z(nz)
INTEGER, INTENT(IN), OPTIONAL :: levels(nz)
INTEGER, INTENT(IN), OPTIONAL :: first_rho_of_constant_height
REAL, INTENT(IN), OPTIONAL :: height_at_top_theta_level
REAL, INTENT(IN), OPTIONAL :: list_eta_theta(:)
REAL, INTENT(IN), OPTIONAL :: list_eta_rho(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_DEFINED_VERTICAL_GRID'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

this%num_levels = nz
  
IF (PRESENT(num_model_levels)) THEN
  this%num_model_levels = num_model_levels
END IF

IF (PRESENT(list_z)) THEN
  ALLOCATE(this%level_values(this%num_levels))
  this%level_values = list_z
END IF

IF (PRESENT(levels)) THEN  
  ALLOCATE(this%levels(this%num_levels))
  this%levels = levels
END IF

IF (PRESENT(first_rho_of_constant_height)) THEN
  this%first_rho_of_constant_height = first_rho_of_constant_height
END IF

IF (PRESENT(height_at_top_theta_level)) THEN
  this%height_at_top_theta_level = height_at_top_theta_level
END IF

IF (PRESENT(list_eta_theta)) THEN
  ! Eta theta values are dimensioned by model levels+1
  ALLOCATE(this%eta_theta(0:this%num_model_levels))
  this%eta_theta = list_eta_theta
END IF

IF (PRESENT(list_eta_rho)) THEN
   ! Eta rho values are dimensioned by model level
  ALLOCATE(this%eta_rho(this%num_model_levels))
  this%eta_rho = list_eta_rho
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_defined_vertical_grid

!-------------------------------------------------------------------------------

SUBROUTINE get_defined_vertical_grid(this, nz, num_model_levels, list_z, levels, &
     first_rho_of_constant_height, height_at_top_theta_level,                    &
     list_eta_theta, list_eta_rho)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(OUT) :: nz
INTEGER, INTENT(OUT), OPTIONAL :: num_model_levels
REAL,    INTENT(OUT), OPTIONAL :: list_z(this%num_levels)
INTEGER, INTENT(OUT), OPTIONAL :: levels(this%num_levels)
INTEGER, INTENT(OUT), OPTIONAL :: first_rho_of_constant_height
REAL,    INTENT(OUT), OPTIONAL :: height_at_top_theta_level
REAL,    INTENT(OUT), OPTIONAL :: list_eta_theta(this%num_levels+1)
REAL,    INTENT(OUT), OPTIONAL :: list_eta_rho(this%num_levels)
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_DEFINED_VERTICAL_GRID'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nz = this%num_levels

IF (PRESENT(list_z)) THEN
  list_z = this%level_values
END IF

IF (PRESENT(num_model_levels)) THEN
  num_model_levels = this%num_model_levels
END IF
  
IF (PRESENT(levels)) THEN
  levels = this%levels
END IF

IF (PRESENT(first_rho_of_constant_height)) THEN
  first_rho_of_constant_height = this%first_rho_of_constant_height
END IF

IF (PRESENT(height_at_top_theta_level)) THEN
  height_at_top_theta_level = this%height_at_top_theta_level
END IF

IF (PRESENT(list_eta_theta)) THEN
  list_eta_theta = this%eta_theta
END IF

IF (PRESENT(list_eta_rho)) THEN
  list_eta_rho = this%eta_rho
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_defined_vertical_grid

!-------------------------------------------------------------------------------

SUBROUTINE calc_first_level(this, grid_stagger)
USE fieldsfile_constants_mod, ONLY: endgame
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_stagger
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_FIRST_LEVEL'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Level code mapping comes from STASHmaster LevelF code. Please see UMDP C4 for
! more details.
IF (this%first_level_code == 1) THEN       ! First atmospheric level
  this%first_level = 1
ELSE IF (this%first_level_code == 38) THEN ! Zeroth atmospheric level (Surface)
  this%first_level = 0
ELSE IF (this%first_level_code == 40) THEN
  IF ( grid_stagger == endgame) THEN   ! ENDGame grid stagger
    this%first_level = 0
  ELSE
    this%first_level = 1
  END IF
ELSE IF (this%first_level_code == -1) THEN  ! Single level field
  this%first_level = 1
ELSE
  this%first_level = imdi  ! Level code not currently supported
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_first_level

!-------------------------------------------------------------------------------

SUBROUTINE calc_last_level(this, grid_stagger)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: grid_stagger
CHARACTER(LEN=*), PARAMETER :: routinename = 'CALC_LAST_LEVEL'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Level code mapping comes from STASHmaster LevelL code. Please see UMDP C4 for
! more details.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (this%last_level_code == 2) THEN          ! Model levels
  this%last_level = this%num_model_levels
ELSE IF (this%last_level_code == 3) THEN     ! Currently wet = model levels
  this%last_level = this%num_model_levels
ELSE IF (this%last_level_code == 19) THEN    ! Model levels + 1
  this%last_level = this%num_model_levels + 1
ELSE IF (this%last_level_code == -1) THEN    ! Single level
  this%last_level = 1
ELSE IF (this%last_level_code == 11) THEN    ! Tracer level (= model levels)
  this%last_level = this%num_model_levels
ELSE
  this%last_level = imdi  ! Level code not currently supported
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calc_last_level

!-------------------------------------------------------------------------------

SUBROUTINE recalc_num_levels(this)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'RECALC_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

this%num_levels = this%last_level - this%first_level + 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE recalc_num_levels

!-------------------------------------------------------------------------------

SUBROUTINE check_num_levels(this)
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'CHECK_NUM_LEVELS'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (this%last_level - this%first_level + 1 /= this%num_levels) THEN
  WRITE(cmessage, '(A,I10,A,I10,A,I10)') "Number of levels: ",this%num_levels, &
       " not consistant with values of first level: ", this%first_level,       &
       " and last level: ", this%last_level
  icode = 10
  CALL ereport(ModuleName//':'//routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_num_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_num_levels(this, num_levels)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

this%num_levels = num_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_levels

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_levels(this)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

get_num_levels = this%num_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_levels

!-------------------------------------------------------------------------------

SUBROUTINE set_num_model_levels(this, num_model_levels)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(INOUT) :: this
INTEGER, INTENT(IN) :: num_model_levels
CHARACTER(LEN=*), PARAMETER :: routinename = 'SET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%num_model_levels = num_model_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_num_model_levels

!-------------------------------------------------------------------------------

INTEGER FUNCTION get_num_model_levels(this)
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(IN) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'GET_NUM_MODEL_LEVELS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
get_num_model_levels = this%num_model_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_num_model_levels

!-------------------------------------------------------------------------------

SUBROUTINE compare_grids(this, compare_grid)
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
CLASS(vertical_grid_type), INTENT(IN) :: this
CLASS(vertical_grid_type), INTENT(IN) :: compare_grid
! Error reporting variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'COMPARE_GRIDS'

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (this%num_levels /= compare_grid%num_levels) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "Number of levels not consistant "//                   &
       "between the two grid objects. Grid 1: ",this%num_levels,                      &
       " Grid 2: ",compare_grid%num_levels
  icode = 10
  CALL ereport(ModuleName//':'//routinename, icode, cmessage)
END IF
IF (this%first_rho_of_constant_height /= compare_grid%first_rho_of_constant_height) THEN
  WRITE(cmessage, '(A,I0,A,I0)') "First rho level of constant height not consistent "// &
       "between the two grid objects. Grid 1: ",this%first_rho_of_constant_height,    &
       " Grid 2: ",compare_grid%first_rho_of_constant_height
  icode = 20
  CALL ereport(ModuleName//':'//routinename, icode, cmessage)
END IF
IF (ABS(this%height_at_top_theta_level - compare_grid%height_at_top_theta_level)    &
     > EPSILON(1.0)) THEN  
  WRITE(cmessage, '(A,F20.10,A,F20.10)') "Height at top theta level not consistent "//  &
       "between the two grid objects. Grid 1: ",this%height_at_top_theta_level,       &
       " Grid 2: ",compare_grid%height_at_top_theta_level
  icode = 30
  CALL ereport(ModuleName//':'//routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE compare_grids

!-------------------------------------------------------------------------------

END MODULE vertical_grid_mod
