! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module containing the file management types and routines
!              which handle the assignment and release of file units
!              and storage of metadata relating to files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: IO Services
!
! Code Description:
!   Language: FORTRAN 90

MODULE file_manager

#if !defined(UTILIO)
USE netcdf, ONLY: nf90_max_dims    ! External netCDF library module
#endif
USE filenamelength_mod, ONLY: filenamelength
USE missing_data_mod, ONLY: imdi
USE ereport_mod, ONLY: ereport
USE umprintmgr, ONLY: umprint, ummessage
USE um_types, ONLY: integer32

IMPLICIT NONE 

PRIVATE

! Every file accessed by the model will be assigned an element of this type
!-------------------------------------------------------------------------------
TYPE um_file_type

  ! The most basic information about a file, for simple files such
  ! as text-output or files to read namelists from this will suffice
  INTEGER :: UNIT = imdi
  CHARACTER(LEN=filenamelength) :: filename = ""

  ! Unique identifier string of file 
  CHARACTER(LEN=20) :: id = ""
  
  ! For more complex model input/output files a variety of extra 
  ! data may be required.  These pointers can be used to add that
  ! information to a given element

  ! Extra data common to all output formats
  TYPE(um_file_meta_type), POINTER :: meta => NULL()

  ! Extra data needed for pp output only
  TYPE(pp_file_meta_type), POINTER :: pp_meta => NULL()

  ! Extra data needed for netCDF output only
  TYPE(nc_file_meta_type), POINTER :: nc_meta => NULL()

  ! Pointers for the linked list
  TYPE(um_file_type), POINTER :: next => NULL()
  TYPE(um_file_type), POINTER :: prev => NULL()

END TYPE um_file_type

! As described above these types stores the various bits of additional
! information required by the more complex model input/output files
!-------------------------------------------------------------------------------

! Information common to all output formats
TYPE um_file_meta_type

  ! Indicates if the file is an output stream file
  LOGICAL :: is_output_file = .FALSE.
  ! Indicates if the file is being provided with partially written
  ! data (for a CRUN) 
  LOGICAL :: partially_written = .FALSE.
  ! Indicates whether this file is due to be initialised this timestep
  LOGICAL :: initialise = .FALSE.
  ! Frequency of initialisation in time-steps
  INTEGER :: init_steps = 0
  INTEGER :: init_start_step = 0    ! Starting only after this timestep
  INTEGER :: init_end_step = 0      ! ... and finishing at this timestep
  ! Indicates type of packing to use
  INTEGER :: packing_code = imdi
  ! Base name for filename (if to be used with appended timestep info)
  CHARACTER(LEN=filenamelength) :: filename_base = ""

END TYPE um_file_meta_type

! Information needed for pp output only
TYPE pp_file_meta_type

  ! Number of field headers reservered (non-mean ouput files)
  ! These may be calculated and/or set by the user
  INTEGER :: reserved_headers = 0
  INTEGER :: reserved_headers_calc = 0
  ! If partially written, gives the position of the final written field 
  INTEGER :: last_written_field = 0
  ! File/stream type indicator
  CHARACTER(LEN=1) :: init_file_type = ""

  ! Attributes taken from model_file (IOS)
  LOGICAL           :: managed = .FALSE.
  INTEGER           :: lookup_address = imdi
  INTEGER, POINTER  :: fixed_header(:)   => NULL() ! Fixed Length Header
  INTEGER, POINTER  :: lookup_table(:,:) => NULL()

END TYPE pp_file_meta_type

! Information needed for netCDF output only
TYPE nc_file_meta_type

  ! NetCDF 4 compression switch
  LOGICAL :: l_compress = .FALSE.
  ! Deflate level, value between 0 and 9
  INTEGER :: nccomp_level = imdi
  ! Switch to turn on netCDF4 shuffle filter
  LOGICAL :: l_shuffle = .FALSE.
  ! NetCDF id
  INTEGER(integer32) :: ncid = -1
  ! NetCDF open/create mode parameter
  INTEGER(integer32) :: mode = 0
  ! Count of the number of fields of a particular STASH section/item 
  INTEGER :: stash_count = 0
#if !defined(UTILIO)
  ! STASH index for a given dimension id
  INTEGER :: stash_index_dim(nf90_max_dims) = -1
  ! longitude dimension id given the corresponding latitude dimension id
  INTEGER :: dimid_lon(nf90_max_dims) = -1
#endif

END TYPE nc_file_meta_type

! The type used to store and manage the list
!-------------------------------------------------------------------------------
TYPE um_file_list_type
  ! The number of file units currently in the list
  INTEGER :: n_file_units = 0
  ! The allowed range limits for the units
  INTEGER :: file_unit_min = 0
  INTEGER :: file_unit_max = 0
  ! A logical array to store the status of the units
  LOGICAL, ALLOCATABLE :: unit_in_use(:)
  INTEGER :: highest_unit_used = 0
  ! Pointers to the first and last elements in the list
  TYPE(um_file_type), POINTER :: head
  TYPE(um_file_type), POINTER :: tail
END TYPE um_file_list_type

! The instances of the type used as the master list objects by this module, 
! there are two distinct file lists depending on the type of file
!-------------------------------------------------------------------------------
TYPE(um_file_list_type), SAVE, TARGET :: um_file_list_fortran
TYPE(um_file_list_type), SAVE, TARGET :: um_file_list_portio
TYPE(um_file_list_type), SAVE, TARGET :: um_file_list_netcdf

! A flag which forces the portio list to use unique unit numbers
LOGICAL, PUBLIC :: portio_unique_units = .FALSE.

PUBLIC :: init_file_manager
PUBLIC :: stop_file_manager
PUBLIC :: assign_file_unit
PUBLIC :: release_file_unit
PUBLIC :: get_file_by_id
PUBLIC :: get_file_by_unit
PUBLIC :: print_file_list
PUBLIC :: get_file_unit_by_id
PUBLIC :: um_file_type
PUBLIC :: init_file_loop

INTERFACE release_file_unit
MODULE PROCEDURE &
    release_file_unit_by_id, &
    release_file_unit_by_unit
END INTERFACE

CONTAINS

! Return target pointer for file list (for looping
!-------------------------------------------------------------------------------
FUNCTION init_file_loop(handler) RESULT(list_head)

IMPLICIT NONE 

CHARACTER(*)                :: handler
TYPE(um_file_type), POINTER :: list_head

TYPE(um_file_list_type), POINTER :: um_file_list

! Choose the correct list and return the head of that list
um_file_list => select_list(handler)
list_head => um_file_list % head

END FUNCTION init_file_loop

! Inserts an item at the end of the list
!-------------------------------------------------------------------------------
SUBROUTINE insert_at_end(list, new_elem)

IMPLICIT NONE 

TYPE(um_file_list_type), INTENT(INOUT)   :: list
TYPE(um_file_type), INTENT(INOUT), POINTER :: new_elem

IF (list % n_file_units == 0) THEN
  ! If this is the first element inserted it becomes both the head and tail
  ! of the list, and its previous link should not point to anything
  list % head => new_elem
  NULLIFY(list % head % prev)
  list % tail => list % head
ELSE
  ! Otherwise make it the next link after the final element in the list, 
  ! its previous link should be the current end of the list and then it
  ! can replace that element as the new end of the list
  list % tail % next => new_elem
  new_elem % prev => list % tail
  list % tail => new_elem
END IF
! Make sure the new elements next link is cleared and increment the counter
NULLIFY(list % tail % next)
list % n_file_units = list % n_file_units + 1

END SUBROUTINE insert_at_end

! Inserts an item at the start of the list
!-------------------------------------------------------------------------------
SUBROUTINE insert_at_start(list, new_elem)

IMPLICIT NONE 

TYPE(um_file_list_type), INTENT(INOUT)   :: list
TYPE(um_file_type), INTENT(INOUT), POINTER :: new_elem

IF (list % n_file_units == 0) THEN
  ! If this is the first element inserted it becomes both the head and tail
  ! of the list, and its next link should not point to anything
  list % tail => new_elem
  NULLIFY(list % tail % next)
  list % head => list % tail
ELSE
  ! Otherwise make it the previous link before the first element in the list,
  ! its next link should be the current head of the list and then it can
  ! replace that element as the new head of the list
  list % head % prev => new_elem
  new_elem % next => list % head
  list % head => new_elem
END IF
! Make sure the new elements previous link is cleared and increment the counter
NULLIFY(list % head % prev)
list % n_file_units = list % n_file_units + 1

END SUBROUTINE insert_at_start

! Removes a specific element from the list
!-------------------------------------------------------------------------------
SUBROUTINE remove_element(list, elem)

IMPLICIT NONE 

TYPE(um_file_list_type), INTENT(INOUT)   :: list
TYPE(um_file_type), INTENT(INOUT), POINTER :: elem

IF (list % n_file_units > 1) THEN
  ! If this is not the only element in the list
  IF (.NOT. ASSOCIATED(elem % next)) THEN
    ! If this element is the tail of the list, break the link to it from the
    ! previous item and make that item the new tail of the list
    NULLIFY(elem % prev % next)
    list % tail => elem % prev
  ELSE IF (.NOT. ASSOCIATED(elem % prev)) THEN
    ! If this element is the head of the list, break the link to it from the
    ! next item and make that item the new head of the list
    NULLIFY(elem % next % prev)
    list % head => elem % next
  ELSE
    ! Otherwise it is an element inside the list somewhere, in which case the
    ! elements either side need to be pointed at each other
    elem % prev % next => elem % next
    elem % next % prev => elem % prev
  END IF
ELSE
  ! If this was the last element in the list, clean up the pointers to the
  ! head and tail of the list
  NULLIFY(list % head)
  NULLIFY(list % tail)
END IF

! Clean up the element and decrement the counter
IF (ASSOCIATED(elem % meta)) THEN
  DEALLOCATE(elem % meta)
  NULLIFY(elem % meta)
END IF
IF (ASSOCIATED(elem % pp_meta)) THEN
  DEALLOCATE(elem % pp_meta)
  NULLIFY(elem % pp_meta)
END IF
IF (ASSOCIATED(elem % nc_meta)) THEN
  DEALLOCATE(elem % nc_meta)
  NULLIFY(elem % nc_meta)
END IF
DEALLOCATE(elem)
NULLIFY(elem)
list % n_file_units = list % n_file_units - 1

END SUBROUTINE remove_element

! Prints out some details of the list elements
!-------------------------------------------------------------------------------
SUBROUTINE list_print(list)

IMPLICIT NONE 

TYPE(um_file_list_type), INTENT(INOUT) :: list

TYPE(um_file_type), POINTER :: elem
INTEGER :: i_elem

CALL umprint("FILE_MANAGER: Printing current list-----", &
              src="file_manager:list_print")

elem => list % head
DO i_elem = 1, list % n_file_units
!$OMP CRITICAL(internal_write)
  WRITE(ummessage,"(A,I4)") "FILE_MANAGER: File     : ",i_elem
!$OMP END CRITICAL(internal_write)
  CALL umprint(ummessage, src="file_manager:list_print")
  WRITE(ummessage,"(A,A)")  "FILE_MANAGER: filename : ",TRIM(elem % filename)
  CALL umprint(ummessage, src="file_manager:list_print")
  WRITE(ummessage,"(A,A)")  "FILE_MANAGER: file id  : ",TRIM(elem % id)
  CALL umprint(ummessage, src="file_manager:list_print")
!$OMP CRITICAL(internal_write)
  WRITE(ummessage,"(A,I4)") "FILE_MANAGER: file unit: ",elem % UNIT
!$OMP END CRITICAL(internal_write)
  CALL umprint(ummessage, src="file_manager:list_print")
  elem => elem % next
END DO

END SUBROUTINE list_print

! Initialises a file list
!-------------------------------------------------------------------------------
SUBROUTINE init_file_list(list, unit_min, unit_max)
IMPLICIT NONE 

TYPE(um_file_list_type), INTENT(INOUT) :: list
INTEGER, INTENT(IN) :: unit_min
INTEGER, INTENT(IN) :: unit_max

list % file_unit_min = unit_min
list % file_unit_max = unit_max
list % highest_unit_used = unit_min - 1
ALLOCATE(list % unit_in_use(unit_max))
list % unit_in_use = .FALSE.

END SUBROUTINE init_file_list

! Initialises the file manager
!-------------------------------------------------------------------------------
SUBROUTINE init_file_manager()
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ERROR_UNIT, INPUT_UNIT, OUTPUT_UNIT
IMPLICIT NONE 

INTEGER :: i_unit

INTEGER, PARAMETER :: start_unit_fortran = 10
INTEGER, PARAMETER :: end_unit_fortran = 300

INTEGER, PARAMETER :: start_unit_portio = 10
INTEGER, PARAMETER :: end_unit_portio = 300

INTEGER, PARAMETER :: start_unit_netcdf = 10
INTEGER, PARAMETER :: end_unit_netcdf = 300

! Setup the list limits
CALL init_file_list(um_file_list_fortran, start_unit_fortran, end_unit_fortran)
CALL init_file_list(um_file_list_portio,  start_unit_portio, end_unit_portio)
CALL init_file_list(um_file_list_netcdf,  start_unit_netcdf, end_unit_netcdf)

! Manually turn off the compiler-reserved units for the native list if
! they appear in the usual range as these should not be used
DO i_unit = start_unit_fortran, end_unit_fortran
  IF (     i_unit == ERROR_UNIT                          &
      .OR. i_unit == INPUT_UNIT                          &
      .OR. i_unit == OUTPUT_UNIT   ) THEN
    um_file_list_fortran % unit_in_use(i_unit) = .TRUE.
  END IF
END DO

CALL umprint("FILE_MANAGER: Initialised---------------", &
             src="file_manager:init_file_manager")

END SUBROUTINE init_file_manager

! Terminates the file manager
!-------------------------------------------------------------------------------
SUBROUTINE stop_file_manager()
IMPLICIT NONE 

CALL umprint("FILE_MANAGER: Stopping------------------", &
             src="file_manager:init_file_manager")

END SUBROUTINE stop_file_manager

! Returns a pointer to the desired list
!-------------------------------------------------------------------------------
FUNCTION select_list(NAME) RESULT(ptr)
IMPLICIT NONE 

CHARACTER(*) :: NAME
TYPE(um_file_list_type), POINTER :: ptr

INTEGER :: errorstatus

SELECT CASE (NAME)
CASE ("fortran")
  ptr => um_file_list_fortran
CASE ("portio")
  ptr => um_file_list_portio
CASE ("netcdf")
  ptr => um_file_list_netcdf
CASE DEFAULT
  errorstatus = 1
  CALL ereport("file_manager:select_handler", errorstatus, &
      "Attempt to select an unrecognised handler: "//TRIM(NAME))
END SELECT

END FUNCTION select_list

! For obtaining a unit for a given file - you must provide a filename
! and can optionally request for the extra metadata to be associated with 
! that unit for UM model type files
!-------------------------------------------------------------------------------
SUBROUTINE assign_file_unit(filename, UNIT, handler, id, force)

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE 

CHARACTER(*),           INTENT(IN)  :: filename
INTEGER,                INTENT(INOUT) :: UNIT
CHARACTER(*),           INTENT(IN)  :: handler
CHARACTER(*), OPTIONAL, INTENT(IN)  :: id
INTEGER,      OPTIONAL, INTENT(IN)  :: force

TYPE(um_file_type),      POINTER :: new_entry
TYPE(um_file_meta_type), POINTER :: new_metadata
TYPE(pp_file_meta_type), POINTER :: new_ppmetadata
TYPE(nc_file_meta_type), POINTER :: new_ncmetadata
TYPE(um_file_type),      POINTER :: ptr
TYPE(um_file_list_type), POINTER :: um_file_list

INTEGER :: i_elem
INTEGER :: errorstatus
CHARACTER(LEN=errormessagelength) :: message

! Choose the correct list
um_file_list => select_list(handler)

! Allow the unit number to be forced by the caller (used for testing or in
! a few special cases such as the IO PEs when using the IO server)
IF (PRESENT(force)) THEN
  ! Don't allow a unit to be forced if it is in use
  IF (um_file_list % unit_in_use(UNIT)) THEN
    errorstatus = 1
!$OMP CRITICAL(internal_write)
    WRITE(message,"(A,I3,A)") &
        "Attempt to force assign an in-use unit (",UNIT,")"
!$OMP END CRITICAL(internal_write)
    CALL ereport("file_manager:assign_file_unit", errorstatus, message)
  END IF
  UNIT = force
ELSE
  ! Retrieve a free unit number, do this first in 
  ! case there aren't any available (to crash earlier rather than later)
  IF (TRIM(handler) == "portio" .AND. portio_unique_units) THEN
    CALL get_file_unit(UNIT, um_file_list, unique=.TRUE.)
  ELSE
    CALL get_file_unit(UNIT, um_file_list, unique=.FALSE.)
  END IF
END IF

! Allocate the memory for this new element and set the basic information
ALLOCATE(new_entry)
new_entry % filename = TRIM(filename)
new_entry % UNIT = UNIT

! If a file ID has been provided check that it isn't in use (must be unique)
IF (PRESENT(id)) THEN
  ptr => um_file_list % head
  DO i_elem = 1, um_file_list % n_file_units
    IF (TRIM(ptr % id) == TRIM(id)) THEN
      errorstatus = 1
      CALL ereport("file_manager:assign_file_unit", errorstatus, &
          "Attempt to assign a file with ID already in use ("//TRIM(id)//")")
    END IF
    ptr => ptr%next
  END DO
  new_entry % id = id
END IF

! Add the element to the correct list, and add the extra metadata for
! the files which require it
ALLOCATE(new_metadata)
new_entry % meta => new_metadata
IF (TRIM(handler) == "portio") THEN
  ALLOCATE(new_ppmetadata)
  new_entry % pp_meta => new_ppmetadata
ELSE IF (TRIM(handler) == "netcdf") THEN
  ALLOCATE(new_ncmetadata)
  new_entry % nc_meta => new_ncmetadata
END IF

CALL insert_at_end(um_file_list, new_entry)

CALL umprint("FILE_MANAGER: Assigned : "//TRIM(filename), &
             src="file_manager:assign_file_unit")
IF (PRESENT(id)) THEN
  CALL umprint("FILE_MANAGER:          : id   : "//TRIM(id), &
               src="file_manager:assign_file_unit")
END IF
!$OMP CRITICAL(internal_write)
WRITE(ummessage,"(A,I4,A)") &
    "FILE_MANAGER:          : Unit :",UNIT," ("//TRIM(handler)//")"
!$OMP END CRITICAL(internal_write)
CALL umprint(ummessage, src="file_manager:assign_file_unit")

END SUBROUTINE assign_file_unit

! Deregister a particular file unit from the list using its ID
!-------------------------------------------------------------------------------
SUBROUTINE release_file_unit_by_id(id, handler)
IMPLICIT NONE 

CHARACTER(*), INTENT(IN) :: id
CHARACTER(*), INTENT(IN) :: handler

TYPE(um_file_type), POINTER :: ptr
TYPE(um_file_list_type), POINTER :: um_file_list

um_file_list => select_list(handler)
ptr => get_file_by_id(id, handler)

CALL umprint("FILE_MANAGER: Released : id   : "//&
              TRIM(id)//" ("//TRIM(handler)//")", &
              src="file_manager:release_file_by_id")

CALL free_file_unit(ptr % UNIT, um_file_list)
CALL remove_element(um_file_list, ptr)

END SUBROUTINE release_file_unit_by_id

! Deregister a particular file unit from the list using its unit
!-------------------------------------------------------------------------------
SUBROUTINE release_file_unit_by_unit(UNIT, handler)
IMPLICIT NONE 

INTEGER,      INTENT(IN) :: UNIT
CHARACTER(*), INTENT(IN) :: handler

TYPE(um_file_type), POINTER :: ptr
TYPE(um_file_list_type), POINTER :: um_file_list

um_file_list => select_list(handler)
ptr => get_file_by_unit(UNIT, handler)

!$OMP CRITICAL(internal_write)
WRITE(ummessage,"(A,I4,A)") &
    "FILE_MANAGER: Released : Unit :",UNIT," ("//TRIM(handler)//")"
!$OMP END CRITICAL(internal_write)
CALL umprint(ummessage, src="file_manager:release_file_unit_by_unit")

CALL free_file_unit(ptr % UNIT, um_file_list)
CALL remove_element(um_file_list, ptr)

END SUBROUTINE release_file_unit_by_unit

! Simply returns the file unit given the file id
!-------------------------------------------------------------------------------
FUNCTION get_file_unit_by_id(id, handler, ignore_missing) RESULT(UNIT)
IMPLICIT NONE 

CHARACTER(*) :: id
CHARACTER(*) :: handler
INTEGER      :: UNIT
LOGICAL, OPTIONAL :: ignore_missing

TYPE(um_file_type), POINTER :: ptr

ptr => get_file_by_id(id, handler,ignore_missing)
IF (ASSOCIATED(ptr)) THEN
  UNIT = ptr % UNIT
ELSE
  UNIT = imdi
END IF

END FUNCTION get_file_unit_by_id

! Assigns a pointer to a file object based on the file unit number
!-------------------------------------------------------------------------------
FUNCTION get_file_by_unit(UNIT, handler, ignore_missing) RESULT(ptr)

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE 

INTEGER :: UNIT
CHARACTER(*) :: handler
LOGICAL, OPTIONAL :: ignore_missing
TYPE(um_file_type), POINTER :: ptr

INTEGER :: errorstatus
INTEGER :: i_elem
TYPE(um_file_list_type), POINTER :: um_file_list
CHARACTER(LEN=errormessagelength) :: message

um_file_list => select_list(handler)

ptr => um_file_list % head

DO i_elem = 1, um_file_list % n_file_units
  IF (ptr % UNIT == UNIT) THEN
    RETURN
  END IF
  ptr => ptr%next
END DO

! If it wasn't found in the loop we have a problem, unless the caller
! is expecting to catch this possibility
IF (PRESENT(ignore_missing)) THEN
  IF (ignore_missing) THEN
    ptr => NULL() 
    RETURN
  END IF
END IF

errorstatus = 1
!$OMP CRITICAL(internal_write)
WRITE(message,"(A,I3,A)") &
    "Unable to find specified unit (",UNIT,") in list ("//TRIM(handler)//")"
!$OMP END CRITICAL(internal_write)
CALL ereport("file_manager:get_file_by_unit", errorstatus, message)

END FUNCTION get_file_by_unit

! Assigns a pointer to a file object based on its id (if it has metadata)
!-------------------------------------------------------------------------------
FUNCTION get_file_by_id(id, handler, ignore_missing) RESULT(ptr)
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE 

CHARACTER(*) :: id
CHARACTER(*) :: handler
LOGICAL, OPTIONAL :: ignore_missing
TYPE(um_file_type), POINTER :: ptr

INTEGER :: errorstatus
INTEGER :: i_elem
TYPE(um_file_list_type), POINTER :: um_file_list
CHARACTER(LEN=errormessagelength) :: message

um_file_list => select_list(handler)

ptr => um_file_list % head

DO i_elem = 1, um_file_list % n_file_units
  IF (TRIM(ptr%id) == TRIM(id)) THEN
    RETURN
  END IF
  ptr => ptr%next
END DO

! If it wasn't found in the loop we have a problem, unless the caller
! is expecting to catch this possibility
IF (PRESENT(ignore_missing)) THEN
  IF (ignore_missing) THEN
    ptr => NULL() 
    RETURN
  END IF
END IF

errorstatus = 1
WRITE(message,"(A)") &
    "Unable to find specified id ("//TRIM(id)//") in list ("//TRIM(handler)//")"
CALL ereport("file_manager:get_file_by_id", errorstatus, message)

END FUNCTION get_file_by_id

! Provides access to the list print function
!-------------------------------------------------------------------------------
SUBROUTINE print_file_list()
IMPLICIT NONE 

CALL list_print(um_file_list_fortran)
CALL list_print(um_file_list_portio)
CALL list_print(um_file_list_netcdf)

END SUBROUTINE print_file_list

! Returns the next free unit number
!-------------------------------------------------------------------------------
SUBROUTINE get_file_unit(UNIT, list, unique)
IMPLICIT NONE 

INTEGER, INTENT(OUT)    :: UNIT
TYPE(um_file_list_type), INTENT(INOUT) :: list
LOGICAL, INTENT(IN) :: unique
INTEGER :: errorstatus

! For debugging purposes this flag will disable unit re-use
IF (unique) THEN
  UNIT = list % highest_unit_used + 1
ELSE
  UNIT = list % file_unit_min
  DO WHILE ( list % unit_in_use(UNIT) .AND. &
             UNIT <= list % file_unit_max )
    UNIT = UNIT + 1
  END DO
END IF

IF (UNIT == list % file_unit_max ) THEN
  errorstatus = 1
  CALL ereport("file_manager:get_file_unit", errorstatus, &
      "All units are in use, cannot provide unit")
ELSE
  list % unit_in_use(UNIT) = .TRUE.
  IF (UNIT > list % highest_unit_used) THEN
    list % highest_unit_used = UNIT
  END IF
END IF

END SUBROUTINE get_file_unit

! Frees the given unit number
!-------------------------------------------------------------------------------
SUBROUTINE free_file_unit(UNIT, list) 
IMPLICIT NONE 

INTEGER, INTENT(IN)    :: UNIT
TYPE(um_file_list_type), INTENT(INOUT) :: list

INTEGER :: errorstatus

IF (UNIT >= list % file_unit_min .AND. &
    UNIT <= list % file_unit_max) THEN
  list % unit_in_use(UNIT) = .FALSE.
ELSE
  errorstatus = 1
  CALL ereport("file_manager:free_file_unit", errorstatus, &
      "Provided with unit outside of allowed range")
END IF

END SUBROUTINE free_file_unit

END MODULE file_manager
