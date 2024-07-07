! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE autodetect_file_mod

USE um_types,           ONLY: integer32
USE io,                 ONLY: setpos, buffin, file_open, file_close
USE io_constants,       ONLY: ioNameProvided, ioOpenReadOnly
USE filenamelength_mod, ONLY: filenamelength
USE file_manager,       ONLY: assign_file_unit, release_file_unit
USE fieldsfile_constants_mod, ONLY: fixed_header_1
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! Description:
!   Subroutine to try and detect the file format of the specified filename.
!   Currently will detect a fieldsfile or a GRIB file. The latter case, 
!   is present only as an example, it is not supported in CreateBC and, 
!   together with unknown file formats will cause the code to fail.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

IMPLICIT NONE

INTEGER, PARAMETER :: start_of_file = 0

INTEGER, PARAMETER :: filetype_unknown = 0
INTEGER, PARAMETER :: filetype_fieldsfile = 1
INTEGER, PARAMETER :: filetype_grib2 = 2

! denotes GRIB format (backwards due to the way the header is converted)
CHARACTER(LEN=4), PARAMETER :: grib_identifier = 'BIRG'

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'AUTODETECT_FILE_MOD'

CONTAINS

INTEGER FUNCTION detect_file_type(input_filename)
! The value of this function will be used to determine the correct subclass of
! datafile_type to instanciate as the list of input files.

IMPLICIT NONE
  
CHARACTER(LEN=filenamelength) :: input_filename

INTEGER :: single_64bit_value(1)
INTEGER(KIND=integer32) :: four_32bit_values(4)
CHARACTER(LEN=4) :: four_character_string
INTEGER :: unit_num
INTEGER :: single_length = 1
CHARACTER(LEN=*), PARAMETER :: routinename = 'DETECT_FILE_TYPE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Open the file
CALL assign_file_unit(input_filename, unit_num, handler="portio")
CALL file_open(unit_num, input_filename, filenamelength,                   &
               read_write=ioOpenReadOnly, name_in_environ=ioNameProvided)

! First, attempt to ascertain if this is a fieldsfile
CALL setpos(unit_num, start_of_file)
CALL buffin(unit_num, single_64bit_value, single_length)

! Second, check for GRIB2  
CALL setpos(unit_num, start_of_file)
CALL buffin(unit_num, four_32bit_values, 4)
four_character_string = TRANSFER(four_32bit_values, four_character_string)
  
! Further file formats can be added here, by looking for whatever marker
! indicates a given file type.

! Analyse the results

IF (single_64bit_value(1) == fixed_header_1) THEN
  ! This is a fieldsfile
  detect_file_type = filetype_fieldsfile
ELSE IF (four_character_string == grib_identifier) THEN
  ! this is a GRIB2 file
  detect_file_type = filetype_grib2
ELSE
  ! unknown file type
  detect_file_type = filetype_unknown
END IF

CALL file_close(unit_num, input_filename)
CALL release_file_unit(unit_num, handler="portio")
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION detect_file_type

END MODULE autodetect_file_mod

