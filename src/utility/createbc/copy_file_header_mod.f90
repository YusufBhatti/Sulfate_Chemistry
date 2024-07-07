! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE copy_file_header_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'COPY_FILE_HEADER_MOD'

CONTAINS

SUBROUTINE copy_file_header(source, target_file)

USE fieldsfile_mod,         ONLY: fieldsfile_type
USE lbcfile_mod,            ONLY: lbcfile_type
USE datafile_mod,           ONLY: datafile_type
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
  
IMPLICIT NONE
  
CLASS(datafile_type), INTENT(IN) :: source
CLASS(datafile_type), INTENT(INOUT) :: target_file

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename='COPY_FILE_HEADER'
INTEGER :: icode
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT TYPE(source)
CLASS IS (fieldsfile_type)
      
  SELECT TYPE(target_file)
  CLASS IS (fieldsfile_type)
  
    CALL target_file%copy_common_metadata(source)

  CLASS IS (lbcfile_type)

    CALL target_file%copy_common_metadata(source)

    CLASS DEFAULT
    cmessage = 'Unknown target class, unable to copy header'
    icode = 10
    CALL ereport(routinename, icode, cmessage)
  END SELECT
    
  CLASS DEFAULT
  cmessage = 'Unknown source class, unable to copy header'
  icode = 9
  CALL ereport(routinename, icode, cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE copy_file_header

END MODULE copy_file_header_mod
