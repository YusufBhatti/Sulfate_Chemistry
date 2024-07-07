! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level


MODULE land_soil_dimensions_mod

USE nlsizes_namelist_mod, ONLY: land_field
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, newline

IMPLICIT NONE

SAVE


INTEGER :: land_points      ! No. of land points  (can be 0)
INTEGER :: land_ice_points  ! Number of land ice points
INTEGER :: soil_points      ! Number of soil points

INTEGER, ALLOCATABLE :: land_index    (:) ! Array of land points.
INTEGER, ALLOCATABLE :: land_ice_index(:) ! Array of land ice points.
INTEGER, ALLOCATABLE :: soil_index    (:) ! Array of soil points.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LAND_SOIL_DIMENSIONS_MOD'

CONTAINS

SUBROUTINE allocate_land_soil_arrays()

! Purpose: Allocate the arrays for land points, land ice points
!          and soil points

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_LAND_SOIL_ARRAYS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0

CALL deallocate_land_soil_arrays()

ALLOCATE( land_index (land_field), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
   CALL ereport( modulename//":"//routinename, icode                   &
                 , 'Failure in allocating array [land_index]'//        &
                 newline// TRIM(umMessage) )
ENDIF

ALLOCATE( land_ice_index (land_field), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
   CALL ereport( modulename//":"//routinename, icode                   &
                 , 'Failure in allocating array [land_ice_index]'//    &
                 newline// TRIM(umMessage) )
ENDIF

ALLOCATE( soil_index (land_field), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
   CALL ereport( modulename//":"//routinename, icode                   &
                 , 'Failure in allocating array [soil_index]'//        &
                 newline// TRIM(umMessage) )
ENDIF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_land_soil_arrays


SUBROUTINE deallocate_land_soil_arrays()

IMPLICIT NONE

IF (ALLOCATED( soil_index    )) DEALLOCATE(soil_index)
IF (ALLOCATED( land_ice_index)) DEALLOCATE(land_ice_index)
IF (ALLOCATED( land_index    )) DEALLOCATE(land_index)

END SUBROUTINE deallocate_land_soil_arrays

END MODULE land_soil_dimensions_mod


