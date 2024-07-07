! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt 
! This file belongs in section: Top Level


MODULE river_routing_sizes_mod


USE nlsizes_namelist_mod, ONLY: aocpl_row_length, aocpl_p_rows
USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, newline


IMPLICIT NONE

SAVE

! Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids 

REAL, ALLOCATABLE :: xpa(:)  ! Atmosphere TP longitude coordinates
REAL, ALLOCATABLE :: xua(:)  ! Atmosphere U longitude coordinates
REAL, ALLOCATABLE :: xva(:)  ! Atmosphere V longitude coordinates
REAL, ALLOCATABLE :: ypa(:)  ! Atmosphere TP latitude coordinates
REAL, ALLOCATABLE :: yua(:)  ! Atmosphere U latitude coordinates
REAL, ALLOCATABLE :: yva(:)  ! Atmosphere V latitude coordinates

CHARACTER(LEN=*), PARAMETER, PRIVATE :: modulename = 'RIVER_ROUTING_SIZES_MOD'

CONTAINS

SUBROUTINE allocate_river_routing_sizes()

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: routinename = 'ALLOCATE_RIVER_ROUTING_SIZES'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

icode=0


IF (ALLOCATED (xpa)) DEALLOCATE(xpa)
IF (ALLOCATED (xua)) DEALLOCATE(xua)
IF (ALLOCATED (xva)) DEALLOCATE(xva)
IF (ALLOCATED (ypa)) DEALLOCATE(ypa)
IF (ALLOCATED (yua)) DEALLOCATE(yua)
IF (ALLOCATED (yva)) DEALLOCATE(yva)


ALLOCATE( xpa(aocpl_row_length+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [xpa]'//          newline//         &
            TRIM(umMessage) )

ALLOCATE( xua(0:aocpl_row_length), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [xua]'//          newline//         &
            TRIM(umMessage) )

ALLOCATE( xva(aocpl_row_length+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [xva]'//         newline//          &
            TRIM(umMessage) )

ALLOCATE( ypa(aocpl_p_rows), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [ypa]'//          newline//         &
            TRIM(umMessage) )

ALLOCATE( yua(aocpl_p_rows), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [yua]'//          newline//         &
            TRIM(umMessage) )

ALLOCATE( yva(0:aocpl_p_rows), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0)                                                                &
CALL ereport( modulename//routinename, icode                                   &
            , 'Failure in allocating array [yva]'//         newline//          &
            TRIM(umMessage) )


IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_river_routing_sizes

END MODULE river_routing_sizes_mod
