! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_dx_out_grid(fld_type,x,y,z,                               &
                          row_length,rows,levels,                       &
                          halo_x, halo_y, offx, offy, mype,             &
                          planet_radius)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod, ONLY: ereport
USE file_manager, ONLY: assign_file_unit, release_file_unit

IMPLICIT NONE
!
! Description:
!
! Output the 3D grid for use by openDX
!
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(8), INTENT(IN) :: row_length, rows, levels
INTEGER(8), INTENT(IN) :: halo_x, halo_y, offx, offy
INTEGER(8), INTENT(IN) :: fld_type, mype

REAL(8),    INTENT(IN) :: z(1-halo_x:row_length+halo_x,                 &
                         1-halo_y:rows+halo_y,levels)
REAL(8),    INTENT(IN) :: x(1-halo_x:row_length+halo_x)
REAL(8),    INTENT(IN) :: y(1-halo_y:rows+halo_y)
REAL(8),    INTENT(IN) :: planet_radius

REAL(4), ALLOCATABLE   :: grid(:)
INTEGER                :: npx, npy, npz
INTEGER                :: i, j, k, l
INTEGER                :: eg_unit
CHARACTER(LEN=15)      :: cgrd
CHARACTER(LEN=4)       :: cpe

INTEGER                :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DX_OUT_GRID'



! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


npx   = row_length
npy   = rows
npz   = levels

WRITE(cpe,"('-',i3.3)")  mype

SELECT CASE(fld_type)
CASE (1)
  cgrd = 'u_grid'//cpe//'.dat'
CASE (2)
  cgrd = 'v_grid'//cpe//'.dat'
CASE (3)
  cgrd = 'w_grid'//cpe//'.dat'
CASE (4)
  cgrd = 'p_grid'//cpe//'.dat'
CASE DEFAULT
  errcode = 1
  CALL Ereport('eg_dx_out_grid',errcode,'Unrecognized field type')
END SELECT

ALLOCATE(grid(3*npx*npy*npz))

l = 1
DO k = 1, npz
  DO j = 1, rows
    DO i = 1, row_length
      grid(l)   = x(i)
      grid(l+1) = y(j)
      grid(l+2) = z(i,j,k) - planet_radius
      l = l + 3
    END DO
  END DO
END DO

CALL assign_file_unit(cgrd, eg_unit, handler="fortran")
OPEN(UNIT=eg_unit, FILE=cgrd, STATUS="unknown", FORM="unformatted")
WRITE(UNIT=eg_unit) grid
CLOSE(UNIT=eg_unit)
CALL release_file_unit(eg_unit, handler="fortran")

DEALLOCATE(grid)


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_dx_out_grid
