! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE eg_dxout_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_DXOUT_MOD'

CONTAINS

SUBROUTINE eg_dxout(fnam,fno,fld,                                        &
                    row_length,rows,levels,halo_x, halo_y,               &
                    fld_type,mype,iter)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE um_parvars, ONLY: nproc_x,nproc_y
USE ereport_mod, ONLY: ereport
USE file_manager, ONLY: assign_file_unit, release_file_unit

IMPLICIT NONE
!
! Description:
!
! This subroutine produces output in a form suitable for use with
! openDX. The header file is in fnam-xxx.dx while the raw data
! file is in fnam-xxx.dat. Here xxx is the integer fno
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics 
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER, INTENT(IN) :: row_length, rows, levels, halo_x, halo_y
INTEGER, INTENT(IN) :: fno, fld_type
REAL,    INTENT(IN) :: fld(1-halo_x:row_length+halo_x,                   &
                           1-halo_y:rows+halo_y,levels)

CHARACTER(LEN=*), INTENT(IN) :: fnam

INTEGER             :: eg_unit
INTEGER             :: npx, npy, npz, npnts, mype
CHARACTER(LEN=30)   :: cfile
CHARACTER(LEN=14)   :: cgrd
CHARACTER(LEN=7)    :: cvar
CHARACTER(LEN=4)    :: cit
CHARACTER(LEN=4)    :: cpe

INTEGER             :: err_code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DXOUT'

INTEGER, OPTIONAL :: iter




! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


npnts = LEN_TRIM(fnam)
IF ( npnts > 12 ) THEN
  err_code = 1
  CALL Ereport('eg_dxout',err_code,'Filename prefix too long!')
END IF

npx   = row_length
npy   = rows
npz   = levels

npnts = npx*npy*npz

! set name of grid data file based on fld_type

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
  err_code = 1
  CALL Ereport('eg_dxout', err_code,'Unrecognized field type')
END SELECT

! Convert the integer fno into a character string

WRITE(cvar,"('-',i6.6)") fno

IF (PRESENT(iter)) THEN
  WRITE(cit,"('-',i3.3)") iter
  cfile = TRIM(fnam)//cvar//cit//cpe
ELSE
  cfile = TRIM(fnam)//cvar//cpe
END IF


! Create header file
CALL assign_file_unit(TRIM(cfile)//".dx", eg_unit, handler="fortran")
OPEN(UNIT=eg_unit, FILE = TRIM(cfile)//".dx", STATUS="unknown")
WRITE(eg_unit,FMT='(A,I5,A,I5)') '# ncpusx: ', nproc_x,' ncpusy: ',nproc_y
WRITE(eg_unit,"('object 1 class array type float rank 1 shape 3 ')")
WRITE(eg_unit,"('items ',i7,' msb binary ')") npnts
WRITE(eg_unit,"('DATA FILE ',a13,' ,4')") cgrd

WRITE(eg_unit,*) 'object 2 class gridconnections counts ', npz, npy, npx
WRITE(eg_unit,"('object 3 class array type double rank 0 items')")
WRITE(eg_unit,"(i8,' msb binary')") npnts
WRITE(eg_unit,*) 'DATA FILE ',TRIM(cfile)//'.dat',',4'

WRITE(eg_unit,*) 'attribute "dep" string "positions"'
WRITE(eg_unit,*) 'object "irreg positions regular connections" class field'
WRITE(eg_unit,*) 'component "positions" value 1'
WRITE(eg_unit,*) 'component "connections" value 2'
WRITE(eg_unit,*) 'component "data" value 3'
WRITE(eg_unit,*) 'end'

CLOSE(eg_unit)

! Now output unformatted binary file

OPEN(UNIT=eg_unit, FILE=TRIM(cfile)//".dat", STATUS="unknown", FORM="unformatted")
WRITE(eg_unit) fld(1:npx,1:npy,:)
CLOSE(eg_unit)
CALL release_file_unit(eg_unit, handler="fortran")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_dxout

END MODULE
