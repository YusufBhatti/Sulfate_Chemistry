! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! ALLOCATE alloc_hires_data
MODULE alloc_hires_data_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOC_HIRES_DATA_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!   Allocate arrays to hold the high resolution data read in from the
!  fieldsfile/pp files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------

SUBROUTINE alloc_hires_data( row_length, rows, nlev)



USE hires_data_mod, ONLY:                                               &
    orog, landsea, precip, rain, snow, zh, lh, sh, pstar, tstar,        &
    theta, thetav, q, qcl, qcf, qrain, qgraup, iclass_col,              &
    p_theta_lev, rh, t,                                                 &
    u, v, w, density, dpdx, dpdy,                                       &
    dt1, dt2, dt4, dt9, dt12, dt30,                                     &
    dq4, dq9, dq12, dq30,                                               &
    dqcl4, dqcl9, dqcl12, dqcl30,                                       &
    dqcf4, dqcf3, dqcf12, dqcf30,                                       &
    drho, dqrain30, dqgr30

USE crmstyle_cntl_mod, ONLY:                                      &
  l_class_col, l_sect30   

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::  &
  row_length            & ! Number of columns
 ,rows                  & ! Number of rows
 ,nlev                    ! Number of levels

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k                  ! loop counters

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_HIRES_DATA'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Allocate 2d arrays for diagnostics

ALLOCATE(orog(row_length,rows) )
ALLOCATE(landsea(row_length,rows) )
ALLOCATE(precip(row_length,rows) )
ALLOCATE(rain(row_length,rows) )
ALLOCATE(snow(row_length,rows) )
ALLOCATE(zh(row_length,rows) )
ALLOCATE(lh(row_length,rows) )
ALLOCATE(sh(row_length,rows) )
ALLOCATE(pstar(row_length,rows) )
ALLOCATE(tstar(row_length,rows) )

ALLOCATE ( theta(row_length,rows,nlev) )
ALLOCATE ( thetav(row_length,rows,nlev) )
ALLOCATE ( q(row_length,rows,nlev) )
ALLOCATE ( qcl(row_length,rows,nlev) )
ALLOCATE ( qcf(row_length,rows,nlev) )
ALLOCATE ( qrain(row_length,rows,nlev) )
ALLOCATE ( qgraup(row_length,rows,nlev) )

ALLOCATE ( p_theta_lev(row_length,rows,nlev) )
ALLOCATE ( density(row_length,rows,nlev) )
ALLOCATE ( dpdx(row_length,rows,nlev) )
ALLOCATE ( dpdy(row_length,rows,nlev) )
ALLOCATE ( rh(row_length,rows,nlev) )
ALLOCATE ( t(row_length,rows,nlev) )

ALLOCATE ( u(row_length,rows,nlev) )
ALLOCATE ( v(row_length,rows,nlev) )
ALLOCATE ( w(row_length,rows,nlev) )

ALLOCATE ( dt1(row_length,rows,nlev) )
ALLOCATE ( dt2(row_length,rows,nlev) )
ALLOCATE ( dt4(row_length,rows,nlev) )
ALLOCATE ( dt9(row_length,rows,nlev) )
ALLOCATE ( dt12(row_length,rows,nlev) )

ALLOCATE ( dq4(row_length,rows,nlev) )
ALLOCATE ( dq9(row_length,rows,nlev) )
ALLOCATE ( dq12(row_length,rows,nlev) )

ALLOCATE ( dqcl4(row_length,rows,nlev) )
ALLOCATE ( dqcl9(row_length,rows,nlev) )
ALLOCATE ( dqcl12(row_length,rows,nlev) )

ALLOCATE ( dqcf4(row_length,rows,nlev) )
ALLOCATE ( dqcf3(row_length,rows,nlev) )
ALLOCATE ( dqcf12(row_length,rows,nlev) )

! Section 30 diagnostics
IF (l_sect30) THEN
  ALLOCATE ( dt30(row_length,rows,nlev) )
  ALLOCATE ( dq30(row_length,rows,nlev) )
  ALLOCATE ( dqcf30(row_length,rows,nlev) )
  ALLOCATE ( dqcl30(row_length,rows,nlev) )
  ALLOCATE ( drho(row_length,rows,nlev) )
  ALLOCATE ( dqrain30(row_length,rows,nlev) )
  ALLOCATE ( dqgr30(row_length,rows,nlev) )
END IF

! column classification

IF (l_class_col) THEN
  ALLOCATE ( iclass_col(row_length,rows) )
END IF

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE alloc_hires_data
END MODULE alloc_hires_data_mod
