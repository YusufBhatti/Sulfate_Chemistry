! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE aspang_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASPANG_MOD'
CONTAINS

SUBROUTINE aspang(row_length, rows, bear_rot_NP,                        &
                  phi, delta_phi, delta_lambda, orog)

USE planet_constants_mod, ONLY: planet_radius
USE conversions_mod, ONLY: pi
USE solinc_data, ONLY: slope_aspect, slope_angle, horiz_limit
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Description:
!   Calculate mean slope aspect and angle in each gridbox and update
!   these variables in the global data module 'solinc_data'
!
! Method:
!   Method is taken from Gallant & Wilson, Computers & Geosciences
!   Vol. 22, No. 7, pp. 713-722, 1996
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.


! Subroutine arguments
INTEGER, INTENT(IN) ::                                                &
  row_length,                                                         &
  rows
!     Grid size

REAL, INTENT(IN) ::                                                   &
  bear_rot_NP(row_length,rows),                                       &
!     Bearing of 'pseudo' N pole (rads)
    phi(1-horiz_limit:rows+horiz_limit),                                &
!     Gridpoint latitude.
    delta_phi(1-horiz_limit:rows+horiz_limit),                          &
!     NS (y) grid spacing in radians.
    delta_lambda(1-horiz_limit:row_length+horiz_limit),                 &
!     EW (x) grid spacing in radians.
    orog(0:row_length+1,0:rows+1)
!     Mean gridpoint heights in m.

! Local variables
INTEGER ::                                                            &
  i, j
REAL ::                                                               &
  z2(0:row_length+1,0:rows+1),                                        &
  z4(0:row_length+1,0:rows+1),                                        &
  z6(0:row_length+1,0:rows+1),                                        &
  z8(0:row_length+1,0:rows+1),                                        &
!     Shifted orography grids
    zx(row_length,rows),                                                &
!     dz/dx slope in the x direction.
    zy(row_length,rows),                                                &
!     dz/dy slope in the y direction.
    xscale(row_length,rows),                                            &
!     Twice width of gridbox (x) in m.
    yscale(row_length,rows),                                            &
!     Twice width of gridbox (y) in m.
    eps
!     Tolerance to avoid division by zero.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASPANG'

! End of header


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

eps=EPSILON(eps)

! Allocate space for arrays in global data module:
IF (ALLOCATED(slope_aspect)) DEALLOCATE(slope_aspect)
IF (ALLOCATED(slope_angle))  DEALLOCATE(slope_angle)
ALLOCATE(slope_aspect (row_length,rows), &
         slope_angle  (row_length,rows) )

! Calculate the horizontal scales of the gridboxes (in metres):
DO j=1,rows
  DO i=1,row_length
    yscale(i,j)=SUM(delta_phi(j-1:j))*planet_radius
    xscale(i,j)=MAX(                                                  &
      SUM(delta_lambda(i-1:i))*planet_radius*COS(phi(j)), eps)
  END DO
END DO

! Calculate x and y slopes, using orography halos at boundaries:
z2 = EOSHIFT( orog, shift =  1, dim=1 )
z6 = EOSHIFT( orog, shift = -1, dim=1 )

z8 = EOSHIFT( orog, shift =  1, dim=2 )
z4 = EOSHIFT( orog, shift = -1, dim=2 )

zx=( z2(1:row_length,1:rows) - z6(1:row_length,1:rows) ) / xscale
zy=( z8(1:row_length,1:rows) - z4(1:row_length,1:rows) ) / yscale

slope_angle  = ATAN( (zx**2 + zy**2)**0.5 )
WHERE (zx == 0.0) zx=eps
slope_aspect = pi - ATAN(zy/zx) + (pi/2.0)*( zx/ABS(zx) )

! Add bearing of 'pseudo' N pole so aspects are relative to
! true North.
slope_aspect = MODULO(slope_aspect + bear_rot_NP, pi*2.0)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE aspang
END MODULE aspang_mod
