! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE aspang_ancil_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASPANG_ANCIL_MOD'
CONTAINS

SUBROUTINE aspang_ancil(row_length, rows, land_points, land_sea_mask, &
                        grad_x, grad_y, bear_rot_NP)

USE conversions_mod, ONLY: pi

USE solinc_data, ONLY: slope_aspect, slope_angle
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Description:
!   Calculate mean slope aspect and angle in each gridbox and update
!   these variables in the global data module 'solinc_data'
!
! Method:
!   Uses X & Y gradients from ancillary fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 95
!
! Declarations:
!
! Global variables (#include statements etc):
!   Data from module 'solinc_data' available by USE association.

! Subroutine arguments

INTEGER, INTENT(IN) :: &
     row_length, rows, &      ! grid size
     land_points              ! number of land points

LOGICAL, INTENT(IN) :: land_sea_mask(row_length,rows)  ! land-sea mask

REAL, INTENT(IN) :: &
     grad_x(land_points),   &         ! orographic X-gradient on land points
     grad_y(land_points)              ! orographic Y-gradient on land points

REAL, INTENT(IN) :: &
     bear_rot_NP(row_length,rows)     ! bearing of 'pseudo' N pole (rads)

! Local variables

REAL :: work (land_points)            ! work array on land points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASPANG_ANCIL'

!- End of header

! Allocate space for arrays in global data module:

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ALLOCATED(slope_angle))  DEALLOCATE(slope_angle)
IF (ALLOCATED(slope_aspect)) DEALLOCATE(slope_aspect)
ALLOCATE(slope_angle (row_length,rows))
ALLOCATE(slope_aspect(row_length,rows))
slope_angle=0.0
slope_aspect=0.0

! Find slope angles and aspects from x & y gradients:

work = ATAN( (grad_x**2 + grad_y**2)**0.5 )
slope_angle = UNPACK(work, land_sea_mask, slope_angle)

work = grad_x
WHERE (grad_x == 0.0) work = EPSILON(grad_x)
work = pi - ATAN(grad_y/work) + SIGN(pi/2.0,work)
slope_aspect = UNPACK(work, land_sea_mask, slope_aspect)

! Add bearing of 'pseudo' N pole so aspects are relative to
! true North:

slope_aspect = MODULO(slope_aspect + bear_rot_NP, pi*2.0)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE aspang_ancil
END MODULE aspang_ancil_mod
