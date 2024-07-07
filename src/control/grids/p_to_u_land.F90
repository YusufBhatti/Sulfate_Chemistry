! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Subroutine P_TO_U_LAND for calculating land variables held at
!   p points at u points.
!
!   This routine does interior points of array not halos,
!   but requires halo information to be set.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! global code has E-W wrap around
! LAM code has most east u point set to zero
!
! Code Description:
! Language: FORTRAN 95

MODULE p_to_u_land_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'P_TO_U_LAND_MOD'

CONTAINS

SUBROUTINE p_to_u_land                                                         &
  ( array_on_p_points, fland_on_p_points, ini_start, ini_end                   &
  , inj_start, inj_end, outi_start, outi_end, outj_start, outj_end             &
  , outk_start, outk_end, array_on_u_points )

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: ini_start, ini_end
INTEGER, INTENT(IN) :: inj_start, inj_end
INTEGER, INTENT(IN) :: outi_start, outi_end
INTEGER, INTENT(IN) :: outj_start, outj_end
INTEGER, INTENT(IN) :: outk_start, outk_end

REAL, INTENT(IN)  :: array_on_p_points( ini_start:ini_end                      &
                                     , inj_start:inj_end                       &
                                     , outk_start:outk_end )

REAL, INTENT(IN)  :: fland_on_p_points( ini_start:ini_end                      &
                                      , inj_start:inj_end )
                                      ! Land fraction on p points.

REAL, INTENT(OUT) :: array_on_u_points( outi_start:outi_end                    &
                                      , outj_start:outj_end                    &
                                      , outk_start:outk_end )

! Local variables
INTEGER :: i, j, k
REAL    :: fland_total    ! Total land fraction on u grid point

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'P_TO_U_LAND'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO k=outk_start, outk_end
  DO j=outj_start, outj_end

    DO i=outi_start, outi_end
      fland_total = fland_on_p_points(i,j)                                   &
                  + fland_on_p_points(i+1,j)

      IF (fland_total >  0.0) THEN
        array_on_u_points(i,j,k) = 1.0/fland_total*                          &
          ( fland_on_p_points(i,j)   * array_on_p_points(i,j,k)              &
          + fland_on_p_points(i+1,j) * array_on_p_points(i+1,j,k) )
      ELSE
        array_on_u_points(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE p_to_u_land

END MODULE p_to_u_land_mod
