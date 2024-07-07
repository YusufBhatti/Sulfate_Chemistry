! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Grids
MODULE latlon_eq_rotation_mod

USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: rotate_eq_to_latlon, rotate_latlon_to_eq,                            &
          vector_eq_to_latlon, vector_latlon_to_eq, eq_latlon_vector_coeffs

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
CHARACTER(LEN=*),   PARAMETER :: ModuleName = 'LATLON_EQ_ROTATION_MOD'

CONTAINS

!------------------------------------------------------------------------------!

SUBROUTINE rotate_latlon_to_eq                                                 &
      (phi, lambda, phi_eq, lambda_eq, phi_pole, lambda_pole, points, cartesian)

USE f_shum_latlon_eq_grids_mod, ONLY: f_shum_latlon_to_eq

IMPLICIT NONE

INTEGER, INTENT(IN) :: points          ! Number of points to be processed

REAL, INTENT(IN)  :: phi(points)       ! Latitude
REAL, INTENT(IN)  :: lambda(points)    ! Longitude
REAL, INTENT(IN)  :: phi_pole          ! Latitude  of equatorial lat-lon pole
REAL, INTENT(IN)  :: lambda_pole       ! Longitude of equatorial lat-lon pole
REAL, INTENT(OUT) :: lambda_eq(points) ! Longitude in equatorial lat-lon coords
REAL, INTENT(OUT) :: phi_eq(points)    ! Latitude  in equatorial lat-lon coords

LOGICAL, INTENT(IN), OPTIONAL :: cartesian
LOGICAL :: l_cartesian
INTEGER :: status
CHARACTER(LEN=errormessagelength) :: message

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ROTATE_LATLON_TO_EQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check optional flag
l_cartesian = .FALSE.
IF (PRESENT(cartesian)) l_cartesian = cartesian

! Nothing to be done for cartesian grids
IF (l_cartesian) THEN
  lambda_eq = lambda
  phi_eq = phi
ELSE
  status = f_shum_latlon_to_eq(                                                &
                 phi, lambda, phi_eq, lambda_eq, phi_pole, lambda_pole, message)

  IF (status /= 0) THEN
    CALL ereport(RoutineName, status, message)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE rotate_latlon_to_eq

!------------------------------------------------------------------------------!

SUBROUTINE rotate_eq_to_latlon                                                 &
      (phi_eq, lambda_eq, phi, lambda, phi_pole, lambda_pole, points, cartesian)

USE f_shum_latlon_eq_grids_mod, ONLY: f_shum_eq_to_latlon

IMPLICIT NONE

INTEGER, INTENT(IN) :: points            ! Number of points to be processed

REAL, INTENT(IN)  :: lambda_eq(points) ! Longitude in equatorial lat-lon coords
REAL, INTENT(IN)  :: phi_eq(points)    ! Latitude in equatorial lat-lon coords
REAL, INTENT(IN)  :: phi_pole          ! Latitude of equatorial lat-lon pole
REAL, INTENT(IN)  :: lambda_pole       ! Longitude of equatorial lat-lon pole
REAL, INTENT(OUT) :: phi(points)       ! Latitude
REAL, INTENT(OUT) :: lambda(points)    ! Longitude (0 =< lon < 360)

LOGICAL, INTENT(IN), OPTIONAL :: cartesian
LOGICAL :: l_cartesian
INTEGER :: status
CHARACTER(LEN=errormessagelength) :: message

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='ROTATE_EQ_TO_LATLON'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check optional flag
l_cartesian = .FALSE.
IF (PRESENT(cartesian)) l_cartesian = cartesian

! Nothing to be done for cartesian grids
IF (l_cartesian) THEN
  phi = phi_eq
  lambda = lambda_eq
ELSE
  status = f_shum_eq_to_latlon(                                                &
                 phi_eq, lambda_eq, phi, lambda, phi_pole, lambda_pole, message)

  IF (status /= 0) THEN
    CALL ereport(RoutineName, status, message)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE rotate_eq_to_latlon

!------------------------------------------------------------------------------!

SUBROUTINE eq_latlon_vector_coeffs                                             &
   (coeff1, coeff2, lambda, lambda_eq, phi_pole, lambda_pole, points, cartesian)

USE f_shum_latlon_eq_grids_mod, ONLY: f_shum_latlon_eq_vector_coeff

IMPLICIT NONE

INTEGER, INTENT(IN) :: points          ! Number of points to be processed

REAL, INTENT(OUT) :: coeff1(points)    ! Coefficient of rotation no 1
REAL, INTENT(OUT) :: coeff2(points)    ! Coefficient of rotation no 2

REAL, INTENT(IN)  :: lambda(points)    ! Longitude
REAL, INTENT(IN)  :: lambda_eq(points) ! Longitude in equatorial lat-lon coords
REAL, INTENT(IN)  :: phi_pole          ! Latitude of equatorial lat-lon pole
REAL, INTENT(IN)  :: lambda_pole       ! Longitude of equatorial lat-lon pole

LOGICAL, INTENT(IN), OPTIONAL :: cartesian
LOGICAL :: l_cartesian
INTEGER :: status
CHARACTER(LEN=errormessagelength) :: message

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'EQ_LATLON_VECTOR_COEFFS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check optional flag
l_cartesian = .FALSE.
IF (PRESENT(cartesian)) l_cartesian = cartesian

IF (l_cartesian) THEN
  coeff1(:) = 1.0
  coeff2(:) = 0.0
ELSE
  status = f_shum_latlon_eq_vector_coeff(                                      &
              coeff1, coeff2, lambda, lambda_eq, phi_pole, lambda_pole, message)

  IF (status /= 0) THEN
    CALL ereport(RoutineName, status, message)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE eq_latlon_vector_coeffs

!------------------------------------------------------------------------------!

SUBROUTINE vector_eq_to_latlon                                                 &
                    (coeff1, coeff2, u_eq, v_eq, u, v, points, l_mdi, cartesian)

USE missing_data_mod, ONLY: rmdi  
USE f_shum_latlon_eq_grids_mod, ONLY: f_shum_eq_to_latlon_vector

IMPLICIT NONE

INTEGER, INTENT(IN) :: points       ! Number of points to be processed

REAL, INTENT(IN)  :: coeff1(points) ! Coefficient of rotation no 1
REAL, INTENT(IN)  :: coeff2(points) ! Coefficient of rotation no 2
REAL, INTENT(IN)  :: u_eq(points)   ! u component of wind on equatorial grid
REAL, INTENT(IN)  :: v_eq(points)   ! v component of wind on equatorial grid
REAL, INTENT(OUT) :: u(points)      ! u component of wind on lat-lon grid
REAL, INTENT(OUT) :: v(points)      ! v component of wind on lat-lon grid

LOGICAL, INTENT(IN) :: l_mdi        ! Shall we check for MDIs?
LOGICAL, INTENT(IN), OPTIONAL :: cartesian
LOGICAL :: l_cartesian
INTEGER :: status
CHARACTER(LEN=errormessagelength) :: message

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'VECTOR_EQ_TO_LATLON'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_cartesian = .FALSE.
IF (PRESENT(cartesian)) l_cartesian = cartesian

IF (l_cartesian) THEN
  u = u_eq
  v = v_eq
ELSE
  IF (l_mdi) THEN  ! supply optional rmdi
    status = f_shum_eq_to_latlon_vector(                                       &
                                coeff1, coeff2, u_eq, v_eq, u, v, message, rmdi)
  ELSE
    status = f_shum_eq_to_latlon_vector(                                       &
                                coeff1, coeff2, u_eq, v_eq, u, v, message)
  END IF

  IF (status /= 0) THEN
    CALL ereport(RoutineName, status, message)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE vector_eq_to_latlon

!------------------------------------------------------------------------------!

SUBROUTINE vector_latlon_to_eq                                                 &
                    (coeff1, coeff2, u, v, u_eq, v_eq, points, l_mdi, cartesian)

USE missing_data_mod, ONLY: rmdi
USE f_shum_latlon_eq_grids_mod, ONLY: f_shum_latlon_to_eq_vector

IMPLICIT NONE

INTEGER, INTENT(IN) :: points       !  Number of points to be processed

REAL, INTENT(IN)  :: coeff1(points) ! Coefficient of rotation no 1
REAL, INTENT(IN)  :: coeff2(points) ! Coefficient of rotation no 2
REAL, INTENT(OUT) :: u_eq(points)   ! u component of wind on equatorial grid
REAL, INTENT(OUT) :: v_eq(points)   ! v component of wind on equatorial grid
REAL, INTENT(IN)  :: u(points)      ! u component of wind on lat-lon grid
REAL, INTENT(IN)  :: v(points)      ! v component of wind on lat-lon grid

LOGICAL, INTENT(IN)           :: l_mdi        ! Shall we check for MDIs?
LOGICAL, INTENT(IN), OPTIONAL :: cartesian
LOGICAL :: l_cartesian
INTEGER :: status
CHARACTER(LEN=errormessagelength) :: message

REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='VECTOR_LATLON_TO_EQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_cartesian = .FALSE.
IF (PRESENT(cartesian)) l_cartesian = cartesian

IF (l_cartesian) THEN
  u_eq = u
  v_eq = v
ELSE
  IF (l_mdi) THEN  ! supply optional rmdi
    status = f_shum_latlon_to_eq_vector(                                       &
                                coeff1, coeff2, u, v, u_eq, v_eq, message, rmdi)
  ELSE
    status = f_shum_latlon_to_eq_vector(                                       &
                                coeff1, coeff2, u, v, u_eq, v_eq, message)
  END IF

  IF (status /= 0) THEN
    CALL ereport(RoutineName, status, message)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE vector_latlon_to_eq

!------------------------------------------------------------------------------!

END MODULE latlon_eq_rotation_mod
