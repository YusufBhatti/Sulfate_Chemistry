! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine profile_increment
!
MODULE profile_increment_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROFILE_INCREMENT_MOD'

CONTAINS

SUBROUTINE profile_increment(profile, ftype, fdims, field_inc)

! Purpose: Apply increment profile to a field
!
! Method: 
!   (1) Increment (tendency) forcing data are interpolated in time to give a
!       forcing profile valid at the current model time. Note that input data
!       are tendencies (SI units per second), while this routine returns an
!       increment scaled by the model timestep.
!   (2) The forcing profile is interpolated onto the model grid. The vertical
!       coordinate used to specify the profile may be non-dimensional height
!       (ie. eta), height above surface or reference pressure.
!

USE planet_constants_mod,     ONLY: planet_radius
USE timestep_mod,             ONLY: timestep
USE field_types,              ONLY: fld_type_u, fld_type_v, fld_type_w
USE atm_fields_bounds_mod,    ONLY: array_dims

USE level_heights_mod,        ONLY: eta_theta_levels, eta_rho_levels,          &
                                    z_ref_theta, z_ref_rho,                    &
                                    r_theta_levels, r_at_u, r_at_v

USE idealise_run_mod,         ONLY: zero_orography
USE profiles_mod,             ONLY: varying_profile, eta_coord, z_coord,       &
                                    p_coord, p_prof_theta, p_prof_rho
USE prof_temporal_interp_mod,                                                  &
                              ONLY: prof_temporal_interp

USE Ereport_mod,              ONLY: Ereport
USE errormessagelength_mod,   ONLY: errormessagelength
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

! Input
! ---------------------------------------------------------------------
! - profile information 
TYPE(varying_profile),      INTENT(IN) :: profile
! - array bounds for field to be incremented
TYPE(array_dims),           INTENT(IN) :: fdims
! - grid-stagger of the field (p, u, v or w)
INTEGER,                    INTENT(IN) :: ftype

! Output
! ---------------------------------------------------------------------
! - increment to the field
REAL, INTENT(OUT)             :: field_inc(fdims%i_start:fdims%i_end,          &
                                           fdims%j_start:fdims%j_end,          &
                                           fdims%k_start:fdims%k_end)
! Local
! ---------------------------------------------------------------------
INTEGER                       :: i, j, k
INTEGER                       :: interp_pt
INTEGER                       :: n_heights
INTEGER                       :: coord_dims
REAL                          :: interp_wt
REAL                          :: height(profile % n_heights)
REAL                          :: prof_now(profile % n_heights)
REAL, POINTER                 :: levels_1D(:)
REAL, POINTER                 :: levels_3D(:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus
CHARACTER(LEN=errormessagelength)                                              &
                              :: Cmessage
CHARACTER(LEN=*), PARAMETER   :: RoutineName='PROFILE_INCREMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Interpolate increment profile to current time
CALL prof_temporal_interp(profile, prof_now)

! Local copy of data levels
n_heights = profile % n_heights
height(:) = profile % height(:)

! If forcing data is specified with 1D vertical coordinate (eta):
SELECT CASE(profile % coord_type)

CASE(eta_coord)

  coord_dims = 1

  ! Select grid appropriate to the field
  SELECT CASE(ftype)
  CASE(fld_type_u, fld_type_v)
    levels_1D => eta_rho_levels
  CASE(fld_type_w)
    levels_1D => eta_theta_levels
  END SELECT

CASE(p_coord)

  coord_dims = 1

  ! Select grid appropriate to the field
  SELECT CASE(ftype)
  CASE(fld_type_u, fld_type_v)
    levels_1D => p_prof_rho
  CASE(fld_type_w)
    levels_1D => p_prof_theta
  END SELECT

CASE(z_coord)

  IF (zero_orography) THEN

    coord_dims = 1

    SELECT CASE(ftype)
    CASE(fld_type_u, fld_type_v)
      levels_1D => z_ref_rho
    CASE(fld_type_w)
      levels_1D => z_ref_theta
    END SELECT

  ELSE

    coord_dims = 3

    ! Select grid appropriate to the field
    SELECT CASE(ftype)
    CASE(fld_type_u)
      levels_3D => r_at_u
    CASE(fld_type_v)
      levels_3D => r_at_v
    CASE(fld_type_w)
      levels_3D => r_theta_levels
    END SELECT

    levels_3D = levels_3D - planet_radius

  END IF

CASE DEFAULT
  ErrorStatus = 1
  WRITE(Cmessage, '(''Unknown profile % coord_type = '', I0)')                 &
                  profile % coord_type
  CALL Ereport(RoutineName, ErrorStatus, Cmessage)

END SELECT

IF (coord_dims == 1) THEN

  DO k = fdims%k_start, fdims%k_end

! Find data point immediately below this grid level and calculate
! linear interpolation weight
    IF (levels_1D(k) <= height(1)) THEN
      interp_wt = 0.0
      interp_pt = 1
    ELSE IF (levels_1D(k) >= height(n_heights)) THEN
      interp_wt = 1.0
      interp_pt = n_heights - 1
    ELSE
      DO interp_pt = 1, n_heights - 1
        IF (height(interp_pt+1) > levels_1D(k)) EXIT
      END DO
      interp_wt = (levels_1D(k)        - height(interp_pt)) /                  &
                  (height(interp_pt+1) - height(interp_pt))
    END IF

    DO j = fdims%j_start, fdims%j_end
      DO i = fdims%i_start, fdims%i_end
        field_inc(i,j,k) = timestep *                                          &
                           ( (1.0 - interp_wt) * prof_now(interp_pt) +         &
                                     interp_wt * prof_now(interp_pt+1) )
      END DO
    END DO

  END DO

ELSE

  DO k = fdims%k_start, fdims%k_end
    DO j = fdims%j_start, fdims%j_end
      DO i = fdims%i_start, fdims%i_end
        IF (levels_3D(i,j,k) <= height(1)) THEN
          field_inc(i,j,k) = timestep * prof_now(1)
        ELSE IF (levels_3D(i,j,k) >= height(n_heights)) THEN
          field_inc(i,j,k) = timestep * prof_now(n_heights)
        ELSE
          DO interp_pt = 1, n_heights - 1
            IF (height(interp_pt+1) > levels_3D(i,j,k)) EXIT
          END DO
          interp_wt = (levels_3D(i,j,k)    - height(interp_pt)) /              &
                      (height(interp_pt+1) - height(interp_pt))
        END IF
        field_inc(i,j,k) = timestep *                                          &
                           ( (1.0 - interp_wt) * prof_now(interp_pt) +         &
                                     interp_wt * prof_now(interp_pt+1) )
      END DO
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE profile_increment

END MODULE profile_increment_mod
