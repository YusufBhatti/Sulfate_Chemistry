! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE profile_uv_geo_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROFILE_UV_GEO_MOD'

CONTAINS

SUBROUTINE profile_uv_geo(profile, fdims, ftype, field_ref)

! Purpose: To calculate the geostrophic wind profile prior to the forcing
!
! Method:  Work out the geostrophic wind profile for this time step.
!

USE atm_fields_bounds_mod,        ONLY: array_dims
USE field_types,                  ONLY: fld_type_u, fld_type_v, fld_type_w

USE level_heights_mod,            ONLY: eta_theta_levels, eta_rho_levels,      &
                                        z_ref_theta, z_ref_rho
USE profiles_mod,                 ONLY: varying_profile,                       &
                                        eta_coord, p_coord, z_coord,           &
                                        p_prof_theta, p_prof_rho
USE prof_temporal_interp_mod,     ONLY: prof_temporal_interp
USE prof_interp_mod,              ONLY: prof_interp, interp_linear,            &
                                        extrap_constant
USE idealise_run_mod,             ONLY: zero_orography

USE ereport_mod,                  ONLY: ereport
USE errormessagelength_mod,       ONLY: errormessagelength
USE yomhook,                      ONLY: lhook, dr_hook
USE parkind1,                     ONLY: jprb, jpim
  
IMPLICIT NONE

! ---------------------------------------------------------------------
TYPE(varying_profile), INTENT(IN) :: profile   ! U or V forcing profiles
TYPE(array_dims),      INTENT(IN) :: fdims     ! array bounds for data

INTEGER,               INTENT(IN) :: ftype

REAL,                INTENT(OUT)  :: field_ref(fdims%k_start:fdims%k_end)
                                      ! U or V profile for this time
! ---------------------------------------------------------------------
! Local
! ---------------------------------------------------------------------
INTEGER                       :: i, j, k
REAL                          :: prof_now(profile % n_heights)
REAL, POINTER                 :: levels(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus
CHARACTER(LEN=errormessagelength)                                              &
                              :: Cmessage
CHARACTER(LEN=*), PARAMETER   :: RoutineName='PROFILE_UV_GEO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that orography is flat
IF (.NOT. zero_orography) THEN
  ErrorStatus = 1
  Cmessage = 'Relaxation not coded for non-zero orography'
  CALL ereport(RoutineName, ErrorStatus, Cmessage)
END IF

! Interpolate reference profile to current time 
CALL prof_temporal_interp(profile, prof_now)

SELECT CASE(profile % coord_type)

CASE(eta_coord)

  ! Select grid appropriate to the field
  SELECT CASE(ftype)
  CASE(fld_type_u, fld_type_v)
    levels => eta_rho_levels
  CASE(fld_type_w)
    levels => eta_theta_levels
  END SELECT

CASE(p_coord)

  ! Select grid appropriate to the field
  SELECT CASE(ftype)
  CASE(fld_type_u, fld_type_v)
    levels => p_prof_rho
  CASE(fld_type_w)
    levels => p_prof_theta
  END SELECT

CASE(z_coord)

    ! Select grid appropriate to the field
    SELECT CASE(ftype)
    CASE(fld_type_u)
      levels => z_ref_rho
    CASE(fld_type_v)
      levels => z_ref_rho
    CASE(fld_type_w)
      levels => z_ref_theta
    END SELECT

CASE DEFAULT

  ErrorStatus = 1
  WRITE(Cmessage, '(''Unknown profile % coord_type = '', I0)')                 &
                  profile % coord_type
  CALL Ereport(RoutineName, ErrorStatus, Cmessage)

END SELECT


! - and onto the reference grid
! NB. This could be made to work with orography - see profile_increment.F90
CALL prof_interp(interp_linear, extrap_constant, profile % n_heights,          &
                 profile % height, prof_now, fdims%k_end-fdims%k_start+1,      &
                 levels, field_ref)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE profile_uv_geo

END MODULE profile_uv_geo_mod
