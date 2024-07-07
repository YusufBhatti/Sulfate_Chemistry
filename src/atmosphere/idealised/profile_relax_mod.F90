! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine profile_relax
!
MODULE profile_relax_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROFILE_RELAX_MOD'

CONTAINS

SUBROUTINE profile_relax(profile, fdims, fdims_h, fdims_inc, ftype, field,   &
                         field_inc)

! Purpose: To apply Newtonian relaxation to a field
!
! Method:  Calculates the domain average of the field for each level.
!          Calculates the difference between the domain avg profile and
!          the 1D forcing profile.
!          Applies Newtonian relaxation of this increment using a
!          defined relaxation timescale.
!

USE timestep_mod,                 ONLY: timestep
USE atm_fields_bounds_mod,        ONLY: array_dims
USE field_types,                  ONLY: fld_type_u, fld_type_v, fld_type_w
USE nlsizes_namelist_mod,         ONLY: global_row_length, global_rows
USE global_2d_sums_mod,           ONLY: global_2d_sums
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

! Input
! ---------------------------------------------------------------------
TYPE(varying_profile), INTENT(IN) :: profile
TYPE(array_dims),      INTENT(IN) :: fdims     ! array bounds for data
TYPE(array_dims),      INTENT(IN) :: fdims_h   ! bounds including halos
TYPE(array_dims),      INTENT(IN) :: fdims_inc ! bounds of field_inc

INTEGER,               INTENT(IN) :: ftype
! Field to relax
REAL,                  INTENT(IN) :: field(fdims_h%i_start:fdims_h%i_end,      &
                                           fdims_h%j_start:fdims_h%j_end,      &
                                           fdims_h%k_start:fdims_h%k_end)

! Changed on output
! ---------------------------------------------------------------------
! Field increment
REAL, INTENT(INOUT)           :: field_inc(fdims_inc%i_start:fdims_inc%i_end,  &
                                           fdims_inc%j_start:fdims_inc%j_end,  &
                                           fdims_inc%k_start:fdims_inc%k_end)
! Local
! ---------------------------------------------------------------------
INTEGER                       :: i, j, k
REAL                          :: relax_coeff, relax_increment
REAL                          :: prof_now(profile % n_heights)
REAL, POINTER                 :: levels(:)
REAL                          :: field_ref(fdims%k_start:fdims%k_end)
REAL                          :: field_2D_ave(fdims%k_start:fdims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus
CHARACTER(LEN=errormessagelength)                                              &
                              :: Cmessage
CHARACTER(LEN=*), PARAMETER   :: RoutineName='PROFILE_RELAX'

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

! Sum the field on each grid-level
CALL global_2d_sums(field, fdims%i_len, fdims%j_len,                           &
                           fdims_h%halo_i, fdims_h%halo_j,                     &
                           fdims%k_len, field_2d_ave)

relax_coeff = timestep / profile % timescale

DO k = fdims%k_start, fdims%k_end
  ! Average of field on each grid-level
  field_2d_ave(k) = field_2d_ave(k) / (global_row_length*global_rows)
  ! Only apply relaxation within specified height range
  ! - this need to be changed to work with orography
  IF (levels(k) < profile % height(1) .OR.                                     &
      levels(k) > profile % height(profile % n_heights)) CYCLE
  relax_increment = relax_coeff * (field_2d_ave(k) - field_ref(k))
  DO j = fdims%j_start, fdims%j_end
    DO i = fdims%i_start, fdims%i_end
      field_inc(i,j,k) = field_inc(i,j,k) - relax_increment
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE profile_relax

END MODULE profile_relax_mod
