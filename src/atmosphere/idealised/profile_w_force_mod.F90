! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine profile_w_force
!
MODULE profile_w_force_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROFILE_W_FORCE_MOD'

CONTAINS

SUBROUTINE profile_w_force(profile, fdims, ftype, w_force)

! Purpose: To calculate the w profile prior to calculating the forcing
!
! Method:  Work out the w profile for this time step.
!

USE timestep_mod,                 ONLY: timestep
USE atm_fields_bounds_mod,        ONLY: array_dims
USE field_types,                  ONLY: fld_type_u, fld_type_v, fld_type_w
USE level_heights_mod,            ONLY: eta_theta_levels, eta_rho_levels,      &
                                        z_ref_theta, z_ref_rho
USE idealise_run_mod,             ONLY: zero_orography
USE profiles_mod,                 ONLY: varying_profile,                       &
                                        eta_coord, p_coord, z_coord,           &
                                        p_prof_theta, p_prof_rho
USE prof_temporal_interp_mod,     ONLY: prof_temporal_interp
USE prof_interp_mod,              ONLY: prof_interp, interp_linear,            &
                                        extrap_constant


USE umPrintMgr,                   ONLY: PrintStatus, PrStatus_Diag,            &
                                        umMessage, umPrint
USE ereport_mod,                  ONLY: ereport
USE errormessagelength_mod,       ONLY: errormessagelength
USE yomhook,                      ONLY: lhook, dr_hook
USE parkind1,                     ONLY: jprb, jpim
  
IMPLICIT NONE

! Input
! ---------------------------------------------------------------------
TYPE(varying_profile), INTENT(IN) :: profile
TYPE(array_dims),      INTENT(IN) :: fdims     ! array bounds for data

INTEGER,               INTENT(IN) :: ftype

! ---------------------------------------------------------------------
! w profile for forcing
REAL, INTENT(INOUT)           :: w_force(fdims%k_start:fdims%k_end)

! Local
! ---------------------------------------------------------------------
INTEGER                       :: k
REAL                          :: prof_now(profile % n_heights)
REAL, POINTER                 :: levels(:)
REAL                          :: dz     ! layer depth (m)
REAL                          :: dzdt   ! layer depth / time step

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus
CHARACTER(LEN=errormessagelength)                                              &
                              :: Cmessage
CHARACTER(LEN=*), PARAMETER   :: RoutineName='PROFILE_W_FORCE'

! ---------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that orography is flat
IF (.NOT. zero_orography) THEN
  ErrorStatus = 1
  Cmessage = 'w forcing not coded for non-zero orography'
  CALL ereport(RoutineName, ErrorStatus, Cmessage)
END IF

! Interpolate reference profile to current time 
CALL prof_temporal_interp(profile, prof_now)

! Optional print out for checking 
IF (PrintStatus == PrStatus_Diag) THEN
  WRITE(umMessage,'(A)')  'Height  w profile now'
  CALL umPrint(umMessage,src=RoutineName)
  DO k=1,profile % n_heights
    WRITE(umMessage,'(2F18.6)') profile % height(k),prof_now(k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

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


! Work out w profile for this time
CALL prof_interp(interp_linear, extrap_constant, profile % n_heights,          &
                 profile % height, prof_now, fdims%k_end-fdims%k_start+1,      &
                 levels, w_force)

! Check w profile will not break CFL limit  i.e. w < dz/dt
! If w value breaks CFL print a warning as w forcing code is not designed
! to cope with w > dz/dt
! Note this check is using layer depths assuming no orography.

DO k = 1, fdims%k_end-1
  dz = z_ref_rho(k+1) - z_ref_rho(k)
  dzdt = dz/timestep
  IF ( w_force(k) >= dzdt) THEN
    ErrorStatus = -1
    WRITE(Cmessage, '(A,I3,F12.6,A,F12.6)')  "WARNING w forcing level ",  &
                  k,w_force(k)," exceeds CFL limit ",dzdt
    CALL Ereport(RoutineName, ErrorStatus, Cmessage)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE profile_w_force

END MODULE profile_w_force_mod
