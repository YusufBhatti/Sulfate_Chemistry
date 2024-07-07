! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module profile_interp_mod
MODULE prof_temporal_interp_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                  ModuleName='PROF_TEMPORAL_INTERP_MOD'

CONTAINS

SUBROUTINE prof_temporal_interp(profile, prof_now)

! Purpose:
!         Linear interpolation in time of profile to cuurent time. 
!
!
USE profiles_mod,        ONLY: varying_profile
USE timestep_mod,        ONLY: timestep, timestep_number

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

! Input
! - time-varying profile to interpolate to current time
TYPE(varying_profile), INTENT(IN) :: profile

! Output
! - input profile interpolated to current time
REAL, INTENT(OUT)                 :: prof_now(profile % n_heights)

! Local variables -
! - loop counters
INTEGER                           :: i
! - current model time in seconds
REAL                              :: t
! - temporal interpolation weight
REAL                              :: tau

INTEGER(KIND=jpim), PARAMETER     :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER     :: zhook_out = 1
REAL(KIND=jprb)                   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PROF_TEMPORAL_INTERP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Temporal interpolation not required, if there's only one input profile
IF (profile % n_times == 1) THEN

  prof_now(:) = profile % vprof(:,1)

ELSE

! Current model time
  t = timestep * timestep_number

! Temporal interplation of input profiles to current time
  IF ( t < profile % tsec(1) ) THEN
    prof_now(:) = profile % vprof(:,1)
  ELSE IF ( t >= profile % tsec(profile % n_times) ) THEN
    prof_now(:) = profile % vprof(:,profile % n_times)
  ELSE
    DO i = 1, profile % n_times - 1
      IF (profile % tsec(i+1) > t) THEN
        tau = (t - profile % tsec(i)) /                                      &
              (profile % tsec(i+1) - profile % tsec(i))
        prof_now(:) = (1.0 - tau) * profile % vprof(:,i) +                   &
                              tau * profile % vprof(:,i+1)
        EXIT
      END IF
    END DO
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE prof_temporal_interp

END MODULE prof_temporal_interp_mod
