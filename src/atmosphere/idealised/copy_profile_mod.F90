! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE copy_profile_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COPY_PROFILE_MOD'

CONTAINS

SUBROUTINE copy_profile(num_times, num_heights, tsec, height,                  &
                        profile_data, profile, field_type, timescale)

! Purpose:
!          Copies profiles read from namelist into derived type
!          that holds time-varying profile data.
!

USE profiles_mod,     ONLY: varying_profile, num_data_max, eta_coord, z_coord, &
                            p_coord, l_p_profile

USE missing_data_mod, ONLY: imdi

USE Ereport_mod,      ONLY: Ereport
USE errormessagelength_mod,                                                    &
                      ONLY: errormessagelength

USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim

IMPLICIT NONE

! Input
INTEGER, INTENT(IN)                :: num_times, num_heights
REAL, INTENT(IN)                   :: tsec(num_data_max)
REAL, INTENT(IN)                   :: height(num_data_max)
REAL, INTENT(IN)                   :: profile_data(num_data_max)

! Optional arguments
INTEGER, OPTIONAL, INTENT(IN)      :: field_type
REAL,    OPTIONAL, INTENT(IN)      :: timescale

! Output
TYPE(varying_profile), INTENT(OUT) :: profile

! Local
INTEGER                            :: i, k, n_minus_k
INTEGER                            :: ind
REAL                               :: tmp_swap

! Error reporting
CHARACTER(LEN=*), PARAMETER        :: RoutineName = 'COPY_PROFILE'
CHARACTER(LEN=errormessagelength)  :: Cmessage
INTEGER                            :: ErrorStatus

! Dr Hook
INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

profile % n_times   = 0
profile % n_heights = 0

IF (PRESENT(field_type)) profile % field_type = field_type
IF (PRESENT(timescale))  profile % timescale  = timescale

IF (num_times /= imdi .AND.  num_times > 0) THEN
  profile % n_times   = num_times
  IF (num_heights <= 0) THEN
    WRITE(Cmessage,'(''Improper number of heights = '', I0)') num_heights
    ErrorStatus = 1
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  ELSE
    profile % n_heights = num_heights
    ALLOCATE(profile % tsec(num_times))
    ALLOCATE(profile % height(num_heights))
    ALLOCATE(profile % vprof(num_heights, num_times))
    profile % tsec(:)   = tsec(1:num_times)
    profile % height(:) = height(1:num_heights)
    ind = 0
    DO i = 1, num_times
      DO k = 1, num_heights
        ind = ind + 1
        profile % vprof(k,i) = profile_data(ind)
      END DO
    END DO

! Determine the coordinate type
    IF (MINVAL(height(1:num_heights)) >= 0.0 .AND.                             &
        MAXVAL(height(1:num_heights)) <= 1.0 ) THEN
      profile % coord_type = eta_coord
    ELSE IF (height(num_heights) < height(1)) THEN
      profile % coord_type = p_coord
      l_p_profile      = .TRUE.
    ELSE
      profile % coord_type = z_coord
    END IF

! For pressure coordinate reverse the order of the height index.
! This allows prof_interp to be used.
    IF (profile % coord_type == p_coord) THEN
      DO k = 1, num_heights/2
        n_minus_k                   = num_heights - k + 1
        tmp_swap                    = profile % height(k)
        profile % height(k)         = profile % height(n_minus_k)
        profile % height(n_minus_k) = tmp_swap
      END DO
      DO i = 1, num_times
        DO k = 1, num_heights/2
          n_minus_k                    = num_heights - k + 1
          tmp_swap                     = profile % vprof(k,i)
          profile % vprof(k,i)         = profile % vprof(n_minus_k,i)
          profile % vprof(n_minus_k,i) = tmp_swap
        END DO
      END DO
    END IF

  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE copy_profile

END MODULE copy_profile_mod
