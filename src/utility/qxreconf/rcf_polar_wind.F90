! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Averaging of Polar wind rows

MODULE Rcf_Polar_Wind_Mod

!  Subroutine Rcf_Polar_Wind - Averaging of Polar wind rows
!
! Description:
!   Performs a vector mean on Polar wind rows.
!
! Method:
!   Magnitude and direction of non-polar rows are calculated and used to
!   set adjecent polar rows.  Normally this is using v to estimate u at pole
!   but ENDGAME requires v at poles so uses u to calculate v.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_POLAR_WIND_MOD'

CONTAINS

SUBROUTINE Rcf_Polar_Wind( polar_wind, non_polar_wind, delta_lon)

USE Rcf_Field_Type_Mod, ONLY: &
field_type

USE UM_ParVars, ONLY: &
atNorth,            &
atSouth,            &
datastart,          &
gc_proc_row_group

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
PrintStatus,            &
PrStatus_Normal

USE conversions_mod, ONLY: &
    pi_over_180

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type), INTENT(INOUT)  :: polar_wind
TYPE( field_type), INTENT(IN)     :: non_polar_wind
REAL,              INTENT(IN)     :: delta_lon

! Local Variables.

INTEGER      :: i
INTEGER      :: k
INTEGER      :: info
INTEGER      :: npw_start_pos, pw_start_pos, pole_dir
INTEGER      :: num_pole, pole

REAL         :: a_p
REAL         :: b_p
REAL         :: longitude
REAL         :: mag_vector
REAL         :: dir_vector

REAL         :: wrk( 2 * non_polar_wind % levels)
REAL         :: rwrk(non_polar_wind % row_len, 2 * non_polar_wind % levels)
REAL         :: delta_lon_rad

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_POLAR_WIND'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Finding vector mean for polar wind rows'
  CALL umPrint(umMessage,src='rcf_polar_wind')
END IF

delta_lon_rad = pi_over_180 * delta_lon
! ----------------------------------------------------------------------
! Section 1.   Calculate magnitude and direction of polar vector wind
!              from non-polar component of wind on row around pole.
! ----------------------------------------------------------------------

! Only need data at poles

IF (atNorth .AND. atSouth) THEN
  num_pole = 2
ELSE
  num_pole = 1
END IF

IF (atNorth .OR. atSouth) THEN
  DO pole = 1, num_pole

    ! Calculate start index in data depending on which pole we are.
    ! Pole_dir gives the sign of the result (cos(t)+sin(t)) after
    ! rotating a vector.
    IF (atNorth .AND. pole == 1 ) THEN
      npw_start_pos = (non_polar_wind % rows - 1) * non_polar_wind % row_len
      pw_start_pos = (polar_wind % rows - 1) * polar_wind % row_len
      pole_dir = 1.0
    ELSE IF (atSouth) THEN
      npw_start_pos = 0
      pw_start_pos = 0
      pole_dir = -1.0
    END IF

    DO k = 1, 2 * non_polar_wind % levels
      wrk(k) = 0.0
    END DO

    DO k = 1, non_polar_wind % levels
      DO i = 1, non_polar_wind % row_len
        ! This strange way of writing things is to fox the Cray
        ! optimiser and thus retain bit comparison (compiler bug)
        longitude = (datastart(1) + i - 2) * delta_lon_rad
        rwrk(i,2*(k-1)+1) = non_polar_wind % DATA( npw_start_pos + i, k) * &
                            COS( longitude )
        rwrk(i,2*k)       = non_polar_wind % DATA( npw_start_pos + i, k) * &
                            SIN( longitude )
      END DO
    END DO

    CALL gcg_rvecsumr(non_polar_wind % row_len, non_polar_wind % row_len, &
                      1, 2 * non_polar_wind % levels, rwrk,               &
                      gc_proc_row_group, info, wrk)

    ! If estimating v using u so rotation is needed in the opposite
    ! direction.
    IF (polar_wind % stashmaster % grid_type == 19) THEN
      pole_dir = -1.0 * pole_dir
    END IF

    DO k = 1, non_polar_wind % levels

      a_p = 2.0 * wrk(2*(k-1)+1) / non_polar_wind % glob_row_len
      b_p = 2.0 * wrk(2*k) / non_polar_wind % glob_row_len

      mag_vector = SQRT(a_p*a_p + b_p*b_p)

      IF (a_p  ==  0.0 .AND. b_p  ==  0.0) THEN
        dir_vector = 0.0
      ELSE
        dir_vector = ATAN2 (b_p, a_p)
      END IF

      DO i = 1, polar_wind % row_len
        polar_wind % DATA( pw_start_pos + i, k) = pole_dir *   &
          ( mag_vector * SIN( (datastart(1) + i - 1 - .5) *    &
          delta_lon_rad - dir_vector ) )
      END DO
    END DO

  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Polar_Wind

END MODULE Rcf_Polar_Wind_Mod
