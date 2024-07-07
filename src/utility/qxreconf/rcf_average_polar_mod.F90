! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  performs polar row averaging for theta fields

MODULE Rcf_average_polar_mod

!  Subroutine Rcf_average_polar_real    - for real data
!  Subroutine Rcf_average_polar_int     - for integer data
!  Subroutine Rcf_average_polar_logical - for logical data
!
! Description:
!   Polar row averaging perfomed on serial (non distributed) fields.
!   Horizontal interpolation done on a level-by-level basis.
!
! Method:
!   A generic interface is used to determine which averaging method
!   is utilised
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

!-------------------------------------------------------------------
! Use a generic (overloaded) interface for different data-types
!-------------------------------------------------------------------
INTERFACE Rcf_average_polar
MODULE PROCEDURE Rcf_average_polar_real, Rcf_average_polar_int, &
                 Rcf_average_polar_logical
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_AVERAGE_POLAR_MOD'

CONTAINS

!------------------------------------------------------------------
! For real data
!------------------------------------------------------------------
SUBROUTINE Rcf_average_polar_real(field, rows, row_length, global, &
                                  average_required)

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)   :: rows
INTEGER, INTENT(IN)   :: row_length
REAL, INTENT(INOUT)   :: field(row_length * rows)
LOGICAL, INTENT(IN)   :: global
LOGICAL, INTENT(OUT)  :: average_required

! Local variables
INTEGER           :: i
INTEGER           :: last_row
REAL              :: sum_start
REAL              :: sum_end
REAL              :: avg_start
REAL              :: avg_end
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_AVERAGE_POLAR_REAL'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
last_row = (rows - 1) * row_length
sum_start = 0.0
sum_end   = 0.0
average_required = .FALSE.

IF (global) THEN
  DO i = 1, row_length

    IF ( field( i ) /= field( 1 ) .OR.                              &
         field( last_row + i ) /= field( last_row + 1 ) ) THEN
      average_required = .TRUE.
    END IF

    sum_start = sum_start + field( i )
    sum_end = sum_end + field( last_row + i )
  END DO

  IF ( average_required ) THEN
    avg_start = sum_start/REAL(row_length)
    avg_end   = sum_end/REAL(row_length)

    DO i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    END DO
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Average_polar_real

!------------------------------------------------------------------
! For integer data
!------------------------------------------------------------------
SUBROUTINE Rcf_average_polar_int(field, rows, row_length, global, &
                                 average_required)

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)        :: rows
INTEGER, INTENT(IN)        :: row_length
INTEGER, INTENT(INOUT)     :: field(row_length * rows)
LOGICAL, INTENT(IN)        :: global
LOGICAL, INTENT(OUT)       :: average_required

! Local variables
INTEGER           :: i
INTEGER           :: last_row
INTEGER           :: sum_start
INTEGER           :: sum_end
INTEGER           :: avg_start
INTEGER           :: avg_end
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_AVERAGE_POLAR_INT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
last_row  = (rows - 1) * row_length
sum_start = 0
sum_end   = 0
average_required = .FALSE.

IF (global) THEN
  DO i = 1, row_length

    IF ( field( i ) /= field( 1 ) .OR.                              &
         field( last_row + i ) /= field( last_row + 1 ) ) THEN
      average_required = .TRUE.
    END IF

    sum_start = sum_start + field( i )
    sum_end = sum_end + field( last_row + i )
  END DO

  IF ( average_required ) THEN
    avg_start = NINT( REAL(sum_start)/REAL(row_length) )
    avg_end   = NINT( REAL(sum_end)  /REAL(row_length) )

    DO i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    END DO
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Average_polar_int

!------------------------------------------------------------------
! For logical data
!------------------------------------------------------------------
SUBROUTINE Rcf_average_polar_logical(field, rows, row_length, global, &
                                     average_required)

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)     :: rows
INTEGER, INTENT(IN)     :: row_length
LOGICAL, INTENT(INOUT)  :: field(row_length * rows)
LOGICAL, INTENT(IN)     :: global
LOGICAL, INTENT(OUT)    :: average_required

! Local variables
INTEGER           :: i
INTEGER           :: last_row
INTEGER           :: sum_start
INTEGER           :: sum_end
LOGICAL           :: avg_start
LOGICAL           :: avg_end
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_AVERAGE_POLAR_LOGICAL'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
last_row = (rows - 1) * row_length
sum_start = 0
sum_end   = 0
average_required = .FALSE.

IF (global) THEN
  DO i = 1, row_length

    IF ( field( i ) .neqv. field( 1 ) .OR.                         &
         field( last_row + i ) .neqv. field( last_row + 1 ) ) THEN
      average_required = .TRUE.
    END IF

    IF ( field(i) ) THEN
      sum_start = sum_start + 1
    END IF

    IF ( field( last_row + i) ) THEN
      sum_end = sum_end + 1
    END IF
  END DO

  IF (average_required) THEN
    IF ( sum_start >= row_length/2 ) THEN
      avg_start = .TRUE.
    ELSE
      avg_start = .FALSE.
    END IF

    IF ( sum_end >= row_length/2 ) THEN
      avg_end = .TRUE.
    ELSE
      avg_end = .FALSE.
    END IF

    DO i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    END DO

  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Average_polar_logical

END MODULE Rcf_Average_Polar_Mod
