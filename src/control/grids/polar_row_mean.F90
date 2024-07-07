! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_Row_Mean
! Purpose:
!         Resets polar values to mean over row to remove any
!         deviations due to rounding error.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
MODULE polar_row_mean_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'POLAR_ROW_MEAN_MOD'

CONTAINS

SUBROUTINE polar_row_mean                                                      &
  ( field, ini_start, ini_end, inj_start, inj_end, ink_start, ink_end          &
  , global_row_length, n_proc, n_procy, proc_row_group, at_extremity )

USE global_2d_sums_mod, ONLY: global_2d_sums
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams

IMPLICIT NONE

INTEGER, INTENT(IN)  ::  &
  ini_start,ini_end      & ! Row_length array bounds
, inj_start,inj_end      & ! Rows array bounds
, ink_start,ink_end      & ! Levels array bounds
, n_proc                 & ! Number of procs
, n_procy                & ! Number of procs N/S
, global_row_length      & ! Number of points per global row
, proc_row_group           ! Group id for processors on the same row

LOGICAL ::                                                                     &
  at_extremity(4)          ! Indicates if this processor is at north,
                           ! south, east or west of the processor grid

REAL, INTENT(INOUT) ::                                                         &
  field (ini_start:ini_end, inj_start:inj_end, ink_start:ink_end)


! Local Variables.
INTEGER  ::              &
  i,k                    & ! Loop indices
, info, levels, rlength


REAL     ::                                                                    &
  polar_row (ini_start:ini_end,ink_start:ink_end)                              &
, tmp_sum (ink_start:ink_end)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'POLAR_ROW_MEAN'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1.   Reset polar values to mean over the row.
! ----------------------------------------------------------------------

rlength = ini_end - ini_start + 1
levels  = ink_end - ink_start + 1

! Variables on all model levels.

IF (n_procy  >   1) THEN

  IF (at_extremity(PNorth)) THEN
    DO k=ink_start, ink_end
      DO i=ini_start, ini_end
        polar_row(i,k) = field(i,inj_end,k)
      END DO
    END DO

  ELSE IF (at_extremity(PSouth)) THEN
    DO k=ink_start, ink_end
      DO i=ini_start, ini_end
        polar_row(i,k) = field(i,inj_start,k)
      END DO
    END DO

  ELSE
    DO k=ink_start, ink_end
      DO i=ini_start, ini_end
        polar_row(i,k) = 0.0
      END DO
    END DO

  END IF

  CALL global_2d_sums( polar_row, rlength, 1, 0, 0, levels                     &
                     , tmp_sum, proc_row_group )

  DO i=ink_start, ink_end
    tmp_sum(i) = tmp_sum(i) / global_row_length
  END DO

  IF (at_extremity(PNorth)) THEN
    DO k=ink_start, ink_end
      DO i=ini_start, ini_end
        field(i,inj_end,k) = tmp_sum(k)
      END DO
    END DO


  ELSE IF (at_extremity(PSouth)) THEN

    DO k=ink_start, ink_end
      DO i=ini_start, ini_end
        field(i,inj_start,k) = tmp_sum(k)
      END DO
    END DO

  END IF

ELSE


  ! One processor N/S so must have both North and South
  ! Do north pole

  DO k=ink_start, ink_end
    DO i=ini_start, ini_end
      polar_row(i,k) = field(i,inj_end,k)
    END DO
  END DO

  CALL global_2d_sums( polar_row, rlength, 1, 0, 0, levels                     &
                     , tmp_sum, proc_row_group )

  DO i=ink_start, ink_end
    tmp_sum(i) = tmp_sum(i) / global_row_length
  END DO

  DO k=ink_start, ink_end
    DO i=ini_start, ini_end
      field(i,inj_end,k) = tmp_sum(k)
    END DO
  END DO



  ! Do south pole

  DO k=ink_start, ink_end
    DO i=ini_start, ini_end
      polar_row(i,k) = field(i,inj_start,k)
    END DO
  END DO


  CALL global_2d_sums( polar_row, rlength, 1, 0, 0, levels                     &
                     , tmp_sum, proc_row_group )

  DO i=ink_start, ink_end
    tmp_sum(i) = tmp_sum(i) / global_row_length
  END DO

  DO k=ink_start, ink_end
    DO i=ini_start, ini_end
      field(i,inj_start,k) = tmp_sum(k)
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Polar_Row_Mean

END MODULE polar_row_mean_mod
