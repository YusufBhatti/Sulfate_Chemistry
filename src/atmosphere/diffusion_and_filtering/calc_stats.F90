! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT*****************************
!
! Subroutine calc_stats

SUBROUTINE calc_stats(                                            &
                      field,                                      &
                      row_length, rows, levels,                   &
                      halo_x, halo_y, n_proc, array_size,         &
                      mean, variance, std_dev, skew, kurtosis)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
IMPLICIT NONE
!
! description:
! Standard statistics of a field
!
! method:
! Do global sums  on processor 0 only.
! Results from different processor configuration will not
! bit-reproduce but model fields are not affected (diagnostic only)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! subroutine arguments:

INTEGER ::                                                        &
  row_length                                                      &
                     ! in columns
, rows                                                            &
                     ! in rows
, levels                                                          &
                     ! no. of levels
, array_size                                                      &
                     ! array size for summing over pe's
, halo_x                                                          &
                     ! in halo in i direction ew
, halo_y                                                          &
                     ! in halo in j direction ns
, n_proc             ! Total number of processors

! inout
REAL ::                                                           &
 field(1-halo_x: row_length+halo_x, 1-halo_y: rows+halo_y         &
       , levels)            ! in array containing fields

! inout
REAL ::                                                           &
  mean(levels)                                                    &
, std_dev(levels)                                                 &
, variance(levels)                                                &
, skew(levels)                                                    &
, kurtosis(levels)
!
! local variables
REAL ::                                                            &
  summ(array_size)                                                 &
             ! various sums
, points                                                          &
                   ! field length (local and global)
, recip_points                                                    &
                   ! reciprocal field length
, diff                                                            &
                 ! difference
, diff_power     ! various powers of difference

INTEGER ::                                                        &
  i, j, k, k1                                                     &
                 ! loop variables
, info           ! return code from gc stuff

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_STATS'


!-----------------------------------------------------------------------
! mpp code for summation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 1. Do the sums
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO k = 1, levels
  summ(k) = 0.0   ! to sum field
  DO j = 1, rows
    DO i = 1, row_length
      summ(k) = summ(k) + field(i,j,k)
    END DO
  END DO
END DO  ! k = 1, levels

!  Put local number of points in sum(levels + 1)
k1 = levels + 1
summ(k1) = row_length * rows  ! to field size

CALL gc_rsumr(k1, n_proc, info, summ)

!  Total number of points is in sum(levels + 1)
points = summ(levels+1)
recip_points = 1.0 / points

DO k = 1, levels
  mean(k) = summ(k) * recip_points
END DO

!   summ(1) for difference from mean, summ(2) for variance
!   summ(3) for skew, summ(4) for kurtosis
DO i = 1, array_size
  summ(i) = 0.0   ! 4 sums required
END DO

DO k = 1, levels
  k1 = 4 * k - 3
  summ(k1) = 0.0
  summ(k1+1) = 0.0
  summ(k1+2) = 0.0
  summ(k1+3) = 0.0
  DO j = 1, rows
    DO i = 1, row_length
      diff  = field(i,j,k) - mean(k)
      summ(k1) = summ(k1) + diff
      diff_power  = diff * diff
      summ(k1+1) = summ(k1+1) + diff_power  !  variance
      diff_power  = diff_power * diff
      summ(k1+2) = summ(k1+2) + diff_power  !  skew
      diff_power  = diff_power * diff
      summ(k1+3) = summ(k1+3) + diff_power  !  kurtosis
    END DO
  END DO
END DO

k1 = 4 * levels
CALL gc_rsumr(k1, n_proc, info, summ)

DO k = 1, levels
  k1 = 4 * k - 3
  variance(k) = (summ(k1+1) - summ(k1) * summ(k1) * recip_points) /  &
                                                    (points - 1.0)
  IF ( variance(k) > 0 ) THEN
    std_dev(k) = SQRT(variance(k))
    skew(k) = summ(k1+2) * recip_points / std_dev(k)**3
    kurtosis(k) = summ(k1+3) * recip_points /                      &
                                   (variance(k) * variance(k)) - 3
  ELSE
    std_dev(k) = 0.0
    skew(k) = 0.0
    kurtosis(k) = 0.0
  END IF  !  variance(k) > 0
END DO !  k = 1, levels

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_stats
