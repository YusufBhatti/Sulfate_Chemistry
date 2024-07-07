! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description: Compute UKCA diagnostics on pressure levels
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Method:      Uses methodology in EOT diags for Q on PRESSURE LEVELS
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
!  Code description:
!    Language: FORTRAN 90
!
! ---------------------------------------------------------------------

MODULE ukca_calc_plev_diag_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CALC_PLEV_DIAG_MOD'

CONTAINS

SUBROUTINE ukca_calc_plev_diag( num_stash_levels, stash_levels,   &
    row_length, rows, model_levels, out_plevs,                    &
    diag_in,                                                      &
    exner_at_theta,                                               &
    diag_out, hvyside_out)

USE planet_constants_mod,  ONLY: p_zero, kappa
USE UM_parvars, ONLY: at_extremity

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: num_stash_levels  
INTEGER, INTENT(IN) :: row_length   ! dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: out_plevs
INTEGER, INTENT(IN) :: stash_levels(num_stash_levels+1 )

! Input arrays
REAL, INTENT(IN):: diag_in(row_length,rows,model_levels) 
REAL, INTENT(IN):: exner_at_theta(row_length,rows,model_levels)

! Output variables
REAL, INTENT(OUT):: diag_out(row_length,rows,out_plevs)
REAL, INTENT(OUT):: hvyside_out(row_length,rows,out_plevs)

! Local variables
! Interpolation type: 1 = linear, 2 = cubic, 3 = quintic
INTEGER, PARAMETER  :: interp_order=1  
INTEGER :: k

REAL    :: p_levels(num_stash_levels) ! Pressure levels requested in STASH
REAL    :: pressure_pa, pressure_ex

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_PLEV_DIAG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Loop over output pressure levels 
DO k = 1, out_plevs

  p_levels(k) = stash_levels(k+1) / 1000.0   ! Get levs from STASH

  pressure_pa = p_levels(k)*100.             ! convert to Pascals
  pressure_ex = ( pressure_pa / p_zero )**kappa   ! exner

  CALL ukca_vert_interp2(diag_in,                               &
       row_length, rows, model_levels,                          &
       pressure_ex,                                             &
! halo sizes
       0, 0, 0, 0,                                              &
       exner_at_theta, interp_order,                            &
       diag_out(:,:,k), hvyside_out(:,:,k) )

END DO ! over output STASH pressure levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


CONTAINS
SUBROUTINE ukca_vert_interp2(                             &
                   data_in, row_length, data_rows,        &
                   data_levels, desired_p,                &
                   halo_x1, halo_y1,                      &
                   halo_x2, halo_y2,                      &
                   p_at_data, interp_order,               &
                   data_out, hvside )

! ---------------------------------------------------------------------
! Purpose:
!          Performs vertical interpolation of a field to a
!          desired p surface given the value of p at each of
!          the data points. Where the desired surface is below/above
!          the bottom/top data point Missing Data value is returned.
!
!          p can be any quantity that is monotonic decreasing with
!          model level, such as pressure, or Exner pressure.
!          It is usual, and recommended, that the p surface be exner.
!          This routine can  do linear, cubic or quintic interpolation,
!          as specified by the interp_order argument.
!
! Method:
!          Lagrange interpolation.
!
! ---------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT (IN) :: row_length     ! number of points on a row
INTEGER, INTENT (IN) :: data_rows      ! number of rows of data
INTEGER, INTENT (IN) :: data_levels    ! number of levels of data
INTEGER, INTENT (IN) :: interp_order   ! 1 = linear, 3 = cubic, 5=quintic
INTEGER, INTENT (IN) :: halo_x1
INTEGER, INTENT (IN) :: halo_x2
INTEGER, INTENT (IN) :: halo_y1
INTEGER, INTENT (IN) :: halo_y2

REAL, INTENT (IN) :: desired_p         ! desired value to which data
                                       ! should be  interpolated to.

REAL, INTENT (IN) :: data_in (1-halo_x1:row_length+halo_x1,       &
                              1-halo_y1:data_rows+halo_y1, data_levels)
REAL, INTENT (IN) :: p_at_data (1-halo_x2:row_length+halo_x2,     &
                                1-halo_y2:data_rows+halo_y2, data_levels)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT (OUT) :: data_out (row_length, data_rows)
REAL, INTENT (OUT) :: hvside (row_length, data_rows)

! Local variables
INTEGER :: i, j, k   ! Loopers
INTEGER :: level_below(row_length, data_rows)
INTEGER :: last  ! level from last loop

REAL :: p_here
REAL :: p_here_plus
REAL :: p_here_plus2
REAL :: p_here_minus
REAL :: p_here_plus3
REAL :: p_here_minus2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_VERT_INTERP2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!   (if requested level is above top of model, set to -1, which will
!    use value from top level)
! ----------------------------------------------------------------------
#if defined(VECTOR)
DO j = 1, data_rows
  DO i = 1, row_length
    level_below(i,j) = -1
  END DO
END DO


DO k = 1, data_levels
  DO j = 1, data_rows
    DO i = 1, row_length
      IF ((p_at_data(i,j,k) < desired_p) .AND.                    &
          (level_below(i,j) == -1)) THEN
        level_below(i,j) = k-1
      END IF
    END DO
  END DO
END DO
#else

! For scalar platforms we can take advantage of the fact that fields
! tend to vary smoothly by storing the last level we found and searching
! up or down from this point until the right level is found. This is not as
! simple to follow, but much quicker.

last=1
DO j = 1, data_rows
  DO i = 1, row_length

! Searching up. Default is -1 (i.e. top level)
    level_below(i,j) = -1
    IF (p_at_data(i,j,last) >= desired_p) THEN
      DO k = last+1, data_levels
        IF (p_at_data(i,j,k) < desired_p) THEN
          level_below(i,j) = k-1
          EXIT
        END IF
      END DO

    ELSE

! Searching down. Default is 0 (i.e. bottom level)
      level_below(i,j) = 0
      DO k = last-1, 1, -1
        IF (p_at_data(i,j,k) >= desired_p) THEN
          level_below(i,j) = k
          EXIT
        END IF
      END DO

    END IF

    last = MAX(level_below(i,j), 1)
  END DO
END DO
#endif

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------
hvside(:,:) = 1.0           ! Heavyside fn to track levels below Psurf
DO j = 1, data_rows
  DO i = 1, row_length

    IF (level_below(i,j)  ==  -1) THEN
      data_out(i,j) = data_in(i, j, data_levels)

    ELSE IF (level_below(i,j)  ==  0) THEN
      data_out(i,j) = 0.0
      hvside(i,j) = 0.0

    ELSE IF (level_below(i,j)  ==  1 .OR.                         &
             level_below(i,j)  ==  data_levels - 1                &
             .OR. interp_order  ==  1 ) THEN
! linearly interpolate
      data_out (i,j) = ( (desired_p -                             &
                          p_at_data(i,j,level_below(i,j)) )       &
                          * data_in (i,j,level_below(i,j)+1)      &
                         -(desired_p -                            &
                           p_at_data(i,j,level_below(i,j)+1)) *   &
                           data_in (i,j,level_below(i,j)) ) /     &
                        ( p_at_data(i,j,level_below(i,j)+1) -     &
                          p_at_data(i,j,level_below(i,j)) )

    ELSE IF (level_below(i,j)  ==  2 .OR.                         &
             level_below(i,j)  ==  data_levels - 2                &
             .OR. interp_order  ==  3 ) THEN
! cubicly interpolate
      p_here_minus = p_at_data(i,j,level_below(i,j)-1)
      p_here = p_at_data(i,j,level_below(i,j))
      p_here_plus = p_at_data(i,j,level_below(i,j)+1)
      p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)

      data_out (i,j) = ( (desired_p - p_here) *                   &
                                (desired_p - p_here_plus )*       &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here_minus - p_here) *         &
                                (p_here_minus - p_here_plus )*    &
                                (p_here_minus - p_here_plus2 ) ) *&
                              data_in (i,j,level_below(i,j)-1) +  &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here_plus )*       &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here - p_here_minus) *         &
                                (p_here - p_here_plus )*          &
                                (p_here - p_here_plus2 ) ) *      &
                              data_in (i,j,level_below(i,j)) +    &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here )*            &
                                (desired_p - p_here_plus2 ) ) /   &
                              ( (p_here_plus - p_here_minus) *    &
                                (p_here_plus - p_here )*          &
                                (p_here_plus - p_here_plus2 ) ) * &
                              data_in (i,j,level_below(i,j)+1) +  &
                              ( (desired_p - p_here_minus) *      &
                                (desired_p - p_here )*            &
                                (desired_p - p_here_plus ) ) /    &
                              ( (p_here_plus2 - p_here_minus) *   &
                                (p_here_plus2 - p_here )*         &
                                (p_here_plus2 - p_here_plus ) ) * &
                              data_in (i,j,level_below(i,j)+2)


    ELSE
! interpolate quinticly
      p_here_minus2 = p_at_data(i,j,level_below(i,j)-2)
      p_here_minus = p_at_data(i,j,level_below(i,j)-1)
      p_here = p_at_data(i,j,level_below(i,j))
      p_here_plus = p_at_data(i,j,level_below(i,j)+1)
      p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)
      p_here_plus3 = p_at_data(i,j,level_below(i,j)+3)

      data_out (i,j) = ((desired_p - p_here_minus) *              &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_minus2 - p_here_minus) *    &
                              (p_here_minus2 - p_here )*          &
                              (p_here_minus2 - p_here_plus )*     &
                              (p_here_minus2 - p_here_plus2 )*    &
                              (p_here_minus2 - p_here_plus3 ) ) * &
                              data_in (i,j,level_below(i,j)-2) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_minus - p_here_minus2) *    &
                              (p_here_minus - p_here )*           &
                              (p_here_minus - p_here_plus )*      &
                              (p_here_minus - p_here_plus2 )*     &
                              (p_here_minus - p_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)-1) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here - p_here_minus2) *          &
                              (p_here - p_here_minus )*           &
                              (p_here - p_here_plus )*            &
                              (p_here - p_here_plus2 )*           &
                              (p_here - p_here_plus3 ) ) *        &
                              data_in (i,j,level_below(i,j)) +    &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus2 )*        &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_plus - p_here_minus2) *     &
                              (p_here_plus - p_here_minus )*      &
                              (p_here_plus - p_here )*            &
                              (p_here_plus - p_here_plus2 )*      &
                              (p_here_plus - p_here_plus3 ) ) *   &
                              data_in (i,j,level_below(i,j)+1) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus3 ))/       &
                            ( (p_here_plus2 - p_here_minus2) *    &
                              (p_here_plus2 - p_here_minus )*     &
                              (p_here_plus2 - p_here )*           &
                              (p_here_plus2 - p_here_plus )*      &
                              (p_here_plus2 - p_here_plus3 ) ) *  &
                              data_in (i,j,level_below(i,j)+2) +  &
                            ((desired_p - p_here_minus2) *        &
                              (desired_p - p_here_minus )*        &
                              (desired_p - p_here )*              &
                              (desired_p - p_here_plus )*         &
                              (desired_p - p_here_plus2 ))/       &
                            ( (p_here_plus3 - p_here_minus2) *    &
                              (p_here_plus3 - p_here_minus )*     &
                              (p_here_plus3 - p_here )*           &
                              (p_here_plus3 - p_here_plus )*      &
                              (p_here_plus3 - p_here_plus2 ) ) *  &
                              data_in (i,j,level_below(i,j)+3)

    END IF
  END DO
END DO

! end of routine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_vert_interp2

END SUBROUTINE UKCA_CALC_PLEV_DIAG

END MODULE ukca_calc_plev_diag_mod

