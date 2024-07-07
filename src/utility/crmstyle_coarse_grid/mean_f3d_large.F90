! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculate means and variances on coarse grid for fields like theta

MODULE mean_f3d_large_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MEAN_F3D_LARGE_MOD'

CONTAINS

SUBROUTINE mean_f3d_large(ncols,nrows,nlevs,num_y,num_x,index_row,index_col, &
                         mask, field_in,                                     &
                         field_prime, mean_out, var_out)

USE missing_data_mod, ONLY: rmdi
USE word_sizes_mod, ONLY:  wp


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description:
!  Calculate means and variances on coarse grid taking account of points
!  with no values as below surface. Special version for fields like theta
!  and pressure which have large values and small variations. 32 bit summation
!  is not accurate enough when there are fields of more than 200x200 values
!  to sum.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

INTEGER, INTENT(IN) ::     &
  ncols                    & ! number of columns
 ,nrows                    & ! number of rows
 ,nlevs                    & ! number of model levels
 ,num_y                    & ! Number of rows coarse grid
 ,num_x                      ! Number of columns coarse grid

INTEGER, INTENT(IN) ::     &
  index_row(nrows)         & ! index for row on coarse grid
 ,index_col(ncols)           ! index for column on coarse grid

LOGICAL, INTENT(IN) ::    &
  mask(ncols,nrows,nlevs)   ! mask true if point above ground

REAL(wp), INTENT(IN)  ::         &
  field_in(ncols,nrows,nlevs)      ! input field

REAL(wp), INTENT(OUT) ::         &
  field_prime(ncols,nrows,nlevs) & ! field minus mean value
 ,mean_out(num_x,num_y,nlevs)    & ! output field mean fields
 ,var_out(num_x,num_y,nlevs)       ! output field mean field

! local variables

INTEGER :: i,j,k, ii,jj          ! loop counters

INTEGER ::                       &
  nvalues(num_x,num_y,nlevs)     & ! number of values
 ,imean(num_x,num_y,nlevs)         ! integer values close to mean

REAL ::                          &
  field_temp(ncols,nrows,nlevs)  & ! field with a estimated mean removed
 ,mean_est(num_x,num_y,nlevs)      ! estimated mean for each sub-region.

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MEAN_F3D_LARGE'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Mean_x  = sum x(i)/N       sum over 1-N
!
!  var x = sqrt ( sum(x(i) -x_mean)**2)/N  )  = sum(x(i)*x(i))/N  - mean_x**2
!
!-------------------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(i,j,ii,jj,k) DEFAULT(NONE)                         &
!$OMP& SHARED(nlevs, nrows, ncols, num_x, num_y, index_row, index_col, imean,&
!$OMP&     nvalues, mean_out, var_out, field_in, mask, mean_est, field_temp)
DO k=1,nlevs

  ! Estimate a mean by taking nearest integer of first value in the field for
  ! each sub-region. This will avoid the problem of parallisation giving
  ! different first value if this done for the full grid on a node.
  ! Another problem that the first point could be below the surface. Need to
  !  use a valid point.

  DO jj = 1,num_y
    DO ii = 1,num_x
      nvalues(ii,jj,k)  = 0
      mean_out(ii,jj,k) = 0.0
      var_out(ii,jj,k)  = 0.0
      imean(ii,jj,k)  = -99
      mean_est(ii,jj,k)  = 0.0
    END DO
  END DO

  DO j=1,nrows
    jj = index_row(j)
    DO i=1,ncols
      ii = index_col(i)
      IF (imean(ii,jj,k) == -99 .AND. mask(i,j,k) ) THEN
        imean(ii,jj,k) = INT(field_in(i,j,k))
        mean_est(ii,jj,k) = REAL(imean(ii,jj,k))
      END IF
    END DO
  END DO

  DO j=1,nrows
    jj = index_row(j)
    DO i=1,ncols
      ii = index_col(i)
      field_temp(i,j,k) =   field_in(i,j,k) - mean_est(ii,jj,k)
      IF (mask(i,j,k)) THEN
        mean_out(ii,jj,k)= mean_out(ii,jj,k)+field_temp(i,j,k)
        nvalues(ii,jj,k) = nvalues(ii,jj,k)+ 1
      END IF
    END DO
  END DO

  DO jj = 1,num_y
    DO ii = 1,num_x
      IF (nvalues(ii,jj,k) > 0) THEN
        mean_out(ii,jj,k) = mean_out(ii,jj,k)/REAL(nvalues(ii,jj,k)) +  &
                            mean_est(ii,jj,k)
      ELSE
        mean_out(ii,jj,k) = rmdi
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


! problem with x**2 method as very big numbers and rounding issues

!$OMP PARALLEL DO PRIVATE(i,j,ii,jj,k) DEFAULT(NONE)                      &
!$OMP& SHARED(nlevs, nrows, ncols, num_x, num_y, index_row, index_col,    &
!$OMP&     nvalues, mean_out, var_out, field_in, field_prime, mask)
DO k=1,nlevs
  DO j=1,nrows
    jj = index_row(j)
    DO i=1,ncols
      ii = index_col(i)
      IF (mask(i,j,k)) THEN
        field_prime(i,j,k) = field_in(i,j,k) - mean_out(ii,jj,k)
        var_out(ii,jj,k)= var_out(ii,jj,k)                               &
                               +(field_prime(i,j,k)* field_prime(i,j,k))
      ELSE
        field_prime(i,j,k) = 0.0     ! or should I set it to rmdi
      END IF
    END DO
  END DO

  DO jj = 1,num_y
    DO ii = 1,num_x
      IF (nvalues(ii,jj,k) > 0) THEN
        IF (var_out(ii,jj,k) > 0.0) THEN
          var_out(ii,jj,k) = SQRT(var_out(ii,jj,k)/REAL(nvalues(ii,jj,k)))
        ELSE
          var_out(ii,jj,k) = 0.0
        END IF
      ELSE
        var_out(ii,jj,k) = rmdi
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE mean_f3d_large
END MODULE mean_f3d_large_mod
