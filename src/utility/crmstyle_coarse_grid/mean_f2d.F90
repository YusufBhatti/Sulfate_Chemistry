! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculate means and variances on coarse grid for a 2d field

MODULE mean_f2d_mod

USE word_sizes_mod, ONLY:  wp

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MEAN_F2D_MOD'

CONTAINS



SUBROUTINE mean_f2d(ncols,nrows,num_y,num_x,new_res,index_row,index_col,  &
                         field_in, mean_out, var_out)

USE missing_data_mod, ONLY: rmdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! Description:
!  Calculate means and variances on coarse grid
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   This code is written to UM programming standards version 8.3.

INTEGER, INTENT(IN) ::     &
  ncols                    & ! number of columns
 ,nrows                    & ! number of rows
 ,num_y                    & ! Number of rows coarse grid
 ,num_x                    & ! Number of columns coarse grid
 ,new_res                    ! Number of rows/columns making new coarse gridbox

INTEGER, INTENT(IN) ::     &
  index_row(nrows)         & ! index for row on coarse grid
 ,index_col(ncols)           ! index for column on coarse grid

REAL(wp), INTENT(IN)  ::   &
  field_in(ncols,nrows)      ! input field

REAL(wp), INTENT(OUT) ::   &
  mean_out(num_x,num_y)    & ! output field mean fields
 ,var_out(num_x,num_y)       ! output field mean field

! local variables

INTEGER :: i,j,k, ii,jj, iii, jjj        ! loop counters

REAL :: nvalues,rnvalues    ! number of values and 1/nvalues

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MEAN_F2D'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Mean_x  = sum x(i)/N       sum over 1-N
!
!  var x = sqrt ( sum(x(i) -x_mean)**2)/N  )  = sum(x(i)*x(i))/N  - mean_x**2
!
!-------------------------------------------------------------------------------
! Regular grid where nvalues depends on resolution so constant

nvalues = REAL(new_res*new_res)
rnvalues = 1.0/nvalues

!$OMP PARALLEL DO PRIVATE(i,j,ii,jj,iii,jjj) DEFAULT(NONE)                 &
!$OMP& SHARED(num_x, num_y, new_res, rnvalues, mean_out, var_out, field_in)

DO jj=1,num_y
  DO ii= 1,num_x
    mean_out(ii,jj) = 0.0
    var_out(ii,jj)  = 0.0

    DO j=1,new_res
      jjj = (jj-1)*new_res+j
      DO i=1,new_res
        iii = (ii-1)*new_res+i
        mean_out(ii,jj)= mean_out(ii,jj)+field_in(iii,jjj)
      END DO
    END DO

    mean_out(ii,jj) = mean_out(ii,jj)*rnvalues

  END DO
END DO

!$OMP END PARALLEL DO

! problem with x**2 method as very big numbers and rounding issues

!$OMP PARALLEL DO PRIVATE(i,j,ii,jj,iii,jjj) DEFAULT(NONE)                 &
!$OMP& SHARED(num_x, num_y, new_res, rnvalues, mean_out, var_out, field_in)

DO jj = 1,num_y
  DO ii = 1,num_x

    DO j=1,new_res
      jjj = (jj-1)*new_res+j
      DO i=1,new_res
        iii = (ii-1)*new_res+i
        var_out(ii,jj)= var_out(ii,jj) + ( field_in(iii,jjj)-mean_out(ii,jj) ) &
                                        *( field_in(iii,jjj)-mean_out(ii,jj) )
      END DO
    END DO
    IF (var_out(ii,jj) > 0.0) THEN
      var_out(ii,jj) = SQRT(var_out(ii,jj)*rnvalues)
    ELSE
      var_out(ii,jj) = 0.0
    END IF

  END DO
END DO
!$OMP END PARALLEL DO


!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE mean_f2d
END MODULE mean_f2d_mod
