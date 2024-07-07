! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_Cubic_Lagrange
!
SUBROUTINE eg_Cubic_Lagrange(                                            &
                          fld,                                           &
                          dim_i_in, dim_j_in, dim_k_in,                  &
                          dim_i_out, dim_j_out, dim_k_out,               &
                          halo_i, halo_j, number_of_inputs,              &
                          weight_lambda, weight_phi,                     &
                          s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,      &
                          i_out, j_out, k_out,                           &
                          coeff_z, lev_ext,                              &
                          data_out)

! Purpose:
!          Performs cubic Lagrange interpolation of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
!
! Method:
!          ENDGame formulation version 3.01
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

USE um_types, ONLY: integer32
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

! Array dimensions
INTEGER, INTENT(IN) :: dim_i_in     ! X dimension of input data
INTEGER, INTENT(IN) :: dim_j_in     ! Y dimension of input data
INTEGER, INTENT(IN) :: dim_k_in     ! Z dimension of input data
INTEGER, INTENT(IN) :: dim_i_out    ! X dimension of output data
INTEGER, INTENT(IN) :: dim_j_out    ! Y dimension of output data
INTEGER, INTENT(IN) :: dim_k_out    ! Z dimension of output data
INTEGER, INTENT(IN) :: halo_i       ! Halo size X
INTEGER, INTENT(IN) :: halo_j       ! Halo size Y
INTEGER, INTENT(IN) :: number_of_inputs    ! last dimension for input/output

INTEGER, INTENT(IN) :: lev_ext      ! vertical levels to loop over in coeff_z

! Input field
REAL, INTENT(IN) :: fld(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j, &
                        -1:dim_k_in+2, number_of_inputs)

! Interpolation weights
REAL, INTENT(IN) :: weight_lambda (dim_i_out*dim_j_out*dim_k_out)
REAL, INTENT(IN) :: weight_phi (dim_i_out*dim_j_out*dim_k_out)

REAL, INTENT(IN) :: s_xi1(2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: s_xi2(2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: t_xi1(2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: t_xi2(2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: q_xi1(-1:2,2-halo_i:dim_i_in+halo_i-2)
REAL, INTENT(IN) :: q_xi2(-1:2,2-halo_j:dim_j_in+halo_j-2)
REAL, INTENT(IN) :: coeff_z(dim_i_out*dim_j_out*dim_k_out,-2:3)

INTEGER(KIND=integer32), INTENT(IN) :: i_out (dim_i_out*dim_j_out*dim_k_out)
INTEGER(KIND=integer32), INTENT(IN) :: j_out (dim_i_out*dim_j_out*dim_k_out)
INTEGER(KIND=integer32), INTENT(IN) :: k_out (dim_i_out*dim_j_out*dim_k_out)

! Arguments with Intent OUT. ie: Output variables.
REAL, INTENT(OUT) :: data_out(dim_i_out*dim_j_out*dim_k_out, number_of_inputs)

! Local Variables.
! Loopers and indexers
INTEGER  :: ind, ko,kk,ll, n, im1,i0,ip1,ip2, jm1,j0,jp1,jp2

REAL :: x, s, t, cij      ! temporaries

REAL :: val_i_minus
REAL :: val_i
REAL :: val_i_plus
REAL :: val_i_plus2
REAL :: coeff_i_minus
REAL :: coeff_i_zero
REAL :: coeff_i_plus
REAL :: coeff_i_plus2
REAL :: coeff_j_minus
REAL :: coeff_j_zero
REAL :: coeff_j_plus
REAL :: coeff_j_plus2

INTEGER :: dim_out

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CUBIC_LAGRANGE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

dim_out = dim_i_out*dim_j_out*dim_k_out

! Interpolate data to get the data at the point for the k interpolation
!! rest of benchmark changes change results
!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP& PRIVATE(ind,n,val_i_plus2,val_i_plus,val_i,val_i_minus,                 &
!$OMP& im1, i0, ip1, ip2, jm1, j0, jp1, jp2, ko, kk, ll,                       &
!$OMP& coeff_j_plus2,coeff_j_plus,coeff_j_zero,coeff_j_minus,                  &
!$OMP& coeff_i_plus2,coeff_i_plus,coeff_i_zero,coeff_i_minus,x,s,t,cij)
!DIR$ NOBLOCKING
DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC) 
!DIR$ VECTOR ALWAYS
  DO ind = 1, dim_out
    i0  = i_out(ind)
    im1 = i0-1
    ip1 = i0+1
    ip2 = i0+2

    j0  = j_out(ind)
    jm1 = j0-1
    jp1 = j0+1
    jp2 = j0+2

    x             = weight_phi(ind)
    s             = s_xi2(j0)
    t             = t_xi2(j0)
    
    coeff_j_minus = x*(x-1.0)*(x-t)    *q_xi2(-1,j0)
    coeff_j_zero  = (x-s)*(x-1.0)*(x-t)*q_xi2(0,j0)
    coeff_j_plus  = (x-s)*x*(x-t)      *q_xi2(1,j0)
    coeff_j_plus2 = (x-s)*x*(x-1.0)    *q_xi2(2,j0)


    x             = weight_lambda(ind)
    s             = s_xi1(i0)
    t             = t_xi1(i0)
    
    coeff_i_minus = x*(x-1.0)*(x-t)    *q_xi1(-1,i0)
    coeff_i_zero  = (x-s)*(x-1.0)*(x-t)*q_xi1(0,i0)
    coeff_i_plus  = (x-s)*x*(x-t)      *q_xi1(1,i0)
    coeff_i_plus2 = (x-s)*x*(x-1.0)    *q_xi1(2,i0)
    
    ko  = k_out(ind)

    data_out(ind,n) = 0.0
!DIR$ loop_info min_trips(2) max_trips(6)
    DO ll = -lev_ext, lev_ext+1
      kk = ko + ll

      val_i_minus = coeff_j_plus2*fld(im1,jp2,kk,n)                &
                   +coeff_j_plus*fld(im1,jp1,kk,n)                 &
                   +coeff_j_zero*fld(im1,j0,kk,n)                  &
                   +coeff_j_minus*fld(im1,jm1,kk,n)
      
      val_i       = coeff_j_plus2*fld(i0,jp2,kk,n)                 &
                   +coeff_j_plus*fld(i0,jp1,kk,n)                  &
                   +coeff_j_zero*fld(i0,j0,kk,n)                   &
                   +coeff_j_minus*fld(i0,jm1,kk,n)
      
      val_i_plus  = coeff_j_plus2*fld(ip1,jp2,kk,n)                &
                   +coeff_j_plus*fld(ip1,jp1,kk,n)                 &
                   +coeff_j_zero*fld(ip1,j0,kk,n)                  &
                   +coeff_j_minus*fld(ip1,jm1,kk,n)

      val_i_plus2 = coeff_j_plus2*fld(ip2,jp2,kk,n)                &
                   +coeff_j_plus*fld(ip2,jp1,kk,n)                 &
                   +coeff_j_zero*fld(ip2,j0,kk,n)                  &
                   +coeff_j_minus*fld(ip2,jm1,kk,n)

      cij         = coeff_i_plus2*val_i_plus2                      &
                   +coeff_i_plus*val_i_plus                        &
                   +coeff_i_zero*val_i                             &
                   +coeff_i_minus*val_i_minus

      data_out(ind,n) = data_out(ind,n)+ coeff_z(ind,ll)*cij
    END DO  ! ll
  END DO ! ind
!$OMP END DO NOWAIT
END DO ! n
!$OMP END PARALLEL


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_Cubic_Lagrange
