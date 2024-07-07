! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Tri_Linear

SUBROUTINE eg_Tri_Linear(                                         &
                      Ext_Data,                                   &
                      dim_i_in, dim_j_in, dim_k_in,               &
                      dim_i_out, dim_j_out, dim_k_out,            &
                      halo_i, halo_j, number_of_inputs,           &
                      weight_lambda, weight_phi,                  &
                      i_out, j_out, k_out,                        &
                      coeff_z_lin,                                &
                      Data_out)

! Purpose:
!          Performs linear interpolation of the input field to a
!          set of points given by weight_lambda,weight_phi,r_out.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE um_types, ONLY: integer32
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.


INTEGER ::                                                        &
  dim_i_in                                                        &
              ! Dimension of Data_in in i direction.
, dim_j_in                                                        &
              ! Dimension of Data_in in j direction.
, dim_k_in                                                        &
              ! Dimension of Data_in in k direction.
, dim_i_out                                                       &
              ! Dimension of Data_out in i direction.
, dim_j_out                                                       &
              ! Dimension of Data_out in j direction.
, dim_k_out                                                       &
              ! Dimension of Data_out in k direction.
, halo_i                                                          &
              ! Size of halo in i direction.
, halo_j                                                          &
              ! Size of halo in j direction.
, number_of_inputs ! number of fields to interpolate

REAL ::                                                           &
  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
            1-halo_j:dim_j_in+halo_j,                             &
                                             !data to be
           -1:dim_k_in+2, number_of_inputs)                       &
                                             ! interpolated
, coeff_z_lin(dim_i_out*dim_j_out*dim_k_out, -2:3)              &
, weight_lambda (dim_i_out*dim_j_out*dim_k_out)                 &
                                                  ! a number
                                                ! between 0 & 1
, weight_phi (dim_i_out*dim_j_out*dim_k_out)    ! a number between
                                                ! 0 & 1

INTEGER(KIND=integer32) ::                                        &
  i_out (dim_i_out*dim_j_out*dim_k_out)                         &
                                                ! point such that
, j_out (dim_i_out*dim_j_out*dim_k_out)                         &
                                                ! the desired
                                                ! output point
, k_out (dim_i_out*dim_j_out*dim_k_out)         ! lies between it
                                                ! and it+1

! Arguments with Intent OUT. ie: Output variables.

REAL ::                                                           &
  Data_out (dim_i_out*dim_j_out*dim_k_out, number_of_inputs)
                ! data interpolated to desired locations.

! Local Variables.

INTEGER ::                                                        &
  i, j, k, n , ii, jj, kk, ind  ! Loop indices

INTEGER :: dim_out

REAL ::                                                           &
  val_here                                                        &
, val_here_plus                                                   &
, a_coeff, b_coeff, c_coeff, d_coeff


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TRI_LINEAR'


! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Perform Linear Interpolation horizontal directions
!              simultaneously.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

dim_out = dim_i_out*dim_j_out*dim_k_out

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(n,ind,val_here,val_here_plus)&
!$OMP& PRIVATE(a_coeff,b_coeff,c_coeff,d_coeff,ii,jj,kk)                &
!$OMP& SHARED(dim_out,weight_lambda,weight_phi,                         &
!$OMP&        number_of_inputs,i_out,j_out,k_out,Ext_data,coeff_z_lin,  &
!$OMP&        data_out) 

!DIR$ NOBLOCKING
DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC) 
!DIR$ VECTOR ALWAYS
  DO ind = 1, dim_out
    d_coeff = weight_lambda(ind)*weight_phi(ind)
    c_coeff = weight_phi(ind) - d_coeff
    b_coeff = weight_lambda(ind) - d_coeff 
    a_coeff = 1.0-weight_lambda(ind)-c_coeff

    ii=i_out(ind)
    jj=j_out(ind)
    kk=k_out(ind)

    val_here = a_coeff   * Ext_Data (ii,jj,kk,n)         &
              +  b_coeff * Ext_Data (ii+1,jj,kk,n)       &
              +  c_coeff * Ext_Data (ii,jj+1,kk,n)       &
              +  d_coeff * Ext_Data (ii+1,jj+1,kk,n)

    val_here_plus = a_coeff * Ext_Data (ii,jj,kk+1,n)    &
              +  b_coeff    * Ext_Data (ii+1,jj,kk+1,n)  &
              +  c_coeff    * Ext_Data (ii,jj+1,kk+1,n)  &
              +  d_coeff    * Ext_Data (ii+1,jj+1,kk+1,n)

    ! interpolate data vertically
    Data_out (ind,n) = coeff_z_lin(ind,1)                &
                   * val_here_plus                       &
                   + coeff_z_lin(ind,0) * val_here

  END DO 
!$OMP END  DO NOWAIT

END DO
!$OMP END PARALLEL

! End of routine.
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_Tri_Linear
