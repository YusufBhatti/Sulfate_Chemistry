! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Quintic_Lagrange

      SUBROUTINE eg_Quintic_Lagrange(                                   &
                                Ext_Data,                               &
                                dim_i_in, dim_j_in, dim_k_in,           &
                                dim_i_out, dim_j_out, dim_k_out,        &
                                halo_i, halo_j, number_of_inputs,       &
                                weight_lambda, weight_phi,              &
                                i_out, j_out, k_out,                    &
                                coeff_z,                                &
                                Data_out)

! Purpose:
!          Performs quintic Lagrange interpolation of the input field to
!          a set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
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
      USE um_types, ONLY: integer32
      USE yomhook, ONLY: lhook, dr_hook
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
      , weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
      , weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

      INTEGER(KIND=integer32) ::                                        &
        i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
      , j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
      , k_out (dim_i_out, dim_j_out, dim_k_out)       ! lies between it
                                                      ! and it+1

! Arguments with Intent IN
      REAL ::                                                           &
        coeff_z(dim_i_out, dim_j_out, dim_k_out, -2:3)

! Arguments with Intent OUT. ie: Output variables.

      REAL ::                                                           &
        Data_out (dim_i_out, dim_j_out, dim_k_out, number_of_inputs)
                      ! data interpolated to desired locations.

! Local Variables.

      INTEGER ::                                                        &
        i, j, k, n, INDEX ! Loop indices

      REAL ::                                                           &
        one_sixth                                                       &
      , one_twelve                                                      &
      , one_twentyfour                                                  &
      , one_hundredtwenty                                               &
      , val_i_minus                                                     &
      , val_i_minus2                                                    &
      , val_i                                                           &
      , val_i_plus                                                      &
      , val_i_plus2                                                     &
      , val_i_plus3                                                     &
      , phi                                                             &
      , phi_cubed                                                       &
      , phi_sq                                                          &
      , phi_fourth                                                      &
      , phi_fifth                                                       &
      , lambda_cubed                                                    &
      , lambda_sq                                                       &
      , lambda                                                          &
      , lambda_fourth                                                   &
      , lambda_fifth                                                    &
      , coeff_minus2(dim_i_out, dim_j_out, dim_k_out)                   &
      , coeff_minus(dim_i_out, dim_j_out, dim_k_out)                    &
      , coeff_zero(dim_i_out, dim_j_out, dim_k_out)                     &
      , coeff_plus(dim_i_out, dim_j_out, dim_k_out)                     &
      , coeff_plus2(dim_i_out, dim_j_out, dim_k_out)                    &
      , coeff_plus3(dim_i_out, dim_j_out, dim_k_out)                    &
      , coeffl_minus2(dim_i_out, dim_j_out, dim_k_out)                  &
      , coeffl_minus(dim_i_out, dim_j_out, dim_k_out)                   &
      , coeffl_zero(dim_i_out, dim_j_out, dim_k_out)                    &
      , coeffl_plus(dim_i_out, dim_j_out, dim_k_out)                    &
      , coeffl_plus2(dim_i_out, dim_j_out, dim_k_out)                   &
      , coeffl_plus3(dim_i_out, dim_j_out, dim_k_out)

      REAL ::                                                           &
        col_data (dim_i_out, dim_j_out, -2:3) ! Horizontally interpolatd
                                              ! data ready for vertical
                                              ! interpolation.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_QUINTIC_LAGRANGE'



! External Routines: None.

! NB: This routine has not been optimised unlike the other SL
!     interpolation routines.

! ----------------------------------------------------------------------
! Section 1.   Perform quintic Lagrange Interpolation in j direction and
!              then in i direction.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
      one_sixth = 1.0/6.0
      one_twelve = 1.0/12.0
      one_twentyfour = 1.0/24.0
      one_hundredtwenty = 1.0/120.0

      DO k = 1, dim_k_out
#if !defined(NEC)
        DO j = 1, dim_j_out
        DO i = 1, dim_i_out
#else
        j=1
        DO i = 1, dim_i_out * dim_j_out
#endif
          phi = weight_phi (i,j,k)
          phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
          phi_cubed = weight_phi (i,j,k) * phi_sq
          phi_fourth = phi_sq * phi_sq
          phi_fifth = phi_sq * phi_cubed

          coeff_plus3(i,j,k) = one_hundredtwenty * (phi_fifth -         &
                            5.0 * phi_cubed + 4.0 * phi)
          coeff_plus2(i,j,k) = - one_twentyfour * (phi_fifth -          &
                            phi_fourth - 7.0 * phi_cubed + phi_sq        &
                             + 6.0 * phi)
          coeff_plus(i,j,k)  = one_twelve * (phi_fifth                  &
                             - 2.0 * phi_fourth                          &
                             - 7.0 * phi_cubed + 8.0 * phi_sq +           &
                            12.0 * phi)
          coeff_zero(i,j,k)  = - one_twelve * (phi_fifth -              &
                            3.0 * phi_fourth - 5.0 * phi_cubed +          &
                            15.0 * phi_sq + 4.0 * phi -12.0 )
          coeff_minus(i,j,k) =  one_twentyfour * (phi_fifth -           &
                             4.0 * phi_fourth - phi_cubed +              &
                             16.0 * phi_sq - 12.0 * phi )
          coeff_minus2(i,j,k) = - one_hundredtwenty * (phi_fifth -      &
                              5.0 * phi_fourth + 5.0 * phi_cubed +        &
                              5.0 * phi_sq - 6.0 * phi )

          lambda = weight_lambda (i,j,k)
          lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
          lambda_cubed = weight_lambda (i,j,k) * lambda_sq
          lambda_fourth = lambda_sq * lambda_sq
          lambda_fifth = lambda_sq * lambda_cubed

          coeffl_plus3(i,j,k) = one_hundredtwenty * (lambda_fifth -     &
                             5.0 * lambda_cubed + 4.0 * lambda)
          coeffl_plus2(i,j,k) = - one_twentyfour * (lambda_fifth -      &
                             lambda_fourth - 7.0 * lambda_cubed +        &
                             lambda_sq + 6.0 * lambda )
          coeffl_plus(i,j,k)  = one_twelve * (lambda_fifth -            &
                             2.0 * lambda_fourth - 7.0 * lambda_cubed     &
                             +  8.0 * lambda_sq + 12.0 * lambda )
          coeffl_zero(i,j,k)  = - one_twelve * (lambda_fifth -          &
                             3.0 * lambda_fourth - 5.0 * lambda_cubed     &
                             + 15.0 * lambda_sq + 4.0 * lambda - 12.0 )
          coeffl_minus(i,j,k) =  one_twentyfour * (lambda_fifth -       &
                              4.0 * lambda_fourth - lambda_cubed         &
                              + 16.0 * lambda_sq  - 12.0 * lambda )
          coeffl_minus2(i,j,k) = -1.0 * one_hundredtwenty *              &
                              (lambda_fifth -                           &
                              5.0 * lambda_fourth + 5.0 * lambda_cubed    &
                              + 5.0 *lambda_sq - 6.0 * lambda )

        END DO ! i loop
#if !defined(NEC)
        END DO ! j loop
#endif
      END DO

! begin loop over number of inputs which ends at end of subroutine

      DO n = 1, number_of_inputs

! begin loop over levels which ends at end of subroutine

      DO k = 1, dim_k_out

        DO INDEX = -2, 3
#if !defined(NEC)
            DO j = 1, dim_j_out
            DO i = 1, dim_i_out
#else
            j=1
            DO i = 1, dim_i_out * dim_j_out
#endif
!CDIR UNROLL=6
! interpolate in j at each index level for the six points needed for
! the i interpolation.


              val_i_minus2 = coeff_plus3(i,j,k) *                       &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)&
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k),k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)&
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)

              val_i_minus = coeff_plus3(i,j,k) *                        &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)&
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)&
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)

              val_i       = coeff_plus3(i,j,k) *                        &
           Ext_Data (i_out(i,j,k),j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)  &
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)  &
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)  &
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+INDEX,n)    &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k),j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)

              val_i_plus  = coeff_plus3(i,j,k) *                        &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)&
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)&
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)

              val_i_plus2 = coeff_plus3(i,j,k) *                        &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)&
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)&
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)


              val_i_plus3 = coeff_plus3(i,j,k) *                        &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+3,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus2(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+2,k_out(i,j,k)+INDEX,n)&
                               + coeff_plus(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+1,k_out(i,j,k)+INDEX,n)&
                               + coeff_zero(i,j,k) *                    &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k),k_out(i,j,k)+INDEX,n)  &
                               + coeff_minus(i,j,k) *                   &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-1,k_out(i,j,k)+INDEX,n)&
                               + coeff_minus2(i,j,k) *                  &
           Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-2,k_out(i,j,k)+INDEX,n)


! interpolate in i and store.

              col_data (i,j,INDEX) = coeffl_plus3(i,j,k) * val_i_plus3  &
                                   + coeffl_plus2(i,j,k) * val_i_plus2  &
                                   + coeffl_plus(i,j,k) * val_i_plus    &
                                   + coeffl_zero(i,j,k) * val_i         &
                                   + coeffl_minus(i,j,k) * val_i_minus  &
                                   + coeffl_minus2(i,j,k) * val_i_minus2

            END DO     ! i loop
#if !defined(NEC)
            END DO     ! j loop
#endif
        END DO    ! index loop

! ----------------------------------------------------------------------
! Section 2.   Perform quintic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        DO j = 1, dim_j_out
          DO i = 1, dim_i_out
            Data_out (i,j,k,n) = coeff_z(i,j,k,-2) *                    &
                                 col_data (i,j,-2) +                    &
                                 coeff_z(i,j,k,-1) *                    &
                                 col_data (i,j,-1) +                    &
                                 coeff_z(i,j,k,0) *                     &
                                 col_data (i,j,0) +                     &
                                 coeff_z(i,j,k,1) *                     &
                                 col_data (i,j,1) +                     &
                                 coeff_z(i,j,k,2) *                     &
                                 col_data (i,j,2) +                     &
                                 coeff_z(i,j,k,3) *                     &
                                 col_data (i,j,3)

          END DO
        END DO

! End loop over levels.
      END DO

! End loop over number of inputs.
      END DO

! End of routine.
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_Quintic_Lagrange
