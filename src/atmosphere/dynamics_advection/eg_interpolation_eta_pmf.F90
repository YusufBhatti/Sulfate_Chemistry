! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_interpolation_eta_pmf_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE ::  &
                  ModuleName='EG_INTERPOLATION_ETA_PMF_MOD'

CONTAINS
!
!=============================================================================
!
!
SUBROUTINE eg_interpolation_eta_pmf(                                         &
    eta_in,  pnt_type,                                                       &
    number_of_inputs,                                                        &
    dim_i_in, dim_j_in, dim_k_in,                                            &
    dim_j_in_w,                                                              &
    dim_i_out, dim_j_out, dim_k_out,                                         &
    high_order_scheme, monotone_scheme,                                      &
    l_high, l_mono,                                                          &
    eta_out, lambda_out, phi_out,                                            &
    me, n_proc, n_procx, n_procy,                                            &
    halo_i, halo_j, g_row_length,                                            &
    datastart, at_extremity, g_i_pe,                                         &
    proc_row_group, proc_col_group,                                          &
    halo_data_out_i, halo_data_out_j,                                        &
    error_code,                                                              &
    data_in, data_out,                                                       &
    k_int_linear_in)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE enforce_mono_mz_mod
USE sl_input_mod, ONLY: pmf_identifier, spmf_identifier 
USE eg_interpolation_eta_mod, ONLY: eg_interpolation_eta

IMPLICIT NONE
!
! Description:
!
!           This routine is wraper for eg_interpolation_eta so different
!           options for monotonicity (including PMF) can be used
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INTERPOLATION_ETA_PMF'


! Arguments with Intent IN. ie: Input variables.
            
INTEGER ::  dim_i_in, dim_j_in, dim_j_in_w, dim_k_in, dim_i_out,dim_j_out, & 
            dim_k_out, me, n_proc, n_procx, n_procy, halo_i, halo_j,       &
            halo_data_out_i, halo_data_out_j, proc_row_group,              &
            proc_col_group, g_row_length, datastart(3),                    &    
            g_i_pe(1-halo_i:g_row_length+halo_i) 

INTEGER ::  pnt_type, high_order_scheme, monotone_scheme, number_of_inputs
INTEGER, OPTIONAL :: k_int_linear_in
LOGICAL ::  l_high, l_mono
REAL    ::  eta_in(dim_k_in) 
REAL    ::  data_in(1-halo_i:dim_i_in+halo_i,  &
                    1-halo_j:dim_j_in+halo_j,  &
                    dim_k_in, number_of_inputs )

REAL ::  lambda_out (dim_i_out, dim_j_out, dim_k_out),        &
            phi_out (dim_i_out, dim_j_out, dim_k_out),        &
            eta_out (dim_i_out, dim_j_out, dim_k_out)

LOGICAL ::   at_extremity(4)

! Arguments with Intent OUT. ie: Output variables.

REAL ::  data_out (1-halo_data_out_i:dim_i_out+halo_data_out_i,           &
                   1-halo_data_out_j:dim_j_out+halo_data_out_j,           &
                   dim_k_out, number_of_inputs)
REAL ::  data_out_linear (1-halo_data_out_i:dim_i_out+halo_data_out_i,    &
                          1-halo_data_out_j:dim_j_out+halo_data_out_j,    &
                          dim_k_out, number_of_inputs)
INTEGER :: error_code, k_int_linear
LOGICAL :: pmf_stringent

! local variables
LOGICAL ::  l_high_in, l_mono_in

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (PRESENT(k_int_linear_in)) THEN
  k_int_linear=k_int_linear_in
ELSE
  k_int_linear=1
END IF

IF ( (l_mono .AND. l_high) .AND.                                 &
                      ((monotone_scheme == pmf_identifier ) .OR. &
                       (monotone_scheme == spmf_identifier))     ) THEN
    !
    ! if monotone_scheme == pmf/spmf then a high-order solution is needed
    ! obtained by the first call of "eg_interpolation_eta" with the same
    ! arguments passed through this subroutine and a lower-order
    ! solution (i.e., linear) using the second call with specific
    ! arguments to return a tri-linear solution.
    !
    ! PMF scheme uses these two solutions and product of slopes to decide
    ! which point, the high-order solution is oscillatory and replace it
    ! locally with the low-order one.
    ! 
     
    pmf_stringent = .FALSE.
    IF ( monotone_scheme == spmf_identifier ) THEN
        pmf_stringent = .TRUE.
    END IF
    
     ! This returns a high-order solution
     l_high_in = .TRUE.
     l_mono_in = .FALSE.
     CALL eg_interpolation_eta(                                         &
                               eta_in,  pnt_type,                       &
                               number_of_inputs,                        &
                               dim_i_in, dim_j_in, dim_k_in,            &
                               dim_j_in_w,                              &
                               dim_i_out, dim_j_out, dim_k_out,         &
                               high_order_scheme, 1,                    &
                               l_high_in, l_mono_in,                    &
                               eta_out, lambda_out, phi_out,            &
                               me, n_proc, n_procx, n_procy,            &
                               halo_i, halo_j, g_row_length,            &
                               datastart, at_extremity, g_i_pe,         &
                               proc_row_group, proc_col_group,          &
                               halo_data_out_i, halo_data_out_j,        &
                               error_code,                              &
                               data_in, data_out,                       &
                               k_int_linear                             )
     
      ! This to return a the low-order (linear) solution
     l_high_in = .FALSE.
     l_mono_in = .TRUE.     
     CALL eg_interpolation_eta(                                         &
                               eta_in,  pnt_type,                       &
                               number_of_inputs,                        &
                               dim_i_in, dim_j_in, dim_k_in,            &
                               dim_j_in_w,                              &
                               dim_i_out, dim_j_out, dim_k_out,         &
                               0, 1,                                    &
                               l_high_in, l_mono_in,                    &
                               eta_out, lambda_out, phi_out,            &
                               me, n_proc, n_procx, n_procy,            &
                               halo_i, halo_j, g_row_length,            &
                               datastart, at_extremity, g_i_pe,         &
                               proc_row_group, proc_col_group,          &
                               halo_data_out_i, halo_data_out_j,        &
                               error_code,                              &
                               data_in, data_out_linear,                &
                               k_int_linear                             )


   ! Enforce monotonicity using the PMF scheme
    
    CALL enforce_mono_PMF(data_out, data_out_linear,                    &
                          1-halo_data_out_i,dim_i_out+halo_data_out_i,  &
                          1-halo_data_out_j,dim_j_out+halo_data_out_j,  &
                          1, dim_k_out, 1, dim_i_out, 1, dim_j_out,     &
                          1,dim_k_out, number_of_inputs, pmf_stringent  )

  ELSE 
 
    ! else use the original eg_interpolation_eta as it is
    ! when monotone_scheme = 0, 1
  
    CALL eg_interpolation_eta(                                          &
                               eta_in,  pnt_type,                       &
                               number_of_inputs,                        &
                               dim_i_in, dim_j_in, dim_k_in,            &
                               dim_j_in_w,                              &
                               dim_i_out, dim_j_out, dim_k_out,         &
                               high_order_scheme, monotone_scheme,      &
                               l_high, l_mono,                          &
                               eta_out, lambda_out, phi_out,            &
                               me, n_proc, n_procx, n_procy,            &
                               halo_i, halo_j, g_row_length,            &
                               datastart, at_extremity, g_i_pe,         &
                               proc_row_group, proc_col_group,          &
                               halo_data_out_i, halo_data_out_j,        &
                               error_code,                              &
                               data_in, data_out,                       &
                               k_int_linear                             )
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE eg_interpolation_eta_pmf
END MODULE eg_interpolation_eta_pmf_mod
