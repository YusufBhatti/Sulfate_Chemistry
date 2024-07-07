! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_correct_thetav_priestley_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='EG_CORRECT_THETAV_PRIESTLEY_MOD'

CONTAINS

SUBROUTINE eg_correct_thetav_priestley(g_i_pe, rho, rho_np1,            &
                                       depart_lambda, depart_phi,       &
                                       depart_r, thetav, s_thetav,      &
                                       thetav_np1)


USE um_parvars,     ONLY: nproc_x, nproc_y,                             &
                          at_extremity,gc_proc_row_group,               &
                          gc_proc_col_group, datastart
USE um_parcore,     ONLY: mype, nproc

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

USE level_heights_mod,       ONLY: eta_theta_levels, eta_rho_levels      
USE atm_fields_bounds_mod,   ONLY: pdims, pdims_s,                      &
                                   wdims, tdims, tdims_s, tdims_l
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE ref_pro_mod
USE metric_terms_mod
USE Field_Types
USE priestley_algorithm_mod, ONLY: priestley_algorithm3
USE eg_helmholtz_mod,        ONLY: ec_vol
USE timestep_mod,            ONLY: timestep_number, timestep
USE parkind1,                ONLY: jpim, jprb       !DrHook
USE yomhook,                 ONLY: lhook, dr_hook   !DrHook  

IMPLICIT NONE             
!
! Description: Apply the Priestley conservation algorithm to
!              potential temperature
!
! Method:
!
! The Priestley algorithm (Monthly Weather Review, 121, 621--629)
! is applied to enforce conservation of dryrho*thetav in ENDGame.
! Proper account is taken of the vertical Charney-Phillips 
! grid-staggering. See UMDP16.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) :: g_i_pe(1-(tdims%i_start-tdims_l%i_start):        &
                              global_row_length +                       &
                              (tdims%i_start-tdims_l%i_start) )

REAL, INTENT(IN)    :: rho(pdims_s%i_start:pdims_s%i_end,               &
                           pdims_s%j_start:pdims_s%j_end,               &
                           pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN)    :: rho_np1(pdims_s%i_start:pdims_s%i_end,           &
                               pdims_s%j_start:pdims_s%j_end,           &
                               pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN)    :: depart_lambda(wdims%i_start:wdims%i_end,         &
                                     wdims%j_start:wdims%j_end,         &
                                     wdims%k_start:wdims%k_end)

REAL, INTENT(IN)    :: depart_phi(wdims%i_start:wdims%i_end,            &
                                  wdims%j_start:wdims%j_end,            &
                                  wdims%k_start:wdims%k_end)

REAL, INTENT(IN)    :: depart_r(wdims%i_start:wdims%i_end,              &
                                wdims%j_start:wdims%j_end,              &
                                wdims%k_start:wdims%k_end)

REAL, INTENT(IN)    :: thetav(tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end,            &
                              tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)    :: s_thetav(tdims_s%i_start:tdims_s%i_end,          &
                                tdims_s%j_start:tdims_s%j_end,          &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: thetav_np1(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

! Local
LOGICAL :: l_high_in   ! local setting for use in CALLs
LOGICAL :: l_mono_in   ! local setting for use in CALLs

INTEGER             :: i, j, k, error_code

REAL                :: corr_factor

REAL                :: psi(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,                   &
                           tdims%k_start:tdims%k_end)

REAL                :: psi_np1(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                               tdims%k_start:tdims%k_end)

REAL                :: alfa_za(pdims%k_start:pdims%k_end)
REAL                :: beta_za(pdims%k_start:pdims%k_end)

REAL                :: ext_data(tdims_l%i_start:tdims_l%i_end,          &
                                tdims_l%j_start:tdims_l%j_end,          &
                                tdims_l%k_start:tdims_l%k_end)

REAL                :: thetav_np1_low(tdims_s%i_start:tdims_s%i_end,    &
                                      tdims_s%j_start:tdims_s%j_end,    &
                                      tdims_s%k_start:tdims_s%k_end)

INTEGER, PARAMETER :: no_high_order_scheme=0
INTEGER, PARAMETER :: tri_linear=1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_THETAV_PRIESTLEY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Vertical averaging weights
DO k = pdims%k_start, pdims%k_end
  alfa_za(k) = (eta_rho_levels(k) - eta_theta_levels(k-1)) /            &
               (eta_theta_levels(k) - eta_theta_levels(k-1))
  beta_za(k) = 1.0 - alfa_za(k)      
END DO
alfa_za(pdims%k_start) = 1.0
beta_za(pdims%k_start) = 0.0

! Calculate integration weights
! psi(i,j,k) = rho*volume*average
! so mass_tracer = SUM[psi*tracer] = SUM[av(tracer)*rho*vol]
psi    (:,:,tdims%k_start) = 0.0
psi_np1(:,:,tdims%k_start) = 0.0

DO k = tdims%k_start + 1, tdims%k_end - 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi(i,j,k)    =  rho(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )      &
                    +  rho(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
      
      psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                     + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
    END DO
  END DO
END DO

k = tdims%k_end
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    psi(i,j,k)     =   rho(i,j,k) * ec_vol(i,j,k) * alfa_za(k)   
    psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
  END DO
END DO

! Copy time t_n solution to extended array; remove
! atmos_physics2 sources from time t_(n+1) solution
DO k=tdims_s%k_start, tdims_s%k_end
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      ext_data(i,j,k)=thetav(i,j,k)
      thetav_np1(i,j,k)=thetav_np1(i,j,k) - s_thetav(i,j,k)
    END DO
  END DO
END DO  

CALL Swap_Bounds(ext_data, tdims%i_len, tdims%j_len, tdims%k_len,       &
                 tdims_l%halo_i, tdims_l%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

l_high_in=.FALSE.
l_mono_in=.TRUE.

! Apply linear interpolation to get low-order solution at time t_(n+1)
CALL eg_interpolation_eta_pmf(                                          &
                     eta_theta_levels,fld_type_w, 1,                    &
                     tdims%i_len, tdims%j_len, tdims%k_len, tdims%j_len,&
                     tdims%i_len, tdims%j_len, tdims%k_len,             &
                     no_high_order_scheme, tri_linear,                  &
                     l_high_in,l_mono_in,                               &
                     depart_r, depart_lambda, depart_phi,               &
                     mype, nproc, nproc_x, nproc_y,                     &
                     tdims_l%halo_i, tdims_l%halo_j,                    &
                     global_row_length, datastart, at_extremity,        &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                     tdims_s%halo_i, tdims_s%halo_j, error_code,        &
                     ext_data, thetav_np1_low)

! Priestley algorithm adjusts thetav_np1 to ensure conservation:
!   sum(thetav_np1 * psi_np1) = sum(thetav * psi)
CALL priestley_algorithm3(tdims, tdims_s, 1, psi, psi_np1,              &
                          thetav, thetav_np1_low, thetav_np1)

! Add back the atmos_physics2 source
DO k=tdims_s%k_start, tdims_s%k_end
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      thetav_np1(i,j,k)=thetav_np1(i,j,k) + s_thetav(i,j,k)
    END DO
  END DO
END DO

CALL Swap_Bounds(thetav_np1, tdims%i_len, tdims%j_len, tdims%k_len,     &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_correct_thetav_priestley

END MODULE eg_correct_thetav_priestley_mod
