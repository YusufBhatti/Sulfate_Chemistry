! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_thermo_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_THERMO_MOD'

CONTAINS
SUBROUTINE eg_sl_thermo(g_i_pe,                                       &
                l_inc_solver,  high_order_scheme,                     &
                monotone_scheme, l_high, l_mono,                      &
                r_theta,  etadot, r_theta_d,                          &
                rho_n, rho_np1, error_code                            )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_parcore,        ONLY: mype, nproc

USE um_parvars,        ONLY: nproc_x,nproc_y,                         &
                             at_extremity,gc_proc_row_group,           &
                             gc_proc_col_group,                        &
                             offx,offy,halo_i,halo_j, datastart

USE nlsizes_namelist_mod, ONLY: global_row_length,                     &
                                row_length, rows, n_rows, model_levels


USE timestep_mod,      ONLY: timestep
USE eg_alpha_mod,      ONLY: tau_theta
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,         &
                              xi3_at_theta=>r_theta_levels,            &
                              xi3_at_rho=>r_rho_levels,                &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE atm_fields_bounds_mod
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE departure_pts_mod
USE Field_Types
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE atm_step_local, ONLY: L_print_L2norms, cycleno
USE sl_thermo_norm_mod
USE lam_config_inputs_mod, ONLY: n_rims_to_do
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE umPrintMgr, ONLY: umMessage, umPrint
USE eg_zlf_conservation_mod, ONLY: eg_zlf_conservation
USE eg_zlf_mod, ONLY: zlf_cfl_top_level_theta
USE dynamics_input_mod, ONLY: l_sl_bc_correction,              &
                              zlf_conservation_theta_option

IMPLICIT NONE
!
! Description:
!  Find timelevel n dependent quantity R_theta_d
!
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_THERMO'

LOGICAL, INTENT(IN) :: l_inc_solver

INTEGER, INTENT(IN) ::                                                &
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction


! Loop index bounds for arrays defined on p, u, v points respectively


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono
                   ! True, if interpolation required to be monotone.


! Timelevel n arrival point quantities

REAL, INTENT(IN) ::                                                   &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  r_theta(1-offx:row_length+offx,1-offy:rows+offy,                    &
          0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.
REAL, INTENT(IN) ::   rho_n(pdims_s%i_start:pdims_s%i_end,            &
                            pdims_s%j_start:pdims_s%j_end,            &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN) :: rho_np1(pdims_s%i_start:pdims_s%i_end,            &
                            pdims_s%j_start:pdims_s%j_end,            &
                            pdims_s%k_start:pdims_s%k_end)


! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_theta_d(1-offx:row_length+offx,1-offy:rows+offy,                  &
            0:model_levels)


! Local variables

INTEGER :: i,j,k, number_of_inputs
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

! tmp arrays

REAL :: rdz
INTEGER :: number_of_inputs_zlf


REAL :: work(1-halo_i:row_length+halo_i,                              &
              1-halo_j:rows+halo_j,0:model_levels)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO SHARED(work,r_theta,pdims,model_levels),              &
!$OMP             PRIVATE(i,j,k) SCHEDULE(STATIC) DEFAULT(NONE)
DO k=0, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      work(i,j,k) = r_theta(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL swap_bounds(work,                                                 &
                 tdims_l%i_len - 2*tdims_l%halo_i,                     &
                 tdims_l%j_len - 2*tdims_l%halo_j,                     &
                 tdims_l%k_len,                                        &
                 tdims_l%halo_i, tdims_l%halo_j,                       &
                 fld_type_p,swap_field_is_scalar)


number_of_inputs = 1

! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

IF ( zlf_conservation_theta_option > 0 ) THEN  

  number_of_inputs_zlf   = number_of_inputs
     
  CALL eg_zlf_conservation(work, r_theta_d,                         &
           number_of_inputs, number_of_inputs_zlf,                  &
           row_length, rows, n_rows, model_levels, halo_i, halo_j,  &
           offx, offy, datastart, g_i_pe, high_order_scheme,        &
           monotone_scheme,  l_high, l_mono,                        &
           zlf_conservation_theta_option, zlf_cfl_top_level_theta,  &
           rho_n, rho_np1, error_code                               )
      
                         
ELSE ! use the old non-conservative advection

  CALL eg_interpolation_eta_pmf(                                    &
                   eta_theta_levels,fld_type_w,                     &
                   number_of_inputs,                                &
                   row_length, rows, model_levels+1, rows,          &
                   row_length, rows, model_levels+1,                &
                   high_order_scheme, monotone_scheme,              &
                   l_high, l_mono, depart_xi3_w, depart_xi1_w,      &
                   depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                   halo_i, halo_j,                                  &
                   global_row_length, datastart, at_extremity,      &
                   g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                   offx, offy, error_code,                          &
                   work, r_theta_d, k_int_linear_in=k_int_linear)

END IF ! endif for use ZLF

! Compute thetav_ref_eta

IF ( .NOT. l_inc_solver ) THEN
  ! NOTE: etadot = 0 at k=0 and k= model_levels.
  !       R_theta_d does not need updating there.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,rdz)                       &
!$OMP          SHARED(model_levels, pdims, eta_rho_levels, r_theta_d, &
!$OMP                 thetav_ref_pro, tau_theta, timestep, etadot,    &
!$OMP                 thetav_ref_eta)
!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels-1
    rdz = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k) )
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_theta_d(i,j,k) = r_theta_d(i,j,k)                           &
                          -thetav_ref_pro(i,j,k)                      &
                          +tau_theta*timestep*etadot(i,j,k)           &
                          *( thetav_ref_eta(i,j,k+1)                  &
                            -thetav_ref_eta(i,j,k)                    &
                           )*rdz
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      k = 0
      r_theta_d(i,j,k) = r_theta_d(i,j,k) - thetav_ref_pro(i,j,k)
      k = model_levels
      r_theta_d(i,j,k) = r_theta_d(i,j,k) - thetav_ref_pro(i,j,k)
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
END IF

IF ( L_print_L2norms ) THEN
  WRITE(umMessage,'(A, I2, A)') ' ** cycleno =  ', cycleno              &
                              , ' **   L2 norms after eg_sl_thermo **'
  CALL umPrint(umMessage,src='EG_SL_THERMO')
  CALL sl_thermo_norm(                                                  &
                      norm_lev_start, norm_lev_end, n_rims_to_do,       &
                      r_theta, r_theta_d, .TRUE., .FALSE., .FALSE. )
END IF !  L_print_L2norms

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_thermo
END MODULE eg_sl_thermo_mod
