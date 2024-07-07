! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!
MODULE eg_correct_tracers_ukca_mod

! Description:
!  This routine carries out mass conservation of UKCA tracers
!  using a modified version of the Priestley (1993) scheme.
!
! Method: ENDGame formulation 4.00 (Feb 2014)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Tracer Advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

USE missing_data_mod,  ONLY: imdi

IMPLICIT NONE

INTEGER :: eg_tracer_ukca_start=imdi
INTEGER :: eg_tracer_ukca_end=imdi
                     ! Indices of UKCA tracers in the super arrays.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CORRECT_TRACERS_UKCA_MOD'

CONTAINS

SUBROUTINE eg_correct_tracers_ukca(                                     &
                            row_length, rows, model_levels,             &
                            halo_i, halo_j, offx, offy,                 &
                            datastart, g_i_pe,                          &
                            super_array_size,                           &
                            super_tracer_phys1, super_tracer_phys2,     &
                            tracer_ukca, rho_n, rho_np1,                &
                            depart_lambda, depart_phi, depart_r,        &
                            tr_ukca                                  )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook


USE um_parvars,     ONLY: nproc_x, nproc_y,                       &
                          at_extremity,gc_proc_row_group,         &
                          gc_proc_col_group
USE um_parcore,     ONLY: mype, nproc

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels
USE atm_fields_bounds_mod
USE sl_input_mod,            ONLY: tracer_sl, high_order_scheme, &
                                   monotone_scheme, l_high, l_mono
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE departure_pts_mod
USE Field_Types
USE priestley_algorithm_mod
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE eg_helmholtz_mod,      ONLY: ec_vol
USE ukca_option_mod,   ONLY: i_ukca_conserve_method,                   &
                             i_ukca_hiorder_scheme,                    &
                             L_ukca_src_in_conservation,               &
                             ukca_conserve_um, priestley_old,          &
                             priestley_optimal, ukca_no_conserve

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length, rows, model_levels,                &
                       halo_i, halo_j, offx, offy, tr_ukca,           &
                       super_array_size

INTEGER, INTENT(IN) :: datastart(3)
INTEGER, INTENT(IN) :: g_i_pe(1-halo_i:global_row_length+halo_i)

REAL, INTENT(IN) :: super_tracer_phys1                                &
                     (tdims_l%i_start:tdims_l%i_end,                  &
                      tdims_l%j_start:tdims_l%j_end,                  &
                      tdims_l%k_start:tdims_l%k_end,                  &
                      super_array_size)

REAL, INTENT(IN) :: super_tracer_phys2                                &
                     (tdims%i_start:tdims%i_end,                      &
                      tdims%j_start:tdims%j_end,                      &
                      tdims%k_start:tdims%k_end,                      &
                      super_array_size)

REAL, INTENT(IN) :: rho_n                                            &
                      ( pdims_s%i_start:pdims_s%i_end,               &
                       pdims_s%j_start:pdims_s%j_end,                &
                       pdims_s%k_start:pdims_s%k_end),               &
!
                    rho_np1                                          &
                      ( pdims_s%i_start:pdims_s%i_end,               &
                       pdims_s%j_start:pdims_s%j_end,                &
                       pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN) ::  depart_lambda ( wdims%i_start:wdims%i_end,      &
                       wdims%j_start:wdims%j_end,                    &
                       wdims%k_start:wdims%k_end),                   &
!
                     depart_phi ( wdims%i_start:wdims%i_end,         &
                       wdims%j_start:wdims%j_end,                    &
                       wdims%k_start:wdims%k_end),                   &
!
                     depart_r  ( wdims%i_start:wdims%i_end,          &
                       wdims%j_start:wdims%j_end,                    &
                       wdims%k_start:wdims%k_end)
!
!in/out

REAL, INTENT(INOUT)  ::  tracer_ukca (                                &
                           tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end, tr_ukca)
! locals

REAL                 :: tracer_ukca_n (                               &
                           tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,             &
                           tr_ukca)

REAL                 :: tracer_ukca_np1_low (                         &
                           tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,             &
                           tr_ukca)

REAL                 :: tracer_ukca_sources (                         &
                           tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,             &
                           tr_ukca)

REAL                 :: tracer_ukca_sources2 (                        &
                           tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,             &
                           tr_ukca)

REAL                 :: tr_temp1 (                                    &
                           tdims_l%i_start:tdims_l%i_end,             &
                           tdims_l%j_start:tdims_l%j_end,             &
                           tdims_l%k_start:tdims_l%k_end,             &
                           tr_ukca)

! locals

LOGICAL :: l_high_in   ! local setting for use in CALLs
LOGICAL :: l_mono_in   ! local setting for use in CALLs

INTEGER   ::  error_code, i,j,k,kk,kksup
REAL      :: temp_real1 (2,tr_ukca), temp_real2 (2,tr_ukca)
REAL      :: mass_tracer (3,tr_ukca)
REAL      :: physics2_mass1 (tr_ukca), physics2_mass2 (tr_ukca)

REAL      :: alfa_za(pdims%k_start:pdims%k_end)
REAL      :: beta_za(pdims%k_start:pdims%k_end)

REAL      :: psi_n (tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,                    &
                       tdims%k_start:tdims%k_end)

REAL      :: psi_np1 (tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,                    &
                       tdims%k_start:tdims%k_end)

REAL      :: tr_min(1:tr_ukca)              ! Minimum allowed value

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_TRACERS_UKCA'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tr_min(:) = 0.0        ! to clip negative values

IF (i_ukca_conserve_method == ukca_no_conserve) THEN

  tracer_ukca = MAX(tracer_ukca, 0.0 )

ELSE

  !------------------------------------------------------------------
  ! Compute the vertical averaging weights
  !------------------------------------------------------------------

  DO k = pdims%k_start, pdims%k_end
    alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /      &
                 ( eta_theta_levels(k) - eta_theta_levels(k-1) )
    beta_za(k) = 1.0 - alfa_za(k)
  END DO

  ! Set values for bottom-most level, assuming tr(:,:,0) = tr(:,:,1)
  alfa_za(pdims%k_start) = 1.0
  beta_za(pdims%k_start) = 0.0
  psi_n  (:,:,tdims%k_start) = 0.0
  psi_np1(:,:,tdims%k_start) = 0.0

  !-------------------------------------------------------------
  ! Compute the 3D array psi(i,j,k) = rho*volume*average
  !   so mass_tracer = SUM[psi*tracer] = SUM[av(tracer)*rho*vol]
  !-------------------------------------------------------------

  DO k = tdims%k_start + 1, tdims%k_end - 1
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                        +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)

        psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                        + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
      END DO
    END DO
  END DO

  k = tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
      psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    END DO
  END DO

  !-------------------------------------------------------------
  ! isolate tracers(n) + phys1 from super_tracer_phys1
  !-------------------------------------------------------------

  tr_temp1(:,:,:,1:tr_ukca) =                                           &
      super_tracer_phys1(:,:,:,eg_tracer_ukca_start:eg_tracer_ukca_end)

  !-------------------------------------------------------------
  ! isolate phys2-increments from super_tracer_phys2
  !-------------------------------------------------------------

  DO kk = 1, tr_ukca
    kksup = kk + eg_tracer_ukca_start - 1
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tracer_ukca_sources(i,j,k,kk) = super_tracer_phys2(i,j,k,kksup)
        END DO
      END DO
    END DO
  END DO

  tracer_ukca_sources(:,:,tdims%k_start,:) =                            &
                        tracer_ukca_sources(:,:,tdims%k_start+1,:)

  tr_temp1(:,:,tdims_l%k_start,:) = tr_temp1(:,:,tdims_l%k_start+1,:)

  tracer_ukca_n(tdims_s%i_start:tdims_s%i_end,                          &
              tdims_s%j_start:tdims_s%j_end,                            &
              tdims_s%k_start:tdims_s%k_end,:)                          &
     = tr_temp1(tdims_s%i_start:tdims_s%i_end,                          &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end,:)

  !------------------------------------------------------------------
  ! Call tri-linear interpolation to compute  tracer_ukca_np1_low
  !------------------------------------------------------------------

  l_high_in=.FALSE.
  l_mono_in=.TRUE.
  CALL eg_interpolation_eta_pmf(                                        &
                     eta_theta_levels,fld_type_w,                       &
                     tr_ukca,                                           &
                     row_length, rows, model_levels+1, rows,            &
                     row_length, rows, model_levels+1,                  &
                     0, 1, l_high_in,l_mono_in,                         &
                     depart_r, depart_lambda, depart_phi,               &
                     mype, nproc, nproc_x, nproc_y,                     &
                     halo_i, halo_j,                                    &
                     global_row_length, datastart, at_extremity,        &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                     offx, offy, error_code,                            &
                     tr_temp1, tracer_ukca_np1_low                    )

  ! Apply separate interpolation method, only if user has requested a 
  ! method other than the one used for Tracers in this run.
  IF ( i_ukca_hiorder_scheme > 0 .AND.                                  &
       i_ukca_hiorder_scheme /= high_order_scheme(tracer_sl) ) THEN

    CALL eg_interpolation_eta_pmf(                                      &
                     eta_theta_levels,fld_type_w,                       &
                     tr_ukca,                                           &
                     row_length, rows, model_levels+1, rows,            &
                     row_length, rows, model_levels+1,                  &
                     i_ukca_hiorder_scheme,                             &
                     monotone_scheme(tracer_sl),                        &
                     l_high(tracer_sl), l_mono(tracer_sl),              &
                     depart_r, depart_lambda, depart_phi,               &
                     mype, nproc, nproc_x, nproc_y,                     &
                     halo_i, halo_j,                                    &
                     global_row_length, datastart, at_extremity,        &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                     offx, offy, error_code,                            &
                     tr_temp1, tracer_ukca                        )

  END IF     ! Use separate interpolation for UKCA

  IF ( .NOT. L_ukca_src_in_conservation ) THEN  

    DO kk = 1, tr_ukca
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            tracer_ukca(i,j,k,kk) = tracer_ukca(i,j,k,kk) -            &
                                        tracer_ukca_sources(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO

    tracer_ukca_sources2 = 0.0

  ELSE
    DO kk = 1, tr_ukca
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            tracer_ukca_np1_low(i,j,k,kk) =                            &
                                      tracer_ukca_np1_low(i,j,k,kk) +  &
                                      tracer_ukca_sources(i,j,k,kk)
            tracer_ukca_sources2(i,j,k,kk) =                           &
                                      tracer_ukca_sources(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO

  END IF

  tracer_ukca_np1_low(:,:,tdims%k_start,:) =                            &
                          tracer_ukca_np1_low(:,:,tdims%k_start+1,:)
  tracer_ukca(:,:,tdims%k_start,:)     =                                &
                          tracer_ukca(:,:,tdims%k_start+1,:)

  !------------------------------------------------------------------
  ! Apply conservation to tracer_ukca
  !------------------------------------------------------------------

  IF ( i_ukca_conserve_method == priestley_old ) THEN

    CALL priestley_algorithm( tracer_ukca_n, tracer_ukca,               &
                              tracer_ukca_np1_low,                      &
                              tracer_ukca_sources2,                     &
                              psi_n, psi_np1,                           &
                              tdims_s%i_start, tdims_s%i_end,           &
                              tdims_s%j_start, tdims_s%j_end,           &
                              tdims_s%k_start, tdims_s%k_end,           &
                              tdims%i_start  , tdims%i_end,             &
                              tdims%j_start  , tdims%j_end,             &
                              tdims%k_start  , tdims%k_end,             &
                              tr_ukca                                 )
  ELSE        ! Priestley optimised

    CALL priestley_algorithm2(tracer_ukca_n, tracer_ukca,              &
                                 tracer_ukca_np1_low,                  &
                                 tracer_ukca_sources2,                 &
                                 psi_n, psi_np1,                       &
                                 tdims_s%i_start, tdims_s%i_end,       &
                                 tdims_s%j_start, tdims_s%j_end,       &
                                 tdims_s%k_start, tdims_s%k_end,       &
                                 tdims%i_start  , tdims%i_end,         &
                                 tdims%j_start  , tdims%j_end,         &
                                 tdims%k_start  , tdims%k_end,         &
                                 tr_ukca                             )

  END IF

  IF ( .NOT. L_ukca_src_in_conservation  ) THEN  

    DO kk = 1, tr_ukca
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            tracer_ukca(i,j,k,kk) = tracer_ukca(i,j,k,kk) +            &
                                      tracer_ukca_sources(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO

  END IF

  !=========================================================================
  ! Impose minimum values on the corrected qs while maintaining conservation
  !=========================================================================

  CALL impose_minima_with_conservation (                                &
                   tracer_ukca, tr_min, psi_np1,                        &
                   tdims_s%i_start, tdims_s%i_end,                      &
                   tdims_s%j_start, tdims_s%j_end,                      &
                   tdims_s%k_start, tdims_s%k_end,                      &
                   tdims%i_start  , tdims%i_end,                        &
                   tdims%j_start  , tdims%j_end,                        &
                   tdims%k_start  , tdims%k_end,                        &
                   tr_ukca                        )

  CALL Swap_Bounds (                                                    &
           tracer_ukca,                                                 &
           tdims%i_len, tdims%j_len,                                    &
           tr_ukca*tdims%k_len,                                         &
           tdims%i_start-tdims_s%i_start,tdims%j_start-tdims_s%j_start, &
           fld_type_p,swap_field_is_scalar)

END IF           ! Conserve / no conserve

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_correct_tracers_ukca

END MODULE eg_correct_tracers_ukca_mod

