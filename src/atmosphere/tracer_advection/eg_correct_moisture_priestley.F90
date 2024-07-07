! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_correct_moisture_priestley_mod

IMPLICIT NONE
! Description:
!            This routine enforces the mass conservation for
!            6 moisture variables using a modified versions of
!            the Priestley scheme. It should work with
!            either mixing ratios and dry density or with
!            specific humdities and wet density.
!            If L_conserve_mass=.FALSE., then the routine make
!            sure that the moisture variables are above the
!            minimum values required only (doesn't impose
!            any mass conservation constraint).
!
! Method: ENDGame formulation version 4.00(Feb 2014)

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: tracer_advection
!
! Code description:
! Language: Fortran 90.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='EG_CORRECT_MOISTURE_PRIESTLEY_MOD'

CONTAINS

SUBROUTINE eg_correct_moisture_priestley(rho_n, rho_np1,               &
               q1_n, q2_n, q3_n, q4_n, q5_n, q6_n,                     &
               q1_np1, q2_np1, q3_np1, q4_np1, q5_np1, q6_np1,         &
               q1_s, q2_s, q3_s, q4_s, q5_s, q6_s,                     &
               q1min, L_q4, L_q5, L_q6, L_conserve_mass,               &
               g_i_pe, depart_lambda, depart_phi, depart_r  )


USE um_parvars,     ONLY: nproc_x, nproc_y,                            &
                          at_extremity,gc_proc_row_group,              &
                          gc_proc_col_group, datastart
USE um_parcore,     ONLY: mype, nproc

USE nlsizes_namelist_mod,  ONLY: global_row_length

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels
USE atm_fields_bounds_mod, ONLY : tdims, tdims_s, tdims_l,             &
                                  pdims, pdims_s, wdims
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE Field_Types
USE priestley_algorithm_mod
USE sl_input_mod, ONLY: moist_priestley_opt, L_moist_src_in_conserve
USE eg_helmholtz_mod,      ONLY: ec_vol
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE timestep_mod, ONLY: timestep_number, timestep
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER, INTENT(IN) :: g_i_pe(1-(tdims%i_start-tdims_l%i_start):       &
                   global_row_length+(tdims%i_start-tdims_l%i_start) )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_MOISTURE_PRIESTLEY'

REAL, INTENT(IN) :: rho_n(pdims_s%i_start:pdims_s%i_end,                &
                          pdims_s%j_start:pdims_s%j_end,                &
                          pdims_s%k_start:pdims_s%k_end),               &

                    rho_np1(pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN) :: q1_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q2_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q3_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q4_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q5_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q6_n(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q1_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q2_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q3_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q4_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q5_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end),                    &

                    q6_s(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         tdims%k_start:tdims%k_end)

REAL, INTENT(INOUT):: q1_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end),               &

                      q2_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end),               &

                      q3_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end),               &

                      q4_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end),               &

                      q5_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end),               &

                      q6_np1(tdims_s%i_start:tdims_s%i_end,             &
                          tdims_s%j_start:tdims_s%j_end,                &
                          tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(IN)     :: q1min
LOGICAL, INTENT(IN) ::  L_q4, L_q5, L_q6, L_conserve_mass

REAL, INTENT(IN)    :: depart_lambda( wdims%i_start:wdims%i_end,        &
                                      wdims%j_start:wdims%j_end,        &
                                      wdims%k_start:wdims%k_end),       &

                       depart_phi( wdims%i_start:wdims%i_end,           &
                                   wdims%j_start:wdims%j_end,           &
                                   wdims%k_start:wdims%k_end),          &

                       depart_r( wdims%i_start:wdims%i_end,             &
                                 wdims%j_start:wdims%j_end,             &
                                 wdims%k_start:wdims%k_end)

! locals

LOGICAL :: l_high_in   ! local setting for use in CALLs
LOGICAL :: l_mono_in   ! local setting for use in CALLs

INTEGER :: i_halo_big, j_halo_big, i_halo_small, j_halo_small

REAL, ALLOCATABLE :: qs_n(:,:,:,:)
REAL, ALLOCATABLE :: qs_np1(:,:,:,:)
REAL, ALLOCATABLE :: qs_s (:,:,:,:)
REAL, ALLOCATABLE :: qs_np1_low(:,:,:,:)
REAL, ALLOCATABLE :: qs_s2(:,:,:,:)
REAL, ALLOCATABLE :: qs_n2(:,:,:,:)

REAL :: alfa_za(pdims%k_start:pdims%k_end)
REAL :: beta_za(pdims%k_start:pdims%k_end)

REAL :: psi_n (tdims%i_start:tdims%i_end,                              &
               tdims%j_start:tdims%j_end,                              &
               tdims%k_start:tdims%k_end)
REAL :: psi_np1(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                              &
               tdims%k_start:tdims%k_end)

INTEGER :: i,j,k,kk,error_code
REAL, ALLOCATABLE :: qsmin(:)

INTEGER :: num_qs, pt_q4, pt_q5, pt_q6


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


i_halo_big   = tdims%i_start-tdims_l%i_start
i_halo_small = tdims%i_start-tdims_s%i_start
j_halo_big   = tdims%j_start-tdims_l%j_start
j_halo_small = tdims%j_start-tdims_s%j_start

num_qs = 3

IF ( L_q4  ) THEN
  num_qs = num_qs + 1
  pt_q4  = num_qs
END IF

IF ( L_q5 ) THEN
  num_qs = num_qs + 1
  pt_q5  = num_qs
END IF

IF ( L_q6   ) THEN
  num_qs = num_qs + 1
  pt_q6  = num_qs
END IF

ALLOCATE (qs_n (tdims_s%i_start:tdims_s%i_end,                        &
                tdims_s%j_start:tdims_s%j_end,                        &
                tdims_s%k_start:tdims_s%k_end,                        &
                num_qs))

ALLOCATE (qs_np1(tdims_s%i_start:tdims_s%i_end,                        &
                 tdims_s%j_start:tdims_s%j_end,                        &
                 tdims_s%k_start:tdims_s%k_end,                        &
                 num_qs))

ALLOCATE (qs_s(tdims_s%i_start:tdims_s%i_end,                         &
               tdims_s%j_start:tdims_s%j_end,                         &
               tdims_s%k_start:tdims_s%k_end,                         &
               num_qs))

ALLOCATE (qs_np1_low(tdims_s%i_start:tdims_s%i_end,                    &
                     tdims_s%j_start:tdims_s%j_end,                    &
                     tdims_s%k_start:tdims_s%k_end,                    &
                     num_qs))

ALLOCATE (qs_s2(tdims_s%i_start:tdims_s%i_end,                         &
                tdims_s%j_start:tdims_s%j_end,                         &
                tdims_s%k_start:tdims_s%k_end,                         &
                num_qs))

ALLOCATE  (qs_n2 (tdims_l%i_start:tdims_l%i_end,                       &
                  tdims_l%j_start:tdims_l%j_end,                       &
                  tdims_l%k_start:tdims_l%k_end,                       &
                  num_qs))

ALLOCATE (qsmin(num_qs))
qsmin(1)               = q1min
qsmin(2:num_qs) = 0.0

IF ( L_conserve_mass ) THEN

  !========================================================
  ! copy the 6 moisture variables  and sources into a
  ! 4-dimensional superarray
  !========================================================

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(k, j, i)                                                 &
!$OMP SHARED(tdims, qs_n, qs_s, qs_np1, q1_n, q1_s, q1_np1, q2_n,      & 
!$OMP        q2_s, q2_np1, q3_n, q3_s, q3_np1, L_q4, L_q5, L_q6,       &
!$OMP        pt_q4, pt_q5, pt_q6, q4_n, q4_s, q4_np1, q5_n, q5_s,      &
!$OMP        q5_np1, q6_n, q6_s, q6_np1)

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,1) = q1_n(i,j,k)
        qs_s(i,j,k,1) = q1_s(i,j,k)
        qs_np1(i,j,k,1) = q1_np1(i,j,k)
      END DO
    END DO
  END DO 
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,2) = q2_n(i,j,k)
        qs_s(i,j,k,2) = q2_s(i,j,k)
        qs_np1(i,j,k,2) = q2_np1(i,j,k)
      END DO
    END DO
  END DO 
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,3) = q3_n(i,j,k)
        qs_s(i,j,k,3) = q3_s(i,j,k)
        qs_np1(i,j,k,3) = q3_np1(i,j,k)
      END DO
    END DO
  END DO 
!$OMP END DO NOWAIT

  IF( L_q4 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_n(i,j,k,pt_q4)   = q4_n(i,j,k)
          qs_s(i,j,k,pt_q4)   = q4_s(i,j,k)
          qs_np1(i,j,k,pt_q4) = q4_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q5 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_n(i,j,k,pt_q5)   = q5_n(i,j,k)
          qs_s(i,j,k,pt_q5)   = q5_s(i,j,k)
          qs_np1(i,j,k,pt_q5) = q5_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q6 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_n(i,j,k,pt_q6)   = q6_n(i,j,k)
          qs_s(i,j,k,pt_q6)   = q6_s(i,j,k)
          qs_np1(i,j,k,pt_q6) = q6_np1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
!$OMP END PARALLEL

  !===========================================================================
  ! make sure the boundary condition field(surface)=field(level 1) is imposed.
  ! This may be not necessary here because it may already be done higher up
  ! but for safety we keep it.
  !=========================================================================

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(tdims, qs_n, qs_s, qs_np1, qs_n2, num_qs)

  k = tdims%k_start
!$OMP DO SCHEDULE(STATIC)
  DO kk = 1, num_qs
    DO j = tdims%j_start, tdims%j_end          
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,kk)   = qs_n(i,j,k+1,kk)
        qs_s(i,j,k,kk)   = qs_s(i,j,k+1,kk)
        qs_np1(i,j,k,kk) = qs_np1(i,j,k+1,kk)
      END DO
    END DO  
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO kk = 1, num_qs
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_n2(i,j,k,kk) = qs_n(i,j,k,kk)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  CALL Swap_Bounds(qs_n2,                                              &
          tdims%i_len,tdims%j_len, num_qs*tdims_l%k_len,               &
          tdims_l%halo_i,tdims_l%halo_j, fld_type_p, swap_field_is_scalar)

  !=========================================================================
  ! Compute the vertical linear averaging weights
  !=========================================================================

  DO k = pdims%k_start, pdims%k_end
    alfa_za(k) = intw_w2rho(k,1)
    beta_za(k) = intw_w2rho(k,2)
  END DO

  ! Set values for bottom-most level, assuming q(:,:,0) = q(:,:,1)
  alfa_za(pdims%k_start) = 1.0
  beta_za(pdims%k_start) = 0.0

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(k, j, i)                                                 &
!$OMP SHARED(tdims, psi_n, psi_np1, rho_n, ec_vol, alfa_za, beta_za,   &
!$OMP        rho_np1)

!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n  (i,j,tdims%k_start) = 0.0
      psi_np1(i,j,tdims%k_start) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

  !=========================================================================
  ! Compute the 3D array psi(i,j,k) = rho*volume*average
  !   so mass_tracer = SUM[psi*tracer] = SUM[av(tracer)*rho*vol]
  !=========================================================================

  DO k = tdims%k_start + 1, tdims%k_end - 1
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        psi_n(i,j,k)  =  rho_n(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )    &
                      +  rho_n(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)

        psi_np1(i,j,k) = rho_np1(i,j,k  )*ec_vol(i,j,k  )*alfa_za(k  )  &
                       + rho_np1(i,j,k+1)*ec_vol(i,j,k+1)*beta_za(k+1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  k = tdims%k_end
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
      psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  !=====================================================================
  ! Call tri-linear interpolation to compute qs_np1_low
  !=====================================================================

  l_high_in=.FALSE.
  l_mono_in=.TRUE.
  CALL eg_interpolation_eta_pmf(                                       &
                    eta_theta_levels,fld_type_w, num_qs,               &
                    tdims%i_len, tdims%j_len, tdims%k_len,             &
                    tdims%j_len, tdims%i_len, tdims%j_len,             &
                    tdims%k_len,0, 1,                                  &
                    l_high_in,l_mono_in,                               &
                    depart_r, depart_lambda, depart_phi,               &
                    mype, nproc, nproc_x, nproc_y,                     &
                    i_halo_big, j_halo_big,                            &
                    global_row_length, datastart, at_extremity,        &
                    g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                    i_halo_small, j_halo_small,                        &
                    error_code,                                        &
                    qs_n2, qs_np1_low                              )

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(L_moist_src_in_conserve, tdims, num_qs, qs_np1, qs_s,     &
!$OMP        qs_s2, qs_np1_low)

  IF ( .NOT. L_moist_src_in_conserve ) THEN  
    DO kk = 1, num_qs
      DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) - qs_s(i,j,k,kk)
            qs_s2(i,j,k,kk) = 0.0
          END DO
        END DO
!$OMP END DO NOWAIT
      END DO
    END DO
  ELSE
    DO kk = 1, num_qs
      DO k = tdims%k_start, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1_low(i,j,k,kk) = qs_np1_low(i,j,k,kk) + qs_s(i,j,k,kk)
            qs_s2(i,j,k,kk) =  qs_s(i,j,k,kk)
          END DO
        END DO
!$OMP END DO NOWAIT
      END DO
    END DO

  END IF

  k = tdims%k_start

  DO kk = 1, num_qs
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_np1_low(i,j,k,kk) = qs_np1_low(i,j,k+1,kk)
        qs_np1(i,j,k,kk)     = qs_np1(i,j,k+1,kk)
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

!$OMP END PARALLEL

  !=====================================================================
  ! Apply conservation to qs_np1
  !=====================================================================

  IF ( moist_priestley_opt == 1 ) THEN    ! Original Priestley method

    CALL priestley_algorithm( qs_n, qs_np1, qs_np1_low, qs_s2,         &
                              psi_n, psi_np1,                          &
                              tdims_s%i_start, tdims_s%i_end,          &
                              tdims_s%j_start, tdims_s%j_end,          &
                              tdims_s%k_start, tdims_s%k_end,          &
                              tdims%i_start  , tdims%i_end,            &
                              tdims%j_start  , tdims%j_end,            &
                              tdims%k_start  , tdims%k_end,            &
                              num_qs                             )
  ELSE                                    ! Optimised Priestley

    CALL priestley_algorithm2(qs_n, qs_np1, qs_np1_low, qs_s2,         &
                              psi_n, psi_np1,                          &
                              tdims_s%i_start, tdims_s%i_end,          &
                              tdims_s%j_start, tdims_s%j_end,          &
                              tdims_s%k_start, tdims_s%k_end,          &
                              tdims%i_start  , tdims%i_end,            &
                              tdims%j_start  , tdims%j_end,            &
                              tdims%k_start  , tdims%k_end,            &
                              num_qs                              )

  END IF



  IF ( .NOT. L_moist_src_in_conserve ) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP PRIVATE(kk, k, j, i)                                             &
!$OMP SHARED(L_moist_src_in_conserve, tdims, num_qs, qs_np1, qs_s)
    DO kk = 1,num_qs
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qs_np1(i,j,k,kk) = qs_np1(i,j,k,kk) + qs_s(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF


  !=========================================================================
  ! Impose minimum values on the corrected qs while maintaining conservation
  !=========================================================================

  CALL impose_minima_with_conservation(qs_np1, qsmin, psi_np1,         &
                                      tdims_s%i_start, tdims_s%i_end,  &
                                      tdims_s%j_start, tdims_s%j_end,  &
                                      tdims_s%k_start, tdims_s%k_end,  &
                                      tdims%i_start  , tdims%i_end,    &
                                      tdims%j_start  , tdims%j_end,    &
                                      tdims%k_start  , tdims%k_end,    &
                                      num_qs                        )

  !=====================================================================
  ! put the corrected moisture variables back into their original arrays
  !=====================================================================

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(k, j, i)                                                 &
!$OMP SHARED(tdims, q1_np1, q2_np1, q3_np1, qs_np1, L_q4,              &
!$OMP        L_q5, L_q6, q4_np1, q5_np1, q6_np1, pt_q4, pt_q5, pt_q6)

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q1_np1(i,j,k) = qs_np1(i,j,k,1)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q2_np1(i,j,k) = qs_np1(i,j,k,2)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q3_np1(i,j,k) = qs_np1(i,j,k,3)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF( L_q4 ) THEN 
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end    
          q4_np1(i,j,k) = qs_np1(i,j,k,pt_q4)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q5 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end   
          q5_np1(i,j,k) = qs_np1(i,j,k,pt_q5)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q6 ) THEN 
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end   
          q6_np1(i,j,k) = qs_np1(i,j,k,pt_q6)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
!$OMP END PARALLEL 

ELSE  ! if L_conserve_mass==.false. then simply make sure that the moisture
     ! variables are above the minimum values required

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP PRIVATE(k, j, i)                                                 &
!$OMP SHARED(q1_np1, q1min, q2_np1, q3_np1, q4_np1, q5_np1, q6_np1,    &
!$OMP        L_q4, L_q5, L_q6, tdims)

!$OMP DO SCHEDULE(STATIC)
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q1_np1(i,j,k) = MAX(q1_np1(i,j,k), q1min)
        q2_np1(i,j,k) = MAX(q2_np1(i,j,k), 0.0  )
        q3_np1(i,j,k) = MAX(q3_np1(i,j,k), 0.0  )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF( L_q4 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end  
          q4_np1(i,j,k) = MAX(q4_np1(i,j,k), 0.0  )
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q5 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end  
          q5_np1(i,j,k) = MAX(q5_np1(i,j,k), 0.0  )
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF( L_q6 ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end 
          q6_np1(i,j,k) = MAX(q6_np1(i,j,k), 0.0  )
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
!$OMP END PARALLEL

END IF  ! L_conserve_mass

DEALLOCATE (qs_n)
DEALLOCATE (qs_np1)
DEALLOCATE (qs_s)
DEALLOCATE (qs_np1_low)
DEALLOCATE (qs_s2)
DEALLOCATE (qs_n2)
DEALLOCATE (qsmin)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_correct_moisture_priestley

END MODULE eg_correct_moisture_priestley_mod
