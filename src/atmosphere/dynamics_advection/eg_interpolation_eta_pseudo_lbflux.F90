! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_interpolation_eta_pseudo_lbflux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='EG_INTERPOLATION_ETA_PSEUDO_LBFLUX_MOD'

CONTAINS
SUBROUTINE eg_interpolation_eta_pseudo_lbflux(                    &
                          eta_in,  pnt_type,                      &
                          number_of_inputs,                       &
                          dim_i_in, dim_j_in, dim_k_in,           &
                          dim_j_in_w,                             &
                          dim_i_out, dim_j_out, dim_k_out,        &
                          high_order_scheme, monotone_scheme,     &
                          l_high, l_mono,                         &
                          eta_out, lambda_out, phi_out,           &
                          me, n_proc, n_procx, n_procy,           &
                          halo_i, halo_j, g_row_length,           &
                          datastart, at_extremity, g_i_pe,        &
                          proc_row_group, proc_col_group,         &
                          halo_data_out_i, halo_data_out_j,       &
                          off_x, off_y,                           &
                          error_code,                             &
                          data_in, data_out,                      &
                          rho_n, rho_np1, l_conserv,              &
                          pseudo_lbflux,                          &
                          k_int_linear)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE umPrintMgr,  ONLY: printstatus, prstatus_diag, &
                        ummessage, umprint

USE interp_grid_const_mod

USE lbc_mod,                      ONLY: rimwidtha
USE rimtypes,                     ONLY: rima_type_norm
USE atm_fields_bounds_mod,        ONLY: pdims, pdims_s, tdims, tdims_s
USE level_heights_mod,            ONLY: eta_theta_levels, eta_rho_levels
USE eg_helmholtz_mod,             ONLY: ec_vol
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE eg_mass_conserv_mod,          ONLY: eg_mass_conservation_fix
USE global_2d_sums_mod,           ONLY: global_2d_sums
USE eg_lam_domain_kind_mod,       ONLY: lam_domain_kind
USE eg_total_mass_region_mod,     ONLY: eg_total_mass_region
USE model_domain_mod,             ONLY: model_type, mt_lam

IMPLICIT NONE
! Description:
!             This routine is the wrapper of eg_interpolation_eta to
!             compute the pseudo lateral boundary flux (PLF).
!
! Method: ENDGame formulation version 4.xx
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INTERPOLATION_ETA_PSEUDO_LBFLUX'


! Arguments with Intent IN. ie: Input variables.

INTEGER ::                                                            &
  dim_i_in                                                            &
              ! Dimension of Data_in in i direction.
, dim_j_in                                                            &
              ! Dimension of Data_in in j direction.
, dim_j_in_w                                                          &
              ! Dimension of Data_in in j direction.
, dim_k_in                                                            &
              ! Dimension of Data_in in k direction.
, dim_i_out                                                           &
              ! Dimension of Data_out in i direction.
, dim_j_out                                                           &
              ! Dimension of Data_out in j direction.
, dim_k_out                                                           &
              ! Dimension of Data_out in k direction.
, me                                                                  &
              ! My processor number
, n_proc                                                              &
              ! Total number of processors
, n_procx                                                             &
              ! Number of processors in longitude
, n_procy                                                             &
              ! Number of processors in latitude
, halo_i                                                              &
              ! Size of halo in i direction.
, halo_j                                                              &
              ! Size of halo in j direction.
, off_x                                                               &
, off_y                                                               &
, halo_data_out_i                                                     &
                  ! size of data out halo in i direction
, halo_data_out_j                                                     &
                  ! size of data out halo in j direction
, proc_row_group                                                      &
                 ! Group id for processors on the same row
, proc_col_group                                                      &
                 ! Group id for processors on the same column
, g_row_length                                                        &
               ! global number of points on a row
, datastart(3)                                                        &
               ! First gridpoints held by this processor.
, g_i_pe(1-halo_i:g_row_length+halo_i) ! processor on my procr-row
                       ! holding a given value in i direction

INTEGER ::                                                            &
  pnt_type                                                            &
              ! Defines via an integer code the nature of the
              ! interpolation in terms of which grid the input
              ! Data is on. The codes are given in
              ! terms of a primary variable that would be held
              ! at that point and are u=1, v=2, w=3.
, high_order_scheme                                                   &
                     ! a code saying which high order scheme to
                     ! use.
, monotone_scheme                                                     &
                  ! a code saying which monotone scheme to use.
, number_of_inputs
                   !number of fields to interpolate.

LOGICAL ::                                                            &
  l_high                                                              &
                 ! True, if high order interpolation required.
, l_mono
                 ! True, if interpolation required to be monotone.

REAL ::                                                               &
  eta_in(dim_k_in) ! eta coordinate levels.

REAL ::                                                               &
  data_in  (1-halo_i:dim_i_in+halo_i,                                 &
            1-halo_j:dim_j_in+halo_j, 0:dim_k_in-1,                   &
            number_of_inputs )
                                                ! data to be
                                                ! interpolated

REAL ::                                                               &
  lambda_out (dim_i_out, dim_j_out, 0:dim_k_out-1)                    &
                                                ! Lambda
                                                ! co-ordinate of
                                                ! output data on
                                                ! input.
, phi_out (dim_i_out, dim_j_out, 0:dim_k_out-1)                       &
                                                ! Phi Co-ordinate
                                                ! of output data
                                                ! on input.
, eta_out (dim_i_out, dim_j_out, 0:dim_k_out-1) ! Vertical
                                                ! co-ordinate
                                                ! of output data.

LOGICAL ::                                                            &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid


REAL,    INTENT(IN)    ::  rho_n(pdims_s%i_start:pdims_s%i_end,       &
                                 pdims_s%j_start:pdims_s%j_end,       &
                                 pdims_s%k_start:pdims_s%k_end)

REAL,    INTENT(IN)    :: rho_np1(pdims_s%i_start:pdims_s%i_end,      &
                                  pdims_s%j_start:pdims_s%j_end,      &
                                  pdims_s%k_start:pdims_s%k_end)

LOGICAL, INTENT(IN)    :: l_conserv

INTEGER, INTENT(IN), OPTIONAL :: k_int_linear

! Arguments with Intent OUT. ie: Output variables.
REAL ::                                                               &
     ! data interpolated to desired locations.
  data_out (1-halo_data_out_i:dim_i_out+halo_data_out_i,              &
            1-halo_data_out_j:dim_j_out+halo_data_out_j,              &
            0:dim_k_out-1, number_of_inputs)

REAL, INTENT(OUT) :: pseudo_lbflux(number_of_inputs)

INTEGER ::                                                            &
  error_code     ! Non-zero on exit if error detected.

!LOCAL VARIABLES
LOGICAL :: l_mono_in ! local setting for use in CALLs

REAL,    ALLOCATABLE :: data_in3(:,:,:,:,:)
REAL,    ALLOCATABLE :: data_out3(:,:,:,:,:)
REAL,    ALLOCATABLE :: qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_s(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE :: psi_n(:,:,:)
REAL,    ALLOCATABLE :: psi_np1(:,:,:)
REAL,    ALLOCATABLE :: qsmin(:)
REAL,    ALLOCATABLE :: local_sums(:,:,:)
REAL,    ALLOCATABLE :: alfa_za(:)
REAL,    ALLOCATABLE :: beta_za(:)

INTEGER :: i, j, k, kk
INTEGER :: nrim
INTEGER :: IS, ie, js, je
REAL    :: mass_n(number_of_inputs)
REAL    :: mass_sl(number_of_inputs)
REAL    :: mass_slmono(number_of_inputs)
REAL    :: errsl(number_of_inputs)
REAL    :: errslmono(number_of_inputs)

LOGICAL :: L_check_sl_mass = .TRUE.
LOGICAL :: L_conserv_smooth_lap
! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (model_type == mt_lam .AND. l_conserv) THEN
  ALLOCATE(data_in3(1-halo_i:dim_i_in+halo_i,                     &
                    1-halo_j:dim_j_in+halo_j, 0:dim_k_in-1,       &
                    number_of_inputs, 2)                          )

  ALLOCATE(data_out3(1-halo_data_out_i:dim_i_out+halo_data_out_i, &
              1-halo_data_out_j:dim_j_out+halo_data_out_j,        &
              0:dim_k_out-1, number_of_inputs, 2)                 )

  ALLOCATE(    qs_n(tdims_s%i_start:tdims_s%i_end,  &
                    tdims_s%j_start:tdims_s%j_end,  &
                    tdims_s%k_start:tdims_s%k_end,  &
                    number_of_inputs)               )

  ALLOCATE(    qs_s(tdims_s%i_start:tdims_s%i_end,  &
                    tdims_s%j_start:tdims_s%j_end,  &
                    tdims_s%k_start:tdims_s%k_end,  &
                    number_of_inputs)               )

  ALLOCATE(  qs_np1(tdims_s%i_start:tdims_s%i_end,  &
                    tdims_s%j_start:tdims_s%j_end,  &
                    tdims_s%k_start:tdims_s%k_end,  &
                    number_of_inputs)               )

  ALLOCATE(   psi_n(tdims%i_start:tdims%i_end,  &
                    tdims%j_start:tdims%j_end,  &
                    tdims%k_start:tdims%k_end)  )

  ALLOCATE( psi_np1(tdims%i_start:tdims%i_end,  &
                    tdims%j_start:tdims%j_end,  &
                    tdims%k_start:tdims%k_end)  )

  ALLOCATE(  qsmin(number_of_inputs)     )
  ALLOCATE(  local_sums(tdims%i_start:tdims%i_end,  &
                        tdims%j_start:tdims%j_end,  &
                        number_of_inputs)           )

  ALLOCATE(alfa_za(pdims%k_start:pdims%k_end))
  ALLOCATE(beta_za(pdims%k_start:pdims%k_end))


  CALL eg_total_mass_region(IS, ie, js, je, l_exclude_rim=.TRUE.)

!$OMP PARALLEL DO DEFAULT(NONE)                               &
!$OMP& PRIVATE(k, kk, j, i)                                   &
!$OMP& SHARED(dim_k_in, number_of_inputs, halo_j, dim_j_in,   &
!$OMP&        halo_i, dim_i_in,                               &
!$OMP&        lam_domain_kind, data_in3, data_in)
  DO k = 0, dim_k_in-1
    DO kk = 1, number_of_inputs
      DO j = 1-halo_j, dim_j_in+halo_j
        DO i = 1-halo_i, dim_i_in+halo_i
          IF (lam_domain_kind(i,j) == 1) THEN
            ! forecast domain close to the boundaries
            data_in3(i,j,k,kk,1) = data_in(i,j,k,kk)
            data_in3(i,j,k,kk,2) = 0.0
          ELSE IF (lam_domain_kind(i,j) == 0) THEN
            ! rim
            data_in3(i,j,k,kk,1) = 0.0
            data_in3(i,j,k,kk,2) = data_in(i,j,k,kk)
          ELSE IF (lam_domain_kind(i,j) == 2) THEN
            ! forecast domain far enough from the boundaries
            data_in3(i,j,k,kk,1) = 0.0
            data_in3(i,j,k,kk,2) = 0.0
          END IF
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  l_mono_in=.FALSE.

  IF (PRESENT(k_int_linear)) THEN
    CALL eg_interpolation_eta_pmf(                                    &
                          eta_in,  pnt_type,                          &
                          2*number_of_inputs,                         &
                          dim_i_in, dim_j_in, dim_k_in,               &
                          dim_j_in_w,                                 &
                          dim_i_out, dim_j_out, dim_k_out,            &
                          high_order_scheme, monotone_scheme,         &
                          l_high, l_mono_in,                          &
                          eta_out, lambda_out, phi_out,               &
                          me, n_proc, n_procx, n_procy,               &
                          halo_i, halo_j, g_row_length,               &
                          datastart, at_extremity, g_i_pe,            &
                          proc_row_group, proc_col_group,             &
                          halo_data_out_i, halo_data_out_j,           &
                          error_code,                                 &
                          data_in3, data_out3,                        &
                          k_int_linear_in = k_int_linear)
  ELSE
    CALL eg_interpolation_eta_pmf(                                    &
                          eta_in,  pnt_type,                          &
                          2*number_of_inputs,                         &
                          dim_i_in, dim_j_in, dim_k_in,               &
                          dim_j_in_w,                                 &
                          dim_i_out, dim_j_out, dim_k_out,            &
                          high_order_scheme, monotone_scheme,         &
                          l_high, l_mono_in,                          &
                          eta_out, lambda_out, phi_out,               &
                          me, n_proc, n_procx, n_procy,               &
                          halo_i, halo_j, g_row_length,               &
                          datastart, at_extremity, g_i_pe,            &
                          proc_row_group, proc_col_group,             &
                          halo_data_out_i, halo_data_out_j,           &
                          error_code,                                 &
                          data_in3, data_out3)
  END IF

!$OMP PARALLEL DO DEFAULT(NONE)                      &
!$OMP& PRIVATE(k, kk, j, i)                          &
!$OMP& SHARED(dim_k_in, number_of_inputs, tdims,     &
!$OMP&        qs_n, data_in3, qs_np1, data_out3, qs_s)
  DO k = 0, dim_k_in-1
    DO kk = 1, number_of_inputs
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qs_n(i,j,k,kk)   =  data_in3(i,j,k,kk,1)
          qs_np1(i,j,k,kk) = data_out3(i,j,k,kk,1)
          qs_s(i,j,k,kk)   = 0.0 ! qs_s is an array for AtmosPhysics2,
                                ! therefore not necessary at this point.
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

    !conservation for the region including rim
  qsmin(:) = -1.0e+30  ! i.e. no qsmin
  CALL eg_mass_conservation_fix(rho_n, rho_np1, qs_n, qs_np1, qs_s,   &
                                qsmin,            number_of_inputs,   &
                                L_conserv_smooth_lap    )


  DO k = pdims%k_start+1, pdims%k_end
    alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) ) /    &
               ( eta_theta_levels(k) - eta_theta_levels(k-1) )
    beta_za(k) = 1.0 - alfa_za(k)
  END DO
  alfa_za(pdims%k_start) = 1.0
  beta_za(pdims%k_start) = 0.0

  k = tdims%k_start

!$OMP PARALLEL DO DEFAULT(NONE)              &
!$OMP& PRIVATE(j, i)                         &
!$OMP& SHARED(k, tdims,                      &
!$OMP&        psi_n, rho_n, ec_vol, beta_za, &
!$OMP&        psi_np1, rho_np1)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k) =   rho_n(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)
      psi_np1(i,j,k) = rho_np1(i,j,k+1) * ec_vol(i,j,k+1) * beta_za(k+1)
    END DO
  END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE)                       &
!$OMP& PRIVATE(k, j, i)                               &
!$OMP& SHARED(tdims,                                  &
!$OMP&        psi_n, rho_n, ec_vol, alfa_za, beta_za, &
!$OMP&        psi_np1,rho_np1)
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
!$OMP END PARALLEL DO

  k = tdims%k_end
!$OMP PARALLEL DO DEFAULT(NONE)              &
!$OMP& PRIVATE(j, i)                         &
!$OMP& SHARED(k, tdims,                      &
!$OMP&        psi_n, rho_n, ec_vol, alfa_za, &
!$OMP&        psi_np1, rho_np1)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      psi_n(i,j,k)   =   rho_n(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
      psi_np1(i,j,k) = rho_np1(i,j,k) * ec_vol(i,j,k) * alfa_za(k)
    END DO
  END DO
!$OMP END PARALLEL DO

    ! calculate pseudo lateral boundary flux such that
    ! the total mass in the forecast region is conserved.
!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(kk, k, j, i)                               &
!$OMP& SHARED(number_of_inputs, dim_k_in, js, je, is, ie, &
!$OMP&        local_sums, data_out3, psi_np1,             &
!$OMP&        qs_np1, data_in3, psi_n)
  DO kk = 1, number_of_inputs
    local_sums(:,:,kk) = 0.0
    DO k = 0, dim_k_in-1
      DO j = js, je
        DO i = IS, ie
          local_sums(i,j,kk) = local_sums(i,j,kk)    &
            + data_out3(i,j,k,kk,2) * psi_np1(i,j,k) & !rim non-conserv @ t=n+1
            + qs_np1(i,j,k,kk)      * psi_np1(i,j,k) & !fcst conserv    @ t=n+1
            - data_in3(i,j,k,kk,1)  * psi_n(i,j,k)     !inside rim      @ t=n
        END DO
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL global_2d_sums(local_sums, tdims%i_len,           &
                                  tdims%j_len,           &
                                  0, 0, number_of_inputs, pseudo_lbflux)

  IF (printstatus >= prstatus_diag) THEN
!$OMP PARALLEL DO DEFAULT(NONE)                   &
!$OMP& PRIVATE(kk, k, j, i)                       &
!$OMP& SHARED(number_of_inputs, dim_k_in,         &
!$OMP&        halo_j, dim_j_in, halo_i, dim_i_in, &
!$OMP&        data_in3, data_in)
    DO k = 0, dim_k_in-1
      DO kk = 1, number_of_inputs
        DO j = 1-halo_j, dim_j_in+halo_j
          DO i = 1-halo_i, dim_i_in+halo_i
            data_in3(i,j,k,kk,1) = data_in(i,j,k,kk)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    l_mono_in=.FALSE.
    CALL eg_interpolation_eta_pmf(                                    &
                          eta_in,  pnt_type,                          &
                          number_of_inputs,                           &
                          dim_i_in, dim_j_in, dim_k_in,               &
                          dim_j_in_w,                                 &
                          dim_i_out, dim_j_out, dim_k_out,            &
                          high_order_scheme, monotone_scheme,         &
                          l_high, l_mono_in,                          &
                          eta_out, lambda_out, phi_out,               &
                          me, n_proc, n_procx, n_procy,               &
                          halo_i, halo_j, g_row_length,               &
                          datastart, at_extremity, g_i_pe,            &
                          proc_row_group, proc_col_group,             &
                          halo_data_out_i, halo_data_out_j,           &
                          error_code,                                 &
                          data_in3, data_out3)
!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(kk, k, j, i)                               &
!$OMP& SHARED(number_of_inputs, dim_k_in, js, je, is, ie, &
!$OMP&        local_sums, data_out3, psi_np1)
    DO kk = 1, number_of_inputs
      local_sums(:,:,kk) = 0.0
      DO k = 0, dim_k_in-1
        DO j = js, je
          DO i = IS, ie
            local_sums(i,j,kk) = local_sums(i,j,kk)    &
              + data_out3(i,j,k,kk,1) * psi_np1(i,j,k)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    CALL global_2d_sums(local_sums, tdims%i_len,   &
                                    tdims%j_len,   &
                                    0, 0, number_of_inputs, mass_sl)
    l_mono_in=.TRUE.
    CALL eg_interpolation_eta_pmf(                                    &
                          eta_in,  pnt_type,                          &
                          number_of_inputs,                           &
                          dim_i_in, dim_j_in, dim_k_in,               &
                          dim_j_in_w,                                 &
                          dim_i_out, dim_j_out, dim_k_out,            &
                          high_order_scheme, monotone_scheme,         &
                          l_high, l_mono_in,                          &
                          eta_out, lambda_out, phi_out,               &
                          me, n_proc, n_procx, n_procy,               &
                          halo_i, halo_j, g_row_length,               &
                          datastart, at_extremity, g_i_pe,            &
                          proc_row_group, proc_col_group,             &
                          halo_data_out_i, halo_data_out_j,           &
                          error_code,                                 &
                          data_in3, data_out3)
!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(kk, k, j, i)                               &
!$OMP& SHARED(number_of_inputs, dim_k_in, js, je, is, ie, &
!$OMP&        local_sums, data_out3, psi_np1)
    DO kk = 1, number_of_inputs
      local_sums(:,:,kk) = 0.0
      DO k = 0, dim_k_in-1
        DO j = js, je
          DO i = IS, ie
            local_sums(i,j,kk) = local_sums(i,j,kk)    &
              + data_out3(i,j,k,kk,1) * psi_np1(i,j,k)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    CALL global_2d_sums(local_sums, tdims%i_len,       &
                                    tdims%j_len,       &
                                    0, 0, number_of_inputs, mass_slmono)

!$OMP PARALLEL DO DEFAULT(NONE)                           &
!$OMP& PRIVATE(kk, k, j, i)                               &
!$OMP& SHARED(number_of_inputs, dim_k_in, js, je, is, ie, &
!$OMP&        local_sums, data_in3, psi_n)
    DO kk = 1, number_of_inputs
      local_sums(:,:,kk) = 0.0
      DO k = 0, dim_k_in-1
        DO j = js, je
          DO i = IS, ie
            local_sums(i,j,kk) = local_sums(i,j,kk)    &
              + data_in3(i,j,k,kk,1) * psi_n(i,j,k)
          END DO
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    CALL global_2d_sums(local_sums, tdims%i_len,       &
                                    tdims%j_len,       &
                                    0, 0, number_of_inputs, mass_n)
    DO kk = 1, number_of_inputs
      IF (mass_n(kk) > 0.0) THEN
        errsl(kk) =      &
             (mass_sl(kk)     - mass_n(kk) - pseudo_lbflux(kk)) / mass_n(kk)
        errslmono(kk) =  &
             (mass_slmono(kk) - mass_n(kk) - pseudo_lbflux(kk)) / mass_n(kk)
      ELSE
        errsl(kk) = 0.0
        errslmono(kk) = 0.0
      END IF
    END DO
    WRITE(umMessage, '(A23,A11,2A30)')         &
             '                      |',        &
             '  species |',                    &
             ' M_SL^{n+1}-(M^{n}+PLF)      |', &
             ' M_SL-QMSL^{n+1}-(M^{n}+PLF) |'
    CALL umPrint(umMessage,src='eg_interpolation_eta_pseudo_lbflux')
    DO kk = 1, number_of_inputs
      WRITE(umMessage, '(A23,I9,A2,2(E28.20,A2))') &
            "Estimated mass errors |",             &
            kk, ' |',                              &
            errsl(kk), ' |',                       &
            errslmono(kk), ' |'
      CALL umPrint(umMessage,src='eg_interpolation_eta_pseudo_lbflux')
    END DO
  END IF ! (printstatus > prstatus_normal)


  DEALLOCATE(data_in3)
  DEALLOCATE(data_out3)
  DEALLOCATE(qs_n)
  DEALLOCATE(qs_s)
  DEALLOCATE(qs_np1)
  DEALLOCATE(psi_n)
  DEALLOCATE(psi_np1)
  DEALLOCATE(qsmin)
  DEALLOCATE(local_sums)
  DEALLOCATE(alfa_za)
  DEALLOCATE(beta_za)
ELSE
  pseudo_lbflux(:) = 0.0
END IF


!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_interpolation_eta_pseudo_lbflux
END MODULE eg_interpolation_eta_pseudo_lbflux_mod
