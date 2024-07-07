! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Print_diag
MODULE print_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PRINT_DIAG_MOD'

CONTAINS

SUBROUTINE Print_diag_4A(                                                    &
    u, v, theta, rho_r2, w, q, qcl, qcf,                                     &
    rows, n_rows, row_length, model_levels,                                  &
    off_x, off_y,timestep_number,                                            &
    rpemax, rpemin, ipesum, rpesum,                                          &
    max_w_run, max_wind_run, min_theta1_run,                                 &
    dtheta1_run, max_div_run, min_div_run,                                   &
    min_lapse_run, max_shear_run, time_max_shear,                            &
    time_div_max, time_div_min, time_lapse_min,                              &
    time_w_max, time_max_wind, time_theta1_min,                              &
    max_KE_run, min_KE_run, max_noise_run,                                   &
    time_KE_max, time_KE_min, time_noise_max )

! Purpose:
!          Diagnostic routine for w, divergence, lapse_rates
!             and level 1 theta
!          v-at-poles version
!
! Method:
!          Is described in ;
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE trignometric_mod, ONLY: sec_theta_latitude,                              &
    cos_theta_latitude, cos_v_latitude,                                      &
    cos_theta_longitude, sin_theta_longitude


USE um_parvars, ONLY: at_extremity,datastart
USE um_parcore, ONLY: mype, nproc
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

USE horiz_grid_mod, ONLY: delta_xi1, delta_xi2
USE level_heights_mod

USE mpp_conf_mod, ONLY: swap_field_is_vector, swap_field_is_scalar

USE atm_fields_bounds_mod, ONLY:                                           &
    udims, vdims, pdims, udims_s, vdims_s, tdims_s, pdims_s,               &
    tdims_l,wdims_s

USE turb_diff_mod, ONLY:                                                   &
    L_print_pe, L_print_w, L_print_wmax, L_print_max_wind, L_print_div,    &
    L_print_lapse, L_print_theta1, L_print_shear, L_diag_wind,             &
    L_diag_noise, print_step, diag_interval, w_print_limit

USE conversions_mod, ONLY: pi, rad_to_deg => recip_pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
USE um_parparams, ONLY: pnorth, psouth
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
USE umPrintMgr, ONLY:                                                      &
    umPrint,                                                               &
    umMessage

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.

INTEGER , INTENT(IN) ::                                                    &
    rows, n_rows, row_length,off_x,off_y, model_levels

INTEGER , INTENT(IN) ::                                                    &
    timestep_number,                                                       &
    rpemax,                                                                &
              ! total size of array needed to find max over pe's
    rpemin,                                                                &
              ! total size of array needed to find min over pe's
    ipesum,                                                                &
              ! total size of integer array needed to sum over pe's
    rpesum    ! total size of real array needed to sum over pe's


INTEGER , INTENT(INOUT) ::                                                 &
    time_theta1_min,                                                       &
    time_ke_max(model_levels + 1),                                         &
    time_ke_min(model_levels + 1),                                         &
    time_w_max(model_levels),                                              &
    time_max_shear(model_levels),                                          &
    time_max_wind(model_levels),                                           &
    time_div_max(model_levels),                                            &
    time_div_min(model_levels),                                            &
    time_lapse_min(model_levels),                                          &
    time_noise_max(model_levels)

REAL , INTENT(INOUT) ::                                                    &
    min_theta1_run,                                                        &
    dtheta1_run,                                                           &
    max_ke_run(model_levels + 1),                                          &
    min_ke_run(model_levels + 1),                                          &
    max_w_run(0:model_levels),                                             &
    max_shear_run(model_levels),                                           &
    max_wind_run(model_levels),                                            &
    max_div_run(model_levels),                                             &
    min_div_run(model_levels),                                             &
    min_lapse_run(model_levels),                                           &
    max_noise_run(model_levels)


REAL , INTENT(IN) ::                                                       &
    u(udims_s%i_start:udims_s%i_end,                                       &
    udims_s%j_start:udims_s%j_end,                                         &
    udims_s%k_start:udims_s%k_end),                                        &
    v(vdims_s%i_start:vdims_s%i_end,                                       &
    vdims_s%j_start:vdims_s%j_end,                                         &
    vdims_s%k_start:vdims_s%k_end),                                        &
    w(wdims_s%i_start:wdims_s%i_end,                                       &
    wdims_s%j_start:wdims_s%j_end,                                         &
    wdims_s%k_start:wdims_s%k_end),                                        &
    theta(tdims_s%i_start:tdims_s%i_end,                                   &
    tdims_s%j_start:tdims_s%j_end,                                         &
    tdims_s%k_start:tdims_s%k_end),                                        &
    rho_r2(pdims_s%i_start:pdims_s%i_end,                                  &
    pdims_s%j_start:pdims_s%j_end,                                         &
    pdims_s%k_start:pdims_s%k_end),                                        &
    q(tdims_l%i_start:tdims_l%i_end,                                       &
    tdims_l%j_start:tdims_l%j_end,                                         &
    tdims_l%k_start:tdims_l%k_end),                                        &
    qcl(tdims_l%i_start:tdims_l%i_end,                                     &
    tdims_l%j_start:tdims_l%j_end,                                         &
    tdims_l%k_start:tdims_l%k_end),                                        &
    qcf(tdims_l%i_start:tdims_l%i_end,                                     &
    tdims_l%j_start:tdims_l%j_end,                                         &
    tdims_l%k_start:tdims_l%k_end)

! Local Variables.

INTEGER ::                                                                 &
    i, j, k,                                                               &
                 ! Loop counters
    gi, gj,                                                                &
                ! global pointers
    ki, kr,                                                                &
                ! pointers for arrays summed over pe's
    kminr,                                                                 &
               ! pointers for arrays summed over pe's
    kmaxr,                                                                 &
               ! pointers for arrays summed over pe's
    j_start_l, j_stop_l,                                                   &
                            ! Loop indices
    info,                                                                  &
    SIZE,                                                                  &
             !  size needed for length for array sums in stats
    level_max,                                                             &
    pe,                                                                    &
    pe_max,                                                                &
    level_max_run,                                                         &
    time_max_run,                                                          &
    i_min_theta,                                                           &
                    ! pointers for theta prints level 1
    j_min_theta     ! pointers for theta prints level 1

REAL ::                                                                    &
    recip_delta_xi1,                                                       &
    recip_delta_xi2,                                                       &
    sum_n,                                                                 &
    sum_s,                                                                 &
    min_theta_pe,                                                          &
                 ! min theta level 1 on this pe
    lambda_max,                                                            &
    phi_max,                                                               &
    max_w,                                                                 &
    max_run,                                                               &
    wind,                                                                  &
    dvdz_at_u,                                                             &
    dudz_at_v,                                                             &
    shear,                                                                 &
    dtheta1,                                                               &
    dtheta1n, dtheta1s, dtheta1e, dtheta1w

! Local arrays

REAL, ALLOCATABLE ::                                                       &
    work1(:,:),                                                            &
    dudz(:,:),                                                             &
    dvdz(:,:)

REAL, ALLOCATABLE ::                                                       &
    u_rho(:,:,:),                                                          &
    v_rho(:,:,:),                                                          &
    w_rho(:,:,:),                                                          &
    rho_dry(:,:,:)

REAL ::                                                                    &
    sumr(rpesum),                                                          &
                 ! array  for summing integers over pe's
    l_n_poles(row_length),                                                 &
    l_s_poles(row_length),                                                 &
    mag_vector_np (model_levels),                                          &
    dir_vector_np (model_levels),                                          &
    mag_vector_sp (model_levels),                                          &
    dir_vector_sp (model_levels),                                          &
    max_wind_pe(model_levels),                                             &
                              ! max wind each level on this pe
    max_w_pe(0:model_levels - 1),                                          &
                               ! max w each level on this pe
    max_shear_pe(model_levels - 1),                                        &
                                   ! max w each level on this pe
    min_lap_pe(model_levels),                                              &
                             ! min lapse rate on this pe
    max_div_pe(model_levels),                                              &
                             ! max div rate on this pe
    min_div_pe(model_levels),                                              &
                             ! min div rate on this pe
    max_real(rpemax),                                                      &
                     ! array  for finding max over pe's
    min_real(rpemin),                                                      &
                     ! array  for finding min over pe's
    w_mean(model_levels - 1),                                              &
    w_variance(model_levels - 1),                                          &
    w_std_dev(model_levels - 1),                                           &
    w_skew(model_levels - 1),                                              &
    w_kurtosis(model_levels - 1)

INTEGER ::                                                                 &
    sumi(ipesum),                                                          &
                 ! array  for summing integers over pe's
    i_min_lap(model_levels),                                               &
                            ! i pointers for prints
    j_min_lap(model_levels),                                               &
                            ! j pointers for prints
    i_max_div(model_levels),                                               &
                            ! i pointers for prints
    j_max_div(model_levels),                                               &
                            ! j pointers for prints
    i_min_div(model_levels),                                               &
                            ! i pointers for prints
    j_min_div(model_levels),                                               &
                            ! j pointers for prints
    i_max_wind(model_levels),                                              &
                             ! i pointers for prints
    j_max_wind(model_levels),                                              &
                             ! j pointers for prints
    i_max_shear(model_levels - 1),                                         &
                                  ! i pointers for prints
    j_max_shear(model_levels - 1),                                         &
                                  ! j pointers for prints
    i_max_w(0:model_levels - 1),                                           &
                              ! i pointers for prints
    j_max_w(0:model_levels - 1),                                           &
                              ! j pointers for prints
    w_count(0:model_levels - 1),                                           &
                              ! counter for w > w_print_limit
    lap_count(model_levels) ! counter for lapse rate < 0

INTEGER :: min_indices(2),max_indices(2)


REAL ::                                                                    &
    recip_row_length,                                                      &
    recip_rows,                                                            &
    weight1,                                                               &
    lambda,                                                                &
                    ! longitude of max
    phi,                                                                   &
                    ! latitude of  max
    sum_levels_ke

CHARACTER(LEN=8) ::                                                        &
    l_string,                                                              &
    p_string,                                                              &
    string_lon,                                                            &
    string_lat

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_DIAG_4A'


! No External Routines:

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( MOD( timestep_number , diag_interval )  ==  0) THEN

  ! ----------------------------------------------------------------------
  ! Section 0. Initialise values
  ! ----------------------------------------------------------------------
  recip_row_length = 1.0 / REAL(global_row_length)
  recip_rows = 1.0 / REAL(global_rows)
  recip_delta_xi1 = 1.0 / delta_xi1
  recip_delta_xi2 = 1.0 / delta_xi2

  IF ( L_print_div ) THEN
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(rho_r2, row_length,rows, model_levels,                &
        off_x, off_y, fld_type_p, swap_field_is_scalar)
  END IF

  IF ( L_print_div .OR. L_print_shear .OR. L_diag_wind                     &
      .OR. L_print_max_wind) THEN
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(v, row_length,n_rows, model_levels,                   &
        off_x, off_y, fld_type_v,swap_field_is_vector)
    ! DEPENDS ON: swap_bounds
    CALL swap_bounds(u, row_length, rows, model_levels,                    &
        off_x,off_y,fld_type_u,swap_field_is_vector)
  END IF

  ! ----------------------------------------------------------------------
  ! Section 1  !  Initialise counters
  ! ----------------------------------------------------------------------
  kminr = 0
  kmaxr = 0

  ! ----------------------------------------------------------------------
  ! Section 2  Find theta min on level 1
  ! ----------------------------------------------------------------------
  IF ( L_print_theta1 ) THEN

    min_indices = MINLOC(theta(1:row_length,1:rows,1))
    i_min_theta = min_indices(1)
    j_min_theta = min_indices(2)
    min_theta_pe = theta(i_min_theta,j_min_theta,1)
    ! copy into min_real(1) which will hold global min later
    kminr = kminr + 1
    min_real(kminr) = min_theta_pe

    ! find dtheta's surrounding minimum (all will be negative)
    dtheta1e = theta(i_min_theta,j_min_theta,1) -                          &
        theta(i_min_theta+1,j_min_theta,1)
    dtheta1w = theta(i_min_theta,j_min_theta,1) -                          &
        theta(i_min_theta-1,j_min_theta,1)
    dtheta1n = theta(i_min_theta,j_min_theta,1) -                          &
        theta(i_min_theta,j_min_theta+1,1)
    dtheta1s = theta(i_min_theta,j_min_theta,1) -                          &
        theta(i_min_theta,j_min_theta-1,1)

  END IF ! L_print_theta1

  ! ----------------------------------------------------------------------
  ! Section 3. Loop over levels > 2 for static stability
  ! ----------------------------------------------------------------------

  IF ( L_print_lapse ) THEN

    ALLOCATE (work1(row_length,rows))

    lap_count(1) = 0
    DO k = 2, model_levels

      lap_count(k) = 0

      DO j = 1, rows
        DO i = 1, row_length
          work1(i,j) = (theta(i,j,k) - theta(i,j,k-1)) /                   &
              (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
          IF ( work1(i,j) < 0.0 ) THEN
            lap_count(k) = lap_count(k) + 1
          END IF ! work1(i,j) < 0.0
        END DO
      END DO

      !Find min lapse_rates on each processor
      min_indices  = MINLOC(work1(1:row_length,1:rows))
      i_min_lap(k) = min_indices(1)
      j_min_lap(k) = min_indices(2)
      min_lap_pe(k) = work1(i_min_lap(k),j_min_lap(k))

      ! Copy local min lap into global max/min_real which will hold
      !          global min after calling gc_rmin
      kminr = kminr + 1
      min_real(kminr) = min_lap_pe(k)

    END DO  ! k = 2, model_levels

    DEALLOCATE ( work1 )

  END IF ! L_print_lapse

  ! ----------------------------------------------------------------------
  ! Section 4. Calculate divergence
  ! ----------------------------------------------------------------------
  IF ( L_print_div ) THEN

    ALLOCATE (work1(row_length,rows))

    DO  k = 1, model_levels

      ! Calculate u derivative at rho points
      DO j = 1, rows
        DO i = 1, row_length
          work1(i,j) = 0.5 * recip_delta_xi1 * ( u(i,j,k) *                &
              ( rho_r2(i,j,k) + rho_r2(i+1,j,k) ) -                        &
              u(i-1,j,k) *                                                 &
              ( rho_r2(i-1,j,k) + rho_r2(i,j,k) ) )*                       &
              sec_theta_latitude(i,j)
        END DO
      END DO

      ! Calculate v derivative at rho points

      DO j = 1, rows
        DO i = 1, row_length
          work1(i,j) = work1(i,j) + 0.5 * recip_delta_xi2 *                &
              ( v(i,j,k) *                                                 &
              ( rho_r2(i,j,k) + rho_r2(i,j+1,k) ) *                        &
              cos_v_latitude(i,j) - v(i,j-1,k) *                           &
              ( rho_r2(i,j,k) + rho_r2(i,j-1,k) ) *                        &
              cos_v_latitude(i,j-1) ) *                                    &
              sec_theta_latitude(i,j)
        END DO
      END DO

      DO j = 1, rows
        DO i = 1, row_length
          work1(i,j) = work1(i,j)                                          &
              / (rho_r2(i,j,k) * r_rho_levels(i,j,k))
        END DO
      END DO

      ! Find maximum absolute divergence on each processor
      min_indices  = MINLOC(work1(1:row_length,1:rows))
      i_min_div(k) = min_indices(1)
      j_min_div(k) = min_indices(2)
      min_div_pe(k) = work1(i_min_div(k),j_min_div(k))

      max_indices  = MAXLOC(work1(1:row_length,1:rows))
      i_max_div(k) = max_indices(1)
      j_max_div(k) = max_indices(2)
      max_div_pe(k) = work1(i_max_div(k),j_max_div(k))

      ! Copy local max/min div into global_max/min_real which will hold
      !    global max/min_div after calling gc_rmax
      kminr = kminr + 1
      kmaxr = kmaxr + 1
      max_real(kmaxr) = max_div_pe(k)
      min_real(kminr) = min_div_pe(k)

    END DO !  k = 1, model_levels

    DEALLOCATE ( work1 )

  END IF ! L_print_div

  ! ----------------------------------------------------------------------
  ! Section 5   Now find max of w for each level
  !               w(model_levels) = 0.0 everywhere so omit
  ! ----------------------------------------------------------------------
  IF (L_print_w .OR. L_print_wmax ) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, w_count, rows, row_length, w, kmaxr,       &
!$OMP         w_print_limit, i_max_w, j_max_w, max_w_pe, max_real )    &
!$OMP PRIVATE( i, j, k, max_indices )
    DO   k = 0, model_levels - 1

      w_count(k) = 0
      DO j = 1, rows
        DO i = 1, row_length
          IF (w(i,j,k) > w_print_limit) THEN
            w_count(k) = w_count(k) + 1
          END IF
        END DO
      END DO

      max_indices = MAXLOC(w(1:row_length,1:rows,k))
      i_max_w(k) = max_indices(1)
      j_max_w(k) = max_indices(2)
      max_w_pe(k) = w(i_max_w(k),j_max_w(k),k)

      ! Copy local max w into w_max which will hold global w_max after gc_rmax
      max_real(k+kmaxr+1) = max_w_pe(k)

    END DO  ! k = 1, model_levels - 1
!$OMP END PARALLEL DO
    kmaxr = kmaxr + model_levels

  END IF ! L_print_w .or. L_print_wmax

  ! ----------------------------------------------------------------------
  ! Section 6   Wind shear diagnostics
  ! ----------------------------------------------------------------------
  IF ( L_print_shear ) THEN

    j_start_l = vdims%j_start
    j_stop_l = vdims%j_end
    IF (model_type == mt_global) THEN
      IF (at_extremity(PSouth)) j_start_l = vdims%j_start +1
      IF (at_extremity(PNorth)) j_stop_l = vdims%j_end -1
    END IF ! model_type  ==  mt_global

    ALLOCATE ( dudz(udims_s%i_start:udims_s%                               &
        i_end,udims_s%j_start:udims_s%j_end) )
    ALLOCATE ( dvdz(vdims_s%i_start:vdims_s%                               &
        i_end,vdims_s%j_start:vdims_s%j_end) )

    DO k = 1, model_levels - 1

      max_shear_pe(k) = 0.0
      i_max_shear(k) = 0
      j_max_shear(k) = 0

      DO j = udims_s%j_start,udims_s%j_end
        DO i = udims_s%i_start,udims_s%i_end
          dudz(i,j) = ( u(i,j,k+1) - u(i,j,k) ) /                          &
              ( r_at_u(i,j,k+1) - r_at_u(i,j,k) )
        END DO
      END DO

      DO j = vdims_s%j_start,vdims_s%j_end
        DO i = vdims_s%i_start,vdims_s%i_end
          dvdz(i,j) = ( v(i,j,k+1) - v(i,j,k) ) /                          &
              ( r_at_v(i,j,k+1) - r_at_v(i,j,k) )
        END DO
      END DO

      DO j = udims%j_start,udims%j_end
        DO i = udims%i_start,udims%i_end
          dvdz_at_u = 0.25 * ( dvdz(i,j) + dvdz(i+1,j) +                   &
              dvdz(i,j-1) + dvdz(i+1,j-1) )
          shear = SQRT( dudz(i,j) * dudz(i,j) +                            &
              dvdz_at_u * dvdz_at_u )
          IF ( shear > max_shear_pe(k) ) THEN
            max_shear_pe(k) = shear
            i_max_shear(k) = i
            j_max_shear(k) = j
          END IF ! shear > max_shear_pe(k)
        END DO
      END DO

      DO j = j_start_l, j_stop_l
        DO i = vdims%i_start,vdims%i_end
          dudz_at_v = 0.25 * ( dudz(i,j  ) + dudz(i-1,j  ) +               &
              dudz(i,j+1) + dudz(i-1,j+1) )
          shear = SQRT( dvdz(i,j) * dvdz(i,j) +                            &
              dudz_at_v * dudz_at_v )
          IF ( shear > max_shear_pe(k) ) THEN
            max_shear_pe(k) = shear
            i_max_shear(k) = i
            j_max_shear(k) = j
          END IF ! shear > max_shear_pe(k)
        END DO
      END DO

      ! Copy local max shear into global_max_real which will hold
      !    global max_shear after calling gc_rmax
      kmaxr = kmaxr + 1
      max_real(kmaxr) = max_shear_pe(k)

    END DO  ! k = 1, model_levels - 1

    DEALLOCATE ( dudz )
    DEALLOCATE ( dvdz )

  END IF ! L_print_shear

  ! ----------------------------------------------------------------------
  ! Section 7.  KE diagnostics and max_wind
  ! ----------------------------------------------------------------------

  !  Initialise pointers in arrays for summing over pe's
  ki = 0
  kr = 0

  IF ( L_diag_wind .OR. L_print_max_wind ) THEN

    ALLOCATE ( u_rho(row_length, rows, model_levels) )
    ALLOCATE ( v_rho(row_length, rows, model_levels) )

    !-------------------------------------------------------------------
    ! 7.1 Interpolate winds to rho points (already on same vertical level)
    !-------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( v, vdims_s, pdims, model_levels, at_extremity, v_rho,    &
!$OMP         u, udims_s, u_rho )
    CALL v_to_p(v,                                                         &
        vdims_s%i_start,vdims_s%i_end,                                     &
        vdims_s%j_start,vdims_s%j_end,                                     &
        pdims%i_start,pdims%i_end,                                         &
        pdims%j_start,pdims%j_end,                                         &
        model_levels, at_extremity,v_rho)

    CALL u_to_p(u,                                                         &
        udims_s%i_start,udims_s%i_end,                                     &
        udims_s%j_start,udims_s%j_end,                                     &
        pdims%i_start,pdims%i_end,                                         &
        pdims%j_start,pdims%j_end,                                         &
        model_levels, at_extremity,u_rho)
!$OMP END PARALLEL

    !-------------------------------------------------------------------
    ! 7.2  Find max wind on this pe
    !-------------------------------------------------------------------
    IF ( L_print_max_wind ) THEN

      ALLOCATE (work1(row_length,rows))

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, max_wind_pe, i_max_wind, j_max_wind, rows, &
!$OMP         row_length, u_rho, v_rho, max_real, kmaxr )              &
!$OMP PRIVATE( i, j, k, work1, max_indices )
      DO k = 1, model_levels
        max_wind_pe(k) = 0.0
        i_max_wind(k) = 0
        j_max_wind(k) = 0
        DO j = 1, rows
          DO i = 1, row_length
            work1(i,j) = SQRT( u_rho(i,j,k) * u_rho(i,j,k) +               &
                v_rho(i,j,k) * v_rho(i,j,k) )
          END DO
        END DO

        max_indices = MAXLOC(work1(1:row_length,1:rows))
        i_max_wind(k) = max_indices(1)
        j_max_wind(k) = max_indices(2)   !! + j_start - 1
        max_wind_pe(k) = work1(i_max_wind(k),j_max_wind(k))

        ! Copy max wind_local into max_real which for global gc_rmax later
        max_real(kmaxr+k) = max_wind_pe(k)

      END DO  !  k = 1, model_levels
!$OMP END PARALLEL DO
      kmaxr = kmaxr + model_levels

      DEALLOCATE (work1)

    END IF ! L_print_max_wind

    !-------------------------------------------------------------------
    ! 7.3  Polar processing for max_wind and wind_diag
    !-------------------------------------------------------------------
    ! Uses v on row next to poles to calculate wind at the polar point

    IF ( .NOT. L_diag_wind ) THEN
      DEALLOCATE ( u_rho )
      DEALLOCATE ( v_rho )
    END IF ! .not. L_diag_wind

  END IF ! L_diag_wind .or. L_print_max_wind

  !-------------------------------------------------------------------
  ! 7.4  Now complete wind_diag
  !-------------------------------------------------------------------
  IF ( L_diag_wind ) THEN

    ALLOCATE ( w_rho(row_length, rows, model_levels) )
    ALLOCATE ( rho_dry(row_length, rows, model_levels) )

    k = 1
    kr = kr + 1
    !  Initialise sumr for ke summing later
    sumr(kr) = 0.0
    DO j = 1, rows
      DO i = 1, row_length
        weight1 = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k)) /          &
            (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
        rho_dry(i,j,k) = rho_r2(i,j,k) *                                   &
            (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
        w_rho(i,j,k) = ( 1.0 -  weight1 ) * w(i,j,k)
      END DO
    END DO

    DO k = 2, model_levels
      kr = kr + 1
      sumr(kr) = 0.0
      DO j = 1, rows
        DO i = 1, row_length
          weight1 = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k)) /        &
              (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
          rho_dry(i,j,k) = rho_r2(i,j,k) * ( (1.0 - weight1) *             &
              (1.0 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k) ) +                &
              weight1 *                                                    &
              (1.0 - q(i,j,k-1) - qcl(i,j,k-1) - qcf(i,j,k-1) ) )
          w_rho(i,j,k) = (1.0 -  weight1) * w(i,j,k) +                     &
              weight1 * w(i,j,k-1)
        END DO
      END DO
    END DO  !  k = 2, model_levels

    ! Problem of u & v values at poles
    ! Uses v on row next to poles to calculate wind at the polar point
    ! reset sum array counter
    kr = kr - model_levels

    !----------------------------------------------------------------------
    ! Integrals over fields using values interpolated to rho grid
    ! At present all integrals done after interpolation to rho points.
    ! Integrals for energy involve using dry mass.
    !----------------------------------------------------------------------
    ! KE terms for u, v & w
    DO k = 1, model_levels

      kr = kr + 1
      DO j = 1, rows
        DO i = 1, row_length
          sumr(kr) = sumr(kr) + 0.5 * rho_dry(i,j,k) *                     &
              delta_xi1 * delta_xi2 *                                      &
              (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) *          &
              ( u_rho(i,j,k) * u_rho(i,j,k) +                              &
              v_rho(i,j,k) * v_rho(i,j,k) +                                &
              w_rho(i,j,k) * w_rho(i,j,k) ) *                              &
              cos_theta_latitude(i,j)
        END DO
      END DO

    END DO  !  k = 1, model_levels

    DEALLOCATE ( u_rho )
    DEALLOCATE ( v_rho )
    DEALLOCATE ( w_rho )
    DEALLOCATE ( rho_dry )

  END IF ! L_diag_wind

  ! Global sum for KE done with all the other sums over processors

  ! ----------------------------------------------------------------------
  ! Section 8. Now find max/mins over all processors
  !            All required fields done at same time
  ! ----------------------------------------------------------------------

  IF ( kmaxr > 0 ) THEN
    CALL gc_rmax(kmaxr, nproc, info, max_real)
  END IF !  kmaxr > 0
  IF ( kminr > 0 ) THEN
    CALL gc_rmin(kminr, nproc, info, min_real)
  END IF !  kminr > 0

  !  Re-initialise pointers in max/min arrays
  kminr = 0
  kmaxr = 0

  ! ----------------------------------------------------------------------
  ! Section 9. Now find time and place of max/mins
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Section 9.1  Obtain max and mins from max_real, min_real
  !              and fill sumi, sumr arrays for summing over pe's
  ! ----------------------------------------------------------------------

  IF ( L_print_theta1 ) THEN
    kminr = kminr + 1
    dtheta1 = MIN(dtheta1e, dtheta1w, dtheta1n, dtheta1s)
    ki = ki + 1
    kr = kr + 3
    sumi(ki) = 0
    sumr(kr-2) = 0.0
    sumr(kr-1) = 0.0
    sumr(kr) = 0.0
    IF ( min_theta_pe <= min_real(kminr) ) THEN
      sumi(ki) = mype ! min is on this processor
      ! store largest negative gradient at min theta1 point
      sumr(kr-2) = dtheta1
      gi = datastart(1) + i_min_theta - 1
      gj = datastart(2) + j_min_theta - 1
      IF (model_type == mt_global) THEN
        sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
        sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
      ELSE ! for LAM output fraction of point
        sumr(kr-1) = (gi-1) * recip_row_length * 100.0
        sumr(kr) = (gj-1) * recip_rows * 100.0
      END IF ! model_type  ==  mt_global
    END IF ! min_theta_pe <= min_real(kminr)

  END IF ! L_print_theta1

  IF ( L_print_lapse ) THEN
    DO k = 2, model_levels
      kminr = kminr + 1
      IF ( min_real(kminr) < min_lapse_run(k) ) THEN
        min_lapse_run(k) = min_real(kminr)
        time_lapse_min(k) = timestep_number
      END IF ! min_real(kminr) > min_lapse_run(k)
      ki = ki + 2
      kr = kr + 2
      sumi(ki - 1) = lap_count(k)
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( min_lap_pe(k) <= min_real(kminr)) THEN
        sumi(ki) = mype ! min is on this processor for this level
        gi = datastart(1) + i_min_lap(k) - 1
        gj = datastart(2) + j_min_lap(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! min_lap_pe(k) <= min_real(kminr)
    END DO ! k = 2, model_levels

  END IF ! L_print_lapse

  IF ( L_print_div ) THEN

    DO k = 1, model_levels
      kmaxr = kmaxr + 1
      IF ( max_real(kmaxr) > max_div_run(k) ) THEN
        max_div_run(k) =  max_real(kmaxr)
        time_div_max(k) = timestep_number
      END IF ! max_real(kminr)  > max_div_run(k)
      ki = ki + 1
      kr = kr + 2
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( max_div_pe(k) >= max_real(kmaxr)) THEN
        sumi(ki) = mype ! max is on this processor for this level
        gi = datastart(1) + i_max_div(k) - 1
        gj = datastart(2) + j_max_div(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! max_div_pe(k) >= max_real(kmaxr)
    END DO  !  k = 1, model_levels

    DO k = 1, model_levels
      kminr = kminr + 1
      IF ( min_real(kminr) < min_div_run(k) ) THEN
        min_div_run(k) =  min_real(kminr)
        time_div_min(k) = timestep_number
      END IF ! min_real(kminr) > min_div_run(k)
      ki = ki + 1
      kr = kr + 2
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( min_div_pe(k) <= min_real(kminr)) THEN
        sumi(ki) = mype ! min is on this processor for this level
        gi = datastart(1) + i_min_div(k) - 1
        gj = datastart(2) + j_min_div(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! min_div_pe(k) <= min_real(kminr)
    END DO  !  k = 1, model_levels

  END IF ! L_print_div

  IF ( L_print_w .OR. L_print_wmax ) THEN
    DO   k = 0, model_levels - 1
      kmaxr = kmaxr + 1
      IF ( max_real(kmaxr) > max_w_run(k) ) THEN
        max_w_run(k) =  max_real(kmaxr)
        time_w_max(k+1) = timestep_number
      END IF ! max_real(kmaxr) > max_w_run(k)
      ki = ki + 2
      kr = kr + 2
      sumi(ki - 1) = w_count(k)
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( max_w_pe(k) >= max_real(kmaxr)) THEN
        sumi(ki) = mype ! max is on this processor for this level
        gi = datastart(1) + i_max_w(k) - 1
        gj = datastart(2) + j_max_w(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-0.5) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! max_w_pe(k) >= max_real(kmaxr)
    END DO  ! k = 1, model_levels - 1
  END IF ! L_print_w .or. L_print_wmax

  IF ( L_print_shear ) THEN
    DO   k = 1, model_levels - 1
      kmaxr = kmaxr + 1
      IF ( max_real(kmaxr) > max_shear_run(k) ) THEN
        max_shear_run(k) =  max_real(kmaxr)
        time_max_shear(k) = timestep_number
      END IF ! max_real(kmaxr) > max_shear_run(k)
      ki = ki + 1
      kr = kr + 2
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( max_shear_pe(k) >= max_real(kmaxr)) THEN
        sumi(ki) = mype ! max is on this processor for this level
        gi = datastart(1) + i_max_shear(k) - 1
        gj = datastart(2) + j_max_shear(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! max_shear_pe(k) >= max_real(kmaxr)
    END DO  ! k = 1, model_levels - 1
  END IF ! L_print_shear

  ! ----------------------------------------------------------------------

  IF ( L_print_max_wind ) THEN
    DO   k = 1, model_levels
      kmaxr = kmaxr + 1
      IF ( max_real(kmaxr) > max_wind_run(k) ) THEN
        max_wind_run(k) =  max_real(kmaxr)
        time_max_wind(k) = timestep_number
      END IF ! max_real(kmaxr) > max_wind_run(k)
      ki = ki + 1
      kr = kr + 2
      sumi(ki) = 0
      sumr(kr-1) = 0.0
      sumr(kr) = 0.0
      IF ( max_wind_pe(k) >= max_real(kmaxr)) THEN
        sumi(ki) = mype ! max is on this processor for this level
        gi = datastart(1) + i_max_wind(k) - 1
        gj = datastart(2) + j_max_wind(k) - 1
        IF (model_type == mt_global) THEN
          sumr(kr-1) = (gi-1) * delta_xi1 * rad_to_deg
          sumr(kr) = ( (gj-1) * delta_xi2 - 0.5 * Pi ) * rad_to_deg
        ELSE ! for LAM output fraction of point
          sumr(kr-1) = (gi-1) * recip_row_length * 100.0
          sumr(kr) = (gj-1) * recip_rows * 100.0
        END IF ! model_type  ==  mt_global
      END IF  ! max_wind_pe(k) <= max_real(kmaxr)
    END DO  ! k = 1, model_levels
  END IF ! L_print_max_wind

  ! ----------------------------------------------------------------------
  !  Printing will be done in section 12
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Section 10  Summing over pe's to obtain sums and location of max/mins
  ! ----------------------------------------------------------------------

  ! We do not need a reproducible sum, as each element of sumr consists of
  ! one non-zero entry and all other entries as zero.
  ! Adding nproc-1 zero's to X give X whether you use a reproducible sum
  ! or not.

  IF ( ki > 0 ) CALL gc_isum(ki, nproc, info, sumi)
  IF ( kr > 0 ) CALL gc_rsum(kr, nproc, info, sumr)

  ! ----------------------------------------------------------------------
  ! Section 10.1  Obtain max/mins for KE
  ! ----------------------------------------------------------------------
  IF ( L_diag_wind ) THEN

    ! Re-initialise pointer
    kr = 0
    sum_levels_ke = 0.0
    DO   k = 1, model_levels
      kr = kr + 1
      IF ( sumr(kr) > max_ke_run(k) ) THEN
        max_ke_run(k) = sumr(kr)
        time_ke_max(k) = timestep_number
      END IF ! sumr(kr) > max_ke_run(k)
      IF ( sumr(kr) < min_ke_run(k) ) THEN
        min_ke_run(k) =  sumr(kr)
        time_ke_min(k) = timestep_number
      END IF ! sumr(kr) > max_ke_run(k)
      sum_levels_ke = sum_levels_ke + sumr(kr)
    END DO   !  k = 1, model_levels

    k = model_levels + 1  ! Put vertical sum in  model_levels + 1
    IF ( sum_levels_ke > max_ke_run(k) ) THEN
      max_ke_run(k) =  sum_levels_ke
      time_ke_max(k) = timestep_number
    END IF ! sum_levels_ke > max_ke_run(k)
    IF ( sum_levels_ke < min_ke_run(k) ) THEN
      min_ke_run(k) =  sum_levels_ke
      time_ke_min(k) = timestep_number
    END IF ! sum_levels_ke < min_ke_run(k)

  END IF ! L_diag_wind

  ! ----------------------------------------------------------------------
  ! Section 11. calculate stats for w field to assess noise
  ! ----------------------------------------------------------------------

  IF ( L_diag_noise ) THEN

    SIZE = 4 * (model_levels - 1)

    ! DEPENDS ON: calc_stats
    CALL calc_stats(                                                       &
        w(1-off_x, 1-off_y, 1),                                            &
        row_length, rows, model_levels - 1,                                &
        off_x, off_y, nproc, SIZE,                                        &
        w_mean, w_variance, w_std_dev,                                     &
        w_skew, w_kurtosis)

  END IF ! L_diag_noise

END IF ! mod( timestep_number , diag_interval )  ==  0

! ----------------------------------------------------------------------
! Section 12. Print diagnostic information
! ----------------------------------------------------------------------

IF ( MOD( timestep_number , print_step ) == 0                               &
    .AND.   (L_print_pe .OR. mype == 0)       ) THEN

  ! Re-initialise pointers
  kminr = 0
  kmaxr = 0
  ki = 0
  kr = 0

  IF ( L_diag_noise ) THEN

    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Vertical velocity characteristics at ',           &
        timestep_number,' time steps'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'Level   Mean      Variance    std_dev   ',         &
        ' Skew      Kurtosis'
    CALL umPrint(umMessage,src='print_diag_4a')

    max_w = -100.0
    max_run = -100.0
    DO k = 1, model_levels - 1

      IF ( w_variance(k) > max_noise_run(k) ) THEN
        max_noise_run(k) = w_variance(k)
        time_noise_max(k) = timestep_number
      END IF ! w_variance(k) > max_noise_run(k)
      IF ( w_variance(k) > max_w ) THEN
        max_w = w_variance(k)
        level_max = k
      END IF ! w_variance(k) > max_noise
      IF ( max_noise_run(k) > max_run ) THEN
        max_run = max_noise_run(k)
        level_max_run = k
        time_max_run = time_noise_max(k)
      END IF ! max_noise_run(k) > max_run

      WRITE(umMessage,FMT='(1x,i4,5e11.3)') k, w_mean(k), w_variance(k),   &
          w_std_dev(k), w_skew(k), w_kurtosis(k)
      CALL umPrint(umMessage,src='print_diag_4a')
    END DO !   k = 1, model_levels - 1

    WRITE(umMessage,*) '    '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Maximum noise (w variance) at timestep ',         &
        timestep_number
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '   Max noise      level',                          &
        '   max this run  level   timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(1x,e12.3,i9,e14.3,2i9)') max_w, level_max,       &
        max_run, level_max_run, time_max_run
    CALL umPrint(umMessage,src='print_diag_4a')

  END IF ! L_diag_noise

  IF ( L_diag_wind ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' KE at ',timestep_number,' time steps'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level   This step     min this run (time step)',   &
        '  max this run  (time step)'
    CALL umPrint(umMessage,src='print_diag_4a')
    DO   k = 1, model_levels
      kr = kr + 1
      WRITE(umMessage,FMT='(1x,i4,2e15.3,i9,e15.3,i9)') k, sumr(kr),       &
          min_ke_run(k), time_ke_min(k), max_ke_run(k), time_ke_max(k)
      CALL umPrint(umMessage,src='print_diag_4a')
    END DO   !  k = 1, model_levels

    k = model_levels + 1  ! Put vertical sum in  model_levels + 1
    WRITE(umMessage,*)' Total KE    Min KE this run at timestep ',         &
        ' Max KE this run at timestep '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(1x,2e15.3,i9,e15.3,i9)') sum_levels_ke,          &
        min_ke_run(k), time_ke_min(k), max_ke_run(k), time_ke_max(k)
    CALL umPrint(umMessage,src='print_diag_4a')
  END IF ! L_diag_wind

  IF ( L_print_theta1 ) THEN
    kminr = kminr + 1
    ki = ki + 1
    kr = kr + 3
    pe = sumi(ki)
    dtheta1 = sumr(kr-2)
    lambda = sumr(kr-1)
    phi = sumr(kr)
    IF ( min_real(kminr) < min_theta1_run ) THEN
      min_theta1_run = min_real(kminr)
      time_theta1_min = timestep_number
      dtheta1_run = dtheta1
    END IF ! min_real(kminr) < min_theta1_run
    IF (model_type == mt_global) THEN
      IF (lambda > 180.0) THEN
        lambda = 360.0 - lambda
        l_string ='deg W'
      ELSE
        l_string ='deg E'
      END IF      !  lambda(k) > 180.0
      IF (phi > 0.0) THEN
        p_string ='deg N'
      ELSE
        p_string ='deg S'
      END IF      !  phi(k) > 0.0
    ELSE ! LAM domains
      p_string ='% North'
      l_string ='% East'
    END IF ! model_type  ==  mt_global
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*)'Minimum theta level 1 for timestep ',timestep_number
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '               This timestep',                     &
        '                         This run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '  Min theta1     proc          position',          &
        '            Min theta1 timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(1x,f11.2,i8,f8.1,a7,f8.1,a7,f11.2,i6)')          &
        min_real(kminr), pe, lambda, l_string, phi, p_string,              &
        min_theta1_run, time_theta1_min
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Largest negative delta theta1 at minimum theta1 '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(a17,f8.2,a20,f8.2,a1)') ' This timestep = ',     &
        dtheta1,'K. At min for run = ', dtheta1_run,'K'
    CALL umPrint(umMessage,src='print_diag_4a')

  END IF ! L_print_theta1

  IF ( L_print_lapse ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' *****  Minimum lapse_rate at timestep ',          &
        timestep_number, '     ****** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '      number          this timestep',              &
        '                      this run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Level <0   dtheta/dz  proc       position',       &
        '            run min    timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    DO k = 2, model_levels
      kminr = kminr + 1
      ki = ki + 2
      kr = kr + 2
      lap_count(1) = sumi(ki - 1)
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda_min(k) > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi_min(k) > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global

      WRITE(umMessage,FMT='(1x,i4,i6,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)')  &
          k, lap_count(1), min_real(kminr), pe, lambda, l_string, phi,     &
          p_string, min_lapse_run(k), time_lapse_min(k)
      CALL umPrint(umMessage,src='print_diag_4a')

    END DO ! k = 2, model_levels

  END IF ! L_print_lapse

  IF ( L_print_div ) THEN

    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' ***** Maximum divergence at timestep ',           &
        timestep_number,' ****** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '         this timestep',                           &
        '                           this run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level  max div   proc     and   position',         &
        '      max div at timestep'
    CALL umPrint(umMessage,src='print_diag_4a')

    DO k = 1, model_levels
      kmaxr = kmaxr + 1
      ki = ki + 1
      kr = kr + 2
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global

      WRITE(umMessage,FMT='(1x,i4,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)')     &
          k, max_real(kmaxr), pe, lambda, l_string, phi, p_string,         &
          max_div_run(k), time_div_max(k)
      CALL umPrint(umMessage,src='print_diag_4a')

    END DO  !  k = 1, model_levels

    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*)                                                     &
        ' **** Minimum divergence (CONVERGENCE) at timestep ',             &
        timestep_number,' ***** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '         this timestep',                           &
        '                           this run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level  min div   proc     and   position',         &
        '      min div at timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    DO k = 1, model_levels
      kminr = kminr + 1
      ki = ki + 1
      kr = kr + 2
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global
      WRITE(umMessage,FMT='(1x,i4,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)')     &
          k, min_real(kminr), pe, lambda, l_string, phi, p_string,         &
          min_div_run(k), time_div_min(k)
      CALL umPrint(umMessage,src='print_diag_4a')
    END DO  !  k = 1, model_levels

  END IF ! L_print_div

  IF ( L_print_w ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' ***** Maximum vertical velocity w at timestep ',  &
        timestep_number,' ****** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(a51,ES10.3,a4)')                                   &
     ' The column below w_limit shows number of points > ',w_print_limit,' m/s'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '         pts >        this timestep',              &
        '                     this run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level w_limit  w_max   proc     and   position',   &
        '        w_max at timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
  END IF ! L_print_w
  IF ( L_print_w .OR. L_print_wmax ) THEN
    max_w = -100.0
    max_run = -100.0
    DO   k = 0, model_levels - 1
      kmaxr = kmaxr + 1
      ki = ki + 2
      kr = kr + 2
      w_count(1) = sumi(ki - 1)
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global
      IF ( max_w < max_real(kmaxr) ) THEN
        max_w = max_real(kmaxr)
        level_max = k
        pe_max = pe
        lambda_max = lambda
        phi_max = phi
        string_lon = l_string
        string_lat = p_string
      END IF ! max_w < max_real(kmaxr)
      IF ( max_run < max_w_run(k) ) THEN
        max_run = max_w_run(k)
        level_max_run = k
        time_max_run = time_w_max(k+1)
      END IF ! max_run < max_w_run(k)
      IF ( L_print_w ) THEN
        WRITE(umMessage,FMT='(1x,i4,i6,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)') &
            k, w_count(1), max_real(kmaxr), pe, lambda, l_string, phi,      &
            p_string, max_w_run(k), time_w_max(k+1)
        CALL umPrint(umMessage,src='print_diag_4a')
      END IF ! L_print_w
    END DO   !  k = 1, model_levels - 1
  END IF ! L_print_w .or. L_print_wmax

  IF ( L_print_wmax ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Maximum vertical velocity at timestep ',          &
        timestep_number,'      Max w this run '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '   w_max   level  proc         position        ',  &
        '     run w_max level timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(1x,e11.3,i4,i7,f7.1,a7,f7.1,a7,e11.3,i5,i6)')    &
        max_w, level_max, pe_max, lambda_max, string_lon, phi_max,         &
        string_lat, max_run, level_max_run, time_max_run
    CALL umPrint(umMessage,src='print_diag_4a')
  END IF ! L_print_wmax

  IF ( L_print_shear ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' ***** Maximum vertical wind shear at timestep ',  &
        timestep_number,' ****** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level limit  max   proc     and   position',       &
        '        max shear at timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    max_w = -100.0
    max_run = -100.0
    DO   k = 1, model_levels - 1
      kmaxr = kmaxr + 1
      ki = ki + 1
      kr = kr + 2
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global
      WRITE(umMessage,FMT='(1x,i4,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)')     &
          k, max_real(kmaxr), pe, lambda, l_string, phi, p_string,         &
          max_shear_run(k), time_max_shear(k)
      CALL umPrint(umMessage,src='print_diag_4a')
    END DO   !  k = 1, model_levels - 1
  END IF ! L_print_shear

  IF ( L_print_max_wind ) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' ***** Maximum wind speed at timestep ',           &
        timestep_number,' ****** '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '                this timestep',                    &
        '                     this run'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) 'level  max wind   proc     and   position',        &
        '      max wind at timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    max_w = 0.0
    max_run = 0.0
    DO   k = 1, model_levels
      kmaxr = kmaxr + 1
      ki = ki + 1
      kr = kr + 2
      pe = sumi(ki)
      lambda = sumr(kr-1)
      phi = sumr(kr)
      IF (model_type == mt_global) THEN
        IF (lambda > 180.0) THEN
          lambda = 360.0 - lambda
          l_string ='deg W'
        ELSE
          l_string ='deg E'
        END IF      !  lambda > 180.0
        IF (phi > 0.0) THEN
          p_string ='deg N'
        ELSE
          p_string ='deg S'
        END IF      !  phi > 0.0
      ELSE ! LAM domains
        p_string ='% North'
        l_string ='% East'
      END IF ! model_type  ==  mt_global
      IF ( max_w < max_real(kmaxr) ) THEN
        max_w = max_real(kmaxr)
        level_max = k
        pe_max = pe
        lambda_max = lambda
        phi_max = phi
        string_lon = l_string
        string_lat = p_string
      END IF ! max_w < max_real(kmaxr)
      IF ( max_run < max_wind_run(k) ) THEN
        max_run = max_wind_run(k)
        level_max_run = k
        time_max_run = time_max_wind(k)
      END IF ! max_run < max_wind_run(k)
      WRITE(umMessage,FMT='(1x,i4,e11.3,i7,f6.1,a7,f6.1,a7,e11.3,i6)')     &
          k,  max_real(kmaxr), pe, lambda, l_string, phi, p_string,        &
          max_wind_run(k), time_max_wind(k)
      CALL umPrint(umMessage,src='print_diag_4a')
    END DO   !  k = 1, model_levels - 1

    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) ' Maximum horizontal wind at timestep ',            &
        timestep_number,'      Max wind this run '
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,*) '   max_wind   level  proc         position   ',    &
        '     run max_wind level timestep'
    CALL umPrint(umMessage,src='print_diag_4a')
    WRITE(umMessage,FMT='(1x,e11.3,i4,i7,f7.1,a7,f7.1,a7,e11.3,i5,i6)')    &
        max_w, level_max, pe_max, lambda_max, string_lon, phi_max,         &
        string_lat, max_run, level_max_run, time_max_run
    CALL umPrint(umMessage,src='print_diag_4a')
  END IF ! L_print_max_wind

END IF !  mod(timestep_number, print_step) == 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Print_diag_4A
END MODULE
