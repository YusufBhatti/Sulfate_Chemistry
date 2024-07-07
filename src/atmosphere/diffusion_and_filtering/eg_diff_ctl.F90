! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE eg_diff_ctl
SUBROUTINE eg_diff_ctl(                                                 &
                       L_Backwards,                                     &
                       pos_timestep, neg_timestep,                      &
                       theta, moist, u, v, w,                           &
                       m_cl_np1, m_cf_np1, m_r_np1,                     &
                       m_cf2_np1, m_gr_np1,                             &
                       r_theta_levels, r_rho_levels,                    &
                       offx, offy, halo_i, halo_j,                      &
                       row_length, rows, n_rows,                        &
                       model_levels, flat_level,                        &
                       theta_star, moist_star, S_u, S_v, S_w,           &
                       s_m_cl, s_m_cf, s_m_r, s_m_cf2, s_m_gr)

! Purpose:
!          Subroutine to interface to diffusion code
!
! Method:
!          Is described in ;
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE atm_fields_bounds_mod, ONLY: udims,vdims, udims_s, vdims_s, pdims,  &
     pdims_s, pdims_l, tdims, tdims_s, tdims_l
USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE timestep_mod, ONLY: timestep
USE metric_terms_mod, ONLY: h2_p_eta, h2_xi1_u, h2_xi2_v
USE horiz_grid_mod, ONLY: xi1_p, xi1_u, xi2_p, xi2_v,                   &
                          csxi2_p, csxi2_v, intw_p2u, intw_p2v
USE turb_diff_mod, ONLY: turb_startlev_horiz,turb_endlev_horiz,         &
                         hyd_mix_opt, liq_only, liq_and_ice, all_hydro
USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, max_diff
USE mphys_inputs_mod, ONLY: l_mcr_qrain, l_mcr_qgraup, l_mcr_qcf2

USE swapable_field_mod, ONLY: swapable_field_pointer_type
USE mpp_conf_mod,       ONLY: swap_field_is_scalar, swap_field_is_vector

USE UM_ParVars, ONLY: at_extremity
USE UM_ParParams, ONLY: pnorth, psouth
USE Field_Types,  ONLY: fld_type_p, fld_type_u, fld_type_v
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

LOGICAL, INTENT(IN) ::                                                  &
  L_Backwards

INTEGER, INTENT(IN) ::                                                  &
  row_length                                                            &
                   ! number of point on a row.
, rows                                                                  &
                   ! number of rows.
, n_rows                                                                &
                   ! number of v rows.
, model_levels                                                          &
                   ! number of model levels.
, flat_level                                                            &
                   ! first flat level
, offx                                                                  &
               ! Size of small halo in i
, offy                                                                  &
               ! Size of small halo in j.
, halo_i                                                                &
                 ! Size of halo in i direction.
, halo_j         ! Size of halo in j direction.

REAL, INTENT(IN) ::                                                     &
 pos_timestep                                                           &
              ! = +timestep.
,neg_timestep
              ! = -timestep.

REAL, INTENT(IN) ::                                                     &
     ! primary model variables
  u(udims_s%i_start:udims_s%i_end,                                      &
    udims_s%j_start:udims_s%j_end, model_levels )                       &
, v(vdims_s%i_start:vdims_s%i_end,                                      &
    vdims_s%j_start:vdims_s%j_end, model_levels )                       &
, w(tdims_s%i_start:tdims_s%i_end,                                      &
    tdims_s%j_start:tdims_s%j_end, 0:model_levels)                      &
, moist (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end, 0:model_levels)                 &
, m_cl_np1 (tdims_s%i_start:tdims_s%i_end,                              &
            tdims_s%j_start:tdims_s%j_end, 0:model_levels)              &
, m_cf_np1 (tdims_s%i_start:tdims_s%i_end,                              &
            tdims_s%j_start:tdims_s%j_end, 0:model_levels)              &
, m_r_np1 (tdims_s%i_start:tdims_s%i_end,                               &
           tdims_s%j_start:tdims_s%j_end, 0:model_levels)               &
, m_gr_np1 (tdims_s%i_start:tdims_s%i_end,                              &
            tdims_s%j_start:tdims_s%j_end, 0:model_levels)              &
, m_cf2_np1 (tdims_s%i_start:tdims_s%i_end,                             &
             tdims_s%j_start:tdims_s%j_end, 0:model_levels)             &
, theta (tdims_s%i_start:tdims_s%i_end,                                 &
         tdims_s%j_start:tdims_s%j_end, 0:model_levels)

REAL, INTENT(IN) ::                                                     &
     ! vertical co-ordinate arrays.
  r_theta_levels (tdims_l%i_start:tdims_l%i_end,                        &
                  tdims_l%j_start:tdims_l%j_end, 0:model_levels)        &
, r_rho_levels (pdims_l%i_start:pdims_l%i_end,                          &
                pdims_l%j_start:pdims_l%j_end, model_levels)

REAL, ALLOCATABLE ::                                                    &
 coeff_u(:,:,:)                                                         &
                       ! visc_m interpolated onto
,coeff_v (:,:,:)                                                        &
                       ! appropriate points
,coeff_th(:,:,:)                                                        &
                       ! for use in the turb_diff_*
,work_th(:,:,:)                                                         &
                       ! for use in the turb_diff_*
,coeff_centre(:,:,:)                                                    &
                       ! diffusion routines.
,w_coeff_u(:,:,:)                                                       &
                       ! horizontal diffusion coefficients for
,w_coeff_v(:,:,:)        ! diffusion of w

! local
INTEGER :: i, j, k, l

!   variables with intent INOUT
!   For ENDGame, these are Fast physics fields
REAL, INTENT(INOUT) ::                                                  &
  S_u(udims_s%i_start:udims_s%i_end,                                    &
      udims_s%j_start:udims_s%j_end, model_levels )                     &
, S_v(vdims_s%i_start:vdims_s%i_end,                                    &
      vdims_s%j_start:vdims_s%j_end, model_levels )                     &
, S_w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:model_levels) &
, moist_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end, 0:model_levels)                 &
, s_m_cl(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end, 0:model_levels)                     &
, s_m_cf(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end, 0:model_levels)                     &
, s_m_r(tdims%i_start:tdims%i_end,                                      &
        tdims%j_start:tdims%j_end, 0:model_levels)                      &
, s_m_gr(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end, 0:model_levels)                     &
, s_m_cf2(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end, 0:model_levels)                    &
, theta_star(tdims_s%i_start:tdims_s%i_end,                             &
             tdims_s%j_start:tdims_s%j_end, 0:model_levels)

! Local variables
INTEGER, PARAMETER :: zero_halo = 0
INTEGER :: j_start, j_stop    ! j start/end points

INTEGER :: i_field

REAL :: weight ! vertical weight

REAL ::  delta_z(pdims_s%i_start:pdims_s%i_end,                         &
                 pdims_s%j_start:pdims_s%j_end, model_levels)
REAL ::  recip_r_squared_delz(pdims_s%i_start:pdims%i_end,              &
                              pdims_s%j_start:pdims%j_end, model_levels)
REAL ::  r_theta_uv(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end, flat_level)

TYPE(swapable_field_pointer_type) :: fields_to_swap(2)
                                     ! multivariate swapbounds

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DIFF_CTL'

! ----------------------------------------------------------------------
! Section 1.0  Horizontal Diffusion
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

j_start = vdims%j_start
j_stop = vdims%j_end
IF (at_extremity(PSouth)) j_start = j_start+1
IF (at_extremity(PNorth)) j_stop = j_stop-1

ALLOCATE( coeff_u(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( coeff_v(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( coeff_th(pdims_s%i_start:pdims_s%i_end,                       &
                   pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( work_th(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( coeff_centre(pdims_s%i_start:pdims_s%i_end,                   &
                       pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( w_coeff_u(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end, model_levels))
ALLOCATE( w_coeff_v(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end, model_levels))

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                           &
!$OMP SHARED(model_levels, pdims, coeff_u, coeff_v,                     &
!$OMP    coeff_th, work_th, coeff_centre, w_coeff_u, w_coeff_v,         &
!$OMP    turb_startlev_horiz, turb_endlev_horiz,  visc_h, visc_m,       &
!$OMP    max_diff)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  coeff_u(:,:,k) = 0.0
  coeff_v(:,:,k) = 0.0
  coeff_th(:,:,k) = 0.0
  work_th(:,:,k) = 0.0
  coeff_centre(:,:,k) = 0.0
  w_coeff_u(:,:,k) = 0.0
  w_coeff_v(:,:,k) = 0.0
END DO
!$OMP END DO NOWAIT

DO k = 1, model_levels - 2

  IF (k >= turb_startlev_horiz .AND. k <= turb_endlev_horiz) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_h(i,j,k) = MIN(visc_h(i,j,k), max_diff(i,j))
        visc_m(i,j,k) = MIN(visc_m(i,j,k), max_diff(i,j))
      END DO
    END DO
!$OMP END DO NOWAIT

  ELSE

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_h(i,j,k) = 0.0
        visc_m(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF

END DO

! Set visc_m and visc_h at the top two levels

DO k = model_levels - 1, model_levels

  IF (turb_startlev_horiz <= k .AND. turb_endlev_horiz >= k) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_h(i,j,k) = visc_h(i,j,model_levels-2)
        visc_m(i,j,k) = visc_m(i,j,model_levels-2)
      END DO
    END DO
!$OMP END DO NOWAIT

  ELSE

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_h(i,j,k) = 0.0
        visc_m(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF

END DO

!$OMP END PARALLEL

i_field = 0

i_field = i_field + 1
fields_to_swap(i_field) % field      => visc_m(:,:,:)
fields_to_swap(i_field) % field_type =  fld_type_p
fields_to_swap(i_field) % levels     =  model_levels
fields_to_swap(i_field) % rows       =  rows
fields_to_swap(i_field) % vector     =  .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field      => visc_h(:,:,:)
fields_to_swap(i_field) % field_type =  fld_type_p
fields_to_swap(i_field) % levels     =  model_levels
fields_to_swap(i_field) % rows       =  rows
fields_to_swap(i_field) % vector     =  .FALSE.

CALL swap_bounds_mv( fields_to_swap, i_field,                           &
                               row_length, halo_i, halo_j)

!        Run with a positive timestep IF integrating backwards.
IF (L_Backwards) timestep = pos_timestep
  !
  !Interpolate visc_m onto appropriate points for the diffusion routines.
  ! At this point the diffusion coefficients are on w (theta) points.

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k, weight)                   &
!$OMP SHARED(model_levels, pdims, pdims_s, flat_level, work_th,         &
!$OMP    visc_m, r_rho_levels, r_theta_levels, coeff_u, visc_h,         &
!$OMP    w_coeff_u, coeff_v, w_coeff_v, coeff_centre, coeff_th,         &
!$OMP    delta_z, recip_r_squared_delz, h2_p_eta, intw_p2u, intw_p2v)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    work_th(i,j,1) = visc_m(i,j,1)
  END DO
END DO
!$OMP END DO NOWAIT

DO k = 2, model_levels
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      weight = ( r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1) ) /      &
               (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
! visc_m on rho-levels
      work_th(i,j,k) = ( weight * visc_m(i,j,k) +                       &
                       (1.0 - weight) * visc_m(i,j,k-1) )
    END DO
  END DO
!$OMP END DO NOWAIT
END DO ! k = 2, model_levels

! Flush work_th
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start-1, pdims%i_end
! visc_h on u points, theta-levels
      coeff_u(i,j,k) = intw_p2u(i,1) * visc_h(i,j,k) +                  &
                       intw_p2u(i,2) * visc_h(i+1,j,k)
! visc_m on u points, theta-levels
      w_coeff_u(i,j,k) = intw_p2u(i,1) * visc_m(i,j,k) +                &
                         intw_p2u(i,2) * visc_m(i+1,j,k)
    END DO
  END DO

  DO j = pdims%j_start-1, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! visc_h on v points, theta-levels
      coeff_v(i,j,k) = intw_p2v(j,1) * visc_h(i,j,k) +                  &
                       intw_p2v(j,2) * visc_h(i,j+1,k)
! visc_m on v points, theta-levels
      w_coeff_v(i,j,k) = intw_p2v(j,1) * visc_m(i,j,k) +                &
                         intw_p2v(j,2) * visc_m(i,j+1,k)
    END DO
  END DO

  DO j = pdims%j_start-1, pdims%j_end
    DO i = pdims%i_start-1, pdims%i_end
! visc_m on rho-levels at centre point (uv point on b-grid)
      coeff_centre(i,j,k) = intw_p2v(j,1) * (                           &
                                intw_p2u(i,1) * work_th(i,j,k) +        &
                                intw_p2u(i,2) * work_th(i+1,j,k) )      &
                            + intw_p2v(j,2) * (                         &
                                intw_p2u(i,1) * work_th(i,j+1,k) +      &
                                intw_p2u(i,2) * work_th(i+1,j+1,k) )
! visc_m on rho-levels
      coeff_th(i,j,k) = work_th(i,j,k)
    END DO
  END DO

END DO !  k = 1, model_levels
!$OMP END DO NOWAIT

! calculate dr/dz about theta levels
DO k = 1, flat_level - 1
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start-1, pdims%j_end
    DO i = pdims%i_start-1, pdims%i_end
      delta_z(i,j,k) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !   k = 1, flat_level - 1
DO k = flat_level, model_levels
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start-1, pdims%j_end
    DO i = pdims%i_start-1, pdims%i_end
      delta_z(i,j,k) = 1.0
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !  k = flat_level, model_levels

! Flush delta_z
!$OMP BARRIER

! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      recip_r_squared_delz(i,j,k) = 1.0 / ( h2_p_eta(i,j,k)             &
                                            * h2_p_eta(i,j,k) *         &
                                              delta_z(i,j,k) )
    END DO
  END DO
END DO !   k = 1, model_levels
!$OMP END DO

!$OMP END PARALLEL

!  CALL diffusion routines with (3D) coefficients from the subgrid
!  turbulence scheme.

! DEPENDS ON: eg_turb_diff
CALL eg_turb_diff(                                                      &
                    theta(tdims_s%i_start,tdims_s%j_start,1),fld_type_p,&
                    row_length, rows,                                   &
                    model_levels,                                       &
                    offx, offy,                                         &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    theta_star(tdims_s%i_start,tdims_s%j_start,1) )

! DEPENDS ON: eg_turb_diff
CALL eg_turb_diff(                                                      &
                    moist(tdims_s%i_start,tdims_s%j_start,1),fld_type_p,&
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    moist_star(:,:,1:) )

IF (hyd_mix_opt >= liq_only) THEN ! mix cloud liquid
  ! DEPENDS ON: eg_turb_diff
  CALL eg_turb_diff(                                                    &
                    m_cl_np1(tdims_s%i_start,tdims_s%j_start,1),        &
                    fld_type_p,                                         &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    s_m_cl(:,:,1:) )
END IF

IF (hyd_mix_opt >= liq_and_ice ) THEN ! mix cloud liq and ice
  ! DEPENDS ON: eg_turb_diff
  CALL eg_turb_diff(                                                    &
                    m_cf_np1(tdims_s%i_start,tdims_s%j_start,1),        &
                    fld_type_p,                                         &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    s_m_cf(:,:,1:) )
END IF

IF (hyd_mix_opt == all_hydro) THEN
  IF (l_mcr_qrain) THEN ! mix rain
  ! DEPENDS ON: eg_turb_diff
    CALL eg_turb_diff(                                                  &
                    m_r_np1(tdims_s%i_start,tdims_s%j_start,1),         &
                    fld_type_p,                                         &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    s_m_r(:,:,1:) )
  END IF
  IF (l_mcr_qgraup) THEN ! mix graupel
  ! DEPENDS ON: eg_turb_diff
    CALL eg_turb_diff(                                                  &
                    m_gr_np1(tdims_s%i_start,tdims_s%j_start,1),        &
                    fld_type_p,                                         &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    s_m_gr(:,:,1:) )
  END IF
  IF (l_mcr_qcf2) THEN ! mix 2nd ice
  ! DEPENDS ON: eg_turb_diff
    CALL eg_turb_diff(                                                  &
                    m_cf2_np1(tdims_s%i_start,tdims_s%j_start,1),       &
                    fld_type_p,                                         &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    zero_halo, zero_halo,                               &
                    coeff_u, coeff_v,                                   &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar,                               &
                    s_m_cf2(:,:,1:) )
  END IF
END IF

! DEPENDS ON: eg_turb_diff
CALL eg_turb_diff(                                                      &
                    w(tdims_s%i_start,tdims_s%j_start,1), fld_type_p,   &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    0, 0,                                               &
                    w_coeff_u, w_coeff_v,                               &
                    delta_z, recip_r_squared_delz ,                     &
                    swap_field_is_scalar, S_w(tdims%i_start,tdims%j_start,1) )

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                           &
!$OMP SHARED(model_levels, udims, flat_level, r_theta_uv,               &
!$OMP r_theta_levels, delta_z, recip_r_squared_delz, intw_p2u,          &
!$OMP h2_xi1_u)

! calculate dr/dz about theta levels
DO k = 1, flat_level
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start-1, udims%j_end
    DO i = udims%i_start, udims%i_end+1
      r_theta_uv(i,j,k) = intw_p2u(i,1) * r_theta_levels(i,j,k) +       &
                          intw_p2u(i,2) * r_theta_levels(i+1,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !  k = 1, flat_level - 1
DO k = 1, flat_level - 1
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start-1, udims%j_end
    DO i = udims%i_start, udims%i_end+1
      delta_z(i,j,k) = r_theta_uv(i,j,k+1) - r_theta_uv(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !   k = 1, flat_level - 1

! Flush delta_z
!$OMP BARRIER

! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      recip_r_squared_delz(i,j,k) = 1.0 / ( h2_xi1_u(i,j,k) *           &
                                              h2_xi1_u(i,j,k) *         &
                                              delta_z(i,j,k) )
    END DO
  END DO
END DO !   k = 1, model_levels
!$OMP END DO

!$OMP END PARALLEL

! DEPENDS ON: eg_turb_diff_u
CALL eg_turb_diff_u(                                                    &
                    u, fld_type_u,                                      &
                    row_length, rows,                                   &
                    model_levels,                                       &
                    offx, offy,                                         &
                    coeff_th, coeff_centre,                             &
                    delta_z, recip_r_squared_delz,                      &
                    swap_field_is_vector, S_u )

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                           &
!$OMP SHARED(model_levels, vdims, flat_level, r_theta_uv,               &
!$OMP r_theta_levels, delta_z, recip_r_squared_delz, intw_p2v,          &
!$OMP h2_xi2_v)

! calculate dr/dz about theta levels
DO k = 1, flat_level
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end+1
    DO i = vdims%i_start-1, vdims%i_end
      r_theta_uv(i,j,k) = intw_p2v(j,1) * r_theta_levels(i,j,k) +       &
                          intw_p2v(j,2) * r_theta_levels(i,j+1,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !  k = 1, flat_level - 1
DO k = 1, flat_level - 1
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end+1
    DO i = vdims%i_start-1, vdims%i_end
      delta_z(i,j,k) = r_theta_uv(i,j,k+1) - r_theta_uv(i,j,k)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO !   k = 1, flat_level - 1

! Flush delta_z
!$OMP BARRIER

! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      recip_r_squared_delz(i,j,k) = 1.0 / ( h2_xi2_v(i,j,k) *           &
                                            h2_xi2_v(i,j,k) *           &
                                              delta_z(i,j,k) )
    END DO
  END DO
END DO !   k = 1, model_levels
!$OMP END DO

!$OMP END PARALLEL


! DEPENDS ON: eg_turb_diff_v
CALL eg_turb_diff_v(                                                    &
                    v, fld_type_v,                                      &
                    row_length, n_rows,                                 &
                    model_levels,                                       &
                    offx, offy,                                         &
                    coeff_centre, coeff_th,                             &
                    delta_z, recip_r_squared_delz,                      &
                    j_start, j_stop, swap_field_is_vector, S_v )

!        Go back to negative timestep IF integrating backwards.
IF (L_Backwards) timestep = neg_timestep

DEALLOCATE (coeff_u)
DEALLOCATE (coeff_v)
DEALLOCATE (coeff_th)
DEALLOCATE (work_th)
DEALLOCATE (coeff_centre)
DEALLOCATE (w_coeff_u)
DEALLOCATE (w_coeff_v)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_diff_ctl
