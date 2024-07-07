! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Find departure point
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics

MODULE departure_point_eta_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEPARTURE_POINT_ETA_MOD'

CONTAINS

SUBROUTINE departure_point_eta(                                        &
                 g_i_pe,                                               &
                 pnt_type, dep_row_len, dep_rows,                      &
                 off_i, off_j, off_k,offz, l_rk_dps,depart_order,      &
                 depart_high_order_scheme,                             &
                 depart_monotone_scheme, first_constant_r_rho_level,   &
                 l_depart_high,l_depart_mono, l_shallow,               &
                 alpha, eta_rho_levels,                                &
                 xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,         &
                 etadot, u_np1, v_np1, etadot_np1,                     &
                 u, v, w, depart_xi1,depart_xi2,depart_xi3)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi
USE planet_constants_mod, ONLY: planet_radius
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE eg_adjust_vert_bound_mod
USE horiz_grid_mod
USE level_heights_mod,            ONLY: eta_theta_levels

USE UM_ParParams
USE UM_parvars, ONLY:  offx,offy,halo_i,halo_j, datastart,              &
                       gc_proc_row_group,gc_proc_col_group,at_extremity,&
                       nproc_x, nproc_y
USE um_parcore, ONLY: mype, nproc
USE nlsizes_namelist_mod,  ONLY: global_row_length,                     &
                                 row_length, rows, n_rows, model_levels

USE Field_Types
USE atm_fields_bounds_mod

USE eg_check_sl_domain_mod
USE eg_parameters_mod, ONLY: l_slice
USE atm_step_local,    ONLY: cycleno
USE timestep_mod,      ONLY: timestep

USE model_domain_mod, ONLY: l_regular, model_type, mt_global

#if !defined (INTEL_FORTRAN)
USE vectlib_mod, ONLY: asin_v, atan2_v
#endif

IMPLICIT NONE
!
! Description:
!   Find u,v,w,rho departure point using a local cartesian transform method.
!
! Method:
!   ENDGame formulation version 3.02 (John Thuburns modification)
!
! Code Description:
!   Language: FORTRAN 90.
!   This code is written to UMDP3 v8.0 programming standards.
!
! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEPARTURE_POINT_ETA'

! Model dimensions

INTEGER, INTENT(IN) :: pnt_type, first_constant_r_rho_level


INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,             &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) ::                                                  &
        g_i_pe(1-halo_i:global_row_length+halo_i)
                           ! processor on my processor-row
                           ! holding a given value in i direction

! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                  &
        depart_order,                                                   &
                           ! for the chosen departure point scheme how
                           ! many iterations/terms to use.
        depart_high_order_scheme,                                       &
                           ! code choosing high order
                           ! interpolation scheme used in Depart routine
        depart_monotone_scheme
                           ! code choosing monotone
                           ! interpolation scheme used in Depart routine

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                  &
        l_depart_high,                                                  &
                         ! True if high order interpolation scheme to
                         ! be used in Depart scheme
        l_depart_mono,                                                  &
                         ! True if monotone interpolation scheme to
                         ! be used in Depart scheme
        l_shallow

LOGICAL, INTENT(IN) :: l_rk_dps ! Flag to enable RK2 scheme

REAL, INTENT(IN) :: alpha

REAL, INTENT(IN) ::                                                     &
        eta_rho_levels(1-offz:model_levels+offz)

REAL, INTENT(IN) ::                                                     &
        xi3_at_theta(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j,0:model_levels),              &
                    ! vertical coordinate at theta points
        xi3_at_rho(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,     &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at rho points. Vert offsetting
                    ! to allow use of extended height when needed
        xi3_at_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,      &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at u points
        xi3_at_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,    &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at v points
        etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
        u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,             &
                   1-offz:model_levels+offz),                           &
                    ! Use vertical offsetting to allow use of extended
                    ! fields for w-point interpolation
        v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,           &
                   1-offz:model_levels+offz),                           &
                    ! Use vertical offsetting to allow use of extended
                    ! fields for w-point interpolation
        w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,              &
                   0:model_levels),                                     &
        u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),   &
        v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels), &
        etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,             &
                   0:model_levels)

! Departure point coordinates

REAL, INTENT(OUT) ::                   &
      depart_xi1(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                1-off_k:model_levels),                                   &
      depart_xi2(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                1-off_k:model_levels),                                   &
      depart_xi3(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                1-off_k:model_levels)

! Local variables

INTEGER :: uv_mlevels_in, uv_mlevels_out, w_mlevels_in,                 &
           w_mlevels_out, uv_first_constant_rho,                        &
           number_of_inputs, i, j, k, kk, extra_i, extra_j,             &
           error_code

REAL :: beta, dxi1, dxi2,                                               &
              max_xi1, min_xi1, max_xi2, min_xi2

REAL ::                                &
        int_wind_u(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                   1-off_k:model_levels),                                  &
        int_wind_v(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                   1-off_k:model_levels),                                  &
        int_wind_w(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                   1-off_k:model_levels)

! Components of rotation matrix

REAL :: rm_11, rm_12, q_33, gam

! Arrival and departure velocities

REAL :: u_a, v_a, w_a, u_d, v_d, w_d
REAL :: snxi2_d, csxi2_d, csxi1_d, snxi1_d
REAL :: x_d, y_d, z_d, IS, rad

! Temporary local

REAL    :: temp_Csalpha
INTEGER :: jv_begin, jv_stop, dep_its

! Vectlib temporary arrays
REAL :: v1_a  (0:row_length)
REAL :: v1_b  (0:row_length)
REAL :: v1_res(0:row_length)
REAL :: v2_a  (0:row_length)
REAL :: v2_res(0:row_length)

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

dep_its = depart_order
IF ( l_rk_dps ) dep_its = 1

jv_begin = vdims%j_start
jv_stop  = vdims%j_end
IF ( .NOT. l_slice ) THEN
  IF (model_type == mt_global) THEN
    IF ( at_extremity(PSouth) ) jv_begin = jv_begin + 1
    IF ( at_extremity(PNorth) ) jv_stop  = jv_stop  - 1
  END IF
END IF

IS = 1.0
IF ( l_shallow ) IS = 0.0

IF ( l_regular ) THEN
  dxi1 = delta_xi1
  dxi2 = delta_xi2
ELSE  ! variable resolution
  dxi1 = xi1_u(udims%i_end)-xi1_u(udims%i_end-1)
  dxi2 = xi2_v(vdims%i_end)-xi2_v(vdims%i_end-1)
END IF ! l_regular

! Set Lateral boundary limits.

min_xi1 = (2-halo_i)*dxi1
max_xi1 = 2.0*pi + (halo_i-2)*dxi1

min_xi2 = -pi*0.5  
max_xi2 = pi*0.5 

number_of_inputs = 1
beta = 1.0 - alpha

!-----------------------------------------------------------------------
! 2.0 Compute a departure point for a u, v, w, rho point.
!-----------------------------------------------------------------------


IF ( l_rk_dps .AND. cycleno < 3 ) THEN
  SELECT CASE(pnt_type)
  CASE (fld_type_u)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, udims, u_np1, intw_p2u, intw_v2p,          &
!$OMP         intw_w2rho, v_np1, etadot_np1, IS, xi3_at_u, timestep,   &
!$OMP         xi1_u, csxi2_p, snxi2_p, planet_radius, eta_rho_levels,  &
!$OMP         depart_xi1, depart_xi2, depart_xi3 )                     &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, rad, x_d, y_d, z_d,             &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end

          ! Arrival points

          u_a    = u_np1(i,j,k)
          v_a    = intw_p2u(i,1)*(intw_v2p(j,1)*v_np1(i,j-1,k)         &
                                 +intw_v2p(j,2)*v_np1(i,j,k) )         &
                  +intw_p2u(i,2)*(intw_v2p(j,1)*v_np1(i+1,j-1,k)       &
                                 +intw_v2p(j,2)*v_np1(i+1,j,k) )
          w_a    = intw_p2u(i,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                 +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                  +intw_p2u(i,2)*(intw_w2rho(k,1)*etadot_np1(i+1,j,k)  &
                                 +intw_w2rho(k,2)*etadot_np1(i+1,j,k-1))

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_u(i,j,k) - planet_radius))
          u_a    = u_a*rad
          v_a    = v_a*rad

          ! Update departure points

          x_d = -timestep*u_a
          y_d = -timestep*v_a
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO
        
        CALL atan2_v(udims%i_len, v1_a  (udims%i_start:udims%i_end), &
                                  v1_b  (udims%i_start:udims%i_end), &
                                  v1_res(udims%i_start:udims%i_end))

        CALL asin_v(udims%i_len, v2_a  (udims%i_start:udims%i_end), &
                                 v2_res(udims%i_start:udims%i_end))

        DO i = udims%i_start, udims%i_end

          depart_xi1(i,j,k) = xi1_u(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)
        
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_v)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, jv_begin, jv_stop, vdims, intw_p2v,        &
!$OMP         intw_u2p, u_np1, intw_w2rho, etadot_np1, Is, xi3_at_v,   &
!$OMP         timestep, v_np1, planet_radius, xi1_p, csxi2_v, snxi2_v, &
!$OMP         depart_xi1, depart_xi2, depart_xi3, eta_rho_levels )     &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, rad, x_d, y_d, z_d,             &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = jv_begin, jv_stop
        DO i = vdims%i_start, vdims%i_end

          ! Arrival points

          u_a    = intw_p2v(j,1)*(intw_u2p(i,1)*u_np1(i-1,j,k)         &
                                 +intw_u2p(i,2)*u_np1(i,j,k))          &
                  +intw_p2v(j,2)*(intw_u2p(i,1)*u_np1(i-1,j+1,k)       &
                                 +intw_u2p(i,2)*u_np1(i,j+1,k))
          v_a    = v_np1(i,j,k)
          w_a    = intw_p2v(j,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                 +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                  +intw_p2v(j,2)*(intw_w2rho(k,1)*etadot_np1(i,j+1,k)  &
                                 +intw_w2rho(k,2)*etadot_np1(i,j+1,k-1))

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_v(i,j,k) - planet_radius))

          u_a = u_a*rad
          v_a = v_a*rad

          ! Update departure points

          x_d = -timestep*u_a
          y_d = -timestep*v_a
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_v(j)-y_d*snxi2_v(j)
          v2_a(i) = y_d*csxi2_v(j)+z_d*snxi2_v(j)

        END DO

        CALL atan2_v(vdims%i_len, v1_a  (vdims%i_start:vdims%i_end), &
                                  v1_b  (vdims%i_start:vdims%i_end), &
                                  v1_res(vdims%i_start:vdims%i_end))

        CALL asin_v(vdims%i_len, v2_a  (vdims%i_start:vdims%i_end), &
                                 v2_res(vdims%i_start:vdims%i_end))

        DO i = vdims%i_start, vdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_w)

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( model_levels, pdims, intw_u2p, u_np1, v_np1, intw_v2p,   &
!$OMP         Is, xi3_at_theta, timestep, csxi2_p, snxi2_p,intw_rho2w, &
!$OMP         etadot_np1, depart_xi3, xi1_p, depart_xi1, depart_xi2,   &
!$OMP         xi2_p, eta_theta_levels, planet_radius )                 &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, rad, x_d, y_d, z_d,             &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          ! Arrival points

          u_a = intw_rho2w(k,1)*(intw_u2p(i,1)*u_np1(i-1,j,k+1)        &
                                +intw_u2p(i,2)*u_np1(i,j,k+1))         &
               +intw_rho2w(k,2)*(intw_u2p(i,1)*u_np1(i-1,j,k)          &
                                +intw_u2p(i,2)*u_np1(i,j,k))
          v_a = intw_rho2w(k,1)*(intw_v2p(j,1)*v_np1(i,j-1,k+1)        &
                                +intw_v2p(j,2)*v_np1(i,j,k+1))         &
               +intw_rho2w(k,2)*(intw_v2p(j,1)*v_np1(i,j-1,k)          &
                                +intw_v2p(j,2)*v_np1(i,j,k))
          w_a = etadot_np1(i,j,k)

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_theta(i,j,k) - planet_radius))

          u_a = u_a*rad
          v_a = v_a*rad

          ! Update departure points

          x_d = -timestep*u_a
          y_d = -timestep*v_a
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_theta_levels(k) - timestep*w_a

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO

        CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                  v1_b  (pdims%i_start:pdims%i_end), &
                                  v1_res(pdims%i_start:pdims%i_end))

        CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                                 v2_res(pdims%i_start:pdims%i_end))

        DO i = pdims%i_start, pdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    k = 0
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        ! Arrival points

        u_a    = intw_u2p(i,1)*u_np1(i-1,j,k+1) +                    &
                 intw_u2p(i,2)*u_np1(i,j,k+1)
        v_a    = intw_v2p(j,1)*v_np1(i,j-1,k+1) +                    &
                 intw_v2p(j,2)*v_np1(i,j,k+1)

        rad    = 1.0/(planet_radius                                  &
                         + IS*(xi3_at_theta(i,j,k) - planet_radius))

        u_a    = u_a*rad
        v_a    = v_a*rad

        ! Update departure points

        x_d = -timestep*u_a
        y_d = -timestep*v_a
        z_d = SQRT(1.0 - x_d**2 - y_d**2)

        depart_xi3(i,j,k) = eta_theta_levels(k)

        v1_a(i) = x_d
        v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
        v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

      END DO

      CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                v1_b  (pdims%i_start:pdims%i_end), &
                                v1_res(pdims%i_start:pdims%i_end))

      CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                               v2_res(pdims%i_start:pdims%i_end))

      DO i = pdims%i_start, pdims%i_end
        
        depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
        depart_xi2(i,j,k) = v2_res(i)
        
      END DO
    END DO
!$OMP END DO NOWAIT

    k = model_levels
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        ! Arrival points

        u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                    &
                 intw_u2p(i,2)*u_np1(i,j,k)
        v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                    &
                 intw_v2p(j,2)*v_np1(i,j,k)

        rad    = 1.0/(planet_radius                                &
                       + IS*(xi3_at_theta(i,j,k) - planet_radius))

        u_a = u_a*rad
        v_a = v_a*rad

        ! Update departure points

        x_d = -timestep*u_a
        y_d = -timestep*v_a
        z_d = SQRT(1.0 - x_d**2 - y_d**2)

        depart_xi3(i,j,k) = eta_theta_levels(k)

        v1_a(i) = x_d
        v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
        v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

      END DO

      CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                v1_b  (pdims%i_start:pdims%i_end), &
                                v1_res(pdims%i_start:pdims%i_end))

      CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                               v2_res(pdims%i_start:pdims%i_end))

      DO i = pdims%i_start, pdims%i_end
        
        depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
        depart_xi2(i,j,k) = v2_res(i)
        
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  CASE (fld_type_p)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, pdims, intw_u2p, u_np1, intw_v2p, v_np1,   &
!$OMP         intw_w2rho, etadot_np1, Is, xi3_at_rho, timestep, xi1_p, &
!$OMP         csxi2_p, planet_radius, snxi2_p, eta_rho_levels,         &
!$OMP         depart_xi1, depart_xi2, depart_xi3 )                     &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, rad, x_d, y_d, z_d,             &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          ! Arrival points

          u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                      &
                   intw_u2p(i,2)*u_np1(i,j,k)
          v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                      &
                   intw_v2p(j,2)*v_np1(i,j,k)
          w_a    = intw_w2rho(k,1)*etadot_np1(i,j,k) +                 &
                   intw_w2rho(k,2)*etadot_np1(i,j,k-1)

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_rho(i,j,k) - planet_radius))

          u_a    = u_a*rad
          v_a    = v_a*rad

          ! Update departure points

          x_d = -timestep*u_a
          y_d = -timestep*v_a
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO
        
        CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                  v1_b  (pdims%i_start:pdims%i_end), &
                                  v1_res(pdims%i_start:pdims%i_end))

        CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                                 v2_res(pdims%i_start:pdims%i_end))

        DO i = pdims%i_start, pdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)        

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! Check the values are in range

  END SELECT

  CALL eg_check_sl_domain( depart_xi2, depart_xi1,                      &
                     dep_row_len, dep_rows, model_levels+off_k,         &
                     max_xi1, min_xi1, max_xi2, min_xi2)

  CALL eg_adjust_vert_bound(                                            &
           row_length, rows, model_levels, g_i_pe,                      &
           pnt_type, dep_row_len,                                       &
           dep_rows, off_i, off_j, off_k, offz,                         &
           etadot, etadot_np1, depart_xi1, depart_xi2,                  &
           depart_xi3)

END IF

IF ( l_rk_dps .AND. cycleno == 1 ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF


SELECT CASE (pnt_type)
CASE (fld_type_u)

  ! Settings for levels used in Interpolation subroutine
  uv_mlevels_in  = model_levels
  uv_mlevels_out = model_levels
  w_mlevels_in   = model_levels + 1
  w_mlevels_out  = model_levels
  uv_first_constant_rho = first_constant_r_rho_level

CASE (fld_type_v)

  ! Settings for levels used in Interpolation subroutine
  uv_mlevels_in  = model_levels
  uv_mlevels_out = model_levels
  w_mlevels_in   = model_levels  + 1
  w_mlevels_out  = model_levels
  uv_first_constant_rho = first_constant_r_rho_level

CASE (fld_type_w)

  ! Settings for levels used in Interpolation subroutine
  uv_mlevels_in  = model_levels  + 2
  uv_mlevels_out = model_levels  + 1
  w_mlevels_in   = model_levels  + 1
  w_mlevels_out  = model_levels  + 1
  uv_first_constant_rho = first_constant_r_rho_level + 2

CASE (fld_type_p)

  ! Settings for levels used in Interpolation subroutine
  uv_mlevels_in  = model_levels
  uv_mlevels_out = model_levels
  w_mlevels_in   = model_levels + 1
  w_mlevels_out  = model_levels
  uv_first_constant_rho = first_constant_r_rho_level

END SELECT


!-----------------------------------------------------------------------
!
! 2.2: Iterate algorithm
!
!-----------------------------------------------------------------------

DO kk=1, dep_its

  CALL eg_interpolation_eta_pmf(                                        &
                          eta_rho_levels, fld_type_u, number_of_inputs, &
                          row_length, rows, uv_mlevels_in, rows,        &
                          dep_row_len, dep_rows, uv_mlevels_out,        &
                          depart_high_order_scheme,                     &
                          depart_monotone_scheme,                       &
                          L_depart_high, L_depart_mono,                 &
                          depart_xi3, depart_xi1, depart_xi2,           &
                          mype, nproc, nproc_x, nproc_y,                &
                          halo_i, halo_j,                               &
                          global_row_length, datastart, at_extremity,   &
                          g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                          0, 0,error_code,                              &
                          u,int_wind_u)

  CALL eg_interpolation_eta_pmf(                                        &
                          eta_rho_levels,fld_type_v, number_of_inputs,  &
                          row_length, n_rows, uv_mlevels_in,            &
                          rows,                                         &
                          dep_row_len, dep_rows, uv_mlevels_out,        &
                          depart_high_order_scheme,                     &
                          depart_monotone_scheme,                       &
                          L_depart_high, L_depart_mono,                 &
                          depart_xi3, depart_xi1, depart_xi2,           &
                          mype, nproc, nproc_x, nproc_y,                &
                          halo_i, halo_j,                               &
                          global_row_length, datastart, at_extremity,   &
                          g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                          0,0,error_code,                               &
                          v,int_wind_v)

  CALL eg_interpolation_eta_pmf(                                        &
                          eta_theta_levels,fld_type_w,                  &
                          number_of_inputs,                             &
                          row_length, rows, w_mlevels_in,               &
                          rows,                                         &
                          dep_row_len, dep_rows, w_mlevels_out,         &
                          depart_high_order_scheme,                     &
                          depart_monotone_scheme,                       &
                          L_depart_high, L_depart_mono,                 &
                          depart_xi3, depart_xi1, depart_xi2,           &
                          mype, nproc, nproc_x, nproc_y,                &
                          halo_i, halo_j,                               &
                          global_row_length, datastart, at_extremity,   &
                          g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                          0, 0,error_code,                              &
                          w,int_wind_w)

  SELECT CASE(pnt_type)
  CASE (fld_type_u)

    ! Compute u-departure point using current iteration estimates

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, udims, depart_xi1, xi1_u, depart_xi2,      &
!$OMP         snxi2_p, csxi2_p, timestep, u_np1, intw_p2u, intw_v2p,   &
!$OMP         v_np1, intw_w2rho, etadot_np1, Is, xi3_at_u, int_wind_u, &
!$OMP         int_wind_v, int_wind_w, depart_xi3, eta_rho_levels,      &
!$OMP         alpha, beta, planet_radius )                             &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, u_d, v_d, w_d, rad, x_d,        &
!$OMP          y_d, z_d, rm_11, rm_12,  gam, q_33, temp_Csalpha,       &
!$OMP          snxi1_d, csxi1_d, snxi2_d, csxi2_d,                     &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end

          snxi1_d = SIN(depart_xi1(i,j,k)-xi1_u(i))
          csxi1_d = COS(depart_xi1(i,j,k)-xi1_u(i))
          snxi2_d = SIN(depart_xi2(i,j,k))
          csxi2_d = COS(depart_xi2(i,j,k))

          q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
          gam          = (1.0 + q_33)
          temp_Csalpha = 1.0/gam
          gam          = 0.5*gam*timestep

          ! Components of rotation matrix

          rm_11 = (csxi2_d*csxi2_p(j)                                  &
                       +(1.0+snxi2_p(j)*snxi2_d)*csxi1_d)              &
                 *temp_Csalpha
          rm_12 =-Snxi1_d*(Snxi2_p(j) + Snxi2_d)*temp_Csalpha

          ! Arrival points

          u_a    = u_np1(i,j,k)
          v_a    = intw_p2u(i,1)*(intw_v2p(j,1)*v_np1(i,j-1,k)         &
                                 +intw_v2p(j,2)*v_np1(i,j,k) )         &
                  +intw_p2u(i,2)*(intw_v2p(j,1)*v_np1(i+1,j-1,k)       &
                                 +intw_v2p(j,2)*v_np1(i+1,j,k) )
          w_a    = intw_p2u(i,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                 +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                  +intw_p2u(i,2)*(intw_w2rho(k,1)*etadot_np1(i+1,j,k)  &
                                 +intw_w2rho(k,2)*etadot_np1(i+1,j,k-1))

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_u(i,j,k) - planet_radius))
          u_a    = u_a *rad
          v_a    = v_a *rad

          ! Rotated velocities at departure points

          u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
          v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
          w_d = int_wind_w(i,j,k)

          ! Update departure points

          x_d = -gam*(alpha*u_a + beta*u_d)
          y_d = -gam*(alpha*v_a + beta*v_d)
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                             -timestep*(alpha*w_a + beta*w_d)

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO
        
        CALL atan2_v(udims%i_len, v1_a  (udims%i_start:udims%i_end), &
                                  v1_b  (udims%i_start:udims%i_end), &
                                  v1_res(udims%i_start:udims%i_end))

        CALL asin_v(udims%i_len, v2_a  (udims%i_start:udims%i_end), &
                                 v2_res(udims%i_start:udims%i_end))

        DO i = udims%i_start, udims%i_end

          depart_xi1(i,j,k) = xi1_u(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)
        
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_v)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, vdims, depart_xi1, depart_xi2, timestep,   &
!$OMP         u_np1, intw_w2rho, etadot_np1, v_np1, Is, int_wind_u,    &
!$OMP         int_wind_v, int_wind_w, xi3_at_v, alpha, beta, intw_u2p, &
!$OMP         eta_rho_levels, depart_xi3, jv_begin, jv_stop, csxi2_v,  &
!$OMP         intw_p2v, xi1_p, snxi2_v, planet_radius )                &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, u_d, v_d, w_d, rad, x_d, y_d,   &
!$OMP          z_d, rm_11, rm_12, gam, q_33, temp_Csalpha, snxi1_d,    &
!$OMP          csxi1_d, snxi2_d, csxi2_d,                              &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = jv_begin, jv_stop
        DO i = vdims%i_start, vdims%i_end

          snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
          csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
          snxi2_d = SIN(depart_xi2(i,j,k))
          csxi2_d = COS(depart_xi2(i,j,k))

          q_33         = snxi2_d*snxi2_v(j)+csxi2_d*csxi2_v(j)*csxi1_d
          gam          = (1.0 + q_33)
          temp_Csalpha = 1.0/gam
          gam          = 0.5*gam*timestep

          ! Components of rotation matrix
          rm_11 = (csxi2_d*csxi2_v(j)                                  &
                    +(1.0+snxi2_v(j)*snxi2_d)*csxi1_d)*temp_Csalpha
          rm_12 =-snxi1_d*(snxi2_v(j) + snxi2_d)*temp_Csalpha

          ! Arrival points

          u_a    = intw_p2v(j,1)*(intw_u2p(i,1)*u_np1(i-1,j,k)         &
                                 +intw_u2p(i,2)*u_np1(i,j,k))          &
                  +intw_p2v(j,2)*(intw_u2p(i,1)*u_np1(i-1,j+1,k)       &
                                 +intw_u2p(i,2)*u_np1(i,j+1,k))
          v_a    = v_np1(i,j,k)
          w_a    = intw_p2v(j,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                 +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                  +intw_p2v(j,2)*(intw_w2rho(k,1)*etadot_np1(i,j+1,k)  &
                                 +intw_w2rho(k,2)*etadot_np1(i,j+1,k-1))

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_v(i,j,k) - planet_radius))

          u_a = u_a *rad
          v_a = v_a *rad

          ! Rotated velocities at departure points

          u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
          v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
          w_d = int_wind_w(i,j,k)

          ! Update departure points

          x_d = -gam*(alpha*u_a + beta*u_d)
          y_d = -gam*(alpha*v_a + beta*v_d)
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                             - timestep*(alpha*w_a + beta*w_d)

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_v(j)-y_d*snxi2_v(j)
          v2_a(i) = y_d*csxi2_v(j)+z_d*snxi2_v(j)

        END DO

        CALL atan2_v(vdims%i_len, v1_a  (vdims%i_start:vdims%i_end), &
                                  v1_b  (vdims%i_start:vdims%i_end), &
                                  v1_res(vdims%i_start:vdims%i_end))

        CALL asin_v(vdims%i_len, v2_a  (vdims%i_start:vdims%i_end), &
                                 v2_res(vdims%i_start:vdims%i_end))

        DO i = vdims%i_start, vdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_w)

    ! Compute w-departure point using current iteration estimates
    ! Use no vertical shear assumption for horizontal momentum

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( model_levels, pdims, depart_xi1, depart_xi2, snxi2_p,    &
!$OMP         csxi2_p, timestep, u_np1, intw_v2p, v_np1, etadot_np1,   &
!$OMP         Is, int_wind_u, int_wind_v, int_wind_w, planet_radius,   &
!$OMP         beta, alpha,xi2_p, eta_theta_levels, depart_xi3,         &
!$OMP         intw_rho2w, xi1_p, xi3_at_theta ,intw_u2p )              &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, u_d, v_d, w_d, rad, x_d, y_d,   &
!$OMP          z_d, rm_11, rm_12, gam, q_33, temp_Csalpha, snxi1_d,    &
!$OMP          csxi1_d, snxi2_d, csxi2_d,                              &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
          csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
          snxi2_d = SIN(depart_xi2(i,j,k))
          csxi2_d = COS(depart_xi2(i,j,k))

          q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
          gam          = (1.0 + q_33)
          temp_Csalpha = 1.0/gam
          gam          = 0.5*gam*timestep

          ! Components of rotation matrix
          rm_11 = (csxi2_d*csxi2_p(j) +                                &
                       (1+snxi2_p(j)*snxi2_d)*csxi1_d)*temp_Csalpha
          rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

          ! Arrival points

          u_a = intw_rho2w(k,1)*(intw_u2p(i,1)*u_np1(i-1,j,k+1)        &
                                +intw_u2p(i,2)*u_np1(i,j,k+1))         &
               +intw_rho2w(k,2)*(intw_u2p(i,1)*u_np1(i-1,j,k)          &
                                +intw_u2p(i,2)*u_np1(i,j,k))
          v_a = intw_rho2w(k,1)*(intw_v2p(j,1)*v_np1(i,j-1,k+1)        &
                                +intw_v2p(j,2)*v_np1(i,j,k+1))         &
               +intw_rho2w(k,2)*(intw_v2p(j,1)*v_np1(i,j-1,k)          &
                                +intw_v2p(j,2)*v_np1(i,j,k))
          w_a = etadot_np1(i,j,k)

          rad    = 1.0/(planet_radius                                  &
                           + IS*(xi3_at_theta(i,j,k) - planet_radius))

          u_a = u_a *rad
          v_a = v_a *rad

          ! Rotated velocities at departure points

          u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
          v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
          w_d = int_wind_w(i,j,k)

          ! Update departure points

          x_d = -gam*(alpha*u_a + beta*u_d)
          y_d = -gam*(alpha*v_a + beta*v_d)
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_theta_levels(k)                      &
                              -timestep*(alpha*w_a + beta*w_d)

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO

        CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                  v1_b  (pdims%i_start:pdims%i_end), &
                                  v1_res(pdims%i_start:pdims%i_end))

        CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                                 v2_res(pdims%i_start:pdims%i_end))

        DO i = pdims%i_start, pdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)

        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    k = 0
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
        csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
        snxi2_d = SIN(depart_xi2(i,j,k))
        csxi2_d = COS(depart_xi2(i,j,k))

        q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
        gam          = (1.0 + q_33)
        temp_Csalpha = 1.0/gam
        gam          = 0.5*gam*timestep

        ! Components of rotation matrix
        rm_11 = (csxi2_d*csxi2_p(j) + (1+snxi2_p(j)*snxi2_d)*csxi1_d)   &
                *temp_Csalpha
        rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

        ! Arrival points

        u_a    = intw_u2p(i,1)*u_np1(i-1,j,k+1) +                       &
                 intw_u2p(i,2)*u_np1(i,j,k+1)
        v_a    = intw_v2p(j,1)*v_np1(i,j-1,k+1) +                       &
                 intw_v2p(j,2)*v_np1(i,j,k+1)

        rad    = 1.0/(planet_radius                                     &
                            + IS*(xi3_at_theta(i,j,k) - planet_radius))

        u_a    = u_a *rad
        v_a    = v_a *rad

        ! Rotated velocities at departure points

        u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
        v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)

        ! Update departure points

        x_d = -gam*(alpha*u_a + beta*u_d)
        y_d = -gam*(alpha*v_a + beta*v_d)
        z_d = SQRT(1.0 - x_d**2 - y_d**2)

        depart_xi3(i,j,k) = eta_theta_levels(k)

        v1_a(i) = x_d
        v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
        v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

      END DO

      CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                v1_b  (pdims%i_start:pdims%i_end), &
                                v1_res(pdims%i_start:pdims%i_end))

      CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                               v2_res(pdims%i_start:pdims%i_end))

      DO i = pdims%i_start, pdims%i_end
        
        depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
        depart_xi2(i,j,k) = v2_res(i)
        
      END DO
    END DO
!$OMP END DO NOWAIT

    ! no vertical momentum shear above top rho level
    k = model_levels
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
        csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
        snxi2_d = SIN(depart_xi2(i,j,k))
        csxi2_d = COS(depart_xi2(i,j,k))

        q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
        gam          = (1.0 + q_33)
        temp_Csalpha = 1.0/gam
        gam          = 0.5*gam*timestep

        ! Components of rotation matrix
        rm_11 = (csxi2_d*csxi2_p(j)                                     &
                        + (1+snxi2_p(j)*snxi2_d)*csxi1_d)*temp_Csalpha
        rm_12 =-Snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

        ! Arrival points

        u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                         &
                 intw_u2p(i,2)*u_np1(i,j,k)
        v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                         &
                 intw_v2p(j,2)*v_np1(i,j,k)

        rad    = 1.0/(planet_radius                                     &
                            + IS*(xi3_at_theta(i,j,k) - planet_radius))

        u_a = u_a *rad
        v_a = v_a *rad

        ! Rotated velocities at departure points

        u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
        v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)

        ! Update departure points

        x_d = -gam*(alpha*u_a + beta*u_d)
        y_d = -gam*(alpha*v_a + beta*v_d)
        z_d = SQRT(1.0 - x_d**2 - y_d**2)

        depart_xi3(i,j,k) = eta_theta_levels(k)

        v1_a(i) = x_d
        v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
        v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

      END DO

      CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                v1_b  (pdims%i_start:pdims%i_end), &
                                v1_res(pdims%i_start:pdims%i_end))

      CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                               v2_res(pdims%i_start:pdims%i_end))

      DO i = pdims%i_start, pdims%i_end
        
        depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
        depart_xi2(i,j,k) = v2_res(i)
        
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  CASE (fld_type_p)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, pdims, depart_xi1, depart_xi2, snxi2_p,    &
!$OMP         csxi2_p, timestep, u_np1, intw_v2p, v_np1, intw_w2rho,   &
!$OMP         etadot_np1, Is, int_wind_u, int_wind_v, int_wind_w,      &
!$OMP         alpha, beta, eta_rho_levels, depart_xi3, xi3_at_rho,     &
!$OMP         xi1_p, intw_u2p, planet_radius)                          &
!$OMP PRIVATE( i, j, k, u_a, v_a, w_a, u_d, v_d, w_d, rad, x_d, y_d,   &
!$OMP          z_d, rm_11, rm_12, gam, q_33, temp_Csalpha, snxi1_d,    &
!$OMP          csxi1_d, snxi2_d, csxi2_d,                              &
!$OMP          v1_a, v1_b, v1_res, v2_a, v2_res)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
          csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
          snxi2_d = SIN(depart_xi2(i,j,k))
          csxi2_d = COS(depart_xi2(i,j,k))

          q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
          gam          = (1.0 + q_33)
          temp_Csalpha = 1.0/gam
          gam          = 0.5*gam*timestep

          ! Components of rotation matrix
          rm_11 = (csxi2_d*csxi2_p(j) + (1.0+snxi2_p(j)*snxi2_d)       &
                  *csxi1_d)*temp_Csalpha
          rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

          ! Arrival points

          u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                      &
                   intw_u2p(i,2)*u_np1(i,j,k)
          v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                      &
                   intw_v2p(j,2)*v_np1(i,j,k)
          w_a    = intw_w2rho(k,1)*etadot_np1(i,j,k) +                 &
                   intw_w2rho(k,2)*etadot_np1(i,j,k-1)

          rad    = 1.0/(planet_radius                                   &
                            + IS*(xi3_at_rho(i,j,k) - planet_radius))

          u_a    = u_a *rad
          v_a    = v_a *rad

          ! Rotated velocities at departure points

          u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
          v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
          w_d = int_wind_w(i,j,k)

          ! Update departure points

          x_d = -gam*(alpha*u_a + beta*u_d)
          y_d = -gam*(alpha*v_a + beta*v_d)
          z_d = SQRT(1.0 - x_d**2 - y_d**2)

          depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                             -timestep*(alpha*w_a + beta*w_d)

          v1_a(i) = x_d
          v1_b(i) = z_d*csxi2_p(j)-y_d*snxi2_p(j)
          v2_a(i) = y_d*csxi2_p(j)+z_d*snxi2_p(j)

        END DO
        
        CALL atan2_v(pdims%i_len, v1_a  (pdims%i_start:pdims%i_end), &
                                  v1_b  (pdims%i_start:pdims%i_end), &
                                  v1_res(pdims%i_start:pdims%i_end))

        CALL asin_v(pdims%i_len, v2_a  (pdims%i_start:pdims%i_end), &
                                 v2_res(pdims%i_start:pdims%i_end))

        DO i = pdims%i_start, pdims%i_end

          depart_xi1(i,j,k) = xi1_p(i) + v1_res(i)
          depart_xi2(i,j,k) = v2_res(i)        

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END SELECT


  CALL eg_check_sl_domain( depart_xi2, depart_xi1,                         &
                        dep_row_len, dep_rows, model_levels+off_k,         &
                        max_xi1, min_xi1, max_xi2, min_xi2)

  CALL eg_adjust_vert_bound(                                               &
              row_length, rows, model_levels, g_i_pe,                      &
              pnt_type, dep_row_len,                                       &
              dep_rows, off_i, off_j, off_k, offz,                         &
              etadot, etadot_np1, depart_xi1, depart_xi2,                  &
              depart_xi3)

END DO ! dep_its iterations

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE departure_point_eta

#if defined (INTEL_FORTRAN)

! The following two routines are defined directly into the module
! instead of taking them from vectlib_mod for performance reason.
! Intel compiler does not inline these two routines when they are in
! the vectlib_mod.

! This routine are too small to have DRHOOK calls.

SUBROUTINE atan2_v(n,a,b,y)
IMPLICIT NONE

INTEGER :: n
REAL :: y(n), a(n),b(n)
INTEGER :: i

DO i=1, n
  y(i) = ATAN2(a(i),b(i))
END DO

RETURN
END SUBROUTINE atan2_v

SUBROUTINE asin_v(n,x,y)
IMPLICIT NONE

INTEGER :: n
REAL :: y(n), x(n)
INTEGER :: i

DO i=1, n
  y(i) = ASIN(x(i))
END DO

RETURN
END SUBROUTINE asin_v

#endif

END MODULE departure_point_eta_mod
