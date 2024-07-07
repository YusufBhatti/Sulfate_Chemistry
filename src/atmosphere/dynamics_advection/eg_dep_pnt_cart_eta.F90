! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_dep_pnt_cart_eta_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_DEP_PNT_CART_ETA_MOD'

CONTAINS
SUBROUTINE eg_dep_pnt_cart_eta(                                       &
            g_i_pe,                                                   &
            pnt_type, dep_row_len, dep_rows, off_i, off_j, off_k,     &
            offz,depart_scheme, depart_order, l_rk_dps,               &
            depart_high_order_scheme,                                 &
            depart_monotone_scheme, first_constant_r_rho_level,       &
            l_depart_high, l_depart_mono, alpha, eta_rho_levels,      &
            etadot, u_np1, v_np1, w_np1,                              &
            etadot_np1, u, v, w, depart_xi1, depart_xi2,              &
            depart_xi3 )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE metric_terms_mod
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE eg_adjust_vert_bound_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE atm_step_local, ONLY: cycleno
USE Field_Types
USE eg_check_sl_domain_mod
USE level_heights_mod,     ONLY: eta_theta_levels
USE timestep_mod,          ONLY: timestep
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows,               &
                                 row_length, rows, n_rows, model_levels

USE UM_parvars, ONLY:  offx,offy,halo_i,halo_j, datastart,                     &
                       gc_proc_row_group,gc_proc_col_group,at_extremity,       &
                       nproc_x, nproc_y
USE um_parcore, ONLY: mype, nproc

USE model_domain_mod, ONLY: model_type, mt_lam, mt_cyclic_lam                  &
                          , mt_bi_cyclic_lam

IMPLICIT NONE
!
! Description:
!   Find u,v,w,rho departure point using a local cartesian transform method.

!  
!
! Method: ENDGame formulation version 1.02,
!         section 10.
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DEP_PNT_CART_ETA'


! Model dimensions

INTEGER, INTENT(IN) :: pnt_type, first_constant_r_rho_level

! row_length, rows and indexing offset for u or v or w type point

INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,           &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) ::                                                &
                     ! global number of points in a column
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction

LOGICAL, INTENT(IN) :: l_rk_dps

! Loop index bounds for arrays defined on p, u, v points respectively


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  depart_scheme,                                                      &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                       &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                           &
                     ! code choosing high order
                     ! interpolation scheme used in Depart routine
  depart_monotone_scheme

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_depart_high,                                                      &
                   ! True if high order interpolation scheme to
                   ! be used in Depart scheme
  l_depart_mono   ! True if monotone interpolation scheme to
                   ! be used in Depart scheme

REAL, INTENT(IN) :: alpha

REAL, INTENT(IN) ::                                                   &
  eta_rho_levels(1-offz:model_levels+offz)

REAL, INTENT(IN) ::                                                   &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
             1-offz:model_levels+offz),                               &
              ! Use vertical offsetting to allow use of extended
              ! fields for w-point interpolation
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
             1-offz:model_levels+offz),                               &
              ! Use vertical offsetting to allow use of extended
              ! fields for w-point interpolation
  w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,                  &
             0:model_levels),                                         &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

! Departure point coordinates

REAL, INTENT(OUT) ::                 &
     depart_xi1(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels),                                &
     depart_xi2(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels),                                &
     depart_xi3(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

! Local variables

INTEGER :: uv_mlevels_in, uv_mlevels_out, w_mlevels_in,               &
           w_mlevels_out, uv_first_constant_rho,                      &
           number_of_inputs, i, j, k, kk, extra_i, extra_j,           &
           error_code

INTEGER :: dep_its

LOGICAL :: l_above_top

REAL :: beta, Lx, Ly,                                                 &
        max_xi1, min_xi1, max_xi2, min_xi2
REAL :: u_a, v_a, w_a

REAL ::                              &
        int_wind_u(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                   1-off_k:model_levels),                                &
        int_wind_v(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                   1-off_k:model_levels),                                &
        int_wind_w(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                   1-off_k:model_levels)


! End of header

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)



!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

! Length and width of domain
Lx = glob_xi1_u(global_row_length) - glob_xi1_u(0)
Ly = glob_xi2_v(global_rows) - glob_xi2_v(0)

! Set Lateral boundary limits.

IF (model_type == mt_lam) THEN

  min_xi1 = glob_xi1_u(0)
  min_xi2 = glob_xi2_v(0)
  max_xi1 = glob_xi1_u(global_row_length)
  max_xi2 = glob_xi2_v(global_rows)

ELSE IF (model_type == mt_cyclic_lam) THEN
  ! check horizontal in the periodic in x LAM.

  SELECT CASE(pnt_type)
  CASE (fld_type_u)
    min_xi1 = glob_xi1_u(1-halo_i)
    max_xi1 = glob_xi1_u(global_row_length+halo_i-3)
        
    min_xi2 = glob_xi2_p(1)
    max_xi2 = glob_xi2_p(global_rows)
  CASE (fld_type_v)
    min_xi1 = glob_xi1_p(2-halo_i)
    max_xi1 = glob_xi1_p(global_row_length+halo_i-2)
        
    min_xi2 = glob_xi2_v(0)
    max_xi2 = glob_xi2_v(global_rows) 
  CASE DEFAULT       ! p-w
    min_xi1 = glob_xi1_p(2-halo_i)
    max_xi1 = glob_xi1_p(global_row_length+halo_i-2)
        
    min_xi2 = glob_xi2_p(1)
    max_xi2 = glob_xi2_p(global_rows)
  END SELECT 

ELSE IF (model_type == mt_bi_cyclic_lam) THEN

  SELECT CASE(pnt_type)
  CASE (fld_type_u)
    min_xi1 = glob_xi1_u(1-halo_i)
    max_xi1 = glob_xi1_u(global_row_length+halo_i-3)
    min_xi2 = glob_xi2_p(2-halo_j)
    max_xi2 = glob_xi2_p(global_rows+halo_j-2)
  CASE (fld_type_v)
    min_xi1 = glob_xi1_p(2-halo_i)
    max_xi1 = glob_xi1_p(global_row_length+halo_i-2)
    min_xi2 = glob_xi2_v(1-halo_j)
    max_xi2 = glob_xi2_v(global_rows+halo_j-3)
  CASE DEFAULT       ! p-w
    min_xi1 = glob_xi1_p(2-halo_i)
    min_xi2 = glob_xi2_p(2-halo_j)
    max_xi1 = glob_xi1_p(global_row_length+halo_i-2)
    max_xi2 = glob_xi2_p(global_rows+halo_j-2)
  END SELECT

END IF  !  model_type == mt_global

number_of_inputs = 1
beta = 1.0 - alpha

!----------------------------------------------------------------------
! 2.0 Compute a departure point for a u, v, w, rho point.
!----------------------------------------------------------------------

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

!----------------------------------------------------------------------
! Predictor step
!----------------------------------------------------------------------

dep_its = depart_order
IF ( l_rk_dps ) dep_its = dep_its - 1

IF ( l_rk_dps .AND. cycleno <3 ) THEN
  SELECT CASE(pnt_type)
  CASE (fld_type_u)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a,w_a)  &
!$OMP&  SHARED(model_levels,udims,u_np1,h1_xi1_u,intw_p2u,intw_v2p,v_np1, &
!$OMP& h2_xi1_u,intw_w2rho,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_u,xi2_p,timestep,eta_rho_levels)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end

          u_a  = u_np1(i,j,k)/h1_xi1_u(i,j,k)

          v_a  = (intw_p2u(i,1)*(                                     &
                                 intw_v2p(j,1)*v_np1(i,j-1,k) +       &
                                 intw_v2p(j,2)*v_np1(i,j,k)           &
                                ) +                                   &
                  intw_p2u(i,2)*(                                     &
                                 intw_v2p(j,1)*v_np1(i+1,j-1,k) +     &
                                 intw_v2p(j,2)*v_np1(i+1,j,k)         &
                                )                                     &
                 )/h2_xi1_u(i,j,k)

          w_a  = intw_p2u(i,1)*(                                      &
                                intw_w2rho(k,1)*etadot_np1(i,j,k) +   &
                                intw_w2rho(k,2)*etadot_np1(i,j,k-1)   &
                               ) +                                    &
                 intw_p2u(i,2)*(                                      &
                                intw_w2rho(k,1)*etadot_np1(i+1,j,k) + &
                                intw_w2rho(k,2)*etadot_np1(i+1,j,k-1) &
                               )

          depart_xi1(i,j,k) = xi1_u(i) - timestep*u_a
          depart_xi2(i,j,k) = xi2_p(j) - timestep*v_a
          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_v)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a,w_a)  &
!$OMP&  SHARED(model_levels,vdims,u_np1,h1_xi2_v,intw_p2v,intw_u2p,v_np1, &
!$OMP& h2_xi2_v,intw_w2rho,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_v,timestep,eta_rho_levels)
    DO k = 1, model_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end

          u_a = (intw_p2v(j,1)*(                                      &
                                intw_u2p(i,1)*u_np1(i-1,j,k) +        &
                                intw_u2p(i,2)*u_np1(i,j,k))  +        &
                 intw_p2v(j,2)*(                                      &
                                intw_u2p(i,1)*u_np1(i-1,j+1,k) +      &
                                intw_u2p(i,2)*u_np1(i,j+1,k)          &
                               )                                      &
                )/h1_xi2_v(i,j,k)

          v_a = v_np1(i,j,k)/h2_xi2_v(i,j,k)

          w_a = intw_p2v(j,1)*(                                       &
                               intw_w2rho(k,1)*etadot_np1(i,j,k) +    &
                               intw_w2rho(k,2)*etadot_np1(i,j,k-1)) + &
                intw_p2v(j,2)*(                                       &
                               intw_w2rho(k,1)*etadot_np1(i,j+1,k) +  &
                               intw_w2rho(k,2)*etadot_np1(i,j+1,k-1)  &
                               )

          depart_xi1(i,j,k) = xi1_p(i) - timestep*u_a
          depart_xi2(i,j,k) = xi2_v(j) - timestep*v_a
          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_w)

    k = 0
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        u_a = (intw_u2p(i,1)*u_np1(i-1,j,k+1) +                       &
               intw_u2p(i,2)*u_np1(i,j,k+1)                           &
              )/h1_p_eta(i,j,k)

        v_a = (intw_v2p(j,1)*v_np1(i,j-1,k+1) +                       &
               intw_v2p(j,2)*v_np1(i,j,k+1)                           &
              )/h2_p_eta(i,j,k)

        depart_xi1(i,j,k) = xi1_p(i) - timestep*u_a
        depart_xi2(i,j,k) = xi2_p(j) - timestep*v_a
        depart_xi3(i,j,k) = eta_theta_levels(k)

      END DO
    END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a)  &
!$OMP&  SHARED(model_levels,pdims,u_np1,h1_p_eta,intw_u2p,intw_rho2w,v_np1, &
!$OMP& h2_p_eta,intw_v2p,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_p,timestep,eta_theta_levels)
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          u_a = (intw_rho2w(k,1)*(                                    &
                                  intw_u2p(i,1)*u_np1(i-1,j,k+1) +    &
                                  intw_u2p(i,2)*u_np1(i,j,k+1)        &
                                 ) +                                  &
                 intw_rho2w(k,2)*                                     &
                                 (intw_u2p(i,1)*u_np1(i-1,j,k) +      &
                                  intw_u2p(i,2)*u_np1(i,j,k))         &
                )/h1_p_eta(i,j,k)

          v_a = (intw_rho2w(k,1)*(                                    &
                                  intw_v2p(j,1)*v_np1(i,j-1,k+1) +    &
                                  intw_v2p(j,2)*v_np1(i,j,k+1)        &
                                 ) +                                  &
                 intw_rho2w(k,2)*(                                    &
                                  intw_v2p(j,1)*v_np1(i,j-1,k) +      &
                                  intw_v2p(j,2)*v_np1(i,j,k))         &
                )/h2_p_eta(i,j,k)

          depart_xi1(i,j,k) = xi1_p(i) - timestep*u_a
          depart_xi2(i,j,k) = xi2_p(j) - timestep*v_a 
  
          depart_xi3(i,j,k) = eta_theta_levels(k) -                   &
                                      timestep*etadot_np1(i,j,k) 
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    k = model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        u_a = (intw_u2p(i,1)*u_np1(i-1,j,k) +                         &
               intw_u2p(i,2)*u_np1(i,j,k)                             &
              )/h1_p_eta(i,j,k)
        v_a =(intw_v2p(j,1)*v_np1(i,j-1,k) +                          &
              intw_v2p(j,2)*v_np1(i,j,k)                              &
             )/h2_p_eta(i,j,k)

        depart_xi1(i,j,k) = xi1_p(i) - timestep*u_a
        depart_xi2(i,j,k) = xi2_p(j) - timestep*v_a
        depart_xi3(i,j,k) = eta_theta_levels(k)

      END DO
    END DO

  CASE (fld_type_p)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a,w_a) &
!$OMP&  SHARED(model_levels,pdims,u_np1,h1_p,intw_u2p,intw_v2p,v_np1, &
!$OMP& h2_p,intw_w2rho,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_p,timestep,eta_rho_levels)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          u_a = (intw_u2p(i,1)*u_np1(i-1,j,k) +                       &
                 intw_u2p(i,2)*u_np1(i,j,k)                           &
                )/h1_p(i,j,k)

          v_a = (intw_v2p(j,1)*v_np1(i,j-1,k) +                       &
                 intw_v2p(j,2)*v_np1(i,j,k)                           &
                )/h2_p(i,j,k)

          w_a = intw_w2rho(k,1)*etadot_np1(i,j,k)  +                  &
                intw_w2rho(k,2)*etadot_np1(i,j,k-1)

          depart_xi1(i,j,k) = xi1_p(i) - timestep*u_a
          depart_xi2(i,j,k) = xi2_p(j) - timestep*v_a
          depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END SELECT

  ! Check the values are in range

  CALL eg_check_sl_domain( depart_xi2, depart_xi1,                  &
                    dep_row_len, dep_rows, model_levels+off_k,      &
                    max_xi1, min_xi1, max_xi2, min_xi2,             &
                    Length = Lx, Width = Ly )

  CALL eg_adjust_vert_bound(                                        &
          row_length, rows, model_levels,  g_i_pe,                  &
          pnt_type, dep_row_len,                                    &
          dep_rows, off_i, off_j, off_k, offz,                      &
          etadot, etadot_np1, depart_xi1, depart_xi2,               &
          depart_xi3 )

END IF

IF ( l_rk_dps .AND. cycleno == 1 ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!
! Iterative/corrector part
!

DO kk=1, dep_its

  CALL eg_interpolation_eta_pmf(                                      &
                     eta_rho_levels,fld_type_u, number_of_inputs,     &
                     row_length, rows, uv_mlevels_in, rows,           &
                     dep_row_len, dep_rows, uv_mlevels_out,           &
                     depart_high_order_scheme,                        &
                     depart_monotone_scheme,                          &
                     l_depart_high, l_depart_mono,                    &
                     depart_xi3, depart_xi1, depart_xi2,              &
                     mype, nproc, nproc_x, nproc_y,                   &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, error_code,                                &
                     u,int_wind_u)

  CALL eg_interpolation_eta_pmf(                                      &
                     eta_rho_levels,fld_type_v,number_of_inputs,      &
                     row_length, n_rows, uv_mlevels_in,               &
                     rows,                                            &
                     dep_row_len, dep_rows, uv_mlevels_out,           &
                     depart_high_order_scheme,                        &
                     depart_monotone_scheme,                          &
                     l_depart_high, l_depart_mono,                    &
                     depart_xi3, depart_xi1, depart_xi2,              &
                     mype, nproc, nproc_x, nproc_y,                   &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, error_code,                                &
                     v,int_wind_v)

  CALL eg_interpolation_eta_pmf(                                      &
                     eta_theta_levels,fld_type_w,                     &
                     number_of_inputs,                                &
                     row_length, rows, w_mlevels_in,                  &
                     rows,                                            &
                     dep_row_len, dep_rows, w_mlevels_out,            &
                     depart_high_order_scheme,                        &
                     depart_monotone_scheme,                          &
                     l_depart_high, l_depart_mono,                    &
                     depart_xi3, depart_xi1, depart_xi2,              &
                     mype, nproc, nproc_x, nproc_y,                   &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, error_code,                                &
                     w,int_wind_w)

  SELECT CASE(pnt_type)
  CASE (fld_type_u)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a,w_a) &
!$OMP&  SHARED(model_levels,udims,u_np1,h1_xi1_u,intw_p2u,intw_v2p,v_np1, &
!$OMP& h2_xi1_u,intw_w2rho,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_u,xi2_p,timestep,eta_rho_levels,alpha,beta, int_wind_u &
!$OMP& ,int_wind_v,int_wind_w)
    DO k = 1, model_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end

          u_a  = u_np1(i,j,k)/h1_xi1_u(i,j,k)
          v_a  = (intw_p2u(i,1)*(                                     &
                                 intw_v2p(j,1)*v_np1(i,j-1,k) +       &
                                 intw_v2p(j,2)*v_np1(i,j,k)           &
                                ) +                                   &
                  intw_p2u(i,2)*(                                     &
                                 intw_v2p(j,1)*v_np1(i+1,j-1,k) +     &
                                 intw_v2p(j,2)*v_np1(i+1,j,k)         &
                                )                                     &
                 )/h2_xi1_u(i,j,k)

          w_a  = intw_p2u(i,1)*(                                      &
                                intw_w2rho(k,1)*etadot_np1(i,j,k) +   &
                                intw_w2rho(k,2)*etadot_np1(i,j,k-1)   &
                               ) +                                    &
                 intw_p2u(i,2)*(                                      &
                                intw_w2rho(k,1)*etadot_np1(i+1,j,k) + &
                                intw_w2rho(k,2)*etadot_np1(i+1,j,k-1) &
                               )

          depart_xi1(i,j,k) = xi1_u(i) -                              &
                     timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

          depart_xi2(i,j,k) = xi2_p(j) -                              &
                     timestep*(alpha*v_a + beta*int_wind_v(i,j,k))

          depart_xi3(i,j,k) = eta_rho_levels(k)  -                     &
                     timestep*(alpha*w_a + beta*int_wind_w(i,j,k))

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_v)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a,w_a) &
!$OMP&  SHARED(model_levels,vdims,u_np1,h1_xi2_v,intw_u2p,intw_p2v,v_np1, &
!$OMP& h2_xi2_v,intw_w2rho,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_v,timestep,eta_rho_levels,alpha,beta, int_wind_u &
!$OMP& ,int_wind_v,int_wind_w)
    DO k = 1, model_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end

          u_a = (intw_p2v(j,1)*(                                           &
                                intw_u2p(i,1)*u_np1(i-1,j,k) +             &
                                intw_u2p(i,2)*u_np1(i,j,k))  +             &
                 intw_p2v(j,2)*(                                           &
                                intw_u2p(i,1)*u_np1(i-1,j+1,k) +           &
                                intw_u2p(i,2)*u_np1(i,j+1,k)               &
                               )                                           &
                )/h1_xi2_v(i,j,k)

          v_a = v_np1(i,j,k)/h2_xi2_v(i,j,k)

          w_a = intw_p2v(j,1)*(                                            &
                               intw_w2rho(k,1)*etadot_np1(i,j,k) +         &
                               intw_w2rho(k,2)*etadot_np1(i,j,k-1)) +      &
                intw_p2v(j,2)*(                                            &
                               intw_w2rho(k,1)*etadot_np1(i,j+1,k) +       &
                               intw_w2rho(k,2)*etadot_np1(i,j+1,k-1)       &
                             )

          depart_xi1(i,j,k) = xi1_p(i) -                                   &
                     timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

          depart_xi2(i,j,k) = xi2_v(j) -                                   &
                     timestep*(alpha*v_a + beta*int_wind_v(i,j,k))

          depart_xi3(i,j,k) = eta_rho_levels(k)  -                          &
                     timestep*(alpha*w_a + beta*int_wind_w(i,j,k))

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (fld_type_w)

    k = 0
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        u_a = (intw_u2p(i,1)*u_np1(i-1,j,k+1) +                     &
               intw_u2p(i,2)*u_np1(i,j,k+1)                         &
              )/h1_p_eta(i,j,k)

        v_a = (intw_v2p(j,1)*v_np1(i,j-1,k+1) +                     &
               intw_v2p(j,2)*v_np1(i,j,k+1)                         &
              )/h2_p_eta(i,j,k)

        depart_xi1(i,j,k) = xi1_p(i) -                              &
                   timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

        depart_xi2(i,j,k) = xi2_p(j) -                              &
                   timestep*(alpha*v_a + beta*int_wind_v(i,j,k))

        depart_xi3(i,j,k) = eta_theta_levels(k)
      END DO
    END DO

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a) &
!$OMP&  SHARED(model_levels,pdims,u_np1,h1_p_eta,intw_u2p,intw_rho2w,v_np1, &
!$OMP& h2_p_eta,intw_v2p,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_p,timestep,eta_theta_levels,alpha,beta, &
!$OMP& int_wind_u ,int_wind_v,int_wind_w)
    DO k = 1, model_levels-1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          u_a = (intw_rho2w(k,1)*(                                    &
                                  intw_u2p(i,1)*u_np1(i-1,j,k+1) +    &
                                  intw_u2p(i,2)*u_np1(i,j,k+1)        &
                                 ) +                                  &
                 intw_rho2w(k,2)*                                     &
                                 (intw_u2p(i,1)*u_np1(i-1,j,k) +      &
                                  intw_u2p(i,2)*u_np1(i,j,k))         &
                )/h1_p_eta(i,j,k)

          v_a = (intw_rho2w(k,1)*(                                    &
                                  intw_v2p(j,1)*v_np1(i,j-1,k+1) +    &
                                  intw_v2p(j,2)*v_np1(i,j,k+1)        &
                                 ) +                                  &
                 intw_rho2w(k,2)*(                                    &
                                  intw_v2p(j,1)*v_np1(i,j-1,k) +      &
                                  intw_v2p(j,2)*v_np1(i,j,k))         &
                )/h2_p_eta(i,j,k)

          depart_xi1(i,j,k) = xi1_p(i) -                              &
                     timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

          depart_xi2(i,j,k) = xi2_p(j) -                              &
                     timestep*(alpha*v_a + beta*int_wind_v(i,j,k))

          depart_xi3(i,j,k) = eta_theta_levels(k)  -                   &
                      timestep*(alpha*etadot_np1(i,j,k) +              &
                                     beta*int_wind_w(i,j,k))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    k = model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        u_a = (intw_u2p(i,1)*u_np1(i-1,j,k) +                       &
               intw_u2p(i,2)*u_np1(i,j,k)                           &
              )/h1_p_eta(i,j,k)
        v_a =(intw_v2p(j,1)*v_np1(i,j-1,k) +                        &
              intw_v2p(j,2)*v_np1(i,j,k)                            &
             )/h2_p_eta(i,j,k)

        depart_xi1(i,j,k) = xi1_p(i) -                              &
                 timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

        depart_xi2(i,j,k) = xi2_p(j) -                              &
                 timestep*(alpha*v_a + beta*int_wind_v(i,j,k))

        depart_xi3(i,j,k) = eta_theta_levels(k)

      END DO
    END DO

  CASE (fld_type_p)
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k,u_a,v_a) &
!$OMP&  SHARED(model_levels,pdims,u_np1,h1_p,intw_u2p,intw_w2rho,v_np1, &
!$OMP& h2_p,intw_v2p,etadot_np1,depart_xi1,depart_xi2,   &
!$OMP& depart_xi3,xi1_p,xi2_p,timestep,alpha,beta, int_wind_u &
!$OMP& ,int_wind_v,int_wind_w, eta_rho_levels)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          u_a = (intw_u2p(i,1)*u_np1(i-1,j,k) +                            &
                 intw_u2p(i,2)*u_np1(i,j,k)                                &
                )/h1_p(i,j,k)

          v_a = (intw_v2p(j,1)*v_np1(i,j-1,k) +                            &
                 intw_v2p(j,2)*v_np1(i,j,k)                                &
                )/h2_p(i,j,k)

          depart_xi1(i,j,k) = xi1_p(i) -                                   &
                     timestep*(alpha*u_a + beta*int_wind_u(i,j,k))

          depart_xi2(i,j,k) = xi2_p(j) -                                   &
                    timestep*(alpha*v_a + beta*int_wind_v(i,j,k))
           
          depart_xi3(i,j,k) = eta_rho_levels(k) -                         &
               timestep*(alpha*(intw_w2rho(k,1)*etadot_np1(i,j,k) +       &
                                intw_w2rho(k,2)*etadot_np1(i,j,k-1)) +    &
                                         beta*int_wind_w(i,j,k))

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END SELECT

  ! Check the values are in range

  CALL eg_check_sl_domain( depart_xi2, depart_xi1,                  &
                    dep_row_len, dep_rows, model_levels+off_k,      &
                    max_xi1, min_xi1, max_xi2, min_xi2,             &
                    Length = Lx, Width = Ly )

  CALL eg_adjust_vert_bound(                                        &
          row_length, rows, model_levels,  g_i_pe,                  &
          pnt_type, dep_row_len,                                    &
          dep_rows, off_i, off_j, off_k, offz,                      &
          etadot, etadot_np1, depart_xi1, depart_xi2,               &
          depart_xi3 )

END DO ! depart_order iterations

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_dep_pnt_cart_eta
END MODULE eg_dep_pnt_cart_eta_mod
