! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE eg_turb_smagorinsky(                                         &
                               u, v, w,                                 &
                               model_levels,                            &
                               r_theta_levels, r_rho_levels)

USE horiz_grid_mod, ONLY: xi1_p, xi1_u, xi2_p, xi2_v,                   &
     csxi2_p, csxi2_v,intw_u2p, intw_v2p, intw_p2u, intw_p2v
USE metric_terms_mod, ONLY: h1_p_eta, h2_p_eta, h1_p, h2_p
USE turb_diff_mod, ONLY: diff_factor, mix_factor
USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, max_diff, shear,           &
                             delta_smag
                             
USE planet_constants_mod, ONLY: vkman
USE bl_option_mod, ONLY: blending_option, off
USE timestep_mod, ONLY: timestep, recip_timestep
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, pdims_l,               &
     udims, udims_s, vdims, vdims_s, tdims_l
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Description: Calculates coefficients for use in subgrid turbulence
!              scheme
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
IMPLICIT NONE

! Variables with Intent (In)
INTEGER, INTENT(IN) ::                                                  &
  model_levels
                         ! number of model levels.

REAL, INTENT(IN) ::                                                     &
  u(udims_s%i_start:udims_s%i_end,                                      &
    udims_s%j_start:udims_s%j_end, model_levels)                        &
, v(vdims_s%i_start:vdims_s%i_end,                                      &
    vdims_s%j_start:vdims_s%j_end, model_levels)                        &
, w(pdims_s%i_start:pdims_s%i_end,                                      &
    pdims_s%j_start:pdims_s%j_end, 0:model_levels)                      &
, r_theta_levels (tdims_l%i_start:tdims_l%i_end,                        &
                  tdims_l%j_start:tdims_l%j_end, 0:model_levels)        &
, r_rho_levels (pdims_l%i_start:pdims_l%i_end,                          &
                pdims_l%j_start:pdims_l%j_end, model_levels)

!Local parameters
INTEGER :: i, j, k              ! Loop counters

REAL ::                                                                 &
  delta_x                                                               &
        ! Sea level longitude spacing in m
, delta_y                                                               &
        ! Sea level latitude spacing in m
, delta_z(pdims%i_start:pdims%i_end,                                    &
          pdims%j_start:pdims%j_end, model_levels)                      &
        ! Backward difference dz at rho points
, delta_zn(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end, model_levels)                 &
        ! Backward difference dz at theta points
        ! Note: k=1 is special - uses surface
        ! Note: delta_zn(k) is correct for theta(k-1)
, delta_zn_u                                                            &
        ! Average dz at theta points to u points (theta levels)
, delta_zn_v                                                            &
        ! Average dz at theta points to v points (theta levels)
, rdz(pdims%i_start:pdims%i_end,                                        &
      pdims%j_start:pdims%j_end, model_levels)                          &
        ! 1/delta_z
, rdzn_u(udims%i_start:udims%i_end+1,                                   &
         udims%j_start:udims%j_end, model_levels)                       &
        ! 1/delta_zn_u
, rdzn_v(vdims%i_start:vdims%i_end,                                     &
         vdims%j_start:vdims%j_end+1, model_levels)                     &
        ! 1/delta_zn_v
, z(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, model_levels)

REAL ::                                                                 &
  dx_rho                                                                &
!               ! dx on rho points
, dx_theta_u                                                            &
!               ! dx above u points on theta levels
, dx_rho_centre                                                         &
!               ! dx in centre of grid box on rho levels
, dy_rho                                                                &
!               ! dy on rho points
, dy_theta_v                                                            &
!               ! dy above v points on theta levels
, dy_rho_centre                                                         &
!               ! dy in centre of grid box on rho levels
, cx_rho                                                                &
!               ! reciprocal of dx_rho
, cy_rho                                                                &
                ! reciprocal of dy_rho
, cx_theta_u(udims%i_start:udims%i_end+1,                               &
             udims%j_start:udims%j_end, model_levels)                   &
                    ! reciprocal of dx_theta_u
, cx_rho_centre(udims%i_start:udims%i_end+1,                            &
                vdims%j_start:vdims%j_end+1, model_levels)              &
                    ! reciprocal of dx_rho_centre
, cy_theta_v(vdims%i_start:vdims%i_end,                                 &
             vdims%j_start:vdims%j_end+1, model_levels)                 &
                    ! reciprocal of dy_theta_v
, cy_rho_centre(udims%i_start:udims%i_end+1,                            &
                vdims%j_start:vdims%j_end+1, model_levels)              &
                    ! reciprocal of dy_rho_centre
, cx2_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          model_levels)                                                 &
                    ! square of cx_rho
, cy2_rho(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          model_levels) ! square of cy_rho

REAL, PARAMETER :: smallp=1.0e-14
                   ! A small number

REAL ::                                                                 &
  weight_pl                                                             &
, weight_min                                                            &
, weight_2pl                                                            &
, weight_2min                                                           &
, ssq11                                                                 &
                      ! ii component of shear stress
, ssq22                                                                 &
                      ! jj component of shear stress
, ssq33                                                                 &
                      ! kk component of shear stress
, ssq13                                                                 &
                      ! ik component of shear stress
, ssq23                                                                 &
                      ! jk component of shear stress
, ssq12                                                                 &
                      ! ij component of shear stress
, ssq12up                                                               &
                      ! part of shear stress
, ssq12k(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)           &
                      ! part of shear stress
, sum_sij
                      ! sum of Sij (ssqij)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_TURB_SMAGORINSKY'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,                             &
!$OMP  dx_theta_u,dy_theta_v,dy_rho_centre,dx_rho,dy_rho,cx_rho,cy_rho, &
!$OMP  dx_rho_centre,delta_zn_u,delta_zn_v,weight_pl,weight_min,        &
!$OMP  weight_2pl,weight_2min,ssq11,ssq22,ssq33,ssq13,ssq23,ssq12up,    &
!$OMP  ssq12,sum_sij,delta_y,delta_x)                                   &
!$OMP SHARED(model_levels,pdims,udims,vdims,r_theta_levels,z,u,v,w,     &
!$OMP  cx_rho_centre,cy_rho_centre,intw_p2u,xi1_p,csxi2_p,cx_theta_u,   &
!$OMP  intw_p2v,xi2_p,cy_theta_v,r_rho_levels,rdz,delta_zn,delta_z,     &
!$OMP  csxi2_v,xi1_u,xi2_v,cx2_rho,cy2_rho,rdzn_u,rdzn_v, ssq12k,       &
!$OMP  intw_v2p, intw_u2p, shear,h2_p_eta,h1_p_eta,delta_smag,          &
!$OMP  mix_factor,blending_option,max_diff,recip_timestep,diff_factor,  &
!$OMP  visc_m,visc_h,h1_p,h2_p,pdims_s)

!$OMP DO SCHEDULE(STATIC)
DO k= 1, model_levels
! delta x on theta-levels at u-points
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end+1
      dx_theta_u = ( intw_p2u(i,1) * h1_p_eta(i,j,k) +                  &
                     intw_p2u(i,2) * h1_p_eta(i+1,j,k) ) *              &
                           ( xi1_p(i+1) - xi1_p(i) )
      cx_theta_u(i,j,k) =  1.0/dx_theta_u
    END DO
  END DO
! delta y on theta-levels at v-points
  DO j = vdims%j_start, vdims%j_end+1
    DO i = vdims%i_start, vdims%i_end
      dy_theta_v = ( intw_p2v(j,1) * h2_p_eta(i,j,k) +                  &
                     intw_p2v(j,2) * h2_p_eta(i,j+1,k) ) *              &
                             ( xi2_p(j+1) - xi2_p(j) )
      cy_theta_v(i,j,k) =  1.0/dy_theta_v
    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end+1
    DO i = udims%i_start, udims%i_end+1
! dx or dy on rho levels at centre of grid (ie b-grid uv point)
      !      cy_rho_centre =  r_ave  !  use as work array
      cy_rho_centre(i,j,k) =                                            &
         intw_p2v(j,2) * (intw_p2u(i,2) * h2_p(i+1,j+1,k) +             &
                          intw_p2u(i,1) * h2_p(i,j+1,k) ) +             &
         intw_p2v(j,1) * (intw_p2u(i,2) * h2_p(i+1,j,k) +               &
                          intw_p2u(i,1) * h2_p(i,j,k) )

      !      dy_rho_centre = r_ave * dphi_p(i,j)
      dy_rho_centre = cy_rho_centre(i,j,k) * (xi2_p(j+1) - xi2_p(j))
      cy_rho_centre(i,j,k) = 1.0 / dy_rho_centre

      !      cx_rho_centre =  r_ave * cos_v_latitude !  use as work array
      cx_rho_centre(i,j,k) =                                            &
         intw_p2v(j,2) * (intw_p2u(i,2) * h1_p(i+1,j+1,k) +             &
                          intw_p2u(i,1) * h1_p(i,j+1,k) ) +             &
         intw_p2v(j,1) * (intw_p2u(i,2) * h1_p(i+1,j,k) +               &
                          intw_p2u(i,1) * h1_p(i,j,k) )

      !      dx_rho_centre =  r_ave * dlambda_p(i)
      dx_rho_centre =  cx_rho_centre(i,j,k) * (xi1_p(i+1) - xi1_p(i))
      cx_rho_centre(i,j,k) =  1.0 / dx_rho_centre
    END DO
  END DO

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      z(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
! dx or dy on rho levels at p-point
      dx_rho = h1_p(i,j,k) * ( xi1_u(i) - xi1_u(i-1) )
      dy_rho = h2_p(i,j,k) * ( xi2_v(j) - xi2_v(j-1) )
      cx_rho = 1.0/dx_rho
      cy_rho =  1.0/dy_rho
      cx2_rho(i,j,k) = cx_rho*cx_rho
      cy2_rho(i,j,k) = cy_rho*cy_rho
    END DO
  END DO

END DO   !  k = 1, model_levels
!$OMP END DO NOWAIT

! Vertical grid spacings used in calculation of rate of strain term
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    delta_z(i,j,1) = r_theta_levels(i,j,1) - r_theta_levels(i,j,0)
    rdz(i,j,1) = 1.0/delta_z(i,j,1)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 2, model_levels
! dz on rho levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      delta_z(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
      rdz(i,j,k) = 1.0/delta_z(i,j,k)
    END DO
  END DO

! dz on theta levels
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      delta_zn(i,j,k-1) = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
    END DO
  END DO

END DO  ! k = 2, model_levels
!$OMP END  DO
! Implicit barrier 

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
! dz on theta levels at u-points
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end+1
      delta_zn_u = intw_p2u(i,1) * delta_zn(i,j,k) +                    &
                   intw_p2u(i,2) * delta_zn(i+1,j,k)
      rdzn_u(i,j,k) = 1.0/delta_zn_u
    END DO
  END DO
! dz on theta levels at v-points
  DO j = vdims%j_start, vdims%j_end+1
    DO i = vdims%i_start, vdims%i_end
      delta_zn_v = intw_p2v(j,2) * delta_zn(i,j+1,k) +                  &
                   intw_p2v(j,1) * delta_zn(i,j,k)
      rdzn_v(i,j,k) = 1.0/delta_zn_v
    END DO
  END DO

END DO
!$OMP END DO
! Implicit barrier 

!--------------------------------------------------------
! As in the LEM code
! _Now calculate half-squared strain rate SSQ on w-points
! _CX=1./DX, CY=1./DY, RDZ(K)=1./DZ(K),
!       RDZN(K) =1./DZN(K)
! _SSQ= 0.5*^^DU_I/DX_J+DU_J/DX_I^^**2
! _SSQIJ= (DU_I/DX_J+DU_J/DX_I)**2
! _Hence SSQ= SUM_sij(I,J) {0.5*(SSQIJ)}
! _Note that for a simple shear S, SSQ=S**2
!   (as in the notation of Mason and Cullen 1986)
!--------------------------------------------------------
k = 1
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! ssq12 = (du/dy+dv/dx)^2 on p-point and 1st rho-level
!                         to be averaged to theta-level later
    ssq12k(i,j) = intw_v2p(j,1) * intw_u2p(i,1) *                       &
        ((u(i-1,j,k)-u(i-1,j-1,k)) * cy_rho_centre(i-1,j-1,k)           &
        +(v(i,j-1,k)-v(i-1,j-1,k)) * cx_rho_centre(i-1,j-1,k) )**2 +    &
                  intw_v2p(j,1) * intw_u2p(i,2) *                       &
        ((u(i,j,k)-u(i,j-1,k)) * cy_rho_centre(i,j-1,k)                 &
        +(v(i+1,j-1,k)-v(i,j-1,k)) * cx_rho_centre(i,j-1,k) )**2 +      &
                  intw_v2p(j,2) * intw_u2p(i,1) *                       &
        ((u(i-1,j+1,k)-u(i-1,j,k)) * cy_rho_centre(i-1,j,k)             &
        +(v(i,j,k)-v(i-1,j,k)) * cx_rho_centre(i-1,j,k) )**2 +          &
                  intw_v2p(j,2) * intw_u2p(i,2) *                       &
        ((u(i,j+1,k)-u(i,j,k)) * cy_rho_centre(i,j,k)                   &
        +(v(i+1,j,k)-v(i,j,k)) * cx_rho_centre(i,j,k) )**2

  END DO   ! on i
END DO    ! on j
!$OMP END DO NOWAIT


DO k = 1, model_levels -1
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      ! Vertical weights
! interpolation of rho to theta levels
      weight_pl = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/        &
                  (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
      weight_min = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/     &
                   (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
! interpolation of theta to theta levels
      weight_2pl = delta_z(i,j,k)/(delta_z(i,j,k)+delta_z(i,j,k+1))
      weight_2min= delta_z(i,j,k+1)/(delta_z(i,j,k)+delta_z(i,j,k+1))

      ! ssq11 = 2 * backward difference (du/dx)^2 averaged to w point
      ssq11 =  2.0*(                                                    &
        weight_pl*cx2_rho(i,j,k+1)*(u(i,j,k+1)-u(i-1,j,k+1))**2 +       &
        weight_min*cx2_rho(i,j,k)*(u(i,j,k)-u(i-1,j,k))**2   )

      ! ssq22 = 2 * backward difference (dv/dy)^2 averaged to w point
      ssq22 = 2.0*(                                                     &
        weight_pl*cy2_rho(i,j,k+1)*(v(i,j,k+1)-v(i,j-1,k+1))**2 +       &
        weight_min*cy2_rho(i,j,k)*(v(i,j,k)-v(i,j-1,k))**2  )

      ! ssq33 =  2 * backward difference (dw/dz)^2 averaged to w point
      ssq33 =  2.0*(                                                    &
        weight_2pl*((w(i,j,k)-w(i,j,k-1))*rdz(i,j,k))**2 +              &
        weight_2min*((w(i,j,k+1)-w(i,j,k))*rdz(i,j,k+1))**2  )

      ! ssq13 =  (du/dz+dw/dx)^2 averaged to w point
      ! ssq31 = ssq13
      ssq13 = intw_u2p(i,1) * ( (u(i-1,j,k+1)-u(i-1,j,k)) *             &
                                                  rdzn_u(i-1,j,k) +     &
                     (w(i,j,k)-w(i-1,j,k)) * cx_theta_u(i-1,j,k) )**2 + &
              intw_u2p(i,2) * ( (u(i,j,k+1) - u(i,j,k)) *               &
                                                  rdzn_u(i,j,k) +       &
                     (w(i+1,j,k)-w(i,j,k)) * cx_theta_u(i,j,k) )**2

      ! ssq23 =  (dw/dy+dv/dz)^2 averaged to w point
      ! ssq32 = ssq23
      ssq23 = intw_v2p(j,2) * ( (w(i,j+1,k)-w(i,j,k)) *                 &
                                                  cy_theta_v(i,j,k) +   &
                  (v(i,j,k+1) - v(i,j,k)) * rdzn_v(i,j,k))**2 +         &
              intw_v2p(j,1) * ( (w(i,j,k)-w(i,j-1,k)) *                 &
                                                  cy_theta_v(i,j-1,k) + &
                  (v(i,j-1,k+1) - v(i,j-1,k)) * rdzn_v(i,j-1,k))**2

      ! ssq12 = (du/dy+dv/dx)^2 on p-point and k+1 rho-level
      ! ssq21 = ssq12
      ssq12up = intw_v2p(j,1) * intw_u2p(i,1) *                         &
              ((u(i-1,j,k+1)-u(i-1,j-1,k+1)) * cy_rho_centre(i-1,j-1,k+1)     &
              +(v(i,j-1,k+1)-v(i-1,j-1,k+1)) * cx_rho_centre(i-1,j-1,k+1) )**2&
                    + intw_v2p(j,1) * intw_u2p(i,2) *                         &
              ((u(i,j,k+1)-u(i,j-1,k+1))     * cy_rho_centre(i,j-1,k+1)       &
              +(v(i+1,j-1,k+1)-v(i,j-1,k+1)) * cx_rho_centre(i,j-1,k+1) )**2  &
                    + intw_v2p(j,2) * intw_u2p(i,1) *                         &
              ((u(i-1,j+1,k+1)-u(i-1,j,k+1)) * cy_rho_centre(i-1,j,k+1)       &
              +(v(i,j,k+1)-v(i-1,j,k+1)) *     cx_rho_centre(i-1,j,k+1) )**2  &
                    + intw_v2p(j,2) * intw_u2p(i,2) *                         &
              ((u(i,j+1,k+1)-u(i,j,k+1)) * cy_rho_centre(i,j,k+1)             &
              +(v(i+1,j,k+1)-v(i,j,k+1)) * cx_rho_centre(i,j,k+1) )**2

      ! average to theta level k
      ssq12 = weight_pl * ssq12up + weight_min * ssq12k(i,j)
      ! save the k+1 value for use at the next iteration of k
      ssq12k(i,j) = ssq12up

      ! sum_ssqij would be 2*ssq11 + 2*ssq22 + 2*ssq33 + ssq13
      !                      + ssq23 + ssq12 + ssq31 + ssq32 + ssq21
      ! hence by the equivalence of the off-diagonal terms (ssq12 etc),
      ! ssq = 0.5 (sum_ssqij) is given by sum_sij below
      !
      ! the extra factor of 2 in the diagonal terms comes from the fact
      ! that they should be (2*du/dx)^2 etc, not 2*(du/dx)^2
      sum_sij = ssq11 + ssq22 + ssq33 +                          &
                                   ssq13 + ssq23 + ssq12 + smallp
      shear(i,j,k) = SQRT(sum_sij)

    END DO   ! on i
  END DO    ! on j
!$OMP END DO NOWAIT
END DO ! k = 1, model_levels -1


! Horizontal grid spacings used in calculation of mixing length and
! maximum diffusivity.
!
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

    !   delta_y = Earth_radius * dphi_v(i,j-1) 
    delta_y = h2_p_eta(i,j,0) * ( xi2_v(j) - xi2_v(j-1))

    !   delta_x = Earth_radius * dlambda_u(i-1) * cos_v_latitude(i,j)
    delta_x = h1_p_eta(i,j,0) * ( xi1_u(i) - xi1_u(i-1))

    IF (blending_option == off) THEN
      delta_smag(i,j) = MAX(delta_x, delta_y)
    ELSE
      delta_smag(i,j) = MIN(delta_x, delta_y)
    END IF
    !
    ! maximum diffusion coefficient allowed in this run
    max_diff(i,j) = 0.25 * recip_timestep * diff_factor /               &
                                      ( 1.0 / (delta_x * delta_x) +     &
                                        1.0 / (delta_y * delta_y) )
  END DO  ! on i
END DO  ! on j
!$OMP END DO NOWAIT

DO k = 1, model_levels-1
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      visc_m(i,j,k)    = shear(i,j,k)
      visc_h(i,j,k)    = shear(i,j,k)
      ! visc are now S and still need to be multiplied by
      ! stability functions and lambda^2 (done in BL)
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP END PARALLEL 

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_turb_smagorinsky
