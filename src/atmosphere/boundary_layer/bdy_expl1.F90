! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BDY_EXPL1----------------------------------------------

!  Purpose: Calculate explicit scalling parameters and additional
!           boundary layer information required by the surface
!           exchange scheme

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer

!---------------------------------------------------------------------
SUBROUTINE bdy_expl1 (                                                  &
! IN values defining vertical grid of model atmosphere :
 bl_levels, p_theta_levels, rho_wet_rsq, rho_wet, rho_dry,              &
! IN U, V and W momentum fields.
 u_px, v_px,                                                            &
! IN cloud data :
 cf_bulk,q,qcf,qcl,t,                                                   &
! OUT
 dtrdz_charney_grid,rdz_charney_grid,dtrdz_u,dtrdz_v,                   &
 rdz_u,rdz_v,rho_mix,rho_mix_tq,dzl_charney,rdz,                        &
 zblend_tq,z1_tq_top,zblend_uv,z1_uv_top,                               &
 k_blend_tq,k_blend_uv,qw,tl,qw_blend,tl_blend,u_blend_px,v_blend_px,   &
 bt,bq,bt_blend,bq_blend,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt   &
 )

USE atm_fields_bounds_mod, ONLY:                                        &
 udims, vdims, pdims, pdims_s, tdims, tdims_l
USE planet_constants_mod, ONLY: lcrcp, lsrcp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE vertical_diffs_mod, ONLY: vertical_diffs
USE buoy_tq_mod, ONLY: buoy_tq

IMPLICIT NONE

!  IN data
! Defining vertical grid of model atmosphere.
INTEGER, INTENT(IN) ::                                                  &
 bl_levels
                                 ! IN Max. no. of "boundary" levels

REAL, INTENT(IN) ::                                                     &
  p_theta_levels(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end, 0:bl_levels+1),             &
  rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,bl_levels+1),               &
                                             ! IN Density * R**2
  rho_wet(pdims%i_start:pdims%i_end,                                    &
          pdims%j_start:pdims%j_end, bl_levels+1),                      &
                      ! IN wet density on rho levels (kg/m3)
  rho_dry(pdims%i_start:pdims%i_end,                                    &
          pdims%j_start:pdims%j_end, bl_levels+1),                      &
                      ! IN dry density on rho levels (kg/m3)
 u_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
     bl_levels),                                                        &
 v_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
     bl_levels)

! Cloud data.
REAL, INTENT(IN) ::                                                     &
 cf_bulk(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end, bl_levels),                         &
                                        ! IN Cloud fraction (decimal).
 qcf(tdims_l%i_start:tdims_l%i_end,                                     &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),          &
                                   ! IN Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,                                     &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),          &
                                   ! IN Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,                                       &
   tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),            &
                                   ! IN specific humidity
 t(tdims%i_start:tdims%i_end,                                           &
   tdims%j_start:tdims%j_end,bl_levels)       ! IN temperature

! OUT data

INTEGER, INTENT(OUT) ::                                                 &

 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                         ! OUT Theta level for blending height.
 k_blend_uv(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
                         ! OUT Rho level for blending height.
REAL, INTENT(OUT) ::                                                    &
 dtrdz_charney_grid(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,bl_levels),               &
                              ! OUT dt/(rho*r*r*dz) for scalar
                              !     flux divergence
 rdz_charney_grid(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,bl_levels),                 &
                              ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                           ! OUT dt/(rho*r*r*dz) for
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! OUT U,V flux divergence
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT 1/(Z_U(K)-Z_U(K-1)) for K > 1
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT 1/(Z_V(K)-Z_V(K-1)) for K > 1
   rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels+1),                                                &
                                ! OUT density on UV (ie. rho) levels;
                                !    used in RHOKH so dry density if
                                !    Lq_mix_bl is true
   rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
             bl_levels),                                                &
                                ! OUT density on TQ (ie. theta) levels;
                                !    used in non-turb flux integration
                                !    so dry density if Lq_mix_bl is true
   dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               bl_levels),                                              &
                                ! OUT DZL(,K) is depth in m of theta
!                                 level K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
   rdz( pdims_s%i_start:pdims_s%i_end,                                  &
        pdims_s%j_start:pdims_s%j_end, bl_levels ),                     &
                                ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   zblend_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                                ! OUT Height of theta level blending height
   zblend_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                ! OUT Height of u,v level blending height
   u_blend_px(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end),                           &
                                ! OUT
   v_blend_px(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end),                           &

   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT Ice/liquid water temperature
   qw_blend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                ! OUT water content at blend height
   tl_blend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                ! OUT ice/liquid water temperature
!                                 at blend height
   bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bt_blend(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                ! OUT buoyancy parameter at blend height
   bq_blend(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                ! OUT buoyancy parameter at blend height
   bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT A buoyancy parameter for cloudy
                                ! air on p,T,q-levels (full levels).
   bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT A buoyancy parameter for cloudy
                                ! air on p,T,q-levels (full levels).
   bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! OUT A grid-box mean buoyancy parameter
                                ! on p,T,q-levels (full levels).
   bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! OUT A grid-box mean buoyancy parameter
                                ! on p,T,q-levels (full levels).
   a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                ! OUT Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                                ! OUT Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                ! OUT Derivative of q_SAT w.r.t. T
REAL, INTENT(OUT) :: z1_uv_top(pdims%i_start:pdims%i_end,               &
                               pdims%j_start:pdims%j_end)
                              ! Height of top of lowest uv-layer
                              ! above the surface
REAL, INTENT(OUT) :: z1_tq_top(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end)
                              ! Height of top of lowest Tq-layer
                              ! above the surface
!-----------------------------------------------------------------------
!  Workspace :-

!  Local scalars :-

INTEGER ::                                                              &
 i,j,                                                                   &
                ! LOCAL Loop counter (horizontal field index).
 k          ! LOCAL Loop counter (vertical level index).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BDY_EXPL1'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Calculate layer depths and heights, and construct wind fields on
!     P-grid.
!-----------------------------------------------------------------------

CALL vertical_diffs(                                                    &
! IN Field dimensions/logicals
    bl_levels,                                                          &
! IN Fields.
    rho_wet_rsq, rho_wet, rho_dry,                                      &
! OUT Vertical differences required by physics.
    rho_mix, rho_mix_tq,                                                &
    dzl_charney, rdz,                                                   &
    zblend_uv, z1_uv_top, zblend_tq, z1_tq_top,                         &
    k_blend_tq, k_blend_uv,                                             &
    rdz_charney_grid, dtrdz_charney_grid,                               &
    dtrdz_u, dtrdz_v, rdz_u, rdz_v                                      &
    )

!-----------------------------------------------------------------------
! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                        &
!$OMP&         PRIVATE(i,j,k)                                           &
!$OMP&         SHARED(bl_levels,tdims, qw,q,qcl,qcf,tl,t,lcrcp,lsrcp)
DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
      tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! Extract momentum, temp and moisture at blending height
!-----------------------------------------------------------------------
DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    u_blend_px(i,j) = u_px(i,j,k_blend_uv(i,j))
    v_blend_px(i,j) = v_px(i,j,k_blend_uv(i,j))
  END DO
END DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    tl_blend(i,j) = tl(i,j,k_blend_tq(i,j))
    qw_blend(i,j) = qw(i,j,k_blend_tq(i,j))
  END DO
END DO

!-----------------------------------------------------------------------
! Calculate buoyancy parameters BT and BQ.
!-----------------------------------------------------------------------

CALL buoy_tq (                                                          &
! IN dimensions/logicals
   bl_levels,                                                           &
! IN fields
   p_theta_levels,t,q,qcf,qcl,cf_bulk,                                  &
! OUT fields
   bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                   &
    )

!-----------------------------------------------------------------------
! Extract buoyancy parameters at blending height
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    bt_blend(i,j) = bt_gb(i,j,k_blend_tq(i,j))
    bq_blend(i,j) = bq_gb(i,j,k_blend_tq(i,j))
  END DO
END DO


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bdy_expl1
