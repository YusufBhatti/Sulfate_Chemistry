! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PURPOSE: Calculate vertical differences required by BL scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Boundary Layer

! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.

MODULE vertical_diffs_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'VERTICAL_DIFFS_MOD'
CONTAINS

SUBROUTINE vertical_diffs (                                             &
! IN Field dimensions/logicals
  bl_levels,                                                            &
! IN Fields.
  rho_wet_rsq, rho_wet, rho_dry,                                        &
! OUT Vertical differences required by physics.
  rho_mix, rho_mix_tq,                                                  &
  dzl_charney, rdz,                                                     &
  zblend_uv, z1_uv_top, zblend_tq, z1_tq_top,                           &
  k_blend_tq,k_blend_uv,                                                &
  rdz_charney_grid,                                                     &
  dtrdz_charney_grid,                                                   &
  dtrdz_u, dtrdz_v, rdz_u, rdz_v                                        &
  )

USE atm_fields_bounds_mod, ONLY:                                        &
 udims, vdims, pdims_s, pdims, pdims_l, tdims
USE bl_option_mod, ONLY:                                                &
   on, blend_height_opt, h_blend_fix, blend_level1, blend_levelk
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE jules_surface_mod, ONLY: i_modiscopt
USE level_heights_mod, ONLY:                                            &
  r_theta_levels, r_rho_levels
USE p_to_t_mod, ONLY: p_to_t
USE p_to_u_mod, ONLY: p_to_u
USE p_to_v_mod, ONLY: p_to_v
USE timestep_mod, ONLY: timestep
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

! LOGICAL SWITCHES
INTEGER, INTENT(IN) :: bl_levels

REAL, INTENT(IN) ::                                                     &
  rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,bl_levels+1),               &
                      ! wet density times r^2 on rho levels (kg/m3)
  rho_wet(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1),                                                 &
                      ! wet density on rho levels (kg/m3)
  rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1)
                      ! dry density on rho levels (kg/m3)

! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

REAL, INTENT(OUT) ::                                                    &
  rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
         bl_levels+1),                                                  &
                              ! OUT density on UV (ie. rho) levels;
                              !    used in RHOKH so dry density if
                              !    l_mr_physics is true

  rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             bl_levels),                                                &
                              ! OUT density on TQ (ie. theta) levels;
                              !    used in non-turb flux integration
                              !    so dry density if l_mr_physics is true
  dtrdz_charney_grid (tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels),             &
  dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              bl_levels),                                               &
  rdz( pdims_s%i_start:pdims_s%i_end,                                   &
       pdims_s%j_start:pdims_s%j_end, bl_levels ),                      &
                              ! OUT  Reciprocal of distance between
  rdz_charney_grid(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,bl_levels),                &
  zblend_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
  zblend_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

INTEGER, INTENT(OUT) ::                                                 &
  k_blend_uv(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end),                            &
  k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(OUT) :: z1_uv_top(pdims%i_start:pdims%i_end,               &
                               pdims%j_start:pdims%j_end)
                              ! Top of lowest uv-layer
REAL, INTENT(OUT) :: z1_tq_top(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end)
                              ! Top of lowest Tq-layer

REAL, INTENT(OUT) ::                                                    &
            ! quantities interpolated to u and v grids.
  rdz_u (udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
        2:bl_levels),                                                   &
  rdz_v (vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
        2:bl_levels),                                                   &
  dtrdz_u (udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
           bl_levels),                                                  &
  dtrdz_v (vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
           bl_levels)

! LOCAL VARIABLES.
INTEGER ::                                                              &
  k,                                                                    &
  i,j

REAL ::                                                                 &
  dz_p,                                                                 &
  dz_t

REAL ::                                                                 &
  work(pdims_s%i_start:pdims_s%i_end,                                   &
       pdims_s%j_start:pdims_s%j_end, bl_levels)

REAL ::                                                                 &
    h_blend(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VERTICAL_DIFFS'

! ----------------------------------------------------------------------
! CALCULATE rho on layer centres and layer boundaries
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP&         PRIVATE(i,j,k,dz_p,dz_t)
IF ( l_mr_physics ) THEN
      ! conservation of mixing ratios requires use of dry density
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels+1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rho_mix(i,j,k) = rho_dry(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels+1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rho_mix(i,j,k) = rho_wet(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO
END IF

! Interpolate Rho_mix to temperature levels for Rho_mix_tq
!  - can't think of a better name but this will be wet or dry
!    depending on l_mr_physics

CALL p_to_t (                                                           &
! IN Field dimensions and pointers.
    pdims%i_end, pdims%j_end, pdims_l%halo_i, pdims_l%halo_j, 0, 0,     &
    bl_levels,                                                          &
! IN Vertical coordinate levels.
    r_theta_levels, r_rho_levels,                                       &
! IN field to be interpolated.
    rho_mix,                                                            &
! OUT Interpolated field.
    rho_mix_tq                                                          &
   )

! Need to barrier here to make sure previous routine has fully completed
! before we use the results
!$OMP BARRIER


!-----------------------------------------------------------------
! Calculate k_blend_tq(i,j) and k_blend_uv(i,j)
! The model levels used for implicitly coupling to the surface
!-----------------------------------------------------------------
SELECT CASE (blend_height_opt)

CASE (blend_level1) ! Use bottom model level
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      k_blend_uv(i,j) = 1
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      k_blend_tq(i,j) = 1
    END DO
  END DO
!$OMP END DO

CASE (blend_levelk) ! Use specified blending height (h_blend_fix)
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ! Set up h_blend(i,j) from h_blend_fix
      h_blend(i,j) = h_blend_fix
      k_blend_tq(i,j) = 1
      DO WHILE ( (r_theta_levels(i,j,k_blend_tq(i,j)) -                 &
                  r_theta_levels (i,j,0)) < h_blend(i,j) )
        ! Pick out the first theta model level above h_blend
        k_blend_tq(i,j) = k_blend_tq(i,j) + 1
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      ! Set up h_blend(i,j) from h_blend_fix
      h_blend(i,j) = h_blend_fix
      k_blend_uv(i,j) = 1
      DO WHILE ( (r_rho_levels(i,j,k_blend_uv(i,j)) -                   &
                  r_theta_levels (i,j,0)) < h_blend(i,j) )
        ! Pick out the first rho model level above h_blend
        k_blend_uv(i,j) = k_blend_uv(i,j) + 1
      END DO
    END DO
  END DO
!$OMP END DO

END SELECT

! ----------------------------------------------------------------------
! CALCULATE DTRDZ, DZL, RDZ
! ----------------------------------------------------------------------
k = 1
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    dz_p = r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)
    dz_t = r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1)
    zblend_uv(i,j) = r_rho_levels(i,j,k_blend_uv(i,j)) -                &
                     r_theta_levels (i,j,k-1)
    zblend_tq(i,j) = r_theta_levels(i,j,k_blend_tq(i,j)) -              &
                     r_theta_levels (i,j,k-1)
    rdz_charney_grid(i,j,k) = 1.0/ dz_p
        ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
        ! equation for scalars - hence dry density if l_mr_physics
        ! is true, wet otherwise
    dtrdz_charney_grid(i,j,k) = timestep/                               &
              ( r_theta_levels (i,j,k)*r_theta_levels (i,j,k)*          &
                                      rho_mix_tq(i,j,k)*dz_t )
    ! Dzl_charney( ,,1) is such that
    ! 0.5 * DZL_Charney = height of first temperature level
    dzl_charney(i,j,k) = 2.0 *                                          &
               (r_theta_levels(i,j,1) - r_theta_levels(i,j,0))
  END DO
END DO
!$OMP END DO NOWAIT

!     Additionally calculate the tops of the lowest layers if using
!     fully conservative discretization.
IF (i_modiscopt == on) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      z1_uv_top(i,j) = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
      z1_tq_top(i,j) = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k-1)
    END DO
  END DO
!$OMP END DO
END IF

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      dz_p = r_theta_levels (i,j,k) - r_theta_levels(i,j,k-1)
      dz_t = r_rho_levels (i,j,k+1) - r_rho_levels(i,j,k)
      dzl_charney(i,j,k) = dz_t
      rdz_charney_grid(i,j,k) = 1.0/ dz_p
        ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
        ! equation for scalars - hence dry density if l_mr_physics
        ! is true, wet otherwise
      dtrdz_charney_grid(i,j,k) = timestep/                             &
              ( r_theta_levels (i,j,k)*r_theta_levels (i,j,k)*          &
                                      rho_mix_tq(i,j,k)*dz_t )
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! ----------------------------------------------------------------------
! horizontal interpolations to momentum points.
! Note the interpolation routines are threadsafe and parallel
! ----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k = 1, bl_levels
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      dz_p = r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)
      work(i,j,k) = timestep/( rho_wet_rsq(i,j,k) * dz_p )
    END DO
  END DO
END DO
!$OMP END DO

CALL p_to_u (work,                                                      &
             pdims_s%i_start,pdims_s%i_end,                             &
             pdims_s%j_start,pdims_s%j_end,                             &
             udims%i_start,udims%i_end,                                 &
             udims%j_start,udims%j_end,                                 &
             1,bl_levels,dtrdz_u)

CALL p_to_v (work,                                                      &
             pdims_s%i_start,pdims_s%i_end,                             &
             pdims_s%j_start,pdims_s%j_end,                             &
             vdims%i_start,vdims%i_end,                                 &
             vdims%j_start,vdims%j_end,                                 &
             1,bl_levels,dtrdz_v)

! ----------------------------------------------------------------------
! Calculation of 1/dz (rdz) arrays on momentum flux levels
! ----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
    ! Temporarily store dz in the array for 1/dz
    rdz(i,j,1) = r_rho_levels(i,j,1) - r_theta_levels(i,j,0)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      ! Temporarily store dz in the array for 1/dz
      rdz(i,j,k) = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END DO

CALL p_to_u (rdz(:,:,2:bl_levels),                                      &
             pdims_s%i_start,pdims_s%i_end,                             &
             pdims_s%j_start,pdims_s%j_end,                             &
             udims%i_start,udims%i_end,                                 &
             udims%j_start,udims%j_end,                                 &
             2,bl_levels,rdz_u)

CALL p_to_v (rdz(:,:,2:bl_levels),                                      &
             pdims_s%i_start,pdims_s%i_end,                             &
             pdims_s%j_start,pdims_s%j_end,                             &
             vdims%i_start,vdims%i_end,                                 &
             vdims%j_start,vdims%j_end,                                 &
             2,bl_levels,rdz_v)
! Note: dz is now also temporarily stored in the 1/dz arrays for u,v points.
! Done this way so that we interpolate as dz, then convert to 1/dz.

! Further note: the above interpolation is needless, as we could have just
! vertically-differenced the arrays of r on u,v points,
! r_at_u, r_at_v stored in level_heights_mod.  Then rdz would not need haloes.
! Differencing these pre-interpolated heights would be mathematically identical
! to interpolating the vertical differences as above, however the changed order
! of calculations would change the rounding error, so would change answers.

! Need to make sure all parallelised calculations in these routines are
! complete before we use the computed values
!$OMP BARRIER

! Take reciprocals to complete calculation of rdz arrays

!$OMP DO SCHEDULE(STATIC)
DO k = 1, bl_levels
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      rdz(i,j,k) = 1.0/rdz(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      rdz_u(i,j,k) = 1.0/rdz_u(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      rdz_v(i,j,k) = 1.0/rdz_v(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vertical_diffs
END MODULE vertical_diffs_mod
