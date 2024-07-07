! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE Restart_Dump
!   PURPOSE:-           To create restart dump for subsequent runs
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!

MODULE restart_dump_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RESTART_DUMP_MOD'
CONTAINS

SUBROUTINE restart_dump                                                       &
  ( row_length, rows, land_points, nbl_levs                                   &
  , nsoilt_levs, nsoilm_levs, n_cca_lev, land_sea_mask, u, v, w               &
  , t, theta, q, qcl, qcf, layer_cloud, p, rho, t_deep_soil, smc, canopy_gb   &
  , snodep, tstar, zh, z0msea, cca, cca_dp, cca_md, cca_sh                    &
  , iccb, icct, smcl, bl_w_var )
! Note: if you add more fields to the dump used here (array resdump),
! you need to increase the declared array size for resdump, set in the
! module scm_utils.
! You also need to modify the routine which reads the dump consistently;
! dumpinit (called from subroutine runinit and from scm_main).

USE scm_utils, ONLY: resdump
USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

INTEGER ::       &
  row_length     &! In x dimension of arrays
, rows           &! In y dimension of arrays
, land_points    &! In no of land_points
, nbl_levs       &! In Number of Boundary layer levels
, nsoilt_levs    &! In Number of soil temperature levels
, nsoilm_levs    &! In Number of soil moisture levels
, n_cca_lev       ! In no of cca levels
!
!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------
!     Primary model variables + T

INTEGER ::                     &
  iccb(row_length, rows)       &! In Convective cloud base
, icct(row_length, rows)        ! In Convective cloud top

REAL ::                                     &
  canopy_gb(land_points)                    &! In Canopy water content (kg/m^2)
, cca(row_length,rows,n_cca_lev)            &! In Convective cloud amount
, cca_dp(row_length,rows,n_cca_lev)         &! In deep Conv cloud amount
, cca_md(row_length,rows,n_cca_lev)         &! In mid-level Conv cloud amount
, cca_sh(row_length,rows,n_cca_lev)         &! In shallow Conv cloud amount
, layer_cloud(row_length,rows,model_levels) &! In Layer cloud amount (decima
, q(row_length,rows,model_levels)           &! In Specific humidity (kg/kg)
, qcf(row_length,rows,model_levels)         &! In Cloud ice content (kg/kg)
, qcl(row_length,rows,model_levels)         &! In Cloud water content (kg/kg)
, smc(land_points)                          &! In Soil moisture content (kg/m^2)
, smcl(land_points,nsoilm_levs)           &! In Soil moisture in levels (kg/m^2)
, snodep(row_length,rows)                   &! In Snow depth (kg/m^2)
, t(row_length,rows,model_levels)           &! In Temperature at each level (K)
, t_deep_soil(land_points,nsoilt_levs)      &! In Deep soil temperatures K
, theta(row_length,rows,model_levels)       &! In Potential temperature (K)
, tstar(row_length,rows)                    &! In Surface temperature (K)
, rho(row_length,rows,model_levels)         &
, p(row_length,rows,model_levels+1)         &! In Pressure (mb)
, u(row_length,rows,model_levels)           &! In Zonal wind (m/s^2)
, v(row_length,rows,model_levels)           &! In Meridional wind (m/s^2)
, w(row_length,rows,0:model_levels)         &! In vertical velocity
, zh(row_length,rows)                       &! In Height above surface of
                                             !    top of boundary layer (m)
, z0msea(row_length,rows)                   &! In Sea surface roughness length
, bl_w_var(row_length,rows,model_levels)     ! In BL W variance 

LOGICAL ::                            &
  Land_sea_mask(row_length,rows)       ! In True if land point

INTEGER ::                            &
  i, j, k                             &! Loop counter
, icount                               ! Counter

INTEGER ::                            &
  land_cnt                             ! land point counter

CHARACTER(LEN=*), PARAMETER :: RoutineName='RESTART_DUMP'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

land_cnt = 0

DO k=1, rows
  DO j=1, row_length

    IF (land_sea_mask(j,k)) THEN
      land_cnt = land_cnt + 1
    END IF

    DO i=1, model_levels
      resdump(j,k, i) = resdump(j,k,i) + u(j,k,i)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + v(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels
      resdump(j,k,i) = resdump(j,k,i) + w(j,k,i-icount)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + t(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + theta(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + q(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + qcl(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(j,k,i) = resdump(j,k,i) + qcf(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount,icount + model_levels-1
      resdump(k, j, i) = resdump(j,k,i) +                           &
                layer_cloud(j,k,i-icount + 1)
    END DO

    icount = i
    IF (land_sea_mask(j,k)) THEN
      DO i=icount,icount + nsoilt_levs-1
        resdump(j,k,i) = resdump(j,k,i)                             &
                       + t_deep_soil(land_cnt,i-icount + 1)
      END DO
      icount = i
    END IF

    DO i=icount,icount + model_levels+1
      resdump(j,k,i) = resdump(j,k,i) + p(j,k,i-icount + 1)
    END DO

    icount = i
    DO i=icount, icount + model_levels
      resdump(j,k,i) = resdump(j,k,i) + rho(j,k,i-icount + 1)
    END DO

    icount = i
    IF (land_sea_mask(j,k)) THEN
      resdump(j,k,icount) = resdump(j,k,icount)                     &
                          + smc(land_cnt)
      icount = icount + 1
      resdump(j,k,icount) = resdump(j,k,icount)                     &
                          + canopy_gb(land_cnt)
      icount = icount + 1
    END IF

    resdump(j,k,icount) = resdump(j,k,icount) + snodep(j,k)
    icount = icount + 1
    resdump(j,k,icount) = resdump(j,k,icount) + tstar(j,k)
    icount = icount + 1
    resdump(j,k,icount) = resdump(j,k,icount) + zh(j,k)
    icount = icount + 1
    resdump(j,k,icount) = resdump(j,k,icount) + z0msea(j,k)
    icount = icount + 1

    DO i=icount,icount + n_cca_lev
      resdump(j,k,i) = resdump(j,k,i) + cca(j,k,i-icount + 1)
    END DO

    icount = icount + 1
    resdump(j,k,icount) = resdump(j,k,icount) + REAL(iccb(j,k))
    icount = icount + 1
    resdump(j,k,icount) = resdump(j,k,icount) + REAL(icct(j,k))

    IF (land_sea_mask(j,k)) THEN
      DO i=1, nsoilm_levs
        resdump(j,k,icount + i) = resdump(j,k,icount + i)           &
                                + smcl(land_cnt,i)
      END DO
    END IF

    icount = icount + 1
    DO i=icount,icount + n_cca_lev
      resdump(j,k,i) = resdump(j,k,i) + cca_dp(j,k,i-icount + 1)
    END DO

    icount = icount + 1
    DO i=icount,icount + n_cca_lev
      resdump(j,k,i) = resdump(j,k,i) + cca_md(j,k,i-icount + 1)
    END DO

    icount = icount + 1
    DO i=icount,icount + n_cca_lev
      resdump(j,k,i) = resdump(j,k,i) + cca_sh(j,k,i-icount + 1)
    END DO

    icount = icount + 1
    DO i=icount,icount + model_levels
      resdump(j,k,i) = resdump(j,k,i) + bl_w_var(j,k,i-icount + 1)
    END DO

  END DO                   ! j
END DO                   ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE restart_dump

END MODULE restart_dump_mod
