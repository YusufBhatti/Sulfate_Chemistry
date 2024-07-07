! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE DumpInit
!   PURPOSE:- To initialise primary variables from RESDUMP read in
!             previously from tape

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

MODULE dumpinit_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DUMPINIT_MOD'
CONTAINS

SUBROUTINE dumpinit                                                           &
  ( row_length, rows, land_points, nbl_levs                                   &
  , nsoilt_levs, nsoilm_levs, n_cca_lev, land_sea_mask, u, v, w               &
  , t, theta, q, qcl, qcf, layer_cloud, p, rho, t_deep_soil, smc, canopy_gb   &
  , snodep, tstar, zh, z0msea, cca, cca_dp, cca_md, cca_sh                    &
  , rccb, rcct, smcl, bl_w_var )
! Note: if you add more fields to the dump used here (array resdump),
! you need to increase the declared array size for resdump, set in the
! module scm_utils.
! You also need to modify the routine which saves the dump consistently;
! restart_dump (called from scm_main).
! You also need to modify both calls to dumpinit, from run_init and scm_main.

USE scm_utils, ONLY: resdump
USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(In)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) :: &
  row_length           &! x dimension
, rows                 &! y dimension
, nbl_levs             &! Number of Boundary layer levels
, nsoilt_levs          &! Number of soil temperature levels
, nsoilm_levs          &! Number of soil moisture levels
, land_points          &! Number of land_points
, n_cca_lev             ! Number of levels of cca


LOGICAL, INTENT(IN) ::                &
  land_sea_mask(row_length,rows)       ! True if land point

!-----------------------------------------------------------------------------
! Arguments with INTENT(Out)
!-----------------------------------------------------------------------------
REAL, INTENT(OUT) ::                  &
  u(row_length,rows,model_levels)     &! Zonal wind (m/s^2)
, v(row_length,rows,model_levels)     &! Meridional wind (m/s^2)
, w(row_length,rows,0:model_levels)   &! vertical velocity
, p( row_length,rows,model_levels+1)  &
, rho(row_length,rows,model_levels)   &
, t(row_length,rows,model_levels)     &! Temperature at each level
, theta(row_length,rows,model_levels) &! Potential temp. (K)
, q(row_length,rows,model_levels)     &! Specific humidity (kg/kg)
, layer_cloud(row_length,rows,model_levels) &! Layer cloud amount (decimal)
, qcl(row_length,rows,model_levels)   &! Cloud water content (kg/kg)
, qcf(row_length,rows,model_levels)   &! Cloud ice content (kg/kg)
, t_deep_soil(land_points,nsoilt_levs)&! Deep soil temperatures
, smc(land_points)                    &! Soil moisture content (kg/m^2)
, smcl(land_points,nsoilm_levs)       &! Soil moisture in layers (kg/m^2)
, rccb(row_length,rows)               &! Convective cloud base
, rcct(row_length,rows)               &! Convective cloud top
, cca(row_length,rows,n_cca_lev)      &! Convective cloud amount
, cca_dp(row_length,rows,n_cca_lev)   &! Deep Conv cloud amount
, cca_md(row_length,rows,n_cca_lev)   &! Mid-level Conv cloud amount
, cca_sh(row_length,rows,n_cca_lev)   &! Shallow Conv cloud amount
, bl_w_var(row_length,rows,model_levels) &! BL w variance
, canopy_gb(land_points)              &! Canopy water content (kg/m^2)
, snodep(row_length,rows)             &! Snow depth (kg/m^2)
, tstar(row_length,rows)              &! Surface temperature (K)
, z0msea(row_length,rows)             &! Sea surface roughness length
, zh(row_length,rows)                  ! Height above surface of to
                                       ! boundary layer (m)


!-----------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------------

INTEGER ::           &
  i,j,k              &! Loop counter
, icount             &! Counter
, land_cnt

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUMPINIT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

land_cnt = 0

DO k=1, rows
  DO j=1, row_length

    IF (land_sea_mask(j,k)) THEN
      land_cnt = land_cnt + 1
    END IF

    DO i=1, model_levels
      u(j,k, i) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      v(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels
      w(j,k, i-icount) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      t(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      theta(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      q(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      qcl(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      qcf(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels-1
      layer_cloud(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    IF (land_sea_mask(j,k)) THEN
      DO i=icount, icount + nsoilt_levs-1
        t_deep_soil(land_cnt, i-icount + 1) = resdump(j,k, i)
      END DO
      icount = i
    END IF

    DO i=icount, icount + model_levels+1
      p(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    DO i=icount, icount + model_levels
      rho(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = i
    IF (land_sea_mask(j,k)) THEN
      smc(land_cnt) = resdump(j,k, icount)
      icount = icount + 1
      canopy_gb(land_cnt) = resdump(j,k, icount)
      icount = icount + 1
    END IF

    snodep(j,k) = resdump(j,k, icount)

    icount = icount + 1
    tstar(j,k) = resdump(j,k, icount)

    icount = icount + 1
    zh(j,k) = resdump(j,k, icount)

    icount = icount + 1
    z0msea(j,k) = resdump(j,k, icount)

    icount = icount + 1
    DO i=icount, icount + n_cca_lev
      cca(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = icount + 1
    rccb(j,k) = resdump(j,k, icount)

    icount = icount + 1
    rcct(j,k) = resdump(j,k, icount)

    IF (land_sea_mask(j,k)) THEN
      DO i=1, nsoilm_levs
        smcl(land_cnt, i) = resdump(j,k, icount + i)
      END DO
    END IF

    icount = icount + 1
    DO i=icount, icount + n_cca_lev
      cca_dp(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = icount + 1
    DO i=icount, icount + n_cca_lev
      cca_md(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = icount + 1
    DO i=icount, icount + n_cca_lev
      cca_sh(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

    icount = icount + 1
    DO i=icount, icount + model_levels
      bl_w_var(j,k, i-icount + 1) = resdump(j,k, i)
    END DO

  END DO                     ! j
END DO                     ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE dumpinit

END MODULE dumpinit_mod
