! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SUB_DATA
! Purpose:-           To be used in test runs when detailed
!                     sub-timestep diagnostics are required after
!                     calls to each subroutine or when budget calcs.
!                     are required and so called at start and
!                     end of meaning period
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model


SUBROUTINE sub_data                                                        &
  ! In leading dimensions of arrays
  ( row_length, rows, nbl_levs, nsoilt_levs                                &
  , nsoilm_levs, title1, istep, ayear, aday, u, v                          &
  , theta, q, qcl, qcf, lca, p,rho, exner_rh, exner_th, t_deep_soil        &
  , smc,canopy, snodep, tstar, zh, z0msea, cca, iccb, icct, smcl )

USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------

INTEGER ::                                                        &
  row_length             &! In x dimension.
, rows                   &! In y dimension.
, nbl_levs               &! In Number of Boundary layer levels
, nsoilt_levs            &! In Number of soil temperature levels
, nsoilm_levs             ! In Number of soil moisture levels


CHARACTER(LEN=35) ::          &
  title1                  ! Table heading

INTEGER ::               &
  iccb(row_length,rows)  &! Convective cloud base
, icct(row_length,rows)  &! Convective cloud top
, aday                   &! Actual day number
, ayear                  &! Actual year number
, istep                   ! Timestep number

REAL ::                      &
  canopy(row_length,rows)    &! Surface/canopy water (kg/m^2)
, cca(row_length,rows)       &! Convective cloud amount
, lca(row_length,rows,model_levels)   ! Layer cloud amount (decimal fraction)

REAL ::                                &
  p(row_length,rows,model_levels)      &! pressure on rho levels (Pa)
, exner_rh(row_length,rows,model_levels) &! exner pressure on rho levels
, exner_th(row_length,rows,model_levels) &! exner pressure on theta levels
, rho(row_length,rows,model_levels)    &! density
, q(row_length,rows,model_levels)      &! Specific humidity (kg/kg)
, qcf(row_length,rows,model_levels)    &! Specific cloud ice (kg/kg)
, qcl(row_length,rows,model_levels)    &! Specific cloud water (kg/kg)
, smc(row_length,rows)                 &! Soil moisture content (kg/m^2)
, smcl(row_length,rows,nsoilm_levs)    &! Soil moisture in layers (kg/m^2)
, snodep(row_length,rows)              &! Snow depth (kg/m^2)
, t(row_length,rows,model_levels)      &! Temperature (K)
, theta(row_length,rows,model_levels)  &! Potential temperature (K)
, tstar(row_length,rows)               &! Surface temp.(K)
, u(row_length,rows,model_levels)      &! Zonal wind (m/s)
, v(row_length,rows,model_levels)      &! Meridional wind (m/s)
, zh(row_length,rows)                  &! Boundary layer depth (m)
, Z0mSea(row_length,rows)               ! Sea surface roughness length (m)

REAL ::                                    &
  t_deep_soil(row_length,rows,nsoilt_levs)  ! Soil layer temps (K)


RETURN

END SUBROUTINE sub_data
