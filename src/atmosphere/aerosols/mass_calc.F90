! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE MASS_CALC
!
! Purpose:
!   To calculate mass of air in a model layer
!
!   Called by mcr_ctl
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran 90.
!  This code is written to UMDP3 v8 programming standards
!
! Documentation: UMDP20
!----------------------------------------------------------------------
MODULE mass_calc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MASS_CALC_MOD'

CONTAINS

SUBROUTINE mass_calc(row_length, rows, r_rho_levels, r_theta_levels,           &
                     timestep, rho_r2, q, qcl, qcf, dm )

! calculates mass of (dry) air per square metre

USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with intent IN:

INTEGER :: row_length
INTEGER :: rows
INTEGER :: ccldbase(row_length,rows)      !convective cloud base
INTEGER :: ccldtop(row_length,rows)       !convective cloud top

REAL :: r_rho_levels(row_length,rows,model_levels)
REAL :: r_theta_levels(row_length,rows,0:model_levels)
REAL :: rho_r2(row_length,rows,model_levels)              ! density*r*r
REAL :: q(row_length,rows,model_levels)                   ! water vapour(kg/kg)
REAL :: qcl(row_length,rows,model_levels)
REAL :: qcf(row_length,rows,model_levels)

REAL ::  timestep                                         ! timestep in secs

! Arguments with intent OUT :
REAL :: dm(row_length, rows, model_levels)                ! mass of air

! Local variables

INTEGER :: i, j, k                                        ! Loop variables

REAL    :: rho1, rho2                                     ! air densities

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MASS_CALC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( model_levels, rows, row_length, rho_r2, r_rho_levels,    &
!$OMP         dm, r_theta_levels, q, qcl, qcf )                        &
!$OMP PRIVATE( i, j, k, rho1, rho2 )
!$OMP DO SCHEDULE(STATIC)
DO k = 1,model_levels-1
  DO j = 1,rows
    DO i = 1,row_length
      ! Remove the r squared factor from rho before interpolation
      rho1 = rho_r2(i,j,k) /                                                   &
        ( r_rho_levels(i,j,k) * r_rho_levels(i,j,k) )
      rho2 = rho_r2(i,j,k+1)/                                                  &
        ( r_rho_levels(i,j,k+1) *  r_rho_levels(i,j,k+1) )
      ! DM = density (interpolated on to theta levels) * delta r
      dm(i,j,k) = rho2 * ( r_theta_levels(i,j,k) -                             &
        r_rho_levels(i,j,k) ) +                                                &
        rho1 * ( r_rho_levels(i,j,k+1) -                                       &
        r_theta_levels(i,j,k) )

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !MODEL_LEVELS
!$OMP END DO

! Special case for lowest layer to get correct mass
!$OMP DO SCHEDULE(STATIC)
DO j = 1,rows
  DO i = 1,row_length
    dm(i,j,1) =                                                                &
      dm(i,j,1) * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0)) /              &
      (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
  END DO !ROW_LENGTH
END DO !ROWS
!$OMP END DO

! Convert DM to DRY density if level is wet
!$OMP DO SCHEDULE(STATIC)
DO k = 1,model_levels-1
  DO j = 1,rows
    DO i = 1,row_length
      dm(i,j,k) = dm (i,j,k) /                                                 &
        (1.0 + q(i,j,k) + qcl(i,j,k) + qcf(i,j,k))
    END DO !ROW_LENGTH
  END DO !ROWS
END DO !MODEL_LEVELS
!$OMP END DO NOWAIT

! Top level
!$OMP DO SCHEDULE(STATIC)
DO j = 1,rows
  DO i = 1,row_length
    dm(i,j,model_levels) = 0.0
  END DO !ROW_LENGTH
END DO !ROWS
!$OMP END DO
!$OMP END PARALLEL


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mass_calc
END MODULE mass_calc_mod
