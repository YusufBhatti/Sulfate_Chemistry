! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE ukca_eg_tracers_total_mass_mod
IMPLICIT NONE

! Description:
!             This routine computes the total mass of species qs,
!             massqs = sum[rho*av(qs)*volume], vertically integrated.
!             Note that qs are averaged to rho-levels.
!
! Note:
! UKCA version of routine from tracer_advection to return vertically integrated
! tracer mass, and using rho_r2 to calculate rho.
!
! Method: ENDGame formulation version 3.02
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'UKCA_EG_TRACERS_TOTAL_MASS_MOD'

CONTAINS

SUBROUTINE ukca_eg_tracers_total_mass_fix(rho_r2, qs, number_qs, mass_2d)

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE level_heights_mod,     ONLY: r_rho_levels
USE eg_helmholtz_mod,      ONLY: ec_vol
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: number_qs     ! No. of fields

! air density * radius^2
REAL, INTENT(IN) :: rho_r2(row_length, rows, model_levels)

! tracer fields
REAL, INTENT(IN) :: qs(1:row_length, 1:rows, 0:model_levels, 1:number_qs)

! vertically integrated tracer mass
REAL, INTENT(OUT) :: mass_2d(row_length, rows, number_qs)

! Local variables
! ===============

! minimum qs for calculation
REAL, PARAMETER :: qsmin = 1.0E-30

! air density on rho_levels
REAL :: rho(row_length, rows, model_levels)

! interpolation weights
REAL :: alfa_za(model_levels)
REAL :: beta_za(model_levels)

! loop counters
INTEGER  :: i,j,k,kk

! partial sum
REAL     :: t1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EG_TRACERS_TOTAL_MASS_FIX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mass_2d(:,:,:) = 0.0

IF (MAXVAL(qs) > qsmin) THEN
! compute column integrated mass of each species in qs array

  DO k = 1, model_levels
    rho(:,:,k) = rho_r2(:,:,k)/(r_rho_levels(1:row_length,1:rows,k) *  &
                                r_rho_levels(1:row_length,1:rows,k))
  END DO

! compute vertical averaging weights
  DO k = pdims%k_start+1, pdims%k_end
    alfa_za(k) = (   eta_rho_levels(k) - eta_theta_levels(k-1) )/ &
                 ( eta_theta_levels(k) - eta_theta_levels(k-1) )
    beta_za(k) = 1.0 - alfa_za(k)
  END DO
  alfa_za(pdims%k_start) = 1.0
  beta_za(pdims%k_start) = 0.0

  DO kk = 1, number_qs
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          t1 =  alfa_za(k)*  qs(i,j,k,kk) +  &
                beta_za(k)*  qs(i,j,k-1,kk)

          mass_2d(i,j,kk) = mass_2d(i,j,kk) + t1*rho(i,j,k)*ec_vol(i,j,k)
        END DO
      END DO
    END DO
  END DO

END IF   ! ANY qs > qsmin

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_eg_tracers_total_mass_fix

END MODULE ukca_eg_tracers_total_mass_mod
