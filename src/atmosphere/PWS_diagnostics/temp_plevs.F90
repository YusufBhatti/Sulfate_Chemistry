! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE temp_plevs_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TEMP_PLEVS_MOD'

CONTAINS

SUBROUTINE temperature_plevs(NumPlevs, plevs, temp_plevs)

! Description: Routine to calculate temperature on p levels
!
! Method: vertical interpolation followed by horiz interp, 
!         using standard utility routines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod
USE atm_fields_real_mod, ONLY: exner_theta_levels, theta, p_theta_levels

USE vert_interp2_mod, ONLY: vert_interp2
USE t_vert_interp_to_p_mod, ONLY: t_vert_interp_to_p

USE planet_constants_mod, ONLY: kappa, p_zero
USE UM_ParVars          , ONLY: offx, offy, halo_i, halo_j
USE interpor_mod        , ONLY: interp_order_linear
USE nlsizes_namelist_mod, ONLY: model_levels, row_length, rows, n_rows, &
                                global_row_length, bl_levels

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumPLevs  ! Number of P levels for output

REAL, INTENT(IN) :: plevs(NumPLevs) 

! Temperature p levs
REAL, INTENT(OUT) :: temp_plevs(tdims%i_start:tdims%i_end, &
                                tdims%j_start:tdims%j_end, NumPLevs) 


! Local variables:
INTEGER :: i, j, k

REAL :: pressure_pa ! Pressure in Pascal
REAL :: pressure_ex ! Exner pressure

! Temperature on model levs
REAL :: temp_mlevs(tdims%i_start:tdims%i_end, &
                   tdims%j_start:tdims%j_end, &
                               1:model_levels) 

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'TEMPERATURE_PLEVS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate temperature at theta points, model levels

DO k = 1, model_levels
  DO j = tdims%j_start, tdims%j_end
    DO i =   tdims%i_start , tdims%i_end
      temp_mlevs(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k

DO  k = 1, NumPLevs

  pressure_pa = PLevs(k)*100.0 ! Convert from hPa to Pa
  pressure_ex = ( pressure_pa /p_zero )**kappa
  
  ! Interpolate temperature to the required p level
  CALL  t_vert_interp_to_p(temp_mlevs, theta, row_length, rows                &
                    ,model_levels, pressure_pa, offx, offy, halo_i, halo_j    &
                    ,p_theta_levels, bl_levels, exner_theta_levels            &
                    ,temp_plevs(1,1,k) )

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE temperature_plevs

END MODULE temp_plevs_mod
