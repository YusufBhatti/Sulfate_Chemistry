! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE curl_at_poles_mod

USE global_2d_sums_mod, ONLY: global_2d_sums

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CURL_AT_POLES_MOD'

CONTAINS
SUBROUTINE curl_at_poles(pole,u,curl_xi3_pole)

USE horiz_grid_mod,        ONLY: xi1_u, xi2_v
USE level_heights_mod,     ONLY: xi3_at_rho => r_rho_levels
USE atm_fields_bounds_mod, ONLY: udims, udims_s
USE nlsizes_namelist_mod,  ONLY: global_row_length
USE um_parvars,            ONLY: gc_proc_row_group
USE conversions_mod,       ONLY: pi
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Description:  Calcualte r-component of curl at either north or south pole.
!              !!!! Only works for zero orography !!!!
!
! Method: Chapter 8, ENDGame formulation version 1.01
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: DYNAMICS
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine arguments
INTEGER, INTENT(IN) :: pole ! 1 = North, -1 = South
REAL, INTENT(IN)    :: u(udims_s%i_start:udims_s%i_end,                        &
                         udims_s%j_start:udims_s%j_end,                        &
                         udims_s%k_start:udims_s%k_end)

REAL, INTENT(OUT)   :: curl_xi3_pole(udims%k_start:udims%k_end)

INTEGER, PARAMETER  :: NP=1, SP=-1
INTEGER             :: i, j_polar, k
REAL                :: delta_lambda, delta_phi, coeff
REAL                :: upole(udims%i_start:udims%i_end,                        &
                             udims%k_start:udims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CURL_AT_POLES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(pole)
CASE(NP)
  j_polar=udims%j_end
CASE(SP)
  j_polar=udims%j_start
END SELECT

DO k=udims%k_start,udims%k_end
  DO i=udims%i_start,udims%i_end
    upole(i,k)=u(i,j_polar,k)
  END DO
END DO

CALL global_2d_sums(upole,udims%i_end-udims%i_start+1,                         &
                    1, 0, 0, udims%k_end-udims%k_start+1,                      &
                    curl_xi3_pole,  gc_proc_row_group)

delta_lambda=xi1_u(1)-xi1_u(0)
delta_phi=xi2_v(1)-xi2_v(0)
coeff=pole*0.5*delta_lambda/(pi*TAN(0.25*delta_phi))

DO k=udims%k_start,udims%k_end
  curl_xi3_pole(k)=coeff*curl_xi3_pole(k)/xi3_at_rho(1,j_polar,k)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE curl_at_poles
END MODULE curl_at_poles_mod
