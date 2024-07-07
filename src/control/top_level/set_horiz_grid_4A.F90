MODULE set_horiz_grid_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_HORIZ_GRID_MOD'

CONTAINS
SUBROUTINE set_horiz_grid()

USE conversions_mod, ONLY: pi_over_180

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE horiz_grid_mod

USE atm_fields_bounds_mod

USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows
USE um_parvars, ONLY: halo_i, halo_j, datastart

USE model_domain_mod, ONLY: model_type, mt_global, mt_cyclic_lam

IMPLICIT NONE
!
! Description: Computes the horizontal ENDGame grid.
!
!
! Method: Two separate cases have been coded: (i) even number of points and
!         (ii) odd number points across the domain. In (ii) the grid is
!         construction ensures symmetry.
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Local variables
INTEGER :: i, j, i_mid, j_mid, gi, gj


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_HORIZ_GRID'

INTEGER :: alloc_stat,ErrorStatus


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Adjust grid if necessary
IF ( .NOT. cartesian_grid .AND.                                      &
             (model_type == mt_global .OR.                           &
              model_type == mt_cyclic_lam) ) THEN

  base_xi1  = 0.0
  delta_xi1 = 360.0/global_row_length

  IF (model_type == mt_global) THEN
    delta_xi2 = 180.0/global_rows
    base_xi2 = -90.0
  END IF
END IF

! Set grid assuming a regular grid

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*) '** Generating standard E-W grid **'
  CALL umPrint(umMessage,src='set_horiz_grid_4A')
END IF

DO i=pdims_l%i_start,pdims_l%i_end
  gi = datastart(1) + i - 1
  xi1_p(i) = base_xi1 + (gi-0.5)*delta_xi1
END DO

DO i=udims_l%i_start,udims_l%i_end
  gi = datastart(1) + i - 1
  xi1_u(i) = base_xi1 + gi*delta_xi1
END DO

IF (PrintStatus > PrStatus_Normal) THEN
  WRITE(umMessage,*) '** Generating standard N-S grid **'
  CALL umPrint(umMessage,src='set_horiz_grid_4A')
END IF
DO j=pdims_l%j_start,pdims_l%j_end
  gj = datastart(2) + j - 1
  xi2_p(j) = base_xi2 + (gj-0.5)*delta_xi2
END DO

DO j=vdims_l%j_start,vdims_l%j_end
  gj = datastart(2) + j - 1
  xi2_v(j) = base_xi2 + gj*delta_xi2
END DO

! Global arrays used by SLICE

DO i = 1-halo_i, global_row_length+halo_i
  glob_xi1_p(i) = base_xi1 + (i-0.5)*delta_xi1
END DO
DO i = -halo_i, global_row_length+halo_i
  glob_xi1_u(i) = base_xi1 + i*delta_xi1
END DO

DO j = 1-halo_j, global_rows+halo_j
  glob_xi2_p(j) = base_xi2 + (j-0.5)*delta_xi2
END DO

DO j = -halo_j, global_rows+halo_j
  glob_xi2_v(j) = base_xi2 + j*delta_xi2
END DO

! If not a cartesian grid convert to Radians

IF ( .NOT. cartesian_grid ) THEN

  xi1_u = pi_over_180*xi1_u
  xi1_p = pi_over_180*xi1_p

  xi2_v = pi_over_180*xi2_v
  xi2_p = pi_over_180*xi2_p

  glob_xi1_p = pi_over_180*glob_xi1_p
  glob_xi1_u = pi_over_180*glob_xi1_u
  glob_xi2_p = pi_over_180*glob_xi2_p
  glob_xi2_v = pi_over_180*glob_xi2_v

  ! Reset grid descriptors (degrees -> radians)

  base_xi1  = pi_over_180*base_xi1
  base_xi2  = pi_over_180*base_xi2
  delta_xi1 = pi_over_180*delta_xi1
  delta_xi2 = pi_over_180*delta_xi2

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_horiz_grid
END MODULE set_horiz_grid_mod
