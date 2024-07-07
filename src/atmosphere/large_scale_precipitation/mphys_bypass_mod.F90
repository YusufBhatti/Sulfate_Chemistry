! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Global bypass module for switches concerned with microphysics

MODULE mphys_bypass_mod

IMPLICIT NONE
SAVE

!  Description:
!   Module containing logical switches used by the microphysics code
!   that are not part of the run_precip or run_cloud namelists

!  Method:
!   Variables declared here are initialised in atm_step and then
!   used in the microphysics scheme.

!   This module is only for switches not in the run_cloud or
!   run_precip namelists. New microphysics or cloud scheme
!   logicals should be put in the run_precip or run_cloud namelist
!   and not here.

!   In theory, this module can be made redundant if all the switches
!   below are put into other modules.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large-scale preciptation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

!-----------------------------------------------------
! 1. Microphysics Switches not part of run_precip
!-----------------------------------------------------

LOGICAL :: l_crystals     = .FALSE.
! Controls whether a crystal mode is turned on or not.

LOGICAL :: l_ref_diag     = .FALSE.
! Controls whether the subroutine mphys_reflec is called
! or not (called if .TRUE.). If no reflectivity diagnostics
! are required, this should save time by skipping this bit
! of code.

LOGICAL :: l_last_iter    = .FALSE.
! Useful for diagnostics called only once per timestep.
! This is .TRUE. when the number of column iterations is exactly
! 1 or when the iteration is the last before the end of the timestep.
! It saves having code run unnecessarily when multiple iterations
! are in operation.

!-----------------------------------------------------
! 2. Top of model;
! required in microphysics
!-----------------------------------------------------

REAL :: mphys_mod_top

!-----------------------------------------------------
! 3. Dimensions of qcf2; required for electric scheme
!-----------------------------------------------------

INTEGER :: qcf2_idims_start
INTEGER :: qcf2_idims_end
INTEGER :: qcf2_jdims_start
INTEGER :: qcf2_jdims_end

INTEGER, PARAMETER :: qcf2_kdims_start = 1

INTEGER :: qcf2_kdims_end

END MODULE mphys_bypass_mod
