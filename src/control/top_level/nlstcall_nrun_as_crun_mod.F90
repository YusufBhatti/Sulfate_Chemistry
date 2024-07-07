! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top level control

!  Description:

! Fortran module to provide l_nrun_as_crun variable. This is in a separate
! file to allow the sstpert library to pick up a default value.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_nrun_as_crun_mod

IMPLICIT NONE

PRIVATE

! Logical controlling whether to treat an NRUN as a CRUN, bypassing the
! timestep one logic for a number of variables:
! run_stochastic - stphseed=2
! run_cloud - l_pc2_check_init=.false.
! run_dyn - alpha_relax_type=1
LOGICAL, PUBLIC :: l_nrun_as_crun = .FALSE.


END MODULE nlstcall_nrun_as_crun_mod
