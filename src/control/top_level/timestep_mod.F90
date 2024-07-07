! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Timestepping values.

MODULE timestep_mod

! Description:
!              Timestep information
!              Updated at start of ATM_STEP every timestep

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

INTEGER    :: timestep_number

REAL       :: timestep             ! atmosphere model timestep
REAL       :: recip_timestep       ! recip model timestep

! These need the SAVE attribute, because they are set from a subroutine
! call outside of atm_step; timestep_mod then goes out of scope before
! timestepping begins.
REAL, SAVE :: radiation_tstep_diag ! timestep of fast radiation (3C)
REAL, SAVE :: radiation_tstep_prog ! timestep of slow radiation (3C)

REAL       :: chemistry_timestep   ! must be  <=  model timestep
REAL       :: pos_timestep         ! = +timestep.
REAL       :: neg_timestep         ! = -timestep.

END MODULE timestep_mod
