! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Control variables for idealised dry planet suite

! Description:
!   
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised

! Method:
!   
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE planet_suite_mod

USE tforce_mod, ONLY: tf_none
USE trelax_mod, ONLY: tr_none

IMPLICIT NONE

! Items read in through IDEALISED namelist
INTEGER :: tforce_number = tf_none    ! Choice of forcing profile
INTEGER :: trelax_number = tr_none    ! Choice of relaxation timescale
INTEGER :: nsteps_consv_print = 0     ! Frequency of printing of AAM and KE

END MODULE planet_suite_mod
