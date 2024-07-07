! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sructure containing diagnostics for the SPT scheme
!  This permits easier addition diagnostics without additional
!  passing of arguments through code tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.
!  Module is used atm_step and stp_main


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE spt_diag_mod

IMPLICIT NONE

TYPE strsptdiag

  ! Need to create a flag and a pointer

  LOGICAL ::  l_spt_forcing_pattern

  LOGICAL ::  l_spt_theta_inc
  LOGICAL ::  l_spt_q_inc
  LOGICAL ::  l_spt_u_inc
  LOGICAL ::  l_spt_v_inc
  LOGICAL ::  l_cfl_br_marker
  LOGICAL ::  l_spt_t_inc

  REAL, POINTER :: spt_forcing_pattern(:, :,:)
  !         23     SPT forcing pattern field
  REAL, POINTER :: spt_theta_inc(:, :, :)
  !         24     Theta SPT increment
  REAL, POINTER :: spt_q_inc(:, :, :)
  !         25     q SPT increment
  REAL, POINTER :: spt_u_inc(:, :, :)
  !         26     u SPT increment
  REAL, POINTER :: spt_v_inc(:, :, :)
  !         27     v SPT increment
  REAL, POINTER :: cfl_br_marker(:, :)
  !         28     Array that indicates where CFL is breached
  REAL, POINTER :: spt_t_inc(:, :, :)
  !         29     T SPT increment


END TYPE strsptdiag

END MODULE spt_diag_mod
