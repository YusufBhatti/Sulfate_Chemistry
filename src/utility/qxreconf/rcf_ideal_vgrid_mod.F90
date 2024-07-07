! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with the idealised models.

MODULE rcf_ideal_vgrid_mod

  ! Description:
  !   Module containing runtime logicals/options used by the idealised code.
  !
  ! Method:
  !   All switches/options which are contained in the &recon_idealised
  !   namelist in the RECONA control file are declared in this module.
  !   Default values have been declared where appropriate.
  !   
  !   Any routine wishing to use these options may do so with the 'USE' 
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Idealised
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v10.0 programming standards.
  !
  ! Declarations:

  IMPLICIT NONE

  INTEGER, PARAMETER :: vert_regular=1
  INTEGER, PARAMETER :: vert_quadratic_theta=21
  INTEGER, PARAMETER :: vert_bi_quadratic=22
  INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
  INTEGER, PARAMETER :: vert_schar=3
  INTEGER, PARAMETER :: vert_dwd=4
  INTEGER, PARAMETER :: vert_stretch_plus_regular=5
  INTEGER, PARAMETER :: vert_quad_stretch_thin=6
  INTEGER, PARAMETER :: vert_regular_thin=7
  INTEGER, PARAMETER :: vert_geometric_theta=8
  INTEGER, PARAMETER :: vert_dump=10
  INTEGER, PARAMETER :: vert_idl_um_grid=11

END MODULE rcf_ideal_vgrid_mod
