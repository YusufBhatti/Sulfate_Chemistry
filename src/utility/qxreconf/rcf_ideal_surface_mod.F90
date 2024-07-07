! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with the idealised models.

MODULE rcf_ideal_surface_mod

  ! Description:
  !   Module containing runtime logicals/options used by the idealised code.
  !
  ! Method:
  !   All switches/options which are contained in the &recon_idealised
  !   namelist in the RECONA control file are declared in this module.
  !   
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Idealised
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v10.0 programming standards.

  IMPLICIT NONE

  INTEGER, PARAMETER :: surface_zero=0
  INTEGER, PARAMETER :: surface_ellipse=1
  INTEGER, PARAMETER :: surface_ridge=2
  INTEGER, PARAMETER :: surface_plateau=3
  INTEGER, PARAMETER :: surface_massif=4
  INTEGER, PARAMETER :: surface_mask=5
  INTEGER, PARAMETER :: surface_gauss=6
  INTEGER, PARAMETER :: surface_ridge_series=7
  ! ENDGAME-only parameters
  INTEGER, PARAMETER :: surface_schar_ridge=8
  INTEGER, PARAMETER :: surface_baroclinic=9
  ! End of ENDGAME-only parameters
  INTEGER, PARAMETER :: surface_dump=10
  ! More ENDGAME-only parameters
  INTEGER, PARAMETER :: surface_wave_sin=11
  INTEGER, PARAMETER :: surface_wave_gauss=12
  INTEGER, PARAMETER :: surface_ancil=13

END MODULE rcf_ideal_surface_mod
