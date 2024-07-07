! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Semi-Lagrangian parameters.
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics

MODULE diag_ctl_mod

  ! Description:
  !          Diagnostic printing control parameters
  !
  ! Method:
  !         The required parameters are read in by NAMELIST
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v6 programming standards.

  IMPLICIT NONE

  INTEGER :: rpemax ! array size needed for diagnostic printing
  INTEGER :: rpemin ! array size needed for diagnostic printing
  INTEGER :: rpesum ! array size needed for diagnostic printing
  INTEGER :: ipesum ! array size needed for diagnostic printing

  INTEGER :: time_theta1_min  ! Timestep of min level 1 theta
  INTEGER, ALLOCATABLE :: time_w_max(:) ! Timestep of max w
  INTEGER, ALLOCATABLE :: time_div_max(:) ! Timestep of max div
  INTEGER, ALLOCATABLE :: time_div_min(:) ! Timestep of min div
  INTEGER, ALLOCATABLE :: time_lapse_min(:) ! Timestep of min
  INTEGER, ALLOCATABLE :: time_max_shear(:) !Timestep max shear
  INTEGER, ALLOCATABLE :: time_max_wind(:) ! Timestep of max wind
  INTEGER, ALLOCATABLE :: time_KE_max(:) ! Timestep of max KE
  INTEGER, ALLOCATABLE :: time_KE_min(:) ! Timestep of min KE
  INTEGER, ALLOCATABLE :: time_noise_max(:) ! Timestep of max

  REAL :: min_theta1_run  ! Min theta level 1
  REAL :: dtheta1_run     ! Largest -ve delta theta at min theta1
  REAL, ALLOCATABLE :: max_w_run(:)  ! Max w at a level
  REAL, ALLOCATABLE :: max_div_run(:) ! Max divergence at a level
  REAL, ALLOCATABLE :: min_div_run(:) ! Min divergence at a level
  REAL, ALLOCATABLE :: min_lapse_run(:) ! Min dtheta/dz at a level
  REAL, ALLOCATABLE :: max_shear_run(:) ! Max shear at a level
  REAL, ALLOCATABLE :: max_wind_run(:) ! Max wind at a level
  REAL, ALLOCATABLE :: max_KE_run(:)   ! Max KE at a level
  REAL, ALLOCATABLE :: min_KE_run(:)   ! Min KE at a level
  REAL, ALLOCATABLE :: max_noise_run(:) ! Max noise at a level

CONTAINS
  SUBROUTINE init_print_diag()
    USE turb_diff_mod, ONLY:                                                  &
       diag_interval, print_step, L_diag_L2norms, L_diag_L2helm,              &
       L_print_theta1, L_print_lapse, L_print_div, L_print_w, L_print_wmax,   &
       L_print_shear, L_print_max_wind, L_diag_wind, L_print_pe, L_flush6,    &
       norm_lev_start, norm_lev_end
    USE diag_print_mod, ONLY:                                                 &
       L_print_allpe, L_flush, norm_start_lev, norm_end_lev
    USE nlsizes_namelist_mod, ONLY: model_levels
    USE umPrintMgr, ONLY : umPrint, umMessage

    IMPLICIT NONE

    IF ( L_diag_L2norms ) THEN
      WRITE(umMessage,FMT='(A)') 'Printing of norms during timestep activated'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF ! L_diag_L2norms
    IF ( L_diag_L2helm ) THEN
      WRITE(umMessage,FMT='(A)') 'Printing of coefficient norms in solver'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF ! L_diag_L2helm

    ! Sampling interval must not be larger than print interval
    IF ( diag_interval > print_step) THEN
      diag_interval = print_step
    END IF !  diag_interval > print_step

    ! Set some New Dynamics only varaibles
    L_print_allpe  = L_print_pe
    L_flush        = L_flush6
    norm_start_lev = norm_lev_start
    norm_end_lev   = norm_lev_end

    ! Set array sizes needed for chosen diagnostic printing
    ! These are needed to do sums and max/mins across processors
    rpemax = 0
    rpemin = 0
    ipesum = 0
    rpesum = 0
    IF (L_print_theta1 ) THEN
      rpemin = rpemin + 1
      ipesum = ipesum + 1
      rpesum = rpesum + 3
    END IF !  L_print_theta1
    IF (L_print_lapse ) THEN
      rpemin = rpemin + (model_levels - 1)
      ipesum = ipesum + 2 * (model_levels - 1)
      rpesum = rpesum + 2 * (model_levels - 1)
    END IF !  L_print_lapse
    IF (L_print_div) THEN
      rpemax = rpemax + model_levels
      rpemin = rpemin + model_levels
      ipesum = ipesum + 2 * model_levels
      rpesum = rpesum + 4 * model_levels
    END IF !  L_print_div
    IF (L_print_w .OR. L_print_wmax) THEN
      rpemax = rpemax + model_levels
      ipesum = ipesum + 2 * (model_levels)
      rpesum = rpesum + 2 * (model_levels)
    END IF !  L_print_w .or. L_print_wmax
    IF (L_print_shear) THEN
      rpemax = rpemax + model_levels - 1
      ipesum = ipesum + model_levels - 1
      rpesum = rpesum + 2 * (model_levels - 1)
    END IF !  L_print_shear
    IF (L_print_max_wind) THEN
      rpemax = rpemax + model_levels
      ipesum = ipesum + model_levels
      rpesum = rpesum + 2 * model_levels
    END IF !  L_print_max_wind
    IF ( L_diag_wind ) THEN
      rpesum = rpesum + 2 * model_levels
    END IF  ! L_diag_wind
    min_theta1_run = 1000.0
    time_theta1_min = 0

    ! With ENDGame need surface level for max_w, even though always zero!
    ALLOCATE(max_w_run(0:model_levels))

    ALLOCATE(time_w_max(model_levels),                                     &
      max_div_run(model_levels), time_div_max(model_levels),               &
      min_div_run(model_levels), time_div_min(model_levels),               &
      min_lapse_run(model_levels), time_lapse_min(model_levels),           &
      max_shear_run(model_levels), time_max_shear(model_levels),           &
      max_wind_run(model_levels), time_max_wind(model_levels),             &
      max_KE_run(model_levels+1), min_KE_run(model_levels+1),              &
      time_KE_max(model_levels+1), time_KE_min(model_levels+1),            &
      time_noise_max(model_levels), max_noise_run(model_levels) )

    max_w_run      = 0.0
    time_w_max     = 0
    max_div_run    = 0.0
    time_div_max   = 0
    min_div_run    = 0.0
    time_div_min   = 0
    min_lapse_run  = 1.0e6
    time_lapse_min = 0
    max_shear_run  = 0.0
    time_max_shear = 0
    max_wind_run   = 0.0
    time_max_wind  = 0
    max_KE_run     = 0.0
    min_KE_run     = 1.0e30
    time_KE_max    = 0
    time_KE_min    = 0
    time_noise_max = 0
    max_noise_run  = 0.0

  END SUBROUTINE init_print_diag

END MODULE diag_ctl_mod
