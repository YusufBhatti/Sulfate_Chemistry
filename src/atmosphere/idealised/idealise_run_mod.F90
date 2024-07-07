! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Global data module for switches/options concerned with the idealised models.

MODULE idealise_run_mod

  ! Description:
  !   Module containing runtime logicals/options used by the idealised code.
  !
  ! Method:
  !   All switches/options which are contained in the &IDEALISED
  !   namelist in the IDEALISE control file are declared in this module.
  !   Default values have been declared where appropriate.
  !
  !   Any routine wishing to use these options may do so with the 'Use'
  !   statement.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Idealised
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.
  !
  ! Declarations:

  USE atmos_max_sizes, ONLY: model_levels_max
  USE model_domain_mod, ONLY: output_grid_stagger,                            &
                              FH_GridStagger_Endgame

  USE horiz_grid_mod, ONLY :  Nxi1L, Nxi1V, Nxi2L, Nxi2V,                     &
                              delta_xi1_H, delta_xi1_L,                       &
                              delta_xi2_H, delta_xi2_L

  USE planet_suite_mod, ONLY: tforce_number, trelax_number,                   &
                              nsteps_consv_print

  USE profiles_mod,     ONLY: num_data_max,                                    &
                              num_theta_relax_times, num_theta_relax_heights,  &
                              num_mv_relax_times, num_mv_relax_heights,        &
                              num_uv_relax_times, num_uv_relax_heights,        &
                              theta_relax_timescale, mv_relax_timescale,       &
                              uv_relax_timescale,                              &
                              theta_relax_height, theta_relax_time,            &
                              theta_relax_data,                                &
                              mv_relax_height, mv_relax_time,                  &
                              mv_relax_data,                                   &
                              uv_relax_height, uv_relax_time,                  &
                              u_relax_data, v_relax_data,                      &
                              num_theta_inc_times, num_theta_inc_heights,      &
                              theta_inc_field_type,                            &
                              num_mv_inc_times,    num_mv_inc_heights,         &
                              num_uv_inc_times,    num_uv_inc_heights,         &
                              theta_inc_time,      theta_inc_height,           &
                              mv_inc_time,         mv_inc_height,              &
                              uv_inc_time,         uv_inc_height,              &
                              theta_inc_data,      mv_inc_data,                &
                              u_inc_data,          v_inc_data,                 &
                              num_w_force_times,   num_w_force_heights,        &
                              w_force_time, w_force_height, w_force_data,      &
                              num_uv_geo_times,   num_uv_geo_heights,          &
                              uv_geo_time, uv_geo_height,                      &
                              u_geo_data, v_geo_data


  USE surface_flux_mod, ONLY: IdlSurfFluxSeaOption, IdlSurfFluxseaParams,      &
                              num_surface_flux_times, surface_flux_time,       &
                              sh_flux, lh_flux, time_varying, num_surf_max


  USE local_heat_mod,   ONLY: local_heat_option, local_heat_xoffset,           &
                              local_heat_yoffset, local_heat_amp,              &
                              local_heat_sigma, local_heat_base,               &
                              local_heat_top, local_heat_period, analytic


  USE missing_data_mod, ONLY: rmdi, imdi
  USE yomhook,  ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE
  SAVE

  ! ROSE guidelines suggest that we should not use default values for
  ! variables defined below but code will not work in places without them.
  ! This needs to be revisited.

  ! Values used to size arrays
  INTEGER,PARAMETER:: max_num_profile_data = 1
  INTEGER,PARAMETER:: max_num_force_times  = 1

  ! Non-namelist variables for the idealised model:
  LOGICAL :: zero_orography = .FALSE.     ! True if all points have zero orog.

  ! ===============================================>>

  ! Options that do not belong here (run_dyntest?):
  LOGICAL :: L_shallow = .FALSE.
  LOGICAL :: L_const_grav = .TRUE. !.FALSE.

  ! General idealised options:
  INTEGER :: tstep_plot_frequency = imdi
  INTEGER :: tstep_plot_start = imdi
  LOGICAL :: L_vert_Coriolis = .FALSE.
  LOGICAL :: L_fixed_lbcs = .FALSE.
  LOGICAL :: L_force_lbc = .FALSE.

  ! Orography specification:
  INTEGER :: grow_steps = imdi

  ! Surface Characteristics:
  REAL    :: T_surface = rmdi
  REAL    :: p_surface = rmdi
  REAL    :: roughlen_z0m = rmdi
  REAL    :: roughlen_z0h = rmdi
  LOGICAL :: L_spec_z0 = .FALSE.

  ! Fixed latitude Coriolis terms:
  REAL    :: f_plane = rmdi
  REAL    :: ff_plane = rmdi

  ! Forcing profiles (New Dynamics only):
  INTEGER :: num_tforce_times = imdi
  INTEGER :: num_qforce_times = imdi
  INTEGER :: num_uvforce_times = imdi
  INTEGER :: num_pforce_times = imdi
  REAL    :: tforce_time_interval = rmdi
  REAL    :: qforce_time_interval = rmdi
  REAL    :: uvforce_time_interval = rmdi
  REAL    :: pforce_time_interval = rmdi
  REAL    :: tforce_data(max_num_profile_data, max_num_force_times) = rmdi
  REAL    :: qforce_data(max_num_profile_data, max_num_force_times) = rmdi
  REAL    :: uforce_data(max_num_profile_data, max_num_force_times) = rmdi
  REAL    :: vforce_data(max_num_profile_data, max_num_force_times) = rmdi
  REAL    :: p_surface_data(max_num_force_times) = rmdi

  ! LBC forcing data
  REAL    :: tforce_data_modlev(model_levels_max, max_num_force_times) = rmdi
  REAL    :: qforce_data_modlev(model_levels_max, max_num_force_times) = rmdi
  REAL    :: uforce_data_modlev(model_levels_max, max_num_force_times) = rmdi
  REAL    :: vforce_data_modlev(model_levels_max, max_num_force_times) = rmdi


  LOGICAL :: L_geo_for = .FALSE.

  ! Held-Suarez (New Dynamics Only):
  INTEGER :: SuHe_fric = imdi
  REAL    :: SuHe_newtonian_timescale_ka = rmdi
  REAL    :: SuHe_newtonian_timescale_ks = rmdi
  REAL    :: SuHe_pole_equ_deltaT = rmdi
  REAL    :: SuHe_static_stab = rmdi
  REAL    :: SuHe_sigma_cutoff = rmdi
  REAL    :: base_frictional_timescale = rmdi
  
  ! Held-Suarez
  ! All the above parameters are hard coded, may want to re-integrate
  ! if want proper control of them again.
  LOGICAL :: L_HeldSuarez = .FALSE.
  LOGICAL :: L_HeldSuarez1_drag = .FALSE.

  ! 2d like simulations - Used for controlling local heating
  !                     - if .TRUE. then only varies in x direction
  !                     - In future may control bubble anomalies as well 
  !                     - indicates a "2d like" simulation.
  LOGICAL :: l_ideal_2d = .FALSE.

  !----------------------------------------------------------------------------
  ! Define namelist &IDEALISED read in from IDEALISE control file.
  !----------------------------------------------------------------------------

  NAMELIST/Idealised/                                                          &

  ! Options that do not belong here:
  L_shallow, L_const_grav,                                                    &
  
  ! How do I get rid of these?
  Nxi1L, Nxi1V, Nxi2L, Nxi2V,                                                 &
  delta_xi1_H, delta_xi1_L, delta_xi2_H, delta_xi2_L,                         &
  
  ! General idealised options:
  L_vert_Coriolis, L_fixed_lbcs, L_force_lbc,                                 &
  tstep_plot_frequency, tstep_plot_start,                                     &

  ! Surface Characteristics:
  T_surface, p_surface, L_spec_z0, roughlen_z0m, roughlen_z0h,                &

  ! Initial profiles:
  f_plane, ff_plane,                                                          &

  ! Other forcing options:
  IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                                 &
  L_geo_for,                                                                  &
  num_uv_geo_times, num_uv_geo_heights, uv_geo_height, uv_geo_time,           &
  u_geo_data, v_geo_data,                                                     &

  ! Idealised forcing
  num_theta_relax_heights, num_theta_relax_times, num_mv_relax_heights,       &
  num_mv_relax_times, num_uv_relax_heights, num_uv_relax_times,               &
  theta_relax_timescale, mv_relax_timescale, uv_relax_timescale,              &
  theta_relax_height, theta_relax_time, theta_relax_data,                     &
  mv_relax_height, mv_relax_time, mv_relax_data,                              &
  uv_relax_height, uv_relax_time, u_relax_data, v_relax_data,                 &
  num_theta_inc_times, num_theta_inc_heights, theta_inc_field_type,           &
  num_mv_inc_times, num_mv_inc_heights, num_uv_inc_times, num_uv_inc_heights, &
  theta_inc_time, theta_inc_height, mv_inc_time, mv_inc_height,               &
  uv_inc_time, uv_inc_height, theta_inc_data, mv_inc_data,                    &
  u_inc_data, v_inc_data,                                                     &
  num_w_force_times, num_w_force_heights, w_force_time, w_force_height,      &
  w_force_data,                                                               &
  num_surface_flux_times, surface_flux_time, sh_flux, lh_flux,                &

  ! 2d like simulation
  l_ideal_2d,                                                                 &

  ! local heating or 2d like
  local_heat_option, local_heat_xoffset, local_heat_yoffset, local_heat_amp,  &
  local_heat_sigma, local_heat_base, local_heat_top, local_heat_period,       &

  ! Planetary modelling options:
  tforce_number, trelax_number,                                               &
  nsteps_consv_print,                                                         &

  ! Held-Suarez (New Dynamics only):
  SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,                   &
  SuHe_pole_equ_deltaT, SuHe_static_stab,                                     &
  SuHe_sigma_cutoff, SuHe_fric, base_frictional_timescale,                    &

  ! Held-Suarez
  L_HeldSuarez, L_HeldSuarez1_drag

  !----------------------------------------------------------------------------

  !DrHook-related parameters
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IDEALISE_RUN_MOD'

CONTAINS
  SUBROUTINE check_nlist_idealised()

    USE dynamics_testing_mod,   ONLY: l_dry, problem_number
    USE problem_mod,            ONLY: idealised_problem
    USE chk_opts_mod, ONLY: chk_var, def_src
    USE bl_option_mod, ONLY: flux_bc_opt,interactive_fluxes

    USE Ereport_mod,            ONLY: Ereport
    USE errormessagelength_mod, ONLY: errormessagelength
    USE umPrintMgr
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'CHECK_NLIST_IDEALISED'
    REAL(KIND=jprb)                   :: zhook_handle

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    def_src = RoutineName

    ! Only check if running with specified fluxes
    IF (flux_bc_opt /= interactive_fluxes) THEN
      ! Surface forcing
      CALL chk_var(IdlSurfFluxSeaOption,'IdlSurfFluxSeaOption','[1:5]')

      IF (IdlSurfFluxSeaOption == time_varying) THEN
        CALL chk_var(num_surface_flux_times,'num_surface_flux_times','[1:100]')
        ! Not checking data arrays as could be a mix of values and missing data
      END IF

    END IF

    ! Only check if running an idealised problem i.e. number 3
    IF (problem_number == idealised_problem) THEN
      ! Local heating option
      CALL chk_var(local_heat_option,'local_heat_option','[0,1]')

      ! Check local heating variables if option 1 (analytic) or above chosen
      ! Limits same as in meta-data
      IF (local_heat_option >= analytic) THEN
        CALL chk_var(local_heat_xoffset,'local_heat_xoffset','[0.0:1.0]')
        IF (.NOT. l_ideal_2d) THEN  ! not checking if 2d like
          CALL chk_var(local_heat_yoffset,'local_heat_yoffset','[0.0:1.0]')
        END IF
        CALL chk_var(local_heat_amp,'local_heat_amp','[0.0:1000.0]')
        CALL chk_var(local_heat_sigma,'local_heat_sigma','[0.0:10000.0]')
        CALL chk_var(local_heat_base,'local_heat_base','[0.0:20000.0]')
        CALL chk_var(local_heat_top,'local_heat_top','[0.0:80000.0]')
        CALL chk_var(local_heat_period,'local_heat_period','[0.0:86400.0]')
      END IF

      ! Subsidence forcing - limited by array dimension
      CALL chk_var(num_w_force_times,'num_w_force_times','[0:100]')

      ! Only check subsidence forcing if it is being used.   
      IF (num_w_force_times > 0) THEN
        ! Number of values restricted by array size 
        CALL chk_var(num_w_force_heights,'num_w_force_heights','[0:100]')
        ! Check that number of heights * number of times is not greater >100
        IF (num_w_force_heights*num_w_force_times > 100) THEN
          WRITE(umMessage,'(A,I0)')                                       &
          'Array size is limited to 100, you need a branch if you want ', &
           num_w_force_heights*num_w_force_times
          CALL umPrint(umMessage,src=RoutineName)
        END IF 
        ! Not checking data arrays as can have valid missing data values if
        ! using only part of the input array.
      END IF 

      ! Time and vertical-varying geostrophic forcing 
      IF (num_uv_geo_times > 0) THEN
        ! Number of values restricted by array size 
        CALL chk_var(num_uv_geo_heights,'num_w_force_heights','[0:100]')
        ! Check that number of heights * number of times is not greater >100
        IF (num_uv_geo_heights*num_uv_geo_times > 100) THEN
          WRITE(umMessage,'(A,I0)')                                       &
          'Array size is limited to 100, you need a branch if you want ', &
           num_w_force_heights*num_w_force_times
          CALL umPrint(umMessage,src=RoutineName)
        END IF 
        ! Not checking data arrays as can have valid missing data values if
        ! using only part of the input array.
      END IF

    END IF

    def_src = ''
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE check_nlist_idealised

  SUBROUTINE print_nlist_idealised()
    USE umPrintMgr, ONLY : umPrint
    IMPLICIT NONE
    INTEGER :: i
    CHARACTER(LEN=50000) :: lineBuffer
    REAL(KIND=jprb) :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_IDEALISED'
    INTEGER :: k

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL umPrint('Contents of namelist idealised', &
        src='idealise_run_mod')

    ! Options that do not belong here:
    WRITE(lineBuffer,*)' L_shallow = ',L_shallow
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' L_const_grav = ',L_const_grav
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    
    ! How do I get rid of these?
    WRITE(lineBuffer,*)' Nxi1L = ',Nxi1L
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' Nxi1V = ',Nxi1V
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' Nxi2L = ',Nxi2L
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' Nxi2V = ',Nxi2V
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' delta_xi1_H = ',delta_xi1_H
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' delta_xi1_L = ',delta_xi1_L
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' delta_xi2_H = ',delta_xi2_H
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' delta_xi2_L = ',delta_xi2_L
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! General idealised options:
    WRITE(lineBuffer,*)' L_vert_Coriolis = ',L_vert_Coriolis
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' L_fixed_lbcs = ',L_fixed_lbcs
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' L_force_lbc = ',L_force_lbc
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' tstep_plot_frequency = ',tstep_plot_frequency
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' tstep_plot_start = ',tstep_plot_start
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! Surface Characteristics:
    WRITE(lineBuffer,*)' T_surface = ',T_surface
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' p_surface = ',p_surface
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' L_spec_z0 = ',L_spec_z0
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' roughlen_z0m = ',roughlen_z0m
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' roughlen_z0h = ',roughlen_z0h
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! Initial profiles:
    WRITE(lineBuffer,*)' f_plane = ',f_plane
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' ff_plane = ',ff_plane
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! Other forcing options:
    WRITE(lineBuffer,*)' IdlSurfFluxSeaOption = ',IdlSurfFluxSeaOption
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' IdlSurfFluxSeaParams = ',IdlSurfFluxSeaParams
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    ! Time varying surface fluxes
    IF (IdlSurfFluxSeaOption ==  time_varying) THEN
      WRITE(lineBuffer,'(A,I0)')' num_surface_flux_times = ',          &
                                  num_surface_flux_times
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      WRITE(lineBuffer,'(A)')' I    time (s)   sh (W/m2)  lh (W/m2)'
      CALL umPrint(lineBuffer,src='idealise_run_mod')

      DO i = 1,  num_surface_flux_times
        WRITE(lineBuffer,'(I4,F16.0,2F16.3)') i,surface_flux_time(i),  &
                                            sh_flux(i),lh_flux(i)
        CALL umPrint(lineBuffer,src='idealise_run_mod')
      END DO
    END IF

    WRITE(lineBuffer,*)' L_geo_for = ',L_geo_for
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! Idealised forcing profiles
    WRITE(lineBuffer,'(A,I0)')' num_theta_relax_heights = ',                   &
                                num_theta_relax_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_theta_relax_times = ',                     &
                                num_theta_relax_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_mv_relax_heights = ',                      &
                                num_mv_relax_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_mv_relax_times = ',                        &
                                num_mv_relax_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_uv_relax_heights = ',                      &
                                num_uv_relax_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_uv_relax_times = ',                        &
                                num_uv_relax_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_theta_inc_heights = ',                     &
                                num_theta_inc_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_theta_inc_times = ',                       &
                                num_theta_inc_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' theta_inc_field_type = ',                      &
                                theta_inc_field_type
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_mv_inc_heights = ',                        &
                                num_mv_inc_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_mv_inc_times = ',                          &
                                num_mv_inc_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_uv_inc_heights = ',                        &
                                num_uv_inc_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_uv_inc_times = ',                          &
                                num_uv_inc_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_w_force_times = ',                         &
                                num_w_force_times
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(A,I0)')' num_w_force_heights = ',                       &
                                num_w_force_heights
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    IF (num_w_force_times > 0) THEN
    WRITE(lineBuffer,'(A)') ' k height w_force_data'
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    DO k=1,num_w_force_heights
      WRITE(lineBuffer,'(I4,F10.1,F12.6)') k, w_force_height(k), w_force_data(k)
      CALL umPrint(lineBuffer,src='idealise_run_mod')     
    END DO
    END IF

    ! Time and vertical- varying geostrophic forcing 
    IF (num_uv_geo_times > 0) THEN
      WRITE(lineBuffer,'(A)') ' k height u_geo_data'
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      DO k=1,num_uv_geo_heights
        WRITE(lineBuffer,'(I4,F10.1,F12.6,F12.6)') k, uv_geo_height(k), &
                           u_geo_data(k), v_geo_data(k)
        CALL umPrint(lineBuffer,src='idealise_run_mod')
      END DO 
    END IF
    ! Idealised forcing data arrays not printed, as they can be very large

    ! 2d like
    WRITE(lineBuffer,'(a,l1)')' l_ideal_2d = ',l_ideal_2d
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,'(a,i12)')' local_heat_option = ',local_heat_option
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    ! Only print local heating settings if in use.
    ! Already checked values within known ranges so can use format statements.
    IF (local_heat_option > 0) THEN 
      WRITE(lineBuffer,'(a,f7.4)')' local_heat_xoffset = ',local_heat_xoffset
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      IF (.NOT. l_ideal_2d) THEN  ! local_heat_yoffset not used
        WRITE(lineBuffer,'(a,f7.4)')' local_heat_yoffset = ',local_heat_yoffset
        CALL umPrint(lineBuffer,src='idealise_run_mod')
      END IF
      WRITE(lineBuffer,'(a,f8.2)')' local_heat_amp = ',local_heat_amp
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      WRITE(lineBuffer,'(a,f10.2)')' local_heat_sigma = ',local_heat_sigma
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      WRITE(lineBuffer,'(a,f10.2)')' local_heat_base = ',local_heat_base
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      WRITE(lineBuffer,'(a,f10.2)')' local_heat_top = ',local_heat_top
      CALL umPrint(lineBuffer,src='idealise_run_mod')
      WRITE(lineBuffer,'(a,f8.0)')' local_heat_period = ',local_heat_period
      CALL umPrint(lineBuffer,src='idealise_run_mod')
    END IF

    ! Planetary modelling options:
    WRITE(lineBuffer,*)'tforce_number = ',tforce_number
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)'trelax_number = ',trelax_number
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)'nsteps_consv_print = ',nsteps_consv_print
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    ! Held-Suarez (New Dynamics only):
    WRITE(lineBuffer,*)' SuHe_newtonian_timescale_ka = ',                     &
                                                    SuHe_newtonian_timescale_ka
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' SuHe_newtonian_timescale_ks = ',                     &
                                                    SuHe_newtonian_timescale_ks
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' SuHe_pole_equ_deltaT = ',SuHe_pole_equ_deltaT
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' SuHe_static_stab = ',SuHe_static_stab
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' SuHe_sigma_cutoff = ',SuHe_sigma_cutoff
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' SuHe_fric = ',SuHe_fric
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' base_frictional_timescale = ',base_frictional_timescale
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    
    ! Held-Suarez
    WRITE(lineBuffer,*)' L_HeldSuarez = ',L_HeldSuarez
    CALL umPrint(lineBuffer,src='idealise_run_mod')
    WRITE(lineBuffer,*)' L_HeldSuarez1_drag = ',L_HeldSuarez1_drag
    CALL umPrint(lineBuffer,src='idealise_run_mod')

    CALL umPrint('- - - - - - end of namelist - - - - - -', &
        src='idealise_run_mod')

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE print_nlist_idealised

  SUBROUTINE read_nml_idealised()

    USE um_parcore, ONLY: mype

    USE check_iostat_mod, ONLY : check_iostat

    USE setup_namelist, ONLY : setup_nml_type

    USE get_env_var_mod, ONLY: get_env_var
    
    USE ereport_mod,            ONLY: ereport
    USE errormessagelength_mod, ONLY: errormessagelength
    USE umPrintMgr

    USE file_manager, ONLY: assign_file_unit, release_file_unit
    USE filenamelength_mod, ONLY: filenamelength
    USE errormessagelength_mod, ONLY: errormessagelength

    IMPLICIT NONE

    INTEGER :: my_comm
    INTEGER :: mpl_nml_type
    INTEGER :: ErrorStatus
    INTEGER :: icode
    INTEGER :: idealise_unit
    CHARACTER(LEN=errormessagelength) :: iomessage
    CHARACTER (LEN=filenamelength)    :: idealise_file
    CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'READ_NML_IDEALISED'
    CHARACTER(LEN=errormessagelength) :: Cmessage
    REAL(KIND=jprb) :: zhook_handle

    ! set number of each type of variable in my_namelist type
    INTEGER, PARAMETER :: no_of_types = 3
    INTEGER, PARAMETER :: n_int   = 30
    INTEGER, PARAMETER :: n_real  = 26 + (1 * 4) + (27 * num_data_max)  &
                                    + (3 * num_surf_max)
    INTEGER, PARAMETER :: n_log   = 10

    TYPE my_namelist
      SEQUENCE
      INTEGER :: nxi1l
      INTEGER :: nxi1v
      INTEGER :: nxi2l
      INTEGER :: nxi2v
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: idlsurffluxseaoption
      INTEGER :: suhe_fric
      INTEGER :: tforce_number
      INTEGER :: trelax_number
      INTEGER :: nsteps_consv_print
      INTEGER :: num_theta_relax_heights
      INTEGER :: num_theta_relax_times
      INTEGER :: num_mv_relax_heights
      INTEGER :: num_mv_relax_times
      INTEGER :: num_uv_relax_heights
      INTEGER :: num_uv_relax_times
      INTEGER :: num_theta_inc_heights
      INTEGER :: num_theta_inc_times
      INTEGER :: theta_inc_field_type
      INTEGER :: num_mv_inc_heights
      INTEGER :: num_mv_inc_times
      INTEGER :: num_uv_inc_heights
      INTEGER :: num_uv_inc_times
      INTEGER :: num_uv_geo_times
      INTEGER :: num_uv_geo_heights
      INTEGER :: num_w_force_times
      INTEGER :: num_w_force_heights
      INTEGER :: num_surface_flux_times
      INTEGER :: local_heat_option
      REAL    :: delta_xi1_h
      REAL    :: delta_xi1_l
      REAL    :: delta_xi2_h
      REAL    :: delta_xi2_l
      REAL    :: t_surface
      REAL    :: p_surface
      REAL    :: roughlen_z0m
      REAL    :: roughlen_z0h
      REAL    :: f_plane
      REAL    :: ff_plane
      REAL    :: idlsurffluxseaparams(4)
      REAL    :: surface_flux_time(num_surf_max)
      REAL    :: sh_flux(num_surf_max)
      REAL    :: lh_flux(num_surf_max)
      REAL    :: suhe_newtonian_timescale_ka
      REAL    :: suhe_newtonian_timescale_ks
      REAL    :: suhe_pole_equ_deltat
      REAL    :: suhe_static_stab
      REAL    :: suhe_sigma_cutoff
      REAL    :: base_frictional_timescale
      REAL    :: theta_relax_height(num_data_max)
      REAL    :: theta_relax_time(num_data_max)
      REAL    :: theta_relax_data(num_data_max)
      REAL    :: mv_relax_height(num_data_max)
      REAL    :: mv_relax_time(num_data_max)
      REAL    :: mv_relax_data(num_data_max)
      REAL    :: uv_relax_height(num_data_max)
      REAL    :: uv_relax_time(num_data_max)
      REAL    :: u_relax_data(num_data_max)
      REAL    :: v_relax_data(num_data_max)
      REAL    :: theta_relax_timescale
      REAL    :: mv_relax_timescale
      REAL    :: uv_relax_timescale
      REAL    :: theta_inc_height(num_data_max)
      REAL    :: theta_inc_time(num_data_max)
      REAL    :: theta_inc_data(num_data_max)
      REAL    :: mv_inc_height(num_data_max)
      REAL    :: mv_inc_time(num_data_max)
      REAL    :: mv_inc_data(num_data_max)
      REAL    :: uv_inc_height(num_data_max)
      REAL    :: uv_inc_time(num_data_max)
      REAL    :: u_inc_data(num_data_max)
      REAL    :: v_inc_data(num_data_max)
      REAL    :: uv_geo_height(num_data_max)
      REAL    :: uv_geo_time(num_data_max)
      REAL    :: u_geo_data(num_data_max)
      REAL    :: v_geo_data(num_data_max)
      REAL    :: w_force_height(num_data_max)
      REAL    :: w_force_time(num_data_max)
      REAL    :: w_force_data(num_data_max)
      REAL    :: local_heat_xoffset
      REAL    :: local_heat_yoffset
      REAL    :: local_heat_amp
      REAL    :: local_heat_sigma
      REAL    :: local_heat_base
      REAL    :: local_heat_top
      REAL    :: local_heat_period
      LOGICAL :: l_shallow
      LOGICAL :: l_const_grav
      LOGICAL :: l_vert_coriolis
      LOGICAL :: l_fixed_lbcs
      LOGICAL :: l_force_lbc
      LOGICAL :: l_spec_z0
      LOGICAL :: l_geo_for
      LOGICAL :: l_heldsuarez
      LOGICAL :: l_heldsuarez1_drag
      LOGICAL :: l_ideal_2d
    END TYPE my_namelist

    TYPE (my_namelist) :: my_nml

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL gc_get_communicator(my_comm, icode)

    CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                        n_real_in=n_real, n_log_in=n_log)

    IF (mype == 0) THEN

      ! Until idealised namelists are part of CNTLATM or SHARED files
      ! we need to manually open the IDEALISE namelist file

      ! Find and open the file containing Idealised model settings
      CALL get_env_var('IDEALISE', idealise_file)

      ! Open the file containing Idealised model settings
      CALL assign_file_unit(idealise_file, idealise_unit, HANDLER="fortran")
      OPEN(UNIT=idealise_unit, FILE=idealise_file, IOSTAT=ErrorStatus,       &
           IOMSG=iomessage)
      IF (errorstatus /= 0) THEN
        cmessage = ' Error opening file:'//TRIM(idealise_file) // ' :' //   &
                     newline  // TRIM(iomessage)
        ErrorStatus = 20
        CALL Ereport( routineName, errorStatus, cmessage )
      END IF

      ! Open the file containing idealised namelist
      IF ( PrintStatus >= PrStatus_Oper ) THEN
        WRITE(umMessage,'(''Idealise file: '',a)') TRIM(idealise_file)
        CALL umPrint(umMessage,src='idealise_run_mod')
      END IF

      ! Read idealised namelist
      READ( UNIT=idealise_unit, NML=Idealised, IOSTAT=ErrorStatus,            &
            IOMSG=iomessage )
      CALL check_iostat(errorstatus, "namelist Idealised", iomessage)

      CLOSE( UNIT=idealise_unit )
      CALL release_file_unit(idealise_unit, HANDLER="fortran")

      my_nml % nxi1l                          = nxi1l
      my_nml % nxi1v                          = nxi1v
      my_nml % nxi2l                          = nxi2l
      my_nml % nxi2v                          = nxi2v
      my_nml % tstep_plot_frequency           = tstep_plot_frequency
      my_nml % tstep_plot_start               = tstep_plot_start
      my_nml % idlsurffluxseaoption           = idlsurffluxseaoption
      my_nml % suhe_fric                      = suhe_fric
      my_nml % tforce_number                  = tforce_number
      my_nml % trelax_number                  = trelax_number
      my_nml % nsteps_consv_print             = nsteps_consv_print
      my_nml % num_theta_relax_heights        = num_theta_relax_heights 
      my_nml % num_theta_relax_times          = num_theta_relax_times  
      my_nml % num_mv_relax_heights           = num_mv_relax_heights    
      my_nml % num_mv_relax_times             = num_mv_relax_times     
      my_nml % num_uv_relax_heights           = num_uv_relax_heights    
      my_nml % num_uv_relax_times             = num_uv_relax_times     
      my_nml % num_theta_inc_heights          = num_theta_inc_heights 
      my_nml % num_theta_inc_times            = num_theta_inc_times  
      my_nml % theta_inc_field_type           = theta_inc_field_type  
      my_nml % num_mv_inc_heights             = num_mv_inc_heights    
      my_nml % num_mv_inc_times               = num_mv_inc_times     
      my_nml % num_uv_inc_heights             = num_uv_inc_heights    
      my_nml % num_uv_inc_times               = num_uv_inc_times     
      my_nml % num_uv_geo_heights             = num_uv_geo_heights    
      my_nml % num_uv_geo_times               = num_uv_geo_times     
      my_nml % num_w_force_heights            = num_w_force_heights    
      my_nml % num_w_force_times              = num_w_force_times     
      my_nml % num_surface_flux_times         = num_surface_flux_times     
      my_nml % local_heat_option              = local_heat_option
! end of integers
      my_nml % delta_xi1_h                    = delta_xi1_h
      my_nml % delta_xi1_l                    = delta_xi1_l
      my_nml % delta_xi2_h                    = delta_xi2_h
      my_nml % delta_xi2_l                    = delta_xi2_l
      my_nml % t_surface                      = t_surface
      my_nml % p_surface                      = p_surface
      my_nml % roughlen_z0m                   = roughlen_z0m
      my_nml % roughlen_z0h                   = roughlen_z0h
      my_nml % f_plane                        = f_plane
      my_nml % ff_plane                       = ff_plane
      my_nml % idlsurffluxseaparams           = idlsurffluxseaparams
      my_nml % surface_flux_time              = surface_flux_time
      my_nml % sh_flux                        = sh_flux
      my_nml % lh_flux                        = lh_flux
      my_nml % suhe_newtonian_timescale_ka    = suhe_newtonian_timescale_ka
      my_nml % suhe_newtonian_timescale_ks    = suhe_newtonian_timescale_ks
      my_nml % suhe_pole_equ_deltat           = suhe_pole_equ_deltat
      my_nml % suhe_static_stab               = suhe_static_stab
      my_nml % suhe_sigma_cutoff              = suhe_sigma_cutoff
      my_nml % base_frictional_timescale      = base_frictional_timescale
      my_nml % theta_relax_height             = theta_relax_height
      my_nml % theta_relax_time               = theta_relax_time
      my_nml % theta_relax_data               = theta_relax_data
      my_nml % mv_relax_height                = mv_relax_height
      my_nml % mv_relax_time                  = mv_relax_time
      my_nml % mv_relax_data                  = mv_relax_data
      my_nml % uv_relax_height                = uv_relax_height
      my_nml % uv_relax_time                  = uv_relax_time
      my_nml % u_relax_data                   = u_relax_data
      my_nml % v_relax_data                   = v_relax_data
      my_nml % theta_relax_timescale          = theta_relax_timescale
      my_nml % mv_relax_timescale             = mv_relax_timescale
      my_nml % uv_relax_timescale             = uv_relax_timescale
      my_nml % theta_inc_height               = theta_inc_height
      my_nml % theta_inc_time                 = theta_inc_time
      my_nml % theta_inc_data                 = theta_inc_data
      my_nml % mv_inc_height                  = mv_inc_height
      my_nml % mv_inc_time                    = mv_inc_time
      my_nml % mv_inc_data                    = mv_inc_data
      my_nml % uv_inc_height                  = uv_inc_height
      my_nml % uv_inc_time                    = uv_inc_time
      my_nml % u_inc_data                     = u_inc_data
      my_nml % v_inc_data                     = v_inc_data
      my_nml % uv_geo_height                  = uv_geo_height
      my_nml % uv_geo_time                    = uv_geo_time
      my_nml % u_geo_data                     = u_geo_data
      my_nml % v_geo_data                     = v_geo_data
      my_nml % w_force_height                 = w_force_height
      my_nml % w_force_time                   = w_force_time
      my_nml % w_force_data                   = w_force_data
      my_nml % local_heat_xoffset             = local_heat_xoffset  
      my_nml % local_heat_yoffset             = local_heat_yoffset  
      my_nml % local_heat_amp                 = local_heat_amp
      my_nml % local_heat_sigma               = local_heat_sigma
      my_nml % local_heat_base                = local_heat_base
      my_nml % local_heat_top                 = local_heat_top
      my_nml % local_heat_period              = local_heat_period
! end of reals
      my_nml % l_shallow                      = l_shallow
      my_nml % l_const_grav                   = l_const_grav
      my_nml % l_vert_coriolis                = l_vert_coriolis
      my_nml % l_fixed_lbcs                   = l_fixed_lbcs
      my_nml % l_force_lbc                    = l_force_lbc
      my_nml % l_spec_z0                      = l_spec_z0
      my_nml % l_geo_for                      = l_geo_for
      my_nml % l_heldsuarez                   = l_heldsuarez
      my_nml % l_heldsuarez1_drag             = l_heldsuarez1_drag
      my_nml % l_ideal_2d                     = l_ideal_2d
! end of logicals

    END IF

    CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

    IF (mype /= 0) THEN

      nxi1l                          = my_nml % nxi1l
      nxi1v                          = my_nml % nxi1v
      nxi2l                          = my_nml % nxi2l
      nxi2v                          = my_nml % nxi2v
      tstep_plot_frequency           = my_nml % tstep_plot_frequency
      tstep_plot_start               = my_nml % tstep_plot_start
      idlsurffluxseaoption           = my_nml % idlsurffluxseaoption
      suhe_fric                      = my_nml % suhe_fric
      tforce_number                  = my_nml % tforce_number
      trelax_number                  = my_nml % trelax_number
      nsteps_consv_print             = my_nml % nsteps_consv_print
      num_theta_relax_heights        = my_nml % num_theta_relax_heights 
      num_theta_relax_times          = my_nml % num_theta_relax_times  
      num_mv_relax_heights           = my_nml % num_mv_relax_heights    
      num_mv_relax_times             = my_nml % num_mv_relax_times     
      num_uv_relax_heights           = my_nml % num_uv_relax_heights    
      num_uv_relax_times             = my_nml % num_uv_relax_times     
      num_theta_inc_heights          = my_nml % num_theta_inc_heights 
      num_theta_inc_times            = my_nml % num_theta_inc_times  
      theta_inc_field_type           = my_nml % theta_inc_field_type  
      num_mv_inc_heights             = my_nml % num_mv_inc_heights    
      num_mv_inc_times               = my_nml % num_mv_inc_times     
      num_uv_inc_heights             = my_nml % num_uv_inc_heights    
      num_uv_inc_times               = my_nml % num_uv_inc_times     
      num_uv_geo_heights             = my_nml % num_uv_geo_heights    
      num_uv_geo_times               = my_nml % num_uv_geo_times     
      num_w_force_heights            = my_nml % num_w_force_heights    
      num_w_force_times              = my_nml % num_w_force_times     
      num_surface_flux_times         = my_nml % num_surface_flux_times     
      local_heat_option              = my_nml % local_heat_option
! end of integers
      delta_xi1_h                    = my_nml % delta_xi1_h
      delta_xi1_l                    = my_nml % delta_xi1_l
      delta_xi2_h                    = my_nml % delta_xi2_h
      delta_xi2_l                    = my_nml % delta_xi2_l
      t_surface                      = my_nml % t_surface
      p_surface                      = my_nml % p_surface
      roughlen_z0m                   = my_nml % roughlen_z0m
      roughlen_z0h                   = my_nml % roughlen_z0h
      f_plane                        = my_nml % f_plane
      ff_plane                       = my_nml % ff_plane
      idlsurffluxseaparams           = my_nml % idlsurffluxseaparams
      surface_flux_time              = my_nml % surface_flux_time
      sh_flux                        = my_nml % sh_flux
      lh_flux                        = my_nml % lh_flux
      suhe_newtonian_timescale_ka    = my_nml % suhe_newtonian_timescale_ka
      suhe_newtonian_timescale_ks    = my_nml % suhe_newtonian_timescale_ks
      suhe_pole_equ_deltat           = my_nml % suhe_pole_equ_deltat
      suhe_static_stab               = my_nml % suhe_static_stab
      suhe_sigma_cutoff              = my_nml % suhe_sigma_cutoff
      base_frictional_timescale      = my_nml % base_frictional_timescale
      theta_relax_height             = my_nml % theta_relax_height
      theta_relax_time               = my_nml % theta_relax_time
      theta_relax_data               = my_nml % theta_relax_data
      mv_relax_height                = my_nml % mv_relax_height
      mv_relax_time                  = my_nml % mv_relax_time
      mv_relax_data                  = my_nml % mv_relax_data
      uv_relax_height                = my_nml % uv_relax_height
      uv_relax_time                  = my_nml % uv_relax_time
      u_relax_data                   = my_nml % u_relax_data
      v_relax_data                   = my_nml % v_relax_data
      theta_relax_timescale          = my_nml % theta_relax_timescale
      mv_relax_timescale             = my_nml % mv_relax_timescale
      uv_relax_timescale             = my_nml % uv_relax_timescale
      theta_inc_height               = my_nml % theta_inc_height
      theta_inc_time                 = my_nml % theta_inc_time
      theta_inc_data                 = my_nml % theta_inc_data
      mv_inc_height                  = my_nml % mv_inc_height
      mv_inc_time                    = my_nml % mv_inc_time
      mv_inc_data                    = my_nml % mv_inc_data
      uv_inc_height                  = my_nml % uv_inc_height
      uv_inc_time                    = my_nml % uv_inc_time
      u_inc_data                     = my_nml % u_inc_data
      v_inc_data                     = my_nml % v_inc_data
      uv_geo_height                  = my_nml % uv_geo_height
      uv_geo_time                    = my_nml % uv_geo_time
      u_geo_data                     = my_nml % u_geo_data
      v_geo_data                     = my_nml % v_geo_data
      w_force_height                 = my_nml % w_force_height
      w_force_time                   = my_nml % w_force_time
      w_force_data                   = my_nml % w_force_data
      local_heat_xoffset             = my_nml % local_heat_xoffset
      local_heat_yoffset             = my_nml % local_heat_yoffset
      local_heat_amp                 = my_nml % local_heat_amp
      local_heat_sigma               = my_nml % local_heat_sigma
      local_heat_base                = my_nml % local_heat_base
      local_heat_top                 = my_nml % local_heat_top
      local_heat_period              = my_nml % local_heat_period
! end of reals
      l_shallow                      = my_nml % l_shallow
      l_const_grav                   = my_nml % l_const_grav
      l_vert_coriolis                = my_nml % l_vert_coriolis
      l_fixed_lbcs                   = my_nml % l_fixed_lbcs
      l_force_lbc                    = my_nml % l_force_lbc
      l_spec_z0                      = my_nml % l_spec_z0
      l_geo_for                      = my_nml % l_geo_for
      l_heldsuarez                   = my_nml % l_heldsuarez
      l_heldsuarez1_drag             = my_nml % l_heldsuarez1_drag
      l_ideal_2d                     = my_nml % l_ideal_2d
! end of logicals

    END IF

    CALL mpl_type_free(mpl_nml_type,icode)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE read_nml_idealised

END MODULE idealise_run_mod
