! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables for turbulent diffusion input
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE turb_diff_mod

! Subgrid turbulence scheme input.

USE umPrintMgr, ONLY: umprint, newline
USE atmos_max_sizes, ONLY: model_levels_max
USE Control_Max_Sizes, ONLY: max_121_rows, max_sponge_width,           &
  max_updiff_levels
USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

LOGICAL :: L_diff_active
LOGICAL :: L_subfilter_horiz = .FALSE. ! activate horiz diffusion
LOGICAL :: L_subfilter_vert  = .FALSE. ! activate vertical diffusion
LOGICAL :: L_blend_isotropic = .TRUE.
                  ! Use same diffusion coefficient in all directions

INTEGER :: turb_startlev_horiz = imdi  ! 1st lev for horiz subgrid turb
INTEGER :: turb_endlev_horiz   = imdi  ! last lev for horiz subgrid turb
INTEGER :: turb_startlev_vert  = imdi  ! 1st lev for vert subgrid turb
INTEGER :: turb_endlev_vert    = imdi  ! last lev for vert subgrid turb
INTEGER :: hyd_mix_opt         = imdi  ! mixing option for hydrometeors
INTEGER, PARAMETER :: liq_only = 1    ! only mix qcl
INTEGER, PARAMETER :: liq_and_ice = 2 ! mix qcl and qcf
INTEGER, PARAMETER :: all_hydro = 3   ! mix qcl, qcf, qr, qgraup, qcf2

REAL :: diff_factor = rmdi
REAL :: mix_factor  = rmdi

LOGICAL :: L_print_w        = .FALSE.
LOGICAL :: L_print_div      = .FALSE.
LOGICAL :: L_diag_print_ops = .FALSE. ! diagnostic prints for ops
LOGICAL :: L_print_pe       = .FALSE. ! print diagnostics on all pe's
LOGICAL :: L_print_shear    = .FALSE. ! wind diagnostic prints
LOGICAL :: L_print_max_wind = .FALSE. ! wind diagnostic prints
LOGICAL :: L_diag_L2norms   = .FALSE. ! l2norm diagnostic prints
LOGICAL :: L_diag_L2helm    = .FALSE. ! l2norm diagnostic prints from solver
LOGICAL :: L_diag_wind      = .FALSE. ! wind diagnostic prints
LOGICAL :: L_diag_noise     = .FALSE. ! w diagnostic prints
LOGICAL :: L_flush6         = .FALSE. ! flush buffers on failure
LOGICAL :: L_diag_print     = .FALSE. ! Print diagnostics
LOGICAL :: L_print_lapse    = .FALSE. ! Print lapse_rate diagnostics
LOGICAL :: L_print_wmax     = .FALSE. ! Print max w diagnostic
LOGICAL :: L_print_theta1   = .FALSE. ! Print level 1 theta diagnostic

LOGICAL :: L_diffusion    = .FALSE. ! Removed from NL, set in code
LOGICAL :: L_upper_ramp   = .FALSE. ! ramp upper-level diffusion
LOGICAL :: L_vdiff_uv     = .FALSE. ! activate targeted diffusion uv
LOGICAL :: L_pftheta      = .FALSE. ! activate polar filter of theta
LOGICAL :: L_pfuv         = .FALSE. ! activate polar filter of u,v
LOGICAL :: L_pfw          = .FALSE. ! activate polar filter of w
LOGICAL :: L_pofil_new    = .FALSE. ! activate new polar filter

LOGICAL :: L_filter      = .FALSE. ! Activate polar filter or diffusion
                                   ! (not in NL)
LOGICAL :: L_filter_incs = .FALSE. ! Activate polar filter of incs (not in NL)

LOGICAL :: L_pfexner  ! Activate polar filter of Exner pressure
                      ! Always equal to L_pofil_hadgem2
                      ! Removed from NL, set in code

LOGICAL :: L_pofil_hadgem2 = .FALSE. ! Use HadGEM2 polar filter setting
LOGICAL :: L_diff_exner              ! activate diffusion of Exner pressure
LOGICAL :: L_pfcomb                  ! combined polar filter/diffusion active
LOGICAL :: L_sponge      = .FALSE.   ! activate lateral bndaries sponge zones
LOGICAL :: L_pfincs                  ! activate polar filter of incs
LOGICAL :: L_diff_thermo = .FALSE.   ! horiz. diffusion of theta
LOGICAL :: L_diff_wind   = .FALSE.   ! horiz. diffusion of u, v 
LOGICAL :: L_diff_wind_ew_only = .FALSE.  ! Omit horiz. diff. in NS dirn. 
LOGICAL :: L_diff_w      = .FALSE.   ! horiz. diffusion of w
LOGICAL :: L_diff_incs   = .FALSE.   ! horiz. diffusion of increments
LOGICAL :: L_diff_auto   = .FALSE.   ! UM calculates diffusion parameters
LOGICAL :: L_tardiff_q   = .FALSE.   ! activate targeted diffusion q

LOGICAL :: L_diff_ctl                ! general diffusion control
                                     ! Removed from NL, set in code

! L_diff_ctl == L_diffusion .or. L_cdiffusion .or. L_vertical_diffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print

LOGICAL :: L_cdiffusion = .FALSE.    ! Removed from NL, set in code

INTEGER :: hdiffopt = imdi ! Horizontal diffusion option
                         ! 0: off, 1: old, 2: conservative, 3: subgrid
INTEGER :: vdiffopt = imdi ! Vertical diffusion option
                         ! 0: off, 1: uniform, 2: ramped

INTEGER :: hdiff_strlev_wind (model_levels_max)=imdi ! Start, end levels for
INTEGER :: hdiff_endlev_wind (model_levels_max)=imdi !  horiz diffn of wind
REAL    :: hdiff_coeff_wind  (model_levels_max)=rmdi ! Horiz diff coeffs wind
INTEGER :: hdiff_order_wind  (model_levels_max)=imdi ! Horiz diff orders wind
INTEGER :: hdiff_strlev_theta(model_levels_max)=imdi ! Start, end levels for
INTEGER :: hdiff_endlev_theta(model_levels_max)=imdi !  horiz diffn of theta
REAL    :: hdiff_coeff_theta (model_levels_max)=rmdi ! Horiz diff coeffs theta
INTEGER :: hdiff_order_theta (model_levels_max)=imdi ! Horiz diff orders theta
INTEGER :: hdiff_strlev_q    (model_levels_max)=imdi ! Start, end levels for
INTEGER :: hdiff_endlev_q    (model_levels_max)=imdi !  horiz diffusion of q
REAL    :: hdiff_coeff_q     (model_levels_max)=rmdi ! Horiz diffn coeffs q
INTEGER :: hdiff_order_q     (model_levels_max)=imdi ! Horiz diffn orders q

INTEGER :: pofil_opt = imdi ! Polar filter option (0:none, 1:combined, 2:old)

INTEGER  :: diffusion_order_thermo(model_levels_max)
INTEGER  :: diffusion_order_wind(model_levels_max)
INTEGER  :: diffusion_order_q(model_levels_max)
INTEGER  :: diffusion_order_w(model_levels_max)
REAL     :: diffusion_coefficient_thermo(model_levels_max)
REAL     :: diffusion_coefficient_wind(model_levels_max)
REAL     :: diffusion_coefficient_q(model_levels_max)
REAL     :: diffusion_coefficient_w(model_levels_max)

INTEGER :: print_step       = imdi ! To control diagnostic printing interval
INTEGER :: diag_interval    = imdi ! diagnostic printing sampling frequency
INTEGER :: norm_lev_start   = imdi ! start level for norm diagnostics
INTEGER :: norm_lev_end     = imdi ! end level for norm diagnostics
INTEGER :: first_norm_print = imdi ! first timestep for norm printing
INTEGER :: dom_w_in = 0     ! define start for block printing
INTEGER :: dom_e_in = 0     ! define end for block printing
INTEGER :: dom_s_in = 0     ! define start for block printing
INTEGER :: dom_n_in = 0     ! define end for block printing
INTEGER :: blockx_in = 0    ! define size for block printing
INTEGER :: blocky_in = 0    ! define size for  block printing

INTEGER :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
INTEGER :: u_end(0:max_121_rows)   ! Sweep control on 121 filter
INTEGER :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
INTEGER :: v_end(0:max_121_rows)   ! Sweep control on 121 filter
INTEGER :: u_sweeps(max_121_rows)  ! Sweep control on 121 filter
INTEGER :: v_sweeps(max_121_rows)  ! Sweep control on 121 filter
INTEGER :: max_sweeps = imdi       ! Max sweeps wanted for 121 filter
INTEGER :: global_u_filter         ! Sweep control on 121 filter
INTEGER :: global_v_filter         ! Sweep control on 121 filter
INTEGER :: diff_order_thermo = imdi      ! diffusion order for theta
INTEGER :: diff_order_wind   = imdi      ! diffusion order for winds
INTEGER :: diff_timescale_thermo = imdi  ! diffusion timescale for theta
INTEGER :: diff_timescale_wind   = imdi  ! diffusion timescale for wind
INTEGER :: vdiffuv_timescale = imdi  ! diffusion e-folding timesteps
INTEGER :: vdiffuv_start     = imdi  ! start level  targeted diffusion
INTEGER :: vdiffuv_end       = imdi  ! end level targeted diffusion
INTEGER :: top_filt_start = imdi ! start level upper-level diffusion
INTEGER :: top_filt_end   = imdi ! end level upper-level diffusion
INTEGER :: sponge_ew = imdi      ! left/right boundaries sponge zone width
INTEGER :: sponge_ns = imdi      ! north/south boundaries sponge zone width
INTEGER :: sponge_power = imdi   ! sponge zone weighting order
INTEGER :: tardiffq_test  = imdi ! test level test w targetted diffusion
INTEGER :: tardiffq_start = imdi ! start level test w targetted diffusion
INTEGER :: tardiffq_end   = imdi ! end level test w targetted diffusion

!  level - assume surfaces are horizontal
INTEGER  :: horizontal_level = imdi
INTEGER  :: tar_horizontal = imdi  ! steep slope test targeted diffusion

REAL :: vdiffuv_test = rmdi ! test to activate shear diffusion of u,v
REAL :: vdiffuv_factor ! vertical diffusn coeff for shear diff (not in NL)
REAL :: scale_ratio = rmdi ! Pass control on 121 filter
REAL :: ref_lat_deg = rmdi ! Reference lat for auto diffusion
REAL :: top_diff = rmdi    !  upper-level diffusion coefficient
REAL :: up_diff_scale = rmdi    !  upper-level diffusion ramping factor
REAL :: adjust_lapse_min = 0.0 !  min dtheta/dz in vertical adjustment

REAL :: up_diff(max_updiff_levels) ! upper-level diffusion coeff
REAL :: sponge_wts_ew(max_sponge_width) ! sponge weights
REAL :: sponge_wts_ns(max_sponge_width) ! sponge weights

REAL :: diff_coeff_ref = rmdi ! EW diffusion coefficient at polar cap
REAL :: diff_coeff_thermo = rmdi  ! NS theta diffusion coeff
REAL :: diff_coeff_wind = rmdi   ! NS u,v diffusion coeff
REAL :: diff_coeff_phi = rmdi   ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
REAL :: polar_cap = rmdi ! Apply 1-2-1 filter polewards
REAL :: tardiffq_factor = rmdi ! targeted diffusion coefficient
REAL :: w_print_limit ! w Threshold for diagnostic printing
REAL :: w_conv_limit = rmdi  ! w Threshold for limiting convection

! Divergence damping control options:
LOGICAL :: L_divdamp = .FALSE.  ! Obsolete option - removed from NL
REAL :: div_damp_coefficient(model_levels_max)

! Polar filter control options:
LOGICAL :: L_polar_filter = .FALSE.
! T: use polar filter to filter increment
LOGICAL :: L_polar_filter_incs = .FALSE.

REAL :: polar_filter_north_lat_limit = rmdi
REAL :: polar_filter_south_lat_limit = rmdi
REAL :: polar_filter_coefficient = rmdi

! amount in radians to increment start latitude by per sweep
REAL :: polar_filter_step_per_sweep = rmdi

! max latitude at which filter can star
REAL :: polar_filter_lat_limit = rmdi

! number of sweeps of filter to do
INTEGER :: polar_filter_n_sweeps = imdi

! Vertical Diffusion control options:
LOGICAL :: L_vertical_diffusion = .FALSE. ! Removed from NL, set in code
LOGICAL :: L_ramp               = .FALSE. ! Removed from NL, set in code

INTEGER :: level_start_wind = imdi
INTEGER :: level_stop_wind = imdi
INTEGER :: level_start_theta = imdi
INTEGER :: level_stop_theta = imdi
INTEGER :: level_start_q = imdi
INTEGER :: level_stop_q = imdi

REAL :: vert_diffusion_coeff_wind = rmdi
REAL :: vert_diffusion_coeff_theta = rmdi
REAL :: vert_diffusion_coeff_q = rmdi
REAL :: ramp_lat_radians = rmdi

! Moisture resetting control options:
LOGICAL :: l_qpos = .FALSE.           ! logical to run qpos code
LOGICAL :: l_qpos_diag_pr = .FALSE.   ! Diagnostic print switch

INTEGER :: q_pos_method = imdi        ! Algorithm choice
INTEGER :: q_pos_tracer_method = imdi ! Algorithm choice for tracers

REAL :: qpos_diag_limit = rmdi        ! Limit for diagnostic prints
REAL :: qlimit = rmdi                 ! lowest allowed value of q

! Inputs for the Leonard terms
LOGICAL :: l_leonard_term = .FALSE.
REAL :: leonard_kl = rmdi

NAMELIST/RUN_Diffusion/                                           &
  hdiffopt, vdiffopt,                                             &
  hdiff_strlev_wind, hdiff_endlev_wind,                           &
  hdiff_coeff_wind, hdiff_order_wind,                             &
  hdiff_strlev_theta, hdiff_endlev_theta,                         &
  hdiff_coeff_theta, hdiff_order_theta,                           &
  hdiff_strlev_q, hdiff_endlev_q,                                 &
  hdiff_coeff_q, hdiff_order_q,                                   &
  pofil_opt,                                                      &
  diff_order_thermo, diff_timescale_thermo,                       &
  diff_order_wind, diff_timescale_wind,                           &
  horizontal_level, tar_horizontal,                               &
  ramp_lat_radians,                                               &
  L_polar_filter,L_polar_filter_incs,                             &
  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
  polar_filter_coefficient,polar_filter_step_per_sweep,           &
  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
  l_qpos,  q_pos_method, q_pos_tracer_method, qlimit,             &
  l_qpos_diag_pr, qpos_diag_limit,                                &
  L_tardiff_q, w_conv_limit,                                      &
  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
  L_subfilter_horiz, L_subfilter_vert, hyd_mix_opt,               &
  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
  turb_startlev_vert, turb_endlev_vert,                           &
  L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pofil_hadgem2,            &
  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_wind_ew_only,   &
  L_diff_w,                                                       &
  level_start_wind, level_stop_wind,                              &
  level_start_theta, level_stop_theta,                            &
  level_start_q, level_stop_q,                                    &
  L_pofil_new, L_diff_auto,                                       &
  diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,            &
  max_sweeps, L_upper_ramp, up_diff_scale, top_diff,              &
  top_filt_start, top_filt_end, L_vdiff_uv, vdiffuv_timescale,    &
  vdiffuv_test, vdiffuv_start, vdiffuv_end,                       &
  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
  vert_diffusion_coeff_q,                                         &
  L_diag_print, L_diag_print_ops, L_print_pe,                     &
  L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,         &
  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
  L_diag_noise, L_diag_L2norms, L_diag_L2helm,                    &
  norm_lev_start, norm_lev_end, first_norm_print,                 &
  print_step, diag_interval, w_print_limit, L_Flush6,             &
  l_leonard_term, leonard_kl

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TURB_DIFF_MOD'

CONTAINS

SUBROUTINE check_run_diffusion()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the run_diffusion namelist.

USE ereport_mod, ONLY: ereport
USE bl_option_mod, ONLY: non_local_bl, blending_option, off,         &
                         i_bl_vn, i_bl_vn_0
USE cv_run_mod, ONLY: cldbase_opt_sh
USE cv_param_mod, ONLY: sh_wstar_closure
USE conversions_mod, ONLY: pi_over_180
USE chk_opts_mod, ONLY: chk_var, def_src

IMPLICIT NONE

INTEGER                       :: ErrorStatus   ! used for ereport
INTEGER                       :: i,j 
CHARACTER (LEN=errormessagelength)    :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CHECK_RUN_DIFFUSION'
    

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! set BL parameters dependent on diffusion options
IF (l_subfilter_vert .AND. blending_option == off) THEN
  non_local_bl = off
  cmessage =                                                 newline//&
  'Non_local_bl switched off since vertical Smagorinsky'//   newline//&
  'chosen without blending.'
  ErrorStatus = -100
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

! set convection parameters dependent on diffusion options
IF (.NOT. l_subfilter_horiz) THEN
  cldbase_opt_sh = sh_wstar_closure
  cmessage =                                                 newline//&
  'cldbase_opt_sh set to sh_wstar_closure since'//           newline//&
  'Smagorinsky diffusion not chosen.'
  ErrorStatus = -100
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

! check horizontal Smagorinsky options
IF (l_subfilter_horiz) THEN
  CALL chk_var(hyd_mix_opt,'hyd_mix_opt',                               &
       [0,liq_only,liq_and_ice,all_hydro])
END IF

! UMUI-to-ROSE transition UM vn9.0.
! Under the UMUI regime, the variable hdiffopt was not in the name list but
! was an internal umui variable used to set a number of logicals.
! To facilitate transition to ROSE, hdiffopt was moved into the name list,
! while the dependent logicals were removed from the it.
! The value of hdiffopt is used here to set dependent logicals.
! Logical l_diffusion switches on the "old diffusion option", while
! l_cdiffusion switches on the "conservative diffusion option".
IF (hdiffopt==1) l_diffusion = .TRUE.
IF (hdiffopt==2) l_cdiffusion = .TRUE.

IF (hdiffopt==1 .OR. hdiffopt==2) THEN
  DO i = 1, model_levels_max
    IF (hdiff_strlev_wind(i) > 0) THEN
      DO j = hdiff_strlev_wind(i), hdiff_endlev_wind(i)
        diffusion_coefficient_wind(j) = hdiff_coeff_wind(i)
        diffusion_order_wind      (j) = hdiff_order_wind(i)
      END DO
    END IF
    IF (hdiff_strlev_theta(i) > 0) THEN
      DO j = hdiff_strlev_theta(i), hdiff_endlev_theta(i)
        diffusion_coefficient_thermo(j) = hdiff_coeff_theta(i)
        diffusion_order_thermo      (j) = hdiff_order_theta(i)
      END DO
    END IF
    IF (hdiff_strlev_q(i) > 0) THEN
      DO j = hdiff_strlev_q(i), hdiff_endlev_q(i)
        diffusion_coefficient_q(j) = hdiff_coeff_q(i)
        diffusion_order_q      (j) = hdiff_order_q(i)
      END DO
    END IF
  END DO
END IF

IF (vdiffopt==1 .OR. vdiffopt==2) l_vertical_diffusion = .TRUE.
IF (vdiffopt==2) l_ramp = .TRUE.

IF (hdiffopt /=0 .OR. vdiffopt /=0 .OR. l_tardiff_q  &
     .OR. l_diag_print) THEN
  L_diff_ctl = .TRUE.
ELSE
  L_diff_ctl = .FALSE.
END IF

! These two logicals were both in the NL pre-UM 9.0, but the umui
!  code gave them the same value. One of them has been removed from
!  the NL under ROSE rationalisation.
l_pfexner = l_pofil_hadgem2


! Convert from degrees to radians
polar_filter_north_lat_limit= polar_filter_north_lat_limit * pi_over_180
polar_filter_south_lat_limit= polar_filter_south_lat_limit * pi_over_180
polar_filter_lat_limit      = polar_filter_lat_limit       * pi_over_180
polar_filter_step_per_sweep = polar_filter_step_per_sweep  * pi_over_180
ramp_lat_radians            = ramp_lat_radians             * pi_over_180

! Check that Leonard term parameter within allowed range, if used.
IF (l_leonard_term)  CALL chk_var( leonard_kl, 'leonard_kl', '[0.0:6.0]' )

! Check not trying to use Leonard terms without BL scheme; won't work as
! Leonard term calculation uses some arrays calculated in the BL scheme
IF ( l_leonard_term .AND. i_bl_vn == i_bl_vn_0 ) THEN
  cmessage =                                                 newline//&
  'Cannot run with Leonard terms switched on'              //newline//&
  '(l_leonard_term = .TRUE.)'                              //newline//&
  'but the boundary-layer scheme switched off'             //newline//&
  '(i_bl_vn == i_bl_vn_0).'                                //newline//&
  'You must switch on the BL scheme if you want to run'    //newline//&
  'with Leonard terms, as the Leonard term calculations'   //newline//&
  'depend on outputs from the BL code.'
  ErrorStatus = 100
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF


def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE check_run_diffusion

SUBROUTINE print_nlist_run_diffusion()
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_DIFFUSION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_diffusion',                         &
    src='nl_cruntimc_h_run_diffusion')

WRITE(lineBuffer,'(A,I0)')' hdiffopt = ',hdiffopt
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,'(A,I0)')' vdiffopt = ',vdiffopt
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,'(A,I0)')' pofil_opt = ',pofil_opt
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diffusion_order_thermo = ',diffusion_order_thermo
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diffusion_order_wind = ',diffusion_order_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diffusion_order_q = ',diffusion_order_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' diffusion_coefficient_thermo = ',diffusion_coefficient_thermo
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' diffusion_coefficient_wind = ',diffusion_coefficient_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diffusion_coefficient_q = ',diffusion_coefficient_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_order_thermo = ',diff_order_thermo
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_timescale_thermo = ',diff_timescale_thermo
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_order_wind = ',diff_order_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_timescale_wind = ',diff_timescale_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' horizontal_level = ',horizontal_level
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' tar_horizontal = ',tar_horizontal
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' ramp_lat_radians = ',ramp_lat_radians
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_polar_filter = ',L_polar_filter
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_polar_filter_incs = ',L_polar_filter_incs
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' polar_filter_north_lat_limit = ',polar_filter_north_lat_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' polar_filter_south_lat_limit = ',polar_filter_south_lat_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' polar_filter_coefficient = ',polar_filter_coefficient
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' polar_filter_step_per_sweep = ',polar_filter_step_per_sweep
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' polar_filter_lat_limit = ',polar_filter_lat_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' polar_filter_n_sweeps = ',polar_filter_n_sweeps
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' l_qpos = ',l_qpos
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' q_pos_method = ',q_pos_method
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' q_pos_tracer_method = ',q_pos_tracer_method
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' qlimit = ',qlimit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' l_qpos_diag_pr = ',l_qpos_diag_pr
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' qpos_diag_limit = ',qpos_diag_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_tardiff_q = ',L_tardiff_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' w_conv_limit = ',w_conv_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' tardiffq_factor = ',tardiffq_factor
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' tardiffq_test = ',tardiffq_test
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' tardiffq_start = ',tardiffq_start
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' tardiffq_end = ',tardiffq_end
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_subfilter_horiz = ',L_subfilter_horiz
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_subfilter_vert = ',L_subfilter_vert
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,'(A,I0)')' hyd_mix_opt = ',hyd_mix_opt
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_factor = ',diff_factor
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' mix_factor = ',mix_factor
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' turb_startlev_horiz = ',turb_startlev_horiz
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' turb_endlev_horiz = ',turb_endlev_horiz
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' turb_startlev_vert = ',turb_startlev_vert
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' turb_endlev_vert = ',turb_endlev_vert
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pftheta = ',L_pftheta
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pfuv = ',L_pfuv
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pfw = ',L_pfw
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pfincs = ',L_pfincs
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pofil_hadgem2 = ',L_pofil_hadgem2
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_incs = ',L_diff_incs
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_thermo = ',L_diff_thermo
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_wind = ',L_diff_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_wind_ew_only = ',L_diff_wind_ew_only
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_w = ',L_diff_w
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_start_wind = ',level_start_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_stop_wind = ',level_stop_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_start_theta = ',level_start_theta
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_stop_theta = ',level_stop_theta
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_start_q = ',level_start_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' level_stop_q = ',level_stop_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_pofil_new = ',L_pofil_new
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diff_auto = ',L_diff_auto
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diff_coeff_ref = ',diff_coeff_ref
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' polar_cap = ',polar_cap
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' scale_ratio = ',scale_ratio
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' ref_lat_deg = ',ref_lat_deg
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' max_sweeps = ',max_sweeps
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_upper_ramp = ',L_upper_ramp
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' up_diff_scale = ',up_diff_scale
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' top_diff = ',top_diff
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' top_filt_start = ',top_filt_start
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' top_filt_end = ',top_filt_end
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_vdiff_uv = ',L_vdiff_uv
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vdiffuv_timescale = ',vdiffuv_timescale
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vdiffuv_test = ',vdiffuv_test
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vdiffuv_factor = ',vdiffuv_factor
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vdiffuv_start = ',vdiffuv_start
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vdiffuv_end = ',vdiffuv_end
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_sponge = ',L_sponge
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' sponge_power = ',sponge_power
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' sponge_ew = ',sponge_ew
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' sponge_ns = ',sponge_ns
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vert_diffusion_coeff_wind = ',vert_diffusion_coeff_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)                                                        &
    ' vert_diffusion_coeff_theta = ',vert_diffusion_coeff_theta
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' vert_diffusion_coeff_q = ',vert_diffusion_coeff_q
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_print = ',L_diag_print
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_print_ops = ',L_diag_print_ops
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_pe = ',L_print_pe
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_w = ',L_print_w
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_wmax = ',L_print_wmax
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_lapse = ',L_print_lapse
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_theta1 = ',L_print_theta1
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_div = ',L_print_div
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_wind = ',L_diag_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_shear = ',L_print_shear
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_print_max_wind = ',L_print_max_wind
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_noise = ',L_diag_noise
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_L2norms = ',L_diag_L2norms
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_diag_L2helm = ',L_diag_L2helm
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' norm_lev_start = ',norm_lev_start
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' norm_lev_end = ',norm_lev_end
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' first_norm_print = ',first_norm_print
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' print_step = ',print_step
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' diag_interval = ',diag_interval
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' w_print_limit = ',w_print_limit
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' L_Flush6 = ',L_Flush6
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' l_leonard_term = ',l_leonard_term
CALL umPrint(lineBuffer,src='turb_diff_mod')
WRITE(lineBuffer,*)' leonard_kl = ',leonard_kl
CALL umPrint(lineBuffer,src='turb_diff_mod')


CALL umPrint('- - - - - - end of namelist - - - - - -',                    &
    src='nl_cruntimc_h_run_diffusion')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_diffusion

#if !defined(LFRIC)
SUBROUTINE read_nml_run_diffusion(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_DIFFUSION'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 40 + 9 * model_levels_max
INTEGER, PARAMETER :: n_real = 24 + 3 * model_levels_max
INTEGER, PARAMETER :: n_log = 38

TYPE my_namelist
  SEQUENCE
  INTEGER :: hdiffopt
  INTEGER :: vdiffopt
  INTEGER :: hdiff_strlev_wind(model_levels_max)
  INTEGER :: hdiff_endlev_wind(model_levels_max)
  INTEGER :: hdiff_order_wind(model_levels_max)
  INTEGER :: hdiff_strlev_theta(model_levels_max)
  INTEGER :: hdiff_endlev_theta(model_levels_max)
  INTEGER :: hdiff_order_theta(model_levels_max)
  INTEGER :: hdiff_strlev_q(model_levels_max)
  INTEGER :: hdiff_endlev_q(model_levels_max)
  INTEGER :: hdiff_order_q(model_levels_max)
  INTEGER :: pofil_opt
  INTEGER :: diff_order_thermo
  INTEGER :: diff_timescale_thermo
  INTEGER :: diff_order_wind
  INTEGER :: diff_timescale_wind
  INTEGER :: horizontal_level
  INTEGER :: tar_horizontal
  INTEGER :: polar_filter_n_sweeps
  INTEGER :: q_pos_method
  INTEGER :: q_pos_tracer_method
  INTEGER :: tardiffq_test
  INTEGER :: tardiffq_start
  INTEGER :: tardiffq_end
  INTEGER :: turb_startlev_horiz
  INTEGER :: turb_endlev_horiz
  INTEGER :: turb_startlev_vert
  INTEGER :: turb_endlev_vert
  INTEGER :: level_start_wind
  INTEGER :: level_stop_wind
  INTEGER :: level_start_theta
  INTEGER :: level_stop_theta
  INTEGER :: level_start_q
  INTEGER :: level_stop_q
  INTEGER :: max_sweeps
  INTEGER :: top_filt_start
  INTEGER :: top_filt_end
  INTEGER :: vdiffuv_timescale
  INTEGER :: vdiffuv_start
  INTEGER :: vdiffuv_end
  INTEGER :: sponge_power
  INTEGER :: sponge_ew
  INTEGER :: sponge_ns
  INTEGER :: norm_lev_start
  INTEGER :: norm_lev_end
  INTEGER :: first_norm_print
  INTEGER :: print_step
  INTEGER :: diag_interval
  INTEGER :: hyd_mix_opt
  REAL :: hdiff_coeff_wind(model_levels_max)
  REAL :: hdiff_coeff_theta(model_levels_max)
  REAL :: hdiff_coeff_q(model_levels_max)
  REAL :: ramp_lat_radians
  REAL :: polar_filter_north_lat_limit
  REAL :: polar_filter_south_lat_limit
  REAL :: polar_filter_coefficient
  REAL :: polar_filter_step_per_sweep
  REAL :: polar_filter_lat_limit
  REAL :: qlimit
  REAL :: qpos_diag_limit
  REAL :: w_conv_limit
  REAL :: tardiffq_factor
  REAL :: diff_factor
  REAL :: mix_factor
  REAL :: diff_coeff_ref
  REAL :: polar_cap
  REAL :: scale_ratio
  REAL :: ref_lat_deg
  REAL :: up_diff_scale
  REAL :: top_diff
  REAL :: vdiffuv_test
  REAL :: vert_diffusion_coeff_wind
  REAL :: vert_diffusion_coeff_theta
  REAL :: vert_diffusion_coeff_q
  REAL :: w_print_limit
  REAL :: leonard_kl
  LOGICAL :: L_polar_filter
  LOGICAL :: L_polar_filter_incs
  LOGICAL :: l_qpos
  LOGICAL :: l_qpos_diag_pr
  LOGICAL :: L_tardiff_q
  LOGICAL :: L_subfilter_horiz
  LOGICAL :: L_subfilter_vert
  LOGICAL :: L_pftheta
  LOGICAL :: L_pfuv
  LOGICAL :: L_pfw
  LOGICAL :: L_pfincs
  LOGICAL :: L_pofil_hadgem2
  LOGICAL :: L_diff_incs
  LOGICAL :: L_diff_thermo
  LOGICAL :: L_diff_wind
  LOGICAL :: L_diff_wind_ew_only
  LOGICAL :: L_diff_w
  LOGICAL :: L_pofil_new
  LOGICAL :: L_diff_auto
  LOGICAL :: L_upper_ramp
  LOGICAL :: L_vdiff_uv
  LOGICAL :: L_sponge
  LOGICAL :: L_diag_print
  LOGICAL :: L_diag_print_ops
  LOGICAL :: L_print_pe
  LOGICAL :: L_print_w
  LOGICAL :: L_print_wmax
  LOGICAL :: L_print_lapse
  LOGICAL :: L_print_theta1
  LOGICAL :: L_print_div
  LOGICAL :: L_diag_wind
  LOGICAL :: L_print_shear
  LOGICAL :: L_print_max_wind
  LOGICAL :: L_diag_noise
  LOGICAL :: L_diag_L2norms
  LOGICAL :: L_diag_L2helm
  LOGICAL :: L_Flush6
  LOGICAL :: l_leonard_term
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Diffusion, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Diffusion", iomessage)

  my_nml % hdiffopt           = hdiffopt
  my_nml % vdiffopt           = vdiffopt
  my_nml % hdiff_strlev_wind  = hdiff_strlev_wind
  my_nml % hdiff_endlev_wind  = hdiff_endlev_wind
  my_nml % hdiff_order_wind   = hdiff_order_wind
  my_nml % hdiff_strlev_theta = hdiff_strlev_theta
  my_nml % hdiff_endlev_theta = hdiff_endlev_theta
  my_nml % hdiff_order_theta  = hdiff_order_theta
  my_nml % hdiff_strlev_q     = hdiff_strlev_q
  my_nml % hdiff_endlev_q     = hdiff_endlev_q
  my_nml % hdiff_order_q      = hdiff_order_q
  my_nml % pofil_opt          = pofil_opt
  my_nml % diff_order_thermo  = diff_order_thermo
  my_nml % diff_timescale_thermo = diff_timescale_thermo
  my_nml % diff_order_wind    = diff_order_wind
  my_nml % diff_timescale_wind = diff_timescale_wind
  my_nml % horizontal_level   = horizontal_level
  my_nml % tar_horizontal     = tar_horizontal
  my_nml % polar_filter_n_sweeps = polar_filter_n_sweeps
  my_nml % q_pos_method       = q_pos_method
  my_nml % q_pos_tracer_method = q_pos_tracer_method
  my_nml % tardiffq_test      = tardiffq_test
  my_nml % tardiffq_start     = tardiffq_start
  my_nml % tardiffq_end       = tardiffq_end
  my_nml % turb_startlev_horiz = turb_startlev_horiz
  my_nml % turb_endlev_horiz  = turb_endlev_horiz
  my_nml % turb_startlev_vert = turb_startlev_vert
  my_nml % turb_endlev_vert   = turb_endlev_vert
  my_nml % level_start_wind   = level_start_wind
  my_nml % level_stop_wind    = level_stop_wind
  my_nml % level_start_theta  = level_start_theta
  my_nml % level_stop_theta   = level_stop_theta
  my_nml % level_start_q      = level_start_q
  my_nml % level_stop_q       = level_stop_q
  my_nml % max_sweeps         = max_sweeps
  my_nml % top_filt_start     = top_filt_start
  my_nml % top_filt_end       = top_filt_end
  my_nml % vdiffuv_timescale  = vdiffuv_timescale
  my_nml % vdiffuv_start      = vdiffuv_start
  my_nml % vdiffuv_end        = vdiffuv_end
  my_nml % sponge_power       = sponge_power
  my_nml % sponge_ew          = sponge_ew
  my_nml % sponge_ns          = sponge_ns
  my_nml % norm_lev_start     = norm_lev_start
  my_nml % norm_lev_end       = norm_lev_end
  my_nml % first_norm_print   = first_norm_print
  my_nml % print_step         = print_step
  my_nml % diag_interval      = diag_interval
  my_nml % hyd_mix_opt        = hyd_mix_opt
  ! end of integers
  my_nml % hdiff_coeff_wind   = hdiff_coeff_wind
  my_nml % hdiff_coeff_theta  = hdiff_coeff_theta
  my_nml % hdiff_coeff_q      = hdiff_coeff_q
  my_nml % ramp_lat_radians   = ramp_lat_radians
  my_nml % polar_filter_north_lat_limit = polar_filter_north_lat_limit
  my_nml % polar_filter_south_lat_limit = polar_filter_south_lat_limit
  my_nml % polar_filter_coefficient     = polar_filter_coefficient
  my_nml % polar_filter_step_per_sweep  = polar_filter_step_per_sweep
  my_nml % polar_filter_lat_limit       = polar_filter_lat_limit
  my_nml % qlimit             = qlimit
  my_nml % qpos_diag_limit    = qpos_diag_limit
  my_nml % w_conv_limit       = w_conv_limit
  my_nml % tardiffq_factor    = tardiffq_factor
  my_nml % diff_factor        = diff_factor
  my_nml % mix_factor         = mix_factor
  my_nml % diff_coeff_ref     = diff_coeff_ref
  my_nml % polar_cap          = polar_cap
  my_nml % scale_ratio        = scale_ratio
  my_nml % ref_lat_deg        = ref_lat_deg
  my_nml % up_diff_scale      = up_diff_scale
  my_nml % top_diff           = top_diff
  my_nml % vdiffuv_test       = vdiffuv_test
  my_nml % vert_diffusion_coeff_wind  = vert_diffusion_coeff_wind
  my_nml % vert_diffusion_coeff_theta = vert_diffusion_coeff_theta
  my_nml % vert_diffusion_coeff_q     = vert_diffusion_coeff_q
  my_nml % w_print_limit      = w_print_limit
  my_nml % leonard_kl         = leonard_kl
  ! end of reals
  my_nml % L_polar_filter      = L_polar_filter
  my_nml % L_polar_filter_incs = L_polar_filter_incs
  my_nml % l_qpos              = l_qpos
  my_nml % l_qpos_diag_pr      = l_qpos_diag_pr
  my_nml % L_tardiff_q         = L_tardiff_q
  my_nml % L_subfilter_horiz   = L_subfilter_horiz
  my_nml % L_subfilter_vert    = L_subfilter_vert
  my_nml % L_pftheta           = L_pftheta
  my_nml % L_pfuv              = L_pfuv
  my_nml % L_pfw               = L_pfw
  my_nml % L_pfincs            = L_pfincs
  my_nml % L_pofil_hadgem2     = L_pofil_hadgem2
  my_nml % L_diff_incs         = L_diff_incs
  my_nml % L_diff_thermo       = L_diff_thermo
  my_nml % L_diff_wind         = L_diff_wind
  my_nml % L_diff_wind_ew_only = L_diff_wind_ew_only
  my_nml % L_diff_w            = L_diff_w
  my_nml % L_pofil_new         = L_pofil_new
  my_nml % L_diff_auto         = L_diff_auto
  my_nml % L_upper_ramp        = L_upper_ramp
  my_nml % L_vdiff_uv          = L_vdiff_uv
  my_nml % L_sponge            = L_sponge
  my_nml % L_diag_print        = L_diag_print
  my_nml % L_diag_print_ops    = L_diag_print_ops
  my_nml % L_print_pe          = L_print_pe
  my_nml % L_print_w           = L_print_w
  my_nml % L_print_wmax        = L_print_wmax
  my_nml % L_print_lapse       = L_print_lapse
  my_nml % L_print_theta1      = L_print_theta1
  my_nml % L_print_div         = L_print_div
  my_nml % L_diag_wind         = L_diag_wind
  my_nml % L_print_shear       = L_print_shear
  my_nml % L_print_max_wind    = L_print_max_wind
  my_nml % L_diag_noise        = L_diag_noise
  my_nml % L_diag_L2norms      = L_diag_L2norms
  my_nml % L_diag_L2helm       = L_diag_L2helm
  my_nml % L_Flush6            = L_Flush6
  my_nml % l_leonard_term      = l_leonard_term

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  hdiffopt          = my_nml % hdiffopt
  vdiffopt          = my_nml % vdiffopt
  hdiff_strlev_wind = my_nml % hdiff_strlev_wind
  hdiff_endlev_wind = my_nml % hdiff_endlev_wind
  hdiff_order_wind  = my_nml % hdiff_order_wind
  hdiff_strlev_theta = my_nml % hdiff_strlev_theta
  hdiff_endlev_theta = my_nml % hdiff_endlev_theta
  hdiff_order_theta  = my_nml % hdiff_order_theta
  hdiff_strlev_q    = my_nml % hdiff_strlev_q
  hdiff_endlev_q    = my_nml % hdiff_endlev_q
  hdiff_order_q     = my_nml % hdiff_order_q
  pofil_opt         = my_nml % pofil_opt
  diff_order_thermo     = my_nml %diff_order_thermo
  diff_timescale_thermo = my_nml % diff_timescale_thermo
  diff_order_wind   = my_nml % diff_order_wind
  diff_timescale_wind = my_nml % diff_timescale_wind
  horizontal_level  = my_nml % horizontal_level
  tar_horizontal    = my_nml % tar_horizontal
  polar_filter_n_sweeps = my_nml % polar_filter_n_sweeps
  q_pos_method      = my_nml % q_pos_method
  q_pos_tracer_method = my_nml % q_pos_tracer_method
  tardiffq_test     = my_nml % tardiffq_test
  tardiffq_start    = my_nml % tardiffq_start
  tardiffq_end      = my_nml % tardiffq_end
  turb_startlev_horiz = my_nml % turb_startlev_horiz
  turb_endlev_horiz   = my_nml % turb_endlev_horiz
  turb_startlev_vert  = my_nml % turb_startlev_vert
  turb_endlev_vert  = my_nml % turb_endlev_vert
  level_start_wind  = my_nml % level_start_wind
  level_stop_wind   = my_nml % level_stop_wind
  level_start_theta = my_nml % level_start_theta
  level_stop_theta  = my_nml % level_stop_theta
  level_start_q     = my_nml % level_start_q
  level_stop_q      = my_nml % level_stop_q
  max_sweeps        = my_nml % max_sweeps
  top_filt_start    = my_nml % top_filt_start
  top_filt_end      = my_nml % top_filt_end
  vdiffuv_timescale = my_nml % vdiffuv_timescale
  vdiffuv_start     = my_nml % vdiffuv_start
  vdiffuv_end       = my_nml % vdiffuv_end
  sponge_power      = my_nml % sponge_power
  sponge_ew         = my_nml % sponge_ew
  sponge_ns         = my_nml % sponge_ns
  norm_lev_start    = my_nml % norm_lev_start
  norm_lev_end      = my_nml % norm_lev_end
  first_norm_print  = my_nml % first_norm_print
  print_step        = my_nml % print_step
  diag_interval     = my_nml % diag_interval
  hyd_mix_opt       = my_nml % hyd_mix_opt
  ! end of integers
  hdiff_coeff_wind  = my_nml % hdiff_coeff_wind
  hdiff_coeff_theta = my_nml % hdiff_coeff_theta
  hdiff_coeff_q     = my_nml % hdiff_coeff_q
  ramp_lat_radians  = my_nml % ramp_lat_radians
  polar_filter_north_lat_limit = my_nml % polar_filter_north_lat_limit
  polar_filter_south_lat_limit = my_nml % polar_filter_south_lat_limit
  polar_filter_coefficient     = my_nml % polar_filter_coefficient
  polar_filter_step_per_sweep  = my_nml % polar_filter_step_per_sweep
  polar_filter_lat_limit       = my_nml % polar_filter_lat_limit
  qlimit            = my_nml % qlimit
  qpos_diag_limit   = my_nml % qpos_diag_limit
  w_conv_limit      = my_nml % w_conv_limit
  tardiffq_factor   = my_nml % tardiffq_factor
  diff_factor       = my_nml % diff_factor
  mix_factor        = my_nml % mix_factor
  diff_coeff_ref    = my_nml % diff_coeff_ref
  polar_cap         = my_nml % polar_cap
  scale_ratio       = my_nml % scale_ratio
  ref_lat_deg       = my_nml % ref_lat_deg
  up_diff_scale     = my_nml % up_diff_scale
  top_diff          = my_nml % top_diff
  vdiffuv_test      = my_nml % vdiffuv_test
  vert_diffusion_coeff_wind  = my_nml % vert_diffusion_coeff_wind
  vert_diffusion_coeff_theta = my_nml % vert_diffusion_coeff_theta
  vert_diffusion_coeff_q     = my_nml % vert_diffusion_coeff_q
  w_print_limit     = my_nml % w_print_limit
  leonard_kl        = my_nml % leonard_kl
  ! end of reals
  L_polar_filter      = my_nml % L_polar_filter
  L_polar_filter_incs = my_nml % L_polar_filter_incs
  l_qpos              = my_nml % l_qpos
  l_qpos_diag_pr      = my_nml % l_qpos_diag_pr
  L_tardiff_q         = my_nml % L_tardiff_q
  L_subfilter_horiz   = my_nml % L_subfilter_horiz
  L_subfilter_vert    = my_nml % L_subfilter_vert
  L_pftheta           = my_nml % L_pftheta
  L_pfuv              = my_nml % L_pfuv
  L_pfw               = my_nml % L_pfw
  L_pfincs            = my_nml % L_pfincs
  L_pofil_hadgem2     = my_nml % L_pofil_hadgem2
  L_diff_incs         = my_nml % L_diff_incs
  L_diff_thermo       = my_nml % L_diff_thermo
  L_diff_wind         = my_nml % L_diff_wind
  L_diff_wind_ew_only = my_nml % L_diff_wind_ew_only
  L_diff_w            = my_nml % L_diff_w
  L_pofil_new         = my_nml % L_pofil_new
  L_diff_auto         = my_nml % L_diff_auto
  L_upper_ramp        = my_nml % L_upper_ramp
  L_vdiff_uv          = my_nml % L_vdiff_uv
  L_sponge            = my_nml % L_sponge
  L_diag_print        = my_nml % L_diag_print
  L_diag_print_ops    = my_nml % L_diag_print_ops
  L_print_pe          = my_nml % L_print_pe
  L_print_w           = my_nml % L_print_w
  L_print_wmax        = my_nml % L_print_wmax
  L_print_lapse       = my_nml % L_print_lapse
  L_print_theta1      = my_nml % L_print_theta1
  L_print_div         = my_nml % L_print_div
  L_diag_wind         = my_nml % L_diag_wind
  L_print_shear       = my_nml % L_print_shear
  L_print_max_wind    = my_nml % L_print_max_wind
  L_diag_noise        = my_nml % L_diag_noise
  L_diag_L2norms      = my_nml % L_diag_L2norms
  L_diag_L2helm       = my_nml % L_diag_L2helm
  L_Flush6            = my_nml % L_Flush6
  l_leonard_term      = my_nml % l_leonard_term

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_diffusion
#endif

END MODULE turb_diff_mod
