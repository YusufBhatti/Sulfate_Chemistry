! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE g_wave_5a_mod

USE Field_Types

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='G_WAVE_5A_MOD'

CONTAINS
!
! Calls components of version 5A of gravity wave drag scheme.
!
SUBROUTINE g_wave_5a(                                           &
  theta, u, v, row_length, rows, nrows,u_rows,v_rows,           &
  off_x, off_y,                                                 &
  global_row_length,n_proc, n_procy, proc_row_group,            &
  at_extremity, levels, rho,                                    &
  delta_lambda, delta_phi, true_latitude,                       &
  sd_orog_land, orog_grad_xx_land, orog_grad_xy_land,           &
  orog_grad_yy_land,                                            &
  land_index, land_points, timestep,                            &
  kay, frc, nsigma, r_u, r_v, T_inc, l_taus_scale,              & 
  l_fix_gwsatn, l_gwd_40km, sat_scheme, fsat, gw_seg_size,      &
  Gsharp,fbcd, l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,   &
! diagnostics
        stress_ud,     stress_ud_on,     stress_ud_p_on,        &
                                         points_stress_ud,      &
        stress_vd,     stress_vd_on,     points_stress_vd,      &
        stress_ud_satn,stress_ud_satn_on,points_stress_ud_satn, &
        stress_vd_satn,stress_vd_satn_on,points_stress_vd_satn, &
        stress_ud_wake,stress_ud_wake_on,points_stress_ud_wake, &
        stress_vd_wake,stress_vd_wake_on,points_stress_vd_wake, &
        du_dt_satn,    du_dt_satn_on,    du_dt_satn_p_on,       &
                                         points_du_dt_satn,     &
        dv_dt_satn,    dv_dt_satn_on,    points_dv_dt_satn,     &
        du_dt_wake,    du_dt_wake_on,    points_du_dt_wake,     &
        dv_dt_wake,    dv_dt_wake_on,    points_dv_dt_wake,     &
        u_s_d,         u_s_d_on,         points_u_s_d,          &
        v_s_d,         v_s_d_on,         points_v_s_d,          &
        nsq_s_d,       nsq_s_d_on,       points_nsq_s_d,        &
        fr_d,          fr_d_on,          points_fr_d,           &
        bld_d,         bld_d_on,         points_bld_d,          &
        bldt_d,        bldt_d_on,        points_bldt_d,         &
        num_lim_d,     num_lim_d_on,     points_num_lim_d,      &
        num_fac_d,     num_fac_d_on,     points_num_fac_d,      &
        tausx_d,       tausx_d_on,       points_tausx_d,        &
        tausy_d,       tausy_d_on,       points_tausy_d,        &
        taus_scale_d,  taus_scale_d_on,  points_taus_scale_d,   &
        orog_slope_d,  orog_slope_d_on,  points_orog_slope_d,   &
        orog_anis_d,   orog_anis_d_on,   points_orog_anis_d,    &
        orog_dir_d,    orog_dir_d_on,    points_orog_dir_d,     &
        iret)

USE um_parparams, ONLY: nodomain, pnorth, peast, psouth, pwest
USE atm_fields_bounds_mod, ONLY:               &
  udims, vdims, wdims, tdims, pdims,           &
  udims_s, vdims_s, wdims_s, tdims_s, pdims_s, array_dims

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels,            &  ! Radii on theta levels (m)
  r_rho_levels                  ! Radii on rho levels (m)
USE planet_constants_mod, ONLY: recip_a2

USE eg_v_at_poles_mod, ONLY: eg_v_at_poles

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
USE p_to_u_mod, ONLY: p_to_u
USE p_to_v_mod, ONLY: p_to_v
USE polar_row_mean_mod, ONLY: polar_row_mean
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p

USE non_blocking_halo_exchange, ONLY: &
  begin_swap_bounds, end_swap_bounds

USE model_domain_mod, ONLY: model_type, mt_global

USE gw_block_mod, ONLY: gw_block
USE gw_setup_mod, ONLY: gw_setup
USE gw_wave_mod, ONLY: gw_wave

IMPLICIT NONE
!
! Description:
! 1) Interpolate winds to theta points
!    gather data for land points only
! 2) Call routine to set-up variables for the gravity wave and
!    flow blocking routines
! 3) Calculate stress profiles due to different components of the
!    scheme. Calculate associated wind increments.
! 4) Interpolate acceleration to wind points and update winds
! 5) Gather diagnostics from land points and interpolate from p to
!    u,v staggering on 'c' grid
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Gravity Wave Drag
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.


!
! SUBROUTINE ARGUMENTS
!
!----------------------
! Intent(in) variables
!----------------------

INTEGER, INTENT(IN)  ::        &!block for intent(in)
  row_length,                  &! number of points per row
  rows,                        &! number of rows on theta grid
  nrows,                       &! number of rows on v grid
  u_rows,                      &! rows on u grid
  v_rows,                      &! rows on v grid
  off_x,                       &! small x-halo
  off_y,                       &! small y-halo
  global_row_length,           &! number of points on a row
  proc_row_group,              &! Group id for processors on the same row
  n_proc,                      &! Total number of processors
  n_procy,                     &! Number of processors in latitude
  levels,                      &! number of model levels
  land_points,                 &! number of land points
  land_index(rows*row_length)   ! index for land points

INTEGER, INTENT(IN) :: gw_seg_size ! segment size for memory/OpenMP optimisation

INTEGER, INTENT(IN) ::         &!intent(in) 4a scheme only
  sat_scheme                    !Switch to determine whether to use
                                !amplitude or stress based saturation test

INTEGER, INTENT(IN) ::         &! intent(in) for diags
  points_stress_ud,            &
  points_stress_vd,            &
  points_stress_ud_satn,       &
  points_stress_vd_satn,       &
  points_stress_ud_wake,       &
  points_stress_vd_wake,       &
  points_du_dt_satn,           &
  points_dv_dt_satn,           &
  points_du_dt_wake,           &
  points_dv_dt_wake,           &
  points_u_s_d,                &
  points_v_s_d,                &
  points_nsq_s_d,              &
  points_fr_d,                 &
  points_bld_d,                &
  points_bldt_d,               &
  points_num_lim_d,            &
  points_num_fac_d,            &
  points_tausx_d,              &
  points_tausy_d,              &
  points_taus_scale_d,         &
  points_orog_slope_d,         &
  points_orog_anis_d,          &
  points_orog_dir_d

REAL, INTENT(IN)  ::                            &!block for intent(in)
  theta(row_length,rows,levels),                &!primary theta field (K)
  rho(pdims_s%i_start:pdims_s%i_end,            &!density *r*r (kg/m)
      pdims_s%j_start:pdims_s%j_end,            &!
      pdims_s%k_start:pdims_s%k_end),           &!
  delta_lambda,                                 &!spacing between pts - i dir
  delta_phi,                                    &!spacing between pts - j dir
  true_latitude(row_length, rows),              &!latitude
  timestep,                                     &!timestep (s)
  kay,                                          &!sfc stress const [4A only]
  frc,                                          &!critical Froude number
  fsat,                                         &!crit Fr num - wave breaking
  sd_orog_land(land_points),                    &!standard deviation orog(m)
  orog_grad_xx_land(land_points),               &!(dh/dx)**2 gradient orog
  orog_grad_xy_land(land_points),               &!(dh/dx)(dh/dy) grad orog
  orog_grad_yy_land(land_points)                 !(dh/dy)**2 gradient orog


REAL, INTENT(IN) ::    &!block for intent(in) 5a scheme only
  Gsharp,              &!function of mtn sharpness
  fbcd,                &!flow blocking drag coefficient
  nsigma                !Scaling factor for sd_orog_land

LOGICAL, INTENT(IN) :: &!intent(in):
  at_extremity(4)       ! indicates if this processor is at north,
                        ! south, east or west of the processor grid

LOGICAL, INTENT(IN) :: &!intent(in) 4a scheme only
  l_taus_scale,        &!if true - sfc stress depends on the low level Fr num
  l_fix_gwsatn,        &!if true - invoke minor bug fixes in gwsatn
  l_gwd_40km            !if true - don't apply GWD above 40km

LOGICAL, INTENT(IN) :: &!intent(in) 5a scheme only
  l_dynbeta,           &!if true:dynamically adjusting beta (group vel angle)
  l_nonhydro,          &!if true:use nonhydro scheme
  l_smooth,            &!if true:lambda_z smoothing of acc
  l_gw_heating(3)       !Calculate heating tendency
                        !l_gw_heating(1)=.true.: due to flow blocking drag
                        !l_gw_heating(2)=.true.: due to gravity wave drag

! Below are the stash flags for calculating diagnostics:
LOGICAL, INTENT(IN) :: &!intent(in) for diags
  stress_ud_on,        &!u stress
  stress_ud_p_on,      &!u stress p points
  stress_vd_on,        &!v stress
  stress_ud_satn_on,   &!u satn stress
  stress_vd_satn_on,   &!v satn stress
  stress_ud_wake_on,   &!u wake stress
  stress_vd_wake_on,   &!v wake stress
  du_dt_satn_on,       &!u accel (saturation)
  du_dt_satn_p_on,     &!u accel (saturation) on p points
  dv_dt_satn_on,       &!v accel (saturation)
  du_dt_wake_on,       &!u accel blocked flow
  dv_dt_wake_on,       &!v accel blocked flow
  u_s_d_on,            &!u_s_d diag switch
  v_s_d_on,            &!v_s_d diag switch
  nsq_s_d_on,          &!nsq_s_d diag switch
  fr_d_on,             &!fr_d switch
  bld_d_on,            &!bld_d switch
  bldt_d_on,           &!bldt_d switch
  num_lim_d_on,        &!num_lim_d switch (4a scheme only)
  num_fac_d_on,        &!num_fac_d switch (4a scheme only)
  tausx_d_on,          &!tausx_d switch
  tausy_d_on,          &!tausy_d switch
  taus_scale_d_on,     &!taus_scale_d switch (4a scheme only)
  orog_slope_d_on,     &!orog_slope_d switch (5a scheme only)
  orog_anis_d_on,      &!orog_anis_d switch (5a scheme only)
  orog_dir_d_on         !orog_dir_d switch (5a scheme only)

!----------------------
! Intent(inout) variables
!----------------------
REAL, INTENT(INOUT) ::                &!block for intent(inout)
  u(udims_s%i_start:udims_s%i_end,    & ! primary u field (ms**-1)
    udims_s%j_start:udims_s%j_end,    &
    udims_s%k_start:udims_s%k_end),   &
  v(vdims_s%i_start:vdims_s%i_end,    & ! primary v field (ms**-1)
    vdims_s%j_start:vdims_s%j_end,    &
    vdims_s%k_start:vdims_s%k_end),   &
  T_inc(tdims%i_start:tdims%i_end,    & !Temperature increment
        tdims%j_start:tdims%j_end,    &
                    1:tdims%k_end)

REAL, INTENT(INOUT) ::                &!intent(inout):
  r_u(udims_s%i_start:udims_s%i_end,  & !u wind increment diagnostic
      udims_s%j_start:udims_s%j_end,  &
      udims_s%k_start:udims_s%k_end), &
  r_v(vdims_s%i_start:vdims_s%i_end,  & !v wind increment diagnostic
      vdims_s%j_start:vdims_s%j_end,  &
      vdims_s%k_start:vdims_s%k_end)

!----------------------
! Intent(out) variables
!----------------------
INTEGER, INTENT(OUT)  ::          &!block for intent(in)
  iret                             ! return code : iret=0 normal exit

REAL, INTENT(OUT) ::                         &!block for intent(out)
  stress_ud(row_length,u_rows,0:levels),     &!u total stress
  stress_vd(row_length,v_rows,0:levels),     &!v total stress
  stress_ud_satn(row_length,u_rows,0:levels),&!u satn stress
  stress_vd_satn(row_length,v_rows,0:levels),&!v satn stress
  stress_ud_wake(row_length,u_rows,0:levels),&!u wake stress
  stress_vd_wake(row_length,v_rows,0:levels),&!v wake stress
  du_dt_satn(row_length,u_rows,levels),      &!u accel (saturation)
  dv_dt_satn(row_length,v_rows,levels),      &!v accel (saturation)
  du_dt_wake(row_length,u_rows,levels),      &!u accel (blocked flow)
  dv_dt_wake(row_length,v_rows,levels)        !v accel (blocked flow)

REAL              ::            &!block for 4a only (do not use INOUT)
  num_lim_d (row_length,  rows),&!% of time numerical limiter used (4a only)
  num_fac_d (row_length,  rows),&!% redn. of flow-blocking stress (4a only)
  taus_scale_d(row_length,rows)  !Factor surface stress scaled by (4a only)

REAL, INTENT(OUT) ::            &!block for intent(out)
  u_s_d  (row_length,  rows),   &!u_s diag at theta pts
  v_s_d  (row_length,  rows),   &!v_s diag at theta pts
  nsq_s_d(row_length,  rows),   &!n_s(not nsq) diag at theta pts
  fr_d   (row_length,  rows),   &!Fr diag at theta pts
  bld_d  (row_length,  rows),   &!blocked layer depth at theta pts
  bldt_d (row_length,  rows),   &!% of time blocked layer diagnosed
                                !after numerical limiter invoked
  tausx_d(row_length,u_rows),   &!x-component of surface stress
  tausy_d(row_length,v_rows),   &!y-component of surface stress
  orog_slope_d(row_length,rows),&!Slope of sub-grid orography (5a only)
  orog_anis_d(row_length,rows), &!Anisotropy of sub-grid orography (5a only)
  orog_dir_d(row_length,rows)    !Orientation of sub-grid orography (5a only)


!--------------------------------------------------------------------
! Local variables
!--------------------------------------------------------------------

INTEGER  ::           &! block for local variables
  k_top(land_points), &! first model level above mountain tops
  i,j,k,l              ! loop counters in routine

! Segmentation is an attempt to 
! a) improve memory accesses by utilising cache better
! b) provide a means by which OpenMP parallelisation may be introduced
! It needs a few control variables
INTEGER :: seg_start                  ! start of segment
INTEGER :: seg_end                    ! end of segment
INTEGER :: seg_len                    ! length of segment

! All the seg_ variables below are for the optional diagnostics
INTEGER :: seg_nsq_s_d_start, seg_nsq_s_d_end, seg_nsq_s_d_len
INTEGER :: seg_u_s_d_start, seg_u_s_d_end, seg_u_s_d_len
INTEGER :: seg_v_s_d_start, seg_v_s_d_end, seg_v_s_d_len
INTEGER :: seg_du_dt_wake_start, seg_du_dt_wake_end, seg_du_dt_wake_len
INTEGER :: seg_dv_dt_wake_start, seg_dv_dt_wake_end, seg_dv_dt_wake_len
INTEGER :: seg_stress_ud_start, seg_stress_ud_end, seg_stress_ud_len
INTEGER :: seg_stress_vd_start, seg_stress_vd_end, seg_stress_vd_len
INTEGER :: seg_stress_ud_wake_start, seg_stress_ud_wake_end
INTEGER :: seg_stress_ud_wake_len
INTEGER :: seg_stress_vd_wake_start, seg_stress_vd_wake_end
INTEGER :: seg_stress_vd_wake_len
INTEGER :: seg_fr_d_start, seg_fr_d_end, seg_fr_d_len
INTEGER :: seg_bld_d_start, seg_bld_d_end, seg_bld_d_len
INTEGER :: seg_bldt_d_start, seg_bldt_d_end, seg_bldt_d_len
INTEGER :: seg_tausx_d_start, seg_tausx_d_end, seg_tausx_d_len
INTEGER :: seg_tausy_d_start, seg_tausy_d_end, seg_tausy_d_len
INTEGER :: seg_du_dt_satn_start, seg_du_dt_satn_end, seg_du_dt_satn_len
INTEGER :: seg_dv_dt_satn_start, seg_dv_dt_satn_end, seg_dv_dt_satn_len
INTEGER :: seg_stress_ud_satn_start, seg_stress_ud_satn_end
INTEGER :: seg_stress_ud_satn_len
INTEGER :: seg_stress_vd_satn_start, seg_stress_vd_satn_end
INTEGER :: seg_stress_vd_satn_len

REAL      ::                                  &!block for work arrays
  work_u(row_length,rows,levels),             &
  work_v(row_length,rows,levels),             &
  work_on_u_grid(udims%i_start:udims%i_end,   &
                 udims%j_start:udims%j_end,   &
                 udims%k_start:udims%k_end),  &
  work_on_v_grid(vdims%i_start:vdims%i_end,   &
                 vdims%j_start:vdims%j_end,   &
                 vdims%k_start:vdims%k_end)

! The scatter of accelerations and of diagnostics and the
! interpolation from p to u,v staggering on 'c' grid is performed
! using two super-arrays, one for u fields and one for v fields. The
! halo exchange is performed twice with overlapping of communications
! with computations. The various fields are separated in the
! super-arrays by different levels.

REAL, ALLOCATABLE :: work_halo_u(:,:,:), work_halo_v(:,:,:)
INTEGER           :: k_total_u, k_total_v, k_w
INTEGER           :: k_du_dt_start, k_du_dt_end
INTEGER           :: k_dv_dt_start, k_dv_dt_end
INTEGER           :: k_stress_ud_start, k_stress_ud_end
INTEGER           :: k_stress_vd_start, k_stress_vd_end
INTEGER           :: k_stress_ud_satn_start, k_stress_ud_satn_end
INTEGER           :: k_stress_vd_satn_start, k_stress_vd_satn_end
INTEGER           :: k_stress_ud_wake_start, k_stress_ud_wake_end
INTEGER           :: k_stress_vd_wake_start, k_stress_vd_wake_end
INTEGER           :: k_du_dt_satn_start, k_du_dt_satn_end
INTEGER           :: k_dv_dt_satn_start, k_dv_dt_satn_end
INTEGER           :: k_du_dt_wake_start, k_du_dt_wake_end
INTEGER           :: k_dv_dt_wake_start, k_dv_dt_wake_end
INTEGER           :: k_tausx_d_start, k_tausx_d_end
INTEGER           :: k_tausy_d_start, k_tausy_d_end


REAL ::               &!currently local variables
  slope(land_points), &!sso params - slope
  anis(land_points),  &!sso params - anisotropy
  mtdir(land_points), &!sso params - angle of major axis
  mt_high(land_points) !sso params - height

REAL    ::                        &!block for local variables
  up_land(land_points,levels),    &!interpolated u on theta grid
  vp_land(land_points,levels),    &!interpolated V on theta grid
  theta_land(land_points,levels), &!land theta field on theta levels
  r_rho_levels_land               &!land field of heights of
  (land_points,levels),           &!rho levels above z=0
  r_theta_levels_land             &!land field of heights of
  (land_points,levels),           &!theta levels above z=0
  rho_land(land_points,levels),   &!density at land points
  du_dt(land_points,levels),      &!total mountain drag du/dt on land/theta
  dv_dt(land_points,levels),      &!total mountain drag dv/dt on land/theta
  dt_dt(land_points,levels),      &!temperature tendency on land/theta
  nsq(land_points,levels),        &!buoyancy frequency
  latitude_land(land_points),     &!true latitude on land points
  banis(land_points),             &!const (function of anisotropy)
  canis(land_points),             &!const (function of anisotropy)
  ulow(land_points),              &!u averaged from 0.5h to h
  vlow(land_points),              &!v averaged from 0.5h to h
  modu(land_points),              &!modulus of horizontal wind(sqrt(u^2+v^2))
  nlow(land_points),              &!N bulk averaged from 0.5h to h
  psilow(land_points),            &!psi averaged from 0.5h to h
  rholow(land_points),            &!rho averaged from 0.5h to h
  psi1(land_points),              &!atan(vlow/ulow)
  zb(land_points)                  !depth of flow blocking layer

! The next set of variables are the segmented versions used in the main
! subroutine calls
REAL, ALLOCATABLE :: du_dt_s(:,:)
REAL, ALLOCATABLE :: dv_dt_s(:,:)
REAL, ALLOCATABLE :: nsq_s(:,:)
REAL, ALLOCATABLE :: up_land_s(:,:)
REAL, ALLOCATABLE :: vp_land_s(:,:)
REAL, ALLOCATABLE :: rho_land_s(:,:)
REAL, ALLOCATABLE :: theta_land_s(:,:)
REAL, ALLOCATABLE :: r_rho_levels_land_s(:,:)
REAL, ALLOCATABLE :: r_theta_levels_land_s(:,:)
REAL, ALLOCATABLE :: dt_dt_s(:,:)
REAL, ALLOCATABLE :: du_dt_wake_land_s(:,:)
REAL, ALLOCATABLE :: dv_dt_wake_land_s(:,:)
REAL, ALLOCATABLE :: stress_ud_land_s(:,:)
REAL, ALLOCATABLE :: stress_vd_land_s(:,:)
REAL, ALLOCATABLE :: stress_ud_wake_land_s(:,:)
REAL, ALLOCATABLE :: stress_vd_wake_land_s(:,:)
REAL, ALLOCATABLE :: du_dt_satn_land_s(:,:)
REAL, ALLOCATABLE :: dv_dt_satn_land_s(:,:)
REAL, ALLOCATABLE :: stress_ud_satn_land_s(:,:)
REAL, ALLOCATABLE :: stress_vd_satn_land_s(:,:)

LOGICAL   ::          &
   l_drag(land_points) !whether point has a non-zero stress or not

!--------------------------------------------------------------------
! Local diagnostic variables
!--------------------------------------------------------------------
!
! Land points arrays below are for the total GWD stress and
! its 4 individual components. (x and y components for each)
!
REAL  ::                                                   &
  stress_ud_land      ( points_stress_ud      , 0:levels ),&
  stress_vd_land      ( points_stress_vd      , 0:levels ),&
  stress_ud_satn_land ( points_stress_ud_satn , 0:levels ),&
  stress_vd_satn_land ( points_stress_vd_satn , 0:levels ),&
  stress_ud_wake_land ( points_stress_ud_wake , 0:levels ),&
  stress_vd_wake_land ( points_stress_vd_wake , 0:levels )
!
!  Land point arrays below are for the 4 individual components of
!  the GWD wind increment.  (x and y components for each)
!
REAL  ::                                         &
  du_dt_satn_land ( points_du_dt_satn , levels ),&
  dv_dt_satn_land ( points_dv_dt_satn , levels ),&
  du_dt_wake_land ( points_du_dt_wake , levels ),&
  dv_dt_wake_land ( points_dv_dt_wake , levels )

!
!  Land point arrays below are for the 9 GWD 'surface' diagnostics.
!
REAL  ::                                 &
  u_s_d_land(points_u_s_d),              &
  v_s_d_land(points_v_s_d),              &
  nsq_s_d_land(points_nsq_s_d),          &
  fr_d_land(points_fr_d),                &
  bld_d_land(points_bld_d),              &
  bldt_d_land(points_bldt_d),            &
  tausx_d_land(points_tausx_d),          &
  tausy_d_land(points_tausy_d)

TYPE (array_dims)  :: gwd_udims          &  ! gwd 0 to top udims for diags
                     ,gwd_vdims          &  ! gwd 0 to top vdims for diags
                     ,gwd_udims1         &  ! gwd single level udims
                     ,gwd_vdims1            ! gwd single level vdims

LOGICAL, PARAMETER   ::  calc_k_levels = .TRUE. 
LOGICAL, PARAMETER   ::  l_vector = .FALSE.   ! scalar fields for swapbounds


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='G_WAVE_5A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------
! set up local dims for gwd diags as required by eg_v_at_poles
! -----------------------------------------------------------------

gwd_udims                   = udims
gwd_udims%k_start           = 0
gwd_vdims                   = vdims
gwd_vdims%k_start           = 0

gwd_udims1                  = udims
gwd_udims1%k_start          = 1
gwd_udims1%k_end            = 1
gwd_vdims1                  = vdims
gwd_vdims1%k_start          = 1
gwd_vdims1%k_end            = 1

!------------------------------------------------------------------
! 0.1 set up the start and end levels for each fields in the
! super-arrays
! -----------------------------------------------------------------

! 0.1.1 u fields

k = 1
CALL set_levels(calc_k_levels, k, levels, k_du_dt_start, k_du_dt_end) 
CALL set_levels((stress_ud_on .OR. stress_ud_p_on), k, levels+1, & 
                k_stress_ud_start, k_stress_ud_end) 
CALL set_levels(stress_ud_satn_on, k, levels+1,                  &
                k_stress_ud_satn_start, k_stress_ud_satn_end) 
CALL set_levels(stress_ud_wake_on, k, levels+1,                  &
                k_stress_ud_wake_start, k_stress_ud_wake_end) 
CALL set_levels((du_dt_satn_on .OR. du_dt_satn_p_on), k, levels, &
                k_du_dt_satn_start, k_du_dt_satn_end) 
CALL set_levels(du_dt_wake_on, k, levels,                        &
                k_du_dt_wake_start, k_du_dt_wake_end) 
CALL set_levels(tausx_d_on, k, 1, k_tausx_d_start, k_tausx_d_end) 

! Get the total number of levels in the u super-array.
k_total_u = k-1

! 0.1.2 v fields

k = 1
CALL set_levels(calc_k_levels, k, levels, k_dv_dt_start, k_dv_dt_end) 
CALL set_levels(stress_vd_on,                       k, levels+1, & 
                k_stress_vd_start, k_stress_vd_end) 
CALL set_levels(stress_vd_satn_on, k, levels+1,                  &
                k_stress_vd_satn_start, k_stress_vd_satn_end) 
CALL set_levels(stress_vd_wake_on, k, levels+1,                  &
                k_stress_vd_wake_start, k_stress_vd_wake_end) 
CALL set_levels(dv_dt_satn_on,                        k, levels, &
                k_dv_dt_satn_start, k_dv_dt_satn_end) 
CALL set_levels(dv_dt_wake_on, k, levels,                        &
                k_dv_dt_wake_start, k_dv_dt_wake_end) 
CALL set_levels(tausy_d_on, k, 1, k_tausy_d_start, k_tausy_d_end) 

! Get the total number of levels in the v super-array.
k_total_v = k-1


! 0.1.3 allocate super-arrays

ALLOCATE(work_halo_u(pdims_s%i_start:pdims_s%i_end,            &
                     pdims_s%j_start:pdims_s%j_end,            &
                     1:k_total_u))

ALLOCATE(work_halo_v(pdims_s%i_start:pdims_s%i_end,            &
                     pdims_s%j_start:pdims_s%j_end,            &
                     1:k_total_v))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,l,j,i,k_w,                     &
!$OMP    seg_start, seg_end, seg_len,                                   &
!$OMP    seg_u_s_d_start, seg_u_s_d_end, seg_u_s_d_len,                 &    
!$OMP    seg_v_s_d_start, seg_v_s_d_end, seg_v_s_d_len,                 &    
!$OMP    seg_du_dt_wake_start, seg_du_dt_wake_end, seg_du_dt_wake_len,  &    
!$OMP    seg_dv_dt_wake_start, seg_dv_dt_wake_end, seg_dv_dt_wake_len,  &    
!$OMP    seg_stress_ud_start, seg_stress_ud_end, seg_stress_ud_len,     &    
!$OMP    seg_stress_vd_start, seg_stress_vd_end, seg_stress_vd_len,     &    
!$OMP    seg_stress_ud_wake_start, seg_stress_ud_wake_end,              &
!$OMP    seg_stress_ud_wake_len,                                        &    
!$OMP    seg_stress_vd_wake_start, seg_stress_vd_wake_end,              &
!$OMP    seg_stress_vd_wake_len,                                        &    
!$OMP    seg_fr_d_start, seg_fr_d_end, seg_fr_d_len,                    &    
!$OMP    seg_bld_d_start, seg_bld_d_end, seg_bld_d_len,                 &    
!$OMP    seg_bldt_d_start, seg_bldt_d_end, seg_bldt_d_len,              &    
!$OMP    seg_tausx_d_start, seg_tausx_d_end, seg_tausx_d_len,           &    
!$OMP    seg_tausy_d_start, seg_tausy_d_end, seg_tausy_d_len,           &    
!$OMP    seg_du_dt_satn_start, seg_du_dt_satn_end, seg_du_dt_satn_len,  &    
!$OMP    seg_dv_dt_satn_start, seg_dv_dt_satn_end, seg_dv_dt_satn_len,  &    
!$OMP    seg_stress_ud_satn_start, seg_stress_ud_satn_end,              &
!$OMP    seg_stress_ud_satn_len,                                        &    
!$OMP    seg_stress_vd_satn_start, seg_stress_vd_satn_end,              &
!$OMP    seg_stress_vd_satn_len,                                        &
!$OMP    seg_nsq_s_d_start, seg_nsq_s_d_end, seg_nsq_s_d_len,           &
!$OMP    du_dt_s, dv_dt_s, nsq_s, up_land_s, vp_land_s, rho_land_s,     &
!$OMP    theta_land_s, r_rho_levels_land_s, r_theta_levels_land_s,      &
!$OMP    dt_dt_s, du_dt_wake_land_s, dv_dt_wake_land_s,                 &
!$OMP    stress_ud_land_s, stress_vd_land_s, stress_ud_wake_land_s,     &
!$OMP    stress_vd_wake_land_s, du_dt_satn_land_s, dv_dt_satn_land_s,   &
!$OMP    stress_ud_satn_land_s, stress_vd_satn_land_s)

!------------------------------------------------------------------
!l    1.1 interpolate winds to p/theta-grid
!------------------------------------------------------------------

CALL u_to_p(u,                                                    &
                  udims_s%i_start,udims_s%i_end,                  &
                  udims_s%j_start,udims_s%j_end,                  &
                  pdims%i_start,pdims%i_end,                      &
                  pdims%j_start,pdims%j_end,                      &
                  levels, at_extremity,work_u)
!
CALL v_to_p(v,                                                    &
                  vdims_s%i_start,vdims_s%i_end,                  &
                  vdims_s%j_start,vdims_s%j_end,                  &
                  pdims%i_start,pdims%i_end,                      &
                  pdims%j_start,pdims%j_end,                      &
                  levels, at_extremity,work_v)
!$OMP BARRIER

!------------------------------------------------------------------
!l    1.2  gather winds at land points
!------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,levels
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    up_land(l,k) =work_u(i,j,k)
    vp_land(l,k) =work_v(i,j,k)
  END DO
END DO
!$OMP END DO

!------------------------------------------------------------------
!l    1.3  gather theta, rho and heights at land points
!------------------------------------------------------------------

IF (land_points  >   0) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      j = (land_index(l)-1)/row_length + 1
      i = land_index(l) - (j-1)*row_length
      r_rho_levels_land(l,k)  = r_rho_levels(i,j,k) -               &
                                r_theta_levels(i,j,0)
      r_theta_levels_land(l,k)= r_theta_levels(i,j,k) -             &
                                r_theta_levels(i,j,0)
      rho_land(l,k)           = rho(i,j,k)*recip_a2
      theta_land(l,k)         = theta(i,j,k)
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
      j = (land_index(l)-1)/row_length + 1
      i = land_index(l) - (j-1)*row_length
      latitude_land(l)       = true_latitude(i,j)
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(DYNAMIC)
  DO l = 1, land_points, gw_seg_size

    ! Calculate the addressing for segments. Note that for diagnostics, 
    ! this is long and tedious, but essentially trivial.
    seg_start = l
    seg_end   = seg_start + gw_seg_size - 1
    IF (seg_end > land_points) seg_end = land_points
    seg_len = seg_end - seg_start + 1

    CALL set_segments(points_u_s_d, seg_start, seg_end, seg_len,              &
             seg_u_s_d_start, seg_u_s_d_end, seg_u_s_d_len)

    CALL set_segments(points_v_s_d, seg_start, seg_end, seg_len,              &
             seg_v_s_d_start, seg_v_s_d_end, seg_v_s_d_len)

    CALL set_segments(points_nsq_s_d, seg_start, seg_end, seg_len,            &
             seg_nsq_s_d_start, seg_nsq_s_d_end, seg_nsq_s_d_len)

    CALL set_segments(points_du_dt_wake, seg_start, seg_end, seg_len,         &
             seg_du_dt_wake_start, seg_du_dt_wake_end, seg_du_dt_wake_len)

    CALL set_segments(points_dv_dt_wake, seg_start, seg_end, seg_len,         &
             seg_dv_dt_wake_start, seg_dv_dt_wake_end, seg_dv_dt_wake_len)
      
    CALL set_segments(points_stress_ud, seg_start, seg_end, seg_len,          &
             seg_stress_ud_start, seg_stress_ud_end, seg_stress_ud_len)

    CALL set_segments(points_stress_vd, seg_start, seg_end, seg_len,          &
             seg_stress_vd_start, seg_stress_vd_end, seg_stress_vd_len)

    CALL set_segments(points_stress_ud_wake, seg_start, seg_end, seg_len,     &
       seg_stress_ud_wake_start, seg_stress_ud_wake_end, seg_stress_ud_wake_len)

    CALL set_segments(points_stress_vd_wake, seg_start, seg_end, seg_len,     &
       seg_stress_vd_wake_start, seg_stress_vd_wake_end, seg_stress_vd_wake_len)

    CALL set_segments(points_fr_d, seg_start, seg_end, seg_len,               &
             seg_fr_d_start, seg_fr_d_end, seg_fr_d_len)

    CALL set_segments(points_bld_d, seg_start, seg_end, seg_len,              &
             seg_bld_d_start, seg_bld_d_end, seg_bld_d_len)

    CALL set_segments(points_bldt_d, seg_start, seg_end, seg_len,             &
             seg_bldt_d_start, seg_bldt_d_end, seg_bldt_d_len)

    CALL set_segments(points_tausx_d, seg_start, seg_end, seg_len,            &
             seg_tausx_d_start, seg_tausx_d_end, seg_tausx_d_len)

    CALL set_segments(points_tausy_d, seg_start, seg_end, seg_len,            &
             seg_tausy_d_start, seg_tausy_d_end, seg_tausy_d_len)

    CALL set_segments(points_du_dt_satn, seg_start, seg_end, seg_len,         &
             seg_du_dt_satn_start, seg_du_dt_satn_end, seg_du_dt_satn_len)

    CALL set_segments(points_dv_dt_satn, seg_start, seg_end, seg_len,         &
             seg_dv_dt_satn_start, seg_dv_dt_satn_end, seg_dv_dt_satn_len)

    CALL set_segments(points_stress_ud_satn, seg_start, seg_end, seg_len,     &
      seg_stress_ud_satn_start, seg_stress_ud_satn_end, seg_stress_ud_satn_len)

    CALL set_segments(points_stress_vd_satn, seg_start, seg_end, seg_len,     &
      seg_stress_vd_satn_start, seg_stress_vd_satn_end, seg_stress_vd_satn_len)

! Allocate the segmented variables
    ALLOCATE(du_dt_s(seg_start:seg_end,levels))
    ALLOCATE(dv_dt_s(seg_start:seg_end,levels))
    ALLOCATE(nsq_s(seg_start:seg_end,levels))
    ALLOCATE(up_land_s(seg_start:seg_end,levels))
    ALLOCATE(vp_land_s(seg_start:seg_end,levels))
    ALLOCATE(rho_land_s(seg_start:seg_end,levels))
    ALLOCATE(theta_land_s(seg_start:seg_end,levels))
    ALLOCATE(r_rho_levels_land_s(seg_start:seg_end,levels))
    ALLOCATE(r_theta_levels_land_s(seg_start:seg_end,levels))
    ALLOCATE(dt_dt_s(seg_start:seg_end,levels))
    ALLOCATE(du_dt_wake_land_s(seg_du_dt_wake_start:seg_du_dt_wake_end,levels))
    ALLOCATE(dv_dt_wake_land_s(seg_dv_dt_wake_start:seg_dv_dt_wake_end,levels))
    ALLOCATE(stress_ud_land_s(seg_stress_ud_start:seg_stress_ud_end,0:levels))
    ALLOCATE(stress_vd_land_s(seg_stress_vd_start:seg_stress_vd_end,0:levels))
    ALLOCATE(stress_ud_wake_land_s(seg_stress_ud_wake_start:             &
                                   seg_stress_ud_wake_end,0:levels))
    ALLOCATE(stress_vd_wake_land_s(seg_stress_vd_wake_start:             &
                                   seg_stress_vd_wake_end,0:levels))
    ALLOCATE(du_dt_satn_land_s(seg_du_dt_satn_start:seg_du_dt_satn_end,levels))
    ALLOCATE(dv_dt_satn_land_s(seg_dv_dt_satn_start:seg_dv_dt_satn_end,levels))
    ALLOCATE(stress_ud_satn_land_s(seg_stress_ud_satn_start:             &
                                   seg_stress_ud_satn_end,0:levels))
    ALLOCATE(stress_vd_satn_land_s(seg_stress_vd_satn_start:             &
                                   seg_stress_vd_satn_end,0:levels))

! Where segmented data needs to be input - copy it in.
    DO j=1,levels
      DO i=seg_start,seg_end
        up_land_s(i,j)             = up_land(i,j)
        vp_land_s(i,j)             = vp_land(i,j)
        rho_land_s(i,j)            = rho_land(i,j)
        theta_land_s(i,j)          = theta_land(i,j)
        r_rho_levels_land_s(i,j)   = r_rho_levels_land(i,j)
        r_theta_levels_land_s(i,j) = r_theta_levels_land(i,j)
      END DO
    END DO

! Note the next three routines use the segmented versions of variables, or
! pass in segments via array notation.

    CALL gw_setup(levels,seg_len,                                              &
                up_land_s,vp_land_s,                                           &
                rho_land_s,theta_land_s,                                       &
                nsq_s,ulow(seg_start:seg_end),                                 &
                vlow(seg_start:seg_end),rholow(seg_start:seg_end),             &
                nlow(seg_start:seg_end),psilow(seg_start:seg_end),             &
                psi1(seg_start:seg_end),modu(seg_start:seg_end),               &
                r_rho_levels_land_s,                                           &
                r_theta_levels_land_s,                                         &
                sd_orog_land(seg_start:seg_end),                               &
                orog_grad_xx_land(seg_start:seg_end),                          &
                orog_grad_xy_land(seg_start:seg_end),                          &
                orog_grad_yy_land(seg_start:seg_end),                          &
                mt_high(seg_start:seg_end),                                    &
                slope(seg_start:seg_end),anis(seg_start:seg_end),              &
                banis(seg_start:seg_end),canis(seg_start:seg_end),             &
                mtdir(seg_start:seg_end), k_top(seg_start:seg_end),            &
                l_drag(seg_start:seg_end),nsigma,                              &
  !diagnostics
                u_s_d_land(seg_u_s_d_start:seg_u_s_d_end),                     &
                u_s_d_on,seg_u_s_d_len,                                        &
                v_s_d_land(seg_v_s_d_start:seg_v_s_d_end),                     &
                v_s_d_on,seg_v_s_d_len,                                        &
                nsq_s_d_land(seg_nsq_s_d_start:seg_nsq_s_d_end),               &
                nsq_s_d_on, seg_nsq_s_d_len)

    CALL gw_block(levels,seg_len,timestep,up_land_s,                           &
       vp_land_s,rho_land_s,                                                   &
       nsq_s,ulow(seg_start:seg_end),                                          &
       vlow(seg_start:seg_end),rholow(seg_start:seg_end),                      &
       psilow(seg_start:seg_end),modu(seg_start:seg_end),                      &
       r_rho_levels_land_s,                                                    &
       r_theta_levels_land_s,mt_high(seg_start:seg_end),                       &
       sd_orog_land(seg_start:seg_end),slope(seg_start:seg_end),               &
       anis(seg_start:seg_end),mtdir(seg_start:seg_end),zb(seg_start:seg_end), &
       banis(seg_start:seg_end),canis(seg_start:seg_end),                      &
       du_dt_s,dv_dt_s,                                                        &
       dt_dt_s,fbcd,frc,l_drag(seg_start:seg_end),                             &
       l_gw_heating(1),                                                        &
  !diagnostics
       du_dt_wake_land_s, seg_du_dt_wake_len,du_dt_wake_on,                    &
       dv_dt_wake_land_s, seg_dv_dt_wake_len,dv_dt_wake_on,                    &
       stress_ud_land_s, stress_ud_on,seg_stress_ud_len,                       &
       stress_ud_p_on,                                                         &
       stress_vd_land_s, stress_vd_on,seg_stress_vd_len,                       &
       stress_ud_wake_land, seg_stress_ud_wake_len, stress_ud_wake_on,         &
       stress_vd_wake_land, seg_stress_vd_wake_len, stress_vd_wake_on,         &
       fr_d_land(seg_fr_d_start:seg_fr_d_end),fr_d_on,seg_fr_d_len,            &
       bld_d_land(seg_bld_d_start:seg_bld_d_end),bld_d_on,                     &
       seg_bld_d_len,                                                          &
       bldt_d_land(seg_bldt_d_start:seg_bldt_d_end),bldt_d_on,                 &
       seg_bldt_d_len,                                                         &
       tausx_d_land(seg_tausx_d_start:seg_tausx_d_end),tausx_d_on,             &
       seg_tausx_d_len,                                                        &
       tausy_d_land(seg_tausy_d_start:seg_tausy_d_end),tausy_d_on,             &
       seg_tausy_d_len)

  !------------------------------------------------------------------
  !l    3. calculate stress profile and accelerations,
  !l       CALL gw_vert
  !------------------------------------------------------------------
    CALL gw_wave(levels,seg_len,up_land_s,                                     &
       vp_land_s,rho_land_s,                                                   &
       nsq_s,ulow(seg_start:seg_end),                                          &
       vlow(seg_start:seg_end),rholow(seg_start:seg_end),                      &
       psi1(seg_start:seg_end),psilow(seg_start:seg_end),                      &
       nlow(seg_start:seg_end),modu(seg_start:seg_end),                        &
       k_top(seg_start:seg_end),r_rho_levels_land_s,                           &
       r_theta_levels_land_s,delta_lambda,delta_phi,                           &
       latitude_land(seg_start:seg_end),mt_high(seg_start:seg_end),            &
       sd_orog_land(seg_start:seg_end),slope(seg_start:seg_end),               &
       zb(seg_start:seg_end),banis(seg_start:seg_end),canis(seg_start:seg_end),&
       du_dt_s,dv_dt_s,                                                        &
       dt_dt_s,timestep,l_dynbeta,l_nonhydro,l_smooth,fsat,                    &
       Gsharp,l_drag(seg_start:seg_end),l_gw_heating(2),                       &
  !diagnostics
       du_dt_satn_land_s, seg_du_dt_satn_len,du_dt_satn_on,du_dt_satn_p_on,    &
       dv_dt_satn_land_s, seg_dv_dt_satn_len,dv_dt_satn_on,                    &
       stress_ud_land_s, seg_stress_ud_len, stress_ud_on, stress_ud_p_on,      &
       stress_vd_land_s, seg_stress_vd_len, stress_vd_on,                      &
       stress_ud_satn_land_s, seg_stress_ud_satn_len, stress_ud_satn_on,       &
       stress_vd_satn_land_s, seg_stress_vd_satn_len, stress_vd_satn_on,       &
       tausx_d_land(seg_tausx_d_start:seg_tausx_d_end), tausx_d_on,            &
       seg_tausx_d_len,                                                        &
       tausy_d_land(seg_tausy_d_start:seg_tausy_d_end), tausy_d_on,            &
       seg_tausy_d_len)

! For segmented data that is output - copy it out
        DO j=1,levels
          DO i=seg_start,seg_end
            du_dt(i,j) = du_dt_s(i,j)
            dv_dt(i,j) = dv_dt_s(i,j)
            dt_dt(i,j) = dt_dt_s(i,j)
            nsq(i,j)   = nsq_s(i,j)
          END DO
        END DO

        IF (du_dt_wake_on) THEN
          DO j=1,levels
            DO i=seg_du_dt_wake_start,seg_du_dt_wake_end
              du_dt_wake_land(i,j) = du_dt_wake_land_s(i,j)
            END DO
          END DO
        END IF

        IF (dv_dt_wake_on) THEN
          DO j=1,levels
            DO i=seg_dv_dt_wake_start,seg_dv_dt_wake_end
              dv_dt_wake_land(i,j) = dv_dt_wake_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_ud_on .OR. stress_ud_p_on) THEN 
          DO j=0,levels
            DO i=seg_stress_ud_start,seg_stress_ud_end
              stress_ud_land(i,j) = stress_ud_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_vd_on) THEN
          DO j=0,levels
            DO i=seg_stress_vd_start,seg_stress_vd_end
              stress_vd_land(i,j) = stress_vd_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_ud_wake_on) THEN 
          DO j=0,levels
            DO i=seg_stress_ud_wake_start,seg_stress_ud_wake_end
              stress_ud_wake_land(i,j) = stress_ud_wake_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_vd_wake_on) THEN 
          DO j=0,levels
            DO i=seg_stress_vd_wake_start,seg_stress_vd_wake_end
              stress_vd_wake_land(i,j) = stress_vd_wake_land_s(i,j)
            END DO
          END DO
        END IF

        IF (du_dt_satn_on .OR. du_dt_satn_p_on) THEN
          DO j=1,levels
            DO i=seg_du_dt_satn_start, seg_du_dt_satn_end
              du_dt_satn_land(i,j) = du_dt_satn_land_s(i,j)
            END DO
          END DO
        END IF

        IF (dv_dt_satn_on) THEN
          DO j=1,levels
            DO i=seg_dv_dt_satn_start, seg_dv_dt_satn_end
              dv_dt_satn_land(i,j) = dv_dt_satn_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_ud_satn_on) THEN 
          DO j=0,levels
            DO i=seg_stress_ud_satn_start,seg_stress_ud_satn_end
              stress_ud_satn_land(i,j) = stress_ud_satn_land_s(i,j)
            END DO
          END DO
        END IF

        IF (stress_vd_satn_on) THEN 
          DO j=0,levels
            DO i=seg_stress_vd_satn_start,seg_stress_vd_satn_end
              stress_vd_satn_land(i,j) = stress_vd_satn_land_s(i,j)
            END DO
          END DO
        END IF

! And deallocate the segmented variables - next segment may be a different size
        DEALLOCATE(stress_vd_satn_land_s)
        DEALLOCATE(stress_ud_satn_land_s)
        DEALLOCATE(dv_dt_satn_land_s)
        DEALLOCATE(du_dt_satn_land_s)
        DEALLOCATE(stress_vd_wake_land_s)
        DEALLOCATE(stress_ud_wake_land_s)
        DEALLOCATE(stress_vd_land_s)
        DEALLOCATE(stress_ud_land_s)
        DEALLOCATE(dv_dt_wake_land_s)
        DEALLOCATE(du_dt_wake_land_s)
        DEALLOCATE(dt_dt_s)
        DEALLOCATE(r_theta_levels_land_s)
        DEALLOCATE(r_rho_levels_land_s)
        DEALLOCATE(theta_land_s)
        DEALLOCATE(rho_land_s)
        DEALLOCATE(vp_land_s)
        DEALLOCATE(up_land_s)
        DEALLOCATE(nsq_s)
        DEALLOCATE(dv_dt_s)
        DEALLOCATE(du_dt_s)
      
  END DO
!$OMP END DO NOWAIT

END IF ! on land_points > 0

!------------------------------------------------------------------
!l    4.1  gather accelerations and diagnostics from land points 
!l         to super-array  (u fields)
!------------------------------------------------------------------
! initialise work array: required to ensure zero values generated by
!   non-land points. 

!$OMP DO SCHEDULE(STATIC)
DO k = 1, k_total_u
  DO j = pdims_s%j_start,pdims_s%j_end
    DO i = pdims_s%i_start,pdims_s%i_end
      work_halo_u(i,j,k)=0.0
    END DO !i
  END DO  !j
END DO   !k
!$OMP END DO
! Barrier 


!$OMP DO SCHEDULE(STATIC)
DO k = 1, levels
  DO l=1,land_points
    k_w = k_du_dt_start + k - 1
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    work_halo_u(i,j,k_w) = du_dt(l,k)
  END DO ! l
END DO ! k
!$OMP END DO NOWAIT

IF (stress_ud_on .OR. stress_ud_p_on) THEN  
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, levels 
    DO l=1,land_points 
      k_w = k_stress_ud_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_u(i,j,k_w) = stress_ud_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (stress_ud_satn_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=0,levels
    DO l=1,land_points
      k_w = k_stress_ud_satn_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_u(i,j,k_w) = stress_ud_satn_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (stress_ud_wake_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=0,levels
    DO l=1,land_points
      k_w = k_stress_ud_wake_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_u(i,j,k_w) = stress_ud_wake_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (du_dt_satn_on .OR. du_dt_satn_p_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      k_w = k_du_dt_satn_start + k - 1
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_u(i,j,k_w) = du_dt_satn_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (du_dt_wake_on) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      k_w = k_du_dt_wake_start + k - 1
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_u(i,j,k_w) = du_dt_wake_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (tausx_d_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    k_w = k_tausx_d_start 
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    work_halo_u(i,j,k_w) = tausx_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

!------------------------------------------------------------------
!l    4.2  start swap bounds of the super array (u fields) 
!------------------------------------------------------------------

! The variable work_halo_u is going to be used in the "p_to_u" routine
! call later which computes the interior points and requires the halo
! information to be set. This routine inside uses the field data at j
! index only (no j+1 or j-1). Thus the North and South halo
! information is not needed. Furthermore, the "p_to_u" routine uses
! the field data at i and i+1 indices with the starting point at zero
! (udims%i_start) which is the West halo. Thus the East halo
! information is not needed.  
! Begin the swap bound of the West halo only.

CALL begin_swap_bounds(                                           &
     work_halo_u,                                                 &
     row_length, rows,                                            &
     k_total_u, off_x, off_y, fld_type_p, l_vector,               &
     do_east_arg=.FALSE.,  do_west_arg=.TRUE.,                    &
     do_south_arg=.FALSE., do_north_arg=.FALSE.,                  &
     do_corners_arg=.FALSE.)


!------------------------------------------------------------------
!l    4.3  gather accelerations and diagnostics from land points 
!l         to super-array  (v fields)
!------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i,l,k_w)                 &
!$OMP SHARED(pdims_s,levels,k_total_v,land_points,land_index,     &
!$OMP    row_length,work_halo_v,dv_dt,stress_vd_land,             &
!$OMP    stress_vd_satn_land,stress_vd_wake_land,dv_dt_satn_land, &
!$OMP    dv_dt_wake_land,tausy_d_land,k_dv_dt_start,              &
!$OMP    k_stress_vd_start,k_stress_vd_satn_start,stress_vd_on,   &
!$OMP    k_stress_vd_wake_start,k_dv_dt_satn_start, tausy_d_on,   &
!$OMP    k_dv_dt_wake_start,k_tausy_d_start,stress_vd_satn_on,    &
!$OMP    stress_vd_wake_on,dv_dt_satn_on,dv_dt_wake_on)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, k_total_v
  DO j = pdims_s%j_start,pdims_s%j_end
    DO i = pdims_s%i_start,pdims_s%i_end
      work_halo_v(i,j,k)=0.0
    END DO !i
  END DO  !j
END DO   !k
!$OMP END DO
! Barrier

!$OMP DO SCHEDULE(STATIC)
DO k= 1, levels
  DO l=1,land_points
    k_w = k_dv_dt_start + k - 1
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    work_halo_v(i,j,k_w) = dv_dt(l,k)
  END DO ! l
END DO ! k
!$OMP END DO NOWAIT

IF (stress_vd_on) THEN  
!$OMP DO SCHEDULE(STATIC)
  DO k = 0, levels 
    DO l=1,land_points
      k_w = k_stress_vd_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_v(i,j,k_w) = stress_vd_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (stress_vd_satn_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=0,levels
    DO l=1,land_points
      k_w = k_stress_vd_satn_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_v(i,j,k_w) = stress_vd_satn_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (stress_vd_wake_on) THEN  
!$OMP DO SCHEDULE(STATIC)
  DO k=0,levels
    DO l=1,land_points
      k_w = k_stress_vd_wake_start + k
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_v(i,j,k_w) = stress_vd_wake_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (dv_dt_satn_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      k_w = k_dv_dt_satn_start + k - 1
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_v(i,j,k_w) = dv_dt_satn_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (dv_dt_wake_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      k_w = k_dv_dt_wake_start + k - 1
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      work_halo_v(i,j,k_w) = dv_dt_wake_land(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF

IF (tausy_d_on) THEN 
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    k_w = k_tausy_d_start 
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    work_halo_v(i,j,k_w) = tausy_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

!------------------------------------------------------------------
!l    4.4  start swap bounds of the super array (v fields) 
!------------------------------------------------------------------

! The variable work_halo_v is going to be used in the "p_to_v" routine
! call later which computes the interior points and requires the halo
! information to be set. This routine inside uses the field data at i
! index only (no i+1 or i-1). Thus the West and East halo information
! is not needed. Furthermore, the "p_to_v" routine uses the field data
! at j and j+1 indices with the starting point at zero (vdims%j_start)
! which is the South halo. Thus the North halo information is not
! needed. 
! Begin the swap bound of the South halo only.

CALL begin_swap_bounds(                                           &
     work_halo_v,                                                 &
     row_length, rows,                                            &
     k_total_v, off_x, off_y, fld_type_p, l_vector,               &
     do_east_arg=.FALSE.,  do_west_arg=.FALSE.,                   &
     do_south_arg=.TRUE., do_north_arg=.FALSE.,                   &
     do_corners_arg=.FALSE.)


!------------------------------------------------------------------
!l    4.5  end swap bounds of the super array (u fields) 
!------------------------------------------------------------------

CALL end_swap_bounds(work_halo_u)


!------------------------------------------------------------------
!l    4.6  interpole from  p to u staggering on 'c' grid 
!------------------------------------------------------------------


!------------------------------------------------------------------
!l    4.6.1  accelerations (u fields)
!------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                       &
!$OMP SHARED(pdims_s,udims,r_u,timestep,work_on_u_grid,           &
!$OMP    work_halo_u,levels,k_du_dt_start,k_du_dt_end,            &
!$OMP    k_stress_ud_start,k_stress_ud_end,stress_ud,             &
!$OMP    k_stress_ud_satn_start,k_stress_ud_satn_end,             &
!$OMP    stress_ud_satn,k_stress_ud_wake_start,                   &
!$OMP    k_stress_ud_wake_end,stress_ud_wake,k_du_dt_satn_start,  &
!$OMP    k_du_dt_satn_end,du_dt_satn,k_du_dt_wake_start,          &
!$OMP    k_du_dt_wake_end,du_dt_wake,stress_ud_on,stress_ud_p_on, &
!$OMP    stress_ud_satn_on,stress_ud_wake_on,du_dt_satn_on,       &
!$OMP    du_dt_satn_p_on,du_dt_wake_on)

CALL p_to_u(work_halo_u(:,:,k_du_dt_start:k_du_dt_end),     &
             pdims_s%i_start,pdims_s%i_end,                 &
             pdims_s%j_start,pdims_s%j_end,                 &
             udims%i_start,udims%i_end,                     &
             udims%j_start,udims%j_end,                     &
             1,levels,work_on_u_grid)

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
DO k=1,levels
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      r_u(i,j,k)=r_u(i,j,k)+ timestep*work_on_u_grid(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!------------------------------------------------------------------
!l    4.6.2 diagnostics (u fields)
!------------------------------------------------------------------


IF (stress_ud_on .OR. stress_ud_p_on) THEN  

  CALL p_to_u(work_halo_u(:,:,k_stress_ud_start:k_stress_ud_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               0,levels,stress_ud)

END IF 

IF (stress_ud_satn_on) THEN 

  CALL p_to_u(work_halo_u(:,:,k_stress_ud_satn_start:k_stress_ud_satn_end), &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               0,levels,stress_ud_satn)


END IF ! (stress_ud_satn_on) on stashflag


IF (stress_ud_wake_on) THEN 

  CALL p_to_u(work_halo_u(:,:,k_stress_ud_wake_start:k_stress_ud_wake_end), &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               0,levels,stress_ud_wake)

END IF ! (stress_ud_wake_on) on stashflag


IF (du_dt_satn_on .OR. du_dt_satn_p_on) THEN 

  CALL p_to_u(work_halo_u(:,:,k_du_dt_satn_start:k_du_dt_satn_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               1,levels,du_dt_satn)

END IF ! (du_dt_satn_on .OR. du_dt_satn_p_on) on stashflag

IF (du_dt_wake_on) THEN 

  CALL p_to_u(work_halo_u(:,:,k_du_dt_wake_start:k_du_dt_wake_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               1,levels,du_dt_wake)

END IF ! (du_dt_wake_on) on stashflag

!$OMP END PARALLEL

IF (tausx_d_on) THEN  

  CALL p_to_u(work_halo_u(:,:,k_tausx_d_start:k_tausx_d_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               udims%i_start,udims%i_end,                     &
               udims%j_start,udims%j_end,                     &
               1,1,tausx_d(udims%i_start,udims%j_start))

END IF ! (tausx_d_on) on stashflag


!------------------------------------------------------------------
!l    4.7  end swap bounds of the super array (v fields) 
!------------------------------------------------------------------

CALL end_swap_bounds(work_halo_v)


!------------------------------------------------------------------
!l    4.8  interpole from  p to v staggering on 'c' grid 
!------------------------------------------------------------------


!------------------------------------------------------------------
!l    4.8.1  accelerations (v fields)
!------------------------------------------------------------------


!$OMP PARALLEL DEFAULT(NONE)                                &
!$OMP SHARED(pdims_s,vdims,work_halo_v,k_dv_dt_start,       &
!$OMP    k_dv_dt_end,work_on_v_grid,levels)

CALL p_to_v(work_halo_v(:,:,k_dv_dt_start:k_dv_dt_end),     &
             pdims_s%i_start,pdims_s%i_end,                 &
             pdims_s%j_start,pdims_s%j_end,                 &
             vdims%i_start,vdims%i_end,                     &
             vdims%j_start,vdims%j_end,                     &
             1,levels,work_on_v_grid)

!$OMP END PARALLEL

! correct dv/dt at poles
IF (model_type == mt_global) THEN
  IF ( at_extremity(psouth) ) THEN
    
    CALL eg_v_at_poles(work_on_u_grid, work_on_v_grid, 1.0, &
                       udims%j_start, vdims%j_start,        &
                       udims,vdims)
    
  END IF

  IF ( at_extremity(pnorth) ) THEN

    CALL eg_v_at_poles(work_on_u_grid, work_on_v_grid, -1.0,&
                       udims%j_end, vdims%j_end,            &
                       udims,vdims)
  END IF
END IF

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(k,j,i)                       &
!$OMP SHARED(vdims,r_v,timestep,work_on_v_grid,levels)

!$OMP DO SCHEDULE(STATIC)
DO k=1,levels
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      r_v(i,j,k)=r_v(i,j,k)+ timestep*work_on_v_grid(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

!------------------------------------------------------------------
!l    4.8.2 diagnostics (v fields)
!------------------------------------------------------------------

IF (stress_vd_on) THEN  

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP SHARED(vdims,pdims_s,k_stress_vd_start,k_stress_vd_end,     &
!$OMP    work_halo_v,levels,stress_vd)

  CALL p_to_v(work_halo_v(:,:,k_stress_vd_start:k_stress_vd_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               0,levels,stress_vd)

!$OMP END PARALLEL

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(stress_ud, stress_vd, 1.0,           &
                         udims%j_start, vdims%j_start,        &
                         gwd_udims,gwd_vdims)

    END IF

    IF ( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(stress_ud, stress_vd, -1.0,          &
                         udims%j_end, vdims%j_end,            &
                         gwd_udims,gwd_vdims)

    END IF
  END IF

END IF ! (stress_vd_on) on stashflag


IF (stress_vd_satn_on) THEN 

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP SHARED(vdims,pdims_s,k_stress_vd_satn_start,                &
!$OMP    k_stress_vd_satn_end,work_halo_v,levels,stress_vd_satn)

  CALL p_to_v(work_halo_v(:,:,k_stress_vd_satn_start:k_stress_vd_satn_end), &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               0,levels,stress_vd_satn)

!$OMP END PARALLEL

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN
      
      CALL eg_v_at_poles(stress_ud_satn, stress_vd_satn, 1.0, &
                         udims%j_start, vdims%j_start,        &
                         gwd_udims,gwd_vdims)

    END IF

    IF ( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(stress_ud_satn, stress_vd_satn, -1.0,&
                         udims%j_end, vdims%j_end,            &
                         gwd_udims,gwd_vdims)
      
    END IF
  END IF

END IF ! (stress_vd_satn_on) on stashflag


IF (stress_vd_wake_on) THEN  

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP SHARED(vdims,pdims_s,k_stress_vd_wake_start,                &
!$OMP    k_stress_vd_wake_end,work_halo_v,levels,stress_vd_wake)

  CALL p_to_v(work_halo_v(:,:,k_stress_vd_wake_start:k_stress_vd_wake_end), &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               0,levels,stress_vd_wake)

!$OMP END PARALLEL

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(stress_ud_wake, stress_vd_wake, 1.0, &
                         udims%j_start, vdims%j_start,        &
                         gwd_udims,gwd_vdims)

    END IF

    IF ( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(stress_ud_wake, stress_vd_wake, -1.0,&
                         udims%j_end, vdims%j_end,            &
                         gwd_udims,gwd_vdims)

    END IF
  END IF

END IF ! (stress_vd_wake_on) on stashflag


IF (dv_dt_satn_on) THEN 

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP SHARED(vdims,pdims_s,k_dv_dt_satn_start,work_halo_v,        &
!$OMP    k_dv_dt_satn_end,levels,dv_dt_satn)

  CALL p_to_v(work_halo_v(:,:,k_dv_dt_satn_start:k_dv_dt_satn_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               1,levels,dv_dt_satn)

!$OMP END PARALLEL

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(du_dt_satn, dv_dt_satn, 1.0,         &
                         udims%j_start, vdims%j_start,        &
                         udims,vdims)
      
    END IF
    
    IF ( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(du_dt_satn, dv_dt_satn, -1.0,        &
                         udims%j_end, vdims%j_end,            &
                         udims,vdims)

    END IF
  END IF

END IF ! (dv_dt_satn_on) on stashflag


IF (dv_dt_wake_on) THEN  

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP SHARED(vdims,pdims_s,k_dv_dt_wake_start,work_halo_v,        &
!$OMP    k_dv_dt_wake_end,levels,dv_dt_wake)

  CALL p_to_v(work_halo_v(:,:,k_dv_dt_wake_start:k_dv_dt_wake_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               1,levels,dv_dt_wake)

!$OMP END PARALLEL

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(du_dt_wake, dv_dt_wake, 1.0,         &
                         udims%j_start, vdims%j_start,        &
                         udims,vdims)

    END IF

    IF ( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(du_dt_wake, dv_dt_wake, -1.0,        &
                         udims%j_end, vdims%j_end,            &
                         udims,vdims)

    END IF
  END IF

END IF ! (dv_dt_wake_on) on stashflag


IF (tausy_d_on) THEN  

  CALL p_to_v(work_halo_v(:,:,k_tausx_d_start:k_tausx_d_end),   &
               pdims_s%i_start,pdims_s%i_end,                 &
               pdims_s%j_start,pdims_s%j_end,                 &
               vdims%i_start,vdims%i_end,                     &
               vdims%j_start,vdims%j_end,                     &
               1,1,tausy_d(vdims%i_start,vdims%j_start))

  ! correct stress at poles
  IF (model_type == mt_global) THEN
    IF ( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(tausx_d,tausy_d, 1.0,                &
                         udims%j_start, vdims%j_start,        &
                         gwd_udims1,gwd_vdims1)

    END IF

    IF ( at_extremity(pnorth) ) THEN
     
      CALL eg_v_at_poles(tausx_d,tausy_d, -1.0,               &
                         udims%j_end, vdims%j_end,            &
                         gwd_udims1,gwd_vdims1)

    END IF
  END IF

END IF ! (tausy_d_on) on stashflag

!------------------------------------------------------------------
!l    5. scatter temperature increments to full area,
!l       and update
!------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( rows, row_length, u_s_d, land_points, land_index,        &
!$OMP         u_s_d_land, v_s_d, v_s_d_land, u_s_d_on, v_s_d_on,       &
!$OMP         nsq_s_d_on, nsq_s_d, bld_d_on, bld_d, bld_d_land,        &
!$OMP         bldt_d_on, bldt_d, bldt_d_land, nsq_s_d_land, fr_d_on,   &
!$OMP         fr_d, fr_d_land, levels, l_gw_heating, T_inc, timestep,  &
!$OMP         dt_dt)                                                   &
!$OMP PRIVATE( k, i, j, l )

IF ( l_gw_heating(1) .OR. l_gw_heating(2)) THEN !Apply heating tendency
!$OMP DO SCHEDULE(STATIC)
  DO k=1,levels
    DO l=1,land_points
      j=(land_index(l)-1)/row_length + 1
      i=land_index(l) - (j-1)*row_length
      T_inc(i,j,k)=T_inc(i,j,k)+ timestep*dt_dt(l,k)
    END DO ! l
  END DO ! k
!$OMP END DO NOWAIT
END IF ! l_gw_heating


IF (u_s_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      u_s_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    u_s_d(i,j) = u_s_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (u_s_d_on) on stashflag

IF (v_s_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      v_s_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    v_s_d(i,j) = v_s_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (v_s_d_on) on stashflag

IF (nsq_s_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      nsq_s_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    nsq_s_d(i,j) = nsq_s_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (nsq_s_d_on) on stashflag

IF (fr_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      fr_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    fr_d(i,j) = fr_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (fr_d_on) on stashflag

IF (bld_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      bld_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO 
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    bld_d(i,j) = bld_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (bld_d_on) on stashflag

IF (bldt_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
!$OMP DO SCHEDULE(STATIC)
  DO j=1,rows
    DO i=1,row_length
      bldt_d(i,j) = 0.0
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    bldt_d(i,j) = bldt_d_land(l)
  END DO ! l
!$OMP END DO NOWAIT
END IF ! (bldt_d_on) on stashflag

!$OMP END PARALLEL

IF (orog_slope_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
  DO j=1,rows
    DO i=1,row_length
      orog_slope_d(i,j) = 0.0
    END DO
  END DO
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    orog_slope_d(i,j) = slope(l)
  END DO ! l
END IF ! (orog_slope_d_on) on stashflag

IF (orog_anis_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
  DO j=1,rows
    DO i=1,row_length
      orog_anis_d(i,j) = 0.0
    END DO
  END DO
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    orog_anis_d(i,j) = anis(l)
  END DO ! l
END IF ! (orog_anis_d_on) on stashflag

IF (orog_dir_d_on) THEN  ! stashflag for this diagnostic
  ! expand from land points
  DO j=1,rows
    DO i=1,row_length
      orog_dir_d(i,j) = 0.0
    END DO
  END DO
  DO l=1,land_points
    j=(land_index(l)-1)/row_length + 1
    i=land_index(l) - (j-1)*row_length
    orog_dir_d(i,j) = mtdir(l)
  END DO ! l
END IF ! (orog_dir_d_on) on stashflag

DEALLOCATE(work_halo_v)
DEALLOCATE(work_halo_u)

iret=0
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE g_wave_5a

! ************************************************************************

SUBROUTINE set_segments(points_in, s_start_in,  s_end_in,  s_len_in, &
                                   s_start_out, s_end_out, s_len_out)

! This is a tiny subroutine, that is purely to improve code readability
! in the main subroutine by replacing lots of identical logic blocks with
! shorter subroutine calls
! It deliberately doesn't use DrHook - costs will be minimal, often
! inlined in the above routine and the DrHook overhead will be greater
! than the cost of the routine itself!

IMPLICIT NONE

INTEGER, INTENT(IN)  :: points_in   ! How many points before segmentation?
INTEGER, INTENT(IN)  :: s_start_in  ! input segment start 
INTEGER, INTENT(IN)  :: s_end_in    ! input segment end
INTEGER, INTENT(IN)  :: s_len_in    ! input segment length
INTEGER, INTENT(OUT) :: s_start_out ! output segment start 
INTEGER, INTENT(OUT) :: s_end_out   ! output segment end
INTEGER, INTENT(OUT) :: s_len_out   ! output segment length

IF (points_in /= 1) THEN
  s_start_out = s_start_in
  s_end_out   = s_end_in
  s_len_out   = s_len_in
ELSE
  s_start_out = 1
  s_end_out   = 1
  s_len_out   = 1
END IF

RETURN
END SUBROUTINE set_segments

! ************************************************************************

SUBROUTINE set_levels(flag, k, levels, k_start,  k_end)

! This is a tiny subroutine, that is purely to improve code readability
! in the main subroutine by replacing lots of identical logic blocks with
! shorter subroutine calls
! It deliberately doesn't use DrHook - costs will be minimal, often
! inlined in the above routine and the DrHook overhead will be greater
! than the cost of the routine itself!

IMPLICIT NONE

LOGICAL, INTENT(IN)     :: flag  
INTEGER, INTENT(INOUT)  :: k     
INTEGER, INTENT(IN)     :: levels
INTEGER, INTENT(OUT)    :: k_start     
INTEGER, INTENT(OUT)    :: k_end     

IF (flag) THEN
  k_start     = k
  k_end       = k_start + levels - 1
  k           = k_end + 1
ELSE
  k_start     = 1
  k_end       = 0 
END IF

RETURN
END SUBROUTINE set_levels

END MODULE
