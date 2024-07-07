! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! purpose: Interface to Atmospheric Physics GWD Schemes.
!         Scheme 1: Orographic flow blocking and gravity wave scheme.
!         Scheme 2: Non-orographic ultra-simple spectral gw scheme.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity wave drag


MODULE NI_gwd_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NI_GWD_CTL_MOD'

CONTAINS

SUBROUTINE NI_gwd_ctl (                                           &

! Parallel variables
        halo_i, halo_j, off_x, off_y                                    &
      ,  global_row_length,n_proc, n_procy, proc_row_group              &
      , at_extremity, neighbour                                         &

! model dimensions.
      , row_length, rows, n_rows, land_points                           &
      , model_levels                                                    &

! trig arrays
      , sin_theta_longitude, sin_theta_latitude                         &

! in coordinate information
      , delta_lambda,delta_phi,true_latitude                            &
      , exner_theta_levels                                              &

! in time stepping information.
      , timestep, timestep_number                                       &

! diagnostic info
      , STASHwork                                                       &

! SCM diagnostics (dummy in full UM)
      , nSCMDpkgs, L_SCMDiags                                           &
! in data fields.
      , u, v                                                            &
      , totalppn, land_sea_mask, p_layer_boundaries                     &
      , rho, theta_latest, sd_orog_land, orog_grad_xx_land              &
      , orog_grad_xy_land, orog_grad_yy_land, land_index                &

! in/out
      , R_u, R_v, T_inc                                                 &

! error information
      , Error_code  )

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                  &
   udims, vdims, wdims, tdims, pdims                              &
,  udims_s, vdims_s, wdims_s, tdims_s, pdims_s


! Model level heights from centre of Earth
USE level_heights_mod, ONLY:  &
   r_theta_levels             &  ! Radii on theta levels (m)
,  r_rho_levels                  ! Radii on rho levels (m)

USE g_wave_input_mod, ONLY: l_gwd,                &
                             L_use_ussp,           &
                             l_gwd_40km,           &
                             l_nonhydro,           &
                             l_dynbeta,            &
                             l_taus_scale,         &
                             l_fix_gwsatn,         &
                             sat_scheme,           &
                             i_gwd_vn,             &
                             i_gwd_vn_4a,          &
                             i_gwd_vn_5a,          &
                             kay_gwave,            &
                             gwd_frc,              &
                             nsigma,               &
                             gwd_fsat,             &
                             gsharp,               &
                             fbcd,                 &
                             l_smooth,             &
                             l_gw_heating,         &
                             i_ussp_vn,            &
                             i_ussp_vn_1a

! ----------------------------------------------------------------------+-------

USE tuning_segments_mod, ONLY : gw_seg_size

USE g_wave_4a_mod, ONLY: g_wave_4a
USE g_wave_5a_mod, ONLY: g_wave_5a
USE gw_ussp_mod, ONLY: gw_ussp

USE ereport_mod, ONLY: Ereport
USE umPrintMgr

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stash_array_mod, ONLY: sf, si
USE um_stashcode_mod, ONLY: stashcode_total_precip
USE physics_tendencies_mod,  ONLY:                                &
    l_retain_gwd_tendencies, dt_gwd, du_gwd, dv_gwd

USE s_scmop_mod,          ONLY: t_inst, d_sl, d_all,              &
                                default_streams, scmdiag_gwd
USE scmoutput_mod,        ONLY: scmoutput
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_single_column

USE diagnostics_gwd_mod, ONLY: diagnostics_gwd

USE pws_diags_mod, ONLY: flag_wafc_cat, flag_mtn_wave_turb

IMPLICIT NONE

!     Fixed starting theta-level (e.g. for P_layer_boundaries)
INTEGER, PARAMETER :: tkfix0start    = 0
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
INTEGER, PARAMETER :: tkfix1start    = 1
!
! arguments with intent in. ie: input variables.

! Parallel setup variables
INTEGER ::                                                        &
  halo_i                                                          &
             ! Size of halo in i direction.
, halo_j                                                          &
             ! Size of halo in j direction.
, off_x                                                           &
             ! Size of small halo in i
, off_y                                                           &
             ! Size of small halo in j.
, global_row_length                                               &
             ! number of points on a row
, proc_row_group                                                  &
             ! Group id for processors on the same row
, n_proc                                                          &
             ! Total number of processors
, n_procy                                                         &
             ! Number of processors in latitude
, neighbour(4)   ! Array with the Ids of the four neighbours
                 ! in the horizontal plane

LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north, south
                   ! east or west of the processor grid


! Model dimensions
INTEGER ::                                                        &
  row_length                                                      &
, rows                                                            &
, n_rows                                                          &
, model_levels                                                    &
, land_points

! model parameters
REAL ::                                                           &
  timestep


! Diagnostics info
REAL ::                                                          &
STASHwork(*) ! STASH workspace for section 6 (GW Drag)

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN)  ::  &
  nSCMDpkgs                ! No of SCM diagnostics packages

LOGICAL, INTENT(IN)  ::  &
  L_SCMDiags(nSCMDpkgs)    ! Logicals for SCM diagnostics packages

! Data arrays

INTEGER ::                                                        &
  land_index (land_points)            ! set from land_sea_mask

REAL ::                                                           &
 u(udims_s%i_start:udims_s%i_end,                                 &
   udims_s%j_start:udims_s%j_end,                                 &
   udims_s%k_start:udims_s%k_end)                                 &
      !primary u field (ms**-1)
,v(vdims_s%i_start:vdims_s%i_end,                                 &
   vdims_s%j_start:vdims_s%j_end,                                 &
   vdims_s%k_start:vdims_s%k_end)                                 &
      !primary v field (ms**-1)
,rho(pdims_s%i_start:pdims_s%i_end,                               &
     pdims_s%j_start:pdims_s%j_end,                               &
     pdims_s%k_start:pdims_s%k_end)
      !density *r*r (kg/m)
REAL ::                                                           &
  theta_latest(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end,tkfix1start:tdims%k_end)
REAL ::                                                           &
  totalppn(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

REAL ::                                                           &
  p_layer_boundaries(tdims%i_start:tdims%i_end,                   &
                     tdims%j_start:tdims%j_end,                   &
                     tkfix0start:tdims%k_end)                     &
        ! pressure at layer boundaries. Same as p except at
        ! bottom level = pstar, and at top = 0.
, sd_orog_land (land_points)                                      &
                             ! orog/qrparm.orog.stdev
, orog_grad_xx_land(land_points)                                  &
                                 ! orog/qrparm.orog.sigmaxx
, orog_grad_xy_land(land_points)                                  &
                                 ! orog/qrparm.orog.sigmaxy
, orog_grad_yy_land(land_points) ! orog/qrparm.orog.sigmayy


REAL ::                                                           &
 exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                &
                    tdims_s%j_start:tdims_s%j_end,                &
                    tdims_s%k_start:tdims_s%k_end)
                                 ! Exner on theta level

LOGICAL ::                                                        &
  land_sea_mask(row_length, rows)

! Co-ordinate arrays
REAL ::                                                           &
! local vertical co-ordinate information
   delta_lambda                                                    &
 , delta_phi                                                       &
 , true_latitude(row_length, rows)

REAL ::                                                           &
  sin_theta_longitude (row_length, rows)                          &
, sin_theta_latitude  (row_length, rows)

! time information for current timestep
INTEGER ::                                                        &
  timestep_number


! arguments with intent in/out. ie: input variables changed on output.

! arguments with intent out. ie: output variables.
REAL       ::                                                     &
  r_u(udims_s%i_start:udims_s%i_end,                              &
      udims_s%j_start:udims_s%j_end,                              &
      udims_s%k_start:udims_s%k_end)                              &
                                     !u wind increment diagnostic
, r_v(vdims_s%i_start:vdims_s%i_end,                              &
      vdims_s%j_start:vdims_s%j_end,                              &
      vdims_s%k_start:vdims_s%k_end)                              &
                                     !v wind increment diagnostic
, T_inc(tdims%i_start:tdims%i_end,                                &
        tdims%j_start:tdims%j_end,                                &
        1:tdims%k_end)               !Temperature increment

INTEGER ::                                                        &
  Error_code

! local variables
LOGICAL ::                                                        &
  stress_ud_on                                                    &
, stress_ud_p_on                                                  &
, stress_vd_on                                                    &
, stress_ud_satn_on                                               &
, stress_vd_satn_on                                               &
, stress_ud_wake_on                                               &
, stress_vd_wake_on                                               &
, du_dt_satn_on                                                   &
, du_dt_satn_p_on                                                 &
, dv_dt_satn_on                                                   &
, du_dt_wake_on                                                   &
, dv_dt_wake_on                                                   &
, gwspec_eflux_on                                                 &
, gwspec_eflux_p_on                                               &
, gwspec_sflux_on                                                 &
, gwspec_wflux_on                                                 &
, gwspec_wflux_p_on                                               &
, gwspec_nflux_on                                                 &
, gwspec_ewacc_on                                                 &
, gwspec_ewacc_p_on                                               &
, gwspec_nsacc_on                                                 &
, u_s_d_on                                                        &
, v_s_d_on                                                        &
, nsq_s_d_on                                                      &
, fr_d_on                                                         &
, bld_d_on                                                        &
, bldt_d_on                                                       &
, num_lim_d_on                                                    &
, num_fac_d_on                                                    &
, tausx_d_on                                                      &
, tausy_d_on                                                      &
, taus_scale_d_on                                                 &
, orog_slope_d_on                                                 &
, orog_anis_d_on                                                  &
, orog_dir_d_on                                                   &
, L_u_incr_gwd                                                    &
, L_v_incr_gwd                                                    &
, L_t_incr_gwd                                                    &
, l_make_pws_prereq

INTEGER ::                                                        &
  i,j,k      ! loop counters

INTEGER ::                                                        &
  errorstatus      ! Return code : 0 Normal Exit : >0 Error

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage         ! Error message if return code >0


! Local data arrays
!

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
REAL,  ALLOCATABLE::                               &
  u_incr_diagnostic(:,:,:)                                               &
                        ! u wind increment for STASH
, v_incr_diagnostic(:,:,:)                                               &
                        ! v wind increment for STASH
, t_incr_diagnostic(:,:,:)                                               &
                        ! temperature increment for STASH
, stress_ud(:,:,:)                                                       &
, stress_vd(:,:,:)                                                       &
, stress_ud_satn (:,:,:)                                                 &
, stress_vd_satn(:,:,:)                                                  &
, stress_ud_wake(:,:,:)                                                  &
, stress_vd_wake(:,:,:)                                                  &
, du_dt_satn (:,:,:)                                                     &
, dv_dt_satn (:,:,:)                                                     &
, du_dt_wake (:,:,:)                                                     &
, dv_dt_wake (:,:,:)                                                     &
, gwspec_eflux(:,:,:)                                                    &
, gwspec_sflux(:,:,:)                                                    &
, gwspec_wflux(:,:,:)                                                    &
, gwspec_nflux(:,:,:)                                                    &
, gwspec_ewacc(:,:,:)                                                    &
, gwspec_nsacc(:,:,:)

REAL,    ALLOCATABLE::                                 &
  u_s_d (:,:)                                                          &
, v_s_d (:,:)                                                          &
, nsq_s_d (:,:)                                                        &
, fr_d  (:,:)                                                          &
, bld_d (:,:)                                                          &
, bldt_d (:,:)                                                         &
, num_lim_d (:,:)                                                      &
, num_fac_d (:,:)                                                      &
, tausx_d (:,:)                                                        &
, tausy_d (:,:)                                                        &
, taus_scale_d (:,:)                                                   &
, orog_slope_d (:,:)                                                   &
, orog_anis_d  (:,:)                                                   &
, orog_dir_d   (:,:)

! Diagnostic land_point array sizes
INTEGER ::                                                        &
 points_stress_ud                                                 &
,points_stress_vd                                                 &
,points_stress_ud_satn                                            &
,points_stress_vd_satn                                            &
,points_stress_ud_wake                                            &
,points_stress_vd_wake                                            &
,points_du_dt_satn                                                &
,points_dv_dt_satn                                                &
,points_du_dt_wake                                                &
,points_dv_dt_wake                                                &
,points_u_s_d                                                     &
,points_v_s_d                                                     &
,points_nsq_s_d                                                   &
,points_fr_d                                                      &
,points_bld_d                                                     &
,points_bldt_d                                                    &
,points_num_lim_d                                                 &
,points_num_fac_d                                                 &
,points_tausx_d                                                   &
,points_tausy_d                                                   &
,points_taus_scale_d                                              &
,points_orog_slope_d                                              &
,points_orog_anis_d                                               &
,points_orog_dir_d

INTEGER ::            &
  u_rows,             & ! rows on u grid
  v_rows                ! rows on v grid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
! Error reporting variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'NI_GWD_CTL'

! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches
!-----------------------------------------------------------------------------
! Section 0 - initialisation and setup
!-----------------------------------------------------------------------------
! Work out number of rows for U and V grids - will depend on whether
! ENDGame or not.

u_rows = udims%j_len
v_rows = vdims%j_len
! ----------------------------------------------------------------------
! Section GWD.1 Set stash diagnostic switches
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! General case of the atmosphere model, ie : with stash.
SELECT CASE (model_type)
CASE DEFAULT

  gwspec_eflux_on = sf(101,6)
  gwspec_eflux_p_on = sf(111,6)
  gwspec_sflux_on = sf(102,6)
  gwspec_wflux_on = sf(103,6)
  gwspec_wflux_p_on = sf(113,6)
  gwspec_nflux_on = sf(104,6)
  gwspec_ewacc_on = sf(105,6)
  gwspec_ewacc_p_on = sf(115,6)
  gwspec_nsacc_on = sf(106,6)
  L_t_incr_gwd      = sf(181,6)
  L_u_incr_gwd      = sf(185,6)
  L_v_incr_gwd      = sf(186,6)
  ! PWS diagnostics 20/007 and 20/047 require the fields in items
  ! 6/201 and 6/202.
  l_make_pws_prereq = (flag_wafc_cat .OR. flag_mtn_wave_turb)
  stress_ud_on      = (sf(201,6) .OR. l_make_pws_prereq)
  stress_ud_p_on    = sf(241,6)
  stress_vd_on      = (sf(202,6) .OR. l_make_pws_prereq)
  du_dt_satn_on     = sf(207,6)
  du_dt_satn_p_on   = sf(247,6)
  dv_dt_satn_on     = sf(208,6)
  u_s_d_on          = sf(214,6)
  v_s_d_on          = sf(215,6)
  nsq_s_d_on        = sf(216,6)
  fr_d_on           = sf(217,6)
  bld_d_on          = sf(218,6)
  bldt_d_on         = sf(222,6)
  stress_ud_satn_on = sf(223,6)
  stress_vd_satn_on = sf(224,6)
  stress_ud_wake_on = sf(227,6)
  stress_vd_wake_on = sf(228,6)
  du_dt_wake_on     = sf(231,6)
  dv_dt_wake_on     = sf(232,6)
  num_lim_d_on      = sf(233,6)
  num_fac_d_on      = sf(234,6)
  tausx_d_on        = sf(235,6)
  tausy_d_on        = sf(236,6)
  taus_scale_d_on   = sf(237,6)
  orog_slope_d_on   = sf(248,6)
  orog_anis_d_on    = sf(249,6)
  orog_dir_d_on     = sf(250,6)
CASE (mt_single_column)
      ! These fields will need to be allocated if the USSP scheme is used
  gwspec_eflux_on   = l_use_ussp
  gwspec_eflux_p_on = l_use_ussp
  gwspec_sflux_on   = l_use_ussp
  gwspec_wflux_on   = l_use_ussp
  gwspec_wflux_p_on = l_use_ussp
  gwspec_nflux_on   = l_use_ussp
  gwspec_ewacc_on   = l_use_ussp
  gwspec_ewacc_p_on = l_use_ussp
  gwspec_nsacc_on   = l_use_ussp
  L_u_incr_gwd      = .FALSE.
  L_v_incr_gwd      = .FALSE.
  L_t_incr_gwd      = .FALSE.
  stress_ud_p_on    = .FALSE.
  du_dt_satn_p_on   = .FALSE.
  bldt_d_on         = .FALSE.
  num_lim_d_on      = .FALSE.
  num_fac_d_on      = .FALSE.
  tausx_d_on        = .FALSE.
  tausy_d_on        = .FALSE.
  taus_scale_d_on   = .FALSE.
  orog_slope_d_on   = .FALSE.
  orog_anis_d_on    = .FALSE.
  orog_dir_d_on     = .FALSE.
  IF (L_SCMDiags(SCMDiag_gwd)) THEN
    stress_ud_on      = .TRUE.
    stress_vd_on      = .TRUE.
    stress_ud_satn_on = .TRUE.
    stress_vd_satn_on = .TRUE.
    stress_ud_wake_on = .TRUE.
    stress_vd_wake_on = .TRUE.
    du_dt_satn_on     = .TRUE.
    dv_dt_satn_on     = .TRUE.
    du_dt_wake_on     = .TRUE.
    dv_dt_wake_on     = .TRUE.
    u_s_d_on          = .TRUE.
    v_s_d_on          = .TRUE.
    nsq_s_d_on        = .TRUE.
    fr_d_on           = .TRUE.
    bld_d_on          = .TRUE.
  ELSE
    stress_ud_on      = .FALSE.
    stress_vd_on      = .FALSE.
    stress_ud_satn_on = .FALSE.
    stress_vd_satn_on = .FALSE.
    stress_ud_wake_on = .FALSE.
    stress_vd_wake_on = .FALSE.
    du_dt_satn_on     = .FALSE.
    dv_dt_satn_on     = .FALSE.
    du_dt_wake_on     = .FALSE.
    dv_dt_wake_on     = .FALSE.
    u_s_d_on          = .FALSE.
    v_s_d_on          = .FALSE.
    nsq_s_d_on        = .FALSE.
    fr_d_on           = .FALSE.
    bld_d_on          = .FALSE.
  END IF
END SELECT ! model_type


! ----------------------------------------------------------------------
! Section GWD.1
! ----------------------------------------------------------------------

IF (error_code  ==  0 ) THEN
  ! Save R_u    before updating
  IF ( l_u_incr_gwd) THEN  ! STASHflag set
    ALLOCATE ( u_incr_diagnostic(udims%i_start:udims%i_end,      &
                                 udims%j_start:udims%j_end,      &
                                 udims%k_start:udims%k_end)  )
    DO k=1,model_levels
      DO j=udims%j_start,udims%j_end
        DO i=udims%i_start,udims%i_end
          u_incr_diagnostic(i,j,k) = R_u(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( u_incr_diagnostic(1,1,1) )

  END IF                    ! on STASHflag

  ! Save R_v    before updating
  IF ( l_v_incr_gwd) THEN  ! STASHflag set
    ALLOCATE ( v_incr_diagnostic(vdims%i_start:vdims%i_end,      &
                                 vdims%j_start:vdims%j_end,      &
                                 vdims%k_start:vdims%k_end)  )
    DO k=1,model_levels
      DO j=vdims%j_start,vdims%j_end
        DO i=vdims%i_start,vdims%i_end
          v_incr_diagnostic(i,j,k) = R_v(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( v_incr_diagnostic(1,1,1) )
  END IF                    ! on STASHflag

  ! Save T_inc   before updating
  IF ( l_t_incr_gwd) THEN  ! STASHflag set
    ALLOCATE ( t_incr_diagnostic(tdims%i_start:tdims%i_end,      &
                                 tdims%j_start:tdims%j_end,      &
                                 tkfix1start:tdims%k_end)  )
    DO k=tkfix1start,model_levels
      DO j=tdims%j_start,tdims%j_end
        DO i=tdims%i_start,tdims%i_end
          t_incr_diagnostic(i,j,k) = T_inc(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( t_incr_diagnostic(1,1,1) )
  END IF                    ! on STASHflag

  IF ( stress_ud_on .OR. stress_ud_p_on) THEN  ! STASHflag set
    ALLOCATE ( stress_ud(row_length,u_rows,0:model_levels) )
    points_stress_ud = land_points
  ELSE
    points_stress_ud = 1
    ALLOCATE ( stress_ud(1,1,1) )
  END IF

  IF ( stress_vd_on ) THEN  ! STASHflag set
    ALLOCATE ( stress_vd(row_length,v_rows,0:model_levels) )
    points_stress_vd = land_points
  ELSE
    points_stress_vd = 1
    ALLOCATE ( stress_vd(1,1,1) )
  END IF

  IF ( du_dt_satn_on .OR. du_dt_satn_p_on ) THEN  ! STASHflag set
    ALLOCATE ( du_dt_satn(row_length,u_rows,model_levels) )
    points_du_dt_satn = land_points
  ELSE
    points_du_dt_satn = 1
    ALLOCATE ( du_dt_satn(1,1,1) )
  END IF

  IF ( dv_dt_satn_on ) THEN  ! STASHflag set
    ALLOCATE ( dv_dt_satn(row_length,v_rows,model_levels) )
    points_dv_dt_satn = land_points
  ELSE
    points_dv_dt_satn = 1
    ALLOCATE ( dv_dt_satn(1,1,1) )
  END IF

  IF ( u_s_d_on ) THEN  ! STASHflag set
    ALLOCATE ( u_s_d(row_length,rows) )
    points_u_s_d = land_points
  ELSE
    points_u_s_d = 1
    ALLOCATE ( u_s_d(1,1) )
  END IF

  IF ( v_s_d_on ) THEN  ! STASHflag set
    ALLOCATE ( v_s_d(row_length,rows) )
    points_v_s_d = land_points
  ELSE
    points_v_s_d = 1
    ALLOCATE ( v_s_d(1,1) )
  END IF

  IF ( nsq_s_d_on ) THEN  ! STASHflag set
    ALLOCATE ( nsq_s_d(row_length,rows) )
    points_nsq_s_d = land_points
  ELSE
    points_nsq_s_d = 1
    ALLOCATE ( nsq_s_d(1,1) )
  END IF

  IF ( fr_d_on ) THEN  ! STASHflag set
    ALLOCATE ( fr_d(row_length,rows) )
    points_fr_d = land_points
  ELSE
    points_fr_d = 1
    ALLOCATE ( fr_d(1,1) )
  END IF

  IF ( bld_d_on ) THEN  ! STASHflag set
    ALLOCATE ( bld_d(row_length,rows) )
    points_bld_d = land_points
  ELSE
    points_bld_d = 1
    ALLOCATE ( bld_d(1,1) )
  END IF

  IF ( bldt_d_on ) THEN  ! STASHflag set
    ALLOCATE ( bldt_d(row_length,rows) )
    points_bldt_d = land_points
  ELSE
    points_bldt_d = 1
    ALLOCATE ( bldt_d(1,1) )
  END IF

  IF ( num_lim_d_on ) THEN  ! STASHflag set
    ALLOCATE ( num_lim_d(row_length,rows) )
    points_num_lim_d = land_points
  ELSE
    points_num_lim_d = 1
    ALLOCATE ( num_lim_d(1,1) )
  END IF

  IF ( num_fac_d_on ) THEN  ! STASHflag set
    ALLOCATE ( num_fac_d(row_length,rows) )
    points_num_fac_d = land_points
  ELSE
    points_num_fac_d = 1
    ALLOCATE ( num_fac_d(1,1) )
  END IF

  IF ( stress_ud_satn_on ) THEN  ! STASHflag set
    ALLOCATE ( stress_ud_satn(row_length,u_rows,0:model_levels) )
    points_stress_ud_satn = land_points
  ELSE
    points_stress_ud_satn = 1
    ALLOCATE ( stress_ud_satn(1,1,1) )
  END IF

  IF ( stress_vd_satn_on ) THEN  ! STASHflag set
    ALLOCATE ( stress_vd_satn(row_length,v_rows,0:model_levels) )
    points_stress_vd_satn = land_points
  ELSE
    points_stress_vd_satn = 1
    ALLOCATE ( stress_vd_satn(1,1,1) )
  END IF

  IF ( stress_ud_wake_on ) THEN  ! STASHflag set
    ALLOCATE ( stress_ud_wake(row_length,u_rows,0:model_levels) )
    points_stress_ud_wake = land_points
  ELSE
    points_stress_ud_wake = 1
    ALLOCATE ( stress_ud_wake(1,1,1) )
  END IF

  IF ( stress_vd_wake_on ) THEN  ! STASHflag set
    ALLOCATE ( stress_vd_wake(row_length,v_rows,0:model_levels) )
    points_stress_vd_wake = land_points
  ELSE
    points_stress_vd_wake = 1
    ALLOCATE ( stress_vd_wake(1,1,1) )
  END IF

  IF ( du_dt_wake_on ) THEN  ! STASHflag set
    ALLOCATE ( du_dt_wake(row_length,u_rows,model_levels) )
    points_du_dt_wake = land_points
  ELSE
    points_du_dt_wake = 1
    ALLOCATE ( du_dt_wake(1,1,1) )
  END IF

  IF ( dv_dt_wake_on ) THEN  ! STASHflag set
    ALLOCATE ( dv_dt_wake(row_length,v_rows,model_levels) )
    points_dv_dt_wake = land_points
  ELSE
    points_dv_dt_wake = 1
    ALLOCATE ( dv_dt_wake(1,1,1) )
  END IF

  IF ( tausx_d_on ) THEN  ! STASHflag set
    ALLOCATE ( tausx_d(row_length,u_rows) )
    points_tausx_d = land_points
  ELSE
    points_tausx_d = 1
    ALLOCATE ( tausx_d(1,1) )
  END IF

  IF ( tausy_d_on ) THEN  ! STASHflag set
    ALLOCATE ( tausy_d(row_length,v_rows) )
    points_tausy_d = land_points
  ELSE
    points_tausy_d = 1
    ALLOCATE ( tausy_d(1,1) )
  END IF

  IF ( taus_scale_d_on ) THEN  ! STASHflag set
    ALLOCATE ( taus_scale_d(row_length,rows) )
    points_taus_scale_d = land_points
  ELSE
    points_taus_scale_d = 1
    ALLOCATE ( taus_scale_d(1,1) )
  END IF

  IF ( orog_slope_d_on ) THEN  ! STASHflag set
    ALLOCATE ( orog_slope_d(row_length,rows) )
    points_orog_slope_d = land_points
  ELSE
    points_orog_slope_d = 1
    ALLOCATE ( orog_slope_d(1,1) )
  END IF

  IF ( orog_anis_d_on ) THEN  ! STASHflag set
    ALLOCATE ( orog_anis_d(row_length,rows) )
    points_orog_anis_d = land_points
  ELSE
    points_orog_anis_d = 1
    ALLOCATE ( orog_anis_d(1,1) )
  END IF

  IF ( orog_dir_d_on ) THEN  ! STASHflag set
    ALLOCATE ( orog_dir_d(row_length,rows) )
    points_orog_dir_d = land_points
  ELSE
    points_orog_dir_d = 1
    ALLOCATE ( orog_dir_d(1,1) )
  END IF

  IF ( gwspec_eflux_on .OR. gwspec_eflux_p_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_eflux(udims%i_start:udims%i_end,               &
                udims%j_start:udims%j_end,tkfix1start:tdims%k_end))
  ELSE
    ALLOCATE (gwspec_eflux(1,1,1))
  END IF
  IF ( gwspec_sflux_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_sflux(vdims%i_start:vdims%i_end,               &
                vdims%j_start:vdims%j_end,tkfix1start:tdims%k_end))
  ELSE
    ALLOCATE (gwspec_sflux(1,1,1))
  END IF
  IF ( gwspec_wflux_on .OR. gwspec_wflux_p_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_wflux(udims%i_start:udims%i_end,               &
                udims%j_start:udims%j_end,tkfix1start:tdims%k_end))
  ELSE
    ALLOCATE (gwspec_wflux(1,1,1))
  END IF
  IF ( gwspec_nflux_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_nflux(vdims%i_start:vdims%i_end,               &
                vdims%j_start:vdims%j_end,tkfix1start:tdims%k_end))
  ELSE
    ALLOCATE (gwspec_nflux(1,1,1))
  END IF
  IF ( gwspec_ewacc_on .OR. gwspec_ewacc_p_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_ewacc(udims%i_start:udims%i_end,               &
               udims%j_start:udims%j_end,udims%k_start:udims%k_end))
  ELSE
    ALLOCATE (gwspec_ewacc(1,1,1))
  END IF
  IF ( gwspec_nsacc_on ) THEN  ! STASHflag set
    ALLOCATE (gwspec_nsacc(vdims%i_start:vdims%i_end,               &
               vdims%j_start:vdims%j_end,vdims%k_start:vdims%k_end))
  ELSE
    ALLOCATE (gwspec_nsacc(1,1,1))
  END IF
  
  ! Save initial values before GWD scheme
  IF (l_retain_gwd_tendencies) THEN

    ! Theta
    DO k = 1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          dt_gwd(i,j,k) = T_inc(i,j,k)
        END DO
      END DO
    END DO
    ! u
    DO k = udims%k_start,udims%k_end
      DO j = udims%j_start,udims%j_end
        DO i = udims%i_start,udims%i_end
          du_gwd(i,j,k) = R_u(i,j,k)
        END DO
      END DO
    END DO
    ! v
    DO k = vdims%k_start,vdims%k_end
      DO j = vdims%j_start,vdims%j_end
        DO i = vdims%i_start,vdims%i_end
          dv_gwd(i,j,k) = R_v(i,j,k)
        END DO
      END DO
    END DO

  END IF ! end if over l_retain_gwd_tendencies
  
  !-------------------------------------------------------------------
  !Section GWD.2b  Call Orographic Gravity Wave Drag Scheme
  !-------------------------------------------------------------------
  IF (L_gwd) THEN

    SELECT CASE ( i_gwd_vn )
    CASE ( i_gwd_vn_4a )

      CALL g_wave_4a(                                                &
      theta_latest, u, v, row_length, rows, n_rows,u_rows,v_rows,    &
      off_x, off_y,                                                  &
      global_row_length,n_proc, n_procy, proc_row_group,at_extremity,&
      model_levels, rho, delta_lambda,delta_phi,                     &
      true_latitude,sd_orog_land,                                    &
      orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
      land_index, land_points, timestep, kay_gwave, gwd_frc, nsigma, &
      r_u, r_v, T_inc, l_taus_scale, l_fix_gwsatn, l_gwd_40km,       &
      sat_scheme, gwd_fsat,                                          &
      Gsharp, fbcd, l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,   &
      stress_ud      , stress_ud_on     , stress_ud_p_on       ,     &
                                          points_stress_ud     ,     &
      stress_vd      , stress_vd_on     , points_stress_vd     ,     &
      stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
      stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
      stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
      stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
      du_dt_satn     , du_dt_satn_on    , du_dt_satn_p_on      ,     &
                                          points_du_dt_satn    ,     &
      dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
      du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
      dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
      u_s_d          , u_s_d_on         , points_u_s_d         ,     &
      v_s_d          , v_s_d_on         , points_v_s_d         ,     &
      nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
      fr_d           , fr_d_on          , points_fr_d          ,     &
      bld_d          , bld_d_on         , points_bld_d         ,     &
      bldt_d         , bldt_d_on        , points_bldt_d        ,     &
      num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
      num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
      tausx_d        , tausx_d_on       , points_tausx_d       ,     &
      tausy_d        , tausy_d_on       , points_tausy_d       ,     &
      taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
      orog_slope_d   , orog_slope_d_on  , points_orog_slope_d  ,     &
      orog_anis_d    , orog_anis_d_on   , points_orog_anis_d   ,     &
      orog_dir_d     , orog_dir_d_on    , points_orog_dir_d    ,     &
      error_code)

    CASE ( i_gwd_vn_5a )

      CALL g_wave_5a(                                                &
      theta_latest, u, v, row_length, rows, n_rows,u_rows,v_rows,    &
      off_x, off_y,                                                  &
      global_row_length,n_proc, n_procy, proc_row_group,at_extremity,&
      model_levels, rho, delta_lambda,delta_phi,                     &
      true_latitude,sd_orog_land,                                    &
      orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
      land_index, land_points, timestep, kay_gwave, gwd_frc, nsigma, &
      r_u, r_v, T_inc, l_taus_scale, l_fix_gwsatn, l_gwd_40km,       &
      sat_scheme, gwd_fsat, gw_seg_size,                             &
      Gsharp, fbcd, l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,   &
      stress_ud      , stress_ud_on     , stress_ud_p_on       ,     &
                                          points_stress_ud     ,     &
      stress_vd      , stress_vd_on     , points_stress_vd     ,     &
      stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
      stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
      stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
      stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
      du_dt_satn     , du_dt_satn_on    , du_dt_satn_p_on      ,     &
                                          points_du_dt_satn    ,     &
      dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
      du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
      dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
      u_s_d          , u_s_d_on         , points_u_s_d         ,     &
      v_s_d          , v_s_d_on         , points_v_s_d         ,     &
      nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
      fr_d           , fr_d_on          , points_fr_d          ,     &
      bld_d          , bld_d_on         , points_bld_d         ,     &
      bldt_d         , bldt_d_on        , points_bldt_d        ,     &
      num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
      num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
      tausx_d        , tausx_d_on       , points_tausx_d       ,     &
      tausy_d        , tausy_d_on       , points_tausy_d       ,     &
      taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
      orog_slope_d   , orog_slope_d_on  , points_orog_slope_d  ,     &
      orog_anis_d    , orog_anis_d_on   , points_orog_anis_d   ,     &
      orog_dir_d     , orog_dir_d_on    , points_orog_dir_d    ,     &
      error_code)

    CASE DEFAULT ! i_gwd_vn

      errorstatus = 10
      WRITE (cmessage,'(A,A,I6)') 'Gravity wave drag version ',     &
      'not recognised. i_gwd_vn = ',i_gwd_vn
      CALL Ereport ( RoutineName, errorstatus, cmessage)

    END SELECT ! i_gwd_vn

  END IF ! l_gwd

  IF ( L_use_ussp ) THEN
!   Check code to fail in absence of totalppn in D1 was tested here, where
!   failure would suggest option code logic for request in tstmask has been 
!   screwed up. Formally wish to know if D1 has Total ppn at Start of TS and 
!   fail if (.NOT. l_totppn_softs) but code as introduced would simply set
!   l_totppn_softs == l_use_ussp and thus makes such a test redundant.
!
!   Presently only one version of ussp scheme, hence CASE not required
!   SELECT CASE ( i_ussp_vn )
!
!   CASE ( i_ussp_vn_1a )
!!
      CALL gw_ussp(model_levels, rows, n_rows,                           &
           off_x, off_y, halo_i, halo_j, row_length,                     &
           global_row_length,n_proc,n_procy,proc_row_group,at_extremity, &
           r_rho_levels, r_theta_levels, p_layer_boundaries,             &
           sin_theta_longitude, sin_theta_latitude,                      &
           theta_latest, rho, u, v, totalppn, timestep,                  &
           r_u, r_v, T_inc,                                              &
           l_gw_heating(3),                                              &
           gwspec_eflux,gwspec_sflux,gwspec_wflux,gwspec_nflux,          &
           gwspec_ewacc,gwspec_nsacc,                                    &
           gwspec_eflux_on,gwspec_eflux_p_on,gwspec_sflux_on,            &
           gwspec_wflux_on,gwspec_wflux_p_on,gwspec_nflux_on,            &
           gwspec_ewacc_on,gwspec_ewacc_p_on,gwspec_nsacc_on)

!!
!   CASE DEFAULT ! i_ussp_vn
!
!     errorstatus = 10
!   Still run the Case Default check as flag for namelist propagation errors
    IF ( i_ussp_vn /= i_ussp_vn_1a )  THEN
      WRITE (cmessage,'(A,A,I6)') 'Non-orographic GW version ',          &
      'not recognised. i_ussp_vn = ',i_ussp_vn
      CALL Ereport ( RoutineName, errorstatus, cmessage)
    END IF
!
!   END SELECT ! i_ussp_vn
  END IF ! l_use_ussp
  
  ! Copy GWD tendencies to physics_tendencies_mod module
  IF (l_retain_gwd_tendencies) THEN

    ! Theta
    DO k = 1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          dt_gwd(i,j,k) = T_inc(i,j,k) - dt_gwd(i,j,k)
        END DO
      END DO
    END DO
    ! u
    DO k = udims%k_start,udims%k_end
      DO j = udims%j_start,udims%j_end
        DO i = udims%i_start,udims%i_end
          du_gwd(i,j,k) = R_u(i,j,k) - du_gwd(i,j,k)
        END DO
      END DO
    END DO
    ! v
    DO k = vdims%k_start,vdims%k_end
      DO j = vdims%j_start,vdims%j_end
        DO i = vdims%i_start,vdims%i_end
          dv_gwd(i,j,k) = R_v(i,j,k) - dv_gwd(i,j,k)
        END DO
      END DO
    END DO

  END IF ! end if over l_retain_gwd_tendencies
  
  ! ----------------------------------------------------------------------
  ! Section GWD.3 Call GWD diagnostics
  ! ----------------------------------------------------------------------
  SELECT CASE (model_type)
  CASE DEFAULT

    IF (sf(0,6) .OR. l_make_pws_prereq) THEN
      ! diagnostics requested this timestep
      CALL diagnostics_gwd(                                         &
                         row_length, rows, model_levels             &
  ,                      n_rows, u_rows, v_rows                     &
  ,                      off_x, off_y                               &
  ,                      at_extremity                               &
  ,                      u, v, R_u, R_v, T_inc                      &
  ,                      u_incr_diagnostic,v_incr_diagnostic        &
  ,                      t_incr_diagnostic                          &
  ,                      exner_theta_levels                         &
  ,                      stress_ud     ,  stress_vd                 &
  ,                      stress_ud_satn,  stress_vd_satn            &
  ,                      stress_ud_wake,  stress_vd_wake            &
  ,                      du_dt_satn    ,  dv_dt_satn                &
  ,                      du_dt_wake   ,  dv_dt_wake                 &
  ,                      u_s_d, v_s_d, nsq_s_d                      &
  ,                      num_lim_d, num_fac_d                       &
  ,                      fr_d, bld_d, bldt_d                        &
  , tausx_d, tausy_d, taus_scale_d                                  &
  , sd_orog_land, land_sea_mask, land_points                        &
  , orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
  , orog_slope_d, orog_anis_d, orog_dir_d                           &
  , gwspec_eflux, gwspec_sflux, gwspec_wflux                        &
  , gwspec_nflux, gwspec_ewacc, gwspec_nsacc,                       &
   STASHwork                                                        &
   )
    END IF            ! on sf(0,6)

  CASE (mt_single_column)

    IF (L_SCMDiags(SCMDiag_gwd)) THEN
      !
      !-----------------------------------------------------------------------
      !     SCM GWD Package
      !-----------------------------------------------------------------------
      !       Stash 6,201
      CALL scmoutput(stress_ud,'gw_stress_u',                         &
           'Mtn Drag: full u stress','N/m^2',                         &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,202
      CALL scmoutput(stress_vd,'gw_stress_v',                         &
           'Mtn Drag: full v stress','N/m^2',                         &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,223
      CALL scmoutput(stress_ud_satn,'gw_satn_stress_u',               &
           'Mtn Drag: gw u stress','N/m^2',                           &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,224
      CALL scmoutput(stress_vd_satn,'gw_satn_stress_v',               &
           'Mtn Drag: gw v stress','N/m^2',                           &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,227
      CALL scmoutput(stress_ud_wake,'gw_wake_stress_u',               &
           'Mtn Drag: wake u stress','N/m^2',                         &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,228
      CALL scmoutput(stress_vd_wake,'gw_wake_stress_v',               &
           'Mtn Drag: wake v stress','N/m^2',                         &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,207
      CALL scmoutput(du_dt_satn,'gw_satn_acc_u',                      &
           'Mtn Drag: gw u acc','m/s^2',                              &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,208
      CALL scmoutput(dv_dt_satn,'gw_satn_acc_v',                      &
           'Mtn Drag: gw v acc','m/s^2',                              &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,231
      CALL scmoutput(du_dt_wake,'gw_wake_acc_u',                      &
           'Mtn Drag: wake u acc','m/s^2',                            &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,232
      CALL scmoutput(dv_dt_wake,'gw_wake_acc_v',                      &
           'Mtn Drag: wake v acc','m/s^2',                            &
           t_inst,d_all,default_streams,'',routinename)
      !       Stash 6,214
      CALL scmoutput(u_s_d,'gw_u_in',                                 &
           'Mtn Drag: low-level u','m/s',                             &
           t_inst,d_sl,default_streams,'',routinename)
      !       Stash 6,215
      CALL scmoutput(v_s_d,'gw_v_in',                                 &
           'Mtn Drag: low-level v','m/s',                             &
           t_inst,d_sl,default_streams,'',routinename)
      !       Stash 6,216
      CALL scmoutput(nsq_s_d,'gw_nsq_in',                             &
           'Mtn Drag: low-level nsq','/m^2',                          &
           t_inst,d_sl,default_streams,'',routinename)
      !       Stash 6,217
      CALL scmoutput(fr_d,'gw_fr',                                    &
           'Mtn Drag: low-level froude num','',                       &
           t_inst,d_sl,default_streams,'',routinename)
      !       Stash 6,218
      CALL scmoutput(bld_d,'gw_bl_depth',                             &
           'Mtn Drag: blocked layer depth','m',                       &
           t_inst,d_sl,default_streams,'',routinename)

    END IF ! SCMDiag_gwd
!
  END SELECT ! model_type

  IF ( ALLOCATED( u_incr_diagnostic ) ) DEALLOCATE ( u_incr_diagnostic )
  IF ( ALLOCATED( v_incr_diagnostic ) ) DEALLOCATE ( v_incr_diagnostic )
  IF ( ALLOCATED( t_incr_diagnostic ) ) DEALLOCATE ( t_incr_diagnostic )
  IF ( ALLOCATED( stress_ud         ) ) DEALLOCATE ( stress_ud         )
  IF ( ALLOCATED( stress_vd         ) ) DEALLOCATE ( stress_vd         )
  IF ( ALLOCATED( du_dt_satn        ) ) DEALLOCATE ( du_dt_satn        )
  IF ( ALLOCATED( dv_dt_satn        ) ) DEALLOCATE ( dv_dt_satn        )
  IF ( ALLOCATED( u_s_d             ) ) DEALLOCATE ( u_s_d             )
  IF ( ALLOCATED( v_s_d             ) ) DEALLOCATE ( v_s_d             )
  IF ( ALLOCATED( nsq_s_d           ) ) DEALLOCATE ( nsq_s_d           )
  IF ( ALLOCATED( fr_d              ) ) DEALLOCATE ( fr_d              )
  IF ( ALLOCATED( bld_d             ) ) DEALLOCATE ( bld_d             )
  IF ( ALLOCATED( bldt_d            ) ) DEALLOCATE ( bldt_d            )
  IF ( ALLOCATED( num_lim_d         ) ) DEALLOCATE ( num_lim_d         )
  IF ( ALLOCATED( num_fac_d         ) ) DEALLOCATE ( num_fac_d         )
  IF ( ALLOCATED( stress_ud_satn    ) ) DEALLOCATE ( stress_ud_satn    )
  IF ( ALLOCATED( stress_vd_satn    ) ) DEALLOCATE ( stress_vd_satn    )
  IF ( ALLOCATED( stress_ud_wake    ) ) DEALLOCATE ( stress_ud_wake    )
  IF ( ALLOCATED( stress_vd_wake    ) ) DEALLOCATE ( stress_vd_wake    )
  IF ( ALLOCATED( du_dt_wake        ) ) DEALLOCATE ( du_dt_wake        )
  IF ( ALLOCATED( dv_dt_wake        ) ) DEALLOCATE ( dv_dt_wake        )
  IF ( ALLOCATED( tausx_d           ) ) DEALLOCATE ( tausx_d           )
  IF ( ALLOCATED( tausy_d           ) ) DEALLOCATE ( tausy_d           )
  IF ( ALLOCATED( taus_scale_d      ) ) DEALLOCATE ( taus_scale_d      )
  IF ( ALLOCATED( orog_slope_d      ) ) DEALLOCATE ( orog_slope_d      )
  IF ( ALLOCATED( orog_anis_d       ) ) DEALLOCATE ( orog_anis_d       )
  IF ( ALLOCATED( orog_dir_d        ) ) DEALLOCATE ( orog_dir_d        )
  IF ( ALLOCATED( gwspec_eflux      ) ) DEALLOCATE ( gwspec_eflux      )
  IF ( ALLOCATED( gwspec_sflux      ) ) DEALLOCATE ( gwspec_sflux      )
  IF ( ALLOCATED( gwspec_wflux      ) ) DEALLOCATE ( gwspec_wflux      )
  IF ( ALLOCATED( gwspec_nflux      ) ) DEALLOCATE ( gwspec_nflux      )
  IF ( ALLOCATED( gwspec_ewacc      ) ) DEALLOCATE ( gwspec_ewacc      )
  IF ( ALLOCATED( gwspec_nsacc      ) ) DEALLOCATE ( gwspec_nsacc      )

END IF ! on error code equal to zero

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE NI_gwd_ctl

END MODULE NI_gwd_ctl_mod
