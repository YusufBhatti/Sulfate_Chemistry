! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Stochastic Perturbations of Tendencies (SPT) main
!               routine, it multiplies the physical tendencies
!               stored in physics tendencies_spt_mod times a
!               stochastic forcing pattern stored in fp_mod
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics

MODULE SPT_main_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SPT_MAIN_MOD'

CONTAINS

SUBROUTINE SPT_main(                                                    &
!in
        row_length, rows, n_rows, land_points,                          &
        exner_theta_levels, p_theta_levels, rho_in,                     &
        first_atmstep_call, theta_star, q_star,                         &
        land_sea_mask, sd_orog_land,                                    &
        stashwork35)

! Use subroutine for conservation if requested
USE SPT_conservation_mod, ONLY: SPT_conservation

! SPT gui settings passed in via NameList READ
USE stochastic_physics_run_mod, ONLY:                             &
    spt_toplev, spt_botlev, spt_bot_tap_lev, spt_top_tap_lev,     &
    l_spt_rain, l_spt_rad, l_spt_gwd, l_spt_conv, l_spt_conv_mom, &
    l_spt_cfl, nsmooth_spt,                                       &
    rain_std, rad_std, gwd_std, conv_std,                         &
    skeb2_up_flux, offx_spt, offy_spt, mask_smooth_spt,           &
    offx_spt, offy_spt, skeb2_up_flux,                            &
    l_spt_qcons, l_spt_mse,                                       &
    sd_orog_thres, psif_orog_thres 
    
! Sructure containing SPT physics diagnostics
USE spt_diag_mod, ONLY:                                           &
    strsptdiag

! SPT Forcing pattern field psif and psif2
USE fp_mod, ONLY: psif, psif2

! To define orography array
USE mask_compression, ONLY : expand_from_mask

! Get Atmos_physics 1 and 2 tendencies
USE physics_tendencies_mod,  ONLY:                                &
    dt_conv, dq_conv, du_conv, dv_conv,                           &
    dt_mic, dq_mic, dt_sw, dq_sw, dt_lw, dq_lw,                   &
    du_gwd,dv_gwd,                                                &
    ! get tendencies for stph
    dtheta_stph, du_stph, dv_stph,                                &
    l_retain_stph_tendencies
    
! Load SPT increment arrays
USE fields_rhs_mod,       ONLY:                                   &
   theta_spt, q_spt, r_u_spt, r_v_spt

! Get timestep in seconds
USE timestep_mod,      ONLY:  timestep

! Bounds of arrays
USE atm_fields_bounds_mod, ONLY:                                  &
    tdims, tdims_s,                                               &
    pdims, pdims_s,                                               &
    udims, vdims, array_dims

USE diagnostics_spt_mod, ONLY: diagnostics_spt

USE nlsizes_namelist_mod, ONLY: model_levels

USE umPrintMgr
USE UM_ParVars
USE Field_Types,     ONLY: fld_type_p, fld_type_u, fld_type_v
USE stash_array_mod, ONLY: sf
USE mpp_conf_mod,    ONLY: swap_field_is_vector, swap_field_is_scalar

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new => qsat_wat,                                 &
                    qsat_wat_mix_new => qsat_wat_mix,                         &
                    l_new_qsat_cntl !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

! DrHook modules
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
     row_length                                                         &
! local number of points on a row
 ,   rows                                                               &
! local number of rows for u
 ,   n_rows
! local number of rows for v

REAL, INTENT (IN) ::                                                    &
    exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                       tdims_s%j_start:tdims_s%j_end,                   &
                       tdims_s%k_start:tdims_s%k_end)                   &
! Exner pressure on theta levels (to convert T -> theta)
 ,   p_theta_levels(pdims_s%i_start:pdims_s%i_end,                      &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)                            &
! Pressure array
 ,   rho_in(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,                              &
            pdims_s%k_start:pdims_s%k_end)                              &
! density * square of radius of earth

 ,   theta_star(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)                              &
! theta after atmos_physics2


 ,   q_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                            &
              tdims%k_start:tdims%k_end)
! q after atmos_physics2

LOGICAL, INTENT(IN) ::                                                  &
     first_atmstep_call
! Is true for first step of: NRUN and each CRUN

REAL ::                                                                 &
     stashwork35(*)
! Array containing requested sec35 STASH diagnostics

REAL ::                                                                 &
     delta_theta(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end,                             &
                 1:tdims%k_end)                                         &
! Matrix holding the SPT increments for theta
 ,   delta_q(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             1:tdims%k_end)                                             &
! Matrix holding the SPT increments for q
 ,   delta_u(udims%i_start:udims%i_end,                                 &
             udims%j_start:udims%j_end,                                 &
             1:udims%k_end)                                             &
! Matrix holding the SPT increments for u
 ,   delta_v(vdims%i_start:vdims%i_end,                                 &
             vdims%j_start:vdims%j_end,                                 &
             1:vdims%k_end)                                             &
! Matrix holding the SPT increments for v
 ,   cfl_marker(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end)                              &
! Marks where the CFL criteria is breached
 ,   psif_halos(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)                          &
! Forcing pattern with haloes for U and V forcing
 ,   psif2_halos(pdims_s%i_start:pdims_s%i_end,                         &
                 pdims_s%j_start:pdims_s%j_end)                         &
! 2nd Forcing pattern with haloes for U and V forcing
,    dt_spt_mic(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                1:tdims%k_end)                                          &
! Internal SPT array for MIC T tend
,    dq_spt_mic(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                1:tdims%k_end)                                          &
! Internal SPT array for MIC q tend
,    dt_spt_lw(tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               1:tdims%k_end)                                           &
! Internal SPT array for LW T tend
,    dq_spt_lw(tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               1:tdims%k_end)                                           &
! Internal SPT array for LW q tend
,    dt_spt_sw(tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               1:tdims%k_end)                                           &
! Internal SPT array for SW T tend
,    dq_spt_sw(tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               1:tdims%k_end)                                           &
! Internal SPT array for SW q tend
,    du_spt_gwd(udims%i_start:udims%i_end,                              &
                udims%j_start:udims%j_end,                              &
                udims%k_start:udims%k_end)                              &
! Internal SPT array for GWD u tend
,    dv_spt_gwd(vdims%i_start:vdims%i_end,                              &
                vdims%j_start:vdims%j_end,                              &
                vdims%k_start:vdims%k_end)                              &
! Internal SPT array for GWD v tend

,    dt_spt_conv(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end,                             &
                 1:tdims%k_end)                                         &
! Internal SPT array for conv T tend
,    dq_spt_conv(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end,                             &
                 1:tdims%k_end)                                         &
! Internal SPT array for conv q tend
,    du_spt_conv(udims%i_start:udims%i_end,                             &
                 udims%j_start:udims%j_end,                             &
                 udims%k_start:udims%k_end)                             &
! Internal SPT array for conv u tend
,    dv_spt_conv(vdims%i_start:vdims%i_end,                             &
                 vdims%j_start:vdims%j_end,                             &
                 vdims%k_start:vdims%k_end)
! Internal SPT array for conv v tend

REAL ,ALLOCATABLE ::                                                    &
     delta_theta0(:,:,:)                                                &
 ! holding variable for smoothing of theta incr
 ,   delta_q0(:,:,:)                                                    &
 ! holding variable for smoothing of q incr
 ,   delta_u0(:,:,:)                                                    &
 ! holding variable for smoothing of u incr
 ,   delta_v0(:,:,:)
 ! holding variable for smoothing of v incr

INTEGER ::                                                              &
     i                                                                  &
 ! loop index over x direction
 ,   j                                                                  &
 ! loop index over y direction
 ,   k                                                                  &
 ! loop index over z direction
 ,   ismooth                                                            &
 ! loop index over smoothing iterations
 ,   ii, jj
 ! Integers to convert conv. mom increments from  P grid to U grid.

INTEGER, PARAMETER ::                                                   &
     zero  = 0
             ! used for identifying zero'th PE

! ++++ Orographic capping variables
INTEGER ::                                                              &

     land_points   ! Number of land points

REAL ::                                                                 &
    sd_orog_land(land_points)        
    ! Standard deviation on land points
    
REAL, ALLOCATABLE, SAVE ::                                              &
    sd_orog_spt(:,:)                                                    &
    ! Standard deviation on model gridpoints
,   sd_orog_spt_haloes(:,:)             
    ! Standard deviation on model gridpoints with haloes

LOGICAL, INTENT(IN) ::                                                  &
    land_sea_mask( pdims%i_start : pdims%i_end,                         &
                   pdims%j_start : pdims%j_end )
         ! land_sea_mask

! ++++ End of orographic capping variables

TYPE (array_dims), SAVE::  t_sptdims_l,q_sptdims_l,                     &
                           u_sptdims_l,v_sptdims_l
     ! array to store vort dimensions


! Non-haloed of cloud arrays for pc2_checks
REAL, ALLOCATABLE ::                                                    &
     t_inc(:,:,:)                                                       &
   ! Temperature (K) before timestep to compute qs.
 ,   p_theta_levels_nh(:,:,:)
! Convert levels into a real domain to be used as the multiplicative
! factor to ramp up (or down) the delta (perturbation added)

REAL, SAVE ::                                                           &
     spt_top_tap_lev_real, spt_bot_tap_lev_real,spt_botlev_real,        &
     spt_toplev_real

REAL :: k_real,spt_amp_factor
 ! for k index (UM levels) and the amplification factor for the ramp up.


REAL ::                                                                 &
     qsl_t(row_length,rows)
     ! Saturated specific humidity for dry bulb temperature T

TYPE (strsptdiag) :: spt_diag
!     Declaration of Stochastic Physics diagnostics.
!     Using 35 (Stochastic diagnostics diagnostics)
!     SPT array defined in spt_diag_mod

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPT_MAIN'

! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Error trap if this is a LAM configuration (SPT should not be called)
IF (model_type /= mt_global) THEN
  WRITE(umMessage,'("**ERROR**: SPT not available in a Limited Area Model")')
  CALL umPrint(umMessage,src='SPT_main')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='SPT_main')
END IF

IF (spt_toplev >= model_levels) THEN
  WRITE(umMessage,'("**ERROR**: SPT TOP LEVEL = OR > MODEL LEVELS")')
  CALL umPrint(umMessage,src='SPT_main')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='SPT_main')
END IF

IF (spt_botlev <= 1) THEN
  WRITE(umMessage,'("ERROR: SPT BOTTOM LEVEL = OR < MODEL BASE")')
  CALL umPrint(umMessage,src='SPT_main')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='SPT_main')
END IF

! ------------------------------------------------------------------
!      Define STASH arrays
! ------------------------------------------------------------------

! Prepare STASH diagnostics where requested
! Allocate STASH diagnostic arrays using Structure TYPE spt_diag

spt_diag%l_spt_forcing_pattern                = sf(23,35)
! Include SPT increments as well
spt_diag%l_spt_theta_inc                      = sf(24,35)
spt_diag%l_spt_q_inc                          = sf(25,35)
spt_diag%l_spt_u_inc                          = sf(26,35)
spt_diag%l_spt_v_inc                          = sf(27,35)
! CFL criteria array
spt_diag%l_cfl_br_marker                      = sf(28,35)

! Temp Inc for budget
spt_diag%l_spt_t_inc                          = sf(29,35)

IF (sf(0,35)) THEN

  ! SPT forcing pattern allocation or re-initialization
  IF (spt_diag%l_spt_forcing_pattern) THEN
    ALLOCATE(spt_diag%spt_forcing_pattern(pdims%i_start:pdims%i_end,    &
                                          pdims%j_start:pdims%j_end,    &
                                          pdims%k_start:pdims%k_end))
  ELSE
    ! Code to allocate unity arrays when not
    ! used (for portability)
    ALLOCATE(spt_diag%spt_forcing_pattern(1,1,1))
    spt_diag%spt_forcing_pattern(:,:,:) = 0.0
  END IF

  ! ----- Increments
   ! Theta Inc SPT forcing allocation or re-initialization
  IF (spt_diag%l_spt_theta_inc) THEN
    ALLOCATE(spt_diag%spt_theta_inc(tdims%i_start:tdims%i_end,          &
                                    tdims%j_start:tdims%j_end,          &
                                    1:tdims%k_end))
  ELSE
    ALLOCATE(spt_diag%spt_theta_inc(1,1,1))
    spt_diag%spt_theta_inc(:,:,:) = 0.0
  END IF

  ! q Inc SPT forcing allocation or re-initialization
  IF (spt_diag%l_spt_q_inc) THEN
    ALLOCATE(spt_diag%spt_q_inc(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                1:tdims%k_end))
  ELSE
    ALLOCATE(spt_diag%spt_q_inc(1,1,1))
    spt_diag%spt_q_inc(:,:,:) = 0.0
  END IF

  ! u Inc SPT forcing allocation or re-initialization
  IF (spt_diag%l_spt_u_inc) THEN
    ALLOCATE(spt_diag%spt_u_inc(udims%i_start:udims%i_end,              &
                                udims%j_start:udims%j_end,              &
                                1:udims%k_end))
  ELSE
    ALLOCATE(spt_diag%spt_u_inc(1,1,1))
    spt_diag%spt_u_inc(:,:,:) = 0.0
  END IF

  ! v Inc SPT forcing allocation or re-initialization
  IF (spt_diag%l_spt_v_inc) THEN
    ALLOCATE(spt_diag%spt_v_inc(vdims%i_start:vdims%i_end,              &
                                vdims%j_start:vdims%j_end,              &
                                1:vdims%k_end))
  ELSE
    ALLOCATE(spt_diag%spt_v_inc(1,1,1))
    spt_diag%spt_v_inc(:,:,:) = 0.0
  END IF

    ! CFL criteria
  IF (spt_diag%l_cfl_br_marker) THEN
    ALLOCATE(spt_diag%cfl_br_marker(tdims%i_start:tdims%i_end,          &
                                    tdims%j_start:tdims%j_end))
  ELSE
    ALLOCATE(spt_diag%cfl_br_marker(1,1))
    spt_diag%cfl_br_marker(:,:) = 0.0
  END IF

  ! T Inc SPT forcing allocation or re-initialization
  IF (spt_diag%l_spt_t_inc) THEN
    ALLOCATE(spt_diag%spt_t_inc(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                1:tdims%k_end))
  ELSE
    ALLOCATE(spt_diag%spt_t_inc(1,1,1))
    spt_diag%spt_t_inc(:,:,:) = 0.0
  END IF

END IF ! SF(0,35)

! --------------------------------------------------------------------
! copy psif and "before SPT forcing" variables on the std_diag array
! --------------------------------------------------------------------

! SPT forcing pattern
IF (spt_diag%l_spt_forcing_pattern) THEN
  DO k =  pdims%k_start, pdims%k_end
    DO j =  pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        spt_diag%spt_forcing_pattern(i,j,k) = psif(i,j,k)
      END DO
    END DO
  END DO
END IF

! --------------------------------------------------------------------
! Save tendencies into internal SPT arrays
! --------------------------------------------------------------------

! rain
dt_spt_mic = 0.0
dq_spt_mic = 0.0

IF (l_spt_rain) THEN
  IF (.NOT. l_spt_mse) THEN
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dt_spt_mic(i,j,k)= dt_mic(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dq_spt_mic(i,j,k)= dq_mic(i,j,k)
      END DO
    END DO
  END DO
END IF

! radiation
dt_spt_lw = 0.0
dq_spt_lw = 0.0
dt_spt_sw = 0.0
dq_spt_sw = 0.0

IF (l_spt_rad) THEN
  IF (.NOT. l_spt_mse) THEN
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dt_spt_sw(i,j,k)= dt_sw(i,j,k)
          dt_spt_lw(i,j,k)= dt_lw(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dq_spt_sw(i,j,k)= dq_sw(i,j,k)
        dq_spt_lw(i,j,k)= dq_lw(i,j,k)
      END DO
    END DO
  END DO
END IF

! GWD
du_spt_gwd = 0.0
dv_spt_gwd = 0.0

IF (l_spt_gwd) THEN
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_spt_gwd(i,j,k)= du_gwd(i,j,k)
      END DO
    END DO
  END DO

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_spt_gwd(i,j,k)= dv_gwd(i,j,k)
      END DO
    END DO
  END DO

END IF

! convection
dt_spt_conv = 0.0
dq_spt_conv = 0.0

IF (l_spt_conv) THEN
  IF (.NOT. l_spt_mse) THEN
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dt_spt_conv(i,j,k)= dt_conv(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dq_spt_conv(i,j,k)= dq_conv(i,j,k)
      END DO
    END DO
  END DO
END IF

! convection momentum
du_spt_conv = 0.0
dv_spt_conv = 0.0

IF (l_spt_conv_mom) THEN
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_spt_conv(i,j,k)= du_conv(i,j,k)
      END DO
    END DO
  END DO

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_spt_conv(i,j,k)= dv_conv(i,j,k)
      END DO
    END DO
  END DO

END IF

! ---------------------------------------------------------------
! Apply CFL criteria to avoid convection instabilities.
! ---------------------------------------------------------------
! Initialize cfl_marker
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    cfl_marker(i,j) = 1.0
  END DO
END DO

IF (l_spt_cfl) THEN

  ! Compute CFL criteria and mark points where it is breached to
  ! 1 in the cfl_marker
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      DO k = spt_bot_tap_lev+1, spt_top_tap_lev-1
        IF ( skeb2_up_flux(i,j,k)*( 1 + conv_std*psif2(i,j) )*timestep  &
             >= ABS(p_theta_levels(i,j,k) - p_theta_levels(i,j,k-1)) )  &
             cfl_marker(i,j)=0.0
      END DO
    END DO
  END DO

  ! Apply CFL mask to the delta_x variables
  DO k= spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dt_spt_conv(i,j,k) = cfl_marker(i,j) * dt_spt_conv(i,j,k)
        dq_spt_conv(i,j,k) = cfl_marker(i,j) * dq_spt_conv(i,j,k)
      END DO
    END DO
  END DO

END IF ! end if l_spt_cfl

  ! CFL criteria outputted to STASH
IF (spt_diag%l_cfl_br_marker) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      spt_diag%cfl_br_marker(i,j) = cfl_marker(i,j)
    END DO
  END DO
END IF

! --------------------------------------------------------------------
! Compute SPT increments (tendencies x forcing pattern)
! --------------------------------------------------------------------

 ! First get psif with haloes if u and v are requested from GWD
IF (l_spt_gwd) THEN

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        psif_halos(i,j,k) = psif(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( psif_halos,                                         &
                 pdims_s%i_len - 2*pdims_s%halo_i,                      &
                 pdims_s%j_len - 2*pdims_s%halo_j,                      &
                 pdims_s%k_len,                                         &
                 pdims_s%halo_i, pdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

ELSE
  psif_halos(:,:,:) = 0.
END IF

 ! First get psif2 if u and v from convection are in.
IF (l_spt_conv_mom) THEN

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      psif2_halos(i,j) = psif2(i,j)
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( psif2_halos,                                        &
                 pdims_s%i_len - 2*pdims_s%halo_i,                      &
                 pdims_s%j_len - 2*pdims_s%halo_j,                      &
                 1, pdims_s%halo_i, pdims_s%halo_j,                     &
                 fld_type_p,swap_field_is_scalar)

ELSE
    psif2_halos(:,:) = 0.
END IF

! +++++++++++ theta
! only apply increments if MSE conservation is false
IF (.NOT. l_spt_mse) THEN 
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        delta_theta(i,j,k)= ( psif(i,j,k)*(                             &
                              rain_std*dt_spt_mic(i,j,k) + rad_std*(    &
                              dt_spt_sw(i,j,k) + dt_spt_lw(i,j,k) ) ) + &
                              psif2(i,j)*conv_std* dt_spt_conv(i,j,k) ) &
                              / exner_theta_levels(i,j,k)


      END DO    
    END DO
  END DO

END IF

! Set field below and above tapering to zero
DO k = 1,spt_bot_tap_lev-1
  DO j =  tdims%j_start, tdims%j_end
    DO i =  tdims%i_start, tdims%i_end
      delta_theta(i,j,k)=0.0
    END DO
  END DO
END DO

DO k = spt_top_tap_lev+1,tdims%k_end
  DO j =  tdims%j_start, tdims%j_end
    DO i =  tdims%i_start, tdims%i_end
      delta_theta(i,j,k)=0.0
    END DO
  END DO
END DO

! ++++++++++ q
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      delta_q(i,j,k)= psif(i,j,k)*(rain_std*dq_spt_mic(i,j,k) +         &
                      rad_std*(dq_spt_sw(i,j,k) + dq_spt_lw(i,j,k))) +  &
                      psif2(i,j)*conv_std*dq_spt_conv(i,j,k)

    END DO    
  END DO
END DO

DO k = 1,spt_bot_tap_lev-1
  DO j = tdims%j_start, tdims%j_end
    DO i =  tdims%i_start, tdims%i_end
      delta_q(i,j,k)=0.0
    END DO
  END DO
END DO

DO k = spt_top_tap_lev+1,tdims%k_end
  DO j =  tdims%j_start, tdims%j_end
    DO i =  tdims%i_start, tdims%i_end
      delta_q(i,j,k)=0.0
    END DO
  END DO
END DO

! ++++++++++ u
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

      delta_u(i,j,k) = ( psif_halos(i,j,k) + psif_halos(i+1,j,k) )*0.5  &
                       *gwd_std*du_spt_gwd(i,j,k)  +                    &
                       ( psif2_halos(i,j) + psif2_halos(i+1,j))*0.5     &
                       *conv_std* du_spt_conv(i,j,k)

    END DO
  END DO
END DO

DO k = 1,spt_bot_tap_lev-1
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      delta_u(i,j,k)=0.0
    END DO
  END DO
END DO

DO k = spt_top_tap_lev+1,udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      delta_u(i,j,k)=0.0
    END DO
  END DO
END DO

! +++++++++++ V
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end

      delta_v(i,j,k)= ( psif_halos(i,j,k) + psif_halos(i,j+1,k) )* 0.5  &
                      *gwd_std*dv_spt_gwd(i,j,k)  +                     &
                      ( psif2_halos(i,j) + psif2_halos(i,j+1) )*0.5     &
                      *conv_std* dv_spt_conv(i,j,k)
    END DO
  END DO
END DO

DO k = 1,spt_bot_tap_lev-1
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      delta_v(i,j,k)=0.0
    END DO
  END DO
END DO

DO k = spt_top_tap_lev+1,vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      delta_v(i,j,k)=0.0
    END DO
  END DO
END DO

! --------------------------------------------------------------------
! RAMP UP/DOWN SPT increments for selected areas
! --------------------------------------------------------------------

! +++++++++ DEFINE THE TOP AND BOTTOM LID AND RAMP UP/DOWN LEVELS

! Settings required only at the first time-step
IF (first_atmstep_call) THEN
  ! Declare all the real numbers for the levels
  spt_top_tap_lev_real=REAL(spt_top_tap_lev)
  spt_bot_tap_lev_real=REAL(spt_bot_tap_lev)
  spt_botlev_real=REAL(spt_botlev)
  spt_toplev_real=REAL(spt_toplev)
END IF


! ++++++ Ramp up from BL

! +++++++++++ theta
DO k = spt_bot_tap_lev, spt_botlev
  k_real=REAL(k)
  spt_amp_factor=(k_real-spt_bot_tap_lev_real)/                         &
  (spt_botlev_real-spt_bot_tap_lev_real)
  
  IF (.NOT. l_spt_mse) THEN 
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        delta_theta(i,j,k)=spt_amp_factor*delta_theta(i,j,k)
      END DO
    END DO
  END IF 
  ! ++++++++++ q
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      delta_q(i,j,k)=spt_amp_factor*delta_q(i,j,k)
    END DO
  END DO

  ! ++++++++++ U
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      delta_u(i,j,k)= spt_amp_factor*delta_u(i,j,k)
    END DO
  END DO

  ! +++++++++++ V
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      delta_v(i,j,k)= spt_amp_factor*delta_v(i,j,k)
    END DO
  END DO
END DO !End loop over levels k

! ++++++ Ramp down when getting close to tropopause

! +++++++++++ theta
DO k = spt_toplev, spt_top_tap_lev
  k_real=REAL(k)
  spt_amp_factor=(spt_top_tap_lev_real-k_real)/                         &
  (spt_top_tap_lev_real-spt_toplev_real)


  IF (.NOT. l_spt_mse) THEN 
    DO j =tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        delta_theta(i,j,k)=spt_amp_factor*delta_theta(i,j,k)
      END DO
    END DO
  END IF
  
  ! ++++++++++ q
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      delta_q(i,j,k)=spt_amp_factor*delta_q(i,j,k)
    END DO
  END DO

  ! ++++++++++ U
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      delta_u(i,j,k)= spt_amp_factor*delta_u(i,j,k)
    END DO
  END DO

  ! +++++++++++ V
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      delta_v(i,j,k)= spt_amp_factor*delta_v(i,j,k)
    END DO
  END DO

END DO ! End loop over levels.

! --------------------------------------------------------------------
! Remove Theta and q if q_star+delta_q > q_sat
! --------------------------------------------------------------------


 ! Set the non-haloed pressure levels matrix
IF (.NOT. ALLOCATED(p_theta_levels_nh)) THEN
  ALLOCATE(p_theta_levels_nh(pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,                 &
                             pdims%k_start:pdims%k_end))
END IF

! Set the non-haloed T increments and copy the
! current values * exner pressure
IF (.NOT. ALLOCATED(t_inc)) THEN
  ALLOCATE(t_inc(tdims%i_start:tdims%i_end,                             &
                 tdims%j_start:tdims%j_end,                             &
                 tdims%k_start:tdims%k_end))
END IF

! Transform theta values onto T values of T_inc
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      t_inc(i,j,k)= theta_star(i,j,k) * exner_theta_levels(i,j,k)
    END DO
  END DO
END DO

!copy pressure values array into non-haloed array "_nh"
DO k = pdims%k_start, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_theta_levels_nh(i,j,k)= p_theta_levels(i,j,k)
    END DO
  END DO
END DO

! (set delta_q and delta_theta to 0.0 at this gridpoint).
DO k = spt_bot_tap_lev, spt_top_tap_lev

  IF (l_new_qsat_cntl) THEN
    IF (l_mr_physics) THEN
      CALL qsat_wat_mix_new(qsl_t, t_inc(:,:,k), p_theta_levels_nh(:,:,k),    &
                    row_length,rows)
    ELSE
      CALL qsat_wat_new(qsl_t, t_inc(:,:,k), p_theta_levels_nh(:,:,k),        &
                    row_length,rows)
    END IF
  ELSE
    ! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl_t, t_inc(1,1,k), p_theta_levels_nh(1,1,k),          &
                      row_length*rows, l_mr_physics)
  END IF

  ! If q >q_sat  or q <0 -> Set the SPT increments in q and theta to 0.
  DO j = 1, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF ( q_star(i,j,k) + delta_q(i,j,k) >= qsl_t(i,j) .OR.            &
           q_star(i,j,k) + delta_q(i,j,k) <0.0) THEN

        delta_q(i,j,k) = 0.0
        delta_theta(i,j,k) = 0.0

      END IF
    END DO
  END DO
END DO  ! End loop over levels k

! Deallocate non-haloed variables
DEALLOCATE(t_inc)
DEALLOCATE(p_theta_levels_nh)

! ----------------------------------------------------------------
! Apply an extra stability constraint: 
! If the standard deviation over orography is higher 
! than sd_orog_thres and abs(psif) > psif_orog_thres 
!
! Then kill the SPT perturbations over these points. As they 
! try to act against the GWD stabilizing tendencies.
! ----------------------------------------------------------------

IF(first_atmstep_call) THEN

  IF (.NOT. ALLOCATED(sd_orog_spt)) THEN
    ALLOCATE(sd_orog_spt(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end))
  ENDIF

! Create orography mask
  CALL expand_from_mask(sd_orog_spt, sd_orog_land, land_sea_mask,       &
                        pdims%i_end*pdims%j_end, land_points)
                        
! Copy orography mask to array with haloes (initialise first)
  IF (.NOT. ALLOCATED(sd_orog_spt_haloes)) THEN
    ALLOCATE(sd_orog_spt_haloes(tdims_s%i_start:tdims_s%i_end,          &
                                tdims_s%j_start:tdims_s%j_end))
  ENDIF
  
  
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Avoid copying the values over the sea (set to -1E9).
      IF (sd_orog_spt(i,j) >= 0.0) THEN
        sd_orog_spt_haloes(i,j) = sd_orog_spt(i,j)
      ELSE
        sd_orog_spt_haloes(i,j) = 0.0
      END IF
    END DO
  END DO
  
! DEPENDS ON: swap_bounds
  CALL swap_bounds(sd_orog_spt_haloes,                                  &
                 tdims_s%i_end - tdims_s%i_start+1 - 2*tdims_s%halo_i,  &
                 tdims_s%j_end - tdims_s%j_start+1 - 2*tdims_s%halo_j,  &
                 1, tdims_s%halo_i, tdims_s%halo_j,                     &
                 fld_type_p,swap_field_is_scalar)

! End if over 1st atm step call.
END IF

! --- Apply stability constraint
! For theta
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
  ! Remove theta increments
    IF ( sd_orog_spt(i,j) >= sd_orog_thres .AND.                        &
         ABS(psif(i,j,spt_botlev)) >= psif_orog_thres) THEN    
      delta_theta(i,j,:) = 0.0
    END IF
    
  END DO
END DO 

! For q
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
  ! Remove humidity increments
    IF ( sd_orog_spt(i,j) >= sd_orog_thres .AND.                        &
         ABS(psif(i,j,spt_botlev)) >= psif_orog_thres) THEN    
      delta_q(i,j,:) = 0.0
    END IF
    
  END DO
END DO  

! For u
DO j = udims%j_start, udims%j_end
  DO i = udims%i_start, udims%i_end
  ! Remove zonal velocity increments
    IF ( sd_orog_spt_haloes(i,j) >= sd_orog_thres .AND.                 &
        ABS(psif_halos(i,j,spt_botlev)) >= psif_orog_thres) THEN    
      delta_u(i,j,:) = 0.0
    END IF
  END DO
END DO  

! For v
DO j = vdims%j_start, vdims%j_end
  DO i = vdims%i_start, vdims%i_end
  ! Remove meridional velocity increments
    IF ( sd_orog_spt_haloes(i,j) >= sd_orog_thres .AND.                 &
        ABS(psif_halos(i,j,spt_botlev)) >= psif_orog_thres) THEN    
      delta_v(i,j,:) = 0.0
    END IF
  END DO
END DO  

! ------------------------------------------------------------------------
! Smooth delta function before going into the star (increments) variables
! -----------------------------------------------------------------------
IF (nsmooth_spt > 0) THEN


  ! ++++++++ theta !!!
  IF (.NOT. l_spt_mse) THEN 
     ! Declare arrays to save the original field
     ! Advanced smoothing option: using pre-defined smoothing array

    !create array for delta_theta0
    IF (first_atmstep_call) THEN
      t_sptdims_l%i_start = tdims%i_start - offx_spt
      t_sptdims_l%i_end   = tdims%i_end   + offx_spt
      t_sptdims_l%j_start = tdims%j_start - offy_spt
      t_sptdims_l%j_end   = tdims%j_end   + offy_spt
      t_sptdims_l%k_start = 1
      t_sptdims_l%k_end   = tdims%k_end
    END IF

    !  Allocate smoothing buffer array
    IF (.NOT. ALLOCATED(delta_theta0)) THEN
      ALLOCATE(delta_theta0(t_sptdims_l%i_start:t_sptdims_l%i_end,      &
                            t_sptdims_l%j_start:t_sptdims_l%j_end,      &
                            t_sptdims_l%k_start:t_sptdims_l%k_end))
    END IF

    ! Pass in the values to the THETA haloed variable to do the smoothing
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      DO j =  tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          delta_theta0(i,j,k) = delta_theta(i,j,k)
        END DO
      END DO
    END DO !end loop over k

    !Initialize delta_theta for smoothing
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          delta_theta(i,j,k)=0.0
        END DO
      END DO
    END DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds( delta_theta0,row_length,  rows,                   &
                      model_levels, offx_spt, offy_spt, fld_type_p,     &
                      swap_field_is_scalar )

    !  2D 1-2-1 spatial filter (applied nsmooth times)
    DO k = spt_bot_tap_lev, spt_top_tap_lev
      !  Apply smoother
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          DO jj = -offy_spt, offy_spt
            DO ii = -offx_spt, offx_spt
              delta_theta(i,j,k) = delta_theta(i,j,k)        +          &
                                   mask_smooth_spt(ii,jj) *             &
                                   delta_theta0(ii+i,jj+j,k)
            END DO  ! ii
          END DO  ! jj
        END DO  ! i
      END DO  ! j
    END DO  ! k

    ! Deallocate delta_theta0
    IF (ALLOCATED(delta_theta0)) DEALLOCATE(delta_theta0)

  END IF    ! End loop if MSE conservation is not selected.

  ! ++++++++ q !!!
     !create array for delta_theta0
  IF (first_atmstep_call) THEN
    q_sptdims_l%i_start = tdims%i_start - offx_spt
    q_sptdims_l%i_end   = tdims%i_end   + offx_spt
    q_sptdims_l%j_start = tdims%j_start - offy_spt
    q_sptdims_l%j_end   = tdims%j_end   + offy_spt
    q_sptdims_l%k_start = 1
    q_sptdims_l%k_end   = tdims%k_end
  END IF

    !  Allocate smoothing buffer array
  IF (.NOT. ALLOCATED(delta_q0)) THEN
    ALLOCATE(delta_q0(q_sptdims_l%i_start:q_sptdims_l%i_end,            &
                      q_sptdims_l%j_start:q_sptdims_l%j_end,            &
                      q_sptdims_l%k_start:q_sptdims_l%k_end))
  END IF


  ! Pass in the values to the THETA haloed variable to do the smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j =  tdims%j_start, tdims%j_end
      DO i =  tdims%i_start,tdims%i_end
        delta_q0(i,j,k) = delta_q(i,j,k)
      END DO
    END DO
  END DO !end loop over k

  !Initialize delta_q for smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j =  tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        delta_q(i,j,k)=0.0
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( delta_q0, row_length,  rows,                        &
                    model_levels, offx_spt, offy_spt, fld_type_p,       &
                    swap_field_is_scalar  )

  !  2D 1-2-1 spatial filter (applied nsmooth times)
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    !  Apply smoother
    DO j =  tdims%j_start, tdims%j_end
      DO i =  tdims%i_start, tdims%i_end
        DO jj = -offy_spt, offy_spt
          DO ii = -offx_spt, offx_spt
            delta_q(i,j,k) = delta_q(i,j,k)       +                     &
                           mask_smooth_spt(ii,jj) * delta_q0(ii+i,jj+j,k)
          END DO  ! ii
        END DO  ! jj
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! Deallocate delta_q0
  IF (ALLOCATED(delta_q0)) DEALLOCATE(delta_q0)

  ! ++++++++ u !!!

     !create array for delta_theta0
  IF (first_atmstep_call) THEN
    u_sptdims_l%i_start = udims%i_start - offx_spt
    u_sptdims_l%i_end   = udims%i_end   + offx_spt
    u_sptdims_l%j_start = udims%j_start - offy_spt
    u_sptdims_l%j_end   = udims%j_end   + offy_spt
    u_sptdims_l%k_start = 1
    u_sptdims_l%k_end   = udims%k_end
  END IF

    !  Allocate smoothing buffer array
  IF (.NOT. ALLOCATED(delta_u0)) THEN
    ALLOCATE(delta_u0(u_sptdims_l%i_start:u_sptdims_l%i_end,            &
                       u_sptdims_l%j_start:u_sptdims_l%j_end,           &
                       u_sptdims_l%k_start:u_sptdims_l%k_end))
  END IF

  ! Pass in the values to the THETA haloed variable to do the smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        delta_u0(i,j,k) = delta_u(i,j,k)
      END DO
    END DO
  END DO !end loop over k

  ! Initialize delta_u for smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        delta_u(i,j,k) =0.0
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( delta_u0, row_length, rows,                         &
                    model_levels, offx_spt, offy_spt, fld_type_u,       &
                    swap_field_is_vector  )

  !  2D 1-2-1 spatial filter (applied nsmooth times)
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    !  Apply smoother
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        DO jj = -offy_spt, offy_spt
          DO ii = -offx_spt, offx_spt
            delta_u(i,j,k) = delta_u(i,j,k)    +                        &
                  mask_smooth_spt(ii,jj) * delta_u0(ii+i,jj+j,k)
          END DO  ! ii
        END DO  ! jj
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! Deallocate delta_q0
  IF (ALLOCATED(delta_u0)) DEALLOCATE(delta_u0)

   ! ++++++++ v !!!

  !create array for delta_theta0
  IF (first_atmstep_call) THEN
    v_sptdims_l%i_start = vdims%i_start - offx_spt
    v_sptdims_l%i_end   = vdims%i_end   + offx_spt
    v_sptdims_l%j_start = vdims%j_start - offy_spt
    v_sptdims_l%j_end   = vdims%j_end   + offy_spt
    v_sptdims_l%k_start = 1
    v_sptdims_l%k_end   = vdims%k_end
  END IF

    !  Allocate smoothing buffer array
  IF (.NOT. ALLOCATED(delta_v0)) THEN
    ALLOCATE(delta_v0(v_sptdims_l%i_start:v_sptdims_l%i_end,            &
                       v_sptdims_l%j_start:v_sptdims_l%j_end,           &
                       v_sptdims_l%k_start:v_sptdims_l%k_end))
  END IF

   ! Pass in the values to the THETA haloed variable to do the smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        delta_v0(i,j,k) = delta_v(i,j,k)
      END DO
    END DO
  END DO !end loop over k

     !Initialize delta_v for smoothing
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        delta_v(i,j,k)=0.0
      END DO
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( delta_v0, row_length, n_rows,                       &
                    model_levels, offx_spt, offy_spt, fld_type_v,       &
                    swap_field_is_vector )

  !  2D 1-2-1 spatial filter (applied nsmooth times)
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    !  Apply smoother
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        DO jj = -offy_spt, offy_spt
          DO ii = -offx_spt, offx_spt
            delta_v(i,j,k) = delta_v(i,j,k)   +                         &
                      mask_smooth_spt(ii,jj) * delta_v0(ii+i,jj+j,k)
          END DO  ! ii
        END DO  ! jj
      END DO  ! i
    END DO  ! j
  END DO  ! k

  ! Deallocate delta_q0
  IF (ALLOCATED(delta_v0)) DEALLOCATE(delta_v0)

END IF !end nsmooth > 0 loop

! --------------------------------------------------------------------
! Apply conservation if requested
! --------------------------------------------------------------------

IF (l_spt_qcons .OR. l_spt_mse)                                        &
  CALL SPT_conservation(delta_q, delta_theta, q_star,                  &
                        exner_theta_levels, rho_in)

! --------------------------------------------------------------------
! Pass in forced fields to star variables
! --------------------------------------------------------------------

! Copy increments into theta
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      theta_spt(i,j,k)= delta_theta(i,j,k)
    END DO
  END DO
END DO

! Copy increments into q
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      q_spt(i,j,k)= delta_q(i,j,k)
    END DO
  END DO
END DO

! Copy increments into U
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_u_spt(i,j,k)= delta_u(i,j,k)
    END DO
  END DO
END DO

! Copy increments into V
DO k = spt_bot_tap_lev, spt_top_tap_lev
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_v_spt(i,j,k)= delta_v(i,j,k)
    END DO
  END DO
END DO

! --------------------------------------------------------------------
! Copy SPT increments for each variable into std_diag array
! --------------------------------------------------------------------

! Theta SPT increment
IF (spt_diag%l_spt_theta_inc) THEN
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        spt_diag%spt_theta_inc(i,j,k) = delta_theta(i,j,k)
      END DO
    END DO
  END DO
END IF

! Theta SPT increment
IF (spt_diag%l_spt_t_inc) THEN
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        spt_diag%spt_t_inc(i,j,k) = delta_theta(i,j,k) *                     &
                                    exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
END IF

! q SPT increment
IF (spt_diag%l_spt_q_inc) THEN
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        spt_diag%spt_q_inc(i,j,k) = delta_q(i,j,k)
      END DO
    END DO
  END DO
END IF

! u SPT increment
IF (spt_diag%l_spt_u_inc) THEN
  DO k = 1, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        spt_diag%spt_u_inc(i,j,k) = delta_u(i,j,k)
      END DO
    END DO
  END DO
END IF

 ! v SPT increment
IF (spt_diag%l_spt_v_inc) THEN
  DO k = 1, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        spt_diag%spt_v_inc(i,j,k) = delta_v(i,j,k)
      END DO
    END DO
  END DO
END IF

! Add tendencies to physics_tendencies_mod
IF (l_retain_stph_tendencies) THEN

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        dtheta_stph(i,j,k) = delta_theta(i,j,k)
      END DO
    END DO
  END DO
  
  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = udims%j_start,udims%j_end
      DO i = udims%i_start,udims%i_end
        du_stph(i,j,k) = du_stph(i,j,k) + delta_u(i,j,k)
      END DO
    END DO
  END DO

  DO k = spt_bot_tap_lev, spt_top_tap_lev
    DO j = vdims%j_start,vdims%j_end
      DO i = vdims%i_start,vdims%i_end
        dv_stph(i,j,k) = dv_stph(i,j,k) + delta_v(i,j,k)
      END DO
    END DO
  END DO

END IF ! end if over l_retain_stph_tendencies

! --------------------------------------------------------------------
! Output Stash Diagnostics
! --------------------------------------------------------------------

IF (sf(0,35)) THEN

  CALL  diagnostics_spt( row_length, rows,                              &
                         n_rows, at_extremity, spt_diag,                &
                            stashwork35)
    ! ------------------------
    ! Tidy allocatable arrays
    ! ------------------------
  DEALLOCATE(spt_diag%spt_forcing_pattern)
  DEALLOCATE(spt_diag%spt_theta_inc)
  DEALLOCATE(spt_diag%spt_q_inc)
  DEALLOCATE(spt_diag%spt_u_inc)
  DEALLOCATE(spt_diag%spt_v_inc)
  DEALLOCATE(spt_diag%cfl_br_marker)
  DEALLOCATE(spt_diag%spt_t_inc)

END IF

! Deallocate forcing patterns
IF (ALLOCATED(psif)) DEALLOCATE(psif)
IF (ALLOCATED(psif2)) DEALLOCATE(psif2)

! Deallocate mass flux upwards
IF (ALLOCATED(skeb2_up_flux)) DEALLOCATE(skeb2_up_flux)
! Deallocate physics tendencies

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! --------------------------------------------------------------------
END SUBROUTINE SPT_main

END MODULE SPT_main_mod
