! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     Subroutine forcing
!
!     Purpose: Called by scm_main (Single Column Model main routine)
!              to apply the appropriate Forcing (code was previously
!              in the main Calling routine scm_main ).
!
!     Code Description:
!     Language - FORTRAN 90
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
!     Documentation: Single Column Model Guide - J. Lean
!=====================================================================
! Options to set initial profiles
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================
!---------------------------------------------------------------------

SUBROUTINE forcing                                                            &

  ! (In)
  ( row_length, rows, nfor, nbl_levs, nsoilt_levs, nsoilm_levs                &
  , ntrop, sec_day, stepcount, daycount, dayno_wint, daysteps, nSCMdpkgs      &
  , ichgf, t, q, qcl, qcf, u, v, w, l_scmdiags, p, exner_theta_levels         &
  , rp, r_theta_levs, ch_tstar_forcing, ch_ustar_forcing                      &
  , ch_flux_h, ch_flux_e, ch_tls,  ch_qls,  ch_uls,  ch_vls,  ch_wls          &
  , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg                               &
  , ad, at, avn, aw, cdbar, ctbar, cvnbar, cwbar, dbar, tbar, vnbar, wbar     &
  , vpbar, cdsd, ctsd, cvnsd, cwsd, dsd, tsd, vnsd, wsd, tdash, ddash         &
  , deltan, px, py                                                            &

  ! (InOut)
  , ilscnt, flux_h_scm, flux_e_scm, ustar_in, tstar, ti_scm, qi_scm, t_inc_scm&
  , q_star_scm, qcl_inc, qcf_inc, u_inc_scm, v_inc_scm, w_inc_scm             &
  , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                          &
  , tls, qls, uls, vls, wls, tr, qr, vnr, vpr, wr                             &

  ! (Out)
  , factor_rhokh, rhokh, dab1, dap1 )

USE scm_utils, ONLY:                                                        &
  old_rlx, old_vertadv, old_nml                                             &
, rlx_none, rlx_init, rlx_bgrd, rlx_inst_init, rlx_inst_bgrd

USE s_main_force, ONLY:                                                     &
  timestep, ui, vi, wi, obs, obs_surf, stats, prindump_obs, prinstat        &
, l_vertadv, rlx_t, rlx_q, rlx_u, rlx_v, rlx_w                              &
, plev_t, plev_q, plev_u, plev_v, plev_w, tau_t, tau_q, tau_u, tau_v, tau_w

USE nlsizes_namelist_mod, ONLY: model_levels

USE ereport_mod, ONLY: ereport
USE s_scmop_mod, ONLY: default_streams                                      &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps                                      &
  , d_sl, d_soilt, d_bl, d_wet, d_all, d_soilm, d_tile, d_vis, d_point      &
  , d_allxtra, d_land, d_cloud, d_mlsnow                                    &
  , SCMDiag_gen, SCMDiag_rad, SCMDiag_bl, SCMDiag_surf, SCMDiag_land        &
  , SCMDiag_sea, SCMDiag_lsp, SCMDiag_conv, SCMDiag_lscld, SCMDiag_pc2      &
  , SCMDiag_forc, SCMDiag_incs, SCMDiag_gwd, SCMDiag_MLSnow
USE scmoutput_mod, ONLY: scmoutput

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------

!==========================
! Arguments with Intent(In)
!==========================
INTEGER, INTENT(IN) ::  &
  row_length            &! Leading dimension of SCM
, rows                  &
, nfor                  &! Number terms for observational forcing
, nbl_levs              &! Number of Boundary layer levels
, nsoilt_levs           &! Number of soil temp levels
, nsoilm_levs           &! Number of soil moisture levels
, ntrop                  ! Max number of levels in the troposphere

INTEGER, INTENT(IN) ::  &
  sec_day               &
, stepcount             &! Timestep counter
, daycount              &! Daynumber (1 represents 1st january)
, dayno_wint            &! Daynumber relative to winter solstice
, daysteps              &! No. of timesteps in 1 day
, nSCMdpkgs             &! No of SCM diagnostics packages
, ichgf                  ! No. of timesteps between change
                         ! in observational forcing

REAL, INTENT(IN) ::                &
  t(row_length,rows,model_levels)  &! Temperature (K)
, q(row_length,rows,model_levels)  &! Specific humidity (kg/kg)
, qcl(row_length,rows,model_levels)&! Cloud water content (kg/kg)
, qcf(row_length,rows,model_levels)&! Cloud ice content (kg/kg)
, u(row_length,rows,model_levels)  &! Zonal wind (m/s)
, v(row_length,rows,model_levels)  &! Meridional wind (m/s)
, w(row_length,rows,0:model_levels)! Vertical velocity (m/s)

LOGICAL, INTENT(IN) ::  &
  l_scmdiags(nscmdpkgs)  ! Logicals for SCM diagnostics packages


REAL, INTENT(IN) ::                         &
  p(row_length,rows,model_levels)           &  ! Pressure coordinates (Pa)
, exner_theta_levels(row_length,rows,model_levels) &
, rp(row_length,rows,model_levels)          &  ! 1/p for rho levels
                                               ! (1/HPa or 1/mb)
, r_theta_levs(row_length,rows,0:model_levels) ! Distance of theta levels
                                             ! from centre of earth (m)


! Rate of change of surface forcing for ...
REAL, INTENT(IN) ::                         &
  ch_tstar_forcing(row_length,rows,nfor-1)  &! ... tstar,  K/s
, ch_ustar_forcing(row_length,rows,nfor-1)  &! ... ustar,  m/s^2
, ch_flux_h(row_length,rows,nfor-1)         &! ... flux_h, W/(m^2*s)
, ch_flux_e(row_length,rows,nfor-1)          ! ... flux_e, W/(m^2*s)

! Rate of change of atmospheric forcing tendencies for ...
REAL, INTENT(IN) ::                            &
  ch_tls(row_length,rows,nfor-1,model_levels)  &! ... Temperature K/(day*s)
, ch_qls(row_length,rows,nfor-1,model_levels)  &! .. Spec. humid (kg/kg)/(day*s)
, ch_uls(row_length,rows,nfor-1,model_levels)  &! ... Zonal  wind m/(day*s^2)
, ch_vls(row_length,rows,nfor-1,model_levels)  &! ... Merid. wind m/(day*s^2)
, ch_wls(row_length,rows,nfor-1,0:model_levels) ! ... Vert.  wind m/(day*s^2)

! Rate of change of background fields for ...
REAL, INTENT(IN) ::                              &
  ch_t_bg(row_length,rows,nfor-1,model_levels)   &! ... Temperature  K/s
, ch_q_bg(row_length,rows,nfor-1,model_levels)   &! ... Spec. humid  (kg/kg)/s
, ch_u_bg(row_length,rows,nfor-1,model_levels)   &! ... Zonal  wind  m/(s^2)
, ch_v_bg(row_length,rows,nfor-1,model_levels)   &! ... Merid. wind  m/(s^2)
, ch_w_bg(row_length,rows,nfor-1,0:model_levels)  ! ... Vert.  wind  m/(s^2)


! Declarations for Statistical forcing
!--------------------------------------
REAL, INTENT(IN) ::                      &
  ad(row_length,rows,model_levels-1)     &! Term a of eqn. 2.22 for
, at(row_length,rows,model_levels-1)     &!   dew pt depression and temperature
, avn(row_length,rows,model_levels-1)    &! Term a of equ. 2.22 for
, aw(row_length,rows,ntrop-1)             !   horiz. and vertical velocities

! Mean of random variable for ...
REAL, INTENT(IN) ::                      &
  cdbar(row_length,rows,model_levels)    &!   ... Dew Point Depression (K)
, ctbar(row_length,rows,model_levels)    &!   ... Temperature (K)
, cvnbar(row_length,rows,model_levels)   &!   ... Velocity, VN (m/s)
, cwbar(row_length,rows,ntrop)            !   ... Vertical Velocity (mb/s)

! Mean at daycount days from winter solstice of ...
REAL, INTENT(IN) ::                      &
  dbar(row_length,rows,model_levels)     &!   ... Dew Point Depression (K)
, tbar(row_length,rows,model_levels)     &!   ... Temperature (K)
, vnbar(row_length, rows,model_levels)   &!   ... Velocity, VN (m/s)
, wbar(row_length,rows,ntrop)            &!   ... Vertical velocity (mb/s)
, vpbar(row_length,rows,model_levels)     !   ... Velocity, VP (m/s)

! Standard Deviation of random variable for ...
REAL, INTENT(IN) ::                      &
  cdsd(row_length,rows,model_levels)     &!   ... Dew Point Depression (K)
, ctsd(row_length,rows,model_levels)     &!   ... Temperature (K)
, cvnsd(row_length,rows,model_levels)    &!   ... Velocity, VN (m/s)
, cwsd(row_length,rows,ntrop)             !   ... Vertical velocity (mb/s)

! Standard Deviation at daycount days from winter solstice of ...
REAL, INTENT(IN) ::                      &
  dsd(row_length,rows,model_levels)      &!   ... Dew Point Depression (K)
, tsd(row_length,rows,model_levels)      &!   ... Temperature (K)
, vnsd(row_length,rows,model_levels)     &!   ... Velocity, VN (m/s)
, wsd(row_length,rows,ntrop)              !   ... Vertical Velocity (mb/s)


REAL, INTENT(IN) ::                    &
  tdash(row_length,rows,model_levels)  &! Temp. corrections (K)
, ddash(row_length,rows,model_levels)  &! Dew pt. corrections
, deltan(row_length,rows)              &! Radius of area (m)
, px(row_length,rows,ntrop)            &! Reciprocal log functions for calc. of
, py(row_length,rows,ntrop-1)           !   advection used in eqns 2.12 and 2.13


!==============================
! Arguments with intent (InOut)
!==============================
INTEGER, INTENT(INOUT) ::             &
  ilscnt                               ! Count for observational forcing

REAL, INTENT(INOUT) ::                &
  flux_h_scm(row_length,rows)         &! Surface sensible heat flux (W/m2)
, flux_e_scm(row_length,rows)         &! Surface latent heat flux   (W/m2)
, tstar(row_length,rows)              &! Surface temperature        (K)
, ustar_in(row_length,rows)            ! Surface friction velocity  (m/s)

REAL, INTENT(INOUT) ::                &
  ti_scm(row_length,rows,model_levels)&! Initial temperature  (K)
, qi_scm(row_length,rows,model_levels) ! Initial spec. humid. (kg/kg)

! Forcing increments of ...
REAL, INTENT(INOUT) ::                       &
  t_inc_scm(row_length,rows,model_levels)    &! ... Temperature  (K)
, q_star_scm(row_length,rows,model_levels)   &! ... Spec. humid. (kg/kg)
, qcl_inc(row_length,rows,model_levels)      &! ... Cloud liquid (kg/kg)
, qcf_inc(row_length,rows,model_levels)      &! ... Cloud ice    (kg/kg)
, u_inc_scm(row_length,rows,model_levels)    &! ... Zonal  wind  (m/s)
, v_inc_scm(row_length,rows,model_levels)    &! ... Merid. wind  (m/s)
, w_inc_scm(row_length,rows,0:model_levels)   ! ... Vert.  wind  (m/s)

! Obs. background field of ...
REAL, INTENT(INOUT) ::                       &
  t_bg_scm(row_length,rows,model_levels)     &! ... Temperature  (K)
, q_bg_scm(row_length,rows,model_levels)     &! ... Spec. humid. (kg/kg)
, u_bg_scm(row_length,rows,model_levels)     &! ... Zonal  wind  (m/s)
, v_bg_scm(row_length,rows,model_levels)     &! ... Merid. wind  (m/s)
, w_bg_scm(row_length,rows,0:model_levels)    ! ... Vert.  wind  (m/s)


! Forcing tendency due to large-scale horizontal/vertical advection of ...
REAL, INTENT(INOUT) ::                &
  tls(row_length,rows,model_levels)   &! ... Temperature (K)/day
, qls(row_length,rows,model_levels)   &! ... Spec. Humid (kg/kg)/day
, uls(row_length,rows,model_levels)   &! ... Zonal wind  (m/s)/day
, vls(row_length,rows,model_levels)   &! ... Merid wind  (m/s)/day
, wls(row_length,rows,0:model_levels)         ! ... Vert. wind  (m/s)/day


! Statistical forcing
!---------------------
REAL, INTENT(INOUT) ::                &! Randomly sampled ...
  tr(row_length,rows,model_levels,2)  &!   ... Temperature (K)
, qr(row_length,rows,model_levels ,2) &!   ... Specific humidity (kg/kg)
, vnr(row_length,rows,model_levels,2) &!   ... Horizontal velocity (m/s)
, vpr(row_length,rows,model_levels,2) &!   ... Horizontal velocity (m/s)
, wr(row_length,rows,ntrop,2)          !   ... Vertical velocity (mb/s)


!==============================
! Arguments with intent (Out)
!==============================

REAL, INTENT(OUT) ::                  &
  factor_rhokh(row_length,rows)       &
, rhokh(row_length,rows,nbl_levs)     &
, dab1(row_length,rows,44)            &! Observational diagnostics
, dap1(row_length,rows,36,model_levels)! Observational diagnostics


!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------
CHARACTER(LEN=errormessagelength) ::  &
  Cmessage                             ! Error message if ErrorStatus > 0

CHARACTER (LEN=*), PARAMETER ::       &
  RoutineName = 'FORCING'

INTEGER ::                            &
  i, j, k, l

INTEGER ::                            &
  ErrorStatus

REAL ::                               &
  tstpfd                              &
, alpha                               &
, a2out_all(model_levels)             &
, a2out_wet(model_levels)

REAL ::                                &
  qclls(row_length,rows,model_levels)  &
, qcfls(row_length,rows,model_levels)

! The increments to t, u, v, q, qcl, qcf on input
REAL ::                                      &
  t_inc_in(row_length,rows,model_levels)     &
, u_inc_in(row_length,rows,model_levels)     &
, v_inc_in(row_length,rows,model_levels)     &
, w_inc_in(row_length,rows,0:model_levels)   &
, q_inc_in(row_length,rows,model_levels)     &
, qcl_inc_in(row_length,rows,model_levels)   &
, qcf_inc_in(row_length,rows,model_levels)

! The increments to T, u, v, q due to observational forcing
REAL ::                                      &
  t_inc_obs(row_length,rows,model_levels)    &
, u_inc_obs(row_length,rows,model_levels)    &
, v_inc_obs(row_length,rows,model_levels)    &
, w_inc_obs(row_length,rows,0:model_levels)  &
, q_inc_obs(row_length,rows,model_levels)

! The increments to T, q, qcl, qcf due to interactive vertical
! advection
REAL ::                                      &
  t_vertadv(row_length,rows,model_levels)    &
, q_vertadv(row_length,rows,model_levels)    &
, qcl_vertadv(row_length,rows,model_levels)  &
, qcf_vertadv(row_length,rows,model_levels)

REAL ::                               &
  th_p                                &
, th_m

REAL ::                               &
  w_vertadv(row_length,rows,0:model_levels)

REAL ::  &
  factor &! Holding variable for calculations
, dth    &! Change in potential temperature
, dz      ! Layer thickness

! Dr Hook
!=============================================================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------
! Control variable
!---------------------------------------------------------------------
tstpfd = timestep / sec_day
alpha  = timestep / 3600.0

!---------------------------------------------------------------------
! Make copies of the increments as they are on input
!---------------------------------------------------------------------

t_inc_in   (:,:,:) = t_inc_scm  (:,:,:)
u_inc_in   (:,:,:) = u_inc_scm  (:,:,:)
v_inc_in   (:,:,:) = v_inc_scm  (:,:,:)
w_inc_in   (:,:,:) = w_inc_scm  (:,:,:)
q_inc_in   (:,:,:) = q_star_scm (:,:,:)
qcl_inc_in (:,:,:) = qcl_inc    (:,:,:)
qcf_inc_in (:,:,:) = qcf_inc    (:,:,:)

!---------------------------------------------------------------------
! Set instantaneous profiles and budgets to zero for OBS forcing
!---------------------------------------------------------------------

IF (obs) THEN
  IF (prindump_obs) THEN
    dap1(:,:,:,:) = 0.0
    dab1(:,:,:)   = 0.0
  END IF ! prindump_obs
END IF  ! obs

!---------------------------------------------------------------------
! If statistical forcing required:-
!   Set up 2 profiles. 1 for start of day plus 1 for start of
!   following day and linearly interpolate between 2 values for
!   all forcing variables.  Increments to T and Q added and U
!   and V calculated
!---------------------------------------------------------------------


! Initialise qcl and qcf keep
qclls(:,:,:) = 0.0
qcfls(:,:,:) = 0.0

IF (stats) THEN
  ! _ls incs = Atmos_physics1 ; _inc = forced incs
  !  We need to keep hold of the increments calculated from physics1 and
  !  then add them on to the ones calculated from the stats forcing

  tls(:,:,:) = t_inc_scm(:,:,:)
  uls(:,:,:) = u_inc_scm(:,:,:)
  vls(:,:,:) = v_inc_scm(:,:,:)
  wls(:,:,:) = w_inc_scm(:,:,:)

  qls  (:,:,:) = q_star_scm(:,:,:)
  qclls(:,:,:) = qcl_inc(:,:,:)
  qcfls(:,:,:) = qcf_inc(:,:,:)

  ! DEPENDS ON: statstep
  CALL statstep                                                             &
    ! In
    ( row_length, rows, ntrop, deltan, px, py, daysteps                     &
    , stepcount, dayno_wint, tr, vnr, vpr, qr, wr, tbar, tsd, tdash, dbar   &
    , dsd, ddash, vnbar, vpbar, vnsd, wbar, wsd, ctbar, ctsd, at, cdbar     &
    , cdsd, ad, cvnbar, cvnsd, avn, cwbar, cwsd, aw, p, rp, u, v, w, t, q   &
    , prinstat, u_inc_scm, v_inc_scm, w_inc_scm, t_inc_scm, q_star_scm      &
    , daycount, timestep)

  ! qcl and qcf temporarily forced using q/100.0
  qcl_inc(:,:,:) = q_star_scm(:,:,:)/100.0
  qcf_inc(:,:,:) = q_star_scm(:,:,:)/100.0

  !  Need to add on increments calculated in physics1
  t_inc_scm(:,:,:) = t_inc_scm(:,:,:) + tls(:,:,:)
  u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)
  v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)

  q_star_scm (:,:,:) = q_star_scm (:,:,:) + qls   (:,:,:)
  qcl_inc    (:,:,:) = qcl_inc    (:,:,:) + qclls (:,:,:)
  qcf_inc    (:,:,:) = qcf_inc    (:,:,:) + qcfls (:,:,:)


  ! W not used in 4.5, added here, but will need someone to use
  ! stats forcing in anger to modify.
  DO k=1, ntrop
    DO j=1, rows
      DO i=1, row_length
        w_inc_scm(i,j,k) = w_inc_scm(i,j,k) + wls(i,j,k)
      END DO
    END DO
  END DO


ELSE IF (obs) THEN
  ! _ls = forced increments ; _inc = Atmos_Physics1 incs


  !---------------------------------------------------------------------
  !       Select forcing value for time of day
  !---------------------------------------------------------------------
  !
  ! NOTE: ilscnt is used to count the number of observational profiles
  !       as the run proceeds. Among other things, it is used to select
  !       the appropriate change/gradient to the forcing tendencies when
  !       observational intervals greater than the model timestep.
  !       THIS FUNCTIONALITY IS FLAWED! It is recommended that users
  !       have ilscnt set to 0 in their SCM namelists and have obs
  !       forcings commence at the same time they wish to begin their
  !       SCM run.                                    Wong/Kerr-Munslow
  IF (MOD((daycount-1) * INT(sec_day)                             &
    + (stepcount-1) * INT(timestep)                               &
    , ichgf * INT(timestep))  ==  0) THEN
    ilscnt = ilscnt + 1
  END IF

  IF (ilscnt  ==  0) THEN
    ilscnt = 1
  ELSE IF (ilscnt  >=  nfor) THEN
    WRITE(Cmessage,*) 'time exceeds forcing period'
    ErrorStatus = 1

    CALL ereport(RoutineName,ErrorStatus,Cmessage)
  END IF

  IF (obs_surf) THEN
    ! Update surface forcings
    ! -----------------------
    flux_h_scm(:,:) = flux_h_scm(:,:) + timestep*ch_flux_h(:,:,ilscnt)
    flux_e_scm(:,:) = flux_e_scm(:,:) + timestep*ch_flux_e(:,:,ilscnt)
    tstar(:,:)      = tstar(:,:)      + timestep*ch_tstar_forcing(:,:,ilscnt)
    ustar_in(:,:)   = ustar_in(:,:)   + timestep*ch_ustar_forcing(:,:,ilscnt)

    rhokh(:,:,1)      = flux_h_scm(:,:)
    factor_rhokh(:,:) = flux_e_scm(:,:)
  END IF

  ! Update current large-scale forcings
  ! -----------------------------------
  tls(:,:,:) = tls(:,:,:) + timestep*ch_tls(:,:,ilscnt,:)
  uls(:,:,:) = uls(:,:,:) + timestep*ch_uls(:,:,ilscnt,:)
  vls(:,:,:) = vls(:,:,:) + timestep*ch_vls(:,:,ilscnt,:)
  wls(:,:,:) = wls(:,:,:) + timestep*ch_wls(:,:,ilscnt,:)
  qls(:,:,:) = qls(:,:,:) + timestep*ch_qls(:,:,ilscnt,:)


  ! Update observed background state (T,q,u,v,w)
  !---------------------------------------------
  t_bg_scm (:,:,:)   = t_bg_scm (:,:,:) + timestep*ch_t_bg(:,:,ilscnt,:)
  q_bg_scm (:,:,:)   = q_bg_scm (:,:,:) + timestep*ch_q_bg(:,:,ilscnt,:)
  u_bg_scm (:,:,:)   = u_bg_scm (:,:,:) + timestep*ch_u_bg(:,:,ilscnt,:)
  v_bg_scm (:,:,:)   = v_bg_scm (:,:,:) + timestep*ch_v_bg(:,:,ilscnt,:)
  w_bg_scm (:,:,:)   = w_bg_scm (:,:,:) + timestep*ch_w_bg(:,:,ilscnt,:)


  ! Add increment due to large-scale forcings
  !------------------------------------------
  t_inc_scm  (:,:,:) = t_inc_scm(:,:,:)  + tstpfd*tls(:,:,:)
  q_star_scm (:,:,:) = q_star_scm(:,:,:) + tstpfd*qls(:,:,:)

  IF (.NOT. old_vertadv) THEN
    w_inc_scm(:,:,:) = w_inc_scm(:,:,:)  + tstpfd*wls(:,:,:)
  END IF

  ! DEPENDS ON: relax_forcing
  CALL relax_forcing                                                        &
    ! In
    ( row_length, rows, model_levels, t, ti_scm, t_bg_scm                   &
    , p, plev_t, tau_t, rlx_t                                               &

    ! InOut
    , t_inc_scm )

  ! DEPENDS ON: relax_forcing
  CALL relax_forcing                                                        &
    ! In
    ( row_length, rows, model_levels, q, qi_scm, q_bg_scm                   &
    , p, plev_q, tau_q, rlx_q                                               &

    ! InOut
    , q_star_scm )

  IF (old_rlx) THEN
    !------------------------------------------------------------------------
    ! OLD BEHAVIOUR, NOT RECOMMENDED
    ! (would like to get rid of this when possible)
    !------------------------------------------------------------------------
    SELECT CASE (rlx_u)

    CASE (rlx_none)
      u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd

    CASE (rlx_init, rlx_inst_init)
      ! Old EUROCS code
      ! This relaxs back to initial conditions.  This is a mess because the
      ! original eurocs logical to relax uv forcings was not affected by the
      ! old l_windrlx logical. The old eurocs relaxation was also done
      ! BEFORE the winds were updated, while the original non-eurocs
      ! wind relaxtion was done AFTER the winds were updated.

      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                    &
        ! In
        ( row_length, rows, model_levels, u, ui, u_bg_scm, p, plev_u        &
        , tau_u, rlx_u                                                      &

        ! InOut
        , u_inc_scm )

      u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd

    CASE (rlx_bgrd, rlx_inst_bgrd)
      ! This relaxs back to obs increments.
      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                    &
        ! In
        ( row_length, rows, model_levels, u, ui, u_bg_scm, p, plev_u        &
        , tau_u, rlx_u                                                      &

        ! InOut
        , u_inc_scm )
    END SELECT


    SELECT CASE (rlx_v)

    CASE (rlx_none)
      v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

    CASE (rlx_init, rlx_inst_init)
      ! Old EUROCS code
      ! This relaxs back to initial conditions.  This is a mess because the
      ! original eurocs logical to relax uv forcings was not affected by the
      ! old l_windrlx logical. The old eurocs relaxation was also done
      ! BEFORE the winds were updated, while the original non-eurocs
      ! wind relaxtion was done AFTER the winds were updated.

      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                    &
        ! In
        ( row_length, rows, model_levels, v, vi, v_bg_scm, p, plev_v        &
        , tau_v, rlx_v                                                      &

        ! InOut
        , v_inc_scm )

      v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

    CASE (rlx_bgrd, rlx_inst_bgrd)
      ! This relaxs back to obs increments
      ! DEPENDS ON: relax_forcing
      CALL relax_forcing                                                    &
        ! In
        ( row_length, rows, model_levels, v, vi, v_bg_scm, p, plev_v        &
        , tau_v, rlx_v                                                      &

        ! InOut
        , v_inc_scm )

    END SELECT


  ELSE ! Use new revised wind relaxation code

    ! Update
    u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uls(:,:,:)*tstpfd
    v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vls(:,:,:)*tstpfd

    ! DEPENDS ON: relax_forcing
    CALL relax_forcing                                                      &
      ! In
      ( row_length, rows, model_levels, u, ui, u_bg_scm, p, plev_u          &
      , tau_u, rlx_u                                                        &

      ! InOut
      , u_inc_scm )

    ! DEPENDS ON: relax_forcing
    CALL relax_forcing                                                      &
      ! In
      ( row_length, rows, model_levels, v, vi, v_bg_scm, p, plev_v          &
      , tau_v, rlx_v                                                        &

      ! InOut
      , v_inc_scm )

  END IF ! (old_rlx)


  ! DEPENDS ON: relax_forcing
  CALL relax_forcing                                                        &
    ! In
    ( row_length, rows, model_levels, w(:,:,1:model_levels)                 &
    , wi(:,:,1:model_levels), w_bg_scm(:,:,1:model_levels), p, plev_w       &
    , tau_w, rlx_w                                                          &

    ! InOut
    , w_inc_scm(:,:,1:model_levels) )

  !-----------------------------------------------------------------------
  ! SCM Forcing OR Increments Diagnostics Package
  !-----------------------------------------------------------------------

  IF (l_SCMdiags(SCMdiag_forc)                                              &
       .OR. l_SCMdiags(SCMdiag_incs)) THEN

    ! These are large-scale forcings and INCLUDE any relaxation effects
    t_inc_obs(:,:,:) = t_inc_scm(:,:,:)  - t_inc_in(:,:,:)
    q_inc_obs(:,:,:) = q_star_scm(:,:,:) - q_inc_in(:,:,:)
    u_inc_obs(:,:,:) = u_inc_scm(:,:,:)  - u_inc_in(:,:,:)
    v_inc_obs(:,:,:) = v_inc_scm(:,:,:)  - v_inc_in(:,:,:)
    w_inc_obs(:,:,:) = w_inc_scm(:,:,:)  - w_inc_in(:,:,:)

  END IF


  !-----------------------------------------------------------------------
  ! Interactive vertical advection
  !-----------------------------------------------------------------------
  ! NOTE: This assumes that the specified forcing does NOT contain
  !       increments resulting from LS vertical advection.
  !       If this is not the case, using this option will result in
  !       double counting.
  !-----------------------------------------------------------------------


  IF (l_vertadv) THEN

    IF (old_vertadv) THEN
      w_vertadv(:,:,:) = wls(:,:,:)
    ELSE
      ! w_inc_scm at this point has allowed for increments due to
      ! LS-tendency and relaxation. The original code calculates
      ! w*d(x)/dz, though uses w = large scale w tendency. Shouldn't
      ! it be using w + w increment from specified LS forcing routine.

      ! Interactive vertical advection should be using the
      ! current wind profile not the specified LS tendency
      w_vertadv(:,:,:) = w(:,:,:) + w_inc_scm(:,:,:)
    END IF

    t_vertadv(:,:,:)   = 0.0
    q_vertadv(:,:,:)   = 0.0
    qcl_vertadv(:,:,:) = 0.0
    qcf_vertadv(:,:,:) = 0.0

    !-----------------------------------------------------------------
    ! Vertical temperature advection
    !-----------------------------------------------------------------
    ! Specified large scale forcing does not include contribution
    ! from vertical advection.  Approximate contributions due to
    ! vertical advection using upstream approximation.
    ! NOTE: tls holds horizontal advection plus other fixed forcing.

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          factor = - timestep*w_vertadv(i,j,k)*exner_theta_levels(i,j,k)

          IF (w_vertadv(i,j,k) < 0.0) THEN

            ! Subsidence
            IF (k < model_levels) THEN
              th_p = t(i,j,k+1) / exner_theta_levels(i,j,k+1)
              th_m = t(i,j,k)   / exner_theta_levels(i,j,k)

              dth  = th_p - th_m
              dz   = r_theta_levs(i,j,k+1) - r_theta_levs(i,j,k)

              t_vertadv(i,j,k)  = factor*dth / dz
            END IF

          ELSE

            ! Ascent
            IF (k > 1) THEN
              th_p = t(i,j,k)   / exner_theta_levels(i,j,k)
              th_m = t(i,j,k-1) / exner_theta_levels(i,j,k-1)

              dth  = th_p - th_m
              dz   = r_theta_levs(i,j,k) - r_theta_levs(i,j,k-1)

              t_vertadv(i,j,k) = factor*dth / dz
            END IF

          END IF ! w_vertadv < 0.0

        END DO
      END DO
    END DO

    !-----------------------------------------------------------------
    ! Vertical moisture advection
    !-----------------------------------------------------------------
    ! Specified large scale forcing does not include contribution
    ! from vertical advection.  Approximate contributions due to
    ! vertical advection using upstream approximation.
    ! NOTE: qls holds horizontal advection plus other fixed forcing.

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          factor = - timestep*w_vertadv(i,j,k)

          IF (w_vertadv(i,j,k) < 0.0) THEN

            ! Subsidence
            IF (k < model_levels) THEN
              dz = r_theta_levs(i,j,k+1) - r_theta_levs(i,j,k)

              q_vertadv(i,j,k)   = factor*(q(i,j,k+1)   - q(i,j,k))   / dz
              qcl_vertadv(i,j,k) = factor*(qcl(i,j,k+1) - qcl(i,j,k)) / dz
              qcf_vertadv(i,j,k) = factor*(qcf(i,j,k+1) - qcf(i,j,k)) / dz
            END IF

          ELSE

            ! Ascent
            IF (k > 1) THEN
              dz = r_theta_levs(i,j,k) - r_theta_levs(i,j,k-1)

              q_vertadv(i,j,k)   = factor*(q(i,j,k)   - q(i,j,k-1))   / dz
              qcl_vertadv(i,j,k) = factor*(qcl(i,j,k) - qcl(i,j,k-1)) / dz
              qcf_vertadv(i,j,k) = factor*(qcf(i,j,k) - qcf(i,j,k-1)) / dz
            END IF

          END IF ! w_vertadv < 0.0

        END DO
      END DO
    END DO

    t_inc_scm(:,:,:)  = t_inc_scm(:,:,:)  + t_vertadv(:,:,:)
    q_star_scm(:,:,:) = q_star_scm(:,:,:) + q_vertadv(:,:,:)
    qcl_inc(:,:,:)    = qcl_inc(:,:,:)    + qcl_vertadv(:,:,:)
    qcf_inc(:,:,:)    = qcf_inc(:,:,:)    + qcf_vertadv(:,:,:)

  END IF ! vert_adv

  IF (prindump_obs) THEN
    dap1(:,:,10,:) = t_inc_scm(:,:,:) / sec_day
    dap1(:,:,20,:) = q_star_scm(:,:,:) * 1000.0 / sec_day
  END IF

END IF                     ! stats or obs

!-----------------------------------------------------------------------
!     SCM Forcing OR Increments Diagnostics Packages
!-----------------------------------------------------------------------
IF (l_SCMdiags(SCMdiag_forc)                                                &
    .OR. l_SCMdiags(SCMdiag_incs)) THEN

  CALL scmoutput(uls, 'ls_u_inc'                                            &
    , 'Large-scale u-wind forcing tendency', '(m/s)/day'                    &
    , t_inst, d_all, default_streams, '', routinename)

  CALL scmoutput(vls, 'ls_v_inc'                                            &
    , 'Large-scale v-wind forcing tendency', '(m/s)/day'                    &
    , t_inst, d_all, default_streams, '', routinename)

  DO k=1, model_levels
    a2out_all(k) = wls(1,1,k)
  END DO

  CALL scmoutput(a2out_all, 'ls_w_inc'                                      &
    , 'Large-scale w-wind forcing tendency', '(m/s)/day'                    &
    , t_inst, d_all, default_streams, '', routinename)

  IF (.NOT. old_nml) THEN

    CALL scmoutput(u_bg_scm, 'bg_u'                                         &
      , 'Background u-wind', 'm/s'                                          &
      , t_inst, d_all, default_streams, '', routinename)

    CALL scmoutput(v_bg_scm, 'bg_v'                                         &
      , 'Background v-wind', 'm/s'                                          &
      , t_inst, d_all, default_streams, '', routinename)

    DO k=1, model_levels
      a2out_all(k) = w_bg_scm(1,1,k)
    END DO

    CALL scmoutput(a2out_all, 'bg_w'                                        &
      , 'Background w-wind', 'm/s'                                          &
      , t_inst, d_all, default_streams, '', routinename)

    CALL scmoutput(q_bg_scm, 'bg_q'                                         &
      , 'Background specific humidity', 'kg/kg'                             &
      , t_inst, d_all, default_streams, '', routinename)

    CALL scmoutput(t_bg_scm, 'bg_t'                                         &
      , 'Background temperature', 'K'                                       &
      , t_inst, d_all, default_streams, '', routinename)

  END IF ! old_nml

  DO k=1, model_levels
    a2out_all(k) = t_inc_scm(1,1,k) - t_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_all, 'dt_totforc'                                    &
    , 'Total temperature increment from s_forcng','K'                       &
    , t_avg, d_all, default_streams, '', RoutineName)

  DO k=1, model_levels
    a2out_all(k) = u_inc_scm(1,1,k) - u_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_all, 'du_totforc'                                    &
    , 'Total u increment from s_forcng', 'm/s'                              &
    , t_avg, d_all, default_streams, '', RoutineName)

  DO k=1, model_levels
    a2out_all(k) = v_inc_scm(1,1,k) - v_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_all, 'dv_totforc'                                    &
    , 'Total v increment from s_forcng', 'm/s'                              &
    , t_avg, d_all, default_streams, '', RoutineName)

  DO k=1, model_levels
    a2out_wet(k) = q_star_scm(1,1,k) - q_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_wet, 'dq_totforc'                                    &
    , 'Total humidity increment from s_forcng', 'kg/kg'                     &
    , t_avg, d_wet, default_streams, '', RoutineName)

  DO k=1, model_levels
    a2out_wet(k) = qcl_inc(1,1,k) - qcl_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_wet, 'dqcl_totforc'                                  &
    , 'Total QCL increment from s_forcng', 'kg/kg'                          &
    , t_avg, d_wet, default_streams, '', RoutineName)

  DO k=1, model_levels
    a2out_wet(k) = qcf_inc(1,1,k) - qcf_inc_in(1,1,k)
  END DO

  CALL scmoutput(a2out_wet, 'dqcf_totforc'                                  &
    , 'Total QCF increment from s_forcng', 'kg/kg'                          &
    , t_avg, d_wet, default_streams, '', RoutineName)

  !-----------------------------------------------------------------------
  !       SCM Forcing OR Increments Diagnostics Packages
  !       when observational based forcing specified
  !-----------------------------------------------------------------------
  IF (obs) THEN

    CALL scmoutput(t_inc_obs, 'dt_obsforc'                                  &
      , 'Temperature increment from observational forcing', 'K'             &
      , t_avg, d_all, default_streams, '', RoutineName)

    CALL scmoutput(u_inc_obs, 'du_obsforc'                                  &
      , 'U increment from observational forcing', 'm/s'                     &
      , t_avg, d_all, default_streams, '', RoutineName)

    CALL scmoutput(v_inc_obs, 'dv_obsforc'                                  &
      , 'V increment from observational forcing', 'm/s'                     &
      , t_avg, d_all, default_streams, '', RoutineName)

    DO k=1, model_levels
      a2out_all(k) = w_inc_obs(1,1,k)
    END DO

    CALL scmoutput(a2out_all, 'dw_obsforc'                                  &
      , 'W increment from observational forcing', 'm/s'                     &
      , t_avg, d_all, default_streams, '', RoutineName)

    CALL scmoutput(q_inc_obs, 'dq_obsforc'                                  &
      , 'Humidity increment from observational forcing', 'kg/kg'            &
      , t_avg, d_wet, default_streams, '', RoutineName)

  END IF ! obs

  !-----------------------------------------------------------------------
  !       SCM Forcing OR Increments Diagnostics Packages
  !       when vertical advection specified
  !-----------------------------------------------------------------------
  IF (l_vertadv) THEN

    CALL scmoutput(t_vertadv, 'dt_vertadv'                                  &
      , 'Temperature increment from vertical advection', 'K'                &
      , t_avg, d_all, default_streams, '', RoutineName)

    CALL scmoutput(q_vertadv, 'dq_vertadv'                                  &
      , 'Humidity increment from vertical advection', 'kg/kg'               &
      , t_avg, d_wet, default_streams, '', RoutineName)

    CALL scmoutput(qcl_vertadv, 'dqcl_vertadv'                              &
      , 'QCL increment from vertical advection', 'kg/kg'                    &
      , t_avg, d_wet, default_streams, '', RoutineName)

    CALL scmoutput(qcf_vertadv, 'dqcf_vertadv'                              &
      , 'QCF increment from vertical advection', 'kg/kg'                    &
      , t_avg, d_wet, default_streams, '', RoutineName)

  END IF ! l_vertadv

END IF ! l_SCMdiags(SCMdiag_forc).OR l_SCMdiags(SCMdiag_incs)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE forcing
