! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Interface to downdraughts and evaporation of precipitation

MODULE dd_evap_ctl_mod

IMPLICIT NONE

! Description :
!  See UMDP 27 for full details of the scheme for treatment of evaporation
!  of precipitation and downdraughts. This routine is the interface to the 
!  scheme. 
!  This code assumes the number of tracer levels is equal to the number of
!  model levels so only passes nlev values down to dd_evap_scheme. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UMDP 3 standards

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DD_EVAP_CTL_MOD'

CONTAINS


SUBROUTINE dd_evap_ctl(npnts, nreal, nlev, trlev, ntra, kmax,              &
                       index_real, iccb, kterm, freeze_lev,                &
                       l_tracer,                                           &
                       timestep,                                           &
                       r_theta, r_rho, z_theta, z_rho,                     &
                       dr_across_th, dr_across_rh,                         &
                       exner_layer_centres,exner_layer_boundaries,         &
                       p_layer_centres, p_layer_boundaries,                &
                       rho, rho_dry, rho_theta, rho_dry_theta,             &
                       th, q, u, v, qse,                                   &
                       precip_rain, precip_snow, ud_flx, sigma_up,         &
                       th_ud, q_ud, u_ud, v_ud, tracer_ud, tracer,         &
! inout                     
                       rain, snow, rain_3d, snow_3d,                       &
                       dthbydt, dqbydt, dubydt, dvbydt, dtrabydt,          & 
! out
                       dd_flux, sigma_dd, dt_dd, dq_dd, du_dd, dv_dd)

USE planet_constants_mod, ONLY: g

USE cv_run_mod, ONLY: l_mom_dd

! Subroutines - Will call
!USE dd_evap_scheme_mod, ONLY: dd_evap_scheme 

USE ereport_mod, ONLY: ereport

USE umPrintMgr                       ! Required to write output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
INTEGER, INTENT(IN) :: &
  npnts                & ! number of columns
 ,nreal                & ! points with real convection
 ,nlev                 & ! Number of model levels
 ,trlev                & ! Number of tracer model levels
 ,ntra                 & ! Number of tracers
 ,kmax                   ! Maximum level reached by convection

INTEGER, INTENT(IN) :: &
  index_real(npnts)    & ! points with real convection
 ,iccb(npnts)          & ! Lowest cloud base level
 ,kterm(npnts)         & ! highest cloud top level
 ,freeze_lev(npnts)      ! freezing level i.e. level with a temperature
                         ! just above freezing

LOGICAL, INTENT(IN) :: &
  l_tracer               ! true if tracers to be moved by convection

REAL, INTENT(IN)    :: &
  timestep               ! Convection timestep (s)

REAL, INTENT(IN)    ::      &
  r_theta(npnts,0:nlev)     & ! radius of theta levels (m)
 ,r_rho(npnts,nlev)         & ! radius of rho levels (m)
 ,z_theta(npnts,nlev)       & ! height of theta levels (m)
 ,z_rho(npnts,nlev)         & ! height of rho levels (m)
 ,dr_across_th(npnts,nlev)  & ! thickness of theta levels (m)
 ,dr_across_rh(npnts,nlev)    ! thickness of rho levels (m)

REAL, INTENT(IN)    ::                 &
  exner_layer_centres(npnts,0:nlev)    & ! Exner pressure
 ,exner_layer_boundaries(npnts,0:nlev) & ! Exner at half level above
                                         ! exner_layer_centres
 ,p_layer_centres(npnts,0:nlev)        & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:nlev)       ! Pressure at half level above
                                         ! p_layer_centres (Pa)

REAL, INTENT(IN)    ::      &
  rho(npnts,nlev)           & ! wet density for rho lev (kg/m3)
 ,rho_dry(npnts,nlev)       & ! dry density for rho lev (kg/m3)
 ,rho_theta(npnts,nlev)     & ! wet density for theta lev (kg/m3)
 ,rho_dry_theta(npnts,nlev)   ! dry density for theta lev (kg/m3)

REAL, INTENT(IN)    ::     &
  th(npnts,nlev)           & ! potential temperature (K)
 ,q(npnts,nlev)            & ! water vapour (kg/kg)
 ,u(npnts,nlev)            & ! Model u field (m/s)
 ,v(npnts,nlev)            & ! Model v field (m/s)
 ,qse(npnts,nlev)            ! Saturation water vapour of environment (kg/kg)

REAL ,INTENT(IN)    ::      &
  precip_rain(npnts,nlev)   & ! rain added when descending from layer k to k-1
                              ! (kg/m**2/s)
 ,precip_snow(npnts,nlev)   & ! snow added when descending from layer k to k-1
                              ! (kg/m**2/s)
 ,ud_flx(npnts,nlev)        & ! updraught mass flux (Pa/s)
 ,sigma_up(npnts,nlev)      & ! fractional area of updraught cores 
 ,th_ud(npnts,nlev)         & ! updraught theta (K)
 ,q_ud(npnts,nlev)          & ! updraught q (kg/kg)
 ,u_ud(npnts,nlev)          & ! updraught u (m/s)
 ,v_ud(npnts,nlev)            ! updraught v (m/s)

REAL, INTENT(IN) ::            &
  tracer_ud(npnts,trlev,ntra)  & ! Tracer field value in updraught (kg/kg)
 ,tracer(npnts,trlev,ntra)       ! Model tracer fields (kg/kg)

REAL, INTENT(INOUT) ::     &
  rain(npnts)              & ! Surface convective rainfall (kg/m2/s)
 ,snow(npnts)              & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(npnts,nlev)      & ! Convective rainfall flux (kg/m2/s)
 ,snow_3d(npnts,nlev)      & ! Convective snowfall flux (kg/m2/s)
 ,dthbydt(npnts,nlev)      & ! Increments to potential temperature due to 
                             ! convection (K/s)
 ,dqbydt(npnts,nlev)       & ! Increments to q due to convection (kg/kg/s)
 ,dubydt(npnts,nlev)       & ! Increments to U due to convection (m/s/s)
 ,dvbydt(npnts,nlev)       & ! Increments to V due to convection (m/s/s)
 ,dtrabydt(npnts,nlev,ntra)  ! Increment to tracer due to convection (kg/kg/s)

REAL, INTENT(OUT)    ::    &
  dd_flux(npnts,nlev)      & ! downdraught mass flux (Pa/s)
 ,sigma_dd(npnts,nlev)     & ! fractional downdraught area
 ,dt_dd(npnts,nlev)        & ! tendency for T due to DD & evap (K/s)
 ,dq_dd(npnts,nlev)        & ! tendency for q due to DD & evap (kg/kg/s)
 ,du_dd(npnts,nlev)        & ! tendency for U due to DD  (m/s/s)
 ,dv_dd(npnts,nlev)          ! tendency for V due to DD  (m/s/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::                   &
  i,j,k,l                      ! loop counters

INTEGER ::                   &
  icode                        ! Error code
! compressed arrays

INTEGER    ::                &
  iccb_c(nreal)              & ! Lowest cloud base level
 ,kterm_c(nreal)             & ! highest cloud top level
 ,freeze_lev_c(nreal)          ! freezing level

REAL     ::                  &
  r_theta_c(nreal,0:nlev)    & ! radius of theta levels (m)
 ,r_rho_c(nreal,nlev)        & ! radius of rho levels (m)
 ,z_theta_c(nreal,nlev)      & ! height of theta levels (m)
 ,z_rho_c(nreal,nlev)        & ! height of rho levels (m)
 ,rho_c(nreal,nlev)          & ! wet density for rho lev (kg/m3)
 ,rho_dry_c(nreal,nlev)      & ! dry density for rho lev (kg/m3)
 ,rho_theta_c(nreal,nlev)    & ! wet density for theta lev (kg/m3)
 ,rho_dry_theta_c(nreal,nlev)& ! dry density for theta lev (kg/m3)
 ,dr_across_th_c(nreal,nlev) & ! thickness of theta levels (m)
 ,dr_across_rh_c(nreal,nlev)   ! thickness of rho levels (m)

REAL    ::                               &
  exner_layer_centres_c(nreal,0:nlev)    & ! Exner pressure
 ,exner_layer_boundaries_c(nreal,0:nlev) & ! Exner at half level above
                                           ! exner_layer_centres
 ,p_layer_centres_c(nreal,0:nlev)        & ! Pressure (Pa)
 ,p_layer_boundaries_c(nreal,0:nlev)       ! Pressure at half level above
                                           ! p_layer_centres (Pa)

REAL                         &
  th_c(nreal,nlev)           & ! potential temperature (K)
 ,q_c(nreal,nlev)            & ! water vapour (kg/kg)
 ,u_c(nreal,nlev)            & ! Model u field (m/s)
 ,v_c(nreal,nlev)            & ! Model v field (m/s)
 ,qse_c(nreal,nlev)            ! Saturation water vapour of environment (kg/kg)

REAL                         &
  precip_rain_c(nreal,nlev)  & ! rain added when descending from layer k to k-1
                               ! (kg/m**2/s)
 ,precip_snow_c(nreal,nlev)  & ! snow added when descending from layer k to k-1
                               ! (kg/m**2/s)
 ,ud_flx_c(nreal,nlev)       & ! updraught mass flux (Pa/s)
 ,sigma_up_c(nreal,nlev)     & ! fractional area of updraught cores 
 ,th_ud_c(nreal,nlev)        & ! updraught theta (K)
 ,q_ud_c(nreal,nlev)         & ! updraught q (kg/kg)
 ,u_ud_c(nreal,nlev)         & ! updraught U (m/s)
 ,v_ud_c(nreal,nlev)           ! updraught V (m/s)

REAL                           &
  tracer_c(nreal,nlev,ntra)    & ! Model tracer fields (kg/kg)
 ,tracer_ud_c(nreal,nlev,ntra) & ! tracer fields in updraught (kg/kg)
 ,dtrabydt_c(nreal,nlev,ntra)    ! Increment to tracer due to convection 
                                 ! (kg/kg/s)

REAL                         &
  rain_c(nreal)              & ! Surface convective rainfall (kg/m2/s)
 ,snow_c(nreal)              & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d_c(nreal,nlev)      & ! Convective rainfall flux (kg/m2/s)
 ,snow_3d_c(nreal,nlev)      & ! Convective snowfall flux (kg/m2/s)
 ,dthbydt_c(nreal,nlev)        ! Increment to theta due to convection (K/s)

REAL                         &
  dd_flux_c(nreal,nlev)      & ! downdraught mass flux (kg/m2/s)
 ,sigma_dd_c(nreal,nlev)     & ! fractional downdraught area
 ,dt_dd_c(nreal,nlev)        & ! tendency for T due to DD & evap (K/s)
 ,dq_dd_c(nreal,nlev)        & ! tendency for q due to DD & evap (kg/kg/s)
 ,lp_ud_c(nreal,nlev)        & ! precip condensate in updraught (kg/kg)
 ,du_dd_c(nreal,nlev)        & ! tendency for U due to DD  (m/s/s)
 ,dv_dd_c(nreal,nlev)        & ! tendency for V due to DD  (m/s/s)
 ,t_dd_c(npnts,nlev)         & ! downdraught temperature (K)
 ,q_dd_c(npnts,nlev)         & ! downdraught water vapour (kg/kg)
 ,qrain_dd_c(npnts,nlev)     & ! downdraught rain/snow condensate (kg/kg)
 ,rain3d_dd_c(npnts,nlev)    & ! downdraught rain rate (kg/m2/s)
 ,snow3d_dd_c(npnts,nlev)    & ! downdraught rain rate (kg/m2/s)
 ,rain3d_ud_c(npnts,nlev)    & ! updraught rain rate (kg/m2/s)
 ,snow3d_ud_c(npnts,nlev)      ! updraught rain rate (kg/m2/s)

CHARACTER(LEN=80) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName ='DD_EVAP_CTL'


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! -----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! -----------------------------------------------------------------------
! 1.0 Compress to just those points where convection occurred and 
!     precipitation was generated by the updraught scheme.
!     The downdraught and evaporation scheme will only work on columns 
!     where convection has occurred and generated precipitation in the 
!     updraught.
! -----------------------------------------------------------------------

DO i = 1,nreal
  j=index_real(i)
  iccb_c(i)  = iccb(j)
  kterm_c(i) = kterm(j)
  freeze_lev_c(i) = freeze_lev(j)
  rain_c(i)   = rain(j)
  snow_c(i)   = snow(j)
END DO

DO k = 0,nlev
  DO i = 1,nreal
    j=index_real(i)
    exner_layer_centres_c(i,k)    = exner_layer_centres(j,k)
    exner_layer_boundaries_c(i,k) = exner_layer_boundaries(j,k)
    p_layer_centres_c(i,k)        = p_layer_centres(j,k)
    p_layer_boundaries_c(i,k)     = p_layer_boundaries(j,k)
    r_theta_c(i,k)                = r_theta(j,k)
  END DO
END DO
DO k = 1,nlev
  DO i = 1,nreal
    j=index_real(i)
    r_rho_c(i,k)         = r_rho(j,k)
    z_theta_c(i,k)       = z_theta(j,k)
    z_rho_c(i,k)         = z_rho(j,k)
    rho_c(i,k)           = rho(j,k)
    rho_dry_c(i,k)       = rho_dry(j,k)
    rho_theta_c(i,k)     = rho_theta(j,k)
    rho_dry_theta_c(i,k) = rho_dry_theta(j,k)
    dr_across_th_c(i,k)  = dr_across_th(j,k)
    dr_across_rh_c(i,k)  = dr_across_rh(j,k)
    th_c(i,k)            = th(j,k)
    q_c(i,k)             = q(j,k)
    qse_c(i,k)           = qse(j,k)
    u_c(i,k)             = u(j,k)
    v_c(i,k)             = v(j,k)
    precip_rain_c(i,k)   = precip_rain(j,k)
    precip_snow_c(i,k)   = precip_snow(j,k)
    ud_flx_c(i,k)        = ud_flx(j,k)
    sigma_up_c(i,k)      = sigma_up(j,k)
    th_ud_c(i,k)         = th_ud(j,k)
    q_ud_c(i,k)          = q_ud(j,k)
    dthbydt_c(i,k)       = dthbydt(j,k)
  END DO
END DO

! Only required if DD CMT calculation otherwise can remain unset
! Currently data unavailable? not stored
IF (l_mom_dd) THEN
  DO k = 1,nlev
    DO i = 1,nreal
      j=index_real(i)
      u_ud_c(i,k) = u_ud(j,k)
      v_ud_c(i,k) = v_ud(j,k)
    END DO
  END DO
END IF


IF (l_tracer) THEN
  DO l = 1,ntra
    DO k = 1,nlev
      DO i = 1,nreal
        j=index_real(i)
        tracer_c(i,k,l)    = tracer(j,k,l) 
        tracer_ud_c(i,k,l) = tracer_ud(j,k,l) 
        dtrabydt_c(i,k,l)  = dtrabydt(j,k,l) 
      END DO
    END DO
  END DO
END IF

! Initialise output arrays
DO k = 1,nlev
  DO i = 1,npnts
    dd_flux(i,k) = 0.0
    sigma_dd(i,k) = 0.0
    dt_dd(i,k) = 0.0
    dq_dd(i,k) = 0.0
    du_dd(i,k) = 0.0
    dv_dd(i,k) = 0.0
  END DO
END DO
! Initialise intermediate arrays
DO k = 1,nlev
  DO i = 1,nreal
    rain_3d_c(i,k)  = 0.0
    snow_3d_c(i,k)  = 0.0
    dt_dd_c(i,k)    = 0.0
    dq_dd_c(i,k)    = 0.0
    dd_flux_c(i,k)  = 0.0
    sigma_dd_c(i,k) = 0.0
    du_dd_c(i,k)    = 0.0
    dv_dd_c(i,k)    = 0.0
  END DO
END DO

! ----------------------------------------------------------------------
! 1.1 Calculate any large-scale fluxes/gradients from environmental profiles.
!     (Gungho may be able to provide these directly.)
!     Note the convective parametrization relies on knowing sub-grid
!     fluxes not mean fluxes.
! ----------------------------------------------------------------------
! Chosen not to do this at present

! ----------------------------------------------------------------------
! 2.0 Call scheme  Downdraughts & evaporation of precipitation outside cloud
! ----------------------------------------------------------------------
WRITE(umMessage,'(a)') 'New downdraught scheme not available yet'
CALL umPrint(umMessage,src='dd_evap_ctl')


!CALL dd_evap_scheme( nreal, nlev, ntra, kmax,                              &
!                     iccb_c, kterm_c, freeze_lev_c,                        &
!                     l_tracer, timestep,                                   &
!                     r_theta_c,r_rho_c,                                    &
!                     z_theta_c, z_rho_c, dr_across_th_c, dr_across_rh_c,   &
!                     exner_layer_centres_c, exner_layer_boundaries_c,      &
!                     p_layer_centres_c, p_layer_boundaries_c,              &
!                     rho_c, rho_dry_c, rho_theta_c, rho_dry_theta_c,       &
!                     th_c, q_c, u_c, v_c, qse_c,                           &
!                     precip_rain_c, precip_snow_c, ud_flx_c, sigma_up_c,   &
!                     th_ud_c, q_ud_c, u_ud_c, v_ud_c, tracer_ud_c,tracer_c,&
!                     dthbydt_c, rain_c, snow_c, rain_3d_c, snow_3d_c,      &
!                     dtrabydt_c,                                           & 
!                     dd_flux_c, sigma_dd_c, dt_dd_c, dq_dd_c,              &
!                     du_dd_c, dv_dd_c,                                     & 
!                     t_dd_c, q_dd_c, lp_ud_c, qrain_dd_c,                  &
!                     rain3d_dd_c, snow3d_dd_c, rain3d_ud_c, snow3d_ud_c )

! Cause an error  
icode = 1
cmessage="l_new_dd = .true. - CONTROL CODE ONLY please do not set to true "

CALL ereport(routinename,icode,cmessage)

! ----------------------------------------------------------------------
! 3.0 Expand output back to full arrays.
! ----------------------------------------------------------------------
DO i = 1,nreal
  j=index_real(i)
  rain(j) = rain_c(i)
  snow(j) = snow_c(i)
END DO

! Only copy back for levels processed by dd_evap_scheme

DO k = 1,kmax
  DO i = 1,nreal
    j=index_real(i)
    ! May be rain/snow from deep/shallow before a mid-level call. Add on
    ! profiles from this call.
    rain_3d(j,k) = rain_3d(j,k) + rain_3d_c(i,k)
    snow_3d(j,k) = snow_3d(j,k) + snow_3d_c(i,k)
    dthbydt(j,k) = dthbydt(j,k) + dt_dd_c(i,k)/exner_layer_centres(j,k)
    dqbydt(j,k)  = dqbydt(j,k)  + dq_dd_c(i,k)
    ! convert to Pa/s from kg/m2/s
    dd_flux(j,k)  = dd_flux_c(i,k) *g
    sigma_dd(j,k) = sigma_dd_c(i,k)
    dt_dd(j,k)    = dt_dd_c(i,k)
    dq_dd(j,k)    = dq_dd_c(i,k)
  END DO
END DO

! Only copy back if CMT being calculated in DD

IF (l_mom_dd) THEN
  DO k = 1,kmax
    DO i = 1,nreal
      j=index_real(i)
      dubydt(j,k) = dubydt(j,k) + du_dd_c(i,k)
      dvbydt(j,k) = dvbydt(j,k) + dv_dd_c(i,k)
      du_dd(j,k)  = du_dd_c(i,k)
      dv_dd(j,k)  = dv_dd_c(i,k)
    END DO
  END DO
END IF

! Tracers - passively transported by scheme

IF (l_tracer) THEN
  DO l = 1,ntra
    DO k = 1,kmax
      DO i = 1,nreal
        j=index_real(i)
        dtrabydt(j,k,l) = dtrabydt_c(i,k,l) 
      END DO
    END DO
  END DO
END IF

! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! ----------------------------------------------------------------------
RETURN
END SUBROUTINE dd_evap_ctl

END MODULE dd_evap_ctl_mod
