! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE conv_diag_5a_mod

USE UM_ParParams

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CONV_DIAG_5A_MOD'

CONTAINS
!
!   to diagnose convective occurrence and type
!
SUBROUTINE conv_diag_5a(                                          &

! IN values defining field dimensions and subset to be processed :
       row_length, rows                                                 &

! IN values defining vertical grid of model atmosphere :
      , bl_levels                                                       &
      , p, P_theta_lev, exner_rho                                       &
      , rho_only, rho_theta, z_full, z_half                             &

! IN Model switches
      , l_extra_call                                                    &
      , no_cumulus                                                      &

! IN Cloud data :
      , qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :
      , pstar, q, theta, exner_theta_levels, u_p, v_p, u_0_p, v_0_p     &
      , tstar_land, tstar_sea, tstar_sice, z0msea                       &
      , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm           &
      , tstar, land_mask, frac_land, ice_fract, timestep                &
      , w_copy, w_max                                                   &
      , deep_flag, past_precip, past_conv_ht                            &

! IN surface fluxes
      , fb_surf, ustar                                                  &

! SCM Diagnostics (dummy values in full UM)
      , nSCMDpkgs, L_SCMDiags                                           &

! OUT data required elsewhere in UM system :

      , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad           &
      , ntml,ntpar,nlcl                                                 &
      , cumulus, L_shallow,l_congestus, l_congestus2, conv_type, cin    &
      , cape, wstar, wthvs, entrain_coef, qsat_lcl                      &
      , error                                                           &
       )

!-----------------------------------------------------------------------
! Purpose:
!    To diagnose convection occurrence and type
!
!   Called by ATMOS_PHYSICS2.
! Method:
! NOTE: this routine works on the model levels and not the surface level
! so all level loops start at 1 regardless of whether new dynamics or endgame.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v6 programming standards
!
!-----------------------------------------------------------------------

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                        &
  pdims, pdims_s, tdims_s, tdims, wdims, ScmRowLen, ScmRow
USE bl_option_mod, ONLY: one_third
! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

! Trig information
USE trignometric_mod,  ONLY:                                      &
  sin_theta_longitude, cos_theta_longitude

USE cv_run_mod, ONLY:                                                   &
  icvdiag, l_jules_flux

USE cv_derived_constants_mod, ONLY:                                     &
  ls, lsrcp, lcrcp, gamma_dry

USE planet_constants_mod, ONLY:                                         &
  vkman, cp, kappa, r, repsilon, g, c_virtual

USE jules_sea_seaice_mod, ONLY: z0sice
USE water_constants_mod, ONLY: lc, lf

! subroutines
USE conv_surf_flux_mod, ONLY: conv_surf_flux

USE model_domain_mod,   ONLY: model_type, mt_single_column
USE s_scmop_mod,        ONLY: default_streams,                          &
                              t_inst,d_wet,d_point,                     &
                              scmdiag_bl,scmdiag_conv
USE scmoutput_mod,      ONLY: scmoutput

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE missing_data_mod, ONLY: rmdi
USE UM_ParParams

USE conv_diag_comp_5a_mod, ONLY: conv_diag_comp_5a
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
! (a) Defining horizontal grid and subset thereof to be processed.

INTEGER, INTENT(IN) ::  &
  row_length            & ! Local number of points on a row
 ,rows                    ! Local number of rows in a theta field

! (b) Defining vertical grid of model atmosphere.

INTEGER, INTENT(IN) :: &
  bl_levels              ! Max. no. of "boundary" levels allowed.

REAL, INTENT(IN) ::                          &
  p(pdims_s%i_start:pdims_s%i_end,           & ! pressure  on rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,           &
    pdims_s%k_start:pdims_s%k_end)           &
 ,p_theta_lev(tdims%i_start:tdims%i_end,     & ! P on theta lev (Pa)
              tdims%j_start:tdims%j_end,     &
                          1:tdims%k_end)     &
 ,exner_rho(pdims_s%i_start:pdims_s%i_end,   & ! Exner on rho level
            pdims_s%j_start:pdims_s%j_end,   & !
            pdims_s%k_start:pdims_s%k_end)   &
 ,rho_only(row_length,rows,1:tdims%k_end)    & ! density (kg/m3)
 ,rho_theta(row_length,rows,1:tdims%k_end-1) & ! rho th lev (kg/m3)
 ,z_full(row_length,rows,1:tdims%k_end)      & ! height th lev (m)
 ,z_half(row_length,rows,1:tdims%k_end)        ! height rho lev (m)


LOGICAL,INTENT(IN) ::   &
  l_spec_z0             & ! true if roughness length has been specified
 ,l_extra_call            ! true this is an additional call to conv_diag
                          ! within a timestep

LOGICAL,INTENT(IN) ::   &
  no_cumulus(row_length,rows)   ! Points overruled by BL

! (c) Cloud data.

REAL, INTENT(IN) ::                          &
  qcf(tdims%i_start:tdims%i_end,             & ! Cloud ice (kg per kg air)
      tdims%j_start:tdims%j_end,             &
                  1:tdims%k_end)             &
 ,qcl(tdims%i_start:tdims%i_end,             & ! Cloud liquid water (kg/kg air)
      tdims%j_start:tdims%j_end,             &
                  1:tdims%k_end)             &
 ,cloud_fraction(tdims%i_start:tdims%i_end,  & ! Cloud fraction
                 tdims%j_start:tdims%j_end,  &
                             1:tdims%k_end)

! (d) Atmospheric + any other data not covered so far, incl control.

REAL, INTENT(IN) ::                             &
  pstar(row_length, rows)                       & ! Surface pressure (Pa)
 ,q(tdims%i_start:tdims%i_end,                  & ! water vapour (kg/kg)
    tdims%j_start:tdims%j_end,                  &
                1:tdims%k_end)                  &
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,              &
                    1:tdims%k_end)              &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)
REAL, INTENT(IN) ::            &
  u_p(row_length, rows)        & ! U(1) on P-grid.
 ,v_p(row_length, rows)        & ! V(1) on P-grid.
 ,u_0_p(row_length,rows)       & ! W'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,v_0_p(row_length,rows)       & ! S'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,flux_e(row_length,rows)      & ! Specified surface
                                 !    latent heat flux (W/m^2)
 ,flux_h(row_length,rows)      & ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
 ,ustar_in(row_length,rows)    & ! Specified surface friction velocity (m/s)
 ,z0msea(row_length,rows)      & ! Sea roughness length for momentum (m)
 ,z0m_scm(row_length,rows)     & ! Namelist input z0m (if >0)
 ,z0h_scm(row_length,rows)       ! Namelist input z0h (if >0)

REAL, INTENT(IN) :: &
  tstar_land(row_length, rows) & ! Surface T on land
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows)   ! Surface T on sea-ice

REAL, INTENT(INOUT) ::    &
  tstar(row_length,rows)    ! Surface temperature (K).
                            ! NOTE only inout for SCM

LOGICAL,INTENT(IN) ::         &
  land_mask(row_length, rows)   ! T If land, F Elsewhere.

REAL, INTENT(IN) ::           &
  frac_land(pdims_s%i_start:pdims_s%i_end, & ! fraction of gridbox that is
            pdims_s%j_start:pdims_s%j_end) & ! land, on all points
 ,ice_fract(row_length,rows)                 ! fraction of sea that has ice

REAL, INTENT(IN) ::                  &
  timestep                           & ! timestep (seconds).
 ,w_copy(wdims%i_start:wdims%i_end,  & ! vertical velocity W (m/s)
         wdims%j_start:wdims%j_end,  &
                     0:wdims%k_end)  & ! Not exact match to module values
 ,w_max(row_length,rows)               ! Column max vertical velocity (m/s)

REAL, INTENT(IN) ::              &
  deep_flag(row_length,rows)     & ! 0-1.0, 1 if deep last time step
 ,past_precip(row_length,rows)   & ! convective precip rate last step
                                   ! or a decayed value.
 ,past_conv_ht(row_length,rows)    ! convective height (m)

REAL, INTENT(IN) ::                                                      &
   ustar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! GBM surface friction velocity

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::  &
  nSCMDpkgs               ! No of diagnostics packages

LOGICAL,INTENT(IN) ::   &
  L_SCMDiags(nSCMDpkgs)   ! Logicals for diagnostics packages

REAL, INTENT(INOUT) ::           &
  zh(row_length,rows)              ! Height above surface of top
                                   ! of boundary layer (metres).
REAL, INTENT(OUT) ::             &
  zhpar(row_length,rows)         & ! Height of max parcel ascent (m)
 ,dzh(row_length,rows)           & ! Height of inversion top (m)
 ,qcl_inv_top(row_length,rows)   & ! Parcel water content at inversion top
 ,z_lcl(row_length,rows)         & ! Height of lifting condensation
                                   ! level  (m)
 ,z_lcl_uv(row_length,rows)      & ! Height of lifting condensation
                                   ! level on uv grid (m)
 ,delthvu(row_length,rows)       & ! Integral of undilute parcel buoyancy
                                   ! over convective cloud layer
                                   ! (for convection scheme)
 ,ql_ad(row_length,rows)         & ! adiabatic liquid water content at
                                   ! inversion or cloud top (kg/kg)
 ,entrain_coef(row_length,rows)  & ! Entrainment coefficient
 ,qsat_lcl(row_length,rows)      & ! qsat at cloud base (kg/kg)
 ,cape(row_length, rows)         & ! CAPE from parcel ascent (m2/s2)
 ,cin(row_length, rows)            ! CIN from parcel ascent (m2/s2)

INTEGER, INTENT(OUT) ::      &
  ntpar(row_length,rows)     & ! Max levels for parcel ascent
 ,nlcl(row_length,rows)        ! No. of model layers below the
                               ! lifting condensation level.

LOGICAL,INTENT(OUT) ::         &
  cumulus(row_length,rows)     & ! Logical indicator for convection
 ,l_shallow(row_length,rows)   & ! Logical indicator for shallow Cu
 ,l_congestus(row_length,rows) & ! Logical indicator for congestus Cu
 ,l_congestus2(row_length,rows)  ! Logical ind 2 for congestus Cu

! Required for extra call to conv_diag as values from original BL call
! may return zero values based on initial timestep profiles.
REAL, INTENT(INOUT) ::    &
  wstar(row_length,rows)  & ! Convective sub-cloud velocity scale (m/s)
 ,wthvs(row_length,rows)    ! surface flux of wthv (K m/s)

INTEGER, INTENT(INOUT) :: &
  error                   & ! 0 - no error in this routine
 ,ntml(row_length,rows)     ! Number of model levels in the
                            ! turbulently mixed layer.
!
! Redundant arguments
! -------------------
! Convective type array ::
INTEGER, INTENT(IN) :: &
  conv_type(row_length, rows)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'CONV_DIAG_5A'

INTEGER ::   &
  i,j        & ! LOCAL Loop counter (horizontal field index).
 , ii        & ! Local compressed array counter.
 , k         & ! LOCAL Loop counter (vertical level index).
 ,nunstable    ! total number of unstable points

! uncompressed arrays - all points

REAL ::                       &
  fb_surf(row_length, rows)   & ! Change in theta_v from surface
                                ! to layer 1 (note diff from BL)
 ,z_lcl_nlcl(row_length,rows)   ! Height of nlcl model level (m)

! Arrays added for extra conv_diag calls
INTEGER ::                     &
  ntml_copy(row_length,rows)

! compressed arrays store only values for unstable columns

INTEGER ::                   &
  index_i(row_length*rows)   & ! column number of unstable points
 ,index_j(row_length*rows)     ! row number of unstable points

REAL  ::                      &
  tv1_sd( row_length*rows)    & ! Approx to standard dev of level
                                ! 1 virtual temperature (K).
 ,bl_vscale2(row_length*rows) & ! Velocity scale squared for
                                ! boundary layer eddies (m2/s2)
 ,frac_land_c(row_length*rows)  ! fraction of land compressed

REAL :: w_m ! velocity scale
!-----------------------------------------------------------------------
! Required for SCM diagnostics
!-----------------------------------------------------------------------
REAL :: TmpScm2d(ScmRowLen,ScmRow)                ! Single level field
REAL :: TmpScm3d(ScmRowLen,ScmRow,1:tdims%k_end)  ! Full field

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
!
! svl = (Tl+gz/cp)*(1+(1/repsilon-1)qt)
!     ~ (Tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T)
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation
!  e = rp/(repsilon+r)   - accurate
!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
error = 0


!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

IF ( bl_levels <  1 .OR. rows <  1 .OR. tdims%k_end <  1 ) THEN
  error = 1
  GO TO 9999

END IF
!-----------------------------------------------------------------------
! 1.1a copy ntml values where do not want to overwrite
!-----------------------------------------------------------------------
IF (l_extra_call) THEN
  DO j=1, rows
    DO i=1,row_length
      IF (no_cumulus(i,j)) THEN
        ntml_copy(i,j) = ntml(i,j)
      END IF
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 1.2 initialisation of output arrays
!-----------------------------------------------------------------------

DO j=1,rows
  DO i=1,row_length

    cumulus(i,j)     = .FALSE.
    L_shallow(i,j)   = .FALSE.
    L_congestus(i,j) = .FALSE.
    L_congestus2(i,j) = .FALSE.
    ntml(i,j) = 1
    nlcl(i,j) = 1
    ntpar(i,j) = 1
    delthvu(i,j)  = 0.0
    ql_ad(i,j)    = 0.0
    cape(i,j)     = 0.0
    cin(i,j)      = 0.0
    z_lcl(i,j)    = 0.0       ! set to zero
    zhpar(i,j)    = 0.0       ! set to zero
    dzh(i,j)      = rmdi      ! initialise to missing data
    qcl_inv_top(i,j)= 0.0
    z_lcl_nlcl(i,j) = z_half(i,j,nlcl(i,j)+1)

    ! Set LCL for UV grid to level below z_lcl (note that (at vn5.1) UV
    ! levels are half levels (below) relative to P-grid. Consistent with
    ! BL scheme's treatment of staggered vertical grid.)

    z_lcl_uv(i,j)= z_full(i,j,nlcl(i,j))

    entrain_coef(i,j) = -99.0    ! indicates not set

  END DO
END DO

!-----------------------------------------------------------------------
! 1.4 Calculation of surface buoyancy flux and level one standard
!  deviation of virtual temperature.
!  Also modifies tstar for some SCM runs.
!-----------------------------------------------------------------------

IF (l_jules_flux) THEN

  ii = 0

  DO j=1, rows
    DO i=1, row_length
      IF (fb_surf(i,j) > 0.0) THEN
        ii = ii+1
        w_m = ( 0.25 * zh(i,j) * fb_surf(i,j) +                         &
                ustar(i,j) * ustar(i,j) * ustar(i,j) ) ** one_third
        tv1_sd(ii) = theta(i,j,1) *                                     &
                    (1.0+c_virtual*q(i,j,1)-qcl(i,j,1)-qcf(i,j,1)) *    &
                    ( 1.93 * fb_surf(i,j) / ( w_m * g) )
                     ! tv1_sd = 1.93 * wthvbar / w_m
                     ! fb_surf = g * wthvbar / thv
        bl_vscale2(ii) = MAX(0.01, 2.5 * 2.52 * w_m * w_m)
                       ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
      END IF
    END DO
  END DO

ELSE

  CALL conv_surf_flux(                                                   &
       row_length, rows                                                  &
       , l_spec_z0                                                       &
       , land_mask                                                       &
       , pstar, tstar_land, tstar_sea, tstar_sice, zh, frac_land         &
       , ice_fract, u_p, v_p, u_0_p, v_0_p                               &
       , flux_e, flux_h, ustar_in, z0msea, z0m_scm, z0h_scm              &
       , z_full, q, theta, exner_theta_levels                            &
       , tstar, fb_surf, tv1_sd, bl_vscale2 )

END IF
!-----------------------------------------------------------------------
! 2.0 Decide on unstable points  ( fb_surf > 0.0 )
!     Only work on these points for the rest of the calculations.
!-----------------------------------------------------------------------

nunstable = 0           ! total number of unstable points
DO j=1,rows
  DO i=1,row_length
    IF ( fb_surf(i,j)  >   0.0 ) THEN
      nunstable = nunstable + 1
      index_i(nunstable) = i
      index_j(nunstable) = j
    END IF
  END DO
END DO

! land fraction on just unstable points
IF (nunstable > 0) THEN
  DO ii=1,nunstable
    i = index_i(ii)
    j = index_j(ii)
    frac_land_c(ii) = frac_land(i,j)
  END DO

  !-----------------------------------------------------------------------
  ! Work on just unstable points to diagnose whether convective.
  !-----------------------------------------------------------------------
  CALL conv_diag_comp_5a(                                         &
    row_length, rows                                              &
  , bl_levels, nunstable                                          &
  , index_i, index_j                                              &
  , p, P_theta_lev, exner_rho                                     &
  , rho_only, rho_theta, z_full, z_half                           &
  , qcf, qcl, cloud_fraction                                      &
  , pstar, q, theta, exner_theta_levels                           &
  , land_mask, frac_land_c, timestep                              &
  , w_copy ,w_max, tv1_sd, bl_vscale2                             &
  , deep_flag, past_precip, past_conv_ht                          &
  , nSCMDpkgs, L_SCMDiags                                         &
  , ntml, ntpar, nlcl                                             &
  , cumulus, L_shallow, l_congestus, l_congestus2                 &
  , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad         &
  , cin, cape,entrain_coef                                        &
  , qsat_lcl                                                      &
    )

ELSE

  ! Reset zh as done in conv_diag_comp for all points
  DO j=1,rows
    DO i=1,row_length
      zh(i,j) = z_half(i,j,2)
    END DO
  END DO

  !-----------------------------------------------------------------------

END IF        ! unstable points only

!=========================================================================
! Extra calculation required if conv_diag called more than once per
! model physics timestep
!=======================================================================
IF (l_extra_call) THEN

  DO j=1, rows
    DO i=1,row_length
      IF (fb_surf(i,j) > 0.0) THEN
        wstar(i,j) = (zh(i,j)*fb_surf(i,j))**(1.0/3.0)
        wthvs(i,j) = fb_surf(i,j)*theta(i,j,1)                        &
                               *exner_theta_levels(i,j,1)/g
      ELSE
        wstar(i,j) = 0.0
        wthvs(i,j) = 0.0
      END IF
      ! BL overruled first sweep diagnosis of cumulus for a good reason
      ! So prevent subsequent sweeps from diagnosing convection
      IF (no_cumulus(i,j)) THEN
        cumulus(i,j)   = .FALSE.
        L_shallow(i,j) = .FALSE.

        ! Copy back orignal ntml value
        ntml(i,j) = ntml_copy(i,j)

      END IF
    END DO
  END DO

END IF        ! test on l_extra_call




!-----------------------------------------------------------------------
!  SCM  Diagnostics for no unstable points
!  Unstable diagnostics are output from conv_diag_comp but this routine
!  is not called for stable points so null (zero values) are output from
!  this routine for stable points.
!-----------------------------------------------------------------------
IF (model_type == mt_single_column) THEN
  ! Zero output
  TmpScm2d(:,:)   = 0.0
  TmpScm3d(:,:,:) = 0.0

  IF (nunstable == 0) THEN
    IF (l_scmdiags(scmdiag_conv)) THEN

      ! Single level fields compressed

      CALL scmoutput(TmpScm2d,'k_plume',                          &
           'model level for parcel start',' ',                    &
           t_inst,d_point,default_streams,'',routinename)

      CALL scmoutput(TmpScm2d,'qw_plume',                         &
           'initial parcel water','kg/kg',                        &
           t_inst,d_point,default_streams,'',routinename)

      CALL scmoutput(TmpScm2d,'sl_plume',                         &
           'initial parcel energy','J',                           &
           t_inst,d_point,default_streams,'',routinename)

      ! model level output - compressed fields

      CALL scmoutput(TmpScm3d,'buoy_undil',                       &
           'parcel buoyancy -undilute','K',                       &
           t_inst,d_wet,default_streams,'',routinename)

      CALL scmoutput(TmpScm3d,'t_parc_undil',                     &
           'Undilute Parcel temperature','K',                     &
           t_inst,d_wet,default_streams,'',routinename)


      IF (icvdiag >= 2) THEN          ! dilute parcel options only

        CALL scmoutput(TmpScm3d,'buoy_dil',                       &
             'parcel buoyancy - dilute','K',                      &
             t_inst,d_wet,default_streams,'',routinename)

        CALL scmoutput(TmpScm3d,'t_parc_dil',                     &
             'dilute parcel temperature','K',                     &
             t_inst,d_wet,default_streams,'',routinename)

        CALL scmoutput(TmpScm3d,'entrain_frac',                   &
             'entrainment fraction',' ',                          &
             t_inst,d_wet,default_streams,'',routinename)

      END IF ! test on icvdiag > 2 i.e. dilute parcel

    END IF ! scmdiag_conv

    IF (l_scmdiags(scmdiag_bl)) THEN

      CALL scmoutput(TmpScm3d,'thv_par',                          &
           'thetav parcel','K',                                   &
           t_inst,d_wet,default_streams,'',routinename)

      CALL scmoutput(TmpScm3d,'thv_env',                          &
           'Environment thetav','K',                              &
           t_inst,d_wet,default_streams,'',routinename)


    END IF ! scmdiag_bl

  END IF !  nunstable == 0

  !----------------------------------------------------------------
  ! SCM Diagnostics always present
  !----------------------------------------------------------------
  IF (l_scmdiags(scmdiag_conv)) THEN

    ! single fields not compressed

    CALL scmoutput(entrain_coef,'ent_coef',                       &
         'entrainment coefficient',' ',                           &
         t_inst,d_point,default_streams,'',routinename)

    CALL scmoutput(delthvu,'delthvu',                             &
         'CAPE from conv_diag','J',                               &
         t_inst,d_point,default_streams,'',routinename)

    CALL scmoutput(z_lcl,'z_lcl',                                 &
         'LCL height','m',                                        &
         t_inst,d_point,default_streams,'',routinename)

    CALL scmoutput(z_lcl_nlcl,'z_lcl_nlcl',                       &
         'LCL height at nlcl','m',                                &
         t_inst,d_point,default_streams,'',routinename)

    CALL scmoutput(cin, 'cin',                                    &
         'undilute parcel CIN','J/kg',                            &
         t_inst,d_point,default_streams,'',routinename)

    CALL scmoutput(cape,'cape',                                   &
         'undilute parcel CAPE','J/kg ',                          &
         t_inst,d_point,default_streams,'',routinename)

    ! model level output -

    DO k=1, tdims%k_end
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d(i,j,k)=theta(i,j,k) * exner_theta_levels(i,j,k)
        END DO
      END DO
    END DO

    CALL scmoutput(TmpScm3d,'t_cdiag',                            &
         'temperature cond_diag','K',                             &
         t_inst,d_wet,default_streams,'',routinename)

    CALL scmoutput(q, 'q_cdiag',                                  &
         'q cond_diag','kg/kg',                                   &
         t_inst,d_wet,default_streams,'',routinename)

  END IF ! scmdiag_conv

END IF ! model_type


9999  CONTINUE  ! Branch for error exit.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE conv_diag_5a

END MODULE conv_diag_5a_mod
