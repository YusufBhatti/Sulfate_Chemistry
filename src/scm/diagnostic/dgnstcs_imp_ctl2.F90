! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Calculate and output some SCM diagnostics

SUBROUTINE dgnstcs_imp_ctl2                                                   &
  ! In
  ( row_length, rows, rhc_row_length, rhc_rows                                &
  , bl_levels, cloud_levels                                                   &
  , rhcrit, combined_cloud, nclds, cumulus, ntml                              &
  , plsp, conv_rain, conv_snow, ls_rain, ls_snow, R_u, R_v                    &
  , zh, zht, bl_type_1, bl_type_2                                             &
  , bl_type_3, bl_type_4, bl_type_5, bl_type_6, bl_type_7, bl_alltypes        &
  , fqt, ftl, t, q, cf, cfl, cff, T_earliest, q_earliest, qcl_earliest        &
  , qcf_earliest, cf_earliest, cfl_earliest, cff_earliest, p_theta_levels     &
  , lwp, iwp, z0m, z0h_eff_gb, z0m_eff_gb, sea_ice_htf                        &
  , taux, tauy, area_cloud_fraction, p_star                                   &
  , cca_2d, rho1, qcl, qcf, surf_ht_flux_sice, rhcpt                          &
  , surf_ht_flux_gb, ustar_imp, aerosol, nSCMDpkgs, L_SCMDiags                &
  , BL_diag, sf_diag)

USE planet_constants_mod, ONLY: lcrcp
USE atm_fields_bounds_mod, ONLY: vdims_s, vdims, udims_s, pdims_s           &
                               , ScmRowLen,ScmRow
USE conversions_mod, ONLY: recip_pi_over_180
USE murk_inputs_mod, ONLY: l_murk_vis
USE bl_diags_mod, ONLY:                                                     &
    strnewbldiag
USE sf_diags_mod, ONLY: strnewsfdiag
USE bl_option_mod, ONLY:                                                    &
    Fric_heating
USE water_constants_mod, ONLY: lc

USE vertnamelist_mod, ONLY:  first_constant_r_rho_level,                    &
    z_top_of_model, eta_theta, eta_rho, vertlevs

USE level_heights_mod, ONLY: eta_theta_levels

USE nlsizes_namelist_mod, ONLY: model_levels

USE visbty_constants_mod, ONLY: n_vis_thresh, vis_thresh, calc_prob_of_vis

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random

USE beta_precip_mod, ONLY: beta_precip

USE umPrintMgr, ONLY: umPrint, umMessage, newline
USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim

USE s_scmop_mod, ONLY: default_streams                                 &
  , t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div       &
  , t_acc_mult, t_const, only_radsteps                                 &
  , d_sl, d_soilt, d_bl, d_wet, d_all, d_soilm, d_tile, d_vis, d_point &
  , d_allxtra, d_land, d_cloud, d_mlsnow                               &
  , SCMDiag_gen, SCMDiag_rad, SCMDiag_bl, SCMDiag_surf, SCMDiag_land   &
  , SCMDiag_sea, SCMDiag_lsp, SCMDiag_conv, SCMDiag_lscld, SCMDiag_pc2 &
  , SCMDiag_forc, SCMDiag_incs, SCMDiag_gwd, SCMDiag_MLSnow
USE scmoutput_mod, ONLY: scmoutput

USE r2_calc_total_cloud_cover_mod, ONLY: r2_calc_total_cloud_cover

USE qsat_mod, ONLY: qsat, qsat_wat
USE dewpnt_mod, ONLY: dewpnt

IMPLICIT NONE


! Description:
! Much of what is in here is (regrettably) duplication of code
! in diagnostics_lscld and diagnostics_bl. These routines can't
! be called in the SCM because of explicit STASH calls, but we
! still want the diagnostics.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran77 with some bits of Fortran90

! ALL PARAMETERS ARE INTENT In APART FROM T1P5M AND Q1P5M
! WHICH ARE InOut.

INTEGER ::          &
  row_length        &
, rows              &
, rhc_row_length    &
, rhc_rows          &
, bl_levels         &
, cloud_levels      &
, nclds

LOGICAL ::                     &
  cumulus(row_length,rows)     

INTEGER :: ntml(row_length,rows)         ! Height of diagnosed BL top

REAL ::                            &
  rhcrit(model_levels)             &
, plsp(row_length,rows)            &
, conv_rain(row_length,rows)       &
, conv_snow(row_length,rows)       &
, ls_rain(row_length,rows)         &
, ls_snow(row_length,rows)         &
, R_u(udims_s%i_start:udims_s%i_end    &
     ,udims_s%j_start:udims_s%j_end    &
     , model_levels)               &
, R_v(vdims_s%i_start:vdims_s%i_end    &
     ,vdims_s%j_start:vdims_s%j_end    &
     , model_levels)

REAL ::                                               &
  zh(row_length,rows)                                 &
, zht(row_length,rows)

REAL ::                        &
  bl_type_1(row_length,rows)   &! In Indicator set to 1.0 if stable b.l.
                                !    diagnosed, 0.0 otherwise.
, bl_type_2(row_length,rows)   &! In Indicator set to 1.0 if Sc over
                                !    stable surface layer diagnosed,
                                !    0.0 otherwise.
, bl_type_3(row_length,rows)   &! In Indicator set to 1.0 if well
                                !    mixed b.l. diagnosed, 0.0 otherwise.
, bl_type_4(row_length,rows)   &! In Indicator set to 1.0 if decoupled
                                !    Sc layer (not over cumulus) diagnosed,
                                !    0.0 otherwise.
, bl_type_5(row_length,rows)   &! In Indicator set to 1.0 if decoupled
                                !    Sc layer over cumulus diagnosed,
                                !    0.0 otherwise.
, bl_type_6(row_length,rows)   &! In Indicator set to 1.0 if a cumulus
                                !    capped b.l. diagnosed, 0.0 otherwise.
, bl_type_7(row_length,rows)   &! In Indicator set to 1.0 if a
                                !    shear-dominated b.l. diagnosed, 0.0
                                !    otherwise.
, bl_alltypes(ScmRowLen,ScmRow)

REAL ::                          &
  fqt(row_length,rows,bl_levels) &
, ftl(row_length,rows,bl_levels) &
, z0m(row_length,rows)           &! In Roughness length for mom
, z0m_eff_gb(row_length,rows)    &! In Orographic roughness length
, z0h_eff_gb(row_length,rows)    &! In Roughness length for heat

! variables needed for PC2 diagnostics
  , t(row_length, rows, model_levels)                           &
  , q(row_length, rows, model_levels)                           &
  , cf(row_length, rows, model_levels)                          &
  , cfl(row_length, rows, model_levels)                         &
  , cff(row_length, rows, model_levels)                         &

! _earliest arrays contain fields at start of imp_ctl
  , T_earliest(row_length, rows, model_levels)                  &
  , q_earliest(row_length, rows, model_levels)                  &
  , qcl_earliest(row_length, rows, model_levels)                &
  , qcf_earliest(row_length, rows, model_levels)                &

! _earliest values contain the values of temperature, water contents
! and cloud fractions before the boundary layer call.
  , cf_earliest(row_length, rows, model_levels)                 &
  , cfl_earliest(row_length, rows, model_levels)                &
  , cff_earliest(row_length, rows, model_levels)                &
  , p_theta_levels(row_length, rows, model_levels)

REAL ::                              &
  lwp(ScmRowLen,ScmRow)              &! Liquid water path  (kg/m2)
, iwp(ScmRowLen,ScmRow)              &! Ice water path     (kg/m2)
, sea_ice_htf(row_length,rows)       &
, taux(row_length,rows,bl_levels)    &
, tauy(row_length,rows,bl_levels)    &

, area_cloud_fraction(row_length,rows,model_levels) &

, p_star(row_length,rows)            &
, ustar_imp(row_length,rows)         &
, cca_2d(row_length,rows)            &
, rho1(row_length,rows)              &! Air density at level 1/kg/m^3
, qcl(row_length,rows,model_levels)  &
, qcf(row_length,rows,model_levels)  &
, surf_ht_flux_sice(row_length,rows) &
, surf_ht_flux_gb(row_length,rows)   &

, rhcpt(rhc_row_length,rhc_rows,model_levels)     &
, aerosol(pdims_s%i_start:pdims_s%i_end    &
         ,pdims_s%j_start:pdims_s%j_end    &
         , model_levels)                          &
, combined_cloud(row_length,rows,model_levels)

! Logicals for SCM Diagnostics packages
INTEGER ::               &
  nSCMDpkgs               ! No of SCM diagnostics packages

LOGICAL ::               &
  L_SCMDiags(nSCMDpkgs)   ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
TYPE (Strnewbldiag) :: BL_diag
TYPE (Strnewsfdiag) :: sf_diag

! Local variables...
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'DGNSTCS_IMP_CTL2'

INTEGER ::             &
  i                    &
, j                    &
, k                    &
, error_code_ignored   &
, kinvert               ! Vertical index for inverted arrays

! Work array
REAL ::                                       &
  interp_data(row_length*rows*model_levels)   &
, work_2d(row_length,rows)                    &
, work_3d_w(row_length,rows,model_levels)     &
, work_3d(row_length,rows,model_levels)

! Global variables

! Variables used to calculate visibility stuff
REAL ::                                          &
  Beta_LS_Rain(row_length,rows)                  &! Scattering in LS Rain.
, Beta_LS_Snow(row_length,rows)                  &! Scattering in LS Snow.
, Beta_C_Rain(row_length,rows)                   &! Scattering in Conv Rain
, Beta_C_Snow(row_length,rows)                   &! Scattering in Conv Snow
, Vis_Threshold(row_length,rows,1,n_vis_thresh)  &! FOG_FR works for n levels,
                                                  ! we want 1
, PVis(row_length,rows,n_vis_thresh)

REAL ::                                   &
  v_1p5m(row_length,rows)                 &! Visibility overall at 1.5m
, v_no_precip(row_length,rows)            &! Visibility without precip.
, v_ls_precip(row_length,rows)            &! Visibility in LS Precip.
, v_c_precip(row_length,rows)             &! Visibility in Conv Precip.
, v_probs(row_length,rows,1,n_vis_thresh)  ! Vis probs at first level
                                           ! (decimal fraction)

! Total cloud amounts.
! (X,X,1) - max. overlap. (X,X,2) - random overlap.
REAL :: tot_cloud(row_length,rows,2)

! More output diagnostics
REAL ::                             &
  layer_cloud1p5m(row_length,rows)  &! Layer cloud at 1.5m
                                     ! (decimal fraction)
, qcl1p5m(row_length,rows)           ! Cloud liquid water content at
                                     ! 1.5m (kg/kg of air)

! Variables used in the calculation of the low_cloud, med_cloud
! and high_cloud diagnostics.
INTEGER, PARAMETER :: num_cloud_types = 3

REAL ::                           &
  h_split(num_cloud_types+1)      &
, cloud_bound(num_cloud_types+1)

INTEGER ::       &
  kk             &
, level          &
, low_bot_level  &
, low_top_level  &
, med_bot_level  &
, med_top_level  &
, high_bot_level &
, high_top_level


REAL ::                         &
  low_cloud(row_length,rows)    &
, med_cloud(row_length,rows)    &
, high_cloud(row_length,rows)   &
, td1p5m(row_length,rows)       &! 1.5m dewpoint temperature
, rh1p5m(row_length,rows)       &! 1.5m relative humidity
, rhw1p5m(row_length,rows)      &! 1.5m relative humidity over liquid water
, wspd10m(row_length,rows)      &! 10m wind speed
, wdrn10m(row_length,rows)      &! 10m wind direction
, bl_gust(row_length,rows)      &! Shear wind gust
, alpha                         &! 'Angle' of wind
, qst                            ! A saturation mixing ratio

! Parameters
LOGICAL, PARAMETER ::           &
  l_mr_physics = .FALSE.       ! Use mixing ratios

LOGICAL :: pct,  & ! T:Cloud amounts are in %
           avg     ! T:Precip =local*prob

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     SCM Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_lscld)) THEN

  ! To calculate low, medium and high cloud amounts, first need to set up
  ! variables low_bot_level, low_top_level, med_bot_level, etc. In global
  ! and limited-area models this is done during initialisation.

  ! Default settings for h_split as set in routine READLSTA
  ! (which is not called in the SCM)
  !   low   : middle :  high   cloud model levels
  ! (1)->(2):(2)->(3):(3)->(4)

  h_split(1) =   111.0    ! ICAO 1000mb height (m)
  h_split(2) =  1949.0    ! ICAO  800mb height (m)
  h_split(3) =  5574.0    ! ICAO  500mb height (m)
  h_split(4) = 13608.0    ! ICAO  150mb height (m)

  ! Set up cloud type boundaries for low/medium/high cloud as in
  ! routine SETCONA (which is not called in the SCM). See this
  ! routine for more comments on what's going on here.
  DO kk=1, num_cloud_types+1

    level = 1
    DO WHILE ((eta_theta_levels(level)*z_top_of_model <= h_split(kk)) &
               .AND. (level <= model_levels))
      level=level+1
    END DO

    IF (level >  model_levels) THEN
      WRITE(umMessage,'(3(A,I3),A)')                                           &
        '=========================================================='//newline//&
        '| ERROR in ni_imp_ctl:'//                                    newline//&
        '| level ', level, ' > ', model_levels, ' model levels for '         //&
        'cloud type ', kk,                                            newline//&
        '=========================================================='
      CALL umPrint(umMessage,src='dgnstcs_imp_ctl2')
      level = model_levels
    END IF

    cloud_bound(kk) = level

  END DO

  low_bot_level  = cloud_bound(1)
  low_top_level  = cloud_bound(2) - 1
  med_bot_level  = cloud_bound(2)
  med_top_level  = cloud_bound(3) - 1
  high_bot_level = cloud_bound(3)
  high_top_level = cloud_bound(4) - 1

  IF (low_top_level >  cloud_levels) THEN
    WRITE(umMessage,'(2(A,I3),A)')                                             &
      '=========================================================='//  newline//&
      '| ERROR in NI_Imp_Ctl:'//                                      newline//&
      '| No of cloud levels less than top of low level cloud'//       newline//&
      '| cloud_levels [',cloud_levels,'] < low_top_level [',                   &
      low_top_level,']'//                                             newline//&
      '=========================================================='
    CALL umPrint(umMessage,src='dgnstcs_imp_ctl2')
  END IF

  IF (med_top_level >  cloud_levels) THEN

    WRITE(umMessage,'(2(A,I3),A)')                                             &
      '=========================================================='//  newline//&
      '| ERROR in NI_Imp_Ctl:'//                                      newline//&
      '| No of cloud levels less than top of med level cloud'//       newline//&
      '| cloud_levels [',cloud_levels,'] < low_top_level [',                   &
      med_top_level,']'//                                             newline//&
      '=========================================================='
    CALL umPrint(umMessage,src='dgnstcs_imp_ctl2')
  END IF

  IF (high_top_level >  cloud_levels) THEN
    WRITE(umMessage,'(2(A,I3),A)')                                             &
      '=========================================================='//  newline//&
      '| ERROR in NI_Imp_Ctl:'//                                      newline//&
      '| No of cloud levels less than top of low level cloud'//       newline//&
      '| cloud_levels [',cloud_levels,'] < low_top_level [',                   &
      high_top_level,']'//                                            newline//&
      '=========================================================='
    CALL umPrint(umMessage,src='dgnstcs_imp_ctl2')
  END IF

  !---------------------------------------------------------------------
  ! In non-SCM runs the following stuff would be done in diagnostics_lscld

      ! Diagnostic 203: Low Cloud Amount
      ! Initialize cloud amount to lowest level value.
  DO j=1, rows
    DO i=1, row_length
      low_cloud(i,j)=area_cloud_fraction(i,j,low_bot_level)
    END DO
  END DO

  ! Cloud amount is calculated under maximum overlap assumption.
  ! Limit cloud amount to maximum of 1.
  DO k=low_bot_level+1, low_top_level
    DO j=1, rows
      DO i=1, row_length
        low_cloud(i,j) = MAX(low_cloud(i,j),area_cloud_fraction(i,j,k))
        low_cloud(i,j) = MIN(low_cloud(i,j),1.0)
      END DO
    END DO
  END DO

  CALL scmoutput(low_cloud, 'lowcld'                                        &
    , 'Low cloud fraction', 'Fraction'                                      &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Diagnostic 204: Medium Cloud Amount
  ! Initialize cloud amount to lowest level value.
  DO j=1, rows
    DO i=1, row_length
      med_cloud(i,j)=area_cloud_fraction(i,j,med_bot_level)
    END DO
  END DO

  ! Cloud amount is calculated under maximum overlap assumption.
  ! Limit cloud amount to maximum of 1.
  DO k=med_bot_level+1, med_top_level
    DO j=1, rows
      DO i=1, row_length
        med_cloud(i,j) = MAX(med_cloud(i,j),area_cloud_fraction(i,j,k))
        med_cloud(i,j) = MIN(med_cloud(i,j),1.0)
      END DO
    END DO
  END DO

  CALL scmoutput(med_cloud, 'medcld'                                        &
    , 'Medium cloud fraction', 'Fraction'                                   &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Diagnostic 205: High Cloud Amount
  ! Initialize cloud amount to lowest level value.
  DO j=1, rows
    DO i=1, row_length
      high_cloud(i,j)=area_cloud_fraction(i,j,high_bot_level)
    END DO
  END DO

  ! Cloud amount is calculated under maximum overlap assumption.
  ! Limit cloud amount to maximum of 1.
  DO k=high_bot_level+1, high_top_level
    DO j=1, rows
      DO i=1, row_length
        high_cloud(i,j) = MAX(high_cloud(i,j),area_cloud_fraction(i,j,k))
        high_cloud(i,j) = MIN(high_cloud(i,j),1.0)
      END DO
    END DO
  END DO

  CALL scmoutput(high_cloud, 'highcld'                                      &
    , 'High cloud fraction', 'Fraction'                                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Total Cloud Amount RANDOM Overlap
  CALL r2_calc_total_cloud_cover                                            &
    ( row_length*rows, model_levels, nclds, ip_cloud_mix_random             &
    , combined_cloud(1,1,1), tot_cloud(1,1,2), row_length*rows, model_levels )

  ! Total Cloud Amount MAX/RANDOM Overlap
  CALL r2_calc_total_cloud_cover                                            &
    ( row_length*rows, model_levels, nclds, ip_cloud_mix_max                &
    , combined_cloud(1,1,1), tot_cloud(1,1,1), row_length*rows              &
    , model_levels )

  CALL scmoutput(tot_cloud(1,1,2), 'tcarndm'                                &
    , 'Total cloud amount (random overlap)', 'Fraction'                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(tot_cloud(1,1,1), 'tcamxrn'                                &
    , 'Total cloud amount (max. random)', 'Fraction'                        &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Note that stash equivalents 9,18[123] are identical to 3,18[123]
  ! and so done below stash 9,226
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        IF (area_cloud_fraction(i, j, k) <= 0.0) THEN  
          work_3d_w(i, j, k) = 0.0
        ELSE
          work_3d_w(i, j, k) = 1.0
        END IF
      END DO ! i
    END DO ! j
  END DO ! k

  !
  CALL scmoutput( work_3d_w, 'lyrcldfreq'                                   &
    , 'Layer cloud frequency indicator', 'Indicator'                        &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 9,228
  CALL scmoutput(rhcpt, 'rhcpt'                                             &
    , 'Critical relative humidity', '%'                                     &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 9,229
  DO k=1, model_levels
    CALL qsat(work_2d, t(:,:,k), p_theta_levels(:,:,k),row_length,rows)

    DO j=1, rows
      DO i=1, row_length

        ! relative humidity in per cent
        work_3d_w(i, j, k) = q(i, j, k) / work_2d(i, j) *100.0

        ! Supersaturation (>100%) can occur with mixed phase
        ! scheme but negative humidity is removed from the
        ! diagnostic:
        IF ( work_3d_w(i, j, k) < 0.0) THEN
          work_3d_w(i, j, k) = 0.0
        END IF

      END DO ! i
    END DO ! j
  END DO ! k

  !
  CALL scmoutput(work_3d_w, 'rhaftcld'                                      &
    , 'Relative humidity after main cloud', '%'                             &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 9,231
  DO k=1, model_levels

    ! NB: Convention in Sect 70 (Radiation) is to invert levels,
    !     1 at top. Combined_cloud is calculated this way but
    !     re-inverted for STASH.

    kinvert = model_levels+1-k

    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i, j, k) = combined_cloud(i,j,kinvert)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'combca_ls'                                     &
    , 'Combined cloud amount in each layer', 'Fraction'                     &
    , t_avg, d_wet, default_streams, '', RoutineName)


  ! End of stuff that would be done in diagnostics_lscld               !
  !---------------------------------------------------------------------


END IF ! L_SCMDiags(SCMDiag_lscld)


!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_surf)) THEN

  ! Stash 3,255
  CALL SCMoutput(sf_diag%q1p5m, 'qt1p5m'                                    &
    , '1.5m total water kg water/kg air', 'kgH20/kgAIR'                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,254
  CALL SCMoutput(sf_diag%t1p5m, 'tl1p5m'                                    &
    , '1.5m liquid temperature', 'K'                                        &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_surf)

!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)) THEN

  CALL SCMoutput(lwp, 'LWP'                                                 &
    , 'Liquid water path', 'kg/m2'                                          &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL SCMoutput(iwp, 'IWP'                                                 &
    , 'Ice water path', 'kg/m2'                                             &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)


!---------------------------------------------------------------------
! In non-SCM runs the following stuff would be done in diagnostics_bl

  ! This call to ls_cld converts 1.5m TL and QT to T and q,
  ! assuming p_star and rhcrit at model level 1 approximates to
  ! 1.5m. This should be done even if ssfm=.false. It also
  ! calculates layer_cloud1p5m and qcl1p5m.
! DEPENDS ON: ls_cld
CALL ls_cld                                                                 &
  ( p_star, rhcpt, 1, bl_levels, rhc_row_length, rhc_rows                   &
  , ntml, cumulus                                                           &
  , l_mr_physics, sf_diag%t1p5m, layer_cloud1p5m, sf_diag%q1p5m, qcf,       &
    qcl1p5m                                                                 &
  , interp_data(2*row_length*rows+1), interp_data(3*row_length*rows+1)      &
  , error_code_ignored)

!-----------------------------------------------------------------------
!     SCM Large Scale Cloud OR Surface Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)) THEN

  CALL scmoutput(layer_cloud1p5m, 'lca1p5m'                                 &
    , '1.5m layer cloud amount', 'Fraction'                                 &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(qcl1p5m, 'qcl1p5m'                                         &
     , '1.5m cloud water', 'kg/kg'                                          &
     , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)


!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_surf)) THEN

  ! 1.5m Fog Fraction and Mist Fraction
  DO i=1, row_length
    DO j=1, rows
      DO k=1, n_vis_thresh
        vis_threshold(i,j,1,k) = vis_thresh(k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: fog_fr
  CALL fog_fr                                                               &
    ( p_star, rhcrit, 1, row_length*rows, sf_diag%t1p5m, aerosol,L_murk_vis,&
      sf_diag%q1p5m                                                         &
    , qcl, qcf, vis_threshold, PVis, n_vis_thresh )

  pct = .FALSE.
  avg = .TRUE.
   ! Calculate scattering coefficients due to precipitation
  CALL beta_precip                                                          &
    ( ls_rain, ls_snow, conv_rain,conv_snow,qcf(1,1,1), rho1, sf_diag%T1p5m,&
      p_star                                                                &
    , plsp, cca_2d, pct, avg, row_length*rows, row_length*rows, 1           &
    , beta_ls_rain, beta_ls_snow, beta_c_rain, beta_c_snow                  &
    , error_code_ignored )

   ! Calculate screen level probability of visibility
   ! less than thresholds
  ! DEPENDS ON: calc_vis_prob
  CALL calc_vis_prob                                                        &
    ( p_star, rhcrit, 1, row_length*rows, row_length*rows, sf_diag%t1p5m,   &
      aerosol                                                               &
    , l_murk_vis, sf_diag%q1p5m, qcl, qcf, vis_thresh, n_vis_thresh, plsp,  &
      cca_2d                                                                &
    , pct, beta_ls_rain, beta_ls_snow, beta_c_rain, beta_c_snow             &
    , v_probs, error_code_ignored )

   ! Visibility at 1.5 m including precipitation
  ! DEPENDS ON: visbty
  CALL visbty                                                               &
    ( p_star, sf_diag%T1p5m, sf_diag%q1p5m, Qcl, Qcf, Aerosol,              &
      calc_prob_of_vis, RHcrit                                              &
    , L_murk_vis, row_length*rows, v_no_precip )

  ! DEPENDS ON: vis_precip
  CALL vis_precip                                                           &
    ( v_no_precip, plsp, cca_2d, pct, Beta_ls_Rain, Beta_ls_Snow            &
    , Beta_C_Rain, Beta_C_Snow, row_length*rows, row_length*rows, 1         &
    , v_1p5m, v_ls_precip, v_C_Precip, error_code_ignored)

  CALL scmoutput(v_probs(1,1,1,1), 'pfog1p5m'                               &
    , 'Probability of fog at 1.5m', '-'                                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(v_probs(1,1,1,2), 'pmist1p5m'                              &
    , 'Probability of mist at 1.5m', '-'                                    &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(v_probs, 'pvisthresh'                                      &
    , 'Probability of vis<=threshold', '-'                                  &
    , t_avg, d_vis, default_streams, '', RoutineName)

  CALL scmoutput(v_no_precip, 'visnop1p5m'                                  &
    , '1.5m visibility outside precip', 'm'                                 &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(v_ls_precip, 'vislsp1p5m'                                  &
    , '1.5m visibility in LS precip', 'm'                                   &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(v_c_precip, 'viscp1pm5'                                    &
    , '1.5m visibility in conv precip', 'm'                                 &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(v_1p5m, 'vis1p5m'                                          &
    , '1.5m visibility', 'm'                                                &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Calculate relative humidities at screen level
  DO j=1, rows
    DO i=1, row_length
      CALL qsat(qst, sf_diag%t1p5m(i,j), p_star(i,j))
      rh1p5m(i,j) = sf_diag%q1p5m(i,j)/qst*100.0
      CALL qsat_wat(qst, sf_diag%t1p5m(i,j), p_star(i,j))
      rhw1p5m(i,j) = sf_diag%q1p5m(i,j)/qst*100.0
    END DO
  END DO

  CALL scmoutput(rh1p5m, 'rh1p5m'                                            &
    , 'Relative humidity at 1.5m', '%'                                       &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(rhw1p5m, 'rhw1p5m'                                          &
    , 'Relative humidity wrt H2O at 1.5m', 'W/m2'                            &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Calculate dewpoint

  CALL dewpnt (sf_diag%q1p5m, p_star, sf_diag%t1p5m, row_length*rows, td1p5m)

  CALL scmoutput(td1p5m, 'td1p5m'                                            &
    , '1.5m dewpoint temperature', 'W/m2'                                    &
    , t_avg, d_sl, default_streams, '', RoutineName)

    ! Calculate the magnitude and direction (by meteorological
    ! convention) of the wind at 10m from the 10m u and v
    ! components
  DO j=1, rows
    DO i=1, row_length
      wspd10m(i,j) = SQRT(sf_diag%u10m(i,j)**2+sf_diag%v10m(i,j)**2)
      IF (wspd10m(i,j) == 0.0) THEN
        wdrn10m(i,j) = 0.0
      ELSE
        alpha = ASIN(sf_diag%v10m(i,j)/wspd10m(i,j)) * recip_pi_over_180
        IF (sf_diag%u10m(i,j) >  0.0) THEN
          wdrn10m(i,j) = 270 - alpha
        ELSE
          wdrn10m(i,j) = 90 + alpha
        END IF
      END IF
    END DO
  END DO

  CALL scmoutput(wspd10m, 'wspd10m'                                        &
    , '10m wind speed', 'm/s'                                              &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(wdrn10m, 'wdrn10m'                                        &
    ,'10m wind direction', 'degs'                                          &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Boundary layer gust diagnostic
  DO j=1, rows
    DO i=1, row_length
      bl_gust(i,j) = wspd10m(i,j)+2.0*2.5                                &
           *SQRT(sqrt(taux(i,j,1)**2+tauy(i,j,1)**2))
    END DO
  END DO

  CALL scmoutput(bl_gust, 'gust10m'                                        &
    , '10m Gust', 'm/s'                                                    &
    , t_avg, d_sl, default_streams, '', RoutineName)


END IF ! L_SCMDiags(SCMDiag_surf)


!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_sea)) THEN

  CALL scmoutput(sf_diag%sice_mlt_htf, 'sice_mlt_htf'                       &
    , 'Heat flux due to melting sea ice', 'W/m2'                            &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sea_ice_htf, 'sea_ice_htf'                                 &
    , 'Heat flux through sea ice', 'W/m2'                                   &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(surf_ht_flux_sice, 'surf_ht_flux_si'                       &
    , 'Net downward heat flux into sea frctn', 'W/m2'                       &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_sea)

!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_surf)) THEN

  CALL scmoutput(sf_diag%latent_heat, 'lat_ht'                              &
    , 'Surface latent heat flux', 'W/m2'                                    &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(surf_ht_flux_gb, 'surf_ht_flux'                            &
    , 'Net downward heat flux at surface', 'W/m2'                           &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(z0m, 'z0m'                                                 &
    , 'Roughness length for momentum', 'm'                                  &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(z0h_eff_gb, 'z0h'                                          &
    , 'Effective roughness length for heat', 'm'                            &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(z0m_eff_gb, 'z0m_eff'                                      &
    , 'Effective roughness length for momentum', 'm'                        &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%q1p5m, 'q1p5m'                                     &
    , '1.5m specific humidity', 'kg/kg'                                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%t1p5m, 't1p5m'                                     &
    , '1.5m temperature', 'K'                                               &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%t1p5m, 't1p5m_max'                                 &
    , 'Max 1.5m temperature', 'K'                                           &
    , t_max, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%t1p5m, 't1p5m_min'                                 &
    , 'Min 1.5m temperature', 'K'                                           &
    , t_min, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%u10m, 'u10m'                                       &
    , 'Zonal 10m wind', 'm/s'                                               &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(sf_diag%v10m, 'v10m'                                       &
    , 'Meridional 10m wind', 'm/s'                                          &
    , t_avg, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_bl)) THEN

  ! Stash 3, 025
  CALL scmoutput(zh, 'zh'                                                   &
    , 'Boundary layer depth after B.layer', 'm'                             &
    , t_avg, d_sl, default_streams, '', RoutineName)

  !       Since model level fields are on rho_levels and the surface
  !       is a theta_level, these are output as separate diagnostics

  ! Stash 3, 217
  DO j=1, rows
    DO i=1, row_length
      work_2d(i, j) = ftl(i, j, 1)
    END DO ! i
  END DO ! j

  CALL scmoutput(work_2d, 'ftl_surf'                                        &
    , 'Surface sensible heat flux from B.layer', 'W/m2'                     &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3, 223
  DO j=1, rows
    DO i=1, row_length
      work_2d(i,j) = fqt(i,j,1)
    END DO ! i
  END DO ! j

  CALL scmoutput(work_2d, 'fqt_surf'                                        &
    , 'Surface sensible moisture flux from B.layer', 'kg/m2/s'              &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(ustar_imp, 'ustar_imp'                                     &
    , 'Implicit friction velocity', 'm/s'                                   &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,222
  CALL scmoutput(fqt, 'fqt_bl'                                              &
    , 'Sensible moisture flux', 'kg/m2/s'                                   &
    , t_avg, d_bl, default_streams, '', RoutineName)

  ! Stash 3,473
  CALL scmoutput(BL_diag%tke, 'TKE'                                         &
    , 'BL diagnostic of Turbulent Kinetic Energy', 'm2/s'                   &
    , t_avg, d_bl, default_streams, '', RoutineName)

  IF (Fric_heating /= 0) THEN
    ! Stash 3,188
    CALL scmoutput(BL_diag%dTfric, 'dt_fric'                                  &
      , 'T increment from turbulence dissipation', 'K'                        &
      , t_avg, d_bl, default_streams, '', RoutineName)
  END IF

  ! Stash 3,304
  CALL scmoutput(zht, 'zht'                                                 &
    , 'Turbulent mixing height after B.layer', 'm'                          &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,305
  CALL scmoutput(bl_type_1, 'bl_type_1'                                     &
    , 'Boundary layer type: stable', 'Indicator'                            &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,306
  CALL scmoutput(bl_type_2, 'bl_type_2'                                     &
    , 'Boundary layer type: Sc over stable', 'Indicator'                    &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,307
  CALL scmoutput(bl_type_3, 'bl_type_3'                                     &
    , 'Boundary layer type: well mixed', 'Indicator'                        &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,308
  CALL scmoutput(bl_type_4, 'bl_type_4'                                     &
    , 'Boundary layer type: decoup Sc not over Cu', 'Indicator'             &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,309
  CALL scmoutput(bl_type_5, 'bl_type_5'                                     &
    , 'Boundary layer type: decoup Sc over Cu', 'Indicator'                 &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,310
  CALL scmoutput(bl_type_6, 'bl_type_6'                                     &
    , 'Boundary layer type: cumulus capped', 'Indicator'                    &
    , t_avg, d_sl, default_streams, '', RoutineName)

  ! Stash 3,340
  CALL scmoutput(bl_type_7, 'bl_type_7'                                     &
    , 'Boundary layer type: shear driven', 'Indicator'                      &
    , t_avg, d_sl, default_streams, '', RoutineName)

  CALL scmoutput(bl_alltypes, 'bl_alltypes'                                 &
    , 'Boundary layer types', ''                                            &
    , t_inst, d_sl, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_bl)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Increments Diagnostics Package
!-----------------------------------------------------------------------
IF (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)) THEN

  ! Stash 3,181 and 9,181
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d(i,j,k) = t(i,j,k) - T_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d, 'dt_pc2blls'                                      &
    , 'Temperature inc PC2+bdy layer+ls cld', 'K'                           &
    ,  t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 3,182 and 9,182
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = q(i,j,k) - q_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dq_pc2blls'                                    &
    , 'Specific humidity inc PC2+bdy layer+ls cld', 'kg/kg'                 &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 3,183 and 9,183
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = qcl(i,j,k) - qcl_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dqcl_pc2blls'                                  &
    , 'QCL increment PC2+bdy layer+ls cld', 'kg/kg'                         &
    , t_avg, d_wet, default_streams,  '',  RoutineName)

  ! Stash 3,184
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = qcf(i,j,k) - qcf_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dqcf_pc2blls'                                  &
    , 'QCF increment PC2+bdy layer+ls cld', 'kg/kg'                         &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 3,185
  CALL scmoutput(bl_diag%u_incr, 'du_bl'                                    &
    , 'U wind increment bdy layer', 'm/s'                                   &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 3,186
  CALL scmoutput(bl_diag%v_incr, 'dv_bl'                                    &
    , 'V wind increment bdy layer', 'm/s'                                   &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 3,192
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = cf(i,j,k) - cf_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dbcf_bl'                                       &
    , 'Bulk cloud fraction increment bdy layer', 'Fraction'                 &
    ,  t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 3,193
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = cfl(i,j,k) - cfl_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dcfl_bl'                                       &
    , 'Liquid cloud fraction increment bdy layer', 'Fraction'               &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 3,194
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) = cff(i,j,k) - cff_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'dcff_bl'                                       &
    , 'Frozen cloud fraction increment bdy layer', 'Fraction'               &
    , t_avg, d_wet, default_streams, '', RoutineName)

  ! Stash 3,189
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d(i,j,k) = t(i,j,k) -   T_earliest(i,j,k)                     &
             - lcrcp * ( qcl(i,j,k) - qcl_earliest(i,j,k) )
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d, 'dlwt_bl'                                         &
    , 'Liquid water temp increment bdy layer', 'K'                          &
    , t_avg, d_all, default_streams, '', RoutineName)

  ! Stash 3,190
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        work_3d_w(i,j,k) =   q(i,j,k) -   q_earliest(i,j,k)                 &
                         + qcl(i,j,k) - qcl_earliest(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(work_3d_w, 'tqinc_bl'                                      &
    , 'Total (liquid) water increment bdy layer', 'kg/kg'                   &
    , t_avg, d_all, default_streams, '', RoutineName)

END IF ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)


! End of stuff that would be done in diagnostics_bl
!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE dgnstcs_imp_ctl2

