! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE parcel_ascent_mod
IMPLICIT NONE
PRIVATE
! Calculates parcel ascent
PUBLIC parcel_ascent
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PARCEL_ASCENT_MOD'

CONTAINS
!
! Subroutine Interface:
SUBROUTINE parcel_ascent (npnts, nSCMDpkgs,                              &
                          nlcl, k_plume,                                 &
                          l_dilute, L_SCMDiags,                          &
                          sl_plume, qw_plume,                            &
                          t, Tl, qcl, qcf, qw, t_dens_env,               &
                          p_theta_lev, exner_theta_levels,               &
                          z_theta, entrain_fraction,                     &
                          T_parc, t_parc_dil, ql_parc,                   &
                          buoyancy, buoyancy_dil,                        &
                          env_svl, par_svl, par_svl_dil,                 &
                          denv_bydz, dpar_bydz, dqsatdz )

USE cv_run_mod, ONLY:                                                    &
  plume_water_load, dil_plume_water_load, qlmin, fac_qsat

USE cv_param_mod, ONLY:                                                  &
  qlcrit

USE cv_derived_constants_mod, ONLY:                                      &
  ls, lsrcp, lcrcp, gamma_dry

USE planet_constants_mod, ONLY:                                          &
  r, repsilon, c_virtual, g

USE water_constants_mod, ONLY: lc, lf, tm

USE model_domain_mod,      ONLY: model_type, mt_single_column
USE atm_fields_bounds_mod, ONLY: ScmRowLen, ScmRow
USE s_scmop_mod,           ONLY: default_streams,                        &
                                 t_inst,d_wet,scmdiag_conv
USE scmoutput_mod,         ONLY: scmoutput
USE nlsizes_namelist_mod,  ONLY: model_levels

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_conv !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates a parcel ascent from k_plume.
!   An undilute ascent is always calculated. If the option l_dilute is set
!   to .true. a dilute ascent is also calculated. The dilute ascent uses
!   the entrainment rate held in entrain_fraction to mix in environmental air.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------

! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,nSCMDpkgs              ! SCM  - No of diagnostics packages

INTEGER, INTENT(IN) :: &
  nlcl(npnts)          & ! Lifting condensation level
 ,k_plume(npnts)         ! Starting model level for plume ascent

LOGICAL, INTENT(IN) :: &
  l_dilute             & ! .TRUE. if a dilute parcel ascent also required.
 ,L_SCMDiags(nSCMDpkgs)  ! SCM - Logicals for diagnostics packages

REAL, INTENT(IN) ::           &
  t(npnts,model_levels)   & ! Temperature on model levels (K)
 ,Tl(npnts,model_levels)  & ! Liquid water temperature on model lev (K)
 ,qcl(npnts,model_levels) & ! cloud liquid water on model levels (kg/kg)
 ,qcf(npnts,model_levels) & ! cloud ice water on model levels (kg/kg)
 ,qw(npnts,model_levels)  & ! total water (kg/kg)
 ,t_dens_env(npnts,model_levels)
                             ! Density potential temperature of environment (K)

REAL, INTENT(IN) ::                          &
  p_theta_lev(npnts,model_levels)        & ! Pressure on theta levels (Pa)
 ,exner_theta_levels(npnts,model_levels) & ! Exner Pressure on theta levels
 ,z_theta(npnts,model_levels)            & ! Height of theta levels  (m)
 ,entrain_fraction(npnts,model_levels)     ! fraction of environmental air
                                               ! to mix with parcel

REAL, INTENT(INOUT) :: &
  sl_plume(npnts)      & ! SL at start of plume
 ,qw_plume(npnts)      & ! total water at plume start (kg/kg)
 ,t_parc_dil(npnts,model_levels)     ! Dilute Parcel temperature (K)

REAL, INTENT(OUT) ::                   &
  t_parc(npnts,model_levels)       & ! Parcel temperature  (K)
 ,ql_parc(npnts,model_levels)      & ! Parcel water content
 ,buoyancy(npnts,model_levels)     & ! Parcel buoyancy  (K)
 ,buoyancy_dil(npnts,model_levels) & ! Dilute parcel buoyancy (K)
 ,env_svl(npnts,model_levels)      & ! Density (virtual) static energy
                                         ! over CP for layer.
 ,par_svl(npnts,model_levels)      & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,par_svl_dil(npnts,model_levels)  & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,denv_bydz(npnts, model_levels)   & ! Gradient of density potential
                                         ! temperature in the environment.
 ,dpar_bydz(npnts, model_levels)   & ! Gradient of density potential
 ,dqsatdz(npnts, model_levels)       ! dqsat/dz along an adiabat undilute
                                         ! parcel

! Local variables

INTEGER ::        &
 ii,k                ! loop counters

REAL ::           &
  q_liq_env       &  ! Condensed water content of environment.
 ,dq_sat_env      &  ! DQSAT/DT for environment
 ,lrcp_const      &  ! lc or lc+lf over cp
 ,lrcp_const_env  &  ! lc or lc+lf over cp
 ,lrcp_const_parc &  ! lc or lc+lf over cp
 ,l_const         &  ! lc or lc+lf
 ,l_const_env     &  ! lc or lc+lf
 ,dz              &  ! layer depth
 ,dtdz            &  ! temperature gradient along undilute parcel ascent
 ,z_pr            &  ! used in estimating th_ref at next level
 ,th_par          &  ! theta value for parcel
 ,dq_sat_par      &  ! dqsat/dT for parcel
 ,dq_sat_par_dil  &  ! dqsat/dT for parcel
 ,temp_parc       &  ! average temperature of parcel after entrainment
 ,q_vap_parc      &  ! Vapour content of undilute parcel
 ,q_liq_parc      &  ! Liquid water content of undilute parcel
 ,qcl_parc        &  ! parcel qcl - dilute parcel cal
 ,qcf_parc        &  ! parcel qcf - dilute parcel cal
 ,ql_remove          ! ql removed from plume


! parcel calculation

REAL ::                                &
  t_ref(npnts)                         & ! reference temperature
 ,th_ref(npnts)                        & ! reference potential temperature
 ,th_par_km1(npnts)                    & ! theta refernce for level below
 ,qsat_lev(npnts, model_levels)        & ! qsat for reference temperature
 ,qsat_env(npnts, model_levels)        & ! qsat for environment temperature
 ,t_dens_parc(npnts, model_levels)       ! Density potential temperature
                                         ! of parcel.
REAL ::                                &
  max_qw_plume                           ! maximum plume water content

! Arrays added for dilute parcel calculation

REAL ::                               &
  t_ref_dil(npnts)                    & ! dilute parcel reference temperature
 ,th_ref_dil(npnts)                   & ! reference potential temperature
 ,th_par_km_dil(npnts)                & ! reference potential temperature 2nd
 ,qsat_lev_dil(npnts)                 & ! qsat for dilute parcel
 ,qw_parc(npnts,model_levels)         & ! parcel total water undilute plume
 ,qv_parc(npnts,model_levels)         & ! parcel water undilute plume (array not
                                        ! essential but may require in future)
 ,sl_parc(npnts,model_levels)         & ! parcel SL undilute
 ,t_dens_parc_dil(npnts,model_levels) & ! dilute parcel t_dens
 ,ql_parc_dil(npnts,model_levels)       ! dilute parcel liquid water

LOGICAL ::     &
  l_keep_water    ! if true keeps water loading in plume
                  ! false removed if water exceeds 1g/kg

REAL :: TmpScm3d(ScmRowLen, ScmRow, model_levels)    ! Full field

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'PARCEL_ASCENT'



! Model constants

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

l_keep_water  = .FALSE. ! water loading restricted to 1g/kg

! Initialise parcel reference theta and value for level below

DO ii=1, npnts
  th_ref(ii) = tl(ii,k_plume(ii))/exner_theta_levels(ii,k_plume(ii))
  th_par_km1(ii) = th_ref(ii)
END DO

IF (l_dilute) THEN     ! dilute parcel ascent
  DO ii=1, npnts
    th_ref_dil(ii)    = th_ref(ii)
    th_par_km_dil(ii) = th_ref_dil(ii)
  END DO
END IF

!-----------------------------------------------------------------------
! 2.0 Parcel ascent
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted upwards to top of model.
! Dilute parcel ascent - mix in environmental air above lifting
! condensation level
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!-----------------------------------------------------------------------

! Level loop calculating parcel ascent

IF ( l_new_qsat_conv ) THEN
  IF ( l_mr_physics ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, qsat_env, t, p_theta_lev, npnts,           &
!$OMP         l_mr_physics )                                           &
!$OMP PRIVATE( k )
    DO  k = 1,model_levels
      CALL qsat_mix_new(qsat_env(:,k),t(:,k),p_theta_lev(:,k),npnts)
    END DO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, qsat_env, t, p_theta_lev, npnts,           &
!$OMP         l_mr_physics )                                           &
!$OMP PRIVATE( k )
    DO  k = 1,model_levels
      CALL qsat_new(qsat_env(:,k),t(:,k),p_theta_lev(:,k),npnts)
    END DO
!$OMP END PARALLEL DO
  END IF
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( model_levels, qsat_env, t, p_theta_lev, npnts,           &
!$OMP         l_mr_physics )                                           &
!$OMP PRIVATE( k )
  DO  k = 1,model_levels
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qsat_env(1,k),t(1,k),p_theta_lev(1,k),npnts,l_mr_physics)
  END DO
!$OMP END PARALLEL DO
END IF

DO  k = 1,model_levels

  ! Require t_ref on all point for qsat call

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(t_ref, th_ref,                    &
!$OMP exner_theta_levels, k, npnts) PRIVATE(ii)
  DO ii=1, npnts
    t_ref(ii)     = th_ref(ii)*exner_theta_levels(ii,k)
  END DO
!$OMP END PARALLEL DO

  IF ( l_new_qsat_conv ) THEN
    IF ( l_mr_physics ) THEN
      CALL qsat_mix_new(qsat_lev(:,k),t_ref,p_theta_lev(:,k),npnts)
    ELSE
      CALL qsat_new(qsat_lev(:,k),t_ref,p_theta_lev(:,k),npnts)
    END IF
  ELSE
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qsat_lev(1,k),t_ref,p_theta_lev(1,k),npnts,l_mr_physics)
  END IF

  IF (l_dilute) THEN     ! dilute parcel ascent

!$OMP PARALLEL DO DEFAULT(NONE) SHARED(t_ref_dil, th_ref_dil,            &
!$OMP exner_theta_levels, k, npnts) PRIVATE(ii)
    DO ii=1, npnts
      t_ref_dil(ii) = th_ref_dil(ii)*exner_theta_levels(ii,k)
    END DO
!$OMP END PARALLEL DO

    IF ( l_new_qsat_conv ) THEN
      IF ( l_mr_physics ) THEN
        CALL qsat_mix_new(qsat_lev_dil,t_ref_dil,p_theta_lev(:,k),npnts)
      ELSE
        CALL qsat_new(qsat_lev_dil,t_ref_dil,p_theta_lev(:,k),npnts)
      END IF
    ELSE
      ! DEPENDS ON: qsat_mix
      CALL qsat_mix(qsat_lev_dil,t_ref_dil,p_theta_lev(1,k),npnts,            &
                    l_mr_physics)
    END IF

  END IF                 ! dilute parcel

  ! Undilute parcel calculation always required

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ii, lrcp_const, lrcp_const_env,   &
!$OMP l_const, dq_sat_env, dq_sat_par, q_liq_parc, q_liq_env, z_pr,      &
!$OMP lrcp_const_parc, ql_remove, max_qw_plume, q_vap_parc, th_par,      &
!$OMP temp_parc, qcl_parc, qcf_parc, l_const_env, dq_sat_par_dil)

!$OMP DO SCHEDULE(STATIC)
  DO ii=1, npnts

    IF (T_ref(ii) >  tm) THEN
      lrcp_const = lcrcp
      l_const    = lc
    ELSE
      lrcp_const = lsrcp
      l_const    = ls
    END IF
    IF (t(ii,k) >  tm) THEN
      lrcp_const_env = lcrcp
      l_const_env    = lc
    ELSE
      lrcp_const_env = lsrcp
      l_const_env    = ls
    END IF

    dq_sat_env = repsilon*l_const_env*qsat_env(ii,k)/(r*t(ii,k)**2)
    dq_sat_par = repsilon*l_const    *qsat_lev(ii,k)/(r*T_ref(ii)**2)

    q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii,k)                   &
              -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,k)-T_ref(ii) ) &
                                 ) / (1.0+lrcp_const*dq_sat_par) )

    q_liq_env  = MAX( 0.0, ( qw(ii,k) - qsat_env(ii,k)                       &
                -dq_sat_env*( tl(ii,k)               - t(ii,k) )             &
                                 ) / (1.0+Lrcp_const_env*dq_sat_env) )
    !
    ! add on the difference in the environment's ql as calculated by the
    ! UM cloud scheme (using some RH_CRIT value) and what it
    ! would be If RH_CRIT=1. This THEN imitates partial condensation
    ! in the parcel.
    !
    ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env

    t_parc(ii,k) = sl_plume(ii)-gamma_dry*z_theta(ii,k) + lrcp_const         &
                  *ql_parc(ii,k)


    ! May need to recalculate if T_parc is > Tm and T_ref < Tm

    IF (T_ref(ii) <= tm .AND. T_parc(ii,k) >  tm) THEN

      ! recalculate using corrected latent heats
      lrcp_const_parc = lcrcp

      q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii,k)                 &
           -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,k)-T_ref(ii) )    &
                          ) / (1.0+lrcp_const_parc*dq_sat_par) )

      ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k)- q_liq_env

      ! revised at parcel calculation

      t_parc(ii,k)=sl_plume(ii)-gamma_dry*z_theta(ii,k)                      &
                                       +lrcp_const_parc*ql_parc(ii,k)

    END IF       ! test on T_ref and t_parc

    ! Add water removal from undilute plume
    IF (plume_water_load == 1) THEN

      ! water removed from parcel after condensation if > 0.001 kg/kg
      IF (ql_parc(ii,k) > 0.001) THEN
        ql_remove = ql_parc(ii,k) -  0.001
        ql_parc(ii,k) = 0.001
        qw_plume(ii) = qw_plume(ii) - ql_remove
        ! Also adjust sl_plume as well as altered energy
        sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove
      END IF

    ELSE IF (plume_water_load == 2) THEN

      ! water remaining depends on qsat environment
      max_qw_plume = fac_qsat*qsat_env(ii,k)
      max_qw_plume = MIN (qlcrit,max_qw_plume)  ! max
      max_qw_plume = MAX (qlmin,max_qw_plume)   ! min value
      IF (ql_parc(ii,k) > max_qw_plume) THEN
        ql_remove = ql_parc(ii,k) - max_qw_plume
        ql_parc(ii,k) = max_qw_plume
        qw_plume(ii) = qw_plume(ii) - ql_remove
        ! Also adjust sl_plume as well as altered energy
        sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove

      END IF

    END IF    ! end test on plume_water_load

    q_vap_parc=qw_plume(ii)-ql_parc(ii,k)
    qv_parc(ii,k) = q_vap_parc

    t_dens_parc(ii,k)=t_parc(ii,k)*(1.0+c_virtual*q_vap_parc-ql_parc(ii,k))


    ! calculate t_ref for next level
    IF (k >  1 .AND. k <   model_levels-1) THEN
      z_pr = (z_theta(ii,k+1)-z_theta(ii,k))/(z_theta(ii,k)-z_theta(ii,k-1))
      th_par = t_parc(ii,k)/exner_theta_levels(ii,k)
      th_ref(ii) = th_par*(1.0+z_pr) - th_par_km1(ii)*z_pr

      ! Check sensible value otherwise set to previous reference value
      ! Problems can occur near top of model where calculation are nolonger
      ! important.
      IF (th_ref(ii) < 0.0) THEN
        th_ref(ii) = th_par_km1(ii)
      END IF
      IF (th_par > 0.0) THEN
        th_par_km1(ii) = th_par
      END IF
    END IF

    !-----------------------------------------------------------------------
    ! Dilute parcel ascent
    !-----------------------------------------------------------------------
    IF (l_dilute) THEN

      IF (k <= nlcl(ii)) THEN  ! Dilute parcel same as undilute
                               ! parcel ascent (no entrainment)

        sl_parc(ii,k) = sl_plume(ii)
        qw_parc(ii,k) = qw_plume(ii)
        t_parc_dil(ii,k)     = t_parc(ii,k)
        ql_parc_dil(ii,k)    = ql_parc(ii,k)
        t_dens_parc_dil(ii,k)= t_dens_parc(ii,k)
        th_ref_dil(ii)       = th_ref(ii)
        th_par_km_dil(ii)    = th_par_km1(ii)

      ELSE                      ! Dilute parcel ascent now required

        IF (t_ref_dil(ii) >  tm) THEN
          lrcp_const = lcrcp
          l_const    = lc
        ELSE
          lrcp_const = lsrcp
          l_const    = ls
        END IF

        !-----------------------------------------------------------------------
        ! Dilute parcel
        !-----------------------------------------------------------------------
        ! Mix in entrain_fraction from environmental air from level below and
        ! raise this to current level.
        ! Assume mix in fraction of mass from environment.
        ! Estimate parcel properties after mixing air from environment with
        ! parcel. Temperature given approximately by average

        temp_parc = (t_parc_dil(ii,k-1)                                   &
                               + entrain_fraction(ii,k)*t(ii,k-1))        &
                        /(1.0+entrain_fraction(ii,k))


        qw_parc(ii,k) = (qw_parc(ii,k-1) +                                &
                              entrain_fraction(ii,k)*qw(ii,k-1))          &
                         /(1.0+entrain_fraction(ii,k))

        qcl_parc = (ql_parc_dil(ii,k-1)   +                               &
                              entrain_fraction(ii,k)*qcl(ii,k-1))         &
                         /(1.0+entrain_fraction(ii,k))

        qcf_parc = (0.0     +                                             &
                              entrain_fraction(ii,k)*qcf(ii,k-1))         &
                         /(1.0+entrain_fraction(ii,k))

        ! All condensed water either ice or liquid based on t_ref
        sl_parc(ii,k) = temp_parc - lrcp_const*(qcl_parc+qcf_parc)        &
                                +gamma_dry*z_theta(ii,k-1)

        dq_sat_par_dil = repsilon*l_const*qsat_lev_dil(ii)                &
                                    /(r*t_ref_dil(ii)**2)

        q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)            &
      -dq_sat_par_dil*( sl_parc(ii,k)-gamma_dry*z_theta(ii,k)-t_ref_dil(ii)) &
                                   ) / (1.0+lrcp_const*dq_sat_par_dil) )

        ! add on the dIfference in the environment's ql as calculated by the
        ! UM cloud scheme (using some RH_CRIT value) and what it
        ! would be If RH_CRIT=1. This THEN imitates partial condensation
        ! in the parcel.
        !
        ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env

        t_parc_dil(ii,k) = sl_parc(ii,k)-gamma_dry*z_theta(ii,k)             &
                                      +lrcp_const*ql_parc_dil(ii,k)

        ! May need to recalculate if T_parc is > Tm and T_ref < Tm

        IF (t_ref_dil(ii) <= tm .AND. t_parc_dil(ii,k) >  tm) THEN

          ! recalculate using corrected latent heats
          lrcp_const_parc = lcrcp

          q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)          &
       -dq_sat_par_dil*(sl_parc(ii,k)-gamma_dry*z_theta(ii,k)-t_ref_dil(ii)) &
                         ) / (1.0+lrcp_const_parc*dq_sat_par_dil) )

          ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k)                         &
                                              + qcf(ii,k)- q_liq_env

          ! revised at parcel calculation

          t_parc_dil(ii,k)=sl_parc(ii,k)-gamma_dry*z_theta(ii,k)             &
                                   +lrcp_const_parc*ql_parc_dil(ii,k)

        END IF   ! test on t_ref

        q_vap_parc=qw_parc(ii,k)-ql_parc_dil(ii,k)

        ! Water loading in plume  ( option 0 water remains in plume)

        IF (dil_plume_water_load == 1) THEN

          ! water removed from parcel after condesation if > 0.001 kg/kg
          IF (ql_parc_dil(ii,k) > 0.001) THEN
            ql_parc_dil(ii,k) = 0.001
            qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
          END IF

        ELSE IF (dil_plume_water_load == 2) THEN

          ! water remaining depends on qsat environment
          max_qw_plume = 0.5*qsat_env(ii,k)
          max_qw_plume = MIN (qlcrit,max_qw_plume)  ! max
          max_qw_plume = MAX (qlmin,max_qw_plume)   ! min value
          IF (ql_parc_dil(ii,k) > max_qw_plume) THEN
            ql_parc_dil(ii,k) = max_qw_plume
            qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
          END IF

        END IF    ! end test on dil_plume_water_load

        t_dens_parc_dil(ii,k)=t_parc_dil(ii,k)*                           &
                        (1.0+c_virtual*q_vap_parc-ql_parc_dil(ii,k))

        ! calculate dilute t_ref for next level
        IF (k >  1 .AND. k <   model_levels-1) THEN
          z_pr = (z_theta(ii,k+1)-z_theta(ii,k))                          &
                                      /(z_theta(ii,k)-z_theta(ii,k-1))

          th_par = t_parc_dil(ii,k)/exner_theta_levels(ii,k)
          th_ref_dil(ii) = th_par*(1.0+z_pr) - th_par_km_dil(ii)*z_pr
          ! Check new reference sensible
          IF (th_ref_dil(ii) < 0.0) THEN
            th_ref_dil(ii) = th_par_km_dil(ii)
          END IF
          IF (th_par > 0.0) THEN
            th_par_km_dil(ii) = th_par
          END IF
        END IF      ! k level test

      END IF   ! test on LCL
    END IF    ! test on L_dilute


    buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

    env_svl(ii,k)  = t_dens_env(ii,k)  + gamma_dry*z_theta(ii,k)

    par_svl(ii,k)  = t_dens_parc(ii,k) + gamma_dry*z_theta(ii,k)

    IF (L_dilute) THEN
      buoyancy_dil(ii,k) = t_dens_parc_dil(ii,k) - t_dens_env(ii,k)
      par_svl_dil(ii,k)  = t_dens_parc_dil(ii,k) + gamma_dry*z_theta(ii,k)
    END IF

  END DO     ! Loop over ii
!$OMP END DO
!$OMP END PARALLEL

END DO      ! level loop

!$OMP PARALLEL DO DEFAULT(NONE)                                        &
!$OMP SHARED( model_levels, qsat_lev, T_parc, p_theta_lev, npnts,      &
!$OMP         l_mr_physics, z_theta, dpar_bydz, par_svl, denv_bydz,    &
!$OMP         env_svl, ls, dqsatdz, repsilon, r, g, l_new_qsat_conv )  &
!$OMP PRIVATE(dz, l_const, dtdz, ii, k)
DO  k = 1,model_levels
  ! Gradient calculations

  IF (k >= 2) THEN

    IF ( l_new_qsat_conv ) THEN
      IF ( l_mr_physics ) THEN
        CALL qsat_mix_new(qsat_lev(:,k),T_parc(:,k),p_theta_lev(:,k),npnts)
      ELSE
        CALL qsat_new(qsat_lev(:,k),T_parc(:,k),p_theta_lev(:,k),npnts)
      END IF
    ELSE
      ! DEPENDS ON: qsat_mix
      CALL qsat_mix(qsat_lev(1,k),T_parc(1,k),p_theta_lev(1,k)           &
                                                  ,npnts,l_mr_physics)
    END IF

    DO ii=1,npnts
      !-------------------------------------------------------------
      ! Find vertical gradients in parcel and environment SVL
      ! (using values from level below (i.e. K-1)).
      !-------------------------------------------------------------

      dz = z_theta(ii,k) - z_theta(ii,k-1)

      dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz

      denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

      !-----------------------------------------------------------------
      ! Temperature gradient and qsat/dz along an undilute parcel ascent
      !-----------------------------------------------------------------

      dtdz = (T_parc(ii,k) - T_parc(ii,k-1))/dz

      IF (T_parc(ii,k) >  tm) THEN
        l_const=lc
      ELSE
        l_const=ls
      END IF

      dqsatdz(ii,k) = (qsat_lev(ii,k)/(r*T_parc(ii,k)))                &
                            * (dtdz*repsilon*l_const/T_parc(ii,k) + g)

    END DO    ! ii loop

  ELSE

    DO ii=1,npnts
      dpar_bydz(ii,k) = 0.0
      denv_bydz(ii,k) = 0.0
      dqsatdz(ii,k) = 0.0
    END DO    ! ii loop

  END IF   ! test on k

END DO      ! level loop
!$OMP END PARALLEL DO


!-----------------------------------------------------------------------
! SCM diagnostics from convective diagnosis - only if required
!-----------------------------------------------------------------------
! Dilute parcel information not passed back to higher level routine
!-----------------------------------------------------------------------
IF (model_type == mt_single_column) THEN
  IF (l_scmdiags(scmdiag_conv) .AND. l_dilute ) THEN
    DO k=1, model_levels
      TmpScm3d(1,1,k)=ql_parc(1,k)
    END DO
    CALL SCMoutput(TmpScm3d,'ql_parc',                      &
         'parcel water                ','kg/kg',            &
         t_inst,d_wet,default_streams,'',routinename)
    DO k= 1,model_levels
      TmpScm3d(1,1,k)=ql_parc_dil(1,k)
    END DO
    CALL SCMoutput(TmpScm3d,'ql_parc_dil',                  &
         'parcel water dilute plume     ','kg/kg',          &
         t_inst,d_wet,default_streams,'',routinename)

    DO k=1, model_levels
      TmpScm3d(1,1,k)=sl_parc(1,k)
    END DO
    CALL SCMoutput(TmpScm3d,'sl_parc',                      &
         'Parcel energy          ','K',                     &
         t_inst,d_wet,default_streams,'',routinename)

    DO k=1, model_levels
      TmpScm3d(1,1,k)=qw_parc(1,k)
    END DO
    CALL SCMoutput(TmpScm3d,'qw_parc',                      &
         'total parcel water          ','kg/kg',            &
         t_inst,d_wet,default_streams,'',routinename)

  END IF ! scmdiag_conv

END IF ! model_type


!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE parcel_ascent
END MODULE parcel_ascent_mod
