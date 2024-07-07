! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


MODULE cosp_init_mod
  USE cosp_constants_mod, ONLY: n_hydro
  USE cosp_types_mod, ONLY: cosp_config, cosp_gridbox, construct_cosp_gridbox, &
                            cosp_subgrid, cosp_sghydro, construct_cosp_subgrid,&
                            construct_cosp_sghydro
  USE cosp_input_mod
  USE planet_constants_mod, ONLY: planet_radius
  USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,pdims_l,tdims, &
                                  tdims_s,tdims_l
  USE ereport_mod
  USE stash_array_mod, ONLY: sf
  USE mcica_mod, ONLY: tot_subcol_gen
  USE rad_input_mod, ONLY: i_inhom
  USE rad_pcf, ONLY: ip_mcica
  USE mphys_inputs_mod, ONLY: l_mcr_qgraup, l_mcr_qrain

  USE qsat_mod, ONLY: qsat

  IMPLICIT NONE

! Description:
!   Routine that allocates and initialises the derived types with
!   the gridbox-mean inputs and configuration information for COSP.
!
! Method:
!   Sets the logical flags in cosp_cfg that control the COSP outputs.
!   It also deals with the dependencies between diagnostics and instrument
!   simulators.
!   It allocates the arrays in cosp_gbx and fills them in with relevant model
!   variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_init(L_cosp,L_Rad_Step_prog,L_radiation,row_length,rows,     &
      model_levels,cosp_crain_3d,cosp_csnow_3d,p,q_n,T_n,                      &
      t_surf,p_star,p_theta_levels,r_theta_levels,r_rho_levels,                &
      L_cosp_call,L_cosp_run,cosp_npoints,cosp_cfg,cosp_gbx,cosp_sgx,cosp_sgh)

  USE errormessagelength_mod, ONLY: errormessagelength

  IMPLICIT NONE

!-----Input arguments
! Switch that tells if COSP is requested
  LOGICAL,INTENT(IN) :: L_cosp
! Is this a prognostic radiation timestep?
  LOGICAL,INTENT(IN) :: L_Rad_Step_prog
! Is this a radiation timestep
  LOGICAL,INTENT(IN) :: L_radiation
! Domain dimensions
  INTEGER,INTENT(IN) :: row_length,rows,model_levels
! Convective rainfall and snowfall 3D fields
  REAL,INTENT(IN) :: cosp_crain_3d(row_length,rows,model_levels)
  REAL,INTENT(IN) :: cosp_csnow_3d(row_length,rows,model_levels)
! Other model fields needed by COSP
  REAL,INTENT(IN) :: p(pdims_s%i_start:pdims_s%i_end, &
                       pdims_s%j_start:pdims_s%j_end, &
                       pdims_s%k_start:pdims_s%k_end)
  REAL,INTENT(IN) :: q_n(tdims%i_start:tdims%i_end, &
                          tdims%j_start:tdims%j_end,&
                          1:tdims%k_end)
  REAL,INTENT(IN) :: p_theta_levels(tdims_s%i_start:tdims_s%i_end, &
                                    tdims_s%j_start:tdims_s%j_end, &
                                    tdims_s%k_start:tdims_s%k_end)
  REAL,INTENT(IN) :: r_theta_levels(tdims_l%i_start:tdims_l%i_end, &
                                    tdims_l%j_start:tdims_l%j_end, &
                                    tdims_l%k_start:tdims_l%k_end)
  REAL,INTENT(IN) :: r_rho_levels(pdims_l%i_start:pdims_l%i_end, &
                                  pdims_l%j_start:pdims_l%j_end, &
                                  pdims_l%k_start:pdims_l%k_end)
  REAL,INTENT(IN) :: T_n(tdims%i_start:tdims%i_end, &
                         tdims%j_start:tdims%j_end, &
                         1:tdims%k_end)
  REAL,INTENT(IN) :: T_surf(row_length, rows)
  REAL,INTENT(IN) :: p_star(pdims%i_start:pdims%i_end, &
                            pdims%j_start:pdims%j_end)

!-----Output arguments
! COSP needs to be called in this timestep
  LOGICAL,INTENT(OUT) :: L_cosp_call
! COSP run, otherwise used in diagnostic mode
  LOGICAL,INTENT(OUT) :: L_cosp_run
! Number of horizontal points in the COSP domain
  INTEGER,INTENT(OUT) :: cosp_npoints
! COSP Configuration options
  TYPE(cosp_config),INTENT(OUT)  :: cosp_cfg
! Gridbox information. Input for COSP
  TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx
! Subgrid information. Input for COSP
  TYPE(cosp_subgrid),INTENT(OUT) :: cosp_sgx
! Subgrid hydrometeors. Input for COSP
  TYPE(cosp_sghydro),INTENT(OUT) :: cosp_sgh

!-----Local variables
  INTEGER :: i, k, icode, cosp_ncolumns
  REAL :: cosp_time,cosp_time_bnds(2)
  CHARACTER(LEN=errormessagelength) :: cmessage = ' ' ! Error/warning message
  CHARACTER(LEN=9) :: routine_name='COSP_INIT'

  icode = 9

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Set diagnostic flags
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF (L_cosp.AND.L_Rad_Step_prog) THEN
!   Copy instrument flags to cfg structure
    cosp_cfg%Lradar_sim = cosp_cloudsat_sim
    cosp_cfg%Llidar_sim = cosp_lidar_sim
    cosp_cfg%Lisccp_sim = cosp_isccp_sim
    cosp_cfg%Lmisr_sim  = cosp_misr_sim
    cosp_cfg%Lmodis_sim = cosp_modis_sim
    cosp_cfg%Lrttov_sim = cosp_rttov_sim
!   Flag to control output to file. Always false in in-line version
    cosp_cfg%Lwrite_output = .FALSE.
!-------Copy diagnostic flags to cfg structure
!   ISCCP
    cosp_cfg%Lalbisccp       = sf(331,2)
    cosp_cfg%Ltauisccp       = sf(332,2)
    cosp_cfg%Lpctisccp       = sf(333,2)
    cosp_cfg%Lcltisccp       = sf(334,2)
    cosp_cfg%Lmeantbisccp    = sf(335,2)
    cosp_cfg%Lmeantbclrisccp = sf(336,2)
    cosp_cfg%Lclisccp        = sf(337,2)
    cosp_cfg%Lboxptopisccp   = .FALSE.
    cosp_cfg%Lboxtauisccp    = .FALSE.
!   CALIPSO
    cosp_cfg%LlidarBetaMol532 = (sf(340,2).OR.sf(357,2))
    cosp_cfg%Latb532          = sf(341,2)
    cosp_cfg%Latb532gbx       = (sf(355,2).OR.sf(356,2))
    cosp_cfg%LcfadLidarsr532  = (sf(342,2).OR.sf(370,2))
    cosp_cfg%Lclcalipso       = (sf(343,2).OR.sf(371,2))
    cosp_cfg%Lcllcalipso      = sf(344,2)
    cosp_cfg%Lclmcalipso      = sf(345,2)
    cosp_cfg%Lclhcalipso      = sf(346,2)
    cosp_cfg%Lcltcalipso      = sf(347,2)
    cosp_cfg%LparasolRefl     = sf(348,2)
    cosp_cfg%Lclcalipsoliq    = sf(473,2)
    cosp_cfg%Lclcalipsoice    = sf(474,2)
    cosp_cfg%Lclcalipsoun     = sf(475,2)
    cosp_cfg%Lcllcalipsoliq   = sf(476,2)
    cosp_cfg%Lclmcalipsoliq   = sf(477,2)
    cosp_cfg%Lclhcalipsoliq   = sf(478,2)
    cosp_cfg%Lcltcalipsoliq   = sf(479,2)
    cosp_cfg%Lcllcalipsoice   = sf(480,2)
    cosp_cfg%Lclmcalipsoice   = sf(481,2)
    cosp_cfg%Lclhcalipsoice   = sf(482,2)
    cosp_cfg%Lcltcalipsoice   = sf(483,2)
    cosp_cfg%Lcllcalipsoun    = sf(484,2)
    cosp_cfg%Lclmcalipsoun    = sf(485,2)
    cosp_cfg%Lclhcalipsoun    = sf(486,2)
    cosp_cfg%Lcltcalipsoun    = sf(487,2)
    cosp_cfg%Lclcalipsotmp    = .FALSE.
    cosp_cfg%Lclcalipsotmpliq = .FALSE.
    cosp_cfg%Lclcalipsotmpice = .FALSE.
    cosp_cfg%Lclcalipsotmpun  = .FALSE.
!   CloudSat
    cosp_cfg%Lcfaddbze94  = (sf(350,2).OR.sf(372,2))
    cosp_cfg%Ldbze94      = sf(351,2)
    cosp_cfg%Ldbze94gbx   = (sf(353,2).OR.sf(354,2))
!   MISR
    cosp_cfg%LclMISR = sf(360,2)
!   MODIS
    cosp_cfg%Lclmodis      = sf(450,2)
    cosp_cfg%Lcltmodis     = sf(451,2)
    cosp_cfg%Lclwmodis     = sf(452,2)
    cosp_cfg%Lclimodis     = sf(453,2)
    cosp_cfg%Lclhmodis     = sf(454,2)
    cosp_cfg%Lclmmodis     = sf(455,2)
    cosp_cfg%Lcllmodis     = sf(456,2)
    cosp_cfg%Ltautmodis    = sf(457,2)
    cosp_cfg%Ltauwmodis    = sf(458,2)
    cosp_cfg%Ltauimodis    = sf(459,2)
    cosp_cfg%Ltautlogmodis = sf(460,2)
    cosp_cfg%Ltauwlogmodis = sf(461,2)
    cosp_cfg%Ltauilogmodis = sf(462,2)
    cosp_cfg%Lreffclwmodis = sf(463,2)
    cosp_cfg%Lreffclimodis = sf(464,2)
    cosp_cfg%Lpctmodis     = sf(465,2)
    cosp_cfg%Llwpmodis     = sf(466,2)
    cosp_cfg%Liwpmodis     = sf(467,2)
!   CloudSat and CALIPSO
    cosp_cfg%Lclcalipso2    = (sf(349,2).OR.sf(374,2))
    cosp_cfg%Lcltlidarradar = .FALSE.
    cosp_cfg%Lcllidarradar  = (sf(358,2).OR.sf(359,2))
!   RTTOV
    cosp_cfg%Ltbrttov = .FALSE.
!   Other
    cosp_cfg%Lfracout = sf(352,2)

!   Deal with dependencies
    IF (.NOT.cosp_cfg%Lradar_sim) THEN
      cosp_cfg%Lcfaddbze94    = .FALSE.
      cosp_cfg%Lclcalipso2    = .FALSE.
      cosp_cfg%Lcltlidarradar = .FALSE.
      cosp_cfg%Lcllidarradar  = .FALSE.
      cosp_cfg%Ldbze94        = .FALSE.
      cosp_cfg%Ldbze94gbx     = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Llidar_sim) THEN
      cosp_cfg%Latb532 = .FALSE.
      cosp_cfg%LcfadLidarsr532  = .FALSE.
      cosp_cfg%Lclcalipso2      = .FALSE.
      cosp_cfg%Lclcalipso       = .FALSE.
      cosp_cfg%Lclhcalipso      = .FALSE.
      cosp_cfg%Lcllcalipso      = .FALSE.
      cosp_cfg%Lclmcalipso      = .FALSE.
      cosp_cfg%Lcltcalipso      = .FALSE.
      cosp_cfg%Lcltlidarradar   = .FALSE.
      cosp_cfg%Lcllidarradar    = .FALSE.
      cosp_cfg%LparasolRefl     = .FALSE.
      cosp_cfg%LlidarBetaMol532 = .FALSE.
      cosp_cfg%Latb532gbx       = .FALSE.
      cosp_cfg%Lclcalipsoliq    = .FALSE.
      cosp_cfg%Lclcalipsoice    = .FALSE.
      cosp_cfg%Lclcalipsoun     = .FALSE.
      cosp_cfg%Lclcalipsotmp    = .FALSE.
      cosp_cfg%Lclcalipsotmpliq = .FALSE.
      cosp_cfg%Lclcalipsotmpice = .FALSE.
      cosp_cfg%Lclcalipsotmpun  = .FALSE.
      cosp_cfg%Lcllcalipsoliq   = .FALSE.
      cosp_cfg%Lclmcalipsoliq   = .FALSE.
      cosp_cfg%Lclhcalipsoliq   = .FALSE.
      cosp_cfg%Lcltcalipsoliq   = .FALSE.
      cosp_cfg%Lcllcalipsoice   = .FALSE.
      cosp_cfg%Lclmcalipsoice   = .FALSE.
      cosp_cfg%Lclhcalipsoice   = .FALSE.
      cosp_cfg%Lcltcalipsoice   = .FALSE.
      cosp_cfg%Lcllcalipsoun    = .FALSE.
      cosp_cfg%Lclmcalipsoun    = .FALSE.
      cosp_cfg%Lclhcalipsoun    = .FALSE.
      cosp_cfg%Lcltcalipsoun    = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lisccp_sim) THEN
      cosp_cfg%Lalbisccp       = .FALSE.
      cosp_cfg%Lboxptopisccp   = .FALSE.
      cosp_cfg%Lboxtauisccp    = .FALSE.
      cosp_cfg%Lclisccp        = .FALSE.
      cosp_cfg%Lpctisccp       = .FALSE.
      cosp_cfg%Ltauisccp       = .FALSE.
      cosp_cfg%Lcltisccp       = .FALSE.
      cosp_cfg%Lmeantbisccp    = .FALSE.
      cosp_cfg%Lmeantbclrisccp = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lmisr_sim) THEN
      cosp_cfg%LclMISR = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lmodis_sim) THEN
      cosp_cfg%Lclmodis      = .FALSE.
      cosp_cfg%Lcltmodis     = .FALSE.
      cosp_cfg%Lclwmodis     = .FALSE.
      cosp_cfg%Lclimodis     = .FALSE.
      cosp_cfg%Lclhmodis     = .FALSE.
      cosp_cfg%Lclmmodis     = .FALSE.
      cosp_cfg%Lcllmodis     = .FALSE.
      cosp_cfg%Ltautmodis    = .FALSE.
      cosp_cfg%Ltauwmodis    = .FALSE.
      cosp_cfg%Ltauimodis    = .FALSE.
      cosp_cfg%Ltautlogmodis = .FALSE.
      cosp_cfg%Ltauwlogmodis = .FALSE.
      cosp_cfg%Ltauilogmodis = .FALSE.
      cosp_cfg%Lreffclwmodis = .FALSE.
      cosp_cfg%Lreffclimodis = .FALSE.
      cosp_cfg%Lpctmodis     = .FALSE.
      cosp_cfg%Llwpmodis     = .FALSE.
      cosp_cfg%Liwpmodis     = .FALSE.
    END IF
    IF (.NOT.cosp_cfg%Lrttov_sim) THEN
      cosp_cfg%Ltbrttov = .FALSE.
    END IF
    IF ((.NOT.cosp_cfg%Lradar_sim).AND.(.NOT.cosp_cfg%Llidar_sim).AND.     &
        (.NOT.cosp_cfg%Lisccp_sim).AND.(.NOT.cosp_cfg%Lmisr_sim).AND.      &
        (.NOT.cosp_cfg%Lmodis_sim)) THEN
      cosp_cfg%Lfracout = .FALSE.
    END IF

!       Are the appropriate masks requested?
    IF ((sf(343,2).OR.sf(349,2)).AND.(.NOT.sf(320,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.320 along any of "//   &
                   "the following: 2.343,2.349")
    IF ((sf(344,2).OR.sf(476,2).OR.sf(480,2).OR.sf(484,2)).AND.                &
        (.NOT.sf(321,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.321 along any of "//   &
                   "the following: 2.344,2.476,2.480,2.484")
    IF ((sf(345,2).OR.sf(477,2).OR.sf(481,2).OR.sf(485,2)).AND.                &
        (.NOT.sf(322,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.322 along any of "//   &
                   "the following: 2.345,2.477,2.481,2.485")
    IF ((sf(346,2).OR.sf(478,2).OR.sf(482,2).OR.sf(486,2)).AND.                &
        (.NOT.sf(323,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.323 along any of "//   &
                   "the following: 2.346,2.478,2.482,2.486")
    IF ((sf(347,2).OR.sf(479,2).OR.sf(483,2).OR.sf(487,2)).AND.                &
        (.NOT.sf(324,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.324 along any of "//   &
                   "the following: 2.347,2.479,2.483,2.487")
    IF ((sf(371,2).OR.sf(374,2).OR.sf(473,2).OR.sf(474,2).OR.                  &
         sf(475,2)).AND.(.NOT.sf(325,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.325 along any of "//   &
                   "the following: 2.371,2.374,2.473,2.474,2.475")
    IF (sf(358,2).AND.(.NOT.sf(326,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.326 along 2.358")
    IF (sf(359,2).AND.(.NOT.sf(327,2))) &
      CALL ereport(routine_name,icode,"You need stash 2.327 along 2.359")
    IF ((sf(331,2).OR.sf(332,2).OR.sf(333,2).OR. &
      sf(334,2).OR.sf(337,2).OR.sf(360,2)).AND.(.NOT.sf(330,2))) THEN
      cmessage = " You need stash 2.330 along any of "// &
                 "the following: 2.331-2.334,3.337,2.360"
      CALL ereport(routine_name,icode,cmessage)
    END IF
!   Does the stats routines need to be called?
    cosp_cfg%Lstats = .FALSE.
    IF ((cosp_cfg%Lradar_sim).OR.(cosp_cfg%Llidar_sim).OR. &
      (cosp_cfg%Lisccp_sim)) cosp_cfg%Lstats = .TRUE.
!   Check that vertical grid and diagnostics requested match
    IF (cosp_use_vgrid) THEN
      IF (sf(320,2).OR.sf(342,2).OR.sf(343,2).OR.sf(349,2).OR.sf(350,2).OR. &
          sf(353,2).OR.sf(355,2).OR.sf(358,2)) THEN
        cmessage = " COSP: Stash 2.320,2.342,2.343,2.349,2.350,2.353,"//    &
                   " 2.355,2.358 cannot be selected"//       &
                   " with use_vgrid=.TRUE."
        CALL Ereport(routine_name,icode,cmessage)
      END IF
    ELSE
      IF (sf(325,2).OR.sf(354,2).OR.sf(356,2).OR.sf(357,2).OR.sf(359,2).OR. &
          sf(370,2).OR.sf(371,2).OR.sf(372,2).OR.sf(374,2).OR.              &
          ANY(sf(473:487,2))) THEN
        cmessage = " COSP: Stash 2.325,2.354,2.356,2.357,2.359,2.370,"// &
                   " 2.371,2.372,2.374,2.473 to 2.487 cannot be selected"// &
                   " with use_vgrid=.FALSE."
        CALL Ereport(routine_name,icode,cmessage)
      END IF
    END IF
!   Consistency of microphysical options
    IF (cosp_use_precipitation_fluxes .AND. l_mcr_qrain) THEN
      cmessage = routine_name//": COSP_USE_PRECIPITATION_FLUXES must be"// &
                 " .FALSE. if prognostic precipitation is used (l_mcr_qrain)."
      CALL ereport(routine_name,icode,cmessage)
    END IF
    IF (cosp_use_precipitation_fluxes .AND. l_mcr_qgraup) THEN
      cmessage = routine_name//": COSP_USE_PRECIPITATION_FLUXES must be"// &
                 " .FALSE. if prognostic graupel is used (l_mcr_qgraup)."
      CALL ereport(routine_name,icode,cmessage)
    END IF
!   Is McICA used?
    IF (i_inhom /= ip_mcica) THEN
      cmessage = routine_name//": COSP is only compatible with McICA."
      CALL ereport(routine_name,icode,cmessage)
    END IF
  END IF ! L_cosp.and.L_Rad_Step_prog

!     Check if COSP needs to be run or is configured in diagnostic mode
  L_cosp_run = .FALSE.
  IF (cosp_cfg%Lradar_sim .OR. cosp_cfg%Llidar_sim .OR. &
      cosp_cfg%Lisccp_sim .OR. cosp_cfg%Lmodis_sim .OR. &
      cosp_cfg%Lmisr_sim  .OR. cosp_cfg%Lrttov_sim) L_cosp_run = .TRUE.

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Allocate and fill in input variables on full radiation timesteps
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  L_cosp_call = .FALSE.
  IF (L_radiation.AND.L_Rad_Step_prog.AND.L_cosp) THEN
    L_cosp_call = .TRUE.
    cosp_npoints=row_length*rows
    cosp_ncolumns = tot_subcol_gen

!   cosp_time* are irrelevant in in-line version. Set to sensible values.
    cosp_time=1.0
    cosp_time_bnds(1)=0.0
    cosp_time_bnds(2)=2.0
!-------Allocate memory for gridbox type
    CALL construct_cosp_gridbox(cosp_time, cosp_time_bnds,                     &
          cosp_radar_freq, cosp_surface_radar, cosp_use_mie_tables,            &
          cosp_use_gas_abs, cosp_do_ray, cosp_melt_lay, cosp_k2,               &
          cosp_npoints, model_levels, cosp_Ncolumns, N_HYDRO,                  &
          cosp_Nprmts_max_hydro, cosp_Naero, cosp_Nprmts_max_aero,             &
          cosp_npoints, cosp_lidar_ice_type, cosp_isccp_topheight,             &
          cosp_isccp_topheight_direction, cosp_overlap, cosp_emsfc_lw,         &
          cosp_use_precipitation_fluxes,cosp_use_reff, cosp_platform,          &
          cosp_satellite, cosp_instrument, cosp_Nchannels, cosp_ZenAng,        &
          cosp_channels(1:cosp_nchannels), cosp_surfem(1:cosp_nchannels),      &
          cosp_co2,cosp_ch4,cosp_n2o,cosp_co,cosp_gbx)
!-------Here code to populate input structure
! Following EG dimension changes assume that
! COSP is not interested in the tdims zeroth level.
    cosp_gbx%sh        = RESHAPE(q_n,(/cosp_npoints,model_levels/))
    cosp_gbx%p         = RESHAPE(p_theta_levels(1:row_length,1:rows,           &
                                                1:tdims_s%k_end),              &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%ph        = RESHAPE(p(1:row_length,1:rows,:),                     &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%psfc      = RESHAPE(p_star,(/cosp_npoints/))
    cosp_gbx%zlev      = RESHAPE(r_theta_levels(1:row_length,                  &
                                  1:rows,1:tdims_l%k_end),                     &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%zlev      = cosp_gbx%zlev - planet_radius
    cosp_gbx%zlev_half = RESHAPE(r_rho_levels(1:row_length,1:rows,:),          &
                                  (/cosp_npoints,model_levels/))
    cosp_gbx%zlev_half = cosp_gbx%zlev_half - planet_radius
    cosp_gbx%T         = RESHAPE(T_n,(/cosp_npoints,model_levels/))
    cosp_gbx%skt       = RESHAPE(t_surf,(/cosp_npoints/))
!-------Compute relative humidity. Using saturated mixing ratio in cosp_gbx%q.
    DO k=1,model_levels
      CALL qsat(cosp_gbx%q(:,k),cosp_gbx%T(:,k),cosp_gbx%p(:,k),cosp_npoints)
      DO i=1,cosp_npoints
        IF (cosp_gbx%q(i,k)/=0.0) &
              cosp_gbx%q(i,k) = cosp_gbx%sh(i,k)/cosp_gbx%q(i,k)*100.0
      END DO
    END DO ! k, model levels

!   LS rain is filled in later in cosp_main, after call to mycrophysics.
    cosp_gbx%rain_cv = RESHAPE(cosp_crain_3d,(/cosp_npoints,model_levels/))
    cosp_gbx%snow_cv = RESHAPE(cosp_csnow_3d,(/cosp_npoints,model_levels/))

!-------Allocate memory for subgrid types. Minimum allocation when COSP
!       is called in diagnostic mode
    IF (L_cosp_run) THEN
      CALL construct_cosp_subgrid(cosp_npoints, cosp_Ncolumns,                 &
                    model_levels, cosp_sgx)
      CALL construct_cosp_sghydro(cosp_npoints, cosp_Ncolumns, model_levels,   &
                                N_HYDRO, cosp_sgh)
    ELSE
      CALL construct_cosp_subgrid(1, 1, 1, cosp_sgx)
      CALL construct_cosp_sghydro(1, 1, 1, 1, cosp_sgh)
    END IF

  END IF !L_radiation.and.L_Rad_Step_prog.and.L_cosp

  RETURN
  END SUBROUTINE cosp_init
END MODULE cosp_init_mod
