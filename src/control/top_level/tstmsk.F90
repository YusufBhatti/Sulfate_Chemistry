! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Checks version mask and option code in a ppx record
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Subroutine Interface:

SUBROUTINE tstmsk(modl,isec,lmask,ladres,errorstatus,cmessage)

#if defined(RECON)
USE rcf_address_vars_mod, ONLY:                                   &
    iopn,                                                         &
    ispace,                                                       &
    vmsk

USE rcf_o3intp_mod

USE rcf_nlist_recon_technical_mod, ONLY: &
    var_recon
#else
! JULES
USE cstash_mod, ONLY: vmsk, ispace, iopn
USE turb_diff_mod,       ONLY: l_leonard_term
#endif

USE stash_model_mod, ONLY: h_global

USE jules_snow_mod, ONLY: nsmax

USE rimtypes
USE lbc_mod

USE free_tracers_inputs_mod, ONLY: i_free_tracer_lbc,              &
    a_max_trvars,        tracer_a,  l_pv_tracer, l_pv_dyn,         &
    l_pv_adv_only, l_diab_tracer, l_diab_tr_bl, l_diab_tr_rad
USE g_wave_input_mod, ONLY:                                       &
    l_gwd,l_use_ussp

USE dust_parameters_mod, ONLY: l_dust, l_dust_diag, l_twobin_dust
USE run_aerosol_mod, ONLY:                                             &
    l_ocff, l_ocff_surem, l_ocff_hilem,                                &
    l_biomass, l_bmass_surem, l_bmass_hilem, l_bmass_hilem_variable,   &
    l_soot, l_soot_surem, l_soot_hilem, l_nitrate,                     &
    l_sulpc_so2, l_so2_surfem, l_so2_natem, l_so2_hilem, l_sulpc_dms,  &
    l_sulpc_ozone, l_sulpc_online_oxidants, l_use_seasalt_pm,          &
    l_sulpc_nh3, l_nh3_em, l_dms_em
USE carbon_options_mod, ONLY: l_co2_interactive
USE coupling_control_mod, ONLY: l_oasis_icecalve, l_oasis_obgc,        &
    l_oasis, l_oasis_ocean

USE model_domain_mod, ONLY: output_grid_stagger,                       &
                            FH_GridStagger_C,                          &
                            FH_GridStagger_Endgame

USE rad_input_mod, ONLY: l_use_sulpc_indirect_sw, l_use_biogenic, &
                         l_use_seasalt_direct, l_use_seasalt_indirect, &
                         l_use_arclbiom, l_use_arclblck, l_use_arcldust, &
                         l_use_arclocff, l_use_arclsslt, l_use_arclsulp, &
                         l_use_arcldlta

USE mphys_inputs_mod,    ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain,         &
                               l_subgrid_qcl_mp, l_casim, l_psd,               &
                               l_orograin,l_orogrime

USE casim_switches, ONLY: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom,    &
                          l_mp_icenumber, l_mp_snownumber, l_mp_snow3mom,      &
                          l_mp_graupnumber, l_mp_graup3mom,                    &
                          l_mp_activesolliquid, l_mp_activesolrain,            &
                          l_mp_activeinsolice, l_mp_activesolice,              &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,        &
                          l_mp_activeinsolnumber

USE electric_inputs_mod, ONLY: l_use_electric, electric_method
USE cloud_inputs_mod,    ONLY: i_cld_area, i_cld_vn, i_rhcpt
USE pc2_constants_mod,   ONLY: rhcpt_off, acf_off, i_cld_pc2
USE river_inputs_mod,    ONLY: l_rivers, l_inland
USE stochastic_physics_run_mod,  ONLY: l_skeb2, l_spt, i_pert_theta,     &
    i_pert_theta_type, pert_theta_correl_seq
USE cv_run_mod, ONLY: l_ccrad, l_3d_cca, l_conv_hist, l_param_conv,      &
    l_conv_prog_group_1, l_conv_prog_group_2, l_conv_prog_group_3,       &
    l_conv_prog_precip, cnv_cold_pools, ccp_off
USE murk_inputs_mod,  ONLY: l_murk, l_murk_source

USE dynamics_testing_mod, ONLY: l_idealised_data

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ancilcta_namelist_mod, ONLY: l_sstanom
USE bl_option_mod,  ONLY: i_bl_vn, i_bl_vn_1a, off
USE rad_input_mod, ONLY: l_use_cariolle, l_use_tpps_ozone, l_use_grad_corr,&
     l_orog_unfilt, i_ozone_int, l_cca_dp_prog, l_cca_md_prog,             &
     l_cca_sh_prog
USE easyaerosol_option_mod, ONLY: l_easyaerosol_sw, l_easyaerosol_lw

! JULES
USE jules_surface_mod, ONLY: l_flake_model, iscrntdiag, l_aggregate,          &
                             i_aggregate_opt, ISrfExCnvGust, formdrag,        &
                             l_vary_z0m_soil, l_elev_land_ice, l_urban2t
USE jules_hydrology_mod, ONLY: l_top
USE jules_radiation_mod, ONLY: l_spec_albedo, l_snow_albedo, l_albedo_obs,   &
                                l_sea_alb_var_chl, l_embedded_snow
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_sice_multilayers, l_ctile, &
                                 l_sice_meltponds_cice, l_saldep_freeze
USE jules_vegetation_mod, ONLY: l_triffid, l_landuse, can_model, l_nitrogen, &
                                l_trif_crop

USE umPrintMgr
USE UM_ParParams
USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon,      &
                      io3_trop_map, io3_trop_map_masscon

USE eng_corr_inputs_mod, ONLY: l_emcorr
USE ukca_option_mod,   ONLY: tc_lbc_ukca, l_ukca_chem_aero,       &
                             l_ukca_dust, l_ukca_prim_moc
USE ukca_tracer_stash, ONLY: a_max_ukcavars
USE ukca_d1_defs, ONLY: ukca_sect, ukca_glomap_diag_sect, ukca_diag_sect, &
                        ukca_s34_plev, ukca_s50_plev
USE tstmsk_ukca_mod, ONLY: tstmsk_ukca
USE um_stashcode_mod, ONLY: stashcode_glomap_clim_sec
USE tstmsk_glomap_clim_mod, ONLY: tstmsk_glomap_clim
USE submodel_mod, ONLY: internal_model_list, n_internal_model_max, atmos_im
USE h_vers_mod, ONLY: h_vers
USE missing_data_mod,     ONLY: imdi
USE g_wave_input_mod,     ONLY: i_gwd_vn
USE ereport_mod, ONLY: ereport

USE nlsizes_namelist_mod, ONLY: &
    tr_ukca, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE

! Description:
!   Determines whether a diagnostic is available to a particular
!   (internal model,section) version. Also checks the option code IOPN
!   Called by INACTR, PRELIM.

! Method:

! The decimal value of the version mask, VMSK, was read from the
! ppxref file by GETPPX, each (model, section, item). The version
! number of the (model, section) - NMASK - is obtained from the
! H_VERS array.

! The procedure for checking diagnostic availability is as follows:
! (1) Check whether the relevant internal model is included in the
!     submodel configuration, by examining the INTERNAL_MODEL_LIST
!     array; if not, return.
! (2) Check whether NMASK=0 - this implies that the diag is
!      unavailable to any version. If so, return.
! (3) Check whether the diag is available to the specified version.
!     The 'version mask' binary code in the STASHmaster specifies
!     which versions the diag is available to. Eg., if it is available
!     to vns. 1 and 3, but not vn.2, the version mask is 101. VMSK
!     is the decimal equivalent of the version mask (5 in this example).
!     TSTMSK therefore calculates the quantity

!             IMOD = MOD(VMSK,2**NMASK)/2**(NMASK-1)

!       If IMOD=1, the diag is available; if IMOD=0, it isn't.

!   The option code is checked in accordance with the ppxref option
!   code definitions.

!   Note for code developers adding new tests: LMASK is initialised to
!   .TRUE. at top of deck. A series of tests is made and if any one
!   test fails LMASK is reset to .FALSE. Therefore do not reinitialise
!   LMASK to .TRUE. anywhere otherwise you may inadvertently overwrite
!   a preceding .FALSE. setting or you may mislead future code writers
!   who may do so.
!    For this reason, all unnecessary LMASK=.TRUE. were removed at 5.1

!  Code description:
!  Language: Fortran 90.
!  This code is written to UM programming standards version 8.3.

!  System component covered:
!  System task:               Sub-Models Project

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER :: modl    ! Internal model number
INTEGER :: isec    ! Section number

!   Scalar arguments with intent(out):
LOGICAL :: lmask   ! T if diag is available to specified version
LOGICAL :: ladres  ! T if diag is available for addressing primary
CHARACTER(LEN=errormessagelength) :: cmessage

! Local scalars
INTEGER :: nmask   ! Version number for (model,section)
INTEGER :: imod    ! Determines whether diag is available
INTEGER :: twonm   ! Used in calculation of IMOD
INTEGER :: twonm1  ! Used in calculation of IMOD
INTEGER :: COUNT
INTEGER :: i
INTEGER :: n0,n1,n2,n2n1,n3,n4,n5,n6,n7,n8,n9,n10
INTEGER :: n11,n12,n13,n14,n15,n16,n17,n18,n19,n20
INTEGER :: n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
INTEGER :: n10n9
INTEGER :: sum_iopn

#if defined(RECON)
! hard wire l_leonard_term to false since it isn't required 
! in the reconfiguration 
LOGICAL, PARAMETER :: l_leonard_term = .FALSE.
#else
! Hard-wire var_recon to false if not the reconfiguration 
LOGICAL, PARAMETER :: var_recon=.FALSE.
#endif
! ErrorStatus
INTEGER :: errorstatus
CHARACTER(LEN=*), PARAMETER  :: RoutineName='TSTMSK'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

sum_iopn=iopn(1)+iopn(2)+iopn(3)+iopn(4)+iopn(5)+iopn(6)
!--------------------------------------------------------------------
! Check whether internal model is included in submodel configuration
!--------------------------------------------------------------------
COUNT = 0
DO i = 1,n_internal_model_max
  IF (internal_model_list(i) == modl) THEN
    lmask =.TRUE.
    ladres=.TRUE.
    COUNT = COUNT + 1
  END IF
END DO
IF (COUNT == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
END IF

!--------------------------------------------------------------------
! Check whether diagnostic is unavailable to any version
!--------------------------------------------------------------------
nmask=h_vers(modl,isec)
IF (nmask == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
END IF

! If a section version is not correctly set, h_vers returns imdi into nmask.
! This would result in twonm1 being zero below and the calculation of imod
! failing with a divide by zero.
IF (nmask == imdi) THEN
  WRITE(cmessage, '(A26,I6,A19,I6)')'UNDEFINED VERSION NUMBER: ',nmask,    &
       ' FOR ATMOS SECTION:',isec
  errorstatus = 10
  CALL ereport(routinename, errorstatus, cmessage)
END IF

! Determine whether the diag is available to the specified version
twonm  = 2**nmask
twonm1 = 2**(nmask-1)
imod   = MOD(vmsk,twonm)/twonm1
IF (imod == 1) THEN
  lmask =.TRUE.
ELSE IF (imod == 0) THEN
  lmask =.FALSE.
  ladres=.FALSE.
  GO TO 9999
ELSE
  WRITE(umMessage,*)'S: TSTMSK INVALID DECODING OF VMSK',vmsk
  CALL umPrint(umMessage,src='tstmsk')
  WRITE(umMessage,*)'s: ... imod=',imod,'     NMASK=',nmask
  CALL umPrint(umMessage,src='tstmsk')
END IF

!-------------------------------------------------------------------
! Check option codes
!-------------------------------------------------------------------
lmask=.TRUE.
!-------------------------------------------------------------------
! Var reconfiguration now is not restricted to primary fields.
! Therefore need to perform on all sections.
!-------------------------------------------------------------------
IF ( var_recon ) THEN
  n17 = MOD((iopn(4)/10),10)
  IF (n17 /= 1) THEN
    lmask = .FALSE.
  END IF
END IF

IF ((isec == 0) .OR. (isec == 31) .OR.                           &
                     (isec == 33) .OR.                           &
   (isec == 36) .OR. (isec == 37)) THEN
  ! Atmosphere primary field
  IF (sum_iopn /= 0) THEN
    n2n1=MOD (iopn(1),100)
    ! Up to A_MAX_TRVARS=150 free tracers now allowed in section 33
    ! Up to 150 free tracer lbcs in section 36
    ! Up to 150 UKCA tracer lbcs in section 37
    IF (isec == 33 .OR. isec == 36 .OR. isec == 37) THEN
      n2n1=MOD (iopn(1),1000)
    END IF
    n3  =MOD((iopn(1)/100),10)
    n4  =MOD((iopn(1)/1000),10)
    n5  =MOD((iopn(1)/10000),10)
    n6  =MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    n8  =MOD((iopn(2)/100),10)
    ! n10n9 is a 2 digit option code for CLASSIC aerosols
    n10n9 =MOD((iopn(2)/1000),100)
    n11 =MOD( iopn(3),10)
    n12 =MOD((iopn(3)/10),10)
    n13 =MOD((iopn(3)/100),10)
    n14 =MOD((iopn(3)/1000),10)
    n15 =MOD((iopn(3)/10000),10)
    n16 =MOD( iopn(4),10)
    n17 =MOD((iopn(4)/10),10)
    n18 =MOD((iopn(4)/100),10)
    n19 =MOD((iopn(4)/1000),10)
    n20 =MOD((iopn(4)/10000),10)
    n21 =MOD( iopn(5),10)
    n22 =MOD((iopn(5)/10),10)
    n23 =MOD((iopn(5)/100),10)
    n24 =MOD((iopn(5)/1000),10)
    n25 =MOD((iopn(5)/10000),10)
    n26 =MOD( iopn(6),10)
    n27 =MOD((iopn(6)/10),10)
    n28 =MOD((iopn(6)/100),10)
    n29 =MOD((iopn(6)/1000),10)
    n30 =MOD((iopn(6)/10000),10)
    IF ((n2n1 == 99) .AND.                                    &
        (isec /= 33) .AND.                                    &
        (isec /= 36) .AND. (isec /= 37)) THEN
      lmask=.FALSE.

      ! Free tracers (s33) and free tracer LBCs (s36)
    ELSE IF (isec == 33 .OR. isec == 36) THEN
      IF (n2n1 > 0 .AND. n2n1 <= a_max_trvars) THEN
        IF ((isec == 33 .AND. .NOT. tracer_a(n2n1)) .OR.      &
            (isec == 36 .AND. i_free_tracer_lbc(n2n1) == 0)) THEN
          lmask=.FALSE.
        END IF
      END IF

      ! UKCA LBCs
    ELSE IF (isec == 37) THEN
      IF (n2n1 > 0 .AND. n2n1 <= a_max_ukcavars) THEN
        IF (tc_lbc_ukca(n2n1) == 0) THEN
          lmask=.FALSE.
        END IF
      END IF

      ! Start of section 0 option codes
      ! Make n30 a metacharacter to allow more option codes
    ELSE IF (n30 == 0) THEN

      !             n3=1 not to be used because of tracers
      !             n3=2 not to be used
      IF ((n3 == 5) .AND. (.NOT. l_top)) THEN
        lmask=.FALSE.   ! TOPMODEL hydrology
      ELSE IF ((n4 == 1) .AND. (.NOT. l_use_biogenic)) THEN
        lmask=.FALSE.   ! biogenic aerosol
      ELSE IF ((n4 == 2) .AND. (.NOT. l_use_arclbiom)) THEN
        lmask=.FALSE.   ! biomass burning aerosol
      ELSE IF ((n4 == 3) .AND. (.NOT. l_use_arclblck)) THEN
        lmask=.FALSE.   ! black carbon aerosol
      ELSE IF ((n4 == 4) .AND. (.NOT. l_use_arclsslt)) THEN
        lmask=.FALSE.   ! sea salt aerosol
      ELSE IF ((n4 == 5) .AND. (.NOT. l_use_arclsulp)) THEN
        lmask=.FALSE.   ! sulphate aerosol
      ELSE IF ((n4 == 6) .AND. (.NOT. l_use_arcldust)) THEN
        lmask=.FALSE.   ! dust aerosol
      ELSE IF ((n4 == 7) .AND. (.NOT. l_use_arclocff)) THEN
        lmask=.FALSE.   ! organic carbon fossil fuel
      ELSE IF ((n4 == 8) .AND. (.NOT. l_use_arcldlta)) THEN
        lmask=.FALSE.   ! delta aerosol
      ELSE IF ((n4 > 8)) THEN
        lmask=.FALSE.  ! these are not yet defined
      ELSE IF ((n5 == 1) .AND. (formdrag == 0)) THEN
        lmask=.FALSE.   ! orographic roughness
      ELSE IF ((n5 == 2) .AND. (i_gwd_vn == imdi)) THEN
        lmask=.FALSE.   ! orographic gradient
      ELSE IF ((n5 == 3) .AND. (.NOT. l_use_grad_corr)) THEN
        lmask=.FALSE.   ! gradient correction for SW radiation
      ELSE IF ((n5 == 4) .AND. (.NOT. l_orog_unfilt)) THEN
        lmask=.FALSE.   ! unfiltered orography for horizon angles
      ELSE IF ((n5 == 7) .AND. (.NOT. l_use_ussp)) THEN
        ! Could define Logical: l_totppn_softs for this option but only if
        ! total ppn required beyond non-orographic spectral gravity wave scheme
        lmask=.FALSE.   ! Total ppn at Start of TS needed in D1
      ELSE IF ((n6 == 1) .AND.                                   &
          (rimwidtha(rima_type_norm) == 0 .OR. h_global(atmos_im) == 'Y')) THEN
        lmask=.FALSE.   ! limited area model boundary condition
      ELSE IF ((n7 == 1) .AND. (.NOT. l_oasis_ocean)) THEN
        lmask=.FALSE.   ! OASIS coupling to ocean model
      ELSE IF ((n7 == 2) .AND. (.NOT. (l_oasis_ocean .AND. l_oasis_obgc    &
                                       .AND. (n15 == 3)))) THEN
        ! MEDUSA coupling for a field solely contained in section 0.
        ! Note that if (n15 == 3) there is a check below that 
        ! l_co2_interactive is true.
        lmask=.FALSE.
      ELSE IF ((n7 == 3) .AND. (.NOT. (l_oasis_ocean .AND.                 &
                                       l_oasis_icecalve))) THEN
        lmask=.FALSE.   ! oasis iceberg calving
      ELSE IF ((n7 == 4)                                                   &
               .AND. (.NOT. (l_oasis_ocean .AND. l_oasis_obgc              &
                             .AND. (l_dust .OR. l_ukca_dust)))) THEN
        lmask=.FALSE.   ! MEDUSA coupling for total dust deposition rate
      ELSE IF (n7 == 7) THEN
        lmask=.FALSE.   ! other coupling fields currently excluded
      ELSE IF ((n8 == 1) .AND. (.NOT. l_sstanom)) THEN
        lmask=.FALSE.   ! SST anomaly
      ELSE IF ((n8 == 2) .AND. (iscrntdiag /= 2) .AND. (iscrntdiag /= 3)) THEN
        lmask=.FALSE.
      ELSE IF ((n8 == 3) .AND. (isrfexcnvgust == 0)) THEN
        lmask=.FALSE.   ! effect of convective downdraughts on surface exchange
      ELSE IF ((n8 == 4) .AND. (.NOT. l_murk)) THEN
        lmask=.FALSE.   ! total aerosol (murk) for visibility
      ELSE IF ((n8 == 5) .AND. (.NOT. l_murk_source)) THEN
        lmask=.FALSE.   ! total aerosol emissions
      ELSE IF ((n8 == 6) .AND. (                                           &
        .NOT. (l_snow_albedo .OR. l_embedded_snow))   &
                       .AND. (nsmax == 0)) THEN
        lmask=.FALSE.   ! snow albedo
      ELSE IF ((n8 == 7) .AND. (i_bl_vn /= i_bl_vn_1a)) THEN
        lmask=.FALSE.   ! tke closure
      ELSE IF ( ( n8  ==  8) .AND. (.NOT. l_emcorr) ) THEN
        lmask=.FALSE.   ! energy adjustment scheme
      ELSE IF ( (n8 == 9) .AND. (.NOT. l_use_electric) ) THEN
        lmask=.FALSE.   ! electric scheme is in use

        !
        ! n10n9 is used for all CLASSIC aerosol related prognostics
        ! Sulphur cycle (1-19)
        ! Some fields also required for UKCA Aerosol scheme
      ELSE IF ((n10n9 == 1) .AND. (.NOT. l_sulpc_so2)) THEN
        lmask=.FALSE.   ! sulphur dioxide cycle
      ELSE IF ((n10n9 == 2) .AND. ((.NOT. l_sulpc_so2)                &
                                   .OR. (.NOT. l_so2_surfem)) ) THEN
        lmask=.FALSE.   ! surface SO2 emissions
      ELSE IF ((n10n9 == 3) .AND. ((.NOT. l_sulpc_so2)                &
                                   .OR. (.NOT. l_so2_hilem)) ) THEN
        lmask=.FALSE.   ! high level SO2 emissions
      ELSE IF ((n10n9 == 4) .AND. ((.NOT. l_sulpc_so2)                &
                                   .OR. (.NOT. l_so2_natem)) ) THEN
        lmask=.FALSE.   ! natural SO2 emissions
      ELSE IF ((n10n9 == 5) .AND. ((.NOT. l_sulpc_so2)                &
                        .OR.  (.NOT. l_sulpc_dms)) ) THEN
        lmask=.FALSE.   ! dimethyl sulphide in SO2 cycle
      ELSE IF ((n10n9 == 6) .AND. ((.NOT. l_sulpc_so2)                &
                                   .OR.  (.NOT. l_sulpc_dms)          &
                                   .OR.  (.NOT. l_dms_em)) ) THEN
        lmask=.FALSE.   ! dimethyl sulphide emissions (excl sea conc)
      ELSE IF ((n10n9 == 7) .AND. ((.NOT. l_sulpc_so2)                &
                        .OR.  (l_sulpc_online_oxidants)               &
                        .OR.  (.NOT. l_sulpc_ozone)) ) THEN
        lmask=.FALSE.   ! offline ozone oxidant in SO2 cycle
      ELSE IF ((n10n9 == 8) .AND. ((.NOT. l_sulpc_so2)                &
                        .OR.  (.NOT. l_sulpc_ozone)                   &
                        .OR.  (.NOT. l_sulpc_nh3)) ) THEN
        lmask=.FALSE.   ! ozone and ammonia in SO2 cycle
      ELSE IF ((n10n9 == 9) .AND. ((.NOT. l_sulpc_so2)                &
                                   .OR. (.NOT. l_sulpc_ozone)         &
                                   .OR. (.NOT. l_sulpc_nh3)           &
                                   .OR. (.NOT. l_nh3_em)) ) THEN
        lmask=.FALSE.   ! ammonia emissions and O3 in SO2 cycle
      ELSE IF ((n10n9 == 10) .AND. ((.NOT. l_sulpc_so2)               &
                        .OR.  (l_sulpc_online_oxidants)) ) THEN
        lmask=.FALSE.   ! offline oxidants (not ozone)
      ELSE IF ((n10n9 == 11) .AND. ((.NOT. l_sulpc_so2)               &
                                    .OR.  (.NOT. l_sulpc_dms)         &
                                    .OR.  (.NOT. l_dms_em))           &
                             .AND. (.NOT. l_ukca_chem_aero) ) THEN
        lmask=.FALSE.   ! dimethyl sulphide sea water concentration
        ! Soot (21 - 23)
      ELSE IF ((n10n9 == 21) .AND. (.NOT. l_soot)) THEN
        lmask=.FALSE.   ! soot scheme
      ELSE IF ((n10n9 == 22) .AND. ((.NOT. l_soot)                    &
                        .OR.  (.NOT. l_soot_surem)) ) THEN
        lmask=.FALSE.   ! soot scheme with surface emissions
      ELSE IF ((n10n9 == 23) .AND. ((.NOT. l_soot)                    &
                       .OR.  (.NOT. l_soot_hilem)) ) THEN
        lmask=.FALSE.   ! soot scheme with high level emissions
        ! Biomass (24 - 26)
      ELSE IF ((n10n9 == 24) .AND. (.NOT. l_biomass)) THEN
        lmask=.FALSE.   ! biomass scheme
      ELSE IF ((n10n9 == 25) .AND. ((.NOT. l_biomass)                 &
                       .OR. (.NOT. l_bmass_surem)) ) THEN
        lmask=.FALSE.   ! biomass scheme with surface emissions
      ELSE IF ((n10n9 == 26) .AND. ((.NOT. l_biomass)                 &
                         .OR. (.NOT. l_bmass_hilem)) ) THEN
        lmask=.FALSE.   ! biomass scheme with high level emissions
        ! 2 bin and 6 bin Mineral dust
      ELSE IF ((n10n9 == 27) .AND. (.NOT. l_dust)) THEN
        lmask=.FALSE.   ! mineral dust scheme (prognostic)
        ! 6 bin (and currently also 2 bin) mineral dust diagnosis
      ELSE IF ((n10n9 == 28) .AND. (.NOT. l_dust) .AND. (.NOT. l_dust_diag)) THEN
        lmask=.FALSE.   ! mineral dust scheme (diagnostic lifting only)
        ! Nitrate (31)
      ELSE IF ((n10n9 == 31) .AND. (.NOT. l_nitrate)) THEN
        lmask=.FALSE.   ! nitrate scheme
        ! OCFF (35-37)
      ELSE IF ((n10n9 == 35) .AND. (.NOT. l_ocff)) THEN
        lmask=.FALSE.   ! organic carbon fossil fuel scheme
      ELSE IF ((n10n9 == 36) .AND. ((.NOT. l_ocff)                    &
                         .OR. (.NOT. l_ocff_surem)) ) THEN
        lmask=.FALSE.   ! OCFF with surface emissions
      ELSE IF ((n10n9 == 37) .AND. ((.NOT. l_ocff)                    &
                         .OR. (.NOT. l_ocff_hilem)) ) THEN
        lmask=.FALSE.   ! OCFF with high level emissions
        ! 6 bin only Mineral dust
      ELSE IF ((n10n9 == 38) .AND. (l_twobin_dust                    &
                         .OR. (.NOT. l_dust)) ) THEN
        lmask=.FALSE.   ! mineral dust scheme (prognostic)
      ELSE IF ((n10n9 == 39) .AND. ((.NOT. l_biomass)                 &
                         .OR. (.NOT. l_bmass_hilem)                   &
                         .OR. (.NOT. l_bmass_hilem_variable)) ) THEN
        lmask=.FALSE.   ! biomass scheme with high level emissions
                        ! and variable injection heights
        ! End of classic aerosol section

      ELSE IF (n11 == 1) THEN
        lmask=.FALSE.   ! basic vegetation scheme, no longer supported
      ELSE IF ((n11 == 2) .AND. .FALSE.) THEN
        lmask=.FALSE.   ! fractional vegetation scheme is always on
      ELSE IF ((n11 == 3) .AND. (.NOT. l_triffid) ) THEN
        lmask=.FALSE.   ! TRIFFID vegetation scheme
      ELSE IF ((n11 == 4) .AND. ((.NOT. l_albedo_obs) .OR. l_spec_albedo)) THEN
        lmask=.FALSE.   ! SW albedo_obs
      ELSE IF ((n11 == 5) .AND. (.NOT. (l_albedo_obs .AND. l_spec_albedo))) THEN
        lmask=.FALSE.   ! VIS and NIR albedo_obs
      ELSE IF ( (n11 == 6) .AND.                                              &
        ( .NOT. (l_aggregate .AND. (i_aggregate_opt==1)) ) ) THEN
        lmask=.FALSE.   ! No separate prognostic for thermal roughness
      ELSE IF ( (n11 == 7) .AND. (.NOT. l_vary_z0m_soil) ) THEN
        lmask=.FALSE.   ! No prgnostic for bare soil (momentum) roughness length
      ELSE IF ( (n11 == 8) .AND. ((.NOT. l_triffid) .OR. (.NOT. l_landuse )) ) THEN
        lmask=.FALSE.   ! L_TRIFFID and L_LANDUSE both need to be TRUE for the land use prognostics 
      ELSE IF ((n11 == 9) .AND. ((.NOT. l_triffid).OR.(.NOT. l_nitrogen))) THEN
        lmask=.FALSE.   ! use both TRIFFID vegetation and NITROGEN schemes  
      ELSE IF ((n12 == 3) .AND. (.NOT. l_mcr_qcf2)) THEN
        lmask=.FALSE.   ! QCF2 is prognostic
      ELSE IF ((n12 == 4) .AND. (.NOT. l_mcr_qrain)) THEN
        lmask=.FALSE.   ! QRAIN is prognostic
      ELSE IF ((n12 == 5) .AND. (.NOT. l_mcr_qgraup)) THEN
        lmask=.FALSE.   ! QGRAUP is prognostic
      ELSE IF ((n12 == 6) .AND. (i_cld_vn /= i_cld_pc2)) THEN
        lmask=.FALSE.   ! PC2 cloud fraction
      ELSE IF ( (n12 ==7) .AND. (.NOT. l_subgrid_qcl_mp) ) THEN
        lmask=.FALSE.   ! No mixed-phase cloud generated by turbulence.
      ELSE IF ( (n12 ==8) .AND. (.NOT. l_casim) ) THEN
        lmask=.FALSE.   ! CASIM microphysics inactive, therefore not available.
      ELSE IF ((n13 == 1) .AND. (l_3d_cca)) THEN
        lmask=.FALSE.    !  convective cloud amount is 2D
      ELSE IF ((n13 == 2) .AND. (.NOT. l_3d_cca)) THEN
        lmask=.FALSE.    !  convective cloud amount is 3D
      ELSE IF ((n13 == 3) .AND. (.NOT. l_ccrad)) THEN
        lmask=.FALSE.   ! CCRad scheme
      ELSE IF ((n13 == 4) .AND. (cnv_cold_pools <= ccp_off)) THEN
        lmask=.FALSE.   ! prognostics for convective cold pools
      ELSE IF ((n13 == 5) .AND. (.NOT. l_conv_hist)) THEN
        lmask=.FALSE.   ! fields for convection history
      ELSE IF ((n13 == 6) .AND. (.NOT. l_conv_prog_group_1)) THEN
        lmask=.FALSE.   ! Convection experimental prognostic fields group 1
      ELSE IF ((n13 == 7) .AND. (.NOT. l_conv_prog_group_2)) THEN
        lmask=.FALSE.   ! Convection experimental prognostic fields group 2
      ELSE IF ((n13 == 8) .AND. (.NOT. l_conv_prog_group_3)) THEN
        lmask=.FALSE.   ! Convection experimental prognostic fields group 3
      ELSE IF ((n13 == 9) .AND. (.NOT. l_conv_prog_precip)) THEN
        lmask=.FALSE.   ! Convection prognostic recent precipitation field
      ELSE IF ( (n14 == 1) .AND. ((.NOT. l_triffid) .OR.                  & 
        (.NOT. l_trif_crop )) ) THEN
        lmask=.FALSE.   ! L_TRIFFID and L_TRIF_CROP both need to be TRUE 
                        ! for the pasture prognostics 
      ELSE IF ( (n14 == 2) .AND. (l_triffid) .AND.                  & 
        (l_trif_crop ) ) THEN
        lmask=.FALSE.   ! Total disturbed veg fraction prognostics should 
                        ! be turned off when L_TRIFFID and L_TRIF_CROP are 
                        ! true - separate crop and pasture fractions are 
                        ! used instead
      ELSE IF ((n15 == 3) .AND. (.NOT. l_co2_interactive)) THEN
        lmask=.FALSE.   ! carbon cycle scheme
      ELSE IF ((n15 == 4) .AND. ((.NOT. l_co2_interactive) .AND.            &
        (.NOT. l_oasis_obgc)) ) THEN
        lmask=.FALSE.   ! interactive CO2 or MEDUSA coupling of CO2 flux 
      ELSE IF ((n16 == 2) .AND. (nsmax==0)) THEN
        lmask=.FALSE.   ! JULES snow scheme with multiple snow layers
      ELSE IF ((n16 == 3) .AND. (                                         &
        .NOT. (l_snow_albedo .OR. l_embedded_snow))) THEN
        lmask=.FALSE.   ! snow soot
      ELSE IF ((n16 == 4) .AND. (.NOT. l_urban2t)) THEN
        lmask=.FALSE.   ! URBAN-2T schemes including MORUSES
      ELSE IF ((n16 == 5) .AND. (.NOT. l_flake_model)) THEN
        lmask=.FALSE.   ! FLake lake scheme
      ELSE IF ((n16 == 6) .AND. (.NOT. l_elev_land_ice)) THEN
        lmask=.FALSE.   ! Tiled subsurface for elevated tiles
        ! n17 is left out here, as it is for all sections not just prognostics.
      ELSE IF ((n18 == 1) .AND. (.NOT. l_ctile)) THEN
        lmask=.FALSE.   ! coastal tiling scheme
      ELSE IF ((n18 == 2) .AND. (                                         &
         .NOT. (l_sea_alb_var_chl .OR.  l_ukca_prim_moc))) THEN
        lmask=.FALSE.   ! variable open sea chlorophyll concentration
      ELSE IF ((n18 == 3) .AND. (.NOT. l_rivers)) THEN
        lmask=.FALSE.   ! river routing scheme
      ELSE IF ((n18 == 4) .AND. (.NOT. l_inland)) THEN
        lmask=.FALSE.   ! inland re-routing scheme
      ELSE IF ((n18 == 6) .AND. (nice_use > 1)) THEN
        lmask=.FALSE.   ! single sea ice category
      ELSE IF ((n18 == 7) .AND. (nice == 1)) THEN
        lmask=.FALSE.   ! sea ice categories
      ELSE IF ((n18 == 8) .AND. (nice_use == 1)) THEN
        lmask=.FALSE.   ! sea ice categories used fully
      ELSE IF ((n18 == 9) .AND. (.NOT. l_sice_multilayers)) THEN
        lmask=.FALSE.   ! multilayer sea ice scheme
      ELSE IF ((n19 == 1) .AND. (.NOT. l_use_tpps_ozone)) THEN
        lmask=.FALSE.   ! tropopause ozone scheme
      ELSE IF ((n20 == 1) .AND. (can_model /= 4)                  &
                         .AND. (nsmax == 0)) THEN
        lmask=.FALSE.   ! snow canopy scheme
      ELSE IF ((n20 == 2) .AND. (.NOT. l_sice_meltponds_cice)) THEN
        lmask=.FALSE.   ! l_sice_meltponds_cice: sea ice meltponds used in
                        ! albedo calculation
      ELSE IF ((n20 == 3) .AND. (.NOT. l_saldep_freeze)) THEN
        lmask=.FALSE.   ! Salinity-dependent sea surface freezing temperature
      ELSE IF ((n21 == 1) .AND. (.NOT. l_cca_dp_prog)) THEN
        lmask=.FALSE.   ! cca_dp used in cloud generator
      ELSE IF ((n21 == 2) .AND. (.NOT. l_cca_md_prog)) THEN
        lmask=.FALSE.   ! cca_md used in cloud generator
      ELSE IF ((n21 == 3) .AND. (.NOT. l_cca_sh_prog)) THEN
        lmask=.FALSE.   ! cca_sh used in cloud generator
      ELSE IF ((n22 == 1) .AND. (.NOT. l_mp_cloudnumber)) THEN
        lmask=.FALSE.   ! Prognostic cloud number not in use
      ELSE IF ((n22 == 2) .AND. (.NOT. l_mp_rainnumber)) THEN
        lmask=.FALSE.   ! Prognostic rain number not in use
      ELSE IF ((n22 == 3) .AND. (.NOT. l_mp_rain3mom)) THEN
        lmask=.FALSE.   ! Prognostic rain 3rd moment not in use
      ELSE IF ((n22 == 4) .AND. (.NOT. l_mp_icenumber)) THEN
        lmask=.FALSE.   ! Prognostic ice number not in use
      ELSE IF ((n22 == 5) .AND. (.NOT. l_mp_snownumber)) THEN
        lmask=.FALSE.   ! prognostic snow number not in use
      ELSE IF ((n22 == 6) .AND. (.NOT. l_mp_snow3mom)) THEN
        lmask=.FALSE.   ! Prognostic snow 3rd moment not in use
      ELSE IF ((n22 == 7) .AND. (.NOT. l_mp_graupnumber)) THEN
        lmask=.FALSE.   ! Prognostic graupel number not in use
      ELSE IF ((n22 == 8) .AND. (.NOT. l_mp_graup3mom)) THEN
        lmask=.FALSE.   ! Prognostic graupel third moment not in use
      ELSE IF ((n23 == 1) .AND. (.NOT. l_mp_activesolliquid)) THEN
        lmask=.FALSE.   ! Prognostic for soluble liquid not in use 
      ELSE IF ((n23 == 2) .AND. (.NOT. l_mp_activesolrain)) THEN
        lmask=.FALSE.   ! Prognostic for soluble rain not in use 
      ELSE IF ((n23 == 3) .AND. (.NOT. l_mp_activeinsolice)) THEN
        lmask=.FALSE.   ! Prognostic for insoluble ice not in use 
      ELSE IF ((n23 == 4) .AND. (.NOT.  l_mp_activesolice)) THEN
        lmask=.FALSE.   ! Prognostic for soluble ice not in use
      ELSE IF ((n23 == 5) .AND. (.NOT. l_mp_activeinsolliquid)) THEN
        lmask=.FALSE.   ! Prognostic for insoluble liquid not in use 
      ELSE IF ((n23 == 6) .AND. (.NOT. l_mp_activesolnumber)) THEN
        lmask=.FALSE.   ! Prognostic for soluble number not in use
      ELSE IF ((n23 == 7) .AND. (.NOT. l_mp_activeinsolnumber)) THEN
        lmask=.FALSE.   ! Prognostic for insoluble number not in use 
      ELSE IF ((n24 == 1) .AND. (.NOT. l_pv_tracer)) THEN
        lmask=.FALSE.   ! Disable PV-tracers prognostics if not used,
      ELSE IF ((n24 == 2) .AND. (.NOT. l_pv_dyn)) THEN
        lmask=.FALSE.   ! Disable dynamics PV-tracers prognostics if not used
      ELSE IF ((n24 == 3) .AND. (.NOT. l_pv_adv_only)) THEN
        lmask=.FALSE.   ! Disable advection only PV-tracer prog if not used
      ELSE IF ((n24 == 4) .AND. (.NOT. l_diab_tracer)) THEN
        lmask=.FALSE.   ! Disable diabatic tracer is scheme is not used
      ELSE IF ((n24 == 5) .AND. (.NOT. l_diab_tr_bl)) THEN
        lmask=.FALSE.   ! Disable MIX and LH BL diabatic tracers
      ELSE IF ((n24 == 6) .AND. (.NOT. l_diab_tr_rad)) THEN
        lmask=.FALSE.   ! Disable SW and LW rad diabatic tracers
      ELSE IF (n25 == 3) THEN
        lmask=.FALSE.   ! Direct PAR prognostic not currently used by any scheme
      ELSE IF ((n26 == 1) .AND. ( i_pert_theta == off .OR.              &
          i_pert_theta_type /= pert_theta_correl_seq )) THEN
        lmask = .FALSE. ! Time correlated stochastic BL perturbations
                        ! switched off
      ELSE IF ((n28 == 1) .AND. (.NOT. l_use_cariolle)) THEN
        lmask=.FALSE.   ! cariolle ozone scheme
      ELSE IF ((n29 == 1) .AND.                                  &
                         (output_grid_stagger == FH_GridStagger_C)) THEN
        lmask=.FALSE.   ! Disable ENDGame prognostics in ND run
      ELSE IF ((n29 == 2) .AND.                                  &
                         (output_grid_stagger == FH_GridStagger_Endgame)) THEN
        lmask=.FALSE.   ! Disable unecessary ND prognostics in EG run
      END IF

    ELSE IF (n30 == 1) THEN
      ! Can be used for a new set of option codes with n30=1
      ! when those with n30=0 are fully used up.
    END IF ! n30
  END IF ! SUM_IOPN
END IF ! isec

! Print out details of any tracer boundary conditions.
IF ( (isec == 36) .AND. printstatus >= prstatus_diag) THEN
  WRITE(umMessage,*) 'TSTMSK:',isec,n2n1,i_free_tracer_lbc(n2n1),lmask
  CALL umPrint(umMessage,src='tstmsk')
END IF
IF ( (isec == 37)  .AND. printstatus >= prstatus_diag) THEN
  WRITE(umMessage,*) 'TSTMSK:',isec,n2n1,tc_lbc_ukca(n2n1),lmask
  CALL umPrint(umMessage,src='tstmsk')
END IF
! UKCA - has own subroutine
IF (isec == ukca_sect .OR. isec == ukca_glomap_diag_sect .OR. &
    isec == ukca_diag_sect .OR. isec == ukca_s34_plev .OR.    &
    isec == ukca_s50_plev) THEN
  ! calculate option codes
  n1 = MOD (iopn(1)    ,10)
  n2 = MOD((iopn(1)/10),10)
  n3 = MOD((iopn(1)/100),10)
  n4  =MOD((iopn(1)/1000),10)
  n5  =MOD((iopn(1)/10000),10)
  n6  =MOD( iopn(2),10)
  n7  =MOD((iopn(2)/10),10)
  n8  =MOD((iopn(2)/100),10)
  n9  =MOD((iopn(2)/1000),10)
  n10 =MOD((iopn(2)/10000),10)
  n11 =MOD( iopn(3),10)
  n12 =MOD((iopn(3)/10),10)
  n13 =MOD((iopn(3)/100),10)
  n14 =MOD((iopn(3)/1000),10)
  n15 =MOD((iopn(3)/10000),10)
  n16 =MOD( iopn(4),10)
  n17 =MOD((iopn(4)/10),10)
  n18 =MOD((iopn(4)/100),10)
  n19 =MOD((iopn(4)/1000),10)
  n20 =MOD((iopn(4)/10000),10)
  n21 =MOD( iopn(5),10)
  n22 =MOD((iopn(5)/10),10)
  n23 =MOD((iopn(5)/100),10)
  n24 =MOD((iopn(5)/1000),10)
  n25 =MOD((iopn(5)/10000),10)
  n26 =MOD( iopn(6),10)
  n27 =MOD((iopn(6)/10),10)
  n28 =MOD((iopn(6)/100),10)
  n29 =MOD((iopn(6)/1000),10)
  n30 =MOD((iopn(6)/10000),10)

  ! Call tstmsk_ukca to set lmask
  CALL tstmsk_ukca(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,            &
  n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,    &
  n29,n30,lmask)
END IF

! GLOMAP_MODE NWP CLIM - has own subroutine
IF ( isec == stashcode_glomap_clim_sec ) THEN
   ! calculate option codes
   n1 = MOD (iopn(1)    ,10)
   n2 = MOD((iopn(1)/10),10)
   n3 = MOD((iopn(1)/100),10)
   n4  =MOD((iopn(1)/1000),10)
   n5  =MOD((iopn(1)/10000),10)
   n6  =MOD( iopn(2),10)
   n7  =MOD((iopn(2)/10),10)
   n8  =MOD((iopn(2)/100),10)
   n9  =MOD((iopn(2)/1000),10)
   n10 =MOD((iopn(2)/10000),10)
   n11 =MOD( iopn(3),10)
   n12 =MOD((iopn(3)/10),10)
   n13 =MOD((iopn(3)/100),10)
   n14 =MOD((iopn(3)/1000),10)
   n15 =MOD((iopn(3)/10000),10)
   n16 =MOD( iopn(4),10)
   n17 =MOD((iopn(4)/10),10)
   n18 =MOD((iopn(4)/100),10)
   n19 =MOD((iopn(4)/1000),10)
   n20 =MOD((iopn(4)/10000),10)
   n21 =MOD( iopn(5),10)
   n22 =MOD((iopn(5)/10),10)
   n23 =MOD((iopn(5)/100),10)
   n24 =MOD((iopn(5)/1000),10)
   n25 =MOD((iopn(5)/10000),10)
   n26 =MOD( iopn(6),10)
   n27 =MOD((iopn(6)/10),10)
   n28 =MOD((iopn(6)/100),10)
   n29 =MOD((iopn(6)/1000),10)
   n30 =MOD((iopn(6)/10000),10)
   
   ! Call tstmsk_glomap_clim to set lmask
   CALL tstmsk_glomap_clim(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,          &
        n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,    &
        n29,n30,lmask)
END IF

! End of atmos primary block

! Atmos diagnostics
IF (isec == 1) THEN
  ! Shortwave radiation
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1)    ,10)
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n4 = MOD((iopn(1)/1000),10)
    n5 = MOD((iopn(1)/10000),10)
    n6 = MOD( iopn(2),10)
    n8  =MOD((iopn(2)/100),10)
    n9 = MOD((iopn(2)/1000),10)
    IF ((n1 == 1) .AND. (h_global(atmos_im) /= 'Y')) THEN
      lmask=.FALSE.
    ELSE IF ((n4  ==  1) .AND. (.NOT. l_use_sulpc_indirect_sw)) THEN
      lmask = .FALSE.
    ELSE IF ((n5  ==  1) .AND.                                     &
            (.NOT. (l_use_seasalt_direct                           &
              .OR. l_use_seasalt_indirect))) THEN
      lmask = .FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n8 == 1) .AND. (nice_use /= nice)) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 1) .AND. (.NOT. l_easyaerosol_sw)) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 2) THEN
  ! Longwave radiation
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1)   ,10)
    IF ((n1 == 1) .AND.                                            &
        (i_ozone_int /= io3_trop_map_masscon)) THEN
      lmask=.FALSE.
    END IF
    n6 = MOD( iopn(2),10)
    n8  =MOD((iopn(2)/100),10)
    n9 = MOD((iopn(2)/1000),10)
    IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n8 == 1) .AND. (nice_use /= nice)) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 1) .AND. (.NOT. l_easyaerosol_lw)) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 3) THEN
  ! Boundary layer
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1)    ,10)
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n4 = MOD((iopn(1)/1000),10)
    n5 = MOD((iopn(1)/10000),10)
    n6 = MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    IF ((n1 == 1) .AND. (formdrag == 0)) THEN
      lmask=.FALSE.
    END IF
    IF ((n2 == 1) .AND. (.NOT. l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2) .AND. (.NOT. l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3) .AND. (.NOT. l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4) .AND. (.NOT. l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5) .AND. (.NOT. l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6) .AND. (.NOT. l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7) .AND. (.NOT. l_nitrate)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 8) .AND. (.NOT. l_dust) .AND. (.NOT. l_dust_diag)) THEN
      lmask=.FALSE.
    END IF
    IF ((n3 == 1) .AND. (.NOT. l_co2_interactive)) THEN
      lmask=.FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    IF ((n4 == 1) .AND. (nice == 1)) THEN
      lmask = .FALSE.
    ELSE IF ((n4 == 2) .AND. (nice_use == 1)) THEN
      lmask = .FALSE.
    END IF
    IF ((n5 == 1) .AND. (.NOT. l_leonard_term)) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust) .AND. (.NOT. l_dust_diag)) THEN
      lmask = .FALSE.
       ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT. l_dust) .AND. (.NOT. l_dust_diag)) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 4) THEN
  ! Large-scale precipitation
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n5 = MOD((iopn(1)/10000),10)
    n6 = MOD( iopn(2),10)
    n7 = MOD((iopn(2)/10),10)
    n8 = MOD((iopn(2)/100),10)
    n9 = MOD((iopn(2)/1000),10)
    IF ((n2 == 1) .AND. (.NOT. l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2) .AND. (.NOT. l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3) .AND. (.NOT. l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4) .AND. (.NOT. l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5) .AND. (.NOT. l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6) .AND. (.NOT. l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7) .AND. (.NOT. l_nitrate)) THEN
      lmask=.FALSE.
    END IF
    ! Mixed phase scheme indicator
    IF ( ( n5  ==  1 ) .AND. ( .NOT. l_subgrid_qcl_mp ) ) THEN
      lmask=.FALSE.
    END IF
    ! PC2 indicator
    IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust)) THEN
      lmask = .FALSE.
      ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT. l_dust)) ) THEN
      lmask = .FALSE.
    END IF
    
    ! CASIM Microphysics option codes
    IF ( ( n8 == 1 ) .AND. ( .NOT. l_casim ) ) THEN
      lmask = .FALSE.
    ELSE IF ( (n8 == 2 ) .AND. ( ( .NOT. l_casim )          .OR. &
                                 ( .NOT. l_mp_cloudnumber ) .OR. &
                                 ( .NOT. l_mp_icenumber ) ) ) THEN
      lmask = .FALSE.
    ELSE IF ( (n8 == 3 ) .AND. ( ( .NOT. l_casim )          .OR. &
                                 ( .NOT. l_mp_rainnumber  ) .OR. &
                                 ( .NOT. l_mp_icenumber ) ) ) THEN
      lmask = .FALSE.
    ELSE IF ( (n8 == 4 ) .AND. ( ( .NOT. l_casim )          .OR. &
                                 ( .NOT. l_mp_icenumber ) ) ) THEN
      lmask = .FALSE.
    ELSE IF ( (n8 == 5 ) .AND. ( ( .NOT. l_casim )          .OR. &
                                 ( .NOT. l_mp_snownumber ) ) ) THEN
      lmask = .FALSE.
    ELSE IF ( (n8 == 6 ) .AND. ( ( .NOT. l_casim )          .OR. &
                                 ( .NOT. l_mp_graupnumber ) ) ) THEN
      lmask = .FALSE. 
    END IF

    ! Option code for STASH section 4 diagnostics which should only be
    ! output when an appropriate prognostic hydrometeor variable is selected
    ! for the precipitation scheme (and is in use throughout the whole model). 
    IF ( ( n9 == 1 ) .AND. ( .NOT. l_mcr_qgraup ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 2 ) .AND. ( .NOT. l_mcr_qrain ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 3 ) .AND. ( .NOT. l_mcr_qcf2 ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 4 ) .AND. ( l_psd .OR. ( .NOT. l_mcr_qcf2 ) ) ) THEN
      ! Large Scale Precipitation has the option to output two sets of 
      ! diagnostics - one for 'aggregates' and one for 'crystals'. Whenever
      ! the Generic Ice Particle Size Distribution (logical l_psd) is in use
      ! or the second ice prognostic (l_mcr_qcf2) is not in use, some crystal
      ! diagnostics will not be produced by the precipitation scheme.
      ! Therefore, turn off these diagnostics if this is the case.
      lmask = .FALSE.
    ELSE IF ( ( n9 == 5 ) .AND. ( .NOT. l_orograin ) ) THEN
      lmask = .FALSE.
    ELSE IF ( ( n9 == 6 ) .AND. ( .NOT. l_orogrime ) ) THEN
      lmask = .FALSE.
    END IF

  END IF

ELSE IF (isec == 5) THEN
  ! Convection
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n6 = MOD( iopn(2),10)
    n7  =MOD((iopn(2)/10),10)
    IF ((n2 == 1) .AND. (.NOT. l_sulpc_so2)) THEN
      lmask=.FALSE.
    ELSE IF ( (n3 == 1) .AND. (l_3d_cca) ) THEN
      lmask=.FALSE.
    ELSE IF ( (n3 == 2) .AND. (.NOT. l_3d_cca) ) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 2) .AND. (.NOT. l_sulpc_nh3)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 3) .AND. (.NOT. l_soot)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 4) .AND. (.NOT. l_biomass)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 5) .AND. (.NOT. l_dust)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 6) .AND. (.NOT. l_ocff)) THEN
      lmask=.FALSE.
    ELSE IF ((n2 == 7) .AND. (.NOT. l_nitrate)) THEN
      lmask=.FALSE.
    END IF
    IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
    ! 2 and 6 bin dust fields
    IF ((n7 == 1) .AND. (.NOT. l_dust)) THEN
      lmask = .FALSE.
      ! 6 bin dust only fields
    ELSE IF ((n7 == 2) .AND. (l_twobin_dust              &
                        .OR. (.NOT. l_dust)) ) THEN
      lmask = .FALSE.
    END IF
  END IF
  n4 = MOD((iopn(1)/1000),10)
  IF ( (n4 /= 1) .AND. (.NOT. l_param_conv) ) THEN
    lmask=.FALSE.
  END IF

ELSE IF (isec == 6) THEN
  ! Gravity Wave Drag parametrization
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    IF ((n2 == 1) .AND. (.NOT. l_gwd)) THEN
      lmask=.FALSE.
    ELSE IF ((n3 == 1) .AND. (.NOT. l_use_ussp)) THEN
      lmask=.FALSE.
    END IF
  END IF

ELSE IF (isec == 8) THEN
  ! Hydrology
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1),10)
    n22 = MOD((iopn(5)/10),10)
    IF ((n1 == 1) .AND. (.NOT. l_top)) THEN
      lmask=.FALSE.
    ELSE IF ((n22 == 1) .AND. (.NOT. l_rivers)) THEN
      lmask=.FALSE.    !River routing switched off
    ELSE IF ((n22 == 2) .AND. (.NOT. l_inland)) THEN
      lmask=.FALSE.    !Inland re-routing switched off
    END IF
  END IF

ELSE IF (isec == 9) THEN
  ! Cloud parametrization
  IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)
    n3 = MOD((iopn(1)/100),10)
    n6 = MOD( iopn(2),10)
    IF ((n2 == 1) .AND. (i_cld_area == acf_off)) THEN
      lmask=.FALSE.
    ELSE IF ((n3 == 1) .AND. (i_rhcpt == rhcpt_off)) THEN
      lmask=.FALSE.
    ELSE IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 12) THEN
  ! Dynamics Advection
  IF (sum_iopn /= 0) THEN
    n6 = MOD( iopn(2),10)
    IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 16) THEN
  ! Extra physics
  IF (sum_iopn /= 0) THEN
    n2n1 = MOD(iopn(1),100)
    n6 = MOD( iopn(2),10)
    IF ((n2n1 /= 0) .AND. (.NOT. tracer_a(n2n1))) THEN
      lmask=.FALSE.
    END IF
    IF ( ( n6  ==  1 ) .AND. ( i_cld_vn /= i_cld_pc2 ) ) THEN
      lmask = .FALSE.
    END IF
  END IF

ELSE IF (isec == 17) THEN
  ! CLASSIC Aerosol section
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1),     10)
    n2 = MOD((iopn(1)/10),10)
    ! Sulphur cycle diagnostics
    IF ((n1 == 1) .AND. (.NOT. l_sulpc_so2)) THEN
      lmask=.FALSE.
      ! Soot diagnostics
    ELSE IF ((n1 == 2) .AND. (.NOT. l_soot)) THEN
      lmask=.FALSE.
      ! Biomass aerosol diagnostics
    ELSE IF ((n1 == 3) .AND. (.NOT. l_biomass)) THEN
      lmask=.FALSE.
      ! Dust aerosol diagnostics
    ELSE IF ((n1 == 4) .AND. (.NOT. l_dust)) THEN
      lmask=.FALSE.
      ! OCFF aerosol diagnostics
    ELSE IF ((n1 == 5) .AND. (.NOT. l_ocff)) THEN
      lmask=.FALSE.
      ! SOA aerosol diagnostics
    ELSE IF ((n1 == 6) .AND. (.NOT. l_use_biogenic)) THEN
      lmask=.FALSE.
      ! Sea-salt PM diagnostics
    ELSE IF ((n1 == 7) .AND. (.NOT. l_use_seasalt_pm)) THEN
      lmask=.FALSE.
      ! Nitrate aerosol diagnostics
    ELSE IF ((n1 == 8) .AND. (.NOT. l_nitrate)) THEN
      lmask=.FALSE.
      ! aerosol PM10/PM2.5 diagnostics - any of the above mean it is available
    ELSE IF ((n1 == 9) .AND. (.NOT.                                   &
           (l_sulpc_so2      .OR. l_soot .OR. l_biomass      .OR.   &
            l_dust           .OR. l_ocff .OR. l_use_biogenic .OR.   &
            l_use_seasalt_pm .OR. l_nitrate) ) ) THEN
      lmask=.FALSE.
    END IF
    ! The following option codes are usually used
    ! in conjunection with n1 ==1
    ! DMS diagnostics
    IF ((n2 == 1) .AND. (.NOT. l_sulpc_dms)) THEN
      lmask=.FALSE.
      ! Diagnostics which depend on oxidation of sulphur by ozone
    ELSE IF ((n2 == 2) .AND. (.NOT. l_sulpc_ozone)) THEN
      lmask=.FALSE.
    END IF
  END IF

ELSE IF (isec == 18) THEN
  ! Data assimilation
  IF (sum_iopn /= 0) THEN
    n1 = MOD (iopn(1),10)
    n2 = MOD((iopn(1)/10),10)
  END IF

ELSE IF (isec == 20) THEN
  IF (sum_iopn /= 0) THEN 
    n1 = MOD (iopn(1),10)
    ! Diagnostic not available when convection is turned off 
    IF ((n1 == 1) .AND. (.NOT. l_param_conv)) THEN
      lmask=.FALSE.
    END IF
  END IF

ELSE IF ((isec >= 1) .AND. (isec <= 20)) THEN  ! BUT NOT 1,3,18
  IF (sum_iopn /= 0) THEN
    WRITE(umMessage,*)                                                    &
   'MESSAGE FROM ROUTINE TSTMSK: UNEXPECTED OPTION CODE',         &
    iopn(1)
    CALL umPrint(umMessage,src='tstmsk')
    WRITE(umMessage,*)'IN ATMOSPHERE SECTION ',isec
    CALL umPrint(umMessage,src='tstmsk')
  END IF

ELSE IF (isec == 21) THEN

  ! Thunderstorm electrification
  IF (.NOT. l_use_electric) THEN
    lmask = .FALSE.
  ELSE IF (sum_iopn /= 0) THEN
    n2 = MOD((iopn(1)/10),10)

    IF ((n2 == 1) .AND. (electric_method /= 2)) THEN
      lmask = .FALSE.
    END IF

  END IF

ELSE IF ((isec >= 22) .AND. (isec <= 24)) THEN
  ! Atmos climate mean diagnostics - redundant

ELSE IF (isec  ==  26) THEN
  !  River Routing Diagnostics
  IF (sum_iopn /= 0) THEN
    n22 = MOD((iopn(5)/10),10)
    IF ((n22 == 1) .AND. (.NOT. l_rivers)) THEN
      lmask=.FALSE.    !  River routing switched off
    ELSE IF ((n22 == 2) .AND. (.NOT. l_inland)) THEN
      lmask=.FALSE.    !Inland re-routing switched off
    END IF
  END IF  
ELSE IF (isec == 30) THEN

  ! FV-TRACK diagnostics and CO2 mass scalar
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1), 10)
    n15 =MOD((iopn(3)/10000),10)
    IF ((n1 == 1) .AND. (h_global(atmos_im) /= "Y"  )) THEN
      lmask = .FALSE.
    ELSE IF ((n15 == 3) .AND. (.NOT. l_co2_interactive)) THEN
      lmask = .FALSE.   ! only when i_co2_opt = 3 
    END IF
  END IF
ELSE IF (isec  ==  35) THEN
  ! Stochastic Physics (sect. 35)
  IF (sum_iopn /= 0) THEN
    n1 = MOD(iopn(1),     10)
    IF ((n1 == 1) .AND. (.NOT. l_skeb2)) THEN
      lmask=.FALSE.    !  l_skeb2 switched off
    ELSE IF ((n1 == 2) .AND. (.NOT. l_spt)) THEN
      lmask=.FALSE.    ! l_spt switched off
    END IF
  END IF

ELSE IF (isec == 53) THEN

  ! Idealised UM
  IF (.NOT. l_idealised_data) THEN
    lmask = .FALSE.
  END IF

END IF
! End of atmos diagnostic block

IF (lmask) THEN
  ladres=.TRUE.
  IF ((ispace == 3) .OR. (ispace == 5)) THEN
    lmask=.FALSE.
  END IF
ELSE
  ladres=.FALSE.
END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE tstmsk
