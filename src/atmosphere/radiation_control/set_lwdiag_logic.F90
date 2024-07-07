! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Set the LW radiation diagnostic flags to true if necessary
!
! Method:
!
! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! LW_diag is set to true.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
MODULE set_lwdiag_logic_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_LWDIAG_LOGIC_MOD'
CONTAINS

SUBROUTINE set_lwdiag_logic(sf, nitems, nsects, j_lw, i_off)

USE lw_diag_mod, ONLY: lw_diag
USE rad_input_mod, ONLY: l_rad_perturb, l_radiance, l_forcing
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE model_domain_mod, ONLY: model_type, mt_single_column
USE spec_sw_lw, ONLY: lw_spectrum
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
USE stash_array_mod, ONLY:                                                   &
      len_stlist, stash_pseudo_levels, num_stash_pseudo, stindex, stlist
USE submodel_mod, ONLY: atmos_im

IMPLICIT NONE

! arguments with intent in

INTEGER, INTENT(IN) :: nitems       ! STASH item number
INTEGER, INTENT(IN) :: nsects       ! STASH section number
INTEGER, INTENT(IN) :: j_lw         ! call to LW radiation
INTEGER, INTENT(IN) :: i_off        ! offset for diagnostics
LOGICAL, INTENT(IN) :: sf(0:,0:)    ! STASH Flags

INTEGER :: im_index ! internal model index
INTEGER :: sect    ! STASH section number
INTEGER :: item    ! STASH item number
INTEGER :: icode                    ! Return code  =0 Normal exit  >1 Error
INTEGER :: nps     ! number of pseudo levels for 3D aerosol-radiation diags
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_LWDIAG_LOGIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0 ! Initialise error status

SELECT CASE (model_type)
CASE DEFAULT

  IF (l_radiance .AND. (j_lw > 1)) THEN

    ! TOA radiances
    lw_diag(j_lw)%l_toa_radiance      = sf(297+i_off,2)

  ELSE IF (l_forcing .AND. (j_lw > 1)) THEN

    ! Fluxes and Heating Rates with an applied forcing
    lw_diag(j_lw)%l_flux_up                 = sf(217+i_off,2)
    lw_diag(j_lw)%l_flux_down               = sf(218+i_off,2)
    lw_diag(j_lw)%l_flux_up_clear           = sf(219+i_off,2)
    lw_diag(j_lw)%l_flux_down_clear         = sf(220+i_off,2)
    lw_diag(j_lw)%l_clear_hr                = sf(233+i_off,2)
    lw_diag(j_lw)%l_net_flux_trop           = sf(237+i_off,2)
    lw_diag(j_lw)%l_down_flux_trop          = sf(238+i_off,2)
    lw_diag(j_lw)%l_flux_up_band            = sf(509+i_off,2)
    lw_diag(j_lw)%l_flux_down_band          = sf(510+i_off,2)
    lw_diag(j_lw)%l_flux_up_clear_band      = sf(511+i_off,2)
    lw_diag(j_lw)%l_flux_down_clear_band    = sf(512+i_off,2)
    lw_diag(j_lw)%l_emission_spectrum       = sf(513+i_off,2)      
    lw_diag(j_lw)%l_emission_spectrum_clear = sf(514+i_off,2)      
    lw_diag(j_lw)%l_emission_spectrum_clean       = sf(515+i_off,2)      
    lw_diag(j_lw)%l_emission_spectrum_clear_clean = sf(516+i_off,2)      

    lw_diag(j_lw)%l_diag_call = lw_diag(j_lw)%l_flux_up                   .OR. &
                                lw_diag(j_lw)%l_flux_down                 .OR. &
                                lw_diag(j_lw)%l_flux_up_clear             .OR. &
                                lw_diag(j_lw)%l_flux_down_clear           .OR. &
                                lw_diag(j_lw)%l_clear_hr                  .OR. &
                                lw_diag(j_lw)%l_net_flux_trop             .OR. &
                                lw_diag(j_lw)%l_down_flux_trop            .OR. &
                                lw_diag(j_lw)%l_flux_up_band              .OR. &
                                lw_diag(j_lw)%l_flux_down_band            .OR. &
                                lw_diag(j_lw)%l_flux_up_clear_band        .OR. &
                                lw_diag(j_lw)%l_flux_down_clear_band      .OR. &
                                lw_diag(j_lw)%l_emission_spectrum         .OR. &
                                lw_diag(j_lw)%l_emission_spectrum_clear   .OR. &
                                lw_diag(j_lw)%l_emission_spectrum_clean   .OR. &
                                lw_diag(j_lw)%l_emission_spectrum_clear_clean

  ELSE

    ! Fluxes
    lw_diag(j_lw)%l_flux_up           = sf(217,2)
    lw_diag(j_lw)%l_flux_down         = sf(218,2)
    lw_diag(j_lw)%l_net_flux_trop     = sf(237,2)
    lw_diag(j_lw)%l_down_flux_trop    = sf(238,2)

    ! Cloud amount
    lw_diag(j_lw)%l_total_cloud_cover     = sf(204,2)
    lw_diag(j_lw)%l_total_cloud_on_levels = sf(261,2)

    ! Grid-box mean cloud diagnostics as seen by radiation:
    lw_diag(j_lw)%l_ls_qcl_rad            = sf(308,2)
    lw_diag(j_lw)%l_ls_qcf_rad            = sf(309,2)
    lw_diag(j_lw)%l_cc_qcl_rad            = sf(310,2)
    lw_diag(j_lw)%l_cc_qcf_rad            = sf(311,2)
    lw_diag(j_lw)%l_ls_cl_rad             = sf(312,2)
    lw_diag(j_lw)%l_ls_cf_rad             = sf(313,2)
    lw_diag(j_lw)%l_cc_cl_rad             = sf(314,2)
    lw_diag(j_lw)%l_cc_cf_rad             = sf(315,2)
    ! Cloud effective dimension diagnostics
    lw_diag(j_lw)%l_ls_del_rad            = sf(397,2)
    lw_diag(j_lw)%l_ls_def_rad            = sf(398,2)
    lw_diag(j_lw)%l_cc_del_rad            = sf(399,2)
    lw_diag(j_lw)%l_cc_def_rad            = sf(400,2)

    ! Cloud water paths in radiation
    lw_diag(j_lw)%l_ls_qcl_rad_path       = sf(391,2)
    lw_diag(j_lw)%l_ls_qcf_rad_path       = sf(392,2)
    lw_diag(j_lw)%l_cc_qcl_rad_path       = sf(393,2)
    lw_diag(j_lw)%l_cc_qcf_rad_path       = sf(394,2)

    IF (.NOT. (l_rad_perturb .AND. (j_lw == 2))) THEN
      ! For the incremental time-stepping scheme (l_rad_perturb) many
      ! of the diagnostics are not calculated on the "cloud only"
      ! radiation calls (j_lw==2).

      ! Clear-sky Fluxes and Heating Rates
      lw_diag(j_lw)%l_flux_up_clear     = sf(219,2)
      lw_diag(j_lw)%l_flux_down_clear   = sf(220,2)
      lw_diag(j_lw)%l_toa_clear_weighted= sf(228,2)
      lw_diag(j_lw)%l_total_clear_area  = sf(229,2)
      lw_diag(j_lw)%l_clear_olr         = sf(206,2)
      lw_diag(j_lw)%l_surf_down_clr     = sf(208,2)
      lw_diag(j_lw)%l_clear_hr          = sf(233,2)
              
      ! Extinction and absorptivity diagnostics
      lw_diag(j_lw)%l_cloud_absorptivity            = sf(262,2)
      lw_diag(j_lw)%l_cloud_weight_absorptivity     = sf(263,2)
      lw_diag(j_lw)%l_ls_cloud_absorptivity         = sf(264,2)
      lw_diag(j_lw)%l_ls_cloud_weight_absorptivity  = sf(265,2)
      lw_diag(j_lw)%l_cnv_cloud_absorptivity        = sf(266,2)
      lw_diag(j_lw)%l_cnv_cloud_weight_absorptivity = sf(267,2)

      ! CLASSIC aerosol optical depth diagnostics
      lw_diag(j_lw)%l_aod_sulphate                  = sf(284,2)
      lw_diag(j_lw)%l_aod_dust                      = sf(285,2)
      lw_diag(j_lw)%l_aod_seasalt                   = sf(286,2)
      lw_diag(j_lw)%l_aod_soot                      = sf(287,2)
      lw_diag(j_lw)%l_aod_biomass                   = sf(288,2)
      lw_diag(j_lw)%l_aod_biogenic                  = sf(289,2)
      lw_diag(j_lw)%l_aod_ocff                      = sf(295,2)
      lw_diag(j_lw)%l_aod_delta                     = sf(296,2)
      lw_diag(j_lw)%l_aod_nitrate                   = sf(297,2)
      lw_diag(j_lw)%l_aod_total_radn                = sf(298,2)
      lw_diag(j_lw)%l_angst_total_radn              = sf(299,2)
      lw_diag(j_lw)%l_aod_prog_sulphate             = sf(421,2)
      lw_diag(j_lw)%l_aod_prog_dust                 = sf(422,2)
      lw_diag(j_lw)%l_aod_prog_seasalt              = sf(423,2)
      lw_diag(j_lw)%l_aod_prog_soot                 = sf(424,2)
      lw_diag(j_lw)%l_aod_prog_biomass              = sf(425,2)
      lw_diag(j_lw)%l_aod_prog_ocff                 = sf(426,2)
      lw_diag(j_lw)%l_aod_prog_nitrate              = sf(427,2)
              
      ! CLASSIC aerosol absorption optical depth diagnostics
      lw_diag(j_lw)%l_aaod_sulphate                 = sf(584,2)
      lw_diag(j_lw)%l_aaod_dust                     = sf(585,2)
      lw_diag(j_lw)%l_aaod_seasalt                  = sf(586,2)
      lw_diag(j_lw)%l_aaod_soot                     = sf(587,2)
      lw_diag(j_lw)%l_aaod_biomass                  = sf(588,2)
      lw_diag(j_lw)%l_aaod_biogenic                 = sf(589,2)
      lw_diag(j_lw)%l_aaod_ocff                     = sf(595,2)
      lw_diag(j_lw)%l_aaod_nitrate                  = sf(597,2)

      ! UKCA aerosol optical depth diagnostics
      lw_diag(j_lw)%l_aod_ukca_ait_sol              = sf(300,2)
      lw_diag(j_lw)%l_aod_ukca_acc_sol              = sf(301,2)
      lw_diag(j_lw)%l_aod_ukca_cor_sol              = sf(302,2)
      lw_diag(j_lw)%l_aod_ukca_ait_ins              = sf(303,2)
      lw_diag(j_lw)%l_aod_ukca_acc_ins              = sf(304,2)
      lw_diag(j_lw)%l_aod_ukca_cor_ins              = sf(305,2)

      ! UKCA stratospheric aerosol optical depth diagnostics
      lw_diag(j_lw)%l_sod_ukca_ait_sol              = sf(251,2)
      lw_diag(j_lw)%l_sod_ukca_acc_sol              = sf(252,2)
      lw_diag(j_lw)%l_sod_ukca_cor_sol              = sf(253,2)
      lw_diag(j_lw)%l_sod_ukca_ait_ins              = sf(254,2)
      lw_diag(j_lw)%l_sod_ukca_acc_ins              = sf(255,2)
      lw_diag(j_lw)%l_sod_ukca_cor_ins              = sf(256,2)

      ! UKCA absorption aerosol optical depth diagnostics
      lw_diag(j_lw)%l_aaod_ukca_ait_sol              = sf(240,2)
      lw_diag(j_lw)%l_aaod_ukca_acc_sol              = sf(241,2)
      lw_diag(j_lw)%l_aaod_ukca_cor_sol              = sf(242,2)
      lw_diag(j_lw)%l_aaod_ukca_ait_ins              = sf(243,2)
      lw_diag(j_lw)%l_aaod_ukca_acc_ins              = sf(244,2)
      lw_diag(j_lw)%l_aaod_ukca_cor_ins              = sf(245,2)

      ! 3D aerosol-radiation diagnostics are defined for AOD wavelengths
      ! using pseudo levels for the wavelength dimension
      !
      nps      = lw_spectrum(j_lw)%aerosol%n_aod_wavel  
      im_index = 1
      sect     = 2
      
      ! CLASSIC 3D aerosol extinction
      item = 540
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_clas_aerosol_ext,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! CLASSIC 3D aerosol absorption
      item = 541
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_clas_aerosol_abs,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! CLASSIC 3D aerosol scattering
      item = 542
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_clas_aerosol_sca,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! UKCA 3D aerosol extinction
      item = 530
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_ukca_aerosol_ext,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! UKCA 3D aerosol absorption
      item = 531
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_ukca_aerosol_abs,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! UKCA 3D aerosol scattering
      item = 532
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_ukca_aerosol_sca,                     &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      ! UKCA 3D aerosol asymmetry * scattering
      item = 533
      IF (sf(item,sect)) THEN
        ! DEPENDS ON: set_pseudo_list
        CALL set_pseudo_list(nps, len_stlist,                      &
             stlist(1,stindex(1,item,sect,im_index)),              &
             lw_diag(j_lw)%l_ukca_aerosol_gsca,                    &
             stash_pseudo_levels,num_stash_pseudo,icode,cmessage)
        IF (icode >  0) THEN
          WRITE(cmessage,'(A,I3,2A)')                              &
               'Error in set_pseudo_list (item ', item, '): ',     &
               cmessage(1:LEN_TRIM(cmessage))
          CALL Ereport(RoutineName,icode,cmessage)
        END IF
      END IF

      lw_diag(j_lw)%n_clas_aerosol_ext = COUNT(lw_diag(j_lw)%l_clas_aerosol_ext)
      lw_diag(j_lw)%n_clas_aerosol_abs = COUNT(lw_diag(j_lw)%l_clas_aerosol_abs)
      lw_diag(j_lw)%n_clas_aerosol_sca = COUNT(lw_diag(j_lw)%l_clas_aerosol_sca)
      lw_diag(j_lw)%n_ukca_aerosol_ext = COUNT(lw_diag(j_lw)%l_ukca_aerosol_ext)
      lw_diag(j_lw)%n_ukca_aerosol_abs = COUNT(lw_diag(j_lw)%l_ukca_aerosol_abs)
      lw_diag(j_lw)%n_ukca_aerosol_sca = COUNT(lw_diag(j_lw)%l_ukca_aerosol_sca)
      lw_diag(j_lw)%n_ukca_aerosol_gsca=COUNT(lw_diag(j_lw)%l_ukca_aerosol_gsca)

      ! Total aerosol Optical Properties
      lw_diag(j_lw)%l_aerosol_optical_depth         = sf(506,2)
      lw_diag(j_lw)%l_aerosol_scat_optical_depth    = sf(507,2)
      lw_diag(j_lw)%l_aerosol_asymmetry_scat        = sf(508,2)

      ! Band-by-band Fluxes
      lw_diag(j_lw)%l_flux_up_band                  = sf(509,2)
      lw_diag(j_lw)%l_flux_down_band                = sf(510,2)
      lw_diag(j_lw)%l_flux_up_clear_band            = sf(511,2)
      lw_diag(j_lw)%l_flux_down_clear_band          = sf(512,2)      
      lw_diag(j_lw)%l_emission_spectrum             = sf(513,2)      
      lw_diag(j_lw)%l_emission_spectrum_clear       = sf(514,2)      

      ! Clean-air Fluxes
      lw_diag(j_lw)%l_flux_up_clean                 = sf(517,2)
      lw_diag(j_lw)%l_flux_down_clean               = sf(518,2)
      lw_diag(j_lw)%l_flux_up_clear_clean           = sf(519,2)
      lw_diag(j_lw)%l_flux_down_clear_clean         = sf(520,2)
      lw_diag(j_lw)%l_flux_up_clean_band            = sf(545,2)
      lw_diag(j_lw)%l_flux_down_clean_band          = sf(546,2)
      lw_diag(j_lw)%l_flux_up_clear_clean_band      = sf(547,2)
      lw_diag(j_lw)%l_flux_down_clear_clean_band    = sf(548,2)
      lw_diag(j_lw)%l_emission_spectrum_clean       = sf(515,2)      
      lw_diag(j_lw)%l_emission_spectrum_clear_clean = sf(516,2)      

      ! Fluxes after green-house gas forcing
      lw_diag(j_lw)%l_flux_up_forc                  = sf(521,2)
      lw_diag(j_lw)%l_flux_down_forc                = sf(522,2)
      lw_diag(j_lw)%l_flux_up_clear_forc            = sf(523,2)
      lw_diag(j_lw)%l_flux_down_clear_forc          = sf(524,2)
      lw_diag(j_lw)%l_flux_up_forc_band             = sf(525,2)
      lw_diag(j_lw)%l_flux_down_forc_band           = sf(526,2)
      lw_diag(j_lw)%l_flux_up_clear_forc_band       = sf(527,2)
      lw_diag(j_lw)%l_flux_down_clear_forc_band     = sf(528,2)

      ! EasyAerosol
      lw_diag(j_lw)%l_easyaerosol_extinction        = sf(550,2)
      lw_diag(j_lw)%l_easyaerosol_absorption        = sf(551,2)
      lw_diag(j_lw)%l_easyaerosol_scattering        = sf(552,2)
      lw_diag(j_lw)%l_easyaerosol_asytimscat        = sf(553,2)

      ! Convective core diagnotics
      lw_diag(j_lw)%l_ccore_clt_rad                 = sf(317,2)
      lw_diag(j_lw)%l_ccore_qcl_rad                 = sf(318,2)
      lw_diag(j_lw)%l_ccore_qcf_rad                 = sf(319,2)
      lw_diag(j_lw)%l_ccore_qcl_rad_path            = sf(395,2)
      lw_diag(j_lw)%l_ccore_qcf_rad_path            = sf(396,2)

    END IF ! .not.(l_rad_perturb.and.(j_lw == 2))

  END IF

CASE (mt_single_column)

  lw_diag(j_lw)%l_total_cloud_cover      = .TRUE.
  lw_diag(j_lw)%l_ls_qcl_rad             = .TRUE.
  lw_diag(j_lw)%l_ls_qcf_rad             = .TRUE.
  lw_diag(j_lw)%l_cc_qcl_rad             = .TRUE.
  lw_diag(j_lw)%l_cc_qcf_rad             = .TRUE.
  lw_diag(j_lw)%l_ls_cl_rad              = .TRUE.
  lw_diag(j_lw)%l_ls_cf_rad              = .TRUE.
  lw_diag(j_lw)%l_cc_cl_rad              = .TRUE.
  lw_diag(j_lw)%l_cc_cf_rad              = .TRUE.

  IF (.NOT. (l_rad_perturb .AND. (j_lw == 2))) THEN
    ! For the incremental time-stepping scheme (l_rad_perturb) many
    ! of the diagnostics are not calculated on the "cloud only"
    ! radiation calls (j_lw==2).
    lw_diag(j_lw)%l_clear_olr              = .TRUE.
    lw_diag(j_lw)%l_surf_down_clr          = .TRUE.
    lw_diag(j_lw)%l_clear_hr               = .TRUE.
  END IF

END SELECT ! model_type


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_lwdiag_logic





END MODULE set_lwdiag_logic_mod
