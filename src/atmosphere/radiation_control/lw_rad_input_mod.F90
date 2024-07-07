! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for lw radiation.

! Description:
!   Module containing control as used by the lw radiation code.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE lw_rad_input_mod

USE sat_opt_mod
USE lw_control_struct, ONLY: lw_control
USE rad_input_mod, ONLY: n_lwcall
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER  (LEN=132) :: spectral_file_lw    ! Spectral file

! Options for gases:
LOGICAL :: l_h2o_lw            ! flag for water vapour
LOGICAL :: l_co2_lw            ! flag for carbon dioxide
LOGICAL :: l_o3_lw             ! flag for ozone
LOGICAL :: l_n2o_lw            ! flag for nitrous oxide
LOGICAL :: l_ch4_lw            ! flag for methane
LOGICAL :: l_so2_lw            ! flag for sulphur dioxide
LOGICAL :: l_cfc11_lw          ! flag for cfc11
LOGICAL :: l_cfc12_lw          ! flag for cfc12
LOGICAL :: l_cfc113_lw         ! flag for cfc113
LOGICAL :: l_cfc114_lw         ! flag for cfc114
LOGICAL :: l_hcfc22_lw         ! flag for hcfc22
LOGICAL :: l_hfc125_lw         ! flag for hfc125
LOGICAL :: l_hfc134a_lw        ! flag for hfc134a
LOGICAL :: l_co_lw             ! flag for carbon monoxide
LOGICAL :: l_nh3_lw            ! flag for ammonia
LOGICAL :: l_tio_lw            ! flag for titanium oxide
LOGICAL :: l_vo_lw             ! flag for vanadium oxide
LOGICAL :: l_h2_lw             ! flag for H2-H2 CIA
LOGICAL :: l_he_lw             ! flag for H2-He CIA
LOGICAL :: l_na_lw             ! flag for sodium
LOGICAL :: l_k_lw              ! flag for potassium
LOGICAL :: l_li_lw             ! flag for lithium
LOGICAL :: l_rb_lw             ! flag for rubidium
LOGICAL :: l_cs_lw             ! flag for cesium

LOGICAL :: l_solar_tail_flux   ! Flag for adding solar tail flux to LW
LOGICAL :: l_microphysics_lw   ! Flag for microphysics in LW

INTEGER :: i_scatter_method_lw ! treatment of scattering
INTEGER :: i_gas_overlap_lw    ! treatment of gaseous overlaps

! types of droplets or ice crystals used for parametrisations
INTEGER :: i_st_water_lw             ! type for stratiform water
INTEGER :: i_cnv_water_lw            ! type for convective water
INTEGER :: i_st_ice_lw               ! type for stratiform ice
INTEGER :: i_cnv_ice_lw              ! type for convective ice

! Second set of control options for the second radiation call:
CHARACTER  (LEN=132) :: spectral_file_lw2
LOGICAL :: l_h2o_lw2
LOGICAL :: l_co2_lw2
LOGICAL :: l_o3_lw2
LOGICAL :: l_n2o_lw2
LOGICAL :: l_ch4_lw2
LOGICAL :: l_so2_lw2
LOGICAL :: l_cfc11_lw2
LOGICAL :: l_cfc12_lw2
LOGICAL :: l_cfc113_lw2
LOGICAL :: l_cfc114_lw2
LOGICAL :: l_hcfc22_lw2
LOGICAL :: l_hfc125_lw2
LOGICAL :: l_hfc134a_lw2
LOGICAL :: l_co_lw2
LOGICAL :: l_nh3_lw2
LOGICAL :: l_tio_lw2
LOGICAL :: l_vo_lw2
LOGICAL :: l_h2_lw2
LOGICAL :: l_he_lw2
LOGICAL :: l_na_lw2
LOGICAL :: l_k_lw2
LOGICAL :: l_li_lw2
LOGICAL :: l_rb_lw2
LOGICAL :: l_cs_lw2
LOGICAL :: l_solar_tail_flux_2
LOGICAL :: l_microphysics_lw2
INTEGER :: i_scatter_method_lw2
INTEGER :: i_gas_overlap_lw2
INTEGER :: i_st_water_lw2
INTEGER :: i_cnv_water_lw2
INTEGER :: i_st_ice_lw2
INTEGER :: i_cnv_ice_lw2

! --------------------

NAMELIST/r2lwclnl/                                                      &
     spectral_file_lw, spectral_file_lw2,                               &
     l_solar_tail_flux, l_solar_tail_flux_2,                            &
     i_gas_overlap_lw, i_gas_overlap_lw2,                               &
     l_h2o_lw, l_h2o_lw2, l_co2_lw, l_co2_lw2, l_o3_lw, l_o3_lw2,       &
     l_n2o_lw, l_n2o_lw2, l_ch4_lw, l_ch4_lw2, l_so2_lw, l_so2_lw2,     &
     l_cfc11_lw, l_cfc11_lw2, l_cfc12_lw, l_cfc12_lw2,                  &
     l_cfc113_lw, l_cfc113_lw2, l_cfc114_lw, l_cfc114_lw2,              &
     l_hcfc22_lw, l_hcfc22_lw2, l_hfc125_lw, l_hfc125_lw2,              &
     l_hfc134a_lw, l_hfc134a_lw2,                                       &
     l_co_lw, l_co_lw2, l_nh3_lw, l_nh3_lw2, l_tio_lw, l_tio_lw2,       &
     l_vo_lw, l_vo_lw2, l_h2_lw, l_h2_lw2, l_he_lw, l_he_lw2,           &
     l_na_lw, l_na_lw2, l_k_lw, l_k_lw2,                                &
     l_li_lw, l_li_lw2, l_rb_lw, l_rb_lw2, l_cs_lw, l_cs_lw2,           &
     i_st_water_lw, i_st_water_lw2, i_cnv_water_lw, i_cnv_water_lw2,    &
     i_st_ice_lw, i_st_ice_lw2, i_cnv_ice_lw, i_cnv_ice_lw2,            &
     i_scatter_method_lw, i_scatter_method_lw2,                         &
     l_microphysics_lw, l_microphysics_lw2

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ----------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LW_RAD_INPUT_MOD'

CONTAINS

! Subroutine to set the input values of the lw control structure.

SUBROUTINE lw_input

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica, ip_cloud_clear,                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan, ip_solver_triple_app_scat,                    &
  ip_solver_mix_app_scat, ip_solver_no_scat,                            &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand,                            &
  ip_scatter_approx, ip_no_scatter_abs, ip_no_scatter_ext
USE ereport_mod, ONLY: ereport
USE rad_input_mod
USE file_manager, ONLY: get_file_unit_by_id

USE model_domain_mod, ONLY: model_type, mt_single_column
USE lw_control_default_mod, ONLY: lw_control_default

IMPLICIT NONE

INTEGER :: j,                        &
     errorstatus      ! Return code : 0 Normal Exit : >0 Error
INTEGER :: atmoscntl_unit

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='LW_INPUT')
CHARACTER(LEN=errormessagelength) :: CMessage  ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialisation: set namelist items to MDI / FALSE.
spectral_file_lw           = ''
l_solar_tail_flux          = .FALSE.
i_gas_overlap_lw           = imdi
l_h2o_lw                   = .FALSE.
l_co2_lw                   = .FALSE.
l_o3_lw                    = .FALSE.
l_n2o_lw                   = .FALSE.
l_ch4_lw                   = .FALSE.
l_so2_lw                   = .FALSE.
l_cfc11_lw                 = .FALSE.
l_cfc12_lw                 = .FALSE.
l_cfc113_lw                = .FALSE.
l_cfc114_lw                = .FALSE.
l_hcfc22_lw                = .FALSE.
l_hfc125_lw                = .FALSE.
l_hfc134a_lw               = .FALSE.
l_co_lw                    = .FALSE.
l_nh3_lw                   = .FALSE.
l_tio_lw                   = .FALSE.
l_vo_lw                    = .FALSE.
l_h2_lw                    = .FALSE.
l_he_lw                    = .FALSE.
l_na_lw                    = .FALSE.
l_k_lw                     = .FALSE.
l_li_lw                    = .FALSE.
l_rb_lw                    = .FALSE.
l_cs_lw                    = .FALSE.
i_st_water_lw              = imdi
i_cnv_water_lw             = imdi
i_st_ice_lw                = imdi
i_cnv_ice_lw               = imdi
i_scatter_method_lw        = imdi
l_microphysics_lw          = .FALSE.

spectral_file_lw2          = ''
l_solar_tail_flux_2        = .FALSE.
i_gas_overlap_lw2          = imdi
l_h2o_lw2                  = .FALSE.
l_co2_lw2                  = .FALSE.
l_o3_lw2                   = .FALSE.
l_n2o_lw2                  = .FALSE.
l_ch4_lw2                  = .FALSE.
l_so2_lw2                  = .FALSE.
l_cfc11_lw2                = .FALSE.
l_cfc12_lw2                = .FALSE.
l_cfc113_lw2               = .FALSE.
l_cfc114_lw2               = .FALSE.
l_hcfc22_lw2               = .FALSE.
l_hfc125_lw2               = .FALSE.
l_hfc134a_lw2              = .FALSE.
l_co_lw2                   = .FALSE.
l_nh3_lw2                  = .FALSE.
l_tio_lw2                  = .FALSE.
l_vo_lw2                   = .FALSE.
l_h2_lw2                   = .FALSE.
l_he_lw2                   = .FALSE.
l_na_lw2                   = .FALSE.
l_k_lw2                    = .FALSE.
l_li_lw2                   = .FALSE.
l_rb_lw2                   = .FALSE.
l_cs_lw2                   = .FALSE.
i_st_water_lw2             = imdi
i_cnv_water_lw2            = imdi
i_st_ice_lw2               = imdi
i_cnv_ice_lw2              = imdi
i_scatter_method_lw2       = imdi
l_microphysics_lw2         = .FALSE.

! Options for the longwave
SELECT CASE (model_type)

CASE (mt_single_column)
  READ (10, r2lwclnl)

CASE DEFAULT
  atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
  CALL read_nml_r2lwclnl(atmoscntl_unit)

END SELECT

DO j=1, n_lwcall
  ! Set default values of control variables.
  ! Note: to allow a simple interface for ROSE only the commonly used
  ! items are now set via the namelist. In order to set further variables
  ! (eg. for radiance diagnostics) a branch must now be used.

  CALL lw_control_default(lw_control(j))

  ! Transfer namelist items to the data structure.
  IF (j==2) THEN
    lw_control(j)%spectral_file          = ADJUSTL(spectral_file_lw2)
    lw_control(j)%l_solar_tail_flux      = l_solar_tail_flux_2
    lw_control(j)%i_gas_overlap          = i_gas_overlap_lw2
    lw_control(j)%i_cloud_representation = i_cloud_representation_2
    lw_control(j)%i_overlap              = i_overlap_2
    lw_control(j)%i_inhom                = i_inhom_2
    lw_control(j)%l_continuum            = l_h2o_lw2
    lw_control(j)%l_h2o                  = l_h2o_lw2
    lw_control(j)%l_co2                  = l_co2_lw2
    lw_control(j)%l_o3                   = l_o3_lw2
    lw_control(j)%l_n2o                  = l_n2o_lw2
    lw_control(j)%l_ch4                  = l_ch4_lw2
    lw_control(j)%l_so2                  = l_so2_lw2
    lw_control(j)%l_cfc11                = l_cfc11_lw2
    lw_control(j)%l_cfc12                = l_cfc12_lw2
    lw_control(j)%l_cfc113               = l_cfc113_lw2
    lw_control(j)%l_cfc114               = l_cfc114_lw2
    lw_control(j)%l_hcfc22               = l_hcfc22_lw2
    lw_control(j)%l_hfc125               = l_hfc125_lw2
    lw_control(j)%l_hfc134a              = l_hfc134a_lw2
    lw_control(j)%l_co                   = l_co_lw2
    lw_control(j)%l_nh3                  = l_nh3_lw2
    lw_control(j)%l_tio                  = l_tio_lw2
    lw_control(j)%l_vo                   = l_vo_lw2
    lw_control(j)%l_h2                   = l_h2_lw2
    lw_control(j)%l_he                   = l_he_lw2
    lw_control(j)%l_na                   = l_na_lw2
    lw_control(j)%l_k                    = l_k_lw2
    lw_control(j)%l_li                   = l_li_lw2
    lw_control(j)%l_rb                   = l_rb_lw2
    lw_control(j)%l_cs                   = l_cs_lw2
    lw_control(j)%i_st_water             = i_st_water_lw2
    lw_control(j)%i_cnv_water            = i_cnv_water_lw2
    lw_control(j)%i_st_ice               = i_st_ice_lw2
    lw_control(j)%i_cnv_ice              = i_cnv_ice_lw2
    lw_control(j)%i_scatter_method       = i_scatter_method_lw2
    lw_control(j)%l_microphysics         = l_microphysics_lw2
  ELSE
    lw_control(j)%spectral_file          = ADJUSTL(spectral_file_lw)
    lw_control(j)%l_solar_tail_flux      = l_solar_tail_flux
    lw_control(j)%i_gas_overlap          = i_gas_overlap_lw
    lw_control(j)%i_cloud_representation = i_cloud_representation
    lw_control(j)%i_overlap              = i_overlap
    lw_control(j)%i_inhom                = i_inhom
    lw_control(j)%l_continuum            = l_h2o_lw
    lw_control(j)%l_h2o                  = l_h2o_lw
    lw_control(j)%l_co2                  = l_co2_lw
    lw_control(j)%l_o3                   = l_o3_lw
    lw_control(j)%l_n2o                  = l_n2o_lw
    lw_control(j)%l_ch4                  = l_ch4_lw
    lw_control(j)%l_so2                  = l_so2_lw
    lw_control(j)%l_cfc11                = l_cfc11_lw
    lw_control(j)%l_cfc12                = l_cfc12_lw
    lw_control(j)%l_cfc113               = l_cfc113_lw
    lw_control(j)%l_cfc114               = l_cfc114_lw
    lw_control(j)%l_hcfc22               = l_hcfc22_lw
    lw_control(j)%l_hfc125               = l_hfc125_lw
    lw_control(j)%l_hfc134a              = l_hfc134a_lw
    lw_control(j)%l_co                   = l_co_lw
    lw_control(j)%l_nh3                  = l_nh3_lw
    lw_control(j)%l_tio                  = l_tio_lw
    lw_control(j)%l_vo                   = l_vo_lw
    lw_control(j)%l_h2                   = l_h2_lw
    lw_control(j)%l_he                   = l_he_lw
    lw_control(j)%l_na                   = l_na_lw
    lw_control(j)%l_k                    = l_k_lw
    lw_control(j)%l_li                   = l_li_lw
    lw_control(j)%l_rb                   = l_rb_lw
    lw_control(j)%l_cs                   = l_cs_lw
    lw_control(j)%i_st_water             = i_st_water_lw
    lw_control(j)%i_cnv_water            = i_cnv_water_lw
    lw_control(j)%i_st_ice               = i_st_ice_lw
    lw_control(j)%i_cnv_ice              = i_cnv_ice_lw
    lw_control(j)%i_scatter_method       = i_scatter_method_lw
    lw_control(j)%l_microphysics         = l_microphysics_lw
  END IF
END DO

! Set i_cloud from combination of i_overlap and i_representation
DO j=1,n_lwcall
  IF (lw_control(j)%i_cloud_representation == ip_cloud_homogen) THEN

    lw_control(j)%i_solver=ip_solver_homogen_direct
    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      lw_control(j)%i_cloud=ip_cloud_mcica
    ELSE
      IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF (lw_control(j)%i_cloud_representation == ip_cloud_ice_water) THEN

    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      lw_control(j)%i_cloud=ip_cloud_mcica
      IF ( (lw_control(j)%i_scatter_method == ip_no_scatter_abs) .OR.   &
           (lw_control(j)%i_scatter_method == ip_no_scatter_ext) ) THEN
        lw_control(j)%i_solver      =ip_solver_no_scat
        lw_control(j)%i_solver_clear=ip_solver_no_scat
      ELSE
        lw_control(j)%i_solver=ip_solver_homogen_direct
      END IF
    ELSE
      IF (lw_control(j)%i_scatter_method == ip_scatter_approx) THEN
        lw_control(j)%i_solver=ip_solver_mix_app_scat
      ELSE
        lw_control(j)%i_solver=ip_solver_mix_direct_hogan
      END IF
      IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF ((lw_control(j)%i_cloud_representation == ip_cloud_conv_strat) .OR.  &
           (lw_control(j)%i_cloud_representation == ip_cloud_csiw)) THEN

    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//                 &
                 ' cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage)
    ELSE
      IF (lw_control(j)%i_scatter_method == ip_scatter_approx) THEN
        lw_control(j)%i_solver=ip_solver_triple_app_scat
      ELSE
        lw_control(j)%i_solver=ip_solver_triple_hogan
      END IF
      IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_triple
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE

    ! No treatment of cloud for LW radiation
    lw_control(j)%l_cloud=.FALSE.
    lw_control(j)%i_cloud=ip_cloud_clear
    IF ( (lw_control(j)%i_scatter_method == ip_no_scatter_abs) .OR.            &
         (lw_control(j)%i_scatter_method == ip_no_scatter_ext) ) THEN
      lw_control(j)%i_solver      =ip_solver_no_scat
      lw_control(j)%i_solver_clear=ip_solver_no_scat
    ELSE
      lw_control(j)%i_solver=ip_solver_homogen_direct
    END IF
    lw_control(j)%l_microphysics=.FALSE.
  END IF

  IF (lw_control(j)%i_inhom == ip_scaling) THEN
    l_inhom_cloud=.TRUE.
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE lw_input

SUBROUTINE read_nml_r2lwclnl(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type


IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_R2LWCLNL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 12
INTEGER, PARAMETER :: n_log = 52
INTEGER, PARAMETER :: n_chars = 2 * 132

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_gas_overlap_lw
  INTEGER :: i_gas_overlap_lw2
  INTEGER :: i_st_water_lw
  INTEGER :: i_st_water_lw2
  INTEGER :: i_cnv_water_lw
  INTEGER :: i_cnv_water_lw2
  INTEGER :: i_st_ice_lw
  INTEGER :: i_st_ice_lw2
  INTEGER :: i_cnv_ice_lw
  INTEGER :: i_cnv_ice_lw2
  INTEGER :: i_scatter_method_lw
  INTEGER :: i_scatter_method_lw2
  LOGICAL :: l_solar_tail_flux
  LOGICAL :: l_solar_tail_flux_2
  LOGICAL :: l_h2o_lw
  LOGICAL :: l_h2o_lw2
  LOGICAL :: l_co2_lw
  LOGICAL :: l_co2_lw2
  LOGICAL :: l_o3_lw
  LOGICAL :: l_o3_lw2
  LOGICAL :: l_n2o_lw
  LOGICAL :: l_n2o_lw2
  LOGICAL :: l_ch4_lw
  LOGICAL :: l_ch4_lw2
  LOGICAL :: l_so2_lw
  LOGICAL :: l_so2_lw2
  LOGICAL :: l_cfc11_lw
  LOGICAL :: l_cfc11_lw2
  LOGICAL :: l_cfc12_lw
  LOGICAL :: l_cfc12_lw2
  LOGICAL :: l_cfc113_lw
  LOGICAL :: l_cfc113_lw2
  LOGICAL :: l_cfc114_lw
  LOGICAL :: l_cfc114_lw2
  LOGICAL :: l_hcfc22_lw
  LOGICAL :: l_hcfc22_lw2
  LOGICAL :: l_hfc125_lw
  LOGICAL :: l_hfc125_lw2
  LOGICAL :: l_hfc134a_lw
  LOGICAL :: l_hfc134a_lw2
  LOGICAL :: l_co_lw
  LOGICAL :: l_co_lw2
  LOGICAL :: l_nh3_lw
  LOGICAL :: l_nh3_lw2
  LOGICAL :: l_tio_lw
  LOGICAL :: l_tio_lw2
  LOGICAL :: l_vo_lw
  LOGICAL :: l_vo_lw2
  LOGICAL :: l_h2_lw
  LOGICAL :: l_h2_lw2
  LOGICAL :: l_he_lw
  LOGICAL :: l_he_lw2
  LOGICAL :: l_na_lw
  LOGICAL :: l_na_lw2
  LOGICAL :: l_k_lw
  LOGICAL :: l_k_lw2
  LOGICAL :: l_li_lw
  LOGICAL :: l_li_lw2
  LOGICAL :: l_rb_lw
  LOGICAL :: l_rb_lw2
  LOGICAL :: l_cs_lw
  LOGICAL :: l_cs_lw2
  LOGICAL :: l_microphysics_lw
  LOGICAL :: l_microphysics_lw2
  CHARACTER (LEN=132) :: spectral_file_lw
  CHARACTER (LEN=132) :: spectral_file_lw2
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=r2lwclnl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist R2LWCLNL", iomessage)

  my_nml % i_gas_overlap_lw     = i_gas_overlap_lw
  my_nml % i_gas_overlap_lw2    = i_gas_overlap_lw2
  my_nml % i_st_water_lw        = i_st_water_lw
  my_nml % i_st_water_lw2       = i_st_water_lw2
  my_nml % i_cnv_water_lw       = i_cnv_water_lw
  my_nml % i_cnv_water_lw2      = i_cnv_water_lw2
  my_nml % i_st_ice_lw          = i_st_ice_lw
  my_nml % i_st_ice_lw2         = i_st_ice_lw2
  my_nml % i_cnv_ice_lw         = i_cnv_ice_lw
  my_nml % i_cnv_ice_lw2        = i_cnv_ice_lw2
  my_nml % i_scatter_method_lw  = i_scatter_method_lw
  my_nml % i_scatter_method_lw2 = i_scatter_method_lw2
  ! end of integers
  my_nml % l_solar_tail_flux   = l_solar_tail_flux
  my_nml % l_solar_tail_flux_2 = l_solar_tail_flux_2
  my_nml % l_h2o_lw      = l_h2o_lw
  my_nml % l_h2o_lw2     = l_h2o_lw2
  my_nml % l_co2_lw      = l_co2_lw
  my_nml % l_co2_lw2     = l_co2_lw2
  my_nml % l_o3_lw       = l_o3_lw
  my_nml % l_o3_lw2      = l_o3_lw2
  my_nml % l_n2o_lw      = l_n2o_lw
  my_nml % l_n2o_lw2     = l_n2o_lw2
  my_nml % l_ch4_lw      = l_ch4_lw
  my_nml % l_ch4_lw2     = l_ch4_lw2
  my_nml % l_so2_lw      = l_so2_lw
  my_nml % l_so2_lw2     = l_so2_lw2
  my_nml % l_cfc11_lw    = l_cfc11_lw
  my_nml % l_cfc11_lw2   = l_cfc11_lw2
  my_nml % l_cfc12_lw    = l_cfc12_lw
  my_nml % l_cfc12_lw2   = l_cfc12_lw2
  my_nml % l_cfc113_lw   = l_cfc113_lw
  my_nml % l_cfc113_lw2  = l_cfc113_lw2
  my_nml % l_cfc114_lw   = l_cfc114_lw
  my_nml % l_cfc114_lw2  = l_cfc114_lw2
  my_nml % l_hcfc22_lw   = l_hcfc22_lw
  my_nml % l_hcfc22_lw2  = l_hcfc22_lw2
  my_nml % l_hfc125_lw   = l_hfc125_lw
  my_nml % l_hfc125_lw2  = l_hfc125_lw2
  my_nml % l_hfc134a_lw  = l_hfc134a_lw
  my_nml % l_hfc134a_lw2 = l_hfc134a_lw2
  my_nml % l_co_lw       = l_co_lw
  my_nml % l_co_lw2      = l_co_lw2
  my_nml % l_nh3_lw      = l_nh3_lw
  my_nml % l_nh3_lw2     = l_nh3_lw2
  my_nml % l_tio_lw      = l_tio_lw
  my_nml % l_tio_lw2     = l_tio_lw2
  my_nml % l_vo_lw       = l_vo_lw
  my_nml % l_vo_lw2      = l_vo_lw2
  my_nml % l_h2_lw       = l_h2_lw
  my_nml % l_h2_lw2      = l_h2_lw2
  my_nml % l_he_lw       = l_he_lw
  my_nml % l_he_lw2      = l_he_lw2
  my_nml % l_na_lw       = l_na_lw
  my_nml % l_na_lw2      = l_na_lw2
  my_nml % l_k_lw        = l_k_lw
  my_nml % l_k_lw2       = l_k_lw2
  my_nml % l_li_lw       = l_li_lw
  my_nml % l_li_lw2      = l_li_lw2
  my_nml % l_rb_lw       = l_rb_lw
  my_nml % l_rb_lw2      = l_rb_lw2
  my_nml % l_cs_lw       = l_cs_lw
  my_nml % l_cs_lw2      = l_cs_lw2
  my_nml % l_microphysics_lw  = l_microphysics_lw
  my_nml % l_microphysics_lw2 = l_microphysics_lw2
  ! end of logicals
  my_nml % spectral_file_lw  = spectral_file_lw
  my_nml % spectral_file_lw2 = spectral_file_lw2

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_gas_overlap_lw  = my_nml % i_gas_overlap_lw
  i_gas_overlap_lw2 = my_nml % i_gas_overlap_lw2
  i_st_water_lw     = my_nml % i_st_water_lw
  i_st_water_lw2    = my_nml % i_st_water_lw2
  i_cnv_water_lw    = my_nml % i_cnv_water_lw
  i_cnv_water_lw2   = my_nml % i_cnv_water_lw2
  i_st_ice_lw       = my_nml % i_st_ice_lw
  i_st_ice_lw2      = my_nml % i_st_ice_lw2
  i_cnv_ice_lw      = my_nml % i_cnv_ice_lw
  i_cnv_ice_lw2     = my_nml % i_cnv_ice_lw2
  i_scatter_method_lw  = my_nml % i_scatter_method_lw
  i_scatter_method_lw2 = my_nml % i_scatter_method_lw2
  ! end of integers
  l_solar_tail_flux   = my_nml % l_solar_tail_flux
  l_solar_tail_flux_2 = my_nml % l_solar_tail_flux_2
  l_h2o_lw      = my_nml % l_h2o_lw
  l_h2o_lw2     = my_nml % l_h2o_lw2
  l_co2_lw      = my_nml % l_co2_lw
  l_co2_lw2     = my_nml % l_co2_lw2
  l_o3_lw       = my_nml % l_o3_lw
  l_o3_lw2      = my_nml % l_o3_lw2
  l_n2o_lw      = my_nml % l_n2o_lw
  l_n2o_lw2     = my_nml % l_n2o_lw2
  l_ch4_lw      = my_nml % l_ch4_lw
  l_ch4_lw2     = my_nml % l_ch4_lw2
  l_so2_lw      = my_nml % l_so2_lw
  l_so2_lw2     = my_nml % l_so2_lw2
  l_cfc11_lw    = my_nml % l_cfc11_lw
  l_cfc11_lw2   = my_nml % l_cfc11_lw2
  l_cfc12_lw    = my_nml % l_cfc12_lw
  l_cfc12_lw2   = my_nml % l_cfc12_lw2
  l_cfc113_lw   = my_nml % l_cfc113_lw
  l_cfc113_lw2  = my_nml % l_cfc113_lw2
  l_cfc114_lw   = my_nml % l_cfc114_lw
  l_cfc114_lw2  = my_nml % l_cfc114_lw2
  l_hcfc22_lw   = my_nml % l_hcfc22_lw
  l_hcfc22_lw2  = my_nml % l_hcfc22_lw2
  l_hfc125_lw   = my_nml % l_hfc125_lw
  l_hfc125_lw2  = my_nml % l_hfc125_lw2
  l_hfc134a_lw  = my_nml % l_hfc134a_lw
  l_hfc134a_lw2 = my_nml % l_hfc134a_lw2
  l_co_lw       = my_nml % l_co_lw
  l_co_lw2      = my_nml % l_co_lw2
  l_nh3_lw      = my_nml % l_nh3_lw
  l_nh3_lw2     = my_nml % l_nh3_lw2
  l_tio_lw      = my_nml % l_tio_lw
  l_tio_lw2     = my_nml % l_tio_lw2
  l_vo_lw       = my_nml % l_vo_lw
  l_vo_lw2      = my_nml % l_vo_lw2
  l_h2_lw       = my_nml % l_h2_lw
  l_h2_lw2      = my_nml % l_h2_lw2
  l_he_lw       = my_nml % l_he_lw
  l_he_lw2      = my_nml % l_he_lw2
  l_na_lw       = my_nml % l_na_lw
  l_na_lw2      = my_nml % l_na_lw2
  l_k_lw        = my_nml % l_k_lw
  l_k_lw2       = my_nml % l_k_lw2
  l_li_lw       = my_nml % l_li_lw
  l_li_lw2      = my_nml % l_li_lw2
  l_rb_lw       = my_nml % l_rb_lw
  l_rb_lw2      = my_nml % l_rb_lw2
  l_cs_lw       = my_nml % l_cs_lw
  l_cs_lw2      = my_nml % l_cs_lw2
  l_microphysics_lw  = my_nml % l_microphysics_lw
  l_microphysics_lw2 = my_nml % l_microphysics_lw2
  ! end of logicals
  spectral_file_lw  = my_nml % spectral_file_lw
  spectral_file_lw2 = my_nml % spectral_file_lw2

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_r2lwclnl

END MODULE lw_rad_input_mod
