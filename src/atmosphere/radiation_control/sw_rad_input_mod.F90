! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for sw radiation.

! Description:
!   Module containing control as used by the sw radiation code.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE sw_rad_input_mod

USE sat_opt_mod
USE sw_control_struct, ONLY: sw_control
USE rad_input_mod, ONLY: n_swcall
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER (LEN=132) :: spectral_file_sw     ! Spectral file

LOGICAL :: l_h2o_sw                  ! flag for water vapour
LOGICAL :: l_co2_sw                  ! flag for carbon dioxide
LOGICAL :: l_o3_sw                   ! flag for ozone
LOGICAL :: l_o2_sw                   ! flag for oxygen
LOGICAL :: l_so2_sw                  ! flag for sulphur dioxide
LOGICAL :: l_ch4_sw                  ! flag for methane
LOGICAL :: l_n2o_sw                  ! flag for N2O
LOGICAL :: l_co_sw                   ! flag for carbon monoxide
LOGICAL :: l_nh3_sw                  ! flag for ammonia
LOGICAL :: l_tio_sw                  ! flag for titanium oxide
LOGICAL :: l_vo_sw                   ! flag for vanadium oxide
LOGICAL :: l_h2_sw                   ! flag for H2-H2 CIA
LOGICAL :: l_he_sw                   ! flag for H2-He CIA
LOGICAL :: l_na_sw                   ! flag for sodium
LOGICAL :: l_k_sw                    ! flag for potassium
LOGICAL :: l_li_sw                   ! flag for lithium
LOGICAL :: l_rb_sw                   ! flag for rubidium
LOGICAL :: l_cs_sw                   ! flag for cesium

LOGICAL :: l_solvar_sw               ! treatment of solar variability
LOGICAL :: l_spherical_solar         ! spherical geometry for direct beam

INTEGER :: i_gas_overlap_sw          ! treatment of gaseous overlaps

! types of droplets or ice crystals used for parametrizations
INTEGER :: i_st_water_sw             ! type for stratiform water
INTEGER :: i_cnv_water_sw            ! type for convective water
INTEGER :: i_st_ice_sw               ! type for stratiform ice
INTEGER :: i_cnv_ice_sw              ! type for convective ice

! Second set of control options for the second radiation call:
CHARACTER (LEN=132) :: spectral_file_sw2
LOGICAL :: l_h2o_sw2
LOGICAL :: l_co2_sw2
LOGICAL :: l_o3_sw2
LOGICAL :: l_o2_sw2
LOGICAL :: l_so2_sw2
LOGICAL :: l_ch4_sw2
LOGICAL :: l_n2o_sw2
LOGICAL :: l_co_sw2
LOGICAL :: l_nh3_sw2
LOGICAL :: l_tio_sw2
LOGICAL :: l_vo_sw2
LOGICAL :: l_h2_sw2
LOGICAL :: l_he_sw2
LOGICAL :: l_na_sw2
LOGICAL :: l_k_sw2
LOGICAL :: l_li_sw2
LOGICAL :: l_rb_sw2
LOGICAL :: l_cs_sw2
LOGICAL :: l_solvar_sw2
LOGICAL :: l_spherical_solar_2
INTEGER :: i_gas_overlap_sw2
INTEGER :: i_st_water_sw2
INTEGER :: i_cnv_water_sw2
INTEGER :: i_st_ice_sw2
INTEGER :: i_cnv_ice_sw2

! -----------------------
!     Elements of namelists for options in the radiation code.
NAMELIST/r2swclnl/                                                      &
     spectral_file_sw, spectral_file_sw2,                               &
     i_gas_overlap_sw, i_gas_overlap_sw2,                               &
     l_h2o_sw, l_h2o_sw2, l_co2_sw, l_co2_sw2, l_o3_sw, l_o3_sw2,       &
     l_o2_sw, l_o2_sw2, l_ch4_sw, l_ch4_sw2, l_n2o_sw, l_n2o_sw2,       &
     l_so2_sw, l_so2_sw2,                                               &
     l_co_sw, l_co_sw2, l_nh3_sw, l_nh3_sw2, l_tio_sw, l_tio_sw2,       &
     l_vo_sw, l_vo_sw2, l_h2_sw, l_h2_sw2, l_he_sw, l_he_sw2,           &
     l_na_sw, l_na_sw2, l_k_sw, l_k_sw2,                                &
     l_li_sw, l_li_sw2, l_rb_sw, l_rb_sw2, l_cs_sw, l_cs_sw2,           &
     i_st_water_sw, i_st_water_sw2, i_cnv_water_sw, i_cnv_water_sw2,    &
     i_st_ice_sw, i_st_ice_sw2, i_cnv_ice_sw, i_cnv_ice_sw2,            &
     l_solvar_sw, l_solvar_sw2, l_spherical_solar, l_spherical_solar_2

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ----------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SW_RAD_INPUT_MOD'

CONTAINS

! Subroutine to set the input values of the sw control structure

SUBROUTINE sw_input

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica, ip_cloud_clear,                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan,                                               &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand
USE ereport_mod, ONLY: ereport
USE rad_input_mod
USE file_manager, ONLY: get_file_unit_by_id

USE model_domain_mod, ONLY: model_type, mt_single_column
USE sw_control_default_mod, ONLY: sw_control_default

IMPLICIT NONE

INTEGER :: j,                                                     &
           errorstatus      ! Return code : 0 Normal Exit : >0 Error
INTEGER :: atmoscntl_unit

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='SW_INPUT')
CHARACTER(LEN=errormessagelength) :: CMessage ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialisation: set namelist items to MDI / FALSE.
spectral_file_sw           = ''
i_gas_overlap_sw           = imdi
l_solvar_sw                = .FALSE.
l_spherical_solar          = .FALSE.
l_h2o_sw                   = .FALSE.
l_co2_sw                   = .FALSE.
l_o3_sw                    = .FALSE.
l_o2_sw                    = .FALSE.
l_so2_sw                   = .FALSE.
l_n2o_sw                   = .FALSE.
l_ch4_sw                   = .FALSE.
l_co_sw                    = .FALSE.
l_nh3_sw                   = .FALSE.
l_tio_sw                   = .FALSE.
l_vo_sw                    = .FALSE.
l_h2_sw                    = .FALSE.
l_he_sw                    = .FALSE.
l_na_sw                    = .FALSE.
l_k_sw                     = .FALSE.
l_li_sw                    = .FALSE.
l_rb_sw                    = .FALSE.
l_cs_sw                    = .FALSE.
i_st_water_sw              = imdi
i_cnv_water_sw             = imdi
i_st_ice_sw                = imdi
i_cnv_ice_sw               = imdi

spectral_file_sw2          = ''
i_gas_overlap_sw2          = imdi
l_solvar_sw2               = .FALSE.
l_spherical_solar_2        = .FALSE.
l_h2o_sw2                  = .FALSE.
l_co2_sw2                  = .FALSE.
l_o3_sw2                   = .FALSE.
l_o2_sw2                   = .FALSE.
l_so2_sw2                  = .FALSE.
l_n2o_sw2                  = .FALSE.
l_ch4_sw2                  = .FALSE.
l_co_sw2                   = .FALSE.
l_nh3_sw2                  = .FALSE.
l_tio_sw2                  = .FALSE.
l_vo_sw2                   = .FALSE.
l_h2_sw2                   = .FALSE.
l_he_sw2                   = .FALSE.
l_na_sw2                   = .FALSE.
l_k_sw2                    = .FALSE.
l_li_sw2                   = .FALSE.
l_rb_sw2                   = .FALSE.
l_cs_sw2                   = .FALSE.
i_st_water_sw2             = imdi
i_cnv_water_sw2            = imdi
i_st_ice_sw2               = imdi
i_cnv_ice_sw2              = imdi

! Options for the shortwave
SELECT CASE (model_type)

CASE (mt_single_column)
  READ (10, r2swclnl)

CASE DEFAULT
  atmoscntl_unit = get_file_unit_by_id("atmoscntl", handler="fortran")
  CALL read_nml_r2swclnl(atmoscntl_unit)

END SELECT

DO j=1,n_swcall
  ! Set default values of control variables.
  ! Note: to allow a simple interface for ROSE only the commonly used
  ! items are now set via the namelist. In order to set further variables
  ! (eg. for radiance diagnostics) a branch must now be used.

  CALL sw_control_default(sw_control(j))

  ! Transfer namelist items to the data structure.
  IF (j==2) THEN
    sw_control(j)%spectral_file          = ADJUSTL(spectral_file_sw2)
    sw_control(j)%i_gas_overlap          = i_gas_overlap_sw2
    sw_control(j)%l_solvar               = l_solvar_sw2
    sw_control(j)%l_spherical_solar      = l_spherical_solar_2
    sw_control(j)%l_continuum            = l_h2o_sw2
    sw_control(j)%l_h2o                  = l_h2o_sw2
    sw_control(j)%l_co2                  = l_co2_sw2
    sw_control(j)%l_o3                   = l_o3_sw2
    sw_control(j)%l_o2                   = l_o2_sw2
    sw_control(j)%l_so2                  = l_so2_sw2
    sw_control(j)%l_n2o                  = l_n2o_sw2
    sw_control(j)%l_ch4                  = l_ch4_sw2
    sw_control(j)%l_co                   = l_co_sw2
    sw_control(j)%l_nh3                  = l_nh3_sw2
    sw_control(j)%l_tio                  = l_tio_sw2
    sw_control(j)%l_vo                   = l_vo_sw2
    sw_control(j)%l_h2                   = l_h2_sw2
    sw_control(j)%l_he                   = l_he_sw2
    sw_control(j)%l_na                   = l_na_sw2
    sw_control(j)%l_k                    = l_k_sw2
    sw_control(j)%l_li                   = l_li_sw2
    sw_control(j)%l_rb                   = l_rb_sw2
    sw_control(j)%l_cs                   = l_cs_sw2
    sw_control(j)%i_st_water             = i_st_water_sw2
    sw_control(j)%i_cnv_water            = i_cnv_water_sw2
    sw_control(j)%i_st_ice               = i_st_ice_sw2
    sw_control(j)%i_cnv_ice              = i_cnv_ice_sw2
    sw_control(j)%i_cloud_representation = i_cloud_representation_2
    sw_control(j)%i_overlap              = i_overlap_2
    sw_control(j)%i_inhom                = i_inhom_2
  ELSE
    sw_control(j)%spectral_file          = ADJUSTL(spectral_file_sw)
    sw_control(j)%i_gas_overlap          = i_gas_overlap_sw
    sw_control(j)%l_solvar               = l_solvar_sw
    sw_control(j)%l_spherical_solar      = l_spherical_solar
    sw_control(j)%l_continuum            = l_h2o_sw
    sw_control(j)%l_h2o                  = l_h2o_sw
    sw_control(j)%l_co2                  = l_co2_sw
    sw_control(j)%l_o3                   = l_o3_sw
    sw_control(j)%l_o2                   = l_o2_sw
    sw_control(j)%l_so2                  = l_so2_sw
    sw_control(j)%l_n2o                  = l_n2o_sw
    sw_control(j)%l_ch4                  = l_ch4_sw
    sw_control(j)%l_co                   = l_co_sw
    sw_control(j)%l_nh3                  = l_nh3_sw
    sw_control(j)%l_tio                  = l_tio_sw
    sw_control(j)%l_vo                   = l_vo_sw
    sw_control(j)%l_h2                   = l_h2_sw
    sw_control(j)%l_he                   = l_he_sw
    sw_control(j)%l_na                   = l_na_sw
    sw_control(j)%l_k                    = l_k_sw
    sw_control(j)%l_li                   = l_li_sw
    sw_control(j)%l_rb                   = l_rb_sw
    sw_control(j)%l_cs                   = l_cs_sw
    sw_control(j)%i_st_water             = i_st_water_sw
    sw_control(j)%i_cnv_water            = i_cnv_water_sw
    sw_control(j)%i_st_ice               = i_st_ice_sw
    sw_control(j)%i_cnv_ice              = i_cnv_ice_sw
    sw_control(j)%i_cloud_representation = i_cloud_representation
    sw_control(j)%i_overlap              = i_overlap
    sw_control(j)%i_inhom                = i_inhom
  END IF
END DO

! Set i_cloud from combination of i_overlap and i_representation
DO j=1,n_swcall
  IF (sw_control(j)%i_cloud_representation == ip_cloud_homogen) THEN

    sw_control(j)%i_solver=ip_solver_homogen_direct
    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      sw_control(j)%i_cloud=ip_cloud_mcica
    ELSE
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF (sw_control(j)%i_cloud_representation == ip_cloud_ice_water) THEN

    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      sw_control(j)%i_cloud=ip_cloud_mcica
      sw_control(j)%i_solver=ip_solver_homogen_direct
    ELSE
      sw_control(j)%i_solver=ip_solver_mix_direct_hogan
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE IF ((sw_control(j)%i_cloud_representation == ip_cloud_conv_strat) .OR.  &
           (sw_control(j)%i_cloud_representation == ip_cloud_csiw)) THEN

    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//                 &
                 ' cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage)
    ELSE
      sw_control(j)%i_solver=ip_solver_triple_hogan
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_triple
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage)
      END IF
    END IF

  ELSE

    ! No treatment of cloud for SW radiation
    sw_control(j)%l_cloud        = .FALSE.
    sw_control(j)%i_cloud        = ip_cloud_clear
    sw_control(j)%i_solver       = ip_solver_homogen_direct
    sw_control(j)%l_microphysics = .FALSE.

  END IF
  IF (sw_control(j)%i_inhom == ip_scaling) THEN
    l_inhom_cloud=.TRUE.
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sw_input

SUBROUTINE read_nml_r2swclnl(unit_in)

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
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_R2SWCLNL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 10
INTEGER, PARAMETER :: n_log = 40
INTEGER, PARAMETER :: n_chars =  2 * 132

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_gas_overlap_sw
  INTEGER :: i_gas_overlap_sw2
  INTEGER :: i_st_water_sw
  INTEGER :: i_st_water_sw2
  INTEGER :: i_cnv_water_sw
  INTEGER :: i_cnv_water_sw2
  INTEGER :: i_st_ice_sw
  INTEGER :: i_st_ice_sw2
  INTEGER :: i_cnv_ice_sw
  INTEGER :: i_cnv_ice_sw2
  LOGICAL :: l_h2o_sw
  LOGICAL :: l_h2o_sw2
  LOGICAL :: l_co2_sw
  LOGICAL :: l_co2_sw2
  LOGICAL :: l_o3_sw
  LOGICAL :: l_o3_sw2
  LOGICAL :: l_o2_sw
  LOGICAL :: l_o2_sw2
  LOGICAL :: l_so2_sw
  LOGICAL :: l_so2_sw2
  LOGICAL :: l_ch4_sw
  LOGICAL :: l_ch4_sw2
  LOGICAL :: l_n2o_sw
  LOGICAL :: l_n2o_sw2
  LOGICAL :: l_co_sw
  LOGICAL :: l_co_sw2
  LOGICAL :: l_nh3_sw
  LOGICAL :: l_nh3_sw2
  LOGICAL :: l_tio_sw
  LOGICAL :: l_tio_sw2
  LOGICAL :: l_vo_sw
  LOGICAL :: l_vo_sw2
  LOGICAL :: l_h2_sw
  LOGICAL :: l_h2_sw2
  LOGICAL :: l_he_sw
  LOGICAL :: l_he_sw2
  LOGICAL :: l_na_sw
  LOGICAL :: l_na_sw2
  LOGICAL :: l_k_sw
  LOGICAL :: l_k_sw2
  LOGICAL :: l_li_sw
  LOGICAL :: l_li_sw2
  LOGICAL :: l_rb_sw
  LOGICAL :: l_rb_sw2
  LOGICAL :: l_cs_sw
  LOGICAL :: l_cs_sw2
  LOGICAL :: l_solvar_sw
  LOGICAL :: l_solvar_sw2
  LOGICAL :: l_spherical_solar
  LOGICAL :: l_spherical_solar_2
  CHARACTER (LEN=132) :: spectral_file_sw
  CHARACTER (LEN=132) :: spectral_file_sw2
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=r2swclnl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist R2SWCLNL", iomessage)

  my_nml % i_gas_overlap_sw  = i_gas_overlap_sw
  my_nml % i_gas_overlap_sw2 = i_gas_overlap_sw2
  my_nml % i_st_water_sw     = i_st_water_sw
  my_nml % i_st_water_sw2    = i_st_water_sw2
  my_nml % i_cnv_water_sw    = i_cnv_water_sw
  my_nml % i_cnv_water_sw2   = i_cnv_water_sw2
  my_nml % i_st_ice_sw       = i_st_ice_sw
  my_nml % i_st_ice_sw2      = i_st_ice_sw2
  my_nml % i_cnv_ice_sw      = i_cnv_ice_sw
  my_nml % i_cnv_ice_sw2     = i_cnv_ice_sw2
  ! end of integers
  my_nml % l_h2o_sw  = l_h2o_sw
  my_nml % l_h2o_sw2 = l_h2o_sw2
  my_nml % l_co2_sw  = l_co2_sw
  my_nml % l_co2_sw2 = l_co2_sw2
  my_nml % l_o3_sw   = l_o3_sw
  my_nml % l_o3_sw2  = l_o3_sw2
  my_nml % l_o2_sw   = l_o2_sw
  my_nml % l_o2_sw2  = l_o2_sw2
  my_nml % l_so2_sw  = l_so2_sw
  my_nml % l_so2_sw2 = l_so2_sw2
  my_nml % l_ch4_sw  = l_ch4_sw
  my_nml % l_ch4_sw2 = l_ch4_sw2
  my_nml % l_n2o_sw  = l_n2o_sw
  my_nml % l_n2o_sw2 = l_n2o_sw2
  my_nml % l_co_sw   = l_co_sw
  my_nml % l_co_sw2  = l_co_sw2
  my_nml % l_nh3_sw  = l_nh3_sw
  my_nml % l_nh3_sw2 = l_nh3_sw2
  my_nml % l_tio_sw  = l_tio_sw
  my_nml % l_tio_sw2 = l_tio_sw2
  my_nml % l_vo_sw   = l_vo_sw
  my_nml % l_vo_sw2  = l_vo_sw2
  my_nml % l_h2_sw   = l_h2_sw
  my_nml % l_h2_sw2  = l_h2_sw2
  my_nml % l_he_sw   = l_he_sw
  my_nml % l_he_sw2  = l_he_sw2
  my_nml % l_na_sw   = l_na_sw
  my_nml % l_na_sw2  = l_na_sw2
  my_nml % l_k_sw    = l_k_sw
  my_nml % l_k_sw2   = l_k_sw2
  my_nml % l_li_sw   = l_li_sw
  my_nml % l_li_sw2  = l_li_sw2
  my_nml % l_rb_sw   = l_rb_sw
  my_nml % l_rb_sw2  = l_rb_sw2
  my_nml % l_cs_sw   = l_cs_sw
  my_nml % l_cs_sw2  = l_cs_sw2
  my_nml % l_solvar_sw       = l_solvar_sw
  my_nml % l_solvar_sw2      = l_solvar_sw2
  my_nml % l_spherical_solar = l_spherical_solar
  my_nml % l_spherical_solar_2 = l_spherical_solar_2
  ! end of logicals
  my_nml % spectral_file_sw  = spectral_file_sw
  my_nml % spectral_file_sw2 = spectral_file_sw2

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_gas_overlap_sw  = my_nml % i_gas_overlap_sw
  i_gas_overlap_sw2 = my_nml % i_gas_overlap_sw2
  i_st_water_sw     = my_nml % i_st_water_sw
  i_st_water_sw2    = my_nml % i_st_water_sw2
  i_cnv_water_sw    = my_nml % i_cnv_water_sw
  i_cnv_water_sw2   = my_nml % i_cnv_water_sw2
  i_st_ice_sw       = my_nml % i_st_ice_sw
  i_st_ice_sw2      = my_nml % i_st_ice_sw2
  i_cnv_ice_sw      = my_nml % i_cnv_ice_sw
  i_cnv_ice_sw2     = my_nml % i_cnv_ice_sw2
  ! end of integers
  l_h2o_sw  = my_nml % l_h2o_sw
  l_h2o_sw2 = my_nml % l_h2o_sw2
  l_co2_sw  = my_nml % l_co2_sw
  l_co2_sw2 = my_nml % l_co2_sw2
  l_o3_sw   = my_nml % l_o3_sw
  l_o3_sw2  = my_nml % l_o3_sw2
  l_o2_sw   = my_nml % l_o2_sw
  l_o2_sw2  = my_nml % l_o2_sw2
  l_so2_sw  = my_nml % l_so2_sw
  l_so2_sw2 = my_nml % l_so2_sw2
  l_ch4_sw  = my_nml % l_ch4_sw
  l_ch4_sw2 = my_nml % l_ch4_sw2
  l_n2o_sw  = my_nml % l_n2o_sw
  l_n2o_sw2 = my_nml % l_n2o_sw2
  l_co_sw   = my_nml % l_co_sw
  l_co_sw2  = my_nml % l_co_sw2
  l_nh3_sw  = my_nml % l_nh3_sw
  l_nh3_sw2 = my_nml % l_nh3_sw2
  l_tio_sw  = my_nml % l_tio_sw
  l_tio_sw2 = my_nml % l_tio_sw2
  l_vo_sw   = my_nml % l_vo_sw
  l_vo_sw2  = my_nml % l_vo_sw2
  l_h2_sw   = my_nml % l_h2_sw
  l_h2_sw2  = my_nml % l_h2_sw2
  l_he_sw   = my_nml % l_he_sw
  l_he_sw2  = my_nml % l_he_sw2
  l_na_sw   = my_nml % l_na_sw
  l_na_sw2  = my_nml % l_na_sw2
  l_k_sw    = my_nml % l_k_sw
  l_k_sw2   = my_nml % l_k_sw2
  l_li_sw   = my_nml % l_li_sw
  l_li_sw2  = my_nml % l_li_sw2
  l_rb_sw   = my_nml % l_rb_sw
  l_rb_sw2  = my_nml % l_rb_sw2
  l_cs_sw   = my_nml % l_cs_sw
  l_cs_sw2  = my_nml % l_cs_sw2
  l_solvar_sw       = my_nml % l_solvar_sw
  l_solvar_sw2      = my_nml % l_solvar_sw2
  l_spherical_solar = my_nml % l_spherical_solar
  l_spherical_solar_2 = my_nml % l_spherical_solar_2
  ! end of logicals
  spectral_file_sw  = my_nml % spectral_file_sw
  spectral_file_sw2 = my_nml % spectral_file_sw2

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_r2swclnl

END MODULE sw_rad_input_mod
