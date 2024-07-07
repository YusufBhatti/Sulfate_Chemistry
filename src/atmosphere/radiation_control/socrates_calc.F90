! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Call the core radiation calculations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE socrates_calc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOCRATES_CALC_MOD'
CONTAINS

SUBROUTINE socrates_calc(                                                      &
  spectrum, control, diag, dimen, atm, cld, aer, bound,                        &
  radout, radout_clean, radout_forc)

USE rad_pcf,      ONLY: i_normal
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm
USE def_cld,      ONLY: StrCld
USE def_aer,      ONLY: StrAer
USE def_bound,    ONLY: StrBound
USE def_out,      ONLY: StrOut
USE def_diag,     ONLY: StrDiag

USE gas_list_pcf, ONLY: ip_co2, ip_n2o, ip_ch4, ip_o2, ip_cfc11, ip_cfc12,     &
                        ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a
USE coradoca,     ONLY: c2c_aerosol,                                           &
                        c2c_co2, co2_mmr_scl, co2_mmr_add,                     &
                        c2c_n2o, n2o_mmr_scl, n2o_mmr_add,                     &
                        c2c_ch4, ch4_mmr_scl, ch4_mmr_add,                     &
                        c2c_o2, o2_mmr_scl, o2_mmr_add,                        &
                        c2c_cfc11, cfc11_mmr_scl, cfc11_mmr_add,               &
                        c2c_cfc12, cfc12_mmr_scl, cfc12_mmr_add,               &
                        c2c_c113, cfc113_mmr_scl, cfc113_mmr_add,              &
                        c2c_hcfc22, hcfc22_mmr_scl, hcfc22_mmr_add,            &
                        c2c_hfc125, hfc125_mmr_scl, hfc125_mmr_add,            &
                        c2c_hfc134, hfc134a_mmr_scl, hfc134a_mmr_add

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim

IMPLICIT NONE


! Spectral data:
TYPE (StrSpecData), INTENT(IN) :: spectrum

! Controlling options:
TYPE (StrCtrl), INTENT(IN) :: control

! Dimensions:
TYPE (StrDim), INTENT(IN) :: dimen

! Atmospheric properties:
TYPE(StrAtm), INTENT(IN) :: atm

! Cloud properties:
TYPE(StrCld), INTENT(IN) :: cld

! Aerosol properties:
TYPE(StrAer), INTENT(IN) :: aer

! Boundary conditions:
TYPE(StrBound), INTENT(IN) :: bound

! Diagnostics:
TYPE (StrDiag), INTENT(IN) :: diag

! Output fields from core radiation code:
TYPE(StrOut), INTENT(INOUT) :: radout

! Output fields from diagnostic calls:
TYPE(StrOut), INTENT(INOUT) :: radout_clean
TYPE(StrOut), INTENT(INOUT) :: radout_forc


! Loop variable
INTEGER ::  k

! Controlling options for diagnostic calls:
TYPE (StrCtrl) :: control_clean
TYPE (StrCtrl) :: control_forc

! Atmospheric properties for diagnostic calls:
TYPE(StrAtm) :: atm_forc

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SOCRATES_CALC'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Separate call to radiance_calc if clean-air diagnostics are required
IF (diag%l_flux_up_clean .OR. diag%l_flux_down_clean .OR.                      &
    diag%l_flux_up_clear_clean .OR. diag%l_flux_down_clear_clean .OR.          &
    diag%l_flux_up_clean_band .OR. diag%l_flux_down_clean_band .OR.            &
    diag%l_flux_up_clear_clean_band .OR.                                       &
    diag%l_flux_down_clear_clean_band .OR.                                     &
    diag%l_flux_direct_clean_sph .OR.                                          &
    diag%l_flux_direct_clean_div .OR.                                          &
    diag%l_flux_direct_clear_clean_sph .OR.                                    &
    diag%l_flux_direct_clear_clean_div .OR.                                    &
    diag%l_flux_direct_clean_sph_band .OR.                                     &
    diag%l_flux_direct_clean_div_band .OR.                                     &
    diag%l_flux_direct_clear_clean_sph_band .OR.                               &
    diag%l_flux_direct_clear_clean_div_band .OR.                               &
    diag%l_emission_spectrum_clean .OR.                                        &
    diag%l_emission_spectrum_clear_clean .OR.                                  &
    diag%l_transmission_spectrum_clean .OR.                                    &
    diag%l_transmission_spectrum_clear_clean) THEN
  
  control_clean = control
  control_clean%l_aerosol = .FALSE.
  control_clean%l_aerosol_mode = .FALSE.
  control_clean%l_clear = diag%l_flux_up_clear_clean .OR.                      &
                          diag%l_flux_down_clear_clean .OR.                    &
                          diag%l_flux_up_clear_clean_band .OR.                 &
                          diag%l_flux_down_clear_clean_band .OR.               &
                          diag%l_flux_direct_clear_clean_sph .OR.              &
                          diag%l_flux_direct_clear_clean_div .OR.              &
                          diag%l_flux_direct_clear_clean_sph_band .OR.         &
                          diag%l_flux_direct_clear_clean_div_band .OR.         &
                          diag%l_emission_spectrum_clear_clean .OR.            &
                          diag%l_transmission_spectrum_clear_clean

  control_clean%l_flux_up_band =                                               &
                          diag%l_flux_up_clean_band .OR.                       &
                          diag%l_emission_spectrum_clean .OR.                  &
                          diag%l_transmission_spectrum_clean
  control_clean%l_flux_down_band =                                             &
                          diag%l_flux_down_clean_band
  control_clean%l_flux_up_clear_band =                                         &
                          diag%l_flux_up_clear_clean_band .OR.                 &
                          diag%l_emission_spectrum_clear_clean .OR.            &
                          diag%l_transmission_spectrum_clear_clean
  control_clean%l_flux_down_clear_band =                                       &
                          diag%l_flux_down_clear_clean_band
  control_clean%l_flux_direct_sph_band =                                       &
                          diag%l_flux_direct_clean_sph_band .OR.               &
                          diag%l_transmission_spectrum_clean
  control_clean%l_flux_direct_clear_sph_band =                                 &
                          diag%l_flux_direct_clear_clean_sph_band .OR.         &
                          diag%l_transmission_spectrum_clear_clean
  control_clean%l_flux_direct_div_band =                                       &
                          diag%l_flux_direct_clean_div_band
  control_clean%l_flux_direct_clear_div_band =                                 &
                          diag%l_flux_direct_clear_clean_div_band

  control_clean%l_aerosol_absorption_band = .FALSE.
  control_clean%l_aerosol_scattering_band = .FALSE.
  control_clean%l_aerosol_asymmetry_band  = .FALSE.
  control_clean%l_cloud_extinction        = .FALSE.
  control_clean%l_ls_cloud_extinction     = .FALSE.
  control_clean%l_cnv_cloud_extinction    = .FALSE.
  control_clean%l_cloud_absorptivity      = .FALSE.
  control_clean%l_ls_cloud_absorptivity   = .FALSE.
  control_clean%l_cnv_cloud_absorptivity  = .FALSE.

! DEPENDS ON: radiance_calc
  CALL radiance_calc(control_clean, dimen, spectrum, atm, cld, aer, bound,     &
                     radout_clean)
END IF


! Separate call to radiance_calc if GHG forcing diagnostics are required
IF (diag%l_flux_up_forc .OR. diag%l_flux_down_forc .OR.                        &
    diag%l_flux_up_clear_forc .OR. diag%l_flux_down_clear_forc .OR.            &
    diag%l_flux_up_forc_band .OR. diag%l_flux_down_forc_band .OR.              &
    diag%l_flux_up_clear_forc_band .OR.                                        &
    diag%l_flux_down_clear_forc_band .OR.                                      &
    diag%l_flux_direct_sph_forc .OR.                                           &
    diag%l_flux_direct_div_forc .OR.                                           &
    diag%l_flux_direct_clear_sph_forc .OR.                                     &
    diag%l_flux_direct_clear_div_forc .OR.                                     &
    diag%l_flux_direct_sph_forc_band .OR.                                      &
    diag%l_flux_direct_div_forc_band .OR.                                      &
    diag%l_flux_direct_clear_sph_forc_band .OR.                                &
    diag%l_flux_direct_clear_div_forc_band) THEN
  control_forc = control
  atm_forc = atm
  IF (c2c_aerosol) THEN
    control_forc%l_aerosol = .FALSE.
    control_forc%l_aerosol_mode = .FALSE.
  END IF
  IF (c2c_co2 .AND. control_forc%l_co2) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_co2) THEN
        atm_forc%gas_mix_ratio(:, :, k) = co2_mmr_add +                        &
          atm_forc%gas_mix_ratio(:, :, k)*co2_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_n2o .AND. control_forc%l_n2o) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_n2o) THEN
        atm_forc%gas_mix_ratio(:, :, k) = n2o_mmr_add +                        &
          atm_forc%gas_mix_ratio(:, :, k)*n2o_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_ch4 .AND. control_forc%l_ch4) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_ch4) THEN
        atm_forc%gas_mix_ratio(:, :, k) = ch4_mmr_add +                        &
          atm_forc%gas_mix_ratio(:, :, k)*ch4_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_o2 .AND. control_forc%l_o2) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_o2) THEN
        atm_forc%gas_mix_ratio(:, :, k) = o2_mmr_add +                         &
          atm_forc%gas_mix_ratio(:, :, k)*o2_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_cfc11 .AND. control_forc%l_cfc11) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_cfc11) THEN
        atm_forc%gas_mix_ratio(:, :, k) = cfc11_mmr_add +                      &
          atm_forc%gas_mix_ratio(:, :, k)*cfc11_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_cfc12 .AND. control_forc%l_cfc12) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_cfc12) THEN
        atm_forc%gas_mix_ratio(:, :, k) = cfc12_mmr_add +                      &
          atm_forc%gas_mix_ratio(:, :, k)*cfc12_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_c113 .AND. control_forc%l_cfc113) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_cfc113) THEN
        atm_forc%gas_mix_ratio(:, :, k) = cfc113_mmr_add +                     &
          atm_forc%gas_mix_ratio(:, :, k)*cfc113_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_hcfc22 .AND. control_forc%l_hcfc22) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_hcfc22) THEN
        atm_forc%gas_mix_ratio(:, :, k) = hcfc22_mmr_add +                     &
          atm_forc%gas_mix_ratio(:, :, k)*hcfc22_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_hfc125 .AND. control_forc%l_hfc125) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_hfc125) THEN
        atm_forc%gas_mix_ratio(:, :, k) = hfc125_mmr_add +                     &
          atm_forc%gas_mix_ratio(:, :, k)*hfc125_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  IF (c2c_hfc134 .AND. control_forc%l_hfc134a) THEN
    DO k=1, spectrum%gas%n_absorb
      IF (spectrum%gas%type_absorb(k) == ip_hfc134a) THEN
        atm_forc%gas_mix_ratio(:, :, k) = hfc134a_mmr_add +                    &
          atm_forc%gas_mix_ratio(:, :, k)*hfc134a_mmr_scl
        EXIT
      END IF
    END DO
  END IF
  control_forc%l_clear = diag%l_flux_up_clear_forc .OR.                        &
                         diag%l_flux_down_clear_forc .OR.                      &
                         diag%l_flux_up_clear_forc_band .OR.                   &
                         diag%l_flux_down_clear_forc_band .OR.                 &
                         diag%l_flux_direct_clear_sph_forc .OR.                &
                         diag%l_flux_direct_clear_div_forc .OR.                &
                         diag%l_flux_direct_clear_sph_forc_band .OR.           &
                         diag%l_flux_direct_clear_div_forc_band

  control_forc%l_flux_up_band = diag%l_flux_up_forc_band
  control_forc%l_flux_down_band = diag%l_flux_down_forc_band
  control_forc%l_flux_up_clear_band = diag%l_flux_up_clear_forc_band
  control_forc%l_flux_down_clear_band = diag%l_flux_down_clear_forc_band
  control_forc%l_flux_direct_sph_band = diag%l_flux_direct_sph_forc_band
  control_forc%l_flux_direct_clear_sph_band =                                 &
                                        diag%l_flux_direct_clear_sph_forc_band
  control_forc%l_flux_direct_div_band = diag%l_flux_direct_div_forc_band
  control_forc%l_flux_direct_clear_div_band =                                 &
                                        diag%l_flux_direct_clear_div_forc_band

  control_forc%l_aerosol_absorption_band = .FALSE.
  control_forc%l_aerosol_scattering_band = .FALSE.
  control_forc%l_aerosol_asymmetry_band  = .FALSE.
  control_forc%l_cloud_extinction        = .FALSE.
  control_forc%l_ls_cloud_extinction     = .FALSE.
  control_forc%l_cnv_cloud_extinction    = .FALSE.
  control_forc%l_cloud_absorptivity      = .FALSE.
  control_forc%l_ls_cloud_absorptivity   = .FALSE.
  control_forc%l_cnv_cloud_absorptivity  = .FALSE.

  CALL radiance_calc(control_forc, dimen, spectrum, atm_forc, cld, aer, bound, &
                     radout_forc)
END IF

! Attention: radout is allocated inside routine.
!------------------------------------------------------------------------------
CALL radiance_calc(control, dimen, spectrum, atm, cld, aer, bound, radout)
!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE socrates_calc
END MODULE socrates_calc_mod
