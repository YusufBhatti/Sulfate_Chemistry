! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE compress_spectrum_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'COMPRESS_SPECTRUM_MOD'
CONTAINS

SUBROUTINE compress_spectrum(con, spec)

! Subroutine to compress the spectral file data
!
! Purpose:
!   Spectral data from the full spectrum is reduced to only
!   those properites required.

USE def_control,         ONLY: StrCtrl
USE def_spectrum,        ONLY: StrSpecData
USE gas_list_pcf
USE run_aerosol_mod,     ONLY: l_soot, l_biomass, l_ocff, l_sulpc_so2
USE dust_parameters_mod, ONLY: l_dust, l_twobin_dust
USE rad_input_mod,       ONLY: l_use_arcldust, l_use_arclsulp, l_use_arclocff, &
                               l_use_arclblck, l_use_arclbiom, l_use_arclsslt, &
                               l_use_seasalt_direct
USE missing_data_mod,    ONLY: rmdi
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE r2_set_690nm_weight_mod, ONLY: r2_set_690nm_weight
IMPLICIT NONE


TYPE (StrCtrl),     INTENT(IN)    :: con
TYPE (StrSpecData), INTENT(INOUT) :: spec

! Local variables
INTEGER :: i, j, n_band_absorb, n_aerosol_mr
LOGICAL :: l_retain_absorb(spec%gas%n_absorb)
!   Flags for the retention of gases in the spectral file

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPRESS_SPECTRUM'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Search the spectrum to find those gases to be retained.
l_retain_absorb=.FALSE.
DO i=1, spec%gas%n_absorb
  IF (((spec%gas%type_absorb(i) == ip_h2o)     .AND. con%l_h2o    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_co2)     .AND. con%l_co2    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_o3)      .AND. con%l_o3     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_o2)      .AND. con%l_o2     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_n2o)     .AND. con%l_n2o    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_ch4)     .AND. con%l_ch4    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_so2)     .AND. con%l_so2    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc11)   .AND. con%l_cfc11  ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc12)   .AND. con%l_cfc12  ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc113)  .AND. con%l_cfc113 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cfc114)  .AND. con%l_cfc114 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hcfc22)  .AND. con%l_hcfc22 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hfc125)  .AND. con%l_hfc125 ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_hfc134a) .AND. con%l_hfc134a) .OR.       &
      ((spec%gas%type_absorb(i) == ip_co)      .AND. con%l_co     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_nh3)     .AND. con%l_nh3    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_tio)     .AND. con%l_tio    ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_vo)      .AND. con%l_vo     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_h2)      .AND. con%l_h2     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_he)      .AND. con%l_he     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_na)      .AND. con%l_na     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_k)       .AND. con%l_k      ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_li)      .AND. con%l_li     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_rb)      .AND. con%l_rb     ) .OR.       &
      ((spec%gas%type_absorb(i) == ip_cs)      .AND. con%l_cs     )) THEN
    l_retain_absorb(i)=.TRUE.
  END IF
END DO

DO i=1, spec%basic%n_band
  n_band_absorb=0
  DO j=1, spec%gas%n_band_absorb(i)
    IF (l_retain_absorb(spec%gas%index_absorb(j, i))) THEN
      n_band_absorb = n_band_absorb + 1
      spec%gas%index_absorb(n_band_absorb, i) = spec%gas%index_absorb(j, i)
    END IF
  END DO
  spec%gas%n_band_absorb(i)=n_band_absorb
END DO


! Increase number of aerosol mixing ratios as required: if both
! the prognostic and climatological aerosols are present then
! both need to be stored in the aerosol_mix_ratio array, depending on the
! number of species used for that aerosol type.
n_aerosol_mr = spec%aerosol%n_aerosol_mr
IF (l_biomass .AND. l_use_arclbiom) n_aerosol_mr = n_aerosol_mr + 2
IF (l_soot .AND. l_use_arclblck) n_aerosol_mr = n_aerosol_mr + 2
IF (l_use_seasalt_direct .AND. l_use_arclsslt) n_aerosol_mr = n_aerosol_mr + 2
IF (l_sulpc_so2 .AND. l_use_arclsulp) n_aerosol_mr = n_aerosol_mr + 2
IF (l_dust .AND. l_use_arcldust .AND. (.NOT. l_twobin_dust)) THEN
    n_aerosol_mr = n_aerosol_mr + 6
END IF
IF (l_ocff .AND. l_use_arclocff) n_aerosol_mr = n_aerosol_mr + 2
spec%aerosol%n_aerosol_mr = n_aerosol_mr
spec%dim%nd_aerosol_mr    = n_aerosol_mr


! Set "blue" weight for each band if not read from spectral file
IF (spec%solar%weight_blue(1) == rmdi) THEN
  CALL r2_set_690nm_weight(spec%basic%n_band,                                  &
    spec%basic%l_present,                                                      &
    spec%basic%n_band_exclude,                                                 &
    spec%basic%index_exclude,                                                  &
    spec%basic%wavelength_short,                                               &
    spec%basic%wavelength_long,                                                &
    spec%solar%weight_blue,                                                    &
    spec%dim%nd_band, spec%dim%nd_exclude,                                     &
    spec%dim%nd_type)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE compress_spectrum
END MODULE compress_spectrum_mod
