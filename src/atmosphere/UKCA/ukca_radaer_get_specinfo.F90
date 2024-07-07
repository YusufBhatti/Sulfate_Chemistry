! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Returns basic information on the spectral files loaded in memory.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_get_specinfo_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_GET_SPECINFO_MOD'

CONTAINS

SUBROUTINE ukca_get_specinfo(                                     &
       nbr_band_sw,                                               &
       nbr_band_lw,                                               &
       nbr_aod_wavel,                                             &
       wavelength_short,                                          &
       wavelength_long                                            &
  )

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE spec_sw_lw, ONLY: sw_spectrum, lw_spectrum
USE ukca_radaer_precalc, ONLY: npd_ukca_spectrum, npd_ukca_band
USE spcrg3a_mod, ONLY: ip_solar, ip_infra_red

IMPLICIT NONE

!     Arguments with intent(out)
INTEGER :: nbr_band_sw
INTEGER :: nbr_band_lw
INTEGER :: nbr_aod_wavel
REAL :: wavelength_short(npd_ukca_spectrum, npd_ukca_band)
REAL :: wavelength_long(npd_ukca_spectrum, npd_ukca_band)

!     Local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_GET_SPECINFO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise
wavelength_short(:,:) = 0.0
wavelength_long(:,:)  = 0.0

!     Copy the required info

nbr_band_sw   = sw_spectrum(1)%basic%n_band
nbr_band_lw   = lw_spectrum(1)%basic%n_band
nbr_aod_wavel = lw_spectrum(1)%aerosol%n_aod_wavel
wavelength_short(ip_solar,1:nbr_band_sw) =                      &
                        sw_spectrum(1)%basic%wavelength_short(:)
wavelength_long(ip_solar,1:nbr_band_sw)  =                      &
                        sw_spectrum(1)%basic%wavelength_long(:)
wavelength_short(ip_infra_red,1:nbr_band_lw)=                   &
                        lw_spectrum(1)%basic%wavelength_short(:)
wavelength_long(ip_infra_red,1:nbr_band_lw) =                   &
                        lw_spectrum(1)%basic%wavelength_long(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_get_specinfo
END MODULE ukca_get_specinfo_mod
