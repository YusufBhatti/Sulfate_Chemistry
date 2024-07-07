! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Re-arrange UKCA-MODE input to match the expectations of
!  routine ukca_radaer_band_average().
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_prepare_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_PREPARE_MOD'

CONTAINS

SUBROUTINE ukca_radaer_prepare(                                         &
      ! Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                     &
      ! UKCA_RADAER structure
      ukca_radaer,                                                      &
      ! Component mass-mixing ratios
      ukca_mix_ratio,                                                   &
      ! Modal mass-mixing ratios
      ukca_modal_mixr,                                                  &
      ! Input modal number concentrations
      ukca_modal_nbr,                                                   &
      ! Output modal number concentrations
      ukca_modal_number,                                                &
      ! Pressure and temperature
      pressure, temperature,                                            &
      ! Fixed array dimensions
      npd_profile, npd_layer, npd_aerosol_mode                          &
  )

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE ukca_radaer_struct_mod

IMPLICIT NONE

!
! Arguments with intent(in)
!
!
! Fixed array dimensions
!
INTEGER :: npd_profile, &
           npd_layer,   &
           npd_aerosol_mode
!
! Actual array dimensions
!
INTEGER :: n_profile,   &
           n_layer,     &
           n_ukca_mode, &
           n_ukca_cpnt
!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer

!
! Component mass-mixing ratios
!
REAL :: ukca_mix_ratio(npd_profile, npd_layer, n_ukca_cpnt)

!
! Modal number concentrations divided by molecular concentration of air
!
REAL :: ukca_modal_nbr(npd_profile, npd_layer, n_ukca_mode)

!
! Pressure and temperature fields.
REAL :: pressure(npd_profile, npd_layer), &
        temperature(npd_profile, npd_layer)
!
!
! Arguments with intent(out)
!
!
!
! Modal mass-mixing ratios
!
REAL :: ukca_modal_mixr(npd_profile, npd_layer, npd_aerosol_mode)

!
! Modal number concentrations (in m-3)
!
REAL :: ukca_modal_number(npd_profile, npd_layer, n_ukca_mode)

!
! Local variables
!
INTEGER :: i, &
        j, &
        k, &
        l
INTEGER :: this_cpnt

!
! Boltzmann constant
!
REAL, PARAMETER :: k_boltzmann = 1.3807e-23 ! J/K

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_PREPARE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!
! Modal mass-mixing ratios.
!
! Simply sum the mixing ratios of all components included in a
! given mode to get the mixing ratio for that mode.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer

    DO l = 1, n_profile

      ukca_modal_mixr(l, k, j) = 0.0

    END DO ! l

  END DO ! k

  DO i = 1, ukca_radaer%n_cpnt_in_mode(j)

    this_cpnt = ukca_radaer%i_cpnt_index(i, j)

    DO k = 1, n_layer

      DO l = 1, n_profile

        ukca_modal_mixr(l, k, j) = &
          ukca_modal_mixr(l, k, j) + ukca_mix_ratio(l, k, this_cpnt)

      END DO ! l

    END DO ! k

  END DO ! i

END DO ! j

!
! Modal number concentrations
!
! Multiply by the molecular concentration of air (p/kT) to obtain
! the acutal aerosol number concentrations.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer

    DO l = 1, n_profile

      ukca_modal_number(l, k, j) = ukca_modal_nbr(l, k, j) * &
                                   pressure(l, k) / &
                                   (k_boltzmann * temperature(l, k))

    END DO ! l

  END DO ! k

END DO ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_prepare
END MODULE ukca_radaer_prepare_mod
