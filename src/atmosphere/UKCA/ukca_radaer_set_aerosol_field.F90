! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Copy UKCA-MODE aerosol fields to arrays that follow the radiation
!  convention (top-to-bottom in the vertical, gathered horizontal
!  points)
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_set_aerosol_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'UKCA_RADAER_SET_AEROSOL_FIELD_MOD'

CONTAINS

SUBROUTINE ukca_radaer_set_aerosol_field(                               &
   i_gather, l_extra_top,                                               &
   n_layer, n_profile,                                                  &
   ukca_dim1, ukca_dim2,                                                &
   ukca_mmr, ukca_cvl,                                                  &
   ukca_dry, ukca_wet, ukca_rho, ukca_vol, ukca_wtv, ukca_nbr, trindx,  &
   ukcaaer_mix_ratio, ukcaaer_comp_vol,                                 &
   ukcaaer_dry_diam, ukcaaer_wet_diam, ukcaaer_modal_rho,               &
   ukcaaer_modal_vol, ukcaaer_modal_wtv, ukcaaer_modal_nbr, trindxrad,  &
   npd_field, npd_profile, npd_layer,                                   &
   n_ukca_cpnt, n_ukca_mode                                             &
   )

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!
! Arguments
!
!
! Fixed array dimensions
!
  ! Field size
INTEGER :: npd_field

! Size of array of profiles
INTEGER :: npd_profile

! Maximum number of layers
INTEGER :: npd_layer

! Total number of UKCA aerosol components
INTEGER :: n_ukca_cpnt

! Total number of UKCA aerosol modes
INTEGER :: n_ukca_mode
!
! Actual array dimensions
!
  ! Number of profiles
INTEGER :: n_profile

! Number of layers seen by radiation
INTEGER :: n_layer

! Dimensions for input UKCA arrays
INTEGER :: ukca_dim1
INTEGER :: ukca_dim2

!
! With intent in
!
  ! model switch to include an extra top layer in the radiation scheme
LOGICAL :: l_extra_top

! gathering array
INTEGER :: i_gather(npd_field)

! UKCA component mass mixing ratios
REAL :: ukca_mmr(ukca_dim1, ukca_dim2, n_ukca_cpnt)

! UKCA component volumes
REAL :: ukca_cvl(ukca_dim1, ukca_dim2, n_ukca_cpnt)

! UKCA modal dry and wet diameters
REAL :: ukca_dry(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL :: ukca_wet(ukca_dim1, ukca_dim2, n_ukca_mode)

! UKCA modal densities
REAL :: ukca_rho(ukca_dim1, ukca_dim2, n_ukca_mode)

! UKCA modal volumes
REAL :: ukca_vol(ukca_dim1, ukca_dim2, n_ukca_mode)

! UKCA modal volume of water
REAL :: ukca_wtv(ukca_dim1, ukca_dim2, n_ukca_mode)

! UKCA modal number concentrations
REAL :: ukca_nbr(ukca_dim1, ukca_dim2, n_ukca_mode)

! Level of tropopause
INTEGER :: trindx(npd_field)

!
! With intent out
!
  ! UKCA component mass mixing ratios on radiation code domain
REAL :: ukcaaer_mix_ratio(npd_profile, npd_layer, n_ukca_cpnt)

! UKCA component volumes on radiation code domain
REAL :: ukcaaer_comp_vol(npd_profile, npd_layer, n_ukca_cpnt)

! UKCA modal dry and wet diameters on radiation code domain
REAL :: ukcaaer_dry_diam(npd_profile, npd_layer, n_ukca_mode)
REAL :: ukcaaer_wet_diam(npd_profile, npd_layer, n_ukca_mode)

! UKCA modal densities on radiation code domain
REAL :: ukcaaer_modal_rho(npd_profile, npd_layer, n_ukca_mode)

! UKCA modal volumes
REAL :: ukcaaer_modal_vol(npd_profile, npd_layer, n_ukca_mode)

! UKCA modal volumes of water
REAL :: ukcaaer_modal_wtv(npd_profile, npd_layer, n_ukca_mode)

! UKCA modal number concentrations
REAL :: ukcaaer_modal_nbr(npd_profile, npd_layer, n_ukca_mode)

! Level of tropopause for use with UKCA-MODE aerosols
INTEGER :: trindxrad(npd_profile)

!
! Local variables
!
INTEGER :: i, &
        j, &
        l
INTEGER :: i_top_copy
INTEGER :: lg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_SET_AEROSOL_FIELD'

!
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (l_extra_top) THEN
  i_top_copy = 2
ELSE
  i_top_copy = 1
END IF

DO j = 1, n_ukca_mode

  ! Here would a good place to switch off the direct effect
  ! of a given mode depending on a model switch.

  DO i = i_top_copy, n_layer
    DO l = 1, n_profile
      lg = i_gather(l)
      ukcaaer_dry_diam(l, i, j)  = ukca_dry(lg, n_layer+1-i, j)
      ukcaaer_wet_diam(l, i, j)  = ukca_wet(lg, n_layer+1-i, j)
      ukcaaer_modal_rho(l, i, j) = ukca_rho(lg, n_layer+1-i, j)
      ukcaaer_modal_vol(l, i, j) = ukca_vol(lg, n_layer+1-i, j)
      ukcaaer_modal_wtv(l, i, j) = ukca_wtv(lg, n_layer+1-i, j)
      ukcaaer_modal_nbr(l, i, j) = ukca_nbr(lg, n_layer+1-i, j)
    END DO ! l
  END DO ! i

  ! If using an extra top layer extrapolate the modal properties
  ! from the adjacent layer.
  IF (l_extra_top) THEN
    DO l = 1, n_profile
      ukcaaer_dry_diam(l, 1, j)  = ukcaaer_dry_diam(l, 2, j)
      ukcaaer_wet_diam(l, 1, j)  = ukcaaer_wet_diam(l, 2, j)
      ukcaaer_modal_rho(l, 1, j) = ukcaaer_modal_rho(l, 2, j)
      ukcaaer_modal_vol(l, 1, j) = ukcaaer_modal_vol(l, 2, j)
      ukcaaer_modal_wtv(l, 1, j) = ukcaaer_modal_wtv(l, 2, j)
      ukcaaer_modal_nbr(l, 1, j) = ukcaaer_modal_nbr(l, 2, j)
    END DO ! l
  END IF

END DO ! j

DO j = 1, n_ukca_cpnt

  ! Here would be a good place to switch off the direct effect
  ! of a given component depending on a model switch.

  DO i = i_top_copy, n_layer
    DO l = 1, n_profile
      lg = i_gather(l)
      ukcaaer_mix_ratio(l, i, j) = ukca_mmr(lg, n_layer+1-i, j)
      ukcaaer_comp_vol(l, i, j)  = ukca_cvl(lg, n_layer+1-i, j)
    END DO ! l
  END DO ! i

  ! If using an extra top layer extrapolate the mixing ratio
  ! from the adjacent layer.
  IF (l_extra_top) THEN
    DO l = 1, n_profile
      ukcaaer_mix_ratio(l, 1, j) = ukcaaer_mix_ratio(l, 2, j)
      ukcaaer_comp_vol(l, 1, j)  = ukcaaer_comp_vol(l, 2, j)
    END DO ! l
  END IF

END DO ! j

!
! Convert the numbering of model level tropopause from the
! model convention to the radiation scheme convention.
!
DO l = 1, n_profile

  lg = i_gather(l)
  trindxrad(l) = n_layer + 1 - trindx(lg)

END DO ! l

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_set_aerosol_field
END MODULE ukca_radaer_set_aerosol_field_mod
