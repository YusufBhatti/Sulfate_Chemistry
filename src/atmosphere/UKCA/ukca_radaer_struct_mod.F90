! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module file ukca_radaer_struct_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

! UKCA_RADAER:
!
! Defines maximum dimensions
! Defines type ukca_radaer_struct, the structure used by UKCA_RADAER

! Contained subroutines:
!      allocate_radaer_struct
!      deallocate_radaer_struct

MODULE ukca_radaer_struct_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_STRUCT_MOD'

! Internal IDs for the types of UKCA aerosol modes.
!
INTEGER, PARAMETER :: ip_ukca_mode_nucleation = 0
INTEGER, PARAMETER :: ip_ukca_mode_aitken     = 1
INTEGER, PARAMETER :: ip_ukca_mode_accum      = 2
INTEGER, PARAMETER :: ip_ukca_mode_coarse     = 3

!
! Internal IDs for the types of UKCA aerosol components.
!
! Note: When adding an aerosol component, also increase
!       npd_ukca_maxcomptype in module ukca_radaer_precalc_mod and
!       update the pre-computed file read by ukca_radaer_read_precalc.
!       Any new aerosol component should be added before ip_ukca_water.
!
! Water is not an aerosol species, but is included here
! as it behaves like one in ukca_radaer_band_average().
! However, no aerosol component should be of type ip_ukca_water.
!
INTEGER, PARAMETER :: ip_ukca_sulphate      = 1
INTEGER, PARAMETER :: ip_ukca_blackcarbon   = 2
INTEGER, PARAMETER :: ip_ukca_organiccarbon = 3
INTEGER, PARAMETER :: ip_ukca_seasalt       = 4
INTEGER, PARAMETER :: ip_ukca_dust          = 5
INTEGER, PARAMETER :: ip_ukca_secondorganic = 6
INTEGER, PARAMETER :: ip_ukca_nitrate       = 7
INTEGER, PARAMETER :: ip_ukca_h2so4         = 8
INTEGER, PARAMETER :: ip_ukca_water         = 9

INTEGER, SAVE :: npd_ukca_cpnt       ! nmodes*ncp

! Main structure holding all the variables needed for
! interacting UKCA aerosols with radiation.

TYPE ukca_radaer_struct

  !
  ! Information about UKCA aerosol modes
  !

  ! Actual number of modes, with a default value for
  ! minimising array dimensions.
  INTEGER :: n_mode = 1

  ! Type of mode (i.e. nucleation, Aitken, accum, or coarse)
  INTEGER, ALLOCATABLE :: i_mode_type(:)

  ! Solubility of mode (soluble if true)
  LOGICAL, ALLOCATABLE :: l_soluble(:)

  ! Lower and upper limits on the geometric mean diameter (m)
  ! in each mode
  REAL, ALLOCATABLE :: d0low(:)
  REAL, ALLOCATABLE :: d0up(:)

  ! Geometric standard deviation in this mode
  REAL, ALLOCATABLE :: sigma(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal dry diameter.
  INTEGER, ALLOCATABLE :: stashcode_dry(:)
  INTEGER, ALLOCATABLE :: d1_address_dry(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_dry(:)
  INTEGER, ALLOCATABLE :: d1_length_dry(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal wet diameter.
  INTEGER, ALLOCATABLE :: stashcode_wet(:)
  INTEGER, ALLOCATABLE :: d1_address_wet(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_wet(:)
  INTEGER, ALLOCATABLE :: d1_length_wet(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to modal density.
  INTEGER, ALLOCATABLE :: stashcode_rho(:)
  INTEGER, ALLOCATABLE :: d1_address_rho(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_rho(:)
  INTEGER, ALLOCATABLE :: d1_length_rho(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to water volume in each mode.
  INTEGER, ALLOCATABLE :: stashcode_wtv(:)
  INTEGER, ALLOCATABLE :: d1_address_wtv(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_wtv(:)
  INTEGER, ALLOCATABLE :: d1_length_wtv(:)

  ! STASH code, D1 address, number of levels, total length
  ! and halo type of the D1 fields corresponding to
  ! modal number concentrations.
  INTEGER, ALLOCATABLE :: stashcode_nbr(:)
  INTEGER, ALLOCATABLE :: d1_address_nbr(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_nbr(:)
  INTEGER, ALLOCATABLE :: d1_length_nbr(:)
  INTEGER, ALLOCATABLE :: d1_halo_type_nbr(:)

  ! Number of components in each mode and index of each
  ! component in array ukca_cpnt_info
  INTEGER, ALLOCATABLE :: n_cpnt_in_mode(:)
  INTEGER, ALLOCATABLE :: i_cpnt_index(:,:)

  ! Modal diameter of the dry aerosol (m)
  REAL, POINTER :: dry_diam(:, :, :, :)

  ! Modal diameter of the wet aerosol (m)
  REAL, POINTER :: wet_diam(:, :, :, :)

  ! Modal densities (kg/m3)
  REAL, POINTER :: modal_rho(:, :, :, :) 

  ! Modal volumes (including water for soluble modes)
  REAL, POINTER :: modal_vol(:, :, :, :) 

  ! Fractional volume of water in each mode
  REAL, POINTER :: modal_wtv(:, :, :, :)

  ! Modal number concentrations
  REAL, POINTER :: modal_nbr(:, :, :, :) 

  !
  ! Information about UKCA aerosol components
  !

  ! Actual number of components, with a default value for
  ! minimising array dimensions.
  INTEGER :: n_cpnt = 1

  ! Type of component (e.g. sulphate)
  INTEGER, ALLOCATABLE :: i_cpnt_type(:)

  ! Mass density of each component (kg/m3)
  REAL, ALLOCATABLE :: density(:)

  ! Array index of the mode this component belongs to
  INTEGER, ALLOCATABLE :: i_mode(:)

  ! STASH code, D1 address, number of levels, total length
  ! and halo type of the D1 fields corresponding to the
  ! mass-mixing ratio of each component.
  INTEGER, ALLOCATABLE :: stashcode_mmr(:)
  INTEGER, ALLOCATABLE :: d1_address_mmr(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_mmr(:)
  INTEGER, ALLOCATABLE :: d1_length_mmr(:)
  INTEGER, ALLOCATABLE :: d1_halo_type_mmr(:)

  ! STASH code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to the fractional volume of
  ! each component.
  INTEGER, ALLOCATABLE :: stashcode_cvl(:)
  INTEGER, ALLOCATABLE :: d1_address_cvl(:)
  INTEGER, ALLOCATABLE :: d1_nlevs_cvl(:)
  INTEGER, ALLOCATABLE :: d1_length_cvl(:)

  ! Component mass-mixing ratio (kg/kg)
  REAL, POINTER :: mix_ratio(:, :, :, :)

  ! Component volumes
  REAL, POINTER :: comp_vol(:, :, :, :) 

  ! Switch: if true, use sulphuric acid optical properties in the
  ! stratosphere, instead of ammonium sulphate (if false).
  ! Has the same default as run_ukca:l_ukca_radaer_sustrat or 
  ! run_glomap_aeroclim:l_glomap_clim_radaer_sustrat
  ! from which is it assigned.
  LOGICAL :: l_sustrat = .FALSE.


END TYPE ukca_radaer_struct


CONTAINS

! #############################################################################

SUBROUTINE allocate_radaer_struct(ukca_radaer)

! To allocate arrays in the structure ukca_radaer using the number of modes
! and number of components configured in the GLOMAP setup routine.

USE ukca_mode_setup, ONLY: nmodes, ncp
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook
IMPLICIT NONE

TYPE(ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_RADAER_STRUCT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

npd_ukca_cpnt = ncp * nmodes
ukca_radaer%n_mode = nmodes
ukca_radaer%n_cpnt = ncp

IF (.NOT. ALLOCATED(ukca_radaer%i_mode_type))               &
      ALLOCATE(ukca_radaer%i_mode_type(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%l_soluble))                 &
      ALLOCATE(ukca_radaer%l_soluble(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d0low))                     &
      ALLOCATE(ukca_radaer%d0low(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d0up))                      &
      ALLOCATE(ukca_radaer%d0up(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%sigma))                     &
      ALLOCATE(ukca_radaer%sigma(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_dry))             &
      ALLOCATE(ukca_radaer%stashcode_dry(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_dry))            &
      ALLOCATE(ukca_radaer%d1_address_dry(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_dry))              &
      ALLOCATE(ukca_radaer%d1_nlevs_dry(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_dry))             &
      ALLOCATE(ukca_radaer%d1_length_dry(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_wet))             &
      ALLOCATE(ukca_radaer%stashcode_wet(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_wet))            &
      ALLOCATE(ukca_radaer%d1_address_wet(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_wet))              &
      ALLOCATE(ukca_radaer%d1_nlevs_wet(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_wet))             &
      ALLOCATE(ukca_radaer%d1_length_wet(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_rho))             &
      ALLOCATE(ukca_radaer%stashcode_rho(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_rho))            &
      ALLOCATE(ukca_radaer%d1_address_rho(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_rho))              &
      ALLOCATE(ukca_radaer%d1_nlevs_rho(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_rho))             &
      ALLOCATE(ukca_radaer%d1_length_rho(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_wet))             &
      ALLOCATE(ukca_radaer%d1_length_wet(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_wtv))             &
      ALLOCATE(ukca_radaer%stashcode_wtv(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_wtv))            &
      ALLOCATE(ukca_radaer%d1_address_wtv(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_wtv))              &
      ALLOCATE(ukca_radaer%d1_nlevs_wtv(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_wtv))             &
      ALLOCATE(ukca_radaer%d1_length_wtv(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_nbr))             &
      ALLOCATE(ukca_radaer%stashcode_nbr(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_nbr))            &
      ALLOCATE(ukca_radaer%d1_address_nbr(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_nbr))              &
      ALLOCATE(ukca_radaer%d1_nlevs_nbr(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_nbr))             &
      ALLOCATE(ukca_radaer%d1_length_nbr(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%d1_halo_type_nbr))          &
      ALLOCATE(ukca_radaer%d1_halo_type_nbr(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%n_cpnt_in_mode))            &
      ALLOCATE(ukca_radaer%n_cpnt_in_mode(nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%i_cpnt_index))              &
      ALLOCATE(ukca_radaer%i_cpnt_index(ncp,nmodes))
IF (.NOT. ALLOCATED(ukca_radaer%i_cpnt_type))               &
      ALLOCATE(ukca_radaer%i_cpnt_type(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%density))                   &
      ALLOCATE(ukca_radaer%density(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%i_mode))                    &
      ALLOCATE(ukca_radaer%i_mode(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_mmr))             &
      ALLOCATE(ukca_radaer%stashcode_mmr(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_mmr))            &
      ALLOCATE(ukca_radaer%d1_address_mmr(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_mmr))              &
      ALLOCATE(ukca_radaer%d1_nlevs_mmr(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_mmr))             &
      ALLOCATE(ukca_radaer%d1_length_mmr(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_halo_type_mmr))          &
      ALLOCATE(ukca_radaer%d1_halo_type_mmr(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%stashcode_cvl))             &
      ALLOCATE(ukca_radaer%stashcode_cvl(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_address_cvl))            &
      ALLOCATE(ukca_radaer%d1_address_cvl(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_nlevs_cvl))              &
      ALLOCATE(ukca_radaer%d1_nlevs_cvl(npd_ukca_cpnt))
IF (.NOT. ALLOCATED(ukca_radaer%d1_length_cvl))             &
      ALLOCATE(ukca_radaer%d1_length_cvl(npd_ukca_cpnt))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE allocate_radaer_struct

! ##########################################################################

SUBROUTINE deallocate_radaer_struct(ukca_radaer)

! To deallocate arrays in the ukca_radaer_struct.

USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook
IMPLICIT NONE

TYPE(ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOCATE_RADAER_STRUCT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(ukca_radaer%i_mode_type))               &
      DEALLOCATE(ukca_radaer%i_mode_type)
IF (ALLOCATED(ukca_radaer%l_soluble))                 &
      DEALLOCATE(ukca_radaer%l_soluble)
IF (ALLOCATED(ukca_radaer%d0low))                     &
      DEALLOCATE(ukca_radaer%d0low)
IF (ALLOCATED(ukca_radaer%d0up))                      &
      DEALLOCATE(ukca_radaer%d0up)
IF (ALLOCATED(ukca_radaer%sigma))                     &
      DEALLOCATE(ukca_radaer%sigma)
IF (ALLOCATED(ukca_radaer%stashcode_dry))             &
      DEALLOCATE(ukca_radaer%stashcode_dry)
IF (ALLOCATED(ukca_radaer%d1_address_dry))            &
      DEALLOCATE(ukca_radaer%d1_address_dry)
IF (ALLOCATED(ukca_radaer%d1_nlevs_dry))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_dry)
IF (ALLOCATED(ukca_radaer%d1_length_dry))             &
      DEALLOCATE(ukca_radaer%d1_length_dry)
IF (ALLOCATED(ukca_radaer%stashcode_wet))             &
      DEALLOCATE(ukca_radaer%stashcode_wet)
IF (ALLOCATED(ukca_radaer%d1_address_wet))            &
      DEALLOCATE(ukca_radaer%d1_address_wet)
IF (ALLOCATED(ukca_radaer%d1_nlevs_wet))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_wet)
IF (ALLOCATED(ukca_radaer%d1_length_wet))             &
      DEALLOCATE(ukca_radaer%d1_length_wet)
IF (ALLOCATED(ukca_radaer%stashcode_rho))             &
      DEALLOCATE(ukca_radaer%stashcode_rho)
IF (ALLOCATED(ukca_radaer%d1_address_rho))            &
      DEALLOCATE(ukca_radaer%d1_address_rho)
IF (ALLOCATED(ukca_radaer%d1_nlevs_rho))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_rho)
IF (ALLOCATED(ukca_radaer%d1_length_rho))             &
      DEALLOCATE(ukca_radaer%d1_length_rho)
IF (ALLOCATED(ukca_radaer%d1_length_wet))             &
      DEALLOCATE(ukca_radaer%d1_length_wet)
IF (ALLOCATED(ukca_radaer%stashcode_wtv))             &
      DEALLOCATE(ukca_radaer%stashcode_wtv)
IF (ALLOCATED(ukca_radaer%d1_address_wtv))            &
      DEALLOCATE(ukca_radaer%d1_address_wtv)
IF (ALLOCATED(ukca_radaer%d1_nlevs_wtv))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_wtv)
IF (ALLOCATED(ukca_radaer%d1_length_wtv))             &
      DEALLOCATE(ukca_radaer%d1_length_wtv)
IF (ALLOCATED(ukca_radaer%stashcode_nbr))             &
      DEALLOCATE(ukca_radaer%stashcode_nbr)
IF (ALLOCATED(ukca_radaer%d1_address_nbr))            &
      DEALLOCATE(ukca_radaer%d1_address_nbr)
IF (ALLOCATED(ukca_radaer%d1_nlevs_nbr))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_nbr)
IF (ALLOCATED(ukca_radaer%d1_length_nbr))             &
      DEALLOCATE(ukca_radaer%d1_length_nbr)
IF (ALLOCATED(ukca_radaer%d1_halo_type_nbr))          &
      DEALLOCATE(ukca_radaer%d1_halo_type_nbr)
IF (ALLOCATED(ukca_radaer%n_cpnt_in_mode))            &
      DEALLOCATE(ukca_radaer%n_cpnt_in_mode)
IF (ALLOCATED(ukca_radaer%i_cpnt_index))              &
      DEALLOCATE(ukca_radaer%i_cpnt_index)
IF (ALLOCATED(ukca_radaer%i_cpnt_type))               &
      DEALLOCATE(ukca_radaer%i_cpnt_type)
IF (ALLOCATED(ukca_radaer%density))                   &
      DEALLOCATE(ukca_radaer%density)
IF (ALLOCATED(ukca_radaer%i_mode))                    &
      DEALLOCATE(ukca_radaer%i_mode)
IF (ALLOCATED(ukca_radaer%stashcode_mmr))             &
      DEALLOCATE(ukca_radaer%stashcode_mmr)
IF (ALLOCATED(ukca_radaer%d1_address_mmr))            &
      DEALLOCATE(ukca_radaer%d1_address_mmr)
IF (ALLOCATED(ukca_radaer%d1_nlevs_mmr))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_mmr)
IF (ALLOCATED(ukca_radaer%d1_length_mmr))             &
      DEALLOCATE(ukca_radaer%d1_length_mmr)
IF (ALLOCATED(ukca_radaer%d1_halo_type_mmr))          &
      DEALLOCATE(ukca_radaer%d1_halo_type_mmr)
IF (ALLOCATED(ukca_radaer%stashcode_cvl))             &
      DEALLOCATE(ukca_radaer%stashcode_cvl)
IF (ALLOCATED(ukca_radaer%d1_address_cvl))            &
      DEALLOCATE(ukca_radaer%d1_address_cvl)
IF (ALLOCATED(ukca_radaer%d1_nlevs_cvl))              &
      DEALLOCATE(ukca_radaer%d1_nlevs_cvl)
IF (ALLOCATED(ukca_radaer%d1_length_cvl))             &
      DEALLOCATE(ukca_radaer%d1_length_cvl)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE deallocate_radaer_struct

END MODULE ukca_radaer_struct_mod
