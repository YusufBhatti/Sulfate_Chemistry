! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_dimen_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_DIMEN_MOD'
CONTAINS

SUBROUTINE set_dimen(control, dimen, spectrum,                                 &
  nd_field, n_points, n_layer, cloud_levels, n_ukca_mode, l_easyaerosol_rad,   &
  nd_field_flux_diag, nd_field_rad_diag)

! Subroutine to set dimensions for the radiation code.
!
! Purpose:
!   To set the dimensions of arrays for use in the radiation code
!   depending on the options chosen.

USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData
USE mcica_mod,    ONLY: tot_subcol_gen, subcol_need
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE jules_sea_seaice_mod, ONLY: l_ctile

IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
  TYPE(StrDim),       INTENT(INOUT) :: dimen

! Spectral data:
  TYPE (StrSpecData), INTENT(IN)    :: spectrum

INTEGER, INTENT(IN) :: nd_field
!   Horizontal size of input fields
INTEGER, INTENT(IN) :: n_points
!   Number of grid points to operate on
INTEGER, INTENT(IN) :: n_layer
!   Number of layers for radiation
INTEGER, INTENT(IN) :: cloud_levels
!   Number of potentially cloudy layers in the model
INTEGER, INTENT(IN) :: n_ukca_mode
!   Number of UKCA aerosol modes
LOGICAL, INTENT(IN) :: l_easyaerosol_rad
!   EasyAerosol interaction with radiation
INTEGER, INTENT(OUT) :: nd_field_flux_diag
!   Size to be allocated for flux diagnostics
INTEGER, INTENT(OUT) :: nd_field_rad_diag
!   Size to be allocated for radiance diagnostics

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_DIMEN'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid dimensions
! (On some architectures an odd size is preferred to avoid memory
! bank conflicts).
dimen%nd_profile = 2*(n_points/2)+1
dimen%nd_layer = n_layer


! Cloud
dimen%nd_cloud_type      = 4
dimen%nd_cloud_component = 4
dimen%nd_cloud_representation = 4
SELECT CASE(control%i_cloud)
CASE (ip_cloud_column_max)
  dimen%nd_column = 3 * cloud_levels + 2
CASE DEFAULT
  dimen%nd_column = 1
END SELECT
dimen%nd_subcol_gen = tot_subcol_gen
dimen%nd_subcol_req = subcol_need

dimen%id_cloud_top = dimen%nd_layer + 1 - cloud_levels
IF (control%l_cloud) THEN
  dimen%nd_layer_clr = dimen%id_cloud_top - 1
ELSE
  dimen%nd_layer_clr = dimen%nd_layer
END IF


! Aerosol
IF (l_easyaerosol_rad) THEN
  dimen%nd_aerosol_mode = MAX(1,n_ukca_mode + 1)
ELSE
  dimen%nd_aerosol_mode = MAX(1,n_ukca_mode)
END IF

! Arrays for prescribed optical properties (not used here)
dimen%nd_profile_aerosol_prsc   = 1
dimen%nd_profile_cloud_prsc     = 1
dimen%nd_opt_level_aerosol_prsc = 1
dimen%nd_opt_level_cloud_prsc   = 1


! Tiled surface.
dimen%nd_tile_type  = 3
dimen%nd_point_tile = MAX(1,n_points)
IF (l_ctile) THEN
  dimen%nd_tile     = 3
ELSE
  dimen%nd_tile     = 2
END IF


! Solver dependent dimensions
IF (control%i_angular_integration == ip_two_stream) THEN

  dimen%nd_viewing_level    = 1
  dimen%nd_radiance_profile = 1
  dimen%nd_j_profile        = 1
  dimen%nd_direction        = 1
  dimen%nd_brdf_basis_fnc   = 2
  dimen%nd_brdf_trunc       = 1
  nd_field_flux_diag        = nd_field
  dimen%nd_flux_profile     = dimen%nd_profile
  nd_field_rad_diag         = 1
  dimen%nd_channel          = 1

  dimen%nd_2sg_profile      = dimen%nd_profile
  dimen%nd_source_coeff     = 2
  dimen%nd_max_order        = 1
  dimen%nd_sph_coeff        = 1

  SELECT CASE(control%i_solver)

  CASE(ip_solver_mix_app_scat, ip_solver_mix_direct, ip_solver_mix_direct_hogan)
    dimen%nd_overlap_coeff=8
    dimen%nd_region=2

  CASE(ip_solver_triple_app_scat, ip_solver_triple, ip_solver_triple_hogan)
    dimen%nd_overlap_coeff=18
    dimen%nd_region=3

  CASE DEFAULT
    dimen%nd_overlap_coeff=1
    dimen%nd_region=2

  END SELECT

ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

  dimen%nd_direction        = 1
  dimen%nd_brdf_basis_fnc   = 2
  dimen%nd_radiance_profile = dimen%nd_profile
  dimen%nd_j_profile        = 1

  IF (control%i_sph_mode == ip_sph_mode_flux) THEN
    dimen%nd_viewing_level  = n_layer+1
    dimen%nd_brdf_trunc     = 1
    dimen%nd_flux_profile   = dimen%nd_profile
    nd_field_flux_diag      = nd_field
    nd_field_rad_diag       = 1
    
    ! When calculating fluxes all spectral bands are mapped into a
    ! single output channel.
    dimen%nd_channel        = 1
  ELSE
    dimen%nd_viewing_level  = 1
    dimen%nd_brdf_trunc     = 1
    dimen%nd_flux_profile   = 1
    nd_field_flux_diag      = 1
    nd_field_rad_diag       = nd_field
    
    ! We assume for now that each band in the spectral file
    ! corresponds to a particular channel.
    dimen%nd_channel        = spectrum%basic%n_band
  END IF

  dimen%nd_2sg_profile      = 1
  dimen%nd_source_coeff     = 1
  dimen%nd_overlap_coeff    = 1
  dimen%nd_region           = 1
  dimen%nd_max_order        = control%ls_global_trunc + 2
  IF (control%i_truncation == ip_trunc_triangular) THEN
    dimen%nd_sph_coeff =                                                       &
      (control%ls_global_trunc+3)*(control%ls_global_trunc+4)/2
  ELSE IF (control%i_truncation == ip_trunc_azim_sym) THEN
    dimen%nd_sph_coeff = control%ls_global_trunc+2
  ELSE
    cmessage = 'Illegal truncation'
    ierr=i_err_fatal
    GO TO 9999
  END IF

END IF


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_dimen
END MODULE set_dimen_mod
