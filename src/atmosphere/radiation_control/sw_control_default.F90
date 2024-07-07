! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE sw_control_default_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SW_CONTROL_DEFAULT_MOD'
CONTAINS

SUBROUTINE sw_control_default(sw_control)

! Subroutine to set the default values of the control structure.

! Description:
!   Performs suitable default initializations of the controlling 
!   structure for SW radiation.

USE def_control, ONLY: StrCtrl
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE rad_pcf,     ONLY: ip_solar, ip_two_stream, ip_pifm80, &
                       ip_solver_homogen_direct, ip_trunc_triangular, &
                       ip_sph_direct, ip_sph_mode_rad, ip_scatter_full

IMPLICIT NONE

TYPE (StrCtrl), INTENT(INOUT) :: sw_control
!   The block of controlling options for the code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SW_CONTROL_DEFAULT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Spectral region and bands
sw_control%isolir     = ip_solar
sw_control%first_band = 1

! Physical processes
sw_control%l_microphysics = .TRUE.
sw_control%l_gas          = .TRUE.
sw_control%l_rayleigh     = .TRUE.
sw_control%l_cloud        = .TRUE.
sw_control%l_drop         = .TRUE.
sw_control%l_ice          = .TRUE.
sw_control%l_aerosol      = .TRUE.
sw_control%l_aerosol_ccn  = .TRUE.

! Properties of clouds
sw_control%l_local_cnv_partition = .TRUE.
sw_control%l_global_cloud_top    = .TRUE.

! Angular integration (including algorithmic options)
sw_control%n_channel              = 1
sw_control%i_angular_integration  = ip_two_stream
sw_control%i_2stream              = ip_pifm80
sw_control%i_solver_clear         = ip_solver_homogen_direct
sw_control%n_order_gauss          = 0
sw_control%i_truncation           = ip_trunc_triangular
sw_control%i_sph_algorithm        = ip_sph_direct
sw_control%n_order_phase_solar    = 1
sw_control%ls_global_trunc        = 9
sw_control%ms_min                 = 0
sw_control%ms_max                 = 0
sw_control%ls_brdf_trunc          = 0
sw_control%accuracy_adaptive      = 1.0e-04
sw_control%l_rescale              = .TRUE.
sw_control%l_henyey_greenstein_pf = .TRUE.
sw_control%i_sph_mode             = ip_sph_mode_rad
sw_control%i_solar_src            = 3
sw_control%i_scatter_method       = ip_scatter_full
sw_control%n_esft_red             = 16
sw_control%gpnt_split             = 0.0

! Switches for diagnostic output
sw_control%l_blue_flux_surf = .TRUE.

! Satellite data
sw_control%sat_hgt      = 0.0
sw_control%sat_lon      = 0.0
sw_control%sat_lat      = 0.0
sw_control%max_view_lon = 0.0
sw_control%min_view_lon = 0.0
sw_control%max_view_lat = 0.0
sw_control%min_view_lat = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sw_control_default
END MODULE sw_control_default_mod
