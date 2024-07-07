! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE lw_control_default_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LW_CONTROL_DEFAULT_MOD'
CONTAINS

SUBROUTINE lw_control_default(lw_control)

! Subroutine to set the default values of the control structure.

! Description:
!   Performs suitable default initializations of the controlling 
!   structure for LW radiation.

USE def_control, ONLY: StrCtrl
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE rad_pcf,     ONLY: ip_infra_red, ip_two_stream, ip_elsasser, &
                       ip_solver_homogen_direct, ip_trunc_azim_sym, &
                       ip_sph_direct, ip_sph_mode_rad

IMPLICIT NONE

TYPE (StrCtrl), INTENT(INOUT) :: lw_control
!   The block of controlling options for the code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LW_CONTROL_DEFAULT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Spectral region and bands
lw_control%isolir     = ip_infra_red
lw_control%first_band = 1

! Physical processes
lw_control%l_gas          = .TRUE.
lw_control%l_cloud        = .TRUE.
lw_control%l_drop         = .TRUE.
lw_control%l_ice          = .TRUE.
lw_control%l_aerosol      = .TRUE.
lw_control%l_aerosol_ccn  = .TRUE.

! Properties of clouds
lw_control%l_local_cnv_partition = .TRUE.
lw_control%l_global_cloud_top    = .TRUE.

! Angular integration (including algorithmic options)
lw_control%n_channel              = 1
lw_control%i_angular_integration  = ip_two_stream
lw_control%i_2stream              = ip_elsasser
lw_control%i_solver_clear         = ip_solver_homogen_direct
lw_control%n_order_gauss          = 0
lw_control%i_truncation           = ip_trunc_azim_sym
lw_control%i_sph_algorithm        = ip_sph_direct
lw_control%n_order_phase_solar    = 1
lw_control%ls_global_trunc        = 9
lw_control%ms_min                 = 0
lw_control%ms_max                 = 0
lw_control%ls_brdf_trunc          = 0
lw_control%accuracy_adaptive      = 1.0e-04
lw_control%l_rescale              = .TRUE.
lw_control%l_henyey_greenstein_pf = .TRUE.
lw_control%i_sph_mode             = ip_sph_mode_rad
lw_control%i_solar_src            = 3
lw_control%l_ir_source_quad       = .TRUE.
lw_control%n_esft_red             = 16
lw_control%gpnt_split             = 0.0

! Satellite data
lw_control%sat_hgt      = 0.0
lw_control%sat_lon      = 0.0
lw_control%sat_lat      = 0.0
lw_control%max_view_lon = 0.0
lw_control%min_view_lon = 0.0
lw_control%max_view_lat = 0.0
lw_control%min_view_lat = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE lw_control_default
END MODULE lw_control_default_mod
