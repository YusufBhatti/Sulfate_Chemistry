! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE h_vers_mod

USE  version_mod,  ONLY: nsectp
USE  submodel_mod, ONLY: N_Internal_Model_Max

IMPLICIT NONE

!  Subroutine h_vers_mod
!
! Description:
!   Identify which 'version' masks are available for each section of UM code.
!   This is done by using the run_time code section variables i_<sec>_vn
!
! Method:
!   Treat each code section in turn to provide a version number for use
!   by tstmsk, based upon the numbers used in previous IFDEFs.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

INTEGER   ::   h_vers  (N_Internal_Model_Max,0:nsectp)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='H_VERS_MOD'

CONTAINS

SUBROUTINE h_vers_init()

USE bl_option_mod,        ONLY: i_bl_vn, i_bl_vn_0,                       &
                                i_bl_vn_1a, i_bl_vn_9b, i_bl_vn_9c
USE cv_run_mod,           ONLY: i_convection_vn, i_convection_vn_5a,      &
                                i_convection_vn_6a
USE dust_parameters_mod,  ONLY: i_dust
USE eng_corr_inputs_mod,  ONLY: l_emcorr
USE g_wave_input_mod,     ONLY: i_gwd_vn, i_gwd_vn_4a , i_gwd_vn_5a
USE iau_mod,              ONLY: l_iau
USE jules_vegetation_mod, ONLY: i_veg_vn
USE missing_data_mod,     ONLY: imdi
USE nudging_input_mod,    ONLY: l_nudging
USE river_inputs_mod
USE ukca_option_mod,      ONLY: l_ukca, l_ukca_mode, l_ukca_chem_plev,    &
                                l_ukca_asad_plev 
USE glomap_clim_option_mod,      ONLY: l_glomap_mode_clim

USE model_domain_mod, ONLY: output_grid_stagger,                          &
                            FH_GridStagger_C

USE jules_hydrology_mod,  ONLY: l_hydrology
USE run_aerosol_mod,      ONLY: l_soot, l_biomass, l_ocff, l_nitrate,     &
                                l_sulpc_so2
USE electric_inputs_mod,  ONLY: l_use_electric
USE mphys_inputs_mod,     ONLY: l_casim
USE dynamics_testing_mod, ONLY: l_idealised_data

USE um_stashcode_mod,     ONLY: stashcode_pws_sec

USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook

IMPLICIT NONE

! Local variables
LOGICAL :: l_aero

INTEGER :: im,IS


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='H_VERS_INIT'


! ----------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------
! Hard-wire to atmos sub model atmos_im=1
im=1

DO IS = 1,nsectp
  h_vers(im,IS)=0
END DO

! The following provides the 'magic numbers' of the versions
! available to each code section.
! Where possible these are set via section runtime variables

! Radiation A01, A02, only version 3 available
h_vers(im,1)=3
h_vers(im,2)=3

! BL A03 has multiple versions, 0A, 1A, 9B and 9C
! problematic as two 9 versions so cannot purely use runtime inputs.
SELECT CASE (i_bl_vn)
CASE (i_bl_vn_0)
  h_vers(im,3)=i_bl_vn_0
CASE (i_bl_vn_1a)
  h_vers(im,3)=i_bl_vn_1a
CASE (i_bl_vn_9b)
  h_vers(im,3)=9
CASE (i_bl_vn_9c)
  h_vers(im,3)=9
END SELECT

! LSP A04 has now two options, version 4 (CASIM) or 
! version 3 (Wilson and Ballard, 1999).
IF (l_casim) THEN
  h_vers(im,4)=4
ELSE
  h_vers(im,4)=3
END IF

! Convection A05 has multiple versions 5A, 6A
h_vers(im,5)=i_convection_vn

! GWD A06 has multiple options  4A 5A
IF (i_gwd_vn /= imdi) THEN
  h_vers(im,6)=i_gwd_vn
END IF

! land surface A08
IF (l_hydrology) THEN
  h_vers(im,8)=8
END IF

! large scale cloud A09
h_vers(im,9)=2

! dynamics solver only valid for ND, 2A or 2B
! EG and section 10 diags are invalid.
IF (output_grid_stagger == FH_GridStagger_C) THEN
  h_vers(im,10)=2    ! no diags are available for ND
ELSE
  h_vers(im,10)=4    ! EG diags are available, see stashmaster.
END IF

! Tracer advection ther are no diags so can set to anything.
! h_vers(im,11)=0

! dynamics advection diags. ND
! Unsure which are valid for EG will need to be revisited.
h_vers(im,12)=2

! Diffusion and Filtering  too complicated at present.
! Very complex assume 2A,2B selected
h_vers(im,13)=2

! Energy Correction A14 single version
IF (l_emcorr) THEN
  h_vers(im,14)=1
END IF

! Dynamics diagnostics A15 assume always available
h_vers(im,15)=1

! Physics diagnostics A16 assume always available
h_vers(im,16)=1

! Aerosols A17 very complex to follow if diags available or not.
l_aero = (l_soot .OR. l_biomass .OR.                         &
         (i_dust /= imdi .AND. i_dust /=0 ) .OR.             &
         l_sulpc_so2 .OR. l_ocff .OR. l_nitrate)

IF (l_aero) THEN
  h_vers(im,17)=2
END IF

! Data assimilation A18 single version
IF (l_iau) THEN
  h_vers(im,18)=2
END IF

! Veg distribution A19
h_vers(im,19)=i_veg_vn

! PWS diagnostics (migrated from FieldCalc)
h_vers(im,stashcode_pws_sec)=1

! Thunderstorm electrification
IF (l_use_electric) THEN
  h_vers(im,21)=1
END IF

! River routing
IF (l_rivers) THEN
  h_vers(im,26)=i_river_vn
END IF

! Climate diagnostics assume always wanted
h_vers(im,30)=1

! LBC INPUT A31
h_vers(im,31)=1

! Free Tracers  set as available should someone use them.
h_vers(im,33)=1

! UKCA A34
IF (l_ukca) THEN
  h_vers(im,34)=1
END IF

! Stochastic physics
h_vers(im,35)=1

! Atmos tracer LBCs
h_vers(im,36)=1

! UKCA LBCs
IF (l_ukca) THEN
  h_vers(im,37)=1
END IF

! UKCA Aerosols
IF (l_ukca_mode) THEN
  h_vers(im,38)=1
END IF

! Nudging scheme
IF (l_nudging) THEN
  h_vers(im,39)=1
END IF

! UKCA ASAD diags
IF (l_ukca) THEN
  h_vers(im,50)=1
END IF

! UKCA Chem on pressure levels
IF (l_ukca_chem_plev) THEN
  h_vers(im,51)=1
END IF

! UKCA ASAD on pressure levels
IF (l_ukca_asad_plev) THEN
  h_vers(im,52)=1
END IF

! Idealised model diagnostics
IF (l_idealised_data) THEN
  h_vers(im,53)=1
END IF

! GLOMAP-mode aerosol climatologies A54
IF (l_glomap_mode_clim) THEN
   h_vers(im,54)=1
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE h_vers_init
END MODULE  h_vers_mod

