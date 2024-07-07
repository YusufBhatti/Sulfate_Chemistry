! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control.
!
! This module declares 'short-term' temporary logicals used to protect
! science bug fixes that lead to significant alterations in science results.
! It is expected that these logicals will be short lived as the preference
! should be for all configurations to use the corrected code. But
! to maintain short term reproducibility of results across UM versions
! the fixes are protected by logicals until the fixes become the default
! in all model configurations and the logical is retired.
!
! All logicals below should have a review period attached to them for future
! retirement of both the logical and the broken code.
!
! ! ticket #xxxx
! LOGICAL :: fix_me = .FALSE.   ! review again MM YYYY
!
! then add logical to namelist /temp_fixes/
! and add code to subroutine warn_temp_fixes to report when fix is not used
! ie when the logical = .FALSE.
!
! -----------------------------------------------------------
! -----------------------------------------------------------

MODULE science_fixes_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, only: newline
USE iau_mod, only: l_iau
USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE

! ticket  5456
LOGICAL :: l_fix_conserv = .FALSE.      ! Review Again April 2015
! ticket #4982
LOGICAL :: l_rm_neg_par = .FALSE.   ! Review Again Jan 2015
! ticket #5401
LOGICAL :: l_roughnesslength_fix = .FALSE.  ! Review Again June 2015
! ticket #5626
LOGICAL :: l_mphys_gr_out = .FALSE. ! Review Again Oct 2014
! ticket #5841
LOGICAL :: l_fix_arcl_eg_levs = .FALSE. ! Review Again April 2015
! ticket #4481
LOGICAL :: l_fix_drop_settle = .FALSE. ! Review Again May 2015
! ticket #5551
LOGICAL :: l_emis_ssi_full = .FALSE.   ! Review again June 2015
! ticket #6498
LOGICAL :: l_iau_pc2check = .FALSE. ! Review Sep 2016
! ticket #6091
LOGICAL :: l_fail_p_layers_inconsis = .FALSE. ! Review again September 2015
! ticket #6514
LOGICAL :: l_pc2_homog_turb_q_neg = .FALSE. ! Review again Sept 2015
! ticket #6357
LOGICAL :: l_methox_fix = .FALSE.       ! Review Again January 2016

! Below this point tickets correspond to the Science Repository Service
!----------------------------------------------------------------------
! ticket #47
LOGICAL :: l_stph_rhcrit_unbias = .FALSE.     ! Review Again January 2016
! ticket #575
LOGICAL :: l_dtcanfix = .FALSE.               ! Review Again May 2016
! ticket #646
LOGICAL :: l_fix_nh4no3_equilibrium = .FALSE. ! Review Again October 2016
! ticket #911
LOGICAL :: l_eg_damp_height_lid = .FALSE.     ! Review Again Oct 2016
! ticket #367
LOGICAL :: l_fix_mphys_diags_iter = .FALSE.   ! Review Again July 2016
! ticket #430
LOGICAL :: l_rm_hardwire_gas360 = .FALSE.   ! Review again Sep 2016
! ticket #1017
LOGICAL :: l_fix_ctile_orog = .FALSE.   ! Review again Sep 2018
! ticket #1167
LOGICAL :: l_fix_conv_precip_evap = .FALSE.   ! Review Again April 2019
! ticket #1421
LOGICAL :: l_fix_ukca_impscav = .FALSE.   ! Review again Jan 2018
! ticket #1638
LOGICAL :: l_fix_rp_shock_amp = .FALSE.   ! Review again Jan 2018
! ticket #1729
LOGICAL :: l_fix_ustar_dust = .FALSE.     ! Review again Jan 2018
! ticket #2077
LOGICAL :: l_fix_dyndiag = .FALSE.        ! Review again Sep 2018
! ticket #2556
LOGICAL :: l_fix_riming = .FALSE.         ! Review in Dec 2018
! ticket #3076
LOGICAL :: l_fix_ccb_cct = .FALSE.        ! Review in July 2017
! ticket #3080  (also JULES ticket #547)
LOGICAL :: l_fix_alb_ice_thick = .FALSE.  ! Review in Nov 2018
! ticket #3005
LOGICAL :: l_fix_zh = .FALSE.             ! Review in Dec 2018
! ticket #2405
LOGICAL :: l_fix_nacl_density = .FALSE.   ! Review in Dec 2018
! ticket #2710
LOGICAL :: l_fix_iau_rim_density = .FALSE. ! Review in Sep 2019
! ticket #3011
LOGICAL :: l_fix_albsnow_ts = .FALSE.     ! Review in Dec 2020
! ticket #3011
LOGICAL :: l_fix_rcf_mlsnow_icefreemax = .FALSE. ! Review in Dec 2020
! ticket #3681
LOGICAL :: l_fix_conv_diags_var = .FALSE. ! Review in Jan 2020
! ticket #2545
LOGICAL :: l_fix_lsp_incs_to_spt = .FALSE. ! Review in Feb 2019
! ticket #2070
LOGICAL :: l_fix_ec_gen_hgt = .FALSE.     ! Review in Dec 2018
! ticket #1250 (review in ticket #3997)
LOGICAL :: l_fix_improve_drydep = .FALSE. ! Review in Jan 2019
! ticket #4038
LOGICAL :: l_fix_wind_snow = .FALSE.      ! Review in May 2019

NAMELIST/temp_fixes/l_roughnesslength_fix,                    &
        l_rm_neg_par,                                         &
        l_mphys_gr_out, l_fix_conserv, l_fix_arcl_eg_levs,    &
        l_fix_drop_settle, l_emis_ssi_full, l_iau_pc2check,   &
        l_fail_p_layers_inconsis, l_pc2_homog_turb_q_neg,     &
        l_methox_fix, l_stph_rhcrit_unbias,                   &
        l_dtcanfix, l_fix_nh4no3_equilibrium,                 &
        l_eg_damp_height_lid, l_fix_mphys_diags_iter,         &
        l_rm_hardwire_gas360, l_fix_ctile_orog,               &
        l_fix_conv_precip_evap, l_fix_ukca_impscav,           &
        l_fix_rp_shock_amp, l_fix_ustar_dust, l_fix_dyndiag,  &
        l_fix_riming, l_fix_alb_ice_thick, l_fix_zh,          &
        l_fix_ccb_cct, l_fix_nacl_density,                    &
        l_fix_iau_rim_density, l_fix_albsnow_ts,              &
        l_fix_rcf_mlsnow_icefreemax, l_fix_conv_diags_var,    &
        l_fix_lsp_incs_to_spt, l_fix_ec_gen_hgt,              &
        l_fix_improve_drydep, l_fix_wind_snow

! -----------------------------------------------------------
! -----------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SCIENCE_FIXES_MOD'

CONTAINS

SUBROUTINE warn_temp_fixes()

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WARN_TEMP_FIXES'

INTEGER :: ErrorStatus            ! Return code : 0 Normal Exit : >0 Error
CHARACTER(LEN=errormessagelength) :: cmessage
                                  ! Error message if Errorstatus /=0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0


! -----------------------------------------------------------
! -----------------------------------------------------------
! define whether the fix is appropriate to RCF, UM RUN or both
! and warn the user if the fix is not used.


#if defined(RECON)

IF (.NOT. l_roughnesslength_fix) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Recon run excludes ticket #5401 as l_roughnesslength_fix=.FALSE.'//newline//&
  'TStep=1 10m coastal winds may be unrealistic.'//                   newline//&
  'Coastal-tiling runs advised to check impact first.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_ec_gen_hgt) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Recon run excludes ticket #2070 as l_fix_ec_gen_hgt=.FALSE.'//     newline//&
  'Level height generation for ECMWF data will not be accurate.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

#else

IF ( .NOT. l_fix_conserv ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5456 as l_fixed_conserv=.FALSE.'//      newline//&
  'This affects the accuracy of tracers conservation'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_rm_neg_par) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #4982 as l_rm_neg_par=.FALSE.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_mphys_gr_out) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5626 as l_mphys_gr_out=.FALSE.'//       newline//&
  'This will affect any model runs which use prognostic graupel.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_arcl_eg_levs) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5841 as l_fix_arcl_eg_levs=.FALSE.'//   newline//&
  'This will affect any ENDGame model runs which use aerosol climatologies.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_drop_settle) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #4481 as l_fix_drop_settle=.FALSE.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_emis_ssi_full) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #5551 as l_emis_ssi_full=.FALSE.'//      newline//&
  'Your results will be WRONG unless you have set the'//              newline//&
  'emissivity of sea and sea-ice to 1.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_iau_pc2check) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6498 as l_iau_pc2check=.FALSE.'//       newline//&
  'This can lead to inconsistent cloud fraction and condensate'//     newline//&
  'being generated by the IAU.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fail_p_layers_inconsis) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6091 as l_fail_p_layers_inconsis=.FALSE.'//      &
                                                                      newline//&
  'This will allow inconsistent pressures between model levels to'//  newline//&
  'persist beyond the radiation scheme (and cause hard-to-diagnose failures!)'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_pc2_homog_turb_q_neg) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6514 as l_pc2_homog_turb_q_neg=.FALSE.'//        &
                                                                      newline//&
  'your results might contain negative humidies especially if'//      newline//&
  'forced_cu setting is >0.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_methox_fix) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #6357 as l_methox_fix=.FALSE.'//         newline//&
  'This will affect any model runs which use the basic'//             newline//&
  'methane oxidation scheme.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_stph_rhcrit_unbias) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #47 as l_stph_rhcrit_unbias=.FALSE.'//   newline//&
  'If using Random Parameters, the RHCRIT perturbations can cause'//  newline//&
  'cloud cover biases.'

  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_dtcanfix) THEN
  ErrorStatus = -100
  CMessage    = 'Model run excludes ticket #575 as'// &
                ' l_dtcanfix=.FALSE.'              // &
                ' This will affect any model run.'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_nh4no3_equilibrium) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #646 as l_fix_nh4no3_equilibrium=.FALSE.'//       &
                                                                      newline//&
  'This will affect any model runs which include the formation of'//  newline//&
  'ammonium nitrate within the CLASSIC aerosol scheme.'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_eg_damp_height_lid) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #911 as l_eg_damp_height_lid=.FALSE.'//  newline//&
  'The user is required to manually input the reference height' //    newline//&
  '(via the variable damp_height) rather than use the model top.' //  newline//&
  'This will affect any ENDGame model run.'
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_mphys_diags_iter) THEN
  ErrorStatus = -100
  CMessage    =                                                       newline//&
  'Model run excludes ticket #367 as '//                              newline//&
  'l_fix_mphys_diags_iter=.FALSE.'    //                              newline//&
  'This will affect any model run with more than one iteration of' // newline//&
  'the large scale precipitation scheme.'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_rm_hardwire_gas360) THEN
  ErrorStatus = -100 
  cmessage    =                                                       newline//&
  'Model run excludes ticket #430 as l_rm_hardwire_gas360=.FALSE.'//  newline//&
  'This will hardwire trace gas interpolation to use 360day'//        newline//&
  'calendar even if gregorian calendar set'

  CALL ereport(RoutineName, ErrorStatus, CMessage) 
END IF

IF (.NOT. l_fix_ctile_orog) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #1017 as'// newline//&
                ' l_fix_ctile_orog=.FALSE.'         // newline//&
                ' This will affect runs with coastal tiling.'
                          
  CALL ereport(RoutineName, ErrorStatus, CMessage) 
END IF

IF (.NOT. l_fix_conv_precip_evap) THEN
  ErrorStatus = -100 
  CMessage    = 'Model run excludes ticket #1167 as'               // newline//&
                ' l_fix_conv_precip_evap=.FALSE.'                  // newline//&
                ' This will affect all runs with parameterised'    // newline//&
                ' precipitating convection; evaporation rates for' // newline//&
                ' convective precipitation will be underestimated.'
  CALL ereport(RoutineName, ErrorStatus, CMessage) 
END IF

IF (.NOT. l_fix_ukca_impscav) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #1421 as l_fix_ukca_impscav=.FALSE.'//   newline//&
  'This will affect any model runs which include UKCA.'
          
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_rp_shock_amp) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #1638 as l_fix_rp_shock_amp=.FALSE.'//   newline//&
  'This will affect any model using the Random Parameters 2b Scheme'//newline//&
  'by doubling the size of the shock amplitude.'
        
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_ustar_dust) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #1729 as l_fix_ustar_dust=.FALSE.'//     newline//&
  'This will affect any model runs which include interactive dust.'
          
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_dyndiag) THEN
  ErrorStatus = -100 
  cmessage    =                                                       newline//&
  'Model run excludes ticket #2077 as l_fix_dyndiag=.FALSE.'//        newline//&
  'This will mess up shear-dominated PBLs slightly'

  CALL ereport(RoutineName, ErrorStatus, CMessage) 
END IF

IF (.NOT. l_fix_riming ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #2556 as'//                newline//&
  ' l_fix_riming=.FALSE.'//                                           newline//&
  'This could affect any model runs where l_shape_rime is .TRUE.   '//newline//&
  'and may produce small numbers in the microphysics leading to    '//newline//&
  'the UM crashing in the solver.'
          
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_alb_ice_thick ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from JULES ticket 547 as'//            newline//&
  ' l_fix_alb_ice_thick=.FALSE.'//                                    newline//&
  'This will affect any model runs where l_sice_multilayers is     '//newline//&
  '.TRUE. and will result in an incorrect sea ice thickness being  '//newline//&
  'used in the calculation of bare ice albedo.'
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_zh) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #3005 as l_fix_zh=.FALSE.'//             newline//&
  'This will particularly affect the boundary layer depth  '//        newline//&
  'diagnostic (25 or 3025) in model runs which have non-zero '//      newline//&
  'setting for forced_cu but also have a small impact on '//          newline//&
  'model evolution in all runs.'
          

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_ccb_cct ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3076 as'//                newline//&
  ' l_fix_ccb_cct=.FALSE.  '//                                        newline//&
  'This affects any model runs using 5A or 6A convection schemes.  '//newline//&
  'Without this switch, the convection scheme may occasionally '//    newline//&
  'diagnose cloud-top and cloud-base levels inconsistently.  '//      newline//&
  'In the 5A scheme it may spuriously leave the cloud-top level '//   newline//&
  'unset even where it has set the cloud-base, leading to '//         newline//&
  'inconsistent ccrad fields.'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_nacl_density ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #2405 as'//                newline//&
  ' l_fix_nacl_density=.FALSE.'//                                     newline//&
  'This affects the sea-salt emissions generated by the model for  '//newline//&
  'the UKCA + GLOMAP aerosol model'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (l_iau .AND. model_type == mt_lam .AND. .NOT.l_fix_iau_rim_density) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #2710 as l_fix_iau_rim_density=.FALSE.'//newline//&
  'This will cause rim dryrho updates from the LBCs to be ignored'//  newline//&
  'when using an IAU scheme that inserts increments after the'//      newline//&
  'model basis time.'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_albsnow_ts ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3011 as'//                newline//&
  ' l_fix_albsnow_ts=.FALSE.'//                                       newline//&
  'This affects the albedo of snow as calculated in the two-stream '//newline//&
  'scheme in JULES.'

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_rcf_mlsnow_icefreemax ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3011 as'//                newline//&
  ' l_fix_rcf_mlsnow_icefreemax=.FALSE.'//                            newline//&
  'This affects the cap imposed on the mass of snow at ice points '// newline//&
  'reconfigured to ice-free points in the multilayer snow scheme. '

  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_conv_diags_var ) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #3681 as'//                newline//&
  ' l_fix_conv_diags_var=.FALSE.'//                                   newline//&
  'This affects the diagnostics that are passed from the'//           newline//&
  'convection scheme to the PF model. Without this switch'//          newline//&
  'the increments will incorrectly include contributions'//           newline//&
  'from cloud erosion, a temperature rather than potential'//         newline//&
  'temperature will be passed, and the mass flux will exclude'//      newline//&
  'points where it is reducing with height.'
  CALL ereport(RoutineName, ErrorStatus, CMessage)
END IF

IF (.NOT. l_fix_lsp_incs_to_spt) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #2545 as l_fix_lsp_incs_to_spt=.FALSE.'//newline//&
  'This will mean that the microphysics changes due to mixed phase'// newline//&
  'turbulence will not be seen in the SPT scheme.'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_improve_drydep) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes a change from ticket #1250 as'//                newline//&
  ' l_fix_conv_diags_var=.FALSE.'//                                   newline//&
  ' This will mean that dry deposition velocities are set to null'//  newline//&
  ' for HCl, HOCl, HBr, HOBr, H2SO4, MeOH and Sec_Org and that dry'// newline//&
  ' deposition velocities for 9 tiles are inconsistant with 13/17/27 tiles.'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF (.NOT. l_fix_wind_snow) THEN
  ErrorStatus = -100
  cmessage    =                                                       newline//&
  'Model run excludes ticket #4038 as l_fix_wind_snow=.FALSE.. '//    newline//&
  'This will mean that a zero wind speed will incorrectly be used '// newline//&
  'in the calculation of wind-dependent unloading of snow from '//    newline//&
  'vegetation on timesteps when 10m wind diagnostics are not requested.'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

#endif

! -----------------------------------------------------------
! -----------------------------------------------------------


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE warn_temp_fixes


SUBROUTINE print_nlist_temp_fixes()

USE umPrintMgr, ONLY: umPrint, maxLineLen

IMPLICIT NONE
CHARACTER(LEN=maxLineLen) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_TEMP_FIXES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist temp_fixes', src=ModuleName)

WRITE(lineBuffer,'(A,L1)') ' l_rm_neg_par = ',l_rm_neg_par
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_roughnesslength_fix = ',l_roughnesslength_fix
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_mphys_gr_out = ',l_mphys_gr_out
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_conserv = ', l_fix_conserv
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_arcl_eg_levs = ', l_fix_arcl_eg_levs
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_drop_settle = ', l_fix_drop_settle
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_emis_ssi_full = ', l_emis_ssi_full
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_iau_pc2check = ',l_iau_pc2check
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fail_p_layers_inconsis = ',      &
                             l_fail_p_layers_inconsis
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_pc2_homog_turb_q_neg = ',l_pc2_homog_turb_q_neg
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_methox_fix = ',l_methox_fix
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_stph_rhcrit_unbias = ', l_stph_rhcrit_unbias
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_dtcanfix = ',l_dtcanfix
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_nh4no3_equilibrium = ',      &
                             l_fix_nh4no3_equilibrium
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_eg_damp_height_lid = ',l_eg_damp_height_lid
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_mphys_diags_iter = ',l_fix_mphys_diags_iter
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_rm_hardwire_gas360 = ',l_rm_hardwire_gas360
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_ctile_orog = ',l_fix_ctile_orog
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_conv_precip_evap = ',l_fix_conv_precip_evap
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_ukca_impscav = ',l_fix_ukca_impscav
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_rp_shock_amp = ',l_fix_rp_shock_amp
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_ustar_dust = ',l_fix_ustar_dust
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_dyndiag = ',l_fix_dyndiag
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_riming = ',l_fix_riming
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_alb_ice_thick = ',l_fix_alb_ice_thick
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_zh = ',l_fix_zh
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_ccb_cct = ',l_fix_ccb_cct
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_nacl_density = ',l_fix_nacl_density
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_iau_rim_density = ',l_fix_iau_rim_density
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_albsnow_ts = ',l_fix_albsnow_ts
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_rcf_mlsnow_icefreemax = ',                  &
                             l_fix_rcf_mlsnow_icefreemax
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_conv_diags_var = ',l_fix_conv_diags_var
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_lsp_incs_to_spt = ',l_fix_lsp_incs_to_spt
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_ec_gen_hgt = ',l_fix_ec_gen_hgt
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_improve_drydep', l_fix_improve_drydep
CALL umPrint(lineBuffer,src=ModuleName)
WRITE(lineBuffer,'(A,L1)') ' l_fix_wind_snow = ',l_fix_wind_snow
CALL umPrint(lineBuffer,src=ModuleName)
CALL umPrint('- - - - - - end of namelist - - - - - -', src=ModuleName)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_temp_fixes

#if !defined(LFRIC)
SUBROUTINE read_nml_temp_fixes(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_TEMP_FIXES'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_log = 36

TYPE my_namelist
  SEQUENCE
  LOGICAL :: l_roughnesslength_fix
  LOGICAL :: l_rm_neg_par
  LOGICAL :: l_mphys_gr_out
  LOGICAL :: l_fix_conserv
  LOGICAL :: l_fix_arcl_eg_levs
  LOGICAL :: l_fix_drop_settle
  LOGICAL :: l_emis_ssi_full
  LOGICAL :: l_iau_pc2check
  LOGICAL :: l_fail_p_layers_inconsis
  LOGICAL :: l_pc2_homog_turb_q_neg
  LOGICAL :: l_methox_fix
  LOGICAL :: l_stph_rhcrit_unbias
  LOGICAL :: l_dtcanfix
  LOGICAL :: l_fix_nh4no3_equilibrium
  LOGICAL :: l_eg_damp_height_lid
  LOGICAL :: l_fix_mphys_diags_iter
  LOGICAL :: l_rm_hardwire_gas360
  LOGICAL :: l_fix_ctile_orog
  LOGICAL :: l_fix_conv_precip_evap
  LOGICAL :: l_fix_ukca_impscav
  LOGICAL :: l_fix_rp_shock_amp
  LOGICAL :: l_fix_ustar_dust
  LOGICAL :: l_fix_dyndiag
  LOGICAL :: l_fix_riming
  LOGICAL :: l_fix_alb_ice_thick
  LOGICAL :: l_fix_zh
  LOGICAL :: l_fix_ccb_cct
  LOGICAL :: l_fix_nacl_density
  LOGICAL :: l_fix_iau_rim_density
  LOGICAL :: l_fix_albsnow_ts
  LOGICAL :: l_fix_rcf_mlsnow_icefreemax
  LOGICAL :: l_fix_conv_diags_var
  LOGICAL :: l_fix_lsp_incs_to_spt
  LOGICAL :: l_fix_ec_gen_hgt
  LOGICAL :: l_fix_improve_drydep
  LOGICAL :: l_fix_wind_snow
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=temp_fixes, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist temp_fixes", iomessage)

  my_nml % l_roughnesslength_fix      = l_roughnesslength_fix
  my_nml % l_rm_neg_par               = l_rm_neg_par
  my_nml % l_mphys_gr_out             = l_mphys_gr_out
  my_nml % l_fix_conserv              = l_fix_conserv
  my_nml % l_fix_arcl_eg_levs         = l_fix_arcl_eg_levs
  my_nml % l_fix_drop_settle          = l_fix_drop_settle
  my_nml % l_emis_ssi_full            = l_emis_ssi_full
  my_nml % l_iau_pc2check             = l_iau_pc2check
  my_nml % l_fail_p_layers_inconsis   = l_fail_p_layers_inconsis
  my_nml % l_pc2_homog_turb_q_neg     = l_pc2_homog_turb_q_neg
  my_nml % l_methox_fix               = l_methox_fix
  my_nml % l_stph_rhcrit_unbias       = l_stph_rhcrit_unbias
  my_nml % l_dtcanfix                 = l_dtcanfix
  my_nml % l_fix_nh4no3_equilibrium   = l_fix_nh4no3_equilibrium
  my_nml % l_eg_damp_height_lid       = l_eg_damp_height_lid
  my_nml % l_fix_mphys_diags_iter     = l_fix_mphys_diags_iter
  my_nml % l_rm_hardwire_gas360       = l_rm_hardwire_gas360
  my_nml % l_fix_ctile_orog           = l_fix_ctile_orog
  my_nml % l_fix_conv_precip_evap     = l_fix_conv_precip_evap
  my_nml % l_fix_ukca_impscav         = l_fix_ukca_impscav
  my_nml % l_fix_rp_shock_amp         = l_fix_rp_shock_amp
  my_nml % l_fix_ustar_dust           = l_fix_ustar_dust
  my_nml % l_fix_dyndiag              = l_fix_dyndiag
  my_nml % l_fix_riming               = l_fix_riming
  my_nml % l_fix_alb_ice_thick        = l_fix_alb_ice_thick
  my_nml % l_fix_zh                   = l_fix_zh
  my_nml % l_fix_ccb_cct              = l_fix_ccb_cct
  my_nml % l_fix_nacl_density         = l_fix_nacl_density
  my_nml % l_fix_iau_rim_density      = l_fix_iau_rim_density
  my_nml % l_fix_albsnow_ts           = l_fix_albsnow_ts
  my_nml % l_fix_rcf_mlsnow_icefreemax = l_fix_rcf_mlsnow_icefreemax
  my_nml % l_fix_conv_diags_var       = l_fix_conv_diags_var
  my_nml % l_fix_lsp_incs_to_spt      = l_fix_lsp_incs_to_spt
  my_nml % l_fix_ec_gen_hgt           = l_fix_ec_gen_hgt
  my_nml % l_fix_improve_drydep       = l_fix_improve_drydep
  my_nml % l_fix_wind_snow            = l_fix_wind_snow
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  l_roughnesslength_fix      = my_nml % l_roughnesslength_fix
  l_rm_neg_par               = my_nml % l_rm_neg_par
  l_mphys_gr_out             = my_nml % l_mphys_gr_out
  l_fix_conserv              = my_nml % l_fix_conserv
  l_fix_arcl_eg_levs         = my_nml % l_fix_arcl_eg_levs
  l_fix_drop_settle          = my_nml % l_fix_drop_settle
  l_emis_ssi_full            = my_nml % l_emis_ssi_full
  l_iau_pc2check             = my_nml % l_iau_pc2check
  l_fail_p_layers_inconsis   = my_nml % l_fail_p_layers_inconsis
  l_pc2_homog_turb_q_neg     = my_nml % l_pc2_homog_turb_q_neg
  l_methox_fix               = my_nml % l_methox_fix
  l_stph_rhcrit_unbias       = my_nml % l_stph_rhcrit_unbias
  l_dtcanfix                 = my_nml % l_dtcanfix
  l_fix_nh4no3_equilibrium   = my_nml % l_fix_nh4no3_equilibrium
  l_eg_damp_height_lid       = my_nml % l_eg_damp_height_lid
  l_fix_mphys_diags_iter     = my_nml % l_fix_mphys_diags_iter
  l_rm_hardwire_gas360       = my_nml % l_rm_hardwire_gas360
  l_fix_ctile_orog           = my_nml % l_fix_ctile_orog
  l_fix_conv_precip_evap     = my_nml % l_fix_conv_precip_evap
  l_fix_ukca_impscav         = my_nml % l_fix_ukca_impscav
  l_fix_rp_shock_amp         = my_nml % l_fix_rp_shock_amp
  l_fix_ustar_dust           = my_nml % l_fix_ustar_dust
  l_fix_dyndiag              = my_nml % l_fix_dyndiag
  l_fix_riming               = my_nml % l_fix_riming
  l_fix_alb_ice_thick        = my_nml % l_fix_alb_ice_thick
  l_fix_zh                   = my_nml % l_fix_zh
  l_fix_ccb_cct              = my_nml % l_fix_ccb_cct
  l_fix_nacl_density         = my_nml % l_fix_nacl_density
  l_fix_iau_rim_density      = my_nml % l_fix_iau_rim_density
  l_fix_albsnow_ts           = my_nml % l_fix_albsnow_ts
  l_fix_rcf_mlsnow_icefreemax = my_nml % l_fix_rcf_mlsnow_icefreemax
  l_fix_conv_diags_var       = my_nml % l_fix_conv_diags_var    
  l_fix_lsp_incs_to_spt      = my_nml % l_fix_lsp_incs_to_spt   
  l_fix_ec_gen_hgt           = my_nml % l_fix_ec_gen_hgt
  l_fix_improve_drydep       = my_nml % l_fix_improve_drydep
  l_fix_wind_snow            = my_nml % l_fix_wind_snow
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_temp_fixes
#endif

END MODULE science_fixes_mod

