! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : LOGIC

MODULE s_logic

USE s_maxdim,  ONLY:                                                        &
  mx_rw_lng, mx_rw

USE scm_cntl_mod,  ONLY: scm_nml

USE scm_utils, ONLY:                                                        &
  imdi, rw_lng, rw

USE s_main_force, ONLY:                                                     &
  dev_test, test, altdat, geoforce, geoinit, noforce, obs, obs_surf, stats  &
, local_time, ancyc, l_flux_bc, l_spec_z0, l_geo_centred, l_qpos_for        &
, grafdump_day, grafdump_days, grafdump_step, prindump_day                  &
, prindump_days, prindump_obs, prindump_step, prinstat                      &
, radcloud_fixed, tapein, tapeout                                           &
, scm_land_sea_mask => land_sea_mask &
, scm_land_ice_mask => land_ice_mask &
, scm_soil_mask     => soil_mask     &
, scm_cumulus       => cumulus

USE bl_option_mod, ONLY: flux_bc_opt, interactive_fluxes

USE umPrintMgr, ONLY:      &
    umPrint, umMessage, newline

USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


!=============================================================================
!
! Description:
!   Allows SCM to read in LOGIC namelist from forcing file, <scm_nml>
!
! Method:
!   Namelist LOGIC is defined in this module and read in by contained
!   subroutine read_logic.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped before being transferred to arrays of the correct size/shape in
!   s_main_force.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

LOGICAL, PRIVATE ::                      &
  land_sea_mask (mx_rw_lng*mx_rw) = .FALSE. &! land point flag (any type)
, land_ice_mask (mx_rw_lng*mx_rw) = .FALSE. &! land point flag (ice)
, soil_mask     (mx_rw_lng*mx_rw) = .FALSE. &! land point flag (soil)
, cumulus       (mx_rw_lng*mx_rw) = .FALSE.  ! bl convection flag


NAMELIST/logic/                                                             &
  dev_test, test, altdat, geoforce, geoinit, noforce, obs, obs_surf, stats  &
, local_time, ancyc, l_flux_bc, l_spec_z0, l_geo_centred, l_qpos_for        &
, grafdump_day, grafdump_days, grafdump_step, prindump_obs, prinstat        &
, prindump_day, prindump_days, prindump_step, radcloud_fixed, tapein        &
, tapeout, land_sea_mask, land_ice_mask, soil_mask, cumulus

PRIVATE :: logic

!=============================================================================
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='S_LOGIC'

CONTAINS

SUBROUTINE read_logic

IMPLICIT NONE


! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

INTEGER :: istatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_LOGIC'

!=========================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), ACTION='READ', IOSTAT=istatus, &
     IOMSG=iomessage)

IF (istatus /= 0) THEN
  icode = 500
  WRITE(umMessage,*) ' Error opening file on unit 10 from '//routinename
  CALL umPrint(umMessage,src='s_logic')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
  CALL umPrint(umMessage,src='s_logic')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='s_logic')
  WRITE(umMessage,*) ' Error: ', TRIM(iomessage)
  CALL umPrint(umMessage,src='s_logic')
END IF

READ  (10, logic)
CLOSE (10)

! Check l_flux_bc setting is consistent with UM switch
IF (l_flux_bc) THEN 
  IF (flux_bc_opt == interactive_fluxes) THEN
    icode = 50
    WRITE(umMessage,'(A,I1)')                                         newline//&
       ' Error: l_flux_bc is true in your SCM namelist but this '//   newline//&
       '        now requires flux_bc_opt to be set in the '//         newline//&
       '        run_bl namelist of the UM'
    CALL ereport(RoutineName,icode,umMessage)
  END IF
ELSE
  IF (flux_bc_opt /= interactive_fluxes) THEN
    WRITE(umMessage,'(A,I1)')                                         newline//&
       ' Warning: l_flux_bc is false in your SCM namelist but '//     newline//&
       '        flux_bc_opt is not zero so the SCM will be '//        newline//&
       '        using the SCM driving data - this may be correct!'
    CALL umPrint(umMessage,src='s_logic')
  END IF
END IF

scm_land_sea_mask = RESHAPE( land_sea_mask, (/rw_lng,rw/) )

scm_land_ice_mask = RESHAPE( land_ice_mask, (/rw_lng,rw/) )

scm_soil_mask     = RESHAPE( soil_mask,     (/rw_lng,rw/) )

scm_cumulus       = RESHAPE( cumulus,       (/rw_lng,rw/) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_logic

!=============================================================================
END MODULE s_logic
