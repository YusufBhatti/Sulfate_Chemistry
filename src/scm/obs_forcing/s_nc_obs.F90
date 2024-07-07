! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : NC_OBS

MODULE s_nc_obs

USE scm_cntl_mod, ONLY: scm_nml

USE s_main_force, ONLY: netcdf_file, obs_t0_prf, source

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


!=============================================================================
!
! Description:
!   Allows SCM to read in NC_OBS namelist from forcing file, scm_nml.
!   NC_OBS contains information relating to the use of NetCDF driver files to
!   force the SCM
!
! Method:
!   Namelist NC_OBS is defined in this module and read in by contained
!   subroutine read_nc_obs.  Scalar variables are read directly to
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

  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
NAMELIST/nc_obs/                                                            &
 netcdf_file, obs_t0_prf, source

PRIVATE :: nc_obs

!=============================================================================
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='S_NC_OBS'

CONTAINS

SUBROUTINE read_nc_obs

IMPLICIT NONE


! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

INTEGER :: istatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NC_OBS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set defaults
netcdf_file = ''
obs_t0_prf  = 1

OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), ACTION='READ', IOSTAT=istatus, &
     IOMSG=iomessage)

IF (istatus /= 0) THEN
  icode = 500
  WRITE(umMessage,*) ' Error opening file on unit 10 from '//routinename
  CALL umPrint(umMessage,src='s_nc_obs')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
  CALL umPrint(umMessage,src='s_nc_obs')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='s_nc_obs')
  WRITE(umMessage,*) ' Error: ', iomessage
  CALL umPrint(umMessage,src='s_nc_obs')
END IF

READ  (10, nc_obs)
CLOSE (10)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_nc_obs

!=============================================================================
END MODULE s_nc_obs
