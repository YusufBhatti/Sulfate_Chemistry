! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : PHYSWITCH

MODULE s_physwitch

USE scm_cntl_mod, ONLY: scm_nml

USE s_main_force, ONLY: conv_mode

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
!   Allows SCM to read in PHYSWITCH namelist from forcing file, scm_nml.
!   PHYSWITCH contains information on the different physics modes.
!
! Method:
!   Namelist PHYSWITCH is defined in this module and read in by contained
!   subroutine read_physwitch.  Scalar variables are read directly to
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

NAMELIST/physwitch/ conv_mode

PRIVATE :: physwitch

!=============================================================================
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='S_PHYSWITCH'

CONTAINS

SUBROUTINE read_physwitch

IMPLICIT NONE


! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

INTEGER :: istatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_PHYSWITCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), ACTION='READ', IOSTAT=istatus, &
     IOMSG=iomessage)

IF (istatus /= 0) THEN
  icode = 500
  WRITE(umMessage,*) ' Error opening file on unit 10 from '//routinename
  CALL umPrint(umMessage,src='s_physwitch')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
  CALL umPrint(umMessage,src='s_physwitch')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='s_physwitch')
  WRITE(umMessage,*) ' Error: ', TRIM(iomessage)
  CALL umPrint(umMessage,src='s_physwitch')
END IF

READ  (10, physwitch)
CLOSE (10)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_physwitch

!=============================================================================
END MODULE s_physwitch
