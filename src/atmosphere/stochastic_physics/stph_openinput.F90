! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_openinput_mod

USE missing_data_mod, ONLY: imdi
USE umPrintMgr, ONLY:                                                   &
    umPrint, umMessage, newline

IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_OPENINPUT_MOD'

CONTAINS


SUBROUTINE stph_openinput(coeftype, basefilename)

! Opens the file stream using the name in the
! basefilename input parameter.
! If the stphseed value from run_stochastic namelist has value 0 or 1, a 
! string from the input parameter "coeftype" is appended to the filename
! Pattern coefficients will only be IO if l_stphseed_file is true

USE filenamelength_mod, ONLY: filenamelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype

USE stochastic_physics_run_mod,  ONLY: stphseed_unit, stphseed
USE timestep_mod, ONLY: timestep_number, timestep

USE file_manager, ONLY: assign_file_unit
USE history, ONLY: blank_file_name
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



! Local variables
INTEGER           :: icode       ! Error return strings

CHARACTER(LEN=filenamelength)  :: filename    ! Name of filename
CHARACTER(LEN=filenamelength)  :: basefilename ! BaseName of filename
CHARACTER(LEN=11) :: fileformat  ! Format of filename's contents
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error string
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*)  :: coeftype    ! String identifying coefs src

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'STPH_OPENINPUT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Default file format is unformatted
fileformat="unformatted"

! Checks If output is required, and if unit is already open
IF (mype == 0) THEN

  IF (TRIM(basefilename) == "") THEN
    cmessage = 'Random seed filename **rpseed** not set in namelist'
    icode = 10
    CALL ereport(RoutineName, icode, cmessage)
  END IF

  ! Append scheme name suffix to filename (e.g. SKEB, RP2, SPT)
  ! Note: If this is a continuation (CRUN) this file should 
  !       have been written during the previous timestep
  IF ((stphseed == 0 .OR. stphseed == 1) .AND. .NOT. l_nrun_as_crun) THEN
    WRITE(filename,'(A,I10.10)')                                        &
          TRIM(basefilename) // "_" // TRIM(coeftype) // "_",           &
          timestep_number - 1
  END IF
  
! For safety close the unit which connects to a previously modified file 
! before opening read only.
  IF ( stphseed_unit /= imdi) CLOSE(UNIT=stphseed_unit)
  CALL assign_file_unit(filename, stphseed_unit, handler="fortran")
  OPEN(UNIT=stphseed_unit,FILE=TRIM(filename),STATUS='old',             &
      ACTION='read', FORM=fileformat,IOSTAT=icode, IOMSG=iomessage)

  IF (icode == 0) THEN  
    WRITE(umMessage,'(A)') "Read StochPhys Fields: " // TRIM(filename)  &
          // " has been opened as " // fileformat
    CALL umPrint(umMessage,src='stph_openinput')
  ELSE
    WRITE(umMessage,'(A,I0,A)')                                         &
          "RPSEED IN: Error opening file " // TRIM(filename)            &
          // " Icode=", icode, newline // 'IoMsg: ' // TRIM(iomessage)
    CALL umPrint(umMessage,src='stph_openinput')
    cmessage = ' Error opening file:'//TRIM(filename) // ' :' //        &
                 newline // 'IoMsg: ' // TRIM(iomessage)
    CALL ereport(RoutineName, icode, cmessage)
  END IF

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stph_openinput
END MODULE stph_openinput_mod
