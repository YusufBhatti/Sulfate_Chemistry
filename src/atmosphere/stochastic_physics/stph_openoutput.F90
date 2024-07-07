! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_openoutput_mod

USE umPrintMgr, ONLY:                                                   &
    umPrint, umMessage, newline
IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_OPENOUTPUT_MOD'

CONTAINS


SUBROUTINE stph_openoutput(io_pos, coeftype)

! Opens the file stream using the name in the
! basefilename input parameter.
! If the stphseed value from run_stochastic namelist has value 0 or 1, a 
! string from the input parameter "coeftype" is appended to the filename
! Pattern coefficients will only be IO if l_stphseed_file is true


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength

USE UM_ParVars
USE UM_ParCore, ONLY: mype

USE nlcfiles_namelist_mod, ONLY: stphseed_file => rp2_seed
USE file_manager, ONLY: assign_file_unit
USE stochastic_physics_run_mod,  ONLY: stphseed_unit
USE timestep_mod, ONLY: timestep_number, timestep
USE history, ONLY: blank_file_name

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



! Local variables

INTEGER           :: icode       ! Error return strings
LOGICAL           :: fileexists  ! Check If file already exists

CHARACTER(LEN=filenamelength)  :: filename     ! Name of filename
CHARACTER(LEN=filenamelength)  :: basefilename ! BaseName of filename
CHARACTER(LEN=10) :: io_pos      ! IO position (REWIND, APPEND, ASIS)
CHARACTER(LEN=15) :: fileformat  ! Name of filename
CHARACTER(LEN=errormessagelength) :: cmessage    ! error string
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*)  :: coeftype    ! String identifying coefs src

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'STPH_OPENOUTPUT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Default file format is formatted
fileformat="unformatted"

! Checks If output file is active, and if unit is already open
! Output mode will be used in CRUNs after initial read
IF (mype == 0) THEN
  basefilename =  stphseed_file

  IF (TRIM(basefilename) == "") THEN
    cmessage = 'Random Field filename **rpseed** not set in namelist'
    icode = 10
    CALL ereport(RoutineName, icode, cmessage)
  END IF

  ! Append scheme name suffix to filename (e.g. SKEB, SPT)
  WRITE(filename,'(A,I10.10)')                                          &
        TRIM(basefilename) // "_" // TRIM(coeftype) // "_",             &
        timestep_number

  !----------------------------------------------------------------
  ! If called with APPEND position, we will be appending the
  ! dpsidtc/s coeffs to the seed file, to be read in at the next CRUN
  !----------------------------------------------------------------
  CALL assign_file_unit(filename, stphseed_unit, handler="fortran")
  OPEN(UNIT=stphseed_unit, FILE=TRIM(filename), ACTION='write',         &
       POSITION=io_pos, FORM=fileformat, IOSTAT=icode, IOMSG=iomessage)
  IF (icode > 0) THEN
    WRITE(umMessage,'(A,I0,A)')                                         &
      "Write Random Field: Failed to open:"//                  newline//&
      TRIM(filename)//                                         newline//&
      "icode: ", icode,                                        newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL umPrint(umMessage,src='stph_openoutput')
    cmessage = 'Cannot open Random Field file for writing'
    CALL ereport(RoutineName, icode, cmessage)
  END IF
  IF (icode == 0) THEN
    WRITE(umMessage,'(A)') "Write Random Fields: " // TRIM(filename)    &
        // " opened for " // coeftype
    CALL umPrint(umMessage,src='stph_openoutput')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE stph_openoutput

END MODULE stph_openoutput_mod
