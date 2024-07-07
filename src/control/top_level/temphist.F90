! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine : TEMPHIST -------------------------------------------------
!
! Purpose :Write current contents of history module to temporary
!          or interim history file - overwriting previous record
!
! Documentation:  Unified Model Documentation Paper
!                  H- History Bricks

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

MODULE temphist_mod

IMPLICIT NONE 

! Identifiers for which type of history file to write
INTEGER, PARAMETER :: write_temporary = 1
INTEGER, PARAMETER :: write_interim = 2

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TEMPHIST_MOD'

CONTAINS

SUBROUTINE temphist                                                           &
    ( write_file,icode,cmessage )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstcall_pp_namelist_mod, ONLY: write_nml_nlstcall_pp_hist
USE nlstcall_nc_namelist_mod, ONLY: write_nml_nlstcall_nc_hist
USE filenamelength_mod,     ONLY:                                             &
    filenamelength
USE MPPIO_file_utils,       ONLY:                                             &
    file_copy,                                                                &
    dir_create
USE umPrintMgr, ONLY: printstatus, prstatus_oper, ummessage, umprint, newline
USE ereport_mod,            ONLY: ereport
USE nlcfiles_namelist_mod, ONLY: nlcfiles
USE file_manager, ONLY: assign_file_unit, release_file_unit, um_file_type,    &
                        init_file_loop
USE missing_data_mod, ONLY: imdi
USE history, ONLY:                                                            &
    nlihistg, nlihisto, nlchistg, nlchisto
USE get_env_var_mod, ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)           :: write_file ! Flag indicating file 
INTEGER, INTENT(INOUT)        :: icode      ! Return code from routine

CHARACTER(LEN=errormessagelength) :: cmessage 


INTEGER, SAVE                 :: sequence_id=1
INTEGER, PARAMETER            :: sequence_id_max=9999
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
INTEGER                       :: i_unit
INTEGER                       :: i_id
INTEGER                       :: history_unit
LOGICAL, SAVE                 :: first=.TRUE.
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEMPHIST'
CHARACTER(LEN=filenamelength) :: filename
CHARACTER(LEN=filenamelength) :: filename_base
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER   :: hist_dir="history_archive"
TYPE(um_file_type), POINTER   :: um_file

!
! 1. Open, rewind and write a new record, assign a versioned name for the
!    file
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (First) THEN
  First=.FALSE.
  CALL dir_create(hist_dir,force_local=.TRUE.)
END IF

! Get the filename of the history file to write to, depending on the flag
! passed into this routine
IF (write_file == write_temporary) THEN
  CALL get_env_var("HISTORY_TEMP",filename_base)
ELSE IF (write_file == write_interim) THEN
  CALL get_env_var("HISTORY",filename_base)
ELSE
  CALL ereport('temphist',icode,                                             &
      'Called with an invalid file indicator')
END IF

WRITE(filename,'(A,A,A,A,I4.4)')                                             &
    hist_dir,                                                                &
    '/',                                                                     &
    'temp_hist',                                                             &
    '.',                                                                     &
    sequence_id

IF (printstatus>=prstatus_oper) THEN
  WRITE(umMessage,'(A,A)')'Temphist: Writing temp history as :',             &
      TRIM(filename)
  CALL umPrint(umMessage,src='temphist')
  WRITE(umMessage,'(A,A)')'Temphist: Final destination       :',             &
      TRIM(filename_base)
  CALL umPrint(umMessage,src='temphist')
END IF

CALL assign_file_unit(filename, history_unit, handler="fortran")

OPEN(history_unit,FILE=filename,FORM='FORMATTED',                            &
    DELIM='APOSTROPHE',IOSTAT=icode, IOMSG=iomessage)

IF (icode <= 0) THEN
  IF (icode < 0) THEN
    WRITE(umMessage,'(A)')                                                   &
      'Temphist: Warning message on OPEN of history file:'//        newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL umPrint(umMessage,src='temphist')
    WRITE(umMessage,*)'iostat= ',icode
    CALL umPrint(umMessage,src='temphist')
  END IF

  REWIND(history_unit)

  WRITE(history_unit,nlihisto,IOSTAT=icode, IOMSG=iomessage)
  CALL ereport('temphist',icode,                                             &
      'Write ERROR on history file for namelist NLIHISTO:'//        newline//&
      'IoMsg: '//TRIM(iomessage))

  WRITE(history_unit,nlchisto,IOSTAT=icode, IOMSG=iomessage)
  CALL ereport('temphist',icode,                                             &
      'Write ERROR on history file for namelist NLCHISTO:'//        newline//&
      'IoMsg: '//TRIM(iomessage))

  WRITE(history_unit,nlihistg,IOSTAT=icode, IOMSG=iomessage)
  CALL ereport('temphist',icode,                                             &
      'Write ERROR on history file for namelist NLIHISTG:'//        newline//&
      'IoMsg: '//TRIM(iomessage))

  WRITE(history_unit,nlchistg,IOSTAT=icode, IOMSG=iomessage)
  CALL ereport('temphist',icode,                                             &
      'Write ERROR on history file for namelist NLCHISTG:'//        newline//&
      'IoMsg: '//TRIM(iomessage))

  WRITE(history_unit,nlcfiles,IOSTAT=icode, IOMSG=iomessage)
  CALL ereport('temphist',icode,                                             &
      'Write ERROR on history file for namelist NLCFILES:'//        newline//&
      'IoMsg: '//TRIM(iomessage))

  ! Write out the nlstcall_pp namelists
  CALL write_nml_nlstcall_pp_hist(history_unit)

  ! Write out the nlstcall_nc namelists
  CALL write_nml_nlstcall_nc_hist(history_unit)
  !
  ! 2. Close and return
  !
  CLOSE(history_unit)
  CALL release_file_unit(history_unit, handler="fortran")

  ! Instigate copy from versioned name to final name.
  CALL file_copy(TRIM(filename),TRIM(filename_base))

  sequence_id=sequence_id+1
  IF (sequence_id>sequence_id_max) THEN
    sequence_id=1
  END IF

ELSE ! Open failed
  cmessage='Temphist: Failed in OPEN of history file'
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE temphist

END MODULE temphist_mod
