! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: READHIST
!
!  Purpose: Initialise full History either from history_mod or
!           from most recent record in history file input,
!           OR initialise filenames only for a new run.
!
!  Code description:
!    Language: FORTRAN 90
!    Programming standard: UMDP3
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

SUBROUTINE readhist( icode,cmessage )

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE filenamelength_mod,    ONLY: filenamelength
USE nlcfiles_namelist_mod, ONLY: read_nml_nlcfiles, print_nlist_nlcfiles
USE nlstcall_pp_namelist_mod, ONLY: read_nml_nlstcall_pp_hist
USE nlstcall_nc_namelist_mod, ONLY: read_nml_nlstcall_nc_hist
USE history,               ONLY: blank_file_name, read_nml_nlihisto,       &
                                 read_nml_nlihistg, read_nml_nlchisto,     &
                                 read_nml_nlchistg, history_file_exists,   &
                                 history_unit
USE ereport_mod,           ONLY: ereport
USE umPrintMgr,            ONLY: umPrint, umMessage, PrStatus_min,         &
                                 printstatus, prstatus_normal, newline
USE um_parcore,            ONLY: mype
USE mpl,                   ONLY: mpl_logical
USE file_manager,          ONLY: get_file_unit_by_id, assign_file_unit,    &
                                 release_file_unit, get_file_by_id, um_file_type
USE get_env_var_mod,       ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
INTEGER                           :: icode     ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage

! Local variables
INTEGER                       :: i    ! Loop counter
INTEGER                       :: UNIT ! File unit

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'READHIST'
CHARACTER(LEN=filenamelength) :: filename
CHARACTER(LEN=8)              :: nlname
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER :: my_comm

TYPE(um_file_type), POINTER :: shared_file ! Object to contain SHARED file data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

icode = 0
!
! 2. Open history file and rewind
!
CALL get_env_var("HISTORY",filename)

IF (mype == 0) THEN
  INQUIRE(FILE=filename, EXIST=history_file_exists)
  IF (history_file_exists) THEN
    ! History available - CRun
    WRITE(umMessage,'(A)') 'Opening history file: '//TRIM(filename)
    CALL umPrint(umMessage,level=PrStatus_min,src=routinename)

    CALL assign_file_unit(filename, history_unit, handler="fortran")
    OPEN(history_unit,FILE=filename,FORM='FORMATTED',ACTION='READ',      &
         IOSTAT=icode, IOMSG=iomessage)

    IF (icode <= 0) THEN
      IF (icode < 0) THEN
        WRITE(umMessage,'(A)')                                             &
          'READHIST: Warning message on OPEN of permanent history file:' //&
                                                                  newline//&
          'IoMsg: '//TRIM(iomessage)

        CALL umPrint(umMessage,level=PrStatus_min,src=routinename)
        WRITE(umMessage,'(A,I6)')'IOSTAT= ',icode
        CALL umPrint(umMessage,level=PrStatus_min,src=routinename)
      END IF
    ELSE ! Open failed
      cmessage='READHIST: Failed in OPEN of history file'
    END IF

  END IF !!! history_file_exists

END IF  !!! mype==0

CALL MPl_bcast(history_file_exists,1,mpl_logical,0,my_comm,icode)

IF (history_file_exists) THEN

  ! 3. Read most recent records
  CALL read_nml_nlihisto(history_unit)

  CALL read_nml_nlchisto(history_unit)

  CALL read_nml_nlihistg(history_unit)

  CALL read_nml_nlchistg(history_unit)

  CALL read_nml_nlcfiles(history_unit)
  IF (printstatus >= prstatus_normal .AND. mype == 0) THEN
    CALL print_nlist_nlcfiles()
  END IF

  CALL read_nml_nlstcall_pp_hist(history_unit)
  IF ( mype == 0 ) REWIND(history_unit)

  CALL read_nml_nlstcall_nc_hist(history_unit)

  ! Close the history unit since we no longer
  ! need it
  IF (mype == 0) THEN
     CLOSE(history_unit)
     CALL release_file_unit(history_unit, handler="fortran")
  END IF

ELSE

  ! No history file available
  ! Initialise filenames from history namelist in SHARED file.
  ! Retreive the SHARED file data:
  shared_file => get_file_by_id("shared", handler="fortran")
  WRITE(umMessage,'(A)')                                                      &
    'History file: '//TRIM(filename)//' not available.'//            newline//&
    'Input/output filenames initialised from file: '//                        &
    TRIM(shared_file % filename)//                                   newline//&
    'All other history namelists initialised to default values.'
  CALL umPrint(umMessage,level=PrStatus_min,src=routinename)

  ! Do not need to open SHARED file (opened in UM shell)
  CALL read_nml_nlcfiles(shared_file % unit)
  IF (printstatus >= prstatus_normal .AND. mype == 0) THEN
    CALL print_nlist_nlcfiles()
  END IF
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE readhist
