! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE chk_opts_mod
! Description:
! This module is the interface to a utility subroutine
! that checks namelist input variables for sensible
! values based on developer input strings or an integer array
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr,             ONLY: umprint, ummessage, newline
USE ereport_mod,            ONLY: ereport

IMPLICIT NONE

PRIVATE
PUBLIC :: chk_var, def_src

INTEGER, PARAMETER :: max_n_conditions = 15
INTEGER, PARAMETER :: condition_len    = 20
INTEGER, PARAMETER :: str_len          = 100

CHARACTER(LEN=*), PARAMETER :: modname ='CHK_OPTS_MOD'

CHARACTER(LEN=2) :: operators(max_n_conditions)
CHARACTER(LEN=2) :: operator_string

CHARACTER(LEN=errormessagelength) :: errstr
CHARACTER(LEN=condition_len) :: values  (max_n_conditions)
CHARACTER(LEN=condition_len) :: lbounds (max_n_conditions)
CHARACTER(LEN=condition_len) :: ubounds (max_n_conditions)
CHARACTER(LEN=condition_len) :: tmp_str

CHARACTER(LEN=str_len) :: srcname
CHARACTER(LEN=str_len) :: loc_str
CHARACTER(LEN=str_len) :: loc_name
CHARACTER(LEN=str_len) :: loc_src

CHARACTER(LEN=30) :: fmt_str
CHARACTER(LEN=8)  :: result_chk

LOGICAL :: pass_chk
INTEGER :: nconditions

! Setup inteface to allow user to call same routinename in order to check
! supported datatypes. In order to pass, the variable must pass ANY of
! the code owner specified checks, these checks should be consistent with
! the UM rose-metadata
INTERFACE chk_var
  MODULE PROCEDURE chk_int_str   ! Integers (values/ranges)
  MODULE PROCEDURE chk_real_str  ! Reals    (values/ranges)
  MODULE PROCEDURE chk_int_arr   ! Integers (values only)
END INTERFACE

CHARACTER(LEN=str_len) :: def_src = '' ! Calling routine if specified

CONTAINS

SUBROUTINE chk_int_str (var, name, string, report_pass, cmessage)

! Checks an integer variable against a list of criteria
! passed as a string in the form '[Check_A, Check_B, ... ]'.
!
! Recognised checks
!  * ==Integer (This is assumed if no operator is provided)
!  * <Integer
!  * >Integer
!  * >=Integer
!  * <=Integer
!  * Integer lower bound : Integer upper bound
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if ANY of the criteria
! are met. If none of the specified checks are met, ereport will be
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If Chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

IMPLICIT NONE

INTEGER,          INTENT(IN) :: var       ! The integer variable to be checked
CHARACTER(LEN=*), INTENT(IN) :: name      ! The name of the test variable
CHARACTER(LEN=*), INTENT(IN) :: string    ! String containg test criteria

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
LOGICAL,          INTENT(INOUT), OPTIONAL :: report_pass

! Additonal comments to af to ereport in event of failure of tests
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmessage


! Local variables
CHARACTER(LEN=*), PARAMETER :: routinename ='CHK_INT_STR'

INTEGER :: i, index_val, icode
INTEGER :: condition
INTEGER :: lowerbound
INTEGER :: upperbound

!==========================================================
IF (.NOT. PRESENT(cmessage)) THEN
  errstr = ''
ELSE
  errstr = newline//TRIM(ADJUSTL(cmessage))
END IF

IF (def_src == '') THEN
  loc_src = ''
  srcname = modname//':'//routinename
ELSE
  loc_src = 'Problem performing defensive checks ('//                          &
            TRIM(ADJUSTL(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//TRIM(ADJUSTL(def_src))
END IF

loc_str  = string
pass_chk   = .FALSE.
icode      = 0

result_chk = 'FAIL'

! Strip off brackets
index_val = INDEX(loc_str, '[')
loc_str(index_val:index_val) = ''
index_val = INDEX(loc_str, ']')
loc_str(index_val:index_val) = ''

loc_str  = TRIM(ADJUSTL(loc_str))
loc_name = TRIM(ADJUSTL(name))

! Check for reals in the condition string
IF ( INDEX(loc_str,'.') /= 0) THEN
  WRITE(errstr,'(A,I0,A)')                                            newline//&
    TRIM(loc_src)//'Attempting to test an integer variable, '//TRIM(loc_name)//&
    ', (',var,')'//newline//'against a condition involving a non-integer, (' //&
    TRIM(loc_str) //&
    ')'
  icode = 10
  CALL ereport(TRIM(srcname), icode, errstr)
END IF

! Separate the condition string into separate tests
CALL split_conditions()

! Check that there are some valid tests
IF (nconditions == 0) THEN
  WRITE(errstr,'(A,I0,A)')                                            newline//&
    TRIM(loc_src)//'No valid checks, '//string//                               &
    ', provided for integer variable, '//TRIM(loc_name)//', (',var,')'
  icode = 30
  CALL ereport(TRIM(srcname), icode, errstr)
END IF


! Loop over the conditions, converting them to integers so they can be tested
DO i=1, nconditions

  operator_string  = operators(i)

  IF (TRIM(operator_string) == ':') THEN

    ! For range tests which have upper and lower bounds
    READ(UNIT=lbounds(i),FMT=*, IOSTAT=icode) lowerbound
    IF (icode /= 0) THEN
      icode = 40
      WRITE(errstr,'(A,I0,A)')                                        newline//&
        TRIM(loc_src)//'Cannot convert lower bound string "'//                 &
        TRIM(lbounds(i))//'" into an integer when'//                  newline//&
        'running defensive checks on '//TRIM(loc_name)//', (',var,').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

    READ(UNIT=ubounds(i),FMT=*, IOSTAT=icode) upperbound
    IF (icode /= 0) THEN
      icode = 50
      WRITE(errstr,'(A,I0,A)')                                        newline//&
        TRIM(loc_src)//'Cannot convert upper bound string "'//                 &
        TRIM(ubounds(i))//'" into an integer when'//                  newline//&
        'running defensive checks on '//TRIM(loc_name)//', (',var,').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

  ELSE

    READ(UNIT=values(i),FMT=*, IOSTAT=icode) condition
    IF (icode /= 0) THEN
      icode = 60
      WRITE(errstr,'(A,I0,A)')                                        newline//&
        TRIM(loc_src)//'Cannot convert string "'//                             &
        TRIM(values(i))//'" into an integer when running '//          newline//&
        'defensive checks on '//TRIM(loc_name)//', (',var,').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

  END IF


  ! Now do the test
  SELECT CASE (TRIM(operator_string))
  CASE ('<')
    IF ( var < condition )  pass_chk = .TRUE.
  CASE ('>')
    IF ( var > condition )  pass_chk = .TRUE.
  CASE ('>=')
    IF ( var >= condition ) pass_chk = .TRUE.
  CASE ('<=')
    IF ( var <= condition ) pass_chk = .TRUE.
  CASE (':')
    IF ( var >= lowerbound .AND. var <= upperbound ) pass_chk = .TRUE.
  CASE ('==')
    IF ( var == condition ) pass_chk = .TRUE.

  CASE DEFAULT
    WRITE(errstr,'(A)')                                               newline//&
      TRIM(loc_src)//'Unrecognised operator, ('//TRIM(operator_string)//       &
      '), when checking'//newline//'integer variable, '//TRIM(loc_name)//      &
      ', (',var,')'
    icode = 70
    CALL ereport(TRIM(srcname), icode, errstr)
    EXIT
  END SELECT

  IF (pass_chk) EXIT

END DO

WRITE(fmt_str,'(A,I0,A)') '(A,I0,A,', nconditions-1, '(I0,","),I0,A)'

IF (PRESENT(report_pass)) THEN
  report_pass = pass_chk
  IF (pass_chk) result_chk = 'pass'
  WRITE(ummessage,'(A,I0,A)')                                                  &
    'Defensive checks on '//TRIM(loc_name)//' (',var,') against '//            &
    TRIM(ADJUSTL(string))//' .... '//TRIM(result_chk)
  CALL umprint(ummessage, src=TRIM(srcname))

ELSE IF (.NOT. pass_chk) THEN

  WRITE(errstr,'(A,I0,A)')                                            newline//&
    TRIM(loc_name)//', (',var,'), has failed to pass any defensive checks, ' //&
    newline//TRIM(ADJUSTL(string))//'.'//TRIM(errstr)
  icode = 80
  CALL ereport(TRIM(srcname), icode, errstr)

END IF

RETURN
END SUBROUTINE chk_int_str



!==============================================================================
!==============================================================================



SUBROUTINE chk_real_str (var, name, string, report_pass, cmessage)

! Checks a real variable against a list of criteria
! passed as a string in the form '[Check_A, Check_B, ... ]'.
!
! Recognised checks
!  * ==Real (This is assumed if no operator is provided)
!  * <Real
!  * >Real
!  * >=Real
!  * <=Real
!  * Real lower bound : Real upper bound
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if ANY of the criteria
! are met. If none of the specified checks are met, ereport will
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If Chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

IMPLICIT NONE

REAL,             INTENT(IN) :: var    ! The real variable to be checked
CHARACTER(LEN=*), INTENT(IN) :: name   ! The name of the test variable
CHARACTER(LEN=*), INTENT(IN) :: string ! String containg test criteria

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
LOGICAL,          INTENT(INOUT), OPTIONAL :: report_pass

! Additonal comments to af to ereport in event of failure of tests
CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: cmessage

! Local variables

CHARACTER(LEN=*), PARAMETER :: routinename ='CHK_REAL_STR'

INTEGER :: i, index_val, icode
REAL :: condition
REAL :: lowerbound
REAL :: upperbound

!==========================================================
IF (.NOT. PRESENT(cmessage)) THEN
  errstr = ''
ELSE
  errstr = newline//TRIM(ADJUSTL(cmessage))
END IF

IF (def_src =='') THEN
  loc_src = ''
  srcname = modname//':'//routinename
ELSE
  loc_src = 'Problem performing defensive checks ('//                          &
            TRIM(ADJUSTL(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//TRIM(ADJUSTL(def_src))
END IF

loc_str    = string
pass_chk   = .FALSE.
icode      = 0
result_chk = 'FAIL'

! Strip off brackets
index_val = INDEX(loc_str, '[')
loc_str(index_val:index_val) = ''
index_val = INDEX(loc_str, ']')
loc_str(index_val:index_val) = ''

WRITE(tmp_str,'(E10.3)') var

loc_str  = ADJUSTL(loc_str)
loc_name = ADJUSTL(name)
tmp_str  = ADJUSTL(tmp_str)

CALL split_conditions ()

IF (nconditions == 0) THEN
  WRITE(errstr,'(A,I0,A)')                                            newline//&
    TRIM(loc_src)//'No valid checks, '//TRIM(ADJUSTL(string))//                &
    ', provided for real variable, '//TRIM(loc_name)//', (',var,')'
  icode = 20
  CALL ereport(TRIM(srcname), icode, errstr)
END IF


DO i=1, nconditions

  operator_string = operators(i)

  IF (TRIM(ADJUSTL(operator_string)) == ':') THEN
    IF ( INDEX(lbounds(i),'.') == 0) THEN

      WRITE(errstr, '(A)')                                            newline//&
        TRIM(loc_src)//'You are attempting to test a real variable, '//        &
        TRIM(loc_name)//'('//TRIM(tmp_str)//'), '//                   newline//&
        'with a non-real lower bound, ('//TRIM(lbounds(i))//')'

      icode = 30
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

    READ(UNIT=lbounds(i),FMT=*, IOSTAT=icode) lowerbound
    IF (icode /= 0) THEN
      icode = 40
      WRITE(errstr,'(A)')                                             newline//&
        TRIM(loc_src)//'Cannot convert lower bound string "'//                 &
        TRIM(lbounds(i))//'" into a real when'//                      newline//&
        'running defensive checks on '//TRIM(loc_name)//', ('//                &
        TRIM(tmp_str)//').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF


    IF ( INDEX(ubounds(i),'.') == 0) THEN
      WRITE(errstr,'(A)')                                             newline//&
        TRIM(loc_src)//'You are attempting to test a real variable, '//        &
        TRIM(loc_name)//'('//TRIM(tmp_str)//'), '//                   newline//&
        'with a non-real upper bound, ('//TRIM(ubounds(i))//')'

      icode = 50
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

    READ(UNIT=ubounds(i),FMT=*, IOSTAT=icode) upperbound
    IF (icode /= 0) THEN
      icode = 60
      WRITE(errstr,'(A)')                                             newline//&
        TRIM(loc_src)//'Cannot convert upper bound string "'//                 &
        TRIM(ubounds(i))//'" into a real when'//                      newline//&
        'running defensive checks on '//TRIM(loc_name)//', ('//                &
        TRIM(tmp_str)//').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

  ELSE
    IF ( INDEX(values(i),'.') == 0) THEN
      WRITE(errstr,'(A)')                                             newline//&
        TRIM(loc_src)//'You are attempting to test real variable, '//          &
        TRIM(loc_name)//', ('//TRIM(tmp_str)//'), '//                 newline//&
        'with a non-real, ('//TRIM(values(i))//')'
      icode = 70
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

    READ(UNIT=values(i),FMT=*, IOSTAT=icode) condition
    IF (icode /= 0) THEN
      icode = 80
      WRITE(errstr,'(A)')                                             newline//&
        TRIM(loc_src)//'Cannot convert string "'//                             &
        TRIM(values(i))//'" into a real when running '//              newline//&
        'defensive checks on '//TRIM(ADJUSTL(name))//', ('//                   &
        TRIM(tmp_str)//').'
      CALL ereport(TRIM(srcname), icode, errstr)
    END IF

  END IF


  SELECT CASE (TRIM(operator_string))
  CASE ('<')
    IF ( var <  condition ) pass_chk = .TRUE.
  CASE ('>')
    IF ( var >  condition ) pass_chk = .TRUE.
  CASE ('>=')
    IF ( var >= condition ) pass_chk = .TRUE.
  CASE ('<=')
    IF ( var <= condition ) pass_chk = .TRUE.
  CASE (':')
    IF ( var >= lowerbound .AND. var <= upperbound ) pass_chk = .TRUE.
  CASE ('==')
    IF ( var == condition ) pass_chk = .TRUE.

  CASE DEFAULT
    WRITE(errstr,'(A)')                                                        &
      TRIM(loc_src)//'Unrecognised operator, ('//TRIM(operator_string)//')'
    icode = 90
    CALL ereport(TRIM(srcname), icode, errstr)
  END SELECT

  IF (pass_chk) EXIT

END DO

IF (PRESENT(report_pass)) THEN
  report_pass = pass_chk

  IF (pass_chk) result_chk = 'pass'
  WRITE(ummessage,'(A)')                                                       &
    'Defensive checks on '//TRIM(loc_name)//' ('//TRIM(tmp_str)//') against '//&
    TRIM(ADJUSTL(string))//' .... '//TRIM(result_chk)
  CALL umprint(ummessage,src=TRIM(srcname))

ELSE IF (.NOT. pass_chk) THEN
  WRITE(errstr,'(A)')                                                 newline//&
    TRIM(loc_src)//TRIM(loc_name)//', ('//TRIM(tmp_str)//                      &
    '), has failed to pass any '//newline//'defensive checks, '//TRIM(string)//&
    '.'//TRIM(errstr)
  icode = 100
  CALL ereport(TRIM(srcname), icode, errstr)
END IF

RETURN
END SUBROUTINE chk_real_str



!==============================================================================
!==============================================================================



SUBROUTINE chk_int_arr (var, name, intarr, report_pass, cmessage)

! Checks an integer variable against a list of criteria
! passed as an integer array in the form '[Check_A, Check_B, ... ]'.
!
! This routine only checks if the integer variable exists in the
! the array of valid values
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if ANY of the criteria
! are met. If none of the specified checks are met, ereport will
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

IMPLICIT NONE

INTEGER,          INTENT(IN) :: var       ! The integer variable to be checked
CHARACTER(LEN=*), INTENT(IN) :: name      ! The name of the test variable
INTEGER,          INTENT(IN) :: IntArr(:) ! Valid integer values

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
LOGICAL,          INTENT(INOUT), OPTIONAL :: report_pass

! Additonal comments to af to ereport in event of failure of tests
CHARACTER(LEN=*), INTENT(IN),    OPTIONAL :: cmessage


! Local variables
CHARACTER(LEN=*), PARAMETER :: routinename ='CHK_INT_ARR'

LOGICAL :: pass_chk  
INTEGER :: icode

!==========================================================
IF (.NOT. PRESENT(cmessage)) THEN
  errstr = ''
ELSE
  errstr = newline//TRIM(ADJUSTL(cmessage))
END IF

IF ( def_src == '' ) THEN
  loc_src = ''
  srcname = modname//':'//routinename
ELSE
  loc_src = 'Problem performing defensive checks ('//                          &
            TRIM(ADJUSTL(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//TRIM(ADJUSTL(def_src))
END IF

result_chk = 'FAIL'
pass_chk   = .FALSE.
loc_name   = ADJUSTL(name)

IF (ANY(var == IntArr(:))) pass_chk = .TRUE.
nconditions = SIZE(IntArr)

WRITE(fmt_str,'(A,I0,A)') '(A,I0,A,', nconditions-1, '(I0,","),I0,A)'

IF (PRESENT(report_pass)) THEN

  report_pass = pass_chk
  IF (pass_chk) result_chk = 'pass'

  WRITE(ummessage,fmt_str)                                                     &
    'Defensive checks on '//TRIM(loc_name)//' (',var,                          &
    '), against [',intarr,'] .... '//TRIM(result_chk)
  CALL umprint(ummessage,src=TRIM(srcname))

ELSE IF (.NOT. pass_chk) THEN

  WRITE(errstr,fmt_str)                                               newline//&
    TRIM(loc_src)//TRIM(loc_name)//',(',var,                                   &
    '), has failed defensive checks, [',intarr,'].'//TRIM(errstr)
  icode = 10
  CALL ereport(TRIM(srcname), icode, errstr)

END IF

RETURN
END SUBROUTINE chk_int_arr



!==============================================================================
!==============================================================================



SUBROUTINE split_conditions ()

! Support routine for chk_var, the provided string is separated and
! used to populate the module arrays:
!
! * operators
! * values
! * lbounds
! * ubounds
!
! String inputs are separated based on "," with the following recognised
! operators
!  * ==
!  * <
!  * >
!  * >=
!  * <=
!  * :
! If none of these are found, equals (==) is assumed 


IMPLICIT NONE

CHARACTER(LEN=condition_len) :: condition_str

INTEGER :: i, index_val
LOGICAL :: last_condition

! Take a local copy of the input string and 
! initialise checking arrays
operators(:) = ''
values(:)    = ''
lbounds(:)   = ''
ubounds(:)   = ''
nconditions  = 0

IF (LEN(TRIM(loc_str)) == 0) RETURN

! 1) Grab the first character,
!    to see if it is > or <, if so check for an >= or <=.
! 2) Check for range character ":"
! 3) Assume the check is for equals, i.e. "=="
DO i=1, max_n_conditions

  index_val = INDEX(loc_str, ',')

  ! Check if this is the last condition, i.e. there 
  ! should be no more "," to separate the checks.
  last_condition = (index_val == 0)

  IF (last_condition) THEN
    condition_str  = TRIM(ADJUSTL(loc_str))
  ELSE
    condition_str  = TRIM(ADJUSTL(loc_str(1:index_val-1)))

    loc_str(1:index_val) = ''
    loc_str = TRIM(ADJUSTL(loc_str))
  END IF

  IF ( condition_str(1:1) == '<' .OR. condition_str(1:1) == '>' ) THEN

    IF (condition_str(2:2) == '=') THEN
      operators(i) = condition_str(1:2)
      values(i)    = TRIM(ADJUSTL(condition_str(3:)))
    ELSE
      operators(i) = condition_str(1:1)
      values(i)    = TRIM(ADJUSTL(condition_str(2:)))
    END IF

  ELSE IF ( INDEX(condition_str,':') /= 0 ) THEN
    operators(i) = ':'
    index_val  = INDEX(condition_str, ':')
    lbounds(i) = TRIM(ADJUSTL(condition_str(:index_val-1)))
    ubounds(i) = TRIM(ADJUSTL(condition_str(index_val+1:)))
  ELSE 
    operators(i) = '=='
    IF ( INDEX(condition_str,'==') /= 0 ) THEN
      values(i) = TRIM(ADJUSTL(condition_str(3:)))
    ELSE
      values(i) = TRIM(ADJUSTL(condition_str))
    END IF
  END IF

  nconditions = nconditions + 1
  IF (last_condition) EXIT
END DO

RETURN
END SUBROUTINE split_conditions


END MODULE chk_opts_mod
