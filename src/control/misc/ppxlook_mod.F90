! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc

! Purpose: Declares stashmaster (ppxref) look-up arrays used by the UM and
!          associated arrays and parameters as well as routines to poplulate,
!          manipulate and query them

MODULE ppxlook_mod

USE version_mod,  ONLY: nitemp, nsectp, ndiagp
USE cppxref_mod,  ONLY: ppxref_codelen, ppxref_charlen
USE get_env_var_mod,        ONLY: get_env_var
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,   ONLY: newline
USE submodel_mod, ONLY: fieldcalc_im, atmos_im

IMPLICIT NONE

PRIVATE ppxi, ppxc, ppxptr

! No. of STASH items per section
INTEGER, PARAMETER :: ppxref_items = nitemp

! No. of STASH sections per internal model
INTEGER, PARAMETER :: ppxref_sections = nsectp

! Max. number of non-null records in ppxref file (=4500 in version_mod)
INTEGER, PARAMETER :: num_diag_max = ndiagp

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
INTEGER :: ppxrecs

! Global arrays:
! ppxref look-up arrays
INTEGER, ALLOCATABLE :: ppxi(:,:)
INTEGER, ALLOCATABLE :: ppxptr(:,:)

CHARACTER, ALLOCATABLE :: ppxc(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PPXLOOK_MOD'

CONTAINS

SUBROUTINE read_atmos_stashmaster( stash_master_path )

USE filenamelength_mod, ONLY: filenamelength
USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen,       &
                       ppx_model_number, ppx_section_number, &
                       ppx_item_number
USE um_parcore, ONLY: mype, nproc
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY: umprint, ummessage, prstatus_normal
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE readstm_mod, ONLY: readstm

IMPLICIT NONE

! Description:
!   Reads the atmosphere STASHmaster
!
! Method:
!    Reads all the elements on PE0 and then broadcasts

  ! Subroutine arguments
CHARACTER(LEN=filenamelength), OPTIONAL, INTENT(IN) :: stash_master_path

! Local constants
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_ATMOS_STASHMASTER'

! Local variables
CHARACTER(LEN=filenamelength)   :: stash_master   ! path for stashmaster
CHARACTER(LEN=errormessagelength) :: cmessage       ! error message
CHARACTER(LEN=errormessagelength) :: iomessage      ! IO error message
CHARACTER(LEN=1)                :: char1          ! first character
CHARACTER                       :: dnam(ppxref_charlen) ! For ppxc record

INTEGER    :: codes(ppxref_codelen)               ! For ppxi record
INTEGER    :: imask(20)                           ! For version mask
INTEGER    :: errorstatus                         ! returned error status
INTEGER    :: ierr                                ! comms error code
INTEGER    :: im_ident                            ! internal model
INTEGER    :: iostatus                            ! status from I/O
INTEGER    :: rownumber                           ! location in table
INTEGER    :: section                             ! STASHmaster section
INTEGER    :: item                                ! STASHmaster item
INTEGER    :: unit_num                            ! STASHmaster file unit


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First allocate space
ALLOCATE(ppxi(ppxref_codelen, num_diag_max))
ALLOCATE(ppxc(ppxref_charlen, num_diag_max))
ALLOCATE(ppxptr(ppxref_items, 0:ppxref_sections))

ppxptr(:,:) = 0
ppxi(:,:) = 0
ppxc(:,:) = ''

! Construct filename
IF (mype == 0) THEN

  IF (PRESENT(stash_master_path)) THEN
    stash_master = stash_master_path
  ELSE
    CALL get_env_var("STASHMASTER",stash_master)
  END IF

  stash_master = TRIM(stash_master) // '/STASHmaster_A'

  ! Open the file
  CALL assign_file_unit(stash_master, unit_num, handler="fortran")
  OPEN (UNIT=unit_num, FILE=stash_master, ACTION='READ',                   &
        IOSTAT=iostatus, IOMSG=iomessage)

  IF (iostatus /= 0) THEN
    ! Error opening file - lets report it
    WRITE (ummessage,'(A,A)') 'Error opening STASHmaster file. ' //          &
               'Error: ',TRIM(iomessage)
    CALL umprint(ummessage, src=routinename)
    WRITE (ummessage, '(A)') 'FILE=' // stash_master
    CALL umprint(ummessage, src=routinename)
    errorstatus = 200
    cmessage = 'Error opening STASHmaster file:'//stash_master // ' :' //  &
                   newline // TRIM(iomessage)
    CALL ereport( routinename, errorstatus, cmessage )
  END IF

  ! Let's report what the STASHmaster is...
  WRITE (ummessage, '(A)') 'STASHmaster: ' // stash_master
  CALL umprint(ummessage, src=routinename, pe=0, level=prstatus_normal)

  ! And read the file in
  im_ident = 0
  rownumber = 0
  DO WHILE(im_ident /= -1)
    READ(unit_num, '(A1)') char1

    IF (char1 == '1') THEN

      ! A block of record to read
      BACKSPACE unit_num
      CALL readstm (imask, dnam, codes, unit_num, errorstatus, cmessage)

      im_ident = codes(ppx_model_number)
      section  = codes(ppx_section_number)
      item     = codes(ppx_item_number)

      IF (im_ident /= -1) THEN   ! not at the end record
        rownumber = rownumber + 1

        ! Check we've space still
        IF (rownumber > ndiagp) THEN
          cmessage = 'ndiagp is too small to read in STASHmaster'
          errorstatus = 300
          CALL ereport( routinename, errorstatus, cmessage )
        END IF

        ! Transfer STASHmaster record to look-up arrays
        ppxc(:,rownumber) = dnam(:)
        ppxi(:,rownumber) = codes(:)

        ! set the pointer to the record
        ppxptr(item, section) = rownumber
      END IF
    END IF ! char1 == 1
  END DO

  CLOSE(unit_num)
  CALL release_file_unit(unit_num, handler="fortran")

END IF   ! mype == 0

! Broadcast data
IF (nproc > 1) THEN
  CALL gc_ibcast(101, 1, 0, nproc, ierr, rownumber)
  CALL gc_ibcast(102, rownumber*ppxref_codelen, 0, nproc, ierr, ppxi)
  CALL gc_ibcast(103, SIZE(ppxptr), 0, nproc, ierr, ppxptr)
  CALL gc_cbcast(104, rownumber*ppxref_charlen, 0, nproc, ierr, ppxc)
END IF

ppxrecs = rownumber           ! store the size of the table

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_atmos_stashmaster

!----------------------------------------------------------------------------

SUBROUTINE compress_atmos_stashmaster( )

USE stextend_mod, ONLY: in_s
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Description:
!   Compresses the atmosphere STASHmaster
!
! Method:
!    Takes the internal STASHmaster representation and compresses it
!    to active members only based on the in_s STASH index
!

! Local constants
CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPRESS_ATMOS_STASHMASTER'


! Local variables
INTEGER, ALLOCATABLE   :: ppxptr_tmp(:,:)
INTEGER, ALLOCATABLE   :: ppxi_tmp(:,:)
CHARACTER, ALLOCATABLE :: ppxc_tmp(:,:)

INTEGER :: rownumber
INTEGER :: im_index
INTEGER :: section
INTEGER :: item
INTEGER :: i
INTEGER :: j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Make a temporary copy of the STASHmaster information
ALLOCATE( ppxc_tmp(ppxref_charlen, ppxrecs) )
ALLOCATE( ppxi_tmp(ppxref_codelen, ppxrecs) )
ALLOCATE( ppxptr_tmp(ppxref_items, 0:ppxref_sections) )

rownumber=0
im_index  = 1
DO section  = 0, ppxref_sections
  DO item = 1, ppxref_items

    ! Check whether there is a stash entry
    IF (in_s(1,im_index,section,item)  /=  0) THEN
      rownumber = rownumber + 1
      ppxptr_tmp(item, section) = rownumber
      ppxc_tmp(:, rownumber) = ppxc(:, ppxptr(item,section))
      ppxi_tmp(:, rownumber) = ppxi(:, ppxptr(item,section))
    ELSE
      ppxptr_tmp(item,section) = 0
    END IF
  END DO
END DO

! The _tmp versions now have the compressed data - reallocate space
! and copy into this.
ppxrecs = rownumber
DEALLOCATE(ppxc)
DEALLOCATE(ppxi)
ALLOCATE(ppxi(ppxref_codelen, ppxrecs))
ALLOCATE(ppxc(ppxref_charlen, ppxrecs))

DO j = 1, ppxrecs
  DO i = 1, ppxref_codelen
    ppxi(i,j) = ppxi_tmp(i,j)
  END DO
END DO

DO j = 1, ppxrecs
  DO i = 1, ppxref_charlen
    ppxc(i,j) = ppxc_tmp(i,j)
  END DO
END DO

ppxptr(:,:) = ppxptr_tmp(:,:)

DEALLOCATE(ppxc_tmp)
DEALLOCATE(ppxi_tmp)
DEALLOCATE(ppxptr_tmp)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE compress_atmos_stashmaster


! ---------------------------------------------------------------------
!  Integer Function to extract data from lookup array PPXI
!
! Function Interface:
INTEGER FUNCTION exppxi(Im_ident,section,item,element,            &
                        ErrorStatus ,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY: ereport
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE

!
! Description:
!   Extracts an individual data value from ppxref lookup array PPXI.
!
! Method:
!   The required data element is identified by the function arguments
!   Im_ident, section, item, element. The appropriate row in PPXI is
!   found from the 3-d pointer array PPXPTR as PPXPTR(i,s). The
!   address of the required element in PPXI is then given by
!   (element, row).
!
! Function arguments:
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: im_ident  ! Internal model identifier (must be atmos_im)
INTEGER, INTENT(IN) :: section   ! STASH section no.
INTEGER, INTENT(IN) :: item      ! STASH item no.
INTEGER, INTENT(IN) :: element   ! Position of required value in PPXI row

! Backwards compatibility - not used
CHARACTER(LEN=80), OPTIONAL, INTENT(INOUT) ::  cmessage
INTEGER, OPTIONAL, INTENT(INOUT) :: errorstatus !+ve = fatal error

! Local scalars
INTEGER :: row         ! Row no. in PPXI array

CHARACTER (LEN=80)            :: my_cmessage
CHARACTER (LEN=*),  PARAMETER :: RoutineName='EXPPXI'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Include file

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (im_ident /= atmos_im .AND. im_ident /= fieldcalc_im) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid model number ', im_ident
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

IF (section  <  0 .OR. section > ppxref_sections) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid section number ', section
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

IF (item <= 0 .OR. item > ppxref_items) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid item number ', item
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

! Obtain row no. in PPXI array
row = ppxptr(item,section)

! Obtain required data value
IF (row >  0) THEN
  exppxi = ppxi(element, row)
ELSE
  ! Invalid record: return an invalid value.
  exppxi = imdi
END IF

! Set a clean error code if it exists
IF ( PRESENT(errorstatus) ) errorstatus = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION exppxi


! ---------------------------------------------------------------------
!  Character Function to extract names from lookup array PPXC
!
! Function Interface:
CHARACTER(LEN=36) FUNCTION exppxc(Im_ident,section,item,               &
                             ErrorStatus ,cmessage)
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Description:
!   Extracts a diagnostic name from ppxref lookup array PPXC.
!
! Method:
!   The required name is identified by the function arguments
!   Im_ident, section, item. The appropriate row in PPXC is found
!   from the 3-d pointer array PPXPTR as PPXPTR(i,s).
!

! Function arguments:
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: im_ident  ! Internal model identifier (must be atmos_im)
INTEGER, INTENT(IN) :: section   ! STASH section no.
INTEGER, INTENT(IN) :: item      ! STASH item no.

!   Scalar arguments with intent(out):
CHARACTER(LEN=80), OPTIONAL, INTENT(INOUT) :: cmessage
INTEGER, OPTIONAL, INTENT(INOUT) :: errorstatus !+ve = fatal error

! Local scalars
INTEGER :: row         ! Row no. in PPXC array
INTEGER :: i           ! Loop counter

CHARACTER(LEN=80) :: my_cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'EXPPXC'
! Error status:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//routinename,zhook_in,zhook_handle)

IF (im_ident /= atmos_im .AND. im_ident /= fieldcalc_im) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid model number ', im_ident
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

IF (section  <  0) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid section number ', section
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

IF (item <= 0) THEN
  WRITE(my_cmessage,'(A,I5)') 'Invalid item number ', item
  errorstatus = 100
  CALL ereport(routinename, errorstatus, my_cmessage)
END IF

! Obtain row no. in PPXC array
row = ppxptr(item,section)

! Obtain required name
IF (row >  0) THEN
  DO i = 1,ppxref_charlen
    exppxc(i:i) = ppxc(i,row)
  END DO
ELSE
  ! Invalid record: return an invalid value.
  exppxc = 'ERROR: invalid STASH record'
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//routinename,zhook_out,zhook_handle)
RETURN
END FUNCTION exppxc

!-----------------------------------------------------------------------------

LOGICAL FUNCTION stashmaster_record_exists(section, item)
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
! Description:
!   Checks if a stashmaster record for a given section and item exists
!
! Method:
!    Is the pointer non-zero?

INTEGER, INTENT(IN)  :: section         ! section number
INTEGER, INTENT(IN)  :: item            ! item number

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'STASHMASTER_RECORD_EXISTS'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ppxptr(item,section) /= 0) THEN
  stashmaster_record_exists = .TRUE.
ELSE
  stashmaster_record_exists = .FALSE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION stashmaster_record_exists

!-----------------------------------------------------------------------------
SUBROUTINE set_ppxi(section, item, element, my_value)
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
! Description:
!   Sets an element of the ppxi array to a value
!
! Method:
!    trivial

INTEGER, INTENT(IN) :: section     ! section number
INTEGER, INTENT(IN) :: item        ! item number
INTEGER, INTENT(IN) :: element     ! STASHmaster element to change
INTEGER, INTENT(IN) :: my_value    ! new value for element

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SET_PPXI'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ppxi(element, ppxptr(item,section) ) = my_value

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_ppxi

END MODULE ppxlook_mod
