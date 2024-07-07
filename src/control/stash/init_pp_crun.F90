! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Initialises buffer space for lookups and headers for PP files
!   on a Crun.
!
! Method:
!    Simply allocates space, reads in the data and then attaches
!    data to the pointers in the PP buffer structures.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE Init_PP_Crun( ftn_unit, filename, &
    len1_lookup, pp_len2_lookup,              &
    filetype)

USE Model_file, ONLY:     &
    InitHeader,            &
    attachHeader,          &
    LoadHeader,            &
    setHeader,             &
    FixedHeader,           &
    InitLookups,           &
    AttachLookups,         &
    model_file_open
USE ios, ONLY:            &
    IOS_MergeLookup
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength
USE io, ONLY:              &
    setpos,                &
    readableWords,         &
    set_unit_bcast_flag,   &
    buffin,                &
    clear_unit_bcast_flag
USE io_constants, ONLY: ioFileTypeMPIIO
USE umPrintMgr
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore,  ONLY: mype
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER :: ftn_unit
INTEGER :: len1_lookup
INTEGER :: pp_len2_lookup
INTEGER :: len_fixhd

CHARACTER(LEN=filenamelength) :: filename

CHARACTER(LEN=1) :: filetype

! Local variables
INTEGER :: icode
INTEGER :: len_io
INTEGER :: isize
INTEGER :: address_ipplook
INTEGER :: dummy
INTEGER :: step
INTEGER :: i
INTEGER :: wordsRemaining

INTEGER, PARAMETER :: current_io_pe=0

REAL :: ierr

! Pointer for the fixed length header
INTEGER, POINTER :: pp_fixhd(:)

! Pointer for the lookup table
INTEGER, POINTER :: ipplook(:,:)

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_PP_CRUN'
CHARACTER(LEN=errormessagelength)           :: Cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! first we need to open the preattached file

CALL Model_File_Open(ftn_unit, filename, LEN_TRIM(filename), 1, 1, icode, &
     fileType=ioFileTypeMPIIO )
IF (icode /= 0) THEN
  cmessage='Error opening file: '//TRIM(filename)
  CALL Ereport(RoutineName, icode, cmessage)
END IF

! Is the file a pp file ?
IF ((filetype == 'p') .OR. (filetype == 'c')) THEN

  ! allocate the arrays
  ALLOCATE( ipplook(len1_lookup,pp_len2_lookup) )

  !  attach the fixed length header
  CALL initHeader(ftn_unit, FixedHeader)
  CALL loadHeader(ftn_unit, FixedHeader, icode=icode)
  pp_fixhd=>attachHeader(ftn_unit, FixedHeader)

  IF (icode /= 0) THEN
    icode = -10
    CALL Ereport(RoutineName, icode, 'Error: Failed to load fixed header.')
    ! If warning error here, skip to the end of this routine.
    GO TO 9999
  END IF

  ! If we are using an async stash server we need to propogate
  ! the fixed length header to the server
  CALL setHeader(ftn_unit,FixedHeader)

  ! init the lookup table
  step=1
  CALL initLookups(ipplook, ftn_unit, len1_lookup,               &
      pp_len2_lookup,address_ipplook,step)

  ! the position of the lookup table is pp_fixhd(150)-1
  IF (mype == 0) THEN
    address_ipplook = pp_fixhd(150)-1
  END IF

  step = 2
  CALL initLookups(ipplook, ftn_unit, dummy, dummy,              &
      address_ipplook, step)

  ! attach the lookup table
  IF (mype == current_io_pe) THEN
    ipplook=>attachLookups(ftn_unit)
  END IF

  ! read the lookup table

  CALL setpos(ftn_unit, address_ipplook, icode )
  IF (icode /= 0) THEN
    Cmessage = 'Error in setpos'
    CALL Ereport(RoutineName, icode, cmessage)
  END IF

  ! Only current_io_pe should do the read as only this PE has
  ! memory allocated for this operation!
  IF (mype == current_io_pe) THEN
    CALL set_unit_bcast_flag(ftn_unit)
    wordsRemaining=readableWords(ftn_unit)
    IF (printstatus>=prstatus_oper) THEN
      WRITE(umMessage,'(A,I3,A,I10,A)')               &
          'init_pp_crun: unit ',ftn_unit,' has ', &
          wordsRemaining,' readable words for lookup read'
      CALL umPrint(umMessage,src='init_pp_crun')
    END IF
    IF (wordsRemaining>= len1_lookup*pp_len2_lookup) THEN
      CALL buffin(ftn_unit, ipplook,                         &
          len1_lookup*pp_len2_lookup, len_io, ierr )

      CALL IOS_MergeLookup(ftn_unit, ipplook, &
          len1_lookup*pp_len2_lookup)
    ELSE
      icode = -20
      Cmessage = 'Error reading lookup table'

      CALL Ereport(RoutineName, icode, cmessage)
      ! If warning error here, skip to the end of this routine.
      GO TO 9999
    END IF
    CALL clear_unit_bcast_flag(ftn_unit)
  END IF ! current_io_pe

  ! Nullify ipplook and pp_fixhd for safety
  NULLIFY(ipplook)
  NULLIFY(pp_fixhd)

END IF ! is_ppfile

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_pp_crun
