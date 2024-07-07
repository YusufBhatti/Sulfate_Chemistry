! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate number of lookup entries in ancillary files

MODULE Calc_nlookups_Mod
IMPLICIT NONE

!  Subroutine Calc_nlookups  - Get no of lookups in ancillary files.
!
! Description:
!    Calculate number of lookup entries (NLOOKUPS) in ancillary files
!
! Method:
!    Loop over ancillary files to be opened and count the number of
!    lookup entries from the fixed headers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_NLOOKUPS_MOD'

CONTAINS

SUBROUTINE Calc_nlookups (nlookups)

USE Ancil_mod, ONLY:        &
    ancf_unitno,            &
    ancil_files,            &
    ancil_requests,         &
    num_ancil_requests,     &
    num_ancil_files

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_UMhead_Mod, ONLY: &
    LenFixHd

USE umPrintMgr, ONLY:      &
    umPrint,               &
    umMessage,             &
    PrintStatus,           &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE io

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(OUT)  :: nlookups  ! No of lookup entries

! Local variables
INTEGER       :: i               ! Loop index
INTEGER       :: j               ! Loop index
INTEGER       :: icode           ! Return code
INTEGER       :: ErrorStatus     ! Error return code
INTEGER       :: IOStatus        ! Return code from 'Inquire'
INTEGER, PARAMETER     :: filenameprovided=1
INTEGER, ALLOCATABLE   :: fixhd(:) ! Fixed Header

CHARACTER (LEN=errormessagelength)  :: Cmessage   
CHARACTER (LEN=20) :: tempFormat
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CALC_NLOOKUPS'


LOGICAL            :: l_exist       ! T if Ancillary File exists.
LOGICAL            :: l_files_found ! F is an ancillary file missing.

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! Allocate space for fixed header
!-----------------------------------------------------------------
ALLOCATE ( fixhd(LenFixHd) )

!-----------------------------------------------------------------
! Initialise nlookups
!-----------------------------------------------------------------
nlookups = 0
l_files_found = .TRUE.

!-----------------------------------------------------------------
! Open ancillary files, read in fixed header and count number
! of entries in lookup-table
!-----------------------------------------------------------------

IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='calc_nlookups_mod')
  WRITE(umMessage,'(A)') 'Ancillary Files to be opened : '
  CALL umPrint(umMessage,src='calc_nlookups_mod')
END IF


DO i=1,num_ancil_files
  IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    WRITE(umMessage, '(A,I0)') "Ancil file num: ", i
    CALL umPrint(umMessage, src='calc_nlookups_mod')
    WRITE(umMessage, '(A,A)') "Ancil filename     : ", ancil_files(i)%filename
    CALL umPrint(umMessage, src='calc_nlookups_mod')
    DO j = 1, ancil_files(i)%num_stash
      WRITE(umMessage, '(A,I0)') "Stash req = ",  ancil_files(i)%stashcodes(j)
      CALL umPrint(umMessage, src='calc_nlookups_mod')
    END DO
    WRITE(umMessage,'(A)') "--------------------------------" // &
         "-------------------------------"
    CALL umPrint(umMessage, src='calc_nlookups_mod')
  END IF

  !------------------------------------
  !    Check that ancillary file exists
  !------------------------------------
  IF (mype == 0) THEN
    INQUIRE ( FILE=ancil_files(i)%filename,         &
         EXIST=l_exist, IOSTAT=IOStatus )
    
    IF ( .NOT. l_exist .OR. IOStatus /= 0) THEN
      WRITE(umMessage,'(A)') 'Ancillary File does not exist.'
      CALL umPrint(umMessage,src='calc_nlookups_mod')
      WRITE(umMessage,'(A,A)') 'File : ',TRIM(ancil_files(i)%filename)
      CALL umPrint(umMessage,src='calc_nlookups_mod')
      ! Create format statement at runtime based on how many stashcodes are
      ! are requested from file
      WRITE(tempFormat, '(A,I0,A)') "(A,", ancil_files(i)%num_stash, "I5)"
      WRITE(umMessage, TRIM(tempFormat)) 'Stashcodes : ', &
           ancil_files(i)%stashcodes(1:ancil_files(i)%num_stash)
      CALL umPrint(umMessage,src='calc_nlookups_mod')
      l_files_found = .FALSE.
      CYCLE
    END IF
  END IF

  !---------------------------
  !    Open the ancillary file
  !---------------------------
  CALL File_Open (AncF_UnitNo,ancil_files(i)%filename ,      &
       LEN_TRIM(ancil_files(i)%filename), 0, filenameprovided,icode)
  
  IF (icode /= 0) THEN
    WRITE(umMessage,'(A)') 'Problem opening Ancillary File.'
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    WRITE(umMessage,'(A,A)') 'File : ',TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    Cmessage = 'Problem opening Ancillary file.'
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF

  !---------------------------
  !    Ensure file is at start
  !---------------------------
  CALL Setpos (AncF_UnitNo, 0, icode)
  
  IF (icode /= 0) THEN
    WRITE(umMessage,'(A)') 'Problem in SETPOS for Ancillary File.'
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    WRITE(umMessage,'(A,A)') 'File : ',TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    Cmessage = 'Error in SETPOS for Ancillary File.'
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF

  !-----------------------------------
  !    Read in the fixed length header
  !-----------------------------------
  ! DEPENDS ON: read_flh
  CALL Read_flh (AncF_UnitNo, fixhd, LenFixHd, icode, cmessage)
  
  IF (icode /= 0) THEN
    WRITE(umMessage,'(A)') 'Problem in reading fixed header from Anc File.'
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    WRITE(umMessage,'(A,A)') 'File : ',TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='calc_nlookups_mod')
    Cmessage='Problem with reading fixed header from Ancillary File.'
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF
  
  
  CALL File_Close (AncF_UnitNo,ancil_files(i)%filename ,    &
       LEN_TRIM(ancil_files(i)%filename),filenameprovided )
  
  !------------------------------------------------------
  !    Accumulate no of lookup entries in ancillary files
  !------------------------------------------------------
  nlookups = nlookups + fixhd (152)

END DO  !  Loop over ancillary files

!---------------------------------------------------------
! If a file has not been found, generate an error
!---------------------------------------------------------
IF ( .NOT. l_files_found ) THEN
  ErrorStatus = 10
  Cmessage = 'Ancillary files have not been found - Check output '//&
        'for details'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-----------------------------------------------------------------
! De-allocate space for fixed header
!-----------------------------------------------------------------
DEALLOCATE ( fixhd )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE calc_nlookups

END MODULE Calc_nlookups_Mod
