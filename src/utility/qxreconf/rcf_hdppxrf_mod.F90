! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reads sizes of stashmaster for later space allocation.

MODULE Rcf_Hdppxrf_Mod

!  Subroutine Rcf_hdppxrf  - finds sizes of stashmaster
!
! Description:
!   Counts the number of records in a stashmaster
!
! Method:
!   Relies on the fixed formatting of a STASHmaster file to count the
!   the number of records.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_HDPPXRF_MOD'

CONTAINS

SUBROUTINE Rcf_Hdppxrf( StmsrNam , EnvDir )

USE Ereport_Mod, ONLY: &
    Ereport

USE file_manager, ONLY: &
    assign_file_unit,    &
    release_file_unit

USE filenamelength_mod, ONLY: &
    filenamelength

USE Rcf_Ppx_Info_Mod, ONLY: &
    ppxrecs

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Min,           &
    newline

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE get_env_var_mod, ONLY: get_env_var

USE um_version_mod, ONLY: &
    um_version_char

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
CHARACTER (LEN=13), INTENT(IN) :: StmsrNam   ! Name of STASHmaster
                                             ! file
CHARACTER (LEN=*), INTENT(IN)  :: EnvDir     ! Name of Env. Var for
                                             ! directory for STASHmaster

! Local parameters:
INTEGER, PARAMETER :: no_version = -9999 ! Value returned by function
                                         ! get_umversion if environment
                                         ! variable VN not set

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_HDPPXRF'
INTEGER            :: icode
INTEGER            :: Int_Model_No
INTEGER            :: FirstBlank
INTEGER            :: IOStatus
INTEGER            :: old_ppxRecs   ! Stores old value of ppxRecs
INTEGER            :: nft           ! Fortran unit number
INTEGER            :: msg           ! GCOM message id
INTEGER            :: info          ! GCOM info flag

LOGICAL            :: found_version   !Indicates presence of STM version

CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=1)  :: char1
CHARACTER (LEN=8)  :: c_stm_version  !STASHmaster version string

CHARACTER (LEN=filenamelength) :: stash_mstr
                                    ! Do. STASH master files
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------
! Most of the processing here only needs to be done by pe0 as it
! involves file manipulations; the results will be broadcast at the end.
!-------------------------------------------------------------------

pe0: IF (mype == 0 ) THEN

  IOStatus = 0

  ! Store current value of ppxRecs
  old_ppxRecs = ppxRecs

  ! Open STASHmaster file for current internal model
  !  Get directory name for STASHmaster & append rest of filename
  CALL get_env_var( EnvDir, stash_mstr)

  FirstBlank = LEN_TRIM( stash_mstr ) + 1
  stash_mstr = TRIM( stash_mstr )

  stash_mstr( FirstBlank:FirstBlank ) = '/'
  stash_mstr( FirstBlank+1:FirstBlank+13 ) = StmsrNam

  !    Open the STASHmaster file
  CALL assign_file_unit ( stash_mstr, nft, handler="fortran" )
  OPEN( UNIT=nft, FILE=stash_mstr,  ACTION='READ',                &
                  IOSTAT=IOStatus, IOMSG=iomessage )
  IF ( IOStatus > 0 ) THEN
    cmessage = 'Error opening STASHmaster file:'//stash_mstr //  &
               ' :' // newline // TRIM(iomessage)
    CALL Ereport( RoutineName, IOStatus, Cmessage )
  END IF

  IF ( PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) '!!!! STASH_MSTR ',stash_mstr
    CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
  END IF

  IF (IOStatus /= 0) THEN
    cmessage=                                                       &
      'Error opening STASHmaster file, routine HDPPXRF:'// newline//&
      'IoMsg: '//TRIM(iomessage)
    IF ( PrintStatus >= PrStatus_Min ) THEN
      WRITE(umMessage,*) &
        'HDPPXRF: Fortran Error Response = ',IOStatus, &
        ' Opening STASHmaster file ',StmsrNam
      CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
    END IF

    CALL Ereport ( RoutineName, icode, cmessage)
  END IF

  ! Check through the header section of the STASHmaster
  ! file looking for H3
  found_version = .FALSE.
  READ( nft, '(A1)') char1

  DO WHILE( char1  ==  'H' .OR. char1  ==  '#')

    IF (char1  ==  'H') THEN
      BACKSPACE nft
      READ (nft, '(1X, A1)') char1

      IF (char1  ==  '3') THEN
        !  This line starts with H3 and should indicate the
        !  STASHmaster version. The line should look like
        !     H3| UM_VERSION=4.3
        found_version = .TRUE.
        BACKSPACE nft
        READ (nft, '(15x,a8)') c_stm_version

        ! Now perform the check against the UM version
        IF (TRIM(c_stm_version)  /=  um_version_char) THEN
          WRITE (cmessage, '(A)') &
              'HDPPXRF : UM version and STASHmaster version differ'
          icode = -90

          IF ( PrintStatus >= PrStatus_Min ) THEN
            WRITE(umMessage,*) 'Version of STASHmaster file (' &
              ,TRIM(c_stm_version), &
              ') does not match UM version (' &
              ,um_version_char,') in file ',StmsrNam
            CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
          END IF

          CALL Ereport ( RoutineName, icode, cmessage)
        END IF  ! version check

      END IF  ! char1 == '3'

    END IF                 ! char1 == 'H'
    READ (nft, '(a1)') char1

  END DO


  IF (.NOT. found_version) THEN
    IF ( PrintStatus >= PrStatus_Min ) THEN
      WRITE(umMessage,*) &
         'HDPPXRF : No STASHmaster version available; Unable to'
      CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
      WRITE(umMessage,*) &
         'check against UM version for file ',StmsrNam
      CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
    END IF

    cmessage = 'HDPPXRF : No STASHmaster version available'
    icode = -100

    CALL Ereport( RoutineName, icode, cmessage )
  END IF

  !     For safety, rewind to the start of the STASH file.
  REWIND (nft)

  ! Count records - ppxRecs is counter
  Int_Model_No = 0
  DO WHILE (Int_Model_No /= -1 )
    READ( nft, '(A1)' ) char1
    IF ( char1 == '1' ) THEN
      BACKSPACE nft
      READ( nft, '(2X,I5)' ) Int_Model_No
      IF ( Int_Model_No /= -1 ) THEN
        ppxRecs = ppxRecs + 1
      END IF
    END IF
  END DO

  IF (PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) 'Hdppxrf : No of STASHmaster records in ', &
                 StmsrNam,' ',ppxRecs - old_ppxRecs
    CALL umPrint(umMessage,src='rcf_hdppxrf_mod')
  END IF

  CLOSE( UNIT = nft )
  CALL release_file_unit( nft, handler="fortran" )

END IF pe0

!-------------------------------------------------------------------
! Now to broadcast the the variables calculated
!-------------------------------------------------------------------
msg = 7001
CALL gc_ibcast( msg, 1, 0, nproc, info, ppxRecs )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_hdppxrf

END MODULE Rcf_Hdppxrf_Mod
