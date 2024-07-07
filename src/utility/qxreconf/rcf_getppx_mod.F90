! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reads stashmaster file into an array of records

MODULE Rcf_Getppx_Mod
IMPLICIT NONE
!  Subroutine Rcf_Getppx - read in stashmaster file into array
!
! Description:
!   Reads a series of records from a stashmaster file into the
!   internal array data-structure and updating the ppxptr reference
!   pointer as required
!
! Method:
!   Normal and User stashmaster entries dealt with seperately
!   Normal entries are assumed to be increasing in section/item no
!   Whereas user stashmaster requires more work to ensure that
!   records are inserted/replaced at the correct point.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Counters for array position
INTEGER, PRIVATE, TARGET, SAVE   :: RowNumber    = 0
INTEGER, PRIVATE, TARGET, SAVE   :: RowNumber_in = 0

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GETPPX_MOD'

CONTAINS

SUBROUTINE Rcf_Getppx( StmsrNam, EnvDir )

USE file_manager, ONLY: &
    assign_file_unit,    &
    release_file_unit

USE filenamelength_mod, ONLY: &
    filenamelength

USE Ereport_Mod, ONLY:   &
    Ereport

USE Rcf_Ppx_Info_Mod, ONLY:  &
    ppxptr,               &
    ppxrecs,              &
    STM_record_type,      &
    STM_record

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    newline

USE Rcf_ReadSTM_Mod, ONLY:   &
    Rcf_ReadSTM

USE Rcf_BcastSTM_Mod, ONLY:  &
    Rcf_BcastSTM

USE UM_ParCore, ONLY:   &
    mype,               &
    nproc

USE version_mod, ONLY:  &
    NDiagP

USE errormessagelength_mod, ONLY: errormessagelength
USE get_env_var_mod, ONLY: get_env_var

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
CHARACTER (LEN=13), INTENT(IN)  :: StmsrNam    ! Names of stash master
                                               ! files
CHARACTER (LEN=*), INTENT(IN)   :: EnvDir      ! Env variable for
                                               ! directory name
! Local Variables
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_GETPPX'
INTEGER            :: i, id         ! Loop counters
INTEGER            :: ErrorStatus   ! Error return code
INTEGER            :: hashcount
INTEGER            :: ifil,irec     ! Do.
INTEGER            :: nftppxref     ! Unit number for STASHmaster
INTEGER            :: nftstmstu     ! Unit number for UserSTASH
INTEGER            :: iostatus
INTEGER            :: Im_index      !
INTEGER            :: Im_ident      !
INTEGER            :: Section       !
INTEGER            :: Item          !
INTEGER            :: LModel  ,dm
INTEGER            :: LSection,ds
INTEGER            :: LItem   ,di
INTEGER            :: USTrow
INTEGER            :: FirstBlank    !Used to append Upsm filename to dir
INTEGER            :: ri            ! Row index
INTEGER            :: Info          ! For GCOM
INTEGER            :: msg           ! For GCOM tag

CHARACTER (LEN=filenamelength) :: stash_mstr
                                    ! Do. STASH master files

CHARACTER (LEN=errormessagelength):: cmessage      ! Error return message
CHARACTER (LEN=errormessagelength):: iomessage
CHARACTER (LEN=1)  :: char1
CHARACTER (LEN=80), POINTER :: ustsfils_ptr(:) ! Pointer for ustsfils/_in

LOGICAL            :: overwrite     ! T if a system stash master record
                                    !is being overwritten by a user rec

TYPE (STM_record_type) :: STM_tmp   ! Temporary record
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!- End of header -------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0
IOStatus    = 0

IF ( .NOT. ASSOCIATED ( STM_record ) ) THEN
  ALLOCATE( STM_record( ppxRecs ) )

  ! Initialise important data
  STM_record(:) % RowIndex = 0
END IF

!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
!----------------------------------------------------------------------
IF ( ppxRecs > NDiagP) THEN
  WRITE(umMessage,*) 'ERROR: no. of diags. requested exceeds max'
  CALL umPrint(umMessage,src='rcf_getppx_mod')
  WRITE(umMessage,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDiagP
  CALL umPrint(umMessage,src='rcf_getppx_mod')
  Errorstatus=104
  Cmessage = 'GETPPX: ppxRecs > NDIAGP '

  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
!----------------------------------------------------------------------

  !---------------------------------------------------------------------
  !Read in records from STASHmaster for current internal model
  !--------------------------------------------------------------------
  ! Get a file unit

IF ( mype == 0 ) THEN

  !Open STASHmaster file for current internal model
  !  Get directory name for STASHmaster & append rest of filename

  CALL get_env_var( EnvDir, stash_mstr )

  FirstBlank = LEN_TRIM( stash_mstr ) + 1
  stash_mstr = TRIM( stash_mstr )

  stash_mstr(FirstBlank:FirstBlank)='/'
  stash_mstr(FirstBlank+1:FirstBlank+13)=StmsrNam

  !   Open the STASHmaster file
  CALL assign_file_unit( stash_mstr, nftppxref, handler="fortran" )
  OPEN(UNIT=nftppxref, FILE=stash_mstr, ACTION='READ',                     &
       IOSTAT=IOStatus, IOMSG=iomessage)
  IF (IOStatus /= 0) THEN
    cmessage = 'Error opening STASHmaster file:'//stash_mstr // &
               ' :' // newline // TRIM(iomessage)
    CALL Ereport( RoutineName, IOStatus, Cmessage )
  END IF

END IF   ! pe0

Im_ident = 0
DO WHILE (Im_ident /= -1)
  IF ( mype == 0 ) THEN
    READ( nftppxref, '(A1)' ) char1
  END IF

  msg = 9001
  CALL gc_cbcast( msg, 1, 0, nproc, info, char1 )

  IF (char1 == '1') THEN
    !Read block of records
    IF ( mype == 0 ) THEN
      BACKSPACE nftppxref
      CALL Rcf_ReadStm( STM_tmp, nftppxref)
    END IF

    CALL Rcf_BcastSTM( STM_tmp, 0 )

    Im_ident = STM_tmp % model
    Section  = STM_tmp % section
    Item     = STM_tmp % item

    IF (Im_ident /= -1) THEN      ! Not end of file
      Im_index= 1
      RowNumber = RowNumber + 1   !   Increment row number
      ppxptr(Im_ident,Section,Item) = RowNumber
      STM_record( RowNumber ) = STM_tmp

      ! Set row index - indicates values of model,sec,item for
      ! this row
      STM_record( RowNumber ) % RowIndex  =  Im_ident*100000 &
                                             + Section *1000 &
                                             + Item

      !   Set flag to indicate record originated from ppxref file
      STM_record( RowNumber ) % OriginFlag = 'P'

      IF (RowNumber > ppxRecs) THEN
        ErrorStatus = 10
        WRITE ( Cmessage , '(A, I8)' ) &
          ' PPXI row number exceeds total no. of ppx records ', &
          RowNumber

        CALL Ereport( RoutineName, ErrorStatus, Cmessage )

      END IF
    ELSE      ! Im_ident /= -1
      IF ( mype == 0 ) THEN
        CLOSE( nftppxref )
        CALL release_file_unit( nftppxref, handler="fortran" )
      END IF  ! mype == 0
    END IF    ! Im_ident /= -1
  END IF      ! Char == '1'
END DO        ! Im_ident /= -1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Getppx

END MODULE Rcf_Getppx_Mod
