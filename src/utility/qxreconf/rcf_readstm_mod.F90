! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read in a STASHmaster record

MODULE Rcf_ReadSTM_Mod

! Subroutine Rcf_ReadSTM    - Reads in a STASHmaster record
!
! Description:
!   The routine reads in a STASHmaster record.
!
! Method:
!   A STASHmaster record (format as in UMDP C4) is read in on
!   unit nft and returned in STM_tmp.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE Ereport_Mod, ONLY:   &
    Ereport

USE umPrintMgr, ONLY:    &
    newline

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READSTM_MOD'

CONTAINS

SUBROUTINE Rcf_ReadSTM( STM_tmp, nft)

USE Rcf_Ppx_Info_mod, ONLY: &
    STM_record_type,     &
    STM_OptCodeLen,      &
    STM_VerCodeLen,      &
    Ppxref_pack_profs

USE check_stm_codes_mod, ONLY: &
    check_stm_codes

IMPLICIT NONE

!     ARGUMENTS
INTEGER, INTENT(IN)             :: nft        ! UNIT NUMBER STMSTS

TYPE (STM_record_type), INTENT(OUT) :: STM_tmp ! returned data

! Local Variables
INTEGER       :: imsk          ! Decimal equivalent of binary IMASK
INTEGER       :: ii            ! loop mark.
INTEGER       :: jj
INTEGER       :: opcod( STM_OptCodeLen )     ! option codes
INTEGER       :: imask( STM_VerCodeLen )     ! version mask
INTEGER       :: icode         ! Error code
! Format spec. for variable length codes
CHARACTER(LEN=27)                   :: CodeFormat
CHARACTER(LEN=errormessagelength)   :: cmessage ! Error message
CHARACTER(LEN=errormessagelength)   :: cmessage2! Error message

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_READSTM'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE( CodeFormat, '(A,I2,A,I2,A)' )                                      &
       '(3X,', STM_OptCodeLen, '(I1),3X,', STM_VerCodeLen, '(I1),2X,I5)'

icode=0
cmessage=' '

READ(nft,'(2x,3(i5,2x),a36)',IOSTAT=icode, IOMSG=cmessage) &
      STM_tmp % model,           &
      STM_tmp % section,         &
      STM_tmp % item,            &
      STM_tmp % NAME
cmessage2 = 'ERROR Occurred while reading STASHMaster file.' // newline // &
            'IoMsg: '//TRIM(cmessage)
CALL Rcf_check_iostat(icode, RoutineName, cmessage2)

IF (STM_tmp % model /= -1) THEN

  READ(nft,'(2x,11(i5,2x))',IOSTAT=icode, IOMSG=cmessage) &
        STM_tmp % space_code,      &
        STM_tmp % ptr_code,        &
        STM_tmp % timavail_code,   &
        STM_tmp % grid_type,       &
        STM_tmp % lv_code,         &
        STM_tmp % lb_code,         &
        STM_tmp % lt_code,         &
        STM_tmp % pt_code,         &
        STM_tmp % pf_code,         &
        STM_tmp % pl_code,         &
        STM_tmp % lev_flag
  cmessage2 = 'ERROR Occurred while reading STASHMaster file.' // newline // &
              'IoMsg: '//TRIM(cmessage)
  CALL Rcf_check_iostat(icode, RoutineName, cmessage2)
  
  READ( nft, CodeFormat, IOSTAT=icode, IOMSG=cmessage ) &
        (opcod(ii) , ii=1,STM_OptCodeLen),    &
        (Imask(ii) , ii=1,STM_VerCodeLen),    &
        STM_tmp % halo_type
  cmessage2 = 'ERROR Occurred while reading STASHMaster file.' // newline // &
              'IoMsg: '//TRIM(cmessage)
  CALL Rcf_check_iostat(icode, RoutineName, cmessage2)
  
  ! Option code converted to 5-digit groups
  DO ii = 1, STM_OptCodeLen / 5
    jj = STM_OptCodeLen - (ii-1)*5
    ! ii = 1, 2, 3, 4, ...
    ! jj = 30, 25, 20, ...
    STM_tmp % opt_code(ii) =                             &
         opcod(jj) + opcod(jj-1)*10   + opcod(jj-2)*100  &
                   + opcod(jj-3)*1000 + opcod(jj-4)*10000
  END DO
  
  !   Binary version mask was read into array IMASK
  !   Convert version mask to decimal form IMSK
  imsk = 0
  DO ii=STM_VerCodeLen,1,-1
    IF (imask(ii) == 1) THEN
      imsk=imsk+2**(STM_VerCodeLen-ii)
    ELSE IF ( imask(ii) /= 0 ) THEN
      icode = 1
      WRITE(cmessage,'(A,3(A,I0))')                                  &
                     'Improper IMASK in STASHMaster' // newline,     &
                     'Model : ', STM_tmp % model,                    &
                     ', Section : ', STM_tmp % section,              &
                     ', Item ', STM_tmp % item
      CALL EReport(ModuleName//':'//RoutineName, icode, cmessage)
    END IF
  END DO
  !     Insert decimal value of version mask
  STM_tmp % version_mask = imsk
  
  READ(nft,'(2x,i5,2x,i5,3x,i3,9(2x,i3))',IOSTAT=icode, IOMSG=cmessage) &
        STM_tmp % data_type,       &
        STM_tmp % dump_packing,    &
       (STM_tmp % packing_acc(ii), &
        ii= 1, ppxref_pack_profs)
  cmessage2 = 'ERROR Occurred while reading STASHMaster file.' // newline // &
              'IoMsg: '//TRIM(cmessage)
  CALL Rcf_check_iostat(icode, RoutineName, cmessage2)
  
  READ(nft,'(2x,9(i5,2x))',IOSTAT=icode, IOMSG=cmessage) &
        STM_tmp % rotate_code,     &
        STM_tmp % field_code,      &
        STM_tmp % user_code,       &
        STM_tmp % lbvc_code,       &
        STM_tmp % base_level,      &
        STM_tmp % top_level,       &
        STM_tmp % ref_lbvc_code,   &
        STM_tmp % cf_levelcode,    &
        STM_tmp % cf_fieldcode
  cmessage2 = 'ERROR Occurred while reading STASHMaster file.' // newline // &
              'IoMsg: '//TRIM(cmessage)
  CALL Rcf_check_iostat(icode, RoutineName, cmessage2)

  CALL check_stm_codes(stm_tmp % model,     &
                       stm_tmp % section,   &
                       stm_tmp % item,      &
                       stm_tmp % grid_type, &
                       stm_tmp % data_type, &
                       stm_tmp % dump_packing)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_ReadSTM

SUBROUTINE Rcf_check_iostat(icode, Reading_Routine, message)
IMPLICIT NONE

CHARACTER(LEN=errormessagelength)   :: cmessage ! Error message


CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CHECK_IOSTAT'
INTEGER, INTENT(INOUT)         :: icode         ! Error code
CHARACTER (LEN=*), INTENT(IN)  :: Reading_Routine
CHARACTER (LEN=*), INTENT(IN)  :: message

INTEGER(KIND=jpim), PARAMETER  :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER  :: zhook_out = 1
REAL(KIND=jprb)                :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (icode /= 0) THEN
  cmessage = TRIM(message)
  CALL EReport(ModuleName//':'//Reading_Routine, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Rcf_check_iostat


END MODULE Rcf_ReadSTM_Mod
