! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE readstm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READSTM_MOD'

CONTAINS
!
! SUBROUTINE READSTM
!   PURPOSE  TO READ A RECORD FROM THE PRE STASH-MASTERS FILE
! AND RETURN THE PPXREF CODES AND NAME OF A GIVEN DIAGNOSTIC
!
! LOGICAL COMPONENT R913
!
! PROJECT TASK: C4
!
! PROGRAMMING STANDARD  UMDP 4
!
! EXTERNAL DOCUMENT C4
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc
!-----------------------------------------------------------------
SUBROUTINE readstm                                                &
          (imask,dnam,codes,nft,icode,cmessage)

USE check_stm_codes_mod, ONLY: check_stm_codes
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cppxref_mod
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!     ARGUMENTS
INTEGER :: nft             !IN:    UNIT NUMBER STMSTS

INTEGER :: icode           !OUT: RETURN CODE
CHARACTER(LEN=errormessagelength) :: cmessage 
                           !OUT: RETURN MESSAGE IF THERE IS A FAILURE

CHARACTER(LEN=1) :: dnam(ppxref_charlen) !OUT: VARIABLE NAME FROM RECORD
INTEGER :: codes(ppxref_codelen)    !OUT: PPXREF CODES FROM RECORD
INTEGER :: imask(20)                !OUT: VERSION MASK

INTEGER :: imsk          ! Decimal equivalent of binary IMASK
INTEGER :: ii                       !LOCAL: loop mark.
INTEGER :: opcod(30)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READSTM'
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=' '

READ(nft,2010,END=3100,err=3200)                                  &
 codes(ppx_model_number)  ,                                       &
 codes(ppx_section_number),                                       &
 codes(ppx_item_number)   , dnam
2010 FORMAT(2x,3(i5,2x),36a1)

IF (codes(ppx_model_number) == -1) GO TO 9999

READ(nft,2110,END=3100,err=3200)                                  &
 codes(ppx_space_code),                                           &
 codes(ppx_ptr_code),                                             &
 codes(ppx_timavail_code),                                        &
 codes(ppx_grid_type),                                            &
 codes(ppx_lv_code),                                              &
 codes(ppx_lb_code),                                              &
 codes(ppx_lt_code),                                              &
 codes(ppx_pt_code),                                              &
 codes(ppx_pf_code),                                              &
 codes(ppx_pl_code),                                              &
 codes(ppx_lev_flag)
2110 FORMAT(2x,11(i5,2x))

! 30-digit option code read in as 6x5 digit groups instead of 4x5
READ(nft,2120,END=3100,err=3200)                                  &
(opcod(ii),ii=1,30),                                              &
(imask(ii),ii=1,20),                                              &
codes(ppx_halo_type)
2120 FORMAT(3x,30(i1),3x,20(i1),2x,i5)

codes(ppx_opt_code  )=                                            &
 opcod(30)+opcod(29)*10  +opcod(28)*100  +                        &
           opcod(27)*1000+opcod(26)*10000
codes(ppx_opt_code+1)=                                            &
 opcod(25)+opcod(24)*10  +opcod(23)*100  +                        &
           opcod(22)*1000+opcod(21)*10000
codes(ppx_opt_code+2)=                                            &
 opcod(20)+opcod(19)*10  +opcod(18)*100  +                        &
           opcod(17)*1000+opcod(16)*10000
codes(ppx_opt_code+3)=                                            &
 opcod(15)+opcod(14)*10  +opcod(13)*100  +                        &
           opcod(12)*1000+opcod(11)*10000
codes(ppx_opt_code+4)=                                            &
 opcod(10)+opcod( 9)*10  +opcod( 8)*100  +                        &
           opcod( 7)*1000+opcod( 6)*10000
codes(ppx_opt_code+5)=                                            &
 opcod( 5)+opcod( 4)*10  +opcod( 3)*100  +                        &
           opcod( 2)*1000+opcod( 1)*10000

!   Binary version mask was read into array IMASK
!   Convert version mask to decimal form IMSK
imsk = 0
DO ii=20,1,-1
  IF ((imask(ii) /= 0) .AND. (imask(ii) /= 1)) THEN
    WRITE(umMessage,*) 'READSTM: improper IMASK in user diag'
    CALL umPrint(umMessage,src='readstm')
    WRITE(umMessage,*) 'Model, Section, Item ',                      &
    codes(ppx_model_number)  ,                                &
    codes(ppx_section_number),                                &
    codes(ppx_item_number)
    CALL umPrint(umMessage,src='readstm')
  ELSE
    IF (imask(ii) == 1) THEN
      imsk=imsk+2**(20-ii)
    END IF
  END IF
END DO
!     Insert decimal value of version mask
codes(ppx_version_mask)=imsk

READ(nft,2130,END=3100,err=3200)                                  &
 codes(ppx_data_type),                                            &
 codes(ppx_dump_packing),                                         &
(codes(ii),                                                       &
 ii= ppx_pack_acc, ppx_pack_acc+ppxref_pack_profs-1)
2130 FORMAT(2x,i5,2x,i5,3x,i3,9(2x,i3))

READ(nft,2140,END=3100,err=3200)                                  &
 codes(ppx_rotate_code),                                          &
 codes(ppx_field_code),                                           &
 codes(ppx_user_code),                                            &
 codes(ppx_lbvc_code),                                            &
 codes(ppx_base_level),                                           &
 codes(ppx_top_level),                                            &
 codes(ppx_ref_lbvc_code),                                        &
 codes(ppx_cf_levelcode),                                         &
 codes(ppx_cf_fieldcode)
2140 FORMAT(2x,9(i5,2x))
3100 GO TO 9999 ! Normal completion
3200 CALL umPrint(' MESSAGE FROM ROUTINE READSTM: ',src='readstm')
WRITE(umMessage,*)' ERROR OCCURRED WHILE READING STASHmaster FILE '
CALL umPrint(umMessage,src='readstm')
cmessage=' READSTM: ERROR READING STASHMASTERS FILE'
icode=2

9999 CONTINUE
     ! If we successfully read in codes lets check them.
IF (icode == 0 ) THEN
  CALL check_stm_codes(codes(ppx_model_number),   &
                       codes(ppx_section_number), &
                       codes(ppx_item_number),    &
                       codes(ppx_grid_type),      &
                       codes(ppx_data_type),      &
                       codes(ppx_dump_packing))
END IF 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE readstm

END MODULE readstm_mod
