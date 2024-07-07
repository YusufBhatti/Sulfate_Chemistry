! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_LOOK----------------------------------------
!
!    Purpose: Prints out Kth 64-word PP header
!
!    Programming standard:  Unified Model Documentation Paper No 3
!
!    Documentation:  Unified Model Documentation Paper No F3
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

MODULE pr_look_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: pr_look

CHARACTER(LEN=*), PARAMETER :: ModuleName='PR_LOOK_MOD'

CONTAINS

SUBROUTINE pr_look(lookup,k)
#if defined(RECON)
USE rcf_exppx_mod, ONLY: rcf_exppx
USE rcf_ppx_info_mod, ONLY: stm_record_type
#endif

USE ppxlook_mod, ONLY: exppxc
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: lookup(:,:) ! IN Integer equivalence of PP LOOKUP
INTEGER, INTENT(IN) :: k           ! IN Field number in Look Up Table

! Local variables:---------------------------------------------
INTEGER :: icode             !Error code
INTEGER :: item              !STASH item number
INTEGER :: section           !STASH section number
INTEGER :: model             !Internal model number
INTEGER :: i                 !Index
INTEGER :: lsec,lsecd        !local seconds values (0 if lbrel<3)

CHARACTER(LEN=36) :: phrase       !Character part of PPXREF record
CHARACTER(LEN=errormessagelength) :: cmessage     !Error message

#if defined(RECON)
TYPE (stm_record_type), POINTER ::  record
#endif

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL, PARAMETER :: rdummy = 0.0

CHARACTER(LEN=*), PARAMETER :: RoutineName='PR_LOOK'

!--------------------------------------------------------------------

!  Write time and field type
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
item=MOD(lookup(42,k),1000)
section=(lookup(42,k)-item)/1000
model=lookup(45,k)

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
IF ( model == 10 ) THEN
  model = 1
END IF
icode = 0
#if defined(RECON)
record => rcf_exppx(model,section,item, NoFindArg = .TRUE.)
IF ( ASSOCIATED ( record ) ) THEN
  phrase = record % NAME
ELSE
  phrase = "Name not known"
END IF
#else
phrase=exppxc(model,section,item,                                 &
              icode,cmessage)
#endif
IF (icode /= 0) THEN
  phrase='NON-STANDARD FIELD'
END IF

WRITE(umMessage,'('' FIELD NO.'', i5,4x,a)') k,phrase
CALL umPrint(umMessage,src='pr_look')

! For lookup release < 3, lbsec(d) held day_number(s), so output zero

IF (lookup(lbrel,k) < 3) THEN
  lsec = 0
  lsecd = 0
ELSE
  lsec = lookup(lbsec,k)
  lsecd = lookup(lbsecd,k)
END IF

WRITE(umMessage,'('' VALID AT: '',2(i2.2,'':''),i2.2,''Z  '',2(i2.2,''/''),'//&
    'i4.4,''    DATA TIME: '',2(i2.2,'':''),i2.2,''Z  '',2(i2.2,''/''),'//&
    'i4.4)') lookup(lbhr,k),lookup(lbmin,k),lsec,  &
    lookup(lbdat,k),lookup(lbmon,k),lookup(lbyr,k),  &
    lookup(lbhrd,k),lookup(lbmind,k),lsecd,          &
    lookup(lbdatd,k),lookup(lbmond,k),lookup(lbyrd,k)
CALL umPrint(umMessage,src='pr_look')

!  Rest of header

CALL umPrint( &
    '   LBTIM   LBFT    LBLREC LBCODE  LBHEM  LBROW  LBNPT  LBEXT LBPACK', &
    src='pr_look')
WRITE(umMessage,'(" ",2i7, i10, 6i7)')(lookup(i,k),i=13,21)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '   LBREL   LBFC  LBCFC LBPROC   LBVC  LBRVC  LBEXP   LBBEGIN    LBNREC',&
    src='pr_look')
WRITE(umMessage,'(" ",7i7, 2i10)')(lookup(i,k),i=22,30)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '  LBPROJ  LBTYP  LBLEV LBRSVD LBRSVD LBRSVD LBRSVD   LBSRCE',&
    src='pr_look')
WRITE(umMessage,'(" ",7i7, i9)')(lookup(i,k),i=31,38)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '  DATA_TYPE     NADDR    LBUSER ITEM_CODE    LBPLEV    LBUSER MODEL_CODE',&
    src='pr_look')
WRITE(umMessage,'(" ",6i10, i11)')(lookup(i,k),i=39,45)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '         BULEV       BHULEV     BRSVD(3)     BRSVD(4)       BDATUM', &
    src='pr_look')
WRITE(umMessage,'(1p," ",5e13.4)')(TRANSFER(lookup(i,k),rdummy),i=46,50)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '          BACC         BLEV        BRLEV        BHLEV       BHRLEV',&
    src='pr_look')
WRITE(umMessage,'(1p," ",5e13.4)')(TRANSFER(lookup(i,k),rdummy),i=51,55)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '         BPLAT        BPLON         BGOR          BZY          BDY', &
    src='pr_look')
WRITE(umMessage,'(1p," ",5e13.4)')(TRANSFER(lookup(i,k),rdummy),i=56,60)
CALL umPrint(umMessage,src='pr_look')

CALL umPrint( &
    '           BZX          BDX         BMDI         BMKS',&
    src='pr_look')
WRITE(umMessage,'(1p," ",4e13.4)')(TRANSFER(lookup(i,k),rdummy),i=61,64)
CALL umPrint(umMessage,src='pr_look')
CALL umPrint('',src='pr_look')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_look

END MODULE pr_look_mod
