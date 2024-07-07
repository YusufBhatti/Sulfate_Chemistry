! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE f_type_mod

IMPLICIT NONE 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='F_TYPE_MOD'

CONTAINS
!
!  SUBROUTINE F_TYPE-------------------------------------------------
!
!  Purpose:  Returns each field code and associated field length from
!            the PP header and a count of the number of fields
!            of each type.
!
! Programming standard : UMDP3
!
! Documentation: None
!
!--------------------------------------------------------------------
!
!    Arguments:------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!--------------------------------------------------------------------
SUBROUTINE f_type(lookup,len2_lookup,pp_num,n_types                &
,pp_len,pp_stash,pp_type,pp_pos,pp_ls,fixhd5,                      &
title, verbose)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
USE ppxlook_mod, ONLY: exppxc
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength
USE lookup_addresses, ONLY: lbpack

USE Packing_Codes_Mod, ONLY: PC_No_Compression, PC_BitMask_CompressType

IMPLICIT NONE


INTEGER ::                                                        &
 len2_lookup                                                      &
                         !IN 2nd dimension of LOOKUP
,n_types                                                          &
                         !IN No of separate field types in file
,lookup(64,len2_lookup)                                           &
                         !IN LOOKUP record
,pp_num(len2_lookup)                                              &
                         !OUT No of successive fields with same co
,pp_len(len2_lookup)                                              &
                         !OUT Length of field
,pp_stash(len2_lookup)                                            &
                         !OUT PP code of field
,pp_type(len2_lookup)                                             &
                         !OUT Integer/real/timeseries
,pp_pos(len2_lookup)                                              &
                         !OUT Pointer to number of PP field
,pp_ls(len2_lookup)                                               &
                         !OUT Data stored on land or sea pts
,fixhd5                  !IN Fixed header item 5 (file type)

CHARACTER(LEN=80)title
LOGICAL :: verbose


! Local variables: -----------------------------------------------------
INTEGER :: model             !Internal model number from LOOKUP

! Local arrays:---------------------------------------------------------
INTEGER ::                                                        &
 pp_xref(ppxref_codelen)  !PPXREF codes for a given section/item

! ----------------------------------------------------------------------
!    Local variables:---------------------------------------------------
INTEGER ::                                                        &
 icode                                                            &
            ! Error code
,item_code                                                        &
            ! STASH item code
,section    ! STASH section number

CHARACTER(LEN=errormessagelength) ::  cmessage
             ! Error message
CHARACTER(LEN=ppxref_charlen) :: phrase ! Name of field

INTEGER :: i,k
CHARACTER(LEN=20) :: valid(len2_lookup)
INTEGER :: fc_time(len2_lookup)
INTEGER :: lbc_levs(len2_lookup)
INTEGER :: lbproc(len2_lookup)
INTEGER :: lbtyp(len2_lookup)
REAL,ALLOCATABLE :: blev(:)
REAL,ALLOCATABLE :: bacc(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='F_TYPE'
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (verbose) THEN
  ALLOCATE(blev(len2_lookup))
  ALLOCATE(bacc(len2_lookup))
  blev(:)=0.0
  bacc(:)=0.0
END IF

! Initialise arrays
DO k=1,len2_lookup
  pp_num(k)=1
  pp_len(k)=0
  pp_stash(k)=0
  pp_type(k)=0
  pp_pos(k)=0
  pp_ls(k)=PC_No_Compression
  lbproc(k)=0
  lbtyp(k)=0
END DO

WRITE(valid(1),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)')              &
                 lookup(1,1),':',lookup(2,1),':',                 &
                 lookup(3,1),':',lookup(4,1),':',                 &
                 lookup(5,1),':',lookup(6,1)
fc_time(1)=lookup(14,1)

n_types=1
pp_len(1)=lookup(15,1)
pp_stash(1)=lookup(42,1)
lbproc(1)=lookup(25,1)
lbtyp(1)=lookup(32,1)
IF (verbose) THEN
  blev(1)=TRANSFER(lookup(52,1),blev(1))
  bacc(1)=TRANSFER(lookup(51,1),bacc(1))
END IF
pp_type(1)=lookup(39,1)
pp_pos(1)=1
IF (MOD(INT(lookup(lbpack,1)/10),10) == PC_BitMask_CompressType) THEN
  pp_ls(1)=MOD(INT(lookup(lbpack,1)/100),10)
END IF

DO k=2,len2_lookup
  IF (lookup(42,k) == lookup(42,k-1) .AND.                          &
     (lookup(18,k)*lookup(19,k)) ==                               &
     (lookup(18,k-1)*lookup(19,k-1))                              &
            .AND.                                                 &
       lookup(1,k) == lookup(1,k-1) .AND.                         &
       lookup(2,k) == lookup(2,k-1) .AND.                         &
       lookup(3,k) == lookup(3,k-1) .AND.                         &
       lookup(4,k) == lookup(4,k-1) .AND.                         &
       lookup(5,k) == lookup(5,k-1) .AND.                         &
       lookup(6,k) == lookup(6,k-1) .AND.                         &
       lookup(25,k) == lookup(25,k-1) .AND.                       &
       lookup(14,k) == lookup(14,k-1) .AND.                       &
       .NOT. verbose ) THEN
    pp_num(n_types)=pp_num(n_types)+1
  ELSE
    n_types=n_types+1
    pp_len(n_types)=lookup(15,k)
    pp_stash(n_types)=lookup(42,k)
    lbproc(n_types)=lookup(25,k)
    lbtyp(n_types)=lookup(32,k)
    IF (verbose) THEN
      blev(n_types)=TRANSFER(lookup(52,k),blev(n_types))
      bacc(n_types)=TRANSFER(lookup(51,k),bacc(n_types))
    END IF
    pp_type(n_types)=lookup(39,k)
    pp_pos(n_types)=pp_pos(n_types-1)+pp_num(n_types-1)
    WRITE(valid(n_types),'(i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)')    &
lookup(1,k),':',lookup(2,k),':',                                  &
lookup(3,k),':',lookup(4,k),':',                                  &
lookup(5,k),':',lookup(6,k)
    fc_time(n_types)=lookup(14,k)
    IF (fixhd5==5) THEN
      lbc_levs(n_types)=lookup(17,k)-100
    END IF
    IF (MOD(INT(lookup(lbpack,k)/10),10) == PC_BitMask_CompressType) THEN
      pp_ls(n_types)=MOD(INT(lookup(lbpack,k)/100),10)
    END IF
  END IF
END DO

! Print out details of fields
CALL umPrint('',src='f_type')
CALL umPrint( '************************************************'// &
    '********************************',src='f_type')
CALL umPrint('',src='f_type')
CALL umPrint(title,src='f_type')
CALL umPrint('',src='f_type')

SELECT CASE(fixhd5)
CASE (3)
  CALL umPrint('Key to output- fields file',src='f_type')
  CALL umPrint('Data stored on land/sea points',src='f_type')
  CALL umPrint('No of fields',src='f_type')
  IF (verbose) THEN
    CALL umPrint('Level value',src='f_type')
    CALL umPrint('Packing accuracy',src='f_type')
  END IF
  CALL umPrint('Length of field before unpacking (==LBLREC)', &
      src='f_type')
  CALL umPrint('Data type',src='f_type')
  CALL umPrint('STASH Code',src='f_type')
  CALL umPrint('Met08 Code',src='f_type')
  CALL umPrint('LBProc',src='f_type')
  CALL umPrint('Start of first field',src='f_type')
  CALL umPrint('Description',src='f_type')
  CALL umPrint('Validity time (yyyy:mm:dd:hh:mn:ss)',src='f_type')
  CALL umPrint('Forecast period',src='f_type')
CASE (5)
  CALL umPrint('Key to output - lbc file',src='f_type')
  CALL umPrint('Data stored on land/sea points',src='f_type')
  CALL umPrint('No of fields',src='f_type')
  CALL umPrint('Length of field before unpacking (==LBLREC)', &
      src='f_type')
  CALL umPrint('Data type',src='f_type')
  CALL umPrint('STASH Code',src='f_type')
  CALL umPrint('Start of first field',src='f_type')
  CALL umPrint('Description',src='f_type')
  CALL umPrint('Validity time (yyyy:mm:dd:hh:mn:ss)', &
      src='f_type')
CASE DEFAULT
  CALL umPrint('Key to output - dump or other',src='f_type')
  CALL umPrint('Data stored on land/sea points',src='f_type')
  CALL umPrint('No of fields',src='f_type')
  CALL umPrint('Length of field before unpacking (==LBLREC)',&
      src='f_type')
  CALL umPrint('Data type',src='f_type')
  CALL umPrint('STASH Code',src='f_type')
  CALL umPrint('Start of first field',src='f_type')
  CALL umPrint('Description',src='f_type')
  CALL umPrint('Validity time (yyyy:mm:dd:hh:mn:ss)',src='f_type')
  CALL umPrint('Forecast period',src='f_type')
END SELECT

i=1
DO k=1,n_types
  IF (lookup(42,i) == -99) THEN
    EXIT
  END IF
  phrase=' '
  item_code=MOD(lookup(42,i),1000)
  section=(lookup(42,i)-item_code)/1000
  model=lookup(45,i)
  icode = 0

  !       All diagnostics under model code of 10 are in section 20
  !       of Atmos StashMaster file.
  IF (model == 10) THEN
    model = 1
  END IF

  phrase=exppxc(model,section,item_code,                          &
              icode,cmessage)
  IF (icode /= 0) THEN
    phrase='NON-STANDARD FIELD'
  END IF
  i=i+pp_num(k)
  SELECT CASE(fixhd5)
  CASE (3)
    IF (verbose) THEN
      WRITE(umMessage, &
          '('' '',I2,I5,F8.2,F6.1,I8,I2,4I6,1x,A36,1x,a22,i3)')    &
          pp_ls(k),pp_num(k),blev(k),bacc(k),pp_len(k),                &
          pp_type(k),pp_stash(k),lbtyp(k),lbproc(k),pp_pos(k),         &
          phrase,valid(k),fc_time(k)
      CALL umPrint(umMessage,src='f_type')
    ELSE
      WRITE(umMessage, &
          '('' '',I2,I5,I8,I4,4I6,1x,A36,1x,a22,i3)')         &
          pp_ls(k),pp_num(k),pp_len(k),                             &
          pp_type(k),pp_stash(k),lbtyp(k),lbproc(k),pp_pos(k),      &
          phrase,valid(k),fc_time(k)
      CALL umPrint(umMessage,src='f_type')
    END IF
  CASE (5)
    WRITE(umMessage,'('' '',I2,I5,I8,I4,2I6,1x,A36,1x,a22)')    &
    pp_ls(k),lbc_levs(k),pp_len(k),                             &
    pp_type(k),pp_stash(k),pp_pos(k),                           &
    phrase,valid(k)
    CALL umPrint(umMessage,src='f_type')
  CASE DEFAULT
    WRITE(umMessage,'('' '',I2,I5,I8,I4,2I6,1x,A36,1x,a22,i3)') &
    pp_ls(k),pp_num(k),pp_len(k),                               &
    pp_type(k),pp_stash(k),pp_pos(k),                           &
    phrase,valid(k),fc_time(k)
    CALL umPrint(umMessage,src='f_type')
  END SELECT

END DO

IF (ALLOCATED(bacc)) DEALLOCATE(bacc)
IF (ALLOCATED(blev)) DEALLOCATE(blev)
CALL umPrint( '************************************************'// &
    '********************************',src='f_type')
CALL umPrint('',src='f_type')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE f_type
END MODULE f_type_mod
