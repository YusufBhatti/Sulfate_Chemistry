! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: STASH
!
!  Purpose: Initialises direct access PP files at the start of
!           the run.  NB: Sequential PP files need no initialisation.
!
!  External documentation: F03 - Input and Output File Formats



SUBROUTINE init_pp ( ftn_unit,file_type_letter,                   &
                     len1_lookup,pp_len2_lookup,fixhd,            &
                     inthd,realhd,levdepc,rowdepc,coldepc,        &
                     len_inthd,                                   &
                     len_realhd,len1_levdepc,len2_levdepc,        &
                     len1_rowdepc,len2_rowdepc,                   &
                     len1_coldepc,len2_coldepc,                   &
                     pp_len_inthd, pp_len_realhd,                 &
                     icode,cmessage)

#if defined(C84_1A)
USE Model_file
#endif
USE ios, ONLY: IOS_Set_Header
USE io_configuration_mod, ONLY:   &
     io_data_alignment
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE io, ONLY:   &
    ioDiskSynchronise,             &
    buffout
USE UM_ParVars
USE missing_data_mod, ONLY: imdi, rmdi
USE errormessagelength_mod, ONLY: errormessagelength

USE poserror_mod, ONLY: poserror
IMPLICIT NONE
!
CHARACTER(LEN=1) ::                                                    &
    file_type_letter    ! IN  - File type (p-PP, b-bndry)
INTEGER ::                                                        &
    ftn_unit                                                      &
                        ! IN  - Fortran unit number
,   len1_lookup                                                   &
                        ! IN  - Size of PP header
,   pp_len2_lookup                                                &
                        ! IN  - Max allowable fields
,   len_inthd                                                     &
                        ! IN    LENGTH OF INTEGER CONSTANTS
,   len_realhd                                                    &
                        ! IN    LENGTH OF REAL CONSTANTS
,   len1_levdepc                                                  &
                        ! IN    LENGTH OF 1st Dim of lev depndt
,   len2_levdepc                                                  &
                        ! IN    LENGTH OF 2nd Dim of lev depndt
,   len1_rowdepc                                                  &
                        ! IN    LENGTH OF 1st Dim of row depndt
,   len2_rowdepc                                                  &
                        ! IN    LENGTH OF 2nd Dim of row depndt
,   len1_coldepc                                                  &
                        ! IN    LENGTH OF 1st Dim of col depndt
,   len2_coldepc                                                  &
                        ! IN    LENGTH OF 2nd Dim of col depndt
,   icode                                                         &
                        ! OUT - Error exit code
,   pp_len_inthd                                                  &
                        ! IN - Length of PP FILE integer header
,   pp_len_realhd       ! IN - Length of PP FILE real header
!
!
INTEGER ::                                                        &
    fixhd(*)                                                      &
                              ! IN    ARRAY OF FIXED CONSTANTS
,   inthd(len_inthd)
                              ! IN    ARRAY OF integer CONSTANTS
REAL ::                                                           &
    levdepc(len1_levdepc*len2_levdepc)                            &
                                        ! IN LEV DEP CONSTANTS
,   rowdepc(len1_rowdepc*len2_rowdepc)                            &
                                        ! IN ROW DEP CONSTANTS
,   coldepc(len1_coldepc*len2_coldepc)
                                        ! IN COL DEP CONSTANT
INTEGER ::                                                        &
    pp_inthd(pp_len_inthd)
                              ! OUT   ARRAY of integer constants
REAL ::                                                           &
    pp_levdepc(len1_levdepc*len2_levdepc)                         &
                                           ! OUT Level dep cts
,   pp_rowdepc(len1_rowdepc*len2_rowdepc)                         &
                                           ! OUT Row dep cts
,   pp_coldepc(len1_coldepc*len2_coldepc)  ! OUT Col dep cts

#if defined(C84_1A)
INTEGER, POINTER :: pp_fixhd(:)
#else
INTEGER ::                                                        &
    pp_fixhd(256)
#endif

!
REAL ::                                                           &
    realhd(len_realhd)                                            &
                              ! IN    ARRAY OF REAL CONSTANTS
,   pp_realhd(pp_len_realhd)  ! OUT   ARRAY OF REAL CONSTANTS
!
CHARACTER(LEN=errormessagelength) ::                              &
    cmessage            ! OUT - Error message
!
! ----------------------------------------------------------------------
!
!
!  Local variables
!
#if defined(C84_1A)
INTEGER, POINTER   :: ipplook(:,:)
INTEGER, PARAMETER :: current_io_pe=0

INTEGER            :: dummy
INTEGER            :: step
#else
INTEGER :: ipplook(len1_lookup,pp_len2_lookup)
#endif
!
!     pp_rowdepc, pp_coldepc, ipplook
INTEGER ::                                                        &
       ii,jj,iwa,ix,len_io,start_block  !
REAL :: a_io

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_PP'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! ----------------------------------------------------------------------
!  1. Reserve space
!
#if defined(C84_1A)
NULLIFY(ipplook)
step = 1
CALL initLookups(ipplook, ftn_unit, len1_lookup,                 &
                  pp_len2_lookup, dummy, step)
#else
DO ii=1,pp_len2_lookup
  DO jj=1,len1_lookup
    ipplook(jj,ii)=-99
  END DO
END DO
#endif

! ----------------------------------------------------------------------
!  1.1 Set up FIXED header record for the PP FILE
!
#if defined(C84_1A)
! Attach fixed length header
CALL initHeader(ftn_unit,FixedHeader)
pp_fixhd=>attachHeader(ftn_unit,FixedHeader)
#endif
DO ii=1,SIZE(pp_fixhd)
  pp_fixhd(ii)=fixhd(ii)
END DO
IF (fixhd(5) == 4) THEN  ! file type letter is unset for ancil.
  pp_fixhd(5)=4
ELSE IF (file_type_letter == 'p' .OR.                                  &
    file_type_letter == 'c') THEN
  pp_fixhd(5)=3
ELSE IF (file_type_letter == 'b') THEN
  pp_fixhd(5)=5
ELSE
  icode=100
  cmessage='INIT_PP  : Unknown output file type letter'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
pp_fixhd(101)=pp_len_inthd
pp_fixhd(105)=pp_fixhd(100)+pp_fixhd(101)
pp_fixhd(106)=pp_len_realhd
pp_fixhd(110)=pp_fixhd(105)+pp_fixhd(106)
pp_fixhd(111)=len1_levdepc
pp_fixhd(112)=len2_levdepc
pp_fixhd(115)=0
pp_fixhd(116)=imdi
pp_fixhd(117)=imdi
pp_fixhd(120)=0
pp_fixhd(121)=imdi
pp_fixhd(122)=imdi
pp_fixhd(125)=0
pp_fixhd(126)=imdi
pp_fixhd(127)=imdi
pp_fixhd(130)=0
pp_fixhd(131)=imdi
pp_fixhd(135)=0
pp_fixhd(136)=imdi
pp_fixhd(140)=0
pp_fixhd(141)=imdi
pp_fixhd(142)=0
pp_fixhd(143)=imdi
pp_fixhd(144)=0
pp_fixhd(145)=imdi
pp_fixhd(150)=pp_fixhd(110)+ pp_fixhd(111)*pp_fixhd(112)
IF (len2_rowdepc > 0) THEN
  pp_fixhd(115)=pp_fixhd(110)+ pp_fixhd(111)*pp_fixhd(112)
  pp_fixhd(116)=len1_rowdepc
  pp_fixhd(117)=len2_rowdepc
  pp_fixhd(150)=pp_fixhd(115)+ pp_fixhd(116)*pp_fixhd(117)
END IF
IF (len2_coldepc > 0) THEN
  pp_fixhd(120)=pp_fixhd(115)+ pp_fixhd(116)*pp_fixhd(117)
  pp_fixhd(121)=len1_coldepc
  pp_fixhd(122)=len2_coldepc
  pp_fixhd(150)=pp_fixhd(120)+ pp_fixhd(121)*pp_fixhd(122)
END IF
pp_fixhd(151)=len1_lookup
pp_fixhd(152)=pp_len2_lookup
pp_fixhd(160)=                                                    &
                   ! make sure the data starts on a sector bndry
 ((pp_fixhd(150)+pp_len2_lookup*len1_lookup-1+io_data_alignment-1)/ &
 io_data_alignment)*io_data_alignment+1

#if defined(C84_1A)
CALL SetHeader(ftn_unit,FixedHeader)
#endif

! ----------------------------------------------------------------------
!  1.2 Set up INTEGER constants record for the PP FILE
!
IF (pp_fixhd(5) <= 2) THEN !  set all values initially to MDI
  DO ii=1,pp_len_inthd
    pp_inthd(ii)=inthd(21)
  END DO
ELSE
  DO ii=1,pp_len_inthd
    pp_inthd(ii)=imdi
  END DO
END IF

pp_inthd(6)=inthd(6)
pp_inthd(7)=inthd(7)
pp_inthd(8)=inthd(8)
pp_inthd(9)=inthd(9)
pp_inthd(10)=inthd(10)
pp_inthd(12)=inthd(12)
pp_inthd(13)=inthd(13)

IF (pp_fixhd(5) /= 4) THEN !  ancils only have 15 elements in INTHD
  pp_inthd(17)=inthd(17)
  pp_inthd(24)=inthd(24)
  pp_inthd(25)=inthd(25)
  pp_inthd(28)=inthd(28)
END IF
! ----------------------------------------------------------------------
!  1.3 Set up REAL constants record for the PP FILE
!
DO ii = 1, pp_len_realhd
  pp_realhd(ii) = rmdi   ! Set all values to RMDI initially
END DO
pp_realhd(1)=realhd(1)
pp_realhd(2)=realhd(2)
pp_realhd(3)=realhd(3)
pp_realhd(4)=realhd(4)
! Set to RMDI for VR
IF (len2_rowdepc > 0 .AND. len2_coldepc > 0) THEN
  pp_realhd(1) = rmdi
  pp_realhd(2) = rmdi
  pp_realhd(3) = rmdi
  pp_realhd(4) = rmdi
END IF
pp_realhd(5)=realhd(5)
pp_realhd(6)=realhd(6)
#if defined(FLDMOD)
IF (pp_len_realhd >= 17) THEN
#endif
  IF (pp_fixhd(5) /= 4) THEN ! ancils only have 6 elements in REALHD.
    pp_realhd(16)=realhd(16)
    pp_realhd(17)=realhd(17)
  END IF
#if defined(FLDMOD)
END IF
#endif

! ----------------------------------------------------------------------
!  1.4 Set up LEVEL/ROW/COL DEPENDANT constants record for the PP FILE
!
DO ii=1,len1_levdepc*len2_levdepc
  pp_levdepc(ii)=levdepc(ii)
END DO

DO ii=1,len1_rowdepc*len2_rowdepc
  pp_rowdepc(ii)=rowdepc(ii)
END DO

DO ii=1,len1_coldepc*len2_coldepc
  pp_coldepc(ii)=coldepc(ii)
END DO

! ----------------------------------------------------------------------
!  2.1 BUFFER OUT Header Records starting with the FIXED LENGTH
!

CALL buffout(ftn_unit,pp_fixhd,SIZE(pp_fixhd),len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= SIZE(pp_fixhd)) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('bufferout of fixed length header',a_io,len_io, &
                 SIZE(pp_fixhd))
  cmessage='INIT_PP:I/O error'
  icode=1
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
start_block=SIZE(pp_fixhd)+1
! ----------------------------------------------------------------------
!  2.2 BUFFER OUT Integer Constants
!

IF (fixhd(100) >  0) THEN  ! Any integer constants to output ?

  ! Check for error in file pointers

  !        WRITE(6,*)  'START_BLOCK FIXHD(100)'
  !        WRITE(6,*)   START_BLOCK
  !        WRITE(6,*)   FIXHD(100)
  !        WRITE(6,*)   FTN_UNIT
  IF (fixhd(100) /= start_block) THEN  ! Check start address
    CALL poserror('integer constants',start_block,100,          &
    pp_fixhd(100))
    cmessage='INIT_PP:  Addressing conflict'
    icode=2
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  CALL buffout (ftn_unit,pp_inthd(1:pp_len_inthd),pp_fixhd(101),len_io,a_io)

  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= pp_fixhd(101)) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of integer constants',a_io,len_io  &
   ,pp_fixhd(101))
    cmessage='INIT_PP: I/O Error'
    icode=3
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+pp_fixhd(101)

END IF

! ----------------------------------------------------------------------
!  2.3 BUFFER OUT Real Constants
!

IF (pp_fixhd(105) >  0) THEN   ! Any real constants to output ?

  ! Check for error in file pointers

  IF (pp_fixhd(105) /= start_block) THEN
    CALL poserror('real constants',start_block,100,pp_fixhd(105))
    cmessage='INIT_PP: Addressing conflict'
    icode=4
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  CALL buffout(ftn_unit,pp_realhd(1:pp_len_realhd),pp_fixhd(106),len_io,a_io)

  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= pp_fixhd(106)) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of real constants',a_io,len_io       &
                 ,pp_fixhd(106))
    cmessage='INIT_PP: I/O Error'
    icode=5
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+pp_fixhd(106)

END IF

! ----------------------------------------------------------------------
!  2.4.1 BUFFER OUT Level Dependant Constants.
!
! Any level dependant constants ?
! note that ancils misleadingly set dep constants arrays to size 1 when not used

IF ( (pp_fixhd(112) >  0 .AND. pp_fixhd(5) /= 4) .OR. &
   (pp_fixhd(5) == 4   .AND. pp_fixhd(111) /= 1)  ) THEN 

  ! Check for error in file pointers

  IF (pp_fixhd(110) /= start_block) THEN
    CALL poserror('real constants',start_block,100,             &
                   pp_fixhd(110))
    cmessage='INIT_PP: Addressing conflict'
    icode=6
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffout (ftn_unit,pp_levdepc(1:)                           &
             ,pp_fixhd(111)*pp_fixhd(112),len_io,a_io)

  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= (pp_fixhd(111)*pp_fixhd(112)      &
       )) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of lev dep constants',a_io,len_io   &
           ,pp_fixhd(111))
    cmessage='INIT_PP: I/O Error'
    icode=7
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+ pp_fixhd(111)*pp_fixhd(112)

END IF
! ----------------------------------------------------------------------
!  2.4.2 BUFFER OUT Row Dependant Constants.
!

IF (pp_fixhd(115) >  0) THEN ! Any row dependant constants ?

  ! Check for error in file pointers

  IF (pp_fixhd(115) /= start_block) THEN
    CALL poserror('real constants',start_block,100,               &
                   pp_fixhd(115))
    cmessage='INIT_PP: Addressing conflict'
    icode=6
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffout (ftn_unit,pp_rowdepc(1:),                           &
                pp_fixhd(116)*pp_fixhd(117),len_io,a_io)

  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= (pp_fixhd(116)*pp_fixhd(117)       &
       )) THEN

    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of row dep constants',a_io,len_io,   &
                  pp_fixhd(116))
    cmessage='INIT_PP: I/O Error'
    icode=7
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+ pp_fixhd(116)*pp_fixhd(117)

END IF
! ----------------------------------------------------------------------
!  2.4.3 BUFFER OUT Col Dependant Constants.
!

IF (pp_fixhd(120) >  0) THEN ! Any col dependant constants ?

  ! Check for error in file pointers

  IF (pp_fixhd(120) /= start_block) THEN
    CALL poserror('real constants',start_block,100,               &
                   pp_fixhd(120))
    cmessage='INIT_PP: Addressing conflict'
    icode=6
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffout (ftn_unit,pp_coldepc(1:),                           &
                pp_fixhd(121)*pp_fixhd(122),len_io,a_io)

  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= (pp_fixhd(121)*pp_fixhd(122)       &
       )) THEN

    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of col dep constants',a_io,len_io ,  &
                  pp_fixhd(121))
    cmessage='INIT_PP: I/O Error'
    icode=7
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+ pp_fixhd(121)*pp_fixhd(122)

END IF
! ----------------------------------------------------------------------
!  2.5 BUFFER OUT Lookup Table
!
!     IWA= 0
!     CALL SETPOS(FTN_UNIT,3,IWA,ICODE)
IF (pp_fixhd(152) >  0) THEN

  ! Check for error in file pointers

  IF (pp_fixhd(150) /= start_block) THEN
    CALL poserror('lookup table',start_block,100,            &
         pp_fixhd(150))
    cmessage='INIT_PP: Addressing conflict'
    icode=8
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffout (ftn_unit,                                           &
                ipplook,len1_lookup*pp_len2_lookup,len_io,a_io)

  !
  ! Check for I/O errors

  IF (a_io /= -1.0 .OR. len_io /= (pp_fixhd(151)*pp_fixhd(152))) &
      THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('buffer out of PP LOOKUP TABLE ',a_io,len_io &
        ,pp_fixhd(152))
    cmessage='INIT_PP: I/O Error'
    icode=9
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Force last write to then file to commit to disk
  ! to avoid problems with continuation runs following
  ! hard failures.
  CALL ioDiskSynchronise(ftn_unit)

  start_block=start_block+(pp_fixhd(151)*pp_fixhd(152))
#if defined(C84_1A)
  !
  ! If we are the current I/O PE, we need to update our copy
  ! of the LOOKUP Table disk address
  !
  step = 2
  CALL initLookups(ipplook, ftn_unit, dummy, dummy,                &
                    pp_fixhd(150) - 1, step)
  NULLIFY(pp_fixhd)
#endif

END IF
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_pp
