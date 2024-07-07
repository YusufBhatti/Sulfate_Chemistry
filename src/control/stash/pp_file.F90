! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE PP_FILE -----------------------------------------
!
!  Purpose:- To output a field to a PP_FILE
!
!  External documentation  C4
!
!--------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE pp_file(ppfield,icurrll,lenbuf,num_words,    &
    rmdi,comp_accrcy,pphoriz_out,unitpp,iwa,n_cols_out, &
    n_rows_out,packing,StashItem,                       &
    packing_type,current_io_pe,len_extra,               &
    srow_out,wcol_out,icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stwork_aux, ONLY: &
    flushPendingStashForUnit
USE io, ONLY: &
    buffout,   &
    setpos
USE ios, ONLY:       &
    IOS_Setpos_Stash, &
    IOS_Write_Stash_PrePacked_Data
USE ios_stash_common, ONLY: &
    isUsingAsyncStash
USE io_configuration_mod, ONLY: &
    io_field_padding
USE wgdos_packing_mod, ONLY: wgdos_compress_field

USE UM_ParVars
USE UM_ParCore, ONLY:       &
    mype
USE umPrintMgr, ONLY:       &
    umPrint,                 &
    umMessage,               &
    printstatus,             &
    prstatus_diag


USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing,                                                             &
    PC_WGDOS_packing,                                                          &
    PC_RunLength_Packing

IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage !OUT OUT MESSAGE FROM ROUTINE

INTEGER, INTENT(IN) :: &
    lenbuf,       &! Length of PP buffer
    unitpp,       &! Output PP unit number
    StashItem,    &! STASHitem=1000*section+item
    current_io_pe,&! PE which will do the I/O
    n_rows_out,   &! PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
    n_cols_out,   &! PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
    comp_accrcy,  &! Packing accuracy in power of 2
    pphoriz_out,  &! Size of output field
    icurrll,      &! Record Number
    len_extra,    &! IN size of expected extra data
    srow_out,     &! 1st southern row to output
    wcol_out       ! 1st western column to output

LOGICAL, INTENT(IN) :: &
    packing        !IN Overall Packing switch (T if pckng reqd)

INTEGER, INTENT(OUT) :: &
    packing_type   ! OUT set to 1 if WGDOS packing else set to zero.

INTEGER, INTENT(INOUT) :: &
    num_words,    &! Number of 64 bit words worth of data
    icode          ! Return code from routine


REAL ::                        &
    ppfield(pphoriz_out),      & !INOUT Array to store PP data
    bufout(lenbuf+len_extra),  & !OUTPUT PP buffer (ROUNDED UP)
    rmdi                         !IN     Missing Data Indicator

INTEGER ::                     &
    int_bufout(lenbuf+len_extra) ! Integer representation of bufout

!   WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
! -------------------------------------------------------------------

! -------------------------------------------------------------------
!   MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!--------------------------------------------------------------------

!    DEFINE LOCAL VARIABLES
INTEGER ::          &
    iwa,            &! record number
    ocode,          &! copy of the input value of icode
    len_buf_words,  &! num_words rounded by 512 and actually
    num_out,        &! number of compressed (32 bit) words
    len_io,         &! not used, but needed for buffout call
    jj               ! Local counter
LOGICAL :: uv
REAL :: ix           ! Return value from unit command

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PP_FILE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!    Check if an error has already been encountered, and get out
!    if it has.
ocode = 0
IF (icode  >   0) THEN
  GO TO 999
ELSE IF (icode  <   0) THEN
  ocode = icode
  icode = 0
END IF

!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!======================================================================
!    At this point packing, if required, will be done using the WGDOS
!    method of packing.
packing_type=PC_No_Packing
! Note the value of -26 corresponds to -15 (F) in ppxref.
! The packing acuracy is scaled to allow greater accuracy.
! Packing will only be attempted if there are at least 2 points per row
! in the PPfield.

IF (packing .AND. n_cols_out  >=  2) THEN
  ! 1. Climate wgdos packing has been selected for current file
  !    stream via namelist/gui
  IF (comp_accrcy  >   -99) THEN
    ! 2. STASH packing profile for the field is set.
    packing_type = PC_WGDOS_packing
  END IF
END IF
!
IF (packing_type == PC_WGDOS_packing) THEN
  CALL wgdos_compress_field(ppfield,pphoriz_out,int_bufout,lenbuf,         &
                            n_cols_out,num_out,comp_accrcy,rmdi,           &
                            stashitem,icode)

  bufout = TRANSFER(int_bufout, bufout)
  num_words=(num_out+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
                          ! wgdos_compress_field returns the number of IBM
                          ! words needed to hold the packed data.
                          !                  ~~~
  len_buf_words=((num_words+io_field_padding-1)/io_field_padding)*    &
   io_field_padding
  IF (printstatus >= prstatus_diag) THEN
    WRITE(umMessage,'(1X,A,I7,A,I7)')'PP_FILE: Field Packed from',pphoriz_out,&
    ' to', num_words
    CALL umPrint(umMessage,src='pp_file')
  END IF
ELSE IF (packing_type == PC_RunLength_Packing) THEN
  IF (mype  ==  current_io_pe) THEN
    ! DEPENDS ON: runlen_encode
    CALL runlen_encode(ppfield,pphoriz_out,bufout,pphoriz_out,    &
                     num_out,rmdi,icode,cmessage)
    ! Size of run length encoded data is greater than unpacked
    ! field therefore leave field unpacked.
    IF (num_out  >=  pphoriz_out) THEN
      packing_type = PC_No_Packing
      DO jj=1,pphoriz_out
        bufout(jj) = ppfield(jj)
      END DO
      num_words=pphoriz_out
      len_buf_words=lenbuf
    ELSE
      num_words=num_out
      len_buf_words=((num_words+io_field_padding-1)/io_field_padding)*&
      io_field_padding
    END IF

  END IF
ELSE  ! No packing required.
  IF (printstatus >= prstatus_diag) THEN
  WRITE(umMessage,'(1X,A,I7,I7)')'PP_FILE: Output field with grid dimensions',&
    n_rows_out,n_cols_out
    CALL umPrint(umMessage,src='pp_file')
  END IF
  DO jj=1,pphoriz_out
    bufout(jj) = ppfield(jj)
  END DO

  num_words=pphoriz_out
  len_buf_words=lenbuf
END IF

IF (isUsingAsyncStash() .AND. icurrll>=0) THEN
  !! If we are here we received a stash field object thats
  !  not using the async-queue method, and we are calling the
  !  fallback code. Hence we should purge anything for the
  !  unit that may be cached and not yet sent before we send
  !  more data.
  CALL flushPendingStashForUnit(unitpp)
END IF

IF (mype  ==  current_io_pe) THEN
  DO jj=num_words+1,len_buf_words
    bufout(jj)= 0.0
  END DO
  IF (isUsingAsyncStash() .AND. icurrll>=0) THEN
    ! This is an alternate route for stash fields, where we use data
    ! on the backend IOS to keep a track of file locations. If we are
    ! executing this code it means that some fields are written via async
    ! and this is a fallback for those that are not - in this case the
    ! file locations are unknown and we need to use the stash aware setpos
    ! with the field number....
    CALL IOS_Setpos_Stash(unitpp,icurrll,icode)
    CALL IOS_Write_Stash_PrePacked_Data(unitpp,icurrll, &
        bufout(1:len_buf_words))
    len_io=len_buf_words
    ix = -1.0
  ELSE
    CALL setpos(unitpp,iwa,icode)
    CALL buffout(unitpp,bufout,len_buf_words,len_io,ix)
  END IF

  IF (ix /= -1.0 .OR. len_io /= len_buf_words) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('Buffer out Data Field',ix,len_io,                 &
                  len_buf_words)
    cmessage='PPFILE  : I/O error - PP Data Field Output'
    icode=7
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF ! (mype  ==  current_io_pe)
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
IF (icode  ==  0 .AND. ocode  /=  0) THEN
  icode = ocode
END IF
999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pp_file
