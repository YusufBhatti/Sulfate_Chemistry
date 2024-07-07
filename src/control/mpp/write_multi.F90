! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM interface to BUFFOUT

! Subroutine Interface:
SUBROUTINE write_multi(nft,d1,isize,len_io,local_len,             &
                       lookup,fixhd12,compbuf,                    &
                       addr_info,icode,cmessage)
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE UM_ParVars
USE Decomp_DB
USE lookup_addresses
USE submodel_mod, ONLY: atmos_im
USE d1_array_mod, ONLY: d1_list_len, d1_imodl, d1_section,        &
                        d1_item, d1_grid_type
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength
USE packing_codes_mod, ONLY: PC_Cray32_Packing
IMPLICIT NONE


! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Unified Model. It is used where each process must write out a
!  local section of a global field.

! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) :: nft             !   IN : FORTRAN unit number
INTEGER, INTENT(IN) :: isize           !   IN : no. of words to write out
INTEGER, INTENT(OUT) :: len_io         !  OUT : no. of words written out
INTEGER, INTENT(OUT) :: local_len      !  OUT : size of the local field written out
INTEGER, INTENT(IN) :: lookup(64)      !   IN : LOOKUP header from dump
INTEGER, INTENT(IN) :: fixhd12         !   IN : 12th. element of fixed length header
                                       !        required for packing fields

INTEGER, INTENT(IN) ::  addr_info(d1_list_len) ! IN addressing info about field
REAL             :: d1(*)                      !   IN : Array to write out
REAL             :: compbuf(*)                 !   IN : Workspace for compressing field

INTEGER, INTENT(OUT) ::  icode                 !  OUT : Return code
CHARACTER(LEN=errormessagelength) ::   cmessage !  OUT : Error message

REAL ::  buf(isize*2)      ! Buffer for holding data to be written.
                           ! Factor of two is incase the data is
                           ! packed (ISIZE is the size on disk)
!DIR$ CACHE_ALIGN buf
INTEGER  ::  i             ! loop counter

REAL  ::  buf_icode           ! return code from buffout

INTEGER ::  orig_decomp    ! decomposition on entry
INTEGER ::  new_decomp     ! decomposition to change to

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WRITE_MULTI'

! ------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
len_io=isize
local_len=0

orig_decomp=current_decomp_type
new_decomp=orig_decomp

#if defined(FLDOP) || defined(MERGE) || defined(PPTOANC)
new_decomp=decomp_smexe
#else
IF ((addr_info(d1_imodl)  ==  atmos_im) .AND.                     &
    (orig_decomp  /=  decomp_standard_atmos)) THEN

  new_decomp=decomp_standard_atmos
END IF
#endif

IF (new_decomp  /=  orig_decomp) THEN

  CALL change_decomposition(new_decomp,icode)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'WTMULT : Error changing to decomposition ',       &
      new_decomp
    CALL umPrint(umMessage,src='write_multi')
    WRITE(umMessage,*)                                                    &
      'Attempting to write field (Model,Section,Item) ',          &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    CALL umPrint(umMessage,src='write_multi')
    icode=1
    cmessage='Failure changing decomposition'
    GO TO 9999
  END IF

END IF

! Gather the field from the local D1 array to buf


#if defined(UTILIO)
DO i=1,lookup(lblrec)
  buf(i)=d1(i)
END DO
#else
! DEPENDS ON: general_gather_field
CALL general_gather_field(                                        &
  d1,buf,local_len,lookup(lblrec),1,                              &
  addr_info,0,                                                    &
  -1,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'WTMULT : Call to GENERAL_GATHER_FIELD failed'
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'Return code was ',icode
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'Error message was ',cmessage
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'Field number ',lookup(item_code)
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'dimensions ',lookup(lbnpt),' x ',lookup(lbrow)
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'Grid type ',addr_info(d1_grid_type)
  CALL umPrint(umMessage,src='write_multi')
  WRITE(umMessage,*) 'Field was not written out'
  CALL umPrint(umMessage,src='write_multi')

  icode=300
  cmessage='Failure to gather field'
  GO TO 9999
END IF
#endif

! ------------------------------------------------------------------
! And finally the code to write the global field in array buf
! out to disk.

IF (mype  ==  0) THEN
  !       Does this field need to be compressed?
  IF (MOD((lookup(lbpack)),10)  ==  PC_Cray32_Packing) THEN
    IF (lookup(data_type)  ==  1) THEN
      ! DEPENDS ON: pack21
      CALL pack21(lookup(lblrec),buf,                             &
                  compbuf)
    END IF
  ELSE ! no compression required - just do a copy
    DO i=1,lookup(lblrec)
      compbuf(i)=buf(i)
    END DO
  END IF

  ! Now write out the global field

  IF (MOD((lookup(lbpack)),10)  ==  PC_Cray32_Packing) THEN
    IF (lookup(data_type)  ==  1) THEN
      ! Data is packed using CRAY 32 bit method - note that we need to write
      ! out 2*ISIZE 32 bit words using BUFFOUT32_F77 (as the array is 64 bit)
      ! DEPENDS ON: buffout32_f77
      CALL buffout32_f77(nft,compbuf(1:2*isize),2*isize,len_io,buf_icode)
      ! And then halve LEN_IO to satisfy tests against ISIZE
      len_io = len_io/2
    END IF
  ELSE
    ! For non-packed data
    CALL buffout(nft,compbuf(1:isize),isize,len_io,buf_icode)
  END IF
  IF ((buf_icode  /=  -1.0) .OR. (len_io  /=  isize)) THEN
    WRITE(umMessage,*) 'WTMULT : Error in call to BUFFOUT'
    CALL umPrint(umMessage,src='write_multi')
    WRITE(umMessage,*) 'LEN_IO : ',len_io
    CALL umPrint(umMessage,src='write_multi')
    WRITE(umMessage,*) 'IOSTAT : ',buf_icode
    CALL umPrint(umMessage,src='write_multi')
    WRITE(umMessage,*)                                                    &
      'Attempting to read field (Model,Section,Item) ',           &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    CALL umPrint(umMessage,src='write_multi')
    icode=400
    cmessage='Failure writing out field'
    GO TO 9999
  END IF

END IF ! am I PE 0 ?


! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

IF (MOD((lookup(lbpack)),10)  ==  PC_Cray32_Packing) THEN
  IF (lookup(data_type)  ==  1) THEN
    ! DEPENDS ON: pack21
    CALL pack21(local_len,d1,compbuf)
    ! DEPENDS ON: expand21
    CALL expand21(local_len,compbuf,d1)
  END IF
END IF

IF (new_decomp  /=  orig_decomp) THEN  ! change back

  CALL change_decomposition(orig_decomp,icode)

END IF
9999 CONTINUE  ! point to jump to if there is a failure

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE write_multi
