! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
!Purpose: To read a direct access PP file and convert it to a sequential
!         file read to be passed across to the IBM
!
SUBROUTINE readpp(len_inthd,len_realhd,len1_levdpc,len2_levdpc,   &
   len1_rowdpc,len2_rowdpc,len1_coldpc,len2_coldpc,               &
   len1_lookup,len2_lookup,len_fixhd,pp_fixhd,lookup,rookup,      &
   ppunit1,ppunit2,icode,cmessage)
USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE io
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE fort2c_interfaces, ONLY: get_file
USE errormessagelength_mod, ONLY: errormessagelength

USE poserror_mod, ONLY: poserror
IMPLICIT NONE

INTEGER ::                                                        &
     len_fixhd                                                    &
    ,len_inthd                                                    &
    ,len_realhd                                                   &
    ,len_levdpc                                                   &
    ,len_rowdpc                                                   &
    ,len_coldpc                                                   &
    ,len1_levdpc                                                  &
    ,len2_levdpc                                                  &
    ,len1_rowdpc                                                  &
    ,len2_rowdpc                                                  &
    ,len1_coldpc                                                  &
    ,len2_coldpc                                                  &
    ,len_lookup                                                   &
    ,len1_lookup                                                  &
    ,len2_lookup                                                  &
    ,len1_looknew                                                 &
    ,len2_looknew                                                 &
    ,lookup(len1_lookup,len2_lookup)                              &
    ,pp_inthd(len_inthd)                                          &
    ,pp_fixhd(len_fixhd)                                          &
    ,len_io                                                       &
    ,icode                                                        &
    ,ppunit1                                                      &
    ,ppunit2
REAL ::                                                           &
     rookup(len1_lookup,len2_lookup)                              &
    ,pp_realhd(len_realhd)                                        &
    ,pp_levdpc(len1_levdpc*len2_levdpc+1)                         &
    ,pp_rowdpc(len1_rowdpc*len2_rowdpc)                          &
    ,pp_coldpc(len1_coldpc*len2_coldpc)                          &
    ,a_io
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=filenamelength) :: outfile
! Local variables
INTEGER ::                                                        &
     start_block                                                  &
    ,nent                                                         &
    ,k                                                            &
    ,Kk                                                           &
    ,iwa                                                          &
    ,RECL                                                         &
    ,ierr                                                         &
    ,pp_len_inthd                                                 &
    ,pp_len_realhd

!---------------------------------------------------------------------
pp_len_realhd=pp_fixhd(106)
pp_len_inthd=pp_fixhd(101)

!---------------------------------------------------------------------
len_levdpc=len1_levdpc*len2_levdpc
len_rowdpc=len1_rowdpc*len2_rowdpc
len_coldpc=len1_coldpc*len2_coldpc
len_lookup=len1_lookup*len2_lookup
! The calculation of LEN_LEVDPC has PLUS 1 which is only true
! for PP headers and not model headers, hopefully the PLUS one will
! be removed as it is inconsistent)
start_block=len_fixhd+1
! ---------------------------------------------------------------
!   Read in the integer constants
! ---------------------------------------------------------------
IF (len_inthd >  0) THEN  ! Integer constants to be read in
  IF (pp_fixhd(100) /= start_block) THEN   ! Address incorrect
    CALL poserror('integer constants',start_block,100,            &
    pp_fixhd(100))
    cmessage=' READPP :  Adressing Conflict'
    icode=2
    RETURN
  END IF
  CALL buffin(ppunit1,pp_inthd,len_inthd,len_io,a_io)
  WRITE(umMessage,*)pp_inthd
  CALL umPrint(umMessage,src='readpp')
  IF (a_io /= -1.0 .OR. len_io /= len_inthd) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror(' Buffer in of Integer constants',a_io,len_io    &
  ,  len_inthd)
    cmessage='READPP : I/O error'
    icode=3
    RETURN
  END IF
  start_block=start_block+len_inthd
END IF
! ---------------------------------------------------------------
!   Read in the real constants
! ---------------------------------------------------------------
IF (len_realhd >  0) THEN  ! Real constants to be read in
  IF (pp_fixhd(105) /= start_block) THEN   ! Address incorrect
    CALL poserror('Real constants',start_block,100,               &
    pp_fixhd(105))
    cmessage=' READPP :  Adressing Conflict'
    icode=4
    RETURN
  END IF
  CALL buffin(ppunit1,pp_realhd,len_realhd,len_io,a_io)
  IF (a_io /= -1.0 .OR. len_io /= len_realhd) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror(' Buffer in of Real constants',a_io,len_io       &
    ,len_realhd)
    cmessage='READPP : I/O error'
    icode=5
    RETURN
  END IF
  start_block=start_block+len_realhd
END IF
! ---------------------------------------------------------------
!   Read in the level dependant constants
! ---------------------------------------------------------------
IF (len_levdpc >  0) THEN  ! Level dep constants to be read in
  IF (pp_fixhd(110) /= start_block) THEN   ! Address incorrect
    CALL poserror('Level depndt constants',start_block,100,       &
    pp_fixhd(110))
    cmessage=' READPP :  Adressing Conflict'
    icode=6
    RETURN
  END IF
  CALL buffin(ppunit1,pp_levdpc,len_levdpc,len_io,a_io)
  IF (a_io /= -1.0 .OR. len_io /= len_levdpc) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror(' Buffer in of Level constants',a_io,len_io      &
    ,len_levdpc)
    cmessage='READPP : I/O error'
    icode=7
    RETURN
  END IF
  start_block=start_block+len_levdpc
END IF
! ---------------------------------------------------------------
!   Read in the Row dependant constants
! ---------------------------------------------------------------
IF (len_rowdpc >  0) THEN  ! row dep constants to be read in
  IF (pp_fixhd(115) /= start_block) THEN   ! Address incorrect
    CALL poserror('Row depndt constants',start_block,100,         &
                   pp_fixhd(115))
    cmessage=' READPP :  Adressing Conflict'
    icode=10
    RETURN
  END IF
  CALL buffin(ppunit1,pp_rowdpc,len_rowdpc,len_io,a_io)
  IF (a_io /= -1.0 .OR. len_io /= len_rowdpc) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('Buffer in of Row constants',a_io,len_io,        &
                  len_rowdpc)
    cmessage='READPP : I/O error'
    icode=11
    RETURN
  END IF
  start_block=start_block+len_rowdpc
END IF
! ---------------------------------------------------------------
!   Read in the Col dependant constants
! ---------------------------------------------------------------
IF (len_coldpc >  0) THEN  ! col dep constants to be read in
  IF (pp_fixhd(120) /= start_block) THEN   ! Address incorrect
    CALL poserror('Col depndt constants',start_block,100,         &
                   pp_fixhd(120))
    cmessage=' READPP :  Adressing Conflict'
    icode=20
    RETURN
  END IF
  CALL buffin(ppunit1,pp_coldpc,len_coldpc,len_io,a_io)
  IF (a_io /= -1.0 .OR. len_io /= len_rowdpc) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror('Buffer in of Col constants',a_io,len_io,        &
                  len_coldpc)
    cmessage='READPP : I/O error'
    icode=21
    RETURN
  END IF
  start_block=start_block+len_coldpc
END IF
! ---------------------------------------------------------------
!   Read in the LOOKUP TABLE
! ---------------------------------------------------------------
IF (len_lookup >  0) THEN  ! Lookup Table to be read in
  IF (pp_fixhd(150) /= start_block) THEN   ! Address incorrect
    WRITE(umMessage,*) 'READPP : WARNING'
    CALL umPrint(umMessage,src='readpp')
    WRITE(umMessage,*) 'Conflict between start position of Lookup table'
    CALL umPrint(umMessage,src='readpp')
    WRITE(umMessage,*) 'block and pointer in fixed length header: ',     &
               'FIXHD(150) = ',pp_fixhd(150)
    CALL umPrint(umMessage,src='readpp')
    WRITE(umMessage,*) 'Current position in file = ',start_block,        &
               ' words in'
    CALL umPrint(umMessage,src='readpp')
    WRITE(umMessage,*) 'Pointer moved to ',pp_fixhd(150),' words in'
    CALL umPrint(umMessage,src='readpp')
    CALL setpos(ppunit1,pp_fixhd(150)-1,ierr)
  END IF
  CALL buffin(ppunit1,lookup,len_lookup,len_io,a_io)
  IF (a_io /= -1.0 .OR. len_io /= len_lookup) THEN
    ! DEPENDS ON: ioerror
    CALL ioerror(' Buffer in of Lookup table   ',a_io,len_io      &
    ,len_lookup)
    cmessage='READPP : I/O error'
    icode=9
    RETURN
  END IF
  start_block=start_block+len_lookup
END IF
WRITE(umMessage,*)' ARRIVED HERE  ',start_block
CALL umPrint(umMessage,src='readpp')
nent=0
DO k=1,len2_lookup
  IF (lookup(1,k) >  0) THEN
    nent=nent+1
  ELSE
    EXIT
  END IF
END DO ! K

WRITE(umMessage,*)' VALUE OF NENT   ',nent
CALL umPrint(umMessage,src='readpp')
DO k=nent-2,nent+1
  WRITE(umMessage,*)'k=',k
  CALL umPrint(umMessage,src='readpp')
  WRITE(umMessage,*) (lookup(kk,k),kk=1,44)
  CALL umPrint(umMessage,src='readpp')
END DO
!-----------------------------------------------------------------
!    OPEN NEW TARGET FIELDSFILE INITIALISING BY CALLING INITPP
!-----------------------------------------------------------------
!
!         Open named file on unit 60
!
WRITE(umMessage,*)"*** Opening new file on unit ",pPUNIT2
CALL umPrint(umMessage,src='readpp')
CALL get_file(ppunit2,outfile,filenamelength,icode)
CALL file_open(ppunit2,outfile,filenamelength,1,1,icode)
!
! DEPENDS ON: init_pp
CALL init_pp(ppunit2,'p',len1_lookup,len2_lookup,pp_fixhd,        &
             pp_inthd,pp_realhd,pp_levdpc,pp_rowdpc,pp_coldpc,    &
                       len_inthd,len_realhd,len1_levdpc,          &
             len2_levdpc,len1_rowdpc,len2_rowdpc,                 &
             len1_coldpc,len2_coldpc,pp_len_inthd,                &
             pp_len_realhd,icode,cmessage)

IF (icode /= 0) THEN
  WRITE(7,'(A,I2)')' ICODE EQUAL TO ',icode
  WRITE(7,'(A80)') cmessage
  !FixMe : ereport here?
  RETURN
END IF
len1_looknew=len1_lookup
len2_looknew=len2_lookup
! DEPENDS ON: control
CALL control(ppunit1,ppunit2,len1_looknew,len2_looknew,           &
             lookup,pp_inthd,len_inthd,                           &
             pp_fixhd,len_fixhd,icode,cmessage,nent)
IF (icode /= 0) THEN
  WRITE(7,'(A,I2)') icode
  WRITE(7,'(A80)') cmessage
  !FixMe : ereport here?
  RETURN
END IF

RETURN
END SUBROUTINE readpp
