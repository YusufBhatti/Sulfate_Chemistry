! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: FIELDS ----------------------------------------------
!
! Purpose: To calculate fields from the Fields File such as those
! normaly derived in the Derived Printfile Program
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------

SUBROUTINE fields(pp_fixhd,len_fixhd,lenbuf,len_field,            &
                  lookup,rookup,len1_lookup,len2_lookup,nent,     &
                  stime_mod,etime_mod,nfields_mod,                &
                                        mtyp_mod,mlevs_mod,amult, &
                  stime_sel,etime_sel,nfields_sel,                &
                                        mtyp_sel,mlevs_sel,       &
                  stime_rej,etime_rej,nfields_rej,                &
                                        mtyp_rej,mlevs_rej,       &
                  stime_thi,etime_thi,nfields_thi,                &
                      mtyp_thi,mlevs_thi,ixxstep_thi,iyystep_thi, &
                  modify,SELECT,reject,thin,output_pack_type,     &
                  wind_10m,wind_10m_orog,wind_10m_scale,          &
                                                   ppunit_orog,   &
                  ppunit1,ppunit2,icode,cmessage)
USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE io
USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE fort2c_interfaces, ONLY: get_file

USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY:                                                   &
    PC_No_Packing

IMPLICIT NONE

!    Stash variables
INTEGER ::                                                        &
      len1_lookup                                                 &
,     len2_lookup                                                 &
,     lenbuf                                                      &
,     len_field                                                   &
,     stime_mod                                                   &
,     etime_mod                                                   &
,     nfields_mod                                                 &
,     mtyp_mod(NFIELDS_mod)                                       &
,     mlevs_mod(NFIELDS_mod)                                      &
,     stime_sel                                                   &
,     etime_sel                                                   &
,     nfields_sel                                                 &
,     mtyp_sel(NFIELDS_sel)                                       &
,     mlevs_sel(NFIELDS_sel)                                      &
,     stime_rej                                                   &
,     etime_rej                                                   &
,     nfields_rej                                                 &
,     mtyp_rej(NFIELDS_rej)                                       &
,     mlevs_rej(NFIELDS_rej)                                      &
,     ppunit_orog                                                 &
,     stime_thi                                                   &
,     etime_thi                                                   &
,     nfields_thi                                                 &
,     mtyp_thi(nfields_thi)                                       &
,     mlevs_thi(nfields_thi)                                      &
,     ixxstep_thi(nfields_thi)                                    &
,     iyystep_thi(nfields_thi)

REAL ::                                                           &
      amult(NFIELDS_mod)                                          &
,     wind_10m_orog                                               &
,     wind_10m_scale

INTEGER ::                                                        &
      lookup(len1_lookup,len2_lookup),                            &
      looknew(len1_lookup,len2_lookup),                           &
      bufout(lenbuf)

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=filenamelength) ::  orogfile
CHARACTER :: output_pack_type*6
LOGICAL :: last                !IN   indicates last record process
LOGICAL ::                                                        &
      modify                                                      &
,     reject                                                      &
,     SELECT                                                      &
,     wind_10m                                                    &
,     thin                                                        &
,     thin_all

INTEGER ::                                                        &
     len_fixhd                                                    &
    ,pp_fixhd(len_fixhd)                                          &
    ,icode                                                        &
    ,ppunit1                                                      &
    ,ppunit2                                                      &
    ,data_addr                                                    &
                            !    start address of data
    ,iextra(10)                                                   &
                            !IN  Used within FFREAD
    ,ier                                                          &
                            !IN  error RETURN CODE from conversion
    ,ilabel(45)                                                   &
                            !IOUT  holds integet part of lookup
    ,nent                                                         &
                            !IN  NO. ENTRIES IN OLD LOOKUP
    ,ilabel_orog(45)
REAL ::                                                           &
     rookup(len1_lookup,len2_lookup)                              &
    ,rooknew(len1_lookup,len2_lookup)                             &
    ,rlabel(19)                                                   &
                            !OUT holds real part of LOOKUP
    ,field(len_field)                                             &
    ,rlabel_orog(19)                                              &
    ,model_orog(lenbuf)
!

LOGICAL ::                                                        &
    packing                                                       &
,   READ                                                          &
,   convert

INTEGER ::                                                        &
     i                                                            &
                            ! local counter
    ,j                                                            &
                            ! local counter
    ,k                                                            &
                            ! local counter
    ,ix                                                           &
                            !
    ,il                                                           &
                            !
    ,bl                                                           &
                            !
    ,tl                                                           &
                            !
    ,iwl                                                          &
                            !
    ,nlev                                                         &
                            !
    ,iwa                                                          &
                            !
    ,iwb                                                          &
                            !
    ,ient                                                         &
                            !
    ,iproj                                                        &
                            !
    ,fct                                                          &
                            !
    ,itype                                                        &
                            !
    ,level                                                        &
                            !
    ,IDIM                                                         &
                            !
    ,len_lookup                                                   &
                            !
    ,len_io                                                       &
                            !
    ,len_buf_words                                                &
                            !
    ,num_words                                                    &
                            !
    ,pack_code                                                    &
    ,ixx                                                          &
                            ! X dimension for THIN_FIELD
    ,iyy                                                          &
                            ! Y dimension for THIN_FIELD
    ,ierr                   ! Error return from SETPOS
REAL ::                                                           &
     a_io                   !

CHARACTER ::                                                      &
     pack_type(6)*6                                               &
    ,input_pack_type*6


pack_type(1)='NONE  '
pack_type(2)='WGDOS '
pack_type(3)='CRAY32'
pack_type(4)='GRIB  '
pack_type(5)='RUNLEN'
pack_type(6)='      '

!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!

len_lookup=len1_lookup*len2_lookup

!----------------------- Section 4 ----------------------------------
!      Write to the PP file . First read in the  LOOKUP table.
!--------------------------------------------------------------------

iwa=0
CALL setpos(ppunit2,iwa,ierr)
CALL buffin(ppunit2,pp_fixhd,len_fixhd,len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= len_fixhd) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('Buffer in fixed length header',a_io,len_io,       &
                                               len_fixhd)
  icode=1
  cmessage='REPLACE: I/O error'
  RETURN
END IF
iwl=pp_fixhd(150)-1
iwa=iwl
data_addr = pp_fixhd(160)

CALL setpos(ppunit2,iwa,ierr)
CALL buffin(ppunit2,looknew,len_lookup,len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= (pp_fixhd(152)*pp_fixhd(151))) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('Buffer in Lookup table   ',a_io,len_io,           &
                                    pp_fixhd(152)*pp_fixhd(151))
  icode=1
  cmessage='Derived: I/O error in reading LOOKUP'
  RETURN
END IF

DO i=1,10
  iextra(i)=0
END DO

!     IF 10M WINDS TO BE FIXED GET MODEL OROGRAPHY FIELD FROM PP0
IF (wind_10m) THEN
  WRITE(umMessage,*) 'open unit',ppunit_orog
  CALL umPrint(umMessage,src='fields')
  CALL get_file(ppunit_orog,orogfile,filenamelength,icode)
  CALL file_open(ppunit_orog,orogfile,filenamelength,0,1,icode)
  iextra(1) = 0
  IDIM  = lenbuf
  fct   = 0
  iproj = lookup(31,1)
  itype = 73
  level = 9999
  WRITE(umMessage,*) ' read orography'
  CALL umPrint(umMessage,src='fields')
  ! DEPENDS ON: ffread
  CALL ffread(iproj,fct,itype,level,ppunit_orog,model_orog,IDIM,  &
               ilabel_orog,rlabel_orog,iextra,icode,cmessage)
  WRITE(umMessage,*) 'close unit',ppunit_orog
  CALL umPrint(umMessage,src='fields')
  CLOSE(ppunit_orog)
END IF
!
!     loop through lookup read/write all fields
iextra(1) = 1      !DO NOT UNPACK
DO ient=1,nent
  READ=.TRUE.
  convert=.FALSE.
  IDIM=lenbuf
  fct=lookup(14,ient)
  iproj=lookup(31,ient)
  itype=lookup(32,ient)
  level=lookup(33,ient)
  pack_code=MOD(lookup(lbpack,ient),10)
  input_pack_type=pack_type(pack_code+1)
  IF (input_pack_type /= output_pack_type .AND.                     &
          pack_code >  PC_No_Packing) THEN ! leave unpacked data unpacked
    convert=.TRUE.
  END IF
  WRITE(umMessage,*)' pack code=',pack_code
  CALL umPrint(umMessage,src='fields')
  WRITE(umMessage,*)input_pack_type,output_pack_type,convert
  CALL umPrint(umMessage,src='fields')
  IF (SELECT) THEN
    READ=.FALSE.
    IF (fct >= stime_sel .AND. fct <= etime_sel) THEN
      DO j=1,nfields_sel
        IF (itype == mtyp_sel(j) .AND. level == mlevs_sel(j)) THEN
          READ=.TRUE.
        END IF
      END DO
    END IF
  END IF
  IF (reject) THEN
    READ=.TRUE.
    IF (fct >= stime_rej .AND. fct <= etime_rej) THEN
      DO j=1,nfields_rej
        IF (itype == mtyp_rej(j) .AND. level == mlevs_rej(j)) THEN
          READ=.FALSE.
        END IF
      END DO
    END IF
  END IF

  WRITE(7,*) ' READ=',READ,iproj,fct,itype,level,ppunit1
  IF (READ) THEN
    ! DEPENDS ON: ffread
    CALL ffread(iproj,fct,itype,level,ppunit1,bufout,IDIM,        &
                ilabel,rlabel,iextra,icode,cmessage)
    num_words = ilabel(15)
    len_buf_words = ilabel(30)

    IF (stime_thi == -9999) THEN
      thin_all=.TRUE.
    ELSE
      thin_all=.FALSE.
    END IF
    IF (thin .OR. thin_all) THEN
      IF ((fct >= stime_thi .AND. fct <= etime_thi)                  &
          .OR. thin_all) THEN
        DO j=1,nfields_thi
          IF ((itype == mtyp_thi(j) .AND. level == mlevs_thi(j))     &
              .OR. thin_all) THEN
            iyy = ilabel(18)
            ixx = ilabel(19)
            WRITE(7,*) ' THINNING FIELD,',itype,level,fct,        &
                                     ixxstep_thi(j),iyystep_thi(j)
            WRITE(umMessage,*) ' THINNING FIELD,',itype,level,fct,       &
                                     ixxstep_thi(j),iyystep_thi(j)
            CALL umPrint(umMessage,src='fields')
            ! DEPENDS ON: thin_field
            CALL thin_field(bufout,bufout,num_words,ixx,iyy,      &
                             ixxstep_thi(j),iyystep_thi(j),       &
                               IDIM,pack_code,rlabel(18),         &
                               ilabel(15),icode,cmessage)
            len_buf_words =((num_words+511)/512)*512
            ilabel(15) = num_words
            ilabel(30) = len_buf_words
            ilabel(18) = iyy
            ilabel(19) = ixx
            rlabel(15) = rlabel(15) * iyystep_thi(j)
            rlabel(17) = rlabel(17) * ixxstep_thi(j)
          END IF
        END DO
      END IF
    END IF

    IF (modify) THEN
      IF (fct >= stime_mod .AND. fct <= etime_mod) THEN
        DO j=1,nfields_mod
          IF (itype == mtyp_mod(j) .AND. level == mlevs_mod(j)) THEN
            WRITE(7,*) ' SCALING FIELD,',itype,level,fct,amult(j)
            WRITE(umMessage,*) ' SCALING FIELD,',itype,level,fct,amult(j)
            CALL umPrint(umMessage,src='fields')
            ! DEPENDS ON: scale_field
            CALL scale_field(bufout,bufout,len_field,amult(j),    &
                            ilabel(15),IDIM,pack_code,rlabel(18), &
                            num_words, icode, cmessage)
            len_buf_words =((num_words+511)/512)*512
            ilabel(15) = num_words
            ilabel(30) = len_buf_words
          END IF
        END DO
      END IF
    END IF
    IF (wind_10m) THEN
      IF (itype == 75 .OR. itype == 76) THEN
        WRITE(umMessage,*) 'call wind fix'
        CALL umPrint(umMessage,src='fields')
        ! DEPENDS ON: wind_10m_fix
        CALL wind_10m_fix(bufout,bufout,num_words,                &
                          fct,itype,level,iproj,ppunit1,          &
                          wind_10m_scale,wind_10m_orog,           &
                          model_orog,ilabel_orog,rlabel_orog,     &
                          IDIM,pack_code,rlabel(18))
        len_buf_words =((num_words+511)/512)*512
        ilabel(15) = num_words
        ilabel(30) = len_buf_words
      END IF
    END IF

    IF (convert) THEN
      ! DEPENDS ON: conv_pack
      CALL conv_pack(ilabel,rlabel,pack_code,                     &
                     input_pack_type,output_pack_type,            &
                     bufout,IDIM,num_words,                       &
                     pp_fixhd,icode,cmessage)
      len_buf_words =((num_words+511)/512)*512
      ilabel(15) = num_words
      ilabel(30) = len_buf_words
    END IF


    ! DEPENDS ON: fld_out
    CALL fld_out(icode,cmessage,bufout,lenbuf,                    &
      len_buf_words,num_words,                                    &
      ppunit2,len1_lookup,len2_lookup,looknew,looknew,            &
      ilabel,rlabel,iwl,data_addr)

    !----------------------- Section 5 ----------------------------------
    !          Output lookup table
    !--------------------------------------------------------------------

    iwa=iwl
    CALL setpos(ppunit2,iwa,ierr)
    CALL buffout(ppunit2,looknew,len_lookup,len_io,a_io)
    !
    IF (a_io /= -1.0 .OR. len_io /=                                  &
                              (pp_fixhd(152)*pp_fixhd(151))) THEN
      ! DEPENDS ON: ioerror
      CALL ioerror('Buffer in fixed length header',a_io,len_io,   &
                    pp_fixhd(151)*pp_fixhd(152))
      icode=1
      cmessage='Derived: I/O error in writing LOOKUP'
      RETURN
    END IF
  END IF
END DO


9999 CONTINUE
RETURN
END SUBROUTINE fields
