! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Routine: READFF
!
!    Purpose: To read a   direct access PP file.
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered:
!
!    Project task:
!
!    External documentation:
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Small execs

SUBROUTINE READFF_convpp(ppunit,field,npts,entry_no,              &
ilabel,rlabel,iextra,pp_len2_lookup,len1_lookup,                  &
pp_fixhd,lookup,rookup,data_add,                                  &
model_flag,max_len_ilabel,max_len_rlabel,                         &
len_ilabel,len_rlabel,                                            &
icode,cmessage)

USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE missing_data_mod, ONLY: rmdi, imdi
USE errormessagelength_mod, ONLY: errormessagelength
USE packing_codes_mod, ONLY: PC_No_Packing, PC_Cray32_Packing
IMPLICIT NONE

!     arguments
CHARACTER ::                                                      &
     cmessage*(errormessagelength)           !OUT error message
LOGICAL ::                                                        &
     model_flag             !IN  True => Dump False =>Fieldsfile
INTEGER ::                                                        &
     len1_lookup                                                  &
                            !IN  first dimension of the lookup
    ,pp_len2_lookup                                               &
                            !IN  secnd dimension of the lookup
    ,ppunit                                                       &
                            !IN  unit no of required fieldsfile
    ,npts                                                         &
                            !IN  dimension of FIELD
    ,max_len_rlabel                                               &
                            !IN  max sixe of RLABEL
    ,max_len_ilabel                                               &
                            !IN  max sixe of ILABEL
    ,iextra(10)                                                   &
                            !IN  spare for future use
    ,data_add                                                     &
                            !IN  The word address of the data.
    ,entry_no                                                     &
                            !IN  Lookup entry no of the Field.
    ,pp_fixhd(*)                                                  &
                            !IN  PPfile fixed header
    ,lookup(len1_lookup,pp_len2_lookup)                           &
                                        !IN integer lookup
    ,len_rlabel                                                   &
                            !OUT actual size of RLABEL
    ,len_ilabel                                                   &
                            !OUT actual size of ILABEL
    ,ilabel(max_len_ilabel)                                       &
                            !OUT integer part of LOOKUP
    ,icode                  !OUT error code
REAL ::                                                           &
     field(npts)                                                  &
                            !OUT array holding final output data.
    ,rookup(len1_lookup,pp_len2_lookup)                           &
                                        !IN real lookup
    ,rlabel(max_len_rlabel) !OUT real part of LOOKUP

!     arguments for called routines
INTEGER ::                                                        &
     pack_type                                                    &
                            ! packing type N1 of LBPACK
    ,num_cray_words                                               &
                            ! number of words for field
    ,nvals                                                        &
                            ! number of points in a data field
    ,iwa                    ! Word address in call SETPOS
! ---------------------------------------------------------------------
!     LOCAL VARIABLES
INTEGER ::                                                        &
     i                                                            &
                            ! Local counter
    ,j                                                            &
                            ! Local counter
    ,length_of_data                                               &
                            ! Length of a particular field
    ,addr                                                         &
                            ! Address of a field in the data store
    ,in_lbvc                                                      &
                            ! Local copy of LBVC required to searc
    ,num_ibm_words                                                &
                            ! No of IBM words used to hold the dat
    ,pos_rlabel                                                   &
                            ! position of first REAL in PPhdr
    ,pack_type_i
REAL ::                                                           &
     amdi                   ! Missing data indicator for lookup


amdi=rookup(bmdi,entry_no)
IF (amdi /= rmdi) THEN
  WRITE(umMessage,'(A)')' NON-STANDARD MISSING DATA USED'
  CALL umPrint(umMessage,src='readff-convpp')
END IF

!
!     CALL PR_LOOK(LOOKUP(1,1),ROOKUP(1,1),ENTRY_NO)
!
!     decode LBPACK
pack_type = MOD(lookup(lbpack,entry_no),10)

!----------------------------------------------------------------------
!=== Reading a model type dump =======================================
!    A model dump has no direct addressing only relative.
!
IF (model_flag) THEN
  ! Old Format dumpfiles
  IF ((lookup(lbnrec,entry_no) == 0) .OR.                          &
    ! Prog lookups in dump before vn3.2:
              ((lookup(lbnrec,entry_no) == imdi) .AND.                      &
                                       (pp_fixhd(12) <= 301))) THEN

    IF (pack_type == PC_Cray32_Packing) THEN            ! 32 bit packing.
      num_cray_words=(lookup(lblrec,entry_no)+1)/2
    ELSE IF (pack_type > PC_No_Packing) THEN
      num_cray_words=lookup(lblrec,entry_no)/2
    ELSE
      num_cray_words=lookup(lblrec,entry_no)
    END IF
    nvals=lookup(lblrec,entry_no) ! No of data points
    addr=data_add
    IF (entry_no >  1) THEN
      DO i=1,entry_no-1
        pack_type_i = MOD(lookup(lbpack,i),10)
        IF (pack_type_i == PC_Cray32_Packing) THEN ! 32 Bit packed
          length_of_data=(lookup(lblrec,i)+1)/2
        ELSE
          length_of_data=lookup(lblrec,i)
        END IF
        addr=addr+length_of_data
      END DO
    ELSE       !  If the first entry.
      addr=data_add  !
      IF (pack_type == PC_Cray32_Packing) THEN ! 32 Bit packed
        length_of_data=(lookup(lblrec,1)+1)/2
      ELSE
        length_of_data=lookup(lblrec,1)
      END IF
      WRITE(umMessage,'(A,I0)')'  LENGTH_OF_DATA  ',length_of_data
      CALL umPrint(umMessage,src='readff-convpp')
    END IF
    iwa=addr  ! Not -1 as this is already done in dump
  ELSE
    ! New format Dumpfiles (vn4.4 onwards)

    IF (pack_type == PC_Cray32_Packing) THEN            ! 32 bit packing.
      num_cray_words=(lookup(lblrec,entry_no)+1)/2
    ELSE IF (pack_type > PC_No_Packing) THEN
      num_cray_words=lookup(lblrec,entry_no)/2
    ELSE
      num_cray_words=lookup(lblrec,entry_no)
    END IF
    iwa = lookup(lbegin,entry_no)
    nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)
  END IF
ELSE
  !=== Reading a PP type file.==========================================
  num_cray_words=lookup(lblrec,entry_no) ! PP type file
  iwa=lookup(lbegin,entry_no)
  nvals=lookup(lbrow,entry_no)*lookup(lbnpt,entry_no)             &
         +lookup(lbext,entry_no)
END IF
!==============================================================
!       WRITE(6,107) ENTRY_NO,NUM_CRAY_WORDS,NVALS
107 FORMAT(' ENTRY NO=',I0,'NUM_CRAY_WORDS= ',I0,'NVALS=',I0)
IF (npts <  num_cray_words) THEN
  icode=num_cray_words
  cmessage='READFF  npts to small ICODE holds correct value'
  GO TO 9999
END IF
icode=0
!     RETURN
! DEPENDS ON: read_rec_convpp
CALL read_rec_convpp(field,num_cray_words,iwa,ppunit,icode,cmessage)
2212 FORMAT('  FIELDS FILE NUMBER ',I0,'  ON UNIT',I0,2x,'BEING READ')
num_ibm_words=num_cray_words*2

!FIX ME WHAT IS UNIT 7 HERE?
WRITE(7,106) entry_no,                                            &
                                          ! Field No
             lookup(lbtyp,entry_no),                              &
                                          ! M08 Type
             lookup(lbfc,entry_no),                               &
                                          ! PP Field Code
             lookup(item_code,entry_no),                          &
                                          ! Stash Code
             lookup(lblev,entry_no),                              &
                                          ! M08 Level
             lookup(lbft,entry_no),                               &
                                          ! Forecast period
             lookup(lbproj,entry_no),                             &
                                          ! M08 Projection no
             num_ibm_words,                                       &
             nvals,                                               &
             pack_type                    ! Packing Code

106 FORMAT(' Field No ',I0,' M08/PP/Stash Code ',I0,'/',I0,'/',I0,&
           ' Level ',I0,' Fcst ',I0,' Proj ',I0,                  &
           ' NWords=',I0,' NVals=',I0,' Pack Type=',I0)

IF (icode == 0) THEN
  pos_rlabel=MOD(lookup(lbrel,entry_no),100)

  ! Treat lookup(45) as an integer to preserve submodel
  ! identifier in PP fields transferred between Cray and IBM.
  pos_rlabel=46


  len_rlabel=1+len1_lookup-pos_rlabel
  len_ilabel=len1_lookup-len_rlabel
  DO i=1,len_ilabel
    ilabel(i)=lookup(i,entry_no)
  END DO

  !         check for valid release number
  IF (ilabel(lbrel) <  1) THEN
    WRITE(umMessage,'(A,I0,A)')' resetting LBREL from ',ilabel(lbrel),' to 3'
    CALL umPrint(umMessage,src='readff-convpp')
    ilabel(lbrel)=3
  END IF

  !  test of header with position of reals

  !         ilabel(lbrel)= 3*1000 + pos_rlabel
  !         ilabel(lbrel)= 3
  !         ilabel(lbsrce)=pos_rlabel

  !  end of test

  DO i=1,len_rlabel
    rlabel(i)=rookup(i+pos_rlabel-1,entry_no)
  END DO
END IF
!=======================================================================
! At this point FIELD holds the data either PACKED or UN-PACKED
! Is the packing indicator set and is un-packing required? If so then
! the data is temp un-packed into a work ARRAY of length npts
IF (pack_type > PC_No_Packing) THEN                ! Is the field packed.
  IF (iextra(1) == 0) THEN  ! unpacking is required
    ! DEPENDS ON: un_pack_convpp
    CALL un_pack_convpp(pack_type,npts,field,num_cray_words,  &
                 ilabel,len_ilabel,aMDI,pp_fixhd,icode,cmessage)
    !           WRITE(7,*) ' NOW UNPACKED INTO ',ILABEL(LBLREC),' WORDS'
  END IF
ELSE IF (lookup(data_type,entry_no) == 3) THEN !Fld is logical
  ! DEPENDS ON: logical_to_real_convpp
  CALL logical_to_real_convpp(npts,field,field,nvals,         &
                       ilabel,icode,cmessage)
ELSE IF (lookup(data_type,entry_no) == 2) THEN !Fld is integer
  ! DEPENDS ON: integer_to_real_convpp
  CALL integer_to_real_convpp(npts,field,field,nvals,         &
                       ilabel,icode,cmessage)
END IF
!=======================================================================
9999 CONTINUE
100 FORMAT(//,32x,'   ARRAY        ',//,32(16f5.0/))
101 FORMAT(//,32x,'   LOOKUP       ',//,32(16i5/))
103 FORMAT('   LENIN  ',i0)
RETURN
END SUBROUTINE READFF_convpp
