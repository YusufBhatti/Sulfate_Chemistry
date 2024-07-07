! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: FFREAD and others (see below) ----------------------------
!
! Purpose: To read a   direct access PP file  and convert it to a
! sequential file read to be passed across to the IBM
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! Logical components covered:
!
! Project task:
!
! External documentation:
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------
!
!    IEXTRA(1) == 0  ! unpacking is required
!    IEXTRA(2) == 0  ! lookup table entry deleted after access.
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
SUBROUTINE ffreadb(iproj,fct,itype,int_level,ppunit,field,IDIM,   &
ilabel,rlabel,iextra,pp_len2_lookup,len1_lookup,                  &
iwa,len1_record,maxff,record,pfno,pp_fixhd,lookup,rookup,data_add,&
icode,cmessage)
USE ereport_mod, ONLY: ereport
USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: rmdi 
USE errormessagelength_mod, ONLY: errormessagelength
USE Packing_Codes_Mod, ONLY: PC_No_Packing, PC_Cray32_Packing
IMPLICIT NONE


CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER ::                                                        &
     maxff                                                        &
                            !IN  Max number of opened files
    ,len1_lookup                                                  &
                            !IN  first dimension of the lookup
    ,len1_record                                                  &
                            !IN  First dimension of record
    ,pfno                                                         &
                            !IN  No of fields files opened
    ,pp_len2_lookup                                               &
                            !IN  secnd dimension of the lookup
    ,iproj                                                        &
                            !IN  map projection of data to read
    ,fct                                                          &
                            !IN  forecast period in hours
    ,itype                                                        &
                            !IN  M08 FIELD field type
    ,int_level                                                    &
                            !IN  LEVEL code (could be real)
    ,ppunit                                                       &
                            !IN  unit no of required fieldsfile
    ,IDIM                                                         &
                            !IN  dimension of FIELD
    ,ilabel(45)                                                   &
                            !OUT holds integer part of LOOKUP
    ,iextra(10)                                                   &
                            !IN  spare for future use
    ,icode                                                        &
                            !OUT return code
    ,lenbuf                                                       &
                            !OUT input buffer length for data
    ,num_cray_words                                               &
                            !OUT no of values in an input field
    ,data_add                                                     &
                            !IN  The word address of the data.
    ,nvals                  !OUT The num of points in a data field
INTEGER ::                                                        &
     pp_fixhd(*),                                                 &
                                        !IN  PPfile fixed header
     lookup(len1_lookup,pp_len2_lookup) !OUT integer lookup
REAL ::                                                           &
     field(IDIM)                                                  &
                            !OUT array holding final output data.
    ,rookup(len1_lookup,pp_len2_lookup)                           &
                                          !OUT real lookup
    ,rlabel(19)                                                   &
                            !OUT holds real part of LOOKUP
    ,real_level             !IN  LEVEL code (could be real)
LOGICAL ::                                                        &
     record(len1_record,maxff) !IN Record of the field no read
!     LOCAL VARIABLES
REAL ::                                                           &
     amdi                   ! Missing data indicator
INTEGER ::                                                        &
     i                                                            &
                            ! local counter
    ,j                      ! local counter
INTEGER ::                                                        &
     ix                                                           &
                            ! used as a dummy variable in UNIT
    ,iwa                                                          &
                            ! Word address in call SETPOS
    ,ik                                                           &
                            ! Word address in call SETPOS
    ,length_of_data                                               &
                            ! Length of a particular field
    ,addr                                                         &
                            ! Address of a field in the data store
    ,in_lbvc                                                      &
                            ! Local copy of LBVC required to searc
    ,pack_type                                                    &
                            ! Packing type N1 of LBPACK
    ,pack_type_i

!
!     DO I=112,112
!       CALL PR_LOOK(LOOKUP(1,1),ROOKUP(1,1),I)
!     ENDDO
!
!----------------------------------------------------------------------
!   Search for the required FIELD
!----------------------------------------------------------------------
IF (iextra(3) == 0) THEN ! Search on LBTYP/LBLEV/LBPROJ/LBFT
  DO  i=1,pp_len2_lookup
    IF (itype == lookup(lbtyp,i)) THEN
      IF (int_level == lookup(lblev,i)) THEN
        IF (iproj == lookup(lbproj,i)) THEN
          IF (fct == lookup(lbft,i)) THEN
            ik=i
            GO TO 3
          END IF
        END IF
      END IF
    END IF
  END DO
ELSE
  in_lbvc=iextra(3)
  IF (iextra(4) == 0) THEN  !  IEXTRA(3) HAS LBVC so search on LBVC
    IF (int_level == 8888 .OR. int_level == 9999) THEN ! special lev
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbtyp,i)) THEN
          IF (int_level == lookup(lblev,i)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                ik=i
                GO TO 3
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE  ! Not a special level so additional search on LBVC
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbtyp,i)) THEN
          IF (int_level == lookup(lblev,i)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                IF (in_lbvc == lookup(lbvc,i)) THEN
                  ik=i
                  GO TO 3
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    END IF   ! End of special level block
    !----------------------------------------------------------------------
    !     Search now on BLEV ie REAL_LEVEL except for special levels      C
    !     and Data on model levels (BLEV would contain BK so would need   C
    !     to search in this case on LBLEV).For special level search on    C
    !     just LBVC LBFT LBPROJ and LBTYP.For model level convert the     C
    !     real input value to integer and search as above plus LBLEV.     C
    !     Note a special level cannot have an LBVC of 9                   C
    !----------------------------------------------------------------------
  ELSE IF (iextra(4) == 1) THEN !  IEXTRA(4) is not zero.
    ! DEPENDS ON: level_rlevel
    CALL level_rlevel(int_level,int_level,real_level)
    IF (real_level <= 0.0) THEN !  Special level indicated.
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbtyp,i)) THEN
          IF (iproj == lookup(lbproj,i)) THEN
            IF (fct == lookup(lbft,i)) THEN
              IF (in_lbvc == lookup(lbvc,i)) THEN
                ik=i
                GO TO 3
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE IF (in_lbvc == 9) THEN  ! model level data
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbtyp,i)) THEN
          int_level=real_level+0.0000001  !
          !               IF(REAL_LEVEL <= (ROOKUP(BLEV,I)+0.0001).AND.
          !    *          REAL_LEVEL >= (ROOKUP(BLEV,I)-0.0001)) THEN
          ! That MOD is only for un-corrected model dumps
          IF (int_level == lookup(lblev,i)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                IF (in_lbvc == lookup(lbvc,i)) THEN
                  ik=i
                  GO TO 3
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE        ! not model level data or a special level
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbtyp,i)) THEN
          IF (real_level <= (rookup(blev,i)+0.0001) .AND.           &
          real_level >= (rookup(blev,i)-0.0001)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                IF (in_lbvc == lookup(lbvc,i)) THEN
                  ik=i
                  GO TO 3
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    END IF
    !----------------------------------------------------------------------
    !     Search now on BLEV ie REAL_LEVEL except for special levels      C
    !     and Data on model levels (BLEV would contain BK so would need   C
    !     to search in this case on LBLEV).For special level search on    C
    !     just LBVC LBFT LBPROJ and LBTYP.For model level convert the     C
    !     real input value to integer and search as above plus LBLEV.     C
    !     Note a special level cannot have an LBVC of 9                   C
    !----------------------------------------------------------------------
  ELSE IF (iextra(4) == 2) THEN !  IEXTRA(4) is not zero.
    ! DEPENDS ON: level_rlevel
    CALL level_rlevel(int_level,int_level,real_level)
    IF (real_level <= 0.0) THEN !  Special level indicated.
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbfc,i)) THEN
          IF (iproj == lookup(lbproj,i)) THEN
            IF (fct == lookup(lbft,i)) THEN
              IF (in_lbvc == lookup(lbvc,i)) THEN
                ik=i
                GO TO 3
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE IF (in_lbvc == 9) THEN  ! model level data
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbfc,i)) THEN
          int_level=real_level+0.0000001  !
          IF (int_level == lookup(lblev,i)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                IF (in_lbvc == lookup(lbvc,i)) THEN
                  ik=i
                  GO TO 3
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    ELSE        ! not model level data or a special level
      DO  i=1,pp_len2_lookup
        IF (itype == lookup(lbfc,i)) THEN
          IF (real_level == rookup(blev,i)) THEN
            IF (iproj == lookup(lbproj,i)) THEN
              IF (fct == lookup(lbft,i)) THEN
                IF (in_lbvc == lookup(lbvc,i)) THEN
                  ik=i
                  GO TO 3
                END IF
              END IF
            END IF
          END IF
        END IF
      END DO
    END IF
  END IF    ! IEXTRA(3) == 0 IF block
END IF      ! IEXTRA(4) == 0 IF block

IF (iextra(4) >= 1 .AND. iextra(3) /= 0) THEN
  WRITE(umMessage,112) itype,real_level,iproj,fct
  CALL umPrint(umMessage,src='ffreadb')
ELSE
  WRITE(umMessage,104) itype,int_level,iproj,fct
  CALL umPrint(umMessage,src='ffreadb')
END IF
104  FORMAT('  FIELD NOT FOUND FOR ITYPE,INT_LEVEL,IPROJ,FCT',4i5)
112  FORMAT('  FIELD NOT FOUND FOR ITYPE,REAL_LEVEL,IPROJ,FCT',i5,     &
           f7.1,2i5)
icode=1
cmessage=' FFREAD  field not found'
GO TO 9999
3 CONTINUE
!=== Decode LBPACK code
pack_type = MOD(lookup(lbpack,ik),10)
!=== Reading a model type dump =======================================
!    A model dump has no direct addressing only relative.
IF (lookup(lbnrec,ik) == 0) THEN ! A model dump
  IF (pack_type == PC_Cray32_Packing) THEN          ! Is the field packed.
    num_cray_words=lookup(lblrec,ik)/2
  ELSE
    num_cray_words=lookup(lblrec,ik)
  END IF
  nvals=lookup(lblrec,ik) ! No of data points
  addr=data_add
  DO i=1,ik-1
    pack_type_i = MOD(lookup(lbpack,i),10)
    IF (pack_type_i == PC_Cray32_Packing) THEN ! 32 Bit packed
      length_of_data=lookup(lblrec,i)/2
    ELSE
      length_of_data=lookup(lblrec,i)
    END IF
    addr=addr+length_of_data
  END DO
  iwa=addr-1
ELSE
  !=== Reading a PP type file.==========================================
  num_cray_words=lookup(lblrec,ik) ! PP type file
  iwa=lookup(29,ik)
  nvals=lookup(44,ik)
END IF
record(ik,pfno)=.TRUE.   ! Record which the no of the field read
lenbuf=lookup(lbnrec,ik) !
!==============================================================
IF (iextra(4) >= 1 .AND. iextra(3) /= 0) THEN
  WRITE(7,110) itype,real_level,iproj,fct,ik,num_cray_words,nvals
ELSE
  WRITE(7,106) itype,int_level,iproj,fct,ik,num_cray_words,nvals
END IF
106 FORMAT(' FIELD ','ITYPE=',i3,' LEVEL=',i5,' PROJ=',i4,' FCST=',   &
    i5,' FIELD NO',i4,' NWORDS=',i5,' NVALS=',i5)
110 FORMAT(' FIELD FOUND','ITYPE=',i4,'LEVEL=',f7.1,'PROJ=',i4,'FCST=' &
    ,i5,'FIELD NO',i4,'NWORDS=',i5,'NVALS=',i5)
IF (IDIM <  num_cray_words) THEN
  icode=num_cray_words
  cmessage='FFREAD  Idim to small ICODE holds correct value'
  GO TO 9999
END IF
icode=0
!     RETURN
! DEPENDS ON: read_rec_ffread1a
CALL READ_REC_ffread1a(field,num_cray_words,iwa,ppunit,icode,cmessage)
2212 FORMAT('  FIELDS FILE NUMBER ',i2,'  ON UNIT',i2,2x,'BEING READ')
!     CLOSE(PPUNIT)
IF (icode == 0) THEN
  DO i=1,45
    ilabel(i)=lookup(i,ik)
  END DO
  DO i=1,19
    rlabel(i)=rookup(i+45,ik)
  END DO
END IF
!=======================================================================
! At this point FIELD holds the data either PACKED or UN-PACKED
! Is the packing indicator set and is un-packing required? If so then
! the data is temp un-packed into a work ARRAY of length IDIM
IF (pack_type > PC_No_Packing) THEN               ! Is the field packed.
  IF (iextra(1) == 0) THEN  ! unpacking is required
    !           get missing data indicator from pp header
    amdi = rookup(bmdi,ik)
    !           compare with MDI
    IF (amdi /= rmdi) THEN
      WRITE(umMessage,*)' WARNING non-standard MDI in use'
      CALL umPrint(umMessage,src='ffreadb')
    END IF
    ! DEPENDS ON: un_pack_ffread1a
    CALL UN_PACK_ffread1a(pack_type,IDIM,field,num_cray_words   &
                ,ilabel,amdi,pp_fixhd,icode,cmessage)
  END IF
ELSE IF (lookup(data_type,ik) == 3) THEN   !Fld is logical
  ! DEPENDS ON: logical_to_real_ffread1a
  CALL LOGICAL_TO_REAL_ffread1a(IDIM,field,field,nvals,         &
                       ilabel,icode,cmessage)
ELSE IF (lookup(data_type,ik) == 2) THEN   !Fld is integer
  ! DEPENDS ON: integer_to_real_ffread1a
  CALL INTEGER_TO_REAL_ffread1a(IDIM,field,field,nvals,         &
                       ilabel,icode,cmessage)
END IF
!=======================================================================
9999 CONTINUE

IF (icode  /=  0) CALL ereport("FFREADB", icode, cmessage)
RETURN
END SUBROUTINE ffreadb
