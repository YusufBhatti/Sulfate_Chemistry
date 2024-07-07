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
SUBROUTINE ffreada(iproj,fct,itype,int_level,ppunit,field,IDIM,   &
ilabel,rlabel,iextra,test,pp_len2_lookup,len1_lookup,pp_fixhd,    &
iwa,len1_record,maxff,record,pfno,data_add,icode,cmessage)
USE io
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER ::                                                        &
     maxff                                                        &
                            !IN  Max number of opened files
    ,len1_lookup                                                  &
                            !IN  first dimension of the lookup
    ,len1_record                                                  &
                            !INOUT First dimension of record
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
    ,maxpp                                                        &
                            !    maximum number of unit number
    ,data_add                                                     &
                            !IN  The word address of the data.
    ,pfno                   !INOUT No of fields files opened
INTEGER ::                                                        &
     pp_fixhd(*),                                                 &
                                        !IN PPfile fixed header
     lookup(len1_lookup,pp_len2_lookup) !OUTinteger lookup
REAL ::                                                           &
     field(IDIM)                                                  &
                            !OUT array holding final output data.
    ,rlabel(19)                                                   &
                            !OUT holds real part of LOOKUP
    ,real_level             !IN  LEVEL code (could be real)
LOGICAL ::                                                        &
     record(len1_record,maxff) !INOUT Record of the field no read
!     LOCAL VARIABLES
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
    ,icount                                                       &
                            ! Counter
    ,len_io                                                       &
                            ! Length of data transferred from BUF
    ,len_io_expected                                              &
                            ! Length od data expected in transfer
    ,length_of_data                                               &
                            ! Length of a particular field
    ,addr                                                         &
                            ! Address of a field in the data store
    ,in_lbvc                ! Local copy of LBVC required to searc
REAL ::                                                           &
     a_io                   ! OUTPUT from UNIT command
LOGICAL ::                                                        &
     test
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
!
!     Read in the LOOKUP table.
!
CALL setpos(ppunit,iwa,icode) ! C coded routine
len_io_expected=pp_len2_lookup*len1_lookup
CALL buffin(ppunit,lookup,len_io_expected,len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= len_io_expected) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('Buffer in lookup table   ',a_io,len_io,         &
  len_io_expected)
  cmessage='FFREADA: I/O error reading LOOKUP TABLE  '
  icode=3
  RETURN
END IF
IF (iextra(2) == 0) THEN  ! Allows duplicate entries to be read
  IF (pp_len2_lookup >  len1_record) THEN
    cmessage='FFREADA: LEN1_RECORD NOT LARGE ENOUGH    '
    icode=4
    RETURN
  END IF
  DO i=1,len1_record
    IF (record(i,pfno)) THEN
      lookup(14,i)=-99
    END IF
  END DO
END IF
! DEPENDS ON: ffreadb
CALL ffreadb      (iproj,fct,itype,int_level,ppunit,field,IDIM,   &
ilabel,rlabel,iextra,pp_len2_lookup,len1_lookup,                  &
iwa,len1_record,maxff,record,pfno,pp_fixhd,lookup,lookup,data_add,&
icode,cmessage)


IF (icode  /=  0) CALL ereport("FFREADA", icode, cmessage)

RETURN
END SUBROUTINE ffreada
