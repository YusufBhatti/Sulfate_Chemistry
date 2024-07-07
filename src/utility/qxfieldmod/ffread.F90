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
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE ffread(iproj,fct,itype,int_level,ppunit,field,IDIM,    &
ilabel,rlabel,iextra,icode,cmessage)
USE io
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER ::                                                        &
     maxff                                                        &
                            !OUT Max number of opened files
    ,len1_lookup                                                  &
                            !IN  first dimension of the lookup
    ,pp_len2_lookup                                               &
                            !OUT secnd dimension of the lookup
    ,len1_record                                                  &
                            !OUT First dimension of record
    ,reclen                                                       &
                            !OUT Total length of record.
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
                            !INOUT Controls certain functions.
    ,icode                                                        &
                            !OUT return code
    ,maxpp                                                        &
                            !    maximum number of unit number
    ,data_add                                                     &
                            !OUT The word address of the data.
    ,pfno                   !OUT No of fields files opened
PARAMETER(maxff=10)
PARAMETER(maxpp=100)
PARAMETER(len1_lookup=64)
PARAMETER(len1_record=30000)!Max size of a lookup table allowed
PARAMETER(reclen=len1_record*maxff) ! Total length of RECORD
REAL ::                                                           &
     field(IDIM)                                                  &
                            !OUT array holding final output data.
    ,rlabel(19)                                                   &
                            !OUT holds real part of LOOKUP
    ,real_level             !IN  LEVEL code (could be real)

! -------------------------------------------------------------------
!     LOCAL VARIABLES
INTEGER ::                                                        &
     table(maxpp)                                                 &
                            ! associates unit no and file no
    ,prev_ppunit(maxff)                                           &
                            ! a record of unit nos already used
    ,i                                                            &
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
    ,len_fixhd                                                    &
                            ! Length of Fixed length header
    ,pp_fixhd(256)                                                &
                            ! Fixed length header
    ,in_lbvc                ! Local copy of LBVC required to searc
REAL ::                                                           &
     a_io                   ! OUTPUT from UNIT command
LOGICAL ::                                                        &
     test                                                         &
    ,record(len1_record,maxff)
SAVE prev_ppunit
SAVE pfno
SAVE table
SAVE record
!
PARAMETER(len_fixhd=256)
DATA prev_ppunit/maxff*0/
DATA table/maxpp*0/
DATA pfno/0/
DATA record/reclen*.FALSE./
!    Remember that BUFFER OUT starts at address 0
IF (ppunit >  maxpp) THEN
  icode=1
  cmessage=' FFREAD   the unit number is too large'
  RETURN
END IF
test=.TRUE.
!  Establish if a completely new FF is being read.
DO i=1,maxff
  IF (ppunit == prev_ppunit(i)) THEN
    test=.FALSE.
  END IF
END DO
!   A TABLE is set up associating FIELDS_FILE NO (1 to MAXFF) with a
!   PP unit number. On succesive calls to the FF this table is used
!   help record which LOOKUP table belongs to which FF (PPUNIT)
IF (test) THEN   ! A FF never read in before.
  pfno=pfno+1
  prev_ppunit(pfno)=ppunit
  table(ppunit)=pfno
END IF
!   Read in the word address of the LOOKUP table (IWA), the length of
!   the LOOKUP table (PP_LEN2_LOOKUP) and the start address of the data
!   DATA_ADD
pfno=table(ppunit)
iwa=0
CALL setpos(ppunit,iwa,icode) ! C coded routine
CALL buffin(ppunit,pp_fixhd,len_fixhd,len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= len_fixhd) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('Buffer in fixed length header',a_io,len_io,       &
                len_fixhd)
  cmessage='  FFREAD : I/O error reading FIXED LENGTH HEADER'
  icode=2
  RETURN
END IF
iwa= pp_fixhd(150)-1  ! NOTE for BUFFER IN the start address
!                             ! is zero for word 1
data_add= pp_fixhd(160)  ! The start address of the data
pp_len2_lookup=pp_fixhd(152)
! DEPENDS ON: ffreada
CALL ffreada(iproj,fct,itype,int_level,ppunit,field,IDIM,         &
ilabel,rlabel,iextra,test,pp_len2_lookup,len1_lookup,pp_fixhd,    &
iwa,len1_record,maxff,record,pfno,data_add,icode,cmessage)

IF (icode  /=  0) CALL ereport("FFREAD", icode, cmessage)

RETURN
END SUBROUTINE ffread
