! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: FLDMOD --------------------------------------------------
!
!    Purpose: To read a   direct access PP file  and convert it to a
!    sequential file read to be passed across to the IBM
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
PROGRAM fldmod
USE filenamelength_mod, ONLY:                                    &
    filenamelength
USE io
USE ereport_mod, ONLY: ereport
USE UM_Config, ONLY: &
    appInit,          &
    appTerminate,     &
    exe_fldmod
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE fort2c_interfaces, ONLY: get_file
USE errormessagelength_mod, ONLY: errormessagelength
USE hostname_mod,           ONLY: get_hostname

IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=filenamelength) :: diagfile
CHARACTER(LEN=filenamelength) :: infile
CHARACTER(LEN=8) :: c_nproc            ! to get nproc_x and nproc_y from
!                                    ! environment variables.
!                         up to an EVEN no for conversion to IBM format
INTEGER ::                                                        &
     len_fixhd                                                    &
                            !    Length of fixed length header
    ,len_inthd                                                    &
    ,len_realhd                                                   &
    ,len1_levdpc                                                  &
    ,len2_levdpc                                                  &
    ,len1_rowdpc                                                  &
    ,len2_rowdpc                                                  &
    ,len1_coldpc                                                  &
    ,len2_coldpc                                                  &
    ,len1_lookup                                                  &
    ,len2_lookup                                                  &
    ,ppunit1                                                      &
                            !OUT unit no of required fieldsfile 1
    ,ppunit2                                                      &
                            !OUT unit no of required fieldsfile 2
    ,icode                                                        &
                            !IN  return code
    ,err
!    LOCAL VARIABLES
PARAMETER(len_fixhd=256)
INTEGER ::                                                        &
     i                                                            &
                            ! local counter
    ,pp_fixhd(len_fixhd)                                          &
                            !IN  Fixed length header
    ,iwa                                                          &
                            !
    ,ix                                                           &
                            !
    ,len_io                 !
REAL ::                                                           &
     a_io                   !
!

CHARACTER(LEN=*),PARAMETER :: RoutineName = 'fldmod'
INTEGER, PARAMETER :: diag_unit = 7

CALL appInit(exe_fldmod)

! Write out hostname
WRITE(umMessage,'(A,A)') 'Host is ',TRIM(get_hostname())
CALL umPrint(umMessage, src='fldmod')

!    OPEN DIAGNOSTIC FILE
CALL get_file(diag_unit,diagfile,filenamelength,icode)
OPEN(UNIT=diag_unit,FILE=diagfile)

! ****************************************************************
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
icode = 0
cmessage= '                                         '
!
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
! ****************************************************************
ppunit1=10
ppunit2=11
! ****************************************************************
!     Buffer in the Fixed Length Header and obtain lengths
! ****************************************************************
CALL get_file(ppunit1,infile,filenamelength,icode)
CALL file_open(ppunit1,infile,filenamelength,0,1,icode)
CALL buffin(ppunit1,pp_fixhd,len_fixhd,len_io,a_io)
IF (a_io /= -1.0 .OR. len_io /= len_fixhd) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('Buffer in fixed length header',a_io,len_io,     &
                len_fixhd)
  cmessage='FFREAD : I/O error reading FIXED LENGTH HEADER'
  icode=2
  WRITE(umMessage,*)'umthin1 - I/O error reading FIXED LENGTH HEADER'
  CALL umPrint(umMessage,src='fldmod')
  CALL ereport(RoutineName,icode,cmessage)
END IF
len_inthd=pp_fixhd(101)
len_realhd=pp_fixhd(106)
len1_levdpc=pp_fixhd(111)
len2_levdpc=pp_fixhd(112)
len1_rowdpc=pp_fixhd(116)
len2_rowdpc=pp_fixhd(117)
len1_coldpc=pp_fixhd(121)
len2_coldpc=pp_fixhd(122)
len1_lookup=pp_fixhd(151)
len2_lookup=pp_fixhd(152)
! DEPENDS ON: dimens1
CALL dimens1(len_inthd,len_realhd,len1_levdpc,len2_levdpc,        &
   len1_rowdpc,len2_rowdpc,len1_coldpc,len2_coldpc,               &
   len1_lookup,len2_lookup,len_fixhd,pp_fixhd,ppunit1,ppunit2,    &
   icode,cmessage)
IF (icode /= 0) THEN
  WRITE(diag_unit,100) icode
  WRITE(diag_unit,110) cmessage
END IF

100  FORMAT(' ICODE EQUAL TO ',i2)
110  FORMAT(a80)

CALL appTerminate()
END PROGRAM fldmod
