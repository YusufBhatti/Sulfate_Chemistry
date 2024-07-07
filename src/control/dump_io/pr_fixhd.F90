! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_FIXHD---------------------------------------
!
!    Purpose: Prints out fixed length header record and checks
!             validity of information.
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!             Code Owner: Please refer to the UM file CodeOwners.txt           
!             This file belongs in section: Dump I/O

SUBROUTINE pr_fixhd                                               &
(fixhd,len_fixhd,len_inthd,len_realhd,len1_levdepc                &
,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc &
,len1_flddepc,len2_flddepc,len_extcnst,len_dumphist,len_cfi1      &
,len_cfi2,len_cfi3,len1_lookup,len2_lookup,len_data               &
,icode,cmessage)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE umPrintMgr
USE missing_data_mod, ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER ::                                                        &
 len_fixhd                                                        &
               !IN Length of fixed length header
,len_inthd                                                        &
               !IN Length of integer header
,len_realhd                                                       &
               !IN Length of real header
,len1_levdepc                                                     &
               !IN 1st dim of level dep consts
,len2_levdepc                                                     &
               !IN 2nd dim of level dep consts
,len1_rowdepc                                                     &
               !IN 1st dim of row dep consts
,len2_rowdepc                                                     &
               !IN 2nd dim of row dep consts
,len1_coldepc                                                     &
               !IN 1st dim of column dep consts
,len2_coldepc                                                     &
               !IN 2nd dim of column dep consts
,len1_flddepc                                                     &
               !IN 1st dim of field dep consts
,len2_flddepc                                                     &
               !IN 2nd dim of field dep consts
,len_extcnst                                                      &
               !IN Length of extra constants
,len_dumphist                                                     &
               !IN Length of history block
,len_cfi1                                                         &
               !IN Length of comp field index 1
,len_cfi2                                                         &
               !IN Length of comp field index 2
,len_cfi3                                                         &
               !IN Length of comp field index 3
,len1_lookup                                                      &
               !IN 1st dim of lookup
,len2_lookup   !IN 2nd dim of lookup

INTEGER ::                                                        &
 fixhd(len_fixhd)                                                 &
                  !IN Fixed length header
,len_data                                                         &
                  !IN Length of real data
,icode          !OUT Return code; successful=0
                !                 error > 0

CHARACTER(LEN=errormessagelength) ::                              &
 cmessage       !OUT Error message if ICODE > 0

! -------------------------------------------------------------
! Local variables:---------------------------------------------
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PR_FIXHD'
!--------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=' '

IF ( printstatus >= prstatus_oper ) THEN
  CALL umPrint('',src='pr_fixhd')
  WRITE(umMessage,'('' FIXED LENGTH HEADER'')')
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' -------------------'')')
  CALL umPrint(umMessage,src='pr_fixhd')

  WRITE(umMessage,'('' Dump format version'',I6)')fixhd(1)
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' UM Version No      '',I6)')fixhd(12)
  CALL umPrint(umMessage,src='pr_fixhd')
END IF

IF (fixhd(2) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Atmospheric data'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(2) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Oceanic data'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(2) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Wave sub-model data'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE
  icode=10
  WRITE(cmessage, '(A,I0,A)')                                                 &
    'PR_FIXHD: Consistency check' // newline //                               &
    'Invalid sub-model indicator: Fixed Length Header(2) = ',fixhd(2),        &
     newline //                                                               &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (fixhd(3) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' On hybrid levels'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(3) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' On sigma levels'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(3) == 3) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' On pressure levels'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(3) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' On depth levels'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(3) == 5) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Charney-Phillips on radius levels'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(3) == 6) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    CALL umPrint(                                                &
    ' Wave model direction levels and frequency pseudo-levels')
  END IF
ELSE IF (fixhd(3) == imdi .AND. fixhd(5) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    CALL umPrint(                                                &
    ' Missing data indicator used for vert coord type')
  END IF
ELSE
  icode=20
  WRITE(cmessage, '(A,I0,A)')                                                &
    'PR_FIXHD: Consistency check' // newline //                              &
    'Invalid vertical coordinate type: Fixed Length Header(3) = ',fixhd(3),  &
     newline //                                                              &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (fixhd(4) == 0) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over global domain'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(4) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over N. Hemispheric domain'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(4) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over S. Hemispheric domain'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(4) == 3) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over LAM domain with no wrap around'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(4) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over LAM domain with wrap around'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(4) == 103) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Over rotated LAM domain'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE
  icode=30
  WRITE(cmessage, '(A,I0,A)')                                               &
    'PR_FIXHD: Consistency check' // newline //                             &
    'Invalid horizontal grid type: Fixed Length Header(4) = ',fixhd(4),     &
     newline //                                                             &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (fixhd(5) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Instantaneous dump'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Meaned dump'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 3) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' FIELDS file'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Ancillary dataset'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 5) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Boundary dataset'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 6) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' AC Observation File'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 7) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Var Observation File'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 8) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Cx file (model columns at ob locations)'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 9) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Covariance File'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(5) == 10) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage, '('' OPS Obstore file '')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE
  icode=40
  WRITE(cmessage, '(A,I0,A)')                                             &
    'PR_FIXHD: Consistency check' // newline //                           &
    'Invalid dataset type: Fixed Length Header(5) = ',fixhd(5),           &
     newline //                                                           &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF ( printstatus >= prstatus_oper ) THEN
  WRITE(umMessage,'('' Exp No ='',I6,'' Run Id ='',I6)') fixhd(7),fixhd(6)
  CALL umPrint(umMessage,src='pr_fixhd')
END IF

IF (fixhd(8) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Gregorian calendar'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(8) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' 360-day calendar'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(8) == imdi .AND. fixhd(5) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Missing data indcator used as calendar'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE
  icode=50
  WRITE(cmessage, '(A,I0,A)')                                                 &
    'PR_FIXHD: Consistency check' // newline //                               &
    'Invalid calendar type: Fixed Length Header(8) = ',fixhd(8),              &
     newline //                                                               &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (fixhd(9) == 1) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Arakawa A grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == 2) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Arakawa B grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == 3) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Arakawa C grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Arakawa D grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == 5) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Arakawa E grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == 6) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' V-AT-POLES/ENDGAME grid'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == imdi .AND. fixhd(5) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Missing data indicator used for grid type'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE IF (fixhd(9) == imdi .AND. fixhd(5) == 5) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(umMessage,'('' Missing data indicator used for grid type'')')
    CALL umPrint(umMessage,src='pr_fixhd')
  END IF
ELSE
  icode=60
  WRITE(cmessage, '(A,I0,A)')                                                 &
    'PR_FIXHD: Consistency check' // newline //                               &
    'Invalid grid staggering indicator: Fixed Length Header(9) = ',fixhd(9),  &
     newline //                                                               &
    'Please see UMDP F03 for details of appropriate values.'
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF ( printstatus >= prstatus_oper ) THEN

  IF ( fixhd(5) == 5 ) THEN    !  Boundary dataset
    WRITE(umMessage, '(A21,7A7)') 'Record:              ',               &
         'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec', 'DayNo'
    CALL umPrint(umMessage, src='pr_fixhd')
    WRITE(umMessage,'(''First Validity time ='',7I7)')(fixhd(i),i=21,27)
    CALL umPrint(umMessage,src='pr_fixhd')
    WRITE(umMessage,'(''Last  Validity time ='',7I7)')(fixhd(i),i=28,34)
    CALL umPrint(umMessage,src='pr_fixhd')
    WRITE(umMessage,'(''Interval            ='',7I7)')(fixhd(i),i=35,41)
    CALL umPrint(umMessage,src='pr_fixhd')

  ELSE
    WRITE(umMessage, '(A15,7A7)') 'Record:        ',                      &
         'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec', 'DayNo'
    CALL umPrint(umMessage, src='pr_fixhd')
    WRITE(umMessage,'(''Data time     ='',7I7)')(fixhd(i),i=21,27)
    CALL umPrint(umMessage,src='pr_fixhd')
    WRITE(umMessage,'(''Validity time ='',7I7)')(fixhd(i),i=28,34)
    CALL umPrint(umMessage,src='pr_fixhd')
    WRITE(umMessage,'(''Creation time ='',7I7)')(fixhd(i),i=35,41)
    CALL umPrint(umMessage,src='pr_fixhd')

  END IF

  CALL umPrint(&
  '                        Start     1st dim    2nd dim   1st parm    2nd parm'&
  ,src='pr_fixhd')
  WRITE(umMessage,'('' Integer Consts   '',2I11,11X,I11)')fixhd(100),     &
  fixhd(101),len_inthd
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Real Consts      '',2I11,11X,I11)')fixhd(105),     &
  fixhd(106),len_realhd
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Level Dep Consts '',5I11)')fixhd(110),             &
  fixhd(111),fixhd(112),len1_levdepc,len2_levdepc
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Row Dep Consts   '',5I11)')fixhd(115),             &
  fixhd(116),fixhd(117),len1_rowdepc,len2_rowdepc
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Column Dep Consts'',5I11)')fixhd(120),             &
  fixhd(121),fixhd(122),len1_coldepc,len2_coldepc
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Fields of Consts '',5I11)')fixhd(125),             &
  fixhd(126),fixhd(127),len1_flddepc,len2_flddepc
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Extra Consts     '',2I11,11X,I11)')fixhd(130),     &
  fixhd(131),len_extcnst
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' History Block    '',2I11,11X,I11)')fixhd(135),     &
  fixhd(136),len_dumphist
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' CFI No 1         '',2I11,11X,I11)')fixhd(140),     &
  fixhd(141),len_cfi1
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' CFI No 2         '',2I11,11X,I11)')fixhd(142),     &
  fixhd(143),len_cfi2
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' CFI No 3         '',2I11,11X,I11)')fixhd(144),     &
  fixhd(145),len_cfi3
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Lookup Tables    '',5I11)')fixhd(150),             &
  fixhd(151),fixhd(152),len1_lookup,len2_lookup
  CALL umPrint(umMessage,src='pr_fixhd')
  WRITE(umMessage,'('' Model Data       '',2I11,11X,I11)')fixhd(160),     &
  fixhd(161),len_data
  CALL umPrint(umMessage,src='pr_fixhd')
END IF

! Check model parameters against header record entries

IF (fixhd(101) >  0) THEN
  IF (len_inthd /= fixhd(101)) THEN
    icode=70
    WRITE(cmessage, '(A,I0,A,I0)')                                            &
    'PR_FIXHD: Consistency check' // newline //                               &
    'Length of integer constants is: ',len_inthd,                             &
     newline //                                                               &
    'This does not match size specified in Fixed Length Header (101): ',      &
    fixhd(101)
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF (fixhd(105) >  0) THEN
  IF (len_realhd /= fixhd(106)) THEN
    icode=80
    WRITE(cmessage, '(A,I0,A,I0)')                                            &
    'PR_FIXHD: Consistency check' // newline //                               &
    'Length of real constants is: ',len_realhd,                               &
     newline //                                                               &
    'This does not match size specified in Fixed Length Header (106): ',      &
    fixhd(106)
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF (fixhd(110) >  0) THEN
  IF (len1_levdepc /= 0) THEN
    IF (len1_levdepc /= fixhd(111) ) THEN
      icode=100
      WRITE(cmessage, '(A,I0,A,I0)')                                         &
      'PR_FIXHD: Consistency check' // newline //                            &
      'Length of first dimension of level dependent constants is: ',         &
       len1_levdepc,  newline //                                             &
      'This does not match size specified in Fixed Length Header (111): ',   &
       fixhd(111)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    IF (len2_levdepc /= fixhd(112) ) THEN
      icode=110
      WRITE(cmessage, '(A,I0,A,I0)')                                          &
      'PR_FIXHD: Consistency check' // newline //                             &
      'Length of second dimension of level dependent constants is: ',         &
       len2_levdepc,  newline //                                              &
      'This does not match size specified in Fixed Length Header (112): ',    &
      fixhd(112)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(115) >  0) THEN
  IF (len1_rowdepc /= 0) THEN
    IF (len1_rowdepc /= fixhd(116) ) THEN
      icode=120
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of first dimension of row dependent constants is: ',             &
       len1_rowdepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (116): ',     &
       fixhd(116)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    IF (len2_rowdepc /= fixhd(117) ) THEN
      icode=130
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of second dimension of row dependent constants is: ',            &
       len2_rowdepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (117): ',     &
      fixhd(117)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(120) >  0) THEN
  IF (len1_coldepc /= 0) THEN
   IF (len1_coldepc /= fixhd(121) ) THEN
      icode=140
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of first dimension of col dependent constants is: ',             &
       len1_coldepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (121): ',     &
       fixhd(121)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    IF (len2_coldepc /= fixhd(122) ) THEN
      icode=150
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of second dimension of col dependent constants is: ',            &
       len2_coldepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (122): ',     &
      fixhd(122)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(125) >  0) THEN
  IF (len1_flddepc /= 0) THEN
   IF (len1_flddepc /= fixhd(126) ) THEN
      icode=160
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of first dimension of fields of constants is: ',                 &
       len1_flddepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (126): ',     &
       fixhd(126)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    IF (len2_flddepc /= fixhd(127) ) THEN
      icode=170
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of second dimension of fields of constants is: ',                &
       len2_flddepc,  newline //                                               &
      'This does not match size specified in Fixed Length Header (127): ',     &
      fixhd(127)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(130) >  0) THEN
  IF (len_extcnst /= 0) THEN
    IF (len_extcnst /= fixhd(131)) THEN
      icode=180
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of extra constants is: ',                                        &
       len_extcnst,  newline //                                                &
      'This does not match size specified in Fixed Length Header (131): ',     &
      fixhd(131)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(135) >  0) THEN
  IF (len_dumphist /= 0) THEN
    IF (len_dumphist /= fixhd(136)) THEN
      icode=190
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of history file is: ',                                           &
       len_dumphist,  newline //                                               &
      'This does not match size specified in Fixed Length Header (136): ',     &
      fixhd(136)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(140) >  0) THEN
  IF (len_cfi1 /= 0) THEN
    IF (len_cfi1 /= fixhd(141)) THEN
      icode=200
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of compressed field index 1: ',                                  &
       len_cfi1,  newline //                                                   &
      'This does not match size specified in Fixed Length Header (141): ',     &
      fixhd(141)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(142) >  0) THEN
  IF (len_cfi2 /= 0) THEN
    IF (len_cfi2 /= fixhd(143)) THEN
      icode=210
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of compressed field index 2: ',                                  &
       len_cfi2,  newline //                                                   &
      'This does not match size specified in Fixed Length Header (143): ',     &
      fixhd(143)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(144) >  0) THEN
  IF (len_cfi3 /= 0) THEN
    IF (len_cfi3 /= fixhd(145) ) THEN
      icode=220
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of compressed field index 3: ',                                  &
       len_cfi2,  newline //                                                   &
      'This does not match size specified in Fixed Length Header (145): ',     &
      fixhd(145)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(150) >  0) THEN
  IF (len2_lookup /= imdi) THEN
    IF (len1_lookup/=fixhd(151)) THEN
      icode=230
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of first dimension of lookup table: ',                           &
       len1_lookup,  newline //                                                &
      'This does not match size specified in Fixed Length Header (151): ',     &
      fixhd(151)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
   END IF
    IF (len2_lookup/=fixhd(152)) THEN
      icode=240
      WRITE(cmessage, '(A,I0,A,I0)')                                           &
      'PR_FIXHD: Consistency check' // newline //                              &
      'Length of second dimension of lookup table: ',                          &
       len2_lookup,  newline //                                                &
      'This does not match size specified in Fixed Length Header (152): ',     &
      fixhd(152)
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_fixhd
