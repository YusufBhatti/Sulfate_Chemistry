! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE FLD_OUT------------------------------------------
!
!    REPLACES THE OUTPUT FROM STASH EITHER ON TO A PP FILE OR
!    BACK TO THE MAIN ARRAY D1
!
!    PURPOSE:   TO PROCESS DIAGNOSTICS CONTROLLED BY STASH
!
!
!

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs


!
!    ARGUMENTS:---------------------------------------------------
SUBROUTINE fld_out                                                &
          (icode,cmessage,bufout,lenbuf,len_buf_words,num_words,  &
           unitpp,len1_lookup,pp_len2_lookup,ipplook,rpplook,     &
           ilabel,rlabel,iwl,data_addr)
USE io
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


CHARACTER(LEN=errormessagelength) :: cmessage !OUT OUT MESSAGE FROM ROUTINE
!

INTEGER ::                                                        &
  icode                                                           &
                     !IN    RETURN CODE FROM ROUTINE
, len1_lookup                                                     &
                     !IN    FIRST DIMENSION OF LOOKUP TABLE
, pp_len2_lookup                                                  &
                     !IN    SECND DIMENSION OF LOOKUP TABLE
, lenbuf                                                          &
                     !IN     LENGTH OFF PP BUFFER
, unitpp                                                          &
                     !IN     OUTPUT PP UNIT NUMBER
, len_buf_words                                                   &
                     !IN
, num_words          !IN
!
INTEGER ::                                                        &
  jj            !IN    ITEM NUMBER
INTEGER ::                                                        &
  ipplook(len1_lookup,pp_len2_lookup)                             &
                                      !IN INTEGER LOOKUP TABLE
, ilabel(45)                                                      &
                ! INTEGER PART OF LOOKUP
, iwl                                                             &
                !IN    Address of the PP LOOKUP Table
, data_addr     !IN    Address of start of data
!
REAL ::                                                           &
  bufout(lenbuf)                                                  &
                         !OUTPUT PP BUFFER (ROUNDED UP)
, rpplook(len1_lookup,pp_len2_lookup)                             &
                                      !IN REAL LOOKUP TABLE
, rlabel(19)    ! REAL PART OF LOOKUP

! ---------------------------------------------------------------------

!    WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
! ---------------------------------------------------------------------
!     EQUIVALENCE(IPPLOOK,RPPLOOK)
!
! ------------------------------------------------------------------
!   MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
! ---------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
INTEGER ::                                                        &
  addr                                                            &
                !
, iwa                                                             &
                !     RECORD NUMBER
, ix                                                              &
                !     RETURN VALUE FROM UNIT COMMAND
, len_io                                                          &
                !
, ii                                                              &
                !     COUNTER
, i                                                               &
                !     COUNTER
,ierr           ! Error return from SETPOS

REAL ::                                                           &
  a_io          !

INTEGER ::                                                        &
  lresid                                                          &
                !
, icurrll                                                         &
                !
, ipast                                                           &
                !
, iproj         !     M08 PROJECTION NUMBER

LOGICAL ::                                                        &
  first              !
DATA first/.TRUE./

LOGICAL :: found
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
first=.TRUE.
icurrll=0

103 FORMAT(//,32x,' ARRAY FROM START OF PPOUT  ',//,32(10f8.0/))
lresid=len_buf_words-num_words
DO jj=num_words+1,lresid
  bufout(jj)= 0.0
END DO
!
found=.TRUE.
IF (first) THEN
  found=.FALSE.
  DO  jj=1,pp_len2_lookup
    IF (ipplook(1,jj) <  0) THEN  ! Search for last entry
      icurrll=jj
      IF (jj == 1) THEN
        iwa=((iwl+511)/512)*512+pp_len2_lookup*len1_lookup
        WRITE(umMessage,*) 'Start data',iwa,data_addr,iwa-1
        CALL umPrint(umMessage,src='fld_out')
        iwa=data_addr
        iwa=iwa-1
      ELSE
        iwa= ipplook(29,jj-1)+ipplook(30,jj-1) !ADDR+LGTH
      END IF
      found=.TRUE.
      EXIT
    END IF
  END DO
  IF (.NOT. found) THEN
    icode=1
    cmessage="FROM PPOUT CANNOT FIND SUITABLE ENTRY IN LOOKUP"
    ! FIXME : Should there be an ereport here?
  END IF
ELSE
  ipast=icurrll-1
  WRITE(7,105) ipast
  105    FORMAT('  FROM PPOUT AND FIRST IS FALSE IPAST=',i8)
  iwa=ipplook(29,ipast) + ipplook(30,ipast) ! ADDR + LENGTH
  WRITE(7,106) iwa
  106    FORMAT('  FROM PPOUT AND FIRST IS FALSE IWA=',i8)
END IF

IF (found) THEN
  !
  !     update lookup for this field
  DO i=1,45
    ipplook(i,icurrll) = ilabel(i)
  END DO
  DO i=1,19
    rpplook(i+45,icurrll) = rlabel(i)
  END DO
  ipplook(29,icurrll)=iwa
  ipplook(30,icurrll)=len_buf_words
  ipplook(40,icurrll)=iwa
  !
  CALL setpos(unitpp,iwa,ierr)
  CALL buffout(unitpp,bufout,LEN_buf_words,len_io,a_io)
  100     FORMAT(//,32x,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10f8.0/))
  101     FORMAT(//,32x,'   IPPLOOK AT END OF  PPOUT   ',//,32(16i5/))
  102     FORMAT('     IWA  LEN_BUF_WORDS ',2i12)
END IF
RETURN
END SUBROUTINE fld_out
