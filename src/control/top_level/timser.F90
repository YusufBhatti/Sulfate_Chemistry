! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Process time-series domain data (if any)
! Subroutine interface:
SUBROUTINE timser(cmessage,ErrorStatus)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE stextend_mod, ONLY: npos_ts, nrecs_ts
USE cstash_mod, ONLY: elim_ts, wlim_ts, nlim_ts, slim_ts, i1_ts,  &
                      i51_ts, tlimr_ts, rlevlst_d, ig_ts, imn_d,  &
                      ndprof, blim_ts, iopl_d, tlim_ts, levlst_d, &
                      blimr_ts, levb_d, levt_d

USE stparam_mod, ONLY: &
  st_levels_model_rho, st_levels_model_theta, st_levels_deep_soil,  &
  st_levels_single, st_levels_cloud_thresh

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE missing_data_mod, ONLY: imdi

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!

! Subroutine arguments
!   Scalar arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage      ! Error return message

! Error status:
INTEGER ::     ErrorStatus ! Error return code

! Local variables:
INTEGER :: BlkId           !Time series block identifier
INTEGER :: BlkSt           !Start position of ts block data
INTEGER :: Nrecs_prev      !No of recs in previous time ser block
INTEGER :: idp             !Domain profile loop counter
INTEGER :: ipos            !Position in ts limits arrays
INTEGER :: isblim,istlim   !Used for converting vertical ts
INTEGER :: ib,it,il,ilvl   !  domain limits to sequence nos.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TIMSER'

!- End of Header ------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!Loop over domain profiles
BlkSt =1
DO idp=1,ndprof
  IF (npos_ts(idp) >  0) THEN
    !  This domain profile has a time series
    !    Identify TS block using pointer array
    BlkId = npos_ts (idp)
    !    Find start position (in LIM_TS arrays) of data for this block
    IF (BlkId >  1) THEN
      BlkSt=BlkSt+Nrecs_prev
    END IF
    !  Loop over records in ts block corresponding to domain profile IDP.
    !  Adjust the TS records for domain profiles with vertical or horiz
    !    averaging.
    !  Convert the ts domain vertical limits to sequence nos.
    !    in the domain profile levels range/levels list.
    DO ipos=BlkSt,BlkSt+nrecs_ts(npos_ts(idp))-1
      !    Vertical levels
      IF (iopl_d(idp) == st_levels_model_rho   .OR.  &
          iopl_d(idp) == st_levels_model_theta .OR.  &
          iopl_d(idp) == st_levels_deep_soil) THEN
        !           Model levels
        IF (imn_d(idp) == 1) THEN
          !             Vertical mean
          blim_ts(ipos)=1
          tlim_ts(ipos)=1
        ELSE
          !             No vertical mean
          IF (levb_d(idp) >= 0) THEN
            !               Range of model levels
            IF (blim_ts(ipos) <  levb_d(idp) .OR.                  &
              tlim_ts(ipos) >  levt_d(idp)) THEN
              WRITE(umMessage,'(A,A,A,I6,A,I6)') 'ERROR, TIMSER: ',      &
            ' TS_DOMAIN LEVEL LIMIT OUT OF RANGE; ',                     &
            ' DOM PROF: ',idp,                                           &
            ' TS RECORD: ',ipos
              CALL umPrint(umMessage,src='timser')
              ErrorStatus=1
              cmessage='TS DOMAIN LEVEL LIMIT OUT OF RANGE'
              GO TO 9999
            END IF
            blim_ts(ipos)=blim_ts(ipos)-levb_d(idp)+1
            tlim_ts(ipos)=tlim_ts(ipos)-levb_d(idp)+1
          ELSE
            !               List of selected model levels;
            !               LEVT_D(IDP)=no. of levels in list
            isblim=imdi
            istlim=imdi
            DO il=1,levt_d(idp)
              IF (blim_ts(ipos) == levlst_d(il,idp)) isblim=il
              IF (tlim_ts(ipos) == levlst_d(il,idp)) istlim=il
            END DO
            IF ((istlim == imdi) .OR.                               &
               (isblim == imdi)) THEN
              WRITE(umMessage,'(A,A,I6,A,I6)')                           &
             'ERROR TIMSER:T-SERIES INTEGER LEVEL NOT IN ',              &
             'LEVELS LIST; DOM PROF: ',idp,' TS RECORD: ',ipos
              CALL umPrint(umMessage,src='timser')
              WRITE(umMessage,'(A,2I6)') 'SPECIFIED TS LEVELS LIMITS: ', &
              blim_ts(ipos),tlim_ts(ipos)
              CALL umPrint(umMessage,src='timser')
              ErrorStatus = 1
              cmessage=                                           &
             'ERROR TIMSER:T-SERIES LEVEL NOT IN LEVELS LIST'
              GO TO 9999
            END IF
            !                 Store seq. nos. of ts domain level limits
            blim_ts(ipos)=isblim
            tlim_ts(ipos)=istlim
          END IF
        END IF
        !           List of specified real levels
      ELSE IF ((iopl_d(idp) /= st_levels_single) .AND.     &
               (iopl_d(idp) <= st_levels_cloud_thresh)) THEN
        IF (imn_d(idp) == 1) THEN
          blim_ts(ipos)=1
          tlim_ts(ipos)=1
        ELSE
          !             Determine sequence nos. of top & bottom ts domain
          !             levels in real levels list (ISBLIM, ISTLIM), by
          !             representing real level values as integers.
          isblim=imdi
          istlim=imdi
          ib=(blimr_ts(ipos)*1000.0+0.5)
          it=(tlimr_ts(ipos)*1000.0+0.5)
          DO il=1,levt_d(idp)
            ilvl=(rlevlst_d(il,idp)*1000.0+0.5)
            IF (ib == ilvl) isblim=il
            IF (it == ilvl) istlim=il
          END DO
          IF ((istlim == imdi) .OR.                                 &
             (isblim == imdi)) THEN
            WRITE(umMessage,'(A,A,I6,A,I6)')                              &
           'ERROR TIMSER:T-SERIES REAL LEVEL NOT IN ',                   &
           'LEVELS LIST; DOM PROF: ',idp,' TS RECORD: ',ipos
            CALL umPrint(umMessage,src='timser')
            WRITE(umMessage,'(A,2I6)') 'SPECIFIED TS LEVELS LIMITS: ',    &
            blimr_ts(ipos),tlimr_ts(ipos)
            CALL umPrint(umMessage,src='timser')
            ErrorStatus = 1
            cmessage=                                             &
           'ERROR TIMSER:T-SERIES LEVEL NOT IN LEVELS LIST'
          END IF
          !               Store seq. nos. of ts domain level limits
          blim_ts(ipos)=isblim
          tlim_ts(ipos)=istlim
        END IF
      ELSE IF (iopl_d(idp) == st_levels_single) THEN
        !           Single level
        blim_ts(ipos)=1
        tlim_ts(ipos)=1
      ELSE
        WRITE(umMessage,'(A,I6)')                                        &
            'ERROR TIMSER: UNEXPECTED LEVEL TYPE CODE',iopl_d(idp)
        CALL umPrint(umMessage,src='timser')
        ErrorStatus=1
        GO TO 9999
      END IF
      !    Horizontal area
      IF (imn_d(idp) == 2) THEN
        elim_ts(ipos)=1
        wlim_ts(ipos)=1
      ELSE IF (imn_d(idp) == 3) THEN
        nlim_ts(ipos)=1
        slim_ts(ipos)=1
      ELSE IF (imn_d(idp) == 4) THEN
        elim_ts(ipos)=1
        wlim_ts(ipos)=1
        nlim_ts(ipos)=1
        slim_ts(ipos)=1
      END IF
      ig_ts =0  ! These constants are left-overs from the
      i1_ts =1  !  pre-vn3.5 TIMSER routine: they are used
      i51_ts=51 !  in the UM time-series routines.
    END DO      ! IPOS loop
    Nrecs_prev=nrecs_ts(npos_ts(idp)) ! For next TS block
  END IF        ! TS(IDP) == 'Y'
END DO          ! IDP loop

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE timser
