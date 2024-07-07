! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE DIAGOPR ------------------------------------------------
!
!    Purpose : Provide Statistics on Precipitation Rates
!
!              Calculate and print out the following :
!              - Contingency tables of Observation vs Forecast data
!              - Various skill scores
!
!

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE diagopr_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGOPR_MOD'

CONTAINS

SUBROUTINE diagopr  (lenobt_dim,                                             &
    intpr,obdata,LENOBT_local,lenob,ndv,                                     &
    lmissd,nobipt,nptobt,no_anal_levs)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc, nproc_max
USE umPrintMgr
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod

IMPLICIT NONE

INTEGER :: lenobt_dim
INTEGER :: LENOBT_local,lenob,ndv,nobipt,nptobt,no_anal_levs
REAL :: intpr  (lenobt_dim),obdata (lenobt_dim,ndv)
LOGICAL :: lmissd (lenob+1,no_anal_levs)

!   INTENT=INOUT-----------------------------------------------
!   INTPR        - interpolated precipitation rates
!   OBDATA       - observation data
!   LENOBT       - number of observations of this type
!   NDV          - number of data values
!   NOBIPT       - pointer to data in OBDATA
!   LMISSD       - logical array to point to obs with missing data
!   LENOB        - number of obs in group of ob types
!   NO_ANAL_LEVS - no of analysis levels
!   NPTOBT       - pointer to present ob type within group
!
!   INTENT=OUT-------------------------------------------------
!  ------------------------------------------------------------
!     Local arrays and variables
INTEGER :: lenobt
INTEGER :: iproc

INTEGER :: obcat,fcat,job,jcat,jcat1,                                         &
    missing,obs_used,                                                      &
    jbdy,counter,counterf,                                                 &
    counterf2
REAL :: avg_fcrate,avg_obrate,                                             &
    obs_pcnt,fc_pcnt,crct_pcnt,                                            &
    rms,rmsf,cc,                                                           &
    rmsf2,                                                                 &
    hr(3),far(3),bias_rate,bias_area,                                      &
    hks(3),ts(3),                                                          &
    tab_rate(5,5),tab_yn(2,2)
INTEGER ::                                                                 &
    missing_rem(0:nproc_max),                                              &
    counter_rem(0:nproc_max),counterf_rem(0:nproc_max),                    &
    counterf2_rem(0:nproc_max),                                            &
    lenobt_rem(0:nproc_max)
INTEGER :: istat  ! for calls to GCOM
REAL :: avg_fcrate_rem(0:nproc_max),avg_obrate_rem(0:nproc_max),           &
    RMS_rem(nproc_max),RMSF_rem(0:nproc_max),                              &
    rmsf2_rem(0:nproc_max),                                                &
    tab_rate_rem(5,5,0:nproc_max)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAGOPR'


!
!   TAB_RATE   - precipitation magnitude contingency table
!   TAB_YN     - 2x2 'yes/no' table compiled from TAB_RATE
!   OBCAT,FCAT - pointers within TAB_RATE, TAB_PHASE
!   OBS_USED   - total number of observations used
!   MISSING    - counter for no of obs with missing data
!   AVG_OBRATE - mean observed precipitation rate
!   AVG_FCRATE - mean forecast precipitation rate at obs points
!   OBS_PCNT   - percentage of observations which give precipitation
!   FC_PCNT    - percentage of forecast values which give precipitation
!   CRCT_PCNT  - percentage of correct forecasts
!   RMS        - Root mean square error of forecast rate v obs
!   RMSF       - Root mean square factor
!   RMSF2      - Root mean square factor, for either value > threshold
!   CC         - Tetrachoric Correlation Coefficient
!   HR         - Hit Rate               }  arrays
!   FAR        - False Alarm Rate       }  with values based
!   HKS        - Hanssen & Kuiper Score }  on light, moderate
!   TS         - Threat Score           }  and heavy thresholds
!   BIAS_RATE  - Bias Ratio for precipitation rates
!   BIAS_AREA  - Bias Ratio for number of points with rain
!   JOB        - loop counter over obs
!   JCAT    }  - loop counters
!   JCAT1   }  -               over categories of obs/fc precip
!   JBDY       - loop counter for precip rate category boundaries
!   COUNTER    - counter for obs contributing to rms rate error
!   COUNTERF   - counter for obs contributing to rms factor score
!   COUNTERF2  - counter for obs contributing to RMSF2 score
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0) THEN
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  PRECIPITATION verification routine DIAGOPR")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Called from VANRAIN")')
  CALL umPrint(umMessage,src='diagopr')
  IF (l_verif_range) THEN
    WRITE(umMessage,'(A,F12.6,A)')                                         &
        'Verification range is ',radar_range,' km from nearest radar'
    CALL umPrint(umMessage,src='diagopr')
  ELSE
    WRITE(umMessage,'(A,F12.6,A)')                                         &
        'Verification range is ',radar_range_max,' km from nearest radar'
    CALL umPrint(umMessage,src='diagopr')
  END IF
  CALL umPrint(' ',src='diagopr')
END IF

IF (ltimer_ac) CALL timer('DIAGOPR ',3)

!  ** 1.0  Initialize variables
DO jcat=1,5
  DO jcat1=1,5
    tab_rate(jcat,jcat1) = 0
  END DO
END DO
lenobt=0
missing=0
avg_fcrate=0.0
avg_obrate=0.0
counter=0
counterf=0
counterf2=0
rmsf2=0.0
rms=0.0
rmsf=0.0
cc =9999999.9

!
!  **2.0   Calculate average forecast and observed rates,
!          rms error in fc rate and rms factor.
!          Ignore points where ob is missing data.
!
DO job=1,LENOBT_local
  IF (lmissd(nptobt+job,1)) THEN
    missing = missing + 1
  ELSE
    avg_fcrate = avg_fcrate + intpr(job)
    avg_obrate = avg_obrate + obdata(job,nobipt)
    IF ( intpr(job)  >   thresh_dl .AND.                                   &
        obdata(job,nobipt)  >   thresh_dl ) THEN
      counter = counter + 1
      rms = rms + ( intpr(job) - obdata(job,nobipt) ) ** 2
    END IF
    IF ( intpr(job)  >   thresh_rmsf .AND.                                 &
        obdata(job,nobipt)  >   thresh_rmsf ) THEN
      counterf = counterf + 1
      rmsf=rmsf+(LOG(ABS( intpr(job)/obdata(job,nobipt) )))**2
    END IF
    IF ( intpr(job)  >= thresh_rmsf .OR.                                   &
        obdata(job,nobipt)  >= thresh_rmsf ) THEN
      counterf2 = counterf2 + 1
      rmsf2=rmsf2 + (LOG( MAX(thresh_rmsf/2.0,intpr(job))                  &
          / MAX(thresh_rmsf/2.0,obdata(job,nobipt)) ))**2

    END IF
  END IF
END DO ! Job

!
!  ** 3.0   calculate 4x4 table for rain magnitudes
!
DO job = 1,LENOBT_local
  IF ( .NOT. lmissd(nptobt+job,1) ) THEN
    fcat = 1
    IF ( intpr(job)  >   thresh_dl ) fcat=2
    IF ( intpr(job)  >   thresh_lm ) fcat=3
    IF ( intpr(job)  >   thresh_mh ) fcat=4
    obcat = 1
    IF ( obdata(job,nobipt)  >   thresh_dl ) obcat=2
    IF ( obdata(job,nobipt)  >   thresh_lm ) obcat=3
    IF ( obdata(job,nobipt)  >   thresh_mh ) obcat=4
    tab_rate(obcat,fcat) = tab_rate(obcat,fcat) + 1
  END IF
END DO ! Job

IF (mype == 0) THEN
  lenobt_rem(mype)=lenobt_local
  missing_rem(mype)=missing
  counter_rem(mype)=counter
  counterf_rem(mype)=counterf
  counterf2_rem(mype)=counterf2
  rms_rem(mype)=rms
  rmsf_rem(mype)=rmsf
  rmsf2_rem(mype)=rmsf2
  avg_fcrate_rem(mype)=avg_fcrate
  avg_obrate_rem(mype)=avg_obrate
  jbdy=0
  DO jcat=1,5
    DO jcat1=1,5
      jbdy=jbdy+1
      tab_rate_rem(1,1,jbdy)=tab_rate(jcat1,jcat)
    END DO
  END DO
END IF

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_irecv(iproc,1,iproc,istat,lenobt_rem(iproc),                   &
        lenobt_local)
  ELSE IF (mype == iproc) THEN
    CALL gc_isend(iproc,1,0,istat,lenobt_rem(iproc),                       &
        lenobt_local)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_irecv(iproc,1,iproc,istat,missing_rem(iproc),                  &
        missing)
  ELSE IF (mype == iproc) THEN
    CALL gc_isend(iproc,1,0,istat,missing_rem(iproc),                      &
        missing)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_irecv(iproc,1,iproc,istat,counter_rem(iproc),                  &
        counter)
  ELSE IF (mype == iproc) THEN
    CALL gc_isend(iproc,1,0,istat,counter_rem(iproc),                      &
        counter)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_irecv(iproc,1,iproc,istat,counterf_rem(iproc),                 &
        counterf)
  ELSE IF (mype == iproc) THEN
    CALL gc_isend(iproc,1,0,istat,counterf_rem(iproc),                     &
        counterf)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_irecv(iproc,1,iproc,istat,counterf2_rem(iproc),                &
        counterf2)
  ELSE IF (mype == iproc) THEN
    CALL gc_isend(iproc,1,0,istat,counterf2_rem(iproc),                    &
        counterf2)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,1,iproc,istat,rms_rem(iproc),rms)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,1,0,istat,rms_rem(iproc),rms)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,1,iproc,istat,rmsf_rem(iproc),rmsf)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,1,0,istat,rmsf_rem(iproc),rmsf)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,1,iproc,istat,rmsf2_rem(iproc),rmsf2)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,1,0,istat,rmsf2_rem(iproc),rmsf2)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,1,iproc,istat,avg_fcrate_rem(iproc),               &
        avg_fcrate)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,1,0,istat,avg_fcrate_rem(iproc),                   &
        avg_fcrate)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,1,iproc,istat,avg_obrate_rem(iproc),               &
        avg_obrate)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,1,0,istat,avg_obrate_rem(iproc),                   &
        avg_obrate)
  END IF
END DO

DO iproc=1,nproc-1
  IF (mype == 0) THEN
    CALL gc_rrecv(iproc,25,iproc,istat,tab_rate_rem(1,1,iproc),            &
        tab_rate)
  ELSE IF (mype == iproc) THEN
    CALL gc_rsend(iproc,25,0,istat,tab_rate_rem,                           &
        tab_rate)
  END IF
END DO

IF (mype == 0) THEN

  lenobt=lenobt_local

  DO iproc=1,nproc-1

    lenobt=lenobt+lenobt_rem(iproc)
    missing=missing+missing_rem(iproc)
    counter=counter+counter_rem(iproc)
    counterf=counterf+counterf_rem(iproc)
    counterf2=counterf2+counterf2_rem(iproc)
    rms=rms+rms_rem(iproc)
    rmsf=rmsf+rmsf_rem(iproc)
    rmsf2=rmsf2+rmsf2_rem(iproc)
    avg_fcrate=avg_fcrate+avg_fcrate_rem(iproc)
    avg_obrate=avg_obrate+avg_obrate_rem(iproc)

    DO jcat=1,4
      DO jcat1=1,4
        tab_rate(jcat,jcat1)= tab_rate(jcat,jcat1)                         &
            +tab_rate_rem(jcat,jcat1,iproc)
      END DO
    END DO

  END DO !iproc=0,nproc-1

  obs_used = lenobt - missing
  IF (obs_used == 0) THEN
    CALL umPrint('  *** all obs are missing data ***',src='diagopr')
    GO TO 9999
  END IF

  ! calculate average rates, rms rate error and rms factor
  ! convert from model units of  kg m-2 s-1  to  mm/hr
  !
  avg_fcrate = avg_fcrate / obs_used * 3600.0
  avg_obrate = avg_obrate / obs_used * 3600.0

  IF (counter  ==  0) THEN
    rms = 99999999.0
  ELSE
    rms = SQRT(rms/counter)
  END IF
  IF (counterf  ==  0) THEN
    rmsf= 99999999.0
  ELSE
    rmsf = EXP( SQRT(rmsf/counterf) )
  END IF
  IF (counterf2  ==  0) THEN
    rmsf2 = 99999999.0
  ELSE
    rmsf2 = EXP( SQRT(rmsf2/counterf2) )
  END IF


  !
  !   Convert tables to percentages, and calculate row totals
  !
  DO jcat=1,4
    DO jcat1=1,4
      tab_rate(jcat,jcat1) = tab_rate(jcat,jcat1)*100.0/obs_used
      tab_rate(5,jcat1) = tab_rate(5,jcat1) + tab_rate(jcat,jcat1)
      tab_rate(jcat,5) = tab_rate(jcat,5) + tab_rate(jcat,jcat1)
    END DO
  END DO

  !
  !  ** 4.0  Condense 4x4 TAB_RATE into various 2x2 tables, TAB_YN
  !
  DO jbdy=1,3
    !
    !  Calculate TAB_YN for each boundary
    !
    IF (jbdy  ==  1) THEN
      !   dry / light  boundary
      tab_yn(1,1) = tab_rate(1,1)
      tab_yn(2,1) = tab_rate(2,1) + tab_rate(3,1) + tab_rate(4,1)
      tab_yn(1,2) = tab_rate(1,2) + tab_rate(1,3) + tab_rate(1,4)
      tab_yn(2,2) = tab_rate(2,2) + tab_rate(3,2) + tab_rate(4,2)          &
          + tab_rate(2,3) + tab_rate(3,3) + tab_rate(4,3)                  &
          + tab_rate(2,4) + tab_rate(3,4) + tab_rate(4,4)
    ELSE IF (jbdy  ==  2) THEN
      !   light / moderate  boundary
      tab_yn(1,1) = tab_rate(1,1) + tab_rate(2,1)                          &
          + tab_rate(1,2) + tab_rate(2,2)
      tab_yn(2,1) = tab_rate(3,1) + tab_rate(3,2)                          &
          + tab_rate(4,1) + tab_rate(4,2)
      tab_yn(1,2) = tab_rate(1,3) + tab_rate(2,3)                          &
          + tab_rate(1,4) + tab_rate(2,4)
      tab_yn(2,2) = tab_rate(3,3) + tab_rate(3,4)                          &
          + tab_rate(4,3) + tab_rate(4,4)
    ELSE
      !   moderate / heavy  boundary
      tab_yn(2,2) = tab_rate(4,4)
      tab_yn(2,1) = tab_rate(4,1) + tab_rate(4,2) + tab_rate(4,3)
      tab_yn(1,2) = tab_rate(1,4) + tab_rate(2,4) + tab_rate(3,4)
      tab_yn(1,1) = tab_rate(1,1) + tab_rate(1,2) + tab_rate(1,3)          &
          + tab_rate(2,1) + tab_rate(2,2) + tab_rate(2,3)                  &
          + tab_rate(3,1) + tab_rate(3,2) + tab_rate(3,3)
    END IF
    !
    !  ** 4.1  Now calculate statistics for each set of TAB_YN
    !

    !  Hit Rate : percentage of observed precipitation correctly forecast
    IF ( (tab_yn(2,1)+tab_yn(2,2))  >   0) THEN
      hr(jbdy) = tab_yn(2,2) * 100.0                                       &
          / ( tab_yn(2,1) + tab_yn(2,2) )
    ELSE
      hr(jbdy) = 999999.9
    END IF

    !  False Alarm Rate : percentage of forecast precipitation observed
    !                     to be dry
    IF ( (tab_yn(1,2)+tab_yn(2,2))  >   0) THEN
      far(jbdy)= tab_yn(1,2) * 100.0                                       &
          / ( tab_yn(1,2) + tab_yn(2,2) )
    ELSE
      far(jbdy)= 9999999.9
    END IF

    !  Hanssen & Kuiper : Hit Rate - percentage of dry observations
    !                              incorrectly forecast
    IF ( (tab_yn(1,1)+tab_yn(1,2))  >   0) THEN
      hks(jbdy) = hr(jbdy) - tab_yn(1,2) * 100.0                           &
          / ( tab_yn(1,1) + tab_yn(1,2) )
    ELSE
      hks(jbdy) = 999999.9
    END IF

    !  Threat Score : Correctly forecast precipitation
    !                    / ( total observed precipitation
    !                      + total forecast precipitation
    !                      - correct forecast precipitation )
    !  threat score is also known as CSI - critical success index
    IF ( (tab_yn(1,2)+tab_yn(2,1)+tab_yn(2,2))  >   0.0) THEN
      ts(jbdy) = tab_yn(2,2) * 100.0                                       &
          / ( tab_yn(1,2) + tab_yn(2,1) + tab_yn(2,2) )
    ELSE
      ts(jbdy) = 999999.9
    END IF

    !   Tetrachoric correlation coefficient
    IF ( (jbdy  ==  1) .AND. ( (tab_yn(2,1)+tab_yn(2,2))  >   0.0)         &
        .AND.  ( (tab_yn(1,2)+tab_yn(2,2))  >   0.0 )                      &
        .AND.  ( (tab_yn(1,1)+tab_yn(1,2))  >   0.0 )                      &
        .AND.  ( (tab_yn(1,1)+tab_yn(2,1))  >   0.0 ) ) THEN
      cc = ( tab_yn(1,1)*tab_yn(2,2) - tab_yn(1,2)*tab_yn(2,1) ) /         &
          SQRT( (tab_yn(2,1)+tab_yn(2,2))*(tab_yn(1,2)+tab_yn(2,2))        &
          * (tab_yn(1,1)+tab_yn(1,2)) * (tab_yn(1,1)+tab_yn(2,1)) )
    END IF
  END DO ! jbdy

  !
  !  ** 5.0 Calculate overall scores
  !
  !  percentage of all obs which give precipitation
  obs_pcnt = 100.0 - ( tab_rate(1,1) + tab_rate(1,2) +                     &
      tab_rate(1,3) + tab_rate(1,4) )

  !  percentage of all forecasts which give precipitation
  fc_pcnt = 100.0 - ( tab_rate(1,1) + tab_rate(2,1) +                      &
      tab_rate(3,1) + tab_rate(4,1) )

  !  percentage of forecasts in correct category
  crct_pcnt = tab_rate(1,1) + tab_rate(2,2) +                              &
      tab_rate(3,3) + tab_rate(4,4)

  !  Bias Ratio : forecast precipitation / observed precipitation
  !
  !  BIAS_RATE : by mean rate
  IF (avg_obrate  >   0.0) THEN
    bias_rate = avg_fcrate / avg_obrate
  ELSE
    bias_rate = 99999.9
  END IF
  !  BIAS_AREA : by area
  IF (obs_pcnt  >   0.0) THEN
    bias_area = fc_pcnt / obs_pcnt
  ELSE
    bias_area = 99999.9
  END IF

  !  ** 6.0 Print results

  CALL umPrint(' ',src='diagopr')
  WRITE(umMessage,'("  * table entries are % of obs")')
  CALL umPrint(umMessage,src='diagopr')
  CALL umPrint(' ',src='diagopr')
  WRITE(umMessage,'("  Precipitation Rate Table")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  ------------------------")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(21X,"OBS")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(16X," Dry    Light  Moderate  Heavy   TOTAL")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(13X,"------------------------------------------")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,                                                         &
      '("  F.C.    Dry| ",F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)')          &
      tab_rate(1,1),tab_rate(2,1),tab_rate(3,1),                           &
      tab_rate(4,1),tab_rate(5,1)
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,                                                         &
      '("        Light| ",F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)')          &
      tab_rate(1,2),tab_rate(2,2),tab_rate(3,2),                           &
      tab_rate(4,2),tab_rate(5,2)
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,                                                         &
      '("     Moderate| ",F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)')          &
      tab_rate(1,3),tab_rate(2,3),tab_rate(3,3),                           &
      tab_rate(4,3),tab_rate(5,3)
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,                                                         &
      '("        Heavy| ",F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)')          &
      tab_rate(1,4),tab_rate(2,4),tab_rate(3,4),                           &
      tab_rate(4,4),tab_rate(5,4)
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,                                                         &
      '("       TOTAL | ",F5.1,3X,F5.1,3X,F5.1,3X,F5.1,4X,F5.1)')          &
      tab_rate(1,5),tab_rate(2,5),tab_rate(3,5),                           &
      tab_rate(4,5),100.0
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  STATISTICS")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  ----------")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("    threshold used   :      d/l     l/m    m/h")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Hit Rate           HR= ",3F7.1)') hr
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  False Alarm Rate  FAR= ",3F7.1)') far
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Hanssen & Kuiper  HKS= ",3F7.1)') hks
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Threat Score       TS= ",3F7.1)') ts
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("    (above scores are %)")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  % of obs with precip   ",F5.1)') obs_pcnt
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  % fc with precip       ",F5.1)') fc_pcnt
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Mean Forecast Rate     ",F6.3,"  mm/hour")')        &
      avg_fcrate
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Mean Observation Rate  ",F6.3,"  mm/hour")')        &
      avg_obrate
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Bias Ratio (rate)      = ",F6.3)') bias_rate
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("             (area)      = ",F6.3)') bias_area
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  % of fc in correct cat = ",F5.1)') crct_pcnt
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  RMS rate error         = ",F7.3,"  mm/hr")')        &
      rms*3600
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  RMS factor (intersection)  RMSF= ",F7.3)') rmsf
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("    (calculated from threshold ",F6.3," mm/hr)")')    &
      thresh_rmsf*3600
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("     No. of obs passing threshold = ",I6)') counterf
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  RMS factor (union)         RMSF2= ",F7.3)') rmsf2
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("     No. of obs for RMSF2    = ",I6)') counterf2
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Correlation Coeff  CC= ",F7.3)') cc
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("    (based on rain/no rain)")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Number of Obs used     ",I6)') obs_used
  CALL umPrint(umMessage,src='diagopr')
  IF (l_latlon_prver) THEN
    WRITE(umMessage,'("  Verification in area :")')
    CALL umPrint(umMessage,src='diagopr')
    WRITE(umMessage,'("    lat ",F6.3," deg","  to  ",F6.3," deg")')       &
        southlat,northlat
    CALL umPrint(umMessage,src='diagopr')
    WRITE(umMessage,'("    lon ",F6.3," deg","  to  ",F6.3," deg")')       &
        westlon,eastlon
    CALL umPrint(umMessage,src='diagopr')
  END IF
  WRITE(umMessage,'(" ")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("  Threshold values used :")')
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("     Dry / Light         ",F6.3,"  mm/hour")')        &
      thresh_dl*3600
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("     Light / Moderate    ",F6.3,"  mm/hour")')        &
      thresh_lm*3600
  CALL umPrint(umMessage,src='diagopr')
  WRITE(umMessage,'("     Moderate / Heavy    ",F6.3,"  mm/hour")')        &
      thresh_mh*3600
  CALL umPrint(umMessage,src='diagopr')
  9999  CONTINUE

END IF ! mype == 0


IF (ltimer_ac) CALL timer('DIAGOPR ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagopr
END MODULE diagopr_mod
