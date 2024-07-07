! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE DIAGO --------------------------------------------------
!
!    Purpose : Provide Statistics on 'Observation-Model' Increments
!
!              Calculate and print out the following :
!              Mean Increment
!              Mean Square Increment
!              Extreme Increment and position (Lat/Long) of obs
!
!              Can be used on Individual or group of AC Obs Types
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE diago_mod
USE umPrintMgr
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGO_MOD'

CONTAINS

SUBROUTINE diago (csubr,ktype,kcall,                              &
           obs_incr,normf,obs_lat,obs_long,                       &
           lmissd,lenob,lenobt,nptobt,                            &
           no_anal_levs,no_anal_var)

USE conversions_mod, ONLY: recip_pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE comobs_mod, ONLY: nobtypmx
USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon, vector_eq_to_latlon, &
                                  eq_latlon_vector_coeffs
USE ac_control_mod

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

CHARACTER(LEN=*) :: csubr   !  Subroutine where DIAGO is called from
INTEGER :: ktype,kcall,lenob,lenobt,nptobt,no_anal_levs
INTEGER :: no_anal_var
REAL ::                                                           &
 obs_incr(lenob+1,no_anal_levs,no_anal_var),                      &
 normf (lenob+1,no_anal_levs), obs_lat (lenobt+1),                &
 obs_long (lenobt+1)
LOGICAL :: lmissd(lenob+1,no_anal_levs)

!     INTENT=IN---------------------------------------------------------
!     KTYPE    : AC Observation Type
!     KCALL    : Indicator to where DIAGO called from
!              : 1 = FROM AC
!                2 = FROM VAN## ; IPASS=1 ; BEFORE VERTICAL FILTERING
!                3 = FROM VAN## ; IPASS=1 ; AFTER  VERTICAL FILTERING
!                4 = FROM VAN## ; IPASS=2 ; BEFORE VERTICAL FILTERING
!                5 = FROM VAN## ; IPASS=2 ; AFTER  VERTICAL FILTERING
!                6 = from DIAGOCC ; IPASS=1 ; for multi-level cloud
!                7 = from DIAGOCC ; IPASS=1 ; for low cloud
!                8 = from DIAGOCC ; IPASS=1 ; for medium cloud
!                9 = from DIAGOCC ; IPASS=1 ; for high cloud
!               10 = from DIAGOCC ; IPASS=1 ; for total cloud
!                VAN## is Vertical Analysis routine calling DIAGO.
!     INTENT=INOUT------------------------------------------------------
!     INC      : Increments - P*, T, U, RH
!     VINC     : Increments - V
!     NORMF    : Normalisation factor
!     OBS_LAT  : Observation co-latitudes
!     OBS_LONG : Observation longitudes
!     LMISSD   : Array to control which obs are used on each level
!     LENOB    : No of obs in group of observation types
!     LENOBT   : No of obs for Observation Type KTYPE
!     NPTOBT   : Offset to first observation for type KTYPE
!     NO_ANAL_LEVS : No of analysis levels
!     INTENT=OUT--------------------------------------------------------
!    ------------------------------------------------------------------
REAL :: a_maxinc(model_levels_max)
INTEGER :: iproc,istat

!     Workspace Usage:--------------------------------------------------
!     -----------------------------------------------------------
!     Local work arrays
!     -----------------------------------------------------------
! *   Dynamic allocation
LOGICAL :: obsuse(lenobt+1)
REAL :: wlt(lenobt+1), wln(lenobt+1), wln2(lenobt+1)
REAL :: coeff1(lenobt+1), coeff2(lenobt+1)
REAL :: uwk(lenobt+1,no_anal_levs), vwk(lenobt+1,no_anal_levs)
!     WLT/WLN  : Real lat/long of wind obs
!     WLN2     : ELF  longitude of wind obs
!     COEFF1/2 : Coefficients to rotate wind components
!     UWK/VWK  : Rotated u/v components on real lat/long grid
!     ------------------------------------------------------------------

!     Local arrays and variables

INTEGER :: jvar,jband,jlev,job,jobt,jnanl
INTEGER :: nband,inobs,inobslv,iobt,anal_var
REAL :: rstat(0:no_anal_levs,8),lat,longitude,ui,vi,nf,sqinc
REAL :: bandlat(4),maxinc,sumui,sumvi,sumnf,sumsqi
CHARACTER :: ns*1,we*1,title1*9,title2*25,lab1*10,lab2*10
CHARACTER(LEN=3) :: bandc(4)
LOGICAL :: lwind
LOGICAL :: l_mdi = .TRUE.  ! check for MDI

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAGO'


!     ANAL_VAR : Analysis varibale (=KTYPE/100)
!     BANDLAT  : Boundaries of latitude bands
!     BANDC    : Labels to print out boundaries of latitude bands
!     INOBS    : Number of observations in latitude band
!     INOBSLV  : Number of observations on level
!     IOBT     : Index to list of obs types in MDIAGTPS
!     JBAND    : Loop variable for latitude bands
!     JLEV     : Loop variable for analysis levels
!     JNANL    : Number of analysis levels to be processed
!     JOB      : Loop variable for observations
!     JOBT     : Loop variable for observation types
!     JVAR     : Loop variable for statistic variable
!     LAB1/2   : Labels for print out
!     LAT      : Absolute value of latitude of maximum increment
!     LONG     : Absolute value of longitude of maximum increment
!     LWIND    : Indicator for working with wind obs type
!     MAXINC   : Maximum increment for level
!     NBAND    : Number of latitude bands
!     NF       : Normalisation Factors
!     NS/WE    : Labels for printing out lat/long of maximum increment
!     SQINC    : Square of increments
!     STAT     : Array storing statistics data
!                STAT(0,#) Overall
!                STAT(JLEV,#) For level JLEV
!                STAT(#,1) No of observations on  level JLEV
!                STAT(#,2) Mean Norm factor
!                STAT(#,3) Mean Increment - p*, t, u, rh
!                STAT(#,4) Mean Increment - v
!                STAT(#,5) Mean square Increment
!                STAT(#,6) Extreme increment
!                STAT(#,7) Lat of extreme increment
!                STAT(#,8) Long of extreme increment
!     SUMNF    : Sum of norm. factors
!     SUMUI    : Sum of increments - p*, t, u, rh
!     SUMVI    : Sum of increments - v
!     SUMSQI   : Sum of squared increments
!     TITLE1/2 : Titles for print out
!     UI       : Increments - p*, t, u, rh
!     VI       : Increments - v

!     ------------------------------------------------------------------
DATA bandlat / 0.0, 1.0472, 2.0944, 3.1416 /
DATA bandc/'90N','30N','30S','90S'/
!     -----------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('DIAGO   ',3)

!     Check that 'Obs-Model' statistics required
IF (.NOT. lldac(2)) GO TO 9999

!     For verification purposes, only get statistics on first pass.
IF (lverif .AND. (kcall /= 2 .AND. kcall /= 3)) GO TO 9999

!     Check that statistics required for calling routine
IF (mdiag >  0 .AND. kcall /= mdiag) GO TO 9999

!     Check that statistics required for observation type KTYPE
IF (mdiagtps(1) >  0) THEN
  iobt=0
  DO jobt=1,nobtypmx
    IF (ktype == mdiagtps(jobt)) iobt=jobt
  END DO
  IF (iobt == 0) GO TO 9999
END IF

!---- Set up titles

lab1='    '
lab2='    '
anal_var = ktype/100
IF (anal_var == 1) THEN
  lab1 = 'P* OBS INC'
ELSE IF (anal_var == 2) THEN
  IF (ltemp) THEN
    lab1 = ' T OBS INC'
  ELSE
    lab1 = 'TH OBS INC'
  END IF
ELSE IF (anal_var == 3) THEN
  lab1 = ' U OBS INC'
  lab2 = ' V OBS INC'
ELSE IF (anal_var == 4) THEN
  IF ( kcall >  5 ) THEN
    lab1 = 'CC OBS INC'
  ELSE
    lab1 = 'RH OBS INC'
  END IF
ELSE IF (anal_var == 9) THEN
  lab1 = 'LOG(VIS) OBS INC'
ELSE
  GO TO 9999
END IF

IF (kcall == 2 .OR. kcall == 3) THEN
  title1 = 'IPASS=1'
ELSE IF (kcall == 4 .OR. kcall == 5) THEN
  title1 = 'IPASS=2'
ELSE
  title1 = '  '
END IF

title2 = ' '
IF (ktype == 201 .OR. ktype == 205 .OR. ktype == 206 .OR.         &
    ktype == 207 .OR. ktype == 208 .OR. ktype == 209 .OR.         &
    ktype == 301 .OR. ktype == 401 ) THEN
  IF (kcall == 2 .OR. kcall == 4) THEN
    title2 = 'BEFORE VERTICAL FILTERING'
  ELSE IF (kcall == 3 .OR. kcall == 5) THEN
    title2 = 'AFTER VERTICAL FILTERING'
  ELSE
    title2 = ' '
  END IF
END IF

lwind = anal_var == 3

jnanl = no_anal_levs

IF (ktype == 202 .OR. ktype == 204 .OR.                           &
    ktype == 302 .OR. ktype == 304 .OR.                           &
    ktype == 305 .OR. ktype == 306 .OR.                           &
    ktype == 402 .OR. ktype == 404                                &
    .OR. ktype  ==  901) THEN
  jnanl=1
END IF

!     ONLY ONE BAND OF STATISTICS ALLOWED FOR ELF GRID
nband=1

!     Loop over latitude bands
DO jband=1,nband

  !     Initialise statistics data
  DO jvar=1,8
    DO jlev=0,no_anal_levs
      rstat(jlev,jvar)=0.0
    END DO
  END DO

  !     1.0 Determine which observations to be used
  IF (lenobt <= 0) THEN
    DO jlev=1,jnanl
      DO jvar=0,8
        s_stat(jlev,jvar)=0.0
      END DO
    END DO
  ELSE

    !     1.1 Initialise so that all observations to be used
    DO job=1,lenobt
      obsuse(job)=.TRUE.
    END DO

    IF (model_type == mt_global) THEN
      !     1.2 Skip observations outside latitude band
      IF (llband) THEN
        DO job=1,lenobt
          obsuse(job) = obs_lat(job) >= bandlat(jband)   .AND.          &
                        obs_lat(job) <  bandlat(jband+1)
        END DO
      END IF
    END IF

    !     1.3 Count number of observations to be used
    inobs = 0
    DO job=1,lenobt
      IF (obsuse(job)) inobs = inobs+1
    END DO

    IF (model_type /= mt_global) THEN

      !       2. Co-ordinate tranformation  (ELF version only)
      !       ------------------------------------------------

      !       Get real Lat/Long of ELF observations

      !       Lat/Lon values in OBS_LAT/OBS_LONG are for ELF grid.
      !       First convert from radians to degrees.
      DO job=1,lenobt
        wlt (job) = 90.0-obs_lat(job)*recip_pi_over_180
        wln (job) = obs_long(job)*recip_pi_over_180
        wln2(job) = wln(job)
      END DO
      !       Get real Lat/Lon of wind obs in WLT/WLN, keep ELF lon in WLN2
      CALL rotate_eq_to_latlon(wlt,wln,wlt,wln,elfplat,elfplon,lenobt)

      IF (lwind) THEN

        !         Rotate wind components on ELF grid to real lat/lon system.

        !         Calculate coefficients of rotation.
        CALL eq_latlon_vector_coeffs(coeff1, coeff2, wln, wln2, &
                                     elfplat, elfplon, lenobt)
        !         Rotation of wind components in loop over levels.
        DO jlev=1,jnanl
          CALL vector_eq_to_latlon(coeff1, coeff2,                    &
                                   obs_incr(nptobt+1,jlev,1),         &
                                   obs_incr(nptobt+1,jlev,2),         &
                                   uwk(1,jlev), vwk(1,jlev), lenobt, l_mdi)
        END DO

      END IF
    END IF


    !       3.0 Calculate statistics required
    !       ---------------------------------

    !-----  Loop over levels
    DO jlev=1,jnanl

      inobslv=0
      sumnf  = 0.0
      sumui  = 0.0
      sumvi  = 0.0
      sumsqi = 0.0
      maxinc = 0.0

      IF (model_type == mt_global) THEN
        !-----  Loop over observations
        DO job=1,lenobt

          nf    = 0.0
          ui    = 0.0
          vi    = 0.0
          sqinc = 0.0

          IF (obsuse(job) .AND. .NOT. lmissd(nptobt+job,jlev)) THEN

            inobslv = inobslv+1

            IF (lnormf) THEN
              !           Use NF as weighting
              nf = normf (nptobt+job,jlev)
            ELSE
              !           Do not use NF as weighting (set to 1.0)
              nf = 1.0
            END IF
            ui = obs_incr(nptobt+job,jlev,1)
            IF (lwind) vi = obs_incr(nptobt+job,jlev,2)


            !-----    3.1 Sum of Normalisation factors used as weights

            !         For IPASS=1, preliminary norm factor (Q1M) calculated in
            !                      VAN### is used to weight the obs increments.
            !         For IPASS=2, final norm factor (Q2M) calculated in
            !                      VAN### is used to weight the obs increments.

            sumnf = sumnf + nf

            !-----    3.2 Weighted sum of increments - P*, T, U, RH
            sumui = sumui + nf*ui

            !-----    3.3 Weighted sum of increments - V
            IF (lwind) sumvi = sumvi + nf*vi

            !-----    3.4 Weighted sum of squared increments
            sqinc = ui*ui
            IF (lwind) sqinc = sqinc + vi*vi
            sumsqi = sumsqi + nf*sqinc

            !-----    3.5 Determine maximum (squared) increment
            IF (sqinc  >   maxinc) THEN
              maxinc = sqinc
              IF (lwind) THEN
                !         For wind, extreme is maximum speed of vector increment
                rstat(jlev,6) = SQRT(sqinc)
              ELSE
                !         Extreme is increment
                rstat(jlev,6) = ui
              END IF
              !           Lat/Long of obs
              rstat(jlev,7) = obs_lat(job)
              rstat(jlev,8) = obs_long(job)
              !-----      Convert to degrees
              rstat(jlev,7) = 90.0 - rstat(jlev,7) * recip_pi_over_180
              rstat(jlev,8) = rstat(jlev,8) * recip_pi_over_180
              IF (rstat(jlev,8) >  180.0) rstat(jlev,8)=rstat(jlev,8)-360.0
            END IF

          END IF

        END DO !JOB
      ELSE
        !-----  Loop over observations
        DO job=1,lenobt

          nf    = 0.0
          ui    = 0.0
          vi    = 0.0
          sqinc = 0.0

          IF (obsuse(job) .AND. .NOT. lmissd(nptobt+job,jlev)) THEN

            inobslv = inobslv+1

            IF (lnormf) THEN
              !           Use NF as weighting
              nf = normf (nptobt+job,jlev)
            ELSE
              !           Do not use NF as weighting (set to 1.0)
              nf = 1.0
            END IF
            IF (lwind) THEN
              ui = uwk(job,jlev)
              vi = vwk(job,jlev)
            ELSE
              ui = obs_incr(nptobt+job,jlev,1)
            END IF

            !-----    3.1 Sum of Normalisation factors used as weights

            !         For IPASS=1, preliminary norm factor (Q1M) calculated in
            !                      VAN### is used to weight the obs increments.
            !         For IPASS=2, final norm factor (Q2M) calculated in
            !                      VAN### is used to weight the obs increments.

            sumnf = sumnf + nf

            !-----    3.2 Weighted sum of increments - P*, T, U, RH
            sumui = sumui + nf*ui

            !-----    3.3 Weighted sum of increments - V
            IF (lwind) sumvi = sumvi + nf*vi

            !-----    3.4 Weighted sum of squared increments
            sqinc = ui*ui
            IF (lwind) sqinc = sqinc + vi*vi
            sumsqi = sumsqi + nf*sqinc

            !-----    3.5 Determine maximum (squared) increment
            IF (sqinc  >   maxinc) THEN
              maxinc = sqinc
              IF (lwind) THEN
                !         For wind, extreme is maximum speed of vector increment
                rstat(jlev,6) = SQRT(sqinc)
              ELSE
                !         Extreme is increment
                rstat(jlev,6) = ui
              END IF
              !           Lat/Long of obs
              rstat(jlev,7) = wlt(job)
              rstat(jlev,8) = wln(job)
              !-----      Already in degrees
            END IF

          END IF

        END DO !JOB
      END IF  !  if GLOBAL
      s_stat(jlev,0)=maxinc
      s_stat(jlev,1)=REAL(inobslv)
      s_stat(jlev,2)=sumnf
      s_stat(jlev,3)=sumui
      IF (lwind)s_stat(jlev,4)=sumvi
      s_stat(jlev,5)=sumsqi
      s_stat(jlev,6)=rstat(jlev,6)
      s_stat(jlev,7)=rstat(jlev,7)
      s_stat(jlev,8)=rstat(jlev,8)
    END DO ! JLEV
  END IF   ! LENOBT <= 0
  IF (mype == 0) THEN
    !
    !   set rstat to zero
    !
    DO jlev=1,jnanl
      rstat(jlev,1)=0.0
      rstat(jlev,2)=0.0
      rstat(jlev,3)=0.0
      rstat(jlev,4)=0.0
      rstat(jlev,5)=0.0
      a_maxinc(jlev)=0.0
    END DO
  END IF  ! mype=0

  ! Gather stats on PE0
  IF (mype == 0) THEN ! stats off PE0
    DO jlev=1,jnanl
      rstat(jlev,1)=s_stat(jlev,1)
      rstat(jlev,2)=s_stat(jlev,2)
      rstat(jlev,3)=s_stat(jlev,3)
      IF (lwind)rstat(jlev,4)=s_stat(jlev,4)
      rstat(jlev,5)=s_stat(jlev,5)
      IF (s_stat(jlev,0) >  a_maxinc(jlev)) THEN
        a_maxinc(jlev)=s_stat(jlev,0)
        rstat(jlev,6)=s_stat(jlev,6)
        rstat(jlev,7)=s_stat(jlev,7)
        rstat(jlev,8)=s_stat(jlev,8)
      END IF
    END DO
  END IF


  DO iproc=1,nproc-1
    IF (mype == 0) THEN  ! PE0 receives
      CALL gc_rrecv(iproc,9*model_levels_max,iproc,           &
                   istat,r_stat,s_stat)
      DO jlev=1,jnanl
        rstat(jlev,1)=rstat(jlev,1)+r_stat(jlev,1)
        rstat(jlev,2)=rstat(jlev,2)+r_stat(jlev,2)
        rstat(jlev,3)=rstat(jlev,3)+r_stat(jlev,3)
        IF (lwind)rstat(jlev,4)=rstat(jlev,4)+r_stat(jlev,4)
        rstat(jlev,5)=rstat(jlev,5)+r_stat(jlev,5)
        IF (r_stat(jlev,0) >  a_maxinc(jlev)) THEN
          a_maxinc(jlev)=r_stat(jlev,0)
          rstat(jlev,6)=r_stat(jlev,6)
          rstat(jlev,7)=r_stat(jlev,7)
          rstat(jlev,8)=r_stat(jlev,8)
        END IF
      END DO
    ELSE IF (mype == iproc) THEN  ! other PEs send
      CALL gc_rsend(iproc,9*model_levels_max,0,               &
                   istat,r_stat,s_stat)
    END IF
  END DO



  IF (mype == 0) THEN

    !-----  3.6 Get means for level JLEV and accumulate for overall means.
    DO jlev=1,jnanl
      inobslv=INT(rstat(jlev,1))
      sumnf=rstat(jlev,2)
      sumui=rstat(jlev,3)
      IF (lwind)sumvi=rstat(jlev,4)
      sumsqi=rstat(jlev,5)

      IF (inobslv >  0) THEN

        rstat(jlev,1) = REAL( inobslv )
        rstat(0,1)    = rstat(0,1)+rstat(jlev,1)
        IF (sumnf >  0.0) THEN
          !           Mean weight
          rstat(jlev,2) = sumnf/rstat(jlev,1)
          rstat(0,2)    = rstat(0,2)+sumnf
          !           Mean observation increment - P*, T, U, RH
          rstat(jlev,3) = sumui/sumnf
          rstat(0,3)    = rstat(0,3)+sumui
          IF (lwind) THEN
            !             Mean observation increment - V
            rstat(jlev,4) = sumvi/sumnf
            rstat(0,4)    = rstat(0,4)+sumvi
          END IF
          !           Mean square observation increment
          rstat(jlev,5) = sumsqi/sumnf
          rstat(0,5)    = rstat(0,5)+sumsqi
          !           Maximum Increment
          maxinc = ABS(rstat(jlev,6))
          IF (maxinc  >=  ABS(rstat(0,6))) THEN
            rstat(0,6) = rstat(jlev,6)
            rstat(0,7) = rstat(jlev,7)
            rstat(0,8) = rstat(jlev,8)
          END IF
        END IF

      END IF

    END DO

    !       4. Get means for all levels
    !       ---------------------------
    IF (rstat(0,1) >  0.0 .AND. rstat(0,2) >  0.0) THEN

      sumnf = rstat(0,2)
      !         Mean weight
      rstat(0,2) = rstat(0,2)/rstat(0,1)
      !         Mean observation increment - P*, T, U, RH
      rstat(0,3) = rstat(0,3)/sumnf
      !         Mean observation increment - V
      IF (lwind) rstat(0,4) = rstat(0,4)/sumnf
      !         Mean square observation increment
      rstat(0,5) = rstat(0,5)/sumnf

    END IF

    !       If p*, convert from pascals to mb
    IF (anal_var == 1) THEN
      rstat(0,3) = rstat(0,3) * 0.01    ! Mean Increment
      rstat(0,5) = rstat(0,5) * 0.0001  ! Mean Square Increment
      rstat(0,6) = rstat(0,6) * 0.01    ! Maximum Increment
    END IF

    !       Get RMS values if required
    IF (lrms) THEN
      DO jlev=0,jnanl
        IF (rstat(jlev,5) >  0.0) THEN
          rstat(jlev,5) = SQRT (rstat(jlev,5))
        END IF
      END DO
    END IF

  END IF  !mype=0

  !     5. Print out results
  !     --------------------

  IF (mype == 0) THEN
    !---- 5.1 Print out titles

    IF (jband == 1) THEN
      CALL umPrint(' ',src='diago')
      IF (kcall >  1 .AND. kcall  <  6) THEN
        WRITE(umMessage,'(A,A,T30,A,I5)')' DIAGO called from ',            &
            csubr,'Obs Type ',ktype
        CALL umPrint(umMessage,src='diago')
        WRITE(umMessage,'(" ",A9,2X,A25)') title1,title2
        CALL umPrint(umMessage,src='diago')
      ELSE IF (kcall >= 6) THEN
        WRITE(umMessage,'(" Fit to ",A,T22,"cloud (oktas)")') csubr
        CALL umPrint(umMessage,src='diago')
      ELSE
        WRITE(umMessage,'(" DIAGO called from ",A)') csubr
        CALL umPrint(umMessage,src='diago')
      END IF
    END IF
    IF (model_type == mt_global) THEN
      CALL umPrint(' ',src='diago')
      IF (llband) THEN
        WRITE(umMessage,'(A,A3," to ",A3,5X,"No of obs ",I6)')             &
            ' LATITUDES ', bandc(jband),bandc(jband+1),inobs
        CALL umPrint(umMessage,src='diago')
      ELSE
        WRITE(umMessage,'(A,A3," to ",A3)')' LATITUDES ',                  &
            bandc(1),bandc(4)
        CALL umPrint(umMessage,src='diago')
        CALL umPrint(' ',src='diago')
      END IF
    ELSE
      CALL umPrint(' ',src='diago')
      CALL umPrint(' WHOLE LIMITED AREA ',src='diago')
      CALL umPrint(' ',src='diago')
    END IF  ! if GLOBAL

    IF (lrms) THEN
      WRITE(umMessage,'(A,3X,A10,3X,A10,2X,A)')                              &
          ' LEVEL     N  NORM FACTOR',lab1,lab2,                             &
          '    RMS       EXTREME INC & POSITION'
      CALL umPrint(umMessage,src='diago')
    ELSE
      WRITE(umMessage,'(A,3X,A10,3X,A10,2X,A)')                              &
          ' LEVEL     N  NORM FACTOR',lab1,lab2,                             &
          'MEAN SQUARE   EXTREME INC & POSITION'
      CALL umPrint(umMessage,src='diago')
    END IF

    !---- 5.2 Print overall means

    ns = ' '
    IF (rstat(0,7) >  0.0) THEN
      ns = 'N'
    ELSE
      ns = 'S'
    END IF
    we = ' '
    IF (rstat(0,8) >  0.0) THEN
      we = 'E'
    ELSE
      we = 'W'
    END IF
    inobslv   = INT(rstat(0,1))
    lat       = ABS(rstat(0,7))
    longitude = ABS(rstat(0,8))
    WRITE(umMessage,'(A,I6,5E13.5,F5.1,A1,F6.1,A1)') '   ALL',inobslv,       &
        (rstat(0,jvar),jvar=2,6),lat,ns,longitude,we
    CALL umPrint(umMessage,src='diago')

    !---- 5.3 Print means for individual levels

    IF (jnanl >  1) THEN
      DO jlev=jnanl,1,-1
        inobslv = INT( rstat(jlev,1) )
        IF (inobslv >  0) THEN
          ns = ' '
          IF (rstat(jlev,7) >  0.0) THEN
            ns = 'N'
          ELSE
            ns = 'S'
          END IF
          we = ' '
          IF (rstat(jlev,8) >  0.0) THEN
            we = 'E'
          ELSE
            we = 'W'
          END IF
          lat       = ABS(rstat(jlev,7))
          longitude = ABS(rstat(jlev,8))
          WRITE(umMessage,'(2I6,5E13.5,F5.1,A1,F6.1,A1)')                    &
              jlev,inobslv,(rstat(jlev,jvar),jvar=2,6),lat,ns,longitude,we
          CALL umPrint(umMessage,src='diago')
        END IF
      END DO
    END IF

  END IF
END DO ! jband

!   if jumping to 9999 a problem has occurred but treat as non fatal
!   since this is a diagnostic routine (Fixme: needs an ereport?)
9999   CONTINUE

IF (ltimer_ac) CALL timer('DIAGO   ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diago
END MODULE diago_mod
