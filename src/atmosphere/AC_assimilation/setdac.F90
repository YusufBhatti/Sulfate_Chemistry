! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE SETDAC -------------------------------------------------
!
!    Purpose : Controls diagnostic output from AC Scheme.
!
!              Diagnostic options are set through ADIAG namelist
!              which is read in INITAC.
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE setdac_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETDAC_MOD'

CONTAINS

SUBROUTINE setdac(obs,tndv)
!     -----------------

!   CALLED BY SUBROUTINE RDOBS
!   THE DIAGNOSTICS ARE CONTROLLED BY VARIABLES IN COMACDG, THESE ARE
!   DOCUMENTED IN DOCACDG

!   THIS SUBROUTINE CALLS NO OTHERS


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx, obs_info
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER :: tndv
REAL :: obs(tndv)
INTEGER :: j,job,jobt,inobs
INTEGER :: ip_lat,ip_long
REAL :: lat,long

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETDAC'

!
!     3. SET UP FOR DIAGNOSTIC TYPE 1:- DETAILS ABOUT SELECTED OBS.
!     -------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ndaco = 0
ndac  = 0
IF (ldiagac) THEN
  IF (lldac(1)) THEN

    IF (modaco == 1) THEN

      IF (ndacprt >  ndacp) ndacprt=ndacp
      IF (mype == 0) THEN
        WRITE(umMessage,*)                                                    &
        'MODACO=1 : Diagnostics on first ',ndacprt,                 &
        ' Observations of each type.'
        CALL umPrint(umMessage,src='setdac')
      END IF

    ELSE IF (modaco == 2) THEN

      !           MDACO SET IN NAMELIST. COUNT THE NUMBER SET & CLEAR REST.
      DO j=1,ndacop
        IF (mdaco(j) >  0) ndaco = ndaco+1
      END DO
      lldac(1) = ndaco >  0
      IF (ndaco <  ndacop) THEN
        DO j=ndaco+1,ndacop
          mdaco(j) = 0
        END DO
      END IF

    ELSE IF (modaco == 3) THEN

      !           MDACO SET UP FROM OBS IN LTD AREA DEFINED IN NAMELIST.
      !           CHECK LATITUDE BAND (REMEMBER COLATITUDES ARE USED)
      !           CHECK LONGITUDE BAND (ALLOW FOR AREAS CROSSING LONGITUDE 0)

      inobs=0

      DO jobt=1,obs_info % nobtyp

        IF (obs_info % nobs(jobt) >  0) THEN

          ip_lat  = obs_info % mdispobt(jobt)
          ip_long = obs_info % mdispobt(jobt) + obs_info % nobs(jobt)

          IF (dlonge >= dlongw) THEN

            DO job=1,obs_info % nobs(jobt)
              lat  = obs(ip_lat +job)
              long = obs(ip_long+job)
              IF ( lat  >= dlatn  .AND. lat <= dlats  .AND.         &
                   long >= dlongw .AND. long <= dlonge ) THEN
                inobs = inobs+1
                IF (inobs <= ndacop) mdaco(inobs) =                   &
                                     obs_info % obs_no_st(jobt)+job
              END IF
            END DO

          ELSE

            !               AREA ASSUMED TO CROSS MERIDIAN (ELSE MISTAKE IN SET UP)
            DO job=1,obs_info % nobs(jobt)
              lat  = obs(ip_lat +job)
              long = obs(ip_long+job)
              IF ( lat >= dlatn   .AND. lat <= dlats  .AND.         &
                 (long >= dlongw  .OR.  long <= dlonge) ) THEN
                inobs = inobs+1
                IF (inobs <= ndacop) mdaco(inobs) =                   &
                                     obs_info % obs_no_st(jobt)+job
              END IF
            END DO

          END IF
        END IF
      END DO

      ndaco    = MIN(ndacop,inobs)
      lldac(1) = ndaco >  0
      IF (mype == 0) THEN
        IF (ndaco <  inobs) THEN
          WRITE(umMessage,*) 'MODACO=3 ',inobs,                               &
          ' Observations found in limited area for diagnostics'
          CALL umPrint(umMessage,src='setdac')
          WRITE(umMessage,*) '         ',ndaco,                               &
          ' of these listed and used.'
          CALL umPrint(umMessage,src='setdac')
        ELSE
          WRITE(umMessage,*) 'MODACO=3 ',inobs,                               &
          ' Observations found in limited area for diagnostics'
          CALL umPrint(umMessage,src='setdac')
        END IF
      END IF

    END IF
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setdac
END MODULE setdac_mod
