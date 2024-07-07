! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE VANRAIN-------------------------------------------------
!
!    Purpose : Performs 'vertical' analysis for precip rate/phase data.
!
!       When IPASS=1,
!       model field is interpolated to ob locations and increments
!       calculated. Preliminary weight normalisation factors are also
!       calculated.
!
!       When IPASS=2,
!       data density is interpolated to ob locations and
!       final weight normalisation factors are calculated.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE vanrain_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VANRAIN_MOD'

CONTAINS

SUBROUTINE vanrain (kact,ipass,ls_rain,ls_snow,                  &
                    conv_rain,conv_snow,lenfld,                  &
                    obdata,cf1pt,cf2pt,cf3pt,cf4pt,              &
                    np1pt,np2pt,np3pt,np4pt,                     &
                    princ,normf,obs_lat,obs_long,                &
                    obs_no,lmissd,nptobt,                        &
                    lenobt,ndv,lenob,                            &
                    no_anal_levs,no_anal_var,                    &
                    icode,cmessage)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx, obs_info
USE diagopr_mod, ONLY: diagopr
USE hintmo_mod, ONLY: hintmo
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! ----------------------------------------------------------------------
INTEGER :: lenfld,lenobt,lenob,ndv
INTEGER :: kact,ipass,nptobt,no_anal_levs,no_anal_var
REAL ::                                                           &
               ls_rain(lenfld),ls_snow(lenfld),                   &
               conv_rain(lenfld),conv_snow(lenfld),               &
               obdata(lenobt,ndv),                                &
               princ(lenob),normf(lenob),                         &
               obs_lat(lenobt),obs_long(lenobt),                  &
               cf1pt(lenobt),cf2pt(lenobt),cf3pt(lenobt),         &
               cf4pt(lenobt)
!
LOGICAL :: lmissd(lenob+1,no_anal_levs)
!
INTEGER :: obs_no(lenobt)
INTEGER :: np1pt(lenobt),np2pt(lenobt),np3pt(lenobt)
INTEGER :: np4pt(lenobt)
!
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
!
!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     P_LEVELS  - no of model level
!     LS_RAIN   - }  model precipitation rates (IPASS=1) for
!     LS_SNOW   - }  large-scale and convective rain and snow
!     CONV_RAIN - }  or data density on model grid (IPASS=2)
!     CONV_SNOW - }
!     OBDATA    - observed values
!     OBS_LAT   - ob co-lats
!     OBS_LONG  - ob longitudes
!     OBS_NO    - observation numbers
!     NP1PT,NP2PT,NP3PT,NP4PT,
!               - pointers to nearest model grid points for each ob
!     CF1PT,CF2PT,CF3PT,CF4PT,
!               - weighting of model points in bilinear interpolation
!
!-INTENT=INOUT-----------------------------------------------------
!     PRINC        -  ob-model increments
!     NORMF        -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          - logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!
!---------------------------------------------------------------------
!     Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
REAL ::        intpr(lenobt),totpr(lenfld)
REAL ::        phsint(lenobt),phsfld(lenfld)
REAL ::        wlat(lenobt),wlon(lenobt)
!     INTPR   - model rate field interpolated to ob locations
!     TOTPR   - total rate field, =sum of LS/CONV_RAIN/SNOW
!     WLAT    - Equatorial latitude of obs in degrees
!     WLON    - Equatorial longitude of obs in degrees
!
!
!----------------------------------------------------------------------
!     Define local variables
!-----------------------------------------------------------------------
INTEGER :: ktype,nlev,nobipt,neropt,job,jpt
INTEGER :: jrad
REAL :: ob_to_radar, ob_to_radar_min(lenobt)
REAL :: f2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VANRAIN'

!     KTYPE    -   Observation type
!     NLEV     -   No of observed levels
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JOB      -   Loop counter in loops over obs
!     JPT      -   Loop counter in loops over model points
!     JRAD     -   Loop counter over radars
!     OB_TO_RADAR      - Distance (km) calculated from ob to radar
!     OB_TO_RADAR_MIN  - Min value of OB_TO_RADAR
!----------------------------------------------------------------------
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('VANRAIN ',3)
!
ktype = lact(kact)

nlev  = obs_info % noblev(kact)
IF (ktype == 506) THEN
  nobipt = 1
END IF
neropt = obs_info % nerlev1(kact) - obs_info % ndvhdr
!
!   COMBINE model precipitation rates to give TOTPR
!

IF (ipass  ==  1) THEN
  DO jpt=1,lenfld
    totpr(jpt)=  ls_rain(jpt)   + ls_snow(jpt)                  &
                 + conv_rain(jpt) + conv_snow(jpt)
  END DO
ELSE IF (ipass  ==  2) THEN   ! TOTPR set to obs density
  DO jpt=1,lenfld
    totpr(jpt) = ls_rain(jpt)
  END DO
END IF
!
!
!  ** 1.   OBSERVATIONS GIVING DATA FOR PRECIP RATE/PHASE
!          ----------------------------------------------
!
!  ** 1.1  HORIZONTAL INTERP OF RATE FIELD ON MODEL GRID
!
!
CALL hintmo(totpr,cf1pt,cf2pt,cf3pt,cf4pt,                      &
            np1pt,np2pt,np3pt,np4pt,                            &
            lenfld,1,lenobt,intpr,                              &
            icode,cmessage)
IF (icode >  0) GO TO 9999

!
!  ** 1.2  VERTICAL INTERP/PROCESSING OF MODEL
!          NONE.
!  ** 1.3  SUBTRACTION, CALCULATION OF PRELIM NORM FACTOR (IPASS=1)
!                       CALCULATION OF FINAL NORM FACTOR (IPASS=2)
!
IF (ipass == 1) THEN
  !       PR Increment = PR Observation - Background PR
  DO job=1,lenobt
    princ(nptobt+job) = obdata(job,nobipt) - intpr(job)
  END DO
END IF

!
!     SET WEIGHT FACTOR TO ZERO WHERE OB. VALUE OR ERROR IS FLAGGED
DO job=1,lenobt
  lmissd(nptobt+job,1) = obdata(job,nobipt) == obs_info % missd .OR.  &
                         obdata(job,neropt) == obs_info % missd
END DO
!
IF (ipass == 1) THEN
  IF (mdatadfn == 1) THEN
    DO job=1,lenobt
      IF (lmissd(nptobt+job,1)) THEN
        normf(nptobt+job) = 0.0
      ELSE
        normf(nptobt+job) = 1.0 /                                 &
                  SQRT(obdata(job,neropt)*obdata(job,neropt)+1.0)
      END IF
    END DO !job
  ELSE IF (mdatadfn == 2) THEN
    DO job=1,lenobt
      IF (lmissd(nptobt+job,1)) THEN
        normf(nptobt+job) = 0.0
      ELSE
        normf(nptobt+job) = 1.0 /                                 &
                        ( obdata(job,neropt)*obdata(job,neropt) )
      END IF
    END DO ! job
  END IF
  !
  !
ELSE IF (ipass == 2) THEN
  IF (mdatadfn == 1) THEN
    DO job=1,lenobt
      IF (lmissd(nptobt+job,1)) THEN
        normf(nptobt+job) = 0.0
      ELSE
        normf(nptobt*job) = normf(nptobt+job)*normf(nptobt+job) / &
         ( intpr(job)*normf(nptobt+job) +1.0 -                    &
           normf(nptobt+job)*normf(nptobt+job) )
      END IF
    END DO ! job
  ELSE IF (mdatadfn == 2) THEN
    DO job=1,lenobt
      IF (lmissd(nptobt+job,1)) THEN
        normf(nptobt+job) = 0.0
      ELSE
        normf(nptobt+job) = normf(nptobt+job)/( intpr(job)+1.0 )
      END IF
    END DO ! job
  END IF
  !
END IF
!
!  ** 1.4  VERTICAL INTERP/PROCESSING OF INCS. CALC'N OF WEIGHT FACTORS.
!          NONE.
!
!     DIAGNOSTICS ON FIT OF MODEL PRECIP RATE/PHASE TO OBS
IF (ldiagac .AND. (ipass == 1) ) THEN

  CALL diagopr (lenobt,intpr,obdata,lenobt,lenob,ndv,             &
                lmissd,nobipt,nptobt,no_anal_levs)
END IF

IF (ltimer_ac) CALL timer('VANRAIN ',4)

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vanrain
END MODULE vanrain_mod
