! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    MODULE GETOB3_MOD
!
!    SUBROUTINE GETOB3
!
!    Purpose : Extract appropriate observations from the OBS array 
!              (passed from RDOBS via argument list)
!              GETOBS called from AC gets list of obs in time window
!              GETOB2 called from AC2 gets lat, long, time and
!                     model obs type for an obs type.
!              GETOB3 called from VERTANL gets data specific to a type
!                     (eg data values and assoc error ratios)
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: AC Assimilation
MODULE getob3_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GETOB3_MOD'

CONTAINS

SUBROUTINE getob3 (kact,obs,obdata,                               &
                   obs_lat,obs_long,                              &
                   lenobt,inobs,ndv,obs_no,icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx, obs_info
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER :: kact                   ! IN  Index to obs type
INTEGER :: lenobt                 ! IN  No of obs to be assimilated
INTEGER :: inobs                  ! IN  No of obs for this type
INTEGER :: ndv                    ! IN  No of data values
!                                    !     excluding header section
REAL :: obs (inobs,*)             ! IN  OBS array from RDOBS
REAL :: obdata (lenobt,ndv)       ! OUT ob values and error ratios
REAL :: obs_lat (lenobt)          ! IN  latitudes
REAL :: obs_long(lenobt)          ! IN  longitudes
INTEGER :: obs_no(lenobt)         ! IN  pointers to obs.
INTEGER :: icode                  ! OUT error code and message
CHARACTER(LEN=errormessagelength) :: cmessage
!-----------------------------------------------------------------------
!  local 
REAL :: zminep
PARAMETER(zminep=1.0e-10)
!     ZMINEP IS THE MINIMUM ALLOWED OBSERVATIONAL ERROR TO AVOID /0.
CHARACTER(LEN=10) :: label
INTEGER :: neropt,ktype,nlev,jdv,job,kdv

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETOB3'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ltimer_ac) CALL timer('GETOB3  ',3)

ktype  = lact (kact)                !  AC Observation Type
nlev   = obs_info % noblev(kact)    !  No of levels for this type

!     Get observation data and error ratios
DO jdv=1,ndv
  DO job=1,lenobt
    obdata(job,jdv) = obs(obs_no(job)-obs_info % obs_no_st(kact), &
                          obs_info % ndvhdr+jdv)
  END DO
END DO

IF (ldiagac) THEN


  IF (ndac >  0) THEN

    !      DO LAT/LON TRANSFORMATIONS AND GET COEFF1,COEFF2 OUTSIDE JDV LOOP

    DO jdv=1,ndv
      IF (ktype == 101) THEN
        IF (jdv == 1) label = 'P* DATA'
        IF (jdv == 2) label = 'EO/EP'
      ELSE IF (ktype == 201) THEN
        IF (jdv <= nlev) THEN
          WRITE(label,21)'T  DATA',jdv
        ELSE
          kdv=jdv-nlev
          WRITE(label,21)'EO/EP  ',kdv
        END IF
      ELSE IF (ktype == 202 .OR. ktype == 204) THEN
        IF (jdv == 1) label = 'T DATA'
        IF (jdv == 2) label = 'EO/EP'
      ELSE IF (ktype == 203) THEN
        IF (jdv == 1) label = 'OBS LEV'
        IF (jdv == 2) label = 'T DATA'
        IF (jdv == 3) label = 'EO/EP'
      ELSE IF (ktype == 205 .OR. ktype == 206 .OR.                       &
              ktype == 207 .OR. ktype == 209 .OR.                       &
              ktype == 211) THEN
        IF (jdv <= nlev) THEN
          WRITE(label,21)'T  DATA',jdv
        ELSE
          kdv=jdv-nlev
          WRITE(label,21)'EO/EP  ',kdv
        END IF
      ELSE IF (ktype == 208) THEN
        IF (jdv <= nlev) THEN
          WRITE(label,21)'T  DATA',jdv
        ELSE IF (jdv >  nlev .AND. jdv <= 2*nlev) THEN
          kdv=jdv-nlev
          WRITE(label,21)'T F/G  ',kdv
        ELSE IF (jdv >  2*nlev .AND. jdv <= 3*nlev) THEN
          kdv=jdv-nlev*2
          WRITE(label,21)'EO/EP  ',kdv
        ELSE
          label='QUALIND'
        END IF
      ELSE IF (ktype == 301 .OR. ktype == 311) THEN
        IF (jdv <= nlev) THEN
          WRITE(label,21)'U  DATA',jdv
        ELSE IF (jdv <= nlev*2) THEN
          kdv=jdv-nlev
          WRITE(label,21)'V  DATA',kdv
        ELSE
          kdv=jdv-nlev*2
          WRITE(label,21)'EO/EP  ',kdv
        END IF
      ELSE IF (ktype == 302 .OR. ktype == 304 .OR.                    &
              ktype == 305 .OR. ktype == 306 ) THEN
        IF (jdv == 1) label = 'U DATA'
        IF (jdv == 2) label = 'V DATA'
        IF (jdv == 3) label = 'EO/EP'
      ELSE IF (ktype == 303) THEN
        IF (jdv == 1) label = 'OBS LEV'
        IF (jdv == 2) label = 'U DATA'
        IF (jdv == 3) label = 'V DATA'
        IF (jdv == 4) label = 'EO/EP'
      ELSE IF (ktype == 401 .OR. ktype == 405 .OR. ktype == 406) THEN
        IF (jdv <= nlev) THEN
          WRITE(label,21)'RH DATA',jdv
        ELSE
          kdv=jdv-nlev
          WRITE(label,21)'EO/EP  ',kdv
        END IF
      ELSE IF (ktype == 402 .OR. ktype == 404) THEN
        IF (jdv == 1) label = 'RH DATA'
        IF (jdv == 2) label = 'EO/EP'
      ELSE IF (ktype == 403) THEN
        IF (jdv == 1) label = 'OBS LEV'
        IF (jdv == 2) label = 'RH DATA'
        IF (jdv == 3) label = 'EO/EP'
      ELSE IF (ktype == 407) THEN
        WRITE(label,21)'CTT CNT',jdv
      ELSE IF (ktype == 506) THEN
        IF (jdv == 1) label = 'PR DATA'
        IF (jdv == 2) label = 'EO/EP'
      ELSE IF (ktype == 901) THEN
        IF (jdv == 1) label = 'VIS DATA'
        IF (jdv == 2) label = 'EO/EP'
      ELSE
        icode=1
        cmessage=' GETOB3 : ILLEGAL OB TYPE'
        GO TO 999
      END IF
      !
      21 FORMAT(a7,i3)
      !
    END DO   !  Loop over JDV
  END IF
END IF


neropt = obs_info % nerlev1(kact) - obs_info % ndvhdr

DO jdv=neropt,neropt+nlev-1

  !       ENSURE ALL OB ERRORS ARE >= ZMINEP TO AVOID
  !       DIVIDING BY ZERO IN VAN## Routines - is this necessary?

  DO job=1,lenobt
    IF (obdata(job,jdv) <  zminep .AND. obdata(job,jdv) /=        &
        obs_info % missd) obdata(job,jdv) = zminep
  END DO

END DO

999   CONTINUE
IF (ltimer_ac) CALL timer('GETOB3  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE getob3
END MODULE getob3_mod
