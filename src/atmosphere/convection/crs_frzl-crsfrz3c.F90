! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change of phase routine where precip crosses a melting or freezing level
!
MODULE crs_frzl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CRS_FRZL_MOD'
CONTAINS

SUBROUTINE crs_frzl (npnts, bddwt_km1, exk, exkm1, flx_dd_km1,             &
                     thdd_k, thdd_km1, rain, snow)

USE cv_run_mod, ONLY:                                                      &
    l_snow_rain, t_melt_snow

USE planet_constants_mod, ONLY: g, cp
USE water_constants_mod, ONLY: lf, tm

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!
! Description: Change of phase routine where precipitation crosses a melting
!              or freezing level (The code of this routine is almost identical
!              to chg_phse. At some future date it would be better to have
!              one common routine if possible.)
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards 8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
INTEGER, INTENT(IN) :: &
  npnts                  ! Vector length

LOGICAL, INTENT(IN) :: &
  bddwt_km1(npnts)       ! Mask for those points in downdraught where
                         ! precipitation is liquid in layer k-1

REAL, INTENT(IN) ::    &
  exk(npnts)           & ! Exner ratio in layer k
 ,exkm1(npnts)         & ! Exner ratio in layer k-1
 ,flx_dd_km1(npnts)      ! Downdraught mass flux in layer k-1 (Pa/s)

REAL, INTENT(IN) ::    &
  thdd_k(npnts)          ! Potential temperature of DD in layer k (K)

REAL, INTENT(INOUT) :: &
  thdd_km1(npnts)      & ! In  Potential temperature of DD in layer k-1(K)
                         ! Out Potential temperature of DD in layer k-1
                         ! updated due to change of phase
 ,rain(npnts)          & ! In  Amount of rain descending from k-1 to k-2
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)            ! In  Amount of snow descending from k-1 to k-2
                         ! Out Updated amount of snowfall (kg/m**2/s)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i          ! Loop counter

REAL ::            &
  factor           &  ! Used in the calculation of the change of phase of
                      ! falling precipitation.
 ,melt_factor      &  ! proportion to try to melt.
 ,rain_temp        &  ! rain amount
 ,precip_fre       &  ! Freezing precipitation
 ,precip_melt      &  ! Melting precipitation
 ,precip_melt2        ! Melting precipitation revised estimate of melting snow
REAL :: melt_tscale   ! e-folding temperature thickness for snow melt
REAL :: r_melt_tscale2! reciprocal squared of e-folding temperature thickness
                      ! for snow melt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CRS_FRZL'
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Section (11), equation (42)
!-----------------------------------------------------------------------

IF (l_snow_rain) THEN
  ! New method to melt snow at at rate proportional to the temperature 
  ! above the melting temperature (usually 0C). Note that freezing of
  ! supercooled rain is not considered under this option.
  melt_tscale   = t_melt_snow - tm
  r_melt_tscale2= 1./melt_tscale**2
  DO i=1,npnts

    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)

    IF ( (thdd_km1(i)*exkm1(i) > tm) .AND. (snow(i) > 0.0) ) THEN
      ! Melt
      IF ( (3.0*melt_tscale) > (thdd_km1(i)*exkm1(i)-tm) ) THEN
        ! If temperature still below emergency melting level
        !then calculate melting rate
        precip_melt = snow(i)*(thdd_km1(i)*exkm1(i)-tm)                     &
                    *(thdd_km1(i)*exkm1(i)-thdd_k(i)*exk(i))*r_melt_tscale2
      ELSE
        ! If above emergency melting level melt all snow
        precip_melt = snow(i)
      END IF

      ! Check that precip_melt is not negative.
      ! to prevent refreezing of rain
      precip_melt     = MAX(0.0,precip_melt)

      ! Limit the amount of melting so that is cannot
      ! lower the temperature to below tm.
      precip_melt2    = (thdd_km1(i)-tm/exkm1(i))*factor
      precip_melt     = MIN(precip_melt2,precip_melt)

      ! Limit the melting to the amount of snow available
      precip_melt     = MIN(snow(i),precip_melt)

      ! Calculate temperature change and update rain and snow
      thdd_km1(i)     = thdd_km1(i)-precip_melt/factor
      rain(i)         = rain(i)+precip_melt
      snow(i)         = snow(i)-precip_melt

    END IF  ! test on melting/freezing

  END DO

ELSE
  ! original code melts snow at freezing level and immediately 
  ! freezes any supercooled rain.
  DO i=1,npnts

    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)

    IF (.NOT. bddwt_km1(i) .AND. rain(i) >  0.0 .AND. thdd_km1(i)       &
         *exkm1(i) <  tm) THEN
      ! Freeze
      precip_fre = (tm/exkm1(i)-thdd_km1(i))* factor
      precip_fre = MIN(rain(i),precip_fre)
      thdd_km1(i) = thdd_km1(i)+precip_fre/factor
      rain(i) = rain(i)-precip_fre
      snow(i) = snow(i)+precip_fre

    ELSE IF (bddwt_km1(i) .AND. snow(i) >  0.0) THEN
      ! Melt
      precip_melt = (thdd_km1(i)-tm/exkm1(i))*factor
      precip_melt = MIN(snow(i),precip_melt)
      thdd_km1(i) = thdd_km1(i)-precip_melt/factor
      rain(i) = rain(i)+precip_melt
      snow(i) = snow(i)-precip_melt

    END IF  ! test on melting/freezing

  END DO

END IF  !l_snow_rain

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE crs_frzl

END MODULE crs_frzl_mod
