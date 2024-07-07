! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change phase for points where no downdraught occurring.
!
MODULE chg_phse_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CHG_PHSE_MOD'
CONTAINS

SUBROUTINE chg_phse (npnts,k,rain,snow,dthbydt_km1,               &
                     exk,exkm1,delpkm1,the_k,the_km1,timestep,    &
                     cca)


USE cv_run_mod, ONLY:                                             &
    l_snow_rain, t_melt_snow
USE cv_param_mod, ONLY:                                           &
     precip_area_fac

USE planet_constants_mod, ONLY: g, cp
USE water_constants_mod, ONLY: lf, tm

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Change phase for points where no downdraught occurring. (The code of this
! routine is almost identical to crs_frzl. At some future date it would be
! better to have one common routine if possible.)
! Updates potential temperature of layer k as precipitation changes phase
! in situ.
! Add latent heating where precipitation crosses a melting or freezing level.
!
! See UM Documentation paper 27.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,k                      ! Model layer

REAL, INTENT(IN) ::    &
  timestep               ! Timestep

REAL, INTENT(IN) :: &
  exk(npnts)        & ! Exner ratio for layer k
 ,exkm1(npnts)      & ! Exner ratio for layer k-1
 ,delpkm1(npnts)    & ! Pressure difference across layer k-1 (Pa)
 ,the_k(npnts)      & ! potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! potential temperature of environment in layer k-1 (K)
 ,cca(npnts)          ! Convective cloud amount


!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
REAL, INTENT(INOUT) :: &
  rain(npnts)          & ! IN Amount of falling rain (kg/m**2/s)
                         ! OUT Updated amount of falling rain (kg/m**2/s)
 ,snow(npnts)          & ! IN Amount of falling Snowfall (kg/m**2/s)
                         ! OUT Updated amount of falling snow (kg/m**2/s)
 ,dthbydt_km1(npnts)     ! IN Increment to model potential temperature in
                         !    layer k-1
                         ! OUT Updated increment to model potential
                         !     temperature in layer k-1 due to change of phase

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        &
  i                 ! loop counter

LOGICAL ::       &
  bppnwt_k       & ! Mask where precipitation is liquid in layer k
 ,bppnwt_km1       ! Mask where precipitation is liquid in layer k-1

REAL ::           &
  factor          & ! Used in the calculation of change of phase of falling
                    ! precipitation.
 ,melt_factor     & ! proportion to try to melt.
 ,rain_temp       & ! rain amount
 ,precip_fre      & ! Freezing precipitation
 ,precip_melt     & ! Amount of precipitation melting
 ,precip_melt2    & ! Revised amount of precip to change phase.
 ,precip_area     & ! Area fraction of precip = cca * precip_area_fac
 ,the_km1_new       ! the_km1 updated with all increments prior to this
                    ! subroutine
REAL :: melt_tscale   ! e-folding temperature thickness for snow melt
REAL :: r_melt_tscale2! reciprocal squared of e-folding temperature thickness
                      ! for snow melt


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHG_PHSE'
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! ----------------------------------------------------------------------
!   Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Section (11), equation (42)
! ----------------------------------------------------------------------
IF (l_snow_rain) THEN
  ! New method to melt snow at at rate proportional to the temperature 
  ! above the melting temperature (usually 0C). Note that freezing of
  ! supercooled rain is not considered under this option.
  melt_tscale   = t_melt_snow - tm
  r_melt_tscale2= 1./melt_tscale**2
  DO i=1,npnts
    the_km1_new = the_km1(i)  + timestep*dthbydt_km1(i)
    bppnwt_km1  = the_km1_new >   tm/exkm1(i)
    factor      = lf*g/(exkm1(i)*cp*delpkm1(i))
    precip_area = precip_area_fac*cca(i)

    IF ( bppnwt_km1 .AND. (snow(i) > 0.0) ) THEN
    
      IF ( (3.0*melt_tscale) > (the_km1_new*exkm1(i)-tm) ) THEN
        ! If temperature still below emergency melting level
        !then calculate melting rate
        precip_melt = snow(i)*(the_km1_new*exkm1(i)-tm)                     &
                    *(the_km1_new*exkm1(i) - the_k(i)*exk(i))*r_melt_tscale2
      ELSE
        ! If above emergency melting level melt all snow
        precip_melt = snow(i)
      END IF

      ! Check that precip_melt is not negative.
      ! to prevent refreezing of rain
      precip_melt     = MAX(0.0,precip_melt)

      ! Limit the amount of melting so that is cannot
      ! lower the temperature to below tm.
      precip_melt2    = precip_area*(the_km1_new - tm/exkm1(i)) / &
                                                        (timestep * factor)
      precip_melt     = MIN(precip_melt2,precip_melt)

      ! Limit the melting to the amount of snow 
      precip_melt     = MIN(snow(i),precip_melt)
      
      !Calculate temperature increment and update rain and snow
      dthbydt_km1(i)  = dthbydt_km1(i) - precip_melt * factor
      rain(i)         = rain(i)+precip_melt
      snow(i)         = snow(i)-precip_melt
      
    END IF
  END DO
  
ELSE
  ! original code melts all snow at freezing level and immediately 
  ! freezes any supercooled rain.
  DO i=1,npnts
    the_km1_new = the_km1(i) + timestep*dthbydt_km1(i)
    bppnwt_k    = the_k(i) >  tm/exk(i)
    bppnwt_km1  = the_km1_new  >   tm/exkm1(i)
    factor      = lf*g/(exkm1(i)*cp*delpkm1(i))
    precip_area = precip_area_fac*cca(i)

    IF (.NOT. bppnwt_km1 .AND. (bppnwt_k .OR. rain(i) >  0.0)) THEN
      ! Freeze
      precip_fre = MIN( rain(i), precip_area*(tm/exkm1(i) - the_km1_new) / &
                                                         (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) + precip_fre * factor
      snow(i) = snow(i) + precip_fre
      rain(i) = rain(i) - precip_fre
    END IF


    IF (bppnwt_km1 .AND. (.NOT. bppnwt_k .OR. snow(i) >  0.0)) THEN
      ! Melt
      precip_melt = MIN( snow(i),                                 &
               precip_area * (the_km1_new - tm/exkm1(i)) /        &
                  (timestep * factor) )
      dthbydt_km1(i) = dthbydt_km1(i) - precip_melt * factor
      rain(i) = rain(i) + precip_melt
      snow(i) = snow(i) - precip_melt

    END IF ! test on melting /freezing
  END DO

END IF  !l_snow_rain

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE chg_phse

END MODULE chg_phse_mod
