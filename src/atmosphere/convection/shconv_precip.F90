! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!   Shallow precipitation
!

MODULE shconv_precip_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SHCONV_PRECIP_MOD'
CONTAINS

SUBROUTINE shconv_precip( n_sh, nlev, max_cldlev                  &
,                     zcld, mb, wstar_up, ql_ad, dz_inv, rho_inv  &
,                     eta_half,eta                                &
,                     rainfall, wqr_inv, precip_product_inv       &
,                     wqr, precip_product_th, precip_product_uv )

!
! Purpose:
!   Calculate the shallow precipitation
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 v6 programming standards
!
USE planet_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

INTEGER, INTENT(IN) ::                                            &
   n_sh                                                           &
                   ! No. of shallow convection points
,  nlev                                                           &
                   ! No. model levels
,  max_cldlev      ! Maximum number of cloud levels

REAL, INTENT(IN) ::                                               &
  zcld(n_sh)                                                      &
                         ! depth of convective cloud (m)
, mb(n_sh)                                                        &
                         ! cloud base mass flux
, wstar_up(n_sh)                                                  &
                         ! convective velocity scale, cloud (m/s)
, ql_ad(n_sh)                                                     &
                         ! adiabatic liquid water content (kg/kg)
, dz_inv(n_sh)                                                    &
                         ! dz across inversion (m)
, rho_inv(n_sh)                                                   &
                         ! density at inversion (kg/m3)
, eta_half(n_sh,nlev)                                             &
                         ! eta  uv levels
, eta(n_sh,nlev)         ! eta
!
! Arguments with intent INOUT:
!
!          NONE

!
! Arguments with intent OUT:
!

REAL, INTENT(OUT) ::                                              &
  rainfall(n_sh)                                                  &
                            ! surface rainfall rate
, wqr_inv(n_sh)                                                   &
                            ! flux of w'qr' at inversion
, precip_product_inv(n_sh)                                        &
                            ! Integral of precip_product over inv
                            ! *kappa_rain/density
, wqr(n_sh,nlev)                                                  &
                            !  w'qr'
, precip_product_th(n_sh,nlev)                                    &
                               ! precipitation production rate
, precip_product_uv(n_sh,nlev) ! precipitation production rate


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
REAL  ::                                                          &
 FRACTION(n_sh)                                                   &
                   ! fraction of ensemble with lwc above threshold
,ql_up_inv(n_sh)                                                  &
                  ! mean liquid water content of ensemble at base
                  ! of inversion
,p0(n_sh)                                                         &
                  ! precip productio p0 value
, fr_inv(n_sh)    ! gravitional flux of precipitation at inversion
REAL  ::                                                          &
 fr(n_sh,nlev) ! gravitional flux of precipitation

!
! Loop counters
!

INTEGER :: i,k


REAL, PARAMETER :: ql_t=0.001       ! threshold for precipitation
REAL, PARAMETER :: epsilon_rain=0.5 ! coefficient for surface
                                    ! precipitation calculation.
REAL, PARAMETER :: beta_rain=0.5    ! coefficient for precip
                                    ! production.
REAL, PARAMETER :: gamma_rain=0.5   ! coefficient for rainfall
                                    ! rate at base of inversion.
REAL, PARAMETER :: kappa_rain=1.0   ! coefficient for rainfall
                                    ! rate at base of inversion.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SHCONV_PRECIP'


!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO k=1,nlev
  DO i = 1,n_sh
    precip_product_th(i,k) = 0.0
    precip_product_uv(i,k) = 0.0
    wqr(i,k) = 0.0
  END DO
END DO

! Calculate liquid water content of buoyant updraughts at the base of
! the inversion
DO i = 1, n_sh
  ql_up_inv(i) = (1.5/g)*((wstar_up(i)/mb(i))**0.5)*              &
                 wstar_up(i)*wstar_up(i)/zcld(i)

END DO

! Does precipitation occur ?

DO i = 1, n_sh

  ! case of no precipitation

  IF (ql_t >= ql_ad(i)) THEN
    FRACTION(i) = 0.0
    rainfall(i) = 0.0
    p0(i) = 0.0
    fr_inv(i) = 0.0
    wqr_inv(i) = 0.0
    precip_product_inv(i) = 0.0
  ELSE
    IF (ql_t > ql_up_inv(i)) THEN
      FRACTION(i) = (1.0-ql_t/ql_ad(i))*(1.0-ql_t/ql_ad(i))/       &
                      (1.0-ql_up_inv(i)/ql_ad(i))
    ELSE
      FRACTION(i) = 1.0 - ql_t*ql_t/(ql_up_inv(i)*ql_ad(i))
    END IF

    ! surface precipitation rate

    rainfall(i) = epsilon_rain*mb(i)*ql_up_inv(i)*FRACTION(i)

    ! p0
    p0(i) = rainfall(i)/(beta_rain*(1.0-0.5*beta_rain)              &
                    *zcld(i)*zcld(i) +                             &
                           0.5*beta_rain*zcld(i)*dz_inv(i))
    ! fr at inversion

    fr_inv(i) = gamma_rain*rainfall(i)
    ! wqr at inversion
    !    wqr = rainfall -0.5 p0 beta zcld dz_inv +(gamma-1)rainfall
    ! simplified to
    wqr_inv(i) = gamma_rain*rainfall(i)                                       &
                         - 0.5*p0(i)*beta_rain*zcld(i)*dz_inv(i)


    precip_product_inv(i) =kappa_rain*0.5*p0(i)*beta_rain          &
                                *zcld(i)*dz_inv(i)/rho_inv(i)

  END IF
END DO

!-----------------------------------------------------------------------
! calculate in cloud profiles of precipitation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Rain water flux wqr in cloud integrate precip_product & dfr/dz
!-----------------------------------------------------------------------

DO k=1,max_cldlev+1
  DO i = 1,n_sh
    IF (eta_half(i,k) <= 1.0) THEN

      IF (eta_half(i,k) < beta_rain) THEN
        precip_product_uv(i,k) = eta_half(i,k)*zcld(i)*p0(i)
        fr(i,k) = rainfall(i)
        wqr(i,k) = fr(i,k) - rainfall(i)+precip_product_uv(i,k)   &
                          *0.5*eta_half(i,k)*zcld(i)
      ELSE
        precip_product_uv(i,k) = p0(i)*beta_rain*zcld(i)
        fr(i,k) = rainfall(i)*                                          &
             (1.0 -gamma_rain*beta_rain+eta_half(i,k)*(gamma_rain-1.0)) &
                         /(1.0-beta_rain)
        wqr(i,k) = fr(i,k) - rainfall(i)                             &
                                +p0(i)*beta_rain*zcld(i)*zcld(i)     &
                                      *(eta_half(i,k)-0.5*beta_rain)
      END IF
    ELSE  ! above  cloud
      precip_product_uv(i,k) = 0.0
      fr(i,k) = 0.0
      wqr(i,k) = 0.0
    END IF
  END DO
END DO

! theta levels

DO k=1,max_cldlev+1
  DO i = 1,n_sh
    IF (eta(i,k) <= 1.0) THEN
      IF (eta(i,k) < beta_rain) THEN
        precip_product_th(i,k) = eta(i,k)*zcld(i)*p0(i)
      ELSE
        precip_product_th(i,k) = p0(i)*beta_rain*zcld(i)
      END IF
    ELSE  ! above  cloud
      precip_product_th(i,k) = 0.0
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
!    End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE shconv_precip
END MODULE shconv_precip_mod
