! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates forced detrainment rate in layer K

MODULE det_rate_6a_mod

IMPLICIT NONE

!
! Description:
!   Calculates forced detrainment rate in layer K
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DET_RATE_6A_MOD'

CONTAINS

SUBROUTINE det_rate_6a (npnts, exkp1,                         &
                        thek, thekp1, thpk, thpkp1, thrk,     &
                        qek, qekp1, qpk, qpkp1, qrk,          &
                        Qlkp1, Qfkp1, Frezkp1,                &
                        ekp14, ekp34,                         &
                        deltak)

USE water_constants_mod, ONLY: lc, lf
USE planet_constants_mod, ONLY: cp
USE cv_derived_constants_mod, ONLY: ls
USE cv_run_mod, ONLY: fdet_opt
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: npnts     ! Number of points

REAL,INTENT(IN) :: exkp1(npnts)  ! Exner ratio at mid-point of layer k+1
REAL,INTENT(IN) :: thek(npnts)   ! Env. potential temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts) ! Env. potential temperature in layer k+1
REAL,INTENT(IN) :: thpk(npnts)   ! Par. potential temperature in layer k (K)
REAL,INTENT(IN) :: thpkp1(npnts) ! Par. potential temperature in layer k+1 (K)
REAL,INTENT(IN) :: thrk(npnts)   ! p. temperature of forced detrained
                                 ! parcel in layer k (K)
REAL,INTENT(IN) :: qek(npnts)    ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)  ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qpk(npnts)    ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qpkp1(npnts)  ! Par. specific humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qrk(npnts)    ! Specific humidity of forced detrained
                                 ! parcel in layer k (kg/kg)
REAL,INTENT(IN) :: Qlkp1(npnts)  ! Amount of condensation to liquid water
                                 ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Qfkp1(npnts)  ! Amount of deposition to ice water
                                 ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Frezkp1(npnts)! Amount of freezing from liquid
                                 ! to ice water in the parcel (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)  ! Entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)  ! Entrainment coefficient at level k+3/4
                                 ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: deltak(npnts)! Parcel forced detrainment rate in
                                 ! layer k multiplied by layer thickness

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
REAL :: deltaqk(npnts)           ! Parcel forced detrainment rate based on
                                 ! specific humidity eqn in
                                 ! layer k multiplied by layer thickness
REAL :: deltathk(npnts)          ! Parcel forced detrainment rate based on
                                 ! potential temperature eqn in
                                 ! layer k multiplied by layer thickness

INTEGER :: i                     ! loop counter
INTEGER :: errorstatus           ! error status

REAL :: Factor      ! factor used to calculate the detrainment rate.
REAL :: Denom       ! denominator used to calculate the detrainment rate.
REAL :: Numer       ! numerator used to calculate the detrainment rate.
REAL, PARAMETER :: eps=EPSILON(Denom) !Smallest allowable denominator

CHARACTER(LEN=errormessagelength) :: cmessage    ! error message
CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'DET_RATE_6A'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------
! Calculate the forced detrainment rates for humidity and theta
! Both versions are required because checks are always made on
! whether either is outside of their 0.0 to 1.0 limits.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Calculate the forced detrainment rate for humidity
!---------------------------------------------------------------------
DO i=1,npnts
  Factor = ekp14(i)*qek(i) + (1.0+ekp14(i))*ekp34(i)*qekp1(i)             &
           - (1.0+ekp14(i))*(1.0+ekp34(i))                                &
           * (qpkp1(i) + Qlkp1(i) + Qfkp1(i))

  Numer  = qpk(i) + Factor
  Denom  = qrk(i) + Factor

  IF (ABS(Denom) <= eps*ABS(Numer)) THEN
    !If the denominator is close to zero then set
    !the detrainment rate to one (plus a bit to ensure
    !that this is taken account of later).
    deltaqk(i) = 1.1
  ELSE
    deltaqk(i) = Numer/Denom
  END IF

END DO

!---------------------------------------------------------------------
! Calculate the forced detrainment rate for theta
!---------------------------------------------------------------------
DO i=1,npnts
  Factor = ekp14(i)*thek(i) + (1.0+ekp14(i))*ekp34(i)*thekp1(i)           &
           - (1.0+ekp14(i))*(1.0+ekp34(i))                                &
           * (thpkp1(i) - (lc*Qlkp1(i) + ls*Qfkp1(i) + lf*Frezkp1(i))     &
           / (cp*exkp1(i)))

  Numer  = thpk(i) + Factor
  Denom  = thrk(i) + Factor

  IF (ABS(Denom) <= eps*ABS(Numer)) THEN
    !If the denominator is close to zero then set
    !the detrainment rate to one (plus a bit to ensure
    !that this is taken account of later).
    deltathk(i) = 1.1
  ELSE
    deltathk(i) = Numer/Denom
  END IF

END DO

!---------------------------------------------------------------------
! Copy the either the humidity or theta based detrainment rate to the
! output detrainment rate
!---------------------------------------------------------------------
IF (fdet_opt == 0) THEN
  DO i=1,npnts
    deltak(i) = deltaqk(i)
  END DO
ELSE IF (fdet_opt == 1 .OR. fdet_opt == 2 .OR. fdet_opt == 3) THEN
  DO i=1,npnts
    deltak(i) = deltathk(i)
  END DO
ELSE
  errorstatus = 1   ! will cause model to fail
  WRITE (cmessage,'(a48)') 'Invalid value for fdet_opt. Valid values=0,1,2,3'
  CALL ereport(routinename, errorstatus, cmessage)
END IF

IF (fdet_opt == 0 .OR. fdet_opt == 1) THEN
!---------------------------------------------------------------------
! But if either deltaqk or deltathk are outside of the acceptable
! range then reset to either 0.0 or 1.0 as appropriate.
!---------------------------------------------------------------------
  DO i=1,npnts
    IF (deltaqk(i) <= 0.0 .OR. deltathk(i) <= 0.0) THEN
      deltak(i) = 0.0
    ELSE IF (deltaqk(i) >= 1.0 .OR. deltathk(i) >= 1.0) THEN
      deltak(i) = 1.0
    END IF
  END DO
ELSE IF (fdet_opt == 2 .OR. fdet_opt == 3) THEN
!---------------------------------------------------------------------
! If deltathk are outside of the acceptable
! range then reset to either 0.0 or 1.0 as appropriate.
!---------------------------------------------------------------------
  DO i=1,npnts
    IF (deltathk(i) <= 0.0) THEN
      deltak(i) = 0.0
    ELSE IF (deltaqk(i) >= 1.0) THEN
      deltak(i) = 1.0
    END IF
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE det_rate_6a
END MODULE det_rate_6a_mod
