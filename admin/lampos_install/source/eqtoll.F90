! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine EQTOLL-------------------------------------------------
!
!    Purpose:  Calculates latitude and longitude on standard grid
!              from input arrays of latitude and longitude on
!              equatorial latitude-longitude (eq) grid used
!              in regional models. Both input and output latitudes
!              and longitudes are in degrees.
!
!
!
!    Documentation: The transformation formulae are described in
!                   unified model on-line documentation paper S1.
!
!   Logical components covered : S131
!
!   Project task :
!
!
!Code Owner: See Unified Model Code Owner's HTML page
!This file belongs in section: Grids

SUBROUTINE eqtoll                                                 &
(phi_eq,lambda_eq,phi,lambda,phi_pole,lambda_pole,points)


IMPLICIT NONE

INTEGER, INTENT(IN) ::  points   !IN  Number of points to be processed

REAL, INTENT(IN) ::     lambda_eq(points)
                                 !IN  Longitude in equatorial lat-lon coords
REAL, INTENT(IN) ::     phi_eq(points)
                                 !IN  Latitude in equatorial lat-lon coords
REAL, INTENT(IN) ::     phi_pole !IN  Latitude of equatorial lat-lon pole 
REAL, INTENT(IN) ::     lambda_pole  !IN  Longitude of equatorial lat-lon pole
         

REAL, INTENT(OUT) ::   phi(points) !OUT Latitude
REAL, INTENT(OUT) ::   lambda(points) !OUT Longitude (0 =< LON < 360)

! Local varables
REAL   ::  e_lambda,e_phi,a_lambda,arg,a_phi,sin_phi_pole,cos_phi_pole
REAL   ::  term1,term2,lambda_zero
REAL, PARAMETER :: small=1.0e-6
INTEGER  ::  i

REAL, PARAMETER :: pi=3.1415926535879323846
                                 ! Pi                                                                      
REAL, PARAMETER :: PI_OVER_180 =PI/180.0
                                 ! Conversion factor degrees to radians                                
REAL, PARAMETER :: RECIP_PI_OVER_180 = 180.0/PI 
                                 ! Conversion factor radians to degrees    


!----------------------------------------------------------------------

!  1. Initialise local constants

! Latitude of zeroth meridian
lambda_zero=lambda_pole+180.
! Sine and cosine of latitude of eq pole
IF (phi_pole >= 0.0) THEN
  sin_phi_pole =  SIN(pi_over_180*phi_pole)
  cos_phi_pole =  COS(pi_over_180*phi_pole)
ELSE
  sin_phi_pole = -SIN(pi_over_180*phi_pole)
  cos_phi_pole = -COS(pi_over_180*phi_pole)
END IF

!  2. Transform from equatorial to standard latitude-longitude

DO i= 1,points

! Scale eq longitude to range -180 to +180 degs

  e_lambda=lambda_eq(i)
  IF(e_lambda >   180.0) e_lambda=e_lambda-360.0
  IF(e_lambda <  -180.0) e_lambda=e_lambda+360.0

! Convert eq latitude & longitude to radians

  e_lambda=pi_over_180*e_lambda
  e_phi=pi_over_180*phi_eq(i)

! Compute latitude using equation (4.7)

  arg=cos_phi_pole*COS(e_lambda)*COS(e_phi)                       &
                     +SIN(e_phi)*sin_phi_pole
  arg=MIN(arg, 1.0)
  arg=MAX(arg,-1.0)
  a_phi=ASIN(arg)
  phi(i)=recip_pi_over_180*a_phi

! Compute longitude using equation (4.8)

  term1 =(COS(e_phi)*COS(e_lambda)*sin_phi_pole                   &
         -SIN(e_phi)*cos_phi_pole)
  term2=COS(a_phi)
  
  IF(term2 <  small) THEN
    a_lambda=0.0
  ELSE
    arg=term1/term2
    arg=MIN(arg, 1.0)
    arg=MAX(arg,-1.0)
    a_lambda=recip_pi_over_180*ACOS(arg)
    a_lambda=SIGN(a_lambda,e_lambda)
    a_lambda=a_lambda+lambda_zero
  END IF

! Scale longitude to range 0 to 360 degs

  IF(a_lambda >= 360.0) a_lambda=a_lambda-360.0
  IF(a_lambda <  0.0) a_lambda=a_lambda+360.0
  lambda(i)=a_lambda

END DO

RETURN
END SUBROUTINE eqtoll

