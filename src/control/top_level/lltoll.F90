! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Use subroutine rotate_latlon_to_eq to convert to equatorial lat/lon for LAM
! Subroutine Interface:

SUBROUTINE lltoll(lnthll,lsthll,lestll,lwstll,                    &
                  phi_pole,lambda_pole)

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE latlon_eq_rotation_mod, ONLY: rotate_latlon_to_eq
IMPLICIT NONE

!  Subroutine arguments:

!    Scalar arguments with intent(in):
REAL :: phi_pole    !  Latitude of pole in equatorial system
REAL :: lambda_pole !  Longitude do.

!    Scalar arguments with intent(inout):
INTEGER :: lnthll
INTEGER :: lsthll
INTEGER :: lestll
INTEGER :: lwstll

!  Local parameters:
INTEGER :: points
PARAMETER(points=9)

!  Local arrays:
REAL :: phi      (points)
REAL :: lambda   (points)
REAL :: lambda_eq(points)
REAL :: phi_eq   (points)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LLTOLL'

REAL :: rnthll
REAL :: rsthll
REAL :: restll
REAL :: rwstll

INTEGER :: i

!- End of Header -----------------------------------------------------


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
phi(1)=lnthll
phi(2)=lnthll
phi(3)=lnthll
phi(4)=(lnthll+lsthll)/2
phi(5)=(lnthll+lsthll)/2
phi(6)=lsthll
phi(7)=lsthll
phi(8)=lsthll
phi(9)=(lnthll+lsthll)/2
lambda(1)=lwstll
IF (lwstll <  lestll) THEN
  lambda(2)=(lwstll+lestll)/2
ELSE
  lambda(2)=(lwstll+lestll-360)/2
  IF (lambda(2) <  0) lambda(2)=lambda(2)+360
END IF
lambda(3)=lestll
lambda(4)=lwstll
lambda(5)=lestll
lambda(6)=lwstll
lambda(7)=lambda(2)
lambda(8)=lestll
lambda(9)=lambda(2)

CALL rotate_latlon_to_eq  &
       (phi,lambda,phi_eq,lambda_eq,phi_pole,lambda_pole,points)

IF (lambda_eq(3) <  lambda_eq(2)) lambda_eq(2)=lambda_eq(2)-360.0
IF (lambda_eq(2) <  lambda_eq(1)) lambda_eq(1)=lambda_eq(1)-360.0
IF (lambda_eq(5) <  lambda_eq(9)) lambda_eq(9)=lambda_eq(9)-360.0
IF (lambda_eq(9) <  lambda_eq(4)) lambda_eq(4)=lambda_eq(4)-360.0
IF (lambda_eq(8) <  lambda_eq(7)) lambda_eq(7)=lambda_eq(7)-360.0
IF (lambda_eq(7) <  lambda_eq(6)) lambda_eq(6)=lambda_eq(6)-360.0

rnthll=phi_eq(1)
rsthll=phi_eq(1)
restll=lambda_eq(1)
rwstll=lambda_eq(1)

DO i=2,8
  rnthll=MAX(rnthll,phi_eq(i))
  rsthll=MIN(rsthll,phi_eq(i))
  rwstll=MIN(rwstll,lambda_eq(i))
  restll=MAX(restll,lambda_eq(i))
END DO

IF (rwstll <  0)  rwstll=rwstll+360.0
IF (restll <  0)  restll=restll+360.0
rnthll=rnthll+.999
lnthll=rnthll
lnthll=MIN(90,lnthll)
lsthll=rsthll
lwstll=rwstll
restll=restll+.999
lestll=restll
lestll=MIN(360,lestll)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lltoll
