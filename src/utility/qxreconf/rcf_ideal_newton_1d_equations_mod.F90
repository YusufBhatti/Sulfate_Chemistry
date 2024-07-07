! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_ideal_newton_1d_equations_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_NEWTON_1D_EQUATIONS_MOD'

CONTAINS

SUBROUTINE newton_1d_equations(f,x,t0,g,z,intw_w2rho,p_0,                    &
                               n,prof_type, u_term)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: r, cp, kappa, p_zero

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
!
! Description: Evaluate the 1D equations for the Newton solver.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.

! Subroutine arguments
INTEGER,                   INTENT(IN)   :: n, prof_type
REAL,                      INTENT(IN)   :: p_0
REAL,                      INTENT(IN)   :: intw_w2rho(n,2)
REAL,                      INTENT(OUT)  :: f(3*n+3)
REAL,                      INTENT(IN)   :: x(3*n+3)
REAL,                      INTENT(IN)   :: z(0:n)
REAL,                      INTENT(IN)   :: t0(0:n)
REAL,                      INTENT(IN)   :: g(0:n)
REAL,                      INTENT(IN)   :: u_term(1:n)

REAL                                    :: gma, tmp, dens, dz
REAL                                    :: kp2
INTEGER                                 :: i, k
INTEGER                                 :: errorstatus
CHARACTER(LEN=errormessagelength)       :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NEWTON_1D_EQUATIONS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

gma = r/p_zero
kp2 = (1.0-kappa)/kappa

! Surface equations

f(1) = x(2)**kp2 - gma*x(1)*x(3)
f(2) = x(2) - (p_0/p_zero)**kappa
SELECT CASE(prof_type)
CASE (1)
  f(3) = x(3) - t0(0)
CASE (2)
  f(3) = x(2)*x(3) - t0(0)
CASE DEFAULT

  errorstatus = 10
  WRITE(cmessage, '(A,I0)') 'Invalid profile type: ',prof_type
  CALL ereport("newton_1d_equations", errorstatus, cmessage)
END SELECT

k = 4
DO i = 1, n
  dz   = z(i) - z(i-1)

  ! Equation of state

  tmp  = intw_w2rho(i,1)*x(3*i+3) + intw_w2rho(i,2)*x(3*i)
  f(k) = x(3*i+2)**kp2 - gma*x(3*i+1)*tmp
  k    = k + 1

  ! Hydrostatic balance

  IF ( i == 1 ) THEN
    dens = x(3)
  ELSE
    dens = x(3*i)
  END IF

  tmp  = cp/(g(i-1))

  f(k) = tmp*( x(3*i+2) - x(3*i-1) )/dz + 1.0/dens
  ! Modification for quasi hydrostatic balance
  f(k) = f(k) - u_term(i)/(g(i-1)*dens)
  k    = k + 1

  ! Constraint equation

  SELECT CASE(prof_type)
  CASE (1)
    f(k) = x(3*i+3) - t0(i)
  CASE (2)
    f(k) = x(3*i+2)*(intw_w2rho(i,1)*x(3*i+3)                    &
                    +intw_w2rho(i,2)*x(3*i))- t0(i)
  END SELECT
  k = k + 1
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE newton_1d_equations
END MODULE rcf_ideal_newton_1d_equations_mod
