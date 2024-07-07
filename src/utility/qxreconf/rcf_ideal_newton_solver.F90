! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_ideal_newton_solver_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_NEWTON_SOLVER_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE newton_solver(p,theta,rho,t0_ref,g,z,                      &
                     intw_rho2w, intw_w2rho, p_0,n,prof_type,u_term,  &
                     tol)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE umPrintMgr

USE ereport_mod, ONLY: ereport
USE rcf_ideal_newton_1d_equations_mod, ONLY: newton_1d_equations
USE rcf_ideal_invert_jacob_mod,        ONLY: invert_jacob
IMPLICIT NONE
!
! Description:
!
! The code solves the 1D initialization problem using
! Newton's method
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments
REAL,                     INTENT(IN)    :: tol  
INTEGER,                  INTENT(IN)    :: n, prof_type
REAL,                     INTENT(IN)    :: intw_rho2w(n,2)
REAL,                     INTENT(IN)    :: intw_w2rho(n,2)
REAL,                     INTENT(IN)    :: p_0
REAL,                     INTENT(IN)    :: u_term(1:n)
REAL,                     INTENT(IN)    :: g(0:n)
REAL,                     INTENT(IN)    :: z(0:n)
REAL,                     INTENT(INOUT) :: t0_ref(0:n)
REAL,                     INTENT(INOUT) :: theta(0:n)
REAL,                     INTENT(INOUT) :: p(0:n+1)
REAL,                     INTENT(INOUT) :: rho(n)

INTEGER                                 :: i, k
REAL                                    :: err
REAL                                    :: x(3*n+3)
REAL                                    :: f(3*n+3)
REAL                                    :: dx(3*n+3)
INTEGER                                 :: it
REAL                                    :: err_low

INTEGER                                 :: errcode

INTEGER, PARAMETER                      :: itmax = 1000
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NEWTON_SOLVER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial guess (interleave the variables to
! reduce the bandwidth of the Jacobian matrix)

k = 1
DO i = 0, n
  IF ( i == 0 ) THEN
    x(k) = rho(1)
  ELSE
    x(k) = rho(i)
  END IF
  k    = k + 1
  x(k) = p(i)
  k    = k + 1
  x(k) = theta(i)
  k    = k + 1
END DO

! Perform Newton iterations until converged


dx = 0.0
err_low=HUGE(err_low)
DO it = 1, itmax

  CALL newton_1d_equations(f,x,t0_ref,g,z,intw_w2rho,                &
                       p_0,n,prof_type, u_term)

  err = MAXVAL(ABS(f))+MAXVAL(ABS(dx))
  IF (err < err_low) err_low=err
  IF ( err < tol ) EXIT

  CALL invert_jacob(dx,f,x,g,z,                                      &
                       intw_rho2w, intw_w2rho,n,prof_type)!, u_term)

  x = ABS(x - dx)
END DO

IF (it >= itmax .AND. err > tol ) THEN
   WRITE(umMessage,'(A)') 'Setup failed'
   CALL umPrint(umMessage,src='newton_solver')
   WRITE(umMessage,'(A,2E16.8,I6)')                                &
        'Lowest Error, Tolerance, Max Iterations:',                &
        err_low, tol, itmax
   CALL umPrint(umMessage,src='newton_solver')
   errcode = 1
   CALL ereport("newton_solver", errcode,                          &
               "Newton iteration failed to converge" )
END IF
! Copy solution back (discarding the surface rho value)

k = 1
DO i = 0, n
  IF ( i /= 0 ) rho(i)   = x(k)
  k        = k + 1
  p(i)     = x(k)
  k        = k + 1
  theta(i) = x(k)
  k        = k + 1
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE newton_solver
END MODULE rcf_ideal_newton_solver_mod
