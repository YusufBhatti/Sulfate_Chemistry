! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE rcf_ideal_invert_jacob_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_INVERT_JACOB_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE invert_jacob(dx,f,x,g,z,intw_rho2w,intw_w2rho,       &
                        n,prof_type)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: r, cp, kappa, p_zero

IMPLICIT NONE
!
! Description:
!
!  Construct the Jacbian of the 1D problem and invert
!  using the matrix of cofactors.
!  NOTE: the prblem is in (3x3) block form with
!        Jd containing the diagonal 3x3 entries
!        Jl containing the lower band of 3x3 entries
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! Subroutine arguments

INTEGER, INTENT(IN)   :: n, prof_type
REAL    ,INTENT(IN)   :: f(3*n+3),                                    &
                         x(3*n+3)
REAL    ,INTENT(IN)   :: g(0:n),                                      &
                         z(0:n)
REAL    ,INTENT(IN)   :: intw_rho2w(n,2),                             &
                         intw_w2rho(n,2)
REAL    ,INTENT(OUT)  :: dx(3*n+3)


! Local variables
REAL ::  kp2

REAL ::  jd(3,3,n+1),                                                 &
         jl(3,3,n+1)
REAL ::  cof(2,2)
REAL ::  u(3,3)
REAL ::  dt(3)
REAL ::  gma,                                                         &
         dz
REAL ::  tmp,                                                         &
         t1,                                                          &
         t2
REAL ::  det,                                                         &
         sgn
INTEGER :: i,                                                           &
         j,                                                           &
         k,                                                           &
         l

INTEGER :: ii,                                                          &
         jj,                                                          &
         kk

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INVERT_JACOB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

gma = r/p_zero
kp2 = (1.0-kappa)/kappa

! First calculate nonzero entries in the Jacobian matrix

jd(1,1,1) = -gma*x(3)
jd(1,2,1) = kp2*x(2)**(kp2-1.0)
jd(1,3,1) = -gma*x(1)
jd(2,1,1) = 0.0
jd(2,2,1) = 1.0
jd(2,3,1) = 0.0
SELECT CASE(prof_type)
CASE (1)
  jd(3,1,1) = 0.0
  jd(3,2,1) = 0.0
  jd(3,3,1) = 1.0
CASE (2)
  jd(3,1,1) = 0.0
  jd(3,2,1) = x(3)
  jd(3,3,1) = x(2)
END SELECT

DO i = 1, n

  jd(1,1,i+1) = -gma*( intw_w2rho(i,1)*x(3*i+3) +                    &
                       intw_w2rho(i,2)*x(3*i) )
  jd(1,2,i+1) = kp2*x(3*i+2)**(kp2-1.0)
  jd(1,3,i+1) = -gma*intw_w2rho(i,1)*x(3*i+1)

  jl(1,1,i+1) = 0.0
  jl(1,2,i+1) = 0.0
  jl(1,3,i+1) = -gma*intw_w2rho(i,2)*x(3*i+1)

  dz   = g(i-1)*( z(i) - z(i-1) )/cp
  IF ( i == 1 ) THEN
    tmp  = 1.0
    jd(2,1,i+1) = 0.0
    jl(2,1,i+1) = 0.0
  ELSE
    t1   = intw_w2rho(i,1)*x(3*i+3) + intw_w2rho(i,2)*x(3*i)
    t2   = intw_w2rho(i-1,1)*x(3*i) + intw_w2rho(i-1,2)*x(3*i-3)
    tmp  = x(3*i)/(intw_rho2w(i-1,1)*t1 + intw_rho2w(i-1,2)*t2)
    tmp  = 1.0
    jd(2,1,i+1) = 0.0
    jl(2,1,i+1) = 0.0
  END IF

  jd(2,2,i+1) = tmp/dz
  jd(2,3,i+1) = -1.0/x(3*i)**2

  jl(2,2,i+1) =-tmp/dz
  jl(2,3,i+1) = 0.0

  SELECT CASE(prof_type)
  CASE (1)
    jd(3,1,i+1) = 0.0
    jd(3,2,i+1) = 0.0
    jd(3,3,i+1) = 1.0
  CASE (2)
    jd(3,1,i+1) = 0.0
    jd(3,2,i+1) = intw_w2rho(i,1)*x(3*i+3)                       &
                + intw_w2rho(i,2)*x(3*i)
    jd(3,3,i+1) = intw_w2rho(i,1)*x(3*i+2)
  END SELECT
  jl(3,1,i+1) = 0.0
  jl(3,2,i+1) = 0.0
  jl(3,3,i+1) = intw_w2rho(i,2)*x(3*i+2)
END DO

! Now invert the D (diagonal) 3x3 matrices

DO kk = 1, n+1

  det = jd(1,1,kk)*jd(2,2,kk)*jd(3,3,kk)                             &
       +jd(1,2,kk)*jd(2,3,kk)*jd(3,1,kk)                             &
       +jd(2,1,kk)*jd(3,2,kk)*jd(1,3,kk)                             &
       -jd(3,1,kk)*jd(2,2,kk)*jd(1,3,kk)                             &
       -jd(2,1,kk)*jd(1,2,kk)*jd(3,3,kk)                             &
       -jd(3,2,kk)*jd(2,3,kk)*jd(1,1,kk)

  sgn = 1.0
  DO i = 1, 3
    DO j = 1, 3
      ii = 1
      DO k = 1, 3
        jj = 1
        IF ( i /= k ) THEN
          DO l = 1, 3
            IF ( j /= l ) THEN

              cof(ii,jj) = jd(k,l,kk)
              jj       = jj + 1

            END IF
          END DO

          ii = ii + 1

        END IF
      END DO

      u(i,j) = sgn*( cof(1,1)*cof(2,2) - cof(1,2)*cof(2,1) )/det
      sgn = -sgn

    END DO
  END DO
  DO i = 1, 3
    DO j = 1, 3

      jd(i,j,kk) = u(j,i)

    END DO
  END DO
END DO


! Now solve J(dx) = F

dx = 0.0
dx(1) = jd(1,1,1)*f(1) + jd(1,2,1)*f(2) + jd(1,3,1)*f(3)
dx(2) = jd(2,1,1)*f(1) + jd(2,2,1)*f(2) + jd(2,3,1)*f(3)
dx(3) = jd(3,1,1)*f(1) + jd(3,2,1)*f(2) + jd(3,3,1)*f(3)
k     = 4
DO i = 2, n+1
  dt(1) = f(k)  -jl(1,1,i)*dx(k-3)-jl(1,2,i)*dx(k-2)-jl(1,3,i)*dx(k-1)
  dt(2) = f(k+1)-jl(2,1,i)*dx(k-3)-jl(2,2,i)*dx(k-2)-jl(2,3,i)*dx(k-1)
  dt(3) = f(k+2)-jl(3,1,i)*dx(k-3)-jl(3,2,i)*dx(k-2)-jl(3,3,i)*dx(k-1)

  dx(k) = jd(1,1,i)*dt(1) + jd(1,2,i)*dt(2) + jd(1,3,i)*dt(3)
  k     = k + 1
  dx(k) = jd(2,1,i)*dt(1) + jd(2,2,i)*dt(2) + jd(2,3,i)*dt(3)
  k     = k + 1
  dx(k) = jd(3,1,i)*dt(1) + jd(3,2,i)*dt(2) + jd(3,3,i)*dt(3)
  k     = k + 1
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE invert_jacob
END MODULE rcf_ideal_invert_jacob_mod
