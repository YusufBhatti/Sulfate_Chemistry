! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   solves tridiagonal matrix -  convection scheme
!
MODULE tridiag_all_n_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TRIDIAG_ALL_N_MOD'
CONTAINS

SUBROUTINE tridiag_all_n(n,nvec,a,b,c,r,u)
!
! Purpose: Solves the equations A.X = Y,  where A is a tridiagnol matrix
!
!          for several matrices A at once.
!          Version where all vetcors same length.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
INTEGER, INTENT(IN) :: &
  n                    & ! maximum size of vectors X and Y
 ,nvec                   ! Number of vectors/matrices to solve


REAL, INTENT(IN) ::    &
  a(nvec,n)            & ! Components of tridiagonal matrix
 ,b(nvec,n)            & ! Components of tridiagonal matrix
 ,c(nvec,n)            & ! Components of tridiagonal matrix
 ,r(nvec,n)              ! vector Y, R.H.S of linear equation

REAL, INTENT(OUT) ::   &
  u(nvec,n)              ! solution vectors


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER :: &
  i,j            ! loop counters

REAL ::         &
  gam(nvec,n)   & ! work array
 ,bet(nvec)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIDIAG_ALL_N'

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1,nvec
  bet(i) = b(i,1)
  u(i,1) = r(i,1)/bet(i)
END DO

DO j=2,n
  DO i=1, nvec

    gam(i,j) = c(i,j-1)/bet(i)
    bet(i)   = b(i,j) - a(i,j)*gam(i,j)
    u(i,j)   = (r(i,j) - a(i,j)*u(i,j-1))/bet(i)

  END DO
END DO

DO j=n-1,1,-1
  DO i=1, nvec

    u(i,j) = u(i,j) - Gam(i,j+1)*u(i,j+1)

  END DO
END DO

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE tridiag_all_n
END MODULE tridiag_all_n_mod
