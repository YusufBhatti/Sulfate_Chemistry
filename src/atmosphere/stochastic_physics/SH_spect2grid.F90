! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  This routine transforms the spherical harmonic
!               coefficients of the Forcing pattern into the
!               pattern grid-space global field, it calls
!               "fourier" function
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics

MODULE SH_spect2grid_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SH_SPECT2GRID_MOD'

CONTAINS

SUBROUTINE SH_spect2grid( coeffc, coeffs, stph_n2, nlim, coeff,         &
                          mu, nlat, first_atmstep_call, ii)
!
! Uses spherical harmonic coefficients that have been evolved in
! time using a Markov process to generate a Fourier series representation
! on each latitude circle, then inverse FFT to get gridpoint values
!

USE c_skeb2_mod, ONLY: nblock, nfacts
USE fourier_mod, ONLY: fourier

USE stochastic_physics_run_mod, ONLY: Ymn

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

LOGICAL :: first_atmstep_call
INTEGER :: stph_n2, nlim,lev
REAL    :: coeffc(0:stph_n2,1:stph_n2), coeffs(0:stph_n2,1:stph_n2)
REAL    :: coeff(0:2*nlim+1)
INTEGER :: m,m1,n,nn,mu
! Latitude pointers
INTEGER :: nlat,ii
REAL,ALLOCATABLE,SAVE :: seqalf(:)
! N+2 equivalent to  2*nlim+2
REAL    :: work((2*nlim+2)*nblock),trigs(2*nlim)
INTEGER :: ifax(nfacts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SH_SPECT2GRID'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate variables (saved for subsequent calls)
IF (.NOT. ALLOCATED(seqalf)) THEN
  ALLOCATE(seqalf(nlim*nlat*stph_n2))
END IF
ifax=0
coeff(:)= 0.0

DO m1 = 0, 2*nlim, 2
  m  = m1/2
  nn = MAX0(1,m)
  DO n = nn, stph_n2
    IF (first_atmstep_call) THEN
      seqalf(ii)= Ymn(mu,n,m)
    END IF
     ! cosine coeffs.
    coeff(m1)=   coeff(m1)   + coeffc(m,n)*seqalf(ii)
    ! sine coeffs.
    coeff(m1+1)= coeff(m1+1) + coeffs(m,n)*seqalf(ii)
    ii=ii+1

  END DO
END DO

! Invert fourier series representation on each latitude
CALL fourier( coeff, 2*nlim+2, trigs, ifax, 1, 2*nlim, 2*nlim,          &
             1, 1, work)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE SH_spect2grid

END MODULE SH_spect2grid_mod
