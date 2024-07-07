! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE skeb_forcing_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SKEB_FORCING_MOD'

CONTAINS

SUBROUTINE skeb_forcing( dpsidtc, dpsidts, stph_n1, stph_n2, nlim,      &
                         dpsidt, ilat, nlat, seqalf)

! Uses cosine and sine Fourier coefficients that have been evolved in
! time using a Markov process to generate the Fourier series and
! convert this into grid space
! The Process uses some arcane functions ALF and EPS to create the
! Fourier series. This code is processed in parallel by distributing the
! latitude loop (in calling routine stph_skeb2) over multiple processors.

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE c_skeb2_mod, ONLY: nblock, nfacts
USE fourier_mod, ONLY: fourier
IMPLICIT NONE

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Factors nblock and nfacts needed for FFT arrays

INTEGER, INTENT(IN) ::                                                    &
     stph_n1                                                              &
,    stph_n2                                                              &
,    nlim                                                                 &
,    ilat                                                                 &
,    nlat

REAL, INTENT(IN) ::                                                       &
     dpsidtc(0:stph_n2,stph_n1:stph_n2)                                   & 
,    dpsidts(0:stph_n2,stph_n1:stph_n2)                                   &
,    seqalf(stph_n1:stph_n2,0:nlim,nlat)

REAL,INTENT(OUT) ::                                                       &
     dpsidt(0:2*nlim+1)

! Local variables
INTEGER :: i,m,m1,n,nn

! N+2 equivalent to  2*nlim+2
REAL    :: work((2*nlim+2)*nblock),trigs(2*nlim)
INTEGER :: ifax(nfacts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SKEB_FORCING'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


ifax=0
dpsidt(:)= 0.0

DO m1 = 0, 2*nlim, 2
  m  = m1/2
  nn = MAX(stph_n1,m)
  DO n = nn, stph_n2
    ! cosine coeffs.
    dpsidt(m1)=   dpsidt(m1)   + dpsidtc(m,n)*seqalf(n,m,ilat)
    ! sine coeffs.
    dpsidt(m1+1)= dpsidt(m1+1) + dpsidts(m,n)*seqalf(n,m,ilat)
  END DO
END DO

! Invert fourier series representation on each latitude
CALL fourier( dpsidt, 2*nlim+2, trigs, ifax, 1, 2*nlim, 2*nlim,         &
             1, 1, work)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE skeb_forcing

END MODULE skeb_forcing_mod
