! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to calculate solar flux, using average solar flux, the shape 
!  of the solar cycle, and the time series of the solar cycle. Data used
!  here are from Lean et al. (2005), as downloaded from
!  www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/index.html
!
! Method: A singular value decomposition has been performed on the spectrally
!  resolved solar irradiance data. The data is very well described by the
!  product of a spectrally varying term and a time series. The spectral 
!  composition of the solar cycle has been fed through the FAST-JX binning
!  algorithm. Here we compute the product of the spectral signature
!  and the solar cycle time series to infer the spectral variation.

! If i_ukca_solcyc is set to (1) the Lean et al data are used for the time 
! period for which they exist (1882-2004), before and after this time period 
! an average cycle is repeated. This average cycle is calculated from the last
! 5 solar cycles in the observations. The period of this repeated average cycle 
! is 10 years 8 months. 
! If i_ukca_solcyc is set to (2) this average cycle is used at all times.
!
! 
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------

MODULE ukca_solflux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SOLFLUX_MOD'

CONTAINS

SUBROUTINE ukca_solflux(timearray, lcal360, solcyc, x)

USE fastjx_data, ONLY : solcyc_av, solcyc_ts   
USE ukca_option_mod, ONLY : i_ukca_solcyc
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: timearray(7)

LOGICAL, INTENT(IN) :: lcal360  ! 360-day calendar indicator

REAL, INTENT(IN)    :: solcyc(:) ! spectral component of sol cyc var

REAL, INTENT(OUT)   :: x(:) ! fl or quanta mod modification


! Local variables

REAL :: realtime        ! time in years
REAL :: reftime         ! reference time for periodic cycle
REAL :: sol_ph          ! phase of periodic cycle
INTEGER :: obs_idx
INTEGER :: cyc_idx
REAL :: frac 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SOLFLUX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! calculate real time in years, and reference time - 
! Dec 2004 (the end of the Lean et al. obs)
IF (lcal360) THEN
  realtime = REAL(timearray(1)) + REAL(timearray(7)-1)/360.        &
           + REAL(timearray(4))/8640. + REAL(timearray(5))/518400.
  reftime = 2004. + 345.0/360.
ELSE
  realtime = REAL(timearray(1)) + REAL(timearray(7)-1)/365.24      &
           + REAL(timearray(4))/8765.76 + REAL(timearray(5))/525946.
  reftime  = 2004. + 349.5/365.24
END IF

! modify solar flux with solar cycle
IF ((i_ukca_solcyc==1) .AND. realtime>=1882. .AND. realtime<=2005.) THEN
  obs_idx = FLOOR((realtime - 1882.042)*12.) + 1
  frac = (realtime - 1882.042)*12. - (obs_idx-1)
  IF (obs_idx < 1) THEN
    ! if at start of obs, interpolate between end of avg. cycle and beginning
    ! of obs
    x = solcyc*((1. - frac)* solcyc_av(128) + frac * solcyc_ts(1))
  ELSE IF (obs_idx >= 1476) THEN
    ! if at end of obs, interpolate between end of obs and beginning
    ! of avg cyc
    x = solcyc*((1. - frac) * solcyc_ts(1476) + frac * solcyc_av(1))
  ELSE
    x = solcyc*((1. - frac)*solcyc_ts(obs_idx) + frac*solcyc_ts(obs_idx+1))  
  END IF
ELSE ! using periodic avg. cycle
  ! Calculate phase of cycle, phase 0 set to Dec 2004, 
  ! cycle period = 10y8m (128m)
  sol_ph = MODULO((realtime-reftime)*12.,128.)
  cyc_idx = FLOOR(sol_ph) + 1
  frac = sol_ph - (cyc_idx-1)
  IF (cyc_idx < 128) THEN  
    x = solcyc*((1. - frac) * solcyc_av(cyc_idx) + frac * solcyc_av(cyc_idx+1))
  ELSE
    x = solcyc*((1. - frac) * solcyc_av(cyc_idx) + frac * solcyc_av(1))
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE ukca_solflux

END MODULE ukca_solflux_mod
