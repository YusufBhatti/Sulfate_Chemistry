! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Compute the spectral power factor for the Forcing Pattern,
!               It is based on a Gaussian power law and its amplitude.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE pattern_spectrum2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PATTERN_SPECTRUM2_MOD'

CONTAINS
SUBROUTINE pattern_spectrum2( pspect, nlim, alpha_stoch_tend)

USE atm_step_local,      ONLY: first_atmstep_call
USE stochastic_physics_run_mod, ONLY: stph_n2

USE planet_constants_mod, ONLY: planet_radius

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
IMPLICIT NONE

INTEGER, INTENT(IN)   :: nlim

REAL, INTENT(OUT)     :: pspect(nlim)
! function controlling the power spectrum of pattern def as g(n)
INTEGER               :: n

REAL :: Lbeta
! Correlation scale factor
REAL, SAVE  :: beta
! Modulator of the Gaussian amplitude.
REAL, INTENT(IN) :: alpha_stoch_tend
REAL :: tmp

REAL :: std_rand_nos
 ! square root of the standard deviation of random numbers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PATTERN_SPECTRUM2'
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! square root of Random number's sigma
! (sqrt to make the computation of F0 easier).
std_rand_nos= 1.0/SQRT(12.0)

DO n = 1, nlim
  pspect(n) = 0.0
END DO
!
! specify the unscaled functional form of the power spectrum
!
IF (first_atmstep_call) THEN
  Lbeta = 500e3 !Original value eq 500e3
  beta  = 0.5*(Lbeta/planet_radius)**2
END IF

DO n = 1, stph_n2
  pspect(n) = EXP(-beta*n*(n+1))
END DO

 ! To compute the amplitude of the power Law (equivalent to SKEB's F0),
 ! for constant tau_spt
tmp = 1.0/(2.0 - alpha_stoch_tend)
DO n = 1, stph_n2
  tmp = tmp + (2*n+1.0)*pspect(n)**2/(1.0 - 0.5*alpha_stoch_tend)
END DO
! Compute the power spectra chi function with the amplitude.
DO n = 1, stph_n2
  pspect(n) = pspect(n)/(std_rand_nos*SQRT(tmp))
END DO

IF (printstatus  >  prstatus_normal) THEN
  WRITE(umMessage,'(A,I0,A,I0)') 'N (nlim=', nlim,                      &
       ' varies from 1 to stph_n2 =', stph_n2
  CALL umPrint("PSPECT: wavenumber-dependent power-spectrum amplitude", &
  src='pattern_spectrum2')

  DO n = 1, nlim
    WRITE(umMessage,'(1ES12.4)') pspect(n)
    CALL umPrint(umMessage,src='pattern_spectrum2')
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE pattern_spectrum2
END MODULE pattern_spectrum2_mod
