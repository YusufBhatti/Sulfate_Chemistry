! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE backscatter_spectrum_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BACKSCATTER_SPECTRUM_MOD'

CONTAINS


SUBROUTINE backscatter_spectrum( gspect, nlim, stph_n1, stph_n2,        &
                                 timestep, tot_backscat, alpha)

USE planet_constants_mod, ONLY: planet_radius

USE conversions_mod, ONLY: pi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
IMPLICIT NONE

  ! This include contains the model PrintStatus

INTEGER,INTENT(IN) :: nlim, stph_n1, stph_n2
REAL, INTENT(IN)   :: timestep, tot_backscat, alpha
REAL, INTENT(OUT)  :: gspect(nlim)

!     Local variables
REAL    :: gamman, p
INTEGER :: n        ! Loop counter

!     Functions
REAL :: chi         ! function

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BACKSCATTER_SPECTRUM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This power law is what Glenn Shutts observed in the CRM simulations
! (see Berner et al., 2009: J. Atmos. Sci, pp 603-626)
! Power of -1.54 based on practice at ECMWF
! e.g.: chi=n^^(2p+1); p = -1.27 => 2p+1 = -1.54
p=-1.27

gamman= 0.0
DO n = stph_n1, stph_n2
  ! DEPENDS ON: chi
  gamman = gamman + (n+1)*(2*n+1)*chi(n)
END DO
gamman = gamman/alpha

DO n = 1, nlim
  gspect(n) = 0.0
END DO

! g(n) = b * n^p
! b = SQRT([4 * Pi * a^2 * tot_backscat]/[timestep * vz * GammaN]
! vz = variance of random numbers [-0.5; 0.5] = 1/12
!      tested by 1000 cases of random arrays of size = 1e9
! note: (4/vz) is pre-calculated = 48
DO n = stph_n1, stph_n2
  gspect(n) = planet_radius * SQRT(48.0*pi*tot_backscat/                &
             (timestep*gamman)) * n**p
END DO

IF (printstatus  >  prstatus_normal) THEN
  WRITE(umMessage,'(A,I5,A,2I5)') 'N (nlim=', nlim,                     &
       ' varies from stph_n1 to stph_n2 =', stph_n1, stph_n2
  CALL umPrint(umMessage,src='backscatter_spectrum')
  CALL umPrint("GSPECT: wavenumber-dependent power-spectrum amplitude", &
      src='backscatter_spectrum')
  DO n=1,nlim
    WRITE(umMessage,'(1ES12.4)') gspect(n)
    CALL umPrint(umMessage,src='backscatter_spectrum')
  END DO
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE backscatter_spectrum
END MODULE backscatter_spectrum_mod
