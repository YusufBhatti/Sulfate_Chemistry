! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_dust_diag ---------------------------------------
!
!   Purpose: Calculates dust amounts for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No ??
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_dust_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_DUST_DIAG_MOD'

! Layer dimensions in metres.
REAL, PARAMETER :: top_of_layer = 1524.0
REAL, PARAMETER :: bottom_of_layer = 609.6

CONTAINS

SUBROUTINE pws_dust_diag(pws_dustconc_tot, pws_dustconc_5000)

USE level_heights_mod, ONLY: r_rho_levels
USE planet_constants_mod, ONLY: planet_radius

USE nlsizes_namelist_mod, ONLY: model_levels
USE atm_fields_bounds_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: pws_dustconc_tot  (tdims%i_start:tdims%i_end,  &
                                       tdims%j_start:tdims%j_end,  &
                                       1:tdims%k_end)

REAL, INTENT(INOUT) :: pws_dustconc_5000 (tdims%i_start:tdims%i_end,  &
                                          tdims%j_start:tdims%j_end)

! Local variables
INTEGER :: i, j, k
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_DUST_DIAG'

REAL, PARAMETER :: layer_thickness = top_of_layer - bottom_of_layer

! thickness of dust layer
REAL :: thickness (tdims%i_len, tdims%j_len)

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thickness         (i,j) = 0.0
      pws_dustconc_5000 (i,j) = 0.0
    END DO
END DO

! compute dust conc for 2000/5000 ft range from 
! aerosol dust conc diagnostic stash code 17,257 (units microgm/m3)

! Note that r_rho_levels(:,:,k+1), r_rho_levels(:,:,k) are the upper and lower
! layer boundaries for r_theta_levels(:,:,k)
! r_rho_levels is of dimensions tdims_l

DO k = 2, model_levels-1
  
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF ( (r_rho_levels(i,j,k+1)                               &
          - r_rho_levels(i,j,1)) < bottom_of_layer ) THEN
        thickness (i,j) = 0.0
      ELSE IF ( (r_rho_levels(i,j,k)                            &
          - r_rho_levels(i,j,1)) > top_of_layer) THEN
        thickness (i,j) = 0.0
      ELSE
        thickness (i,j) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      END IF
    END DO
  END DO
  
  
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      pws_dustconc_5000(i,j) = pws_dustconc_5000(i,j) + &
                                  thickness(i,j) * pws_dustconc_tot(i,j,k)
    END DO
  END DO

END DO

! convert from microg/m3 to g/m3 and scale by layer_thickness
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    pws_dustconc_5000(i,j) = 1.0e-6 * pws_dustconc_5000(i,j) / layer_thickness

  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_dust_diag

END MODULE pws_dust_diag_mod
