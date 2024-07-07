! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE zenith_delay_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ZENITH_DELAY_MOD'

CONTAINS

SUBROUTINE zenith_delay(zen_tot_delay)
!
! Description: Routine to calculate Zenith Total Delay
!
! Method: 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod
USE atm_fields_real_mod, ONLY: orography, exner_theta_levels, &
                               theta, p_theta_levels, q ,&
                               exner_rho_levels  ! 1:model_levels+1

USE level_heights_mod, ONLY: r_rho_levels, r_theta_levels

USE planet_constants_mod, ONLY: planet_radius, g, r, cp, pref,            &
                                c_virtual, repsilon, grcp
USE nlsizes_namelist_mod, ONLY: model_levels, row_length, rows

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Subroutine Arguments:
REAL, INTENT(OUT) :: zen_tot_delay &
             (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end) 

! Local Constants:
REAL, PARAMETER :: nalpha=77.6
REAL, PARAMETER :: nbeta=3.73e5

! Local variables:
INTEGER           :: i, j, k
REAL              :: temp1, temp2
REAL              :: exner_on_theta
REAL              :: t
REAL              :: Tv
REAL              :: Nwet
REAL              :: Ndry
REAL              :: c
REAL, ALLOCATABLE :: Refrac(:,:)
REAL, ALLOCATABLE :: Refrac_kp1(:,:)
REAL, ALLOCATABLE :: h_rho_levels(:,:,:)
REAL, ALLOCATABLE :: h_theta_levels(:,:,:)
REAL              :: h_diff

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ZENITH_DELAY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(refrac    (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
ALLOCATE(refrac_kp1(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))

ALLOCATE(h_rho_levels(pdims%i_start:pdims%i_end,  &
                      pdims%j_start:pdims%j_end,  &
                      pdims%k_start:pdims%k_end+1))  ! 1:model_levels+1

ALLOCATE(h_theta_levels(tdims%i_start:tdims%i_end,  &
                        tdims%j_start:tdims%j_end,  &
                                    1:tdims%k_end))  ! 1:model_levels

exner_on_theta           =  0.0
t                        =  0.0
Tv                       =  0.0
Nwet                     =  0.0
Ndry                     =  0.0
refrac                   =  0.0
refrac_kp1               =  0.0
h_diff                   =  0.0

! Set up levels heights for theta and rho.
DO k = pdims%k_start, pdims%k_end
  DO i = pdims%i_start, pdims%i_end
    DO j = pdims%j_start, pdims%j_end
      h_rho_levels(i,j,k) = r_rho_levels(i,j,k) - planet_radius
    END DO
  END DO
END DO

DO k = 1, tdims%k_end
  DO i = tdims%i_start, tdims%i_end
    DO j = tdims%j_start, tdims%j_end
      h_theta_levels(i,j,k) = r_theta_levels(i,j,k) - planet_radius
    END DO
  END DO
END DO

! r_rho_levels extends only to model_levels, but we require model_levels+1,
! so we insert top level (83933.328m for 70 level model). Note that the top 
! two p levels (i.e. model_levels and model_levels+1) are equally spaced
! above and below the top theta level.
DO i = pdims%i_start, pdims%i_end
  DO j = pdims%j_start, pdims%j_end
    h_diff = h_theta_levels(i,j,tdims%k_end) - h_rho_levels(i,j,pdims%k_end)
    h_rho_levels(i,j,pdims%k_end+1) = h_theta_levels(i,j,tdims%k_end) + h_diff
  END DO
END DO

! Calculate the correction for delay above top of model
DO i = tdims%i_start, tdims%i_end
  DO j = tdims%j_start, tdims%j_end
    zen_tot_delay(i,j) =  1.0e-6 * Nalpha /                                 &
                            100.0 *r *(p_theta_levels(i,j,tdims%k_end)/g)
  END DO
END DO

DO k = model_levels, 1, -1
  DO j= tdims%j_start, tdims%j_end
    DO i= tdims%i_start, tdims%i_end

      Nwet  = 0.0

      ! Compute mean layer virtual temperature
      temp1 = grcp *(h_rho_levels(i,j,k+1) - h_rho_levels(i,j,k))
      temp2 = exner_rho_levels(i,j,k) - exner_rho_levels(i,j,k+1)
      Tv    = exner_theta_levels(i,j,k) *temp1/temp2

      ! Compute wet refractivity
      t = Tv / (1.0 + C_virtual *q(i,j,k))
      temp1 = 0.01 *nbeta *p_theta_levels(i,j,k) *q(i,j,k)
      temp2 = ( t*t *(repsilon + (1-repsilon) *q(i,j,k)) )
      Nwet  = temp1/temp2

      ! Compute dry refractivity for this level
      Ndry = 0.01 *Nalpha *p_theta_levels(i,j,k) / t

      ! Total refractivity
      Refrac_kp1(i,j) = Refrac(i,j)
      Refrac(i,j)     = Ndry + Nwet      

      IF (k == model_levels) CYCLE

      c  = (LOG (Refrac_kp1(i,j) / Refrac(i,j))) /                         &
                         (h_theta_levels(i,j,k) - h_theta_levels(i,j,k+1))

      ! Compute zenith delay for this level
      zen_tot_delay(i,j) = zen_tot_delay(i,j) +                            &
                             (-1.0e-6 * Refrac(i,j) *                      &
                             EXP(c *h_theta_levels(i,j,k))*                &
                            (EXP(-1.0 *c *h_theta_levels(i,j,k+1)) -       &
                             EXP(-1.0 *c *h_theta_levels(i,j,k)))  / c) 

    END DO
  END DO
END DO

DEALLOCATE(refrac)
DEALLOCATE(refrac_kp1)
DEALLOCATE(h_rho_levels)
DEALLOCATE(h_theta_levels)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE zenith_delay

END MODULE zenith_delay_mod
