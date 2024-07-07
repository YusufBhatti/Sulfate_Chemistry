! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_contrail_mod ---------------------------------------
!
!   Purpose: Calculates contrail diagnostics, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_contrail_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_CONTRAIL_MOD'

CONTAINS

SUBROUTINE pws_contrail(theta, p_theta_levels,exner_theta_levels,pstar)

! Description
!   This subroutine uses pressure, temperature and height to
!   calculate an upper and lower limit within which contrails are
!   expected to form.
!
! Method:
!   The aim is to find the intersection of the temperature profile and
!   the Mintra -14 deg 'environment curve' at each point.  At the
!   the intersection: Pressure (p_m) = (( Temperature)/A0)**(1/A1)
!   Also,             Temperature    = t(k)*(p_m/p(k))**(R*lapse/g)
!   So  p_m = ( (t(k)/A0)*p(k)**(-R*lapse/g) )**(1/(A1-R*lapse/g))
!   The routine loops through the input temperature fields to find the
!   2 intervals where the curve is crossed.  Interpolation is done
!   within these intervals to find the 2 intersections with the curve,
!   which are the upper and lower contrail limits.
!   Some information on Mintra lines can be found in the Forecasters
!   Reference Book.


USE pws_diags_mod, ONLY: pws_contrail_bot,    &
                         pws_contrail_top,    &
                         flag_contrail_bot, &
                         flag_contrail_top

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s

USE planet_constants_mod, ONLY: g_over_r,planet_radius
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels
USE level_heights_mod, ONLY: r_theta_levels

USE icao_ht_fc_mod, ONLY:ICAOHeight_FC

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


REAL, INTENT(IN)  :: theta(tdims_s%i_start:tdims_s%i_end,                    &
                           tdims_s%j_start:tdims_s%j_end,                    &
                           tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
                                        tdims_s%j_start:tdims_s%j_end,       &
                                        tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: pstar(tdims%i_start:tdims%i_end,           &
                           tdims%j_start:tdims%j_end)

REAL, PARAMETER :: a0 = 137.14816 ! Constants for approximating the
REAL, PARAMETER :: a1 = 0.046822  ! Mintra - 14 deg "environment curve"
REAL, PARAMETER :: heightcut_top = 18000.0 !  arbritary limits for high
REAL, PARAMETER :: low_pressure_cut = 5000.0 ! arbritary minimum press for test

REAL    :: temp(row_length,rows,model_levels) ! temperature

REAL    :: Plwr(tdims%i_start:tdims%i_end,           &
                tdims%j_start:tdims%j_end) ! Max press for contrail formation
REAL    :: Pupr(tdims%i_start:tdims%i_end,           &
                tdims%j_start:tdims%j_end) ! Min press for contrail formation

REAL :: mintra_k  (row_length,rows) ! Mintra line
REAL :: mintra_km1(row_length,rows) !   - 14 degrees
REAL :: p_m       (row_length,rows) ! P at intersection
REAL :: lapse               ! Lapse rate in layer containing
                            ! Intersection

REAL :: work , expon      ! temporary scalars

INTEGER :: LevLwr(row_length,rows)
INTEGER :: LevUpr(row_length,rows)

INTEGER :: i, j, k, lev     ! Loop counters

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_CONTRAIL'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL Timer( RoutineName, 3 )

! calculate the temperature.
DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      temp(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
    END DO 
  END DO 
END DO 


! Initialise level and pressure arrays
LevLwr(:,:) = 0
LevUpr(:,:) = 0
Plwr(:,:)   = 0.0
Pupr(:,:)   = 0.0

! Check lowest level for contrails 

DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    Mintra_k(i,j) = a0*(p_theta_levels(i,j,1) **a1) - temp(i,j,1)
    IF ( Mintra_k(i,j) > 0 ) THEN
      LevLwr(i,j) = 1              ! Contrails forming at surface
      PLwr(i,j) = pstar(i,j)  
    END IF
  END DO
END DO
  
! Find the 2 intervals containing the Mintra -14deg curve
DO k = 2,model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      Mintra_km1(i,j) = Mintra_k(i,j)
      Mintra_k(i,j) = a0*(p_theta_levels(i,j,k)**a1) - temp(i,j,k)

      IF ( (p_theta_levels(i,j,k) > low_pressure_cut)               .AND. &
         (r_theta_levels(i,j,k)-planet_radius < heightcut_top )     .AND. &
         (Mintra_k(i,j) * Mintra_km1(i,j) <= 0.0) ) THEN

         IF (LevLwr(i,j) > 0)  THEN
           LevUpr(i,j) = k   ! Top of contrail layer (bottom already found)
         ELSE IF (LevLwr(i,j) == 0)     THEN     
           LevLwr(i,j) = k            ! Bottom of contrail layer
         END IF
      END IF

    END DO
  END DO
END DO

! Interpolate within the intervals to find the curve exactly
DO lev = 1,2              ! Lwr then Upr
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end

      k = LevLwr(i,j)
      IF ( k > 1 ) THEN
        lapse = (temp(i,j,(k-1)) - temp(i,j,k)) &
                / (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1) )
        work  = (p_theta_levels(i,j,k) **lapse) &
                              * ((a0/temp(i,j,k))**G_over_R)
        expon = 1.0/(lapse - a1*G_over_R)
        p_m(i,j) = work**expon

      END IF

      IF ( lev == 1 ) THEN
        IF ( LevLwr(i,j) > 1 ) THEN
          PLwr(i,j) = p_m(i,j)
        END IF
      ELSE
        IF ( LevUpr(i,j) > 1 ) THEN
          PUpr(i,j) = p_m(i,j)
        END IF

      END IF
    LevLwr(i,j) = LevUpr(i,j)
    END DO
  END DO

END DO

! Convert pressure to ICAO height
CALL ICAOHeight_FC( PLwr, pws_contrail_bot, row_length, rows )
CALL ICAOHeight_FC( PUpr, pws_contrail_top, row_length, rows )

CALL Timer( RoutineName, 4 )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_contrail


END MODULE pws_contrail_mod
