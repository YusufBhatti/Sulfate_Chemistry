! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Module pws_isotherm_mod ---------------------------------------
MODULE pws_isotherm_mod

!   Purpose: Calculates isotherm diags for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

USE pws_diags_mod, ONLY: pws_freezing_ht, pws_freezing_press,                &
                         pws_freezing_icao

USE icao_ht_fc_mod, ONLY:ICAOHeight_FC

IMPLICIT NONE

CONTAINS

SUBROUTINE pws_isotherm(tempref, orog, pstar, theta, p_theta_levels,         &
                         exner_theta_levels ,flag_icao)

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s

USE water_constants_mod, ONLY: tm  ! freezing/melting temp of water.
USE planet_constants_mod, ONLY: g_over_r, planet_radius
USE missing_data_mod, ONLY: rmdi
USE level_heights_mod, ONLY: r_theta_levels
USE nlsizes_namelist_mod, ONLY:   row_length, rows, model_levels

IMPLICIT NONE
REAL, INTENT(IN)  :: tempref           ! ref temp to search for.
REAL, INTENT(IN)  :: orog(row_length,rows)
REAL, INTENT(IN)  :: pstar(row_length,rows)
REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: theta(tdims_s%i_start:tdims_s%i_end,                    &
                           tdims_s%j_start:tdims_s%j_end,                    &
                           tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
                                        tdims_s%j_start:tdims_s%j_end,       &
                                        tdims_s%k_start:tdims_s%k_end)
LOGICAL, INTENT(IN) :: flag_icao ! are icao hts requested?

! Local Variables:
INTEGER :: i, j, k
INTEGER :: IsoLev(row_length,rows)
INTEGER :: icode
REAL    :: lapse_rate        ! lapse rate of layer containing TempRef
REAL    :: temp(row_length,rows,model_levels)

! For dry-bulb freezing level height correction below model orography
REAL, PARAMETER :: dblr = 10.0e-3  ! Dry-bulb lapse rate
! For dry-bulb freezing level pressure correction below model orography
REAL, PARAMETER :: simple_dpdz = -10.0 !
! Cut off height above which no search is required for sub tempref temps.
REAL, PARAMETER :: heightcut = 23000.0 ! 23km.

icode = 0
! calc local temperature.
DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      temp(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k

IsoLev(:,:) = 0

IF (tempref == tm) THEN  ! if looking for 273.15 freezing level

  WHERE ( temp(:,:,1) < tempref )
    IsoLev(:,:) = 1
    ! Approximate correction to height of freezing level below model orography.
    ! The calculation is designed to allow the return of values below model
    ! orography, but is capped to not return anything less than sea level.
    ! Perhaps the cap need not be there.
    pws_freezing_ht(:,:) = MAX(orog(:,:) + (temp(:,:,1) - tempref)/dblr,0.0)
    ! Approx correction to pressure of freezing level below model orography.
    ! NB. This one is not being capped.
    pws_freezing_press(:,:) = pstar(:,:) +                                   &
                             (simple_dpdz* (temp(:,:,1) - tempref)/dblr)
  END WHERE

ELSE  ! where temp is less than supplied threshold set 
      ! pressure and height to surf values.

  WHERE ( temp(:,:,1) < tempref )
    IsoLev(:,:) = 1
    pws_freezing_ht(:,:) = orog(:,:)
    pws_freezing_press(:,:) = pstar(:,:)
  END WHERE
END IF


! Identify the kth model level that is less than threshold and the
! level below was above the threshold. 
! This level must be less than 23km so to avoid odd inversions.
DO k = 2, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF (temp(i,j,k) < tempref .AND.                                        &
            (r_theta_levels(i,j,k)-planet_radius) < heightcut) THEN
        IF (temp(i,j,k-1) > tempref ) THEN
          IsoLev(i,j) = k
        END IF
      END IF
    END DO
  END DO
END DO

! calc the isotherm ht and pressure.

DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    k = IsoLev(i,j)
    IF ( k >= 2 ) THEN

      lapse_rate = ( temp(i,j,k-1) - temp(i,j,k) )                           &
                  / ( r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1) )
      pws_freezing_ht(i,j) = (r_theta_levels(i,j,k-1)-planet_radius) +       &
                      (temp(i,j,k-1)-tempref)/lapse_rate
      pws_freezing_press(i,j) = p_theta_levels(i,j,k-1) *                    &
              (tempref/temp(i,j,k-1))**(g_over_r/lapse_rate)
    END IF
  END DO
END DO

! where no level has been identified need to set output diag to mdi

WHERE ( IsoLev(:,:) == 0 )
  pws_freezing_ht(:,:) = rmdi
  pws_freezing_press(:,:) = rmdi
END WHERE

! if icao ht requested calc from the presuure level.
IF (flag_icao) THEN
  CALL ICAOHeight_FC(pws_freezing_press, pws_freezing_icao,                 &
                      tdims%i_len,tdims%j_len )
END IF


END SUBROUTINE pws_isotherm

END MODULE pws_isotherm_mod
