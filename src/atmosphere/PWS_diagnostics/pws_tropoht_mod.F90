! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_tropoht_mod ---------------------------------------
!
!   Purpose: Calculates PWS tropopause diags, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_tropoht_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_TROPOHT_MOD'

CONTAINS

SUBROUTINE pws_tropoht(theta, p_theta_levels, exner_theta_levels)

! Description
! Find interval (k to k+1) in which tropopause lies by checking the lapse
! rate of surrounding intervals, and then calculate the corresponding
! height, pressure and temperature at the tropopause by using the lapse rates
! above and below the interval identified.

! Method:
! 1) Loop over model_levels; assuming that
!    1,2 and model_levels-1,model_levels will never be required as they
!    lie outside the height range search criteria for the tropopause. 
!    Thus we will not go out of bounds.
!
!    i)  Calculate lapse rate field for interval k~k+1
!        and the lapse rate for the interval below k-1~k
!
!    ii) Where the lapse rate of k~k+1 < 2C and the lapse rate of the interval
!         below is +ve  -- we may have found the tropopause interval -- search 
!         upwards for the next model level that at least 2km above k (k2km)
!         When (if found) check that the lapse rate k~k2km < 2C.
!         If it does then k is the bottom model level of the WMO tropopause.  
!         Otherwise we cycle around again looking at the next k level.   
!
! 2) Having found the model level k that is close to the bottom of the 
!    tropopause, we shall then use the lapse rates of the intervals above and
!    below the k level to calculate the tropopause height at the intersection
!    of these lapse rate lines.     
!
!       Equation of first line: 
!         in form y=mx+c where m=-lapselwr and c = templwr + lapselwr*heightlwr
!       i)  TropTemp = - lapselwr*Tropheight + (templwr + lapselwr*heightlwr)
!       Equation of second line:
!         in form y=mx+c where m=-LapseUpr and c = TempUpr + LapseUpr*HeightUpr
!       ii) TropTemp = - lapseupr*Tropheight + (tempupr + lapseupr*heightupr)

!       Point of intersection is given by:
!         - lapseupr*tropheight + (tempupr + lapseupr*heightupr) =
!                         - lapselwr*tropheight + (templwr + lapselwr*heightlwr)
!       => (lapselwr - lapseypr)*tropheight =
!               (templwr + lapselwr*heightlwr) - (tempupr + lapseupr*heightupr)
!       => TropHeight = (templwr + lapselwr*heightlwr) -
!               (tempupr + lapseupr*heightupr)/(lapselwr - lapseupr)
!
!     iii) Calculate tropopause temperature by substuting trop height into
!          equation of first line:
!          Troptemp = templwr - lapselwr * (Tropheight - heightlwr)
!
!     iv) Calculate tropopause pressure:
!          TropPress = Presslwr * (TropTemp/templwr)**(g_over_r/lapselwr)
!

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s
USE level_heights_mod,     ONLY: r_theta_levels
USE planet_constants_mod,  ONLY: planet_radius, g_over_r, lapse_trop
USE missing_data_mod,      ONLY: imdi, rmdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels

USE pws_diags_mod, ONLY: pws_tropopause_ht,    &
                         pws_tropopause_temp,  &
                         pws_tropopause_press, &
                         pws_tropopause_icao,  &
                         flag_tropopause_icao, &
                         heightcut_top,        &
                         heightcut_bot,        &
                         tempcut

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

! Local variables

INTEGER :: i, j, k, k2km
REAL, PARAMETER :: vsmall = 1.0e-6
REAL    :: temp(row_length,rows,model_levels) ! temperature

INTEGER :: tlev(row_length, rows)
INTEGER :: missing_count
! Lapse rates below and above current layer
REAL :: lapseupr, lapselwr
! Lapse rates below, at and also 2km layer
REAL :: lapse_below(row_length, rows)
REAL :: lapse(row_length, rows)
REAL :: lapse_2km(row_length, rows)
REAL :: delta_lapse        ! lapse_lwr-lapse_upr

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_TROPOHT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -----------------------------------
! calculate the temperature.
DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      temp(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
    END DO 
  END DO 
END DO 

! -----------------------------------

! Initialise tropopause model level to not found, imdi 
tlev (:,:) = imdi


! The following assumes that model levels 1-3 are lower than 4500m and thus
! we will always be able to find a k-2 level and not go out of bounds.
! Likewise model_levels is above 22km. 
! It may be worth adding a check for this and abort.

! For each model point between 4500m and 20km we look for a lapse
! rate that may signify we have found the nearest tropopause model level

DO k=1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ( tlev(i,j)==imdi                    .AND.                           &
        temp(i,j,k) < tempcut .AND.                                          &  
        (r_theta_levels(i,j,k)-planet_radius >  heightcut_bot)  .AND.        &
        (r_theta_levels(i,j,k)-planet_radius <  heightcut_top) ) THEN

         ! lapse rate for k~k+1 interval
         Lapse(i,j)   = (temp(i,j,k) - temp(i,j,k+1)) /                       &
                         (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))

         ! lapse rate for k-1~k interval
         Lapse_below(i,j)  = (temp(i,j,k-1) - temp(i,j,k)) /                  &
                         (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))

         ! lapse rate is < 2C we may have found tropopause interval.
         IF ( (Lapse(i,j) < lapse_trop)    .AND.                              &
              (Lapse_below(i,j) > 0.0)   )      THEN

           ! now look for lapse rate for corresponding 2km interval.
           DO k2km=k,model_levels
             IF (r_theta_levels(i,j,k2km)-planet_radius > heightcut_top) EXIT  
             IF ( (r_theta_levels(i,j,k2km)-r_theta_levels(i,j,k)) >=        &
                                           2000.0 ) THEN
             
               Lapse_2km(i,j)  = (temp(i,j,k) - temp(i,j,k2km)) /            &
                         (r_theta_levels(i,j,k2km) - r_theta_levels(i,j,k))

               ! if 2km interval also < 2 then we have the tropopause level
               IF (Lapse_2km(i,j) < lapse_trop) THEN
                 tlev(i,j) = k 
               END IF

               EXIT

             END IF
                 
           END DO

         END IF
         
      END IF
    END DO
  END DO
END DO


! Assuming the tropopause was found in the steps above, retrieve the
! level index for it here...
! We shall calculate the lapse rates above and below expected level
! and assume linear cross over point is nearer to the true height of the
! tropopause level -- it need not correspond to an exact model level --
! if no model level was identified return mdi

DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end

    ! if level found.
    IF (tlev(i,j) /= imdi) THEN
 
      k = tlev(i,j) 
      ! lapse rate for interval above
      lapseupr = (temp(i,j,k+1) - temp(i,j,k+2)) /                          &
                 (r_theta_levels(i,j,k+2) - r_theta_levels(i,j,k+1))
      ! lapse rate for interval below
      lapselwr = (temp(i,j,k-1) - temp(i,j,k)) /                            &
                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))


      delta_lapse = lapselwr - lapseupr
      IF ( ABS(delta_lapse) < vsmall ) THEN
        IF ( delta_lapse >= 0 ) delta_lapse =  vsmall
        IF ( delta_lapse <  0 ) delta_lapse = -vsmall
      END IF

      ! height of tropopause is
      pws_tropopause_ht(i,j) = ( (temp(i,j,k)+                               &
               (lapselwr*(r_theta_levels(i,j,k)-planet_radius))) -           &
                               (temp(i,j,k+1)+                               &
               (lapseupr*(r_theta_levels(i,j,k+1)-planet_radius))) )         &
                               / delta_lapse

      IF (pws_tropopause_ht(i,j) < (r_theta_levels(i,j,k)-planet_radius) ) THEN
        ! ensure trop level doesn't undershoot
        pws_tropopause_ht(i,j)  = r_theta_levels(i,j,k)-planet_radius 
      END IF
      IF (pws_tropopause_ht(i,j)> (r_theta_levels(i,j,k+1)-planet_radius)) THEN
        ! or overshoot
        pws_tropopause_ht(i,j)  = r_theta_levels(i,j,k+1)-planet_radius
      END IF

      ! Temperature at tropopause
      pws_tropopause_temp(i,j) = temp(i,j,k) - lapselwr *                    &
          (pws_tropopause_ht(i,j)-(r_theta_levels(i,j,k)-planet_radius))

      IF ( ABS(lapselwr) < vsmall ) THEN
        IF ( lapseLwr >= 0 ) lapselwr =  vsmall
        IF ( lapselwr <  0 ) lapselwr = -vsmall
      END IF

      ! P at tropopause is derived from hydrostatic equation.
      pws_tropopause_press(i,j) =  p_theta_levels(i,j,k)*                    &
                   (pws_tropopause_temp(i,j)/temp(i,j,k))**(g_over_r/lapselwr)
    ELSE
      pws_tropopause_press(i,j) = rmdi
      pws_tropopause_temp(i,j)  = rmdi
      pws_tropopause_ht(i,j)    = rmdi
      pws_tropopause_icao(i,j)  = rmdi
    END IF

  END DO
END DO

IF (flag_tropopause_icao) THEN
  CALL ICAOHeight_FC(pws_tropopause_press, pws_tropopause_icao,              &
                      tdims%i_len,tdims%j_len)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_tropoht
END MODULE pws_tropoht_mod
