! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   module pws_icing_pot_mod ---------------------------------------
!
!   Purpose: Calculates PWS icing potential, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_icing_pot_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_ICING_POT_MOD'

REAL , PARAMETER :: thickness = 10000.0 ! 100hPa slabs of atmosphere.

CONTAINS

SUBROUTINE pws_icing_pot(p_theta_levels,theta,exner_theta_levels,cf_bulk,q)

! Description
!-----------------------------------------------------------------------------!
!                 Relative Humidity based icing diagnostic                    !
!-----------------------------------------------------------------------------!
! Method:                                                                     !
! 1. Calculate saturated mixing ratio from P & T and then combine this with   !
!    specific humidity to get RH. Then determine RH on pressure levels.       !
! 2. Create mask where cloud is                                               !
!    present. Determine cloud fraction on pressure levels.                    !
! 3. Create -20 to 0C temperature mask. Convert on to pressure levels.        !
! 4. Multiply the three pressure level arrays to give RH where cloud is       !
!    present and temperature is between 0 and -20C                            !
!-----------------------------------------------------------------------------!

USE pws_diags_mod, ONLY: icing_pot_press, icing_pot_press_levs,              &
                         pws_icing_pot_diag, flag_icing_pot

USE lionplyr_mod, ONLY: lionplyr

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, tdims_l, pdims
USE level_heights_mod,     ONLY: r_theta_levels
USE planet_constants_mod,  ONLY: planet_radius, repsilon
USE conversions_mod,       ONLY: zerodegc
USE missing_data_mod,      ONLY: imdi, rmdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels

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

REAL, INTENT(IN)  :: cf_bulk(tdims_l%i_start:tdims_l%i_end,                  &
                             tdims_l%j_start:tdims_l%j_end,                  &
                             tdims_l%k_start:tdims_l%k_end)

REAL, INTENT(IN)  :: q      (tdims_l%i_start:tdims_l%i_end,                  &
                             tdims_l%j_start:tdims_l%j_end,                  &
                             tdims_l%k_start:tdims_l%k_end)


! Local variables
REAL    :: temp(row_length,rows,model_levels) ! temperature
REAL    :: SatMR(row_length,rows,model_levels) ! saturated mixing ratio
REAL    :: RH(row_length,rows,model_levels) ! relative humidity

! RH on requested pressure levels
REAL    :: RHpress(row_length,rows,icing_pot_press_levs)
! temperature on requested pressure levels
REAL    :: temp_press(row_length,rows,icing_pot_press_levs)
REAL    :: temp_mask(row_length,rows,icing_pot_press_levs)
! bulk cloud on requested pressure levels
REAL    :: cloud_press(row_length,rows,icing_pot_press_levs)
REAL    :: cloud_mask(row_length,rows,icing_pot_press_levs)

INTEGER :: i,j,k

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_ICING_POT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!   calculate the temperature.
!-------
!   Calculate the saturated mixing ratio (for use in relative humidity
!   calculation) given temperature and pressure as input.

!   x_{sat} is calculated as:
!   x_{sat} = ( 0.622 * e_{s} ) / ( p - e_{s} )
!   where e_{s} is the saturated vapour pressure (SVP):
!         e_{s} = 6.11 * 10^( 7.5*T / (237.7 + T) )
!   using the Murray formulation, with T in deg C and p in hPa

!   Temperature is converted into
!   degrees C and pressure into hPa, hence in the equation
!   Temp = temp - 273.15
!   p = p / 100.0

!   Note also the part of the equation (237.7+T) has been reduced to
!   T - 35.45
!   which is equivalent to 237.7 + T - 273.15
!-------
!   Calculate the RH given the specific humidity and saturated mixing ratio
!   RH = x/x_{sat} = (SH/(1-SH))/x_{sat}
!   where RH = relative humidity, SH = specific humidity, and
!         x_{sat} = saturated mixing ratio


DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      temp(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
      SatMR(i,j,k) = ( repsilon *                                          &
                     ( 6.11 *                                              &
                     ( 10.0**( (7.5 * (temp(i,j,k) - zerodegc)) &
                          / (temp(i,j,k) - 35.45) ) )  )                   &
    !
                            /                                              &
    !
                     ( (p_theta_levels(i,j,k) / 100.0) -                   &
                     ( 6.11 * &
                     ( 10.0**( (7.5 * (temp(i,j,k) - zerodegc)) &
                       / (temp(i,j,k) - 35.45) ) ) ) ) )
      RH(i,j,k) = (q(i,j,k) / (1.0 - q(i,j,k))) / SatMR(i,j,k)
    END DO
  END DO
END DO 


! Calculate relative humidity on pressure levels
! using method as previously employed in fieldcalc.

CALL lionplyr( icing_pot_press_levs,                            &     
                     icing_pot_press,                           &  ! in
                     thickness,                                 &  ! in
                     RH,                                        &  ! in
                     p_theta_levels ,                           &  ! in
                     RHpress)          

! Use MIN to make sure relative humidity is in range 0->1 otherwise we get
! unrealistic ICING potential.

DO k=1, icing_pot_press_levs
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
       RHpress(i,j,k)=MIN(RHpress(i,j,k),1.0)
    END DO
  END DO
END DO

!----------------------------------------------------------------
!---- Step 2: Determine temperature on pressure levels        ---
!---- calculating mean temperature across atmospheric layer   ---
!----------------------------------------------------------------

CALL lionplyr( icing_pot_press_levs,                            &     
                     icing_pot_press,                           &  ! in
                     thickness,                                 &  ! in
                     temp,                                      &  ! in
                     p_theta_levels ,                           &  ! in
                     temp_press)    

! Calculate which gridboxes have temperatures between -20 and 0C    
DO k = 1, icing_pot_press_levs
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      IF ( (temp_press(i,j,k) > (zerodegc-20.0))     &
                 .AND. (temp_press(i,j,k) < zerodegc ) )THEN

        temp_mask(i,j,k) = 1.0
      ELSE
        temp_mask(i,j, k) = 0.0
      END IF

    END DO
  END DO
END DO

!----------------------------------------------------------------
!---- Step 3: Determine cloud fraction on pressure levels    ----
!----------------------------------------------------------------
    
CALL lionplyr(icing_pot_press_levs ,                            &     
                     icing_pot_press,                           &  ! in
                     thickness,                                 &  ! in
                     cf_bulk(tdims%i_start:tdims%i_end,         &
                           tdims%j_start:tdims%j_end,           &
                           1:tdims%k_end),                      &  ! in
                     p_theta_levels,                            &  ! in
                     cloud_press)    

! Calculate which gridboxes have cloud present

DO k = 1,icing_pot_press_levs
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      IF ( ( cloud_press(i,j,k) > 0.0)   &
                 .AND. (cloud_press(i,j,k) < 1.1 ) ) THEN

        cloud_mask(i,j,k) = 1.0
      ELSE
        cloud_mask(i,j,k) = 0.0
      END IF

    END DO
  END DO
END DO

!-----------------------------------------------------------------
!---- Step 4: Apply temperature and cloud masks to             ---
!----         relative humidity                                ---
!-----------------------------------------------------------------

DO k = 1, icing_pot_press_levs
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      IF (RHpress(i,j,k) /= rmdi) THEN
        pws_icing_pot_diag(i,j,k)=RHpress(i,j,k)*cloud_mask(i,j,k) &
                                                *temp_mask(i,j,k)
      ELSE 
        pws_icing_pot_diag(i,j,k)=rmdi
      END IF

    END DO
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_icing_pot

END MODULE pws_icing_pot_mod
