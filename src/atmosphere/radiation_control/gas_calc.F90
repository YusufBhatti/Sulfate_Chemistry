! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE gas_calc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'GAS_CALC_MOD'
CONTAINS

SUBROUTINE gas_calc(gas_now                                       &
                   ,gas_index_max                                 &
                   ,gas_year                                      &
                   ,gas_conc                                      &
                   ,gas_rate                                      &
                   ,max_scenario_pts                              &
                   ,icode)

!   Subroutine GAS_CALC ----------------------------------------------

!   Purpose :
!     Calculates the trace gas mixing ratio (or weighting factor for
!   aerosol forcing fields.  Rates of increase (yearly compound factors)
!   can be supplied, or spot values (which will be linearly
!   interpolated) or a mixture of these. It is designed so it can be
!   called each time step, but when rates of increase are being used,
!   values are in fact only updated at New Year.
!   The rules are:
!     If rates exist (i.e. are positive) for the first & current years
!   then all concentrations are ignored, except for the initial value.
!     If there is a positive rate for the current year but not for the
!   start, the current rate & most recent concentration are used.
!     If rates do not exist for the current year then the concentration
!   is calculated by linear interpolation between the concentrations at
!   the given years.
!     The mixing ratios calculated after the last given year use the
!   rate for the final given year.
!   CARE should be taken if solitary rates are specified, as this can
!   result in discontinuities in the concentration time profile at
!   the next given year without a corresponding given rate.
!   ------------------------------------------------------------------

! Use temporary science fix to hardwire code to 360 day calendar
USE science_fixes_mod, ONLY: l_rm_hardwire_gas360

USE clmchfcg_scenario_mod, ONLY: L_Cts_Fcg_Rates
USE conversions_mod, ONLY: isec_per_day
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstcall_mod,  ONLY: lcal360
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE model_time_mod, ONLY: &
    i_day, i_hour, i_minute, i_month, i_second, i_year

IMPLICIT NONE

! Submodel parameters for array sizes

REAL, INTENT(OUT)    :: gas_now      ! Gas concentration at time step
INTEGER, INTENT(IN)  :: gas_index_max, max_scenario_pts
INTEGER, INTENT(IN)  :: gas_year(max_scenario_pts)
REAL, INTENT(IN)     :: gas_conc(max_scenario_pts)
REAL, INTENT(IN)     :: gas_rate(max_scenario_pts)
INTEGER, INTENT(OUT) :: icode        ! Return code: successful=0

! Local variables
LOGICAL :: l_loc_cal360   ! Local calendar logical for "science fix"

INTEGER :: gas_index      ! to subscript gas concs for NOW_TIME
INTEGER :: i              ! Loop over indices
INTEGER :: now_in_secs    ! Current time in seconds from current GAS_YEAR
INTEGER :: year_in_secs   ! Year length in seconds from current GAS_YEAR
INTEGER :: now_time_day, now_time_sec
                          ! Time now in days/secs from base time
INTEGER :: gas_yr_day1,  gas_yr_sec1
                          ! Time in days/secs of current GAS_YEAR
INTEGER :: gas_yr_day2,  gas_yr_sec2
                          ! Time in days/secs of next GAS_YEAR

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GAS_CALC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that GASCNST namelist is defined for this year
IF ( i_year < gas_year(1) ) THEN
  icode = 8325
  WRITE(umMessage, FMT='(A)') 'GAS_CALC: no gas data for this year'
  CALL umPrint(umMessage,src='gas_calc')
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set local logical for calendar to get unfixed behaviour where
! routine is hardwired for 360 day calendar. This should go with removal
! of the temporary logical and all further references of l_loc_cal360
! should be returned to lcal360.
IF (l_rm_hardwire_gas360) THEN
  l_loc_cal360 = lcal360
ELSE
  l_loc_cal360 = .TRUE.
END IF

! Loop over I to find correct index for current NOW_TIME
gas_index = 0
DO i=1, gas_index_max
  IF ( i_year >= gas_year(i) ) gas_index = gas_index+1
END DO

! Calculate base time of current GAS_YEAR
! DEPENDS ON: time2sec
CALL time2sec (gas_year(gas_index), 1, 1, 0, 0, 0,                            &
              0, 0, gas_yr_day1, gas_yr_sec1, l_loc_cal360)

! Calculate time now since base time in seconds
! DEPENDS ON: time2sec
CALL time2sec (i_year, i_month, i_day, i_hour, i_minute, i_second,            &
         gas_yr_day1, gas_yr_sec1, now_time_day, now_time_sec, l_loc_cal360)
now_in_secs = now_time_day * isec_per_day + now_time_sec

! If gas rate at current year is non zero calculate new GAS_NOW
! by considering compound increases of GAS_RATE(1:GAS_INDEX)
IF ( gas_rate(gas_index) > 0.0 ) THEN
  gas_now = gas_conc(1)
  ! Find last concentration, denoted by last year with negative gas_rate
  ! Then apply gas_rates prior to current year
  DO i=1, gas_index-1
    IF ( gas_rate(i) < 0.0 ) THEN
      gas_now = gas_conc(i+1)
    ELSE
      gas_now = gas_now * gas_rate(i) ** ( gas_year(i+1) - gas_year(i) )
    END IF
  END DO

  ! Calculate time of GAS_YEAR+1 from base time in seconds
  ! Note that we need to do this to get year_frac > 1
  ! This will be a limitation of gregorian calendar if going from or to
  ! a leap year.
! DEPENDS ON: time2sec
  CALL time2sec (gas_year(gas_index)+1, 1, 1, 0, 0, 0,                        &
           gas_yr_day1, gas_yr_sec1, gas_yr_day2, gas_yr_sec2, l_loc_cal360)
  year_in_secs = gas_yr_day2 * isec_per_day + gas_yr_sec2
  ! GAS_NOW now holds the concentration in year GAS_INDEX
  ! - need only update it to the current year.
  IF (L_Cts_Fcg_Rates) THEN
    ! Provides a continously varying gas concentration
    gas_now = gas_now * (gas_rate(gas_index) **                               &
              (REAL(now_in_secs) / REAL(year_in_secs)))
  ELSE
    ! Provides concentration that updates at beginning of year
    gas_now = gas_now * (gas_rate(gas_index) **                               &
              REAL(now_in_secs / year_in_secs))
  END IF

! Otherwise calculate by linear interpolation between respective
! GAS concentrations of given years.
ELSE
  ! Calculate time of next GAS_YEAR from base time in seconds
! DEPENDS ON: time2sec
  CALL time2sec (gas_year(gas_index+1), 1, 1, 0, 0, 0,                        &
           gas_yr_day1, gas_yr_sec1, gas_yr_day2, gas_yr_sec2, l_loc_cal360)
  year_in_secs = gas_yr_day2 * isec_per_day + gas_yr_sec2
  gas_now = gas_conc(gas_index) +                                             &
          ( gas_conc(gas_index+1) - gas_conc(gas_index) )                     &
          * REAL(now_in_secs) / REAL(year_in_secs)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE gas_calc
END MODULE gas_calc_mod
