! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Apply time variations to the solar constant and solar spectrum
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE solvar_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOLVAR_MOD'

CONTAINS

SUBROUTINE solvar(model_time, Sp, solar_constant)

USE def_spectrum,           ONLY: StrSpecData
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE umPrintMgr,             ONLY: umPrint, umMessage
USE nlstcall_mod,           ONLY: lcal360

IMPLICIT NONE


! Model time at start of timestep:
! year, month, day of month, hour, minute, second, day number in year
INTEGER, INTENT(IN) :: model_time(7)

! Spectral data:
TYPE (StrSpecData), INTENT(INOUT) :: Sp

! Total solar irradiance at 1AU for given time
REAL, INTENT(OUT) :: solar_constant

! Local variables
INTEGER :: i, j, sub_band, band, term, next_band
!   Loop variables
INTEGER :: i_time
!   Look-up table index for current time
INTEGER :: p_year, p_month, p_day, p_hour, p_minute, p_second, p_day_number
!   Equivalent model time for periodic repeat of the last solar cycle
INTEGER :: base_day, base_sec, repeat_day, repeat_sec
!   Variables for calculating repeat period in days and seconds

INTEGER :: errcode 

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'SOLVAR'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find current time in look-up table
i_time=0
DO i=1, Sp%Var%n_times
  IF ( ( Sp%Var%time(1, i) >  model_time(1) ) .OR. &
       ( Sp%Var%time(1, i) == model_time(1) .AND.  &
         Sp%Var%time(2, i) >  model_time(2) ) .OR. &
       ( Sp%Var%time(1, i) == model_time(1) .AND.  &
         Sp%Var%time(2, i) == model_time(2) .AND.  &
         Sp%Var%time(3, i) >  model_time(3) ) .OR. &
       ( Sp%Var%time(1, i) == model_time(1) .AND.  &
         Sp%Var%time(2, i) == model_time(2) .AND.  &
         Sp%Var%time(3, i) == model_time(3) .AND.  &
         Sp%Var%time(4, i) >  model_time(4)*3600 + &
                              model_time(5)*60 +   &
                              model_time(6) ) ) EXIT
  i_time=i
END DO

IF (i_time == 0) THEN
  cmessage = 'No solar variability data for model time.'
  errcode = 1
  CALL ereport(RoutineName, errcode, cmessage)
ELSE IF (i_time == Sp%Var%n_times .AND. Sp%Var%n_repeat_times > 1) THEN
  ! Repeat times periodically once look-up table ends
  ! Find base time for start of repeat period
  CALL time2sec (Sp%Var%time(1,Sp%Var%n_times - Sp%Var%n_repeat_times), &
                 Sp%Var%time(2,Sp%Var%n_times - Sp%Var%n_repeat_times), &
                 Sp%Var%time(3,Sp%Var%n_times - Sp%Var%n_repeat_times), 0, 0, &
                 Sp%Var%time(4,Sp%Var%n_times - Sp%Var%n_repeat_times), &
                 0, 0, base_day, base_sec, lcal360)
  ! Find length of repeat period in days and seconds
  CALL time2sec (Sp%Var%time(1,Sp%Var%n_times), &
                 Sp%Var%time(2,Sp%Var%n_times), &
                 Sp%Var%time(3,Sp%Var%n_times), 0, 0, &
                 Sp%Var%time(4,Sp%Var%n_times), &
                 base_day, base_sec, repeat_day, repeat_sec, lcal360)
  ! Find base time (days and seconds) for current model time
  CALL time2sec (model_time(1), model_time(2), model_time(3), &
                 model_time(4), model_time(5), model_time(6), &
                 0, 0, base_day, base_sec, lcal360)
  ! Find equivalent look-up time in repeat period
  j = -1
  outer: DO
    ! Subtract required number of repeat periods from model time
    CALL sec2time(j*repeat_day, j*repeat_sec, base_day, base_sec, &
      p_year, p_month, p_day, p_hour, p_minute, p_second, p_day_number, &
      lcal360)
    ! Search through repeated values for equivalent model time
    DO i=Sp%Var%n_times - Sp%Var%n_repeat_times + 1, Sp%Var%n_times
      IF ( ( Sp%Var%time(1, i) >  p_year )  .OR. &
           ( Sp%Var%time(1, i) == p_year  .AND.  &
             Sp%Var%time(2, i) >  p_month ) .OR. &
           ( Sp%Var%time(1, i) == p_year  .AND.  &
             Sp%Var%time(2, i) == p_month .AND.  &
             Sp%Var%time(3, i) >  p_day )   .OR. &
           ( Sp%Var%time(1, i) == p_year  .AND.  &
             Sp%Var%time(2, i) == p_month .AND.  &
             Sp%Var%time(3, i) == p_day   .AND.  &
             Sp%Var%time(4, i) >  p_hour*3600 +  &
                                  p_minute*60 +  &
                                  p_second ) ) EXIT outer
      i_time=i
    END DO
    j = j - 1
  END DO outer
END IF

WRITE(umMessage, '(a,i0.4,a,i2.2,a,i2.2,a,i0,a)') &
  'Setting solar spectral variation at look-up time: ', &
  Sp%Var%time(1, i_time),'-', Sp%Var%time(2, i_time),'-', &
  Sp%Var%time(3, i_time),' : ',Sp%Var%time(4, i_time),' sec'
CALL umPrint(umMessage, src='solvar')

! Set solar constant
solar_constant = Sp%Var%total_solar_flux(i_time)
WRITE(umMessage, '(a,1pe16.9)') 'Solar constant: ', solar_constant
CALL umPrint(umMessage, src='solvar')

! Set normalised flux in each spectral band
Sp%Solar%solar_flux_band = 0.0
DO sub_band=1, Sp%Var%n_sub_band
  band = Sp%Var%index_sub_band(1, sub_band)
  Sp%Solar%solar_flux_band(band) = Sp%Solar%solar_flux_band(band) &
    + Sp%Var%solar_flux_sub_band(sub_band, i_time)
END DO

! Set k-term weights (only given for special k-terms that
! represent the absorption in the sub-band) and Rayleigh coefficients
CALL umPrint('Band Fraction   : Rayleigh coeff  : Fraction per k-term', &
  src='solvar')
Sp%Rayleigh%rayleigh_coeff = 0.0
DO sub_band=1, Sp%Var%n_sub_band
  band = Sp%Var%index_sub_band(1, sub_band)
  term = Sp%Var%index_sub_band(2, sub_band)
  IF (term > 0) THEN
    Sp%Gas%w(term, band, 1) = Sp%Var%solar_flux_sub_band(sub_band, i_time) &
      / Sp%Solar%solar_flux_band(band)
    Sp%Rayleigh%rayleigh_coeff(band) = Sp%Rayleigh%rayleigh_coeff(band) &
      + Sp%Var%rayleigh_coeff(sub_band, i_time) * Sp%Gas%w(term, band, 1)
  ELSE
    Sp%Rayleigh%rayleigh_coeff(band)=Sp%Var%rayleigh_coeff(sub_band, i_time)
  END IF
  IF (sub_band == Sp%Var%n_sub_band) THEN
    next_band = band + 1
  ELSE
    next_band = Sp%Var%index_sub_band(1, sub_band+1)
  END IF
  IF (next_band > band) THEN
    IF (term > 0) THEN
      WRITE(umMessage, '(2(1pe16.9,a),61(1pe16.9))') &
        Sp%Solar%solar_flux_band(band), ' :', &
        Sp%Rayleigh%rayleigh_coeff(band), ' :', &
        Sp%Gas%w(1:term, band, 1)
    ELSE
      WRITE(umMessage, '(1pe16.9,a,1pe16.9)') &
        Sp%Solar%solar_flux_band(band), ' :', &
        Sp%Rayleigh%rayleigh_coeff(band)
    END IF
    CALL umPrint(umMessage, src='solvar')
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solvar

END MODULE solvar_mod
