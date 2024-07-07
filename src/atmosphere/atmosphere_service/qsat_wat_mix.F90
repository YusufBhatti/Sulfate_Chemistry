! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Saturation Specific Humidity Scheme (Qsat_Wat): Vapour to Liquid.
! Subroutine Interface:
SUBROUTINE qsat_wat_mix (                                         &
!      Output field
        QmixS                                                           &
!      Input fields
      , t, p                                                            &
!      Array dimensions
      , npnts                                                           &
!      logical control
      , lq_mix                                                          &
        )

USE qsat_wat_data, ONLY:                                          &
   T_low,                                                         &
   T_high,                                                        &
   delta_T,                                                       &
   es

USE conversions_mod, ONLY: zerodegc

USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon

USE vectlib_mod, ONLY: oneover_v
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Purpose:
!   Returns a saturation specific humidity or mixing ratio given a
!   temperature and pressure using the saturation vapour pressure
!   calculated using the Goff-Gratch formulae, adopted by the WMO as
!   taken from Landolt-Bornstein, 1987 Numerical Data and Functional
!   Relationships in Science and Technolgy. Group V/vol 4B meteorology.
!   Phyiscal and Chemical properties or air, P35.
!
!   Values in the lookup table are over water above and below 0 deg C.
!
!   Note : For vapour pressure over water this formula is valid for
!   temperatures between 373K and 223K. The values for saturated vapour
!   over water in the lookup table below are out of the lower end of
!   this range. However it is standard WMO practice to use the formula
!   below its accepted range for use with the calculation of dew points
!   in the upper atmosphere
!
! Method:
!   Uses lookup tables to find eSAT, calculates qSAT directly from that.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service
!
! Code Description:
!   Language: FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP No.29
!
! Declarations:
!
!  Global Variables:----------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
!
! arguments with intent in. ie: input variables.

INTEGER, INTENT(IN) :: npnts
! Points (=horizontal dimensions) being processed by qSAT scheme.

REAL, INTENT(IN)  :: t(npnts) !  Temperature (K).
REAL, INTENT(IN)  :: p(npnts) !  Pressure (Pa).

LOGICAL, INTENT(IN)  :: lq_mix
              !  .true. return qsat as a mixing ratio
              !  .false. return qsat as a specific humidity

! arguments with intent out

REAL, INTENT(OUT)   ::  QmixS(npnts)
       ! Output Saturation mixing ratio or saturation specific
       ! humidity at temperature T and pressure P (kg/kg).

!-----------------------------------------------------------------------
!  Local scalars
!-----------------------------------------------------------------------

INTEGER :: itable    ! Work variable

REAL :: atable       ! Work variable

REAL :: fsubw
      ! FACTOR THAT CONVERTS FROM SAT VAPOUR PRESSURE IN A PURE
      ! WATER SYSTEM TO SAT VAPOUR PRESSURE IN AIR.

REAL :: tt

REAL :: vect_in(npnts), vect_out(npnts)
      ! temp array for calculating reciprocals

REAL, PARAMETER :: R_delta_T = 1.0/delta_T

INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='QSAT_WAT_MIX'

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


! loop over points
DO i=1,npnts

  tt = MAX(T_low,t(i))
  tt = MIN(T_high,tt)
  atable = (tt - T_low + delta_T) * R_delta_T
  itable = atable
  atable = atable - itable

    !      Use the lookup table to find saturated vapour pressure, and store
    !      it in qmixs.

  QmixS(i)   = (1.0 - atable)    * es(itable)    +              &
       atable*es(itable+1)

  IF (lq_mix) THEN

      !      Compute the factor that converts from sat vapour pressure in a
      !      pure water system to sat vapour pressure in air, fsubw.
      !      This formula is taken from equation A4.7 of Adrian Gill's book:
      !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
      !      of pressure in mb and temperature in Celsius, so conversion of
      !      units leads to the slightly different equation used here.

    fsubw    = 1.0 + 1.0e-8*p(i)   * ( 4.5 + 6.0e-4*    &
               ( t(i) - zerodegc ) * ( t(i) - zerodegc) )

      !      Multiply by fsubw to convert to saturated vapour pressure in air
      !      (equation A4.6 of Adrian Gill's book).

    QmixS(i)   = QmixS(i)   * fsubw

      !      Now form the accurate expression for qmixs, which is a rearranged
      !      version of equation A4.3 of Gill's book.

      !-----------------------------------------------------------------------
      ! For mixing ratio,  rsat = epsilon *e/(p-e)
      ! e - saturation vapour pressure
      ! Note applying the fix to qsat for specific humidity at low pressures
      ! is not possible, this implies mixing ratio qsat tends to infinity.
      ! If the pressure is very low then the mixing ratio value may become
      ! very large.
      !-----------------------------------------------------------------------

    vect_in(i) = ( MAX(p(i),  1.1*QmixS(i))   - QmixS(i) )

  ELSE

      !      Compute the factor that converts from sat vapour pressure in a
      !      pure water system to sat vapour pressure in air, fsubw.
      !      This formula is taken from equation A4.7 of Adrian Gill's book:
      !      Atmosphere-Ocean Dynamics. Note that his formula works in terms
      !      of pressure in mb and temperature in Celsius, so conversion of
      !      units leads to the slightly different equation used here.

    fsubw    = 1.0 + 1.0e-8*p(i) * ( 4.5 + 6.0e-4*                            &
               ( t(i) - zerodegc ) * ( t(i) - zerodegc ))

      !      Multiply by fsubw to convert to saturated vapour pressure in air
      !      (equation A4.6 of Adrian Gill's book).

    QmixS(i)   = QmixS(i)   * fsubw

      !      Now form the accurate expression for qmixs, which is a rearranged
      !      version of equation A4.3 of Gill's book.

      !-----------------------------------------------------------------------
      ! For specific humidity,   qsat = epsilon*e/(p-(1-epsilon)e)
      !
      ! Note that at very low pressure we apply a fix, to prevent a
      ! singularity (qsat tends to 1. kg/kg).
      !-----------------------------------------------------------------------
    vect_in(i) = (MAX(p(i), QmixS(i))-one_minus_epsilon*QmixS(i))

  END IF
END DO

CALL oneover_v(npnts, vect_in, vect_out)

DO i=1, npnts
  QMixS(i) = repsilon*QmixS(i) * vect_out(i)
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE qsat_wat_mix
