! Subroutine Arguments:
INTEGER, INTENT(IN)              :: npnts
REAL (KIND=prec),    INTENT(IN)  :: t(npnts), p(npnts)
REAL (KIND=prec),    INTENT(OUT) :: qs(npnts)

! Local scalars
INTEGER                          :: itable, i
REAL (KIND=prec)                 :: atable, fsubw, tt

! Local parameters allows simple templating of KINDs
REAL (KIND=prec), PARAMETER      :: one   = 1.0,                              &
                                    pconv = 1.0e-8,                           &
                                    term1 = 4.5,                              &
                                    term2 = 6.0e-4,                           &
                                    term3 = 1.1

DO i=1, npnts
  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, fsubw.
  ! This formula is taken from equation A4.7 of Adrian Gill's book:
  ! atmosphere-ocean dynamics. Note that his formula works in terms
  ! of pressure in mb and temperature in celsius, so conversion of
  ! units leads to the slightly different equation used here.
  fsubw = one + pconv * p(i) *                                                &
  (term1 + term2 * (t(i) - zerodegc) * (t(i) - zerodegc))

  ! Use the lookup table to find saturated vapour pressure. Store it in qs.
  tt     = MAX(t_low,t(i))
  tt     = MIN(t_high,tt)
  atable = (tt - t_low + delta_t) / delta_t
  itable = atable
  atable = atable - itable
  qs(i)  = (one - atable) * es(itable) + atable * es(itable+1)

  ! Multiply by fsubw to convert to saturated vapour pressure in air
  ! (equation A4.6 OF Adrian Gill's book).
  qs(i)  = qs(i) * fsubw

  ! Now form the accurate expression for qs, which is a rearranged
  ! version of equation A4.3 of Gill's book.
  ! Note that at very low pressures we apply a fix, to prevent a
  ! singularity (qsat tends to 1. kg/kg).
  qs(i)  = (repsilon * qs(i)) / (MAX(p(i), term3 * qs(i)) - qs(i))
END DO
RETURN