! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme.
!
!  Contains the following subroutines:
!
!  lsp_subgrid (Subbgrid-scale set ups and checks)
!  lsp_qclear  (Calculates mean water vapour of clear-sky portion of gridbox)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

MODULE lsp_subgrid_mod

IMPLICIT NONE


! NOTE: This interface is also called from outside the lsp scheme by:
!       atmosphere/UKCA/ukca_main1-ukca_main1.F90
!       atmosphere/GLOMAP_CLIM/glomap_clim_fields_mod.F90
INTERFACE lsp_qclear
  MODULE PROCEDURE lsp_qclear_64b, lsp_qclear_32b
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_SUBGRID_MOD'

CONTAINS

SUBROUTINE lsp_subgrid(                                                       &
  points,                                                                     &
                                          ! Number of points
  q, qcf_cry, qcf_agg, qcftot, t,                                             &
                                          ! Water contents and temp
  qsl, qs,                                                                    &
                                          ! Saturated water contents
  q_ice, q_clear, q_ice_1, q_ice_2,                                           &
                                          ! Local vapour contents
  area_liq,area_mix,area_ice,area_clear,                                      &
                                              ! Cloud frac partitions
  area_ice_1, area_ice_2,                                                     &
                                          ! Subdivision of area_ice
  areamix_over_cfliq,                                                         &
                                          ! area_mix/cfliq
  rain_liq,rain_mix,rain_ice,rain_clear,                                      &
                                              ! Rain overlap partitions
  cftot, cfliq, cfice, cficei,                                                &
                                          ! Cloud fractions for
  frac_ice_above,                                                             &
                                          ! partition calculations
  cf, cff, rainfrac,                                                          &
                                          ! Cloud and rain fractions
                                          ! for updating
  lsrcp,                                                                      &
                                          ! Latent heat of sublim./cp
  rhcpt                                                                       &
                                          ! RH crit values
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod, ONLY: qcfmin, ice_width, zerodegc,                            &
                      zero, half, one, two

! Cloud modules- logicals and integers
USE cloud_inputs_mod,  ONLY: l_subgrid_qv

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

! Dr Hook Modules
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE

! Purpose:
!   Perform the subgrid-scale setting up calculations

! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!   Calculates the overlaps within each gridbox between  the cloud
!   fraction prognostics and rainfraction diagnostic.
!

! Description of Code:
!   Fortran95.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! The subgrid calculations are a necessary step to calculating the
! subsequent deposition and sublimation transfers and for setting up
! the partition information that is used by the subsequent transfers.

! Subroutine Arguments

INTEGER, INTENT(IN) ::                                                        &
  points            ! Number of points to calculate

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  qs(points),                                                                 &
                        ! Saturated mixing ratio wrt ice
  qsl(points),                                                                &
                        ! Saturated mixing ratio wrt liquid
  cfliq(points),                                                              &
                        ! Fraction of gridbox with liquid cloud
  rainfrac(points),                                                           &
                        ! Fraction of gridbox containing rain
  frac_ice_above(points),                                                     &
                             ! Ice cloud in level above this one
  lsrcp,                                                                      &
                        ! Latent heat of sublimation
                        ! / heat capacity of air / K
  rhcpt(points)     ! RH crit values

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  q(points),                                                                  &
                        ! Vapour content / kg kg-1
  qcf_cry(points),                                                            &
                        ! Ice crystal content / kg kg-1
  qcf_agg(points),                                                            &
                        ! Ice aggregate content / kg kg-1
  qcftot(points),                                                             &
                        ! Total ice content before advection / kg kg-1
  t(points),                                                                  &
                        ! Temperature / K
  cf(points),                                                                 &
                        ! Current cloud fraction
  cff(points)       ! Current ice cloud fraction

REAL (KIND=real_lsprec), INTENT(OUT) ::                                       &
  q_clear(points),                                                            &
                        ! Local vapour in clear-sky region / kg kg-1
  q_ice(points),                                                              &
                        ! Local vapour in ice-only region  / kg kg-1
  q_ice_1(points),                                                            &
                        ! Local vapour in ice-only regions that are:
  q_ice_2(points),                                                            &
                        !   1, depositing; and 2, subliming.
  cftot(points),                                                              &
                        ! Modified cloud fraction for partition calc.
  cfice(points),                                                              &
                        ! Modified ice cloud frac. for partition calc.
  cficei(points),                                                             &
                        ! 1/cfice
  area_liq(points),                                                           &
                        ! Frac of gridbox with liquid cloud but no ice
  area_mix(points),                                                           &
                        ! Frac of gridbox with liquid and ice cloud
  area_ice(points),                                                           &
                        ! Frac of gridbox with ice cloud but no liquid
  area_clear(points),                                                         &
                        ! Frac of gridbox with no cloud
  area_ice_1(points),                                                         &
                        ! Frac of gridbox where ice-only cloud is:
  area_ice_2(points),                                                         &
                        !  1, depositing; and 2, subliming.
  areamix_over_cfliq(points),                                                 &
                        ! area_mix/cfliq
  rain_liq(points),                                                           &
                        ! Frac of gbox with rain and liquid but no ice
  rain_mix(points),                                                           &
                        ! Frac of gbox with rain and liquid and ice
  rain_ice(points),                                                           &
                        ! Frac of gbox with rain and ice but no liquid
  rain_clear(points)! Frac of gbox with rain but no condensate

! Local Variables

INTEGER ::                                                                    &
  i                 ! Loop counter

REAL (KIND=real_lsprec) ::                                                    &
  temp7(points),                                                              &
                        ! Temporary in width of PDF calculation
  width(points)     ! Full width of vapour distribution in ice and
                        ! clear sky.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_SUBGRID'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1, points

      !-----------------------------------------------
      ! Check that ice cloud fraction is sensible.
      !-----------------------------------------------
      ! Difference between the way PC2 and non-PC2 code operates
      ! is kept here in order to be tracable across model versions.
      ! However, perhaps the code ought to be the same.
        ! 0.001 is to avoid divide by zero problems

  cfice(i)  = MAX( cff(i), 0.001_real_lsprec )
  cficei(i) = one/cfice(i)
  cftot(i)  = cf(i)
  cftot(i)  = MIN( MAX( cftot(i), cfice(i) ),( cfice(i) + cfliq(i)) )

    ! -----------------------------------------------
    ! Calculate overlaps of liquid, ice and rain fractions
    ! -----------------------------------------------
  area_liq(i) = MAX(cftot(i)-cfice(i),zero)
  area_mix(i) = MAX(cfice(i)+cfliq(i)-cftot(i),zero)

  IF (cfliq(i) /= zero) THEN
    areamix_over_cfliq(i) = area_mix(i)/cfliq(i)
  END IF

  IF (cfice(i) == cftot(i)) THEN
    area_mix(i)           = cfliq(i)
    areamix_over_cfliq(i) = one
  END IF

      ! Remove tiny mixed phase areas
  IF (area_mix(i) < 0.0001_real_lsprec) THEN
    area_mix(i)           = zero
    areamix_over_cfliq(i) = zero
  END IF

  area_ice(i) = MAX(cftot(i)-cfliq(i),zero)
  area_clear(i) = MAX(one-cftot(i),zero)
  rain_liq(i) = MAX(MIN(area_liq(i),rainfrac(i)),zero)
  rain_mix(i) = MAX(MIN(area_mix(i),rainfrac(i)-rain_liq(i)),zero)
  rain_ice(i) =                                                               &
    MAX(MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i)),zero)
  rain_clear(i) =                                                             &
    MAX(rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i),zero)

END DO

CALL lsp_qclear(                                                              &
!     Input fields
     q, qs, qsl, qcftot, cfliq, cftot, rhcpt,                                 &
!     Output field
     q_clear,                                                                 &
!     Array dimensions
     points)

DO i = 1, points
  IF (cfliq(i)  <   one) THEN

    IF (area_ice(i)  >   zero) THEN
      IF (l_subgrid_qv) THEN
        width(i) = two *(one-rhcpt(i))*qsl(i)                                 &
            *MAX((one-half*qcftot(i)/(ice_width * qsl(i))), 0.001_real_lsprec)

        ! The full width cannot be greater than 2q because otherwise
        ! part of the gridbox would have negative q. Also ensure that
        ! the full width is not zero (possible if rhcpt is 1).
        ! 0.001 is to avoid divide by zero problems
        width(i) = MIN(width(i), MAX(two*q(i),0.001_real_lsprec*qs(i)))

      ! If there is no partitioning of qv between clear-sky and
      ! ice-cloud, then use the same value of qv everywhere outside
      ! of the liquid-cloud
        q_ice(i) = (q(i)-cfliq(i)*qsl(i)-area_clear(i)*q_clear(i))            &
                   / area_ice(i)
      ELSE
        q_ice(i) = q_clear(i)
      END IF !l_subgrid_qv
    ELSE
      q_ice(i) = zero               ! q_ice is a dummy value here
    END IF  ! area_ice gt 0

  ELSE ! cf_liq lt 1

        ! -----------------------------------------------
        ! Specify dummy values for q_clear and q_ice
        ! -----------------------------------------------
    width(i)   = one
    q_clear(i) = zero
    q_ice(i)   = zero

  END IF ! cf_liq lt 1

      ! -------------------------------------------------
      ! Remove any small amount of ice to be tidy.
      ! -------------------------------------------------
      ! If QCF is less than QCFMIN and isn't growing by deposition
      ! (assumed to be given by RHCPT) then evaporate it.
  IF ((qcf_cry(i)+qcf_agg(i)) <  qcfmin) THEN
    IF (t(i) >  zerodegc .OR.                                                 &
       (q_ice(i)  <=  qs(i) .AND. area_mix(i)  <=  zero)                      &
       .OR. (qcf_cry(i)+qcf_agg(i)) <  zero) THEN
      q(i) = q(i) +qcf_cry(i)+qcf_agg(i)
      t(i) = t(i) - lsrcp * (qcf_cry(i)+qcf_agg(i))
      qcf_cry(i)=zero
      qcf_agg(i)=zero
    END IF ! T gt 0 etc.
  END IF ! qcf_cry+qcf_agg lt qcfmin

      ! -------------------------------------------------
      ! First estimate of partition sizes for ice sublimation
      ! and deposition and vapour contents within these partitions
      ! -------------------------------------------------
  IF (q_ice(i)  >   qs(i)) THEN
        ! First estimate is to use a deposition process
    area_ice_1(i) = area_ice(i)
    area_ice_2(i) = zero
    q_ice_1(i)    = q_ice(i)
    q_ice_2(i)    = qs(i)       ! Dummy value

  ELSE ! q_ice gt qs
        ! First estimate is to use a sublimation process
    area_ice_1(i) = zero
    area_ice_2(i) = area_ice(i)
    q_ice_1(i)    = qs(i)       ! Dummy value
    q_ice_2(i)    = q_ice(i)

  END IF ! q_ice gt qs


  ! -----------------------------------------------------
  ! If the water vapor is being partitioned between sublimating
  ! and deposition ice-regions, then calculate the detailed
  ! partition. Else there is no further
  ! subgrid partitioning of qv and the first-estimate(above) is used
  ! -----------------------------------------------------
  IF ( l_subgrid_qv ) THEN
      ! -------------------------------------------------
      ! Detailed estimate of partition sizes for ice sublimation
      ! and deposition and vapour contents within these partitions
      ! -------------------------------------------------
    IF (area_ice(i)  >   zero) THEN
      ! Temp7 is the estimate of the proportion of the gridbox
      ! which contains ice and has local q > than qs (wrt ice)
      temp7(i) = half*area_ice(i) + (q_ice(i)-qs(i)) / width(i)

      IF (temp7(i) >  zero .AND. temp7(i) <  area_ice(i)) THEN
          ! Calculate sizes of regions and q in each region
          ! These overwrite previous estimates
        area_ice_1(i) = temp7(i)
        area_ice_2(i) = area_ice(i) - area_ice_1(i)
        q_ice_1(i) = qs(i) + half * area_ice_1(i) * width(i)
        q_ice_2(i) = qs(i) - half * area_ice_2(i) * width(i)
      END IF ! temp7 gt 0 etc.

    END IF ! area_ice gt 0

  END IF   ! l_subgrid_qv=.true.


END DO  ! Points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_subgrid

!Note regarding variable precision:
!This routine is used in various places beyond the lsp scheme which may require
!either 64 or 32 bit calculation. Therefore, we can't use the real_lsprec 
!parameter as seen generally in the lsp_* modules. Instead, we create an 
!INTERFACE with both 32 and 64 bit versions available.

SUBROUTINE lsp_qclear_64b(                                                    &
!     Input fields
     q, qsmr, qsmr_wat, qcf, cloud_liq_frac, cloud_frac, rhcrit,              &
!     Output field
     q_clear,                                                                 &
!     Array dimensions
     npnts)

  ! Cloud modules
USE cloud_inputs_mod, ONLY: l_subgrid_qv
USE cloud_inputs_mod, ONLY: ice_width
USE um_types,         ONLY: real64

! Dr Hook Modules
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE
INTEGER, PARAMETER :: prec = real64
CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_QCLEAR_64B'
#include "lsp_subgrid_lsp_qclear.h"
RETURN
END SUBROUTINE lsp_qclear_64b

SUBROUTINE lsp_qclear_32b(                                                    &
!     Input fields
     q, qsmr, qsmr_wat, qcf, cloud_liq_frac, cloud_frac, rhcrit,              &
!     Output field
     q_clear,                                                                 &
!     Array dimensions
     npnts)

  ! Cloud modules
USE cloud_inputs_mod, ONLY: l_subgrid_qv
USE cloud_inputs_mod, ONLY: ice_width
USE um_types,         ONLY: real32

! Dr Hook Modules
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE
INTEGER, PARAMETER :: prec = real32
CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_QCLEAR_32B'
#include "lsp_subgrid_lsp_qclear.h"
RETURN
END SUBROUTINE lsp_qclear_32b

END MODULE lsp_subgrid_mod

