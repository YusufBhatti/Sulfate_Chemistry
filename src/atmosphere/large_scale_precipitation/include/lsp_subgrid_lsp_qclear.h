! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Conversion between moments of PSD

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! ######################################################################
! Purpose:
! Returns the average relative humidity in the cloud-free (clear-sky)
! portion of gridboxes.
!
! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!
! Description of Code:
!   Fortran95
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! ######################################################################

! Subroutine Arguments
!
!-----------------------------------------------------------------------
! arguments with intent in. ie: input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN) :: npnts
! Points (dimensions) being processed

REAL (KIND=prec),   INTENT(IN)  :: q(npnts)
! water vapour mixing ratio (kg/kg)

REAL (KIND=prec),   INTENT(IN)  :: qsmr(npnts)
! Saturation mixing ratio, w.r.t. ice where T < 0C, (kg/kg)

REAL (KIND=prec),   INTENT(IN)  :: qsmr_wat(npnts)
! Saturation mixing ratio, w.r.t. water, (kg/kg)

REAL (KIND=prec),   INTENT(IN)  :: qcf(npnts)
! Frozen cloud condensate (kg/kg)

REAL (KIND=prec),   INTENT(IN)  :: cloud_liq_frac(npnts)
! Fraction of gridbox occupied by liquid or mixed phase cloud

REAL (KIND=prec),   INTENT(IN)  :: cloud_frac(npnts) !
! Fraction of gridbox occupied by cloud or any type

REAL (KIND=prec),   INTENT(IN)  :: rhcrit(npnts)
! Critical relative humidity for large-scale cloud formation (fraction)

!-----------------------------------------------------------------------
! arguments with INTENT out. ie: output variables.
!-----------------------------------------------------------------------

! relative humidity of the cloud-free portion of the gridboxes (fraction)
REAL (KIND=prec), INTENT(OUT) :: q_clear(npnts)

!-----------------------------------------------------------------------

!  Local variables

REAL (KIND=prec) :: area_ice        ! Ice-only cloud fraction
REAL (KIND=prec) :: width      ! Width of water vapour PDF
REAL (KIND=prec) :: qa         ! q in the portion of gridbox not containing 
                               ! liquid cloud
INTEGER :: k       ! loop integer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! loop over points
DO k = 1, npnts

  area_ice = MAX(cloud_frac(k) - cloud_liq_frac(k), 0.0_prec)

  IF (cloud_liq_frac(k) < 1.0_prec) THEN

    ! Derived from equation 34 in UMDP 26.
    qa = (q(k) - cloud_liq_frac(k)*qsmr_wat(k)) /                             &
         (1.0_prec - cloud_liq_frac(k))

    width = 2.0_prec *(1.0_prec-rhcrit(k))*qsmr_wat(k)                        &
            *MAX((1.0_prec-0.5_prec*qcf(k)/                                   &
                  (REAL(ice_width,prec) * qsmr_wat(k))),                      &
                  0.001_prec)

    ! The full width cannot be greater than 2q because otherwise
    ! part of the gridbox would have negative q. Also ensure that
    ! the full width is not zero (possible if rhcpt is 1).
    ! 0.001 is to avoid divide by zero problems
    width = MIN(width , MAX(2.0_prec*q(k),0.001_prec*qsmr(k)))

    ! Adjust for ice-only clouds, if subgrid partition
    ! of qv is selected. Else, use the same value of qv
    ! everywhere outside the liquid-cloud

    IF (area_ice  >  0.0_prec) THEN
      IF ( l_subgrid_qv ) THEN
        q_clear(k) = qa - 0.5_prec*width * area_ice
      ELSE
        q_clear(k) = qa
      END IF
    ELSE
      q_clear(k) = qa
    END IF ! area_ice gt 0

    ELSE ! cf_liq lt 1

      ! -----------------------------------------------
      ! If cloud-free or overcast then set q_clear = q
      ! -----------------------------------------------
      q_clear(k) = q(k)
  END IF ! cf_liq lt 1

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)