! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Ensure conservation of energy and optionally water

MODULE cor_engy_6a_mod

IMPLICIT NONE

! ------------------------------------------------------------------------------
!
! Description:
! Adjust the potential temperature increments and humidity increments
! to ensure the conservation of energy and water.
! By default the energy and water is corrected in a manner consistent
! with the global energy correction (see UMDP 85)
!
! It is possible to apply the correction in a manner consistent with
! the assumptions made in the convection scheme - this is useful for
! scheme development but is not recommended for general use.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COR_ENGY_6A_MOD'

CONTAINS

SUBROUTINE cor_engy_6a(npnts, nconv, nlev, index1,                      &
                       p_layer_boundaries, exner_layer_centres,         &
                       exner_rho, r_theta, r_rho, rho,                  &
                       r2rho_th, r2rho,                                 &
                       dr_across_th, dr_across_rh,                      &
                       dubydt, dvbydt, dqclbydt, dqcfbydt,              &
                       rain, snow, th, q, qcl, qcf, u, v,               &
                       !In/Out
                       dthbydt, dqbydt)

USE cv_derived_constants_mod, ONLY: ra2
USE planet_constants_mod, ONLY: g, cp
USE water_constants_mod, ONLY: lc, lf
USE cv_param_mod, ONLY: method_en_rho, method_en_mx_rho, method_en_qx_p
USE cv_dependent_switch_mod, ONLY: cor_method
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER,INTENT(IN) :: npnts         ! Full vector length
INTEGER,INTENT(IN) :: nconv         ! Number of convecting points
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations

INTEGER,INTENT(IN) :: index1(npnts) ! index of points with convection
REAL,INTENT(IN) :: p_layer_boundaries(npnts,0:nlev)   ! Pressure at
                                                      ! layer boundary (Pa)
REAL,INTENT(IN) :: exner_layer_centres(npnts,0:nlev)  ! Exner function
                                                      ! at layer centre
REAL,INTENT(IN) :: exner_rho(npnts,nlev)    ! Exner on rho levels
REAL,INTENT(IN) :: r_theta(npnts,0:nlev)    ! radius of theta levels (m)
REAL,INTENT(IN) :: r_rho(npnts,nlev)        ! radius of rho levels (m)
REAL,INTENT(IN) :: rho(npnts,nlev)          ! wet density for rho lev (kg/m3)
REAL,INTENT(IN) :: r2rho_th(npnts,nlev)     ! radius**2 wet density for
                                            ! theta lev (kg/m)
REAL,INTENT(IN) :: r2rho(npnts,nlev)        ! radius**2 wet density for
                                            ! rho lev (kg/m)
REAL,INTENT(IN) :: dr_across_th(npnts,nlev) ! thickness of theta levels (m)
REAL,INTENT(IN) :: dr_across_rh(npnts,nlev) ! thickness of rho levels (m)

REAL,INTENT(IN) :: dubydt(npnts,nlev)       ! Increment to u-wind on rho levels
                                            ! (m/s/s)
REAL,INTENT(IN) :: dvbydt(npnts,nlev)       ! Increment to u-wind on rho levels
                                            ! (m/s/s)
REAL,INTENT(IN) :: dqclbydt(npnts,nlev)     ! Increment to specific liquid water
                                            ! (kg/kg/s)
REAL,INTENT(IN) :: dqcfbydt(npnts,nlev)     ! Increment to specific frozen water
                                            ! (kg/kg/s)

REAL,INTENT(IN) :: rain(npnts)              ! rain at surface (kg/m**2/s)
REAL,INTENT(IN) :: snow(npnts)              ! Snow at surface (kg/m**2/s)

REAL,INTENT(IN) :: th(npnts,nlev)           ! potential temperature
                                            ! on theta levels (K)
! Specific quantities required to convert the specific increments into
! mixing ratio increments.
REAL,INTENT(IN) :: q(npnts,nlev)            ! Specific humidity on theta levels
                                            ! in kg/kg
REAL,INTENT(IN) :: qcl(npnts,nlev)          ! Specific liquid condensate on
                                            ! theta levels in kg/kg
REAL,INTENT(IN) :: qcf(npnts,nlev)          ! Specific frozen condensate on
                                            ! theta levels in kg/kg
REAL,INTENT(IN) :: u(npnts,nlev)            ! U wind (m/s)
REAL,INTENT(IN) :: v(npnts,nlev)            ! V wind (m/s)

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

REAL,INTENT(INOUT) :: dthbydt(npnts,nlev)   ! Increment to P. temperature
                                            ! due to convection (K/s)
                                            ! In:uncorrected
                                            ! Out:Energy Corrected
REAL,INTENT(INOUT) :: dqbydt(npnts,nlev)    ! Increment to specific
                                            ! water vapour (kg/kg/s)
                                            ! In:uncorrected
                                            ! Out:Mass corrected

!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER :: i,j,k, m               ! loop counters

REAL :: factor                    ! factor used to convert from specific
                                  ! to mixing ratios
REAL :: sumdqxbydt                ! Sum of increments to specific water
                                  ! quantities on theta levels (kg/kg/s)
REAL :: mv(nconv,nlev)            ! mixing ratio of water vapour
                                  ! on theta levels in kg/kg
REAL :: mcl(nconv,nlev)           ! mixing ratio of liquid condensate on
                                  ! theta levels in kg/kg
REAL :: mcf(nconv,nlev)           ! mixing ratio of frozen condensate on
                                  ! theta levels in kg/kg
REAL :: mv_rho(nconv,nlev)        ! mixing ratio of water vapour
                                  ! on theta levels in kg/kg
REAL :: mcl_rho(nconv,nlev)       ! mixing ratio of liquid condensate on
                                  ! rho levels in kg/kg
REAL :: mcf_rho(nconv,nlev)       ! mixing ratio of frozen condensate on
                                  ! rho levels in kg/kg
REAL :: dmvbydt(nconv,nlev)       ! Increment to mixing ratio of water vapour
                                  ! on theta levels (kg/kg/s)
REAL :: dmclbydt(nconv,nlev)      ! Increment to mixing ratio of liquid water
                                  ! on theta levels (kg/kg/s)
REAL :: dmvbydt_rho(nconv,nlev)   ! Increment to mixing ratio of water vapour
                                  ! on rho levels (kg/kg/s)
REAL :: dmclbydt_rho(nconv,nlev)  ! Increment to mixing ratio of liquid water
                                  ! on rho levels (kg/kg/s)
REAL :: dqbydt_rho(nconv,nlev)    ! Increment to specific water vapour
                                  ! on rho levels (kg/kg/s)
REAL :: dqclbydt_rho(nconv,nlev)  ! Increment to specific liquid water
                                  ! on rho levels (kg/kg/s)
REAL :: dqcfbydt_rho(nconv,nlev)  ! Increment to specific frozen water
                                  ! on rho levels (kg/kg/s)
REAL :: dthbydt_rho(nconv,nlev)   ! Increment to potential temperature
                                  ! on rho levels (kg/kg/s)
REAL :: dEKbydt(nconv)            ! Vertically integrated rate of change of
                                  ! kinetic energy (W/m2)
REAL :: dEIbydt(nconv)            ! Vertically integrated rate of change of
                                  ! internal energy (W/m2)
REAL :: dEMbydt(nconv)            ! Vertically integrated rate of change of
                                  ! moist energy (W/m2)
REAL :: flux_in(nconv)            ! Total energy flux into the column (W/m2)
REAL :: dEbydt(nconv)             ! Vertically integrated rate of change of
                                  ! total energy (W/m2)
REAL :: dTCWbydt(nconv)           ! Vertically integrated rate of change of
                                  ! total water (kg/m2/s)
REAL :: dTCVbydt(nconv)           ! Vertically integrated rate of change of
                                  ! water vapour (kg/m2/s)
REAL :: corr(nconv)               ! Energy correction required for the column
                                  ! (W/m2)
REAL :: TCWcorr(nconv)            ! Total water correction required for the
                                  ! column (kg/m2/s)
REAL :: r2_rhod_dr                ! r**2 rho_dry * dr on rho levels (kg)
REAL :: r2_rhow_dr                ! r**2 rho_wet * dr on rho levels (kg)
REAL :: weight1                   ! Weighting factor used in interpolation
                                  ! from theta to rho levels
REAL :: weight2                   ! Weighting factor used in interpolation
                                  ! from theta to rho levels
REAL :: weight3                   ! Weighting factor used in interpolation
                                  ! from theta to rho levels
REAL :: ww1                       ! Weighting factor used in interpolation
                                  ! from theta to rho levels
REAL :: ww2                       ! Weighting factor used in interpolation
                                  ! from theta to rho levels
REAL :: deltap                    ! Pressure difference between layer boundaries
                                  ! Only used for pressure based correction.

REAL :: dtscale(nconv) ! scaling factor used to correct the temperature
                       ! increments
REAL :: dmvscale(nconv)! scaling factor used to correct the water increments


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COR_ENGY_6A'

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)



!-------------------------------------------------------------------------
! WATER CORRECTION
!-------------------------------------------------------------------------

IF (cor_method == method_en_mx_rho) THEN
  !-------------------------------------------------------------------------
  ! Interpolate water increments to the rho grid for use in wate correction
  !-------------------------------------------------------------------------
  ! k=1
  DO j = 1, nconv
    ! arrays defined on all points, npnts, need to be indexed with i
    ! arrays defined on convecting points only, nconv, need to be indexed with j
    i                   = index1(j)
    ! assume bottom rho level value equal to bottom theta level value
    dqbydt_rho(j,1)     = dqbydt(i,1)
    dqclbydt_rho(j,1)   = dqclbydt(i,1)
    dqcfbydt_rho(j,1)   = dqcfbydt(i,1)
  END DO

  DO k=2,nlev
    DO j = 1, nconv
      i = index1(j)
      weight1             = r_theta(i,k) - r_rho(i,k)
      weight2             = r_rho(i,k)   - r_theta(i,k-1)
      weight3             = r_theta(i,k) - r_theta(i,k-1)
      ww1                 = weight1/weight3
      ww2                 = weight2/weight3

      dqbydt_rho(j,k)     = ww2 * dqbydt(i,k)     + ww1 * dqbydt(i,k-1)
      dqclbydt_rho(j,k)   = ww2 * dqclbydt(i,k)   + ww1 * dqclbydt(i,k-1)
      dqcfbydt_rho(j,k)   = ww2 * dqcfbydt(i,k)   + ww1 * dqcfbydt(i,k-1)
    END DO
  END DO
END IF !(cor_method == method_en_mx_rho)

IF (cor_method == method_en_mx_rho .OR. cor_method == method_en_qx_p) THEN
  !-------------------------------------------------------------------------
  ! Water correction: Initialise column integrals
  !-------------------------------------------------------------------------
  DO j=1,nconv
    dmvscale(j) = 0.0
    dTCVbydt(j) = 0.0
    dTCWbydt(j) = 0.0
  END DO
END IF !(cor_method == method_en_mx_rho .OR. cor_method == method_en_qx_p)


IF (cor_method == method_en_mx_rho) THEN
  !-------------------------------------------------------------------------
  ! Water correction on rho levels: Calculate column integrals
  !-------------------------------------------------------------------------
  DO k=1,nlev
    DO j=1,nconv
      i           = index1(j)
      r2_rhow_dr  = r2rho(i,k)*dr_across_rh(i,k)

      dTCVbydt(j) = dTCVbydt(j) + ra2 * r2_rhow_dr * dqbydt_rho(j,k)

      dTCWbydt(j) = dTCWbydt(j) + ra2 * r2_rhow_dr *                          &
                    (dqbydt_rho(j,k) + dqclbydt_rho(j,k) + dqcfbydt_rho(j,k))

    END DO
  END DO

ELSE IF (cor_method == method_en_qx_p) THEN
  !-------------------------------------------------------------------------
  ! Water correction on pressure levels: Calculate column integrals
  !-------------------------------------------------------------------------
  DO k=1,nlev
    DO j=1,nconv
      i           = index1(j)
      deltap      = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)

      dTCVbydt(j) = dTCVbydt(j) + deltap/g * dqbydt(i,k)

      dTCWbydt(j) = dTCWbydt(j) + deltap/g *                                  &
                    (dqbydt(i,k) + dqclbydt(i,k) + dqcfbydt(i,k))
    END DO
  END DO

END IF !cor_method

IF (cor_method == method_en_mx_rho .OR. cor_method == method_en_qx_p) THEN
  !-------------------------------------------------------------------------
  ! Water correction: calculate scaling
  !-------------------------------------------------------------------------
  DO j=1,nconv
    i = index1(j)

    TCWcorr(j)  = - rain(i) - snow(i) - dTCWbydt(j)

    IF (ABS(dTCVbydt(j)) > 1.0e-6*ABS(TCWcorr(j))) THEN
      dmvscale(j)  = 1.0 + TCWcorr(j)/dTCVbydt(j)
    ELSE
      dmvscale(j)  = 1.0
    END IF

    IF (dmvscale(j) > 1.25 .OR. dmvscale(j) < 0.75) THEN
      ! If the scaling is too large then do not bother.
      dmvscale(j) = 1.0
    END IF

  END DO

  !-------------------------------------------------------------------------
  ! Water correction: Apply correction
  !-------------------------------------------------------------------------
  DO k=1,nlev
    DO j=1,nconv
      i           = index1(j)
      dqbydt(i,k) = dqbydt(i,k)*dmvscale(j)
    END DO  ! nconv
  END DO ! nlev

END IF !cor_method

!-------------------------------------------------------------------------
! ENERGY CORRECTION
!-------------------------------------------------------------------------

IF (cor_method == method_en_mx_rho .OR. cor_method == method_en_rho) THEN
  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Convert the specific quantities and
  ! specific water increments to mass mixing ratio and mass mixing ratio
  ! increments for use in the energy correction
  !-------------------------------------------------------------------------
  DO k = 1, nlev
    DO j = 1, nconv
      i             = index1(j)
      factor        = 1.0/(1.0-q(i,k)-qcl(i,k)-qcf(i,k))
      sumdqxbydt    = dqbydt(i,k)+dqclbydt(i,k)+dqcfbydt(i,k)

      mv(j,k)       = factor*q(i,k)
      mcl(j,k)      = factor*qcl(i,k)
      mcf(j,k)      = factor*qcf(i,k)
      dmvbydt(j,k)  = factor*(dqbydt(i,k)   + q(i,k)  * factor * sumdqxbydt)
      dmclbydt(j,k) = factor*(dqclbydt(i,k) + qcl(i,k)* factor * sumdqxbydt)
    END DO
  END DO

  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Interpolate mixing ratios,
  ! mixing ratio increments and theta increment to the rho grid
  ! to be used in the energy correction
  !-------------------------------------------------------------------------

  ! k=1
  DO j = 1, nconv
    ! arrays defined on all points, npnts, need to be indexed with i
    ! arrays defined on convecting points only, nconv, need to be indexed with j
    i = index1(j)
    ! assume bottom rho level value equal to bottom theta level value
    mv_rho(j,1)         = mv(j,1)
    mcl_rho(j,1)        = mcl(j,1)
    mcf_rho(j,1)        = mcf(j,1)
    dthbydt_rho(j,1)    = dthbydt(i,1)
    dmvbydt_rho(j,1)    = dmvbydt(j,1)
    dmclbydt_rho(j,1)   = dmclbydt(j,1)
  END DO

  DO k=2,nlev
    DO j = 1, nconv
      i                   = index1(j)
      weight1             = r_theta(i,k) - r_rho(i,k)
      weight2             = r_rho(i,k)   - r_theta(i,k-1)
      weight3             = r_theta(i,k) - r_theta(i,k-1)
      ww1                 = weight1/weight3
      ww2                 = weight2/weight3

      mv_rho(j,k)         = ww2 * mv(j,k)         + ww1 * mv(j,k-1)
      mcl_rho(j,k)        = ww2 * mcl(j,k)        + ww1 * mcl(j,k-1)
      mcf_rho(j,k)        = ww2 * mcf(j,k)        + ww1 * mcf(j,k-1)
      dthbydt_rho(j,k)    = ww2 * dthbydt(i,k)    + ww1 * dthbydt(i,k-1)
      dmvbydt_rho(j,k)    = ww2 * dmvbydt(j,k)    + ww1 * dmvbydt(j,k-1)
      dmclbydt_rho(j,k)   = ww2 * dmclbydt(j,k)   + ww1 * dmclbydt(j,k-1)
    END DO
  END DO

END IF  !(cor_method == method_en_mx_rho .OR. cor_method == method_en_rho)

!-------------------------------------------------------------------------
! Energy correction: Initialise column integrals
!-------------------------------------------------------------------------
DO j=1,nconv
  dtscale(j)  = 0.0
  dEKbydt(j)  = 0.0
  dEIbydt(j)  = 0.0
  dEMbydt(j)  = 0.0
END DO

IF (cor_method == method_en_mx_rho .OR. cor_method == method_en_rho) THEN
  !-------------------------------------------------------------------------
  ! Energy correction on rho levels: Calculate column integrals
  !-------------------------------------------------------------------------
  DO k=1,nlev
    DO j=1,nconv
      i           = index1(j)
      !-------------------------------------------------------------------------
      ! This expression uses an approximation to dry rho and therefore
      ! is not quite as accurate as it ideally should be.
      ! It inconsistently uses start of timestep wet density combined with the
      ! latest mixing ratios. It also ignores prognostic rain. A more accurate
      ! expression should be considered in future.
      !-------------------------------------------------------------------------

      r2_rhod_dr  = r2rho(i,k)*dr_across_rh(i,k) /                            &
                    (1.0+mv_rho(j,k)+mcl_rho(j,k)+mcf_rho(j,k))

      dEKbydt(j)  = dEKbydt(j) + ra2 * r2_rhod_dr *                           &
                    ( u(i,k)*dubydt(i,k) + v(i,k)*dvbydt(i,k) )

      dEMbydt(j)  = dEMbydt(j) + ra2 * r2_rhod_dr *                           &
                    ( (lc+lf)*dmvbydt_rho(j,k) + lf*dmclbydt_rho(j,k) )

      dEIbydt(j)  = dEIbydt(j) + ra2 * r2_rhod_dr *                           &
                    cp*exner_rho(i,k)*dthbydt_rho(j,k)

    END DO
  END DO
ELSE IF (cor_method == method_en_qx_p) THEN
  !-------------------------------------------------------------------------
  !Energy correction on pressure levels: Calculate column integrals
  !-------------------------------------------------------------------------
  DO k=1,nlev
    DO j=1,nconv
      i           = index1(j)
      deltap      = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)

      dEKbydt(j)  = dEKbydt(j) + deltap/g*                                    &
                    ( u(i,k)*dubydt(i,k) + v(i,k)*dvbydt(i,k) )

      dEMbydt(j)  = dEMbydt(j) + deltap/g*                                    &
                    ( (lc+lf)*dqbydt(i,k) + lf*dqclbydt(i,k) )

      dEIbydt(j)  = dEIbydt(j) + deltap/g*                                    &
                    cp*exner_layer_centres(i,k)*dthbydt(i,k)

    END DO
  END DO
END IF !cor_method

!-------------------------------------------------------------------------
! Energy correction: calculate scaling
!-------------------------------------------------------------------------
DO j=1,nconv
  i           = index1(j)
  flux_in(j)  = -lf*rain(i)
  dEbydt(j)   = dEKbydt(j) + dEMbydt(j) + dEIbydt(j)

  corr(j)     = flux_in(j) - dEbydt(j)

  IF (ABS(dEIbydt(j)) > 1.0e-6*ABS(corr(j))) THEN
    dtscale(j)  = 1.0 + corr(j)/dEIbydt(j)
  ELSE
    dtscale(j)  = 1.0
  END IF

  IF (dtscale(j) > 1.25 .OR. dtscale(j) < 0.75) THEN
    !   If the scaling is too large then do not bother.
    dtscale(j) = 1.0
  END IF

END DO

!-------------------------------------------------------------------------
! Energy correction: Apply correction
!-------------------------------------------------------------------------
DO k=1,nlev
  DO j=1,nconv
    i             = index1(j)
    dthbydt(i,k)  = dthbydt(i,k)*dtscale(j)
  END DO  ! nconv
END DO ! nlev


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE cor_engy_6a

END MODULE cor_engy_6a_mod
