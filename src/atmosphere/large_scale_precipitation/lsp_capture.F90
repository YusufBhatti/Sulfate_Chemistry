! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Capture of raindrops by ice
! Subroutine Interface:
MODULE lsp_capture_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_CAPTURE_MOD'

CONTAINS

SUBROUTINE lsp_capture(                                                       &
  points, timestep,                                                           &
                                          ! Number of points and tstep
  qcf, qrain, qgraup, t,                                                      &
                                          ! Water contents and temp
  area_liq, area_mix, area_ice,                                               &
  cficei,                                                                     &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                               &
                                          ! Rain fraction information
  rain_ice, rain_clear,                                                       &
! cf, cff                                 ! Current cloud fractions for
!                                         ! updating
  rho, rhor, m0, tcg, tcgi,                                                   &
                                          ! Parametrization information
  corr, dhir,                                                                 &
  lfrcp , ice_type,                                                           &
                                          ! Microphysical information
  l_psd,                                                                      &
                                          ! Code options
  ptransfer,                                                                  &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                               &
                                          ! 1/(timestep*iterations)
!    &, cftransfer,cfftransfer            ! Cloud transfer diagnostics
  rftransfer,                                                                 &
                                          ! Rain fraction transfer
  l_use_agg_vt                                                                &
                                          ! Fallspeed branch choice
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod, ONLY: zerodegc, cx, constp, bi,                               &
                      zero, one, two

! Microphysics modules
USE mphys_constants_mod,  ONLY: ice_type_offset
USE mphys_inputs_mod,     ONLY: l_rainfall_as, l_sr2graup,                    &
                                l_mcr_qrain, l_mcr_qgraup,                    &
                                l_diff_icevt

! Dr Hook modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

! Large scale precipitation modules
USE lsp_moments_mod,      ONLY: lsp_moments

IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of capture of raindrops
!   by ice particles

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out a specified
!   distribution of raindrops.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

INTEGER, INTENT(IN) ::                                                        &
  points,                                                                     &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  timestep,                                                                   &
                        ! Timestep / s
  area_liq(points),                                                           &
                        ! Fraction of gridbox with liquid-only cloud
  area_mix(points),                                                           &
                        ! Fraction of gridbox with mixed phase cloud
  area_ice(points),                                                           &
                        ! Fraction of gridbox with ice-only cloud
  cficei(points),                                                             &
                        ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
    rho(points),                                                              &
                          ! Air density / kg m-3
    rhor(points),                                                             &
                          ! 1/Air density / m3 kg-1
    m0,                                                                       &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                              &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                             &
                          ! 1/tcg (no units)
    corr(points),                                                             &
                          ! Fall velocity correction factor (no units)
    dhir(points),                                                             &
                          ! Depth of layer / timestep  / m s-1
    lfrcp,                                                                    &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi          ! 1/(timestep*iterations)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  qcf(points),                                                                &
                        ! Ice water content    / kg kg-1
  qrain(points),                                                              &
                        ! Rain mixing ratio / kg kg-1
  qgraup(points),                                                             &
                        ! Graupel mixing ration / kg kg-1
  t(points),                                                                  &
                        ! Temperature / K
  rainfrac(points),                                                           &
                        ! Rain fraction
  rain_liq(points),                                                           &
                        ! Overlap fraction of rain and liquid
  rain_mix(points),                                                           &
                        ! Overlap frac of rain with mixed phase cloud
  rain_ice(points),                                                           &
                        ! Overlap fraction of rain and ice
  rain_clear(points),                                                         &
                        ! Overlap fraction of rain and clear-sky
  ptransfer(points),                                                          &
                        ! Mass rimed in this timestep / kg kg-1
  rftransfer(points)! Transfer of rain fraction this timestep
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfftransfer(points)! Cumulative ice cloud fraction increment

LOGICAL, INTENT(IN) ::                                                        &
  l_psd,                                                                      &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

! Local Variables

INTEGER ::                                                                    &
  i,                                                                          &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

REAL (KIND=real_lsprec) ::                                                    &
  dpr(points),                                                                &
                        ! Transfer of mixing ratio  / kg kg-1
  vr1(points),                                                                &
                        ! Average fall speed of raindrops  / m s-1
  vi1(points),                                                                &
                        ! Average fall speed of ice particles  / m s-1
  fv1(points),                                                                &
                        ! Average fall speed difference between
                        ! ice particles and raindrops  / m s-1
  lamr1(points),                                                              &
                        ! Reciprocal of slope parameter in raindrop
                        ! size distribution  / m
  lami1(points),                                                              &
                        ! Reciprocal of slope parameter in ice
                        ! particle size distribution  / m
  lam4,                                                                       &
                        ! 4.0 * lamr1
  m_bi_rat,                                                                   &
                        ! Ratio of moments m_bi_1 and m_bi
  lamfac1(points),                                                            &
                        ! Combination of lamr1 and lamr2
  rf_final(points),                                                           &
                        ! Rain fraction at end of the timestep
  m_0(points), m_1(points), m_2(points),                                      &
                        ! Moments of the ice particle size distributn.
  m_bi(points), m_bi_1(points),                                               &
                        ! bi and bi+1 moment of the generic ice psd.
  m_bi_di(points),                                                            &
                        ! bi+di moment of the generic ice size distn
  m_bic_dic(points)
                        ! bic+dic moment of the generic ice size distn

INTEGER, PARAMETER :: it_cry = 0 ! defines ice_type to be crystals

REAL (KIND=real_lsprec), PARAMETER    :: gr_thr = 1.0e-4_real_lsprec
                        ! threshold for converting rain water to
                        ! graupel

! Local compression variable
INTEGER ::                                                            &
  npts,                                                               &
                        ! Number of compressed points to compute
  c,                                                                  &
                        ! Compressed point index 
  ix(points)
                        ! Original index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_CAPTURE'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAY_FORTRAN) && (CRAY_FORTRAN <8004000)
!DIR$ NOVECTOR
#endif


cry_offset = ice_type * ice_type_offset

    !-----------------------------------------------
    ! Calculate moments of size distribution if appropriate
    !-----------------------------------------------
IF (l_psd) THEN
      ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
      ! ice particle size distribution
  CALL lsp_moments(points,rho,t,qcf,cficei,zero,m_0)
  CALL lsp_moments(points,rho,t,qcf,cficei,one,m_1)
  CALL lsp_moments(points,rho,t,qcf,cficei,two,m_2)
  CALL lsp_moments(points,rho,t,qcf,cficei,bi,m_bi)
  CALL lsp_moments(points,rho,t,qcf,cficei,cx(83),m_bi_1)
  IF (.NOT. l_diff_icevt) THEN
          ! Only one set of vt-D parameters
    CALL lsp_moments(points,rho,t,qcf,cficei,cx(82),m_bi_di)
  ELSE
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    CALL lsp_moments(points,rho,t,qcf,cficei,cx(182),m_bic_dic)
                        ! psd moment with crystal parameters
    CALL lsp_moments(points,rho,t,qcf,cficei,cx(82),m_bi_di)
                        ! psd moment with aggregate parameters
  END IF
END IF

! Identify the points that need to be calculated.
npts = 0
DO i = 1, points

  IF (qcf(i)  >   m0 .AND. qrain(i)  >   zero                                 &
      .AND. (rain_mix(i)+rain_ice(i))  >   zero                               &
      .AND. t(i)  <   zerodegc) THEN

    npts = npts + 1
    ix(npts) = i

  END IF ! qcf(i) >  m0 etc.

END DO

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)


      !-----------------------------------------------
      ! Calculate reciprocal of lambda for rain
      !-----------------------------------------------
  IF (l_mcr_qrain) THEN
        ! rain is a mixing ratio (kg kg-1)
    lamr1(i) = (rho(i) * constp(50) * qrain(i) / rainfrac(i))               &
               **(cx(52))

  ELSE
        ! rain is a flux (kg m-2 s-1)
    lamr1(i) = (qrain(i) * rho(i) * dhir(i)                                 &
               /( rainfrac(i)*constp(42)*corr(i)) )**(cx(42))

  END IF  ! l_mcr_qrain

      !-----------------------------------------------
      ! Calculate reciprocal of lambda for ice crystals
      !-----------------------------------------------
  IF (.NOT. l_psd) THEN
    lami1(i) = (rho(i) * qcf(i) * cficei(i)                                 &
           *constp(5+cry_offset)*tcgi(i))**(-cx(7+cry_offset))
  END IF


END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Calculate rain mass-weighted fallspeed
      !-----------------------------------------------
  !          If (l_mcr_qrain) then
  !            ! rain is a mixing ratio (kg kg-1)
  !            vr1(i) = (constp(41)/6.0) * corr(i) *                       &
  !     &             ( rho(i)*constp(50)*qrain(i)/rainfrac(i) )**cx(51)

  IF (l_rainfall_as) THEN ! Include Abel & Shipway terms.

    IF (l_mcr_qrain) THEN

      vr1(i) = ((lamr1(i)**cx(45))/(constp(53)*rainfrac(i))) *              &
       ( ( constp(54) / ((lamr1(i)+cx(56))**cx(59))) +                      &
         ( constp(55) / ((lamr1(i)+cx(57))**cx(60))))
    ELSE

      vr1(i) = ((lamr1(i)**cx(45))/constp(53))   *                          &
       ( ( constp(54) / ((lamr1(i)+cx(56))**cx(59))) +                      &
          ( constp(55) / ((lamr1(i)+cx(57))**cx(60))))

    END IF !l_mcr_qrain


  ELSE

        !Standard UM Rainfall parameters

    IF (l_mcr_qrain) THEN
          ! rain is a mixing ratio (kg kg-1)

      vr1(i) = (constp(41)/6.0_real_lsprec) * corr(i) *                     &
            ( rho(i)*constp(50)*qrain(i)/rainfrac(i) )**cx(51)

    ELSE
          ! rain is a flux (kg m-2 s-1)
      vr1(i) = corr(i) * constp(41)/6.0_real_lsprec *                       &
               (qrain(i) * rho(i) * dhir(i)                                 &
               /( rainfrac(i)*constp(42)*corr(i)) )**cx(41)

    END IF  ! l_mcr_qrain

  END IF !L_rainfall_as

      !-----------------------------------------------
      ! Calculate ice mass-weighted fallspeed
      !-----------------------------------------------
  IF (l_psd) THEN
        ! Use the generic PSD
    IF (.NOT. l_diff_icevt) THEN
       ! Only one set of vt-D parameters
       ! constp(82) = ci*ai
      vi1(i) = constp(82) * corr(i) * m_bi_di(i)                            &
                   / (rho(i) * qcf(i) * cficei(i))
    ELSE
      IF (l_use_agg_vt(i)) THEN
          ! Use aggregate parameters
        vi1(i) = constp(82) * corr(i) * m_bi_di(i)                          &
                   / (rho(i) * qcf(i) * cficei(i))
      ELSE
          ! Use crystal parameters
        vi1(i) = constp(182) * corr(i) * m_bic_dic(i)                       &
                   / (rho(i) * qcf(i) * cficei(i))
      END IF
    END IF

  ELSE
    vi1(i) = constp(4+cry_offset) * corr(i) *                               &
           (rho(i) * qcf(i) * cficei(i)                                     &
           * constp(5+cry_offset) * tcgi(i))**cx(3+cry_offset)
  END IF  ! l_psd

      !-----------------------------------------------
      ! Estimate the mean absolute differences in velocities
      !-----------------------------------------------
  fv1(i) = MAX(ABS(vr1(i)-vi1(i)),(vr1(i)+vi1(i))/8.0_real_lsprec)

END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

      !------------------------------------------------
      ! Calculate transfer
      !------------------------------------------------
  IF (l_psd) THEN

        ! Use the generic ice particle size distribution
        ! constp(86)=pi**2/24 x1r rho_water
        ! constp(87)=gamma(6+x4r)
        ! constp(88)=2 gamma(5+x4r)
        ! constp(89)=gamma(4+x4r)

    dpr(i) = constp(86) * fv1(i) *  timestep * rhor(i) *                    &
            (rain_mix(i)+rain_ice(i)) * lamr1(i)**(-cx(46)) *               &
            (lamr1(i)**cx(43) * constp(87) * m_0(i) +                       &
             lamr1(i)**cx(44) * constp(88) * m_1(i) +                       &
             lamr1(i)**cx(45) * constp(89) * m_2(i) )

  ELSE

    lamfac1(i) = constp(10+cry_offset) * constp(43) *                       &
             (lamr1(i)**cx(43) * lami1(i)**cx(8+cry_offset)) +              &
             constp(11+cry_offset) * constp(44) *                           &
             (lamr1(i)**cx(44) * lami1(i)**cx(9+cry_offset)) +              &
             constp(12+cry_offset) * constp(45) *                           &
             (lamr1(i)**cx(45) * lami1(i)**cx(10+cry_offset))

    dpr(i) = tcg(i) * constp(13+cry_offset) *                               &
           lami1(i)**(-cx(11+cry_offset)) *                                 &
           lamr1(i)**(-cx(46)) * fv1(i) * lamfac1(i) *                      &
           timestep * rhor(i) * (rain_mix(i)+rain_ice(i))

  END IF  ! l_psd

END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

        ! Limit transfer to the mass of rain that is available

  IF (rain_liq(i) == zero .AND. rain_clear(i) == zero) THEN

        ! All rain is through ice-only or mixed-phase cloud, so
        ! (rain_mix + rain_ice) / rainfrac must be 1.

    dpr(i) = MIN(dpr(i),qrain(i))

  ELSE

    dpr(i) = MIN(dpr(i),qrain(i) *                                          &
                      (rain_mix(i)+rain_ice(i))/rainfrac(i))

  END IF

  ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi


END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

      !------------------------------------------------
      ! Adjust ice and rain contents
      !------------------------------------------------

  IF ( l_mcr_qgraup .AND. l_sr2graup ) THEN

    ! Move contents to graupel rather than snow/ice

    IF (l_psd) THEN

      ! With the generic ice PSD active, we examine the
      ! mass-weighted mean diameters of the ice and rain.
      ! If the ice is bigger, we form snow. If the rain is
      ! bigger we form graupel.

      ! Calculate m_bi+1 / m_bi (mass-weighted diameter of ice)

      m_bi_rat = m_bi_1(i) / m_bi(i)

      ! Calculate 4 /lambda  (mass-weighted diameter of rain)
      ! As lamr1 = 1/lambda can just multiply this by 4:

      lam4 = 4.0_real_lsprec * lamr1(i)

      IF ( m_bi_rat >= lam4 ) THEN

        ! Snow m.m.d. bigger than rain
        ! Put transfer quantities as snow

        qcf(i) = qcf(i) + dpr(i)

      ELSE

        ! Rain m.m.d. bigger than snow
        ! Put transfer quantities as graupel

        qgraup(i) = qgraup(i) + dpr(i)

      END IF

    ELSE ! l_psd

      ! Without the generic ice PSD, we have two rain
      ! categories. In this case, we will follow a
      ! similar method to the Met Office LEM, transferring
      ! crystals colliding with small rain amounts to
      ! crystals and anything else to graupel.

      IF ( ice_type == it_cry .AND. qrain(i) < gr_thr ) THEN

        ! For crystals capturing low amounts of liquid
        ! water we shall move the result of the capture
        ! process to snow.

        ! If qr >= gr_thr (heavy rain) the quantities are
        ! assumed to freeze rapidly and produce graupel.

        qcf(i) = qcf(i) + dpr(i)

      ELSE

        ! The capture process is one of
        ! i  ) Crystals capturing large rain amounts or
        ! ii ) aggregates capturing any rain amount
        ! so we shall send the result of this capture process
        ! to graupel.

        qgraup(i) = qgraup(i) + dpr(i)

      END IF ! ice_type eq 0

    END IF ! l_psd

  ELSE    ! l_mcr_qgraup / l_sr2graup

    ! Graupel is not active, so move the amounts to
    ! the appropriate ice category.

    qcf(i) = qcf(i) + dpr(i)

  END IF  ! l_mcr_qgraup / l_sr2graup


END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

  ! Now remove the appropriate transfer amount
  ! from the rain category and adjust the
  ! temperature

  qrain(i) = qrain(i) - dpr(i)
  t(i)     = t(i)     + dpr(i) * lfrcp

  !------------------------------------------------
  ! Update cloud fractions
  !------------------------------------------------
  !         These are commented out since there is currently no
  !         cloud fraction update associated with the capture term.
  !          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
  !          cff_transfer_rate(i) = 0.0 / (timestep*iterations)

  !            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
  !            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)

  !            cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
  !            cff(i) = cff(i) +cff_transfer_rate(i)*timestep*iterations

             !-----------------------------------------------
             ! Update rain fractions
             !-----------------------------------------------
  !          These are commented out to ensure that rainfraction tends
  !          to a 0 or 1 behaviour as RHcrit tends to 1.
  !           rf_final(i) = rainfrac(i) * qrain(i) / (qrain(i)+dpr(i))
  !           rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))

  !             rainfrac(i)= rf_final(i)
  !             rain_liq(i)= min(area_liq(i),rainfrac(i))
  !             rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
  !             rain_ice(i)=
  !     &          min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
  !             rain_clear(i)=
  !     &            rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)


END DO  ! Points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_capture
END MODULE lsp_capture_mod
