! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Deposition of ice particles
! Subroutine Interface:
MODULE lsp_deposition_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_DEPOSITION_MOD'

CONTAINS

SUBROUTINE lsp_deposition(                                                    &
  points, timestep,                                                           &
                                          ! Number of points and tstep
  q, qcl, qcf, qcft, t, p, frac_in_category,                                  &
                                          ! Water contents, temp, pres
  q_ice_1, q_ice_2,                                                           &
                                          ! Subgrid-scale water c'tents
  area_ice_1, area_ice_2,                                                     &
                                          ! Subgrid-scale areas
  esi, qs, qsl,                                                               &
                                          ! Saturated quantities
  area_mix, cfliq, cfice, cficei, areamix_over_cfliq,                         &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  cf, cfl, cff,                                                               &
                                          ! Current cloud fractions for
                                          ! updating
  rho, tcg, tcgi,                                                             &
                                          ! Parametrization information
  corr2, rocor, lheat_correc_ice,                                             &
  lfrcp, lsrcp, ice_type,                                                     &
                                          ! Microphysical information
  l_psd,                                                                      &
                                          ! Code options
  ptransfer,                                                                  &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                               &
                                          ! 1/(timestep*iterations)
  cftransfer,cfltransfer,cfftransfer,                                         &
                                          ! Cloud transfer diagnostics
  l_use_agg_vt                                                                &
                                          ! Vt-D branch to use
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod, ONLY: apb1, apb2, apb3, hm_t_min, hm_t_max, hm_decay,         &
                      hm_rqcl, m0, cx, constp, zerodegc
USE lsprec_mod, ONLY: zero, half, one

  ! Microphysics modules- logicals and integers
USE mphys_constants_mod, ONLY: l_hallett_mossop, ice_type_offset
USE mphys_inputs_mod,    ONLY: l_diff_icevt

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

! Dr Hook Modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

! Large scale precipitation modules
USE lsp_moments_mod,     ONLY: lsp_moments

IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of ice deposition and
!   sublimation

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles in a distribution of vapour
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Source for ice. Sink for liquid water and vapour.
! Hallett Mossop process enhances growth when turned on -
! Check value of HM_T_MAX

! Subroutine Arguments


INTEGER, INTENT(IN) ::                                                        &
  points,                                                                     &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates)

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  timestep,                                                                   &
                        ! Timestep / s
  p(points),                                                                  &
                        ! Air pressure / N m-2
  esi(points),                                                                &
                        ! Saturated vapour pressure over ice / N m-2
  qs(points),                                                                 &
                        ! Saturated humidity wrt ice / kg kg-1
  qsl(points),                                                                &
                        ! Saturated humidity wrt liquid / kg kg-1
  q_ice_1(points),                                                            &
                        ! Mean vapour in ice only region / kg kg-1
  q_ice_2(points),                                                            &
                        ! Mean vapour in clear region / kg kg-1
  area_ice_1(points),                                                         &
                        ! Ice only area that is growing by deposition
  area_ice_2(points),                                                         &
                        ! Ice only area that is subliming
  area_mix(points),                                                           &
                        ! Fraction of gridbox with mixed phase cloud
  cfliq(points),                                                              &
                        ! Fraction of gridbox with liquid cloud
  cfice(points),                                                              &
                        ! Fraction of gridbox with ice cloud
  cficei(points),                                                             &
                        ! 1/Fraction of gridbox with ice cloud
  areamix_over_cfliq(points),                                                 &
                        ! area_mix(points)/cfliq(points)
  rho(points),                                                                &
                        ! Air density / kg m-3
  tcg(points),                                                                &
                        ! T dependent function in ice size dist'n
  tcgi(points),                                                               &
                        ! 1/tcg (no units)
  corr2(points),                                                              &
                        ! Temperature correction factor (no units)
  rocor(points),                                                              &
                        ! Combined fall and corr2 correction factor
  lheat_correc_ice(points),                                                   &
                        ! Ice latent heat correction factor
  frac_in_category(points),                                                   &
                        ! Fraction of the ice that is in this category
  lfrcp,                                                                      &
                        ! Latent heat of fusion
                        ! /heat capacity of air (cP) / K
  lsrcp,                                                                      &
                        ! Latent heat of sublimation/cP / K
  one_over_tsi          ! 1/(timestep*iterations)
LOGICAL, INTENT(IN) ::                                                        &
   l_use_agg_vt(points)
                         ! Determines which vt-D parameters to use

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  q(points),                                                                  &
                        ! Vapour content / kg kg-1
  qcl(points),                                                                &
                        ! Liquid water content / kg kg-1
  qcf(points),                                                                &
                        ! Ice water content in ice category to be
!                           updated    / kg kg-1
    qcft(points),                                                             &
                          ! Ice water in all ice categories
!                           (for cloud fraction calculations)
    t(points),                                                                &
                          ! Temperature / K
    cf(points),                                                               &
                          ! Current cloud fraction
    cfl(points),                                                              &
                          ! Current liquid cloud fraction
    cff(points),                                                              &
                          ! Current ice cloud fraction
    ptransfer(points)  ! Mass deposited in this timestep / kg kg-1

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  cftransfer(points),                                                         &
                         ! Cloud fraction increment this tstep
  cfltransfer(points),                                                        &
                         ! Liquid cloud fraction inc this tstep
  cfftransfer(points)! Ice cloud fraction inc this tstep

LOGICAL, INTENT(IN) ::                                                        &
  l_psd
                        ! Use generic ice particle size distribution

! Local Variables

INTEGER ::                                                                    &
  i,                                                                          &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

REAL (KIND=real_lsprec) ::                                                    &
  aplusb(points),                                                             &
                        ! A+B terms in diffusional growth equation
  tempw_dep(points),                                                          &
                        ! Temporary for available moisture
  tempw_sub(points),                                                          &
                        ! Temporary for saturation deficit
  pr02(points),                                                               &
                        ! Temporary in calculation of PSD slopes
  lamr1(points),                                                              &
                        ! Power of PSD slope
  lamr2(points),                                                              &
                        ! Power of PSD slope
  dqi(points),                                                                &
                        ! Temporary in calculating ice transfers
  dqil(points),                                                               &
                        ! Temporary in calculating liquid transfers
  dqi_dep(points),                                                            &
                        ! Deposition amount / kg kg-1
  dqi_sub(points),                                                            &
                        ! Sublimation amount / kg kg-1
  cfltemp(points),                                                            &
                        ! Temporary in calc of cfl and qcl change
  deltacf(points),                                                            &
                        ! Change in cf across timestep
  deltacfl(points)  ! Change in cfl across timestep

REAL (KIND=real_lsprec) ::                                                    &
  hm_rate,                                                                    &
                        ! Hallett-Mossop rate multiplier
  hm_normalize,                                                               &
                        ! Normalization for T func. in HM process
  m_1(points),                                                                &
                        ! 1st moment of generic particle size dist.
  m_0p5_dicp3(points),                                                        &
                        ! 1+(dic+1)/2 moment of generic PSD
  m_0p5_dip3(points)
                        ! 1+(di+1)/2 moment of generic PSD

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_DEPOSITION'


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


cry_offset=ice_type*ice_type_offset

    ! Use the generic ice particle size distribution
IF (l_psd) THEN

   ! Calculate the 1st (cx(84)) moment of the
   ! ice particle size distribution

  CALL lsp_moments(points,rho,t,qcf,cficei,cx(84),m_1)

  IF (.NOT. l_diff_icevt ) THEN
    ! Use only one vt-D relation
    ! Calculate the 1+0.5(di+1) (cx(85)) moment of the
    ! ice particle size distribution

    CALL lsp_moments(points,rho,t,qcf,cficei,cx(85),m_0p5_dip3)
  ELSE
     ! The vt-D relation which gives the
     ! least mass weighted mean fallspeed will be used so
     ! calculate both required moments

    CALL lsp_moments(points,rho,t,qcf,cficei,cx(185),m_0p5_dicp3)
             ! ice ventilation moment with crystal parameters

    CALL lsp_moments(points,rho,t,qcf,cficei,cx(85),m_0p5_dip3)
             ! ice ventilation moment with aggregate parameters
  END IF ! l_diff_icevt


END IF

    !-----------------------------------------------
    ! Hallett Mossop process normalisation
    !-----------------------------------------------
hm_normalize=one/(one-EXP((hm_t_min-hm_t_max)*hm_decay))


! Identify the points that need to be calculated.
npts = 0
DO i = 1, points

  IF (qcf(i) >  m0 .AND. t(i) <  zerodegc) THEN

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
      ! Diffusional growth parameters
      !-----------------------------------------------
  aplusb(i) = (apb1-apb2*t(i)) * esi(i)
  aplusb(i) = aplusb(i) + (t(i)**3) * p(i) * apb3

      !-----------------------------------------------
      ! Moisture available from subgrid scale calculation
      !-----------------------------------------------
  tempw_dep(i) = qsl(i) * area_mix(i)                                       &
               + MIN(q_ice_1(i),qsl(i)) * area_ice_1(i)                     &
               - qs(i) * (area_mix(i) + area_ice_1(i))
  tempw_sub(i) = (q_ice_2(i) - qs(i)) * area_ice_2(i)

      ! Only allow the sub saturation to be reduced by the
      ! fraction of ice in this category
  tempw_sub(i) = tempw_sub(i)*frac_in_category(i)

      !-----------------------------------------------
      ! Calculate transfer rates
      !-----------------------------------------------
  IF (l_psd) THEN
        ! Use generic particle size distribution
    IF (.NOT. l_diff_icevt) THEN
       ! Only one choice of vt-D parameters for ventilation
       ! constp(83) = 2 pi axial_ratio_correction
       ! constp(84) = ventilation coefficient 1
       ! constp(85) = ventilation coefficient 2
       !        * Sc^(1/3)*ci^0.5/viscosity0^0.5
      dqi(i) = constp(83) * t(i)**2 * esi(i)                                &
                          * (constp(84)*m_1(i)*corr2(i)                     &
                          +  constp(85)*rocor(i)*m_0p5_dip3(i))             &
                          / (qs(i) * aplusb(i) * rho(i))
    ELSE
      IF (l_use_agg_vt(i)) THEN
         ! Use aggregate parameters for ventilation
        dqi(i) = constp(83) * t(i)**2 * esi(i)                              &
                         * (constp(84)*m_1(i)*corr2(i)                      &
                         +  constp(85)*rocor(i)*m_0p5_dip3(i))              &
                         / (qs(i) * aplusb(i) * rho(i))
      ELSE
         ! Use crystal parameters for ventilation
        dqi(i) = constp(83) * t(i)**2 * esi(i)                              &
                         * (constp(84)*m_1(i)*corr2(i)                      &
                         +  constp(185)*rocor(i)*m_0p5_dicp3(i))            &
                         / (qs(i) * aplusb(i) * rho(i))
      END IF

    END IF ! l_diff_icevt

  ELSE
        ! Use particle size distribution based on intercepts
    pr02(i)=rho(i)*qcf(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
    lamr1(i) = pr02(i)**cx(4+cry_offset)
    lamr2(i) = pr02(i)**cx(5+cry_offset)
    dqi(i) = tcg(i) * constp(6+cry_offset) * t(i)**2 *                      &
             esi(i) * (constp(7+cry_offset) * corr2(i) *                    &
             lamr1(i) + constp(8+cry_offset) * rocor(i) *                   &
             lamr2(i)) / (qs(i) * aplusb(i) * rho(i))

  END IF  ! l_psd

  dqi_dep(i) = dqi(i) * tempw_dep(i)
  dqi_sub(i) = dqi(i) * tempw_sub(i)

END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)


  IF (dqi_dep(i) >  zero) THEN  ! Limits depend on whether
                                   ! deposition or sublimation

        !-----------------------------------------------
        ! Deposition is occuring
        !-----------------------------------------------

    IF (l_hallett_mossop .AND. ice_type  ==  0) THEN

          ! Only for crystals when Hallett-Mossop process
          ! is active.

          !-----------------------------------------------
          ! Hallett Mossop Enhancement
          !-----------------------------------------------

      IF ( (t(i)-zerodegc)  >=  hm_t_max) THEN  ! no enhancement
        hm_rate=zero
      ELSE IF ((t(i)-zerodegc)  <   hm_t_max                                &
         .AND. (t(i)-zerodegc)  >   hm_t_min) THEN
            ! Some enhancement occurs between temperature limits
        hm_rate = (one-EXP( (t(i)-zerodegc-hm_t_max)*hm_decay) )            &
                  * hm_normalize
      ELSE  ! Some enhancement at lower temperatures
        hm_rate = EXP( (t(i)-zerodegc - hm_t_min) * hm_decay )
      END IF

          ! Calculate enhancement factor for HM process
      hm_rate = one + hm_rate * qcl(i) * hm_rqcl
      dqi_dep(i) = dqi_dep(i) * hm_rate


    END IF  ! l_hallett_mossop .and. ice_type  ==  0

        !-----------------------------------------------
        ! Molecular diffusion to/from a surface is more efficient
        ! when a particle is at a molecular step. This is more
        ! likely for sublimation. For growth, reduce rate by 10%
        ! ----------------------------------------------
    dqi_dep(i)=0.9_real_lsprec*dqi_dep(i)

        !-----------------------------------------------
        ! Latent heat correction (equivalent to aL in LS cloud)
        !-----------------------------------------------
    tempw_dep(i)=tempw_dep(i)*lheat_correc_ice(i)

        ! Only allow the supersaturation to be reduced by
        ! the fraction of ice in this category
    tempw_dep(i) = tempw_dep(i)*frac_in_category(i)

        !-----------------------------------------------
        ! Calculate available moisture and transfer amount.
        !-----------------------------------------------
    IF (cfliq(i) >  zero) THEN
        ! Include liquid water contribution. Ignore latent
        ! heat correction from freezing liquid.
      tempw_dep(i)=tempw_dep(i)+qcl(i)*areamix_over_cfliq(i)                &
                  +MAX((q_ice_1(i)-qsl(i)),zero)*area_ice_1(i)
    END IF
    dqi_dep(i)=MIN(dqi_dep(i)*timestep,tempw_dep(i))

  END IF  ! dqi_dep gt 0.0

END DO


!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

  IF (dqi_sub(i)  <   zero) THEN
        !-----------------------------------------------
        ! Sublimation is occuring
        !-----------------------------------------------
        ! Limits are spare moisture capacity and QCF
        ! outside liquid cloud
    dqi_sub(i) = MAX( MAX( dqi_sub(i) * timestep,                           &
                           tempw_sub(i) * lheat_correc_ice(i) ),            &
                 -( qcf(i) * area_ice_2(i) * cficei(i) ))

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------
    deltacf(i) = area_ice_2(i) * ( SQRT( MAX(                               &
               one+dqi_sub(i)*cfice(i)/(qcft(i)*area_ice_2(i)),             &
                 zero) ) - one )

    ! The above calculation is ill-conditioned when the square
    ! root argument is slightly non-zero, which often occurs
    ! spuriously due to rounding error. The following test uses
    ! tolerances to detect when this situation has occured, and
    ! revises the calculation of deltacf accordingly:
    IF (ABS( dqi_sub(i)+qcft(i) )     < 1.0e-16_real_lsprec .AND.           &
        ABS( cfice(i)-area_ice_2(i) ) < 1.0e-12_real_lsprec) THEN
      deltacf(i) = -area_ice_2(i)
    END IF

    cff(i) = cff(i) + deltacf(i)
    cf (i) = cf (i) + deltacf(i)

    cftransfer(i) = cftransfer(i)  + deltacf(i) * one_over_tsi
    cfftransfer(i)= cfftransfer(i) + deltacf(i) * one_over_tsi

  END IF  ! dqi_sub lt 0.0

END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

    !-----------------------------------------------
    ! Calculate liquid water change
    !-----------------------------------------------
  IF (cfliq(i) >  zero .AND. area_mix(i) >  zero                            &
          .AND. qcl(i) >  zero) THEN

        ! Deposition removes some liquid water content
        ! First estimate of the liquid water removed is explicit

    dqil(i) = MAX (MIN ( dqi_dep(i)*area_mix(i)                             &
              /(area_mix(i)+area_ice_1(i)),                                 &
              qcl(i)*areamix_over_cfliq(i)) ,zero)

      !-----------------------------------------------
      ! Update liquid cloud fraction (and new liquid water est.)
      !-----------------------------------------------
        ! First estimate of the liquid cloud fraction is based
        ! on this explicit estimate of liquid lost by deposition

    cfltemp(i)=cfl(i)*SQRT(MAX(one-dqil(i)/qcl(i),                          &
                               zero))

        ! Now form a half timestep estimate of the proportion
        ! of the depositing volume which contains liquid cloud

    cfltemp(i)=half*MAX(area_mix(i)-(cfl(i)-cfltemp(i)),zero)               &
                       /(area_mix(i)+area_ice_1(i))                         &
               + half*area_mix(i)/ (area_mix(i)+area_ice_1(i))

        ! Recalculate an estimate of the liquid water removed
    dqil(i)=MAX( MIN( dqi_dep(i)*cfltemp(i),                                &
                        qcl(i)*areamix_over_cfliq(i)),zero)

          ! Update liquid cloud fraction and transfer rate
    deltacfl(i) = cfl(i)

    cfl(i) = cfl(i) * SQRT(MAX(one-dqil(i)/qcl(i),                          &
                               zero))

    cfltransfer(i) = cfltransfer(i) + (cfl(i) - deltacfl(i))                &
                       * one_over_tsi

  ELSE
        ! Deposition does not remove any liquid water content
    dqil(i)=zero
  END IF ! cfliq gt 0.0 etc.


END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Adjust ice content
      !-----------------------------------------------
  qcf(i) = qcf(i) + dqi_dep(i) + dqi_sub(i)


      !-----------------------------------------------
      ! Adjust liquid and vapour contents (liquid adjusts first)
      !-----------------------------------------------

  qcl(i) = qcl(i) - dqil(i)  ! Bergeron Findeisen acts first
  t(i) = t(i) + lfrcp * dqil(i)
  dqi(i) = dqi_dep(i) + dqi_sub(i)- dqil(i)

  q(i) = q(i) - dqi(i)
  t(i) = t(i) + lsrcp * dqi(i)

    !-----------------------------------------------
    ! Store depostion/sublimation rate
    !-----------------------------------------------
  ptransfer(i) = ptransfer(i) + (dqi_dep(i) + dqi_sub(i))                   &
                  * one_over_tsi

END DO  

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_deposition
END MODULE lsp_deposition_mod
