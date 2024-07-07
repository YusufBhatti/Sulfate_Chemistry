! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme.
!
!  Contains the following subroutines:
!
!  lsp_riming (Riming of ice particles including shape-dependence)
!  lsp_riming_sphere  (Riming of spherical ice particles. Used for Graupel)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation
!
! Subroutine Interface:
MODULE lsp_riming_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_RIMING_MOD'

CONTAINS

! Subroutine lsp_riming only called if l_shape_rime is TRUE
SUBROUTINE lsp_riming(                                                        &
  points, timestep,                                                           &
                                          ! Number of points and tstep
  qcl, qcf, t,                                                                &
                                          ! Water contents and temp
  area_liq, area_mix, cfliq, cficei,                                          &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl, cff                      ! Current cloud fractions for
!                                         ! updating
  rho, m0, tcg, tcgi, corr,                                                   &
                                          ! Parametrization information
  lfrcp , ice_type,                                                           &
                                          ! Microphysical information
  l_psd,                                                                      &
                                          ! Code options
  ptransfer,                                                                  &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                               &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt                                                                &
!    &, cftransfer,cfltransfer,cfftransfer! Cloud transfer diagnostics
  )

! Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod, ONLY: cx, constp,qclmin_rime, area_ratio_prefac,              &
                               area_ratio_expn, zerodegc,                     &
                               zero, one

  ! Microphysics modules
USE mphys_constants_mod, ONLY: ice_type_offset
USE mphys_inputs_mod,    ONLY: l_diff_icevt

USE pc2_constants_mod,   ONLY: cloud_pc2_tol_2

USE science_fixes_mod,   ONLY: l_fix_riming

USE um_types,            ONLY: real_lsprec

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

! Large scale precipitation modules
USE lsp_moments_mod, ONLY: lsp_moments
IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of ice particle riming

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.
! This code allows for possible dependence of the the riming rate
! on the cross-sectional area of the ice/snow crystals.
! The parametrization of area in terms of particle size from Heymsfield and
! Miloshevich (J. Atmos. Sci., 60, 936 (2003)) is used.
! Also included is a minimum qcl, above which liquid water
! is available for riming. This follows the findings of
! Marimaya (J. Meteor. Japan, 53(6), 384-392) that
! there is a minimum droplet size for riming to occur.

!
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
  cfliq(points),                                                              &
                        ! Fraction of gridbox with liquid cloud
  cficei(points),                                                             &
                        ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
    rho(points),                                                              &
                          ! Air density / kg m-3
    m0,                                                                       &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                              &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                             &
                          ! 1/tcg (no units)
    corr(points),                                                             &
                          ! Fall velocity correction factor (no units)
    lfrcp,                                                                    &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi
                          ! 1/(timestep*iterations)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  qcl(points),                                                                &
                        ! Liquid water content / kg kg-1
  qcf(points),                                                                &
                        ! Ice water content    / kg kg-1
  t(points),                                                                  &
                        ! Temperature / K
  ptransfer(points) ! Mass rimed in this timestep / kg kg-1

!      Real, Intent(InOut) ::
!    &, cf_transfer_rate(points) ! Cloud fraction increment this tstep
!    &, cfl_transfer_rate(points)! Liquid cloud fraction inc this tstep
!    &, cff_transfer_rate(points)! Ice cloud fraction inc this tstep

LOGICAL, INTENT(IN) ::                                                        &
  l_psd,                                                                      &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines fallspeed branch to use

! Local Variables

INTEGER ::                                                                    &
  i,                                                                          &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

REAL (KIND=real_lsprec) ::                                                    &
  qclnew(points),                                                             &
                        ! For value of qcl after riming  / kg kg-1
  dqi(points),                                                                &
                        ! Amount of ice rimed  / kg kg-1
  m_2_di(points),                                                             &
                        ! riming moment of particle size distribution
  m_2_dic(points),                                                            &
                        ! riming moment of particle size distribution
  ic_qcl(points),                                                             &
                        ! In-cloud liquid water (qcl)
  qcl_rime,                                                                   &
                        ! LWC available for riming
  rime_eff_agg,                                                               &
                        ! Shape-dependent riming efficency of aggregates
  rime_eff_cry,                                                               &
                        ! Shape-dependent riming efficency of crystals
  rime_expn_agg,                                                              &
                        ! Order of riming moment for aggregate,
                        ! modified to account for cross-sectional area
  rime_expn_cry,                                                              &
                        ! Order of riming moment for crystals,
                        ! modified to account for cross-sectional area
  rime_fac
                        ! Factor in the implicit solution used for
                        ! qcl after riming

LOGICAL :: l_riming_conditions_met
REAL (KIND=real_lsprec) :: qcl_min_threshold
  ! Minimum threshold of qcl for riming

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_RIMING'

    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set minimum threshold for qcl. This needs to be such that the in-cloud
! qcl ( i.e. qcl / cloud fraction) will be greater than qclmin_rime.
! To obtain this, multiply the minimum in-cloud value by the minimum
! acceptable cloud fraction
qcl_min_threshold = cloud_pc2_tol_2 * qclmin_rime

! Set other properties
cry_offset        = ice_type * ice_type_offset
rime_expn_agg     = cx(81) + area_ratio_expn
rime_expn_cry     = cx(181) + area_ratio_expn
rime_eff_agg      = area_ratio_prefac
rime_eff_cry      = area_ratio_prefac

IF (l_psd) THEN
      ! Use the generic ice particle size distribution
      ! Calculate the 2+di (cx(81)) moment of the
      ! ice particle size distribution.

  IF (.NOT. l_diff_icevt) THEN
        ! Only one set of vt-D parameters
    CALL lsp_moments(points,rho,t,qcf,cficei,rime_expn_agg,m_2_di)
  ELSE

    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    CALL lsp_moments(points,rho,t,qcf,cficei,rime_expn_cry,m_2_dic)
                      ! psd moment with crystal parameters
    CALL lsp_moments(points,rho,t,qcf,cficei,rime_expn_agg,m_2_di)
                      ! psd moment with aggregate parameters
  END IF
END IF

! May be required for 32-bit PROC comparability for some configs
! m_2_di(:)  = MAX(m_2_di,EPSILON(one))
! m_2_dic(:) = MAX(m_2_dic,EPSILON(one))

DO i = 1, points

  ! Set riming conditions to False and then check if they are being met
  ! due to a different set of conditions. The first set should not produce
  ! very small values of qcl, the second set should.
  l_riming_conditions_met = .FALSE.

  IF ( l_fix_riming .AND. qcf(i) > m0 .AND. qcl(i) > qcl_min_threshold .AND.  &
       t(i) < zerodegc .AND. area_mix(i) > cloud_pc2_tol_2             .AND.  &
       cfliq(i) > cloud_pc2_tol_2                                     ) THEN
    ! Only rime if the in-cloud qcl is bigger than the threshold set
    ! and cloud fractions are bigger than the minimum tolerance set by the
    ! cloud scheme. By definition qcl(i) must be positive if ic_qcl(i) is
    ! above qcl_rime and the cloud fractions are positive and above the
    ! tolerances.

    l_riming_conditions_met = .TRUE.

  ELSE IF ( qcf(i) > m0 .AND. qcl(i) > 0.0 .AND. t(i) < zerodegc              &
           .AND. area_mix(i) > 0.0 .AND. cfliq(i) > 0.0 ) THEN

    ! l_fix_riming turned off: rime only if original conditions have been met
    ! which may lead to small numbers being passed through the model

    l_riming_conditions_met = .TRUE.

  END IF

  IF ( l_riming_conditions_met ) THEN

    ! Determine what the in-cloud qcl value is - used in a few places
    ic_qcl(i) = qcl(i) / cfliq(i)

        !-----------------------------------------------
        ! Calculate water content of mixed phase region
        !-----------------------------------------------
    IF (l_psd) THEN

          ! Calculate the riming rate using the generic PSD
      IF (.NOT. l_diff_icevt) THEN
         ! Only one set of vt-D parameters
         ! constp(81) = (pi/4) ci
        rime_fac = constp(81)                                                 &
                   *corr(i)*timestep*rime_eff_agg*m_2_di(i)

      ELSE
        IF (l_use_agg_vt(i)) THEN
           ! Use aggregate parameters
          rime_fac = constp(81)                                               &
                   *corr(i)*timestep*rime_eff_agg*m_2_di(i)

        ELSE
           ! Use crystal parameters
          rime_fac = constp(181)                                              &
                   *corr(i)*timestep*rime_eff_cry*m_2_dic(i)

        END IF
      END IF
    ELSE
          ! Use the defined gamma distribution

      rime_fac = rime_eff_agg*constp(9+cry_offset)*tcg(i)                     &
                *corr(i)*timestep*(rho(i)*qcf(i)*cficei(i)                    &
                *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset)

    END IF

    ! Only the fraction of qcl above qclmin is
    ! available for riming.
    ! Note that if in-liquid-cloud qcl < qclmin_rime
    ! then it is unchanged by riming, i.e., qclnew=qcl/cfliq
    !
      qcl_rime  = ic_qcl(i) + rime_fac * qclmin_rime

      qclnew(i) = qcl_rime / (one + rime_fac)

      qclnew(i) = MIN(ic_qcl(i), qclnew(i))

        !-----------------------------------------------
        ! Convert to new grid box total water content
        !-----------------------------------------------
        ! The mixed-phase cloud now contains the result
        ! of the riming process, which may be the same
        ! qcl as before riming if qclmin is not reached
    qclnew(i) = qcl(i)*area_liq(i)/cfliq(i)+qclnew(i)*area_mix(i)
    dqi(i)    = ( qcl(i) - qclnew(i) )


        !-----------------------------------------------
        ! Store process rate / kg kg-1 s-1 and cloud fraction changes
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
    qcf(i) = qcf(i) + dqi(i)
    t(i)   = t(i) + lfrcp * dqi(i)
    qcl(i) = qclnew(i)


  END IF ! qcf(i) >  m0 etc.

END DO  ! Points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_riming


! Subroutine lsp_riming only called if l_shape_rime is FALSE
SUBROUTINE lsp_riming_sphere(                                                 &
  points, timestep,                                                           &
                                          ! Number of points and tstep
  qcl, qcf, t,                                                                &
                                          ! Water contents and temp
  area_liq, area_mix, cfliq, cficei,                                          &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl, cff                      ! Current cloud fractions for
!                                         ! updating
  rho, m0, tcg, tcgi, corr,                                                   &
                                          ! Parametrization information
  lfrcp , ice_type,                                                           &
                                          ! Microphysical information
  l_psd,                                                                      &
                                          ! Code options
  ptransfer,                                                                  &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                               &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt                                                                &
!    &, cftransfer,cfltransfer,cfftransfer! Cloud transfer diagnostics
  )

USE lsprec_mod,         ONLY: cx, constp, zerodegc, zero

  ! Microphysics modules
USE mphys_constants_mod, ONLY: ice_type_offset
USE mphys_inputs_mod,    ONLY: l_diff_icevt

! General atmosphere modules

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

! Dr Hook modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

! Large scale precipitation modules
USE lsp_moments_mod, ONLY: lsp_moments
IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of ice particle riming
!   assuming particles are spherical

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.
! In this rountine the ice particles are assumed to be spherical,
! which is appropriate when accumulation of cloud droplets by
! graupel calculated
!

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
  cfliq(points),                                                              &
                        ! Fraction of gridbox with liquid cloud
  cficei(points),                                                             &
                        ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
    rho(points),                                                              &
                          ! Air density / kg m-3
    m0,                                                                       &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                              &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                             &
                          ! 1/tcg (no units)
    corr(points),                                                             &
                          ! Fall velocity correction factor (no units)
    lfrcp,                                                                    &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi
                          ! 1/(timestep*iterations)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  qcl(points),                                                                &
                        ! Liquid water content / kg kg-1
  qcf(points),                                                                &
                        ! Ice water content    / kg kg-1
  t(points),                                                                  &
                        ! Temperature / K
  ptransfer(points) ! Mass rimed in this timestep / kg kg-1

!      Real, Intent(InOut) ::
!    &, cf_transfer_rate(points) ! Cloud fraction increment this tstep
!    &, cfl_transfer_rate(points)! Liquid cloud fraction inc this tstep
!    &, cff_transfer_rate(points)! Ice cloud fraction inc this tstep

LOGICAL, INTENT(IN) ::                                                        &
  l_psd,                                                                      &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines fallspeed branch to use

! Local Variables

INTEGER ::                                                                    &
  i,                                                                          &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

REAL (KIND=real_lsprec) ::                                                    &
  qclnew(points),                                                             &
                        ! For value of qcl after riming  / kg kg-1
  dqi(points),                                                                &
                        ! Amount of ice rimed  / kg kg-1
  m_2_di(points),                                                             &
                        ! 2+DI moment of particle size distribution
  m_2_dic(points)
                        ! 2+DIC moment of particle size distribution

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_RIMING_SPHERE'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cry_offset=ice_type*ice_type_offset

IF (l_psd) THEN
      ! Use the generic ice particle size distribution
      ! Calculate the 2+di (cx(81)) moment of the
      ! ice particle size distribution.

  IF (.NOT. l_diff_icevt) THEN
        ! Only one set of vt-D parameters

    CALL lsp_moments(points,rho,t,qcf,cficei,                                 &
                            cx(81),m_2_di)
  ELSE
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    CALL lsp_moments(points,rho,t,qcf,cficei,                                 &
                            cx(181),m_2_dic)
                      ! psd moment with crystal parameters
    CALL lsp_moments(points,rho,t,qcf,cficei,                                 &
                            cx(81),m_2_di)
                      ! psd moment with aggregate parameters
  END IF
END IF

DO i = 1, points

  IF (qcf(i) >  m0 .AND. qcl(i) >  zero .AND. t(i)  <   zerodegc              &
      .AND. area_mix(i) >  zero .AND. cfliq(i) >  zero) THEN

        !-----------------------------------------------
        ! Calculate water content of mixed phase region
        !-----------------------------------------------
    IF (l_psd) THEN

          ! Calculate the riming rate using the generic PSD
      IF (.NOT. l_diff_icevt) THEN
         ! Only one set of vt-D parameters
         ! constp(81) = (pi/4) ci
        qclnew(i) = qcl(i) /                                                  &
                  (cfliq(i)+cfliq(i)*constp(81)                               &
                   *corr(i)*timestep*m_2_di(i))
      ELSE
        IF (l_use_agg_vt(i)) THEN
           ! Use aggregate parameters
          qclnew(i) = qcl(i) /                                                &
                  (cfliq(i)+cfliq(i)*constp(81)                               &
                   *corr(i)*timestep*m_2_di(i))
        ELSE
           ! Use crystal parameters
          qclnew(i) = qcl(i) /                                                &
                  (cfliq(i)+cfliq(i)*constp(181)                              &
                   *corr(i)*timestep*m_2_dic(i))
        END IF
      END IF
    ELSE
          ! Use the defined gamma distribution

      qclnew(i) = qcl(i) /                                                    &
                (cfliq(i)+cfliq(i)*constp(9+cry_offset)*tcg(i)                &
                *corr(i)*timestep*(rho(i)*qcf(i)*cficei(i)                    &
                *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset))
    END IF

        !-----------------------------------------------
        ! Convert to new grid box total water content
        !-----------------------------------------------
    qclnew(i)=qcl(i)*area_liq(i)/cfliq(i)+qclnew(i)*area_mix(i)
    dqi(i)=(qcl(i)-qclnew(i))


        !-----------------------------------------------
        ! Store process rate / kg kg-1 s-1 and cloud fraction changes
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

    !         There are no cloud fraction updates associated
    !         with the riming term. They will have already been set to
    !         zero on input to this subroutine, which is why these
    !         are commented out.
    !          cf_transfer_rate(i)  = zero / (timestep*iterations)
    !          cfl_transfer_rate(i) = zero / (timestep*iterations)
    !          cff_transfer_rate(i) = zero / (timestep*iterations)

              !-----------------------------------------------
              ! Update water contents
              !-----------------------------------------------
    qcf(i) = qcf(i) + dqi(i)
    t(i)   = t(i) + lfrcp * dqi(i)
    qcl(i) = qclnew(i)

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------
    !          IF (ice_type /= 3) THEN
    !            These are commented out since there is currently no
    !            cloud fraction update associated with the riming term.

    !            cf(i)  = cf(i) + cf_transfer_rate(i) *timestep*iterations
    !            cfl(i) = cfl(i)+ cfl_transfer_rate(i)*timestep*iterations
    !            cff(i) = cff(i)+ cff_transfer_rate(i)*timestep*iterations

    !            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
    !            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)
    !            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)

    !          END IF  ! ice_type /= 3

  END IF ! qcf(i) >  m0 etc.

END DO  ! Points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_riming_sphere
END MODULE lsp_riming_mod
