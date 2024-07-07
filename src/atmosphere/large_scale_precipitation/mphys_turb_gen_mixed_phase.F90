! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Generation of mixed-phase cloud
! by turbulent processes
MODULE mphys_turb_gen_mixed_phase_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName='MPHYS_TURB_GEN_MIXED_PHASE_MOD'

CONTAINS

SUBROUTINE mphys_turb_gen_mixed_phase( q_work, t_work, qcl_work, qcf_work,     &
                                       q_inc, qcl_inc, cfl_inc,  cf_inc,       &
                                       t_inc,  dqcl_mp, bl_levels,             &
                                       bl_w_var, cff_work, cfl_work, cf_work,  &
                                       q_n, cfl_n, cf_n, p_layer_centres,      &
                                       rhodz_dry, rhodz_moist, deltaz,         &
                                       qcl_mpt, tau_d, inv_prt, disprate,      &
                                       inv_mt, si_avg, dcfl_mp, sigma2_s   )

!Microphysics modules
USE mphys_inputs_mod,      ONLY: l_subgrid_cfl_mp_by_erosion,                 &
                                 l_mixed_phase_t_limit, nbins_mp, mp_t_limit, &
                                 mp_dz_scal, mp_tau_d_lim, mp_czero
USE mphys_constants_mod,   ONLY: cx, constp
USE lsp_moments_mod,       ONLY: lsp_moments
USE mphys_ice_mod,         ONLY: rhoi
USE lsp_dif_mod,           ONLY: air_conductivity0, air_diffusivity0, tcor1,  &
                                 tcor2, cpwr

!General and constants modules
USE gen_phys_inputs_mod,   ONLY: l_mr_physics
USE conversions_mod,       ONLY: pi, zerodegc
USE water_constants_mod,   ONLY: lc, lf
USE planet_constants_mod,  ONLY: cp, r, repsilon, pref, rv, g

! Grid bounds module
USE atm_fields_bounds_mod, ONLY: tdims

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_mix_new     => qsat_mix,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_mphys !Currently defaults to FALSE

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! Purpose:
! Produces mixed phase cloud by subgrid turbulence, in order to
! improve the radiation budget and account for missing microphysical
! processes

! Method and paper reference:
! Uses the stochastic model from Field et al (2013), Q. J. R. Met. Soc
! "Mixed-phase clouds in a turbulent environment. Part 2: Analytic treatment"

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!  Fortran90
!  This code is written to UMDP3 programming standards.

! Documentation: UMDP 29A.

! Subroutine arguments

INTEGER, INTENT(IN) :: bl_levels

REAL, INTENT(INOUT) :: q_work(   tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) ::  qcl_work(tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) ::  qcf_work(tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: t_work(  tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: cff_work( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: cfl_work( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(INOUT) :: cf_work( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                            1 : tdims%k_end )

REAL, INTENT(IN) :: q_n(         tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(IN) :: cfl_n(       tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(IN) :: cf_n(        tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(IN) :: bl_w_var(    tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )

REAL, INTENT(IN) :: p_layer_centres( tdims%i_start : tdims%i_end,             &
                                     tdims%j_start : tdims%j_end,             &
                                                 0 : tdims%k_end )

REAL, INTENT(IN) :: rhodz_dry(  tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                            1 : tdims%k_end )

REAL, INTENT(IN) ::   rhodz_moist(tdims%i_start : tdims%i_end,                &
                                  tdims%j_start : tdims%j_end,                &
                                              1 : tdims%k_end )

REAL, INTENT(IN) ::   deltaz( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: q_inc( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: qcl_inc( tdims%i_start : tdims%i_end,                  &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: cfl_inc( tdims%i_start : tdims%i_end,                  &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: cf_inc( tdims%i_start : tdims%i_end,                   &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: t_inc( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT(INOUT) :: dqcl_mp( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                            1 : tdims%k_end )

! Diagnostics
!                    qcl generate by turbulent mixed_phase
REAL, INTENT(OUT) :: qcl_mpt( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Turbulent decorrelation timescale [s]
REAL, INTENT(OUT) :: tau_d(   tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Inverse Phase-Relaxation Timescale [s-1]
REAL, INTENT(OUT) :: inv_prt( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Diagnosed turbulent dissipation rate [m2 s-3]
REAL, INTENT(OUT) :: disprate(tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Inverse Mixing Timescale [s-1]
REAL, INTENT(OUT) :: inv_mt(  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Mean of subgrid supersaturation PDF
REAL, INTENT(OUT) :: si_avg(  tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Liquid cloud fraction diagnosed by mixed phase scheme
REAL, INTENT(OUT) :: dcfl_mp( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

!                    Variance of subgrid supersaturation PDF
REAL, INTENT(OUT) :: sigma2_s(tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

! Local variables
REAL, PARAMETER :: smallnum = 2.2e-14
                   ! Small value used in if tests; taken to be the same as
                   ! lsp_moments uses for consistency.

REAL, PARAMETER :: one_third = 1.0/3.0

REAL :: Ei               ! Saturated vapour pressure wrt ice [Pa]
REAL :: rhice            ! Relative humidity of ice [0-1]
REAL :: siw              ! The value of ice supersaturation at
                         ! water saturation []
REAL :: siw_lim          ! Limit of ice supersaturation
REAL :: aa               ! Thermodynamic term (for definition see Field (2013))
REAL :: b0               ! Thermodynamic term (for definition see Field (2013))
REAL :: Ai               ! Thermodynamic term (for definition see Field (2013))
REAL :: bi               ! Thermodynamic term (for definition see Field (2013))
REAL :: ka               ! temperature-corrected conductivity
REAL :: dv               ! temperature and pressure-corrected diffusivity
REAL :: fdist            ! frequency distribution of
REAL :: qv_excess        ! Excess moisture mixing ratio
REAL :: Sice             ! Saturation
REAL :: deltas           ! Local change of sigma_s
REAL :: fac              ! Local factors ...
REAL :: fac2             ! ... used in calculations
REAL :: four_root_sigmas ! local variable 4.0 * SQRT(sigma_s)
REAL :: dz_scal          ! Scaling factor applied to layer thickness
REAL :: t_limit          ! local temperature limit applied to mixed phase
                         ! calculation
REAL :: t_corr           ! temperature correction for constants
REAL :: p_corr           ! pressure correction for diffusivity
REAL :: tau_d_work       ! local working value of the turbulent decorrelation
                         ! timescale

INTEGER :: grd_pts       ! Number of points in decomposed domain

INTEGER :: ibin          ! Bin number

INTEGER :: i             ! Loop counter in x direction
INTEGER :: j             ! Loop counter in y direction
INTEGER :: k             ! Loop counter in z direction

REAL :: t_local2d(tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end)
! Local value of temperature (K)

REAL :: q_local2d(tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end)
! Local value of humidity (kg/kg)

REAL :: qsi_2d( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end )

REAL :: qsw_2d( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end )

REAL :: mom1(   tdims%i_start : tdims%i_end)
! First moment of the ice particle size distribution

REAL :: rho_air( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end )
! Air density (kg m-3)

REAL :: rho_dry( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end )
! Dry air density (kg m-3)

REAL :: cff_inv(  tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end,                                &
                              1 : tdims%k_end )
! 1/frozen cloud fraction

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MPHYS_TURB_GEN_MIXED_PHASE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=============================================================================
! START OF PHYSICS
!=============================================================================

! Set the temperature limit

IF (l_mixed_phase_t_limit) THEN
  t_limit = mp_t_limit
  siw_lim = -1.0
ELSE
  t_limit = zerodegc - 0.01
  siw_lim = smallnum
END IF

dqcl_mp(:,:,:)  =  0.0
qcl_mpt(:,:,:)  =  0.0
tau_d(:,:,:)    = -2.0 ! To isolate points where the scheme is not used.
inv_prt(:,:,:)  =  0.0
disprate(:,:,:) =  0.0
inv_mt(:,:,:)   =  0.0
si_avg(:,:,:)   =  0.0
dcfl_mp(:,:,:)  =  0.0
sigma2_s(:,:,:) =  0.0

grd_pts =  tdims%j_len  *                                                     &
           tdims%i_len

! Loop over boundary layer levels
DO k = 1, bl_levels-1

  ! Store and modify values of q and T for mixed phase calculation
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      q_local2d(i,j) = q_work(i,j,k)
      t_local2d(i,j) = t_work(i,j,k)

      ! dry air density
      rho_dry(i,j) = rhodz_dry(i,j,k) / deltaz(i,j,k)

      IF (l_mr_physics) THEN
         ! rho is the dry density
        rho_air(i,j) = rho_dry(i,j)
      ELSE
         ! rho is the moist density
        rho_air(i,j) = rhodz_moist(i,j,k) / deltaz(i,j,k)

        ! Subsequent liquid cloud calculation uses mixing
        ! ratios. Convert moisture variable to a mixing ratio
        q_local2d(i,j) = rho_air(i,j)*q_local2d(i,j)/rho_dry(i,j)

      END IF

      cff_inv(i,j,k) = 1.0 / MAX( cff_work(i,j,k), 0.001 )

    END DO
  END DO

  !  Calls to qsat_mix, qsat_wat_mix and lsp_moments.
  !  Done level-by-level for efficiency and generates the values
  !  of mom1, qsi_2d and qsw_2d

  ! Always request saturated values as mixing ratios
  ! for subsequent calculation of liquid cloud
  IF (l_new_qsat_mphys) THEN
    CALL qsat_mix_new(qsi_2d, t_local2d, p_layer_centres(:,:,k),              &
                      tdims%j_len,tdims%i_len)

    CALL qsat_wat_mix_new(qsw_2d, t_local2d, p_layer_centres(:,:,k),          &
                          tdims%j_len,tdims%i_len)
  ELSE
    CALL qsat_mix(qsi_2d, t_local2d, p_layer_centres(1,1,k), grd_pts,         &
                  .TRUE. )

    CALL qsat_wat_mix(qsw_2d, t_local2d, p_layer_centres(1,1,k), grd_pts,     &
                      .TRUE. )
  END IF

  ! Assuming the Field 'Generic' PSD and no ventilation effects here,
  ! to keep things simple.
  ! Note: this may not be consistent with assumptions in the microphysical
  ! deposition rate calculation in lsp_deposition

  DO j = tdims%j_start, tdims%j_end

    !Call lsp_moments row by row to help it play nicely with the INTERFACE
    CALL lsp_moments( tdims%i_len, rho_air(:,j), t_local2d(:,j),              &
                      qcf_work(:,j,k), cff_inv(:,j,k), cx(84), mom1 )

    DO i = tdims%i_start, tdims%i_end

       ! First statement in IF test catches any -ve qv from pc2
      IF (q_local2d(i,j)    > smallnum .AND.                                  &
          t_local2d(i,j)    < t_limit  .AND.                                  &
          bl_w_var(i,j,k)   > 1e-12 ) THEN

          ! N.B. Using approximate mixing ratio of
          ! vapour to get an RH
        rhice = q_local2d(i,j) / qsi_2d(i,j)

        ! N.B. This way of obtaining Siw and Ei
        ! is for qsw, qsi mixing ratios
        siw   = (qsw_2d(i,j) / qsi_2d(i,j)) - 1.0

        ei    = qsi_2d(i,j) * p_layer_centres(i,j,k) /                        &
              ( qsi_2d(i,j) + repsilon )

        ! Pressure and temperature corrections:
        t_corr = ( (t_local2d(i,j) / zerodegc)**cpwr) *                       &
                 ( tcor1 / ( t_local2d(i,j) + tcor2 ) )

        p_corr = ( pref / p_layer_centres(i,j,k) )

        dv = air_diffusivity0 * t_corr * p_corr

        ka = air_conductivity0 * t_corr

        bi = 1.0 / q_local2d(i,j) + (Lc+Lf)**2 /                              &
             (cp * rv * t_local2d(i,j) ** 2)

        Ai = 1.0 / (rhoi *(Lc+Lf)**2 / (ka*rv*t_local2d(i,j)**2) +            &
                  rhoi * rv * t_local2d(i,j) / (Ei*Dv))

        b0 = 4.0 * pi * constp(35) * rhoi * Ai / rho_air(i,j)

        aa = ( g / (r*t_local2d(i,j) ) *                                      &
             ( (Lc+Lf)*r / (cp*rv*t_local2d(i,j))-1.0))

        dz_scal = mp_dz_scal * deltaz(i,j,k)

        disprate(i,j,k) = 2.0 * bl_w_var(i,j,k)**1.5 /                        &
                          (mp_czero * dz_scal)

        ! Limit tau_d to specified limit
        tau_d_work = 2.0 * bl_w_var(i,j,k) / (disprate(i,j,k)*mp_czero)

        inv_prt(i,j,k)  = bi * b0 * mom1(i)

        inv_mt(i,j,k) = (disprate(i,j,k) / dz_scal**2)**one_third

        ! Add factor to speed up calculations
        fac = 1.0 / (inv_prt(i,j,k) + inv_mt(i,j,k))

        sigma2_s(i,j,k) = 0.5 *aa**2 * SQRT(bl_w_var(i,j,k))                  &
                        * dz_scal * fac

        si_avg(i,j,k) = (rhice-1.0) * inv_mt(i,j,k) * fac

        IF ( siw > siw_lim .AND. aa > smallnum .AND.                          &
             tau_d_work <=  mp_tau_d_lim   ) THEN

          tau_d(i,j,k) = tau_d_work

          four_root_sigmas = 4 * SQRT(sigma2_s(i,j,k))
          fac2   = SQRT(2 * pi * sigma2_s(i,j,k) )
          deltas = four_root_sigmas / nbins_mp

          DO ibin = 0, nbins_mp-1

            sice      = siw + ibin * deltas
            qv_excess = qsi_2d(i,j) * ibin * deltas
            fdist     = 0.0

            IF ( ABS(sice-si_avg(i,j,k)) <= four_root_sigmas ) THEN

              fdist = EXP(-(sice-si_avg(i,j,k))**2/(2*sigma2_s(i,j,k)))       &
                      / fac2

            END IF

            qcl_mpt(i,j,k) = qcl_mpt(i,j,k) + qv_excess * fdist * deltas
            dcfl_mp(i,j,k) = dcfl_mp(i,j,k) + fdist * deltas

          END DO ! ibin

          ! The calculated qcl_mpt is a mixing ratio.
          ! Convert to specific quantity if necessary
          IF (.NOT. l_mr_physics) THEN
            qcl_mpt(i,j,k) = rho_dry(i,j)*qcl_mpt(i,j,k)/rho_air(i,j)
          END IF

          IF ( cfl_work(i,j,k) < smallnum .OR.                                &
              .NOT. l_subgrid_cfl_mp_by_erosion ) THEN

            ! Ensure that qcl_inc + qcl_mpt < q_n
            qcl_mpt(i,j,k) = MIN( qcl_mpt(i,j,k),                             &
                                  q_n(i,j,k) - qcl_inc(i,j,k) )

            qcl_inc(i,j,k) = qcl_inc(i,j,k) + qcl_mpt(i,j,k)
            q_inc(i,j,k)   = q_inc(i,j,k)   - qcl_mpt(i,j,k)
            t_inc(i,j,k)   = t_inc(i,j,k)   + Lc * qcl_mpt(i,j,k) / cp

            cfl_inc(i,j,k) = MIN( cfl_inc(i,j,k) + dcfl_mp(i,j,k),            &
                                       1.0 - cfl_n(i,j,k)       )

            cf_inc(i,j,k)  = MIN( cf_inc(i,j,k) + dcfl_mp(i,j,k),             &
                                       1.0 - cf_n(i,j,k))

            qcl_work(i,j,k) = qcl_work(i,j,k) + qcl_mpt(i,j,k)
            q_work(i,j,k)   = q_work(i,j,k)   - qcl_mpt(i,j,k)
            t_work(i,j,k)   = t_work(i,j,k)   + Lc * qcl_mpt(i,j,k) / cp

            cfl_work(i,j,k) = MIN( cfl_work(i,j,k) + dcfl_mp(i,j,k), 1.0)

            cf_work(i,j,k)  = MIN( cf_work(i,j,k) + dcfl_mp(i,j,k), 1.0)

          ELSE ! cfl_work == 0 etc

            dqcl_mp(i,j,k) = qcl_mpt(i,j,k) - qcl_work(i,j,k)

          END IF  ! cfl_work == 0 etc

        ELSE IF (tau_d_work >  mp_tau_d_lim) THEN

          tau_d(i,j,k) = -1.0 ! Flag for tau_d limit exceeded

        END IF    ! Siw > 0 etc

      END IF       ! qlocal2d > 0

    END DO         ! i
  END DO           ! j
END DO             ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE mphys_turb_gen_mixed_phase
END MODULE mphys_turb_gen_mixed_phase_mod
