! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_ideal_1d_profiles_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_1D_PROFILES_MOD'

CONTAINS

! Subroutine rcf_ideal_1d_profiles
!
! Description:
!   Calculate 1D vertical profiles for exner, (dry) theta and (dry) rho.
!   This subroutine acts as an interface to rcf_1d_newton.
!
! Method:
!   Solve discrete equations for hydrostatic balance and equation of state,
!   for given temperature or potential temperature profile. Equations are
!   solved using the Newton-Raphson method.
!
!   Returns potential temperature, exner pressure and density.
!
!   Generates thermodynamic profiles at height levels given by:
!      z = (coord - surface_offset) * depth_factor
!   for coord = xi3_at_theta and xi3_at_rho.
!
!   This routine is used in two ways:
!     (i)  Generate semi-implicit reference profiles with:
!          xi3_at_theta   = eta_theta_levels
!          xi3_at_rho     = eta_rho_levels
!          depth_factor   = height_domain
!          surface_offset = 0.0
!
!     (ii) Generate profiles on a column of the grid:
!          xi3_at_theta   = xi3_at_theta
!          xi3_at_rho     = xi3_at_rho
!          depth_factor   = 1.0
!          surface_offset = planet_radius
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_1d_profiles (exner_profile_ext, rho_profile,             &
                                  theta_profile, tprofile_number,             &
                                  t_surface, p_surface,                       &
                                  xi2_p, xi3_at_theta, xi3_at_rho,            &
                                  depth_factor, surface_offset,               &
                                  dtheta_dz1, height_dz1, l_constant_dz,      &
                                  g_theta_profile, model_levels,              &
                                  print_profile)

USE rcf_interp_weights_mod, ONLY: &
    intw_rho2w,                   &
    intw_w2rho

USE planet_constants_mod, ONLY: &
    cp,                         &
    g,                          &
    kappa,                      &
    planet_radius,              &
    omega,                      &
    p_zero,                     &
    r,                          &
    recip_kappa

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_dthetadz,                  &
    tp_isothermal,                &
    tp_HD209458b_Heng_mid,        &
    tp_HD209458b_Heng_smooth_mid, &
    tp_HD209458b_iro,             &
    tp_GJ1214b,                    &
    tp_Y_Dwarf,                   &
    tp_file,                      &
    tp_HD209458b_Heng_smooth_day, &
    tp_BRuntV,                    &
    tp_BV_isoth,                  &
    tp_dump

USE tforce_mod, ONLY:         &
    tf_HD209458b_Heng,        &
    tf_HD209458b_Heng_smooth, &
    tf_HD209458b_iro,         &
    tf_GJ1214b,               &
    tf_Y_Dwarf,               &
    tf_file

USE rcf_nlist_recon_idealised_mod, ONLY: &
    initial_tolerance,                   &
    k_const,                             &
    profile_filename

USE idealise_run_mod, ONLY: &
    l_shallow

USE rcf_ideal_baro_eta_conv_mod, ONLY: &
    baro_eta_conv

USE rcf_ideal_baro_T_p_mod, ONLY: &
    baro_t_p

USE rcf_ideal_newton_solver_mod, ONLY: &
    newton_solver

USE rcf_ideal_file_forcing_mod, ONLY: &
    read_file_temp

USE rcf_ideal_deep_baroclinic_constants_mod, ONLY: &
    tau1,                                          &
    tau2,                                          &
    tau1_int,                                      &
    tau2_int

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    PrStatus_Normal

USE ereport_mod, ONLY: &
    ereport

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN) :: tprofile_number
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN)    :: xi2_p
REAL, INTENT(IN)    :: xi3_at_theta(0:model_levels+1)
REAL, INTENT(IN)    :: xi3_at_rho(0:model_levels+1)
REAL, INTENT(IN)    :: depth_factor
REAL, INTENT(IN)    :: surface_offset
REAL, INTENT(IN)    :: dtheta_dz1(3)
REAL, INTENT(IN)    :: height_dz1(2)
REAL, INTENT(IN)    :: g_theta_profile(0:model_levels)
LOGICAL, INTENT(IN) :: l_constant_dz

REAL, INTENT(INOUT) :: t_surface
REAL, INTENT(INOUT) :: p_surface

REAL, INTENT(OUT)   :: exner_profile_ext(0:model_levels+1) ! Additional 0 level
REAL, INTENT(OUT)   :: rho_profile(model_levels)
REAL, INTENT(OUT)   :: theta_profile(0:model_levels)

LOGICAL, OPTIONAL, INTENT(IN) :: print_profile

! Local variables:
INTEGER :: k
INTEGER :: profile_type
INTEGER :: tforce_local   ! Initial profile
INTEGER :: icode          ! error reporting

REAL :: t0_profile(0:model_levels)
REAL :: geopotential_factor(0:model_levels)
REAL :: levs(0:model_levels)
REAL :: z_rho(0:model_levels)   ! Height of rho levels above surface
REAL :: z_theta(0:model_levels) ! Height of theta levels above surface
REAL :: u_term(1:model_levels)
REAL :: theta_surface           ! Potential temperature at lower boundary
REAL :: theta_height1
REAL :: theta_height2
REAL :: t_theta(0:model_levels)      ! Temperature on theta levels
REAL :: exner_theta(0:model_levels)  ! Exner on theta levels
REAL :: extrap1                 ! } Weights for extrapolating exner
REAL :: extrap2                 ! } to model_levels+1
REAL :: exner_top               ! exner at top boundary
REAL :: eta(0:model_levels)
REAL :: r_height                ! For deep baroclinic wave calcs
REAL :: tau2_term               ! For deep baroclinic wave calcs
REAL :: u_temp

LOGICAL :: local_print_profile  ! Local value of print_profile
LOGICAL :: profile_available    ! T-p profile is/isn't for printing

CHARACTER(LEN=*), PARAMETER :: routinename = 'RCF_IDEAL_1D_PROFILES'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Process optional argument:
IF (PRESENT(print_profile)) THEN
  local_print_profile = print_profile
ELSE
  local_print_profile = .FALSE.
END IF

! Make initial guesses:
exner_profile_ext(0) = (p_surface / p_zero)**kappa
DO k = 1, model_levels
  exner_profile_ext(k) = exner_profile_ext(0)
  rho_profile(k)   = 1.0
  u_term(k)        = 0.0
END DO
exner_profile_ext(model_levels+1) = exner_profile_ext(0)

DO k = 0, model_levels
  theta_profile(k) = t_surface / exner_profile_ext(0)
  t0_profile(k)    = t_surface

! When using eta-levels for generating reference profile, the
! height_domain factor is carried with the gravitational acceleration.
  geopotential_factor(k) = depth_factor*g_theta_profile(k)
END DO

! Put level data into levs array:
!   for atmospheric state these are z levels in metres.
!   for reference state these are eta levels.

levs(0) = xi3_at_theta(0) - surface_offset
DO k = 1, model_levels
  levs(k) = xi3_at_rho(k) - surface_offset
END DO

z_rho(0) = 0.0
DO k = 1, model_levels
  z_rho(k) = (xi3_at_rho(k) - surface_offset)*depth_factor
END DO
DO k = 0, model_levels
  z_theta(k) = (xi3_at_theta(k) - surface_offset)*depth_factor
END DO

profile_available = .FALSE.

! Check if tprofile_number is valid and set up t0_profile if required.
SELECT CASE(tprofile_number)

CASE(tp_dthetadz)

  ! t0_profile is a theta point
  profile_type = 1    ! theta specified
  theta_surface = t_surface / exner_profile_ext(0)

  IF (l_constant_dz) THEN
    DO k = 0, model_levels
      t0_profile(k) = theta_surface + dtheta_dz1(1) + z_theta(k)
    END DO
  ELSE
    theta_height1 = theta_surface + dtheta_dz1(1)*height_dz1(1)
    theta_height2 = theta_height1 + dtheta_dz1(2) *              &
                                  (height_dz1(2) - height_dz1(1))
    DO k = 0, model_levels
      IF (z_theta(k) <= height_dz1(1)) THEN
        t0_profile(k) = theta_surface + dtheta_dz1(1) * z_theta(k)
      ELSE IF (z_theta(k) <= height_dz1(2)) THEN
        t0_profile(k) = theta_height1 + dtheta_dz1(2) *          &
                                      (z_theta(k) - height_dz1(1))
      ELSE
        t0_profile(k) = theta_height2 + dtheta_dz1(3) *          &
                                      (z_theta(k) - height_dz1(2))
      END IF
    END DO
  END IF

CASE(tp_isothermal)

  ! t0_profile is a rho point
  profile_type = 2   ! T specified
  t0_profile(0) = t_surface
  DO k = 1, model_levels
    t0_profile(k) = t0_profile(k-1) - dtheta_dz1(1) * depth_factor *    &
                                      (levs(k)-levs(k-1))
  END DO

  DO k = 1, model_levels
    exner_profile_ext(k) = exner_profile_ext(0)*EXP(-g*z_rho(k)/(cp*t_surface))
    rho_profile(k) = p_zero*exner_profile_ext(k)**(recip_kappa)/(r*t_surface)
  END DO

  ! Fill theta levels independently
  theta_profile(0) = t_surface/exner_profile_ext(0)
  exner_theta(0)   = exner_profile_ext(0)
  DO k = 1, model_levels
    theta_profile(k) = theta_profile(0)**EXP(g*z_theta(k)/(cp*t_surface))
    ! u_term needs adding in later.
    exner_theta(k) = exner_theta(0)*EXP(-g*z_theta(k)/(cp*t_surface))
  END DO

  ! Exponential extrapolation to model top
  extrap1 = (z_theta(model_levels)-z_rho(model_levels-2)) /                  &
       (z_rho(model_levels-1)-z_rho(model_levels-2))
  extrap2 = (z_theta(model_levels)-z_rho(model_levels-1)) /                  &
       (z_rho(model_levels-1)-z_rho(model_levels-2))
  exner_top = exner_profile_ext(model_levels)**extrap1 *                     &
       exner_profile_ext(model_levels-1)**(-extrap2)

  ! Set exner_profile_ext(model_levels+1) so that extrapolated exner_top is
  ! obtained using 0.5*(exner_profile_ext(model_levels+1) +
  !                     exner_profile_ext(model_levels))
  exner_profile_ext(model_levels+1) = 2.0*exner_top -                        &
                                       exner_profile_ext(model_levels)

CASE(tp_HD209458b_Heng_mid,        &
     tp_HD209458b_Heng_smooth_mid, &
     tp_HD209458b_iro,             &
     tp_Y_Dwarf,                   &
     tp_file,                      &
     tp_HD209458b_Heng_smooth_day, &
     tp_GJ1214b)

  profile_available = .TRUE.   ! T-p profile can be sent to stdout

  ! The initial T-P profile for a Hot Jupiter
  ! based on the midway point between the day and night-side
  ! profiles, or on the IRO profile
  ! Also the intial T-P for a Y-Dwarf and a generic input profile
  profile_type = 2     ! T specified

  ! Set the number used to identify the T-P profile 
  ! in the forcing
  SELECT CASE(tprofile_number)

  CASE(tp_HD209458b_Heng_mid)
    tforce_local = tf_HD209458b_Heng

  CASE(tp_HD209458b_Heng_smooth_mid, tp_HD209458b_Heng_smooth_day)
    tforce_local = tf_HD209458b_Heng_smooth

  CASE(tp_HD209458b_iro)
    tforce_local = tf_HD209458b_iro

  CASE(tp_GJ1214b)
    tforce_local = tf_GJ1214b

  CASE(tp_Y_Dwarf)
    tforce_local = tf_Y_Dwarf

  CASE(tp_file)
    tforce_local = tf_file
    CALL read_file_temp(profile_filename)

  END SELECT

  ! Set surface values.
  ! Set temperature as the mid-points:
  t0_profile(0) = planet_profile_temp(tprofile_number, exner_profile_ext(0), &
                                      tforce_local)
  exner_theta(0) = exner_profile_ext(0)
  t_theta(0) = t0_profile(0)

  ! Build the potential temperature on rho levels
  DO k = 1, model_levels
    ! Use hydrostatic balance to create the pressure at the next level
    exner_profile_ext(k) = exner_profile_ext(k-1) *                          &
        EXP(-g_theta_profile(k-1) * (levs(k) - levs(k-1)) /                  &
           (cp * t0_profile(k-1)) )

    ! Calculate pot. temperature on rho level
    t0_profile(k) = planet_profile_temp(tprofile_number,                     &
                        exner_profile_ext(k), tforce_local)

    ! Calculate density
    rho_profile(k) = p_zero*(exner_profile_ext(k)**(recip_kappa)) /          &
                    (r*t0_profile(k))
  END DO

  ! Exponential extrapolation to model top
  extrap1 = (z_theta(model_levels)-z_rho(model_levels-2)) /                  &
       (z_rho(model_levels-1)-z_rho(model_levels-2))
  extrap2 = (z_theta(model_levels)-z_rho(model_levels-1)) /                  &
       (z_rho(model_levels-1)-z_rho(model_levels-2))
  exner_top = exner_profile_ext(model_levels)**extrap1 *                     &
       exner_profile_ext(model_levels-1)**(-extrap2)

  ! Set exner_profile_ext(model_levels+1) so that extrapolated exner_top is
  ! obtained using 0.5*(exner_profile_ext(model_levels+1) +
  !                     exner_profile_ext(model_levels))
  exner_profile_ext(model_levels+1) = 2.0*exner_top -                        &
                                       exner_profile_ext(model_levels)

  ! Build the potential temperature on theta levels,
  ! done separately to avoid cumulative errors.
  theta_profile(0) = t_theta(0) / exner_theta(0)

  DO k = 1, model_levels
    ! Use hydrostatic balance to create the pressure at the next level
    exner_theta(k) = exner_theta(k-1) *                                      &
        EXP(-g_theta_profile(k-1) * (z_theta(k) - z_theta(k-1)) /            &
           (cp * t_theta(k-1)) )

    ! Calculate temperature on theta level
    t_theta(k) = planet_profile_temp(tprofile_number,                        &
                     exner_theta(k), tforce_local)
    
    ! Derive potential temperature on theta level
    theta_profile(k) = t_theta(k) / exner_theta(k)
  END DO


CASE(tp_BRuntV)

  ! Jablonowski & Williamson 2006 baroclinic wave test
  profile_type = 2   ! T specified
  DO k = 0, model_levels
    eta(k) = baro_eta_conv(xi2_p, depth_factor*levs(k))
    t0_profile(k) = baro_T_p(xi2_p, eta(k))

    ! Make a first guess at exner, theta...
    exner_profile_ext(k) = (p_zero*eta(k))**kappa
    theta_profile(k) = t0_profile(k) / exner_profile_ext(k)
  END DO
  exner_profile_ext(model_levels+1)=SQRT(exner_profile_ext(model_levels)**2 / &
                                         exner_profile_ext(model_levels-1))
  DO k = 1, model_levels
    ! ...and rho
    rho_profile(k) = p_zero*eta(k) / (r*t0_profile(k))
  END DO

  ! Return actual surface temperature and pressure
  p_surface = p_zero*eta(0)  ! Boundary condition for Newton solver
  t_surface = t0_profile(0)

CASE(tp_BV_isoth)

  ! Staniforth Deep atmosphere baroclinic wave test
  IF (.NOT. ALLOCATED(tau1)) THEN
    icode = 10
    CALL ereport(routinename, icode,                                     &
      'Attempted to use tprofile for deep baroclinic wave test ' //      &
      'without setting initial_profile = 3 (deep_baro_inst)')
  END IF

  profile_type = 2   ! T specified
  DO k = 0, model_levels
    IF (l_shallow) THEN
      r_height = planet_radius
    ELSE
      IF (k == 0) THEN
        r_height = surface_offset
      ELSE
        r_height = xi3_at_theta(k)
        r_height = xi3_at_rho(k)
      END IF
    END IF

    tau2_term = (r_height/planet_radius*COS(xi2_p)) ** k_const         &
              - (k_const/(k_const + 2.0)) *                            &
                (r_height/planet_radius*COS(xi2_p)) **(k_const + 2.0)

    t0_profile(k) = (planet_radius/r_height)**2 / (tau1(k) - tau2(k)*tau2_term)

    exner_profile_ext(k) = (EXP(-g/r*(tau1_int(k)-tau2_int(k)*tau2_term))) &
                              **kappa

    theta_profile(k) = t0_profile(k) / exner_profile_ext(k)

    IF (k >= 1) THEN
      rho_profile(k) = p_zero*exner_profile_ext(k)**(recip_kappa) / &
                        (r*t0_profile(k))

      IF (.NOT. l_shallow) THEN
        ! U coriolis term
        u_temp = (g/planet_radius)*k_const*tau2_int(k)                  &
                  *((r_height/planet_radius*COS(xi2_p))**(k_const-1)    &
                  - (r_height/planet_radius*COS(xi2_p))**(k_const+1))   &
                   * t0_profile(k)

        u_temp  = - omega*r_height*COS(xi2_p)                                &
                   + SQRT( (omega*r_height*COS(xi2_p))**2                    &
                          + r_height*COS(xi2_p)*u_temp )

        u_term(k) = u_temp**2/r_height+2.0*omega*u_temp*COS(xi2_p)
      END IF
    END IF
  END DO

  ! Better guess for topmost level
  exner_profile_ext(model_levels+1) = exner_profile_ext(model_levels)        &
           + 2.0 * (xi3_at_theta(model_levels) - xi3_at_rho(model_levels))   &
           / ((r/kappa)*theta_profile(model_levels))                         &
           * (u_term(model_levels)-g_theta_profile(model_levels))

  p_surface = p_zero * exner_profile_ext(0)**recip_kappa
  t_surface = t0_profile(0)

CASE (tp_dump)
  ! Do nothing; fields are taken from the input dump.

CASE DEFAULT

  icode = 20
  WRITE(umMessage,'(A,I0)') &
   'Selected idealised setup does not support tprofile_number:', tprofile_number
  CALL ereport(routinename, icode, umMessage)

END SELECT

! Print profile if requested:
IF (print_profile) THEN
  IF (profile_available) THEN
    WRITE(umMessage,'(A)') 'Initial T-p profile (before Newton-Raphson):'
    CALL umprint(umMessage, src=routinename, level=PrStatus_Normal, pe=0)
    DO k=1, model_levels
      WRITE(umMessage,'(I4,2E16.8)') k, theta_profile(k)*exner_theta(k),     &
                             (exner_theta(k)**recip_kappa) * p_zero
      CALL umprint(umMessage, src=routinename, level=PrStatus_Normal, pe=0)
    END DO
  END IF
END IF

CALL newton_solver(exner_profile_ext, theta_profile, rho_profile, t0_profile,&
       geopotential_factor, levs, intw_rho2w, intw_w2rho, p_surface,         &
       model_levels, profile_type, u_term, initial_tolerance)

! Now re-interpolate for the Exner pressure at theta levels                   
! In all but top Level
DO k=1, model_levels-1
  exner_theta(k) = (intw_rho2w(k,2)*exner_profile_ext(k)+                    &
                    intw_rho2w(k,1)*exner_profile_ext(k+1))
END DO

! The Newton solver doesn't calculate exner_ref_pro(model_levels+1).
! Use exponential extrapolation of exner to model top; then set
! exner_ref_pro(model_levels+1) such that this extrapolated value is
! returned when vertical averaging is done.

! Exponential extrapolation to model top:
extrap1 = (z_theta(model_levels)-z_rho(model_levels-2)) /                    &
          (z_rho(model_levels-1)-z_rho(model_levels-2))
extrap2 = (z_theta(model_levels)-z_rho(model_levels-1)) /                    &
          (z_rho(model_levels-1)-z_rho(model_levels-2))
exner_top = exner_profile_ext(model_levels)**extrap1 *                       &
            exner_profile_ext(model_levels-1)**(-extrap2)

! Set exner_profile_ext(model_levels+1) so that extrapolated exner_top is
! obtained using 0.5*(exner_profile_ext(model_levels+1) +
!                     exner_profile_ext(model_levels))
exner_profile_ext(model_levels+1) = 2.0*exner_top -                          &
                                     exner_profile_ext(model_levels)

! Print final profile if requested:
IF (print_profile) THEN
  WRITE(umMessage,'(A)') &
    'Initial T-p profile (after Newton-Raphson) and g_theta profile:'
  CALL umprint(umMessage, src=routinename, level=PrStatus_Normal, pe=0)
  DO k=1, model_levels
    WRITE(umMessage,'(I4,2E16.8,F13.7)') k, theta_profile(k)*exner_theta(k), &
          (exner_theta(k)**recip_kappa) * p_zero, g_theta_profile(k)
    CALL umprint(umMessage, src=routinename, level=PrStatus_Normal, pe=0)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

FUNCTION planet_profile_temp(tprofile_number, exner_pressure, tforce_local)  &
         RESULT(temperature)

! Description:
!      Wrapper for call to calculation of initial temperature profile

!   Variables to setup the Hot Jupiter Starting T-P profile
USE hd209458b_forcing_mod, ONLY:                                             &
    hd209458b_night_temp,                                                    &
    hd209458b_day_temp,                                                      &
    iro_HD209458b_temp

!  Variables to setup GJ1214 (Warm Neptune) Starting T-P profile
USE gj1214b_forcing_mod, ONLY:                                               &
    gj1214b_mean_temp

!   Variables to setup Y Dwarf starting T-P
USE y_dwarf_forcing_mod, ONLY:                                               &
    y_dwarf_temperature

!   Variables to setup starting profile from file
USE rcf_ideal_file_forcing_mod, ONLY:                                        &
    file_temp

IMPLICIT NONE

INTEGER, INTENT(IN)  :: tprofile_number, tforce_local
REAL,    INTENT(IN)  :: exner_pressure

REAL                 :: temperature

CHARACTER (LEN=*), PARAMETER  :: routinename='PLANET_PROFILE_TEMP'
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(tprofile_number)

CASE(tp_HD209458b_Heng_mid, tp_HD209458b_Heng_smooth_mid)
  temperature = 0.5*(                                                        &
      hd209458b_night_temp(exner_pressure, tforce_local) +                   &
      hd209458b_day_temp(exner_pressure, tforce_local))

CASE(tp_HD209458b_iro)
   temperature = iro_HD209458b_temp(exner_pressure, tforce_local)

CASE(tp_GJ1214b)
   temperature = gj1214b_mean_temp(exner_pressure)

CASE(tp_Y_Dwarf)
   temperature = y_dwarf_temperature(exner_pressure)

CASE(tp_file)
   temperature = file_temp(exner_pressure)

CASE(tp_HD209458b_Heng_smooth_day)
  temperature = hd209458b_day_temp(exner_pressure, tforce_local)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION planet_profile_temp

END SUBROUTINE rcf_ideal_1d_profiles

END MODULE rcf_ideal_1d_profiles_mod
