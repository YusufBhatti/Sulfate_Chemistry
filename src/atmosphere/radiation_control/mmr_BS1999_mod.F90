! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to calculate gas abundances for a hot Jupiter.

! Description:
!   Calculates gas abundance as a function of pressure and temperature
!   from relations in Burrows & Sharp, ApJ, 1999 (for hot Jupiters).

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE mmr_BS1999_mod

USE realtype_rd, ONLY: RealK
USE rad_ccf, ONLY: r_gas, A_H, A_He
USE gas_list_pcf, ONLY: molar_weight

IMPLICIT NONE

! Ideal gas constant in cal/(mol*K)
REAL(RealK), PARAMETER :: r_gas_cal=1.9858775_RealK
! Number of atmospheres per Pa
REAL(RealK), PARAMETER :: atmprpa=1.0_RealK/1.01325e+05_RealK
! Average number of oxygen atoms removed per silicon atom
REAL(RealK), PARAMETER :: x_Si=3.28_RealK
! Temperature below which oxygen is removed by silicon atoms
REAL(RealK), PARAMETER :: T_remove_O=1500_RealK
! Characteristic transition scale in T over which O is removed
REAL(RealK), PARAMETER :: T_trans_scale_remove_O=10.0_RealK
! Characteristic transition scale in T over which the abundance changes
REAL(RealK) ::                                                                 &
    T_trans_scale_TiO = 10.0_RealK,                                            &
    T_trans_scale_VO  = 10.0_RealK,                                            &
    T_trans_scale_Na  = 20.0_RealK,                                            &
    T_trans_scale_K   = 20.0_RealK,                                            &
    T_trans_scale_Li  = 20.0_RealK,                                            &
    T_trans_scale_Rb  = 20.0_RealK,                                            &
    T_trans_scale_Cs  = 20.0_RealK

! Elemental abundances as log(epsilon_z) = log(N_z/N_H) + 12, where N_z is the
! number of atoms of species z
REAL(RealK), PARAMETER ::                                                      &
    log_eps_C  = 8.50_RealK,                                                   &
    log_eps_N  = 7.86_RealK,                                                   &
    log_eps_O  = 8.76_RealK,                                                   &
    log_eps_Na = 6.24_RealK,                                                   &
    log_eps_Si = 7.51_RealK,                                                   &
    log_eps_K  = 5.03_RealK,                                                   &
    log_eps_Ti = 4.95_RealK,                                                   &
    log_eps_V  = 3.93_RealK,                                                   &
    log_eps_Li = 3.26_RealK,                                                   &
    log_eps_Rb = 2.52_RealK,                                                   &
    log_eps_Cs = 1.08_RealK

! Convert elemental abundances to number of atoms relative to H:
! APrim_z = N_z/N_H
REAL(RealK), PARAMETER ::                                                      &
    APrim_C  = (10.0_RealK)**(log_eps_C  - 12.0_RealK),                        &
    APrim_N  = (10.0_RealK)**(log_eps_N  - 12.0_RealK),                        &
    APrim_O  = (10.0_RealK)**(log_eps_O  - 12.0_RealK),                        &
    APrim_Na = (10.0_RealK)**(log_eps_Na - 12.0_RealK),                        &
    APrim_Si = (10.0_RealK)**(log_eps_Si - 12.0_RealK),                        &
    APrim_K  = (10.0_RealK)**(log_eps_K  - 12.0_RealK),                        &
    APrim_Ti = (10.0_RealK)**(log_eps_Ti - 12.0_RealK),                        &
    APrim_V  = (10.0_RealK)**(log_eps_V  - 12.0_RealK),                        &
    APrim_Li = (10.0_RealK)**(log_eps_Li - 12.0_RealK),                        &
    APrim_Rb = (10.0_RealK)**(log_eps_Rb - 12.0_RealK),                        &
    APrim_Cs = (10.0_RealK)**(log_eps_Cs - 12.0_RealK)

! Equilibrium coefficients from Burrows & Sharp, ApJ, 1999:
REAL(RealK), PARAMETER ::                                                      &
    eqcoeff_1(5) = (/                                                          &
         1.106131e+06_RealK,                                                   &
        -5.6895e+04_RealK,                                                     &
         62.565_RealK,                                                         &
        -5.81396e-04_RealK,                                                    &
         2.346515e-08_RealK /),                                                &
    eqcoeff_2(5) = (/                                                          &
         8.16413e+05_RealK,                                                    &
        -2.9109e+04_RealK,                                                     &
         58.5878_RealK,                                                        &
        -7.8284e-04_RealK,                                                     &
         4.729048e-08_RealK /)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MMR_BS1999_MOD'
CONTAINS

! ****************************************************************************

! Calculates the mass mixing ratio
SUBROUTINE calc_gas_mixing_ratio_BS1999(p, t, iump, n_profile, n_layer, mmr)

USE gas_list_pcf, ONLY: ip_h2o, ip_co, ip_ch4, ip_nh3,                         &
                        ip_tio, ip_vo, ip_h2, ip_he, ip_na, ip_k,              &
                        ip_li, ip_rb, ip_cs
USE ereport_mod, ONLY: ereport
USE rad_pcf, ONLY: i_err_fatal
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Pressure
REAL(RealK), INTENT(IN) :: p(:,:)
! Temperature
REAL(RealK), INTENT(IN) :: t(:,:)
! Identifier of molecule
INTEGER, INTENT(IN) :: iump
! Number of profiles
INTEGER, INTENT(IN) :: n_profile
! Number of atmospheric layers
INTEGER, INTENT(IN) :: n_layer
! Mass mixing ratio
REAL(RealK), INTENT(OUT) :: mmr(:,:)

! Array for storage of equilibrium constants
REAL(RealK) :: k(n_profile, n_layer)
! Fraction of O available to deplete that is depleted
REAL(RealK) :: frac_remove_O(n_profile, n_layer)
! Total mass density
REAL(RealK) :: rho_tot(n_profile, n_layer)
! H2 and He partial pressures densities
REAL(RealK) :: pp_H2(n_profile, n_layer)
REAL(RealK) :: pp_He(n_profile, n_layer)
! Log (base 10) of pressure in bars
REAL(RealK) :: log_p_bar(n_profile, n_layer)
! Condensation temperature
REAL(RealK) :: t_cond(n_profile, n_layer)
! Number of oxygen atoms once depleted
REAL(RealK) :: APrim_O_depleted(n_profile, n_layer)
! Working variable in the calculation of CO equilibrium abundance
REAL(RealK) :: CO_work(n_profile, n_layer)
! Working variable in the calculation of NH3 equilibrium abundance
REAL(RealK) :: NH3_work(n_profile, n_layer)

! Loop integers
INTEGER :: i, j

INTEGER :: ierr
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'calc_gas_mixing_ratio_BS1999'


! Calculate fraction of available O to deplete that have depleted given
! temperature in layer
DO j = 1, n_layer
  DO i = 1, n_profile
    frac_remove_O(i,j) = 1.0_RealK /                                           &
      (EXP((t(i,j) - T_remove_O)/T_trans_scale_remove_O) + 1.0_RealK)
  END DO
END DO

! Calculate H2 and He number densities
CALL calc_H2_He_pp(p, n_profile, n_layer, pp_H2, pp_He)

! Total atmospheric mass density
! NB: All mass is assumed to be held in H2 and He. Should possibly use
! a mean molecular mass for all molecules instead to find total mass
! of atmosphere. pp_H2 converted from Pa to Ba
DO j = 1, n_layer
  DO i = 1, n_profile
    rho_tot(i,j) =                                                             &
      molar_weight(ip_h2)*1.0e-03_RealK*pp_H2(i,j)/(r_gas*t(i,j)) +            &
      molar_weight(ip_he)*1.0e-03_RealK*pp_He(i,j)/(r_gas*t(i,j))
  END DO
END DO

! Calculate partial pressures of gases relative to pp_H2. The unit of
! pp_H2 in Burrows & Sharp, ApJ, 1999 are in atmospheres, must convert
! pp_H2 to atmospheres before use in these formulas.
SELECT CASE(iump)

CASE(ip_h2)
  ! H2 partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 1.0_RealK ! pp_H2/pp_H2
    END DO
  END DO

CASE(ip_he)
  ! He partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = pp_He(i,j)/pp_H2(i,j)
    END DO
  END DO

CASE(ip_co, ip_ch4, ip_h2o)
  ! CO partial pressure relative to H2 partial pressure
  CALL calc_K1(t, n_profile, n_layer, k)
  DO j = 1, n_layer
    DO i = 1, n_profile
      APrim_O_depleted(i,j) = APrim_O - x_Si*APrim_Si*frac_remove_O(i,j)

      CO_work(i,j) = APrim_C + APrim_O_depleted(i,j) +                         &
        (pp_H2(i,j)*atmprpa)**2 / (2.0_RealK*k(i,j))

      mmr(i,j) = CO_work(i,j) -                                                &
        SQRT( CO_work(i,j)**2 - 4.0_RealK*APrim_C*APrim_O_depleted(i,j) )
    END DO
  END DO

  ! CH4 partial pressure relative to H2 partial pressure
  IF (iump == ip_ch4) THEN
    DO j = 1, n_layer
      DO i = 1, n_profile
        mmr(i,j) = 2.0_RealK*APrim_C - mmr(i,j)
      END DO
    END DO
  END IF

  ! H2O partial pressure relative to H2 partial pressure
  IF (iump == ip_h2o) THEN
    DO j = 1, n_layer
      DO i = 1, n_profile
        mmr(i,j) = 2.0_RealK*APrim_O_depleted(i,j) - mmr(i,j)
      END DO
    END DO
  END IF

CASE(ip_nh3)
  ! NH3 partial pressure relative to H2 partial pressure
  CALL calc_K2(t, n_profile, n_layer, k)
  DO j = 1, n_layer
    DO i = 1, n_profile
      NH3_work(i,j) = (pp_H2(i,j)*atmprpa)**2/(8.0_RealK*k(i,j))

      IF ( APrim_N / NH3_work(i,j) < 1.0e-10_RealK ) THEN
        mmr(i,j) = 2.0_RealK*APrim_N
      ELSE
        mmr(i,j) = 2.0_RealK*(                                                 &
          SQRT( NH3_work(i,j)*(2.0_RealK*APrim_N + NH3_work(i,j)) )            &
          - NH3_work(i,j) )
      END IF
    END DO
  END DO

CASE(ip_tio)
  ! TiO partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_tio, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_Ti/                                           &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_TiO) + 1.0_RealK)
    END DO
  END DO

CASE(ip_vo)
  ! VO partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_vo, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_V/                                            &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_VO) + 1.0_RealK)
    END DO
  END DO

CASE(ip_na)
  ! Sodium partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_na, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_Na/                                           &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_Na) + 1.0_RealK)
    END DO
  END DO

CASE(ip_k)
  ! Potassium partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_k, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_K/                                            &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_Na) + 1.0_RealK)
    END DO
  END DO

CASE(ip_li)
  ! Lithium partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_li, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_Li/                                           &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_Li) + 1.0_RealK)
    END DO
  END DO

CASE(ip_rb)
  ! Rubidium partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_rb, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_Rb/                                           &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_Rb) + 1.0_RealK)
    END DO
  END DO

CASE(ip_cs)
  ! Cesium partial pressure relative to H2 partial pressure
  DO j = 1, n_layer
    DO i = 1, n_profile
      log_p_bar(i,j) = LOG10(p(i,j)/1.0e+05_RealK)
    END DO
  END DO
  CALL cond_temp_curve(log_p_bar, ip_cs, n_profile, n_layer, t_cond)
  DO j = 1, n_layer
    DO i = 1, n_profile
      mmr(i,j) = 2.0_RealK*APrim_Cs/                                           &
        (EXP(-(t(i,j) - t_cond(i,j))/T_trans_scale_Cs) + 1.0_RealK)
    END DO
  END DO

CASE DEFAULT
  WRITE(cmessage,'(A,I3,A)') 'Gas ',iump,' is not implemented.'
  ierr=i_err_fatal
  CALL ereport(RoutineName, ierr, cmessage)

END SELECT

DO j = 1, n_layer
  DO i = 1, n_profile
    ! Convert to actual pressure (not relative to H2 partial pressure), and
    ! convert to mass density
    mmr(i,j) = mmr(i,j)                                                        &
      *pp_H2(i,j)*molar_weight(iump)*1.0e-03_RealK/(r_gas*t(i,j))

    ! Molecular mass fraction
    mmr(i,j) = mmr(i,j)/rho_tot(i,j)
  END DO
END DO

END SUBROUTINE calc_gas_mixing_ratio_BS1999

! ****************************************************************************

! Calculates the equilibrium coefficients K
SUBROUTINE calc_K1(t, n_profile, n_layer, k)

IMPLICIT NONE

! Number of profiles
INTEGER, INTENT(IN) :: n_profile
! Number of atmospheric layers
INTEGER, INTENT(IN) :: n_layer
! Temperature
REAL(RealK), INTENT(IN) :: t(:,:)
! Equilibrium constants
REAL(RealK), INTENT(OUT) :: k(n_profile, n_layer)

! Loop integers
INTEGER :: i, j

DO j = 1, n_layer
  DO i = 1, n_profile
    k(i,j) = EXP( (eqcoeff_1(1)/t(i,j) + eqcoeff_1(2)                          &
      + eqcoeff_1(3)*t(i,j)                                                    &
      + eqcoeff_1(4)*t(i,j)**2                                                 &
      + eqcoeff_1(5)*t(i,j)**3)                                                &
      / (r_gas_cal*t(i,j)) )
  END DO
END DO

END SUBROUTINE calc_K1

SUBROUTINE calc_K2(t, n_profile, n_layer, k)

IMPLICIT NONE

! Number of profiles
INTEGER, INTENT(IN) :: n_profile
! Number of atmospheric layers
INTEGER, INTENT(IN) :: n_layer
! Temperature
REAL(RealK), INTENT(IN) :: t(:,:)
! Equilibrium coefficients
REAL(RealK), INTENT(OUT) :: k(n_profile, n_layer)

! Loop integers
INTEGER :: i, j

DO j = 1, n_layer
  DO i = 1, n_profile
    k(i,j) = EXP( (eqcoeff_2(1)/t(i,j) + eqcoeff_2(2)                          &
      + eqcoeff_2(3)*t(i,j)                                                    &
      + eqcoeff_2(4)*t(i,j)**2                                                 &
      + eqcoeff_2(5)*t(i,j)**3)                                                &
      / (r_gas_cal*t(i,j)) )
  END DO
END DO

END SUBROUTINE calc_K2

! ****************************************************************************

! Calculates H2 and He partial pressures
SUBROUTINE calc_H2_He_pp(p, n_profile, n_layer, pp_H2, pp_He)

IMPLICIT NONE

! Number of profiles
INTEGER, INTENT(IN) :: n_profile
! Number of atmospheric layers
INTEGER, INTENT(IN) :: n_layer
! Total gas pressure
REAL(RealK), INTENT(IN) :: p(:,:)
! Partial pressures of H2 and He
REAL(RealK), INTENT(OUT) :: pp_H2(n_profile, n_layer)
REAL(RealK), INTENT(OUT) :: pp_He(n_profile, n_layer)

! Loop integers
INTEGER :: i, j

! Calculate partial pressures (weighing total pressure by abundances,
! assume all H is in H2)
DO j = 1, n_layer
  DO i = 1, n_profile
    pp_H2(i,j) = p(i,j)*(A_H/2.0_RealK)/(A_H/2.0_RealK + A_He)
    pp_He(i,j) = p(i,j)*(A_He)/(A_H/2.0_RealK + A_He)
  END DO
END DO

END SUBROUTINE calc_H2_He_pp

SUBROUTINE cond_temp_curve(log_p_bar, iump, n_profile, n_layer, t_cond)

USE gas_list_pcf, ONLY: ip_tio, ip_vo, ip_na, ip_k, ip_li, ip_rb, ip_cs

IMPLICIT NONE

! Number of profiles
INTEGER, INTENT(IN) :: n_profile
! Number of atmospheric layers
INTEGER, INTENT(IN) :: n_layer
! Identifier of molecule
INTEGER, INTENT(IN) :: iump
! Total gas pressure in bar
REAL(RealK), INTENT(IN) :: log_p_bar(n_profile, n_layer)
! Condensation temperature
REAL(RealK), INTENT(OUT) :: t_cond(n_profile, n_layer)

! Polynomial fits to condensation curves
REAL(RealK), PARAMETER ::                                                      &
  poly_highp_tio_vo(2) =  (/                                                   &
    -3.96274342e-05_RealK,                                                     &
     5.20741797e-04_RealK /),                                                  &
!   TiO and VO
  poly_highp_na(2) = (/                                                        &
    -5.58256919e-05_RealK,                                                     &
     8.81851644e-04_RealK /),                                                  &
!   Na
  poly_highp_k(2) = (/                                                         &
    -5.46977180e-05_RealK,                                                     &
     8.19104478e-04_RealK /),                                                  &
!   K
  poly_highp_li(2) = (/                                                        &
    -3.50995394e-05_RealK,                                                     &
     6.51993843e-04_RealK /),                                                  &
!   Li
  poly_highp_rb(2) = (/                                                        &
    -6.06654087e-05_RealK,                                                     &
     8.09569974e-04_RealK /),                                                  &
!   Rb
  poly_highp_cs(2) = (/                                                        &
    -5.29210264e-05_RealK,                                                     &
     7.71577097e-04_RealK /)
!   Cs
REAL(RealK), PARAMETER ::                                                      &
  poly_lowp_tio_vo(2)  = (/                                                    &
    -2.91788901e-05_RealK,                                                     &
     5.11801742e-04_RealK /),                                                  &
!   TiO and VO
  poly_lowp_na(2) = (/                                                         &
    -6.69921629e-05_RealK,                                                     &
     8.90116829e-04_RealK /),                                                  &
!   Na
  poly_lowp_k(2) = (/                                                          &
    -6.46633572e-05_RealK,                                                     &
     8.29549449e-04_RealK /),                                                  &
!   K
  poly_lowp_li(2) = (/                                                         &
    -3.55469185e-05_RealK,                                                     &
     6.52116945e-04_RealK /),                                                  &
!   Li
  poly_lowp_rb(2) = (/                                                         &
    -3.19328287e-05_RealK,                                                     &
    8.69542964e-04_RealK /),                                                   &
!   Rb
  poly_lowp_cs(2) = (/                                                         &
    -3.85306167e-05_RealK,                                                     &
     7.63040762e-04_RealK /)
!   Cs

! Log (base 10) of transition point between polynomial fits
REAL(RealK), PARAMETER ::                                                      &
  log_p_trans_tio_vo =  1.0_RealK,                                             &
!   TiO and VO
  log_p_trans_na =      1.0_RealK,                                             &
!   Na
  log_p_trans_k =       1.0_RealK,                                             &
!   K
  log_p_trans_li =      1.0_RealK,                                             &
!   Li
  log_p_trans_rb =     -2.0_RealK,                                             &
!   Rb
  log_p_trans_cs =      1.0_RealK
!   Cs

! Polynomial fit data for current gas
REAL(RealK) ::                                                                 &
  poly_highp(2),                                                               &
  poly_lowp(2),                                                                &
  log_p_trans

! Working variable
REAL(RealK) :: work(n_profile, n_layer)

! Loop integers
INTEGER :: i, j

! Set polynomial fit data for current absorber
SELECT CASE(iump)
CASE(ip_tio, ip_vo)
  poly_highp = poly_highp_tio_vo
  poly_lowp = poly_lowp_tio_vo
  log_p_trans = log_p_trans_tio_vo
CASE(ip_na)
  poly_highp = poly_highp_na
  poly_lowp = poly_lowp_na
  log_p_trans = log_p_trans_na
CASE(ip_k)
  poly_highp = poly_highp_k
  poly_lowp = poly_lowp_k
  log_p_trans = log_p_trans_k
CASE(ip_li)
  poly_highp = poly_highp_li
  poly_lowp = poly_lowp_li
  log_p_trans = log_p_trans_li
CASE(ip_rb)
  poly_highp = poly_highp_rb
  poly_lowp = poly_lowp_rb
  log_p_trans = log_p_trans_rb
CASE(ip_cs)
  poly_highp = poly_highp_cs
  poly_lowp = poly_lowp_cs
  log_p_trans = log_p_trans_cs
END SELECT

DO j = 1, n_layer
  DO i = 1, n_profile
    work(i,j) = 1.0_RealK /                                                    &
      (EXP(-(log_p_bar(i,j)-log_p_trans)/0.1_RealK) + 1.0_RealK)

    t_cond(i,j) =                                                              &
      1.0_RealK / ( (poly_lowp(1)*log_p_bar(i,j) + poly_lowp(2))) *            &
      (1.0_RealK - work(i,j)) +                                                &
      1.0_RealK / ( (poly_highp(1)*log_p_bar(i,j) + poly_highp(2))) *          &
      work(i,j)
  END DO
END DO

END SUBROUTINE cond_temp_curve

END MODULE mmr_BS1999_mod
