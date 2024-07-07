! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module containing the subroutine CALC_SURF_AREA, which returns
!   surface areas (considering hygroscopic growth) and wet radii
!   for some CLASSIC aerosols.
!
! Method:
!   Identify if the following CLASSIC aerosol types are being modelled:
!   ammonium sulphate, soot (BC), OCFF, biogenic secondary organic aerosol
!   and sea-salt aerosol. If so do the following calculations separately
!   for their respective aerosol modes:
!
!   * Get the total number concentration (m-3) for the aerosol types which
!     are transported as mass mixing ratios (mmr). This is done as follows
!                            rho_air           1
!      N (m-3) = m.m.r. * ------------  * -------------, with
!                         rho_particle     V_particle
!
!      V_particle = (4*PI/3) * (r_bar**3) * exp (4.5*(ln (sigma))**2),
!      2*r_bar    = geometric diameter
!      sigma      = geometric standard deviation
!
!      Rather than doing the whole calculation here, we use
!                          3 * rho_air
!      N (m-3) = m.m.r. * -------------, where
!                          denom_X
!      denom_X = 3 * rho_particle * V_particle 
!
!     All variables whose name starts with 'denom_' are then specific to  
!     given aerosol types and modes. They were calculated offline by the
!     standalone program calc_pm_params.f90 and the values used here are
!     taken from the module calc_pm_diags_mod.
!
!     Note that this calculation is needed for all aerosol types except for
!     sea-salt aerosol, which is already expressed as total number
!     concentration (m-3).
!
!   * Call the routine GROW_PARTICLES to calculate the parameters
!     of the log-normal distribution after hygroscopic growth (i.e.
!     "wet" radii and "wet" standard deviation).
!
!   * Use the parameters derived as indicated above to finally
!     calculate surface areas for each aerosol type and mode.
!     Return surface areas as well as wet radii.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: aerosols
!
! Code Description:
!   Language:  FORTRAN 2003.
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE calc_surf_area_mod

IMPLICIT NONE

! All variables and functions private by default
PRIVATE
PUBLIC :: calc_surf_area

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_SURF_AREA_MOD'

CONTAINS

SUBROUTINE calc_surf_area (n_pnts, rho_air, rh_frac,            &
                           so4_aitken,    so4_accum,            &
                           soot_fresh,    soot_aged,            &
                           ocff_fresh,    ocff_aged, biogenic,  &
                           sea_salt_film, sea_salt_jet,         &
                           sa_so4_ait,    sa_so4_acc,           &
                           sa_bc_fresh,   sa_bc_aged,           &
                           sa_ocff_fresh, sa_ocff_aged, sa_soa, &
                           sa_ss_film,    sa_ss_jet,            &
                           wr_so4_ait,    wr_so4_acc,           &
                           wr_bc_fresh,   wr_bc_aged,           &
                           wr_ocff_fresh, wr_ocff_aged, wr_soa, &
                           wr_ss_film,    wr_ss_jet)

USE run_aerosol_mod,    ONLY: l_sulpc_so2, l_soot, l_ocff,                  &
                              l_use_seasalt_sulpc
USE rad_input_mod,      ONLY: l_use_seasalt_direct, l_use_seasalt_indirect, &
                              l_use_biogenic
USE mphys_inputs_mod,   ONLY: l_use_seasalt_autoconv
USE aero_params_mod,    ONLY:                                            &
  r_bar_su_ait,  sigma_su_ait,  r_bar_su_acc,  sigma_su_acc,             &
  r_bar_bc_fr,   sigma_bc_fr,   r_bar_bc_ag,   sigma_bc_ag,              &
  r_bar_ocff_fr, sigma_ocff_fr, r_bar_ocff_ag, sigma_ocff_ag,            &
  r_bar_ss_fi,   sigma_ss_fi,   r_bar_ss_je,   sigma_ss_je,              &
  r_bar_soa,     sigma_soa
USE calc_pm_diags_mod,  ONLY: denom_su_ait,   denom_su_acc,  s_conv_fac, &
                              denom_bc_fr,    denom_bc_ag,               &
                              denom_ocff_fr,  denom_ocff_ag,             &
                              denom_soa
USE rad_pcf,            ONLY: ip_accum_sulphate, ip_aitken_sulphate,     &
                              ip_fresh_soot,     ip_aged_soot,           &
                              ip_ocff_fresh,     ip_ocff_aged,           &
                              ip_seasalt_film,   ip_seasalt_jet,         &
                              ip_biogenic
USE grow_particles_mod, ONLY: grow_particles
USE conversions_mod,    ONLY: pi

USE parkind1,           ONLY: jpim,  jprb     ! DrHook
USE yomhook,            ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Subroutine arguments
!
INTEGER, INTENT(IN) :: n_pnts           ! number of grid points

REAL, INTENT(IN)    :: rho_air (n_pnts) ! kg m-3
REAL, INTENT(IN)    :: rh_frac (n_pnts) ! RH (fraction: 0.0 - 0.999)

! Mass mixing ratios (kg kg-1) or aerosol number (m-3) from CLASSIC aerosols.
REAL, INTENT(IN)    :: so4_aitken    (n_pnts) ! kg kg-1
REAL, INTENT(IN)    :: so4_accum     (n_pnts)
REAL, INTENT(IN)    :: soot_fresh    (n_pnts)
REAL, INTENT(IN)    :: soot_aged     (n_pnts)
REAL, INTENT(IN)    :: ocff_fresh    (n_pnts)
REAL, INTENT(IN)    :: ocff_aged     (n_pnts)
REAL, INTENT(IN)    :: biogenic      (n_pnts)
REAL, INTENT(IN)    :: sea_salt_film (n_pnts) ! m-3
REAL, INTENT(IN)    :: sea_salt_jet  (n_pnts)

! Aerosol surface areas (cm2 cm-3) 
REAL, INTENT(OUT)   :: sa_so4_ait    (n_pnts)
REAL, INTENT(OUT)   :: sa_so4_acc    (n_pnts)
REAL, INTENT(OUT)   :: sa_bc_fresh   (n_pnts)
REAL, INTENT(OUT)   :: sa_bc_aged    (n_pnts)
REAL, INTENT(OUT)   :: sa_ocff_fresh (n_pnts)
REAL, INTENT(OUT)   :: sa_ocff_aged  (n_pnts)
REAL, INTENT(OUT)   :: sa_soa        (n_pnts)
REAL, INTENT(OUT)   :: sa_ss_film    (n_pnts)
REAL, INTENT(OUT)   :: sa_ss_jet     (n_pnts)

! Aerosol wet radii (cm) 
REAL, INTENT(OUT)   :: wr_so4_ait    (n_pnts)
REAL, INTENT(OUT)   :: wr_so4_acc    (n_pnts)
REAL, INTENT(OUT)   :: wr_bc_fresh   (n_pnts)
REAL, INTENT(OUT)   :: wr_bc_aged    (n_pnts)
REAL, INTENT(OUT)   :: wr_ocff_fresh (n_pnts)
REAL, INTENT(OUT)   :: wr_ocff_aged  (n_pnts)
REAL, INTENT(OUT)   :: wr_soa        (n_pnts)
REAL, INTENT(OUT)   :: wr_ss_film    (n_pnts)
REAL, INTENT(OUT)   :: wr_ss_jet     (n_pnts)

! Local variables
INTEGER :: i            ! loop counter for grid points
INTEGER :: aerosol_type ! input argument to GROW_PARTICLES

REAL :: n_tot      ! total number concentration (m-3)
 
! Parameters of the log-normal size distribution:
!
! Geometric mean radius (also called median radius):
! * r0  = dry radius (not declared because its calculation is not needed)
! * wr0 = wet radius (considering growth with humidity)
REAL :: wr0

! Natural logarithms of dry (r0) and wet (wr0) median radii
REAL :: ln_r0    ! dry
REAL :: ln_wr0   ! wet
!
! Natural logarithms of twice r0 (dry radius) & wr0 (wet radius)
REAL :: ln_2r0   ! dry
REAL :: ln_2wr0  ! wet
!
! Natural logarithms of the geometric standard deviation
REAL :: ln_sigma    ! dry
REAL :: ln_wsigma   ! wet

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_SURF_AREA'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise surface areas (SA)
sa_so4_ait    (:) = 0.0
sa_so4_acc    (:) = 0.0
sa_bc_fresh   (:) = 0.0
sa_bc_aged    (:) = 0.0
sa_ocff_fresh (:) = 0.0
sa_ocff_aged  (:) = 0.0
sa_soa        (:) = 0.0
sa_ss_film    (:) = 0.0
sa_ss_jet     (:) = 0.0

! Initialise wet radii (wr)
wr_so4_ait    (:) = 0.0
wr_so4_acc    (:) = 0.0
wr_bc_fresh   (:) = 0.0
wr_bc_aged    (:) = 0.0
wr_ocff_fresh (:) = 0.0
wr_ocff_aged  (:) = 0.0
wr_soa        (:) = 0.0
wr_ss_film    (:) = 0.0
wr_ss_jet     (:) = 0.0

DO i = 1,n_pnts

  ! The IF blocks below are used to check for the presence of different
  ! CLASSIC aerosol types which can be used.
  ! 
  ! Calculations are commented in detail for the first aerosol type evaluated.
  ! Then they are repeated for the remaining aerosol types considered.

  IF (l_sulpc_so2) THEN
    !-----------------------------------------------------------------------
    ! Ammonium sulphate from CLASSIC - Aitken mode
    !
    ! Identifier of aerosol type in CLASSIC.
    aerosol_type = ip_aitken_sulphate

    ! Calculate the total number concentration (m-3) of the given aerosol type
    ! using the parameter denom_X as indicated in the Method section.
    ! In the case of ammonium sulphate the mmr of sulphur is advected. Hence,
    ! the factor s_conv_fac is needed to convert to mmr of ammonium sulphate.
    n_tot = 3.0 * ( so4_aitken (i) * s_conv_fac ) * rho_air (i) / denom_su_ait

    ! Natural logarithms of the dry parameters of the log-normal distribution
    ! (r0 is the geometric mean radius, also called median radius;
    !  sigma is the geometric standard deviation).
    ln_r0    = LOG (r_bar_su_ait)
    ln_sigma = LOG (sigma_su_ait)

    ! Get natural logarithm for wet parameters of the log-normal distribution
    ! (radius and standard deviation) after growth with humidity.
    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    ! 'Wet' parameters used in the calculation of the surface area    
    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    ! Calculation of the aerosol surface area. Note that these two
    ! formulations are often used and both of them are equivalent:
    !   S_t = N_t * pi * (2r)**2 * exp {2* [ln(sigma)]^2}
    !   S_t = N_t * pi * exp {2 * ln(2r) + 2[ln(sigma)]^2}
    ! We use the second formula, with the factor 1e-2 to convert the resulting
    ! surface area from m2 m-3 to cm2 cm-3.
    sa_so4_ait (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)

    ! Store wet radius too (the factor 1e2 converts from m to cm).
    wr_so4_ait (i) = 1.0e2 * wr0

    !-----------------------------------------------------------------------
    ! Ammonium sulphate from CLASSIC - accumulation mode
    aerosol_type = ip_accum_sulphate

    n_tot = 3.0 * ( so4_accum (i) * s_conv_fac ) * rho_air (i) / denom_su_acc

    ln_r0    = LOG (r_bar_su_acc)
    ln_sigma = LOG (sigma_su_acc)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_so4_acc (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)
    wr_so4_acc (i) = 1.0e2  * wr0
  END IF

  IF (l_soot) THEN
    !-----------------------------------------------------------------------
    ! Black carbon (or soot) from CLASSIC - fresh  mode
    aerosol_type = ip_fresh_soot

    n_tot = 3.0 * soot_fresh (i) * rho_air (i) / denom_bc_fr

    ! Soot is considered to be hydrophobic - no need to calculate growth
    ln_2r0   = LOG (2*r_bar_bc_fr)
    ln_sigma = LOG (sigma_bc_fr)

    sa_bc_fresh (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2r0 + 2*(ln_sigma)**2)

    ! Wet radius is equal to dry radius, but need to convert from m to cm.
    wr_bc_fresh (i) =  1.0e2 * r_bar_bc_fr

    !-----------------------------------------------------------------------
    ! Black carbon (or soot) from CLASSIC - aged  mode
    aerosol_type = ip_aged_soot

    n_tot = 3.0 * soot_aged (i) * rho_air (i) / denom_bc_ag

    ln_2r0   = LOG (2*r_bar_bc_ag)
    ln_sigma = LOG (sigma_bc_ag)

    sa_bc_aged (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2r0 + 2*(ln_sigma)**2)
    wr_bc_aged (i) = 1.0e2  * r_bar_bc_ag    ! wet radius = dry radius
  END IF

  IF (l_ocff) THEN
    !-----------------------------------------------------------------------
    ! Organic carbon from fossil fuel (OCFF) in CLASSIC - fresh  mode
    aerosol_type = ip_ocff_fresh

    n_tot = 3.0 * ocff_fresh (i) * rho_air (i) / denom_ocff_fr

    ln_r0    = LOG (r_bar_ocff_fr)
    ln_sigma = LOG (sigma_ocff_fr)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_ocff_fresh(i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)
    wr_ocff_fresh(i) = 1.0e2  * wr0

    !-----------------------------------------------------------------------
    ! Organic carbon from fossil fuel (OCFF) in CLASSIC - aged  mode
    aerosol_type = ip_ocff_aged

    n_tot = 3.0 * ocff_aged (i) * rho_air (i) / denom_ocff_ag

    ln_r0    = LOG (r_bar_ocff_ag)
    ln_sigma = LOG (sigma_ocff_ag)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_ocff_aged(i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0  + 2*(ln_wsigma)**2)
    wr_ocff_aged(i) = 1.0e2  * wr0
  END IF

  IF (l_use_biogenic) THEN
    !-----------------------------------------------------------------------
    ! Climatology of biogenic organic aerosols (based on oxidation 
    ! from terpenes) used in CLASSIC
    aerosol_type = ip_biogenic

    n_tot = 3.0 * biogenic(i) * rho_air (i) / denom_soa

    ln_r0    = LOG (r_bar_soa)
    ln_sigma = LOG (sigma_soa)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_soa (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)
    wr_soa (i) = 1.0e2  * wr0
  END IF

  IF (l_use_seasalt_autoconv .OR. l_use_seasalt_sulpc .OR.              &
      l_use_seasalt_indirect .OR. l_use_seasalt_direct) THEN
    !-----------------------------------------------------------------------
    ! Sea-salt aerosol from CLASSIC - film mode
    aerosol_type = ip_seasalt_film

    ! Sea-salt is the only aerosol type in CLASSIC which is transported as
    ! total number (m-3), so no need to do any conversion here:
    n_tot = sea_salt_film (i)

    ln_r0    = LOG (r_bar_ss_fi)
    ln_sigma = LOG (sigma_ss_fi)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_ss_film (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)
    wr_ss_film (i) = 1.0e2  * wr0  

    !-----------------------------------------------------------------------
    ! Sea-salt aerosol from CLASSIC - jet mode
    aerosol_type = ip_seasalt_jet

    n_tot = sea_salt_jet (i)

    ln_r0    = LOG (r_bar_ss_je)
    ln_sigma = LOG (sigma_ss_je)

    CALL grow_particles (aerosol_type, rh_frac (i), ln_r0, ln_sigma,    &
                         ln_wr0, ln_wsigma)

    wr0     = EXP (ln_wr0)
    ln_2wr0 = LOG (2*wr0)

    sa_ss_jet (i) = 1.0e-2 * pi * n_tot * EXP (2*ln_2wr0 + 2*(ln_wsigma)**2)
    wr_ss_jet (i) = 1.0e2  * wr0
  END IF

END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE calc_surf_area 

END MODULE calc_surf_area_mod
