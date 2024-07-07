! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
PROGRAM calc_pm_params
!
! This is a standalone program which was used to
! calculate the numbers which need to be parameters
! in the pm diagnostics subroutine and output statements
! to declare them

! MMethod:
! Calculation of the PM10 and PM2.5 mass concentrations (in
! microgram/m3) by using the cumulative volume (from 0 to
! cut-off diameter) of a log-normal, with median radius
! (or gemetric mean radius) r_bar, geometric standard deviation
! sigma and total number Ntot.
! The contributions by the different aerosol species and modes
! are added up to get the total PM10 & PM2.5 concs (microgram/m3).
!
! Formulation derived from eqs. (7.51) & (7.51a) of Seinfeld
! and Pandis (1998) using method of (7.39). Also, checked with
! "The Mechanics of Inhaled Pharmaceutical Aerosols: An
! Introduction", by Warren A. Finlay, Academic Press, 2001
!
! (1) Cumulative_volume for particles with diameter between
! 0 and d_cut (assuming spherical particles):
!        V(d_cut) = (A*Ntot/2.0) + (A*Ntot/2.0) *
!               erf((alog(d_cut)-B)/(sqrt(2.0)*alog(sigma)))
! with
!        A = (pi/6.0)*exp(3.0*alog(dbar)+4.5*(alog(sigma))**2)
!        B = alog (volume median diameter) =
!          = alog(d_bar)+3.0*(alog(sigma))**2


IMPLICIT NONE

! Pi
REAL, PARAMETER :: Pi                 = 3.14159265358979323846
! R IS GAS CONSTANT FOR DRY AIR
REAL, PARAMETER  :: r      = 287.05


! air density, total nr. aerosol conc. & cut-off diams.
REAL :: rho_air, n_tot, d_cut1, d_cut2
! Variables for sulphate
REAL :: rho_su,r_bar_su_ait, sigma_su_ait, denom_su_ait,                       &
        a_param_su_ait, b_param_su_ait,                                        &
        r_bar_su_acc, sigma_su_acc, denom_su_acc,                              &
        a_param_su_acc, b_param_su_acc,                                        &
        erf_su_ait_d1, erf_su_ait_d2, erf_su_acc_d1, erf_su_acc_d2

! Variables for black carbon
REAL :: rho_bc, r_bar_bc_fr, sigma_bc_fr, denom_bc_fr,                         &
        a_param_bc_fr, b_param_bc_fr,                                          &
        r_bar_bc_ag, sigma_bc_ag, denom_bc_ag,                                 &
        a_param_bc_ag, b_param_bc_ag,                                          &
        erf_bc_fr_d1, erf_bc_fr_d2, erf_bc_ag_d1, erf_bc_ag_d2

! Variables for biomass burning aerosol
REAL :: rho_bb, r_bar_bb_fr, sigma_bb_fr, denom_bb_fr,                         &
        a_param_bb_fr, b_param_bb_fr,                                          &
        r_bar_bb_ag, sigma_bb_ag, denom_bb_ag,                                 &
        a_param_bb_ag, b_param_bb_ag,                                          &
        erf_bb_fr_d1, erf_bb_fr_d2, erf_bb_ag_d1, erf_bb_ag_d2

! Variables for OCFF
REAL :: rho_ocff, r_bar_ocff_fr, sigma_ocff_fr, denom_ocff_fr,                 &
        a_param_ocff_fr, b_param_ocff_fr,                                      &
        r_bar_ocff_ag, sigma_ocff_ag, denom_ocff_ag,                           &
        a_param_ocff_ag, b_param_ocff_ag,                                      &
        erf_ocff_fr_d1, erf_ocff_fr_d2, erf_ocff_ag_d1, erf_ocff_ag_d2

! Variables for SOA
REAL :: rho_soa, r_bar_soa, sigma_soa, denom_soa,                              &
        a_param_soa, b_param_soa,                                              &
        erf_soa_d1, erf_soa_d2

! Variables for sea salt aerosol
REAL :: rho_ss, r_bar_ss_fi, sigma_ss_fi,                                      &
        a_param_ss_fi, b_param_ss_fi,                                          &
        r_bar_ss_je, sigma_ss_je,                                              &
        a_param_ss_je, b_param_ss_je,                                          &
        erf_ss_fi_d1, erf_ss_fi_d2, erf_ss_je_d1, erf_ss_je_d2
! Variables for nitrate
REAL ::   rho_ni, r_bar_ni_acc, sigma_ni_acc, denom_ni_acc,                    &
        a_param_ni_acc, b_param_ni_acc,                                        &
        erf_ni_acc_d1, erf_ni_acc_d2

!     Cut-off aerodyn. diameters for PM10 and PM2.5 (metres)
d_cut1 = 10.0e-06
d_cut2 =  2.5e-06

! Calculate parameters for each aerosol type
! (these could be done once and some are parameters not variables)
! Sulphate

rho_su         = 1769.0                !all densities in kg/m3
r_bar_su_ait   = 0.0065e-6             !all radii in metres
sigma_su_ait   = 1.3
r_bar_su_acc   = 0.095e-6
sigma_su_acc   = 1.4
denom_su_ait   = 4.0*pi*rho_su*(r_bar_su_ait**3)                               &
  *EXP(4.5*(ALOG(sigma_su_ait))**2)
a_param_su_ait = pi*(1.0/6.0)                                                  &
  *EXP(3.0*ALOG(r_bar_su_ait*2.0)+4.5*(ALOG(sigma_su_ait))**2)
b_param_su_ait = ALOG(r_bar_su_ait*2.0)+3.0*(ALOG(sigma_su_ait))**2
denom_su_acc   = 4.0*pi*rho_su*(r_bar_su_acc**3)                               &
  *EXP(4.5*(ALOG(sigma_su_acc))**2)
a_param_su_acc = pi*(1.0/6.0)                                                  &
  *EXP(3.0*ALOG(r_bar_su_acc*2.0)+4.5*(ALOG(sigma_su_acc))**2)
b_param_su_acc = ALOG(r_bar_su_acc*2.0)+3.0*(ALOG(sigma_su_acc))**2

erf_su_ait_d1  = ERF ((ALOG(d_cut1) - b_param_su_ait)                          &
  / (SQRT(2.0) * ALOG(sigma_su_ait)))
erf_su_ait_d2  = ERF ((ALOG(d_cut2) - b_param_su_ait)                          &
  / (SQRT(2.0) * ALOG(sigma_su_ait)))

erf_su_acc_d1  = ERF ((ALOG(d_cut1) - b_param_su_acc)                          &
  / (SQRT(2.0) * ALOG(sigma_su_acc)))
erf_su_acc_d2  = ERF ((ALOG(d_cut2) - b_param_su_acc)                          &
  / (SQRT(2.0) * ALOG(sigma_su_acc)))


! Black carbon

rho_bc        = 1900.0
r_bar_bc_fr   = 0.04e-6
sigma_bc_fr   = 2.0
r_bar_bc_ag   = 0.04e-6
sigma_bc_ag   = 2.0
denom_bc_fr   = 4.0*pi*rho_bc*(r_bar_bc_fr**3)                                 &
  *EXP(4.5*(ALOG(sigma_bc_fr))**2)
a_param_bc_fr = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_bc_fr*2.0)+4.5*(ALOG(sigma_bc_fr))**2)
b_param_bc_fr = ALOG(r_bar_bc_fr*2.0)+3.0*(ALOG(sigma_bc_fr))**2
denom_bc_ag   = 4.0*pi*rho_bc*(r_bar_bc_ag**3)                                 &
  *EXP(4.5*(ALOG(sigma_bc_ag))**2)
a_param_bc_ag = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_bc_ag*2.0)+4.5*(ALOG(sigma_bc_ag))**2)
b_param_bc_ag = ALOG(r_bar_bc_ag*2.0)+3.0*(ALOG(sigma_bc_ag))**2

erf_bc_fr_d1  = ERF ((ALOG(d_cut1) - b_param_bc_fr)                            &
  / (SQRT(2.0) * ALOG(sigma_bc_fr)))
erf_bc_fr_d2  = ERF ((ALOG(d_cut2) - b_param_bc_fr)                            &
  / (SQRT(2.0) * ALOG(sigma_bc_fr)))

erf_bc_ag_d1  = ERF ((ALOG(d_cut1) - b_param_bc_ag)                            &
  / (SQRT(2.0) * ALOG(sigma_bc_ag)))
erf_bc_ag_d2  = ERF ((ALOG(d_cut2) - b_param_bc_ag)                            &
  / (SQRT(2.0) * ALOG(sigma_bc_ag)))


! Biomass-burning aerosol

rho_bb        = 1350.0
r_bar_bb_fr   = 0.1e-6
sigma_bb_fr   = 1.3
r_bar_bb_ag   = 0.12e-6
sigma_bb_ag   = 1.3
denom_bb_fr   = 4.0*pi*rho_bb*(r_bar_bb_fr**3)                                 &
  *EXP(4.5*(ALOG(sigma_bb_fr))**2)
a_param_bb_fr = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_bb_fr*2.0)+4.5*(ALOG(sigma_bb_fr))**2)
b_param_bb_fr = ALOG(r_bar_bb_fr*2.0)+3.0*(ALOG(sigma_bb_fr))**2
denom_bb_ag   = 4.0*pi*rho_bb*(r_bar_bb_ag**3)                                 &
  *EXP(4.5*(ALOG(sigma_bb_ag))**2)
a_param_bb_ag = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_bb_ag*2.0)+4.5*(ALOG(sigma_bb_ag))**2)
b_param_bb_ag = ALOG(r_bar_bb_ag*2.0)+3.0*(ALOG(sigma_bb_ag))**2

erf_bb_fr_d1  = ERF ((ALOG(d_cut1) - b_param_bb_fr)                            &
  / (SQRT(2.0) * ALOG(sigma_bb_fr)))
erf_bb_fr_d2  = ERF ((ALOG(d_cut2) - b_param_bb_fr)                            &
  / (SQRT(2.0) * ALOG(sigma_bb_fr)))

erf_bb_ag_d1  = ERF ((ALOG(d_cut1) - b_param_bb_ag)                            &
  / (SQRT(2.0) * ALOG(sigma_bb_ag)))
erf_bb_ag_d2  = ERF ((ALOG(d_cut2) - b_param_bb_ag)                            &
  / (SQRT(2.0) * ALOG(sigma_bb_ag)))


! Organic Carbon (from fossil fuels) aerosol

rho_ocff        = 1350.0
r_bar_ocff_fr   = 0.1e-6
sigma_ocff_fr   = 1.3
r_bar_ocff_ag   = 0.12e-6
sigma_ocff_ag   = 1.3
denom_ocff_fr   = 4.0*pi*rho_ocff*(r_bar_ocff_fr**3)                           &
  *EXP(4.5*(ALOG(sigma_ocff_fr))**2)
a_param_ocff_fr = pi*(1.0/6.0)                                                 &
  *EXP(3.0*ALOG(r_bar_ocff_fr*2.0)+4.5*(ALOG(sigma_ocff_fr))**2)
b_param_ocff_fr = ALOG(r_bar_ocff_fr*2.0)                                      &
  +3.0*(ALOG(sigma_ocff_fr))**2
denom_ocff_ag   = 4.0*pi*rho_ocff*(r_bar_ocff_ag**3)                           &
  *EXP(4.5*(ALOG(sigma_ocff_ag))**2)
a_param_ocff_ag = pi*(1.0/6.0)                                                 &
  *EXP(3.0*ALOG(r_bar_ocff_ag*2.0)+4.5*(ALOG(sigma_ocff_ag))**2)
b_param_ocff_ag = ALOG(r_bar_ocff_ag*2.0)                                      &
  +3.0*(ALOG(sigma_ocff_ag))**2

erf_ocff_fr_d1  = ERF ((ALOG(d_cut1) - b_param_ocff_fr)                        &
  / (SQRT(2.0) * ALOG(sigma_ocff_fr)))
erf_ocff_fr_d2  = ERF ((ALOG(d_cut2) - b_param_ocff_fr)                        &
  / (SQRT(2.0) * ALOG(sigma_ocff_fr)))

erf_ocff_ag_d1  = ERF ((ALOG(d_cut1) - b_param_ocff_ag)                        &
  / (SQRT(2.0) * ALOG(sigma_ocff_ag)))
erf_ocff_ag_d2  = ERF ((ALOG(d_cut2) - b_param_ocff_ag)                        &
  / (SQRT(2.0) * ALOG(sigma_ocff_ag)))


! Biogenic Secondary Organic Aerosol

rho_soa     = 1300.0
r_bar_soa   = 0.095e-6
sigma_soa   = 1.5
denom_soa   = 4.0*pi*rho_soa*(r_bar_soa**3)                                    &
  *EXP(4.5*(ALOG(sigma_soa))**2)
a_param_soa = pi*(1.0/6.0)                                                     &
  *EXP(3.0*ALOG(r_bar_soa*2.0)+4.5*(ALOG(sigma_soa))**2)
b_param_soa = ALOG(r_bar_soa*2.0)+3.0*(ALOG(sigma_soa))**2

erf_soa_d1  = ERF ((ALOG(d_cut1) - b_param_soa)                                &
  / (SQRT(2.0) * ALOG(sigma_soa)))
erf_soa_d2  = ERF ((ALOG(d_cut2) - b_param_soa)                                &
  / (SQRT(2.0) * ALOG(sigma_soa)))


! Sea-salt aerosol

rho_ss        = 2165.0
r_bar_ss_fi   = 0.1e-6
sigma_ss_fi   = 1.9
r_bar_ss_je   = 1.0e-6
sigma_ss_je   = 2.0
a_param_ss_fi = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_ss_fi*2.0)+4.5*(ALOG(sigma_ss_fi))**2)
b_param_ss_fi = ALOG(r_bar_ss_fi*2.0)+3.0*(ALOG(sigma_ss_fi))**2
a_param_ss_je = pi*(1.0/6.0)                                                   &
  *EXP(3.0*ALOG(r_bar_ss_je*2.0)+4.5*(ALOG(sigma_ss_je))**2)
b_param_ss_je = ALOG(r_bar_ss_je*2.0)+3.0*(ALOG(sigma_ss_je))**2

erf_ss_fi_d1  = ERF ((ALOG(d_cut1) - b_param_ss_fi)                            &
  / (SQRT(2.0) * ALOG(sigma_ss_fi)))
erf_ss_fi_d2  = ERF ((ALOG(d_cut2) - b_param_ss_fi)                            &
  / (SQRT(2.0) * ALOG(sigma_ss_fi)))

erf_ss_je_d1  = ERF ((ALOG(d_cut1) - b_param_ss_je)                            &
  / (SQRT(2.0) * ALOG(sigma_ss_je)))
erf_ss_je_d2  = ERF ((ALOG(d_cut2) - b_param_ss_je)                            &
  / (SQRT(2.0) * ALOG(sigma_ss_je)))


! Ammonium nitrate (there is only an accumulation mode)

rho_ni          = 1725.0
r_bar_ni_acc    = 0.095e-6
sigma_ni_acc    = 1.4
denom_ni_acc    = 4.0*pi*rho_ni*(r_bar_ni_acc**3)                              &
  *EXP(4.5*(ALOG(sigma_ni_acc))**2)
a_param_ni_acc  = pi*(1.0/6.0)                                                 &
  *EXP(3.0*ALOG(r_bar_ni_acc*2.0)+4.5*(ALOG(sigma_ni_acc))**2)
b_param_ni_acc  = ALOG(r_bar_ni_acc*2.0)+3.0*(ALOG(sigma_ni_acc))**2

erf_ni_acc_d1  = ERF ((ALOG(d_cut1) - b_param_ni_acc)                          &
  / (SQRT(2.0) * ALOG(sigma_ni_acc)))
erf_ni_acc_d2  = ERF ((ALOG(d_cut2) - b_param_ni_acc)                          &
  / (SQRT(2.0) * ALOG(sigma_ni_acc)))

WRITE(6,*) '! Parameters for sulphate '
WRITE(6,*) 'REAL, PARAMETER :: denom_su_ait    =',denom_su_ait
WRITE(6,*) 'REAL, PARAMETER :: a_param_su_ait  =',a_param_su_ait
WRITE(6,*) 'REAL, PARAMETER :: erf_su_ait_d1   =',erf_su_ait_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_su_ait_d2   =',erf_su_ait_d2
WRITE(6,*) 'REAL, PARAMETER :: denom_su_acc    =',denom_su_acc
WRITE(6,*) 'REAL, PARAMETER :: a_param_su_acc  =',a_param_su_acc
WRITE(6,*) 'REAL, PARAMETER :: erf_su_acc_d1   =',erf_su_acc_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_su_acc_d2   =',erf_su_acc_d2
WRITE(6,*) '! Parameters for Black Carbon (soot)'
WRITE(6,*) 'REAL, PARAMETER :: denom_bc_fr     =',denom_bc_fr
WRITE(6,*) 'REAL, PARAMETER :: a_param_bc_fr   =',a_param_bc_fr
WRITE(6,*) 'REAL, PARAMETER :: erf_bc_fr_d1    =',erf_bc_fr_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_bc_fr_d2    =',erf_bc_fr_d2
WRITE(6,*) 'REAL, PARAMETER :: denom_bc_ag     =',denom_bc_ag
WRITE(6,*) 'REAL, PARAMETER :: a_param_bc_ag   =',a_param_bc_ag
WRITE(6,*) 'REAL, PARAMETER :: erf_bc_ag_d1    =',erf_bc_ag_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_bc_ag_d2    =',erf_bc_ag_d2
WRITE(6,*) '! Parameters for biomass burning aerosol'
WRITE(6,*) 'REAL, PARAMETER :: denom_bb_fr     =',denom_bb_fr
WRITE(6,*) 'REAL, PARAMETER :: a_param_bb_fr   =',a_param_bb_fr
WRITE(6,*) 'REAL, PARAMETER :: erf_bb_fr_d1    =',erf_bb_fr_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_bb_fr_d2    =',erf_bb_fr_d2
WRITE(6,*) 'REAL, PARAMETER :: denom_bb_ag     =',denom_bb_ag
WRITE(6,*) 'REAL, PARAMETER :: a_param_bb_ag   =',a_param_bb_ag
WRITE(6,*) 'REAL, PARAMETER :: erf_bb_ag_d1    =',erf_bb_ag_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_bb_ag_d2    =',erf_bb_ag_d2
WRITE(6,*) '! Parameters for OCFF'
WRITE(6,*) 'REAL, PARAMETER :: denom_ocff_fr   =',denom_ocff_fr
WRITE(6,*) 'REAL, PARAMETER :: a_param_ocff_fr =',a_param_ocff_fr
WRITE(6,*) 'REAL, PARAMETER :: erf_ocff_fr_d1  =',erf_ocff_fr_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_ocff_fr_d2  =',erf_ocff_fr_d2
WRITE(6,*) 'REAL, PARAMETER :: denom_ocff_ag   =',denom_ocff_ag
WRITE(6,*) 'REAL, PARAMETER :: a_param_ocff_ag =',a_param_ocff_ag
WRITE(6,*) 'REAL, PARAMETER :: erf_ocff_ag_d1  =',erf_ocff_ag_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_ocff_ag_d2  =',erf_ocff_ag_d2
WRITE(6,*) '! Parameters for biogenic aerosol (SOA)'
WRITE(6,*) 'REAL, PARAMETER :: denom_soa       =',denom_soa
WRITE(6,*) 'REAL, PARAMETER :: a_param_soa     =',a_param_soa
WRITE(6,*) 'REAL, PARAMETER :: erf_soa_d1      =',erf_soa_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_soa_d2      =',erf_soa_d2
WRITE(6,*) '! Parameters for sea-salt'
WRITE(6,*) 'REAL, PARAMETER :: a_param_ss_fi   =',a_param_ss_fi
WRITE(6,*) 'REAL, PARAMETER :: erf_ss_fi_d1    =',erf_ss_fi_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_ss_fi_d2    =',erf_ss_fi_d2
WRITE(6,*) 'REAL, PARAMETER :: a_param_ss_je   =',a_param_ss_je
WRITE(6,*) 'REAL, PARAMETER :: erf_ss_je_d1    =',erf_ss_je_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_ss_je_d2    =',erf_ss_je_d2
WRITE(6,*) '! Parameters for nitrate'
WRITE(6,*) 'REAL, PARAMETER :: denom_ni_acc    =',denom_ni_acc
WRITE(6,*) 'REAL, PARAMETER :: a_param_ni_acc  =',a_param_ni_acc
WRITE(6,*) 'REAL, PARAMETER :: erf_ni_acc_d1   =',erf_ni_acc_d1
WRITE(6,*) 'REAL, PARAMETER :: erf_ni_acc_d2   =',erf_ni_acc_d2

END PROGRAM calc_pm_params
