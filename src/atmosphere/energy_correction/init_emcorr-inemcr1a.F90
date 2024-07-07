! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Energy Correction

MODULE init_emcorr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INIT_EMCORR_MOD'
CONTAINS

SUBROUTINE init_emcorr(                                           &
              icode,cmessage)
!
!  Subroutine: INIT_EMCORR
!
!  Purpose: Interface routine required to pass super arrays down into
!           ENG_MASS_CAL which initialises the energy correction.
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
! ----------------------------------------------------------------------
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook
USE atm_fields_mod,        ONLY: u, v, w, rho, theta, q, qcl,          &
                                 qcf, qcf2, qrain, qgraup, p,          &
                                 exner_theta_levels, pstar
USE eng_mass_diag_mod,     ONLY: eng_mass_diag
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE wet_to_dry_n_calc_mod
USE eng_corr_inputs_mod,   ONLY: lemq_print, a_energysteps
USE dump_headers_mod, ONLY: rh_tot_mass_init, rh_tot_energy_init,      &
                            rh_energy_corr, rh_tot_m_init, a_realhd
USE missing_data_mod, ONLY: rmdi  
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
!  Arguments
!
!  Configuration-dependent sizes and arrays


INTEGER :: icode             ! Work - Internal return code
CHARACTER(LEN=errormessagelength) :: cmessage  ! Work - Internal error message
!  NOTE icode and cmessage are not altered by this subroutine
!       returned unchanged.
!-------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------
LOGICAL ::                                                        &
  lmass_cor0,                                                     &
                  ! switch for mass correction at T+0
  lmoist_cor0     ! switch for moist correction at T+0

REAL ::                                                           &
  dummy           ! dummy initial mass

REAL ::                                                           &
 weight1,                                                         &
                ! weight arrays for wet to dry
                ! conversion calculation
 weight2,                                                         &
 weight3


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_EMCORR'

!-------------------------------------------------------------------
! Only set initial energy and mass if the header values in the start
! dump are set to missing data otherwise the model will start from
! the values in the dump.
! Values in the dump will come from a previous climate integration.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (a_realhd(rh_tot_mass_init  ) == rmdi .OR.                     &
    a_realhd(rh_tot_energy_init) == rmdi .OR.                     &
    a_realhd(rh_energy_corr)     == rmdi) THEN

  a_realhd(rh_tot_mass_init) = 0.0
  a_realhd(rh_tot_energy_init) = 0.0

  lmass_cor0 =.FALSE.       ! do not attempt mass correction
  lmoist_cor0=.FALSE.       ! do not attempt moist correction
  dummy = 0.0

  CALL wet_to_dry_n_calc( q,                  &
                          qcl,                &
                          qcf,                &
                          qcf2,               &
                          qrain,              &
                          qgraup )

  CALL eng_mass_diag (                                            &
                    delta_lambda,delta_phi,                       &
                    theta,                                        &
                    u,                                            &
                    v,                                            &
                    w,                                            &
                    rho,                                          &
                    q,                                            &
                    qcl,                                          &
                    qcf,                                          &
                    wet_to_dry_n,                                 &
                    exner_theta_levels,                           &
  !  passing pstar as a dummy array
                            pstar,                                &
                            dummy,dummy,lmass_cor0,lmoist_cor0,   &
                            lemq_print,                           &
                            a_energysteps,                        &
                            a_realhd(rh_tot_energy_init),         &
                            a_realhd(rh_tot_mass_init),           &
                            a_realhd(rh_tot_m_init) )


  ! INITIAL RATE OF ENERGY CORRECTION TO ZERO

  a_realhd(rh_energy_corr) = 0.0

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_emcorr

END MODULE init_emcorr_mod
