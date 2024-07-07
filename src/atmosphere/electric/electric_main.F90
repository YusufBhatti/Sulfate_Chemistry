! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Control routine for the electric scheme.
! Calls routines to do the following steps:
!--------------------------------------------
! Initialise all variables
! Define where storms exist in the model
! Calculate flash rate for each storm
! Distribute the flashes as required in the model
! Calculate any diagnostics
!--------------------------------------------

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

MODULE electric_main_mod

USE atm_fields_bounds_mod,    ONLY: tdims, wdims_s
USE electric_init_mod,        ONLY: electric_init
USE define_storm_mod,         ONLY: define_storm
USE distribute_flash_mod,     ONLY: distribute_flash
#if !defined(LFRIC)
USE diagnostics_electric_mod, ONLY: diagnostics_electric
#endif
USE flash_rate_mod,           ONLY: flash_rate
USE calc_wp_mod,              ONLY: calc_wp
USE mphys_bypass_mod,         ONLY: qcf2_idims_start, qcf2_idims_end,         &
                                    qcf2_jdims_start, qcf2_jdims_end,         &
                                    qcf2_kdims_start, qcf2_kdims_end

! Dr Hook modules
USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ELECTRIC_MAIN_MOD'

CONTAINS

SUBROUTINE electric_main( qcf, qcf2, qgraup, rhodz_dry, rhodz_moist, t, w,     &
                          at_extremity, stashwork, flash_pot )

IMPLICIT NONE

REAL, INTENT(IN) :: qcf(         tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
REAL, INTENT(IN) :: qcf2(     qcf2_idims_start : qcf2_idims_end,               &
                              qcf2_jdims_start : qcf2_jdims_end,               &
                              qcf2_kdims_start : qcf2_kdims_end )
REAL, INTENT(IN) :: qgraup(      tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
REAL, INTENT(IN) :: rhodz_dry(   tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
REAL, INTENT(IN) :: rhodz_moist( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
REAL, INTENT(IN) :: t(           tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
REAL, INTENT(IN) :: w(         wdims_s%i_start : wdims_s%i_end,                &
                               wdims_s%j_start : wdims_s%j_end,                &
                               wdims_s%k_start : wdims_s%k_end )

LOGICAL ::  at_extremity(4)
! Indicates if this processor is at north,
! south, east or west of the processor grid)


REAL, INTENT(INOUT) :: stashwork(*)
! STASH workspace for electric section

REAL, INTENT(INOUT) :: flash_pot  ( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )

!------------------------------------------------------------------------------
!Local variables
REAL :: rhodz(              tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end,                       &
                                        1 : tdims%k_end )
REAL :: gwp  (              tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: tiwp (              tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: qcf_tot(            tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end,                       &
                                        1 : tdims%k_end )
REAL :: flash(              tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: flash_pot1d(        tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: qgtot(              tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: fr1_mc(             tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: fr2_mc(             tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
LOGICAL :: storm_field(     tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
REAL :: num_flashes(        tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
INTEGER :: nspts ! number of storm points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ELECTRIC_MAIN'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=====================================================================
! 1. Initialise all parameters and determine the values of parameters
!    required for lightning flash rate
!=====================================================================

CALL electric_init( storm_field, qcf, qcf2, qgraup, qcf_tot, rhodz_dry,       &
                    rhodz_moist, flash_pot, rhodz, flash, qgtot, flash_pot1d, &
                    fr1_mc, fr2_mc )

!Graupel water path (GWP)
CALL calc_wp( qgraup,  rhodz, gwp  )

!Ice water path (IWP)
CALL calc_wp( qcf_tot, rhodz, tiwp )

!=====================================================================
! 2. Locate where storms are present within the model fields
!=====================================================================
! Define where storm points are
CALL define_storm( gwp, storm_field, nspts )

!=====================================================================
! 3. Calculate the flash rate only when there are storm points
!=====================================================================
IF (nspts > 0) CALL flash_rate( storm_field, qgraup, t, w, rhodz, gwp, tiwp,   &
                                flash, fr1_mc, fr2_mc)

!=====================================================================
! 4. Calculate number of lightning flashes this timestep and
!    distribute excess potential for prognosis
!=====================================================================
CALL distribute_flash( qgraup, qgtot, flash, flash_pot1d, flash_pot,           &
                       num_flashes )

!=====================================================================
! 5. Output any diagnostics requested
!=====================================================================
#if !defined(LFRIC)
CALL diagnostics_electric( flash, storm_field, gwp, tiwp, num_flashes,         &
                           fr1_mc, fr2_mc, stashwork, at_extremity )
#endif

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE electric_main

END MODULE electric_main_mod
