#if defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: oasis_initialise_2
!
!  Purpose: Second part of initialisation for OASIS
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Coupling
!  -------------------------------------------------------------------

SUBROUTINE oasis_initialise_2 (cmessage )

USE atm_fields_mod
USE atm_fields_bounds_Mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE oasis_atm_data_mod,   ONLY: oasis_couple_ts_ao, oasis_couple_ts_oa,    & 
                          oasis_cpl_ts_get_hyb, oasis_cpl_ts_put_hyb,      & 
                          oasis_cpl_ts_get_hyb_stats,                      &
                          oasis_cpl_ts_put_hyb_stats,                      &
                          oasis_couple_ts_aw, oasis_couple_ts_wa,          & 
                          put_step_ao, get_step_oa,                        &
                          put_step_hyb, get_step_hyb,                      &
                          put_step_hyb_stats, get_step_hyb_stats,          &
                          put_step_aw, get_step_wa,                        &
                          cpl_update_step, prism_timestep

USE coupling_control_mod, ONLY: oasis_couple_freq_ao, oasis_couple_freq_oa, & 
                          oasis_couple_freq_ac, oasis_couple_freq_ca,       & 
                          oasis_couple_freq_ac_stats,                       &
                          oasis_couple_freq_ca_stats,                       &
                          oasis_couple_freq_aw, oasis_couple_freq_wa,       &
                          l_oasis_ocean, l_senior, l_junior

USE nlstgen_mod,          ONLY: secs_per_periodim, steps_per_periodim

! Declarations:
USE submodel_mod, ONLY: atmos_sm
USE io
USE um_types
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE Decomp_DB

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,   &
    n_cca_lev, n_obj_d1_max, sm_levels, st_levels, tpps_ozone_levels,  &
    tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars

USE model_time_mod,         ONLY: secs_per_stepim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!
! ----------------------------------------------------------------------
!
!

! Subroutine arguments:
CHARACTER(LEN=errormessagelength) :: cmessage    ! OUT - Error return message

! Local variables:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_INITIALISE_2'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  5.1.1  Initialise coupling info for OASIS3-MCT coupler
! ----------------------------------------------------------------------

!     Work out the frequency of coupling timesteps

! Atmos -> Ocean exchanges
oasis_couple_ts_ao = steps_per_periodim(1) *             &
      ((oasis_couple_freq_ao(1)*3600)+(oasis_couple_freq_ao(2)*60)) &
       /secs_per_periodim(1)

! Ocean -> Atmos exchanges
oasis_couple_ts_oa = steps_per_periodim(1) *             &
      ((oasis_couple_freq_oa(1)*3600)+(oasis_couple_freq_oa(2)*60)) & 
      /secs_per_periodim(1)

IF (l_senior) THEN
  ! Set switches to indicate whether a coupling exchange is required
  ! in the senior component of an Atmos-Chemistry hybrid model.  

  ! Atmos -> Chemistry exchanges
  oasis_cpl_ts_put_hyb = steps_per_periodim(1) *                    &
      ((oasis_couple_freq_ac(1)*3600)+(oasis_couple_freq_ac(2)*60)) &
      /secs_per_periodim(1)

  ! Atmos -> Chemistry exchanges for stats
  oasis_cpl_ts_put_hyb_stats = steps_per_periodim(1) *              &
      ((oasis_couple_freq_ac_stats(1)*3600)+                        &
       (oasis_couple_freq_ac_stats(2)*60))                          &
       /secs_per_periodim(1)

  ! Chemistry  -> Atmos exchanges
  oasis_cpl_ts_get_hyb = steps_per_periodim(1) *                    &
      ((oasis_couple_freq_ca(1)*3600)+(oasis_couple_freq_ca(2)*60)) &
      /secs_per_periodim(1)

  ! Chemistry  -> Atmos exchanges for stats
  oasis_cpl_ts_get_hyb_stats = steps_per_periodim(1) *              &
      ((oasis_couple_freq_ca_stats(1)*3600)+                        &
       (oasis_couple_freq_ca_stats(2)*60))                          &
      /secs_per_periodim(1)

ELSE IF (l_junior) THEN
  ! Set switches to indicate whether a coupling exchange is required
  ! in the junior component of an Atmos-Chemistry hybrid model.  

  ! Atmos -> Chemistry exchanges
  oasis_cpl_ts_get_hyb = steps_per_periodim(1) *                    &
      ((oasis_couple_freq_ac(1)*3600)+(oasis_couple_freq_ac(2)*60)) &
      /secs_per_periodim(1)

  ! Atmos -> Chemistry exchanges for stats
  oasis_cpl_ts_get_hyb_stats = steps_per_periodim(1) *              &
      ((oasis_couple_freq_ac_stats(1)*3600)+                        &
       (oasis_couple_freq_ac_stats(2)*60))                          &
       /secs_per_periodim(1)

  ! Chemistry  -> Atmos exchanges
  oasis_cpl_ts_put_hyb = steps_per_periodim(1) *                    &
      ((oasis_couple_freq_ca(1)*3600)+(oasis_couple_freq_ca(2)*60)) &
      /secs_per_periodim(1)

  ! Chemistry  -> Atmos exchanges for stats
  oasis_cpl_ts_put_hyb_stats = steps_per_periodim(1) *              &
      ((oasis_couple_freq_ca_stats(1)*3600)+                        &
       (oasis_couple_freq_ca_stats(2)*60))                          &
      /secs_per_periodim(1)

END IF

! Atmos -> Wave exchanges
oasis_couple_ts_aw = steps_per_periodim(1) *             &
      ((oasis_couple_freq_aw(1)*3600)+(oasis_couple_freq_aw(2)*60)) &
      /secs_per_periodim(1)

! Wave  -> Atmos exchanges
oasis_couple_ts_wa = steps_per_periodim(1) *             &
      ((oasis_couple_freq_wa(1)*3600)+(oasis_couple_freq_wa(2)*60)) &
      /secs_per_periodim(1)

! Initialise the coupling logicals
put_step_ao       =.FALSE.
get_step_oa       =.FALSE.
put_step_hyb      =.FALSE.
get_step_hyb      =.FALSE.
put_step_hyb_stats=.FALSE.
get_step_hyb_stats=.FALSE.
put_step_aw       =.FALSE.
get_step_wa       =.FALSE.
cpl_update_step   =.FALSE.

!     Set timestep for use in OASIS calls to be
!     real version of the standard atmos timestep

prism_timestep = 1.0*secs_per_stepim(atmos_sm)


! ----------------------------------------------------------------------
!  8. If coupled model, initialise addresses of coupling fields.
! ----------------------------------------------------------------------
IF (l_oasis_ocean) THEN

  ! DEPENDS ON: oasis_inita2o
  CALL oasis_inita2o(cmessage)

END IF

IF ( (l_senior) .OR. (l_junior) ) THEN

  ! If stats coupling frequencies are not set, set them to standard 
  ! coupling frequencies
  IF (oasis_cpl_ts_get_hyb_stats == 0)                              &
    oasis_cpl_ts_get_hyb_stats = oasis_cpl_ts_get_hyb
  IF (oasis_cpl_ts_put_hyb_stats == 0)                              &
    oasis_cpl_ts_put_hyb_stats = oasis_cpl_ts_put_hyb

END IF

! DEPENDS ON: OASIS3_GRID
CALL oasis3_grid(.TRUE.)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis_initialise_2
#endif
