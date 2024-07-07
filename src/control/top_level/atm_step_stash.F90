! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Contains various chunks of code from atm_step - the purpose of each
! section is indicated at the head of the section - mainly related to
! stash calls
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Subroutine Interface:
SUBROUTINE atm_step_stash( &
errorstatus, flag)

USE atm_step_local
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, tdims_l,        &
                                 udims, udims_s, vdims,          &
                                 vdims_s, wdims

USE rad_input_mod, ONLY: l_radiation

USE jules_snow_mod, ONLY: nsmax
USE jules_surface_types_mod, ONLY: ntype, npft
USE turb_diff_mod, ONLY: l_polar_filter_incs, l_filter_incs
USE g_wave_input_mod, ONLY: l_gwd, l_use_ussp
USE bl_option_mod, ONLY: i_bl_vn, i_bl_vn_0
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY: sf
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE dynamics_input_mod, ONLY: NumCycles

USE jules_hydrology_mod, ONLY: l_hydrology
USE ukca_option_mod, ONLY: l_ukca
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain, l_rain
USE electric_inputs_mod, ONLY: l_use_electric
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE cv_run_mod, ONLY: l_param_conv ! Convection scheme switch
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, land_field, len1_lookup,                &
    len_dumphist, len_fixhd, len_tot, model_levels,                    &
    mpp_len1_lookup, n_cca_lev, n_obj_d1_max, n_rows,                  &
    ntiles, river_row_length, river_rows, row_length, rows, sm_levels, &
    theta_off_size, tr_levels, tr_ukca, tr_vars

USE atm_fields_mod, ONLY: q, qcl, qcf, qcf2, qrain, qgraup, theta,     &
                          exner_theta_levels, u, v, w
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


! Subroutine arguments

INTEGER :: errorstatus, flag

! Local variables

REAL :: stashwork0_dummy(1)
            ! STASHwork not defined for sec 0, but required as dummy argument.
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if return code >0
CHARACTER(LEN=*) :: RoutineName

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
PARAMETER (   RoutineName='ATM_STEP_STASH')

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (flag == 1) THEN

  ! Clear workspace for radiation
  DEALLOCATE ( ozone3d )
  ! Deallocate space for tropopause diagnostics
  DEALLOCATE (o3_trop_level)
  DEALLOCATE (o3_trop_height)
  DEALLOCATE (t_trop_level)
  DEALLOCATE (t_trop_height)

  ! Diagnostics STASHed for each section in Atmos_Physics1:
  IF (l_radiation) THEN
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,1,stashwork1,                          &
          errorstatus,cmessage)

    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,2,stashwork2,                          &
            errorstatus,cmessage)
  END IF
  DEALLOCATE (stashwork1)
  DEALLOCATE (stashwork2)

  IF (l_rain) THEN
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,4,stashwork4,                          &
            errorstatus,cmessage)
  END IF
  DEALLOCATE (stashwork4)

  IF (l_use_electric) THEN
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,21,stashwork21,                        &
            errorstatus,cmessage)
  END IF
  DEALLOCATE (stashwork21)

  IF (l_gwd .OR. l_use_ussp) THEN
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,6,stashwork6,                          &
            errorstatus,cmessage)
  END IF
  DEALLOCATE (stashwork6)

  IF (l_emcorr) THEN
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,14,stashwork14,                      &
          errorstatus,cmessage)
  END IF   ! l_emcorr
  DEALLOCATE (stashwork14)


  ! ----------------------------------------------------------------------


ELSE IF (flag == 2) THEN

  ! Apply diagnostics only at last cycle.
  IF ( cycleno == numcycles ) THEN

    ! Diagnostics STASHed for each section in Atmos_Physics2:

    IF ( i_bl_vn /= i_bl_vn_0 ) THEN
      ! DEPENDS ON: stash
      CALL stash(atmos_sm,atmos_im,3,stashwork3,                        &
          errorstatus,cmessage)
    END IF ! i_bl_vn
    DEALLOCATE (stashwork3)

    IF (l_param_conv) THEN
      ! DEPENDS ON: stash
      CALL stash(atmos_sm,atmos_im,5,stashwork5,                      &
          errorstatus,cmessage)
    END IF   ! on l_param_conv
    DEALLOCATE (stashwork5)

    IF (l_hydrology) THEN
      ! DEPENDS ON: stash
      CALL stash(atmos_sm,atmos_im,8,stashwork8,                        &
            errorstatus,cmessage)
    END IF
    DEALLOCATE (stashwork8)

    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,9,stashwork9,                        &
          errorstatus,cmessage)
    DEALLOCATE (stashwork9)

    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,19,stashwork19,                      &
          errorstatus,cmessage)
    DEALLOCATE (stashwork19)

    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,26,stashwork26,                  &
          errorstatus,cmessage)
    DEALLOCATE (stashwork26)

  END IF  ! CycleNo == NumCycles


  ! ----------------------------------------------------------------------

ELSE IF (flag == 3) THEN

  ! Section 13 stash call

  IF ( .NOT. ( (model_type == mt_global) .AND.               &
                    (l_polar_filter_incs .OR. l_filter_incs) ) ) THEN

    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,13,stashwork13,                &
            errorstatus,cmessage)
     ! We also deallocate in atm_step for the else statement of this block.
    DEALLOCATE(stashwork13)

  END IF ! .NOT. ( (model_type == mt_global ) etc

  ! ----------------------------------------------------------------------


ELSE IF (flag == 4) THEN

  ! section 0: extraction of primary variables

  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,0,stashwork0_dummy,                      &
                 errorstatus,cmessage)
  ! Check error condition
  IF (errorstatus >  0) THEN

    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,33,stashwork0_dummy,                     &
                 errorstatus,cmessage)
  ! Check error condition
  IF (errorstatus >  0) THEN

    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  ! Section 34: extraction of UKCA tracer variables

  IF (.NOT. l_ukca) THEN      ! UKCA has its own calls to STASH
    ! DEPENDS ON: stash
    CALL stash(atmos_sm,atmos_im,34,stashwork0_dummy,                   &
                   errorstatus,cmessage)
  END IF    ! L_ukca

  ! Check error condition
  IF (errorstatus >  0) THEN

    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  ! -------------------------------------------------------------------------

  ! Section 39: Extraction of Nudging variables
ELSE IF (flag == 39) THEN
  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,39,stashwork39,                      &
                 errorstatus,cmessage)

  ! Check error condition
  IF (errorstatus >  0) THEN

    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  DEALLOCATE (stashwork39)

END IF ! flag
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_stash
