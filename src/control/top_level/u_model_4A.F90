! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: U_MODEL_4A (ENDGAME VERSION) -------------------------
!
!    Purpose: High level control program for the Unified Model
!             (master routine).  Calls lower level control routines
!             according to top level switch settings. Called by
!             top level routine UMSHELL which provides dimension sizes
!             for dynamic allocation of data arrays.
!
!    Programming standard: UM Doc Paper 3, version 8.3
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Top Level

SUBROUTINE u_model_4a(dump_unit)

! D1 replacement module

USE atm_fields_mod
USE atm_fields_bounds_mod

! OASIS Modules
USE oasis_atm_data_mod


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE io, ONLY: io_timestep
USE model_file, ONLY: storeAllLookups, SynchroniseAll
USE um_types

! ensure that module variables cannot go out of scope
USE horiz_grid_mod
USE ref_pro_mod
USE departure_pts_mod
USE fields_rhs_mod
USE metric_terms_mod
USE helmholtz_const_matrix_mod
USE gravity_mod
USE eg_parameters_mod
USE wet_to_dry_n_calc_mod, ONLY: wet_to_dry_n

USE Control_Max_Sizes
USE decomp_DB
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE get_env_var_mod, ONLY: get_env_var

USE coupling_control_mod, ONLY: l_oasis, l_oasis_ocean, l_senior, l_junior
USE river_inputs_mod, ONLY: l_rivers
USE ppxlook_mod
USE acp_namel_mod, ONLY: l_ac

USE submodel_mod, ONLY: atmos_im
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim

USE nlstcall_mod, ONLY: model_basis_time, &
                         lpp, &
                         lnc, &
                         ldump, &
                         lmean, &
                         lprint, &
                         lexit, &
                         lancillary, &
                         lboundary, &
                         ltimer, &
                         lstashdumptimer

USE up_ancil_mod, ONLY: up_ancil

USE history, ONLY: checkpoint_dump_im, blank_file_name

USE temphist_mod, ONLY: temphist, write_interim, write_temporary

USE model_time_mod, ONLY: secs_per_stepim, stepim

USE num_obs_mod, ONLY: ac_num_obs_max, ac_tot_obs_size_max

USE land_soil_dimensions_mod, ONLY: deallocate_land_soil_arrays

USE close_unneeded_stash_files_mod, ONLY: close_unneeded_stash_files,          &
    close_unneeded_stash_files_freq, check_close_unneeded_stash_files

USE run_info, ONLY: start_time 
USE io_configuration_mod, ONLY: print_runtime_info
USE turb_diff_mod, ONLY: l_print_pe
USE timestep_mod, ONLY: timestep_number

USE ncfile_reinit_mod, ONLY: ncfile_reinit

USE del_hist_mod, ONLY: del_hist

! ----------------------------------------------------------------------+-------

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
  
REAL :: time_end_run

!    Interface and arguments: ------------------------------------------

INTEGER, INTENT(IN) :: dump_unit  ! unit attached to input data file 
                                  ! (astart or checkpoint_dump_im)
! ----------------------------------------------------------------------
!
INTEGER :: obs_flag_len,obs_len
INTEGER, ALLOCATABLE :: obs_flag(:)
REAL,    ALLOCATABLE :: obs(:)
!
!  Local variables
!
INTEGER :: i_unit         ! Work - File unit number
INTEGER :: internal_model    ! Work - Internal model identifier
INTEGER :: submodel          ! Work - Submodel id for dump partition
INTEGER :: ngroup            ! Work - Number of steps in "group"
INTEGER :: meanlev           ! Work - Mean level indicator
INTEGER :: iabort            ! Work - Internal return code
INTEGER :: i_step            ! Work - Loop counter over timesteps
INTEGER :: g_theta_field_size                                        &
                             ! Sizes for MPP dynamic allocation
       ,g_imtjmt          ! in A-O coupling routines
!
! River routing
INTEGER :: g_river_field_size   ! Sizes for MPP dynamic allocation
!
LOGICAL :: lexitNOW          ! Work - Immediate exit indicator
LOGICAL, PARAMETER :: called_from_u_model = .FALSE.
INTEGER :: co2_dima,                                                 &
                             ! CO2 array dimensions
        co2_dimo,                                                 &
        co2_dimo2
INTEGER :: dms_dima,                                                 &
                             ! DMS array dimensions
        dms_dimo,                                                 &
        dms_dimo2
INTEGER :: info   ! Return code from GCom routines

! 3-D fields of species to be passed down to radiation
INTEGER, PARAMETER :: ngrgas = 8
INTEGER, SAVE :: grgas_addr(ngrgas)

#if defined(IBM_XL_FORTRAN)
INTEGER                      :: pos               ! string position
INTEGER                      :: length            ! length of string contents
CHARACTER (LEN=256)          :: c_xlfrteopts      ! XLFRTEOPTS content
CHARACTER (LEN=*), PARAMETER :: c_disable='buffering=disable_all'
#endif
! Error reporting
INTEGER ::    icode       ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*) :: RoutineName

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
PARAMETER (   RoutineName='U_MODEL_4A')

! ENDGAME-only declarations


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

icode=0
cmessage=''

! ----------------------------------------------------------------------
!  0. Start Timer call for U_MODEL_4A 
!
IF (ltimer) CALL timer('U_MODEL_4A ',5)


! Compress stashmaster to only the active items
CALL compress_atmos_stashmaster()

! ----------------------------------------------------------------------
!  1. General initialisation of control and physical data blocks
!

icode=0

CALL timer('INITIAL ',5)

! DEPENDS ON: initial_4A
CALL initial_4a(                                                  &
     ngrgas,grgas_addr,                                           &
     internal_model,submodel,ngroup,meanlev, dump_unit)
IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)

CALL timer('INITIAL ',6)

!  Allocate AC scheme arrays using sizes from AC_INIT
IF (l_ac) THEN
  obs_flag_len = ac_num_obs_max
  ALLOCATE (obs_flag(obs_flag_len))
  ! 2048 gives enough space in WORK in RDOBS2 when no or very few obs.
  obs_len = ac_tot_obs_size_max+2048
  ALLOCATE (obs(obs_len))
  WRITE(umMessage,*)'U_MODEL_4A - OBS arrays allocated with sizes ',        &
   ac_num_obs_max,ac_tot_obs_size_max+2048
  CALL umPrint(umMessage,src='u_model_4A')
ELSE
  ac_num_obs_max = 1
  ac_tot_obs_size_max = 1
  obs_flag_len = ac_num_obs_max
  ALLOCATE (obs_flag(obs_flag_len))
  obs_len = ac_tot_obs_size_max
  ALLOCATE (obs(obs_len))
  WRITE(umMessage,*)'U_MODEL_4A - OBS arrays allocated with length 1'
  CALL umPrint(umMessage,src='u_model_4A')
END IF

! ----------------------------------------------------------------------
!  2. Check for nothing-to-do
!

! DEPENDS ON: exitchek
CALL exitchek( internal_model, lexitNOW)

! ----------------------------------------------------------------------
!  3. Start group of timesteps

IF (.NOT. lexitnow) THEN

  WRITE(umMessage,'(A,F10.2,A)') 'Model running with timestep ',          &
                                 secs_per_stepim(atmos_im),' seconds'
  CALL umPrint(umMessage,src='u_model_4A')

#if defined(IBM_XL_FORTRAN)
  ! On IBM turn off buffering of Fortran I/O once initialisation complete
  ! provided the disabling of buffering has been asked for. This is done
  ! by reading  the XLFRTEOPTS environment variable and checking is for the
  ! 'buffering=disable_all' substring (note case and space requirements)
  CALL get_env_var('XLFRTEOPTS', c_xlfrteopts, allow_missing=.TRUE.,      &
                   length=length)

  ! If length > 0, the env var has been read
  IF (length > 0) THEN
    pos = INDEX(c_xlfrteopts, c_disable)

    IF (pos /= 0) THEN    ! Substring found
      IF (printstatus  >=  prstatus_oper .AND. mype == 0) THEN
        WRITE(umMessage,*) 'Disabling Fortran I/O Buffering'
        CALL umPrint(umMessage,src='u_model_4A')
      END IF

      CALL setrteopts(c_disable)
    END IF
  END IF
#endif
END IF
! ----------------------------------------------------------------------
!  3. Start group of timesteps
!
DO WHILE (.NOT. lexitnow) ! Keep looping until the exit flag is set

  ! ----------------------------------------------------------------------
  !  3.1. Start main timestep loop
  !
  !  3.1.1 Increment model time ..

  ! DEPENDS ON: incrtime
  CALL incrtime (                                                   &
      internal_model,icode,cmessage)
  IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)


  ! Keep tabs on PRISM PUT/GET timesteps.
  IF (l_oasis) THEN

    ! Is this a genuine exchange timestep for the components
    ! potentially involved in coupling? 

    IF ( (l_senior) .OR. (l_junior) ) THEN
      ! The location and modulo calculation details for these switches
      ! will need to be reviewed when the full hybrid interface becomes
      ! available.  
      put_step_hyb = (MOD(stepim(atmos_im),oasis_cpl_ts_put_hyb) == 0)
      get_step_hyb = (MOD(stepim(atmos_im),oasis_cpl_ts_get_hyb) == 0)
      put_step_hyb_stats =                                                     &
                     (MOD(stepim(atmos_im),oasis_cpl_ts_put_hyb_stats) == 0)
      get_step_hyb_stats =                                                     &
                     (MOD(stepim(atmos_im),oasis_cpl_ts_get_hyb_stats) == 0)
    END IF ! Senior/Junior atmos-chemistry coupling active

    IF (oasis_couple_ts_aw > 0)                     &
      put_step_aw = (oasis_couple_ts_aw == 1) .OR.  &
        (MOD(stepim(atmos_im),oasis_couple_ts_aw) == 1)

    IF (oasis_couple_ts_wa > 0)                     &
      get_step_wa = (oasis_couple_ts_wa == 1) .OR.  &
       (MOD(stepim(atmos_im),oasis_couple_ts_wa) == 1)

    IF (l_oasis_ocean) THEN 

      IF (oasis_couple_ts_ao > 0)                     &
        put_step_ao = (oasis_couple_ts_ao == 1) .OR.  &
          (MOD(stepim(atmos_im),oasis_couple_ts_ao) == 1)

      IF (oasis_couple_ts_oa > 0)                     &
        get_step_oa = (oasis_couple_ts_oa == 1) .OR.  & 
          (MOD(stepim(atmos_im),oasis_couple_ts_oa) == 1)

      IF (oasis_couple_ts_ao > 0) &
        cpl_update_step=(MOD(stepim(atmos_im),oasis_couple_ts_ao) == 0)

      ! Perform coupling exchanges relating to TS-1
      ! DEPENDS ON: oasis3_geto2a
      CALL oasis3_geto2a()

      ! DEPENDS ON: oasis3_puta2o
      CALL oasis3_puta2o()

    END IF ! Ocean coupling active

    ! Advance date ready for next timestep (if there is one)
    ! DEPENDS ON: OASIS3_ADVANCE_DATE
    CALL oasis3_advance_date()

  END IF ! l_oasis=true

  !  3.1.2 .. set timestep control switches
  ! DEPENDS ON: settsctl
  CALL settsctl (                                                   &
           internal_model,called_from_u_model,meanlev,icode,cmessage)
  IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)

  !  3.1.3 If PPfile initialisation time call PP control routine
  !           for instantaneous data (MEANLEV=0)
  IF (lpp) THEN

    ! DEPENDS ON: ppctl_reinit
    CALL ppctl_reinit(                                              &
       internal_model,icode,cmessage)
    IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)
  END IF

  !  If NetCDF file initialisation time call NetCDF control routine
  !  for instantaneous data
  IF (lnc) THEN

    IF (ltimer) CALL timer('NCFILE_REINIT',5)

    CALL ncfile_reinit()

    IF (ltimer) CALL timer('NCFILE_REINIT',6)
  END IF

  ! Send an EOT message to IOS
  ! He may or may not act on this to purge outstanding items
  CALL io_timestep(stepim(atmos_im))

  ! Close files which are no longer required
  IF (check_close_unneeded_stash_files) THEN
    IF (MOD(stepim(atmos_im),close_unneeded_stash_files_freq)==0) THEN
      CALL close_unneeded_stash_files()
    END IF
  END IF

  !        Integrate atmosphere or ocean by 1 timestep
  IF (internal_model == atmos_im) THEN

    ! River routing
    IF (l_rivers) THEN
      ! Get 'global' atmos horizontal domain sizes from database
      ! in DECOMPDB to set dynamic allocation in ATM_STEP_4A for River routing
      ! on PE 0                                                   .
      g_theta_field_size=                                             &
      decompDB(decomp_standard_atmos)%glsize(1,fld_type_p) *          &
      decompDB(decomp_standard_atmos)%glsize(2,fld_type_p)
      g_river_field_size=                                             &
      decompDB(decomp_standard_atmos)%glsize(1,fld_type_r) *          &
      decompDB(decomp_standard_atmos)%glsize(2,fld_type_r)
    ELSE
      g_theta_field_size=1
      g_river_field_size=1
    END IF

    ! DEPENDS ON: atm_step_4A
    CALL atm_step_4a (                                                      &
           land, cumulus, nbdsc, ntdsc, ntml,                               &
           ccb,cct,                                                         &
           g_theta_field_size,                                              &
           g_river_field_size,                                              &
           obs_flag,obs,obs_flag_len,obs_len,                               &
           ngrgas,grgas_addr)

  END IF        ! internal_model = atmos_im

  IF (l_oasis) THEN

    IF (cpl_update_step) THEN
      !---------------------------------------------------------------------
      ! Ensure atmos coupling data in prognostic areas of D1 are up to date.
      ! The logic here is that cpl_update_step will be TRUE on the timestep
      ! BEFORE coupling is due to take place.

      ! Newly generated coupling data is intercepted after ATM_STEP and
      ! copied to  D1 prognostics.

      ! On the next timestep i.e. a coupling timestep, the prognostic
      ! contents will be sent to the other components in the coupling process.

      ! If there is no subsequent timestep (i.e. if this is the last model
      ! timestep then the D1 contents will be written to the dump
      ! ready for any future restart). There is no "end-of-model"
      ! coupling exchange.
      !---------------------------------------------------------------------

      ! DEPENDS ON: oasis_updatecpl
      CALL oasis_updatecpl(cmessage)

    END IF
  END IF


  !  3.1.4 If dump time, call dump control routine
  IF (ldump) THEN
    IF (ltimer .OR. lstashdumptimer) THEN
      CALL timer('DUMPCTL',5)
    END IF

    ! DEPENDS ON: dumpctl
    CALL dumpctl (meanlev == -1) ! Indicates analysis

    IF (ltimer .OR. lstashdumptimer) THEN
      CALL timer('DUMPCTL ',6)
    END IF
  END IF
  !  3.1.5 If printed output time, call print control routine
  IF (lprint) THEN

    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,*) RoutineName,':Warning, Printing of climate global ' &
                ,'and zonal diagnostics no longer supported'
      CALL umPrint(umMessage,src='u_model_4A')
    END IF  ! PrintStatus test

  END IF

  !  3.1.7 If partial sum/mean creation time, call means control routine
  !        (calls mean PPfield and diagnostic print routines internally)
  IF (lmean) THEN
    IF (ltimer .OR. lstashdumptimer) THEN
      CALL timer('MEANCTL',5)
    END IF

    ! DEPENDS ON: meanctl
    CALL meanctl (submodel,meanlev,icode,cmessage)

    IF (ltimer .OR. lstashdumptimer) THEN
      CALL timer('MEANCTL',6)
    END IF

    IF (icode >  0) THEN
      CALL Ereport(RoutineName,icode,Cmessage)
    END IF
  END IF

  ! 3.1.8 Update history file once dumping and any meaning is complete
  IF (ldump) THEN
    ! checkpoint_dump_im now successfully written. 
    ! Update history file to restart from this
    ! Write history to backup location before overwriting main file
    IF (mype  ==  0) THEN
      CALL temphist(write_temporary,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(umMessage,*)routinename,':Failure writing temporary restart file'
        CALL umPrint(umMessage,src='u_model_4A')
        WRITE(umMessage,*)'Check for problems and restart from main file'
        CALL umPrint(umMessage,src='u_model_4A')
        CALL ereport(routinename,icode,cmessage)
      END IF
      CALL temphist(write_interim,icode,cmessage)
      IF (icode /= 0) THEN
        WRITE(umMessage,*)routinename,':Failure writing main restart file'
        CALL umPrint(umMessage,src='u_model_4A')
        WRITE(umMessage,*)'Check for problems and restart from temporary file'
        CALL umPrint(umMessage,src='u_model_4A')
        WRITE(umMessage,*)'by overwriting main file with temporary file'
        CALL umPrint(umMessage,src='u_model_4A')
        CALL ereport(routinename,icode,cmessage)
      ELSE
        ! Main restart file successfully written, so delete backup
        CALL del_hist()
      END IF

    END IF ! IF (mype == 0)
  END IF ! IF (ldump)

  !  3.1.9 If exit check time, check for immediate exit
  IF (lexit) THEN

    ! DEPENDS ON: exitchek
    CALL exitchek(internal_model, lexitNOW)

    IF (lexitNOW) THEN
      IF (.NOT. ldump) THEN

        WRITE(umMessage,*)routinename,                                         &
           ': Warning: exiting at a period that is not a dump period'
        CALL umPrint(umMessage,src='u_model_4A')
        WRITE(umMessage,*)'Therefore continuing the run will rerun preceding timesteps'
        CALL umPrint(umMessage,src='u_model_4A')
        WRITE(umMessage,*)'This is inefficient and can cause restart problems'
        CALL umPrint(umMessage,src='u_model_4A')

      END IF
      CYCLE ! Restart the loop now (this will exit since lexitnow is true)
    END IF
  END IF
  !  3.1.10 Update ancillary fields if necessary
  IF (lancillary) THEN

    CALL up_ancil (submodel, icode, cmessage)

    IF (icode  /=  0) CALL Ereport(RoutineName,icode,Cmessage)
  END IF
  !  3.1.11 Update boundary fields if necessary
  IF (lboundary) THEN
     ! DEPENDS ON: lbc_update
    CALL lbc_update(submodel, icode, cmessage)
  END IF
  !
  !       End main timestep loop
  ! ----------------------------------------------------------------------

END DO

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    CALL umPrint( '',src='u_model_4A')
    time_end_run = get_wallclock_time()
    WRITE(umMessage,'(A,A,F10.3,A)')                         &
      'u_model_4A: Info: Exiting last atmosphere timestep ',   &
      ' at time=',time_end_run - Start_time,' seconds'
    CALL umPrint(umMessage,src='u_model_4A')
    CALL umPrint( '',src='u_model_4A')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

! ----------------------------------------------------------------------
!  4. Exit processing: Output error messages and perform tidy-up
!

IF (l_oasis) THEN
  ! DEPENDS ON: oasis_tidy
  CALL oasis_tidy()
END IF

!  4.1 Exit processing: If abnormal completion, output error message
iabort = icode
IF (icode /= 0) THEN

  CALL Ereport(RoutineName,icode,Cmessage)

END IF

CALL deallocate_land_soil_arrays()

! ----------------------------------------------------------------------
!  5. Complete Timer call and return
!
icode=iabort
IF (ltimer) CALL timer('U_MODEL_4A ',6)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE u_model_4a
