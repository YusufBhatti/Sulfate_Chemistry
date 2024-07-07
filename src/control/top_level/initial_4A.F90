! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: INITIAL (ENDGAME VERSION) -------------------------------
!
!    Purpose: Initialises the model ready for integration/assimilation.
!   This involves reading the model control files and setting up STASH,
!   reading the initial or restart dump,
!   initialising the ancillary, boundary and interface
!   field control routines and updating the ancillary fields on restart
!   if time to do so, exchanging coupling fields and swapping dumps (if
!   a coupled model), and initialising the assimilation package if
!   necessary.  Subsidiary control routines are called to perform these
!   functions.
!
!    Programming standard: UM Doc Paper 3, version 8 (01/06/2007)
!
!    Logical components covered: C0
!
!    Project task: C0
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE initial_4A(                                            &
     ngrgas,grgas_addr,                                           &
     internal_model,submodel,ngroup,meanlev, dump_unit)

USE atm_fields_mod
USE atm_fields_bounds_Mod
USE atm_d1_indices_mod, ONLY: allocate_d1_indices

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength
USE filenamelength_mod,     ONLY: filenamelength

USE planet_constants_mod, ONLY: recip_epsilon

USE init_casim_run_mod, ONLY: init_casim_run
USE mphys_inputs_mod,   ONLY: l_casim

! IAU scheme:
USE IAU_mod, ONLY:    &
    l_iau,             &
    IAU_FirstCallTS,   &
    IAU_LastCallTS,    &
    L_IAU_DumpTS0State

USE setup_iau_mod, ONLY: setup_iau

USE init_exner_star_mod

! Stochastic physics
USE stph_setup_mod, ONLY: stph_read_ensmem
USE stochastic_physics_run_mod, ONLY: l_rp2, l_skeb2, l_spt,            &
                                      i_pert_theta
USE bl_option_mod, ONLY: off

! Declarations:
USE dump_headers_mod, ONLY: a_inthd, a_realhd, allocate_dump_headers
USE submodel_mod, ONLY: submodel_partition_list, atmos_im
USE d1_array_mod, ONLY: d1
USE stash_array_mod, ONLY: allocate_stash_arrays
USE io
USE um_types
USE eg_alpha_mod
USE eg_alpha_ramp_mod
USE eg_q_to_mix_mod,  ONLY: eg_q_to_mix
USE Control_Max_Sizes
USE UM_ParVars
USE decomp_DB
USE ereport_mod,      ONLY: ereport
USE umPrintMgr
! Make sure um_sleep is picked up from portio.
USE io_dependencies

USE dynamics_testing_mod, ONLY: L_Physics, L_perturb_IC_theta, IntRand_seed
USE coupling_control_mod, ONLY: l_oasis
USE cloud_inputs_mod, ONLY: i_cld_vn, l_pc2_check_init
USE pc2_constants_mod, ONLY: i_cld_pc2
USE mphys_inputs_mod, ONLY:                                           &
     l_mcr_qcf2,          l_mcr_qrain,     l_mcr_qgraup,              &
     l_mcr_qcf2_lbc,      l_mcr_qrain_lbc, l_mcr_qgraup_lbc
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE lbc_read_data_mod, ONLY: albc_num, albc2_starttime_steps, albc_swapstep
USE ukca_option_mod, ONLY: l_ukca

USE nlstcall_mod, ONLY: ncpu, &
                         Num_ALBCs, &
                         ALBC2_StartTime_mins, &
                         ltimer
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun

USE atm_step_local, ONLY: iau_in_initial
USE history, ONLY: h_stepim

USE nlsizes_namelist_mod, ONLY:                                           &
    a_len1_coldepc, a_len1_levdepc, a_len1_rowdepc, a_len2_coldepc,       &
    a_len2_levdepc, a_len2_rowdepc, global_row_length, global_rows,       &
    land_field, model_levels, n_rows, row_length, rows

USE temphist_mod, ONLY: temphist, write_interim

USE model_time_mod, ONLY:                                                    &
    ancillary_stepsim, assim_extrastepsim, assim_firststepim, assim_stepsim, &
    boundary_stepsim, secs_per_stepim, stepim

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam

USE atm_boundary_headers_mod, ONLY: allocate_atm_boundary_headers
USE inancctl_mod, ONLY: inancctl
USE up_ancil_mod, ONLY: up_ancil
USE river_routing_sizes_mod, ONLY: xpa,xua,xva,                             &
                                   ypa,yua,yva

USE land_soil_dimensions_mod, ONLY: allocate_land_soil_arrays, land_index
USE init_radukca_mod, ONLY: init_radukca

USE disturb_veg_category_mod, ONLY: set_disturb_veg_category

!Module code is commented out for now due to variable typing issue caused
!by variables usually passed round by include files
! USE surf_couple_initialise_mod, ONLY: surf_couple_initialise

USE ncfile_init_mod, ONLY: ncfile_init

! ----------------------------------------------------------------------+-------

USE init_emcorr_mod, ONLY: init_emcorr
USE initial_pc2_check_mod, ONLY: initial_pc2_check
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! Subroutine arguments:

INTEGER :: internal_model ! OUT internal model identifier:
!                            !   1:Atmos; 2:Ocean; 3:Slab ; etc
INTEGER :: submodel       ! OUT submodel partition (dump) identifier:
!                            !   1:Atmos; 2:Ocean; etc
INTEGER :: ngroup         ! OUT   - No of steps in "group"n
INTEGER :: meanlev        ! OUT - Mean level indicator
INTEGER, INTENT(IN) :: dump_unit   ! unit attached to input data file 
                                   ! (astart or checkpoint_dump_im)

! 3-D fields of species to be passed down to radiation
INTEGER, INTENT(IN) :: ngrgas
INTEGER, INTENT(OUT) :: grgas_addr(ngrgas)

! Local variables:

INTEGER :: imean      ! Loop index over mean periods
INTEGER :: i, j       ! Loop indices
INTEGER :: ism        ! Loop index over submodels
INTEGER :: ll        ! Counter
INTEGER :: len_filename ! Length of FILENAME array
INTEGER :: Dummyarg  ! Not used, needed to end arg list
INTEGER :: secs_per_stepa ! Atmos timestep length in seconds

CHARACTER(LEN=filenamelength) :: filename
CHARACTER(LEN=2) :: envalue
!
INTEGER :: len_wait_tot     ! Total wait time for boundary data
CHARACTER(LEN=8) :: ch_date2     ! Date from date_and_time
CHARACTER(LEN=10) :: ch_time2    ! Time from date_and_time
INTEGER :: info             ! Return Code from GCOM routine.
INTEGER :: lbc_ntimes       ! No of BC's in communication file
INTEGER :: ms_ntimes        ! No of BC's required in mesoscale
INTEGER :: um_lbc_wait_usec ! No of microseconds to wait
INTEGER :: length           ! Length of ENS_MEMBER's contents

LOGICAL ::                                                        &
   not_enough_lbcs,                                                 &
                    ! True if more LBCs are required
   read_more_values  ! True if we are to read more values in


! Number of dust bins in LBC generation, disabled if not using CreateBC
INTEGER :: ndustbin_in, ndustbin_out


! Error reporting
INTEGER ::    icode       ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='INITIAL_4A')
! Mixing ratio code
LOGICAL, PARAMETER ::                                             &
   l_mr_iau = .FALSE.  ! Use mixing ratio code (if available)

LOGICAL, PARAMETER  ::  called_from_initial = .TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

icode=0
Cmessage=''

!
!  1.3 Option to write RADINCS array to file cache2 removed
!
IF (PrintStatus  >=  PrStatus_Normal) THEN
  WRITE(umMessage,*) 'Fast i/o of radincs directly to core memory'
  CALL umPrint(umMessage,src='initial_4A')
END IF
!
!  Number of CPUs attached to program (ncpu) is hard-wired to 1.
!
ncpu = 1
!
! ---------------------------------------------------------------------
!  2. Initialise STASH control arrays from STASH control file.
! ---------------------------------------------------------------------

! Note that NSECTS=NSECTP, NITEMS=NITEMP : set in WSTLST

CALL allocate_d1_indices()
CALL allocate_stash_arrays()
CALL allocate_dump_headers()
! allocate arrays for the land_soil_dimensions
CALL allocate_land_soil_arrays()


! DEPENDS ON: initctl
CALL initctl(                                                     &
   icode,cmessage )

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to INITCTL'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF


!
! ----------------------------------------------------------------------
!  3. Read appropriate submodel partition dump to memory.  If coupled,
!     page out the D1 part of each dump to its 'swap' file and read the
!     other dump(s) into memory.  Write temporary history file if dumps
!     read successfully on timestep 0.
!
!
!  3.1  Loop over submodel partition dumps
!
submodel=submodel_partition_list(1)

! DEPENDS ON: initdump
CALL initdump(submodel,icode,cmessage, dump_unit)

IF (icode  /=  0) THEN
  WRITE(umMessage,'(A)') 'Failure in call to INITDUMP'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF


! SET_ATM_FIELDS points the fields in atm_fields_mod to the correct parts of D1
! It should be called after INITDUMP (which calls SET_ATM_POINTERS)
! and before INITDIAG (which calls ST_DIAG<n>).

! DEPENDS ON: set_atm_fields
CALL Set_Atm_Fields(d1)

CALL set_disturb_veg_category()

!       Set RUN indicator in atmosphere dump header
!
! DEPENDS ON: set_run_indic_op
CALL set_run_indic_op(                                            &
   icode,cmessage)
IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to SET_RUN_INDIC_OP'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF
!
!  3.3  On NRUN initialise dump LOOKUP headers associated with
!       diagnostic fields with the bare essentials needed to read and
!       write dumps - the rest to be filled in by STASH during the run.
!
IF (h_stepim(atmos_im)  == 0) THEN
  ! DEPENDS ON: inithdrs
  CALL inithdrs(                                                  &
                icode,cmessage)
  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to INITHDRS'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF

END IF    ! End test for NRUN
!

! Read ensemble member from environment variable ENS_MEMBER when
! stochastic physics is active (needed to initialise T+0 lookup hdr)
IF (l_rp2 .OR. l_skeb2 .OR. l_spt .OR. i_pert_theta /= off) THEN
  CALL stph_read_ensmem (length)
  ! Warning if ENS_MEMBER is not used when stochastic physics is on
  IF (length < 1) THEN
    WRITE (cmessage,'(A)') 'Stoch Phys active, but ENS_MEMBER not set or empty'
    CALL Ereport(RoutineName, length, Cmessage)
  END IF
ENDIF

! ----------------------------------------------------------------------
!  6.  Initialise means program control block
!
submodel=submodel_partition_list(1)

! DEPENDS ON: initmean
CALL initmean(submodel,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,'(A)') 'Failure in call to INITMEAN'
               
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

! ----------------------------------------------------------------------
!  4. Set up other control blocks and physical constants
!
!  4.1  Initialise the model time and check that history file data time
!       matches dump(s); set derived time/step information
!
! DEPENDS ON: inittime
CALL inittime(                                                    &
   submodel)

!
!  4.2  Write up temporary history file after successfully reading
!       initial dumps on timestep 0 and setting model_data_time if
!       assimilation run, to allow CRUN from initial dumps.
!
IF (mype  ==  0) THEN
  IF (h_stepim(atmos_im) == 0) THEN
    CALL temphist(write_interim,icode,cmessage)
    IF (icode  /=  0) THEN
      WRITE(umMessage,*) 'Failure in call to TEMPHIST'
      CALL umPrint(umMessage,src='initial_4A')
      CALL Ereport(RoutineName,icode,Cmessage)
    END IF
  END IF
END IF
!
!  4.3  Set up control block for updating ancillary fields
!

CALL inancctl(                                                    &
           icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to INANCCTL'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF
!
!  4.4  Set up control block for updating boundary fields
CALL allocate_atm_boundary_headers()

IF (model_type /= mt_global) THEN

  ALBC_num = 1

  ! Setup for use of two atmos boundary files:
  IF (Num_ALBCs == 2) THEN

    ! Check namelist variable ALBC2_StartTime_mins:
    secs_per_stepa = NINT(SECS_PER_STEPim(atmos_im))

    IF (MOD(ALBC2_StartTime_mins*60, secs_per_stepa) /= 0) THEN
      icode = 1
      WRITE (CMessage,*)                                            &
         'ALBC2_StartTime_mins (', ALBC2_StartTime_mins,             &
         ') does not coincide with the start of a timestep'
      CALL EReport (RoutineName, icode, CMessage)
    ELSE
      ! Convert into steps:
      ALBC2_StartTime_steps = ALBC2_StartTime_mins*60/secs_per_stepa
    END IF

    ! If on a continuation run, we may be able to go straight to
    ! the second boundary file:
    IF (stepim(atmos_im) >= ALBC2_StartTime_steps) ALBC_num = 2

  END IF

END IF  ! not global

IF (model_type /= mt_global) THEN


  ! DEPENDS ON: in_bound
  CALL in_bound(                                                    &
     a_len1_levdepc,a_len2_levdepc,                                 &
     a_len1_rowdepc,a_len2_rowdepc,                                 &
     a_len1_coldepc,a_len2_coldepc,                                 &
                                      ! for dynamic array
                     icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to IN_BOUND'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF

  !
  !    4.4.1  Update atmosphere boundary fields at step zero
  !

  IF (boundary_stepsim(atmos_im)  /=  0) THEN ! If there are BCs

    ! DEPENDS ON: up_bound
    CALL up_bound(submodel,                                         &
                icode,cmessage)

    IF (icode  /=  0) THEN
      WRITE(umMessage,*) 'INITIAL: Failure in call to atmosphere UP_BOUND'
      CALL umPrint(umMessage,src='initial_4A')
      CALL Ereport(RoutineName,icode,Cmessage)
    END IF

    IF (Num_ALBCs == 2) THEN

      ! If the data for the start and end of the current boundary
      ! data interval is to come from different boundary files, we
      ! need to make additional calls to INBOUNDA/UP_BOUND:
      IF (ALBC_num == 1 .AND. stepim(atmos_im) >= ALBC_SwapStep) THEN

        IF (PrintStatus >= PrStatus_Normal) THEN
          CALL umPrint('',src='initial_4A')
          CALL umPrint('INITIAL: Swapping to 2nd atmos boundary file', &
              src='initial_4A')
          CALL umPrint('',src='initial_4A')
        END IF

        ALBC_num = 2


        ! DEPENDS ON: inbounda
        CALL inbounda(                                            &
          a_len1_levdepc,a_len2_levdepc,                          &
          a_len1_rowdepc,a_len2_rowdepc,                          &
          a_len1_coldepc,a_len2_coldepc)

        ! DEPENDS ON: up_bound
        CALL up_bound(submodel,                                   &
          icode,cmessage)


        IF (icode /= 0) THEN
          WRITE(umMessage,*)                                              &
            'INITIAL: Failure in call to atmosphere UP_BOUND'
          CALL umPrint(umMessage,src='initial_4A')
          CALL EReport (RoutineName, icode, cmessage)
        END IF

      END IF ! (ALBC_num == 1 .AND. stepim(atmos_im) >= ALBC_SwapStep)

    END IF ! (Num_ALBCs == 2)

  END IF ! IF (boundary_stepsim(atmos_im)  /=  0)

END IF  ! .NOT. GLOBAL

!
!    4.4.2  Call SETCONA_4A
!
! endgame settings
!
! Set loop index bounds for arrays defined on p, u, v points
!
!      Select Case(model_type)
!         Case(mt_global)
!            i_p_start = 1
!            i_p_end   = row_length
!            j_p_start = 1
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0
!            j_v_end   = n_rows - 1
!         Case(mt_cyclic_lam)
!            i_p_start = 1
!            i_p_end   = row_length
!            j_p_start = 1
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0
!            j_v_end   = n_rows - 1
!         Case(mt_bi_cyclic_lam)
!            i_p_start = 1
!            i_p_end   = row_length
!            j_p_start = 1
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length
!            j_v_start = 0
!            j_v_end   = rows - 1
!         Case(mt_lam)
!            i_p_start = 1
!            i_p_end   = row_length - 1
!            j_p_start = 1
!            j_p_end   = rows
!            i_v_start = 1
!            i_v_end   = row_length - 1
!            j_v_start = 0
!            j_v_end   = n_rows - 1
!         Case Default
!            Print*,'Invalid domain type :',model_type
!            Stop
!         End Select
!
!      i_u_start = 0
!      i_u_end   = row_length-1
!      j_u_start = 1
!      j_u_end   = rows
!      k_p_start = 1
!      k_p_end   = model_levels
!      k_w_start = 0
!      k_w_end   = model_levels
!
! End of endgame settings

! DEPENDS ON: setcona_ctl_4A
CALL setcona_ctl_4a(                                                 &
  icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'INITIAL: Failure in call to SETCONA_4A'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

! 4.5 LBC generation from the UM removed, please use CreateBC

!
!  4.6  Initialise physical constants used in main physics
!       packages - includes radiation lookup tables
!

IF (L_Physics) THEN

  ! DEPENDS ON: initphys
  CALL initphys(icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to INITPHYS'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF

  IF (l_casim) THEN
    CALL init_casim_run()
  END IF

END IF ! on L_Physics

IF (L_emcorr) THEN ! Energy correction enabled
  ! --------------------------------------------------------------------
  !  4.7  Initialise total atmospheric energy & energy correction
  ! --------------------------------------------------------------------
  !   Only done for a new run and then only if the header values in the
  !   dump are set to missing data. If the arrays were set at the
  !   beginning of all NRUNS (new runs) then bit reproducibility
  !   would be lost for short reruns.
  !   The energy correction applied in any day comes from the
  !   change calculated for the previous day. Initialisation for a NRUN
  !   sets the energy correction to zero for the first day.
  !
  !     A_REALHD(18) - total energy flux into the atmosphere
  !                    (used to accumulate change in energy throughout
  !                    a day from physics). Value at start of run zero
  !     A_REALHD(19) - total mass of the atmosphere (wet + dry)
  !                    Note this is not conserved by the dynamics
  !                    The model from UM 5.0 only conserves dry mass.
  !     A_REALHD(20) - total energy of the atmosphere calculated from
  !                    fields in dump.
  !     A_REALHD(21) - energy correction evaluated from previous day of
  !                    run (ie previous run or needs setting to zero).

  IF (stepim(atmos_im) == 0 ) THEN


    CALL init_emcorr(                                              &
       icode,cmessage)


    IF (icode  /=  0) THEN
      WRITE(umMessage,*) 'Failure in call to INIT_EMCORR'
      CALL umPrint(umMessage,src='initial_4A')
      CALL Ereport(RoutineName,icode,Cmessage)
    END IF

  END IF      ! end test on T+0
  !
END IF    !    LEMCORR

!------------------------------------------------------------------------
!  Section 4.8 - Initialise external halos
!                Some T+0 diagnostics and IAU calculations require
!                initialised values of external halos. LBCs will not
!                have been applied yet, so a simple fill of the appropriate
!                variables will suffice to ensure these calculations
!                use sensible values at the boundaries
!------------------------------------------------------------------------
IF (model_type == mt_lam) THEN
  CALL fill_external_halos(u, row_length, rows, model_levels,       &
                           offx, offy)
  CALL fill_external_halos(v, row_length, n_rows, model_levels,     &
                           offx, offy)
  CALL fill_external_halos(theta, row_length, rows, model_levels+1, &
                           offx, offy)
  CALL fill_external_halos(rho, row_length, rows, model_levels,     &
                           offx, offy)
END IF

!------------------------------------------------------------------------
!  Section 4.9 - Initialise external halos
!                In LAM's there is the possibility that the extra moisture
!                fields are used but the LBC files do not hold values for
!                them.
!                This code will initialise the fields to sensible values
!                on the external halos in this case.
!------------------------------------------------------------------------
IF (model_type == mt_lam) THEN
  IF (l_mcr_qrain .AND. .NOT. l_mcr_qrain_lbc)                     &
  CALL set_external_halos(qrain, row_length, rows, model_levels+1, &
                          halo_i, halo_j, 0.0)
  IF (l_mcr_qgraup .AND. .NOT. l_mcr_qgraup_lbc)                   &
  CALL set_external_halos(qgraup, row_length, rows, model_levels+1,&
                          halo_i, halo_j, 0.0)
  IF (l_mcr_qcf2 .AND. .NOT. l_mcr_qcf2_lbc)                       &
  CALL set_external_halos(qcf2 ,row_length, rows, model_levels+1,  &
                          halo_i, halo_j, 0.0)
END IF

! ----------------------------------------------------------------------
!    If coupled model, initialise addresses of coupling fields,
!     and if model has restarted at the end of a coupling period
!     exchange coupling fields and swap data (full ocean model)
!     or both models are at step 0, exchange coupling fields and
!     swap data (in sense O-A at step 0).
!

IF (l_oasis ) THEN
  ! DEPENDS ON: oasis_initialise_2
  CALL oasis_initialise_2 (cmessage)

END IF

! ----------------------------------------------------------------------
!  5. Set timestep group control switches for initial step
!
!  5.1 Set timestep control switches for initial step
!

! DEPENDS ON: settsctl
CALL settsctl (                                                   &
             internal_model,called_from_initial,meanlev,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to SETTSCTL'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

!
!  5.2 Initialise PP files at step 0
!

! DEPENDS ON: ppctl_init
CALL ppctl_init(                            &
   submodel,icode,cmessage)

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to PPCTL'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF
!
!  5.3  Initialise assimilation package (not if assimilation completed)
!
IF ( l_iau .AND. (assim_stepsim(atmos_im)+assim_extrastepsim(atmos_im) >  0 ) &
  .OR.                                                               &
   l_iau .AND. (   stepim(atmos_im)  <   assim_firststepim(atmos_im) + &
                   assim_stepsim(atmos_im) + assim_extrastepsim(atmos_im)) &
   ) THEN
  ! DEPENDS ON: in_acctl
  CALL in_acctl(                                                    &
                    icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to IN_ACCTL'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF
END IF
!
! 5.4 Initialise NetCDF files at step 0
!
IF (ltimer) CALL timer('NCFILE_INIT',5)

CALL ncfile_init()

IF (ltimer) CALL timer('NCFILE_INIT',6)

!----------------------------------------------------------------------
! [6.1]: Incremental Analysis Update (IAU).
!----------------------------------------------------------------------

IF (l_iau) CALL Setup_IAU

IF (l_iau .AND. stepim(atmos_im) >= IAU_FirstCallTS .AND. &
                stepim(atmos_im) <= IAU_LastCallTS) THEN

  IF (ltimer) CALL timer('IAU',3)

  iau_in_initial = .TRUE.

  ! DEPENDS ON: iau
  CALL iau (                                      &
             l_mr_iau,                            & ! in
             frac_typ,                            & ! in
             snowdepth,          nsnow,           & ! in
             vol_smc_wilt,       vol_smc_crit,    & ! in
             vol_smc_sat,                         & ! in
             u,      v,          w,               & ! inout
             theta,  rho,        murk,            & ! inout
             q,      qCL,        qCF,             & ! inout
             TStar,  TStar_tile, Deep_soil_temp,  & ! inout
             smcl,   tsnowlayer, tstar_anom,      & ! inout
             Pstar,  p,                           & ! inout
             p_theta_levels,                      & ! inout
             exner_rho_levels,                    & ! inout
             exner_theta_levels,                  & ! inout
             snodep,                              & ! inout
             cf_area,                             & ! inout
             cf_bulk,                             & ! inout
             cf_liquid,                           & ! inout
             cf_frozen,                           & ! inout
             dust_div1, dust_div2, dust_div3,     & ! inout
             dust_div4, dust_div5, dust_div6,     & ! inout
             ozone_tracer, SO2, tracer_ukca )       ! inout

  CALL init_exner_star(exner_rho_levels,thetav,                   &
                       m_v, m_cl, m_cf, m_r, m_gr, m_cf2,         &
                       exner_surf)

  IF (ltimer) CALL timer('IAU',4)

  ! If required, write out TS0 dump, including any IAU increments
  !  which may have just been added:
  IF (stepim(atmos_im) == 0 .AND. L_IAU_DumpTS0State) THEN

    ! DEPENDS ON: dumpctl
    CALL dumpctl (                                                &
        meanlev == -1) ! Indicates analysis

  END IF

ELSE IF (i_cld_vn == i_cld_pc2 .AND.                                    &
        (l_pc2_check_init  .AND. .NOT. l_nrun_as_crun) .AND.            &
        stepim(atmos_im) == 0) THEN

  IF (.NOT. l_nrun_as_crun) THEN

    CALL initial_pc2_check(exner_theta_levels, p_theta_levels, theta,   &
                           cf_bulk, cf_liquid, cf_frozen, q, qcl, qcf, l_mr_iau)

! Update ENDGame prognostics which may have changed
    CALL eg_q_to_mix(tdims_l, tdims_s, q, qcl, qcf, qcf2, qrain, qgraup,&
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   m_v, m_cl, m_cf, m_cf2, m_r, m_gr, swap_in = .FALSE.)

    thetav = theta*(1.0+m_v*recip_epsilon)

  ELSE
!   l_nrun_as_crun so we don't apply PC2 checks
    cmessage = 'l_nrun_as_crun is TRUE' //  &
               'Ignoring value of l_pc2_check_init (treating as false)'
    icode = -1
    CALL ereport(RoutineName, icode, cmessage)
  END IF
END IF

!
! By this stage, all of the prognostic atmospheric variables have been set to
! their final values. Add bit-level perturbations to the initial theta
! (and thetav(d)) value if requested. This cannot be done earlier in
! case the optional call to dumpctl above packs the dumps to lower
! precision and wipes out these perturbations.
!
IF (L_perturb_IC_theta .AND. h_stepim(atmos_im) == 0) THEN
  ! DEPENDS ON: perturb_theta_ctl
  CALL perturb_theta_ctl(                                          &
     row_length, rows, model_levels,          &
     global_row_length, global_rows,          &
     at_extremity, offx, offy, IntRand_Seed,  &
     datastart                                &
     )
END IF

! ----------------------------------------------------------------------
!
!  7.1  Get derived diagnostics from start fields (atmosphere)
!
IF (stepim(atmos_im) == 0) THEN

  ! DEPENDS ON: initdiag
  CALL initdiag(Dummyarg)


  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to INITDIAG'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF
END IF
!
! 7.2 Code to update the boundaries now moved forwards to 4.4.1
!
!
!  7.3 Update ancillary fields in dump if start time corresponds to
!      an ancillary field update time. Also done at T+0 with values
!      updated to half a period back from first standard update time
!      to ensure reproducibility between long runs and new runs
!      started from dump at any time.
!

IF (ancillary_stepsim(atmos_im) >  0) THEN
  IF (stepim(atmos_im) == 0 .OR.                                  &
     MOD(stepim(atmos_im),ancillary_stepsim(atmos_im)) == 0) THEN 
    CALL up_ancil (submodel, icode, cmessage)
  END IF
END IF

IF (icode  /=  0) THEN
  WRITE(umMessage,*) 'Failure in call to UP_ANCIL'
  CALL umPrint(umMessage,src='initial_4A')
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

!  7.3.3 Ensure that convective cloud cover and liquid water path
!        are consistent with convective cloud base & top. (Corrects
!        for occasional problems caused by reconfiguration.)
IF (stepim(atmos_im) == 0) THEN

  ! DEPENDS ON: init_cnv
  CALL init_cnv(                                                 &
     icode,cmessage)
  IF (icode  /=  0) THEN
    WRITE(umMessage,*) 'Failure in call to INIT_CNV'
    CALL umPrint(umMessage,src='initial_4A')
    CALL Ereport(RoutineName,icode,Cmessage)
  END IF

END IF
!

!  7.3.4 ! initialization of radiative feedback
!   Initialise and check address array to feed chemical tracers from
!   UKCA into radiation scheme.  This needs to be done even for a CRUN
!   so no need to check for stepim(atmos_im) == 0.
grgas_addr = -1
IF (l_ukca) CALL init_radukca(ngrgas,grgas_addr)

! 7.4 LBC generation from the UM removed, please use CreateBC

!------------------------------------------------------------------------------
!Reordering of subroutine calls to put all calls to land surface together
!Dependence on include files removed to keep them out of the JULES repo
!------------------------------------------------------------------------------
!  7.3.1 Initialize tiled prognostics, gridbox mean vegetation
!        parameters and TRIFFID accumulation prognostics.
! DEPENDS ON: surf_couple_initialise
CALL surf_couple_initialise(stepim(atmos_im),                                 &
              a_inthd,               &
              land_index,            &
              !Arguments for river routing
              a_realhd,              &
              xpa, xua, xva, ypa, yua, yva)


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE initial_4A
