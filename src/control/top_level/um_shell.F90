! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: UM_SHELL -------------------------------------------------
!
!  Purpose: Control routine for the Atm Model.
!           Acquires size information needed for dynamic allocation of
!           configuration-dependent arrays and calls U_MODEL (the
!           master control routine) to allocate the arrays and perform
!           the top-level control functions and timestepping.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

SUBROUTINE um_shell

USE UM_Config, ONLY:                                                        &
    appInit,                                                                 &
    appTerminate,                                                            &
    exe_UM
USE mpl, ONLY:                                                              &
    mpl_max_processor_name,                                                  &
    mpl_thread_multiple,                                                     &
    mpl_thread_serialized,                                                   &
    mpl_thread_funneled,                                                     &
    mpl_thread_single
USE affinity_mod, ONLY:                                                      &
    affinity_write
USE wait_policy_mod, ONLY: init_wait_policy
USE windmax_mod, ONLY: eg_unit, windmax_file

USE run_info, ONLY: set_start_time,start_time
USE timestep_mod, ONLY: timestep_number

USE oasis_atm_data_mod, ONLY: comm_in

#if defined(MCT)
USE OASIS3_split_comm_mod, ONLY: oasis3_split_comm
#endif

!$ USE omp_lib           ! Note OpenMP sentinel

USE filenamelength_mod, ONLY:                                               &
    filenamelength

USE nlcfiles_namelist_mod, ONLY: &
    stashreq_log_file => streqlog

USE file_manager, ONLY: &
    assign_file_unit, release_file_unit, um_file_type, init_file_loop

USE atm_fields_bounds_Mod
USE d1_array_mod, ONLY: allocate_d1_array
USE Halo_Exchange, ONLY:                                                    &
    Halo_Exchange_Init
USE Non_Blocking_Halo_Exchange, ONLY:                                       &
    Non_Blocking_Halo_Exchange_Init
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE drhook_control_mod, ONLY: drhook_control_enable, drhook_control_disable
USE io,                 ONLY:                                               &
    ioInit,                                                                  &
    file_close,                                                              &
    is_Unit_Open
USE ios
USE IOS_Constants,      ONLY:                                               &
    IOS_maxServers
USE IOS_Init,           ONLY:                                               &
    IOS_Setup,                                                               &
    IOS_Run
USE IOS_Stash_Common,   ONLY:                                               &
    isUsingAsyncStash,                                                       &
    isUSingAsyncDumps
USE IOS_Stash, ONLY:                                                        &
    IOS_stash_client_fini
USE IOS_Model_Geometry, ONLY:                                               &
    IOS_Client_Geometry_Init
USE mppio_file_utils, ONLY:                                                 &
    mppio_file_utils_init
USE model_file,     ONLY:                                                   &
    model_file_close
USE mpp_conf_mod,       ONLY:                                               &
    extended_halo_size_ew,                                                   &
    extended_halo_size_ns
USE ereport_mod, ONLY: ereport,ereport_finalise
USE UM_ParVars, ONLY: change_decomposition, g_datastart_f, g_blsize
USE UM_ParCore, ONLY: mype, nproc_max, nproc_um_npes=>nproc
USE Decomp_DB
USE umPrintMgr

USE gcom_mod, ONLY: gc_alltoall_multi, gc_alltoall_version
USE decomp_params, ONLY:                                                    &
    decomp_standard_atmos,                                                   &
    decomp_unset
USE rimtypes
USE lbc_mod
USE coupling_control_mod, ONLY: l_oasis
USE version_mod, ONLY:                                                      &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,                    &
    nlevp, npslevp, npslistp

USE um_version_mod, ONLY: um_version_char

USE nlsizes_namelist_mod, ONLY:                             &
    global_row_length, global_rows, max_intf_model_levels,  &
    max_lbcrow_length, max_lbcrows, model_levels, n_intf_a, &
    river_row_length, river_rows, row_length, rows

USE errormessagelength_mod, ONLY: errormessagelength
USE get_env_var_mod,   ONLY: get_env_var
USE model_domain_mod,  ONLY: model_type
USE io_configuration_mod, ONLY: print_runtime_info
USE turb_diff_mod, ONLY: l_print_pe

USE umnetcdf_mod, ONLY: nc_file_close, nc_is_file_open

USE um_submodel_init_mod, ONLY: um_submodel_init

USE get_wallclock_time_mod, ONLY: get_wallclock_time
 
USE timer_mod, ONLY: timer

USE nlstcall_mod, ONLY: lstashdumptimer
 
IMPLICIT NONE
 
!
!  Local parameters
!
CHARACTER(LEN=*) :: RoutineName
PARAMETER (RoutineName = 'UM_SHELL')
!
!  Local variables
!
INTEGER :: icode       ! Work - Internal return code
INTEGER :: istatus     ! RETURN STATUS FROM OPEN
INTEGER :: loop_pe
INTEGER :: loop_pe_start
INTEGER :: loop_pe_stop

CHARACTER(LEN=filenamelength) :: shared_filename    = "dummy filename" !Namelist
CHARACTER(LEN=filenamelength) :: atmoscntl_filename = "dummy filename" !Namelist
CHARACTER(LEN=errormessagelength) :: cmessage ! Work - Internal error message
CHARACTER(LEN=errormessagelength) :: iomessage
INTEGER :: atm_nprocx          ! number of procs EW for atmosphere
INTEGER :: atm_nprocy          ! number of procs NS for atmosphere
INTEGER :: length              ! length of env var contents
INTEGER :: err                 ! error return from subroutine calls
INTEGER :: stash_unit          ! Unit to use for STASH requests log file
INTEGER :: atmoscntl_unit      ! Unit to use for ATMOSCNTL file
INTEGER :: shared_unit         ! Unit to use for SHARED file

CHARACTER(LEN=10) :: c_thread          ! to get nproc_x and nproc_y from
CHARACTER(LEN=8) :: c_nproc            ! to get nproc_x and nproc_y from
CHARACTER(LEN=10) :: c_coupl            ! to get whether coupled run from env.

TYPE(um_file_type), POINTER :: um_file, um_file_next


INTEGER :: um_nam_max_seconds
!
CHARACTER(LEN=8) :: c_nam_max_seconds
!
!  Localized sizes for ocean decomposition:
INTEGER ::                                                                   &
    row_length_oce , rows_oce

INTEGER :: sect_err, rnl_err, um_rnl_skip

CHARACTER(LEN=8) :: c_um_rnl_skip
CHARACTER(LEN=8) :: ch_date2   !  Date returned from date_and_time
CHARACTER(LEN=10) :: ch_time2  !  Time returned from date_and_time

! Variables needed to close all the files
INTEGER                      :: i

! Variables for IO Server setup
LOGICAL                      :: isIOServer
INTEGER                      :: numIOServers
INTEGER                      :: errorCode
CHARACTER(LEN=32)            :: c_io_pes
CHARACTER(LEN=10)            :: thread_level_setc
INTEGER                      :: thread_level_set

INTEGER :: dummy_comm     ! Dummy communicator for OASIS

CHARACTER(LEN=mpl_max_processor_name) :: env_myhost
INTEGER                               :: env_myhost_len

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
REAL(KIND=jprb)               :: zhook_rendez_vous
LOGICAL(KIND=jpim)            :: luser_comm

INTEGER :: dump_unit      ! unit attached to input data file

REAL :: time_end_run

! Prevent any calls to DrHook before MPI has been initialised.
CALL drhook_control_disable()

!-----------------------------------------------------------------------

cmessage = ' '
numIOServers = 0  ! Initialise IO server variable

! ----------------------------------------------------------------------
!----------------------------------------------------------------------
! 1.0 Initialise Message Passing Libraries
!

! Get the atmosphere decomposition
! Cannot call get_env_var before GCOM initialisation.
CALL GET_ENVIRONMENT_VARIABLE('UM_THREAD_LEVEL',c_thread,length,err)
IF (err  /=  0 .OR. length == 0) THEN
  CALL umPrint('Warning: Environment variable UM_THREAD_LEVEL has ' //       &
      'not been set.',src='um_shell')
  CALL umPrint('Setting thread_level to multiple',src='um_shell')
  thread_level_setc = 'MULTIPLE'
ELSE
  READ(c_thread,'(A10)') thread_level_setc
END IF

SELECT CASE (thread_level_setc)
CASE ('MULTIPLE')
  thread_level_set = mpl_thread_multiple
CASE ('SERIALIZED')
  thread_level_set = mpl_thread_serialized
CASE ('FUNNELED')
  thread_level_set = mpl_thread_funneled
CASE ('SINGLE')
  thread_level_set = mpl_thread_single
CASE DEFAULT
  WRITE(umMessage,'(A,A,A)') 'Warning: Thread level ', thread_level_setc,    &
      ' not recognised, setting to MULTIPLE.'
  CALL umPrint(umMessage,src='um_shell')
  thread_level_set = mpl_thread_multiple
END SELECT

! The total number of processors required (nproc_max) is determined
!  by a call to gc_init/gc_init_thread:

! Cannot call get_env_var before GCOM initialisation.
CALL GET_ENVIRONMENT_VARIABLE('COUPLER',c_coupl,STATUS=err)
IF (err  /=  0 .OR. c_coupl == "none" .OR. c_coupl == "") THEN
  l_oasis = .FALSE.
ELSE
  l_oasis = .TRUE.
END IF

IF ( l_oasis ) THEN

  !   Call routine to initialise OASIS

  ! DEPENDS ON: oasis_initialise
  CALL oasis_initialise(mype,nproc_max,comm_in,thread_level_set)

  !   Permit calls to DrHook, then call the first top-level caliper.
  CALL drhook_control_enable()

  luser_comm = .TRUE.
  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle,luser_comm, &
                          INT(comm_in,KIND=jpim))

ELSE ! l_oasis

  !   Standard UM GCOM initialisation when OASIS is not used

  CALL gc_init_thread(mype,nproc_max, thread_level_set)

  !   Permit calls to DrHook, then call the first top-level caliper.
  CALL drhook_control_enable()

  IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

END IF ! l_oasis

! initialise timers after MPI initialisation
CALL set_start_time()

CALL appInit(exe_UM)

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    CALL umPrint( '',src='um_shell')
    CALL umPrint('um_shell: Info: Start time set ',src='um_shell')
    CALL umPrint( '',src='um_shell')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

! Determine the number of processors the model has been configured to run on:
CALL get_env_var('UM_ATM_NPROCX',c_nproc,allow_missing=.TRUE.,length=length)
IF (length <  0) THEN
  cmessage = 'Warning: Environment variable UM_ATM_NPROCX has '              &
      //'not been set. Setting nproc_x to 1.'
  CALL ereport(RoutineName, length, cmessage)
  atm_nprocx=1
ELSE
  READ(c_nproc,'(I4)') atm_nprocx
END IF

CALL get_env_var('UM_ATM_NPROCY',c_nproc,allow_missing=.TRUE.,length=length)
IF (length < 0) THEN
  cmessage = 'Warning: Environment variable UM_ATM_NPROCY has '              &
      //'not been set. Setting nproc_y to 1.'
  CALL ereport(RoutineName, length, cmessage)
  atm_nprocy=1
ELSE
  READ(c_nproc,'(I4)') atm_nprocy
END IF

! Calculate total number of atmos processors required:
nproc_um_npes = atm_nprocx * atm_nprocy


! Get number of I/O PEs from the environment
CALL get_env_var('FLUME_IOS_NPROC',c_io_pes,allow_missing=.TRUE.,length=length)
IF (length < 0) THEN
      ! If not specified, try to work out a valid number
  numIOServers=nproc_max-nproc_um_npes
  icode=-10
  WRITE(cmessage,'(A,A,I4)')                                                 &
      'FLUME_IOS_NPROC environment variable not set',                        &
      ', I/O PE count set to ',numIOServers
  CALL ereport(routinename,icode,cmessage)
ELSE
  READ (c_io_pes,'(I5)') numIOServers
END IF

IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
  icode=-10
  WRITE(cmessage,'(A,I4)')                                                   &
      'I/O PE count is outside allowed range: ',numIOServers
  CALL ereport(routinename,icode,cmessage)

      ! try to work out a valid number
  numIOServers=nproc_max-nproc_um_npes
  IF ( numIOServers < 0 .OR. numIOServers > IOS_maxServers ) THEN
    icode=10
    WRITE(cmessage,'(A)')                                                    &
        'A valid I/O PE count could not be set'
    CALL ereport(routinename,icode,cmessage)
  END IF
ELSE
  CALL umPrint('Enabling '//TRIM(str(numIOServers))//' I/O PEs',             &
      pe=0,src='um_shell')
END IF

! Check the breakdown of processors requested by environment variables
! (UM_ATM_NPROCX * UM_ATM_NPROCY + FLUME_IOS_NPROC) matches the total
! number obtained by GCOM
IF (nproc_max /= nproc_um_npes + numIOServers) THEN
  icode = 100
  WRITE(cmessage,'(A,i7,A,i7,A)')                                            &
      'UM started on ', nproc_max,                                           &
      ' PEs but ',                                                           &
      nproc_um_npes+numIOServers,                                            &
      ' asked for. Please adjust decomposition'
  CALL ereport(routinename,icode,cmessage)
END IF

IF (nproc_max  <   0) THEN
  CALL umPrint( 'Parallel initialisation failed',src='um_shell')
  GO TO 9999
ELSE
   
  ! Set GCOM to use the alternative version of RALLTOALLE
  ! throughout the run
  CALL gc_setopt(gc_alltoall_version, gc_alltoall_multi, err)

  CALL umPrint('',src='um_shell')
  CALL umPrint('**************************** PROCESSOR '//                   &
      'INFORMATION ****************************',src='um_shell')
  CALL umPrint('',src='um_shell')
  CALL umPrint(TRIM(str(nproc_max))//' Processors initialised.',             &
      src='um_shell')
  CALL MPL_Get_processor_name(env_myhost, env_myhost_len, err)
  IF (err /= 0) THEN
    CALL umPrint('I am PE '//TRIM(str(mype)),src='um_shell')
  ELSE
    CALL umPrint('I am PE '//TRIM(str(mype))//' on '//TRIM(env_myhost),      &
        src='um_shell')
  END IF
  ! Only want OpenMP section executing if OpenMP is compiled in,
  ! so protect by sentinal
!$OMP PARALLEL DEFAULT(NONE)
!$OMP MASTER
!$  WRITE(umMessage,'(A,I2,A)') 'I am running with ',                          &
!$      omp_get_num_threads(),' thread(s).'
!$  CALL umPrint(umMessage,src='um_shell')
!$  WRITE(umMessage,'(A,I6)') 'OpenMP Specification: ',openmp_version
!$  CALL umPrint(umMessage,src='um_shell')
!$OMP END MASTER
!$OMP END PARALLEL
END IF

#if defined(IBM_XL_FORTRAN)
! On IBM force buffering of Fortran I/O for initialisation
CALL setrteopts('buffering=enable')
#endif
!
CALL timer(RoutineName,1)

!----------------------------------------------------------------------
!
!    Open the two main namelist files on PE 0 only. 
!    All runtime control variables are subsequently read in from here.
!
CALL get_env_var("ATMOSCNTL",atmoscntl_filename)
CALL assign_file_unit(atmoscntl_filename, atmoscntl_unit, handler="fortran", &
                      id="atmoscntl")

IF (mype == 0) THEN
  OPEN(UNIT=atmoscntl_unit,FILE=atmoscntl_filename, ACTION='READ',           &
       IOSTAT=istatus, IOMSG=iomessage)

  IF (istatus /= 0) THEN
    icode=500
    CALL umPrint( ' ERROR OPENING ATMOSPHERE-ONLY NAMELIST FILE',src='um_shell')
    WRITE(umMessage,'(A,A)') ' FILENAME =',TRIM(atmoscntl_filename)
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A,I0)') ' IOSTAT =',istatus
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A,A)') ' IOMSG =',TRIM(iomessage)
    CALL umPrint(umMessage,src='um_shell')
    GO TO 9999
  END IF
END IF

CALL get_env_var("SHARED_NLIST",shared_filename)
CALL assign_file_unit(shared_filename, shared_unit, handler="fortran", &
                      id="shared")

IF (mype == 0) THEN
  OPEN(UNIT=shared_unit,FILE=shared_filename, ACTION='READ', IOSTAT=istatus,  &
       IOMSG=iomessage)

  IF (istatus /= 0) THEN
    icode=510
    CALL umPrint( ' ERROR OPENING SHARED NAMELIST FILE',src='um_shell')
    WRITE(umMessage,'(A,A)') ' FILENAME =',TRIM(shared_filename)
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A,I0)') ' IOSTAT =',istatus
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A,A)') ' IOMSG =',TRIM(iomessage)
    CALL umPrint(umMessage,src='um_shell')
    GO TO 9999
  END IF
END IF
! ------------------------------------------------------------------
!  0.1 Get submodel/internal model components of model run.
!
icode=0
CALL UM_Submodel_Init(icode)
IF (icode  /=  0) THEN
  cmessage = 'Error calling UM_Submodel_init'
  GO TO 9999
END IF

CALL um_shell_banner('Start')

! Initialise I/O subsystem
CALL ioInit()

! Start I/O Server if required
      ! On exit the IO slave PEs have finished the whole job so can
      ! go to the end of the routine. FIXME check cpu counts against deco
isIOServer=IOS_Setup(numIOServers)

! Set mype to the local model rank
mype=model_rank

! Initialise the OMP wait policy 
CALL init_wait_policy()

!Report affinity
CALL affinity_write(global_comm, writer_rank=0,     &
             long_form=(PrintStatus >= PrStatus_Diag))

#if defined(MCT)
! If we're employing OASIS3-MCT then we have to tell OASIS
! which processes are IO ones and which are atmos ones.
CALL oasis3_split_comm(.NOT. isIOServer)

! Certain collective calls which need performing on all PEs
! regardless of whether they're actually involved in
! coupling or not.
! DEPENDS ON: OASIS3_grid
CALL oasis3_grid(.FALSE.)
#endif

IF (isIOServer) THEN
  CALL IOS_Run()
ELSE
  CALL umPrint('Running Atmospheric code as pe '//TRIM(str(mype)),           &
      src='um_shell')

  !----------------------------------------------------------------------

  CALL mppio_file_utils_init()

  ! ----------------------------------------------------------------------
  !  Read history files for NRUN or CRUN.
  ! DEPENDS ON: readhist
  CALL readhist ( icode,cmessage )
  IF (icode > 0) GO TO 9999

  !  Read Control file on standard input.
  !
  ! DEPENDS ON: readcntl
  CALL readcntl ( icode,cmessage )
  IF (icode  >   0) GO TO 9999

  !----------------------------------------------------------------------
  ! Operational settings:

  ! Allow Override of namelist input FASTRUN in operational environment.
  ! DEPENDS ON: set_fastrun
  CALL set_fastrun

  !----------------------------------------------------------------------
  !  Call READLSTA to read namelists to control atmosphere integration
  !  and diagnostic point print.
  ! DEPENDS ON: readlsta
  CALL readlsta()

  !----------------------------------------------------------------------
  !  1.1 Get configuration-dependent sizes needed for dynamic allocation.
  !
  ! DEPENDS ON: readsize
  CALL readsize(dump_unit)

  ! Decompose atmosphere data and find new local data size

  CALL decompose(decomp_standard_atmos,                                      &
      global_row_length,global_rows,model_levels,                            &
      river_rows, river_row_length,                                          &
      model_type,                                                            &
      atm_nprocx, atm_nprocy,                                                &
      extended_halo_size_ew,                                                 &
      extended_halo_size_ns,                                                 &
      rimwidtha, nrima_max,row_length,rows )

  ! Set up the atmosphere decomposition in PARVARS
  CALL change_decomposition(decomp_standard_atmos,icode)

  ! Now we have a decomposition initialise the blocking and
  ! non-blocking halo swap modules
  CALL Halo_Exchange_Init()
  CALL Non_blocking_Halo_Exchange_Init()

  IF (icode  /=  0) GO TO 9999

  ! Output range of gridpoints for each PE
  IF (PrintStatus >= PrStatus_Diag) THEN
    loop_pe_start=mype
    loop_pe_stop =mype
    IF (mype  ==  0) THEN
      loop_pe_start=0
      loop_pe_stop=atm_nprocx*atm_nprocy-1
    END IF

    CALL umPrint('',src='um_shell')
    CALL umPrint('Range of gridpoints for processing elements:',             &
        src='um_shell')
    CALL umPrint('',src='um_shell')
    WRITE(umMessage,'(A7,A2,A15,A2,A15)')                                    &
        '     PE',' |','  West -   East',' |',' South -  North'
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A7,A2,A15,A2,A15)')                                    &
        '-------','-+','---------------','-+','---------------'
    CALL umPrint(umMessage,src='um_shell')
    DO loop_pe = loop_pe_start,loop_pe_stop
      WRITE(umMessage,'(I7,A2,I6,A3,I6,A2,I6,A3,I6)')                        &
          loop_pe,' |',                                                      &
          g_datastart_f(1,1,loop_pe),' - ',                                  &
          g_datastart_f(1,1,loop_pe)+g_blsize(1,1,loop_pe)-1,' |',           &
          g_datastart_f(2,1,loop_pe),' - ',                                  &
          g_datastart_f(2,1,loop_pe)+g_blsize(2,1,loop_pe)-1
      CALL umPrint(umMessage,src='um_shell')
    END DO
    WRITE(umMessage,'(A7,A2,A15,A2,A15)')                                    &
        '-------','-+','---------------','-+','---------------'
    CALL umPrint(umMessage,src='um_shell')

  END IF

  ! Call DERVSIZE (the call in READSIZE has been deleted)

  icode=0
  ! DEPENDS ON: dervsize
  CALL dervsize(icode,cmessage)
  IF (icode  /=  0) GO TO 9999

  !     Ensure that domain decomposition is set for Atmosphere
  CALL change_decomposition (decomp_standard_atmos,icode)
  IF (icode /= 0) THEN
    WRITE(umMessage,'(A,A)') ' Error returned in CHANGE_DECOMPOSITION',      &
        ' before DERV_LAND_FIELD.'
    CALL umPrint(umMessage,src='um_shell')
    WRITE(umMessage,'(A,I5)') ' Error code ',icode
    CALL umPrint(umMessage,src='um_shell')
    WRITE(cmessage,'(A)') 'UM_SHELL : Error in CHANGE_DECOMPOSITION.'
    GO TO 9999   !  Exit
  END IF

  !     For MPP jobs, calculate the no of land-points on each PE.
  ! DEPENDS ON: derv_land_field
  CALL derv_land_field (dump_unit, icode, cmessage)
  IF (icode >  0) THEN
    CALL umPrint( 'Error returned from DERV_LAND_FIELD.',src='um_shell')
    WRITE(umMessage,'(A,I5)') 'Error code ',icode
    CALL umPrint(umMessage,src='um_shell')
    GO TO 9999   !  Exit
  END IF

  IF (icode >  0) GO TO 9999

  n_intf_a = 1
  max_intf_model_levels = 1
  max_lbcrow_length = 1
  max_lbcrows = 1

  CALL umPrint('********************************************'//            &
      '***********************************',src='um_shell')
  !-----------------------------------------------------------------------
  ! 1.2 Call STASH_PROC: top level control routine for processing of
  !                      STASH requests and STASH addressing.

  ! Open the stash requests file
  IF (mype == 0) THEN
    CALL assign_file_unit(stashreq_log_file, stash_unit, handler="fortran", &
                          id="stash_req_log")
    OPEN(UNIT=stash_unit, FILE=stashreq_log_file, STATUS='REPLACE')
  END IF

  IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
    CALL assign_file_unit(windmax_file, eg_unit, handler="fortran")
    OPEN(UNIT=eg_unit,FILE=windmax_file)
  END IF


  ! DEPENDS ON: stash_proc
  CALL stash_proc(icode, cmessage  )
  IF (icode >  0) GO TO 9999

  ! ----------------------------------------------------------------------
  !  1.3 Calculate addresses of super arrays passed down for dynamic
  !      allocation.
  !
  CALL allocate_d1_array()
  icode=0
  ! DEPENDS ON: um_index
  CALL um_index(icode,cmessage)

  IF (icode >  0) GO TO 9999

  ! ----------------------------------------------------------------------
  !  1.5 Set up geometry of the model on the IO servers if needed
  !
  !  * Note there is a global barrier (INCLUDING THE IO SERVERS) here *
  IF (isUsingAsyncStash() .OR. isUsingAsyncDumps()) THEN
        ! We need to tell the submodels about our geometry
    CALL IOS_Client_Geometry_Init(decomp_standard_atmos)
  END IF

  ! ----------------------------------------------------------------------
  !  2. Call U_MODEL_4A master routine to allocate the main data arrays
  !     and do the calculations.
  !

  ! DEPENDS ON: u_model_4A
  CALL u_model_4A (dump_unit)

  ! Make sure all the Files are properly Closed. We call model_file_close not
  ! file_close to ensure that any open files which are PP and hence have
  ! cached lookups, have the lookups commited to disk before the close.
  um_file => init_file_loop(handler="portio")
  DO WHILE (ASSOCIATED(um_file))
    ! Store the next file so we can still access it once we've freed this one!
    um_file_next => um_file % next
    IF (um_file % pp_meta % managed) THEN
      WRITE(umMessage,'(A,I3,A)')'Managed unit ', &
          um_file % UNIT,' is open at end of run, closing'
      CALL umPrint(umMessage,src='um_shell')
      CALL model_file_close(um_file % UNIT, um_file % filename)
      CALL release_file_unit(um_file % UNIT, handler="portio")
    ELSE IF (is_unit_open(um_file % UNIT)) THEN
      WRITE(umMessage,'(A,I3,A)')'Unmanaged unit ', &
          um_file % UNIT,' is open at end of run, closing'
      CALL umPrint(umMessage,src='um_shell')
      CALL file_close(um_file % UNIT, um_file % filename)
      CALL release_file_unit(um_file % UNIT, handler="portio")
    ELSE
      WRITE(umMessage,'(A,I3,A)')'Unit ',um_file % UNIT,' is closed'
      CALL umPrint(umMessage,src='um_shell')
    END IF
    um_file => um_file_next
  END DO

  IF (mype == 0) THEN
    um_file => init_file_loop(handler="netcdf")
    DO WHILE (ASSOCIATED(um_file))
      ! Store the next file so we can still access it once we've freed this one!
      um_file_next => um_file % next
      IF (nc_is_file_open(um_file)) THEN
        WRITE(umMessage,'(A,I3,A)')'NetCDF unit ', &
            um_file % UNIT,' is open at end of run, closing'
        CALL umPrint(umMessage,src='um_shell')
        CALL nc_file_close(um_file)
        CALL release_file_unit(um_file % UNIT, handler="netcdf")
      ELSE
        WRITE(umMessage,'(A,I3,A)')'NetCDF unit ',um_file % UNIT,' is closed'
        CALL umPrint(umMessage,src='um_shell')
      END IF
      um_file => um_file_next
    END DO
  END IF

  ! Close STASH requests log file
  IF (mype == 0) THEN
    CLOSE(stash_unit)
    CALL release_file_unit(stash_unit, handler="fortran")
  END IF

  IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
    CLOSE(UNIT=eg_unit)
    CALL release_file_unit(eg_unit, handler="fortran")
  END IF


  !
  ! Close IO Server
  !
  IF (L_IOS_Active()) THEN

    IF (lstashdumptimer) THEN
      CALL timer('IOS_Shutdown', 5)
    END IF
    CALL IOS_Shutdown()
    IF (lstashdumptimer) THEN
      CALL timer('IOS_Shutdown', 6)
    END IF

    ! Shut down async stash on all ranks
    IF (isUsingAsyncStash() .OR. isUsingAsyncDumps()) THEN
      IF (lstashdumptimer) THEN
        CALL timer('IOS_stash_client_fini', 5)
      END IF
      CALL IOS_stash_client_fini()
      IF (lstashdumptimer) THEN
        CALL timer('IOS_stash_client_fini', 6)
      END IF
    END IF
  END IF

END IF ! (am an io server, atmos and io rejoin here)
9999  CONTINUE

! Namelist files were only open on PE 0, so should only be closed on such.
IF (mype == 0) THEN
  CLOSE(atmoscntl_unit)
  CLOSE(shared_unit)
END IF
CALL release_file_unit(atmoscntl_unit, handler="fortran")
CALL release_file_unit(shared_unit, handler="fortran")

CALL um_shell_banner('End')

IF (icode /= 0) THEN
  CALL Ereport(RoutineName,icode,Cmessage)
END IF

!Time a barrier to ensure dr_hook sees
!  the time for the IOS and atmos to sync up
IF (lhook) CALL dr_hook(RoutineName//':RENDEZ_VOUS',zhook_in,zhook_rendez_vous)
CALL mpl_barrier(global_comm,icode)
IF (lhook) CALL dr_hook(RoutineName//':RENDEZ_VOUS',zhook_out,zhook_rendez_vous)

IF (print_runtime_info) THEN
  IF ( L_print_pe .OR. mype == 0 ) THEN
    CALL umPrint( '',src='um_shell')
    time_end_run = get_wallclock_time()
    WRITE(umMessage,'(A,A,F10.3,A)')                         &
      'um_shell: Info: End model run',     &
      ' at time=',time_end_run - Start_time,' seconds'
    CALL umPrint(umMessage,src='um_shell')
    CALL umPrint( '',src='um_shell')
  END IF ! L_print_pe .or. mype == 0
END IF  ! print_runtime_info

CALL timer(RoutineName,2)

CALL appTerminate()

! Final top-level DrHook caliper, then prevent any further calls to DrHook.
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle) 
CALL drhook_control_disable()

! AppTerminate does not shut down GCOM for the UM (yet)...

IF ( l_oasis ) THEN
  ! DEPENDS ON: oasis_finalise
  CALL oasis_finalise
ELSE
  CALL gc_exit()
END IF

RETURN

CONTAINS

SUBROUTINE um_shell_banner(stampname)

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) ::  stampname
IF (mype == 0) THEN
  CALL DATE_AND_TIME(ch_date2, ch_time2)
  CALL umPrint('',src='um_shell')
  CALL umPrint('********************************************'//            &
                '***********************************',src='um_shell')
  WRITE(umMessage,'(23A)')                                                 &
         '**************** ',stampname,' of UM RUN Job : ',                &
         ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),                &
          ' on ',                                                          &
          ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),               &
          ' *****************'
  CALL umPrint(umMessage,src='um_shell')
  WRITE(umMessage,'(3A)')                                                  &
        '**************** Based upon UM release vn', um_version_char,      &
        '             *****************'
  CALL umPrint(umMessage,src='um_shell')
  CALL umPrint('*****************************************'//               &
               '**************************************',src='um_shell')
  CALL umPrint('',src='um_shell')
END IF

END SUBROUTINE um_Shell_banner

END SUBROUTINE um_shell
