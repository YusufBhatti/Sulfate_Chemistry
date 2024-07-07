! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module holding procedures for New UKCA emission system
!
! Method:
!   Declares the emissions data structure where the emission
!   values and other relevant data are stored. Most of the
!   emission fields are read from NetCDF files, but others
!   are calculated online.
!   Description of each subroutine contained here:
!    1. UKCA_EMISS_INIT: Identify fields in NetCDF emission files
!         and ititialise emissions structure accordingly.
!    2. UKCA_EM_STRUCT_INIT: Initialise emissions structure to
!         default values.
!    3. UKCA_EMISS_UPDATE: Check if time to update values in emissions
!         structure and do it when needed.
!    4. GET_EMFILE_REC: Get time records from a NetCDF file just
!         before and after the model time.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_emiss_mod
!
USE emiss_io_mod
USE um_types,              ONLY: integer32
USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE Control_max_sizes
USE UM_ParCore,            ONLY: mype
USE umPrintMgr,            ONLY: umPrint, umMessage, PrintStatus,         &
                                 PrStatus_Normal, PrStatus_Diag
USE t_int_mod,             ONLY: t_int
USE ukca_option_mod,       ONLY: l_ukca_qch4inter, l_ukca_prescribech4,   &
                                 ukca_em_dir, ukca_em_files, l_ukca_ibvoc,&
                                 l_ukca_primss, l_ukca_primdu,            &
                                 l_ukca_offline, l_ukca_offline_be,       &
                                 l_ukca_primbcoc, l_ukca_prim_moc,        &
                                 l_ukca_mode, nr_cdf_files
USE ukca_constants,        ONLY: L_ukca_diurnal_isopems, avc, ra, zboltz
USE ukca_d1_defs,          ONLY: n_chem_emissions, n_3d_emissions, em_chem_spec
USE nlstcall_mod,          ONLY: lcal360, ancil_reftime
USE missing_data_mod,      ONLY: rmdi, imdi
USE conversions_mod,       ONLY: pi,                             &
                                 rhour_per_day,                  &
                                 ihour_per_day,                  &
                                 rhour_per_sec,                  &
                                 rsec_per_hour

USE jules_vegetation_mod,  ONLY: l_bvoc_emis
USE asad_mod,              ONLY: advt
USE ukca_emiss_mode_mod,   ONLY: aero_ems_species, ukca_def_mode_emiss
USE ukca_mode_setup,       ONLY: nmodes, mm, rhocomp, mode, component, cp_cl, &
                                 cp_du, cp_oc
USE run_aerosol_mod,       ONLY: l_sulpc_dms

USE parkind1,              ONLY: jpim, jprb      ! DrHook
USE yomhook,               ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Default private
PRIVATE

! Emission Namelist Items
INTEGER             :: em_update_freq (nr_cdf_files) ! update freq. (hours)
INTEGER             :: em_update_type (nr_cdf_files) ! 1 serial, 2 periodic
CHARACTER (LEN=256) :: em_filename    (nr_cdf_files) ! names of source files

INTEGER         :: num_file_nml      ! Number of files in namelist
INTEGER         :: num_cdf_em_flds   ! Number of NetCDF emission fields
INTEGER         :: num_onln_em_flds  ! Number of online emission fields
INTEGER, PUBLIC :: num_em_flds       ! Number of total emission fields
INTEGER         :: num_seasalt_modes ! Number of modes used by seasalt
INTEGER         :: num_pmoc_modes    ! Number of modes used by primary marine OC
INTEGER         :: num_dust_modes    ! Number of modes used by dust

! Indices for any online emissions which are not input as NetCDF
INTEGER, PUBLIC :: inox_light          ! Index for lightning emiss of NOx
INTEGER, PUBLIC :: ich4_wetl           ! Index for wetland emiss of CH4
INTEGER, PUBLIC :: ic5h8_ibvoc         ! Index for iBVOC emiss of isoprene
INTEGER, PUBLIC :: ic10h16_ibvoc       ! Index for iBVOC emiss of monoterpenes
INTEGER, PUBLIC :: ich3oh_ibvoc        ! Index for iBVOC emiss of methanol
INTEGER, PUBLIC :: ich3coch3_ibvoc     ! Index for iBVOC emiss of acetone
INTEGER, PUBLIC :: idms_seaflux        ! Index for marine DMS emissions
INTEGER, PUBLIC :: iseasalt_first      ! Index for seasalt emiss 1st mode
INTEGER, PUBLIC :: ipmoc_first         ! Index for PMOC emiss 1st mode
INTEGER, PUBLIC :: idust_first         ! Index for dust emiss 1st mode

! Emission Data structure
TYPE, PUBLIC :: ukca_em_struct
  CHARACTER (LEN=256)  :: file_name    ! Name of source file
  CHARACTER (LEN=80)   :: var_name     ! Name of variable in file
  CHARACTER (LEN=10)   :: tracer_name  ! Emitted species
  CHARACTER (LEN=256)  :: std_name     ! standard_name attrib in NetCDF files
  CHARACTER (LEN=256)  :: lng_name     ! long_name     attrib in NetCDF files
  CHARACTER (LEN=30)   :: units        ! Units of the emission field

  INTEGER              :: update_freq  ! Update frequency (hours)
  INTEGER              :: update_type  ! 1 serial, 2 periodic, ...
  INTEGER              :: last_update  ! Num anc update intervals at last update

  LOGICAL              :: l_update     ! True if field updated in a given tstep
  ! (to indicate if conversion factors need to be applied in UKCA_NEW_EMISS_CTL)
  LOGICAL              :: l_online     ! True if field is calculated online
                                       ! (not read from file)
  LOGICAL              :: three_dim    ! True if 3D emiss field

  REAL                 :: base_fact     ! Base conv factor
  REAL, POINTER        :: vert_scaling_3d (:,:,:) ! Vertical conv factor

  CHARACTER (LEN=20)   :: hourly_fact  ! hourly_scaling: traffic,
                                                ! none, diurnal_isopems
  CHARACTER (LEN=20)   :: daily_fact   ! daily_scaling:  traffic,
                                                ! none
  CHARACTER (LEN=30)   :: vert_fact    ! vertical_scaling: surface,
                                       ! all_levels, etc.
  INTEGER              :: lowest_lev   ! Lowest and highest level where
  INTEGER              :: highest_lev  ! emiss can be injected

  ! Variables for aerosol emissions...
  LOGICAL              :: l_mode         ! Field is a num/mass mode emiss into
                                         ! which offline emissions are mapped.
  LOGICAL              :: l_mode_so2     ! Mode emission from SO2
  LOGICAL              :: l_mode_biom    ! Mode emission from biomass burning

  ! If l_mode = .TRUE. the following define the moment, mode, component the
  ! field applies to, and the name of the tracer which supplied the source
  ! emissions field.
  INTEGER              :: moment         ! Modal emission moment (0 or 3)
  INTEGER              :: mode           ! Modal emission mode
  INTEGER              :: component      ! Modal emission component
  INTEGER              :: from_emiss     ! Modal source emission index


  ! Allocatable fields of derived types are not allowed in F95 (they are
  ! an extension of Standard F95). As a consequence we use pointers for
  ! variables that need to be allocated within the emissions structure.
  REAL, POINTER:: values (:,:,:)      ! emission data
  REAL, POINTER:: diags  (:,:,:)      ! emission diagnostics
END TYPE ukca_em_struct

! Super array of emissions
TYPE (ukca_em_struct), ALLOCATABLE, PUBLIC :: emissions (:)

! 2-D field to hold biogenic isoprene emissions, with base scaling
! applied (in case units are not kg m-2 s-1), but without diurnal cycle
REAL, ALLOCATABLE, PUBLIC :: biogenic_isop (:,:)

! Subroutines available outside this module
PUBLIC :: ukca_emiss_init, ukca_emiss_update, get_emfile_rec, &
          ukca_emiss_update_mode, ukca_emiss_mode_map

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Description:
!  Get emission NetCDF files from namelist, check fields in them
!  and initialise the emissions structure.
!
! Method:
! 1) Go through all possible NetCDF emission files in the namelist and get
!    the current number of files to read (num_file_nml).
! 2) First loop through all files found in (1) to get all the emission fields
!    present and their update frequencies. Note that one file can contain
!    several emission fields but only one global attribute with the update
!    frequency, which means this freq is the same for all fields in a file.
! 3) Get:
!    * number of emission fields found in the NetCDF files: num_cdf_em_flds
!    * total number of emission fields: num_em_flds, which is equal to
!      num_cdf_em_flds + nr of online emiss (online emissions can include e.g.
!      NOx from lightning or CH4 from wetlands or mode emissions)
! 4) Allocate the emissions structure, with upper bounds equal to num_em_flds
!    and call UKCA_EM_STRUCT_INIT to initialise some items there.
! 5) Repeat loop like in (2) and check also the online emissions again
!    to fill in (or to allocate) some items in the emissions structure.
! ---------------------------------------------------------------------
SUBROUTINE ukca_emiss_init (row_length, rows,  model_levels,  &
                            global_row_length, global_rows)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN) :: row_length       ! model dimensions
INTEGER,           INTENT(IN) :: rows
INTEGER,           INTENT(IN) :: model_levels
INTEGER,           INTENT(IN) :: global_row_length
INTEGER,           INTENT(IN) :: global_rows

! Local variables
LOGICAL             :: l_nml_exis = .FALSE.
INTEGER             :: n, l, fileid, ecount
INTEGER             :: n_em_files     ! Max nr of emiss files in namelist
INTEGER, PARAMETER  :: max_dims = 5           ! Max nr of dims per var
INTEGER             :: dimen_size (max_dims)

INTEGER                    :: nvars, varid
CHARACTER (LEN=256)        :: variable_name
CHARACTER (LEN=errormessagelength) :: cmessage

! Attributes read from the NetCDF files
CHARACTER (LEN=10) :: emitted_tracer

INTEGER            :: ierror    ! Error code

LOGICAL       :: l_exist                     ! True if attribute present
LOGICAL, SAVE :: test_single_emiss = .FALSE. ! True if we intend to have one
                                             ! single emiss field per tracer

INTEGER :: icp  ! Component which a given emission field emits into
INTEGER :: imode
INTEGER :: imoment
INTEGER :: from_emiss  ! Emissions index from which a model emission derives
INTEGER :: num_modes_em
INTEGER :: num_ems_this_file  ! Number of emissions vars found in a single file
LOGICAL :: lmode_emiss(nmodes) ! True for modes which are emitted into for a
                               ! given emission field

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror   = 0

! Initialise indices for interative BVOC emissions
ic5h8_ibvoc     = -99
ic10h16_ibvoc   = -99
ich3oh_ibvoc    = -99
ich3coch3_ibvoc = -99

! Initialise variables to default values before
! reading them from namelist
em_filename    (:) = ' '       ! File names
em_update_freq (:) = imdi      ! Update frequencies
em_update_type (:) = imdi      ! Emission type: serial, periodic, ...

! Go through all possible files in the namelist and find out
! how many of them are not empty. That will be the current
! number of files to read.
num_file_nml = 0
n_em_files   = SIZE (ukca_em_files) ! max nr of emiss files in namelist
DO n = 1, n_em_files
  ! Exit when first empty file found (note that the array ukca_em_files
  ! was initialised as 'ukca_em_files is unset' in UKCA_OPTION_MOD)
  IF (ukca_em_files(n) == ' '                       .OR.   &
      ukca_em_files(n) == 'ukca_em_files is unset') EXIT

  em_filename (n) = TRIM (ukca_em_dir) // '/' // TRIM (ukca_em_files(n))
  num_file_nml    = num_file_nml + 1
END DO

! Stop with error message if no files found (i.e. if user did
! not introduce the name of any NetCDF emission file)
IF (num_file_nml == 0) THEN
  ierror   = 1
  cmessage = 'Need to provide name of at least one NetCDF emission file'
  CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
END IF

IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (umMessage,'(I3,1X,A)') num_file_nml,'files found in emission namelist'
  CALL umPrint(umMessage,src='ukca_emiss_init')
END IF

! -----------------------------------------------------------------------------
! First loop through all emission files to get their update
! frequencies (stored as global attribute) as well as the
! total number of NetCDF emission fields (num_cdf_em_flds).
num_cdf_em_flds  = 0
num_onln_em_flds = 0
DO n = 1, num_file_nml
  CALL em_fopen (TRIM(em_filename(n)),fileid)

  !  Get the update frequency of the file (in hours). Should be an integer but
  !  older files have characters. Request an integer, and lower-level routines
  !  will convert char to int if required. Integer version of routine expects an
  !  array (not scalar), so use array slice of length 1.
  CALL em_get_var_att (fileid, 1, 'GLOBAL', 'update_freq_in_hours',        &
                       em_update_freq(n:n))

  ! Make sure that update frequency is filled with an integer greater
  ! or equal to 1
  IF (em_update_freq(n)  <   1) THEN
    ierror   = 1
    cmessage = 'Update frequency needs to be multiple of one hour for ' // &
               TRIM(em_filename(n))
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  !  Similarly read the indicator for times at which data is provided.
  !  We adopt the same conventions as for ancillary files:
  !    0 Single time
  !    1 Time series
  !    2 Periodic time series
  ! This attribute should be called 'update_type', but 'emission_type' is
  ! allowed for backward compatibility. Include l_exist optional argument to
  ! allow update_type to be missing. If it is missing, 'emission_type' must be
  ! present.
  CALL em_get_var_att (fileid, 1, 'GLOBAL', 'update_type',               &
                       em_update_type(n:n), l_exist)
  IF (.NOT. l_exist) THEN
    CALL em_get_var_att (fileid, 1, 'GLOBAL', 'emission_type',               &
                         em_update_type(n:n))
    ! If no failure in the above call then emission_type exists but 
    ! update_type doesn't: issue deprecation warning.
    ierror   = -1
    cmessage = 'Attribute "emission_type" is deprecated. ' // &
        'Please use "update_type" instead. File: ' // TRIM(em_filename(n))
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  !  Check the update indicator is valid
  IF (em_update_type(n) /= 0 .AND. em_update_type(n) /= 1 .AND. &
      em_update_type(n) /= 2) THEN
    ierror   = 1
    cmessage = 'Single time (0), time series (1) and periodic (2) are ' // &
               'the only options allowed for ' // TRIM(em_filename(n))
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  !  Get number of variables in an emission file
  CALL em_get_file_info (fileid, num_vars = nvars)

  !  Go through all variables in each emission file to calculate the
  !  total number of emission fields (considering all files)
  num_ems_this_file = 0
  DO varid = 1, nvars
    ! Set variable_name to ' ' so that EM_GET_VAR_INFO returns
    ! the variable name associated with the variable ID.
    variable_name = ' '
    l = varid     ! cannot directly pass loop counter as arg declared INOUT
    CALL em_get_var_info (fileid, l, variable_name)

    ! Some variables (e.g. dimension related) will not have all the standard
    ! attributes required for emission fields. To prevent NetCDF errors for
    ! missing attributes, first check whether the variable has attribute 
    ! 'tracer_name', else continue to the next variable.
    ! This is done by passing the optional argument l_exist, which will cause
    ! em_get_var_att to return .TRUE. if attribute exists, otherwise .FALSE.
    l_exist = .FALSE.
    CALL em_get_var_att (fileid, 1, variable_name, 'tracer_name', &
                         emitted_tracer, l_exist)    

    IF (.NOT. l_exist) CYCLE ! 'tracer_name' not found, continue to next varid
    
    ! Some tracers can have their emissions calculated interactively, e.g., 
    ! isoprene (C5H8) or (mono-)terpenes. However, they can also be
    ! prescribed from ancillaries. We need to make sure that for those
    ! tracers we do not double count the emissions. Emissions can be
    ! _either_ prescribed _or_ interactive, but not both. This rule 
    ! applies only to isoprene and (mono-)terpen emissions (we introduce
    ! a separate rule for methanol and acetone later on).
    
    IF ((emitted_tracer == 'C5H8      ') .AND. L_ukca_ibvoc) THEN
      ierror   = 1
      cmessage =                                                              &
         'A biogenic emiss field for C5H8 has been read from NetCDF files ' // &
         'but interactive emissions for C5H8 have also been requested - '  // &
         'Please either deselect interactive emissions or remove the '     // &
         'ancillary file from the list of prescribed emissions.'
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    ELSE IF ((emitted_tracer == 'Monoterp  ') .AND. L_ukca_ibvoc) THEN
      ierror   = 2
      cmessage =                                                           &
         'A biogenic emiss field for Monoterp has been read from '   // &
         'NetCDF files but interactive emissions have also been '       // &
         'requested - Please either deselect interactive emissions or ' // &
         'remove the ancillary file from the list of prescribed emissions.'
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF

    !    Select only variables which are consistent with the emission definition
    !    of the chemistry scheme. For that the attribute 'tracer_name' needs to
    !    be present in the NetCDF field and the variable where that attribute is
    !    stored (i.e. emitted_tracer) needs to be equal to one of the elements
    !    in the array of emitted species (i.e. array em_chem_spec).
    IF (ANY (em_chem_spec == emitted_tracer)) THEN  
                                                       ! emitted in chem scheme
      num_ems_this_file = num_ems_this_file + 1
      num_cdf_em_flds = num_cdf_em_flds + 1
      IF (PrintStatus >= PrStatus_Normal) THEN
        WRITE (umMessage,'(A,A,A)') emitted_tracer, ' found in ',            &
                              TRIM (em_filename (n))
        CALL umPrint(umMessage,src='ukca_emiss_init')
      END IF

      ! If this is an aerosol species, create space for the number and mass it
      ! will be mapped to.  Assumes only one aerosol component for each
      ! emission field, so only one mass field required for each mode (as well
      ! as one number field per mode so 2 fields per mode).
      ! Tell ukca_def_mode_emiss to print a warning if there is a mismatch in
      ! the emission settings which could lead to unexpected behaviour.
      IF (ANY(aero_ems_species == emitted_tracer) .AND. l_ukca_mode) THEN
        CALL ukca_def_mode_emiss(emitted_tracer, lmode_emiss, icp, &
                                 lwarn_mismatch=.TRUE.)
        num_modes_em = COUNT(lmode_emiss)
        num_onln_em_flds = num_onln_em_flds + 2*num_modes_em
      END IF

    ELSE
      ! Report a warning and ignore the emission field if it is not present
      ! in the current chemistry scheme. Note that this error will appear if any
      ! fields in the netcdf file have a tracer_name atribute containing an
      ! empty string.
      ierror   = -1
      cmessage = 'emission field ' // TRIM (emitted_tracer) //       &
        ' not present in the chemistry scheme and therefore ignored'
      CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
    END IF       ! emitted in chem scheme

  END DO         ! loop trhough variables in an emission file

  IF (num_ems_this_file == 0) THEN
    ! File contains no emissions for this chemistry scheme, i.e. no
    ! variable contains tracer_name attribute which matches em_chem_spec.
    ierror   = n
    cmessage = 'Emission file ' // TRIM(em_filename(n)) //       &
      ' contains no emission variables for this chemistry scheme.'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  ELSE
    IF (PrintStatus >= PrStatus_Normal) THEN
      WRITE (umMessage,'(I4,2A)') num_ems_this_file, &
          " emission fields will be used from ", TRIM(em_filename(n))
      CALL umPrint(umMessage,src='ukca_emiss_init')
    END IF
  END IF

  CALL em_fclose(fileid)
END DO    ! loop through all NetCDF emiss files

! -----------------------------------------------------------------------------
! The total nr of emission fields will be equal to the number of fields
! found in the NetCDF files (i.e. num_cdf_em_flds) plus any other fields
! that account for online emissions.
!
! Add first an index for online NOx emiss from lightning,
! which are always present except in the offline oxidants mechanism
IF (.NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
  inox_light = num_cdf_em_flds + num_onln_em_flds + 1   
  num_onln_em_flds = num_onln_em_flds + 1
END IF
!
! Add index for online CH4 from wetland emissions if these are
! active and prescribed surface CH4 concentrations are not used.
! Also increment the index for total number of emissions
IF (L_ukca_qch4inter .AND. (.NOT. L_ukca_prescribech4)) THEN
  ich4_wetl        = num_cdf_em_flds  + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 1
ELSE
  ich4_wetl = -99
END IF

! Interactive emissions of biogenic volatile organic compounds (iBVOC)
IF (l_ukca_ibvoc .AND. l_bvoc_emis) THEN

  IF (ANY (em_chem_spec == 'C5H8      ')) THEN
    ic5h8_ibvoc      = num_cdf_em_flds  + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ! Report a warning to indicate that it will not be used
    ! because species not emitted in the chemistry scheme
    ierror      = -2
    cmessage    = 'Biogenic emission field of C5H8 '        // &
                  'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'Monoterp  ')) THEN
    ic10h16_ibvoc    = num_cdf_em_flds  + num_onln_em_flds + 1
    num_onln_em_flds = num_onln_em_flds + 1
  ELSE
    ierror        = -3
    cmessage      = 'Biogenic emission field of Monoterp '  // &
                    'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  ! Note that methanol can be referred to with two different names
  ! depending on the chemistry scheme.
  IF (ANY (em_chem_spec == 'MeOH      ') .OR.    &
      ANY (em_chem_spec == 'CH3OH     ')) THEN
!     ich3oh_ibvoc     = num_cdf_em_flds  + num_onln_em_flds + 1
!     num_onln_em_flds = num_onln_em_flds + 1
    ich3oh_ibvoc  = -99
  ELSE
    ierror        = -4
    cmessage      = 'Biogenic emission field of methanol '  // &
                    'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF

  IF (ANY (em_chem_spec == 'Me2CO     ')) THEN
!     ich3coch3_ibvoc  = num_cdf_em_flds  + num_onln_em_flds + 1
!     num_onln_em_flds = num_onln_em_flds + 1
    ich3coch3_ibvoc  = -99
  ELSE
    ierror          = -5
    cmessage        = 'Biogenic emission field of Me2CO '   // &
                      'ignored for this chemistry scheme'
    CALL ereport('UKCA_EMISS_INIT', ierror, cmessage)
  END IF
END IF

! Add index for online marine DMS emissions
IF (ANY(em_chem_spec == 'DMS     ') .AND. (.NOT. l_sulpc_dms)) THEN
  idms_seaflux = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 1
ELSE
  idms_seaflux = -99
END IF

! Add index for online seasalt emissions
! Make space for mass and number for every seasalt mode
IF (l_ukca_primss) THEN
  num_seasalt_modes = COUNT(component(:,cp_cl))
  iseasalt_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_seasalt_modes
ELSE
  iseasalt_first = -99
END IF
! Add index for online marine primary organic emissions
! Make space for number and mass for all OC modes
IF (l_ukca_prim_moc) THEN
  IF (.NOT. l_ukca_primss) THEN
    cmessage = 'Prim. marine OC requires sea-salt'          // &
               'emissions to be turned on'
    ierror = 1
    CALL ereport('UKCA_EMISS_INIT',ierror,cmessage)
  END IF
  num_pmoc_modes = COUNT(component(:,cp_oc))
  ipmoc_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_pmoc_modes
ELSE
  ipmoc_first = -99
END IF
! Add index for online dust emissions
IF (l_ukca_primdu) THEN
  num_dust_modes = COUNT(component(:,cp_du))
  idust_first = num_cdf_em_flds + num_onln_em_flds + 1
  num_onln_em_flds = num_onln_em_flds + 2*num_dust_modes
ELSE
  idust_first = -99
END IF


num_em_flds = num_cdf_em_flds + num_onln_em_flds
IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE (umMessage,'(A,1X,I3)') 'Number of NetCDF emission fields:', &
                                num_cdf_em_flds
  CALL umPrint(umMessage,src='ukca_emiss_init')
  WRITE (umMessage,'(A,1X,I3)') 'Total number of emission fields:',  &
                                 num_em_flds
  CALL umPrint(umMessage,src='ukca_emiss_init')
END IF

! Allocate the emission structure and initialise some
! of its elements to default values
ALLOCATE (emissions(num_em_flds))
CALL ukca_em_struct_init

! -----------------------------------------------------------------------------
! Now repeat loops through emission files and fields, with the  objective
! to fill in (or allocate) some items in the emissions structure

ecount = 0   ! initialise counter for elements in emissions tructure
DO n = 1, num_file_nml
  CALL em_fopen (TRIM(em_filename(n)),fileid)

  !  Get number of variables in each emission file
  CALL em_get_file_info (fileid, num_vars = nvars)

  !  Select the right emission fields as before and store
  !  some of their attributes in the emissions structure
  DO varid = 1, nvars
    ! Set variable_name to ' ' so that EM_GET_VAR_INFO returns
    ! the variable name associated with the variable ID.
    variable_name = ' '
    l = varid       ! cannot pass loop counter as arg declared INOUT
    CALL em_get_var_info (fileid, l, variable_name)

    ! As before, do not look for attributes for dimension variables
    ! (this will be the case when attribute 'tracer_name' is absent)
    l_exist = .FALSE.
    CALL em_get_var_att (fileid, 1, variable_name, 'tracer_name', &
                         emitted_tracer, l_exist)

    IF (.NOT. l_exist) CYCLE  ! 'tracer_name' not found, continue to next varid
    
    !    If an emission field is found then do some updates for it
    !    in the emissions structure
    IF (ANY (em_chem_spec == emitted_tracer)) THEN  ! emitted in chem scheme
      ecount = ecount + 1 ! increase counter for field in emissions struct

      ! Populate the structure from the netCDF attributes
      CALL ukca_em_struct_set(ecount, varid, n, fileid,         &
                              row_length, rows,  model_levels,  &
                              global_row_length, global_rows)

      emissions(ecount)%tracer_name = emitted_tracer

      ! For aerosol emissions, additionally populate structures for individual
      ! mode mass/num emissions into which the offline emission field will be
      ! mapped.
      ! We assume only one aerosol component per emission field.
      IF (ANY(aero_ems_species == emitted_tracer) .AND. l_ukca_mode) THEN
        CALL ukca_def_mode_emiss(emitted_tracer, lmode_emiss, icp)
        from_emiss = ecount
        DO imode = 1, nmodes
          IF (lmode_emiss(imode)) THEN
            DO imoment = 0, 3, 3  ! imoment = 0 and 3 only
              ecount = ecount + 1
              ! Set struct as for normal emissions, then set mode attrs.
              ! For mass emissions units are same as in file, but for number
              ! they need overwritten. std_name/lng_name will be wrong for 
              ! number but these are not used in the code.
              ! Number flux is multiplied by MM_DA/AVC to give units of
              ! kg(air) m-2 s-1
              ! If emitted_tracer contains substring SO2 or biomass, set
              ! logicals which control scaling in ukca_add_emiss
              CALL ukca_em_struct_set(ecount, varid, n, fileid,         &
                                      row_length, rows,  model_levels,  &
                                      global_row_length, global_rows)
              emissions(ecount)%l_mode      = .TRUE.
              emissions(ecount)%l_online    = .TRUE.
              emissions(ecount)%mode        = imode
              emissions(ecount)%moment      = imoment
              emissions(ecount)%component   = icp
              emissions(ecount)%tracer_name = 'mode_emiss'
              emissions(ecount)%from_emiss  = from_emiss
              IF (imoment == 0) emissions(ecount)%units = 'kg(air) m-2 s-1'
              IF (INDEX(emitted_tracer,'SO2') /= 0) THEN
                emissions(ecount)%l_mode_so2 = .TRUE.
              ELSE IF (INDEX(emitted_tracer,'biomass') /= 0 ) THEN
                emissions(ecount)%l_mode_biom = .TRUE.
              END IF
            END DO
          END IF
        END DO
      END IF  ! offline aerosol emission

    END IF    ! emitted in chem scheme
  END DO      ! loop through variables in an emission file

  CALL em_fclose(fileid)
END DO        ! loop through all NetCDF emission fields


! Initialise fields for lightning NOx in the emissions structure
! (Note that there is no need to update file_name, update_freq
! and update_type for online emissions)
IF (.NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
  emissions(inox_light)%var_name    = 'NOx_lightning_emissions'
  emissions(inox_light)%tracer_name = 'NO_lightng'
  emissions(inox_light)%l_online    = .TRUE.

  ! Units should be consistent with the emission field we will take from
  ! the call to ukca_light_ctl within ukca_new_emiss_ctl, i.e.
  ! 'kg(N) gridbox-1 s-1'. However, to follow the same approach as
  ! for NetCDF emission fields, we do not indicate 'N' in the units
  ! but use the substring 'expressed as nitrogen' in the long_name.
  emissions(inox_light)%units       = 'kg gridbox-1 s-1'
  emissions(inox_light)%lng_name    = 'tendency of atmosphere mass content' // &
    ' of nitrogen monoxide expressed as nitrogen due to emission from lightning'

  emissions(inox_light)%hourly_fact = 'none'         ! No need to apply
  emissions(inox_light)%daily_fact  = 'none'         ! time profiles

  emissions(inox_light)%vert_fact   = 'all_levels'   ! 3D emissions
  emissions(inox_light)%three_dim   = .TRUE.
  emissions(inox_light)%lowest_lev  = 1
  emissions(inox_light)%highest_lev = model_levels

  ! Emission values and diagostics allocated and initialised here.
  ! They will be updated in ukca_new_emiss_ctl and ukca_add_emiss,
  ! respectively
  ALLOCATE (emissions(inox_light)%values(row_length, rows, model_levels))
  ALLOCATE (emissions(inox_light)%diags (row_length, rows, model_levels))
  emissions(inox_light)%values(:,:,:) = 0.0
  emissions(inox_light)%diags (:,:,:) = 0.0
END IF

! Initialise fields for online wetland emiss of CH4 in emissions structure
! (done in a similar way as done above for lightning NOx, but only if
! these emissions exist)
IF (L_ukca_qch4inter .AND. (.NOT. L_ukca_prescribech4)) THEN
  emissions(ich4_wetl)%var_name    = 'CH4_wetland_emissions'
  emissions(ich4_wetl)%tracer_name = 'CH4_wetlnd'
  emissions(ich4_wetl)%l_online    = .TRUE.

  !  Units should be consistent with CH4 wetland emiss from JULES,
  !  given in 'kg(C) m-2 s-1'. Similarly as seen above for 'NO_lightng',
  !  we do not indicate 'C' in the units but use the substring
  !  'expressed as carbon' in the long_name.
  emissions(ich4_wetl)%units      = 'kg m-2 s-1'
  emissions(ich4_wetl)%lng_name   = 'tendency of atmosphere mass content' // &
     ' of methane expressed as carbon due to emission from wetlands'

  emissions(ich4_wetl)%hourly_fact = 'none'      ! No need to apply
  emissions(ich4_wetl)%daily_fact  = 'none'      ! time profiles

  emissions(ich4_wetl)%vert_fact   = 'surface'   ! surface emissions
  emissions(ich4_wetl)%three_dim   = .FALSE.
  emissions(ich4_wetl)%lowest_lev  = 1
  emissions(ich4_wetl)%highest_lev = 1

  ALLOCATE (emissions(ich4_wetl)%values(row_length, rows, 1))
  ALLOCATE (emissions(ich4_wetl)%diags (row_length, rows, 1))
  emissions(ich4_wetl)%values(:,:,:) = 0.0
  emissions(ich4_wetl)%diags (:,:,:) = 0.0
END IF

! Initialise fields for iBVOC in emissions structure
IF (l_ukca_ibvoc .AND. l_bvoc_emis) THEN
  IF (ic5h8_ibvoc > 0) THEN
    emissions(ic5h8_ibvoc)%var_name    = 'interactive_C5H8'
    emissions(ic5h8_ibvoc)%tracer_name = 'C5H8'
    emissions(ic5h8_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ic5h8_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ic5h8_ibvoc)%lng_name   =                                     &
     'tendency of atmosphere mass content'                               // &
     ' of isoprene expressed as carbon due to emission from vegetation'

    emissions(ic5h8_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ic5h8_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ic5h8_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ic5h8_ibvoc)%three_dim   = .FALSE.
    emissions(ic5h8_ibvoc)%lowest_lev  = 1
    emissions(ic5h8_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ic5h8_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic5h8_ibvoc)%diags (row_length, rows, 1))
    emissions(ic5h8_ibvoc)%values(:,:,:) = 0.0
    emissions(ic5h8_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ic10h16_ibvoc > 0) THEN
    emissions(ic10h16_ibvoc)%var_name    = 'interactive_terpene'
    emissions(ic10h16_ibvoc)%tracer_name = 'Monoterp'
    emissions(ic10h16_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ic10h16_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ic10h16_ibvoc)%lng_name   =                                   &
     'tendency of atmosphere mass content'                               // &
     ' of terpenes expressed as carbon due to emission from vegetation'

    emissions(ic10h16_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ic10h16_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ic10h16_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ic10h16_ibvoc)%three_dim   = .FALSE.
    emissions(ic10h16_ibvoc)%lowest_lev  = 1
    emissions(ic10h16_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ic10h16_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ic10h16_ibvoc)%diags (row_length, rows, 1))
    emissions(ic10h16_ibvoc)%values(:,:,:) = 0.0
    emissions(ic10h16_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3oh_ibvoc > 0) THEN
    emissions(ich3oh_ibvoc)%var_name    = 'interactive_CH3OH'

    ! Note that methanol can be referred to with two different names
    ! depending on the chemistry scheme.
    IF (ANY (em_chem_spec == 'MeOH      ')) THEN
      emissions(ich3oh_ibvoc)%tracer_name = 'MeOH'
    ELSE IF (ANY (em_chem_spec == 'CH3OH     ')) THEN
      emissions(ich3oh_ibvoc)%tracer_name = 'CH3OH'
    ELSE
      ierror   = 1
      cmessage = 'iBVOC field for methanol found in emission structure ' // &
                 ' but not found in em_chem_spec'
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF

    emissions(ich3oh_ibvoc)%l_online   = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ich3oh_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ich3oh_ibvoc)%lng_name   =                                    &
     'tendency of atmosphere mass content'                               // &
     ' of methanol expressed as carbon due to emission from vegetation'

    emissions(ich3oh_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ich3oh_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ich3oh_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ich3oh_ibvoc)%three_dim   = .FALSE.
    emissions(ich3oh_ibvoc)%lowest_lev  = 1
    emissions(ich3oh_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ich3oh_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3oh_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3oh_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3oh_ibvoc)%diags (:,:,:) = 0.0
  END IF

  IF (ich3coch3_ibvoc > 0) THEN
    emissions(ich3coch3_ibvoc)%var_name    = 'interactive_Me2CO'
    emissions(ich3coch3_ibvoc)%tracer_name = 'Me2CO'
    emissions(ich3coch3_ibvoc)%l_online    = .TRUE.

    ! iBVOC emissions from JULES are in 'kg(C) m-2 s-1'. Don't indicate this
    ! in %units but use substring 'expressed as carbon' in %lng_name.
    emissions(ich3coch3_ibvoc)%units      = 'kg m-2 s-1'
    emissions(ich3coch3_ibvoc)%lng_name   =                                 &
     'tendency of atmosphere mass content'                               // &
     ' of acetone expressed as carbon due to emission from vegetation'

    emissions(ich3coch3_ibvoc)%hourly_fact = 'none'      ! No need to apply
    emissions(ich3coch3_ibvoc)%daily_fact  = 'none'      ! time profiles

    emissions(ich3coch3_ibvoc)%vert_fact   = 'surface'   ! surface emissions
    emissions(ich3coch3_ibvoc)%three_dim   = .FALSE.
    emissions(ich3coch3_ibvoc)%lowest_lev  = 1
    emissions(ich3coch3_ibvoc)%highest_lev = 1

    ALLOCATE (emissions(ich3coch3_ibvoc)%values(row_length, rows, 1))
    ALLOCATE (emissions(ich3coch3_ibvoc)%diags (row_length, rows, 1))
    emissions(ich3coch3_ibvoc)%values(:,:,:) = 0.0
    emissions(ich3coch3_ibvoc)%diags (:,:,:) = 0.0
  END IF
END IF

! Initialise online marine DMS emissions 
IF (ANY(em_chem_spec == 'DMS     ') .AND. (.NOT. l_sulpc_dms)) THEN
  emissions(idms_seaflux)%var_name    = 'DMS_marine_emissions'
  emissions(idms_seaflux)%tracer_name = 'DMS       '
  emissions(idms_seaflux)%l_online  = .TRUE.

  !  Units should be consistent with DMS emiss from DMS_FLUX_4A,
  !  given in 'kg(S) m-2 s-1'. Similarly as seen above for 'NO_lightng',
  !  we do not indicate 'S' in the units but use the substring
  !  'expressed as sulfur' in the long_name.
  emissions(idms_seaflux)%units      = 'kg m-2 s-1'
  emissions(idms_seaflux)%lng_name   = 'tendency of atmosphere mass ' // &
     'content of dimethyl sulfide expressed as sulfur due to emission from sea'

  emissions(idms_seaflux)%hourly_fact = 'none'      ! No need to apply
  emissions(idms_seaflux)%daily_fact  = 'none'      ! time profiles

  emissions(idms_seaflux)%vert_fact   = 'surface'   ! surface emissions
  emissions(idms_seaflux)%three_dim   = .FALSE.
  emissions(idms_seaflux)%lowest_lev  = 1
  emissions(idms_seaflux)%highest_lev = 1

  ALLOCATE (emissions(idms_seaflux)%values(row_length, rows, 1))
  ALLOCATE (emissions(idms_seaflux)%diags (row_length, rows, 1))
  emissions(idms_seaflux)%values(:,:,:) = 0.0
  emissions(idms_seaflux)%diags (:,:,:) = 0.0
END IF

! Initialise fields for seasalt emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_ss.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=iseasalt_first to identify as seasalt for ukca_emiss_mode_map
IF (l_ukca_primss) THEN
  ecount = iseasalt_first - 1
  DO imode = 1, nmodes
    DO imoment = 0, 3, 3  ! imoment = 0 or 3
      IF (mode(imode) .AND. component(imode,cp_cl)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'seasalt_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = iseasalt_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_cl

        IF (imoment == 0) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' // &
                     'content_of_seasalt_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primss

! Initialise fields for primary marine organic carbon emissions (PMOC).
! Allow emission into any mode for flexibility. Emissions are 2D and emit into
! OC component.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=ipmao_first to identify as pmao for ukca_emiss_mode_map
IF (l_ukca_prim_moc) THEN
  ecount = ipmoc_first - 1
  DO imode = 1, nmodes
    DO imoment = 0, 3, 3  ! imoment = 0 or 3
      IF (mode(imode) .AND. component(imode,cp_oc)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'pmoc_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = ipmoc_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_oc

        IF (imoment == 0) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' // &
           'content_of_particulate_organic_matter_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_prim_moc

! Initialise fields for dust emissions.
! There is a 2D mass and number field for every mode, consistent with the
! shape of the arrays returned by ukca_prim_du.
! Number flux is multiplied by MM_DA/AVC to give units of kg(air) m-2 s-1
! Set from_emiss=idust_first to identify as dust for ukca_emiss_mode_map
IF (l_ukca_primdu) THEN
  ecount = idust_first - 1
  DO imode = 1, nmodes
    DO imoment = 0, 3, 3  ! imoment = 0 or 3
      IF (mode(imode) .AND. component(imode,cp_du)) THEN
        ecount = ecount + 1
        emissions(ecount)%var_name    = 'dust_online_emissions'
        emissions(ecount)%tracer_name = 'mode_emiss'
        emissions(ecount)%from_emiss  = idust_first

        emissions(ecount)%l_mode    = .TRUE.
        emissions(ecount)%l_online  = .TRUE.
        emissions(ecount)%mode      = imode
        emissions(ecount)%moment    = imoment
        emissions(ecount)%component = cp_du

        IF (imoment == 0) THEN
          emissions(ecount)%units      = 'kg(air) m-2 s-1'
        ELSE
          emissions(ecount)%units      = 'kg m-2 s-1'
        END IF
        emissions(ecount)%std_name  = 'tendency_of_atmosphere_mass_' // &
                     'content_of_dust_dry_aerosol_due_to_emission'

        emissions(ecount)%hourly_fact = 'none'      ! No need to apply
        emissions(ecount)%daily_fact  = 'none'      ! time profiles

        emissions(ecount)%vert_fact   = 'surface'   ! surface emissions
        emissions(ecount)%three_dim   = .FALSE.
        emissions(ecount)%lowest_lev  = 1
        emissions(ecount)%highest_lev = 1

        ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
        ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
        emissions(ecount)%values(:,:,:) = 0.0
        emissions(ecount)%diags (:,:,:) = 0.0
      END IF  ! mode and component in use
    END DO  ! imoment
  END DO  ! imode
END IF  ! l_ukca_primdu





! If the emissions structure contains an offline biogenic isoprene field
! that needs a diurnal cycle to be applied online (via the routine
! UKCA_DIURNAL_ISOP_EMS) then create a 2-D field which will hold
! such emissions before applying the diurnal cycle.
! That field will be updated only when emission values are
! updated, depending on the updating frequency, but will not
! include the diurnal correction to avoid re-applying it
! (see code in UKCA_NEW_EMISS_CTL).
! Also avoid double-counting emissions: Do not allow the use of
! iBVOC emissions for isoprene (from JULES) when an offline
! biogenic emission field for the same species has been read
! from the NetCDF files.
IF (ANY (emissions(:)%tracer_name == 'C5H8      '  .AND. &
    emissions(:)%hourly_fact == 'diurnal_isopems') .AND. &
    L_ukca_diurnal_isopems                         .AND. &
    .NOT. L_ukca_ibvoc) THEN

  ALLOCATE (biogenic_isop (row_length, rows))
  biogenic_isop(:,:) = 0.0

END IF

! Make sure that there are no missing fields in the emission files,
! i.e. that all names in em_chem_spec are present in emissions structure
! An exception is that OC/BC emissions can be missing if l_ukca_primbcoc is 
! false. The same exception doesn't apply to SO2 with l_ukca_primsu, since 
! SO2 is emitted into the chemistry.
! NVOC is also an exception if interactive BVOCs are used - In that case
! biogenic emissions of some species such as methanol and acetone are
! provided by the iBVOC code. Therefore they do not need to come from the
! NVOC (GEIA) emission file which also represents biogenic emissions.
DO n = 1, n_chem_emissions + n_3d_emissions
  IF (.NOT. ANY(emissions(:)%tracer_name == em_chem_spec (n)) ) THEN 
    IF (.NOT. (                                                    &
        ((em_chem_spec(n)(1:2) == 'BC') .AND. .NOT. l_ukca_primbcoc) .OR.  &
        ((em_chem_spec(n)(1:2) == 'OC') .AND. .NOT. l_ukca_primbcoc) .OR.  &
        ((em_chem_spec(n) == 'NVOC      ') .AND. l_ukca_ibvoc) ) ) THEN
      cmessage = TRIM (em_chem_spec (n)) // ' missing in the emission files'
      ierror = n
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF
  END IF
END DO

! The switch test_single_emiss can be set to TRUE if for whatever reason
! it is intended that there is only one emission field for each tracer.
! The code below checks that there is no duplication of names in the final
! emissions structure for that particular case.
! Note however that the switch will be FALSE by default, because some
! chemistry schemes will use several emission fields for most tracers
! in the future to account for different emission source sectors.
IF (test_single_emiss) THEN
  DO n = 1, num_em_flds-1
    IF ( ANY (emissions(n+1:num_em_flds)%tracer_name ==                      &
              emissions(n)%tracer_name) ) THEN
      cmessage = 'Duplicate emission fields for ' // emissions(n)%tracer_name
      ierror = n
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF
  END DO
END IF

! Check that all specified NetCDF files are valid and that dimensions match
! (just double-checking, but not strictly needed because similar checks were
! done before when looking for the nr of vertical levels in each emiss field).
! Also check that all elements of the emissions array have been filled.
IF ( PrintStatus >= PrStatus_Diag ) THEN
  DO n = 1, num_em_flds
    IF (emissions(n)%tracer_name == ' ' .OR. &
        emissions(n)%units == ' ') THEN
      cmessage = 'Unfilled emission field (index is errcode)'
      ierror = n 
      CALL ereport ('UKCA_EMISS_INIT', ierror, cmessage)
    END IF
    IF (emissions(n)%l_online) CYCLE  ! only do checks for offline ems

    WRITE (umMessage,'(A,A,A,A)')  'Verifying File: ',                    &
        TRIM(emissions(n)%file_name), ' for ',                            &
        TRIM(emissions(n)%var_name)
    CALL umPrint(umMessage,src='ukca_emiss_init')

    ! Aborts internally if problems
    CALL em_fopen (TRIM (emissions(n)%file_name), fileid)
    CALL em_var_check_dims (global_row_length, global_rows, model_levels, &
                           fileid, emissions(n)%var_name, dimen_size)

    CALL em_fclose(fileid)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END SUBROUTINE ukca_emiss_init

!---------------------------------------------------------------------------
! Description:
!   Initialisation of emissions structure to default values.
! Method:
!   Initialise elements of the emissions structure, which has been
!   previously allocated with upper bound equal to 'num_em_flds'.
!   Note that as indicated below some elements cannot be
!   allocated yet within this subroutine.
!---------------------------------------------------------------------------
SUBROUTINE ukca_em_struct_init

IMPLICIT NONE

INTEGER                        :: ierror         ! Error code
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! For Dr Hook
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EM_STRUCT_INIT'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierror = 0

IF ( .NOT. ALLOCATED(emissions)) THEN
  ierror = 1
  CALL ereport('UKCA_EM_STRUCT_INIT', ierror,                       &
    'Emissions super array should be allocated before this routine')
END IF

emissions(:)%file_name   = ' '
emissions(:)%var_name    = ' '
emissions(:)%tracer_name = ' '
emissions(:)%std_name    = ' '
emissions(:)%lng_name    = ' '
emissions(:)%units       = ' '

emissions(:)%update_freq = 0
emissions(:)%update_type = 0
emissions(:)%l_update    = .FALSE.
emissions(:)%l_online    = .FALSE.
emissions(:)%three_dim   = .FALSE.

emissions(:)%base_fact   = 1.0 ! by default do not apply any conversion factor

emissions(:)%hourly_fact = ' '
emissions(:)%daily_fact  = ' '
emissions(:)%vert_fact   = ' '
emissions(:)%lowest_lev  = 1   ! Initially assume surface emiss (changed later
emissions(:)%highest_lev = 1   ! in the code for both 3D and high level emiss)

emissions(:)%l_mode      = .FALSE.
emissions(:)%l_mode_so2  = .FALSE.
emissions(:)%l_mode_biom = .FALSE.
emissions(:)%moment      = -1
emissions(:)%mode        = -1
emissions(:)%component   = -1
emissions(:)%from_emiss  = -1

! Note that the pointers emissions(:)%values and emissions(:)%diags
! cannot be initialised yet. Before that they need to be allocated
! independently for each element of emissions(:). This is done in
! ukca_emiss_init after calling this routine.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_em_struct_init


!---------------------------------------------------------------------------
! Description:
!   Get values of core attributes for the emission field at index ecount
!
! Method:
!   Read attributes from file and then set emissions structure attributes.
!   Check the consistency of the vertical scaling information.
!   Allocate arrays for storing emission values and diagnostics.
!---------------------------------------------------------------------------
SUBROUTINE ukca_em_struct_set(ecount, varid, n, fileid,         &
                              row_length, rows,  model_levels,  &
                              global_row_length, global_rows)

IMPLICIT NONE

! Subroutine arguments
INTEGER,           INTENT(IN) :: fileid
INTEGER,           INTENT(IN) :: row_length       ! model dimensions
INTEGER,           INTENT(IN) :: rows
INTEGER,           INTENT(IN) :: model_levels
INTEGER,           INTENT(IN) :: global_row_length
INTEGER,           INTENT(IN) :: global_rows

INTEGER, INTENT(IN) :: ecount  ! index of emission field to populate
INTEGER, INTENT(IN) :: varid  ! netCDF variable ID

INTEGER, INTENT(IN) :: n ! current index in em_filename array

! Attributes read from the netCDF files
CHARACTER (LEN=256)        :: variable_name
CHARACTER (LEN=256)        :: standard_name, long_name
CHARACTER (LEN=30)         :: emiss_units
CHARACTER (LEN=20)         :: hourly_scaling, daily_scaling
CHARACTER (LEN=30)         :: vert_scaling

INTEGER                    :: lowest_level,  highest_level
INTEGER                    :: levelarr(1) ! array to read int attrib
INTEGER                    :: ierror      ! Error code

INTEGER                    :: varid_local
INTEGER (KIND=integer32)   :: var_type, num_atts, iatt
INTEGER, PARAMETER         :: max_dims = 5         ! Max nr of dims per var
INTEGER                    :: dimen_size (max_dims)

! Flag to check presence of an attrib
LOGICAL                    :: l_exist  
! Flags for presence of standard_name & long_name attributes
LOGICAL                    :: l_exist_sn, l_exist_ln 

CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EM_STRUCT_SET'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! The argument varid is INTENT(IN). Need to use local variable
! to pass it as INOUT argument.
varid_local   = varid
! Set variable_name to ' ' so that EM_GET_VAR_INFO returns
! the variable name associated with the variable ID.
variable_name = ' '
CALL em_get_var_info (fileid, varid_local, variable_name)

! Read mandatory attributes from this variable.
CALL em_get_var_att (fileid, 1, variable_name, 'units', emiss_units)

! Read some optional attributes from the variable. This is done by
! passing the INOUT argument l_exist_x, which has input value = .FALSE.
! The output value will be changed to .TRUE. only if att exists.
! Atts which are missing will remain with the initialised value of ' '.
standard_name  = ' '
long_name      = ' '
hourly_scaling = ' '
daily_scaling  = ' '
vert_scaling   = ' '

! First look for standard_name and long_name and make sure that
! at least one of them is present (to comply with CF conventions).
l_exist_sn = .FALSE.
CALL em_get_var_att (fileid, 1, variable_name, 'standard_name',           &
                                                standard_name, l_exist_sn)

l_exist_ln = .FALSE.
CALL em_get_var_att (fileid, 1, variable_name, 'long_name',               &
                                                long_name,     l_exist_ln)

IF ( (.NOT. l_exist_sn) .AND. (.NOT. l_exist_ln) ) THEN
    ierror   = 1
    cmessage = 'standard_name and long_name missing in '      //          &
               TRIM (variable_name)                           //          &
               ': at least one of them needs to be present'
    CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
END IF

! Look for attributes related to vertical / temporal profiles of emissions.
! If not present in the files they will remain as ' ' and no profiles
! will be applied.
l_exist = .FALSE.
CALL em_get_var_att (fileid, 1, variable_name, 'hourly_scaling',          &
                                                hourly_scaling,  l_exist)
l_exist = .FALSE.
CALL em_get_var_att (fileid, 1, variable_name, 'daily_scaling',           &
                                                daily_scaling,   l_exist)
l_exist = .FALSE.
CALL em_get_var_att (fileid, 1, variable_name, 'vertical_scaling',        &
                                                vert_scaling,    l_exist)

! Fill in file_name and update_freq with information read from the
! namelist (these two variables are the same for all fields within
! the same emision file).
emissions(ecount)%file_name   = em_filename    (n)
emissions(ecount)%update_freq = em_update_freq (n) ! hours
emissions(ecount)%update_type = em_update_type (n)

! The other elements of the structure can be read
! from the variables/attributes in the emission files
emissions(ecount)%var_name    = variable_name
emissions(ecount)%std_name    = standard_name
emissions(ecount)%lng_name    = long_name
emissions(ecount)%units       = emiss_units
emissions(ecount)%hourly_fact = hourly_scaling
emissions(ecount)%daily_fact  = daily_scaling
emissions(ecount)%vert_fact   = vert_scaling

! Note that by default ukca_em_struct_init has set:
!    emissions(:)%base_fact   = 1.
! The routine ukca_emiss_factors will update this
! based on the units of the emission field

! If it is a high or a single level emission field then the
! NetCDF file should contain two attributes indicating the
! highest and lowest level where emissions should be injected.
! em_get_var_att_int requires an array as the argument, so the
! length 1 array levelarr is used as a dummy argument, then read
! into the corresponding scalar integer.
IF (vert_scaling == 'high_level'    .OR.                             &
   vert_scaling == 'single_level') THEN

  CALL em_get_var_att (fileid, 1, variable_name, 'lowest_level',     & 
                                                 levelarr(1:1))
  lowest_level = levelarr(1)

  CALL em_get_var_att (fileid, 1, variable_name, 'highest_level',    &
                                                 levelarr(1:1))
  highest_level = levelarr(1)

  IF (lowest_level  < 1 .OR. lowest_level  > model_levels) THEN
    ierror   = 1
    cmessage = 'lowest level for ' // TRIM (variable_name) //         &
               ' outside the range of model vertical levels'
    CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
  END IF

  IF (highest_level < 1 .OR. highest_level > model_levels) THEN
    ierror   = 1
    cmessage = 'highest level for ' // TRIM (variable_name) //        &
               ' outside the range of model vertical levels'
    CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
  END IF

  IF (lowest_level > highest_level) THEN
    ierror   = 1
    cmessage = 'lowest level has to be less than or equal '//         &
               'to highest level for ' // TRIM (variable_name)
    CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
  END IF

  IF (vert_scaling == 'single_level' .AND.                            &
      lowest_level /= highest_level) THEN
    ierror   = 1
    cmessage = 'highest and lowest level should coincide for ' //     &
               'single level emissions ('                      //     &
               TRIM (variable_name) // ')'
    CALL ereport('UKCA_EM_STRUCT_SET', ierror, cmessage)
  END IF

  emissions(ecount)%lowest_lev  = lowest_level
  emissions(ecount)%highest_lev = highest_level
END IF  ! 'high_level' OR 'single_level'

! Checking number of dimensions of the emission field
CALL em_var_check_dims (global_row_length, global_rows, model_levels, &
                        fileid, variable_name, dimen_size)

! Different allocations depending on the nr of elements
! in the vertical dimension
IF (dimen_size (3) == 1) THEN                ! 2D emiss field
  ALLOCATE (emissions(ecount)%values(row_length, rows, 1))
  ALLOCATE (emissions(ecount)%diags (row_length, rows, 1))
  !   If 'high_level' or 'single_level' emission then the highest and
  !   lowest levels where to inject the emissions should have been
  !   updated above. If not then assume that it is a surface emiss
  !   field and therefore set both of them to 1 (just for completeness)
  IF (vert_scaling /= 'high_level'    .AND.                         &
      vert_scaling /= 'single_level') THEN
    emissions(ecount)%lowest_lev  = 1
    emissions(ecount)%highest_lev = 1
  END IF
ELSE
  IF (dimen_size (3) == model_levels) THEN   ! 3D emiss field
    ALLOCATE (emissions(ecount)%values(row_length, rows, model_levels))
    ALLOCATE (emissions(ecount)%diags (row_length, rows, model_levels))
    ! Update the lowest and highest levels where emiss are injected
    ! (just for consistency) and indicate that it is a 3D emiss field
    emissions(ecount)%lowest_lev  = 1
    emissions(ecount)%highest_lev = model_levels
    emissions(ecount)%three_dim   = .TRUE.
  ELSE
    ! Report error just to double-check. However em_var_check_dims
    ! should have failed if the number of vertical levels is not
    ! equal to 1 or to model_levels
    ierror = fileid
    CALL ereport('UKCA_EM_STRUCT_SET', ierror,             &
      'Dimension mismatch for: ' // TRIM(variable_name))
  END IF
END IF  ! dimen_size

emissions(ecount)%values(:,:,:) = 0.0  ! will be filled when we
emissions(ecount)%diags (:,:,:) = 0.0  ! update emissions

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_em_struct_set


! ---------------------------------------------------------------------
! Description:
!   Check if it is time to update the values in the emissions structure,
!   and do it when needed.
! Method:
! 1) Check first if full hour. If so go through all emission fields and
!    decide to update their values only when the hours from start of the
!    run are multiple of their update frequency (in hours).
! 2) When it is time to update then take the current model time
!    (i_year, i_month, i_day, i_hour, i_minute, i_second) from model_time_mod,
!    and pass it to GET_EMFILE_REC to get the previous and next register.
!    Then call EM_GET_DATA to get the emission values for those two registers,
!    and interpolate the values if none of the registers coincides exactly
!    with the model time.
! ---------------------------------------------------------------------
SUBROUTINE ukca_emiss_update (row_length, rows, model_levels, &
                              global_row_length, global_rows, &
                              timestep_number, timestep,      &
                              l_first)

USE model_time_mod, ONLY: &
    i_day, i_hour, i_minute, i_month, i_second, i_year

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length       ! model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: global_row_length
INTEGER, INTENT(IN) :: global_rows
INTEGER, INTENT(IN) :: timestep_number  ! Nr of timesteps since start

REAL,    INTENT(IN) :: timestep         ! Timestep length (sec)

LOGICAL, INTENT(IN) :: l_first          ! True if first call

! Local variables
INTEGER :: ancil_ref_days, ancil_ref_secs  ! Ancil reference time in days and
                                           ! seconds since 0,0 basis time

INTEGER :: i_year_trg, i_month_trg, i_day_trg      ! Target interpolation time:
INTEGER :: i_hour_trg, i_minute_trg, i_second_trg  ! year, month, etc

INTEGER :: day, sec          ! Current time in days/seconds since ancil ref time
REAL    :: hours_from_ancref ! Hours from ancil ref time to current time
REAL    :: fhr               ! Hours from ancil ref time to mid-point of update
                             ! interval
INTEGER :: fdaynum           ! Day number returned by sec2time (not used)
INTEGER :: anc_steps         ! Number of ancil update intervals between ancil
                             ! ref time and current time

INTEGER             :: fid              ! NetCDF file id

INTEGER             :: ireco(2)      ! NetCDF time record indices: prev & next
INTEGER             :: creco(2,6)    ! Time records: 1st dim is for prev & next
                                     ! 2nd dim is for yr, mon, day, hr, min, sec

REAL                :: ftreco(2)     ! Time of prev/next time records (fract
                                     ! hours since 0:00 h on 1st Jan)

REAL                :: fhournow      ! Target time (fractional hours since
                                     ! 0:00 h on 1st Jan)

INTEGER, PARAMETER  :: max_dims = 5            ! Max nr of dims per var
INTEGER             :: dimsize (max_dims)      ! dimension sizes
REAL, ALLOCATABLE   :: tmpdat1 (:,:,:)         ! Temporal data arrays
REAL, ALLOCATABLE   :: tmpdat2 (:,:,:)
INTEGER             :: l                       ! Index
LOGICAL             :: l_noreco                ! True if record not found (EOF)
INTEGER             :: ierror

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_UPDATE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Reset l_update to false for all emission fields.
! Note that this is changed to TRUE below separately for each
! NetCDF emission field if it is updated in the current timestep.
! UKCA_NEW_EMISS_CTL also updates this field to TRUE for all online
! emissions (e.g. lightning NOx and wetland CH4) at each timestep.
emissions(:)%l_update = .FALSE.

! Get ancil reference time as days/secs since 0/0 basis time and use this to
! calculate hours from ancil reference time to now.
CALL time2sec(ancil_reftime(1), ancil_reftime(2), ancil_reftime(3), &
              ancil_reftime(4), ancil_reftime(5), ancil_reftime(6), &
              0,0,ancil_ref_days,ancil_ref_secs,lcal360)

CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,         &
              ancil_ref_days,ancil_ref_secs,day,sec ,lcal360)
hours_from_ancref = day*rhour_per_day + sec*rhour_per_sec


! Loop over NetCDF emission fields
DO l = 1, num_em_flds
! Online emissions are not in files so ignore in this loop
  IF (emissions(l)%l_online) CYCLE

  ! Decide if updating, and set target time for time-interpolation
  IF (emissions(l)%update_type == iupd_single) THEN
  ! Single-time emissions read on first timestep and no other
  ! Set target time to now for get_emfile_rec (not actually used)
    IF (l_first) THEN
      emissions(l)%l_update = .TRUE.
      i_year_trg = i_year
      i_month_trg = i_month
      i_day_trg = i_day
      i_hour_trg = i_hour
      i_minute_trg = i_minute
      i_second_trg = i_second
    END IF

  ELSE

    ! Period or timeseries emissions.
    ! Compute integer number of ancil steps from ref time to the begining of the
    ! present update period (=anc_steps). Use FLOOR to round down (even if 
    ! hours_from_ancref is negative) to get time which is at or before now.
    ! Force update when anc_steps increases from that stored at last update.
    ! Always update on first timestep.
    ! Set target time for time-interp to be mid-point of update interval:
    ! start of this update interval plus half the update frequency
    anc_steps = FLOOR(hours_from_ancref / emissions(l)%update_freq)

    IF (l_first .OR. anc_steps /= emissions(l)%last_update) THEN
      emissions(l)%l_update = .TRUE.
      emissions(l)%last_update = anc_steps

      fhr = (REAL(anc_steps) + 0.5) * REAL(emissions(l)%update_freq)
      day = INT(fhr/ rhour_per_day)
      sec = INT((fhr - day*ihour_per_day) * rsec_per_hour)
      CALL sec2time (day, sec, ancil_ref_days, ancil_ref_secs, &
        i_year_trg, i_month_trg, i_day_trg, &
        i_hour_trg, i_minute_trg, i_second_trg, &
        fdaynum, lcal360 )
      END IF
    END IF


    IF (emissions(l)%l_update) THEN

      !     Debug update times only for first emission field
      IF (l ==1 .AND. PrintStatus >= PrStatus_Diag) THEN
        WRITE (umMessage, '(A,6(1x,I4))')                                  &
         'UKCA_EMISS_UPDATE: Update first emission field for target time', &
         i_year_trg, i_month_trg, i_day_trg,                               &
         i_hour_trg, i_minute_trg, i_second_trg
        CALL umPrint(umMessage,src='ukca_emiss_init')
      END IF

      !     Open emission file to read an emision field
      CALL em_fopen(emissions(l)%file_name, fid)

      !     Pass current model time (i_year, ..., i_second) and other data
      !     to GET_EMFILE_REC in order to get:
      !     * current hour since basis time: fhournow,
      !     * indices of the matching NetCDF records used for
      !       interpolation: ireco (1:2)
      !     * yr, mon, day, hr, min, sec for those two records: creco (1:2,1:6)
      !     * hour since basis time for those records: ftreco (1:2)
      CALL get_emfile_rec (row_length, rows, fid, l,                    &
               emissions(l)%update_type,                                &
             (/i_year_trg, i_month_trg, i_day_trg,                      &
               i_hour_trg, i_minute_trg, i_second_trg /),               &
               l_first, ireco, creco, fhournow, ftreco)

      !     Later we will call EM_GET_DATA to get the emission values for registers
      !     ireco(1) & ireco(2), at hours ftreco(1) and ftreco(2), so that we can
      !     interpolate them to the current model time in hours (i.e. to hournow).

      !     Check that matching time records are found
      IF ( ANY(ireco(:) == -1) ) THEN  
        WRITE (umMessage,'(A,6(1x,I4),1x,F15.5)') 'At time:', i_year, &
          i_month, i_day, i_hour, i_minute, i_second, fhournow
        CALL umPrint(umMessage,src='ukca_emiss_update')
        CALL em_fclose(fid)
        ierror = l
        CALL ereport('UKCA_EMISS_UPDATE', ierror,                             &
          'Required time record not in file ' // emissions(l)%file_name)
      END IF

      !     Go ahead with reading values and interpolations

      !     First get dimensions, in particular the vertical levels, so
      !     that we can allocate the temp data arrays tmpdat1 & tmpdat2.
      CALL em_var_check_dims (global_row_length, global_rows, model_levels,   &
                              fid, emissions(l)%var_name, dimsize)

      l_noreco = .FALSE.  ! initially assume that records are found

      IF ( ireco(1) == ireco(2) ) THEN
        !     No interpolation if exact time match. Get only emiss values
        !     for register ireco(1).
        CALL em_get_data (row_length, rows, fid, ireco(1),          &
                          emissions(l)%var_name,                    &
                          emissions(l)%values (:, :, 1:dimsize(3)), &
                          l_noreco )
        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          CALL ereport('UKCA_EMISS_UPDATE',ierror,                       &
            'Error reading data1: ' // emissions(l)%var_name)
        END IF

      ELSE
        !       Get two emission values tmpdat1 and tmpdat2, for registers
        !       ireco(1) and ireco(2), to later interpolate.
        ALLOCATE ( tmpdat1 (row_length, rows, dimsize(3)) )
        ALLOCATE ( tmpdat2 (row_length, rows, dimsize(3)) )

        !       Read in data for previous matching time. Stop if record not found.
        CALL em_get_data (row_length, rows, fid, ireco(1),          &
                          emissions(l)%var_name,                    &
                          tmpdat1 (:,:,:), l_noreco)
        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          CALL ereport('UKCA_EMISS_UPDATE',ierror,                        &
             'Error reading data1: ' // emissions(l)%var_name)
        END IF

        !       Read in data for next matching time.
        CALL em_get_data (row_length, rows, fid, ireco(2),          &
                          emissions(l)%var_name,                    &
                          tmpdat2 (:,:,:), l_noreco)
        IF (l_noreco) THEN  
          CALL em_fclose(fid)
          ierror = l
          CALL ereport('UKCA_EMISS_UPDATE',ierror,                        &
             'Error reading data2: ' // emissions(l)%var_name)
        END IF

        ! Interpolate using same UM routine as REPLANCA
        CALL t_int (                                  &
             tmpdat1 (:,:,:), ftreco(1),              &  ! Data1 & time1
             tmpdat2 (:,:,:), ftreco(2),              &  ! Data2 & time2
             emissions(l)%values(:, :, 1:dimsize(3)), &  ! Data interpolated
             fhournow,                                &  ! to current time
             row_length*rows*dimsize(3))                 ! Nr of points

        DEALLOCATE(tmpdat2)
        DEALLOCATE(tmpdat1)
      END IF      ! Interpolation required?

      !     Close emission file after reading an emision field
      CALL em_fclose(fid)

  END IF        ! Update time?
END DO          ! loop over emission species

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_update


!---------------------------------------------------------------------------
! Description:
!   Convert total aerosol emission mass into mass and number for each mode.
!
! Method:
!   For each mode and component this emission field maps to, compute number and
!   mass, and add to appropriate emissions structures.
!   This routine doesn't care whether the emissions are 2d or 3d: it just
!   converts the source mass into number and mass for each mode. The (level and
!   time) attributes of the modal emissions are defined in ukca_emiss_init.
!   The source emission is normally read from a netcdf file (offline emission),
!   but could be calculated online by the model.
!---------------------------------------------------------------------------
SUBROUTINE ukca_emiss_update_mode(row_length, rows, model_levels, &
                                  tracer_name, three_dim, emiss_values, &
                                  emmas, emnum, lmode_emiss, &
                                  nlev, icp)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length    ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

CHARACTER (LEN=10), INTENT(IN) :: tracer_name  ! name of emission tracer
LOGICAL, INTENT(IN) :: three_dim  ! source emissions are 3D

! Source emission mass (kg(component)/m2/s) which needs converting to modal
! emissions. Assumed-shape array because 3rd dim can be 1 or model_levels.
REAL, INTENT(IN) :: emiss_values(:,:,:)

! 3D number (equivalent kg/m2/s) and mass (kg(component)/m2/s) emissions for
! each mode.
REAL, INTENT(OUT) :: emnum(row_length, rows, model_levels, nmodes)
REAL, INTENT(OUT) :: emmas(row_length, rows, model_levels, nmodes)

! Which modes receive some emission for this species (TRUE for those that do)
LOGICAL, INTENT(OUT) :: lmode_emiss(nmodes)

INTEGER, INTENT(OUT) :: nlev  ! Number of levels emission is applied to
INTEGER, INTENT(OUT) :: icp   ! Aerosol component emission is applied to

! Mass (kg/m2/s) and particle volume (nm3/m2/s) emissions for an individual mode
REAL :: modemass(row_length, rows, model_levels)
REAL :: modevol(row_length, rows, model_levels)

INTEGER :: imode            ! Loop index for mode number.
REAL :: mode_frac(nmodes)   ! Fraction of emission in each mode (0->1).
REAL :: mode_diam(nmodes)   ! Geometric mean diameter for each mode emiss (nm)
REAL :: mode_stdev(nmodes)  ! Geometric stdev for each mode emiss (units=1)

REAL    :: mm_da  !  Molar mass of dry air (kg/mol) =avc*zboltz/ra
REAL    :: factor  !  converts from partcls per m2/s to kg/m2/s
REAL    :: lstdev  !  log(stdev)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_UPDATE_MODE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Call ukca_def_mode_emiss to obtain emission mask, mode fraction, diameter
! and stdev for each mode, as well as component emissions applied to (=icp).
CALL ukca_def_mode_emiss(tracer_name, lmode_emiss, icp, &
                         mode_frac, mode_diam, mode_stdev)

! emiss%values might only have 1 level, so ensure consistency when assiging to
! modemass. All other arrays used below have all levels, so don't need to use 
! nlev after that.
IF (three_dim) THEN
  nlev = model_levels
ELSE
  nlev = 1
END IF

mm_da = avc*zboltz/ra  ! Molar mass of dry air (kg/mol)
factor = mm_da / avc   ! Converts from partcls per m2/s to kg/m2/s
modemass(:,:,:) = 0.0
modevol(:,:,:) = 0.0
emnum(:,:,:,:) = 0.0
emmas(:,:,:,:) = 0.0

! For each mode, calculate mass, volume and hence number emissions and save into
! emmas and emnum arrays.
DO imode = 1, nmodes
  IF (lmode_emiss(imode)) THEN
    lstdev = LOG(mode_stdev(imode))

    ! Particulate mass emissions (kg of component /m2/s) in mode.
    ! Note that units conversions (C to OC, S to SO2) are done in 
    ! ukca_new_emiss_ctl using base_scaling, and scaling for particulate 
    ! fraction and conversion from SO2 to H2SO4 is done in ukca_add_emiss using 
    ! aero_scaling.
    modemass(:,:,1:nlev) = emiss_values(:,:,1:nlev) * mode_frac(imode)
    emmas(:,:,:,imode) = modemass(:,:,:)

    ! Total particle volume emission rate (nm3 per m2 per s) in this mode is 
    ! the rate of mass emission (kg/m2/s) divided by the density of this mode
    ! (kg/m3) multiplied by 1e27 (1e9^3) to go from m3/m2/s to nm3/m2/s.
    modevol(:,:,:) = 1e27*modemass(:,:,:)/rhocomp(icp)

    ! Total particle number (particles per m2 per s) * factor.
    ! factor converts from ptcls/m2/s to equivalent kg/m2/s as dry air.
    emnum(:,:,:,imode) = modevol(:,:,:) * factor /                           &
      ( (pi/6.0) * (mode_diam(imode)**3) * EXP(4.5*lstdev*lstdev) )

  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_update_mode


!---------------------------------------------------------------------------
! Description:
!   Map number and mass emissions arrays to appropriate mode emissions
!   structures
!    
! Method:
!   For a given emission field (from file or online), take the provided mass and
!   number and find the emissions structures with the corresponding moment,
!   mode, component and source species. Copy the emissions into 
!   emissions(l)%values.
!---------------------------------------------------------------------------
SUBROUTINE ukca_emiss_mode_map(row_length, rows, model_levels, emmas, emnum, &
                               emiss_levs, icp, lmode_emiss, from_emiss)

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length    ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! 3D number (equivalent kg/m2/s) and mass (kg(component)/m2/s) emissions for
! each mode.
REAL, INTENT(IN) :: emmas(row_length, rows, model_levels, nmodes)
REAL, INTENT(IN) :: emnum(row_length, rows, model_levels, nmodes)

! Number levels of emission field. To ensure that emnum/emmas slice is correct
! shape for emissions(l)%values.
INTEGER, INTENT(IN) :: emiss_levs

INTEGER, INTENT(IN) :: icp         ! Aerosol component for this emission
INTEGER, INTENT(IN) :: from_emiss  ! Index of source emission
LOGICAL, INTENT(IN) :: lmode_emiss(nmodes)   ! Which modes receive emission?

! Local variables
INTEGER :: nlev            ! Short name from emiss_levs argument
INTEGER :: imode, imoment  ! Loop indices for mode and moment number
INTEGER :: l               ! Loop index for emissions structure array
LOGICAL :: lfound          ! Flag for whether matching emission found

CHARACTER (LEN=errormessagelength) :: cmessage    ! Error message
INTEGER                            :: ierror      ! Error code

nlev = emiss_levs

! For each mode and for both number and mass, loop through all emissions
! structures until we find a mode emission with matching mode, component,
! moment and name of source tracer. Copy into the %values tag of the matching
! structure.
DO imode = 1, nmodes
  IF (lmode_emiss(imode)) THEN
    DO imoment = 0, 3, 3  ! imoment = 0 and 3 only
      lfound = .FALSE.
      DO l = 1, num_em_flds
        ! Find the corresponding mode emissions structures for mass and number
        IF (emissions(l)%l_mode .AND. &
            emissions(l)%mode == imode .AND. &
            emissions(l)%component  == icp .AND. &
            emissions(l)%moment == imoment .AND. &
            emissions(l)%from_emiss == from_emiss) THEN
          SELECT CASE(imoment)
          CASE (0) ! number
            emissions(l)%values(:,:,1:nlev) = emnum(:,:,1:nlev,imode)
          CASE (3) ! mass
            emissions(l)%values(:,:,1:nlev) = emmas(:,:,1:nlev,imode)
          END SELECT

          lfound = .TRUE.
          EXIT
        END IF  ! mode==imode and component==icp
      END DO  ! num_em_flds
    END DO  ! imoment

    IF (.NOT. lfound) THEN
      ! If this happens something has gone wrong in ukca_emiss_init or
      ! ukca_def_mode_emiss
      WRITE(cmessage,'(2A,3I3)') 'Missing mode emission structure: ',  &
        from_emiss, imode, icp, imoment
      ierror = icp
      CALL ereport('UKCA_EMISS_MODE_MAP', ierror, cmessage)
    END IF
  END IF  ! lmode_emiss(imode) and component(imode,icp)
END DO  ! nmodes

RETURN
END SUBROUTINE ukca_emiss_mode_map


! ---------------------------------------------------------------------
! Description:
!   Get time records (indices as well as values) from a NetCDF
!   emission or oxidant file that correspond to the closest times
!   before and after the requested target_time (usually the mid-point of an
!   ancillary update period)
!
! Method:
! 1) Call time2sec to get  target time in hours and seconds.
! 2) Loop over time records in the file, calling EM_GET_TIME_REC
!    to get time for the previous and next record.
!    When fhr_now is between both records then get the indices and
!    current times (both as a date and in hours since base time)
!    of those two records  and exit the loop.
! ---------------------------------------------------------------------
SUBROUTINE get_emfile_rec ( row_length, rows, fid, nfield,             &
                            iupd_type, target_time, L_first,           &
                            irec, crec, fhr_now, ftrec )

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: row_length     ! model horiz dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: fid            ! NetCDF file ID (already open)
INTEGER, INTENT(IN) :: nfield         ! nr of NetCDF emiss field
INTEGER, INTENT(IN) :: iupd_type      ! iupd_serial=1, iupd_cyclic=2
INTEGER, INTENT(IN) :: target_time(6) ! Target time:
                                      ! yr, mon, day, hr, min, sec

LOGICAL, INTENT(IN) :: L_first

INTEGER, INTENT(OUT):: irec(2)        ! Time record indices: prev & next
INTEGER, INTENT(OUT):: crec(2,6)      ! Readable time of records
                                      ! 1st dim is for prev & next
                                      ! 2nd dim for yr, mn, day, hr, mon, sec

REAL,    INTENT(OUT):: fhr_now        ! Target time in fract hours
                                      ! since base time (0:00 h on 1st Jan).
                                      ! Here 'now' refers to target time not
                                      ! present model time.
REAL,    INTENT(OUT):: ftrec(2)       ! Prev/next time in hours since
                                      ! the same base time

! Local variables
INTEGER          :: daynow, secnow    ! target time in days and secs
                                      ! elapsed since 0:00 on 1st Jan

INTEGER          :: tday1, tday2, tsec1, tsec2, fdaynum ! prev/next times
INTEGER          :: exact_match       ! Null if interpol required, and
                                      ! 1 or 2 if a NetCDF reg is matched

INTEGER          :: target_yr             ! target year
INTEGER          :: itime, itime1, itime2 ! time rec counters for NetCDF files

REAL               :: ftime(2)              ! time record value in file
REAL, PARAMETER    :: small_number = 0.00001   ! For comparing reals numbers

LOGICAL          :: L_notime         ! True if end-of-file
LOGICAL          :: L_new_cdf_yr     ! True if periodic data and need to read
                                     ! Dec of year i and Jan of year i+1
LOGICAL          :: L_new_target_yr  ! True if target year is i+1

INTEGER          :: ierror

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_EMFILE_REC'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Start with negat indices for prev/next time records
! until they are found in the TIME_REC loop
irec(:) = -1

! Time in NetCDF files is in days since the base time (0:00 h on 1st Jan
! of a given year). The code below gets target days/secs since
! the same base time so that we can interpolate.
! When file is periodic then ignore year. Set to cyclic_year to ensure
! consistency with em_get_time_rec.
IF (iupd_type == iupd_cyclic) THEN  ! periodic time series
  target_yr = cyclic_year
ELSE
  target_yr = target_time(1)
END IF

! Get current time in days and seconds (daynow & secnow)
! elapsed since base time (0:00 h on 1st Jan)
! DEPENDS ON: time2sec
CALL time2sec ( target_yr,      target_time(2),  &  ! current time: yr,  mon,
                target_time(3), target_time(4),  &  !               day, hour
                target_time(5), target_time(6),  &  !               min, sec
                0, 0,                          &  ! basis time: 0 h on 1st Jan
                daynow, secnow, lcal360 )

! Note that the code below interpolates to the start time of the
! given timestep (integer hours) plus half the update frequency.
!
! Calculate target time as hours from ref time
fhr_now =  daynow * ihour_per_day + secnow * rhour_per_sec

! If extra diags messg then print time info only for first NetCDF field
IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
  WRITE (umMessage,'(A,I4,5(1x,I2),1x,I6,1x,I6,1x,F12.3)')       &
    'GET_EMFILE_REC: Get data field for ',                       &
     target_time(1:6), daynow, secnow, fhr_now
  CALL umPrint(umMessage,src='ukca_emiss_init')
END IF

! Loop over time records in file until we find two consecutive records
! before and after the target time. Note that in the case of cyclic
! emissions we assume 12 NetCDF registers with monthly means and
! some conditions will be needed for particular cases when the target
! time should fall between register 12 (December) and register 1
! (January) of the NetCDF files.
itime          = 1         ! Initialise time counter
L_notime       = .FALSE.   ! Initially assume that time records found

time_rec: DO    ! loop over time records

  ! Initially assume that no need to change year either in
  ! the target or in the NetCDF files
  L_new_target_yr = .FALSE.
  L_new_cdf_yr   = .FALSE.

  ! Get previous time (from a register in NetCDF file)
  ! in days and secs: tday1, tsec1
  ! This first register will usually be itime1 = itime
  !
  ! However when emissions are periodic/cyclic and the target time is in
  ! early January that would be before the time of the first register in the
  ! NetCDF files (mid-January assuming 12 registers with monthly emissions).
  ! In that particular case no previous register will be found within this loop
  ! and itime will reach 13, which will cause an end-of-file. For this
  ! particular case wee need the condition below to indicate that:
  ! * 1st NetCDF reg. is in Dec of year i   (itime1=12, L_new_cdf_yr = .FALSE.)
  ! * 2nd NetCDF reg. is in Jan of year i+1 (itime2= 1, L_new_cdf_yr = .TRUE.)
  ! * target time is in year i+1 (L_new_target_yr = .TRUE.). That will be
  !   between both NetCDF regs and interpolations should work.
  IF (iupd_type == iupd_cyclic .AND. itime == 13) THEN
    ! 1st NetCDF register in Dec of year i
    itime1     = 12
    CALL em_get_time_rec (row_length, rows, fid, itime1, iupd_type,     &
                          L_first, L_new_cdf_yr, tday1, tsec1, L_notime)

    ! 2nd NetCDF register in Jan of year i+1
    itime2         = 1
    L_new_cdf_yr   = .TRUE.

    ! Recalculate the target time 'fhr_now' as done before this DO loop
    ! but indicating that year is i+1 so that interpolations can be made.
    L_new_target_yr = .TRUE.
    ! DEPENDS ON: time2sec
    CALL time2sec ( target_yr + 1,  target_time(2),  &
                    target_time(3), target_time(4),  &
                    target_time(5), target_time(6),  &
                    0, 0,                          &
                    daynow, secnow, lcal360 )

    fhr_now = daynow * ihour_per_day + secnow * rhour_per_sec

    IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
      WRITE (umMessage,'(A,I4,5(1x,I2),1x,I6,1x,I6,1x,F12.3,1x,A)')   &
        'GET_EMFILE_REC: Get emiss for ',                     &
        target_time(1:6), daynow, secnow, fhr_now,             &
        '(after shifting target time one year forward)'
      CALL umPrint(umMessage,src='ukca_emiss_init')
    END IF

  ELSE
    ! Simple case where the max nr of NetCDF registers has not been exceeded
    itime1     = itime
    CALL em_get_time_rec (row_length, rows, fid, itime1, iupd_type,     &
                          L_first, L_new_cdf_yr, tday1, tsec1, L_notime)
  END IF

  ! The line below checks that the file has not run out of time records.
  ! In the future the code might allow to automatically open the
  ! next file using prefix+timestamp
  IF ( L_notime ) THEN
       ierror = itime1 
       CALL ereport ('GET_EMFILE_REC', ierror,              &
      'Invalid time record 1, End-of-file reached ?')
  END IF

  IF (iupd_type == iupd_single) THEN
    ! Single time in file - get the first time record and return.
    ftrec(1) = REAL(tday1) * rhour_per_day + REAL(tsec1) * rhour_per_sec
    irec(1) = itime1

    ! Get time yr-mon-day-hr-min-sec and store them in crec (1,1:6) for
    ! prev time record (given as tday1-tsec1). Useful for debugging.
    ! DEPENDS ON: sec2time
    CALL sec2time (tday1, tsec1, 0, 0,                                     &
         crec(1,1), crec(1,2), crec(1,3), crec(1,4), crec(1,5), crec(1,6), &
         fdaynum, lcal360 )

    irec  (2)   = irec  (1)       ! record indices
    crec  (2,:) = crec  (1,:)     ! record times: yr, mon, day, hr, min, sec
    ftrec (2)   = ftrec (1)       ! hours since base time

    EXIT  time_rec          ! exit loop over time records
  END IF

  ! Get next time in days and secs: tday2, tsec2
  ! This second register will usually be itime2 = itime1 + 1.
  !
  ! However if emissions are periodic and the first NetCDF register is in
  ! December (12) then the second reg to read is in January (1) of the
  ! following year. L_new_cdf_yr is set to TRUE to inform EM_GET_TIME_REC
  ! that the year should be increased one unit for that second register.
  !
  ! We exclude checks for the case when the target year has also changed
  ! since that has been dealt with above
  IF (.NOT. L_new_target_yr) THEN
    IF (iupd_type == iupd_cyclic .AND. itime1 == 12) THEN
      itime2       = 1
      L_new_cdf_yr = .TRUE.
    ELSE
      itime2       = itime1 + 1
    END IF
  END IF
  CALL em_get_time_rec (row_length, rows, fid, itime2, iupd_type,      &
                        L_first, L_new_cdf_yr, tday2, tsec2, L_notime)

  IF ( L_notime ) THEN
      ierror = itime2
       CALL ereport ('GET_EMFILE_REC', ierror,              &
      'Invalid time record 2, End-of-file reached ?')
  END IF

  ! Convert record times (in days and secs) to fractional hours
  ! to compare with fhr_now
  ftrec(1) = REAL(tday1)*rhour_per_day + REAL(tsec1)*rhour_per_sec ! prev rec
  ftrec(2) = REAL(tday2)*rhour_per_day + REAL(tsec2)*rhour_per_sec ! next rec

  ! Check if current time is between two time records
  IF ( ftrec(1) <= fhr_now .AND. fhr_now <= ftrec(2) ) THEN  

    irec(1) = itime1             ! Prev time record for interpolation
    irec(2) = itime2             ! Next time record for interpolation

    ! Get time yr-mon-day-hr-min-sec and store them in crec (1,1:6) for
    ! prev time record (given as tday1-tsec1). Useful for debugging.
    ! DEPENDS ON: sec2time
    CALL sec2time (tday1, tsec1, 0, 0,                                     &
         crec(1,1), crec(1,2), crec(1,3), crec(1,4), crec(1,5), crec(1,6), &
         fdaynum, lcal360 )

    ! Similarly, get time yr-mon-day-hr-min-sec for next time record
    ! DEPENDS ON: sec2time
    CALL sec2time (tday2, tsec2, 0, 0,                                     &
         crec(2,1), crec(2,2), crec(2,3), crec(2,4), crec(2,5), crec(2,6), &
         fdaynum, lcal360 )

    ! Assume that not exact match, i.e. interpolation required
    exact_match = 0

    ! If exact time match, no interpolation required. Therefore
    ! set both recs to same value
    IF ( ABS (ftrec(1) - fhr_now) < small_number ) THEN
      irec  (2)   = irec  (1)       ! record indices
      crec  (2,:) = crec  (1,:)     ! record times: yr, mon, day, hr, min, sec
      ftrec (2)   = ftrec (1)       ! hours since base time
      exact_match = 1
    END IF

    IF ( ABS (ftrec(2) - fhr_now) < small_number ) THEN
      irec  (1)   = irec  (2)
      crec  (1,:) = crec  (2,:)
      ftrec (1)   = ftrec (2)
      exact_match = 2
    END IF

    ! If extra diag mssg then print time info for first NetCDF field
    IF (PrintStatus >= PrStatus_Diag .AND. nfield == 1) THEN
      SELECT CASE  (exact_match)
      CASE (0)
        WRITE (umMessage, '(A,1X,I5,1X,I5,1X,F12.3,1X,F12.3)')               &
          "GET_EMFILE_REC: Time found between these regs. in NetCDF file:",  &
          irec(1), irec(2), ftrec(1), ftrec(2)
        CALL umPrint(umMessage,src='ukca_emiss_init')
      CASE (1)
        WRITE (umMessage, '(A,1X,I5,1X,F12.3)')                              &
          "GET_EMFILE_REC: Time found for this reg. in NetCDF file:",        &
          irec(1), ftrec(1)
        CALL umPrint(umMessage,src='ukca_emiss_init')
      CASE (2)
        WRITE (umMessage, '(A,1X,I5,1X,F12.3)')                              &
          "GET_EMFILE_REC: Time found for this reg. in NetCDF file:",        &
          irec(2), ftrec(2)
        CALL umPrint(umMessage,src='ukca_emiss_init')
      END SELECT
    END IF

    EXIT  time_rec          ! exit loop over time records

  ELSE
    itime = itime + 1       ! increment time rec counter => next rec in file

  END IF                    ! if t1 <= current <= t2

END DO  time_rec            ! loop over time records

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE get_emfile_rec
! ---------------------------------------------------------------------

END MODULE ukca_emiss_mod
