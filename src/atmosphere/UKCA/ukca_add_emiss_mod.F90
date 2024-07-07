! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module containing the subroutine UKCA_ADD_EMISS to
!    add emission fields to UKCA tracers and do boundary
!    layer mixing of tracers in the new UKCA emission system.
!    Called from the top level emission routine UKCA_NEW_EMISS_CTL
!    at each timestep.
!
!  Method:
!  1) Go through all emission fields in the emissions structure and
!     apply temporal/vertical factors.
!  2) The corrected emission values are stored as diagnostics and
!     added to the position of the array em_field that correspond
!     to the emitted tracer.
!  3) Loop through all tracers: call TRSRCE to add em_field to the
!     corresponding tracer (for all model levels except for surface),
!     and call TR_MIX to do tracer mixing and add surface emissions.
!  4) Finally ukca_volcanic_so2 is called for SO2 emissions from explosive
!     volcanic eruptions (only relevant for the stratosphere)
!  5) The array with updated tracers is sent back to UKCA_MAIN1
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
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
MODULE ukca_add_emiss_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_ADD_EMISS_MOD'

CONTAINS

SUBROUTINE ukca_add_emiss (                                           &
    row_length, rows, bl_levels,                                      &
    n_tracers, iyear, imonth, iday, ihour, timestep,                  &
    delta_lambda,       delta_phi,                                    &
    sin_theta_latitude, true_longitude, tropopause_height,            &
    r_theta_levels,     theta,          q,  qcl,  qcf,                &
    exner_rho_levels,   rho_r2,                                       &
    kent,               kent_dsc,                                     &
    rhokh_mix,          dtrdz_charney_grid,                           &
    we_lim,             t_frac,                 zrzi,                 &
    we_lim_dsc,         t_frac_dsc,             zrzi_dsc,             &
    ml_depth,           zhsc,                   z_half,               &
    surf_area,          totnodens,  volume,  mass,  tracers,          &
    len_stashwork50,stashwork50)

USE bl_option_mod,           ONLY: alpha_cd
USE tr_mix_mod,              ONLY: tr_mix
USE trsrce_mod,              ONLY: trsrce
USE asad_mod,                ONLY: advt, method
USE ukca_chem_schemes_mod,   ONLY: int_method_nr
USE asad_chem_flux_diags
USE conversions_mod, ONLY: pi
USE ukca_constants
USE ukca_trace_gas_mixratio, ONLY: um_ch4_for_ukca
USE dust_parameters_mod,     ONLY: ndiv
USE ukca_option_mod,         ONLY: L_ukca_stratcfc, L_ukca_strat,     &
                                   L_ukca_prescribech4,               &
                                   L_ukca_strattrop,                  &
                                   L_ukca_achem, L_ukca_aerchem,      &
                                   l_ukca_chem,  jpctr,               &
                                   l_ukca_offline, l_ukca_offline_be, &
                                   l_ukca_nr_aqchem,                  &
                                   l_ukca_scale_biom_aer_ems,         &
                                   biom_aer_ems_scaling, mode_parfrac,&
                                   l_ukca_so2ems_expvolc
USE ukca_d1_defs,            ONLY: n_use_emissions, em_chem_spec,     &
                                   n_boundary_vals, lbc_mmr, lbc_spec
USE get_molmass_mod,         ONLY: get_molmass
USE ukca_mode_setup,         ONLY: nmodes, ncp, mode, component, cp_su, mm
USE ukca_mode_tracer_maps_mod, ONLY:  nmr_index, mmr_index
USE ukca_setup_indices,      ONLY: msotwo, mm_gas

! Use data from emissions structure in new emiss system
USE ukca_emiss_mod,          ONLY: num_em_flds, emissions
USE ukca_emiss_factors
USE ukca_emiss_mode_mod,     ONLY: aero_ems_species

USE nlstcall_mod,            ONLY: lcal360
USE day_of_week_mod,         ONLY: day_of_week    ! Calculate day of week
USE parkind1,                ONLY: jprb, jpim     ! DrHook
USE yomhook,                 ONLY: lhook, dr_hook ! DrHook
USE ereport_mod,             ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE UM_ParVars
USE Control_Max_Sizes
USE nlsizes_namelist_mod, ONLY:model_levels
USE umPrintMgr, ONLY: umPrint, umMessage, PrintStatus, PrStatus_Normal

! Still need to decide what to do with some emissions when
! the new UKCA emission system is working with MODE
! USE run_aerosol_mod,         ONLY: SO2_high_level

USE ukca_volcanic_so2_mod, ONLY: ukca_volcanic_so2

USE stash_array_mod,    ONLY: sf, si

IMPLICIT NONE

! Subroutine arguments

! Input arguments with info on model dimensions
INTEGER, INTENT(IN)            :: row_length
INTEGER, INTENT(IN)            :: rows
INTEGER, INTENT(IN)            :: bl_levels

! Input arguments to get nr tracers and current model time
INTEGER, INTENT(IN) :: n_tracers   ! nr traces: chem + mode
INTEGER, INTENT(IN) :: iyear       ! model yr, mo, day, hr
INTEGER, INTENT(IN) :: imonth
INTEGER, INTENT(IN) :: iday
INTEGER, INTENT(IN) :: ihour
REAL,    INTENT(IN) :: timestep    ! timestep lenght (sec)

! Input arguments neded to get SO2 emiss from explosive volcanic eruptions
REAL, INTENT(IN) :: delta_lambda
REAL, INTENT(IN) :: delta_phi
REAL, INTENT(IN) :: sin_theta_latitude    (1:row_length,1:rows)
REAL, INTENT(IN) :: true_longitude        (1:row_length,1:rows)
REAL, INTENT(IN) :: tropopause_height     (1:row_length,1:rows)
REAL, INTENT(IN) :: r_theta_levels        (1:row_length,1:rows,0:model_levels)

! Input arguments needed to call TSRCE
REAL, INTENT(IN) :: theta             (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: q                 (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcl               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: qcf               (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: exner_rho_levels  (1:row_length, 1:rows, 1:model_levels+1)
REAL, INTENT(IN) :: rho_r2            (1:row_length, 1:rows, 1:model_levels)

! Input arguments needed to call TR_MIX
INTEGER, INTENT(IN) :: kent       (1:row_length, 1:rows)
INTEGER, INTENT(IN) :: kent_dsc   (1:row_length, 1:rows)

REAL, INTENT(IN) :: rhokh_mix          (1:row_length, 1:rows, 1:bl_levels)
REAL, INTENT(IN) :: dtrdz_charney_grid (1:row_length, 1:rows, 1:bl_levels)

REAL, INTENT(IN) :: we_lim     (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: t_frac     (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: zrzi       (1:row_length, 1:rows, 1:3)

REAL, INTENT(IN) :: we_lim_dsc (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: t_frac_dsc (1:row_length, 1:rows, 1:3)
REAL, INTENT(IN) :: zrzi_dsc   (1:row_length, 1:rows, 1:3)

REAL, INTENT(IN) :: ml_depth   (1:row_length, 1:rows)
REAL, INTENT(IN) :: zhsc       (1:row_length, 1:rows)
REAL, INTENT(IN) :: z_half     (1:row_length, 1:rows,1:bl_levels)

! Input arguments (surf area, dens, vol and mass) mainly for ASAD emiss diags
REAL, INTENT(IN) :: surf_area  (1:row_length, 1:rows)
REAL, INTENT(IN) :: totnodens  (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: volume     (1:row_length, 1:rows, 1:model_levels)
REAL, INTENT(IN) :: mass       (1:row_length, 1:rows, 1:model_levels)

INTEGER, INTENT(IN) :: len_stashwork50
REAL, INTENT(INOUT) :: stashwork50(len_stashwork50)

! Tracer mass mixing ratios
REAL, INTENT(INOUT) :: tracers (1:row_length, 1:rows, 1:model_levels,   &
                                1:n_tracers)

! Local variables
INTEGER, PARAMETER  :: surface_level = 1

INTEGER, SAVE :: inox                 ! Index for NO/NOx tracer
INTEGER, SAVE :: iso2                 ! Index for SO2 tracer
INTEGER :: i, j, k, l, m, ilev, p, pp ! Loop variables
INTEGER :: imode, icp                 ! Loop variables
INTEGER :: n                          ! counter
INTEGER :: errcode                    ! Error code for ereport
INTEGER :: ierr                       ! Error code from asad diags routines

INTEGER :: iday_week            ! Day of week (for time profiling of emiss
                                ! in regional air quality modelling)
INTEGER :: t_local              ! Local time based on longitude
INTEGER, SAVE:: num_hourly_profs  ! Number of different hourly profiles
                                  ! of emissions used

REAL :: res_factor       (row_length, rows)              ! dry dep
REAL :: tr_flux          (row_length, rows, bl_levels, n_tracers)
                                                         ! tracer flux
REAL :: surf_dep_flux    (row_length, rows, n_tracers)
REAL :: em_field         (row_length, rows, model_levels, n_tracers)
                                                   ! 3D emiss for all tracers

REAL :: em_field_2d      (row_length, rows)    ! 2D emiss field for a tracer
                                               ! at a given model level

REAL :: rho_aresist      (row_length, rows)

REAL, SAVE, ALLOCATABLE :: theta_latitude (:,:)
REAL, SAVE, ALLOCATABLE :: molmass        (:)    ! molar masses
REAL, SAVE, ALLOCATABLE :: lbc_molmass    (:)    ! molar masses for lower
                                                 ! boundary conditions

! Scalar to store daily scaling accounting for weekly cycle of emissions
REAL :: fdaily_scaling

! Scalar to store scaling to remove particulate fraction of SO2 (default=1.0)
REAL :: aero_scaling

! 1-D arrays to store daily and hourly scaling factors of emissions
REAL :: daily_scaling  (7)             ! Sunday=1, ..., Saturday=7
REAL :: hourly_scaling (0:23)          ! from 0 to 23 h

! 2-D array to hold the hourly scaling as a function of the longitude
! (used as an approximation considering a local time based on the longitude)
REAL :: hourly_scaling_2d  (row_length, rows)

! Arrays to hold 1-D (function of local time) and 2-D (function of lon-lat)
! hourly scaling for a number of emission profiles. The last dimension
! is for the profile number.
REAL, SAVE, ALLOCATABLE :: hourly_scaling_all    (:,:)    ! 1D + profile
REAL, SAVE, ALLOCATABLE :: hourly_scaling_all_2d (:,:,:)  ! 2D + profile

! Character variables that will be passed to HOURLY_EMISS_FACTORS,
! i.e. to the subroutine dealing with hourly profiles
CHARACTER (LEN=20), SAVE, ALLOCATABLE :: hourly_prof_name (:)
CHARACTER (LEN=80), SAVE, ALLOCATABLE :: field_name       (:)

! Required for ASAD 3D emiss diags (this replicates code in UKCA_EMISSION_CTL).
REAL  :: conv_aircraftems (row_length,rows,model_levels) ! aircraft emiss
REAL, ALLOCATABLE :: tmp3dems (:,:,:)      ! explosive volcanic SO2 emiss
! Note that NOx lightning diags are in UKCA_NEW_EMISS_CTL and that additional
! diagnostics (e.g. other 3D SO2 emiss) can be included here in the future.

LOGICAL, SAVE :: firstcall = .TRUE.

CHARACTER (LEN=errormessagelength) :: cmessage  ! error message
CHARACTER (LEN=10)  :: mapped_tracer
CHARACTER (LEN=10), SAVE, ALLOCATABLE  :: nm_tracer(:)   
                    ! temporarily holds tracer name from advt
REAL, PARAMETER     :: small_number  = 0.00001 ! to compare real numbers

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_ADD_EMISS'

INTEGER                            :: section   ! stash section
INTEGER                            :: item      ! stash item
INTEGER                            :: icode = 0 ! local error status

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise variables
rho_aresist    (:,:)     = 0.0     ! Set dry deposition
res_factor     (:,:)     = 0.0     ! to zero

em_field       (:,:,:,:) = 0.0     ! 3D Initial ems field for all tracers

em_field_2d    (:,:)     = 0.0     ! 2D emission field for a tracer
                                   ! at a given model level

hourly_scaling_2d (:,:)  = 0.0     ! 2D hourly scaling as a function of
                                   ! lon/lat (related to local time)
                                   

! Get day of week (may be needed for time profiling of emissions in air
! quality modelling). Note that this would not have any meaning and
! therefore is not allowed in the case of a 360-day calendar.
IF (.NOT. lcal360) THEN
  iday_week = day_of_week (iday, imonth, iyear)
END IF

! Set some values in the first call to this routine
IF (firstcall) THEN
  inox             = -99      ! Initial index for nox tracer
  iso2             = -99      ! Initial index for so2 tracer

  ALLOCATE (theta_latitude (row_length, rows))
  theta_latitude = ASIN (sin_theta_latitude)

  !   Find index for SO2 and NOx tracer
  DO k = 1, jpctr
    SELECT CASE (advt(k))
    CASE ('NOx       ')
      inox = k
    CASE ('NO        ')
      inox = k
    CASE ('SO2       ')
      iso2 = k
    END SELECT
  END DO

  IF (inox == -99 .AND. .NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
    cmessage = 'Did not find NO or NOx tracer'
    errcode  = 1
    CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
  END IF

  IF ((L_ukca_aerchem .OR. L_ukca_nr_aqchem .OR. l_ukca_offline_be)            &
      .AND. iso2 == -99) THEN
    cmessage = 'Did not find SO2 tracer'
    errcode  = 1
    CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
  END IF

  ! Temporarily hold tracer names from advt (to check for 'H2O ' before call
  ! to TR_MIX in the loop below. Advt is dimensioned to 1:jpctr, while the 
  ! loop is over chem+aero_tracers, so fill rest of array with dummy values
  ALLOCATE(nm_tracer(n_tracers))
  nm_tracer(:) = 'XXXXXXXXXX'
  nm_tracer(1:jpctr) = advt(1:jpctr)

  !   This block replicates the code in UKCA_EMISSION_CTL. It is only used by
  !   the Newton-Raphson (N-R) solver, which needs to get molmass and lbc_molmass
  !   before passing them to ASAD_EMISSIONS_DIAGNOSTICS for the calculation
  !   of ASAD diagnostics in the right units.
  IF (L_ukca_chem .AND. method == int_method_nr) THEN    ! N-R solver

    !     Allocate and initialise array with molar masses of emitted species
    IF (ANY(em_chem_spec == 'NO_aircrft')) THEN
      ALLOCATE (molmass (n_use_emissions-1))
    ELSE
      ALLOCATE (molmass (n_use_emissions))
    END IF
    molmass(:) = 0.0

    !     Get the molar masses of emitted species
    DO k = 1, SIZE(molmass)
      molmass (k) = get_molmass (em_chem_spec(k))
    END DO

    ! Check if all the emitted species have a valid molecular weight
    IF (ANY(molmass(:) < small_number)) THEN
      n = 0
      DO k = 1, SIZE(molmass)
        IF (molmass(k) < small_number) THEN
          cmessage = ' Species: ' // TRIM (em_chem_spec(k)) //     &
                     ' missing from molmass list.'
          errcode  = -k
          CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
          n = n + 1
        END IF
      END DO
      IF (n > 0) THEN
        cmessage = ' Species missing from molmass list'
        errcode  = n
        CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
      END IF
    END IF  ! is any molmass < small_number ?

    IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc) THEN

      ! Set MMRs of N2O, CFCs and halons to specified global constants,
      ! which may follow a time dependence. Set values of inorganic
      ! Cl and Br compounds to 0 at lower boundary.
      ! Adjust this if new stratospheric species are introduced.
      ! The diagnostic (surface emission) contains the additional
      ! amount of tracer in mols added globally, which is negative
      ! in case of sink gases.

      ALLOCATE (lbc_molmass(n_boundary_vals))
      lbc_molmass(:) = 0.0

      !       Get the molar masses of species with lower boundary condition
      DO k = 1, n_boundary_vals
        lbc_molmass (k) = get_molmass (lbc_spec(k))
      END DO

      n = 0
      DO k = 1, n_boundary_vals
        IF (lbc_molmass(k) < small_number) THEN
          cmessage = 'Species ' // TRIM (lbc_spec(k)) //                  &
                    ' missing from lower boundary condition molmass list.'
          errcode  = -k
          CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
          n = n + 1
        END IF
      END DO
      IF (n > 0) THEN
        cmessage = 'Species missing from lower boundary condition ' //    &
                   'molmass list'
        errcode  = n
        CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
      END IF

    END IF             ! L_ukca_strat, etc
  END IF               ! N-R solver

  IF (PrintStatus > PrStatus_Normal) THEN
    DO k = 1, jpctr
      WRITE (umMessage, '(A,I3,A,A)') 'UKCA_ADD_EMISS - tracer ', k, ': ',&
                                       advt(k)
      CALL umPrint(umMessage,src='ukca_add_emiss')
    END DO
    DO imode = 1, nmodes
      IF (mode(imode)) THEN
        k = nmr_index(imode)
        WRITE (umMessage, '(A,I3,A,I2)') 'UKCA_ADD_EMISS - tracer ', &
          k, ': aerosol number for mode ', imode
        CALL umPrint(umMessage,src='ukca_add_emiss')
        DO icp = 1, ncp
          IF (component(imode,icp)) THEN
            k = mmr_index(imode, icp)
            WRITE (umMessage, '(A,I3,A,I2,A,I2)') 'UKCA_ADD_EMISS - tracer ',&
              k, ': aerosol mass for mode ', imode, ', component ', icp
            CALL umPrint(umMessage,src='ukca_add_emiss')
          END IF  ! component(icp)
        END DO  ! icp
      END IF  ! mode(imode)
    END DO  ! imode

    DO l = 1, n_use_emissions
      WRITE(umMessage,'(A,I3,A,A)')'UKCA_ADD_EMISS - non-interactive emiss ',&
                              l, ': ', em_chem_spec(l)
      CALL umPrint(umMessage,src='ukca_add_emiss')
    END DO
  END IF

  !-----------------------------------------------------------------------
  ! Identify and get all the different scaling factors used
  ! to account for hour to hour variability of emissions.
  !
  ! First find total number of different hourly profiles used.
  ! We must always count the first entry. Then we loop over all
  ! the other hourly factors in the emissions structure (starting
  ! with the second one). If none of the previous entries match
  ! this one, then we have another new type of hourly profile and
  ! so we increment num_hourly_profs.
  num_hourly_profs = 1
  DO l = 2, num_em_flds
    IF (.NOT. ANY (emissions(1:l-1)%hourly_fact ==   &
                   emissions(l)%hourly_fact) ) THEN
      num_hourly_profs = num_hourly_profs + 1
    END IF
  END DO

  WRITE (umMessage,'(A,1x,I2)')                                             &
   'Number of independent hourly emission profiles found: ', num_hourly_profs
  CALL umPrint(umMessage,src='ukca_add_emiss')

  ! Allocate arrays which will contain info on all the different
  ! hourly profiles that will be used
  ALLOCATE (hourly_prof_name      (num_hourly_profs))
  ALLOCATE (field_name            (num_hourly_profs))
  ALLOCATE (hourly_scaling_all    (0:23, num_hourly_profs))
  ALLOCATE (hourly_scaling_all_2d (row_length, rows, num_hourly_profs))

  ! The arrays with hourly scaling will be filled at each time step
  ! because they will depend on the local time. Here we only initialise
  ! and fill in two arrays:
  ! * array with names of independent hourly profiles found
  ! * array with names of the first field in the list of all emissions
  ! for which each hourly profile that was found.
  hourly_prof_name (:) = ''
  field_name       (:) = ''

  ! At least one type of hourly profile should have been found
  n = 1
  hourly_prof_name (1) = emissions(1)%hourly_fact
  field_name       (1) = emissions(1)%var_name

  ! Continue filling in if more than one type of hourly profile
  IF (num_hourly_profs > 1) THEN
    DO l = 2, num_em_flds
      IF (.NOT. ANY (emissions(1:l-1)%hourly_fact ==   &
                   emissions(l)%hourly_fact) ) THEN
        n = n + 1
        hourly_prof_name (n) = emissions(l)%hourly_fact
        field_name       (n) = emissions(l)%var_name
      END IF
    END DO
  END IF
  firstcall = .FALSE.
END IF                  ! firstcall

!-----------------------------------------------------------------------
! Continue with hourly profiling of emissions, now to fill in the
! arrays with hourly scaling which will change for each timestep.
!
! Initialise the arrays that were allocated in the first call
hourly_scaling_all    (:,:)   = 0.0
hourly_scaling_all_2d (:,:,:) = 0.0
!
! Fill in array 'hourly_scaling_all' which contains all the different
! 1-D hourly scalings (from 0 to 23 h)
DO p = 1, num_hourly_profs
  CALL hourly_emiss_factors (hourly_prof_name (p), field_name (p) , &
                             hourly_scaling (0:23))
  hourly_scaling_all (0:23, p) = hourly_scaling (0:23)
END DO
!
! Calculate local hour as a function of the true longitude and from
! it a set of grid specific 2-D hourly factors: hourly_scaling_all_2d
!
! The calculation of the local hour is an approximation. With the loops below
! we simply create up to 24 time zones of 15 degrees longitude each one
! covering the whole globe (or fewer zones in the case of a limited area
! model configuration) and calculate a local time for each of them.
DO j = 1, rows
  DO i = 1, row_length
    ! Take nearest integer so that at 0Z all cells of -7.5 to +7.5
    ! degrees use the hourly scaling indexed at 0
    t_local = NINT (true_longitude(i,j)*24.0/(2.0*pi)) + ihour
    !
    ! Ensure that t_local is bounded to be between 0 and 23
    IF (t_local >= 24) t_local = t_local - 24
    IF (t_local <   0) t_local = t_local + 24

    DO p = 1, num_hourly_profs
      ! calculate scaling for diurnal cycle as function of local time
      hourly_scaling_all_2d (i, j, p) = hourly_scaling_all (t_local, p)
    END DO
  END DO
END DO

! --------------------------------------------------------------------
! Go through all emission fields and add them to em_field (:,:,:,k)
! for the corresponding tracer k. This is done by identifying
! the emitted tracer and adding the corresponding emissions%values,
! which is already in units 'kg tracer m-2 s-1'.
!
! Here we also scale the emission field to account for temporal
! variability (diurnal and weekly) and for the vertical
! spread/distribution of the emissions, as well as for removal of 
! particulate fraction in the case of SO2.
DO l = 1, num_em_flds       ! loop over emission fields
  k = -1

  mapped_tracer = ' '
  IF (emissions(l)%l_mode) THEN
    ! Map mode emission to appropriate index in tracers array
    SELECT CASE(emissions(l)%moment)
    CASE (0)  ! number
      k = nmr_index(emissions(l)%mode)
    CASE (3)  ! mass
      k = mmr_index(emissions(l)%mode, emissions(l)%component)
    CASE DEFAULT
      cmessage = 'Moment must equal 0 or 3 (number or mass).'
      errcode = l
      CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
    END SELECT

  ELSE
    ! Map chemistry emissions to appropriate tracers
    SELECT CASE (TRIM(emissions(l)%tracer_name))
    CASE ('NO_aircrft', 'NO_lightng')
      mapped_tracer =   'NO        '

    CASE ('SO2_low', 'SO2_high', 'SO2_nat')
      mapped_tracer =   'SO2       '

    CASE ('CH4_wetlnd')
      mapped_tracer =   'CH4       '

    CASE ('NVOC')
      mapped_tracer =   'MeOH      '

    CASE DEFAULT
      mapped_tracer = emissions(l)%tracer_name
    END SELECT

    !   Loop over tracers to find the position k in em_field (:,:,:,k)
    !   to add there the new emission field
    DO m = 1, jpctr          ! loop over tracers
      IF (mapped_tracer == advt (m)) THEN
        k = m
        EXIT
      END IF
    END DO

  END IF  ! emissions(l)%l_mode


  !   If k not set then field is neither chemical or mode emission (i.e. is
  !   offline aerosol emission) and the rest of this loop should be skipped.  
  !   If it is not an offline aerosol emission, there is an error.
  IF (k < 0) THEN
    IF (ANY(aero_ems_species == emissions(l)%tracer_name)) THEN
      CYCLE
    ELSE
      IF (emissions(l)%l_mode) THEN
        WRITE(cmessage, '(a, 3i3)') 'Mode emission not mapped to tracer' // &
                   ' array: ' // TRIM(emissions(l)%var_name) //          &
                   'Mode, moment, cpt: ', emissions(l)%mode,                &
                   emissions(l)%moment, emissions(l)%component
      ELSE
        cmessage = 'Emitted species ' // TRIM (emissions(l)%tracer_name) // &
                   ' cannot be mapped to any atmospheric tracer'
      END IF
      errcode  = l
      CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
    END IF  ! ANY(aero_ems_species == mapped_tracer)
  END IF    ! k < 0

  ! For modal H2S04 emissions from SO2, set scaling to get particulate fraction
  ! and to include convertion from kg(SO2) to kg(H2SO4).
  ! For SO2 species set scaling factor to remove particulate fraction.
  ! Scale biomass burning aerosol emissions if required.
  aero_scaling = 1.0
  IF (emissions(l)%l_mode) THEN
    IF (emissions(l)%l_mode_so2) THEN  

      aero_scaling = mode_parfrac/100.0 * mm(cp_su)/mm_gas(msotwo)

    ELSE IF (l_ukca_scale_biom_aer_ems .AND. emissions(l)%l_mode_biom) THEN

      aero_scaling = biom_aer_ems_scaling

    END IF
  ELSE IF (TRIM(mapped_tracer) == 'SO2') THEN

    aero_scaling = (1.0 - mode_parfrac/100.0)

  END IF

  !   Get scaling factor to account for possible hour to hour variability
  !
  !   First look for the index containing the right hourly factor
  !   for the emission field
  pp = -1
  DO p = 1, num_hourly_profs
    IF (emissions(l)%hourly_fact == hourly_prof_name (p)) THEN
      pp = p
      EXIT
    END IF
  END DO
  !
  !   Stop with error message if hourly profile not found
  IF (pp < 0) THEN
    cmessage = 'Emitted field ' // TRIM (emissions(l)%var_name) // &
               ' cannot be mapped to any hourly profile'
    errcode  = -pp
    CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
  END IF
  !
  !   Get the 2-D hourly scaling (function of lon/lat) for the index found
  hourly_scaling_2d (:,:) = hourly_scaling_all_2d (:,:,pp)

  !   Get scaling factor to account for possible day to day variability, i.e.
  !   for a weekely cycle in the emissions. Note that this factor is set to 1
  !   in the case of a 360-day calendar because a weekly cycle does not have
  !   any meaning in that context.
  IF (lcal360) THEN
    fdaily_scaling = 1.0

    ! Do not allow attributes indicating a specific weekly cycle
    IF (emissions(l)%daily_fact  /=  'none' .AND.                         &
        emissions(l)%daily_fact  /=  '') THEN
      cmessage = 'Weekly cycle '     // TRIM (emissions(l)%daily_fact) // &
                 ' not allowed for ' // TRIM (emissions(l)%var_name)   // &
                 '  in the case of a 360-day calendar'
      errcode = 1
      CALL ereport ('UKCA_ADD_EMISS', errcode, cmessage)
    END IF
  ELSE
    ! Get factor for all days in the week
    CALL daily_emiss_factors (emissions(l)%daily_fact,  &
                              emissions(l)%var_name,    &
                              daily_scaling(1:7))
    ! Get factor for the current day of the week
    fdaily_scaling = daily_scaling (iday_week)
  END IF

  !   Note that the units conversion (applying the so called base scaling)
  !   is done in ukca_new_emiss_ctl. The vertical factors are stored
  !   in the emissions structure (from the very first time step of
  !   ukca_new_emiss_ctl) and will be applied below.

  !   When L_ukca_prescribech4 is true then make sure that the corresponding
  !   emission field is zero and prescribe the surface CH4 mass mixing ratio.
  !   For other tracers and for non-prescribed CH4 add the values from
  !   the emission structure (with the appropriate corrections) to
  !   the emission field.
  IF (mapped_tracer == 'CH4       ' .AND. L_ukca_prescribech4) THEN
    em_field (:,:,:,k)             = 0.0
    tracers  (:,:,surface_level,k) = um_ch4_for_ukca
  ELSE


    !     Set emiss diagnostics to zero to ensure that no double-counting
    !     in the case of surface emissions
    emissions(l)%diags (:,:,:) = 0.0

    DO ilev = 1, model_levels
      !       First create 2D field with values of emissions (already
      !       corrected for units) to which we apply temporal profiles
      !       as well as the vertical profile for the given model level
      !       and scalings specific to aerosols if applicable.
      !       Also fill in emission diagnostics.
      IF (emissions(l)%three_dim) THEN
        em_field_2d (:,:) = emissions(l)%values(:,:,ilev)    *      &
                            hourly_scaling_2d (:,:)          *      &
                            fdaily_scaling                   *      &
                            aero_scaling                     *      &
                            emissions(l)%vert_scaling_3d (:,:,ilev)
        emissions(l)%diags (:,:,ilev) = em_field_2d (:,:)
      ELSE
        em_field_2d (:,:) = emissions(l)%values(:,:,1)       *      &
                            hourly_scaling_2d (:,:)          *      &
                            fdaily_scaling                   *      &
                            aero_scaling                     *      &
                            emissions(l)%vert_scaling_3d (:,:,ilev)

        !         Column-integrated diagnostics for 2-D surface emiss:
        !         Add vertically- and time-scaled emissions at all levels
        !         and get a total in the column.
        emissions(l)%diags (:,:,1) = emissions(l)%diags (:,:,1) + &
                                     em_field_2d (:,:)
      END IF

      !       Add the 2D emission field (with vertical and temporal profiles
      !       already applied) to the total emission field
      em_field (:,:,ilev,k) = em_field (:,:,ilev,k) + em_field_2d (:,:)
    END DO  ! loop over model_levels
  END IF    ! CH4 & L_ukca_prescribech4 ?

END DO      ! loop over emission fields

! --------------------------------------------------------------------
! After filling in all emission fields now loop over tracers to
! first add emissions and later call boundary layer mixing
tr_flux       = 0.0
surf_dep_flux = 0.0

DO k = 1, jpctr

  !   For stratospheric chemistry schemes some species should
  !   be set to lower boundary conditions. When that is the
  !   case then set emissions to difference of tracer at surface
  !   and intended value, scaled with mass in gridbox,
  !   area and timestep, to turn it into a surface emission rate.
  !   No emissions for other model levels
  IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop) THEN
    DO l = 1, n_boundary_vals
      IF (advt(k) == lbc_spec(l)) THEN
        em_field (:,:,:,            k) = 0
        em_field (:,:,surface_level,k) =                               &
                 (lbc_mmr(l) - tracers(:,:,surface_level,k))        *  &
                 mass(:,:,surface_level) / surf_area(:,:) / timestep
      END IF
    END DO
  END IF   ! L_ukca_strat, etc

END DO  ! jpctr

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                     &
!$OMP PRIVATE(ilev, k)                                               &
!$OMP SHARED(alpha_cd, bl_levels, dtrdz_charney_grid,                &
!$OMP        em_field, exner_rho_levels, kent, kent_dsc,             &
!$OMP        ml_depth, model_levels, n_tracers, nm_tracer,           &
!$OMP        q, qcf, qcl,                                            &
!$OMP        res_factor, rho_aresist, rho_r2, rhokh_mix,             &
!$OMP        row_length, rows, surf_dep_flux,                        &
!$OMP        t_frac, t_frac_dsc, theta, timestep, tr_flux, tracers,  &
!$OMP        we_lim, we_lim_dsc, z_half, zhsc, zrzi, zrzi_dsc)
DO k = 1, n_tracers
  !   Add emissions over all model layers except at surface
  DO ilev = 2, model_levels
      CALL trsrce(                                                   &
        rows, row_length, 0, 0, 0, 0,                                &
        theta, q, qcl, qcf , exner_rho_levels, rho_r2,               &
        tracers(:,:,ilev,k), em_field(:,:,ilev,k), ilev,             &
        timestep, 1, 1, 0.0)
  END DO

  !   Call boundary layer mixing and add surface emissions.
  !   Exclude H2O tracer here if advected.
  IF (nm_tracer(k) /= 'H2O       ') THEN
    CALL tr_mix (                                               &
          bl_levels, alpha_cd,                                  &
          rhokh_mix (:,:,2:), rho_aresist, dtrdz_charney_grid,  &
          em_field (:,:,surface_level,k), res_factor (:,:),     &
          kent,     we_lim,     t_frac,     zrzi,               &
          kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,           &
          ml_depth, zhsc, z_half,                               &
          tracers       (:, :, 1:bl_levels, k),                 &
          tr_flux       (:, :, 1:bl_levels, k),                 &
          surf_dep_flux (:, :, k)                               &
          )
  END IF

END DO
!$OMP END PARALLEL DO

! --------------------------------------------------------------------
! --------------------------------------------------------------------
! Diagnose emissions; in-situ diagnostics that take samples directly
! at the point in the code where emissions are added to the tracers
! (see above). Whatever is currently stored in "em_field", the array
! that stores the sum of all emissions for each emitted species, is
! *directly* written to a corresponding diagnostic in section 50. By
! probing emissions in this way any ambiguity is removed and any
! potential double counting of emissions will be spotted.
DO k = 1, n_tracers

  IF (nm_tracer(k) == 'C5H8      ') THEN
  
! In-situ emission diagnostics for C5H8 (isoprene) --- STASHitem 50300
    section =  50
    item    = 300
    
    IF (sf(item,section)) THEN
    
! DEPENDS ON: copydiag 
      CALL copydiag (stashwork50(si(item,section,1)),           & 
        em_field (:,:,surface_level,k),                         & 
        row_length,rows,0,0,0,0, at_extremity,                  & 
        1,section,item,icode,cmessage)
    
    END IF

    IF (icode >  0) THEN
      errcode = section*1000+item
      CALL ereport ('In-situ diagnostics', errcode, cmessage)
    END IF

  ELSE IF(nm_tracer(k) == 'Monoterp  ') THEN
  
! In-situ emission diagnostics for C10H16 ((mono-)terpenes) --- STASHitem 50301
    section =  50
    item    = 301
    
    IF (sf(item,section)) THEN 
    
! DEPENDS ON: copydiag 
      CALL copydiag (stashwork50(si(item,section,1)),           & 
        em_field (:,:,surface_level,k),                         & 
        row_length,rows,0,0,0,0, at_extremity,                  & 
        1,section,item,icode,cmessage)
    
    END IF

    IF (icode >  0) THEN
      errcode = section*1000+item
      CALL ereport ('In-situ diagnostics', errcode, cmessage)
    END IF
  
  ELSE IF(nm_tracer(k) == 'MeOH      ' .OR.                     &
          nm_tracer(k) == 'CH3OH     ') THEN
  
! In-situ emission diagnostics for MeOH (methanol) --- STASHitem 50302
! Note that methanol can be referred to with two different names
! depending on the chemistry scheme.
    section =  50
    item    = 302
    
    IF (sf(item,section)) THEN 
    
! DEPENDS ON: copydiag 
      CALL copydiag (stashwork50(si(item,section,1)),           & 
        em_field (:,:,surface_level,k),                         & 
        row_length,rows,0,0,0,0, at_extremity,                  & 
        1,section,item,icode,cmessage)
    
    END IF

    IF (icode >  0) THEN
      errcode = section*1000+item
      CALL ereport ('In-situ diagnostics', errcode, cmessage)
    END IF
    
  ELSE IF(nm_tracer(k) == 'Me2CO     ') THEN
  
! In-situ emission diagnostics for Me2CO (acetone) --- STASHitem 50303
    section =  50
    item    = 303
    
    IF (sf(item,section)) THEN
    
! DEPENDS ON: copydiag 
      CALL copydiag (stashwork50(si(item,section,1)),           & 
        em_field (:,:,surface_level,k),                         & 
        row_length,rows,0,0,0,0, at_extremity,                  & 
        1,section,item,icode,cmessage)
    
    END IF

    IF (icode >  0) THEN
      errcode = section*1000+item
      CALL ereport ('In-situ diagnostics', errcode, cmessage)
    END IF
  
  END IF
  
! End of In-situ Diagnostics
! --------------------------------------------------------------------
! --------------------------------------------------------------------

END DO                       ! end of loop over tracers

! --------------------------------------------------------------------
! Calculate emissions diagnostics for ASAD.
! Note that we are only passing emiss at lowest model level.
IF (L_asad_use_chem_diags .AND. L_asad_use_surf_ems) THEN
  CALL asad_emissions_diagnostics (          &
       row_length,                           &
       rows,                                 &
       model_levels,                         &
       jpctr,                                &
       em_field (:,:,surface_level,1:jpctr), &  ! 3-D sfc field to replicate
                                                ! call in ukca_emission_ctl
       surf_area,                            &
       timestep,                             &
       n_use_emissions,                      &
       n_boundary_vals,                      &
       em_chem_spec,                         &
       lbc_spec,                             &
       molmass,                              &
       lbc_molmass,                          &
       ierr)
END IF

! --------------------------------------------------------------------
! ASAD diagnostics - aircraft NOx emissions.
IF (L_asad_use_chem_diags .AND. L_asad_use_air_ems &
    .AND. .NOT. (l_ukca_offline .OR. l_ukca_offline_be)) THEN
  DO l = 1, num_em_flds
    IF (emissions(l)%tracer_name == 'NO_aircrft') THEN
      ! Note that emissions(l)%values (:,:,:) is already expressed as
      ! kg(NO) m-2 s-1 because it was multiplied by base_scaling
      ! in UKCA_NEW_EMISS_CTL. We will therefore use that field
      ! for ASAD diagnostics without any further conversions.
      conv_aircraftems = emissions(l)%values (:,:,:)

      ! conv_aircraftems    ==> aircraft NOx emiss, in kg(NO)/m^2/s
      ! aircraft_emissions  ==> type flag used in ASAD diags routines
      CALL asad_3D_emissions_diagnostics(                             &
               row_length, rows, model_levels,                        &
               inox, conv_aircraftems,  surf_area,                    &
               totnodens, volume, mass,                               &
               m_no, timestep, aircraft_emissions,                    &
               ierr)
      EXIT
    END IF
  END DO
END IF

! --------------------------------------------------------------------

! Emission of volcanic SO2 from explosive volcanic eruptions
! into stratosphere
IF ((L_ukca_strat .OR. L_ukca_stratcfc .OR.                       &
     L_ukca_strattrop) .AND. iso2 > 0  .AND.                      &
     l_ukca_so2ems_expvolc ) THEN  

  !   Diagnostics - volcanic SO2 emissions
  IF (L_asad_use_chem_diags .AND. L_asad_use_volc_ems) THEN
    ALLOCATE (tmp3dems (row_length, rows, model_levels))
    tmp3dems (:,:,:) = tracers(1:row_length, 1:rows, 1:model_levels, iso2)
  END IF

  CALL ukca_volcanic_so2 (                                        &
           tracers (1:row_length, 1:rows, :, iso2),               &
           mass, theta_latitude, true_longitude,                  &
           delta_phi, delta_lambda,                               &
           row_length, rows, model_levels,                        &
           iyear, imonth, iday, timestep,                         &
           tropopause_height,                                     &
           r_theta_levels (1:row_length, 1:rows, 1:model_levels))

  !   Diagnostics - SO2 emissions from explosive volcanic eruptions
  IF (L_asad_use_chem_diags .AND. L_asad_use_volc_ems) THEN
    tmp3dems (:,:,:) =                                            &
     tracers (1:row_length, 1:rows, :, iso2) - tmp3dems(:,:,:)

    !     tmp3dems               ==> SO2 emission field
    !     volcanic_emissions     ==> type flag used in ASAD diagnostic routines
    CALL asad_3D_emissions_diagnostics(                           &
         row_length, rows, model_levels,                          &
         iso2, tmp3dems, surf_area,                               &
         totnodens, volume, mass,                                 &
         m_so2, timestep, volcanic_emissions,                     &
         ierr)
    DEALLOCATE (tmp3dems)
  END IF

END IF    ! L_ukca_strat, etc.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_add_emiss

END MODULE ukca_add_emiss_mod
