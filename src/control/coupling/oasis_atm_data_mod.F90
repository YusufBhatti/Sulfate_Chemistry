! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE oasis_atm_data_mod
  ! Description:
  ! Useful data and arrays for use with OASIS3-MCT coupling.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !=====================================================================

USE um_types
USE field_types, ONLY: fld_type_u, fld_type_v, fld_type_p
USE conversions_mod, ONLY: zerodegc

IMPLICIT NONE

!=====================================================
! Things needed for all models
!=====================================================

INTEGER (KIND=integer64) :: oasis_couple_ts_ao ! Coupling tstep len (atm->ocn)
INTEGER (KIND=integer64) :: oasis_couple_ts_oa ! Coupling tstep len (ocn->atm)

! Cater for future developments in atmos <-> chemistry coupling
INTEGER (KIND=integer64) :: oasis_cpl_ts_put_hyb ! Coupling tstep len
                                                 ! (hybrid put)
INTEGER (KIND=integer64) :: oasis_cpl_ts_get_hyb ! Coupling tstep len 
                                                 ! (hybrid get)
INTEGER (KIND=integer64) :: oasis_cpl_ts_put_hyb_stats ! Coupling tstep len
                                                 ! for stats (hybrid put)
INTEGER (KIND=integer64) :: oasis_cpl_ts_get_hyb_stats ! Coupling tstep len
                                                 ! for stats (hybrid get)

! Cater for future developments in atmos <-> wave coupling
INTEGER (KIND=integer64) :: oasis_couple_ts_aw ! Coupling tstep len (atm->wav)
INTEGER (KIND=integer64) :: oasis_couple_ts_wa ! Coupling tstep len (wav->atm)

! Logical switches to indicate whether a coupling exchange is needed
LOGICAL :: get_step_oa        ! True on a time step when we need to put data
                              ! to OASIS for ocean coupling, false if not.
LOGICAL :: put_step_ao        ! True on a time step when we need to get ocean
                              ! data from OASIS, false if not.
LOGICAL :: put_step_hyb       ! True on a time step when we need to put data
                              ! to OASIS for chemistry coupling between hybrid
                              ! components
LOGICAL :: get_step_hyb       ! True on a time step when we need to get data
                              ! from OASIS for chemistry coupling between 
                              ! hybrid components
LOGICAL :: put_step_hyb_stats ! True on a time step when we need to put stats 
                              ! data to OASIS for chemistry coupling between 
                              ! hybrid components
LOGICAL :: get_step_hyb_stats ! True on a time step when we need to get stats 
                              ! data from OASIS for chemistry coupling between
                              ! hybrid components
LOGICAL :: get_step_wa        ! True on a timestep when we need to put to OASIS
                              ! for wave coupling, false if not.
LOGICAL :: put_step_aw        ! True on a timestep when we need to get wave 
                              ! data from OASIS, false if not.
LOGICAL :: cpl_update_step    ! True when we need to make sure
                              ! the D1 prognostic data is updated
                              ! in preparation for dump creation or
                              ! coupling actions.

! Address pointers for atmosphere-to-ocean (or any other
! potentially coupled component) fields which may be exchanged through
! any internal or external coupling mechanism.

INTEGER:: ja_taux       ! Surface x-windstress on atmos grid
INTEGER:: ja_tauy       ! Surface y-windstress on atmos grid
INTEGER:: ja_w10        ! 10m wind on atmos grid
INTEGER:: ja_solar      ! Net downward SW at surf on atmos grid
INTEGER:: ja_blue       ! Net blueband SW at surf on atmos grid
INTEGER:: ja_evap       ! Net evaporation over sea on atmos grid
INTEGER:: ja_longwave   ! Net downward LW at surf on atmos grid
INTEGER:: ja_sensible   ! Sensible heat flux over sea on atmos grid
INTEGER:: ja_lssnow     ! Large-scale snowfall rate on atmos grid
INTEGER:: ja_cvsnow     ! Convective snowfall rate on atmos grid
INTEGER:: ja_lsrain     ! Large-scale rainfall rate on atmos grid
INTEGER:: ja_cvrain     ! Convective rainfall rate on atmos grid
INTEGER:: ja_riverout   ! Total river outflow on atmos grid
INTEGER:: ja_sublim     ! Sublimation on atmos grid
INTEGER:: ja_fcondtopn  ! Diff heat thro ice (ncat)
INTEGER:: ja_topmeltn   ! Seaice top melting flux (ncat)
INTEGER:: ja_tstar_sicen ! Sea ice surface skin temp on categories (K)
INTEGER:: ja_mslp       ! Mean Sea Level Pressure 

! STASH indecies for the MEDUSA to UKCA coupling
INTEGER:: ja_CO2        ! It is used for pCO2 coupling, but points at CO2 
INTEGER:: ja_dust_dep   ! Total dust deposition flux

! OASIS-defined or related variables
INTEGER (KIND=integer32) :: oasis_comp_id

INTEGER (KIND=integer64) :: comm_in  ! OASIS defined communicator

#if defined(MCT)
!======================================================
! Things only needed when OASIS3-MCT is active.
!======================================================

!=================================================
! General coupling related variables 
!=================================================

! Set up the names of the various model components which may be
! involved in coupling. These are used to compare against the contents of
! the transient namelist.
CHARACTER (LEN=3), PARAMETER :: atm_component="ATM"
CHARACTER (LEN=3), PARAMETER :: jnr_component="JNR"
CHARACTER (LEN=3), PARAMETER :: ocn_component="OCN"

! The number of coupling fields sent in each direction
INTEGER (KIND=integer64) :: transient_o2a_count
INTEGER (KIND=integer64) :: transient_a2o_count
INTEGER (KIND=integer64) :: transient_c2a_count
INTEGER (KIND=integer64) :: transient_a2c_count

INTEGER (KIND=integer32) :: icpl  ! Indicates to OASIS whether this PE
                                  ! is involved in coupling or not

INTEGER (KIND=integer32), PARAMETER :: icpl_atm = 1  ! Indicates an
                                  ! atmos coupling process

! Standard OASIS partition types likely to be useful
INTEGER (KIND=integer32), PARAMETER :: PartitionSerial = 0
INTEGER (KIND=integer32), PARAMETER :: PartitionBox = 2

! The current prism timestep.
REAL :: prism_timestep
INTEGER (KIND=integer64) :: prism_nsec  ! number of model seconds

! Full set of i and j indices to be used when calculating
! 2nd order conservative field gradients under O3-MCT if
! required.
! Note that we refer to things in terms of i and j rather than
! (W)est-(E)ast and (S)outh-(N)orth because this is intended
! to work on rotated grids w/o any special processing to take
! account of rotation (since we expect that remapping weights
! are generated  based on the same ROTATED grids if such grids are
! employed by the component.)
TYPE gradient_index
  INTEGER (KIND=integer32) :: im1  ! Lower i bounds
  INTEGER (KIND=integer32) :: ip1  ! Upper i bounds
  INTEGER (KIND=integer32) :: jm1  ! Lower j bounds
  INTEGER (KIND=integer32) :: jp1  ! Upper j bounds
  REAL                     :: di   ! spacing following i indexing
  REAL                     :: dj   ! spacing following j indexing
  REAL                     :: Rdi  ! Reciprocal i spacing
  REAL                     :: Rdj  ! Reciprocal j spacing
  REAL                     :: dlat_row ! True lat spacing following i (row)
  REAL                     :: dlat_col ! True lat spacing following j (column)
  REAL                     :: dlon_row ! True lon spacing following i (row)
  REAL                     :: dlon_col ! True lon spacing following j (column) 

END TYPE gradient_index

TYPE (gradient_index), ALLOCATABLE :: GradIndex(:,:)

INTEGER (KIND=integer64) :: max_transients

! Transient field description components
TYPE transient
  INTEGER (KIND=integer32) :: indx        ! Unique identifier
  INTEGER (KIND=integer32) :: field_id    ! Unique field identifier
  INTEGER (KIND=integer32) :: order       ! 1st or 2nd order mapping flag
  INTEGER (KIND=integer32), ALLOCATABLE :: oasis_id(:)  ! OASIS defined ID 
  INTEGER (KIND=integer32) :: oasis_info  ! Return code from OASIS ops
  CHARACTER (LEN=3) :: c_from             ! Source component of field
  CHARACTER (LEN=3) :: c_to               ! Target component of field
  CHARACTER (LEN=1) :: grid               ! Grid cell type
  CHARACTER (LEN=8), ALLOCATABLE :: name(:)  ! Transient name
  LOGICAL (KIND=logical64):: ice_min      ! Ice minima processing requd.
  REAL (KIND=real64), POINTER :: field(:,:)      ! Pointer to 2D field
  REAL (KIND=real64), POINTER :: field3d(:,:,:)  ! Pointer to 3D field.
  INTEGER (KIND=integer32) :: n_slices    ! Number of horizontal levels
END TYPE transient

CHARACTER (LEN=3), PARAMETER :: my_component="ATM"

! Distinguish between different source/target combinations
! for transients.
TYPE (transient), ALLOCATABLE :: transient_a2o(:)
TYPE (transient), ALLOCATABLE :: transient_o2a(:)

TYPE (transient), ALLOCATABLE :: transient_a2c(:)
TYPE (transient), ALLOCATABLE :: transient_c2a(:)

! The index to deal specifically with orography.
INTEGER :: trans_indx_orog

CHARACTER (LEN=4) :: c_out
CHARACTER (LEN=4) :: c_in
CHARACTER (LEN=8) :: n_out
CHARACTER (LEN=8) :: n_in
INTEGER  ::          i_number
INTEGER  ::          f_id
INTEGER  ::          map_type
INTEGER  ::          max_lev

! Limits applied to top layer ice temperature field
REAL, PARAMETER :: TI_max = zerodegc
REAL, PARAMETER :: TI_min = 200.0

! Limit applied to ice thickness - ice thinner than this value will be
! - if multilayers set, vanished by setting ice area = 0
! - if zerolayers set, squished so that it has thickness hi_min
REAL            :: hi_min 

! Limit applied to ice concentration
REAL            :: aicenmin

NAMELIST/transcount/ max_transients, &
         transient_a2c_count,        &     
         transient_c2a_count,        &     
         transient_a2o_count,        &     
         transient_o2a_count 

NAMELIST/transfld/          &
         c_out,             &
         c_in,              &
         n_out,             &
         n_in,              &
         i_number,          &
         f_id,              &
         map_type,          &
         max_lev

! Logical switches indicate whether various coupling
! fields are active or not.

! Define flags for outgoing coupling fields
LOGICAL :: l_heatflux
LOGICAL :: l_bluerad
LOGICAL :: l_icesheet
LOGICAL :: l_solar
LOGICAL :: l_runoff
LOGICAL :: l_w10
LOGICAL :: l_train
LOGICAL :: l_tsnow
LOGICAL :: l_evap2d
LOGICAL :: l_taux
LOGICAL :: l_tauy
LOGICAL :: l_sublim(5)
LOGICAL :: l_topmeltn(5)
LOGICAL :: l_fcondtopn(5)
LOGICAL :: l_mslp 
LOGICAL :: l_tstar_sicen(5)
LOGICAL :: l_pCO2
LOGICAL :: l_dust_dep

! Define flags for incoming coupling fields
LOGICAL :: l_ocn_sst
LOGICAL :: l_ocn_u
LOGICAL :: l_ocn_v
LOGICAL :: l_ocn_freeze
LOGICAL :: l_ocn_sstfrz
LOGICAL :: l_ocn_icetn(5)
LOGICAL :: l_ocn_icekn(5)
LOGICAL :: l_ocn_freezen(5)
LOGICAL :: l_ocn_freezen_first(5)
LOGICAL :: l_ocn_snowthickn(5)
LOGICAL :: l_ocn_hicen(5)
LOGICAL :: l_pond_frac_n(5)
LOGICAL :: l_pond_depth_n(5)
LOGICAL :: l_dms_conc
LOGICAL :: l_CO2flux
LOGICAL :: l_chloro_conc   ! Surface chlorophyll conc.

! Arrays containing the "vind" numbers of the
! transients we want to use. These are used to
! define our transients using a simple loop.

! Indices used to point to various coupling fields.
INTEGER :: vind_sublim(5)
INTEGER :: vind_ocn_freezen(5)
INTEGER :: vind_ocn_freezen_first(5)
INTEGER :: vind_ocn_snowthickn(5)
INTEGER :: vind_ocn_hicen(5)
INTEGER :: vind_topmeltn(5)
INTEGER :: vind_fcondtopn(5)
INTEGER :: vind_ocn_icetn(5)
INTEGER :: vind_ocn_icekn(5)
INTEGER :: vind_pond_frac_n(5)
INTEGER :: vind_pond_depth_n(5)
INTEGER :: vind_tstar_sicen(5)

! Define indexes for outgoing coupling fields

INTEGER , PARAMETER :: vind_heatflux = 1
INTEGER , PARAMETER :: vind_bluerad = 2
INTEGER , PARAMETER :: vind_solar = 54
INTEGER , PARAMETER :: vind_runoff = 3
INTEGER , PARAMETER :: vind_w10 = 4
INTEGER , PARAMETER :: vind_train = 5
INTEGER , PARAMETER :: vind_tsnow = 6
INTEGER , PARAMETER :: vind_evap2d = 7
INTEGER , PARAMETER :: vind_icesheet_n = 72
INTEGER , PARAMETER :: vind_icesheet_s = 73


DATA vind_topmeltn       / 8, 9,10,11,12/
DATA vind_fcondtopn      /13,14,15,16,17/
DATA vind_sublim         /18,19,20,21,22/
DATA vind_tstar_sicen    /66,67,68,69,70/

INTEGER , PARAMETER :: vind_taux = 23
INTEGER , PARAMETER :: vind_tauy = 24

INTEGER , PARAMETER :: vind_pCO2 = 92
INTEGER , PARAMETER :: vind_dust_dep = 93

! Define indexes for incoming coupling fields

INTEGER , PARAMETER :: vind_ocn_sst = 25

! DATA statements are used to initialise the following because
! potentially we might have gaps under which circumstances the
! use of PARAMETER is no good.
DATA vind_ocn_freezen    /26,27,28,29,30/
DATA vind_ocn_snowthickn /31,32,33,34,35/
DATA vind_ocn_hicen      /36,37,38,39,40/
DATA vind_ocn_icetn      /41,42,43,44,45/
DATA vind_ocn_icekn      /46,47,48,49,50/

INTEGER , PARAMETER :: vind_ocn_u = 51
INTEGER , PARAMETER :: vind_ocn_v = 52
INTEGER , PARAMETER :: vind_ocn_freeze = 53
INTEGER , PARAMETER :: vind_mslp = 55 
DATA vind_pond_frac_n    /56,57,58,59,60/
DATA vind_pond_depth_n    /61,62,63,64,65/
INTEGER , PARAMETER :: vind_ocn_sstfrz = 71
DATA vind_ocn_freezen_first /81,82,83,84,85/

INTEGER , PARAMETER :: vind_dms_conc = 90
INTEGER , PARAMETER :: vind_CO2flux  = 91
INTEGER , PARAMETER :: vind_chloro_conc = 94

! Dimensions for coupling transient arrays.
! These will simply hold copies of lasize, but make
! code less long winded.
INTEGER :: oasis_jmt
INTEGER :: oasis_jmt_u
INTEGER :: oasis_jmt_v
INTEGER :: oasis_imt

! Outgoing coupling arrays or components thereof
! Anything defined without a "target" is simply
! a potential constituent value, not a direct coupling field.

REAL, ALLOCATABLE, TARGET :: taux(:,:)
REAL, ALLOCATABLE, TARGET :: tauy(:,:)

REAL, ALLOCATABLE, TARGET :: solar2d(:,:)
REAL, ALLOCATABLE, TARGET :: blue2d(:,:)
REAL, ALLOCATABLE, TARGET :: evap2d(:,:)
REAL, ALLOCATABLE :: longwave2d(:,:)
REAL, ALLOCATABLE :: sensible2d(:,:)
REAL, ALLOCATABLE, TARGET :: sublim(:,:,:)
REAL, ALLOCATABLE, TARGET :: heatflux(:,:)

REAL, ALLOCATABLE :: rainls(:,:)
REAL, ALLOCATABLE :: snowls(:,:)
REAL, ALLOCATABLE :: rainconv(:,:)
REAL, ALLOCATABLE :: snowconv(:,:)
REAL, ALLOCATABLE, TARGET :: totalrain(:,:)
REAL, ALLOCATABLE, TARGET :: totalsnow(:,:)
REAL, ALLOCATABLE, TARGET :: riverout(:,:)
REAL, ALLOCATABLE, TARGET :: calving(:,:)
REAL, ALLOCATABLE, TARGET :: w10(:,:)

REAL, ALLOCATABLE, TARGET :: topmeltn(:,:,:)
REAL, ALLOCATABLE, TARGET :: fcondtopn(:,:,:)
REAL, ALLOCATABLE, TARGET :: mslp(:,:)
REAL, ALLOCATABLE, TARGET :: tstar_sicen(:,:,:)

REAL, ALLOCATABLE, TARGET :: icenorth(:,:)
REAL, ALLOCATABLE, TARGET :: icesouth(:,:)

! Outgoing fields to MEDUSA
! Note that pCO2 is derived from CO2
REAL, ALLOCATABLE, TARGET :: pCO2_out(:,:)
REAL, ALLOCATABLE, TARGET :: dust_dep_out(:,:)

! Incoming coupling transients or components thereof
! Anything defined without a "target" is simply
! a potential constituent value, not a direct coupling field.
REAL, ALLOCATABLE, TARGET :: ocn_sst(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_sstfrz(:,:) 
REAL, ALLOCATABLE :: ocn_sst_orig(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_freeze(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_freezen(:,:,:)
REAL, ALLOCATABLE, TARGET :: ocn_freezen_first(:,:,:)
REAL, ALLOCATABLE :: ocn_hice(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_hicen(:,:,:)
REAL, ALLOCATABLE :: ocn_snowthick(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_snowthickn(:,:,:)
REAL, ALLOCATABLE, TARGET :: ocn_u(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_v(:,:)
REAL, ALLOCATABLE, TARGET :: ocn_icetn(:,:,:)
REAL, ALLOCATABLE :: ocn_icet_gb (:,:)
REAL, ALLOCATABLE, TARGET :: ocn_icekn(:,:,:)
REAL, ALLOCATABLE :: tstar_local(:,:)
REAL, ALLOCATABLE :: tstar_ssi(:,:)
REAL, ALLOCATABLE :: tstar_sice_local(:,:)
REAL, ALLOCATABLE :: tstar_sice_cat_local(:,:,:)
REAL, ALLOCATABLE :: tstar_land_local(:,:)
REAL, ALLOCATABLE, TARGET :: pond_frac_n(:,:,:)
REAL, ALLOCATABLE, TARGET :: pond_depth_n(:,:,:)

! Local land-fraction
REAL, ALLOCATABLE :: fland_loc(:,:)

! This variable is used to read in 'future' ice concentration from the ocean,
! for dividing through atmosphere-to-ice fluxes in oasis3_puta2o
REAL, ALLOCATABLE :: ice_fract_cat_future(:,:,:)

! Incoming fields from MEDUSA
REAL, ALLOCATABLE, TARGET :: dms_conc_in(:,:)
REAL, ALLOCATABLE, TARGET :: CO2flux_in(:,:)
REAL, ALLOCATABLE, TARGET :: chloro_conc_in(:,:) ! Surf. chlorophyll conc. Kg/M3

#endif

END MODULE oasis_atm_data_mod

