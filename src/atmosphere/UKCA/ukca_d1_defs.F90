! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to define fields from D1 needed for UKCA
!  these are set in UKCA_SETD1DEFS
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!  Language:  Fortran 2003
!  This code is written to UMDP3 programming standards.
!
! ######################################################################
!
MODULE ukca_d1_defs
IMPLICIT NONE

!     No of Prognostics and Diagnostics required
INTEGER, SAVE :: Nukca_D1items      ! Size of UkcaD1codes array
INTEGER, SAVE :: n_use_tracers      ! no of tracers used
INTEGER, SAVE :: n_use_emissions    ! no of emissions used
INTEGER, SAVE :: n_ntp              ! No of non-transported prognostics
INTEGER, SAVE :: n_3d_emissions     ! no of 3-D emissions used
INTEGER, SAVE :: n_in_progs         ! No of prognostics required
INTEGER, SAVE :: n_in_diags0        ! No of diagnostics (sect 0)
INTEGER, SAVE :: n_in_diags1        ! No of diagnostics (sect 1)
INTEGER, SAVE :: n_in_diags2        ! No of diagnostics (sect 2)
INTEGER, SAVE :: n_in_diags3        ! No of diagnostics (sect 3)
INTEGER, SAVE :: n_in_diags4        ! No of diagnostics (sect 4)
INTEGER, SAVE :: n_in_diags5        ! No of diagnostics (sect 5)
INTEGER, SAVE :: n_in_diags8        ! No of diagnostics (sect 8)
INTEGER, SAVE :: n_in_diags15       ! No of diagnostics (sect 15)
INTEGER, SAVE :: n_in_diags30       ! No of diagnostics (sect 30)
INTEGER, SAVE :: n_in_diags33       ! No of diagnostics (sect 33)
INTEGER, SAVE :: n_in_diags38       ! No of diagnostics (sect 38)
INTEGER, SAVE :: n_in_diags50       ! No of diagnostics (sect 50)
INTEGER, SAVE :: n_in_diags51       ! No of diagnostics (sect 51)
INTEGER, SAVE :: n_in_diags52       ! No of diagnostics (sect 52)
INTEGER, SAVE :: iemiss_first       ! item for 1st emission diag
INTEGER, SAVE :: iemiss_last        ! item for last emission diag
INTEGER, SAVE :: istrat_first       ! item for 1st strat flux diag
INTEGER, SAVE :: istrat_last        ! item for last strat flux diag
INTEGER, SAVE :: imode_first        ! item for 1st MODE diag
INTEGER, SAVE :: imode_last         ! item for last MODE diag
INTEGER, PARAMETER :: ukca_sect=34       ! stash section for UKCA
INTEGER, PARAMETER :: ukca_glomap_diag_sect=38
                                         ! stash section for UKCA GLOMAP diags
INTEGER, PARAMETER :: ukca_diag_sect=50  ! stash section for UKCA diags
INTEGER, PARAMETER :: ukca_s34_plev=51   ! stash section 34 on pressure levels
INTEGER, PARAMETER :: ukca_s50_plev=52   ! stash section 50 on pressure levels
INTEGER, PARAMETER :: item1_ukca_diags = 1    ! 1st item UKCA diags (s50)
INTEGER, PARAMETER :: item1_stratflux  = 1    ! 1st item No. Strat fluxes (s38)
INTEGER, PARAMETER :: item1_mode_diags = 201  ! 1st item No. for MODE diags (s38)
INTEGER, PARAMETER :: item1_plume_diags = 900 ! 1st item No. for plume
                                              ! scavenging diags (s38)
INTEGER, PARAMETER :: n_boundary_vals  = 7    ! No. lower boundary vals for
                                              ! stratospheric chemistry
INTEGER, PARAMETER :: nfirst_mode_tracer = 101 ! first mode tracer number
INTEGER, PARAMETER :: nmax_mode_tracers  = 39  ! max no. of mode tracers

TYPE code
  INTEGER :: section        ! section code
  INTEGER :: item           ! item code
  INTEGER :: n_levels       ! number of levels
  INTEGER :: address        ! address in D1
  INTEGER :: length         ! length of field
  INTEGER :: halo_type      ! field halo type
  INTEGER :: grid_type      ! grid type
  INTEGER :: field_type     ! field grid type
  INTEGER :: len_dim1       ! length of array dim1
  INTEGER :: len_dim2       ! length of array dim2
  INTEGER :: len_dim3       ! length of array dim3
  LOGICAL :: prognostic     ! prognostic t/f
  LOGICAL :: required       ! t/f
  CHARACTER(LEN=10) :: NAME ! species name
END TYPE code

!     Number of tracers, emissions and diagnostics, set according
!     to GUI choices in UKCA_SETD1DEFS

INTEGER, SAVE   :: n_chem_tracers     ! No. tracers for chemistry
INTEGER, SAVE   :: n_aero_tracers     ! No. tracers for chemistry
INTEGER, SAVE   :: n_nonchem_tracers  ! No. tracers non-chemistry
INTEGER, SAVE   :: n_chem_emissions   ! "   emissions        "
INTEGER, SAVE   :: nmax_strat_fluxdiags  ! Max no strat flux diags
INTEGER, SAVE   :: nmax_mode_diags    ! Max no MODE diags
INTEGER, SAVE   :: nmax_plume_diags = 33 ! Max no plume scavenging diags
INTEGER, SAVE   :: n_strat_fluxdiags  ! No. strat flux diags
INTEGER, SAVE   :: n_MODE_tracers     ! No. tracers for MODE
INTEGER, SAVE   :: n_MODE_emissions   ! "   emissions   "
INTEGER, SAVE   :: n_MODE_diags       ! "   diagnostics "
INTEGER, SAVE   :: n_dust_tracers     ! No. tracers for dust
INTEGER, SAVE   :: n_dust_emissions   ! "   emissions   "
INTEGER, SAVE   :: n_dust_diags       ! "   diagnostics "
INTEGER, SAVE   :: nr_therm           ! "   thermal reactions
INTEGER, SAVE   :: nr_phot            ! "   photolytic reactions

!     list of prognostic/diagnostics from/to D1

TYPE(code), ALLOCATABLE, SAVE :: UkcaD1Codes(:)

!     Names for tracers which have surface emissions

CHARACTER(LEN=10), ALLOCATABLE, SAVE :: em_chem_spec(:)
CHARACTER(LEN=10), ALLOCATABLE, SAVE :: em_MODE_spec(:)
CHARACTER(LEN=10), ALLOCATABLE, SAVE :: em_dust_spec(:)
! Lower BCs for stratospheric species
CHARACTER(LEN=10), ALLOCATABLE, SAVE :: lbc_spec(:)   ! species names
REAL, ALLOCATABLE, SAVE :: lbc_mmr(:)    ! mixing ratios

!     Index arrays for emissions
INTEGER, ALLOCATABLE, SAVE           :: em_index(:)

!      The flux logicals are now turned on via STASH panel requests.
LOGICAL, SAVE :: L_ukca_stratflux = .FALSE.  ! T for stratospheric fluxes
LOGICAL, SAVE :: L_ukca_mode_diags= .FALSE.  ! T for MODE diags.
LOGICAL, SAVE :: L_ukca_emiss_diags= .FALSE. ! T for emission diags.
LOGICAL, SAVE :: l_ukca_plume_diags= .FALSE. ! T for plume scavenging diags.

! STASH Codes of species used for  coupling of UKCA with radiation scheme
! and sulphur cycle of UM

!! Note: TO BE UPDATED FOR ANY SEC 34 ITEM CHANGES     !!

! species used for feedback to radiation scheme
INTEGER, PARAMETER :: i_ukca_grg_o3  =  1, i_ukca_grg_ch4 = 9 ,           &
                      i_ukca_grg_n2o = 49, i_ukca_grg_f11 = 55,           &
                      i_ukca_grg_f12 = 56, i_ukca_grg_f113= 64,           &
                      i_ukca_grg_h22 = 65, i_ukca_grg_cf3chf2 = -1,       &
                      i_ukca_grg_chf2chf2 = -1

! Item numbers for UKCA oxidants used in the CLASSIC sulphur cycle.
! Set in ukca_setup_chem - address of OH and HO2 varies with solver

!       001 : O3 MASS MIXING RATIO AFTER TSTEP
!       007 : HONO2 MASS MIXING RATIO AFTER TSTEP
!       008 : H2O2 MASS MIXING RATIO AFTER TSTEP
!           : OH MASS MIXING RATIO  AFTER TSTEP
!           : HO2 MASS MIXING RATIO AFTER TSTEP
INTEGER, SAVE :: ukca_item_sulpc(5)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_D1_DEFS'

CONTAINS

!------------------------------------------------------------------
!  Description:
!   Given a section and item, find the index in UkcaD1Codes
!
!  Method:
!   Go through all entries in ukca_d1_defs checking section and 
!   item. Set value to return when these match the requested value.
!   If value is not found return -99
!------------------------------------------------------------------

INTEGER FUNCTION ukca_d1_myindex(section, item)

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

! input arguments
INTEGER, INTENT(IN) :: section, item

! local variables
INTEGER :: i, d1index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_D1_MYINDEX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ukca_d1_myindex = -99

DO i=LBOUND(UkcaD1Codes,1), UBOUND(UkcaD1Codes,1)
  IF (UkcaD1Codes(i)%section == section .AND. UkcaD1Codes(i)%item == item) THEN
     ukca_d1_myindex = i
     EXIT
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN 

END FUNCTION ukca_d1_myindex


END MODULE ukca_d1_defs
