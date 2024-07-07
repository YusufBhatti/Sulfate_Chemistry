! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module defining ASAD arrays, variables, and parameters
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

MODULE asad_mod

USE ukca_option_mod, ONLY: jpctr, jpspec, jpdw, jpdd, jpnr, jpbk, jptk,  &
                           jphk, jppj, L_ukca_strattrop, L_ukca_strat,   &
                           L_ukca_stratcfc
USE ukca_chem_schemes_mod, ONLY: int_method_nr
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE
PUBLIC

REAL, ALLOCATABLE :: wp(:)           ! water vapour field (vmr)
REAL, ALLOCATABLE :: co2(:)          ! CO2 field (vmr)
REAL, ALLOCATABLE :: dpd(:,:)
REAL, ALLOCATABLE :: dpw(:,:)
REAL, ALLOCATABLE :: emr(:,:)
REAL, ALLOCATABLE :: fj(:,:,:)       ! Full jacobian
REAL, ALLOCATABLE :: qa(:,:)
REAL, ALLOCATABLE :: ratio(:,:)
REAL, ALLOCATABLE :: p(:)
REAL, ALLOCATABLE :: t(:)
REAL, ALLOCATABLE :: t300(:)
REAL, ALLOCATABLE :: tnd(:)          ! total number density
REAL, ALLOCATABLE :: pmintnd(:)
REAL, ALLOCATABLE :: f(:,:)          ! tracer concentrations
REAL, ALLOCATABLE :: fdot(:,:)
REAL, TARGET, ALLOCATABLE :: pd(:,:)
REAL, POINTER     :: prod(:,:)
REAL, POINTER     :: slos(:,:)
REAL, ALLOCATABLE :: y(:,:)
REAL, ALLOCATABLE :: ydot(:,:)
REAL, ALLOCATABLE :: ftilde(:,:)     ! lower order solution
REAL, ALLOCATABLE :: ej(:,:)
REAL, ALLOCATABLE :: rk(:,:)
REAL, ALLOCATABLE :: prk(:,:)
REAL, ALLOCATABLE :: deriv(:,:,:)
REAL, ALLOCATABLE :: za(:)           ! Aerosol surface area
REAL, ALLOCATABLE :: co3(:)          ! Column ozone
REAL, ALLOCATABLE :: lati(:)         ! Latitude
REAL, ALLOCATABLE :: sphno3(:)       ! Amount of HNO3 in solid phase
REAL, ALLOCATABLE :: sph2o(:)        ! Amount of H2O in solid phase
REAL, ALLOCATABLE :: depvel(:,:,:)
REAL, ALLOCATABLE :: k298(:)         ! K(298) for Henry law (M/atm)
REAL, ALLOCATABLE :: dhr(:)          ! deltaH/R (K^-1)
REAL, ALLOCATABLE :: kd298(:,:)      ! dissociation constant (M)
REAL, ALLOCATABLE :: ddhr(:,:)       ! deltaH/R
REAL, ALLOCATABLE :: ct_k298(:)      ! As above, but for constant species
REAL, ALLOCATABLE :: ct_dhr(:)       ! Allocated in ukca_chem_offline
REAL, ALLOCATABLE :: ct_kd298(:,:)   !
REAL, ALLOCATABLE :: ct_ddhr(:,:)    !
REAL, ALLOCATABLE :: ab(:,:)
REAL, ALLOCATABLE :: at(:,:)
REAL, ALLOCATABLE :: aj(:,:)
REAL, ALLOCATABLE :: ah(:,:)
REAL, ALLOCATABLE :: ztabpd(:,:)
REAL, ALLOCATABLE :: shno3(:)        ! No. density type 1 psc solid phase hno3
REAL, ALLOCATABLE :: sh2o(:)         ! No. density type 2 psc solid phase h2o
REAL, ALLOCATABLE :: fpsc1(:)        ! 1.0 if type 1 psc's are present, else 0
REAL, ALLOCATABLE :: fpsc2(:)        ! 1.0 if type 2 psc's are present, else 0
REAL, ALLOCATABLE :: spfj(:,:)       ! Sparse full Jacobian

INTEGER, ALLOCATABLE :: madvtr(:)
INTEGER, ALLOCATABLE :: majors(:)
INTEGER, ALLOCATABLE :: moffam(:)
INTEGER, ALLOCATABLE :: nodd(:)
INTEGER, ALLOCATABLE :: nltrf(:)
INTEGER, ALLOCATABLE :: nltr3(:)
INTEGER, ALLOCATABLE :: ipa(:,:)     ! Pivot information for solving jacobian
INTEGER, ALLOCATABLE :: ipa2(:)
INTEGER, ALLOCATABLE :: nltrim(:,:)
INTEGER, ALLOCATABLE :: nlpdv(:,:)
INTEGER, ALLOCATABLE :: nfrpx(:)     ! index to fractional product array nfrpx:
! one entry for each reaction. If zero, there are no fractional products. If nonzero,
! contains the array element in frpx for the first coefficient for that reaction.
INTEGER, ALLOCATABLE :: ntabfp(:,:)  ! Table used for indexing the fractional
!    products:  ntabfp(i,1) contains the species no.,
!               ntabfp(i,2) contains the reaction no., and
!               ntabfp(i,3) contains the array location in the frpx array.
INTEGER, ALLOCATABLE :: ntabpd(:,:)
INTEGER, ALLOCATABLE :: npdfr(:,:)
INTEGER, ALLOCATABLE :: ngrp(:,:)
INTEGER, ALLOCATABLE :: njcgrp(:,:)
INTEGER, ALLOCATABLE :: nprdx3(:,:,:)
INTEGER, ALLOCATABLE :: nprdx2(:,:)
INTEGER, ALLOCATABLE :: nprdx1(:)
INTEGER, ALLOCATABLE :: njacx3(:,:,:)
INTEGER, ALLOCATABLE :: njacx2(:,:)
INTEGER, ALLOCATABLE :: njacx1(:)
INTEGER, ALLOCATABLE :: nmpjac(:)
INTEGER, ALLOCATABLE :: npjac1(:,:)
INTEGER, ALLOCATABLE :: nbrkx(:)
INTEGER, ALLOCATABLE :: ntrkx(:)
INTEGER, ALLOCATABLE :: nprkx(:)
INTEGER, ALLOCATABLE :: nhrkx(:)
INTEGER, ALLOCATABLE :: nlall(:)
INTEGER, ALLOCATABLE :: nlstst(:)
INTEGER, ALLOCATABLE :: nlf(:)
INTEGER, ALLOCATABLE :: nlmajmin(:)
INTEGER, ALLOCATABLE :: nldepd(:)
INTEGER, ALLOCATABLE :: nldepw(:)
INTEGER, ALLOCATABLE :: nlemit(:)
INTEGER, ALLOCATABLE :: nldepx(:)
INTEGER, ALLOCATABLE :: njcoth(:,:)
INTEGER, ALLOCATABLE :: nmzjac(:)
INTEGER, ALLOCATABLE :: nzjac1(:,:)
INTEGER, ALLOCATABLE :: njcoss(:,:)
INTEGER, ALLOCATABLE :: nmsjac(:)
INTEGER, ALLOCATABLE :: nsjac1(:,:)
INTEGER, ALLOCATABLE :: nsspt(:)
INTEGER, ALLOCATABLE :: nspi(:,:)
INTEGER, ALLOCATABLE :: nsspi(:,:)
INTEGER, ALLOCATABLE :: nssi(:)
INTEGER, ALLOCATABLE :: nssrt(:)
INTEGER, ALLOCATABLE :: nssri(:,:)
INTEGER, ALLOCATABLE :: nssrx(:,:)
REAL, ALLOCATABLE :: frpb(:)         ! fractional product array (bimol)
REAL, ALLOCATABLE :: frpt(:)         ! fractional product array (trimol)
REAL, ALLOCATABLE :: frpj(:)         ! fractional product array (phot)
REAL, ALLOCATABLE :: frph(:)         ! fractional product array (het)
REAL, ALLOCATABLE :: frpx(:)         ! fractional product array (total)
! sparse algebra
INTEGER, ALLOCATABLE :: nonzero_map_unordered(:,:)  ! Map of nonzero entries
INTEGER, ALLOCATABLE :: modified_map(:,:) ! modified map (after decomposition)
INTEGER, ALLOCATABLE :: nonzero_map(:,:) 
                                    ! Map of nonzero entries, before reordering

INTEGER, ALLOCATABLE :: ro(:)       ! reordering of tracers to minimize fill-in

INTEGER, ALLOCATABLE :: ilcf(:)
INTEGER, ALLOCATABLE :: ilss(:)
INTEGER, ALLOCATABLE :: ilct(:)
INTEGER, ALLOCATABLE :: ilftr(:)
INTEGER, ALLOCATABLE :: ilft(:)
INTEGER, ALLOCATABLE :: ilstmin(:)

LOGICAL, ALLOCATABLE :: linfam(:,:)
LOGICAL, ALLOCATABLE :: ldepd(:)     ! T for dry deposition
LOGICAL, ALLOCATABLE :: ldepw(:)     ! T for wet deposition
LOGICAL, ALLOCATABLE :: lemit(:)     ! T for emission

CHARACTER(LEN=10), ALLOCATABLE :: advt(:)       ! advected tracers
CHARACTER(LEN=10), ALLOCATABLE :: nadvt(:)      ! non-advected species
CHARACTER(LEN=10), ALLOCATABLE :: family(:)     ! family
CHARACTER(LEN=10), ALLOCATABLE :: speci(:)      ! species names
CHARACTER(LEN=2),  ALLOCATABLE :: ctype(:)      ! species type
CHARACTER(LEN=10), ALLOCATABLE :: spb(:,:)      ! species from bimolecular rates
CHARACTER(LEN=10), ALLOCATABLE :: spt(:,:)      ! species from termolecular rates
CHARACTER(LEN=10), ALLOCATABLE :: spj(:,:)      ! species from photolysis rates
CHARACTER(LEN=10), ALLOCATABLE :: sph(:,:)      ! species from heterogenous rates

REAL, PARAMETER    :: pmin = 1.0e-20
REAL, PARAMETER    :: ptol = 1.0e-5            ! tolerance for time integration
REAL, PARAMETER    :: ftol = 1.0e-3            ! tolerance in family member iteration

INTEGER, PARAMETER :: kfphot=0
INTEGER, PARAMETER :: jpss = 16
INTEGER, PARAMETER :: jpssr = 51
INTEGER, PARAMETER :: nvar = 17
INTEGER, PARAMETER :: nllv = 17
INTEGER, PARAMETER :: ninv = 23
INTEGER, PARAMETER :: nout=71
!     nout:        Fortran channel for output in subroutine OUTVMR
INTEGER, PARAMETER :: jpem = 9      ! IS THIS A CONSTANT/USED?
INTEGER, PARAMETER :: jpeq = 2      ! dimension for dissociation arrays
INTEGER, PARAMETER :: jddept = 6
!     jddept:      Number of time periods used in dry deposition i.e.
!                  summer(day,night,24h ave), winter(day,night,24h ave)
INTEGER, PARAMETER :: jddepc = 5
!     jddepc:      Number of land use categories used in dry dep.
INTEGER, PARAMETER :: jpdwio = 56
!     jpdwio       Fortran i/o unit to read/write anything to do with
!                  wet/dry deposition
INTEGER, PARAMETER :: jpemio = 57
!     jpemio       Fortran i/o unit to read in emissions
INTEGER, PARAMETER :: jpfrpd=100
INTEGER, PARAMETER :: jpab    = 3
INTEGER, PARAMETER :: jpat    = 7
INTEGER, PARAMETER :: jpaj    = 3
INTEGER, PARAMETER :: jpah    = 3
INTEGER, PARAMETER :: jpspb   = 6
INTEGER, PARAMETER :: jpspt   = 4
INTEGER, PARAMETER :: jpspj   = 6
INTEGER, PARAMETER :: jpsph   = 6
INTEGER, PARAMETER :: jpmsp   = jpspb
INTEGER, PARAMETER :: jppjac  = 10
INTEGER, PARAMETER :: jpkargs = 10
INTEGER, PARAMETER :: jprargs = 10
INTEGER, PARAMETER :: jpcargs = 1


! Fractional product parameters - these are initialsed to jp values in asad_mod_init
INTEGER :: jpfrpb
INTEGER :: jpfrpt
INTEGER :: jpfrpj
INTEGER :: jpfrph
INTEGER :: jpfrpx

INTEGER, PARAMETER :: jpcio   = 55
INTEGER, PARAMETER :: spfjsize_max =1000 ! maximum number of
                                         ! nonzero matrix elements


CHARACTER(LEN=2),  PARAMETER :: jpfm = 'FM'     ! Family member
CHARACTER(LEN=2),  PARAMETER :: jpif = 'FT'     ! Family member depending in timestep
CHARACTER(LEN=2),  PARAMETER :: jpsp = 'TR'     ! Independent tracer
CHARACTER(LEN=2),  PARAMETER :: jpna = 'SS'     ! Steady-state species
CHARACTER(LEN=2),  PARAMETER :: jpco = 'CT'     ! Constant
CHARACTER(LEN=2),  PARAMETER :: jpcf = 'CF'     ! Constant with spatial field

LOGICAL, PARAMETER :: lvmr=.TRUE.    ! T for volume mixing ratio
LOGICAL      :: o1d_in_ss      ! T for steady state,
LOGICAL      :: o3p_in_ss      ! these are set in routine:
LOGICAL      :: n_in_ss        ! asad_mod_init
LOGICAL      :: h_in_ss        !
INTEGER, PARAMETER :: nss_o1d=1      ! indicies of deriv array
INTEGER, PARAMETER :: nss_o3p=2      !    "         "
INTEGER, PARAMETER :: nss_n=3        !    "         "
INTEGER, PARAMETER :: nss_h=4        !    "         "

REAL :: fch4,fco2,fh2,fn2,fo2
REAL    :: cdt                       ! chemistry timestep
REAL    :: peps                      !
REAL, PARAMETER :: tslimit = 1200.0  ! timestep limit for some solvers

INTEGER :: nrsteps
INTEGER :: nitnr          ! Iterations in ftoy for IMPACT solver
INTEGER :: nitfg          ! Max no of iterations in ftoy
INTEGER :: ntrf           ! Counter for tracers
INTEGER :: ntr3
INTEGER :: nnaf           ! Counter for non-advected tracers
INTEGER :: nuni
INTEGER :: nsst           ! No of steady-state species
INTEGER :: ncsteps        ! No of chemical steps
INTEGER :: nit0=20        ! ftoy iterations with method=0
INTEGER :: nfphot
INTEGER :: jsubs
INTEGER :: method          ! chemistry integration method
INTEGER :: interval = imdi ! interval in timesteps between calls to chemistry
INTEGER :: nnfrp           ! Total number of fractional products
INTEGER :: nstst           ! No of steady state species
INTEGER :: nf
INTEGER :: ndepd           ! No of dry deposited species
INTEGER :: ndepw           ! No of wet deposited species
INTEGER :: nemit           ! No of emitted species
INTEGER :: ntro3, ntroh, ntrho2, ntrno
INTEGER :: nspo1d, nspo3p, nspo3, nspoh
INTEGER :: nspho2, nspno, nspn, nsph
INTEGER :: ih_o3, ih_h2o2, ih_so2, ih_hno3  ! index for soluble species
INTEGER :: ih_dms, ih_msia, ih_hobr !LER
INTEGER :: ih_o3_const                      ! index for soluble species
                                            !  as constant species
INTEGER :: ihso3_h2o2                       ! Index for HSO3- + H2O2(aq) reaction
INTEGER :: ihso3_o3                         ! Index for HSO3- + O3(aq) reaction
INTEGER :: iso3_o3                          ! Index for SO3-- + O3(aq) reaction
INTEGER :: ih2so4_hv                        ! Index for H2SO4 + hv reaction
INTEGER :: iso2_oh                          ! Index for SO2 + OH reaction
INTEGER :: ih2o2_oh                         ! Index for H2O2 + OH reaction
INTEGER :: ihno3_oh                         ! Index for HNO3 + OH reaction
INTEGER :: in2o5_h                          ! Index for N2O5 => HONO2 heterog. reaction
INTEGER :: iho2_h                           ! Index for HO2 + HO2 => H2O2 heterog. "
INTEGER :: idms_o3                          ! Index for aqueous-phase DMS + O3 reaction, LER, Mar2019
INTEGER :: imsia_o3                         ! Index for aqueous-phase MSIA + O3 reaction, LER, Mar2019
INTEGER :: imsi_o3                          ! Index for aqueous-phase MSI + O3 reaction, LER, Mar2019
INTEGER :: ihso3_hobr                       ! Index for aqueous-phase HSO3 + HOBr reaction, LER, Mar2019
INTEGER :: iso3_hobr                        ! Index for aqueous-phase SO3 + HOBr reaction, LER, Mar2019

LOGICAL :: lsvjac         ! Flag for saving jacobian if not recalculated
LOGICAL :: ljacx

! ltrig set to debug slow convergence systems
! shared between asad_spimpmjp and asad_spmjpdriv
LOGICAL :: ltrig

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_MOD'

! Variables which should be stored separately on each thread
!$OMP THREADPRIVATE(cdt, co3, deriv, dpd, dpw, ej, emr,                  &
!$OMP               f, fdot, fj, fpsc1, fpsc2, ftilde,                   &
!$OMP               interval, ipa, jsubs,                                &
!$OMP               lati, linfam, lsvjac, ltrig, modified_map,           &
!$OMP               ncsteps, p, pd, pmintnd, prk, prod,                  &
!$OMP               qa, ratio, rk,                                       &
!$OMP               sh2o, shno3, slos, spfj, sph2o, sphno3,              &
!$OMP               t, t300, tnd, wp, co2, y, ydot, za)

CONTAINS

! ######################################################################
SUBROUTINE asad_mod_init(n_points)

! To allocate and initialise ASAD arrays and variables

USE ukca_chem_defs_mod,       ONLY:  chch_t, chch_defs
USE UM_ParVars
USE Ereport_mod,              ONLY: ereport
USE ukca_option_mod,          ONLY:l_ukca, ukca_int_method

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

INTEGER :: i
INTEGER :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage

LOGICAL, SAVE :: firstcall=.TRUE.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MOD_INIT'

!$OMP THREADPRIVATE(firstcall)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! nullify prod and slos on firstcall to give DISSASSOCIATED attribute
IF (firstcall) THEN
  NULLIFY(prod)
  NULLIFY(slos)
  firstcall = .FALSE.
END IF

! variables and shared arrays should only be set and allocated on
! one thread
!$OMP SINGLE

! Fractional product parameters - set total number allowed to be
! total number of reactions * 2 as each fractional product has 4 potentials!!
jpfrpb  = (jpspb-2)*jpbk
jpfrpt  = (jpspt-2)*jptk
jpfrpj  = (jpspj-2)*jppj
jpfrph  = (jpsph-2)*jphk
jpfrpx  = jpfrpb + jpfrpt + jpfrpj + jpfrph  ! Total FPs



IF (.NOT. ALLOCATED(madvtr)) ALLOCATE(madvtr(jpspec))
IF (.NOT. ALLOCATED(majors)) ALLOCATE(majors(jpctr))
IF (.NOT. ALLOCATED(moffam)) ALLOCATE(moffam(jpspec))
IF (.NOT. ALLOCATED(nodd)) ALLOCATE(nodd(jpspec))
IF (.NOT. ALLOCATED(nltrf)) ALLOCATE(nltrf(jpctr))
IF (.NOT. ALLOCATED(nltr3)) ALLOCATE(nltr3(jpctr))
IF (.NOT. ALLOCATED(advt)) ALLOCATE(advt(jpctr))
IF (.NOT. ALLOCATED(nadvt)) ALLOCATE(nadvt(jpspec-jpctr))
IF (.NOT. ALLOCATED(family)) ALLOCATE(family(jpspec))
IF (.NOT. ALLOCATED(speci)) ALLOCATE(speci(jpspec))
IF (.NOT. ALLOCATED(ctype)) ALLOCATE(ctype(jpspec))
IF (.NOT. ALLOCATED(ipa2)) ALLOCATE(ipa2(jpctr))
IF (.NOT. ALLOCATED(nspi)) ALLOCATE(nspi(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nsspt)) ALLOCATE(nsspt(jpss))
IF (.NOT. ALLOCATED(nsspi)) ALLOCATE(nsspi(jpss,jpssr))
IF (.NOT. ALLOCATED(nssi)) ALLOCATE(nssi(jpss))
IF (.NOT. ALLOCATED(nssrt)) ALLOCATE(nssrt(jpss))
IF (.NOT. ALLOCATED(nssri)) ALLOCATE(nssri(jpss,jpssr))
IF (.NOT. ALLOCATED(nssrx)) ALLOCATE(nssrx(jpss,jpssr))
IF (.NOT. ALLOCATED(depvel)) ALLOCATE(depvel(jddept,jddepc,jpdd))
IF (.NOT. ALLOCATED(k298)) ALLOCATE(k298(jpdw))
IF (.NOT. ALLOCATED(dhr)) ALLOCATE(dhr(jpdw))
IF (.NOT. ALLOCATED(kd298)) ALLOCATE(kd298(jpdw,jpeq))
IF (.NOT. ALLOCATED(ddhr)) ALLOCATE(ddhr(jpdw,jpeq))
IF (.NOT. ALLOCATED(ldepd)) ALLOCATE(ldepd(jpspec))
IF (.NOT. ALLOCATED(ldepw)) ALLOCATE(ldepw(jpspec))
IF (.NOT. ALLOCATED(lemit)) ALLOCATE(lemit(jpspec))
IF (.NOT. ALLOCATED(nltrim)) ALLOCATE(nltrim(0:jpctr,3))
IF (.NOT. ALLOCATED(nlpdv)) ALLOCATE(nlpdv((jpspj-2)*jppj,2))
IF (.NOT. ALLOCATED(ab)) ALLOCATE(ab(jpbk,jpab))
IF (.NOT. ALLOCATED(at)) ALLOCATE(at(jptk,jpat))
IF (.NOT. ALLOCATED(aj)) ALLOCATE(aj(jppj,jpaj))
IF (.NOT. ALLOCATED(ah)) ALLOCATE(ah(jphk,jpah))
IF (.NOT. ALLOCATED(spb)) ALLOCATE(spb(jpbk+1,jpspb))
IF (.NOT. ALLOCATED(spt)) ALLOCATE(spt(jptk+1,jpspt))
IF (.NOT. ALLOCATED(spj)) ALLOCATE(spj(jppj+1,jpspj))
IF (.NOT. ALLOCATED(sph)) ALLOCATE(sph(jphk+1,jpsph))
IF (.NOT. ALLOCATED(frpb)) ALLOCATE(frpb(jpfrpb))
IF (.NOT. ALLOCATED(frpt)) ALLOCATE(frpt(jpfrpt))
IF (.NOT. ALLOCATED(frpj)) ALLOCATE(frpj(jpfrpj))
IF (.NOT. ALLOCATED(frph)) ALLOCATE(frph(jpfrph))
IF (.NOT. ALLOCATED(frpx)) ALLOCATE(frpx(jpfrpx))
IF (.NOT. ALLOCATED(ztabpd)) ALLOCATE(ztabpd(jpfrpd,2))
IF (.NOT. ALLOCATED(nfrpx)) ALLOCATE(nfrpx(jpnr))
IF (.NOT. ALLOCATED(ntabfp)) ALLOCATE(ntabfp(jpfrpx,3))
IF (.NOT. ALLOCATED(ntabpd)) ALLOCATE(ntabpd(jpfrpd,3))
IF (.NOT. ALLOCATED(npdfr)) ALLOCATE(npdfr(jpnr,2))
IF (.NOT. ALLOCATED(ngrp)) ALLOCATE(ngrp(2*jpspec,3))
IF (.NOT. ALLOCATED(njcgrp)) ALLOCATE(njcgrp(jpctr,3))
IF (.NOT. ALLOCATED(nprdx3)) ALLOCATE(nprdx3(3,(jpnr/(3*3))+3*3,2*jpspec))
IF (.NOT. ALLOCATED(nprdx2)) ALLOCATE(nprdx2(2,2*jpspec))
IF (.NOT. ALLOCATED(nprdx1)) ALLOCATE(nprdx1(2*jpspec))
IF (.NOT. ALLOCATED(njacx3)) ALLOCATE(njacx3(3,(jpnr/(3*3))+3*3,jpctr))
IF (.NOT. ALLOCATED(njacx2)) ALLOCATE(njacx2(2,jpctr))
IF (.NOT. ALLOCATED(njacx1)) ALLOCATE(njacx1(jpctr))
IF (.NOT. ALLOCATED(nmpjac)) ALLOCATE(nmpjac(jpctr))
IF (.NOT. ALLOCATED(npjac1)) ALLOCATE(npjac1(jppjac,jpctr))
IF (.NOT. ALLOCATED(nbrkx)) ALLOCATE(nbrkx(jpbk+1))
IF (.NOT. ALLOCATED(ntrkx)) ALLOCATE(ntrkx(jptk+1))
IF (.NOT. ALLOCATED(nprkx)) ALLOCATE(nprkx(jppj+1))
IF (.NOT. ALLOCATED(nhrkx)) ALLOCATE(nhrkx(jphk+1))
IF (.NOT. ALLOCATED(nlall)) ALLOCATE(nlall(jpspec))
IF (.NOT. ALLOCATED(nlstst)) ALLOCATE(nlstst(jpspec))
IF (.NOT. ALLOCATED(nlf)) ALLOCATE(nlf(jpspec))
IF (.NOT. ALLOCATED(nlmajmin)) ALLOCATE(nlmajmin(jpspec))
IF (.NOT. ALLOCATED(nldepd)) ALLOCATE(nldepd(jpspec))
IF (.NOT. ALLOCATED(nldepw)) ALLOCATE(nldepw(jpspec))
IF (.NOT. ALLOCATED(nlemit)) ALLOCATE(nlemit(jpspec))
IF (.NOT. ALLOCATED(nldepx)) ALLOCATE(nldepx(jpspec+6))
! nldepx(1:2) = start, end indices (in nldepx) of dry+wet deposited species
! nldepx(3:4) = start, end indices of species undergoing dry deposition only
! nldepx(5:6) = start, end indices of wet deposited species
IF (.NOT. ALLOCATED(njcoth)) ALLOCATE(njcoth(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmzjac)) ALLOCATE(nmzjac(jpctr))
IF (.NOT. ALLOCATED(nzjac1)) ALLOCATE(nzjac1(jpnr,jpctr))
IF (.NOT. ALLOCATED(njcoss)) ALLOCATE(njcoss(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmsjac)) ALLOCATE(nmsjac(jpctr))
IF (.NOT. ALLOCATED(nsjac1)) ALLOCATE(nsjac1(jpnr,jpctr))

! the following had save attribs.
IF (.NOT. ALLOCATED(ilcf)) ALLOCATE(ilcf(jpspec))
IF (.NOT. ALLOCATED(ilss)) ALLOCATE(ilss(jpspec))
IF (.NOT. ALLOCATED(ilct)) ALLOCATE(ilct(jpspec))
IF (.NOT. ALLOCATED(ilftr)) ALLOCATE(ilftr(jpspec))
IF (.NOT. ALLOCATED(ilft)) ALLOCATE(ilft(jpspec))
IF (.NOT. ALLOCATED(ilstmin)) ALLOCATE(ilstmin(jpspec))


! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_int_method

! Initialize variables that may be changed in cinit
nrsteps = 45            ! No of N-R steps, set > 50 to debug convergence failures
nitnr   = 10            ! Iterations in ftoy
nitfg   = 10            ! Max number of iterations in ftoy

! Initialise arrays
njcoth(:,:) = 0

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(nonzero_map_unordered))   &
      ALLOCATE(nonzero_map_unordered(jpctr, jpctr))
  IF (.NOT. ALLOCATED(nonzero_map))  ALLOCATE(nonzero_map(jpctr, jpctr))
  IF (.NOT. ALLOCATED(ro))  ALLOCATE(ro(jpctr))
END IF

! Find out which species are in steady state (for N-R solver)
o1d_in_ss = .FALSE.
o3p_in_ss = .FALSE.
n_in_ss   = .FALSE.
h_in_ss   = .FALSE.
DO i=1,jpspec
  IF (chch_defs(i)%speci=='O(1D)     ' .AND.                             &
      chch_defs(i)%ctype(1:2)==jpna) o1d_in_ss=.TRUE.
  IF (chch_defs(i)%speci=='O(3P)     ' .AND.                             &
      chch_defs(i)%ctype(1:2)==jpna) o3p_in_ss=.TRUE.
  IF (chch_defs(i)%speci=='N         ' .AND.                             &
      chch_defs(i)%ctype(1:2)==jpna) n_in_ss=.TRUE.
  IF (chch_defs(i)%speci=='H         ' .AND.                             &
      chch_defs(i)%ctype(1:2)==jpna) h_in_ss=.TRUE.
END DO

! Currently N-R solver assumes that O(1D) is in steady-state
IF (.NOT. o1d_in_ss .AND. (L_ukca_strat .OR. &   ! L_ukca_wachem .OR.  &
    L_ukca_stratcfc .OR. L_ukca_strattrop)) THEN
  cmessage=' O(1D) is not a Steady-State species'
  errcode = 1
  CALL ereport('ASAD_MOD_INIT',errcode,cmessage)
END IF

!$OMP END SINGLE

! The arrays which have copies on all threads

! pd is a TARGET
IF (.NOT. ALLOCATED(pd)) ALLOCATE(pd(n_points,2*jpspec))
! prod and slos are pointers
IF (.NOT. ASSOCIATED(prod)) ALLOCATE(prod(n_points,jpspec))
IF (.NOT. ASSOCIATED(slos)) ALLOCATE(slos(n_points,jpspec))

IF (.NOT. ALLOCATED(co3)) ALLOCATE(co3(n_points))
IF (.NOT. ALLOCATED(deriv)) ALLOCATE(deriv(n_points,4,4))
IF (.NOT. ALLOCATED(dpd)) ALLOCATE(dpd(n_points,jpspec))
IF (.NOT. ALLOCATED(dpw)) ALLOCATE(dpw(n_points,jpspec))
IF (.NOT. ALLOCATED(ej)) ALLOCATE(ej(n_points,jpctr))
IF (.NOT. ALLOCATED(emr)) ALLOCATE(emr(n_points,jpspec))
IF (.NOT. ALLOCATED(f)) ALLOCATE(f(n_points,jpctr))
IF (.NOT. ALLOCATED(fdot)) ALLOCATE(fdot(n_points,jpctr))
IF (.NOT. ALLOCATED(fj)) ALLOCATE(fj(n_points,jpctr,jpctr))
IF (.NOT. ALLOCATED(fpsc1)) ALLOCATE(fpsc1(n_points))
IF (.NOT. ALLOCATED(fpsc2)) ALLOCATE(fpsc2(n_points))
IF (.NOT. ALLOCATED(ftilde)) ALLOCATE(ftilde(n_points, jpctr))
IF (.NOT. ALLOCATED(ipa)) ALLOCATE(ipa(n_points,jpctr))
IF (.NOT. ALLOCATED(lati)) ALLOCATE(lati(n_points))
IF (.NOT. ALLOCATED(linfam)) ALLOCATE(linfam(n_points,0:jpctr))
IF (.NOT. ALLOCATED(p)) ALLOCATE(p(n_points))
IF (.NOT. ALLOCATED(pmintnd)) ALLOCATE(pmintnd(n_points))
IF (.NOT. ALLOCATED(prk)) ALLOCATE(prk(n_points,jpnr))
IF (.NOT. ALLOCATED(qa)) ALLOCATE(qa(n_points,jpspec))
IF (.NOT. ALLOCATED(ratio)) ALLOCATE(ratio(n_points,jpspec))
IF (.NOT. ALLOCATED(rk)) ALLOCATE(rk(n_points,jpnr))
IF (.NOT. ALLOCATED(sh2o)) ALLOCATE(sh2o(n_points))
IF (.NOT. ALLOCATED(shno3)) ALLOCATE(shno3(n_points))
IF (.NOT. ALLOCATED(sph2o)) ALLOCATE(sph2o(n_points))
IF (.NOT. ALLOCATED(sphno3)) ALLOCATE(sphno3(n_points))
IF (.NOT. ALLOCATED(t)) ALLOCATE(t(n_points))
IF (.NOT. ALLOCATED(t300)) ALLOCATE(t300(n_points))
IF (.NOT. ALLOCATED(tnd)) ALLOCATE(tnd(n_points))
IF (.NOT. ALLOCATED(wp)) ALLOCATE(wp(n_points))
IF (.NOT. ALLOCATED(co2)) ALLOCATE(co2(n_points))
IF (.NOT. ALLOCATED(y)) ALLOCATE(y(n_points,jpspec))
IF (.NOT. ALLOCATED(ydot)) ALLOCATE(ydot(n_points,jpspec))
IF (.NOT. ALLOCATED(za)) ALLOCATE(za(n_points))

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(modified_map))  ALLOCATE(modified_map(jpctr, jpctr))
  IF (.NOT. ALLOCATED(spfj))  ALLOCATE(spfj(n_points,spfjsize_max))
END IF

! EQUIVALENCE ( pd(1,1), prod(1,1) )
! EQUIVALENCE ( pd(1,jpspec+1), slos(1,1) )
prod => pd(:,1:jpspec)
slos => pd(:,jpspec+1:2*jpspec)

! Initialise arrays
deriv(:,:,:) = 1.0      ! Temp fix for deriv being uninitialised in first
                        ! solver iteration

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_init

! ######################################################################
SUBROUTINE asad_mod_final

! To deallocate ASAD arrays

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MOD_FINAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(nltrf)) DEALLOCATE(nltrf)
IF (ALLOCATED(nltr3)) DEALLOCATE(nltr3)
IF (ALLOCATED(advt)) DEALLOCATE(advt)
IF (ALLOCATED(nadvt)) DEALLOCATE(nadvt)
IF (ALLOCATED(family)) DEALLOCATE(family)
IF (ALLOCATED(speci)) DEALLOCATE(speci)
IF (ALLOCATED(ctype)) DEALLOCATE(ctype)
IF (ALLOCATED(ipa2)) DEALLOCATE(ipa2)
IF (ALLOCATED(nspi)) DEALLOCATE(nspi)
IF (ALLOCATED(nsspt)) DEALLOCATE(nsspt)
IF (ALLOCATED(nsspi)) DEALLOCATE(nsspi)
IF (ALLOCATED(nssi)) DEALLOCATE(nssi)
IF (ALLOCATED(nssrt)) DEALLOCATE(nssrt)
IF (ALLOCATED(nssri)) DEALLOCATE(nssri)
IF (ALLOCATED(nssrx)) DEALLOCATE(nssrx)
IF (ALLOCATED(depvel)) DEALLOCATE(depvel)
IF (ALLOCATED(k298)) DEALLOCATE(k298)
IF (ALLOCATED(dhr)) DEALLOCATE(dhr)
IF (ALLOCATED(kd298)) DEALLOCATE(kd298)
IF (ALLOCATED(ddhr)) DEALLOCATE(ddhr)
IF (ALLOCATED(ldepd)) DEALLOCATE(ldepd)
IF (ALLOCATED(ldepw)) DEALLOCATE(ldepw)
IF (ALLOCATED(lemit)) DEALLOCATE(lemit)
IF (ALLOCATED(nltrim)) DEALLOCATE(nltrim)
IF (ALLOCATED(nlpdv)) DEALLOCATE(nlpdv)
IF (ALLOCATED(ab)) DEALLOCATE(ab)
IF (ALLOCATED(at)) DEALLOCATE(at)
IF (ALLOCATED(aj)) DEALLOCATE(aj)
IF (ALLOCATED(ah)) DEALLOCATE(ah)
IF (ALLOCATED(spb)) DEALLOCATE(spb)
IF (ALLOCATED(spt)) DEALLOCATE(spt)
IF (ALLOCATED(spj)) DEALLOCATE(spj)
IF (ALLOCATED(sph)) DEALLOCATE(sph)
IF (ALLOCATED(frpb)) DEALLOCATE(frpb)
IF (ALLOCATED(frpt)) DEALLOCATE(frpt)
IF (ALLOCATED(frpj)) DEALLOCATE(frpj)
IF (ALLOCATED(frph)) DEALLOCATE(frph)
IF (ALLOCATED(frpx)) DEALLOCATE(frpx)
IF (ALLOCATED(ztabpd)) DEALLOCATE(ztabpd)
IF (ALLOCATED(nfrpx)) DEALLOCATE(nfrpx)
IF (ALLOCATED(ntabfp)) DEALLOCATE(ntabfp)
IF (ALLOCATED(ntabpd)) DEALLOCATE(ntabpd)
IF (ALLOCATED(npdfr)) DEALLOCATE(npdfr)
IF (ALLOCATED(ngrp)) DEALLOCATE(ngrp)
IF (ALLOCATED(njcgrp)) DEALLOCATE(njcgrp)
IF (ALLOCATED(nprdx3)) DEALLOCATE(nprdx3)
IF (ALLOCATED(nprdx2)) DEALLOCATE(nprdx2)
IF (ALLOCATED(nprdx1)) DEALLOCATE(nprdx1)
IF (ALLOCATED(njacx3)) DEALLOCATE(njacx3)
IF (ALLOCATED(njacx2)) DEALLOCATE(njacx2)
IF (ALLOCATED(njacx1)) DEALLOCATE(njacx1)
IF (ALLOCATED(nmpjac)) DEALLOCATE(nmpjac)
IF (ALLOCATED(npjac1)) DEALLOCATE(npjac1)
IF (ALLOCATED(nbrkx)) DEALLOCATE(nbrkx)
IF (ALLOCATED(ntrkx)) DEALLOCATE(ntrkx)
IF (ALLOCATED(nprkx)) DEALLOCATE(nprkx)
IF (ALLOCATED(nhrkx)) DEALLOCATE(nhrkx)
IF (ALLOCATED(nlall)) DEALLOCATE(nlall)
IF (ALLOCATED(nlstst)) DEALLOCATE(nlstst)
IF (ALLOCATED(nlf)) DEALLOCATE(nlf)
IF (ALLOCATED(nlmajmin)) DEALLOCATE(nlmajmin)
IF (ALLOCATED(nldepd)) DEALLOCATE(nldepd)
IF (ALLOCATED(nldepw)) DEALLOCATE(nldepw)
IF (ALLOCATED(nlemit)) DEALLOCATE(nlemit)
IF (ALLOCATED(nldepx)) DEALLOCATE(nldepx)
IF (ALLOCATED(njcoth)) DEALLOCATE(njcoth)
IF (ALLOCATED(nmzjac)) DEALLOCATE(nmzjac)
IF (ALLOCATED(nzjac1)) DEALLOCATE(nzjac1)
IF (ALLOCATED(njcoss)) DEALLOCATE(njcoss)
IF (ALLOCATED(nmsjac)) DEALLOCATE(nmsjac)
IF (ALLOCATED(nsjac1)) DEALLOCATE(nsjac1)

IF (method == int_method_NR) THEN       ! sparse vars
  IF (ALLOCATED(nonzero_map_unordered))   DEALLOCATE(nonzero_map_unordered)
  IF (ALLOCATED(nonzero_map))   DEALLOCATE(nonzero_map)
  IF (ALLOCATED(ro))   DEALLOCATE(ro)
END IF

!$OMP PARALLEL

! pd is a TARGET
IF (ALLOCATED(pd)) DEALLOCATE(pd)
! prod and slos are POINTERs
IF (ASSOCIATED(prod)) DEALLOCATE(prod)
IF (ASSOCIATED(slos)) DEALLOCATE(slos)

IF (ALLOCATED(co3)) DEALLOCATE(co3)
IF (ALLOCATED(deriv)) DEALLOCATE(deriv)
IF (ALLOCATED(dpd)) DEALLOCATE(dpd)
IF (ALLOCATED(dpw)) DEALLOCATE(dpw)
IF (ALLOCATED(ej)) DEALLOCATE(ej)
IF (ALLOCATED(emr)) DEALLOCATE(emr)
IF (ALLOCATED(f)) DEALLOCATE(f)
IF (ALLOCATED(fdot)) DEALLOCATE(fdot)
IF (ALLOCATED(fj)) DEALLOCATE(fj)
IF (ALLOCATED(fpsc1)) DEALLOCATE(fpsc1)
IF (ALLOCATED(fpsc2)) DEALLOCATE(fpsc2)
IF (ALLOCATED(ftilde)) DEALLOCATE(ftilde)
IF (ALLOCATED(ipa)) DEALLOCATE(ipa)
IF (ALLOCATED(lati)) DEALLOCATE(lati)
IF (ALLOCATED(linfam)) DEALLOCATE(linfam)
IF (ALLOCATED(p)) DEALLOCATE(p)
IF (ALLOCATED(pmintnd)) DEALLOCATE(pmintnd)
IF (ALLOCATED(prk)) DEALLOCATE(prk)
IF (ALLOCATED(qa)) DEALLOCATE(qa)
IF (ALLOCATED(ratio)) DEALLOCATE(ratio)
IF (ALLOCATED(rk)) DEALLOCATE(rk)
IF (ALLOCATED(sh2o)) DEALLOCATE(sh2o)
IF (ALLOCATED(shno3)) DEALLOCATE(shno3)
IF (ALLOCATED(sph2o)) DEALLOCATE(sph2o)
IF (ALLOCATED(sphno3)) DEALLOCATE(sphno3)
IF (ALLOCATED(t)) DEALLOCATE(t)
IF (ALLOCATED(t300)) DEALLOCATE(t300)
IF (ALLOCATED(tnd)) DEALLOCATE(tnd)
IF (ALLOCATED(wp)) DEALLOCATE(wp)
IF (ALLOCATED(co2)) DEALLOCATE(co2)
IF (ALLOCATED(y)) DEALLOCATE(y)
IF (ALLOCATED(ydot)) DEALLOCATE(ydot)
IF (ALLOCATED(za)) DEALLOCATE(za)

IF (method == int_method_NR) THEN       ! sparse vars
  IF (ALLOCATED(modified_map))   DEALLOCATE(modified_map)
  IF (ALLOCATED(spfj))   DEALLOCATE(spfj)
END IF

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_final

END MODULE
