! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates aerosol dry deposition and sedimentation.
!    Based on the parameterisation of Zhang et al (2001) which
!    uses the method in the model of Slinn (1982).
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_ddepaer_incl_sedi_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
                  ModuleName='UKCA_DDEPAER_INCL_SEDI_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_ddepaer_incl_sedi(nbox,nd,md,mdt,rhopar,znot,                 &
 dtc,wetdp,ustr,pmid,pupper,plower,t,surtp,seaice,                            &
 rhoa,mfpa,dvisc,bud_aer_mas,jlabove,ilscat,sedi_on,sm)

!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates aerosol dry deposition and sedimentation.
! Based on the parameterisation of Zhang et al (2001) which
! uses the method in the model of Slinn (1982).
!
! Sedimentation is done using a simple explicit discretization
! which should be adequate for this process-split method.
!
! Evaluate deposition velocity in lowest level as:
!
! V_dep = V_g + 1/(A_r + S_r)
!
! where V_dep is the deposition velocity
!       V_g   is the gravitational velocity = rho_p*Dp^2*g*CF/(18*DVISC)
!       CF    is the Cunningham slip correction
!       DVISC is the dynamic viscosity
!       A_r   is the aerodynamic resitance
!       S_r   is the surface resistance
!       Dp    is the particle diameter
!       rho_p is the particle density
!       g     is the gravitational acceleration
!
! Evaluate S_r=1/{ 3 * ustar * (EB + EIM + EIN) }
!
! following parameterization by Zhang et al (2001) where
!
! EB,EIM,EIN are collection efficiencies for Brownian diffusion,
! impaction and interception respectively.
!
! EB = Sc^-YR where Sc is the particle Schmidt number = nu/D
!                                where nu = kinematic viscosity of air
!                                      D =  particle diffusion coeff.
!
!  and YR is surface-dependent constant, values as in Table 3 (Zhang01)
!         0.50 over water          (Land use category 13-14)
!         0.56 over forest         (Land use category  1- 5)
!         0.54 over grass and ice  (Land use category  6,12)
!
! EIM = { St/(ALPHA+St) }^2
!
!    where St is the Stokes number = V_g * ustar^2 / DVISC  (z0<1mm)
!                                  = V_g * ustar   / (g*CR) (z0>1mm)
!
!                                    [smooth & rough flow regimes]
!
!      and ALPHA,CR are surface-dependent constant, values as Table 3:
!         ALPHA=100.0, CR=0.0 over water [only divide by CR for veg]
!         ALPHA= 50.0, CR=0.0 over ice   [only divide by CR for veg]
!         ALPHA=  1.0, CR=0.005 over grass
!         ALPHA=  1.2, CR=0.002 over forest
!
! EIN = 0.5*Dp/CR
!
! Evaluates drydep & sedimentation for number & mass using 0th & 3rd
! order moment specific coefficients for modal aerosol as in Appendix 4
! Binkowski & Shankar (1995) JGR, vol 100, no D12, pp. 26,191--26,209.
!
! Note --- in this routine, sedimentation is included at all levels.
!
! Inputs :
! ------
! NBOX      : Number of grid boxes
! ND        : Initial no. concentration of aerosol mode (ptcls/cc)
! MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
! MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
! RHOPAR    : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! ZNOT      : Roughness length (m)
! DTC       : Chemical timestep (s)
! WETDP     : Wet diameter for ptcl with dry diameter DRYDP (m)
! USTR      : Friction velocity(ms-1)
! PMID      : Centre level pressure (Pa)
! PUPPER    : Pressure at box upper interface (Pa)
! PLOWER    : Pressure at box lower interface (Pa)
! T         : Centre level temperature (K)
! SURTP     : Surface type [0=sea-surf,1=land-surf,2=above-surf]
! SEAICE    : Fraction of horizontal gridbox area containing seaice
! RHOA      : Air density (kg/m3)
! MFPA      : Mean free path of air (m)
! DVISC     : Dynamic viscosity of air (kg m-1 s-1)
! JLABOVE   : Index of box directly above this grid box
! ILSCAT    : Land-use category (1-9 based on UM landsurf types)
! SEDI_ON   : Switch for whether aerosol sedimentation is on/off
! SM        : Grid box mass of air (kg)
!
! Calls subroutine GETROUGH to read roughness length
!
! Outputs
! -------
! Updated particle number density ND (/cm3)
! Updated particle avg mass MD (molecules/particle)
! Updated Avg tot mass of aerosol ptcl in size mode MDT (particle^-1)
! BUD_AER_MAS: Aerosol mass budgets (mlcls/cc/tstep)
!
! Local Variables
! ---------------
! PS_AV_0    : 0th moment avg particle Schmidt Number
! PS_AV_3    : 3rd moment avg particle Schmidt Number
! KVISC      : Kinematic viscosity of air (m2 s-1)
! VGRAV_AV_0 : 0th moment avg grav. settling vel. (m s^-1)
! VGRAV_AV_3 : 3rd moment avg grav. settling vel. (m s^-1)
! VDEP_AV_0  : 0th moment avg deposition velocity (m s^-1)
! VDEP_AV_3  : 3rd moment avg deposition velocity (m s^-1)
! DCOEF_AV_0 : 0th moment avg particle diffusion coefficient(m2 s-1)
! DCOEF_AV_3 : 3rd moment avg particle diffusion coefficient(m2 s-1)
! SN_AV_0    : 0th moment avg Stokes number
! SN_AV_3    : 3rd moment avg Stokes number
! SR_AV_0    : 0th moment avg surface resistance
! SR_AV_3    : 3rd moment avg surface resistance
! EB_AV_0    : 0th moment avg collection eff. for Brownian diffusion
! EB_AV_3    : 3rd moment avg collection eff. for Brownian diffusion
! EIM_AV_0   : 0th moment avg collection eff. for impaction
! EIM_AV_3   : 3rd moment avg collection eff. for impaction
! EIN        : Collection eff. for interception
! AR         : Aerodynamic resistance
! MTOT       : Total aerosol mass conc [all cpts] (molecules/cm3)
! NEWN       : Updated number concentration (/cm3)
! DZ         : Ht difference between box vertical interfaces (m)
! DZMID      : Ht difference between box lower interface & mid-level (m)
! CR,Y,ALPHA: aerosol deposition coefficients
!     [vary with land category & input via DATA statements]
! CR        : Characteristic radius of collectors (m)
! Y         : Parameter for calculating Brownian diffusion
! ALPHA     : Parameter for calculating EIM
! MASK1     : Logical to define regions of domain for roughness categories
! MASK2     : Logical to define regions of domain for JLABOVE
! MASK3     : Logical to define regions with some aerosol in box or 1above
! MASK4     : Logical to define surface regions with aerosol in box/1above
! MASK_SMOO :Logical to define regsions over "smooth"  surface categories
! MASK_VEGE :Logical to define regsions over vegetated surface categories
! MASK_ABOVE_LIM: Logical to define boxes with VGRAV > allowed max value
! cfl_fraction: Safety margin by which to stay within CFL limits
! dtsedi    : Desired sedimentation timestep for each mode
! dtmode    : Actual sedimentation timestep for IMODE (integer
!             fraction of DTC to ensure nsedi is an integer)
! nsedi     : Number of sedimentation substeps for IMODE
! isedi     : Sedimentation substep interation index
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! GG        : Gravitational acceleration = 9.80665 ms^-2
! VKARMN    : Von Karman's constant = 0.4
! RA        : Dry air gas constant = 287.05 Jkg^-1 K^-1
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! MODE      : Defines which modes are set
! COMPONENT : Defines which cpts are allowed in each mode
! SIGMAG    : Geometric standard deviation of mode
! MM        : Molar masses of components (kg/mole)
! NUM_EPS   : Value of ND_0 below which do not recalculate MD (per cc)
!                                              or carry out process
! CP_SU     : Index of component in which H2SO4 is stored
! CP_BC     : Index of component in which BC is stored
! CP_OC     : Index of component in which 1st OC cpt is stored
! CP_CL     : Index of component in which NaCl is stored
! CP_SO     : Index of component in which 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
! References
! ----------
! Slinn, Atmos. En., 1982, 16, 1785-1794
! Zhang et al, Atmos. En., 2001, 35, 549-560
!
!----------------------------------------------------------------------

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE jules_surface_types_mod, ONLY: &
    npft,                          &
    ntype

USE parkind1,                ONLY: &
    jpim,                          &
    jprb

USE ukca_constants,          ONLY: &
    gg,                            &
    vkarmn,                        &
    ra

USE ukca_dcoff_par_av_k_mod, ONLY: &
    ukca_dcoff_par_av_k

USE ukca_mode_setup,         ONLY: &
    nmodes,                        &
    ncp,                           &
    mode,                          &
    component,                     &
    sigmag,                        &
    num_eps,                       &
    cp_su,                         &
    cp_bc,                         &
    cp_oc,                         &
    cp_cl,                         &
    cp_du,                         &
    cp_so

USE ukca_setup_indices

USE ukca_vgrav_av_k_mod,     ONLY: &
    ukca_vgrav_av_k

USE um_types,                ONLY: &
    log_small

USE umPrintMgr,              ONLY: &
    umPrint,                       &
    umMessage,                     &
    PrintStatus,                   &
    PrStatus_oper

USE vectlib_mod,             ONLY: &
    log_v

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook


IMPLICIT NONE

! .. Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: sedi_on
INTEGER, INTENT(IN) :: jlabove(nbox)
INTEGER, INTENT(IN) :: ilscat(nbox)
REAL, INTENT(IN)    :: rhopar(nbox,nmodes)
REAL, INTENT(IN)    :: znot(nbox)
REAL, INTENT(IN)    :: dtc
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: ustr(nbox)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: pupper(nbox)
REAL, INTENT(IN)    :: plower(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: surtp(nbox)
REAL, INTENT(IN)    :: seaice(nbox)
REAL, INTENT(IN)    :: rhoa(nbox)
REAL, INTENT(IN)    :: mfpa(nbox)
REAL, INTENT(IN)    :: dvisc(nbox)
REAL, INTENT(IN)    :: sm(nbox)
REAL, INTENT(INOUT) :: bud_aer_mas(nbox,0:nbudaer)
REAL, INTENT(INOUT) :: nd(nbox,nmodes)
REAL, INTENT(INOUT) :: md(nbox,nmodes,ncp)
REAL, INTENT(INOUT) :: mdt(nbox,nmodes)

! .. Local Variables
INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
INTEGER :: isedi
INTEGER :: nsedi
LOGICAL (KIND=log_small) :: mask2(nbox)
LOGICAL (KIND=log_small) :: mask3(nbox)
LOGICAL (KIND=log_small) :: mask4(nbox)
LOGICAL (KIND=log_small) :: masksurf(nbox)
LOGICAL (KIND=log_small) :: mask_smoo(nbox)
LOGICAL (KIND=log_small) :: mask_vege(nbox)
LOGICAL (KIND=log_small) :: mask_above_lim(nbox)
REAL    :: dtmode
REAL    :: ps_av_0(nbox)
REAL    :: ps_av_3(nbox)
REAL    :: kvisc(nbox)
REAL    :: vgrav_av_0(nbox)
REAL    :: vgrav_av_3(nbox)
REAL    :: vgrav_av_0_up(nbox)
REAL    :: vgrav_av_3_up(nbox)
REAL    :: dcoef_av_0(nbox)
REAL    :: dcoef_av_3(nbox)
REAL    :: eb_av_0(nbox)
REAL    :: eb_av_3(nbox)
REAL    :: eim_av_0(nbox)
REAL    :: eim_av_3(nbox)
REAL    :: ein(nbox)
REAL    :: sn_av_0(nbox)
REAL    :: sn_av_3(nbox)
REAL    :: ar(nbox)
REAL    :: sr_av_0(nbox)
REAL    :: sr_av_3(nbox)
REAL    :: vdep_av_0(nbox)
REAL    :: vdep_av_3(nbox)
REAL    :: dzmid(nbox)
REAL    :: dz(nbox)
REAL    :: dz_up(nbox)
REAL    :: t_up(nbox)
REAL    :: plo_up(nbox)
REAL    :: pup_up(nbox)
REAL    :: sm_up(nbox)
REAL    :: rhoa_up(nbox)
REAL    :: nd0(nbox,nmodes)
REAL    :: md0(nbox,nmodes,ncp)
REAL    :: nd0_up(nbox,nmodes)
REAL    :: md0_up(nbox,nmodes,ncp)
REAL    :: ndnew(nbox)
REAL    :: termin_1(nbox)
REAL    :: termin_2(nbox)
REAL    :: termin_n(nbox)
REAL    :: termout_n(nbox)
REAL    :: termin_m(nbox,ncp)
REAL    :: termout_m(nbox,ncp)
REAL    :: delnsedi(nbox)
REAL    :: delmsedi(nbox,ncp)
REAL    :: delmddep(nbox)
REAL    :: vgrav_lim(nbox)
REAL    :: vgrav_lim_up(nbox)
REAL, PARAMETER :: cfl_fraction = 0.9 
REAL, PARAMETER :: dtsedi(nmodes) =                        &
                (/ 3600.0, 3600.0, 1800.0, 900.0, 3600.0, 1800.0, 900.0 /)

CHARACTER(LEN=errormessagelength) :: cmessage


!--------------------------------------------------------------------

REAL, ALLOCATABLE :: alpha(:)
REAL, ALLOCATABLE :: cr(:)
REAL, ALLOCATABLE :: yr(:)

!! here changed land-use categories 1/2 (trees) to match values used
!! for forest in 5-LS-category representation (should match then).
!
! Now set to match 9 UM Land-use types (ILSCAT)      YR  ALPHA  CR
!! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.57  0.70 0.005
!! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.05 0.002
! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.56  1.00 0.005
! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.00 0.005
! 3=C3 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 4=C4 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 5=Shrub    (Zhang cat=10)         [shrub i. wood] 0.54  1.30 0.010
! 6=Urban    (Zhang cat=15)         [urban        ] 0.56  1.50 1.500
! 7=Water    (Zhang cat=13/14)      [inl wat/ocean] 0.50 100.0 0.000
! 8=Soil     (Zhang cat=8)          [desert       ] 0.54 50.00 0.000
! 9=Ice      (Zhang cat=12)         [ice cap/glac.] 0.54 50.00 0.000
!
! ILSCAT is now an input to this routine and passed in to UKCA_AERO_STEP
!--------------------------------------------------------------------
!!!      REAL, PARAMETER :: YR(5) = (/  0.5, 0.56, 0.54,0.0,0.54/)
!!!      REAL, PARAMETER :: CR(5) = (/ 0.00,5.E-3,2.E-3,0.0,0.00/)
!!!      REAL, PARAMETER :: ALPHA(5) = (/100.0,  1.0,  1.2,0.0,50.0/)
!
! Previously used code below to set ILSCAT based on roughness length
! as either 1 = water/sea, 2= forest, 3 = all other land
!           4 = desert (not used), 5 = sea ice
! and used values from Zhang et al (2001) for YR (gamma in paper),
!                                             ALPHA
!                                             CR (A in paper)
! BELOW IS OLD TREATMENT                            YR   ALPHA    CR
! 1=water  (Zhang cat=13/14)[inland water/ocean  ] 0.50  100.0   0.000
! 2=forest (Zhang cat=1    )[evergreen needleleaf] 0.56    1.0   0.005
! 3=o land (Zhang cat=6/7  )[grass/crops         ] 0.54    1.2   0.002
! 4=desert (Zhang cat=8    )[desert              ] 0.54   50.0   0.000
! 5=seaice (Zhang cat=12   )[ice cap & glacier   ] 0.54   50.0   0.000

REAL :: p_ratio(nbox)
REAL :: p_ratio_out(nbox)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='UKCA_DDEPAER_INCL_SEDI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(alpha(ntype))
ALLOCATE(cr(ntype))
ALLOCATE(yr(ntype))

SELECT CASE (NTYPE)
CASE (9)
  yr    = (/0.56,0.56,0.54,0.54,0.54,0.56,0.50,0.54,0.54/)
  cr    = (/5.0e-3,5.0e-3,2.0e-3,2.0e-3,0.01,1.5e0,0.0e0,0.0e0,0.0e0/)
  alpha = (/1.00,1.00,1.20,1.02,1.30,1.50,100.0,50.0,50.0/)
CASE (13, 17, 27)
  yr(1:6)    = (/   0.56,   0.58,   0.58,   0.56,   0.56,   0.54/)
  cr(1:6)    = (/ 7.0e-3, 5.0e-3, 5.0e-3, 3.2e-3, 2.0e-3, 3.2e-3/)
  alpha(1:6) = (/   0.80,   0.60,   0.60,   1.10,   1.00,   1.20/)
CASE DEFAULT
  WRITE(umMessage,'(A)') 'NTYPE must equal 9 or 13 or 17 or 27'
  CALL umPrint(umMessage,src=RoutineName)
  cmessage = 'Unexpected value of NTYPE'
  errcode = 1000 + ntype
  CALL ereport(RoutineName,errcode,cmessage)
END SELECT

SELECT CASE (NTYPE)
CASE (13)
  yr(7:13)    = (/   0.54,   0.54,   0.54,   0.56,   0.50,   0.54,  0.54/)
  cr(7:13)    = (/ 3.2e-3, 1.0e-2, 1.0e-2, 1.0e-2,  0.0e0,  0.0e0, 0.0e0/)
  alpha(7:13) = (/   1.20,   1.30,   1.30,   1.50,  100.0,   50.0,  50.0/)
CASE (17, 27)
  yr(7:17)    = (/   0.54,   0.54,   0.54,   0.54,   0.54,   0.54,            &
                     0.54,   0.56,   0.50,   0.54,   0.54/)
  cr(7:17)    = (/ 3.2e-3, 3.2e-3, 3.2e-3, 3.2e-3, 3.2e-3, 1.0e-2,            &
                   1.0e-2, 1.0e-2,  0.0e0,  0.0e0,  0.0e0/)
  alpha(7:17) = (/   1.20,   1.20,   1.20,   1.20,   1.20,   1.30,            &
                     1.30,   1.50,  100.0,   50.0,   50.0/)
END SELECT

IF ( ntype == 27 ) THEN
  yr(18:27)    = (/  0.54,  0.54,  0.54,  0.54,  0.54,                        &
                     0.54,  0.54,  0.54,  0.54,  0.54 /)
  cr(18:27)    = (/ 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0,                        &
                    0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0 /)
  alpha(18:27) = (/  50.0,  50.0,  50.0,  50.0,  50.0,                        &
                     50.0,  50.0,  50.0,  50.0,  50.0 /)
END IF

! Find out what category (water,forest,grass,desert)
! based on roughness length znot (desert not used at present)
! This should be updated in later version to read land type
! category directly. Desert not represented here.
!
!!!! water/sea - z0<0.001m
!!!      MASK1(:)=(ZNOT(:) < 1.0E-3) ! water/sea
!!!      WHERE(MASK1(:)) ILSCAT(:)=1
!!!
!!!! forests - z0>0.1m
!!!      MASK1(:)=(ZNOT(:) > 1.0E-1) ! forest
!!!      WHERE(MASK1(:)) ILSCAT(:)=2
!!!
!!!! all other lands, grass 0.001<z0<0.1m
!!!      MASK1(:)=((ZNOT(:) >= 1.0E-3).AND.(ZNOT(:) <= 1.0E-1)) ! grass
!!!      WHERE(MASK1(:)) ILSCAT(:)=3
!!!
!!!! If sea ice covers > 50% of sea surface, treat as sea ice
!!!      MASK1(:)=(SEAICE(:) > 0.5) ! seaice
!!!      WHERE(MASK1(:)) ILSCAT(:)=5
!--------------------------------------------------------------------
!
mask2(:)=(jlabove(:) > 0) ! where using JLABOVE

WHERE (mask2(:))
  t_up(:)=t(jlabove(:))
  plo_up(:)=plower(jlabove(:))
  pup_up(:)=pupper(jlabove(:))
  sm_up(:)=sm(jlabove(:))
  rhoa_up(:)=rhoa(jlabove(:))
ELSEWHERE
  t_up(:)=t(:)
  plo_up(:)=plower(:)
  pup_up(:)=pupper(:)
  sm_up(:)=sm(:)
  rhoa_up(:)=rhoa(:)
END WHERE

p_ratio(:)=plower(:)/pupper(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dz(:)=ra*t(:)* p_ratio_out(:)/gg

p_ratio(:)=plower(:)/pmid(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dzmid(:)=ra*t(:)* p_ratio_out(:)/gg

p_ratio(:)=plo_up(:)/pup_up(:)
CALL log_v(nbox,p_ratio,p_ratio_out)
dz_up(:)=ra*t_up(:)* p_ratio_out(:)/gg

kvisc(:)=dvisc(:)/rhoa(:) ! Calculate kinematic viscosity of air

masksurf(:)=(surtp(:) < 2.0) ! create mask for boxes at surface.

! .. Calculate aerodynamic resistance (only in gridboxes at surface)
WHERE (masksurf(:)) ar(:)=LOG(dzmid(:)/znot(:))/(vkarmn*ustr(:))

DO imode=1,nmodes
  IF (mode(imode)) THEN

     ! .. Compute num iters from timestep dtmode, then recompute dtmode to
     ! .. force it to be an exact divisor of DTC.  Ensure dtmode is not larger
     ! .. than DTC.  Round nsedi up (instead of down or nearest) to ensure
     ! .. timestep is no longer than dtsedi(IMODE)
    dtmode = MIN(dtsedi(imode), dtc)
    nsedi = CEILING(dtc/dtmode)
    dtmode = dtc / REAL(nsedi)

    DO isedi = 1, nsedi

      ! .. First copy original values of ND,MD to ND0,MD0
      nd0(:,imode)=nd(:,imode)
      md0(:,imode,:)=md(:,imode,:)

      WHERE (mask2(:))
        nd0_up(:,imode)=nd0(jlabove(:),imode)
      ELSEWHERE
        nd0_up(:,imode)=0.0
      END WHERE
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask2(:))
            md0_up(:,imode,icp)=md0(jlabove(:),imode,icp)
          ELSEWHERE
            md0_up(:,imode,icp)=md0(:,imode,icp)
          END WHERE
        END IF
      END DO

      ! .. Calculate 0th moment avg. grav. settling velocities
      CALL ukca_vgrav_av_k(nbox,0,wetdp(:,imode),sigmag(imode),               &
             dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_0(:))
      ! .. Calculate 3rd moment avg. grav. settling velocities
      CALL ukca_vgrav_av_k(nbox,3,wetdp(:,imode),sigmag(imode),               &
             dvisc(:),mfpa(:),rhopar(:,imode),vgrav_av_3(:))
      !
      ! .. Store values from box above so that can get flux into box
      WHERE (mask2(:))
        vgrav_av_0_up(:)=vgrav_av_0(jlabove(:))
        vgrav_av_3_up(:)=vgrav_av_3(jlabove(:))
      ELSEWHERE
        vgrav_av_0_up(:)=0.0
        vgrav_av_3_up(:)=0.0
      END WHERE

      ! .. If sedimentation switched off, set all calculated VGRAV to zero
      IF (sedi_on == 0) THEN
        vgrav_av_0(:)=0.0
        vgrav_av_3(:)=0.0
        vgrav_av_0_up(:)=0.0
        vgrav_av_3_up(:)=0.0
      END IF

      !       Calculate 0th moment avg particle diffusion coeffs
      CALL ukca_dcoff_par_av_k(nbox,0,wetdp(:,imode),sigmag(imode),           &
                t(:),dvisc(:),mfpa(:),dcoef_av_0(:))

      !       Calculate 3rd moment avg particle diffusion coeffs
      CALL ukca_dcoff_par_av_k(nbox,3,wetdp(:,imode),sigmag(imode),           &
                t(:),dvisc(:),mfpa(:),dcoef_av_3(:))

      ! .. only calculate Schmidt number and collection eff at surface
      WHERE (masksurf(:))

        !       Calculate 0th and 3rd moment avg. particle Schmidt number
        ps_av_0(:)=kvisc(:)/dcoef_av_0(:)
        ps_av_3(:)=kvisc(:)/dcoef_av_3(:)
        !       Calculate particle collection efficiencies
        !       -- For Brownian Diffusion
        eb_av_0(:)=ps_av_0(:)**(-yr(ilscat(:)))
        eb_av_3(:)=ps_av_3(:)**(-yr(ilscat(:)))
        !
      END WHERE

      ! In new version, using 9 UM landsurf types,
      ! Set smooth surfaces to be water (7), soil (8) or ice (9)
      ! All other surfaces are vegetated (have CR>0)
      mask_smoo=( (ilscat(:) >= 7) .AND. (ilscat(:) <= 9) )
      mask_vege= .NOT. mask_smoo(:)

      !!! Below is as in previous version using 5 landsurf types
      !!        MASK_SMOO=( (ILSCAT(:) == 1).OR.(ILSCAT(:) == 5) )
      !!        MASK_VEGE=( (ILSCAT(:) == 2).OR.(ILSCAT(:) == 3) )

      ! .. only calculate Stokes number at surface
      WHERE (mask_smoo(:) .AND. masksurf(:))
        !        Calculate stokes number for smooth surfaces
        sn_av_0(:)=vgrav_av_0(:)*ustr(:)*ustr(:)/dvisc(:)
        sn_av_3(:)=vgrav_av_3(:)*ustr(:)*ustr(:)/dvisc(:)
      END WHERE
      WHERE (mask_vege(:) .AND. masksurf(:))
        !        Calculate stokes number for vegetated surfcaes
        sn_av_0(:)=vgrav_av_0(:)*ustr(:)/(gg*cr(ilscat(:)))
        sn_av_3(:)=vgrav_av_3(:)*ustr(:)/(gg*cr(ilscat(:)))
      END WHERE

      ! .. only calculate impaction collection efficiency at surface
      WHERE (masksurf(:))

        !       -- For Impaction
        eim_av_0(:)=(sn_av_0(:)/(alpha(ilscat(:))+sn_av_0(:)))**2
        eim_av_3(:)=(sn_av_3(:)/(alpha(ilscat(:))+sn_av_3(:)))**2

      END WHERE

      ! .. only calculate interception collection eff (smooth) at surface
      WHERE (mask_smoo(:) .AND. masksurf(:))

        !       -- For Interception (smooth surfaces)
        ein(:)=0.0

      END WHERE

      ! .. only calculate interception collection eff (vegetd) at surface
      WHERE (mask_vege(:) .AND. masksurf(:))

        !       -- For Interception (vegetd surfaces)
        ein(:)=0.5*(wetdp(:,imode)*wetdp(:,imode)                             &
                             /cr(ilscat(:))/cr(ilscat(:)))

      END WHERE

      ! .. only calculate surface resistance and deposition vel. at surface
      ! .. this section also increases VGRAV due to dry dep (at surface)
      WHERE (masksurf(:))

        !  Calculate surface resistance
        sr_av_0(:)=1.0/(3.0*ustr(:)*(eb_av_0(:)+eim_av_0(:)+ein(:)))
        sr_av_3(:)=1.0/(3.0*ustr(:)*(eb_av_3(:)+eim_av_3(:)+ein(:)))

        !  Calculate deposition velocity
        vdep_av_0(:)=vgrav_av_0(:)+1.0/(ar(:)+sr_av_0(:))
        vdep_av_3(:)=vgrav_av_3(:)+1.0/(ar(:)+sr_av_3(:))

        !  Set gravitational velocity to deposition velocity if in lowest box
        vgrav_av_0(:)=vdep_av_0(:)
        vgrav_av_3(:)=vdep_av_3(:)

        !  VGRAV_AV_UP never at surface so no need to set to dep. vel.

      END WHERE

      ! .. limit 0&3 grav. settling vel so only falls 
      ! .. cfl_fraction*box max [numerical]
      vgrav_lim(:)=cfl_fraction*dz(:)/dtmode

      mask_above_lim(:)=(vgrav_av_0(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_0(:)=vgrav_lim(:)

      mask_above_lim(:)=(vgrav_av_3(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_3(:)=vgrav_lim(:)

      ! .. limit 0&3 grav. settling vel so only falls 
      ! .. cfl_fraction*box max [box above]
      vgrav_lim(:)=cfl_fraction*dz_up(:)/dtmode

      mask_above_lim(:)=(vgrav_av_0_up(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_0_up(:)=vgrav_lim(:)

      mask_above_lim(:)=(vgrav_av_3_up(:) > vgrav_lim(:))
      WHERE (mask_above_lim(:)) vgrav_av_3_up(:)=vgrav_lim(:)

      ! .. Calculate sedimenting in term (number) using 0th moment VGRAV_AV_0
      WHERE (mask2(:)) ! where not top model level
        termin_1(:)=nd0_up(:,imode)/dz_up(:)
        termin_2(:)=sm_up(:)*rhoa(:)/sm(:)/rhoa_up(:)
        termin_n(:)=termin_1(:)*termin_2(:)*dtmode*vgrav_av_0_up(:)
      ELSEWHERE
        !    If top model level TERMIN_N(:)=0.0
        termin_n(:)=0.0
      END WHERE

      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask2(:))
            termin_1(:)=nd0_up(:,imode)*md0_up(:,imode,icp)/dz_up(:)
            termin_m(:,icp)=termin_1(:)*termin_2(:)*dtmode*vgrav_av_3_up(:)
          ELSEWHERE
            !    If top model level TERMIN_M(:,ICP)=0.0
            termin_m(:,icp)=0.0
          END WHERE
        END IF
      END DO

      ! .. Calculate sedimenting out term (number) using 0th moment VGRAV_AV_0
      termout_n(:)=(nd0(:,imode)/dz(:))*dtmode*vgrav_av_0(:)

      ! .. Calculate sedimenting in term (mass  ) using 3rd moment VGRAV_AV_3
      ! .. Calculate sedimenting out term (mass  ) using 3rd moment VGRAV_AV_3
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          termout_m(:,icp)=(nd0(:,imode)*md0(:,imode,icp)/dz(:))              &
                        *dtmode*vgrav_av_3(:)
        END IF
      END DO

      ! .. below calculates net change in number and mass concentration
      delnsedi(:)=-(termin_n(:)-termout_n(:))
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          delmsedi(:,icp)=-(termin_m(:,icp)-termout_m(:,icp))
        END IF
      END DO

      ! .. below masks boxes with some updates to apply
      mask3(:)=( (nd0   (:,imode) > num_eps(imode)) .OR.                      &
                 (nd0_up(:,imode) > num_eps(imode)) )
      mask4(:)=mask3(:) .AND. masksurf(:) ! some ptcls & also at surface

      ! .. below sets NDNEW to new value and re-sets MDT in boxes to update
      WHERE (mask3(:)) ! only do where some particles in/out
        ndnew(:)=nd0(:,imode)-delnsedi(:)
        mdt(:,imode)=0.0
      END WHERE

      ! .. below updates component masses and MDT for each mode
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask3(:)) ! only do where some particles in/out
            md(:,imode,icp)=                                                  &
              (md0(:,imode,icp)*nd0(:,imode)-delmsedi(:,icp))/ndnew(:)
            mdt(:,imode)=mdt(:,imode)+md(:,imode,icp)
          END WHERE
        END IF ! IF COMPONENT(ICP)
      END DO ! loop over cpts

      ! .. below updates number concentration to NDNEW
      ! ..  (only do where some particles in/out)
      WHERE (mask3(:)) nd(:,imode)=ndnew(:)

      ! .. below stores ddep/sedi fluxes to BUD_AER_MAS
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          WHERE (mask4(:)) delmddep(:)=termout_m(:,icp)
          IF (icp == cp_su) THEN
            IF ((imode == 1) .AND. (nmasddepsunucsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsunucsol)=                &
                           bud_aer_mas(:,nmasddepsunucsol)+delmddep(:)
            IF ((imode == 2) .AND. (nmasddepsuaitsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsuaitsol)=                &
                           bud_aer_mas(:,nmasddepsuaitsol)+delmddep(:)
            IF ((imode == 3) .AND. (nmasddepsuaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsuaccsol)=                &
                           bud_aer_mas(:,nmasddepsuaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepsucorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsucorsol)=                &
                           bud_aer_mas(:,nmasddepsucorsol)+delmddep(:)
          END IF
          IF (icp == cp_bc) THEN
            IF ((imode == 2) .AND. (nmasddepbcaitsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaitsol)=                &
                        bud_aer_mas(:,nmasddepbcaitsol)+delmddep(:)
            IF ((imode == 3) .AND. (nmasddepbcaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaccsol)=                &
                        bud_aer_mas(:,nmasddepbcaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepbccorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbccorsol)=                &
                        bud_aer_mas(:,nmasddepbccorsol)+delmddep(:)
            IF ((imode == 5) .AND. (nmasddepbcaitins > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepbcaitins)=                &
                        bud_aer_mas(:,nmasddepbcaitins)+delmddep(:)
          END IF
          IF (icp == cp_oc) THEN
            IF ((imode == 1) .AND. (nmasddepocnucsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocnucsol)=                &
                        bud_aer_mas(:,nmasddepocnucsol)+delmddep(:)
            IF ((imode == 2) .AND. (nmasddepocaitsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaitsol)=                &
                        bud_aer_mas(:,nmasddepocaitsol)+delmddep(:)
            IF ((imode == 3) .AND. (nmasddepocaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaccsol)=                &
                        bud_aer_mas(:,nmasddepocaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepoccorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepoccorsol)=                &
                        bud_aer_mas(:,nmasddepoccorsol)+delmddep(:)
            IF ((imode == 5) .AND. (nmasddepocaitins > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepocaitins)=                &
                        bud_aer_mas(:,nmasddepocaitins)+delmddep(:)
          END IF
          IF (icp == cp_cl) THEN
            IF ((imode == 3) .AND. (nmasddepssaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepssaccsol)=                &
                        bud_aer_mas(:,nmasddepssaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepsscorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsscorsol)=                &
                        bud_aer_mas(:,nmasddepsscorsol)+delmddep(:)
          END IF
          IF (icp == cp_so) THEN
            IF ((imode == 1) .AND. (nmasddepsonucsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsonucsol)=                &
                        bud_aer_mas(:,nmasddepsonucsol)+delmddep(:)
            IF ((imode == 2) .AND. (nmasddepsoaitsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsoaitsol)=                &
                        bud_aer_mas(:,nmasddepsoaitsol)+delmddep(:)
            IF ((imode == 3) .AND. (nmasddepsoaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsoaccsol)=                &
                        bud_aer_mas(:,nmasddepsoaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepsocorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepsocorsol)=                &
                        bud_aer_mas(:,nmasddepsocorsol)+delmddep(:)
          END IF
          IF (icp == cp_du) THEN
            IF ((imode == 3) .AND. (nmasddepduaccsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepduaccsol)=                &
                        bud_aer_mas(:,nmasddepduaccsol)+delmddep(:)
            IF ((imode == 4) .AND. (nmasddepducorsol > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepducorsol)=                &
                        bud_aer_mas(:,nmasddepducorsol)+delmddep(:)
            IF ((imode == 6) .AND. (nmasddepduaccins > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepduaccins)=                &
                        bud_aer_mas(:,nmasddepduaccins)+delmddep(:)
            IF ((imode == 7) .AND. (nmasddepducorins > 0))                    &
             WHERE (mask4(:)) bud_aer_mas(:,nmasddepducorins)=                &
                        bud_aer_mas(:,nmasddepducorins)+delmddep(:)
          END IF
        END IF ! if component present in mode
      END DO ! loop over components
    END DO ! loop over sedimentation timesteps
  END IF ! IF MODE(IMODE)
END DO ! loop over modes

DEALLOCATE(yr)
DEALLOCATE(cr)
DEALLOCATE(alpha)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_ddepaer_incl_sedi

END MODULE ukca_ddepaer_incl_sedi_mod
