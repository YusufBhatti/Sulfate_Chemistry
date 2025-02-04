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
!    Calculates condensation of condensable cpt vapours onto pre-existing
!    aerosol particles. Includes switch for using either Fuchs (1964) or
!    modified Fuchs and Sutugin (1971) calculation of CC.
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
! Subroutine Interface:
MODULE ukca_conden_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CONDEN_MOD'

CONTAINS

SUBROUTINE ukca_conden(nbox,gc,nd,md,mdt,                         &
 dtz,drydp,wetdp,tsqrt,rhoa,airdm3,delgc_cond,                    &
 ifuchs,ageterm1,bud_aer_mas,s_cond_s,pmid,t,s,aird,              &
 idcmfp,icondiam)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates condensation of condensable cpt vapours onto pre-existing
! aerosol particles. Includes switch for using either Fuchs (1964) or
! modified Fuchs and Sutugin (1971) calculation of CC.
!
! Parameters
! ----------
! SE_SOL : Sticking efficiency for   soluble modes [set to 1.0 as in M7]
! SE_INS : Sticking efficiency for insoluble modes [set to 0.3 as in M7]
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! GC       : Condensable cpt number density (molecules cm-3)
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DTZ      : Time Step for nucl/cond competition (s)
! DRYDP    : Avg dry diameter for each aerosol mode (m)
! WETDP    : Avg wet diameter for each aerosol mode (m)
! TSQRT    : Square-root of mid-level temperature (K)
! RHOA     : Air density (kg/m3)
! AIRDM3   : Number density of air (per m3)
! IFUCHS   : Switch for Fuchs (1964) or Fuchs-Sutugin (1971) for CC
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
! S_COND_S :  Condensation sink
! PMID        : Centre level pressure (Pa)
! T           : Centre level temperature (K)
! S        : Specific humidity (kg/kg)
! AIRD     : Number density of air (per cm3)
! IDCMFP      : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
! ICONDIAM : Switch : wet diam in UKCA_CONDEN (1=g.mean,2=condiam.)
!
! Outputs:
! -------
! MD    : Avg aerosol ptcl mass in mode (molecules per ptcl)
! GC    : Condensable cpt number density (molecules per cc)
! AGETERM1: stores mass of soluble material which has condensed onto
!           each of the insoluble modes for use in UKCA_AGEING
!           (molecules per cc)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DELGC_COND : Change in vapour conc due to cond (molecules cpt/cm3)
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
!
! Local variables:
! ---------------
! DMOL    : Molecular diameter of condensable cpt (m)
! CC      : Conden. coeff. for condensable cpt onto particle (m^3s^-1)
! RP      : Radius of aerosol particle (m)
! SE      : Sticking efficiency (accomodation coeff)
! SE_SOL  : Sticking efficiency (accomodation coeff) for soluble mode
! SE_INS  : Sticking efficiency (accomodation coeff) for insoluble mode
! MMCG    : Molar mass of condensing gas (kg/mole)
! NC      : Product of number conc and condensation coefficient
! SUMNC   : Sum of NC over all modes
! DELTAMS : Mass of condensing gas taken up by this   soluble mode
! DELTAMI : Mass of condensing gas taken up by this insoluble mode
! DELTAM  : Mass of condensing gas taken up by both modes (-->soluble)
!    n.b. DELTAMS,DELTAMI,DELTAM all in molecules per cc.
! MASK1-3 : Logical array to define regions of domain to work on
! VPKEL    : vapour pressure of H2SO4 includes kelvin effect (Pa)
! GC_EQ    : equilibrium vapour pressure (mol cm-3)
! PTCLCONC : Product of MDH2SO4 and ND - aerosol conc in mol cm-3
! SURFTEN : Surface tension of H2SO4 aerosol = (0.0728 J m-2)
!
! References
! ----------
! Gong et al, JGR, 108(D1), 4007, doi:10.1029/2001JD002002, 2003.
! Raes et al, J. Aerosol Sci., 23 (7), pp. 759--771, 1992.
! Fuchs & Sutugin, Topics in aerosol research, 1971.
! Fuchs, "Mechanics of aerosols", Pergamon Press, 1964.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! RA       : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RR       : Universal gas constant = 8.314 Jmol^-1 K^-1
! AVC      : Avogadros constant (mol-1)
! CONC_EPS : Threshold for condensable conc (molecules per cc)
! ZBOLTZ   : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Defines which modes are set
! COMPONENT: Defines which cpts are allowed in each mode
! CONDENSABLE : Logical variable defining which cpts are condensable
! MODESOL  : Defines whether the mode is soluble or not (=1 or 0)
! MM       : Molar masses of components (kg/mole)
! DIMEN    : Molecular diamters of condensable components (m)
! NUM_EPS  : Value of NEWN below which do not recalculate MD (per cc)
!                                             or carry out process
! CP_SU    : Index of component in which sulfate is stored
! CP_OC    : Index of component in which 1st OC cpt is stored
! CP_SO    : Index of component in which 2nd OC cpt is stored
! SIGMAG   : Geometric standard deviation of mode
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! MH2SO4   : Index of MM_GAS, WTRATC and S0G for H2SO4
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
! DIMEN    : Molecular diamters of condensable components (m)
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_constants,       ONLY: ra, rr, avc, conc_eps, zboltz
USE ukca_mode_setup,      ONLY: nmodes, ncp, mode, component,     &
                                sigmag, modesol, mm, num_eps,     &
                                cp_su, cp_oc, cp_so
USE ukca_setup_indices
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE um_types,             ONLY: log_small
USE umPrintMgr
USE ukca_cond_coff_v_mod, ONLY: ukca_cond_coff_v
IMPLICIT NONE
!
! Subroutine interface:
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: ifuchs
INTEGER, INTENT(IN) :: idcmfp
INTEGER, INTENT(IN) :: icondiam
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(IN)    :: tsqrt(nbox)
REAL, INTENT(IN)    :: rhoa(nbox)
REAL, INTENT(IN)    :: airdm3(nbox)
REAL, INTENT(IN)    :: dtz
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: drydp(nbox,nmodes)
REAL, INTENT(IN)    :: pmid(nbox)
REAL, INTENT(IN)    :: t(nbox)
REAL, INTENT(IN)    :: s(nbox)
REAL, INTENT(IN)    :: aird(nbox)
REAL, INTENT(INOUT) :: md(nbox,nmodes,ncp)
REAL, INTENT(INOUT) :: mdt(nbox,nmodes)
REAL, INTENT(INOUT) :: gc(nbox,nchemg)
REAL, INTENT(INOUT) :: bud_aer_mas(nbox,0:nbudaer)
REAL, INTENT(OUT)   :: delgc_cond(nbox,nchemg)
REAL, INTENT(OUT)   :: ageterm1(nbox,3,nchemg)
REAL, INTENT(OUT)   :: s_cond_s(nbox)
!
! Local variables
INTEGER :: icp
INTEGER :: imode
INTEGER :: jv
INTEGER :: jl
INTEGER :: ierr        ! Error code
LOGICAL (KIND=log_small) :: mask1(nbox)
LOGICAL (KIND=log_small) :: mask2(nbox)
LOGICAL (KIND=log_small) :: mask3(nbox)
LOGICAL (KIND=log_small) :: mask3i(nbox)
REAL    :: dmol
REAL    :: mmcg
REAL    :: cc(nbox)
REAL    :: rp(nbox)
REAL    :: sumnc(nbox)
REAL    :: nc(nbox,nmodes)
REAL    :: deltam(nbox)
REAL    :: deltams(nbox)
REAL    :: deltami(nbox)
REAL    :: se
REAL    :: y2
REAL    :: aa
REAL    :: aa_modes(nmodes)   ! AA specific to each mode
REAL, PARAMETER :: se_sol=1.0
!!      REAL, PARAMETER :: SE_INS=0.3
REAL, PARAMETER :: se_ins=1.0
!
REAL, PARAMETER :: surften=0.0728
!
REAL :: sinkarr(nbox)
REAL :: difvol
!
!! added below for stratosphere
REAL :: vpkel(nbox)
REAL :: difkel(nbox)
REAL :: ptclconc(nbox)
REAL :: delgc_evap(nbox,nchemg)
REAL :: delgc_evap_mode(nbox)
REAL :: frac_evap_mode(nbox,nmodes)
REAL :: wts(nbox)
REAL :: rhosol_strat(nbox)
REAL :: gc_eq(nbox)
REAL :: diff_gc_mode(nbox)
REAL :: maxevap(nbox)
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CONDEN'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ageterm1(:,:,:)=0.0
s_cond_s(:)=0.0
sinkarr(:)=0.0
!
! Set value of AA for each mode according to ICONDIAM
! .. these values are taken from Figure 1 Lehtinen et al (2003)
SELECT CASE (icondiam)
CASE (1)
  aa_modes(:) = 0.0     ! as v1_gm4c, take g.m. number radius
CASE (2)
  aa_modes(1) = 2.0 ! continuum regime -- 2nd radial moment
  aa_modes(2) = 1.9
  aa_modes(3) = 1.5
  aa_modes(4) = 1.1 ! molecular regime -- 1st radial moment
  aa_modes(5) = 1.9
  aa_modes(6) = 1.5
  aa_modes(7) = 1.1
CASE DEFAULT
  ierr = 1
  WRITE(umMessage,'(A,I5)')'Unexpected Value of ICONDIAM ',icondiam
  CALL umPrint(umMessage,src='ukca_conden')
  CALL ereport('UKCA_CONDEN',ierr,'Unexpected ICONDIAM value')
END SELECT

DO jv=1,nchemg
  IF (condensable(jv)) THEN

    ! .. Set component into which component will condense
    icp=condensable_choice(jv)

    dmol=dimen(jv)
    mmcg=mm_gas(jv)
    delgc_cond(:,jv)=0.0

    mask1(:) = gc(:,jv) > conc_eps

    sumnc(:)=0.0
    DO imode=1,nmodes
      IF (mode(imode)) THEN

        aa = aa_modes(imode)

        !!          Y2=EXP(2.0*LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE)))
        !! original runs tried just setting equivalent to AA=2.0 (above)
        y2=EXP(0.5*aa*aa*LOG(sigmag(imode))*LOG(sigmag(imode)))
        ! now use values of AA for each mode from Lehtinen et al (2003)

        nc(:,imode)=0.0
        mask2(:) = mask1(:) .AND. (nd(:,imode) > num_eps(imode))

        rp(:)=wetdp(:,imode)*0.5*y2 ! radius to use for condensation

        SELECT CASE (modesol(imode))
        CASE (1)
          se=se_sol
        CASE (0)
          se=se_ins
        END SELECT

        !!          IF(JV == MH2SO4  ) DIFVOL=DH2SO4 ! as H2SO4+H2O hydrate
        !!          IF(JV == msec_org) DIFVOL=DSECOR ! as OH-a-pinene radical
        IF (jv == mh2so4) THEN
          difvol = 51.96          ! values from dist_data (BIN)
        ELSE IF (jv == msec_org) THEN
          difvol = 204.14         ! values from dist_data (BIN)
        ELSE            ! trap if DIFVOL is being used w/o being defined
          ierr = 100+jv
          CALL ereport('UKCA_CONDEN',ierr,'DIFVOL remains undefined')
        END IF

        ! ..  Calculate change in condensable cpt conc (molecules cm^-3)
        CALL ukca_cond_coff_v(nbox,mask2,rp,tsqrt,airdm3,rhoa,        &
               mmcg,se,dmol,ifuchs,cc,sinkarr,pmid,t,difvol,idcmfp)
        WHERE (mask2(:))
          nc(:,imode)=nd(:,imode)*cc(:)
          sumnc(:)=sumnc(:)+nc(:,imode)
        END WHERE
        !!          IF(CONDENSABLE_CHOICE(JV).EQ.1) THEN ! if H2SO4
        IF (jv == mh2so4) THEN ! if H2SO4
          s_cond_s(:)=s_cond_s(:)+nd(:,imode)*sinkarr(:)
        END IF

      END IF ! if mode is present
    END DO ! Over modes
    !
    WHERE (mask1(:))                                                 &
     delgc_cond(:,jv)=gc(:,jv)*(1.0-EXP(-sumnc(:)*dtz))

    ! .. Update condensable cpt concentration (molecules cm^-3)
    mask2(:) = mask1(:) .AND. (delgc_cond(:,jv) > conc_eps)
    WHERE (mask2(:) .AND. delgc_cond(:,jv) > gc(:,jv))
      delgc_cond(:,jv)=delgc_cond(:,jv)/gc(:,jv) ! make sure no -ves
    END WHERE
    WHERE (mask2(:)) gc(:,jv)=gc(:,jv)-delgc_cond(:,jv)

    DO imode=1,4 ! loop over sol modes (do cond sol -> ins here too)
      IF (mode(imode)) THEN

        !         Calculate increase in total & cpt masses in each soluble mode
        deltams(:)=0.0
        deltami(:)=0.0

        mask3 (:) = mask2(:) .AND. (nd(:,imode  ) > num_eps(imode))
        mask3i(:) = mask2(:) .AND. (nd(:,imode+3) > num_eps(imode))

        IF (imode == 1) THEN

          IF ((icp == cp_su) .AND. (nmascondsunucsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsunucsol)=                           &
              bud_aer_mas(:,nmascondsunucsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocnucsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondocnucsol)=                           &
              bud_aer_mas(:,nmascondocnucsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsonucsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsonucsol)=                           &
              bud_aer_mas(:,nmascondsonucsol)+deltams(:)

            END WHERE
          END IF

        END IF ! if IMODE=1

        IF (imode == 2) THEN

          IF ((icp == cp_su) .AND. (nmascondsuaitsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsuaitsol)=                           &
              bud_aer_mas(:,nmascondsuaitsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsuaitins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsuaitins)=                           &
              bud_aer_mas(:,nmascondsuaitins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocaitsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondocaitsol)=                           &
              bud_aer_mas(:,nmascondocaitsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocaitins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondocaitins)=                           &
              bud_aer_mas(:,nmascondocaitins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaitsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsoaitsol)=                           &
              bud_aer_mas(:,nmascondsoaitsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaitins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsoaitins)=                           &
              bud_aer_mas(:,nmascondsoaitins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! if IMODE=2

        IF (imode == 3) THEN

          IF ((icp == cp_su) .AND. (nmascondsuaccsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsuaccsol)=                           &
              bud_aer_mas(:,nmascondsuaccsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsuaccins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsuaccins)=                           &
              bud_aer_mas(:,nmascondsuaccins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocaccsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondocaccsol)=                           &
              bud_aer_mas(:,nmascondocaccsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondocaccins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondocaccins)=                           &
              bud_aer_mas(:,nmascondocaccins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaccsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsoaccsol)=                           &
              bud_aer_mas(:,nmascondsoaccsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsoaccins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsoaccins)=                           &
              bud_aer_mas(:,nmascondsoaccins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! if IMODE=3

        IF (imode == 4) THEN

          IF ((icp == cp_su) .AND. (nmascondsucorsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsucorsol)=                           &
              bud_aer_mas(:,nmascondsucorsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_su) .AND. (nmascondsucorins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsucorins)=                           &
              bud_aer_mas(:,nmascondsucorins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondoccorsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondoccorsol)=                           &
              bud_aer_mas(:,nmascondoccorsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_oc) .AND. (nmascondoccorins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondoccorins)=                           &
              bud_aer_mas(:,nmascondoccorins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsocorsol > 0)) THEN
            WHERE (mask3(:))

              deltams(:)=delgc_cond(:,jv)*nc(:,imode)/sumnc(:)

              bud_aer_mas(:,nmascondsocorsol)=                           &
              bud_aer_mas(:,nmascondsocorsol)+deltams(:)

            END WHERE
          END IF

          IF ((icp == cp_so) .AND. (nmascondsocorins > 0)) THEN
            WHERE (mask3i(:))

              deltami(:)=delgc_cond(:,jv)*nc(:,imode+3)/sumnc(:)

              bud_aer_mas(:,nmascondsocorins)=                           &
              bud_aer_mas(:,nmascondsocorins)+deltami(:)

              ageterm1(:,imode-1,jv)=deltami(:)

            END WHERE
          END IF

        END IF ! if IMODE=4

        !         All mass condensed onto sol. & ins. goes to sol. mode
        deltam(:)=deltams(:)+deltami(:)

        WHERE (mask3(:))
          md(:,imode,icp)=                                             &
           (md(:,imode,icp)*nd(:,imode)+deltams(:))/nd(:,imode)
          mdt(:,imode)=                                                &
           (mdt(:,imode)*nd(:,imode)+deltams(:))/nd(:,imode)
        END WHERE

        !!          WHERE(MASK3I(:))
        !!           MD(:,IMODE,ICP)=                                             &
        !!            (MD(:,IMODE,ICP)*ND(:,IMODE)+DELTAMI(:))/ND(:,IMODE)
        !!           MDT(:,IMODE)=                                                &
        !!            (MDT(:,IMODE)*ND(:,IMODE)+DELTAMI(:))/ND(:,IMODE)
        !!          ENDWHERE
        !!
        !! **** HERE HAVE COMMENTED OUT THE UPDATING OF THE SOLUBLE MODE
        !! **** MD,MDT FOR THE CONDENSATION ONTO INSOLUBLE MODES BECAUSE IN
        !! **** THE CASE WHERE THERE ARE NO PARTICLES IN THE AITKEN SOLUBLE
        !! **** MODE BUT THERE ARE IN THE AITKEN INSOLUBLE MODE, THE UPDATING
        !! **** OF THE SOLUBLE MODE MD,MDT WILL NOT BE POSSIBLE UNTIL THE
        !! **** SOLUBLE MODE ND HAS BEEN UPDATED DUE TO THE TRANSFER OF AGED
        !! **** PARTICLES. THIS IS DONE IN THE UKCA_AGEING ROUTINE -- SO THIS IS
        !! **** WHERE THE UPDATE OF THE SOLUBLE MODE MD,MDT FOR THE CONDENSATION
        !! **** ONTO THE SHOULD BE DONE ALSO.
        !!
      END IF ! if mode is present

    END DO ! IMODE=1,4 (soluble modes)

  END IF ! IF CONDENSABLE(JV)

END DO ! DO JV=1,NCHEMG

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_conden
END MODULE ukca_conden_mod
