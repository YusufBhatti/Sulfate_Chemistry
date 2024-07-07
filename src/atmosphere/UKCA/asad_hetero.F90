! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Heterogeneous chemistry routine (includes aqueous phase reactions).
!
!     The purpose of this routine is to set and return the heterogeneous
!     reaction rates. If the user has heterogeneous chemistry turned on
!     then this subroutine will be called. The user must supply their
!     own version of this routine to compute the heterogeneous rates.
!
!     Note that this subroutine is called repeatedly. It should not
!     therefore be used to do any I/O unless absolutely necessary. The
!     routine inihet is provided to initialise the heterogeneous chemist
!     by reading in files etc.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_hetero_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_HETERO_MOD'

CONTAINS

SUBROUTINE asad_hetero(n_points, cld_f, cld_l, rc_het)

USE asad_findreaction_mod, ONLY: asad_findreaction
USE asad_mod,        ONLY: t, p, tnd, rk, ih_o3, ih_h2o2, ih_so2, &
                           ih_hno3, ihso3_h2o2, iho2_h, in2o5_h,  &
                           iso3_o3, ihso3_o3, iso2_oh, ih2o2_oh,  &
                           ihno3_oh, spb, sph, nbrkx, nhrkx,      &
                           jpspb, jpsph, jpeq, ih_o3_const,       &
                           idms_o3, imsia_o3, imsi_o3,            &!LER
                           ihso3_hobr, iso3_hobr,                 &!LER
                           ih_dms, ih_msia, ih_hobr !LER
USE ukca_option_mod, ONLY: L_ukca_achem, L_ukca_trophet,          &
                           L_ukca_tropisop, L_ukca_mode,          &
                           jpbk, jphk, jpdw, l_ukca_offline,      &
                           l_ukca_offline_be, l_ukca_nr_aqchem
USE ukca_chem_offline, ONLY: nwet_constant
USE ukca_fdiss_constant_mod, ONLY: ukca_fdiss_constant
USE ukca_constants,  ONLY: avc, rhow, m_air, H_plus
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook
USE ereport_mod,     ONLY: ereport
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_fdiss_mod, ONLY: ukca_fdiss
IMPLICIT NONE


INTEGER, INTENT(IN) :: n_points         ! No of spatial points

REAL, INTENT(IN) :: cld_f(n_points)     ! Cloud fraction
REAL, INTENT(IN) :: cld_l(n_points)     ! Cloud liquid water (kg/kg)
REAL, INTENT(IN) :: rc_het(n_points,2)  ! Heterog. Chem. Rates (tropospheric)

!     Local variables

REAL, PARAMETER    :: qcl_min = 1.0e-12 ! do calcs when qcl > qcl_min
REAL               :: vr(n_points)      ! volume ratio
REAL               :: fdiss(n_points, jpdw, jpeq+1)
                                        ! fractional dissociation array
                                        ! final index: 1) dissolved
                                        !              2) 1st dissociation
                                        !              3) 2nd dissociation
REAL, ALLOCATABLE  :: fdiss_constant(:,:,:)
                                    ! As fdiss, but for constant species
REAL               :: fdiss_o3(n_points) ! fractional dissociation for O3

INTEGER            :: icode = 0         ! Error code
CHARACTER (LEN=errormessagelength) :: cmessage          ! Error message
CHARACTER(LEN=10)  :: prods(2)          ! Products
LOGICAL, SAVE      :: first = .TRUE.    ! Identifies first call
LOGICAL, SAVE      :: first_pass = .TRUE.    ! Identifies if thread has
                                             ! been through CRITICAL region
LOGICAL            :: todo(n_points)    ! T where cloud frac above threshold

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_HETERO'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!       1. Identify equations and calculate fractional dissociation
!          --------------------------------------------------------

IF (first_pass) THEN
  ! OMP CRITICAL will only allow one thread through this code at a time,
  ! while the other threads are held until completion.
!$OMP CRITICAL (asad_hetero_init)
  IF (first) THEN
    IF (L_ukca_achem) THEN
      ! Check that the indicies of the aqueous arrays are identified
      IF (ih_o3 == 0 .OR. ih_h2o2 == 0 .OR. ih_so2 == 0 .OR. ih_hno3 == 0 &
          .OR. ih_dms == 0 .OR. ih_msia == 0 .OR. ih_hobr == 0) THEN !LER
        cmessage=' Indicies for Aqueous chemistry uninitialised'//  &
                 ' - O3, H2O2, SO2, and HNO3 etc must be made '//       &
                 ' soluble species in CHCH_DEFS array'
        WRITE(umMessage,'(A9,I5)') 'ih_o3:   ',ih_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_h2o2: ',ih_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_hno3: ',ih_hno3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_so2:  ',ih_so2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_dms:  ',ih_dms !LER Mar2019
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_hobr:  ',ih_hobr !LER Mar2019
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'ih_msia:  ',ih_msia !LER Mar2019
        CALL umPrint(umMessage,src='asad_hetero')
        icode=1
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF

    IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
      ! Check that the indices of the aqueous arrays are identified
      IF (ih_h2o2 == 0 .OR. ih_so2 == 0 ) THEN
        cmessage=' Indices for Aqueous chemistry uninitialised'//  &
                 ' - H2O2, and SO2, must be made '//                &
                 ' soluble species in CHCH_DEFS array'
        icode=1
      END IF

      IF (icode /= 0) THEN
        WRITE(umMessage,'(A10,I5)') 'ih_h2o2: ',ih_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A10,I5)') 'ih_so2:  ',ih_so2
        CALL umPrint(umMessage,src='asad_hetero')
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF

    IF (L_ukca_trophet .AND. .NOT. L_ukca_mode) THEN
      cmessage=' Tropospheric heterogeneous chemistry is flagged'// &
             ' but MODE aerosol scheme is not in use'
      icode=1
      CALL ereport('ASAD_HETERO',icode,cmessage)
    END IF

    ! Find reaction locations
    ihso3_h2o2 = 0
    iso3_o3    = 0
    ihso3_o3   = 0
    ih2o2_oh   = 0
    ihno3_oh   = 0
    in2o5_h    = 0
    iho2_h     = 0
    idms_o3    = 0 !LER Mar2019
    imsia_o3   = 0 !LER Mar2019
    imsi_o3    = 0 !LER Mar2019
    ihso3_hobr = 0 !LER Mar2019
    iso3_hobr  = 0 !LER Mar2019


    IF (l_ukca_nr_aqchem .OR. l_ukca_offline_be) THEN

      prods = (/'NULL0     ','          '/)
      ihso3_h2o2 = asad_findreaction( 'SO2       ', 'H2O2      ',   &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL1     ','          '/)   ! Identifies HSO3- + O3(aq)
      ihso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',     &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL2     ','          '/)   ! Identifies SO3-- + O3(aq)
      iso3_o3 = asad_findreaction( 'SO2       ', 'O3        ',      &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL3     ','          '/)   ! Identifies DMS(aq) + O3(aq); LER Mar2019
      idms_o3 = asad_findreaction( 'DMS       ', 'O3        ',      &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL4     ','          '/)   ! Identifies MSIA(aq) + O3(aq); LER Mar2019
      imsia_o3 = asad_findreaction( 'MSIA      ', 'O3        ',      &
                               prods, 1, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL5     ','          '/)   ! Identifies MSI-(aq) + O3(aq); LER Mar2019
      imsi_o3 = asad_findreaction( 'MSIA      ', 'O3        ',      &
                               prods, 1, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL6     ','          '/)   ! Identifies HSO3- + HOBr(aq); LER Mar2019
      ihso3_hobr = asad_findreaction( 'SO2       ', 'HOBr      ',      &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'NULL7     ','          '/)   ! Identifies SO3-- + HOBr(aq); LER Mar2019
      iso3_hobr = asad_findreaction( 'SO2       ', 'HOBr       ',      &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )

      IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
        prods = (/'H2O       ','          '/)
      ELSE
        prods = (/'H2O       ','HO2       '/)
      END IF
      ih2o2_oh = asad_findreaction( 'H2O2      ', 'OH        ',     &
                               prods, 2, spb, nbrkx, jpbk+1, jpspb )
      prods = (/'H2O       ','NO3       '/)
      ihno3_oh = asad_findreaction( 'HONO2     ', 'OH        ',     &
                               prods, 2, spb, nbrkx, jpbk+1, jpspb )

      icode = 0
      IF (ihso3_h2o2 == 0 .OR. iso3_o3 == 0 .OR. ih2o2_oh == 0 .OR. &
          ihso3_o3 == 0 ) THEN
        WRITE(umMessage,'(A12,I5)') 'ihso3_h2o2: ',ihso3_h2o2
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'ihso3_o3: ',ihso3_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'iso3_o3: ',iso3_o3
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A12,I5)') 'ih2o2_oh: ',ih2o2_oh
        CALL umPrint(umMessage,src='asad_hetero')
        icode = 1
      END IF
      IF (l_ukca_achem .AND. ihno3_oh == 0) THEN
        icode = 1
        WRITE(umMessage,'(A12,I5)') 'ihno3_oh: ',ihno3_oh
        CALL umPrint(umMessage,src='asad_hetero')
      END IF
      IF (icode > 0) THEN
        cmessage=' Heterogenous chemistry called, but eqns'//       &
                  ' not found - see output'
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF
    END IF   ! L_ukca_achem.....

    ! Search for tropospheric heterogeneous reactions
    IF (l_ukca_trophet) THEN
      prods = (/'HONO2     ','          '/)
      in2o5_h = asad_findreaction( 'N2O5      ', '          ',        &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )
      prods = (/'H2O2      ','          '/)
      iho2_h = asad_findreaction( 'HO2       ', '          ',         &
                               prods, 2, sph, nhrkx, jphk+1, jpsph )

      IF (iho2_h == 0 .OR. in2o5_h == 0) THEN
        WRITE(umMessage,'(A9,I5)') 'in2o5_h: ',in2o5_h
        CALL umPrint(umMessage,src='asad_hetero')
        WRITE(umMessage,'(A9,I5)') 'iho2_h: ',iho2_h
        CALL umPrint(umMessage,src='asad_hetero')
        cmessage=' Tropospheric heterogenous chemistry is flagged,'// &
                 ' but equations not found - see output'
        icode = 1
        CALL ereport('ASAD_HETERO',icode,cmessage)
      END IF   ! iho3_h=0 etc

    END IF      ! l_ukca_trophet

    first = .FALSE.

  END IF      ! first
!$OMP END CRITICAL (asad_hetero_init)
END IF        ! first_pass

IF ((l_ukca_nr_aqchem .OR. l_ukca_offline_be) .AND. ANY(cld_l > qcl_min)) THEN
  CALL ukca_fdiss(n_points, qcl_min, t, p, cld_l, fdiss)
  todo(:) = cld_l(:) > qcl_min
ELSE
  fdiss(:,:,:) = 0.0
  todo(:) = .FALSE.
END IF

IF (ANY(cld_l > qcl_min)) THEN
  IF ((l_ukca_offline .OR. l_ukca_offline_be) .AND. nwet_constant > 0 ) THEN
    ALLOCATE(fdiss_constant(n_points, nwet_constant, jpeq+1))
    CALL ukca_fdiss_constant(n_points, qcl_min, t, p, cld_l,        &
                           fdiss_constant)
    fdiss_o3(:) = fdiss_constant(:,ih_o3_const,1)
    DEALLOCATE(fdiss_constant)
  ELSE
    fdiss_o3(:) = fdiss(:,ih_o3,1)
  END IF
END IF

!       2. Calculate heterogeneous rates and reduce rates due to aqueous fraction
!          ----------------------------------------------------------------------

IF (l_ukca_nr_aqchem .OR. l_ukca_offline_be) THEN
  WHERE (todo(:))

    ! Convert clw in kg/kg to volume ratio
    vr(:) = cld_l(:)*tnd(:)*m_air*(1e6/avc)/rhow

!     HSO3- + H2O2(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
    rk(:,ihso3_h2o2) = 2.1295e+14*EXP(-4430.0/t(:))*                &
       (H_plus/(1.0 + 13.0*H_plus))*cld_f(:)*fdiss(:,ih_so2,2)*     &
       fdiss(:,ih_h2o2,1)*1000.0/(avc*vr(:))

!     HSO3- + O3(aq) => SO4--  [Kreidenweis et al. (2003), optimised]
    rk(:,ihso3_o3) = 4.0113e+13*EXP(-5530.0/t(:))*                  &
                     cld_f(:)*fdiss(:,ih_so2,2)*fdiss_o3(:)*        &
                     1000.0/(avc*vr(:))

    ! SO3-- + O3(aq) => SO4-- [Kreidenweis et al. (2003), optimised]
    rk(:,iso3_o3) = 7.43e+16*EXP(-5280.0/t(:))*cld_f(:)*            &
                       fdiss(:,ih_so2,3)*fdiss_o3(:)*               &
                       1000.0/(avc*vr(:))

    ! H2O2 + OH: reduce to take account of dissolved fraction
    rk(:,ih2o2_oh) = rk(:,ih2o2_oh)*                                &
          (1.0 - (fdiss(:,ih_h2o2,1)+fdiss(:,ih_h2o2,2))*cld_f(:))
    
    !Chen et al (2018), LER, UC, Mar2019:
    ! DMS(aq) + O3(aq) => DMSO [Chen et al. (2018)], LER
    rk(:,idms_o3) = 5.3e+12*EXP(-2600.0/t(:))*cld_f(:)*            &
                       fdiss(:,ih_dms,1)*fdiss_o3(:)*               &
                       1000.0/(avc*vr(:))

    !Chen et al (2018), LER, UC, Jun2019:
    ! DMS(aq) + O3(aq) => DMSO [Chen et al. (2018)], LER
    !rk(:,idms_o3) = 5.3e+12*EXP(-2600.0*((1.0/t(:))-(1.0/298.0)))*cld_f(:)*& # YAB changed from 8.61e+08 to 5.3e+12
    !                   fdiss(:,ih_dms,1)*fdiss_o3(:)*               &
    !                   1000.0/(avc*vr(:))

    ! MSIA(aq) + O3(aq) => MSA [Chen et al. (2018)], LER
    rk(:,imsia_o3) = 3.50e+07*cld_f(:)*            &
                       fdiss(:,ih_msia,1)*fdiss_o3(:)*               &
                       1000.0/(avc*vr(:))

    ! MSI- + O3(aq) => MS [Chen et al. (2018)], LER
    rk(:,imsi_o3) = 2.00e+06*cld_f(:)*            &
                       fdiss(:,ih_msia,2)*fdiss_o3(:)*               &
                       1000.0/(avc*vr(:))

    ! HSO3- + H2O2(aq) => SO4  [Chen et al. (2018)], LER
!    rk(:,ihso3_h2o2) = 2.36e+03*EXP(-4760.0/t(:))*                &
!       (H_plus/(1.0 + 13.0*H_plus))*cld_f(:)*fdiss(:,ih_so2,2)*     &
!       fdiss(:,ih_h2o2,1)*1000.0/(avc*vr(:))

    ! HSO3- + H2O2(aq) => SO4  [Chen (2018)], LER, Jun2019
!    rk(:,ihso3_h2o2) = 2.36e+03*EXP(-4760.0* &
!       ((1.0/t(:))-(1.0/298.0)))* &
!       (H_plus/(1.0 + 13.0*H_plus))*cld_f(:)*fdiss(:,ih_so2,2)*     &
!       fdiss(:,ih_h2o2,1)*1000.0/(avc*vr(:))

!Check with k value recommended by Kreidenweis(2003) and Jacob(1986)
    ! HSO3- + H2O2(aq) => SO4  [Chen (2018)], LER, Jun2019
    !rk(:,ihso3_h2o2) = 7.45e+07*EXP(-4760.0* &
    !   ((1.0/t(:))-(1.0/298.0)))* &
    !   (H_plus/(1.0 + 13.0*H_plus))*cld_f(:)*fdiss(:,ih_so2,2)*     &
    !   fdiss(:,ih_h2o2,1)*1000.0/(avc*vr(:))

    ! HSO3- + O3(aq) => SO4  [Chen et al. (2018)], LER
!    rk(:,ihso3_o3) = 3.2e+05*EXP(-4830.0/t(:))*                  &
!                     cld_f(:)*fdiss(:,ih_so2,2)*fdiss_o3(:)*        &
!                     1000.0/(avc*vr(:))

    ! HSO3- + O3(aq) => SO4  [Chen et al. (2018)], LER, UC, Jun2019
!    rk(:,ihso3_o3) = 3.2e+05*EXP(-4830.0*((1.0/t(:))-(1.0/298.0)))* &
!                     cld_f(:)*fdiss(:,ih_so2,2)*fdiss_o3(:)*        &
!                     1000.0/(avc*vr(:))


    ! SO3-- + O3(aq) => SO4 [Chen et al. (2018)], LER
!    rk(:,iso3_o3) = 1.00e+09*EXP(-4030.0/t(:))*cld_f(:)*            &
!                       fdiss(:,ih_so2,3)*fdiss_o3(:)*               &
!                       1000.0/(avc*vr(:))

    ! SO3-- + O3(aq) => SO4 [Chen et al. (2018)], LER, UC, Jun2019
!    rk(:,iso3_o3) = 1.00e+09*EXP(-4030.0*((1.0/t(:))-(1.0/298.0)))*cld_f(:)* &
!                       fdiss(:,ih_so2,3)*fdiss_o3(:)*               &
!                       1000.0/(avc*vr(:))

    ! HSO3- + HOBr(aq) => SO4 + 2H + Br [Chen et al. (2018)], LER
    rk(:,ihso3_hobr) = 3.20e+09*cld_f(:)*            &
                       fdiss(:,ih_so2,2)*fdiss(:,ih_hobr,1)*               &
                       1000.0/(avc*vr(:))

    ! SO3 + HOBr(aq) => SO4 + H + Br [Chen et al. (2018)], LER
    rk(:,iso3_hobr) = 5.00e+09*cld_f(:)*            &
                       fdiss(:,ih_so2,3)*fdiss(:,ih_hobr,1)*               &
                       1000.0/(avc*vr(:))

  ELSEWHERE
    rk(:,idms_o3) = 0.0 ! LER
    rk(:,imsia_o3) = 0.0 ! LER
    rk(:,imsi_o3) = 0.0 ! LER
    rk(:,ihso3_h2o2) = 0.0
    rk(:,ihso3_o3) = 0.0
    rk(:,iso3_o3) = 0.0
    rk(:,ihso3_hobr) = 0.0 ! LER
    rk(:,iso3_hobr) = 0.0 ! LER
  END WHERE

END IF      ! L_ukca_achem .....

IF (l_ukca_achem) THEN
  !  HNO3 + OH : reduce to take account of dissolved fraction
  WHERE (todo(:))
    rk(:,ihno3_oh) = rk(:,ihno3_oh)*                              &
        (1.0 - (fdiss(:,ih_hno3,1)+fdiss(:,ih_hno3,2))*cld_f(:))
  END WHERE
END IF

IF (L_ukca_trophet) THEN
  ! N2O5 => HNO3 (heterogenous)
  rk(:,in2o5_h) = rc_het(:,1)

  ! HO2 + HO2 => H2O2 (heterogenous)
  rk(:,iho2_h) = rc_het(:,2)
ELSE
  IF (in2o5_h > 0) rk(:,in2o5_h) = 0.0
  IF (iho2_h > 0)  rk(:,iho2_h) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_hetero
END MODULE asad_hetero_mod
