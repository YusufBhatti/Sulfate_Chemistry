! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to define details of mode emissions for species provided in netCDF
!   files
!
! Method:
!   Contains list of allowed offline aerosol emissions in netCDF files:
!   aero_ems_species.
!   Contains routine to set details of mode emissions for these species:
!   components, mode fractions, ...: ukca_def_mode_emiss
!
!   To add new aerosol emissions species follow these steps:
!    * add the tracer_name to aero_ems_species
!    * add the tracer_name to CASE statements in ukca_def_mode_emiss and
!      specifiy component as well as emitted fraction, diameter and stdev for
!      each mode.
!
! Part of the UKCA model, a community model supported by
! The Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------------
!
MODULE ukca_emiss_mode_mod

USE ukca_mode_setup,     ONLY: nmodes, mode, cp_bc, cp_oc, cp_su, &
                               fracbcem, fracocem, ncp,           &
                               component
USE ukca_option_mod,     ONLY: i_mode_setup, l_ukca_primsu, l_ukca_primbcoc, &
                               l_bcoc_bf, l_bcoc_bm, l_bcoc_ff
USE ukca_constants,      ONLY: m_s, m_so2

USE ereport_mod,         ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,          ONLY: umPrint, umMessage, PrintStatus, PrStatus_Normal

USE missing_data_mod,    ONLY: rmdi, imdi

USE parkind1,            ONLY: jpim, jprb      ! DrHook
USE yomhook,             ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! List of all allowable offline emissions for testing if a given emission
! in a netCDF file is an aerosol emission field.
CHARACTER(LEN=10), PARAMETER :: aero_ems_species(9) = (/"BC_biomass", &
                                                        "OC_biomass", &
                                                        "BC_fossil ", &
                                                        "OC_fossil ", &
                                                        "BC_biofuel", &
                                                        "OC_biofuel", &
                                                        "SO2_low   ", &
                                                        "SO2_high  ", &
                                                        "SO2_nat   "   /)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_MODE_MOD'

CONTAINS

!---------------------------------------------------------------------------
! Description:
!   Routine to define the following details of mode emissions for a given tracer
!   name:
!    * component
!    * mode diameter and stdev for each mode
!    * fraction of emission allocated to each mode
!    * which modes receive some emission
!
! Method:
!   Component and mode frac, diam, stdev are provided in CASE statements.
!   lmode_emiss array of modes receiving emission is derived from mode_frac
!   array, plus arrays of modes and components active for the current scheme.
!
!   Mode emissions choices for sulphate, BC and OC are controlled by the 
!   integer ISO2EMS:
!  
!   ISO2EMS == 1 : primary SO4 ems --> 15%/85% to 10/70 nm g.m.diam. modes as
!                  in Spracklen (2005) and Binkowski & Shankar (1995).
!
!   ISO2EMS == 2 : primary SO4 ems:
!                  road/off-road/domestic      all -->   30nm gm.diam mode
!                  industrial/power-plant/ship all --> 1000nm gm.diam mode
!                  as for original AEROCOM size recommendations.
!  
!   ISO2EMS >= 3 : primary SO4 ems --> 50%/50% to 150/1500nm  g.m.diam. modes
!                  as for Stier et al (2005) modified AEROCOM sizdis
!                  recommendations.
!
!   Note that to be completely consistent with these references, ISO2EMS=1
!   requires mode_parfrac=3.0, while ISO2EMS>=2 requires mode_parfrac=2.5.
!   mode_parfrac is set in namelist run_ukca.
!  
!   ISO2EMS also controls size assumptions for primary carbonaceous aerosol:
!  
!   ISO2EMS /= 3 : biofuel & biomass BC/OC emissions -->  80nm g.m.diam.
!                  fossil-fuel BC/OC emissions --> 30nm g.m.diam.
!  
!   ISO2EMS == 3 : biofuel & biomass BC/OC emissions --> 150nm g.m.diam.
!                  fossil-fuel BC/OC emissions --> 120nm g.m.diam.
! 
!   We use ISO2EMS=1 when chosen only sulphate and sea-salt in 4 modes
!                     (no BC/OC) -- i_mode_setup=1 -- small sizes
!                     make up for lack of primary BC/OC emissions.
!  
!   For all other options use ISO2EMS=3 as standard as in GLOMAP.
!---------------------------------------------------------------------------
SUBROUTINE ukca_def_mode_emiss(tracer_name, lmode_emiss, icp,        &
                               mode_frac_out, mode_diam_out, mode_stdev_out, &
                               lwarn_mismatch)

IMPLICIT NONE

! Subroutine arguments
CHARACTER (LEN=10), INTENT(IN) :: tracer_name  ! name of emission tracer

! Which modes receive some emission? TRUE for modes which do.
LOGICAL, INTENT(OUT) :: lmode_emiss(nmodes)

! Aerosol component for this species (cp_su, cp_bc, etc)
INTEGER, INTENT(OUT) :: icp           

! Fraction of emission applied to each mode (0->1)
REAL, INTENT(OUT), OPTIONAL :: mode_frac_out(nmodes)

! Geometric mean diameter and standard deviation of emiss applied to each mode
REAL, INTENT(OUT), OPTIONAL :: mode_diam_out(nmodes)   ! units = nm
REAL, INTENT(OUT), OPTIONAL :: mode_stdev_out(nmodes)  ! units = 1

! Print warning if this emission will be ignored due to namelist switches or if
! there is a mismatch between active modes and the modes expected for this
! emission. The purpose of the warnings is to alert the user that some
! emissions files provided will be (partially) ignored.
! This routine is called from multiple places and on every timestep, so this
! logical avoids duplicating the warnings: it is only true for the first call
! from ukca_emiss_init.
LOGICAL, INTENT(IN), OPTIONAL :: lwarn_mismatch

! Local variables

! Local versions of optional vars mode_frac_out, mode_diam_out, mode_stdev_out,
! lwarn_mismatch.
! These are used locally and copied to/from subroutine arguments if PRESENT()
REAL :: mode_frac(nmodes)
REAL :: mode_diam(nmodes)
REAL :: mode_stdev(nmodes)
LOGICAL :: lwarn_mismatch_local

! Indicators for whether emissions contain biofuel, biomass burning, or
! fossil fuel BC/OC, for comparison with namelists switches l_bcoc_bf etc.
LOGICAL :: l_emfile_bcoc_bf
LOGICAL :: l_emfile_bcoc_bm
LOGICAL :: l_emfile_bcoc_ff

INTEGER :: iso2ems        ! determines emission assumptions (see below)
INTEGER :: jmode          ! loop index

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_DEF_MODE_EMISS'

CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER                            :: ierror      ! Error code
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise the arrays to enable check that all are set
mode_frac(:) = -1.0
mode_diam(:) = -1.0
mode_stdev(:) = -1.0
icp = 0

! Copy lwarn_mismatch to local variable, or set false if not present
IF (PRESENT(lwarn_mismatch)) THEN
  lwarn_mismatch_local = lwarn_mismatch
ELSE
  lwarn_mismatch_local = .FALSE.
END IF

! Set ISO2EMS according to i_mode_setup (see subroutine header).
IF (i_mode_setup == 1) iso2ems=1       ! SUSS_4mode
IF (i_mode_setup == 2) iso2ems=3       ! SUSSBCOC_5mode
IF (i_mode_setup == 3) iso2ems=3       ! SUSSBCOC_4mode
IF (i_mode_setup == 4) iso2ems=3       ! SUSSBCOCSO_5mode
IF (i_mode_setup == 5) iso2ems=3       ! SUSSBCOCSO_4mode
IF (i_mode_setup == 6) iso2ems=3       ! DUonly_2mode 
IF (i_mode_setup == 7) iso2ems=0       ! DUonly_3mode
IF (i_mode_setup == 8) iso2ems=3       ! SUSSBCOCDU_7mode
IF (i_mode_setup == 9) iso2ems=3       ! SUSSBCOCDU_4mode
IF (i_mode_setup == 10) iso2ems=3      ! SUSSBCOCNTNH_5mode_8cpt

! Initialise bcoc indicators
l_emfile_bcoc_bf = .FALSE.
l_emfile_bcoc_bm = .FALSE.
l_emfile_bcoc_ff = .FALSE.

SELECT CASE(iso2ems)
CASE (3)

  SELECT CASE(TRIM(tracer_name))
  CASE ('BC_biofuel','OC_biofuel')
    ! Component and mode_frac set below for OC and BC.
    l_emfile_bcoc_bf = .TRUE.
    mode_diam(:) = 150.0
    mode_stdev(:) = 1.59
  CASE ('BC_biomass','OC_biomass')
    ! Component and mode_frac set below for OC and BC.
    l_emfile_bcoc_bm = .TRUE.
    mode_diam(:) = 150.0
    mode_stdev(:) = 1.59
  CASE ('BC_fossil','OC_fossil')
    ! FF emission diam is set to 60nm based on Stier et al. 2005
    ! Acceptable range for this parameter  is 30-120nm 
    ! with lower limit from the uncertainty analysis of Lee et al. 2011
    ! and  upper limit of 120nm based on aircraft
    ! measurements from "fresh pollution" flights during TARFOX, ACE2, and
    ! ADRIEX (Osborne et al. 2005; Osborne et al. 2007). 
    ! Component and mode_frac set below for OC and BC.
    l_emfile_bcoc_ff = .TRUE.
    mode_diam(:) = 60.0
    mode_stdev(:) = 1.59
  CASE ('SO2_low','SO2_high')
    icp = cp_su
    mode_frac(:) = (/ 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0 /)
    mode_diam(:) = (/rmdi,rmdi,150.0,1500.0,rmdi,rmdi,rmdi/)
    mode_stdev(:) = (/rmdi,rmdi,1.59, 2.0, rmdi,rmdi,rmdi/)
  CASE ('SO2_nat')
    icp = cp_su
    mode_frac(:) = (/ 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0 /)
    mode_diam(:) = (/rmdi,60.0,150.0,rmdi,rmdi,rmdi,rmdi/)
    mode_stdev(:) = (/rmdi,1.59,1.59,rmdi,rmdi,rmdi,rmdi/)
  CASE ('PMOC')
    icp = cp_oc
    !these settings may change in the future
    mode_frac(:) = (/ 0.0, 0.25, 0.0, 0.0, 0.75, 0.0, 0.0 /)
    mode_diam(:) = 160.0
    mode_stdev(:) = 2.0
  CASE DEFAULT
    cmessage = 'No emission information coded for ' // TRIM(tracer_name)
    CALL ereport ('UKCA_DEF_MODE_EMISS', iso2ems, cmessage)
END SELECT

CASE DEFAULT
  cmessage = 'This value of iso2ems has not been tested with new emissions'
  CALL ereport ('UKCA_DEF_MODE_EMISS', iso2ems, cmessage)
END SELECT

! Set mode_frac for all OC, BC species using central oc/bcfracem.
SELECT CASE(tracer_name(1:3))
CASE ('OC_')
  icp = cp_oc
  mode_frac(:) = fracocem(:)
CASE ('BC_')
  icp = cp_bc
  mode_frac(:) = fracbcem(:)
END SELECT

! Check that the arrays are now all set
IF (ALL(mode_frac(:) < 0.0) .OR. &
    ALL(mode_diam(:) < 0.0) .OR. &
    ALL(mode_stdev(:) < 0.0)) THEN
  cmessage = 'Mode information has not been set for ' // TRIM(tracer_name)
  CALL ereport ('UKCA_DEF_MODE_EMISS', iso2ems, cmessage)
END IF

! Set lmode_emiss to indicate which modes will receive emissions.
! lmode_emiss is only true for a mode if the mode is active, the component is
! present in the model and the emission into that mode is non-zero.
lmode_emiss(:) = (mode(:) .AND. &
                  component(:,icp) .AND. &
                  (mode_frac(:) > 0.0))

! There are individual logicals which can switch off particular emissions - use
! these to override lmode_emiss and issue a warning so that the user knows the
! file is not being used.
IF ((.NOT. l_ukca_primsu .AND. icp == cp_su) .OR. &
    (.NOT. l_ukca_primbcoc .AND. icp == cp_bc) .OR. &
    (.NOT. l_ukca_primbcoc .AND. icp == cp_oc) .OR. &
    (.NOT. l_bcoc_bf .AND. l_emfile_bcoc_bf) .OR. &
    (.NOT. l_bcoc_bm .AND. l_emfile_bcoc_bm) .OR. &
    (.NOT. l_bcoc_ff .AND. l_emfile_bcoc_ff) ) THEN

  lmode_emiss(:) = .FALSE.

  IF (lwarn_mismatch_local) THEN
    cmessage = "Emission file provided for " // TRIM(tracer_name) // " but" //&
               " emissions for this component/species switched off by namelist"
    ierror = -icp
    CALL ereport('UKCA_DEF_MODE_EMISS', ierror, cmessage)
  END IF
END IF

! Warn if we have defined emissions as entering a mode which is not used for
! this component.
IF (lwarn_mismatch_local) THEN
  DO jmode = 1, nmodes
    IF ((mode_frac(jmode) > 0.0) .AND. .NOT. component(jmode,icp)) THEN
      WRITE(umMessage,'(4A,7F6.3,A,7L2)') "Emissions for ", TRIM(tracer_name),&
         " defined as entering a mode which is not used for this component.", &
         "mode_frac:", mode_frac(:), "component:", component(:,icp)
      CALL umPrint(umMessage,src='ukca_def_mode_emiss')

      cmessage = "Mismatch between emission and active mode/component. " // &
                 "Total mass emitted will not agree with mass in emission file."
      ierror = -icp
      CALL ereport('UKCA_DEF_MODE_EMISS', ierror, cmessage)
      EXIT
    END IF
  END DO
END IF

IF (PRESENT(mode_frac_out)) mode_frac_out = mode_frac
IF (PRESENT(mode_diam_out)) mode_diam_out = mode_diam
IF (PRESENT(mode_stdev_out)) mode_stdev_out = mode_stdev

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_def_mode_emiss


END MODULE ukca_emiss_mode_mod
