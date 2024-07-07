! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Specify surface boundary conditions for chemical tracers
!
!  This routine is designed to read a "Representative Concentration
!  Pathways" (RCP) file which can be downloaded via
!
!    http://cmip-pcmdi.llnl.gov/cmip5/forcing.html
!
!  It should be noted that the data in these files is
!
!    VALID FOR THE MIDDLE OF THE YEAR ONLY
!
!  If you are making up a new file to be read-in by this routine
!  it should follow the same time-convention.
!
!  Files provided by this page are:
!
!    PICNTRL_MIDYR_CONC.DAT
!    PRE2005_MIDYR_CONC.DAT
!    RCP3PD_MIDYR_CONC.DAT
!    RCP45_MIDYR_CONC.DAT
!    RCP6_MIDYR_CONC.DAT
!    RCP6SCP6TO45_MIDYR_CONC.DAT
!    RCP85_MIDYR_CONC.DAT
!
!  NOTE: These files were created on a windows computer, and so
!        they may need to be passed through dos2unix to be able
!        to be read correctly.
!
!  These files contain useful information about the file in a
!  'fortran-type' namelist in the file header. It is important that
!  this namelist is there, but due to the non-fortran format of this
!  namelist in some circumstances it will use a '/' character, which
!  will cause errors in a normal read line statement. This is solved
!  buy using the GET_RCP_NML subroutine contained within this module.
!
!  If perform_lumping is .TRUE. then the concentrations of several Cl
!  species will be added to either CFC11 or CFC12 (which goes where
!  is loosly judged by the lifetime of the species), and the
!  concetrations of several Br species will be added to MeBr (CH3Br).
!  This is done in exactly the same way as in the WMOA1 routine.
!  It should be noted that this may cause issues if radiative feedback
!  from UKCA CFC fields is requested.
!
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
!  Called from UKCA_SCENARIO_CTL
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_scenario_rcp_mod

USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus, &
                      PrStatus_Normal, PrStatus_Oper, newline

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_scenario_rcp
PUBLIC :: test_scenario_rcp_ctl

! logical to turn on testing of RCP code (default F)
LOGICAL, PUBLIC, PARAMETER :: L_ukca_test_scenario_rcp=.FALSE.

! values required from namelist in RCP file. Set to initial values which are
! trapped against in GET_RCP_NML routine
INTEGER, SAVE           :: thisfile_datacolumns  = -1
INTEGER, SAVE           :: thisfile_firstyear    = -1
INTEGER, SAVE           :: thisfile_lastyear     = -1
REAL   , SAVE           :: thisfile_annualsteps  = -1.0
INTEGER, SAVE           :: thisfile_firstdatarow = -1
INTEGER, SAVE           :: thisfile_units        = -1
CHARACTER(LEN=12), SAVE :: thisfile_dattype      = 'unset'


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SCENARIO_RCP_MOD'

CONTAINS

SUBROUTINE ukca_scenario_rcp(n_lbc_specs, lbc_specs,              &
                      i_year, i_day_number,                       &
                      lbc_mmr, perform_lumping)

USE ukca_constants,    ONLY: c_cf2cl2,c_cf2clcfcl2, c_cf2clcf2cl, &
    c_cf2clcf3, c_cfcl3, c_ccl4, c_meccl3, c_chf2cl, c_mecfcl2,   &
    c_mecf2cl, c_cf2clbr, c_mecl, c_mebr, c_cf2br2, c_cf3br,      &
    c_cf2brcf2br, c_n2o, c_ch4, c_co2, c_cf3chf2, c_ch2fcf3
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE ukca_option_mod,   ONLY: ukca_RCPdir, ukca_RCPfile
USE ereport_mod,       ONLY: ereport
USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs           ! Number of LBC species
REAL, INTENT(INOUT)      :: lbc_mmr(n_lbc_specs)  ! Lower BC MMRs
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs(n_lbc_specs)! LBC species
LOGICAL, INTENT(IN)      :: perform_lumping       ! T for lumping of LBCs
INTEGER, INTENT(IN)      :: i_year                ! year
INTEGER, INTENT(IN)      :: i_day_number          ! day

LOGICAL, SAVE :: first = .TRUE.

! cater for species below

! integers holding location of possible species that could be passed-in into
! the lbc_spec array
INTEGER,SAVE :: icfcl3=0, icf2cl2=0, icf2clcfcl2=0, icf2clcf2cl=0
INTEGER,SAVE :: icf2clcf3=0, iccl4=0, imeccl3=0, ichf2cl=0
INTEGER,SAVE :: imecfcl2=0, imecf2cl=0, icf2clbr=0, icf2br2=0
INTEGER,SAVE :: icf3br=0, icf2brcf2br=0, imecl=0, imebr=0
INTEGER,SAVE :: in2o=0, ich4=0, ih2=0, ico2=0, ich2br2=0, in2=0
INTEGER,SAVE :: icf3chf2=0, ich2fcf3=0
INTEGER,SAVE :: ichbr3=0
INTEGER,SAVE :: icos=0

! integers holding location of a species in the species_array array which is
! used to work out which column of data read-in from the RCP file contains
! which species
INTEGER,SAVE :: jcfcl3=0, jcf2cl2=0, jcf2clcfcl2=0, jcf2clcf2cl=0
INTEGER,SAVE :: jcf2clcf3=0, jccl4=0, jmeccl3=0, jchf2cl=0
INTEGER,SAVE :: jmecfcl2=0, jmecf2cl=0, jcf2clbr=0, jcf2br2=0
INTEGER,SAVE :: jcf3br=0, jcf2brcf2br=0, jmecl=0, jmebr=0
INTEGER,SAVE :: jn2o=0, jch4=0, jco2=0
INTEGER,SAVE :: jcf3chf2=0, jch2fcf3=0

! name of RCP file, generated from namelist input
CHARACTER(LEN=256) :: rcpfile

! number of species considered for lumping purposes
! if more species that need to be lumped are added, then
! this value will need to increase, and the replace_Cl/Br and
! convfac_Cl/Br arrays will need to be extended and values given
INTEGER, PARAMETER :: n_full_spec=20

REAL :: full_lbc(n_full_spec)

INTEGER :: file_unit

! List which species replaces which in case of lumping. 1 = F11, 2 = F12
! 16 = CH3Br
INTEGER, SAVE :: replace_Cl(n_full_spec) =                        &
                                  (/0,0,2,2,                      &
                                    2,1,1,2,                      &
                                    1,1,1,0,                      &
                                    0,0,1,0,                      &
                                    0,0,0,0/)

INTEGER, SAVE :: replace_Br(n_full_spec) =                        &
                                  (/0,0,0,0,                      &
                                    0,0,0,0,                      &
                                    0,0,16,16,                    &
                                    16,16,0,0,                    &
                                    0,0,0,0/)

! conversions factors to use when lumping:
! Cl
REAL, SAVE :: convfac_Cl(n_full_spec)
! Br
REAL, SAVE :: convfac_Br(n_full_spec)
! these are filled in the IF (first) ... block

! local variables
      ! for time interpolation
REAL :: ftime
REAL :: frac
INTEGER :: yindex ! year index
! loop variables
INTEGER :: i
INTEGER :: j

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage

! allocated on first TS, read-in from file, and then interpolated in time
CHARACTER(LEN=255) :: temp_string
CHARACTER(LEN=3),  ALLOCATABLE, SAVE :: units_array(:)
CHARACTER(LEN=16), ALLOCATABLE, SAVE :: species_array(:)

! the year value that is read-in, before being converted to REAL (and increased
! by 0.5 as the file defines the mid-year values
INTEGER :: int_year

! number of years
INTEGER, SAVE :: nyear
! year array, needed for interpolation
REAL, ALLOCATABLE, SAVE :: year(:)
! these define the order of the species in the file, and what units they are
! (i.e. ppm, ppb, ppt - all in volume mixing ratio)
REAL, ALLOCATABLE :: unit_conversion(:)
REAL, ALLOCATABLE :: rcp_specs(:)

! species considered - these are saved and not deallocated
REAL, ALLOCATABLE, SAVE :: cfc11(:)
REAL, ALLOCATABLE, SAVE :: cfc12(:)
REAL, ALLOCATABLE, SAVE :: cfc113(:)
REAL, ALLOCATABLE, SAVE :: cfc114(:)
REAL, ALLOCATABLE, SAVE :: cfc115(:)
REAL, ALLOCATABLE, SAVE :: ccl4(:)
REAL, ALLOCATABLE, SAVE :: meccl3(:)
REAL, ALLOCATABLE, SAVE :: hcfc22(:)
REAL, ALLOCATABLE, SAVE :: hcfc141b(:)
REAL, ALLOCATABLE, SAVE :: hcfc142b(:)
REAL, ALLOCATABLE, SAVE :: h1211(:)
REAL, ALLOCATABLE, SAVE :: h1202(:)
REAL, ALLOCATABLE, SAVE :: h1301(:)
REAL, ALLOCATABLE, SAVE :: h2402(:)
REAL, ALLOCATABLE, SAVE :: mecl(:)
REAL, ALLOCATABLE, SAVE :: mebr(:)
REAL, ALLOCATABLE, SAVE :: n2o(:)
REAL, ALLOCATABLE, SAVE :: ch4(:)
REAL, ALLOCATABLE, SAVE :: co2(:)
REAL, ALLOCATABLE, SAVE :: cf3chf2(:)
REAL, ALLOCATABLE, SAVE :: ch2fcf3(:)

! time interpolated values
REAL :: n2o_interp
REAL :: ch4_interp
REAL :: co2_interp

REAL, SAVE :: CH2Br2=0.0  ! invariant contribution
REAL, SAVE :: COS_mmr=0.0 ! invariant contribution
REAL, SAVE :: h2=0.0      ! invariant contribution
REAL, SAVE :: n2=0.0      ! invariant contribution
REAL, SAVE :: CHBr3=0.0   ! invariant contribution

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SCENARIO_RCP'


! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set-up conversion factor arrays and read-in file to get mid-year data
IF (first) THEN
   ! SETUP REQUIRED ON FIRST TIMESTEP ONLY
  first = .FALSE.
  convfac_Cl = (/                                                &
       0.0, 0.0, 1.5*(c_cf2cl2/c_cf2clcfcl2),                    &
       (c_cf2cl2/c_cf2clcf2cl),                                  &
       0.5*(c_cf2cl2/c_cf2clcf3),                                &
       1.3333*(c_cfcl3/c_ccl4), c_cfcl3/c_meccl3,                &
       0.5*(c_cf2cl2/c_chf2cl),                                  &
       0.6667*(c_cfcl3/c_mecfcl2), 0.3333*(c_cfcl3/c_mecf2cl),   &
       0.3333*(c_cfcl3/c_cf2clbr), 0.0,                          &
       0.0,0.0,0.3333*(c_cfcl3/c_mecl),0.0,0.0,0.0,0.0,0.0/)

  convfac_Br = (/                                                &
       0.0, 0.0, 0.0, 0.0,                                       &
       0.0, 0.0, 0.0, 0.0,                                       &
       0.0, 0.0, c_mebr/c_cf2clbr, 2.0*(c_mebr/c_cf2br2),        &
       c_mebr/c_cf3br,2.0*(c_mebr/c_cf2brcf2br), &
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)


  ! create file
  rcpfile = TRIM(ADJUSTL(ukca_RCPdir))//'/'                      &
       //TRIM(ADJUSTL(ukca_RCPfile))

  ! diagnostic output if requested
  IF (PrintStatus >= Prstatus_Oper) THEN
    WRITE(umMessage,'(A,A)')                                    &
               'Taking UKCA Lower BC values from RCP file ',    &
               TRIM(ADJUSTL(rcpfile))
    CALL umPrint(umMessage,src='ukca_scenario_rcp')
  END IF

  ! file will contain a 'fortran' namelist which gives values
  ! needed to correctly allocate arrays etc. Read this in here
  CALL get_rcp_nml(rcpfile)

  ! set up array containing years for RCP data - value is for the
  ! middle of the year
  nyear = INT(REAL(thisfile_lastyear - thisfile_firstyear)       &
       *thisfile_annualsteps) + 1
  ! allocate the year array - is not deallocated
  ALLOCATE(year(1:nyear))

  ALLOCATE(units_array(1:thisfile_datacolumns))
  ALLOCATE(species_array(1:thisfile_datacolumns))

  ! get file unit to open file on
  CALL assign_file_unit(TRIM(ADJUSTL(rcpfile)), &
                        file_unit, handler="fortran")

  ! check that file exists as well when opening
  ! reset ierr to 0
  ierr = 0
  OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(rcpfile)), ACTION='READ', &
       POSITION='rewind',STATUS='old',IOSTAT=ierr, IOMSG=iomessage)
  IF (ierr /= 0) THEN
    cmessage =                                           newline//&
      'Error opening file:'//                            newline//&
      TRIM(ADJUSTL(rcpfile))//                           newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL ereport('UKCA_SCENARIO_RCP',ierr,cmessage)
  END IF

  ! read header and throw this away
  DO i=1,thisfile_firstdatarow-3
    READ(file_unit,FMT='(A)') temp_string
  END DO

  ! read in units and species
  READ(file_unit,*) temp_string,units_array(:)
  READ(file_unit,*) temp_string,species_array(:)

  ALLOCATE(unit_conversion(1:thisfile_datacolumns))
  unit_conversion(:) = 0.0

  ! work out where species are in array, and what the unit
  ! conversion is
  ! rhs column is name in file, if different from actual species
  ! formula
  DO i=1,thisfile_datacolumns
    SELECT CASE (species_array(i))
    CASE ('CFCl3           ','CFC_11          ')
      jcfcl3=i
    CASE ('CF2Cl2          ','CFC_12          ')
      jcf2cl2=i
    CASE ('CF2ClCFCl2      ','CFC_113         ')
      jcf2clcfcl2=i
    CASE ('CF2ClCF2Cl      ','CFC_114         ')
      jcf2clcf2cl=i
    CASE ('CF2ClCF3        ','CFC_115         ')
      jcf2clcf3=i
    CASE ('CCl4            ','CARB_TET        ')
      jccl4=i
    CASE ('MeCCl3          ','MCF             ')
      jmeccl3=i
    CASE ('CHF2Cl          ','HCFC_22         ')
      jchf2cl=i
    CASE ('MeCF2Cl         ','HCFC_142B       ')
      jmecf2cl=i
    CASE ('MeCFCl2         ','HCFC_141B       ')
      jmecfcl2=i
    CASE ('CF2ClBr         ','HALON1211       ')
      jcf2clbr=i
    CASE ('CF2Br2          ','HALON1202       ')
      jcf2br2=i
    CASE ('CF3Br           ','HALON1301       ')
      jcf3br=i
    CASE ('CF2BrCF2Br      ','HALON2402       ')
      jcf2brcf2br=i
    CASE ('MeBr            ','CH3BR           ')
      jmebr=i
    CASE ('MeCl            ','CH3CL           ')
      jmecl=i
    CASE ('CF3CHF2         ','HFC125          ')
      jcf3chf2=i
    CASE ('CH2FCF3         ','HFC134a         ')
      jch2fcf3=i
    CASE ('N2O             ')
      jn2o=i
    CASE ('CH4             ')
      jch4=i
    CASE ('CO2             ')
      jco2=i
    END SELECT
    ! units
    SELECT CASE (units_array(i))
    CASE ('ppm')
      unit_conversion(i) = REAL(1.0e-6)
    CASE ('ppb')
      unit_conversion(i) = REAL(1.0e-9)
    CASE ('ppt')
      unit_conversion(i) = REAL(1.0e-12)
    CASE DEFAULT
      WRITE(umMessage,'(a,a,a,a)') ' INCORRECT UNITS: ',       &
                               units_array(i), ' FOR ',        &
                               species_array(i)
      CALL umPrint(umMessage,src='ukca_scenario_rcp')
      cmessage = 'Incorrect Units in RCP scenario'
      ierr     = i
      CALL ereport('UKCA_SCENARIO_RCP',ierr,cmessage)
    END SELECT
  END DO

  ! species_array and units_array no longer needed
  IF (ALLOCATED(species_array)) DEALLOCATE(species_array)
  IF (ALLOCATED(units_array)) DEALLOCATE(units_array)

  ! now get data from file
  ! allocate arrays
  ALLOCATE(rcp_specs(1:thisfile_datacolumns))
  ALLOCATE(cfc11(1:nyear))
  ALLOCATE(cfc12(1:nyear))
  ALLOCATE(cfc113(1:nyear))
  ALLOCATE(cfc114(1:nyear))
  ALLOCATE(cfc115(1:nyear))
  ALLOCATE(ccl4(1:nyear))
  ALLOCATE(meccl3(1:nyear))
  ALLOCATE(hcfc22(1:nyear))
  ALLOCATE(hcfc141b(1:nyear))
  ALLOCATE(hcfc142b(1:nyear))
  ALLOCATE(h1211(1:nyear))
  ALLOCATE(h1202(1:nyear))
  ALLOCATE(h1301(1:nyear))
  ALLOCATE(h2402(1:nyear))
  ALLOCATE(mecl(1:nyear))
  ALLOCATE(mebr(1:nyear))
  ALLOCATE(n2o(1:nyear))
  ALLOCATE(ch4(1:nyear))
  ALLOCATE(co2(1:nyear))
  ALLOCATE(cf3chf2(1:nyear))
  ALLOCATE(ch2fcf3(1:nyear))
  ! initialise to zero
  cfc11 = 0.0
  cfc12 = 0.0
  cfc113 = 0.0
  cfc114 = 0.0
  cfc115 = 0.0
  ccl4 = 0.0
  meccl3 = 0.0
  hcfc22 = 0.0
  hcfc141b = 0.0
  hcfc142b = 0.0
  h1211 = 0.0
  h1202 = 0.0
  h1301 = 0.0
  h2402 = 0.0
  mecl = 0.0
  mebr = 0.0
  n2o = 0.0
  ch4 = 0.0
  co2 = 0.0
  cf3chf2 = 0.0
  ch2fcf3 = 0.0

  ! read in RCP specifications of trace gas concentrations
  DO i=1,nyear
    rcp_specs(:) = 0.0
    READ(file_unit,*) int_year,rcp_specs(:)
    ! Year - VALUES ARE DEFINED FOR THE *MIDDLE* OF THE YEAR
    year(i)  = REAL(int_year) + 0.5
    ! F11
    cfc11(i) = rcp_specs(jcfcl3)*unit_conversion(jcfcl3)*C_CFCl3
    ! F12
    cfc12(i) = rcp_specs(jcf2cl2)*unit_conversion(jcf2cl2)      &
               *C_CF2Cl2
    ! F113
    cfc113(i) = rcp_specs(jcf2clcfcl2)                          &
                *unit_conversion(jcf2clcfcl2)*C_CF2ClCFCl2
    ! F114
    cfc114(i) = rcp_specs(jcf2clcf2cl)                          &
                *unit_conversion(jcf2clcf2cl)*c_cf2clcf2cl
    ! F115
    cfc115(i) = rcp_specs(jcf2clcf3)*unit_conversion(jcf2clcf3) &
                *c_cf2clcf3
    ! CCl4
    ccl4(i) = rcp_specs(jccl4)*unit_conversion(jccl4)*c_ccl4
    ! MeCCl3
    meccl3(i) = rcp_specs(jmeccl3)*unit_conversion(jmeccl3)     &
                *c_meccl3
    ! HCFC-22
    hcfc22(i) = rcp_specs(jchf2cl)*unit_conversion(jchf2cl)     &
                *c_chf2cl
    ! HCFC-141b
    hcfc141b(i) = rcp_specs(jmecfcl2)*unit_conversion(jmecfcl2) &
                  *c_mecfcl2
    ! HCFC-142b
    hcfc142b(i) = rcp_specs(jmecf2cl)*unit_conversion(jmecf2cl) &
                  *c_mecf2cl
    ! H-1211
    h1211(i) = rcp_specs(jcf2clbr)*unit_conversion(jcf2clbr)    &
               *c_cf2clbr
    ! H-1202
    h1202(i) = rcp_specs(jcf2br2)*unit_conversion(jcf2br2)      &
               *c_cf2br2
    ! H-1301
    h1301(i) = rcp_specs(jcf3br)*unit_conversion(jcf3br)*c_cf3br
    ! H-2402
    h2402(i) = rcp_specs(jcf2brcf2br)                           &
               *unit_conversion(jcf2brcf2br)*c_cf2brcf2br
    ! MeCl
    mecl(i) = rcp_specs(jmecl)*unit_conversion(jmecl)*c_mecl
    ! MeBr
    mebr(i) = rcp_specs(jmebr)*unit_conversion(jmebr)*c_mebr
    ! N2O
    n2o(i) = rcp_specs(jn2o)*unit_conversion(jn2o)*c_n2o
    ! CH4
    ch4(i) = rcp_specs(jch4)*unit_conversion(jch4)*c_ch4
    ! CO2
    co2(i) = rcp_specs(jco2)*unit_conversion(jco2)*c_co2
    ! HFC125
    cf3chf2(i) = rcp_specs(jcf3chf2)*unit_conversion(jcf3chf2)  &
                 *c_cf3chf2
    ! HFC134a
    ch2fcf3(i) = rcp_specs(jch2fcf3)*unit_conversion(jch2fcf3)  &
                 *c_ch2fcf3
  END DO

  ! close file
  CLOSE(file_unit)
  CALL release_file_unit(file_unit, handler="fortran")

  ! these arrays are no longer needed, so can be deallocated
  DEALLOCATE(rcp_specs)
  DEALLOCATE(unit_conversion)

  ! Bromoform - constant
  CHBr3  = 6.97904e-12

  ! CH2Br2 (constant)
  !      CH2Br2 = real(3.0e-12)*C_CH2Br2    ! 3 pptv of CH2Br2
  ! as SCENARIO_WMOA1:
  CH2Br2 = 18.0186e-12    ! 3 pptv of CH2Br2

  ! COS (constant)
  COS_mmr = 1.0e-9 ! lbc for COS set as 1.0 ppb MMR

  ! H2 (constant)
  h2 = 3.4528e-8
  ! N2 (constant)
  n2 = 0.754682

END IF ! first


! work out values to interpolate to in time
! currently the routine will give values according to the model
! year and day - however, if the model year is outside the range
! in the file then the output from the routine will be as the
! first/last point would be, i.e.
! if the file goes from 1750-2500 then years <1750 would have the
! 1750 value, and years >2500 would have the 2500 value
ierr = 0
IF (nyear  >   1) THEN
  ftime = REAL(i_year) + REAL(i_day_number - 1)/360.0
  IF (ftime < year(1)) THEN
    yindex = 1
    frac = 0.0
    IF (PrintStatus >= PrStatus_Diag) THEN
       ! print warning message to this effect
      ierr = -1*INT(year(1))
      WRITE(cmessage,FMT='(A,F12.6)')                          &
           'Warning, model year is less than ', year(1)
      CALL ereport('UKCA_SCENARIO_RCP',ierr,cmessage)
    END IF
  ELSE IF (ftime >= year(nyear)) THEN
    yindex = nyear - 1
    frac = 1.0
    IF (PrintStatus >= PrStatus_Diag) THEN
       ! print warning message to this effect
      ierr = -1*INT(year(nyear))
      WRITE(cmessage,FMT='(A,F12.6)')                          &
           'Warning, model year is greater than ', year(nyear)
      CALL ereport('UKCA_SCENARIO_RCP',ierr,cmessage)
    END IF
  ELSE
    yindex = 1
    DO WHILE (ftime > year(yindex + 1))
      yindex = yindex + 1
    END DO
    frac = (ftime-year(yindex))/(year(yindex+1)-year(yindex))
  END IF

  ! INTERPOLATE IN TIME
  ! F11
  full_lbc(1) = (1.0 - frac)* cfc11(yindex)+frac* cfc11(yindex+1)
  ! F12
  full_lbc(2) = (1.0 - frac)* cfc12(yindex)+frac* cfc12(yindex+1)
  ! F113
  full_lbc(3) = (1.0 - frac)*cfc113(yindex)+frac*cfc113(yindex+1)
  ! F114
  full_lbc(4) = (1.0 - frac)*cfc114(yindex)+frac*cfc114(yindex+1)
  ! F115
  full_lbc(5) = (1.0 - frac)*cfc115(yindex)+frac*cfc115(yindex+1)
  ! CCl4
  full_lbc(6) = (1.0 - frac)*  ccl4(yindex)+frac*  ccl4(yindex+1)
  ! MeCCl3
  full_lbc(7) = (1.0 - frac)*meccl3(yindex)+frac*meccl3(yindex+1)
  ! HCFC-22
  full_lbc(8) = (1.0 - frac)*hcfc22(yindex)+frac*hcfc22(yindex+1)
  ! HCFC-141b
  full_lbc(9) = (1.0 - frac)*hcfc141b(yindex)+frac               &
                *hcfc141b(yindex+1)
  ! HCFC-142b
  full_lbc(10)= (1.0 - frac)*hcfc142b(yindex)+frac               &
                *hcfc142b(yindex+1)
  ! H-1211
  full_lbc(11)= (1.0 - frac)* h1211(yindex)+frac* h1211(yindex+1)
  ! H-1202
  full_lbc(12)= (1.0 - frac)* h1202(yindex)+frac* h1202(yindex+1)
  ! H-1301
  full_lbc(13)= (1.0 - frac)* h1301(yindex)+frac* h1301(yindex+1)
  ! H-2402
  full_lbc(14)= (1.0 - frac)* h2402(yindex)+frac* h2402(yindex+1)
  ! MeCl
  full_lbc(15)= (1.0 - frac)*  mecl(yindex)+frac*  mecl(yindex+1)
  ! MeBr
  full_lbc(16)= (1.0 - frac)*  mebr(yindex)+frac*  mebr(yindex+1)

  ! Add 5 pptv of bromine to MeBr to account for the effect of very
  ! short-lived bromine gases
  IF (ich2br2 == 0) THEN
    full_lbc(16) = full_lbc(16) + 5.0e-12*c_mebr
  END IF
  ! CH2Br2
  full_lbc(17)= ch2br2
  ! CHBr3
  full_lbc(18)= chbr3
  ! HFC125
  full_lbc(19)= (1.0 - frac)*cf3chf2(yindex)+frac*cf3chf2(yindex+1)
  ! HFC134a
  full_lbc(20)= (1.0 - frac)*ch2fcf3(yindex)+frac*ch2fcf3(yindex+1)
  ! N2O
  n2o_interp = (1.0 - frac)* n2o(yindex) +frac* n2o(yindex+1)
  ! CH4
  ch4_interp = (1.0 - frac)* ch4(yindex) +frac* ch4(yindex+1)
  ! CO2
  co2_interp = (1.0 - frac)* co2(yindex) +frac* co2(yindex+1)

ELSE ! only one year supplied
  ! F11
  full_lbc(1) = cfc11(1)
  ! F12
  full_lbc(2) = cfc12(1)
  ! F113
  full_lbc(3) = cfc113(1)
  ! F114
  full_lbc(4) = cfc114(1)
  ! F115
  full_lbc(5) = cfc115(1)
  ! CCl4
  full_lbc(6) = ccl4(1)
  ! MeCCl3
  full_lbc(7) = meccl3(1)
  ! HCFC-22
  full_lbc(8) = hcfc22(1)
  ! HCFC-141b
  full_lbc(9) = hcfc141b(1)
  ! HCFC-142b
  full_lbc(10)= hcfc142b(1)
  ! H-1211
  full_lbc(11)= h1211(1)
  ! H-1202
  full_lbc(12)= h1202(1)
  ! H-1301
  full_lbc(13)= h1301(1)
  ! H-2402
  full_lbc(14)= h2402(1)
  ! MeCl
  full_lbc(15)= mecl(1)
  ! MeBr
  full_lbc(16)= mebr(1)

  ! Add 5 pptv of bromine to MeBr to account for the effect of very
  ! short-lived bromine gases
  IF (ich2br2 == 0) THEN
    full_lbc(16) = full_lbc(16) + 5.0e-12*c_mebr
  END IF
  ! CH2Br2
  full_lbc(17)= ch2br2
  ! CHBr3
  full_lbc(18)= chbr3
  ! HFC125
  full_lbc(19)= cf3chf2(1)
  ! HFC134a
  full_lbc(20)= ch2fcf3(1)
  ! N2O
  n2o_interp = n2o(1)
  ! CH4
  ch4_interp = ch4(1)
  ! CO2
  co2_interp = co2(1)
END IF ! nyear>1


! Perform lumping
IF (perform_lumping) THEN
  DO i=1,n_full_spec
    IF (replace_Cl(i) > 0) THEN
      full_lbc(replace_cl(i)) = full_lbc(replace_cl(i)) +         &
                                convfac_cl(i) * full_lbc(i)
    END IF
    IF (replace_Br(i) > 0) THEN
      full_lbc(replace_Br(i)) = full_lbc(replace_Br(i)) +         &
                                convfac_Br(i) * full_lbc(i)
    END IF
  END DO
END IF

icfcl3 = 0
icf2cl2 = 0
icf2clcfcl2 = 0
icf2clcf2cl = 0
icf2clcf3 = 0
iccl4 = 0
imeccl3 = 0
ichf2cl = 0
imecf2cl = 0
imecfcl2 = 0
icf2clbr = 0
icf2br2 = 0
icf3br = 0
icf2brcf2br = 0
imebr = 0
imecl = 0
in2o = 0
ico2 = 0
ich4 = 0
ih2 = 0
ich2br2 = 0
ichbr3 = 0
DO i=1,n_lbc_specs
  SELECT CASE (lbc_specs(i))
  CASE ('CFCl3     ')
    icfcl3=i
  CASE ('CF2Cl2    ')
    icf2cl2=i
  CASE ('CF2ClCFCl2')
    icf2clcfcl2=i
  CASE ('CF2ClCF2Cl')
    icf2clcf2cl=i
  CASE ('CF2ClCF3  ')
    icf2clcf3=i
  CASE ('CCl4      ')
    iccl4=i
  CASE ('MeCCl3    ')
    imeccl3=i
  CASE ('CHF2Cl    ')
    ichf2cl=i
  CASE ('MeCF2Cl   ')
    imecf2cl=i
  CASE ('MeCFCl2   ')
    imecfcl2=i
  CASE ('CF2ClBr   ')
    icf2clbr=i
  CASE ('CF2Br2    ')
    icf2br2=i
  CASE ('CF3Br     ')
    icf3br=i
  CASE ('CF2BrCF2Br')
    icf2brcf2br=i
  CASE ('MeBr      ')
    imebr=i
  CASE ('MeCl      ')
    imecl=i
  CASE ('N2O       ')
    in2o=i
  CASE ('CH4       ')
    ich4=i
  CASE ('H2        ')
    ih2=i
  CASE ('N2        ')
    in2=i
  CASE ('CO2       ')
    ico2=i
  CASE ('CH2Br2    ')
    ich2br2=i
  CASE ('CHBr3     ')
    ichbr3=i
  CASE ('CF3CHF2   ')
    icf3chf2=i
  CASE ('CH2FCF3   ')
    ich2fcf3=i
  CASE ('COS       ')
    icos=i
  END SELECT
END DO


! Assign values to LBC array
! Default values to 0
lbc_mmr = 0.0
IF (icfcl3      > 0) lbc_mmr( icfcl3    ) = full_lbc( 1) ! F11
IF (icf2cl2     > 0) lbc_mmr(icf2cl2    ) = full_lbc( 2) ! F12
IF (icf2clcfcl2 > 0) lbc_mmr(icf2clcfcl2) = full_lbc( 3) ! F113
IF (icf2clcf2cl > 0) lbc_mmr(icf2clcf2cl) = full_lbc( 4) ! F114
IF (icf2clcf3   > 0) lbc_mmr(icf2clcf3  ) = full_lbc( 5) ! F115
IF (iccl4       > 0) lbc_mmr(iccl4      ) = full_lbc( 6) ! CCl4
IF (imeccl3     > 0) lbc_mmr(imeccl3    ) = full_lbc( 7) ! MeCCl3
IF (ichf2cl     > 0) lbc_mmr(ichf2cl    ) = full_lbc( 8) ! HFCF-22
IF (imecfcl2    > 0) lbc_mmr(imecfcl2   ) = full_lbc( 9) ! HCFC-141b
IF (imecf2cl    > 0) lbc_mmr(imecf2cl   ) = full_lbc(10) ! HCFC-142b
IF (icf2clbr    > 0) lbc_mmr(icf2clbr   ) = full_lbc(11) ! H-1211
IF (icf2br2     > 0) lbc_mmr(icf2br2    ) = full_lbc(12) ! H-1202
IF (icf3br      > 0) lbc_mmr(icf3br     ) = full_lbc(13) ! H-1301
IF (icf2brcf2br > 0) lbc_mmr(icf2brcf2br) = full_lbc(14) ! H-2402
IF (imecl       > 0) lbc_mmr(imecl      ) = full_lbc(15) ! MeCl
IF (imebr       > 0) lbc_mmr(imebr      ) = full_lbc(16) ! MeBr
IF (ich2br2     > 0) lbc_mmr(ich2br2    ) = full_lbc(17) ! CH2Br2
IF (ichbr3      > 0) lbc_mmr(ichbr3     ) = full_lbc(18) ! CHBr3
IF (icf3chf2    > 0) lbc_mmr(icf3chf2   ) = full_lbc(19) ! HFC125
IF (ich2fcf3    > 0) lbc_mmr(ich2fcf3   ) = full_lbc(20) ! HFC134a

IF (in2o        > 0) lbc_mmr(in2o       ) = n2o_interp   ! N2O
IF (ich4        > 0) lbc_mmr(ich4       ) = ch4_interp   ! CH4
IF (ico2        > 0) lbc_mmr(ico2       ) = co2_interp   ! CO2

! COS, H2 and N2 not time varying
IF (icos        > 0) lbc_mmr(icos       ) = COS_mmr      ! COS
IF (ih2         > 0) lbc_mmr(ih2        ) = h2           ! H2
IF (in2         > 0) lbc_mmr(in2        ) = n2           ! N2


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_scenario_rcp

SUBROUTINE get_rcp_nml(rcpfile)
USE ereport_mod,      ONLY: ereport
USE file_manager,     ONLY: assign_file_unit, release_file_unit
USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=256), INTENT(IN) :: rcpfile

! length of rcpfile
INTEGER :: rcplen
! temporary strings required when reading in the file and namelist block
CHARACTER(LEN=255) :: c_line
CHARACTER(LEN=255) :: c_tmp01
CHARACTER(LEN=255) :: c_tmp02
CHARACTER(LEN=255) :: c_tmp03
CHARACTER(LEN=255) :: temp_string
! is the annualstep amount defined as a fraction?
LOGICAL :: L_annualstep_fraction
INTEGER :: frac01
INTEGER :: frac02
INTEGER :: asteps

! unit used to read-in file
INTEGER :: file_unit
INTEGER :: i
INTEGER :: j

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_RCP_NML'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! open the file to read the namelist
CALL assign_file_unit(TRIM(ADJUSTL(rcpfile)), &
                      file_unit, handler="fortran")
      
! check that file exists as well when opening
! reset ierr to 0
ierr = 0
OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(rcpfile)),ACTION='READ',&
     POSITION='rewind',STATUS='old',IOSTAT=ierr, IOMSG=iomessage)
IF (ierr /= 0) THEN
  cmessage='Error opening file '//TRIM(ADJUSTL(rcpfile)) // ' : ' // &
  TRIM(iomessage)
  CALL ereport('GET_RCP_NML',ierr,cmessage)
END IF

! work out how long the file is
rcplen=0
ierr=0
len_loop: DO
  rcplen=rcplen+1
  READ(file_unit,FMT='(A)',IOSTAT=ierr) temp_string
  IF (ierr /= 0) EXIT len_loop
END DO len_loop

! have got length of file, so now rewind
REWIND file_unit

! now go through and read-in the values from the 'namelist'
! (but not using a namelist read)
i=0
read_loop: DO
  i=i+1
  IF (i > rcplen) THEN
    cmessage='End of '//TRIM(ADJUSTL(rcpfile))
    ierr    = i
    CALL ereport('GET_RCP_NML',ierr,cmessage)
  END IF
  READ(file_unit,'(A)') c_line
  IF (TRIM(ADJUSTL(c_line)) == &
       '&THISFILE_SPECIFICATIONS') THEN
    DO
      temp_string = ' '
      READ(file_unit,'(A)') c_line
      IF (LEN(TRIM(ADJUSTL(c_line))) > 0) THEN
        READ(c_line,*) temp_string
      ELSE
        READ(c_line,*)
      END IF
      SELECT CASE (TRIM(ADJUSTL(temp_string)))
      CASE ('THISFILE_DATACOLUMNS')
        READ(c_line,FMT='(A27,I11,A1)')                       &
                    c_tmp01,thisfile_datacolumns,c_tmp02
      CASE ('THISFILE_FIRSTYEAR')
        READ(c_line,FMT='(A27,I11,A1)')                       &
                    c_tmp01,thisfile_firstyear,c_tmp02
      CASE ('THISFILE_LASTYEAR')
        READ(c_line,FMT='(A27,I11,A1)')                       &
                    c_tmp01,thisfile_lastyear,c_tmp02
      CASE ('THISFILE_ANNUALSTEPS')
         ! more complicated here, since namelist isn't actually
         ! properly formatted if the number of annual steps is
         ! <1 (i.e. gap between lines is >1 year
         ! check is '/' character is present, and if it is, calculate
         ! the value of the fraction
        L_annualstep_fraction = .FALSE.
        READ(c_line,FMT='(A35,A11,A1)') c_tmp01,c_tmp02,c_tmp03
        frac_loop: DO j=1,LEN(c_tmp02)
          IF (c_tmp02(j:j) == '/') THEN
            L_annualstep_fraction = .TRUE.
            EXIT frac_loop
          END IF
        END DO frac_loop
        IF (L_annualstep_fraction) THEN
          READ(c_line,FMT='(A27,I8,A1,I2,A1)')               &
                      c_tmp01,frac01,c_tmp02,frac02,c_tmp03
          thisfile_annualsteps = REAL(frac01)/REAL(frac02)
        ELSE
           ! '/' character not present, so can read normally
          READ(c_line,FMT='(A27,I11,A1)') c_tmp01,asteps,c_tmp02
          thisfile_annualsteps = REAL(asteps)
        END IF
      CASE ('THISFILE_FIRSTDATAROW')
        READ(c_line,FMT='(A27,I11,A1)')                       &
                    c_tmp01,thisfile_firstdatarow,c_tmp02
      CASE ('THISFILE_UNITS')
        READ(c_line,FMT='(A35,I3,A1)')                        &
                    c_tmp01,thisfile_units,c_tmp02
      CASE ('THISFILE_DATTYPE')
        READ(c_line,FMT='(A27,A)')                            &
                    c_tmp01,thisfile_dattype
        thisfile_dattype = TRIM(ADJUSTL(thisfile_dattype))
      CASE DEFAULT
        EXIT read_loop
      END SELECT
    END DO
  END IF
END DO read_loop

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A,I4)') 'UKCA_RCP:  THISFILE_DATACOLUMNS  = ',&
             thisfile_datacolumns
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,I4)') 'UKCA_RCP:  THISFILE_FIRSTYEAR    = ',&
             thisfile_firstyear
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,I4)') 'UKCA_RCP:  THISFILE_LASTYEAR     = ',&
             thisfile_lastyear
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,F12.6)') 'UKCA_RCP:  THISFILE_ANNUALSTEPS  = ',&
             thisfile_annualsteps
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,I4)') 'UKCA_RCP:  THISFILE_FIRSTDATAROW = ',&
             thisfile_firstdatarow
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,I4)') 'UKCA_RCP:  THISFILE_UNITS        = ',&
             thisfile_units
  CALL umPrint(umMessage,src='get_rcp_nml')
  WRITE(umMessage,'(A,A)') 'UKCA_RCP:  THISFILE_DATTYPE      = ', &
             thisfile_dattype
  CALL umPrint(umMessage,src='get_rcp_nml')
END IF

! we have all the information we need from the namelist now, so close
! the file. It will be opened later to read-in the species concetrations
CLOSE(file_unit)
CALL release_file_unit(file_unit, handler="fortran")

! check for errors
ierr = 0
temp_string = ' '
IF (thisfile_datacolumns < 0) THEN
  ierr = ierr + 1
  temp_string = TRIM(ADJUSTL(temp_string))//' DATACOLUMNS'
END IF
IF (thisfile_firstyear < 0) THEN
  ierr = ierr + 10
  temp_string = TRIM(ADJUSTL(temp_string))//' FIRSTYEAR'
END IF
IF (thisfile_lastyear < 0) THEN
  ierr = ierr + 100
  temp_string = TRIM(ADJUSTL(temp_string))//' LASTYEAR'
END IF
IF (thisfile_annualsteps < 0.0) THEN
  ierr = ierr + 1000
  temp_string = TRIM(ADJUSTL(temp_string))//' ANNUALSTEPS'
END IF
IF (thisfile_firstdatarow < 0) THEN
  ierr = ierr + 10000
  temp_string = TRIM(ADJUSTL(temp_string))//' FIRSTDATAROW'
END IF
IF (thisfile_units < 0) THEN
  ierr = ierr + 100000
  temp_string = TRIM(ADJUSTL(temp_string))//' UNITS'
END IF
IF (TRIM(ADJUSTL(thisfile_dattype)) == 'unset') THEN
  ierr = ierr + 1000000
  temp_string = TRIM(ADJUSTL(temp_string))//' DATTYPE'
END IF
IF (ierr > 0) THEN
  cmessage= 'THISFILE '//TRIM(ADJUSTL(temp_string))//' not filled'
  CALL ereport('GET_RCP_NML',ierr,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE get_rcp_nml

SUBROUTINE test_scenario_rcp_ctl(n_lbc_specs, lbc_specs,            &
                                 perform_lumping)
USE ereport_mod,      ONLY: ereport
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs           ! Number of LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs(n_lbc_specs)! LBC species
LOGICAL, INTENT(IN)      :: perform_lumping       ! T for lumping of LBCs

! loop counters
INTEGER :: i
! to simulate many years to calculate the lbc values
INTEGER :: jday
INTEGER :: jyear

LOGICAL, PARAMETER ::  lumping_off = .FALSE.

REAL :: lbc_c_species(n_lbc_specs)

CHARACTER(LEN=256), PARAMETER :: ukca_test_file=                  &
                                 'Test_RCP_UKCA.dat'
CHARACTER(LEN=256), PARAMETER :: full_test_file=                  &
                                 'Test_RCP_Full.dat'

CHARACTER(LEN=256) :: file_header

INTEGER, PARAMETER :: n_full_specs=26
CHARACTER(LEN=10), ALLOCATABLE :: full_specs(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEST_SCENARIO_RCP_CTL'

! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! this routine calls the routine which tests the
! the scenario_rcp routine by
! cycling through all 1700-2550 and writing the
! lbc_mmr values out to a file so that they can be
! plotted and compared to the ones in the original
! RCP specification file(s)
! this has been separated off, as we want to call
! twice - once for the UKCA lbc_specs and once
! for the full listing of lbc_specs

! check UKCA values, do lumping as UKCA does lumping
file_header = 'LBC_MMR VALUES REQUESTED FOR USE IN UKCA'
CALL test_scenario_rcp(n_lbc_specs, lbc_specs,                    &
                       perform_lumping, ukca_test_file,           &
                       file_header)


! only write out full values if on diag print status, as well as
! having the logical trap on the call of this routine
IF (PrintStatus >= Prstatus_Diag) THEN
  ALLOCATE(full_specs(1:n_full_specs))
  full_specs(:) =   '          '
  full_specs    = (/'CFCl3     ',                                &
                    'CF2Cl2    ',                                &
                    'CF2ClCFCl2',                                &
                    'CF2ClCF2Cl',                                &
                    'CF2ClCF3  ',                                &
                    'CCl4      ',                                &
                    'MeCCl3    ',                                &
                    'CHF2Cl    ',                                &
                    'MeCF2Cl   ',                                &
                    'MeCFCl2   ',                                &
                    'CF2ClBr   ',                                &
                    'CF2Br2    ',                                &
                    'CF3Br     ',                                &
                    'CF2BrCF2Br',                                &
                    'MeBr      ',                                &
                    'MeCl      ',                                &
                    'N2O       ',                                &
                    'CH4       ',                                &
                    'H2        ',                                &
                    'N2        ',                                &
                    'CO2       ',                                &
                    'CH2Br2    ',                                &
                    'CHBr3     ',                                &
                    'CF3CHF2   ',                                &
                    'CH2FCF3   ',                                &
                    'COS       '/)


  ! check all values - DO NOT LUMP
  file_header = 'LBC_MMR VALUES FOR ALL CONSIDERED SPECIES'
  CALL test_scenario_rcp(n_full_specs, full_specs,               &
                         lumping_off, full_test_file,            &
                         file_header)

  IF (ALLOCATED(full_specs)) DEALLOCATE(full_specs)
END IF ! Prstatus_Diag


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE test_scenario_rcp_ctl


SUBROUTINE test_scenario_rcp(n_lbc_specs, lbc_specs,                &
                             perform_lumping, test_file,            &
                             file_header)
USE ukca_constants,   ONLY: c_cf2cl2,c_cf2clcfcl2, c_cf2clcf2cl,  &
    c_cf2clcf3, c_cfcl3, c_ccl4, c_meccl3, c_chf2cl, c_mecfcl2,   &
    c_mecf2cl, c_cf2clbr, c_mecl, c_mebr, c_cf2br2, c_cf3br,      &
    c_cf2brcf2br, c_n2o, c_ch4, c_co2, c_cf3chf2, c_ch2fcf3,      &
    c_cos, c_ch2br2, c_chbr3, c_h2, c_n2
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE ereport_mod,      ONLY: ereport
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs           ! Number of LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs(n_lbc_specs)! LBC species
LOGICAL, INTENT(IN)      :: perform_lumping       ! T for lumping of lbc's
CHARACTER(LEN=256), INTENT(IN) :: test_file       ! file to write to
CHARACTER(LEN=256), INTENT(IN) :: file_header     ! header info for file

REAL, ALLOCATABLE :: lbc_mmr(:)  ! Lower BC mass mixing ratios

INTEGER :: file_unit
! loop counters
INTEGER :: i
INTEGER :: j
! to simulate many years to calculate the lbc values
INTEGER :: jday
INTEGER :: jyear

REAL, ALLOCATABLE :: lbc_c_species(:)

! FORMAT STRINGS FOR OUTPUT STATEMENTS - NEED TO BE VARIABLE LENGTH
! format of the name of the strings
CHARACTER(LEN=80)  :: tr_name_fmt_str
! format of the conversion factors
CHARACTER(LEN=80)  :: tr_conv_fmt_str
! format of the mmr values
CHARACTER(LEN=80)  :: tr_float_fmt_str
! format to number the lines for ease of reading
CHARACTER(LEN=80)  :: tr_num_fmt_str


! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEST_SCENARIO_RCP'

! Body of subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! this routine test the scenario_rcp routine by
! cycling through all 1700-2550 and writing the
! lbc_mmr values out to a file so that they can be
! plotted and compared to the ones in the original
! RCP specification file(s)

! allocate arrays
ALLOCATE(lbc_mmr(1:n_lbc_specs))
ALLOCATE(lbc_c_species(1:n_lbc_specs))

! populate lbc_c_species for printing purposes
lbc_c_species(:) = 0.0
WHERE (lbc_specs == 'CFCl3     ') lbc_c_species = c_cfcl3
WHERE (lbc_specs == 'CF2Cl2    ') lbc_c_species = c_cf2cl2
WHERE (lbc_specs == 'CF2ClCFCl2') lbc_c_species = c_cf2clcfcl2
WHERE (lbc_specs == 'CF2ClCF2Cl') lbc_c_species = c_cf2clcf2cl
WHERE (lbc_specs == 'CF2ClCF3  ') lbc_c_species = c_cf2clcf3
WHERE (lbc_specs == 'CCl4      ') lbc_c_species = c_ccl4
WHERE (lbc_specs == 'MeCCl3    ') lbc_c_species = c_meccl3
WHERE (lbc_specs == 'CHF2Cl    ') lbc_c_species = c_chf2cl
WHERE (lbc_specs == 'MeCF2Cl   ') lbc_c_species = c_mecf2cl
WHERE (lbc_specs == 'MeCFCl2   ') lbc_c_species = c_mecfcl2
WHERE (lbc_specs == 'CF2ClBr   ') lbc_c_species = c_cf2clbr
WHERE (lbc_specs == 'CF2Br2    ') lbc_c_species = c_cf2br2
WHERE (lbc_specs == 'CF3Br     ') lbc_c_species = c_cf3br
WHERE (lbc_specs == 'CF2BrCF2Br') lbc_c_species = c_cf2brcf2br
WHERE (lbc_specs == 'MeBr      ') lbc_c_species = c_mebr
WHERE (lbc_specs == 'MeCl      ') lbc_c_species = c_mecl
WHERE (lbc_specs == 'N2O       ') lbc_c_species = c_n2o
WHERE (lbc_specs == 'CH4       ') lbc_c_species = c_ch4
WHERE (lbc_specs == 'H2        ') lbc_c_species = c_h2
WHERE (lbc_specs == 'N2        ') lbc_c_species = c_n2
WHERE (lbc_specs == 'CO2       ') lbc_c_species = c_co2
WHERE (lbc_specs == 'CH2Br2    ') lbc_c_species = c_ch2br2
WHERE (lbc_specs == 'CHBr3     ') lbc_c_species = c_chbr3
WHERE (lbc_specs == 'CF3CHF2   ') lbc_c_species = c_cf3chf2
WHERE (lbc_specs == 'CH2FCF3   ') lbc_c_species = c_ch2fcf3
WHERE (lbc_specs == 'COS       ') lbc_c_species = c_cos

! find a file unit to write to
CALL assign_file_unit(TRIM(ADJUSTL(test_file)), &
                      file_unit, handler="fortran")

! open file to write test output to
ierr = 0
OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(test_file)),                    &
     POSITION='rewind',STATUS='unknown',IOSTAT=ierr, IOMSG=iomessage)
IF (ierr /= 0) THEN
  cmessage='Error opening file '//TRIM(ADJUSTL(test_file)) // ' : '// &
           TRIM(iomessage)
  CALL ereport('TEST_SCENARIO_RCP',ierr,cmessage)
END IF

! write some header to the file
! start header lines with #
! generate what the format of the header lines should be first
! e.g. "#",X,A10,X,7(A20,X) etc.
! it is necessary to generate these 'on the fly' as there will be a
! variable number of columns dependent on scheme used, or if the output
! is 'full' or just the values UKCA will see (i.e. with lumping).
WRITE(tr_name_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs,'(A20,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_conv_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs,'(E20.10,1X))'
! e.g. 2X,F12.6,7(E20.10,X) etc.
WRITE(tr_float_fmt_str,FMT='(A,I2,A)')                            &
                     '(2X,F12.6,1X,',n_lbc_specs,'(E20.10,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_num_fmt_str,FMT='(A,I2,A)')                              &
                     '(A1,1X,A12,1X,',n_lbc_specs,'(I20,1X))'

! now write header information to the file
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(test_file))
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(file_header))
WRITE(file_unit,FMT=tr_name_fmt_str) '#','SPECIES',               &
                    (TRIM(ADJUSTL(lbc_specs(i))),i=1,n_lbc_specs)
WRITE(file_unit,FMT=tr_conv_fmt_str) '#','CONV FAC',              &
                    (lbc_c_species(i),i=1,n_lbc_specs)
WRITE(file_unit,FMT=tr_num_fmt_str) '#','YEAR',(i,i=1,n_lbc_specs)

IF (PrintStatus >= Prstatus_Oper) THEN
  WRITE(umMessage,'(A,A)') 'Outputting RCP lower BC data to ',   &
       TRIM(ADJUSTL(test_file))
  CALL umPrint(umMessage,src='test_scenario_rcp')
END IF

! call the scenario_rcp routine, as if cycling through many centuries
! do up to 2550 (i.e. 2549+360 days)
DO jyear=1700,2549
  DO jday=1,360
    CALL  ukca_scenario_rcp(n_lbc_specs, lbc_specs,             &
                    jyear, jday,                                &
                    lbc_mmr, perform_lumping)
    ! output mmr values to our output file
    WRITE(file_unit,FMT=tr_float_fmt_str)                       &
                        (REAL(jyear)+REAL(jday)/360.0),         &
                        (lbc_mmr(i),i=1,n_lbc_specs)
  END DO
END DO

! all finished now - close the file
CLOSE(file_unit)
CALL release_file_unit(file_unit, handler="fortran")

! deallocate arrays
IF (ALLOCATED(lbc_c_species)) DEALLOCATE(lbc_c_species)
IF (ALLOCATED(lbc_mmr)) DEALLOCATE(lbc_mmr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE test_scenario_rcp


END MODULE ukca_scenario_rcp_mod
