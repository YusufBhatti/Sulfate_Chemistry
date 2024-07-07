! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
!  Specify surface boundary conditions
!  - This routine will call one of three others, depending on the value
!    of I_UKCA_SCENARIO which is a namelist input:
!    0 == UKCA_SCENARIO_PRESCRIBED which uses UM values
!    1 == UKCA_SCENARIO_WMOA1 which follows the WMOA1(b) scenario
!    2 == UKCA_SCENARIO_RCP which reads-in the scenario from an input
!         datafile, in a format which corresponds to those which describe
!         the represenative concentration pathways defined for CMIP5.
!  - It is also possible to test the UKCA_SCENARIO_RCP code by setting the
!    logical L_UKCA_TEST_SCENARIO_RCP to .TRUE. in UKCA_SCENARIO_CTL_MOD
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
!  Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
MODULE ukca_scenario_ctl_mod
IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_scenario_ctl

! take LBC values from namelist
INTEGER, PARAMETER, PUBLIC :: i_ukca_scenario_um    = 0
! take LBC values from WMOAl routine
INTEGER, PARAMETER, PUBLIC :: i_ukca_scenario_wmoa1 = 1
! take LBC values from RCP routine
INTEGER, PARAMETER, PUBLIC :: i_ukca_scenario_rcp   = 2

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SCENARIO_CTL_MOD'

CONTAINS

SUBROUTINE ukca_scenario_ctl(n_boundary_vals, lbc_spec,           &
                             i_year, i_day_number,                &
                             lbc_mmr, perform_lumping)

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE ereport_mod,           ONLY: ereport
USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus
USE Control_Max_Sizes
USE UM_ParCore,            ONLY: mype
USE UKCA_Scenario_RCP_mod, ONLY: ukca_scenario_rcp,               &
                                 test_scenario_rcp_ctl,           &
                                 l_ukca_test_scenario_rcp
USE ukca_option_mod,       ONLY: i_ukca_scenario
USE errormessagelength_mod, ONLY: errormessagelength

USE ukca_scenario_prescribed_mod, ONLY: ukca_scenario_prescribed
USE ukca_scenario_wmoa1_mod, ONLY: ukca_scenario_wmoa1
IMPLICIT NONE

! Number of lower BC species
INTEGER, INTENT(IN) :: n_boundary_vals
! LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_spec(n_boundary_vals)
INTEGER, INTENT(IN) :: i_year                    ! model year
INTEGER, INTENT(IN) :: i_day_number              ! model day
! Lower BC mass mixing ratios
REAL, INTENT(INOUT) :: lbc_mmr(n_boundary_vals)
! T for lumping of lower boundary conditions
LOGICAL, INTENT(IN) :: perform_lumping

LOGICAL, SAVE :: L_first=.TRUE.

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage

! loop variable for diagnostic output
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SCENARIO_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (i_ukca_scenario)
CASE (i_ukca_scenario_um)
  IF (PrintStatus >= Prstatus_Diag) THEN
    WRITE(umMessage,'(A)') 'Taking UKCA Lower BC values from UM'
    CALL umPrint(umMessage,src='ukca_scenario_ctl')
  END IF
  ! take the LBC values from how they are set in the UM radiation
  ! routines, but assume these numbers are only valid for the surface
  ! don't need to pass in i_year or i_day_number as the UM deals with
  ! the time evolution of these values
  CALL ukca_scenario_prescribed(n_boundary_vals, lbc_spec,       &
       lbc_mmr, perform_lumping)
CASE (i_ukca_scenario_wmoa1)
  IF (PrintStatus >= Prstatus_Diag) THEN
    WRITE(umMessage,'(A)')'Taking UKCA Lower BC values from WMOA1 routine'
    CALL umPrint(umMessage,src='ukca_scenario_ctl')
  END IF
  ! take LBC values from the WMOA1b specification (2006)
  CALL ukca_scenario_wmoa1(n_boundary_vals, lbc_spec,            &
       i_year, i_day_number,                                     &
       lbc_mmr, perform_lumping)
CASE (i_ukca_scenario_rcp)
  IF (PrintStatus >= Prstatus_Diag) THEN
    WRITE(umMessage,'(A)')'Taking UKCA Lower BC values from RCP routine'
    CALL umPrint(umMessage,src='ukca_scenario_ctl')
  END IF
  ! take LBC values from those contained in a file formatted to
  ! provide this data for CMIP5, which follow various RCP scenarios
  CALL ukca_scenario_rcp(n_boundary_vals, lbc_spec,              &
       i_year, i_day_number,                                     &
       lbc_mmr, perform_lumping)
  ! output diagnostic printout of values if requested (recompilation
  ! will be required due to logical trap)
  IF (L_first) THEN
    IF (mype == 0) THEN
      IF (L_ukca_test_scenario_rcp) THEN
        CALL test_scenario_rcp_ctl(n_boundary_vals, lbc_spec, &
                                   perform_lumping)
      END IF
    END IF
  END IF
CASE DEFAULT
   ! should not have come here, so exit with an error
  IF (i_ukca_scenario < 0) THEN
    cmessage = "ERROR: i_ukca_scenario is unset"
  ELSE
    cmessage = "ERROR: incorrect scenario choice for "//        &
               "i_ukca_scenario"
  END IF
  ! make sure number is +ve (cannot be 0 as this is an option in the
  ! CASE) as this needs to be an error and not a warning
  ierr = ABS(i_ukca_scenario)
  CALL ereport('UKCA_SCENARIO_CTL',ierr,cmessage)
END SELECT

! diagnostic output of the LBC values at this timestep (if requested)
IF (PrintStatus >= Prstatus_Diag) THEN
  DO i=1,n_boundary_vals
    WRITE(umMessage,'(A15,1X,A10,1X,A1,1X,E18.8)') ' UKCA Lower BC:',&
         TRIM(ADJUSTL(lbc_spec(i))),'=',lbc_mmr(i)
    CALL umPrint(umMessage,src='ukca_scenario_ctl')
  END DO
END IF

L_first = .FALSE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_scenario_ctl

END MODULE ukca_scenario_ctl_mod
