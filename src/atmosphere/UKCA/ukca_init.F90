! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: To initialize internal values, addresses and
! other information needed for UKCA
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
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_INIT_MOD'

CONTAINS

SUBROUTINE ukca_init

USE asad_mod,              ONLY: cdt, interval, ncsteps, tslimit
USE ukca_cdnc_mod,         ONLY: ukca_cdnc_init
USE ukca_option_mod,       ONLY: l_ukca, l_ukca_aie1,             &
                                 l_ukca_aie2, check_run_ukca,     &
                                 l_ukca_raq, i_mode_setup,        &
                                 l_ukca_mode, ukca_int_method,    &
                                 l_ukca_plume_scav,chem_timestep, &
                                 i_ukca_chem, l_ukca_ageair
USE ukca_chem_schemes_mod, ONLY: i_ukca_chem_offline_be,          &
                                 i_ukca_chem_off, int_method_nr,  &
                                 int_method_be_explicit,          &
                                 int_method_impact
USE ukca_setup_chem_mod,   ONLY: ukca_setup_chem
USE ukca_d1_defs,          ONLY: ukca_item_sulpc,                 &
                                 n_mode_tracers
USE ukca_mode_setup,       ONLY: nmodes,                          &
                                 mode_choice, ncp, mode,          &
                                 component, ukca_mode_suss_4mode, &
                                 ukca_mode_sussbcoc_5mode,        &
                                 ukca_mode_sussbcoc_4mode,        &
                                 ukca_mode_sussbcocso_5mode,      &
                                 ukca_mode_sussbcocso_4mode,      &
                                 ukca_mode_duonly_2mode,          &
                                 ukca_mode_sussbcocdu_7mode,      &
                                 ukca_mode_sussbcocntnh_5mode_8cpt
USE ukca_setup_indices,    ONLY: ukca_indices_sv1,                &
                                 ukca_indices_suss_4mode,         &
                                 ukca_indices_orgv1_soto3,        &
                                 ukca_indices_sussbcoc_5mode,     &
                                 ukca_indices_orgv1_soto3,        &
                                 ukca_indices_sussbcoc_4mode,     &
                                 ukca_indices_orgv1_soto6,        &
                                 ukca_indices_sussbcocso_5mode,   &
                                 ukca_indices_orgv1_soto6,        &
                                 ukca_indices_sussbcocso_4mode,   &
                                 ukca_indices_nochem,             &
                                 ukca_indices_duonly_2mode,       &
                                 ukca_indices_sussbcocdu_7mode
USE ukca_scavenging_mod, ONLY:   ukca_mode_scavcoeff
USE cv_run_mod,            ONLY: i_convection_vn,                 &
                                 i_convection_vn_5a,              &
                                 i_convection_vn_6a
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE Submodel_Mod
USE umPrintMgr,            ONLY: umPrint, umMessage,              &
                                 PrintStatus, PrStatus_Oper
USE ereport_mod,           ONLY: ereport
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Local variables
INTEGER                       :: imode     ! loop counter for modes
INTEGER                       :: icp       ! loop counter for components
INTEGER                       :: i         ! loop counter
INTEGER                       :: j         ! loop counter
INTEGER                       :: jend      ! loop end counter
INTEGER                       :: k1        ! to hold 100s
INTEGER                       :: k2        ! to hold 10s
INTEGER                       :: k3        ! to hold 1s
INTEGER                       :: n_reqd_tracers ! no. of required tracers
INTEGER                       :: errcode=0     ! Error code: ereport
INTEGER                       :: timestep      ! Dynamical timestep
CHARACTER(LEN=errormessagelength)     :: cmessage=' '  ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check the logical switches
CALL check_run_ukca()

! Set internal UKCA values from UKCA namelists
IF (l_ukca) CALL ukca_setup_chem()

! Set item numbers for CDNC etc
IF (l_ukca_aie1 .OR. l_ukca_aie2) THEN
  CALL ukca_cdnc_init(errcode,cmessage)
END IF

! Set up timestep counting. Interval depends on solver: Backward-Euler
! and other solvers are run every dynamical timestep. Newton-Raphson and
! Offline oxidants (backward-Euler) may be run every 2 / 3 timesteps for
! a 30/20 minutes dynamical timestep. Chemistry timestep (chem_timestep)
! is set in the namelist for these solvers.

! dynamical timestep
timestep = secs_per_periodim(atmos_im) / steps_per_periodim(atmos_im)

! Do not check for solver type and timestep if none of the chemistry 
! schemes is selected e.g. for an Age-of-air-only configuration. 
! In that case, set values to default in case they are used elsewhere
!$OMP PARALLEL
L_ukca_chem_off: IF ( i_ukca_chem == i_ukca_chem_off ) THEN
  interval = 1
  cdt = REAL(timestep)
  ncsteps = 1
ELSE   
  ! Calculate interval depending on chemistry scheme
  ukca_method_type:IF (ukca_int_method == int_method_nr) THEN
    ! Newton-Raphson solver
    interval = chem_timestep/timestep
    cdt = REAL(chem_timestep)
    ncsteps = 1
  ELSE IF (ukca_int_method == int_method_impact) THEN     
    ! IMPACT solver use about 15 or 10 minutes, depending on dynamical timestep
    IF (timestep < tslimit) THEN
      cdt = REAL(timestep)
      ncsteps = 1
      interval = 1
    ELSE
      cdt = 0.5*REAL(timestep)
      ncsteps = 2
      interval = 1
    END IF
  ELSE IF (ukca_int_method == int_method_be_explicit) THEN  
    ! explicit Backward-Euler solver
    ! solver interval derived from namelist value of chemical timestep
    interval = chem_timestep/timestep
    cdt = REAL(chem_timestep)
    ncsteps = 1
  ELSE
    ! Unknown solver type
    WRITE(cmessage, '(A,I0,A)') &
     'Type of solver (ukca_int_method = ',ukca_int_method,') not recognised.'
    errcode = 100
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF ukca_method_type

  IF (printstatus >= prstatus_oper) THEN
    WRITE(umMessage,'(A40,I6)') 'Interval for chemical solver set to: ',    &
                                                               interval
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A40,E12.4)') 'Timestep for chemical solver set to: ', &
                                                               cdt
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A40,I6)') 'No. steps for chemical solver set to: ',   &
                                                               ncsteps
    CALL umPrint(umMessage,src='ukca_init')
  END IF

  ! Verify that the interval and timestep values have been set correctly
  IF (ABS(cdt*ncsteps - REAL(timestep*interval)) > 1e-4) THEN
    cmessage=' chemical timestep does not fit dynamical timestep'
    WRITE(umMessage,'(A)') cmessage
    CALL umPrint(umMessage,src='ukca_init')
    WRITE(umMessage,'(A,I6,A,I6)') ' timestep: ',timestep, &
                                   ' interval: ',interval
    CALL umPrint(umMessage,src='ukca_init')
    errcode = chem_timestep
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF

END IF L_ukca_chem_off 
!$OMP END PARALLEL

IF (l_ukca_mode) THEN

  ! Check that convection scheme version supports plume scavenging
  IF ((i_convection_vn /= i_convection_vn_5a .AND.                 &
       i_convection_vn /= i_convection_vn_6a) .AND. l_ukca_plume_scav) THEN
    cmessage = ' Convective plume scavenging not available, '//    &
               'in this convection Vn., but required for GLOMAP-mode'
    WRITE(umMessage,'(A,A,I6)') 'Convective plume scavenging not available ', &
                             ' in version: ', i_convection_vn
    CALL umPrint(umMessage,src='ukca_init')
    errcode = 1
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF

  ! Call appropriate MODE setup routine
  IF (i_mode_setup == 1) THEN
    CALL ukca_indices_sv1
    CALL ukca_indices_suss_4mode
    CALL ukca_mode_suss_4mode
  ELSE IF (i_mode_setup == 2) THEN
    CALL ukca_indices_orgv1_soto3
    CALL ukca_indices_sussbcoc_5mode
    CALL ukca_mode_sussbcoc_5mode
  ELSE IF (i_mode_setup == 3) THEN
    CALL ukca_indices_orgv1_soto3
    CALL ukca_indices_sussbcoc_4mode
    CALL ukca_mode_sussbcoc_4mode
  ELSE IF (i_mode_setup == 4) THEN
    CALL ukca_indices_orgv1_soto6
    CALL ukca_indices_sussbcocso_5mode
    CALL ukca_mode_sussbcocso_5mode
  ELSE IF (i_mode_setup == 5) THEN
    CALL ukca_indices_orgv1_soto6
    CALL ukca_indices_sussbcocso_4mode
    CALL ukca_mode_sussbcocso_4mode
  ELSE IF (i_mode_setup == 6) THEN
!!    CALL ukca_indices_nochem
!! temporarily run 2-mode dust only with chemistry, though it's not needed
    CALL ukca_indices_orgv1_soto3
    CALL ukca_indices_duonly_2mode
    CALL ukca_mode_duonly_2mode
    !!      ELSE IF(i_mode_setup == 7) THEN
    !!        CALL ukca_indices_nochem
    !!        CALL ukca_indices_duonly_3mode
    !!        CALL ukca_mode_duONLY_3mode
  ELSE IF(i_mode_setup == 8) THEN
    CALL ukca_indices_orgv1_soto3
    CALL ukca_indices_sussbcocdu_7mode
    CALL ukca_mode_sussbcocdu_7mode
    !!      ELSE IF(i_mode_setup == 9) THEN
    !!        CALL ukca_indices_orgv1_soto3
    !!        CALL ukca_indices_sussbcocdu_4mode
    !!        CALL ukca_mode_sussbcocdu_4mode
  ELSE IF(i_mode_setup == 10) THEN
    cmessage=' i_mode_setup 10 needs further development, see UMDP84'
    WRITE(umMessage,'(A,I4)') cmessage,i_mode_setup
    CALL umPrint(umMessage,src='ukca_init')
    errcode = 1
    CALL ereport('UKCA_INIT',errcode,cmessage)
!    CALL ukca_mode_sussbcocntnh_5mode_8cpt
  ELSE
    cmessage=' i_mode_setup has unrecognised value'
    WRITE(umMessage,'(A,I4)') cmessage,i_mode_setup
    CALL umPrint(umMessage,src='ukca_init')
    errcode = 1
    CALL ereport('UKCA_INIT',errcode,cmessage)
  END IF       ! i_mode_setup

  ! Calculate number of aerosol tracers required for components and number
  n_reqd_tracers = 0
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      DO icp=1,ncp
        IF (component(imode,icp)) n_reqd_tracers = n_reqd_tracers + 1
      END DO
    END IF
  END DO
  n_mode_tracers = n_reqd_tracers + SUM(mode_choice)
  
  ! Set up the scavenging coefficients for plume scavenging
  CALL ukca_mode_scavcoeff()
ELSE
  ! Allocate arrays that are referred to outside GLOMAP
  ! to avoid compiler errors
  IF (.NOT. ALLOCATED(component)) ALLOCATE(component(1,1))

END IF    ! L_ukca_mode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_init
END MODULE ukca_init_mod
