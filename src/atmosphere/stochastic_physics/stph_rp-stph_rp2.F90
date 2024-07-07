! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Stochastic Physics (sect35) Random Parameters Ver. 2
MODULE stph_rp2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_RP2_MOD'

CONTAINS


SUBROUTINE stph_rp2(model_levels,                                       &
                    rhcrit, rhcrit_max, rhcrit_min,                     &
                    m_ci, m_ci_max, m_ci_min,                           &
                    Charnock)

USE bl_option_mod, ONLY: i_bl_vn,                                       &
                         i_bl_vn_9b,                                    &
                         i_bl_vn_9c

USE stochastic_physics_run_mod, ONLY:                                   &
    stphseed, ran_max, rhcrit_ref_level, ran_count,                     &
    rp2_callfreq, i_rp_scheme, i_rp2, i_rp2b,                           &
    l_rp2_cycle_out, l_rp2_cycle_in, rp2_cycle_tm,                      &
    x1r_rp, x1r_rp_min, x1r_rp_max,                                     &
    ndrop_surf_rp, ndrop_surf_rp_min, ndrop_surf_rp_max,                &
    ec_auto_rp, ec_auto_rp_min, ec_auto_rp_max,                         &
    g0_rp, g0_rp_max, g0_rp_min,                                        &
    par_mezcla, par_mezcla_max, par_mezcla_min,                         &
    charnock_max, charnock_min,                                         &
    lambda_min_rp_max, lambda_min_rp, lambda_min_rp_min,                &
    ricrit_rp, ricrit_rp_max, ricrit_rp_min,                            &
    a_ent_1_rp, a_ent_1_rp_max, a_ent_1_rp_min,                         &
    g1_rp, g1_rp_max, g1_rp_min,                                        &
    lam_meta_rp, lam_meta_rp_max, lam_meta_rp_min,                      &
    stph_rp2_data_present, stph_rp2_data_check, rp_max,                 &
    lai_mult_rp, lai_mult_rp_max, lai_mult_rp_min,                      &
    dz0v_dh_rp, dz0v_dh_rp_max, dz0v_dh_rp_min,                         &
    z0hm_pft_rp, z0hm_pft_rp_max, z0hm_pft_rp_min


USE nlcfiles_namelist_mod, ONLY: rp2_seed

USE dump_headers_mod, ONLY: fdc_rp2_coef_start, a_flddepc


USE timestep_mod,      ONLY: timestep_number, timestep

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,           ONLY: mype, nproc

USE nlstcall_mod, ONLY: ldump
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE stph_rp_pert_mod, ONLY: stph_rp_pert
USE stph_closeinput_mod,  ONLY: stph_closeinput
USE stph_closeoutput_mod, ONLY: stph_closeoutput
USE stph_openinput_mod,   ONLY: stph_openinput
USE stph_openoutput_mod,  ONLY: stph_openoutput
USE stph_readentry_mod,   ONLY: stph_readentry
USE stph_writeentry_mod,  ONLY: stph_writeentry

! Include switch for unbiased RHCRIT perturbations
USE science_fixes_mod,    ONLY: l_stph_rhcrit_unbias

USE jules_surface_types_mod, ONLY: npft
USE nlsizes_namelist_mod, ONLY: land_field
USE max_dimensions, ONLY: npft_max

IMPLICIT NONE

!
! Description:  Introduce stochastic perturbation in some physics params
!               to take into account model error
!
! Method:       The value of a given PARAMETER is chosen randomly
!               between a maximum and minimum values. PARAMETER's values
!               are temporally correlated (1st order markov process)
!
! Code Owner:   Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
!
! Code Description:
!   Language:   FORTRAN 90
!   This code is written to UMDP3 version 8.3 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable  !Description of variable
!
! Global variables (#include statements etc):

!-------------------------------------------------------------
! Variable definition
!-------------------------------------------------------------
!
! IN variables
!
INTEGER, INTENT(IN) :: model_levels
!
! Default, Maximum and minimum values for the STPH_RP scheme
! Large Scale Precipitation
! Note: rhcrit scalars refer to reference level=rhcrit_ref_level
!
REAL, INTENT(IN) :: rhcrit_max        ! Max value critical rh
REAL, INTENT(IN) :: rhcrit_min        ! Min value critical rh
REAL, SAVE       :: rhcrit_0          ! Def value critical rh
REAL, SAVE       :: rhcrit_p          ! Perturbed critical rh
REAL, INTENT(IN) :: m_ci_max          ! Max value of multiplication
                                      ! factor for CI
REAL, INTENT(IN) :: m_ci_min          ! Min value of multiplication
                                      ! factor for CI
REAL, SAVE       :: m_ci_0            ! Def value of multiplication
                                      ! factor for CI
REAL, SAVE       :: x1r_rp_0          ! Def value of rain particle size
                                      ! distribution
REAL, SAVE       :: ndrop_surf_rp_0   ! Def value of surface droplet
                                      ! number 
REAL, SAVE       :: ec_auto_rp_0      ! Def value for autoconversion of
                                      ! cloud water to rain

!
! Boundary Layer
!
REAL, SAVE       :: par_mezcla_0      ! Def value neutral mixing length
                                      ! parameter
REAL, SAVE       :: g0_rp_0           ! Def value of stability functions
REAL, SAVE       :: lambda_min_rp_0   ! Def value of min mixing length
REAL, SAVE       :: ricrit_rp_0       ! Def value of critical Ri
REAL, SAVE       :: a_ent_1_rp_0      ! Def value of entrainment A1
REAL, SAVE       :: g1_rp_0           ! Def value of velocity scale
REAL, SAVE       :: charnock_0        ! Def value Charnock param.
REAL, SAVE       :: lam_meta_rp_0     ! Def value of lam_meta_rp, replaces
                                      ! par_mezcla and lambda_min in RP2b
!
! Land-sfc
!
REAL, SAVE       :: lai_mult_rp_0(npft_max)   ! Def value of LAI multiplier
REAL, SAVE       :: rnumb_lai_mult_rp         ! Random number for LAI multiplier
REAL, SAVE       :: rnumb_lai_mult_rp_0 = 0.0 ! Def value of random number 
                                              ! for LAI multiplier
REAL, SAVE       :: dz0v_dh_rp_0(npft_max)    ! Def value of JULES parameter 
                                              ! dz0v_dh
REAL, SAVE       :: rnumb_dz0v_dh_rp          ! Random number for JULES 
                                              ! parameter dz0v_dh
REAL, SAVE       :: rnumb_dz0v_dh_rp_0 = 0.0  ! Def value of random number for 
                                              ! JULES parameter dz0v_dh
REAL, SAVE       :: z0hm_pft_rp_0(npft_max)   ! Def value of JULES parameter
                                              ! z0hm_pft
REAL, SAVE       :: rnumb_z0hm_pft_rp         ! Random number for JULES 
                                              ! parameter z0hm_pft
REAL, SAVE       :: rnumb_z0hm_pft_rp_0 = 0.0 ! Def value of random number
                                              ! for JULES parameter z0hm_pft

!
! IN/OUT variables
!
REAL, INTENT(INOUT) :: rhcrit(model_levels)
                                      ! Crit relative humidity for
                                      ! layer cloud formation
REAL, INTENT(INOUT) :: m_ci           ! Fall-speed in 3C/D microphysics
REAL, INTENT(INOUT) :: Charnock       ! Charnock Parameter


! LOCAL VARIABLES
!
! Local variables to assign values to rhcrit, which is level
! dependant
INTEGER, ALLOCATABLE,SAVE :: interval(:)
REAL,    ALLOCATABLE,SAVE :: drhcrit(:)
                                      ! used to store difference of
                                      ! rhcrit to reference level value
! Local array parameters

! Other local variables
INTEGER, SAVE   :: max_interval
LOGICAL, SAVE   :: l_firstcall     = .TRUE.
LOGICAL, SAVE   :: l_crunfirstcall = .FALSE.
                                      ! 1st call to RP?
LOGICAL, PARAMETER  :: lp1stcll_off = .FALSE.
LOGICAL         :: local_switch       ! Logical to calculate the first
                                      ! new rhcrit value.
INTEGER         :: i, j, rp_count     ! loop variables
INTEGER         :: icode              ! Return code for error reporting
INTEGER         :: lead_time          ! Current forecast lead-time in seconds
REAL, SAVE      :: rp_coef(rp_max,1)  ! Current RP values
                                      ! (2nd dim for compatibility
                                      !  with stph_readentry)
! True if the coefficients have been read in
LOGICAL         :: coefficients_input 
                                      

! Variables associated with random number generator
REAL            :: rp_rand(ran_max)   ! random number for each parameter
REAL            :: rp_rand0(ran_max)  ! random number for initial values
INTEGER         :: istat              ! To record call in the synch.

CHARACTER(LEN=errormessagelength) :: cmessage
             ! OUT error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STPH_RP2'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lead_time = timestep * (timestep_number - 1)
IF (MOD(lead_time, rp2_callfreq) == 0) THEN

  IF (printstatus  >=   prstatus_normal) THEN
    WRITE(umMessage,FMT='(A)') 'CALLING RANDOM PARAMETERS2'
    CALL umPrint(umMessage,src='stph_rp2')
  END IF

ELSE

  IF (printstatus  >=  prstatus_normal) THEN
    WRITE(umMessage,FMT='(A)') 'NOT CALLING RANDOM PARAMETERS2'
    CALL umPrint(umMessage,src='stph_rp2')
    WRITE(umMessage,FMT='(A,I6,A)')                                     &
      'This routine is only called every ', rp2_callfreq, ' seconds'
    CALL umPrint(umMessage,src='stph_rp2')
  END IF
  ! Write current RP perturbed values to dump even if
  ! RP2 not active this cycle before exiting the subroutine.
  IF (ldump) THEN
    IF (mype == 0) THEN
      ! copy rp_coef to dump header array for writing out in dump file
      WRITE(umMessage,FMT='(A)') 'Copying rp2 to dump headers'
      CALL umPrint(umMessage,src='stph_rp2')
      DO i=1,rp_max
        a_flddepc(fdc_rp2_coef_start+i)=rp_coef(i,1)
      END DO
    END IF
    
  END IF
  
  IF (l_rp2_cycle_out) THEN
  ! If l_rp2_cycle_out = TRUE, write current RP perturbed values
  ! to a file at lead_time == rp_cycle_tm seconds, and tag with 
  ! RP2_CYCLE, even if RP2 is not active this cycle.
  ! These RP2_CYCLE files are used to cycle the RP values between 
  ! forecasts of the same ensemble member set, e.g. if members 1 to 11
  ! are run at 03z and 15z, then T+12 RP values of the 03z run can be 
  ! used to intialise the 15z run.
    IF (lead_time == rp2_cycle_tm ) THEN
      IF (mype == 0) THEN
        CALL stph_openoutput("rewind    ", "RP2_CYCLE")
        CALL stph_writeentry(rp_coef, SIZE(rp_coef))
        CALL stph_closeoutput()
      END IF
      l_rp2_cycle_out = .FALSE.  ! Terminate cycle out condition
    END IF
  END IF
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN

END IF

! initialise value to false, if the coefficients are read in this will be 
! set to true
coefficients_input = .FALSE.

IF (l_firstcall) THEN

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Copy RP2 values from dump header:
  ! ---------------------------------
  ! This section of code is called the first time the UM executes. Depending on 
  ! the timestep value, we could be starting an NRUN or CRUN. If the timestep 
  ! value is 1, then we are in an NRUN. If the timestep value is greater than 1,
  ! then we are in a CRUN. There are 3 cases where we want to retrieve the 
  ! stochastic physics values from the dump header. The first is if this is an
  ! NRUN and the stphseed value is set to 2. The second is all
  ! CRUNs. The third is if you are running an NRUN as if it were a
  ! CRUN (using the l_nrun_as_crun logical). If any of those conditions is met, 
  ! then we also check that the data is actually in the dump. If these checks
  ! are passed, then we proceed with copying the RP2 values from the additional 
  ! parameters section of the dump header to the relevant location here.
  
  IF ( ((stphseed == 2 .OR. l_nrun_as_crun) .AND.                       &
        stph_rp2_data_check == stph_rp2_data_present .AND.              &
        timestep_number == 1)                                           &
        .OR. ( timestep_number > 1 .AND.                                &
        stph_rp2_data_check == stph_rp2_data_present ) ) THEN
    l_crunfirstcall = timestep_number > 1 .AND.                         &
        stph_rp2_data_check == stph_rp2_data_present
    coefficients_input = .TRUE.  

    IF (l_nrun_as_crun) THEN
        cmessage = 'l_nrun_as_crun is TRUE:' //  &
           'Reading random numbers from the dump (equivalent of stphseed=2)'
        icode = -1
        CALL ereport(RoutineName, icode, cmessage)
    END IF
                    
    IF (mype == 0) THEN
      DO i=1,rp_max
        rp_coef(i,1) = a_flddepc(fdc_rp2_coef_start+i)
      END DO
    END IF
    ! Broadcast restart coeffs to all processors
    CALL gc_rbcast(3249, rp_max, 0, nproc, istat, rp_coef)

  ELSE IF (l_rp2_cycle_in) THEN
    coefficients_input = .TRUE.                      
    ! If l_rp2_cycle_in = TRUE then the random parameters are read 
    ! in from a file the first time the routine is called.
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("*** CYCLING IN RP PARAMETERS ***")')
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF
    IF (mype == 0) THEN
      WRITE(umMessage, '(A,A)')                                       &
        'calling stph_openinput with filename: ',rp2_seed
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
      CALL stph_openinput("RP2_CYCLE_IN", rp2_seed)
      CALL stph_readentry(rp_coef, SIZE(rp_coef))
      CALL stph_closeinput()
    END IF

    ! Broadcast restart coeffs to all processors
    CALL gc_rbcast(3249, rp_max, 0, nproc, istat, rp_coef)

  END IF
  
  IF (l_stph_rhcrit_unbias) THEN
    ! ---------------------------------------
    ! Allocate drhcrit the first time STPH_RP is called only
    ! This will be used when perturbing rhcrit.
    ! It is a simplified method of perturbing the column
    !  and removes a bias introduced by the original method
    ! ---------------------------------------
    IF (.NOT. ALLOCATED(drhcrit)) THEN
      ALLOCATE (drhcrit(model_levels))
      drhcrit(:) = 0
    END IF

    DO i = rhcrit_ref_level + 1, model_levels
      drhcrit(i) = rhcrit(i) - rhcrit(rhcrit_ref_level)
    END DO
    rhcrit_p = rhcrit(rhcrit_ref_level)       ! critical RH ref value 

  ELSE
    ! ---------------------------------------
    ! Allocate interval the first time STPH_RP is called only
    ! This will be used when perturbing rhcrit.
    ! ---------------------------------------
    IF (.NOT.ALLOCATED(interval)) THEN
      ALLOCATE (interval(model_levels))
      interval(:) = 0
    END IF
    IF (.NOT.ALLOCATED(drhcrit)) THEN
      ALLOCATE (drhcrit(model_levels))
      drhcrit(:) = 0
    END IF
    j=1
    DO i = rhcrit_ref_level, model_levels
      IF(rhcrit(i) <  rhcrit(i-1)) THEN
        drhcrit(i) = rhcrit(i-1) - rhcrit(i)
        interval(i) = j
        j = j + 1
      END IF
    ENDDO
    max_interval = j-1
!
    DO i = rhcrit_ref_level + 1, model_levels
      IF(interval(i) <  interval(i-1)) THEN
        interval(i)=interval(i-1)
      END IF
      IF(drhcrit(i) <  drhcrit(i-1)) THEN
        drhcrit(i)=drhcrit(i-1)
      END IF
    ENDDO
    rhcrit_0 = rhcrit(rhcrit_ref_level)       ! critical RH mean value
  END IF

END IF

! If this is a CRUN, or an NRUN where we have read in the coefficients from 
! a file  then l_firstcall is no longer needed because we use the values
!  read in.
IF (coefficients_input) THEN
  l_firstcall = .FALSE.
END IF

! End of setup at start of NRUN and CRUN
! --------------------------------------------------------------------


! Setup random numbers for initial perturbed values at start of NRUN
IF (l_firstcall) THEN
  IF (mype == 0) THEN
    CALL RANDOM_NUMBER(rp_rand0)
  END IF
  CALL gc_rbcast(3145, ran_max, 0, nproc, istat, rp_rand0)
END IF

!----------------------------------------------------------------------
! Each time random parameters is called we broadcast the random number
! to all processors
!----------------------------------------------------------------------
IF (mype == 0) THEN
  CALL RANDOM_NUMBER(rp_rand)
END IF
CALL gc_rbcast(3145, ran_max, 0, nproc, istat, rp_rand)
IF (printstatus  >  prstatus_normal) THEN
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  WRITE(umMessage,'("*** STPH_RP2 NEW PARAM VALUES ***")')
  CALL umPrint(umMessage,src='stph_rp-stph_rp2')
END IF

!----------------------------------------------------------------------
! Method: Each parameter is individually controlled, and is updated
!         by an AR1 process that uses a random number and the params
!         previous value.
!         At the start of an NRUN a start value is randomly set using
!         rp_rand0 and the default is set inside "stph_rp_pert".
!         At the start of a CRUN, the previous value is read in from
!         a file and the default is set (outside of "stph_rp_pert").
!----------------------------------------------------------------------

! Initialise rp_rand and rp_coef indices
ran_count = 1
rp_count = 1


!
! Large-scale Precip
!
! Perturb m_ci  (Ice fall speed) - 3C/D Microphysics LSPCON3C
IF ( coefficients_input ) THEN
  m_ci_0 = m_ci
  m_ci = rp_coef(rp_count,1)
END IF
CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, m_ci_0,              &
                   m_ci_max, m_ci_min, m_ci)
rp_coef(rp_count,1) = m_ci
rp_count = rp_count + 1
IF (printstatus  >  prstatus_normal) THEN
  WRITE(umMessage,'("M_CI ............. ",2ES16.8)') m_ci, m_ci_0
  CALL umPrint(umMessage,src='stph_rp-stph_rp2')
END IF
! 
! Additional large-scale precip parameters for RP2b
!
IF ( i_rp_scheme == i_rp2b ) THEN

  ! Perturb x1r_rp (rain particle size dist) - 3C Microphysics LSPCON3C
  IF ( coefficients_input ) THEN
    x1r_rp_0 = x1r_rp
    x1r_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, x1r_rp_0,          &
                     x1r_rp_max, x1r_rp_min, x1r_rp)
  rp_coef(rp_count,1) = x1r_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("X1R_RP ........... ",2ES16.8)')                  &
          x1r_rp, x1r_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

  !  Perturb ndrop_surf_rp (surface droplet number) 
  !  - Microphysics lsp_taper_ndrop
  IF ( coefficients_input ) THEN
  
    ndrop_surf_rp_0 = ndrop_surf_rp
    ndrop_surf_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, ndrop_surf_rp_0,   &
                     ndrop_surf_rp_max, ndrop_surf_rp_min,              &
                     ndrop_surf_rp)
  rp_coef(rp_count,1) = ndrop_surf_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("NDROP_SURF_RP .... ",2ES16.8)')                  &
          ndrop_surf_rp, ndrop_surf_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

  !  Perturb ec_auto_rp  (autoconversion cloud water to rain) 
  !  - 3C Microphysics LSPCON3C
  IF ( coefficients_input ) THEN
    ec_auto_rp_0 = ec_auto_rp
    ec_auto_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, ec_auto_rp_0,      &
                      ec_auto_rp_max, ec_auto_rp_min, ec_auto_rp)
  rp_coef(rp_count,1) = ec_auto_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("EC_AUTO_RP ....... ",2ES16.8)')                  &
          ec_auto_rp, ec_auto_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

END IF ! End of i_rp_scheme == i_rp2b

!
! Boundary Layer
!
IF ( i_bl_vn == i_bl_vn_9b .OR. &
     i_bl_vn == i_bl_vn_9c      ) THEN

  IF ( i_rp_scheme == i_rp2 ) THEN
    !   Perturb par_mezcla (neutral mixing length) - EXCOEF (PBL)
    IF ( coefficients_input ) THEN
      par_mezcla_0 = par_mezcla
      par_mezcla = rp_coef(rp_count,1)
    END IF
    CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, par_mezcla_0,    &
                       par_mezcla_max, par_mezcla_min, par_mezcla)
    rp_coef(rp_count,1) = par_mezcla
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("PAR_MEZCLA ....... ",2ES16.8)')                &
                        par_mezcla, par_mezcla_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF
  END IF

  !   Perturb g0_rp (stability function parameter) - EXCOEF (PBL)
  IF ( coefficients_input ) THEN
    g0_rp_0 = g0_rp
    g0_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, g0_rp_0,           &
                     g0_rp_max, g0_rp_min, g0_rp)
  rp_coef(rp_count,1) = g0_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("G0_RP ............ ",2ES16.8)') g0_rp, g0_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

  IF ( i_rp_scheme == i_rp2 ) THEN

    !   Perturb lambda_min_rp (minimum mixing length) - (PBL)
    IF ( coefficients_input ) THEN
      lambda_min_rp_0 = lambda_min_rp
      lambda_min_rp = rp_coef(rp_count,1)
    END IF
    CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, lambda_min_rp_0, &
              lambda_min_rp_max, lambda_min_rp_min, lambda_min_rp)
    rp_coef(rp_count,1) = lambda_min_rp
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("LAMBDA_MIN_RP..... ",2ES16.8)') lambda_min_rp, &
                                                         lambda_min_rp_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF

    !   Perturb ricrit_rp (critical Ri) - (PBL)
    IF ( coefficients_input ) THEN
      ricrit_rp_0 = ricrit_rp
      ricrit_rp = rp_coef(rp_count,1)
    END IF
    CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, ricrit_rp_0,     &
                       ricrit_rp_max, ricrit_rp_min, ricrit_rp)
    rp_coef(rp_count,1) = ricrit_rp
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("RICRIT_RP ........ ",2ES16.8)')                &
                                       ricrit_rp, ricrit_rp_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF

  ELSE ! Combine parameters to vary together for i_rp_scheme == i_rp2b
    
    ! Perturb new parameter lam_meta_rp (mixing length) - PBL
    IF ( coefficients_input ) THEN
      lam_meta_rp_0 = lam_meta_rp
      lam_meta_rp = rp_coef(rp_count,1)
    END IF
    CALL stph_rp_pert(l_firstcall, rp_rand0, rp_rand, lam_meta_rp_0,    &
                       lam_meta_rp_max, lam_meta_rp_min, lam_meta_rp)
    rp_coef(rp_count,1) = lam_meta_rp
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("LAM_META_RP ...... ",2ES16.8)')                &
                        lam_meta_rp, lam_meta_rp_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF 

    ! Set par_mezcla (neutral mixing length) as a function of 
    ! lam_meta_rp - EXCOEF (PBL)
    IF ( coefficients_input ) THEN
      par_mezcla_0 = par_mezcla
      par_mezcla = rp_coef(rp_count,1)
    END IF
    IF (l_firstcall) par_mezcla_0 = par_mezcla 
    par_mezcla = lam_meta_rp * par_mezcla_0 
    rp_coef(rp_count,1) = par_mezcla
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("PAR_MEZCLA ....... ",2ES16.8)')                &
                        par_mezcla, par_mezcla_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF

    ! Set lambda_min_rp (minimum mixing length) as a function of
    ! lam_meta_rp - (PBL)
    IF ( coefficients_input ) THEN
      lambda_min_rp_0 = lambda_min_rp
      lambda_min_rp = rp_coef(rp_count,1)
    END IF
    IF (l_firstcall) lambda_min_rp_0 = lambda_min_rp 
    lambda_min_rp = lam_meta_rp * lambda_min_rp_0 
    rp_coef(rp_count,1) = lambda_min_rp
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("LAMBDA_MIN_RP..... ",2ES16.8)')                &
            lambda_min_rp, lambda_min_rp_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF       

    ! Set ricrit_rp (critical Ri) as a function of g0_rp - (PBL)
    IF ( coefficients_input ) THEN
      ricrit_rp_0 = ricrit_rp
      ricrit_rp = rp_coef(rp_count,1)
    END IF
    IF ( l_firstcall) ricrit_rp_0 = ricrit_rp
    ricrit_rp = 10.0*ricrit_rp_0/g0_rp
    rp_coef(rp_count,1) = ricrit_rp
    rp_count = rp_count + 1
    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'("RICRIT_RP ........ ",2ES16.8)')                &
                                       ricrit_rp, ricrit_rp_0
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END IF      

  END IF ! i_rp_scheme

  !  Perturb a_ent_1_rp (entrainment parameter A1) - (PBL)
  IF ( coefficients_input ) THEN
    a_ent_1_rp_0 = a_ent_1_rp
    a_ent_1_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, a_ent_1_rp_0,      &
                     a_ent_1_rp_max, a_ent_1_rp_min, a_ent_1_rp)
  rp_coef(rp_count,1) = a_ent_1_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("A_ENT_1_RP ....... ",2ES16.8)')                  &
            a_ent_1_rp, a_ent_1_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

  !   Perturb g1_rp (velocity scale parameter g1) - (PBL)
  IF ( coefficients_input ) THEN
    g1_rp_0 = g1_rp
    g1_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, g1_rp_0,           &
                     g1_rp_max, g1_rp_min, g1_rp)
  rp_coef(rp_count,1) = g1_rp
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("G1_RP ............ ",2ES16.8)') g1_rp, g1_rp_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF

  !   Perturb charnock (Charnock parameter) - (PBL)
  IF ( coefficients_input ) THEN
    charnock_0 = charnock
    charnock = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, charnock_0,        &
                     charnock_max, charnock_min, charnock)
  rp_coef(rp_count,1) = charnock
  rp_count = rp_count + 1
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("CHARNOCK ......... ",2ES16.8)')                  &
                                   charnock, charnock_0
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END IF
END IF

!
! Land-surface parameters
!
! Parameter values are input as an array for each plant function 
! type (1:npft).  To ensure the value for each plant function
! type (PFT) is changing in the same way (for example, making the
! forecast smoother /  rougher overall) we use the same random number
! to perturb each PFT for a given parameter.  This random number, 
! rather than the parameter values, is stored in the rp_coef array.
IF ( i_rp_scheme == i_rp2b ) THEN

  ! Perturb LAI multiplier
  IF ( coefficients_input ) THEN
    rnumb_lai_mult_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert(l_firstcall, rp_rand0, rp_rand, rnumb_lai_mult_rp_0,&
                    1.0, -1.0, rnumb_lai_mult_rp)
  rp_coef(rp_count,1) = rnumb_lai_mult_rp
  rp_count = rp_count + 1

  ! Set default values
  IF ( coefficients_input .OR. l_firstcall) THEN
    lai_mult_rp_0(:) = lai_mult_rp(:)
  END IF

  ! Update value of lai_mult_rp
  IF ( rnumb_lai_mult_rp < 0.0 ) THEN
    DO i = 1, npft
      lai_mult_rp(i) = lai_mult_rp_0(i)                                 &
        + (lai_mult_rp_0(i) - lai_mult_rp_min(i)) * rnumb_lai_mult_rp
    END DO
  ELSE
    DO i = 1, npft
      lai_mult_rp(i) = lai_mult_rp_0(i)                                 &
        + (lai_mult_rp_max(i) - lai_mult_rp_0(i)) * rnumb_lai_mult_rp
    END DO
  END IF

  IF (printstatus  >=  prstatus_oper) THEN
    DO i = 1, npft
      WRITE(umMessage,'("LAI_MULT_RP(PFT) ....... ",I0,2ES16.8)')        &
                        i,lai_mult_rp(i), lai_mult_rp_0(i)
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END DO
  END IF 

  ! Perturb JULES parameter dz0v_dh
  IF ( coefficients_input ) THEN
    rnumb_dz0v_dh_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert(l_firstcall, rp_rand0, rp_rand, rnumb_dz0v_dh_rp_0, &
                    1.0, -1.0, rnumb_dz0v_dh_rp)
  rp_coef(rp_count,1) = rnumb_dz0v_dh_rp
  rp_count = rp_count + 1

  ! Set default values
  IF ( coefficients_input .OR. l_firstcall) THEN
    dz0v_dh_rp_0(:) = dz0v_dh_rp(:)
  END IF

  ! Update value of dz0v_dh_rp
  IF ( rnumb_dz0v_dh_rp < 0.0 ) THEN
    DO i = 1, npft
      dz0v_dh_rp(i) = dz0v_dh_rp_0(i)                                   &
        + (dz0v_dh_rp_0(i) - dz0v_dh_rp_min(i)) * rnumb_dz0v_dh_rp
    END DO
  ELSE
    DO i = 1, npft
      dz0v_dh_rp(i) = dz0v_dh_rp_0(i)                                   &
        + (dz0v_dh_rp_max(i) - dz0v_dh_rp_0(i)) * rnumb_dz0v_dh_rp
    END DO
  END IF

  IF (printstatus  >=  prstatus_oper) THEN
    DO i=1, npft
      WRITE(umMessage,'("DZ0V_DH_RP(PFT) ....... ",I0,2ES16.8)')         &
                        i,dz0v_dh_rp(i), dz0v_dh_rp_0(i)
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END DO
  END IF 

  ! Perturb JULES parameter z0hm_pft
  IF ( coefficients_input ) THEN
    rnumb_z0hm_pft_rp = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert(l_firstcall, rp_rand0, rp_rand, rnumb_z0hm_pft_rp_0,&
                    1.0, -1.0, rnumb_z0hm_pft_rp)
  rp_coef(rp_count,1) = rnumb_z0hm_pft_rp
  rp_count = rp_count + 1

  ! Set default values
  IF ( coefficients_input .OR. l_firstcall) THEN
    z0hm_pft_rp_0(:) = z0hm_pft_rp(:)
  END IF

  ! Update value of z0hm_pft_rp
  IF ( rnumb_z0hm_pft_rp < 0.0 ) THEN
    DO i = 1, npft
      z0hm_pft_rp(i) = z0hm_pft_rp_0(i)                                 &
        + (z0hm_pft_rp_0(i) - z0hm_pft_rp_min(i)) * rnumb_z0hm_pft_rp
    END DO
  ELSE
    DO i = 1, npft
      z0hm_pft_rp(i) = z0hm_pft_rp_0(i)                                 &
        + (z0hm_pft_rp_max(i) - z0hm_pft_rp_0(i)) * rnumb_z0hm_pft_rp
    END DO
  END IF

  IF (printstatus  >=  prstatus_oper) THEN
    DO i = 1, npft
      WRITE(umMessage,'("Z0HM_PFT_RP(PFT) ....... ",I0,2ES16.8)')        &
                        i,z0hm_pft_rp(i), z0hm_pft_rp_0(i)
      CALL umPrint(umMessage,src='stph_rp-stph_rp2')
    END DO
  END IF 

END IF ! i_rp_scheme

! Perturb rhcrit (critical relative humidity) - (LSP)
IF (l_stph_rhcrit_unbias) THEN
  IF ( coefficients_input ) THEN
    rhcrit_0 = rhcrit_p
    rhcrit_p = rp_coef(rp_count,1)
  END IF
  CALL stph_rp_pert( l_firstcall, rp_rand0, rp_rand, rhcrit_0,          &
                     rhcrit_max, rhcrit_min, rhcrit_p)
  ! Prevent perturbations setting RHCRIT(ref_lev) > RHCRIT(ref_lev - 1)
  ! The profile should decrease monotonically with height
  rhcrit_p = MIN(rhcrit_p, rhcrit(rhcrit_ref_level - 1))
  DO i = rhcrit_ref_level + 1, model_levels
    rhcrit(i) = rhcrit_p + drhcrit(i)
  END DO
  rp_coef(rp_count,1) = rhcrit_p
  rp_count = rp_count + 1
  rhcrit(rhcrit_ref_level) = rhcrit_p

ELSE
  IF (l_firstcall) rhcrit_p = rp_rand0(ran_count) *                     &
                             (rhcrit_max-rhcrit_min) + rhcrit_min
!
  local_switch = .TRUE.
  DO j=1,max_interval
    DO i=1,model_levels
      IF (j == 1) THEN
        IF (interval(i) == j) THEN
          IF (local_switch) THEN
            CALL stph_rp_pert( lp1stcll_off, rp_rand, rp_rand, rhcrit_0,     &
                               rhcrit_max, rhcrit_min, rhcrit_p)
            rhcrit(i) = rhcrit_p
            local_switch = .FALSE.
          ELSE
            rhcrit(i) = rhcrit(i-1)
          END IF
        END IF
      ELSE
        IF (interval(i) == j) THEN
          rhcrit(i) = rhcrit_p - (drhcrit(i)*(j-1))
        END IF
      END IF
    ENDDO
  ENDDO
END IF
IF (printstatus  >  prstatus_normal) THEN
  WRITE(umMessage,'("RHCRIT(4,7,10) ... ",4ES16.8)') rhcrit(4),         &
      rhcrit(7), rhcrit(10), rhcrit_0
  CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  DO i=1,model_levels
    WRITE(umMessage,'("RHCRIT(all) ... ",ES16.8)') rhcrit(i)
    CALL umPrint(umMessage,src='stph_rp-stph_rp2')
  END DO
END IF

! Terminate all firstcall conditions
l_firstcall = .FALSE.
l_crunfirstcall = .FALSE.
l_rp2_cycle_in = .FALSE.

! Write current RP perturbed values to file at dump times.
! Generally used for CRUNs.
IF (ldump) THEN
  IF (mype == 0) THEN
    ! copy rp_coef to dump header array for writing out in dump file
    CALL umPrint('Copying rp2 to dump headers', src='stph_rp2')
    DO i=1,rp_max
      a_flddepc(fdc_rp2_coef_start+i)=rp_coef(i,1)
    END DO
  END IF
  
END IF

IF (l_rp2_cycle_out) THEN
  ! Write current RP perturbed values to a file at 
  ! lead_time == rp2_cycle_tm seconds, and tag with RP2_CYCLE. 
  ! These RP2_CYCLE files are used to cycle the RP values between 
  ! forecasts of the same ensemble member set, e.g. if members 1 to 11
  ! are run at 03z and 15z, then T+12 RP values of the 03z run can be 
  ! used to intialise the 15z run.
  IF (lead_time == rp2_cycle_tm) THEN
    IF (mype == 0) THEN
      CALL stph_openoutput("rewind    ", "RP2_CYCLE")
      CALL stph_writeentry(rp_coef, SIZE(rp_coef))
      CALL stph_closeoutput()
    END IF
    l_rp2_cycle_out = .FALSE.  ! Terminate cycle out condition
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stph_rp2
END MODULE stph_rp2_mod
