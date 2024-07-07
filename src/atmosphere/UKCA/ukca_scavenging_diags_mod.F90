! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module for plume scavenging diagnostics to capture before and after state
!    of UKCA tracers in convection, then calculate the rate of change of each
!    aerosol tracer in units of mol/gridbox/s . Contained routines:
!    ukca_plume_scav_initial   (called from ni_conv_ctl)
!    ukca_plume_scav_diags     (called from ukca_main1)
!
!  Method:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office. See
!  www.ukca.ac.uk
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90 (free-form)
!
! ######################################################################
MODULE ukca_scavenging_diags_mod

USE ukca_mode_setup,  ONLY: nmodes, ncp, ncp_max
USE um_stashcode_mod, ONLY:  stashcode_so4_nuc_sol_ps,   &
  stashcode_oc_nuc_sol_ps,  stashcode_so_nuc_sol_ps,    &
  stashcode_so4_ait_sol_ps, stashcode_bc_ait_sol_ps,    &
  stashcode_oc_ait_sol_ps,  stashcode_so_ait_sol_ps,    &
  stashcode_so4_acc_sol_ps, stashcode_bc_acc_sol_ps,    &
  stashcode_oc_acc_sol_ps,  stashcode_ss_acc_sol_ps,    &
  stashcode_du_acc_sol_ps,  stashcode_so_acc_sol_ps,    &
  stashcode_no3_acc_sol_ps, stashcode_nh4_acc_sol_ps,   &
  stashcode_so4_cor_sol_ps, stashcode_bc_cor_sol_ps,    &
  stashcode_oc_cor_sol_ps,  stashcode_ss_cor_sol_ps,    &
  stashcode_du_cor_sol_ps,  stashcode_so_cor_sol_ps,    &
  stashcode_no3_cor_sol_ps, stashcode_nh4_cor_sol_ps,   &
  stashcode_oc_ait_insol_ps, stashcode_bc_ait_insol_ps, &
  stashcode_du_acc_insol_ps, stashcode_du_cor_insol_ps, &
  stashcode_glomap_sec
IMPLICIT NONE
PRIVATE
SAVE

! Tracer arrays for 2D diagnostic calculations
REAL, ALLOCATABLE :: pre_tracers(:,:,:,:)     ! To hold tracer mixing ratios
                                              ! before scavenging
! To hold tracer mixing ratio differences
REAL, ALLOCATABLE, PUBLIC :: d_tracers(:,:,:,:)

! Rate of change due to plume scavenging (2D) (mol/gridbox/s)
REAL, ALLOCATABLE, PUBLIC :: plume_scav_diag_2d(:,:)

! Table of section 38 stash item numbers for 2D plume scavenging diagnostics
INTEGER, PUBLIC :: item_no_plume(nmodes,ncp_max) =          &
  RESHAPE( (/                                               &
! SO4
  stashcode_so4_nuc_sol_ps,   stashcode_so4_ait_sol_ps,     &
  stashcode_so4_acc_sol_ps,   stashcode_so4_cor_sol_ps,     &
  -1,                         -1,                           &
  -1,                                                       &
! BC
  -1,                         stashcode_bc_ait_sol_ps,      &
  stashcode_bc_acc_sol_ps,    stashcode_bc_cor_sol_ps,      &
  stashcode_bc_ait_insol_ps,  -1,                           &
  -1,                                                       &
! OC
  stashcode_oc_nuc_sol_ps,    stashcode_oc_ait_sol_ps,      &
  stashcode_oc_acc_sol_ps,    stashcode_oc_cor_sol_ps,      &
  stashcode_oc_ait_insol_ps,  -1,                           &
  -1,                                                       &
! NaCl
  -1,                         -1,                           &
  stashcode_ss_acc_sol_ps,    stashcode_ss_cor_sol_ps,      &
  -1,                         -1,                           &
  -1,                                                       &
! Dust
  -1,                         -1,                           &
  stashcode_du_acc_sol_ps,    stashcode_du_cor_sol_ps,      &
  -1,                         stashcode_du_acc_insol_ps,    &
  stashcode_du_cor_insol_ps,                                &
! Secondary Organic
  stashcode_so_nuc_sol_ps,    stashcode_so_ait_sol_ps,      &
  stashcode_so_acc_sol_ps,    stashcode_so_cor_sol_ps,      &
  -1,                         -1,                           &
  -1,                                                       &
! NO3
                        -1,                         -1,     &
  stashcode_no3_acc_sol_ps,   stashcode_no3_cor_sol_ps,     &
  -1,                         -1,                           &
  -1,                                                       &
! NH4
                        -1,                         -1,     &
  stashcode_nh4_acc_sol_ps,   stashcode_nh4_cor_sol_ps,     &
  -1,                         -1,                           &
  -1                                                        &
  /),         (/nmodes, ncp_max/))

INTEGER, PARAMETER, PUBLIC :: msect38 = 1000*stashcode_glomap_sec

LOGICAL, PUBLIC :: plume_scav_diag_ping = .FALSE.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SCAVENGING_DIAGS_MOD'

PUBLIC :: ukca_plume_scav_diags_2d
PUBLIC :: ukca_plume_scav_initial
PUBLIC :: ukca_plume_scav_final

CONTAINS

! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_plume_scav_initial(tr_ukca, tracer_ukca)

!   Captures the tracer state before plume scavenging

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY: tdims_s
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook

IMPLICIT NONE

! Number of ukca tracers
INTEGER, INTENT(IN) :: tr_ukca

! UKCA tracer mixing ratios
REAL,    INTENT(IN) :: tracer_ukca(tdims_s%i_start:tdims_s%i_end,          &
                                   tdims_s%j_start:tdims_s%j_end,          &
                                   tdims_s%k_start:tdims_s%k_end,          &
                                   tr_ukca)

CHARACTER (LEN=* ), PARAMETER :: RoutineName = 'UKCA_PLUME_SCAV_INITIAL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (.NOT. ALLOCATED(pre_tracers))                                          &
               ALLOCATE(pre_tracers(tdims_s%i_start:tdims_s%i_end,         &
                                    tdims_s%j_start:tdims_s%j_end,         &
                                    tdims_s%k_start:tdims_s%k_end,         &
                                    tr_ukca) )

pre_tracers(:,:,:,:) = tracer_ukca(:,:,:,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_plume_scav_initial
! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_plume_scav_final(tr_ukca, tracer_ukca)

!   Captures the tracer state after plume scavenging

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY: tdims_s
USE parkind1,              ONLY: jprb, jpim
USE yomhook,               ONLY: lhook, dr_hook

IMPLICIT NONE

! Number of ukca tracers
INTEGER, INTENT(IN) :: tr_ukca

! UKCA tracer mixing ratios
REAL,    INTENT(IN) :: tracer_ukca(tdims_s%i_start:tdims_s%i_end,          &
                                   tdims_s%j_start:tdims_s%j_end,          &
                                   tdims_s%k_start:tdims_s%k_end,          &
                                   tr_ukca)

CHARACTER (LEN=* ), PARAMETER :: RoutineName = 'UKCA_PLUME_SCAV_FINAL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(d_tracers))                                            &
               ALLOCATE(d_tracers(tdims_s%i_start:tdims_s%i_end,           &
                                  tdims_s%j_start:tdims_s%j_end,           &
                                  tdims_s%k_start:tdims_s%k_end,           &
                                  tr_ukca) )

! Calculate ukca mixing ratio differences from plume scavenging
d_tracers(:,:,:,:) = pre_tracers(:,:,:,:) - tracer_ukca(:,:,:,:)

DEALLOCATE(pre_tracers)

! Set indicator when this routine is called from ni_conv_ctl
plume_scav_diag_ping = .TRUE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_plume_scav_final

! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_plume_scav_diags_2d(imode_req, icpt_req, row_length, rows, &
                                    model_levels, tr_ukca)

!   Convert aerosol tracer mass mixing ratio differences into plume scavenging
!    diagnostic fluxes in units of mol/gridbox/s .

! Definitions of prognostic variable array sizes
USE model_time_mod,         ONLY: secs_per_stepim
USE submodel_mod,           ONLY: atmos_im
USE ukca_d1_defs,           ONLY: n_mode_tracers, n_aero_tracers
USE ukca_um_interf_mod,     ONLY: rho_r2
USE ukca_mode_setup,        ONLY: nmodes, ncp, mm
USE ukca_scavenging_mod, ONLY:  mmr_index_um
USE ukca_eg_tracers_total_mass_mod, ONLY: ukca_eg_tracers_total_mass_fix
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: umPrint, umMessage
USE parkind1,               ONLY: jprb, jpim
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: imode_req   ! required mode
INTEGER, INTENT(IN) :: icpt_req    ! required component
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
! Number of ukca tracers
INTEGER, INTENT(IN) :: tr_ukca

! Local variables

! Vertically integrated tracer mass. Need rank 3 array 
! for call to ukca_eg_tracers_total_mass_fix, so 3rd 
! dimension set to 1
REAL :: mass_2d(row_length,rows,1)

REAL     :: timestep                  ! model timestep

INTEGER  :: j                         ! tracer identifier 
INTEGER  :: errcode                   ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage           ! Error message
CHARACTER (LEN=* ), PARAMETER :: RoutineName = 'UKCA_PLUME_SCAV_DIAGS_2D'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

timestep = secs_per_stepim(atmos_im)    ! timestep in seconds

IF (.NOT. ALLOCATED(plume_scav_diag_2d)) THEN
    ALLOCATE(plume_scav_diag_2d(row_length, rows))
    plume_scav_diag_2d(:,:) = 0.0
END IF

IF (plume_scav_diag_ping) THEN

    ! set the index in d_tracers of the component needed
    j = mmr_index_um(imode_req,icpt_req)

    ! check that we do not have an out of bounds issue
    ! with the index
    IF (j < 1.OR. j > tr_ukca) THEN
        cmessage = ' Index j is out of range'
        errcode = ABS(j)
        WRITE(umMessage,'(A10,2I5)') 'errcode: ',errcode,j
        CALL umPrint(umMessage,src=RoutineName)
        CALL ereport(RoutineName,errcode,cmessage)
    END IF

    IF (ALLOCATED(d_tracers)) THEN

    !   Fill 2-D tracer mass array for this tracer from the array holding
    !   differences before and after convection for all tracers.
    !   Uses ukca_eg_tracers_total_mass_fix to return vertically integrated
    !   tracer mass.
        CALL ukca_eg_tracers_total_mass_fix(                                &
            rho_r2(1:row_length,1:rows,1:model_levels),                   &
            d_tracers(1:row_length,1:rows,0:model_levels,j),              &
            1, mass_2d)
    END IF

    ! Convert mass mixing ratio differences in kg/kg per timestep 
    ! to mol/s/gridbox
    plume_scav_diag_2d(:,:) = mass_2d(:,:,1) / (timestep * mm(icpt_req))

ELSE
    plume_scav_diag_2d(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE ukca_plume_scav_diags_2d

END MODULE ukca_scavenging_diags_mod
