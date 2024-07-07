! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Purpose:
!  Gather D1 information about UKCA-MODE inputs required by UKCA_RADAER
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Description:
!  To define and allocate the structure UKCA_RADAER.
!  Note that RADAER assumes that nitrate is in the form NH4NO3, and that
!  sulphate is (NH4)2SO4, so that explicit NH4 information is not required,
!  although space is available for the NH4 component.
!
MODULE ukca_radaer_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_INIT_MOD'

CONTAINS

SUBROUTINE ukca_radaer_init( ierr, cmessage, l_ukca_radaer_sustrat,        &
                             ukca_radaer )

USE parkind1,               ONLY: jpim, jprb
USE yomhook,                ONLY: lhook, dr_hook
USE ukca_mode_setup,        ONLY: nmodes, ncp, mode, component,   &
                                  mode_names, modesol, ddplim0,   &
                                  ddplim1, sigmag, rhocomp,       &
                                  component_names, ncp_max

USE ukca_radaer_struct_mod,  ONLY: &
     ip_ukca_mode_nucleation,      &
     ip_ukca_mode_aitken,          &
     ip_ukca_mode_accum,           &
     ip_ukca_mode_coarse,          &
     ip_ukca_sulphate,             &
     ip_ukca_blackcarbon,          &
     ip_ukca_organiccarbon,        &
     ip_ukca_seasalt,              &
     ip_ukca_dust,                 &
     ip_ukca_secondorganic,        &
     ukca_radaer_struct

USE um_stashcode_mod,        ONLY: &
    stashcode_ukca_sec,            &
    stashcode_nucsol_no,           &
    stashcode_Aitsol_no,           &
    stashcode_accsol_no,           &
    stashcode_corsol_no,           &
    stashcode_Aitinsol_no,         &
    stashcode_accinsol_no,         &
    stashcode_corinsol_no,         &
    stashcode_nucsol_so4,          &
    stashcode_Aitsol_so4,          &
    stashcode_accsol_so4,          &
    stashcode_corsol_so4,          &
    stashcode_Aitsol_bc,           &
    stashcode_accsol_bc,           &
    stashcode_corsol_bc,           &
    stashcode_Aitinsol_bc,         &
    stashcode_nucsol_oc,           &
    stashcode_Aitsol_oc,           &
    stashcode_accsol_oc,           &
    stashcode_corsol_oc,           &
    stashcode_Aitinsol_oc,         &
    stashcode_accsol_ss,           &
    stashcode_corsol_ss,           &
    stashcode_accsol_du,           &
    stashcode_corsol_du,           &
    stashcode_accinsol_du,         &
    stashcode_corinsol_du,         &
    stashcode_nucsol_so,           &
    stashcode_Aitsol_so,           &
    stashcode_accsol_so,           &
    stashcode_corsol_so,           &
    stashcode_dryd_ait_sol,        &
    stashcode_dryd_acc_sol,        &
    stashcode_dryd_cor_sol,        &
    stashcode_dryd_ait_insol,      &
    stashcode_dryd_acc_insol,      &
    stashcode_dryd_cor_insol,      &
    stashcode_wetd_ait_sol,        &
    stashcode_wetd_acc_sol,        &
    stashcode_wetd_cor_sol,        &
    stashcode_rho_ait_sol,         &
    stashcode_rho_acc_sol,         &
    stashcode_rho_cor_sol,         &
    stashcode_rho_ait_insol,       &
    stashcode_rho_acc_insol,       &
    stashcode_rho_cor_insol,       &
    stashcode_pvol_ait_su_sol,     &
    stashcode_pvol_ait_bc_sol,     &
    stashcode_pvol_ait_oc_sol,     &
    stashcode_pvol_ait_so_sol,     &
    stashcode_pvol_ait_h2o_sol,    &
    stashcode_pvol_acc_su_sol,     &
    stashcode_pvol_acc_bc_sol,     &
    stashcode_pvol_acc_oc_sol,     &
    stashcode_pvol_acc_so_sol,     &
    stashcode_pvol_acc_du_sol,     &
    stashcode_pvol_acc_ss_sol,     &
    stashcode_pvol_acc_no3_sol,    &
    stashcode_pvol_acc_h2o_sol,    &
    stashcode_pvol_cor_su_sol,     &
    stashcode_pvol_cor_bc_sol,     &
    stashcode_pvol_cor_oc_sol,     &
    stashcode_pvol_cor_so_sol,     &
    stashcode_pvol_cor_du_sol,     &
    stashcode_pvol_cor_ss_sol,     &
    stashcode_pvol_cor_no3_sol,    &
    stashcode_pvol_cor_h2o_sol,    &
    stashcode_pvol_ait_bc_insol,   &
    stashcode_pvol_ait_oc_insol,   &
    stashcode_pvol_acc_du_insol,   &
    stashcode_pvol_cor_du_insol,   & 
    stashcode_glomap_clim_sec,     &
    stashcode_gc_nd_nuc_sol,       &
    stashcode_gc_nd_ait_sol,       &
    stashcode_gc_nd_acc_sol,       &
    stashcode_gc_nd_cor_sol,       &
    stashcode_gc_nd_ait_ins,       &
    stashcode_gc_nd_acc_ins,       &
    stashcode_gc_nd_cor_ins,       &
    stashcode_gc_nuc_sol_su,       &
    stashcode_gc_ait_sol_su,       &
    stashcode_gc_acc_sol_su,       &
    stashcode_gc_cor_sol_su,       &
    stashcode_gc_ait_sol_bc,       &
    stashcode_gc_acc_sol_bc,       &
    stashcode_gc_cor_sol_bc,       &
    stashcode_gc_ait_ins_bc,       &
    stashcode_gc_nuc_sol_oc,       &
    stashcode_gc_ait_sol_oc,       &
    stashcode_gc_acc_sol_oc,       &
    stashcode_gc_cor_sol_oc,       &
    stashcode_gc_ait_ins_oc,       &
    stashcode_gc_ait_sol_ss,       &
    stashcode_gc_acc_sol_ss,       &
    stashcode_gc_cor_sol_ss,       &
    stashcode_gc_acc_sol_du,       &
    stashcode_gc_cor_sol_du,       &
    stashcode_gc_acc_ins_du,       &
    stashcode_gc_cor_ins_du,       &
    stashcode_gc_nuc_sol_so,       &
    stashcode_gc_ait_sol_so,       &
    stashcode_gc_acc_sol_so,       &
    stashcode_gc_cor_sol_so,       &
    stashcode_gc_acc_sol_no3,      &
    stashcode_gc_cor_sol_no3,      &
    stashcode_gc_dryd_ait_sol,     &
    stashcode_gc_dryd_acc_sol,     &
    stashcode_gc_dryd_cor_sol,     &
    stashcode_gc_dryd_ait_ins,     &
    stashcode_gc_dryd_acc_ins,     &
    stashcode_gc_dryd_cor_ins,     &
    stashcode_gc_wetd_ait_sol,     &
    stashcode_gc_wetd_acc_sol,     &
    stashcode_gc_wetd_cor_sol,     &
    stashcode_gc_rho_ait_sol,      &
    stashcode_gc_rho_acc_sol,      &
    stashcode_gc_rho_cor_sol,      &
    stashcode_gc_rho_ait_ins,      &
    stashcode_gc_rho_acc_ins,      &
    stashcode_gc_rho_cor_ins,      &
    stashcode_gc_pvol_ait_su_sol,  &
    stashcode_gc_pvol_ait_bc_sol,  &
    stashcode_gc_pvol_ait_oc_sol,  &
    stashcode_gc_pvol_ait_so_sol,  &
    stashcode_gc_pvol_ait_h2o_sol, &
    stashcode_gc_pvol_acc_su_sol,  &
    stashcode_gc_pvol_acc_bc_sol,  &
    stashcode_gc_pvol_acc_oc_sol,  &
    stashcode_gc_pvol_acc_ss_sol,  &
    stashcode_gc_pvol_acc_du_sol,  &
    stashcode_gc_pvol_acc_so_sol,  &
    stashcode_gc_pvol_acc_no3_sol, &
    stashcode_gc_pvol_acc_h2o_sol, &
    stashcode_gc_pvol_cor_su_sol,  &
    stashcode_gc_pvol_cor_bc_sol,  &
    stashcode_gc_pvol_cor_oc_sol,  &
    stashcode_gc_pvol_cor_ss_sol,  &
    stashcode_gc_pvol_cor_du_sol,  &
    stashcode_gc_pvol_cor_so_sol,  &
    stashcode_gc_pvol_cor_no3_sol, &
    stashcode_gc_pvol_cor_h2o_sol, &
    stashcode_gc_pvol_ait_bc_ins,  &
    stashcode_gc_pvol_ait_oc_ins,  &
    stashcode_gc_pvol_acc_du_ins,  &
    stashcode_gc_pvol_cor_du_ins,  &
    stashcode_pvol_acc_no3_sol,    &
    stashcode_pvol_cor_no3_sol

USE errormessagelength_mod, ONLY: errormessagelength

USE glomap_clim_option_mod,  ONLY: &
    l_glomap_clim_radaer

USE ukca_option_mod,         ONLY: &
    l_ukca_radaer

USE ukca_radaer_struct_mod,  ONLY: &
    allocate_radaer_struct

IMPLICIT NONE

!
! Arguments
!
!
! Error indicator (0 is OK, >0 error)
!
INTEGER, INTENT(INOUT) :: ierr

!
! Error message if ierr is larger than 0
!
CHARACTER (LEN=errormessagelength), INTENT(INOUT) :: cmessage

!
! Switch: Use sulphuric acid optical properties in the stratosphere,
! instead of ammonium sulphate.
!
LOGICAL, INTENT(IN) :: l_ukca_radaer_sustrat

!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer
!
! Local variables
!
! Loop variables
!
INTEGER :: i, j

!
! In-loop copy of mode names
!
CHARACTER(LEN=7) :: this_name

!
! In-loop mode type
!
INTEGER :: this_type

!
! Local counters for number of modes and components
!
INTEGER :: n_loc_mode
INTEGER :: n_loc_cpnt

INTEGER, PARAMETER :: msect_gc = 1000 * stashcode_glomap_clim_sec
INTEGER, PARAMETER :: msect    = 1000 * stashcode_ukca_sec



!
! STASH codes for the diagnostics given the modal dry and wet
! diameters, the component mass-mixing ratios, and the modal
! densities.
! Those arrays should be obtained through UKCA_MODE_SETUP for
! maximum flexibility.


!
! The following prognostics are expected in section 34.
!
INTEGER, PARAMETER :: stashc_nbr(nmodes) = (/                                  &
  stashcode_nucsol_no, stashcode_Aitsol_no,   stashcode_accsol_no,             &
  stashcode_corsol_no, stashcode_Aitinsol_no, stashcode_accinsol_no,           &
  stashcode_corinsol_no /)

! Mass mixing ratios
INTEGER, PARAMETER :: stashc_mmr(nmodes,ncp_max) = RESHAPE( SOURCE = (/        &
    ! nuc_sol                ait_sol                acc_sol
    ! cor_sol                ait_insol              acc_insol
    ! cor_insol
! SO4
      stashcode_nucsol_so4,  stashcode_Aitsol_so4,  stashcode_accsol_so4,      &
      stashcode_corsol_so4,                    -1,                    -1,      &
                        -1,                                                    &
! BC
                        -1,   stashcode_Aitsol_bc,   stashcode_accsol_bc,      &
       stashcode_corsol_bc, stashcode_Aitinsol_bc,                    -1,      &
                        -1,                                                    &
! OC
       stashcode_nucsol_oc,   stashcode_Aitsol_oc,   stashcode_accsol_oc,      &
       stashcode_corsol_oc, stashcode_Aitinsol_oc,                    -1,      &
                        -1,                                                    &
! NaCl
                        -1,                    -1,   stashcode_accsol_ss,      &
       stashcode_corsol_ss,                    -1,                    -1,      &
                        -1,                                                    &
! Dust
                        -1,                    -1,   stashcode_accsol_du,      &
       stashcode_corsol_du,                    -1, stashcode_accinsol_du,      &
     stashcode_corinsol_du,                                                    &
! Sec_Org
       stashcode_nucsol_so,   stashcode_Aitsol_so,   stashcode_accsol_so,      &
       stashcode_corsol_so,                    -1,                    -1,      &
                        -1,                                                    &
! NO3
                        -1,                    -1,                    -1,      &
                        -1,                    -1,                    -1,      &
                        -1,                                                    &
! NH4  
                        -1,                    -1,                    -1,      &
                        -1,                    -1,                    -1,      &
                        -1                                                     &
/), SHAPE = (/ nmodes, ncp_max /) )

! Component volumes (including water as a component)
INTEGER, PARAMETER :: stashc_cvl(nmodes,ncp_max+1) = RESHAPE( SOURCE = (/   &
    ! nuc_sol                ait_sol                acc_sol
    ! cor_sol                ait_insol
    ! acc_insol              cor_insol
! SO4
    -1, stashcode_pvol_ait_su_sol, stashcode_pvol_acc_su_sol,               &
        stashcode_pvol_cor_su_sol, -1,                                      &
        -1,                        -1,                                      &
! BC
    -1, stashcode_pvol_ait_bc_sol,   stashcode_pvol_acc_bc_sol,             &
        stashcode_pvol_cor_bc_sol,   stashcode_pvol_ait_bc_insol,           &
        -1,                         -1,                                     &
! OC
    -1, stashcode_pvol_ait_oc_sol,   stashcode_pvol_acc_oc_sol,             &
        stashcode_pvol_cor_oc_sol,   stashcode_pvol_ait_oc_insol,           &
        -1,                          -1,                                    &
! NaCl
    -1, -1,                          stashcode_pvol_acc_ss_sol,             &
        stashcode_pvol_cor_ss_sol,   -1,                                    &
        -1,                          -1,                                    &
! Dust
    -1, -1,                          stashcode_pvol_acc_du_sol,             &
        stashcode_pvol_cor_du_sol,   -1,                                    &
        stashcode_pvol_acc_du_insol, stashcode_pvol_cor_du_insol,           &
! Sec_Org
    -1, stashcode_pvol_ait_so_sol,   stashcode_pvol_acc_so_sol,             &
        stashcode_pvol_cor_so_sol,   -1,                                    &
        -1,                          -1,                                    &
! NO3     
    -1,                         -1,  stashcode_pvol_acc_no3_sol,            &
        stashcode_pvol_cor_no3_sol,  -1,                                    &
        -1,                          -1,                                    &
! NH4
    -1, -1,                          -1,                                    &
        -1,                          -1,                                    &
        -1,                          -1,                                    &
! Water
    -1, stashcode_pvol_ait_h2o_sol,  stashcode_pvol_acc_h2o_sol,            &
        stashcode_pvol_cor_h2o_sol,  -1,                                    &
        -1,                          -1                                     &
    /), SHAPE = (/ nmodes, ncp_max+1 /) )

! Number mass mixing ratios
! 7 modes : nuc_sol, ait_sol, acc_sol, cor_sol, ait_ins, acc_ins, cor_ins
INTEGER, PARAMETER :: stashc_gc_nbr(nmodes) = (/ stashcode_gc_nd_nuc_sol,      &
                                                 stashcode_gc_nd_ait_sol,      &
                                                 stashcode_gc_nd_acc_sol,      &
                                                 stashcode_gc_nd_cor_sol,      &
                                                 stashcode_gc_nd_ait_ins,      &
                                                 stashcode_gc_nd_acc_ins,      &
                                                 stashcode_gc_nd_cor_ins  /)

! Component mass mixing ratios
! 7 modes : nuc_sol, ait_sol, acc_sol, cor_sol, ait_ins, acc_ins, cor_ins
! 8 components : sulphate, black_carbon, organic_carbon, salt, dust, 2nd_organic
!                nitrate,  ammonium
INTEGER, PARAMETER :: stashc_gc_mmr(nmodes, ncp_max) =                         &
                           RESHAPE( (/ stashcode_gc_nuc_sol_su,                &
                                       stashcode_gc_ait_sol_su,                &
                                       stashcode_gc_acc_sol_su,                &
                                       stashcode_gc_cor_sol_su,                &
                                       msect_gc - 1,    &      ! su ait_ins
                                       msect_gc - 1,    &      ! su acc_ins
                                       msect_gc - 1,    &      ! su cor_ins
                                       msect_gc - 1,    &      ! bc nuc_sol
                                       stashcode_gc_ait_sol_bc,                &
                                       stashcode_gc_acc_sol_bc,                &
                                       stashcode_gc_cor_sol_bc,                &
                                       stashcode_gc_ait_ins_bc,                &
                                       msect_gc - 1,    &      ! bc acc_ins
                                       msect_gc - 1,    &      ! bc cor_ins
                                       stashcode_gc_nuc_sol_oc,                &
                                       stashcode_gc_ait_sol_oc,                &
                                       stashcode_gc_acc_sol_oc,                &
                                       stashcode_gc_cor_sol_oc,                &
                                       stashcode_gc_ait_ins_oc,                &
                                       msect_gc - 1,    &      ! oc acc_ins
                                       msect_gc - 1,    &      ! oc cor_ins
                                       msect_gc - 1,    &      ! ss nuc_sol
                                       stashcode_gc_ait_sol_ss,                &
                                       stashcode_gc_acc_sol_ss,                &
                                       stashcode_gc_cor_sol_ss,                &
                                       msect_gc - 1,    &      ! ss ait_ins
                                       msect_gc - 1,    &      ! ss acc_ins
                                       msect_gc - 1,    &      ! ss cor_ins
                                       msect_gc - 1,    &      ! du nuc_sol
                                       msect_gc - 1,    &      ! du ait_sol
                                       stashcode_gc_acc_sol_du,                &
                                       stashcode_gc_cor_sol_du,                &
                                       msect_gc - 1,    &      ! du ait_ins
                                       stashcode_gc_acc_ins_du,                &
                                       stashcode_gc_cor_ins_du,                &
                                       stashcode_gc_nuc_sol_so,                &
                                       stashcode_gc_ait_sol_so,                &
                                       stashcode_gc_acc_sol_so,                &
                                       stashcode_gc_cor_sol_so,                &
                                       msect_gc - 1,    &      ! so ait_ins
                                       msect_gc - 1,    &      ! so acc_ins
                                       msect_gc - 1,    &      ! so cor_ins
                                       msect_gc - 1,    &      ! no3 nuc_sol
                                       msect_gc - 1,    &      ! no3 ait_sol
                                       stashcode_gc_acc_sol_no3,               &
                                       stashcode_gc_cor_sol_no3,               &
                                       msect_gc - 1,    &      ! no3 ait_ins
                                       msect_gc - 1,    &      ! no3 acc_ins
                                       msect_gc - 1,    &      ! no3 cor_ins
                                       msect_gc - 1,    &      ! nh4 nuc_sol
                                       msect_gc - 1,    &      ! nh4 ait_sol
                                       msect_gc - 1,    &      ! nh4_acc_sol 
                                       msect_gc - 1,    &      ! nh4_cor_sol
                                       msect_gc - 1,    &      ! nh4 ait_ins
                                       msect_gc - 1,    &      ! nh4 acc_ins
                                       msect_gc - 1 /)  &      ! nh4 cor_ins
                                       , (/ nmodes, ncp_max /) )


!
! The following diagnostics are expected in section 34.
!
INTEGER, PARAMETER :: stashc_dry(nmodes) = (/ -1,                    &
           stashcode_dryd_ait_sol,   stashcode_dryd_acc_sol,         &
           stashcode_dryd_cor_sol,   stashcode_dryd_ait_insol,       &
           stashcode_dryd_acc_insol, stashcode_dryd_cor_insol        &
           /)
INTEGER, PARAMETER :: stashc_wet(nmodes) = (/ -1,                    &
           stashcode_wetd_ait_sol,   stashcode_wetd_acc_sol,         &
           stashcode_wetd_cor_sol,     -1,     -1,     -1            &
           /)
INTEGER, PARAMETER :: stashc_rho(nmodes) = (/ -1,                    &
           stashcode_rho_ait_sol,   stashcode_rho_acc_sol,           &
           stashcode_rho_cor_sol,   stashcode_rho_ait_insol,         &
           stashcode_rho_acc_insol, stashcode_rho_cor_insol          &
           /)

!
! Diagnostics that will have to be taken in also include:
!
!
!
! The following diagnostics are expected in Section 54
!

! Dry diameter stashcodes
INTEGER, PARAMETER :: stashc_gc_dry(nmodes) = (/ -1,                & ! nuc_sol
                                                 stashcode_gc_dryd_ait_sol,    &
                                                 stashcode_gc_dryd_acc_sol,    &
                                                 stashcode_gc_dryd_cor_sol,    &
                                                 stashcode_gc_dryd_ait_ins,    &
                                                 stashcode_gc_dryd_acc_ins,    &
                                                 stashcode_gc_dryd_cor_ins /)

! Wet diameter stashcodes
INTEGER, PARAMETER :: stashc_gc_wet(nmodes) = (/ -1,                & ! nuc_sol
                                                 stashcode_gc_wetd_ait_sol,    &
                                                 stashcode_gc_wetd_acc_sol,    &
                                                 stashcode_gc_wetd_cor_sol,    &
                                                 -1, -1, -1 /)

! Density stashcodes
INTEGER, PARAMETER :: stashc_gc_rho(nmodes) = (/ -1,                & ! nuc_sol
                                                 stashcode_gc_rho_ait_sol,     &
                                                 stashcode_gc_rho_acc_sol,     &
                                                 stashcode_gc_rho_cor_sol,     &
                                                 stashcode_gc_rho_ait_ins,     &
                                                 stashcode_gc_rho_acc_ins,     &
                                                 stashcode_gc_rho_cor_ins /)

! Component volumes (including water as a component)
! 7 modes : nuc_sol, ait_sol, acc_sol, cor_sol, ait_ins, acc_ins, cor_ins
! 7 components : su, bc, oc, ss, du, so, h20
INTEGER, PARAMETER :: stashc_gc_cvl(nmodes, ncp_max+1) =                       &
                           RESHAPE( (/          msect_gc - 1,   & ! su  nuc_sol
                                                stashcode_gc_pvol_ait_su_sol,  &
                                                stashcode_gc_pvol_acc_su_sol,  &
                                                stashcode_gc_pvol_cor_su_sol,  &
                                                msect_gc - 1,   & ! su  ait_ins
                                                msect_gc - 1,   & ! su  acc_ins
                                                msect_gc - 1,   & ! su  cor_ins
                                                msect_gc - 1,   & ! bc  nuc_sol
                                                stashcode_gc_pvol_ait_bc_sol,  &
                                                stashcode_gc_pvol_acc_bc_sol,  &
                                                stashcode_gc_pvol_cor_bc_sol,  &
                                                stashcode_gc_pvol_ait_bc_ins,  &
                                                msect_gc - 1,   & ! bc  acc_ins
                                                msect_gc - 1,   & ! bc  cor_ins
                                                msect_gc - 1,   & ! oc  nuc_sol
                                                stashcode_gc_pvol_ait_oc_sol,  &
                                                stashcode_gc_pvol_acc_oc_sol,  &
                                                stashcode_gc_pvol_cor_oc_sol,  &
                                                stashcode_gc_pvol_ait_oc_ins,  &
                                                msect_gc - 1,   & ! oc  acc_ins
                                                msect_gc - 1,   & ! oc  cor_ins
                                                msect_gc - 1,   & ! ss  nuc_sol
                                                msect_gc - 1,   & ! ss  ait_sol
                                                stashcode_gc_pvol_acc_ss_sol,  &
                                                stashcode_gc_pvol_cor_ss_sol,  &
                                                msect_gc - 1,   & ! ss  ait_ins
                                                msect_gc - 1,   & ! ss  acc_ins
                                                msect_gc - 1,   & ! ss  cor_ins
                                                msect_gc - 1,   & ! du  nuc_sol
                                                msect_gc - 1,   & ! du  ait_sol
                                                stashcode_gc_pvol_acc_du_sol,  &
                                                stashcode_gc_pvol_cor_du_sol,  &
                                                msect_gc - 1,   & ! du  ait_ins
                                                stashcode_gc_pvol_acc_du_ins,  &
                                                stashcode_gc_pvol_cor_du_ins,  &
                                                msect_gc - 1,   & ! so  nuc_sol
                                                stashcode_gc_pvol_ait_so_sol,  &
                                                stashcode_gc_pvol_acc_so_sol,  &
                                                stashcode_gc_pvol_cor_so_sol,  &
                                                msect_gc - 1,   & ! so  ait_ins
                                                msect_gc - 1,   & ! so  acc_ins
                                                msect_gc - 1,   & ! so  cor_ins
                                                msect_gc - 1,   & ! no3 nuc_sol
                                                msect_gc - 1,   & ! no3 ait_sol
                                                stashcode_gc_pvol_acc_no3_sol, &
                                                stashcode_gc_pvol_cor_no3_sol, &
                                                msect_gc - 1,   & ! no3 ait_ins
                                                msect_gc - 1,   & ! no3 acc_ins
                                                msect_gc - 1,   & ! no3 cor_ins
                                                msect_gc - 1,   & ! nh4 nuc_sol
                                                msect_gc - 1,   & ! nh4 ait_sol
                                                msect_gc - 1,   & ! nh4 acc_sol
                                                msect_gc - 1,   & ! nh4 cor_sol
                                                msect_gc - 1,   & ! nh4 ait_ins
                                                msect_gc - 1,   & ! nh4 acc_ins
                                                msect_gc - 1,   & ! nh4 cor_ins
                                                msect_gc - 1,   & ! h2o nuc_sol
                                                stashcode_gc_pvol_ait_h2o_sol, &
                                                stashcode_gc_pvol_acc_h2o_sol, &
                                                stashcode_gc_pvol_cor_h2o_sol, &
                                                msect_gc - 1,   & ! h2o ait_ins
                                                msect_gc - 1,   & ! h2o acc_ins
                                                msect_gc - 1 /) & ! h2o cor_ins
                                                , (/ nmodes, ncp_max+1 /) )


! Index of the "water" component in array stashc_cvl above.
INTEGER, PARAMETER :: ip_water_index =  ncp_max + 1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

! Allocate elements that depend on ncp or nmodes
CALL allocate_radaer_struct(ukca_radaer)

! Loop on all modes and components per mode, and only retain
! included modes and components.
!
n_loc_mode = 0
n_loc_cpnt = 0
DO i = 1, nmodes

  IF (mode(i)) THEN

    this_name = mode_names(i)
    !
    ! Get the mode type. Since there is no direct information,
    ! it is obtained from the mode names.
    !
    SELECT CASE (this_name(1:3))

    CASE ('Nuc')
      this_type = ip_ukca_mode_nucleation

    CASE ('Ait')
      this_type = ip_ukca_mode_aitken

    CASE ('Acc')
      this_type = ip_ukca_mode_accum

    CASE ('Cor')
      this_type = ip_ukca_mode_coarse

    CASE DEFAULT
      ierr = 701
      cmessage = 'ukca_radaer_init: Unexpected mode name.' // this_name
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,         &
                                                            zhook_handle)
      RETURN

    END SELECT

    !
    ! Interaction of nucleation modes with radiation is
    ! neglected.
    !
    IF (this_type /= ip_ukca_mode_nucleation) THEN

      n_loc_mode = n_loc_mode + 1

      ukca_radaer%i_mode_type(n_loc_mode) = this_type

      ukca_radaer%l_soluble(n_loc_mode) = modesol(i) == 1

      ukca_radaer%d0low(n_loc_mode) = ddplim0(i)

      ukca_radaer%d0up(n_loc_mode) = ddplim1(i)

      ukca_radaer%sigma(n_loc_mode) = sigmag(i)

      ukca_radaer%n_cpnt_in_mode(n_loc_mode) = 0
      
      IF (l_glomap_clim_radaer) THEN
        ukca_radaer%stashcode_dry(n_loc_mode) = stashc_gc_dry(i) - msect_gc
        ukca_radaer%stashcode_wet(n_loc_mode) = stashc_gc_wet(i) - msect_gc
        ukca_radaer%stashcode_rho(n_loc_mode) = stashc_gc_rho(i) - msect_gc
        ukca_radaer%stashcode_wtv(n_loc_mode) = stashc_gc_cvl(i,               &
                                                 ip_water_index) - msect_gc
        ukca_radaer%stashcode_nbr(n_loc_mode) = stashc_gc_nbr(i) - msect_gc
      ELSE IF (l_ukca_radaer) THEN
        ukca_radaer%stashcode_dry(n_loc_mode) = stashc_dry(i) - msect
        ukca_radaer%stashcode_wet(n_loc_mode) = stashc_wet(i) - msect
        ukca_radaer%stashcode_rho(n_loc_mode) = stashc_rho(i) - msect
        ukca_radaer%stashcode_wtv(n_loc_mode) = stashc_cvl(i, ip_water_index)  &
                                                              - msect
        ukca_radaer%stashcode_nbr(n_loc_mode) = stashc_nbr(i) - msect
      ELSE
        ierr = 703
        cmessage = &
             'ukca_radaer_init: Neither l_glomap_clim_radaer nor l_ukca_radaer'
        IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,       &
                                                              zhook_handle)
        RETURN
      
      END IF
      !
      ! Loop on components within that mode.
      !
      DO j = 1, ncp

        IF (component(i, j)) THEN

          n_loc_cpnt = n_loc_cpnt + 1

          !
          ! Update the number of components in that mode and
          ! retain the array index of the current component.
          !
          ukca_radaer%n_cpnt_in_mode(n_loc_mode) =                             &
                                     ukca_radaer%n_cpnt_in_mode(n_loc_mode) + 1

          ukca_radaer%i_cpnt_index(ukca_radaer%n_cpnt_in_mode(n_loc_mode),     &
                                   n_loc_mode) = n_loc_cpnt

          !
          ! Get the component type. Since there is no direct
          ! information, it is obtained from the component
          ! names.
          !
          SELECT CASE (component_names(j))

          CASE ('h2so4  ')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_sulphate

          CASE ('bcarbon')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_blackcarbon

          CASE ('ocarbon')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_organiccarbon

          CASE ('nacl   ')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_seasalt

          CASE ('dust   ')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_dust

          CASE ('sec_org')
            ukca_radaer%i_cpnt_type(n_loc_cpnt) = ip_ukca_secondorganic

          CASE DEFAULT
            ierr = 702
            cmessage = 'ukca_radaer_init: Unexpected component name: ' //      &
                                                             component_names(j)
            IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,   &
                                                                  zhook_handle)
            RETURN

          END SELECT


          ukca_radaer%density(n_loc_cpnt) = rhocomp(j)

          ukca_radaer%i_mode(n_loc_cpnt) = n_loc_mode

          IF (l_glomap_clim_radaer) THEN
            
            ukca_radaer%stashcode_mmr(n_loc_cpnt) = stashc_gc_mmr(i,j)-msect_gc
            
            ukca_radaer%stashcode_cvl(n_loc_cpnt) = stashc_gc_cvl(i,j)-msect_gc

          ELSE IF (l_ukca_radaer) THEN
            
            ukca_radaer%stashcode_mmr(n_loc_cpnt) = stashc_mmr(i, j) - msect
            
            ukca_radaer%stashcode_cvl(n_loc_cpnt) = stashc_cvl(i, j) - msect
            
          ELSE
            ierr = 704
            cmessage = &
             'ukca_radaer_init: Neither l_glomap_clim_radaer nor l_ukca_radaer'
            IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,   &
                                                                  zhook_handle)
            RETURN
          END IF

        END IF

      END DO ! j

    END IF

  END IF

END DO ! i

ukca_radaer%n_mode = n_loc_mode
ukca_radaer%n_cpnt = n_loc_cpnt
ukca_radaer%l_sustrat = l_ukca_radaer_sustrat

IF (ukca_radaer%n_mode == 0 .OR. ukca_radaer%n_cpnt == 0) THEN
  ierr = 703
  cmessage = 'ukca_radaer_init: Setup includes no UKCA aerosols.'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_init
END MODULE ukca_radaer_init_mod
