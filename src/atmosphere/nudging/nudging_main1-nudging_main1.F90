! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Description:
!  The code `nudges' model variables towards meteorological analyses
!  in order to reproduce `real weather'.
!  Top level routine for Nudging. Calls routines to read in ECMWF data
!  fields and nudge model variables (T,U,V) towards ECMWF.
!  Diagnostics related to nudging are written into STASHwork to be updated
!  to STASH by ATM_STEP

!  The nudged model was developed as part of the UKCA project.
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

!  Available for non-commerical use at present. Please cite
!
!  Telford, P. J., Braesicke, P., Morgenstern, O., and Pyle, J. A.:
!  Technical Note: Description and assessment of a nudged version of
!  the new dynamics Unified Model, Atmos. Chem. Phys., 8, 1701-1712, 2008.

!  Method:
! 1) Interface with atmosphere model:
!    Prognostic fields are passed in from Atm_Step. The nudged variables
!    are passed back to Atm_Step as arguments.
! 2) Sets up various quantities (eg filenames of analyses )
!    before calling routines to nudge T/theta, U & V
! 3) The diagnostics from the nudging are put into STASH using
!    copy_diag() and stashwork array.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! CONTAINED subroutines:
!     NUDGING_NETCDF_LOADER
!     NUDGING_ALLOC_DATA
!     NUDGING_GETD1_DATA
!     NUDGING_UPDATE_STASH
!     NUDGING_DEALLOC_DATA

!  Code Description:
!    Language:  FORTRAN 90 (formatted)

! ######################################################################

MODULE nudging_main1_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'NUDGING_MAIN1_MOD'
CONTAINS

! Subroutine Interface:
SUBROUTINE nudging_main1(                                        &
! in time stepping information.
  i_year, i_month, i_day, i_hour, i_minute,                      &
  i_second, timestep_number,                                     &

! in data fields.
  theta, u, v, p_rho_levels, exner_theta_levels, p_theta_levels, &

  st_work)

USE timestep_mod, ONLY: timestep
USE nudging_control              ! Include control module for nudging
USE nudging_d1_defs              ! Use own way of accessing D1 array
USE nudging_filename_mod
USE atm_fields_bounds_mod, ONLY:  udims, vdims, tdims, pdims,&
      udims_s, vdims_s, tdims_s, pdims_s

USE trignometric_mod             ! Include latitude and longitude information
USE cderived_mod, ONLY: delta_lambda, delta_phi, base_phi, base_lambda
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE Field_Types
USE UM_ParVars
USE UM_ParParams, ONLY: pnorth, psouth
USE control_max_sizes
USE submodel_mod, ONLY: submodel_for_sm, atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE model_domain_mod

USE p_to_u_mod
USE p_to_v_mod
USE eg_v_at_poles_mod

USE nlsizes_namelist_mod, ONLY: &
    len_tot, model_levels, n_obj_d1_max, n_rows, row_length, rows
USE d1_array_mod, ONLY: d1_section, d1_item, d1_address, d1_length,            &
                        d1_grid_type, d1, d1_addr, no_obj_d1
USE mpp_conf_mod, ONLY: swap_field_is_scalar

USE nudging_getfrac_mod, ONLY: nudging_getfrac
USE nudging_nudging_control_mod, ONLY: nudging_nudging_control

! Get tendencies for nudging
USE physics_tendencies_mod,  ONLY:                                &
    init_nud_tendencies, dtheta_nud, du_nud, dv_nud,              &
    l_retain_nud_tendencies

IMPLICIT NONE

! ************************************************************************
! arguments with intent in. ie: input variables.

! Parameters

! Submodel indices
INTEGER :: m_atm_modl, im_index

! model parameters
! time information for current timestep
INTEGER ::                                                        &
  i_year, i_month, i_day                                          &
, i_hour, i_minute, i_second                                      &
, timestep_number

! Diagnostics info
REAL :: st_work(*)             ! STASH workspace for section

! Data arrays
REAL, INTENT(INOUT) ::                                 &
!
 u( udims_s%i_start : udims_s%i_end,                   &
    udims_s%j_start : udims_s%j_end,                   &
    udims_s%k_start : udims_s%k_end ),                 &
!
 v( vdims_s%i_start : vdims_s%i_end,                   &
   vdims_s%j_start : vdims_s%j_end,                    &
   vdims_s%k_start : vdims_s%k_end ),                  &
!
 theta( tdims_s%i_start : tdims_s%i_end,               &
        tdims_s%j_start : tdims_s%j_end,               &
        1:model_levels ),                              &
!
 p_rho_levels( pdims_s%i_start : pdims_s%i_end,        &
               pdims_s%j_start : pdims_s%j_end,        &
               pdims_s%k_start : pdims_s%k_end + 1 )

REAL, INTENT(IN) ::                                    &
!
 p_theta_levels( tdims_s%i_start : tdims_s%i_end,      &
                 tdims_s%j_start : tdims_s%j_end,      &
                 tdims_s%k_start : tdims_s%k_end ),    &
!
 exner_theta_levels( tdims_s%i_start : tdims_s%i_end,  &
                     tdims_s%j_start : tdims_s%j_end,  &
                     tdims_s%k_start : tdims_s%k_end )

! *************************************************************************
! Begin local Header

! Miscellaneous variables
LOGICAL, SAVE  :: l_first = .TRUE.        ! Is this the first time
INTEGER        :: a_steps_per_hr          ! no. of atm. steps per hour
INTEGER        :: i,j,k,ii,jj,kk,jnext     ! Loop variables

! **********************************************************************************
! Data_type information from D1
INTEGER                   :: grid_type_theta,field_type_theta
INTEGER                   :: grid_type_u,field_type_u
INTEGER                   :: grid_type_v,field_type_v
! Derived variables
REAL, ALLOCATABLE  :: p_ugrid_rho_levels(:,:,:)  ! pressure on u grid
REAL, ALLOCATABLE  :: p_vgrid_rho_levels(:,:,:)  ! pressure on v grid

! Tropopause pressure (and derivations)
REAL, ALLOCATABLE  :: trop_pressure(:,:)       ! tropopause pressure on T grid
REAL, ALLOCATABLE  :: trop_pressure_ugrid(:,:) ! tropopause pressure on u grid
REAL, ALLOCATABLE  :: trop_pressure_vgrid(:,:) ! tropopause pressure on v grid

! Variables associated with the above
INTEGER   :: theta_row_length_min   ! Min theta column (no halo)
INTEGER   :: theta_row_length_max   ! Max theta column
INTEGER   :: theta_rows_min         ! Min theta row (no halo)
INTEGER   :: theta_rows_max         ! Max theta row

INTEGER   :: u_row_length_min       ! Min U column
INTEGER   :: u_row_length_max       ! Max U column
INTEGER   :: u_rows_min             ! Min U row
INTEGER   :: u_rows_max             ! Max U row

INTEGER   :: v_row_length_min       ! Min V column
INTEGER   :: v_row_length_max       ! Max V column
INTEGER   :: v_rows_min             ! Min V row
INTEGER   :: v_rows_max             ! Max V row

REAL      :: frac                   ! Frac between ECMWF timesteps
INTEGER   :: proc_row_length_min    ! min column in PE
INTEGER   :: proc_trow_length_max   ! max column T in PE
INTEGER   :: proc_urow_length_max   ! max column U in PE
INTEGER   :: proc_vrow_length_max   ! max column V in PE
INTEGER   :: proc_rows_min          ! min row in PE
INTEGER   :: proc_trows_max         ! max row in PE
INTEGER   :: proc_urows_max         ! max row in PE
INTEGER   :: proc_vrows_max         ! max row in PE

! ***********************************************************
! Diagnostics

REAL, ALLOCATABLE   :: diag_theta_data   (:,:,:)     ! T diag 1: ECMWf T on mod levs
REAL, ALLOCATABLE   :: diag_theta_model  (:,:,:)     ! T diag 2: model T on mod levs
REAL, ALLOCATABLE   :: diag_theta_ntend  (:,:,:)     ! T diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_theta_mtend  (:,:,:)     ! T diag 4: model tend tendency
REAL, ALLOCATABLE   :: diag_theta_relax  (:,:,:)     ! T diag 5: relax par
REAL, ALLOCATABLE   :: diagsq_theta_ntend(:,:,:)     ! T diag 6: nudging tend^2
REAL, ALLOCATABLE   :: diagsq_theta_mtend(:,:,:)     ! T diag 7: model tend^2

REAL, ALLOCATABLE   :: diag_u_data       (:,:,:)     ! u diag 1: ECMWf U on mod levs
REAL, ALLOCATABLE   :: diag_u_model      (:,:,:)     ! u diag 2: model U on mod levs
REAL, ALLOCATABLE   :: diag_u_ntend      (:,:,:)     ! u diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_u_mtend      (:,:,:)     ! u diag 4: model tendency
REAL, ALLOCATABLE   :: diag_u_relax      (:,:,:)     ! u diag 5: relax par
REAL, ALLOCATABLE   :: diagsq_u_ntend    (:,:,:)     ! u diag 6: nudging tend^2
REAL, ALLOCATABLE   :: diagsq_u_mtend    (:,:,:)     ! u diag 7: model tend^2

REAL, ALLOCATABLE   :: diag_v_data       (:,:,:)     ! v diag 1: ECMWf V on mod levs
REAL, ALLOCATABLE   :: diag_v_model      (:,:,:)     ! v diag 2: model V on mod levs
REAL, ALLOCATABLE   :: diag_v_ntend      (:,:,:)     ! v diag 3: nudging tendency
REAL, ALLOCATABLE   :: diag_v_mtend      (:,:,:)     ! v diag 4: model tendency
REAL, ALLOCATABLE   :: diag_v_relax      (:,:,:)     ! v diag 5: relax par

! Logicals associated with diagnostics
LOGICAL, SAVE :: l_t_data_diag, l_t_model_diag, l_t_ntend_diag,      &
                 l_t_mtend_diag, l_t_relax_diag     ! Temperature

LOGICAL, SAVE :: l_u_data_diag, l_u_model_diag, l_u_ntend_diag,      &
                 l_u_mtend_diag, l_u_relax_diag     ! U-wind

LOGICAL, SAVE :: l_v_data_diag, l_v_model_diag, l_v_ntend_diag,      &
                 l_v_mtend_diag, l_v_relax_diag     ! V-wind

! *********************************************************************
! Variables associated with the netcdf files

REAL, SAVE, ALLOCATABLE   :: temp_file_data        (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: temp_file_data2       (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: u_file_data           (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: u_file_data2          (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: v_file_data           (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: v_file_data2          (:,:,:,:)

REAL, SAVE, ALLOCATABLE   :: surf_logp_file_data   (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logp_file_data2  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpu_file_data  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpu_file_data2 (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpv_file_data  (:,:,:,:)
REAL, SAVE, ALLOCATABLE   :: surf_logpv_file_data2 (:,:,:,:)

!     Dimensions of Netcdf Grid
INTEGER                   :: file_dims(totdims)   ! Array of all dim lengths
INTEGER                   :: file_dims2(totdims)
INTEGER, SAVE             :: file_row_length      ! no. of (T) pts in a row
INTEGER, SAVE             :: file_rows            ! no of (T) rows
INTEGER, SAVE             :: file_levels          ! no. of model levels
INTEGER, SAVE             :: file_timesteps       ! no. of timesteps
INTEGER, SAVE             :: file_urow_length     ! no. of (U) pts in a row
INTEGER, SAVE             :: file_vrows           ! no of (V) rows
INTEGER, SAVE             :: timestep1            ! Timestep of first data point
INTEGER, SAVE             :: timestep2            ! Timestep of second data poin
CHARACTER(LEN=256)        :: dataname1            ! first file+pathname
CHARACTER(LEN=256)        :: dataname2            ! second file+pathname

REAL                      :: base_phi_theta       ! base phi for theta grid
REAL                      :: base_phi_v           ! base phi for v grid

INTEGER                   :: icode                ! return code

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_MAIN1'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ************************************************************************
! End of the Header

! Standard subroutine entry comment
IF (PrintStatus > PrStatus_Normal) THEN
  CALL umPrint('NUDGING_MAIN: Entering routine ', &
      src='nudging_main1_nudging_main1')
END IF

! Save tendencies into physics_tendencies_mod
IF (l_retain_nud_tendencies) THEN

  ! Allocate variables
  CALL init_nud_tendencies()
  
  ! create increments
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_nud(i,j,k) = theta(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_nud(i,j,k) = u(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_nud(i,j,k) = v(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

END IF ! end if loop over l_retain_nud_tendencies

! *********************************************************

!     Calculate various miscellanea before we start

! *********************************************************

IF ( l_first )  THEN
! Confirm that the namelist Run_nudging has been specified properly & read
 IF ( ndg_relax_uvalue < 0.0 .OR. ndg_relax_vvalue < 0.0 .OR.          &
    ndg_relax_tvalue < 0.0 .OR. ndg_lev_bottom   <  1  .OR.            &
    ndg_lev_top < 1 .OR. ndg_on_lev_bottom < 1         .OR.            &
    ndg_on_lev_top < 1 .OR.  ndg_datapath == ' '       .OR.            &
    ndg_lev_top > model_levels .OR. ndg_lev_bottom > model_levels .OR. &
    ndg_strat_fac < 0.0 .OR. ndg_strat_fac > 1.0         .OR.            &
    ndg_hours_perdata < 1 .OR. ndg_analysis_source > 3 ) THEN   

  icode = 999

  CALL ereport ('NUDGING_MAIN',icode,                                  &
     'Improper values found in namelist: RUN_Nudging')
 END IF

! Check whether 'Pressure at Tropopause Level' (30-451) has been requested   
! in the STASH. This is required for determination of tropopause to apply
! the Stratospheric Nudging Factor'.
! trop_pres_index,etc are set in Nudging_control module.
 IF ( .NOT. sf(trop_pres_index,sect_tropht) ) THEN
   icode = 451
   WRITE (umMessage,'(A)') 'Pressure at Tropopause Level needs to be '//&
   'requested in STASH with TALLTS profile'
   CALL umPrint(umMessage,src='nudging_main1_nudging_main1')
   CALL ereport ('NUDGING_MAIN',icode,                                  &
      'Required Diagnostic missing from STASH')
 END IF

 CALL nudging_set_diag_logic

END IF   ! l_first

! Load the submodel code
m_atm_modl   = submodel_for_sm(atmos_im)
im_index     = 1

! Load the number of dynamic timesteps per hour
a_steps_per_hr  = INT(3600.0/timestep)

! ******************************************************
! Load where the local (processor) indices are in the global array
! using Datastart and Blsize from ParVars
proc_row_length_min  = datastart(1)
proc_trow_length_max = proc_row_length_min + blsize(1,fld_type_p) -1
proc_urow_length_max = proc_row_length_min + blsize(1,fld_type_u) -1
proc_vrow_length_max = proc_row_length_min + blsize(1,fld_type_v) -1

proc_rows_min  = datastart(2)
proc_trows_max = proc_rows_min + blsize(2,fld_type_p) -1
proc_urows_max = proc_rows_min + blsize(2,fld_type_u) -1
proc_vrows_max = proc_rows_min + blsize(2,fld_type_v) -1

! If in debug mode write out where we are on the grid
IF (l_first .AND. PrintStatus > PrStatus_Oper) THEN
  WRITE(umMessage,'(A,8(1X,I4))' )                                             &
   ' NUDGING_MAIN: Minimum and Maximum Indices are',                           &
   proc_row_length_min, proc_trow_length_max, proc_urow_length_max,            &
   proc_vrow_length_max, proc_rows_min, proc_trows_max, proc_urows_max,        &
   proc_vrows_max
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')
END IF

! Define location for debug prints --approx centre of local domain
dbg_x = (proc_trow_length_max - proc_row_length_min)/2
dbg_y = (proc_trows_max - proc_rows_min)/2
dbg_z = (ndg_lev_top - ndg_lev_bottom)/2

! Define the size of each grid type we are using (no halo).

u_row_length_min = udims%i_start
u_row_length_max = udims%i_end
u_rows_min       = udims%j_start
u_rows_max       = udims%j_end

v_row_length_min = vdims%i_start
v_row_length_max = vdims%i_end
v_rows_min       = vdims%j_start
v_rows_max       = vdims%j_end

theta_row_length_min = tdims%i_start
theta_row_length_max = tdims%i_end
theta_rows_min       = tdims%j_start
theta_rows_max       = tdims%j_end

! Define shifted base_phi for endgame or new dynamics

base_phi_theta = base_phi + (0.5*delta_phi)
base_phi_v = base_phi

! If in debug mode provide information about dimensions
IF (l_first .AND. PrintStatus > PrStatus_Oper) THEN
  WRITE(umMessage,'(A,4(1X,I4))')                                              &
    ': NUDGING_MAIN: Local u array bounds = ',                                 &
    u_row_length_min,  u_row_length_max,                                       &
    u_rows_min,  u_rows_max
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')

  WRITE(umMessage,'(A,4(1X,I4))')                                             &
   ': NUDGING_MAIN: Local v array bounds = ',                                 &
   v_row_length_min,  v_row_length_max,                                       &
   v_rows_min,  v_rows_max
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')

  WRITE(umMessage, '(A,4(1X,I4))')                                             &
    ': NUDGING_MAIN: Local theta array bounds = ',                             &
    theta_row_length_min,  theta_row_length_max,                               &
    theta_rows_min,  theta_rows_max
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')

END IF ! debug output

! ******************************************************

! Allocate required variables - diagnostics & local arrays
CALL nudging_alloc_data

! Load required diagnostic & u,v,theta grid information
! Need to access D1 for this
CALL nudging_getd1_data

! *******************************************************************
! Migrate Pressure values from P to U,V grids

p_ugrid_rho_levels(:,:,:) = 0.0
p_vgrid_rho_levels(:,:,:) = 0.0

CALL swap_bounds(p_rho_levels,pdims%i_end,pdims%j_end,          &
                   pdims%k_end+1,offx,offy,fld_type_p,swap_field_is_scalar)

CALL p_to_u(p_rho_levels,                                       &
            pdims_s%i_start,pdims_s%i_end,                      &
            pdims_s%j_start,pdims_s%j_end,                      &
            udims%i_start,udims%i_end,                          &
            udims%j_start,udims%j_end,                          &
            udims%k_start,udims%k_end+1,                        &
            p_ugrid_rho_levels)

CALL p_to_v(p_rho_levels,                                       &
            pdims_s%i_start,pdims_s%i_end,                      &
            pdims_s%j_start,pdims_s%j_end,                      &
            vdims%i_start,vdims%i_end,                          &
            vdims%j_start,vdims%j_end,                          &
            vdims%k_start,vdims%k_end+1,                        &
            p_vgrid_rho_levels)

CALL swap_bounds(trop_pressure,pdims%i_end,pdims%j_end,          &
                   0,offx,offy,fld_type_p,swap_field_is_scalar)

CALL p_to_u(trop_pressure,                                      &
            pdims_s%i_start,pdims_s%i_end,                      &
            pdims_s%j_start,pdims_s%j_end,                      &
            udims%i_start,udims%i_end,                          &
            udims%j_start,udims%j_end,                          &
            0,0,                                                &
            trop_pressure_ugrid)

CALL p_to_v(trop_pressure,                                      &
            pdims_s%i_start,pdims_s%i_end,                      &
            pdims_s%j_start,pdims_s%j_end,                      &
            vdims%i_start,vdims%i_end,                          &
            vdims%j_start,vdims%j_end,                          &
            0,0,                                                &
            trop_pressure_vgrid)

! ***************************************************************
! Load netcdf files and obtain basic information

! Get fraction between the two ECMWF timesteps of the model timestep
CALL nudging_getfrac(                 &
 a_steps_per_hr,                      & ! No. of model timesteps per hour
 ndg_hours_perdata,                   & ! No. of hours per data timestep
 i_hour,                              & ! Model hour
 i_minute,                            & ! Model minute
 frac)                                  ! Fraction between data timesteps

! If in debug mode then write out the value of frac
IF (PrintStatus > PrStatus_Oper) THEN
  WRITE(umMessage,'(A,F12.4)')                                                 &
   'NUDGING_MAIN: Fraction between data steps equals ', frac
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')
END IF

! If we have started a new data timestep then reload data
! NB this assumes that the data timestep is made up of an
! integral no. of atm steps
IF (frac == 0 .OR. l_first) CALL nudging_netcdf_loader

! *********************************************************
! Call the routine that controls the nudging -only if parameter
! is to be nudged
! Potential temperature
IF ( ndg_relax_tvalue > 0.0 )             &
 CALL nudging_nudging_control(            &
  temp_name,                              &
  grid_type_theta,                        & ! Grid type
  ana_row_length,                         & ! Global row length
  ana_rows,                               & ! Global rows
  file_levels,                            & ! Netcdf file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_trow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_trows_max,                         & ! Maximum row
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  base_lambda,                            &
  delta_lambda,                           &
  base_phi_theta,                         &
  delta_phi,                              &
  ana_base_lambda,                        &
  ana_delta_lambda,                       &
  ana_base_phi,                           &
  ana_delta_phi,                          &
  temp_file_data       (:,:,:,timestep1), &
  temp_file_data2      (:,:,:,timestep2), &
  surf_logp_file_data  (:,:,1,timestep1), &
  surf_logp_file_data2 (:,:,1,timestep2), &
  frac,                                   & ! Fraction between data
  model_levels,                           & ! model levels
  p_theta_levels (                        & ! model pressure
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max,          &
    1:model_levels ),                     &
  trop_pressure(                          & ! Tropopause pressure
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max),         &
  theta (                                 & ! nudged variable (theta)
  theta_row_length_min:theta_row_length_max, &
  theta_rows_min:theta_rows_max,          &
    1:model_levels ),                     &
  l_t_data_diag, diag_theta_data,         & ! theta diag 1
  l_t_model_diag, diag_theta_model,       & ! theta diag 2
  l_t_ntend_diag, diag_theta_ntend,       & ! theta diag 3
  l_t_mtend_diag, diag_theta_mtend,       & ! theta diag 4
  l_t_relax_diag, diag_theta_relax)         ! theta diag 5

! *******************************************************************
! Zonal Wind
IF ( ndg_relax_uvalue > 0.0 )             &
CALL nudging_nudging_control(             &
  u_name,                                 &
  grid_type_u,                            & ! Grid type
  ana_row_length,                         & ! Global row length
  ana_rows,                               & ! Global rows
  file_levels,                            & ! Netcfd file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_urow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_urows_max,                         & ! Maximum row
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  (base_lambda+0.5*delta_lambda),         &
  delta_lambda,                           &
  base_phi_theta,                         &
  delta_phi,                              &
  ana_base_lambda_u,                      &
  ana_delta_lambda,                       &
  ana_base_phi,                           &
  ana_delta_phi,                          &
  u_file_data(:,:,:,timestep1),           &
  u_file_data2(:,:,:,timestep2),          &
  surf_logpu_file_data(:,:,1,timestep1),  &
  surf_logpu_file_data2(:,:,1,timestep2), &
  frac,                                   & ! fraction between data
  model_levels,                           & ! model levels
  p_ugrid_rho_levels                      & ! model pressure
   ( u_row_length_min:u_row_length_max,   &
     u_rows_min:u_rows_max,               &
     1:model_levels ),                    &
  trop_pressure_ugrid(                    & ! Tropopause pressure
     u_row_length_min:u_row_length_max,   &
     u_rows_min:u_rows_max),              &
  u (                                     & ! nudged variable (u)
   u_row_length_min:u_row_length_max,     &
   u_rows_min:u_rows_max,                 &
   1:model_levels ),                      &
  l_u_data_diag, diag_u_data,             & ! u diag 1
  l_u_model_diag, diag_u_model,           & ! u diag 2
  l_u_ntend_diag, diag_u_ntend,           & ! u diag 3
  l_u_mtend_diag, diag_u_mtend,           & ! u diag 4
  l_u_relax_diag, diag_u_relax)             ! u diag 5

! *******************************************************************
! Meridional Wind
IF ( ndg_relax_vvalue > 0.0 )             &
CALL nudging_nudging_control(             &
  v_name,                                 &
  grid_type_v,                            & ! Grid type
  ana_row_length,                         & ! global row length
  (ana_rows-1),                           & ! global rows
  file_levels,                            & ! Netcdf file levels
  proc_row_length_min,                    & ! Mimimum column
  proc_vrow_length_max,                   & ! Maximum column
  proc_rows_min,                          & ! Minimum row
  proc_vrows_max,                         & ! Maximum row
  sin_theta_latitude(1,1),                & ! sine of min latitude
  sin_theta_latitude(row_length, rows),   & ! sine of max latitude
  base_lambda,                            &
  delta_lambda,                           &
  base_phi_v,                             &
  delta_phi,                              &
  ana_base_lambda,                        &
  ana_delta_lambda,                       &
  ana_base_phi_v,                         &
  ana_delta_phi,                          &
  v_file_data(:,:,:,timestep1),           &
  v_file_data2(:,:,:,timestep2),          &
  surf_logpv_file_data(:,:,1,timestep1),  &
  surf_logpv_file_data2(:,:,1,timestep2), &
  frac,                                   & ! fraction between data
  model_levels,                           & ! Model levels
  p_vgrid_rho_levels                      & ! Model pressure
   ( v_row_length_min:v_row_length_max,   &
     v_rows_min:v_rows_max,               &
     1:model_levels ),                    &
  trop_pressure_vgrid(                    & ! Tropopause pressure
     v_row_length_min:v_row_length_max,   &
     v_rows_min:v_rows_max),              &
  v (                                     & ! nudged variable (v)
    v_row_length_min:v_row_length_max,    &
    v_rows_min:v_rows_max,                &
    1:model_levels ),                     &
  l_v_data_diag, diag_v_data,             & ! v diag 1
  l_v_model_diag, diag_v_model,           & ! v diag 2
  l_v_ntend_diag, diag_v_ntend,           & ! v diag 3
  l_v_mtend_diag, diag_v_mtend,           & ! v diag 4
  l_v_relax_diag, diag_v_relax)             ! v diag 5

! ***************************************************
! Wrap up at the end

! *****************************************************************


! Synchronize V for the polar rows, call as in ni_filter_sl
IF( model_type == mt_global .AND. ndg_relax_vvalue > 0.0 ) THEN 
  IF( at_extremity(psouth) ) THEN 
 
    CALL eg_v_at_poles(u, v, 1.0, udims%j_start, vdims%j_start,    & 
                       udims_s,vdims_s) 
  END IF 
 
  IF( at_extremity(pnorth) ) THEN 
 
    CALL eg_v_at_poles(u, v, -1.0, udims%j_end, vdims%j_end,       & 
                        udims_s,vdims_s) 
  END IF 
END IF      ! Model = global

! Upload variables to the STASHWORK array
CALL nudging_update_stash

! Deallocate arrays used
CALL nudging_dealloc_data

! Copy tendencies to physics_tendencies
IF (l_retain_nud_tendencies) THEN
  
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_nud(i,j,k) = theta(i,j,k) - dtheta_nud(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_nud(i,j,k) = u(i,j,k) - du_nud(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_nud(i,j,k) = v(i,j,k) - dv_nud(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

END IF ! end loop over l_retain_nud_tendencies

! No longer the first time so set to false
l_first = .FALSE.

! Standard subroutine exit comment
IF (PrintStatus > PrStatus_Normal) CALL umPrint(' Leaving NUDGING_MAIN', &
    src='nudging_main1_nudging_main1')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
! ******************************************************

RETURN

! *****************************************************
! Contained subroutines

CONTAINS

! ********************************************************************!
!                                                                     !
!     This routine sets logicals that control calculation of the      !
!     Nudging diagnostics depending on the STASH requests.            !
!                                                                     !
! ********************************************************************!
SUBROUTINE nudging_set_diag_logic

IMPLICIT NONE

REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_SET_DIAG_LOGIC'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set value of logicals depending on STASH Flag. 
! Indices and sect_nudge are set in Nudging_control module

l_t_data_diag = sf(tdiag_data_index,sect_nudge)
l_t_model_diag = sf(tdiag_model_index,sect_nudge)
l_t_ntend_diag = sf(tdiag_ntend_index,sect_nudge)
l_t_mtend_diag = sf(tdiag_mtend_index,sect_nudge)
l_t_relax_diag = sf(tdiag_relax_index,sect_nudge)

l_u_data_diag = sf(udiag_data_index,sect_nudge)
l_u_model_diag = sf(udiag_model_index,sect_nudge)
l_u_ntend_diag = sf(udiag_ntend_index,sect_nudge)
l_u_mtend_diag = sf(udiag_mtend_index,sect_nudge)
l_u_relax_diag = sf(udiag_relax_index,sect_nudge)

l_v_data_diag = sf(vdiag_data_index,sect_nudge)
l_v_model_diag = sf(vdiag_model_index,sect_nudge)
l_v_ntend_diag = sf(vdiag_ntend_index,sect_nudge)
l_v_mtend_diag = sf(vdiag_mtend_index,sect_nudge)
l_v_relax_diag = sf(vdiag_relax_index,sect_nudge)

! Check whether certain Nudging diagnostics have been requested in STASH.
! The only items required are the 'XXXX AFTER NUDGING' fields
! (Xdiag_model_index), if the corresponding 'XXXX INCREMENT DUE TO OTHER' &
! diagnostics (Xdiag_mtend_index) are requested.
IF ( (l_t_mtend_diag .AND. .NOT. l_t_model_diag ) .OR.       &
     (l_u_mtend_diag .AND. .NOT. l_u_model_diag ) .OR.       &
     (l_v_mtend_diag .AND. .NOT. l_v_model_diag ) ) THEN
  
  icode = sect_nudge
  
  WRITE (umMessage,'(A)') 'A "XXXX INCREMENT DUE TO OTHER" diag '// &
   'has been requested without the "XXXX AFTER NUDGING" diagnostic'
  CALL umPrint(umMessage,src='nudging_main1_nudging_main1')

  CALL ereport('NUDGING_SET_DIAG_LOGIC',icode,                              &
    'Required Nudging diagnostic not in STASH' )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nudging_set_diag_logic

! ********************************************************************!
!                                                                     !
!     This routine loads data from netcdf files, mostly using         !
!     procedures from the UKCA emissions I/O module.                  !
!     Separated from the main routine for the sake of tidiness        !
!                                                                     !
! ********************************************************************!
SUBROUTINE nudging_netcdf_loader

USE emiss_io_mod,   ONLY:em_fopen,em_fclose,em_get_var_info
USE nudging_io_mod

IMPLICIT NONE

INTEGER              :: ndfile1, ndfile2      ! File indices
INTEGER              :: varid,ndimen          ! Variable id, n-dimensions
INTEGER              :: f_dims(max_ncdims)    ! File Variable dimensions
INTEGER              :: file_levp             ! Levels for surface press
                                              ! can be 2-D or 3-D  
CHARACTER(LEN=20)    :: varname               ! local copy of variable name
 
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_NETCDF_LOADER'

! **************************************
! End of Header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Load filenames for file immediately before and after model timestep
! Don't worry about ECMWF naming convention, only necessary to change
! if you get data from another source with another naming convention
CALL nudging_getfilename(               &
  i_year,                               &  ! Model year
  i_month,                              &  ! Model month
  i_day,                                &  ! Model day
  i_hour,                               &  ! Model hour
  dataname1,                            &  ! Return datafile name 1
  dataname2,                            &  ! Return datafile name 2
  timestep1,                            &  ! Return timestep 1
  timestep2,                            &  ! Return timestep 2
  0)                                       ! Select ECMWF naming convention

! If in debug mode drop a line giving the filenames
IF (PrintStatus > PrStatus_Oper) THEN
  WRITE(umMessage,'(A,A)')'NUDGING_NETCDF_LOADER: filename 1: ',dataname1
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')
  WRITE(umMessage,'(A,A)')'NUDGING_NETCDF_LOADER: filename 2: ',dataname2
  CALL umPrint(umMessage,src='nudging_main1-nudging_main1')
END IF

! Open netcdf files and load dimensions, using UKCA module procedures.
! Use variable names as defined in nudging_control module.
! --assume x,y,z dims are identical in both files. The sizes are anyway
! checked when variable is read in by NDG_GET_DATA. 

! 'ignore_time=True' prevents the routine from trying to get name of
!  time variable, which is specific to UKCA files only.
CALL em_fopen(dataname1,ndfile1,ignore_time=.TRUE.)

! Clear data arrays
IF (ALLOCATED(temp_file_data))        DEALLOCATE(temp_file_data)
IF (ALLOCATED(temp_file_data2))       DEALLOCATE(temp_file_data2)
IF (ALLOCATED(u_file_data))           DEALLOCATE(u_file_data)
IF (ALLOCATED(u_file_data2))          DEALLOCATE(u_file_data2)
IF (ALLOCATED(v_file_data))           DEALLOCATE(v_file_data)
IF (ALLOCATED(v_file_data2))          DEALLOCATE(v_file_data2)
IF (ALLOCATED(surf_logp_file_data))   DEALLOCATE(surf_logp_file_data)
IF (ALLOCATED(surf_logp_file_data2))  DEALLOCATE(surf_logp_file_data2)
IF (ALLOCATED(surf_logpu_file_data))  DEALLOCATE(surf_logpu_file_data)
IF (ALLOCATED(surf_logpu_file_data2)) DEALLOCATE(surf_logpu_file_data2)
IF (ALLOCATED(surf_logpv_file_data))  DEALLOCATE(surf_logpv_file_data)
IF (ALLOCATED(surf_logpv_file_data2)) DEALLOCATE(surf_logpv_file_data2)

varname = temp_name ! local copy since decl as INTENT(INOUT) in routine
varid = -1   ! Set to negative so ignored by routine
f_dims(:) = 0

CALL em_get_var_info(ndfile1, varid, varname, ndims=ndimen, vdimsize=f_dims) 

file_levels = f_dims(3)  ! Variable used in press-to-model_levs routine

! Compare sizes against those expected for Analysis data
! Levels can be different as they are interpolated onto model levels anyway
IF ( f_dims(1) /= ana_row_length .OR. f_dims(2) /= ana_rows ) THEN
   WRITE(umMessage,'(A,2I6,A,2I6)')                                    &
    'Unexpected T dimensions in file '//TRIM(dataname1)//'. Expected:',&
    ana_row_length, ana_rows,' found: ',f_dims(1), f_dims(2)
   CALL umPrint(umMessage,src='NUDGING_NETCDF_LOADER')
   icode = 999
   CALL ereport ('NUDGING_NETCDF_LOADER',icode,                        &
     'Unexpected T dimensions in file '//TRIM(dataname1))
END IF

! Allocate T/P grid variables - size = 1 for time 
ALLOCATE(temp_file_data(f_dims(1),f_dims(2),f_dims(3),1))
ALLOCATE(temp_file_data2(f_dims(1),f_dims(2),f_dims(3),1))
! Get number of levels for Log(Surface_press) - can be 2-D or 3-D depending 
! in source of dataset
varname = surfp_name 
varid = -1
CALL em_get_var_info(ndfile1, varid, varname, ndims=ndimen, vdimsize=f_dims) 

file_levp = f_dims(3)   ! Store levels for psurf, use for surfpu,surfpv
ALLOCATE(surf_logp_file_data(f_dims(1),f_dims(2),file_levp,1))
ALLOCATE(surf_logp_file_data2(f_dims(1),f_dims(2),file_levp,1))

! Obtain and allocate U-grid variables (same as T currently)
varname = u_name 
varid = -1
CALL em_get_var_info(ndfile1, varid, varname, ndims=ndimen, vdimsize=f_dims) 

ALLOCATE(u_file_data(f_dims(1),f_dims(2),f_dims(3),1))
ALLOCATE(u_file_data2(f_dims(1),f_dims(2),f_dims(3),1))
ALLOCATE(surf_logpu_file_data(f_dims(1),f_dims(2),file_levp,1))
ALLOCATE(surf_logpu_file_data2(f_dims(1),f_dims(2),file_levp,1))

! Obtain and allocate V-grid variables (rows-1 currently)
varname = v_name 
varid = -1
CALL em_get_var_info(ndfile1, varid, varname, ndims=ndimen, vdimsize=f_dims) 

! Compare sizes against those expected for Analysis data
IF ( f_dims(1) /= ana_row_length .OR. f_dims(2) /= (ana_rows-1) ) THEN
   WRITE(umMessage,'(A,2I6,A,2I6)')                                    &
    'Unexpected V dimensions in file '//TRIM(dataname1)//'. Expected:',&
    ana_row_length, (ana_rows-1),' found: ',f_dims(1), f_dims(2)
   CALL umPrint(umMessage,src='NUDGING_NETCDF_LOADER')
   icode = 999
   CALL ereport ('NUDGING_NETCDF_LOADER',icode,                        &
     'Unexpected V dimensions in file '//TRIM(dataname1))
END IF

ALLOCATE(v_file_data(f_dims(1),f_dims(2),f_dims(3),1))
ALLOCATE(v_file_data2(f_dims(1),f_dims(2),f_dims(3),1))
ALLOCATE(surf_logpv_file_data(f_dims(1),f_dims(2),file_levp,1))
ALLOCATE(surf_logpv_file_data2(f_dims(1),f_dims(2),file_levp,1))

! Read values from file1 
CALL ndg_get_data(ana_row_length, ana_rows, f_dims(3), 1, ndfile1,      &
  timestep1, temp_name, temp_file_data )          ! Temperature
CALL ndg_get_data(ana_row_length, ana_rows, f_dims(3), 1, ndfile1,      &
  timestep1, u_name, u_file_data )                ! U wind
CALL ndg_get_data(ana_row_length, (ana_rows-1), f_dims(3), 1, ndfile1,  &
  timestep1, v_name, v_file_data )                ! V wind
CALL ndg_get_data(ana_row_length, ana_rows, file_levp, 1, ndfile1,      &
  timestep1, surfp_name, surf_logp_file_data )    ! Log(Surfp)
CALL ndg_get_data(ana_row_length, ana_rows, file_levp, 1, ndfile1,      &
  timestep1, usurfp_name, surf_logpu_file_data )  ! Log(Surfp)-U grid
CALL ndg_get_data(ana_row_length, (ana_rows-1), file_levp, 1, ndfile1,  &
  timestep1, vsurfp_name, surf_logpv_file_data )  ! Log(Surfp)-V grid
CALL em_fclose(ndfile1)

! Read values from file2 -- same dimensions for all files from a dataset
CALL em_fopen(dataname2,ndfile2,ignore_time=.TRUE.)

CALL ndg_get_data(ana_row_length, ana_rows, f_dims(3), 1, ndfile2,      &
  timestep2, temp_name, temp_file_data2 )          ! Temperature
CALL ndg_get_data(ana_row_length, ana_rows, f_dims(3), 1, ndfile2,      &
  timestep2, u_name, u_file_data2 )                ! U wind
CALL ndg_get_data(ana_row_length, (ana_rows-1), f_dims(3), 1, ndfile2,  &
  timestep2, v_name, v_file_data2 )                ! V wind
CALL ndg_get_data(ana_row_length, ana_rows, file_levp, 1, ndfile2,      &
  timestep2, surfp_name, surf_logp_file_data2 )    ! Log(Surfp)
CALL ndg_get_data(ana_row_length, ana_rows, file_levp, 1, ndfile2,      &
  timestep2, usurfp_name, surf_logpu_file_data2 )  ! Log(Surfp)-U grid
CALL ndg_get_data(ana_row_length, (ana_rows-1), file_levp, 1, ndfile2,  &
  timestep2, vsurfp_name, surf_logpv_file_data2 )  ! Log(Surfp)-V grid
CALL em_fclose(ndfile2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nudging_netcdf_loader

! ***********************************************************!
!                                                            !
!     This routine allocates required diagnotics and         !
!     temporary / local arrays                               !
!                                                            !
! ******************************************** **************!
SUBROUTINE nudging_alloc_data

USE ereport_mod, ONLY: ereport
USE Field_Types
USE UM_ParVars
IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_ALLOC_DATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate the theta diag1 array
ALLOCATE(diag_theta_data(                                               &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_data(:,:,:) = 0.0

! Allocate the theta diag2 array
ALLOCATE(diag_theta_model(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_model(:,:,:) = 0.0

! Allocate the theta diag3 array
ALLOCATE(diag_theta_mtend(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_mtend(:,:,:) = 0.0

! Allocate the theta diag4 array
ALLOCATE(diag_theta_ntend(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_ntend(:,:,:) = 0.0

! Allocate the theta diag5 array
ALLOCATE(diag_theta_relax(                                              &
      theta_row_length_min:theta_row_length_max,                        &
      theta_rows_min:theta_rows_max,                                    &
      1:model_levels))

diag_theta_relax(:,:,:) = 0.0

! Allocate the u diag1 array
ALLOCATE(diag_u_data(                                                   &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_data(:,:,:) = 0.0

! Allocate the u diag2 array
ALLOCATE(diag_u_model(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_model(:,:,:) = 0.0

! Allocate the u diag3 array
ALLOCATE(diag_u_mtend(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_mtend(:,:,:) = 0.0

! Allocate the u diag4 array
ALLOCATE(diag_u_ntend(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_ntend(:,:,:) = 0.0

! Allocate the u diag5 array
ALLOCATE(diag_u_relax(                                                  &
      u_row_length_min:u_row_length_max,                                &
      u_rows_min:u_rows_max,                                            &
      1:model_levels))

diag_u_relax(:,:,:) = 0.0

! Allocate the v diag1 array
ALLOCATE(diag_v_data(                                                   &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_data(:,:,:) = 0.0

! Allocate the v diag2 array
ALLOCATE(diag_v_model(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_model(:,:,:) = 0.0

! Allocate the v diag3 array
ALLOCATE(diag_v_mtend(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_mtend(:,:,:) = 0.0

! Allocate the v diag4 array
ALLOCATE(diag_v_ntend(                                                  &
      v_row_length_min:v_row_length_max,                                &
      v_rows_min:v_rows_max,                                            &
      1:model_levels))

diag_v_ntend(:,:,:) = 0.0

! Allocate the v diag5 array
ALLOCATE(diag_v_relax(                                                 &
      v_row_length_min:v_row_length_max,                               &
      v_rows_min:v_rows_max,                                           &
      1:model_levels))

diag_v_relax(:,:,:) = 0.0

! Allocate pressure on rho levels on the u and v grids
ALLOCATE(p_ugrid_rho_levels(                                           &
      u_row_length_min : u_row_length_max,                             &
      u_rows_min : u_rows_max,                                         &
      udims%k_start:udims%k_end+1))

p_ugrid_rho_levels(:,:,:) = 0.0

ALLOCATE(p_vgrid_rho_levels(                                           &
     v_row_length_min : v_row_length_max,                              &
     v_rows_min : v_rows_max,                                          &
     vdims%k_start:vdims%k_end+1))

p_vgrid_rho_levels(:,:,:) = 0.0

 ! Allocate the tropopause pressure array
ALLOCATE(trop_pressure(                                                &
     pdims_s%i_start:pdims_s%i_end,                                    &
     pdims_s%j_start:pdims_s%j_end))

trop_pressure(:,:) = 0.0

! Allocate the tropopause pressure array
ALLOCATE(trop_pressure_ugrid(                                          &
     u_row_length_min:u_row_length_max,                                &
     u_rows_min:u_rows_max))

trop_pressure_ugrid(:,:) = 0.0

! Allocate the tropopause pressure array
ALLOCATE(trop_pressure_vgrid(                                          &
     v_row_length_min:v_row_length_max,                                &
     v_rows_min:v_rows_max))

trop_pressure_vgrid(:,:) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nudging_alloc_data

! ***********************************************************!
!                                                            !
!     This routine downloads specific variables from the     !
!     D1 array into specific arrays                          !
!     Section & Item numbers are currently set in            !
!     NUDGING_D1_DEFS module                                 !
!                                                            !
! ******************************************** **************!
SUBROUTINE nudging_getd1_data

USE ereport_mod, ONLY: ereport
USE Field_Types
USE UM_ParVars
IMPLICIT NONE

INTEGER           :: d1_addr_st,d1_var_len
                    ! start & length of variable in D1
INTEGER           :: icode
REAL,ALLOCATABLE  :: buff_temp(:)
REAL(KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_GETD1_DATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get grid type information for prognostic variables
! Currently required in some nudging subroutines

! Loop over sections & Items
DO i = 1, no_obj_d1(m_atm_modl)

  IF ( d1_addr(d1_section,i,m_atm_modl) == sect_prog ) THEN
    ! Section for prognostic
    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
      ! Item no.s
    CASE (u_index)
      grid_type_u = d1_addr(d1_grid_type,i,m_atm_modl)
      field_type_u = fld_type_u

    CASE (v_index)
      grid_type_v = d1_addr(d1_grid_type,i,m_atm_modl)
      field_type_v = fld_type_v

    CASE (theta_index)
      grid_type_theta = d1_addr(d1_grid_type,i,m_atm_modl)
      field_type_theta = fld_type_p

    CASE DEFAULT

      ! Do nothing
    END SELECT

  ELSE IF ( d1_addr(d1_section,i,m_atm_modl) == sect_tropht ) THEN
    !Section for climate diagnostics (for tropopause pressure)
    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
      ! Item no. for tropopause pressure (set in nudging_d1_defs module)

    CASE ( trop_pres_index)
      d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
      d1_var_len = d1_addr(d1_length,i,m_atm_modl)

      IF ( d1_var_len /= row_length*rows ) THEN
        WRITE(nmessage,'(3(A,1X,I9))')                        &
         'Mismatch in variable size :STASH vs trop pressure', &
          trop_pres_index,'. Expected',                       &
          row_length*rows,': Actual', d1_var_len
        icode = trop_pres_index
        CALL ereport('NUDGE_GET_D1',icode,nmessage)
      END IF

      ALLOCATE(buff_temp(row_length*rows))

      buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

      trop_pressure(pdims%i_start:pdims%i_end,               &
                    pdims%j_start:pdims%j_end) =             &
      RESHAPE(buff_temp(:),(/row_length,rows/))

      DEALLOCATE(buff_temp)

    END SELECT

    ! Extract values for post-nudged variables of previous timestep
    ! This should ensure bit comparability for CRUNs. Otherwise    &
    ! these would be initialised to zero and diag_model_tendency   &
    ! ( pre-nudging - previous_post-nudged) will get a different value
  ELSE IF ( d1_addr(d1_section,i,m_atm_modl) == sect_nudge ) THEN
    ! Section for nudging

    SELECT CASE( d1_addr(d1_item,i,m_atm_modl) )
      ! Item no.s (set in nudging_d1_defs module)

      ! Theta diag 2
    CASE ( tdiag_model_index )

      d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
      d1_var_len = d1_addr(d1_length,i,m_atm_modl)

      IF ( d1_var_len /= row_length*rows*model_levels ) THEN
        WRITE(nmessage,'(3(A,1X,I9))')                        &
        'Mismatch in variable size :STASH vs nudging diag',   &
         tdiag_model_index,'. Expected',                      &
         row_length*rows*model_levels,': Actual', d1_var_len
        icode = tdiag_model_index
        CALL ereport('NUDGE_GET_D1',icode,nmessage)
      END IF

      ALLOCATE(buff_temp(row_length*rows*model_levels))

      buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

      diag_theta_model(:,:,:) =                                &
       RESHAPE(buff_temp(:),(/row_length,rows,model_levels/))

      DEALLOCATE(buff_temp)

      ! u diag 2
    CASE ( udiag_model_index )

      d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
      d1_var_len = d1_addr(d1_length,i,m_atm_modl)

      IF ( d1_var_len /= row_length*rows*model_levels ) THEN
        WRITE(nmessage,'(3(A,1X,I9))')                         &
        'Mismatch in variable size :STASH vs nudging diag',    &
         udiag_model_index,'. Expected',                       &
         row_length*rows*model_levels,': Actual', d1_var_len
        icode = udiag_model_index
        CALL ereport('NUDGE_GET_D1',icode,nmessage)
      END IF

      ALLOCATE(buff_temp(row_length*rows*model_levels))

      buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

      diag_u_model(:,:,:) =                                &
       RESHAPE(buff_temp(:),(/row_length,rows,model_levels/))

      DEALLOCATE(buff_temp)

      ! v diag 2
    CASE ( vdiag_model_index )

      d1_addr_st = d1_addr(d1_address,i,m_atm_modl)
      d1_var_len = d1_addr(d1_length,i,m_atm_modl)

      IF ( d1_var_len /= row_length*n_rows*model_levels ) THEN
        WRITE(nmessage,'(3(A,1X,I9))')                         &
        'Mismatch in variable size :STASH vs nudging diag',    &
         vdiag_model_index,'. Expected',                       &
         row_length*rows*model_levels,': Actual', d1_var_len
        icode = vdiag_model_index
        CALL ereport('NUDGE_GET_D1',icode,nmessage)
      END IF

      ALLOCATE(buff_temp(row_length*n_rows*model_levels))

      buff_temp(:) = d1(d1_addr_st:d1_addr_st+d1_var_len-1)

      diag_v_model(:,:,:) =                                &
       RESHAPE(buff_temp(:),(/row_length,n_rows,model_levels/))

      DEALLOCATE(buff_temp)

    END SELECT ! Items in Nudging Section

  END IF       ! If Prognostics/Clim Diags/ Nudging Section

END DO         ! Loop over D1 items

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nudging_getd1_data

! ***********************************************************!
!                                                            !
!     This routine copies the nudging diagnostics into       !
!     STASH using the copy_diag routine                      !
!                                                            !
! ******************************************** **************!
SUBROUTINE nudging_update_stash

USE ereport_mod, ONLY: ereport
USE Field_Types
USE UM_ParVars
IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_UPDATE_STASH'
INTEGER :: errcode        ! Error code for ereport
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

! Analysis Theta on Model Grid
IF ( l_t_data_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(tdiag_data_index,sect_nudge,im_index)),        &
      diag_theta_data(:,:,:),                                   &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_data_index,sect_nudge,im_index)),&
      len_stlist, stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_data_index,icode,nmessage)

  IF ( icode > 0 ) THEN
    errcode = tdiag_data_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Theta After Nudging
IF ( l_t_model_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(tdiag_model_index,sect_nudge,im_index)),       &
      diag_theta_model(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_model_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,tdiag_model_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = tdiag_model_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Theta increment due to Nudging
IF ( l_t_ntend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(tdiag_ntend_index,sect_nudge,im_index)),       &
      diag_theta_ntend(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_ntend_index,icode,nmessage)

  IF ( icode > 0 ) THEN
    errcode = tdiag_ntend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Theta increment due to Other
IF ( l_t_mtend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(tdiag_mtend_index,sect_nudge,im_index)),       &
      diag_theta_mtend(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,tdiag_mtend_index,icode,nmessage)

  IF ( icode > 0 ) THEN
    errcode =   tdiag_mtend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Theta Relaxation parameter
IF ( l_t_relax_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(tdiag_relax_index,sect_nudge,im_index)),       &
      diag_theta_relax(:,:,:),                                  &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,tdiag_relax_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,tdiag_relax_index,icode,nmessage)

  IF ( icode > 0 ) THEN
   errcode = tdiag_relax_index
   CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Analysis U on Model Grid
IF ( l_u_data_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(udiag_data_index,sect_nudge,im_index)),        &
      diag_u_data(:,:,:),                                       &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_data_index,sect_nudge,im_index)),&
      len_stlist,stash_levels,num_stash_levels+1,               &
      atmos_im,sect_nudge,udiag_data_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = udiag_data_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! U After Nudging
IF ( l_u_model_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(udiag_model_index,sect_nudge,im_index)),       &
      diag_u_model(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_model_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_model_index,icode,nmessage)

  IF ( icode > 0 )  THEN
   errcode = udiag_model_index
   CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! U increment due to Nudging
IF ( l_u_ntend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(udiag_ntend_index,sect_nudge,im_index)),       &
      diag_u_ntend(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_ntend_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = udiag_ntend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! U increment due to Other
IF ( l_u_mtend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(udiag_mtend_index,sect_nudge,im_index)),       &
      diag_u_mtend(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_mtend_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = udiag_mtend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! U Relaxation parameter
IF ( l_u_relax_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(udiag_relax_index,sect_nudge,im_index)),       &
      diag_u_relax(:,:,:),                                      &
      row_length,rows,model_levels,0,0,0,0, at_extremity,       &
      stlist(1,stindex(1,udiag_relax_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,udiag_relax_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = udiag_relax_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! Analysis V on Model Grid
IF ( l_v_data_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(vdiag_data_index,sect_nudge,im_index)),        &
      diag_v_data(:,:,:),                                       &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_data_index,sect_nudge,im_index)),&
      len_stlist, stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,vdiag_data_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode  = vdiag_data_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! V after nudging
IF ( l_v_model_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(vdiag_model_index,sect_nudge,im_index)),       &
      diag_v_model(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_model_index,sect_nudge,im_index))&
      ,len_stlist,stash_levels,num_stash_levels+1,              &
      atmos_im,sect_nudge,vdiag_model_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = vdiag_model_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! V increment due to Nudging
IF ( l_v_ntend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(vdiag_ntend_index,sect_nudge,im_index)),       &
      diag_v_ntend(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_ntend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_ntend_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = vdiag_ntend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! V increment due to 'other'
IF ( l_v_mtend_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(vdiag_mtend_index,sect_nudge,im_index)),       &
      diag_v_mtend(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_mtend_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_mtend_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = vdiag_mtend_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF

! V Relaxation parameter
IF ( l_v_relax_diag )  THEN
! DEPENDS ON: copydiag_3d
 CALL copydiag_3d (                                             &
      st_work(si(vdiag_relax_index,sect_nudge,im_index)),       &
      diag_v_relax(:,:,:),                                      &
      row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
      stlist(1,stindex(1,vdiag_relax_index,sect_nudge,im_index))&
      ,len_stlist, stash_levels,num_stash_levels+1,             &
      atmos_im,sect_nudge,vdiag_relax_index,icode,nmessage)

  IF ( icode > 0 )  THEN
    errcode = vdiag_relax_index
    CALL ereport('NUDGE COPY DIAG',errcode,nmessage)
  END IF
END IF 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_update_stash

! ***********************************************************!
!                                                            !
!     This routine deallocates allocated diagnotics and      !
!     temporary / local arrays                               !
!                                                            !
! ******************************************** **************!
SUBROUTINE nudging_dealloc_data

USE ereport_mod, ONLY: ereport
USE Field_Types
USE UM_ParVars
IMPLICIT NONE
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUDGING_DEALLOC_DATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Deallocate derived arrays
IF (ALLOCATED(p_ugrid_rho_levels))  DEALLOCATE(p_ugrid_rho_levels)
IF (ALLOCATED(p_vgrid_rho_levels))  DEALLOCATE(p_vgrid_rho_levels)
IF (ALLOCATED(trop_pressure))       DEALLOCATE(trop_pressure)
IF (ALLOCATED(trop_pressure_ugrid)) DEALLOCATE(trop_pressure_ugrid)
IF (ALLOCATED(trop_pressure_vgrid)) DEALLOCATE(trop_pressure_vgrid)

! Deallocate diagnostic arrays
IF (ALLOCATED(diag_theta_data))     DEALLOCATE(diag_theta_data)
IF (ALLOCATED(diag_theta_model))    DEALLOCATE(diag_theta_model)
IF (ALLOCATED(diag_theta_ntend))    DEALLOCATE(diag_theta_ntend)
IF (ALLOCATED(diag_theta_mtend))    DEALLOCATE(diag_theta_mtend)
IF (ALLOCATED(diag_theta_relax))    DEALLOCATE(diag_theta_relax)

IF (ALLOCATED(diag_u_data))         DEALLOCATE(diag_u_data)
IF (ALLOCATED(diag_u_model))        DEALLOCATE(diag_u_model)
IF (ALLOCATED(diag_u_ntend))        DEALLOCATE(diag_u_ntend)
IF (ALLOCATED(diag_u_mtend))        DEALLOCATE(diag_u_mtend)
IF (ALLOCATED(diag_u_relax))        DEALLOCATE(diag_u_relax)

IF (ALLOCATED(diag_v_data))         DEALLOCATE(diag_v_data)
IF (ALLOCATED(diag_v_model))        DEALLOCATE(diag_v_model)
IF (ALLOCATED(diag_v_ntend))        DEALLOCATE(diag_v_ntend)
IF (ALLOCATED(diag_v_mtend))        DEALLOCATE(diag_v_mtend)
IF (ALLOCATED(diag_v_relax))        DEALLOCATE(diag_v_relax)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nudging_dealloc_data

!##########################################################################

END SUBROUTINE nudging_main1
END MODULE nudging_main1_mod
