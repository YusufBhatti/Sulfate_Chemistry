! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Program for taking fieldsfile/pp output and calculating convective info

PROGRAM crmstyle_coarse_grid

!-----------------------------------------------------------------------
! Description:
!   A research utility for processing high resolution (convection permitting)
! simulation output, averaging to create coarser resolution means and calculate
! quantities like buoyant cloudy updraughts etc typically output from CRM
! (cloud resolving model) simulations.
!
! Makes use of UM code from various parts of the model
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 95
!   Software Standards: UMDP3 standards.
!-----------------------------------------------------------------------

! Modules only used by this program

USE crmstyle_grid_info_mod, ONLY:                                 &
  nprocs, nproc_x, nproc_y, thread_level_set,                      &
  full_new_x,  full_new_y,                                         &
  local_row_len, local_rows, local_new_x, local_new_y,             &
  new_res_x, new_res_y, mask_full_x, mask_full_y

USE crmstyle_pp_data_mod, ONLY:                                   &
  bzy,bzx,bdy,bdx, new_bzy, new_bzx, new_bdy, new_bdx, origin_x1,  &
  origin_y1,pseudo_lat,pseudo_lon,                                 &
  mask_bzy, mask_bzx, mask_bdy, mask_bdx

USE crmstyle_cntl_mod, ONLY:                                          &
  nx_start ,ny_start ,nres ,new_res ,ntimes, in_cols, in_rows         &
 ,num_x ,num_y, model_levels, mlevs, height_gen_method, iprint        &
 ,l_bcu, l_wg1, l_acc, l_acu, l_ppd, l_nbd, l_bcu_mask, l_whole_grid  &
 ,l_no_orog, l_all_sea, l_class_col, l_cape, l_pcape, l_bcw, l_sect30 &
 ,l_plume_size, l_pp, num_ff, l_qgraup, run_coarse_grid, bl_levels    &
 ,run_date_info, start_date, modtstep_date, datastep_date             &
 ,nbins_diam                                                          &
 ,adx, ady, azy, azx, apseudo_lat, apseudo_lon, lev_list, num_want    &
 ,timestep                                                            &
! Subroutines to print namelists
 ,print_nlist_run_coarse_grid ,print_nlist_run_date_info

USE crmwork_arrays_mod, ONLY:                                            &
   orog_full, prec_full, r_theta_levels_local, r_rho_levels_local,       &
   r_theta_sea2, h_theta_sea, h_rho_sea, dz_theta_sea, mask,             &
  uv_km_level, th_km_level, th_km_level_full, xcoslat_full,              &
  uv_weight, th_weight, th_weight_full, r_theta_sea, r_rho_sea

USE hires_data_mod, ONLY: orog, landsea

USE crmstyle_filenames_mod, ONLY:                                    &
  ProgName, input_file, orogfile, landseafile, lev_nl_file, pp_file  &
 ,all_file, acc_file, acu_file, bcu_file, wg1_file, ppd_file         &
 ,nbd_file, nid_file, max_ff

! Modules from UM code
! ------------------------
! UM constants

USE conversions_mod,     ONLY: pi_over_180

! Planet constants

USE planet_constants_mod, ONLY:                                          &
    i_planet, ip_earth, set_planet_constants, planet_radius

! Vertical levels namelist

USE vertnamelist_mod, ONLY:                                              &
    first_constant_r_rho_level, z_top_of_model, eta_theta, eta_rho

USE level_heights_mod, ONLY:                                             &
    r_theta_levels, r_rho_levels

USE filenamelength_mod, ONLY:                                           &
  filenamelength

USE file_manager, ONLY: assign_file_unit, release_file_unit

! Info for parallel processing

USE UM_ParCore, ONLY: mype, nproc_max
USE UM_ParVars, ONLY: gc_all_proc_group, change_decomposition

USE mpl, ONLY:                                                              &
    mpl_thread_multiple,                                                     &
    mpl_thread_serialized,                                                   &
    mpl_thread_funneled,                                                     &
    mpl_thread_single

! For Open MP
!$ USE omp_lib           ! Note OpenMP sentinel

! UM IO routines etc

USE io
USE IO_Mod

USE Err_Mod, ONLY:                                                      &
  StatusOK, StatusFatal, StatusWarning

! subroutine and exec number for UM application
USE UM_Config,               ONLY: appInit,appTerminate
USE application_description, ONLY: exe_crmstyle_coarse_grid

USE hostname_mod,            ONLY: get_hostname

! Subroutines used

USE read_lev_info_mod
USE crmstyle_decompose_grid_mod
USE alloc_hires_data_mod
USE crmstyle_get_env_mod
USE crmstyle_read_landsea_mod
USE alloc_crmwork_arrays_mod
USE alloc_sample_arrays_mod
USE crmstyle_cal_weights_mod
USE calc_heights_sea_mod
USE crmstyle_read_orog_mod
USE crmstyle_read_pp_input_mod
USE crmstyle_read_ff_input_mod
USE crmstyle_all_means_mod
USE crmstyle_sample_mod
USE crmstyle_write_ff_out_mod
USE crmstyle_write_ff_mask_mod
USE update_time_mod
USE crmstyle_read_all_ffhdr_mod
USE crmstyle_class_col_mod
USE crmstyle_pcape_mod
USE cape_cin_from_mean_mod
USE crmstyle_hstat_balance_mod

USE ereport_mod, ONLY: ereport, ereport_finalise

! Referenced by umPrintMgr now
!USE printstatus_mod, ONLY:                                              &
!  PrintStatus
USE umPrintMgr                        ! Required for writing output
! Required for Dr Hook
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
 
IMPLICIT NONE

!-----------------------------------------------------------------------
! Declarations
!-----------------------------------------------------------------------

! Info for grid decomposition requried by UM routines
INTEGER, PARAMETER :: model_type = 2        ! indicates Limited area
INTEGER, PARAMETER :: extended_halo_ew = 1
INTEGER, PARAMETER :: extended_halo_ns = 1

! Needed for GCOM
INTEGER, PARAMETER :: gc_alltoall_version = 2
INTEGER, PARAMETER :: gc_alltoall_multi   = 2

! Unit number for namelist input
INTEGER  :: name_unit

INTEGER ::            &
  i,j,k, it,          & ! loop counters
  k1, k2, kloop         ! printing loops

INTEGER ::            &
  coarse_deco         & ! decomposition for coarse grid
 ,hires_deco          & ! decomposition for high resolution input grid
 ,mask_deco           & ! decomposition for high resolution output grid - mask
 ,icode               & ! error return code
 ,errorstatus         & ! error code
 ,length              & ! length of returned string
 ,imultiple_out         ! set to 1 for output in different files

INTEGER ::            &
  date_required(6)    & ! next date required
 ,date_requiredp1(6)  & ! next date required plus one time step
 ,date_current(6)       ! copy of current date

CHARACTER(LEN=10) :: c_thread            ! requested threading type
CHARACTER(LEN=10) :: thread_level_setc   ! actual threading type
CHARACTER(LEN=filenamelength)             ::   &
  stdout_filename = "dummy stdout"             & ! file to WRITE stdout to
 ,stdout_basename = "dummy stdout"               ! base of filename

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CRMSTYLE_COARSE_GRID'

TYPE(UM_Header_type) ::   &
  mask_hdr                & ! UM Headers:  land sea
 ,ff_hdr(max_ff)            ! UM Headers:  fieldsfiles

REAL, ALLOCATABLE :: &
  latitude(:)               ! If no landsea ancillary

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
! Program running on MPP - some things done on all PE/nodes
!------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! Get information from environment variables
! These hold information on processors and I/O filenames
!-----------------------------------------------------------------------

WRITE(umMessage,'(A)') ' Call to get environment variables '
CALL umPrint(umMessage,src=RoutineName)

! Set type of threading
! Cannot call get_env_var before GCOM initialisation.
CALL GET_ENVIRONMENT_VARIABLE('THREAD_LEVEL',c_thread,length,errorstatus)
IF (ErrorStatus  /=  0 .OR. length == 0) THEN
  WRITE(umMessage,'(A,A)') 'Warning: Environment variable THREAD_LEVEL has ',&
      'not been set.'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)') 'Setting thread_level to multiple'
  CALL umPrint(umMessage,src=RoutineName)
  thread_level_setc = 'MULTIPLE'
ELSE
  READ(c_thread,'(A10)') thread_level_setc
END IF
SELECT CASE (thread_level_setc)
CASE ('MULTIPLE')
  thread_level_set = mpl_thread_multiple
CASE ('SERIALIZED')
  thread_level_set = mpl_thread_serialized
CASE ('FUNNELED')
  thread_level_set = mpl_thread_funneled
CASE ('SINGLE')
  thread_level_set = mpl_thread_single
CASE DEFAULT
  WRITE(umMessage,'(A,A,A)') 'Warning: Thread level ', thread_level_setc,  &
      ' not recognised, setting to MULTIPLE.'
  CALL umPrint(umMessage,src=RoutineName)
  thread_level_set = mpl_thread_multiple
END SELECT

!-----------------------------------------------------------------------
! Initialise GCOM need this to start MPP - copied what is being done in
! UM_shell & reconfiguration
!-----------------------------------------------------------------------

CALL gc_init_thread(mype, nprocs, thread_level_set)

!------------------------------------------------------------------------------
! Output separate for each PE/node - Copied from UM_shell
!------------------------------------------------------------------------------
IF (nprocs > 1) THEN
  imultiple_out = 1
ELSE
  imultiple_out = 0
END IF

CALL crmstyle_get_env(1)

nproc_max = nprocs

IF (mype == 0) THEN
  WRITE(umMessage,'(3(A,I3))') 'Processor decomposition:  nproc_x: ',nproc_x, &
                               '  nproc_y: ',nproc_y,'  Total nprocs: ',nprocs
END IF

CALL umPrint(umMessage,src=RoutineName)

! Need a call to appInit  - Now sets up output PEs
CALL appInit(exe_crmstyle_coarse_grid)

! Write out hostname
WRITE(umMessage,'(A,A,A,I0)') 'Host is ',TRIM(get_hostname()),' on PE ',mype
CALL umPrint(umMessage, src=ProgName)

! Set GCOM to use the alternative version of RALLTOALLE
! throughout the run  (copied from UM_SHELL and reconfiguration
CALL Gc_Setopt(gc_alltoall_version, gc_alltoall_multi, Errorstatus)

! Taken from UM_SHELL
! Only want OpenMP section executing if OpenMP is compiled in,
! so protect by sentinal
!$OMP PARALLEL DEFAULT(NONE)
!$OMP MASTER
!$  WRITE(umMessage,'(A,I2,A)') 'I am running with ',                       &
!$     omp_get_num_threads(),' thread(s).'
!$  CALL umPrint(umMessage,src=RoutineName)
!$  WRITE(umMessage,'(A,I6)') 'OpenMP Specification: ',openmp_version
!$  CALL umPrint(umMessage,src=RoutineName)
!$OMP END MASTER
!$OMP END PARALLEL

#if defined(IBM_XL_FORTRAN)
! On IBM force buffering of Fortran I/O for initialisation
CALL setrteopts('buffering=enable')
#endif


WRITE(umMessage,'(A,I3)') ' I am on mype ',mype
CALL umPrint(umMessage,src=RoutineName)
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Read in control namelist - required area of full grid and resolutions
!-----------------------------------------------------------------------
! First need to open unit

CALL assign_file_unit(input_file, name_unit, handler="fortran")

OPEN( UNIT   = name_unit, ACCESS = "SEQUENTIAL", ACTION = "READ",        &
      FILE   = TRIM(input_file), FORM   = "FORMATTED",                   &
      IOSTAT = ErrorStatus,  STATUS = "OLD" )

READ(name_unit,run_coarse_grid)

CALL print_nlist_run_coarse_grid()
WRITE(umMessage,'(A60)')                                           &
       '============================================================'
CALL umPrint(umMessage,src=RoutineName)

! Calculate timestep of increments from namelist input 
! Assumes less than an hour
timestep = REAL(modtstep_date(5))*60. + REAL(modtstep_date(6)) 



full_new_x = num_x
full_new_y = num_y

! only one new output resolution coded at present

new_res_x = new_res(1)
new_res_y = new_res(1)

IF (iprint >= 1) THEN
  PrintStatus = 4  ! All extra print diagnostics (used in UM IO routines)
ELSE
  PrintStatus = 2  ! normal UM output from routines
END IF

! Need to alter lev_list to reflect actual number of levels

DO i=1,num_want
  IF (lev_list(i) /= 1) THEN
    lev_list(i) = mlevs
  END IF
END DO

!-----------------------------------------------------------------------
! Read in namelist - for date and timestep info
!-----------------------------------------------------------------------

READ(name_unit,run_date_info)

CALL print_nlist_run_date_info()
WRITE(umMessage,'(A60)')                                             &
        '============================================================'
CALL umPrint(umMessage,src=RoutineName)

! close namelist input file
CLOSE(name_unit)

CALL release_file_unit(name_unit, handler="fortran")

!-----------------------------------------------------------------------
! Second call to get the rest of environment variables now read namelist
! Needs to be before model namelist read as gets filename for this
!-----------------------------------------------------------------------

CALL crmstyle_get_env(2)

!-----------------------------------------------------------------------
! Read in model level namelist
!-----------------------------------------------------------------------

CALL read_lev_info(lev_nl_file)

!-----------------------------------------------------------------------
! Planet constants: set planet as Earth
!-----------------------------------------------------------------------
i_planet = ip_earth
CALL set_planet_constants()

!-----------------------------------------------------------------------
! Now call decompose grid for MPP - input grid from fieldsfiles
!-----------------------------------------------------------------------
! Portion of new grid on this processor
! Note expect number of processors chosen to give sensible answers here

local_new_x = full_new_x / nproc_x
local_new_y = full_new_y / nproc_y

! proportion of old grid required by this processor

local_row_len = new_res_x * local_new_x
local_rows    = new_res_y * local_new_y

hires_deco = 1   ! first decomposition

CALL crmstyle_decompose_grid(hires_deco, in_cols, in_rows, mlevs,  &
                             model_type,                           &
                             nproc_x, nproc_y, nprocs,             &
                             extended_halo_ew, extended_halo_ns,   &
                             nx_start, ny_start,                   &
                             local_row_len, local_rows)

! Set up the atmosphere decomposition in PARVARS (Needed for information to
! be correctly setup for IO calls)

icode = 0
CALL change_decomposition(hires_deco,icode)

WRITE(umMessage,'(A60)')                                      &
     '============================================================'
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,5I6)') ' I am MYPE : ',mype,local_row_len, local_rows, &
             local_new_x, local_new_y
CALL umPrint(umMessage,src=RoutineName)

!-----------------------------------------------------------------------
! Now call decompose grid for coarse output grid
!-----------------------------------------------------------------------

coarse_deco = 2     ! 2nd decomposition

CALL crmstyle_decompose_grid(coarse_deco, num_x, num_y, mlevs,      &
                           model_type,                              &
                           nproc_x, nproc_y, nprocs,                &
                           extended_halo_ew, extended_halo_ns,      &
                           1,  1,                                   &
                           local_new_x, local_new_y)

!-----------------------------------------------------------------------
! Now call decompose grid for high resolution output grid for Bcu mask
! As input grid but without outer rim not processed by this program.
!-----------------------------------------------------------------------
IF (l_bcu_mask .OR. l_class_col) THEN
  mask_deco = 3     ! 3rd decomposition

  mask_full_x = num_x*new_res_x
  mask_full_y = num_y*new_res_y

  CALL crmstyle_decompose_grid(mask_deco, mask_full_x, mask_full_y, mlevs,  &
                           model_type,                                      &
                           nproc_x, nproc_y, nprocs,                        &
                           extended_halo_ew, extended_halo_ns,              &
                           1,  1,                                           &
                           local_row_len, local_rows)
END IF
!-----------------------------------------------------------------------
! Allocate full arrays required
!-----------------------------------------------------------------------
! Some of these are only needed on mype = 0

CALL alloc_crmwork_arrays(in_cols, in_rows, model_levels)

!-----------------------------------------------------------------------
! Allocate arrays for holding input data on reduced region to be processed
! on this PE
!-----------------------------------------------------------------------

CALL alloc_hires_data( local_row_len, local_rows, mlevs)

IF (l_all_sea) THEN
  ! Problem if idealised model with no land sea grid
  ! Need to provide grid info

  DO j= 1,local_rows
    DO i= 1, local_row_len
      landsea(i,j) =  0.0
    END DO
  END DO

  ! Set grid info from additional namelist info
  bdx = adx
  bdy = ady
  bzy = azy
  bzx = azx
  pseudo_lat = apseudo_lat
  pseudo_lon = apseudo_lon

  ALLOCATE(latitude(in_rows))
  DO j=1,in_rows
    ! latitude of grid points - full grid
    latitude(j) = bzy + REAL(j)*bdy
    DO i= 1,in_cols
      xcoslat_full(i,j) = COS(latitude(j)*pi_over_180)
    END DO
  END DO
  DEALLOCATE(latitude)

ELSE
  !-----------------------------------------------------------------------
  ! Read land sea mask for model grid - get original grid sizes from this
  !-----------------------------------------------------------------------

  CALL crmstyle_read_landsea(in_cols, in_rows, mype, gc_all_proc_group,    &
                             local_row_len,local_rows,                     &
                             bzy,bzx,bdy,bdx,pseudo_lat, pseudo_lon,mask_hdr)

END IF
!-----------------------------------------------------------------------
! Read in orography
!-----------------------------------------------------------------------
IF (l_no_orog) THEN
  ! No orogaphy ancillary as whole grid flat
  DO j= 1,local_rows
    DO i= 1, local_row_len
      orog(i,j) =  0.0
    END DO
  END DO

  DO j= 1,in_rows
    DO i= 1, in_cols
      orog_full(i,j) =  0.0
    END DO
  END DO

ELSE
  CALL crmstyle_read_orog(in_cols, in_rows, mype, gc_all_proc_group,    &
                          local_row_len,local_rows)
END IF
!-----------------------------------------------------------------------
! Work out details of full new output grid
!-----------------------------------------------------------------------

new_bdy = new_res(1)*bdy
new_bdx = new_res(1)*bdx
! first lat of nres(1)
origin_y1 = bzy+ny_start*bdy+ 0.5*REAL(new_res(1)-1)*bdy
origin_x1 = bzx+nx_start*bdx+ 0.5*REAL(new_res(1)-1)*bdx
new_bzy = origin_y1 -new_bdy
new_bzx = origin_x1 -new_bdx

WRITE(umMessage,'(A,6F16.10)') ' New grid details ',new_bzy,new_bdy,&
         new_bzx,new_bdx,pseudo_lat, pseudo_lon
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') '===================================================='
CALL umPrint(umMessage,src=RoutineName)

WRITE(umMessage,'(A,I8)') ' full_new_x :',full_new_x
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' full_new_y :',full_new_y
CALL umPrint(umMessage,src=RoutineName)

WRITE(umMessage,'(A,I8)') ' new_res_x  :',new_res_x
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' new_res_y  :',new_res_y
CALL umPrint(umMessage,src=RoutineName)

WRITE(umMessage,'(A,I8)') ' local_row_len :',local_row_len
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' local_rows    :',local_rows
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' local_new_x :',local_new_x
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' local_new_y :',local_new_y
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' num_x :',num_x
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' num_y :',num_y
CALL umPrint(umMessage,src=RoutineName)

WRITE(umMessage,'(A,I8)') ' in_cols :',in_cols
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,I8)') ' in_rows :',in_rows
CALL umPrint(umMessage,src=RoutineName)

! Details for BCu mask grid or column classification
IF (l_bcu_mask .OR. l_class_col) THEN
  mask_bdy = bdy              ! as input grid
  mask_bdx = bdx
  mask_bzy = bzy+(ny_start-1)*bdy     ! start of grid
  mask_bzx = bzx+(nx_start-1)*bdx
END IF
!-----------------------------------------------------------------------
! Calculate heights of model levels - using central subroutine
! Stores levels in  level_heights_Mod as r_theta_levels & r_rho_levels
!-----------------------------------------------------------------------
ALLOCATE (r_theta_sea(model_levels) )
ALLOCATE (r_rho_sea(model_levels) )

CALL calc_heights_sea(height_gen_method,model_levels, bl_levels,           &
                      local_rows, local_row_len,                           &
                      in_rows, in_cols )

IF (iprint >= 1) THEN
  i = new_res(1)/2
  WRITE(umMessage,'(A,I4)') ' r_theta_levels, r_rho_levels for ',i
  CALL umPrint(umMessage,src=RoutineName)
  DO k=1,model_levels
    WRITE(umMessage,'(2F20.11)') r_theta_levels_local(i,i,k),    &
                          r_rho_levels_local(i,i,k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

!-----------------------------------------------------------------------
! Calculate weight for interpolating from input fields on model
! theta and UV levels to a set of fixed height levels
!-----------------------------------------------------------------------

ALLOCATE (h_theta_sea(mlevs) )
ALLOCATE (h_rho_sea(mlevs) )
ALLOCATE (dz_theta_sea(mlevs) )
ALLOCATE (r_theta_sea2(mlevs) )

! Required heights
DO k=1,mlevs
  h_theta_sea(k) = r_theta_sea(k)-planet_radius
  r_theta_sea2(k)= r_theta_sea(k)
  h_rho_sea(k) = r_rho_sea(k)-planet_radius
END DO

! Thickness of model levels  - used for LWP * IWP only
dz_theta_sea(1) = r_rho_sea(2) - planet_radius
DO k=2,mlevs-1
  dz_theta_sea(k) = r_rho_sea(k+1) - r_rho_sea(k)
END DO
dz_theta_sea(mlevs) = r_theta_sea(mlevs) - r_rho_sea(mlevs)

WRITE(umMessage,'(A)') ' heights required '
CALL umPrint(umMessage,src=RoutineName)

! Problem as number of levels unknown and want all levels. umPrint will
! only print the first 1024 characters of any buffer

kloop=(mlevs+19)/20

DO k=1,kloop
  k1=(k-1)*20+1
  k2=k*20
  IF (k2 > mlevs) THEN
    k2=mlevs
  END IF
  WRITE(umMessage,'(20F10.1)') (h_theta_sea(j),j=k1,k2)
  CALL umPrint(umMessage,src=RoutineName)
END DO

! Work out weights for uw each PE

CALL crmstyle_cal_weights(local_row_len,local_rows,mlevs,                   &
                          r_rho_levels_local,orog,                          &
                          r_theta_sea2, uv_km_level,uv_weight)

! Work out theta weights for whole grid only really need on PE/node 0 but
! currently done on all.

CALL crmstyle_cal_weights(in_cols,in_rows,mlevs,r_theta_levels(1,1,1),      &
                          orog_full,r_theta_sea, th_km_level_full,          &
                          th_weight_full)

! Work out theta weights for local grid

CALL crmstyle_cal_weights(local_row_len,local_rows,mlevs,                   &
                          r_theta_levels_local(1,1,1),                      &
                          orog, r_theta_sea, th_km_level,th_weight)


! Now finished using full model level info so free space

DEALLOCATE (r_rho_sea)
DEALLOCATE (r_theta_sea)
DEALLOCATE (orog_full)
DEALLOCATE (r_rho_levels)
DEALLOCATE (r_theta_levels)
DEALLOCATE (r_theta_levels_local)
DEALLOCATE (r_rho_levels_local)

! Mask for just this area

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                  &
!$OMP& SHARED(mlevs, local_rows, local_row_len, mask, th_km_level)
DO k=1,mlevs
  DO j=1,local_rows
    DO i=1,local_row_len
      IF (th_km_level(i,j,k) >= 0) THEN
        mask(i,j,k) = .TRUE.
      ELSE
        mask(i,j,k) = .FALSE.
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (iprint >= 1) THEN
  WRITE(umMessage,'(A)') ' All theta and UV weights have been calculated'
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)') ' mask th_km_level th_weight, uv_km_level, uv_weight '
  CALL umPrint(umMessage,src=RoutineName)
  i = new_res(1)/2
  DO k=1,mlevs
    WRITE(umMessage,'(L2,2(I5,F13.10))') mask(i,i,k),th_km_level(i,i,k),     &
               th_weight(i,i,k), uv_km_level(i,i,k),uv_weight(i,i,k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

! Arrays for CRM sampling & meaning

CALL alloc_sample_arrays(local_new_x,local_new_y,mlevs,local_row_len,        &
                         local_rows, nbins_diam)

IF (.NOT. l_pp) THEN

  ! Open and Read all input fieldsfile headers
  CALL crmstyle_read_all_ffhdr(num_ff,ff_hdr)

END IF

!-----------------------------------------------------------------------
! Main loop over time
! Read in fieldfiles/ pp files  looping over times
! Note run dkek converted output to pp from fieldsfile so need different
! input
!-----------------------------------------------------------------------

date_required(:) = start_date(:)

CALL update_time(start_date, modtstep_date, date_requiredp1)

WRITE(umMessage,'(A,6i5)') ' start date        ',date_required
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A,6i5)') ' start date plus 1 ',date_requiredp1
CALL umPrint(umMessage,src=RoutineName)

DO it = 1, ntimes

  ! Column classification need to reset surface precipitation to zero
  IF (l_class_col .AND. mype == 0) THEN 
    DO j=1,in_rows
      DO i=1,in_cols
        prec_full(i,j) = 0.0
      END DO
    END DO
  END IF 

  ! Read in data from pp file /fieldsfile ?
  IF (l_pp) THEN
    ! pp files contains data for just the one required time
    ! NO date time checking
    ! Enables use of data from old runs eg EMBRACE

    CALL crmstyle_read_pp_input(pp_file(it),icode)
    IF (icode /= 0) THEN
      CALL EReport("CRMstyle_coarse_grid", icode,                          &
                       "Failed to read data from pp file" )
    END IF
  ELSE

    ! Fieldsfiles may contain more than one date and time
    ! Need date time checking. Expects fields like precip, sensible heat
    ! tendencies offset by one timestep from model prognostics.

    CALL crmstyle_read_ff_input(date_required,date_requiredp1,num_ff,ff_hdr)

  END IF

  WRITE(umMessage,'(A)')    '=========================='
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I4)') ' read in data time ',it
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A)')    '=========================='
  CALL umPrint(umMessage,src=RoutineName)

  ! Column classification if required
  IF (l_class_col) THEN 
    CALL crmstyle_class_col( ) 
  END IF 
  ! Create mean fields on new grid

  CALL crmstyle_all_means(mask)
  IF (iprint >= 1) THEN
    WRITE(umMessage,'(A)') ' After meaning'
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  ! From means profiles derive a hydrostatic pressure
  ! This will not be the same as the actual pressure on theta levels as
  ! the model is non-hydrostatic

  CALL crmstyle_hstat_balance(mlevs)

  ! Create special sampling of fluxes etc

  CALL crmstyle_sample(mask)
  IF (iprint >= 1) THEN
    WRITE(umMessage,'(A)')  ' CRM sampling'
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  ! Calculate z LCL plus undilute CAPE & CIN (zlcl need by PCAPE)
  IF ((l_cape .OR. l_pcape) .AND. l_sect30 .AND. (l_bcu .OR. l_bcw)) THEN 
    CALL cape_cin_from_mean()
  END IF

  ! Calculate dilute PCAPE plus rate of change of PCAPE
  IF (l_pcape) THEN 
    CALL crmstyle_pcape()
  END IF


  ! change decomposition for output 
  CALL change_decomposition(coarse_deco,icode)

  ! Write all fields created out
  CALL crmstyle_write_ff_out(it,ntimes)

  ! Write out fields still on full grid but just for area required.
  IF (l_bcu_mask .OR. l_class_col) THEN    ! BCu mask output, column class
    ! change decomposition for  this grid
    CALL change_decomposition(mask_deco,icode)

    CALL crmstyle_write_ff_mask(it,ntimes)
  END IF

  ! change decomposition back for next time
  CALL change_decomposition(hires_deco,icode)

  ! update time required ready for next iteration
  date_current(:) = date_required(:)

  CALL update_time(date_current,  datastep_date, date_required  )
  CALL update_time(date_required, modtstep_date, date_requiredp1)

  WRITE(umMessage,'(A,6I6)') ' new date        ',date_required
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,6I6)') ' new date plus 1 ',date_requiredp1
  CALL umPrint(umMessage,src=RoutineName)

END DO   ! end time loop

CALL appTerminate()

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! end  - close every thing?
!-----------------------------------------------------------------------

END PROGRAM crmstyle_coarse_grid
