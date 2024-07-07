! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: MULTI_SPATIAL --------------------------------------------
!
! Purpose: Control routine for spatial processing when extracting a
!          timeseries within STASH.  Calls SPATIAL to extract global
!          mean samples from subdomains pointed to by the mother
!          STASHlist record using weighting and masking codes from
!          the STASHlist record within each subdomain.  The list of
!          subdomains is held as part of the STASH control file.
!          They may be in terms of gridpoints, or latitude/longitude
!          relative to the grid coordinates.  All timeseries samples
!          at a given step are appended to the output field.
!          On the first timestep it fills the entire output
!          vector to missing data (apart from values for the
!          first timestep and computes extra data for the timeseries
!          This prevents time meaning routines failing due to 
!          uninitialised data.  However as a result of this the output 
!          vector length will change from timestep to timestep.
!
! External documentation:
!   Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                system (STASH)
!
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE multi_spatial(fieldin,vx,vy,vz,gr,st_grid,             &
                         fld_type,halo_type,                      &
                         halo_x,halo_y,                           &
                         lcyclic,                                 &
     lmasswt,horiz_size,num_vert_levels,                          &
     no_of_levels_masswt,                                         &
     p,pstar,                                                     &
     cos_v_latitude,cos_theta_latitude,land,sea,                  &
     row_length,rows,n_rows,model_levels,                         &
     fieldout,lenout,amdi,                                        &
     control,control_size,                                        &
     stash_series,series_entry_size,no_records,                   &
     num_stash_levels,index_lev,level_list,start_ts,extraw,       &
     n_rows_out,n_cols_out,                                       &
     real_hd,len_realhd,int_hd,len_inthd,                         &
     icode,errmssg)

USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod, ONLY: gc_none
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE atm_fields_bounds_mod, ONLY: vdims_s, tdims, tdims_s, pdims,  &
                                 pdims_s
USE sterr_mod, ONLY: st_not_supported, st_bad_array_param,        &
                     st_unknown,       unknown_processing
USE stparam_mod, ONLY: st_proc_no_code, st_time_series_code,      &
    st_time_series_mean, st_gridpoint_code, block_size,           &
    st_output_length, st_input_bottom, st_input_top,              &
    st_weight_code, st_north_code, st_south_code, st_west_code,   &
    st_east_code, st_period_code, st_freq_code, extract_base,     &
    series_size, series_grid_type, series_grid_code,              &
    series_proc_code, series_north, series_south, series_west,    &
    series_east, series_list_start, series_list_end,              &
    merid_mean_base, vert_mean_base, field_mean_base, zonal_mean_base

IMPLICIT NONE


INTEGER :: vx,vy,vz          ! IN size of fieldin
INTEGER :: gr                ! IN ppxref gridtype code
INTEGER :: st_grid           ! IN STASH internal gridtype code
INTEGER ::                                                        &
 fld_type                                                         &
                         ! IN field type (u/v/p) of input field
,halo_type                                                        &
                         ! IN halo type of input field
,halo_x                                                           &
                         ! IN halo size East-West
,halo_y                  ! IN halo size North-South
LOGICAL :: lcyclic           ! IN TRUE if cyclic EW BCs
LOGICAL :: lmasswt           ! IN TRUE if level-dep mass-wts OK
INTEGER :: row_length,rows,n_rows,model_levels ! IN size parameters
REAL :: fieldin(vx*vy*vz)    ! IN data field
INTEGER :: horiz_size        ! OUT no. of points in horizontal slice
INTEGER :: num_vert_levels   ! OUT no of horizontal slices.
INTEGER :: no_of_levels_masswt ! IN: levels for mass weighting array
INTEGER :: num_stash_levels  ! IN size of vertical levels list.
INTEGER :: index_lev(num_stash_levels) ! IN offsets for each horiz fi
INTEGER :: level_list(num_stash_levels) ! IN level for each horiz. fi

REAL ::                                                           &
 p(pdims_s%i_start:pdims_s%i_end,                                 &
   pdims_s%j_start:pdims_s%j_end,                                 &
   pdims_s%k_start:pdims_s%k_end)                                 &
                                        ! IN pressure (rho levels)
,pstar(pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end)      &
                                        ! IN surface pressure
,cos_v_latitude(vdims_s%i_start:vdims_s%i_end,                    &
                vdims_s%j_start:vdims_s%j_end)                    &
                                        ! IN v-grid area fn
,cos_theta_latitude(tdims_s%i_start:tdims_s%i_end,                &
                    tdims_s%j_start:tdims_s%j_end)
                                        ! IN T-grid area fn

LOGICAL :: land(tdims%i_start:tdims%i_end,                        &
                tdims%j_start:tdims%j_end)          ! IN land mask
LOGICAL :: sea(tdims%i_start:tdims%i_end,                         &
               tdims%j_start:tdims%j_end)           ! IN sea mask
INTEGER :: lenout               ! IN max size of output field
REAL :: fieldout(lenout)        ! OUT output field
REAL :: amdi                    ! IN missing data indicator
INTEGER :: len_realhd           ! IN size of real header
INTEGER :: len_inthd            ! IN size of integer header
REAL :: real_hd(len_realhd)     ! IN real header
INTEGER :: int_hd(len_inthd)    ! IN integer header
INTEGER :: control_size         ! IN size of control array
INTEGER :: control(control_size)! IN control array (mostly not used)
INTEGER :: series_entry_size    ! IN no of entries in each record.
INTEGER :: no_records           ! IN no of records to process.
INTEGER :: extraw               ! OUT no of words required by extra d
INTEGER :: n_rows_out,n_cols_out! OUT data-set size and extent
LOGICAL :: start_ts             ! IN true if first time-series timest
INTEGER :: stash_series(series_entry_size,no_records) ! IN
! IN control data for calls to spatial
INTEGER :: icode                       ! OUT error code
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: errmssg
                ! OUT error message
! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER :: fake_record(control_size) ! fake_record for SPATIAL call
INTEGER :: fake_record_extradata(control_size) ! fake_record to form extra data
INTEGER :: stash_list_start_extradata ! start address in index_levs for levels
INTEGER :: stash_list_end_extradata   ! end address in index_levs for levels
INTEGER :: pp_ptr ! ptr to pp_field for where output from spatial goe
INTEGER :: i                         ! loop variable
INTEGER :: data_size                 ! size of data.
INTEGER :: this_index_lev   ! index of input level
INTEGER :: stash_list_start ! the start address in index_levs for lev
INTEGER :: stash_list_end ! the end address in index_levs for levels
INTEGER :: what_process ! what kind of processing is requested
INTEGER :: what_mask ! what mask is wanted.
INTEGER :: extra_start ! start address in fieldout for extra data
REAL :: bdx,bzx,bdy,bzy ! grid descriptors

INTEGER ::                                                        &
  proc_start_x                                                    &
                  ! processor co-ords of processor at origin
, proc_start_y                                                    &
                  ! (SW corner) of subarea
, start_pe                                                        &
                  ! processor id of this processor
, dummy1,dummy2                                                   &
                  ! unused return variables from
!                         GLOBAL_TO_LOCAL_RC
      , info            ! GCOM return variable

REAL :: global_mean ! global mean returned by call to SPATIAL
REAL :: fieldout_col(num_stash_levels) ! single column values returned 
                                       ! by call to SPATIAL

LOGICAL :: flag_column  ! Indicates column output for a single ts point

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MULTI_SPATIAL'

! ----------------------------------------------------------------------
!  1. Error checking
!
!  Check we are in fact doing a time series. Error if not.
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF ((control(st_proc_no_code) /= st_time_series_code) .AND.        &
   (control(st_proc_no_code) /= st_time_series_mean)) THEN
  icode=st_unknown
  WRITE(errmssg,'(3X,A,I5)') & 
  'MULTI_SP : unexpected call to extract timeseries',control(st_proc_no_code)
  GO TO 9999            ! jump to return
END IF
! ----------------------------------------------------------------------
!  2. Workout what kind of processing we are doing and what mask is used
!
what_mask=MOD(control(st_gridpoint_code),block_size)
what_process=control(st_gridpoint_code)-what_mask
! Set flag_column to generate extra data for time series at 1stt step
IF ( what_process == extract_base ) THEN
  flag_column = .TRUE.
ELSE
  flag_column = .FALSE.
END IF

extraw=0
!  note that the first word in fieldout is assumed to be
!  the word where the first spatial domain mean for this timeseries
!  will be stored
!
!  Next we compute the grid discriptors -- as used by extra data
IF (start_ts) THEN ! compute grid descrip
  extraw=6*(no_records+1)
  extra_start=control(st_output_length)-extraw+1
  ! DEPENDS ON: stash_comp_grid
  CALL stash_comp_grid(bzx,bzy,bdx,bdy,0,st_grid,                &
    1,1,real_hd,len_realhd,int_hd,len_inthd,extract_base+1,      &
    icode,errmssg)
END IF
!  3. Set up pseudo STASH record to be passed to SPATIAL on each call
!     to extract a sample from the input field.
!
fake_record(st_input_bottom)=control(st_input_bottom)
fake_record(st_input_top)=control(st_input_top)
fake_record(st_weight_code)=control(st_weight_code)
! doing a field mean on this sub-domain with mask specified by control.
pp_ptr=1
! ----------------------------------------------------------------------
!  4. Loop over samples and extract global mean within subdomain for
!     each, appending to output field
!
DO i=1,no_records ! loop over sub domains
  data_size=stash_series(series_size,i)
  !  4.1 Do preliminary verifications on stash_series
  !  4.1.1 Gridtype code
  IF (stash_series(series_grid_type,i) /= series_grid_code) THEN
    !
    ! NB: Latitude/longitude range conversion to gridpoint range needs to
    !     be added
    !
    icode=st_not_supported
    errmssg='MULTI_SP : only support grid point processing'
    GO TO 9999
  END IF
  !  ---------------------------------------------------------------------
  !  5. Set up the fake record domain info depending on what kind of
  !     "primary" processing is requested.
  !     As far as stash is concerned everything looks like a global mean
  !     here and it is just a question of setting up the fake record
  !     correctly.
  !
  fake_record(st_gridpoint_code)=                                 &
    (stash_series(series_proc_code,i)/block_size)*block_size      &
      +what_mask

  ! Old style output for each point on each level
  IF (.NOT. flag_column) THEN 
  IF (what_process == extract_base) THEN ! an extract
    fake_record(st_north_code)=stash_series(series_north,i)
    fake_record(st_south_code)=stash_series(series_south,i)
    fake_record(st_west_code)= stash_series(series_west,i)
    fake_record(st_east_code)= stash_series(series_east,i)
    stash_list_start=stash_series(series_list_start,i)
    stash_list_end=stash_series(series_list_end,i)
    fake_record(st_input_bottom)=stash_list_start
    fake_record(st_input_top)=stash_list_end
  ELSE IF (what_process == zonal_mean_base) THEN ! a zonal_mean
    fake_record(st_north_code)=stash_series(series_north,i)
    fake_record(st_south_code)=stash_series(series_south,i)
    fake_record(st_west_code)= control(st_west_code)
    fake_record(st_east_code)= control(st_east_code)
    stash_list_start=stash_series(series_list_start,i)
    stash_list_end=stash_series(series_list_end,i)
    fake_record(st_input_bottom)=stash_list_start
    fake_record(st_input_top)=stash_list_end
  ELSE IF (what_process == merid_mean_base) THEN ! a merid_mean
    fake_record(st_north_code)= control(st_north_code)
    fake_record(st_south_code)= control(st_south_code)
    fake_record(st_east_code)=stash_series(series_east,i)
    fake_record(st_west_code)=stash_series(series_west,i)
    stash_list_start=stash_series(series_list_start,i)
    stash_list_end=stash_series(series_list_end,i)
    fake_record(st_input_bottom)=stash_list_start
    fake_record(st_input_top)=stash_list_end
  ELSE IF (what_process == vert_mean_base) THEN ! a vert_mean
    fake_record(st_north_code)=stash_series(series_north,i)
    fake_record(st_south_code)=stash_series(series_south,i)
    fake_record(st_east_code)=stash_series(series_east,i)
    fake_record(st_west_code)=stash_series(series_west,i)
    stash_list_start=1
    stash_list_end=num_stash_levels
    fake_record(st_input_bottom)=stash_list_start
    fake_record(st_input_top)=stash_list_end
  ELSE IF (what_process == field_mean_base) THEN ! a field_mean
    fake_record(st_north_code)=control(st_north_code)
    fake_record(st_south_code)=control(st_south_code)
    fake_record(st_east_code)=control(st_east_code)
    fake_record(st_west_code)=control(st_west_code)
    stash_list_start=1
    stash_list_end=num_stash_levels
    fake_record(st_input_bottom)=stash_list_start
    fake_record(st_input_top)=stash_list_end
  ELSE ! error code...
    icode=unknown_processing
    WRITE(errmssg,'(A,I5)') &
    'MULTI_SP :  >>>FATAL ERROR << unknown processing option',what_process
    GO TO 9999 ! jump to error return
  END IF

  !  Check record (south > north and west < east)
  IF (fake_record(st_north_code) <                                &
    fake_record(st_south_code)) THEN
    WRITE(errmssg,'(A,2I5,A,I5)') &
    'MULTI_SP : north < south',fake_record(st_north_code), & 
    fake_record(st_south_code),' in record ',i
    icode=st_bad_array_param
    GO TO 9999 ! error exit
  END IF
  IF (fake_record(st_west_code) >                                 &
    fake_record(st_east_code)) THEN
    WRITE(errmssg,'(A,2I5,A,I5)') & 
    'MULTI_SP : west > east',fake_record(st_west_code), & 
    fake_record(st_east_code),'in record ',i
    icode=st_bad_array_param
    GO TO 9999 ! error exit
  END IF

  ! Determine which pe holds the first point of the subdomain, ie the
  ! bottom left hand corner [=SW], that SPATIAL will process.
  ! This is the pe used to gather all points of the subdomain and
  ! calculate their global mean in STGLOM before sending these to pe0
  ! for storage and output.

  ! DEPENDS ON: global_to_local_rc
  CALL global_to_local_rc(gr,halo_type,                           &
    fake_record(st_west_code),fake_record(st_south_code),         &
    proc_start_x, proc_start_y,                                   &
    dummy1,dummy2)

  start_pe=proc_start_x+proc_start_y*nproc_x

  !
  ! NB: At present timeseries samples are global (ie. 3D) means, so
  !     there is no levels loop outside the call to SPATIAL here -
  !     this may be extended at some point to allow multi-level
  !     timeseries sampling inside a levels loop
  !
  !     n_cols_out and n_rows_out are recalculated within SPATIAL but are
  !     now appropriate for an individual timeseries sample, not the whole
  !     field.  They are reset for the whole field after subdomain loop.
  !
  lcyclic=.FALSE.
  this_index_lev = 1     ! no (none) levels loop
  ! DEPENDS ON: spatial
  CALL spatial(fieldin,vx,vy,vz,gr,st_grid,                        &
        fld_type,halo_type,                                        &
        halo_x,halo_y,                                             &
        lcyclic,lmasswt,                                           &
        n_cols_out,n_rows_out,this_index_lev,                      &
        level_list(stash_list_start),                              &
        index_lev(stash_list_start),                               &
        (stash_list_end+1-stash_list_start),                       &
        no_of_levels_masswt,                                       &
        p,pstar,                                                   &
        cos_v_latitude,cos_theta_latitude,land,sea,                &
        row_length,rows,n_rows,blsize(2,fld_type),model_levels,    &
        global_mean,1,                                             &
        fake_record,control_size,amdi,                             &
        icode,errmssg)
  IF (icode /= 0) GO TO 9999 ! got some error so jump to return

  ! Must move the global_mean data to PE 0 which stores all timeseries
  ! data
  ! (NB. This assumes that the output from SPATIAL is just a
  !      single number)

  IF (mype  ==  start_pe) THEN
    info=gc_none
    CALL gc_rsend(100,1,0,info,fieldout(pp_ptr),global_mean)
  END IF

  IF (mype  ==  0) THEN
    info=gc_none
    CALL gc_rrecv(100,1,start_pe,info,                            &
                  fieldout(pp_ptr),global_mean)
  ELSE
    fieldout(pp_ptr)=0.0
  END IF


  pp_ptr=pp_ptr+(n_cols_out*n_rows_out)  ! increment the pp_ptr
  !
  ! NB: n_cols_out and n_rows_out should both be 1 as timeseries samples
  !       are currently designed to be scalar quantities only.
  !
  !  check on n_cols_out and n_rows_out
  IF (n_cols_out /= 1) THEN
    errmssg='MULTI_SP : n_cols_out <> 1'
    icode=st_not_supported
    GO TO 9999
  END IF
  IF (n_rows_out /= 1) THEN
    errmssg='MULTI_SP : n_rows_out <> 1'
    icode=st_not_supported
    GO TO 9999
  END IF

  IF (mype  ==  0) THEN
    IF (start_ts) THEN ! put the descriptive info for this record
      ! DEPENDS ON: extra_make_vector
      CALL extra_make_vector(fake_record,control_size,              &
        i,no_records,fieldout(extra_start),extraw,bzx,bzy,bdx,bdy)
    END IF
  END IF

  !
  ! Extract data by column for each time series geographical point
  !
  ELSE ! flag_column = True
    ! Create fake record for extra data (grid, level descriptions)
    fake_record_extradata(st_north_code)=stash_series(series_north,i)
    fake_record_extradata(st_south_code)=stash_series(series_south,i)
    fake_record_extradata(st_west_code)= stash_series(series_west,i)
    fake_record_extradata(st_east_code)= stash_series(series_east,i)
    stash_list_start_extradata=stash_series(series_list_start,i)
    stash_list_end_extradata=stash_series(series_list_end,i)
    fake_record_extradata(st_input_bottom)=stash_list_start_extradata
    fake_record_extradata(st_input_top)=stash_list_end_extradata

    fieldout(i) = 0.0

    IF (MOD(i,num_stash_levels) == 0) THEN    
      IF (what_process == extract_base) THEN ! an extract == 0

        fake_record(st_north_code)=stash_series(series_north,i)
        fake_record(st_south_code)=stash_series(series_south,i)
        fake_record(st_west_code)= stash_series(series_west,i)
        fake_record(st_east_code)= stash_series(series_east,i)
        stash_list_start=1
        stash_list_end=stash_series(series_list_end,i)
        fake_record(st_input_bottom)=stash_list_start
        fake_record(st_input_top)=stash_list_end

      ELSE ! error code...
        icode=unknown_processing
        WRITE(errmssg,'(A,I5)') &
        'MULTI_SP : >>FATAL ERROR << unknown processing option Column output' &
        ,what_process
         GO TO 9999 ! jump to error return
      END IF

      ! Check record (south > north and west < east)
      IF (fake_record(st_north_code) <                                &
        fake_record(st_south_code)) THEN
        WRITE(errmssg,'(A,2I5,A,I5)') &
        'MULTI_SP : north < south',fake_record(st_north_code), & 
        fake_record(st_south_code),' in record ',i
        icode=st_bad_array_param
        GO TO 9999 ! error exit
      END IF
      IF (fake_record(st_west_code) >                                 &
        fake_record(st_east_code)) THEN
        WRITE(errmssg,'(A,2I5,A,I5)') & 
        'MULTI_SP : west > east',fake_record(st_west_code), & 
        fake_record(st_east_code),'in record ',i
        icode=st_bad_array_param
        GO TO 9999 ! error exit
      END IF

      ! Determine which pe holds the first point of the subdomain
      ! DEPENDS ON: global_to_local_rc
      CALL global_to_local_rc(gr,halo_type,                           &
        fake_record(st_west_code),fake_record(st_south_code),         &
        proc_start_x, proc_start_y,                                   &
        dummy1,dummy2)

      start_pe = proc_start_x + proc_start_y * nproc_x

      lcyclic = .FALSE.
      this_index_lev = 1

      ! DEPENDS ON: spatial
      CALL spatial(fieldin,vx,vy,vz,gr,st_grid,                       &
           fld_type,halo_type,                                        &
           halo_x,halo_y,                                             &
           lcyclic,lmasswt,                                           &
           n_cols_out,n_rows_out,this_index_lev,                      &
           level_list(stash_list_start),                              &
           index_lev(stash_list_start),                               &
           (stash_list_end+1-stash_list_start),                       &
           no_of_levels_masswt,                                       &
           p,pstar,                                                   &
           cos_v_latitude,cos_theta_latitude,land,sea,                &
           row_length,rows,n_rows,blsize(2,fld_type),model_levels,    &
           fieldout_col,num_stash_levels,                             &
           fake_record,control_size,amdi,                             &
           icode,errmssg)

      IF (icode /= 0) GO TO 9999 ! got some error so jump to return

      IF (mype  ==  start_pe) THEN
        info = gc_none
        CALL gc_rsend(100,num_stash_levels,0,info,fieldout(pp_ptr),fieldout_col)
      END IF

      IF (mype  ==  0) THEN
        info = gc_none
        CALL gc_rrecv(100,num_stash_levels,start_pe,info,             &
                  fieldout(pp_ptr),fieldout_col)
      END IF
 
      pp_ptr = pp_ptr + num_stash_levels  ! increment the pp_ptr
    END IF  ! end of  MOD(i,num_stash_leveles) single point

    IF (mype  ==  0) THEN
      IF (start_ts) THEN ! put the descriptive info for this record
        ! DEPENDS ON: extra_make_vector
        CALL extra_make_vector(fake_record_extradata,control_size,           &
          i,no_records,fieldout(extra_start),extraw,bzx,bzy,bdx,bdy)
      END IF
    END IF

  END IF    ! end for column extract
END DO      ! end the loop over time series domain points 
!
horiz_size=pp_ptr-1
num_vert_levels=1
!  --------------------------------------------------------------------
!  7. If this is the first time in a time-series then
!      put the codes describing the extra data into the extra data fld
!      In addition set pphoriz out to the total length
!        as well as setting the input vetor to missing
!        where no values are set
! ----------------------------------------------------------------------

n_cols_out=no_records
n_rows_out=control(st_period_code)/control(st_freq_code)
horiz_size=n_cols_out
IF (start_ts) THEN  ! on start timestep we have entire vector
  horiz_size=n_cols_out*n_rows_out+extraw

  IF (mype  ==  0) THEN

    ! DEPENDS ON: extra_ts_info
    CALL extra_ts_info(fieldout(extra_start),extraw,no_records)
    DO i=no_records+1,extra_start-1
      fieldout(i)=amdi
    END DO

  ELSE
    DO i=no_records+1,lenout
      fieldout(i)=0.0
    END DO
  END IF

END IF
!
9999   CONTINUE ! jump here for error exit
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE multi_spatial
