! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! Takes a general decomposed field on many processors and gathers it
! to a single processor. Also supports the use case of a deferred gather
! to an IO server, when the final target of the gather is output to a dump
! file.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

SUBROUTINE general_gather_field (                               &
  local_field, global_field,                                    &
  local_size, global_size,                                      &
  levels,                                                       &
  addr_info,                                                    &
  gather_pe,                                                    &
  asyncOutputHandle,                                            &
  icode, cmessage)

USE IOS_Constants
USE IOS_stash_common
USE IOS_model_geometry, ONLY:                                &
    ios_dump_init
USE IOS_dump
USE mask_compression, ONLY:                                  &
    compress_to_mask,                                         &
    expand_from_mask
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod, ONLY: gc_none
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE atmos_max_sizes,   ONLY: Max2DFieldSize
USE atm_land_sea_mask, ONLY: atmos_landmask, atmos_landmask_local
USE umPrintMgr     , ONLY: &
    umMessage,              &
    umPrint,                &
    printstatus,            &
    prstatus_diag

USE cppxref_mod, ONLY:      &
    ppx_atm_compressed,     &
    ppx_atm_rim,            &
    ppx_atm_tall,           &
    ppx_atm_tland,          &
    ppx_atm_tsea,           &
    ppx_atm_uall,           &
    ppx_atm_cuall,          &
    ppx_atm_cvall,          &
    ppx_atm_tzonal,         &
    ppx_atm_uzonal,         &
    ppx_atm_tmerid,         &
    ppx_atm_ozone,          &
    ppx_atm_river,          &
    ppx_atm_scalar

USE d1_array_mod, ONLY: d1_list_len, d1_object_type, d1_length,   &
                        d1_grid_type, d1_no_levels, d1_north_code,&
                        d1_south_code, d1_east_code, d1_west_code,&
                        d1_gridpoint_code, d1_proc_no_code,       &
                        d1_halo_type, diagnostic
USE stparam_mod, ONLY: st_time_series_code, st_time_series_mean
USE stash_gather_field_mod, ONLY: stash_gather_field
USE nlsizes_namelist_mod, ONLY: &
    local_land_field

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments:

INTEGER, INTENT(IN)  ::  global_size         ! Size of GLOBAL FIELD
INTEGER, INTENT(IN)  ::  levels              ! How many levels of data to do
INTEGER, INTENT(IN)  ::  gather_pe(levels)   ! Which PE to gather each level to
INTEGER, INTENT(IN)  ::  asyncOutputHandle   ! Indicates this gather should
                                             ! go to an IOS directly via
                                             ! an async op
                                             ! -1  : no IOS op
                                             ! >=0 : handle of the async op
INTEGER, INTENT(OUT) ::  local_size          ! Size of LOCAL_FIELD
INTEGER, INTENT(OUT) ::  icode               ! Return code, 0=OK

INTEGER, INTENT(IN)  ::  addr_info(d1_list_len) ! Addressing info about field

REAL, INTENT(IN)     ::  local_field(*)         ! My local part of field
REAL, INTENT(OUT)    ::  global_field(global_size,*)
                                                ! Array to gather field to

CHARACTER(LEN=errormessagelength)    ::  cmessage               ! Error message

! Local variables
INTEGER   ::                                                      &
  grid_type                                                       &
            ! grid type of field being gathered
, grid_code                                                       &
            ! ppx grid code of field being gathered
, halo_type                                                       &
            ! halo type of field being gathered
, k                                                               &
       ! loop over levels
, my_k                                                            &
       ! value of k for GLOBAL_FIELD on this PE
, my_k_temp                                                       &
       ! Temp. my_k for safe references
, info                                                            &
       ! return code from GCOM routines
, dummy                                                           &
       ! dummy variables - ignored return values
, north,south,east,west                                           &
            ! domain limits for STASH output
, mean_type                                                       &
            ! spatial meaning type on diagnostic
, iproc                                                           &
            ! Loop variable over processors
, i         ! Loop variable for debugging

INTEGER       :: get_fld_type  ! function for finding field type
INTEGER       :: ios_subdomain
LOGICAL, SAVE :: haveInitialisedIOSLandMask=.FALSE.
REAL   ::                                                         &
  buf_expand(max2dfieldsize)                                      &
, buf_expand_local(max2dfieldsize)

LOGICAL, PARAMETER :: data_extracted = .TRUE.
            ! if the data in LOCAL_FIELD has been extracted

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GENERAL_GATHER_FIELD'

!===================================================================

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

grid_code=addr_info(d1_grid_type)
! DEPENDS ON: get_fld_type
grid_type=get_fld_type(grid_code)
halo_type=addr_info(d1_halo_type)

!-------------------------------------------------------------------

! Timeseries data
IF ((addr_info(d1_object_type)  ==  diagnostic) .AND.             &
    ((addr_info(d1_proc_no_code)  ==  st_time_series_code) .OR.    &
    (addr_info(d1_proc_no_code)  ==  st_time_series_mean))) THEN

  ! Multi-level fields not supported here
  IF (levels  >   1) THEN
    WRITE(umMessage,*) 'GENERAL_GATHER_FIELD : Cannot have more than ',   &
        '1 level for gathering timeseries data.'
    CALL umPrint(umMessage,src='general_gather_field')
    icode=1
    cmessage='GENERAL_GATHER_FIELD : Multi-level timeseries field'
    GO TO 9999
  END IF

  IF (asyncOutputHandle /= -1) THEN
    WRITE(umMessage,*) 'GENERAL_GATHER_FIELD : Cannot have time series ', &
        'in async dumps using IOS'
    CALL umPrint(umMessage,src='general_gather_field')
    icode=1
    cmessage= &
        'GENERAL_GATHER_FIELD : timeseries field with IOS not allowed'
    GO TO 9999
  END IF

  ! Copy the data to GLOBAL_FIELD from PE 0
  IF (mype  ==  0) THEN
    info=gc_none
    CALL gc_rsend(99,global_size,gather_pe,info,global_field,     &
        local_field)
  END IF

  IF (mype  ==  gather_pe(1)) THEN
    info=gc_none
    CALL gc_rrecv(99,global_size,0,info,global_field,local_field)
  END IF

  local_size=global_size

  !-------------------------------------------------------------------

  ! Surface (land points only) fields

ELSE IF (grid_code  ==  ppx_atm_compressed) THEN

  IF (.NOT. haveInitialisedIOSLandMask .AND. &
      (asyncOutputHandle /= -1)) THEN

    IF (mype == 0 .AND. printstatus >= prstatus_diag) THEN
      WRITE(umMessage,*)'Initialise landmask on IO Server...'
      CALL umPrint(umMessage,src='general_gather_field')
    END IF
    CALL ios_dump_init(atmos_landmask,&
        glsize(1,grid_type)*glsize(2,grid_type)&
        )
    haveInitialisedIOSLandMask=.TRUE.
  END IF

  my_k=0
  DO k=1,levels
    ! Unpack the local field out to full (local) field size and
    ! put this into the array buf_expand_local

    CALL expand_from_mask(buf_expand_local,          &
        local_field(1+(k-1)*local_land_field),    &
        atmos_landmask_local,                     &
        lasize(1,grid_type,halo_type)*            &
        lasize(2,grid_type,halo_type),            &
        dummy)

    ! Now gather in all the processors local fields into the global
    ! field (array buf_expand) or send to IO server

    IF (asyncOutputHandle /= -1) THEN
      CALL ios_dump_pack_data ( &
       asyncOutputHandle, &                 ! Communicatiosn handle
       buf_expand_local, &                  ! Data In
       arg_grid_type=grid_code, &           ! field type
       arg_halo_type=halo_type, &           ! halo type
       arg_preprocess_flag=ios_stash_preprocess, &! preprocess active
       arg_subdomain_flag=IOS_Full_Field, & ! subdomain control
       arg_landmask_compress=.TRUE. &
       )
    ELSE
      ! DEPENDS ON: gather_field
      CALL gather_field(buf_expand_local,buf_expand,    &
          lasize(1,grid_type,halo_type),                &
          lasize(2,grid_type,halo_type),                &
          glsize(1,grid_type),glsize(2,grid_type),      &
          grid_type,halo_type,                          &
          gather_pe(k),gc_all_proc_group )

      ! And now pack the global field (buf_expand) back to land points
      ! and put into the array GLOBAL_FIELD.

      IF (mype  ==  gather_pe(k)) THEN
        my_k=my_k+1

        CALL compress_to_mask(buf_expand,global_field(1,my_k), &
            atmos_landmask,                                   &
            glsize(1,grid_type)*glsize(2,grid_type),          &
            dummy)
      END IF
    END IF
  END DO ! k : loop over levels

  local_size = local_land_field

  !-------------------------------------------------------------------

  ! Atmosphere/ocean Lateral boundary fields

ELSE IF  (grid_code  ==  ppx_atm_rim) THEN

  WRITE(umMessage,*)                                                    &
    'GENERAL_GATHER_FIELD : Error attempting to gather ',       &
    'boundary data. Data not stored in dump from vn5.3'
  CALL umPrint(umMessage,src='general_gather_field')
  icode=3
  cmessage='GENERAL_GATHER_FIELD : Error gathering LBCs'
  GO TO 9999

  !-------------------------------------------------------------------

  ! "Normal" fields

ELSE IF                                                           &
  !     atmosphere grids
    ((grid_code  ==  ppx_atm_tall)  .OR.                            &
     (grid_code  ==  ppx_atm_tland) .OR.                            &
     (grid_code  ==  ppx_atm_tsea)  .OR.                            &
     (grid_code  ==  ppx_atm_uall) .OR.                             &
     (grid_code  ==  ppx_atm_cuall) .OR.                            &
     (grid_code  ==  ppx_atm_cvall) .OR.                            &
     (grid_code  ==  ppx_atm_tzonal) .OR.                           &
     (grid_code  ==  ppx_atm_uzonal) .OR.                           &
     (grid_code  ==  ppx_atm_tmerid) .OR.                           &
     (grid_code  ==  ppx_atm_compressed) .OR.                       &
     (grid_code  ==  ppx_atm_ozone) .OR.                            &
     (grid_code  ==  ppx_atm_river) .OR.                            &
     (grid_code  ==  ppx_atm_scalar)  )                             &
     THEN

  local_size=addr_info(d1_length)/addr_info(d1_no_levels)
  IOS_Subdomain=IOS_Full_Field

  IF (addr_info(d1_object_type)  ==  diagnostic) THEN
    north=addr_info(d1_north_code)
    south=addr_info(d1_south_code)
    east=addr_info(d1_east_code)
    west=addr_info(d1_west_code)

    mean_type=addr_info(d1_gridpoint_code)/10
    IF (mean_type  ==  2) THEN ! zonal mean
      IOS_Subdomain=IOS_Partial_Field
      east=west
    ELSE IF (mean_type  ==  3) THEN ! meridional mean
      IOS_Subdomain=IOS_Partial_Field
      north=south
    ELSE IF (mean_type  >=  4) THEN ! field/global mean
      IOS_Subdomain=IOS_Partial_Field
      east=west
      north=south
    ELSE IF (mean_type  ==  0 .AND. grid_code  ==  ppx_atm_scalar) THEN
      IOS_Subdomain=IOS_Partial_Field
      east=west
      north=south        
    END IF

  ELSE
    north=glsize(2,grid_type)
    west=1
    east=glsize(1,grid_type)
    south=1
  END IF

  my_k=0
  DO k=1,levels

    IF (mype  ==  gather_pe(k)) THEN
      my_k      = my_k+1
      my_k_temp = my_k
    ELSE
      my_k_temp = 1
    END IF

    IF (asyncOutputHandle /= -1) THEN

      CALL ios_dump_pack_data (            &
          asyncOutputHandle,               & ! Communications handle
          local_field(1+(k-1)*local_size : &
          (k  )*local_size),               & ! Data In
          arg_grid_type=grid_code,         & ! field type
          arg_halo_type=halo_type,         & ! halo type
          arg_preprocess_flag=             &
          ios_stash_preprocess,            & ! preprocess active
          arg_subdomain_flag=IOS_Subdomain,& ! subdomain control
          arg_S_boundary=south,            &
          arg_N_boundary=north,            &
          arg_W_boundary=west,             &
          arg_E_boundary=east,             &
          arg_landmask_compress=.FALSE.    &
          )
    ELSE
      !       STASH_GATHER_FIELD can distribute whole fields, or subarea
      !       fields
      CALL stash_gather_field(                                          &
          local_field(1+(k-1)*local_size),global_field(1,my_k_temp),    &
          local_size,global_size,1,                                     &
          north,east,south,west,                                        &
          grid_code,halo_type,gather_pe(k),data_extracted,              &
          icode=icode, cmessage=cmessage)

      IF (icode  /=  0) THEN
        WRITE(umMessage,*)                                                      &
            'GENERAL_GATHER_FIELD : Error detected in call to ',        &
            'STASH_GATHER_FIELD'
        CALL umPrint(umMessage,src='general_gather_field')
        WRITE(umMessage,*) 'Return code : ',icode
        CALL umPrint(umMessage,src='general_gather_field')
        WRITE(umMessage,*) 'Error message : ',cmessage
        CALL umPrint(umMessage,src='general_gather_field')

        icode=4
        cmessage='GENERAL_GATHER_FIELD : Error gathering field'
        GO TO 9999
      END IF
    END IF
  END DO ! k : loop over levels

  !-------------------------------------------------------------------
  ! Any other type of field
ELSE

  icode=10
  cmessage='GENERAL_GATHER_FIELD : Field type not recognized'

END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE general_gather_field

