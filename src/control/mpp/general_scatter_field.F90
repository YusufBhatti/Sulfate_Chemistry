! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Scatters any type of field from one processor to many processors

! Subroutine interface:

SUBROUTINE general_scatter_field (                                             &
  local_field, global_field, local_size, global_size, levels, grid_code,       &
  halo_type, object_type, proc_no_code, length, no_levels, north_code,         &
  south_code, east_code, west_code, gridpoint_code, scatter_pe, icode, cmessage)

USE mask_compression, ONLY: &
    compress_to_mask,           &
    expand_from_mask
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE d1_array_mod, ONLY: diagnostic
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod, ONLY: gc_none
USE stash_scatter_field_mod, ONLY: stash_scatter_field
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE atmos_max_sizes,   ONLY: Max2DFieldSize
USE atm_land_sea_mask, ONLY: atmos_landmask, atmos_landmask_local
USE rimtypes
USE lbc_mod
USE cppxref_mod, ONLY:                                         &
    ppx_atm_lbc_u, ppx_atm_lbc_v,                              &
    ppx_atm_lbc_theta, ppx_atm_lbc_orog,                       &
    ppx_atm_tall, ppx_atm_tland, ppx_atm_tsea,                 &
    ppx_atm_uall, ppx_atm_cuall, ppx_atm_cvall,                &
    ppx_atm_tzonal, ppx_atm_uzonal, ppx_atm_tmerid,            &
    ppx_atm_ozone, ppx_atm_river, ppx_atm_compressed,          &
    ppx_atm_scalar
USE stparam_mod, ONLY: st_time_series_code, st_time_series_mean
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    local_land_field

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Description:
! Takes a general field on a single processor and decomposes (scatters)
! it over many processors

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine arguments:

INTEGER, INTENT(IN) ::  global_size        ! IN:  size of GLOBAL FIELD
INTEGER, INTENT(IN) ::  levels             ! IN:  How many levels of data to do
INTEGER, INTENT(IN) ::  scatter_pe(levels) ! IN: Which PE to scatter each level from
INTEGER, INTENT(OUT) :: local_size         ! OUT: size of LOCAL_FIELD
INTEGER, INTENT(OUT) :: icode              ! OUT: return code, 0=OK

INTEGER, INTENT(IN) ::  grid_code      ! IN: ppx grid code of field
INTEGER, INTENT(IN) ::  halo_type      ! IN: halo type of field
INTEGER, INTENT(IN) ::  object_type    ! IN: field type
INTEGER, INTENT(IN) ::  proc_no_code   ! IN: Processing Code
INTEGER, INTENT(IN) ::  length         ! IN: Record length of field
INTEGER, INTENT(IN) ::  no_levels      ! IN: number of levels in field
INTEGER, INTENT(IN) ::  north_code     ! IN: northern row of field
INTEGER, INTENT(IN) ::  south_code     ! IN: southern row of field
INTEGER, INTENT(IN) ::  east_code      ! IN: eastern row of field
INTEGER, INTENT(IN) ::  west_code      ! IN: western row of field
INTEGER, INTENT(IN) ::  gridpoint_code ! IN: gridpoint info of field

REAL, INTENT(IN) :: global_field(global_size,*) ! IN:  field to scatter
REAL, INTENT(OUT) :: local_field(*)             ! OUT: my local part of fld

CHARACTER(LEN=errormessagelength) :: cmessage   ! OUT : Error message

! Local variables

INTEGER  ::                                                       &
  grid_type                                                       &
            ! grid type of field being scattered


, info                                                            &
       ! return code from GCOM routines
, dummy                                                           &
        ! dummy variables - ignored return values
, north,south,east,west                                           &
                         ! domain limits for STASH output
, mean_type                                                       &
            ! spatial meaning type on diagnostic
, k                                                               &
     ! loop over levels
, my_k                                                            &
       ! value of k for GLOBAL_FIELD on this PE
, my_k_temp                                                       &
            ! Temp. my_k for safe references.
, rim_type                                                        &
            ! RIMWIDTH type for LBC field
, loc_len_rim                                                     &
              ! length of local rimdata for single level
, iproc                                                           &
            ! Loop variable over processors
, i         ! Loop variable for debugging

INTEGER  ::                                                       &
  get_fld_type  ! function for finding field type

REAL     ::                                                       &
  buf_expand(max2dfieldsize)                                      &
, buf_expand_local(max2dfieldsize)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GENERAL_SCATTER_FIELD'

!===================================================================

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! DEPENDS ON: get_fld_type
grid_type=get_fld_type(grid_code)

!-------------------------------------------------------------------

! Timeseries data

IF ((object_type  ==  diagnostic) .AND.             &
   ((proc_no_code  ==  st_time_series_code) .OR.     &
   (proc_no_code  ==  st_time_series_mean))) THEN
  ! Copy the data from GLOBAL_FIELD to PE 0

  ! Multi-level fields not supported here
  IF (levels  >   1) THEN
    WRITE(umMessage,*) 'GENERAL_SCATTER_FIELD : Cannot have more than ',  &
               '1 level for scattering timeseries data.'
    CALL umPrint(umMessage,src='general_scatter_field')
    icode=1
    cmessage='GENERAL_SCATTER_FIELD : Multi-level timeseries field'
    GO TO 9999
  END IF

  IF (mype  ==  scatter_pe(1)) THEN

    info=gc_none
    CALL gc_rsend(99,global_size,0,info,local_field,              &
                  global_field)

  END IF

  IF (mype  ==  0) THEN

    info=gc_none
    CALL gc_rrecv(99,global_size,scatter_pe,info,                 &
                  local_field,global_field)

  END IF

  local_size=global_size

  !-------------------------------------------------------------------

  ! Surface (land points only) fields

ELSE IF                                                           &
  (grid_code  ==  ppx_atm_compressed) THEN

  my_k=0
  DO k=1,levels

    IF (mype  ==  scatter_pe(k)) THEN

      my_k=my_k+1

      CALL expand_from_mask(buf_expand,global_field(1,my_k),      &
          atmos_landmask,                                      &
          glsize(1,grid_type)*glsize(2,grid_type),             &
          dummy)

      ! SCATTER_PE now contains the expanded version of the full field

    END IF

    ! Now scatter this to all the other processors, putting the local
    ! part of the field into the array buf_expand_local

    ! DEPENDS ON: scatter_field
    CALL scatter_field(buf_expand_local , buf_expand,             &
                       lasize(1,grid_type,halo_type),             &
                       lasize(2,grid_type,halo_type),             &
                       glsize(1,grid_type),glsize(2,grid_type),   &
                       grid_type,halo_type,                       &
                       scatter_pe(k),gc_all_proc_group)

    ! No need to call swap_bounds as field will be recompressed

    ! Pack the local field down to local land points and put
    ! the packed field into LOCAL_FIELD

#if defined (MERGE)
    local_land_field=dummy
#endif

    CALL compress_to_mask(buf_expand_local,          &
        local_field(1+(k-1)*local_land_field),      &
        atmos_landmask_local,                       &
        lasize(1,grid_type,halo_type)*              &
        lasize(2,grid_type,halo_type),              &
        dummy)
  END DO ! k : loop over levels

  local_size = local_land_field

  !-------------------------------------------------------------------

  ! Atmosphere Lateral boundary fields

ELSE IF                                                           &
  ((grid_code  ==  ppx_atm_lbc_theta) .OR.                        &
   (grid_code  ==  ppx_atm_lbc_u) .OR.                            &
   (grid_code  ==  ppx_atm_lbc_v) .OR.                            &
   (grid_code  ==  ppx_atm_lbc_orog)) THEN

  IF (grid_code  ==  ppx_atm_lbc_orog) THEN
    rim_type=rima_type_orog
  ELSE
    rim_type=rima_type_norm
  END IF

  ! DEPENDS ON: scatter_atmos_lbcs
  CALL scatter_atmos_lbcs(global_field,local_field,               &
                          global_lenrima(grid_type,halo_type,     &
                                         rim_type),               &
                          lenrima(grid_type,halo_type,rim_type),  &
                          levels,levels,                          &
                          grid_type,halo_type,rim_type,           &
                          scatter_pe,                             &
                          icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*)                                                    &
      'GENERAL_SCATTER_FIELD : Error detected in call to ',       &
      'SCATTER_ATMOS_LBCS '
    CALL umPrint(umMessage,src='general_scatter_field')
    WRITE(umMessage,*) 'Return code : ',icode
    CALL umPrint(umMessage,src='general_scatter_field')
    WRITE(umMessage,*) 'Error message : ',cmessage
    CALL umPrint(umMessage,src='general_scatter_field')

    icode=2
    cmessage='GENERAL_SCATTER_FIELD : Error scattering LBCs'
    GO TO 9999
  END IF


  local_size=lenrimdata_a

  !-------------------------------------------------------------------

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
     (grid_code  ==  ppx_atm_ozone) .OR.                            &
     (grid_code  ==  ppx_atm_river) .OR.                            &
     (grid_code  ==  ppx_atm_scalar)  )                             &
     THEN

  local_size=length/no_levels


  IF (object_type  ==  diagnostic) THEN
    north=north_code
    south=south_code
    east=east_code
    west=west_code

    mean_type=gridpoint_code/10
    IF (mean_type  ==  2) THEN ! zonal mean
      east=west
    ELSE IF (mean_type  ==  3) THEN ! meridional mean
      south=north
    ELSE IF (mean_type  >=  4) THEN ! field/global mean
      east=west
      south=north
    ELSE IF (mean_type  ==  0 .AND. grid_code  ==  ppx_atm_scalar) THEN
      east=west
      south=north 
    END IF

  ELSE
    north=glsize(2,grid_type)
    west=1
    east=glsize(1,grid_type)
    south=1
  END IF

  my_k=0
  DO k=1,levels

    IF (mype  ==  scatter_pe(k)) THEN
      my_k=my_k+1
      my_k_temp = my_k
    ELSE
      my_k_temp=1
    END IF

    !       STASH_SCATTER_FIELD can distribute whole fields, or subarea
    !       fields
    CALL stash_scatter_field(                                     &
      local_field(1+(k-1)*local_size),global_field(1,my_k_temp),  &
      local_size,global_size,1,                                   &
      north,east,south,west,                                      &
      grid_code,halo_type,scatter_pe(k),icode,cmessage)

    IF (icode  /=  0) THEN
      WRITE(umMessage,*)                                                  &
        'GENERAL_SCATTER_FIELD : Error detected in call to ',     &
        'STASH_SCATTER_FIELD '
      CALL umPrint(umMessage,src='general_scatter_field')
      WRITE(umMessage,*) 'Return code : ',icode
      CALL umPrint(umMessage,src='general_scatter_field')
      WRITE(umMessage,*) 'Error message : ',cmessage
      CALL umPrint(umMessage,src='general_scatter_field')

      icode=4
      cmessage='GENERAL_SCATTER_FIELD : Error scattering field'
      GO TO 9999
    END IF
  END DO ! k : loop over levels

  !-------------------------------------------------------------------
  ! Any other type of field
ELSE

  icode=10
  cmessage='GENERAL_SCATTER_FIELD : Field type not recognized'

END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE general_scatter_field
