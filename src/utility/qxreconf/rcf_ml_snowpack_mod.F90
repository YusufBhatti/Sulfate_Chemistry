! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets variables of multi-layer snow scheme on snow levels.

MODULE rcf_ml_snowpack_mod
IMPLICIT NONE

! Description:
!     This subroutine sets variables of the multi-layer snow scheme 
!     that are defined on tiles. The snowpack is divided into layers
!     in a manner similar to that in the snow scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ML_SNOWPACK_MOD'

CONTAINS

SUBROUTINE rcf_ml_snowpack(  fields_in, field_count_in, hdr_in,    &
                             fields_out, field_count_out, hdr_out, &
                             data_source )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_grid_type_mod, ONLY: &
    grid_type, &
    Input_Grid, &
    Output_Grid

USE rcf_set_dominant_tile_mod, ONLY:       &
    rcf_set_dominant_tile, i_dominant

USE Rcf_horizontal_mod, ONLY: &
    Rcf_horizontal

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE um_stashcode_mod, ONLY:    &
    stashcode_prog_sec          , &
    stashcode_soil_temp         , &
    stashcode_snow_tile         , &
    stashcode_snow_grnd         , &
    stashcode_snowdep_grd_tile  , &
    stashcode_snowpack_bk_dens  , &
    stashcode_nsnow_layrs_tiles , &
    stashcode_snow_laythk_tiles , &
    stashcode_snow_ice_tile     , &
    stashcode_snow_liq_tile     , &
    stashcode_snow_t_tile       , &
    stashcode_snow_laydns_tiles , &
    stashcode_snow_grnsiz_tiles

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_init_snow_bk_dens_mod, ONLY: &
    rcf_init_snow_bk_dens

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE decomp_params, ONLY: &
    decomp_rcf_input,    &
    decomp_rcf_output

USE rcf_data_source_Mod, ONLY: &
    data_source_type

USE rcf_set_interp_flags_mod, ONLY: &
    interp_h_only

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE rcf_nlist_recon_science_mod, ONLY: &
    l_regularize_landice,              &
    snow_icefree_max,                  &
    snow_landice_min

USE rcf_nlist_recon_technical_mod, ONLY: &
    l_rcf_init_flexi

USE rcf_lsm_mod, ONLY: &
    local_land_out, local_land_in

USE rcf_calc_tile_map_mod, ONLY:           &
    rcf_calc_tile_map

USE jules_surface_types_mod, ONLY: tile_map_ids, tile_map_pslevs, ice

USE water_constants_mod,     ONLY: tm, lf, hcapw, hcapi

USE nlsizes_namelist_mod, ONLY: &
    ntiles

USE jules_vegetation_mod, ONLY: can_model

USE jules_snow_mod, ONLY: nsmax, rho_snow_const, rho_snow_fresh,    &
                          cansnowtile, dzsnow, r0, rmax

USE science_fixes_mod, ONLY: l_fix_rcf_mlsnow_icefreemax

USE umPrintMgr, ONLY: umPrint, umMessage, printstatus, prstatus_diag

USE um_parcore, ONLY: mype

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_in(:), fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_in, hdr_out
INTEGER, INTENT(IN)               :: field_count_in, field_count_out
TYPE( data_source_type ), POINTER :: data_source( : )


! Internal variables
TYPE( field_type ), POINTER  :: soil_temp

TYPE( field_type ), POINTER  :: snowdepth
TYPE( field_type ), POINTER  :: snowbulkdens
TYPE( field_type ), POINTER  :: snow_tile
TYPE( field_type ), POINTER  :: snow_grnd

TYPE( field_type ), POINTER  :: nsnow
TYPE( field_type ), POINTER  :: ds
TYPE( field_type ), POINTER  :: sice
TYPE( field_type ), POINTER  :: sliq
TYPE( field_type ), POINTER  :: rhosnow
TYPE( field_type ), POINTER  :: tsnow
TYPE( field_type ), POINTER  :: rgrainl

TYPE (grid_type)                 :: grid_middle
TYPE (field_type)                :: field_middle

TYPE( field_type ), POINTER      :: input_field

INTEGER                          :: tile_map_pslevs_tmp(ntiles)

INTEGER                          :: ntiles_in

INTEGER                          :: pos_in, pos_out, pos

LOGICAL                          :: l_ml_snow_in
LOGICAL                          :: l_from_input

REAL, ALLOCATABLE                :: nsnow_in(:,:)

REAL, ALLOCATABLE                :: ds_in(:,:)
REAL, ALLOCATABLE                :: sice_in(:,:)
REAL, ALLOCATABLE                :: sliq_in(:,:)
REAL, ALLOCATABLE                :: rhosnow_in(:,:)
REAL, ALLOCATABLE                :: tsnow_in(:,:)
REAL, ALLOCATABLE                :: rgrainl_in(:,:)

INTEGER                          :: i, j, n
INTEGER                          :: jin
INTEGER                          :: kin, kout
INTEGER                          :: pseudo_in, pseudo_out

REAL                             :: remains
REAL                             :: depth_in
REAL                             :: depth_out
REAL                             :: depth_done
REAL                             :: dz
REAL                             :: water_mass
REAL                             :: ice_mass
REAL                             :: dsnowmass

REAL                :: total_snow_ground( local_land_out, ntiles )

REAL, PARAMETER                  :: rgrain_icesheet = 150.0
!                                   ! Typical grain size on ice sheets
!                                   ! Used only to initialize the
!                                   ! grain size at points converted
!                                   ! from ice-free to land ice

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_ML_SNOWPACK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
  WRITE(umMessage,'(a,i3)') 'Calculating multilayer snow fields'
  CALL umPrint(umMessage,src='rcf_ml_snowpack_mod')
END IF

!-----------------------------------------------------------------------
! If multilayer fields are in the input dump, make use of them by
! interpolating to the output grid on the original tiles.
!-----------------------------------------------------------------------

! Check for the presence of these fields by looking for snow depth
! in snow layers.
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles,      &
           fields_in, field_count_in, pos_in, zero_ok_arg = .TRUE.)
l_ml_snow_in = (pos_in > 0)

! Set the number of tiles in the input from the snow amount on tiles.
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_tile,                 &
           fields_in, field_count_in, pos, .TRUE.)
IF (pos > 0) THEN
  ntiles_in = fields_in(pos) % levels
ELSE
! The tiled field is not present, so the input dump can only
! have a single value for the whole grid-box.
  ntiles_in = 1
END IF

IF (l_ml_snow_in) THEN

  !---------------------------------------------------------------------
  ! Set up the template for interpolation. This requires the horizontal
  ! size of the output grid, but the number of snow levels of the
  ! input grid to ensure that information is retained.
  !--------------------------------------------------------------------

  IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
    WRITE(umMessage,'(a,i3)') 'Multilayer snow fields present in input dump'
    CALL umPrint(umMessage,src='rcf_ml_snowpack_mod')
  END IF

  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles,      &
             fields_out, field_count_out, pos_out )

  grid_middle = Output_Grid

  field_middle % levels          =  fields_in(pos_in)   % levels
  field_middle % bottom_level    =  fields_in(pos_in)   % bottom_level
  field_middle % top_level       =  fields_in(pos_in)   % top_level
  field_middle % interp          =  fields_in(pos_in)   % interp
  field_middle % rows            =  fields_out(pos_out) % rows
  field_middle % row_len         =  fields_out(pos_out) % row_len
  field_middle % level_size      =  fields_out(pos_out) % level_size
  field_middle % glob_rows       =  fields_out(pos_out) % glob_rows
  field_middle % glob_row_len    =  fields_out(pos_out) % glob_row_len
  field_middle % glob_level_size =  fields_out(pos_out) % glob_level_size
  field_middle % stashmaster     => fields_in(pos_in)   % stashmaster

! Allocate field_middle Data space
  CALL Rcf_Alloc_Field( field_middle )

!------------------------------------------------------------------
! Use field_middle to interpolate all the multilayer fields in the
! horizontal.
!------------------------------------------------------------------

  ALLOCATE( nsnow_in(local_land_out,   fields_in(pos_in) % levels) )
  ALLOCATE( ds_in(local_land_out,      fields_in(pos_in) % levels) )
  ALLOCATE( sice_in(local_land_out,    fields_in(pos_in) % levels) )
  ALLOCATE( sliq_in(local_land_out,    fields_in(pos_in) % levels) )
  ALLOCATE( tsnow_in(local_land_out,   fields_in(pos_in) % levels) )
  ALLOCATE( rhosnow_in(local_land_out, fields_in(pos_in) % levels) )
  ALLOCATE( rgrainl_in(local_land_out, fields_in(pos_in) % levels) )

! Number of snow layers
  CALL rcf_locate( stashcode_prog_sec, stashcode_nsnow_layrs_tiles,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  nsnow_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

! Thicknesses of layers
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  ds_in = field_middle % DATA

! May have arrived here directly because of incorrect indexing of snow layers
! so need to check that the tile_map_ids is correct for the input configuration
! before using it. Use the snow layer thickness to check this as it triggered
! the call to this routine. rcf_calc_tile_map also needs to be called to
! generate the corresponding tile_map_pslevs containing pseudo level numbers.
  IF ( ntiles_in > 1 .AND. ntiles > 1 ) THEN
    CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles,      &
                     fields_out, field_count_out, pos)
    CALL rcf_calc_tile_map( input_field, hdr_in, ntiles_in,                &
                            fields_out(pos), stashcode_snow_laythk_tiles )
  END IF
  CALL rcf_dealloc_field( input_field )

! Frozen content of snow
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_ice_tile,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  sice_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

! Liquid content of snow
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_liq_tile,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  sliq_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

! Temperature of snow
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_t_tile,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  tsnow_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

! Density of snow
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laydns_tiles,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  rhosnow_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

! Snow grain size
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_grnsiz_tiles,      &
             fields_in, field_count_in, pos)
  input_field => fields_in(pos)
  CALL rcf_alloc_field( input_field )
  CALL rcf_read_field( input_field, hdr_in, decomp_rcf_input )
  input_field % interp = interp_h_only
  CALL Rcf_horizontal( input_field, field_middle, Input_Grid, grid_middle )
  rgrainl_in = field_middle % DATA
  CALL rcf_dealloc_field( input_field )

  IF ( ntiles_in > 1 .AND. ntiles == 1 ) THEN
    IF ( .NOT. ALLOCATED(i_dominant) ) THEN
      ALLOCATE( i_dominant(fields_in(pos) % level_size) )

      !     Find the dominant tile at each point.
      CALL Rcf_Set_Dominant_Tile( fields_in, field_count_in, hdr_in, &
                                  fields_out, field_count_out, hdr_out )
    END IF
  END IF

END IF


! Get the snow stores, recalling that they will have been checked
! for consistency.
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_tile,                   &
           fields_out, field_count_out, pos)
snow_tile => fields_out(pos)
CALL rcf_alloc_field( snow_tile )
CALL rcf_read_field( snow_tile, hdr_out, decomp_rcf_output )
total_snow_ground(:,:) = snow_tile % DATA(:,:)
!
IF (can_model == 4) THEN
! Snow under canopy
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_grnd,                 &
             fields_out, field_count_out, pos)
  snow_grnd => fields_out(pos)
  CALL rcf_alloc_field( snow_grnd )
  CALL rcf_read_field( snow_grnd, hdr_out, decomp_rcf_output )
  DO n=1,ntiles
    IF ( cansnowtile(n) ) &
      total_snow_ground(:,n) = snow_grnd % DATA(:,n)
  END DO
END IF

! Soil temperature will be required to set the temperature 
! of the snow pack where there is no input data.
CALL rcf_locate( stashcode_prog_sec, stashcode_soil_temp,                   &
           fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
CALL rcf_alloc_field( soil_temp )
CALL rcf_read_field( soil_temp, hdr_out, decomp_rcf_output )

!-------------------------------------------------------------------
! Read the single-layer fields which will exist on 
! the output grid by now.
!-------------------------------------------------------------------

! Snow depth
CALL rcf_locate( stashcode_prog_sec, stashcode_snowdep_grd_tile,      &
            fields_out, field_count_out, pos)
snowdepth => fields_out(pos)
CALL rcf_alloc_field( snowdepth )
CALL rcf_read_field( snowdepth, hdr_out, decomp_rcf_output )

! Bulk density
CALL rcf_locate( stashcode_prog_sec, stashcode_snowpack_bk_dens,      &
            fields_out, field_count_out, pos)
snowbulkdens => fields_out(pos)
CALL rcf_alloc_field( snowbulkdens )
CALL rcf_read_field( snowbulkdens, hdr_out, decomp_rcf_output )

!---------------------------------------------------------------
! Set the tile map from the input to the output.
!---------------------------------------------------------------
! Compatability with m-to-n reconfiguration:
! tile_map_pslevs not set therefore called from either 1-to-n or n-to-1
! rcf_calc_tiles. Set tile_map_pslevs to map output tiles to single tile in
! input. ( Not needed if called from n-to-1, but does no harm )
! This is not compatible with l_rcf_init_flexi (fail at rcf_reset_data_source).
tile_map_pslevs_tmp(:)=tile_map_pslevs(1:ntiles)
IF ( MAXVAL(tile_map_pslevs) < 0 ) tile_map_pslevs_tmp(:) = 1

!---------------------------------------------------------------
! Allocate fields for output.
!---------------------------------------------------------------
CALL rcf_locate( stashcode_prog_sec, stashcode_nsnow_layrs_tiles,      &
            fields_out, field_count_out, pos)
nsnow => fields_out(pos)
CALL rcf_alloc_field( nsnow )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles,      &
            fields_out, field_count_out, pos)
ds => fields_out(pos)
CALL rcf_alloc_field( ds )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_ice_tile,          &
            fields_out, field_count_out, pos)
sice => fields_out(pos)
CALL rcf_alloc_field( sice )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_liq_tile,          &
            fields_out, field_count_out, pos)
sliq => fields_out(pos)
CALL rcf_alloc_field( sliq )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_t_tile,            &
            fields_out, field_count_out, pos)
tsnow => fields_out(pos)
CALL rcf_alloc_field( tsnow )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laydns_tiles,      &
            fields_out, field_count_out, pos)
rhosnow => fields_out(pos)
CALL rcf_alloc_field( rhosnow )
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_grnsiz_tiles,      &
            fields_out, field_count_out, pos)
rgrainl => fields_out(pos)
CALL rcf_alloc_field( rgrainl )

DO j = 1, ntiles

  DO i = 1, local_land_out

    IF ( (ntiles == 1) .AND. (ntiles_in > 1) ) THEN
      jin = i_dominant(i)
    ELSE
      jin = tile_map_pslevs_tmp(j)
    END IF

!   ---------------------------------------------------------------
!   Divide the snow pack into layers: this code follows the
!   routine layersnow in JULES.
!   ---------------------------------------------------------------
    IF ( snowdepth % DATA(i,j) >= dzsnow(1) ) THEN

      remains = snowdepth % DATA(i,j)
      DO n=1,nsmax
        ! pseudo-level index for NTILES*NSMAX dimension
        pseudo_out = j + (n-1) * ntiles
        ds % DATA(i,pseudo_out) = dzsnow(n)
        remains = remains - dzsnow(n)
        IF ( remains <= dzsnow(n) .OR. n == nsmax ) THEN
          ds % DATA(i,pseudo_out) = ds % DATA(i,pseudo_out) + remains
          EXIT
        END IF
      END DO

      ! Set number of layers.
      nsnow % DATA(i,j) = REAL(n)

!     Decide whether the fields at this point can be set from the input,
!     or whether they must be initialized. A two-stage test is required
!     in case ds_in has not been allocated.
      l_from_input = l_ml_snow_in
      IF ( l_from_input ) l_from_input = (ds_in(i, jin) > 0.0)

      IF (l_from_input ) THEN

!       ---------------------------------------------------------------
!       If there are data on multiple layers at this point,
!       remap the input fields on the output grid to the new layers:
!       this code does not quite need to follow relayersnow in JULES
!       since we are not adding any new material at the top of the
!       snow pack. Note that the input fields are now on the right
!       horizontal grid, but the original tiles.
!       ---------------------------------------------------------------
!
!       Initialize the depth of the snow pack assigned (depth_done)
!       and the depth of the top of the current layer in the new
!       decomposition.
        depth_out  = 0.0
        depth_done = 0.0
!
!       Set up the current layer in the input snow pack.
        kin        = 1
        pseudo_in  = jin + (kin-1) * ntiles_in
        depth_in   = ds_in(i, pseudo_in)
!
!       Loop through output layers, pulling in contributions from the 
!       input layers.
        DO kout = 1, NINT(nsnow % DATA(i,j))
!         Set depth to the bottom of the output layer.
          pseudo_out = j + (kout-1) * ntiles
          depth_out = depth_out + ds % DATA(i, pseudo_out)
!
!         Accumulate conserved quantities in current new layer.
          water_mass                    = 0.0
          ice_mass                      = 0.0
          tsnow   % DATA(i, pseudo_out) = 0.0
          rgrainl % DATA(i, pseudo_out) = 0.0
          DO
            dz = MIN(depth_in, depth_out) - depth_done
            depth_done = depth_done + dz
            water_mass = water_mass + (dz / ds_in(i, pseudo_in) ) * &
              (sice_in(i, pseudo_in) + sliq_in(i, pseudo_in))
            ice_mass = ice_mass + (dz / ds_in(i, pseudo_in) ) * &
              sice_in(i, pseudo_in)
!           Average the inverse of the grain size to preserve the 
!           snow specific area.
            rgrainl % DATA(i, pseudo_out) = rgrainl % DATA(i, pseudo_out) + &
              sice_in(i, pseudo_in) * ( dz / ds_in(i, pseudo_in) ) / &
              rgrainl_in(i, pseudo_in)
!           In the first instance calculate the enthalpy divided by hcapi
!           to deal with melting.
            tsnow % DATA(i, pseudo_out) = tsnow % DATA(i, pseudo_out) + &
              (sice_in(i, pseudo_in) * (dz / ds_in(i, pseudo_in) ) ) * &
              tsnow_in(i, pseudo_in) + &
              (sliq_in(i, pseudo_in) * (dz / ds_in(i, pseudo_in) ) ) * &
              (lf / hcapi + tm + (hcapw / hcapi) * &
              (tsnow_in(i, pseudo_in) - tm) )

!           Note. Use EPSILON with a margin of safety.
            IF ( depth_done >= depth_out - 64.0 * EPSILON(depth_out) ) EXIT
!           Move on to the next layer in the input if it's valid.
            IF (kin < NINT(nsnow_in(i, jin))) THEN
              kin=kin+1
              pseudo_in  = jin + (kin-1) * ntiles_in
              depth_in  = depth_in  + ds_in(i, pseudo_in)
            ELSE
              depth_in  = snowdepth % DATA(i,j)
            END IF
          END DO
!
          rhosnow % DATA(i, pseudo_out) = water_mass / ds % DATA(i, pseudo_out)
          rgrainl % DATA(i, pseudo_out) =   ice_mass / &
                                     rgrainl % DATA(i, pseudo_out)
          rgrainl % DATA(i, pseudo_out) = MAX(r0, MIN(rmax, &
            rgrainl % DATA(i, pseudo_out)) )
!         We now have the enthalpy of each new layer. If this exceeds that
!         which it would have if all the mass were frozen at the freezing
!         point, then some water must be in the liquid state and the
!         temperature of the layer must be the freezing point. If it is
!         lower then all the water will be frozen.
          IF (tsnow % DATA(i, pseudo_out) < water_mass * tm) THEN
!           The temperature will be below freezing, so assume all liquid
!           is frozen.
            sliq % DATA(i, pseudo_out)  = 0.0
            sice % DATA(i, pseudo_out)  = water_mass
            tsnow % DATA(i, pseudo_out) = tsnow % DATA(i, pseudo_out) / &
              water_mass
          ELSE
!           Set the temperature to freezing, melting enough ice to do this.
            sliq % DATA(i, pseudo_out) = (tsnow % DATA(i, pseudo_out) - &
              water_mass * tm) / (lf / hcapi)
            sice % DATA(i, pseudo_out) =  water_mass - &
              sliq % DATA(i, pseudo_out)
            tsnow % DATA(i, pseudo_out) = tm
          END IF

        END DO

      ELSE

        DO n = 1, NINT(nsnow % DATA(i,j))
          pseudo_out = j + (n-1) * ntiles
          IF ( snowdepth % DATA(i,j) > EPSILON(snowdepth % DATA) ) THEN
            sice % DATA(i, pseudo_out)    = total_snow_ground(i,j) * &
              ds % DATA(i, pseudo_out) / snowdepth % DATA(i,j)
            rhosnow % DATA(i, pseudo_out) = rho_snow_const
          ELSE
            sice % DATA(i, pseudo_out)    = 0.0
            rhosnow % DATA(i, pseudo_out) = rho_snow_fresh
          END IF
          sliq % DATA(i, pseudo_out)    = 0.0
          rgrainl % DATA(i, pseudo_out) = r0
          tsnow % DATA(i, pseudo_out)   = &
            MIN( soil_temp % DATA(i,1), tm )
        END DO

      END IF

    ELSE

      nsnow % DATA(i,j) = 0.0

    END IF    !  >dzSnow(1)
!
!   Set sensible values in empty layers
    DO n = NINT(nsnow % DATA(i,j)) + 1, nsmax
      pseudo_out = j + (n-1) * ntiles
      ds % DATA(i, pseudo_out)         = 0.0
      sice % DATA(i, pseudo_out)       = 0.0
      sliq % DATA(i, pseudo_out)       = 0.0
      tsnow % DATA(i, pseudo_out)      = &
        MIN( soil_temp % DATA(i,1), tm )
      rhosnow % DATA(i, pseudo_out)    = rho_snow_fresh
      rgrainl % DATA(i, pseudo_out)    = r0
    END DO

  END DO

END DO

IF (l_regularize_landice) THEN
!
  DO j=1, ntiles
!
    IF (j==ice) THEN
!     
!     If this is an ice point, ensure that the mass of snow is large
!     enough: if it has been reconfigured from an ice-free point this
!     may not be the case.
      DO i=1, local_land_out
        IF ( snow_tile % DATA(i,j) < snow_landice_min ) THEN
          IF ( nsnow % DATA(i,j) < REAL(nsmax) ) THEN
!           The ice is so thin that the snow pack should be reset.
            snow_tile % DATA(i,j) = snow_landice_min
            IF (can_model == 4) snow_grnd % DATA(i,j) = 0.0
            snowbulkdens % DATA(i,j) = rho_snow_const
            snowdepth % DATA(i,j)    = snow_landice_min / rho_snow_const
            nsnow % DATA(i,j)        = REAL(nsmax)
            DO n=1, nsmax
              pseudo_out = j + (n-1) * ntiles
              ds % DATA(i, pseudo_out) = dzsnow(n)
              IF (n==nsmax) THEN
                ds % DATA(i, pseudo_out) = &
                  snowdepth % DATA(i,j) - SUM(dzsnow(1:nsmax-1))
              END IF
              sice % DATA(i, pseudo_out) = &
                rho_snow_const * ds % DATA(i, pseudo_out)
              sliq % DATA(i, pseudo_out) = 0.0
!             Keep it frozen.
              tsnow % DATA(i, pseudo_out) = &
                MIN( soil_temp % DATA(i,1), tm-1.0 )
              rhosnow % DATA(i, pseudo_out) = rho_snow_const
!             A typical value for ice sheets.
              rgrainl % DATA(i, pseudo_out) = rgrain_icesheet
            END DO
          ELSE
!           The snow pack is essentially correct. Simply reset the
!           lowest layer.
            dsnowmass = snow_landice_min - snow_tile % DATA(i,j)
            pseudo_out = j + (nsmax-1) * ntiles
            sice % DATA(i, pseudo_out) = &
              sice % DATA(i, pseudo_out) + dsnowmass
            ds % DATA(i, pseudo_out) = &
              ds % DATA(i, pseudo_out) + &
              dsnowmass / rhosnow % DATA(i, pseudo_out)
            snow_tile % DATA(i,j) = snow_landice_min
            snowdepth % DATA(i,j) = snowdepth % DATA(i,j) + &
              dsnowmass / rhosnow % DATA(i, pseudo_out)
            snowbulkdens % DATA(i,j) = &
              snow_landice_min / snowdepth % DATA(i,j)
          END IF
        END IF
      END DO
!
    ELSE
!     
!     If this is an ice-free point, ensure that the mass of snow is not
!     too large: if it has been reconfigured from a land ice point this
!     may not be the case.
      DO i=1, local_land_out
        IF ( snow_tile % DATA(i,j) > snow_icefree_max ) THEN
          IF ( (can_model == 4) .AND. cansnowtile(n) ) THEN
            snow_tile % DATA(i,j) = 0.0
            IF (l_fix_rcf_mlsnow_icefreemax) THEN
              snow_grnd % DATA(i,j) = snow_icefree_max
            ELSE
              snow_grnd % DATA(i,j) = snow_landice_min
            END IF
          ELSE
            IF (l_fix_rcf_mlsnow_icefreemax) THEN
              snow_tile % DATA(i,j) = snow_icefree_max
            ELSE
              snow_tile % DATA(i,j) = snow_landice_min
            END IF
          END IF
          snowbulkdens % DATA(i,j) = rho_snow_const
          snowdepth % DATA(i,j)    = snow_icefree_max / rho_snow_const
          nsnow % DATA(i,j)        = REAL(nsmax)
          DO n=1, nsmax
            pseudo_out = j + (n-1) * ntiles
            ds % DATA(i, pseudo_out) = dzsnow(n)
            IF (n==nsmax) THEN
              ds % DATA(i, pseudo_out) = &
                snowdepth % DATA(i,j) - SUM(dzsnow(1:nsmax-1))
            END IF
            sice % DATA(i, pseudo_out) = &
              rho_snow_const * ds % DATA(i, pseudo_out)
            sliq % DATA(i, pseudo_out) = 0.0
!           Keep it frozen.
            tsnow % DATA(i, pseudo_out) = &
              MIN( soil_temp % DATA(i,1), tm-1.0 )
            rhosnow % DATA(i, pseudo_out) = rho_snow_const
!           Use a typical value for ice sheets. The point is indeed
!           not land ice, but will be close to land ice, so we assume
!           that surface snow will have aged in much the same way.
            rgrainl % DATA(i, pseudo_out) = rgrain_icesheet
          END DO
        END IF
      END DO
!
    END IF
!
  END DO
!
END IF


!-------------------------------------------------------------------
! Write out the multilayer fields and deallocate the storage.
!-------------------------------------------------------------------
CALL rcf_write_field( snowdepth,    hdr_out, decomp_rcf_output )
CALL rcf_write_field( snowbulkdens, hdr_out, decomp_rcf_output )
CALL rcf_write_field( nsnow,        hdr_out, decomp_rcf_output )
CALL rcf_write_field( ds,           hdr_out, decomp_rcf_output )
CALL rcf_write_field( sice,         hdr_out, decomp_rcf_output )
CALL rcf_write_field( sliq,         hdr_out, decomp_rcf_output )
CALL rcf_write_field( tsnow,        hdr_out, decomp_rcf_output )
CALL rcf_write_field( rhosnow,      hdr_out, decomp_rcf_output )
CALL rcf_write_field( rgrainl,      hdr_out, decomp_rcf_output )

IF(l_ml_snow_in) THEN
  DEALLOCATE( nsnow_in )
  DEALLOCATE( ds_in )
  DEALLOCATE( sice_in )
  DEALLOCATE( sliq_in )
  DEALLOCATE( tsnow_in )
  DEALLOCATE( rhosnow_in )
  DEALLOCATE( rgrainl_in )
  CALL rcf_dealloc_field( field_middle )
  CALL rcf_dealloc_field( snow_tile )
  IF (can_model == 4) CALL rcf_dealloc_field( snow_grnd )
END IF

CALL rcf_dealloc_field( soil_temp )
CALL rcf_dealloc_field( nsnow )
CALL rcf_dealloc_field( ds )
CALL rcf_dealloc_field( sice )
CALL rcf_dealloc_field( sliq )
CALL rcf_dealloc_field( tsnow )
CALL rcf_dealloc_field( rhosnow )
CALL rcf_dealloc_field( rgrainl )
CALL rcf_dealloc_field( snowdepth )
CALL rcf_dealloc_field( snowbulkdens )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_ml_snowpack
END MODULE rcf_ml_snowpack_mod
