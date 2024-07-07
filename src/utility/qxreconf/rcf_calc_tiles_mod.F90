! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reconfiguration for different numbers or types of tiles.

MODULE rcf_calc_tiles_mod

!  Subroutine rcf_calc_tiles

! Description:
!   This routine either initialises multiple tiles from
!   single-tile input, or matches tile ids between multiple-tile
!   input and output, depending on what is required.

! Method:
!   1-to-N, M-to-N and N-to-1.
!   The 1-to-N case calls routines which were previously
!   called directly from rcf_field_calcs.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_TILES_MOD'

CONTAINS

SUBROUTINE rcf_calc_tiles( fields_in, field_count_in, hdr_in,    &
                           fields_out, field_count_out, hdr_out, &
                           tiled_field_out, field_stashcode,     &
                           data_source )

USE rcf_umhead_mod, ONLY:                  &
    um_header_type

USE ereport_mod, ONLY:                     &
    ereport

USE umPrintMgr, ONLY:                      &
   umPrint,                                 &
   umMessage,                               &
   printstatus,                             &
   prstatus_diag

USE um_parcore, ONLY:                      &
    mype

USE rcf_field_type_mod, ONLY:              &
    field_type

USE rcf_data_source_mod, ONLY: &
    data_source_type,          &
    already_processed

USE items_nml_mod, ONLY: &
    Field_Dependent_Calcs

USE rcf_locate_mod, ONLY:                  &
    rcf_locate

USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_out

USE rcf_alloc_field_mod, ONLY:             &
    rcf_alloc_field,                        &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY:              &
    rcf_read_field

USE rcf_init_field_on_tiles_mod, ONLY:     &
    rcf_init_field_on_tiles

USE rcf_init_canopy_water_mod, ONLY:       &
    rcf_init_canopy_water

USE rcf_init_tile_snow_mod, ONLY:          &
    rcf_init_tile_snow

USE rcf_init_tile_t_mod, ONLY:             &
    rcf_init_tile_t

USE rcf_init_snowdep_mod, ONLY:       &
    rcf_init_snowdep

USE rcf_init_snow_bk_dens_mod, ONLY:       &
    rcf_init_snow_bk_dens

USE rcf_set_dominant_tile_mod, ONLY:       &
    rcf_set_dominant_tile, i_dominant

USE rcf_calc_tile_map_mod, ONLY:           &
    rcf_calc_tile_map

USE decomp_params, ONLY:                   &
    decomp_rcf_input,                      &
    decomp_rcf_output

USE um_stashcode_mod, ONLY:                &
    stashcode_prog_sec,                    &
    stashcode_tstar_tile,                  &
    stashcode_tsurf_elev_surft,            &
    stashcode_can_water_tile,              &
    stashcode_snow_tile,                   &
    stashcode_snow_grnd,                   &
    stashcode_catch_snow,                  &
    stashcode_rgrain,                      &
    stashcode_snowpack_bk_dens,            &
    stashcode_nsnow_layrs_tiles,           &
    stashcode_snow_laythk_tiles,           &
    stashcode_snow_ice_tile,               &
    stashcode_snow_liq_tile,               &
    stashcode_snow_t_tile,                 &
    stashcode_snow_laydns_tiles,           &
    stashcode_snow_grnsiz_tiles,           &
    stashcode_snowdep_grd_tile

USE rcf_nlist_recon_technical_mod, ONLY: &
    l_rcf_init_flexi

USE errormessagelength_mod, ONLY: errormessagelength

USE nlsizes_namelist_mod, ONLY: ntiles

USE jules_snow_mod, ONLY: nsmax

USE missing_data_mod, ONLY: imdi

USE Rcf_Interp_Weights_Mod, ONLY: h_int_active

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type), POINTER           :: fields_in(:),fields_out(:)
TYPE( field_type ), INTENT(INOUT)    :: tiled_field_out
TYPE( data_source_type ), POINTER    :: data_source( : )
TYPE( um_header_type ), INTENT(IN)   :: hdr_in, hdr_out
INTEGER, INTENT(IN)                  :: field_count_in, field_count_out, &
                                        field_stashcode

! Local vars.

TYPE( field_type ), POINTER  :: tiled_field_in

INTEGER :: ntiles_in ! Number of tiles in input dump
INTEGER :: nsmax_in  ! Maximum number of snow layers in input dump

LOGICAL :: l_pseudo_ok             ! logical to record a match in level ids
LOGICAL :: l_init_tile_t_zerofrac

INTEGER                :: i, j ! Loopers
INTEGER                :: pos  ! field position (DIFFERENT from dump_pos)
INTEGER                :: posd ! dependent field position
INTEGER                :: snow_tile_pos

INTEGER                :: errorstatus
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CALC_TILES'
CHARACTER (LEN=errormessagelength) :: cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
  WRITE(umMessage,'(a,i3)') 'Populating tile data for field ',    &
     field_stashcode
  CALL umPrint(umMessage,src='rcf_calc_tiles_mod')
END IF

!----------------------------------------------------------------------
! Find input field.
!----------------------------------------------------------------------
CALL rcf_locate( stashcode_prog_sec, field_stashcode,  &
   fields_in, field_count_in, pos, zero_ok_arg = .TRUE. )

IF ( pos /= 0 ) THEN
  tiled_field_in => fields_in(pos)
ELSE ! field_stashcode is not in input dump
  SELECT CASE ( field_stashcode )

  CASE ( stashcode_tstar_tile, &
         stashcode_tsurf_elev_surft, &
     stashcode_snowdep_grd_tile, &
     stashcode_snowpack_bk_dens, &
     stashcode_snow_tile )
    ! These are not expected to be in input dump; rcf_init routines use another
    !   stashcode to initialise
  CASE DEFAULT
    ! Shouldn't have got as far as this, but fatal error in case it has
    errorstatus = 6173
    WRITE(cmessage, '(A, I3, I5)') &
       'rcf_locate did not find field with STASHcode ', &
       stashcode_prog_sec, field_stashcode
    CALL ereport( routinename, errorstatus, cmessage )
  END SELECT
END IF

! Set ntiles for input dump. Required as ML snow fields have ntiles*nsmax
! levels so aggregate tile will have nsmax levels not 1

SELECT CASE ( field_stashcode )

CASE ( stashcode_catch_snow,        &
       stashcode_snow_grnd,         &
       stashcode_rgrain,            &
       stashcode_can_water_tile,    &
       stashcode_tstar_tile,        &
       stashcode_tsurf_elev_surft,  &
       stashcode_snow_tile,         &
       stashcode_snowdep_grd_tile,  &
       stashcode_snowpack_bk_dens,  &
       stashcode_nsnow_layrs_tiles )

  IF ( pos == 0 ) THEN
    ntiles_in  = imdi
  ELSE
    ntiles_in  = tiled_field_in % levels
  END IF

CASE ( stashcode_snow_laythk_tiles, &
       stashcode_snow_ice_tile,     &
       stashcode_snow_liq_tile,     &
       stashcode_snow_t_tile,       &
       stashcode_snow_laydns_tiles, &
       stashcode_snow_grnsiz_tiles )
  ! ML snow fields

  IF ( pos == 0 ) THEN
    ntiles_in  = imdi
  ELSE
    ! Check that nsmax is the same in the output dump as in the 
    ! input dump. If not, the snowpack will need to be relayered.
    CALL rcf_locate( stashcode_prog_sec, stashcode_snow_tile,  &
       fields_in, field_count_in, snow_tile_pos )
    nsmax_in = tiled_field_in % levels / fields_in(snow_tile_pos) % levels
    IF ( nsmax_in /= nsmax ) THEN
      data_source( pos ) % source = Field_Dependent_Calcs
    END IF
    ntiles_in  = tiled_field_in % levels / nsmax_in
  END IF

CASE DEFAULT
  errorstatus = 6174
  WRITE (cmessage, '(A, I5)')                                &
     'STASHcode not recognised for setting ntiles_in: ', &
     tiled_field_out % stashmaster % item
  CALL ereport( routinename, errorstatus, cmessage )
END SELECT

!----------------------------------------------------------------------
! 1 to N (expansion) including N=1 or general initialisation of fields
!    with introduced tiles (l_rcf_init_flexi)
!----------------------------------------------------------------------
IF ( ( ( ntiles_in == 1 .OR. pos == 0 ) .AND. ntiles >= 1 ) &
   .OR. l_rcf_init_flexi ) THEN

  SELECT CASE( field_stashcode )

  CASE ( stashcode_catch_snow,        &
         stashcode_snow_grnd,         &
         stashcode_rgrain,            &
         stashcode_nsnow_layrs_tiles, &
         stashcode_snow_laythk_tiles, &
         stashcode_snow_ice_tile,     &
         stashcode_snow_liq_tile,     &
         stashcode_snow_t_tile,       &
         stashcode_snow_laydns_tiles, &
         stashcode_snow_grnsiz_tiles )
    CALL rcf_init_field_on_tiles( fields_in, field_count_in,        &
                                hdr_in, ntiles_in, tiled_field_out, &
                                tiled_field_out % stashmaster % item )

  CASE ( stashcode_can_water_tile )
    CALL rcf_init_canopy_water( fields_in, field_count_in,       &
                                hdr_in, tiled_field_out )

  CASE ( stashcode_tstar_tile, stashcode_tsurf_elev_surft  )
    l_init_tile_t_zerofrac = .FALSE.
    CALL rcf_init_tile_t( fields_out, field_count_out, hdr_out,         &
                          decomp_rcf_output, local_lsm_out,             &
                          tiled_field_out,l_init_tile_t_zerofrac )

  CASE ( stashcode_snow_tile )
    CALL rcf_init_tile_snow( fields_out, field_count_out,        &
                             hdr_out, tiled_field_out )

  CASE ( stashcode_snowdep_grd_tile )
    CALL rcf_init_snowdep( fields_in, field_count_in, hdr_in,  &
                           fields_out, field_count_out, hdr_out,   &
                           data_source, tiled_field_out )

  CASE ( stashcode_snowpack_bk_dens )
    CALL rcf_init_snow_bk_dens( fields_in, field_count_in, hdr_in,  &
                           fields_out, field_count_out, hdr_out,   &
                           tiled_field_out )

  CASE DEFAULT
    errorstatus = 6174
    WRITE (cmessage, '(A, I5)')                                     &
          '1-to-N: No Field Calculations specified for STASHcode ', &
          tiled_field_out % stashmaster % item
    CALL ereport( routinename, errorstatus, cmessage )

  END SELECT

  !----------------------------------------------------------------------
  ! M to N or M to m (different numbers or types)
  !----------------------------------------------------------------------
ELSE IF ( ntiles_in > 1 .AND. ntiles > 1 ) THEN

  SELECT CASE( field_stashcode )
    CASE ( stashcode_snowdep_grd_tile )

      ! If initializing multilayer snow fields from a dump which does
      ! contain total snow depth, but not other multilayer fields
      ! (ie. starting from a zero-layer start dump), the total
      ! snowdepth is ignored and it is recalculated from the snow mass
      ! to ensure consistency.
      ! Check that the multilayer fields are being initialized by
      ! looking for the absence of thickness of individual snow layers.

      CALL rcf_locate( stashcode_prog_sec, stashcode_snow_laythk_tiles, &
                       fields_in, field_count_in, posd,zero_ok_arg = .TRUE. )
      IF ( (posd == 0) .AND. (nsmax > 0) ) THEN
        ! Initialize from the masses of snow.
        CALL rcf_init_snowdep( fields_in, field_count_in, hdr_in,       &
                               fields_out, field_count_out, hdr_out,    &
                               data_source, tiled_field_out )
      ELSE
        ! Set the fields from the input in the normal way.
        CALL rcf_calc_tile_map( tiled_field_in, hdr_in, ntiles_in,      &
                                tiled_field_out, field_stashcode )
        CALL rcf_init_field_on_tiles( fields_in, field_count_in,        &
                                      hdr_in, ntiles_in,                &
                                      tiled_field_out,                  &
                                      tiled_field_out % stashmaster % item )
      END IF
  CASE DEFAULT
!   Set the field from the input.
    CALL rcf_calc_tile_map( tiled_field_in, hdr_in, ntiles_in,          &
                            tiled_field_out, field_stashcode )
    CALL rcf_init_field_on_tiles( fields_in, field_count_in,            &
                                  hdr_in, ntiles_in, tiled_field_out,   &
                                  tiled_field_out % stashmaster % item )
  END SELECT

ELSE IF ( ntiles_in > 1 .AND. ntiles == 1 ) THEN

  !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !   Note on Canopy Snow: If the dominant tile in the input has a 
  !   canopy snow model, there will be a ground store. Although it is
  !   incompatible with conservation, we ignore the ground store here.
  !   If it were included we would end up transferring snow from the
  !   ground to the canopy, massively increasing the albedo late in
  !   the season: this is probably more serious than a loss of snow
  !   mass.
  !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  SELECT CASE( field_stashcode )

  CASE ( stashcode_can_water_tile )
    CALL rcf_init_canopy_water( fields_in, field_count_in,       &
                                hdr_in, tiled_field_out )

  CASE ( stashcode_tstar_tile, stashcode_tsurf_elev_surft )
    l_init_tile_t_zerofrac = .FALSE.
    CALL rcf_init_tile_t( fields_out, field_count_out, hdr_out,  &
                          decomp_rcf_output, local_lsm_out,      &
                          tiled_field_out, l_init_tile_t_zerofrac)

  CASE ( stashcode_snow_tile )
    CALL rcf_init_tile_snow( fields_out, field_count_out,        &
                             hdr_out, tiled_field_out )

  CASE ( stashcode_catch_snow,        &
         stashcode_snow_grnd,         &
         stashcode_rgrain,            &
         stashcode_snowpack_bk_dens,  &
         stashcode_nsnow_layrs_tiles, &
         stashcode_snow_laythk_tiles, &
         stashcode_snow_ice_tile,     &
         stashcode_snow_liq_tile,     &
         stashcode_snow_t_tile,       &
         stashcode_snow_laydns_tiles, &
         stashcode_snow_grnsiz_tiles, &
         stashcode_snowdep_grd_tile )

! While rcf_init_field_on_tiles is used to set the DATA field from the dominant
! tile. rcf_set_dominant_tile selects the dominant tile from the input grid, so
! any ancillary information or change in grid size is currently not supported.
  IF ( h_int_active ) THEN
    errorstatus = 6174
    WRITE (cmessage, '(A, I5)')                               &
     'Try two-step recon: ' //                                &
     'Change in grid for n-to-1 not currently supported for STASHcode ', &
     tiled_field_out % stashmaster % item
    CALL ereport( routinename, errorstatus, cmessage )
  END IF

  IF ( .NOT. ALLOCATED(i_dominant) ) THEN
    ALLOCATE( i_dominant(tiled_field_in % level_size) )

    !     Find the dominant tile at each point.
    CALL Rcf_Set_Dominant_Tile( fields_in, field_count_in, hdr_in, &
                                fields_out, field_count_out, hdr_out )
  END IF

  CALL rcf_init_field_on_tiles( fields_in, field_count_in,        &
                              hdr_in, ntiles_in, tiled_field_out, &
                              tiled_field_out % stashmaster % item )

  CASE DEFAULT
    errorstatus = 6174
    WRITE (cmessage, '(A, I5)')                                  &
       'N-to-1: No Field Calculations specified for STASHcode ', &
       tiled_field_out % stashmaster % item
    CALL ereport( routinename, errorstatus, cmessage )

  END SELECT

END IF

! These fields need to be processed before multilayer snow fields are
! set.
SELECT CASE ( field_stashcode )
CASE ( stashcode_snow_tile,  &
       stashcode_snow_grnd )
  CALL rcf_locate( stashcode_prog_sec, field_stashcode,  &
     fields_out, field_count_out, pos )
  IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
    WRITE(umMessage,'(a,i5)') &
       'Setting data_source to already_processed: stashcode', &
       field_stashcode
    CALL umPrint(umMessage,src='rcf_calc_tiles_mod')
  END IF
  data_source( pos ) % source = already_processed
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_calc_tiles
END MODULE rcf_calc_tiles_mod

