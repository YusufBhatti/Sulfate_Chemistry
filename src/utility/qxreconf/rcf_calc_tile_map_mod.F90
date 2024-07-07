! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reconfiguration for different numbers or types of tiles.

MODULE rcf_calc_tile_map_mod

!  Subroutine rcf_calc_tile_map

! Description:
!   This routine calculates the tile mapping array, if it is not already set,
!   and performs checks.

! Method:
!   This routine either calculates a tile mapping array given the
!   configuration in the input and output dumps or uses a user specified
!   array. Checks the mapping array whenever there is a change in input
!   surface type configuration. If there is no change in the input field
!   surface type configuration between calls then nothing is done.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_TILE_MAP_MOD'

CONTAINS

SUBROUTINE rcf_calc_tile_map( tiled_field_in, hdr_in, ntiles_in,              &
                              tiled_field_out, field_stashcode )

USE rcf_umhead_mod, ONLY:                  &
    um_header_type

USE ereport_mod, ONLY:                     &
    ereport

USE umPrintMgr, ONLY:                      &
   umPrint,                                &
   umMessage,                              &
   printstatus,                            &
   prstatus_diag,                          &
   maxLineLen

USE um_parcore, ONLY:                      &
    mype

USE rcf_field_type_mod, ONLY:              &
    field_type


USE um_stashcode_mod, ONLY:                &
              stashcode_snow_laythk_tiles, &
              stashcode_snow_ice_tile,     &
              stashcode_snow_liq_tile,     &
              stashcode_snow_t_tile,       &
              stashcode_snow_laydns_tiles, &
              stashcode_snow_grnsiz_tiles

USE lookup_addresses, ONLY: lbplev

USE errormessagelength_mod, ONLY: errormessagelength

USE land_tile_ids, ONLY: surface_type_ids, ml_snow_type_ids, tile_ids_in

USE max_dimensions, ONLY: ntype_max, snow_layers_max

USE nlsizes_namelist_mod, ONLY: ntiles

USE jules_snow_mod, ONLY: nsmax

USE jules_surface_types_mod, ONLY: tile_map_ids, tile_map_pslevs

USE missing_data_mod, ONLY: imdi

USE rcf_nlist_recon_technical_mod, ONLY: &
    l_rcf_init_flexi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: tiled_field_in
TYPE( field_type ), INTENT(IN)       :: tiled_field_out
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
INTEGER, INTENT(IN)                  :: ntiles_in ! ntiles in input dump
INTEGER, INTENT(IN)                  :: field_stashcode

! Local vars.

! tile ids
INTEGER, ALLOCATABLE :: tile_ids_in_tmp(:)
             ! Working array of tile IDs in the input header
INTEGER, ALLOCATABLE :: tile_ids_out_hdr(:)
             ! Working array of expected tile IDs of output configuration
INTEGER, ALLOCATABLE :: tile_ids_out_calc(:)
             ! Working array of calculated tile IDs of output configuration
INTEGER :: max_levels_calc,    & ! Number of tile IDs out sensible to calculate
           nsmax_in              ! Maximum number of snow levels in input dump

LOGICAL :: l_usr_tile_map = .FALSE., &
           l_t1044        = .FALSE. ! Indicates a pre-#1044 indexing in ML field

LOGICAL :: tile_ids_out_check ! Check tile IDs out are correct

INTEGER                :: i, j ! Loopers

INTEGER                :: errorstatus
! Stop umFormat overwriting umMessage length, based on format statement
! (A,n(1X,I6)) < maxLineLen
INTEGER, PARAMETER     :: maxLevelsPrint = ( maxLineLen - 50 ) / 7
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CALC_TILE_MAP'
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=20) :: umFormat
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
ALLOCATE(tile_ids_in_tmp(tiled_field_in % levels))
ALLOCATE(tile_ids_out_calc(tiled_field_out % levels))
ALLOCATE(tile_ids_out_hdr(tiled_field_out % levels))

IF ( tiled_field_out % levels == nsmax*ntiles ) THEN
  tile_ids_out_hdr(:) = ml_snow_type_ids(1:tiled_field_out % levels)
ELSE
  tile_ids_out_hdr(:) = surface_type_ids(1:tiled_field_out % levels)
END IF

! Read INPUT tile configuration from INPUT header
tile_ids_in_tmp(:) = hdr_in % lookup( lbplev, tiled_field_in % dump_pos:  &
   tiled_field_in % dump_pos + tiled_field_in % levels - 1 )
! Check INPUT tile configuration = previous INPUT tile configuration
!   (tile_ids_in) stored in land_tile_ids (initially imdi). If they are not
!   identical then update tile_ids_in, calculate tile_map_pslevs
!   and perform checks.

IF ( MAXVAL(tile_map_ids) > 0 ) THEN
  l_usr_tile_map = .TRUE.
  IF ( l_rcf_init_flexi ) THEN
    errorstatus = 5
    WRITE(cmessage, '(A,I6)') &
       'Both l_rcf_init_flexi = T and tile_map_ids is specified. ' //     &
       'These are mutually exclusive so please use one or the other.'
    CALL ereport ( routinename, errorstatus, cmessage )
  END IF
END IF

IF ( .NOT. compare_configs ( tile_ids_in_tmp(:), &
                             tile_ids_in(1:tiled_field_in % levels) ) ) THEN
  tile_ids_in(1:tiled_field_in % levels) = tile_ids_in_tmp(:)
  ! Should only have to calculate the mapping array once; snow levels map as
  !      pseudo = i + (j-1) * ntiles, where i = 1, ntiles and j = 2, nsmax
  ! Following IF block is only done on first call
  IF ( .NOT. l_usr_tile_map ) THEN
    IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
      WRITE(umMessage,'(A,I3)') 'Calculating tile mapping array for field ',  &
         field_stashcode
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
    END IF
    ! For robustness although currently there are no fields of size npft
    ! indicated.
    IF ( tiled_field_out % levels /= ntiles ) THEN
      SELECT CASE ( field_stashcode )
      CASE ( stashcode_snow_laythk_tiles )
        ! This will be the first field if the multilayer snow scheme has a
        ! difference in indexing (see #1044) or the number of snow layers
        ! differs between the input and output dumps, in which case the recon
        ! path will lead directly here from rcf_ml_snowpack. These scenarios
        ! have now been dealt with, provided that the actual surface type
        ! configuration is the same.
        IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
          CALL umPrint( routinename //                   &
             ': Checking tile_map_ids from rcf_ml_snowpack', &
             src='rcf_calc_tile_map_ids_mod')
        END IF
      CASE DEFAULT
        errorstatus = 10
        WRITE(cmessage, '(A,I6)') &
           'Routine is not compatible with size of first field /= ntiles', &
           field_stashcode
        CALL ereport ( routinename, errorstatus, cmessage )
      END SELECT
    END IF
    DO i = 1, ntiles ! Loop over output configuration
      ! Check required tile types are present in input dump and if so create
      ! tile map. If they are not present check for ML snow reindexing (#1044).
      IF ( COUNT( tile_ids_out_hdr(i) == tile_ids_in ) == 1 ) THEN
        DO j = 1, ntiles_in
          IF ( tile_ids_out_hdr(i) == tile_ids_in(j) ) THEN
            tile_map_pslevs(i) = j
            EXIT
          END IF
        END DO
      ELSE IF ( ( tile_ids_out_hdr(i) - 1 ) / 1000 == tile_ids_in(i) ) THEN
        l_t1044            = .TRUE.
        tile_map_pslevs(i) = i
      ELSE
        errorstatus = 20
        WRITE(cmessage, '(A,I6,A,I6)') &
           'Consider l_rcf_init_flexi = T or set tile_map_ids: '// &
           'Cannot map tiled field ', field_stashcode,         &
           ', tile id ', tile_ids_out_hdr(i)
        CALL ereport ( routinename, errorstatus, cmessage )
      END IF
    END DO
  ELSE IF ( MAXVAL( tile_map_pslevs ) < 0 ) THEN
    ! Routine will need to be altered in future if the need arises to deal
    ! with the first field not having SIZE = ntiles. Fail if this occurs.
    SELECT CASE ( field_stashcode )
    CASE ( stashcode_snow_laythk_tiles )
      ! This ML snow field will be the first if the multilayer snow scheme has a
      ! difference in indexing (see #1044) or the number of snow layers differs
      ! between the input and output dumps, in which case the recon path will
      ! lead directly here from rcf_ml_snowpack. An exception is made for ML
      ! snow.
      ! Do nothing for the minute. In either case that I can think of
      ! tile_ids_in(1:ntiles_in) = tile IDs not ML snow indices
      !   1. pre um:#578 has 1-27 for an old 9T dump i.e. 1-9 in first layer
      !   2. pre um:#1044 has e.g. 1,2,3,4,5,7,8,9,601,602 for ps37 dump
      IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
        CALL umPrint( routinename //                   &
           ': Checking tile_map_ids from rcf_ml_snowpack', &
           src='rcf_calc_tile_map_ids_mod')
      END IF
    CASE DEFAULT
      IF ( tiled_field_out % levels /= ntiles ) THEN
        errorstatus = 30
        WRITE(cmessage, '(A,I6)') &
           'Routine is not compatible with size of first field /= ntiles', &
           field_stashcode
        CALL ereport ( routinename, errorstatus, cmessage )
      END IF
    END SELECT
    ! Check that tile_map_ids is consistent with the surface types in the input
    ! dump.
    errorstatus = 0
    IF ( COUNT ( tile_map_ids(1:ntiles) == imdi ) > 0 ) errorstatus = 40
    ! Currently not called for a npft field, but in the future if the first
    ! field has size npft then it will have to be tweaked as currently, if
    ! the check above is removed it will fail here when checking that
    ! tile_map_ids is consistent.
    DO i = 1, ntiles ! Loop over output configuration
      IF ( COUNT( tile_map_ids(i) == tile_ids_in ) == 0 ) errorstatus = 50
    END DO
    IF ( errorstatus > 0 ) THEN
      WRITE(umFormat,'(A,I3,A)') '(A,',                                      &
         MIN( maxLevelsPrint, tiled_field_in % levels ),'(1X,I6))'
      WRITE(umMessage,umFormat) 'Tile IDs in  ',                             &
         tile_ids_in(1:MIN( maxLevelsPrint, tiled_field_in % levels ))
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
      WRITE(umFormat,'(A,I3,A)') '(A,',                                      &
         MIN( maxLevelsPrint, tiled_field_out % levels ),'(1X,I6))'
      WRITE(umMessage,umFormat) 'Tile IDs out ',                             &
         tile_ids_out_hdr(1:MIN( maxLevelsPrint, tiled_field_out % levels ))
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
      WRITE(umFormat,'(A,I3,A)') '(A,',ntiles,'(1X,I6))'
      WRITE(umMessage,umFormat) 'Mapping array', tile_map_ids(1:ntiles)
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
      WRITE(cmessage,'(A,I3)')                                               &
         'User specified mapping array inconsistent with tile '//            &
         'configuration ', field_stashcode
      CALL ereport( routinename, errorstatus, cmessage )
    END IF

    ! User specified tile_map_ids gives the input surface type ID that the
    ! output surface type should be initialised from, not which pseudo level
    ! it is in so we need to convert it
    DO i = 1, ntiles ! Loop over output configuration
      DO j = 1, ntiles_in
        IF ( tile_map_ids(i) == tile_ids_in(j) ) THEN
          tile_map_pslevs(i) = j
          EXIT
        END IF
      END DO
    END DO
  END IF

  IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
    WRITE(umMessage,'(A,I3)') 'Checking tile mapping array for field ', &
       field_stashcode
    CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
  END IF

  ! Check that the mapping array, tile_ids_in and tile_ids_out_hdr are
  ! consistent including the ML snow levels
  tile_ids_out_calc(:) = imdi
  ! No fields currently indicated that have size < ntiles, but in future field
  ! of size npft might be added so avoid out of bounds
  DO i = 1, MIN( ntiles, tiled_field_out % levels )
    tile_ids_out_calc(i) = tile_ids_in(tile_map_pslevs(i))
    SELECT CASE ( field_stashcode )
    CASE ( stashcode_snow_laythk_tiles, &
           stashcode_snow_ice_tile,     &
           stashcode_snow_liq_tile,     &
           stashcode_snow_t_tile,       &
           stashcode_snow_laydns_tiles, &
           stashcode_snow_grnsiz_tiles)
      IF ( l_t1044 ) tile_ids_out_calc(i) = tile_ids_out_calc(i) * 1000 + 1
      nsmax_in = tiled_field_in % levels / ntiles_in
      IF ( nsmax > 1 ) THEN
        ! Don't calculate tile IDs for more than nsmax_in as they won't exist
        DO j = 2, MIN( nsmax, nsmax_in )
          tile_ids_out_calc(i + (j-1) * ntiles) = &
             tile_ids_in(tile_map_pslevs(i) + (j-1) * ntiles_in)
        END DO
      END IF
    END SELECT
  END DO

  SELECT CASE ( field_stashcode )
  CASE ( stashcode_snow_laythk_tiles )
    ! tile_ids_out_calc won't exist for more than nsmax_in * ntiles
    max_levels_calc = ntiles * MIN( nsmax, nsmax_in )
    tile_ids_out_check = &
       compare_configs ( tile_ids_out_calc(1:max_levels_calc), &
                         tile_ids_out_hdr(1:max_levels_calc) )
  CASE DEFAULT
    tile_ids_out_check = &
       compare_configs ( tile_ids_out_calc(:), tile_ids_out_hdr(:) )
  END SELECT

  IF ( .NOT. tile_ids_out_check ) THEN
    WRITE(umFormat,'(A,I3,A)') '(A,',                                        &
       MIN( maxLevelsPrint, tiled_field_in % levels ), '(1X,I6))'
    WRITE(umMessage,umFormat) 'Tile IDs in                   ',              &
       tile_ids_in(1:MIN( maxLevelsPrint, tiled_field_in % levels ))
    CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
    WRITE(umFormat,'(A,I3,A)') '(A,',                                        &
       MIN( maxLevelsPrint, tiled_field_out % levels ), '(1X,I6))'
    WRITE(umMessage,umFormat) 'Tile IDs out                  ',              &
       tile_ids_out_hdr(1:MIN( maxLevelsPrint, tiled_field_out % levels ))
    CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
    WRITE(umMessage,umFormat) 'Tile IDs out (calculated)     ',              &
       tile_ids_out_calc(1:MIN( maxLevelsPrint, tiled_field_out % levels ))
    CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
    WRITE(umFormat,'(A,I3,A)') '(A,',ntiles,'(1X,I6))'
    IF ( l_usr_tile_map ) THEN
      errorstatus = -60
      WRITE(umMessage,umFormat) 'Mapping array (pseudo levels) ',            &
         tile_map_pslevs(1:ntiles)
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
      WRITE(cmessage,'(A,I3)')                                               &
         'Mapping array user specified; please check that it is correct! ',  &
         field_stashcode
      CALL ereport( routinename, errorstatus, cmessage )
    ELSE
      errorstatus = 70
      WRITE(umMessage,umFormat) 'Mapping array                 ',            &
         tile_map_pslevs(1:ntiles)
      CALL umPrint(umMessage,src='rcf_calc_tile_map_mod')
      WRITE(cmessage,'(A,I3)')                                               &
         'Mapping array inconsistent with tile configuration ',              &
         field_stashcode
      CALL ereport( routinename, errorstatus, cmessage )
    END IF
  END IF

END IF

DEALLOCATE(tile_ids_in_tmp)
DEALLOCATE(tile_ids_out_calc)
DEALLOCATE(tile_ids_out_hdr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_calc_tile_map


FUNCTION compare_configs ( config1, config2 )
IMPLICIT NONE
INTEGER :: config1(:), config2(:)
INTEGER :: i, npslevs ! Loop counters
LOGICAL :: compare_configs

compare_configs = .TRUE.
npslevs = SIZE( config1 )
IF ( npslevs /= SIZE( config2 ) ) compare_configs = .FALSE.
IF ( compare_configs ) THEN
  DO i = 1, npslevs
    IF ( config1(i) /= config2(i) ) compare_configs = .FALSE.
  END DO
END IF

END FUNCTION compare_configs


END MODULE rcf_calc_tile_map_mod

