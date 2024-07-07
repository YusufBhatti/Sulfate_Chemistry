! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Checks the required land-surface types are in the input dump.

MODULE rcf_compare_tiles_mod

!  Subroutine rcf_compare_tiles

! Description:
!   A mismatch in number of tiles or tile types results in
!   the setting of a logical flag, l_match, to FALSE.

! Method:
!   The number of tiles input and output is first compared. If
!   they are the same, a second test compares the tile pseudo-level
!   ids at identical positions. Failure of either test results
!   in the l_match logical changing from TRUE to FALSE.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_COMPARE_TILES_MOD'

CONTAINS

SUBROUTINE rcf_compare_tiles( fields_in, field_count_in, hdr_in, &
   tiled_field_out, hdr_out, field_stashcode, l_match )

USE rcf_umhead_mod, ONLY:                  &
   um_header_type

USE ereport_mod, ONLY:                     &
    ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE umPrintMgr, ONLY:                      &
   umPrint,                                 &
   umMessage,                               &
   printstatus,                             &
   prstatus_diag

USE um_parcore, ONLY:                      &
   mype

USE rcf_field_type_mod, ONLY:              &
   field_type

USE rcf_locate_mod, ONLY:                  &
   rcf_locate

USE rcf_alloc_field_mod, ONLY:             &
   rcf_alloc_field,                         &
   rcf_dealloc_field

USE rcf_read_field_mod, ONLY:              &
   rcf_read_field

USE decomp_params, ONLY:                   &
   decomp_rcf_input

USE um_stashcode_mod, ONLY:              &
   stashcode_prog_sec

USE jules_surface_types_mod, ONLY: tile_map_ids

USE lookup_addresses, ONLY: lbplev

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type), POINTER         :: fields_in(:)
TYPE( field_type ), INTENT(IN)     :: tiled_field_out
TYPE( um_header_type ), INTENT(IN) :: hdr_in, hdr_out
INTEGER, INTENT(IN)                :: field_count_in, field_stashcode
LOGICAL, INTENT(OUT)               :: l_match

! Local vars.

INTEGER  :: tile_id_in, tile_id_out  ! tile ids
INTEGER  :: i    ! Looper
INTEGER  :: pos  ! field position (DIFFERENT from dump_pos)
INTEGER  :: errorstatus

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_COMPARE_TILES'
CHARACTER (LEN=errormessagelength) :: cmessage

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
  WRITE(umMessage,'(a,i3)') 'Comparing tile data for field ',         &
     field_stashcode
  CALL umPrint(umMessage,src='rcf_compare_tiles_mod')
END IF

!----------------------------------------------------------------------
! Find input field.
!----------------------------------------------------------------------
! A call from rcf_snowstores exists to check that stashcode_snow_tile and
! stashcode_snow_grnd have the correct configuration before proceeding.
! However, it is conceivable that a situation may exist when they are not in
! the input dump. This will result in a fatal error. An option may need to
! be introduced to check against e.g. stashcode_frac_surf_type instead.
CALL rcf_locate( stashcode_prog_sec, field_stashcode,  &
   fields_in, field_count_in, pos )

!----------------------------------------------------------------------
! Starting value which testing may change.
!----------------------------------------------------------------------
l_match = .TRUE.

!----------------------------------------------------------------------
! Gross comparison of number of levels.
!----------------------------------------------------------------------
IF (tiled_field_out % levels /= fields_in(pos) % levels) THEN

  l_match=.FALSE.

  IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
    WRITE(umMessage,'(a,i3)') 'Different number of pseudo levels for ', &
       field_stashcode
    CALL umPrint(umMessage,src='rcf_compare_tiles_mod')
  END IF

END IF

!----------------------------------------------------------------------
! Only do the next test if the last test passed.
!----------------------------------------------------------------------
IF (l_match) THEN ! test #2

  !----------------------------------------------------------------------
  ! Futher comparison by tile / pseudo-level id
  !----------------------------------------------------------------------
  DO i = 1, tiled_field_out % levels

    !-----------------------------------------------------
    ! pseudo-level ids stored in the LOOKUP
    ! in the LBPLEV row
    ! in positions DUMP_POS : DUMP_POS+LEVELS-1
    !-----------------------------------------------------
    tile_id_in  = hdr_in  % lookup(lbplev,fields_in(pos)  % dump_pos+i-1)
    tile_id_out = hdr_out % lookup(lbplev,tiled_field_out % dump_pos+i-1)

    !---------------------------------------------------------
    ! Do the IN and OUT fields match *at the same position* ?
    !---------------------------------------------------------
    IF (tile_id_in /= tile_id_out) THEN
      l_match = .FALSE.
      ! Check for ML snow reindexing of first layer (#1044)
      IF ( (tile_id_out - 1)/1000 == tile_id_in ) l_match = .TRUE.
    END IF

  END DO ! i

  IF (.NOT. l_match) THEN

    IF (mype == 0 .AND. printstatus >= prstatus_diag ) THEN
      WRITE(umMessage,'(a,i3)') 'Non-identical pseudo levels for ', &
         field_stashcode
      CALL umPrint(umMessage,src='rcf_compare_tiles_mod')
    END IF
  END IF

END IF ! test #2

IF ( MAXVAL(tile_map_ids) > 0 .AND. l_match ) THEN
  errorstatus = -10
  WRITE(cmessage, '(A,I6)') &
       'tile_map_ids is specified, but it will be ignored as there is no '// &
       'difference in surface type configuration for ', field_stashcode
  CALL ereport ( routinename, errorstatus, cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_compare_tiles
END MODULE rcf_compare_tiles_mod

