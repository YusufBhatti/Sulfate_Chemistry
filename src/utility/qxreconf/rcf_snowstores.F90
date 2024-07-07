! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Processes snow stores.

MODULE rcf_snowstores_mod
IMPLICIT NONE

! Description:
!     This subroutine processes the snow stores,
!     creating fields that are consistent with JULES.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SNOWSTORES_MOD'

CONTAINS

SUBROUTINE rcf_snowstores( fields_in, field_count_in, hdr_in,    &
                             fields_out, field_count_out, hdr_out, &
                             data_source )

USE ereport_mod, ONLY: &
    ereport

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_grid_type_mod, ONLY: &
    grid_type

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE um_stashcode_mod, ONLY:    &
    stashcode_prog_sec          , &
    stashcode_snow_tile         , &
    stashcode_snow_grnd

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_data_source_mod, ONLY: &
    data_source_type
 
USE items_nml_mod, ONLY: &
    field_calcs

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parcore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE rcf_lsm_mod, ONLY: &
    local_land_out

USE nlsizes_namelist_mod, ONLY: &
    ntiles

USE jules_snow_mod, ONLY: cansnowtile

USE jules_vegetation_mod, ONLY: can_model

USE rcf_nlist_recon_science_mod, ONLY: &
    l_canopy_snow_throughfall

USE rcf_compare_tiles_mod, ONLY:                                          &
    rcf_compare_tiles

USE errormessagelength_mod, ONLY: errormessagelength
 
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_in(:), fields_out(:)
TYPE( data_source_type ), POINTER :: data_source( : )
TYPE( um_header_type), INTENT(IN) :: hdr_in, hdr_out
INTEGER, INTENT(IN)               :: field_count_in, field_count_out

! Internal variables
TYPE( field_type ), POINTER  :: snow_tile
TYPE( field_type ), POINTER  :: snow_grnd

REAL         :: total_snow_ground( local_land_out, ntiles )
REAL         :: snow_ground_store( local_land_out, ntiles )

INTEGER      ::  errorstatus = 0
INTEGER      ::  pos         ! position in array
INTEGER      ::  i,n         ! loop counters

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_SNOWSTORES'
CHARACTER (LEN=errormessagelength) :: Cmessage

LOGICAL :: l_match

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE(umMessage,'(a,i4)') &
      ' Processing snow stores '
  CALL umPrint(umMessage,src='rcf_Snowstores')
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------

! snow on tile
CALL rcf_locate( stashcode_prog_sec, stashcode_snow_tile,                   &
           fields_out, field_count_out, pos)
snow_tile => fields_out(pos)
CALL rcf_alloc_field( snow_tile )
CALL rcf_read_field( snow_tile, hdr_out, decomp_rcf_output )
! Check surface type configuration if snow_tile has not already been
! processed i.e. if it is still field_calcs
IF ( data_source( pos ) % source == field_calcs ) THEN
  CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
     fields_out(pos), hdr_out, stashcode_snow_tile, l_match )
  IF ( .NOT. l_match ) THEN
    errorstatus = 1
    WRITE (cmessage, '(A,I5)') &
       'Input field /= jules_surface_types configuration and has not yet '// &
       'passed through field_calcs: stashcode', stashcode_snow_tile
    CALL ereport( routinename, errorstatus, cmessage )
  END IF
END IF

IF (can_model == 4) THEN
  ! snow under canopy
  CALL rcf_locate( stashcode_prog_sec, stashcode_snow_grnd,                 &
             fields_out, field_count_out, pos)
  snow_grnd => fields_out(pos)
  CALL rcf_alloc_field( snow_grnd )
  CALL rcf_read_field( snow_grnd, hdr_out, decomp_rcf_output )
  ! Check surface type configuration if snow_grnd has not already been
  ! processed i.e. if it is still field_calcs
  IF ( data_source( pos ) % source == field_calcs ) THEN
    CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
       fields_out(pos), hdr_out, stashcode_snow_grnd, l_match )
    IF ( .NOT. l_match ) THEN
      errorstatus = 2
      WRITE (cmessage, '(A,I5)') &
         'Input field /= jules_surface_types configuration and has not '// &
         'yet passed through field_calcs: stashcode', stashcode_snow_grnd
      CALL ereport( routinename, errorstatus, cmessage )
    END IF
  END IF
END IF

!-----------------------------------------------------------------------------
! Do sanity checks because these fields are used interchangeably.
!-----------------------------------------------------------------------------
IF ( snow_tile % level_size /= local_land_out ) THEN
  errorstatus = 3
  WRITE (cmessage, '(A, I8, A, I8)')                &
        'Field sizes do not agree: level_size: ',   &
        snow_tile % level_size,                     &
        ' local_land_out: ',                            &
        local_land_out
  CALL ereport( routinename, errorstatus, cmessage )
END IF
IF ( snow_tile % levels /= ntiles ) THEN
  errorstatus = 4
  WRITE (cmessage, '(A, I3, A, I3)')                &
        'Field sizes do not agree: levels: ',       &
        snow_tile % levels,                         &
        ' ntiles: ',                                &
        ntiles
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!-----------------------------------------------------------------------------
! Record the snow amount on the ground,
! putting all snow onto the ground and zeroing canopy snow if required.
!-----------------------------------------------------------------------------
IF (can_model == 4) THEN
  snow_ground_store = snow_grnd % DATA
ELSE
  snow_ground_store = 0.0
END IF

IF (l_canopy_snow_throughfall) THEN

  IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
    CALL umPrint( ' Moving snow from canopy to ground', &
        src='rcf_snowstores')
  END IF

  DO n=1,ntiles
    IF ( cansnowtile(n) .AND. (can_model == 4) ) THEN
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % DATA(i,n), 0.0 ) &
                               + MAX( snow_ground_store(i,n), 0.0 )
      END DO
    ELSE
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % DATA(i,n), 0.0 )
      END DO
    END IF
  END DO
  ! Initialise stores to zero.
  snow_ground_store = 0.0
  snow_tile % DATA  = 0.0

  DO n=1,ntiles
    IF ( canSnowTile(n) .AND. (can_model == 4) ) THEN
      DO i=1,local_land_out
        snow_ground_store(i,n) = total_snow_ground(i,n)
      END DO
    ELSE
      DO i=1,local_land_out
        snow_tile % DATA(i,n) = total_snow_ground(i,n)
      END DO
    END IF
  END DO

ELSE ! no throughfall

  DO n=1,ntiles
    IF ( cansnowtile(n) .AND. (can_model == 4) ) THEN
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX( snow_ground_store(i,n), 0.0 )
      END DO
    ELSE
      DO i=1,local_land_out
        total_snow_ground(i,n) = MAX(  snow_tile % DATA(i,n), 0.0 )
      END DO
    END IF
  END DO

END IF ! throughfall test


!----------------------------------------------------------------------
! Write out the calculated snow fields.
!----------------------------------------------------------------------
CALL rcf_write_field( snow_tile , hdr_out, decomp_rcf_output )
IF (can_model == 4) THEN
  snow_grnd % DATA = snow_ground_store
  CALL rcf_write_field( snow_grnd , hdr_out, decomp_rcf_output )
END IF

!----------------------------------------------------------------------
! Clear up dynamic memory used.
!----------------------------------------------------------------------
CALL rcf_dealloc_field( snow_tile )
IF (can_model == 4) THEN
  CALL rcf_dealloc_field( snow_grnd )
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_snowstores
END MODULE rcf_snowstores_mod
