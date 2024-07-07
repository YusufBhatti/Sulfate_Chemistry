! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialises snow depth on the ground.

MODULE rcf_init_snowdep_mod
IMPLICIT NONE

! Description:
!     This subroutine initilaizes the snow depth on the ground
!     from the ground and canopy snow stores.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_SNOWDEP_MOD'

CONTAINS

SUBROUTINE rcf_init_snowdep( fields_in, field_count_in, hdr_in,    &
                             fields_out, field_count_out, hdr_out, &
                             data_source, snowdepth )

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
    stashcode_snow_grnd         , &
    stashcode_snowdep_grd_tile

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

USE jules_snow_mod, ONLY: rho_snow_const, cansnowtile

USE jules_vegetation_mod, ONLY: can_model

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
TYPE( field_type ), INTENT(INOUT) :: snowdepth

! Internal variables
TYPE( field_type ), POINTER  :: snow_tile
TYPE( field_type ), POINTER  :: snow_grnd

INTEGER      ::  errorstatus = 0
INTEGER      ::  pos         ! position in array
INTEGER      ::  i,n         ! loop counters

CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_INIT_SNOWDEP'
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
      ' Initialising snow depth on the ground.'
  CALL umPrint(umMessage,src='rcf_init_snowdep')
END IF

!----------------------------------------------------------------------
! Find the masses of snow on the ground or the canopy. The detailed
! checking of sizes may be superfluous, but it is retained for
! safety.
!----------------------------------------------------------------------

! Snow on tile
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
    errorstatus = 43
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
      errorstatus = 43
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
  errorstatus = 43
  WRITE (cmessage, '(A, I8, A, I8)')                &
        'Field sizes do not agree: level_size: ',   &
        snow_tile % level_size,                     &
        ' local_land_out: ',                            &
        local_land_out
  CALL ereport( routinename, errorstatus, cmessage )
END IF
IF ( snow_tile % levels /= ntiles ) THEN
  errorstatus = 43
  WRITE (cmessage, '(A, I3, A, I3)')                &
        'Field sizes do not agree: levels: ',       &
        snow_tile % levels,                         &
        ' ntiles: ',                                &
        ntiles
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!-----------------------------------------------------------------------------
! Determine the total amount of snow on the ground, allowing for
! different canopy models.
!-----------------------------------------------------------------------------

DO n=1,ntiles
  IF ( cansnowtile(n) .AND. (can_model == 4) ) THEN
    DO i=1,local_land_out
      snowdepth % DATA(i,n) = MAX(0.0, &
                  snow_grnd % DATA(i,n) / rho_snow_const)
    END DO
  ELSE
    DO i=1,local_land_out
      snowdepth % DATA(i,n) = MAX(0.0, &
                  snow_tile % DATA(i,n) / rho_snow_const)
    END DO
  END IF
END DO

!----------------------------------------------------------------------
! Clear up dynamic memory used.
!----------------------------------------------------------------------
CALL rcf_dealloc_field( snow_tile )
IF (can_model == 4) THEN
  CALL rcf_dealloc_field( snow_grnd )
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_init_snowdep
END MODULE rcf_init_snowdep_mod
