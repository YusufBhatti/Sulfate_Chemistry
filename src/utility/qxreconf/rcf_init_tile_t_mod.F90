! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Performs initiation of surface temperature on tiles/elevated tiles

MODULE Rcf_Init_Tile_T_Mod

!  Subroutine Rcf_Init_Tile_T   Initialise the surface temperature
!                               on tiles
!
! Description:
!   The surface temperature on tiles/elevated tiles needs initialising
!   in a number of situations, for example when the number of surface
!   tiles has changed, or the surface land use dataset has changed.
!
!   The simplest method is to initialise the surface temperature
!   on all tiles to equal the grid-box mean tstar.
!
!   However, if the number of tiles is unchanged, then this method
!   discards useful information held within the input field. Therefore
!   an option is provided to only initialise the temperature
!   on those tiles with zero fraction.
!
!   Further options may be added here.
!
! Method:
!   Initialise the surface temperature (on tiles with zero fraction)
!   to a sensible value. The current best-guess for this is the grid-box
!   mean surface temperature.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_TILE_T_MOD'

CONTAINS

SUBROUTINE Rcf_Init_Tile_T( fields, field_count, hdr,                   &
                            decomp_rcf, local_lsm,                      &
                            field_to_init_surft, l_init_tile_t_zerofrac )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY:    &
    stashcode_tstar,           &
    stashcode_tstar_land,      &
    stashcode_frac_surf_type,  &
    stashcode_prog_sec

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY:  &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE jules_sea_seaice_mod, ONLY: l_ctile

USE mask_compression, ONLY: compress_to_mask

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)     ! input fields
TYPE( field_type ), INTENT(INOUT)  :: field_to_init_surft    ! tstar on tiles
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp_rcf
LOGICAL, INTENT(IN)              :: local_lsm(field_to_init_surft % level_size)
LOGICAL, INTENT(IN)              :: l_init_tile_t_zerofrac

! Local variables
TYPE( field_type ), POINTER          :: tstar_gbm      ! grid-box mean tstar
TYPE( field_type ), POINTER          :: frac_tile      ! tile fraction
INTEGER                              :: size           ! for mask compression
INTEGER                              :: pos, pos_land  ! field positions
INTEGER                              :: i,k            ! looper

REAL, ALLOCATABLE :: tstar_gbm_landpts(:,:)  ! grid-box mean tstar 
                                             ! on land points only

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_TILE_T'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,'(A)') 'Initialising '// &
    TRIM(field_to_init_surft%stashmaster%name) // ' from T*'
  CALL umPrint(umMessage,src='rcf_init_tile_t_mod')
END IF

!----------------------------------------------------------------------
! Get grid-box mean tstar
!----------------------------------------------------------------------
! Find GBM Tstar
CALL Rcf_Locate( stashcode_prog_sec, stashcode_tstar,                   &
                 fields, field_count, pos)
IF (l_ctile) THEN
  ! If coastal tiling use land GBM Tstar if available
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_tstar_land,            &
                   fields, field_count, pos_land, zero_ok_arg = .TRUE.)
  IF (pos_land /= 0) pos = pos_land
END IF
tstar_gbm => fields(pos)
CALL Rcf_Alloc_Field( tstar_gbm )
CALL Rcf_Read_Field( tstar_gbm, hdr, decomp_rcf )

!----------------------------------------------------------------------
! Compress this to land points only
!----------------------------------------------------------------------

ALLOCATE( tstar_gbm_landpts( field_to_init_surft % level_size, 1) )
CALL compress_to_mask( tstar_gbm % Data,                                &
                       tstar_gbm_landpts,                               &
                       local_lsm,                                       &
                       tstar_gbm % level_size,                          &
                       size )

IF (l_init_tile_t_zerofrac) THEN
!----------------------------------------------------------------------
! Get tile fractions
!----------------------------------------------------------------------

  CALL Rcf_Locate ( stashcode_prog_sec,                                 &
                    stashcode_frac_surf_type,                           &
                    fields, field_count, pos)
  frac_tile => fields( pos )
  CALL Rcf_Alloc_Field( frac_tile )
  CALL Rcf_Read_Field( frac_tile, hdr, decomp_rcf )

!----------------------------------------------------------------------
! Over-write field_to_init_surft with mean tstar where fraction is 0
!----------------------------------------------------------------------

  DO k = 1, field_to_init_surft % levels
    DO i = 1, field_to_init_surft % level_size
      IF ( frac_tile % Data(i,k) == 0.0) THEN
        field_to_init_surft % Data(i,k) = tstar_gbm_landpts(i,1)
      END IF
    END DO
  END DO

  CALL Rcf_Dealloc_Field( frac_tile )

ELSE ! l_init_tile_t_zerofrac

  DO i = 1, field_to_init_surft % levels
    field_to_init_surft % Data(:,i) = tstar_gbm_landpts(:,1)
  END DO

END IF ! l_init_tile_t_zerofrac

!----------------------------------------------------------------------
! Deallocate memory
!----------------------------------------------------------------------

DEALLOCATE( tstar_gbm_landpts )
CALL Rcf_Dealloc_Field( tstar_gbm )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_Tile_T

END MODULE Rcf_Init_Tile_T_Mod
