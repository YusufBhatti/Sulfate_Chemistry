! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialisations for snow amount on tiles

MODULE Rcf_Init_Tile_Snow_Mod

! Subroutine Rcf_Init_Tile_Snow
!
! Description:
!   Initialises snow amount on tiles.
!
! Method:
!   Initialises snow amount on tiles from gridbox mean snow amount.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_TILE_SNOW_MOD'

CONTAINS

SUBROUTINE Rcf_Init_Tile_Snow( fields_out, field_count_out, hdr_out, &
                               snow_tile )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_mean_snow,       &
    stashcode_prog_sec

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field
USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_out


USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE mask_compression, ONLY: compress_to_mask

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)  :: snow_tile
TYPE( um_header_type ), INTENT(IN) :: hdr_out

INTEGER, INTENT(IN)                :: field_count_out

! Local variables
TYPE( field_type ), POINTER          :: mean_snow
INTEGER                              :: SIZE
INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper
REAL                                 :: mean_snow_land &
                                       ( snow_tile % level_size, 1)

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_TILE_SNOW'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Initialising snow amount on tiles from gb mean snow'
  CALL umPrint(umMessage,src='rcf_init_tile_snow_mod')
END IF

!----------------------------------------------------------------------
! Find snow amount in output fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_mean_snow,            &
                 fields_out, field_count_out, pos)
mean_snow => fields_out(pos)
CALL Rcf_Alloc_Field( mean_snow )
CALL Rcf_Read_Field( mean_snow, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Copy global mean-snow values into land only variable
!----------------------------------------------------------------------
CALL compress_to_mask( mean_snow % DATA,                 &
                     mean_snow_land,                    &
                     local_lsm_out,                     &
                     mean_snow % level_size,            &
                     SIZE )

! Copy snow amount into tile snow amount (for all psuedo levels) and
! write it out
!----------------------------------------------------------------------
DO i = 1, snow_tile % levels
  snow_tile % DATA(:,i) = mean_snow_land(:,1)
END DO

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( mean_snow )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_Tile_Snow

END MODULE Rcf_Init_Tile_Snow_Mod
