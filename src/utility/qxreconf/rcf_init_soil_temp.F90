! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! initialises soil temperature to surface temperature value

MODULE Rcf_init_soil_temp_Mod

!  Subroutine Rcf_init_soil_temp
!
! Description:
!   initialise soil temperature to surface temperature value
!
! Method:
!   Copy grid-box mean surface temperature to all soil levels
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_SOIL_TEMP_MOD'

CONTAINS

SUBROUTINE Rcf_init_soil_temp(fields_out, field_count_out,              &
                              decomp_rcf_output, hdr_out,               &
                              current_field_out)

USE Rcf_Alloc_Field_Mod, ONLY: Rcf_Alloc_Field, Rcf_Dealloc_Field

USE Rcf_Field_Type_Mod, ONLY: field_type

USE Rcf_Locate_Mod, ONLY: Rcf_Locate

USE Rcf_Lsm_Mod, ONLY: local_lsm_out

USE Rcf_Read_Field_Mod, ONLY: Rcf_Read_Field

USE Rcf_UMhead_Mod, ONLY: um_header_type

USE mask_compression, ONLY: compress_to_mask

USE um_stashcode_mod, ONLY: stashcode_prog_sec, stashcode_tstar

USE yomhook,   ONLY: lhook, dr_hook

USE parkind1,  ONLY: jprb, jpim

IMPLICIT NONE

! Arguments
TYPE( Field_Type ), POINTER         :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)   :: current_field_out
TYPE( um_header_type ), INTENT(IN)  :: hdr_out
INTEGER, INTENT(IN)                 :: field_count_out
INTEGER, INTENT(IN)                 :: decomp_rcf_output

! Local Variables/Parameters
INTEGER               :: i                   ! Looper
INTEGER               :: k                   ! Looper 2
INTEGER               :: pos                 ! position in array
INTEGER               :: field_size          ! for mask compression

REAL, ALLOCATABLE :: tstar_gbm_landpts(:)  ! grid-box mean tstar 
                                           ! on land points only

TYPE( field_type ), POINTER :: tstar_gbm      ! grid-box mean tstar

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_SOIL_TEMP'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Get grid-box mean tstar
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_tstar,                   &
     fields_out, field_count_out, pos)
tstar_gbm => fields_out(pos)
CALL Rcf_Alloc_Field( tstar_gbm )
CALL Rcf_Read_Field( tstar_gbm, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Compress this to land points only
!----------------------------------------------------------------------
ALLOCATE( tstar_gbm_landpts( current_field_out % level_size ) )
CALL compress_to_mask( tstar_gbm % Data,                                &
     tstar_gbm_landpts,                                                 &
     local_lsm_out,                                                     &
     tstar_gbm % level_size,                                            &
     field_size )

!--------------------------------------------------------------------
! Set the soil temperature to match the surface
!--------------------------------------------------------------------
DO k = 1, current_field_out % levels
  DO i = 1, current_field_out % level_size

    current_field_out % Data( i, k ) = tstar_gbm_landpts(i)

  END DO
END DO

!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
DEALLOCATE( tstar_gbm_landpts )
CALL Rcf_Dealloc_Field( tstar_gbm )
      
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_init_soil_temp

END MODULE Rcf_init_soil_temp_Mod
