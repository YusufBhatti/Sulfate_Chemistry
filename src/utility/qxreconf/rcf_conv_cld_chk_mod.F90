! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Performs a sanity check on convective cloud variables

MODULE Rcf_Conv_Cld_Chk_Mod

IMPLICIT NONE

!  Subroutine Rcf_Conv_Cld_Chk - sanity checks.
!
! Description:
!   Ensures that convective cloud top is at least as high as
!   convective cloud base!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CONV_CLD_CHK_MOD'

CONTAINS

SUBROUTINE Rcf_Conv_Cld_Chk( fields, field_count, grid, decomp, hdr )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_done

USE um_stashcode_mod, ONLY: &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_prog_sec

USE UM_ParCore, ONLY: &
    nproc

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)
TYPE( grid_type ), INTENT(IN)      :: grid
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp

! Local variables
INTEGER                            :: pos_cct
INTEGER                            :: pos_ccb
INTEGER                            :: i
INTEGER                            :: base_changed
INTEGER                            :: top_changed
INTEGER                            :: istat
TYPE( field_type ), POINTER        :: cct
TYPE( field_type ), POINTER        :: ccb

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CONV_CLD_CHK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Find and read conv_cloud_top and conv_cloud_base
!-------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields, field_count, pos_cct, zero_ok_arg = .TRUE. )
CALL Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields, field_count, pos_ccb, zero_ok_arg = .TRUE. )

IF (pos_ccb /= 0 .AND. pos_cct /= 0 ) THEN
  cct => fields( pos_cct )
  CALL Rcf_Alloc_Field( cct )
  CALL Rcf_Read_Field( cct, hdr, decomp )

  ccb => fields( pos_ccb )
  CALL Rcf_Alloc_Field( ccb )
  CALL Rcf_Read_Field( ccb, hdr, decomp )
  !--------------------------------------------------------------------
  ! Compare all points
  !--------------------------------------------------------------------
  IF (ccb % interp == interp_done .OR. cct % interp == interp_done ) THEN
    base_changed = 0
    top_changed  = 0

    DO i = 1, cct % level_size
      IF ( ccb % Data_Int(i,1) == 0 .AND. cct % Data_Int(i,1) > 0 ) THEN
        ccb % Data_Int(i,1) = 1
        base_changed = 1
      END IF

      IF ( cct % Data_Int(i,1) == ccb % Data_Int(i,1) .AND.  &
           cct % Data_Int(i,1) /= 0) THEN

        IF ( cct % Data_Int(i,1) /= grid % model_levels ) THEN
          ! Convective cloud top and base the same
          ! increment top by one
          cct % Data_Int(i,1) = cct % Data_Int(i,1) + 1
          top_changed = 1

        ELSE
          ! Convective Cloud top and base the same
          ! decrement base by 1
          ccb % Data_Int(i,1) = ccb % Data_Int(i,1) - 1
          base_changed = 1
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------
    ! Synchronise `changed' flags
    !---------------------------------------------------------------------
    CALL GC_Imax( 1, nproc, istat, base_changed )
    CALL GC_Imax( 1, nproc, istat, top_changed  )

    !---------------------------------------------------------------------
    ! Write out changed fields
    !---------------------------------------------------------------------
    IF (top_changed == 1) THEN
      CALL Rcf_Write_Field( cct, hdr, decomp )
    END IF

    IF (base_changed == 1) THEN
      CALL Rcf_Write_Field( ccb, hdr, decomp )
    END IF

  END IF

  CALL Rcf_Dealloc_Field( cct )
  CALL Rcf_Dealloc_Field( ccb )

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Conv_Cld_Chk
END MODULE Rcf_Conv_Cld_Chk_Mod
