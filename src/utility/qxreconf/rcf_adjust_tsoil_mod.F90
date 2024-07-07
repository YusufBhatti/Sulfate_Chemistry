! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Performs orographic adjustment of soil temperatures

MODULE Rcf_Adjust_Tsoil_Mod

!  Subroutine Rcf_Adjust_Tsoil    Adjusts Tsoil when an ancillary
!                                 orography has been used with an
!                                 interpolated Tsoil
!
! Description:
!   A similar method to that used in rcf_adjust_pstar
!
! Method:
!   Assuming a constant lapse rate, each soil level temperature (Ts)
!   is calculated based on a temperature at a reference height (Tr)
!              Ts = Tr - lapse * (Zr - Z0) 
!   
!   This code was developed for use in the COLPEX project, and assumes
!   the 'magic number' lapse rate of dT/dZ = -6K/km, based on analysis
!   of soil temperature versus height close to the COLPEX domain
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ADJUST_TSOIL_MOD'

CONTAINS

SUBROUTINE Rcf_Adjust_Tsoil(fields_out, field_count_out,                &
                            decomp_rcf_output, hdr_out,                 &
                            orog_in, orog_out,                          &
                            input_grid, output_grid)

USE Rcf_Alloc_Field_Mod, ONLY:                                          &
    Rcf_Alloc_Field,                                                    &
    Rcf_Dealloc_Field

USE Rcf_Field_Equals_Mod, ONLY:                                         &
    Rcf_Field_Equals

USE Rcf_Field_Type_Mod, ONLY:                                           &
    field_type

USE Rcf_Grid_Type_Mod, ONLY:                                            &
    Grid_Type

USE Rcf_Interpolate_Mod, ONLY:                                          &
    Rcf_Interpolate

USE Rcf_Interp_Weights_Mod, ONLY:                                       &
    h_int_active

USE Rcf_Locate_Mod, ONLY:                                               &
    Rcf_Locate

USE Rcf_Lsm_Mod, ONLY:                                                  &
    local_lsm_out

USE Rcf_Read_Field_Mod, ONLY:                                           &
    Rcf_Read_Field

USE Rcf_Set_Interp_Flags_Mod, ONLY:                                     &
    interp_h_only,                                                      &
    interp_copy,                                                        &
    interp_no_op

USE Rcf_UMhead_Mod, ONLY:                                               &
    um_header_type

USE Rcf_Write_Field_Mod, ONLY:                                          &
    Rcf_Write_Field

USE mask_compression, ONLY: compress_to_mask

USE um_stashcode_mod, ONLY:                                             &
    stashcode_soil_temp,                                                &
    stashcode_prog_sec

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( Grid_Type ), INTENT(IN)       :: input_grid
TYPE( Grid_Type ), INTENT(IN)       :: output_grid
TYPE( Field_Type ), INTENT(INOUT)   :: orog_in
TYPE( Field_Type ), INTENT(IN)      :: orog_out
TYPE( Field_Type ), POINTER         :: fields_out(:)
TYPE( um_header_type ), INTENT(IN)  :: hdr_out
INTEGER, INTENT(IN)                 :: field_count_out
INTEGER, INTENT(IN)                 :: decomp_rcf_output

! Local Variables/Parameters
TYPE( Field_Type )    :: orog_interp         ! interpolated input
                                             ! orography
TYPE( Field_Type )    :: dummy               ! dummy field
REAL, PARAMETER       :: soil_lapse = -0.006 ! soil lapse rate assume
REAL                  :: Z0o                 ! point rcf output height
REAL                  :: Z0i                 ! point interpolated height
INTEGER               :: i                   ! Looper
INTEGER               :: k                   ! Looper 2
INTEGER               :: pos                 ! position in array
INTEGER               :: size                ! for mask compression

TYPE( Field_Type ), POINTER         :: soil_temp

REAL, ALLOCATABLE :: orog_out_landpts(:,:)
REAL, ALLOCATABLE :: orog_interp_landpts(:,:)

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_ADJUST_TSOIL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------

! Soil temperature
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,               &
                 fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
CALL Rcf_Alloc_Field( soil_temp )
CALL Rcf_Read_Field( soil_temp, hdr_out, decomp_rcf_output)

!--------------------------------------------------------------------
! Set up the orography fields to obtain interpolated orography
!--------------------------------------------------------------------

CALL Rcf_Field_Equals( orog_interp, orog_out)

IF (h_int_active) THEN
  orog_in % interp   = interp_h_only
ELSE
  orog_in % interp   = interp_copy
END IF

orog_interp % interp = interp_no_op

CALL Rcf_Alloc_Field( orog_interp )
CALL Rcf_Interpolate( orog_in, orog_interp, input_grid, output_grid,    &
                      dummy, dummy )

!--------------------------------------------------------------------
! Compress the orography to land points
!--------------------------------------------------------------------
ALLOCATE( orog_out_landpts( soil_temp % level_size, 1) )
ALLOCATE( orog_interp_landpts( soil_temp % level_size, 1) )

CALL compress_to_mask (orog_out % Data,                                 &
                       orog_out_landpts,                                &
                       local_lsm_out,                                   &
                       orog_out % level_size,                           &
                       size )
CALL compress_to_mask (orog_interp % Data,                              &
                       orog_interp_landpts,                             &
                       local_lsm_out,                                   &
                       orog_interp % level_size,                        &
                       size )

!--------------------------------------------------------------------
! Can now do the tsoil adjustment calculation as specified above
!--------------------------------------------------------------------
DO i = 1, soil_temp % level_size

  Z0o = orog_out_landpts( i, 1 )
  Z0i = orog_interp_landpts( i, 1 )

  DO k = 1, output_grid % st_levels

    soil_temp % Data( i, k ) = soil_temp % Data( i, k ) -               &
                        (Z0i - Z0o)*soil_lapse

  END DO
END DO

!----------------------------------------------------------------------
! Write out corrected soil temperature
!----------------------------------------------------------------------
CALL Rcf_Write_Field( soil_temp, hdr_out, decomp_rcf_output)

!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
CALL Rcf_Dealloc_Field( orog_interp )
CALL Rcf_Dealloc_Field( soil_temp )

DEALLOCATE( orog_out_landpts )
DEALLOCATE( orog_interp_landpts )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Adjust_Tsoil

END MODULE Rcf_Adjust_Tsoil_Mod
