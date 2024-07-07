! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE field_length_mod

USE atm_fields_bounds_mod
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


!  interface
!    function field_length(grid_type,halo_type,levels)
!      integer :: field_length ! function output
!      integer, intent(in) :: grid_type
!      integer, intent(in) :: halo_type
!      integer, intent(in) :: levels
!    end function field_length
!  end interface
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

      ! the set of grid types that can be specified in a stashmaster file
INTEGER, PARAMETER :: theta_points = 1 !1
INTEGER, PARAMETER :: theta_points_land_only = 2 !2
INTEGER, PARAMETER :: theta_points_sea_only = 3 !3
INTEGER, PARAMETER :: zonal_theta_points = 4 !4
INTEGER, PARAMETER :: merid_theta_points = 5 !5
INTEGER, PARAMETER :: uv_points = 6 !11
INTEGER, PARAMETER :: uv_points_land_only = 7 !12
INTEGER, PARAMETER :: uv_points_sea_only = 8 !13
INTEGER, PARAMETER :: zonal_uv_points = 9 !14
INTEGER, PARAMETER :: merid_uv_points = 10 !15
INTEGER, PARAMETER :: scalar_points = 11 !17
INTEGER, PARAMETER :: u_points = 12 !18
INTEGER, PARAMETER :: v_points = 13 !19
INTEGER, PARAMETER :: land_points = 14 !21
INTEGER, PARAMETER :: ozone_points = 15 !22
INTEGER, PARAMETER :: river_points = 16 !23
INTEGER, PARAMETER :: lbc_points = 17 !25
INTEGER, PARAMETER :: lbc_theta_points = 18 !26
INTEGER, PARAMETER :: lbc_u_points = 19 !27
INTEGER, PARAMETER :: lbc_v_points = 20 !28
INTEGER, PARAMETER :: lbc_orography_points = 21 !29
INTEGER, PARAMETER :: non_standard_ocn_points = 22 !30
INTEGER, PARAMETER :: comp_ocn_mass_points = 23 !31
INTEGER, PARAMETER :: comp_ocn_velocity_points = 24 !32
INTEGER, PARAMETER :: cyclic_ocn_mass_points = 25 !36
INTEGER, PARAMETER :: cyclic_ocn_velocity_points = 26 !37
INTEGER, PARAMETER :: cyclic_ocn_u_points = 27 !38
INTEGER, PARAMETER :: cyclic_ocn_v_points = 28 !39
INTEGER, PARAMETER :: ocn_mass_points = 29 !41
INTEGER, PARAMETER :: ocn_velocity_points = 30 !42
INTEGER, PARAMETER :: zonal_ocn_mass_points = 31 !43
INTEGER, PARAMETER :: zonal_ocn_velocity_points = 32 !44
INTEGER, PARAMETER :: merid_ocn_mass_points = 33 !45
INTEGER, PARAMETER :: merid_ocn_velocity_points = 34 !46
INTEGER, PARAMETER :: scalar_ocn_points = 35 !47
INTEGER, PARAMETER :: lbc_ocn_points = 36 !51
INTEGER, PARAMETER :: lbc_ocn_t_points = 37 !52
INTEGER, PARAMETER :: lbc_ocn_uv_points = 38 !53
INTEGER, PARAMETER :: wave_points = 39 !60
INTEGER, PARAMETER :: wave_points_sea_only = 40 !62
INTEGER, PARAMETER :: lbc_wave_points = 41 !65

! the set of halo types that can be specified in a stashmater file
INTEGER, PARAMETER :: single_halo = 1
INTEGER, PARAMETER :: extended_halo = 2
INTEGER, PARAMETER :: no_halo = 3

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FIELD_LENGTH_MOD'

CONTAINS

INTEGER FUNCTION field_length(grid_type,halo_type,levels)

USE UM_ParVars
USE nlsizes_namelist_mod, ONLY:                                        &
    land_field, model_levels, n_cca_lev, river_row_length, river_rows, &
    sm_levels, st_levels, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, &
    tr_ukca, tr_vars

IMPLICIT NONE

! Function arguments:
INTEGER, INTENT(IN) :: grid_type
INTEGER, INTENT(IN) :: halo_type
INTEGER, INTENT(IN) :: levels

! Local variables:
INTEGER :: length
INTEGER :: x, y, halo_x, halo_y

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FIELD_LENGTH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
SELECT CASE (grid_type)
CASE (theta_points                                              &
      ,theta_points_land_only ,theta_points_sea_only            &
      ,uv_points ,uv_points_land_only ,uv_points_sea_only       &
      ,ozone_points                                             &
     )
  x = tdims%i_end
  y = tdims%j_end
CASE (u_points)
  x = udims%i_len
  y = udims%j_end
CASE (v_points)
  x = vdims%i_end
  y = vdims%j_len
CASE (land_points)
  ! this is a special case (assume there are no halos)
  IF (halo_type  /=  no_halo) THEN
    WRITE(umMessage,*) "WARNING MSG FROM FIELD_LENGTH()"
    CALL umPrint(umMessage,src='field_length_mod')
    WRITE(umMessage,*) "POSSIBLE ERROR DETERMINING FIELD LENGTH"
    CALL umPrint(umMessage,src='field_length_mod')
    WRITE(umMessage,*) "UNHANDLED HALO_TYPE (FOR LAND_FIELDS):",halo_type
    CALL umPrint(umMessage,src='field_length_mod')
    WRITE(umMessage,*) "ASSUMING 'no_halo'"
    CALL umPrint(umMessage,src='field_length_mod')
    WRITE(umMessage,*) "END WARNING MSG"
    CALL umPrint(umMessage,src='field_length_mod')
  END IF
  length = land_field * levels
  field_length = length
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
CASE (river_points)
  x = river_row_length
  y = river_rows
  !! OCEAN and WAVE grid types are not handled by this routine.
CASE DEFAULT !! unspecified/unhandled grid_type
  WRITE(umMessage,*) "WARNING MSG FROM FIELD_LENGTH()"
  CALL umPrint(umMessage,src='field_length_mod')
  WRITE(umMessage,*) "POSSIBLE ERROR DETERMINING FIELD LENGTH"
  CALL umPrint(umMessage,src='field_length_mod')
  WRITE(umMessage,*) "UNHANDLED GRID_TYPE: ",grid_type
  CALL umPrint(umMessage,src='field_length_mod')
  WRITE(umMessage,*) "ASSUMING (ROW_LENGTH*ROWS)"
  CALL umPrint(umMessage,src='field_length_mod')
  WRITE(umMessage,*) "END WARNING MSG"
  CALL umPrint(umMessage,src='field_length_mod')
  x = tdims%i_end
  y = tdims%j_end
END SELECT

SELECT CASE (halo_type)
CASE (extended_halo)
  halo_x = halo_i
  halo_y = halo_j
CASE (single_halo)
  halo_x = offx
  halo_y = offy
CASE DEFAULT ! (no_halo)
  halo_x = 0
  halo_y = 0
END SELECT

length =  (((2 * halo_x) + x) * ((2 * halo_y) + y)) * levels
field_length = length
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION field_length

!- End of Function code

END MODULE field_length_mod
