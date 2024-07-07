! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up horizontal interpolation weights

MODULE Rcf_init_h_interp_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!  Subroutine Rcf_init_h_interp - initialises horizontal interp. weights
!
! Description:
!   Sets up interpolation weights based on choice of scheme.
!
! Method:
!   Allocates space for and sets up weights
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_H_INTERP_MOD'

CONTAINS

SUBROUTINE Rcf_Init_H_Interp( grid_in, grid_out, hdr_in, hdr_out )

USE Rcf_Interp_Weights_Mod  ! Almost all of this used.

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_UMhead_Mod, ONLY:    &
    um_header_type

USE Ereport_mod, ONLY: &
    Ereport

USE Rcf_Lsm_Mod, ONLY: &
    cyclic

USE rcf_h_int_init_bl_mod, ONLY: &
    rcf_h_int_init_bl
USE rcf_h_int_init_cart_bl_mod, ONLY: &
    rcf_h_int_init_cart_bl

USE Rcf_HeadAddress_Mod, ONLY:&
    FH_GridStagger_Endgame, FH_GridStagger,  &
    FH_CoordSystem,                          &
    FH_CoordSystem_Cartesian

USE nlstcall_mod, ONLY: &
    LTimer

USE missing_data_mod, ONLY: imdi

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (grid_type), INTENT(IN)      :: grid_in
TYPE (grid_type), INTENT(IN)      :: grid_out
TYPE (um_header_type), INTENT(IN) :: hdr_in
TYPE (um_header_type), INTENT(IN) :: hdr_out

! Local data
INTEGER                       :: i,j   ! loop counters
INTEGER                       :: interp_arr_sz
INTEGER                       :: interp_rows_sz
INTEGER                       :: max_rows_out, max_field_out
INTEGER                       :: max_rows_in
INTEGER                       :: ErrorStatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RCF_INIT_H_INTERP'
CHARACTER (LEN=errormessagelength)    :: Cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (LTimer) CALL Timer( RoutineName, 3)

! Set the cyclic flag for coastal adjustment
IF ( Hdr_Out % Fixhd(4) /= 3 .AND. Hdr_Out % FixHd(4) /= 103 ) THEN
  cyclic=.TRUE.
ELSE
  cyclic=.FALSE.
END IF

!----------------------------------------------------------------
! Test if we know which interpolation scheme we are doing!
!----------------------------------------------------------------
IF ( h_int_method == imdi ) THEN
  Cmessage = 'No interpolation method specified - cannot set weights'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! To deal with Endgame need to find maximum grid size.
! IN EG v has more rows than p and thus glob_v_field can be larger.

interp_arr_sz = MAX(grid_out % glob_v_field,                      &
                    grid_out % glob_p_field)
interp_rows_sz = MAX(grid_out % glob_v_rows + 1,                  &
                     grid_out % glob_p_rows + 1)

SELECT CASE ( h_int_method )
CASE ( bilinear, nearest_neighbour )

  ALLOCATE( bl_index_b_l( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_b_r( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_t_l( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_t_r( interp_arr_sz, IDIM ) )
  ALLOCATE( coeff1( interp_arr_sz ) )
  ALLOCATE( coeff2( interp_arr_sz ) )
  ALLOCATE( coeff3( grid_in % glob_p_field ) )
  ALLOCATE( coeff4( grid_in % glob_p_field ) )
  ALLOCATE( weight_t_r( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_b_r( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_t_l( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_b_l( interp_arr_sz, IDIM ) )


  IF (hdr_out% FixHd( FH_CoordSystem ) == FH_CoordSystem_Cartesian .AND.  &
      hdr_in% FixHd( FH_CoordSystem ) == FH_CoordSystem_Cartesian) THEN

    CALL rcf_h_int_init_cart_bl( grid_in, grid_out, hdr_in, hdr_out )

  ELSE

    CALL rcf_h_int_init_bl( grid_in, grid_out, hdr_in, hdr_out )

  END IF

CASE ( area_weighted )

  ALLOCATE( aw_index_targ_lhs( grid_out % glob_p_row_length + 1, IDIM ))
  ALLOCATE( aw_index_targ_top( interp_rows_sz, IDIM ))
  ALLOCATE( aw_colat_t( interp_rows_sz, IDIM ))
  ALLOCATE( aw_long_l( grid_Out % glob_p_row_length + 1, IDIM ))
  ALLOCATE( bl_index_b_l( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_b_r( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_t_l( interp_arr_sz, IDIM ) )
  ALLOCATE( bl_index_t_r( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_t_r( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_b_r( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_t_l( interp_arr_sz, IDIM ) )
  ALLOCATE( weight_b_l( interp_arr_sz, IDIM ) )

  IF (hdr_out % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    max_rows_out = grid_out % glob_v_rows
    max_field_out = grid_out % glob_v_field
  ELSE
    max_rows_out = grid_out % glob_p_rows
    max_field_out = grid_out % glob_p_field
  END IF

  IF (hdr_in % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    max_rows_in = grid_in % glob_v_rows
  ELSE
    max_rows_in = grid_in % glob_p_rows
  END IF

  ! DEPENDS ON: h_int_init_aw
  CALL h_int_init_aw( icof, IDIM, max_field_out,                          &
               max_rows_in, max_rows_out,                                 &
               grid_in % glob_p_row_length, grid_out % glob_p_row_length, &
               grid_out % global, hdr_in % FixHd,                         &
               hdr_out % FixHd, hdr_in % RealC, hdr_out % RealC,          &
               aw_area_box, aw_index_targ_lhs, aw_index_targ_top,         &
               bl_index_b_l, bl_index_b_r,                                &
               bl_index_t_l, bl_index_t_r,                                &
               aw_colat_t, aw_long_l, weight_t_r, weight_b_r,             &
               weight_t_l, weight_b_l )

CASE DEFAULT
  Cmessage = 'Interpolation method not supported'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_H_Interp

END MODULE Rcf_init_h_interp_mod
