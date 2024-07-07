! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Choose weights for horizontal interpolation.

MODULE Rcf_select_weights_mod

!  Subroutine Rcf_Select_Weights - choose weights for horiz. interp.
!
! Description:
! This module calculates which part of the weights array is
! required by a given horizontal interpolation method, based
! on the type of field being used.
!
! Method:
! The correct weights are chosen on grid-code and interpolation
! method chosen. Pointers are set to relevant weight array slices.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SELECT_WEIGHTS_MOD'

CONTAINS

SUBROUTINE Rcf_select_weights( ptr_bl_index_b_l, ptr_bl_index_b_r, &
                               ptr_bl_index_t_l, ptr_bl_index_t_r, &
                               ptr_weight_b_l, ptr_weight_b_r,     &
                               ptr_weight_t_l, ptr_weight_t_r,     &
                               ptr_aw_index_targ_lhs,              &
                               ptr_aw_index_targ_top,              &
                               ptr_aw_colat_t, ptr_aw_long_l,      &
                               gridcode, section, item )

USE Rcf_Interp_Weights_Mod, ONLY: &
    bilinear,                  &
    area_weighted,             &
    nearest_neighbour,         &
    h_int_method,              &
    weight_b_l, weight_b_r,    &
    weight_t_l, weight_t_r,    &
    bl_index_b_l, bl_index_b_r,&
    bl_index_t_l, bl_index_t_r,&
    aw_index_targ_lhs,         &
    aw_index_targ_top,         &
    aw_colat_t, aw_long_l

USE um_stashcode_mod, ONLY:                        &
    stashcode_prog_sec,      stashcode_tracer_sec,    &
    stashcode_ukca_sec,                               &
    stashcode_u,             stashcode_surf_z_curr,   &
    stashcode_w,             stashcode_3d_cca,        &
    stashcode_v,             stashcode_surf_m_curr,   &
    stashcode_riv_sequence,  stashcode_3d_ccw,        &
    stashcode_u_compnt_pert, stashcode_v_compnt_pert, &
    stashcode_3d_nat_so2_em, stashcode_3d_oh_conc

USE model_domain_mod, ONLY: &
    l_cartesian

USE Ereport_mod, ONLY: &
    Ereport

USE cppxref_mod, ONLY:                               &
    ppx_atm_tall,            ppx_atm_uall,            &
    ppx_atm_cuall,           ppx_atm_cvall,           &
    ppx_atm_tsea,            ppx_atm_usea,            &
    ppx_atm_ozone,           ppx_atm_compressed

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)           :: gridcode
INTEGER, INTENT(IN)           :: section
INTEGER, INTENT(IN)           :: item
INTEGER, POINTER              :: ptr_bl_index_b_l(:)
INTEGER, POINTER              :: ptr_bl_index_b_r(:)
INTEGER, POINTER              :: ptr_bl_index_t_l(:)
INTEGER, POINTER              :: ptr_bl_index_t_r(:)
INTEGER, POINTER              :: ptr_aw_index_targ_lhs(:)
INTEGER, POINTER              :: ptr_aw_index_targ_top(:)
REAL, POINTER                 :: ptr_aw_colat_t(:)
REAL, POINTER                 :: ptr_aw_long_l(:)
REAL, POINTER                 :: ptr_weight_b_l(:)
REAL, POINTER                 :: ptr_weight_b_r(:)
REAL, POINTER                 :: ptr_weight_t_l(:)
REAL, POINTER                 :: ptr_weight_t_r(:)

! Local Data
INTEGER                      :: INDEX
INTEGER                      :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SELECT_WEIGHTS'
CHARACTER (LEN=errormessagelength)  :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Choose the index based on the grid code
!------------------------------------------------------------------
SELECT CASE ( gridcode )
CASE ( ppx_atm_tall, ppx_atm_tsea, ppx_atm_compressed )
  INDEX = 1

CASE ( ppx_atm_uall, ppx_atm_cuall, ppx_atm_usea )

  SELECT CASE ( section )

  CASE ( stashcode_prog_sec )
    SELECT CASE ( item )
      ! Removed items 135, 139, 148 at vn6.6 since no STASHmaster record exists
    CASE (stashcode_u,             &
          stashcode_surf_z_curr,   &
          stashcode_3d_nat_so2_em, &
          stashcode_w,             &
          stashcode_u_compnt_pert, &
          stashcode_3d_cca )
      INDEX = 2

      ! Removed items 136, 140 and 149 at vn6.6 since no STASHmaster record exists
    CASE (stashcode_v,             &
          stashcode_surf_m_curr,   &
          stashcode_3d_oh_conc,    &
          stashcode_riv_sequence,  &
          stashcode_v_compnt_pert, &
          stashcode_3d_ccw )
      INDEX = 3

    CASE DEFAULT
      INDEX = 2
    END SELECT
  CASE DEFAULT
    ! Set other sections to default.
    INDEX = 2
  END SELECT

CASE ( ppx_atm_cvall )
  INDEX = 3

CASE ( ppx_atm_ozone )
  ! This assumes zonal ozone - not sure how to distinguish otherwise
  INDEX = 4

  ! Either output grid is Cartesian or both input and output are Cartesian
  IF (l_cartesian) THEN
    Cmessage = 'Interpolation from a zonal input field NOT TESTED' 
    ErrorStatus = 30
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
CASE DEFAULT

END SELECT

!------------------------------------------------------------------
! Set pointers for relevant weight slices
!------------------------------------------------------------------


SELECT CASE ( h_int_method )
CASE ( bilinear, nearest_neighbour )
  ptr_bl_index_b_l => bl_index_b_l( :, INDEX )
  ptr_bl_index_b_r => bl_index_b_r( :, INDEX )
  ptr_bl_index_t_l => bl_index_t_l( :, INDEX )
  ptr_bl_index_t_r => bl_index_t_r( :, INDEX )
  ptr_weight_b_l   => weight_b_l  ( :, INDEX )
  ptr_weight_b_r   => weight_b_r  ( :, INDEX )
  ptr_weight_t_l   => weight_t_l  ( :, INDEX )
  ptr_weight_t_r   => weight_t_r  ( :, INDEX )
  NULLIFY( ptr_aw_index_targ_lhs )
  NULLIFY( ptr_aw_index_targ_top )
  NULLIFY( ptr_aw_colat_t )
  NULLIFY( ptr_aw_long_l )


CASE ( area_weighted )
  ptr_aw_index_targ_lhs => aw_index_targ_lhs( :, INDEX )
  ptr_aw_index_targ_top => aw_index_targ_top( :, INDEX )
  ptr_aw_colat_t        => aw_colat_t       ( :, INDEX )
  ptr_aw_long_l         => aw_long_l        ( :, INDEX )
  NULLIFY( ptr_bl_index_b_l )
  NULLIFY( ptr_bl_index_b_r )
  NULLIFY( ptr_bl_index_t_l )
  NULLIFY( ptr_bl_index_t_r )
  NULLIFY( ptr_weight_b_l)
  NULLIFY( ptr_weight_b_r)
  NULLIFY( ptr_weight_t_l)
  NULLIFY( ptr_weight_t_r)
CASE DEFAULT
  Cmessage = 'Interpolation method not yet supported'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_select_weights

END MODULE Rcf_select_weights_mod
