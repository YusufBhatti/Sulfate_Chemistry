! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialises bilinear interpolation weights

MODULE rcf_h_int_init_cart_bl_mod

!  Subroutine rcf_h_int_init_cart_bl
!
! Description:
!   Sets up the interpolation and weights for bilinear
!   horizontal interpolation for a cartesian to cartesian grid
!   Cartesian grids can be:
!   1. bicylic  (wrap around X and Y) most common setup
!   2. a channel in X direction not Y (wrap around in X)
!   3. no wrap around will require boundary conditions.
!
! Method:
!   Weights stored in the interp_weights_mod module
!   Seperate calculations for P, U, V and P zonal points.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_H_INT_INIT_CART_BL_MOD'

CONTAINS

SUBROUTINE rcf_h_int_init_cart_bl( grid_in, grid_out, hdr_in, hdr_out )


USE Rcf_Interp_Weights_Mod, ONLY:                                &
    bl_index_b_l, bl_index_b_r, bl_index_t_l, bl_index_t_r,      &
    weight_b_l, weight_b_r, weight_t_l, weight_t_r

USE Ereport_Mod, ONLY: &
    Ereport

USE lam_config_inputs_mod, ONLY: &
    target_delta_lat,            &
    target_delta_lon,            &
    target_frstlata,             &
    target_frstlona

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Headaddress_Mod, ONLY:        &
    FH_RowDepCStart, FH_HorizGrid,    &
    RC_LongSpacing,  RC_LatSpacing,   &
    RC_FirstLat,     RC_FirstLong,    &
    RC_PoleLat,      RC_PoleLong,     &
    FH_GridStagger,  FH_GridStagger_C,&
    FH_GridStagger_Endgame,           &
    FH_GridStagger_A,                 &
    FH_CoordSystem,                   &
    FH_CoordSystem_Cartesian


USE umPrintMgr, ONLY:      &
    umPrint,  &
    umMessage,&
    PrStatus_Oper

USE f_shum_horizontal_field_interp_mod, ONLY :                                 &
                    f_shum_cart_horizontal_field_bi_lin_interp_get_coeffs

USE missing_data_mod, ONLY: rmdi
USE rcf_nlist_recon_technical_mod, ONLY: &
    l_basic_interp

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim


IMPLICIT NONE

! Arguments
TYPE( grid_type ),      INTENT(IN)   :: grid_in
TYPE( grid_type ),      INTENT(IN)   :: grid_out
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
TYPE( um_header_type ), INTENT(IN)   :: hdr_out

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
INTEGER            :: i                ! loop counters
INTEGER            :: j                ! 
INTEGER            :: ij               ! 
INTEGER            :: ErrorStatus
INTEGER            :: glob_rows_in     ! Number of rows in source grid
INTEGER            :: glob_field_in    ! Number of points in source grid
INTEGER            :: glob_rows_out    ! Number of rows in target grid
INTEGER            :: glob_field_out   ! Number of points in target grid

INTEGER            :: source_grid_type ! source grid type
INTEGER            :: icode            ! error code

! Allowed source grid types - should be in a module
INTEGER, PARAMETER :: default_grid_type = 0  ! no wrap around
INTEGER, PARAMETER :: bicylic_grid_type = 1  ! wrap X and Y
INTEGER, PARAMETER :: channel_grid_type = 2  ! wrap X only

REAL               :: grid_tol      ! Grid tolerance between grids.
REAL               :: max_y_in      ! Max size of source grid in Y (m)    
REAL               :: max_x_in      ! Max size of source grid in X (m)    
REAL               :: dx            ! grid spacing (m)
REAL               :: dy            ! grid spacing (m)
REAL               :: s_x_in        ! Source grid bottom left X (m)
REAL               :: s_y_in        ! Source grid bottom left Y (m)
REAL               :: e_x_in        ! Source grid top right X (m)
REAL               :: e_y_in        ! Source grid top right Y (m)
REAL               :: s_x_out       ! Target grid bottom left X (m)
REAL               :: s_y_out       ! Target grid bottom left Y (m)
REAL               :: e_x_out       ! Target grid top right X (m)
REAL               :: e_y_out       ! Target grid top right Y (m)

REAL, ALLOCATABLE  :: x_in( : )     ! X coords of source grid (m)
REAL, ALLOCATABLE  :: y_in( : )     ! Y coords of source grid (m)
REAL, ALLOCATABLE  :: x_out( : )    ! X coords of target grid (m)
REAL, ALLOCATABLE  :: y_out( : )    ! Y coords of target grid (m)

LOGICAL            :: L_endgame_out
LOGICAL            :: L_endgame_in
LOGICAL            :: l_outside

CHARACTER (LEN=*), PARAMETER     :: RoutineName = 'RCF_H_INT_INIT_CART_BL'
CHARACTER (LEN=errormessagelength)  :: Cmessage
TYPE( grid_type )  :: target_grid

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! 0: Allocate arrays depending on grid used.
!-----------------------------------------------------------------------------

L_endgame_in  = Hdr_In  % FixHd(FH_GridStagger) == FH_GridStagger_Endgame
L_endgame_out = Hdr_Out % FixHd(FH_GridStagger) == FH_GridStagger_Endgame

source_grid_type = Hdr_In  % FixHd(FH_HorizGrid)

! Note Cartesian dumps only recongised from UM10.9 onwards so all 
! cartesian dumps are ENDGame

IF ( L_endgame_in ) THEN
  glob_rows_in  = grid_in % glob_v_rows
  glob_field_in = grid_in % glob_v_field
ELSE
  glob_rows_in  = grid_in % glob_p_rows
  glob_field_in = grid_in % glob_p_field
END IF

ALLOCATE(x_in(glob_field_in), &
         y_in(glob_rows_in ))

IF ( L_endgame_out ) THEN
  glob_rows_out  = grid_out % glob_v_rows
  glob_field_out = grid_out % glob_v_field
ELSE
  glob_rows_out  = grid_out % glob_p_rows
  glob_field_out = grid_out % glob_p_field
END IF

ALLOCATE(x_out(glob_field_out), &
         y_out(glob_field_out))


! Initialise error tolerance for grid boundaries. For basic interpolating of 
! input grid only we should not set a grid_tol and instead indicate to 
! only warn if output grid is not inside input grid.
IF (l_basic_interp) THEN
  grid_tol = 0.0
ELSE
  grid_tol = SQRT(EPSILON(1.0))
END IF

!--------------------------------------------------------------
! 1: Test if requested horizontal interpolation is sensible
!--------------------------------------------------------------

! 1.1: Hemispheric or LAM -> global not allowed
IF ( Hdr_Out % FixHd(FH_HorizGrid) == 0 .AND.                   &
     Hdr_In % FixHd(FH_HorizGrid) > 0) THEN
  WRITE(umMessage,'(A,A)') ' *ERROR* Trying to interpolate from a ',          &
      'hemispheric or LAM to a global domain'
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  Cmessage = 'Cannot interpolation from hemisphere/LAM to global'
  ErrorStatus = 2
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! 1.2: LAM -> hemispheric not allowed
IF ( Hdr_Out % FixHd(FH_HorizGrid) < 3 .AND.                     &
     Hdr_In % FixHd(FH_HorizGrid) > 2) THEN
  WRITE(umMessage,'(A,A)')' *ERROR* Trying to interpolate from a ',         &
     'limited area domain to a global or hemispheric domain'
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  Cmessage = 'Cannot interpolate from LAM to global or hemisphere'
  ErrorStatus = 3
  CALL ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!--------------------------------------------------------------
! Begin setting up the target grid for interpolation:
!--------------------------------------------------------------
target_grid % glob_p_rows       = grid_out % glob_p_rows
target_grid % glob_p_row_length = grid_out % glob_p_row_length
target_grid % glob_p_field      = grid_out % glob_p_field
target_grid % glob_u_rows       = grid_out % glob_u_rows
target_grid % glob_u_row_length = grid_out % glob_u_row_length
target_grid % glob_u_field      = grid_out % glob_u_field
target_grid % glob_v_rows       = grid_out % glob_v_rows
target_grid % glob_v_row_length = grid_out % glob_v_row_length
target_grid % glob_v_field      = grid_out % glob_v_field

target_grid % phi_npole         = grid_out % phi_npole
target_grid % lambda_npole      = grid_out % lambda_npole


ALLOCATE(target_grid % lambda_p(target_grid % glob_p_row_length))
ALLOCATE(target_grid % lambda_u(target_grid % glob_u_row_length))
ALLOCATE(target_grid % phi_p(target_grid % glob_p_rows))
ALLOCATE(target_grid % phi_v(target_grid % glob_v_rows))


DO i = 1, target_grid % glob_p_row_length
  target_grid % lambda_p(i) = grid_out % lambda_p(i)
END DO
DO i = 1, target_grid % glob_u_row_length
  target_grid % lambda_u(i) = grid_out % lambda_u(i)
END DO
DO i = 1, target_grid % glob_p_rows
  target_grid % phi_p(i)    = grid_out % phi_p(i)
END DO
DO i = 1, target_grid % glob_v_rows
  target_grid % phi_v(i)    = grid_out % phi_v(i)
END DO

!-------------------------------------------------------------------
! Y and X  of P grid 
!-------------------------------------------------------------------
! New target grid create full field version of coordinates.
ij = 0
DO j = 1, target_grid % glob_p_rows
  DO i = 1, target_grid % glob_p_row_length
    ij = ij + 1
    x_out(ij) = target_grid % lambda_p(i)
    y_out(ij) = target_grid % phi_p(j)
  END DO
END DO

! Y and X of old grid just for rows and columns respectively

DO j = 1, grid_in % glob_p_rows
  y_in(j) = grid_in % phi_p(j)
END DO
DO j = 1, grid_in % glob_p_row_length
  x_in(j) = grid_in % lambda_p(j)
END DO

!-------------------------------------------------------------------
! 3.7: Check if target area contained within source area when
!      LAM->LAM.
!-------------------------------------------------------------------
IF (Hdr_Out % FixHd(FH_HorizGrid) /= 0 .AND.                                   &
    Hdr_In  % FixHd(FH_HorizGrid) /= 0) THEN


!  ENDGame grid    s - bottom left hand corner
!
!     -----V---------V---------V---------V--e  <- Top right hand corner of grid
!                                                 call this e for end. V if not
!     U    P    U    P    U    P    U    P        bicyclic
!                                           
!          V         V         V         V
!
! |   U    P    U    P    U    P    U    P
! |   
! |        V         V         V         V
! |
! |   U    P    U    P    U    P    U    P
! |   
! |   s    V         V         V         V
! y, X ----->

  ! Want out s and e locations
  s_y_out = target_grid%phi_v(1)
  s_x_out = target_grid%lambda_u(1)
  e_y_out = s_y_out + target_grid%glob_p_rows*target_delta_lat
  e_x_out = s_x_out + target_grid%glob_p_row_length*target_delta_lon

  ! Want in and out s and e locations
  s_y_in = grid_in%phi_v(1)
  s_x_in = grid_in%lambda_u(1)
  dy = grid_in%phi_p(2)-grid_in%phi_p(1)
  dx = grid_in%lambda_p(2)-grid_in%lambda_p(1)
  e_y_in = s_y_in + grid_in%glob_p_rows*dy
  e_x_in = s_x_in + grid_in%glob_p_row_length*dx

  l_outside = .FALSE.

  IF ( s_y_out < s_y_in - grid_tol  ) THEN
    l_outside = .TRUE.
    WRITE(umMessage,'(A,f12.1,A,f12.1)') "Target grid bottom ", s_y_out,&
               " below source grid bottom", s_y_in
    CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')
    
  END IF 
  IF ( s_x_out < s_x_in - grid_tol  ) THEN
    l_outside = .TRUE.
    WRITE(umMessage,'(A,f12.1,A,f12.1)') "Target grid left ", s_x_out,&
               " to left of source grid ", s_x_in
    CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')
  END IF 
  IF ( e_y_out > e_y_in + grid_tol  ) THEN
    l_outside = .TRUE.
    WRITE(umMessage,'(A,f12.1,A,f12.1)') "Target grid top ", e_y_out,&
               " above top of source grid ", e_y_in
    CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')
  END IF 
  IF ( e_x_out > e_x_in + grid_tol  ) THEN
    l_outside = .TRUE.
    WRITE(umMessage,'(A,f12.1,A,f12.1)') "Target grid right ", e_x_out,&
               " to right of source grid ", e_x_in
    CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')
  END IF 

  ! New grid not inside old grid
  IF (l_outside) THEN
    cmessage = "target grid outside source grid see output for more details "
    errorstatus = 1
    CALL ereport( routinename, errorstatus, cmessage )
  END IF
 
END IF

!----------------------------------------------------------------------------
! P grid points  - weights and indices required for bi-linear interpolation
! and coastal adjustment
!----------------------------------------------------------------------------
WRITE(umMessage,'(A)') "P grid cartesian "
CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')

CALL f_shum_cart_horizontal_field_bi_lin_interp_get_coeffs                    &
 (bl_index_b_l(1,1), bl_index_b_r(1,1), bl_index_t_l(1,1), bl_index_t_r(1,1), &
  weight_t_r(1,1), weight_b_r(1,1), weight_t_l(1,1), weight_b_l(1,1),         &
  x_in, y_in, x_out, y_out,                                                   &
  grid_in % glob_p_row_length, grid_in % glob_p_rows,                         &
  target_grid % glob_p_field, source_grid_type, icode, cmessage)

IF (icode > 0) THEN
  errorstatus = icode
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!----------------------------------------------------------------------------
! U grid reset grid info passed to h_int_co_cart
!----------------------------------------------------------------------------
ij=0
DO j=1,target_grid % glob_u_rows
  DO i=1,target_grid % glob_u_row_length
    ij=ij+1
    x_out(ij) = target_grid % lambda_u(i)
    y_out(ij) = target_grid % phi_p(j)
  END DO
END DO

! Y and X of source grid
DO j = 1, grid_in % glob_u_rows
  y_in(j) = grid_in % phi_p(j)
END DO
DO i = 1, grid_in % glob_u_row_length
  x_in(i) = grid_in % lambda_u(i)
END DO
WRITE(umMessage,'(A)') "U grid cartesian "
CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')

CALL f_shum_cart_horizontal_field_bi_lin_interp_get_coeffs                    &
 (bl_index_b_l(1,2), bl_index_b_r(1,2), bl_index_t_l(1,2), bl_index_t_r(1,2), &
  weight_t_r(1,2), weight_b_r(1,2), weight_t_l(1,2), weight_b_l(1,2),         &
  x_in, y_in, x_out, y_out,                                                   &
  grid_in % glob_p_row_length, grid_in % glob_p_rows,                         &
  target_grid % glob_p_field, source_grid_type, icode, cmessage)

IF (icode > 0) THEN
  errorstatus = icode
  CALL ereport( routinename, errorstatus, cmessage )
END IF


!----------------------------------------------------------------------------
!  V grid  - bicyclic has same number of rows as U and P
!            channel and oridinary LAM have an extra row
!----------------------------------------------------------------------------

ij=0
DO j=1,target_grid % glob_v_rows
  DO i=1,target_grid % glob_v_row_length
    ij=ij+1
    x_out(ij) = target_grid % lambda_p(i)
    y_out(ij) = target_grid % phi_v(j)
  END DO
END DO

! Y and X of source grid
DO j = 1, grid_in % glob_v_rows
  y_in(j) = grid_in % phi_v(j)
END DO
DO i = 1, grid_in % glob_p_row_length
  x_in(i) = grid_in % lambda_p(i)
END DO

WRITE(umMessage,'(A)') "V grid cartesian "
CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')

CALL f_shum_cart_horizontal_field_bi_lin_interp_get_coeffs                    &
 (bl_index_b_l(1,3), bl_index_b_r(1,3), bl_index_t_l(1,3), bl_index_t_r(1,3), &
  weight_t_r(1,3), weight_b_r(1,3), weight_t_l(1,3), weight_b_l(1,3),         &
  x_in, y_in, x_out, y_out,                                                   &
  grid_in % glob_p_row_length, grid_in % glob_v_rows,                         &
  target_grid % glob_v_field, source_grid_type, icode, cmessage)

IF (icode > 0) THEN
  errorstatus = icode
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!----------------------------------------------------------------------------
! Weights and indices for zonal mean P points:
! Code as for lat-lon case, expected to work but NOT been tested.
! Appears to only be used if grid type is that set for ozone
!----------------------------------------------------------------------------

! Y and X of target grid
x_out(1) = target_grid % lambda_p(1)
DO j=1, target_grid % glob_p_rows
  y_out(j)  = target_grid % phi_p(j)
END DO


! Y and X of source grid
x_in(1) = grid_in % lambda_p(1)
DO j=1, grid_in % glob_p_rows
  y_in(j) = grid_in % phi_p(j)
END DO

WRITE(umMessage,'(A)') "Zonal P grid cartesian "
CALL umPrint(umMessage,src='rcf_h_int_init_cart_bl')

CALL f_shum_cart_horizontal_field_bi_lin_interp_get_coeffs                    &
 (bl_index_b_l(1,4), bl_index_b_r(1,4), bl_index_t_l(1,4), bl_index_t_r(1,4), &
  weight_t_r(1,4), weight_b_r(1,4), weight_t_l(1,4), weight_b_l(1,4),         &
  x_in, y_in, x_out, y_out,                                                   &
  grid_in % glob_p_row_length, grid_in % glob_p_rows,                         &
  target_grid % glob_p_field, source_grid_type, icode, cmessage)

IF (icode > 0) THEN
  errorstatus = icode
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!--------------------------------------------------------------
! Cleanup - deallocate arrays
!--------------------------------------------------------------
IF (ASSOCIATED(target_grid % phi_v))    DEALLOCATE(target_grid % phi_v)
IF (ASSOCIATED(target_grid % phi_p))    DEALLOCATE(target_grid % phi_p)
IF (ASSOCIATED(target_grid % lambda_u)) DEALLOCATE(target_grid % lambda_u)
IF (ASSOCIATED(target_grid % lambda_p)) DEALLOCATE(target_grid % lambda_p)


IF (ALLOCATED(y_in))     DEALLOCATE(y_in)
IF (ALLOCATED(y_out))    DEALLOCATE(y_out)

IF (ALLOCATED(x_in))     DEALLOCATE(x_in)
IF (ALLOCATED(x_out))    DEALLOCATE(x_out)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_h_int_init_cart_bl
END MODULE rcf_h_int_init_cart_bl_mod
