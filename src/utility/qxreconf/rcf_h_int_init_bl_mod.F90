! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialises bilinear interpolation weights

MODULE Rcf_H_Int_Init_BL_Mod

!  Subroutine Rcf_H_Int_Init_BL
!
! Description:
!   Sets up the interpolation and rotation weights for bilinear
!   horizontal interpolation.
!
! Method:
!   Weights stored in the interp_weights_mod module
!   Seperate calculations for P, U, V and P zonal points.
!
!   Based (with many changes) on UM4.5 code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_H_INT_INIT_BL_MOD'

CONTAINS

SUBROUTINE Rcf_H_Int_Init_BL( grid_in, grid_out, hdr_in, hdr_out )

USE Rcf_Interp_Weights_Mod            ! All of it (well almost)

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Headaddress_Mod, ONLY:       &
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
                               f_shum_horizontal_field_bi_lin_interp_get_coeffs
USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon, rotate_latlon_to_eq, &
                                  eq_latlon_vector_coeffs

USE lam_inclusion_mod, ONLY: lam_inclusion
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

USE lam_config_inputs_mod, ONLY: &
    target_delta_lat,            &
    target_delta_lon,            &
    target_frstlata,             &
    target_frstlona

USE fort2c_interfaces, ONLY: logical_to_bool

IMPLICIT NONE

! Arguments
TYPE( grid_type ),      INTENT(IN)   :: grid_in
TYPE( grid_type ),      INTENT(IN)   :: grid_out
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
TYPE( um_header_type ), INTENT(IN)   :: hdr_out

! Local variables
INTEGER            :: i                ! looper
INTEGER            :: j                ! looper
INTEGER            :: ij               ! looper
INTEGER            :: p_points_out     ! Total points in output grid
INTEGER            :: ErrorStatus
REAL               :: grid_tol         ! Grid tolerance between grids.
REAL               :: max_phi_out      ! Maximum output phi
REAL               :: min_phi_out      ! Minimum output phi
REAL               :: lambda_out_row( grid_out % glob_p_row_length )
REAL               :: phi_out_col   ( grid_out % glob_p_rows       )
REAL               :: lambda_in( grid_in % glob_p_row_length )
                                    ! Longitude coords of source p-grid
REAL, ALLOCATABLE  :: lambda_inn( : )
                                    ! Longitude coords of source p-grid
REAL, ALLOCATABLE  :: lambda_out( : )
                                    ! Longitude coords of target grid
REAL, ALLOCATABLE  :: lambda_rot( : )
                                    ! Standard lon coords of source grid
REAL, ALLOCATABLE  :: lambda_tmp( : )
                                    ! Standard lon coords of target grid
REAL, ALLOCATABLE  :: phi_in( : )
                                    ! Latitude coords of source p-grid
REAL, ALLOCATABLE  :: phi_inn( : )
                                    ! Latitude coords of source p-grid
REAL, ALLOCATABLE  :: phi_out( : )
                                    ! Latitude coords of target grid
REAL, ALLOCATABLE  :: phi_tmp( : )
                                    ! Standard lat coords of target grid
LOGICAL            :: Rot_In
LOGICAL            :: Rot_Out
LOGICAL            :: Cyclic
LOGICAL            :: L_endgame_out
LOGICAL            :: L_endgame_in
INTEGER            :: glob_rows_in, glob_field_in
INTEGER            :: glob_rows_out, glob_field_out
CHARACTER (LEN=*), PARAMETER     :: RoutineName = 'RCF_H_INT_INIT_BL'
CHARACTER (LEN=errormessagelength)  :: Cmessage
TYPE( grid_type )  :: target_grid
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!---
! 0: Allocate arrays depending on grid used.
!---
L_endgame_in  = Hdr_In  % FixHd(FH_GridStagger) == FH_GridStagger_Endgame
L_endgame_out = Hdr_Out % FixHd(FH_GridStagger) == FH_GridStagger_Endgame

IF ( L_endgame_in ) THEN
  glob_rows_in  = grid_in % glob_v_rows
  glob_field_in = grid_in % glob_v_field
ELSE
  glob_rows_in  = grid_in % glob_p_rows
  glob_field_in = grid_in % glob_p_field
END IF

ALLOCATE(lambda_inn(glob_field_in), &
         lambda_rot(glob_field_in), &
         phi_in    (glob_rows_in ),  &
         phi_inn   (glob_field_in))

IF ( L_endgame_out ) THEN
  glob_rows_out  = grid_out % glob_v_rows
  glob_field_out = grid_out % glob_v_field
ELSE
  glob_rows_out  = grid_out % glob_p_rows
  glob_field_out = grid_out % glob_p_field
END IF


ALLOCATE(lambda_out(glob_field_out), &
         lambda_tmp(glob_field_out), &
         phi_out   (glob_field_out), &
         phi_tmp   (glob_field_out))


! Initialise error tolerance for grid boundaries. For basic interpolating of 
! input grid only we should not set a grid_tol and instead indicate to 
! only warn if output grid is not inside input grid.
IF (l_basic_interp) THEN
  grid_tol = rmdi
ELSE
  grid_tol = SQRT(EPSILON(1.0))
END IF

!--------------------------------------------------------------
! 1: Test if requested horizontal interpolation is sensible
!--------------------------------------------------------------

! 1.1: Hemispheric or LAM -> global not allowed
IF ( Hdr_Out % FixHd(FH_HorizGrid) == 0 .AND.                   &
     Hdr_In % FixHd(FH_HorizGrid) > 0) THEN
  WRITE(umMessage,'(A,A)')' *ERROR* Trying to interpolate from a hemispheric ',&
      ' or LAM to a global domain'
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  Cmessage = 'Cannot interpolation from hemisphere/LAM to global'
  ErrorStatus = 2
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! 1.2: LAM -> hemispheric not allowed
IF ( Hdr_Out % FixHd(FH_HorizGrid) < 3 .AND.                     &
     Hdr_In % FixHd(FH_HorizGrid) > 2) THEN
  WRITE(umMessage,'(A,A)')' *ERROR* Trying to interpolate from a limited area ',&
      'domain to a global or hemispheric domain'
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  CALL umPrint(umMessage,src='rcf_h_int_init_bl_mod')
  Cmessage = 'Cannot interpolate from LAM to global or hemisphere'
  ErrorStatus = 3
  CALL ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! 2.4: Logical to indicate output grid is rotated
! 2.5: Logical to indicate input grid is rotated
! If rotations are the same no need to rotate grid if told to do so with
! l_limit_rotations.
IF ( l_same_rotation ) THEN
  rot_out = .FALSE.
  rot_in  = .FALSE.
ELSE
  rot_out = grid_out % phi_npole /= 90.0 .OR. grid_out % lambda_npole /= 0.0
  rot_in  = grid_in  % phi_npole /= 90.0 .OR. grid_in  % lambda_npole /= 0.0
END IF

! 2.6: Logical to indicate if input data cyclic
cyclic = Hdr_In % FixHd( FH_HorizGrid ) < 3

! To interpolate from lat-lon to Cartesian we construct a lat-lon grid
! over the domain specified by the original lam_config namelist input, and at 
! the resolution of the intended Cartesian output grid. This represents a
! projection of the intended Cartesian grid onto specified area of the globe.
! This works best for smaller/equatorial domains.

! Begin setting up the target grid for interpolation:
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

IF (hdr_out % fixhd( FH_CoordSystem ) == FH_CoordSystem_Cartesian) THEN
  ! Cartesian output grid: set up a lat-lon grid using the original output
  ! co-ordinates and at the intended resolution.
  DO i = 1, target_grid % glob_p_row_length
    target_grid % lambda_p(i) = target_frstlona + (i-0.5)*target_delta_lon
  END DO
  DO i = 1, target_grid % glob_u_row_length
    target_grid % lambda_u(i) = target_frstlona + (i-1)*target_delta_lon
  END DO
  DO i = 1, target_grid % glob_p_rows
    target_grid % phi_p(i)    = target_frstlata + (i-0.5)*target_delta_lat
  END DO
  DO i = 1, target_grid % glob_v_rows
    target_grid % phi_v(i)    = target_frstlata + (i-1)*target_delta_lat
  END DO
ELSE
  ! Non-Cartesian grid: copy the existing output grid.
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
END IF

!-------------------------------------------------------------------
! 3. Weights and indices for P points
!-------------------------------------------------------------------
! 3.1: Lat and lon of new grid

DO j=1,target_grid % glob_p_rows
  phi_out_col(j)    = target_grid % phi_p(j)
END DO
DO i=1,target_grid % glob_p_row_length
  lambda_out_row(i) = target_grid % lambda_p(i)
END DO
lambda_out_row(1) = MODULO(lambda_out_row(1), 360.0)
WHERE (lambda_out_row < lambda_out_row(1)) &
  lambda_out_row = lambda_out_row + 360.0

! Now create full field version of coordinates.
ij = 0
DO j = 1, target_grid % glob_p_rows
  DO i = 1, target_grid % glob_p_row_length
    ij = ij + 1
    lambda_out(ij) = target_grid % lambda_p(i)
    phi_out(ij)    = target_grid % phi_p(j)
  END DO
END DO

! 3.2: Convert output grid to standard lat-lon if output grid rotated
IF (rot_out) THEN
  CALL rotate_eq_to_latlon(phi_out, lambda_out, phi_tmp, lambda_tmp,      &
                 target_grid % phi_npole, target_grid % lambda_npole,     &
                 target_grid % glob_p_field )

  ! 3.2.1: Calculate coefficients for rotating wind
  CALL eq_latlon_vector_coeffs(coeff1, coeff2, lambda_tmp, lambda_out,    &
                 target_grid % phi_npole, target_grid % lambda_npole,     &
                 target_grid % glob_p_field )

  ! 3.2.2: Copy across latitude and longitude of target grid in standard
  !        coords for next step
  DO i = 1, target_grid % glob_p_field
    lambda_out(i) = lambda_tmp(i)
    phi_out(i)    = phi_tmp(i)
  END DO
END IF

! 3.3: Convert output grid to input lat-lon if input grid rotated
! This already makes sure output is between 0 and 360.
IF (rot_in) THEN
  CALL rotate_latlon_to_eq(phi_out, lambda_out, phi_out, lambda_out,    &
                           grid_in % phi_npole, grid_in % lambda_npole, &
                           target_grid % glob_p_field )
END IF

! 3.4: rotate_latlon_to_eq ensures that lambda_out is between 0 and 360

! 3.5: Lat and lon of old grid
phi_in   (1:grid_in % glob_p_rows)          =                    &
  grid_in % phi_p(1:grid_in % glob_p_rows)
lambda_in(1:grid_in % glob_p_row_length)    =                    &
  grid_in % lambda_p(1:grid_in % glob_p_row_length)

! 3.6 rotate_latlon_to_eq already makes sure lambda_out is between 0->360 but
! lambda_in is left alone since we shouldn't need to modify the source grid.

! 3.7: Check if target area contained within source area when
!      LAM->LAM.
IF (Hdr_Out % FixHd(FH_HorizGrid) /= 0 .AND.                                   &
    Hdr_In  % FixHd(FH_HorizGrid) /= 0) THEN

  ! Loop through output grid on previously populated points only
  ! and not full allocated size
  min_phi_out = CEILING(phi_out(1))
  max_phi_out = FLOOR(phi_out(1))
  DO i = 1, target_grid % glob_p_field
    min_phi_out = MIN(min_phi_out, phi_out(i))
    max_phi_out = MAX(max_phi_out, phi_out(i))
  END DO

  CALL lam_inclusion(grid_in % glob_p_rows, grid_in % glob_p_row_length, 0, 0, &
    lambda_in, phi_in, grid_in  % phi_npole, grid_in  % lambda_npole,          &
    target_grid % glob_p_rows, target_grid % glob_p_row_length, 0, 0,          &
    lambda_out_row, phi_out_col, min_phi_out, max_phi_out,                     &
    target_grid % phi_npole, target_grid % lambda_npole,                       &
    l_same_rotation, rot_in, grid_tol)
END IF

! Calculate coefficients.  Need to create temporary arrays since we need to
! rotate grid to an unrotated pole.
IF (ROT_in) THEN
  ij=0
  DO j=1,grid_in % glob_p_rows
    DO i=1,grid_in % glob_p_row_length
      ij=ij+1
      lambda_inn(ij) = grid_in % lambda_p(i)
      phi_inn(ij)    = grid_in % phi_p(j)
    END DO
  END DO
  CALL rotate_eq_to_latlon(phi_inn, lambda_inn, phi_inn,lambda_rot,     &
                           grid_in % phi_npole, grid_in % lambda_npole, &
                           grid_in % glob_p_field)

  ! Calculate coefficients for rotating wind
  CALL eq_latlon_vector_coeffs(coeff3, coeff4, lambda_rot, lambda_inn,      &
                               grid_in % phi_npole, grid_in % lambda_npole, &
                               grid_in % glob_p_field)
END IF

! 3.8: Initialise Indices and weights for Bi-linear interpolation
!        As these are required for coastal adjustment, this must be
!        called whether bi-linear interpolation is requested or not.

! 3.8.1: Call H_INT_CO

CALL f_shum_horizontal_field_bi_lin_interp_get_coeffs                      &
                      (bl_index_b_l(1,1), bl_index_b_r(1,1),               &
                       bl_index_t_l(1,1), bl_index_t_r(1,1),               &
                       weight_t_r(1,1), weight_b_r(1,1),                   &
                       weight_t_l(1,1), weight_b_l(1,1),                   &
                       lambda_in,phi_in, lambda_out, phi_out,              &
                       grid_in % glob_p_row_length, grid_in % glob_p_rows, &
                       target_grid % glob_p_field, logical_to_bool(cyclic))

!-------------------------------------------------------------------
! U first
!-------------------------------------------------------------------
ij=0
DO j=1,target_grid % glob_u_rows
  DO i=1,target_grid % glob_u_row_length
    ij=ij+1
    lambda_out(ij)=target_grid % lambda_u(i)
    phi_out(ij)   =target_grid % phi_p(j)
  END DO
END DO

! 4.2: Lat and lon of source grid
DO j = 1, grid_in % glob_u_rows
  phi_in(j)    = grid_in % phi_p(j)
END DO
DO i = 1, grid_in % glob_u_row_length
  lambda_in(i) = grid_in % lambda_u(i)
END DO

! 4.3: Convert to standard lat-lon if new grid rotated
IF (rot_out) THEN
  CALL rotate_eq_to_latlon(phi_out, lambda_out, phi_out, lambda_out,    &
              target_grid % phi_npole, target_grid % lambda_npole,      &
              target_grid % glob_u_field )
END IF

! 4.6: Convert target grid to source lat-lon if source grid rotated
IF (rot_in) THEN
  CALL rotate_latlon_to_eq(phi_out, lambda_out, phi_out, lambda_out,    &
                           grid_in % phi_npole, grid_in % lambda_npole, &
                           target_grid % glob_u_field )
END IF

! 4.7 rotate_ routines already make sure lambda_out is between 0->360 but
! lambda_in is left alone since we shouldn't need to modify the source grid.

! 4.9: Indices and weights for horizontal interpolation
CALL f_shum_horizontal_field_bi_lin_interp_get_coeffs                      &
                              ( bl_index_b_l(1,2), bl_index_b_r(1,2),      &
                                bl_index_t_l(1,2), bl_index_t_r(1,2),      &
                                weight_t_r(1,2), weight_b_r(1,2),          &
                                weight_t_l(1,2), weight_b_l(1,2),          &
                                lambda_in, phi_in, lambda_out, phi_out,    &
                                grid_in % glob_u_row_length,               &
                                grid_in % glob_u_rows,                     &
                                target_grid % glob_u_field,                &
                                logical_to_bool(cyclic))

!-------------------------------------------------------------------
!  Then V
!-------------------------------------------------------------------
ij=0
DO j=1,target_grid % glob_v_rows
  DO i=1,target_grid % glob_v_row_length
    ij=ij+1
    lambda_out(ij) = target_grid % lambda_p(i)
    phi_out(ij)    = target_grid % phi_v(j)
  END DO
END DO

! 4.2: Lat and lon of source grid
phi_in   (1:grid_in % glob_v_rows)             =             &
  grid_in % phi_v(1:grid_in % glob_v_rows)
lambda_in(1:grid_in % glob_p_row_length)       =             &
  grid_in % lambda_p(1:grid_in % glob_v_row_length)

! 4.3: Convert to standard lat-lon if new grid rotated
IF (rot_out) THEN
  CALL rotate_eq_to_latlon(phi_out, lambda_out, phi_out, lambda_out,    &
              target_grid % phi_npole, target_grid % lambda_npole,      &
              target_grid % glob_v_field )
END IF

! 4.6: Convert target grid to source lat-lon if source grid rotated
IF (rot_in) THEN
  CALL rotate_latlon_to_eq(phi_out, lambda_out, phi_out, lambda_out,    &
                           grid_in % phi_npole, grid_in % lambda_npole, &
                           target_grid % glob_v_field )
END IF

! 4.7 rotate_ routines already make sure lambda_out is between 0->360 but
! lambda_in is left alone since we shouldn't need to modify the source grid.
! 4.9: Indices and weights for horizontal interpolation
CALL f_shum_horizontal_field_bi_lin_interp_get_coeffs                      &
                              ( bl_index_b_l(1,3), bl_index_b_r(1,3),      &
                                bl_index_t_l(1,3), bl_index_t_r(1,3),      &
                                weight_t_r(1,3), weight_b_r(1,3),          &
                                weight_t_l(1,3), weight_b_l(1,3),          &
                                lambda_in, phi_in, lambda_out, phi_out,    &
                                grid_in % glob_v_row_length,               &
                                grid_in % glob_v_rows,                     &
                                target_grid % glob_v_field,                &
                                logical_to_bool(cyclic))

!--------------------------------------------------------------
! 5: Weights and indices for zonal mean P points:
!--------------------------------------------------------------
! 5.1: Lat and lon of target grid
lambda_out(1) = target_grid % lambda_p(1)
DO j=1, target_grid % glob_p_rows
  phi_out(j)  = target_grid % phi_p(j)
END DO

p_points_out = target_grid % glob_p_rows

! 5.2: Lat and lon of source grid
lambda_in(1) = grid_in % lambda_p(1)
phi_in(1:grid_in%glob_p_rows)             = &
  grid_in % phi_p(1:grid_in%glob_p_rows)

! 5.3: Convert output grid to standard lat-lon if target grid rotated
IF (rot_out) THEN
  CALL rotate_eq_to_latlon(phi_out, lambda_out, phi_out, lambda_out,    &
              target_grid % phi_npole, target_grid % lambda_npole,      &
              p_points_out )
END IF

! 5.4: Convert target grid to source lat-lon if source grid rotated
IF (rot_in) THEN
  CALL rotate_latlon_to_eq(phi_out, lambda_out, phi_out, lambda_out,    &
                           grid_in % phi_npole, grid_in % lambda_npole, &
                           p_points_out )
END IF

! 5.5: Initialise indices and weights for horizontal interpolation
CALL f_shum_horizontal_field_bi_lin_interp_get_coeffs                      &
                        ( bl_index_b_l(1,4), bl_index_b_r(1,4),            &
                          bl_index_t_l(1,4), bl_index_t_r(1,4),            &
                          weight_t_r(1,4), weight_b_r(1,4),                &
                          weight_t_l(1,4), weight_b_l(1,4),                &
                          lambda_in, phi_in, lambda_out, phi_out,          &
                          1, grid_in % glob_p_rows, p_points_out,          &
                          logical_to_bool(cyclic))

IF (ASSOCIATED(target_grid % phi_v))    DEALLOCATE(target_grid % phi_v)
IF (ASSOCIATED(target_grid % phi_p))    DEALLOCATE(target_grid % phi_p)
IF (ASSOCIATED(target_grid % lambda_u)) DEALLOCATE(target_grid % lambda_u)
IF (ASSOCIATED(target_grid % lambda_p)) DEALLOCATE(target_grid % lambda_p)

IF (ALLOCATED(lambda_inn)) DEALLOCATE(lambda_inn)
IF (ALLOCATED(lambda_rot)) DEALLOCATE(lambda_rot)
IF (ALLOCATED(lambda_tmp)) DEALLOCATE(lambda_tmp)
IF (ALLOCATED(phi_in))     DEALLOCATE(phi_in)
IF (ALLOCATED(phi_inn))    DEALLOCATE(phi_inn)
IF (ALLOCATED(phi_out))    DEALLOCATE(phi_out)
IF (ALLOCATED(phi_tmp))    DEALLOCATE(phi_tmp)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_h_int_init_bl
END MODULE rcf_h_int_init_bl_mod
