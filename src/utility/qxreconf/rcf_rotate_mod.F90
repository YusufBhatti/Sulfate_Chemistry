! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Rotates winds for LAM dumps

MODULE Rcf_Rotate_Mod

!  Subroutine Rcf_Rotate - rotates winds
!
! Description:
! This routine reads u/v (and u_adv/v_adv if available) and
! rotates them to the standard orientation (Mode = ToStandard)
! or from the standard orientation (Mode = FromStandard).
!
! Method:
! The algorithm is a little different from UM <= 4.5,
! rather interpolating u/v onto theta points, then calculating
! the rotated winds before interpolating back.
! Uses same gather/scatter as horizontal interpolation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE latlon_eq_rotation_mod, ONLY: vector_eq_to_latlon, vector_latlon_to_eq

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! magic numbers
INTEGER, PARAMETER        :: ToStandard   = 0
INTEGER, PARAMETER        :: FromStandard = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ROTATE_MOD'

CONTAINS

SUBROUTINE Rcf_Rotate( fields, field_count, grid, hdr, decomp, mode )

USE Rcf_Gather_Field_Mod, ONLY: &
    Rcf_Gather_Field_Real

USE Rcf_Scatter_Field_Mod, ONLY: &
    Rcf_Scatter_Field_Real

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE UM_ParVars, ONLY:  &
    gc_all_proc_group,  &
    current_decomp_type,&
    change_decomposition

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE Rcf_Interp_Weights_Mod, ONLY: &
    Coeff1,                    &
    Coeff2,                    &
    Coeff3,                    &
    Coeff4

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_read_field_mod, ONLY: &
    Rcf_read_field

USE Rcf_write_field_mod, ONLY: &
    Rcf_write_field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE rcf_AC_interp_mod, ONLY: &
    rcf_p_to_v,               &
    rcf_p_to_u,               &
    rcf_u_to_p,               &
    rcf_v_to_p

USE nlstcall_mod, ONLY:   &
    LTimer

USE um_stashcode_mod, ONLY: &
    stashcode_u,              stashcode_v,    &
    stashcode_u_adv,          stashcode_v_adv,&
    stashcode_prog_sec

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
INTEGER, INTENT(IN)                :: mode

! Local variables
INTEGER                            :: i
INTEGER                            :: j
INTEGER                            :: upos
INTEGER                            :: vpos
INTEGER                            :: pe
INTEGER                            :: rem
INTEGER                            :: div
INTEGER                            :: p_size
INTEGER                            :: outer
INTEGER                            :: decomp_original
INTEGER                            :: ErrorStatus
LOGICAL                            :: l_mdi = .FALSE. ! check for MDI
REAL, ALLOCATABLE                  :: u_level(:,:)
REAL, ALLOCATABLE                  :: v_level(:,:)
REAL, ALLOCATABLE                  :: u_at_p (:,:)
REAL, ALLOCATABLE                  :: v_at_p (:,:)
REAL, ALLOCATABLE                  :: u_at_p_out (:,:)
REAL, ALLOCATABLE                  :: v_at_p_out (:,:)

TYPE( field_type ), POINTER        :: u_field
TYPE( field_type ), POINTER        :: v_field


CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_ROTATE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( 'Rcf_Rotate', 3)

p_size = grid % glob_p_rows * grid % glob_p_row_length;

decomp_original = current_decomp_type
CALL Change_Decomposition( decomp )

! We need to do this all twice - once for u/v and then for
! u_adv/u_adv
DO outer = 1, 2

  !-----------------------------------------------------------------
  ! Locate the u and v fields, allocate space for the data and
  ! read in the field from file
  !-----------------------------------------------------------------
    ! find the u-field (u/u_adv)
  IF ( outer == 1) THEN
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_u,                &
                     fields, field_count, upos, zero_ok_arg = .TRUE. )
  ELSE     ! outer = 2
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_u_adv,            &
                     fields, field_count, upos, zero_ok_arg = .TRUE.)
  END IF

  IF ( upos == 0 ) THEN
    CYCLE;
  END IF

  ! find the v-field (v/v_adv)
  IF ( outer == 1) THEN
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_v,                &
                     fields, field_count, vpos, zero_ok_arg = .TRUE. )
  ELSE     ! outer = 2
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_v_adv,            &
                     fields, field_count, vpos, zero_ok_arg = .TRUE.)
  END IF

  IF ( vpos == 0 ) THEN
    CYCLE;
  END IF

  ! Allocate space and read in fields
  u_field => fields( upos )
  CALL Rcf_Alloc_Field( u_field )
  CALL Rcf_Read_Field( u_field, hdr, decomp )

  v_field => fields( vpos )
  CALL Rcf_Alloc_Field( v_field )
  CALL Rcf_Read_Field( v_field, hdr, decomp )

  !-----------------------------------------------------------------
  ! Gather the fields - one per pe - for the interpolation/rotation
  !-----------------------------------------------------------------
  ALLOCATE( u_level( grid % glob_u_row_length, grid % glob_u_rows) );
  ALLOCATE( v_level( grid % glob_v_row_length, grid % glob_v_rows) );
  ALLOCATE( u_at_p ( grid % glob_p_row_length, grid % glob_p_rows) );
  ALLOCATE( v_at_p ( grid % glob_p_row_length, grid % glob_p_rows) );
  ALLOCATE( u_at_p_out( grid % glob_p_row_length, grid % glob_p_rows) );
  ALLOCATE( v_at_p_out( grid % glob_p_row_length, grid % glob_p_rows) );

  div = u_field % levels / nproc
  rem = MOD( u_field % levels, nproc )
  pe = 0

  DO i = 1, div
    DO j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on pe
      ! Cannot use generic subroutine as u_level is 2D

      CALL Rcf_Gather_Field_Real( u_field % DATA(:,j), u_level,    &
                                  u_field % row_len,               &
                                  u_field % rows,                  &
                                  u_field % glob_row_len,          &
                                  u_field % glob_rows, pe,         &
                                  gc_all_proc_group )

      CALL Rcf_Gather_Field_Real( v_field % DATA(:,j), v_level,    &
                                  v_field % row_len,               &
                                  v_field % rows,                  &
                                  v_field % glob_row_len,          &
                                  v_field % glob_rows, pe,         &
                                  gc_all_proc_group )

      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO

    CALL Rcf_U_to_P( u_level, u_at_p, grid )
    CALL Rcf_V_to_P( v_level, v_at_p, grid )

    IF ( Mode == ToStandard ) THEN
      CALL vector_eq_to_latlon( Coeff3, Coeff4, u_at_p, v_at_p,  &
                                u_at_p_out, v_at_p_out, p_size, l_mdi)
    ELSE IF ( Mode == FromStandard ) THEN
      CALL vector_latlon_to_eq( Coeff1, Coeff2, u_at_p, v_at_p,  &
                                u_at_p_out, v_at_p_out, p_size, l_mdi)
    END IF

    CALL Rcf_P_to_U( u_at_p_out, u_level, grid )
    CALL Rcf_P_to_V( v_at_p_out, v_level, grid )

    DO j = ((i-1) * nproc) + 1, i * nproc

      CALL Rcf_Scatter_Field_Real( u_field % DATA(:,j), u_level,   &
                                   u_field % row_len,              &
                                   u_field % rows,                 &
                                   u_field % glob_row_len,         &
                                   u_field % glob_rows,pe,         &
                                   gc_all_proc_group )

      CALL Rcf_Scatter_Field_Real( v_field % DATA(:,j), v_level,   &
                                   v_field % row_len,              &
                                   v_field % rows,                 &
                                   v_field % glob_row_len,         &
                                   v_field % glob_rows,pe,         &
                                   gc_all_proc_group )

      pe = pe + 1

      IF (pe == nproc) pe = 0
    END DO
  END DO

  !-------------------------------------------------------------------
  ! There are rem levels left to process. Will do these now.
  !-------------------------------------------------------------------
  pe = 0
  DO i = 1, rem
    j = nproc * div + i

    CALL Rcf_Gather_Field_Real( u_field % DATA(:,j), u_level,         &
                                u_field % row_len,                    &
                                u_field % rows,                       &
                                u_field % glob_row_len,               &
                                u_field % glob_rows, pe,              &
                                gc_all_proc_group )

    CALL Rcf_Gather_Field_Real( v_field % DATA(:,j), v_level,         &
                                v_field % row_len,                    &
                                v_field % rows,                       &
                                v_field % glob_row_len,               &
                                v_field % glob_rows, pe,              &
                                gc_all_proc_group )

    pe = pe + 1
  END DO

  IF (mype < pe) THEN
    CALL Rcf_U_to_P( u_level, u_at_p, grid )
    CALL Rcf_V_to_P( v_level, v_at_p, grid )

    IF ( Mode == ToStandard ) THEN
      CALL vector_eq_to_latlon( Coeff3, Coeff4, u_at_p, v_at_p,  &
                                u_at_p_out, v_at_p_out, p_size, l_mdi )
    ELSE IF ( Mode == FromStandard ) THEN
      CALL vector_latlon_to_eq( Coeff1, Coeff2, u_at_p, v_at_p,  &
                                u_at_p_out, v_at_p_out, p_size, l_mdi )
    END IF

    CALL Rcf_P_to_U( u_at_p_out, u_level, grid )
    CALL Rcf_P_to_V( v_at_p_out, v_level, grid )
  END IF

  pe = 0
  DO i = 1, rem
    j = nproc * div + i

    CALL Rcf_Scatter_Field_Real( u_field % DATA(:,j), u_level,         &
                                 u_field % row_len,                    &
                                 u_field % rows,                       &
                                 u_field % glob_row_len,               &
                                 u_field % glob_rows,pe,               &
                                 gc_all_proc_group )

    CALL Rcf_Scatter_Field_Real( v_field % DATA(:,j), v_level,         &
                                 v_field % row_len,                    &
                                 v_field % rows,                       &
                                 v_field % glob_row_len,               &
                                 v_field % glob_rows,pe,               &
                                 gc_all_proc_group )

    pe = pe + 1

    IF (pe == nproc) pe = 0
  END DO


  DEALLOCATE( u_level );
  DEALLOCATE( v_level );
  DEALLOCATE( u_at_p  );
  DEALLOCATE( v_at_p  );
  DEALLOCATE( u_at_p_out );
  DEALLOCATE( v_at_p_out );

  !------------------------------------------------------------------
  ! Write the rotated data back to file and clear up
  !------------------------------------------------------------------
  CALL Rcf_Write_Field( u_field, hdr, decomp )
  CALL Rcf_Write_Field( v_field, hdr, decomp )

  CALL Rcf_DeAlloc_Field( u_field )
  CALL Rcf_DeAlloc_Field( v_field )

END DO

!-----------------------------------------------------------------
! Back to the original decomposition
!-----------------------------------------------------------------
IF ( decomp_original /= current_decomp_type ) THEN
  CALL Change_Decomposition( decomp_original )
END IF

IF (LTimer) CALL Timer( 'Rcf_Rotate', 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Rotate
END MODULE Rcf_Rotate_Mod
