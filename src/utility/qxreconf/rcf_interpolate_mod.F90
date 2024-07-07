! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  top level wrapper for interpolation

MODULE Rcf_interpolate_mod
IMPLICIT NONE

!  Subroutine Rcf_Interpolate - top level interpolation wrapper
!
! Description:
! This module contains a top-level wrapper subroutine for
! interpolation (both horizontal and vertical). It handles
! conversion between data-types (log,int and real) .
! It is worth noting that the orography fields that are fed in
! need not exist if only horizontal interpolation is done.
!
! Method:
!  Intermediate temporary fields are set up - heights are generated,
!  horizontal and then vertical interpolation is done - and then
!  data is reconverted as required.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INTERPOLATE_MOD'

CONTAINS

SUBROUTINE Rcf_interpolate( field_in, field_out, grid_in, grid_out, &
                            interp_orography_in, orography_out )

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_type

USE Rcf_horizontal_mod, ONLY: &
    Rcf_horizontal

USE Rcf_vertical_mod, ONLY: &
    Rcf_vertical

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,             &
    Rcf_Dealloc_Field

USE Rcf_generate_heights_mod, ONLY: &
    Rcf_generate_heights

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Vert_Cloud_Mod, ONLY: &
    Rcf_Vert_Cloud

USE Rcf_Set_Interp_Flags_Mod, ONLY:                         &
    interp_v_only,          interp_h_only,                   &
    interp_all,             interp_done,                     &
    interp_copied,          interp_copy

USE um_stashcode_mod, ONLY: &
    stashcode_cct,             &
    stashcode_ccb,             &
    stashcode_lvl_bse_dp_sc,   &
    stashcode_lvl_top_dp_sc,   &
    stashcode_prog_sec

USE UM_Parvars, ONLY:         &
    change_decomposition

USE stparam_mod, ONLY:  &
    st_levels_single,   &
    st_levels_model_rho

USE cppxref_mod, ONLY:             &
    ppx_type_int,                   &
    ppx_type_log

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (field_type), INTENT(INOUT) :: field_in
TYPE (field_type), INTENT(INOUT) :: field_out
TYPE (grid_type), INTENT(IN)     :: grid_in
TYPE (grid_type), INTENT(IN)     :: grid_out

! This argument (pos dummy) is for the *FINAL OUTPUT* orography
TYPE (field_type), INTENT(IN)    :: orography_out

! This argument (pos dummy) is for the *INTERPOLATED INPUT* orography
TYPE (field_type), INTENT(IN)    :: interp_orography_in

! Local variables/parameters
TYPE (grid_type)                 :: grid_middle
TYPE (field_type)                :: field_middle
CHARACTER (LEN=*), PARAMETER     :: RoutineName = 'RCF_INTERPOLATE'
CHARACTER (LEN=errormessagelength)  :: Cmessage
INTEGER                          :: ErrorStatus
INTEGER                          :: i
INTEGER                          :: j
INTEGER                          :: log_to_int
REAL, ALLOCATABLE                :: heights_middle( :, : )
REAL, ALLOCATABLE                :: heights_out( :, : )
LOGICAL                          :: vertical_required
LOGICAL                          :: horizontal_required
LOGICAL                          :: interp_required

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

vertical_required   = (field_in % interp == interp_v_only .OR. &
                       field_in % interp == interp_all )
horizontal_required = (field_in % interp == interp_h_only .OR. &
                       field_in % interp == interp_all )
interp_required     = (vertical_required .OR. horizontal_required)

!------------------------------------------------------------
! Tests to see if orography is set should it be required
!------------------------------------------------------------
IF ( vertical_required ) THEN
  IF ( .NOT. ALLOCATED( interp_orography_in % data ) .OR. &
       .NOT. ALLOCATED( orography_out % data ) ) THEN
    Cmessage = 'Orography data not present when vertical '//&
               'interpolation required'
    ErrorStatus = 10
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

!------------------------------------------------------------
! Test to see if Data is allocated where necessary
!------------------------------------------------------------
IF ( (.NOT. ALLOCATED( field_in % data  )       .AND. &
      .NOT. ALLOCATED( field_in % Data_Int  )   .AND. &
      .NOT. ALLOCATED( field_in % Data_Log  ) ) .OR.  &
     (.NOT. ALLOCATED( field_out % data )       .AND. &
      .NOT. ALLOCATED( field_out % Data_Int )   .AND. &
      .NOT. ALLOCATED( field_out % Data_Log ) ) ) THEN
  Cmessage = 'An interpolation field has not had Data space allocated!'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF


!------------------------------------------------------------------
! Do the interpolation in the relevant way.
!------------------------------------------------------------------
! Force convective cloud base and convective cloud top to have
! lv_code = st_levels_model_rho for the duration of this routine to
! generate correct height field
! Note that have to check section as well as item number
! THIS IS DANGEROUS AS IT ALTERS THE INTERNAL STASHMASTER!!!!
IF ( (field_out % stashmaster % item    == stashcode_cct  .OR.  &
     field_out % stashmaster % item    == stashcode_ccb  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_bse_dp_sc  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_top_dp_sc) .AND. &
     field_out % stashmaster % section == stashcode_prog_sec ) THEN
  field_in  % stashmaster % lv_code = st_levels_model_rho
  field_out % stashmaster % lv_code = st_levels_model_rho
END IF

! middle grid/field have out horizont resolution and in vert.
grid_middle = grid_out     ! horizontal resolutions correct
grid_middle % model_levels = grid_in % model_levels
grid_middle % cloud_levels = grid_in % cloud_levels
grid_middle % st_levels    = grid_in % st_levels
grid_middle % sm_levels    = grid_in % sm_levels
grid_middle % bl_levels    = grid_in % bl_levels
grid_middle % ozone_levels = grid_in % ozone_levels
grid_middle % tr_levels    = grid_in % tr_levels
grid_middle % z_top_of_model = grid_in % z_top_of_model
grid_middle % first_constant_r_rho_level =                          &
                               grid_in % first_constant_r_rho_level
grid_middle % height_gen_method = grid_in % height_gen_method

grid_middle % eta_theta_levels => grid_in % eta_theta_levels
grid_middle % eta_rho_levels   => grid_in % eta_rho_levels
grid_middle % soil_depths      => grid_in % soil_depths

field_middle % levels          = field_in  % levels
field_middle % bottom_level    = field_in  % bottom_level
field_middle % top_level       = field_in  % top_level
field_middle % interp          = field_in  % interp
field_middle % rows            = field_out % rows
field_middle % row_len         = field_out % row_len
field_middle % level_size      = field_out % level_size
field_middle % glob_rows       = field_out % glob_rows
field_middle % glob_row_len    = field_out % glob_row_len
field_middle % glob_level_size = field_out % glob_level_size
field_middle % stashmaster     => field_in % stashmaster

! Allocate field_middle Data space
CALL Rcf_Alloc_Field( field_middle )

! Allocate heights - need for vertical function call
ALLOCATE( heights_middle( field_middle % level_size,            &
                      0 :  grid_middle % model_levels + 1) )


ALLOCATE( heights_out( field_out % level_size,                  &
                        0 :  grid_out % model_levels + 1) )


! Only need to generate heights if vertical interpolation is actually
! done.  Check at start of routine is done so we can safely use orography data.
IF (vertical_required) THEN

  CALL rcf_generate_heights( grid_middle, interp_orography_in,        &
                         field_middle % stashmaster % grid_type,      &
                         field_middle % stashmaster % lv_code,        &
                         heights_middle,                              &
                         field_middle % level_size )

  CALL rcf_generate_heights( grid_out, orography_out,             &
                         field_out % stashmaster % grid_type,     &
                         field_out % stashmaster % lv_code,       &
                         heights_out,                             &
                         field_out % level_size )
END IF

!------------------------------------------------------------------
! If data is integer or logical, we will need to convert it to
! Real to do any interpolation on it. Only need to do the conversion
! if some interpolation is actually required!
!------------------------------------------------------------------
IF ( interp_required ) THEN

  IF (field_in % stashmaster % data_type == ppx_type_int ) THEN
    ALLOCATE( field_in % DATA( field_in % level_size,       &
                               field_in % levels ) )
    ALLOCATE( field_out % DATA( field_out % level_size,     &
                                field_out % levels ) )
    ALLOCATE( field_middle % DATA( field_middle % level_size,  &
                                   field_middle % levels ) )

    field_in % DATA(:,:) = REAL( field_in % Data_Int(:,:) )

  ELSE IF (field_in % stashmaster % data_type == ppx_type_log ) THEN

    ALLOCATE( field_in % DATA( field_in % level_size,       &
                               field_in % levels ) )
    ALLOCATE( field_out % DATA( field_out % level_size,     &
                                field_out % levels ) )
    ALLOCATE( field_middle % DATA( field_middle % level_size, &
                                   field_middle % levels ) )
    DO j = 1, field_in % levels
      DO i = 1, field_in % level_size
        IF ( field_in % Data_Log(i,j) ) THEN
          field_in % DATA(i,j) = 1.0
        ELSE
          field_in % DATA(i,j) = 0.0
        END IF
      END DO
    END DO
  END IF
END IF

!-------------------------------------------------------------------
! Now can interpolate
!-------------------------------------------------------------------
CALL Rcf_horizontal( field_in, field_middle, grid_in, grid_middle )

! Ensure we are on the correct (ie output) decomposition before
! we do vertical interpolation
CALL Change_Decomposition( decomp_rcf_output )

IF ( (field_out % stashmaster % item    == stashcode_cct  .OR.  &
     field_out % stashmaster % item    == stashcode_ccb  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_bse_dp_sc  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_top_dp_sc) .AND. &
     field_out % stashmaster % section == stashcode_prog_sec ) THEN

  CALL Rcf_Vert_Cloud( field_middle, field_out, grid_middle, grid_out, &
                       heights_middle, heights_out)

  ! convert lv_code in stashmaster back to the correct value
  field_in  % stashmaster % lv_code = st_levels_single
  field_out % stashmaster % lv_code = st_levels_single
ELSE
  CALL Rcf_vertical( field_middle, field_out, grid_middle, grid_out,&
                     heights_middle, heights_out)
END IF

! Deallocate heights
DEALLOCATE( heights_middle )
DEALLOCATE( heights_out )

!-------------------------------------------------------------------
! Deallocate field_middle Data
!-------------------------------------------------------------------
CALL Rcf_Dealloc_Field( field_middle )

!-------------------------------------------------------------------
! Need to convert back to integer or logical if field is such. Only
! need to do this if interpolation active!
!-------------------------------------------------------------------
IF ( interp_required ) THEN
  IF (field_out % stashmaster % data_type == ppx_type_int ) THEN
    field_out % Data_Int(:,:) = NINT( field_out % DATA(:,:) )

    DEALLOCATE( field_out % DATA )
    DEALLOCATE( field_in % DATA )
    ! middle data already deallocated by dealloc_field above

  ELSE IF (field_out % stashmaster % data_type == ppx_type_log ) THEN
    DO i = 1, field_out % level_size
      DO j = 1, field_out % levels
        log_to_int = INT( field_out % DATA(i,j) + 0.5 )
        IF ( log_to_int == 1 ) THEN
          field_out % Data_Log(i,j) =  .TRUE.
        ELSE IF (log_to_int == 0 ) THEN
          field_out % Data_Log(i,j) =  .FALSE.
        ELSE
          Cmessage = 'Conversion from Logical to Integer: Illegal value'
          ErrorStatus = 30
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF
      END DO
    END DO

    DEALLOCATE( field_out % DATA )
    DEALLOCATE( field_in % DATA )
  END IF
END IF

!--------------------------------------------------------------------
! Set the interp flag for the output field
!--------------------------------------------------------------------
IF (interp_required) THEN
  field_out % interp = interp_done

ELSE IF (field_in % interp == interp_copy) THEN
  field_out % interp = interp_copied

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_interpolate

END MODULE Rcf_interpolate_mod
