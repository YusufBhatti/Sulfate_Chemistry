! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Horizontally interpolates T and Pstar for Height generation

MODULE Rcf_GRIB_Interp_TnPstar_Mod
IMPLICIT NONE

!  Subroutine Rcf_GRIB_Interp_TnPstar
!
! Description: A routine to allocate space and call rcf_interpolate
!              to provide horizontally interpolated versions of
!              input T and Pstar on the output grid's horizontal
!              resolution. (copied and modified rcf_set_orography)
!
! Method: Allocate appropriate space
!         Set interp flag to horizontal only
!         Call rcf_interpolate
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_INTERP_TNPSTAR_MOD'

CONTAINS

SUBROUTINE Rcf_GRIB_Interp_TnPstar( fields_in, fields_out, Hdr_In,    &
                                    field_count_in, field_count_out)

USE science_fixes_mod, ONLY: l_fix_ec_gen_hgt

USE planet_constants_mod, ONLY: c_virtual

USE umPrintMgr, ONLY:      &
    umPrint,                &
    PrintStatus,                &
    PrStatus_Normal

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,               &
    Output_Grid,              &
    grid_type

USE decomp_params, ONLY: &
    decomp_rcf_input

USE Rcf_Set_Interp_Flags_Mod, ONLY:  &
    interp_copy,                      &
    interp_h_only

USE um_stashcode_mod, ONLY: &
    stashcode_orog,          &
    stashcode_pstar,         &
    stashcode_theta,         &
    stashcode_t,             &
    stashcode_q,             &
    stashcode_qcl,           &
    stashcode_qcf,           &
    stashcode_qrain,         &
    stashcode_qgraup,        &
    stashcode_prog_sec

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_field_equals_mod, ONLY: &
    Rcf_field_equals

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE UM_ParCore, ONLY: &
    mype

USE Rcf_GRIB_T_n_Pstar_H_Interp_Mod, ONLY: &
    grib_tv,           &
    GRIB_Pstar,        &
    GRIB_Levels

USE Rcf_Headaddress_Mod, ONLY: &
    ldc_mlindex

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE Ereport_Mod, ONLY: &
    Ereport

USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, only: newline

IMPLICIT NONE

TYPE( field_type ), POINTER          :: Pstar_in

! Arguments
TYPE( field_type ), POINTER          :: fields_in(:)
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
INTEGER, INTENT(IN)                  :: field_count_in
INTEGER, INTENT(IN)                  :: field_count_out

! Local variables

TYPE( grid_type )                    :: grid_middle
TYPE( field_type ), POINTER          :: orog_out
TYPE( field_type ), POINTER          :: T_in
TYPE( field_type ), POINTER          :: T_out
TYPE( field_type ), POINTER          :: qx_in
INTEGER                              :: pos,icode
INTEGER                              :: i,k
REAL,ALLOCATABLE                     :: mfac(:,:)
TYPE (field_type)                    :: dummy

CHARACTER (LEN=errormessagelength) :: CMessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_GRIB_INTERP_TNPSTAR'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  IF (l_fix_ec_gen_hgt) THEN
    CALL umPrint( 'Horizontally Interpolating T, qX and Pstar', &
        src='rcf_grib_interp_tnpstar_mod')
  ELSE
    CALL umPrint( 'Horizontally Interpolating T and Pstar',     &
        src='rcf_grib_interp_tnpstar_mod')
  END IF
END IF

! Setup grid_middle which has the vertical information from input grid but
! horizontal information from output grid.

grid_middle                              = output_grid
grid_middle % model_levels               = input_grid % model_levels
grid_middle % cloud_levels               = input_grid % cloud_levels
grid_middle % st_levels                  = input_grid % st_levels
grid_middle % sm_levels                  = input_grid % sm_levels
grid_middle % bl_levels                  = input_grid % bl_levels
grid_middle % ozone_levels               = input_grid % ozone_levels
grid_middle % tr_levels                  = input_grid % tr_levels
grid_middle % conv_levels                = input_grid % conv_levels
grid_middle % height_gen_method          = input_grid % height_gen_method
grid_middle % first_constant_r_rho_level = &
  input_grid % first_constant_r_rho_level
grid_middle % z_top_of_model             = input_grid % z_top_of_model
grid_middle % eta_theta_levels          => input_grid % eta_theta_levels
grid_middle % eta_rho_levels            => input_grid % eta_rho_levels
grid_middle % rhcrit                    => input_grid % rhcrit
grid_middle % soil_depths               => input_grid % soil_depths

!------------------------------------------------------------------
! find and read input Pstar
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_pstar,                &
                 fields_in, field_count_in, pos )
Pstar_in  => fields_in(pos)

CALL Rcf_Alloc_Field( Pstar_in )
CALL Rcf_Read_Field( Pstar_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Setup Pstar field with output horizontal resolution
!------------------------------------------------------------------
! Pstar is (probably) not in the output dump, it does however need
! to be interpolated. Use the output orography for the 'donor' level
! descriptor.
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos )
orog_out => fields_out(pos)
CALL Rcf_Alloc_Field( orog_out )

! Set the sizes of GRIB_Pstar to be those of orog_out
CALL rcf_field_equals(GRIB_Pstar, orog_out)

! Now need to reset some values which aren't correct
GRIB_Pstar % dump_pos        = Pstar_in % dump_pos
GRIB_Pstar % interp          = interp_h_only
GRIB_Pstar % stashmaster    => Pstar_in % stashmaster

!------------------------------------------------------------------
! Interpolate Pstar
!------------------------------------------------------------------
IF (h_int_active) THEN
  Pstar_in % interp = interp_h_only
ELSE
  Pstar_in % interp = interp_copy
END IF


CALL Rcf_Alloc_Field( GRIB_Pstar )
CALL rcf_interpolate( Pstar_in, GRIB_Pstar, Input_Grid,       &
                      grid_middle, dummy, dummy )

!------------------------------------------------------------------
! find and read input T
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_t/1000, MOD(stashcode_t,1000),              &
                 fields_in, field_count_in, pos )
T_in  => fields_in(pos)
CALL Rcf_Alloc_Field( T_in )
CALL Rcf_Read_Field( T_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Apply moist species fix
!------------------------------------------------------------------
! If using height generation fix, the temperature that goes into
! the height generation algorithm is actually virtual temperature. This
! makes sure that the mass of moist species is accounted for. Note that
! all moist species should really be present (q, qcl, qcf, qrain and
! qsnow), however only q is insisted upon. If using qrain and qsnow
! then note that qrain and qgraup (which the GRIB reconfiguration
! hijacks to make this work) will be explicitly set to zero afterwards
! via rcf_set_data_source as these are not known to be the same
! quantities.
IF (l_fix_ec_gen_hgt) THEN

  ! Set up array to hold moist factor to calculate virtual temperature
  ALLOCATE(mfac(T_in % level_size, T_in % levels))
  DO k=1, T_in % levels
    DO i=1, T_in % level_size
      mfac(i,k) = 1.0
    END DO
  END DO

  ! Set up standard warning message and code
  WRITE(Cmessage,FMT='(A)') ' not in input GRIB file.'//newline//   &
                  'This will compromise the height level calculation'

  !------------------------------------------------------------------
  ! find and read input q (compulsory)
  !------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_q/1000, MOD(stashcode_q,1000),         &
                    fields_in, field_count_in, pos, zero_ok_arg=.TRUE. )
  IF (pos /= 0) THEN
    qx_in  => fields_in(pos)
    CALL Rcf_Alloc_Field( qx_in )
    CALL Rcf_Read_Field( qx_in, Hdr_In, decomp_rcf_input )
    DO k=1, T_in % levels
      DO i=1, T_in % level_size
        mfac(i,k) = mfac(i,k) + c_virtual * qx_in % DATA(i,k)
      END DO
    END DO
    CALL Rcf_DeAlloc_Field( qx_in )
  ELSE
    icode = 10   ! Definitely fail here as impact is huge!
    CALL Ereport( RoutineName, icode, 'Q'//Cmessage )
  END IF

  !------------------------------------------------------------------
  ! find and read input qcl (optional)
  !------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_qcl/1000, MOD(stashcode_qcl,1000),     &
                    fields_in, field_count_in, pos, zero_ok_arg=.TRUE. )
  IF (pos /= 0) THEN
    qx_in  => fields_in(pos)
    CALL Rcf_Alloc_Field( qx_in )
    CALL Rcf_Read_Field( qx_in, Hdr_In, decomp_rcf_input )
    DO k=1, T_in % levels
      DO i=1, T_in % level_size
        mfac(i,k) = mfac(i,k) - qx_in % DATA(i,k)
      END DO
    END DO
    CALL Rcf_DeAlloc_Field( qx_in )
  ELSE
    icode = -20   ! Only need to warn here due to small impact
    CALL Ereport( RoutineName, icode, 'QCL'//Cmessage )
  END IF

  !------------------------------------------------------------------
  ! find and read input qcf (optional)
  !------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_qcf/1000, MOD(stashcode_qcf,1000),     &
                    fields_in, field_count_in, pos, zero_ok_arg=.TRUE. )
  IF (pos /= 0) THEN
    qx_in  => fields_in(pos)
    CALL Rcf_Alloc_Field( qx_in )
    CALL Rcf_Read_Field( qx_in, Hdr_In, decomp_rcf_input )
    DO k=1, T_in % levels
      DO i=1, T_in % level_size
        mfac(i,k) = mfac(i,k) - qx_in % DATA(i,k)
      END DO
    END DO
    CALL Rcf_DeAlloc_Field( qx_in )
  ELSE
    icode = -30   ! Only need to warn here due to small impact
    CALL Ereport( RoutineName, icode, 'QCF'//Cmessage )
  END IF

  !------------------------------------------------------------------
  ! find and read input qrain (optional)
  !------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_qrain/1000, MOD(stashcode_qrain,1000), &
                    fields_in, field_count_in, pos, zero_ok_arg=.TRUE. )
  IF (pos /= 0) THEN
    qx_in  => fields_in(pos)
    CALL Rcf_Alloc_Field( qx_in )
    CALL Rcf_Read_Field( qx_in, Hdr_In, decomp_rcf_input )
    DO k=1, T_in % levels
      DO i=1, T_in % level_size
        mfac(i,k) = mfac(i,k) - qx_in % DATA(i,k)
      END DO
    END DO
    CALL Rcf_DeAlloc_Field( qx_in )
  ELSE
    icode = -40   ! Only need to warn here due to small impact
    CALL Ereport( RoutineName, icode, 'QRAIN'//Cmessage )
  END IF

  !------------------------------------------------------------------
  ! find and read input qsnow/qgraup (optional)
  !------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_qgraup/1000, MOD(stashcode_qgraup,1000), &
                    fields_in, field_count_in, pos, zero_ok_arg=.TRUE. )
  IF (pos /= 0) THEN
    qx_in  => fields_in(pos)
    CALL Rcf_Alloc_Field( qx_in )
    CALL Rcf_Read_Field( qx_in, Hdr_In, decomp_rcf_input )
    DO k=1, T_in % levels
      DO i=1, T_in % level_size
        mfac(i,k) = mfac(i,k) - qx_in % DATA(i,k)
      END DO
    END DO
    CALL Rcf_DeAlloc_Field( qx_in )
  ELSE
    icode = -50   ! Only need to warn here due to small impact
    CALL Ereport( RoutineName, icode, 'QSNOW'//Cmessage )
  END IF

  !------------------------------------------------------------------
  ! Calculate virtual temperature
  !------------------------------------------------------------------
  DO k=1, T_in % levels
    DO i=1, T_in % level_size
      T_in % DATA(i,k) = T_in % DATA(i,k) * mfac(i,k)
    END DO
  END DO

  DEALLOCATE(mfac)
END IF

!------------------------------------------------------------------
! Setup T field with output horizontal resolution
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos )
T_out => fields_out(pos)
CALL Rcf_Alloc_Field( T_out )

! Set the sizes of GRIB_Tv to be those of T_out
CALL rcf_field_equals(grib_tv, T_out)

! Now need to reset some values to maintain same vertical levels
grib_tv % levels          = T_in % levels
grib_tv % bottom_level    = T_in % bottom_level
grib_tv % top_level       = T_in % top_level
! This is probably not required but consistent with what is happening here.
! To support Endgame we used to have 2 different STASHmasters but this has
! been retired so we only need one which means this is redundant but correct.
grib_tv % stashmaster     => T_in % stashmaster

!------------------------------------------------------------------
! Interpolate T
!------------------------------------------------------------------
IF (h_int_active) THEN
  T_in % interp = interp_h_only
ELSE
  T_in % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( grib_tv )
CALL rcf_interpolate( T_in, grib_tv, Input_Grid,       &
                      grid_middle, dummy, dummy )

!------------------------------------------------------------------
! Save level definitons
!------------------------------------------------------------------
ALLOCATE ( GRIB_Levels(T_in % levels) )
GRIB_Levels(:) = Hdr_In % LevDepC (1:T_in % levels,ldc_mlindex)

!------------------------------------------------------------------
! Clean up
!------------------------------------------------------------------
CALL Rcf_DeAlloc_Field( Pstar_in )
CALL Rcf_DeAlloc_Field( T_in )
CALL Rcf_DeAlloc_Field( T_out )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_GRIB_Interp_TnPstar
END MODULE Rcf_GRIB_Interp_TnPstar_Mod
