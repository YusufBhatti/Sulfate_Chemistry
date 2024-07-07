! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Dump creation `main loop'

MODULE Rcf_Create_Dump_Mod

IMPLICIT NONE

!  Subroutine Rcf_Create_Dump - main dump creation loop.
!
! Description:
!   Creates the output dump based on the input choices.
!
! Method:
!   Data sources are setup and then fields in the output dump are
!   looped over. Each one is set appropriately ( interpolated from
!   input dump, set to 0/MDI/constant/from file/etc). Pre and Post
!   processing is performed as required.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CREATE_DUMP_MOD'

CONTAINS

SUBROUTINE Rcf_Create_Dump ( Hdr_In, Hdr_Out, fields_in, fields_out,          &
                             field_count_in, field_count_out, data_source )

USE Rcf_Set_Orography_Mod, ONLY: &
    Rcf_Set_Orography

USE rcf_set_data_source_mod, ONLY: rcf_reset_data_source

USE Rcf_Pre_Interp_Transform_Mod, ONLY: &
    Rcf_Pre_Interp_Transform

USE Rcf_Post_Interp_Transform_Mod, ONLY: &
    Rcf_Post_Interp_Transform

USE Rcf_Post_Process_Mod, ONLY:&
    Rcf_Post_Process_Atmos

USE rcf_nlist_recon_technical_mod, ONLY: &
    l_basic_interp,                      &
    l_trans,                             &
    transp

USE file_manager, ONLY: &
    assign_file_unit,    &
    release_file_unit

USE filenamelength_mod, ONLY: &
    filenamelength

USE Rcf_FreeUMhdr_mod, ONLY: &
    Rcf_FreeUMhdr

USE Rcf_UMhead_Mod, ONLY: &
    UM_header_type

USE Ereport_mod, ONLY: &
    Ereport

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type,          &
    other_field,               &
    already_processed

USE items_nml_mod, ONLY:      &
    input_dump,               &
    ancillary_file,           &
    set_to_zero,              &
    set_to_mdi,               &
    tracer_file,              &
    set_to_const,             &
    external_dump,            &
    field_calcs,              &
    field_dependent_calcs,    &
    netcdf_file

USE rcf_read_field_mod, ONLY: &
    Rcf_Read_Field

USE rcf_write_field_mod, ONLY: &
    Rcf_Write_Field

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE decomp_params, ONLY: &
    Decomp_rcf_input,    &
    Decomp_rcf_output

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,           &
    Output_Grid

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE UM_ParCore, ONLY: &
    mype

USE Rcf_field_equals_mod, ONLY: &
    Rcf_field_equals

USE Rcf_Interp_Weights_Mod,        ONLY: &
    h_int_method,                        &
    nearest_neighbour,                   &
    h_int_active_u,                      &
    h_int_active_v,                      &
    l_same_rotation

USE rcf_interpolate_mod, ONLY: &
    rcf_interpolate

USE Rcf_aux_file_mod, ONLY:     &
    tracers,       user_prog,    &
    transplant,    rcf_aux_file

USE Rcf_Field_Calcs_Mod, ONLY: &
    Rcf_Field_Calcs

USE Rcf_Field_Dependent_Calcs_Mod, ONLY: &
    Rcf_Field_Dependent_Calcs

USE Rcf_Pre_Process_Calcs_Mod, ONLY: &
    Rcf_Pre_Process_Calcs

USE Rcf_Set_Interp_Flags_Mod, ONLY:               &
    Rcf_Set_Interp_Flags,                          &
    interp_v_only,            interp_h_only,       &
    interp_all,               interp_no_op,        &
    interp_copy,              interp_done

USE Rcf_Rotate_Mod, ONLY: &
    Rcf_Rotate,            &
    ToStandard,            &
    FromStandard

USE um_stashcode_mod, ONLY: &
    stashcode_u,               &
    stashcode_v,               &
    stashcode_lsm,             &
    stashcode_orog,            &
    stashcode_icefrac,         &
    stashcode_tstar,           &
    stashcode_tstar_land,      &
    stashcode_tstar_sea,       &
    stashcode_tstar_sice,      &
    stashcode_prog_sec

USE Rcf_GRIB_Interp_TnPstar_Mod, ONLY: &
    Rcf_GRIB_Interp_TnPstar

USE Rcf_Generate_Heights_Mod, ONLY: &
    height_gen_ecmwf_hybrd,          &
    height_gen_ecmwf_press

USE Rcf_locate_alt_field_mod, ONLY: &
    rcf_locate_alt_field

USE cppxref_mod, ONLY: &
    ppx_type_real,     &
    ppx_type_int,      &
    ppx_type_log

USE rcf_polar_wind_mod, ONLY: &
    rcf_polar_wind

USE rcf_headaddress_mod, ONLY: &
    RC_LongSpacing,             &
    FH_GridStagger,             &
    FH_GridStagger_Endgame

USE rcf_ideal_initialisation_mod, ONLY: &
    rcf_ideal_initialisation

USE rcf_nlist_recon_idealised_mod, ONLY: &
    l_init_idealised

USE io_constants, ONLY: ioNoDelete, ioOpenReadOnly
USE model_file, ONLY: model_file_open, model_file_close
USE missing_data_mod, ONLY:  rmdi, imdi

USE get_env_var_mod, ONLY: get_env_var

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)                  :: field_count_in
INTEGER, INTENT(IN)                  :: field_count_out
TYPE (UM_Header_Type), INTENT(INOUT) :: Hdr_In
TYPE (UM_Header_Type), INTENT(INOUT) :: Hdr_Out
TYPE (field_type), POINTER           :: fields_in(:)
TYPE (field_type), POINTER           :: fields_out(:)
TYPE (data_source_type), POINTER     :: data_source(:)


! Local vars
INTEGER                            :: pos,pos1  ! position in fields array
INTEGER                            :: i         ! looper
INTEGER                            :: err       ! for file open/close
INTEGER                            :: IOStatus
INTEGER                            :: ErrorStatus
INTEGER                            :: stashcode_onpole
INTEGER                            :: stashcode_offpole
LOGICAL                            :: l_exist
CHARACTER (LEN=errormessagelength) :: Cmessage
CHARACTER (LEN=filenamelength)     :: TracerFileName
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_CREATE_DUMP'

TYPE (field_type), TARGET  :: interp_orog    ! interpolated input orog
TYPE (field_type), POINTER :: orog_out       ! ptr to output orog
TYPE (field_type), POINTER :: orog_in        ! ptr to intput orog
TYPE (field_type), POINTER :: input_field    ! ptr to input field
TYPE (field_type), POINTER :: output_field   ! ptr to output field
TYPE (field_type), POINTER :: polar_wind     ! ptr to polar wind
TYPE (field_type), POINTER :: non_polar_wind ! ptr to non-polar wind
TYPE (UM_Header_Type)      :: Hdr_Aux        ! auxillary file header

! Formatting
CHARACTER (LEN=*), PARAMETER :: &
    FORM="(a20,1x,i3, ' ( Section', i4, ' )', "//&
    "' ( Stashcode', i4, ' )', ' ', a36)"
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------
! Set up the input, output and interpolated orography fields
!-----------------------------------------------------------
CALL Rcf_Set_Orography( fields_in, fields_out, field_count_in,                &
                        field_count_out, hdr_in, hdr_out,                     &
                        data_source, orog_in, orog_out, interp_orog )

! Reset data_source for output fields which either do not appear
! in input dump or require special processing
CALL rcf_reset_data_source( data_source, fields_in, fields_out,               &
                            field_count_in, field_count_out,                  &
                            hdr_in, hdr_out )

! If recon from GRIB T & Pstar need to interp'd horizontally
! for height generation
IF ( ( input_grid % height_gen_method == height_gen_ecmwf_hybrd )   .OR.      &
     ( input_grid % height_gen_method == height_gen_ecmwf_press ) ) THEN
  CALL Rcf_GRIB_Interp_TnPstar( fields_in, fields_out, Hdr_In,                &
                                field_count_in, field_count_out)
END IF

!-------------------------------------------------------------
! Setup interpolation flags
!-------------------------------------------------------------
CALL Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in,             &
                           field_count_out, data_source )

!------------------------------------------------------------
! Rotate input winds if so required
!------------------------------------------------------------
IF ( (h_int_active_u .OR. h_int_active_v) .AND. Input_Grid % Rotated .AND.    &
     .NOT. l_same_rotation ) THEN
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint( 'Rotating input winds', src='rcf_create_dump_mod')
  END IF

  CALL Rcf_Rotate( fields_in, field_count_in, Input_Grid, Hdr_In,             &
                   decomp_rcf_input, ToStandard)
END IF

!-----------------------------------------------------
! Main loop
!-----------------------------------------------------
DO i = 1, field_count_out ! run through the output fields

  !--------------------------------------------
  ! Which fields have we already dealt with? For example.
  ! Orography
  ! LSM
  ! Fields in rcf_ancil_atmos which have be treated specially. (e.g. tstar)
  !--------------------------------------------
  IF ( data_source( i ) % source == Already_Processed ) THEN

    IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Already processed', i,                           &
                      fields_out( i ) % stashmaster % section,                &
                      fields_out( i ) % stashmaster % item,                   &
                      fields_out( i ) % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    CYCLE    ! just skip this iteration of the loop.

  END IF

  ! To simplify the specification of fields to process we need to make sure we
  ! dont revisit ancillary files since the orography will be rewritten if
  ! orograhy is read in via ancillary.
  output_field => fields_out( i )
  IF ( ( Data_Source( i ) % Source /= Ancillary_File ) .AND.                  &
       ( Data_Source( i ) % Source /= NetCDF_File) )   THEN
    CALL Rcf_Alloc_Field( output_field )
  END IF

  SELECT CASE( Data_Source( i ) % Source )

    !---------------------------------------------------------------
    ! Data interpolated (or copied) from input dump
    !---------------------------------------------------------------
  CASE ( Input_Dump, Other_field )

    ! set up input fields
    CALL Rcf_Locate( fields_out( i ) % stashmaster % section,                 &
                     fields_out( i ) % stashmaster % item,                    &
                     fields_in, field_count_in, pos, zero_ok_arg = .TRUE. )

    IF (pos == 0 .AND. Data_Source( i ) % Source == Other_field ) THEN
      ! Lets try finding another field we can use.
      CALL Rcf_locate_alt_field(fields_out( i ), fields_in,                   &
                                field_count_in, pos)
    END IF
    ! Read in field
    input_field => fields_in( pos )
    CALL Rcf_Alloc_Field( input_field )

    CALL Rcf_Read_Field( input_field, Hdr_In, decomp_rcf_input )

    ! Write out appropriate message
    SELECT CASE( input_field % interp )

    CASE ( interp_h_only, interp_v_only, interp_all)
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,FORM) 'Interpolating Field', i,                       &
                       input_field % stashmaster % section,                   &
                       input_field % stashmaster % item,                      &
                       input_field % stashmaster % NAME
        CALL umPrint(umMessage,src='rcf_create_dump_mod')
      END IF

    CASE ( interp_copy )
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,FORM) 'Copying Field', i,                             &
                       input_field % stashmaster % section,                   &
                       input_field % stashmaster % item,                      &
                       input_field % stashmaster % NAME
        CALL umPrint(umMessage,src='rcf_create_dump_mod')
      END IF

    CASE ( interp_no_op )
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,FORM) 'Skipping Field', i,                            &
                       input_field % stashmaster % section,                   &
                       input_field % stashmaster % item,                      &
                       input_field % stashmaster % NAME
        CALL umPrint(umMessage,src='rcf_create_dump_mod')
      END IF

      CALL Rcf_DeAlloc_Field( input_field )
      CALL Rcf_DeAlloc_Field( output_field )
      CYCLE

    END SELECT

    ! If using nearest neighbour we dont really care about transforming data.
    IF (h_int_method /= nearest_neighbour  .AND. &
        .NOT. l_basic_interp ) THEN
      ! convert fields approriately for interpolation
      CALL Rcf_Pre_Interp_Transform( input_field, fields_in,                  &
                                     field_count_in, hdr_in,                  &
                                     fields_out, field_count_out,             &
                                     hdr_out, data_source,                    &
                                     orog_in )
    END IF

    CALL Rcf_Interpolate( input_field, output_field, Input_Grid,              &
                          Output_Grid, interp_orog, orog_out )

    ! If using nearest neighbour we dont really care about transforming data.
    IF (h_int_method /= nearest_neighbour .AND. &
        .NOT. l_basic_interp ) THEN
      ! Convert fields back to original form (if possible) and
      ! perform any simple post-processing required
      CALL Rcf_Post_Interp_Transform( output_field, fields_out,               &
                                      field_count_out )
    END IF

    CALL Rcf_DeAlloc_Field( input_field )

    !------------------------------------------------------------------
    ! Data from Ancillary File
    !------------------------------------------------------------------
  CASE ( Ancillary_File )

    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Ancillary field', i,                             &
                  fields_out( i ) % stashmaster % section,                    &
                  fields_out( i ) % stashmaster % item,                       &
                  fields_out( i ) % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    !------------------------------------------------------------------
    ! Data to be set to zero
    !------------------------------------------------------------------
  CASE ( Set_To_Zero )

    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Set to Zero, field', i,                          &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    SELECT CASE ( output_field % stashmaster % data_type )
    CASE ( ppx_type_real )
      output_field % DATA( :, : ) = 0.0

    CASE ( ppx_type_int )
      output_field % Data_Int( : , : ) = 0

    CASE ( ppx_type_log )
      output_field % Data_Log( : , : ) = .FALSE.
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        CALL umPrint( 'Assuming zero means FALSE for logicals',               &
            src='rcf_create_dump_mod')
      END IF

    CASE DEFAULT
      CALL umPrint( 'Cannot set fields of this type to zero!',                &
      src='rcf_create_dump_mod')
      WRITE(umMessage,'(A,I7)') 'ppx_type_? =',&
        output_field % stashmaster % data_type
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END SELECT

    !------------------------------------------------------------------
    ! Data to be set as missing
    !------------------------------------------------------------------
  CASE ( Set_To_MDI )

    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Set to MDI ', i,                                 &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    SELECT CASE ( output_field % stashmaster % data_type )
    CASE ( ppx_type_real )
      output_field % DATA( :, : ) = rmdi

    CASE ( ppx_type_int )
      output_field % Data_Int( : , : ) = imdi

    CASE DEFAULT
      CALL umPrint( 'Cannot set fields of this type to MDI!',                 &
      src='rcf_create_dump_mod')
    END SELECT

    !----------------------------------------------------------------
    ! Tracer data
    !----------------------------------------------------------------
  CASE ( Tracer_File )

    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Tracer file data ', i,                           &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    ! Get Tracer File name from env var ATRACER
    CALL get_env_var( 'ATRACER', TracerFileName)

    !     Check that Tracer File exists
    INQUIRE ( FILE=TracerFileName, EXIST=l_exist, IOSTAT=IOStatus )

    IF ( .NOT. l_exist ) THEN
      CALL umPrint( 'Tracer File does not exist.',src='rcf_create_dump_mod')
      WRITE(umMessage,'(2A)') 'File : ',TRIM(TracerFileName)
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
      ErrorStatus=10
      Cmessage = 'Tracer File does not exist.'
      CALL Ereport ( RoutineName, ErrorStatus, Cmessage )
    END IF

    CALL assign_file_unit( TracerFileName, Hdr_Aux % UnitNum,handler="portio")

    CALL Model_File_Open( Hdr_Aux % UnitNum, TracerFileName,                  &
                          read_write=ioOpenReadOnly, error=err )

    CALL Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,         &
                       tracers, output_field % stashmaster % section,         &
                       output_field % stashmaster % item, imdi, imdi)

    CALL Model_File_Close( Hdr_Aux % UnitNum, TracerFileName,                 &
                           delete=ioNoDelete, error=err )

    CALL release_file_unit ( Hdr_Aux % UnitNum, handler="portio" )
    CALL Rcf_FreeUMhdr( Hdr_Aux )

    !----------------------------------------------------------------
    ! Data to be set to constant from the namelist
    !----------------------------------------------------------------
  CASE ( Set_To_Const )
    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Set to user const', i,                           &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    SELECT CASE ( fields_out( i ) % stashmaster % data_type )
    CASE ( ppx_type_real )
      output_field % DATA( :, : ) = data_source( i ) % RConst

    CASE ( ppx_type_int )
      output_field % Data_Int( : , : ) =                                      &
                                 NINT( data_source( i ) % RConst )
    CASE ( ppx_type_log )
      IF ( data_source( i ) % Rconst > 0.5 ) THEN
        output_field % Data_Log( : ,: ) = .TRUE.
      ELSE
        output_field % Data_Log( :, : ) = .FALSE.
      END IF

    END SELECT

    !----------------------------------------------------------------
    ! User prognostics from external dump
    !----------------------------------------------------------------
  CASE ( External_Dump )

    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'User Prognostic', i,                             &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    !     Check that ancillary file exists
    INQUIRE ( FILE=data_source( i ) % Ancil_File,                             &
             EXIST=l_exist, IOSTAT=IOStatus )

    IF ( .NOT. l_exist ) THEN
      CALL umPrint(  'User Prognostic File does not exist.',                  &
          src='rcf_create_dump_mod')
      WRITE(umMessage,'(2A)') 'File : ',TRIM(data_source( i ) % Ancil_File)
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
      ErrorStatus=10
      Cmessage = 'User Prognostic File does not exist.'
      CALL Ereport ( RoutineName, ErrorStatus, Cmessage )
    END IF

    CALL assign_file_unit ( data_source(i) % Ancil_File, Hdr_Aux % UnitNum,   &
                            handler="portio" )

    CALL Model_File_Open( Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,   &
                          read_write=ioOpenReadOnly, error=err )

    CALL Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,         &
                   user_prog, data_source( i ) % Ancil_SctnC,                 &
                   data_source( i ) % Ancil_ItemC,                            &
                   output_field % stashmaster % section,                      &
                   output_field % stashmaster % item)


    CALL Model_File_Close(Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,   &
                          delete=ioNoDelete, error=err )
    CALL release_file_unit( Hdr_Aux % UnitNum, handler="portio" )
    CALL Rcf_FreeUMhdr( Hdr_Aux )

    !-----------------------------------------------------------------
    ! Calculations for fields that are missing in the input dump.
    ! These are skipped now and done in a later loop
    !-----------------------------------------------------------------
  CASE ( Field_Calcs, Field_Dependent_Calcs)
    IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Calcs. needed, skip', i,                         &
                  output_field % stashmaster % section,                       &
                  output_field % stashmaster % item,                          &
                  output_field % stashmaster % NAME
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    !-----------------------------------------------------------------
    ! Data from NetCDF File
    !-----------------------------------------------------------------
  CASE ( NetCDF_File)

     IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
        WRITE(umMessage,FORM) 'NetCDF field', i,                              &
                               output_field % stashmaster % section,          &
                               output_field % stashmaster % item,             &
                               output_field % stashmaster % NAME
        CALL umPrint(umMessage,src='rcf_create_dump_mod')
     END IF

    !-----------------------------------------------------------------
    ! Unrecognised source for data - set it to missing
    !-----------------------------------------------------------------
  CASE DEFAULT

    ErrorStatus = -10
    Cmessage = 'Source code not recognised - will set field to MDI'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,FORM) 'Unknown source', i,                              &
          output_field % stashmaster % section,                               &
          output_field % stashmaster % item, 'Set to MDI!'
      CALL umPrint(umMessage,src='rcf_create_dump_mod')
    END IF

    SELECT CASE ( fields_out(i) % stashmaster % data_type )
    CASE ( ppx_type_real )
      output_field % DATA( :, : ) = rmdi

    CASE ( ppx_type_int )
      output_field % Data_Int( : , : ) = imdi

    CASE ( ppx_type_log )
      output_field % Data_Log( :, : ) = .FALSE.
    END SELECT


  END SELECT

  !-------------------------------------------------------------------
  ! Write out data if required
  !-------------------------------------------------------------------
  IF ( Data_Source( i ) % Source /= External_Dump            .AND.            &
       Data_Source( i ) % Source /= Tracer_File              .AND.            &
       Data_Source( i ) % Source /= Field_Calcs              .AND.            &
       Data_Source( i ) % Source /= Field_Dependent_Calcs    .AND.            &
       Data_Source( i ) % Source /= NetCDF_File              .AND.            &
       Data_Source( i ) % Source /= Ancillary_File )         THEN

    CALL Rcf_Write_Field( output_field, Hdr_Out, decomp_rcf_output )

  END IF

  IF ( ( Data_Source( i ) % Source /= Ancillary_File )   .AND.                &
       ( Data_Source( i ) % Source /= NetCDF_File    ) ) THEN
    CALL Rcf_DeAlloc_Field( output_field )
  END IF
END DO

! If using nearest neighbour we dont really care about fixing data.
IF (h_int_method /= nearest_neighbour) THEN

  ! We need to always fix the winds even if we are just interpolating so lets
  ! do that now.
  IF ( output_grid % global ) THEN
    ! Find the polar and non-polar winds.  Depends on staggering.
    IF (Hdr_Out % Fixhd( FH_GridStagger ) == FH_gridstagger_endgame) THEN
      stashcode_onpole  = stashcode_v
      stashcode_offpole = stashcode_u
    ELSE
      stashcode_onpole  = stashcode_u
      stashcode_offpole = stashcode_v
    END IF

    ! Locate polar wind.
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_onpole,                    &
                     fields_out, field_count_out, pos, zero_ok_arg = .TRUE.)
    ! Locate non polar wind
    CALL Rcf_Locate( stashcode_prog_sec, stashcode_offpole,                   &
                       fields_out, field_count_out, pos1,zero_ok_arg = .TRUE. )
    ! Only continue if both wind fields are present in the dump
    IF ( pos > 0 .AND. pos1 > 0) THEN
      polar_wind => fields_out(pos)
      non_polar_wind => fields_out(pos1)
      ! Only need to perform average of polar wind if we are interpolating
      IF ( polar_wind % interp == interp_done ) THEN
        ! Allocate and read polar wind.
        CALL Rcf_Alloc_Field( polar_wind )
        CALL Rcf_Read_Field( polar_wind, hdr_out, decomp_rcf_output )

        ! Allocate and read non-polar wind.
        CALL Rcf_Alloc_Field( non_polar_wind )
        CALL Rcf_Read_Field( non_polar_wind, hdr_out, decomp_rcf_output )

        ! Perform calculation of polar wind.
        CALL rcf_polar_wind( polar_wind, non_polar_wind,                      &
             hdr_out % RealC( RC_LongSpacing) )

        ! Write out field.
        CALL Rcf_Write_Field( polar_wind, hdr_out, decomp_rcf_output )

        ! Deallocate.
        CALL Rcf_DeAlloc_Field( polar_wind )
        CALL Rcf_DeAlloc_Field( non_polar_wind )
      END IF
    END IF
  END IF


  !------------------------------------------------------------------
  ! Need to perform some field_calcs (8) before Post-Processing
  !------------------------------------------------------------------

  CALL Rcf_Pre_Process_Calcs( fields_in, fields_out, field_count_in,          &
                              field_count_out, data_source, hdr_in, hdr_out )

  IF (.NOT. l_basic_interp) THEN
    !------------------------------------------------------------------
    ! Perform post-processing
    !------------------------------------------------------------------
    CALL Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,          &
                                 fields_out, field_count_out, orog_out,       &
                                 hdr_in, hdr_out, data_source)
  END IF
END IF

!------------------------------------------------------------------
! Tidy up orography fields
!------------------------------------------------------------------
CALL Rcf_DeAlloc_Field( orog_in )
CALL Rcf_DeAlloc_Field( orog_out )
CALL Rcf_DeAlloc_Field( interp_orog )

! If using nearest neighbour we dont really care about transforming data.
! We make sure field calcs is not performed in rcf_set_data_source for
! nearest_neighbour interpolation.
!------------------------------------------------------------------
! Need to revisit the source = field_calcs (8) fields
! and the             source = field_dependent_calcs (9) fields
!------------------------------------------------------------------

CALL Rcf_Field_Calcs( fields_in, fields_out, field_count_in,                  &
                      field_count_out, data_source, hdr_in, hdr_out )

CALL Rcf_Field_Dependent_Calcs( fields_in, fields_out, field_count_in,        &
                                field_count_out, data_source, hdr_in, hdr_out )

!------------------------------------------------------------------------
! When initialising an idealised setup some fields need to be populated
! together in a non-trivial, interconnected way. These are handled
! separately, overwriting any previous data if necessary.
! -----------------------------------------------------------------------
IF (l_init_idealised) THEN
  CALL rcf_ideal_initialisation( fields_out, field_count_out, hdr_out)
END IF

!-------------------------------------------------------------------
! Rotate the winds if so required
!-------------------------------------------------------------------
IF ((h_int_active_u .OR. h_int_active_v) .AND. Output_Grid % Rotated .AND.    &
    .NOT. l_same_rotation .AND. .NOT. l_init_idealised) THEN
  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint( 'Rotating Output Winds',src='rcf_create_dump_mod')
  END IF

  CALL Rcf_Rotate( fields_out, field_count_out, Output_Grid, Hdr_Out,         &
                   decomp_rcf_output, FromStandard)
END IF

!-------------------------------------------------------------------
! Transplant Data
!-------------------------------------------------------------------
IF ( l_trans ) THEN

  CALL assign_file_unit ( transp, Hdr_Aux % UnitNum, handler="portio" )

  CALL Model_File_Open( Hdr_Aux % UnitNum, transp, read_write=ioOpenReadOnly )

  CALL Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,           &
                     transplant, imdi, imdi, imdi, imdi)


  CALL Model_File_Close( Hdr_Aux % UnitNum, transp,                           &
                         delete=ioNoDelete, error=err )
  CALL release_file_unit( Hdr_Aux % UnitNum, handler="portio" )
  CALL Rcf_FreeUMhdr( Hdr_Aux )

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Create_Dump

END MODULE Rcf_Create_Dump_Mod
