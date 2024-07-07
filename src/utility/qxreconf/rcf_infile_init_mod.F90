! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_Infile_Init_Mod

!  Subroutine Rcf_Files_Init - initialisation of input dump
!
! Description:
!   Open input dump, setup input header. Set up Input_Grid datatype 
!   and copy input dump if rotation required.
!
! Method:
!   Input dumps opened with File_Open. Input header setup and read in.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INFILE_INIT_MOD'

CONTAINS

SUBROUTINE Rcf_Infile_Init( hdr_in )

USE Ereport_Mod, ONLY: Ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE file_manager, ONLY: assign_file_unit

USE filenamelength_mod, ONLY: filenamelength

USE get_env_var_mod, ONLY: get_env_var

USE io

USE io_constants, ONLY: &
    ioNameProvided, ioNoDelete, ioOpenReadOnly, ioOpenReadWrite

USE lam_config_inputs_mod, ONLY: &
    polelona, polelata

USE missing_data_mod, ONLY: imdi

USE Rcf_Generate_Heights_Mod, ONLY: &
    height_gen_original,             &
    height_gen_smooth,               &
    height_gen_linear,               &
    height_gen_ECMWF_Press,          &
    height_gen_ECMWF_Hybrd

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid

USE Rcf_HeadAddress_Mod, ONLY: &
    IC_1stConstRho,         IC_BLevels,            IC_TracerLevs,    &
    IC_XLen,                IC_YLen,               IC_NumLandPoints, &
    IC_PLevels,             IC_WetLevels,          IC_NoCloudLevels, &
    IC_SoilTLevels,         IC_SoilMoistLevs,      IC_NumOzoneLevs,  &
    RC_ModelTop,                                                     &
    FH_HorizGrid_Global,    FH_HorizGrid,          IC_ConvectLevs,   &
    FH_HorizGrid_LamWrap,   FH_HorizGrid_LamWrapEQ,                  &
    IC_HeightMethod,        IC_RiverRows,          IC_RiverRowLength,&
    FH_GridStagger,         FH_GridStagger_Endgame,FH_GridStagger_A, &
    RC_PoleLat,             RC_PoleLong,                             &
    RC_FirstLat,            RC_FirstLong,                            &
    RC_LatSpacing,          RC_LongSpacing

USE Rcf_interp_weights_mod, ONLY: &
    l_limit_rotations, l_same_rotation

USE Rcf_readumhdr_mod, ONLY: rcf_readumhdr

USE rcf_nlist_recon_technical_mod, ONLY:                             &
    ainitial,                                                        &
    grib_input_dump,                                                 &
    input_dump_type

USE Rcf_UMhead_Mod, ONLY: um_header_type

USE UM_ParCore, ONLY: mype

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    PrintStatus,      &
    PrStatus_Normal

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_Header_type), INTENT(INOUT)    :: hdr_in

! Local variables
CHARACTER (LEN=*), PARAMETER            :: RoutineName='RCF_INFILE_INIT'
CHARACTER (LEN=errormessagelength)      :: Cmessage
CHARACTER (LEN=filenamelength)          :: DumpName
INTEGER                                 :: ErrorStatus
INTEGER                                 :: err
INTEGER                                 :: i
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Initialise errorstatus
errorstatus = 0
!------------------------------------------------------------------
! First open and setup Input Dump
!------------------------------------------------------------------

IF ( input_dump_type == grib_input_dump ) THEN
  CALL get_env_var('RECONTMP', DumpName)
  CALL assign_file_unit(DumpName, hdr_in % UnitNum, handler="portio")
  CALL File_Open(hdr_in % UnitNum, DumpName, &
                 read_write=ioOpenReadOnly, error=err)
ELSE
  CALL assign_file_unit(ainitial, hdr_in % UnitNum, handler="portio")
  CALL File_Open(hdr_in % UnitNum, ainitial, &
                 read_write=ioOpenReadOnly, error=err)
END IF

IF ( err /= 0) THEN
  Cmessage    = 'Failed to Open Start Dump'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN

  CALL umPrint( '',src='rcf_files_init_mod')
  IF ( input_dump_type == grib_input_dump ) THEN
    WRITE(umMessage,'(2A)') 'Input GRIB data : ',TRIM(DumpName)
    CALL umPrint(umMessage,src='rcf_files_init_mod')
  ELSE
    WRITE(umMessage,'(2A)') 'Input dump : ', TRIM(ainitial)
    CALL umPrint(umMessage,src='rcf_files_init_mod')
  END IF

END IF
CALL Rcf_ReadUMhdr( hdr_in )

!------------------------------------------------------------------
! Set the input grid data type with info from header
! Assume C grid throughout - Ocean however is B grid (?)
!------------------------------------------------------------------
IF ( hdr_in % FixHd( FH_HorizGrid ) == FH_HorizGrid_Global ) THEN
  Input_Grid % global   =  .TRUE.
  Input_Grid % Rotated  =  .FALSE.
ELSE
  Input_Grid % global   = .FALSE.
  IF (  hdr_in % RealC( RC_PoleLat ) == polelata    .AND. &
        hdr_in % RealC( RC_PoleLong) == polelona    .AND. &
        l_limit_rotations ) THEN
    IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
      CALL umPrint( 'LAM Poles are the same. Will limit rotations performed.', &
          src='rcf_files_init_mod')
    END IF
    l_same_rotation = .TRUE.
  END IF
  IF ( hdr_in % RealC( RC_PoleLat ) /= 90 .OR. &
            hdr_in % RealC( RC_PoleLong) /= 0 ) THEN
    Input_Grid % Rotated = .TRUE.
  ELSE
    Input_Grid % Rotated = .FALSE.
  END IF
END IF

Input_Grid % glob_p_row_length = hdr_in % IntC( IC_XLen )
Input_Grid % glob_p_rows       = hdr_in % IntC( IC_YLen )
! Set grid stagger information
Input_Grid % grid_stagger      = Hdr_in % fixhd(FH_GridStagger)

IF ( Hdr_in % fixhd(FH_GridStagger) == FH_GridStagger_A ) THEN
  ! GRIB data on A grid ( v rows not 1 short)
  Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
  Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )


  ! The Fixhd has no means of identifying whether a grid is a cyclic
  ! or bicyclic lam. A wrapping LAM could be E-W N-S or both and the code
  ! treats wrapping as seperate to cyclic.
  ! The RCF does not support cyclic LAMs and those users appear
  ! thankfully not to use this code.

  ! Switching logic for the LAM/global for use by ND and EG.

ELSE     ! C grid (global or LAM) Atmos
  IF (Hdr_In % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
    Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen ) + 1
  ELSE   ! New Dynamics
    Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
    Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen ) - 1
  END IF
  IF (Hdr_In % FixHd(FH_HorizGrid) == FH_HorizGrid_LamWrap .OR. &
      Hdr_In % FixHd(FH_HorizGrid) == FH_HorizGrid_LamWrapEq) THEN
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )
  END IF
END IF


Input_Grid % glob_land_field   = hdr_in % IntC( IC_NumLandPoints )

Input_Grid % model_levels      = hdr_in % IntC( IC_PLevels )

Input_Grid % cloud_levels      = hdr_in % IntC( IC_NoCloudLevels )
Input_Grid % st_levels         = hdr_in % IntC( IC_SoilTLevels )
Input_Grid % sm_levels         = hdr_in % IntC( IC_SoilMoistLevs )
Input_Grid % bl_levels         = hdr_in % IntC( IC_BLevels )
Input_Grid % ozone_levels      = hdr_in % IntC( IC_NumOzoneLevs )
Input_Grid % tr_levels         = hdr_in % IntC( IC_TracerLevs )
Input_Grid % conv_levels       = hdr_in % IntC( IC_ConvectLevs )
Input_Grid % z_top_of_model    = hdr_in % RealC( RC_ModelTop )
Input_Grid % first_constant_r_rho_level = &
                                 hdr_in % IntC( IC_1stConstRho)

! Original height generation method is set if Integer Const is
! either 1 or MDI (ie unset)

IF ( hdr_in % IntC( IC_HeightMethod ) == imdi .OR.                 &
     hdr_in % IntC( IC_HeightMethod ) == height_gen_original ) THEN

  Input_Grid % height_gen_method = height_gen_original

ELSE IF ( hdr_in % IntC( IC_HeightMethod ) == height_gen_smooth ) THEN
  Input_Grid % height_gen_method = height_gen_smooth

ELSE IF ( hdr_in % IntC( IC_HeightMethod ) == height_gen_linear ) THEN
  Input_Grid % height_gen_method = height_gen_linear

ELSE IF (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Press ) THEN
  Input_Grid % height_gen_method = height_gen_ECMWF_Press

ELSE IF (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Hybrd ) THEN
  Input_Grid % height_gen_method = height_gen_ECMWF_Hybrd

ELSE
  ErrorStatus = 40
  Cmessage = 'Input dump has unknown height generation method'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! No. of wet levels is hardwired to be the same as no. of model levels.
! Abort if dump has different no. of wet levels.
IF ( hdr_in % IntC( IC_WetLevels ) /= Input_Grid % model_levels ) THEN
  ErrorStatus = 50
  WRITE (Cmessage, '(A, I3, A, I3)')                                          &
          'No. of wet levels in dump not equal to no. of model levels'        &
          , hdr_in % IntC( IC_WetLevels ) , " != ", Input_Grid % model_levels
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! If available set the river routing rows/row_length from header.
! Otherwise set to 0 to allow decomposition to work correctly.
IF ( hdr_in % IntC( IC_RiverRows ) /= imdi ) THEN
  Input_Grid % glob_r_rows       = hdr_in % IntC( IC_RiverRows )
  Input_Grid % glob_r_row_length = hdr_in % IntC( IC_RiverRowLength )
ELSE
  Input_Grid % glob_r_rows       = 0
  Input_Grid % glob_r_row_length = 0
END IF

ALLOCATE( Input_grid % eta_theta_levels( 0:Input_Grid % model_levels))
ALLOCATE( Input_grid % eta_rho_levels( Input_Grid % model_levels))
ALLOCATE( Input_grid % rhcrit( Input_Grid % model_levels ) )
ALLOCATE( Input_grid % soil_depths( Input_Grid % sm_levels ) )

Input_grid % eta_theta_levels( 0 : Input_Grid % model_levels ) =  &
             hdr_in % LevDepC( 1 : Input_Grid % model_levels + 1, 1)

Input_grid % eta_rho_levels( 1 : Input_Grid % model_levels ) =    &
           hdr_in % LevDepC( 1 : Input_Grid % model_levels, 2)

Input_Grid % rhcrit(  1 : Input_Grid % model_levels ) =           &
    hdr_in % LevDepC( 1 : Input_Grid % model_levels, 3 )

Input_Grid % soil_depths( 1 : Input_Grid % sm_levels ) =          &
        hdr_in % LevDepC( 1 : Input_Grid % sm_levels, 4 )

Input_Grid % glob_p_field = Input_Grid % glob_p_row_length *        &
                            Input_Grid % glob_p_rows
Input_Grid % glob_u_field = Input_Grid % glob_u_row_length *        &
                            Input_Grid % glob_u_rows
Input_Grid % glob_v_field = Input_Grid % glob_v_row_length *        &
                            Input_Grid % glob_v_rows

ALLOCATE(input_grid % lambda_p(Input_Grid % glob_p_row_length))
ALLOCATE(input_grid % lambda_u(Input_Grid % glob_u_row_length))
ALLOCATE(input_grid % phi_p   (Input_Grid % glob_p_rows))
ALLOCATE(input_grid % phi_v   (Input_Grid % glob_v_rows))

! If input variable grid we can copy the dimensions across.
IF (hdr_in % Len2RowDepc == 2 .AND. hdr_in % Len2ColDepc == 2) THEN
  input_grid % lambda_p(:) = &
    hdr_in % coldepc(1:input_grid % glob_p_row_length,1)
  input_grid % lambda_u(:) = &
    hdr_in % coldepc(1:input_grid % glob_u_row_length,2)
  input_grid % phi_p   (:) = hdr_in % rowdepc(1:input_grid % glob_p_rows,1)
  input_grid % phi_v   (:) = hdr_in % rowdepc(1:input_grid % glob_v_rows,2)
ELSE
  ! Assume we need to create the grid if we are not a variable grid.
  IF ( Hdr_in % fixhd(FH_GridStagger) == FH_GridStagger_A ) THEN
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  ELSE IF (Hdr_In % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-0.5)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-0.5)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  ELSE
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-0.5)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-0.5)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  END IF
END IF

! Set the pole coordinates.
input_grid % phi_npole    = hdr_in % RealC( RC_PoleLat  )
input_grid % lambda_npole = hdr_in % RealC( RC_PoleLong )

!--------------------------------------------------------------------
! If input grid is rotated, the winds will be unrotated and
! written back to file, so will copy the file here and
! make sure the input is not overwritten
!--------------------------------------------------------------------
IF ( Input_Grid % Rotated .AND. &
     input_dump_type /= grib_input_dump  .AND. &
     .NOT. l_same_rotation) THEN

  CALL File_Close( hdr_in % UnitNum, ainitial, &
                   name_in_environ=ioNameProvided, delete=ioNoDelete, error=err)

  IF ( err /= 0 ) THEN
    Cmessage    = 'Failed to Close Start Dump'
    ErrorStatus = 50
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  IF (mype == 0 ) THEN
    CALL Shell('cp '//TRIM(ainitial)//' $RECONTMP; chmod +rw $RECONTMP',&
               34+LEN_TRIM(ainitial))
  END IF

  CALL get_env_var('RECONTMP', DumpName)

  CALL File_Open( hdr_in % UnitNum, DumpName, &
                  read_write=ioOpenReadWrite, error=err)

  IF ( err /= 0 ) THEN
    Cmessage    = 'Failed to Open Copied Start Dump'
    ErrorStatus = 50
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Infile_Init
END MODULE Rcf_Infile_Init_Mod
