! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up the orography for the reconfiguration

MODULE Rcf_Set_Orography_Mod
IMPLICIT NONE

!  Subroutine Rcf_Set_Orography - sets the orographies
!
! Description:
! This module sets up input, output and interpolated orography
! fields for interpolation (if required ) (particularly for height
! field generation.
!
! Method:
!  Read input and output orographies - if required interpolate the
!  input orography. Do orogrphic blending for LAM and check if
!  orographic changes force vertical interpolation.
!
!  Orographic blending blends from "outside to inside" of the domain
!  with the specified weight. The blending zone *includes* the
!  Rim - thus a weight of 1 should be specified to leave the Rim
!  untouched.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_OROGRAPHY_MOD'

CONTAINS

SUBROUTINE Rcf_Set_Orography( fields_in, fields_out, field_count_in, &
                              field_count_out, hdr_in, hdr_out,      &
                              data_source, orog_in, orog_out,        &
                              interp_orog )

USE Rcf_Gather_Field_Mod, ONLY: &
    Rcf_Gather_Field_Real

USE Rcf_Scatter_Field_Mod, ONLY: &
    Rcf_Scatter_Field_Real

USE umPrintMgr, ONLY:      &
    umPrint,                &
    newline,                &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParVars, ONLY:      &
    current_decomp_type,    &
    gc_all_proc_group,      &
    change_decomposition

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,               &
    Output_Grid

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE decomp_params, ONLY: &
    decomp_rcf_input,        &
    decomp_rcf_output

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active

USE Rcf_field_equals_mod, ONLY: &
    Rcf_field_equals

USE Rcf_Set_Interp_Flags_Mod, ONLY:  &
    interp_copy,                     &
    interp_h_only

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type,          &
    already_processed

USE items_nml_mod, ONLY:        &
    Input_Dump,                 &
    Ancillary_File

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_orog,            &
    stashcode_prog_sec

USE Rcf_ReadNL_Horizont_Mod, ONLY: &
    orog_blend_width,               &
    blend_weights

USE Rcf_Ideal_Set_Orography_Mod, ONLY: &
    Rcf_Ideal_Set_Orography

USE ereport_mod, ONLY: &
    ereport

USE rcf_nlist_recon_idealised_mod, ONLY: &
    l_init_idealised,                    &
    surface_type

USE rcf_ideal_surface_mod, ONLY: &
    surface_dump,                &
    surface_ancil

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER             :: fields_in(:)
TYPE( field_type ), POINTER             :: fields_out(:)
TYPE( field_type ), POINTER             :: orog_in
TYPE( field_type ), POINTER             :: orog_out
TYPE( field_type ), TARGET              :: interp_orog
TYPE( data_source_type ), INTENT(INOUT) :: data_source(:)
TYPE( um_header_type ), INTENT(IN)      :: hdr_in
TYPE( um_header_type ), INTENT(IN)      :: hdr_out
INTEGER, INTENT(IN)                     :: field_count_in
INTEGER, INTENT(IN)                     :: field_count_out

! Local variables
INTEGER                              :: i
INTEGER                              :: j
INTEGER                              :: k
INTEGER                              :: pos_in
INTEGER                              :: pos_out
INTEGER                              :: decomp_old   ! old decomposition
INTEGER                              :: stat         ! gcom status
INTEGER                              :: v_on = 0     ! vert interp flag
INTEGER                              :: icode = 0
INTEGER                              :: orog_source = 0  ! Source for orography
REAL                                 :: weight       ! blending weight
REAL, ALLOCATABLE                    :: orog_out_fullfield   ( :, : )
REAL, ALLOCATABLE                    :: orog_interp_fullfield( :, : )
TYPE (field_type)                    :: dummy

! Local parameters
INTEGER, PARAMETER                   :: ancil_orog = 1 ! Orography from ancil
INTEGER, PARAMETER                   :: dump_orog  = 2 ! Orography from dump
INTEGER, PARAMETER                   :: ideal_orog = 3 ! Idealised orography

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SET_OROGRAPHY'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  CALL umPrint( 'Processing Orography (stashcode 33) ', &
      src='rcf_set_orography_mod')
END IF

!------------------------------------------------------------------
! Find output orography
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos_out,zero_ok_arg =  .TRUE. )

! Check we have found it.
IF (pos_out == 0) THEN
  CALL umPrint( 'Output orography not found.  Attempting to continue.', &
      src='rcf_set_orography_mod')
  ALLOCATE(orog_in)
  ALLOCATE(orog_out)
  IF (ALLOCATED(blend_weights)) DEALLOCATE( blend_weights )
  GO TO 9999
END IF

!---------------------------------------------------------------
! Find origin of orography
!----------------------------------------------------------------
IF (l_init_idealised) THEN
  ! Idealised model.
  ! Check which surface type has been requested in the recon_idealised namelist.
  SELECT CASE(surface_type)

  CASE(surface_ancil)
    ! Use orography from ancillary.
    orog_source = ancil_orog
    IF (data_source( pos_out ) % source /= Ancillary_File) THEN
      icode = 10
      CALL ereport('rcf_set_orography_mod', icode,                           &
        'Idealised surface_type set to orography from ancil but' //newline// &
        'no items namelist providing the ancillary orography has been found.')
    END IF

  CASE(surface_dump)
    ! Use orography from input dump.
    orog_source = dump_orog
    IF ( data_source( pos_out ) % source /= Input_Dump) THEN
      icode = 20
      CALL ereport('rcf_set_orography_mod', icode,                            &
        'Idealised surface_type set to orography from input dump but'         &
        // newline //                                                         &
        'an items namelist providing an alternative orography has been found.'&
        // newline //                                                         &
        'Please ensure recon_idealised and items namelist settings '          & 
        // 'are consistent.')
    END IF

  CASE DEFAULT
    ! Any other preset orography type.
    orog_source = ideal_orog
    IF ( data_source( pos_out ) % source /= Input_Dump) THEN
      icode = 30
      CALL ereport('rcf_set_orography_mod', icode,                            &
        'Idealised surface_type set to a predefined ideal orography but'      &
        //newline//                                                           &
        'an items namelist providing an alternative orography has been found.'&
        // newline //                                                         &
        'Please ensure recon_idealised and items namelist settings '          &
        // 'are consistent.')
    END IF

  END SELECT

ELSE

  ! Not an idealised model.
  ! Check whether orography is from ancillary or not.
  IF ( data_source( pos_out ) % source == Ancillary_File ) THEN
    orog_source = ancil_orog
  ELSE
    orog_source = dump_orog
  END IF

END IF

! Inform the user which orography has been selected.
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  SELECT CASE (orog_source)

  CASE (ancil_orog)
    CALL umPrint( 'Using Ancillary Orography',src='rcf_set_orography_mod')

  CASE (dump_orog)
    IF (h_int_active) THEN
      CALL umPrint( 'Using interpolated Orography',src='rcf_set_orography_mod')
    ELSE
      CALL umPrint( 'Copying input Orography',src='rcf_set_orography_mod')
    END IF

  CASE (ideal_orog)
    CALL umPrint( 'Using Idealised Orography',src='rcf_set_orography_mod')

  END SELECT
END IF

!------------------------------------------------------------------
! find and read input orography
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_in, field_count_in, pos_in )
orog_in  => fields_in(pos_in)
CALL Rcf_Alloc_Field( orog_in )
CALL Rcf_Read_Field( orog_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Setup output orography
!------------------------------------------------------------------
orog_out => fields_out(pos_out)
CALL Rcf_Alloc_Field( orog_out )

! Set the sizes of interp_orog to be those of orog_out
CALL rcf_field_equals(interp_orog, orog_out)

!------------------------------------------------------------------
! Setup interpolated orography
!------------------------------------------------------------------
IF (h_int_active) THEN
  orog_in % interp = interp_h_only
ELSE
  orog_in % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( interp_orog )
CALL rcf_interpolate( orog_in, interp_orog, Input_Grid,       &
                      Output_Grid, dummy, dummy )

!-------------------------------------------------------------------
! Check for which output orography is required
!-------------------------------------------------------------------
SELECT CASE (orog_source)

CASE(ancil_orog)
  ! read the ancillary back in from output dump
  CALL Rcf_Read_Field( orog_out, Hdr_Out, decomp_rcf_output )

  !----------------------------------------------------------------
  ! Perform the Topog masking
  !----------------------------------------------------------------
  decomp_old = decomp_rcf_output
  IF (current_decomp_type /= decomp_rcf_output) THEN
    decomp_old = current_decomp_type
    CALL Change_Decomposition( decomp_rcf_output )
  END IF

  IF ( orog_blend_width > 0 ) THEN

    ALLOCATE( orog_out_fullfield(    orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )
    ALLOCATE( orog_interp_fullfield( orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )

    ! Gather orogrophies on PE 0
    ! Cannot use generic routine as fullfield is 2D so doesn't match
    ! in rank!
    CALL Rcf_Gather_Field_Real( orog_out % DATA(:,1),                &
                                orog_out_fullfield,                  &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    CALL Rcf_Gather_Field_Real( interp_orog % DATA(:,1),             &
                                orog_interp_fullfield,               &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    ! Do the orography blending on PE 0
    IF (mype == 0) THEN

      ! Northern and Southern Strips (including corners)
      DO i = 1, orog_out % glob_row_len
        DO j = 1, orog_blend_width

          ! First determine which weight to use
          ! Western corners
          IF ( i < orog_blend_width ) THEN
            weight = blend_weights( MIN(i,j) )

            ! Eastern corners
          ELSE IF ( i > orog_out % glob_row_len - orog_blend_width + 1)&
                                                                   THEN
            weight = blend_weights(                                 &
                               MIN( orog_out % glob_row_len - i + 1, j))

            ! Middle section
          ELSE
            weight = blend_weights( j )

          END IF

          ! Set the blended field for the Southern strip
          k = j
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )

          ! Set the blended field for the Northern strip
          k = orog_out % glob_rows - j + 1
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )
        END DO
      END DO

      ! Western and Eastern Strips (excluding corners)
      DO i = 1, orog_blend_width
        DO j = orog_blend_width + 1, orog_out % glob_rows -           &
                                     orog_blend_width

          ! Set the weight used
          weight = blend_weights( i )

          ! Set the blended field for the Western Strip
          k = i
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

          ! Set the blended field for the Eastern Strip
          k = orog_out % glob_row_len - i + 1
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

        END DO
      END DO

    END IF

    ! Only need to scatter the orog_out_fullfield
    CALL Rcf_Scatter_Field_Real( orog_out % DATA, orog_out_fullfield, &
                                 orog_out % row_len, orog_out % rows, &
                                 orog_out % glob_row_len,             &
                                 orog_out % glob_rows, 0,             &
                                 gc_all_proc_group )

    CALL Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )

    DEALLOCATE( orog_out_fullfield )
    DEALLOCATE( orog_interp_fullfield )

  END IF

  ! Change decomposition back
  IF ( current_decomp_type /= decomp_old ) THEN
    CALL Change_Decomposition( decomp_old )
  END IF

CASE (dump_orog)
  ! Use interpolated orography for output
  orog_out => interp_orog
  CALL Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )
  ! We have now handled the orography.
  data_source( pos_out ) % source = already_processed

CASE (ideal_orog)
  ! Use a preset idealised orography
  CALL Rcf_Ideal_Set_Orography( orog_out, Output_Grid, Hdr_Out )
  CALL Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )
  data_source( pos_out ) % source = already_processed

END SELECT

!------------------------------------------------------------------
! If there is a difference between input and output orographies,
! we need to turn on vertical interpolation
!------------------------------------------------------------------
IF ( .NOT. v_int_active ) THEN
  IF ( orog_in % glob_level_size /= orog_out % glob_level_size ) THEN
    IF (orog_source == ancil_orog) THEN
      ! Only switch on interpolation if ancillary
      v_int_active = .TRUE.
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
        CALL umPrint( 'Vertical interpolation has been switched on '//&
        'due to a change in orography',src='rcf_set_orography_mod')
      END IF
    END IF
  ELSE
    DO i = 1, orog_in % level_size
      IF ( ABS( orog_in % DATA( i, 1) - orog_out % DATA( i, 1 ) ) > &
           EPSILON( 1.0 ) ) THEN
        v_on = 1
        EXIT
      END IF
    END DO

    ! Need to make sure all PEs turn on v interp if 1 does
    CALL gc_imax( 1, nproc, stat, v_on )
    IF ( v_on == 1 ) THEN
      v_int_active = .TRUE.
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
        CALL umPrint( 'Vertical interpolation has been switched on '//&
        'due to a change in orography',src='rcf_set_orography_mod')
      END IF
    END IF
  END IF
END IF

IF (ALLOCATED(blend_weights)) DEALLOCATE( blend_weights )

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Set_Orography
END MODULE Rcf_Set_Orography_Mod
