! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Recalculates Dust bin contents.

MODULE Rcf_Change_Dust_Bins_Mod

! Description:
!     Calculates the contents of the bins for a 2 bin dust scheme
!     from the bin contents in a 6 bin dust scheme. Or from a 6 to
!     a 2 bin scheme.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CHANGE_DUST_BINS_MOD'

CONTAINS

SUBROUTINE Rcf_Change_Dust_Bins( fields_in, field_count_in,    &
                                 fields_out, field_count_out,  &
                                 field_stash_code,             &
                                 hdr_in,hdr_out,               &
                                 dust_bin )

USE Rcf_Locate_Mod, ONLY:                       &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY:                  &
    Rcf_Alloc_Field,                             &
    Rcf_Dealloc_Field

USE um_stashcode_mod, ONLY:                   &
    stashcode_prog_sec,                          &
    stashcode_dust1_mmr, stashcode_dust2_mmr,    &
    stashcode_dust3_mmr, stashcode_dust4_mmr,    &
    stashcode_dust5_mmr, stashcode_dust6_mmr,    &
    stashcode_orog

USE Rcf_Field_Type_Mod, ONLY:                   &
    Field_Type

USE Rcf_UMhead_Mod, ONLY:                       &
    um_header_type

USE umPrintMgr, ONLY:                           &
    umPrint,                                     &
    umMessage,                                   &
    PrintStatus,                                 &
    PrStatus_Normal

USE UM_ParCore, ONLY:                              &
    mype

USE decomp_params, ONLY:                     &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE Rcf_Read_Field_Mod, ONLY:                   &
    Rcf_Read_Field

USE Rcf_Interp_Weights_Mod, ONLY:               &
    h_int_active

USE Rcf_V_Int_Ctl_Mod, ONLY:                    &
    v_int_active

USE Rcf_Set_Interp_Flags_Mod, ONLY:             &
    interp_h_only,                               &
    interp_v_only,                               &
    interp_copy,                                 &
    interp_all,                                  &
    interp_no_op

USE dustbin_conversion_mod, ONLY:                &
    convert_dust_six_to_two,                     &
    convert_dust_two_to_six

USE Rcf_Grid_Type_Mod, ONLY:                    &
    Grid_Type,                                   &
    Input_Grid, Output_Grid

USE Ereport_Mod, ONLY:                          &
    Ereport

USE Rcf_Interpolate_Mod, ONLY:                  &
    Rcf_Interpolate

USE Rcf_Field_Equals_Mod, ONLY:                 &
    Rcf_Field_Equals

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( field_type ), POINTER       :: fields_in(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( um_header_type), INTENT(IN) :: hdr_in
TYPE( field_type ), INTENT(INOUT) :: dust_bin
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_count_in
INTEGER, INTENT(IN)               :: field_stash_code

! Internal variables
TYPE( field_type ), POINTER       ::  dust_bin_in
TYPE( field_type )                ::  dust_bin_1_out, dust_bin_2_out
TYPE( field_type )                ::  dust_bin_3_out, dust_bin_4_out
TYPE( field_type )                ::  dust_bin_5_out, dust_bin_6_out
TYPE( field_type ), POINTER       ::  orog_in, orog_out
TYPE( field_type )                ::  orog_interp
TYPE( field_type )                ::  dummy


INTEGER                           ::  pos, pos_in   ! position in array
INTEGER                           ::  interp_option ! interpolation to perform
INTEGER                           ::  ErrorStatus   ! Error no. to pass to ereport
INTEGER                           ::  bin           ! bin no. to recalculate

LOGICAL                           :: l_2_to_6_bin   ! logical indicating conversion

CHARACTER (LEN=errormessagelength)   :: Cmessage       ! used for EReport

! Parameters to stipulate if a field is required to be found in the field list
! of a dump.
LOGICAL, PARAMETER                :: l_not_required = .TRUE.
LOGICAL, PARAMETER                :: l_required = .FALSE.

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_CHANGE_DUST_BINS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
SELECT CASE ( field_stash_code )
CASE (stashcode_dust1_mmr)
  bin = 1

CASE (stashcode_dust2_mmr)
  bin = 2

CASE (stashcode_dust3_mmr)
  bin = 3

CASE (stashcode_dust4_mmr)
  bin = 4

CASE (stashcode_dust5_mmr)
  bin = 5

CASE (stashcode_dust6_mmr)
  bin = 6

CASE DEFAULT
  ErrorStatus = 7
  WRITE (Cmessage, '(A, I4)')                                &
        'Unrecognised STASH Code : ',                        &
        field_stash_code
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,'(A,I1,A)') 'Recalculating contents of dust bin ',bin,'.'
  CALL umPrint(umMessage,src='rcf_change_dust_bins_mod')
END IF

!----------------------------------------------------------------------
! Setup Orographies
!----------------------------------------------------------------------
! orogoraphy in input dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,       &
                 fields_in,field_count_in,pos,l_required)
orog_in => fields_in(pos)
CALL Rcf_Alloc_Field( orog_in )
CALL Rcf_Read_Field( orog_in, hdr_in, decomp_rcf_input )

! orogoraphy in output dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,       &
                 fields_out,field_count_out,pos,l_required)
orog_out => fields_out(pos)
CALL Rcf_Alloc_Field( orog_out )
CALL Rcf_Read_Field( orog_out, hdr_out, decomp_rcf_output )

! Set up input orograohy interploated onto output grid.
CALL Rcf_Field_Equals( orog_interp, orog_out )

IF (h_int_active) THEN
  orog_in % interp   = interp_h_only
ELSE
  orog_in % interp   = interp_copy
END IF

orog_interp % interp = interp_no_op

CALL Rcf_Alloc_Field( orog_interp )
CALL Rcf_Interpolate( orog_in, orog_interp, Input_Grid, Output_Grid, &
                      dummy, dummy )

!----------------------------------------------------------------------
! Set up the interpolation for the incoming fields
!----------------------------------------------------------------------

IF (v_int_active .AND. h_int_active ) THEN
  interp_option = interp_all

ELSE IF (v_int_active) THEN
  interp_option = interp_v_only

ELSE IF (h_int_active) THEN
  interp_option = interp_h_only

ELSE
  interp_option = interp_copy

END IF

!----------------------------------------------------------------------
! Find required fields in the dumps and read them in and interploate.
!----------------------------------------------------------------------

! Bin 1 :
! bin 1 in input dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,       &
                 fields_in,field_count_in,pos,l_required)
dust_bin_in => fields_in(pos)
CALL Rcf_Alloc_Field( dust_bin_in )
CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

! bin 1 in output dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,       &
                 fields_out,field_count_out,pos,l_required)
CALL Rcf_Field_Equals( dust_bin_1_out, fields_out(pos) )
CALL Rcf_Alloc_Field( dust_bin_1_out )

dust_bin_in % interp = interp_option
CALL Rcf_Interpolate( dust_bin_in, dust_bin_1_out, input_grid,  &
                      output_grid, orog_interp, orog_out )
CALL Rcf_Dealloc_Field( dust_bin_in )

! Bin 2 :
! bin 2 in input dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust2_mmr,       &
                 fields_in,field_count_in,pos,l_required)
dust_bin_in => fields_in(pos)
CALL Rcf_Alloc_Field( dust_bin_in )
CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

! bin 2 in output dump, must exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust2_mmr,       &
                 fields_out,field_count_out,pos,l_required)
CALL Rcf_Field_Equals( dust_bin_2_out, fields_out(pos) )
CALL Rcf_Alloc_Field( dust_bin_2_out )

dust_bin_in % interp = interp_option
CALL Rcf_Interpolate( dust_bin_in, dust_bin_2_out, input_grid,  &
                      output_grid, orog_interp, orog_out )
CALL Rcf_Dealloc_Field( dust_bin_in )

! Bin 3 :
! bin 3 in input dump, may not exist...
CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,       &
                 fields_in,field_count_in,pos,l_not_required)
IF ( pos == 0 ) THEN
  ! Must be moving from 2 to 6 bins
  l_2_to_6_bin = .TRUE.
ELSE
  l_2_to_6_bin = .FALSE.
END IF

! If changing from 2 to 6 bins : locate out bin, allocate and set to zero
! Otherwise : locate in bin, alloc, read, interpolate to dummy out.

IF ( l_2_to_6_bin ) THEN
  ! bin 3 in output dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,       &
                   fields_out,field_count_out,pos,l_required)
  CALL Rcf_Field_Equals( dust_bin_3_out, fields_out(pos) )
  CALL Rcf_Alloc_Field( dust_bin_3_out )
  dust_bin_3_out % DATA (:,:) = 0.0

  ! bin 4 in input dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust4_mmr,       &
                   fields_in,field_count_in,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 50
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 2 to 6 bins yet found bin 4 in input dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! bin 4 in output dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust4_mmr,       &
                   fields_out,field_count_out,pos,l_required)
  CALL Rcf_Field_Equals( dust_bin_4_out, fields_out(pos) )
  CALL Rcf_Alloc_Field( dust_bin_4_out )
  dust_bin_4_out % DATA (:,:) = 0.0

  ! bin 5 in input dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust5_mmr,       &
                   fields_in,field_count_in,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 60
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 2 to 6 bins yet found bin 5 in input dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! bin 5 in output dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust5_mmr,       &
                   fields_out,field_count_out,pos,l_required)
  CALL Rcf_Field_Equals( dust_bin_5_out, fields_out(pos) )
  CALL Rcf_Alloc_Field( dust_bin_5_out )
  dust_bin_5_out % DATA (:,:) = 0.0

  ! bin 6 in input dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust6_mmr,       &
                   fields_in,field_count_in,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 70
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 2 to 6 bins yet found bin 6 in input dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  ! bin 6 in output dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust6_mmr,       &
                   fields_out,field_count_out,pos,l_required)
  CALL Rcf_Field_Equals( dust_bin_6_out, fields_out(pos) )
  CALL Rcf_Alloc_Field( dust_bin_6_out )
  dust_bin_6_out % DATA (:,:) = 0.0

  CALL convert_dust_two_to_six  (                                       &
                 dust_bin_1_out % level_size * dust_bin_1_out % levels, &
                 dust_bin_1_out % DATA, dust_bin_2_out % DATA,          &
                 dust_bin_3_out % DATA, dust_bin_4_out % DATA,          &
                 dust_bin_5_out % DATA, dust_bin_6_out % DATA           &
                                )

ELSE ! going from 6 to 2 bins.
  ! bin 3 in input dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,       &
                   fields_in,field_count_in,pos_in,l_required)
  dust_bin_in => fields_in(pos_in)
  CALL Rcf_Alloc_Field( dust_bin_in )
  CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

  ! bin 3 in output dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,       &
                   fields_out,field_count_out,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 100
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 6 to 2 bins yet found bin 3 in output dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  CALL Rcf_Field_Equals( dust_bin_3_out, dust_bin_1_out )
  dust_bin_3_out % interp = interp_no_op
  dust_bin_3_out % dump_pos = 0
  dust_bin_3_out % stashmaster => fields_in(pos_in) % stashmaster
  CALL Rcf_Alloc_Field( dust_bin_3_out )

  dust_bin_in % interp = interp_option
  CALL Rcf_Interpolate( dust_bin_in, dust_bin_3_out, input_grid,  &
                        output_grid, orog_interp, orog_out )
  CALL Rcf_Dealloc_Field( dust_bin_in )

  ! bin 4 in input dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust4_mmr,       &
                   fields_in,field_count_in,pos_in,l_required)
  dust_bin_in => fields_in(pos_in)
  CALL Rcf_Alloc_Field( dust_bin_in )
  CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

  ! bin 4 in output dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust4_mmr,       &
                   fields_out,field_count_out,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 110
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 6 to 2 bins yet found bin 4 in output dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  CALL Rcf_Field_Equals( dust_bin_4_out, dust_bin_1_out )
  dust_bin_4_out % interp = interp_no_op
  dust_bin_4_out % dump_pos = 0
  dust_bin_4_out % stashmaster => fields_in(pos_in) % stashmaster
  CALL Rcf_Alloc_Field( dust_bin_4_out )

  dust_bin_in % interp = interp_option
  CALL Rcf_Interpolate( dust_bin_in, dust_bin_4_out, input_grid,  &
                        output_grid, orog_interp, orog_out )
  CALL Rcf_Dealloc_Field( dust_bin_in )

  ! bin 5 in input dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust5_mmr,       &
                   fields_in,field_count_in,pos_in,l_required)
  dust_bin_in => fields_in(pos_in)
  CALL Rcf_Alloc_Field( dust_bin_in )
  CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

  ! bin 5 in output dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust5_mmr,       &
                   fields_out,field_count_out,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 120
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 6 to 2 bins yet found bin 5 in output dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  CALL Rcf_Field_Equals( dust_bin_5_out, dust_bin_1_out )
  dust_bin_5_out % interp = interp_no_op
  dust_bin_5_out % dump_pos = 0
  dust_bin_5_out % stashmaster => fields_in(pos_in) % stashmaster
  CALL Rcf_Alloc_Field( dust_bin_5_out )

  dust_bin_in % interp = interp_option
  CALL Rcf_Interpolate( dust_bin_in, dust_bin_5_out, input_grid,  &
                        output_grid, orog_interp, orog_out )
  CALL Rcf_Dealloc_Field( dust_bin_in )

  ! bin 6 in input dump, must exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust6_mmr,       &
                   fields_in,field_count_in,pos_in,l_required)
  dust_bin_in => fields_in(pos_in)
  CALL Rcf_Alloc_Field( dust_bin_in )
  CALL Rcf_Read_Field( dust_bin_in, hdr_in, decomp_rcf_input )

  ! bin 6 in output dump, must not exist...
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust6_mmr,       &
                   fields_out,field_count_out,pos,l_not_required)
  IF ( pos /= 0 ) THEN
    ErrorStatus = 130
    WRITE (Cmessage, '(A)')                                          &
          'Changing from 6 to 2 bins yet found bin 6 in output dump'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  CALL Rcf_Field_Equals( dust_bin_6_out, dust_bin_1_out )
  dust_bin_6_out % interp = interp_no_op
  dust_bin_6_out % dump_pos = 0
  dust_bin_6_out % stashmaster => fields_in(pos_in) % stashmaster
  CALL Rcf_Alloc_Field( dust_bin_6_out )

  dust_bin_in % interp = interp_option
  CALL Rcf_Interpolate( dust_bin_in, dust_bin_6_out, input_grid,  &
                        output_grid, orog_interp, orog_out )
  CALL Rcf_Dealloc_Field( dust_bin_in )

  CALL convert_dust_six_to_two  (                                       &
                 dust_bin_1_out % level_size * dust_bin_1_out % levels, &
                 dust_bin_1_out % DATA, dust_bin_2_out % DATA,          &
                 dust_bin_3_out % DATA, dust_bin_4_out % DATA,          &
                 dust_bin_5_out % DATA, dust_bin_6_out % DATA           &
                                )

END IF ! IF ( l_2_to_6_bin )

!----------------------------------------------------------------------
! Return the bin to be written to dump.
!----------------------------------------------------------------------
SELECT CASE ( field_stash_code )
CASE (stashcode_dust1_mmr)
  dust_bin % DATA (:,:) = dust_bin_1_out % DATA (:,:)

CASE (stashcode_dust2_mmr)
  dust_bin % DATA (:,:) = dust_bin_2_out % DATA (:,:)

CASE (stashcode_dust3_mmr)
  dust_bin % DATA (:,:) = dust_bin_3_out % DATA (:,:)

CASE (stashcode_dust4_mmr)
  dust_bin % DATA (:,:) = dust_bin_4_out % DATA (:,:)

CASE (stashcode_dust5_mmr)
  dust_bin % DATA (:,:) = dust_bin_5_out % DATA (:,:)

CASE (stashcode_dust6_mmr)
  dust_bin % DATA (:,:) = dust_bin_6_out % DATA (:,:)

CASE DEFAULT
  ErrorStatus = 300
  WRITE (Cmessage, '(A, I4)')                                &
        'Unrecognised STASH Code : ',                        &
        field_stash_code
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( dust_bin_6_out )
CALL Rcf_Dealloc_Field( dust_bin_5_out )
CALL Rcf_Dealloc_Field( dust_bin_4_out )
CALL Rcf_Dealloc_Field( dust_bin_3_out )
CALL Rcf_Dealloc_Field( dust_bin_2_out )
CALL Rcf_Dealloc_Field( dust_bin_1_out )
CALL Rcf_Dealloc_Field( orog_interp )
CALL Rcf_Dealloc_Field( orog_out )
CALL Rcf_Dealloc_Field( orog_in )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Change_Dust_Bins
END MODULE Rcf_Change_Dust_Bins_Mod
