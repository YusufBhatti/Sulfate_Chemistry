! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  sets up the output grid land-sea mask

MODULE Rcf_Setup_LSM_out_mod

USE coast_aj_mod, ONLY: coast_aj

IMPLICIT NONE

!  Subroutine Rcf_Setup_LSM_Out - output land-sea mask setup
!
! Description:
! This subroutine sets up the Land-Sea Mask for the output grid
! (if not ancillary) and computes gather indexes etc for coastal
! adjustment should these be required. This needs to be done for
! each PE (with whole grid) as for the weights for horizontal
! interpolation as this is where the coastal adjustment occurs.
!
! Method:
!   All pes work out the global mask (if not ancillary).
!   This is scattered to create the local mask.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_LSM_OUT_MOD'

CONTAINS


SUBROUTINE Rcf_Setup_LSM_Out( hdr_out, fields_in, field_count_in, fields_out,  &
                              field_count_out, grid_out )

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE um_stashcode_mod, ONLY: &
    stashcode_lsm,             &
    stashcode_prog_sec

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_FreeUMhdr_Mod, ONLY: &
    Rcf_FreeUMhdr

USE UM_ParVars, ONLY:   &
    current_decomp_type,&
    gc_all_proc_group,  &
    change_decomposition

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE UM_ParParams, ONLY:  &
    halo_type_single

USE Field_Types, ONLY: &
    fld_type_p

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Open_LSM_Ancil_mod, ONLY: &
    Rcf_Open_LSM_Ancil

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Ereport_Mod, ONLY: &
    Ereport

USE rcf_nlist_recon_science_mod, ONLY: &
    coast_adj_circle_method,           &
    coast_adj_method,                  &
    coast_adj_standard_method

USE file_manager, ONLY: &
    release_file_unit

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid,             &
    grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Interp_Weights_Mod, ONLY: &
    bl_index_b_l,               bl_index_b_r,              &
    bl_index_t_l,               bl_index_t_r,              &
    weight_b_l,                 weight_t_r,                &
    weight_t_l,                 weight_b_r,                &
    h_int_active

USE Rcf_Lsm_Mod, ONLY: &
    N_Coastal_Points,           coast_index_in,            &
    coast_index_out,            index_targ_land,           &
    index_targ_sea,             land_unres_index,          &
    sea_unres_index,            n_land_points_unres,       &
    n_sea_points_unres,         l_lsm_out_present,         &
    lsm_source,                 glob_lsm_in,               &
    glob_lsm_out,               local_lsm_out,             &
    glob_land_out,              local_land_out,            &
    cyclic,                     land_unres_constrain_index,&
    lsm_fixed_value,                                       &
    land_unres_notconstrain_index,                         &
    sea_unres_notconstrain_index

USE cppxref_mod, ONLY: &
    ppx_atm_compressed

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE rcf_spiral_circle_s_mod, ONLY : &
    rcf_spiral_circle_s

USE io

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE items_nml_mod, ONLY: &
    Input_Dump,          &
    Ancillary_File,      &
    Set_To_Zero,         &
    Set_To_Const

USE filenamelength_mod, ONLY:          &
    filenamelength

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)              :: field_count_in
INTEGER, INTENT(IN)              :: field_count_out
TYPE (field_type), POINTER       :: fields_in( : )
TYPE (field_type), POINTER       :: fields_out( : )
TYPE (um_header_type), INTENT(INOUT)   :: hdr_out
TYPE (grid_type), INTENT(IN)     :: grid_out

! Note that the two fields for lsm contain sizes but *NOT* data!!!
! Data is held in the lsm module.

! Local Data
INTEGER                :: i           ! Looper
INTEGER                :: ij          ! Looper
INTEGER                :: k           ! Looper
INTEGER                :: kk          ! Looper
INTEGER                :: ijmin       ! used for output regulation
INTEGER                :: pos_in      ! position in input fields array
INTEGER                :: pos_out     ! position in output fields array
INTEGER                :: msg         ! tag for comms
INTEGER                :: dump_pos_tmp
INTEGER                :: ErrorStatus
INTEGER                :: Land_Points ! local counter
INTEGER                :: sea_Points  ! local counter
INTEGER                :: Orig_Decomp ! Temporary for decomp change
INTEGER, PARAMETER     :: filenameprovided=1
INTEGER, ALLOCATABLE   :: int_lsm_in(:)  ! integer version
INTEGER, ALLOCATABLE   :: int_lsm_out(:) ! integer version
INTEGER, ALLOCATABLE   :: land_index(:) ! temp indexing
INTEGER, ALLOCATABLE   :: sea_index(:) ! temp indexing
LOGICAL, ALLOCATABLE   :: Land_Mask_temp(:) ! temp space
LOGICAL, ALLOCATABLE   :: Unres_Mask(:) ! temp space
LOGICAL, ALLOCATABLE   :: sea_Mask_temp(:) ! temp space
LOGICAL                :: is_land_point  ! temp variable to denote a land point
LOGICAL                :: l_lsm_from_ancil  ! .TRUE. = lsm from ancil file
LOGICAL                :: is_land_field
LOGICAL                :: constrained ! if the 200km constraint is to be applied
TYPE (Um_header_type)  :: hdr_anc     ! Header for ancillary
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_SETUP_LSM_OUT'
CHARACTER (LEN=errormessagelength)           :: Cmessage
TYPE (field_type), POINTER   :: lsm_in      ! input lsm field
TYPE (field_type), POINTER   :: lsm_out     ! output lsm field

CHARACTER (LEN=filenamelength) ::  lsm_filename

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------
! Locate the lsm fields from fields arrays
!----------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_in, field_count_in, pos_in, zero_ok_arg = .TRUE. )

IF (pos_in == 0) THEN
  ErrorStatus = -10
  Cmessage = 'Land-Sea Mask is not in input file'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  NULLIFY( lsm_in )
ELSE
  lsm_in => fields_in( pos_in )
END IF

CALL Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_out, field_count_out, pos_out, zero_ok_arg = .TRUE.)

IF (pos_out == 0 ) THEN
  NULLIFY( lsm_out )
  WRITE(cMessage, '(A)') "No land sea mask in output dump. Unable to " &
       // "initialise coastal adjustment"
  ErrorStatus = -20
  CALL ereport(routinename, errorstatus, cmessage)
  l_lsm_out_present = .FALSE.
ELSE
  lsm_out => fields_out( pos_out )
  l_lsm_out_present = .TRUE.
END IF

!------------------------------------------------------------------
! Set decomposition to rcf_output as all computation will be for
! output lsm
!------------------------------------------------------------------
orig_decomp = current_decomp_type
IF ( orig_decomp /= decomp_rcf_output ) THEN
  CALL Change_Decomposition( decomp_rcf_output )
END IF

!------------------------------------------------------------------
! Only need to set the output LSM et al if it is required in the
! output dump!
!-----------------------------------------------------------------
IF ( ASSOCIATED( lsm_out ) .AND. l_lsm_out_present ) THEN

  !------------------------------------------------------------------
  ! Check that output LSM memory is allocated - rather belt and braces
  !-----------------------------------------------------------------
  IF (.NOT. ALLOCATED( Local_LSM_Out ) ) THEN
    ErrorStatus = 30
    Cmessage = 'Local output LSM space not allocated for output!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  IF (.NOT. ALLOCATED( Glob_LSM_Out ) ) THEN
    ErrorStatus = 40
    Cmessage = 'Global output LSM space not allocated for output!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  !---------------------------------------------------------------
  ! Do we have to read in ancillary LSM? If so, do it here
  !---------------------------------------------------------------
  IF (lsm_source == ancillary_file) THEN

    CALL Rcf_Open_LSM_Ancil( Hdr_Anc, Output_Grid, lsm_filename )

    ! Read Local LSM
    ! Note that a temporary overwrite of the dump position is done
    ! as the lsm_out field doesn't correspond to the ancillary dump...
    dump_pos_tmp = lsm_out % dump_pos
    lsm_out % dump_pos = 1
    CALL Rcf_Alloc_Field( lsm_out )

    CALL Rcf_Read_Field( lsm_out, Hdr_Anc,  decomp_rcf_output )
    Local_Lsm_Out(:) = lsm_out % Data_Log(:,1)

    CALL Rcf_Dealloc_Field( lsm_out )
    lsm_out % dump_pos = dump_pos_tmp

    ! Do the comms to get it as a global LSM on all PEs
!DEPENDS ON: gather_field
    CALL Gather_Field( local_lsm_out, glob_lsm_out,               &
                       lsm_out % row_len, lsm_out % rows,         &
                       lsm_out % glob_row_len,                    &
                       lsm_out % glob_rows,                       &
                       fld_type_p, halo_type_single,              &
                       0, gc_all_proc_group )

    ! Note that this communication should ideally use GC_BBcast and
    ! use Kind to determine number of bytes, but we'll not bother.
    ! After all Rcf_Gather_Field (above) and the I/O assume that a
    ! Logical is the same size as a Real as is an Integer. We will do
    ! likewise.

    msg = 801
    CALL GC_IBcast( msg, lsm_out % glob_level_size, 0, nproc,  &
                    ErrorStatus, glob_lsm_out )

    IF ( ErrorStatus /= 0 ) THEN
      Cmessage = 'Problem broadcasting global land-sea mask from PE0'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    CALL File_Close( Hdr_Anc % UnitNum,lsm_filename,       &
          LEN_TRIM(lsm_filename), filenameprovided, ErrorStatus)
    IF ( ErrorStatus /= 0 ) THEN
      Cmessage = 'Problem closing Land-Sea ancillary'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    CALL release_file_unit( Hdr_Anc % UnitNum, handler="portio")
    CALL Rcf_FreeUMhdr( Hdr_Anc )

  END IF

  !----------------------------------------------------------------
  ! Have we changed resolution or LSM? If so, turn on horizontal
  ! interpolation.
  !----------------------------------------------------------------
  IF (lsm_source == ancillary_file .AND. (.NOT. h_int_active) ) THEN

    IF ( ASSOCIATED( lsm_in ) ) THEN
      IF ( lsm_out % glob_level_size /= lsm_in % glob_level_size ) THEN

        h_int_active = .TRUE.
        IF (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
          WRITE(umMessage,'(A)')'Horizontal interpolation is switched on '//&
              'because of a change in resolution'
          CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
        END IF

      ELSE

        DO i = 1, lsm_out % glob_level_size
          IF ( glob_lsm_out( i ) .neqv. glob_lsm_in( i ) ) THEN
            IF (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
              WRITE(umMessage,'(A)') 'Horizontal interpolation is '// &
                  'switched on because of a change in the Land-Sea Mask'
              CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
            END IF
            h_int_active = .TRUE.
            EXIT
          END IF
        END DO

      END IF
    END IF
  END IF


  !-------------------------------------------------------------
  ! Temporary space for integer version of lsms
  !-------------------------------------------------------------
  IF ( ASSOCIATED( lsm_in ) ) THEN
    ALLOCATE( int_lsm_in( lsm_in % glob_level_size ) )
  ELSE
    ALLOCATE( int_lsm_in( lsm_out % glob_level_size ) )
  END IF
  ALLOCATE( int_lsm_out   ( lsm_out % glob_level_size ) )
  ALLOCATE( land_mask_temp( lsm_out % glob_level_size ) )
  ALLOCATE( land_index    ( lsm_out % glob_level_size ) )
  ALLOCATE( sea_mask_temp ( lsm_out % glob_level_size ) )
  ALLOCATE( sea_index     ( lsm_out % glob_level_size ) )
  ALLOCATE( unres_mask( lsm_out % glob_level_size ) )

  !------------------------------------------------------------
  ! Set up integer versions of LSM (in and out)
  !------------------------------------------------------------

  IF ( ASSOCIATED (lsm_in ) ) THEN
    DO i = 1, lsm_in % glob_level_size
      IF ( glob_lsm_in( i ) ) THEN
        int_lsm_in( i ) = 1
      ELSE
        int_lsm_in( i ) = 0
      END IF
    END DO
  ELSE
    int_lsm_in( : ) = 0
  END IF

  ! Only possible if we had ancillary above.
  IF (lsm_source == ancillary_file) THEN
    DO i = 1, lsm_out % glob_level_size
      IF ( glob_lsm_out( i ) ) THEN
        int_lsm_out( i ) = 1
      ELSE
        int_lsm_out( i ) = 0
      END IF
    END DO
  ELSE
    int_lsm_out( : ) = 0
  END IF

  !--------------------------------------------------------------
  ! If interpolation is turned on we need to calculate gather
  ! indexes etc. for coastal adjustment
  !--------------------------------------------------------------
  SELECT CASE (lsm_source)

  CASE(ancillary_file, input_dump)

    IF ( h_int_active .AND. ASSOCIATED( lsm_in) ) THEN

      !-------------------------------------------------------------
      ! Allocate space required for gather indexes
      !-------------------------------------------------------------
          ! Global
      ALLOCATE( Coast_Index_In( lsm_out % glob_level_size ) )
      ALLOCATE( Coast_Index_Out( lsm_out % glob_level_size ) )
      ALLOCATE( Index_Targ_Land( lsm_out % glob_level_size ) )
      ALLOCATE( Index_Targ_Sea( lsm_out % glob_level_size ) )

      ! Note that the following 2 arrays are only used for non-spiral
      ! adjustment - thus spiral saves memory...
      IF ( coast_adj_method == coast_adj_standard_method ) THEN
        ALLOCATE( Land_Unres_Index( lsm_out % glob_level_size ) )
        ALLOCATE( Sea_Unres_Index( lsm_out % glob_level_size ) )
      ELSE IF ( coast_adj_method == coast_adj_circle_method ) THEN
        ALLOCATE( sea_unres_notconstrain_index( lsm_out % glob_level_size ) )
        ALLOCATE( land_unres_notconstrain_index( lsm_out % glob_level_size ) )
        ALLOCATE( land_unres_constrain_index( lsm_out % glob_level_size ) )
      END IF

      ! Set whether mask is estimated (false) or read (true)
      l_lsm_from_ancil = (lsm_source == ancillary_file)
      CALL Coast_AJ( bl_index_b_l, bl_index_b_r, bl_index_t_l,           &
                     bl_index_t_r, weight_t_r, weight_b_r,               &
                     weight_t_l, weight_b_l, lsm_in % glob_row_len,      &
                     lsm_in % glob_rows, lsm_out % glob_level_size,      &
                     int_lsm_in, int_lsm_out, coast_index_out,           &
                     coast_index_in, n_coastal_points, l_lsm_from_ancil, &
                     index_targ_sea, n_sea_points_unres, index_targ_land,&
                     n_land_points_unres )

      IF ( .NOT. l_lsm_from_ancil ) THEN  
        ! Coast_aj will have estimated a lsm for us
        DO i = 1, lsm_out % glob_level_size
          IF ( int_lsm_out( i ) == 0 ) THEN
            glob_lsm_out( i ) = .FALSE.
          ELSE
            glob_lsm_out( i ) = .TRUE.
          END IF
        END DO
      END IF

      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
        WRITE(umMessage,'('' COASTAL PTS ='',I10)')n_coastal_points
        CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
        WRITE(umMessage,'('' UNRES SEA PTS ='',I10)')n_sea_points_unres
        CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
        WRITE(umMessage,'('' UNRES LAND PTS ='',I10)')n_land_points_unres
        CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
      END IF

      ! Print out the input land-sea mask if appropriate.
      IF ( PrintStatus >= PrStatus_Diag .AND.  mype == 0) THEN
        CALL umPrint('',src='rcf_setup_lsm_out_mod')
        CALL umPrint(' Input Land Sea Mask.',src='rcf_setup_lsm_out_mod')
        ij = lsm_in % glob_row_len * (lsm_in % glob_rows - 1 ) + 1
        DO k = lsm_in % glob_rows, 1, -1
          ijmin = MIN( ij+149, ij+ lsm_in % glob_row_len - 1)
          WRITE(umMessage,'('' '',150I1)')( int_lsm_in(i),i = ij,ijmin)
          CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
          ij = ij - lsm_in % glob_row_len
        END DO
      END IF

      !------------------------------------------------------------------
      ! Set up gather indices to satify unresolved land and sea points
      !------------------------------------------------------------------

      ! Compute gather index for sea points minus unresolved points

      IF (coast_adj_method == coast_adj_standard_method) THEN
        DO i=1,lsm_out % glob_level_size
          sea_mask_temp(i)= .NOT. glob_lsm_out(i)
        END DO
        DO i=1,n_sea_points_unres
          IF (.NOT. glob_lsm_out(index_targ_sea(i))) THEN
            sea_mask_temp(index_targ_sea(i))=.FALSE.
          END IF
        END DO

        sea_points = 0
        DO i=1,lsm_out % glob_level_size
          IF (sea_mask_temp(i)) THEN
            sea_points=sea_points + 1
            sea_index(sea_points) = i
          END IF
        END DO

        IF (sea_points < 0 .OR. sea_points > lsm_out % glob_level_size)&
                                                       THEN
          WRITE(umMessage,'(A,I0)') 'sea_points = ', sea_points
          CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
          Cmessage = 'Value for sea_points is not valid'
          ErrorStatus =  45
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF

        ! Check that we have some sea points to use for unresolved points.
        IF (sea_points == 0 .AND. n_sea_points_unres > 0) THEN
          cmessage = 'Cannot resolve sea points.  Using input land points.'
          errorStatus =  -46
          CALL ereport( routinename, errorstatus, cmessage )

          ! Set n_sea_points_unres to zero and use nearest land values.
          n_sea_points_unres = 0
        END IF

        ! Assign each unresolved sea pt to nearest non-unresolved sea pt

        DO i=1,n_sea_points_unres

          IF (index_targ_sea(i) <= sea_index(1)) THEN
            sea_unres_index(i)=sea_index(1)

          ELSE IF (index_targ_sea(i) >  sea_index(sea_points)) THEN

            sea_unres_index(i)=sea_index(sea_points)

          ELSE

            DO kk=1,sea_points-1
              IF (index_targ_sea(i)  >=  sea_index(kk) .AND.        &
                  index_targ_sea(i)  <   sea_index(kk+1)) THEN
                sea_unres_index(i)=sea_index(kk)
              END IF
            END DO

          END IF
        END DO

        ! Compute gather index for land points minus unresolved points

        DO i=1,lsm_out % glob_level_size
          land_mask_temp(i)=glob_lsm_out(i)
        END DO
        DO i=1,n_land_points_unres
          IF (glob_lsm_out(index_targ_land(i))) THEN
            land_mask_temp(index_targ_land(i))=.FALSE.
          END IF
        END DO

        land_points = 0
        DO i=1,lsm_out % glob_level_size
          IF (land_mask_temp(i)) THEN
            land_points=land_points + 1
            land_index(land_points) = i
          END IF
        END DO

        IF (land_points < 0 .OR. land_points > lsm_out % glob_level_size)&
                                                       THEN
          WRITE(umMessage,'(A,I0)') 'land_points = ', land_points
          CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
          Cmessage = 'Value for land_points is not valid'
          ErrorStatus =  47
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF

        ! Since we can have land packed fields the input field maybe empty if 
        ! all sea points.  What values do we use?  Assume this should be decided
        ! in rcf_horizontal.
        IF (land_points == 0 .AND. n_land_points_unres > 0) THEN
          cmessage = 'Cannot resolve land points.  Using input sea points.'
          errorStatus =  -48
          CALL ereport( routinename, errorstatus, cmessage )

          ! Set n_land_points_unres to zero and use nearest sea values.
          n_land_points_unres = 0
        END IF

        ! Assign each unresolved land pt to nearest non-unresolved land pt

        DO i=1,n_land_points_unres

          IF (index_targ_land(i) <= land_index(1)) THEN
            land_unres_index(i)=land_index(1)
          ELSE IF (index_targ_land(i) >  land_index(land_points)) THEN
            land_unres_index(i)=land_index(land_points)
          ELSE

            DO kk=1,land_points-1
              IF (index_targ_land(i) >= land_index(kk) .AND.        &
                index_targ_land(i) <  land_index(kk+1)) THEN
                land_unres_index(i)=land_index(kk)
              END IF
            END DO

          END IF
        END DO

      ELSE IF (coast_adj_method == coast_adj_circle_method) THEN
                     ! new spiral method

        ! Make a mask so that all unresolved points on this are .TRUE. else
        ! .FALSE.
        unres_mask(:)=.FALSE.
        DO I=1,N_SEA_POINTS_UNRES
          unres_mask(INDEX_TARG_SEA(I))=.TRUE.
        END DO
        DO I=1,N_LAND_POINTS_UNRES
          unres_mask(INDEX_TARG_LAND(I))=.TRUE.
        END DO

  ! Do sea points
        is_land_field = .FALSE.
        constrained = .FALSE.
        CALL rcf_spiral_circle_s(glob_lsm_out, index_targ_sea,        &
               n_sea_points_unres, grid_out % glob_p_rows,            &
               grid_out % glob_p_row_length, output_grid % phi_p,     &
               output_grid % lambda_p, is_land_field, constrained,    &
               cyclic, unres_mask,sea_unres_notconstrain_index)

  ! Do land points unconstrained
        is_land_field = .TRUE.
        constrained = .FALSE.
        CALL rcf_spiral_circle_s(glob_lsm_out, index_targ_land,       &
               n_land_points_unres, grid_out % glob_p_rows,           &
               grid_out % glob_p_row_length, output_grid % phi_p,     &
               output_grid % lambda_p, is_land_field, constrained,    &
               cyclic, unres_mask,land_unres_notconstrain_index)

  ! Do land points constrained
        is_land_field = .TRUE.
        constrained = .TRUE.
        CALL rcf_spiral_circle_s(glob_lsm_out, index_targ_land,       &
               n_land_points_unres, grid_out % glob_p_rows,           &
               grid_out % glob_p_row_length, output_grid % phi_p,     &
               output_grid % lambda_p, is_land_field, constrained,    &
               cyclic, unres_mask,land_unres_constrain_index)

      END IF

    ELSE IF (lsm_source /= ancillary_file) THEN  ! h_int_active (or no input)

      ! Reuse the input Land-Sea Mask
      DO i = 1, lsm_out % glob_level_size
        IF ( int_lsm_in( i ) == 0 ) THEN
          glob_lsm_out( i ) = .FALSE.
        ELSE
          glob_lsm_out( i ) = .TRUE.
        END IF
      END DO

    END IF     ! h_int_active

  CASE (set_to_zero)

    DO i = 1, lsm_out % glob_level_size
      glob_lsm_out( i ) = .FALSE.
    END DO

  CASE (set_to_const)

    IF (lsm_fixed_value < 0.5) THEN
      is_land_point = .FALSE.
    ELSE
      is_land_point = .TRUE.
    END IF
    DO i = 1, lsm_out % glob_level_size
      glob_lsm_out( i ) = is_land_point
    END DO

  END SELECT


  !--------------------------------------------------------------
  ! Check number of land points against that specified in namelist
  !--------------------------------------------------------------
  glob_land_out = 0
  DO i = 1, lsm_out % glob_level_size
    IF ( glob_lsm_out( i ) ) THEN
      glob_land_out = glob_land_out + 1
    END IF
  END DO

  IF ( glob_land_out /= Output_Grid % glob_land_field) THEN

    WRITE(umMessage,*) 'Reconfiguration Error'
    CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
    WRITE(umMessage,*) 'No of land points in output land_sea mask     = ',     &
                glob_land_out
    CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
    WRITE(umMessage,*) 'No of land points specified in namelist RECON = ',     &
                Output_Grid % glob_land_field
    CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
    CALL umPrint( 'Please correct the number of land points, '//&
        'via the gui, updating LAND_FIELD within SIZES', &
        src='rcf_setup_lsm_out_mod')

    ErrorStatus = 50
    Cmessage='Number of land points does not agree with input namelist!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END IF

  !----------------------------------------------------------------
  ! Print out the output land-sea mask if appropriate.
  !----------------------------------------------------------------
  IF ( PrintStatus >= PrStatus_Diag .AND. ASSOCIATED( lsm_out) &
       .AND. mype == 0) THEN
    CALL umPrint('',src='rcf_setup_lsm_out_mod')
    CALL umPrint(' Output Land Sea Mask.',src='rcf_setup_lsm_out_mod')
    ij = lsm_out % glob_row_len * (lsm_out % glob_rows - 1 ) + 1
    DO k = lsm_out % glob_rows, 1, -1
      ijmin = MIN( ij+149, ij+ lsm_out % glob_row_len - 1)
      IF ( h_int_active ) THEN
        WRITE(umMessage,'('' '',150I1)')( int_lsm_out(i),i = ij,ijmin)
        CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
      ELSE
        WRITE(umMessage,'('' '',150I1)')( int_lsm_in(i),i = ij,ijmin)
        CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
      END IF
      ij = ij - lsm_out % glob_row_len
    END DO
  END IF

  !-------------------------------------------------------------
  ! Scatter the output grid from PE 0 (global should be same on
  ! all PEs) across LPG and setup local size.
  ! Only need to do this if NOT ancillary lsm.
  !-------------------------------------------------------------
  IF ( lsm_source /= ancillary_file ) THEN
    ! Strictly need a real field to do comms, so convert LSM to real
!DEPENDS ON: scatter_field
    CALL Scatter_Field( local_lsm_out, glob_lsm_out,              &
                        lsm_out % row_len, lsm_out % rows,        &
                        lsm_out % glob_row_len,                   &
                        lsm_out % glob_rows,                      &
                        fld_type_p, halo_type_single,             &
                        0, gc_all_proc_group )

  END IF

  !--------------------------------------------------------------
  ! Clear up allocated space etc
  !--------------------------------------------------------------
  DEALLOCATE( int_lsm_in )
  DEALLOCATE( int_lsm_out )
  DEALLOCATE( land_mask_temp )
  DEALLOCATE( land_index )
  DEALLOCATE( sea_mask_temp )
  DEALLOCATE( sea_index )
  DEALLOCATE( unres_mask )

END IF

!--------------------------------------------------------------
! Count local land-field size
!--------------------------------------------------------------
local_land_out = 0
IF ( ASSOCIATED( lsm_out ) ) THEN
  DO i = 1, lsm_out % level_size
    IF ( local_lsm_out( i ) ) THEN
      local_land_out = local_land_out + 1
    END IF
  END DO
ELSE
  local_lsm_out( : ) = .FALSE.
  glob_lsm_out ( : ) = .FALSE.
  local_land_out = 0
  glob_land_out  = 0

  IF ( glob_land_out /= Output_Grid % glob_land_field) THEN
    ErrorStatus = 51
    Cmessage='No land points calculated but non-zero in input namelist!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

IF ( PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'No. of land points is ', glob_land_out
  CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
  WRITE(umMessage,*) 'Local no. of land points is ',local_land_out,            &
                ' on PE ', mype
  CALL umPrint(umMessage,src='rcf_setup_lsm_out_mod')
END IF

Output_Grid % loc_land_field = local_land_out

!---------------------------------------------------------------
! Set local sizes for land-only fields in fields_out
!---------------------------------------------------------------
DO i = 1, field_count_out
  IF ( fields_out( i ) % stashmaster % grid_type ==                 &
                                       ppx_atm_compressed) THEN
    IF ( .NOT. ASSOCIATED(lsm_out) ) THEN
      ErrorStatus = 60
      Cmessage = 'Output field requested is land-packed but there is'//&
        ' no output Land-Sea Mask'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    fields_out( i ) % level_size = local_land_out
  END IF
END DO

!-------------------------------------------------------------
! Land Sea Mask - need to deal with this a little differently
! as data already read in and processed.
!-------------------------------------------------------------
IF ( pos_out /= 0 ) THEN
  CALL Rcf_Alloc_Field( fields_out ( pos_out ) )
  fields_out( pos_out ) % Data_Log(:,1) = local_lsm_out(:)

  CALL Rcf_Write_Field( fields_out( pos_out ), Hdr_Out, decomp_rcf_output )
  CALL Rcf_DeAlloc_Field( fields_out( pos_out ) )
END IF

!----------------------------------------------------------------
! Change decomposition back to original one
!----------------------------------------------------------------
IF ( orig_decomp /= current_decomp_type ) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_LSM_Out
END MODULE Rcf_Setup_LSM_Out_Mod
