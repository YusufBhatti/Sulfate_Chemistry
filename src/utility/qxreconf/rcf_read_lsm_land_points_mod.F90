! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads land sea mask ancillary/dump to count the number of land points
!
MODULE Rcf_Read_LSM_Land_Points_mod

IMPLICIT NONE

! Description:
!   This subroutine reads the land sea mask either from an ancillary 
!   (if used) or the input dump (if the mask comes from there and
!   there is no horizontal interpolation). This is used to set
!   the number of land points in the output dump etc. 
!   Note: this is not where the land sea mask is read/calculated for 
!   putting into the dump - this is done later in Rcf_Setup_LSM_Out.
!
! Method:
!   Read land mask from ancillary/dump and count number of points
!   which are .true. (i.e. contain at least some land)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READ_LSM_LAND_POINTS_MOD'

CONTAINS

! Subroutine interface:
SUBROUTINE Rcf_Read_LSM_Land_Points( Input_Grid, Output_Grid, Hdr_in)

USE decomp_params, ONLY:               &
    decomp_rcf_input,                  &
    decomp_rcf_output

USE Ereport_Mod, ONLY:                 &
    Ereport

USE errormessagelength_mod, ONLY:      &
    errormessagelength

USE file_manager, ONLY:                &
    release_file_unit

USE filenamelength_mod, ONLY:          &
    filenamelength

USE io, ONLY:                          &
    file_close

USE io_constants, ONLY:                &
    ioNoDelete

USE nlsizes_namelist_mod, ONLY:        &
    land_field

USE Rcf_Address_Mod, ONLY:             &
    Rcf_Address

USE Rcf_Alloc_Field_Mod, ONLY:         &
    Rcf_Alloc_Field,                   &
    Rcf_Dealloc_Field

USE Rcf_Field_Type_Mod, ONLY:          &
    field_type

USE Rcf_FreeUMhdr_Mod, ONLY:           &
    Rcf_FreeUMhdr

USE Rcf_Grid_Type_Mod, ONLY:           &
    grid_type

USE Rcf_Open_LSM_Ancil_mod, ONLY: &
    Rcf_Open_LSM_Ancil

USE Rcf_Locate_Mod, ONLY:              &
    Rcf_Locate

USE Rcf_Lsm_Mod, ONLY:                 &
    glob_land_out,                     &
    lsm_source

USE Rcf_Read_Field_Mod, ONLY:          &
    Rcf_Read_Field

USE Rcf_Setup_Field_mod, ONLY:         &
    Rcf_Setup_Field

USE Rcf_UMhead_Mod, ONLY:              &
    um_header_type

USE UM_ParVars, ONLY:                  &
    change_decomposition,              &
    current_decomp_type

USE UM_ParCore, ONLY:                  &
    mype,                              &
    nproc

USE um_stashcode_mod, ONLY:            &
    stashcode_lsm,                     &
    stashcode_prog_sec

USE items_nml_mod, ONLY :              &
    Ancillary_File

USE umPrintMgr, ONLY:                  &
    umPrint,                           &
    umMessage,                         &
    PrintStatus,                       &
    PrStatus_Normal

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (grid_type), TARGET, INTENT(IN)      :: Input_Grid
TYPE (um_header_type), TARGET, INTENT(IN) :: Hdr_In
TYPE (grid_type), TARGET, INTENT(INOUT)   :: Output_Grid

! Local Data
INTEGER                      :: i           ! Looper
INTEGER                      :: pos         ! Dump/anc position of mask
INTEGER                      :: errorstatus
INTEGER                      :: icode       ! return code for GC routine
INTEGER                      :: orig_decomp ! Domain decomp at start
INTEGER                      :: decomp_lsm  ! Domain decomp of dump/anc
INTEGER, PARAMETER           :: filenameprovided=1
TYPE (um_header_type),TARGET :: hdr_anc     ! Header for ancillary  
CHARACTER (LEN=*), PARAMETER ::                                   &
                   RoutineName = 'RCF_READ_LSM_LAND_POINTS'
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=20)           :: grid_title

TYPE (field_type), POINTER      :: fields_in_file (:)  ! dump/anc    
TYPE (grid_type), POINTER       :: LSM_Grid    ! grid in dump/anc
TYPE (um_header_type), POINTER  :: Hdr_LSM  ! headers of dump/anc

INTEGER              :: glob_land_points    ! land_points on all PEs
INTEGER              :: field_count

CHARACTER (LEN=filenamelength) :: lsm_filename
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------
!   Set variables for mask headers/grid
!----------------------------------------------------------------
IF (lsm_source == Ancillary_File) THEN
  
!----------------------------------------------------------------
!   If reading mask from an ancillary, we need to open the file
!----------------------------------------------------------------

  CALL Rcf_Open_LSM_Ancil( Hdr_Anc, Output_Grid, lsm_filename )
  Hdr_LSM => Hdr_Anc
  decomp_lsm = decomp_rcf_output
  LSM_Grid => Output_Grid
  grid_title = 'Output grid'

ELSE
!----------------------------------------------------------------
!   If reading from the dump, the file is already open
!----------------------------------------------------------------

  IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
    WRITE(umMessage,'(A)')                                        &
     'Reading in Land-Sea Mask from input dump to count land_field'
    CALL umPrint(umMessage,src=RoutineName)
  END IF
  Hdr_LSM => Hdr_In
  decomp_lsm = decomp_rcf_input
  LSM_Grid => Input_Grid
  grid_title = 'Input grid'
  
END IF

!----------------------------------------------------------------
!   Set decomposition to that of the relevant file 
!----------------------------------------------------------------

orig_decomp = current_decomp_type
IF ( orig_decomp /= decomp_lsm ) THEN
  CALL Change_Decomposition( decomp_lsm )
END IF

!----------------------------------------------------------------
! Read mask from file
!----------------------------------------------------------------

NULLIFY(fields_in_file)

CALL Rcf_Setup_Field( fields_in_file, Hdr_LSM, LSM_Grid,          &
                      field_count,grid_title)

! Initialise the number of land points to zero
glob_land_points = 0

!----------------------------------------------------------------
!   Determine the position of the LSM in the file and 
!   allocate the field
!----------------------------------------------------------------

CALL Rcf_Locate( stashcode_prog_sec, stashcode_lsm,               &
                   fields_in_file, field_count, pos, zero_ok_arg = .TRUE. )

IF (pos == 0) THEN ! Behaviour with missing mask depends on 
                   ! where the mask is coming from

  IF (lsm_source == Ancillary_File ) THEN ! If ancillary, then this is 
                                          ! clearly an error
      
    errorstatus = 100
    cmessage = 'Ancillary file does not contain land-sea mask'
    CALL Ereport( RoutineName, errorstatus, cmessage )
      
  END IF

ELSE ! Mask is in file, so read it.

  CALL Rcf_Alloc_Field( fields_in_file(pos) )

!---------------------------------------------------------------
!   Read field on all PEs (using relevant decompositon)
!---------------------------------------------------------------
  
  CALL Rcf_Read_Field( fields_in_file(pos), Hdr_LSM, decomp_lsm )

!--------------------------------------------------------------
! Count the number of land points in local land sea mask
!--------------------------------------------------------------

  DO i = 1, SIZE(fields_in_file(pos) % Data_Log(:,1))
    IF ( fields_in_file(pos) % Data_Log(i,1) ) THEN
      glob_land_points = glob_land_points + 1
    END IF
  END DO
  
  ! global sum of the number of land points
  CALL gc_isum(1, nproc, icode, glob_land_points)  
  IF ( icode /= 0 ) THEN
    ! return error if there are problems in gc_isum
    errorstatus = 101
    WRITE(cmessage, '(A,A,I0)') 'Problem when computing global sum over ', &
                                'land points. Return code from gc_isum = ',&
                                icode
    CALL Ereport( RoutineName, errorstatus, cmessage )
  END IF

  CALL Rcf_Dealloc_Field( fields_in_file(pos) )

END IF ! (pos == 0)

IF (lsm_source == Ancillary_file ) THEN ! Close ancil and dealloc hdr

  CALL File_Close( Hdr_Anc % UnitNum,lsm_filename,                &
       LEN_TRIM(lsm_filename), filenameprovided, ioNoDelete,      &
       errorstatus)
  IF ( errorstatus /= 0 ) THEN
    cmessage = 'Problem closing Land-Sea ancillary'
    CALL Ereport( RoutineName, errorstatus, cmessage )
  END IF

  CALL release_file_unit ( Hdr_Anc % UnitNum, handler="portio" )

  CALL Rcf_FreeUMhdr( Hdr_Anc )

END IF

!-----------------------------------------------------------------
! Reset to original decomposition
!-----------------------------------------------------------------
IF ( orig_decomp /= current_decomp_type ) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
  WRITE(umMessage,'(A)')                                          &
    'Setting number of land points to value from that file'
  CALL umPrint(umMessage,src=RoutineName)
  
  WRITE(umMessage,'(A,I12)')                                      &
     'No. of land points is ', glob_land_points
  CALL umPrint(umMessage,src=RoutineName)
END IF

!--------------------------------------------------------------
! Use this to override the value set in the nlsizes 
! namelist and the mask/grid values set from this
!--------------------------------------------------------------

land_field = glob_land_points
glob_land_out = glob_land_points
Output_Grid % glob_land_field = glob_land_points

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Read_LSM_Land_Points
END MODULE Rcf_Read_LSM_Land_Points_mod
