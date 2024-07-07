! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  initialisation for reconfiguration

MODULE Rcf_Initialise_Mod
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!  Subroutine Rcf_Initialise - initialisation tasks
!
! Description:
! This module initialises the reconfiguration calculation,
! including Gcom, Namelists, Decompositions and File information
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INITIALISE_MOD'

CONTAINS

SUBROUTINE Rcf_Initialise( hdr_in, hdr_out )

USE mpl, ONLY:             &
    mpl_max_processor_name

USE Rcf_Decompose_Mod, ONLY: &
    Rcf_Decompose

USE UM_ParVars, ONLY: &
    nproc_x,            nproc_y,            &
    change_decomposition

USE UM_ParCore, ONLY: &
    mype,               nproc,              &
    nproc_max

USE Ereport_mod, ONLY: &
    Ereport

USE Rcf_Read_Namelists_Mod, ONLY: &
    Rcf_Read_Namelists

USE Rcf_Stash_Init_Mod, ONLY: &
    Rcf_Stash_Init

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Infile_Init_Mod, ONLY: &
    Rcf_Infile_Init

USE Rcf_Outfile_Init_Mod, ONLY: &
    Rcf_Outfile_Init

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,           &
    Output_Grid

USE Decomp_DB

USE decomp_params, ONLY: &
    decomp_rcf_input,    &
    decomp_rcf_output

USE umPrintMgr, ONLY:      &
    newline,                &
    umPrint,                &
    umMessage,              &
    umPrintSetLevel,        &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Min

USE nlstcall_mod, ONLY: &
    LTimer

USE Rcf_Set_Interp_Logicals_Mod, ONLY: &
    Rcf_Set_Interp_Logicals

USE rcf_nlist_recon_technical_mod, ONLY: &
    grib_input_dump,                     &
    input_dump_type
USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Grib_Control_Mod, ONLY: &
    Rcf_Grib_Control

USE Rcf_Grib2ff_Init_Mod, ONLY: &
    Rcf_Grib2ff_Init

USE coupling_control_mod,  ONLY:  l_oasis

USE get_env_var_mod, ONLY: get_env_var

USE rcf_lsm_mod, ONLY: &
    lsm_source

USE rcf_read_lsm_land_points_mod, ONLY: rcf_read_lsm_land_points

USE rcf_set_lsm_land_points_mod, ONLY: &
    rcf_set_lsm_land_points

USE items_nml_mod, ONLY: &
    Input_Dump,          &
    Ancillary_File,      &
    Set_To_Zero,         &
    Set_To_Const

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (um_header_type), INTENT(INOUT)  :: hdr_in
TYPE (um_header_type), INTENT(INOUT)  :: hdr_out

! Local Vars/Params
CHARACTER (LEN=*), PARAMETER    :: RoutineName = 'RCF_INITIALISE'
CHARACTER (LEN=errormessagelength)    :: Cmessage
CHARACTER (LEN=8)               :: c_nprocx
CHARACTER (LEN=8)               :: c_nprocy
CHARACTER (LEN=8)               :: c_printst
CHARACTER (LEN=8)               :: c_timer
CHARACTER (LEN=10)               :: c_coupl
INTEGER                         :: ErrorStatus
INTEGER                         :: err
INTEGER                         :: length       ! length of returned string

CHARACTER(LEN=mpl_max_processor_name) :: env_myhost
INTEGER                               :: env_myhost_len

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------
! is the RCF related to the coupled model?
!------------------------------------

CALL get_env_var('COUPLER',c_coupl, allow_missing=.TRUE., allow_empty=.TRUE.)
IF (c_coupl == "none" .OR. c_coupl == "") THEN
  l_oasis = .FALSE.
ELSE
  l_oasis = .TRUE.
END IF

!----------------------------------------------------------------
! Initialise parallel variables
!----------------------------------------------------------------
! Get the x/y decomposition
CALL get_env_var( 'RCF_NPROCX', c_nprocx, allow_missing=.TRUE., length=length )
IF ( length < 0 .OR.  c_nprocx == 'UNSET' ) THEN
  ErrorStatus = -30
  Cmessage = 'RCF_NPROCX not set: Defaults apply!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  IF ( MOD( nproc, 2 ) == 0 ) THEN
    nproc_x = 2
  ELSE
    nproc_x = 1
  END IF
ELSE
  READ(c_nprocx,'(I4)') nproc_x
END IF

CALL get_env_var( 'RCF_NPROCY', c_nprocy, allow_missing=.TRUE., length=length )
IF ( length < 0 .OR. c_nprocy == 'UNSET' ) THEN
  ErrorStatus = -40
  Cmessage = 'RCF_NPROCY not set: Defaults apply!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  nproc_y = nproc/nproc_x
ELSE
  READ(c_nprocy,'(I4)') nproc_y
END IF


! Check the x/y decomp works ok
IF ( nproc_x * nproc_y /= nproc ) THEN
  ErrorStatus = 50
  Cmessage = 'Total number of processors does not fit EW/NS LPG'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Output the decomposition to STDOUT note that PrintStatus is not
! yet set
IF (mype == 0) THEN
  CALL umPrint( '',src='rcf_initialise_mod')
  WRITE(umMessage,'(A,I7,A)') 'Parallel Reconfiguration using ', nproc, &
      ' procesor(s)'
  CALL umPrint(umMessage,src='rcf_initialise_mod')
  WRITE(umMessage,'(A,I7,A,I7)')'divided into a LPG with nproc_x=',nproc_x,   &
                            'and nproc_y=',nproc_y
  CALL umPrint(umMessage,src='rcf_initialise_mod')
  CALL umPrint( '',src='rcf_initialise_mod')
END IF

! Check on maximum size
nproc_max = nproc

WRITE(umMessage,'(I7,A)') nproc,' Processors initialised.'
CALL umPrint(umMessage,src='rcf_initialise_mod')

#if defined(RECON_SERIAL)
WRITE(umMessage,'(A,I5)') 'I am PE ',mype
CALL umPrint(umMessage,src='rcf_initialise_mod')
#else
CALL MPL_Get_processor_name(env_myhost, env_myhost_len, errorstatus)
IF (errorstatus /= 0) THEN
  WRITE(umMessage,'(A,I5)') 'I am PE ',mype
  CALL umPrint(umMessage,src='rcf_initialise_mod')
ELSE
  WRITE(umMessage,'(A,I5,A,A)') 'I am PE ',mype,' on ', TRIM(env_myhost)
  CALL umPrint(umMessage,src='rcf_initialise_mod')
END IF
#endif


!----------------------------------------------------------------
! Initialise LTimer (logical for timer)
!----------------------------------------------------------------
c_timer = 'true'
CALL get_env_var( 'RCF_TIMER', c_timer, allow_missing=.TRUE., &
                  allow_empty=.TRUE.)

! Errors are ignored and default taken
IF ( c_timer(1:4) == 'TRUE' .OR. c_timer(1:4) == 'true' ) THEN
  LTimer = .TRUE.
  IF (PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint( 'Timer is ON',src='rcf_initialise_mod')
  END IF
ELSE
  LTimer = .FALSE.
  IF (PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint( 'Timer is OFF',src='rcf_initialise_mod')
  END IF
END IF

IF (LTimer) CALL Timer( 'Reconfigure', 1)
IF (LTimer) CALL Timer( 'Initialise', 3)

!------------------------------------------------------------------
! Set Namelist information
!------------------------------------------------------------------

CALL Rcf_Read_Namelists( )

!-----------------------------------------------------------------
! Set STASH information
!-----------------------------------------------------------------

CALL Rcf_Stash_Init()

!-----------------------------------------------------------------
! Handle GRIB1 data
!-----------------------------------------------------------------

! Setup original GRIB conversion.
IF (input_dump_type == grib_input_dump) THEN
  CALL Rcf_Grib_Control( )
END IF

!-----------------------------------------------------------------
! Set input file file headers 
!-----------------------------------------------------------------

CALL Rcf_Infile_Init( hdr_in )

!-----------------------------------------------------------------
! Set Common feature between original GRIB and fieldsfile created
! from GRIB2FF utility.
!-----------------------------------------------------------------

CALL Rcf_Grib2ff_Init( hdr_in, input_grid )

!-----------------------------------------------------------------
! Set up Decompositions
!-----------------------------------------------------------------

CALL Rcf_Decompose( Output_Grid, nproc_x, nproc_y, decomp_rcf_output )

CALL Rcf_Decompose( Input_Grid, nproc_x, nproc_y, decomp_rcf_input )

CALL Change_Decomposition( decomp_rcf_input )

!-----------------------------------------------------------------
! Set interpolation logical flags
!-----------------------------------------------------------------

CALL Rcf_Set_Interp_Logicals( Input_Grid, Output_Grid, Hdr_In)

!-----------------------------------------------------------------
! Count number of land points in land sea mask if using either 
! an ancillary or just copying the mask from the input dump without interp.
! If setting the land sea mask to some other value, calculate according to
! the requirements of the items namelist.
!-----------------------------------------------------------------
SELECT CASE (lsm_source)

CASE (ancillary_file)
  CALL Rcf_Read_LSM_Land_Points( Input_Grid, Output_Grid,Hdr_in)

CASE (input_dump)
  IF (.NOT. h_int_active) THEN
    CALL Rcf_Read_LSM_Land_Points( Input_Grid, Output_Grid,Hdr_in)
  END IF

CASE (set_to_zero, set_to_const)
  CALL Rcf_Set_LSM_Land_Points( Output_Grid )

CASE DEFAULT
  ErrorStatus = 100
  Cmessage = 'Selected method for setting the land-sea mask is not supported.' &
             // newline //                                                     &
             'Supported methods are: from ancillary, from input dump,'         &
             // newline // 'set to zero, set to constant.'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!-----------------------------------------------------------------
! Set output file file headers 
!-----------------------------------------------------------------

CALL Rcf_Outfile_Init( hdr_in, hdr_out )

IF (LTimer) CALL Timer( 'Initialise', 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Initialise
END MODULE Rcf_Initialise_Mod
