! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sizes for many of the UM's main, dynamic data arrays
!
MODULE nlsizes_namelist_mod

! Description:
!  Module-based interface to the nlsizes namelist and associated declarations.
!  Contains the sizes needed for the dynamic allocation of the main data arrays
!  within the model.
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v10.3 programming standards.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Grids
!

USE missing_data_mod,   ONLY: imdi
USE filenamelength_mod, ONLY: filenamelength

IMPLICIT NONE

PRIVATE :: imdi

CHARACTER (LEN=*), PRIVATE, PARAMETER :: ModuleName  = "NLSIZES_NAMELIST_MOD"

SAVE

! Main sizes of fields for each submodel

! Start of variables read in from nlsizes
INTEGER :: global_row_length = imdi ! Points per global row
INTEGER :: global_rows       = imdi ! Number of global (theta) rows
INTEGER :: model_levels      = imdi ! Number of model levels
INTEGER :: land_field        = imdi ! Number of land points in field
INTEGER :: cloud_levels      = imdi ! Number of cloud-levels
INTEGER :: st_levels         = imdi ! Number of soil temperature levels
INTEGER :: bl_levels         = imdi ! Number of boundary-layer-levels
INTEGER :: ozone_levels      = imdi ! Number of ozone-levels

CHARACTER(LEN=filenamelength) :: vert_lev  = 'unset' ! Vertical levels file
CHARACTER(LEN=filenamelength) :: var_grid  = 'unset' ! Variable resolution
                                                     ! grid file
! End of variables read in from nlsizes

NAMELIST /nlsizes/                                                             &
  global_row_length, global_rows, land_field,                                  &
  model_levels, cloud_levels, st_levels,                                       &
  bl_levels, ozone_levels, vert_lev, var_grid

INTEGER :: row_length           ! Number of points per local row
INTEGER :: rows                 ! Number of local (theta) rows

INTEGER :: ntiles               ! Number of land surface tiles


! Physics-related sizes for atmosphere submodel
INTEGER :: sm_levels            ! Number of soil moisture levels
INTEGER :: tpps_ozone_levels    ! Number of tropopause-ozone-levels
INTEGER :: river_rows           ! Number of rows for river routing
INTEGER :: river_row_length     ! Row length for river routing


! Dynamics-related sizes for atmosphere submodel
INTEGER :: tr_levels            ! Number of tracer-levels
INTEGER :: tr_vars              ! Number of passive tracers
INTEGER :: tr_lbc_vars          ! Number of tracers in lbcs
INTEGER :: tr_ukca              ! Number of UKCA tracers
INTEGER :: tr_lbc_ukca          ! Number of UKCA tracer lbcs


! For Small executables
INTEGER :: tot_levels


! Grid related sizes for data structure
! Data structure sizes for atmosphere submodel
INTEGER :: a_prog_lookup        ! Number of prognostic fields
INTEGER :: a_prog_len           ! Total length of prognostic fields
INTEGER :: a_len_inthd          ! Length of integer header
INTEGER :: a_len_realhd         ! Length of real header
INTEGER :: a_len2_levdepc       ! Number of level-dependent arrays
INTEGER :: a_len2_rowdepc       ! Number of row-dependent arrays
INTEGER :: a_len2_coldepc       ! Number of column-dependent arrays
INTEGER :: a_len2_flddepc       ! Number of field arrays
INTEGER :: a_len_extcnst        ! Number of extra scalar constants
INTEGER :: a_len_cfi1           ! Length of compressed field index 1
INTEGER :: a_len_cfi2           ! Length of compressed field index 2
INTEGER :: a_len_cfi3           ! Length of compressed field index 3


! Size of super array holding all tracers
INTEGER :: super_array_size
INTEGER :: moisture_array_size
INTEGER :: turb_array_size


! Data structure sizes for atmosphere interface file control routines
INTEGER :: n_intf_a              ! Number of atmosphere interface areas
INTEGER :: max_intf_model_levels ! Max. number of model levs in all areas
INTEGER :: max_lbcrow_length     ! Max. number of lbc row length in all areas
INTEGER :: max_lbcrows           ! Max. number of lbc rows in all areas


! Data structure sizes for atmosphere boundary file control routines:
! Sizes applicable to all configurations (dumps/fieldsfile)
INTEGER :: pp_len_inthd          ! Length of PP file integer header
INTEGER :: pp_len_realhd         ! Length of PP file real    header

INTEGER :: aocpl_row_length      ! Atmos row length
INTEGER :: aocpl_p_rows          ! Atmos number of p rows


! Data structure sizes for atmosphere submodel (config dependent)
INTEGER :: a_len2_lookup         ! Total no of fields (incl diags)
INTEGER :: a_len_data            ! Total no of words of data
INTEGER :: a_len_d1              ! Total no of words in atmos D1


! Size of main data array for this configuration
INTEGER :: len_tot               ! Length of D1 array
INTEGER :: n_obj_d1_max          ! Number of objects in D1 array

! global (ie. dump version) of *_len_data
INTEGER :: global_a_len_data

INTEGER :: global_land_field     ! Global number of land points
INTEGER :: local_land_field      ! Local  number of land points

! Fundamental data sizes :
! Fundamental parameter  sizes of data structure
! Sizes applicable to all configurations (history file)

! Length of history file in dump
INTEGER, PARAMETER :: len_dumphist = 0

! Sizes applicable to all configurations (dumps/fieldsfile)
! Length of dump fixed header
INTEGER, PARAMETER :: len_fixhd = 256

! Size of a single LOOKUP header
INTEGER, PARAMETER :: len1_lookup     = 64
INTEGER, PARAMETER :: mpp_len1_lookup = 2

! Size of compressed LBC LOOKUP (only used internally and
! contains just the items which change between each set of LBCs
INTEGER, PARAMETER :: len1_lbc_comp_lookup = 8

! Sizes applicable to all configurations (stash)
! Moved to typstsz.h

INTEGER :: intf_len2_levdepc    ! 1st dim of interface out lev dep constants
INTEGER :: intf_len2_rowdepc    ! 2nd dim of interface out row dep constants
INTEGER :: intf_len2_coldepc    ! 2nd dim of interface out col dep constants

! sub-model atmosphere   :
! Data structure sizes derived from grid size
INTEGER :: a_len1_levdepc       ! 1st dim of level  dep constants
INTEGER :: a_len1_rowdepc       ! 1st dim of row    dep constants
INTEGER :: a_len1_coldepc       ! 1st dim of column dep constants
INTEGER :: a_len1_flddepc       ! 1st dim of field  dep constants

! Data structure sizes for atmosphere interface file control routines
INTEGER :: intf_lookupsa        ! No of interface lookups.

! sub-model atmosphere   : derived sizes
! derived from model grid/levels. Arakawa B-grid

! Size of fields on THETA grid ...
INTEGER :: theta_field_size     ! ... with no halos
INTEGER :: theta_off_size       ! ... with simple halos
INTEGER :: theta_halo_size      ! ... with extended halos

! Size of fields on u grid ...
INTEGER :: u_field_size         ! ... with no halos
INTEGER :: u_off_size           ! ... with simple halos
INTEGER :: u_halo_size          ! ... with extended halos

! Size of fields on v grid ...
INTEGER :: v_field_size         ! ... with no halos
INTEGER :: v_off_size           ! ... with simple halos
INTEGER :: v_halo_size          ! ... with extended halos

INTEGER :: n_rows               ! Number of V-rows
INTEGER :: n_cca_lev            ! Number of CCA levels


CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE print_nlist_nlsizes()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist nlsizes',                                   &
  src='nlsizes_namelist_mod')

WRITE(lineBuffer,*)' global_row_length = ', global_row_length
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' global_rows       = ', global_rows
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' land_field        = ', land_field
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' model_levels      = ', model_levels
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' cloud_levels      = ', cloud_levels
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' st_levels         = ', st_levels
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' bl_levels         = ', bl_levels
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' ozone_levels      = ', ozone_levels
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' vert_lev          = ', TRIM(vert_lev)
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')
WRITE(lineBuffer,*)' var_grid          = ', TRIM(var_grid)
CALL umPrint(lineBuffer,src='nlsizes_namelist_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',                        &
  src='nlsizes_namelist_mod')

END SUBROUTINE print_nlist_nlsizes

!-----------------------------------------------------------------------
#if !defined(LFRIC)
SUBROUTINE read_nml_nlsizes(unit_in)

USE um_parcore, ONLY: mype

USE setup_namelist, ONLY: setup_nml_type

USE errormessagelength_mod, ONLY: errormessagelength

USE check_iostat_mod, ONLY: check_iostat

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
CHARACTER (LEN=errormessagelength)  :: iomessage

INTEGER :: icode

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 8
INTEGER, PARAMETER :: n_chars = 2 * filenamelength

TYPE my_namelist
  SEQUENCE
  INTEGER :: global_row_length
  INTEGER :: global_rows
  INTEGER :: model_levels
  INTEGER :: land_field
  INTEGER :: cloud_levels
  INTEGER :: st_levels
  INTEGER :: bl_levels
  INTEGER :: ozone_levels
  CHARACTER(LEN=filenamelength) :: vert_lev
  CHARACTER(LEN=filenamelength) :: var_grid
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,                 &
  n_chars_in=n_chars)

IF ( mype == 0 ) THEN

  READ(UNIT=unit_in, NML=nlsizes, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist nlsizes", iomessage)

  my_nml % global_row_length = global_row_length
  my_nml % global_rows       = global_rows
  my_nml % model_levels      = model_levels
  my_nml % land_field        = land_field
  my_nml % cloud_levels      = cloud_levels
  my_nml % st_levels         = st_levels
  my_nml % bl_levels         = bl_levels
  my_nml % ozone_levels      = ozone_levels
  my_nml % vert_lev          = vert_lev
  my_nml % var_grid          = var_grid

END IF

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF ( mype /= 0 ) THEN

  global_row_length = my_nml % global_row_length
  global_rows       = my_nml % global_rows
  model_levels      = my_nml % model_levels
  land_field        = my_nml % land_field
  cloud_levels      = my_nml % cloud_levels
  st_levels         = my_nml % st_levels
  bl_levels         = my_nml % bl_levels
  ozone_levels      = my_nml % ozone_levels
  vert_lev          = my_nml % vert_lev
  var_grid          = my_nml % var_grid

END IF

CALL mpl_type_free(mpl_nml_type, icode)

END SUBROUTINE read_nml_nlsizes
#endif
!-----------------------------------------------------------------------

SUBROUTINE check_nlsizes()
! Description:
!   Subroutine to apply logic and range checking on variables based on 
!   the options set in the nlsizes namelist

USE atmos_max_sizes,        ONLY: row_length_max, rows_max, &
                                  model_levels_max
USE chk_opts_mod,           ONLY: chk_var, def_src
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,            ONLY: ereport
USE umPrintMgr,             ONLY: newline
USE parkind1,               ONLY: jprb, jpim
USE yomhook,                ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: chk_string
CHARACTER (LEN=*), PARAMETER :: RoutineName = "CHECK_NLSIZES"
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER :: icode
LOGICAL :: pass_log(8) = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

def_src = ModuleName//':'//RoutineName

WRITE(chk_string, '(A3,I0,A1)')'[2:',row_length_max,']'
CALL chk_var(global_row_length, "global_row_length", TRIM(chk_string),   &
             report_pass=pass_log(1) )
IF (MODULO(global_row_length, 2) /= 0) THEN
  WRITE(cmessage, '(A,I0,A)') "global_row_length = ", global_row_length, &
       newline // 'global_row_length must be even'
  icode = 10
  CALL ereport (def_src, icode, cmessage)
END IF

WRITE(chk_string, '(A3,I0,A1)')'[1:',rows_max,']'
CALL chk_var(global_rows, "global_rows", TRIM(chk_string),               &
             report_pass=pass_log(2) )

WRITE(chk_string, '(A3,I0,A1)')'[1:',model_levels_max,']'
CALL chk_var(model_levels, "model_levels", TRIM(chk_string) ,            &
              report_pass=pass_log(3) )

CALL chk_var(land_field, "land_field", "[-99,>=0]", report_pass=pass_log(4) )

WRITE(chk_string, '(A,I0,A)')'[1:',model_levels,']'
CALL chk_var(bl_levels, "bl_levels", chk_string, report_pass=pass_log(5) )
CALL chk_var(cloud_levels, "cloud_levels", chk_string, report_pass=pass_log(6) )
CALL chk_var(ozone_levels, "ozone_levels", chk_string, report_pass=pass_log(7) )

CALL chk_var(st_levels, "st_levels", "[>=1]", report_pass=pass_log(8) )

! var_grid not yet checked here. A string detailing the path to the Variable
! resolution grid file. Only used by recon for variable resolution dumps
! If required, value is checked before use

! vert_lev not yet checked here. A string detailing the path to the vertical
! levels namelist file. Only used by recon
! If required, value is checked before use

IF (ANY(.NOT. pass_log)) THEN
  icode = 20
  WRITE(cmessage, '(I0,A)') COUNT(.NOT. pass_log),                       &
       ' errors have been found in the nlsizes namelist'                 &
       // newline //                                                     &
       'Please see output for details and then correct the namelist'
  CALL ereport(def_src, icode, cmessage)
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_nlsizes

END MODULE nlsizes_namelist_mod
