! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE mpp_conf_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Purpose:
!    This module contains configuration variables for MPP code.
!
! Code Description:
!   Language: FORTRAN 90

USE global_2d_sums_mod, ONLY: global_sum_method
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Parameters for the different possible way of performing blocking and
! non-blocking swap bounds.

INTEGER, PARAMETER :: nb_swap_bounds_off   = 0  ! Blocking swap bounds

! Use MPI point-to-point technique
INTEGER, PARAMETER :: mpi_aao      = 10
INTEGER, PARAMETER :: mpi_nsew_wa  = 11

! Use derived-data-type (DDT) and point-to-point technique of MPI
INTEGER, PARAMETER :: ddt_aao      = 20
INTEGER, PARAMETER :: ddt_nsew_wa  = 21

! Variable that controls the way that swap bounds is performed:

! On 3D a field with large halos, i.e. greater than one
INTEGER :: swap_bounds_method_3D_large    = imdi
! On 3D a field with small halos, i.e. equals 1.
INTEGER :: swap_bounds_method_3D_small    = imdi
! On 2D a field
INTEGER :: swap_bounds_method_2D          = imdi

! Parameters to choose whether to included EW and/or NW haloes when
! calculating local subdomain boundaries
LOGICAL, PARAMETER :: include_halos_ew = .TRUE.
LOGICAL, PARAMETER :: exclude_halos_ew = .FALSE.

LOGICAL, PARAMETER :: include_halos_ns = .TRUE.
LOGICAL, PARAMETER :: exclude_halos_ns = .FALSE.

! Set up readable names for SWAP_BOUNDS argument logical.
LOGICAL, PARAMETER :: swap_field_is_vector     = .TRUE.
LOGICAL, PARAMETER :: swap_field_is_scalar     = .FALSE.

! Variable that controls the way the non-blocking is performed. 
INTEGER :: nb_swap_bounds_method    = imdi

! Variables declared with a sensible default
INTEGER :: extended_halo_size_ew    = imdi
INTEGER :: extended_halo_size_ns    = imdi
INTEGER :: gcom_coll_limit          = 72     ! MetO Broadwell - used as a
                                             ! sensible default for recon.
                                             ! Model has compulsory namelist.
LOGICAL :: gather_sync              = .TRUE. ! Do we sync in gathers?
LOGICAL :: scatter_sync             = .TRUE. ! Do we sync in scatters?

NAMELIST / nlst_mpp / extended_halo_size_ew,            &
                      extended_halo_size_ns,            &
                      gcom_coll_limit,                  &
                      gather_sync,                      &
                      scatter_sync,                     &
                      global_sum_method,                &
                      swap_bounds_method_3D_large,      &
                      swap_bounds_method_3D_small,      &
                      nb_swap_bounds_method

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPP_CONF_MOD'

CONTAINS

SUBROUTINE check_nml_nlst_mpp()

! Description:
!   Subroutine to check control variables based on options nlst_mpp namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_NML_NLST_MPP'
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check sensible values chosen
CALL chk_var( extended_halo_size_ew, 'extended_halo_size_ew', '[1:10]' )

CALL chk_var( extended_halo_size_ns, 'extended_halo_size_ns', '[1:10]' )

CALL chk_var( gcom_coll_limit, 'gcom_coll_limit', '[>=0]' )

CALL chk_var( global_sum_method, 'global_sum_method', '[1,2,3]' )

CALL chk_var( swap_bounds_method_3D_large, 'swap_bounds_method_3D_large', &
  '[0,10:11,20:21]' )
CALL chk_var( swap_bounds_method_3D_small, 'swap_bounds_method_3D_small', &
  '[0,10:11,20:21]' )
CALL chk_var( nb_swap_bounds_method,       'nb_swap_bounds_method',       &
  '[0,20,21]' )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_nlst_mpp

SUBROUTINE read_nml_nlst_mpp(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLST_MPP'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 7
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: extended_halo_size_ew
  INTEGER :: extended_halo_size_ns
  INTEGER :: gcom_coll_limit
  INTEGER :: global_sum_method
  INTEGER :: swap_bounds_method_3D_large
  INTEGER :: swap_bounds_method_3D_small
  INTEGER :: nb_swap_bounds_method
  LOGICAL :: gather_sync
  LOGICAL :: scatter_sync
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=nlst_mpp, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLST_MPP", iomessage)

  my_nml % extended_halo_size_ew       = extended_halo_size_ew
  my_nml % extended_halo_size_ns       = extended_halo_size_ns
  my_nml % gcom_coll_limit             = gcom_coll_limit
  my_nml % global_sum_method           = global_sum_method
  my_nml % swap_bounds_method_3D_large = swap_bounds_method_3D_large
  my_nml % swap_bounds_method_3D_small = swap_bounds_method_3D_small
  my_nml % nb_swap_bounds_method       = nb_swap_bounds_method
  my_nml % gather_sync                 = gather_sync
  my_nml % scatter_sync                = scatter_sync

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  extended_halo_size_ew       = my_nml % extended_halo_size_ew
  extended_halo_size_ns       = my_nml % extended_halo_size_ns
  gcom_coll_limit             = my_nml % gcom_coll_limit
  global_sum_method           = my_nml % global_sum_method
  swap_bounds_method_3D_large = my_nml % swap_bounds_method_3D_large 
  swap_bounds_method_3D_small = my_nml % swap_bounds_method_3D_small
  nb_swap_bounds_method       = my_nml % nb_swap_bounds_method
  gather_sync                 = my_nml % gather_sync
  scatter_sync                = my_nml % scatter_sync
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlst_mpp

! Return the name of the swap bounds method.
FUNCTION name_sb_method (method) RESULT (str)

IMPLICIT NONE
INTEGER, INTENT(IN) :: method
CHARACTER (len=40) :: str

SELECT CASE(method)
CASE(nb_swap_bounds_off)
  str = "off (blocking)"
CASE(mpi_aao)
  str = "mpi_aao"
CASE(mpi_nsew_wa)
  str = "mpi_nsew_wa"
CASE(ddt_aao)
  str = "ddt_aao"
CASE(ddt_nsew_wa)
  str = "ddt_nsew_wa"
CASE DEFAULT
  str = "unknown"
END SELECT

END FUNCTION name_sb_method

END MODULE mpp_conf_mod
