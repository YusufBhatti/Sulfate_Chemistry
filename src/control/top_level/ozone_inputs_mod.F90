! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for ozone related switches and settings
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
MODULE ozone_inputs_mod

USE missing_data_mod, ONLY: rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------
LOGICAL :: zon_av_ozone = .FALSE.
! Ozone is held as:
! zonal average value => TRUE
! full field => FALSE

NAMELIST /run_ozone/ zon_av_ozone

!-----------------------------------------------------
! Items set in module
!-----------------------------------------------------
LOGICAL :: zon_av_tpps_ozone = .FALSE.

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OZONE_INPUTS_MOD'

CONTAINS

SUBROUTINE check_run_ozone()

IMPLICIT NONE

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
zon_av_tpps_ozone=zon_av_ozone

END SUBROUTINE check_run_ozone

SUBROUTINE print_nlist_run_ozone()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_OZONE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_ozone', &
     src='ozone_inputs_mod')

WRITE(lineBuffer,'(A,L1)') ' zon_av_ozone = ',zon_av_ozone
CALL umPrint(lineBuffer,src='ozone_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_ozone

SUBROUTINE read_nml_run_ozone(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_OZONE'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_log = 1

TYPE my_namelist
  SEQUENCE
  LOGICAL :: zon_av_ozone
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_ozone, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist run_ozone", iomessage)

  my_nml % zon_av_ozone = zon_av_ozone

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  zon_av_ozone = my_nml % zon_av_ozone

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_ozone

END MODULE ozone_inputs_mod
