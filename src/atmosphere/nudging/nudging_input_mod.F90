! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE nudging_input_mod

! Description:
!       Input/namelist control for Nudging scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Nudging

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi,imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-------------------- Nudging control parameters ----------------------
!
REAL      :: ndg_relax_uvalue = rmdi
REAL      :: ndg_relax_vvalue = rmdi
REAL      :: ndg_relax_tvalue = rmdi
! Relaxation parameters for u, v, T

INTEGER   :: ndg_lev_bottom = imdi
INTEGER   :: ndg_lev_top    = imdi
! Levels that we nudge upon from lowest to highest

INTEGER   :: ndg_on_lev_bottom = imdi
INTEGER   :: ndg_on_lev_top    = imdi
! Number of levels to go from none to full strength nudging.

CHARACTER(LEN=256)   :: ndg_datapath = ' '
! Path to analyses files
INTEGER   :: ndg_hours_perdata = imdi
! Numbers of hours per analyses data timestep
INTEGER   :: ndg_analysis_source = imdi
! Source of Analyses data being used
!  0 = ECMWF on hybrid levels,
!  1 = ECMWF on pressure levels,
!  2 = UM analyses
!  3 = JRA Analyses
REAL      :: ndg_strat_fac = rmdi
! Factor to weaken nudging in stratosphere
LOGICAL   :: l_nudging
! True when Nudging is switched on

NAMELIST /Run_Nudging/                                              &
    ndg_relax_uvalue,ndg_relax_vvalue,ndg_relax_tvalue,             &
    ndg_lev_bottom,ndg_lev_top,ndg_on_lev_bottom,ndg_on_lev_top,    &
    ndg_datapath,ndg_hours_perdata,ndg_analysis_source,             &
    ndg_strat_fac, l_nudging

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NUDGING_INPUT_MOD'

CONTAINS

SUBROUTINE print_nlist_run_nudging()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_NUDGING'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_nudging', &
    src='nudging_input_mod')

WRITE(lineBuffer,*)'ndg_relax_uvalue = ',ndg_relax_uvalue
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_relax_vvalue = ',ndg_relax_vvalue
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_relax_tvalue = ',ndg_relax_tvalue
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_lev_bottom = ',ndg_lev_bottom
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_lev_top = ',ndg_lev_top
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_on_lev_bottom = ',ndg_on_lev_bottom
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_on_lev_top = ',ndg_on_lev_top
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_datapath = ',ndg_datapath
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_hours_perdata = ',ndg_hours_perdata
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_analysis_source = ',ndg_analysis_source
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'ndg_strat_fac = ',ndg_strat_fac
CALL umPrint(lineBuffer,src='nudging_input_mod')
WRITE(lineBuffer,*)'l_nudging = ',l_nudging
CALL umPrint(lineBuffer,src='nudging_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='nudging_input_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_nudging

#if !defined(LFRIC)
SUBROUTINE read_nml_run_nudging(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_NUDGING'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 4
INTEGER, PARAMETER :: n_int = 6
INTEGER, PARAMETER :: n_real = 4
INTEGER, PARAMETER :: n_log = 1
INTEGER, PARAMETER :: n_chars = 256

TYPE my_namelist
  SEQUENCE
  INTEGER :: ndg_lev_bottom
  INTEGER :: ndg_lev_top
  INTEGER :: ndg_on_lev_bottom
  INTEGER :: ndg_on_lev_top
  INTEGER :: ndg_hours_perdata
  INTEGER :: ndg_analysis_source
  REAL :: ndg_relax_uvalue
  REAL :: ndg_relax_vvalue
  REAL :: ndg_relax_tvalue
  REAL :: ndg_strat_fac
  LOGICAL :: l_nudging
  CHARACTER (LEN=256) :: ndg_datapath
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_real_in=n_real, n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Nudging, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Nudging", iomessage)

  my_nml % ndg_lev_bottom      = ndg_lev_bottom
  my_nml % ndg_lev_top         = ndg_lev_top
  my_nml % ndg_on_lev_bottom   = ndg_on_lev_bottom
  my_nml % ndg_on_lev_top      = ndg_on_lev_top
  my_nml % ndg_hours_perdata   = ndg_hours_perdata
  my_nml % ndg_analysis_source = ndg_analysis_source
  ! end of integers
  my_nml % ndg_relax_uvalue = ndg_relax_uvalue
  my_nml % ndg_relax_vvalue = ndg_relax_vvalue
  my_nml % ndg_relax_tvalue = ndg_relax_tvalue
  my_nml % ndg_strat_fac    = ndg_strat_fac
  ! end of reals
  my_nml % l_nudging    = l_nudging
  my_nml % ndg_datapath = ndg_datapath

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  ndg_lev_bottom      = my_nml % ndg_lev_bottom
  ndg_lev_top         = my_nml % ndg_lev_top
  ndg_on_lev_bottom   = my_nml % ndg_on_lev_bottom
  ndg_on_lev_top      = my_nml % ndg_on_lev_top
  ndg_hours_perdata   = my_nml % ndg_hours_perdata
  ndg_analysis_source = my_nml % ndg_analysis_source
  ! end of integers
  ndg_relax_uvalue = my_nml % ndg_relax_uvalue
  ndg_relax_vvalue = my_nml % ndg_relax_vvalue
  ndg_relax_tvalue = my_nml % ndg_relax_tvalue
  ndg_strat_fac    = my_nml % ndg_strat_fac
  ! end of reals
  l_nudging    = my_nml % l_nudging
  ndg_datapath = my_nml % ndg_datapath

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_nudging
#endif

END MODULE nudging_input_mod
