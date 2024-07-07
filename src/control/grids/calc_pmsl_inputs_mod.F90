! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
! Description:
!   Module containing options used to calculate pmsl
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! Code Description: Language: FORTRAN 95

MODULE calc_pmsl_inputs_mod

USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod,       ONLY: rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE
PUBLIC :: npmsl_height, l_pmsl_sor, run_calc_pmsl, print_nlist_run_calc_pmsl,  &
          check_nml_run_calc_pmsl, read_nml_run_calc_pmsl

! PMSL diagnostic
! Orographic height threshold for new pmsl calculation
REAL    :: npmsl_height = rmdi

! PMSL smoothing
! True:  Use SOR algorithm
! False: Use Jacobi algorithm
LOGICAL :: l_pmsl_sor   = .FALSE.

NAMELIST /run_calc_pmsl/ npmsl_height, l_pmsl_sor

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_PMSL_INPUTS_MOD'

CONTAINS

SUBROUTINE check_nml_run_calc_pmsl()
! Description:
!   Subroutine to check variables based on run_calc_pmsl namelist.

USE chk_opts_mod, ONLY: chk_var, def_src

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CHECK_NML_RUN_CALC_PMSL'
REAL(KIND=jprb)              :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

CALL chk_var( npmsl_height, 'npmsl_height', '[0.0:99999.0]' )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_run_calc_pmsl

SUBROUTINE print_nlist_run_calc_pmsl()

USE umPrintMgr, ONLY: umPrint

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_CALC_PMSL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_calc_pmsl',                             &
  src='calc_pmsl_inputs_mod')

WRITE(lineBuffer,'(A,F0.2)') ' npmsl_height = ',npmsl_height
CALL umPrint(lineBuffer,src='calc_pmsl_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_pmsl_sor   = ',l_pmsl_sor
CALL umPrint(lineBuffer,src='calc_pmsl_inputs_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',                        &
  src='calc_pmsl_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_calc_pmsl



SUBROUTINE read_nml_run_calc_pmsl(unit_in)

USE um_parcore,       ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist,   ONLY: setup_nml_type

IMPLICIT NONE

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_CALC_PMSL'

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage

! Set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_real = 1
INTEGER, PARAMETER :: n_log = 1

TYPE my_namelist
  SEQUENCE
  REAL    :: npmsl_height
  LOGICAL :: l_pmsl_sor
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in=n_real,               &
  n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_calc_pmsl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist run_calc_pmsl", iomessage)

  my_nml % npmsl_height = npmsl_height
  my_nml % l_pmsl_sor   = l_pmsl_sor

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  npmsl_height = my_nml % npmsl_height
  l_pmsl_sor   = my_nml % l_pmsl_sor

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_calc_pmsl

END MODULE calc_pmsl_inputs_mod
