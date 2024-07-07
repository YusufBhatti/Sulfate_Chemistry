! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Control variables

!  Description:

! Migrated from include file cntlall.h

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstcall_mod

USE missing_data_mod, ONLY: imdi, rmdi
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Array holding original data time (prior to assimilation)
INTEGER:: model_basis_time(6) = imdi

! Model analysis time in hours since Basis Time
! UM6.5 - Replace model_analysis_hrs by model_analysis_mins
!         model_analysis_hrs changed to REAL
REAL :: model_analysis_hrs     = rmdi
INTEGER :: model_analysis_mins = imdi

INTEGER :: ncpu              = imdi  ! No of CPUs assigned to the program
INTEGER :: ancil_reftime(6)  = imdi  ! Ref. time for updating ancillaries
INTEGER :: run_target_end(6) = imdi  ! Target end time for this run

INTEGER :: Num_ALBCs            = imdi ! Number of atmos boundary files
INTEGER :: ALBC2_StartTime_mins = imdi ! VT of first block of data in 2nd
                                       ! atmos boundary file, in minutes
                                       ! from start of run

LOGICAL :: latmosnext              = .FALSE.
LOGICAL :: loceannext              = .FALSE.
LOGICAL :: lpp                     = .FALSE. ! Activate PPCTL
LOGICAL :: lnc                     = .FALSE. ! Activate NCFILE
LOGICAL :: ldump                   = .FALSE. ! Activate DUMPCTL
LOGICAL :: lmean                   = .FALSE. ! Activate MEANCTL
LOGICAL :: lhistory                = .FALSE. ! Update TEMP history file
LOGICAL :: lprint                  = .FALSE. ! Activate PRINTCTL
LOGICAL :: lexit                   = .FALSE. ! Activate EXITCHEK

! Select printed diags from means
LOGICAL :: lmeanPR(4)              = .FALSE.

LOGICAL :: lancillary              = .FALSE. ! Activate UP_ANCIL
LOGICAL :: lboundary               = .FALSE. ! Activate UP_BOUND
LOGICAL :: lassimilation           = .FALSE. ! Activate assimilation

!Activate detailed TIMER routine
LOGICAL :: ltimer                  = .FALSE.
LOGICAL :: ltimers_user            = .FALSE.

!Activate stash/dump/mean timer for coupled model profiling
LOGICAL :: lstashdumptimer         = .FALSE.

! T : D1 copied to memory for AO coupling
LOGICAL :: l_ao_d1_memory          = .FALSE.

! Real-period climate means
LOGICAL :: lclimrealyr             = .FALSE.

LOGICAL :: l_fastrun =.FALSE. ! do not run IAU but 
                              ! maintain time conventions of DA
! 360-day calendar.
LOGICAL :: lcal360 = .FALSE.

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

NAMELIST /nlstcall/                                             &
   model_basis_time, ancil_reftime, run_target_end,             &
   lclimrealyr, ltimer, model_analysis_mins, num_albcs,         &
   albc2_starttime_mins, lcal360, ltimers_user, l_nrun_as_crun, &
   lstashdumptimer

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLSTCALL_MOD'

CONTAINS

SUBROUTINE print_nlist_nlstcall()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
INTEGER :: i

CALL umPrint('Contents of namelist nlstcall', &
    src='nlstcall_mod')

DO i=1,6
  WRITE(lineBuffer,'(A,I0,A,I0)')'  model_basis_time(',i,') = ', &
       model_basis_time(i)
  CALL umPrint(lineBuffer,src='nlstcall_mod')
END DO
DO i=1,6
  WRITE(lineBuffer,'(A,I0,A,I0)')' ancil_reftime(',i,') = ', &
       ancil_reftime(i)
  CALL umPrint(lineBuffer,src='nlstcall_mod')
END DO
DO i=1,6
  WRITE(lineBuffer,'(A,I0,A,I0)')' run_target_end(',i,') = ', &
       run_target_end(i)
  CALL umPrint(lineBuffer,src='nlstcall_mod')
END DO
WRITE(lineBuffer,'(A,L1)')' lclimrealyr = ', lclimrealyr
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,L1)')' ltimer = ',ltimer 
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,I0)')' model_analysis_mins = ', model_analysis_mins
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,I0)')' num_albcs = ',  num_albcs     
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,I0)')' albc2_starttime_mins = ', albc2_starttime_mins
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,L1)')' lcal360 = ',lcal360
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,L1)')' ltimers_user = ',ltimers_user
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,L1)')' l_nrun_as_crun = ',l_nrun_as_crun
CALL umPrint(lineBuffer,src='nlstcall_mod')
WRITE(lineBuffer,'(A,L1)')' lstashdumptimer = ',lstashdumptimer
CALL umPrint(lineBuffer,src='nlstcall_mod')
CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='nlstcall_mod')

END SUBROUTINE print_nlist_nlstcall

SUBROUTINE read_nml_nlstcall(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCALL'
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 3 + 3 * 6
INTEGER, PARAMETER :: n_log = 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: model_basis_time(6)
  INTEGER :: ancil_reftime(6)
  INTEGER :: run_target_end(6)
  INTEGER :: model_analysis_mins
  INTEGER :: num_albcs
  INTEGER :: albc2_starttime_mins
  LOGICAL :: lclimrealyr
  LOGICAL :: ltimer
  LOGICAL :: lcal360
  LOGICAL :: ltimers_user
  LOGICAL :: l_nrun_as_crun
  LOGICAL :: lstashdumptimer
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,        &
                    n_log_in= n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=nlstcall, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLSTCALL", iomessage)

  my_nml % model_basis_time     = model_basis_time
  my_nml % ancil_reftime        = ancil_reftime
  my_nml % run_target_end       = run_target_end
  my_nml % model_analysis_mins  = model_analysis_mins
  my_nml % num_albcs            = num_albcs
  my_nml % albc2_starttime_mins = albc2_starttime_mins
  my_nml % lclimrealyr          = lclimrealyr
  my_nml % ltimer               = ltimer
  my_nml % lcal360              = lcal360
  my_nml % ltimers_user         = ltimers_user
  my_nml % l_nrun_as_crun       = l_nrun_as_crun
  my_nml % lstashdumptimer      = lstashdumptimer

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  model_basis_time     = my_nml % model_basis_time
  ancil_reftime        = my_nml % ancil_reftime
  run_target_end       = my_nml % run_target_end
  model_analysis_mins  = my_nml % model_analysis_mins
  num_albcs            = my_nml % num_albcs
  albc2_starttime_mins = my_nml % albc2_starttime_mins
  lclimrealyr          = my_nml % lclimrealyr
  ltimer               = my_nml % ltimer
  lcal360              = my_nml % lcal360
  ltimers_user         = my_nml % ltimers_user
  l_nrun_as_crun       = my_nml % l_nrun_as_crun
  lstashdumptimer      = my_nml % lstashdumptimer
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcall

END MODULE nlstcall_mod
