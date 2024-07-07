! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Purpose: History data needed by top level routines
!          and passed from run to run.

MODULE history

USE Submodel_Mod, ONLY: n_internal_model_max
USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Flag to indicate if history file has been found for this run
LOGICAL, SAVE :: history_file_exists = .FALSE.

! Variable to store the history unit if a history file is used
INTEGER, SAVE :: history_unit = imdi
  
! Integer history variables

! Array containing model data time (Set from dump data time,
! or same as MODEL_ANALYSIS_HRS after assimilation)
INTEGER :: model_data_time(6) = imdi
! Array to store original actual model basis time as opposed to the model
! basis time set on input which may be set to 0. The first run sets it for
! all runs, but only continuation runs for which the input model basis time
! is zero read its value.
INTEGER :: original_basis_time(6) = imdi


NAMELIST /nlihisto/ model_data_time, original_basis_time


! Integer History variables for generic aspects of internal models
! Generic means values likely to be common to the control of any
! sub-model/internal model.

! History block copy of A_STEP held in model_time_mod
INTEGER :: h_stepim(n_internal_model_max) = 0
! No of means activated
INTEGER :: mean_offsetim(n_internal_model_max) = imdi
! Offset between MEAN_REFTIME and model basis time(in model dumps)
INTEGER :: offset_dumpsim(n_internal_model_max) = imdi
! No of mean periods chosen
INTEGER :: mean_numberim(n_internal_model_max) = imdi
! Indicators used to correct logical units are used for
! atmos partial sum dump I/O
INTEGER :: run_meanctl_indicim(4,n_internal_model_max) = 1

NAMELIST /nlihistg/ h_stepim, mean_offsetim, offset_dumpsim, &
   run_meanctl_indicim

! Character History variables for the model overall.
CHARACTER(LEN=8)  ::  run_type = 'Reconfig'       ! Type of run

LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

NAMELIST /nlchisto/ run_type

! Character History variables for managing dump names
! Blank name
CHARACTER(LEN=15), PARAMETER :: blank_file_name = 'blank_file_name'
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max) &
    = blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HISTORY'

CONTAINS

SUBROUTINE read_nml_nlihisto(unithist)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unithist
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLIHISTO'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 6 + 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: model_data_time(6)
  INTEGER :: original_basis_time(6)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  READ(UNIT=unithist, NML=nlihisto, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLIHISTO", iomessage)

  my_nml % model_data_time = model_data_time
  my_nml % original_basis_time = original_basis_time

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  model_data_time = my_nml % model_data_time
  original_basis_time = my_nml % original_basis_time

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlihisto

SUBROUTINE read_nml_nlchisto(unithist)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unithist
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLCHISTO'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_chars = 8

TYPE my_namelist
  SEQUENCE
  CHARACTER(LEN=8)  ::  run_type
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unithist, NML=nlchisto, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLCHISTO", iomessage)

  my_nml % run_type  = run_type

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  run_type  = my_nml % run_type

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlchisto

SUBROUTINE read_nml_nlihistg(unithist)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unithist
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLIHISTG'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 7 * n_internal_model_max

TYPE my_namelist
  SEQUENCE
  INTEGER :: h_stepim(n_internal_model_max)
  INTEGER :: mean_offsetim(n_internal_model_max)
  INTEGER :: offset_dumpsim(n_internal_model_max)
  INTEGER :: run_meanctl_indicim(4,n_internal_model_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  READ(UNIT=unithist, NML=nlihistg, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLIHISTG", iomessage)

  my_nml % h_stepim            = h_stepim
  my_nml % mean_offsetim       = mean_offsetim
  my_nml % offset_dumpsim      = offset_dumpsim
  my_nml % run_meanctl_indicim = run_meanctl_indicim

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  h_stepim            = my_nml % h_stepim
  mean_offsetim       = my_nml % mean_offsetim
  offset_dumpsim      = my_nml % offset_dumpsim
  run_meanctl_indicim = my_nml % run_meanctl_indicim

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlihistg

SUBROUTINE read_nml_nlchistg(unithist)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unithist
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLCHISTG'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_chars =  256 * n_internal_model_max

TYPE my_namelist
  SEQUENCE
  CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unithist, NML=nlchistg, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLCHISTG", iomessage)

  my_nml % checkpoint_dump_im = checkpoint_dump_im

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  checkpoint_dump_im = my_nml % checkpoint_dump_im

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlchistg

END MODULE history
