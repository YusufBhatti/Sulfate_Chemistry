! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Control variables

!  Description:  Control variables for generic aspects of internal models.
!                Generic means values likely to be common to the control
!                of any sub-model/internal model.

! Migrated from include file cntlgen.h

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE nlstgen_mod

USE Submodel_Mod, ONLY: N_Internal_Model_Max
USE missing_data_mod, ONLY:  rmdi, imdi
USE filenamelength_mod, ONLY: filenamelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Max no. of irregular times for dumps
INTEGER, PARAMETER :: Dumptimes_Len1 = 40

! No. of time intervals for climate meaning
INTEGER, PARAMETER :: Meanfreq_Len1 = 4

INTEGER       :: steps_per_periodim(N_Internal_Model_Max) = imdi
INTEGER       :: secs_per_periodim(N_Internal_Model_Max)  = imdi

! Number of steps between atmosphere restart dumps
INTEGER       :: dumpfreqim(N_Internal_Model_Max) = 0

! If a Gregorian climate run is used then the following variable is used to
! control the dumping
INTEGER       :: greg_dump_freq = 0

! If Gregorian climate run is used then the following variable is chosen to
! select monthly dumps
LOGICAL       :: greg_monthly_dump = .FALSE.

! If Gregorian climate run is used then the following variable is chosen to
! produce a dump at the end
LOGICAL       :: end_of_run_dump = .FALSE.

! Archiving frequency  for atmos dumps
INTEGER       :: archdump_freqim(N_Internal_Model_Max) = imdi

! Timesteps (from start of run) at which restart dumps are written
INTEGER       :: dumptimesim(Dumptimes_Len1,N_Internal_Model_Max) = 0

! Indicators for mean dump frequency
INTEGER       :: meanfreqim(Meanfreq_Len1,N_Internal_Model_Max) = 0

! PP field selectors
INTEGER       :: ppselectim(Meanfreq_Len1,N_Internal_Model_Max) = 0

! Number of field headers to reserve for internal model mean  PPfiles
INTEGER       :: pp_len2_meanim(Meanfreq_Len1,N_Internal_Model_Max) = 30000

! Reference time for production of means
INTEGER       :: mean_reftimeim(6,N_Internal_Model_Max) = 0

! Unit reserved for mean PPs
INTEGER       :: ft_meanim(N_Internal_Model_Max) = 27

! Packing indicator for dumps
INTEGER       :: dump_packim(N_Internal_Model_Max) = imdi

! Select dumping and meaning option
INTEGER       :: i_dump_output = imdi
! 1 = No dumping or climate meaning
! 2 = Regular frequency dumps with possible meaning sequence
! 3 = Irregular dump times - no climate meaning possible
! Used by GUI to control whether dumptimesim or dumpfreqim are output
! to namelist.

! Units in which the dump frequency is specified in the GUI.  This
! is then used to convert the relevant arrays into timesteps in
! the UM code.
! 1 = Hours
! 2 = Days
! 3 = Timesteps
! These integers match the parameters defined in nlstcall_subs_mod
INTEGER       :: dump_frequency_units = imdi

INTEGER       :: ppxm = imdi ! Packing profile to use for climate mean files

LOGICAL   :: llboutim(N_Internal_Model_Max) = .FALSE.  ! Lateral b.c.'s

! Turn on climate meaning of output fields
! Used by GUI to control whether meanfreqim, mean_reftimeim and ppselectim
! are output to namelist.
LOGICAL   :: l_meaning_sequence = .FALSE.

CHARACTER(LEN=filenamelength) :: dump_filename_base = ""
CHARACTER(LEN=filenamelength) :: psum_filename_base = ""
CHARACTER(LEN=filenamelength) :: mean_1_filename_base = ""
CHARACTER(LEN=filenamelength) :: mean_2_filename_base = ""
CHARACTER(LEN=filenamelength) :: mean_3_filename_base = ""
CHARACTER(LEN=filenamelength) :: mean_4_filename_base = ""

CHARACTER(LEN=filenamelength) :: mean_filename_base(4) = ""

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

NAMELIST /nlstcgen/                                                   &
      steps_per_periodim, secs_per_periodim, dumpfreqim,              &
      dumptimesim, ppselectim,                                        &
      meanfreqim, mean_reftimeim,                                     &
      dump_packim, i_dump_output, l_meaning_sequence,                 &
      dump_frequency_units, greg_monthly_dump, end_of_run_dump,       &
      dump_filename_base, psum_filename_base,                         &
      mean_1_filename_base, mean_2_filename_base,                     &
      mean_3_filename_base, mean_4_filename_base, ppxm

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NLSTGEN_MOD'

CONTAINS

#if !defined(LFRIC)
SUBROUTINE read_nml_nlstcgen(unit_in)

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
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_NLSTCGEN'
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int =  3 + N_Internal_Model_Max *           &
                      (4 + 6 + Dumptimes_Len1 + 2 * Meanfreq_Len1) 
INTEGER, PARAMETER :: n_chars = 6*filenamelength
INTEGER, PARAMETER :: n_log = 3

TYPE my_namelist
  SEQUENCE
  INTEGER :: steps_per_periodim(N_Internal_Model_Max)
  INTEGER :: secs_per_periodim(N_Internal_Model_Max)
  INTEGER :: dumpfreqim(N_Internal_Model_Max)
  INTEGER :: dumptimesim(Dumptimes_Len1,N_Internal_Model_Max)
  INTEGER :: ppselectim(Meanfreq_Len1,N_Internal_Model_Max)
  INTEGER :: meanfreqim(Meanfreq_Len1,N_Internal_Model_Max)
  INTEGER :: mean_reftimeim(6,N_Internal_Model_Max)
  INTEGER :: dump_packim(N_Internal_Model_Max)
  INTEGER :: i_dump_output 
  INTEGER :: dump_frequency_units
  INTEGER :: ppxm
  LOGICAL :: l_meaning_sequence
  LOGICAL :: greg_monthly_dump
  LOGICAL :: end_of_run_dump
  CHARACTER(LEN=filenamelength) :: dump_filename_base
  CHARACTER(LEN=filenamelength) :: psum_filename_base
  CHARACTER(LEN=filenamelength) :: mean_1_filename_base
  CHARACTER(LEN=filenamelength) :: mean_2_filename_base
  CHARACTER(LEN=filenamelength) :: mean_3_filename_base
  CHARACTER(LEN=filenamelength) :: mean_4_filename_base
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=nlstcgen, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist NLSTCGEN", iomessage)

  my_nml % steps_per_periodim   = steps_per_periodim
  my_nml % secs_per_periodim    = secs_per_periodim
  my_nml % dumpfreqim           = dumpfreqim
  my_nml % dumptimesim          = dumptimesim
  my_nml % ppselectim           = ppselectim
  my_nml % meanfreqim           = meanfreqim
  my_nml % mean_reftimeim       = mean_reftimeim
  my_nml % dump_packim          = dump_packim
  my_nml % i_dump_output        = i_dump_output
  my_nml % dump_frequency_units = dump_frequency_units
  my_nml % ppxm                 = ppxm
  my_nml % l_meaning_sequence   = l_meaning_sequence
  my_nml % greg_monthly_dump    = greg_monthly_dump
  my_nml % end_of_run_dump      = end_of_run_dump
  my_nml % dump_filename_base   = dump_filename_base
  my_nml % psum_filename_base   = psum_filename_base
  my_nml % mean_1_filename_base = mean_1_filename_base
  my_nml % mean_2_filename_base = mean_2_filename_base
  my_nml % mean_3_filename_base = mean_3_filename_base
  my_nml % mean_4_filename_base = mean_4_filename_base
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  steps_per_periodim   = my_nml % steps_per_periodim
  secs_per_periodim    = my_nml % secs_per_periodim
  dumpfreqim           = my_nml % dumpfreqim
  dumptimesim          = my_nml % dumptimesim
  ppselectim           = my_nml % ppselectim
  meanfreqim           = my_nml % meanfreqim
  mean_reftimeim       = my_nml % mean_reftimeim
  dump_packim          = my_nml % dump_packim
  i_dump_output        = my_nml % i_dump_output
  dump_frequency_units = my_nml % dump_frequency_units
  ppxm                 = my_nml % ppxm
  l_meaning_sequence   = my_nml % l_meaning_sequence
  greg_monthly_dump    = my_nml % greg_monthly_dump
  end_of_run_dump      = my_nml % end_of_run_dump
  dump_filename_base   = my_nml % dump_filename_base
  psum_filename_base   = my_nml % psum_filename_base
  mean_1_filename_base = my_nml % mean_1_filename_base
  mean_2_filename_base = my_nml % mean_2_filename_base
  mean_3_filename_base = my_nml % mean_3_filename_base
  mean_4_filename_base = my_nml % mean_4_filename_base
END IF

CALL mpl_type_free(mpl_nml_type,icode)

mean_filename_base(1) = mean_1_filename_base
mean_filename_base(2) = mean_2_filename_base
mean_filename_base(3) = mean_3_filename_base
mean_filename_base(4) = mean_4_filename_base

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_nlstcgen
#endif

END MODULE nlstgen_mod
