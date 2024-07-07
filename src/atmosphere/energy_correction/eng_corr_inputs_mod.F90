! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the energy correction scheme
!
! Method:
!   Switches and associated data values used by the energy correction scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Energy Correction

MODULE eng_corr_inputs_mod

USE missing_data_mod, ONLY: imdi

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!===========================================================================
! LOGICAL options set from run_eng_corr namelist
!===========================================================================

LOGICAL :: l_emcorr    = .FALSE.  ! Turns on energy correction code

LOGICAL :: lmass_corr  = .FALSE.  ! Apply mass correction

LOGICAL :: lemq_print  = .FALSE.  ! Print additional info from em code

!===========================================================================
! INTEGER values set from run_eng_corr namelist
!===========================================================================

INTEGER :: a_energyhours = imdi   ! Number of hours per energy
                                  ! correction period.

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Define the run_eng_corr namelist

NAMELIST /run_eng_corr/ l_emcorr, lmass_corr, lemq_print, a_energyhours

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================

LOGICAL :: lqt_corr    = .FALSE.  ! Apply total moisture correction

LOGICAL :: lenergy     = .FALSE.  ! Timestep to cal energy correction

LOGICAL :: lflux_reset = .FALSE.  ! Timestep to reset flux array in D1

!===========================================================================
! INTEGER values not set in namelist
!===========================================================================

INTEGER :: a_energysteps = imdi   ! Number of model timesteps per energy
                                  ! correction period 
                                  ! (calculated from a_energyhours)

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ENG_CORR_INPUTS_MOD'

CONTAINS


SUBROUTINE print_nlist_run_eng_corr()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_ENG_CORR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_eng_corr', &
    src='eng_corr_inputs_mod')

WRITE(lineBuffer,*)'l_emcorr = ',l_emcorr
CALL umPrint(lineBuffer,src='eng_corr_inputs_mod')
WRITE(lineBuffer,*)'lmass_corr = ',lmass_corr
CALL umPrint(lineBuffer,src='eng_corr_inputs_mod')
WRITE(lineBuffer,*)'lemq_print = ',lemq_print
CALL umPrint(lineBuffer,src='eng_corr_inputs_mod')
WRITE(lineBuffer,*)'a_energyhours = ',a_energyhours
CALL umPrint(lineBuffer,src='eng_corr_inputs_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='eng_corr_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_eng_corr

#if !defined(LFRIC)
SUBROUTINE read_nml_run_eng_corr(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_ENG_CORR'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_log = 3

TYPE my_namelist
  SEQUENCE
  INTEGER :: a_energyhours
  LOGICAL :: l_emcorr
  LOGICAL :: lmass_corr
  LOGICAL :: lemq_print
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Eng_Corr, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Eng_Corr", iomessage)

  my_nml % a_energyhours = a_energyhours
  my_nml % l_emcorr      = l_emcorr
  my_nml % lmass_corr    = lmass_corr
  my_nml % lemq_print    = lemq_print

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  a_energyhours = my_nml % a_energyhours
  l_emcorr      = my_nml % l_emcorr
  lmass_corr    = my_nml % lmass_corr
  lemq_print    = my_nml % lemq_print

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_eng_corr
#endif

SUBROUTINE check_run_eng_corr()

! Description:
!   Subroutine to set control variables based on the
!   options selected in the run_eng_corr.

USE conversions_mod,       ONLY: isec_per_hour
!   Number of seconds per hour
USE nlstgen_mod,            ONLY: secs_per_periodim, steps_per_periodim
!   Number of seconds per period, typically 1 day (86400 sec)
!   and the number of time steps per period
!
!   The "im" in these variables denotes that these are arrays
!   with one entry for each different internal model (or sub-model)
USE Submodel_Mod,           ONLY: atmos_im
!   The "im" array entry for the atmosphere internal model

!   Printing and error reporting:
USE um_parcore,             ONLY: mype
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: newline

IMPLICIT NONE

INTEGER :: errorstatus
CHARACTER(LEN=errormessagelength) :: cmessage

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_ENG_CORR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate number of timesteps per energy correction step from number
! of hours per energy correction step
a_energysteps = a_energyhours * isec_per_hour                                  &
                 * steps_per_periodim(atmos_im) / secs_per_periodim(atmos_im)

! Also check that the model timestep divides exactly into this period.
! We do this by inverting the above calculation to see if one recovers
! the same value of the input variable "a_energyhours".
! This should usually be the case (as the UM timestep will almost
! always divide into the hour) but it is worth trapping this as an error

IF (a_energyhours /= ( (a_energysteps * secs_per_periodim(atmos_im))           &
                      / (steps_per_periodim(atmos_im) * isec_per_hour) ) ) THEN
  
  IF (mype == 0) THEN
    
    errorstatus = 1
    WRITE(cmessage,FMT='(A,I4,A)')                                             &
         'Time period for energy correction steps is ', a_energyhours,  &
         ' hours.' // newline //                                               &
         'This is not a multiple of the model timestep.'
         
    CALL ereport('check_run_eng_corr', errorstatus, cmessage)
  END IF ! mype == 0

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE check_run_eng_corr

END MODULE eng_corr_inputs_mod
