! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the murk scheme
!
! Method:
!   Switches and associated data values used by the murk scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Atmosphere Service

MODULE murk_inputs_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE missing_data_mod, ONLY: rmdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!===========================================================================
! LOGICAL options set from run_murk namelist
!===========================================================================

LOGICAL :: l_murk        = .FALSE.  ! Use total aerosol field

LOGICAL :: l_murk_advect = .FALSE.  ! Aerosol advection

LOGICAL :: l_murk_source = .FALSE.  ! Aerosol source & sink terms

LOGICAL :: l_murk_lbc    = .FALSE.  ! Murk aerosol lbcs active

LOGICAL :: l_murk_vis    = .FALSE.  ! Murk aerosol used for visibility

!===========================================================================
! REAL options set from run_murk namelist
!===========================================================================

REAL :: murk_source_scale   = rmdi     ! Aerosol emissions scaling factor

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Define the run_murk namelist

NAMELIST /run_murk/ l_murk, l_murk_advect, l_murk_source, l_murk_lbc,   &
                    l_murk_vis, murk_source_scale

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================

LOGICAL :: l_murk_bdry   = .FALSE.  ! UK Mes boundary model

LOGICAL :: l_murk_rad    = .FALSE.  ! Include radiative effects of aerosol

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MURK_INPUTS_MOD'

CONTAINS

SUBROUTINE print_nlist_run_murk()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_MURK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_murk', &
    src='murk_inputs_mod')

WRITE(lineBuffer,'(A,L1)')' l_murk = ', l_murk
CALL umPrint(lineBuffer,src='murk_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_murk_advect  = ', l_murk_advect
CALL umPrint(lineBuffer,src='murk_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_murk_source = ', l_murk_source
CALL umPrint(lineBuffer,src='murk_inputs_mod')
WRITE(lineBuffer,'(A,F16.4)')' murk_source_scale = ',murk_source_scale
CALL umPrint(lineBuffer,src='murk_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_murk_lbc = ', l_murk_lbc
CALL umPrint(lineBuffer,src='murk_inputs_mod')
WRITE(lineBuffer,'(A,L1)')' l_murk_vis = ', l_murk_vis
CALL umPrint(lineBuffer,src='murk_inputs_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='murk_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_murk

#if !defined(LFRIC)
SUBROUTINE read_nml_run_murk(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_MURK'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_real = 1
INTEGER, PARAMETER :: n_log = 5

TYPE my_namelist
  SEQUENCE
  REAL    :: murk_source_scale
  LOGICAL :: l_murk
  LOGICAL :: l_murk_advect
  LOGICAL :: l_murk_source
  LOGICAL :: l_murk_lbc
  LOGICAL :: l_murk_vis
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in=n_real,           &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Murk, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Murk", iomessage)

  my_nml % murk_source_scale = murk_source_scale
  ! end of reals
  my_nml % l_murk        = l_murk
  my_nml % l_murk_advect = l_murk_advect
  my_nml % l_murk_source = l_murk_source
  my_nml % l_murk_lbc    = l_murk_lbc
  my_nml % l_murk_vis    = l_murk_vis
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  murk_source_scale = my_nml % murk_source_scale 
  ! end of reals
  l_murk         =  my_nml % l_murk
  l_murk_advect  =  my_nml % l_murk_advect
  l_murk_source  =  my_nml % l_murk_source
  l_murk_lbc     =  my_nml % l_murk_lbc
  l_murk_vis     =  my_nml % l_murk_vis

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_murk
#endif

END MODULE murk_inputs_mod
