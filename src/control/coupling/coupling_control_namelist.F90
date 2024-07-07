! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!  
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
! Code Description:
!   Language: FORTRAN 95
!
MODULE coupling_control_mod

USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

PRIVATE :: imdi

! OASIS coupling switch - set in um_shell/rcf_intialise based upon COUPLER 
! env_var
LOGICAL :: l_oasis = .FALSE.

!===========================================================================
! Component model run time switches 
!===========================================================================

LOGICAL :: l_senior = .FALSE.       ! Indicates this component is a main 
                                    ! UM atmosphere (senior) 

LOGICAL :: l_junior = .FALSE.       ! Indicates this component is a subordinate
                                    ! UM chemistry model controller (junior)

LOGICAL :: l_oasis_ocean = .FALSE.  ! Indicates this component is intended 
                                    ! to couple to an ocean component 

!===========================================================================
! LOGICAL Top Level items set from coupling_control namelist
!===========================================================================

! Ocean contains biogeochemistry
LOGICAL :: l_oasis_obgc = .FALSE.

! OASIS includes iceberg calving ancillary data
LOGICAL :: l_oasis_icecalve = .FALSE.

! Time Coupling Operations
LOGICAL :: l_oasis_timers = .FALSE.

!===========================================================================
! INTEGER Top Level items set from coupling_control namelist
!===========================================================================

! Coupling frequency in hours
INTEGER :: oasis_couple_freq_ao(2) = imdi
INTEGER :: oasis_couple_freq_oa(2) = imdi

INTEGER :: oasis_couple_freq_ac(2) = imdi
INTEGER :: oasis_couple_freq_ca(2) = imdi
INTEGER :: oasis_couple_freq_ac_stats(2) = imdi
INTEGER :: oasis_couple_freq_ca_stats(2) = imdi

INTEGER :: oasis_couple_freq_aw(2) = imdi
INTEGER :: oasis_couple_freq_wa(2) = imdi

!===========================================================================
! Define the coupling_control namelist
!===========================================================================

NAMELIST/coupling_control/  &
     oasis_couple_freq_ao, oasis_couple_freq_oa,             &
     oasis_couple_freq_ac, oasis_couple_freq_ca,             &
     oasis_couple_freq_ac_stats, oasis_couple_freq_ca_stats, &
     oasis_couple_freq_aw, oasis_couple_freq_wa,             &
     l_oasis_icecalve,                                       &
     l_oasis_obgc, l_oasis_timers,                           &
     l_oasis_ocean

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COUPLING_CONTROL_MOD'

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE check_nml_coupling_control()

! Description:
!   Subroutine to check variables based on options coupling_control namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CHECK_NML_COUPLING_CONTROL'
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check sensible values chosen
IF (l_oasis_ocean) THEN
  CALL chk_var(oasis_couple_freq_ao(1), 'oasis_couple_freq_ao(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_ao(2), 'oasis_couple_freq_ao(min)', '[0:60]' )

  CALL chk_var(oasis_couple_freq_oa(1), 'oasis_couple_freq_oa(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_oa(2), 'oasis_couple_freq_oa(min)', '[0:60]' )
END IF

IF (ANY(oasis_couple_freq_ac /= imdi).OR.       &
    ANY(oasis_couple_freq_ca /= imdi)) THEN
  CALL chk_var(oasis_couple_freq_ac(1), 'oasis_couple_freq_ac(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_ac(2), 'oasis_couple_freq_ac(min)', '[0:60]' )

  CALL chk_var(oasis_couple_freq_ca(1), 'oasis_couple_freq_ca(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_ca(2), 'oasis_couple_freq_ca(min)', '[0:60]' )
END IF

IF (ANY(oasis_couple_freq_ca_stats /= imdi).OR. &
    ANY(oasis_couple_freq_ac_stats /= imdi)) THEN
  CALL chk_var( oasis_couple_freq_ac_stats(1), &
               'oasis_couple_freq_ac_stats(hr)', '[>=0]' )
  CALL chk_var( oasis_couple_freq_ac_stats(2), &
               'oasis_couple_freq_ac_stats(min)', '[0:60]' )

  CALL chk_var( oasis_couple_freq_ca_stats(1), &
               'oasis_couple_freq_ca_stats(hr)', '[>=0]' )
  CALL chk_var( oasis_couple_freq_ca_stats(2), &
               'oasis_couple_freq_ca_stats(min)', '[0:60]' )
END IF

IF (ANY(oasis_couple_freq_aw /= imdi).OR. &
    ANY(oasis_couple_freq_aw /= imdi)) THEN
  CALL chk_var(oasis_couple_freq_aw(1), 'oasis_couple_freq_aw(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_aw(2), 'oasis_couple_freq_aw(min)','[0:60]' )
  CALL chk_var(oasis_couple_freq_wa(1), 'oasis_couple_freq_wa(hr)', '[>=0]' )
  CALL chk_var(oasis_couple_freq_wa(2), 'oasis_couple_freq_wa(min)','[0:60]' )
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_coupling_control

!------------------------------------------------------------------------------

SUBROUTINE print_nlist_coupling_control()

USE umPrintMgr, ONLY: umPrint
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_COUPLING_CONTROL'

! Dr Hook calls are retained for purely debugging purposes, since this routine
! has minimal impact on performance. 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist coupling_control', &
    src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_oasis = ',l_oasis
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_oasis_timers = ',l_oasis_timers
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_ao = ',  &
                                oasis_couple_freq_ao(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_oa = ',  & 
                                oasis_couple_freq_oa(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_ac = ',  & 
                                oasis_couple_freq_ac(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_ca = ',  &
                                oasis_couple_freq_ca(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_ac_stats = ',  & 
                                oasis_couple_freq_ac_stats(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_ca_stats = ',  &
                                oasis_couple_freq_ca_stats(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_aw = ',  &
                                oasis_couple_freq_aw(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,I0,I0)') ' oasis_couple_freq_wa = ',  &
                                oasis_couple_freq_wa(1:2)
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_oasis_ocean = ',l_oasis_ocean
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_oasis_obgc = ',l_oasis_obgc
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_oasis_icecalve = ',l_oasis_icecalve
CALL umPrint(lineBuffer,src='coupling_control_mod')
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_senior = ',l_senior
CALL umPrint(lineBuffer,src='coupling_control_mod')
WRITE(lineBuffer,'(A,L1)') ' l_junior = ',l_junior

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_coupling_control

!------------------------------------------------------------------------------

SUBROUTINE read_nml_coupling_control(unit_in)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_COUPLING_CONTROL'
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 16
INTEGER, PARAMETER :: n_log = 4

TYPE my_namelist
  SEQUENCE
  INTEGER :: oasis_couple_freq_ao(2)
  INTEGER :: oasis_couple_freq_oa(2)
  INTEGER :: oasis_couple_freq_ac(2)
  INTEGER :: oasis_couple_freq_ca(2)
  INTEGER :: oasis_couple_freq_ac_stats(2)
  INTEGER :: oasis_couple_freq_ca_stats(2)
  INTEGER :: oasis_couple_freq_aw(2)
  INTEGER :: oasis_couple_freq_wa(2)
  LOGICAL :: l_oasis_timers
  LOGICAL :: l_oasis_ocean
  LOGICAL :: l_oasis_obgc
  LOGICAL :: l_oasis_icecalve
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! Dr Hook calls are retained for purely debugging purposes, since this routine
! has minimal impact on performance. 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_log_in=n_log)

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=coupling_control, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist UM_COUPLING_CONTROL", iomessage)

  my_nml % oasis_couple_freq_ao       = oasis_couple_freq_ao
  my_nml % oasis_couple_freq_oa       = oasis_couple_freq_oa
  my_nml % oasis_couple_freq_ac       = oasis_couple_freq_ac
  my_nml % oasis_couple_freq_ca       = oasis_couple_freq_ca
  my_nml % oasis_couple_freq_ac_stats = oasis_couple_freq_ac_stats
  my_nml % oasis_couple_freq_ca_stats = oasis_couple_freq_ca_stats
  my_nml % oasis_couple_freq_aw       = oasis_couple_freq_aw
  my_nml % oasis_couple_freq_wa       = oasis_couple_freq_wa
  my_nml % l_oasis_timers    = l_oasis_timers
  my_nml % l_oasis_ocean     = l_oasis_ocean
  my_nml % l_oasis_obgc      = l_oasis_obgc
  my_nml % l_oasis_icecalve  = l_oasis_icecalve
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  oasis_couple_freq_ao       = my_nml % oasis_couple_freq_ao
  oasis_couple_freq_oa       = my_nml % oasis_couple_freq_oa
  oasis_couple_freq_ac       = my_nml % oasis_couple_freq_ac
  oasis_couple_freq_ca       = my_nml % oasis_couple_freq_ca
  oasis_couple_freq_ac_stats = my_nml % oasis_couple_freq_ac_stats
  oasis_couple_freq_ca_stats = my_nml % oasis_couple_freq_ca_stats
  oasis_couple_freq_aw       = my_nml % oasis_couple_freq_aw
  oasis_couple_freq_wa       = my_nml % oasis_couple_freq_wa
  l_oasis_timers    = my_nml % l_oasis_timers
  l_oasis_ocean     = my_nml % l_oasis_ocean
  l_oasis_obgc      = my_nml % l_oasis_obgc
  l_oasis_icecalve  = my_nml % l_oasis_icecalve
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_coupling_control

!------------------------------------------------------------------------------

END MODULE coupling_control_mod
