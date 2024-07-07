#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
MODULE OASIS3_atmos_init_mod
  !
  ! Description:
  !     Define modules from OASIS3-MCT system containing variables 
  !     and routines which we want to use.
  !
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !=====================================================================

USE NetCDF

USE mod_prism, ONLY: prism_init_comp_proto &
                    ,prism_get_localcomm_proto

USE oasis_atm_data_mod
USE coupling_control_mod, ONLY: l_junior

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE UM_ParVars
USE um_types
USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CONTAINS
SUBROUTINE OASIS3_atmos_init(comm64,ierror,cmessage)

IMPLICIT NONE
!======================================================================
! Description: Initialise OASIS3/OASIS3-MCT for UM model components
!======================================================================

INTEGER (KIND=integer64), INTENT(OUT) :: ierror
INTEGER (KIND=integer64), INTENT(OUT) :: comm64
CHARACTER(LEN=errormessagelength),INTENT(OUT) :: cmessage

! Local 32 bit integers for oasis interfacing
INTEGER (KIND=integer32)              :: ierr32
INTEGER (KIND=integer32)              :: comm32
CHARACTER (LEN=6)                     :: cpl_component ! Component name

!
! Initialise Component
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!

!   DrHook call commented out: comms not yet initialised.
!   IF (lhook) CALL dr_hook('ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr32=0
IF (l_junior) THEN
  cpl_component = 'junior'
ELSE
  cpl_component = 'toyatm'
END IF

CALL prism_init_comp_proto(OASIS_comp_id, cpl_component, ierr32)
ierror=ierr32

IF (ierror /= 0) THEN
  cmessage = cpl_component // ": PRISM_init_comp_proto Error"
ELSE
  cmessage = cpl_component // ": PRISM_init_comp_proto OK"
END IF

! Get the communicator allocated to this model component.
! For MPI2 OASIS options, this just returns MPI_COMM_WORLD (because
! under those circumstances the MPI2 spawn approach is used), but for
! the MPI1 case we must use this to override the UM's own internal
! communicator
CALL PRISM_get_localcomm_proto(comm32,ierr32 )
ierror=ierr32
comm64=comm32

!   DrHook call commented out: comms not yet initialised.
!   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE OASIS3_atmos_init

END MODULE OASIS3_atmos_init_mod
#endif
