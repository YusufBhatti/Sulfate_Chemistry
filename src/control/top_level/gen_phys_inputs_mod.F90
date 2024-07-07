! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing general physics runtime options
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 95

MODULE gen_phys_inputs_mod

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Methane oxidation
LOGICAL :: l_use_methox = .FALSE.

! Use mixing ratio in atmos_physics1 and atmos_physics2
LOGICAL :: l_mr_physics = .FALSE.

! Switch for leads temperature - 
! If FALSE, leads are set to TFS (freezing temperature of sea water).
! If TRUE, leads temperatures are prognostic values.
LOGICAL :: l_leads_temp_prog = .FALSE. 

! Switch for using the new qsat family
! If FALSE, use the old qsat routines
! If TRUE, use the new family, which fixes a small inconsistency in results
LOGICAL :: l_new_qsat = .FALSE.

NAMELIST /gen_phys_inputs/ l_use_methox, l_mr_physics, l_leads_temp_prog,     &
                           l_new_qsat

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GEN_PHYS_INPUTS_MOD'

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE print_nlist_gen_phys_inputs()

USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist gen_phys_inputs', &
     src='gen_phys_inputs_mod')

WRITE(lineBuffer,'(A,l1)') 'l_use_methox= ',l_use_methox
CALL umPrint(lineBuffer,src='gen_phys_inputs_mod')
WRITE(lineBuffer,'(A,l1)') 'l_mr_physics  = ',l_mr_physics
CALL umPrint(lineBuffer,src='gen_phys_inputs_mod')
WRITE(lineBuffer,'(A,l1)') 'l_leads_temp_prog  = ',l_leads_temp_prog
CALL umPrint(lineBuffer,src='gen_phys_inputs_mod')
WRITE(lineBuffer,'(A,l1)') 'l_new_qsat  = ',l_new_qsat
CALL umPrint(lineBuffer,src='gen_phys_inputs_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
     src='gen_phys_inputs_mod')

END SUBROUTINE print_nlist_gen_phys_inputs

!------------------------------------------------------------------------------
#if !defined(LFRIC)
SUBROUTINE read_nml_gen_phys_inputs(unit_in)

USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim
USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist,   ONLY: setup_nml_type

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_GEN_PHYS_INPUTS'

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_log = 4

TYPE my_namelist
  SEQUENCE
  LOGICAL :: l_use_methox
  LOGICAL :: l_mr_physics
  LOGICAL :: l_leads_temp_prog
  LOGICAL :: l_new_qsat
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log)

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=gen_phys_inputs, IOSTAT=ErrorStatus, &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist gen_phys_inputs", iomessage)
  my_nml % l_use_methox      = l_use_methox
  my_nml % l_mr_physics      = l_mr_physics
  my_nml % l_leads_temp_prog = l_leads_temp_prog
  my_nml % l_new_qsat        = l_new_qsat
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  l_use_methox      = my_nml % l_use_methox
  l_mr_physics      = my_nml % l_mr_physics
  l_leads_temp_prog = my_nml % l_leads_temp_prog
  l_new_qsat        = my_nml % l_new_qsat
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_gen_phys_inputs
#endif
!------------------------------------------------------------------------------

END MODULE gen_phys_inputs_mod
