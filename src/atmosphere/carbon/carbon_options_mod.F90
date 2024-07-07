! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
!
!  Data module for switches for UM inputs,carbon options namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Carbon

MODULE carbon_options_mod

USE missing_data_mod, ONLY: imdi, rmdi
IMPLICIT NONE

! Carbon cycle
! Include surface emissions
LOGICAL :: l_co2_emits        = .FALSE.

! Specification of CO2 absorption
! 1 => Simple method with fixed value.
! 2 => Complex method allowing linear and/or exponential variation.
! 3 => From the interactive carbon cycle.
INTEGER :: i_co2_opt = imdi

NAMELIST /carbon_options/  i_co2_opt, l_co2_emits

!===========================================================================
! items not set in namelist
!===========================================================================
! Carbon cycle
! Interactive 3D CO2 field for use with carbon cycle model
! This is set by the check_carbon_options routine based on the value of
! i_co2_opt
LOGICAL :: l_co2_interactive  = .FALSE.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CARBON_OPTIONS_MOD'

CONTAINS

SUBROUTINE print_nlist_carbon_options()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist carbon_options', &
    src='carbon_options_mod')

WRITE(lineBuffer,'(A,I0)') ' i_co2_opt = ',i_co2_opt
CALL umPrint(lineBuffer,src='carbon_options_mod')
WRITE(lineBuffer,'(A,L1)') ' l_co2_emits = ',l_co2_emits
CALL umPrint(lineBuffer,src='carbon_options_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='carbon_options_mod')

END SUBROUTINE print_nlist_carbon_options

!-----------------------------------------------------------------------

SUBROUTINE check_carbon_options()
! Description:
!   Subroutine to apply logic controls and set control variables based on the
!   options selected in the carbon_options namelist.
IMPLICIT NONE
! Set logicals based on integer choice in namelist
IF (i_co2_opt == 3) THEN
  l_co2_interactive = .TRUE.
END IF
RETURN
END SUBROUTINE check_carbon_options

!-----------------------------------------------------------------------

#if !defined(LFRIC)
SUBROUTINE read_nml_carbon_options(unit_in)

USE um_parcore, ONLY: mype
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength
USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_CARBON_OPTIONS'

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_log = 1

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_co2_opt
  LOGICAL :: l_co2_emits
END TYPE my_namelist
TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)
CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, n_log_in=n_log)

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=carbon_options, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist CARBON_OPTIONS", iomessage)
  my_nml % i_co2_opt         = i_co2_opt
  my_nml % l_co2_emits       = l_co2_emits
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  i_co2_opt         = my_nml % i_co2_opt
  l_co2_emits       = my_nml % l_co2_emits
END IF

CALL mpl_type_free(mpl_nml_type,icode)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_nml_carbon_options
#endif

!-----------------------------------------------------------------------

END MODULE carbon_options_mod
