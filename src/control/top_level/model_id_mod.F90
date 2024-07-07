! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

! Module for the model configuration identifier
! this can be used to explicitly set LBEXP in UMDP F3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

MODULE model_id_mod

USE missing_data_mod, ONLY: imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, ONLY: newline

IMPLICIT NONE

INTEGER   ::    itab = imdi ! model configuration identifier
                            ! if set then LBEXP will be set to this.

NAMELIST/configid/itab

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MODEL_ID_MOD'

CONTAINS

SUBROUTINE read_nml_configid(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_CONFIGID'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 1

TYPE my_namelist
  SEQUENCE
  INTEGER :: itab
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=configid, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist configid", iomessage)

  my_nml % itab = itab

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  itab = my_nml % itab

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_configid


SUBROUTINE check_configid()

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=30), PARAMETER :: routinename = 'check_configid'
INTEGER :: icode
    
IF (itab == imdi) THEN 
  icode    = -10
  cmessage =                                                          newline//&
  'itab is set to missing data, pp_head(LBEXP) will be set to zero'
  CALL ereport( routinename, icode, cmessage )
END IF

END SUBROUTINE check_configid

END MODULE model_id_mod
