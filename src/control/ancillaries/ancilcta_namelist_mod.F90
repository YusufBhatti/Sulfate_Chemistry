! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Ancil-related logicals and values read in by namelist
!
MODULE ancilcta_namelist_mod

! Description:
!   Module containing the ANCILCTA namelist and associated declarations.
!
! Method:
!  Contents of the ANCILCTA namelist are read in and stored here.
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Ancillaries

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: nancil_lookupsa  ! Max no of fields to be read

LOGICAL :: l_sstanom = .FALSE. ! Indicator if SST anom to be formed
                               !  (RECON=T) or used (-DEF,RECON)

LOGICAL :: l_amipii_ice_processing   = .FALSE. 
! True if special AMIP II updating

LOGICAL :: use_lookup_dates_anc_time_interp   = .FALSE. 
! True if using date information from ancil field pp lookup rather 
! than from the fixed header of ancil file

NAMELIST/ancilcta/                                                 &
  l_sstanom, use_lookup_dates_anc_time_interp,                     &
  l_amipii_ice_processing, nancil_lookupsa

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ANCILCTA_NAMELIST_MOD'

CONTAINS

SUBROUTINE print_nlist_ancilcta()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)      :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_ANCILCTA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist ancilcta', &
    src='ancilcta_namelist_mod')

WRITE(lineBuffer,'(A,L1)')' l_sstanom = ',l_sstanom
CALL umPrint(lineBuffer,src='ancilcta_namelist_mod')
WRITE(lineBuffer,'(A,L1)')' l_amipii_ice_processing = ',l_amipii_ice_processing
CALL umPrint(lineBuffer,src='ancilcta_namelist_mod')
WRITE(lineBuffer,'(A,L1)')' use_lookup_dates_anc_time_interp = ', &
     use_lookup_dates_anc_time_interp
CALL umPrint(lineBuffer,src='ancilcta_namelist_mod')
WRITE(lineBuffer,'(A,I0)')' nancil_lookupsa = ',nancil_lookupsa
CALL umPrint(lineBuffer,src='ancilcta_namelist_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='ancilcta_namelist_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_ancilcta

SUBROUTINE read_nml_ancilcta(unit_in)

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
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_ANCILCTA'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_log = 3
INTEGER, PARAMETER :: n_int = 1

TYPE my_namelist
  SEQUENCE
  LOGICAL :: l_sstanom
  LOGICAL :: l_amipii_ice_processing
  LOGICAL :: use_lookup_dates_anc_time_interp
  INTEGER :: nancil_lookupsa
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=ancilcta, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist ANCILCTA", iomessage)

  my_nml % l_sstanom                = l_sstanom
  my_nml % use_lookup_dates_anc_time_interp = use_lookup_dates_anc_time_interp
  my_nml % l_amipii_ice_processing   = l_amipii_ice_processing
  my_nml % nancil_lookupsa          = nancil_lookupsa

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  l_sstanom                = my_nml % l_sstanom
  use_lookup_dates_anc_time_interp = my_nml % use_lookup_dates_anc_time_interp
  l_amipii_ice_processing   = my_nml % l_amipii_ice_processing
  nancil_lookupsa          = my_nml % nancil_lookupsa

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_ancilcta

SUBROUTINE check_nml_ancilcta()

! Description:
!   Subroutine to check control variable in ancilcta namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE parkind1, ONLY: jprb,jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER (LEN=*),  PARAMETER :: RoutineName = 'CHECK_NML_ANCILCTA'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

CALL chk_var( nancil_lookupsa,    'nancil_lookupsa',    '[1:50000]' )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_ancilcta

END MODULE ancilcta_namelist_mod
