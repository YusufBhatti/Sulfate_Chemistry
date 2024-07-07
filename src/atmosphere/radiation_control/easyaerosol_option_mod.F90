! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold all EasyAerosol variables in easyaerosol 
!          namelist
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE easyaerosol_option_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength
USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE

! Declaration for EasyAerosol
!----------------------------
! Namelist items

! Use EasyAerosol in the SW and/or LW spectrum
LOGICAL :: l_easyaerosol_sw = .FALSE.
LOGICAL :: l_easyaerosol_lw = .FALSE.

! Use EasyAerosol Cloud Droplet Number Concentrations
! in the radiation scheme (for cloud albedo)
LOGICAL :: l_easyaerosol_cdnc = .FALSE.

! Use EasyAerosol Cloud Droplet Number Concentrations
! in the calculation of autoconversion rates.
LOGICAL :: l_easyaerosol_autoconv = .FALSE.

! EasyAerosol climatology can be read in as zonal means
LOGICAL :: l_easyaerosol_zonal     = .FALSE.

! Number of distributions in each spectrum: extinction, absorption,
! and asymmetry
INTEGER, PARAMETER :: n_easy = 3

! Information on netCDF files for EasyAerosol distributions
INTEGER, PARAMETER :: n_easyaerosol_files = n_easy*2+1 ! = 7 netCDF files
CHARACTER (LEN=filenamelength) :: easyaerosol_dir = 'unset'
                                                    ! Directory
CHARACTER (LEN=filenamelength) :: easyaerosol_files(n_easyaerosol_files) = (/  &
  'unset', 'unset', &
  'unset', 'unset', &
  'unset', 'unset', &
  'unset' /)
                                                    ! Name of EasyAerosol files

! Define the easyaerosol namelist 

NAMELIST /easyaerosol/ l_easyaerosol_sw, l_easyaerosol_lw, &
                       l_easyaerosol_cdnc, l_easyaerosol_autoconv, &
                       l_easyaerosol_zonal, &
                       easyaerosol_dir, easyaerosol_files

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EASYAEROSOL_OPTION_MOD'

CONTAINS

SUBROUTINE print_nlist_easyaerosol
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
INTEGER :: i ! loop counter
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_EASYAEROSOL'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist easyaerosol', &
    src='easyaerosol_option_mod')

WRITE(lineBuffer,'(A,L1)')' l_easyaerosol_sw = ',l_easyaerosol_sw
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_easyaerosol_lw = ',l_easyaerosol_lw
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_easyaerosol_cdnc = ',l_easyaerosol_cdnc
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_easyaerosol_autoconv = ',l_easyaerosol_autoconv
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,L1)')' l_easyaerosol_zonal = ',l_easyaerosol_zonal
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,A)')' easyaerosol_dir = ',TRIM(easyaerosol_dir)
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
WRITE(lineBuffer,'(A,I0,A,A)')' easyaerosol_files(',1,') = ',                  &
  easyaerosol_files(1)
CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
DO i = 2, n_easyaerosol_files
  IF (easyaerosol_files(i) /= 'unset') THEN
    WRITE(lineBuffer,'(A,I0,A,A)')' easyaerosol_files(',i,') = ',              &
      TRIM(easyaerosol_files(i))
    CALL umPrint(lineBuffer,src='easyaerosol_option_mod')
  END IF
END DO

CALL umPrint('- - - - - - end of namelist - - - - - -',                        &
             src='easyaerosol_option_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_easyaerosol

#if !defined(LFRIC)
SUBROUTINE read_nml_easyaerosol(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER (LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_EASYAEROSOL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_log = 5
INTEGER, PARAMETER :: n_chars = filenamelength * (1+ n_easyaerosol_files)

TYPE my_namelist
  SEQUENCE
  LOGICAL :: l_easyaerosol_sw
  LOGICAL :: l_easyaerosol_lw
  LOGICAL :: l_easyaerosol_cdnc
  LOGICAL :: l_easyaerosol_autoconv
  LOGICAL :: l_easyaerosol_zonal
  CHARACTER (LEN=filenamelength) :: easyaerosol_dir
  CHARACTER (LEN=filenamelength) :: easyaerosol_files(n_easyaerosol_files)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=easyaerosol, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist EASYAEROSOL", iomessage)
  
  my_nml % l_easyaerosol_sw       = l_easyaerosol_sw
  my_nml % l_easyaerosol_lw       = l_easyaerosol_lw
  my_nml % l_easyaerosol_cdnc     = l_easyaerosol_cdnc
  my_nml % l_easyaerosol_autoconv = l_easyaerosol_autoconv
  my_nml % l_easyaerosol_zonal    = l_easyaerosol_zonal
  ! end of logicals
  my_nml % easyaerosol_dir        = easyaerosol_dir
  my_nml % easyaerosol_files      = easyaerosol_files


END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  l_easyaerosol_sw       = my_nml % l_easyaerosol_sw
  l_easyaerosol_lw       = my_nml % l_easyaerosol_lw
  l_easyaerosol_cdnc     = my_nml % l_easyaerosol_cdnc
  l_easyaerosol_autoconv = my_nml % l_easyaerosol_autoconv
  l_easyaerosol_zonal    = my_nml % l_easyaerosol_zonal
  ! end of logicals
  easyaerosol_dir        = my_nml % easyaerosol_dir
  easyaerosol_files      = my_nml % easyaerosol_files

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_easyaerosol
#endif

END MODULE easyaerosol_option_mod
