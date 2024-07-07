! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for switches for inputs to configure LAMs in the
! lam_config namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
MODULE lam_config_inputs_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

REAL ::  delta_lon = rmdi
! Column Spacing

REAL ::  delta_lat = rmdi
! Row Spacing

REAL ::  frstlata = rmdi
! First Latitude

REAL ::  frstlona = rmdi
!  First Longitude

REAL ::  polelata = rmdi
! North Pole Latitude

REAL ::  polelona  = rmdi
! North Pole Longitude

INTEGER :: n_rims_to_do  = 1
! rim size for LAM

NAMELIST /lam_config/ delta_lon, delta_lat, frstlata, frstlona, &
     polelata, polelona, n_rims_to_do

!-----------------------------------------------------
! Related, non-namelist variables
!-----------------------------------------------------

! Variables for storing LAM pseudo-grid co-ordinates
! when interpolating between lat-long/Cartesian grids:
REAL ::  target_delta_lon = rmdi ! Column Spacing 

REAL ::  target_delta_lat = rmdi ! Row Spacing 

REAL ::  target_frstlata = rmdi  ! First Latitude 

REAL ::  target_frstlona = rmdi  ! First Longitude

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LAM_CONFIG_INPUTS_MOD'

CONTAINS

SUBROUTINE check_nml_lam_config()
! Description:
!   Subroutine to check variables based on options lam_config namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE ereport_mod,  ONLY: ereport
USE model_domain_mod, ONLY: model_type, mt_global, l_cartesian

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CHECK_NML_LAM_CONFIG'
REAL(KIND=jprb)              :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check sensible values chosen for LAMs
! These variables set to rmdi for Global runs
IF (model_type /= mt_global) THEN

  IF (.NOT. l_cartesian) THEN  ! Don't check lat-lon for Cartesian grid case
                               ! Cartesian uses info from recon_idealised
    CALL chk_var( delta_lat, 'delta_lat', '[0.0:100.0]' )
    CALL chk_var( delta_lon, 'delta_lon', '[0.0:100.0]' )
    CALL chk_var( frstlata, 'frstlata', '[-90.0:90.0]' )
    CALL chk_var( frstlona, 'frstlona', '[0.0:360.0]' )
    CALL chk_var( polelata, 'polelata', '[-90.0:180.0]' )
    CALL chk_var( polelona, 'polelona', '[0.0:360.0]' )
  END IF

  ! Need to check whether Lat-Lon or Cartesian
  CALL chk_var( n_rims_to_do, 'n_rims_to_do', '[>=1]' )
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_lam_config

SUBROUTINE print_nlist_lam_config()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_LAM_CONFIG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist lam_config', &
     src='lam_config_inputs_mod')

WRITE(lineBuffer,'(A,F0.6)') ' delta_lon = ',delta_lon
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)') ' delta_lat = ',delta_lat
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)') ' frstlata = ',frstlata
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)') ' frstlona = ',frstlona
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)') ' polelata = ',polelata
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,F0.6)') ' polelona = ',polelona
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
WRITE(lineBuffer,'(A,I0)') ' n_rims_to_do = ',n_rims_to_do
CALL umPrint(lineBuffer,src='lam_config_inputs_mod')
CALL umPrint('- - - - - - end of namelist - - - - - -', &
     src='lam_config_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_lam_config

SUBROUTINE read_nml_lam_config(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_LAM_CONFIG'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_real = 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: n_rims_to_do
  REAL :: delta_lon
  REAL :: delta_lat
  REAL :: frstlata
  REAL :: frstlona
  REAL :: polelata
  REAL :: polelona
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
                    n_real_in=n_real)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=lam_config, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist lam_config", iomessage)

  my_nml % n_rims_to_do = n_rims_to_do
  my_nml % delta_lon    = delta_lon
  my_nml % delta_lat    = delta_lat
  my_nml % frstlata     = frstlata
  my_nml % frstlona     = frstlona
  my_nml % polelata     = polelata
  my_nml % polelona     = polelona

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  n_rims_to_do = my_nml % n_rims_to_do
  delta_lon    = my_nml % delta_lon
  delta_lat    = my_nml % delta_lat
  frstlata     = my_nml % frstlata
  frstlona     = my_nml % frstlona
  polelata     = my_nml % polelata
  polelona     = my_nml % polelona

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_lam_config

END MODULE lam_config_inputs_mod
