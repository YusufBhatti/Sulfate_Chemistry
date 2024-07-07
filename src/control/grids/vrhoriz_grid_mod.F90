! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Horizontal grid for variable resolution

MODULE vrhoriz_grid_mod

USE Atmos_Max_Sizes
USE UM_ParParams
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE
!
! Description:
!   Contains the horizontal grid namelist for variable resolution, shared
!   between CreateBC and the reconfiguration
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UMDP3 v10.3 programming standards.
!
REAL :: lambda_input_p(row_length_max)      ! Longitude points on p-grid
REAL :: lambda_input_u(row_length_max)      ! Longitude points on u-grid
REAL :: phi_input_p(rows_max)               ! Latitude  points on p-grid
REAL :: phi_input_v(rows_max)               ! Latitude  points on v-grid

NAMELIST /horizgrid/                                                       &
  lambda_input_p, lambda_input_u,                                          &
  phi_input_p, phi_input_v

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VRHORIZ_GRID_MOD'

CONTAINS
SUBROUTINE print_nlist_horizgrid()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer = ''
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'PRINT_NLIST_HORIZGRID'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

CALL umPrint('Contents of namelist horizgrid', &
    src='vrhoriz_grid_mod')

! Print out first and last values, rather than the entire namelist
! (which is not useful and umprint will truncate anyway)
WRITE(lineBuffer,'(2(A,F16.10))') ' lambda_input_p = ',lambda_input_p(1),     &
                                  ' to ', last_element(lambda_input_p)
CALL umPrint(TRIM(lineBuffer),src='vrhoriz_grid_mod')

WRITE(lineBuffer,'(2(A,F16.10))') ' lambda_input_u = ',lambda_input_u(1),     &
                                  ' to ', last_element(lambda_input_u)
CALL umPrint(TRIM(lineBuffer),src='vrhoriz_grid_mod')

WRITE(lineBuffer,'(2(A,F16.10))') ' phi_input_p = ',phi_input_p(1),           &
                                  ' to ', last_element(phi_input_p)
CALL umPrint(TRIM(lineBuffer),src='vrhoriz_grid_mod')

WRITE(lineBuffer,'(2(A,F16.10))') ' phi_input_v = ',phi_input_v(1),            &
                                  ' to ', last_element(phi_input_v)
CALL umPrint(TRIM(lineBuffer),src='vrhoriz_grid_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='vrhoriz_grid_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

CONTAINS

REAL FUNCTION last_element(array)

! Returns the value of the last element that's not a missing data indicator
! of a real 1-d array

IMPLICIT NONE

REAL :: array(:)
INTEGER :: i

DO i = SIZE(array), 1, -1
  IF (array(i) /= rmdi) THEN
    last_element = array(i)
    RETURN
  END IF
END DO

last_element = rmdi
END FUNCTION last_element

END SUBROUTINE print_nlist_horizgrid


SUBROUTINE read_nml_horizgrid(unit_in)

USE check_iostat_mod,       ONLY: check_iostat
USE setup_namelist,         ONLY: setup_nml_type
USE um_parcore,             ONLY: mype

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'READ_NML_HORIZGRID'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_real = 2 * (row_length_max + rows_max)

TYPE my_namelist
  SEQUENCE
  REAL :: lambda_input_p(row_length_max)
  REAL :: lambda_input_u(row_length_max)
  REAL :: phi_input_p(rows_max)
  REAL :: phi_input_v(rows_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_real_in=n_real)

lambda_input_p(:)      = rmdi
lambda_input_u(:)      = rmdi
phi_input_p(:)         = rmdi
phi_input_v(:)         = rmdi

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=horizgrid, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist horizgrid", iomessage)

  my_nml % lambda_input_p = lambda_input_p
  my_nml % lambda_input_u = lambda_input_u
  my_nml % phi_input_p    = phi_input_p
  my_nml % phi_input_v    = phi_input_v

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  lambda_input_p  = my_nml % lambda_input_p
  lambda_input_u  = my_nml % lambda_input_u
  phi_input_p     = my_nml % phi_input_p
  phi_input_v     = my_nml % phi_input_v

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE read_nml_horizgrid

SUBROUTINE check_nml_horizgrid()

USE chk_opts_mod, ONLY:                chk_var, def_src

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_NML_HORIZGRID'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
def_src = RoutineName

! Are there any sensible checks that can be performed here ?
!CALL chk_var( lambda_input_p, 'lambda_input_p', '[0.0:9999.9]'  )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE check_nml_horizgrid

END MODULE vrhoriz_grid_mod
