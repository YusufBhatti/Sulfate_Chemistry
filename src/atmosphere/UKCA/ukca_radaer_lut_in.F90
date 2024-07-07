! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Reads in a UKCA_RADAER look-up table (namelist format) and copy
!  its contents in the structure given in argument.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_lut_read_in

USE umprintmgr,             ONLY: newline
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,            ONLY: ereport
USE parkind1,               ONLY: jpim, jprb
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

!
! Local copy of the contents of a UKCA look-up table.
! Here, fixed array sizes are needed.
! Note that array indexing starts at 0 in the namelist.
INTEGER, PARAMETER :: npd_x  = 50
INTEGER, PARAMETER :: npd_nr = 50
INTEGER, PARAMETER :: npd_ni = 50

REAL :: stdev
INTEGER :: n_x
INTEGER :: n_nr
INTEGER :: n_ni
INTEGER :: funit
REAL :: x_min
REAL :: x_max
REAL :: nr_min
REAL :: nr_max
REAL :: ni_min
REAL :: ni_max
REAL :: ni_c
REAL :: ukca_absorption(0:npd_x,0:npd_ni,0:npd_nr)
REAL :: ukca_scattering(0:npd_x,0:npd_ni,0:npd_nr)
REAL :: ukca_asymmetry(0:npd_x,0:npd_ni,0:npd_nr)
REAL :: volume_fraction(0:npd_x)

NAMELIST /ukcanml/                                                      &
          stdev,                                                        &
          n_x, n_nr, n_ni,                                              &
          x_min, x_max,                                                 &
          nr_min, nr_max,                                               &
          ni_min, ni_max, ni_c,                                         &
          ukca_absorption, ukca_scattering, ukca_asymmetry,             &
          volume_fraction

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_LUT_READ_IN'

! Overload interface
INTERFACE ukca_radaer_get_lut_index
  MODULE PROCEDURE ukca_radaer_get_lut_index_scalar
  MODULE PROCEDURE ukca_radaer_get_lut_index_array
END INTERFACE

CONTAINS

SUBROUTINE ukca_radaer_lut_in(icode, cmessage, filename, lut)

USE um_parcore, ONLY: mype
USE ukca_radaer_tlut_mod
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE filenamelength_mod, ONLY: filenamelength

IMPLICIT NONE

!
! Arguments
!

! error indicator (0 is OK, <O is KO)
INTEGER :: icode

! error message
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=errormessagelength) :: iomessage

! look-up table filename
CHARACTER (LEN=filenamelength) :: filename

! Structure representing a look-up table
TYPE (ukca_radaer_tlut) :: lut


INTEGER :: ios
INTEGER :: i
INTEGER :: j
INTEGER :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_LUT_IN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!
! Get the file name, open the file and perform a namelist read.
!
IF (mype==0) THEN
  IF (LEN_TRIM(filename) == 0) THEN
    icode = -1
    cmessage = 'UKCA look-up table filename is blank'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, &
                            zhook_handle)
    RETURN
  END IF
  CALL assign_file_unit(filename, funit, handler="fortran")
  OPEN(UNIT=funit, FILE=filename, ACTION='READ', IOSTAT=ios, IOMSG=iomessage)
  IF (ios /= 0) THEN
    icode = -1
    cmessage =                                                     newline//&
      'Error opening UKCA look-up table file:'//                   newline//&
      TRIM(filename) //                                            newline//&
      'IoMsg: '//TRIM(iomessage)
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, &
                            zhook_handle)
    RETURN
  END IF

END IF  ! mype==0

CALL read_nml_ukcanml(funit)
IF (mype==0) THEN
  CLOSE(funit)
  CALL release_file_unit(funit, handler="fortran")
END IF

!
! Check that actual array dimensions do not exceed
! those that were fixed at compile time.
!
IF (n_x > npd_x) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (X).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

IF (n_nr > npd_nr) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (NR).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

IF (n_ni > npd_ni) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (NI).'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! Copy the look-up table into the structure passed as an argument
lut%stdev = stdev

!
! Namelist arrays indices start at 0. For the sake of ease of use,
! online arrays start at 1 -> simply shift the content by 1 and add
! 1 to each dimension.
! All calculations within UKCA_RADAER expect array indexing starting
! at 1.
!
n_x = n_x + 1
n_nr = n_nr + 1
n_ni = n_ni + 1

lut%n_x  = n_x
lut%n_nr = n_nr
lut%n_ni = n_ni

lut%x_min  = x_min
lut%x_max  = x_max
lut%nr_min = nr_min
lut%nr_max = nr_max
lut%incr_nr = (nr_max - nr_min) / REAL(n_nr-1)
lut%ni_min = ni_min
lut%ni_max = ni_max
lut%ni_c   = ni_c

! Allocate the dynamic arrays
ALLOCATE(lut%ukca_absorption(n_x, n_ni, n_nr))
ALLOCATE(lut%ukca_scattering(n_x, n_ni, n_nr))
ALLOCATE(lut%ukca_asymmetry(n_x, n_ni, n_nr))
ALLOCATE(lut%volume_fraction(n_x))

DO k = 1, n_nr
  DO j = 1, n_ni
    DO i = 1, n_x
      lut%ukca_absorption(i, j, k) = ukca_absorption(i-1, j-1, k-1)
      lut%ukca_scattering(i, j, k) = ukca_scattering(i-1, j-1, k-1)
      lut%ukca_asymmetry(i, j, k)  = ukca_asymmetry(i-1, j-1, k-1)
    END DO
  END DO
END DO

DO i = 1, n_x
  lut%volume_fraction(i) = volume_fraction(i-1)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_lut_in

SUBROUTINE read_nml_ukcanml(unit_in)
USE um_parcore, ONLY: mype
USE setup_namelist, ONLY: setup_nml_type
IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_UKCANML'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 3
INTEGER, PARAMETER :: n_real = 8 + (npd_x+1) +                      &
                    3 * ( (npd_x+1) * (npd_ni+1) * (npd_nr+1) )

TYPE my_namelist
  SEQUENCE
  INTEGER ::  n_x
  INTEGER ::  n_nr
  INTEGER ::  n_ni
  REAL :: stdev
  REAL :: x_min
  REAL :: x_max
  REAL :: nr_min
  REAL :: nr_max
  REAL :: ni_min
  REAL :: ni_max
  REAL :: ni_c
  REAL :: ukca_absorption(0:npd_x, 0:npd_ni, 0:npd_nr)
  REAL :: ukca_scattering(0:npd_x, 0:npd_ni, 0:npd_nr)
  REAL :: ukca_asymmetry(0:npd_x, 0:npd_ni, 0:npd_nr)
  REAL :: volume_fraction(0:npd_x)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
                    n_real_in=n_real) 

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=ukcanml, IOSTAT=ErrorStatus)

  my_nml % n_x = n_x
  my_nml % n_nr = n_nr
  my_nml % n_ni = n_ni
  my_nml % stdev = stdev
  my_nml % x_min = x_min
  my_nml % x_max = x_max
  my_nml % nr_min = nr_min
  my_nml % nr_max = nr_max
  my_nml % ni_min = ni_min
  my_nml % ni_max = ni_max
  my_nml % ni_c   = ni_c
  my_nml % ukca_absorption = ukca_absorption
  my_nml % ukca_scattering = ukca_scattering
  my_nml % ukca_asymmetry = ukca_asymmetry
  my_nml % volume_fraction =  volume_fraction

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  n_x  = my_nml % n_x
  n_nr  = my_nml % n_nr
  n_ni  = my_nml % n_ni
  stdev  = my_nml % stdev
  x_min  = my_nml % x_min
  x_max  = my_nml % x_max
  nr_min  = my_nml % nr_min
  nr_max  = my_nml % nr_max
  ni_min  = my_nml % ni_min
  ni_max  = my_nml % ni_max
  ni_c    = my_nml % ni_c
  ukca_absorption  = my_nml % ukca_absorption
  ukca_scattering  = my_nml % ukca_scattering
  ukca_asymmetry  = my_nml % ukca_asymmetry
  volume_fraction  = my_nml % volume_fraction

END IF
   
CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_nml_ukcanml

SUBROUTINE ukca_radaer_get_lut_index_scalar(Nni, im_m, ni_min, ni_max, ni_c,  &
                                            ni_ind, ni_c_power)
IMPLICIT NONE

 INTEGER, INTENT(IN)  :: Nni       ! number of ni values in LUT
 REAL,    INTENT(IN)  :: im_m      ! calculated value of ni
 REAL,    INTENT(IN)  :: ni_min    ! minimum value of ni in the LUT
 REAL,    INTENT(IN)  :: ni_max    ! maximum value of ni in the LUT
 REAL,    INTENT(IN)  :: ni_c   ! parameter required to calculate ni increments
 INTEGER, INTENT(OUT) :: ni_ind    ! nearest-neighbour index for ni
 REAL,    INTENT(IN), OPTIONAL :: ni_c_power

 REAL    :: im_m_array(1)
 INTEGER :: ni_ind_array(1)

 im_m_array(1) = im_m
 CALL ukca_radaer_get_lut_index_array(Nni, im_m_array, ni_min, ni_max, ni_c, &
                                      ni_ind_array, 1, ni_c_power=ni_c_power)
 ni_ind = ni_ind_array(1)

RETURN
END SUBROUTINE ukca_radaer_get_lut_index_scalar

SUBROUTINE ukca_radaer_get_lut_index_array(Nni, im_m, ni_min, ni_max, ni_c,  &
                                           ni_ind, length, ni_c_power)
USE vectlib_mod, only: log_v
IMPLICIT NONE

 INTEGER, INTENT(IN)  :: length           ! Array length
 INTEGER, INTENT(IN)  :: Nni              ! number of ni values in LUT
 REAL,    INTENT(IN)  :: im_m(1:length)   ! calculated value of ni
 REAL,    INTENT(IN)  :: ni_min           ! minimum value of ni in the LUT
 REAL,    INTENT(IN)  :: ni_max           ! maximum value of ni in the LUT
 REAL,    INTENT(IN)  :: ni_c    ! parameter required to calculate ni increments
 INTEGER, INTENT(OUT) :: ni_ind(1:length) ! nearest-neighbour index for ni
 REAL,    INTENT(IN), OPTIONAL :: ni_c_power

! local variables
REAL :: a, b, incr_ni
REAL :: ni_c_power_local
REAL, PARAMETER :: min_ni_c = 0.001 ! Lowest value of ni_c to accept
REAL, PARAMETER :: max_ni_c = 5.0   ! Highest value of ni_c to accept

REAL, PARAMETER :: inv_ln_10 = 1.0 / LOG(10.0)

REAL :: logs_array_in(1:length)
REAL :: logs_array_out(1:length)

INTEGER :: icode  
CHARACTER(LEN=errormessagelength) :: cmessage

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_GET_LUT_INDEX'

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

icode = 0

! If ni_c > max_ni_c then produce a fatal error
! 
IF (ni_c > max_ni_c) THEN 
  icode = 1
  cmessage='UKCA RADAER Look-up table'//                   newline//&
           'NI_C exceeds upper limit'

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END IF

! As ni_c tends towards zero the function below tends towards a linear scale
! With ni_c = 0, a and b become infinity. So we will protect this
! from happening by using a linear scale if ni_c < min_ni_c
!
IF (ni_c > min_ni_c) THEN

   IF(PRESENT(ni_c_power)) THEN
     ni_c_power_local = ni_c_power
   ELSE
     ni_c_power_local = 10.0**ni_c
   END IF

   a = ni_max / ((ni_c_power_local) - 1.)
   b = REAL(Nni) / ni_c

   logs_array_in(:) = (iM_m(:)/a) + 1.0

   CALL log_v(length, logs_array_in, logs_array_out)

   ni_ind(:) = NINT(b * logs_array_out(:) * inv_ln_10 ) + 1

ELSE

   incr_ni   = (ni_max - ni_min) / REAL(Nni-1)
   ni_ind(:) = NINT( (im_m(:) - ni_min) / incr_ni ) + 1

END IF

ni_ind(:) = MIN(Nni, MAX(1, ni_ind(:)))

IF (icode /= 0) THEN
   CALL Ereport(RoutineName,icode,cmessage)
END IF

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_get_lut_index_array

END MODULE ukca_radaer_lut_read_in
