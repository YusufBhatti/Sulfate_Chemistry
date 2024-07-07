! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Read routine for the CCMVal-2 SPARC aerosol climatology.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!   Called from UKCA_FASTJX
!
!
!     Purpose.
!     --------
!     To read the SPARC aerosol climatology at the start of a run and to
!     provide spatial and temporal interpolates during runtime.
!
!     Interface
!     ---------
!     Called from *ukca_main*.
!
!
!     Local variables
!     ---------------
!
!      n_lat = 36          number of latitudes in clim.
!      dlat = 5.          latitude spacing in deg.
!      first_lat = -87.5  first latitude
!      lat(n_lat)         latitudes in climatology

!      n_alt = 56         number of altitudes
!      dalt = 500.        vertical spacing
!      first_alt = 12000. first altitude
!      alt(n_alt)         altitudes in climatology

!      n_months = 1812    number of months
!      sad_clim_3D(:,:,:) allocated and saved array, containing the climatology
!      month_index        time index
!      istat              error indicator
!
!      first              set to FALSE after first call

!      filename = 'Sulfate_SAD_SPARC_1950-2100.asc' Name of clim. file
!      i,j                integer counters
!      one_line           string variable for reading
!      sad_clim(n_lat, n_alt) temporary interpolated fields
!      sad_clim1(n_lat, model_levels)
!      index, frac        temporary variables

!      latitude           allocated model latitude

! Unit conversion: Revert from mum^2/cm^3 to cm^2/cm^3
!      scalefac = 1.0e-8
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_read_reff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_READ_REFF_MOD'

CONTAINS

SUBROUTINE ukca_read_reff(year,month,row_length,rows,             &
                        model_levels,sin_theta_latitude,          &
                        height, aerosol_reff)

USE ukca_option_mod,      ONLY: dir_reff_sulp, file_reff_sulp
USE UM_ParVars
USE UM_ParCore,           ONLY: mype, nproc
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE conversions_mod,      ONLY: recip_pi_over_180
USE file_manager,         ONLY: assign_file_unit, release_file_unit
USE filenamelength_mod,   ONLY: filenamelength
USE vectlib_mod, ONLY: asin_v
IMPLICIT NONE



INTEGER, INTENT(IN) :: month                     ! date of aerosol density
INTEGER, INTENT(IN) :: year                      ! date of aerosol density

INTEGER, INTENT(IN) :: row_length                ! model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

REAL, INTENT(IN)    :: sin_theta_latitude(row_length, rows)
REAL, INTENT(IN)    :: height(model_levels)

REAL, INTENT(OUT)   :: aerosol_reff(row_length, rows, model_levels)

! local variables

INTEGER, PARAMETER :: n_lat = 36        ! number of latitudes in clim.
REAL, PARAMETER :: dlat = 5.0            ! latitude spacing in deg.
REAL, PARAMETER :: first_lat = -87.5    ! first latitude
REAL, SAVE :: lat(n_lat)

INTEGER, PARAMETER :: n_alt = 56        ! number of altitudes
REAL, PARAMETER :: dalt = 500.0          ! m; vertical spacing
REAL, PARAMETER :: first_alt = 12000.0   ! first altitude
REAL, SAVE :: alt(n_alt)

INTEGER, PARAMETER :: n_months = 1812   ! number of months
REAL, ALLOCATABLE, SAVE :: sad_clim_3D(:,:,:)
INTEGER :: month_index
INTEGER :: istat

LOGICAL, SAVE :: first = .TRUE.


INTEGER :: i,j
INTEGER :: ukcafjar_unit

CHARACTER(LEN=80) :: one_line
CHARACTER(LEN=filenamelength) :: filename

REAL :: sad_clim(n_lat, n_alt)

REAL :: sad_clim1(n_lat, model_levels)

INTEGER :: INDEX
REAL :: frac

REAL, ALLOCATABLE, SAVE :: latitude(:,:)

! Unit conversion: Revert from mum^2/cm^3 to cm^2/cm^3
REAL, PARAMETER :: scalefac = 1.0e-8

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)         :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_READ_REFF'

! main program body
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first) THEN
  ALLOCATE (sad_clim_3d(n_lat, n_alt, n_months))

  DO i=1,n_lat
    lat(i) = first_lat + dlat * REAL(i-1)
  END DO
  DO i=1,n_alt
    alt(i) = first_alt + dalt * REAL(i-1)
  END DO
  first = .FALSE.

  ALLOCATE(latitude(row_length, rows))
  CALL asin_v(row_length*rows, sin_theta_latitude,latitude)

  latitude = latitude * recip_pi_over_180

  IF (mype == 0) THEN
    ! read and keep climatology
    filename = TRIM(ADJUSTL(dir_reff_sulp))//         &
               TRIM(ADJUSTL(file_reff_sulp))
    CALL assign_file_unit(filename, ukcafjar_unit, handler="fortran")
    OPEN(ukcafjar_unit,FILE=filename,ACTION='READ')

    ! read data
    READ(ukcafjar_unit,*) sad_clim_3d
    CLOSE(ukcafjar_unit)
    CALL release_file_unit(ukcafjar_unit, handler="fortran")

    WHERE (sad_clim_3d < 0) sad_clim_3d = 0.0

    ! revert to cm^2/cm^3
    sad_clim_3d = sad_clim_3d * scalefac
  END IF
  CALL gc_rbcast(1,n_lat*n_alt*n_months,0,nproc,istat,            &
                 sad_clim_3d)
END IF  ! first

! main body, executed every timestep

! calculate monthly index
month_index = (year - 1950)*12 + month
DO WHILE (month_index < 1)
  month_index = month_index + 12
END DO
DO WHILE (month_index > n_months)
  month_index = month_index - 12
END DO

! recover climatlogy for chosen month
sad_clim = sad_clim_3d(:,:,month_index)

! interpolate to UM grid
! interpolate in height

! Generate climatlogy on UM height levels
DO i=1,model_levels
  IF ((height(i) < first_alt) .OR.                                 &
      (height(i) > alt(n_alt))) THEN
    sad_clim1(:,i) = 0.0
  ELSE IF (height(i) == alt(n_alt)) THEN
    sad_clim1(:,i) = sad_clim(:,n_alt)
  ELSE
    INDEX = INT ((height(i) - first_alt)/dalt) + 1
    frac = (height(i) - alt(INDEX))/dalt
    sad_clim1(:,i) = (1.0 - frac) * sad_clim(:,INDEX) +             &
                           frac  * sad_clim(:,INDEX + 1)
  END IF
END DO


! interpolate in latitude

DO i=1,rows
  DO j=1,row_length
    IF (latitude(j,i) <= first_lat) THEN
      aerosol_reff(j,i,:) = sad_clim1(1,:)
    ELSE IF (latitude(j,i) >= lat(n_lat)) THEN
      aerosol_reff(j,i,:) = sad_clim1(n_lat,:)
    ELSE
      INDEX = INT ((latitude(j,i) - first_lat)/dlat + 0.001)+1
      frac = (latitude(j,i) - lat(INDEX))/dlat
      aerosol_reff(j,i,:) = (1.0 - frac) * sad_clim1(INDEX    ,:)+    &
                                frac  * sad_clim1(INDEX + 1,:)
    END IF
  END DO
END DO

WHERE (aerosol_reff < 0.0) aerosol_reff = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_read_reff

END MODULE ukca_read_reff_mod
