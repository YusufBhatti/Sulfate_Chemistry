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
!   Called from UKCA_MAIN
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
!      n_lat = 36          number of latitudes in climatology
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
MODULE ukca_read_aerosol_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_READ_AEROSOL_MOD'

CONTAINS

SUBROUTINE ukca_read_aerosol(year,month,row_length,rows,           &
                        model_levels,sin_theta_latitude,           &
                        height, use_background, aerosol_sd)

USE ukca_option_mod,        ONLY: dir_strat_aer, file_strat_aer,     &
                            i_ukca_sad_months, i_ukca_sad_start_year
USE UM_ParVars
USE UM_ParCore,             ONLY: mype, nproc
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE ereport_mod,            ONLY: ereport
USE conversions_mod,        ONLY: recip_pi_over_180
USE file_manager,           ONLY: assign_file_unit, release_file_unit
USE filenamelength_mod,     ONLY: filenamelength
USE vectlib_mod,            ONLY: asin_v
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



INTEGER, INTENT(IN) :: month   ! date required
INTEGER, INTENT(IN) :: year    !
! Model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

REAL, INTENT(IN) :: sin_theta_latitude(row_length, rows)
REAL, INTENT(IN) :: height(model_levels)
LOGICAL, INTENT(IN) :: use_background ! T: use background file (equiv.
                                      ! to 2100 value from full file)

REAL, INTENT(OUT) :: aerosol_sd(row_length, rows, model_levels)

! local variables

INTEGER, PARAMETER :: n_lat = 36 ! number of latitudes in clim.
REAL, PARAMETER :: dlat = 5.0    ! latitude spacing in deg.
REAL, PARAMETER :: first_lat = -87.5 ! first latitude
REAL, SAVE :: lat(n_lat)

INTEGER, PARAMETER :: n_alt = 56 ! number of altitudes
REAL, PARAMETER :: dalt = 500.0  ! m; vertical spacing
REAL, PARAMETER :: first_alt = 12000.0 ! first altitude
REAL, SAVE :: alt(n_alt)

INTEGER, PARAMETER :: n_months_bg   = 12   ! number of months for background
INTEGER, SAVE :: n_months

REAL, ALLOCATABLE, SAVE :: sad_clim_3D(:,:,:)
INTEGER :: month_index
INTEGER :: istat
INTEGER :: ukcastar_unit
LOGICAL, SAVE :: first = .TRUE.

INTEGER :: i,j

CHARACTER(LEN=80) :: one_line
CHARACTER(LEN=filenamelength) :: filename

REAL :: sad_clim(n_lat, n_alt)

REAL :: sad_clim1(n_lat, model_levels)

INTEGER :: INDEX
REAL :: frac

REAL, ALLOCATABLE, SAVE :: latitude(:,:)

! Unit conversion: Revert from mum^2/cm^3 to cm^2/cm^3
REAL, PARAMETER :: scalefac = 1.0e-8

! For errors
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: ierr
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)         :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_READ_AEROSOL'

! main program body
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ((.NOT. use_background) .AND. (year < i_ukca_sad_start_year)) THEN
  cmessage='MODEL YEAR IS INCOMPATABLE WITH AEROSOL'
  ierr    = 34999
  CALL ereport('UKCA_READ_AEROSOL',ierr,cmessage)
END IF


IF (first) THEN
  ! get correct number of months
  IF (use_background) THEN
    n_months = n_months_bg
  ELSE
    n_months = i_ukca_sad_months
  END IF

  ALLOCATE (sad_clim_3d(n_lat, n_alt, n_months))

  DO i=1,n_lat
    lat(i) = first_lat + dlat * REAL(i-1)
  END DO
  DO i=1,n_alt
    alt(i) = first_alt + dalt * REAL(i-1)
  END DO
  first = .FALSE.

  ALLOCATE(latitude(row_length, rows))

  CALL asin_v(row_length*rows, sin_theta_latitude, latitude)

  latitude = latitude * recip_pi_over_180

  IF (mype == 0) THEN
    ! read and keep climatology
    filename = TRIM(ADJUSTL(dir_strat_aer))//      &
                    '/'//TRIM(ADJUSTL(file_strat_aer))
    CALL assign_file_unit(filename, ukcastar_unit, handler="fortran")
    OPEN(ukcastar_unit,FILE=filename, ACTION='READ')
    ! read data
    READ(ukcastar_unit,*) sad_clim_3d
    CLOSE(ukcastar_unit)
    CALL release_file_unit(ukcastar_unit, handler="fortran")
    WHERE (sad_clim_3d < 0) sad_clim_3d = 0.0

    ! revert to cm^2/cm^3
    sad_clim_3d = sad_clim_3d * scalefac
  END IF
  CALL gc_rbcast(1,n_lat*n_alt*n_months,0,nproc,istat,            &
                 sad_clim_3d)
END IF  ! first

! main body, executed every timestep

! calculate monthly index
IF (use_background) THEN
  month_index = month
ELSE
  month_index = (year - i_ukca_sad_start_year)*12 + month
  DO WHILE (month_index < 1)
    month_index = month_index + 12
  END DO
  DO WHILE (month_index > n_months)
    month_index = month_index - 12
  END DO
END IF

! recover climatology for chosen month
sad_clim = sad_clim_3d(:,:,month_index)

aerosol_sd(:,:,:) = 0.0

! interpolate to UM grid
! interpolate in height

! Generate climatology on UM height levels
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
      aerosol_sd(j,i,:) = sad_clim1(1,:)
    ELSE IF (latitude(j,i) >= lat(n_lat)) THEN
      aerosol_sd(j,i,:) = sad_clim1(n_lat,:)
    ELSE
      INDEX = MAX(MIN(FLOOR((latitude(j,i) -                     &
           first_lat)/dlat)+1,n_lat-1),1)
      frac = MAX(MIN((latitude(j,i) - lat(INDEX))/dlat,1.0),0.0)
      aerosol_sd(j,i,:) = (1.0 - frac) * sad_clim1(INDEX    ,:)+  &
           frac  * sad_clim1(INDEX + 1,:)
    END IF
  END DO
END DO

WHERE (aerosol_sd < 0.0) aerosol_sd = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_read_aerosol
END MODULE ukca_read_aerosol_mod
