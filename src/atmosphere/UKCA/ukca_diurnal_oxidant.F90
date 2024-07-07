! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description:
!  The routine UKCA_SET_DIURNAL_OX distributes the mean oxidant concentrations
!  read in from external files to give a diurnally varying concentration.
!  Part of the offline oxidant chemistry scheme. A routine to calculate the time
!  integral of cos(zenith angle) is also included (UKCA_INT_COSZ).
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds, University of Oxford,
!  and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines:
!   ukca_set_diurnal_ox
!   ukca_int_cosz
!
!-----------------------------------------------------------------------
!
MODULE ukca_diurnal_oxidant

USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

IMPLICIT NONE
PRIVATE

! cos (zenith angle) and integral of cos(zenith angle) for diurnal variation
REAL, PUBLIC, ALLOCATABLE, SAVE :: cosza(:)          ! cos(za)
REAL, PUBLIC, ALLOCATABLE, SAVE :: intcosza(:)       ! integral cos(za)
REAL, PUBLIC, ALLOCATABLE, SAVE :: daylength1d(:)    ! radians

PUBLIC ukca_set_diurnal_ox                ! called from asad_fyinit
PUBLIC ukca_set_diurnal_ox_col            ! called from asad_fyinit
PUBLIC ukca_int_cosz                      ! called from ukca_main1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_DIURNAL_OXIDANT'

CONTAINS

! #############################################################################

SUBROUTINE ukca_set_diurnal_ox(species, field, n_points, nlev)

! This routine provides code to add a diurnal profile to oxidants provided
! as mean fields. This varies with each oxidant. In the case of OH and HO2
! the diurnal profile is proportional to cos(zenith angle).

USE ukca_chem_offline,        ONLY: o3_offline_diag, oh_offline_diag,          &
                                    no3_offline_diag, ho2_offline_diag
USE conversions_mod,          ONLY: rsec_per_day, pi
USE ereport_mod,              ONLY: ereport

USE errormessagelength_mod,   ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_points                ! field dimension
INTEGER, INTENT(IN)  :: nlev                    ! model level
CHARACTER (LEN=10), INTENT(IN)  :: species      ! Species char strng

REAL,    INTENT(INOUT) :: field(n_points)       ! field to be scaled

! Local variables

REAL, PARAMETER :: daylen_cutoff = 3600.0
REAL            :: daylen_secs(n_points)        ! daylength (s)

CHARACTER (LEN=errormessagelength) :: cmessage                  ! Error message
INTEGER            :: errcode                   ! Variable passed to ereport

REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_DIURNAL_OX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(species(1:4))
CASE ('OH  ','HO2 ')

  ! Scale the field according to cos(zenith angle)/integral(cos(zenith_angle))
  ! Night time values are zero
  WHERE ((cosza > 0.0) .AND. (intcosza > 1.0e-1))
    field(:) = field(:)*rsec_per_day*(cosza(:)/intcosza(:))
  ELSE WHERE
    field(:) = 0.0
  END WHERE

CASE ('NO3 ')

  ! convert from radians to seconds
  daylen_secs = daylength1d*rsec_per_day/(2.0*pi)

  ! With short day length ensures no NO3 diurnal variation done
  WHERE (daylen_secs < daylen_cutoff)
    daylen_secs = rsec_per_day
  END WHERE

  ! Restrict NO3 to nightime
  WHERE ((cosza > 0.0) .AND. (intcosza > 1.0e-1) .AND.                       &
    (daylen_secs < rsec_per_day))
    field(:) = 0.0
  END WHERE

CASE ('O3  ')
  ! No diurnal profile for O3

CASE DEFAULT

  cmessage = 'ukca_set_diurnal_ox: Species: '//species//' not found'
  errcode = 1
  CALL ereport('ukca_set_diurnal_ox',errcode,cmessage)

END SELECT

! Get rid of any negative values
WHERE (field < 0.0)
  field = 0.0
END WHERE

! Set diagnostic fields
SELECT CASE(species(1:4))
CASE ('O3  ')
  o3_offline_diag(:,nlev) = field(:)
CASE ('OH  ')
  oh_offline_diag(:,nlev) = field(:)
CASE ('NO3 ')
  no3_offline_diag(:,nlev) = field(:)
CASE ('HO2 ')
  ho2_offline_diag(:,nlev) = field(:)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_diurnal_ox

! ######################################################################

SUBROUTINE ukca_set_diurnal_ox_col(species, field, n_points, ix, jy)

! This routine provides code to add a diurnal profile to oxidants provided
! as mean fields. This varies with each oxidant. In the case of OH and HO2
! the diurnal profile is proportional to cos(zenith angle).

USE ukca_chem_offline,      ONLY: o3_offline_diag, oh_offline_diag,           &
                                  no3_offline_diag, ho2_offline_diag
USE conversions_mod,        ONLY: pi, rsec_per_day
USE ereport_mod,            ONLY: ereport

USE nlsizes_namelist_mod,   ONLY: row_length
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n_points                ! field dimension
INTEGER, INTENT(IN) :: ix                       ! i counter
INTEGER, INTENT(IN) :: jy                       ! j counter
CHARACTER (LEN=10), INTENT(IN)  :: species      ! Species char strng

REAL,    INTENT(INOUT) :: field(n_points)       ! field to be scaled

! Local variables

REAL, PARAMETER :: daylen_cutoff = 3600.0
REAL            :: daylen_secs        ! daylength (s)

CHARACTER (LEN=errormessagelength) :: cmessage                  ! Error message
INTEGER            :: errcode                   ! Variable passed to ereport

INTEGER :: location ! mapping to theta_field

REAL(KIND=jprb)    :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_DIURNAL_OX_COL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! get mapping from theta_field to proper x*y. Field is 1 column here!
location = ix + ((jy - 1)*row_length)

SELECT CASE(species(1:4))
CASE ('OH  ','HO2 ')

  ! Scale the field according to cos(zenith angle)/integral(cos(zenith_angle))
  ! Night time values are zero
  IF ((cosza(location) > 0.0) .AND. (intcosza(location) > 1.0e-1)) THEN
    field(:) = field(:)*rsec_per_day*(cosza(location) / intcosza(location))
  ELSE
    field(:) = 0.0
  END IF

CASE ('NO3 ')

  ! convert from radians to seconds
  daylen_secs = daylength1d(location)*rsec_per_day / (2.0*pi)

  ! With short day length ensures no NO3 diurnal variation done
  IF (daylen_secs < daylen_cutoff) daylen_secs = rsec_per_day

  ! Restrict NO3 to nightime
  IF (((cosza(location) > 0.0) .AND. (intcosza(location) > 1.0e-1)) .AND.     &
       (daylen_secs < rsec_per_day)) THEN
    field(:) = 0.0
  END IF

CASE ('O3  ')
  ! No diurnal profile for O3

CASE DEFAULT

  cmessage = 'ukca_set_diurnal_ox: Species: '//species//' not found'
  errcode = 1
  CALL ereport('ukca_set_diurnal_ox',errcode,cmessage)

END SELECT

! Get rid of any negative values
WHERE (field < 0.0)
  field = 0.0
END WHERE

! Set diagnostic fields
SELECT CASE(species(1:4))
CASE ('O3  ')
  o3_offline_diag(location,:) = field(:)
CASE ('OH  ')
  oh_offline_diag(location,:) = field(:)
CASE ('NO3 ')
  no3_offline_diag(location,:) = field(:)
CASE ('HO2 ')
  ho2_offline_diag(location,:) = field(:)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_diurnal_ox_col

! ######################################################################

SUBROUTINE ukca_int_cosz(row_length,rows,true_latitude,sindec,                 &
  cos_zenith_angle, int_zenith_angle, daylength)

!  Description:
!   Calculate integral of cos(solar zenith angle) dt over a period of
!   one day multiplied by seconds per day.
!   See Gates, D.M. (1980) Biophysical Ecology, Springer-Verlag, p103-4
!   Only depends on latitude and solar declination, original expression:
!   int_za = (daysecs/pi)*(hs*sin(phi)*sin(decls)+cos(phi)*cos(decls)*sin(hs)),
!   was simplified to:
!   int_za = (daysecs/pi)*0.5*(cos(phi-decls)*(hs + sin(hs)) +
!                              cos(phi+decls)*(sin(hs) - hs)).
!   where hs is day length / 2.0.
!
! ######################################################################
!

USE ukca_option_mod,  ONLY: l_ukca_offline, l_ukca_offline_be
USE parcons_mod,      ONLY: pi
IMPLICIT NONE

INTEGER, INTENT(IN) ::  row_length       ! field dimension
INTEGER, INTENT(IN) ::  rows             ! field dimension

REAL, INTENT(IN)    :: true_latitude(row_length, rows)    ! latitude
REAL, INTENT(IN)    :: sindec            ! Sin(solar declination)
REAL, INTENT(IN)    :: cos_zenith_angle(row_length, rows)
! Cosine of zenith angle

REAL, INTENT(OUT)   :: int_zenith_angle(row_length, rows)
! Integral of zenith angle
REAL, INTENT(OUT)   :: daylength(row_length, rows)
! Daylength

! Local variables

REAL, PARAMETER :: daysecs=24.0*3600.0         ! seconds in day
REAL            :: decls                       ! solar declination
REAL            :: hs(row_length,rows)         ! half day length (radians)
REAL            :: cos_hs(row_length,rows)     ! cos of hs
REAL            :: intcosza2d(row_length,rows) ! integral of cos(za)

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_INT_COSZ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

decls = ASIN(sindec)
cos_hs(:,:) = -TAN(true_latitude(:,:))*TAN(decls)

! Polar summer, set cos_hs ~ -1
WHERE (ABS(cos_hs) > 1.0 .AND. true_latitude*decls > 0.0)
  cos_hs = -0.99999999
END WHERE

WHERE (ABS(cos_hs) < 1.0)
  hs = ACOS(cos_hs)
  int_zenith_angle = (daysecs/pi)*0.5*(COS(true_latitude-decls)*               &
    (hs+SIN(hs)) + COS(true_latitude+decls)*(SIN(hs)-hs))
  daylength = 2.0*hs
ELSE WHERE
  ! Polar night
  int_zenith_angle = 0.0
  daylength = 0.0
END WHERE

IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
  ! Make 1-D arays for set_ukca_diurnal_ox routine
  IF (.NOT. ALLOCATED(intcosza)) ALLOCATE(intcosza(row_length*rows))
  IF (.NOT. ALLOCATED(cosza)) ALLOCATE(cosza(row_length*rows))
  IF (.NOT. ALLOCATED(daylength1d)) ALLOCATE(daylength1d(row_length*rows))

  intcosza = RESHAPE(int_zenith_angle,(/row_length*rows/))
  cosza = RESHAPE(cos_zenith_angle,(/row_length*rows/))
  daylength1d = RESHAPE(daylength,(/row_length*rows/))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_int_cosz
END MODULE ukca_diurnal_oxidant
