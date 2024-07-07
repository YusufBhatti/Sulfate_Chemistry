! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Reflectivity calculation

MODULE max_hail_size_mod
! Description:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

INTEGER, PARAMETER :: nbins = 50  ! Number of bins

REAL     :: hail_sizes(nbins) ! Selection of hail size bins [m] - largest first
REAL     :: gwc_min(nbins)    ! Minimum graupel water content value
                              ! for each bin [kg m-3] - largest first

REAL, PARAMETER :: m2mm = 1.0E3 ! Conversion factor from [m] to [mm]

REAL, PARAMETER :: size_largest  = 0.075    ! Size of the largest hail bin [m]
REAL, PARAMETER :: size_smallest = 500.0E-6 ! Size of the smallest hail bin [m]

REAL, PARAMETER :: gwc_abs = 1.0E-10
! An absolute, minimum value of graupel mass [kg m-3] that is considered
! as having some physical size - below this threshold maximum diameter
! will be zero.

REAL, PARAMETER :: threshold_concentration = 0.005
! A threshold concentration [m-3]. The maximum hail size is deemed to be the
! diameter at which this threshold is exceeded, starting with the largest
! particles first.

LOGICAL  :: l_hail_bins_set = .FALSE. ! Flag for hail bins set
                                      ! Initially false

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MAX_HAIL_SIZE_MOD'

CONTAINS

!==============================================================================
SUBROUTINE generate_hail_bins()

! Generates a lookup table of minimum hail size for each bin

USE mphys_psd_mod,       ONLY: x1g, x2g, x4g
USE mphys_constants_mod, ONLY: cx, constp
USE missing_data_mod,    ONLY: rmdi
USE umprintmgr,          ONLY: umprint, newline
USE ereport_mod,         ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL :: bin_edges(nbins + 1) ! Temporary working edge of each bin [m]
REAL :: bin_widths(nbins)    ! Width of each size bin [m]
REAL :: ls_frac              ! Fraction: (size_largest / size_smallest)
REAL :: log_smallest         ! Natural Log (ln) of size_smallest
REAL :: binfrac              ! Fraction: bin number / number of bins
REAL :: working_gwc          ! Local working graupel water content [kg m-3]
REAL :: working_gwc_inc      ! Local working increment of graupel water content
                             ! [kg m-3]
REAL :: lamg                 ! Lambda (Slope parameter) for graupel particle
                             ! size distribution [m-1]
REAL :: n0g                  ! Intercept parameter for the graupel particle
                             ! size distribution [m-3]
REAL :: diam_g               ! Local graupel diameter [m]
REAL :: sum_in_bin           ! Local number of particles in bin [m-3]
REAL :: sum_n0g              ! Sum of the number of particles when integrating
                             ! from the tail of the graupel particle size
                             ! distribution [m-3]

! Factor for use in increment calculations. The smaller this value is, the
! more accurate the solution will be, but the longer and more expensive the
! computational calculations will be.
REAL, PARAMETER :: inc_fact = 1.0E-2

INTEGER :: l                    ! hail bin loop counter
INTEGER :: ErrorStatus          ! Error Status, set initially to normal
INTEGER :: previous_bin         ! Previous bin number used in loop calculations
INTEGER :: working_bin          ! Current bin number used within loop
                                ! calculations
INTEGER :: iter_count           ! Number of interations counted in DO WHILE
                                ! loop

CHARACTER(LEN=errormessagelength) :: emessage    ! Error message text
CHARACTER(LEN=1000)               :: line_buffer ! Text for UM print

LOGICAL :: l_bins_incorrect     ! Flag for incorrectly set hail bins

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='GENERATE_HAIL_BINS'

! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initially set the flag for incorrectly set bins to False and error status
! to normal:
l_bins_incorrect = .FALSE.
ErrorStatus      = 0

! Set hail sizes initially to missing data and check for negative hail sizes
! later as this can be used as a flag in case hail sizes are incorrectly set.
hail_sizes(:) = rmdi

! Set gwc_min for each bin this to zero - this means that any unset values
! of gwc_min should just produce a hail size of zero.
gwc_min(:) = 0.0

! Initialise the smallest hail and largest bin edges - biggest first
bin_edges(1)       = size_smallest
bin_edges(nbins+1) = size_largest

ls_frac      = LOG(size_largest / size_smallest)
log_smallest = LOG(size_smallest)

! Fill in the missing bin edges
DO l = 2, nbins
  binfrac      = REAL(l-1) / REAL(nbins)
  bin_edges(l) = EXP( binfrac * ls_frac + log_smallest)
END DO

DO l = 1, nbins
  hail_sizes(l) = SQRT( bin_edges(l)*bin_edges(l+1) )
  bin_widths(l) = bin_edges(l+1) - bin_edges(l)

  ! If this bin has been set incorrectly, it will still be zero.
  ! Raise the flag ready to call ereport as this can only be due
  ! to someone doing something weird by editing the hail size
  ! code.
  IF ( hail_sizes(l) < 0.0 ) l_bins_incorrect = .TRUE.

END DO

! Check the maximum hail size is not negative (this will indicate
! missing data) and so need to call ereport

IF (l_bins_incorrect ) THEN
  ErrorStatus = 100 ! Flag to cause an issue

  emessage = 'Negative values found in hail size bins. '// newline//          &
     'This suggests that they have been incorrectly set'// newline//          &
     'or remain as missing data. Please check setup of '// newline//          &
     'hail bins in generate_hail_bins subroutine'

  CALL ereport(RoutineName, ErrorStatus, emessage)

END IF

WRITE(line_buffer, '(A)')                                                     &
' About to perform look-up calculation for maximum hail size'
CALL umprint(line_buffer, src=RoutineName)

! Starting with the smallest possible value (gwc_abs) work out the
! bin this will lie in and therefore a minimum value of

working_gwc  = gwc_abs
previous_bin = 0
iter_count   = 0

DO WHILE ( previous_bin < nbins )

  ! First increment the iteration count and set the graupel increment
  iter_count      = iter_count + 1
  working_gwc_inc = working_gwc * inc_fact

  ! Next, determine the hail size for this value of working_gwc:

  ! Determine slope and intercept parameters
  lamg = ( constp(100) / working_gwc ) ** cx(93)
  n0g  = x1g * ( lamg**x2g )

  ! Loop over the bins to try and find the bin which this value of
  ! graupel fits into.
  sum_n0g      = 0.0
  working_bin  = 0

  DO l = nbins, 1, -1
    diam_g     = hail_sizes(l)
    sum_in_bin = n0g * ( diam_g ** x4g ) * EXP( -1.0 * lamg * diam_g )        &
                 * bin_widths(l)

    sum_n0g = sum_n0g + sum_in_bin

    IF ( sum_n0g > threshold_concentration ) THEN
      working_bin = l
      EXIT ! This should exit the 'l' loop over nbins
    END IF

  END DO ! l

  IF (working_bin == previous_bin + 1) THEN
    ! GWC has moved to the next bin along. This is expected behaviour.
    ! Make this bin the new previous bin and set the current value
    ! of graupel water as the minimum value.
    previous_bin = working_bin

    gwc_min(working_bin) = working_gwc

  ELSE IF (previous_bin > 0 .AND. (previous_bin - working_bin) > 1) THEN

    ! Whoops - moved more than one bin at once. The only time
    ! when this would be acceptable is if no bins have been set
    ! (previous_bin == 0). In this case, the graupel increments
    ! are too big, causing the rapid jump between bins. If this is
    ! the case, reduce inc_fact. Throw an error message to alert
    ! the user to try this.

    ErrorStatus = 100 ! Flag to cause an issue

    emessage = 'Graupel increment too large when setting min '//newline//     &
               'values of GWC for lookup table in routine '   //newline//     &
               'generate_hail_bins. Try reducing the value ' //newline//      &
               'of inc_fact in generate_hail_bins '

    CALL ereport(RoutineName, ErrorStatus, emessage)

  END IF

  ! Finally, move the working_gwc on by the increment
  working_gwc = working_gwc + working_gwc_inc

END DO ! WHILE working_bin < nbins

WRITE(line_buffer, '(A,I6)')                                                  &
'Hail size look-up table found. Number of iterations =', iter_count
CALL umprint(line_buffer, src=RoutineName)

! Reverse the order so the hail_sizes array is largest first
! This makes it in the right direction for checking from the largest
! mass first.
hail_sizes(:) = hail_sizes(nbins:1:-1)
gwc_min(:)    = gwc_min(nbins:1:-1)

! Finally after this code has completed, set the flag for hail bins to be True.
l_hail_bins_set = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE generate_hail_bins
!==============================================================================
SUBROUTINE max_hail_size_2d(gwc_2d, max_hail_diam2d)

USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, INTENT(IN) :: gwc_2d( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end )
! An input value of graupel water content in [kg m-3]
! This can either be the maximum value in each column or for a specified
! level (e.g. the surface)

REAL, INTENT(OUT) :: max_hail_diam2d( tdims%i_start : tdims%i_end,            &
                                      tdims%j_start : tdims%j_end )
! An output value of hail size in units of [mm].


INTEGER :: i, j, l
! Loop counters

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='MAX_HAIL_SIZE_2D'

! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. l_hail_bins_set ) CALL generate_hail_bins()

DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    ! Initialise size to zero to ensure that something
    ! is always output
    max_hail_diam2d(i,j) = 0.0

    ! For the vast majority of grid boxes, there will be no graupel
    ! or hail present. Therefore having an IF test here will ensure
    ! that the maximum size is zero for no graupel and avoid lots
    ! of unnecessary IF tests in the 'l' loop below.
    IF ( gwc_2d(i,j) > gwc_abs ) THEN

      DO l = 1, nbins

        IF (gwc_2d(i,j) >= gwc_min(l) ) THEN
          max_hail_diam2d(i,j) = hail_sizes(l) * m2mm
          EXIT ! exits the 'l' DO loop
        END IF ! gwc_2d

      END DO ! l

    END IF ! gwc above mimimum

  END DO   ! i
END DO     ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE max_hail_size_2d
!==============================================================================
SUBROUTINE max_hail_size_3d(gwc_3d, max_hail_diam3d)

USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, INTENT(IN) :: gwc_3d( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end,                      &
                                        1 : tdims%k_end )

REAL, INTENT(OUT) :: max_hail_diam3d( tdims%i_start : tdims%i_end,            &
                                      tdims%j_start : tdims%j_end,            &
                                                  1 : tdims%k_end )

INTEGER :: i, j, k, l
! Loop counters

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='MAX_HAIL_SIZE_3D'

! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. l_hail_bins_set ) CALL generate_hail_bins()

DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Initialise size to zero to ensure that something
      ! is always output
      max_hail_diam3d(i,j,k) = 0.0

      ! For the vast majority of grid boxes, there will be no graupel
      ! or hail present. Therefore having an IF test here will ensure
      ! that the maximum size is zero for no graupel and avoid lots
      ! of unnecessary IF tests in the 'l' loop below.
      IF ( gwc_3d(i,j,k) > gwc_abs ) THEN

        DO l = 1, nbins

          IF (gwc_3d(i,j,k) >= gwc_min(l) ) THEN
            max_hail_diam3d(i,j,k) = hail_sizes(l) * m2mm
            EXIT ! exits the 'l' DO loop
          END IF ! gwc_3d

        END DO ! l

      END IF ! gwc_3d above minimum

    END DO   ! i
  END DO     ! j
END DO       ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE max_hail_size_3d
!==============================================================================
END MODULE max_hail_size_mod
