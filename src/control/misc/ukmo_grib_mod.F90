! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module UKMO_GRIB_MOD --------------------------------------
!
! Purpose : Read GRIB 1 & 2 files
!
! Coding Standard : UM documentation paper no. 3
!
! Documentation : None
!
!----------------------------------------------------------------
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
MODULE ukmo_grib_mod

USE missing_data_mod, ONLY: rmdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    nl => newline
USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE

PRIVATE

PUBLIC :: ukmo_decode, grib_api_computed_keys

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKMO_GRIB_MOD'

! Type to contain a selection of computed keys from 
! the GRIB API. Computed keys are not directly linked to 
! values in the GRIB message (as coded keys are) instead they 
! are calculated from other coded and computed keys.
TYPE grib_api_computed_keys
  INTEGER :: validityDate  ! Date of validity of the forecast 
  INTEGER :: validityTime  ! Time of validity of the forecast 
END TYPE grib_api_computed_keys

CONTAINS
SUBROUTINE ukmo_decode(fp_data, len_fp, num_fp,      &
                  vert_coords, len_vert, num_vert,   &
                  bitmap, len_bitmap, num_bitmap,    &
                  quasi, len_q, num_q,               &
                  block_0, block_1, block_2,         &
                  block_4, block_r,                  &
                  computed_keys,                     &
                  mesg, len_mesg)

#if defined(GRIB_API) 
USE grib_api
#endif

IMPLICIT NONE

! Scalar arguments
INTEGER, INTENT(IN) :: len_fp
INTEGER, INTENT(IN) :: len_vert
INTEGER, INTENT(IN) :: len_bitmap
INTEGER, INTENT(IN) :: len_q
INTEGER, INTENT(IN) :: len_mesg

! Array arguments
INTEGER, INTENT(IN) :: mesg(len_mesg)

! Array arguments
INTEGER, INTENT(INOUT) :: block_4(2)

! Scalar arguments
INTEGER, INTENT(OUT) :: num_fp
INTEGER, INTENT(OUT) :: num_vert
INTEGER, INTENT(OUT) :: num_bitmap
INTEGER, INTENT(OUT) :: num_q

! Array arguments
REAL,    INTENT(OUT) :: fp_data(len_fp)
REAL,    INTENT(OUT) :: vert_coords(len_vert)
INTEGER, INTENT(OUT) :: bitmap(len_bitmap)
INTEGER, INTENT(OUT) :: quasi(len_q)
INTEGER, INTENT(OUT) :: block_0(4)
INTEGER, INTENT(OUT) :: block_1(*)
INTEGER, INTENT(OUT) :: block_2(20)
REAL,    INTENT(OUT) :: block_r(20)

! Derived type arguments
TYPE(grib_api_computed_keys), INTENT(OUT) :: computed_keys

! Local variables

CHARACTER (LEN=*), PARAMETER     :: RoutineName='UKMO_DECODE'
CHARACTER (LEN=errormessagelength)   :: Cmessage      ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Main part of routine needs to be protected with GRIB API ifdef
#if defined(GRIB_API)

INTEGER(KIND=KindOfInt ) :: igrib_out
CHARACTER(LEN=1), ALLOCATABLE   :: message(:)
INTEGER(KIND=KindOfLong) :: intval
INTEGER(KIND=KindOfLong), ALLOCATABLE :: intarray(:)
REAL(KIND=KindOfDouble) :: realval
REAL(KIND=KindOfDouble), ALLOCATABLE :: realarray(:)
CHARACTER(LEN=32)       :: charval
INTEGER(KIND=KindOfInt)  :: missing_val

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! mesg is GRIB message
ALLOCATE(message(len_mesg*BIT_SIZE(mesg)/8))

message = TRANSFER(mesg, message)

CALL grib_new_from_message(igrib_out, message)

DEALLOCATE(message)

! block_0(1)  GRIB edition number
CALL grib_get(igrib_out, 'editionNumber', intval)
block_0(1) = intval
IF (block_0(1) /= 1) THEN
  Cmessage    = 'Attempting to convert GRIB2 to GRIB1'
  ErrorStatus = -1
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Lets set the level number before converting GRIB2 to GRIB1.
! block_1( 7) Zero if level type requires no parameters.
!             Height, pressure, etc. of level, if level type requires a single
!             descriptive parameter.  First height, pressure, etc. if level
!             requires two descriptive parameters.
! block_1( 8) Zero if level type requires no descriptive parameter or one
!             descriptive parameter.  Second height, pressure, etc. if level
!             requires two descriptive parameters.
CALL grib_get(igrib_out,'typeOfLevel',charval)
SELECT CASE (TRIM(charval))

CASE ('surface')
  block_1(7:8) = 0

CASE ('isobaricInhPa','heightAboveSea','heightAboveGround','depthBelowLand')
  CALL grib_get(igrib_out,'level',intval)
  block_1(7) = intval
  block_1(8) = 0

CASE ('hybrid','hybridHeight','hybridPressure')
  CALL grib_get(igrib_out,'level',intval)
  block_1(7) = intval
  block_1(8) = 0

CASE ('depthBelowLandLayer')
  CALL grib_get(igrib_out,'topLevel',intval)
  block_1(7) = intval

  ! Check if missing
  CALL grib_is_missing(igrib_out, 'bottomLevel',missing_val)
  IF (missing_val == 1) THEN
    ! It is missing
    intval = 289
    Cmessage    = 'bottomLevel is missing, setting to 289'
    ErrorStatus = -2
    CALL EReport( RoutineName, ErrorStatus, Cmessage )
  ELSE
    CALL grib_get(igrib_out,'bottomLevel',intval)
  END IF
  block_1(8) = intval

CASE DEFAULT
  Cmessage    = 'Level type not supported'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage )

END SELECT

! Lets make sure we store the length of the PV array.
! block_2( 1) Number of vertical coordinate parameters if GRIB 1, undefined
!             otherwise.
CALL grib_get(igrib_out,'PVPresent',intval)
IF (intval == 1) THEN
  CALL grib_get_size(igrib_out,'pv',intval)
  block_2(1) = intval
END IF

! Lets get the PV array before we try converting GRIB2 to GRIB1.
! num_vert    Number of elements in vert_coords.
CALL grib_get(igrib_out,'PVPresent',intval)
IF (intval == 1) THEN
  CALL grib_get_size(igrib_out,'pv',intval)
  num_vert = intval
ELSE
  num_vert = 0
END IF

IF (num_vert > 0 .AND. len_vert > num_vert) THEN
  ! vert_coords If message is hybrid.  Undefined for other types.
  ALLOCATE(realarray(num_vert))
  CALL grib_get(igrib_out,'pv',realarray)
  vert_coords(1:num_vert) = realarray(1:num_vert)
  DEALLOCATE(realarray)
ELSE IF (len_vert < num_vert) THEN
  Cmessage    = 'Vertical coordinate is too small'
  ErrorStatus = 20
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Now if we are not GRIB1 lets set some values to safe GRIB1 values.
IF (block_0(1) /= 1) THEN
  CALL grib_get(igrib_out,'typeOfLevel',charval)
  SELECT CASE (TRIM(charval))
  CASE ('hybrid','hybridHeight','hybridPressure')
    intval = 1
    CALL grib_set(igrib_out,'level',intval)
    realval = 0.0
    CALL grib_set(igrib_out,'pv',(/realval,realval/))
  END SELECT

  ! Now lets convert to GRIB1.
  intval = 1
  CALL grib_set(igrib_out,'editionNumber',intval)
END IF

! From this point on we should have GRIB1 message.

! block_0(2)  Tables version number if GRIB 1.  Undefined otherwise.
CALL grib_get(igrib_out, 'table2Version', intval)
block_0(2) = intval

! block_0(3)  Total length of message in octets, if GRIB 1.  Undefined
!             otherwise.
CALL grib_get(igrib_out, 'totalLength', intval)
block_0(3) = intval

! block_0(4)  Undefined.
block_0(4) = 0

! block_1( 1) Originating centre from GRIB code table C-1
CALL grib_get(igrib_out, 'centre', intval)
block_1(1) = intval

! block_1( 2) Generating process
CALL grib_get(igrib_out, 'generatingProcessIdentifier', intval)
block_1(2) = intval

! block_1( 3) Grid identification
CALL grib_get(igrib_out,'gridDefinition',intval)
block_1(3) = intval

! block_1( 4) Block indicator flags (GRIB table 1)
CALL grib_get(igrib_out,'section1Flags',intval)
block_1(4) = intval

! block_1( 5) Identification of parameter (GRIB table 2)
CALL grib_get(igrib_out,'indicatorOfParameter',intval)
block_1(5) = intval

! block_1( 6) Indicator of type of level (GRIB table 3)
CALL grib_get(igrib_out,'indicatorOfTypeOfLevel',intval)
block_1(6) = intval

! block_1( 9) Year of century of reference/data time
CALL grib_get(igrib_out,'yearOfCentury',intval)
block_1(9) = intval

! block_1(10) Month of reference/data time
CALL grib_get(igrib_out,'month',intval)
block_1(10) = intval

! block_1(11) Day of reference/data time
CALL grib_get(igrib_out,'day',intval)
block_1(11) = intval

! block_1(12) Hour of reference/data time
CALL grib_get(igrib_out,'hour',intval)
block_1(12) = intval

! block_1(13) Minute of reference time
CALL grib_get(igrib_out,'minute',intval)
block_1(13) = intval

! block_1(14) Indicator of time unit, as defined in GRIB table 4
CALL grib_get(igrib_out,'unitOfTimeRange',intval)
block_1(14) = intval

! block_1(15) First time interval (P1), if required for the type of time range
!             used. Undefined otherwise.
CALL grib_get(igrib_out,'P1',intval)
block_1(15) = intval

! block_1(16) Second time interval (P2), if required for the type of time range
!             used.  Undefined otherwise.
CALL grib_get(igrib_out,'P2',intval)
block_1(16) = intval

! block_1(17) Time range indicator, as defined in GRIB table 5
CALL grib_get(igrib_out,'timeRangeIndicator',intval)
block_1(17) = intval

! block_1(18) Number included in average or accumulation, if required for the
!             type of time range used.  Zero otherwise.
CALL grib_get(igrib_out,'numberIncludedInAverage',intval)
block_1(18) = intval

! block_1(19) Century of reference time, if GRIB 1.  Undefined otherwise.
CALL grib_get(igrib_out,'centuryOfReferenceTimeOfData',intval)
block_1(19) = intval

! block_1(20) Decimal scale factor, if GRIB 1. Undefined otherwise.
CALL grib_get(igrib_out, 'decimalScaleFactor', intval)
block_1(20) = intval

! block_1(21) Length of block 1 in octets, if GRIB 1, Undefined otherwise.
block_1(21) = 34

! block_1(22..33) Undefined (reserved for future use)
block_1(22:33) = 0

! block_1(34) Sub centre number ID
CALL grib_get(igrib_out, 'subCentre', intval)
block_1(34) = intval

! block_1(35..) Centre defined octets in each element.  Length is (length of
! block 1 - 40).  If centre is 74, sub-centre 2, then super STASH pp info.
! Not set if not used since size of array is unknown.

! block_2( 2) PV or PL if GRIB 1 and hybrid/quasi regular.  255 otherwise or
!             undefined for other editions.
CALL grib_get(igrib_out,'pvlLocation',intval)
block_2(2) = intval

! block_2( 3) Data representation type, GRIB table 6.
CALL grib_get(igrib_out,'dataRepresentationType',intval)
block_2(3) = intval

CALL grib_get(igrib_out,'gridType',charval)
SELECT CASE(TRIM(charval))

  ! For lat/long grid
CASE ('regular_ll')
  ! block_2( 4) Number of points along a parallel.
  CALL grib_get(igrib_out,'Ni',intval)
  block_2(4) = intval

  ! block_2( 5) Number of points along a meridian.
  CALL grib_get(igrib_out,'Nj',intval)
  block_2(5) = intval

  ! block_2( 6) Latitude of first grid point.
  CALL grib_get(igrib_out,'latitudeOfFirstGridPoint',intval)
  block_2(6) = intval

  ! block_2( 7) Longitude of first grid point.
  CALL grib_get(igrib_out,'longitudeOfFirstGridPoint',intval)
  block_2(7) = intval

  ! blocK_2( 8) Resolution/component flag, GRIB table 7.
  CALL grib_get(igrib_out,'resolutionAndComponentFlags',intval)
  block_2(8) = intval

  ! block_2( 9) Latitude of extreme point.
  CALL grib_get(igrib_out,'latitudeOfLastGridPoint',intval)
  block_2(9) = intval

  ! block_2(10) Longitude of extreme point.
  CALL grib_get(igrib_out,'longitudeOfLastGridPoint',intval)
  block_2(10) = intval

  ! block_2(11) Increment along a parallel.
  CALL grib_get(igrib_out,'jDirectionIncrement',intval)
  block_2(11) = intval

  ! block_2(12) Increment along a meridian.
  CALL grib_get(igrib_out,'iDirectionIncrement',intval)
  block_2(12) = intval

  ! block_2(13) Scanning mode flags, as defined in GRIB table 8.
  CALL grib_get(igrib_out,'scanningMode',intval)
  block_2(13) = intval

  ! block_2(14) Latin, if grid type is Mercator, or latitude of southern pole if
  !             rotated pole. Undefined for other GRIB types.
  block_2(14) = 0.0

  ! block_2(15) Longitude of southern pole if rotated pole.  Undefined otherwise.
  block_2(15) = 0.0

  ! block_2(16..20) Undefined.

  ! For rotated lat/long grid
CASE ('rotated_ll')
  ! block_2( 4) Number of points along a parallel.
  CALL grib_get(igrib_out,'Ni',intval)
  block_2(4) = intval

  ! block_2( 5) Number of points along a meridian.
  CALL grib_get(igrib_out,'Nj',intval)
  block_2(5) = intval

  ! block_2( 6) Latitude of first grid point.
  CALL grib_get(igrib_out,'latitudeOfFirstGridPoint',intval)
  block_2(6) = intval

  ! block_2( 7) Longitude of first grid point.
  CALL grib_get(igrib_out,'longitudeOfFirstGridPoint',intval)
  block_2(7) = intval

  ! blocK_2( 8) Resolution/component flag, GRIB table 7.
  CALL grib_get(igrib_out,'resolutionAndComponentFlags',intval)
  block_2(8) = intval

  ! block_2( 9) Latitude of extreme point.
  CALL grib_get(igrib_out,'latitudeOfLastGridPoint',intval)
  block_2(9) = intval

  ! block_2(10) Longitude of extreme point.
  CALL grib_get(igrib_out,'longitudeOfLastGridPoint',intval)
  block_2(10) = intval

  ! block_2(11) Increment along a parallel.
  CALL grib_get(igrib_out,'jDirectionIncrement',intval)
  block_2(11) = intval

  ! block_2(12) Increment along a meridian.
  CALL grib_get(igrib_out,'iDirectionIncrement',intval)
  block_2(12) = intval

  ! block_2(13) Scanning mode flags, as defined in GRIB table 8.
  CALL grib_get(igrib_out,'scanningMode',intval)
  block_2(13) = intval

  ! block_2(14) Latin, if grid type is Mercator, or latitude of southern pole if
  !             rotated pole.  Undefined for other GRIB types.
  CALL grib_get(igrib_out,'latitudeOfSouthernPole',intval)
  block_2(14) = intval

  ! block_2(15) Longitude of southern pole if rotated pole.  Undefined otherwise.
  CALL grib_get(igrib_out,'longitudeOfSouthernPole',intval)
  block_2(15) = intval

  ! block_2(16..20) Undefined.

CASE DEFAULT
  Cmessage    = 'Grid type not supported'
  ErrorStatus = 30
  CALL EReport( RoutineName, ErrorStatus, Cmessage )

END SELECT


! block_3     Undefined.
!
! block_4( 1) Flag (packing type), if GRIB 1.  Undefined otherwise.
CALL grib_get(igrib_out,'dataFlag',intval)
block_4(1) = intval

! block_4( 2) Input value, unchanged.
! Do nothing

! block_r( 1) Angle of rotation, if grid type rotated.  Undefined otherwise.

CALL grib_get(igrib_out,'gridType',charval)
IF (TRIM(charval) == 'rotated_ll') THEN
  CALL grib_get(igrib_out,'angleOfRotationInDegrees',realval)
  block_r(1) = realval
END IF

! block_r( 2) Missing data indicator, if bitmapped grid used.  Undefined
!             otherwise.
! This is up to the user to define.
CALL grib_get(igrib_out,'bitmapPresent',intval)
IF (intval == 1) THEN
  realval = rmdi
  block_r(2) = realval
  ! Set the missing value.
  CALL grib_set(igrib_out,'missingValue',realval)
END IF

! block_r(3..20) Undefined.

! Get computed keys

! Validity/forecast date
CALL grib_get(igrib_out, 'validityDate', intval)
computed_keys%validityDate = intval

! Validity/forecast time
CALL grib_get(igrib_out, 'validityTime', intval)
computed_keys%validityTime = intval

! num_fp      Number of elements of fp_data
CALL grib_get(igrib_out,'numberOfPoints',intval)
num_fp = intval
IF (num_fp > len_fp) THEN
  WRITE(umMessage,'(A,I10,A,I10)') 'len_fp is ', len_fp,  &
                      ' but numberOfPoints is ', num_fp
  CALL umPrint(umMessage,src='ukmo_grib_mod')
  Cmessage    = 'fp_data too small'
  ErrorStatus = 40
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! fp_data     Field of data contained in message
ALLOCATE(realarray(num_fp))
CALL grib_get(igrib_out,'values',realarray);

! fp_data contains full field.  missingValue above should have set the mdi.
fp_data(1:num_fp) = realarray(1:num_fp)

DEALLOCATE(realarray)

! num_bitmap  Number of elements of bitmap.
CALL grib_get(igrib_out,'bitmapPresent',intval)
IF (intval == 1) THEN
  num_bitmap = num_fp
ELSE
  num_bitmap = 0
END IF

IF (num_bitmap > 0 .AND. len_bitmap > num_bitmap) THEN
  ! bitmap      If bitmap is used.
  bitmap(1:num_fp) = MERGE(1, 0, fp_data(1:num_fp) == block_r(2))
ELSE IF (len_bitmap < num_bitmap) THEN
  WRITE(umMessage,'(A,I10,A,I10)') 'len_bitmap is ', len_bitmap,  &
                              ' but num_bitmap is ', num_bitmap
  CALL umPrint(umMessage,src='ukmo_grib_mod')
  Cmessage    = 'Bitmap array is too small'
  ErrorStatus = 40
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Quasi-regular grid unsupported.
! num_q       Number of elements in quasi
num_q = 0

! quasi       Number of points along each parallel.
quasi = 0

! Release memory
CALL grib_release(igrib_out)


#else
! GRIB_API ifdef
Cmessage = 'Attempting to read GRIB data, but GRIB_API ifdef not present'//nl//&
     ' in build. Please note the GRIB API is currently the only method'  //nl//&
     ' for decoding GRIB'
ErrorStatus = 100
CALL EReport( RoutineName, ErrorStatus, Cmessage )
#endif
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukmo_decode
END MODULE ukmo_grib_mod
