! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Module containing a routine which processes filename strings
!              containing special control characters and replaces these with
!              date/time information
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!             This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!
MODULE filename_generation_mod

USE filenamelength_mod, ONLY: filenamelength
USE submodel_mod, ONLY: atmos_im

IMPLICIT NONE

PRIVATE

! Mappings for month and season names
CHARACTER(LEN=3), PARAMETER :: &
    month_names_3char(12)  = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', &
                              'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
CHARACTER(LEN=3), PARAMETER :: &
    season_names_3char(12) = ['jfm', 'fma', 'mam', 'amj', 'mjj', 'jja', &
                              'jas', 'aso', 'son', 'ond', 'ndj', 'djf']

! Arrays used to store previously generated filenames
CHARACTER(LEN=filenamelength), ALLOCATABLE, SAVE :: &
    used_filenames(:)
CHARACTER(LEN=filenamelength), ALLOCATABLE :: &
    temp_used_filenames(:)

PUBLIC :: get_filename

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FILENAME_GENERATION_MOD'

CONTAINS

SUBROUTINE get_filename(base_filename, filename, absolute_offset_steps, &
                        relative_offset_steps, reinit_steps, meanlev, dump)

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE filenamelength_mod,   ONLY: filenamelength
USE nlstcall_mod,         ONLY: model_basis_time, model_analysis_mins, lcal360
USE umprintmgr,           ONLY: umprint, ummessage, newline
USE ereport_mod,          ONLY: ereport
USE conversions_mod,      ONLY: isec_per_day

USE model_time_mod, ONLY:                               &
    i_second, i_minute, i_hour, i_day, i_month, i_year, &
    i_day_number, secs_per_stepim, stepim

IMPLICIT NONE

CHARACTER(LEN=*), &
         INTENT(IN)  :: base_filename   ! Base part of filename to generate
CHARACTER(LEN=filenamelength), &
         INTENT(OUT) :: filename        ! Generated file name

! If provided, this provides an offset in timesteps which any absolute time
! in the generated name should be adjusted by (if the name is generated at an
! awkward point such as at the end of a meaning period)
INTEGER, INTENT(IN), OPTIONAL :: absolute_offset_steps
! If provided, this provides an offset in timesteps which any relative time
! in the generated name should be adjusted by (to allow the reference time to
! be different to the basis time, typically when using an analysis)
INTEGER, INTENT(IN), OPTIONAL :: relative_offset_steps
! The frequency between re-initialisations in whole timesteps; this is only
! needed if using the Legacy style absolute naming for non-mean output files
INTEGER, INTENT(IN), OPTIONAL :: reinit_steps
! The meaning period level.  This is only needed if using the Legacy style
! absolute naming for climate mean output files
INTEGER, INTENT(IN), OPTIONAL :: meanlev
! This flag is included for legacy compatibility to deal with the slightly
! altered method for naming dump names as part of the absolute naming
LOGICAL, INTENT(IN), OPTIONAL :: dump

INTEGER :: mean_level ! Meaning level
INTEGER :: init_steps ! Re-initialising step frequency

! Total time relative to the analysis recorded in each unit (not the components
! i.e. total_seconds will be the same as total_minutes*60, total_hours*3600,etc)
! If no analysis is in use this will instead relate to the basis time
INTEGER :: total_steps
INTEGER :: total_seconds
INTEGER :: total_minutes
INTEGER :: total_hours
INTEGER :: total_days
INTEGER :: total_months
INTEGER :: total_years

! Components of total time relative to the analysis (not the total times as
! above, i.e. the total time is rel_years, rel_months, ... , rel_minutes and
! rel_seconds). If no analysis is in use this will relate to the basis time
INTEGER :: rel_seconds
INTEGER :: rel_minutes
INTEGER :: rel_hours
INTEGER :: rel_days
INTEGER :: rel_months
INTEGER :: rel_years

! Components of the analysis time (or the basis time if no analysis is in use)
INTEGER :: analysis_seconds
INTEGER :: analysis_minutes
INTEGER :: analysis_hours
INTEGER :: analysis_days
INTEGER :: analysis_months
INTEGER :: analysis_years
INTEGER :: analysis_dayno

! Components of the absolute time (the exact calendar date/time)
INTEGER :: abs_seconds
INTEGER :: abs_minutes
INTEGER :: abs_hours
INTEGER :: abs_days
INTEGER :: abs_months
INTEGER :: abs_years
INTEGER :: abs_dayno

! Working variables to aid calculation of the above
INTEGER :: total_analysis_days
INTEGER :: total_analysis_secs
INTEGER :: total_abs_days
INTEGER :: total_abs_secs

! Variables to store a copy of the selected value for output into the filename
INTEGER :: seconds
INTEGER :: minutes
INTEGER :: hours
INTEGER :: days
INTEGER :: months
INTEGER :: years

INTEGER :: month           ! Month index for season/month name array
INTEGER :: i_char          ! Loop counter
INTEGER :: i_file          ! Loop counter
INTEGER :: icode           ! Error status indicator
LOGICAL :: relative        ! Switch to indicate relative time selected
LOGICAL :: relative_total  ! Switch to indicate total time selected
LOGICAL :: analysis        ! Switch to indicate analysis time selected
LOGICAL :: after_marker    ! Switch to indicate if marker was encountered
LOGICAL :: legacy_dumpmode ! Switch to indicate legacy dump naming

INTEGER, PARAMETER :: max_suffix_len = 64  ! Length of below suffix
CHARACTER(LEN=max_suffix_len)  &           ! Suffix to be appended to the name
                   :: filename_suffix      ! after each pass of the loop

CHARACTER(LEN=1), &                  ! Marker character used to denote special
    PARAMETER     :: marker = "%"    ! special control sequences in template
CHARACTER(LEN=1)  :: chr             ! Character of filename in loop
CHARACTER(LEN=1)  :: rel_sign        ! Sign to indicate if relative time is
                                     ! before or after the analysis time

CHARACTER(LEN=10) :: pad      ! Padding
INTEGER           :: int_pad  ! Padding as an integer for comparisons

CHARACTER(LEN=*), PARAMETER :: RoutineName="GET_FILENAME"

! Dr-Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Calculate "total_" variables; these give the total amount of time relative
! to the start of the model run or analysis time (if used) in each unit
! ------------------------------------------------------------------------------

! Get the total number of timesteps (and adjust if an offset was provided)
total_steps = stepim(atmos_im)
IF (PRESENT(relative_offset_steps)) THEN
  total_steps = total_steps + relative_offset_steps
END IF

total_seconds = total_steps * secs_per_stepim(atmos_im)

! Negative relative time means we are currently before the analysis time; set up
! a special character to indicate this before making the seconds positive
IF (total_seconds < 0) THEN
  total_seconds = ABS(total_seconds)
  rel_sign = "-"
ELSE
  rel_sign = "+"
END IF

total_minutes = total_seconds/60
total_hours   = total_minutes/60
total_days    = total_hours/24
IF (lcal360) THEN
  total_months = total_days/30
  total_years  = total_months/12
ELSE
  total_months = 0
  total_years  = 0
END IF

! Calculate "rel_" variables; these give the components of the time relative
! to the start of the model run or analysis time (if used).
! ------------------------------------------------------------------------------
IF (lcal360) THEN
  rel_years  = total_years
  rel_months = total_months - total_years*12
  rel_days = total_days - total_months*30
ELSE
  rel_months = 0
  rel_years  = 0
  rel_days = total_days
END IF
rel_hours = total_hours - total_days*24
rel_minutes = total_minutes - total_hours*60
rel_seconds = total_seconds - total_minutes*60

! Calculate analysis variables; these give the components of the analysis time
! if model_analysis_mins wasn't set this is equal to the basis time
!-------------------------------------------------------------------------------
! First get the basis time in whole days + seconds since calendar ref. date
! DEPENDS ON: time2sec
CALL time2sec(model_basis_time(1), model_basis_time(2),             &
              model_basis_time(3), model_basis_time(4),             &
              model_basis_time(5), model_basis_time(6),             &
              0, 0, total_analysis_days, total_analysis_secs, lcal360)

! Then adjust this based on the analysis time, allowing for day boundary
total_analysis_secs = total_analysis_secs + model_analysis_mins*60
DO WHILE (total_analysis_secs >= isec_per_day)
  total_analysis_days = total_analysis_days + 1
  total_analysis_secs = total_analysis_secs - isec_per_day
END DO

! DEPENDS ON: sec2time
CALL sec2time(0, 0, total_analysis_days, total_analysis_secs, &
              analysis_years,   analysis_months,              &
              analysis_days,    analysis_hours,               &
              analysis_minutes, analysis_seconds,             &
              analysis_dayno, &
              lcal360)

! Calculate absolute variables; these give the exact date/time of the output,
! but may be offset if required
! ------------------------------------------------------------------------------
IF (PRESENT(absolute_offset_steps)) THEN
  ! First get the absolute time in whole days + seconds since calendar ref. date
  ! DEPENDS ON: time2sec
  CALL time2sec(i_year, i_month, i_day, i_hour, i_minute, i_second,            &
                0, 0, total_abs_days, total_abs_secs, lcal360)

  ! Then adjust this based on the given offset, allowing for day boundary
  total_abs_secs = total_abs_secs                                              &
                 + absolute_offset_steps*secs_per_stepim(atmos_im)

  DO WHILE (total_abs_secs >= isec_per_day)
    total_abs_days = total_abs_days + 1
    total_abs_secs = total_abs_secs - isec_per_day
  END DO

  ! DEPENDS ON: sec2time
  CALL sec2time(0, 0, total_abs_days, total_abs_secs,           &
                abs_years, abs_months, abs_days, abs_hours,     &
                abs_minutes, abs_seconds, abs_dayno, lcal360)
ELSE
  abs_years   = i_year
  abs_months  = i_month
  abs_days    = i_day
  abs_hours   = i_hour
  abs_minutes = i_minute
  abs_seconds = i_second
END IF

! Iterate through the filename template and convert specially marked control
! characters with equivalent name information
! ------------------------------------------------------------------------------
!
! Initialise the output filename
filename = ""
after_marker = .FALSE.

! Iterate through each character in the base filename
DO i_char = 1, LEN_TRIM(base_filename)
  chr = base_filename(i_char:i_char)

  ! If the marker character is encountered, activate the flag and go onto the
  ! next character for processing, also initialise/reset the padding string and
  ! other indicator flags
  IF (chr == marker .AND. .NOT. after_marker) THEN
    after_marker   = .TRUE.
    relative       = .FALSE.
    relative_total = .FALSE.
    analysis       = .FALSE.
    pad            = ""
    CYCLE
  END IF

  ! The suffix is stored in the below variable and appended to the generated
  ! filename after each pass
  filename_suffix = ""

  ! If the flag was set on the previous character then this character is a
  ! special character to be interpreted
  IF (after_marker) THEN
    ! First switch off the flag for the following character
    after_marker = .FALSE.
    SELECT CASE(chr)

      CASE(marker)
        ! If the marker is followed by another marker symbol, this is
        ! an escape to print the literal marker symbol
        filename_suffix = chr

      CASE("A")
        ! This specifier indicates that the value should be taken from the
        ! analysis time (or the basis time if not analysis is in use) instead
        ! of the validity time
        after_marker = .TRUE.
        analysis     = .TRUE.

      CASE("Q")
        ! This specifier indicates the start of a "relative" value, and will be
        ! followed by (optionally) some padding digits (see below), and one of
        ! the other specifiers above - the difference between this and the "R"
        ! specifier below is that this returns the relative time component
        ! in the given unit (i.e. RH gives the total relative hours from the
        ! analysis/basis but QH gives only the hour component of the relative
        ! time from the analysis/basis)
        after_marker = .TRUE.
        relative     = .TRUE.

      CASE("R")
        ! This specifier indicates the start of a "relative total" value, and
        ! will be followed by (optionally) some padding digits (see below), and
        ! one of the other specifiers above - the difference between this and
        ! the "Q" specifier above is that this returns the *total* relative time
        ! in the given unit (i.e. RH gives the total relative hours from the
        ! analysis/basis but QH gives only the hour component of the relative
        ! time from the analysis/basis)
        after_marker   = .TRUE.
        relative_total = .TRUE.

      CASE("0","1","2","3","4","5","6","7","8","9")
        ! If we encounter a number, ensure we carry on reading the next char
        after_marker = .TRUE.
        ! And if it followed a "R" or "Q" (relative) specifier, begin recording
        ! it to the padding string (otherwise skip it, as padding isn't needed
        ! for the absolute specifiers).
        IF (relative .OR. relative_total) THEN
          pad = TRIM(pad) // chr
          READ(pad, "(I10)") int_pad
          IF (int_pad > max_suffix_len) THEN
            icode = 1
            CALL ereport(routinename, icode,                      &
                "Maximum padding size (64) exceeded:"//  newline//&
                TRIM(filename))
          END IF
        END IF

      CASE("Y")
        ! 4 digit year
        !-------------
        IF (relative_total) THEN
          ! Use total years relative to analysis/basis time
          years = total_years
        ELSE IF (relative) THEN
          ! Use the years component of the relative time
          years = rel_years
        ELSE IF (analysis) THEN
          ! Use the years component of the analysis/basis time
          years = analysis_years
        ELSE
          ! Otherwise use the actual validity year
          years = abs_years
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "4"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((years / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") years

      CASE("m")
        ! 2 digit month number
        !---------------------
        IF (relative_total) THEN
          ! Use total months relative to analysis/basis time
          months = total_months
        ELSE IF (relative) THEN
          ! Use the months component of the relative time
          months = rel_months
        ELSE IF (analysis) THEN
          ! Use the months component of the analysis/basis time
          months = analysis_months
        ELSE
          ! Otherwise use the actual validity month
          months = abs_months
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "2"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((months / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") months

      CASE("d")
        ! 2 digit day number
        !---------------------
        IF (relative_total) THEN
          ! Use total days relative to analysis/basis time
          days = total_days
        ELSE IF (relative) THEN
          ! Use the days component of the relative time
          days = rel_days
        ELSE IF (analysis) THEN
          ! Use the days component of the analysis/basis time
          days = analysis_days
        ELSE
          ! Otherwise use the actual validity day
          days = abs_days
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "2"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((days / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") days

      CASE("H")
        ! Hour using 24-hour clock
        !-------------------------
        IF (relative_total) THEN
          ! Use total hours relative to analysis/basis time
          hours = total_hours
        ELSE IF (relative) THEN
          ! Use the hours component of the relative time
          hours = rel_hours
        ELSE IF (analysis) THEN
          ! Use the hours component of the analysis/basis time
          hours = analysis_hours
        ELSE
          ! Otherwise use the actual validity hour
          hours = abs_hours
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "2"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((hours / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") hours

      CASE("M")
        ! 2 digit minutes
        !----------------
        IF (relative_total) THEN
          ! Use total minutes relative to analysis/basis time
          minutes = total_minutes
        ELSE IF (relative) THEN
          ! Use the minutes component of the relative time
          minutes = rel_minutes
        ELSE IF (analysis) THEN
          ! Use the minutes component of the analysis/basis time
          minutes = analysis_minutes
        ELSE
          ! Otherwise use the actual validity minutes
          minutes = abs_minutes
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "2"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((minutes / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") minutes

      CASE("S")
        ! 2 digit seconds
        !----------------
        IF (relative_total) THEN
          ! Use total seconds relative to analysis/basis time
          seconds = total_seconds
        ELSE IF (relative) THEN
          ! Use the seconds component of the relative time
          seconds = rel_seconds
        ELSE IF (analysis) THEN
          ! Use the seconds component of the analysis/basis time
          seconds = analysis_seconds
        ELSE
          ! Otherwise use the actual validity seconds
          seconds = abs_seconds
        END IF

        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "2"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((seconds / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") seconds

      CASE("T")
        ! Timestep
        !---------
        ! If no padding specified, set to the default for this specifier
        IF (TRIM(pad) == "") pad = "4"

        ! Detect if the final padding size will fail to contain the value,
        ! and if it will reset the padding to zero (which will make it expand
        ! to fit the value exactly)
        READ(pad, "(I10)") int_pad
        IF ((total_steps / (10**int_pad)) /= 0) THEN
          pad = "0"
        END IF

        ! Write the result into the filename suffix
        WRITE(filename_suffix, "(i"//TRIM(pad)//"."//TRIM(pad)//")") total_steps

      CASE("b")
        ! 3 Character month name
        !-----------------------
        IF (analysis) THEN
          ! Use the month component of the analysis/basis time
          month = analysis_months
        ELSE
          ! Otherwise use the actual validity month
          month = abs_months
        END IF

        ! The suffix is taken from the array of month names
        filename_suffix = month_names_3char(month)

      CASE("s")
        ! 3 Character season_name
        !------------------------
        IF (analysis) THEN
          ! Use the month component of the analysis/basis time
          month = analysis_months
        ELSE
          ! Otherwise use the actual validity month
          month = abs_months
        END IF

        ! The suffix is taken from the array of season names
        filename_suffix = season_names_3char(month)

      CASE("N")
        ! Legacy option #1: reproduces the old "Relative" naming convention
        ! for hourly output.  This gives the hours since the basis time up
        ! to a maximum of 3600 hours, but strictly fit into 3-characters by
        ! making the "hundreds" digit use a base-36 notation with the letters
        ! of the alphabet (so once an hour of 999 is reached the next value
        ! is a00 (1,000) then a01 (1,001) up to the maximum z99 (3,599))

        ! The legacy relative naming did not correctly apply the offset for
        ! output stream files, so we will have to undo any offset applied
        ! already here, but not for dump files (the effect of this is that
        ! output streams will always be relative to the *basis* time, never
        ! the analysis)
        legacy_dumpmode = .FALSE.
        IF (PRESENT(dump)) legacy_dumpmode = dump

        IF (PRESENT(relative_offset_steps) .AND. (.NOT. legacy_dumpmode)) THEN
          filename_suffix = legacy_relative_naming(                           &
              INT((total_steps - relative_offset_steps)*                      &
                                   secs_per_stepim(atmos_im)/3600))
        ELSE
          filename_suffix = legacy_relative_naming(total_hours)
        END IF

      CASE("n")
        ! Legacy option #2: reproduces the old "Sub Hourly" naming convention
        ! for sub-hourly output.  This is similar to the above "N" option but
        ! is in a unit of hours and minutes. This goes up to the same maximum
        ! of 3600 hours using 4-characters, with the first digit in the same
        ! base-36 notation using the letters of the alphabet (so once
        ! 9959 (99 hours, 59 minutes) is reached the next value is
        ! a000 (100 hours, 0 minutes) then a001 (100 hours, 1 minute) up to
        ! the maximum value z959 (359 hours, 59 minutes))

        ! The legacy relative naming did not correctly apply the offset for
        ! output stream files, so we will have to undo any offset applied
        ! already here, but not for dump files (the effect of this is that
        ! output streams will always be relative to the *basis* time, never
        ! the analysis)
        legacy_dumpmode = .FALSE.
        IF (PRESENT(dump)) legacy_dumpmode = dump

        IF (PRESENT(relative_offset_steps) .AND. (.NOT. legacy_dumpmode)) THEN
          filename_suffix = legacy_subhourly_naming(                          &
              INT( (total_steps - relative_offset_steps) *                    &
                   secs_per_stepim(atmos_im) ))
        ELSE
          filename_suffix = legacy_subhourly_naming(total_seconds)
        END IF

      CASE("z")
        ! Legacy option #3: similar to the new "+" option but using the old
        ! convention that relative times after the analysis would use the letter
        ! "a" and names before the analysis time would use "z"
        IF (rel_sign == "+") THEN
          filename_suffix = "a"
        ELSE
          filename_suffix = "z"
        END IF

      CASE("C")
        ! Legacy option #4: reproduces the old "Absolute" naming conventions
        ! for climate style files.  This attempts to guess what the name should
        ! be based on the periodicity of the output frequency.  For example if
        ! the frequency is a multiple of 30 days it'll go for a name of the
        ! form "yyyy{mon}" where mon is the name of the *previous* month, if
        ! it fits into 90 days it'll be "yyyy{ssn}" where ssn is the name of
        ! the *previous* season, and other periods are labelled with a special
        ! letter or number followed by the date in a yyyymmdd[_hh] form

        ! This naming style treats dumps slightly differently to other types
        ! of file, so check if the optional flag was provided
        legacy_dumpmode = .FALSE.
        IF (PRESENT(dump)) legacy_dumpmode = dump

        ! This style also requires the mean period level index
        mean_level = 0
        IF (PRESENT(meanlev)) mean_level = meanlev

        ! And the reinit frequency
        init_steps = 0
        IF (PRESENT(reinit_steps)) init_steps = reinit_steps

        filename_suffix =                                                &
                   TRIM(legacy_absolute_naming(mean_level, init_steps,   &
                                               legacy_dumpmode))

      CASE("+")
        ! This specifier will be evaluated to a plus or minus sign depending
        ! on if the current time is before or after the basis time
        filename_suffix = rel_sign

      CASE DEFAULT
        ! If the character didn't match the above just reprint it including
        ! the marker to make it obvious which part of the name didn't expand
        filename_suffix = marker // chr

    END SELECT
  ELSE
    ! If this character wasn't a special marked character just add it to the
    ! output name unmodified
    filename_suffix = chr
  END IF

  ! After exiting the above the suffix contains the result of any expansions.
  ! To prevent output overflows test the length here and abort if it exceeds
  ! the filename size
  IF (LEN_TRIM(filename) + LEN_TRIM(filename_suffix) > filenamelength) THEN
    icode = 1
    CALL ereport(routinename, icode,                                          &
        "Generated filename has exceeded maximum allowed size:"//    newline//&
        TRIM(filename))
  END IF

  ! Otherwise append the suffix and move onto the next character
  filename = TRIM(filename) // TRIM(filename_suffix)

END DO  ! Loop over characters in input filename

! Print a message explaining the results of the generation
WRITE(ummessage, "(A)") &
    routinename//": Generated filename:"//TRIM(filename)
CALL umprint(ummessage,src=routinename)
WRITE(ummessage, "(A)") &
    routinename//":             (From): "//TRIM(base_filename)
CALL umprint(ummessage,src=routinename)

! Now ensure that the generated name hasn't been generated before
IF (ALLOCATED(used_filenames)) THEN
  ! Check to see if the filename clashes with any previously generated filename
  DO i_file = 1, SIZE(used_filenames)
    IF (filename == used_filenames(i_file)) THEN
      icode = 1
      CALL ereport(routinename, icode,                                        &
          "Filename clashes with previously generated name:"//       newline//&
          TRIM(filename))
    END IF
  END DO
  ! Add the filename to the list of previously generated filenames
  ALLOCATE(temp_used_filenames(SIZE(used_filenames)+1))
  temp_used_filenames(1:SIZE(used_filenames)) = used_filenames
  DEALLOCATE(used_filenames)
  CALL MOVE_ALLOC(temp_used_filenames, used_filenames)
  used_filenames(SIZE(used_filenames)) = filename
ELSE
  ALLOCATE(used_filenames(1))
  used_filenames = filename
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

END SUBROUTINE get_filename

!-------------------------------------------------------------------------------
FUNCTION legacy_relative_naming(total_hours) RESULT(extension)

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER          :: total_hours
CHARACTER(LEN=3) :: extension
CHARACTER(LEN=1), &
    PARAMETER    :: char_id(36) = ['0','1','2','3','4','5','6','7','8','9',   &
                                   'a','b','c','d','e','f','g','h','i','j',   &
                                   'k','l','m','n','o','p','q','r','s','t',   &
                                   'u', 'v', 'w', 'x', 'y', 'z']
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode

IF (total_hours >= 3600) THEN
  icode = 1
  WRITE(cmessage, '(A,I0,A)') &
      'Cannot use legacy naming with hours = ', total_hours, ', must be < 3600'
  CALL ereport('filename_generation:legacy_relative_naming', icode, cmessage)
ELSE
  extension = char_id(MOD(total_hours/100,36)+1) // &
              char_id(MOD(total_hours/10 ,10)+1) // &
              char_id(MOD(total_hours    ,10)+1)
END IF

END FUNCTION legacy_relative_naming

!-------------------------------------------------------------------------------
FUNCTION legacy_subhourly_naming(total_seconds) RESULT(extension)

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER          :: total_seconds
CHARACTER(LEN=4) :: extension
INTEGER          :: hours_and_mins
CHARACTER(LEN=1), &
    PARAMETER    :: char_id(36) = ['0','1','2','3','4','5','6','7','8','9',   &
                                   'a','b','c','d','e','f','g','h','i','j',   &
                                   'k','l','m','n','o','p','q','r','s','t',   &
                                   'u', 'v', 'w', 'x', 'y', 'z']
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode

! Turn the seconds into hours multiplied by 100
hours_and_mins = 100*(REAL(total_seconds)/3600.0)
! Now convert the decimal part of this into base 60 (minutes)
hours_and_mins = 100*(hours_and_mins/100) + &
                 NINT((REAL(MOD(hours_and_mins,100))/100.0)*60)

IF (hours_and_mins >= 36000) THEN
  icode = 1
  WRITE(cmessage, '(A,I0,A)') &
      'Cannot use legacy naming with seconds = ', total_seconds, &
      ', must be < 216000 (3600 hours)'
  CALL ereport('filename_generation:legacy_relative_naming', icode, cmessage)
ELSE

  extension = char_id(MOD(hours_and_mins/1000,36)+1) // &
              char_id(MOD(hours_and_mins/100 ,10)+1) // &
              char_id(MOD(hours_and_mins/10  ,10)+1) // &
              char_id(MOD(hours_and_mins     ,10)+1)
END IF

END FUNCTION legacy_subhourly_naming

!-------------------------------------------------------------------------------
FUNCTION legacy_absolute_naming(mean_period, reinit_steps, dump) &
                                RESULT(extension)

USE model_time_mod, ONLY: i_day, i_hour, i_month, i_year
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim, &
                       meanfreqim, dumpfreqim

IMPLICIT NONE

INTEGER :: mean_period
INTEGER :: reinit_steps
LOGICAL :: dump
CHARACTER(LEN=12) :: extension

INTEGER :: steps_per_hour
INTEGER :: i
INTEGER :: hours
INTEGER :: year
INTEGER :: month
CHARACTER(LEN=1) :: special_char

extension = ""

steps_per_hour = 3600*steps_per_periodim(atmos_im)/secs_per_periodim(atmos_im)

! Check to see if this is an instantaneous file
IF (mean_period == 0) THEN
  ! This section was always skipped for dump files in the previous code
  IF (dump) THEN
    ! Use the "full" style with hours
    WRITE(extension, "(i4.4, i2.2, i2.2, A1, i2.2)") &
        i_year, i_month, i_day, "_", i_hour
  ELSE
    ! Check for Gregorian reinit
    IF (reinit_steps < 0) THEN
      IF (i_day == 1) THEN
        ! If first day of month use month naming style
        WRITE(extension, "(i4.4, A3)") i_year, month_names_3char(i_month)
      ELSE
        ! Otherwise use the "full" style without hours
        WRITE(extension, "(i4.4, i2.2, i2.2)") i_year, i_month, i_day
      END IF
    ELSE
      ! If not using Gregorian re-init, work out reinit frequency in hours
      hours = reinit_steps/steps_per_hour
      ! See if it's a multiple of whole days
      IF (MOD(hours, 24) == 0) THEN
        ! ...and if it's the first day in what looks like monthly intervals
        IF (i_day == 1 .AND. hours/24 == 30) THEN
          ! Use the year + month name style
          WRITE(extension, "(i4.4, A3)") i_year, month_names_3char(i_month)
        ELSE
          ! Use the "full" style without hours
          WRITE(extension, "(i4.4, i2.2, i2.2)") i_year, i_month, i_day
        END IF
      ELSE
        ! Use the "full" style with hours
        WRITE(extension, "(i4.4, i2.2, i2.2, A1, i2.2)") &
            i_year, i_month, i_day, "_", i_hour
      END IF
    END IF
  END IF
ELSE
  ! If these are mean files, untangle the cascading mean syntax to work out
  ! how many hours there are in this meaning period (each period is a multiple
  ! of the last, with the first being the regular dumping period)
  hours = dumpfreqim(atmos_im)/steps_per_hour
  DO i = 1, mean_period
    hours = hours*meanfreqim(i,atmos_im)
  END DO
  ! Check if it's a special multiple of years
  special_char = ""
  IF (MOD(hours, 8640) == 0) THEN
    IF (hours == 8640) THEN
      special_char = "y"             !    1 Year
    ELSE IF (hours == 43200) THEN
      special_char = "v"             !    5 Years
    ELSE IF (hours == 86400) THEN
      special_char = "x"             !   10 Years
    ELSE IF (hours == 432000) THEN
      special_char = "l"             !   50 Years
    ELSE IF (hours == 864000) THEN
      special_char = "u"             !  100 Years
    ELSE IF (hours == 8640000) THEN
      special_char = "z"             ! 1000 Years
    END IF
  ELSE IF (MOD(hours, 24) == 0) THEN
    ! Or a special multiple of days
    IF (hours == 120) THEN
      special_char = "p"             !  5 Days
    ELSE IF (hours == 168) THEN
      special_char = "w"             !  1 Week
    ELSE IF (hours == 240) THEN
      special_char = "t"             ! 10 Days
    ELSE IF (hours == 336) THEN
      special_char = "r"             !  1 Fortnight

    ELSE IF (hours == 720 .AND. i_day == 1) THEN
      ! For monthly meaning, use the month array but note that since the mean
      ! files are produced at the end of a timestep they should be labelled with
      ! the month name from the previous month (and this is only done if the
      ! request happens to have been made on the 1st day of the month)
      month = i_month - 1
      IF (month == 0) THEN
        month = 12
        year = i_year - 1
      ELSE
        year = i_year
      END IF
      ! Use the year + month name style
      WRITE(extension, "(A1, i4.4, A3)") "m", year, month_names_3char(month)

    ELSE IF (hours == 2160 .AND. i_day == 1) THEN
      ! Similarly for seasonal (3-monthly) meaning the month index must be
      ! offset back to correspond to the beginning of the 3-month period.
      ! The year in the name will correspond to the *middle* month of this
      ! 3-month period in each case
      month = i_month - 3
      year  = i_year
      IF (month <= 0) THEN
        month = 12 + (i_month - 3)
        ! Roll back the year unless this is "djf" (the only season where the 
        ! middle month is not in the same year as the start month)
        IF (month /= 12) THEN
          year = i_year - 1
        END IF
      END IF

      ! Use the year + season name style
      WRITE(extension, "(A1, i4.4, A3)") "s", year, season_names_3char(month)

    END IF

  END IF

  ! If one of the special character branches above was set, the output name
  ! should be the "full" format, not including hours, and including that
  ! character at the front
  IF (special_char /= "") THEN
    WRITE(extension, "(A1, i4.4, i2.2, i2.2)") &
        special_char, i_year, i_month, i_day
  END IF

  ! Finally, if the extension hasn't been created by this point it defaults
  ! to a numerical name describing the meaning period and including the hours
  IF (TRIM(extension) == "") THEN
    WRITE(extension, "(i1.1, i4.4, i2.2, i2.2, A1, i2.2)") &
        mean_period, i_year, i_month, i_day, "_", i_hour
  END IF
END IF

END FUNCTION legacy_absolute_naming

END MODULE filename_generation_mod
