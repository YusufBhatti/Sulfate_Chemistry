! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   module lionplyr_mod ---------------------------------------
!
!   Purpose: Calculates pressure layer fields from model layer fields
!   Used in the calculation of icing potential diagnostics.
!   Level Interpolation ONto Pressure LaYeR. lionplyr.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE lionplyr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LIONPLYR_MOD'

CONTAINS

SUBROUTINE lionplyr(NumLyrs,           &  ! in
                     req_plevels,       &  ! in
                     thickness,         &  ! in
                     Field_on_mls,      &  ! in
                     PFields ,          &  ! in
                     Field_on_press)       ! inout

! Description:
!   Calculates a vector of pressure level fields from a vector of model
!   level fields.
!
! Method:
!   Two vectors of pressures are supplied.  These contain the lower and
!   upper boundaries of a set of layers of the atmosphere.  For each of
!   these layers the corresponding data on the model level fields is
!   found using pressures on model levels supplied.  To each gridpoint
!   in the fields comprising the output, the mean of the corresponding
!   data is assigned.

USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels
USE missing_data_mod,      ONLY: imdi, rmdi
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, pdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Number of output (std level) levels
INTEGER, INTENT(IN) :: NumLyrs
! requested output std_levels
REAL, INTENT(IN)    :: req_plevels(NumLyrs)
! thickness of std levels to work with
REAL, INTENT(IN)    :: thickness
! Model level fields of parameter to be determined at std levels (no_halo)
REAL, INTENT(IN)    :: Field_on_mls(pdims%i_start:pdims%i_end,           &
                           pdims%j_start:pdims%j_end,                    &
                           pdims%k_start:pdims%k_end)
! Model level Pressure
REAL, INTENT(IN)    :: PFields(tdims_s%i_start:tdims_s%i_end,            &
                           tdims_s%j_start:tdims_s%j_end,                &
                           tdims_s%k_start:tdims_s%k_end)
! Output fields on std levels
REAL, INTENT(INOUT) :: Field_on_press(tdims%i_start:tdims%i_end,         &
                           tdims%j_start:tdims%j_end,                    &
                           NumLyrs)

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "LIONPLYR"

REAL    :: RData(model_levels)  ! Holds data part of LevFields(k)
REAL    :: Press(model_levels)  ! Holds pressure on model levels in hPa
LOGICAL :: NotMDI(model_levels) ! Mask showing where LevFields and corresponding
                                !   pressure fields do not have missing data
LOGICAL :: Mask(model_levels)   ! Holds mask showing which model level pressure
                                !   fields are within the slab of atmosphere

INTEGER :: i,j,k,n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Convert model level fields to pressure levels by determining the
! number of model pressure levels which are in the surrounding slab
! of atmosphere and calculating the mean of the model level fields
! in that slab of atmosphere

! First loop over gridpoints
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end

    DO k = 1, model_levels
      RData(k) = Field_on_mls(i,j,k)
      Press(k) = PFields(i,j,k)
    END DO

    NotMDI = (RData /= rmdi) .AND. (Press /= rmdi)

    ! Now loop over desired atmosphere layers
    DO k = 1, NumLyrs

      Mask = NotMDI .AND.                                                &
             (Press <= (req_plevels(k)+thickness*0.5)) .AND.             &
             (Press >= (req_plevels(k)-thickness*0.5))
      n = COUNT( Mask )

      IF (n > 0) THEN
        Field_on_press(i,j,k) = SUM( RData, Mask = Mask )
        Field_on_press(i,j,k) =  Field_on_press(i,j,k) / REAL( n )
      ELSE
        Field_on_press(i,j,k) = rmdi
      END IF

    END DO

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE LIOnPLyr

END MODULE lionplyr_mod
