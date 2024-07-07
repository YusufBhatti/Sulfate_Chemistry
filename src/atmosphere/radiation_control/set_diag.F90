! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the diagnostic fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_DIAG_MOD'
CONTAINS

SUBROUTINE set_diag(n_points, n_layer, list, col_list, row_list, model_levels, &
  nd_field, obs_solid_angle, trans_solid_angle, dir_flux_to_trans, &
  control, atm, spectrum, radout, radout_clean, radout_forc, diag)

USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm
USE def_spectrum, ONLY: StrSpecData
USE def_diag,     ONLY: StrDiag
USE def_out,      ONLY: StrOut, deallocate_out

USE rad_pcf, ONLY: ip_two_stream, ip_spherical_harmonic, ip_sph_mode_flux
USE conversions_mod, ONLY: pi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER, INTENT(IN) :: n_points
!   Number of atmospheric points
INTEGER, INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations
INTEGER, INTENT(IN) :: list(n_points)
!   List of points where radiation is to be calculated
INTEGER, INTENT(IN) :: col_list(n_points)
!   List of column indices
INTEGER, INTENT(IN) :: row_list(n_points)
!   List of row indices
INTEGER, INTENT(IN) :: model_levels
!   Number of layers in full model
INTEGER, INTENT(IN) :: nd_field
!   Field size in full model
REAL, INTENT(IN) :: obs_solid_angle(nd_field)
!   Solid angle subtended by grid-box for an observer at 1 AU
REAL, INTENT(IN) :: trans_solid_angle(nd_field)
!   Solid angle subtended by grid-box for a transit observer at 1 AU
REAL, INTENT(IN) :: dir_flux_to_trans(nd_field)
!   Conversion factor from direct flux to flux at transit observer

! Controlling options:
TYPE (StrCtrl), INTENT(IN) :: control

! Atmospheric properties:
TYPE(StrAtm), INTENT(IN) :: atm

! Spectral data:
TYPE (StrSpecData), INTENT(IN) :: spectrum

! Output fields from core radiation code:
TYPE(StrOut), INTENT(INOUT) :: radout

! Output fields from diagnostic calls:
TYPE(StrOut), INTENT(INOUT) :: radout_clean
TYPE(StrOut), INTENT(INOUT) :: radout_forc

! Diagnostics:
TYPE (StrDiag), INTENT(INOUT) :: diag


! Local variables.
INTEGER :: i, l, k
!   Loop variables

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_DIAG'
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER :: ierr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Clean-air diagnostics
IF (diag%l_flux_up_clean) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_up_clean(col_list(l), row_list(l),i) &
        = radout_clean%flux_up(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_down_clean) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_down_clean(col_list(l), row_list(l),i) &
        = radout_clean%flux_down(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_up_clear_clean) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_up_clear_clean(col_list(l), row_list(l),i) &
        = radout_clean%flux_up_clear(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_down_clear_clean) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_down_clear_clean(col_list(l), row_list(l),i) &
        = radout_clean%flux_down_clear(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_up_clean_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up_clean_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_up_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_down_clean_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down_clean_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_down_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_up_clear_clean_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up_clear_clean_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_up_clear_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_down_clear_clean_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down_clear_clean_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_down_clear_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_emission_spectrum_clean) THEN
  DO k=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%emission_spectrum_clean(col_list(l), row_list(l), k) &
        = radout_clean%flux_up_band(l, 0, k) * obs_solid_angle(list(l)) / pi
    END DO
  END DO
END IF
IF (diag%l_emission_spectrum_clear_clean) THEN
  DO k=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%emission_spectrum_clear_clean(col_list(l), row_list(l), k) &
        = radout_clean%flux_up_clear_band(l, 0, k)*obs_solid_angle(list(l))/pi
    END DO
  END DO
END IF

! Clean-air direct fluxes for spherical geometry
IF (diag%l_flux_direct_clean_sph) THEN
  DO i=0, model_levels
    DO l=1, n_points
      diag%flux_direct_clean_sph(col_list(l), row_list(l), i) &
        = radout_clean%flux_direct_sph(l, n_layer+1-i, 1)
    END DO
  END DO
  i=model_levels+1
  DO l=1, n_points
    diag%flux_direct_clean_sph(col_list(l), row_list(l), i) &
      = radout_clean%flux_direct_sph(l, 0, 1)
  END DO
END IF
IF (diag%l_flux_direct_clean_div) THEN
  DO i=1, model_levels
    DO l=1, n_points
      diag%flux_direct_clean_div(col_list(l), row_list(l),i) &
        = radout_clean%flux_direct_div(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_clean_sph) THEN
  DO i=0, model_levels
    DO l=1, n_points
      diag%flux_direct_clear_clean_sph(col_list(l), row_list(l), i) &
        = radout_clean%flux_direct_clear_sph(l, n_layer+1-i, 1)
    END DO
  END DO
  i=model_levels+1
  DO l=1, n_points
    diag%flux_direct_clear_clean_sph(col_list(l), row_list(l), i) &
      = radout_clean%flux_direct_clear_sph(l, 0, 1)
  END DO
END IF
IF (diag%l_flux_direct_clear_clean_div) THEN
  DO i=1, model_levels
    DO l=1, n_points
      diag%flux_direct_clear_clean_div(col_list(l), row_list(l),i) &
        = radout_clean%flux_direct_clear_div(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clean_sph_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_clean_sph_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_direct_sph_band(l, n_layer+1-i, k)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_clean_sph_band(col_list(l), row_list(l), i, k) &
        = radout_clean%flux_direct_sph_band(l, 0, k)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clean_div_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_clean_div_band(col_list(l), row_list(l), i, k) &
          = radout_clean%flux_direct_div_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_clean_sph_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_clean_sph_band(col_list(l), row_list(l),i,k) &
          = radout_clean%flux_direct_clear_sph_band(l, n_layer+1-i, k)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_clear_clean_sph_band(col_list(l), row_list(l),i,k) &
        = radout_clean%flux_direct_clear_sph_band(l, 0, k)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_clean_div_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_clean_div_band(col_list(l),row_list(l),i,k) &
          = radout_clean%flux_direct_clear_div_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_transmission_spectrum_clean) THEN
  DO k=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%transmission_spectrum_clean(col_list(l), row_list(l), k) &
        = - radout_clean%flux_direct_sph_band(l, 0, k) &
          * dir_flux_to_trans(list(l)) &
          + radout_clean%flux_up_band(l, 0, k) &
          * trans_solid_angle(list(l)) / pi
    END DO
  END DO
END IF
IF (diag%l_transmission_spectrum_clear_clean) THEN
  DO k=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%transmission_spectrum_clear_clean(col_list(l), row_list(l), k) &
        = - radout_clean%flux_direct_clear_sph_band(l, 0, k) &
          * dir_flux_to_trans(list(l)) &
          + radout_clean%flux_up_clear_band(l, 0, k) &
          * trans_solid_angle(list(l)) / pi  
    END DO
  END DO
END IF

CALL deallocate_out(radout_clean)


! GHG forcing diagnostics
IF (diag%l_flux_up_forc) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_up_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_up(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_down_forc) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_down_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_down(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_up_clear_forc) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_up_clear_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_up_clear(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_down_clear_forc) THEN
  DO i=1, model_levels+1
    DO l=1, n_points
      diag%flux_down_clear_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_down_clear(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_up_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_up_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_down_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_down_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_up_clear_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up_clear_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_up_clear_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_down_clear_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down_clear_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_down_clear_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_direct_sph_forc) THEN
  DO i=0, model_levels
    DO l=1, n_points
      diag%flux_direct_sph_forc(col_list(l), row_list(l), i) &
        = radout_forc%flux_direct_sph(l, n_layer+1-i, 1)
    END DO
  END DO
  i=model_levels+1
  DO l=1, n_points
    diag%flux_direct_sph_forc(col_list(l), row_list(l), i) &
      = radout_forc%flux_direct_sph(l, 0, 1)
  END DO
END IF
IF (diag%l_flux_direct_div_forc) THEN
  DO i=1, model_levels
    DO l=1, n_points
      diag%flux_direct_div_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_direct_div(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_sph_forc) THEN
  DO i=0, model_levels
    DO l=1, n_points
      diag%flux_direct_clear_sph_forc(col_list(l), row_list(l), i) &
        = radout_forc%flux_direct_clear_sph(l, n_layer+1-i, 1)
    END DO
  END DO
  i=model_levels+1
  DO l=1, n_points
    diag%flux_direct_clear_sph_forc(col_list(l), row_list(l), i) &
      = radout_forc%flux_direct_clear_sph(l, 0, 1)
  END DO
END IF
IF (diag%l_flux_direct_clear_div_forc) THEN
  DO i=1, model_levels
    DO l=1, n_points
      diag%flux_direct_clear_div_forc(col_list(l), row_list(l),i) &
        = radout_forc%flux_direct_clear_div(l,n_layer+1-i, 1)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_sph_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_sph_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_direct_sph_band(l, n_layer+1-i, k)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_sph_forc_band(col_list(l), row_list(l), i, k) &
        = radout_forc%flux_direct_sph_band(l, 0, k)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_div_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_div_forc_band(col_list(l), row_list(l), i, k) &
          = radout_forc%flux_direct_div_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_sph_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_sph_forc_band(col_list(l), row_list(l), i,k) &
          = radout_forc%flux_direct_clear_sph_band(l, n_layer+1-i, k)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_clear_sph_forc_band(col_list(l), row_list(l), i,k) &
        = radout_forc%flux_direct_clear_sph_band(l, 0, k)
    END DO
  END DO
END IF
IF (diag%l_flux_direct_clear_div_forc_band) THEN
  DO k=1, spectrum%basic%n_band
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_div_forc_band(col_list(l),row_list(l),i,k) &
          = radout_forc%flux_direct_clear_div_band(l, n_layer+1-i, k)
      END DO
    END DO
  END DO
END IF

CALL deallocate_out(radout_forc)


! Processing depends on whether the code has been invoked to
! calculate radiances or fluxes.
IF ( (control%i_angular_integration == ip_two_stream) .OR. &
     ( (control%i_angular_integration == ip_spherical_harmonic) .AND. &
        (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

  ! Mask for radiation calculations
  IF (diag%l_rad_mask) THEN
    DO l=1, n_points
      diag%rad_mask(col_list(l), row_list(l)) = 1.0
    END DO
  END IF

  ! Total upward & downward all-sky & clear-sky flux at all levels:
  IF (diag%l_flux_up) THEN
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up(col_list(l), row_list(l),i) &
          = radout%flux_up(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_down) THEN
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down(col_list(l), row_list(l),i) &
          = radout%flux_down(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_up_clear) THEN
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_up_clear(col_list(l), row_list(l),i) &
          = radout%flux_up_clear(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_down_clear) THEN
    DO i=1, model_levels+1
      DO l=1, n_points
        diag%flux_down_clear(col_list(l), row_list(l),i) &
          = radout%flux_down_clear(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_up_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels+1
        DO l=1, n_points
          diag%flux_up_band(col_list(l), row_list(l), i, k) &
            = radout%flux_up_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_flux_down_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels+1
        DO l=1, n_points
          diag%flux_down_band(col_list(l), row_list(l), i, k) &
            = radout%flux_down_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_flux_up_clear_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels+1
        DO l=1, n_points
          diag%flux_up_clear_band(col_list(l), row_list(l), i, k) &
            = radout%flux_up_clear_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_flux_down_clear_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels+1
        DO l=1, n_points
          diag%flux_down_clear_band(col_list(l), row_list(l), i, k) &
            = radout%flux_down_clear_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_emission_spectrum) THEN
    DO k=1, spectrum%basic%n_band
      DO l=1, n_points
        diag%emission_spectrum(col_list(l), row_list(l), k) &
          = radout%flux_up_band(l, 0, k) * obs_solid_angle(list(l)) / pi
      END DO
    END DO
  END IF
  IF (diag%l_emission_spectrum_clear) THEN
    DO k=1, spectrum%basic%n_band
      DO l=1, n_points
        diag%emission_spectrum_clear(col_list(l), row_list(l), k) &
          = radout%flux_up_clear_band(l, 0, k) * obs_solid_angle(list(l)) / pi
      END DO
    END DO
  END IF

  ! Direct fluxes for spherical geometry
  IF (diag%l_flux_direct_sph) THEN
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_sph(col_list(l), row_list(l), i) &
          = radout%flux_direct_sph(l, n_layer+1-i, 1)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_sph(col_list(l), row_list(l), i) &
        = radout%flux_direct_sph(l, 0, 1)
    END DO
  END IF
  IF (diag%l_flux_direct_div) THEN
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_div(col_list(l), row_list(l),i) &
          = radout%flux_direct_div(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_direct_clear_sph) THEN
    DO i=0, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_sph(col_list(l), row_list(l), i) &
          = radout%flux_direct_clear_sph(l, n_layer+1-i, 1)
      END DO
    END DO
    i=model_levels+1
    DO l=1, n_points
      diag%flux_direct_clear_sph(col_list(l), row_list(l), i) &
        = radout%flux_direct_clear_sph(l, 0, 1)
    END DO
  END IF
  IF (diag%l_flux_direct_clear_div) THEN
    DO i=1, model_levels
      DO l=1, n_points
        diag%flux_direct_clear_div(col_list(l), row_list(l),i) &
          = radout%flux_direct_clear_div(l,n_layer+1-i, 1)
      END DO
    END DO
  END IF
  IF (diag%l_flux_direct_sph_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=0, model_levels
        DO l=1, n_points
          diag%flux_direct_sph_band(col_list(l), row_list(l), i, k) &
            = radout%flux_direct_sph_band(l, n_layer+1-i, k)
        END DO
      END DO
      i=model_levels+1
      DO l=1, n_points
        diag%flux_direct_sph_band(col_list(l), row_list(l), i, k) &
          = radout%flux_direct_sph_band(l, 0, k)
      END DO
    END DO
  END IF
  IF (diag%l_flux_direct_div_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels
        DO l=1, n_points
          diag%flux_direct_div_band(col_list(l), row_list(l), i, k) &
            = radout%flux_direct_div_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_flux_direct_clear_sph_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=0, model_levels
        DO l=1, n_points
          diag%flux_direct_clear_sph_band(col_list(l), row_list(l), i, k) &
            = radout%flux_direct_clear_sph_band(l, n_layer+1-i, k)
        END DO
      END DO
      i=model_levels+1
      DO l=1, n_points
        diag%flux_direct_clear_sph_band(col_list(l), row_list(l), i, k) &
          = radout%flux_direct_clear_sph_band(l, 0, k)
      END DO
    END DO
  END IF
  IF (diag%l_flux_direct_clear_div_band) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels
        DO l=1, n_points
          diag%flux_direct_clear_div_band(col_list(l), row_list(l), i, k) &
            = radout%flux_direct_clear_div_band(l, n_layer+1-i, k)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_transmission_spectrum) THEN
    DO k=1, spectrum%basic%n_band
      DO l=1, n_points
        diag%transmission_spectrum(col_list(l), row_list(l), k) &
          = - radout%flux_direct_sph_band(l, 0, k) &
            * dir_flux_to_trans(list(l)) &
            + radout%flux_up_band(l, 0, k) &
            * trans_solid_angle(list(l)) / pi 
      END DO
    END DO
  END IF
  IF (diag%l_transmission_spectrum_clear) THEN
    DO k=1, spectrum%basic%n_band
      DO l=1, n_points
        diag%transmission_spectrum_clear(col_list(l), row_list(l), k) &
          = - radout%flux_direct_clear_sph_band(l, 0, k) &
            * dir_flux_to_trans(list(l)) &
            + radout%flux_up_clear_band(l, 0, k) &
            * trans_solid_angle(list(l)) / pi
      END DO
    END DO
  END IF


  ! Aerosol optical properties
  IF (diag%l_aerosol_optical_depth) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels
        DO l=1, n_points
          diag%aerosol_optical_depth(col_list(l), row_list(l), i, k) &
            = (radout%aerosol_absorption_band(l, n_layer+1-i, k) &
             + radout%aerosol_scattering_band(l, n_layer+1-i, k)) &
             * atm%mass(l, n_layer+1-i)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_aerosol_scat_optical_depth) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels
        DO l=1, n_points
          diag%aerosol_scat_optical_depth(col_list(l), row_list(l), i, k) &
            = radout%aerosol_scattering_band(l, n_layer+1-i, k) &
            * atm%mass(l, n_layer+1-i)
        END DO
      END DO
    END DO
  END IF
  IF (diag%l_aerosol_asymmetry_scat) THEN
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels
        DO l=1, n_points
          diag%aerosol_asymmetry_scat(col_list(l), row_list(l), i, k) &
            = radout%aerosol_asymmetry_band(l, n_layer+1-i, k) &
            * atm%mass(l, n_layer+1-i)
        END DO
      END DO
    END DO
  END IF

ELSE

  ! Radiances are being calculated. This code is written on the
  ! assumption that only one radiance is calculated at a point:
  ! this can of course be generalized.

  IF (diag%l_toa_radiance) THEN
    DO k=1, control%n_channel
      DO l=1, n_points
        diag%toa_radiance(col_list(l), row_list(l), k) &
          = radout%radiance(l, 1, 1, k)
      END DO
    END DO
  END IF

END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_diag
END MODULE set_diag_mod
