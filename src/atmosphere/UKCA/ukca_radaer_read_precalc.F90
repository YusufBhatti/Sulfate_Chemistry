! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Reads in the pre-computed variables for use in the
!  interaction between UKCA aerosols and the radiation
!  code.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_read_precalc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'UKCA_RADAER_READ_PRECALC_MOD'

CONTAINS

SUBROUTINE ukca_radaer_read_precalc(                                    &
                                    ierr,                               &
                                    cmessage                            &
)

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE ereport_mod,     ONLY: ereport
USE umprintmgr, ONLY: newline
USE ukca_radaer_precalc, ONLY: &
    npd_integ_pts,             &
    npd_ukca_band,             &
    npd_ukca_spectrum,         &
    npd_ukca_maxcomptype,      &
    npd_ukca_aod_wavel,        &
    precalc

USE filenamelength_mod, ONLY: filenamelength

USE ukca_option_mod,        ONLY: &
    l_ukca_radaer,                &
    ukcaprec_file => ukcaprec

USE glomap_clim_option_mod, ONLY: &
    l_glomap_clim_radaer,         &
    gclmprec_file => gclmprec

USE file_manager, ONLY: assign_file_unit, release_file_unit
USE errormessagelength_mod, ONLY: errormessagelength

USE um_parcore, ONLY: mype
USE setup_namelist, ONLY: setup_nml_type

USE ukca_get_specinfo_mod, ONLY: ukca_get_specinfo
IMPLICIT NONE

!
! Arguments
!
!
! Error indicator
!
INTEGER :: ierr
!
! Error message
!
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=errormessagelength) :: iomessage
CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'UKCA_RADAER_READ_PRECALC'
!
! Local variables
!
CHARACTER (LEN=filenamelength) :: filename
INTEGER :: ios


CHARACTER (LEN=9) :: spectrum_name (npd_ukca_spectrum)
!
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode
!

INTEGER :: n_band

INTEGER :: aerosol_index

REAL :: wl_short
REAL :: wl_long
!
! Number of bands and AOD wavelengths in the spectral files.
! These must match the values read from the precalc file.
!
INTEGER :: specf_n_band(npd_ukca_spectrum)
INTEGER :: specf_n_aod_wavel
REAL :: wavelength_short(npd_ukca_spectrum, npd_ukca_band)
REAL :: wavelength_long(npd_ukca_spectrum, npd_ukca_band)

INTEGER :: funit

CHARACTER (LEN=80) :: line ! buffer for neglected lines (headers,...)

!
! Loop indices
!
INTEGER :: i, &
        j, &
        k, &
        n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = 2 * npd_integ_pts * npd_ukca_band *         &
                                                  npd_ukca_spectrum +      &
                      npd_integ_pts + npd_ukca_band * npd_ukca_spectrum +  &
 2*npd_ukca_maxcomptype*npd_integ_pts* npd_ukca_band* npd_ukca_spectrum +  &
                      npd_ukca_aod_wavel +                                 &
                      2* npd_ukca_maxcomptype* npd_ukca_aod_wavel

TYPE my_namelist
  SEQUENCE
  INTEGER :: n_integ_pts
  INTEGER :: n_ukca_aod_wavel
  REAL :: wavelength(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)
  REAL :: irrad(npd_integ_pts, npd_ukca_band, npd_ukca_spectrum)
  REAL :: weight(npd_integ_pts)
  REAL :: flux(npd_ukca_band, npd_ukca_spectrum)
  REAL :: realrefr(npd_ukca_maxcomptype, npd_integ_pts, npd_ukca_band,  &
                   npd_ukca_spectrum)
  REAL :: imagrefr(npd_ukca_maxcomptype, npd_integ_pts, npd_ukca_band,  &
                   npd_ukca_spectrum)
  REAL :: aod_wavel(npd_ukca_aod_wavel)
  REAL :: aod_realrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)
  REAL :: aod_imagrefr(npd_ukca_maxcomptype, npd_ukca_aod_wavel)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_int_in=n_int, n_real_in=n_real)
ierr = 0

!
! Get the number of wavebands as defined in the spectral files
!
CALL ukca_get_specinfo(  &
  specf_n_band(1),       &
  specf_n_band(2),       &
  specf_n_aod_wavel,     &
  wavelength_short,      &
  wavelength_long        &
  )


IF ( mype == 0 ) THEN
  !
  ! Get the name of the file to open and open it.
  !
  IF (l_glomap_clim_radaer) THEN
    filename = gclmprec_file
  ELSE IF(l_ukca_radaer) THEN
    filename = ukcaprec_file
  ELSE
    ierr = 1
    cmessage = 'l_glomap_clim_radaer and l_ukca_radaer must not both be false'
    CALL ereport(routinename, ierr, cmessage)
  END IF 
  IF (LEN_TRIM(filename) == 0) THEN
    ierr = -1
    cmessage = 'UKCA pre-computed filename is blank'
    CALL ereport(routinename, ierr, cmessage)
  END IF

  CALL assign_file_unit(filename, funit, handler="fortran")
  OPEN(UNIT=funit, FILE=filename, ACTION='READ', IOSTAT=ios, IOMSG=iomessage)

  IF (ios /= 0) THEN
    ierr = 1
    cmessage =                                                       newline//&
      'Error opening UKCA pre-computed file:'//                      newline//&
      TRIM(filename)//                                               newline//&
      'IoMsg: '//TRIM(iomessage)
    CALL ereport(routinename,ierr,cmessage)
  END IF

  !
  ! Read the formatted file.
  !
  !
  ! Number of integration points and weight of each point in the
  ! integration.
  !
  READ(funit, '(29x,i3)') precalc%n_integ_pts
  IF (precalc%n_integ_pts > npd_integ_pts) THEN
    ierr = 1
    cmessage = 'Increase the maximum number of integration points in ' // &
         'UKCA_RADAER pcalc file ' // filename
    CLOSE(funit)
    CALL release_file_unit(funit, handler="fortran")
    CALL ereport(routinename,ierr,cmessage)
  END IF

  READ(funit, '(a)') line ! column header

  DO n = 1, precalc%n_integ_pts

    READ(funit, '(4x,E12.5)') precalc%weight(n)

  END DO ! n

  !
  ! For each spectrum
  !
  DO i = 1, npd_ukca_spectrum

    READ(funit, '(a)') line ! spectrum header

    !
    ! Number of wavebands
    !

    spectrum_name = (/ 'shortwave', 'longwave ' /)

    READ(funit, '(20x,i3)') n_band
    IF (n_band /= specf_n_band(i)) THEN
      ierr = 1
      cmessage = 'UKCA RADAER pcalc file ' // filename // &
           ' is not compatible with the ' // &
                 spectrum_name(i) // ' spectral file.'
      CLOSE(funit)
      CALL release_file_unit(funit, handler="fortran")
      CALL ereport(routinename,ierr,cmessage)
    END IF

    IF (n_band > npd_ukca_band) THEN
      ierr = 1
      cmessage = &
        'Too many wavebands. Increase the maximum number of wavebands.'
      CLOSE(funit)
      CALL release_file_unit(funit, handler="fortran")
      CALL ereport(routinename,ierr,cmessage)
    END IF

    !
    ! Loop on all wavebands
    !
    DO j = 1, n_band

      READ(funit, '(14x,E12.5,3x,E12.5,2x)') wl_short, wl_long
      IF ( (ABS(wl_short-wavelength_short(i,j)) / wl_short  >=  1.0e-5) .OR. &
           (ABS(wl_long -wavelength_long(i,j) ) / wl_long   >=  1.0e-5)) THEN
        ierr = 1
        cmessage = &
          'Wavelength limits in RADAER PCALC.ukca file do not match those in ' // &
                   spectrum_name(i) // ' spectral file.'
        CLOSE(funit)
        CALL release_file_unit(funit, handler="fortran")
        CALL ereport(routinename,ierr,cmessage)
      END IF

      READ(funit, '(a)') line ! column header

      !
      ! Wavelength at integration point and monochromatic solar
      ! irradiance or Planckian at this wavelength
      !
      DO n = 1, precalc%n_integ_pts

        READ(funit, '(6x,2(E12.5,1x))') &
          precalc%wavelength(n, j, i), precalc%irrad(n, j, i)

      END DO ! n (integration points)

      !
      ! Waveband-integrated irradiance or Planckian
      !
      READ(funit, '(27x,E12.5)') precalc%flux(j, i)

      !
      ! Complex refractive index of aerosols at integration points.
      !
      READ(funit, '(a)') line ! header

      DO k = 1, npd_ukca_maxcomptype

        READ(funit, '(24x,i2)') aerosol_index

        READ(funit, '(a)') line ! column header

        DO n = 1, precalc%n_integ_pts

          READ(funit, '(8x,2(E12.5,1x))') &
            precalc%realrefr(aerosol_index, n, j, i), &
            precalc%imagrefr(aerosol_index, n, j, i)

        END DO ! n

      END DO ! k

    END DO ! j (wavebands)

  END DO ! i (spectrum)

  !
  ! Number and values of wavelengths for the aerosol optical depth
  !
  READ(funit, '(a)') line ! header

  READ(funit, '(22x,i3)') precalc%n_ukca_aod_wavel
  IF (precalc%n_ukca_aod_wavel > npd_ukca_aod_wavel) THEN
    ierr = 1
    cmessage = 'Increase npd_ukca_aod_wavel in UKCAPCALC.'
    CLOSE(funit)
    CALL release_file_unit(funit, handler="fortran")
    CALL ereport(routinename,ierr,cmessage)
  END IF

  READ(funit, '(a)') line ! column header

  DO i = 1, precalc%n_ukca_aod_wavel

    READ(funit, '(4x,E12.5)') precalc%aod_wavel(i)

  END DO ! i

  !
  ! Complex refractive index of aerosols at AOD wavelengths.
  !
  READ(funit, '(a)') line ! header

  DO k = 1, npd_ukca_maxcomptype

    READ(funit, '(24x,i2)') aerosol_index

    READ(funit, '(a)') line ! column header

    DO n = 1, precalc%n_ukca_aod_wavel

      READ(funit, '(8x,2(E12.5,1x))') &
        precalc%aod_realrefr(aerosol_index, n), &
        precalc%aod_imagrefr(aerosol_index, n)

    END DO ! n

  END DO ! k

  CLOSE(funit)
  CALL release_file_unit(funit, handler="fortran")

  my_nml % n_integ_pts      = precalc % n_integ_pts
  my_nml % n_ukca_aod_wavel = precalc % n_ukca_aod_wavel
  my_nml % wavelength       = precalc % wavelength
  my_nml % irrad            = precalc % irrad
  my_nml % weight           = precalc % weight
  my_nml % flux             = precalc % flux
  my_nml % realrefr         = precalc % realrefr
  my_nml % imagrefr         = precalc % imagrefr
  my_nml % aod_wavel        = precalc % aod_wavel
  my_nml % aod_realrefr     = precalc % aod_realrefr
  my_nml % aod_imagrefr     = precalc % aod_imagrefr

END IF !!mype=0

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF (mype /= 0) THEN

  precalc % n_integ_pts      = my_nml % n_integ_pts
  precalc % n_ukca_aod_wavel = my_nml % n_ukca_aod_wavel
  precalc % wavelength       = my_nml % wavelength
  precalc % irrad            = my_nml % irrad
  precalc % weight           = my_nml % weight
  precalc % flux             = my_nml % flux
  precalc % realrefr         = my_nml % realrefr
  precalc % imagrefr         = my_nml % imagrefr
  precalc % aod_wavel        = my_nml % aod_wavel
  precalc % aod_realrefr     = my_nml % aod_realrefr
  precalc % aod_imagrefr     = my_nml % aod_imagrefr

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_read_precalc
END MODULE ukca_radaer_read_precalc_mod
