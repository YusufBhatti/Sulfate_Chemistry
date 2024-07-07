! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine INITPHYS
!
! Purpose:
!   Calls the routines required to read in spectral files and other
!   data files for the radiation scheme.
!
! Documentation:
!   UM documentation paper no 23
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!- ---------------------------------------------------------------------
SUBROUTINE initphys(icode, cmessage)

USE rad_input_mod, ONLY: n_swcall, n_lwcall
USE sw_control_struct, ONLY: sw_control
USE lw_control_struct, ONLY: lw_control
USE spec_sw_lw, ONLY: sw_spectrum, lw_spectrum
USE mcica_mod, ONLY: read_mcica_data
USE rad_pcf, ONLY: i_err_io, ip_cloud_mcica
USE filenamelength_mod, ONLY: filenamelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_option_mod, ONLY: l_ukca_radaer
USE glomap_clim_option_mod, ONLY: l_glomap_clim_radaer
USE umPrintMgr, ONLY: umPrint, umMessage
USE um_parcore, ONLY: mype
USE setup_spectra_mod, ONLY: setup_spectra
USE errormessagelength_mod, ONLY: errormessagelength
USE get_env_var_mod, ONLY: get_env_var

USE ukca_radaer_read_luts_mod, ONLY: ukca_radaer_read_luts
USE ukca_radaer_read_precalc_mod, ONLY: ukca_radaer_read_precalc
USE compress_spectrum_mod, ONLY: compress_spectrum
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: icode             ! Return code
CHARACTER(LEN=errormessagelength), INTENT(INOUT) :: cmessage ! Error message

! Local variables
INTEGER :: j, path_end
CHARACTER (LEN=filenamelength) :: spectral_file_dir
!   Common directory for spectral files
CHARACTER (LEN=filenamelength) :: mcica_data
!   Path to McICA data file

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITPHYS'


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

CALL umPrint('', src='initphys')
CALL umPrint('*********************** Reading spectral files '//      &
  '***********************', src='initphys')

! Get the base directory for the spectral files from the
! environment variable $SPECTRAL_FILE_DIR.
CALL get_env_var("SPECTRAL_FILE_DIR", spectral_file_dir, length=path_end)

! Add a trailing / if not present.
IF (spectral_file_dir(path_end:path_end) /= '/') THEN
  path_end = path_end + 1
  spectral_file_dir(path_end:path_end) = '/'
END IF

! ------------- Shortwave Radiation -------------------------
IF (mype == 0) THEN
  DO j=1, n_swcall
    ! Relative paths are appended onto spectral_file_dir
    IF (sw_control(j)%spectral_file(1:1) /= '/') THEN
      sw_control(j)%spectral_file(path_end+1:) =                      &
        sw_control(j)%spectral_file(1:)
      sw_control(j)%spectral_file(:path_end) =                        &
        spectral_file_dir(:path_end)
    END IF
    CALL umPrint('SW spectral file: '//                               &
      TRIM(sw_control(j)%spectral_file), src='initphys')
    ! DEPENDS ON: read_spectrum
    CALL read_spectrum( sw_control(j)%spectral_file, sw_spectrum(j) )
    CALL compress_spectrum( sw_control(j), sw_spectrum(j) )
  END DO
END IF  !! mype==0
IF (n_swcall > 0) THEN
! Broadcast spectral file to all other PEs
  CALL setup_spectra(sw_control, sw_spectrum, n_swcall)
END IF


! ------------- Longwave Radiation --------------------------
IF (mype == 0) THEN
  DO j=1,n_lwcall
    ! Relative paths are appended onto spectral_file_dir
    IF (lw_control(j)%spectral_file(1:1) /= '/') THEN
      lw_control(j)%spectral_file(path_end+1:) =                      &
        lw_control(j)%spectral_file(1:)
      lw_control(j)%spectral_file(:path_end) =                        &
        spectral_file_dir(:path_end)
    END IF
    CALL umPrint('LW spectral file: '//                               &
      TRIM(lw_control(j)%spectral_file), src='initphys')
    CALL read_spectrum( lw_control(j)%spectral_file, lw_spectrum(j) )
    CALL compress_spectrum( lw_control(j), lw_spectrum(j) )
  END DO
END IF  !! mype==0
IF (n_lwcall > 0) THEN
! Broadcast spectral file to all other PEs
  CALL setup_spectra(lw_control, lw_spectrum, n_lwcall)
END IF


! ------------- MCICA data file -----------------------------
IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                    &
    (lw_control(1)%i_cloud == ip_cloud_mcica)) THEN

  IF (mype == 0) THEN
    path_end=SCAN(sw_control(1)%spectral_file,'/',.TRUE.)
    mcica_data(:path_end)=sw_control(1)%spectral_file(:path_end)
    mcica_data(path_end+1:)='mcica_data'

    CALL umPrint('',src='initphys')
    WRITE(umMessage,'(A,A)') 'Reading McICA data file: ',mcica_data
    CALL umPrint(umMessage,src='initphys')
  END IF  !! mype==0
  CALL read_mcica_data(mcica_data)
END IF


! ------------- UKCA_RADAER ---------------------------------
IF ( l_ukca_radaer .OR. l_glomap_clim_radaer ) THEN
  CALL ukca_radaer_read_luts(icode, cmessage)
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
  CALL ukca_radaer_read_precalc(icode, cmessage)
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF


CALL umPrint( '**********************************************'//      &
  '************************', src='initphys')
CALL umPrint('', src='initphys')
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE initphys
