! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Coordinates the input of UKCA_RADAER look-up tables for each
!  spectra and aerosol modes.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_read_luts_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_READ_LUTS_MOD'

CONTAINS

SUBROUTINE ukca_radaer_read_luts(icode, cmessage)

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE glomap_clim_option_mod,  ONLY: &
    l_glomap_clim_radaer,          &
    gclmacsw_file => gclmacsw,     &
    gclmaclw_file => gclmaclw,     &
    gclmcrsw_file => gclmcrsw,     &
    gclmcrlw_file => gclmcrlw,     &
    gclmansw_file => gclmansw,     &
    gclmanlw_file => gclmanlw

USE parkind1,                ONLY: &
    jpim,                          &
    jprb

USE ukca_option_mod,         ONLY: &
    l_ukca_radaer,                 &
    ukcaacsw_file => ukcaacsw,     &
    ukcaaclw_file => ukcaaclw,     &
    ukcacrsw_file => ukcacrsw,     &
    ukcacrlw_file => ukcacrlw,     &
    ukcaansw_file => ukcaansw,     &
    ukcaanlw_file => ukcaanlw

USE ukca_radaer_lut,         ONLY: &
    ip_ukca_lut_accum,             &
    ip_ukca_lut_coarse,            &
    ip_ukca_lut_accnarrow,         &
    ip_ukca_lut_sw,                &
    ip_ukca_lut_lw,                &
    ukca_lut

USE ukca_radaer_lut_read_in, ONLY: &
    ukca_radaer_lut_in

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

!
! Arguments
!

! error indicator (0 is OK, <0 is KO)
INTEGER :: icode

! error message
CHARACTER (LEN=errormessagelength) :: cmessage

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_READ_LUTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

icode = 0

IF (l_glomap_clim_radaer) THEN
  IF (l_ukca_radaer) THEN
    icode = -2
    cmessage = 'l_glomap_clim_radaer and l_ukca_radaer must not both be true'
    CALL ereport(RoutineName, icode, cmessage)
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmacsw_file,                      &
                          ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
     IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,          &
                                                           zhook_handle)
     RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmcrsw_file, &
                          ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmansw_file,                      &
                          ukca_lut(ip_ukca_lut_accnarrow, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmaclw_file,                      &
                          ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmcrlw_file,                      &
                          ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, gclmanlw_file,                      &
                          ukca_lut(ip_ukca_lut_accnarrow, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

ELSE IF(l_ukca_radaer) THEN

  CALL ukca_radaer_lut_in(icode, cmessage, ukcaacsw_file,                      &
                          ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
     IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,          &
                                                           zhook_handle)
     RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, ukcacrsw_file, &
                          ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, ukcaansw_file,                      &
                          ukca_lut(ip_ukca_lut_accnarrow, ip_ukca_lut_sw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, ukcaaclw_file,                      &
                          ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, ukcacrlw_file,                      &
                          ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

  CALL ukca_radaer_lut_in(icode, cmessage, ukcaanlw_file,                      &
                          ukca_lut(ip_ukca_lut_accnarrow, ip_ukca_lut_lw))
  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF

ELSE
  icode    = 1
  cmessage = 'l_glomap_clim_radaer and l_ukca_radaer must not both be false'
  CALL ereport(RoutineName, icode, cmessage)
END IF 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_read_luts
END MODULE ukca_radaer_read_luts_mod
