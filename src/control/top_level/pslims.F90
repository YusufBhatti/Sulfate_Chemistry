! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns the limits of the pseudo levels
! Subroutine Interface:
SUBROUTINE pslims(ipfirst,iplast,ifirst,ilast)

USE jules_snow_mod, ONLY: nsmax
USE clmchfcg_scenario_mod, ONLY: nsulpat

USE jules_surface_types_mod
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE um_ParParams
USE cosp_constants_mod, ONLY: sr_bins,dbze_bins,n_hydro
USE cosp_input_mod, ONLY: cosp_ncolumns_max
USE aodband_mod, ONLY: aod_bands
USE rad_input_mod, ONLY: h_swbands,h_lwbands
USE jules_sea_seaice_mod, ONLY: nice
USE nlsizes_namelist_mod, ONLY: ntiles
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Subroutine arguments:
!   Input: First and last pseudo-level codes from STASHmaster
INTEGER, INTENT(IN) :: ipfirst
INTEGER, INTENT(IN) :: iplast

!   Output: First and last level numbers that the codes refer to
INTEGER, INTENT(OUT) :: ifirst
INTEGER, INTENT(OUT) :: ilast

! Local variables
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER  :: err_code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PSLIMS'

!- End of Header --------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (ipfirst == 1) THEN
  ifirst=1
ELSE
  err_code = 2
  WRITE(cmessage,'(A,I5)') 'pslims: Invalid first pseudo-level code ',ipfirst
  CALL ereport('pslevcod',err_code,cmessage)
END IF

IF (iplast == 1) THEN
  ilast=h_swbands
ELSE IF (iplast == 2) THEN
  ilast=h_lwbands
ELSE IF (iplast == 3) THEN
  !Aerosol Optical Depth diagnostics
  ilast = aod_bands
ELSE IF (iplast == 6) THEN
  !Sulphate loading patterns
  ilast=nsulpat
ELSE IF (iplast == 7) THEN
  !Direct vegetation parametrization: all surface types
  ilast=ntype
ELSE IF (iplast == 8) THEN
  !Direct vegetation parametrization: only plant functional types
  ilast=npft
ELSE IF (iplast == 9) THEN
  !Direct vegetation parametrization: all tiles
  ilast=ntiles
ELSE IF (iplast == 10) THEN
  !All sea ice categories
  ilast=nice
ELSE IF (iplast == 11) THEN
  !New snow scheme: all tiles <times> max no. of snow levels
  ilast=ntiles*nsmax
ELSE IF (iplast == 12) THEN
  !COSP: Radar reflectivity bins
  ilast=dbze_bins
ELSE IF (iplast == 13) THEN
  !COSP: Number of hydrometeors
  ilast=n_hydro
ELSE IF (iplast == 14) THEN
  !COSP: Backscattering ratio bins
  ilast=sr_bins
ELSE IF (iplast == 15) THEN
  !COSP: ISCCP tau bins
  ilast=7
ELSE IF (iplast == 16) THEN
  !COSP: Max number of subcolumns allowed to be output
  ilast=cosp_ncolumns_max
ELSE
  err_code = 2
  WRITE(cmessage,'(A,I5)') 'pslims: Invalid last pseudo-level code ',iplast
  CALL ereport('pslevcod',err_code,cmessage)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pslims

