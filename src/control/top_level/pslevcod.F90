! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Decode the STASH pseudo level code
! Subroutine Interface:
SUBROUTINE pslevcod(ilin,ilout,swtch,ErrorStatus,cmessage)

USE jules_surface_types_mod
USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice
USE clmchfcg_scenario_mod, ONLY: nsulpat
USE rad_input_mod, ONLY: h_swbands,h_lwbands

USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE um_ParParams
USE cosp_constants_mod, ONLY: sr_bins,dbze_bins,n_hydro
USE cosp_input_mod, ONLY: cosp_ncolumns_max
USE aodband_mod, ONLY: aod_bands
USE nlsizes_namelist_mod, ONLY: &
    ntiles

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Description:
!   Sets ILOUT to an appropriate pseudo level size according
!    to the value of IL
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!

! Subroutine arguments:
INTEGER,          INTENT(IN) :: ilin  ! Model pseudo level code
CHARACTER(LEN=1) ,INTENT(IN) :: swtch ! Set to 'F' or 'L' to indicate if
                                      ! first or last level code to be checked

INTEGER,          INTENT(OUT):: ilout   ! An actual pseudo level
CHARACTER(LEN=errormessagelength),INTENT(OUT):: cmessage

! Local scalars:
INTEGER :: i
INTEGER :: j

! Error Status:
INTEGER :: ErrorStatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PSLEVCOD'

!- End of Header --------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (swtch == 'F') THEN
  IF (ilin == 1) THEN
    ilout=1
    ! Ocean assimilation groups - code removed
  ELSE
    WRITE(cmessage,'(A,I5)')                                        &
       'pslevcod: Inappropriate first pseudo level code found ',ilin
    ErrorStatus=2
    CALL ereport('pslevcod',ErrorStatus,cmessage)
  END IF
ELSE IF (swtch == 'L') THEN
  IF (ilin == 1) THEN
    ilout=h_swbands
  ELSE IF (ilin == 2) THEN
    ilout=h_lwbands
  ELSE IF (ilin == 3) THEN
    ! Radiation bands for measuring aerosol optical depth
    ilout=aod_bands
  ELSE IF ( ilin  ==  6 ) THEN
    ! Last index for HadCM2 sulphate loading patterns.
    ilout = nsulpat
  ELSE IF ( ilin  ==  7 ) THEN
    ! All surface types
    ilout = ntype
  ELSE IF ( ilin  ==  8 ) THEN
    ! Plant functional types only
    ilout = npft
  ELSE IF ( ilin  ==  9 ) THEN
    ! All tiles
    ilout = ntiles
  ELSE IF ( ilin  ==  10 ) THEN
    ! All sea ice categories
    ilout = nice
  ELSE IF ( ilin  ==  11 ) THEN
    !New snow scheme: all tiles <times> max no. of snow levels
    ilout = ntiles*nsmax
  ELSE IF ( ilin  ==  12 ) THEN
    !COSP: Radar reflectivity bins
    ilout = dbze_bins
  ELSE IF ( ilin  ==  13 ) THEN
    !COSP: Number of hydrometeors
    ilout = n_hydro
  ELSE IF ( ilin  ==  14 ) THEN
    !COSP: Backscattering ratio bins
    ilout = sr_bins
  ELSE IF ( ilin  ==  15 ) THEN
    !COSP: ISCCP tau bins
    ilout = 7
  ELSE IF ( ilin  ==  16 ) THEN
    !COSP: Max number of subcolumns allowed to be output
    ilout = cosp_ncolumns_max
  ELSE
    WRITE(cmessage,'(A,I5)')                                       &
       'pslevcod: Inappropriate last pseudo level code found ',ilin
    ErrorStatus=2
    CALL ereport('pslevcod',ErrorStatus,cmessage)
  END IF

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pslevcod

