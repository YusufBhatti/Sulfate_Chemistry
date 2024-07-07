! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module profile_interp_mod
MODULE prof_interp_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

! Interpolation methods
INTEGER, PARAMETER :: interp_constant=0
INTEGER, PARAMETER :: interp_linear=1

! Extrapolation methods
INTEGER, PARAMETER :: extrap_constant=0
INTEGER, PARAMETER :: extrap_linear=1
INTEGER, PARAMETER :: extrap_zero=2
INTEGER, PARAMETER :: extrap_linear_to_zero=3

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROF_INTERP_MOD'

CONTAINS

SUBROUTINE prof_interp(interp_method, extrap_method,                           &
                       ndata, zdata, fdata,                                    &
                       n_out, z_out, f_out)

! Purpose:
!          Sets up initial data for idealised problems.
!          Interpolate data provided in idealised namelist onto
!          given height levels.
!
!
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE Ereport_mod,            ONLY: Ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Input
INTEGER, INTENT(IN)  :: interp_method
INTEGER, INTENT(IN)  :: extrap_method
INTEGER, INTENT(IN)  :: ndata
INTEGER, INTENT(IN)  :: n_out

REAL, INTENT(IN)     :: zdata(ndata)
REAL, INTENT(IN)     :: fdata(ndata)
REAL, INTENT(IN)     :: z_out(n_out)

! Output
REAL, INTENT(OUT)    :: f_out(n_out)

! Local
INTEGER              :: k, kdata, k_ex1, k_ex2
LOGICAL              :: extrap_up(n_out), extrap_down(n_out)
REAL                 :: xi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                           :: ErrorStatus
CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='PROF_INTERP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If there's only one data point then this is the interpolated value
! everywhere. NB. algorithm below would still work in this case but
! that isn't immediately obvious.
IF (ndata == 1) THEN

  f_out(:) = fdata(1)

ELSE

! Interpolate fdata in the region zdata(1) <= z_out < zdata(ndata).
! Points outside this region are flagged for extrapolation.
  extrap_up(:)=.FALSE.
  extrap_down(:)=.FALSE.

  DO k=1,n_out
    IF (z_out(k) < zdata(1)) THEN
      extrap_down(k)=.TRUE.
      CYCLE
    ELSE IF (z_out(k) >= zdata(ndata)) THEN
      extrap_up(k)=.TRUE.
      CYCLE
    ELSE
      DO kdata=1,ndata-1
        IF (z_out(k) < zdata(kdata+1)) EXIT
      END DO
    END IF
    SELECT CASE(interp_method)
    CASE(interp_constant)
      f_out(k)=fdata(kdata)
    CASE(interp_linear)
      xi=(z_out(k)-zdata(kdata))/(zdata(kdata+1)-zdata(kdata))
      f_out(k)=(1.0-xi)*fdata(kdata)+xi*fdata(kdata+1)
    CASE DEFAULT
      ErrorStatus = 1
      Cmessage    = 'Unknown profile interpolation method.'
      CALL Ereport(RoutineName, ErrorStatus, Cmessage)
    END SELECT
  END DO

! Extrapolate where necessary
  DO k=1,n_out

    IF (extrap_down(k)) THEN
      k_ex1=1
      k_ex2=2
    ELSE IF (extrap_up(k)) THEN
      k_ex1=ndata
      k_ex2=ndata-1
    ELSE
      CYCLE
    END IF

    SELECT CASE(extrap_method)
    CASE(extrap_constant)
      f_out(k)=fdata(k_ex1)
    CASE(extrap_linear)
      xi=(z_out(k)-zdata(k_ex1))/(zdata(k_ex2)-zdata(k_ex1))
      f_out(k)=(1.0-xi)*fdata(k_ex1)+xi*fdata(k_ex2)
    CASE(extrap_zero)
      f_out(k)=0.0
    CASE(extrap_linear_to_zero)
      IF (extrap_up(k)) THEN
        k_ex2=n_out
      ELSE
        k_ex2=1
      END IF
      xi=(z_out(k)-zdata(k_ex1))/(z_out(k_ex2)-zdata(k_ex1))
      f_out(k)=(1.0-xi)*fdata(k_ex1)
    CASE DEFAULT
      ErrorStatus = 1
      WRITE(Cmessage, '(''Unknown extrap_method = '', I0)') extrap_method
      CALL Ereport(RoutineName, ErrorStatus, Cmessage)
    END SELECT

  END DO

END IF ! check if only single data point

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE prof_interp

END MODULE prof_interp_mod
