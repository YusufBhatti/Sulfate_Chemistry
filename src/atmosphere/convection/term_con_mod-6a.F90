! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns a mask for points at which convection is terminating
!

MODULE term_con_6a_mod

IMPLICIT NONE

!
! Description:
!   Returns a mask for points at which convection is terminating!
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TERM_CON_6A_MOD'

CONTAINS

SUBROUTINE term_con_6a(npnts,nlev,k,flxkp1,flx_init,deltak,bterm)

USE water_constants_mod, ONLY: lc, lf
USE planet_constants_mod, ONLY: cp, c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations
INTEGER,INTENT(IN) :: k             ! present model layer

REAL,INTENT(IN) :: flxkp1(npnts)    ! parcel massflux in layer k+1 (Pa/s)
REAL,INTENT(IN) :: flx_init(npnts)  ! Initial par. massflux at cloud base (Pa/s)
REAL,INTENT(IN) :: deltak(npnts)    ! Parcel forced detrainment rate in
                                    ! layer k multiplied by layer thickness

!---------------------------------------------------------------------
! Variables which are input and output
!---------------------------------------------------------------------
LOGICAL,INTENT(INOUT) :: bterm(npnts)   ! Mask for parcels which terminate
                                        ! in layer k+1

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::   i                 ! loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TERM_CON_6A'

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
!  Convection will terminate when the mass flux falls below a small
!  fraction of the cloud base mass flux or if the forced detrainment
!  rate is large or if approaching the top of the model.
!  It will also terminate if bterm has been set to true at a higher
!  level.
!----------------------------------------------------------------------

DO i=1,npnts
  bterm(i) =  bterm(i)                         .OR. &
             (flxkp1(i)  <   0.05*flx_init(i)) .OR. &
             (deltak(i)  >=  0.95)             .OR. &
             (k+1        ==  nlev)
END DO  ! i

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE term_con_6a

END MODULE term_con_6a_mod
