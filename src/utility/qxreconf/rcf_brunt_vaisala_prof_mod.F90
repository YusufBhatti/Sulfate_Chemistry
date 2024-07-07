! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Idealised Brunt Vaisala profile setup

MODULE rcf_brunt_vaisala_prof_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = 'RCF_BRUNT_VAISALA_PROF_MOD'

CONTAINS

! Description:
!   Sets up an initial theta profile using Brunt Vaisala frequencies 
!   The profile can have a isoentropic boundary layer above the surface 
!   with up to two different stabliilties above this given by Brunt-Vaisala 
!   frequencies.  
!
! Method:
!       Uses data in theta_init % height(:) and theta_init % vprof(:,1)
!       where the assumptions is that:
!          theta_init % height(1) is the height up to which a constant 
!          isentropic layer is required.
!          theta_init % height(2:k) give the heights up to which the
!          Brunt-Vaisala frequencies are to be applied.
!          theta_init % vprof(1,1) gives the surface theta
!          theta_init % vprof(2:k,1) give Brunt-Vaisala frequencies.
!          Note if theta_init % height(1)=0.0 then there is no isentropic layer.
!   See diagram below
!
!    ---------------- Model top 
!               |     Another isothermal layer if theta_init % height(3)< top
!    ---------------- theta_init % height(3)  
!               /
!              /     Brunt-Vaisala frequency = theta_init % vprof(3:1) (in /s)
!             /
!    ---------------- theta_init % height(2)   
!             \
!              \
!               \     Brunt-Vaisala frequency = theta_init % vprof(2:1) (in /s)
!                \
!    ---------------- theta_init % height(1)   
!                |
!                |   Isothermal layer  theta_init % vprof(1,1) (in K)
!    ----------------   
!    ///////////////   surface
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 2003
!   This code is written to UMDP3 programming standards.


SUBROUTINE rcf_brunt_vaisala_prof(model_levels, ndata,                   &
                                  zdata, fdata, z_theta,                 &
                                  theta)

USE planet_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage

IMPLICIT NONE

INTEGER, INTENT(IN)  ::   &
  model_levels            & ! Number of model levels
 ,ndata                     ! Number of input information levels

REAL, INTENT(IN)     ::   &
  zdata(ndata)            & ! Heights at which profile changes
 ,fdata(ndata)            & ! Information required to construct profile
                            ! between heights (m). 
                            ! First value theta at surface (K)
                            ! rest different Brunt-Vaisala frequency (N) in (/s)
 ,z_theta(0:model_levels)   ! Height above surface of theta levels (m)

REAL, INTENT(OUT)   ::    &
  theta(0:model_levels)     ! Initial Theta profile on model levels (K)


!---------------------------------------------------------------------------
! Local variables

INTEGER :: j,k              ! loop counter

REAL ::                  &
  dz                        ! distnace between theta level (m)

REAL ::                  &
  BV_squared_over_g(ndata-1) ! Brunt_Vaisala * Brunt_Vaisala / g  (/m)


CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_BRUNT_VAISALA_PROF'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! Construct initial profile based on namelist input data.

DO k=1,ndata-1
  BV_squared_over_g(k) = fdata(k+1)*fdata(k+1)/g
END DO

WRITE(umMessage,'(A)') 'Brunt-Vaisala initial profile option'
CALL umPrint(umMessage,src=RoutineName)

DO k= 0, model_levels

  ! Isentropic part of profile - none if zdata(0)=0.0

  IF (z_theta(k) <= zdata(1) ) THEN
    theta(k) = fdata(1)                 ! constant theta
  ELSE IF  (z_theta(k) <= zdata(ndata)) THEN
    dz = z_theta(k) - z_theta(k-1)
    ! Constant Brunt-Vaisala  dtheta/dz = theta *N*N/g
    ! Expect in most cases the loop to be over no more than two values
    ! i.e. a tropospheric value and a stratospheric value.
    DO j= 1,ndata-1
      IF (z_theta(k) > zdata(j) .AND. z_theta(k) <= zdata(j+1)) THEN            
        theta(k) = theta(k-1) * EXP( BV_squared_over_g(j)*dz)
      END IF
    END DO

  ELSE
    ! Allow another isentropic layer at the top if zdata(ndata) is not
    ! greater than or equal to the model top.  
    theta(k) = theta(k-1)
  END IF

END DO

WRITE(umMessage,'(A)')    ' k z_theta   theta' 
CALL umPrint(umMessage,src=RoutineName)
DO k= 0, model_levels
  WRITE(umMessage,'(i3,f10.2, f10.2)') k, z_theta(k), theta(k)
  CALL umPrint(umMessage,src=RoutineName)
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_brunt_vaisala_prof
END MODULE rcf_brunt_vaisala_prof_mod
