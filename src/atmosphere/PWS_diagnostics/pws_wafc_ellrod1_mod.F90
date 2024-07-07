! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_wafc_ellrod1_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_WAFC_ELLROD1_MOD'

CONTAINS

SUBROUTINE pws_wafc_ellrod1(UStdLev,     &  ! in 1,2
                         VStdLev,     &  ! in
                         dUdX,        &  ! in
                         dUdY,        &  ! in
                         dVdX,        &  ! in
                         dVdY,        &  ! in
                         dWdZ,        &  ! in
                         CATPred)
! Description
!  Routine to calculate Ellrods 1st CAT index TI1.
! Method:
!   1) Calculate deformation and vertical wind shear.
!   2) Use Ellrod's empirical formula to calculate the index:
!        TI1 = vert. wind shear * Deformation
!        (see Ellrod and Knapp, 1991)
!   3) Convert TI1 values to percentage probability using  equation
!         prob_TI1 = (1383149.3*TI1) + 2.7258455
!        which was obtained by correlating TI1 with CAT probability.
!   4) Lastly, set all values < 3.5%  to zero.

USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims
USE trignometric_mod, ONLY: sin_v_latitude
USE missing_data_mod, ONLY: rmdi


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:
REAL , INTENT(IN) :: UStdLev(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end)   ! U & V wind comps on
REAL , INTENT(IN) :: VStdLev(udims%i_start:udims%i_end, &
                             vdims%j_start:vdims%j_end)   ! std level
REAL , INTENT(IN) :: dUdX(udims%i_start:udims%i_end, &
                          vdims%j_start:vdims%j_end)
REAL , INTENT(IN) :: dUdY(udims%i_start:udims%i_end, &
                          vdims%j_start:vdims%j_end)
REAL , INTENT(IN) :: dVdX(udims%i_start:udims%i_end, &
                          vdims%j_start:vdims%j_end)
REAL , INTENT(IN) :: dVdY(udims%i_start:udims%i_end, &
                          vdims%j_start:vdims%j_end)
REAL , INTENT(IN) :: dWdZ(udims%i_start:udims%i_end, &
                          vdims%j_start:vdims%j_end)

REAL , INTENT(OUT) :: catpred(udims%i_start:udims%i_end, &
                              vdims%j_start:vdims%j_end)

! Local variables
INTEGER :: i, j                  ! Loop counters
REAL    :: ShearV                ! Vertical wind shear   (VWS)
REAL    :: St_Def                ! Stretching deformation
REAL    :: Sh_Def                ! Shearing deformation
REAL    :: Def                   ! Total deformation
REAL    :: ti1                   ! Ellrod turbulence indicator TI1
REAL    :: prob_TI1              ! TI1 converted to percentage probability
                                 !   of encountering turbulence


CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_WAFC_ELLROD1'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------
! Calculate shearing & stretching deformation and TI1 for each point
! Then convert TI1 to percentage probability.

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end

    IF ( (dUdX(i,j) == rmdi) .OR.  &
         (dUdY(i,j) == rmdi) .OR.  &
         (dVdX(i,j) == rmdi) .OR.  &
         (dVdY(i,j) == rmdi) ) THEN
      CATPred(i,j) = rmdi
    ELSE

      ShearV = dWdZ(i,j)

      !Calculate shearing & stretching deformation, then TI1
      St_Def = dUdX(i,j) - dVdY(i,j)
      Sh_Def = dVdX(i,j) + dUdY(i,j)
      Def = SQRT( (St_Def**2) + (Sh_Def**2) )
      ti1 = ShearV * Def

      !Convert to percentage probability
      prob_TI1 = (1383149.3*ti1) + 2.7258455
      IF (prob_TI1 < 3.5) THEN
        prob_TI1=0.0
      END IF

      CATPred(i,j) = prob_TI1

    END IF
  END DO
END DO



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_wafc_ellrod1

END MODULE pws_wafc_ellrod1_mod
