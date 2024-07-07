! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected AND a full boundary layer treatment is not
!           performed.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer

MODULE bl_lsp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BL_LSP_MOD'
CONTAINS

SUBROUTINE bl_lsp( bl_levels,qcf,q,t )

USE atm_fields_bounds_mod, ONLY: tdims
USE planet_constants_mod, ONLY: lsrcp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
  bl_levels             ! IN   Number of boundary layer levels

REAL, INTENT(INOUT) ::                                                  &
  qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      bl_levels),                                                       &
                                 ! INOUT Ice water content
  q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
    bl_levels),                                                         &
                                 ! INOUT
!                                  IN    Vapour+liquid+ice content
!                                  OUT   Vapour+liquid content
    t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      bl_levels)                   ! INOUT
!                                  IN    Liquid ice temperature
!                                  OUT   Liquid temperature
! Temporary Space
INTEGER ::                                                              &
        i,                                                              &
                               ! Counter over points
        j,                                                              &
                               ! Counter over points
        k                ! Counter over boundary layer levels
REAL :: newqcf              ! Temporary variable for QCF

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BL_LSP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                        &
!$OMP&         PRIVATE(i,j,k,newqcf)                                    &
!$OMP&         SHARED(bl_levels,tdims,q,qcf,t,lsrcp)
DO k = 1, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      ! Convert Q (vapour+liquid+ice) to (vapour+liquid)
      q(i,j,k)=q(i,j,k)-qcf(i,j,k)
      ! Check that Q is not negative
      IF (q(i,j,k)  <   0.0) THEN
        ! Evaporate ice to keep Q positive, but don't let ice go negative
        ! itself
        newqcf=MAX(qcf(i,j,k)+q(i,j,k),0.0)
        q(i,j,k)=q(i,j,k)+(qcf(i,j,k)-newqcf)
        qcf(i,j,k)=newqcf
      END IF
      ! Adjust T from T liquid ice to T liquid
      t(i,j,k)=t(i,j,k)+lsrcp*qcf(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
! End the subroutine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bl_lsp
END MODULE bl_lsp_mod
