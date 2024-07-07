! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
! Allocate space for idealised output diagnostic arrays

MODULE alloc_ideal_diag_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   ALLOCATE and initialise diagnostic arrays required by the idealised UM
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UMDP standards.
!------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOC_IDEAL_DIAG_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE alloc_ideal_diag(iflag)

USE idealised_diag_mod, ONLY:                                             &
  dt_inc_ideal_um, dq_inc_ideal_um, du_inc_ideal_um, dv_inc_ideal_um,     &
  dtheta_inc_ideal_um, dcolqdt_ideal_um, l_stored_ref, diag_theta_ref,    &
  diag_exner_ref, diag_rho_ref, diag_u_ref, diag_v_ref, diag_q_ref,       &
  de_cvt_ideal_um, de_u2_ideal_um, de_v2_ideal_um

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY: udims, vdims, tdims, pdims

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: iflag     ! indicates where called from and what to do

!------------------------------------------------------------------------------
! Local variables

INTEGER ::           &
  i, j, k            & ! loop counters
 ,icode                ! Return code  =0 Normal, exit  >1 Error

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_IDEAL_DIAG'


! Variables required for Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
IF (iflag == 0) THEN

  ! Called from initialisation routines - allocate arrays to hold
  ! reference info for diagnostics
  ! Arrays allocated for complete run
  ! ENDGame versions will hold level 0

  l_stored_ref = .FALSE.    ! indicator of whether reference arrays set

  ALLOCATE (diag_theta_ref(tdims%k_start:tdims%k_end))
  ALLOCATE (diag_exner_ref(pdims%k_start:pdims%k_end))
  ALLOCATE (diag_rho_ref(pdims%k_start:pdims%k_end))
  ALLOCATE (diag_u_ref(udims%k_start:udims%k_end))
  ALLOCATE (diag_v_ref(vdims%k_start:vdims%k_end))
  ALLOCATE (diag_q_ref(tdims%k_start:tdims%k_end))

ELSE IF (iflag == 1) THEN
  
  ! Normal calling every timestep before forcing call

  ALLOCATE ( dcolqdt_ideal_um(tdims%i_start:tdims%i_end,    &
                              tdims%j_start:tdims%j_end))
  ALLOCATE ( de_cvt_ideal_um(tdims%i_start:tdims%i_end,    &
                             tdims%j_start:tdims%j_end))
  ALLOCATE ( de_u2_ideal_um(udims%i_start:udims%i_end,    &
                             udims%j_start:udims%j_end))
  ALLOCATE ( de_v2_ideal_um(vdims%i_start:vdims%i_end,    &
                             vdims%j_start:vdims%j_end))

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      de_cvt_ideal_um(i,j)  = 0.0
      dcolqdt_ideal_um(i,j) = 0.0
    END DO
  END DO
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      de_u2_ideal_um(i,j) = 0.0
    END DO
  END DO
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      de_v2_ideal_um(i,j) = 0.0
    END DO
  END DO

  ALLOCATE ( dt_inc_ideal_um(tdims%i_start:tdims%i_end,    &
                             tdims%j_start:tdims%j_end,    &
                             tdims%k_start:tdims%k_end))
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dt_inc_ideal_um(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ALLOCATE ( dq_inc_ideal_um(tdims%i_start:tdims%i_end,    &
                             tdims%j_start:tdims%j_end,    &
                             tdims%k_start:tdims%k_end))
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dq_inc_ideal_um(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ALLOCATE ( du_inc_ideal_um(udims%i_start:udims%i_end,    &
                             udims%j_start:udims%j_end,    &
                             udims%k_start:udims%k_end))

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_inc_ideal_um(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ALLOCATE ( dv_inc_ideal_um(vdims%i_start:vdims%i_end,    &
                             vdims%j_start:vdims%j_end,    &
                             vdims%k_start:vdims%k_end))
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_inc_ideal_um(i,j,k) = 0.0
      END DO
    END DO
  END DO

  ALLOCATE ( dtheta_inc_ideal_um(tdims%i_start:tdims%i_end,    &
                                 tdims%j_start:tdims%j_end,    &
                                 tdims%k_start:tdims%k_end))
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_inc_ideal_um(i,j,k) = 0.0
      END DO
    END DO
  END DO

ELSE

  ! Should never reach here - will cause model to stop
  WRITE(cmessage,'(a,i6)')                                               &
        'Alloc_ideal_diag called with invalid iflag value : ',iflag
  icode = 1
  CALL ereport(routinename,icode,cmessage)
  
END IF ! test on iflag
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE alloc_ideal_diag

END MODULE alloc_ideal_diag_mod
