! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Ensure conv.cloud cover & liquid water path zero if no conv.cloud.
!
! Subroutine Interface:
SUBROUTINE init_cnv(                                              &
                    icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: tdims
USE atm_fields_mod, ONLY: cca, ccb, cct, cclwp, lcbase
USE cv_run_mod, ONLY: l_param_conv, l_ccrad
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: n_cca_lev
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Description:
!   Sets conv.cloud cover and liquid water path to zero at timestep 0
!   if conv.cloud base and top are zero (no conv.cloud is present).
!
! Method:
!   For full model field tests conv.cloud base and top, and if either
!   is zero sets conv.cloud cover and liquid water path to zero.
!   This consistency check at timestep 0 is needed as interpolation in
!   the reconfiguration can give rise to inconsistent fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations: these are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!

! Subroutine arguments
!   Array  arguments with Intent(Out):

!   ErrorStatus <Delete this & the next 2 lines if ErrorStatus not used>
INTEGER ::   icode                ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: cmessage        ! Error message if ICODE>0

! Local parameters:
INTEGER :: int_val = 1      ! used for transfers
! Local scalars:
INTEGER :: i,k,j ! Loop counters over tdims, n_cca_lev

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_CNV'

!- End of header

! 1.0 Ensure that conv.cloud cover and liquid water path are zero
!      when there is zero conv.cloud base and top.
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
WRITE(umMessage,'(A)')'INIT_CNV:resets conv.cld cover zero for base/top zero'
CALL umPrint(umMessage,src='init_cnv')
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    IF ( TRANSFER(ccb(i,j), int_val) == 0 .OR.       &
         TRANSFER(cct(i,j), int_val) == 0 ) THEN
      DO k = 1, n_cca_lev
        cca(i,j,k) = 0.0
      END DO
      cclwp(i,j) = 0.0
    END IF
  END DO ! i
END DO ! j

! Make sure lowest cloud base zeroed if no convection but L_CCrad is true
! Note the radiation scheme uses the lcbase array in this case but the
! convection scheme is not setting a value.

IF (.NOT.l_param_conv .AND. l_ccrad) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      lcbase(i,j) = 0
    END DO ! i
  END DO ! j
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_cnv

