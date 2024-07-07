! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE diag_R_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAG_R_MOD'

CONTAINS

SUBROUTINE diag_R()

USE atm_fields_bounds_mod
USE um_parcore, ONLY: mype, nproc
USE fields_rhs_mod,  ONLY:  r_u ,r_v
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE umPrintMgr

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


REAL :: max_r_u,min_r_u,max_r_v,min_r_v

INTEGER :: ierr
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAG_R'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)


max_r_u     = MAXVAL(R_u(udims%i_start:udims%i_end,                   &
                         udims%j_start:udims%j_end,                   &
                         udims%k_start:udims%k_end))
min_r_u     = MINVAL(R_u(udims%i_start:udims%i_end,                   &
                         udims%j_start:udims%j_end,                   &
                         udims%k_start:udims%k_end))
max_r_v     = MAXVAL(R_v(vdims%i_start:vdims%i_end,                   &
                         vdims%j_start:vdims%j_end,                   &
                         vdims%k_start:vdims%k_end))
min_r_v     = MINVAL(R_v(vdims%i_start:vdims%i_end,                   &
                         vdims%j_start:vdims%j_end,                   &
                         vdims%k_start:vdims%k_end))

CALL gc_rmax(1,nproc,ierr,max_r_u)
CALL gc_rmax(1,nproc,ierr,max_r_v)
CALL gc_rmin(1,nproc,ierr,min_r_u)
CALL gc_rmin(1,nproc,ierr,min_r_v)

IF (mype==0) THEN
  WRITE(umMessage,FMT='(2A)') '============================================',  &
                        '========================================'
  CALL umPrint(umMessage,src='diag_R')

  WRITE(umMessage,FMT='(A)') 'Slow physics source terms from atmos_physics1:'
  CALL umPrint(umMessage,src='diag_R')
  WRITE(umMessage,FMT='(A,2E32.16)') 'r_u      :',min_r_u,max_r_u
  CALL umPrint(umMessage,src='diag_R')
  WRITE(umMessage,FMT='(A,2E32.16)') 'r_v      :',min_r_v,max_r_v
  CALL umPrint(umMessage,src='diag_R')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE
END MODULE
