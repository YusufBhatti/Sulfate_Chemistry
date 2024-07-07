
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE SET_RELAX_CF -------------------------------------------
!
!    Purpose :
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE set_relax_cf_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_RELAX_CF_MOD'

CONTAINS

SUBROUTINE set_relax_cf (jgroup, n_rows_local,                    &
                       relax_cf_local, lwind,                     &
           timestep, timestep_no, iter_no,                        &
           icode, cmessage)

USE umPrintMgr
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE UM_ParParams, ONLY: pnorth, psouth
USE field_types, ONLY: fld_type_p
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: jgroup               !IN identifier for group
INTEGER :: n_rows_local
REAL :: relax_cf_local(n_rows_local)
LOGICAL :: lwind                !IN switch to identify wind grid
REAL :: timestep             !IN timestep
INTEGER :: timestep_no          !IN timestep number
INTEGER :: iter_no              !IN iteration count
INTEGER :: icode                !OUT error code and message
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: n_rows
REAL :: relax_cf(rows_max)
INTEGER :: jrow,jrow_print,jrow_print_end
INTEGER :: ir1,ir2,ir3,ir4,ir5,inr1,inr2,inr3,inr4,inr5
REAL :: troplatn, troplats
REAL :: gradn,grads,nudgec
REAL :: nudge_nh
REAL :: nudge_tr
REAL :: nudge_sh
REAL :: nudge_lam

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_RELAX_CF'


!-----------------------------------------------------------------------
!     RELAXATION COEFFICIENTS
!     -----------------------
!
!     these are now input as nudging (N in 3.23 of TN27) and
!     converted to relaxation coefficients in RELAX_CF (using 3.23)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_rows=glsize(2,fld_type_p)

IF (model_type == mt_global) THEN
  nudge_nh  = def_nudge_nh (group_index(jgroup))
  nudge_tr  = def_nudge_tr (group_index(jgroup))
  nudge_sh  = def_nudge_sh (group_index(jgroup))

  DO jrow=1,n_rows
    relax_cf(jrow) = nudge_sh
  END DO
ELSE

  nudge_lam = def_nudge_lam(group_index(jgroup))
  DO jrow=1,n_rows
    relax_cf(jrow) = nudge_lam
  END DO
END IF  ! if GLOBAL

!     convert the nudging coefficients to relaxation coeffs
!     using eqn 3.23 of TN27 (lamda=N*dT/(1+N*dT) )

DO jrow=1,n_rows

  relax_cf(jrow) =                                                &
  relax_cf(jrow)*timestep/(1+relax_cf(jrow)*timestep)

END DO
IF (n_rows+1  <=  rows_max) THEN
  relax_cf(n_rows+1)=relax_cf(n_rows)
END IF

IF (ldiagac .AND. timestep_no == 1 .AND. iter_no == 1) THEN
  IF (mype == 0) THEN

    CALL umPrint(' RELAX_CF for Group No'//TRIM(str(jgroup)),       &
        src='set_relax_cf')
    DO jrow_print=1,n_rows,10
      jrow_print_end=MIN(jrow_print+9,n_rows)
      WRITE(umMessage,'(10F10.3)') (relax_cf(jrow),                &
       jrow=jrow_print,jrow_print_end)
      CALL umPrint(umMessage,src='set_relax_cf')
    END DO
  END IF
END IF
IF (at_extremity(PSouth)) THEN
  relax_cf_local(1)=relax_cf(datastart(2))
ELSE
  relax_cf_local(1)=relax_cf(datastart(2)-1)
END IF
DO jrow=2,n_rows_local
  relax_cf_local(jrow)=relax_cf(jrow+datastart(2)-2)
END DO
IF (at_extremity(PNorth)) THEN
  relax_cf_local(n_rows_local)=                                   &
    relax_cf(n_rows_local+datastart(2)-3)
ELSE
  relax_cf_local(n_rows_local)=                                   &
    relax_cf(n_rows_local+datastart(2)-2)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_relax_cf
END MODULE set_relax_cf_mod
