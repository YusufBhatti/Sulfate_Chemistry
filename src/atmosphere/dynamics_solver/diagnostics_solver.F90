! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE diagnostics_solver_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_SOLVER_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE Diagnostics_solver(row_length,rows,n_rows,STASHwork)

USE atm_fields_bounds_mod, ONLY : udims, vdims, tdims, wdims

USE solver_increments_mod, ONLY : solv_inc_t, solv_inc_u,               &
                                  solv_inc_v, solv_inc_w

USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                              &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
!   Diagnostics_solver extracts diagnostics of (N+1) time-level estimates
!   of primary fields after solver has been called, to be processed
!   by STASH routines for UM section 10 (solver)).
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays from stash_array_mod
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (185,10) u increment          = delta(R_u) across solver
!   (186,10) v increment          = delta(R_v) across solver
!   (187,10) w increment          = delta(w  ) across solver
!   (181,10) T increment          = delta(theta)/exner
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics solver
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::                                                        &
 row_length,rows                                                  &
                  ! horizontal dimensions
,n_rows
                  ! rows for last (N) row of pes

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
REAL :: STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAGNOSTICS_SOLVER'

INTEGER, PARAMETER :: sect = 10      ! STASH section for diagnostics

! Local scalars:
INTEGER ::                                                        &
 im_index                                                         &
             !  internal model index for STASH arrays
,item                                                             &
             !  STASH item of diagnostic
,Errorstatus !  Error status

CHARACTER(LEN=errormessagelength) :: CMessage !  Error message

! Local dynamic arrays:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header

!
! 1. Initialisation
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
im_index    = 1
Errorstatus = 0
Cmessage    = ''

!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! u wind increment
item = 185
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        solv_inc_u,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! v wind increment
item = 186
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        solv_inc_v,                                               &
        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! w wind increment
item = 187
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        solv_inc_w,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! T increment
item = 181
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        solv_inc_t,                                               &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! 3. Error handling
!
IF (Errorstatus /= 0) CALL Ereport(RoutineName,Errorstatus,Cmessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_solver

END MODULE diagnostics_solver_mod
