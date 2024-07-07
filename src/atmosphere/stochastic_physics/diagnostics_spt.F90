
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic physics
MODULE diagnostics_spt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_SPT_MOD'

CONTAINS


SUBROUTINE diagnostics_spt( row_length, rows,                           &
                             n_rows, at_extremity, spt_diag,            &
                             stashwork35)

! Purpose:
!  Calculates diagnostics generated from SPT on UM section 35.

! Method:
! Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the stochastic
! physics routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.

!  Diagnostics currently available:

USE um_parparams,  ONLY: nodomain, pnorth, peast, psouth, pwest
USE spt_diag_mod,  ONLY: strsptdiag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE Field_Types
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE missing_data_mod
USE nlsizes_namelist_mod, ONLY: model_levels

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

LOGICAL ::                                                              &
  at_extremity(4)
                   ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

INTEGER ::                                                              &
  row_length                                                            &
                   ! number of points on a row
, rows                                                                  &
                   ! number of rows in a theta field
, n_rows
                   ! number of rows in a v field

!     Declaration of Stochastic Physics diagnostics.
TYPE (strsptdiag) :: spt_diag


!  Global Variables:----------------------------------------------------

! Diagnostics info
REAL ::                                                                 &
 stashwork35(*)
                  ! STASH workspace
INTEGER ::                                                              &
  im_index        ! internal model index

! Local variables

INTEGER ::                                                              &
  icode                                                                 &
                  ! Return code  =0 Normal exit  >1 Error
 ,item                                                                  &
                  ! STASH item
 ,sect            ! STASH section
PARAMETER( sect = 35 ) ! for stochastic physics

CHARACTER(LEN=errormessagelength) :: cmessage

CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName='DIAGNOSTICS_SPT')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

! Initialise error status
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0
im_index = 1

! ----------------------------------------------------------------------
! DIAG.35023 SPT Forcing pattern
! ----------------------------------------------------------------------

item = 23  ! Forcing pattern for the SPT scheme.
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_forcing_pattern,                                         &
  row_length,rows,model_levels,0,0,0,0, at_extremity,                   &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 23)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35024 Theta SPT increment
! ----------------------------------------------------------------------

item = 24  ! Theta SPT increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_theta_inc,                                               &
  row_length,rows,model_levels,0,0,0,0, at_extremity,                   &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 24)"//cmessage
  END IF
END IF !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35025 q SPT increment
! ----------------------------------------------------------------------

item = 25  ! q SPT increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_q_inc,                                                   &
  row_length,rows,model_levels,0,0,0,0, at_extremity,                   &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 25)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35026 u SPT increment
! ----------------------------------------------------------------------

item = 26  !u SPT increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_u_inc,                                                   &
  row_length,rows,model_levels,0,0,0,0, at_extremity,                   &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 26)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35027 v SPT increment
! ----------------------------------------------------------------------

item = 27  ! v SPT increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_v_inc,                                                   &
  row_length,n_rows,model_levels,0,0,0,0, at_extremity,                 &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 27)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35028 CFL criteria exceeded
! ----------------------------------------------------------------------
item = 28  ! CFL criteria exceeded
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       spt_diag%cfl_br_marker,                                          &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 28)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35029 T SPT increment
! ----------------------------------------------------------------------

item = 29  ! T SPT increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
  spt_diag%spt_t_inc,                                                   &
  row_length,rows,model_levels,0,0,0,0, at_extremity,                   &
  stlist(1,stindex(1,item,sect,im_index)),len_stlist,                   &
  stash_levels,num_stash_levels+1,                                      &
  atmos_im,sect,item,                                                   &
  icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 29)"//cmessage
  END IF
END IF !  sf(item,sect)
! ----------------------------------------------------------------------
!  single point error handling.
! ----------------------------------------------------------------------

IF (icode /= 0) THEN

  CALL ereport(routinename,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_spt

END MODULE diagnostics_spt_mod
