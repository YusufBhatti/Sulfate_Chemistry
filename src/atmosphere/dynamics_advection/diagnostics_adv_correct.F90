! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE diag_adv_correct_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAG_ADV_CORRECT_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE Diagnostics_adv_correct(row_length, rows, stashwork)

USE atm_fields_bounds_mod, ONLY : tdims

USE adv_correct_incs_mod, ONLY :                                         &
    advcor_inc_t, advcor_inc_m_v, advcor_inc_m_cl, advcor_inc_m_cf,      &
    advcor_inc_m_r, advcor_inc_m_gr, advcor_inc_m_cf2, advcor_inc_q,     &
    advcor_inc_qcl, advcor_inc_qcf, advcor_inc_qrain, advcor_inc_qgraup, &
    advcor_inc_qcf2

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE nlsizes_namelist_mod, ONLY: model_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE submodel_mod, ONLY: atmos_im, atmos_sm
USE stash_array_mod, ONLY:                                               &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf, &
    stash_maxlen
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
!   Diagnostics_adv_correct extracts 
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays from stash_array_mod
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (381,12) T increment          = delta(T)   across adv correction
!   (382,12) q  increment         = delta(q)   across adv correction
!   (383,12) qcl increment        = delta(qcl) across adv correction
!   (384,12) qcf increment        = delta(qcf) across adv correction
!   (389,12) qrain increment      = delta(qrain) across adv correction
!   (390,12) qgraup increment     = delta(qgraup) across adv correction
!   (391,12) qcf2 increment       = delta(q)   across adv correction
!   (395,12) mv  increment        = delta(mv)  across adv correction
!   (396,12) mcl increment        = delta(mcl) across adv correction
!   (397,12) mcf increment        = delta(mcf) across adv correction
!   (398,12) mrain increment      = delta(mrain) across adv correction
!   (399,12) mgruap increment     = delta(mgraup) across adv correction
!   (400,12) mcf2 increment       = delta(mcf2) across adv correction
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to the UM programming standards in UMDP3
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: row_length,rows

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
REAL, INTENT(INOUT) :: STASHwork(*)

!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):

! Local parameters:
CHARACTER(LEN=*), PARAMETER :: RoutineName='DIAGNOSTICS_ADV_CORRECT'

INTEGER, PARAMETER :: sect = 12      ! STASH section for diagnostics

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local scalars:
INTEGER ::                                                        &
 im_index                                                         &
             !  internal model index for STASH arrays
,item                                                             &
             !  STASH item of diagnostic
,Errorstatus !  Error status

CHARACTER(LEN=errormessagelength) :: CMessage !  Error message

! Local dynamic arrays:

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

! Temperature increment
item = 381
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_t,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! q increment
item = 382
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_q,                                             &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcl increment
item = 383
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_qcl,                                           &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf increment
item = 384
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_qcf,                                           &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qrain increment
item = 389
IF (l_mcr_qrain .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_qrain,                                         &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qgraup increment
item = 390
IF (l_mcr_qgraup .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_qgraup,                                        &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! qcf2 increment
item = 391
IF (l_mcr_qcf2 .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_qcf2,                                          &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_v increment
item = 395
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_v,                                           &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cl increment
item = 396
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_cl,                                          &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cf increment
item = 397
IF (sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_cf,                                          &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_r increment
item = 398
IF (l_mcr_qrain .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_r,                                           &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_gr increment
item = 399
IF (l_mcr_qgraup .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_gr,                                          &
        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,sect,item,                                       &
        Errorstatus,cmessage)

END IF ! sf(item,sect)

! m_cf2 increment
item = 400
IF (l_mcr_qcf2 .AND. sf(item,sect) .AND. Errorstatus == 0) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
        advcor_inc_m_cf2,                                         &
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
END SUBROUTINE diagnostics_adv_correct

END MODULE diag_adv_correct_mod
