! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large scale precip
MODULE diagnostics_pc2checks_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_PC2CHECKS_MOD'

CONTAINS

SUBROUTINE diagnostics_pc2checks(                                             &
                       row_length, rows                                       &
,                      me                                                     &
,                      at_extremity                                           &
,                      t_inc, q_inc, qcl_inc, qcf_inc                         &
,                      cf_inc, cfl_inc, cff_inc                               &
,                      stashwork                                  )

! Purpose:
!          Calculates diagnostics and outputs them.

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                    &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE nlsizes_namelist_mod, ONLY: model_levels
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

LOGICAL ::                                                                    &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Parameters

INTEGER ::                                                                    &
  row_length      & ! number of points on a row
, rows              ! number of rows in a theta field

INTEGER ::                                                                    &
  me                  ! IN. Processor number

! Primary Arrays used in all models
REAL ::                                                                       &
  t_inc(row_length, rows, model_levels)                                       &
, q_inc(row_length,rows, model_levels)                                        &
, qcl_inc(row_length, rows, model_levels)                                     &
, qcf_inc(row_length, rows, model_levels)                                     &
, cf_inc(row_length,rows, model_levels)                                       &
, cfl_inc(row_length, rows, model_levels)                                     &
, cff_inc(row_length, rows, model_levels)

! 3d work array for calculating separate +/- PC2 increments
! only to be allocated if it is needed.
REAL, ALLOCATABLE ::  work3d(:,:,:)


! Diagnostic variables
REAL ::                                                                       &
 stashwork(*)        ! STASH workspace for section 4 (PC2 checks)

! Local variables
INTEGER ::                                                                    &
 i, j, k                                                                      &
,    icode           ! Return code  =0 Normal exit  >1 Error

INTEGER :: sect,item    ! STASH section, item no.s
PARAMETER (sect = 4) !  for microphysics - large scale rain

CHARACTER(LEN=errormessagelength) :: cmessage

CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName='DIAGNOSTICS_PC2CHECKS')

INTEGER ::                                                                    &
  im_index        ! internal model index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0 ! Initialise error status
im_index = 1

! Allocate work3d array if calculating +/- increments for cfl,cff,qcl,qcf
IF (sf(130,sect) .OR. sf(131,sect) .OR.                                       &
    sf(132,sect) .OR. sf(133,sect) .OR.                                       &
    sf(136,sect) .OR. sf(137,sect) .OR.                                       &
    sf(138,sect) .OR. sf(139,sect) ) THEN
  ALLOCATE ( work3d(row_length,rows,model_levels) )
END IF

! Copy diagnostic information to STASHwork for STASH processing

! increment diagnostics= modified - previous

item = 141  ! temperature increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       t_inc,                                                                 &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 141)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 142  ! vapour increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       q_inc,                                                                 &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 142)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 143  ! liquid increment net
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       qcl_inc,                                                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 143)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 130  ! liquid increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,qcl_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MAX(0.0,qcl_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 130)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 131  ! liquid increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,qcl_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MIN(0.0,qcl_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 131)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 144  ! ice increment: net
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       qcf_inc,                                                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 144)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 132  ! ice increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,qcf_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MAX(0.0,qcf_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 132)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 133  ! ice increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,qcf_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MIN(0.0,qcf_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 133)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 152  ! total cloud fraction increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       cf_inc,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 152)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 153  ! liquid cloud fraction increment: net
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       cfl_inc,                                                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 153)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 136  ! liquid cloud fraction increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,cfl_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MAX(0.0,cfl_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 136)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 137  ! liquid cloud fraction increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,cfl_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MIN(0.0,cfl_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 137)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 154  ! ice cloud fraction increment: net
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       cff_inc,                                                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 154)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 138  ! ice cloud fraction increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,cff_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MAX(0.0,cff_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 138)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 139  ! ice cloud fraction increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work3d,cff_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work3d(i,j,k) = MIN(0.0,cff_inc(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                        &
       work3d,                                                                &
       row_length,rows,model_levels,0,0,0,0, at_extremity,                    &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,                    &
       stash_levels,num_stash_levels+1,                                       &
       atmos_im,sect,item,                                                    &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 139)"//cmessage
  END IF

END IF  !  sf(item,sect)

! Allocate work3 array if calculating +/- increments for cfl,cff,qcl,qcf.
IF (sf(130,sect) .OR. sf(131,sect) .OR.                                       &
    sf(132,sect) .OR. sf(133,sect) .OR.                                       &
    sf(136,sect) .OR. sf(137,sect) .OR.                                       &
    sf(138,sect) .OR. sf(139,sect) ) THEN
  DEALLOCATE ( work3d )
END IF

9999 CONTINUE
IF (icode /= 0) THEN
  CALL ereport(RoutineName,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_pc2checks
END MODULE diagnostics_pc2checks_mod
