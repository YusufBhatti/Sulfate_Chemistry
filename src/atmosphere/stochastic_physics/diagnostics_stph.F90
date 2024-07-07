! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic physics
MODULE diagnostics_stph_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_STPH_MOD'

CONTAINS


SUBROUTINE diagnostics_stph( row_length, rows, model_levels,            &
                             n_rows, at_extremity, stph_diag,           &
                             stashwork35)

! Purpose:
!  Calculates diagnostics generated from stochastic physics routines
!  (UM section 35).

! Method:
! Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the stochastic
! physics routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.

!  Diagnostics currently available:

USE stph_diag_mod,  ONLY: strstphdiag

USE um_parparams, ONLY: nodomain, pnorth, peast, psouth, pwest
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE Field_Types
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                              &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
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
, n_rows                                                                &
                   ! number of rows in a v field
, model_levels
                   ! number of model levels

!     Declaration of Stochastic Physics diagnostics.
TYPE (strstphdiag) :: stph_diag


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
PARAMETER ( RoutineName='DIAGNOSTICS_STPH')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------

! Initialise error status
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0
im_index = 1

! ----------------------------------------------------------------------
! DIAG.35001 Copy U wind after SKEB2 to stashwork
! ----------------------------------------------------------------------
item = 1  ! U wind after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_u,                                               &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 1)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35002 Copy V wind after SKEB2 to stashwork
! ----------------------------------------------------------------------

item = 2  ! V wind after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_v,                                               &
       row_length,n_rows,model_levels,0,0,0,0, at_extremity,            &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 2)"//cmessage
  END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35003 Copy SKEB2: Full u increment to stashwork
! ----------------------------------------------------------------------
item = 3  ! skeb2 full u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_u_incr,                                          &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 3)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35004 Copy SKEB2: full V INCR to stashwork
! ----------------------------------------------------------------------

item = 4  ! skeb2 full v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_v_incr,                                          &
       row_length,n_rows,model_levels,0,0,0,0, at_extremity,            &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 4)"//cmessage
  END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35005 Copy SKEB2: Rotational U INCR to stashwork
! ----------------------------------------------------------------------
item = 5  ! skeb2 rotational u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_u_rot,                                           &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 1)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35006 Copy SKEB2: Rotational V INCR to stashwork
! ----------------------------------------------------------------------

item = 6  ! skeb2 rotational  v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_v_rot,                                           &
       row_length,n_rows,model_levels,0,0,0,0, at_extremity,            &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 2)"//cmessage
  END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35007 Copy SKEB2: Divergent U INCR to stashwork
! ----------------------------------------------------------------------
item = 7  ! skeb2 divergent u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_u_div,                                           &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 3)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35008 Copy SKEB2: Divergent V INCR to stashwork
! ----------------------------------------------------------------------

item = 8  ! skeb2 divergent v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_v_div,                                           &
       row_length,n_rows,model_levels,0,0,0,0, at_extremity,            &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 4)"//cmessage
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.35009 Copy SKEB2: dissipation field from smagorinsky code to stashwork
! ----------------------------------------------------------------------
item = 9  ! skeb2 dissipation field from smagorinsky
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_disp_smag,                                       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 1)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35010 Copy SKEB2: dissipation field from convection to stashwork
! ----------------------------------------------------------------------

item = 10  ! skeb2 dissipation field from convection
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_disp_conv,                                       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 2)"//cmessage
  END IF

END IF  !  sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.35011 Copy SKEB2: dissipation field from SKEB1 to stashwork
! ----------------------------------------------------------------------
item = 11  ! skeb2 dissipation field from SKEB1-type
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_disp_skeb1,                                      &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 3)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35012 Copy SKEB2: smoothed modulating field to stashwork
! ----------------------------------------------------------------------

item = 12  ! skeb2 smoothed modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_smodfield,                                       &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35013 Copy SKEB2: final stream function field
! ----------------------------------------------------------------------
item = 13  ! skeb2 raw modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_streamfunction,                                  &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 3)"//cmessage
  END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35014 Copy SKEB2: intial random pattern
! ----------------------------------------------------------------------

item = 14  ! skeb2 smoothed modulating field
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork35(si(item,sect,im_index)),                &
       stph_diag%skeb2_random_pattern,                                  &
       row_length,rows,model_levels,0,0,0,0, at_extremity,              &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35015 Copy SKEB2: Vert Integ. KE of initial SF forcing
! ----------------------------------------------------------------------

item = 15  ! Vert Integ. KE of initial SF forcing
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_psif,                                         &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35016 Copy SKEB2: Vert Integ. KE of numerical diss
! ----------------------------------------------------------------------

item = 16  ! Vert Integ. KE of numerical diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_sdisp,                                        &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35017 Copy SKEB2: Vert Integ. KE of convection diss
! ----------------------------------------------------------------------

item = 17  ! Vert Integ. KE of convection diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_cdisp,                                        &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35018 Copy SKEB2: Vert Integ. KE of 2nd convection diss
! ----------------------------------------------------------------------

item = 18  ! Vert Integ. KE of mflx-based "w" convection diss
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_kdisp,                                        &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35019 Copy SKEB2: Vert Integ. KE of modulated SF forcing
! ----------------------------------------------------------------------

item = 19  ! Vert Integ. KE of modulated SF forcing
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_m_psif,                                       &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35020 Copy SKEB2: Vert Integ. KE of total wind incr before SKEB2
! ----------------------------------------------------------------------

item = 20  ! Vert Integ. KE of total wind incr before SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_prewindincr,                                  &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35021 Copy SKEB2: Vert Integ. KE of wind incr from SKEB2
! ----------------------------------------------------------------------

item = 21  ! Vert Integ. KE of wind incr from SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_windincr,                                     &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.35022 Copy SKEB2: Vert Integ. KE of total wind incr after SKEB2
! ----------------------------------------------------------------------

item = 22  ! Vert Integ. KE of total wind incr after SKEB2
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork35(si(item,sect,im_index)),                   &
       stph_diag%skeb2_ke_postwindincr,                                 &
       row_length,rows,0,0,0,0, at_extremity,                           &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag(item 4)"//cmessage
  END IF
END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
!  single point error handling.
! ----------------------------------------------------------------------

IF (icode /= 0) THEN

  CALL ereport(routinename,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_stph
END MODULE diagnostics_stph_mod
