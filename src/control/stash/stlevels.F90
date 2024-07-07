! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine STLEVELS -----------------------------------------------
!
! Purpose: Generate a level index from STASHrecord and level_lists
!          and number of levels tailored to a particular diagnostic.
!          Also set levels and pseudo-levels information for encoding
!          PPheader details. 
!          New subroutine STLEVELS is based on GEN_INDEX and
!          PP_COMPUTE_LEVEL with merged functionality.
!          A general note as levels list is an integer
!          real values are multiplied by a 1000.0.
!          When computing the real value of the level for the
!          pp header it is necessary to divide by a 1000.0.
!          Levels that are affected by this are theta, pressure and
!          height.
!
! Programming standard: UM Doc Paper 3
!
! Project task: C4
!
! External documentation : UMDP no C4
!
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE stlevels(stash_control,stash_control_size,             &
     stash_levels,num_stash_levels,num_level_lists,               &
     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,       &
     max_stash_levs,num_levs_in,num_levs_out,num_pseudo_out,      &
     index_size,index_lev,level_list,                             &
     lbvcl,level,pseudo_level,                                    &
     icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE sterr_mod, ONLY: nonsense
USE stparam_mod, ONLY: st_input_bottom, st_output_bottom, st_input_top,&
    st_output_top, st_special_code, st_proc_no_code, st_gridpoint_code,&
    block_size, st_time_series_code, st_time_series_mean,              &
    st_append_traj_code, st_pseudo_in, st_pseudo_out, vert_mean_base,  &
    global_mean_base
                       
USE cppxref_mod, ONLY: ppx_lbvc_height, ppx_lbvc_pressure,        &
                       ppx_lbvc_theta, ppx_lbvc_PV
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER ::                                                        &
       stash_control_size,                                        &
                           ! IN size of stash control record
       stash_control(stash_control_size),                         &
                                         ! IN  stash control
       num_stash_levels,                                          &
                           ! IN max. no of hts for a levels list
       num_level_lists,                                           &
                           ! IN max. no of level lists
       stash_levels(num_stash_levels+1,num_level_lists),          &
                           ! IN lookup table for level lists
             num_stash_pseudo,num_pseudo_lists,                         &
                                               ! IN dims of pseudo_levs
             stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists),  &
                                 ! IN lookup table for pseudo-lev lists
             max_stash_levs,                                            &
                                 ! IN max. no of output levels
             num_levs_in,                                               &
                                 ! OUT no of levels in input data
             num_levs_out,                                              &
                                 ! OUT no of levels in output data
             num_pseudo_out,                                            &
                                 ! OUT no of pseudo levels in output data
             index_size,                                                &
                                 ! OUT no of levels in levels index
             index_lev(max_stash_levs),                                 &
                                 ! OUT index of output level rel to input level
             level_list(max_stash_levs),                                &
                                         ! OUT value of model level
             pseudo_level(max_stash_levs),                              &
                                         ! OUT Value of pseudo levels
             lbvcl,                                                     &
                                 ! IN  vertical coordinate PP code
             icode               ! OUT error code
REAL ::                                                           &
       level(max_stash_levs)  ! OUT Value of output levels (real)
CHARACTER(LEN=errormessagelength) ::                                    &
       cmessage            ! OUT error message

! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER ::                                                        &
       index_pseudo_lev(max_stash_levs),                          &
                                         ! Pseudo-level 1D index
       num_pseudo_in,                                             &
                                 ! Number of pseudo levels in input data
       k2,ml,kl,                                                  &
                                 ! loop counts
       ni,no,                                                     &
                                 ! Number In/Out
       indx1,                                                     &
                                 ! index count
       ilev,                                                      &
                                 ! Integer level/pseudo-level
       what_mean,what_proc       ! Meaning and processing code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STLEVELS'
!
! First compute the index for physical levels
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (stash_control(st_input_bottom) <  0) THEN ! Input LEVELS list
  ni=-stash_control(st_input_bottom)
  num_levs_in=stash_levels(1,ni)
  IF (stash_control(st_output_bottom) <  0) THEN ! LEVELS LIST out
    no=-stash_control(st_output_bottom)
    num_levs_out=stash_levels(1,no)
    indx1=0
    DO ml=1,num_levs_out
      ilev=stash_levels(ml+1,no)    !  Level required
      DO kl=1,num_levs_in
        IF (stash_levels(kl+1,ni) == ilev) THEN
          indx1=indx1+1
          index_lev(indx1)=kl   ! Relative position of Input to Ou
          level_list(indx1)=ilev
          GO TO 400
        END IF
      END DO

      icode=nonsense
      WRITE(cmessage,'(A,I6,A)') 'STLEVELS : Output level ',ilev, &
                                 ' not found in input levels list'
      GO TO 999
      400        CONTINUE
    END DO

  ELSE           !  Output as a Level range
    num_levs_out=stash_control(st_output_top)-                    &
                 stash_control(st_output_bottom)+1
    ilev=stash_control(st_output_bottom) !1st output model level
    DO kl=1,num_levs_in
      IF (stash_levels(kl+1,ni) == ilev) THEN
        index_lev(1)=kl ! Relative posn of Input to the 1st level
        level_list(1)=ilev
        GO TO 401
      END IF
    END DO
    icode=nonsense
    WRITE(cmessage,'(A,I6,A)') 'STLEVELS : Output bottom model level ', &
                               ilev, ' not found in input levels list'
    GO TO 999
    401      CONTINUE
    DO kl=2,num_levs_out
      index_lev(kl)=index_lev(kl-1)+1
      level_list(kl)=level_list(kl-1)+1
    END DO
  END IF
ELSE IF (stash_control(st_input_bottom) == 100) THEN !Special level
  num_levs_in=1
  num_levs_out=1
  index_lev(1)=1
  level_list(1)=1 ! could be worth setting to some nonsense no.
ELSE     !  Input is Model level range
  num_levs_in=stash_control(st_input_top)-                        &
              stash_control(st_input_bottom)+1
  IF (stash_control(st_output_bottom) <  0) THEN ! LEVELS LIST out
    no=-stash_control(st_output_bottom)
    num_levs_out=stash_levels(1,no)
    indx1=0
    DO ml=1,num_levs_out
      ilev=stash_levels(ml+1,no)    ! Output level reqd
      DO kl=1,num_levs_in
        IF ((stash_control(st_input_bottom)+kl-1) == ilev) THEN
          indx1=indx1+1
          index_lev(indx1)=kl   ! Relative posn of output to inpt
          level_list(indx1)=ilev
          GO TO 402
        END IF
      END DO
      icode=nonsense
      WRITE(cmessage,'(A,I6,A)') 'STLEVELS : Output model level ', &
                          ilev, ' not in input model level range'
      GO TO 999
      402        CONTINUE
    END DO
  ELSE     !   Output as model level range
    ! Do some consistency checks here to ensure valid processing request
    ! output bottom should be greater or equal to input bottom
    IF (stash_control(st_output_bottom) <                         &
       stash_control(st_input_bottom)) THEN
      icode=nonsense
      WRITE(cmessage,'(A,A,2I5)') 'STLEVELS : >> FATAL ERROR << ', & 
      'bad level spec, bot input>output', stash_control(st_input_bottom), &
       stash_control(st_output_bottom)
      GO TO 999 ! jump to error
    ELSE IF (stash_control(st_output_top) >                        &
         stash_control(st_input_top)) THEN
      icode=nonsense
      WRITE(cmessage,'(A,A,2I5)') 'STLEVELS : >> FATAL ERROR << ', &
      'bad level spec, top input<output', stash_control(st_input_top), &
       stash_control(st_output_top)
      GO TO 999 ! jump to error
    END IF
    num_levs_out=stash_control(st_output_top)-                    &
                 stash_control(st_output_bottom)+1
    index_lev(1)=stash_control(st_output_bottom)-                 &
                 stash_control(st_input_bottom)+1
    level_list(1)=stash_control(st_output_bottom)
    DO kl=2,num_levs_out
      index_lev(kl)=index_lev(kl-1)+1
      level_list(kl)=level_list(kl-1)+1
    END DO
  END IF
END IF
index_size=num_levs_out
IF (num_levs_out >  num_levs_in) THEN   ! things very badly wrong
  icode=nonsense
  WRITE(cmessage,'(A,A,2I5)') 'STLEVELS : >> FATAL ERROR << ', & 
  'asking for num_levs_out>num_levs_in', num_levs_out,num_levs_in
  GO TO 999 ! jump to return
END IF
!
! Next, compute actual (physical) levels for encoding PPheaders
!
IF (stash_control(st_output_bottom) <  0) THEN ! Levels List ?
  no=-stash_control(st_output_bottom)     ! Index of Levels list

    ! Remove scaling (by factor 1000) of vertical level coord
    ! for certain types of STASH output [originally needed to
    ! store in an intermediary integer array]
  IF ( lbvcl  ==  ppx_lbvc_height   .OR.                         &
                                        !  height levels
      lbvcl  ==  ppx_lbvc_pressure .OR.                         &
                                        ! pressure levels
      lbvcl  ==  ppx_lbvc_theta    .OR.                         &
                                        ! theta levels
      lbvcl  ==  ppx_lbvc_PV ) THEN     ! potential vorticity


    DO ml=1,num_levs_out
      level(ml)=REAL(stash_levels(ml+1,no))*0.001+1.0e-10
    END DO
  ELSE
    DO ml=1,num_levs_out
      level(ml)=REAL(stash_levels(ml+1,no))
    END DO
  END IF
ELSE IF (stash_control(st_output_bottom) == st_special_code) THEN
  ! Special level.
  ! The LEVEL array is not used by the model except to construct pp
  ! header items at output. The value of -1.0 is set as a flag for
  ! special levels so that routine PP_HEAD will insert the lbvc
  ! item in STASHmaster record.
  DO ml=1,num_levs_out
    level(ml)=-1.0
  END DO
ELSE
  DO ml=1,num_levs_out
    level(ml)=REAL(stash_control(st_output_bottom)+ml-1)
  END DO
END IF
!
!
! Now reset the number of output levels to 1 if vertical compression is
! to be done in SPATIAL.  NB: index_lev and level_list need to be filled
! with values corresponding to the full range of levels processed.
!
what_proc=stash_control(st_proc_no_code)
what_mean=(stash_control(st_gridpoint_code)/block_size)*block_size
IF (what_mean == vert_mean_base .OR. what_mean == global_mean_base &
   .OR. what_proc == st_time_series_code                          &
   .OR. what_proc == st_time_series_mean                          &
   .OR. what_proc == st_append_traj_code) num_levs_out=1
!
! Next compute the index for pseudo levels, if there are any
!
IF (stash_control(st_pseudo_in) >  0) THEN ! Input PSEUDO_LEVELS
  ni=stash_control(st_pseudo_in)
  num_pseudo_in=stash_pseudo_levels(1,ni)
  IF (stash_control(st_pseudo_out) >  0) THEN ! Output PSEUDO_LEVS
    no=stash_control(st_pseudo_out)
    num_pseudo_out=stash_pseudo_levels(1,no)
    indx1=0
    DO ml=1,num_pseudo_out
      ilev=stash_pseudo_levels(ml+1,no)   !  Level required
      DO kl=1,num_pseudo_in
        IF (stash_pseudo_levels(kl+1,ni) == ilev) THEN
          indx1=indx1+1
          index_pseudo_lev(indx1)=kl
          pseudo_level(indx1)=ilev
          GO TO 500
        END IF
      END DO
      icode=nonsense
      WRITE(cmessage,'(A,I6,A)') 'STLEVELS : Output pseudo level ', &
                           ilev, ' not found in input levels list'
      GO TO 999
      500        CONTINUE
    END DO
  ELSE  ! Illegal combination
    icode=nonsense
    WRITE(cmessage,'(A,I6,A)') 'STLEVELS : Input pseudo level list ', &
                        ni, ' has illegal output pseudo levels list'
    GO TO 999
  END IF
ELSE  ! Only levels lists are supported for pseudo levels
  num_pseudo_out=0
END IF
!
! Next expand the separate indexes and physical levels arrays into
! combined arrays if necessary, taking care not to overwrite earlier
! parts of the arrays.  If no pseudo-levels, set pseudo-level to 0.
!
IF (num_pseudo_out >  0) THEN
  DO k2=num_pseudo_out,1,-1
    DO ml=1,num_levs_out
      index_lev(ml+(k2-1)*num_levs_out)=                          &
        (index_pseudo_lev(k2)-1)*num_levs_in+index_lev(ml)
      level(ml+(k2-1)*num_levs_out)=level(ml)
    END DO
    DO ml=num_levs_out,1,-1
      pseudo_level(ml+(k2-1)*num_levs_out)=pseudo_level(k2)
    END DO
  END DO
  num_levs_out=num_levs_out*num_pseudo_out
ELSE
  DO ml=1,num_levs_out
    pseudo_level(ml)=0
  END DO
END IF
!
999   CONTINUE ! jump here for error return
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stlevels

