! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Idealised UM diagnostics routine

MODULE diagnostics_ideal_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Calculates Idealised UM diagnostics (held in STASH section 53).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   Idealised routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.

!   Diagnostics currently available (in numerical order):
!   Item  Description
!    181  Temperature increment
!    182  Water vapour increment
!    183  Reserved for Liquid Cloud Condensate increment if inc possible
!    184  Reserved for Frozen Cloud Condensate increment if inc possible
!    185  u wind increment
!    186  v wind increment
!
!    190  Potential temperature increment
!
!    201  Rate of change of column water vapour (kg/m2/s)
!    202  theta reference profile (K)
!    203  q reference profile (kg/kg)
!    204  u wind reference profile (m/s)
!    205  v wind reference profile (m/s)
!    206  Column cvT energy change due to forcing (Js/m2)
!    207  Column 0.5u*u energy change due to U forcing (Js/m2)
!    208  Column 0.5v*v energy change due to V forcing (Js/m2)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UMDP standards.
!------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_IDEAL_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE diagnostics_ideal(row_length, rows,    &
     stashwork                                    &
     )

USE idealised_diag_mod, ONLY:                                             &
  dt_inc_ideal_um, dq_inc_ideal_um, du_inc_ideal_um, dv_inc_ideal_um,     &
  dtheta_inc_ideal_um, dcolqdt_ideal_um, l_stored_ref, diag_theta_ref,    &
  diag_q_ref, diag_u_ref, diag_v_ref, de_cvt_ideal_um, de_u2_ideal_um,    &
  de_v2_ideal_um

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   udims, vdims, tdims, pdims

USE UM_ParVars, ONLY: &
  at_extremity          ! Indicates if this processor is at north,
                        ! south, east or west of the processor grid


USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf

USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr    ! allows printing

USE nlsizes_namelist_mod, ONLY: model_levels

! Subroutines
USE ereport_mod, ONLY: ereport

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::  &
  row_length            & ! model row_length
 ,rows                    ! number of rows for theta grid

REAL, INTENT(INOUT) ::                                            &
 stashwork(*)     ! STASH workspace

!------------------------------------------------------------------------------
! Local variables

INTEGER ::           &
  i, j, k            & ! loop counters
 ,u_rows             & ! number of rows for u grid
 ,v_rows             & ! number of rows for v grid
 ,icode              & ! error code
 ,sect               & ! Model section number
 ,im_index           & ! model index
 ,item                 ! diagnostic stash item number

REAL  ::                                 &
  work_nozero(tdims%i_start:tdims%i_end, &
       tdims%j_start:tdims%j_end,        &
       1:tdims%k_end)                    & 
 ,work(tdims%i_start:tdims%i_end,        &
       tdims%j_start:tdims%j_end,        &
       tdims%k_start:tdims%k_end)        & 
 ,uwork(udims%i_start:udims%i_end,       &
        udims%j_start:udims%j_end,       &
        udims%k_start:udims%k_end)       &  
 ,vwork(vdims%i_start:vdims%i_end,       &
        vdims%j_start:vdims%j_end,       &
        vdims%k_start:vdims%k_end)  


CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'DIAGNOSTICS_IDEAL'


! Variables required for Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Initialise variables

sect = 53         ! stash section for idealised UM

icode = 0         ! set error code to zero

im_index = 1

! Work out number of rows for U and V grids - will depend on whether
! ENDGame or not.

u_rows = udims%j_end - udims%j_start + 1
v_rows = vdims%j_end - vdims%j_start + 1

!-----------------------------------------------------------------------------
! Increments to prognostics
!-----------------------------------------------------------------------------

item = 181  ! temperature increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! Required for ENDGame as copydiag_3d not able to cope with zeroth level

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_nozero,dt_inc_ideal_um)
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_nozero(i,j,k) = dt_inc_ideal_um(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       work_nozero,                                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag_3d(item 181)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 182  ! humidity increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! Required for ENDGame as copydiag_3d not able to cope with zeroth level

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_nozero,dq_inc_ideal_um)
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_nozero(i,j,k) = dq_inc_ideal_um(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       work_nozero,                                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag_3d(item 182)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF

item = 185  ! u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       du_inc_ideal_um,                                          &
       row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag_3d(item 185)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 186  ! v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       dv_inc_ideal_um,                                          &
       row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag_3d(item 186)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)


item = 190  ! theta (potential temperature) increment
IF (icode <= 0 .AND. sf(item,sect)) THEN
  ! Required for ENDGame as copydiag_3d not able to cope with zeroth level

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work_nozero,dtheta_inc_ideal_um)
  DO k =             1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_nozero(i,j,k) = dtheta_inc_ideal_um(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       work_nozero,                                              &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag_3d(item 190)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!------------------------------------------------------------------------------
! rate of change of column water vapour
!------------------------------------------------------------------------------

item = 201
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                   &
       dcolqdt_ideal_um,row_length,rows,                             &
       0,0,0,0, at_extremity,atmos_im,sect,item,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 201)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!------------------------------------------------------------------------------
! Reference profiles on model levels
!------------------------------------------------------------------------------

item = 202  ! theta reference profile
IF (icode <= 0 .AND. sf(item,sect)) THEN

  IF (l_stored_ref) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work,diag_theta_ref)
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work(i,j,k) = diag_theta_ref(k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE         ! zero array output
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work)
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       work,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 202)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 203  ! q reference profile
IF (icode <= 0 .AND. sf(item,sect)) THEN
  IF (l_stored_ref) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work,diag_q_ref)
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work(i,j,k) = diag_q_ref(k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE         ! zero array output
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(tdims,work)
    DO k =             1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          work(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       work,                                                     &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 203)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 204  ! u reference profile
IF (icode <= 0 .AND. sf(item,sect)) THEN
  IF (l_stored_ref) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,uwork,diag_u_ref)
    DO k =             1, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          uwork(i,j,k) = diag_u_ref(k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE         ! zero array output
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(udims,uwork)
    DO k =             1, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          uwork(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       uwork,                                                    &
       row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 204)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

item = 205  ! v reference profile
IF (icode <= 0 .AND. sf(item,sect)) THEN
  IF (l_stored_ref) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,vwork,diag_v_ref)
    DO k =             1, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          vwork(i,j,k) = diag_v_ref(k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  ELSE         ! zero array output
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(vdims,vwork)
    DO k =             1, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          vwork(i,j,k) = 0.0
        END DO
     END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       vwork,                                                    &
       row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 205)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!------------------------------------------------------------------------------
! Change in energy due to dT forcing   integral of cvdT
!------------------------------------------------------------------------------

item = 206
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                   &
       de_cvt_ideal_um,row_length,rows,                              &
       0,0,0,0, at_extremity,atmos_im,sect,item,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 201)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)
!------------------------------------------------------------------------------
! Change in energy due to dU forcing   integral of udu+0.5du*du
!------------------------------------------------------------------------------

item = 207
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                   &
       de_u2_ideal_um,row_length,u_rows,                             &
       0,0,0,0, at_extremity,atmos_im,sect,item,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 207)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)
!------------------------------------------------------------------------------
! Change in energy due to dV forcing   integral of vdv+0.5dv*dv
!------------------------------------------------------------------------------

item = 208
IF (icode <= 0 .AND. sf(item,sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,sect,im_index)),                   &
       de_v2_ideal_um,row_length,v_rows,                             &
       0,0,0,0, at_extremity,atmos_im,sect,item,                     &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Idealised: error in copydiag(item 208)"
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_ideal

END MODULE diagnostics_ideal_mod
