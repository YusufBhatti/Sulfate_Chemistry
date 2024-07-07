! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   SUBROUTINE SPIRAL_S-----------------------------------------------
!
!
!   Programming standard:
!      Unified Model Documentation Paper No 3
!
!   System component: S121
!
!   System task: S1
!
!   Purpose:
!     Attempts to set a value at points which are unresolved when
!     interpolating between one grid and another.  A value is set
!     by finding the mean of surrounding points which do have data
!     set within a search radius determined by NSEARCH.
!
!   Documentation:
!     UMDP S1
!
!   -------------------------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration
SUBROUTINE spiral_s(land_sea_mask,index_unres,no_point_unres,     &
           points_phi,points_lambda,data_field,nsearch,sea_land,  &
           cyclic)

USE missing_data_mod, ONLY: rmdi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

!   ARGUMENTS:---------------------------------------------------

INTEGER ::                                                        &
 points_phi                                                       &
                  !IN number of rows in grid
,points_lambda                                                    &
                  !IN number of columns in grid
,nsearch                                                          &
                  !IN number of points in each direction to search
,no_point_unres                                                   &
                  !INOUT number of unresolved points
,land_sea_mask(points_lambda*points_phi)                          &
                  !IN land sea mask
,index_unres(points_lambda*points_phi)                            &
                  !INOUT index to unresolved pts
,sea_land         !IN =0 for sea field  =1/-1 for land field

REAL ::                                                           &
 data_field(points_lambda*points_phi) !INOUT field

LOGICAL ::                                                        &
 cyclic           ! IN =T if data covers complete latitude circle

!   Parameters

!   LOCAL VARIABLES

INTEGER ::                                                        &
 i,j,jj,jjj,k,jk                                                  &
                  ! indices
,ipt,irow,icol                                                    &
                ! coordinate of unresolved point
,ipoint,iunres                                                    &
                ! do loop variable
,npoints                                                          &
                ! number of points in serach box
,ir((1+2*nsearch)*(1+2*nsearch))                                  &
                   ! row numbers of points to serach
,ic((1+2*nsearch)*(1+2*nsearch))                                  &
                   ! col numbers of points to search
,ind_search((1+2*nsearch)*(1+2*nsearch))                          &
                   ! index to points to search
!
      ,not_yet_set                                                      &
                                    ! number of points still to set
      ,ind_yet_set(points_lambda*points_phi)                            &
                                             ! index of points
                     ! still unresolved after calling this subroutine
      ,isum_mask     ! number of surrounding points which have data

REAL ::                                                           &
 sum_data                                                         &
               ! sum of data surrounding unresolved points
, rmdi_tol ! values within this tolerance counted as missing

CHARACTER (LEN=*), PARAMETER :: RoutineName='SPIRAL_S'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!     EXTERNAL ROUTINES
!     None
!---------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

rmdi_tol  = ABS (rmdi) * 0.0001
!

! Calculate number of points in search box
npoints=(1+2*nsearch)**2 ! number of grid points in search box

! Loop around unresolved points
not_yet_set=0
DO iunres=1,no_point_unres

  ! find unresolved point coordinate in terms of rows and cols
  ipt=index_unres(iunres)
  irow= (ipt - 1)/points_lambda   +   1
  icol=ipt-(irow-1)*points_lambda

  ! calculate surrounding points' coords in terms of rows and cols
  jjj=1
  DO j=-nsearch,nsearch
    DO jj=jjj,jjj+2*nsearch
      ir(jj)=irow+j
    END DO
    jjj=jjj+1+2*nsearch
  END DO

  jjj=1+2*nsearch
  jk=1
  DO j=-nsearch,nsearch
    DO jj=0,2*nsearch
      ic(jk+jj*jjj)=icol+j
    END DO
    jk=jk+1
  END DO

  ! Check that col and rows are in range of grid
  DO ipoint=1,npoints
    IF (ic(ipoint) >  points_lambda) THEN
      IF (cyclic) THEN
        ic(ipoint)=ic(ipoint)-points_lambda
      ELSE
        ic(ipoint)=points_lambda
      END IF
    END IF
    IF (ic(ipoint) <  1) THEN
      IF (cyclic) THEN
        ic(ipoint)=ic(ipoint)+points_lambda
      ELSE
        ic(ipoint)=1
      END IF
    END IF
    IF (ir(ipoint) <  1) ir(ipoint)=1
    IF (ir(ipoint) >  points_phi) ir(ipoint)=points_phi
  END DO

  ! Form index search array
  DO ipoint=1,npoints
    ind_search(ipoint)=(ir(ipoint)-1)*points_lambda+ic(ipoint)
  END DO

  ! search for data around this point. If no data is found the point
  ! remains unresolved

  isum_mask=0   ! number of points with data found
  sum_data=0.0  ! sum of data of surrounding grid points

  DO ipoint=1,npoints
    IF (IABS(land_sea_mask(ind_search(ipoint))) == IABS(sea_land)  &
    .AND. data_field(ind_search(ipoint))  >   rmdi+rmdi_tol) THEN
      sum_data=sum_data+data_field(ind_search(ipoint))
      isum_mask=isum_mask+1
    END IF
  END DO

  IF (isum_mask  >   0) THEN
    ! data found - take mean
    data_field(ipt)=sum_data/REAL(isum_mask)
  ELSE
    ! data not found - point remains unresolved
    not_yet_set=not_yet_set+1
    ind_yet_set(not_yet_set)=ipt
  END IF

END DO


! amend output array with points remaining unresolved
IF (not_yet_set >  0) THEN
  DO ipoint=1,not_yet_set
    index_unres(ipoint)=ind_yet_set(ipoint)
  END DO
  no_point_unres=not_yet_set
ELSE
  no_point_unres=0
END IF
IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE spiral_s
