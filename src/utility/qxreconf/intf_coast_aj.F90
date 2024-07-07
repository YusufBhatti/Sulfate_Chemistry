! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: INTF_COAST_AJ------------------------------------------
!
!    Purpose: To calculate a suitable value of NSEARCH for call to
!             SPIRAL_S and then call SPIRAL_S
!
!    Reviewer:  Date of review:
!
!    Tested under compiler: cft77
!    Tested under OS version: UNICOS 7
!
!    Code version no: 1       Date: 17 November 1993
!
!   Programming standard: UM Doc Paper 3, version
!
!   Logucal component number:
!
!   Project task: S1
!
!
!   Documentation:
!     UMDP S1
!
!   -------------------------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration
SUBROUTINE intf_coast_aj                                          &
           (land_sea_mask,index_unres,no_point_unres,             &
           points_phi,points_lambda,data_field,sea_land,          &
           cyclic,maxdim)

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
,no_point_unres                                                   &
                  !INOUT number of unresolved points
,land_sea_mask(points_lambda*points_phi)                          &
                  !IN land sea mask
,index_unres(points_lambda*points_phi)                            &
                  !INOUT index to unresolved pts
,sea_land         !IN =0 for sea field  =1/-1 for land field

REAL ::                                                           &
 data_field(points_lambda*points_phi) !IN field

LOGICAL ::                                                        &
 cyclic           ! IN =T if data covers complete latitude circle

!   LOCAL VARIABLES

INTEGER ::                                                        &
 i,j,jj,jjj,k,jk,nsch                                             &
                      ! indices
,ipt,irow,icol                                                    &
                ! coordinate of unresolved point
,ipoint,iunres                                                    &
                ! do loop variable
,npoints                                                          &
                ! number of points in serach box
,maxdim                                                           &
               ! largest dimension of field
,ir((1+2*maxdim)*(1+2*maxdim))                                    &
                   ! row numbers of points to serach
,ic((1+2*maxdim)*(1+2*maxdim))                                    &
                   ! col numbers of points to search
,ind_search((1+2*maxdim)*(1+2*maxdim))                            &
                   ! index to points to search
               ! still unresolved after calling this subroutine
,isum_mask                                                        &
               ! number of surrounding points which have data
,isearch                                                          &
               ! largest dimension of field
,nsearch                                                          &
               ! minimum search radius required.
,land_sea_temp(points_lambda*points_phi) ! local copy of mask

CHARACTER (LEN=*), PARAMETER  :: RoutineName='INTF_COAST_AJ'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
nsearch=0

! Take local copy of land sea mask to work with
DO ipoint=1,points_lambda*points_phi
  land_sea_temp(ipoint)=land_sea_mask(ipoint)
END DO

! toggle land sea mask to exclude unresolved points from meaning process
DO iunres=1,no_point_unres
  IF (sea_land == 0) THEN
    land_sea_temp(index_unres(iunres))=1
  ELSE
    land_sea_temp(index_unres(iunres))=0
  END IF
END DO


! Loop around unresolved points
DO iunres=1,no_point_unres

  ! find unresolved point coordinate in terms of rows and cols
  ipt=index_unres(iunres)
  irow=INT(REAL(ipt)/REAL(points_lambda)+1)
  icol=ipt-(irow-1)*points_lambda

  isearch=maxdim
  DO i=1,maxdim

    ! Calculate number of points in search box
    npoints=(1+2*i)**2 ! number of grid points in search box

    ! calculate surrounding points' coords in terms of rows and cols
    jjj=1
    DO j=-i,i
      DO jj=jjj,jjj+2*i
        ir(jj)=irow+j
      END DO
      jjj=jjj+1+2*i
    END DO

    jjj=1+2*i
    jk=1
    DO j=-i,i
      DO jj=0,2*i
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

    DO ipoint=1,npoints
      IF (IABS(land_sea_temp(ind_search(ipoint))) == IABS(sea_land) &
        .AND. data_field(ind_search(ipoint)) >= 0.0) THEN
        isum_mask=isum_mask+1
      END IF
    END DO

    IF (isum_mask >  0) THEN
      isearch=MIN0(isearch,i)
      GO TO 100
    END IF

  END DO

  100    CONTINUE
  nsearch=MAX0(isearch,nsearch)

  ! If the search radius reaches the maximum, we have no need to
  ! do testing for the rest of the unresolved points
  IF (nsearch >= maxdim) EXIT
END DO


! DEPENDS ON: spiral_s
CALL spiral_s(land_sea_temp,index_unres,no_point_unres,           &
                 points_phi,points_lambda,data_field,nsearch,     &
                 sea_land,cyclic)


IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE intf_coast_aj
