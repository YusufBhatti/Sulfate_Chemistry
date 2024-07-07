! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: STMERM ---------------------------------------------------
!
!    Purpose: Calculate weighted meridional mean within a region
!             specified by lower left hand and upper right hand corner.
!             (STASH service routine).
!
!
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!
!    Project task: D7
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------


!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE stmerm(fieldin,vx,vy,fld_type,gr,halo_type,            &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  global_xstart,global_ystart,                    &
                  global_xend,global_yend,                        &
                  fieldout,                                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
                  icode,cmessage)


USE um_parvars
USE um_parcore, ONLY: mype
USE global_2d_sums_mod, ONLY: global_2d_sums
USE mpp_conf_mod, ONLY: exclude_halos_ew, exclude_halos_ns
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE sterr_mod, ONLY: st_no_data

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) ::  vx                ! input field x size
INTEGER, INTENT(IN) ::  vy                ! input field y size
INTEGER, INTENT(IN) ::  fld_type          ! field type (u/v/p)
INTEGER, INTENT(IN) ::  gr                ! input fld grid
INTEGER, INTENT(IN) ::  halo_type         ! halo type
INTEGER, INTENT(IN) ::  xstart            ! lower LH corner
INTEGER, INTENT(IN) ::  ystart            ! lower LH corner
INTEGER, INTENT(IN) ::  xend              ! upper RH corener
INTEGER, INTENT(IN) ::  yend              ! upper RH corner
INTEGER, INTENT(IN) ::  global_xstart     ! global version
INTEGER, INTENT(IN) ::  global_ystart     ! global version
INTEGER, INTENT(IN) ::  global_xend       ! global version
INTEGER, INTENT(IN) ::  global_yend       ! global version
INTEGER, INTENT(IN) ::  level_code        ! input level code
INTEGER, INTENT(IN) ::  mask_code         ! masking code
INTEGER, INTENT(IN) ::  weight_code       ! weighting code
INTEGER, INTENT(INOUT) :: icode           ! error return code

CHARACTER(LEN=errormessagelength),INTENT(INOUT) :: cmessage  ! error return msg

LOGICAL, INTENT(IN) :: lwrap              ! TRUE if wraparound
LOGICAL, INTENT(IN) :: lmasswt            ! TRUE if masswts OK
LOGICAL, INTENT(IN) :: mask(vx+1,vy)      ! mask array

REAL, INTENT(IN) ::  rmdi                   ! missing data indicator
REAL, INTENT(IN) ::  fieldin(vx,vy)         ! input field
REAL, INTENT(IN) ::  pstar_weight(vx+1,vy)  ! mass weight factor
REAL, INTENT(IN) ::  area_weight(vx+1,vy)   ! area weight factor

REAL, INTENT(OUT) :: fieldout(xstart:xend)  ! output field

! ----------------------------------------------------------------------

! Local variables

INTEGER :: i              ! Looper
INTEGER :: ii             ! Looper
INTEGER :: j              ! Looper



! size of subarea on this processor, not including halo areas
INTEGER :: local_sum_xstart
INTEGER :: local_sum_xend
INTEGER :: local_sum_ystart
INTEGER :: local_sum_yend

! number of columns of meridional data to sum on this processor
! (xend-xstart+1)
INTEGER :: partial_sum_data_sizex
INTEGER :: partial_sum_data_sizey

! Sums and partial sums
REAL, ALLOCATABLE :: partial_sumtop_2d(:,:)
REAL, ALLOCATABLE :: partial_sumbot_2d(:,:)
REAL, ALLOCATABLE :: sumtop(:)
REAL, ALLOCATABLE :: sumbot(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STMERM'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking


! 1.0 Find the bounds of the actual data required in the summation
!    (ie. excluding the halos, contained within
!    xstart,xend,ystart,yend.

! DEPENDS ON: global_to_local_subdomain
CALL global_to_local_subdomain(exclude_halos_ew, exclude_halos_ns,&
  gr,halo_type,mype,                                              &
  global_ystart,global_xend,                                      &
  global_yend,global_xstart,                                      &
  local_sum_ystart,local_sum_xend,                                &
  local_sum_yend,local_sum_xstart)

IF (local_sum_xstart  >   local_sum_xend)                         &
  local_sum_xend=local_sum_xend+vx-                               &
                 2*halosize(1,halo_type)

! 1.1 And the number of partial sums on this processor
!     (If there is a wrap around, with the subarea starting
!      and ending on this processor, the halo cols are included
!      in this number - although no sums will actually be
!      carried out there)

partial_sum_data_sizex=                                           &
  local_sum_xend-local_sum_xstart+1
partial_sum_data_sizey=                                           &
  local_sum_yend-local_sum_ystart+1


! 1.2 Initialise the sum arrays
!     Note that partial arrays are transposed to make global sums
!     possible and straightforward later...

ALLOCATE (partial_sumtop_2d(local_sum_ystart:local_sum_yend, &
                      local_sum_xstart:local_sum_xend))
ALLOCATE (partial_sumbot_2d(local_sum_ystart:local_sum_yend, &
                      local_sum_xstart:local_sum_xend))
ALLOCATE (sumtop(local_sum_xstart:local_sum_xend))
ALLOCATE (sumbot(local_sum_xstart:local_sum_xend))
partial_sumtop_2d(:,:) = 0.0
partial_sumbot_2d(:,:) = 0.0
sumtop(:) = 0.0
sumbot(:) = 0.0

! 2.0 Calculate the partial sums

! Only do the calculations if some of the subdomain exists on this
! processor

IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

  ! 2.2 Do the actual sum

  DO i=local_sum_xstart,local_sum_xend

    IF ( lwrap .AND.                                            &
        (i  >   (lasize(1,fld_type,halo_type)-                  &
                 halosize(1,halo_type)))) THEN
      ! miss halos on wrap around
      ii=i-blsize(1,fld_type)
    ELSE
      ii=i
    END IF

    DO j=local_sum_ystart,local_sum_yend
      IF (mask(ii,j)) THEN

        partial_sumbot_2d(j,i)=                  &
          pstar_weight(ii,j)*area_weight(ii,j)
        partial_sumtop_2d(j,i)=                  &
          fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
      END IF ! if this point is to be processed
    END DO ! j : loop over rows
  END DO ! i : loop over columns
END IF ! if subdomain covers this processor

! 3.0 Sums up the partial sums down each column to make a full sum

! So a sum down the processor column if the subdomain covers any
! processor(s) along the column

IF ((local_sum_xstart  /=  st_no_data) .AND.                      &
    (local_sum_xend  /=  st_no_data)) THEN

  CALL global_2d_sums(partial_sumbot_2d, partial_sum_data_sizey,  &
                      1, 0, 0, partial_sum_data_sizex,            &
                      sumbot, gc_proc_col_group)
  CALL global_2d_sums(partial_sumtop_2d, partial_sum_data_sizey,  &
                      1, 0, 0, partial_sum_data_sizex,            &
                      sumtop, gc_proc_col_group)

  ! So now the partial_* arrays actually contain the full sums
  ! along the column

  ! 3.1 And put the mean meridional values into the fieldout array

  ! Only processors in the subdomain area need to record the
  ! results
  IF ((local_sum_ystart  /=  st_no_data) .AND.                    &
      (local_sum_yend  /=  st_no_data)) THEN

    DO i=local_sum_xstart,local_sum_xend
      ii = (local_sum_xstart - xstart) + i
      IF (sumbot(i)  ==  0.0) THEN
        fieldout(ii)=rmdi
      ELSE
        fieldout(ii)=sumtop(i)/sumbot(i)
      END IF
    END DO

  END IF ! is this processor in the subdomain

END IF ! does the subdomain intersect with this processor column

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stmerm
