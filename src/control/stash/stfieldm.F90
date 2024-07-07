! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: STFIELDM -------------------------------------------------
!
!    Purpose: Calculate weighted field mean within a region specified
!             by a lower left hand and upper right hand corner.
!             Single level fields only.
!             (STASH service routine).
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE stfieldm(fieldin,vx,vy,fld_type,gr,halo_type,          &
                    lwrap,lmasswt,                                &
                    xstart,ystart,xend,yend,                      &
                    global_xstart,global_ystart,                  &
                    global_xend,global_yend,                      &
                    fieldout,                                     &
                    pstar_weight,                                 &
                    area_weight,mask,                             &
                    level_code,mask_code,weight_code,rmdi,        &
                    icode,cmessage)


USE UM_parvars
USE um_parcore, ONLY: mype
USE global_2d_sums_mod, ONLY: global_2d_sums
USE mpp_conf_mod, ONLY: exclude_halos_ew, exclude_halos_ns
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE sterr_mod, ONLY: st_no_data
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: vx                ! input x size
INTEGER, INTENT(IN) :: vy                ! input y size
INTEGER, INTENT(IN) :: fld_type          ! field type (u/v/p)
INTEGER, INTENT(IN) :: gr                ! input fld grid
INTEGER, INTENT(IN) :: halo_type         ! halo type
INTEGER, INTENT(IN) :: xstart            ! lower LH corner
INTEGER, INTENT(IN) :: ystart            ! lower LH corner
INTEGER, INTENT(IN) :: xend              ! upper RH corner
INTEGER, INTENT(IN) :: yend              ! upper RH corner
INTEGER, INTENT(IN) :: global_xstart     ! global version of xstart
INTEGER, INTENT(IN) :: global_ystart     ! global version of ystart
INTEGER, INTENT(IN) :: global_xend       ! global version of xend
INTEGER, INTENT(IN) :: global_yend       ! global version of yend
INTEGER, INTENT(IN) :: level_code        ! input level code
INTEGER, INTENT(IN) :: mask_code         ! masking code
INTEGER, INTENT(IN) :: weight_code       ! weighting code

INTEGER, INTENT(OUT) :: icode            ! error return code
CHARACTER(LEN=errormessagelength)     :: cmessage         ! OUT error return msg

LOGICAL, INTENT(IN) :: lwrap             ! TRUE if wraparound
LOGICAL, INTENT(IN) :: lmasswt           ! TRUE if masswts OK
LOGICAL, INTENT(IN) :: mask(vx+1,vy)     ! mask array

REAL, INTENT(IN) :: rmdi                 ! missing data indicator
REAL, INTENT(IN) :: fieldin(vx,vy)       ! input field
REAL, INTENT(IN) :: pstar_weight(vx+1,vy) ! mass weight factor
REAL, INTENT(IN) :: area_weight(vx+1,vy) ! area weight factor
REAL, INTENT(OUT) :: fieldout            ! output scalar

!----------------------------------------------------------------------
! Local variables

INTEGER :: i                             ! Looper
INTEGER :: ii                            ! Looper
INTEGER :: j                             ! Looper

! limits of local data to be summed
INTEGER :: local_sum_xstart
INTEGER :: local_sum_xend
INTEGER :: local_sum_ystart
INTEGER :: local_sum_yend

REAL :: my_sum(2)                       ! result of sums
REAL :: partial(vx+1,vy,2)              ! partial sums

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STFIELDM'

! ----------------------------------------------------------------------
!  0. Initialise sums
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

icode          = 0
my_sum(:)      = 0.0
partial(:,:,:) = 0.0

! ----------------------------------------------------------------------

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
  local_sum_xend=local_sum_xend+vx-2*offx

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

        partial(ii,j,1) =                                       &
          pstar_weight(ii,j)*area_weight(ii,j)
        partial(ii,j,2) =                                       &
          fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)

      END IF ! if this point is to be processed
    END DO ! j : loop over rows
  END DO ! i : loop over columns
END IF ! if subdomain covers this processor

! 3.0  add all the partial sums together, and store
CALL global_2d_sums(partial, vx+1, vy, 0, 0, 2, my_sum,           &
                    gc_all_proc_group)

IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

  IF (my_sum(1) ==  0.0) THEN
    fieldout=rmdi
  ELSE
    fieldout=my_sum(2)/my_sum(1)
  END IF
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stfieldm
