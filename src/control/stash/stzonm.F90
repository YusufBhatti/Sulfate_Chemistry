! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: STZONM ---------------------------------------------------
!
!    Purpose: Calculate weighted zonal mean within a region specified
!             by a lower left hand and upper right hand corner.
!             (STASH service routine).
!
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    Project task: D7
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: STASH

SUBROUTINE stzonm(fieldin,vx,vy,fld_type,gr,halo_type,            &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  global_xstart,global_ystart,                    &
                  global_xend,global_yend,                        &
                  fieldout,                                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
                  icode,cmessage)

USE global_2d_sums_mod, ONLY: global_2d_sums
USE mpp_conf_mod, ONLY: exclude_halos_ew, exclude_halos_ns
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE sterr_mod, ONLY: st_no_data

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)  :: vx                ! input field size
INTEGER, INTENT(IN)  :: vy                ! input field size
INTEGER, INTENT(IN)  :: fld_type          ! field type (u/v/p)
INTEGER, INTENT(IN)  :: gr                ! input fld grid
INTEGER, INTENT(IN)  :: halo_type         ! halo type
INTEGER, INTENT(IN)  :: xstart            ! lower LH corner
INTEGER, INTENT(IN)  :: ystart            ! lower LH corner
INTEGER, INTENT(IN)  :: xend              ! upper RH corner
INTEGER, INTENT(IN)  :: yend              ! upper RH corner
INTEGER, INTENT(IN)  :: global_xstart     ! global version of xstart
INTEGER, INTENT(IN)  :: global_ystart     ! global version of ystart
INTEGER, INTENT(IN)  :: global_xend       ! global version of xend
INTEGER, INTENT(IN)  :: global_yend       ! global version of yend
INTEGER, INTENT(IN)  :: level_code        ! input level code
INTEGER, INTENT(IN)  :: mask_code         ! masking code
INTEGER, INTENT(IN)  :: weight_code       ! weighting code
REAL,    INTENT(IN)  :: rmdi              ! missing data indicator
INTEGER, INTENT(INOUT) :: icode             ! error return code

LOGICAL, INTENT(IN)  :: lwrap             ! TRUE if wraparound
LOGICAL, INTENT(IN)  :: lmasswt           ! TRUE if masswts OK
LOGICAL, INTENT(IN)  :: mask(vx+1,vy)     ! mask array

CHARACTER (LEN=errormessagelength), INTENT(INOUT) :: cmessage ! error return msg

REAL, INTENT(IN)     :: fieldin(vx,vy)        ! input field
REAL, INTENT(OUT)    :: fieldout(ystart:yend) ! output field
REAL, INTENT(IN)     :: pstar_weight(vx+1,vy) ! mass weight factor
REAL, INTENT(IN)     :: area_weight(vx+1,vy)  ! area weight factor

! ----------------------------------------------------------------------
! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking
!-----------------------------------------------------------------------

! Local variables
INTEGER :: i    ! Looper
INTEGER :: ii   ! Looper
INTEGER :: j    ! Looper
INTEGER :: jj   ! Looper

! sizes of subarea on this processor, not including halo areas
INTEGER :: local_sum_xstart
INTEGER :: local_sum_xend
INTEGER :: local_sum_ystart
INTEGER :: local_sum_yend

! number of rows of zonal data to sum on this processor (yend-ystart+1)
INTEGER :: partial_sum_data_sizex
INTEGER :: partial_sum_data_sizey

! return code from GCOM routines
INTEGER :: info

! Sums and partial sums
REAL, ALLOCATABLE :: partial_sumz(:,:,:)
REAL, ALLOCATABLE :: sumz(:,:)

INTEGER, PARAMETER :: bot=1     ! Index for bottom sum
INTEGER, PARAMETER :: top=2     ! Index for top sum

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STZONM'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! Initialise fieldout array - so all PE's have valid data
! (Only PEs on left of subdomain get the zonal means)
fieldout(:) = 0

! 1.0 Find the bounds of the actual data required in the summation
!    (ie. excluding the halos, contained within
!    xstart,xend,ystart,yend.

! DEPENDS ON: global_to_local_subdomain
CALL global_to_local_subdomain(exclude_halos_ew,exclude_halos_ns, &
                               gr, halo_type, mype,               &
                               global_ystart, global_xend,        &
                               global_yend, global_xstart,        &
                               local_sum_ystart, local_sum_xend,  &
                               local_sum_yend, local_sum_xstart)

! 1.1 And the number of partial sums on this processor

IF (local_sum_xstart  >   local_sum_xend)                         &
  local_sum_xend=local_sum_xend+vx-2*halosize(1,halo_type)

partial_sum_data_sizex=local_sum_xend-local_sum_xstart+1
partial_sum_data_sizey=local_sum_yend-local_sum_ystart+1

! 2.0 Calculate the partial sums

! Only do calculations if some of the subdomain exists on this
! processor
IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

  ! 2.2 Add do the actual sum

  ALLOCATE(partial_sumz(local_sum_xstart:local_sum_xend,    &
                 local_sum_ystart:local_sum_yend,2))
  ALLOCATE(sumz(local_sum_ystart:local_sum_yend,2))

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( local_sum_ystart, local_sum_yend, local_sum_xstart,      &
!$OMP         local_sum_xend, sumz, partial_sumz,                      &
!$OMP         lwrap, lasize, halosize, fld_type, halo_type, blsize,    &
!$OMP         mask, pstar_weight, area_weight, fieldin)                &
!$OMP PRIVATE ( i, j, ii )
  DO j=local_sum_ystart,local_sum_yend
    sumz(j,bot) = 0.0
    sumz(j,top) = 0.0
    DO i=local_sum_xstart,local_sum_xend
      partial_sumz(i,j,bot)=0.0
      partial_sumz(i,j,top)=0.0
    END DO

    DO i=local_sum_xstart,local_sum_xend

      IF ( lwrap .AND.                                              &
            (i  >   (lasize(1,fld_type,halo_type)-                  &
                     halosize(1,halo_type)))) THEN
        ! miss halos on wrap around
        ii=i-blsize(1,fld_type)
        ! ii > local_sum_xstart
      ELSE
        ii=i
      END IF

      ! Only do the sum if this is not a halo point

      IF (mask(ii,j)) THEN

        partial_sumz(ii,j,bot) = pstar_weight(ii,j)                 &
                               * area_weight(ii,j)
        partial_sumz(ii,j,top) = fieldin(ii,j) * pstar_weight(ii,j) &
                               * area_weight(ii,j)

      END IF ! if this point is to be processed
    END DO ! i : loop over columns
  END DO ! j : loop over rows
!$OMP END PARALLEL DO
END IF ! if subdomain covers this processor

! 3.0 Sums up the partial sums along each row to make a full sum

! So a sum along the processor row if the subdomain covers any
! processor(s) along the row

IF ((local_sum_ystart  /=  st_no_data) .AND.                      &
    (local_sum_yend  /=  st_no_data)) THEN

  CALL global_2d_sums(partial_sumz, partial_sum_data_sizex,      &
                      1, 0, 0, 2*partial_sum_data_sizey, sumz,   &
                      gc_proc_row_group)

  ! 3.1 And put the mean zonal values into the fieldout array when in
  !     the subdomain
  IF ((local_sum_xstart  /=  st_no_data) .AND.                    &
      (local_sum_xend  /=  st_no_data)) THEN

    DO j=local_sum_ystart,local_sum_yend
      jj=j+1-local_sum_ystart
      IF (sumz(j,bot)  ==  0.0) THEN
        fieldout(jj)=rmdi
      ELSE
        fieldout(jj)=sumz(j,top)/sumz(j,bot)
      END IF
    END DO

  END IF ! is this processor in the subdomain

  DEALLOCATE(partial_sumz)
  DEALLOCATE(sumz)

END IF ! does the subdomain intersect with this processor row

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE stzonm
