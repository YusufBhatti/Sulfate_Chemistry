! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: STGLOM ---------------------------------------------------
!
!    Purpose: Calculate weighted global mean within a region specified
!             by a lower left hand and upper right hand corner.
!             Multiple level fields.
!             (STASH service routine).
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------

!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH
SUBROUTINE stglom(fieldin,vx,vy,vz,fld_type,gr,halo_type,         &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  global_xstart,global_ystart,                    &
                  global_xend,global_yend,                        &
                  fieldout,index_lev,zsize,                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
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

INTEGER, INTENT(IN) :: vx                 ! input x size
INTEGER, INTENT(IN) :: vy                 ! input y size
INTEGER, INTENT(IN) :: vz                 ! input z size
INTEGER, INTENT(IN) :: fld_type           ! field type (u/v/p)
INTEGER, INTENT(IN) :: gr                 ! input fld grid
INTEGER, INTENT(IN) :: halo_type          ! halo type
INTEGER, INTENT(IN) :: xstart             ! lower LH corner
INTEGER, INTENT(IN) :: ystart             ! lower LH corner
INTEGER, INTENT(IN) :: xend               ! upper RH corner
INTEGER, INTENT(IN) :: yend               ! upper RH corner
INTEGER, INTENT(IN) :: global_xstart      ! global version
INTEGER, INTENT(IN) :: global_ystart      ! global version
INTEGER, INTENT(IN) :: global_xend        ! global version
INTEGER, INTENT(IN) :: global_yend        ! global version
INTEGER, INTENT(IN) :: zsize              ! no of horiz levels to process
INTEGER, INTENT(IN) :: index_lev(zsize)   ! offset for each horiz lev
INTEGER, INTENT(IN) :: level_code         ! input level code
INTEGER, INTENT(IN) :: mask_code          ! masking code
INTEGER, INTENT(IN) :: weight_code        ! weighting code
INTEGER, INTENT(INOUT) :: icode           ! error return code

CHARACTER (LEN=errormessagelength), INTENT(INOUT) :: cmessage ! error return msg

LOGICAL, INTENT(IN) :: lwrap              ! TRUE if wraparound
LOGICAL, INTENT(IN) :: lmasswt            ! TRUE if masswts OK
LOGICAL, INTENT(IN) :: mask(vx+1,vy)      ! mask array
LOGICAL :: ts_column_flag                 ! time series column output

REAL, INTENT(IN) :: rmdi                  ! missing data indic
REAL, INTENT(IN) :: fieldin(vx,vy,vz)     ! input field
REAL, INTENT(IN) :: pstar_weight(vx+1,vy,zsize) ! mass weight factor
REAL, INTENT(IN) :: area_weight(vx+1,vy)  ! area weight factor
REAL, INTENT(OUT) :: fieldout(zsize)      ! output field
! ----------------------------------------------------------------------
! Local variables

! Loopers
INTEGER :: i
INTEGER :: j
INTEGER :: k
INTEGER :: ii
INTEGER :: jj
INTEGER :: kk
INTEGER :: kkk

! limits of local data to be summed
INTEGER :: local_sum_xstart
INTEGER :: local_sum_xend
INTEGER :: local_sum_ystart
INTEGER :: local_sum_yend

! Data sizes
INTEGER :: xsize
INTEGER :: ysize

! Sums/partial sums
REAL :: sumgtop                           ! top sum
REAL :: sumgbot                           ! bottom sum
REAL, ALLOCATABLE :: partial(:,:,:,:)     ! Array for partial sums
REAL :: level_sums(zsize,2)               ! level sums

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STGLOM'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! ----------------------------------------------------------------------
!  0. Initialise sums
sumgtop            = 0.0
sumgbot            = 0.0
ts_column_flag     = .FALSE.

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

! Set some sizes and allocate the partial array
xsize = local_sum_xend - local_sum_xstart + 1
ysize = local_sum_yend - local_sum_ystart + 1
ALLOCATE( partial(xsize, ysize, zsize, 2) )
partial(:,:,:,:) = 0.0

! 2.0 Calculate the partial sums

! Only do the calculations if some of the subdomain exists on this
! processor

IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

  ! 2.2 Do the actual sum
  DO kk=1,zsize
    k=index_lev(kk)


    IF (lmasswt) THEN  ! mass weighting?
      kkk=kk
    ELSE              ! no: therefore only 1 level of pstar_weight
      kkk=1          !     initialised (=1.0)
    END IF

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

          partial(i-local_sum_xstart+1,j-local_sum_ystart+1,kk,1) = &
                  pstar_weight(ii,j,kkk)*area_weight(ii,j)
          partial(i-local_sum_xstart+1,j-local_sum_ystart+1,kk,2) = &
                  fieldin(ii,j,k)*pstar_weight(ii,j,kkk)*area_weight(ii,j)
        END IF ! if this point is to be processed
      END DO ! j : loop over rows
    END DO ! i : loop over columns
  END DO ! kk : loop over levels
END IF ! if subdomain covers this processor


IF (xsize == 1 .AND. ysize == 1 .AND. zsize >= 1) THEN
  ! If only single point, assume timeseries column output and skip global sum
  ts_column_flag = .TRUE.
ELSE
  ! 3.0  add all the partial sums together, and store
  CALL global_2d_sums(partial, xsize, ysize, 0, 0, zsize*2, level_sums, &
                        gc_all_proc_group)

  ! Sum up bottom and tops across levels
  DO k = 1, zsize
    sumgbot = sumgbot + level_sums(k,1)
    sumgtop = sumgtop + level_sums(k,2)
  END DO
END IF

! Assign fieldout value(s)
IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

  IF (ts_column_flag) THEN ! ts column partial summs
    DO k = 1, zsize
      IF (partial(1,1,k,1) == 0.0) THEN
        fieldout(k) = rmdi
      ELSE
        fieldout(k) = partial(1,1,k,2) / partial(1,1,k,1)
      END IF
    END DO
  ELSE ! Global field mean
    IF (sumgbot  ==  0.0) THEN
      fieldout(1)=rmdi
    ELSE
      fieldout(1)=sumgtop/sumgbot
    END IF
  END IF
ELSE ! Put some value
  DO k = 1, zsize
    fieldout(k) = rmdi
  END DO
END IF  ! end of assign fieldout values 

DEALLOCATE( partial )

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stglom
