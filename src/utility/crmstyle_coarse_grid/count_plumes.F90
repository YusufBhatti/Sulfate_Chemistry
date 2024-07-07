! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Number plumes

MODULE count_plumes_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COUNT_PLUMES_MOD'

CONTAINS

SUBROUTINE count_plumes(ncols,nrows,iprint,nn_mask,           &
                        plume_mask,                           &
                        nups, diam_pdf, size_pdf, fract_pdf,dxdy_pdf)


USE crmstyle_cntl_mod, ONLY:                                               &
   nbins_diam, nbins_size, nbins_fract, nbins_dxdy,                        &
   diam_bin, size_bin, fract_bin

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words

USE umPrintMgr                       ! Required to write output

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! ------------------------------------------------------------------------------
! Description:
! For a given mask of 1's and 0's identify separate plumes, numbering them
! and working out their sizes
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

INTEGER, INTENT(IN) ::  &
  ncols                 &  ! Number of columns in mask
 ,nrows                 &  ! Number of rows in mask
 ,iprint                &  ! printing switch
 ,nn_mask                  ! total number of buoyant points in mask

LOGICAL, INTENT(IN) ::     &
  plume_mask(ncols,nrows)    ! A mask identifying plumes

INTEGER, INTENT(OUT) ::    &
  nups                       ! number of separate plumes

! May want at some point
!INTEGER, INTENT(OUT) ::    &
!  plume_id(ncols,nrows)      ! Id numbers for plumes each true point in
!                             ! plume_mask is given a number indicating its
!                             ! plume identifier from 1 to nups

REAL(wp),INTENT(OUT) ::    &
  diam_pdf(nbins_diam)     & ! PDF of diameters
 ,size_pdf(nbins_size)     & ! PDF of sizes       (NOT USED)
 ,fract_pdf(nbins_fract)   & ! PDF of fractions   (NOT USED)
 ,dxdy_pdf(nbins_dxdy)       ! PDF of directions  (NOT USED)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! The total number of buoyant points can exceed a iwp number in some cases
! so it is possible for the plume number of exceed this.
!----------------------------------------------------------------------------
! dxdy (10)
!   <-1.5,-1.5 <-1.0, -1.0, -1.0<-0.5, -0.5=<0.0,
!   0.0<=0.5, 0.5<1.0, 1.0, 1.<=1.5, <1.5
!
!
!----------------------------------------------------------------------------
! Local arrays
INTEGER ::            &
  i, j, i2, j2                 ! loop counters
INTEGER ::    kk               ! loop counters all points - big numbers

INTEGER ::           &
 im1,ip1,jm1,jp1          & ! indicies for points surrounding a buoyant point
                            ! safe up to values of ~32768
 ,iproblem                  ! Used to indicate a problem when identifying plumes

INTEGER ::                &
  n1,n2,n3,n4,n5,n6,n7,n8   ! Used in identifying plumes


INTEGER   :: dx ,dy
INTEGER      ::     &
  nn                & ! label a plume
 ,imatch            & ! Used to indicate plumes touch
 ,iplume_match        ! Holds number of matching plume


INTEGER :: ntmp, nfuture, nf2
INTEGER :: idx, idy
INTEGER ::          &
  ncnt              & ! Number of bouyant points making up plumes
 ,nups2             & ! Revised number of plumes after check on touching
 ,max_plume           ! Maximum number of points making a plume

INTEGER    ::        &
  bcu2(ncols,nrows)      ! Each plume has a different number

INTEGER        ::             &
  coordx(nn_mask)             &
 ,coordy(nn_mask)             &
 ,coordp(nn_mask)


INTEGER, ALLOCATABLE ::                     &
  num_pl(:),coordx_pl(:,:),coordy_pl(:,:),  &
  plume_no(:),coordx_min(:),coordy_min(:),  &
  coordx_max(:),coordy_max(:),index_pl(:)

INTEGER                :: &
  plume_id(ncols,nrows)      ! Id numbers for plumes each true point in
                             ! plume_mask is given a number indicating its
                             ! plume identifier from 1 to nups

! Used in diameter of a plume calculation
INTEGER ::      &
 i_xmax         &  ! x max
,i_ymax         &  ! y max
,i_xmin         &  ! x min
,i_ymin         &  ! y min
,i_st           &  ! start coordinate for x
,dia1           &
,dia2

INTEGER   :: &
  sum_idx          ! x - x start   - get a sign for dx/dy

REAL :: dxdy, fract, diam, num_box

REAL ::          &
  old_to_p_sq    &
 ,old_to_p       &
 ,cenx           &
 ,ceny           &
 ,rad            &
 ,rad_sq         &
 ,old_to_new     &
 ,xspan          &
 ,yspan          &
 ,maxspan

REAL :: &
 xmax   &  ! x  max
,ymax   &  ! y max
,xmin   &  ! x  min
,ymin      ! y min

REAL(wp), ALLOCATABLE  :: &
  x(:,:)                  &  ! x  - coordinate real value
 ,y(:,:)                     ! y  - coordinate real value

CHARACTER(LEN=1) ::            &
  clet(0:503)                  &  ! Letters for printing plume maps
 ,scut(ncols,nrows)               ! Holds plume map to be printed

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'COUNT_PLUMES'

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------------
DATA  (clet(i),i=0,503) /'-','0','1','2','3','4','5','6','7','8','9', &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','@',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','%',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','!',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','*',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','?',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','#',                              &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        '0','1','2','3','4','5','6','7','8','9',  &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z','+',                              &
                        'A','B','C','D','E','F','G','H','I','J',  &
                        'K','L','M','N','O','P','Q','R','S','T',  &
                        'U','V','W','X','Y','Z','a','b','c','d',  &
                        'e','f','g','h','i','j','k','l','m','n',  &
                        'o','p','q','r','s','t','u','v','w','x',  &
                        'y','z'/

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
!----------------------------------------------------------------------------
! Initialise arrays
!----------------------------------------------------------------------------
bcu2(:,:) = 0

nups = 0
nups2 = 0
ncnt = 1

diam_pdf(:) = 0.0
dxdy_pdf(:) = 0.0
size_pdf(:) = 0.0
fract_pdf(:) = 0.0

!----------------------------------------------------------------------------
! First sweep through field identifying different plumes
! Results from this sweep often have plumes touching.
!----------------------------------------------------------------------------

DO j=1,nrows
  jm1 = j-1
  IF (jm1 == 0) jm1=1
  jp1 = j+1
  IF (jp1 > nrows) jp1=nrows
  DO i= 1,ncols
    im1 = i-1
    IF (im1 == 0)  im1=1
    ip1 = i+1
    IF (ip1 > ncols) ip1=ncols

    n1 =1
    n2 =1
    n3 =1
    n4 =1
    n5 =1
    n6 =1
    n7 =1
    n8 =1

    IF (plume_mask(i,j)) THEN

      coordx(ncnt) = i
      coordy(ncnt) = j
      ncnt = ncnt+ 1

      ! Are there any squares that have been counted before?
      ! Go through squares that qualified last time around, .AND. set ns to
      ! zero=> this 'thing' has already been counted

      !        i-1 i  i+1
      !  j-1   n3  n2  n4
      !  j     n1  X   n5
      !  j+1   n8  n7  n6

      IF (plume_mask(im1, j) .AND. i /= 1 )                   n1 = 0
      IF (plume_mask(i, jm1) .AND. j /= 1 )                   n2 = 0
      IF (plume_mask(im1, jm1) .AND. (i /= 0 .AND. j /= 0))   n3 = 0
      IF (plume_mask(ip1, jm1) .AND. j /= 1 .AND. i /= ncols) n4 = 0

      IF (plume_mask(ip1, j) .AND. i /= ncols )                    n5 = 0
      IF (plume_mask(ip1, jp1) .AND. i /= ncols .AND. j /= nrows ) n6 = 0
      IF (plume_mask(i  , jp1) .AND. j /= nrows)                   n7 = 0
      IF (plume_mask(im1, jp1) .AND. i /= 1 .AND. j /= nrows)      n8 = 0

      ! If there's a touched square before, then don't count as new
      ntmp = n1*n2*n3*n4

      ! Do future squares touch
      nfuture = n5*n6*n7*n8

      ! Does one of the future neighbours already have a name
      nf2=bcu2(ip1,j)+bcu2(ip1,jp1)+bcu2(i,jp1)+bcu2(im1,jp1)

      IF (ntmp == 1) THEN
        ! This is a new thing  does this plume already or
        ! a future neighbour have a name ?
        IF (bcu2(i,j) == 0 .AND. nf2 == 0) THEN
          nups = nups + 1
          bcu2(i,j) = nups      ! plume number
        ELSE
          ! A future neighbour has a name
          IF (n5 == 0 .AND. bcu2(ip1,  j) /= 0)  bcu2(i,j)= bcu2(ip1, j)
          IF (n6 == 0 .AND. bcu2(ip1,jp1) /= 0 ) bcu2(i,j)= bcu2(ip1,jp1)
          IF (n7 == 0 .AND. bcu2(i , jp1) /= 0 ) bcu2(i,j)= bcu2(i, jp1)
          IF (n8 == 0 .AND. bcu2(im1,jp1) /= 0 ) bcu2(i,j)= bcu2(im1,jp1)
        END IF  ! test on bcu2 & nf2
        IF (nfuture == 0) THEN
          ! set bcu2 for future points to same plume
          IF (n5 == 0) bcu2(ip1, j) = bcu2(i,j)
          IF (n6 == 0) bcu2(ip1, jp1) = bcu2(i,j)
          IF (n7 == 0) bcu2(i, jp1) = bcu2(i,j)
          IF (n8 == 0) bcu2(im1, jp1) = bcu2(i,j)
        END IF
      ELSE       ! Not a new plume  ntmp = 0

        !        i-1 i  i+1
        !  j-1   n3  n2  n4
        !  j     n1  X   n5
        !  j+1   n8  n7  n6

        ! Specail cases
        IF (j == 1 .AND. i == 1) THEN    ! First row .AND. col

          nups = nups + 1
          bcu2(i,j) = nups      ! plume number

          IF (nfuture == 0) THEN
            ! set bcu2 for future points to same plume
            IF (n5 == 0) bcu2(ip1, j)  = nups
            IF (n6 == 0) bcu2(ip1,jp1) = nups
            IF (n7 == 0) bcu2(i, jp1)  = nups
            ! Note n8  - no gridpoint
          END IF

        ELSE IF (i == 1 .AND. n1 == 0 .AND. n2*n4 == 1) THEN   ! first column

          nups = nups + 1
          bcu2(i,j) = nups      ! plume number  number 1 to n
          IF (nfuture == 0) THEN
            ! set bcu2 for future points to same plume
            IF (n5 == 0) bcu2(ip1, j)  = nups
            IF (n6 == 0) bcu2(ip1,jp1) = nups
            IF (n7 == 0) bcu2(i, jp1)  = nups
            ! Note n8  - no gridpoint
          END IF

        ELSE IF (j == 1 .AND. n2 == 0 .AND. n1 == 1 ) THEN   ! first row

          nups = nups + 1
          bcu2(i,j) = nups      ! plume number  number 1 to n
          IF (nfuture == 0) THEN
            ! set bcu2 for future points to same plume
            IF (n5 == 0) bcu2(ip1, j)  = nups
            IF (n6 == 0) bcu2(ip1,jp1) = nups
            IF (n7 == 0) bcu2(i, jp1)  = nups
            IF (n8 == 0) bcu2(im1, jp1) = nups
          END IF

        ELSE

          ! Point already has neighbours
          ! Do any of the neighbours have a plume number ? Should be true.

          IF (n1 == 0) bcu2(i,j)=bcu2(im1, j)
          IF (n2 == 0) bcu2(i,j)=bcu2(i, jm1)
          IF (n3 == 0) bcu2(i,j)=bcu2(im1, jm1)
          IF (n4 == 0) bcu2(i,j)=bcu2(ip1, jm1)

          IF (nfuture == 0) THEN

            ! set bcu2 for future points to same plume
            IF (n5 == 0) bcu2(ip1, j)  = bcu2(i,j)
            IF (n6 == 0) bcu2(ip1,jp1) = bcu2(i,j)
            IF (n7 == 0) bcu2(i, jp1)  = bcu2(i,j)
            IF (n8 == 0) bcu2(im1,jp1) = bcu2(i,j)

          END IF
        END IF

      END IF  ! ntmp

      ! store plume number of point

      coordp(ncnt-1)= bcu2(i,j)
      IF (coordp(ncnt-1) == 0 ) THEN   ! problem
        PRINT*, ' PROBLEM ',i,j,ntmp,nfuture,nf2, ' n ',n1,n2,n3,n4,n5,n6,n7,n8
      END IF

    END IF   ! test on mask
  END DO
END DO

! take one away to get correct number of buoyant points
ncnt = ncnt -1

!----------------------------------------------------------------------------
! Second sweep through data - required if more than 1 plume to check not
! touching
!----------------------------------------------------------------------------

IF (nups > 1 .AND. ncnt > 0) THEN

  ALLOCATE(num_pl(nups))
  ALLOCATE(index_pl(nups))
  ALLOCATE(coordx_pl(ncnt,nups))
  ALLOCATE(coordy_pl(ncnt,nups))

  index_pl(:)=1
  num_pl(:) = 0

  ! Resort data to make checking of touching plumes easier?
  iproblem = 0
  ! I think there are now no problem points left in the code but still has a
  ! check.
  DO i =1,ncnt
    IF ( coordp(i) == 0 .OR. coordp(i) > nups) THEN
      WRITE(umMessage,'(A,4I10)') ' PROBLEM ',i,coordp(i),coordx(i),coordy(i)
      CALL umPrint(umMessage)
      iproblem = 1
    ELSE
      num_pl(coordp(i)) = num_pl(coordp(i)) + 1
      coordx_pl(num_pl(coordp(i)),coordp(i)) = coordx(i)
      coordy_pl(num_pl(coordp(i)),coordp(i)) = coordy(i)
    END IF
  END DO

  ! Check through plumes starting with the last looking for touching plumes
  nups2= nups

  DO i = nups,2,-1

    imatch=0

    DO j=1,num_pl(i)    ! members of plume i

      DO i2 = i-1,1,-1   ! other plumes

        DO j2=1,num_pl(i2)
          ! Does plume i2 touch  plume i ?
          IF (((coordx_pl(j2,i2) == coordx_pl(j,i)+1 .AND.                &
              coordy_pl(j2,i2) == coordy_pl(j,i))     .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)+1 .AND.                 &
              coordy_pl(j2,i2) == coordy_pl(j,i)+1)   .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)+1 .AND.                 &
              coordy_pl(j2,i2) == coordy_pl(j,i)-1)   .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)   .AND.                 &
              coordy_pl(j2,i2) == coordy_pl(j,i)+1)   .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)   .AND.                 &
              coordy_pl(j2,i2) == coordy_pl(j,i)-1)   .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)-1  .AND.                &
              coordy_pl(j2,i2) == coordy_pl(j,i)+1)   .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)-1  .AND.                &
              coordy_pl(j2,i2) == coordy_pl(j,i))     .OR.                &
              (coordx_pl(j2,i2) == coordx_pl(j,i)-1  .AND.                &
              coordy_pl(j2,i2) == coordy_pl(j,i)-1) ) .AND. imatch == 0 ) THEN

            ! relabel all of plume
            nups2 = nups2 -1
            imatch=1
            iplume_match=i2
          END IF   ! plumes touch
        END DO !j2

      END DO   ! i2

    END DO  ! j

    IF (imatch == 1) THEN
      index_pl(i) = 0      ! plume copied
      ! relabel all last plume as matching plume
      nn=num_pl(iplume_match)
      num_pl(iplume_match)=num_pl(iplume_match)+ num_pl(i)
      DO j=1,num_pl(i)
        nn=nn+1
        coordx_pl(nn,iplume_match) = coordx_pl(j,i)
        coordy_pl(nn,iplume_match) = coordy_pl(j,i)
        bcu2(coordx_pl(j,i), coordy_pl(j,i)) = iplume_match
      END DO
    END IF

  END DO    !i


  ALLOCATE(plume_no(nups))

  nn=0
  DO i=1,nups
    IF (index_pl(i) == 1) THEN
      nn=nn+1
      plume_no(i) = nn        ! maps to new number

      ! Set plume number in plume_id

      DO j=1,num_pl(i)
        plume_id(coordx_pl(j,i), coordy_pl(j,i)) = nn ! plume number
      END DO

    END IF
  END DO

  DEALLOCATE(index_pl)
  DEALLOCATE(num_pl)
  DEALLOCATE(coordx_pl)
  DEALLOCATE(coordy_pl)

ELSE   ! Only one plume so no need to check whether touching

  nups2 = nups
  ALLOCATE(plume_no(nups2))
  plume_no(:) = 1

END IF  ! End of test for touching plumes

!----------------------------------------------------------------------------
! Analysis of plumes using bcu2
!----------------------------------------------------------------------------

IF (nups2 > 0 .AND. ncnt > 0) THEN

  ALLOCATE(num_pl(nups2))
  ALLOCATE(coordx_max(nups2))
  ALLOCATE(coordx_min(nups2))
  ALLOCATE(coordy_max(nups2))
  ALLOCATE(coordy_min(nups2))

  DO i=1,nups2
    num_pl(i) = 0
    coordx_max(i) = 0
    coordx_min(i) = 0
    coordy_max(i) = 0
    coordy_min(i) = 0
  END DO

  DO j=1,nrows
    DO i= 1,ncols
      IF (bcu2(i,j) > 0) THEN    ! plume point
        kk=plume_no(bcu2(i,j))
        num_pl(kk)=num_pl(kk)+1
        IF (j >  coordy_max(kk)) THEN
          coordy_max(kk) = j
        END IF
        IF (i >  coordx_max(kk)) THEN
          coordx_max(kk) = i
        END IF
        IF ( coordy_min(kk) == 0) THEN
          coordy_min(kk) = j
        ELSE IF (j < coordy_min(kk)) THEN
          coordy_min(kk) = j
        END IF
        IF ( coordx_min(kk) == 0) THEN
          coordx_min(kk) = i
        ELSE IF ( i < coordx_min(kk)) THEN
          coordx_min(kk) = i
        END IF
      END IF
    END DO
  END DO


  !-------------------------------------------------------------------------
  ! min circle diameter calculation from Paul Field
  ! -------------------------------------------------------------------------
  !  Max number of points per plume ?
  max_plume = num_pl(1)

  IF (nups2 > 1) THEN
    DO i=2,nups2
      IF (num_pl(i) > max_plume) THEN
        max_plume = num_pl(i)
      END IF
    END DO
  END IF

  ALLOCATE(x(max_plume,nups2))      ! max and min x & y
  ALLOCATE(y(max_plume,nups2))
  num_pl(:) = 0

  DO j=1,nrows
    DO i= 1,ncols
      IF (bcu2(i,j) > 0) THEN    ! plume point
        kk=plume_no(bcu2(i,j))
        num_pl(kk)=num_pl(kk)+1
        x(num_pl(kk),kk) = REAL(i)
        y(num_pl(kk),kk) = REAL(j)
      END IF
    END DO
  END DO

  DO i=1,nups2
    idx = coordx_max(i) - coordx_min(i) +1
    idy = coordy_max(i) - coordy_min(i) +1
    dxdy = REAL(idx)/REAL(idy)               ! orentation of plume
    fract = REAL(num_pl(i))/REAL(idx*idy)    ! fraction of square plume

    ! Tests to put in correct bins for fraction
    IF (fract == 1.0) THEN
      fract_pdf(11)= fract_pdf(11)+1.0
    ELSE
      DO j=1,10
        IF (fract >= fract_bin(j) .AND. fract < fract_bin(j+1)) THEN
          fract_pdf(j)= fract_pdf(j)+1.0
        END IF
      END DO
    END IF

    IF (num_pl(i) == 1) THEN

      diam_pdf(1) = diam_pdf(1) +1.0    ! one grid point
      size_pdf(1) = size_pdf(1) +1.0    ! one grid point

    ELSE IF (num_pl(i) == 2) THEN

      diam_pdf(2) = diam_pdf(2) +1.0    ! one grid point
      size_pdf(2) = size_pdf(2) +1.0    ! one grid point

      IF (x(2,i) < x(1,i)) THEN
        dxdy = -1.0* dxdy
      END IF

    ELSE IF (num_pl(i) > 2) THEN   ! all other plumes
      num_box = REAL(num_pl(i))
      ymin = y(1,i)
      ymax = y(1,i)
      xmin = x(1,i)
      xmax = x(1,i)
      i_xmax = 1
      i_xmin = 1
      i_ymax = 1
      i_ymin = 1
      i_st =  INT(x(1,i))
      sum_idx = 0

      DO j=2,num_pl(i)
        IF (y(j,i) < ymin) THEN
          ymin=y(j,i)
          i_ymin = j
        END IF
        IF (y(j,i) > ymax) THEN
          ymax=y(j,i)
          i_ymax = j
        END IF
        IF (x(j,i) > xmin) THEN
          xmin=x(j,i)
          i_xmin = j
        END IF
        IF (x(j,i) > xmax) THEN
          xmax=x(j,i)
          i_xmax = j
        END IF
        sum_idx = sum_idx + INT(x(j,i)) - i_st

      END DO
      ! sign for dx/dy
      IF (sum_idx < 0) THEN
        dxdy = -1.0* dxdy
      END IF

      ! Set xspan = distance between the 2 points xmin & xmax (squared) *
      dx=x(i_xmax,i)-x(i_xmin,i)
      dy=y(i_xmax,i)-y(i_xmin,i)
      xspan=dx*dx+dy*dy

      dx=x(i_ymax,i)-x(i_ymin,i)
      dy=y(i_ymax,i)-y(i_ymin,i)
      yspan=dx*dx+dy*dy
      !  Set points dia1 & dia2 to the maximally separated pair
      dia1=i_xmin  ! assume xspan biggest
      dia2=i_xmax
      maxspan = xspan
      IF (yspan > maxspan) THEN
        maxspan=yspan
        dia1=i_ymin
        dia2=i_ymax
      END IF

      ! dia1,dia2 is a diameter of initial sphere
      ! calc initial center
      cenx=(x(dia1,i)+x(dia2,i))/2.0
      ceny=(y(dia1,i)+y(dia2,i))/2.0
      ! calc initial radius
      dx= x(dia2,i) -cenx
      dy= y(dia2,i) -ceny
      rad_sq=dx*dx + dy*dy
      rad=SQRT(rad_sq)

      ! SECOND PASS: increment current disk

      DO j=1,num_pl(i)
        dx=x(j,i)-cenx
        dy=y(j,i)-ceny
        old_to_p_sq= dx*dx +dy*dy
        IF (old_to_p_sq > rad_sq) THEN
          ! point is outside current disk
          old_to_p = SQRT(old_to_p_sq);
          ! calc radius of new sphere
          rad = (rad + old_to_p) / 2.0;
          rad_sq = rad*rad       ! for next r**2 compare
          old_to_new = old_to_p - rad
          !calc center of new sphere
          cenx = (rad*cenx + old_to_new*x(j,i)) / old_to_p
          ceny = (rad*ceny + old_to_new*y(j,i)) / old_to_p
        END IF
      END DO
      diam = rad*2.0
      ! tests on bins for diameter PDF
      DO j=1,nbins_diam
        IF (diam > diam_bin(j) .AND. diam <= diam_bin(j+1)) THEN
          diam_pdf(j) =  diam_pdf(j)+1
        END IF
      END DO
      ! Size based on num_box
      DO j=1,nbins_size
        IF (num_box > size_bin(j) .AND. num_box <= size_bin(j+1)) THEN
          size_pdf(j) =  size_pdf(j)+1
        END IF
      END DO

    END IF
    ! Tests for dxdy
    IF ( dxdy <= -1.5) THEN
      dxdy_pdf(1)= dxdy_pdf(1)+1.0
    ELSE IF ( dxdy > -1.5 .AND. dxdy < -1.0) THEN
      dxdy_pdf(2)= dxdy_pdf(2)+1.0
    ELSE IF ( dxdy == -1.0) THEN
      dxdy_pdf(3)= dxdy_pdf(3)+1.0
    ELSE IF ( dxdy > -1.0 .AND. dxdy <= -0.5) THEN
      dxdy_pdf(4)= dxdy_pdf(4)+1.0
    ELSE IF ( dxdy > -0.5 .AND. dxdy < 0.0) THEN
      dxdy_pdf(5)= dxdy_pdf(5)+1.0
    ELSE IF ( dxdy > 0.0 .AND. dxdy < 0.5) THEN
      dxdy_pdf(6)= dxdy_pdf(6)+1.0
    ELSE IF ( dxdy >= 0.5 .AND. dxdy < 1.0) THEN
      dxdy_pdf(7)= dxdy_pdf(7)+1.0
    ELSE IF ( dxdy ==  1.0) THEN
      dxdy_pdf(8)= dxdy_pdf(8)+1.0
    ELSE IF ( dxdy > 1.0 .AND. dxdy <= 1.5) THEN
      dxdy_pdf(9)= dxdy_pdf(9)+1.0
    ELSE IF ( dxdy > 1.5 ) THEN
      dxdy_pdf(10)= dxdy_pdf(10)+1.0
    END IF

  END DO    ! loop over plumes

  DEALLOCATE(coordx_max)
  DEALLOCATE(coordx_min)
  DEALLOCATE(coordy_max)
  DEALLOCATE(coordy_min)

  DEALLOCATE(y)
  DEALLOCATE(x)

  DEALLOCATE(num_pl)

END IF

DEALLOCATE(plume_no)
!----------------------------------------------------------------------------
! Print info if required - helps check program
!----------------------------------------------------------------------------
IF (iprint == 2) THEN

  scut(:,:) = ' '     ! Initialise

  WRITE(umMessage,'(A,I10,A,2I6)') 'Number of buoyant points ',ncnt,        &
                           ' Number of plumes ',nups,nups2
  CALL umPrint(umMessage,src=RoutineName)

  ! Problem printing if more than 503 plumes as clet not setup with enough
  ! letters/symbols

  IF (ncnt > 1 .AND. nups < 504) THEN

    DO j=1,nrows
      DO i=1,ncols
        scut(i,j)=clet(bcu2(i,j))
      END DO
    END DO
    ! Choose format based on region (setup for values most likely).

    WRITE(umMessage,'(A)') 'Map of plumes'
    CALL umPrint(umMessage,src=RoutineName)
    SELECT CASE (ncols)
    CASE (600)
      DO j=1,nrows
        WRITE(umMessage,'(600A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (200)
      DO j=1,nrows
        WRITE(umMessage,'(200A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (120)
      DO j=1,nrows
        WRITE(umMessage,'(120A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (100)
      DO j=1,nrows
        WRITE(umMessage,'(100A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (50)
      DO j=1,nrows
        WRITE(umMessage,'(50A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (30)
      DO j=1,nrows
        WRITE(umMessage,'(30A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (20)
      DO j=1,nrows
        WRITE(umMessage,'(20A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE (10)
      DO j=1,nrows
        WRITE(umMessage,'(10A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    CASE DEFAULT
      DO j=1,nrows
        WRITE(umMessage,'(100A1)') (scut(i,j),i=1,ncols)
        CALL umPrint(umMessage,src=RoutineName)
      END DO
    END SELECT

  END IF    ! if possible to print plume maps

END IF

! Reset total number of plumes for output to the number of non-touching plumes.

nups= nups2

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE count_plumes

END MODULE count_plumes_mod
