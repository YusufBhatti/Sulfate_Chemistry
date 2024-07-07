! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE HINTCF -------------------------------------------------
!
!    Purpose :
!
!   THIS CALCULATES INTERPOLATION COEFFICIENTS FOR BI-LINEAR
!   INTERPOLATION OF VALUES FROM MODEL GRID POINTS TO OBSERVATION
!   POINTS.IT STORES THE COEFFICIENTS IN ANALYSIS WORK ARRAY ANWORK.
!   IN CALCULATING THE COEFFICIENTS,DISTANCES ARE MEASURED RELATIVE
!   TO THE NEAREST MODEL POINT TO THE OBSERVATION.INCREMENT VECTORS
!   (OVER OBSERVATIONS) ARE CALCULATED TO DEFINE THE POSITION
!   RELATIVE TO THE NEAREST POINT,OF THE FOUR POINTS SURROUNDING
!   THE OBSERVATION.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!   WORKS FOR ARAKAWA 'B' GRID ONLY
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE hintcf_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HINTCF_MOD'

CONTAINS

SUBROUTINE hintcf (lwind,lenob,obs_lat,obs_long,obs_no,           &
                   row_length,p_rows,                             &
                   cf1pt,cf2pt,cf3pt,cf4pt,                       &
                   np1pt,np2pt,np3pt,np4pt,                       &
                   icode,lambda_p,phi_p,cmessage)
!     --------------------------------------------------
!
!
USE conversions_mod, ONLY: pi_over_180, pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: l_regular
USE umprintmgr, ONLY: ummessage, umprint

IMPLICIT NONE
!     -----------------------------------------------------------
LOGICAL :: lwind              !IN switch to identify wind grid
INTEGER :: lenob              !IN number of obs
INTEGER :: row_length,p_rows  !IN size of model grid
!     -----------------------------------------------------------
REAL ::  obs_lat(lenob)    !IN ob co-latitudes between 0 & PI
REAL ::  obs_long(lenob)   !IN ob longitudes between 0 & 2*PI
!     NP#PT are pointers to the 4 gridpoints surrounding each ob.
INTEGER :: np1pt(lenob)      !OUT pointer to nearest model point
INTEGER :: np2pt(lenob)      !OUT pointer to same row, +or-1 pt.
INTEGER :: np3pt(lenob)      !OUT pointer to +or-1 row, same pt.
INTEGER :: np4pt(lenob)      !OUT pointer to +or-1 row, +or-1 pt.
!     CF#PT are the weights given to each pt pointed to by NP#PT.
REAL ::  cf1pt(lenob)      !OUT interpolation coeffs
REAL ::  cf2pt(lenob)      !OUT interpolation coeffs
REAL ::  cf3pt(lenob)      !OUT interpolation coeffs
REAL ::  cf4pt(lenob)      !OUT interpolation coeffs
INTEGER :: obs_no(lenob)     !IN pointers to obs

REAL :: lambda_p(1-halo_i:row_length+halo_i)
REAL :: phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)


!     -----------------------------------------------------------
INTEGER :: icode              !OUT error code and message
CHARACTER(LEN=errormessagelength) :: cmessage

!     -----------------------------------------------------------
! DECLARATIONS FOR VARIABLE GRID
REAL :: grid_lats_thispe(0:p_rows)
REAL :: grid_lons_thispe(0:row_length)
REAL :: r_delta_lats_thispe(0:p_rows)
REAL :: r_delta_lons_thispe(0:row_length)
REAL :: delta_lats_thispe(0:p_rows)
REAL :: delta_lons_thispe(0:row_length)
LOGICAL :: variable
INTEGER :: i,k
INTEGER :: locpos(1)
INTEGER :: n_row_test(lenob) ! nearest row to ob     minus 1
INTEGER :: n_pnt_test(lenob) ! nearest point to ob   minus 1
INTEGER :: iobs_temp(lenob)
REAL :: obs_long_test,obs_lat_test
!     -----------------------------------------------------------
! *   DYNAMIC ALLOCATION WITH LENOB
REAL ::        wklon(lenob)
REAL ::        work5(lenob)
!     CF#PT and NP1PT are used for workspace during the calculation
INTEGER :: n_row(lenob) ! nearest row to ob     minus 1
INTEGER :: n_pnt(lenob) ! nearest point to ob   minus 1
INTEGER :: i_row(lenob) ! increment to row other side of ob   (+or-1)
INTEGER :: i_pnt(lenob) ! increment to point other side of ob (+or-1)
!     -----------------------------------------------------------
REAL ::                                                           &
               zdlat,zdlong,zlatn,zlongw,r_zdlat,r_zdlong,zlonge
INTEGER :: job
!     -----------------------------------------------------------
REAL :: pi2p
PARAMETER (pi2p = 2.0*pi)
!     -----------------------------------------------------------
REAL :: TINY

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HINTCF'

PARAMETER (TINY = 1.6e-7)
REAL :: Tiny1
PARAMETER (Tiny1 = 9.9e-13)
!     -----------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('HINTCF  ',3)

!
! PART 1
! ------
! SET UP LATITUDE OF NORTHERN ROW AND LONGITUDE OF WESTERN POINT
! FOR EACH MODEL GRID
!
! PARAMETERS FOR MODEL GRID
!
!     GET MODEL GRID INFORMATION FROM COMMG AND CONVERT TO RADIANS
!     CONVERT XLATN TO CO-LATITUDE FIRST
!
zdlat  = dlat         * pi_over_180
r_zdlat = 1.0/zdlat
zdlong = dlong        * pi_over_180
r_zdlong = 1.0/zdlong
zlatn  = (90.0-lat_n) * pi_over_180
zlongw = long_w_model * pi_over_180
zlonge = long_e_model * pi_over_180

!     MAKE SURE ZLONGW,ZLONGE ARE IN RANGE 0-->2*PI
IF (zlongw <  0.0) THEN  
  zlongw=zlongw+pi2p
ELSE IF (zlongw >  pi2p) THEN
  zlongw=zlongw-pi2p
END IF
IF (zlonge <  0.0) THEN  
  zlonge=zlonge+pi2p
ELSE IF (zlonge >  pi2p) THEN
  zlonge=zlonge-pi2p
END IF

! Set up cooridinates of grid (in radians)
! If a regular grid these are calculated from input parameters otherwise
! copy from arrays phi_p and lambda_p
! Also note that variable grids have monotonically increasing longitudes

IF (l_regular) THEN
  DO i=1,p_rows
    grid_lats_thispe(i)=(90.0-lat_n)+REAL(i-1)*dlat
    grid_lats_thispe(i)=grid_lats_thispe(i) * pi_over_180

    ! spacing and reciprocal between row
    delta_lats_thispe(i)= dlat * pi_over_180
    r_delta_lats_thispe(i)=1.0/delta_lats_thispe(i)
  END DO

  DO i=1,row_length
    grid_lons_thispe(i)=long_w_model+(i-1)*dlong
    ! allow for Meridian
    IF (grid_lons_thispe(i) >  360.0) THEN
      grid_lons_thispe(i)=grid_lons_thispe(i)-360.0
    END IF
    grid_lons_thispe(i)=grid_lons_thispe(i) * pi_over_180

    ! spacing and reciprocal between cols
    delta_lons_thispe(i)=dlong * pi_over_180
    r_delta_lons_thispe(i)=1.0/delta_lons_thispe(i)
  END DO

ELSE ! variable grid

  DO i=1,row_length
    grid_lons_thispe(i)=lambda_p(i)
  END DO
  DO i=1,p_rows
    ! colatitudes in radians are required
    grid_lats_thispe(i)=phi_p(1,p_rows-i+1)/pi_over_180
    grid_lats_thispe(i)=90.0-grid_lats_thispe(i)
    grid_lats_thispe(i)=grid_lats_thispe(i)*pi_over_180
  END DO

  ! Calculate reciprocal latitude spacings between rows
  DO i=1,p_rows-1
    delta_lats_thispe(i)=grid_lats_thispe(i+1)-grid_lats_thispe(i)
    r_delta_lats_thispe(i)=1.0/delta_lats_thispe(i)
  END DO
  delta_lats_thispe(p_rows)=delta_lats_thispe(p_rows-1)
  r_delta_lats_thispe(p_rows)=r_delta_lats_thispe(p_rows-1)

  ! Calculate reciprocal of longitude spacings between columns
  DO i=1,row_length-1
    delta_lons_thispe(i)=grid_lons_thispe(i+1)-grid_lons_thispe(i)
    r_delta_lons_thispe(i)=1.0/delta_lons_thispe(i)
  END DO
  delta_lons_thispe(row_length)=delta_lons_thispe(row_length-1)
  r_delta_lons_thispe(row_length)=r_delta_lons_thispe(row_length-1)

END IF ! l_regular

! Calculate zeroth row/col details
delta_lats_thispe(0)=delta_lats_thispe(1)
delta_lons_thispe(0)=delta_lons_thispe(1)
r_delta_lats_thispe(0)=r_delta_lats_thispe(1)
r_delta_lons_thispe(0)=r_delta_lons_thispe(1)
grid_lats_thispe(0)=grid_lats_thispe(1)-delta_lats_thispe(0)
grid_lons_thispe(0)=grid_lons_thispe(1)-delta_lons_thispe(0)
!
!
! PART 2
! ------
! CALCULATE (NEAREST-1) MODEL GRID PT TO OBSERVATION
!
!     INPUT  VECTOR IN OBS_LAT  - COLATITUDE (RADIANS) SET UP IN AC
!                   IN OBS_LONG - LONGITUDE  (RADIANS) SET UP IN AC
!
!     OUTPUT VECTORS
!
!     OUTPUT VECTORS OF INTERPOLATION COEFFICIENTS ARE IN
!     CF1PT CF2PT CF3PT CF4PT NP1PT NP2PT NP3PT NP4PT

! The TINY1 offset is required to prevent invalid
! NP1PT NP2PT NP3PT NP4PT on variable grids (-ve values)
!
! 2.1
! ---
!
DO job=1,lenob

  !       FIND NEAREST MODEL ROW minus 1.  in N_ROW
  IF (l_regular) THEN
    n_row(job) = NINT((obs_lat(job)-zlatn)*r_zdlat)
  ELSE
    obs_lat_test=obs_lat(job)-tiny1

    n_row(job)=p_rows-1
    IF (obs_lat_test <  grid_lats_thispe(1)) THEN
      n_row(job)=0
    ELSE

      DO k=1,p_rows-1
        IF (obs_lat_test  >=  grid_lats_thispe(k) .AND.                 &
           obs_lat_test  <   grid_lats_thispe(k+1)) THEN
          obs_lat_test=grid_lats_thispe(k+1)-obs_lat_test
          ! To mimic regular grid check which half of grid box point lies in
          IF (obs_lat_test <  delta_lats_thispe(k)/2.0) THEN
            n_row(job)=k
          ELSE
            n_row(job)=k-1
          END IF
          EXIT
        END IF
      END DO
    END IF
  END IF
  !
  ! 2.2
  ! ---
  ! FIND NEAREST POINT minus 1 ON ROW.   in N_PNT
  !
  !
  !       SPECIAL TREATMENT WHEN LIMITED AREA SPANS MERIDIAN.
  !       OBS EAST OF THE MERIDIAN HAVE LONGITUDES < WESTERN BOUNDARY.
  !       CALCULATE THEIR LONGITUDE DISPLACEMENT WITH EXTRA 2*PI SHIFT.
  !
  !       INITIALISE WORK ARRAY OF LONGITUDE SHIFTS
  wklon(job)=0.0
  IF (l_regular) THEN
    IF (zlongw  >   zlonge) THEN
      !         ASSUME AREA STRADDLES MERIDIAN
      !         SET UP LONGITUDE SHIFT FOR OBS BETWEEN MERIDIAN
      !                                      AND EASTERN BOUNDARY.
      IF (obs_long(job) >= 0.0 .AND. obs_long(job)                  &
                    <= (zlonge+zdlong/2.0))   wklon(job)= pi2p
    END IF
    n_pnt(job) = NINT((obs_long(job)+wklon(job)-zlongw)*r_zdlong)
  ELSE
    obs_long_test=obs_long(job)-tiny1

    ! Check to see if any part of this PE crosses Meridan
    ! (Variable grids work in monotonically increasing longitude)
    IF (MAXVAL(grid_lons_thispe) > pi2p) THEN
      ! Observations east of the Meridian need 2PI added to longitude
      IF ((obs_long_test-grid_lons_thispe(0)) <= 0.0) THEN
        wklon(job)= pi2p
        obs_long_test=obs_long_test+wklon(job)
      END IF
    END IF

    n_pnt(job)=row_length-1
    IF (obs_long_test <  grid_lons_thispe(1)) THEN
      n_pnt(job)=0
    ELSE

      DO k=1,row_length-1
        IF (obs_long_test  >=  grid_lons_thispe(k) .AND.                 &
           obs_long_test  <   grid_lons_thispe(k+1)) THEN
          obs_long_test=grid_lons_thispe(k+1)-obs_long_test
          ! To mimic regular grid check which half of grid box point lies in
          IF (obs_long_test  <   delta_lons_thispe(k)/2.0) THEN
            n_pnt(job)=k
          ELSE
            n_pnt(job)=k-1
          END IF
          EXIT
        END IF
      END DO
    END IF
  END IF
  !
  ! PART 3
  ! ------
  ! CALCULATE INTERPOLATION COEFFS. In lat in CF1PT, in long in CF2PT.
  !
  IF (l_regular) THEN
    cf1pt(job) = obs_lat(job) -(zlatn +n_row(job)*zdlat )
    cf2pt(job) = obs_long(job)-(zlongw+n_pnt(job)*zdlong)         &
                 + wklon(job)
  ELSE
    IF (n_row(job) <  p_rows) THEN
      cf1pt(job)=obs_lat(job)-grid_lats_thispe(n_row(job)+1)
    ELSE
      cf1pt(job)=1.0  ! can this happen?
    END IF


    IF (n_pnt(job) <  row_length) THEN
      cf2pt(job)=obs_long(job)-grid_lons_thispe(n_pnt(job)+1)     &
                + wklon(job)
    ELSE
      cf2pt(job)=1.0   ! can this happen?
    END IF
  END IF
  !
  ! GET INCREMENT VECTORS TO OBTAIN FOUR SURROUNDING GRID POINTS
  ! INCLUDE POSSIBILITY OF OBS BEING ON SOUTHERN OR EASTERN HALO
  ! BOUNDARIES. (Cannot be on northern or western because of the way
  ! obs are distributed in RDOBS)
  !
  !   N.B Value of Tiny is chosen so that observations within
  !   approx 1 metre of grid edge are affected.

  !
  IF (cf1pt(job) >= 0.0) THEN
    i_row(job)=1
    !         Set I_ROW to 0 if CF1PT is very small.
    IF (cf1pt(job)  <   TINY) THEN
      i_row(job) = 0
    END IF
  ELSE
    i_row(job)=-1
  END IF
  !
  IF (cf2pt(job) >= 0.0) THEN
    i_pnt(job) = 1
    !         Set I_PNT to 0 if CF2PT is very small.
    IF (cf2pt(job)  <   TINY) THEN
      i_pnt(job) = 0
    END IF
  ELSE
    i_pnt(job) = -1
  END IF
  !
  cf1pt(job) = ABS(cf1pt(job))
  cf2pt(job) = ABS(cf2pt(job))

  !
  !       Make sure that all pointers are within grid, and different
  np1pt(job)=n_row(job)+i_row(job) ! use NP1PT as workspace
  IF (np1pt(job) <  0) THEN         ! row-1 is outside grid
    i_row(job)=1                   ! point instead to row+1
    cf1pt(job)=0.0                  ! but do not use row+1.
  END IF

  IF (np1pt(job) >= p_rows) THEN    ! row+1 is outside grid
    i_row(job)=-1                  ! point instead to row-1
    cf1pt(job)=0.0                  ! but do not use row-1.
  END IF
  IF ((n_pnt(job) + i_pnt(job)) <  1) THEN
    i_pnt(job) = 1
    cf2pt(job) = 0.0
  END IF

  IF ((n_pnt(job) + i_pnt(job)) >=  row_length) THEN
    i_pnt(job) = -1
    cf2pt(job) = 0.0
  END IF

  n_row(job)=MAX(n_row(job),0)
  n_pnt(job)=MAX(n_pnt(job),0)
  !
  !       Convert to 2-D coeffs for 4 surrounding gridpoints
  IF (l_regular) THEN
    cf4pt(job)  = cf1pt(job)*r_zdlat   ! lat coeff for row +or-1
    work5 (job) = cf2pt(job)*r_zdlong  ! long coeff for pt +or-1
  ELSE
    cf4pt(job) = cf1pt(job)*r_delta_lats_thispe(n_row(job))
    work5(job) = cf2pt(job)*r_delta_lons_thispe(n_row(job))
  END IF
  cf2pt(job)  = 1.0-cf4pt(job)       ! lat  coeff for row
  cf3pt(job)  = 1.0-work5(job)       ! long coeff for pt
  !
  cf1pt(job) = cf2pt(job)*cf3pt(job) ! row, pt
  cf2pt(job) = cf2pt(job)*work5(job) ! row, pt +or-1
  cf3pt(job) = cf4pt(job)*cf3pt(job) ! row +or-1, pt
  cf4pt(job) = cf4pt(job)*work5(job) ! row +or-1, pt +or-1
  !
  !  up to now have assumed north most row is first
  !  but ND order is reversed, so adjust use of N_ROW
  np1pt(job) =n_pnt(job) + 1 +                                    &
                       (p_rows -1 - n_row(job) )*row_length
  np2pt(job) =np1pt(job) + i_pnt(job)

  !  respecting S->N order of ND fields,
  !  use -I_ROW to point to other row adjacent to ob
  np3pt(job) =np1pt(job) - i_row(job)*row_length
  np4pt(job) =np3pt(job) + i_pnt(job)

  IF (np1pt(job) >  row_length*p_rows) THEN
    cmessage = 'NP1PT(JOB) too big!!'
    WRITE(umMessage,*)'NP1PT(JOB) too big!!',np1pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF
  IF (np1pt(job) <  1) THEN
    cmessage = 'NP1PT(JOB) too small!!'
    WRITE(umMessage,*)'NP1PT(JOB) too small!!',np1pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF

  IF (np2pt(job) >  row_length*p_rows) THEN
    cmessage = 'NP2PT(JOB) too big!!'
    WRITE(umMessage,*)'NP2PT(JOB) too big!!',np2pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF
  IF (np2pt(job) <  1) THEN
    cmessage = 'NP2PT(JOB) too small!!'
    WRITE(umMessage,*)'NP2PT(JOB) too small!!',np2pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF

  IF (np3pt(job) >  row_length*p_rows) THEN
    cmessage = 'NP3PT(JOB) too big!!'
    WRITE(umMessage,*)'NP3PT(JOB) too big!!',np3pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF
  IF (np3pt(job) <  1) THEN
    cmessage = 'NP3PT(JOB) too small!!'
    WRITE(umMessage,*)'NP3PT(JOB) too small!!',np3pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF

  IF (np4pt(job) >  row_length*p_rows) THEN
    cmessage = 'NP4PT(JOB) too big!!'
    WRITE(umMessage,*)'NP4PT(JOB) too big!!',np4pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF
  IF (np4pt(job) <  1) THEN
    cmessage = 'NP4PT(JOB) too small!!'
    WRITE(umMessage,*)'NP4PT(JOB) too small!!',np4pt(job),job,mype
    CALL umPrint(umMessage,src='hintcf')
    icode=1
  END IF

END DO ! JOB
!
!
IF (ltimer_ac) CALL timer('HINTCF  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hintcf
END MODULE hintcf_mod
