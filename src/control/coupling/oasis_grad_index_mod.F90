MODULE OASIS_Grad_Index_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

IMPLICIT NONE

CONTAINS

SUBROUTINE OASIS_Grad_Index()
!
! Description:
! Set up indexing arrays for use in calculation of field gradients
! to be used when employing 2nd order conservative remapping
! with the OASIS-MCT coupler.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!=======================================================================

USE UM_ParVars,   ONLY: at_extremity
USE UM_ParParams, ONLY: PNorth, PSouth, PWest, PEast
USE Field_Types,  ONLY: fld_type_p
USE um_types, ONLY: integer32, real64
USE oasis_atm_data_mod, ONLY: GradIndex, oasis_imt, oasis_jmt, fland_loc
USE conversions_mod, ONLY: pi
USE model_domain_mod, ONLY : mt_global, model_type
USE mpp_conf_mod, ONLY: swap_field_is_scalar

USE trignometric_mod, ONLY: true_longitude, true_latitude

IMPLICIT NONE

REAL (KIND=real64) :: lat_im1, lat_ip1      ! Temp values holding from and to
REAL (KIND=real64) :: lon_im1, lon_ip1      ! lats and longs
REAL (KIND=real64) :: lat_jm1, lat_jp1      
REAL (KIND=real64) :: lon_jm1, lon_jp1      

REAL (KIND=real64) :: temp            ! Holds intermediate result

INTEGER (KIND=integer32) :: i, j

REAL (KIND=real64) :: radian_lat(0:oasis_imt+1,0:oasis_jmt+1)
REAL (KIND=real64) :: radian_lon(0:oasis_imt+1,0:oasis_jmt+1)

! The true_latitude arrays etc. don't contain halo rows and on global grids 
! don't even contain usable longitudes in the top and bottom rows. 
! So we have to do extra work here to make our own copy
! with halos and to get sensible longitudes in ALL true domain grid points. 

radian_lat(1:oasis_imt,1:oasis_jmt) =  true_latitude(1:oasis_imt,1:oasis_jmt) 
radian_lon(1:oasis_imt,1:oasis_jmt) =  true_longitude(1:oasis_imt,1:oasis_jmt) 

! Ensure our radian latitude halos are up to date.
CALL Swap_Bounds(radian_lat, oasis_imt, oasis_jmt, 1,        &
          1, 1, fld_type_p,swap_field_is_scalar)

radian_lon(1:oasis_imt,1:oasis_jmt) =  true_longitude(1:oasis_imt,1:oasis_jmt)

! Ensure our radian longitude halos are up to date.
CALL Swap_Bounds(radian_lon, oasis_imt, oasis_jmt, 1,        &
          1, 1, fld_type_p,swap_field_is_scalar)

! Then if it's a global model we have to get sensible longitudes in 
! the top and bottom rows!
IF (model_type  ==  mt_Global ) THEN

  IF (at_extremity(PSouth)) THEN

    DO i = 0,oasis_imt+1
      radian_lon(i,1) = radian_lon(i,2)
    END DO 

  END IF

  IF (at_extremity(PNorth)) THEN

    DO i = 0,oasis_imt+1
      radian_lon(i,oasis_jmt) = radian_lon(i,oasis_jmt-1)
    END DO 

  END IF

END IF

! Adjust longitude values at extreme grid edges to make gradient calculations 
! sensible for distances between grid points in the same vicinity.
IF (at_extremity(PWest)) THEN
  DO j = 0,  oasis_jmt+1
    IF (radian_lon(0,j) > radian_lon(1,j)) THEN
      radian_lon(0,j) = radian_lon(0,j) - (2.0 * pi)
    END IF
  END DO
END IF

IF (at_extremity(PEast)) THEN
  DO j = 0,  oasis_jmt+1
    IF (radian_lon(oasis_imt+1,j) < radian_lon(oasis_imt,j)) THEN
      radian_lon(oasis_imt+1,j) = radian_lon(oasis_imt+1,j) + (2.0 * pi)
    END IF
  END DO
END IF

! Our gradient arrays must be large enough to cater for halos
ALLOCATE(GradIndex(0:oasis_imt+1,0:oasis_jmt+1))

DO j = 0, oasis_jmt+1
  DO i = 0, oasis_imt+1

    GradIndex(i,j)%ip1 = i
    GradIndex(i,j)%im1 = i
    GradIndex(i,j)%jp1 = j
    GradIndex(i,j)%jm1 = j

    GradIndex(i,j)%dlat_row = 0.0
    GradIndex(i,j)%dlat_col = 0.0
    GradIndex(i,j)%dlon_row = 0.0
    GradIndex(i,j)%dlon_col = 0.0

    GradIndex(i,j)%di = 0.0
    GradIndex(i,j)%dj = 0.0
    GradIndex(i,j)%Rdi = 0.0
    GradIndex(i,j)%Rdj = 0.0

  END DO
END DO

DO j = 1, oasis_jmt
  DO i = 1, oasis_imt
    ! If this point is a land point then don't bother with
    ! the gradient - it has no meaning for A-O coupling and we certainly
    ! don't want to calculate it through land from adjacent sea points.
    ! We also need to cater for any surrounding land
    ! points. 
    ! Unlike original OASIS3, where gradient calculations were calculated
    ! by the coupler, we deliberately don't take the "best we can do"
    ! approach (i.e. we must span cells at i-1 to i+1 or from j-1 to j+1
    ! to bother calculating the gradient), rather, we explicitly use 1st 
    ! order, i.e. no gradient terms at all, on points bordering any land 
    ! and at the N and S poles since this is generally thought to 
    ! produce more realistic results.   
    IF ((fland_loc(i,j) == 1.0)   .OR. &
        (fland_loc(i-1,j) == 1.0) .OR. & 
        (fland_loc(i+1,j) == 1.0) .OR. &
        (fland_loc(i,j-1) == 1.0) .OR. & 
        (fland_loc(i,j+1) == 1.0)) THEN
      GradIndex(i,j)%im1 = i
      GradIndex(i,j)%ip1 = i
      GradIndex(i,j)%jm1 = j
      GradIndex(i,j)%jp1 = j
    ELSE
      ! Set start and end indices for
      ! use in gradient calculations.
      GradIndex(i,j)%im1 = i-1
      GradIndex(i,j)%ip1 = i+1

      GradIndex(i,j)%jm1 = j-1
      GradIndex(i,j)%jp1 = j+1

    END IF ! This is not on or adjacent to a land point
  END DO  ! over i
END DO ! over j

! We need to cater for the fact that grids may be rotated by finding the
! true "as-the-crow-flies" distance between points in the i and j
! orientation (in terms radians on the arc.)
! This won't be quick because of all the trig functions etc but we only
! do this once at startup so it's not a major issue.
DO j = 1, oasis_jmt
  DO i = 1, oasis_imt
      
    lat_im1 = radian_lat(GradIndex(i,j)%im1,j)
    lat_ip1 = radian_lat(GradIndex(i,j)%ip1,j)

    lat_jm1 = radian_lat(i,GradIndex(i,j)%jm1)
    lat_jp1 = radian_lat(i,GradIndex(i,j)%jp1)

    lon_im1 = radian_lon(GradIndex(i,j)%im1,j)
    lon_ip1 = radian_lon(GradIndex(i,j)%ip1,j)

    lon_jm1 = radian_lon(i,GradIndex(i,j)%jm1)
    lon_jp1 = radian_lon(i,GradIndex(i,j)%jp1)

    ! We need to cater for two potential components of the latitude 
    ! differences since our i and j indexing won't necessarily follow
    ! latitude and longitude. So we take lat and long differences 
    ! across both i points and j points. 
    GradIndex(i,j)%dlat_row = lat_ip1 - lat_im1
    GradIndex(i,j)%dlat_col = lat_jp1 - lat_jm1

    ! Similarly we need to cater for two components of the 
    ! longitude differences. We also need to scale longitude
    ! separation by a cosine latitude factor.  
    GradIndex(i,j)%dlon_row = (lon_ip1 - lon_im1) * COS(radian_lat(i,j))
    GradIndex(i,j)%dlon_col = (lon_jp1 - lon_jm1) * COS(radian_lat(i,j))

    ! Work out the distance between our start and end points 
    ! following true longitude.
    ! Haversine formula should be OK for these purposes.
    temp = SQRT( ((SIN((lat_ip1-lat_im1)/2.0))**2) +  &
                   COS(lat_im1)*COS(lat_ip1)*         &
                 ((SIN((lon_ip1-lon_im1)/2.0))**2) )

    ! Calculate W-E distance in radians
    GradIndex(i,j)%di = 2 * ASIN(temp)

    ! Work out the distance between our start and end points 
    ! following true latitude.
    ! Haversine formula should be OK for these purposes.
    temp = SQRT( ((SIN((lat_jp1-lat_jm1)/2.0))**2) +  &
                   COS(lat_jm1)*COS(lat_jp1)*         &
                 ((SIN((lon_jp1-lon_jm1)/2.0))**2) )

    GradIndex(i,j)%dj = 2 * ASIN(temp)

    ! Calculate reciprocal distances for performance purposes because we
    ! dont want to do repeated divisions later on.
    ! If there's a zero distance, reciprocal remains set to zero!
    IF (ABS(GradIndex(i,j)%di) > 0.0) THEN
      GradIndex(i,j)%Rdi=1/GradIndex(i,j)%di
    END IF

    IF (ABS(GradIndex(i,j)%dj) > 0.0) THEN
      GradIndex(i,j)%Rdj=1/GradIndex(i,j)%dj
    END IF
 
  END DO  ! over i

END DO ! over j

END SUBROUTINE OASIS_Grad_Index

END MODULE OASIS_Grad_Index_mod
