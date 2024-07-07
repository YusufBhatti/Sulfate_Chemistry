! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculate simple test diagnostics based on a simple analytic formula
!
! Subroutine Interface:
SUBROUTINE Testdiag(                                              &
  p_field,v_field,rows,n_rows,row_length                          &
 ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole     &
 ,elf                                                             &
 ,press_levels_list,no_press_levels                               &
 ,model_levels_list,no_model_levels,forecast_hrs                  &
 ,diag1,diag2,diag3,diag4                                         &
 ,qdia1,qdia2,qdia3,qdia4)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon

USE horiz_grid_mod, ONLY: cartesian_grid
IMPLICIT NONE
!
! Description:

! Calculate simple test diagnostics based on a simple analytic formula:
!
!    VALUE=A*(LATITUDE+90.)+B*LONGITUDE+C*LEVEL+D*FORECAST_HRS
!    where A=1.0, B=1.0E2, C=1.0E3, D=1.0E4
!    and (LAT,LONG) are in degrees, actual position (rotated for LAM),
!    LEVEL is either the model level (real number) or
!                        pressure level (mb),
!    FORECAST_HRS in T+hours after assimilation time.
!    Theses diagnostics are to be used for checking output procedures
!    for various post-processing routes.
!    Four diagnostics are supported:
!    1. single-level FIELD (LEVEL=0.) at V points of Arakawa C grid
!    2. single-level FIELD (LEVEL=0.) at P points of Arakawa C grid
!    3. multi -level FIELD (LEVEL=press level) at P points of C grid
!    4. multi -level FIELD (LEVEL=model level) at P points of C grid
!
! Method:
!      Calculate lat,long offsets for MPP
!   1. Calculate FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!   2. Calculate ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!   3. Calculate SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!   4. Calculate THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!   5. Calculate FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!
!    DOCUMENTATION:  UM Doc Paper D7
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::                                                        &
  p_field                                                         &
                   !IN   Horizontal field size p points
, v_field                                                         &
                   !IN   Horizontal field size v points
, rows                                                            &
                   !IN   No. of rows for p field
, n_rows                                                          &
                   !IN   No. of rows for v field
, row_length                                                      &
                   !IN   No. of points per row
, no_model_levels                                                 &
                   !IN   model levels for output
, no_press_levels
                   !IN   press levels for output

! UM6.5 - MODEL_ANALYSIS_HRS chnged to real -
!                   requires FORECAST_HRS changed to REAL also
REAL :: forecast_hrs     !IN   FORECAST HOURS T+0, etc

LOGICAL ::                                                        &
  elf                                                             &
                   !IN  TRUE IF MODEL IS LAM WITH ROTATED GRID
 ,qdia1                                                           &
                   !IN  STASHflag for DIAG1
 ,qdia2                                                           &
                   !IN  STASHflag for DIAG2
 ,qdia3                                                           &
                   !IN  STASHflag for DIAG3
 ,qdia4            !IN  STASHflag for DIAG4

REAL ::                                                           &
  ew_space                                                        &
                   !IN  DELTA LONGITUDE (DEGREES)
, ns_space                                                        &
                   !IN  DELTA  LATITUDE (DEGREES)
, first_lat                                                       &
                   !IN  latitude  of first p row (degrees)
, first_long                                                      &
                   !IN  longitude of first p col (degrees)
, phi_pole                                                        &
                   !IN  latitude  of the pseudo pole
, lambda_pole      !IN  longitude of the pseudo pole

!   Array  arguments with intent(in):
REAL ::                                                           &
  model_levels_list(no_model_levels)                              &
                                     !IN LEVELS list (for DIAG3)
, press_levels_list(no_press_levels) !IN LEVELS list (for DIAG4)

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
REAL ::                                                           &
  diag1(v_field)                                                  &
                                  !OUT DIAGNOSTIC 1
, diag2(p_field)                                                  &
                                  !OUT DIAGNOSTIC 2
, diag3(p_field,no_press_levels)                                  &
                                  !OUT DIAGNOSTIC 3
, diag4(p_field,no_model_levels)  !OUT DIAGNOSTIC 4


! Local parameters:
REAL ::                                                           &
  a,b,c,d  ! COEFFICIENTS FOR CALCULATING VALUES OF FIELD
!
PARAMETER(a=1.0,b=1.0e2,c=1.0e3,d=1.0e4)

! Local scalars:
INTEGER ::                                                        &
  i,j,k                                                           &
               !  LOOP COUNTERS
 ,l                                                               &
               !  LOOP INDEX
 ,offsetx                                                         &
               !  INDEX OFFSETs (for calculating MPP lat,longs)
 ,offsety

! Local dynamic arrays:
REAL ::                                                           &
  latitude(p_field)                                               &
                    ! latitude  in degrees
, longitude(p_field)                                              &
                    ! longitude in degrees
, lat(p_field)                                                    &
                    ! latitude  in degrees on equatorial grid
, lon(p_field)      ! longitude in degrees on equatorial grid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TESTDIAG'


!- End of header


!
!  Calculate lat,long offsets for MPP
!

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
offsetx=datastart(1)- Offx - 1
offsety=datastart(2)- Offy - 1

!-----------------------------------------------------------------------
!   1. CALCULATE FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
IF (qdia1) THEN
  !
  !   1a. FIND EQUATORIAL LATITUDES,LONGITUDES
  !
  DO j=1,n_rows
    DO i=1,row_length
      l= i + (j-1)*row_length
      lat(l)=  first_lat  - ns_space*(j+offsety-0.5)
      lon(l)=  first_long + ew_space*(i+offsetx-0.5)
    END DO
  END DO
  !
  !   1b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
  !
  IF (elf .AND. .NOT. cartesian_grid) THEN
    CALL rotate_eq_to_latlon(lat, lon, latitude, longitude,  &
                             phi_pole, lambda_pole, v_field)
  ELSE
    DO i=1,v_field
      latitude(i) =lat(i)
      longitude(i)=lon(i)
    END DO
  END IF
  !
  !   1c. CALCULATE VALUE FROM ANALYTIC FUNCTION
  !

  DO i=1,v_field
    diag1(i)=a*(latitude(i)+90.0) + b*longitude(i) +             &
             d*forecast_hrs
  END DO

END IF               ! END OF qdia1 TEST

!-----------------------------------------------------------------------
!   2. CALCULATE ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!-----------------------------------------------------------------------
IF (qdia2 .OR. qdia3 .OR. qdia4) THEN
  !
  !   2a. FIND EQUATORIAL LATITUDES,LONGITUDES
  !
  DO j=1,rows
    DO i=1,row_length
      l= i + (j-1)*row_length
      lat(l)=  first_lat  - ns_space*(j+offsety-1)
      lon(l)=  first_long + ew_space*(i+offsetx-1)
    END DO
  END DO
  !
  !   2b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
  !
  IF (elf .AND. .NOT. cartesian_grid) THEN
    CALL rotate_eq_to_latlon(lat, lon, latitude, longitude,  &
                             phi_pole, lambda_pole, p_field)
  ELSE
    DO i=1,p_field
      latitude(i) =lat(i)
      longitude(i)=lon(i)
    END DO
  END IF
END IF                      ! END OF qdia2-4 TEST
!-----------------------------------------------------------------------
!   3. CALCULATE SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
IF (qdia2) THEN

  DO i=1,p_field
    diag2(i)=a*(latitude(i)+90.0) + b*longitude(i) +             &
             d*forecast_hrs
  END DO

END IF               ! END OF qdia2 TEST
!-----------------------------------------------------------------------
!   4. CALCULATE THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!-----------------------------------------------------------------------
IF (qdia3) THEN

  DO k=1,no_press_levels
    DO i=1,p_field
      diag3(i,k)=a*(latitude(i)+90.0)   + b*longitude(i) +       &
                 c*press_levels_list(k) + d*forecast_hrs
    END DO
  END DO

END IF               ! END OF qdia3 TEST
!-----------------------------------------------------------------------
!   5. CALCULATE FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!-----------------------------------------------------------------------
IF (qdia4) THEN

  DO k=1,no_model_levels
    DO i=1,p_field
      diag4(i,k)=a*(latitude(i)+90.0)   + b*longitude(i) +       &
                 c*model_levels_list(k) + d*forecast_hrs
    END DO
  END DO

END IF               ! END OF qdia4 TEST

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Testdiag
