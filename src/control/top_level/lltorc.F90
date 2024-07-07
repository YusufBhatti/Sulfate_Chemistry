! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convert lat/long specification to row/column specification
! Subroutine Interface:

SUBROUTINE lltorc(grid_code,                                      &
                  north_lat_in,south_lat_in,                      &
                  west_long_in,east_long_in,                      &
                  start_row,end_row,start_col,end_col             &
#if defined(RECON)
                 ,recon_grid                                      &
#endif
                 )
#if defined(RECON)
USE Rcf_Grid_Type_Mod, ONLY:                                     &
    Grid_Type

USE Rcf_Global_To_Local_Mod, ONLY:                               &
    Rcf_Get_Fld_Type
#endif

USE stash_model_mod, ONLY:                                                    &
    h_a_firstlat, h_a_firstlong, h_a_nsspace, h_a_ewspace, h_global,          &
    h_a_polelat, h_a_polelong


USE ozone_inputs_mod, ONLY: zon_av_ozone
USE submodel_mod, ONLY: atmos_im
USE Ereport_Mod, ONLY: Ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE Decomp_DB
USE decomp_params, ONLY: decomp_standard_atmos

USE model_domain_mod, ONLY: output_grid_stagger,                   &
                            FH_GridStagger_C

USE grdtypes_mod, ONLY: gt_unset, gt_atmos, gt_ocean, gt_wave,     &
                        gt_thetamass, gt_velocity, gt_u_c, gt_v_c, &
                        gt_hybrid, gt_river, gt_allpts, gt_land,   &
                        gt_sea, gt_full, gt_zonal, gt_meridional,  &
                        gt_ozone, gt_scalar, gt_compressed, gt_lbc,&
                        gt_nocyclic, gt_optcyclic, gt_cyclic

IMPLICIT NONE

! Description:
!   Uses the gridpoint code, the lat/long spec of the required area,
!   and the model area sizes to calculate the model row/column
!   numbers for the area.
!   Called by addrln, inputl, outptl, prelim, rcf_address_length.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 95
!    Written to UMDP 003 standards.
!

! Subroutine arguments:


INTEGER ::                                                        &
  grid_code                                                       &
              ! IN : Grid type code from STASHmaster

, north_lat_in                                                    &
                 ! IN : Latitude of Northern boundary
, south_lat_in                                                    &
                 ! IN : Latitude of Southern boundary
, west_long_in                                                    &
                 ! IN : Longitude of Western boundary
, east_long_in                                                    &
                 ! IN : Longitude of Eastern boundary

, start_row                                                       &
              ! OUT : Row number of start of area
, end_row                                                         &
              ! OUT : Row number of end of area
, start_col                                                       &
              ! OUT : Column number of start of area
, end_col     ! OUT : Column number of end of area

#if defined(RECON)
TYPE (grid_type) recon_grid
#endif

! Local variables
INTEGER ::                                                        &
  model_type                                                      &
              ! model type of grid
, content                                                         &
              ! content type of grid
, coverage                                                        &
              ! coverage of grid
, domain                                                          &
              ! domain of grid
, cyclic                                                          &
              ! does grid contain cyclic wrap columns
, fld_type    ! P, U or V field type

INTEGER ::                                                        &
  north_lat                                                       &
             ! Modifiable copies
, south_lat                                                       &
             ! of the input arguments
, west_long                                                       &
             ! so that the format can
, east_long                                                       &
             ! be changed if necessary

, nrows                                                           &
             ! number of rows
, ncols      ! number of columns

LOGICAL ::                                                        &
  lam_model  ! True if this is a LAM configuration

REAL ::                                                           &
  pole_lat                                                        &
             ! Latitude of rotated pole (LAM only)
, pole_long  ! Longitude of rotated pole (LAM only)

INTEGER ::                                                        &
  row_ordering                                                    &
                ! ordering North->South or South->North
, North_to_South                                                  &
                 ! indicates North->South ordering
, South_to_North ! indicates South->North ordering

PARAMETER                                                         &
  (North_to_South=1,South_to_North=2)

REAL ::                                                           &
  start_lat                                                       &
             ! Starting latitude
, start_long                                                      &
             ! Starting longitude
, del_lat                                                         &
             ! Latitude grid spacing
, del_long                                                        &
             ! Longitude grid spacing
, r_start_row                                                     &
              ! first row number
, r_end_row                                                       &
              ! last row number
, r_start_col                                                     &
              ! first column number
, r_end_col   ! last column number

! Function call
#if !defined(RECON)
INTEGER :: get_fld_type
#endif

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LLTORC'


!-----------------------------------------------------------------------

! Copy input arguments
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
north_lat=north_lat_in
south_lat=south_lat_in
west_long=west_long_in
east_long=east_long_in

IF (west_long  <   0) west_long=west_long+360
IF (east_long  <=  0) east_long=east_long+360

! Get information about grid type

! DEPENDS ON: gt_decode
CALL gt_decode(grid_code,                                         &
               model_type,content,coverage,domain,cyclic)

! Get information of field type for PARVARS variables

#if defined(RECON)
fld_type=rcf_get_fld_type(grid_code)
#else
! DEPENDS ON: get_fld_type
fld_type=get_fld_type(grid_code)
#endif

! Find start latitude and longitude of full grid, grid spacing and
! size of grid

IF (model_type  ==  gt_atmos) THEN

  start_lat=h_a_firstlat
  start_long=h_a_firstlong
  del_lat=h_a_nsspace
  del_long=h_a_ewspace
#if defined(RECON)
  IF (fld_type == fld_type_p) THEN
    nrows=recon_grid % glob_p_rows
    ncols=recon_grid % glob_p_row_length
  ELSE IF (fld_type == fld_type_u) THEN
    nrows=recon_grid % glob_u_rows
    ncols=recon_grid % glob_u_row_length
  ELSE IF (fld_type == fld_type_v) THEN
    nrows=recon_grid % glob_v_rows
    ncols=recon_grid % glob_v_row_length
  ELSE IF (fld_type == fld_type_r) THEN
    nrows=recon_grid % glob_r_rows
    ncols=recon_grid % glob_r_row_length
  END IF
#else
  nrows=decompDB(decomp_standard_atmos)%glsize(2,fld_type)
  ncols=decompDB(decomp_standard_atmos)%glsize(1,fld_type)
#endif
  row_ordering=South_to_North

  IF (h_global(atmos_im) == 'N') THEN ! LAM Configuration
    lam_model=.TRUE.
    pole_lat=h_a_polelat
    pole_long=h_a_polelong
  ELSE
    lam_model=.FALSE.
  END IF

END IF

! Ensure that DEL_LAT is always positive
! (ie. doesn't imply a direction)
IF (del_lat  <   0) del_lat=-del_lat

! This assumes the mass grid. The start latitude and longitude
! may be offset for velocity/U/V fields

IF (content  ==  gt_velocity) THEN  ! B grid U/V
  start_lat=start_lat+del_lat/2
  start_long=start_long+del_long/2
ELSE IF (content  ==  gt_U_C) THEN   ! C grid U
  IF (output_grid_stagger == FH_GridStagger_C) THEN
    start_long=start_long+del_long/2
  ELSE
    start_long=start_long-del_long/2
  END IF
  start_long=start_long+del_long/2
ELSE IF (content  ==  gt_V_C) THEN   ! C grid V
  IF (output_grid_stagger == FH_GridStagger_C) THEN
    start_lat=start_lat+del_lat/2
  ELSE
    start_lat=start_lat-del_lat/2
  END IF
END IF

! This assumes full domain. Now take account of zonal,meridional
! and scalar fields

IF (domain  ==  gt_zonal) THEN
  ncols=1
ELSE IF (domain  ==  gt_meridional) THEN
  nrows=1
ELSE IF (domain  ==  gt_scalar) THEN
  ncols=1
  nrows=1
ELSE IF (domain  ==  gt_ozone) THEN
  IF (zon_av_ozone) ncols=1
END IF

! This is the global sizes - this may be all we need

IF ((north_lat  ==  90) .AND. (south_lat  ==  -90) .AND.          &
    (west_long  ==  0) .AND. (east_long  ==  360)) THEN

  start_row=1
  end_row=nrows
  start_col=1
  end_col=ncols

ELSE ! Not a simple case

  ! If this is a LAM configuration we need to transform the latitudes
  ! and longitudes relative to the rotated pole

  IF (lam_model) THEN
    ! DEPENDS ON: lltoll
    CALL lltoll(                                                  &
        north_lat,south_lat,east_long,west_long,                  &
        pole_lat,pole_long)
  END IF

  ! Make sure that DEL_LAT has a sign consistent with the
  ! grid ordering. Further up we have ensured that DEL_LAT
  ! is positive. Now we change it negative if necessary.
  !               ( North->South => -ive
  !                 South->North => +ive )

  IF (row_ordering  ==  North_to_South) del_lat=-del_lat

  ! Calculate the start and end rows

  IF (row_ordering  ==  North_to_South) THEN

    r_start_row=1.0   + (north_lat-start_lat)/del_lat
    r_end_row  =1.999 + (start_lat-south_lat)/del_lat

  ELSE ! South->North

    r_start_row=1.0   + (south_lat-start_lat)/del_lat
    r_end_row  =1.999 + (north_lat-start_lat)/del_lat

  END IF

  start_row=MAX(1.0,r_start_row)
  end_row=MIN(REAL(nrows),r_end_row)

  IF (start_long  <   0) start_long=start_long+360.0

  IF ((REAL(start_long) + del_long*(ncols-1))  <=  360.0)         &
  THEN  ! If the total model domain doesn't cross
        ! the zero meridian

    r_start_col=1.0 +   (west_long-start_long)/del_long
    r_end_col=  1.999 + (east_long-start_long)/del_long

  ELSE ! model domain crosses the zero meridian

    IF (REAL(west_long)  <   start_long) THEN
      ! If the start is before the zero meridian

      r_start_col=1.0  + (west_long + 360.0 - start_long)/        &
                         del_long

    ELSE ! the start lies after the meridian

      r_start_col=1.0 + (REAL(west_long)-start_long)/             &
                        del_long

    END IF

    IF (REAL(east_long)  <   start_long) THEN
      ! If the end is after the zero meridian

      r_end_col=1.0 + (east_long + 360.0 - start_long)/           &
                      del_long

    ELSE ! the end lies before the zero meridian

      r_end_col=1.0 + (REAL(east_long) - start_long)/             &
                      del_long

    END IF

  END IF

  start_col=MIN(REAL(ncols),MAX(1.0,r_start_col))
  end_col=MIN(REAL(ncols),MAX(1.0,r_end_col))

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lltorc
