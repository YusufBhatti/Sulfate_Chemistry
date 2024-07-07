! *****************************COPYRIGHT*******************************
! 
! Copyright 2017-2018 University of Reading
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation 
! and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without 
! specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: NetCDF output
!
!  Purpose: Create the netCDF dimensions, coordinate variables,
!           auxiliary coordinate variables and bounds variables,
!           plus associated attributes needed for latitude/longitude 
!           coordinates. Update coordinates and cell_methods attributes
!           for STASH variable.

MODULE ncfile_write_horiz_dim_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE um_parvars,             ONLY: glsize
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,           ONLY: um_file_type
USE cderived_mod,           ONLY: elf
USE stash_model_mod,        ONLY: h_a_polelat,h_a_polelong
USE cstash_mod,             ONLY: idom_b,dompro
USE profilename_length_mod, ONLY: profilename_length
USE stash_array_mod,        ONLY: stlist
USE nlsizes_namelist_mod,   ONLY: global_row_length,river_row_length
USE stparam_mod,            ONLY: &
    st_diag_address,st_gridpoint_code,st_weight_code,st_domain_code, &
    st_west_code,st_east_code,st_north_code,st_south_code, &
    zonal_mean_base,merid_mean_base,field_mean_base, &
    global_mean_base,extract_base,vert_mean_base,block_size, &
    stash_land_mask_code,stash_sea_mask_code, &
    stash_weight_area_code,stash_weight_volume_code,stash_weight_mass_code, &
    st_domain_n_hemisphere, st_domain_s_hemisphere, st_domain_30_to_90_N, &
    st_domain_30_to_90_S, st_domain_0_to_30_N, st_domain_0_to_30_S, &
    st_domain_30_S_to_30_N, st_domain_whole_degrees, st_domain_gridpoints
USE cppxref_mod,            ONLY: &
    ppx_atm_tall,ppx_atm_tland,ppx_atm_tsea,ppx_atm_tzonal,ppx_atm_tmerid, &
    ppx_atm_uall,ppx_atm_uland,ppx_atm_usea,ppx_atm_uzonal,ppx_atm_umerid, &
    ppx_atm_compressed,ppx_atm_ozone,ppx_atm_cuall,ppx_atm_cvall,ppx_atm_river
USE nc_dimension_id_mod,    ONLY: &
    nc_horiz_id_theta,nc_horiz_id_u,nc_horiz_id_v, &
    nc_horiz_id_uv,nc_horiz_id_river
USE umnetcdf_mod,           ONLY: &
    nc_create_dim,nc_create_var,nc_put_att,check_cf_name,nc_dim_type_real

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_HORIZ_DIM_MOD'

CONTAINS

SUBROUTINE ncfile_write_horiz_dim(um_file,grid_type_code,ix, &
                                  coordinates,cell_methods, &
                                  dimid_lon,dimid_lat)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER           , INTENT(IN)    :: grid_type_code ! STASH gridtype
INTEGER           , INTENT(IN)    :: ix             ! STASH item index
CHARACTER(LEN=*)  , INTENT(INOUT) :: coordinates    ! CF coordinates attribute
CHARACTER(LEN=*)  , INTENT(INOUT) :: cell_methods   ! CF cell_methods attribute
INTEGER           , INTENT(OUT)   :: dimid_lon,dimid_lat ! NetCDF dimensions IDs

!
!  Local variables
!
INTEGER :: n_rows_out              ! No of rows used for a diagnostic
INTEGER :: n_cols_out              ! No of cols used
INTEGER :: srow_in,srow_out        ! North, South, East and West
INTEGER :: nrow_in,nrow_out        ! subdomain limits in the horizontal sense
INTEGER :: wcol_in,wcol_out        ! corresponding to the subarea before being
INTEGER :: ecol_in,ecol_out        ! processed (IN) and after processing (OUT)
INTEGER :: domain_opt              ! domain option
INTEGER :: weight_opt              ! weighting option
INTEGER :: mean_opt                ! meaning option
INTEGER :: mask_opt                ! masking option (land+sea/land/sea)
INTEGER :: idiag                   ! STASH diagnostic address
INTEGER :: fld_type                ! field type: u,v or t location on C grid
INTEGER :: dimid_array(1)          ! used to convert dimid value into 1-D array
INTEGER :: dimid_truell(2)         ! dimension IDs for true lat/lon 
                                   ! auxiliary coordinates
INTEGER :: dimid_bnds              ! dimension ID for bounds variable
INTEGER :: dimid_lon_bnds(2)       ! dimension IDs for longitude bounds variable
INTEGER :: dimid_lat_bnds(2)       ! dimension IDs for latitude bounds variable
INTEGER :: varid_lon               ! longitude coordinate variable ID
INTEGER :: varid_lat               ! latitude coordinate variable ID
INTEGER :: varid_true_lon          ! true longitude coordinate variable ID
INTEGER :: varid_true_lat          ! true latitude coordinate variable ID
INTEGER :: varid_bnds              ! lat/lon bounds variable ID
INTEGER :: num_long_pts            ! Number of longitude points for the grid
                                   ! type associated with this STASH item

LOGICAL :: lexists_lon,lexists_lat ! TRUE if lat/lon dimension already exists
LOGICAL :: lcyclic                 ! TRUE if cyclic EW BCs
LOGICAL :: lrotgrid                ! TRUE if rotated grid
LOGICAL :: ltruell                 ! TRUE if true lat/lon values to be written
LOGICAL :: l_unique_lon            ! TRUE if unique longitude dimension name 
                                   ! needed for domain options 2-8

CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'NCFILE_WRITE_HORIZ_DIM'
CHARACTER(LEN=profilename_length) :: domname ! Domain profile name
                                   ! lat/lon dimension names
CHARACTER(LEN=nf90_max_name)      :: dimname_lon,dimname_lat
                                   ! true lat/lon dimension names
CHARACTER(LEN=nf90_max_name)      :: dimname_true_lon,dimname_true_lat

! External functions:
INTEGER :: get_fld_type

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get domain name for this STASH diagnostic
idiag = stlist(st_diag_address,ix)
domname = dompro(idom_b(idiag))

! Make sure domname has no characters disallowed by the CF standard
CALL check_cf_name(domname)
lcyclic = .NOT. elf
lrotgrid = .NOT. (h_a_polelat == 90.0 .AND. h_a_polelong == 0.0)

! Get domain, meaning and weighting options
domain_opt = stlist(st_domain_code,ix)
mean_opt = (stlist(st_gridpoint_code,ix)/block_size)*block_size
mask_opt = MOD(stlist(st_gridpoint_code,ix),block_size)
weight_opt = stlist(st_weight_code,ix)
IF (printstatus >= prdiag) THEN
  CALL umPrint(TRIM(RoutineName)//' variables:',src=RoutineName)
  WRITE(umMessage,'(A,3(I3))') &
             'domain_opt,mean_opt,weight_opt = ', &
              domain_opt,mean_opt,weight_opt
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Find the horizontal dimensions for the output grid
wcol_in = stlist(st_west_code,ix)    ! Input subdomain limits
ecol_in = stlist(st_east_code,ix)
nrow_in = stlist(st_north_code,ix)
srow_in = stlist(st_south_code,ix)

SELECT CASE (mean_opt)
CASE (extract_base,vert_mean_base)
  wcol_out = wcol_in
  ecol_out = ecol_in
  nrow_out = nrow_in
  srow_out = srow_in
CASE (zonal_mean_base)
  wcol_out = 1
  ecol_out = 1
  nrow_out = nrow_in
  srow_out = srow_in
CASE (merid_mean_base)
  wcol_out = wcol_in
  ecol_out = ecol_in
  nrow_out = 1
  srow_out = 1
CASE (field_mean_base,global_mean_base)
  wcol_out = 1
  ecol_out = 1
  nrow_out = 1
  srow_out = 1
END SELECT

! Get field type, ie u or v or p location in Arakawa C grid staggering
! DEPENDS ON: get_fld_type
fld_type = get_fld_type(grid_type_code)

! Adjust easternmost column if field wraps east-west
IF ((wcol_in > ecol_in) .AND. lcyclic) &
  ecol_in = ecol_in + glsize(1,fld_type)

IF ((wcol_out > ecol_out) .AND. lcyclic) &
  ecol_out = ecol_out + glsize(1,fld_type)

! Calculate number of rows and columns of output field
n_rows_out = nrow_out - srow_out + 1
n_cols_out = ecol_out - wcol_out + 1

IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(2(A,I6))')'wcol_in  = ',wcol_in,' ecol_in  = ',ecol_in
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(2(A,I6))')'nrow_in  = ',nrow_in,' srow_in  = ',srow_in
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(2(A,I6))')'wcol_out = ',wcol_out,' ecol_out = ',ecol_out
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(2(A,I6))')'nrow_out = ',nrow_out,' srow_out = ',srow_out
  CALL umPrint(umMessage,src=RoutineName)
END IF

! IF ltruell is .TRUE. then true lat/lon grid is also included in NetCDF file
ltruell = lrotgrid .AND. n_cols_out > 1 .AND. n_rows_out > 1

! Set horizontal grid dimension names
IF (lrotgrid) THEN
  dimname_lon = 'grid_longitude_'
  dimname_lat = 'grid_latitude_'
ELSE
  dimname_lon = 'longitude_'
  dimname_lat = 'latitude_'
END IF

! Number of longitude points for grid type unless river routing grid
num_long_pts = global_row_length

! Make unique latitude/longitude dimension names for grid type
SELECT CASE (grid_type_code)
CASE (ppx_atm_tall, &
      ppx_atm_tland, &
      ppx_atm_tsea, &
      ppx_atm_tzonal, &
      ppx_atm_tmerid, &
      ppx_atm_compressed, &
      ppx_atm_ozone)
  dimname_lon = TRIM(dimname_lon)//nc_horiz_id_theta
  dimname_lat = TRIM(dimname_lat)//nc_horiz_id_theta
CASE (ppx_atm_cuall)
  dimname_lon = TRIM(dimname_lon)//nc_horiz_id_u
  dimname_lat = TRIM(dimname_lat)//nc_horiz_id_u
CASE (ppx_atm_cvall)
  dimname_lon = TRIM(dimname_lon)//nc_horiz_id_v
  dimname_lat = TRIM(dimname_lat)//nc_horiz_id_v
CASE (ppx_atm_uall, &
      ppx_atm_uland, &
      ppx_atm_usea, &
      ppx_atm_uzonal, &
      ppx_atm_umerid)
  dimname_lon = TRIM(dimname_lon)//nc_horiz_id_uv
  dimname_lat = TRIM(dimname_lat)//nc_horiz_id_uv
CASE (ppx_atm_river)
  dimname_lon = TRIM(dimname_lon)//nc_horiz_id_river
  dimname_lat = TRIM(dimname_lat)//nc_horiz_id_river
  num_long_pts = river_row_length
END SELECT

! For domain options 2-8 longitude coordinate should be unchanged but for UKV
! grids this is not the case, so check whether unique dimension name needed
l_unique_lon = wcol_in /= 1 .OR. ecol_in /= num_long_pts

! Make unique latitude/longitude dimension name for selected domain option
SELECT CASE (domain_opt)
CASE (st_domain_n_hemisphere)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_90N_0'
  dimname_lat = TRIM(dimname_lat)//'_90N_0'
CASE (st_domain_s_hemisphere)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_0_90S'
  dimname_lat = TRIM(dimname_lat)//'_0_90S'
CASE (st_domain_30_to_90_N)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_90N_30N'
  dimname_lat = TRIM(dimname_lat)//'_90N_30N'
CASE (st_domain_30_to_90_S)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_30S_90S'
  dimname_lat = TRIM(dimname_lat)//'_30S_90S'
CASE (st_domain_0_to_30_N)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_30N_0'
  dimname_lat = TRIM(dimname_lat)//'_30N_0'
CASE (st_domain_0_to_30_S)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_0_30S'
  dimname_lat = TRIM(dimname_lat)//'_0_30S'
CASE (st_domain_30_S_to_30_N)
  IF (l_unique_lon) dimname_lon = TRIM(dimname_lon)//'_30N_30S'
  dimname_lat = TRIM(dimname_lat)//'_30N_30S'
CASE (st_domain_whole_degrees, st_domain_gridpoints)
  ! Use domname to make latitude/longitude dimension names  
  ! unique for a particular domain profile
  dimname_lon = TRIM(dimname_lon)//'_'//TRIM(domname)
  dimname_lat = TRIM(dimname_lat)//'_'//TRIM(domname)
END SELECT

! Make unique latitude/longitude dimension names for selected meaning option
SELECT CASE (mean_opt)
CASE (zonal_mean_base)
  dimname_lon = TRIM(dimname_lon)//'_zonal_mean'
CASE (merid_mean_base)
  dimname_lat = TRIM(dimname_lat)//'_meridional_mean'
CASE (field_mean_base)
  dimname_lon = TRIM(dimname_lon)//'_field_mean'
  dimname_lat = TRIM(dimname_lat)//'_field_mean'
CASE (global_mean_base)
  dimname_lon = TRIM(dimname_lon)//'_global_mean'
  dimname_lat = TRIM(dimname_lat)//'_global_mean'
END SELECT

IF (ltruell) THEN
  ! If ltruell is TRUE then lrotgrid must also be TRUE so remove 'grid_' 
  ! from dimension names to get true latitude/logitude dimension names
  dimname_true_lon = dimname_lon(6:)
  dimname_true_lat = dimname_lat(6:)
END IF

IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(4(A))') &
        'dimname_lon = ',TRIM(dimname_lon), &
        ' : dimname_lat = ',TRIM(dimname_lat)
  CALL umPrint(umMessage,src=RoutineName)
  IF (ltruell) THEN
    WRITE(umMessage,'(4(A))') &
          'dimname_true_lon = ',TRIM(dimname_true_lon), &
          ' : dimname_true_lat = ',TRIM(dimname_true_lat)
    CALL umPrint(umMessage,src=RoutineName)
  END IF
END IF

! Create NetCDF dimensions if they don't already exist
CALL nc_create_dim(um_file,dimname_lon,n_cols_out,dimid_lon,lexists_lon)
CALL nc_create_dim(um_file,dimname_lat,n_rows_out,dimid_lat,lexists_lat)

IF (ltruell) THEN
  ! Add true latitude/longitude dimensions to coordinates attribute
  IF (TRIM(coordinates) == "") THEN
    coordinates = TRIM(dimname_true_lon)
  ELSE
    coordinates = TRIM(coordinates)//' '//TRIM(dimname_true_lon)
  END IF
  coordinates = TRIM(coordinates)//' '//TRIM(dimname_true_lat)

  ! Associate longitude dimension id with latitude dimension id
  ! Needed when writing out true latitude/longitude data
  um_file % nc_meta % dimid_lon(dimid_lat) = dimid_lon
END IF

! Add any meaning/weighting options to the cell_methods attribute
! for the horizontal dimensions
SELECT CASE (mean_opt)
CASE (zonal_mean_base) ! Zonal mean
  cell_methods = TRIM(dimname_lon)//': mean'
CASE (merid_mean_base) ! Meridional mean
  cell_methods = TRIM(dimname_lat)//': mean'
CASE (field_mean_base) ! Horizontal mean
  cell_methods = 'area: mean'
CASE (global_mean_base) ! Global mean
  cell_methods = TRIM(dimname_lon)//': '//TRIM(dimname_lat)//':'
END SELECT

IF (mean_opt == zonal_mean_base .OR. mean_opt == merid_mean_base .OR. &
    mean_opt == field_mean_base) THEN

  SELECT CASE (mask_opt)
  CASE (stash_land_mask_code)
    cell_methods = TRIM(cell_methods)//' where land'
  CASE (stash_sea_mask_code)
    cell_methods = TRIM(cell_methods)//' where sea'
  END SELECT

  SELECT CASE (weight_opt)
  CASE (stash_weight_area_code)
    cell_methods = TRIM(cell_methods)//' (area-weighted)'
  CASE (stash_weight_volume_code)
    cell_methods = TRIM(cell_methods)//' (volume-weighted)'
  CASE (stash_weight_mass_code)
    cell_methods = TRIM(cell_methods)//' (mass-weighted)'
  END SELECT
END IF
IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(A,A)')'cell_methods = ',TRIM(cell_methods)
  CALL umPrint(umMessage,src=RoutineName)
END IF

IF (lexists_lon .AND. lexists_lat) THEN
  ! Both dimensions already exist so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Create bounds dimension if it doesn't already exist
CALL nc_create_dim(um_file,'bounds2',2,dimid_bnds)

IF (ltruell) THEN
  ! True latitude/longitude variables have 
  ! latitude/longitude dimensions for coordinates
  dimid_truell(1) = dimid_lon
  dimid_truell(2) = dimid_lat
END IF

IF (.NOT. lexists_lon) THEN

  ! Create longitude coordinate variable
  dimid_array(1) = dimid_lon
  CALL nc_create_var(um_file,dimname_lon,nc_dim_type_real,varid_lon,dimid_array)

  IF (ltruell) THEN
    ! Create the true longitude auxiliary coordinate variable
    CALL nc_create_var(um_file,dimname_true_lon,nc_dim_type_real, &
                       varid_true_lon,dimid_truell)
  END IF

  ! Bounds variable dimension IDs for longitude coordinate
  dimid_lon_bnds(1) = dimid_bnds
  dimid_lon_bnds(2) = dimid_lon

  ! Create bounds variable for the longitude coordinate
  CALL nc_create_var(um_file,TRIM(dimname_lon)//'_bounds',nc_dim_type_real, &
                     varid_bnds,dimid_lon_bnds)

  ! Create longitude coordinate variable attributes
  IF (lrotgrid) THEN
    CALL nc_put_att(um_file,varid_lon,'standard_name','grid_longitude')
    CALL nc_put_att(um_file,varid_lon,'long_name', &
                                      'longitude in rotated pole grid')
    CALL nc_put_att(um_file,varid_lon,'units','degrees')
  ELSE
    CALL nc_put_att(um_file,varid_lon,'standard_name','longitude')
    CALL nc_put_att(um_file,varid_lon,'long_name','longitude')
    CALL nc_put_att(um_file,varid_lon,'units','degrees_east')
  END IF
  IF (ltruell) THEN
    CALL nc_put_att(um_file,varid_true_lon,'standard_name','longitude')
    CALL nc_put_att(um_file,varid_true_lon,'long_name','longitude')
    CALL nc_put_att(um_file,varid_true_lon,'units','degrees_east')
  END IF
  CALL nc_put_att(um_file,varid_lon,'axis','X')
  CALL nc_put_att(um_file,varid_lon,'bounds',TRIM(dimname_lon)//'_bounds')

  ! Save stash index value for later use
  um_file % nc_meta % stash_index_dim(dimid_lon) = ix
END IF

IF (.NOT. lexists_lat) THEN

  ! Create latitude coordinate variable
  dimid_array(1) = dimid_lat
  CALL nc_create_var(um_file,dimname_lat,nc_dim_type_real,varid_lat,dimid_array)

  IF (ltruell) THEN
    ! Create the true latitude auxiliary coordinate variable
    CALL nc_create_var(um_file,dimname_true_lat,nc_dim_type_real, &
                       varid_true_lat,dimid_truell)
  END IF

  ! Bounds variable dimension IDs for latitude coordinate
  dimid_lat_bnds(1) = dimid_bnds
  dimid_lat_bnds(2) = dimid_lat

  ! Create bounds variable for the latitude coordinate
  CALL nc_create_var(um_file,TRIM(dimname_lat)//'_bounds',nc_dim_type_real, &
                     varid_bnds,dimid_lat_bnds)

  ! Create latitude coordinate variable attributes
  IF (lrotgrid) THEN
    CALL nc_put_att(um_file,varid_lat,'standard_name','grid_latitude')
    CALL nc_put_att(um_file,varid_lat,'long_name', &
                                      'latitude in rotated pole grid')
    CALL nc_put_att(um_file,varid_lat,'units','degrees')
  ELSE
    CALL nc_put_att(um_file,varid_lat,'standard_name','latitude')
    CALL nc_put_att(um_file,varid_lat,'long_name','latitude')
    CALL nc_put_att(um_file,varid_lat,'units','degrees_north')
  END IF
  IF (ltruell) THEN
    CALL nc_put_att(um_file,varid_true_lat,'standard_name','latitude')
    CALL nc_put_att(um_file,varid_true_lat,'long_name','latitude')
    CALL nc_put_att(um_file,varid_true_lat,'units','degrees_north')
  END IF
  CALL nc_put_att(um_file,varid_lat,'axis','Y')
  CALL nc_put_att(um_file,varid_lat,'bounds',TRIM(dimname_lat)//'_bounds')

  ! Save stash index value for later use
  um_file % nc_meta % stash_index_dim(dimid_lat) = ix
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_horiz_dim

END MODULE ncfile_write_horiz_dim_mod
