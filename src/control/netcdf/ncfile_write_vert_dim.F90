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
!           plus associated attributes needed for vertical coordinates. 
!           Update coordinates and cell_methods attributes for STASH variable.


MODULE ncfile_write_vert_dim_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,           ONLY: um_file_type
USE cstash_mod,             ONLY: idom_b,dompro
USE profilename_length_mod, ONLY: profilename_length
USE stash_array_mod,        ONLY: stlist,stash_levels
USE stparam_mod,            ONLY:                                              &
    st_diag_address,st_output_bottom,st_output_top,                            &
    st_gridpoint_code,st_weight_code,                                          &
    global_mean_base,vert_mean_base,block_size,                                &
    stash_weight_area_code,stash_weight_volume_code,stash_weight_mass_code,    &
    st_levels_model_rho, st_levels_model_theta, st_levels_pressure,            &
    st_levels_geometric_height, st_levels_single, st_levels_deep_soil,         &
    st_levels_potential_temp, st_levels_potential_vort, st_levels_cloud_thresh
USE nc_dimension_id_mod,    ONLY: &
    nc_vert_id_model_level, &
    nc_vert_id_eta_rho,nc_vert_id_zsea_rho,nc_vert_id_C_rho, &
    nc_vert_id_eta_theta,nc_vert_id_zsea_theta,nc_vert_id_C_theta
USE umnetcdf_mod,           ONLY: &
    nc_create_dim,nc_create_var,nc_put_att,check_cf_name, &
    nc_dim_type_real,nc_type_int,nc_max_attlen

IMPLICIT NONE

! Parameters for state of the 'positive' dimension attribute
INTEGER, PARAMETER, PRIVATE :: positive_unset = 0
INTEGER, PARAMETER, PRIVATE :: positive_up = 1
INTEGER, PARAMETER, PRIVATE :: positive_down = 2

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_VERT_DIM_MOD'

CONTAINS

SUBROUTINE ncfile_write_vert_dim(um_file,level_code,ix,cf_extra, &
                                 coordinates,cell_methods,dimid)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN) :: level_code               ! STASH level type code
INTEGER, INTENT(IN) :: ix                       ! STASH item index
CHARACTER(LEN=*), INTENT(IN   ) :: cf_extra     ! Contains height info
CHARACTER(LEN=*), INTENT(INOUT) :: coordinates  ! CF coordinates attribute
CHARACTER(LEN=*), INTENT(INOUT) :: cell_methods ! CF cell_methods attribute
INTEGER, INTENT(OUT) :: dimid                   ! NetCDF dimension ID

!
!  Local variables
!
INTEGER :: idiag              ! STASH diagnostic address
INTEGER :: nlevs              ! number of vertical levels
INTEGER :: levbot             ! < 0 points to levels list -levbot
                              ! > 0 start level for model level output
INTEGER :: levtop             ! levbot > 0 last level for model level output
INTEGER :: dimid_array(1)     ! used to convert dimid value into 1-D array
INTEGER :: dimid_bnds         ! dimension ID for bounds variable
INTEGER :: dimid_vert_bnds(2) ! dimension IDs for vertical bounds variable
INTEGER :: varid              ! variable ID for vertical coordinates
INTEGER :: varid_bnds         ! variable ID for bounds variable
INTEGER :: datatype           ! data type for main vertical coordinate
INTEGER :: mean_opt           ! meaning option
INTEGER :: weight_opt         ! weighting option
INTEGER :: positive           ! direction of positive for dimension
INTEGER :: idx                ! Position of . in dimname

LOGICAL :: lexists            ! TRUE if vertical dimension already exists
LOGICAL :: lusebounds         ! TRUE if vertical dimension has a 
                              ! associated bounds variable

CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'NCFILE_WRITE_VERT_DIM'
CHARACTER(LEN=profilename_length) :: domname      ! Domain profile name
CHARACTER(LEN=nf90_max_name)      :: dimname      ! Dimension name
CHARACTER(LEN=nf90_max_name)      :: varname1     ! Auxiliary coordinate
CHARACTER(LEN=nf90_max_name)      :: varname2     ! variable names
CHARACTER(LEN=nc_max_attlen)      :: longname     ! Long name for dimension
CHARACTER(LEN=nc_max_attlen)      :: stdname      ! Standard name for dimension
CHARACTER(LEN=nc_max_attlen)      :: units        ! Units name for dimension

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dimid = 0

! Get domain name for this STASH diagnostic
idiag = stlist(st_diag_address,ix)
domname = dompro(idom_b(idiag))

! Make sure domname has no characters disallowed by the CF standard
CALL check_cf_name(domname)

! Get levels information
levbot = stlist(st_output_bottom,ix)
levtop = stlist(st_output_top,ix)

IF (printstatus >= prdiag) THEN
  CALL umPrint(TRIM(RoutineName)//' variables:',src=RoutineName)
  WRITE(umMessage,'(A,2(I4))') &
       'level_code = ',level_code
  WRITE(umMessage,'(A,I3,I3)') &
       'levbot,levtop = ',levbot,levtop
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Get meaning and weighting options
mean_opt = (stlist(st_gridpoint_code,ix)/block_size)*block_size
weight_opt = stlist(st_weight_code,ix)

datatype = nc_dim_type_real ! default data type
lusebounds = .FALSE.

! Get number of levels for this vertical dimension variable
IF (level_code == st_levels_single) THEN
  IF (LEN(cf_extra) > 10 .AND. cf_extra(1:7) == 'height=') THEN
    ! Single height level, height value specified in STASH to CF metadata
    nlevs = 1
  ELSE
    dimid = -1
    ! Single level with no value. No dimension variable needed, so finish
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE IF (mean_opt == vert_mean_base .OR. mean_opt == global_mean_base) THEN
  ! vertical or global mean
  nlevs = 1
  ! Always have bounds for means
  lusebounds = .TRUE.
ELSE IF (levbot < 0) THEN                   ! levels defined using STASH list
  nlevs = stash_levels(1,-levbot)
ELSE                                        ! model level range
  nlevs = levtop-levbot+1
END IF

! Default values of vertical dimension attributes
stdname = ' '
longname = ' '
units = ' '
positive = positive_unset
dimname=domname

! Depending on type of level specify various vertical dimension attributes
SELECT CASE (level_code)
CASE (st_levels_model_rho)
  stdname = "atmosphere_hybrid_height_coordinate"
  longname = "eta value of rho levels"
  dimname=TRIM(domname)//'_'//nc_vert_id_eta_rho
  positive = positive_up
  ! Add auxiliary vertical coordinate variables to coordinates attribute
  IF (TRIM(coordinates) == "") THEN
    coordinates = TRIM(domname)//'_'//nc_vert_id_zsea_rho
  ELSE
    coordinates = TRIM(coordinates)//' '// &
                  TRIM(domname)//'_'//nc_vert_id_zsea_rho
  END IF
  coordinates = TRIM(coordinates)//' '// &
                TRIM(domname)//'_'//nc_vert_id_C_rho
  coordinates = TRIM(coordinates)//' '// &
                TRIM(domname)//'_'//nc_vert_id_model_level
  lusebounds = .TRUE.
CASE (st_levels_model_theta)
  stdname = "atmosphere_hybrid_height_coordinate"
  longname = "eta value of theta levels"
  dimname=TRIM(domname)//'_'//nc_vert_id_eta_theta
  positive = positive_up
  ! Add auxiliary vertical coordinate variables to coordinates attribute
  IF (TRIM(coordinates) == "") THEN
    coordinates = TRIM(domname)//'_'//nc_vert_id_zsea_theta
  ELSE
    coordinates = TRIM(coordinates)//' '// &
                  TRIM(domname)//'_'//nc_vert_id_zsea_theta
  END IF
  coordinates = TRIM(coordinates)//' '// &
                TRIM(domname)//'_'//nc_vert_id_C_theta
  coordinates = TRIM(coordinates)//' '// &
                TRIM(domname)//'_'//nc_vert_id_model_level
  lusebounds = .TRUE.
CASE (st_levels_deep_soil)
  stdname = "depth"
  longname = "soil layer mid-point depth"
  units = "m"
  positive = positive_down
  ! Add auxiliary vertical coordinate variable to coordinates attribute
  IF (TRIM(coordinates) == "") THEN
    coordinates = TRIM(domname)//'_'//nc_vert_id_model_level
  ELSE
    coordinates = TRIM(coordinates)//' '// &
                  TRIM(domname)//'_'//nc_vert_id_model_level
  END IF
  lusebounds = .TRUE.
CASE (st_levels_pressure)
  stdname = "air_pressure"
  longname = "pressure levels"
  units = "hPa"
  positive = positive_down
CASE (st_levels_geometric_height)
  stdname = "altitude"
  longname = "geometric height levels"
  units = "m"
  positive = positive_up
CASE (st_levels_potential_temp)
  stdname = "air_potential_temperature"
  longname = "constant theta surfaces"
  units = "K"
  positive = positive_down
CASE (st_levels_potential_vort)
  stdname = "ertel_potential_vorticity"
  longname = "potential vorticity levels"
  units = "K m2 kg-1 s-1"
  positive = positive_up
CASE (st_levels_cloud_thresh)
  longname = "cloud threshold levels"
CASE (st_levels_single)
  ! Only single level with a dimension defined is height
  stdname = "height"
  longname = cf_extra(8:LEN_TRIM(cf_extra))//' height'
  dimname='height_'//cf_extra(8:LEN_TRIM(cf_extra))
  ! If dimname contains a decimal point it is an invalid CF name
  ! Replace . with _
  idx = INDEX(dimname,".")
  IF (idx > 0) THEN
    dimname(idx:idx) = "_"
  END IF
  units = "m"
CASE DEFAULT
  WRITE(umMessage,'(A,I3)') 'Unknown level_code: ',level_code
  CALL umPrint(umMessage,src=RoutineName)
END SELECT

! Get cell_methods attribute information for the vertical dimension
IF (mean_opt == vert_mean_base .OR. mean_opt == global_mean_base) THEN
  ! Vertical mean
  IF (cell_methods == ' ') THEN
    cell_methods = TRIM(dimname)//': mean'
  ELSE
    cell_methods = TRIM(cell_methods)//' '//TRIM(dimname)//': mean'
  END IF

  ! Add any weighting options to cell_methods
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
  WRITE(umMessage,'(8(A))') &
            'dimname = ',TRIM(dimname), &
            ' : stdname = ',TRIM(stdname), &
            ' : longname = ',TRIM(longname), &
            ' : units = ',TRIM(units)
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(3(A),I2,A,I3)') &
            'cell_methods = ',TRIM(cell_methods), &
            ' : positive = ',positive, &
            ' : nlevs = ',nlevs
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Create or read NetCDF dimension
CALL nc_create_dim(um_file,dimname,nlevs,dimid,lexists)

IF (lexists) THEN
  ! Vertical dimension already exists so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

dimid_array(1) = dimid

! Create NetCDF dimension variables
CALL nc_create_var(um_file,dimname,datatype,varid,dimid_array)

! Save stash index value for later use
um_file % nc_meta % stash_index_dim(dimid) = ix

! Create or read bounds dimension if needed
IF (lusebounds) THEN
  CALL nc_create_dim(um_file,'bounds2',2,dimid_bnds)
  dimid_vert_bnds(1) = dimid_bnds
  dimid_vert_bnds(2) = dimid
END IF

! Create NetCDF dimension variable attributes
IF (stdname(1:1) /= ' ') THEN
  CALL nc_put_att(um_file,varid,'standard_name',stdname)
END IF
IF (longname(1:1) /= ' ') THEN
  CALL nc_put_att(um_file,varid,'long_name',longname)
END IF
IF (units(1:1) /= ' ') THEN
  CALL nc_put_att(um_file,varid,'units',units)
END IF
CALL nc_put_att(um_file,varid,'axis','Z')
IF (positive == positive_up) THEN
  CALL nc_put_att(um_file,varid,'positive','up')
ELSE IF (positive == positive_down) THEN
  CALL nc_put_att(um_file,varid,'positive','down')
END IF
IF (lusebounds) THEN
  CALL nc_put_att(um_file,varid,'bounds',TRIM(dimname)//'_bounds')
  CALL nc_create_var(um_file,TRIM(dimname)//'_bounds',datatype, &
                     varid_bnds,dimid_vert_bnds)
END IF

! Extra variables/attributes for model levels and soil levels
SELECT CASE (level_code)
CASE (st_levels_model_rho, st_levels_model_theta)

  IF (level_code == st_levels_model_rho) THEN
    varname1 = TRIM(domname)//'_'//nc_vert_id_zsea_rho
    varname2 = TRIM(domname)//'_'//nc_vert_id_C_rho
  ELSE IF (level_code == st_levels_model_theta) THEN
    varname1 = TRIM(domname)//'_'//nc_vert_id_zsea_theta
    varname2 = TRIM(domname)//'_'//nc_vert_id_C_theta
  END IF

  ! Add zsea variable
  CALL nc_create_var(um_file,varname1,nc_dim_type_real,varid,dimid_array)
  CALL nc_put_att(um_file,varid,'standard_name', &
                                'height_above_reference_ellipsoid')
  CALL nc_put_att(um_file,varid,'long_name','Height above mean sea level')
  CALL nc_put_att(um_file,varid,'units','m')
  CALL nc_put_att(um_file,varid,'positive','up')
  IF (lusebounds) THEN
    CALL nc_put_att(um_file,varid,'bounds',TRIM(varname1)//'_bounds')
    CALL nc_create_var(um_file,TRIM(varname1)//'_bounds',nc_dim_type_real, &
                       varid_bnds,dimid_vert_bnds)
  END IF

  ! Add C variable
  CALL nc_create_var(um_file,varname2,nc_dim_type_real,varid,dimid_array)
  CALL nc_put_att(um_file,varid,'long_name','Fraction of orographic height')
  CALL nc_put_att(um_file,varid,'units','1')
  IF (lusebounds) THEN
    CALL nc_put_att(um_file,varid,'bounds',TRIM(varname2)//'_bounds')
    CALL nc_create_var(um_file,TRIM(varname2)//'_bounds',nc_dim_type_real, &
                       varid_bnds,dimid_vert_bnds)
  END IF

  ! Add model level number as auxiliary coordinate variable
  varname1 = TRIM(domname)//'_'//nc_vert_id_model_level
  CALL nc_create_var(um_file,varname1,nc_type_int,varid,dimid_array)
  CALL nc_put_att(um_file,varid,'standard_name','model_level_number')
  IF (level_code == st_levels_model_rho) THEN
    longname = 'model rho levels (Charney-Phillips grid)'
  ELSE IF (level_code == st_levels_model_theta) THEN
    longname = 'model theta levels (Charney-Phillips grid)'
  END IF
  CALL nc_put_att(um_file,varid,'long_name',longname)
  CALL nc_put_att(um_file,varid,'units','1')
  CALL nc_put_att(um_file,varid,'positive','up')
  IF (mean_opt == vert_mean_base .OR. mean_opt == global_mean_base) THEN
    CALL nc_put_att(um_file,varid,'bounds',TRIM(varname1)//'_bounds')
    CALL nc_create_var(um_file,TRIM(varname1)//'_bounds',nc_type_int, &
                       varid_bnds,dimid_vert_bnds)
  END IF

CASE (st_levels_deep_soil)

  ! Add soil model level number as auxiliary coordinate variable
  varname1 = TRIM(domname)//'_'//nc_vert_id_model_level
  CALL nc_create_var(um_file,varname1,nc_type_int,varid,dimid_array)
  CALL nc_put_att(um_file,varid,'standard_name','model_level_number')
  CALL nc_put_att(um_file,varid,'long_name','soil model level number')
  CALL nc_put_att(um_file,varid,'units','1')
  CALL nc_put_att(um_file,varid,'positive','down')
  IF (mean_opt == vert_mean_base .OR. mean_opt == global_mean_base) THEN
    CALL nc_put_att(um_file,varid,'bounds',TRIM(varname1)//'_bounds')
    CALL nc_create_var(um_file,TRIM(varname1)//'_bounds',nc_type_int, &
                       varid_bnds,dimid_vert_bnds)
  END IF

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_vert_dim

END MODULE ncfile_write_vert_dim_mod
