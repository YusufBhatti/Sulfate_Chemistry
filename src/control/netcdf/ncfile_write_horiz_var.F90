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
!  Purpose: Writes data for coordinate variables, auxiliary coordinate variables
!           and bounds variables created in ncfile_write_horiz_dim

MODULE ncfile_write_horiz_var_mod

USE netcdf,                ONLY: nf90_max_name ! External netCDF library module
USE yomhook,               ONLY: lhook,dr_hook
USE parkind1,              ONLY: jprb,jpim
USE umprintmgr,            ONLY: umprint,ummessage,printstatus,prdiag
USE ereport_mod,           ONLY: ereport
USE file_manager,          ONLY: um_file_type
USE cderived_mod,          ONLY: elf
USE model_domain_mod,      ONLY: l_regular
USE stash_array_mod,       ONLY: stlist
USE latlon_eq_rotation_mod,ONLY: rotate_eq_to_latlon
USE nlsizes_namelist_mod,  ONLY: &
    global_row_length,global_rows,river_row_length,river_rows
USE dump_headers_mod,      ONLY: a_rowdepc,a_coldepc
USE stash_model_mod,       ONLY: &
    h_a_ewspace,h_a_nsspace,h_a_firstlat,h_a_firstlong
USE stparam_mod,           ONLY: &
    st_west_code,st_east_code,st_north_code,st_south_code, &
    zonal_mean_base,merid_mean_base,field_mean_base,global_mean_base, &
    st_gridpoint_code,block_size
USE atm_d1_indices_mod,    ONLY: &
    jlambda_input_p,jlambda_input_u,jphi_input_p,jphi_input_v
USE nc_dimension_id_mod,   ONLY: &
    nc_horiz_id_theta,nc_horiz_id_u,nc_horiz_id_v, &
    nc_horiz_id_uv,nc_horiz_id_river
USE umnetcdf_mod,          ONLY: &
    nc_inq_ndims,nc_inq_dim,nc_get_varid,nc_get_var,nc_get_att,nc_put_var

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_HORIZ_VAR_MOD'

PRIVATE ncfile_write_truell

CONTAINS

SUBROUTINE ncfile_write_horiz_var(um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
!  Local variables
!
INTEGER :: ndims           ! Number of dimensions in NetCDF file
INTEGER :: dimlen1,dimlen2 ! Size of dimensions
INTEGER :: dimid1          ! NetCDF Dimension id
INTEGER :: varid           ! Variable id of dimension
INTEGER :: nlen            ! Length of dimension before any processing
INTEGER :: idim1           ! First dimension value before any processing
INTEGER :: idim2           ! Last dimension value before any processing
INTEGER :: ix              ! STASH item index
INTEGER :: i,j             ! Loop indices
INTEGER :: start(1)        ! Start position of dimension output data
INTEGER :: count1(1)       ! Number of items of dimension output data
INTEGER :: start_bnds(2)   ! Start position of dimension bounds output data
INTEGER :: count_bnds(2)   ! Number of items of dimension bounds output data
INTEGER :: lbnd            ! Column/Row array lower bound
INTEGER :: index1,index2   ! Indices of grid type name in dimension name
INTEGER :: mean_opt        ! Meaning option
INTEGER :: ErrorStatus     ! Error code

REAL :: dim0                    ! First dimension value
REAL :: diminc                  ! Dimension increment value
REAL :: pole_lon,pole_lat       ! North pole lat/lon values
REAL, ALLOCATABLE :: dimdata(:) ! Variable to contain dimension values
REAL, ALLOCATABLE :: dimdata_bnds(:,:)
                                ! Variable to contain dimension bounds values

LOGICAL :: lexists  ! TRUE if selected attribute exists
LOGICAL :: lrotgrid ! TRUE if rotated grid
LOGICAL :: lcyclic  ! TRUE if cyclic EW BCs

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NCFILE_WRITE_HORIZ_VAR'
CHARACTER(LEN=nf90_max_name) :: dimname1,dimname2 ! lat/lon dimension names
CHARACTER(LEN=nf90_max_name) :: dimname3,dimname4 ! true lat/lon dimension names
CHARACTER(LEN=1)             :: axis              ! Dimension axis type
CHARACTER(LEN=2)             :: grid_type         ! STASH gridtype

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0
lcyclic = .NOT. elf

! Check if output is on a rotated grid and get north pole lat/lon
CALL nc_get_varid(um_file,'rotated_latitude_longitude',varid)
IF (varid == -1) THEN
  lrotgrid = .FALSE.
ELSE
  lrotgrid = .TRUE.
  CALL nc_get_att(um_file,varid,'grid_north_pole_longitude',pole_lon,lexists)
  CALL nc_get_att(um_file,varid,'grid_north_pole_latitude',pole_lat,lexists)
END IF

! Get number of dimensions in NetCDF file
CALL nc_inq_ndims(um_file,ndims)
IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(A,A,A,I4)') &
            'Number of dimensions in NetCDF file ',TRIM(um_file % filename), &
            ' is ',ndims
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Loop over each dimension, i is netCDF dimid
DO i=1,ndims

  ! Get dimension name and length for this dimid
  CALL nc_inq_dim(um_file,i,dimname1,dimlen1)

  ! Not interested in bounds dimension
  IF (dimname1(1:6) == 'bounds') CYCLE

  ! Get id of variable with same name as dimension
  CALL nc_get_varid(um_file,dimname1,varid)
 
  ! Get axis attribute
  CALL nc_get_att(um_file,varid,'axis',axis,lexists)
  IF (.NOT. lexists) CYCLE

  ! This routine is for processing lat/long dimensions, so only select those
  IF (axis == 'X' .OR. axis == 'Y') THEN

    ! Get stash index value
    ix = um_file % nc_meta % stash_index_dim(i)

    IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
            'Horizontal dimension name ',TRIM(dimname1), &
            ' of size ',dimlen1,' on unit ',um_file % UNIT,' found'
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,A,A,I5)') &
                'STASH index for horizontal dimension ',TRIM(dimname1),' is ',ix
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Intialise start and count variables needed to write netCDF data
    start(1) = 1
    count1(1) = dimlen1
    start_bnds(1) = 1
    start_bnds(2) = 1
    count_bnds(1) = 2
    count_bnds(2) = dimlen1

    ! Extract grid type from dimension name
    index1 = INDEX(dimname1,'itude_')
    IF (index1 > 1) THEN
      index1 = index1 + 6
      index2 = INDEX(dimname1(index1:),'_')
      IF (index2 == 0) THEN
        index2 = LEN_TRIM(dimname1)
      ELSE
        index2 = index1 + index2 - 2
      END IF
      IF (index2-index1 == 0 .OR. index2-index1 == 1) THEN
        grid_type = dimname1(index1:index2)
      ELSE
        ErrorStatus = 99
      END IF
    ELSE
      ErrorStatus = 99
    END IF
    IF (ErrorStatus == 99) THEN
      CALL ereport(RoutineName,ErrorStatus, &
                  'Dimension name '//TRIM(dimname1)//' malformed')
    ELSE IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A)') 'Grid type is '//TRIM(grid_type)
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Get the first and last point for this dimension
    ! Values are input values before any meaning options
    IF (axis == 'X') THEN
      idim1 = stlist(st_west_code,ix)
      idim2 = stlist(st_east_code,ix)
    ELSE IF (axis == 'Y') THEN
      idim1 = stlist(st_south_code,ix)
      idim2 = stlist(st_north_code,ix)
    END IF
    IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A,A1,A,I4,I4)') &
            'Input domain limits for ',axis,' axis are ',idim1,idim2
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    mean_opt = (stlist(st_gridpoint_code,ix)/block_size)*block_size

    ! Allocate arrays to hold dimension and bounds data values
    ALLOCATE(dimdata(dimlen1))
    ALLOCATE(dimdata_bnds(2,dimlen1))

    IF (l_regular) THEN ! Regular grid, global or LAM

      ! Calculate first grid point, number of grid points and grid increment
      ! for a given regular grid type and axis
      IF (axis == 'X') THEN
        diminc=h_a_ewspace
        nlen=global_row_length
        SELECT CASE (grid_type)
        CASE (nc_horiz_id_theta,nc_horiz_id_v)
          dim0=h_a_firstlong+diminc/2.0
        CASE (nc_horiz_id_u,nc_horiz_id_uv)
          dim0=h_a_firstlong
        CASE (nc_horiz_id_river)
          diminc=river_row_length/dimlen1
          nlen=river_row_length
          dim0=h_a_firstlong+diminc/2.0
        CASE DEFAULT
          CYCLE
        END SELECT
      ELSE IF (axis == 'Y') THEN
        diminc=h_a_nsspace
        SELECT CASE (grid_type)
        CASE (nc_horiz_id_theta,nc_horiz_id_u)
          nlen=global_rows
          dim0=h_a_firstlat+diminc/2.0
        CASE (nc_horiz_id_v,nc_horiz_id_uv)
          nlen=global_rows+1
          dim0=h_a_firstlat
        CASE (nc_horiz_id_river)
          diminc=river_rows/dimlen1
          nlen=river_rows
          dim0=h_a_firstlat+diminc/2.0
        CASE DEFAULT
          CYCLE
        END SELECT
      END IF

      ! Calculate values for this dimension
      IF (axis == 'X' .AND. idim1 > idim2 .AND. lcyclic) THEN
        ! Field wraps east-west so make dimension values continuous
        ! East portion of dimension
        DO j=1,nlen-idim1+1
          dimdata(j) = dim0 + (j+idim1-2)*diminc
        END DO
        ! West portion of dimension shifted 360 degrees
        DO j=nlen-idim1+2,dimlen1
          dimdata(j) = dim0 + (j+idim1-nlen-2)*diminc + 360.0
        END DO
      ELSE
        ! Calculate dimension values from whole grid
        DO j=1,dimlen1
          dimdata(j) = dim0 + (j+idim1-2)*diminc
        END DO
      END IF

      ! Calculate bounds for the dimension in the same way as above
      ! dimdata_bnds(1,:) - is lower bound, half a grid point below grid value
      ! dimdata_bnds(2,:) - is upper bound, half a grid point above grid value
      IF (dimlen1 == 1 .AND.                                          &
         ((mean_opt == zonal_mean_base .AND. axis == 'X') .OR.        &
          (mean_opt == merid_mean_base .AND. axis == 'Y') .OR.        &
           mean_opt == field_mean_base .OR.                           &
           mean_opt == global_mean_base)) THEN
        ! For means, bounds are first and last points of input grid
        dimdata_bnds(1,1) = dimdata(1)
        dimdata_bnds(2,1) = dim0 + (idim2-1)*diminc
        IF (axis == 'X' .AND. idim1 > idim2 .AND. lcyclic) THEN
          ! Shift upper bound 360 degrees for east-west wrapped dimension
          dimdata_bnds(2,1) = dimdata_bnds(2,1) + 360.0
        END IF
      ELSE
        DO j=1,dimlen1
          dimdata_bnds(1,j) = dimdata(j) - 0.5*diminc
          dimdata_bnds(2,j) = dimdata(j) + 0.5*diminc
        END DO
      END IF

    ELSE ! Variable resolution grid

      ! Get values for this dimension
      IF (axis == 'X') THEN
        SELECT CASE (grid_type)
        CASE (nc_horiz_id_theta,nc_horiz_id_v)
          lbnd = jlambda_input_p+idim1-1
        CASE (nc_horiz_id_u,nc_horiz_id_uv)
          lbnd = jlambda_input_u+idim1-1
        CASE DEFAULT
          CYCLE
        END SELECT
        dimdata(1:dimlen1) = a_coldepc(lbnd:lbnd+dimlen1-1)
      ELSE IF (axis == 'Y') THEN
        SELECT CASE (grid_type)
        CASE (nc_horiz_id_theta,nc_horiz_id_u)
          lbnd = jphi_input_p+idim1-1
        CASE (nc_horiz_id_v,nc_horiz_id_uv)
          lbnd = jphi_input_v+idim1-1
        CASE DEFAULT
          CYCLE
        END SELECT
        dimdata(1:dimlen1) = a_rowdepc(lbnd:lbnd+dimlen1-1)
      END IF

      ! Get bounds for the dimension
      ! dimdata_bnds(1,:) - is lower bound
      ! dimdata_bnds(2,:) - is upper bound
      IF (dimlen1 == 1 .AND.                                          &
         ((mean_opt == zonal_mean_base .AND. axis == 'X') .OR.        &
          (mean_opt == merid_mean_base .AND. axis == 'Y') .OR.        &
           mean_opt == field_mean_base .OR.                           &
           mean_opt == global_mean_base)) THEN
        ! For means, bounds are first and last points of input grid
        dimdata_bnds(1,1) = dimdata(1)
        IF (axis == 'X') THEN
          SELECT CASE (grid_type)
          CASE (nc_horiz_id_theta,nc_horiz_id_v)
            dimdata_bnds(2,1) = a_coldepc(jlambda_input_p+idim2-1)
          CASE (nc_horiz_id_u,nc_horiz_id_uv)
            dimdata_bnds(2,1) = a_coldepc(jlambda_input_u+idim2-1)
          END SELECT
        ELSE IF (axis == 'Y') THEN
          SELECT CASE (grid_type)
          CASE (nc_horiz_id_theta,nc_horiz_id_u)
            dimdata_bnds(2,1) = a_rowdepc(jphi_input_p+idim2-1)
          CASE (nc_horiz_id_v,nc_horiz_id_uv)
            dimdata_bnds(2,1) = a_rowdepc(jphi_input_v+idim2-1)
          END SELECT
        END IF
      ELSE
        IF (axis == 'X') THEN
          SELECT CASE (grid_type)
          CASE (nc_horiz_id_theta,nc_horiz_id_v)
            ! U grid longitude values are bounds to Theta grid longitude values
            lbnd = jlambda_input_u+idim1-1
            dimdata_bnds(1,1:dimlen1) = a_coldepc(lbnd:lbnd+dimlen1-1)
            dimdata_bnds(2,1:dimlen1-1) = a_coldepc(lbnd+1:lbnd+dimlen1-1)
            IF (idim2 == global_row_length) THEN
              ! Last dimension value has no upper bound array value
              ! so extrapolate
              dimdata_bnds(2,dimlen1) = &
                2*a_coldepc(jlambda_input_p+global_row_length-1) - &
                  a_coldepc(jlambda_input_u+global_row_length-1)
            ELSE
              dimdata_bnds(2,dimlen1) = a_coldepc(lbnd+dimlen1)
            END IF
          CASE (nc_horiz_id_u,nc_horiz_id_uv)
            ! Theta grid longitude values are bounds to U grid longitude values
            lbnd = jlambda_input_p+idim1-2
            IF (idim1 == 1) THEN
              ! First dimension value has no lower bound array value
              ! so extrapolate
              dimdata_bnds(1,1) = &
                2*a_coldepc(jlambda_input_u) - a_coldepc(jlambda_input_p)
            ELSE
              dimdata_bnds(1,1) = a_coldepc(lbnd)
            END IF
            dimdata_bnds(1,2:dimlen1) = a_coldepc(lbnd+1:lbnd+dimlen1-1)
            dimdata_bnds(2,1:dimlen1) = a_coldepc(lbnd+1:lbnd+dimlen1)
          END SELECT
        ELSE IF (axis == 'Y') THEN
          SELECT CASE (grid_type)
          CASE (nc_horiz_id_theta,nc_horiz_id_u)
            ! V grid latitude values are bounds to Theta grid latitude values
            lbnd = jphi_input_v+idim1-1
            dimdata_bnds(1,1:dimlen1) = a_rowdepc(lbnd:lbnd+dimlen1-1)
            dimdata_bnds(2,1:dimlen1) = a_rowdepc(lbnd+1:lbnd+dimlen1)
          CASE (nc_horiz_id_v,nc_horiz_id_uv)
            ! Theta grid latitude values are bounds to V grid latitude values
            lbnd = jphi_input_p+idim1-2
            IF (idim1 == 1) THEN
              ! First dimension value has no lower bound array value
              ! so extrapolate
              dimdata_bnds(1,1) = &
                2*a_rowdepc(jphi_input_v) - a_rowdepc(jphi_input_p)
            ELSE
              dimdata_bnds(1,1) = a_rowdepc(lbnd)
            END IF
            dimdata_bnds(1,2:dimlen1) = a_rowdepc(lbnd+1:lbnd+dimlen1-1)
            dimdata_bnds(2,1:dimlen1-1) = a_rowdepc(lbnd+1:lbnd+dimlen1-1)
            IF (idim2 == global_rows+1) THEN
              ! Last dimension value has no upper bound array value
              ! so extrapolate
              dimdata_bnds(2,dimlen1) = &
                2*a_rowdepc(jphi_input_v+global_rows) - &
                  a_rowdepc(jphi_input_p+global_rows-1)
            ELSE
              dimdata_bnds(2,dimlen1) = a_rowdepc(lbnd+dimlen1)
            END IF
          END SELECT
        END IF
      END IF
    END IF

    ! Write dimension information to NetCDF file
    CALL nc_put_var(um_file,varid,dimdata,start,count1)

    ! Get dimension bounds variable id
    CALL nc_get_varid(um_file,TRIM(dimname1)//'_bounds',varid)

    ! Write Dimension bounds information to NetCDF file
    CALL nc_put_var(um_file,varid,dimdata_bnds,start_bnds,count_bnds)

    DEALLOCATE(dimdata)
    DEALLOCATE(dimdata_bnds)

  END IF
END DO

IF (lrotgrid) THEN

  ! Write out 2D true lat/lon fields to netCDF file
  DO i=1,ndims
  
    ! Get dimension name and length for this dimid
    CALL nc_inq_dim(um_file,i,dimname2,dimlen2)
  
    ! Not interested in bounds dimension
    IF (dimname2(1:6) == 'bounds') CYCLE
  
    ! Get id of variable with same name as dimension
    CALL nc_get_varid(um_file,dimname2,varid)
   
    ! Get axis attribute
    CALL nc_get_att(um_file,varid,'axis',axis,lexists)
    IF (.NOT. lexists) CYCLE

    ! Use latitude dimension to find matching longitude
    ! Only write true lat/lon fields when number of latitudes is > 1
    IF (axis == 'Y' .AND. dimlen2 > 1) THEN

      IF (printstatus >= prdiag) THEN
        WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
              'Latitude dimension name ',TRIM(dimname2), &
              ' of size ',dimlen2,' on unit ',um_file % UNIT,' found'
        CALL umPrint(umMessage,src=RoutineName)
      END IF

      ! Find longitude variable which goes with this latitude variable
      dimid1 = um_file % nc_meta % dimid_lon(i)
      IF (dimid1 == -1) CYCLE ! No matching longitude variable

      ! Get dimension name and length for this longitude dimid
      CALL nc_inq_dim(um_file,dimid1,dimname1,dimlen1)
      IF (printstatus >= prdiag) THEN
        WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
              'Matching longitude dimension name ',TRIM(dimname1), &
              ' of size ',dimlen1,' on unit ',um_file % UNIT,' found'
        CALL umPrint(umMessage,src=RoutineName)
      END IF
      ! Only write true lat/lon fields when number of longitudes is > 1
      IF (dimlen1 <= 1) CYCLE

      ! True longitude dimension name is rotated longitude dimension name
      ! with 'grid_' removed from start of string
      IF (dimname1(1:5) /= 'grid_') THEN
        ErrorStatus = -100
        CALL ereport(RoutineName,ErrorStatus, &
                     TRIM(dimname1)//" does not start with 'grid_'")
        CYCLE
      END IF
      dimname3 = dimname1(6:)

      ! True latitude dimension name is rotated latitude dimension name
      ! with 'grid_' removed from start of string
      IF (dimname2(1:5) /= 'grid_') THEN
        ErrorStatus = -101
        CALL ereport(RoutineName,ErrorStatus, &
                     TRIM(dimname2)//" does not start with 'grid_'")
        CYCLE
      END IF
      dimname4 = dimname2(6:)

      CALL ncfile_write_truell(um_file, &
                               dimname1,dimname2, &
                               dimname3,dimname4, &
                               dimlen1,dimlen2, &
                               pole_lon,pole_lat)

    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_horiz_var

! Take rotated latitude/longitude grid and 
! calculate true latitude/longitude grid.
! Write true latitude/longitude grid to netCDF file

SUBROUTINE ncfile_write_truell(um_file, &
                               lonnamein,latnamein, &
                               lonnameout,latnameout, &
                               nlon,nlat, &
                               lon_pole,lat_pole)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
CHARACTER(LEN=*), INTENT(IN)  :: lonnamein  ! Input longitude NetCDF name 
CHARACTER(LEN=*), INTENT(IN)  :: latnamein  ! Input latitude NetCDF name 
CHARACTER(LEN=*), INTENT(IN)  :: lonnameout ! Output longitude NetCDF name 
CHARACTER(LEN=*), INTENT(IN)  :: latnameout ! Output latitude NetCDF name 
INTEGER         , INTENT(IN)  :: nlon       ! Number of longitude points
INTEGER         , INTENT(IN)  :: nlat       ! Number of latitude points
REAL            , INTENT(IN)  :: lon_pole   ! North pole longitude
REAL            , INTENT(IN)  :: lat_pole   ! North pole latitude

!
!  Local variables
!
INTEGER :: start(2)  ! Start position of NetCDF output data for each dimension
INTEGER :: count1(2) ! Number of items of NetCDF output data for each dimension
INTEGER :: varid     ! Latitude/Longitude variable ID
INTEGER :: i,j       ! Loop indices

REAL, ALLOCATABLE :: alonin(:,:)  ! Input longitude values on rotated grid
REAL, ALLOCATABLE :: alatin(:,:)  ! Input latitude values on rotated grid
REAL, ALLOCATABLE :: alonout(:,:) ! Output longitude values on true lat/lon grid
REAL, ALLOCATABLE :: alatout(:,:) ! Output latitude values on true lat/lon grid

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NCFILE_WRITE_TRUELL'

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(alonin(nlon,nlat))
ALLOCATE(alatin(nlon,nlat))
ALLOCATE(alonout(nlon,nlat))
ALLOCATE(alatout(nlon,nlat))

! Intialise start and count variables needed to write netCDF data
start(1) = 1
start(2) = 1
count1(1) = nlon
count1(2) = nlat

! Get Input longitude values
CALL nc_get_varid(um_file,lonnamein,varid)
CALL nc_get_var(um_file,varid,alonin(:,1),start(1:1),count1(1:1))

! Expand to 2D array for rotate_eq_to_latlon
DO j=2,nlat
  DO i=1,nlon
    alonin(i,j) = alonin(i,1)
  END DO
END DO

! Get Input latitude values
CALL nc_get_varid(um_file,latnamein,varid)
CALL nc_get_var(um_file,varid,alatin(1,:),start(2:2),count1(2:2))

! Expand to 2D array for rotate_eq_to_latlon
DO j=1,nlat
  DO i=2,nlon
    alatin(i,j) = alatin(1,j)
  END DO
END DO

! Calculate true lat/lon values
CALL rotate_eq_to_latlon  &
            (alatin,alonin,alatout,alonout,lat_pole,lon_pole,nlon*nlat)

! Convert longitude range to -180,180
alonout = MODULO(alonout+180.0,360.0)-180.0

! Write true lat/lon values to NetCDF file
CALL nc_get_varid(um_file,lonnameout,varid)
CALL nc_put_var(um_file,varid,alonout,start,count1)
CALL nc_get_varid(um_file,latnameout,varid)
CALL nc_put_var(um_file,varid,alatout,start,count1)

DEALLOCATE(alonin)
DEALLOCATE(alatin)
DEALLOCATE(alonout)
DEALLOCATE(alatout)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_truell

END MODULE ncfile_write_horiz_var_mod
