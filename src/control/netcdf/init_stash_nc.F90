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
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: NetCDF output
!
!  Purpose: Add netCDF attributes defined by STASH to
!           each netCDF file

MODULE init_stash_nc_mod

! External netCDF library module
USE netcdf,                     ONLY: nf90_max_name
USE parkind1,                   ONLY: jprb,jpim
USE yomhook,                    ONLY: lhook,dr_hook
USE file_manager,               ONLY: &
    um_file_type,get_file_by_unit,init_file_loop
USE umprintmgr,                 ONLY: umprint,ummessage,printstatus,prdiag
USE missing_data_mod,           ONLY: rmdi,imdi
USE stash_model_mod,            ONLY: h_a_polelat,h_a_polelong
USE stextend_mod,               ONLY: in_s
USE stparam_mod,                ONLY: &
    st_output_code,st_diag_address,st_output_type,st_netcdf
USE cstash_mod,                 ONLY: itim_b,idom_b,timpro,dompro
USE ppxlook_mod,                ONLY: exppxi,exppxc
USE cppxref_mod,                ONLY: ppxref_charlen
USE submodel_mod,               ONLY: atmos_im
USE cppxref_mod,                ONLY: &
    ppx_data_type,ppx_grid_type,ppx_lv_code,ppx_pt_code,ppx_timavail_code, &
    ppx_type_real,ppx_type_int,ppx_type_log,ppx_packing_acc
USE stash_array_mod,            ONLY: stlist,stindex,nitems,nsects
USE nlstcall_nc_namelist_mod,   ONLY: ncformat,l_checksum
USE um_version_mod,             ONLY: um_version_char
USE umnetcdf_mod,               ONLY: &
    ncformat_netcdf4_classic_model,ncformat_netcdf4, &
    nc_time_index,nc_dim_used,nc_inq_dimlen, &
    nc_create_var,nc_var_type_real,nc_type_int, &
    nc_def_var_chunking,nc_def_var_deflate,nc_def_var_fletcher32, &
    nc_put_att,set_varid_val,nc_file_enddef,nc_file_sync, &
    nc_max_attlen,nc_max_var_dims
USE cf_metadata_mod,            ONLY: standard_name,units,cf_extra_info
USE ncfile_write_horiz_dim_mod, ONLY: ncfile_write_horiz_dim
USE ncfile_write_vert_dim_mod,  ONLY: ncfile_write_vert_dim
USE ncfile_write_pseudo_dim_mod,ONLY: ncfile_write_pseudo_dim
USE ncfile_write_time_dim_mod,  ONLY: ncfile_write_time_dim
USE ncfile_write_horiz_var_mod, ONLY: ncfile_write_horiz_var
USE ncfile_write_vert_var_mod,  ONLY: ncfile_write_vert_var
USE ncfile_write_pseudo_var_mod,ONLY: ncfile_write_pseudo_var
USE ncfile_write_time_var_mod,  ONLY: ncfile_write_time_var

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_STASH_NC_MOD'

CONTAINS

SUBROUTINE init_stash_nc ()

IMPLICIT NONE

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle
!
! Local variables
!
! Length of STASH variable in MSI format (mXXsYYiZZZ)
INTEGER, PARAMETER :: stash_msi_len = 10
INTEGER :: dimids(nc_max_var_dims)     ! NetCDF dimension IDs array
INTEGER :: chunksizes(nc_max_var_dims) ! Chunk size of each dimension
INTEGER :: im_index            ! Internal model index
INTEGER :: data_type           ! STASH datatype
INTEGER :: grid_type           ! STASH gridtype
INTEGER :: level_code          ! STASH level type code
INTEGER :: pseudo_level_code   ! Pseudo level type code
INTEGER :: time_avail          ! Time availability code
INTEGER :: ppunit              ! Output unit
INTEGER :: pptype              ! Output file type
INTEGER :: packing_code        ! Indicates type of packing to use
INTEGER :: nccomp_level        ! Deflate level, value between 0 and 9
INTEGER :: comp_accrcy         ! Packing accuracy in power of 2
INTEGER :: idim1               ! Number of NetCDF dimensions
INTEGER :: dimid1              ! First netCDF dimension ID
INTEGER :: dimid2              ! Second netCDF dimension ID
INTEGER :: varid               ! NetCDF variable id
INTEGER :: nc_data_type        ! Netcdf datatype
INTEGER :: idiag               ! STASH diagnostic address
INTEGER :: ix                  ! STASH item index
INTEGER :: i                   ! loop index over STASH item codes
INTEGER :: j                   ! loop index over STASH section codes

LOGICAL :: l_compress   ! .TRUE. if netCDF file uses compression
LOGICAL :: l_shuffle    ! Switch to turn on netCDF4 shuffle filter
LOGICAL :: l_chunked    ! .TRUE. if data uses netCDF4 chunking

CHARACTER(LEN=*), PARAMETER   :: RoutineName = "INIT_STASH_NC"
CHARACTER(LEN=ppxref_charlen) :: stash_name   ! STASH long_name attribute
CHARACTER(LEN=stash_msi_len)  :: stash_msi    ! STASH variable in MSI format
CHARACTER(LEN=nf90_max_name)  :: varname      ! STASH variable name
CHARACTER(LEN=nc_max_attlen)  :: coordinates  ! CF coordinates attribute
CHARACTER(LEN=nc_max_attlen)  :: cell_methods ! CF cell_methods attribute

TYPE(um_file_type), POINTER :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

im_index = 1

! Loop over all STASH entries
DO j=0,nsects
  DO i=1,nitems
    IF (in_s(1,atmos_im,j,i) >= 1) THEN

      stash_name        = exppxc(atmos_im,j,i)
      data_type         = exppxi(atmos_im,j,i,ppx_data_type)
      grid_type         = exppxi(atmos_im,j,i,ppx_grid_type)
      level_code        = exppxi(atmos_im,j,i,ppx_lv_code)
      pseudo_level_code = exppxi(atmos_im,j,i,ppx_pt_code)
      time_avail        = exppxi(atmos_im,j,i,ppx_timavail_code)

      ! Reset all values of stash_count
      NULLIFY(um_file)
      um_file => init_file_loop(handler="netcdf")
      DO WHILE (ASSOCIATED(um_file))
        um_file % nc_meta % stash_count = 0
        ! Increment the pointer to the next file
        um_file => um_file % next
      END DO

      DO ix=stindex(1,i,j,im_index), &
            stindex(1,i,j,im_index)+stindex(2,i,j,im_index)-1

        ppunit = -stlist(st_output_code,ix)
        pptype =  stlist(st_output_type,ix)
        IF (ppunit >  0 .AND. pptype == st_netcdf) THEN ! NetCDF output unit

          NULLIFY(um_file) 
          um_file => get_file_by_unit(ppunit,handler="netcdf")

          IF (um_file % meta % initialise .AND. &
              um_file % meta % init_steps /= 0 .AND. &
              .NOT. um_file % meta % partially_written) THEN

            um_file % nc_meta % stash_count = &
              um_file % nc_meta % stash_count + 1
            packing_code = um_file % meta % packing_code
            l_compress = um_file % nc_meta % l_compress

            SELECT CASE(data_type)
            CASE (ppx_type_real)
              nc_data_type = nc_var_type_real
            CASE (ppx_type_int,ppx_type_log)
              nc_data_type = nc_type_int
            END SELECT

            IF (packing_code >= 100) THEN
              packing_code = packing_code - 100
            END IF
            IF (packing_code > 0) THEN
              comp_accrcy = exppxi(atmos_im,j,i, &
                                   ppx_packing_acc+packing_code-1)
            ELSE
              comp_accrcy = -99
            END IF
            idiag = stlist(st_diag_address,ix)
            IF (printstatus >= prdiag) THEN
              WRITE(umMessage,'(A,A,I4,A,I4,I4,A,I6)') &
                        TRIM(stash_name),': on unit ',ppunit, &
                        ' STASH section,item = ',j,i, &
                        ' STASH index = ',ix
              CALL umPrint(umMessage,src=RoutineName)
              WRITE(umMessage,'(2X,3(A,I4))') &
                        'STASH count = ',um_file % nc_meta % stash_count, &
                        ' Accuracy = ',comp_accrcy, &
                        ' Diagnostic = ',idiag
              CALL umPrint(umMessage,src=RoutineName)
              WRITE(umMessage,'(2X,4(A))') &
                        'Domain profile = ',TRIM(dompro(idom_b(idiag))), &
                        ' Time profile = ',TRIM(timpro(itim_b(idiag)))
              CALL umPrint(umMessage,src=RoutineName)
            END IF

            coordinates = ' '
            cell_methods = ' '
            idim1=0
            nc_time_index(ix) = 0
            nc_dim_used(:,ix) = .FALSE.

            ! Create netCDF variable name from the model ID, 
            ! STASH section, STASH item and STASH count
            IF (um_file % nc_meta % stash_count == 1) THEN
              WRITE(varname,'(SS,A,I2.2,A,I2.2,A,I3.3)') &
                    'STASH_m',atmos_im,'s',j,'i',i
            ELSE
              WRITE(varname,'(SS,A,I2.2,A,I2.2,A,I3.3,A,I0)') &
                    'STASH_m',atmos_im,'s',j,'i',i,"_", &
                    um_file % nc_meta % stash_count
            END IF

            ! Write horizontal grid dimension information to NetCDF file
            CALL ncfile_write_horiz_dim(um_file,grid_type,ix,         &
                                        coordinates,cell_methods,     &
                                        dimid1,dimid2)
            IF (dimid1 > 0) THEN
              nc_dim_used(1,ix) = .TRUE.
              idim1=idim1+1
              dimids(idim1) = dimid1
              CALL nc_inq_dimlen(um_file,dimid1,chunksizes(idim1))
            END IF
            IF (dimid2 > 0) THEN
              nc_dim_used(2,ix) = .TRUE.
              idim1=idim1+1
              dimids(idim1) = dimid2
              CALL nc_inq_dimlen(um_file,dimid2,chunksizes(idim1))
            END IF

            ! Write vertical dimension information to NetCDF file
            CALL ncfile_write_vert_dim(um_file,level_code,ix,         &
                                       cf_extra_info(j,i),            &
                                       coordinates,cell_methods,      &
                                       dimid1)
            IF (dimid1 > 0) THEN
              nc_dim_used(3,ix) = .TRUE.
              idim1=idim1+1
              dimids(idim1) = dimid1
              chunksizes(idim1) = 1
            END IF

            ! Write pseudo dimension information to NetCDF file
            CALL ncfile_write_pseudo_dim(um_file,pseudo_level_code,   &
                                         ix,coordinates,dimid1)
            IF (dimid1 > 0) THEN
              nc_dim_used(4,ix) = .TRUE.
              idim1=idim1+1
              dimids(idim1) = dimid1
              chunksizes(idim1) = 1
            END IF

            ! Write time dimension information to NetCDF file
            CALL ncfile_write_time_dim(um_file,time_avail,j,ix,       &
                                       cell_methods,dimid1)
            IF (dimid1 > 0) THEN
              nc_dim_used(5,ix) = .TRUE.
              idim1=idim1+1
              dimids(idim1) = dimid1
              chunksizes(idim1) = 1
            ELSE
              CYCLE
            END IF

            ! Create NetCDF variable for this STASH item
            CALL nc_create_var(um_file,varname,nc_data_type,           &
                               varid,dimids(1:idim1))

            ! Save varid value for later use
            CALL set_varid_val(ix,varid)

            ! If using netCDF4 output
            IF (ncformat == ncformat_netcdf4_classic_model .OR.       &
                ncformat == ncformat_netcdf4) THEN

              nccomp_level = um_file % nc_meta % nccomp_level
              l_shuffle = um_file % nc_meta % l_shuffle
              l_chunked = .TRUE.

              ! Set chunking information
              CALL nc_def_var_chunking(um_file,varid,                 &
                                       l_chunked,                     &
                                       chunksizes(1:idim1))

              IF (printstatus >= prdiag) THEN
                WRITE(umMessage,'(A,A,20(I8,1X))')                    &
                      TRIM(varname),' has chunk sizes ',chunksizes(1:idim1)
                CALL umPrint(umMessage,src=RoutineName)
              END IF

              ! Compress data if required
              IF (l_compress) THEN
                CALL nc_def_var_deflate(um_file,varid,l_shuffle,nccomp_level)
              END IF

              ! Compute checksum value
              IF (l_checksum) THEN
                CALL nc_def_var_fletcher32(um_file,varid,l_checksum)
              END IF
            END IF

            ! Write STASH name as long_name attribute for this STASH item
            CALL nc_put_att(um_file,varid,'long_name',stash_name)

            ! Write standard name attribute for this STASH item
            IF (standard_name(j,i)(1:9) /= 'undefined') THEN
              CALL nc_put_att(um_file,varid,'standard_name',          &
                              TRIM(standard_name(j,i)))
            END IF

            ! Write units attribute for this STASH item
            IF (units(j,i)(1:9) /= 'undefined') THEN
              CALL nc_put_att(um_file,varid,'units',TRIM(units(j,i)))
            END IF

            ! Write coordinates attribute for this STASH item
            IF (coordinates /= ' ') THEN
              CALL nc_put_att(um_file,varid,'coordinates',TRIM(coordinates))
            END IF

            ! Write cell_methods attribute for this STASH item
            IF (cell_methods /= ' ') THEN
                CALL nc_put_att(um_file,varid,'cell_methods',         &
                                TRIM(cell_methods))
            END IF

            ! Write grid_mapping attribute for this STASH item
            IF (h_a_polelat == 90.0 .AND. h_a_polelong == 0.0) THEN
              CALL nc_put_att(um_file,varid,'grid_mapping','grid_crs')
            ELSE
              CALL nc_put_att(um_file,varid,'grid_mapping',           &
                              'rotated_latitude_longitude')
            END IF

            ! Write _FillValue attribute for this STASH item
            SELECT CASE (data_type)
            CASE (ppx_type_real)
              CALL nc_put_att(um_file,varid,'_FillValue',rmdi)
            CASE (ppx_type_int)
              CALL nc_put_att(um_file,varid,'_FillValue',imdi)
            END SELECT

            ! Write um_version, um_stash_source attributes for this STASH item
            CALL nc_put_att(um_file,varid,'um_version',um_version_char)
            WRITE(stash_msi,'(A1,I2.2,A1,I2.2,A1,I3.3)')              &
                            'm',atmos_im,'s',j,'i',i
            CALL nc_put_att(um_file,varid,'um_stash_source',stash_msi)

            ! Write packing_method attribute for this STASH item
            IF (l_compress .AND. comp_accrcy > -99) THEN
              CALL nc_put_att(um_file,varid,'packing_method','quantization')
              CALL nc_put_att(um_file,varid,'precision_measure','binary')
              CALL nc_put_att(um_file,varid,'precision_value',comp_accrcy)
            ELSE
              CALL nc_put_att(um_file,varid,'packing_method','none')
            END IF

          END IF

        END IF

      END DO

    END IF
  END DO
END DO

NULLIFY(um_file)
um_file => init_file_loop(handler="netcdf")
DO WHILE (ASSOCIATED(um_file))

  IF (um_file % meta % initialise .AND. &
      um_file % meta % init_steps /= 0 .AND. &
      .NOT. um_file % meta % partially_written) THEN

    ! End define mode for netCDF file
    CALL nc_file_enddef(um_file)

    ! Write horizontal dimension data
    CALL ncfile_write_horiz_var(um_file)

    ! Write vertical dimension data
    CALL ncfile_write_vert_var(um_file)

    ! Write pseudo dimension data
    CALL ncfile_write_pseudo_var(um_file)

    ! Write time dimension data
    CALL ncfile_write_time_var(um_file)

    ! Sync netCDF file data to disk
    CALL nc_file_sync(um_file)

    ! Set the flag which indicates the file is active
    um_file % meta % partially_written = .TRUE.

  END IF
  
  ! Increment the pointer to the next file
  um_file => um_file % next
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_stash_nc

END MODULE init_stash_nc_mod
