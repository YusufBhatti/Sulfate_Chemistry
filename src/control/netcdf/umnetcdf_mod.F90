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

! A module to output files in NetCDF format
MODULE umnetcdf_mod

USE nlstcall_nc_namelist_mod, ONLY: ncformat,ncdim_prec,ncvar_prec
USE yomhook,                  ONLY: lhook,dr_hook
USE parkind1,                 ONLY: jprb,jpim
USE um_types,                 ONLY: integer32,real32,real64
USE filenamelength_mod,       ONLY: filenamelength
USE version_mod,              ONLY: nrecdp
USE missing_data_mod,         ONLY: rmdi,imdi
USE file_manager,             ONLY: get_file_by_unit,um_file_type
USE errormessagelength_mod,   ONLY: errormessagelength
USE ereport_mod,              ONLY: ereport
USE umPrintMgr,               ONLY: umPrint,umMessage,PrintStatus,PrNorm,PrDiag
USE io_configuration_mod,     ONLY: io_external_control
USE MPPIO_file_utils,         ONLY: file_touch,file_delete
    
IMPLICIT NONE

PRIVATE

! Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

INTEGER, PARAMETER :: ndiag_um = nrecdp   ! Maximum number of STASH list 
                                          ! records (prognostic + diagnostic)

! Maximum number of variable dimensions, lat, lon, lev, pseudo and time
INTEGER, PARAMETER, PUBLIC :: nc_max_var_dims = 5

! Parameters defining netCDF file open mode
INTEGER, PARAMETER, PUBLIC :: nc_mode_create = 0  ! nc_file_open mode is create
INTEGER, PARAMETER, PUBLIC :: nc_mode_open = 1    ! nc_file_open mode is open 
                                                  ! read/write
INTEGER, PARAMETER, PUBLIC :: nc_mode_open_ro = 2 ! nc_file_open mode is open 
                                                  ! read only

! Parameters defining netCDF format
INTEGER, PARAMETER, PUBLIC :: ncformat_netcdf3 = 1
INTEGER, PARAMETER, PUBLIC :: ncformat_netcdf3_64bit_offset = 2
INTEGER, PARAMETER, PUBLIC :: ncformat_netcdf4 = 3
INTEGER, PARAMETER, PUBLIC :: ncformat_netcdf4_classic_model = 4

! Parameters defining netCDF floating point precision
INTEGER, PARAMETER :: ncprec_32bit = 1
INTEGER, PARAMETER :: ncprec_64bit = 2

! Maximum length of a netCDF attribute character string
INTEGER, PARAMETER, PUBLIC :: nc_max_attlen = 1024

INTEGER, PUBLIC :: nc_dim_type_real        ! Real data type for dimensions
INTEGER, PUBLIC :: nc_var_type_real        ! Real data type for variables
INTEGER, PUBLIC :: nc_type_int             ! Integer data type
INTEGER, PUBLIC :: nc_time_index(ndiag_um) ! Index value for time dimension
INTEGER         :: varid_val(ndiag_um) = -1! variable id for a given STASH index

! TRUE if dimension number used for STASH diagnostic
LOGICAL, PUBLIC :: nc_dim_used(nc_max_var_dims,ndiag_um) 

CHARACTER(LEN=errormessagelength) :: nc_error_message ! Error message string
!$OMP THREADPRIVATE (nc_error_message)

! Interfaces to overload subroutines

! Write data to Netcdf file
INTERFACE nc_put_var
MODULE PROCEDURE      &
  nc_put_var_real_1d, &
  nc_put_var_real_2d, &
  nc_put_var_real_3d, &
  nc_put_var_real_4d, &
  nc_put_var_int_1d,  &
  nc_put_var_int_2d,  &
  nc_put_var_int_3d,  &
  nc_put_var_int_4d
END INTERFACE

! Read data from Netcdf file
INTERFACE nc_get_var
MODULE PROCEDURE      &
  nc_get_var_real_1d, &
  nc_get_var_real_2d, &
  nc_get_var_real_3d, &
  nc_get_var_real_4d, &
  nc_get_var_int_1d,  &
  nc_get_var_int_2d,  &
  nc_get_var_int_3d,  &
  nc_get_var_int_4d
END INTERFACE

! Write attributes to Netcdf file
INTERFACE nc_put_att
MODULE PROCEDURE      &
  nc_put_att_str,     &
  nc_put_att_int,     &
  nc_put_att_real
END INTERFACE

! Read attributes from Netcdf file
INTERFACE nc_get_att
MODULE PROCEDURE      &
  nc_get_att_str,     &
  nc_get_att_int,     &
  nc_get_att_real
END INTERFACE

! Only make the interfaces available publicly
PUBLIC nc_put_var
PUBLIC nc_get_var
PUBLIC nc_put_att
PUBLIC nc_get_att

! Other publicly available subroutines/functions
PUBLIC nc_file_open
PUBLIC nc_file_close
PUBLIC nc_file_sync
PUBLIC nc_file_redef
PUBLIC nc_file_enddef
PUBLIC nc_create_dim
PUBLIC nc_create_var
PUBLIC nc_inq_ndims
PUBLIC nc_inq_dim
PUBLIC nc_inq_dimname
PUBLIC nc_inq_dimlen
PUBLIC nc_get_dimid
PUBLIC nc_get_varid
PUBLIC nc_inq_var_dimid
PUBLIC nc_def_var_chunking
PUBLIC nc_inq_var_chunking
PUBLIC nc_def_var_deflate
PUBLIC nc_def_var_fletcher32
PUBLIC nc_is_file_open
PUBLIC set_varid_val
PUBLIC get_varid_val
PUBLIC check_cf_name
PUBLIC nc_set_data_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UMNETCDF_MOD'

CONTAINS

!----------------------------------------------------------------------
! Subroutine: nc_file_open
! Open NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_file_open(um_file,mode)

USE netcdf, ONLY: &
    nf90_clobber,nf90_write,nf90_nowrite, &
    nf90_64bit_offset,nf90_netcdf4,nf90_classic_model, &
    nf90_create,nf90_open,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)                     :: mode

INTEGER                       :: errcode
INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: mode_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_FILE_OPEN'
CHARACTER(LEN=filenamelength) :: filename 

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (nc_is_file_open(um_file)) THEN
  errcode = 100
  CALL nc_ereport(RoutineName, errcode ,                                      &
                  'Unit already open', um_file % UNIT)
END IF

ncid = -1
filename = um_file % filename

SELECT CASE (mode)
CASE (nc_mode_create)
  mode_nc = nf90_clobber
  SELECT CASE (ncformat)
  CASE (ncformat_netcdf3)
    WRITE(umMessage,'(A,A,A,I4)')'Creating netCDF3 file ',                    &
                                 TRIM(filename),' on unit ',um_file % UNIT
    CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  CASE (ncformat_netcdf3_64bit_offset)
    mode_nc = IOR(mode_nc,nf90_64bit_offset)
    WRITE(umMessage,'(A,A,A,I4)')'Creating netCDF3 64bit offset file ',       &
                                 TRIM(filename),' on unit ',um_file % UNIT
    CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  CASE (ncformat_netcdf4_classic_model)
    mode_nc = IOR(mode_nc,nf90_netcdf4)
    mode_nc = IOR(mode_nc,nf90_classic_model)
    WRITE(umMessage,'(A,A,A,I4)')'Creating netCDF4 classic model file ',      &
                                 TRIM(filename),' on unit ',um_file % UNIT
    CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  CASE (ncformat_netcdf4)
    mode_nc = IOR(mode_nc,nf90_netcdf4)
    WRITE(umMessage,'(A,A,A,I4)')'Creating netCDF4 file ',                    &
                                 TRIM(filename),' on unit ',um_file % UNIT
    CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  CASE DEFAULT
    WRITE(nc_error_message,'(A47,I3)')                                        &
      'Cannot open netCDF file, unknown format number ',ncformat
    errcode = 101
    CALL nc_ereport(RoutineName,errcode,nc_error_message,um_file % UNIT)
  END SELECT
  status_nc = nf90_create(filename, mode_nc, ncid)
CASE (nc_mode_open)
  mode_nc = nf90_write
  WRITE(umMessage,'(A,A,A,I4,A)')'Opening netCDF file ',                      &
                                 TRIM(filename),' on unit ',um_file % UNIT,   &
                                 ' with read/write access'
  CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  status_nc = nf90_open(filename, mode_nc, ncid)
CASE (nc_mode_open_ro)
  mode_nc = nf90_nowrite
  WRITE(umMessage,'(A,A,A,I4,A)')'Opening netCDF file ',                      &
                                 TRIM(filename),' on unit ',um_file % UNIT,   &
                                 ' with read only access'
  CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)
  status_nc = nf90_open(filename, mode_nc, ncid)
CASE DEFAULT
  WRITE(nc_error_message,'(A,I3)')                                            &
    'Cannot open netCDF file, unknown mode ',mode
  errcode = 102
  CALL nc_ereport(RoutineName,errcode,nc_error_message,um_file % UNIT)
END SELECT

IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

um_file % nc_meta % ncid = ncid
um_file % nc_meta % mode = mode

! Remove any .done files for this UNIT
IF (io_external_control .AND. mode /= nc_mode_open_ro) THEN
  CALL file_delete(TRIM(filename)//'.done',silent=.TRUE.)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_file_open

!----------------------------------------------------------------------
! Subroutine: nc_file_close
! Close NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_file_close(um_file)

USE netcdf, ONLY: nf90_close,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_FILE_CLOSE'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

WRITE(umMessage,'(A,I4)')'Closing NetCDF file '//TRIM(um_file % filename)// &
                         ' on unit = ',um_file % UNIT
CALL umPrint(umMessage,src='umnetcdf_mod',level=PrNorm)

status_nc = nf90_close(ncid)

IF (status_nc == nf90_noerr) THEN
  um_file % nc_meta % ncid = -1
  um_file % meta % partially_written = .FALSE.
ELSE
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

! Create a .done files for this UNIT
IF (io_external_control .AND. um_file % nc_meta % mode /= nc_mode_open_ro) THEN
  CALL file_touch(TRIM(um_file % filename)//'.done')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_file_close

!----------------------------------------------------------------------
! Subroutine: nc_file_sync
! Sync NetCDF file to disk
!----------------------------------------------------------------------

SUBROUTINE nc_file_sync(um_file)

USE netcdf, ONLY: nf90_sync,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_FILE_SYNC'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,I4)')'Syncing NetCDF file on unit = ',um_file % UNIT
CALL umPrint(umMessage,src='umnetcdf_mod',level=PrDiag)

ncid = um_file % nc_meta % ncid

status_nc = nf90_sync(ncid)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_file_sync

!----------------------------------------------------------------------
! Subroutine: nc_file_redef
! Start define mode for NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_file_redef(um_file)

USE netcdf, ONLY: nf90_redef,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_FILE_REDEF'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,I4)')'Start define mode for NetCDF file on unit = ', &
                         um_file % UNIT
CALL umPrint(umMessage,src='umnetcdf_mod',level=PrDiag)

ncid = um_file % nc_meta % ncid

status_nc = nf90_redef(ncid)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_file_redef

!----------------------------------------------------------------------
! Subroutine: nc_file_enddef
! End define mode for NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_file_enddef(um_file)

USE netcdf, ONLY: nf90_enddef,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_FILE_ENDDEF'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,I4)')'End define mode for NetCDF file on unit = ', &
                         um_file % UNIT
CALL umPrint(umMessage,src='umnetcdf_mod',level=PrDiag)

ncid = um_file % nc_meta % ncid

status_nc = nf90_enddef(ncid)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_file_enddef

!----------------------------------------------------------------------
! Subroutine: nc_create_dim
! Define NetCDF dimension
!----------------------------------------------------------------------

SUBROUTINE nc_create_dim(um_file,dimname,dimlen,dimid,lexists)

USE netcdf, ONLY: &
    nf90_def_dim,nf90_inq_dimid,nf90_unlimited,nf90_enameinuse,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
CHARACTER(LEN=*), INTENT(IN)  :: dimname  ! NetCDF dimension name
INTEGER, INTENT(IN)           :: dimlen
INTEGER, INTENT(OUT)          :: dimid
LOGICAL, INTENT(OUT), OPTIONAL:: lexists  ! Set to .TRUE. if 
                                          ! dimension already exists

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: dimid_nc, dimlen_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_CREATE_DIM'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (dimlen == -1) THEN
  dimlen_nc = nf90_unlimited
ELSE
  dimlen_nc = dimlen
END IF

ncid = um_file % nc_meta % ncid

IF (PrintStatus>=PrDiag) THEN
  WRITE(umMessage,'(A,A,A,I6)')'Creating dimension ',   &
                               TRIM(dimname),' dimlen = ',dimlen
  CALL umPrint(umMessage,src='umnetcdf_mod')
END IF
! Define dimension id
status_nc = nf90_def_dim(ncid, dimname, dimlen_nc, dimid_nc)
IF (status_nc == nf90_noerr) THEN
  IF (PRESENT(lexists)) lexists = .FALSE.
ELSE IF (status_nc == nf90_enameinuse) THEN
  IF (PRESENT(lexists)) lexists = .TRUE.
  IF (PrintStatus>=PrDiag) THEN
    WRITE(umMessage,'(A,A,A)')'Dimension ',TRIM(dimname),' already exists'
    CALL umPrint(umMessage,src='umnetcdf_mod')
  END IF
  ! Get existing dimension id
  status_nc = nf90_inq_dimid(ncid, dimname, dimid_nc)
END IF

IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

dimid = dimid_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_create_dim

!----------------------------------------------------------------------
! Subroutine: nc_create_var
! Define NetCDF variable
!----------------------------------------------------------------------

SUBROUTINE nc_create_var(um_file,varname,data_type,varid,dimids,lexists)

USE netcdf, ONLY: nf90_def_var,nf90_inq_varid,nf90_enameinuse,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
CHARACTER(LEN=*), INTENT(IN)  :: varname  ! NetCDF variable name
INTEGER, INTENT(IN)           :: data_type
INTEGER, INTENT(OUT)          :: varid
INTEGER, INTENT(IN), OPTIONAL :: dimids(:)
LOGICAL, INTENT(OUT), OPTIONAL:: lexists

INTEGER                       :: errcode
INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: xtype
INTEGER(integer32)            :: dimids_nc(nc_max_var_dims)
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_CREATE_VAR'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

xtype = data_type
ncid = um_file % nc_meta % ncid

IF (PrintStatus>=PrDiag) THEN
  WRITE(umMessage,'(A,A)')'Creating variable ',TRIM(varname)
  CALL umPrint(umMessage,src='umnetcdf_mod')
END IF
! Define variable id
IF (PRESENT(dimids)) THEN
  dimids_nc = dimids
  status_nc = nf90_def_var(ncid, varname, xtype,           &
                           dimids_nc(1:SIZE(dimids)), varid_nc)
ELSE
  status_nc = nf90_def_var(ncid, varname, xtype, varid_nc)
END IF
IF (status_nc == nf90_noerr) THEN
  IF (PRESENT(lexists)) lexists = .FALSE.
ELSE IF (status_nc == nf90_enameinuse) THEN
  IF (PRESENT(lexists)) lexists = .TRUE.
  IF (PrintStatus>=PrDiag) THEN
    WRITE(umMessage,'(A,A,A)')'Variable ',TRIM(varname),' already exists'
    CALL umPrint(umMessage,src='umnetcdf_mod')
  END IF
  ! Get existing variable id
  status_nc = nf90_inq_varid(ncid, varname, varid_nc)
END IF

IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

varid = varid_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_create_var

!----------------------------------------------------------------------
! Subroutine: nc_inq_ndims
! Return number of dimensions in NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_inq_ndims(um_file,ndims)

USE netcdf, ONLY: nf90_inquire,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(OUT)          :: ndims

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: ndims_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_NDIMS'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

status_nc = nf90_inquire(ncid, nDimensions=ndims_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

ndims = ndims_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_ndims

!----------------------------------------------------------------------
! Subroutine: nc_inq_dim
! Return dimension name and length for given dimension id
!----------------------------------------------------------------------

SUBROUTINE nc_inq_dim(um_file,dimid,dimname,dimlen)

USE netcdf, ONLY: nf90_inquire_dimension,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: dimid
CHARACTER(LEN=*), INTENT(OUT) :: dimname
INTEGER, INTENT(OUT)          :: dimlen

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: dimid_nc
INTEGER(integer32)            :: dimlen_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_DIM'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dimid_nc = dimid

ncid = um_file % nc_meta % ncid

status_nc = nf90_inquire_dimension(ncid, dimid_nc, &
                                   dimname, dimlen_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

dimlen = dimlen_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_dim

!----------------------------------------------------------------------
! Subroutine: nc_inq_dimname
! Return dimension name for given dimension id
!----------------------------------------------------------------------

SUBROUTINE nc_inq_dimname(um_file,dimid,dimname)

USE netcdf, ONLY:  nf90_inquire_dimension,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: dimid
CHARACTER(LEN=*), INTENT(OUT) :: dimname

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: dimid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_DIMNAME'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dimid_nc = dimid

ncid = um_file % nc_meta % ncid

status_nc = nf90_inquire_dimension(ncid, dimid_nc, NAME=dimname)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_dimname

!----------------------------------------------------------------------
! Subroutine: nc_inq_dimlen
! Return dimension length for given dimension id
!----------------------------------------------------------------------

SUBROUTINE nc_inq_dimlen(um_file,dimid,dimlen)

USE netcdf, ONLY: nf90_inquire_dimension,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: dimid
INTEGER, INTENT(OUT)          :: dimlen

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: dimid_nc
INTEGER(integer32)            :: dimlen_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_DIMLEN'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

dimid_nc = dimid

ncid = um_file % nc_meta % ncid

status_nc = nf90_inquire_dimension(ncid, dimid_nc, LEN=dimlen_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,UNIT=um_file % UNIT,nc_status=status_nc)
END IF

dimlen = dimlen_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_dimlen

!----------------------------------------------------------------------
! Subroutine: nc_get_dimid
! Return dimension id for given dimension name
!----------------------------------------------------------------------

SUBROUTINE nc_get_dimid(um_file,dimname,dimid)

USE netcdf, ONLY: nf90_inq_dimid,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
CHARACTER(LEN=*), INTENT(IN)  :: dimname  ! NetCDF dimension name
INTEGER, INTENT(OUT)          :: dimid    ! NetCDF dimension id

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: dimid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_DIMID'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

status_nc = nf90_inq_dimid(ncid, dimname, dimid_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_INQ_DIMID', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

dimid = dimid_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_dimid

!----------------------------------------------------------------------
! Subroutine: nc_get_varid
! Return variable id for given variable name
!----------------------------------------------------------------------

SUBROUTINE nc_get_varid(um_file,varname,varid)

USE netcdf, ONLY: nf90_inq_varid,nf90_enotvar,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
CHARACTER(LEN=*), INTENT(IN)  :: varname  ! NetCDF variable name
INTEGER, INTENT(OUT)          :: varid    ! NetCDF variable id

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VARID'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

status_nc = nf90_inq_varid(ncid, varname, varid_nc)
IF (status_nc == nf90_enotvar) THEN ! variable does not exist
  IF (PrintStatus>=PrDiag) THEN
    WRITE(umMessage,'(A,A,A)')'Varname ',TRIM(varname),' does not exist'
    CALL umPrint(umMessage,src='umnetcdf_mod')
  END IF
  varid = -1
ELSE IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_INQ_VARID', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
ELSE
  varid = varid_nc
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_varid

!----------------------------------------------------------------------
! Subroutine: nc_inq_var_dimid
! Return netCDF variable dimids for variable varid
!----------------------------------------------------------------------

SUBROUTINE nc_inq_var_dimid(um_file,varid,dimids,ndims)

USE netcdf, ONLY: nf90_inquire_variable,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
INTEGER, INTENT(OUT)          :: dimids(:)
INTEGER, INTENT(OUT)          :: ndims

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: dimids_nc(nc_max_var_dims)
INTEGER(integer32)            :: ndims_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_VAR_DIMID'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

status_nc = nf90_inquire_variable(ncid, varid_nc,     &
                                  ndims=ndims_nc,     &
                                  dimids=dimids_nc(1:SIZE(dimids)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_INQUIRE_VARIABLE', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

dimids(:) = dimids_nc(1:SIZE(dimids))
ndims = ndims_nc

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_var_dimid

!----------------------------------------------------------------------
! Subroutine: nc_put_var_real_1d
! Write real 1-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_real_1d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid      ! NetCDF variable id
REAL, INTENT(IN)              :: values(:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_REAL_1D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,          &
                         start_nc(1:SIZE(start_in)),      &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_real_1d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_real_2d
! Write real 2-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_real_2d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid       ! NetCDF variable id
REAL, INTENT(IN)              :: values(:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_REAL_2D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_real_2d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_real_3d
! Write real 3-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_real_3d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid          ! NetCDF variable id
REAL, INTENT(IN)              :: values (:,:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_REAL_3D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_real_3d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_real_4d
! Write real 2-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_real_4d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid           ! NetCDF variable id
REAL, INTENT(IN)              :: values(:,:,:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_REAL_4D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_real_4d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_int_1d
! Write integer 1-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_int_1d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid      ! NetCDF variable id
INTEGER, INTENT(IN)           :: values(:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_INT_1D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_int_1d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_int_2d
! Write integer 2-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_int_2d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid       ! NetCDF variable id
INTEGER, INTENT(IN)           :: values(:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_INT_2D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_int_2d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_int_3d
! Write integer 3-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_int_3d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid         ! NetCDF variable id
INTEGER, INTENT(IN)           :: values(:,:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_INT_3D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_int_3d

!----------------------------------------------------------------------
! Subroutine: nc_put_var_int_4d
! Write integer 4-dimensional variable to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_var_int_4d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_put_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid           ! NetCDF variable id
INTEGER, INTENT(IN)           :: values(:,:,:,:) ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_VAR_INT_4D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_put_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_var_int_4d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_real_1d
! Read real 1-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_real_1d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid      ! NetCDF variable id
REAL, INTENT(OUT)             :: values(:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_REAL_1D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_real_1d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_real_2d
! Read real 2-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_real_2d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid        ! NetCDF variable id
REAL, INTENT(OUT)             :: values(:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_REAL_2D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_real_2d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_real_3d
! Read real 3-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_real_3d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid          ! NetCDF variable id
REAL, INTENT(OUT)             :: values(:,:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_REAL_3D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_real_3d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_real_4d
! Read real 4-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_real_4d(um_file,varid,values, &
                              start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid            ! NetCDF variable id
REAL, INTENT(OUT)             :: values(:,:,:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_REAL_4D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,        &
                         start_nc(1:SIZE(start_in)),    &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_real_4d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_int_1d
! Read integer 1-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_int_1d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid      ! NetCDF variable id
INTEGER, INTENT(OUT)          :: values(:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_INT_1D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_int_1d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_int_2d
! Read integer 2-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_int_2d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid        ! NetCDF variable id
INTEGER, INTENT(OUT)          :: values(:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_INT_2D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_int_2d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_int_3d
! Read integer 3-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_int_3d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid          ! NetCDF variable id
INTEGER, INTENT(OUT)          :: values(:,:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_INT_3D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_int_3d

!----------------------------------------------------------------------
! Subroutine: nc_get_var_int_4d
! Read integer 4-dimensional variable from NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_var_int_4d(um_file,varid,values, &
                             start_in,count_in)

USE netcdf, ONLY: nf90_get_var,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid            ! NetCDF variable id
INTEGER, INTENT(OUT)          :: values(:,:,:,:)  ! Variable data values
INTEGER, INTENT(IN)           :: start_in(:)
INTEGER, INTENT(IN)           :: count_in(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: start_nc(nc_max_var_dims)
INTEGER(integer32)            :: count_nc(nc_max_var_dims)
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_VAR_INT_4D'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

start_nc(1:SIZE(start_in)) = start_in(:)
count_nc(1:SIZE(count_in)) = count_in(:)
status_nc = nf90_get_var(ncid, varid_nc, values,       &
                         start_nc(1:SIZE(start_in)),   &
                         count_nc(1:SIZE(count_in)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_VAR', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_var_int_4d

!----------------------------------------------------------------------
! Subroutine: nc_put_att_str
! Add string attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_att_str(um_file,varid,attname,values)

USE netcdf, ONLY: nf90_put_att,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
CHARACTER(LEN=*), INTENT(IN)  :: values   ! Attribute value

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_ATT_STR'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

status_nc = nf90_put_att(ncid, varid_nc, attname, values)

IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_att_str

!----------------------------------------------------------------------
! Subroutine: nc_put_att_int
! Add integer attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_att_int(um_file,varid,attname,values)

USE netcdf, ONLY: nf90_put_att,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
INTEGER, INTENT(IN)           :: values   ! Attribute value

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_ATT_INT'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

status_nc = nf90_put_att(ncid, varid_nc, attname, values)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_att_int

!----------------------------------------------------------------------
! Subroutine: nc_put_att_real
! Add real attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_put_att_real(um_file,varid,attname,values)

USE netcdf, ONLY: &
    nf90_inquire_variable,nf90_put_att,nf90_float,nf90_double,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
REAL, INTENT(IN)              :: values   ! Attribute value

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: nctype
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_PUT_ATT_REAL'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

status_nc = nf90_inquire_variable(ncid, varid_nc, xtype=nctype)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_INQUIRE_VARIABLE', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (nctype == nf90_float) THEN
  status_nc = nf90_put_att(ncid, varid_nc, attname, REAL(values,real32))
ELSE IF (nctype == nf90_double) THEN
  status_nc = nf90_put_att(ncid, varid_nc, attname, REAL(values,real64))
ELSE
  status_nc = nf90_put_att(ncid, varid_nc, attname, values)
END IF
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_PUT_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF
 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_put_att_real

!----------------------------------------------------------------------
! Subroutine: nc_get_att_str
! Add string attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_att_str(um_file,varid,attname,values,lexists)

USE netcdf, ONLY: nf90_get_att,nf90_enotatt,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
CHARACTER(LEN=*), INTENT(OUT) :: values   ! Attribute value
LOGICAL, INTENT(OUT)          :: lexists  ! .TRUE. if attribute does not exist

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_ATT_STR'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
lexists = .TRUE.

status_nc = nf90_get_att(ncid, varid_nc, attname, values)
IF (status_nc == nf90_enotatt) THEN ! attribute does not exist
  values = ''
  lexists = .FALSE.
ELSE IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_att_str

!----------------------------------------------------------------------
! Subroutine: nc_get_att_int
! Add integer attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_att_int(um_file,varid,attname,values,lexists)

USE netcdf, ONLY: nf90_get_att,nf90_enotatt,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
INTEGER, INTENT(OUT)          :: values   ! Attribute value
LOGICAL, INTENT(OUT)          :: lexists  ! .TRUE. if attribute does not exist

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_ATT_INT'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
lexists = .TRUE.

status_nc = nf90_get_att(ncid, varid_nc, attname, values)
IF (status_nc == nf90_enotatt) THEN ! attribute does not exist
  values = imdi
  lexists = .FALSE.
ELSE IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_att_int

!----------------------------------------------------------------------
! Subroutine: nc_get_att_real
! Add real attribute to NetCDF file
!----------------------------------------------------------------------

SUBROUTINE nc_get_att_real(um_file,varid,attname,values,lexists)

USE netcdf, ONLY: nf90_get_att,nf90_enotatt,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
CHARACTER(LEN=*), INTENT(IN)  :: attname  ! Attribute name
REAL, INTENT(OUT)             :: values   ! Attribute value
LOGICAL, INTENT(OUT)          :: lexists  ! .TRUE. if attribute does not exist

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_GET_ATT_REAL'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
lexists = .TRUE.

status_nc = nf90_get_att(ncid, varid_nc, attname, values)
IF (status_nc == nf90_enotatt) THEN ! attribute does not exist
  values = rmdi
  lexists = .FALSE.
ELSE IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_GET_ATT', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_get_att_real

!----------------------------------------------------------------------
! Subroutine: nc_def_var_chunking
! Set netCDF variable chunking parameters
!----------------------------------------------------------------------

SUBROUTINE nc_def_var_chunking(um_file,varid,lchunked,chunksizes)

USE netcdf, ONLY: nf90_def_var_chunking,nf90_chunked,nf90_contiguous,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
LOGICAL, INTENT(IN)           :: lchunked
INTEGER, INTENT(IN)           :: chunksizes(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: storage_nc
INTEGER(integer32)            :: chunksizes_nc(nc_max_var_dims)
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_DEF_VAR_CHUNKING'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
IF (lchunked) THEN
  storage_nc = nf90_chunked
ELSE
  storage_nc = nf90_contiguous
END IF
chunksizes_nc(1:SIZE(chunksizes)) = chunksizes(:)

status_nc = nf90_def_var_chunking(ncid, varid_nc, storage_nc,      &
                                  chunksizes_nc(1:SIZE(chunksizes)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_DEF_VAR_CHUNKING', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_def_var_chunking

!----------------------------------------------------------------------
! Subroutine: nc_inq_var_chunking
! Return netCDF variable chunking parameters
!----------------------------------------------------------------------

SUBROUTINE nc_inq_var_chunking(um_file,varid,lchunked,chunksizes)

USE netcdf, ONLY: nf90_inq_var_chunking,nf90_chunked,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
LOGICAL, INTENT(OUT)          :: lchunked
INTEGER, INTENT(OUT)          :: chunksizes(:)

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: storage_nc
INTEGER(integer32)            :: chunksizes_nc(nc_max_var_dims)
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_INQ_VAR_CHUNKING'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid

status_nc = nf90_inq_var_chunking(ncid, varid_nc, storage_nc,      &
                                  chunksizes_nc(1:SIZE(chunksizes)))
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_INQ_VAR_CHUNKING', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

lchunked = storage_nc == nf90_chunked
chunksizes(:) = chunksizes_nc(1:SIZE(chunksizes))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_inq_var_chunking

!----------------------------------------------------------------------
! Subroutine: nc_def_var_deflate
! Set netCDF variable compression options
!----------------------------------------------------------------------

SUBROUTINE nc_def_var_deflate(um_file,varid,lshuffle,deflate_level)

USE netcdf, ONLY: nf90_def_var_deflate,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
LOGICAL, INTENT(IN)           :: lshuffle
INTEGER, INTENT(IN)           :: deflate_level

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: shuffle_nc
INTEGER(integer32)            :: deflate_nc
INTEGER(integer32)            :: deflate_level_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_DEF_VAR_DEFLATE'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
IF (lshuffle) THEN
  shuffle_nc = 1
ELSE
  shuffle_nc = 0
END IF
deflate_nc = 1
deflate_level_nc = deflate_level

status_nc = nf90_def_var_deflate(ncid, varid_nc, shuffle_nc,      &
                                 deflate_nc, deflate_level_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_DEF_VAR_DEFLATE', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_def_var_deflate

!----------------------------------------------------------------------
! Subroutine: nc_def_var_fletcher32
! Set netCDF variable checksum option for compression data
!----------------------------------------------------------------------

SUBROUTINE nc_def_var_fletcher32(um_file,varid,lchecksum)

USE netcdf, ONLY: nf90_def_var_fletcher32,nf90_fletcher32,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN)           :: varid    ! NetCDF variable id
LOGICAL, INTENT(IN)           :: lchecksum

INTEGER(integer32)            :: ncid
INTEGER(integer32)            :: status_nc
INTEGER(integer32)            :: varid_nc
INTEGER(integer32)            :: checksum_nc
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_DEF_VAR_FLETCHER32'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

varid_nc = varid
IF (lchecksum) THEN
  checksum_nc = nf90_fletcher32
ELSE
  checksum_nc = 0
END IF

status_nc = nf90_def_var_fletcher32(ncid, varid_nc, checksum_nc)
IF (status_nc /= nf90_noerr) THEN
  CALL nc_ereport(RoutineName,cmessage='NF90_DEF_VAR_FLETCHER32', &
                  UNIT=um_file % UNIT,nc_status=status_nc)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_def_var_fletcher32

!----------------------------------------------------------------------
! Subroutine: nc_ereport
! Error handling routine for netCDF code
!----------------------------------------------------------------------

SUBROUTINE nc_ereport(RoutineName_in,errorstatus,cmessage,UNIT,nc_status)

USE netcdf, ONLY: nf90_strerror

IMPLICIT NONE

CHARACTER(LEN=*),INTENT(IN)  :: RoutineName_in
INTEGER, INTENT(IN),                                                       &
    OPTIONAL                 :: errorstatus ! code to pass through to ereport
CHARACTER(LEN=*),INTENT(IN),                                               &
    OPTIONAL                 :: cmessage
INTEGER, INTENT(IN),                                                       &
    OPTIONAL                 :: UNIT        ! unit responsible for the error
INTEGER(integer32), INTENT(IN),                                            &
    OPTIONAL                 :: nc_status   ! netCDF library error status

INTEGER                      :: code

CHARACTER(LEN=*), PARAMETER  :: RoutineName='NC_EREPORT'

REAL(KIND=jprb)              :: zhook_handle

TYPE(um_file_type), POINTER  :: um_file

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A,A)')                                                   &
    '****************** NetCDF_File Error Report ',                        &
    '***************************'
CALL umPrint(umMessage,src='umnetcdf_mod')

IF (PRESENT(UNIT)) THEN

  NULLIFY(um_file) 
  um_file => get_file_by_unit(UNIT,handler="netcdf")

  WRITE(umMessage,'(A,I5,A,A)')'Problem with unit ',UNIT,' filename is ',  &
                           TRIM(um_file % filename)
  CALL umPrint(umMessage,src='umnetcdf_mod')

END IF

IF (PRESENT(nc_status)) THEN
  ! Make error code positive for ereport.
  ! NetCDF library error codes are negative, system error codes are positive
  code = ABS(nc_status)
  IF (PRESENT(cmessage)) THEN
    CALL ereport(RoutineName_in,code,                                      &
                    TRIM(nf90_strerror(nc_status))//' : '//TRIM(cmessage))
  ELSE
    CALL ereport(RoutineName_in,code,TRIM(nf90_strerror(nc_status)))
  END IF
ELSE
  IF (PRESENT(errorstatus)) THEN
    code = errorstatus
  ELSE
    code = 99
  END IF
  IF (PRESENT(cmessage)) THEN
    CALL ereport(RoutineName_in,code,cmessage)
  ELSE
    CALL ereport(RoutineName_in,code,'Error from umnetcdf_mod')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nc_ereport

!----------------------------------------------------------------------
! Function: nc_is_file_open
! Check if netCDF file is open using netCDF API routine nf90_inquire
!----------------------------------------------------------------------

LOGICAL FUNCTION nc_is_file_open(um_file) RESULT(lopen)

USE netcdf, ONLY: nf90_inquire,nf90_noerr

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

INTEGER                       :: idx
INTEGER(integer32)            :: ncid
CHARACTER(LEN=*), PARAMETER   :: RoutineName='NC_IS_FILE_OPEN'

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ncid = um_file % nc_meta % ncid

lopen = nf90_inquire(ncid) == nf90_noerr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION nc_is_file_open

!----------------------------------------------------------------------
! Subroutine: set_varid_var
! Set varid array value for a given STASH index
!----------------------------------------------------------------------

SUBROUTINE set_varid_val(stash_index,varid)

IMPLICIT NONE

INTEGER, INTENT(IN)               :: stash_index
INTEGER, INTENT(IN)               :: varid

INTEGER                           :: errorstatus
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='SET_VARID_VAL'

REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (stash_index < 1 .OR. stash_index > ndiag_um) THEN
  WRITE(cmessage,'(A27,I4,A13)') &
       'set_varid_val: stash_index ',stash_index,' out of range'
  errorstatus = 1
  CALL ereport(RoutineName,errorstatus,cmessage)
END IF

varid_val(stash_index) = varid

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_varid_val

!----------------------------------------------------------------------
! Subroutine: get_varid_val
! Get varid array value from a given STASH index
!----------------------------------------------------------------------

SUBROUTINE get_varid_val(stash_index,varid)

IMPLICIT NONE

INTEGER, INTENT(IN)               :: stash_index
INTEGER, INTENT(OUT)              :: varid

INTEGER                           :: errorstatus
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='GET_VARID_VAL'

REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

varid = -1

IF (stash_index < 1 .OR. stash_index > ndiag_um) THEN
  WRITE(cmessage,'(A27,I4,A13)') &
       'get_varid_val: stash_index ',stash_index,' out of range'
  errorstatus = 1
  CALL ereport(RoutineName,errorstatus,cmessage)
END IF

varid = varid_val(stash_index)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE get_varid_val

!----------------------------------------------------------------------
! Subroutine: check_cf_name
! Check string corresponds to the CF naming convention for 
! dimension,variable and attribute names
!----------------------------------------------------------------------

SUBROUTINE check_cf_name(cf_name)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(INOUT) :: cf_name

INTEGER :: idx         ! Location of cf_name character 
                       ! in allowed/disallowed1 string
INTEGER :: i           ! Loop index
INTEGER :: errorstatus ! Error status for ereport

CHARACTER(LEN=*), PARAMETER   :: RoutineName='CHECK_CF_NAME'
! Characters allowed in names according to CF standard
CHARACTER(LEN=*), PARAMETER   :: allowed = 'abcdefghijklmnopqrstuvwxyz'// &
                                           'ABCDEFGHIJKLMNOPQRSTUVWXYZ'// &
                                           '_0123456789'
! Characters from "allowed" not allowed to be the first character in names
CHARACTER(LEN=*), PARAMETER   :: disallowed1 = '_0123456789'
CHARACTER(LEN=errormessagelength) :: cmessage

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check all characters are in the allowed string
DO i=1,LEN_TRIM(cf_name)
  idx = INDEX(allowed,cf_name(i:i))
  IF (idx == 0) THEN
    ! cf_name character not in allowed character string
    WRITE(cmessage,'(A)') &
      'Invalid character '''//cf_name(i:i)// &
      ''' in dimension name '//TRIM(cf_name)
    errorstatus = 1
    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF
END DO

! Check first character is not in the disallowed string
idx = INDEX(disallowed1,cf_name(1:1))
IF (idx > 0) THEN
  ! First cf_name character is in disallowed1 character string
  WRITE(cmessage,'(A)') &
    'Invalid first character '''//cf_name(1:1)// &
    ''' in dimension name '//TRIM(cf_name)
  errorstatus = 2
  CALL ereport(RoutineName,errorstatus,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_cf_name

!----------------------------------------------------------------------
! Subroutine: nc_set_data_type
! Select netCDF real precison for dimensions and variables 
! depending on namelist inputs
!----------------------------------------------------------------------

SUBROUTINE nc_set_data_type()

USE netcdf, ONLY: nf90_float,nf90_double,nf90_int

IMPLICIT NONE

IF (ncdim_prec == ncprec_32bit) THEN
  nc_dim_type_real = nf90_float
ELSE
  nc_dim_type_real = nf90_double
END IF
IF (ncvar_prec == ncprec_32bit) THEN
  nc_var_type_real = nf90_float
ELSE
  nc_var_type_real = nf90_double
END IF
nc_type_int = nf90_int

END SUBROUTINE nc_set_data_type

END MODULE umnetcdf_mod
