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
!  Purpose: Add global attributes and grid mapping variable to
!           netCDF file

MODULE init_nc_mod

! External netCDF library module
USE netcdf,               ONLY: nf90_global,nf90_char,nf90_max_name
USE parkind1,             ONLY: jprb,jpim
USE yomhook,              ONLY: lhook,dr_hook
USE file_manager,         ONLY: um_file_type
USE stash_model_mod,      ONLY: h_a_polelat,h_a_polelong
USE planet_constants_mod, ONLY: planet_radius
USE um_version_mod,       ONLY: um_version_char
USE umnetcdf_mod,         ONLY: nc_put_att,nc_create_var
USE cf_metadata_mod,      ONLY: conventions

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_NC_MOD'

CONTAINS

SUBROUTINE init_nc (um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
! Local variables
!
INTEGER           :: varid            ! NetCDF variable id
INTEGER           :: varid_glob       ! NetCDF global attribute variable id
INTEGER           :: data_type        ! Data type of grid_mapping variable

CHARACTER(LEN=*), PARAMETER  :: RoutineName = "INIT_NC"
CHARACTER(LEN=nf90_max_name) :: varname ! NetCDF variable name

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

varid_glob = nf90_global
data_type = nf90_char

! Write some global attributes
CALL nc_put_att(um_file,varid_glob,'Conventions',conventions)
CALL nc_put_att(um_file,varid_glob,'source', &
                'Met Office Unified Model v'//um_version_char)

! write grid_mapping variable
IF (h_a_polelat == 90.0 .AND. h_a_polelong == 0.0) THEN
  varname = 'grid_crs'
  CALL nc_create_var(um_file,varname,data_type,varid)
  CALL nc_put_att(um_file,varid,'grid_mapping_name','latitude_longitude')
  CALL nc_put_att(um_file,varid,'earth_radius',planet_radius)
ELSE
  varname = 'rotated_latitude_longitude'
  CALL nc_create_var(um_file,varname,data_type,varid)
  CALL nc_put_att(um_file,varid,'grid_mapping_name', &
                                 'rotated_latitude_longitude')
  CALL nc_put_att(um_file,varid,'grid_north_pole_latitude',h_a_polelat)
  CALL nc_put_att(um_file,varid,'grid_north_pole_longitude',h_a_polelong)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_nc

END MODULE init_nc_mod
