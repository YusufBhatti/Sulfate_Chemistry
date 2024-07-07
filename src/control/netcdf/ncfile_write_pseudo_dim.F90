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
!           and auxiliary coordinate variables,
!           plus associated attributes needed for pseudo coordinates.

MODULE ncfile_write_pseudo_dim_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,           ONLY: um_file_type
USE cstash_mod,             ONLY: idom_b,dompro
USE profilename_length_mod, ONLY: profilename_length
USE stash_array_mod,        ONLY: stlist,stash_pseudo_levels
USE stparam_mod,            ONLY: st_diag_address,st_pseudo_out
USE nc_dimension_id_mod,    ONLY: nc_pseudo_id_surf_type
USE umnetcdf_mod,           ONLY: &
    nc_create_dim,nc_create_var,nc_put_att,check_cf_name, &
    nc_type_int,nc_max_attlen

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_PSEUDO_DIM_MOD'

CONTAINS

SUBROUTINE ncfile_write_pseudo_dim(um_file,pseudo_level_code,ix, &
                                   coordinates,dimid)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN) :: pseudo_level_code        ! Pseudo dimension level code
INTEGER, INTENT(IN) :: ix                       ! STASH item index
CHARACTER(LEN=*), INTENT(INOUT) :: coordinates  ! CF coordinates attribute
INTEGER, INTENT(OUT) :: dimid                   ! NetCDF dimension ID

!
!  Local variables
!
INTEGER :: idiag              ! STASH diagnostic address
INTEGER :: nlevs              ! number of pseudo dimension levels
INTEGER :: dimid_array(1)     ! used to convert dimid value into 1-D array
INTEGER :: varid              ! variable ID for pseudo coordinates

LOGICAL :: lexists            ! TRUE if pseudo dimension already exists

CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'NCFILE_WRITE_PSEUDO_DIM'
CHARACTER(LEN=profilename_length) :: domname      ! Domain profile name
CHARACTER(LEN=nf90_max_name)      :: dimname      ! Dimension name
CHARACTER(LEN=nf90_max_name)      :: varname      ! Auxiliary coordinate 
                                                  ! variable name
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

IF (printstatus >= prdiag) THEN
  CALL umPrint(TRIM(RoutineName)//' variables:',src=RoutineName)
  WRITE(umMessage,'(A,2(I4))') &
            'pseudo_level_code = ',pseudo_level_code
END IF

! Get number of levels for this pseudo dimension variable
IF (pseudo_level_code > 0) THEN
  nlevs = stash_pseudo_levels(1,stlist(st_pseudo_out,ix))
ELSE
  dimid = -1
  ! No Pseudo dimension, so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Default values of pseudo dimension attributes
stdname = ' '
longname = ' '
units = '1'
dimname=TRIM(domname)//'_pseudo'

SELECT CASE (pseudo_level_code)
CASE (1)
  longname = "SW radiation bands"
CASE (2)
  longname = "LW radiation bands"
CASE (3)
  longname = "Atmospheric assimilation groups"
CASE (4)
  longname = "Radiation bands for calculating aerosol optical depth"
CASE (8)
  longname = "HadCM2 Sulphate Loading Pattern Index"
CASE (9)
  longname = "Land and Vegetation Surface types"
  IF (TRIM(coordinates) == "") THEN
    coordinates = TRIM(domname)//'_'//nc_pseudo_id_surf_type
  ELSE
    coordinates = TRIM(coordinates)//' '// &
                  TRIM(domname)//'_'//nc_pseudo_id_surf_type
  END IF
CASE (10)
  longname = "Sea ice categories"
CASE (11)
  longname = "Number of land surface tiles x maximum number of snow layers"
CASE (12)
  longname = "COSP radar reflectivity bins"
CASE (13)
  longname = "COSP number of hydrometeors"
CASE (14)
  longname = "COSP backscattering ratio bins"
CASE (15)
  longname = "COSP ISCCP tau bins"
CASE (16)
  longname = "COSP max. num. of subcolumns output"
CASE DEFAULT
  WRITE(umMessage,'(A,I3)') 'Unknown pseudo level code: ',pseudo_level_code
  CALL umPrint(umMessage,src=RoutineName)
END SELECT

IF (printstatus >= prdiag) THEN
  WRITE(umMessage,'(9(A),I4)') &
            'dimname = ',TRIM(dimname), &
            ' : stdname = ',TRIM(stdname), &
            ' : longname = ',TRIM(longname), &
            ' : units = ',TRIM(units), &
            ' : nlevs = ',nlevs
  CALL umPrint(umMessage,src=RoutineName)
END IF

! Create or read NetCDF dimension
CALL nc_create_dim(um_file,dimname,nlevs,dimid,lexists)

IF (lexists) THEN
  ! Pseudo dimension already exists so finish
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

dimid_array(1) = dimid

! Create NetCDF dimension variables
CALL nc_create_var(um_file,dimname,nc_type_int,varid,dimid_array)

! Save stash index value for later use
um_file % nc_meta % stash_index_dim(dimid) = ix

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

! Create any auxiliary coordinate variables
SELECT CASE (pseudo_level_code)
CASE (9)
  ! Surface type IDs auxiliary coordinate
  varname = TRIM(domname)//'_'//nc_pseudo_id_surf_type
  CALL nc_create_var(um_file,varname,nc_type_int,varid,dimid_array)
  CALL nc_put_att(um_file,varid,'long_name', &
                  'Land and Vegetation Surface type ID')
  CALL nc_put_att(um_file,varid,'units','1')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_pseudo_dim

END MODULE ncfile_write_pseudo_dim_mod
