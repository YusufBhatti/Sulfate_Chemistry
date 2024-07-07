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
!  Purpose: Write STASH output data in netCDF format
!  for given level and time value

MODULE ncfile_write_data_mod

USE yomhook,      ONLY: lhook,dr_hook
USE parkind1,     ONLY: jprb,jpim
USE file_manager, ONLY: um_file_type
USE cppxref_mod,  ONLY: ppx_type_real
USE umnetcdf_mod, ONLY: nc_time_index,nc_max_var_dims,nc_dim_used, &
                        nc_put_var,get_varid_val

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_DATA_MOD'

CONTAINS

SUBROUTINE ncfile_write_data(um_file,stash_index,data_type_code, &
                             n_cols,n_rows,level,nlevels,npseudo,ppdata)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file
INTEGER, INTENT(IN) :: stash_index    ! STASH index
INTEGER, INTENT(IN) :: data_type_code ! STASH datatype
INTEGER, INTENT(IN) :: n_cols         ! Number of columns for ppdata
INTEGER, INTENT(IN) :: n_rows         ! Number of rows for ppdata
INTEGER, INTENT(IN) :: level          ! Level number for ppdata 
                                      ! (includes pseudo levels if any)
INTEGER, INTENT(IN) :: nlevels        ! Number of levels
                                      ! (includes pseudo levels if any)
INTEGER, INTENT(IN) :: npseudo        ! Number of pseudo dimension levels

REAL, INTENT(IN)    :: ppdata(n_rows*n_cols) ! Output data to be written

!
!  Local variables
!
INTEGER :: level_num                 ! level number part of output level
INTEGER :: pseudo_level_num          ! pseudo level number part of output level
INTEGER :: nlev                      ! Number of vertical levels
INTEGER :: start1(nc_max_var_dims)   ! Start index of output data
                                     ! for each dimension
INTEGER :: count1(nc_max_var_dims)   ! Number of items of output data
                                     ! for each dimension
INTEGER :: ndim                      ! Number of output data dimensions
INTEGER :: varid                     ! Variable id of output data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NCFILE_WRITE_DATA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If first level then this a new time value for output data
IF (level == 1) THEN
  nc_time_index(stash_index) = nc_time_index(stash_index) + 1
END IF

! Extract level number and pseudo level number from level
IF (npseudo == 0) THEN
   level_num = level
   pseudo_level_num = 0
ELSE
   nlev = nlevels/npseudo
   level_num = 1 + MOD(level-1,nlev)
   pseudo_level_num = 1 + (level-1)/nlev
END IF

! Get saved varid from stash_index
CALL get_varid_val(stash_index,varid)

! Set netCDF output's dimension index/size values
ndim = 0
IF (nc_dim_used(1,stash_index)) THEN
  ndim = ndim + 1
  start1(ndim) = 1
  count1(ndim) = n_cols
END IF
IF (nc_dim_used(2,stash_index)) THEN
  ndim = ndim + 1
  start1(ndim) = 1
  count1(ndim) = n_rows
END IF
IF (nc_dim_used(3,stash_index)) THEN
  ndim = ndim + 1
  start1(ndim) = level_num
  count1(ndim) = 1
END IF
IF (nc_dim_used(4,stash_index)) THEN
  ndim = ndim + 1
  start1(ndim) = pseudo_level_num
  count1(ndim) = 1
END IF
IF (nc_dim_used(5,stash_index)) THEN
  ndim = ndim + 1
  start1(ndim) = nc_time_index(stash_index)
  count1(ndim) = 1
END IF

! Write ppdata to netCDF file
IF (data_type_code == ppx_type_real) THEN
  CALL nc_put_var(um_file,varid,ppdata,start1(1:ndim),count1(1:ndim))
ELSE
  CALL nc_put_var(um_file,varid,TRANSFER(ppdata,ndim,n_rows*n_cols), &
                  start1(1:ndim),count1(1:ndim))
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_data

END MODULE ncfile_write_data_mod
