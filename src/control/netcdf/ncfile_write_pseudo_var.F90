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
!  Purpose: Writes data for coordinate variables
!  created in ncfile_write_pseudo_dim.

MODULE ncfile_write_pseudo_var_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,           ONLY: um_file_type
USE cstash_mod,             ONLY: idom_b,dompro
USE ppxlook_mod,            ONLY: exppxi
USE submodel_mod,           ONLY: atmos_im
USE profilename_length_mod, ONLY: profilename_length
USE version_mod,            ONLY: npslevp
USE stash_array_mod,        ONLY: stlist,stash_pseudo_levels
USE cppxref_mod,            ONLY: ppx_pt_code
USE land_tile_ids,          ONLY: surface_type_ids
USE stparam_mod,            ONLY: &
    st_diag_address,st_sect_no_code,st_item_code,st_pseudo_out
USE nc_dimension_id_mod,    ONLY: nc_pseudo_id_surf_type
USE umnetcdf_mod,           ONLY: &
    nc_inq_ndims,nc_inq_dim,nc_get_varid,nc_get_att,nc_put_var

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_PSEUDO_VAR_MOD'

CONTAINS

SUBROUTINE ncfile_write_pseudo_var(um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
!  Local variables
!
INTEGER :: ndims                 ! Number of dimensions in NetCDF file
INTEGER :: dimlen                ! Dimensions size
INTEGER :: varid                 ! Variable id of dimension
INTEGER :: ix                    ! STASH item index
INTEGER :: i                     ! NetCDF dimension index
INTEGER :: j                     ! Pointer to output pseudo level list
INTEGER :: k                     ! Pseudo dimension index
INTEGER :: idiag                 ! STASH diagnostic address
INTEGER :: section               ! STASH section number
INTEGER :: item                  ! STASH item number
INTEGER :: start(1)              ! Start position of dimension output data
INTEGER :: count1(1)             ! Number of items of dimension output data
INTEGER :: pseudo_level_code     ! PP pseudo level code
INTEGER :: dim_strlen            ! Length of dimension name string
INTEGER :: surf_type_id(npslevp) ! Surface type ID output array

LOGICAL :: lexists    ! TRUE if selected attribute exists

CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'NCFILE_WRITE_PSEUDO_VAR'
CHARACTER(LEN=profilename_length) :: domname ! Domain profile name
CHARACTER(LEN=nf90_max_name)      :: dimname ! NetCDF dimension name
CHARACTER(LEN=nf90_max_name)      :: varname ! NetCDF variable name
CHARACTER(LEN=1)                  :: axis    ! Dimension axis type

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get numnber of dimensions in NetCDF file
CALL nc_inq_ndims(um_file,ndims)

! Loop over each dimension, i is netCDF dimid
DO i=1,ndims

  ! Get dimension name and length for this dimid
  CALL nc_inq_dim(um_file,i,dimname,dimlen)

  ! Not interested in bounds dimension
  IF (dimname(1:6) == 'bounds') CYCLE

  ! Get id of variable with same name as dimension
  CALL nc_get_varid(um_file,dimname,varid)

  ! Get axis attribute
  CALL nc_get_att(um_file,varid,'axis',axis,lexists)
  ! Pseudo dimensions won't have an axis attribute 
  ! so don't go any further if it does
  IF (lexists) CYCLE

  dim_strlen = LEN_TRIM(dimname)

  ! This routine is for processing pseudo dimensions, so only select those
  IF (dimname(dim_strlen-5:dim_strlen) == 'pseudo') THEN

    ! Get stash index value
    ix = um_file % nc_meta % stash_index_dim(i)

    IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
            'Pseudo dimension name ',TRIM(dimname), &
            ' of size ',dimlen,' on unit ',um_file % UNIT,' found'
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,A,A,I5)') &
                'STASH index for pseudo dimension ',TRIM(dimname),' is ',ix
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Intialise start and count variables needed to write netCDF data
    start(1) = 1
    count1(1) = dimlen

    ! Get STASH and pseudo levels information
    section = stlist(st_sect_no_code,ix)
    item = stlist(st_item_code,ix)
    pseudo_level_code = exppxi(atmos_im,section,item,ppx_pt_code)

    varname = dimname

    ! Get variable id and write pseudo dimension values
    CALL nc_get_varid(um_file,varname,varid)

    IF (pseudo_level_code > 0) THEN
      j = stlist(st_pseudo_out,ix)
      CALL nc_put_var(um_file,varid,stash_pseudo_levels(2:dimlen+1,j), &
                      start,count1)
    END IF

    idiag = stlist(st_diag_address,ix)
    domname = dompro(idom_b(idiag))

    ! Write out any auxiliary coordinate variables
    SELECT CASE (pseudo_level_code)
    CASE (9) ! Land and Vegetation Surface types

      ! Get Surface type IDs
      DO k=1,dimlen
        surf_type_id(k) = surface_type_ids(stash_pseudo_levels(k+1,j))
      END DO

      ! Variable name for auxilliary coordinate
      varname = TRIM(domname)//'_'//nc_pseudo_id_surf_type

      ! Get variable id and write auxilliary coordinates for pseudo dimension
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,surf_type_id,start,count1)

    END SELECT

  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_pseudo_var

END MODULE ncfile_write_pseudo_var_mod
