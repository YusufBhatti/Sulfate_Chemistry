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
! Description: Module contains short descriptors of a UM field grid type.
!              To be used to form part of the netCDF dimension name.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: NetCDF output
!

MODULE nc_dimension_id_mod

IMPLICIT NONE

! Horizontal grid types, derived from the STASHmaster Grid parameter

! NetCDF dimension name identifier for theta grid points
CHARACTER(LEN=*), PARAMETER :: nc_horiz_id_theta = "t"

! NetCDF dimension name identifier for U points on 'c' grid
CHARACTER(LEN=*), PARAMETER :: nc_horiz_id_u = "cu"

! NetCDF dimension name identifier for V points on 'c' grid
CHARACTER(LEN=*), PARAMETER :: nc_horiz_id_v = "cv"

! NetCDF dimension name identifier for UV points on 'b' grid
CHARACTER(LEN=*), PARAMETER :: nc_horiz_id_uv = "uv"

! NetCDF dimension name identifier for river-routing grid
CHARACTER(LEN=*), PARAMETER :: nc_horiz_id_river = "r"

! Vertical level types, derived from the STASHmaster LevelT parameter

! NetCDF dimension name identifier for model levels
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_model_level = "model_level_number"

! NetCDF dimension name identifiers for model rho levels
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_eta_rho = "eta_rho"
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_zsea_rho = "zsea_rho"
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_C_rho = "C_rho"

! NetCDF dimension name identifiers for model theta levels
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_eta_theta = "eta_theta"
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_zsea_theta = "zsea_theta"
CHARACTER(LEN=*), PARAMETER :: nc_vert_id_C_theta = "C_theta"
 
! Pseudo dimension types, derived from the STASHmaster PseudT parameter

! NetCDF dimension name identifier for surface type id
CHARACTER(LEN=*), PARAMETER :: nc_pseudo_id_surf_type = "surface_type_id"

END MODULE nc_dimension_id_mod
