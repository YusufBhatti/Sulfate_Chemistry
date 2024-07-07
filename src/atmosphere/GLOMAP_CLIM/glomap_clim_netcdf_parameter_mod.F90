! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
! Define structures to hold NetCDF FIELD related information
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

MODULE glomap_clim_netcdf_parameter_mod

USE ancil_mod,          ONLY: &
    stash_num_max

USE filenamelength_mod, ONLY: &
    filenamelength

USE items_nml_mod,      ONLY: &
    netcdf_varname_len

IMPLICIT NONE

INTEGER, SAVE      :: ncdf_file_unit_no      ! Unit No for NetCDF files

INTEGER, SAVE      :: number_ncdf_files = 0  ! number of NetCDF files to open

! Expected max No. of values in any attribute, just in case some attributes
! may be arrays (could be changed to higher value if needed)
INTEGER, PARAMETER :: max_att_len    = 10

INTEGER, PARAMETER :: max_dims_per_var = 5   ! Max nr of dimens per variable

INTEGER, PARAMETER :: max_ncdf_files = 20    ! Maximum number of netcdf files

INTEGER, PARAMETER :: stash_name_len = 36    ! Fixed length of stashname

! Length for the strings indicating variable type (e.g. 'CHAR', 'REAL', etc.)
INTEGER, PARAMETER :: vartype_len = 4

! Expected maximum length of the calendar attribute in the time variable
INTEGER, PARAMETER :: t_calendar_len = 20

! Expected maximum length of the units attribute in the time variable
INTEGER, PARAMETER :: t_units_len    = 40

! Structure used to hold NetCDF information from ITEMS.
TYPE netcdf_fileinfo
  CHARACTER (LEN=filenamelength)     :: ncdf_file
  CHARACTER (LEN=netcdf_varname_len) :: varname    ( stash_num_max )
  CHARACTER (LEN=stash_name_len)     :: stash_name ( stash_num_max )
  INTEGER                            :: stashcode  ( stash_num_max )
  INTEGER                            :: section    ( stash_num_max )
  INTEGER                            :: item       ( stash_num_max )
  INTEGER                            :: fields_in_file
END TYPE  netcdf_fileinfo

! ncdf_list is an array of type netcdf_fileinfo
TYPE (netcdf_fileinfo)               :: ncdf_list ( max_ncdf_files )

END MODULE glomap_clim_netcdf_parameter_mod
