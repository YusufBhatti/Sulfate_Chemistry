! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module structure and procedures for EasyAerosol prescriptions
!
! Method:
!   Declares the EasyAerosol data structure where the EasyAerosol
!   distributions and other relevant data are stored. EasyAerosol
!   distributions are read from netCDF files.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE easyaerosol_mod

USE filenamelength_mod, ONLY: &
    filenamelength

IMPLICIT NONE

! EasyAerosol Data structure
TYPE, PUBLIC :: easyaerosol_struct
  CHARACTER (LEN=filenamelength)  :: file_name    ! Name of source file
  CHARACTER (LEN=80)   :: var_name     ! Name of variable in file
  INTEGER              :: varid        ! ID of variable in file
  CHARACTER (LEN=256)  :: std_name     ! standard_name attrib in NetCDF files
  CHARACTER (LEN=256)  :: long_name    ! long_name     attrib in NetCDF files
  CHARACTER (LEN=30)   :: units        ! Units of the field

  INTEGER              :: update_freq  ! Update frequency (hours)
  INTEGER              :: update_type  ! 1 serial, 2 periodic, ...
  INTEGER              :: last_update  ! Num anc update intervals at last update
  LOGICAL              :: l_update     ! True if field updated in a given tstep

  INTEGER              :: ndims        ! Number of dimensions (excluding time)
                                       ! (default: 3)

  INTEGER              :: n_specbands  ! Additional dimension for optical data  
                                       ! that is defined on spectral bands
                                       ! (default: 0 meaning none)

  REAL, ALLOCATABLE:: values_3d(:,:,:)   ! EasyAerosol 3D data
  REAL, ALLOCATABLE:: values_4d(:,:,:,:) ! EasyAerosol 4D data

END TYPE easyaerosol_struct

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EASYAEROSOL_MOD'

CONTAINS

SUBROUTINE set_easyaerosol_field(dimen, n_band, atm, aer, &
                                 col_list, row_list, &
                                 easyaerosol_rad)

  !
  ! Modules used
  !
  USE def_dimen,       ONLY: StrDim
  USE def_atm,         ONLY: StrAtm
  USE def_aer,         ONLY: StrAer
  USE def_easyaerosol, ONLY: t_easyaerosol_rad 

  USE rad_input_mod,   ONLY: l_extra_top

  USE yomhook,         ONLY: lhook, dr_hook
  USE parkind1,        ONLY: jprb, jpim

  IMPLICIT NONE

  ! 
  ! Arguments
  !
  ! Dimensions
  TYPE(StrDim), INTENT(IN) :: dimen

  ! Number of wavebands
  INTEGER, INTENT(IN) :: n_band

  ! Atmospheric properties
  TYPE(StrAtm), INTENT(IN) :: atm

  ! Aerosol properties
  TYPE (StrAer), INTENT(INOUT) :: aer

  ! List of column and row indices on the full grid
  INTEGER, INTENT(IN) :: col_list(dimen%nd_profile)
  INTEGER, INTENT(IN) :: row_list(dimen%nd_profile)

  ! EasyAerosol optical properties
  TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_rad

  !
  ! Local variables
  !
  INTEGER i, j, k
  INTEGER i_top_copy

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
   
  CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_EASYAEROSOL_FIELD'
 
  aer%n_mode = aer%n_mode + 1

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_in, zhook_handle)

  IF (l_extra_top) THEN
    i_top_copy = 2
  ELSE
    i_top_copy = 1
  END IF
    
  !
  ! Set EasyAerosol mass-mixing ratio to 1 so that the
  ! multiplication in opt_prop_ukca_aerosol() resolves
  ! to the EasyAerosol optical properties, unscaled.
  !
  DO j = i_top_copy, atm%n_layer

    DO i = 1, atm%n_profile

      aer%mode_mix_ratio(i, j, aer%n_mode) = 1.0

    END DO ! i

  END DO ! j

  !
  ! Copy the waveband-dependent optical properties,
  ! converting absorption and scattering coefficients
  ! from m-1 to m2 kg-1.
  !
  DO k = 1, n_band

    DO j = i_top_copy, atm%n_layer 
      
      DO i = 1, atm%n_profile
        
        aer%mode_absorption(i, j, aer%n_mode, k) = &
          easyaerosol_rad%absorption(col_list(i), row_list(i), &
                                       atm%n_layer+1-j, k) / &
          atm%density(i, j) 
        aer%mode_scattering(i, j, aer%n_mode, k) = &
          MAX(easyaerosol_rad%extinction(col_list(i), row_list(i), &
                                         atm%n_layer+1-j, k) - &
              easyaerosol_rad%absorption(col_list(i), row_list(i), &
                                         atm%n_layer+1-j, k), 0.0) / &
          atm%density(i, j) 
        aer%mode_asymmetry(i, j, aer%n_mode, k)  = &
          easyaerosol_rad%asymmetry(col_list(i), row_list(i), &
                                    atm%n_layer+1-j, k) 

      END DO ! i
    
    END DO ! j
  
  END DO ! k

  !
  ! If using an extra top layer, set its EasyAerosol properties
  ! to zero.
  ! 
  IF (l_extra_top) THEN

    j = 1

    DO i = 1, atm%n_profile

      aer%mode_mix_ratio(i, j, aer%n_mode) = 0.0

    END DO ! i

    DO k = 1, n_band
      
      DO i = 1, atm%n_profile

        aer%mode_absorption(i, j, aer%n_mode, k) = 0.0
        aer%mode_scattering(i, j, aer%n_mode, k) = 0.0
        aer%mode_asymmetry(i, j, aer%n_mode, k)  = 0.0

      END DO ! i

    END DO ! k

  END IF ! l_extra_top

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_out, zhook_handle)
 
  RETURN
END SUBROUTINE set_easyaerosol_field

END MODULE easyaerosol_mod
