! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
SUBROUTINE ice_sheet_mass( land_ice_index,land_ice_points, ntiles,      &
land_points, land_index, fland, snodep_tile,                            &
greenland_volume, antarctica_volume               )

USE cderived_mod,         ONLY: delta_lambda, delta_phi
USE column_rh_mod,        ONLY: r_theta_levels
USE trignometric_mod,     ONLY: cos_theta_latitude,true_latitude
USE conversions_mod,      ONLY: recip_pi_over_180
USE global_2d_sums_mod,   ONLY: global_2d_sums
USE nlsizes_namelist_mod, ONLY: row_length, rows
USE jules_surface_types_mod, ONLY: ice
USE parkind1,             ONLY: jpim, jprb
USE yomhook,              ONLY: lhook, dr_hook

IMPLICIT NONE

!
! Description:
! Calculate the landice masses of the Northern and Southern Hemisphere land-ice
! volumes. This will be split in the ocean model into ice-shelf melt and
! iceberg flux - being used to weight spatial and temporal distributions. This
! will produce a spatially unchanging melt-water distribution but should
! conserve coupled model fresh water
! 
! N.B. This methodology is only applicable for the AO version of the model, be
! independent of model resolution and replaces the need to use the waterfix
! ancillary file. It will be replaced by a more direct methodology when
! interactive ice sheets and ice shelf cavities are operational.
!
! Input varaibles
INTEGER, INTENT(IN) :: land_points              ! Number of land points
INTEGER, INTENT(IN) :: land_ice_points          ! Number of land ice points
INTEGER, INTENT(IN) :: ntiles                   ! Number of land tiles
INTEGER, INTENT(IN) :: land_index(land_points) ! Lookup global points
                                               ! from land points
INTEGER, INTENT(IN) :: land_ice_index(land_ice_points) ! Lookup land points
                                                       ! from land ice points
REAL, INTENT(IN)    :: fland(land_points)       ! Fraction of land in grid box
REAL, INTENT(IN)    :: snodep_tile(land_points,ntiles) ! Snow depth on tiles

! Output variables
REAL, INTENT(OUT)   :: greenland_volume         ! Snow volume over Greenland
REAL, INTENT(OUT)   :: antarctica_volume        ! Snow volume over Antarctica

! Local varaibles
REAL :: icevol(row_length,rows,2)  ! Snow volume on the global grid for
                                   ! Greenland (1) and Antarctica (2)
REAL :: alat                       ! The latitude
REAL :: cell_area                  ! Area of each grid box
REAL :: snowvol                    ! Snow volume (only used in loop)
REAL :: sum_all(2)                 ! Results of global sums.
                                   ! sum_all(1) = Snow volume over Greenland
                                   ! sum_all(2) = Snow volume over Antarctica
                                   
INTEGER :: i , j, l, m ! loop counter

! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ICE_SHEET_MASS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Ultimately we want a single value  of total snow mass to pass through the
! coupler to the ocean.

! Initialise ice volume as we want it to have zero at non land-ice points
icevol(:,:,:)=0.0

DO m=1,land_ice_points
  l = land_ice_index(m)
  j = (land_index(l)-1)/row_length + 1
  i = land_index(l) - (j-1)*row_length
  
  ! Calculate the snow volume in each grid box
  cell_area = r_theta_levels(i,j,0)                       &
              * r_theta_levels(i,j,0)                     &
              * delta_lambda * delta_phi                  &
              * cos_theta_latitude(i,j)
  snowvol = snodep_tile(l,ice)*fland(l)*cell_area
  
  ! Separate into hemispheres
  alat= recip_pi_over_180 * true_latitude(i,j)   
  IF (alat >= 0.0) THEN 
    icevol(i,j,1)=snowvol
  ELSE
    icevol(i,j,2)=snowvol
  END IF
  
END DO 
        
! Conduct global sum from local arrays and output single floats for the
! ice sheet volumes.
CALL global_2d_sums(icevol,row_length,rows,0,0,2,sum_all)
greenland_volume=sum_all(1)
antarctica_volume=sum_all(2)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE
