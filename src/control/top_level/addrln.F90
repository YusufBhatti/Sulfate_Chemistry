! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Determine STASH input length per vertical level for prog var
! Subroutine Interface:

SUBROUTINE addrln(igpl,halo_type,LEN,size_type)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore,   ONLY: mype
USE UM_ParParams, ONLY: local_data, peast, pwest
USE rimtypes
USE lbc_mod
USE grdtypes_mod, ONLY: gt_unset, gt_atmos, gt_ocean, gt_wave,    &
                        gt_thetamass, gt_velocity, gt_u_c, gt_v_c,&
                        gt_hybrid, gt_river, gt_allpts, gt_land,  &
                        gt_sea, gt_full, gt_zonal, gt_meridional, &
                        gt_ozone, gt_scalar, gt_compressed,       &
                        gt_lbc, gt_nocyclic, gt_optcyclic, gt_cyclic
USE cppxref_mod, ONLY: ppx_atm_lbc_orog
USE mpp_conf_mod, ONLY: include_halos_ew, include_halos_ns
USE ozone_inputs_mod, ONLY: zon_av_ozone
USE nlsizes_namelist_mod, ONLY: &
    global_land_field, local_land_field

IMPLICIT NONE
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!

! Subroutine arguments:
INTEGER ::                                                        &
  igpl                                                            &
                     ! IN : grid code
, halo_type                                                       &
                     ! IN : type of halo (none, single
                     !                    or extended)
, size_type                                                       &
                     ! IN : "local_data" or "global_data"
, LEN                ! OUT : length of field

! Local variables

INTEGER ::                                                        &
  ix1                                                             &
                  ! Column number of start of area
, ix2                                                             &
                  ! Column number of end of area
, iy1                                                             &
                  ! Row number at start of area
, iy2                                                             &
                  ! Row number at end of area

, local_IX1                                                       &
                  ! local IX1 for this processor
, local_IX2                                                       &
                  ! local IX2 for this processor
, local_IY1                                                       &
                  ! local IY1 for this processor
, local_IY2                                                       &
                  ! local IY2 for this processor

, fld_type                                                        &
                  ! Which kind of variable (theta u or v?)

                  ! Variables from GT_DECODE:
, model_type                                                      &
                  ! model type of grid
, content                                                         &
                  ! content type of grid
, coverage                                                        &
                  ! coverage of grid
, domain                                                          &
                  ! domain of grid
, cyclic                                                          &
                  ! does grid contain cyclic wrap columns

, rim_type                                                        &
                  ! Which type of rimwidth - normal or orography?
, ocean_extra_pts ! extra points added to ocean row length
!                       ! to allow for the wrap around pts at the
!                       ! start and end of each global row

! Functions
INTEGER ::                                                        &
  get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADDRLN'

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Get information about grid type
! DEPENDS ON: gt_decode
CALL gt_decode(igpl,                                              &
               model_type,content,coverage,domain,cyclic)

! Determine row/column nos. for global domain
! DEPENDS ON: lltorc
CALL lltorc(igpl,90,-90,0,360,iy1,iy2,ix1,ix2)

IF ((size_type  ==  local_data) .AND.                             &
    ((domain  /=  gt_compressed) .AND.                            &
     (domain  /=  gt_LBC))) THEN

  ! Convert the global subdomain limits to local subdomain limits
  ! DEPENDS ON: global_to_local_subdomain
  CALL global_to_local_subdomain( include_halos_ew,include_halos_ns,&
                                  igpl,halo_type,mype,            &
                                  iy1,ix2,iy2,ix1,                &
                                  local_IY1,local_IX2,            &
                                  local_IY2,local_IX1)


  ix1=local_IX1
  ix2=local_IX2
  iy1=local_IY1
  iy2=local_IY2

  ocean_extra_pts=0
  IF (at_extremity(PWest)) ocean_extra_pts=ocean_extra_pts+1
  IF (at_extremity(PEast)) ocean_extra_pts=ocean_extra_pts+1
ELSE
  ocean_extra_pts=2 ! extra pt at start and end of row for
  !                         ! wrap around
END IF

IF (domain  ==  gt_compressed) THEN

  IF (model_type  ==  gt_atmos) THEN

    IF (size_type  ==  local_data) THEN
      LEN=local_land_field
    ELSE
      LEN=global_land_field
    END IF

  ELSE IF (model_type  ==  gt_ocean) THEN

    LEN=-1 ! Set flag at this stage for a multi-level
           ! compress

  END IF ! MODEL_TYPE

ELSE IF (domain  ==  gt_ozone) THEN

  IF (zon_av_ozone) THEN ! zonal

    IF (size_type  ==  local_data) THEN
      LEN=iy2-iy1+1+2*halosize(2,halo_type)
    ELSE
      LEN=iy2-iy1+1
    END IF

  ELSE ! Full fields

    LEN=(ix2-ix1+1)*(iy2-iy1+1)

  END IF

ELSE IF (domain  ==  gt_zonal) THEN
  IF (size_type  ==  local_data) THEN
    LEN=iy2-iy1+1+2*halosize(2,halo_type)
  ELSE
    LEN=iy2-iy1+1
  END IF

ELSE IF (domain  ==  gt_LBC) THEN

  IF (model_type  ==  gt_atmos) THEN

    ! DEPENDS ON: get_fld_type
    fld_type=get_fld_type(igpl)
    IF (igpl  ==  ppx_atm_lbc_orog) THEN
      rim_type=rima_type_orog
    ELSE
      rim_type=rima_type_norm
    END IF

    ! NB - these sizes are just for a single level. The LBCs are stored
    ! with all vertical levels in one field. This number will need to
    ! be multiplied by the number of levels by the calling routine

    IF (size_type  ==  local_data) THEN
      LEN=lenrima(fld_type,halo_type,rim_type)
    ELSE
      LEN=global_LENRIMA(fld_type,halo_type,rim_type)
    END IF

  END IF

ELSE

  LEN =(ix2-ix1+1)*(iy2-iy1+1)

END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE addrln
