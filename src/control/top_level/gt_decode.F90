! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Returns information about a given grid type

SUBROUTINE gt_decode(                                             &
                      grid_type,                                  &
                      model_type,content,coverage,domain,cyclic)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE grdtypes_mod, ONLY: gt_unset, gt_atmos, gt_ocean, gt_wave,    &
                        gt_thetamass, gt_velocity, gt_u_c, gt_v_c,&
                        gt_hybrid, gt_river, gt_allpts, gt_land,  &
                        gt_sea, gt_full, gt_zonal, gt_meridional, &
                        gt_ozone, gt_scalar, gt_compressed,       &
                        gt_lbc, gt_nocyclic, gt_optcyclic, gt_cyclic

USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Given an input GRIDTYPE, this routine will return a value
! for each of MODEL_TYPE, CONTENT, COVERAGE, DOMAIN and CYCLIC
! from the gt_* variables defined in grdtypes_mod


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Arguments:

INTEGER ::                                                        &
  grid_type                                                       &
             ! IN  : Grid type to investigate
, model_type                                                      &
             ! OUT : What model type grid is for
, content                                                         &
             ! OUT : What type of data the grid is for
, coverage                                                        &
             ! OUT : What type of points are on the grid
, domain                                                          &
             ! OUT : What subset of points are on the grid
, cyclic     ! OUT : If the grid includes extra cyclic wrap
             !       around points at the start and end of
             !       each row

! Local variables

INTEGER :: max_grid_types
PARAMETER (max_grid_types=65)

INTEGER :: grid_data(5,max_grid_types)
!             Array holding all the descriptions of the different
!             grid types.5 descriptors for each grid type.

INTEGER ::                                                        &
  icode      ! Error code

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage   ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GT_DECODE'

! Definition of the field types
DATA                                                              &
  grid_data(1:5,1)                                                &
   /gt_atmos,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/          &
, grid_data(1:5,2)                                                &
   /gt_atmos,gt_thetamass,gt_land,gt_full,gt_nocyclic/            &
, grid_data(1:5,3)                                                &
   /gt_atmos,gt_thetamass,gt_sea,gt_full,gt_nocyclic/             &
, grid_data(1:5,4)                                                &
   /gt_atmos,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
, grid_data(1:5,5)                                                &
   /gt_atmos,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
, grid_data(1:5,6)  /5*gt_unset/                                  &
, grid_data(1:5,7)  /5*gt_unset/                                  &
, grid_data(1:5,8)  /5*gt_unset/                                  &
, grid_data(1:5,9)  /5*gt_unset/                                  &
, grid_data(1:5,10) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,11)                                               &
   /gt_atmos,gt_velocity,gt_allpts,gt_full,gt_nocyclic/           &
, grid_data(1:5,12)                                               &
   /gt_atmos,gt_velocity,gt_land,gt_full,gt_nocyclic/             &
, grid_data(1:5,13)                                               &
   /gt_atmos,gt_velocity,gt_sea,gt_full,gt_nocyclic/              &
, grid_data(1:5,14)                                               &
   /gt_atmos,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
, grid_data(1:5,15)                                               &
   /gt_atmos,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
, grid_data(1:5,16) /5*gt_unset/                                  &
, grid_data(1:5,17)                                               &
   /gt_atmos,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
, grid_data(1:5,18)                                               &
   /gt_atmos,gt_U_C,gt_allpts,gt_full,gt_nocyclic/                &
, grid_data(1:5,19)                                               &
   /gt_atmos,gt_V_C,gt_allpts,gt_full,gt_nocyclic/                &
, grid_data(1:5,20) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,21)                                               &
   /gt_atmos,gt_thetamass,gt_land,gt_compressed,gt_nocyclic/      &
, grid_data(1:5,22)                                               &
   /gt_atmos,gt_thetamass,gt_allpts,gt_ozone,gt_nocyclic/         &
, grid_data(1:5,23)                                               &
   /gt_atmos,gt_river,gt_allpts,gt_full,gt_nocyclic/              &
, grid_data(1:5,24) /5*gt_unset/                                  &
, grid_data(1:5,25)                                               &
   /gt_atmos,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
, grid_data(1:5,26)                                               &
   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
, grid_data(1:5,27)                                               &
   /gt_atmos,gt_U_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
, grid_data(1:5,28)                                               &
   /gt_atmos,gt_V_C,gt_allpts,gt_LBC,gt_nocyclic/                 &
, grid_data(1:5,29)                                               &
   /gt_atmos,gt_thetamass,gt_allpts,gt_LBC,gt_nocyclic/           &
, grid_data(1:5,30) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,31)                                               &
   /gt_ocean,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/       &
, grid_data(1:5,32)                                               &
   /gt_ocean,gt_velocity,gt_sea,gt_compressed,gt_nocyclic/        &
, grid_data(1:5,33) /5*gt_unset/                                  &
, grid_data(1:5,34) /5*gt_unset/                                  &
, grid_data(1:5,35) /5*gt_unset/                                  &
, grid_data(1:5,36)                                               &
   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_optcyclic/         &
, grid_data(1:5,37)                                               &
   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_optcyclic/          &
, grid_data(1:5,38)                                               &
   /gt_ocean,gt_U_C,gt_allpts,gt_full,gt_optcyclic/               &
, grid_data(1:5,39)                                               &
   /gt_ocean,gt_V_C,gt_allpts,gt_full,gt_optcyclic/               &
, grid_data(1:5,40) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,41)                                               &
   /gt_ocean,gt_thetamass,gt_allpts,gt_full,gt_cyclic/            &
, grid_data(1:5,42)                                               &
   /gt_ocean,gt_velocity,gt_allpts,gt_full,gt_cyclic/             &
, grid_data(1:5,43)                                               &
   /gt_ocean,gt_thetamass,gt_allpts,gt_zonal,gt_nocyclic/         &
, grid_data(1:5,44)                                               &
   /gt_ocean,gt_velocity,gt_allpts,gt_zonal,gt_nocyclic/          &
, grid_data(1:5,45)                                               &
   /gt_ocean,gt_thetamass,gt_allpts,gt_meridional,gt_nocyclic/    &
, grid_data(1:5,46)                                               &
   /gt_ocean,gt_velocity,gt_allpts,gt_meridional,gt_nocyclic/     &
, grid_data(1:5,47)                                               &
   /gt_ocean,gt_hybrid,gt_allpts,gt_scalar,gt_nocyclic/           &
, grid_data(1:5,48) /5*gt_unset/                                  &
, grid_data(1:5,49) /5*gt_unset/                                  &
, grid_data(1:5,50) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,51)                                               &
   /gt_ocean,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/              &
, grid_data(1:5,52) /5*gt_unset/                                  &
, grid_data(1:5,53) /5*gt_unset/                                  &
, grid_data(1:5,54) /5*gt_unset/                                  &
, grid_data(1:5,55) /5*gt_unset/                                  &
, grid_data(1:5,56) /5*gt_unset/                                  &
, grid_data(1:5,57) /5*gt_unset/                                  &
, grid_data(1:5,58) /5*gt_unset/                                  &
, grid_data(1:5,59) /5*gt_unset/

DATA                                                              &
  grid_data(1:5,60)                                               &
   /gt_wave,gt_thetamass,gt_allpts,gt_full,gt_nocyclic/           &
, grid_data(1:5,61) /5*gt_unset/                                  &
, grid_data(1:5,62)                                               &
   /gt_wave,gt_thetamass,gt_sea,gt_compressed,gt_nocyclic/        &
, grid_data(1:5,63) /5*gt_unset/                                  &
, grid_data(1:5,64) /5*gt_unset/                                  &
, grid_data(1:5,65)                                               &
   /gt_wave,gt_hybrid,gt_allpts,gt_LBC,gt_nocyclic/

! And the code

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF ((grid_type  <   1) .OR.                                       &
    (grid_type  >   max_grid_types)) THEN

  icode=1
  WRITE(cmessage,*) 'Invalid GRID_TYPE ',grid_type,               &
                    ' passed to GT_DECODE'

  CALL ereport('GT_DECODE',icode,cmessage)

END IF

model_type = grid_data(1,grid_type)
content    = grid_data(2,grid_type)
coverage   = grid_data(3,grid_type)
domain     = grid_data(4,grid_type)
cyclic     = grid_data(5,grid_type)

IF (model_type  ==  gt_unset) THEN

  icode=2
  WRITE(cmessage,*) 'GRID_TYPE ',grid_type,                       &
                    ' undefined in GT_DECODE'
  CALL ereport('GT_DECODE',icode,cmessage)

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE gt_decode
