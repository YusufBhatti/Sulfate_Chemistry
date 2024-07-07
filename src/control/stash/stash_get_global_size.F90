! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculates the "global" size of STASHed data.
!
! Subroutine Interface:
SUBROUTINE stash_get_global_size(                                 &
  global_north_in , global_east_in ,                              &
  global_south_in , global_west_in ,                              &
  levels_in ,                                                     &
  gridpoint_code , processing_code ,                              &
  global_field_size ,                                             &
  icode , cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE field_types, ONLY: fld_type_p
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE stparam_mod, ONLY: st_replace_code, st_accum_code, st_max_code,&
                    st_min_code, st_time_mean_code, vert_mean_base,&
                 zonal_mean_base, merid_mean_base, field_mean_base,&
                 global_mean_base, global_mean_top

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Description:
! Calculates the global (ie. total size on disk) size of a
! STASH request.
!
! Method:
! Using the PROCESSING_CODE to indicate the type of STASH request,
! the GRIDPOINT_CODE to indicate the grid type of the data,
! and the subdomain limits, the total size of the data is calculated.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stash
!
! Subroutine arguments:

INTEGER ::                                                        &
  global_north_in                                                 &
                    ! IN: specification of subdomain boundaries
, global_east_in                                                  &
                    ! IN: ""
, global_south_in                                                 &
                    ! IN: ""
, global_west_in                                                  &
                    ! IN: ""
, levels_in                                                       &
                    ! IN: number of levels
, gridpoint_code                                                  &
                    ! IN: indicates the output grid type
, processing_code   ! IN: indicates the type of STASH processing

INTEGER ::                                                        &
  global_field_size ! OUT: size of STASH data on disk

INTEGER ::                                                        &
  icode             ! OUT: Return code (0=OK)

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage          ! OUT: Error message

! Local variables

INTEGER ::                                                        &
! copies of input arguments, which get modified according the
! type of output grid
        global_north,global_east,global_south,global_west               &
      , levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH_GET_GLOBAL_SIZE'

! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
global_north = global_north_in
global_east  = global_east_in
global_south = global_south_in
global_west  = global_west_in
levels       = levels_in

! Fix wrap-arounds s.t. east > west

IF (global_west  >   global_east)                                 &
  global_east=global_east+glsize(1,fld_type_p)

! Full field or subdomain output:

IF ((processing_code  ==  st_replace_code) .OR.                   &
    (processing_code  ==  st_accum_code) .OR.                     &
    (processing_code  ==  st_time_mean_code) .OR.                 &
    (processing_code  ==  st_max_code) .OR.                       &
    (processing_code  ==  st_min_code)) THEN

  IF ((gridpoint_code  >=  vert_mean_base) .AND.                  &
                                                   ! vertical
      (gridpoint_code  <   zonal_mean_base)) THEN ! mean
    levels=1

  ELSE IF ((gridpoint_code  >=  zonal_mean_base) .AND.             &
                                                       ! zonal
          (gridpoint_code  <   merid_mean_base)) THEN  ! mean
    global_east=global_west

  ELSE IF ((gridpoint_code  >=  merid_mean_base) .AND.             &
                                                      ! merid.
          (gridpoint_code  <   field_mean_base)) THEN ! mean
    global_south=global_north

  ELSE IF ((gridpoint_code  >=  field_mean_base)  .AND.            &
                                                       ! field
          (gridpoint_code  <   global_mean_base)) THEN ! fmean
    global_east=global_west
    global_south=global_north

  ELSE IF (gridpoint_code  >=  global_mean_base) THEN
    global_east=global_west
    global_south=global_north
    levels=1

  ELSE IF (gridpoint_code  >   global_mean_top) THEN
    icode=1
    WRITE(umMessage,*) 'Grid type ',gridpoint_code,                       &
               ' not yet supported by MPP STASH'
    CALL umPrint(umMessage,src='stash_get_global_size')
    cmessage='Unsupported grid type'
    GO TO 9999
  END IF

  global_field_size=(global_east-global_west+1)*                  &
                    (global_north-global_south+1)*                &
                    levels

ELSE
  icode=2
  WRITE(umMessage,*) 'Processing code ',processing_code,                  &
             ' not yet supported by MPP STASH'
  CALL umPrint(umMessage,src='stash_get_global_size')
  cmessage='Unsupported processing code'
  GO TO 9999
END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE stash_get_global_size

