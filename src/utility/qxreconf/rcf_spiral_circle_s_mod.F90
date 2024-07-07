! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: RCF_SPIRAL_CIRCLE_S--------------------------------------

MODULE rcf_spiral_circle_s_mod

IMPLICIT NONE

! Description:
! Calculate the values for unresolved points.
! Uses a new spiral circle method, which finds the closest point
! by distance (in m)
!
! Method:
! Uses shumlib f_shum_spiral_search_algorithm function
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code description:
! Language: Fortran
! This code is written to UMDP3 standards.

PRIVATE

PUBLIC :: rcf_spiral_circle_s

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SPIRAL_CIRCLE_S_MOD'

CONTAINS
SUBROUTINE rcf_spiral_circle_s                                    &
           (lsm,index_unres,no_point_unres,                       &
            points_phi,points_lambda,lats,lons,                   &
            is_land_field, constrained, cyclic_domain, unres_mask,&
            indices)

USE f_shum_spiral_search_mod, ONLY: &
    f_shum_spiral_search_algorithm
USE planet_constants_mod, ONLY : & 
    planet_radius          ! Earth radius in m
USE ereport_mod,         ONLY : &
    ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook
USE parkind1,  ONLY: &
    jprb,            &
    jpim
USE fort2c_interfaces, ONLY: logical_to_bool

IMPLICIT NONE


! Subroutine arguments

INTEGER, INTENT(IN) :: points_phi ! number of rows in grid
INTEGER, INTENT(IN) :: points_lambda ! number of columns in grid
INTEGER, INTENT(IN) :: no_point_unres ! number of unresolved points
! index to unresolved pts
INTEGER, INTENT(INOUT) :: index_unres(points_lambda*points_phi)
! index of resolved points
INTEGER, INTENT(INOUT) :: indices(points_lambda*points_phi)

REAL, INTENT(IN) :: lats(points_phi) ! field
REAL, INTENT(IN) :: lons(points_lambda) ! field

LOGICAL, INTENT(IN) :: cyclic_domain ! =T if data wraps E/W
LOGICAL, INTENT(IN) :: constrained ! =T if constraint is to be applied
LOGICAL, INTENT(IN) :: lsm(points_lambda*points_phi) ! land sea mask
LOGICAL, INTENT(IN) :: is_land_field ! =F for sea field  =T for land field
! =F for a point that is resolved  =T for an unresolved point
LOGICAL, INTENT(IN) :: Unres_Mask(points_lambda*points_phi)

! LOCAL VARIABLES

! Shorter copies of the inputs to use
INTEGER :: index_unres_tmp(no_point_unres)
INTEGER :: indices_tmp(no_point_unres)

INTEGER :: errorstatus

CHARACTER (Len=*), PARAMETER  :: RoutineName = 'RCF_SPIRAL_CIRCLE_S'
CHARACTER (LEN=errormessagelength) :: cmessage

! Maximum distance when using constrained search (in m)
REAL, PARAMETER :: constrained_max_dist = 200000.0
! Step size coefficient for spiral search
REAL, PARAMETER :: dist_step = 3.0

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header
! Dr Hook is included in this routine instead of spiral_circle_search as that
! routine is used by other codes and should have no dependencies.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Populate short form Shumlib input variable
index_unres_tmp(:) = index_unres(1:no_point_unres)

errorstatus = 0
errorstatus = f_shum_spiral_search_algorithm(                &
               logical_to_bool(lsm),                         &
               index_unres_tmp, no_point_unres, points_phi,  &
               points_lambda, lats, lons,                    &
               logical_to_bool(is_land_field),               &
               logical_to_bool(constrained),                 &
               constrained_max_dist, dist_step,              &
               logical_to_bool(cyclic_domain),               & 
               logical_to_bool(unres_mask),                  &
               indices_tmp, planet_radius, cmessage)

! Populate long form output variable
indices(1:no_point_unres) = indices_tmp(:)

IF ( errorstatus /= 0 ) THEN
    CALL ereport(RoutineName, errorstatus, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_spiral_circle_s
END MODULE rcf_spiral_circle_s_mod
