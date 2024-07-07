! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Provides utility functionality for manipulating decomposition data

MODULE IOS_Geometry_utils


USE IOS_Model_geometry
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IOS_print_mgr, ONLY:                                                    &
    IOS_print,                                                               &
    IOS_Verbosity,                                                           &
    IOS_message

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! A simple derived type for conveying a boundary.
TYPE box
  INTEGER :: n
  INTEGER :: s
  INTEGER :: e
  INTEGER :: w
END TYPE box

CHARACTER (LEN=errormessagelength), PRIVATE :: geom_util_message

! Classification values
INTEGER, PARAMETER  :: not_computed=-3
INTEGER, PARAMETER  :: no_intersection=-2
INTEGER, PARAMETER  :: complete_intersection=-1
! >0 implies the number of 'patches' from the local domain
! =0 implies 'not relevent'

INTEGER, PARAMETER  :: WesterlyZone=1
INTEGER, PARAMETER  :: EasterlyZone=2

! params/vars  for dr_hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IOS_GEOMETRY_UTILS'

CONTAINS


SUBROUTINE printBox(b)

IMPLICIT NONE
TYPE(box), INTENT(IN) :: b
!$OMP CRITICAL(internal_write)
WRITE(IOS_message,'(A,i5,A,i5,A,i5,A,i5,A,i5,A)')                          &
    'Box: S=',b%s,' N=',b%n,' W=',b%w,' E=',b%e,                           &
    '(data_size=',(b%n-b%s+1)*(b%e-b%w+1),')'
!$OMP END CRITICAL(internal_write)
CALL IOS_print(IOS_message,src='ios_geometry_utils')
END SUBROUTINE printBox

! provide coordinates in the local domain on processor 'cpu', corresponding
! to the bit of the local domain included in the global region
! defined by 'box'. This routine ASSUMES that cpu intersects the box_in
! so you need to call classify_subdomain first
SUBROUTINE getLocalSubdomainBounds(cpu,fld,box_in,localSubdomain)

IMPLICIT NONE

TYPE(box), INTENT(IN)        :: box_in   ! A global box
TYPE(box), INTENT(OUT)       :: localSubdomain  ! A local box
INTEGER, INTENT(IN)          :: fld
INTEGER, INTENT(IN)          :: cpu
TYPE(box)                    :: globalSubdomain
TYPE(box)                    :: processorDomain
TYPE(box)                    :: modelDomain
LOGICAL                      :: ew_wrap
INTEGER                      :: East
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETLOCALSUBDOMAINBOUNDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Normalize the box
globalSubdomain=box_in
CALL Normalize_Box(globalSubdomain,fld,ew_wrap)

! Set the whole domain in local coords first
localSubdomain%s=1
localSubdomain%n=size_map(2,fld,cpu)
localSubdomain%w=1
localSubdomain%e=size_map(1,fld,cpu)

! my domain in global coords
processorDomain%s=offset_map(2,fld,cpu)
processorDomain%n=offset_map(2,fld,cpu)+size_map(2,fld,cpu)-1
processorDomain%w=offset_map(1,fld,cpu)
processorDomain%e=offset_map(1,fld,cpu)+size_map(1,fld,cpu)-1

modelDomain%s=1
modelDomain%n=atm_global_points(2,fld)
modelDomain%w=1
modelDomain%e=atm_global_points(1,fld)

! and then pull the sides in

! South
IF (globalSubdomain%s >= processorDomain%s) THEN
  IF (globalSubdomain%s <= processorDomain%n) THEN
    localSubdomain%s=globalSubdomain%s-processorDomain%s+1
  ELSE
    CALL IOS_Ereport('getLocalSubdomainBounds',99,                         &
        'South boundary is North of my domain')
  END IF
END IF

!North
IF (globalSubdomain%n <= processorDomain%n) THEN
  IF (globalSubdomain%n >= processorDomain%s ) THEN
    localSubdomain%n=globalSubdomain%n-processorDomain%s+1
  ELSE
    CALL IOS_Ereport('getLocalSubdomainBounds',99,                         &
        'North boundary is South of my domain')
  END IF
END IF


IF (ew_wrap) THEN
  East=globalSubdomain%e-modelDomain%e

  !East
  IF (East >= processorDomain%w) THEN
    IF (East <= processorDomain%e) THEN
      localSubdomain%e=East-processorDomain%w+1
    END IF
  END IF
  !West
  IF (globalSubdomain%w <= processorDomain%e) THEN
    IF (globalSubdomain%w >= processorDomain%w) THEN
      localSubdomain%w=globalSubdomain%w-processorDomain%w+1
    END IF
  END IF
ELSE
  East=globalSubdomain%e
  IF (East <= processorDomain%e) THEN
    IF (East >= processorDomain%w) THEN
      localSubdomain%e=East-processorDomain%w+1
    ELSE
      CALL printBox(box_in)
      CALL printBox(processorDomain)
      CALL IOS_Ereport('getLocalSubdomainBounds',99,                       &
          'East boundary is west of my domain')
    END IF
  END IF
  !West
  IF (globalSubdomain%w >= processorDomain%w) THEN
    IF (globalSubdomain%w <= processorDomain%e) THEN
      localSubdomain%w=globalSubdomain%w-processorDomain%w+1
    ELSE
      CALL printBox(box_in)
      CALL printBox(processorDomain)
      CALL IOS_Ereport('getLocalSubdomainBounds',99,                       &
          'West boundary is east of my domain')
    END IF
  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE getLocalSubdomainBounds

SUBROUTINE getLocalSubdomainBoundsCyclic(cpu,fld,box_in,localSubdomain,      &
    region)

IMPLICIT NONE

TYPE(box), INTENT(IN)        :: box_in   ! A global box
TYPE(box), INTENT(OUT)       :: localSubdomain  ! A local box
INTEGER, INTENT(IN)          :: fld
INTEGER, INTENT(IN)          :: cpu
INTEGER, INTENT(IN)          :: region
TYPE(box)                    :: globalSubdomain
TYPE(box)                    :: processorDomain
TYPE(box)                    :: modelDomain
LOGICAL                      :: ew_wrap
INTEGER                      :: East
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETLOCALSUBDOMAINBOUNDSCYCLIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Normalize the box
globalSubdomain=box_in
CALL Normalize_Box(globalSubdomain,fld,ew_wrap)

! Set the whole domain in local coords first
localSubdomain%s=1
localSubdomain%n=size_map(2,fld,cpu)
localSubdomain%w=1
localSubdomain%e=size_map(1,fld,cpu)

! my domain in global coords
processorDomain%s=offset_map(2,fld,cpu)
processorDomain%n=offset_map(2,fld,cpu)+size_map(2,fld,cpu)-1
processorDomain%w=offset_map(1,fld,cpu)
processorDomain%e=offset_map(1,fld,cpu)+size_map(1,fld,cpu)-1

modelDomain%s=1
modelDomain%n=atm_global_points(2,fld)
modelDomain%w=1
modelDomain%e=atm_global_points(1,fld)

! and then pull the sides in

! South
IF (globalSubdomain%s >= processorDomain%s) THEN
  IF (globalSubdomain%s <= processorDomain%n) THEN
    localSubdomain%s=globalSubdomain%s-processorDomain%s+1
  ELSE
    CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99,                   &
        'South boundary is North of my domain')
  END IF
END IF

!North
IF (globalSubdomain%n <= processorDomain%n) THEN
  IF (globalSubdomain%n >= processorDomain%s ) THEN
    localSubdomain%n=globalSubdomain%n-processorDomain%s+1
  ELSE
    CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99,                   &
        'North boundary is South of my domain')
  END IF
END IF

! East west wraping is a given for the cyclic case

IF (region==WesterlyZone) THEN
  ! The westerly zone is the easterly bit of the local cpu's
  ! data

  ! East : The  boundary must not move, otherwise
  !        we aren't cyclic at all!

  ! West :
  localSubdomain%w=globalSubdomain%w-processorDomain%w+1
ELSE IF (region==EasterlyZone) THEN
  ! The easterly zone is the westerly bit of the local cpu's
  ! data

  ! West : The  boundary must not move, otherwise
  !        we aren't cyclic at all!

  ! East : compute the on processor column of the easterly edge
  East=globalSubdomain%e-modelDomain%e
  localSubdomain%e=East-processorDomain%w+1
ELSE
  CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99,                     &
      'Bad region: must provide WesterlyZone or EasterlyZone')

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE getLocalSubdomainBoundsCyclic


! Describe how the box, box_in intersects with my local domain,
! the options are 'it doesnt','completely', or some number
! of overlaps (with 1 processor, a global model the most overlaps
! is 4. For simplicity we will assume that only EW wrapping will occur.
FUNCTION classify_subdomain(cpu,fld,box_in) RESULT(r)

IMPLICIT NONE
TYPE(box), INTENT(IN)  :: box_in   ! The domain
INTEGER, INTENT(IN)    :: fld
INTEGER, INTENT(IN)    :: cpu
INTEGER                :: r
LOGICAL                :: ew_wrap  ! Does domain wraps the meridian
TYPE(box)              :: bx       ! Local copy of box_in
TYPE(box)              :: me       ! My processor domain
INTEGER                :: bxWrapEast
REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CLASSIFY_SUBDOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

r=not_computed
bx=box_in

me%n=offset_map(2,fld,cpu) + size_map(2,fld,cpu) - 1
me%s=offset_map(2,fld,cpu)
me%e=offset_map(1,fld,cpu) + size_map(1,fld,cpu) - 1
me%w=offset_map(1,fld,cpu)

CALL Normalize_Box(bx,fld,ew_wrap)

IF (ew_wrap) THEN
  bxWrapEast=bx%e-atm_global_points(1,fld)

  IF (                                                                     &
      bx%s <= me%s .AND. bx%n >= me%n .AND.                                &
      (bx%w <= me%w .OR. bxWrapEast  >= me%e)) THEN
    ! The box encloses this processor
    r=complete_intersection
  ELSE IF (                                                                &
      bx%s > me%n .OR. bx%n < me%s .OR.                                    &
      (bx%w > me%e .AND. bxWrapEast < me%w ) ) THEN
    r=no_intersection
  ELSE IF (                                                                &
      bx%w <= me%e .AND. bxWrapEast >= me%w ) THEN
    r=2
  ELSE
    r=1
  END IF
ELSE ! No wrapping
  IF (                                                                     &
      bx%s <= me%s .AND. bx%n >= me%n .AND.                                &
      bx%w <= me%w .AND. bx%e >= me%e) THEN
    ! The box encloses this processor
    r=complete_intersection
  ELSE IF (                                                                &
      bx%s > me%n .OR. bx%n < me%s .OR.                                    &
      bx%w > me%e .OR. bx%e < me%w ) THEN
    r=no_intersection
  ELSE ! Without wrapping there can be at most one patch
    r=1
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION classify_subdomain

! Put the coordinates in a standard form such that
! the east bound is off-map and east of west and check for
! errors
SUBROUTINE Normalize_Box(bx,fld,ew_wrap)

IMPLICIT NONE

TYPE(box), INTENT(INOUT) :: bx
INTEGER, INTENT(IN)      :: fld
LOGICAL, INTENT(OUT)     :: ew_wrap
REAL(KIND=jprb)          :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NORMALIZE_BOX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ew_wrap=.FALSE.

IF (bx%w > bx%e) THEN
  bx%e=bx%e+atm_global_points(1,fld)
END IF

! If east is off-map, set the flag
IF (bx%e > atm_global_points(1,fld))                                       &
    ew_wrap=.TRUE.

! Check for other wierd conditions (we don't handle these yet)

IF (bx%s > bx%n) THEN
!$OMP CRITICAL(internal_write)
  WRITE(geom_util_message,'(A,I6,A,I6,A)')                                 &
      'Normalize Box: South(',bx%s,                                        &
      ') greater than North(',bx%n,                                        &
      ') condition not supported yet'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
END IF

IF (bx%s < 1) THEN
!$OMP CRITICAL(internal_write)
  WRITE(geom_util_message,'(A,I6,A)')                                      &
      'Normalize Box: South(',bx%s,                                        &
      ') is off map'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
END IF

IF (bx%n > atm_global_points(2,fld)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(geom_util_message,'(A,I6,A,I6,A)')                                 &
      'Normalize Box: North(',bx%n,                                        &
      ') is off map (',atm_global_points(2,fld),')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
END IF

IF (bx%w > atm_global_points(1,fld)) THEN
!$OMP CRITICAL(internal_write)
  WRITE(geom_util_message,'(A,I6,A,I6,A)')                                 &
      'Normalize Box: West(',bx%w,                                         &
      ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
END IF

IF (bx%w < 1) THEN
!$OMP CRITICAL(internal_write)
  WRITE(geom_util_message,'(A,I6,A,I6,A)')                                 &
      'Normalize Box: West(',bx%w,                                         &
      ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
  CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
END IF

IF (ew_wrap) THEN
  IF (bx%e <= atm_global_points(1,fld)) THEN
!$OMP CRITICAL(internal_write)
    WRITE(geom_util_message,'(A,I6,A,I6,A)')                               &
        'Normalize Box: East WRAP(',bx%e,                                  &
        ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
  END IF
  IF (bx%e > 2*atm_global_points(1,fld)) THEN
!$OMP CRITICAL(internal_write)
    WRITE(geom_util_message,'(A,I6,A,I6,A)')                               &
        'Normalize Box: East WRAP(',bx%e,                                  &
        ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
  END IF

ELSE
  IF (bx%e < 1) THEN
!$OMP CRITICAL(internal_write)
    WRITE(geom_util_message,'(A,I6,A,I6,A)')                               &
        'Normalize Box: East(',bx%e,                                       &
        ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
  END IF
  IF (bx%e > atm_global_points(1,fld)) THEN
!$OMP CRITICAL(internal_write)
    WRITE(geom_util_message,'(A,I6,A,I6,A)')                               &
        'Normalize Box: East(',bx%e,                                       &
        ') is off map (',atm_global_points(1,fld),')'
!$OMP END CRITICAL(internal_write)
    CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
  END IF

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE Normalize_Box

END MODULE IOS_Geometry_utils

