! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

MODULE d1_array_mod

IMPLICIT NONE

! Information for accessing D1 addressing array including the number of 
! items of info needed for each object and maximum number of objects in D1

INTEGER :: alt_n_submodel_partition
INTEGER, PARAMETER :: alt_n_submodel_partition_max=1

! Number of items of information in D1 addressing array
INTEGER, PARAMETER :: d1_list_len = 17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

! Prognostic, Diagnostic, Secondary or other
INTEGER, PARAMETER :: d1_object_type    = 1 ! Internal model id
INTEGER, PARAMETER :: d1_imodl          = 2 ! Internal model id
INTEGER, PARAMETER :: d1_section        = 3 ! Section
INTEGER, PARAMETER :: d1_item           = 4 ! Item
INTEGER, PARAMETER :: d1_address        = 5 ! Address in D1
INTEGER, PARAMETER :: d1_length         = 6 ! Record length
INTEGER, PARAMETER :: d1_grid_type      = 7 ! Grid type
INTEGER, PARAMETER :: d1_no_levels      = 8 ! Number of levels

! Stash list number for diags. -1 for progs
INTEGER, PARAMETER :: d1_stlist_no      = 9

! Pointer to dump header lookup table
INTEGER, PARAMETER :: d1_lookup_ptr     = 10
INTEGER, PARAMETER :: d1_north_code     = 11 ! Northern row
INTEGER, PARAMETER :: d1_south_code     = 12 ! Southern row
INTEGER, PARAMETER :: d1_east_code      = 13 ! Eastern row
INTEGER, PARAMETER :: d1_west_code      = 14 ! Western row
INTEGER, PARAMETER :: d1_gridpoint_code = 15 ! gridpoint info
INTEGER, PARAMETER :: d1_proc_no_code   = 16 ! Processing Code
INTEGER, PARAMETER :: d1_halo_type      = 17 ! Halo width type

! Types of items for d1_type
INTEGER, PARAMETER :: prognostic = 0
INTEGER, PARAMETER :: diagnostic = 1
INTEGER, PARAMETER :: secondary  = 2
INTEGER, PARAMETER :: other      = 3

! Main data array
REAL, ALLOCATABLE    ::  d1 (:)
! D1 addressing array and number of objects in each submodel
INTEGER, ALLOCATABLE :: d1_addr (:,:,:)

INTEGER :: no_obj_d1 (alt_n_submodel_partition_max)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'D1_ARRAY_MOD'

CONTAINS

SUBROUTINE allocate_d1_array()

USE nlsizes_namelist_mod, ONLY: len_tot, n_obj_d1_max

USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, newline

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_D1_ARRAY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( d1 (len_tot), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
   CALL ereport( modulename//":"//routinename, icode                   &
                 , 'Failure in allocating array [d1]'//                &
                 newline// TRIM(umMessage) )
END IF

ALLOCATE( d1_addr (d1_list_len,n_obj_d1_max,alt_n_submodel_partition), &
          STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
   CALL ereport( modulename//":"//routinename, icode                   &
                 , 'Failure in allocating array [d1_addr]'//           &
                 newline// TRIM(umMessage) )
END IF



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE allocate_d1_array

END MODULE d1_array_mod
