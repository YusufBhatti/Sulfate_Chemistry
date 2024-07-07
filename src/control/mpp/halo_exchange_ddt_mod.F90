! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! This module provides the low level swap bounds routines which use
! derived data types and the point-to-point techinque of MPI.
!
! ATTENTION: these routines should not be called directly; use the
! top-level routines in halo_exchange and non_blocking_halo_exchange
! modules.


MODULE halo_exchange_ddt_mod

USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim

USE um_types, ONLY:                  &
  umFortranIntegerSize,              &
  umFortranLogicalSize,              &
  real32,                            &
  real64

USE model_domain_mod, ONLY:          &
  mt_global, mt_cyclic_lam

USE um_parcore,   ONLY:              &
  mype

USE um_parparams, ONLY:              &
  bc_cyclic,                         &
  nodomain,                          &
  peast, pwest, pnorth, psouth,      &
  pnorthwest, pnortheast,            &
  psouthwest, psoutheast

USE um_parvars,   ONLY:              &
  nproc_x, nproc_y,                  &
  bound,                             &
  at_extremity,                      &
  sb_model_type

USE field_types,  ONLY:              &
  fld_type_p, fld_type_v

USE mpl, ONLY:                       &
  mpl_real4,                         &
  mpl_real8,                         &
  mpl_logical,                       &
  mpl_integer,                       &
  mpl_status_size,                   &
  mpl_request_null,                  &
  mpl_address_kind,                  &
  mpl_tag

USE halo_exchange_base_mod, ONLY:    & 
  halos_comm,                        &
  he_neighbour

USE tools_halo_exchange_mod, ONLY:   &
  copy_ns_halo_single,               &
  copy_edge_to_ns_halo

USE tags_params, ONLY:               &
  tag_to_n, tag_to_s,                &
  tag_to_w, tag_to_e,                &
  tag_to_n_p, tag_to_s_p,            &
  tag_from_n, tag_from_s,            &
  tag_from_w, tag_from_e,            &
  tag_from_n_p, tag_from_s_p,        &
  tag_to_ne, tag_to_se,              &
  tag_to_sw, tag_to_nw,              &
  tag_to_ne_p, tag_to_se_p,          &
  tag_to_sw_p, tag_to_nw_p,          &
  tag_from_ne, tag_from_se,          &
  tag_from_sw, tag_from_nw,          &
  tag_from_ne_p, tag_from_se_p,      &
  tag_from_sw_p, tag_from_nw_p

USE tools_halo_exchange_mod, ONLY:   &
  normal_halo,                       &
  full_halo 

USE missing_data_mod, ONLY: imdi

USE ereport_mod,   ONLY: ereport

IMPLICIT NONE

PRIVATE

! LIMITATIONS
!
! The swap_bounds_ddt_??? routines are not thread-safe. This mean that
! they need to be called outiside an OpenMP parallel region. If they
! are called inside a parallel region, they need to be protected by a
! MASTER or CRITICAL region.

! Contains all the derived data type variables.
TYPE, PUBLIC :: ddt_type

  INTEGER :: ew_halo_type 
  ! Data type that represents the e/w halo, i.e columns, it can
  ! only be used to access the data on a single level.

  INTEGER :: ew_halo_ext_type
  ! Data type that represents the extended e/w halo, it can be used
  ! to access the data on multiple levels.
  
  INTEGER :: ns_halo_type
  ! Data type that represents the n/s halo, i.e. rows, it can only
  ! be used to access the data on a single level.

  INTEGER :: ns_halo_ext_type
  ! Data type that represents the extended n/s halo, it can
  ! be used to access the data on multiple levels.

  INTEGER :: rev_ns_halo_type          
  ! Data type that represents a reverse n/s halo, i.e. rows accessed
  ! in reverse Y order, it can only be used to access the data on a
  ! single level. This is used for the over poles communication.

  INTEGER :: rev_ns_halo_ext_type
  ! Data type that represents the extended reverse n/s halo, it can
  ! be used to access the data on multiple levels.

  INTEGER :: cr_halo_type
  ! Data type that represents the corner halo, it can
  ! only be used to access the data on a single level.

  INTEGER :: cr_halo_ext_type
  ! Data type that represents the extended corner halo, it can
  ! be used to access the data on multiple levels.

  INTEGER :: rev_cr_halo_type
  ! Data type that represents the reverse corner halo, i.e. rows
  ! accessed in reverse Y order. it can only be used to access the
  ! data on a single level. This is used for the over poles
  ! communication.

  INTEGER :: rev_cr_halo_ext_type
  ! Data type that represents the extended reverse corner halo, it can
  ! be used to access the data on multiple levels.
    
END TYPE ddt_type

! Manage the derive data type structure.
PUBLIC :: free_ddt
PUBLIC :: create_ddt_ew
PUBLIC :: create_ddt_ns
PUBLIC :: create_ddt_cr

! Begin the swapping of NS halos and then EW halos using DDT and
! wait-any technique.
PUBLIC :: begin_swap_bounds_ddt_nsew_wa
INTERFACE begin_swap_bounds_ddt_nsew_wa
MODULE PROCEDURE begin_swap_bounds_ddt_nsew_wa_dp,  &
                 begin_swap_bounds_ddt_nsew_wa_sp,  &
                 begin_swap_bounds_ddt_nsew_wa_int, &
                 begin_swap_bounds_ddt_nsew_wa_log
END INTERFACE

! End the swapping of NS halos and then EW halos using DDT and
! wait-any technique.
PUBLIC :: end_swap_bounds_ddt_nsew_wa
INTERFACE end_swap_bounds_ddt_nsew_wa
MODULE PROCEDURE end_swap_bounds_ddt_nsew_wa_dp,   &
                 end_swap_bounds_ddt_nsew_wa_sp,   &
                 end_swap_bounds_ddt_nsew_wa_int,  &
                 end_swap_bounds_ddt_nsew_wa_log
END INTERFACE

! Swap NS and then EW bounds using DDT and wait-any technique.
PUBLIC :: swap_bounds_ddt_nsew_wa
INTERFACE swap_bounds_ddt_nsew_wa
MODULE PROCEDURE swap_bounds_ddt_nsew_wa_dp,       &
                 swap_bounds_ddt_nsew_wa_sp,       &
                 swap_bounds_ddt_nsew_wa_int,      &
                 swap_bounds_ddt_nsew_wa_log
END INTERFACE

! Begin the swapping of the bounds all-at-once using DDT.
PUBLIC :: begin_swap_bounds_ddt_aao
INTERFACE begin_swap_bounds_ddt_aao
MODULE PROCEDURE begin_swap_bounds_ddt_aao_dp,   &
                 begin_swap_bounds_ddt_aao_sp,   &
                 begin_swap_bounds_ddt_aao_int,  &
                 begin_swap_bounds_ddt_aao_log 
END INTERFACE

! End the swapping of the bounds all-at-once using DDT
PUBLIC :: end_swap_bounds_ddt_aao
INTERFACE end_swap_bounds_ddt_aao
MODULE PROCEDURE end_swap_bounds_ddt_aao_dp,     &
                 end_swap_bounds_ddt_aao_sp,     &
                 end_swap_bounds_ddt_aao_int,    &
                 end_swap_bounds_ddt_aao_log
END INTERFACE

! Swap the bounds all-at-once using DDT.
PUBLIC :: swap_bounds_ddt_aao
INTERFACE swap_bounds_ddt_aao
MODULE PROCEDURE swap_bounds_ddt_aao_dp,         &
                 swap_bounds_ddt_aao_sp,         &
                 swap_bounds_ddt_aao_int,        &
                 swap_bounds_ddt_aao_log
END INTERFACE


INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HALO_EXCHANGE_DDT_MOD'

CONTAINS

!****************************!
!         FREE DDT           !
!****************************!

! ---------------------------------------------------------------------------- !

! Free previously created derived data type
SUBROUTINE free_ddt ( ddt )

IMPLICIT NONE

TYPE(ddt_type), INTENT(INOUT) :: ddt

! Local variables

INTEGER  :: ierr                    ! MPI error code

REAL(KIND=jprb)  :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FREE_DDT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ddt%ew_halo_type /= imdi) THEN
  CALL mpl_type_free(ddt%ew_halo_type, ierr)
END IF

IF (ddt%ew_halo_ext_type /= imdi) THEN
  CALL mpl_type_free(ddt%ew_halo_ext_type, ierr)
END IF

IF (ddt%ns_halo_type /= imdi) THEN
  CALL mpl_type_free(ddt%ns_halo_type, ierr)
END IF

IF (ddt%ns_halo_ext_type /= imdi) THEN
  CALL mpl_type_free(ddt%ns_halo_ext_type, ierr)
END IF

IF (ddt%rev_ns_halo_type /= imdi) THEN
  CALL mpl_type_free(ddt%rev_ns_halo_type, ierr)
END IF

IF (ddt%rev_ns_halo_ext_type /= imdi) THEN
  CALL mpl_type_free(ddt%rev_ns_halo_ext_type, ierr)
END IF

IF (ddt%cr_halo_type /= imdi) THEN
  CALL mpl_type_free(ddt%cr_halo_type, ierr)
END IF

IF (ddt%cr_halo_ext_type /= imdi) THEN
  CALL mpl_type_free(ddt%cr_halo_ext_type, ierr)
END IF

IF (ddt%rev_cr_halo_type /= imdi) THEN
  CALL mpl_type_free(ddt%rev_cr_halo_type, ierr)
END IF

IF (ddt%rev_cr_halo_ext_type /= imdi) THEN
  CALL mpl_type_free(ddt%rev_cr_halo_ext_type, ierr)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE free_ddt

! ---------------------------------------------------------------------------- !

!****************************!
!      CREATE DDT EW         !
!****************************!

! ---------------------------------------------------------------------------- !

! Create the derived data type needed to perform the E-W halos
SUBROUTINE create_ddt_ew (              &
  row_length, rows, levels,             & ! field values
  halo_x,halo_y,                        & ! halos
  mpl_type, byte_width,                 & ! precision values
  halo_type,                            & 
  ddt)

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length ! number of points on a row
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! number of model levels
INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction
INTEGER, INTENT(IN) :: mpl_type   ! Defines the mpl base type
INTEGER, INTENT(IN) :: byte_width ! Defines the field precision in bytes
INTEGER, INTENT(IN) :: halo_type  ! Full or normal halo

TYPE(ddt_type), INTENT(INOUT) :: ddt 

! Local variables

INTEGER :: ierr  ! MPL error code
INTEGER :: error_code

INTEGER ::                   &
  rows_size,                 &                 
  full_row_length,           & ! length of row including halos
  full_rows                    ! number of rows including halo rows

INTEGER(KIND=mpl_address_kind) lb, extent 

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CREATE_DDT_EW'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

full_row_length = row_length + 2*halo_x
full_rows       = rows + 2*halo_y

SELECT CASE(halo_type)
CASE(normal_halo)
  rows_size = rows
CASE(full_halo)
  rows_size = full_rows
CASE DEFAULT
  error_code = 10
  CALL ereport(RoutineName, error_code, 'Unknown halo type')
END SELECT

! - EW -      

! Create a type that represents the ew halo on a level.
CALL mpl_type_vector(rows_size, halo_x,full_row_length, mpl_type, &
  ddt%ew_halo_type, ierr)
CALL mpl_type_commit(ddt%ew_halo_type, ierr)

! Increase the extent of the ew halo type to include also the
! halo_y area. The new type can be used to retrieve ew halos in
! multiple levels. 
CALL mpl_type_create_resized(ddt%ew_halo_type,0,   & 
  (full_row_length)*(full_rows)*byte_width,        & 
  ddt%ew_halo_ext_type,ierr)      
CALL mpl_type_commit(ddt%ew_halo_ext_type, ierr)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE create_ddt_ew

! ---------------------------------------------------------------------------- !

!****************************!
!      CREATE DDT NS         !
!****************************!

! ---------------------------------------------------------------------------- !

! Create the derived data type needed to perform the N-S halos
SUBROUTINE create_ddt_ns (              &
  row_length, rows, levels,             & ! field values
  halo_x,halo_y,                        & ! halos
  mpl_type, byte_width,                 & ! precision values
  halo_type,                            & 
  ddt)

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length ! number of points on a row
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! number of model levels
INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction
INTEGER, INTENT(IN) :: mpl_type   ! Defines the mpl base type
INTEGER, INTENT(IN) :: byte_width ! Defines the field precision in bytes
INTEGER, INTENT(IN) :: halo_type  ! Full or normal halo

TYPE(ddt_type), INTENT(INOUT) :: ddt 

! Local variables

INTEGER :: ierr  ! MPL error code
INTEGER :: error_code

INTEGER ::                   &
  row_length_size,           &                 
  full_row_length,           & ! length of row including halos
  full_rows                    ! number of rows including halo rows

INTEGER(KIND=mpl_address_kind) lb, extent 

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CREATE_DDT_NS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

full_row_length = row_length + 2*halo_x
full_rows       = rows + 2*halo_y

SELECT CASE(halo_type)
CASE(normal_halo)
  row_length_size = row_length
CASE(full_halo)
  row_length_size = full_row_length
CASE DEFAULT
  error_code = 10
  CALL ereport(RoutineName, error_code, 'Unknown halo type')
END SELECT

! - NS -

! Create a type that represents the ns halo on a level.
CALL mpl_type_vector(halo_y,row_length_size,full_row_length, mpl_type,     &
  ddt%ns_halo_type, ierr)
CALL mpl_type_commit(ddt%ns_halo_type, ierr)

! Increase the extent of the ns halo type to include also the halo_y
! area. The new type can be used to retrieve rows in multiple
! levels. 
CALL mpl_type_create_resized(ddt%ns_halo_type,0,   & 
  (full_row_length)*(full_rows)*byte_width,        & 
  ddt%ns_halo_ext_type,ierr)
CALL mpl_type_commit(ddt%ns_halo_ext_type, ierr)

! - NS REVERSE -

! Create a type that represents the reverse ns halo on a level.
! Attention, notice the reverse stride.
CALL mpl_type_vector(halo_y,row_length_size,-full_row_length, mpl_type, &
  ddt%rev_ns_halo_type, ierr)
CALL mpl_type_commit(ddt%rev_ns_halo_type, ierr)

! Get the extent and lower bound of the reverse halo data type.
CALL mpl_type_get_extent(ddt%rev_ns_halo_type, lb, extent, ierr)

! Increase the extent of the reverse ns halo type to include
! also the halo_y area, keep the same lower bound. The new type
! can be used to retrieve rows in multiple levels.
CALL mpl_type_create_resized(ddt%rev_ns_halo_type,lb, & 
  (full_row_length)*(full_rows)*byte_width,           & 
  ddt%rev_ns_halo_ext_type,ierr)
CALL mpl_type_commit(ddt%rev_ns_halo_ext_type, ierr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE create_ddt_ns

! ---------------------------------------------------------------------------- !

!****************************!
!      CREATE DDT CR        !
!****************************!

! ---------------------------------------------------------------------------- !

! Create the derived data type needed to perform the corners halos
SUBROUTINE create_ddt_cr (              &
  row_length, rows, levels,             & ! field values
  halo_x,halo_y,                        & ! halos
  mpl_type, byte_width,                 & ! precision values
  ddt)

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length ! number of points on a row
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! number of model levels
INTEGER, INTENT(IN) :: halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! size of halo in "j" direction
INTEGER, INTENT(IN) :: mpl_type   ! Defines the mpl base type
INTEGER, INTENT(IN) :: byte_width ! Defines the field precision in bytes

TYPE(ddt_type), INTENT(INOUT) :: ddt 

! Local variables

INTEGER :: ierr  ! MPL error code

INTEGER ::                   &
  full_row_length,           & ! length of row including halos
  full_rows                    ! number of rows including halo rows

INTEGER(KIND=mpl_address_kind) lb, extent 

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CREATE_DDT_CR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

full_row_length = row_length + 2*halo_x
full_rows       = rows + 2*halo_y

! - CORNER -

! Create a type that represents the corner on a level.
CALL mpl_type_vector(halo_y,halo_x,full_row_length, mpl_type,     &
  ddt%cr_halo_type, ierr)
CALL mpl_type_commit(ddt%cr_halo_type, ierr)

! Increase the extent of the corner halo type. The new type can
! be used to retrieve rows in multiple levels.
CALL mpl_type_create_resized(ddt%cr_halo_type,0,   & 
  (full_row_length)*(full_rows)*byte_width,        & 
  ddt%cr_halo_ext_type,ierr)
CALL mpl_type_commit(ddt%cr_halo_ext_type, ierr)

! - CORNER REVERSE -

! Create a type that represents the reverse corner halo on a level.
! Attention, notice the reverse stride.
CALL mpl_type_vector(halo_y,halo_x,-full_row_length, mpl_type, &
  ddt%rev_cr_halo_type, ierr)
CALL mpl_type_commit(ddt%rev_cr_halo_type, ierr)

! Get the extent and lower bound of the reverse halo data type/
CALL mpl_type_get_extent(ddt%rev_cr_halo_type, lb, extent, ierr)

! Increase the extent of the reverse corner halo, keep the same
! lower bound. The new type can be used to retrieve rows in
! multiple levels.
CALL mpl_type_create_resized(ddt%rev_cr_halo_type,lb, & 
  (full_row_length)*(full_rows)*byte_width,           & 
  ddt%rev_cr_halo_ext_type,ierr)
CALL mpl_type_commit(ddt%rev_cr_halo_ext_type, ierr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE create_ddt_cr

! ---------------------------------------------------------------------------- !

!****************************!
!     BEGIN DDT NSEW WA      !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_nsew_wa_dp(      &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_NSEW_WA_DP'

#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE begin_swap_bounds_ddt_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_nsew_wa_sp(      &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_NSEW_WA_SP'

#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE begin_swap_bounds_ddt_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_nsew_wa_int(     &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_NSEW_WA_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_INTEGER 

RETURN
END SUBROUTINE begin_swap_bounds_ddt_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_nsew_wa_log(     &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_NSEW_WA_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE begin_swap_bounds_ddt_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!      END DDT NSEW WA       !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_nsew_wa_dp(        &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_NSEW_WA_DP'

#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE end_swap_bounds_ddt_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_nsew_wa_sp(        &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_NSEW_WA_SP'

#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE end_swap_bounds_ddt_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_nsew_wa_int(       &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_NSEW_WA_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_INTEGER 

RETURN
END SUBROUTINE end_swap_bounds_ddt_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_nsew_wa_log(       &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_NSEW_WA_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE end_swap_bounds_ddt_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!      SWAP DDT NSEW WA      !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_nsew_wa_dp(      &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8
INTEGER, PARAMETER :: byte_width = 8

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_NSEW_WA_DP'

#include "halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE swap_bounds_ddt_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_nsew_wa_sp(      &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4
INTEGER, PARAMETER :: byte_width = 4

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_NSEW_WA_SP'

#include "halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h"

RETURN
END SUBROUTINE swap_bounds_ddt_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_nsew_wa_int(     &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer
INTEGER            :: byte_width ! Calculated within the include file

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_NSEW_WA_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_INTEGER

RETURN
END SUBROUTINE swap_bounds_ddt_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_nsew_wa_log(     &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical
INTEGER            :: byte_width ! Calculated within the include file

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_NSEW_WA_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_swap_bounds_ddt_nsew_wa.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE swap_bounds_ddt_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!       BEGIN DDT AAO        !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_aao_dp(          &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_AAO_DP'

#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE begin_swap_bounds_ddt_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_aao_sp(          &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_AAO_SP'

#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE begin_swap_bounds_ddt_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_aao_int(         &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_AAO_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_INTEGER 

RETURN
END SUBROUTINE begin_swap_bounds_ddt_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_ddt_aao_log(         &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_DDT_AAO_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_begin_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE begin_swap_bounds_ddt_aao_log

! ---------------------------------------------------------------------------- !

!****************************!
!         END DDT AAO        !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_aao_dp(            &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_AAO_DP'

#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE end_swap_bounds_ddt_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_aao_sp(            &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_AAO_SP'

#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE end_swap_bounds_ddt_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_aao_int(           &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_AAO_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_INTEGER

RETURN
END SUBROUTINE end_swap_bounds_ddt_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_ddt_aao_log(           &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  ddt,                                            &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_DDT_AAO_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_end_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE end_swap_bounds_ddt_aao_log

! ---------------------------------------------------------------------------- !

!****************************!
!        SWAP DDT AAO        !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_aao_dp(          &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8
INTEGER, PARAMETER :: byte_width = 8

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_AAO_DP'

#include "halo_exchange_ddt_mod_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE swap_bounds_ddt_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_aao_sp(          &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4
INTEGER, PARAMETER :: byte_width = 4

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_AAO_SP'

#include "halo_exchange_ddt_mod_swap_bounds_ddt_aao.h"

RETURN
END SUBROUTINE swap_bounds_ddt_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_aao_int(         &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer
INTEGER            :: byte_width ! Calculated within the include file

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_AAO_INT'

#define   HALO_EXCHANGE_DDT_INTEGER 1
#include "halo_exchange_ddt_mod_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_INTEGER

RETURN
END SUBROUTINE swap_bounds_ddt_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_ddt_aao_log(         &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical
INTEGER            :: byte_width ! Calculated within the include file

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_DDT_AAO_LOG'

#define   HALO_EXCHANGE_DDT_LOGICAL 1
#include "halo_exchange_ddt_mod_swap_bounds_ddt_aao.h"
#undef    HALO_EXCHANGE_DDT_LOGICAL

RETURN
END SUBROUTINE swap_bounds_ddt_aao_log

! ---------------------------------------------------------------------------- !


END MODULE halo_exchange_ddt_mod
