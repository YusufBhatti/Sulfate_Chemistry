! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! This module provides the low level swap bounds routines which use
! local buffers and the point-to-point techinque of MPI.
!
! ATTENTION: these routines should not be called directly; use the
! top-level routines in halo_exchange and non_blocking_halo_exchange
! modules.


MODULE halo_exchange_mpi_mod

USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim

USE um_types, ONLY:                  &
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
  g_datastart,                       &
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
  mpl_tag

USE halo_exchange_base_mod, ONLY:    & 
  halos_comm,                        &
  he_neighbour

USE tools_halo_exchange_mod, ONLY:   &
  copy_ew_halo_to_buffer,            &
  copy_ns_halo_to_buffer,            &
  copy_cr_halo_to_buffer,            &
  copy_buffer_to_ew_halo,            &
  copy_buffer_to_ns_halo,            &
  copy_buffer_to_cr_halo,            &
  copy_ns_halo_single,               &
  copy_edge_to_ns_halo,              &
  normal_halo,                       &
  full_halo 

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


USE ereport_mod,   ONLY: ereport

IMPLICIT NONE

PRIVATE

! LIMITATIONS
!
! The red/black swap bounds routine works only on specifc situations;
! see halo_exchange.F90.


! Begin the swapping of NS halos and then EW halos using buffered MPI
! and wait-any technique.
PUBLIC :: begin_swap_bounds_mpi_nsew_wa
INTERFACE begin_swap_bounds_mpi_nsew_wa
MODULE PROCEDURE begin_swap_bounds_mpi_nsew_wa_dp,  &
                 begin_swap_bounds_mpi_nsew_wa_sp,  &
                 begin_swap_bounds_mpi_nsew_wa_int, &
                 begin_swap_bounds_mpi_nsew_wa_log
END INTERFACE

! End the swapping of NS halos and then EW halos using buffered MPI
! and wait-any technique.
PUBLIC :: end_swap_bounds_mpi_nsew_wa
INTERFACE end_swap_bounds_mpi_nsew_wa
MODULE PROCEDURE end_swap_bounds_mpi_nsew_wa_dp,  &
                 end_swap_bounds_mpi_nsew_wa_sp,  &
                 end_swap_bounds_mpi_nsew_wa_int, &
                 end_swap_bounds_mpi_nsew_wa_log
END INTERFACE


! Swap NS and then EW bounds using buffered MPI and wait-any technique.
PUBLIC :: swap_bounds_mpi_nsew_wa
INTERFACE swap_bounds_mpi_nsew_wa
MODULE PROCEDURE swap_bounds_mpi_nsew_wa_dp,  &
                 swap_bounds_mpi_nsew_wa_sp,  &
                 swap_bounds_mpi_nsew_wa_int, &
                 swap_bounds_mpi_nsew_wa_log
END INTERFACE


! Begin the swapping of the bounds all-at-once using buffered MPI technique.
PUBLIC :: begin_swap_bounds_mpi_aao
INTERFACE begin_swap_bounds_mpi_aao
MODULE PROCEDURE begin_swap_bounds_mpi_aao_dp,  &
                 begin_swap_bounds_mpi_aao_sp,  &
                 begin_swap_bounds_mpi_aao_int, &
                 begin_swap_bounds_mpi_aao_log
END INTERFACE


! End the swapping of the bounds all-at-once using buffered MPI technique.
PUBLIC :: end_swap_bounds_mpi_aao
INTERFACE end_swap_bounds_mpi_aao
MODULE PROCEDURE end_swap_bounds_mpi_aao_dp,  &
                 end_swap_bounds_mpi_aao_sp,  &
                 end_swap_bounds_mpi_aao_int, &
                 end_swap_bounds_mpi_aao_log
END INTERFACE


! Swap the bounds all-at-once using buffered MPI technique.
PUBLIC :: swap_bounds_mpi_aao
INTERFACE swap_bounds_mpi_aao
MODULE PROCEDURE swap_bounds_mpi_aao_dp, &
                 swap_bounds_mpi_aao_sp,  &
                 swap_bounds_mpi_aao_int, &
                 swap_bounds_mpi_aao_log
END INTERFACE


! Special red/black swap bound using buffered MPI technique.
PUBLIC :: swap_bounds_mpi_rb
INTERFACE swap_bounds_mpi_rb
MODULE PROCEDURE swap_bounds_mpi_rb_dp, &
                 swap_bounds_mpi_rb_sp
END INTERFACE


INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HALO_EXCHANGE_MPI_MOD'

CONTAINS

!****************************!
!     BEGIN MPI NSWE WA      !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_nsew_wa_dp(      &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_NSEW_WA_DP'

#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE begin_swap_bounds_mpi_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_nsew_wa_sp(      &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_NSEW_WA_SP'

#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE begin_swap_bounds_mpi_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_nsew_wa_int(     &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_NSEW_WA_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE begin_swap_bounds_mpi_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_nsew_wa_log(     &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_NSEW_WA_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE begin_swap_bounds_mpi_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!      END MPI NSWE WA       !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_nsew_wa_dp(        &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_NSEW_WA_DP'

#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE end_swap_bounds_mpi_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_nsew_wa_sp(      &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_NSEW_WA_SP'

#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE end_swap_bounds_mpi_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_nsew_wa_int(       &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_NSEW_WA_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE end_swap_bounds_mpi_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_nsew_wa_log(       &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  ns_halos_recv,                                  &
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_NSEW_WA_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE end_swap_bounds_mpi_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!     SWAP MPI NSWE WA       !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_nsew_wa_dp(      &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_NSEW_WA_DP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE swap_bounds_mpi_nsew_wa_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_nsew_wa_sp(      &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_NSEW_WA_SP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_nsew_wa.h"

RETURN
END SUBROUTINE swap_bounds_mpi_nsew_wa_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_nsew_wa_int(     &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_NSEW_WA_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE swap_bounds_mpi_nsew_wa_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_nsew_wa_log(     &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_NSEW_WA_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_swap_bounds_mpi_nsew_wa.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE swap_bounds_mpi_nsew_wa_log

! ---------------------------------------------------------------------------- !

!****************************!
!       BEGIN MPI AAO        !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_aao_dp(          &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_AAO_DP'

#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE begin_swap_bounds_mpi_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_aao_sp(          &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_AAO_SP'

#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE begin_swap_bounds_mpi_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_aao_int(         &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_AAO_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE begin_swap_bounds_mpi_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_mpi_aao_log(         &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   


IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_MPI_AAO_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_begin_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE begin_swap_bounds_mpi_aao_log

! ---------------------------------------------------------------------------- !

!****************************!
!        END MPI AAO         !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_aao_dp(            &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_AAO_DP'

#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE end_swap_bounds_mpi_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_aao_sp(            &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_AAO_SP'

#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE end_swap_bounds_mpi_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_aao_int(           &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_integer

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_AAO_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE end_swap_bounds_mpi_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_mpi_aao_log(           &
  field, row_length, rows, levels,                & 
  halo_x,halo_y,                                  & 
  field_type, l_vector,                           & 
  do_east_arg, do_west_arg,                       & 
  do_south_arg, do_north_arg,                     & 
  do_corners_arg,                                 &
  recv_buffer_e, recv_buffer_w,                   &  
  recv_buffer_n, recv_buffer_s,                   &  
  recv_buffer_ne, recv_buffer_se,                 &  
  recv_buffer_nw, recv_buffer_sw,                 &  
  nreq_recv, ireq_recv,                           & 
  tag_base)                                   

IMPLICIT NONE

INTEGER, PARAMETER :: mpl_type   = mpl_logical

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_MPI_AAO_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_end_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE end_swap_bounds_mpi_aao_log

! ---------------------------------------------------------------------------- !

!****************************!
!        SWAP MPI AAO        !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_aao_dp(          &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_AAO_DP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE swap_bounds_mpi_aao_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_aao_sp(          &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_AAO_SP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_aao.h"

RETURN
END SUBROUTINE swap_bounds_mpi_aao_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_aao_int(         &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_AAO_INT'

#define   HALO_EXCHANGE_MPI_INTEGER 1
#include "halo_exchange_mpi_mod_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_INTEGER 

RETURN
END SUBROUTINE swap_bounds_mpi_aao_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_aao_log(         &
  field, row_length, rows, levels,          & 
  halo_x,halo_y,                            & 
  field_type, l_vector,                     & 
  do_east_arg, do_west_arg,                 & 
  do_south_arg, do_north_arg,               & 
  do_corners_arg)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_AAO_LOG'

#define   HALO_EXCHANGE_MPI_LOGICAL 1
#include "halo_exchange_mpi_mod_swap_bounds_mpi_aao.h"
#undef    HALO_EXCHANGE_MPI_LOGICAL 

RETURN
END SUBROUTINE swap_bounds_mpi_aao_log

! ---------------------------------------------------------------------------- !

!****************************!
!        SWAP MPI RB         !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_rb_dp(           &
    field, row_length, rows, levels,        & ! field
    red)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_RB_DP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_rb.h"

RETURN
END SUBROUTINE swap_bounds_mpi_rb_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_mpi_rb_sp(           &
    field, row_length, rows, levels,        & ! field
    red)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MPI_RB_SP'

#include "halo_exchange_mpi_mod_swap_bounds_mpi_rb.h"

RETURN
END SUBROUTINE swap_bounds_mpi_rb_sp


END MODULE halo_exchange_mpi_mod
