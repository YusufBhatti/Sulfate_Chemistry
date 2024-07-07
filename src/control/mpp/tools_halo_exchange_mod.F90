! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! This module provides tools to reduce code duplication in swap_bounds
! routines.
!
! - Routines that copy data from halo/borders to buffers and vice-versa

MODULE tools_halo_exchange_mod

USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE um_types,     ONLY: real32, real64
USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

PRIVATE

! There are two types of halos that can be copied to/from buffers. A
! 'normal halo' which does NOT consider the data of the corners and a
! 'full halo' which DOES consider the data of the corners.
INTEGER, PARAMETER, PUBLIC ::  normal_halo = 0
INTEGER, PARAMETER, PUBLIC ::  full_halo   = 1

! Copy data from the East or West halo to a buffer
PUBLIC :: copy_ew_halo_to_buffer
INTERFACE copy_ew_halo_to_buffer
MODULE PROCEDURE copy_ew_halo_to_buffer_dp,  &
                 copy_ew_halo_to_buffer_sp,  &
                 copy_ew_halo_to_buffer_int, &
                 copy_ew_halo_to_buffer_log
END INTERFACE

! Copy data from buffer to the East or West halo
PUBLIC :: copy_buffer_to_ew_halo
INTERFACE copy_buffer_to_ew_halo
MODULE PROCEDURE copy_buffer_to_ew_halo_dp,  &
                 copy_buffer_to_ew_halo_sp,  &
                 copy_buffer_to_ew_halo_int, &
                 copy_buffer_to_ew_halo_log
END INTERFACE

! Copy data from the North or South halo to a buffer
PUBLIC :: copy_ns_halo_to_buffer
INTERFACE copy_ns_halo_to_buffer
MODULE PROCEDURE copy_ns_halo_to_buffer_dp,  &
                 copy_ns_halo_to_buffer_sp,  &
                 copy_ns_halo_to_buffer_int, &
                 copy_ns_halo_to_buffer_log
END INTERFACE

! Copy data from buffer to the North or South halo
PUBLIC :: copy_buffer_to_ns_halo
INTERFACE copy_buffer_to_ns_halo
MODULE PROCEDURE copy_buffer_to_ns_halo_dp,  &
                 copy_buffer_to_ns_halo_sp,  &
                 copy_buffer_to_ns_halo_int, &
                 copy_buffer_to_ns_halo_log
END INTERFACE

! Copy data from the Corner halo to a buffer
PUBLIC :: copy_cr_halo_to_buffer
INTERFACE copy_cr_halo_to_buffer
MODULE PROCEDURE copy_cr_halo_to_buffer_dp,  &
                 copy_cr_halo_to_buffer_sp,  &
                 copy_cr_halo_to_buffer_int, &
                 copy_cr_halo_to_buffer_log
END INTERFACE

! Copy data from buffer to the Corner halo
PUBLIC :: copy_buffer_to_cr_halo
INTERFACE copy_buffer_to_cr_halo
MODULE PROCEDURE copy_buffer_to_cr_halo_dp,  &
                 copy_buffer_to_cr_halo_sp,  &
                 copy_buffer_to_cr_halo_int, &
                 copy_buffer_to_cr_halo_log
END INTERFACE

! Copy data into the North or South halo in the following cases:
! a) Overpole with a single E/W PE 
! a) Not overpole with a single N/S PE 
PUBLIC :: copy_ns_halo_single
INTERFACE copy_ns_halo_single
MODULE PROCEDURE copy_ns_halo_single_dp,     &
                 copy_ns_halo_single_sp,     &
                 copy_ns_halo_single_int,    &
                 copy_ns_halo_single_log
END INTERFACE

! Copy North/South edge data to the North or South halo
PUBLIC :: copy_edge_to_ns_halo
INTERFACE copy_edge_to_ns_halo
MODULE PROCEDURE copy_edge_to_ns_halo_dp,  &
                 copy_edge_to_ns_halo_sp,  &
                 copy_edge_to_ns_halo_int, &
                 copy_edge_to_ns_halo_log
END INTERFACE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TOOLS_HALO_EXCHANGE_MOD'

CONTAINS


!****************************!
!   From EW halo to buffer   !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ew_halo_to_buffer_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EW_HALO_TO_BUFFER_DP'

#include "tools_halo_exchange_mod_copy_ew_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_ew_halo_to_buffer_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ew_halo_to_buffer_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EW_HALO_TO_BUFFER_SP'

#include "tools_halo_exchange_mod_copy_ew_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_ew_halo_to_buffer_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ew_halo_to_buffer_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EW_HALO_TO_BUFFER_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_ew_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_ew_halo_to_buffer_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ew_halo_to_buffer_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EW_HALO_TO_BUFFER_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_ew_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_ew_halo_to_buffer_log


!****************************!
!   From buffer to EW halo   !  
!****************************!


! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ew_halo_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_EW_HALO_DP'

#include "tools_halo_exchange_mod_copy_buffer_to_ew_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_ew_halo_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ew_halo_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_EW_HALO_SP'

#include "tools_halo_exchange_mod_copy_buffer_to_ew_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_ew_halo_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ew_halo_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_EW_HALO_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_buffer_to_ew_halo.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_buffer_to_ew_halo_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ew_halo_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_EW_HALO_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_buffer_to_ew_halo.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_buffer_to_ew_halo_log


!****************************!
!   From NS halo to buffer   !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_to_buffer_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type)
 
IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_TO_BUFFER_DP'

#include "tools_halo_exchange_mod_copy_ns_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_ns_halo_to_buffer_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_to_buffer_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_TO_BUFFER_SP'

#include "tools_halo_exchange_mod_copy_ns_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_ns_halo_to_buffer_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_to_buffer_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_TO_BUFFER_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_ns_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_ns_halo_to_buffer_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_to_buffer_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_TO_BUFFER_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_ns_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_ns_halo_to_buffer_log


!****************************!
!   From buffer to NS halo   !  
!****************************!


! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ns_halo_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type,                                    &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_NS_HALO_DP'

#include "tools_halo_exchange_mod_copy_buffer_to_ns_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_ns_halo_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ns_halo_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type,                                    &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_NS_HALO_SP'

#include "tools_halo_exchange_mod_copy_buffer_to_ns_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_ns_halo_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ns_halo_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type,                                    &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_NS_HALO_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_buffer_to_ns_halo.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_buffer_to_ns_halo_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_ns_halo_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_y, halo_type,                                    &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_NS_HALO_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_buffer_to_ns_halo.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_buffer_to_ns_halo_log


!****************************!
!   From CR halo to buffer   !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_cr_halo_to_buffer_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y)
 
IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_CR_HALO_TO_BUFFER_DP'

#include "tools_halo_exchange_mod_copy_cr_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_cr_halo_to_buffer_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_cr_halo_to_buffer_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_CR_HALO_TO_BUFFER_SP'

#include "tools_halo_exchange_mod_copy_cr_halo_to_buffer.h"

RETURN
END SUBROUTINE copy_cr_halo_to_buffer_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_cr_halo_to_buffer_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_CR_HALO_TO_BUFFER_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_cr_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_cr_halo_to_buffer_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_cr_halo_to_buffer_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_CR_HALO_TO_BUFFER_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_cr_halo_to_buffer.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_cr_halo_to_buffer_log


!****************************!
!   From buffer to CR halo   !  
!****************************!


! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_cr_halo_dp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y,                                       &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_CR_HALO_DP'

#include "tools_halo_exchange_mod_copy_buffer_to_cr_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_cr_halo_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_cr_halo_sp(                             &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y,                                       &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_CR_HALO_SP'

#include "tools_halo_exchange_mod_copy_buffer_to_cr_halo.h"

RETURN
END SUBROUTINE copy_buffer_to_cr_halo_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_cr_halo_int(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y,                                       &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_CR_HALO_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_buffer_to_cr_halo.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_buffer_to_cr_halo_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_buffer_to_cr_halo_log(                            &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    buffer, disp_x, disp_y,                                       &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_BUFFER_TO_CR_HALO_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_buffer_to_cr_halo.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_buffer_to_cr_halo_log


!****************************!
!       NS halo single       !  
!****************************!


! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_single_dp(                                &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    north_off, south_off, overpoles,                              &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_SINGLE_DP'

#include "tools_halo_exchange_mod_copy_ns_halo_single.h"

RETURN
END SUBROUTINE copy_ns_halo_single_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_single_sp(                                &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    north_off, south_off, overpoles,                              &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_SINGLE_SP'

#include "tools_halo_exchange_mod_copy_ns_halo_single.h"

RETURN
END SUBROUTINE copy_ns_halo_single_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_single_int(                               &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    north_off, south_off, overpoles,                              &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_SINGLE_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_ns_halo_single.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_ns_halo_single_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_ns_halo_single_log(                               &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    north_off, south_off, overpoles,                              &
    overpole_north, overpole_south, change_sign)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_NS_HALO_SINGLE_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_ns_halo_single.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_ns_halo_single_log

! ---------------------------------------------------------------------------- !


!****************************!
!   From edge to NS halo   !  
!****************************!


! ---------------------------------------------------------------------------- !

SUBROUTINE copy_edge_to_ns_halo_dp(                               &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    edge_north, edge_south)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EDGE_TO_NS_HALO_DP'

#include "tools_halo_exchange_mod_copy_edge_to_ns_halo.h"

RETURN
END SUBROUTINE copy_edge_to_ns_halo_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_edge_to_ns_halo_sp(                               &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    edge_north, edge_south)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EDGE_TO_NS_HALO_SP'

#include "tools_halo_exchange_mod_copy_edge_to_ns_halo.h"

RETURN
END SUBROUTINE copy_edge_to_ns_halo_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_edge_to_ns_halo_int(                              &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    edge_north, edge_south)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EDGE_TO_NS_HALO_INT'

#define   TOOLS_HALO_EXCHANGE_INTEGER 1
#include "tools_halo_exchange_mod_copy_edge_to_ns_halo.h"
#undef    TOOLS_HALO_EXCHANGE_INTEGER 

RETURN
END SUBROUTINE copy_edge_to_ns_halo_int

! ---------------------------------------------------------------------------- !

SUBROUTINE copy_edge_to_ns_halo_log(                              &
    field, row_length, rows, levels,                              & ! field
    halo_x, halo_y,                                               & ! halos
    edge_north, edge_south)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPY_EDGE_TO_NS_HALO_LOG'

#define   TOOLS_HALO_EXCHANGE_LOGICAL 1
#include "tools_halo_exchange_mod_copy_edge_to_ns_halo.h"
#undef    TOOLS_HALO_EXCHANGE_LOGICAL 

RETURN
END SUBROUTINE copy_edge_to_ns_halo_log


END MODULE tools_halo_exchange_mod
