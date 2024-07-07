! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

!
! This module provides methods for updating boundaries (halos)
! on fields, either by extending the field or exchanging data
! with neighbouring processors depending on the model.
!
! Communications take place in private communicators which musch be
! initialised with a call to halo_exchange_init() before other methods
! can be used.
!
! Logic and communicator setup is non-trivial for global models.

MODULE Halo_Exchange

USE mpl, ONLY:       &
    mpl_real,        &
    mpl_real4,       &
    mpl_real8,       &
    mpl_logical,     &
    mpl_status_size, &
    mpl_request_null,&
    mpl_tag,         &
    mpl_undefined,   &
    mpl_address_kind

USE um_types, ONLY:  &
     real32,         &
     real64

USE ereport_mod,  ONLY: ereport

USE missing_data_mod, ONLY: imdi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE mpp_conf_mod, ONLY:             &
    name_sb_method,                 &
    mpi_aao,                        &
    mpi_nsew_wa,                    &
    ddt_aao,                        &
    ddt_nsew_wa,                    &
    swap_bounds_method_3D_small,    &
    swap_bounds_method_3D_large,    &
    swap_bounds_method_2D

USE UM_ParVars
USE UM_ParCore,   ONLY: mype, nproc
USE UM_ParParams, ONLY: bc_cyclic, bc_overpole, nodomain, &
                        peast, pwest, pnorth, psouth, &
                        pnortheast, pnorthwest, psoutheast, psouthwest
USE Field_Types,  ONLY: fld_type_p, fld_type_v
USE model_domain_mod, ONLY: mt_global, mt_lam

USE fill_external_halos_mod, ONLY:  &
  fill_external_halos

USE halo_exchange_base_mod, ONLY:   &
  halo_exchange_base_init

USE halo_exchange_mpi_mod, ONLY:    &
  swap_bounds_mpi_rb,               &
  swap_bounds_mpi_aao,              &
  swap_bounds_mpi_nsew_wa  

USE halo_exchange_ddt_mod, ONLY:    &
  swap_bounds_ddt_aao,              &
  swap_bounds_ddt_nsew_wa  

IMPLICIT NONE

INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PRIVATE, PARAMETER :: zhook_out = 1

! Communicators to use for comms
INTEGER, PRIVATE :: global_comm
INTEGER, PRIVATE :: ew_comm
INTEGER, PRIVATE :: ns_comm
INTEGER, PRIVATE :: diag_comm 
INTEGER, PRIVATE :: diag_p_comm ! Diagonal polar, one for North and South poles.

INTEGER, PRIVATE, PARAMETER :: all_neighbours = 8

INTEGER          :: he_neighbour(all_neighbours)

LOGICAL          :: HE_Initialised=.FALSE.

! Tags for some of the communications
INTEGER, PARAMETER :: tag_to_e     = 501       
INTEGER, PARAMETER :: tag_to_w     = 502       
INTEGER, PARAMETER :: tag_from_e   = tag_to_w
INTEGER, PARAMETER :: tag_from_w   = tag_to_e
INTEGER, PARAMETER :: tag_to_s     = 503         ! not polar  
INTEGER, PARAMETER :: tag_to_n     = 504         ! not polar
INTEGER, PARAMETER :: tag_from_s   = tag_to_n    ! not polar
INTEGER, PARAMETER :: tag_from_n   = tag_to_s    ! not polar
INTEGER, PARAMETER :: tag_to_s_p   = 505         ! polar  
INTEGER, PARAMETER :: tag_to_n_p   = 506         ! polar
INTEGER, PARAMETER :: tag_from_s_p = tag_to_s_p  ! polar
INTEGER, PARAMETER :: tag_from_n_p = tag_to_n_p  ! polar
INTEGER, PARAMETER :: tag_to_ne    = 507         
INTEGER, PARAMETER :: tag_to_se    = 508         
INTEGER, PARAMETER :: tag_to_sw    = 509         
INTEGER, PARAMETER :: tag_to_nw    = 510         
INTEGER, PARAMETER :: tag_from_ne  = tag_to_sw   
INTEGER, PARAMETER :: tag_from_se  = tag_to_nw   
INTEGER, PARAMETER :: tag_from_sw  = tag_to_ne   
INTEGER, PARAMETER :: tag_from_nw  = tag_to_se   
INTEGER, PARAMETER :: tag_to_ne_p  = 511         ! polar
INTEGER, PARAMETER :: tag_to_se_p  = 512         ! polar
INTEGER, PARAMETER :: tag_to_sw_p  = 513         ! polar
INTEGER, PARAMETER :: tag_to_nw_p  = 514         ! polar
INTEGER, PARAMETER :: tag_from_ne_p= tag_to_nw_p ! polar
INTEGER, PARAMETER :: tag_from_nw_p= tag_to_ne_p ! polar  
INTEGER, PARAMETER :: tag_from_se_p= tag_to_sw_p ! polar
INTEGER, PARAMETER :: tag_from_sw_p= tag_to_se_p ! polar  

! Interface declaration for the swap_bound routines

INTERFACE swap_bounds_hub
MODULE PROCEDURE swap_bounds_hub_dp,    &
                 swap_bounds_hub_sp,    &
                 swap_bounds_hub_int,   &
                 swap_bounds_hub_log,   &
                 swap_bounds_hub_2D_dp, &
                 swap_bounds_hub_2D_sp, &
                 swap_bounds_hub_2D_int,&
                 swap_bounds_hub_2D_log
END INTERFACE

INTERFACE swap_bounds
MODULE PROCEDURE swap_bounds_hub_dp,    &
                 swap_bounds_hub_sp,    &
                 swap_bounds_hub_int,   &
                 swap_bounds_hub_log,   &
                 swap_bounds_hub_2D_dp, &
                 swap_bounds_hub_2D_sp, &
                 swap_bounds_hub_2D_int,&
                 swap_bounds_hub_2D_log
END INTERFACE

INTERFACE swap_bounds_RB
MODULE PROCEDURE swap_bounds_RB_hub_dp,    &
                 swap_bounds_RB_hub_sp
END INTERFACE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HALO_EXCHANGE'

CONTAINS

!==============================================================================
! Purpose:
! Initialise the module, must be called before other routines
! Set up communicators for halo exchange ops
!==============================================================================
SUBROUTINE Halo_Exchange_Init()
USE umPrintMgr, ONLY: umPrint, umMessage, PrNorm

IMPLICIT NONE

INTEGER             :: ierror
INTEGER             :: row_start
INTEGER             :: ns_comm_size
INTEGER             :: ns_index_delta
INTEGER             :: diag_p_comm_size
INTEGER             :: key
INTEGER             :: colour
INTEGER             :: my_proc_col
INTEGER             :: my_proc_row
INTEGER             :: proc_id_north
INTEGER             :: proc_id_south
INTEGER             :: proc_id_east
INTEGER             :: proc_id_west
INTEGER             :: proc_id_north_east
INTEGER             :: proc_id_south_east
INTEGER             :: proc_id_south_west
INTEGER             :: proc_id_north_west
INTEGER             :: half_nproc_x
INTEGER             :: half_nproc_x_is_odd
INTEGER             :: mycol_is_odd


LOGICAL             :: isOverpole

REAL(KIND=jprb)     :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HALO_EXCHANGE_INIT'



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up convenience indexing of processor neighbours
proc_id_north      =nodomain
proc_id_south      =nodomain
proc_id_east       =nodomain
proc_id_west       =nodomain
proc_id_north_east =nodomain
proc_id_south_east =nodomain
proc_id_south_west =nodomain
proc_id_north_west =nodomain

CALL gc_get_communicator(global_comm,ierror)
row_start      = (mype/nproc_x) * nproc_x
my_proc_col    = mype - row_start
my_proc_row    = (mype/nproc_x)

!------------------------------------------------------------------
! E/W communicator

! Define a row communicator in natural order, moving from west (0)
! to east.
colour         =my_proc_row
key            =my_proc_col
CALL MPL_Comm_Split(global_comm, colour, key, ew_comm, ierror)

! EW neighbours are trivial due to the natural ordering of
! the communicator
IF (neighbour(peast) /= nodomain) THEN
  proc_id_east=my_proc_col+1
  IF (proc_id_east==nproc_x) proc_id_east=0
END IF

IF (neighbour(pwest) /= nodomain) THEN
  proc_id_west=my_proc_col-1
  IF (proc_id_west==-1) proc_id_west=nproc_x-1
END IF

!------------------------------------------------------------------
! N/S communicator

! Define a column communicator in natural order, moving from south (0), in
! an initially northerly direction.

! For the global model a column is more complex. A column is a complete
! south->north->south circumnavigation. Columns start at the south pole
! in the west hemisphere traverse the north pole and end at the south
! pole in the east hemisphere, thus the set of 'columns' forms a set of
! nproc_x/2 rings.
!
! LAM and Torus2D global models form more conventional columns without/with
! wrapping boundary conditions respectively.

colour         =my_proc_col
key            =my_proc_row
isOverpole     =.FALSE.
IF (nproc_x>1) THEN  ! it is an even number by definition
  IF (sb_model_type == mt_global .AND. bound(2) /= bc_cyclic) THEN
    IF (my_proc_col>=nproc_x/2) THEN
      ! Easter hemisphere processors want to join the communicator
      ! that loops over the north pole from the
      ! western hemisphere
      colour=my_proc_col-nproc_x/2
      ! and I want to index such that the +1
      ! direction is south (the continued direction of north
      ! as you move over the pole)
      key   =2*nproc_y-my_proc_row-1
      isOverPole=.TRUE.
    END IF
  END IF
END IF

CALL MPL_Comm_Split(global_comm, colour, key, ns_comm, ierror)
CALL MPL_Comm_Size(ns_comm ,ns_comm_size  ,ierror)


! For NS the direction of northerly flips at the pole
! so we change our 'increment' to be negative.
ns_index_delta=1                 ! normally north is +1
IF (isOverPole)ns_index_delta=-1 ! except in the e. hemisphere
!                                ! in global models

! If there is only me, then my neighbour must be me, also.
IF (nproc_y==1 .AND. nproc_x==1)ns_index_delta=0

IF (neighbour(pnorth) /= nodomain) THEN
  proc_id_north=key+ns_index_delta

  ! If we are at maximal processor in the south->north
  ! direction, then we are either
  ! - a global model, and hence at the southerly pole
  !   in east hemisphere, so traversing the south pole
  !   takes us back to rank 0 (southerly in the western
  !   hemisphere)
  ! - a cyclic model on the northerly edge, and north
  !   wraps around to the south edge, and hence back
  !   to rank 0
  IF (proc_id_north==ns_comm_size) proc_id_north=0

END IF

IF (neighbour(psouth) /= nodomain) THEN
  proc_id_south=key-ns_index_delta

  ! adjust s. hem southerly row south to max in communicator
  ! wrapping under the pole for global etc etc
  IF (proc_id_south==-1)           proc_id_south=ns_comm_size-1
  IF (proc_id_south==ns_comm_size) proc_id_south=0
END IF


!------------------------------------------------------------------
! Corners communicator

! Define a SW to NE and SE to NW diagonal communicator. This is a copy
! of the global communicator. Attention all the following code
! regarding PEs uses indices starting from zero.
CALL MPL_Comm_Dup(global_comm, diag_comm, ierror)

! NE neighbour 
IF (neighbour(pnorth) /= nodomain .AND. neighbour(peast) /= nodomain) THEN
  proc_id_north_east=MODULO((my_proc_row+1)*nproc_x,nproc) +     & 
    MODULO(my_proc_col+1,nproc_x)
END IF

IF (neighbour(psouth) /= nodomain .AND. neighbour(peast) /= nodomain) THEN
  proc_id_south_east=MODULO((my_proc_row-1)*nproc_x,nproc) +     &
    MODULO(my_proc_col+1,nproc_x)
END IF

IF (neighbour(psouth) /= nodomain .AND. neighbour(pwest) /= nodomain) THEN
  proc_id_south_west=MODULO((my_proc_row-1)*nproc_x,nproc) +     &
    MODULO(my_proc_col-1,nproc_x)
END IF

IF (neighbour(pnorth) /= nodomain .AND. neighbour(pwest) /= nodomain) THEN
  proc_id_north_west=MODULO((my_proc_row+1)*nproc_x,nproc) +     &
    MODULO(my_proc_col-1,nproc_x)
END IF

!------------------------------------------------------------------
! Corners Over Poles communicator

! Define Overpoles diagonal communicator. Attention all the
! following code and comments regarding PEs uses indices starting
! from zero.

! The idea is to set up the communicator in order that the
! neighbours' values follow a natural W->E ordering, i.e. rank+1
! for east PE neighbour and rank-1 for west PE neighbour.

! All the PEs that are at the poles need to communicate with the
! PEs that are on the same pole, i.e. on the same row. Thus, the
! diagonal overpole communicator is a subgroup of the E/W
! communicator, since it groups PEs by row.

! All the PEs that are not at the extremeties, i.e. not at the poles,
! should never send messages using this communicator. Thus their
! colour is undentified and the communicator null.
IF (.NOT. at_extremity(pnorth) .AND. .NOT. at_extremity(psouth) ) THEN 
  colour           = mpl_undefined      
  key              = 0
  diag_p_comm_size = 1
ELSE
  ! The way that PEs are connected to communicate the corners overpole
  ! depends on the number of PEs in the X direction and the column
  ! number of each PE.

  ! Consider the number of PEs in the X direction (nproc_x). If
  ! this is halved (nproc_x/2) and the resulting number is even
  ! (MODULE(nproc_x/2,2)==0), all the PEs on a pole are on the
  ! same communicator, i.e. same colour.
  ! If nproc_x/2 is odd, then column even PEs and column odd PEs
  ! on a pole are on different communicators, i.e. one colour
  ! for the odd and one colour for the even.

  half_nproc_x         = nproc_x/2
  half_nproc_x_is_odd  = MODULO(half_nproc_x,2)
  mycol_is_odd         = MODULO(my_proc_col,2)

  IF (half_nproc_x_is_odd==0) THEN
    ! Ex with nproc_x = 8 and  half_nproc_x = 4
    ! mycol   0  1  2  3  4  5  6  7  
    ! colour  0  0  0  0  0  0  0  0  
    ! Key     0  5  2  7  4  1  6  3
    colour = 0        
    diag_p_comm_size = nproc_x
  ELSE
    ! Ex with nproc_x = 10 and  half_nproc_x = 5
    ! mycol   0  1  2  3  4  5  6  7  8  9
    ! colour  0  1  0  1  0  1  0  1  0  1
    ! Key     0  0  2  2  4  4  1  1  3  3
    colour = mycol_is_odd
    diag_p_comm_size = half_nproc_x
  END IF

  ! About the key (rank of the PE in the same colour)....

  ! If the column is even, then the key value is the remainder of
  ! the column value divided by the size of the communicator, i.e
  ! MODULO(mycol,diag_p_comm_size)

  ! If the column is odd and the half nproc_x/2 is odd, then the
  ! key value is the remainder of the column value minus 1 divided
  ! by the size of the communicator, i.e
  ! MODULO(mycol-1,diag_p_comm_size)

  ! If the column is odd and the half nproc_x/2 is even, then the
  ! key value is the remainder of the column value plus nproc_x/2
  ! divided by the size of the communicator, i.e
  ! MODULO(mycol+half_nproc_x,diag_p_comm_size)

  key = MODULO((my_proc_col-colour)+(half_nproc_x*mycol_is_odd), & 
    diag_p_comm_size)

END IF

! Notice the ew_comm is the parent communicator
CALL MPL_Comm_Split(ew_comm, colour, key, diag_p_comm, ierror)


IF (bound(2) == bc_overpole) THEN 
  ! The polar SW NE SE and NW neighbours values follow the natural
  ! W -> E ordering of the diagonal polar communicator.
  IF (at_extremity(pnorth)) THEN
    IF (neighbour(pnorth) /= nodomain .AND. neighbour(peast) /= nodomain) THEN
      proc_id_north_east=key+1
      IF (proc_id_north_east==diag_p_comm_size) proc_id_north_east=0
    END IF
    IF (neighbour(pnorth) /= nodomain .AND. neighbour(pwest) /= nodomain) THEN
      proc_id_north_west=key-1
      IF (proc_id_north_west==-1) proc_id_north_west=diag_p_comm_size-1
    END IF
  END IF
  IF (at_extremity(psouth)) THEN
    IF (neighbour(psouth) /= nodomain .AND. neighbour(peast) /= nodomain) THEN
      proc_id_south_east=key+1
      IF (proc_id_south_east==diag_p_comm_size) proc_id_south_east=0
    END IF
    IF (neighbour(psouth) /= nodomain .AND. neighbour(pwest) /= nodomain) THEN
      proc_id_south_west=key-1
      IF (proc_id_south_west==-1) proc_id_south_west=diag_p_comm_size-1
    END IF
  END IF
END IF

!------------------------------------------------------------------
! Store indexing

he_neighbour(pEast)=proc_id_east
he_neighbour(pWest)=proc_id_west
he_neighbour(pNorth)=proc_id_north
he_neighbour(pSouth)=proc_id_south
he_neighbour(pNorthEast)=proc_id_north_east
he_neighbour(pSouthEast)=proc_id_south_east
he_neighbour(pSouthWest)=proc_id_south_west
he_neighbour(pNorthWest)=proc_id_north_west


HE_Initialised=.TRUE.

! The method for 2D fields is the same one for 3D fields with small halos.
! This is to reduce the number of input variables.
swap_bounds_method_2D = swap_bounds_method_3D_small

CALL halo_exchange_base_init()

WRITE(umMessage,'(A,A40)') 'swap_bounds_method_3D_large  = ', &
  name_sb_method(swap_bounds_method_3D_large)
CALL umPrint(umMessage,src='halo_exchange',level=PrNorm)

WRITE(umMessage,'(A,A40)') 'swap_bounds_method_3D_small  = ', &
  name_sb_method(swap_bounds_method_3D_small)
CALL umPrint(umMessage,src='halo_exchange',level=PrNorm)

WRITE(umMessage,'(A,A40)') 'swap_bounds_method_2D        = ', &
  name_sb_method(swap_bounds_method_2D)
CALL umPrint(umMessage,src='halo_exchange',level=PrNorm)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Halo_Exchange_Init

!****************************!
!     SWAP BOUNDS HUB        !
!****************************!

! Call the specific swap bounds routine depending on the size of the
! field and halos.

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_dp(          &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_DP'

#include "halo_exchange_mod_swap_bounds_hub.h"

RETURN
END SUBROUTINE swap_bounds_hub_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_sp(          &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_SP'

#include "halo_exchange_mod_swap_bounds_hub.h"

RETURN
END SUBROUTINE swap_bounds_hub_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_int(         &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_INT'

#define   HALO_EXCHANGE_INTEGER 1
#include "halo_exchange_mod_swap_bounds_hub.h"
#undef    HALO_EXCHANGE_INTEGER

RETURN
END SUBROUTINE swap_bounds_hub_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_log(         &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_LOG'

#define   HALO_EXCHANGE_LOGICAL 1
#include "halo_exchange_mod_swap_bounds_hub.h"
#undef    HALO_EXCHANGE_LOGICAL

RETURN
END SUBROUTINE swap_bounds_hub_log

! ---------------------------------------------------------------------------- !

!****************************!
!     SWAP BOUNDS HUB 2D     !
!****************************!

! Forward the routine call to the 3D version.

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_2D_dp(       &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_2D_DP'

#include "halo_exchange_mod_swap_bounds_hub_2D.h"

RETURN
END SUBROUTINE swap_bounds_hub_2D_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_2D_sp(       &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_2D_SP'

#define   HALO_EXCHANGE_SINGLE 1
#include "halo_exchange_mod_swap_bounds_hub_2D.h"
#undef    HALO_EXCHANGE_SINGLE

RETURN
END SUBROUTINE swap_bounds_hub_2D_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_2D_int(      &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_2D_INT'

#define   HALO_EXCHANGE_INTEGER 1
#include "halo_exchange_mod_swap_bounds_hub_2D.h"
#undef    HALO_EXCHANGE_INTEGER

RETURN
END SUBROUTINE swap_bounds_hub_2D_int

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_hub_2D_log(      &
    field, row_length, rows, levels,    & 
    halo_x, halo_y,                     & 
    field_type, l_vector,               & 
    do_east_arg, do_west_arg,           & 
    do_south_arg, do_north_arg,         & 
    do_corners_arg                      & 
    )

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_HUB_2D_LOG'

#define   HALO_EXCHANGE_LOGICAL 1
#include "halo_exchange_mod_swap_bounds_hub_2D.h"
#undef    HALO_EXCHANGE_LOGICAL

RETURN
END SUBROUTINE swap_bounds_hub_2D_log

! ---------------------------------------------------------------------------- !

!****************************!
!     SWAP BOUNDS RB HUB     !
!****************************!

!   This is a special subroutine that performs boundary swapping in EW
!   and NS directions when all the following are true:
!   1. The data is computed in a red or black order (like a chessboard)
!   2. The halos have fixed size of halo_x=halo_y=1.
!   3. The corners of the halos are NOT swapped
!   4. The field type is fld_type_p
!   5. The data does not have a horizontal vector component
!   This subroutine has the following limitations:
!   1. Does NOT work if nproc_x or nproc_y  == 1
!   2. Does NOT work if the model type is NOT (mt_global or mt_lam)
!   3. Does NOT work if the model type is mt_global and
!      global_row_length/2 is not even
!   In these cases the routine calls the generic swap_bounds
!
!   If red is .TRUE. then the red elements are swapped. The element in position
!   (1,1) in the global grid is red.
!    

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_RB_hub_dp(       &
    field, row_length, rows, levels,    & 
    red                                 & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_RB_HUB_DP'

#include "halo_exchange_mod_swap_bounds_RB_hub.h"

RETURN
END SUBROUTINE swap_bounds_RB_hub_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE swap_bounds_RB_hub_sp(       &
    field, row_length, rows, levels,    & 
    red                                 & 
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_RB_HUB_SP'

#include "halo_exchange_mod_swap_bounds_RB_hub.h"

RETURN
END SUBROUTINE swap_bounds_RB_hub_sp

! ---------------------------------------------------------------------------- !

END MODULE Halo_Exchange
