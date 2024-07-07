! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

!
! This module provides methods for updating boundaries (halos) on
! real fields using non-blocking swap-bounds, i.e. the routine that start
! the swap bounds returns without waiting for the communication to
! complete so that computation can continue while the halo exchange is
! happening on the background.
!
! This module is an extension of the halo_exchange module and the
! halo_exchange_init() method must be called before other methods of
! this module can be used.
!
MODULE Non_Blocking_Halo_Exchange

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_LOC, C_INTPTR_T

USE um_types, ONLY:  &
     real32,         &
     real64

USE mpl, ONLY:       &
  mpl_real4,                        &
  mpl_real8,                        &
  mpl_request_null 

USE tags_params, ONLY:              &
  nb_tag_base,                      &
  nb_max_tag_index

USE mpp_conf_mod, ONLY:             &
  name_sb_method,                   &
  ddt_aao,                          &
  ddt_nsew_wa,                      &
  nb_swap_bounds_off,               &
  nb_swap_bounds_method          

USE fill_external_halos_mod, ONLY:  &
  fill_external_halos

USE tools_halo_exchange_mod, ONLY:   &
  normal_halo,                       &
  full_halo 

USE halo_exchange, ONLY:            & 
  he_initialised,                   &
  swap_bounds

USE halo_exchange_ddt_mod, ONLY:    &
  ddt_type,                         &
  free_ddt,                         &
  create_ddt_ew,                    &
  create_ddt_ns,                    &
  create_ddt_cr,                    &
  begin_swap_bounds_ddt_aao,        &
  begin_swap_bounds_ddt_nsew_wa,    &  
  end_swap_bounds_ddt_aao,          &
  end_swap_bounds_ddt_nsew_wa

USE missing_data_mod, ONLY: imdi

USE ereport_mod,   ONLY: ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParVars, ONLY: halo_j, halo_i,sb_model_type
USE model_domain_mod, ONLY: mt_lam

USE fort2c_addr_mod, ONLY: um_addr_value

IMPLICIT NONE

PRIVATE 

! ----------------------------------------------------------------------
! ATTENTION: 
!
! Multiple non-blocking swap bounds can be active at the same time
! on different fields of different type and size.
!
! Furthermore, a real field can have multiple non-blocking swap bounds
! active at different levels of the field.
!
! Linked list are used to manage all these independent active swap
! bounds executions which are represented as nodes in the list. Each
! non-blocking communication is uniquely identified in the list by a
! KEY. This key is the memory address of the field passed to the
! non-blocking routines as returned by the internal get_key function
! which uses the C_LOC intrinsic and um_addr_value C routine. Thus if
! the same field is passed to the non-blocking routines at different
! levels, different independent simultaneously active swap-bounds can
! be executed.
! 
! It is very important to always pass the field at the
! same address/index into the routines of this module. 
!
! WARNING: 
!
! In the Fortran standard 2003 a compiler is free to
! implement the C_PTR type any way it feels as long as it a
! derived type that contains a compatible pointer address
! component, which is what gets passed to any BIND(c)
! procedure. This mean that in some compilers what is returned by
! C_LOC(field) call may not actually be the address of the
! field. Thus in order to be portable the get_key routine passes the
! C_LOC(field) value to the um_addr_value C routine which returns the
! integer address of the field.
!
! LIMITATIONS:
!
! The begin and end swap bounds routines are not thread-safe. This
! means that they need to be called outiside an OpenMP parallel
! region. If they are called inside a parallel region, they need to be
! protected by a MASTER or CRITICAL region.
!
! ----------------------------------------------------------------------

! The type used to store a key. Since a key is a memory address which
! is retrieved using C_LOC and um_addr_value C function, the type is
! the largest integer that can store a C pointer.
INTEGER, PARAMETER :: keytype = C_INTPTR_T

! The status of the non-blocking swap bounds of a field.
INTEGER, PARAMETER :: nb_ok            = 0 ! The communication is finished
                                           ! or/and the field is ready 
INTEGER, PARAMETER :: nb_started       = 1 ! Performing a swap bound


! The maximum number of MPL requests needed by a non-blocking
! swap bounds for each action, (8 four send and 8 for recv).
INTEGER, PARAMETER :: nb_requests_max = 8


! Contains the information about a field associated with a
! non-blocking swap-bounds. This type acts as a node into the linked
! list.
TYPE nb_field_node_type

  INTEGER(KIND=keytype) :: key          ! Identified the field in the
                                        ! list. This is the address
                                        ! returned by the LOC
                                        ! intrinsic.

  INTEGER :: tag_index                  ! Unique index that is used
                                        ! for the unique tags in the
                                        ! swap bounds

  INTEGER :: status                     ! The status of the communication.

  INTEGER :: rows                       ! number of rows in a field
                                        ! (not including halos)
  INTEGER :: row_length                 ! number of points on a row
                                        ! (not including halos)
  INTEGER :: levels                     ! number of model levels
  INTEGER :: halo_x                     ! size of halo in "i" direction
  INTEGER :: halo_y                     ! size of halo in "j" direction

  INTEGER :: field_type                 ! Defines the grid interpolation type
                                        ! of the input FIELD (u,v or w)
  LOGICAL :: l_vector                   ! TRUE:  Data is a horizontal vector
                                        !     component
                                        ! FALSE: Data is a scalar

  LOGICAL :: do_east                    ! Perform east swap
  LOGICAL :: do_west                    ! Perform west swap
  LOGICAL :: do_south                   ! Perform south swap
  LOGICAL :: do_north                   ! Perform north swap
  LOGICAL :: do_corners                 ! Perform corners swap

  ! Variables used by specific implementations

  ! Common to DDT 
  INTEGER :: nreq_recv                  ! The number of MPL receive requests
  INTEGER :: ireq_recv(nb_requests_max) ! The MPL receive requests

  ! Common to NSEW_WA
  INTEGER  :: ns_halos_recv             ! Number of N/S halos received

  ! Pointers for the linked list
  TYPE(nb_field_node_type), POINTER :: prev
  TYPE(nb_field_node_type), POINTER :: next 
    
END TYPE nb_field_node_type

! The type used to store and manage the list
TYPE nb_field_list_type

  ! The global tag_index for this list
  INTEGER :: tag_index

  ! Pointers to the first elements in the list
  TYPE(nb_field_node_type), POINTER :: head 
END TYPE nb_field_list_type

! Master field list object.
TYPE(nb_field_list_type), SAVE :: nb_field_list


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER :: ModuleName='NON_BLOCKING_HALO_EXCHANGE'

! Private routines used to retrieve the key, i.e. the memory address.
INTERFACE  get_key
MODULE PROCEDURE get_key_3d_dp, get_key_3d_sp, get_key_2d_dp
END INTERFACE

! Public routines

! Initlaise the data needed for non blocking swap bounds
PUBLIC :: non_blocking_halo_exchange_init

! Start the non_blocking swap bounds of the given field.
PUBLIC ::  begin_swap_bounds
INTERFACE  begin_swap_bounds
MODULE PROCEDURE begin_swap_bounds_hub_dp,     &
                 begin_swap_bounds_hub_sp,     &
                 begin_swap_bounds_hub_2D_dp  
END INTERFACE

! End the non_blocking swap bounds of the given field.
PUBLIC ::  end_swap_bounds
INTERFACE  end_swap_bounds
MODULE PROCEDURE end_swap_bounds_hub_dp,       &
                 end_swap_bounds_hub_sp,       &
                 end_swap_bounds_hub_2D_dp
END INTERFACE

CONTAINS

!==============================================================================
! Purpose:
! Return the the key, i.e the memory address, of a double precision 3D array 
!==============================================================================
FUNCTION get_key_3d_dp ( field ) RESULT (key)

IMPLICIT NONE

REAL(KIND=real64), TARGET, INTENT(IN) :: field(:,:,:) 
INTEGER(KIND=keytype) :: key

key = um_addr_value(C_LOC(field))

END FUNCTION get_key_3d_dp

!==============================================================================
! Purpose:
! Return the the key, i.e the memory address, of a single precision 3D array 
!==============================================================================
FUNCTION get_key_3d_sp ( field ) RESULT (key)

IMPLICIT NONE

REAL(KIND=real32), TARGET, INTENT(IN) :: field(:,:,:) 
INTEGER(KIND=keytype) :: key

key = um_addr_value(C_LOC(field))

END FUNCTION get_key_3d_sp

!==============================================================================
! Purpose:
! Return the the key, i.e the memory address, of a double precision 2D array 
!==============================================================================
FUNCTION get_key_2d_dp ( field ) RESULT (key)

IMPLICIT NONE

REAL(KIND=real64), TARGET, INTENT(IN) :: field(:,:) 
INTEGER(KIND=keytype) :: key

key = um_addr_value(C_LOC(field))

END FUNCTION get_key_2d_dp

!==============================================================================
! Purpose:
! Initialise the linked list.
! Attention, this routine is not thread-safe .
!==============================================================================
SUBROUTINE init_nb_field_list ( list )

IMPLICIT NONE

TYPE(nb_field_list_type), INTENT(INOUT) :: list

! Make sure the list pointer is not associated with any nodes.
NULLIFY(list % head)

! Initialise the global tag_index
list % tag_index = 1

END SUBROUTINE init_nb_field_list

!==============================================================================
! Purpose:
! Return the existing or new field node from the list with the given key
! Attention, this routine is not thread-safe .
!==============================================================================
FUNCTION get_nb_field_node ( list, key) RESULT (node)

IMPLICIT NONE

TYPE(nb_field_list_type), INTENT(INOUT) :: list
INTEGER(KIND=keytype),    INTENT(IN)    :: key

TYPE(nb_field_node_type), POINTER       :: node

! Set the node to be the first element in the list, if the list is
! empty the node will not be associated
node => list % head

! Loop through the nodes in the list until the node is not a null
! pointer or a matching key was found. 
DO WHILE ( ASSOCIATED(node) )
  IF (node % key == key) THEN 
    RETURN 
  END IF
  node => node % next  
END DO

! If the routine reach this point, then no matching key was found and
! thus we need to create a new node. 
ALLOCATE(node)

! Set the key,  status and all the pointers to not associated.
node % key    = key
node % status = nb_ok
NULLIFY(node % next)
NULLIFY(node % prev)

! Set the tag_index from the global one and increase it. If the global
! index values is too high, reset it.
node % tag_index = list % tag_index
IF(list % tag_index < nb_max_tag_index) THEN
  list % tag_index = list % tag_index +1
ELSE
  list % tag_index = 1
END IF

! Make sure that the eventual existing head become the next of the new
! node and its previous is linked to the new node.
IF(ASSOCIATED( list % head )) THEN
  node % next         => list % head
  list % head % prev  => node
END IF

! Set the new node as the head of the list
list % head   => node

RETURN
END FUNCTION get_nb_field_node

!==============================================================================
! Purpose:
! Remove the given node from the list
! Attention, this routine is not thread-safe .
!==============================================================================
SUBROUTINE remove_nb_field_node ( list, node)

IMPLICIT NONE

TYPE(nb_field_list_type), INTENT(INOUT) :: list
TYPE(nb_field_node_type), POINTER       :: node

! If the node is not associated, then do nothing.
IF( .NOT. ASSOCIATED(node)) THEN
  RETURN
END IF

! If the previous node is associated, then set its next to be this
! node next. Do the same for the next node.
IF( ASSOCIATED(node % prev) ) THEN 
  NULLIFY(node % prev % next)
  IF (ASSOCIATED(node % next) ) THEN
    node % prev % next => node % next
  END IF
END IF
IF( ASSOCIATED(node % next) ) THEN 
  NULLIFY(node % next % prev)
  IF (ASSOCIATED(node % prev) ) THEN
    node % next % prev => node % prev
  END IF
END IF

! If the head of the list is the node to remove, then set the next
! node to be the new head.
IF (ASSOCIATED(list % head, target=node)) THEN
  NULLIFY(list % head)
  IF (ASSOCIATED(node % next) ) THEN
    list % head => node % next
  END IF
END IF

! Deallocate node
DEALLOCATE(node)

RETURN
END SUBROUTINE remove_nb_field_node


!==============================================================================
! Purpose:
! Initialise the module, must be called before other routines
!==============================================================================
SUBROUTINE Non_Blocking_Halo_Exchange_Init()

USE umPrintMgr, ONLY: & 
  umPrint,            &
  umMessage,          &
  PrNorm

USE nlsizes_namelist_mod, ONLY: &
  row_length,                   &
  rows,                         &
  model_levels

IMPLICIT NONE

! Local variables
INTEGER :: f
INTEGER :: errorstatus

! Used for computing message size statistics
INTEGER :: size_small
INTEGER :: size_large


! Read env
CHARACTER(LEN=10) :: env_msg  ! To get environmental variable

REAL(KIND=jprb)   :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NON_BLOCKING_HALO_EXCHANGE_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. he_initialised) THEN
  errorstatus = 99
  CALL ereport(RoutineName,errorstatus,                                    &
               ' halo exchange not initialised')
END IF

! Initialise field list
CALL init_nb_field_list(nb_field_list)

! Check that the non-blocking swap bounds method selected is one of
! the working ones.
IF(nb_swap_bounds_method /= nb_swap_bounds_off .AND.                       &
   nb_swap_bounds_method /= ddt_aao .AND.                                  &
   nb_swap_bounds_method /= ddt_nsew_wa)  THEN

  errorstatus = -99
  CALL ereport(RoutineName,errorstatus,                                    &
               ' nb_swap_bounds_method is unknown, set to Off (blocking)')  

  ! If an unknown swap bounds method, then use a blocking method.
  nb_swap_bounds_method = nb_swap_bounds_off
END IF
    

! Print information

! Type of non-blocking swap-bounds
WRITE(umMessage,'(A,A40)') 'nb_swap_bounds_method        = ', &
  name_sb_method(nb_swap_bounds_method)
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

WRITE(umMessage,"(A52)")  "----------------------------------------------------"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)
WRITE(umMessage,"(A52)") "|        Typical non blocking message sizes        |"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)
WRITE(umMessage,"(A52)")  "| Dir Type  Small Halo Size(B)  Large Halo Size(B) |"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

! Compute E/W size (no corners)
size_small = rows * model_levels
size_large = rows * halo_i * model_levels

WRITE(umMessage,"(A10,I19,I20,A3)") "| E/W  SP ",size_small*4,size_large*4,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)
WRITE(umMessage,"(A10,I19,I20,A3)") "| E/W  DP ",size_small*8,size_large*8,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

! Compute N/S size (no corners)
size_small = row_length * model_levels
size_large = row_length * halo_j * model_levels

WRITE(umMessage,"(A10,I19,I20,A3)") "| N/S  SP ",size_small*4,size_large*4,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)
WRITE(umMessage,"(A10,I19,I20,A3)") "| N/S  DP ",size_small*8,size_large*8,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

! Compute Corners size 
size_small = model_levels
size_large = halo_i * halo_j * model_levels

WRITE(umMessage,"(A10,I19,I20,A3)") "| Crns SP ",size_small*4,size_large*4,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)
WRITE(umMessage,"(A10,I19,I20,A3)") "| Crns DP ",size_small*8,size_large*8,"|"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

WRITE(umMessage,"(A52)")  "----------------------------------------------------"
CALL umPrint(umMessage,src='non_blocking_halo_exchange',level=PrNorm)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE Non_Blocking_Halo_Exchange_Init


!****************************!
!         BEGIN HUB          !
!****************************!

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_hub_dp(    &
    field, row_length, rows, levels,    & ! field
    halo_x, halo_y,                     & ! halos
    field_type, l_vector,               & ! supporting information
    do_east_arg, do_west_arg,           & ! Optionals
    do_south_arg, do_north_arg,         & ! Optionals
    do_corners_arg,                     & ! Optional
    field_node_arg                      & ! Optional
    )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8
INTEGER, PARAMETER :: byte_width = 8

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_HUB_DP'

#include "non_blocking_halo_exchange_mod_begin_swap_bounds_hub.h"

RETURN
END SUBROUTINE begin_swap_bounds_hub_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_hub_sp(    &
    field, row_length, rows, levels,    & ! field
    halo_x, halo_y,                     & ! halos
    field_type, l_vector,               & ! supporting information
    do_east_arg, do_west_arg,           & ! Optionals
    do_south_arg, do_north_arg,         & ! Optionals
    do_corners_arg,                     & ! Optional
    field_node_arg                      & ! Optional
    )    

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4
INTEGER, PARAMETER :: byte_width = 4

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_HUB_SP'

#include "non_blocking_halo_exchange_mod_begin_swap_bounds_hub.h"

RETURN
END SUBROUTINE begin_swap_bounds_hub_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE begin_swap_bounds_hub_2D_dp( &
    field, row_length, rows, levels,    & ! field
    halo_x, halo_y,                     & ! halos
    field_type, l_vector,               & ! supporting information
    do_east_arg, do_west_arg,           & ! Optionals
    do_south_arg, do_north_arg,         & ! Optionals
    do_corners_arg                      & ! Optional
    )

IMPLICIT NONE

INTEGER, INTENT(IN) ::  row_length ! number of points on a row
!                                  ! (not including halos)
INTEGER, INTENT(IN) ::  rows       ! number of rows in a theta field
!                                  ! (not including halos)
INTEGER, INTENT(IN) ::  levels     ! number of model levels
INTEGER, INTENT(IN) ::  halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) ::  halo_y     ! size of halo in "j" direction
INTEGER, INTENT(IN) ::  field_type ! Defines the grid interpolation type
!                                  ! of the input FIELD (u,v or w)
LOGICAL, INTENT(IN) :: l_vector    ! TRUE:  Data is a horizontal vector
!                                  !     component
!                                  !     FALSE: Data is a scalar

! The next set of logicals control whether swaps are performed for a
! given direction. If none are present, all are done. If any are
! present, only the ones given TRUE values will be done.
LOGICAL, INTENT(IN), OPTIONAL :: do_east_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_west_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_south_arg
LOGICAL, INTENT(IN), OPTIONAL :: do_north_arg

! The next logical controls if corners are swapped or not. 
! If not, and single point halos are used, an optimised version can
! be used.
LOGICAL, INTENT(IN), OPTIONAL :: do_corners_arg

! Field to have its halos updated
REAL(KIND=real64), INTENT(INOUT) :: field(1-halo_x:row_length+halo_x, &
                                          1-halo_y:rows+halo_y) 
! Local variables

! The field node with all the information
TYPE(nb_field_node_type), POINTER       :: field_node

INTEGER(KIND=keytype) :: key

! Local variables needed to fixup for optional args.
LOGICAL :: do_east    ! perform east swap
LOGICAL :: do_west    ! perform west swap
LOGICAL :: do_south   ! perform south swap
LOGICAL :: do_north   ! perform north swap
LOGICAL :: do_corners ! do corner swaps

INTEGER :: error_code

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BEGIN_SWAP_BOUNDS_HUB_2D_DP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Setup the logicals

do_east    = .TRUE.
do_west    = .TRUE.
do_south   = .TRUE.
do_north   = .TRUE.
do_corners = .TRUE.

IF (PRESENT(do_corners_arg)) do_corners=do_corners_arg


IF (PRESENT(do_east_arg)  .OR. PRESENT(do_west_arg) .OR.  &
    PRESENT(do_south_arg) .OR. PRESENT(do_north_arg) ) THEN
  do_east  = .FALSE.
  do_west  = .FALSE.
  do_south = .FALSE.
  do_north = .FALSE.

  IF (PRESENT(do_east_arg))  do_east  = do_east_arg
  IF (PRESENT(do_west_arg))  do_west  = do_west_arg
  IF (PRESENT(do_south_arg)) do_south = do_south_arg
  IF (PRESENT(do_north_arg)) do_north = do_north_arg
END IF

!------------------------------------------------------------------
! Do blocking swap-bounds? 
! Call the default swap bounds routine and exit

IF (nb_swap_bounds_method == nb_swap_bounds_off) THEN
  
  CALL swap_bounds( field, row_length, rows, levels,    & 
                    halo_x, halo_y,                     &
                    field_type, l_vector,               &
                    do_east, do_west,                   & 
                    do_south, do_north,                 & 
                    do_corners                          & 
                    )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)  
  RETURN
END IF

!------------------------------------------------------------------
! Retrieve node

! Retrieve the node of the given key, i.e the addess location.
key = get_key( field )
field_node => get_nb_field_node(nb_field_list, key )

! Check that the swap bounds was started
IF(field_node % status /= nb_ok) THEN
  error_code = 99
  CALL ereport(RoutineName,error_code,                                     &
               ' wrong staus of non-blocking communication')  
END IF

!------------------------------------------------------------------
! Call the 3D field routine and pass the field node already found.

CALL begin_swap_bounds_hub_dp(field, row_length, rows, levels, &
                              halo_x, halo_y,                  &
                              field_type,l_vector,             &
                              do_east, do_west,                &
                              do_south,do_north,               &
                              do_corners, field_node)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE begin_swap_bounds_hub_2D_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_hub_dp ( field, field_node_arg )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64
INTEGER, PARAMETER :: mpl_type   = mpl_real8
INTEGER, PARAMETER :: byte_width = 8

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_HUB_DP'

#include "non_blocking_halo_exchange_mod_end_swap_bounds_hub.h"

RETURN
END SUBROUTINE end_swap_bounds_hub_dp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_hub_sp ( field, field_node_arg )

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32
INTEGER, PARAMETER :: mpl_type   = mpl_real4
INTEGER, PARAMETER :: byte_width = 4

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_HUB_SP'

#include "non_blocking_halo_exchange_mod_end_swap_bounds_hub.h"

RETURN 
END SUBROUTINE end_swap_bounds_hub_sp

! ---------------------------------------------------------------------------- !

SUBROUTINE end_swap_bounds_hub_2D_dp ( field )

IMPLICIT NONE

! Assumed shape field to have its halos updated
REAL(KIND=real64), INTENT(INOUT) :: field(:,:) 

! Local variables

INTEGER :: row_length ! Number of points on a row 
!                                 ! (not including halos)
INTEGER :: rows       ! Number of rows in a theta field
!                     ! (not including halos)
INTEGER :: levels     ! Number of model levels
INTEGER :: halo_x     ! Size of halo in "i" direction
INTEGER :: halo_y     ! Size of halo in "j" direction


! The field node with all the information
TYPE(nb_field_node_type), POINTER       :: field_node

INTEGER(KIND=keytype) :: key

INTEGER :: error_code

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_HUB_2D_DP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Do blocking swap-bounds? 
! This routine does nothing.
IF (nb_swap_bounds_method == nb_swap_bounds_off) THEN  
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)  
  RETURN
END IF

!------------------------------------------------------------------
! Retrieve node

! Retrieve the node of the given key, i.e the addess location.
key = get_key( field )
field_node => get_nb_field_node(nb_field_list, key )

! Check that the swap bounds was started
IF(field_node % status /= nb_started) THEN
  error_code = 99
  CALL ereport(RoutineName,error_code,                                     &
               ' wrong staus of non-blocking communication')  
END IF

!------------------------------------------------------------------
! Retrieve the field information

rows       = field_node % rows       
row_length = field_node % row_length 
levels     = field_node % levels     
halo_x     = field_node % halo_x     
halo_y     = field_node % halo_y     

!------------------------------------------------------------------
! Call the wrap routine and pass the field node already found.

CALL end_swap_bounds_wrap_dp(        & 
  field, row_length, rows, levels,   & 
  halo_x, halo_y, field_node)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE end_swap_bounds_hub_2D_dp


!==============================================================================
! Purpose:
! This is used by the 2D end swap bounds routine in order to call the
! end_swap_bounds_hub routine which uses the 3D array.
! ==============================================================================
SUBROUTINE end_swap_bounds_wrap_dp (    &
    field, row_length, rows, levels,    & ! field
    halo_x, halo_y, field_node)                     
  
IMPLICIT NONE

INTEGER, INTENT(IN) ::  row_length ! number of points on a row
!                                  ! (not including halos)
INTEGER, INTENT(IN) ::  rows       ! number of rows in a theta field
!                                  ! (not including halos)
INTEGER, INTENT(IN) ::  levels     ! number of model levels
INTEGER, INTENT(IN) ::  halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) ::  halo_y     ! size of halo in "j" direction

REAL(KIND=real64), INTENT(INOUT) :: field(1-halo_x:row_length+halo_x, &
                                          1-halo_y:rows+halo_y,levels) 

! The field node with all the information
TYPE(nb_field_node_type), INTENT(IN), POINTER  :: field_node

! Local variables

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='END_SWAP_BOUNDS_WRAP_DP'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL end_swap_bounds_hub_dp(field, field_node)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE end_swap_bounds_wrap_dp


END MODULE Non_Blocking_Halo_Exchange

