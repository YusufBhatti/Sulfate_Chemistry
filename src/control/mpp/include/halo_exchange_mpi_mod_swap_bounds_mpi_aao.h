! -- begin halo_exchange_mpi_mod_swap_bounds_mpi_aao.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP


INTEGER, INTENT(IN) :: row_length ! Number of points on a row 
!                                 ! (not including halos)
INTEGER, INTENT(IN) :: rows       ! Number of rows in a theta field
!                                 ! (not including halos)
INTEGER, INTENT(IN) :: levels     ! Number of model levels
INTEGER, INTENT(IN) :: halo_x     ! Size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! Size of halo in "j" direction
INTEGER, INTENT(IN) :: field_type ! Defines the grid interpolation type
!                                 !     of the input FIELD (u,v or w)
LOGICAL, INTENT(IN) :: l_vector   ! TRUE:  Horizontal vector
!                                 ! FALSE: Scalar

#if   defined(HALO_EXCHANGE_MPI_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_MPI_LOGICAL)
LOGICAL, INTENT(INOUT) ::                                                  &
#else
REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
#endif
  field(1-halo_x:row_length+halo_x,                                        &
        1-halo_y:rows+halo_y,                                              &
        levels)                   
  
! The next set of logicals control whether swaps are performed for a
! given direction. 
LOGICAL, INTENT(IN) :: do_east_arg
LOGICAL, INTENT(IN) :: do_west_arg
LOGICAL, INTENT(IN) :: do_south_arg
LOGICAL, INTENT(IN) :: do_north_arg
LOGICAL, INTENT(IN) :: do_corners_arg

! Local variables

! The buffers used to receive the data
#if   defined(HALO_EXCHANGE_MPI_INTEGER)
INTEGER ::                                                                 &
#elif defined(HALO_EXCHANGE_MPI_LOGICAL)
LOGICAL ::                                                                 &
#else
REAL(KIND=field_kind) ::                                                   &
#endif
  recv_buffer_e(halo_x, rows, levels),                                     &
  recv_buffer_w(halo_x, rows, levels),                                     &
  recv_buffer_n(row_length, halo_y, levels),                               &
  recv_buffer_s(row_length, halo_y, levels),                               &
  recv_buffer_ne(halo_x, halo_y, levels),                                  &
  recv_buffer_se(halo_x, halo_y, levels),                                  &
  recv_buffer_nw(halo_x, halo_y, levels),                                  &
  recv_buffer_sw(halo_x, halo_y, levels)

INTEGER  :: ireq_recv(8)
INTEGER  :: nreq_recv
INTEGER  :: tag_base

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use always the same tag
tag_base = 0

CALL begin_swap_bounds_mpi_aao(                   &
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


CALL end_swap_bounds_mpi_aao(                     &
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


! -- end halo_exchange_mpi_mod_swap_bounds_mpi_aao.h --
