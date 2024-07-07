! -- begin halo_exchange_mpi_mod_swap_bounds_mpi_rb.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

INTEGER, INTENT(IN) :: row_length ! number of points on a row
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! number of rows in a theta field
!                                 !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! number of model levels

LOGICAL, INTENT(IN) :: red        ! If .TRUE. the red elements need to
                                  ! be swapped

REAL(KIND=field_kind), INTENT(INOUT) ::                                    &
    field(0:row_length+1,0:rows+1,levels)  ! Field to have its halos updated

! Local variables

! The send and receive buffers for EW halos
REAL(KIND=field_kind) ::                                                   &
  send_buffer_e(rows/2+1, levels),                                         &
  send_buffer_w(rows/2+1, levels),                                         &
  recv_buffer_e(rows/2+1, levels),                                         &
  recv_buffer_w(rows/2+1, levels)

! The send and receive buffers for NS halos
REAL(KIND=field_kind) ::                                                   &
  send_buffer_n(row_length/2+1, levels),                                   &
  send_buffer_s(row_length/2+1, levels),                                   &
  recv_buffer_n(row_length/2+1, levels),                                   &
  recv_buffer_s(row_length/2+1, levels)


INTEGER  ::                &
  ns_halo_size,            & ! size of NS halo exchanged
  ew_halo_size               ! size of EW halo exchanged    

INTEGER  :: ierror                         ! MPL error code
INTEGER  :: istat(mpl_status_size)         ! status from MPL
INTEGER  :: ireq(8)                        ! MPL requests
INTEGER  :: ireq_recv(4)                   ! MPL requests
INTEGER  :: ireq_send(4)                   ! MPL requests
INTEGER  :: nreq                           ! number of requests
INTEGER  :: nreq_recv                      ! number of requests
INTEGER  :: nreq_send                      ! number of requests
INTEGER  :: i,j,k                          ! Spatial loop counters
INTEGER  :: ib, jb                         ! Buffer loop counters
INTEGER  :: recv_loop                      ! recv/unpack loop counter
INTEGER  :: idx                            ! Index of recv message
INTEGER  :: tag                            ! tag from status

! If .TRUE. the North/South halos need to be communicated over the poles
LOGICAL  :: overpoles

! Tags of the north/south messagges
INTEGER :: tag_send_pnorth
INTEGER :: tag_send_psouth
INTEGER :: tag_recv_pnorth
INTEGER :: tag_recv_psouth


! If .TRUE. the element in the corner need to be swapped
LOGICAL  :: sw_cr                          ! South-West corner
LOGICAL  :: se_cr                          ! South-East corner
LOGICAL  :: ne_cr                          ! North-East corner
LOGICAL  :: nw_cr                          ! North-West corner

LOGICAL  :: row_length_even                ! Row length is even
LOGICAL  :: rows_even                      ! The number of rows is even

LOGICAL  :: l_vector

! Offset used to read or write the red or black only elements.
INTEGER  :: send_off_east
INTEGER  :: send_off_west
INTEGER  :: recv_off_east
INTEGER  :: recv_off_west
INTEGER  :: send_off_north
INTEGER  :: send_off_south
INTEGER  :: recv_off_north
INTEGER  :: recv_off_south


INTEGER  :: half_rows
INTEGER  :: half_length

! Number of Red or Black elements in a level of the E/W halo and N/S halo
INTEGER  :: send_num_east    
INTEGER  :: recv_num_east
INTEGER  :: send_num_west    
INTEGER  :: recv_num_west
INTEGER  :: send_num_north
INTEGER  :: recv_num_north
INTEGER  :: send_num_south   
INTEGER  :: recv_num_south

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------
! Set up  

! MPL communication
nreq         = 0
nreq_send    = 0
nreq_recv    = 0
ireq(:)      = mpl_request_null
ireq_send(:) = mpl_request_null
ireq_recv(:) = mpl_request_null

! Sizes
half_rows      = rows/2
half_length    = row_length/2
ew_halo_size   = (half_rows+1)  *levels
ns_halo_size   = (half_length+1)*levels

!---------------------------------------
! Red/black set up

! These are true if the respective dimensions are even
row_length_even = (MOD(row_length,2)==0)
rows_even       = (MOD(rows,2)==0)

! Find if the element in each corner needs to be swapped.
sw_cr = (MOD(g_datastart(1,mype) + g_datastart(2,mype), 2)==0) .EQV. red
se_cr = (sw_cr .NEQV. row_length_even)
nw_cr = (sw_cr .NEQV. rows_even)
ne_cr = (nw_cr .NEQV. row_length_even)

! Compute the offsets for reading and writing the red only or black
! only elements. 
!
! This routine work only with field of type field_type_p. 
!
! In a UM global model, the number of points in the E/W direction,
! global_row_length, is always even. Furthermore, when this number is
! halved the resulting number must be even. 
!
! The receive offsets for the East/West halos are always the opposite
! of the send offsets also in case of cycling E/W boundary. This is
! the case since the global_row_length in UM is always even in the
! global model. The lam model does not have cycling E/W boundary
! so the receive offset is not used at the boundary.
!
! The receive offsets for the North/South halos are the same ones as
! of the send offsets when the data is sent/received over the
! poles. This is because in the global model, the global_row_length is
! always even when it is halved.


send_off_east  = MERGE(0,1,       se_cr)
recv_off_east  = MERGE(0,1, .NOT. se_cr)

send_off_west  = MERGE(0,1,       sw_cr)
recv_off_west  = MERGE(0,1, .NOT. sw_cr)

send_off_south = MERGE(0,1,       sw_cr)
recv_off_south = MERGE(0,1, .NOT. sw_cr)

send_off_north = MERGE(0,1,       nw_cr)
recv_off_north = MERGE(0,1, .NOT. nw_cr)

! Compute the number of red only or black only elements that need to
! be sent or received

send_num_west  = half_rows   + MERGE(1,0,       sw_cr .AND.       nw_cr)
recv_num_west  = half_rows   + MERGE(1,0, .NOT. sw_cr .AND. .NOT. nw_cr)

send_num_east  = half_rows   + MERGE(1,0,       se_cr .AND.       ne_cr)
recv_num_east  = half_rows   + MERGE(1,0, .NOT. se_cr .AND. .NOT. ne_cr)

send_num_south = half_length + MERGE(1,0,       sw_cr .AND.       se_cr)
recv_num_south = half_length + MERGE(1,0, .NOT. sw_cr .AND. .NOT. se_cr)

send_num_north = half_length + MERGE(1,0,       nw_cr .AND.       ne_cr)
recv_num_north = half_length + MERGE(1,0, .NOT. nw_cr .AND. .NOT. ne_cr)


!---------------------------------------
! Over the poles set up 
IF(bound(2)  /=  bc_cyclic .AND. sb_model_type == mt_global) THEN
  overpoles = .TRUE.
ELSE
  overpoles = .FALSE.
END IF

tag_send_pnorth = tag_to_n
tag_recv_pnorth = tag_from_n

tag_send_psouth = tag_to_s
tag_recv_psouth = tag_from_s

IF(overpoles)THEN

  IF (at_extremity(pnorth)) THEN
    
    tag_send_pnorth = tag_to_n_p
    tag_recv_pnorth = tag_from_n_p
    
    ! Same offset and number as the sender
    recv_off_north = send_off_north
    recv_num_north = send_num_north
  END IF

  IF (at_extremity(psouth)) THEN

    tag_send_psouth = tag_to_s_p
    tag_recv_psouth = tag_from_s_p    

    ! Same offset and number as the sender
    recv_off_south = send_off_south
    recv_num_south = send_num_south
  END IF

END IF

!---------------------------------------
! Post MPL receive 

IF (he_neighbour(pwest)  /=  nodomain) THEN
  ! Receive from West
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_w, ew_halo_size,                   &
                 mpl_type,he_neighbour(pwest),tag_from_w,       &
                 halos_comm,ireq_recv(nreq_recv),ierror)
END IF

IF (he_neighbour(peast)  /=  nodomain) THEN
  ! Receive from East
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_e, ew_halo_size,                   &
                 mpl_type,he_neighbour(peast),tag_from_e,       &
                 halos_comm,ireq_recv(nreq_recv),ierror)
END IF

IF (he_neighbour(psouth)  /=  nodomain) THEN
  ! Receive from South
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_s, ns_halo_size,                   &
                 mpl_type,he_neighbour(psouth),tag_recv_psouth, &
                 halos_comm,ireq_recv(nreq_recv),ierror)
END IF

IF (he_neighbour(pnorth)  /=  nodomain) THEN
  ! Receive from North
  nreq_recv=nreq_recv+1
  CALL mpl_irecv(recv_buffer_n, ns_halo_size,                   &
                 mpl_type,he_neighbour(pnorth),tag_recv_pnorth, &
                 halos_comm,ireq_recv(nreq_recv),ierror)
END IF

!---------------------------------------
! East/West halo packing and MPL send

IF (he_neighbour(peast)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, send_num_east, send_off_east, send_buffer_e,     &
!$OMP         field, row_length)                                       &
!$OMP PRIVATE( k, jb, j )
  DO k=1,levels
    DO jb = 1, send_num_east
      j = jb*2 - 1 + send_off_east
      ! Copy data from the Eastern side of the grid
      send_buffer_e(jb,k)=field(row_length,j,k)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Send East
  nreq_send=nreq_send+1
  CALL mpl_isend(send_buffer_e, ew_halo_size,                    &
                 mpl_type,he_neighbour(peast),tag_to_e,          &
                 halos_comm,ireq_send(nreq_send),ierror)
END IF

IF (he_neighbour(pwest)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, send_num_west, send_off_west, send_buffer_w,     &
!$OMP         field)                                                   &
!$OMP PRIVATE( k, jb, j )
  DO k=1,levels
    DO jb = 1, send_num_west
      j = jb*2 - 1 + send_off_west
      ! Copy data from the Western side of the grid
      send_buffer_w(jb,k)=field(1,j,k)
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO

  ! Send West
  nreq_send=nreq_send+1
  CALL mpl_isend(send_buffer_w,ew_halo_size,                    &
                 mpl_type,he_neighbour(pwest),tag_to_w,         &
                 halos_comm,ireq_send(nreq_send),ierror)
END IF


!---------------------------------------
! North/South halo packing and MPL send

IF (he_neighbour(psouth)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, send_num_south, send_off_south, send_buffer_s,   &
!$OMP         field)                                                   &
!$OMP PRIVATE( k, ib, i )
  DO k=1,levels
    DO ib = 1, send_num_south
      i = ib*2 - 1 + send_off_south
      ! Copy data from the South side of the grid
      send_buffer_s(ib,k) = field(i,1,k)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Send South
  nreq_send=nreq_send+1
  CALL mpl_isend(send_buffer_s, ns_halo_size,                   &
                 mpl_type,he_neighbour(psouth),tag_send_psouth, &
                 halos_comm,ireq_send(nreq_send),ierror)
END IF

IF (he_neighbour(pnorth)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, send_num_north, send_off_north, send_buffer_n,   &
!$OMP         field, rows)                                             &
!$OMP PRIVATE( k, ib, i )
  DO k=1,levels
    DO ib = 1, send_num_north
      i = ib*2 - 1 + send_off_north
      ! Copy data from the North side of the grid
      send_buffer_n(ib,k) = field(i,rows,k)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Send North
  nreq_send=nreq_send+1
  CALL mpl_isend(send_buffer_n, ns_halo_size,                   &
                 mpl_type,he_neighbour(pnorth),tag_send_pnorth, &
                 halos_comm,ireq_send(nreq_send),ierror)
END IF

!---------------------------------------
! Copy send and receive request to the general request array
! This is necessary if MPI uses a Rendezvous protocol. In this case a
! mpl_isend might only be started when a mpl_test or mpl_wait is
! called.

nreq = nreq_recv + nreq_send
DO i = 1, nreq_recv
  ireq(i) = ireq_recv(i)
END DO
DO i = 1, nreq_send
  ireq(i+nreq_recv) = ireq_send(i)
END DO

!---------------------------------------
! Wait one-by-one for the halos to arrive and wait also for the sends

DO recv_loop = 1, nreq
  CALL mpl_waitany( nreq, ireq, idx, istat, ierror )

  ! Continue loop if it was a send that completed
  IF (idx>nreq_recv) THEN
    CYCLE
  END IF

  tag = istat(mpl_tag)

  !---------------------------------------
  ! Fill the East/West halos with data

  IF (tag == tag_from_e) THEN
    IF (he_neighbour(peast)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, recv_num_east, recv_off_east, field,             &
!$OMP         recv_buffer_e, row_length)                               &
!$OMP PRIVATE( k, jb, j )
      DO k=1,levels
        DO jb = 1, recv_num_east
          j = jb*2 - 1 + recv_off_east
          field(row_length+1,j,k) = recv_buffer_e(jb,k)
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF 
  END IF

  IF (tag == tag_from_w) THEN
    IF (he_neighbour(pwest)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, recv_num_west, recv_off_west, field,             &
!$OMP         recv_buffer_w)                                           &
!$OMP PRIVATE( k, jb, j )
      DO k=1,levels
        DO jb = 1, recv_num_west
          j = jb*2 - 1 + recv_off_west
          field(0,j,k)=recv_buffer_w(jb,k)
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF 
  END IF

  !---------------------------------------
  ! Fill the North/South halos with data

  IF (tag == tag_from_s .OR. tag == tag_from_s_p) THEN
    IF (he_neighbour(psouth)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, recv_num_south, recv_off_south, field,           &
!$OMP         recv_buffer_s)                                           &
!$OMP PRIVATE( k, ib, i )
      DO k=1,levels
        DO ib = 1, recv_num_south
          i = ib*2 - 1 + recv_off_south
          field(i,0,k) = recv_buffer_s(ib,k)
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF 
  END IF

  IF (tag == tag_from_n .OR. tag == tag_from_n_p) THEN
    IF (he_neighbour(pnorth)  /=  nodomain) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP SHARED( levels, recv_num_north, recv_off_north, field,           &
!$OMP         recv_buffer_n,rows)                                      &
!$OMP PRIVATE( k, ib, i )
      DO k=1,levels
        DO ib = 1, recv_num_north
          i = ib*2 - 1 + recv_off_north
          field(i,rows+1,k) = recv_buffer_n(ib,k)
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF
  END IF

END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end halo_exchange_mpi_mod_swap_bounds_mpi_rb.h --
