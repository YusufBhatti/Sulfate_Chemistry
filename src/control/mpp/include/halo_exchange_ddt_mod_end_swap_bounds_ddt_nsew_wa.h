! -- begin halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h --
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

#if   defined(HALO_EXCHANGE_DDT_INTEGER)
INTEGER, INTENT(INOUT) ::                                                  &
#elif defined(HALO_EXCHANGE_DDT_LOGICAL)
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

! The derived that types. This is expected to be already initialised.
TYPE(ddt_type), INTENT(IN) :: ddt 

INTEGER, INTENT(INOUT) :: ns_halos_recv          ! Number of N/S halos received
INTEGER, INTENT(INOUT) :: nreq_recv              ! number of requests
INTEGER, INTENT(INOUT) :: ireq_recv(4)           ! MPL request

INTEGER, INTENT(IN) :: tag_base ! Used to identify unique tags


! Local variables

INTEGER  :: ierr                    ! MPI error code
INTEGER  :: istat_recv(mpl_status_size) 
INTEGER  :: m,idx


INTEGER  :: i,j,k                     ! Spatial loop counters
LOGICAL  :: change_sign               ! .TRUE. if sign change across pole
LOGICAL  :: overpole_south            ! halos over the south pole
LOGICAL  :: overpole_north            ! halos over the north pole    
LOGICAL  :: edge_north                ! copy edge into north halo
LOGICAL  :: edge_south                ! copy edge into south halo

LOGICAL  :: send_south, send_north, send_east, send_west
LOGICAL  :: recv_south, recv_north, recv_east, recv_west

INTEGER  ::                &
  north_off,               & ! Offsets to use when copying data
  south_off                  ! to send around poles

! Tags of the north/south messagges
INTEGER :: tag_send_pnorth
INTEGER :: tag_send_psouth
INTEGER :: tag_recv_pnorth
INTEGER :: tag_recv_psouth

! If .TRUE. the North/South halos need to be communicated over the poles
LOGICAL  :: overpoles

REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Over the poles set up

IF(bound(2)  /=  bc_cyclic .AND. sb_model_type == mt_global) THEN
  overpoles = .TRUE.
ELSE
  overpoles = .FALSE.
END IF

recv_south = do_south_arg
send_north = do_south_arg

recv_north = do_north_arg
send_south = do_north_arg

recv_east  = do_east_arg
send_west  = do_east_arg

recv_west  = do_west_arg
send_east  = do_west_arg

north_off=0
south_off=0
  
tag_send_pnorth = tag_to_n
tag_recv_pnorth = tag_from_n

tag_send_psouth = tag_to_s
tag_recv_psouth = tag_from_s

overpole_south = .FALSE.
overpole_north = .FALSE.
change_sign    = .FALSE.

IF (overpoles) THEN

  ! Set up the offsets. When copying data that is to be passed over
  ! the pole, on wind (u or v) grid, then copy data one row away
  ! from the pole
  
  IF ( field_type  ==  fld_type_v) THEN
    IF (at_extremity(pnorth)) north_off=1
    IF (at_extremity(psouth)) south_off=1
  END IF
  
  ! Set up the sign factor. If L_VECTOR is true and data has been passed
  ! over the poles in a global model, then the variables must change
  ! sign
  
  IF (.NOT. ((l_vector))) THEN
    change_sign=.FALSE.
  ELSE
    change_sign=.TRUE.
  END IF
  
  ! If the processor is at the poles, change the tags

  IF (at_extremity(pnorth)) THEN
    
    tag_send_pnorth = tag_to_n_p
    tag_recv_pnorth = tag_from_n_p
    
    overpole_north = .TRUE.
    
    recv_north = do_north_arg
    send_north = do_north_arg

  END IF

  IF (at_extremity(psouth)) THEN

    tag_send_psouth = tag_to_s_p
    tag_recv_psouth = tag_from_s_p
    
    overpole_south = .TRUE.

    recv_south = do_south_arg
    send_south = do_south_arg
    
  END IF

END IF

!--------------------------------------------------------------------
! Wait for any message and fill the corresponding halo. Then copy and
! send East and West halos when both North and South are received.

DO m = 1, nreq_recv
  CALL mpl_waitany(nreq_recv, ireq_recv, idx, istat_recv, ierr)

  ! If received from East or West do nothing

  ! Received from North
  IF (istat_recv(mpl_tag) == tag_base + tag_recv_pnorth) THEN

    IF (overpole_north .AND. change_sign) THEN
#if !defined(HALO_EXCHANGE_DDT_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, rows, row_length)
      DO k=1,levels
        DO j=1,halo_y
          DO i=1,row_length
            field(i,rows+j,k)= -field(i,rows+j,k)                
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO 
#endif
    END IF 

    ! Increase the number of N/S halos received
    ns_halos_recv = ns_halos_recv + 1

  END IF

  ! Received from South
  IF (istat_recv(mpl_tag) == tag_base + tag_recv_psouth) THEN

    IF (overpole_south .AND. change_sign) THEN
#if !defined(HALO_EXCHANGE_DDT_LOGICAL)
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP SHARED(levels, halo_y, halo_x, field, rows, row_length)
      DO k=1,levels
        DO j=1,halo_y
          DO i=1,row_length
            field(i,1-j,k)= -field(i,1-j,k)                      
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO 
#endif
    END IF 

    ! Increase the number of N/S halos received
    ns_halos_recv = ns_halos_recv + 1

  END IF

  ! If we have received both North and South halos. Then copy and
  ! send the East and West halos.
  IF (ns_halos_recv==2) THEN

    !------------------------------------------------------------------
    ! East-West communications
    
    IF (he_neighbour(pwest)  /=  nodomain .AND. send_west) THEN

      ! Send West
      CALL mpl_send(field(1,1-halo_y,1), levels,                               &
        ddt%ew_halo_ext_type, he_neighbour(pwest), tag_base + tag_to_w,        &
        halos_comm, ierr)
      

    END IF

    IF (he_neighbour(peast)  /=  nodomain .AND. send_east) THEN

      ! Send East
      CALL mpl_send(field(row_length-halo_x+1,1-halo_y,1), levels,             &
        ddt%ew_halo_ext_type, he_neighbour(peast), tag_base + tag_to_e,        &
        halos_comm, ierr)

    END IF

    ! We want to do this only once
    ns_halos_recv = -1

  END IF

END DO

!--------------------------------------------------------------------
! Special cases:

! he_neighbour(psouth)==mype and he_neighbour(pnorth)==mype
! 
! Copy the data of the North/South halos directly if:
!
! a) Overpoles communication with only one single PE in the E/W
! b) Cyclic (no overpoles) communication on the Y direction with only
!    one single PE in the N/S

IF(                                                               &
  (nproc_x == 1 .AND. overpoles) .OR.                             &
  (nproc_y == 1 .AND. bound(2) == bc_cyclic )                     &
  ) THEN

  CALL copy_ns_halo_single(field,row_length,rows,levels,          &
    halo_x,halo_y,                                                &
    north_off, south_off, overpoles,                              &
    overpole_north, overpole_south, change_sign)

END IF

! EW cyclic LAMs
! Copy adjacent edge data at the North and South extremities into halos

IF (sb_model_type == mt_cyclic_lam) THEN

    edge_north = .FALSE.
    edge_south = .FALSE.

  IF(at_extremity(pnorth) .AND. he_neighbour(pnorth)  ==  nodomain) THEN

    edge_north = .TRUE.
    
  END IF

  IF(at_extremity(psouth) .AND. he_neighbour(psouth)  ==  nodomain) THEN

    edge_south = .TRUE.

  END IF

  CALL copy_edge_to_ns_halo(field,row_length,rows,levels,         &
    halo_x,halo_y, edge_north, edge_south)
  
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

! -- end halo_exchange_ddt_mod_end_swap_bounds_ddt_nsew_wa.h --
