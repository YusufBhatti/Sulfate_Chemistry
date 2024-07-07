! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Direct MPI version of swap_bounds to deal with multiple variables
! at once

MODULE swap_bounds_mv_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SWAP_BOUNDS_MV_MOD'

CONTAINS

SUBROUTINE swap_bounds_mv(  &
   input_fields,            &        ! Fields to be swapped
   n_multi,                 &        ! The number of Fields to swap
   row_length,              &        ! field size
   halo_x, halo_y)                   ! halos

USE mpl, ONLY:              &
         mpl_request_null,  &
         mpl_real,          &
         mpl_status_size

USE swapable_field_mod, ONLY: &
    swapable_field_pointer_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParParams, ONLY: bc_cyclic, nodomain, peast, pwest, pnorth, psouth
USE UM_ParCore, ONLY: mype
USE field_types, ONLY: fld_type_v
USE model_domain_mod, ONLY: mt_global, mt_lam

IMPLICIT NONE

!  This code performs the same process as swap_bounds but
!  can swap multiple variables at once and hence utilise bandwidth
!  better. This is the MPL version.

!  Note that it can only deal with multiple variables with the same
!  row_length and same halo size.

! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.

! Implementation
!   The logic flow is non-trivial!
!   The across pole differencing in particular must be handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.
!   The North/South halos are done first, then the East/West
!   halos.

!   Note that due to the fact that pointers to the data are used,
!   addressing is (1:row_length+2*halo_x, 1:rows+2*halo_y) rather
!   than (1-halo_x:row_length+halo_x, 1-halo_y:rows+halo_y) as used
!   elsewhere. This is unavoidable as pointers change the addressing
!   mode.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Arguments:

INTEGER, INTENT(IN)  :: n_multi      ! number of fields to be swapped
INTEGER, INTENT(IN)  :: row_length   ! number of points on row (no halos)
INTEGER, INTENT(IN)  :: halo_x       ! sixze of "i" halo
INTEGER, INTENT(IN)  :: halo_y       ! sixze of "j" halo

! Fields to swap
TYPE(swapable_field_pointer_type), TARGET :: input_fields(n_multi)

! Local scalar variables
INTEGER :: levels                ! number of levels in field
INTEGER :: rows                  ! number of rows in field (no halos)
INTEGER :: field_type            ! The grid type of the field (u,v,p)
INTEGER :: i,j,k                 ! Spatial loop counters
INTEGER :: info                  ! GCOM return code
INTEGER :: length
INTEGER :: i_field
INTEGER :: ierror                ! MPI return code
INTEGER :: nreq_r_ew             ! number of MPI recv requests ew
INTEGER :: nreq_s_ew             ! number of MPI send requests ew
INTEGER :: nreq_r_ns             ! number of MPI recv requests ns
INTEGER :: nreq_s_ns             ! number of MPI send requests ns
INTEGER :: full_row_length       ! length of row including halos
INTEGER :: full_rows             ! number of rows including halos
INTEGER :: ew_halo_size          ! size of EW halo region
INTEGER :: ns_halo_size          ! size of NS halo region
INTEGER :: west_halo_source      ! first column of data for west halo
INTEGER :: east_halo_source      ! first column of data for east halo
INTEGER :: south_halo_source
INTEGER :: north_halo_source
INTEGER :: half_full_row_length  ! half of the total EW dimension
INTEGER :: half_row_length       ! half of the data (no halo) EW dim
INTEGER :: max_levels            ! maximum levels used over all fields
INTEGER :: max_rows              ! maximum rows used over all fields
INTEGER :: my_comm               ! Communicator

LOGICAL :: l_vector              ! TRUE if a vector field


! Local arrays
INTEGER :: istat(mpl_status_size,4)  ! MPI status
INTEGER :: ireq_r_ew(4)              ! MPI requests
INTEGER :: ireq_s_ew(4)              ! MPI requests
INTEGER :: ireq_r_ns(4)              ! MPI requests
INTEGER :: ireq_s_ns(4)              ! MPI requests
INTEGER :: north_off(n_multi)        ! Offsets to use when copying data
INTEGER :: south_off(n_multi)        ! to send around poles

REAL, POINTER :: field(:, :, :)

! Send and receive buffers to be allocated dynamically
REAL, ALLOCATABLE :: send_buffer_e(:,:,:,:)
REAL, ALLOCATABLE :: send_buffer_w(:,:,:,:)
REAL, ALLOCATABLE :: recv_buffer_e(:,:,:,:)
REAL, ALLOCATABLE :: recv_buffer_w(:,:,:,:)
REAL, ALLOCATABLE :: send_buffer_n(:,:,:,:)
REAL, ALLOCATABLE :: send_buffer_s(:,:,:,:)
REAL, ALLOCATABLE :: recv_buffer_n(:,:,:,:)
REAL, ALLOCATABLE :: recv_buffer_s(:,:,:,:)

LOGICAL :: change_sign(n_multi)     ! .TRUE. if sign change across pole



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS_MV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Check if there is anything to do
IF (((halo_x == 0) .AND. (halo_y == 0)) .OR. n_multi == 0) GO TO 9999

!------------------------------------------------------------------
! Initialise variables

ireq_r_ns(:) = mpl_request_null
ireq_s_ns(:) = mpl_request_null
ireq_r_ew(:) = mpl_request_null
ireq_s_ew(:) = mpl_request_null

! Maximum rows and levels
max_rows   = 0
max_levels = 0
DO i_field=1, n_multi
  max_levels= MAX(input_fields(i_field) % levels, max_levels)
  max_rows  = MAX(input_fields(i_field) % rows,   max_rows)
END DO

full_row_length = row_length + 2*halo_x
full_rows       = max_rows   + 2*halo_y

half_full_row_length = full_row_length/2
half_row_length      = row_length/2

ew_halo_size         = full_rows * halo_x
ns_halo_size         = row_length * halo_y



! Get communicator we will be using from GCOM
CALL gc_get_communicator(my_comm,ierror)


!------------------------------------------------------------------
!  North-South communications
!  section of code for cyclic N-S boundary conditions
IF (halo_y > 0) THEN

  ALLOCATE( send_buffer_n(row_length, halo_y, max_levels, n_multi) )
  ALLOCATE( send_buffer_s(row_length, halo_y, max_levels, n_multi) )
  ALLOCATE( recv_buffer_n(row_length, halo_y, max_levels, n_multi) )
  ALLOCATE( recv_buffer_s(row_length, halo_y, max_levels, n_multi) )

  IF (bound(2) == bc_cyclic) THEN

    IF (nproc_y == 1) THEN ! only 1 processor north-south
      ! 2.1 Simple case of only one processor

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, i, j,  &
!$OMP& k, i_field, south_halo_source, north_halo_source, levels, rows)  &
!$OMP& SHARED(halo_y, halo_x, input_fields, n_multi, row_length)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        south_halo_source=rows+1              ! copy from opposite end
        north_halo_source=halo_y + 1          ! of each column

        DO k=1,levels
          DO j=1,halo_y
            DO i=1 ,row_length

              ! Fill southern halo
              field(i+halo_x,j,k)=field(i+halo_x,south_halo_source+j-1,k)

              ! Fill northern halo
              field(i+halo_x,rows+halo_y+j,k)=&
                field(i+halo_x,north_halo_source+j-1,k)

            END DO ! I
          END DO ! J
        END DO ! K
      END DO ! loop over fields
!$OMP END PARALLEL DO

        !---------------------------------------
        ! Now the more common case of having
        ! a number of processors in the
        ! North-South direction

    ELSE ! If there is more than 1 processor north_south

      !---------------------------------------
      ! Copy the data into buf_send

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, i_field,  &
!$OMP& k, levels, rows, i, j) SHARED(row_length, input_fields, halo_y,    &
!$OMP& n_multi, send_buffer_s, send_buffer_n, halo_x)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        DO k=1,levels
          DO j=1,halo_y
            DO i=1, row_length

              ! Copy stuff from the southern side of the grid
              send_buffer_s(i,j,k,i_field) = field(i+halo_x,halo_y+j,k)

              ! Copy stuff from the northern side of the grid
              send_buffer_n(i,j,k,i_field) = field(i+halo_x,rows+j,k)

            END DO ! I
          END DO ! J
        END DO ! K
      END DO ! loop over fields
!$OMP END PARALLEL DO

        !---------------------------------------
        ! Send and receive the data

      length = n_multi * ns_halo_size * max_levels

      IF (neighbour(psouth) /= nodomain) THEN

        ! Send south
        CALL gc_rsend(2,length,neighbour(psouth),info,                    &
                       recv_buffer_n, send_buffer_s)
      END IF

      IF (neighbour(pnorth) /= nodomain) THEN

        ! Send north
        CALL gc_rsend(3,length,neighbour(pnorth),info,                    &
                       recv_buffer_s, send_buffer_n)
      END IF

      IF (neighbour(pnorth) /= nodomain) THEN

        ! Receive from north
        CALL gc_rrecv(2,length,neighbour(pnorth),info,                    &
                       recv_buffer_n, send_buffer_s)
      END IF

      IF (neighbour(psouth) /= nodomain) THEN

        ! Receive from south
        CALL gc_rrecv(3,length,neighbour(psouth),info,                    &
                       recv_buffer_s, send_buffer_n)

      END IF

      !---------------------------------------
      ! Fill the halos with data

      IF (neighbour(pnorth) /= nodomain) THEN

        ! unpack data from receive_buffer into field

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, i_field, &
!$OMP& k, levels, rows, i, j) SHARED(input_fields, n_multi, halo_x,       &
!$OMP& halo_y, recv_buffer_n, row_length)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          rows        = input_fields(i_field) % rows

          DO k=1,levels
            DO j=1,halo_y
              DO i=1,row_length
                field(i+halo_x,j+halo_y+rows,k) = recv_buffer_n(i,j,k,i_field)
              END DO
            END DO
          END DO
        END DO ! loop over fields
!$OMP END PARALLEL DO

      END IF ! IF (neighbour(Pnorth) /= NoDomain)

      IF (neighbour(psouth) /= nodomain) THEN

        ! unpack data from receive_buffer into field

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,  &
!$OMP& i, j, k, i_field) SHARED(row_length, input_fields, n_multi,        &
!$OMP& halo_x, halo_y, recv_buffer_s)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          DO k=1,levels
            DO j=1,halo_y
              DO i=1,row_length
                field(i+halo_x,j,k) = recv_buffer_s(i,j,k,i_field)
              END DO
            END DO
          END DO
        END DO ! loop over fields
!$OMP END PARALLEL DO


      END IF ! IF (neighbour(Psouth) /= NoDomain)

    END IF ! IF (nproc_y == 1)

  ELSE                 !!! bc_cyclic in NS

    ! Set up some variables

    ! Set up the offsets. When copying data that is to be passed over
    ! the pole, on wind (u or v) grid, then copy data one row away
    ! from the pole

    ! parameters: fld_type_v, fld_type_p, fld_type_u, mt_global
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i_field,        &
!$OMP& l_vector, field_type) SHARED(n_multi, north_off, south_off,        &
!$OMP& change_sign, sb_model_type, at_extremity, input_fields)

    DO i_field=1, n_multi
      north_off(i_field)=0
      south_off(i_field)=0
      field_type = input_fields(i_field) % field_type
      IF ((sb_model_type == mt_global) .AND.                                &
           (field_type  ==  fld_type_v)) THEN
        IF (at_extremity(pnorth)) north_off(i_field)=1
        IF (at_extremity(psouth)) south_off(i_field)=1
      END IF

      ! Set up the sign factor. If l_vector is true and data has been passed
      ! over the poles in a global model, then the variables must change
      ! sign

      l_vector = input_fields(i_field) % vector
      IF (.NOT. ((l_vector) .AND. (sb_model_type == mt_global))) THEN
        change_sign(i_field)=.FALSE.
      ELSE
        change_sign(i_field)=.TRUE.
      END IF
    END DO
!$OMP END PARALLEL DO

      !---------------------------------------
      ! Copy data into the send_buffer
      ! But not if:
      !   - Not a global model and the data is at the North/South edge
      !   - A global model at the North/South edge but only 1 processor
      !     in the East-West direction.

    IF (.NOT. (at_extremity(psouth) .AND.                                   &
                ((nproc_x == 1)  .OR.                                       &
                (sb_model_type /= mt_global))))                             &
       THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k) SHARED(input_fields, halo_y, halo_x, send_buffer_s,&
!$OMP& n_multi, south_off, row_length)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        DO k=1,levels
          DO j=1,halo_y
            DO i=1,row_length

              send_buffer_s(i,j,k,i_field) =                    &
                 field(i+halo_x,j+halo_y+south_off(i_field),k)

            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF

    IF (.NOT. (at_extremity(pnorth) .AND.                                   &
                ((nproc_x == 1) .OR.                                        &
                (sb_model_type /= mt_global))))                             &
       THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k, rows) SHARED(row_length, halo_x, halo_y,         &
!$OMP& send_buffer_n, input_fields, n_multi, north_off)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        DO k=1,levels
          DO j=1,halo_y
            DO i=1,row_length

              send_buffer_n(i,j,k,i_field) =                    &
                 field(i+halo_x,rows-north_off(i_field)+j,k)

            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF

    !---------------------------------------
    ! Send and receive the data

    !---------------------------------------
    !---------------------------------------
    ! The more general case, where
    ! each buffer is sent to a
    ! different processor

    nreq_r_ns=0

    length = n_multi * ns_halo_size * max_levels

    IF (at_extremity(psouth)) THEN

      IF ((neighbour(psouth) /= nodomain) .AND.                             &
           (neighbour(psouth) /= mype)) THEN

        nreq_r_ns=nreq_r_ns+1
        CALL mpl_irecv(recv_buffer_s,                                      &
          length, mpl_real, neighbour(psouth), 11, my_comm,                &
          ireq_r_ns(nreq_r_ns), ierror)

      END IF

    ELSE ! not at the South

      nreq_r_ns=nreq_r_ns+1
      CALL mpl_irecv(recv_buffer_s,                                        &
        length, mpl_real, neighbour(psouth), 14, my_comm,                  &
        ireq_r_ns(nreq_r_ns), ierror)

    END IF

    IF (at_extremity(pnorth)) THEN

      IF ((neighbour(pnorth) /= nodomain) .AND.                             &
           (neighbour(pnorth) /= mype)) THEN

        nreq_r_ns=nreq_r_ns+1
        CALL mpl_irecv(recv_buffer_n,                                      &
          length, mpl_real, neighbour(pnorth), 13, my_comm,                &
          ireq_r_ns(nreq_r_ns), ierror)

      END IF

    ELSE

      nreq_r_ns=nreq_r_ns+1
      CALL mpl_irecv(recv_buffer_n,                                      &
        length, mpl_real, neighbour(pnorth), 12, my_comm,                &
        ireq_r_ns(nreq_r_ns), ierror)

    END IF

    nreq_s_ns=0
    IF (at_extremity(psouth)) THEN

      IF ((neighbour(psouth) /= nodomain) .AND.                             &
           (neighbour(psouth) /= mype)) THEN

        nreq_s_ns=nreq_s_ns+1
        CALL mpl_isend(send_buffer_s,                                      &
          length, mpl_real, neighbour(psouth), 11, my_comm,                &
          ireq_s_ns(nreq_s_ns),ierror)

      END IF

    ELSE ! not at the South

      nreq_s_ns=nreq_s_ns+1
      CALL mpl_isend(send_buffer_s,                                      &
        length, mpl_real, neighbour(psouth), 12, my_comm,                &
        ireq_s_ns(nreq_s_ns),ierror)

    END IF

    IF (at_extremity(pnorth)) THEN

      IF ((neighbour(pnorth) /= nodomain) .AND.                             &
           (neighbour(pnorth) /= mype)) THEN

        nreq_s_ns=nreq_s_ns+1
        CALL mpl_isend(send_buffer_n,                                        &
          length, mpl_real, neighbour(pnorth), 13, my_comm,                  &
          ireq_s_ns(nreq_s_ns),ierror)

      END IF

    ELSE ! not at the North

      nreq_s_ns=nreq_s_ns+1
      CALL mpl_isend(send_buffer_n,                                        &
        length, mpl_real, neighbour(pnorth), 14, my_comm,                  &
        ireq_s_ns(nreq_s_ns),ierror)

    END IF

    CALL mpl_waitall( nreq_r_ns, ireq_r_ns, istat, ierror )
    !---------------------------------------
    ! Fill the halos with data

    !---------------------------------------
    ! Southern halo

    IF (at_extremity(psouth)) THEN

      IF (neighbour(psouth)  ==  nodomain) THEN
        IF (sb_model_type /= mt_lam) THEN
          ! Only see this section in EW cyclic LAMs
          ! Just copy adjacent rows into halo area
          DO i_field = 1, n_multi
            field   => input_fields(i_field) % field
            levels  =  input_fields(i_field) % levels

            DO k = 1, levels
              DO j = 1, halo_y
                DO i = 1,row_length+(2*halo_x)
                  field(i,j,k) = field(i,halo_y+1,k)
                END DO
              END DO
            END DO
          END DO
        END IF

      ELSE IF (neighbour(psouth) == mype) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k) SHARED(input_fields, change_sign, halo_x,        &
!$OMP& halo_y, n_multi, half_full_row_length, half_row_length, south_off)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          IF (change_sign(i_field)) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length
                  field(half_row_length+halo_x+i,1-j+halo_y,k)=             &
                     -field(i+halo_x,j+halo_y+south_off(i_field),k)
                  field(i,1-j+halo_y,k)=                                    &
                     -field(half_row_length+i,j+halo_y+                     &
                     south_off(i_field),k)
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length
                  field(half_row_length+halo_x+i,1-j+halo_y,k)=             &
                     field(i+halo_x,j+halo_y+south_off(i_field),k)
                  field(i,1-j+halo_y,k)=                                    &
                     field(half_row_length+i,j+halo_y+                      &
                     south_off(i_field),k)
                END DO
              END DO
            END DO
          END IF ! IF (change_sign(i_field))
        END DO ! loop over field
!$OMP END PARALLEL DO

      ELSE IF (neighbour(psouth) /= nodomain) THEN ! not LAM

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k) SHARED(input_fields, halo_x, halo_y, row_length, &
!$OMP& n_multi, recv_buffer_s, change_sign)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          IF (change_sign(i_field)) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,row_length
                  field(i+halo_x,1-j+halo_y,k)=                             &
                  -recv_buffer_s(i,j,k,i_field)
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,row_length
                  field(i+halo_x,1-j+halo_y,k)=                             &
                   recv_buffer_s(i,j,k,i_field)
                END DO
              END DO
            END DO
          END IF !  IF (change_sign(i_field))
        END DO ! loop over fields
!$OMP END PARALLEL DO

      END IF ! What type of South extremity

    ELSE ! IF (at_extremity(PSouth)

      ! not at a South extremity

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels, &
!$OMP& i_field, i, j, k) SHARED(recv_buffer_s, halo_x, input_fields,   &
!$OMP& n_multi, row_length, halo_y)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        DO k=1,levels
          DO j=1,halo_y
            DO i=1,row_length
              field(i+halo_x,j,k) = recv_buffer_s(i,j,k,i_field)
            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF ! IF (at_extremity(PSouth)


    !---------------------------------------
    ! Northern halo

    IF (at_extremity(pnorth)) THEN

      IF (neighbour(pnorth)  ==  nodomain) THEN
        IF (sb_model_type /= mt_lam) THEN
          ! Only see this section in EW cyclic LAMs
          ! Just copy adjacent rows into halo area
          DO i_field = 1, n_multi
            field   => input_fields(i_field) % field
            levels  =  input_fields(i_field) % levels
            rows    =  input_fields(i_field) % rows

            DO k = 1, levels
              DO j = 1, halo_y
                DO i = 1,row_length+(2*halo_x)
                  field(i,rows+halo_y+j,k)=field(i,rows+halo_y,k)
                END DO
              END DO
            END DO
          END DO
        END IF

      ELSE IF (neighbour(pnorth) == mype) THEN
        ! Local across pole difference

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k, rows) SHARED(input_fields, half_full_row_length, &
!$OMP& half_row_length, n_multi, north_off, halo_y, halo_x, change_sign)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          rows        = input_fields(i_field) % rows
          IF (change_sign(i_field)) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length

                  field(half_row_length+halo_x+i,rows+halo_y+j,k)=          &
                     -field(i+halo_x,rows+halo_y-j+1-north_off(i_field),k)
                  field(i,rows+j+halo_y,k)=                                 &
                     -field(half_row_length+i,                              &
                          rows-j+halo_y+1-north_off(i_field),k)

                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,half_full_row_length

                  field(half_row_length+halo_x+i,rows+halo_y+j,k)=          &
                     field(i+halo_x,rows+halo_y-j+1-north_off(i_field),k)
                  field(i,rows+j+halo_y,k)=                                 &
                     field(half_row_length+i,                               &
                         rows-j+1+halo_y-north_off(i_field),k)

                END DO
              END DO
            END DO
          END IF ! IF (change_sign(i_field))
        END DO ! loop over fields
!$OMP END PARALLEL DO

      ELSE IF (neighbour(pnorth) /= nodomain) THEN ! not LAM

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels, &
!$OMP& i_field, i, j, k, rows) SHARED(input_fields, change_sign, halo_y,&
!$OMP& halo_x, recv_buffer_n, row_length, n_multi)
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          rows        = input_fields(i_field) % rows
          IF (change_sign(i_field)) THEN
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,row_length
                  field(i+halo_x,rows+halo_y+j,k)=                          &
                     -recv_buffer_n(i,halo_y-j+1,k,i_field)
                END DO
              END DO
            END DO
          ELSE ! don't change sign
            DO k=1,levels
              DO j=1,halo_y
                DO i=1,row_length
                  field(i+halo_x,rows+halo_y+j,k)=                          &
                     recv_buffer_n(i,halo_y-j+1,k,i_field)
                END DO
              END DO
            END DO
          END IF ! IF (change_sign(i_field))
        END DO ! loop over fields
!$OMP END PARALLEL DO

      END IF ! What type of North extremity

    ELSE ! IF (at_extremity(PNorth)

      ! not at a North extremity

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, levels,   &
!$OMP& i_field, i, j, k, rows) SHARED(n_multi, halo_y, halo_x, row_length,&
!$OMP& recv_buffer_n, input_fields)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        DO k=1,levels
          DO j=1,halo_y
            DO i=1,row_length
              field(i+halo_x,rows+halo_y+j,k) = recv_buffer_n(i,j,k,i_field)
            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF ! IF (at_extremity(PNorth)

    ! This need to be inside the bound(2) /= bc_cyclic branch
    CALL mpl_waitall( nreq_s_ns, ireq_s_ns, istat, ierror )

  END IF ! BC_CYCLIC  for NS

  IF (ALLOCATED(send_buffer_n))  DEALLOCATE(send_buffer_n)
  IF (ALLOCATED(send_buffer_s))  DEALLOCATE(send_buffer_s)
  IF (ALLOCATED(recv_buffer_n))  DEALLOCATE(recv_buffer_n)
  IF (ALLOCATED(recv_buffer_s))  DEALLOCATE(recv_buffer_s)
END IF ! halo_y > 0

!------------------------------------------------------------------
! East-West communications
!---------------------------------------
! Simple case of only one processor
! in the East-West direction

IF (halo_x > 0) THEN
  ! Allocate buffers for communication
  ALLOCATE( send_buffer_e(halo_x, full_rows, max_levels, n_multi) )
  ALLOCATE( send_buffer_w(halo_x, full_rows, max_levels, n_multi) )
  ALLOCATE( recv_buffer_e(halo_x, full_rows, max_levels, n_multi) )
  ALLOCATE( recv_buffer_w(halo_x, full_rows, max_levels, n_multi) )

  nreq_s_ew = 0
  nreq_s_ns = 0
  IF (nproc_x == 1) THEN ! only 1 processor East-West

    IF (bound(1) == bc_cyclic) THEN        ! cyclic boundary conditions
      west_halo_source=row_length+1        ! copy from opposite end
      east_halo_source=halo_x+1            ! of each row

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(field,         &
!$OMP& k, i_field, levels, rows, i, j) SHARED(input_fields, halo_x,      &
!$OMP& n_multi, halo_y, east_halo_source, west_halo_source, row_length)
      DO i_field=1, n_multi
        field       => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        !CDIR NOVECTOR
        DO i=1,halo_x
          DO k=1,levels
            DO j=1,rows + (2 * halo_y)

              ! Fill Western halo
              field(i,j,k) = field(west_halo_source+i-1,j,k)

              ! Fill Eastern halo
              field(row_length+halo_x+i,j,k) =                              &
                              field(east_halo_source+i-1,j,k)

            END DO ! J
          END DO ! K
        END DO ! I

      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF   !  bound(1) == BC_CYCLIC

    !---------------------------------------
    ! Now the more common case of having
    ! a number of processors in the
    ! East-West direction

  ELSE ! If there is more than 1 processor East-West

    !---------------------------------------
    ! Copy the data into send_buffer

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, k,       &
!$OMP& i_field, levels, rows, i, j) SHARED(input_fields, row_length,      &
!$OMP& send_buffer_e, send_buffer_w, halo_x, halo_y, n_multi)
    DO i_field=1, n_multi
      field       => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO k=1,levels
        DO j=1,rows+2*halo_y
          DO i=1,halo_x

            ! Copy stuff from the Western side of the grid
            send_buffer_w(i,j,k,i_field) = field(halo_x+i,j,k)

            ! Copy stuff from the Eastern side of the grid
            send_buffer_e(i,j,k,i_field) = field(row_length+i,j,k)

          END DO ! I
        END DO ! J
      END DO ! K
    END DO ! loop over fields
!$OMP END PARALLEL DO

    !---------------------------------------
    ! Send and receive the data

    !---------------------------------------
    length = n_multi * ew_halo_size * max_levels

    nreq_r_ew=0
    IF (neighbour(peast) /= nodomain) THEN

      ! Receive from East
      nreq_r_ew=nreq_r_ew+1
      CALL mpl_irecv(recv_buffer_e,                                       &
                     length, mpl_real, neighbour(peast), 2,               &
                     my_comm, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    IF (neighbour(pwest) /= nodomain) THEN

      ! Receive from West
      nreq_r_ew=nreq_r_ew+1
      CALL mpl_irecv(recv_buffer_w,                                       &
                     length, mpl_real, neighbour(pwest), 3,               &
                     my_comm, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    nreq_s_ew=0
    IF (neighbour(peast) /= nodomain) THEN

      ! Send East
      nreq_s_ew=nreq_s_ew+1
      CALL mpl_isend(send_buffer_e,                                       &
                     length, mpl_real, neighbour(peast), 3,               &
                     my_comm, ireq_s_ew(nreq_s_ew), ierror)
    END IF

    IF (neighbour(pwest) /= nodomain) THEN

      ! Send West
      nreq_s_ew=nreq_s_ew+1
      CALL mpl_isend(send_buffer_w,                                       &
                     length, mpl_real, neighbour(pwest), 2,               &
                     my_comm, ireq_s_ew(nreq_s_ew), ierror)

    END IF

    CALL mpl_waitall( nreq_r_ew, ireq_r_ew, istat, ierror )


    CALL mpl_waitall( nreq_s_ew, ireq_s_ew, istat, ierror )

    !---------------------------------------
    ! Fill the halos with data

    IF (neighbour(peast) /= nodomain) THEN

      ! unpack data from receive_buffer into field

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field, i_field, &
!$OMP& k, levels, rows, i, j) SHARED(recv_buffer_e, halo_x,               &
!$OMP& row_length, halo_y, input_fields, n_multi)
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        DO k=1,levels
          DO j=1,rows+2*halo_y
            DO i=1,halo_x
              field(row_length+halo_x+i,j,k) = recv_buffer_e(i,j,k,i_field)
            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF ! IF (neighbour(PEast) /= NoDomain)

    IF (neighbour(pwest) /= nodomain) THEN

      ! unpack data from receive_buffer into field

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(field,        &
!$OMP& k, i_field,  levels, rows, i, j) SHARED(input_fields, halo_y,    &
!$OMP& n_multi, recv_buffer_w, halo_x )
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        DO k=1,levels
          DO j=1,rows+2*halo_y
            DO i=1,halo_x
              field(i,j,k) = recv_buffer_w(i,j,k,i_field)
            END DO
          END DO
        END DO
      END DO ! loop over fields
!$OMP END PARALLEL DO

    END IF ! IF (neighbour(PWest) /= NoDomain)

  END IF ! IF (nproc_x == 1)

  IF (ALLOCATED(send_buffer_e))  DEALLOCATE(send_buffer_e)
  IF (ALLOCATED(send_buffer_w))  DEALLOCATE(send_buffer_w)
  IF (ALLOCATED(recv_buffer_e))  DEALLOCATE(recv_buffer_e)
  IF (ALLOCATED(recv_buffer_w))  DEALLOCATE(recv_buffer_w)
END IF ! halo_x > 0

9999 CONTINUE


! Fill external halos for LAMs
IF (sb_model_type == mt_lam) THEN
  ! DEPENDS ON: fill_external_halos
  DO i_field = 1, n_multi
    field  => input_fields(i_field) % field
    levels =  input_fields(i_field) % levels
    rows   =  input_fields(i_field) % rows
    CALL fill_external_halos(field, row_length, rows, levels, &
                             halo_x, halo_y)
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE swap_bounds_mv
END MODULE swap_bounds_mv_mod
