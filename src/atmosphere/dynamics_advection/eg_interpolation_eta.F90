! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_interpolation_eta_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_INTERPOLATION_ETA_MOD'

CONTAINS
SUBROUTINE eg_interpolation_eta(                                             &
    eta_in,  pnt_type,                                                       &
    number_of_inputs,                                                        &
    dim_i_in, dim_j_in, dim_k_in,                                            &
    dim_j_in_w,                                                              &
    dim_i_out, dim_j_out, dim_k_out,                                         &
    high_order_scheme, monotone_scheme,                                      &
    l_high, l_mono,                                                          &
    eta_out, lambda_out, phi_out,                                            &
    me, n_proc, n_procx, n_procy,                                            &
    halo_i, halo_j, g_row_length,                                            &
    datastart, at_extremity, g_i_pe,                                         &
    proc_row_group, proc_col_group,                                          &
    halo_data_out_i, halo_data_out_j,                                        &
    error_code,                                                              &
    data_in, data_out,                                                       &
    k_int_linear_in)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE locate_hdps_mod

USE interp_grid_const_mod

USE mpl, ONLY:                                                             &
    mpl_integer,                                                           &
    mpl_real,                                                              &
    mpl_address_kind,                                                      &
    mpl_request_null,                                                      &
    mpl_source,                                                            &
    mpl_any_source,                                                        &
    mpl_status_size

USE missing_data_mod, ONLY:                                                &
    imdi, rmdi

USE um_types,   ONLY: integer32
USE highos_mod, ONLY: cubicLagrange, quinticLagrange,                      &
    ECMWF_quasiCubic, ECMWF_mono_quasiCubic,                               &
    hCubic_vLin, hCubic_vQuintic,                                          &
    hLag3_vHerm3_d2, hLag3_vHerm3_d4

USE monots_mod,  ONLY: mono_quasiCubic

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE
!
! Description:
!
!          Performs interpolation of a field or fields defined on one
!          grid to another grid.
!          Requested output points must
!          lie inside the region defined for the input Data.
!
! Method:
!          For the global model a 'compute-on-demand' method is used
!          for points outside the local grid in the X direction.
!          The 'compute-on-demand' points data is communicated using a
!          Non-blocking Barrier eXchange (NBX) like algorithm from
!          Hoefler et al "Scalable Communication Protocols for Dynamic
!          Sparse Data Exchange". A description of the algorithm (in C) 
!          is also available in the first edition of the Gropp et all 
!          book "Using Advanced MPI" in section 2.1.6 page 29 
!          (ISBN 978-0-262-52763-7).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INTERPOLATION_ETA'


! Arguments with Intent IN. ie: Input variables.

INTEGER ::                                                                 &
    dim_i_in,                                                              &
        ! Dimension of Data_in in i direction.
    dim_j_in,                                                              &
        ! Dimension of Data_in in j direction.
    dim_j_in_w,                                                            &
        ! Dimension of Data_in in j direction.
    dim_k_in,                                                              &
        ! Dimension of Data_in in k direction.
    dim_i_out,                                                             &
        ! Dimension of Data_out in i direction.
    dim_j_out,                                                             &
        ! Dimension of Data_out in j direction.
    dim_k_out,                                                             &
        ! Dimension of Data_out in k direction.
    me,                                                                    &
        ! My processor number
    n_proc,                                                                &
        ! Total number of processors
    n_procx,                                                               &
        ! Number of processors in longitude
    n_procy,                                                               &
        ! Number of processors in latitude
    halo_i,                                                                &
        ! Size of halo in i direction.
    halo_j,                                                                &
        ! Size of halo in j direction.
    halo_data_out_i,                                                       &
        ! size of data out halo in i direction
    halo_data_out_j,                                                       &
        ! size of data out halo in j direction
    proc_row_group,                                                        &
        ! Group id for processors on the same row
    proc_col_group,                                                        &
        ! Group id for processors on the same column
    g_row_length,                                                          &
        ! global number of points on a row
    datastart(3),                                                          &
        ! First gridpoints held by this processor.
    g_i_pe(1-halo_i:g_row_length+halo_i) ! processor on my procr-row
        ! holding a given value in i direction

INTEGER ::                                                                 &
    pnt_type,                                                              &
        ! Defines via an integer code the nature of the
        ! interpolation in terms of which grid the input
        ! Data is on. The codes are given in
        ! terms of a primary variable that would be held
        ! at that point and are u=1, v=2, w=3.
    high_order_scheme,                                                     &
        ! a code saying which high order scheme to use.
    monotone_scheme,                                                       &
        ! a code saying which monotone scheme to use.
    number_of_inputs
        ! number of fields to interpolate.

INTEGER :: k_int_linear_in
        ! Level up to which linear interpolation is used

LOGICAL ::                                                                 &
    l_high,                                                                &
        ! True, if high order interpolation required.
    l_mono
        ! True, if interpolation required to be monotone.


REAL ::                                                                    &
    eta_in(dim_k_in) ! eta coordinate levels.

REAL ::                                                                    &
    data_in  (1-halo_i:dim_i_in+halo_i,                                    &
    1-halo_j:dim_j_in+halo_j, dim_k_in,                                    &
    number_of_inputs )
        ! data to be interpolated

REAL ::                                                                    &
    lambda_out (dim_i_out, dim_j_out, dim_k_out),                          &
        ! Lambda co-ordinate of output data on input.
    phi_out (dim_i_out, dim_j_out, dim_k_out),                             &
        ! Phi Co-ordinate of output data on input.
    eta_out (dim_i_out, dim_j_out, dim_k_out)
        ! Vertical co-ordinate of output data.

LOGICAL ::                                                                 &
    at_extremity(4)
        ! Indicates if this processor is at north,
        ! south, east or west of the processor grid


! Arguments with Intent OUT. ie: Output variables.

REAL ::                                                                    &
    data_out (1-halo_data_out_i:dim_i_out+halo_data_out_i,                 &
    1-halo_data_out_j:dim_j_out+halo_data_out_j,                           &
    dim_k_out, number_of_inputs)
        ! data interpolated to desired locations.

INTEGER ::                                                                 &
    error_code     ! Non-zero on exit if error detected.

! Local Variables.

! scalars

INTEGER :: i, j, k, n, id, jd
INTEGER :: k_int_linear
        ! Level up to which linear interpolation is used
        ! The value of this may be changed, optionally,
        ! by k_int_linear_in

REAL ::   rdel, temp

! arrays

INTEGER (KIND=integer32) ::                                                &
    i_out (dim_i_out, dim_j_out, dim_k_out),                               &
    j_out (dim_i_out, dim_j_out, dim_k_out),                               &
    k_out (dim_i_out, dim_j_out, dim_k_out)


REAL ::                                                                    &
    ext_data(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,          &
    -1:dim_k_in+2,number_of_inputs),                                       &
    data_out_intp (dim_i_out, dim_j_out, dim_k_out,                        &
    number_of_inputs),                                                     &
    weight_lambda (dim_i_out, dim_j_out, dim_k_out),                       &
    weight_phi (dim_i_out, dim_j_out, dim_k_out),                          &
    coeff_z(dim_i_out, dim_j_out, dim_k_out, -2:3)


! Variables applied in the "compute-on-demand" strategy
INTEGER ::                                                                 &
    irecv,                                                                 &
    my_imin,                                                               &
    my_imax,                                                               &
    dim_e_out,                                                             &
    h_factor,                                                              &
    nrecv,                                                                 &
    info,                                                                  &
    itmp,                                                                  &
    kk,                                                                    &
    dsm1,                                                                  &
    itmpt


! Variables needed for the NBX like algorithm
INTEGER, PARAMETER :: tag_in   = 11
INTEGER, PARAMETER :: tag_out  = 41
INTEGER, PARAMETER :: tag_loop = 2
INTEGER, SAVE      :: tag_call = 0


INTEGER :: barrier_req
INTEGER :: reqs_send_in (n_procx)
INTEGER :: reqs_send_out(n_procx)
INTEGER :: reqs_recv_out(n_procx)
INTEGER :: reqs(n_procx*2)

INTEGER :: idx
INTEGER :: pe 
INTEGER :: num_reqs_send_in
INTEGER :: num_reqs_send_out
INTEGER :: num_reqs_recv_out
INTEGER :: num_reqs


INTEGER :: status(mpl_status_size)
INTEGER :: statuses(mpl_status_size,n_procx)

LOGICAL :: barrier_done
LOGICAL :: barrier_active
LOGICAL :: flag

INTEGER :: n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)
INTEGER :: n_recv_out

INTEGER (KIND=integer32), ALLOCATABLE ::                                   &
    i_out_e(:),j_out_e(:),k_out_e(:)

! The following arrays are defined as 3/4D but treated locally as 1D
REAL, ALLOCATABLE ::                                                       &
    weight_lambda_e(:,:,:),                                                &
    weight_phi_e(:,:,:),                                                   &
    coeff_z_e(:,:,:,:),                                                    &
    eta_out_e(:,:,:)

! Variables used for index vectors

INTEGER ::                                                                 &
    ic,                                                                    &
    itmpa(dim_i_out, dim_j_out, dim_k_out),                                &
    icnt(dim_j_out, dim_k_out),                                            &
    inx(dim_i_out, dim_j_out, dim_k_out)

! pointer for stretched grid data
REAL, POINTER :: s_xi1(:), t_xi1(:), s_xi2(:), t_xi2(:),                   &
    q_xi1(:,:), q_xi2(:,:)

! Variables needed to create the MPI derived data type which is used
! to send and receive the input data of the 'compute-on-demand' points
INTEGER (KIND=mpl_address_kind) :: extent
INTEGER (KIND=mpl_address_kind) :: lb
INTEGER, SAVE :: mpl_send_type = imdi
INTEGER :: oldtypes(0:1)
INTEGER :: blockcounts(0:1)
INTEGER (KIND=mpl_address_kind) :: offsets(0:1)

TYPE sendrecv_type
  SEQUENCE
  INTEGER  :: i_out
  INTEGER  :: j_out
  INTEGER  :: k_out
  REAL     :: weight_lambda
  REAL     :: weight_phi
  REAL     :: eta_out
END TYPE sendrecv_type

! Arrays used to store the indices of the
! 'compute-on-demands' points in the sender
TYPE indices_type
  INTEGER (KIND=integer32), ALLOCATABLE :: i(:)
  INTEGER (KIND=integer32), ALLOCATABLE :: j(:)
  INTEGER (KIND=integer32), ALLOCATABLE :: k(:)
END TYPE indices_type

TYPE (indices_type) :: store(0:n_procx-1)

! Arrays used to store the input data of the
! 'compute-on-demands' points

TYPE  data_in_type
  TYPE (sendrecv_type), ALLOCATABLE :: points(:)  
END TYPE data_in_type

TYPE (data_in_type) :: recv_in(0:n_procx-1)
TYPE (data_in_type) :: send_in(0:n_procx-1)  

! Arrays used to store the output data of the 'compute-on-demand'
! points, i.e. the interpolation results

TYPE  data_out_type
  REAL, ALLOCATABLE :: intp_e(:,:,:,:)
END TYPE data_out_type

TYPE (data_out_type) :: recv_out(0:n_procx-1)  
TYPE (data_out_type) :: send_out(0:n_procx-1)  

! Extra variables

INTEGER :: d_imin,d_imax, lev_ext
INTEGER :: k_found

LOGICAL                     :: l_cubic_interp


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Check if linear interpolation is to be used on more than one layer
! near the lower boundary.

 k_int_linear=k_int_linear_in

! ---------------------------------------------------------------------
! Section 2.   Call appropriate routine to extend data.
!              Minimum amount of extending done to cope with highest
!              order interpolation scheme requested. This can leave
!              unset values in Ext_Data and Ext_r_in, but
!              these values will not be used.

!              Parallel version: No extension in r required.
! ---------------------------------------------------------------------

IF ( high_order_scheme  ==  ecmwf_quasicubic      .OR.                      &
    high_order_scheme  ==  ecmwf_mono_quasicubic .OR.                      &
    monotone_scheme    ==  mono_quasiCubic           ) THEN
  
  error_code = 1
  CALL ereport(RoutineName, error_code,                                     &
      "Interpolation options not available within EG" )
END IF

! This is used for the comms-on-demand and wind limiter
! - basically the size of the interpolation stencil

h_factor = 1    ! linear
IF (high_order_scheme  ==  quinticlagrange ) h_factor = 3

l_cubic_interp = .FALSE.
IF ( high_order_scheme == cubiclagrange     .OR.                            &
    high_order_scheme == hcubic_vquintic   .OR.                            &
    high_order_scheme == hcubic_vlin       .OR.                            &
    high_order_scheme == hLag3_vHerm3_d2   .OR.                            &
    high_order_scheme == hLag3_vHerm3_d4        ) THEN

  h_factor       = 2
  l_cubic_interp = .TRUE.
END IF

lev_ext = 0
IF ( high_order_scheme  ==  quinticlagrange .OR.                           &
    high_order_scheme  ==  hcubic_vquintic .OR.                            &
    high_order_scheme  ==  hLag3_vHerm3_d4      ) THEN

  lev_ext = 2

ELSE IF ( high_order_scheme  ==  cubiclagrange .OR.                        &
    high_order_scheme  ==  hLag3_vHerm3_d2    ) THEN

  lev_ext = 1

END IF

! Copy core data

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( number_of_inputs, dim_k_in, halo_j, halo_i, dim_i_in,    &
!$OMP         dim_j_in, Ext_Data, Data_in, lev_ext )                   &
!$OMP PRIVATE( i, j, k, n )
DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, dim_k_in
    DO j = 1-halo_j, dim_j_in+halo_j
      DO i = 1-halo_i, dim_i_in+halo_i
        Ext_Data(i,j,k,n) = Data_in(i,j,k,n)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END DO

! Extend bottom

DO n = 1, number_of_inputs
  DO k = 1-lev_ext, 0
!$OMP DO SCHEDULE(STATIC)
    DO j = 1-halo_j, dim_j_in+halo_j
      DO i = 1-halo_i, dim_i_in+halo_i
        Ext_Data(i,j,k,n) = rmdi
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  ! Extend Top

  DO k = 1+dim_k_in, dim_k_in + lev_ext
!$OMP DO SCHEDULE(STATIC)
    DO j = 1-halo_j, dim_j_in+halo_j
      DO i = 1-halo_i, dim_i_in+halo_i
        Ext_Data (i,j,k,n) = rmdi
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END DO
!$OMP END PARALLEL

! ---------------------------------------------------------------------
! Section 3.  For each output point find i,j,k so that the point on the
!             output grid lies between i and i+1, j and j+1, k and k+1
! ---------------------------------------------------------------------

SELECT CASE(pnt_type)
CASE (fld_type_u)
  s_xi1 => sig_xi1_u
  t_xi1 => tau_xi1_u
  s_xi2 => sig_xi2_p
  t_xi2 => tau_xi2_p
  q_xi1 => q1_u
  q_xi2 => q2_p
CASE (fld_type_v)
  s_xi1 => sig_xi1_p
  t_xi1 => tau_xi1_p
  s_xi2 => sig_xi2_v
  t_xi2 => tau_xi2_v
  q_xi1 => q1_p
  q_xi2 => q2_v
CASE DEFAULT
  s_xi1 => sig_xi1_p
  t_xi1 => tau_xi1_p
  s_xi2 => sig_xi2_p
  t_xi2 => tau_xi2_p
  q_xi1 => q1_p
  q_xi2 => q2_p
END SELECT

dsm1 = datastart(2)-1

! Find i and j point.

CALL locate_hdps(i_out, j_out, weight_lambda, weight_phi,                  &
    lambda_out, phi_out, dsm1, datastart(1), h_factor,                     &
    pnt_type, dim_i_out, dim_j_out, dim_k_out, halo_i, halo_j)

! New search algorithm using the continuous nature of the search
! space to optimise

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(static)                              &
!$OMP     PRIVATE(k,j,i,kk, k_found) SHARED(dim_k_out, dim_j_out,              &
!$OMP     dim_i_out, dim_k_in, k_out, eta_out, eta_in)
DO k = 1, dim_k_out
  DO j = 1, dim_j_out
    DO i = 1, dim_i_out

      ! Start searching from the input data point
      k_found = k
      IF (eta_out(i,j,k) <= eta_in(k_found)) THEN

        !search down
        DO kk = k_found, 1, -1
          k_found = kk
          IF (eta_out(i,j,k) > eta_in(kk)) EXIT
        END DO

      ELSE

        ! search up
        IF (k_found == dim_k_in) THEN
          k_found = k_found-1
        ELSE

          DO kk = k_found+1, dim_k_in-1
            IF (eta_out(i,j,k) < eta_in(kk)) THEN
              k_found = kk-1
              EXIT
            END IF
          END DO
        END IF
      END IF

      k_out(i,j,k) = k_found

    END DO
  END DO
END DO
!$OMP END PARALLEL DO


IF ( n_procx > 1 ) THEN
  IF (model_type == mt_global) THEN

    ! Send the points outside my region to the appropriate processor for
    ! interpolation. Only performed if the domain is decomposed in the
    ! i direction.

    ! The first and last point I can interpolate in, based on available
    ! data on this processor (minus/plus one to avoid use of ge/le)

    my_imin = datastart(1) - halo_i + h_factor - 1
    my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1

    DO pe = 0, n_procx-1
      n_sendto(pe)   = 0
      n_recvfrom(pe) = 0
    END DO

    ! And now decide where a point should be evaluated

    d_imin=my_imin-datastart(1) + 1
    d_imax=my_imax-datastart(1) + 1

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC)                          &
!$OMP SHARED( dim_k_out, dim_j_out, dim_i_out, i_out, datastart, d_imin,   &
!$OMP         d_imax, itmpa, icnt, inx, g_row_length, halo_i, h_factor,    &
!$OMP         n_sendto, g_i_pe)                                            &
!$OMP PRIVATE( i, j, k, ic, irecv )
    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        END DO
        icnt(j,k) = 0

        DO i = 1, dim_i_out
          IF (i_out(i,j,k)  <=  d_imin .OR.                                &
              i_out(i,j,k)  >=  d_imax) THEN
            icnt(j,k) = icnt(j,k)+1
            inx(icnt(j,k),j,k) = i
          END IF
        END DO
        DO ic=1, icnt(j,k)
          i = inx(ic,j,k)
          i_out(i,j,k) = i_out(i,j,k) + datastart(1) - 1

          ! Send to a remote processor, given by the array g_i_pe

          !   CODE TO STOP BIT NON-REPRODUCIBILITY
          IF (i_out(i,j,k) >  g_row_length+halo_i-h_factor) THEN
            i_out(i,j,k) = i_out(i,j,k)-g_row_length
          END IF
          IF (i_out(i,j,k) <  1-halo_i+h_factor) THEN
            i_out(i,j,k) = i_out(i,j,k)+g_row_length
          END IF
          !   END CODE TO STOP BIT NON-REPRODUCIBILITY
        END DO
      END DO      
    END DO
!$OMP END PARALLEL DO

    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        ! Process locally, so find the local destination
        DO ic=1, icnt(j,k)
          i = inx(ic,j,k)
          irecv = g_i_pe(i_out(i,j,k))
          n_sendto(irecv) = n_sendto(irecv) + 1
          itmpa(ic,j,k)  = n_sendto(irecv)
        END DO
      END DO
    END DO

    ! Allocate the array to send the input data and the arrays with
    ! the indices needed to send and receive the 'compute-on-demand'
    ! points
    DO pe = 0, n_procx-1
      IF (n_sendto(pe)  >   0) THEN
        ALLOCATE( send_in(pe) % points( n_sendto(pe) ) )
        ALLOCATE( store  (pe) % i     ( n_sendto(pe) ) )
        ALLOCATE( store  (pe) % j     ( n_sendto(pe) ) )
        ALLOCATE( store  (pe) % k     ( n_sendto(pe) ) )
      END IF
    END DO

    ! Prepare the input data
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE)                          &
!$OMP SHARED( dim_k_out, dim_j_out, icnt, inx, g_i_pe, i_out, j_out, k_out,&
!$OMP         send_in, weight_lambda, weight_phi, eta_out, store, itmpa)   &
!$OMP PRIVATE( i, j, k, ic, irecv, itmp )
    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        ! Process locally, so find the local destination
        DO ic=1, icnt(j,k)
          i = inx(ic,j,k)
          irecv = g_i_pe(i_out(i,j,k))
          itmp = itmpa(ic,j,k)

          send_in(irecv) % points (itmp) % i_out         = i_out(i,j,k)
          send_in(irecv) % points (itmp) % j_out         = j_out(i,j,k)
          send_in(irecv) % points (itmp) % k_out         = k_out(i,j,k)
          send_in(irecv) % points (itmp) % weight_lambda = weight_lambda(i,j,k)
          send_in(irecv) % points (itmp) % weight_phi    = weight_phi(i,j,k)
          send_in(irecv) % points (itmp) % eta_out       = eta_out(i,j,k)
          store(irecv) % i(itmp) = i
          store(irecv) % j(itmp) = j
          store(irecv) % k(itmp) = k
          i_out(i,j,k) = i

        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! Get types setup if not done
    IF (mpl_send_type == imdi) THEN
      offsets    (0) = 0
      oldtypes   (0) = mpl_integer
      blockcounts(0) = 3

      CALL mpl_type_get_extent(mpl_integer, lb, extent, info)

      offsets    (1) = 3 * extent
      oldtypes   (1) = mpl_real
      blockcounts(1) = 3

      CALL mpl_type_create_struct(2, blockcounts, offsets, oldtypes,       &
          mpl_send_type, info)
      CALL mpl_type_commit(mpl_send_type, info)
    END IF

    ! Initialise NBX protocol

    barrier_req      = mpl_request_null
    reqs_send_in(:)  = mpl_request_null
    reqs_send_out(:) = mpl_request_null
    reqs_recv_out(:) = mpl_request_null

    barrier_done   = .FALSE.
    barrier_active = .FALSE.

    ! Each call of this routine should use unique tags values that
    ! cycle between a range.
    !
    ! Rationale:    
    !
    ! There is a possible race condition, which was observed in XC40,
    ! when the routine (and thus the NBX protocol) uses the same tag
    ! for each call and the routine is called twice in
    ! succession. There is an unknown situation when all the PEs have
    ! activated the non-blocking barrier but one/some PEs do an extra
    ! iprobe and start to receive the data from the PEs that are
    ! already in the second call of the routine. The use of different
    ! tags for each call stop this race condition.
    ! 
    tag_call = tag_call + 1
    tag_call = MOD(tag_call,tag_loop)

    ! Send the input data using a non-blocking synchronised send,
    ! i.e. the request is completed only when the receiver has started
    ! executing the matching receive and the send buffer can be safely
    ! re-used. 

    num_reqs_send_in  = 0

    DO pe = 0,n_procx-1
      IF (n_sendto(pe)  >   0) THEN

        num_reqs_send_in  = num_reqs_send_in  + 1
                        
        CALL mpl_issend(send_in(pe) % points, n_sendto(pe),                  &
          mpl_send_type, pe, tag_in + tag_call, proc_row_group,              &
          reqs_send_in(num_reqs_send_in), info )

      END IF
    END DO

    ! Loop until all neighbour processes pass the non-blocking
    ! barrier, i.e. have finished sending and receiving the input data.
    DO WHILE(.NOT. barrier_done)

      CALL mpl_iprobe(mpl_any_source, tag_in + tag_call, proc_row_group,     &
        flag, status, info)

      ! If a message with the input data is ready, receive and process
      ! the input data.
      IF (flag) THEN

        ! Get the process sending the data.
        pe = status(mpl_source)

        ! Get the number of elements in the data.
        CALL mpl_get_count(status, mpl_send_type, n_recvfrom(pe), info)

        ! Allocate the memory needed and receive the input data.
        ALLOCATE( recv_in(pe) % points (n_recvfrom(pe)) )

        CALL mpl_recv(recv_in(pe) % points, n_recvfrom(pe), mpl_send_type,  &
          pe, tag_in + tag_call, proc_row_group, status, info )

      END IF

      IF (.NOT. barrier_active) THEN

        ! This process is still waiting for its input data to be
        ! received by the neighbour processes. Check if it is still
        ! the case.
        CALL mpl_testall (num_reqs_send_in, reqs_send_in, flag, statuses, info)
        
        IF (flag) THEN
          
          ! The non-blocking synchronised send requests are completed,
          ! i.e. the neighbour processes started to execute the
          ! matching receive.  Pass the barrier, i.e. signal that
          ! all the input data was sent.
          CALL mpl_ibarrier(proc_row_group, barrier_req, info)
          barrier_active = .TRUE.

        END IF

      ELSE 

        ! This process is waiting for all the processes to pass the
        ! non-blocking barrier
        CALL mpl_test(barrier_req, barrier_done, status, info)

      END IF

    END DO ! Barrier done and all input data sent and received   
    
  ELSE  !  model_type is NOT GLOBAL

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( dim_k_out, dim_j_out, dim_i_out, i_out, datastart )      &
!$OMP PRIVATE( i, j, k )
    DO k = 1, dim_k_out
      DO j = 1, dim_j_out
        DO i = 1, dim_i_out
          i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO


  END IF ! model_type == mt_global
END IF ! n_procx > 1

! CODE TO STOP BIT NON-REPRODUCIBILITY
IF (model_type == mt_global .AND. n_procx == 1) THEN
  h_factor = 2
  IF (high_order_scheme  ==  quinticlagrange ) THEN
    h_factor = 3
  END IF
  my_imin = datastart(1) - halo_i + h_factor - 1
  my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
  DO k = 1, dim_k_out
    DO j = 1,dim_j_out
      DO i = 1, dim_i_out
        IF (i_out(i,j,k) >= my_imax) THEN
          i_out(i,j,k)=i_out(i,j,k) - g_row_length
        END IF
        IF (i_out(i,j,k) <= my_imin ) THEN
          i_out(i,j,k)=i_out(i,j,k)+g_row_length
        END IF
      END DO
    END DO
  END DO
END IF   ! model_type == mt_global .and. n_procx == 1
!  END CODE TO STOP BIT NON-REPRODUCIBILITY

! ----------------------------------------------------------------------
! Section 4.   Perform required Interpolations.
! ----------------------------------------------------------------------

! Compute interpolation procedure for the "compute-on-demand" points

IF (n_procx > 1 .AND. model_type == mt_global) THEN

  ! Allocate the array to receive and send the output data
  DO pe = 0, n_procx-1
    IF (n_sendto(pe)  >   0) THEN
      ALLOCATE( recv_out(pe) % intp_e( 1, 1, n_sendto(pe),   number_of_inputs ))
    END IF
    IF (n_recvfrom(pe)  >   0) THEN
      ALLOCATE( send_out(pe) % intp_e( 1, 1, n_recvfrom(pe), number_of_inputs ))
    END IF
  END DO

  num_reqs_recv_out = 0
  num_reqs_send_out = 0

  ! Post the command to receive the output data

  DO pe = 0,n_procx-1
    IF (n_sendto(pe)  >   0) THEN
      
      num_reqs_recv_out = num_reqs_recv_out + 1

      CALL mpl_irecv(recv_out(pe) % intp_e, n_sendto(pe)*number_of_inputs, & 
        mpl_real, pe, tag_out + tag_call, proc_row_group,                  &
        reqs_recv_out(num_reqs_recv_out), info )
          
    END IF
  END DO

  ! Interpolate points

  DO pe = 0, n_procx-1

    IF (n_recvfrom(pe)  >   0) THEN
      
      dim_e_out = n_recvfrom(pe)

      ALLOCATE( i_out_e (dim_e_out) )
      ALLOCATE( j_out_e (dim_e_out) )
      ALLOCATE( k_out_e (dim_e_out) )
      ALLOCATE( weight_lambda_e(1,1,dim_e_out) )
      ALLOCATE( weight_phi_e(1,1,dim_e_out) )
      ALLOCATE( eta_out_e(1,1,dim_e_out) )
      ALLOCATE( coeff_z_e(1,1,dim_e_out,-2:3) )
  
      ! Prepare the "compute-on-demand" points data
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( dim_e_out, i_out_e, j_out_e, k_out_e, recv_in, pe,       &
!$OMP         weight_lambda_e, weight_phi_e, eta_out_e, datastart)     &
!$OMP PRIVATE( j )
      DO j = 1, dim_e_out
        i_out_e(j) = recv_in(pe) % points(j) % i_out - datastart(1) + 1
        j_out_e(j) = recv_in(pe) % points(j) % j_out
        k_out_e(j) = recv_in(pe) % points(j) % k_out
        weight_lambda_e(1,1,j) = recv_in(pe) % points(j) % weight_lambda
        weight_phi_e(1,1,j)    = recv_in(pe) % points(j) % weight_phi
        eta_out_e(1,1,j)       = recv_in(pe) % points(j) % eta_out
      END DO
!$OMP END PARALLEL DO

      ! DEPENDS ON: eg_vert_weights_eta
      CALL eg_vert_weights_eta (                                           &
        dim_k_in, 1, 1, dim_e_out,                                         &
        high_order_scheme, monotone_scheme,                                &
        l_high, k_int_linear,                                              &
        k_out_e, eta_in, eta_out_e,                                        &
        coeff_z_e)

      ! Call high order scheme if required.
      
      IF (l_cubic_interp) THEN
        ! DEPENDS ON: eg_cubic_lagrange
        CALL eg_cubic_lagrange(                                            &
          ext_data,                                                        &
          dim_i_in, dim_j_in, dim_k_in,                                    &
          1, 1, dim_e_out,                                                 &
          halo_i, halo_j, number_of_inputs,                                &
          weight_lambda_e,weight_phi_e,                                    &
          s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,                        &
          i_out_e, j_out_e, k_out_e,                                       &
          coeff_z_e, lev_ext,                                              &
          send_out(pe) % intp_e)
        
      ELSE IF (high_order_scheme  ==  quinticlagrange) THEN
        
        ! DEPENDS ON: eg_quintic_lagrange
        CALL eg_quintic_lagrange(                                          &
          ext_data,                                                        &
          dim_i_in, dim_j_in, dim_k_in,                                    &
          1, 1, dim_e_out,                                                 &
          halo_i, halo_j, number_of_inputs,                                &
          weight_lambda_e,weight_phi_e,                                    &
          i_out_e, j_out_e, k_out_e,                                       &
          coeff_z_e,                                                       &
          send_out(pe) % intp_e)
        
      END IF
      
      IF ( l_mono ) THEN
        IF ( l_high ) THEN
          
          ! DEPENDS ON: mono_enforce
          CALL mono_enforce(                                               &
            ext_data, number_of_inputs,                                    &
            dim_i_in, dim_j_in, dim_k_in,                                  &
            1, 1, dim_e_out,                                               &
            halo_i, halo_j,                                                &
            i_out_e, j_out_e, k_out_e,                                     &
            send_out(pe) % intp_e)
        ELSE
          
          ! DEPENDS ON: eg_tri_linear
          CALL eg_tri_linear(                                              &
            ext_data,                                                      &
            dim_i_in, dim_j_in, dim_k_in,                                  &
            1, 1, dim_e_out,                                               &
            halo_i, halo_j, number_of_inputs,                              &
            weight_lambda_e,weight_phi_e,                                  &
            i_out_e, j_out_e, k_out_e,                                     &
            coeff_z_e, send_out(pe) % intp_e)
        END IF
      END IF
      
      num_reqs_send_out  = num_reqs_send_out  + 1

      ! Send the output data of the 'compute-on-demand' points back
      ! to the original process
      CALL mpl_isend(send_out(pe) % intp_e, dim_e_out*number_of_inputs,    &
        mpl_real, pe, tag_out + tag_call, proc_row_group,                  &
        reqs_send_out(num_reqs_send_out), info )
            
      ! Deallocate array
      DEALLOCATE( coeff_z_e )
      DEALLOCATE( eta_out_e )
      DEALLOCATE( weight_phi_e )
      DEALLOCATE( weight_lambda_e )
      DEALLOCATE( k_out_e )
      DEALLOCATE( j_out_e )
      DEALLOCATE( i_out_e )

    END IF 
    
  END DO 

END IF 


! Compute interpolation procedure for the "local" points

! DEPENDS ON: eg_vert_weights_eta
CALL eg_vert_weights_eta(                                                  &
    dim_k_in, dim_i_out, dim_j_out, dim_k_out,                             &
    high_order_scheme, monotone_scheme,                                    &
    l_high, k_int_linear,                                                  &
    k_out, eta_in, eta_out,                                                &
    coeff_z )

! Call high order scheme if required.

IF (l_cubic_interp) THEN
  ! DEPENDS ON: eg_cubic_lagrange
  CALL eg_cubic_lagrange(                                                  &
      ext_data,                                                            &
      dim_i_in, dim_j_in, dim_k_in,                                        &
      dim_i_out, dim_j_out, dim_k_out,                                     &
      halo_i, halo_j, number_of_inputs,                                    &
      weight_lambda, weight_phi,                                           &
      s_xi1, s_xi2, t_xi1, t_xi2, q_xi1, q_xi2,                            &
      i_out, j_out, k_out,                                                 &
      coeff_z, lev_ext,                                                    &
      data_out_intp)

ELSE IF (high_order_scheme  ==  quinticlagrange) THEN

  ! DEPENDS ON: eg_quintic_lagrange
  CALL eg_quintic_lagrange(                                                &
      ext_data,                                                            &
      dim_i_in, dim_j_in, dim_k_in,                                        &
      dim_i_out, dim_j_out, dim_k_out,                                     &
      halo_i, halo_j, number_of_inputs,                                    &
      weight_lambda, weight_phi,                                           &
      i_out, j_out, k_out,                                                 &
      coeff_z,                                                             &
      data_out_intp )
END IF

IF ( l_mono ) THEN
  IF (l_high ) THEN
    ! DEPENDS ON: mono_enforce
    CALL mono_enforce(                                                     &
        ext_data, number_of_inputs,                                        &
        dim_i_in, dim_j_in, dim_k_in,                                      &
        dim_i_out, dim_j_out, dim_k_out,                                   &
        halo_i, halo_j,                                                    &
        i_out, j_out, k_out,                                               &
        data_out_intp)
  ELSE
    ! DEPENDS ON: eg_tri_linear
    CALL eg_tri_linear(                                                    &
        ext_data,                                                          &
        dim_i_in, dim_j_in, dim_k_in,                                      &
        dim_i_out, dim_j_out, dim_k_out,                                   &
        halo_i, halo_j, number_of_inputs,                                  &
        weight_lambda, weight_phi,                                         &
        i_out, j_out, k_out,                                               &
        coeff_z,data_out_intp )
  END IF
END IF


IF (n_procx > 1 .AND. model_type == mt_global) THEN

  ! Receive the data output of the "compute-on-demand" points

  ! Merge and wait the requests for the receiving and sending of
  ! output data.
  num_reqs = 0
  reqs(:)  = mpl_request_null

  num_reqs = num_reqs_recv_out + num_reqs_send_out
  DO i = 1, num_reqs_recv_out
    reqs(i) = reqs_recv_out(i)
  END DO
  DO i = 1, num_reqs_send_out
    reqs(i+num_reqs_recv_out) = reqs_send_out(i)
  END DO

  IF (num_reqs > 0) THEN 
    
    ! At the end of this loop, all send and receive requests must have
    ! completed, since the waitany will have been called the requisite
    ! number of times
    DO i = 1, num_reqs

      CALL mpl_waitany(num_reqs, reqs, idx, status, info)
    
      ! Continue loop if it was a "send" that completed
      IF (idx>num_reqs_recv_out) THEN
        CYCLE
      END IF

      ! It was a "receive" that completed

      ! Get the the process sending the data.
      pe = status(mpl_source)
                
      ! Get the number of output elements received.
      CALL mpl_get_count(status, mpl_real, n_recv_out, info)

      ! Check that the number of output element is equal to the number
      ! of input elements (that was communicated to the pe) multiplied
      ! by the number of inputs. If this check fails something wrong
      ! happen during the communication
      IF(n_recv_out /= n_sendto(pe)*number_of_inputs) THEN
        error_code = 1
        CALL ereport(RoutineName, error_code,                           &
          "Wrong number of output elements received" )
      END IF
      
      ! Copy the output data 
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( pe, number_of_inputs, n_sendto, data_out_intp, store,    &
!$OMP         recv_out )                                               &
!$OMP PRIVATE( n, j )

      DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC)
        DO j = 1,n_sendto(pe)
          data_out_intp(store(pe)%i(j), store(pe)%j(j), store(pe)%k(j), n) = &
            recv_out(pe) % intp_e(1,1,j,n)
        END DO
!$OMP END DO NOWAIT
      END DO

!$OMP END PARALLEL

    END DO
  END IF

  ! Deallocate the array used to send and receive the input data
  DO pe = 0, n_procx-1
    IF ( ALLOCATED( send_in(pe) % points ) ) DEALLOCATE (send_in(pe) % points)
    IF ( ALLOCATED( recv_in(pe) % points ) ) DEALLOCATE (recv_in(pe) % points)
    IF ( ALLOCATED( store  (pe) % i      ) ) DEALLOCATE (store  (pe) % i)
    IF ( ALLOCATED( store  (pe) % j      ) ) DEALLOCATE (store  (pe) % j)
    IF ( ALLOCATED( store  (pe) % k      ) ) DEALLOCATE (store  (pe) % k)
  END DO
  
  ! Deallocate the array used to send and receive the output data
  DO pe = 0, n_procx-1
    IF ( ALLOCATED( recv_out(pe) % intp_e ) ) DEALLOCATE (recv_out(pe) % intp_e)
    IF ( ALLOCATED( send_out(pe) % intp_e ) ) DEALLOCATE (send_out(pe) % intp_e)
  END DO
  
END IF 


! ----------------------------------------------------------------------
! Section 5.  Put interpolated field into Data_out.
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                         &
!$OMP& SHARED( number_of_inputs, dim_k_out, dim_j_out, dim_i_out,         &
!$OMP&         data_out, data_out_intp )                                  &
!$OMP& PRIVATE( i, j, k, n ) COLLAPSE(2)
DO n = 1, number_of_inputs
  DO k = 1, dim_k_out
    DO j = 1, dim_j_out
      DO i = 1, dim_i_out
        data_out (i,j,k,n) = data_out_intp (i,j,k,n)
      END DO
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_interpolation_eta
END MODULE eg_interpolation_eta_mod
