! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

!     SUBROUTINE RFCSL-1A
!     RECURSIVE FILTER - CONSTANT COEFF - FOR A SCALAR LTD.AREA FIELD
!
!     THIS VERSION FOR A VARIABLE RESOLUTION GRID AND DLAT AND DLONG
!     ARE VECTORS OF THE GRID SPACING
!
!     VERSION 1A: 
!     This version works by dividing the global domain into horizontal
!     and vertical stripes. Then the field data is continously moved
!     from the logical processed grid (LPG) decomposition to the
!     horizontal stripes decomposition when the E/W filter is
!     performed and to the vertical stripes decomposition when the N/S
!     filter is performed. The filter is applied to the global field
!     at multiple levels by continously moving data between domains,
!     i.e. global communication between PEs.
!
!
!     ARGUMENTS:
MODULE rfcsl_1a_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RFCSL_1A_MOD'

CONTAINS

SUBROUTINE rfcsl_1a(f_l,rows_l,row_length_l,model_levels,             &
  rows_g, row_length_g, m_grid, f_boundary ,                          &
  coslat_g, dlat_g, dlong_g, fi_scale, npass)                            
  

USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc, nproc_max
USE field_types, ONLY: fld_type_p
USE mpl, ONLY: mpl_real
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!---- INPUT
INTEGER, INTENT(IN) :: rows_l             ! Rows in local grid
INTEGER, INTENT(IN) :: row_length_l       ! Points in each row of local grid
INTEGER, INTENT(IN) :: model_levels       ! Levels
INTEGER, INTENT(IN) :: rows_g             ! Rows in global grid
INTEGER, INTENT(IN) :: row_length_g       ! Points in each row of global grid
INTEGER, INTENT(IN) :: m_grid             ! 0=use f_boundary, -1=no  boundary
REAL, INTENT(IN) :: f_boundary            ! f assumed outside ltd.area
REAL, INTENT(IN) :: coslat_g(rows_g)      ! Cosine(latitude)
REAL, INTENT(IN) :: dlat_g(rows_g)        ! Row spacing (radians)
REAL, INTENT(IN) :: dlong_g(row_length_g) ! Point spacing (radians)
REAL, INTENT(IN) :: fi_scale              ! Filter scale (radians)
INTEGER, INTENT(IN) :: npass              ! Number of passes in each direction
                                 ! npass=2  equivalent of soar function
                                 ! npass=large equivalent of gaussian function

REAL, INTENT(INOUT) :: f_l(row_length_l,rows_l,model_levels)     
! Local field to be filtered

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RFCSL_1A'

!-----------------------------------------------------------------------
! **        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS F
!     FIRST PASS IN ONE DIRECTION GIVES G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES H (OVERWRITING ARRAY F):
!      H(I)=A(I)*H(I+1)+B(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!     FOR FILTER OF LIMITED AREA THE GRID IS
!      ASSUMED REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!      BECAUSE OF THIS APPROX, THE GRID SHOULD NOT GO NEAR POLE.
!-----------------------------------------------------------------------
!     NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------

! The global domain is divided in nproc_x X nproc_y PEs. The field array
! passed to this routine is a local field, i.e. the data only of this PE.
! 
! Since the filtering is recursive a full global row is needed for a
! E/W pass of the filter and a full global column is needed for a N/S
! pass of the filter. 
!
! Thus in order to parallelise the computation, two new domains are
! created: 
! 1) multiple horizontal stripes domain where the the E/W pass is performed.
! 2) multiple vertical   stripes domain  where the the N/S pass is performed.
! 
! The data is communicated between the various domains using MPI alltoall.
!
! Example:
!
! * Normal logical processed grid (LPG) decomposition 
!
! LPG decomp   2x3
! global_row_length = 9
! global_rows       = 9

!      4444      55555
! P4   4444      55555   P5
!      4444      55555

!      2222      33333
! P2   2222      33333   P3
!      2222      33333

!      0000      11111
! P0   0000      11111   P1
!      0000      11111
!
! * Horizontal stripes - substripes
!
! The horizonantal stripes decomposition is done so that when data is
! communicated from LPG decomp to Horiz stripe decomp, each PE needs only
! to commmunicate with the PEs on the same LPG row, gc_proc_row_group.
!
! The number of horizontal stripes is the same as nproc_y (3)
! The number of horizontal substripes in a stripe is the same as
! nproc_x (2)
!
! P5   444455555
!                     Stripe 1
!      444455555       
! P4   444455555
!
! -----------------
!
! P3   222233333
!                     Stripe 2
!      222233333      
! P2   222233333
!
! -----------------
!
! P1   000011111
!                     Stripe 3
!      000011111
! P0   000011111
!
!
! * Vertical stripes - substripes
!
! The vertical stripes decomposition is done so that when data is
! communicated from vertical stripes decomp to LPG decomp, each PE needs only
! to commmunicate with the PEs on the same LPG column, gc_proc_col_group.
!
! The number of vertical stripes is the same as nproc_x (2)
! The number of vertical substripes in a stripe is the same as
! nproc_y (3)
!
!   Stripe 1       Stripe 2
!              | 
!  P0  P2  P4  |  P1  P3  P5
!
!  44  4   4   |  55  55  5
!  44  4   4   |  55  55  5
!  44  4   4   |  55  55  5
!  22  2   2   |  33  33  3
!  22  2   2   |  33  33  3
!  22  2   2   |  33  33  3
!  00  0   0   |  11  11  1
!  00  0   0   |  11  11  1
!  00  0   0   |  11  11  1
!              |
!
! When data needs to be communicated from Horizontal stripes decomp to
! Vertical stripes decomp, and vice-versa, each PE needs to commuicate
! wirh all the PEs.
!
!
! The variables used in this routine can have different subscripts:
! '_g_'  for global variables
! '_l_'  for local LPG variables 
! '_h_'  for horizontal stripe/sub-stripe variables
! '_v_'  for vertical   stripe/sub-stripe variables

! Work space
REAL :: a_ns(rows_g)               ! Filter coefficients N-S-N
REAL :: a_we(row_length_g,rows_g)  ! Filter coefficients W-E-W
REAL :: e,z                        ! Intermediate values for filter coefficient

INTEGER :: i,j,k, dstart
INTEGER :: pe, pe_i, pe_j, prev_pe  
INTEGER :: ierr
INTEGER :: jpass

LOGICAL, SAVE :: first_call = .TRUE.

REAL, ALLOCATABLE ::                                                    &
    f_h  (:,:,:)                                                        &
             ! Horizontal sub-stripe of the field on which to perform the E/W
             ! pass of the filter

,   f_v  (:,:,:)                                               
             ! Vertical sub-stripe of the field on which to perform the N/S
             ! pass of the filter

! Local variables for sub-stripe decompositions 
INTEGER, SAVE ::                                                        & 
     rows_h                                                             &
             ! The rows in the horizontal sub-stripe for the PE
,    row_length_v                                                       &
             ! The row length in the vertical sub-stripe for the PE
,    max_rows_h                                                         &
             ! Maximum rows amongst all the horizontal sub-stripes
,    max_row_length_v                                                   &
             ! Maximum row length amongst all the vertical sub-stripes
,    max_rows_l                                                         &
             ! Maximum rows amongst all the LPG  
,    max_row_length_l                                                   
             ! Maximum row length amongst all the lpg

INTEGER, ALLOCATABLE, SAVE ::                                           & 
     pe_rows_h(:)                                                       &
             ! The rows in the horizontal sub-stripe of all PEs
,    pe_row_length_v(:)                                                 &     
             ! The row length in the vertical sub-stripe of all PEs
,    pe_rows_offset_h(:)                                                &
             ! The offset in rows of the sub-stripe inside a
             ! horizontal stripe for each PE
,    pe_row_length_offset_v(:)                                          &      
             ! The offset in row length of the sub-stripe inside a
             ! vertical stripe for each PE
,    pe_rows_displ_h(:)                                                 &
             ! The global displacement in rows of each sub-stripe
,    pe_row_length_displ_v(:)                                      
             ! The global displacement in row length of each sub-stripe


! Buffer for all-to-all commuication
REAL, ALLOCATABLE ::                                                    &
    sendbuf_l2h (:,:,:,:)                                               &
             ! Buffer used to send the local LPG data to horizontal
             ! sub-stripe
,   recvbuf_l2h (:,:,:,:)                                               &
             ! Buffer used to recv the local LPG data from horizontal
             ! sub-stripe                                                 
,   sendbuf_h2v (:,:,:,:)                                               &
             ! Buffer used to send the horizontal sub-stripe data to
             ! the vertical sub-stripe, and vice-versa
,   recvbuf_h2v (:,:,:,:)                                               &
             ! Buffer used to recv the horizontal sub-stripe data from
             ! the vertical sub-stripe, and vice-versa
,   sendbuf_v2l (:,:,:,:)                                               &
             ! Buffer used to send the vertical sub-stripe data to
             ! the local LPG domain
,   recvbuf_v2l (:,:,:,:)
             ! Buffer used to recv the vertical sub-stripe data from
             ! the local LPG domain


! Arguments for all-to-all communication
INTEGER ::                                                              &
     sendcount_l2h                                                      &
             ! The number of elements to send to each processor when
             ! transfering from LPG to horizontal sub-stripe
,    recvcount_l2h                                                      &
             ! The maximum number of elements that can be received
             ! from each processor when transfering from LPG to
             ! horizontal sub-stripe
,    sendcount_h2v                                                      &
             ! The number of elements to send to each processor when
             ! transfering from horizontal sub-stripe to vertical one,
             ! and vice-versa
,    recvcount_h2v                                                      &
             ! The maximum number of elements that can be received
             ! from each processor when transfering from horizontal
             ! sub-stripe to vertical one, and vice-versa
,    sendcount_v2l                                                      &
             ! The number of elements to send to each processor when
             ! transfering from vertical sub-stripe to LPG domain
,    recvcount_v2l                                                     
             ! The maximum number of elements that can be received
             ! from each processor when transfering from vertical
             ! sub-stripe to LPG domain.

! Local variable for stripes decompositions
INTEGER ::                                                              &
    stripe_rows_h(0:nproc_max)                                          &
             ! The rows in the  horizontal stripes 
,   stripe_row_length_v(0:nproc_max)                                    
             ! The row length in the vertical stripes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!------------------------------------------------------------------
! 1.0 Set up horizontal and vertical stripes domains information if it
! is the first call to the routine.
!
! ATTENTION: this routine has not been tested if:
!
! rows_l       /= g_blsize(2, fld_type_p, mype)
! row_length_l /= g_blsize(1, fld_type_p, mype)
! rows_g       /= glsize(2, fld_type_p)
! row_length_g /= glsize(1, fld_type_p)

IF(first_call) THEN

  ! Allocate saved arrays
  ALLOCATE(pe_rows_h(0:nproc-1))
  ALLOCATE(pe_row_length_v(0:nproc-1))
  ALLOCATE(pe_rows_offset_h(0:nproc-1))
  ALLOCATE(pe_row_length_offset_v(0:nproc-1))
  ALLOCATE(pe_rows_displ_h(0:nproc-1))
  ALLOCATE(pe_row_length_displ_v(0:nproc-1))

  max_rows_l       = 0
  max_row_length_l = 0

  DO pe = 0, nproc-1
    ! Number of rows of the horizontal stripe where the specific PE is
    ! located, this is the same number of rows that the PE has in the
    ! LPG decomposition
    stripe_rows_h(pe) = g_blsize(2, fld_type_p, pe)

    ! Row length of the vertical stripe where the specific PE is
    ! located, this is the same row_length that the PE has in the
    ! LPG decomposition
    stripe_row_length_v(pe) = g_blsize(1, fld_type_p, pe)
    
    max_rows_l       = MAX(max_rows_l,       g_blsize(2, fld_type_p, pe))
    max_row_length_l = MAX(max_row_length_l, g_blsize(1, fld_type_p, pe))
    
  END DO

  max_rows_h       = 0
  max_row_length_v = 0

  ! Compute the size  of the horiz/vertical sub-stripes
  DO pe_j = 0, nproc_y-1
    DO pe_i = 0, nproc_x-1
      
      ! The rank 
      pe=pe_i+pe_j*nproc_x
      
      ! Compute the number of rows of the horizontal sub-stripe and
      ! the row length of the vertical sub-stripe  on each PE. 
      pe_rows_h(pe)       = stripe_rows_h(pe)      /nproc_x  
      pe_row_length_v(pe) = stripe_row_length_v(pe)/nproc_y  
      
      ! If the number of rows in a horizontal stripe cannot be evenly
      ! divided between all the sub-stripes, add one extra row in the
      ! last few sub-stripes as needed. Thus the sum of rows between
      ! all the sub-stripes is equal to the number of rows of the
      ! stripe. Do the same for row length on vertical
      ! stripes/sub-stripes.
      IF (pe_i >= nproc_x - MOD(stripe_rows_h(pe),nproc_x)) THEN
        pe_rows_h(pe) = pe_rows_h(pe) + 1
      END IF
      IF (pe_j >= nproc_y - MOD(stripe_row_length_v(pe),nproc_y)) THEN
        pe_row_length_v(pe) = pe_row_length_v(pe) + 1
      END IF
      
      max_rows_h       = MAX(max_rows_h,       pe_rows_h(pe))
      max_row_length_v = MAX(max_row_length_v, pe_row_length_v(pe))
      
    END DO
  END DO
  
  ! Compute offset
  DO pe_j = 0, nproc_y-1
    DO pe_i = 0, nproc_x-1
      
      ! The rank 
      pe=pe_i+pe_j*nproc_x
      
      IF (pe_i == 0) THEN
        pe_rows_offset_h(pe) = 0
      ELSE
        ! Previous PE on the same row
        prev_pe = (pe_i-1)+pe_j*nproc_x
        pe_rows_offset_h(pe) = pe_rows_offset_h(prev_pe) +                 &
          pe_rows_h(prev_pe)
      END IF
      
      IF (pe_j == 0) THEN
        pe_row_length_offset_v(pe) = 0
      ELSE
        ! Previous PE on the same column
        prev_pe = pe_i+(pe_j-1)*nproc_x
        pe_row_length_offset_v(pe) = pe_row_length_offset_v(prev_pe) +     &
          pe_row_length_v(prev_pe)
      END IF
    END DO
  END DO
  
  ! Compute horizontal sub stripes global displacement
  prev_pe = 0
  DO pe_j = 0, nproc_y-1
    DO pe_i = 0, nproc_x-1
      
      ! The rank 
      pe=pe_i+pe_j*nproc_x
      
      IF (pe == 0) THEN
        pe_rows_displ_h(pe) = 0
      ELSE
        pe_rows_displ_h(pe) = pe_rows_displ_h(prev_pe) +                 &
          pe_rows_h(prev_pe)
      END IF
      prev_pe = pe
    END DO
  END DO
  
  ! Compute vertical sub stripes global displacement
  prev_pe = 0
  DO pe_i = 0, nproc_x-1
    DO pe_j = 0, nproc_y-1
      ! The rank 
      pe=pe_i+pe_j*nproc_x
      
      IF (pe == 0) THEN
        pe_row_length_displ_v(pe) = 0
      ELSE
        pe_row_length_displ_v(pe) = pe_row_length_displ_v(prev_pe) +     &
          pe_row_length_v(prev_pe)
      END IF
      prev_pe = pe
    END DO
  END DO
  
  ! Local sizes
  rows_h       = pe_rows_h(mype)
  row_length_v = pe_row_length_v(mype)

  first_call = .FALSE.

END IF

!------------------------------------------------------------------
! 1.1 Allocate horizontal/vertical stripes and communication buffers


! Allocate horizonatal and vertical sub-stripes field
ALLOCATE(f_h(row_length_g, rows_h, model_levels))
ALLOCATE(f_v(row_length_v, rows_g, model_levels))

! Allocate communication buffers
ALLOCATE(sendbuf_l2h(max_row_length_l,max_rows_h,model_levels,0:nproc_x-1))
ALLOCATE(recvbuf_l2h(max_row_length_l,max_rows_h,model_levels,0:nproc_x-1))

ALLOCATE(sendbuf_h2v(max_row_length_v,max_rows_h,model_levels,0:nproc-1))
ALLOCATE(recvbuf_h2v(max_row_length_v,max_rows_h,model_levels,0:nproc-1))

ALLOCATE(sendbuf_v2l(max_row_length_v,max_rows_l,model_levels,0:nproc_y-1))
ALLOCATE(recvbuf_v2l(max_row_length_v,max_rows_l,model_levels,0:nproc_y-1))


! The number of elements to send and receive
sendcount_l2h=max_row_length_l*max_rows_h*model_levels
recvcount_l2h=max_row_length_l*max_rows_h*model_levels
sendcount_h2v=max_row_length_v*max_rows_h*model_levels
recvcount_h2v=max_row_length_v*max_rows_h*model_levels
sendcount_v2l=max_row_length_v*max_rows_l*model_levels
recvcount_v2l=max_row_length_v*max_rows_l*model_levels  

!------------------------------------------------------------------
! 2.0 Calculate coeffs for E->W sweep (also used for w->e)
!   

DO j=1,rows_g  
  DO i=1,row_length_g
    z=npass*(dlong_g(i)*coslat_g(j))**2*0.25
    e=z/fi_scale**2
    a_we(i,j) = 1.0+e-SQRT(e*(e+2.0))
  END DO
END DO

!------------------------------------------------------------------
! 2.1 Calculate coeffs for N->S sweep (also used for S->N)
!

DO j=1,rows_g
  z=npass*dlat_g(j)**2*0.25
  e=z/fi_scale**2
  a_ns(j)     = 1.0+e-SQRT(e*(e+2.0))
END DO ! jrow


!------------------------------------------------------------------
! 3.0 Perform the filter passes


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, jpass, pe, pe_j, pe_i,  & 
!$OMP dstart)

DO jpass = 1, npass     

  !------------------------------------------------------------------
  ! 3.1 Transfer data to horizontal stripes depending where it is
  ! located at the moment, if jpass>1 then the data is in the vertical
  ! stripe
  IF ( jpass == 1 ) THEN
    ! ----------- LPG  -> Horizontal -----------------------------    

    ! Copy data from local LPG field to the send buffer 
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels           
      ! For each PE in a LPG row
      pe_j = gridpos(2)
      DO pe_i = 0, nproc_x-1   
        pe     = pe_i + pe_j*nproc_x          ! global rank
        dstart = pe_rows_offset_h(pe)         ! rows offset
        DO j = 1, pe_rows_h(pe)            
          DO i = 1, row_length_l
            sendbuf_l2h(i,j,k,pe_i) = f_l(i,dstart+j,k)      
          END DO ! i
        END DO ! j
      END DO ! pe
    END DO ! k  
!$OMP END DO 

    ! Transfer from local LPG to horizontal stripes
!$OMP MASTER
    CALL mpl_alltoall(sendbuf_l2h, sendcount_l2h, mpl_real,  &
      recvbuf_l2h, recvcount_l2h, mpl_real,                  &
      gc_proc_row_group, ierr)
!$OMP END MASTER
!$OMP BARRIER

    ! Copy data from receive buffer to horizontal sub-stripe
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels 

      ! For each PE in a LPG row
      pe_j = gridpos(2)
      DO pe_i = 0, nproc_x-1    
        pe = pe_i + pe_j*nproc_x              ! global rank
        dstart = g_datastart(1, pe) - 1       ! row length offset 
        DO j = 1, rows_h       
          DO i = 1, g_blsize(1, fld_type_p, pe)
            f_h(dstart+i,j,k) = recvbuf_l2h(i,j,k,pe_i) 
          END DO ! i
        END DO ! j
      END DO ! pe
    END DO ! k  
!$OMP END DO NOWAIT

  ELSE ! jpass > 1

    ! ----------- Vertical  -> Horizontal -----------------------------

    ! Copy data from vertical sub-stripe field to the send buffer 
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels 

      ! For each PE 
      DO pe = 0, nproc-1
        dstart = pe_rows_displ_h(pe)         ! rows offset global displacement
        DO j = 1, pe_rows_h(pe)
          DO i = 1, row_length_v                       
            sendbuf_h2v(i,j,k,pe)  = f_v(i,dstart+j,k)
          END DO ! i
        END DO ! j
      END DO ! pe
    END DO ! k  
!$OMP END DO

    ! Transfer from vertical sub-stripes to horizontal sub-stripes
!$OMP MASTER
    CALL mpl_alltoall(sendbuf_h2v, sendcount_h2v, mpl_real,            &
      recvbuf_h2v, recvcount_h2v, mpl_real,            &
      gc_all_proc_group, ierr)
!$OMP END MASTER
!$OMP BARRIER

    ! Copy data from receive buffer to horizontal sub-stripe
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels 

      ! For each PE 
      DO pe = 0, nproc-1
        dstart = pe_row_length_displ_v(pe)    ! row length global displacement
        DO j = 1, rows_h
          DO i = 1, pe_row_length_v(pe)                  
            f_h(dstart+i,j,k) = recvbuf_h2v(i,j,k,pe)      
          END DO ! i
        END DO ! j
      END DO ! pe
    END DO ! k  
!$OMP END DO NOWAIT

  END IF ! jpass

  ! Rows offset global displacement
  dstart = pe_rows_displ_h(mype)  

  ! The filtering of the levels is independent from each other.
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels           

    !------------------------------------------------------------------
    ! 3.2 Filter W->E
    
    !------------------------------------------------------------------
    ! 3.2.1 Boundary
    IF (m_grid == 0) THEN
      ! Perform W BC
      ! The boundary conditions depend on the number of previous
      ! passes in the opposite direction.
      
      ! IF f_boundary=0 the formulae for 0,1,2 previous passes are:-
      ! G1=(1-A)F1                          (6.5.26)
      ! G1=(1/(1+A))F1                      (6.5.27)
      ! G1=(1/(1+A)(1-A**2))(F1-A**3*F2)    (6.5.28)
      ! If f_boundary /= 0 it must be subtracted before using formulae.
      
      IF (jpass == 1) THEN              ! No previous E->W pass        
        
        DO j = 1, rows_h
          f_h(1,j,k)=f_boundary+(f_h(1,j,k)-f_boundary)                  &
            *(1.0-a_we(1,dstart+j))
        END DO ! j

      ELSE IF (jpass == 2) THEN         ! One previous E->W pass
      
        DO j= 1, rows_h
          f_h(1,j,k)=f_boundary+(f_h(1,j,k)-f_boundary)                  &
            /(1.0+a_we(1,dstart+j))
        END DO ! j
      
      ELSE IF (jpass >= 3) THEN         ! Two previous E->W passes
                                        ! BC for 2 is also used for >2   
        DO j = 1, rows_h
          f_h(1,j,k)=f_boundary+(f_h(1,j,k)-f_boundary-                  &
            a_we(1,dstart+j)**3*(f_h(2,j,k)-f_boundary))                  &
            /((1.0+a_we(1,dstart+j))*(1.0-a_we(1,dstart+j)**2))
        END DO ! j
      END IF ! ! Check jpass
    END IF ! m_grid

    !------------------------------------------------------------------
    ! 3.2.2 Filter W->E EQN 6.5.1
    
    DO j= 1, rows_h
      DO i= 2, row_length_g
        f_h(i,j,k) = f_h(i,j,k)+                                         &
          (f_h(i-1,j,k)-f_h(i,j,k))*a_we(i,dstart+j)
      END DO ! i
    END DO ! 

    !------------------------------------------------------------------
    ! 3.3 Filter E->W
    
    
    !------------------------------------------------------------------
    ! 3.3.1 Boundary
    IF (m_grid == 0) THEN
      ! Perform E BC

      IF (jpass == 1) THEN              ! One previous W->E pass

        DO j= 1, rows_h
          f_h(row_length_g,j,k)=f_boundary+(f_h(row_length_g,j,k)-f_boundary) &
                              /(1.0+a_we(row_length_g,dstart+j))
        END DO ! j
        
      ELSE IF (jpass >= 2) THEN         ! Two previous W->E passes
                                        ! BC for 2 is also used for >2

        DO j = 1, rows_h
          f_h(row_length_g,j,k)=f_boundary+(f_h(row_length_g,j,k) -        &
            f_boundary - a_we(row_length_g,dstart+j)**3 *                  &
            (f_h(row_length_g-1,j,k)-f_boundary))                          &
            /((1.0+a_we(row_length_g,dstart+j))                            &
            *(1.0-a_we(row_length_g,dstart+j)**2))
        END DO ! j

      END IF ! Check jpass
    END IF ! m_grid
  
    !------------------------------------------------------------------
    ! 3.3.2 Filter E->W EQN 6.5.1
    
    DO j = 1, rows_h
      DO i = row_length_g - 1, 1, -1
        f_h(i,j,k)=f_h(i,j,k)+                                             &
          (f_h(i+1,j,k)-f_h(i,j,k))*a_we(i,dstart+j)
      END DO ! i
    END DO ! j

  END DO ! k  ! End W-E-W filter
!$OMP END DO NOWAIT

  !------------------------------------------------------------------
  ! 3.4 Transfer data to vertical stripes 

  ! ----------- Horizontal -> Vertical -----------------------------

  ! Copy data from horizontal sub-stripe field to the send buffer 
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels 

    ! For each PE 
    DO pe = 0, nproc-1
      dstart = pe_row_length_displ_v(pe)    ! row length global displacement
      DO j = 1, rows_h
        DO i = 1, pe_row_length_v(pe)                  
          sendbuf_h2v(i,j,k,pe) = f_h(dstart+i,j,k)      
        END DO ! i
      END DO ! j
    END DO ! pe
  END DO ! k  
!$OMP END DO 

  ! Transfer from horizontal sub-stripes to vertical sub-stripes
!$OMP MASTER
  CALL mpl_alltoall(sendbuf_h2v, sendcount_h2v, mpl_real,            &
    recvbuf_h2v, recvcount_h2v, mpl_real,            &
    gc_all_proc_group, ierr)
!$OMP END MASTER
!$OMP BARRIER

  ! Copy data from receive buffer to vertical sub-stripe
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels 

    ! For each PE 
    DO pe = 0, nproc-1
      dstart = pe_rows_displ_h(pe)         ! rows offset global displacement
      DO j = 1, pe_rows_h(pe)
        DO i = 1, row_length_v                       
          f_v(i,dstart+j,k) = recvbuf_h2v(i,j,k,pe) 
        END DO ! i
      END DO ! j
    END DO ! pe
  END DO ! k  
!$OMP END DO NOWAIT

  ! The filtering of the levels is independent from each other.
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels           

    !------------------------------------------------------------------
    ! 3.5 Filter S->N
    
    !------------------------------------------------------------------
    ! 3.5.1 Boundary
    IF (m_grid == 0) THEN
      ! Perform S BC

      IF (jpass == 1) THEN              ! No previous S->N pass
        
        DO i = 1, row_length_v
          f_v(i,1,k)=f_boundary+(f_v(i,1,k)-f_boundary)                    &
            *(1.0-a_ns(1))
        END DO ! i

      ELSE IF (jpass == 2) THEN         ! One previous S->N pass

        DO i = 1, row_length_v
          f_v(i,1,k)=f_boundary+(f_v(i,1,k)-f_boundary)                    &
            /(1.0+a_ns(1))
        END DO ! i

      ELSE IF (jpass >= 3) THEN         ! Two previous S->N passes
                                        ! BC for 2 is also used for >2
        DO i = 1, row_length_v
          f_v(i,1,k)=f_boundary+(f_v(i,1,k)-f_boundary-                    &
            a_ns(1)**3*(f_v(i,2,k)-f_boundary))                            &
            /((1.0+a_ns(1))*(1.0-a_ns(1)**2))
        END DO ! i

      END IF ! Check  jpass
    END IF ! m_grid

    !------------------------------------------------------------------
    ! 3.5.2 Filter 
    DO j= 2, rows_g
      DO i= 1, row_length_v
        f_v(i,j,k)=f_v(i,j,k)+                                             &
          (f_v(i,j-1,k)-f_v(i,j,k))*a_ns(j)
      END DO ! i
    END DO ! j

    !------------------------------------------------------------------
    ! 3.6 Filter N->S
    
    !------------------------------------------------------------------
    ! 3.6.1 Boundary
    IF (m_grid == 0) THEN
      ! Perform N BC

      IF (jpass == 1) THEN              ! One previous N->S pass
                                        ! BC for after a single sweep is
                                        ! used for all later sweeps
        DO i = 1, row_length_v
          f_v(i,rows_g,k)=f_boundary+(f_v(i,rows_g,k)-f_boundary)          &
                          /(1.0+a_ns(rows_g))
        END DO ! i

      ELSE IF (jpass >= 2) THEN         ! Two previous S->N passes
                                        ! BC for 2 is also used for >2
        DO i = 1, row_length_v
          f_v(i,rows_g,k)=f_boundary+(f_v(i,rows_g,k)-f_boundary-          &
            a_ns(rows_g)**3*(f_v(i,rows_g-1,k)-f_boundary))                &
            /((1.0+a_ns(rows_g))*(1.0-a_ns(rows_g)**2))
        END DO ! i

      END IF ! Check  jpass
    END IF ! m_grid
    
    !------------------------------------------------------------------
    ! 3.6.2 Filter 
    DO j= rows_g-1, 1, -1
      DO i= 1, row_length_v
        f_v(i,j,k)=f_v(i,j,k)+                                             &
          (f_v(i,j+1,k)-f_v(i,j,k))*a_ns(j)
      END DO ! i
    END DO ! j


  END DO ! k  ! End S-N-S filter
!$OMP END DO NOWAIT

END DO ! Loop jpass

!------------------------------------------------------------------
! 4.0 Filter passes finished, transfer data back to LPG domain


! ----------- Vertical -> LPG -----------------------------    

! Copy data from local vertical sub stripe to the send buffer 
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels 

  ! For each PE in a LPG column
  pe_i = gridpos(1)
  DO pe_j = 0, nproc_y-1   
    pe     = pe_i + pe_j*nproc_x          ! global rank
    dstart = g_datastart(2, pe) - 1       ! rows offset 
    DO j = 1, g_blsize(2, fld_type_p, pe)
      DO i = 1, row_length_v
        sendbuf_v2l(i,j,k,pe_j) = f_v(i,dstart+j,k)      
      END DO ! i
    END DO ! j
  END DO ! pe
END DO ! k  
!$OMP END DO

! Transfer from vertical stripes to local LPG. 
!$OMP MASTER
CALL mpl_alltoall(sendbuf_v2l, sendcount_v2l, mpl_real,            &
  recvbuf_v2l, recvcount_v2l, mpl_real,            &
  gc_proc_col_group, ierr)
!$OMP END MASTER
!$OMP BARRIER

! Copy data from receive buffer to local LPG domain
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels 

  ! For each PE in a LPG column
  pe_i = gridpos(1)
  DO pe_j = 0, nproc_y-1   
    pe = pe_i + pe_j*nproc_x              ! global rank
    dstart = pe_row_length_offset_v(pe)   ! row length offset 
    DO j = 1, rows_l
      DO i = 1, pe_row_length_v(pe)
        f_l(dstart+i,j,k) = recvbuf_v2l(i,j,k,pe_j)
      END DO ! i
    END DO ! j
  END DO ! pe
END DO ! k  
!$OMP END DO NOWAIT

!$OMP END PARALLEL

DEALLOCATE(recvbuf_v2l)
DEALLOCATE(sendbuf_v2l)

DEALLOCATE(recvbuf_h2v)
DEALLOCATE(sendbuf_h2v)

DEALLOCATE(recvbuf_l2h)
DEALLOCATE(sendbuf_l2h)


DEALLOCATE(f_v)
DEALLOCATE(f_h)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rfcsl_1a

END MODULE rfcsl_1a_mod
