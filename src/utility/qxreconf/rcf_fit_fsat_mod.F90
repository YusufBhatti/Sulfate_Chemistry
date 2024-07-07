! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up Fitting parameters for large-scale hydrology

MODULE Rcf_Fit_Fsat_Mod

! Subroutine Rcf_Fit_Fsat
!
! Description:
!   Calls Calc_Fit_Fsat which calculates the fitting parameters for
!   LSH model and stores them.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FIT_FSAT_MOD'

CONTAINS

SUBROUTINE Rcf_Fit_Fsat( fields_out, field_count_out,        &
                         field_pos, hdr_out )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY:  &
    stashcode_sthzw,            &
    stashcode_vol_smc_sat,      &
    stashcode_unfrozen_soil,    &
    stashcode_ti_mean,          &
    stashcode_ti_sig,           &
    stashcode_gamtot,           &
    stashcode_fexp,             &
    stashcode_a_fsat,           &
    stashcode_c_fsat,           &
    stashcode_a_fwet,           &
    stashcode_c_fwet,           &
    stashcode_prog_sec

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE nlstcall_mod, ONLY:   &
    LTimer

USE UM_ParCore, ONLY: &
    mype,             &
    nproc

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Ereport_mod, ONLY: &
    Ereport

USE jules_hydrology_mod, ONLY: l_top

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_out

INTEGER, INTENT(IN)                :: field_pos
INTEGER, INTENT(IN)                :: field_count_out

! Local variables
TYPE( field_type ), POINTER          :: vol_smc_sat
TYPE( field_type ), POINTER          :: unfrozen_soil
TYPE( field_type ), POINTER          :: ti_mean
TYPE( field_type ), POINTER          :: ti_sig
TYPE( field_type ), POINTER          :: gamtot
TYPE( field_type ), POINTER          :: fexp

INTEGER                              :: soil_index  &
                                       (fields_out(field_pos) % level_size)
INTEGER                              :: soil_pts

INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper
INTEGER                              :: n    ! looper
INTEGER                              :: j    ! looper


! These fields will be saved
REAL, ALLOCATABLE, SAVE  :: a_fsat (:,:)
REAL, ALLOCATABLE, SAVE  :: c_fsat (:,:)
REAL, ALLOCATABLE, SAVE  :: a_fwet (:,:)
REAL, ALLOCATABLE, SAVE  :: c_fwet (:,:)

INTEGER      :: ErrorStatus

CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_FIT_FSAT'

REAL                                 :: soil_depth
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( RoutineName, 3)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (.NOT. l_top) THEN
  ErrorStatus = 30
  WRITE(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
  CALL Ereport ( RoutineName, ErrorStatus, CMessage)
END IF

! Only calculate when all fields are not allocated.
IF (.NOT. ALLOCATED(a_fsat) .AND. &
    .NOT. ALLOCATED(c_fsat) .AND. &
    .NOT. ALLOCATED(a_fwet) .AND. &
    .NOT. ALLOCATED(c_fwet) ) THEN

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) 'Setting up call for calc_fit_fsat in rcf_fit_fsat'
    CALL umPrint(umMessage,src='rcf_fit_fsat_mod')
  END IF

  ! Allocate all fields.
  ALLOCATE(a_fsat( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(c_fsat( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(a_fwet( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(c_fwet( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )

  !----------------------------------------------------------------------
  ! Find required fields in output dump and read them in:
  !----------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,            &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  CALL Rcf_Alloc_Field( vol_smc_sat )
  CALL Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,          &
                   fields_out, field_count_out, pos)
  unfrozen_soil => fields_out(pos)
  CALL Rcf_Alloc_Field( unfrozen_soil )
  CALL Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,                &
                   fields_out, field_count_out, pos)
  ti_mean => fields_out(pos)
  CALL Rcf_Alloc_Field( ti_mean )
  CALL Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,                 &
                   fields_out, field_count_out, pos)
  ti_sig => fields_out(pos)
  CALL Rcf_Alloc_Field( ti_sig )
  CALL Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,                 &
                   fields_out, field_count_out, pos)
  gamtot => fields_out(pos)
  CALL Rcf_Alloc_Field( gamtot )
  CALL Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_fexp,                   &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  CALL Rcf_Alloc_Field( fexp )
  CALL Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

  !----------------------------------------------------------------------
  ! Initialise fitting variables to zero:
  !----------------------------------------------------------------------
  a_fsat(:,:) = 0.0
  c_fsat(:,:) = 0.0
  a_fwet(:,:) = 0.0
  c_fwet(:,:) = 0.0

  !----------------------------------------------------------------------
  ! Set up soil index:
  !----------------------------------------------------------------------

  soil_pts=0
  soil_index(:)=0
  DO i=1 , vol_smc_sat % level_size
    IF (vol_smc_sat % DATA(i,1) > 0.0) THEN
      soil_pts = soil_pts + 1
      soil_index(soil_pts) = i
    END IF
  END DO

  Soil_depth=0.0
  DO n = 1,unfrozen_soil % levels
    soil_depth=soil_depth+Output_Grid % soil_depths(n)
  END DO

  !----------------------------------------------------------------------
  ! Set up variables needed for call to calc_fit_fsat
  !----------------------------------------------------------------------

  DO j = 1, soil_pts
    i = soil_index(j)
    IF (gamtot % DATA(i,1) <= 0.0) THEN
      ErrorStatus = 20
      WRITE (CMessage, '(A)')'Gamma integral should be positive and non-zero.'
      CALL Ereport ( RoutineName, ErrorStatus, CMessage)
    END IF

    IF (fexp % DATA(i,1) <= 0.0) THEN
      ErrorStatus = 21
      WRITE (CMessage, '(A)')'Exp. decay in Vol_Smc_Sat wrongly/not set'
      CALL Ereport ( RoutineName, ErrorStatus, CMessage)
    END IF

  END DO

  IF (nproc > 1) THEN
! Note RECON_SERIAL is compiled without MPI support, so all MPI handling
! needs to be deffed out
#if !defined(RECON_SERIAL)
    ! Do load-balancing
    CALL rcf_balance_fit_fsat(   soil_pts,   soil_index,         &
                                 ti_mean % level_size,           &
                                 a_fsat,     c_fsat,             &
                                 a_fwet,     c_fwet,             &
                                 ti_mean,    fexp,               &
                                 ti_sig,     gamtot,             &
                                 soil_depth)
#else
    
    ErrorStatus = 22
    WRITE (CMessage, '(A)')'Reconfiguration compiled with RECON_SERIAL def' //&
                           'but using > 1 processes'
    CALL Ereport ( RoutineName, ErrorStatus, CMessage)
#endif

  ELSE
    
  !----------------------------------------------------------------------
  ! Get fitting parameters:
  !----------------------------------------------------------------------

    ! DEPENDS ON: calc_fit_fsat
    CALL Calc_Fit_Fsat(                       &
                       soil_pts,              &
                       soil_index,            &
                       ti_mean % level_size,  &
                       fexp % DATA(:,1),      &
                       ti_mean % DATA(:,1),   &
                       ti_sig % DATA(:,1),    &
                       gamtot % DATA(:,1),    &
                       soil_depth,            &
                       a_fsat,                &
                       c_fsat,                &
                       a_fwet,                &
                       c_fwet)

  END IF   ! nproc > 1
  !----------------------------------------------------------------------
  ! Tidy Up
  !----------------------------------------------------------------------
  CALL Rcf_Dealloc_Field( vol_smc_sat )
  CALL Rcf_Dealloc_Field( unfrozen_soil )
  CALL Rcf_Dealloc_Field( ti_mean )
  CALL Rcf_Dealloc_Field( ti_sig )
  CALL Rcf_Dealloc_Field( gamtot )
  CALL Rcf_Dealloc_Field( fexp )
END IF

SELECT CASE( fields_out( field_pos ) % stashmaster % item )

CASE ( stashcode_a_fsat)
  fields_out(field_pos) % DATA (:,:) = a_fsat (:,:)
  DEALLOCATE(a_fsat)

CASE ( stashcode_c_fsat )
  fields_out(field_pos) % DATA (:,:) = c_fsat (:,:)
  DEALLOCATE(c_fsat)

CASE ( stashcode_a_fwet)
  fields_out(field_pos) % DATA (:,:) = a_fwet (:,:)
  DEALLOCATE(a_fwet)

CASE ( stashcode_c_fwet )
  fields_out(field_pos) % DATA (:,:) = c_fwet (:,:)
  DEALLOCATE(c_fwet)

END SELECT

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Fit_Fsat

!------------------------------------------------------------------------------

#if !defined(RECON_SERIAL)
SUBROUTINE rcf_balance_fit_fsat( soil_pts,   soil_index,         &
                                 level_size,                     &
                                 a_fsat,     c_fsat,             &
                                 a_fwet,     c_fwet,             &
                                 ti_mean,    fexp,               &
                                 ti_sig,     gamtot,             &
                                 soil_depth)

! Subroutine Rcf_balance_fit_fsat
!
! Description:
!   Performs load-balancing of the expensive fit-fsat calaculations
!
! Method:
!  - Compress soil fields to just the soil points
!  - Calculate a new equal distribution of the soil points
!  - Calculate maps of how many points to distribute to/from each process
!  - Fill buffers with data to send
!  - Communicate via alltoall
!  - Unpack buffers and existing data into new contiguous arrays
!  - Do the main calculation
!  - Re-pack buffers
!  - Communicate results back to original processes
!  - Unpack buffers back into the compressed fields
!  - Uncompress fields and continue

USE um_parcore, ONLY: mype, nproc

USE rcf_field_type_mod, ONLY: field_type

USE mpl, ONLY : mpl_integer, mpl_real

IMPLICIT NONE
! Arguments
INTEGER, INTENT(IN)  :: soil_pts               ! number of soil points
INTEGER, INTENT(IN)  :: level_size             ! 2d level size
REAL, INTENT(IN)     :: soil_depth             ! soil depth
INTEGER, INTENT(IN)  :: soil_index(level_size) ! mapping of xy to soil points
TYPE( field_type ), INTENT(IN) :: ti_mean
TYPE( field_type ), INTENT(IN) :: fexp
TYPE( field_type ), INTENT(IN) :: ti_sig
TYPE( field_type ), INTENT(IN) :: gamtot

REAL, INTENT(INOUT)  :: a_fsat (level_size,1)
REAL, INTENT(INOUT)  :: c_fsat (level_size,1)
REAL, INTENT(INOUT)  :: a_fwet (level_size,1)
REAL, INTENT(INOUT)  :: c_fwet (level_size,1)


! Local variables
INTEGER              :: i               ! looper
INTEGER              :: j               ! looper
INTEGER              :: pos             ! position in buffer
INTEGER              :: comm            ! Communicator
INTEGER              :: ierr            ! Return from MPL
INTEGER              :: mpl_a2a_type    ! MPI type for alltoall
INTEGER, ALLOCATABLE :: sizes(:)        ! How many soil pts on each process
INTEGER, ALLOCATABLE :: new_sizes(:)    ! Distributed numbers of soil pts
INTEGER, ALLOCATABLE :: map(:,:)        ! map of send/recv values
INTEGER              :: offsets         ! offsets in MPI data type
INTEGER              :: oldtypes        ! old data types used in composite
INTEGER              :: blockcounts     ! number of each old datatype used
INTEGER              :: total_size      ! global number of soil points
INTEGER              :: current_size    ! intermediate size value
INTEGER              :: mysize          ! size of data once re-distributed
INTEGER              :: excess          ! excess in redistribute
INTEGER              :: diff            ! difference between current and needed
INTEGER              :: orig_disps(0:nproc-1)  ! original alltoall displacements
INTEGER              :: rem_disps(0:nproc-1)   ! remote alltoall displacements
INTEGER              :: orig_counts(0:nproc-1) ! original alltoall counts
INTEGER              :: rem_counts(0:nproc-1)  ! remote alltoall counts

REAL, ALLOCATABLE :: fexp_r(:)         !
REAL, ALLOCATABLE :: ti_mean_r(:)      !
REAL, ALLOCATABLE :: ti_sig_r(:)       !  The _r versions are the 
REAL, ALLOCATABLE :: gamtot_r(:)       !  redistributed variables used
REAL, ALLOCATABLE :: a_fsat_r(:)       !  in the load-balanced case
REAL, ALLOCATABLE :: c_fsat_r(:)       !
REAL, ALLOCATABLE :: a_fwet_r(:)       !
REAL, ALLOCATABLE :: c_fwet_r(:)       !

REAL, ALLOCATABLE :: fexp_c(:)         !
REAL, ALLOCATABLE :: ti_mean_c(:)      !
REAL, ALLOCATABLE :: ti_sig_c(:)       ! The _c versions are the variables
REAL, ALLOCATABLE :: gamtot_c(:)       ! compressed onto just the soil points
REAL, ALLOCATABLE :: a_fsat_c(:)       !
REAL, ALLOCATABLE :: c_fsat_c(:)       !
REAL, ALLOCATABLE :: a_fwet_c(:)       !
REAL, ALLOCATABLE :: c_fwet_c(:)       !

INTEGER, ALLOCATABLE :: s_index_r(:)   ! The compressed, redistributed
                                       ! soil index

! Type for temporary store of alltoall data. Component names are generic
! fields (f1, f2, ...) rather than being related to the data sent/received
! as different fields are distributed at the start than are gathered back 
! at the end.
TYPE a2a_type
  SEQUENCE
  REAL :: f1    ! Field 1
  REAL :: f2    ! Field 2 
  REAL :: f3    ! Field 3
  REAL :: f4    ! Field 4
END TYPE a2a_type

TYPE (a2a_type), ALLOCATABLE :: orig_buffer(:)  ! Buffer on original state
TYPE (a2a_type), ALLOCATABLE :: rem_buffer(:)   ! Buffer for remote state

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_BALANCE_FIT_FSAT'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays for soil-compressed data
ALLOCATE ( fexp_c(soil_pts) )
ALLOCATE ( ti_mean_c(soil_pts) )
ALLOCATE ( ti_sig_c(soil_pts) )
ALLOCATE ( gamtot_c(soil_pts) )
ALLOCATE ( a_fsat_c(soil_pts) )
ALLOCATE ( c_fsat_c(soil_pts) )
ALLOCATE ( a_fwet_c(soil_pts) )
ALLOCATE ( c_fwet_c(soil_pts) )

! Compress my soil point fields
DO j = 1, soil_pts
  i = soil_index(j)
  fexp_c(j)    = fexp    % data(i,1)
  ti_mean_c(j) = ti_mean % data(i,1)
  ti_sig_c(j)  = ti_sig  % data(i,1)
  gamtot_c(j)  = gamtot  % data(i,1)
END DO

CALL gc_get_communicator( comm, ierr )
ALLOCATE( sizes(0:nproc-1) )
ALLOCATE( new_sizes(0:nproc-1) )
ALLOCATE( map(0:nproc-1, 0:nproc-1) )

! Setup MPI type for alltoall
offsets     = 0
oldtypes    = mpl_real
blockcounts = 4
CALL mpl_type_create_struct(1, blockcounts, offsets, oldtypes,       &
                            mpl_a2a_type, ierr)
CALL mpl_type_commit(mpl_a2a_type, ierr)


! Gather values of soil_pts from all processes and put into array "sizes"
! so that all pes can work out which data needs to be communicated later
CALL mpl_allgather( soil_pts, 1, mpl_integer,  &
                    sizes,    1, mpl_integer,  &
                    comm,     ierr )

! Calculate new distribution. "new_sizes" will contain the number of points
! for an equally distributed set of fields
total_size = SUM(sizes)
DO i = 0, nproc-1
  new_sizes(i) = total_size/nproc
  IF (i < MOD(total_size, nproc)) new_sizes(i) = new_sizes(i) + 1
END DO

! Each row    map(:,i) is a sendmap    i.e. how much do pes send to pe i
! Each column map(j,:) is a receivemap i.e. how much do pes receive from pe j
! All process calculate the map for all processes.
map(:,:) = 0
DO  i = 0, nproc-1
  ! Does pe i have more data than it needs?
  IF (sizes(i) > new_sizes(i)) THEN
    excess = sizes(i) - new_sizes(i)
    ! PE i has "excess" more points than it needs. We'll now loop around
    ! the processes again to find those that need this excess data
    DO j = 0, nproc-1

      ! Calculate the amount of data pe j will have once all currently 
      ! known redistribution is taken into account
      current_size = sizes(j)+SUM(map(j,:)) 

      ! Does this process need any more data?
      IF (current_size < new_sizes(j) ) THEN
        ! The amount of data I can take will either be the amount of excess
        ! available from pe i, or pe j needs to to take it up to the new size
        diff = MIN(excess, new_sizes(j)-current_size)
        map(j,i) = map(j,i) + diff
        excess = excess - diff
      END IF
    END DO 
  END IF
END DO

! Setup comms buffers, counts and displacements
ALLOCATE( orig_buffer(SUM(map(:,mype))) )
ALLOCATE( rem_buffer(SUM(map(mype,:))) )
orig_disps(0) = 0
rem_disps(0)  = 0


! Counts are the number of points being sent to each process
! Displacements are where in the send buffer this data is (0 addressed)
DO  i = 0, nproc-1
  rem_counts(i)  = map(mype, i)
  orig_counts(i) = map(i, mype)
  IF (i>0) THEN
    orig_disps(i) = SUM(map(:i-1, mype))
    rem_disps(i)  = SUM(map(mype, :i-1))
  END IF
END DO

! Put data in send buffer
mysize = new_sizes(mype)
IF (mysize < soil_pts) THEN
  pos=1
  DO j = 0, nproc-1
    IF (map(j,mype) > 0) THEN   ! I've got data to send
      orig_buffer(pos : pos+map(j,mype)-1 ) % f1 =                        &
                  fexp_c   (mysize+pos : mysize+pos+map(j,mype)-1)
      orig_buffer(pos : pos+map(j,mype)-1 ) % f2 =                        &
                  ti_mean_c(mysize+pos : mysize+pos+map(j,mype)-1)
      orig_buffer(pos : pos+map(j,mype)-1 ) % f3 =                        &
                  ti_sig_c (mysize+pos : mysize+pos+map(j,mype)-1)
      orig_buffer(pos : pos+map(j,mype)-1 ) % f4 =                        &
                  gamtot_c (mysize+pos : mysize+pos+map(j,mype)-1)
      pos = pos + map(j,mype)
    END IF
  END DO
END IF

! Do the data exchange.
CALL mpl_alltoallv(orig_buffer, orig_counts, orig_disps, mpl_a2a_type, &
                   rem_buffer,  rem_counts,  rem_disps,  mpl_a2a_type, &
                   comm,        ierr)

! Allocate and Unpack the data into standard arrays
ALLOCATE( fexp_r   (mysize) )
ALLOCATE( ti_mean_r(mysize) )
ALLOCATE( ti_sig_r (mysize) )
ALLOCATE( gamtot_r (mysize) )
ALLOCATE( a_fsat_r (mysize) )
ALLOCATE( c_fsat_r (mysize) )
ALLOCATE( a_fwet_r (mysize) )
ALLOCATE( c_fwet_r (mysize) )
ALLOCATE( s_index_r(mysize) )

a_fsat_r(:) = 0.0
c_fsat_r(:) = 0.0
a_fwet_r(:) = 0.0
c_fwet_r(:) = 0.0

DO i = 1, mysize
  s_index_r(i) = i      ! we have a compressed set of data
END DO

IF ( mysize <= soil_pts ) THEN  ! Just a subset of original data
  fexp_r   (1:mysize) = fexp_c   (1:mysize)
  ti_mean_r(1:mysize) = ti_mean_c(1:mysize)
  ti_sig_r (1:mysize) = ti_sig_c (1:mysize)
  gamtot_r (1:mysize) = gamtot_c (1:mysize)
ELSE                           ! Mix of original and new data
  fexp_r   (1:soil_pts) = fexp_c   (1:soil_pts)
  ti_mean_r(1:soil_pts) = ti_mean_c(1:soil_pts)
  ti_sig_r (1:soil_pts) = ti_sig_c (1:soil_pts)
  gamtot_r (1:soil_pts) = gamtot_c (1:soil_pts)

  fexp_r   (soil_pts+1 : mysize) = rem_buffer(:) % f1
  ti_mean_r(soil_pts+1 : mysize) = rem_buffer(:) % f2
  ti_sig_r (soil_pts+1 : mysize) = rem_buffer(:) % f3
  gamtot_r (soil_pts+1 : mysize) = rem_buffer(:) % f4
END IF
 

  ! DEPENDS ON: calc_fit_fsat
CALL Calc_Fit_Fsat(                     &
                   mysize,                &
                   s_index_r,             &
                   mysize,                &
                   fexp_r,                &
                   ti_mean_r,             &
                   ti_sig_r,              &
                   gamtot_r,              &
                   soil_depth,            &
                   a_fsat_r,              &
                   c_fsat_r,              &
                   a_fwet_r,              &
                   c_fwet_r)

! Now we have to repeat most of the above to get the results back to 
! the original locations
IF ( mysize <= soil_pts ) THEN  ! Just a subset of original data
  a_fsat_c (1:mysize) = a_fsat_r (1:mysize)
  c_fsat_c (1:mysize) = c_fsat_r (1:mysize)
  a_fwet_c (1:mysize) = a_fwet_r (1:mysize)
  c_fwet_c (1:mysize) = c_fwet_r (1:mysize)
ELSE                           ! Mix of original and new data
  a_fsat_c (1:soil_pts) = a_fsat_r (1:soil_pts)
  c_fsat_c (1:soil_pts) = c_fsat_r (1:soil_pts)
  a_fwet_c (1:soil_pts) = a_fwet_r (1:soil_pts)
  c_fwet_c (1:soil_pts) = c_fwet_r (1:soil_pts)

  ! Note we'll be sending the remote buffer as this is the reverse
  ! journey
  rem_buffer(:) % f1 = a_fsat_r (soil_pts+1 : mysize)
  rem_buffer(:) % f2 = c_fsat_r (soil_pts+1 : mysize)
  rem_buffer(:) % f3 = a_fwet_r (soil_pts+1 : mysize)
  rem_buffer(:) % f4 = c_fwet_r (soil_pts+1 : mysize)
END IF

! Now the comms - everything is backwards from before as we are 
! returning the result
CALL mpl_alltoallv(rem_buffer,  rem_counts,  rem_disps,  mpl_a2a_type, &
                   orig_buffer, orig_counts, orig_disps, mpl_a2a_type, &
                   comm,        ierr)

! Unpack the results
IF (mysize < soil_pts) THEN
  pos=1
  DO j = 0, nproc-1
    IF (map(j,mype) > 0) THEN   ! I've received data in orig_buffer
      a_fsat_c(mysize+pos : mysize+pos+map(j,mype)-1) =                   &
        orig_buffer(pos : pos+map(j,mype)-1) % f1
      c_fsat_c(mysize+pos : mysize+pos+map(j,mype)-1) =                   &
        orig_buffer(pos : pos+map(j,mype)-1) % f2
      a_fwet_c(mysize+pos : mysize+pos+map(j,mype)-1) =                   &
        orig_buffer(pos : pos+map(j,mype)-1) % f3
      c_fwet_c(mysize+pos : mysize+pos+map(j,mype)-1) =                   &
        orig_buffer(pos : pos+map(j,mype)-1) % f4
      pos = pos + map(j,mype)
    END IF
  END DO
END IF

! And uncompress
DO j = 1, soil_pts
  i = soil_index(j)
  a_fsat(i,1) = a_fsat_c(j)
  c_fsat(i,1) = c_fsat_c(j)
  a_fwet(i,1) = a_fwet_c(j)
  c_fwet(i,1) = c_fwet_c(j)
END DO

! And deallocate my allocatables, including MPI types
DEALLOCATE( s_index_r )
DEALLOCATE( c_fwet_r )
DEALLOCATE( a_fwet_r )
DEALLOCATE( c_fsat_r )
DEALLOCATE( a_fsat_r )
DEALLOCATE( gamtot_r )
DEALLOCATE( ti_sig_r )
DEALLOCATE( ti_mean_r )
DEALLOCATE( fexp_r )
DEALLOCATE( rem_buffer )
DEALLOCATE( orig_buffer )
DEALLOCATE( map )
DEALLOCATE( new_sizes )
DEALLOCATE( sizes )
DEALLOCATE( c_fwet_c )
DEALLOCATE( a_fwet_c )
DEALLOCATE( c_fsat_c )
DEALLOCATE( a_fsat_c )
DEALLOCATE( gamtot_c )
DEALLOCATE( ti_sig_c )
DEALLOCATE( ti_mean_c )
DEALLOCATE( fexp_c )
  
CALL mpl_type_free(mpl_a2a_type, ierr)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_balance_fit_fsat
#endif

END MODULE Rcf_Fit_Fsat_Mod
