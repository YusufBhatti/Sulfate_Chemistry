! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT Met Office. ALL RIGHTS RESERVED.
! For further details please refer to the file copyright.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  scatters a field from one processor to many processors

! SUBROUTINE interface:
SUBROUTINE scatter_field_mpl32(local_field,     global_field,           &
                               local_row_len,   local_rows,             &
                               global_row_len,  global_rows,            &
                               grid_type,       halo_type,              &
                               scatter_pe,      proc_group )

USE mpl, ONLY: mpl_real4

USE um_types, ONLY: real32

USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE field_types, ONLY: fld_type_unknown

USE ereport_mod, ONLY: ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



! Description:
!  Takes a model field which is stored entirely on one processor
!  and distributes it over a group of processors using the
!  standard UM decomposition.

! Method:
!   Copies each processors data into contiguous space in a 1D
!   array and then uses mpl_scatterv to do the collective operation.


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! SUBROUTINE arguments:

INTEGER, INTENT(IN) :: local_row_len  ! length of rows in local part of field
INTEGER, INTENT(IN) :: local_rows     ! number of rows in local part of field
INTEGER, INTENT(IN) :: global_row_len ! length of rows in global field
INTEGER, INTENT(IN) :: global_rows    ! number of rows in global field
INTEGER, INTENT(IN) :: grid_type      ! type (p,u or v) of grid
INTEGER, INTENT(IN) :: halo_type      ! halo type (hence width) of grid
INTEGER, INTENT(IN) :: scatter_pe     ! processor to scatter global
                                      ! field from
INTEGER, INTENT(IN) :: proc_group     ! group id of processors involved
                                      ! here

REAL(KIND=real32), INTENT(OUT) :: local_field(local_row_len*local_rows)
                                      ! local part of field
REAL(KIND=real32), INTENT(INOUT) ::  global_field(global_row_len*global_rows)
                                      ! (on pe scatter_pe) global field



! local variables

REAL(KIND=real32) :: global_field_copy( global_row_len * global_rows )
                         ! Space for re-ordered copy of global field

INTEGER :: info          ! Status return from MPL/GCOM
INTEGER :: i,j, iproc    ! Loopers
INTEGER :: igpos         ! Position pointer
INTEGER :: ipos          ! Position pointer
INTEGER :: my_comm       ! Current communicator
INTEGER :: icode         ! return code

INTEGER :: counts(0:nproc-1)          ! Counts for Scatter
INTEGER :: displacements(0:nproc-1)   ! Displacements for Scatter

INTEGER, SAVE :: old_global_row_len = -1234
INTEGER, SAVE :: old_global_rows    = -1234
INTEGER, SAVE :: old_proc_group     = -1234
INTEGER, SAVE :: old_scatter_pe     = -1234
INTEGER, SAVE :: old_decomp         = -1234
INTEGER, SAVE :: old_grid_type      = -1234
INTEGER, SAVE :: old_halo_type      = -1234
INTEGER, SAVE :: recv_type          = -1234        ! receive data type

INTEGER, ALLOCATABLE :: iposa(:,:)
INTEGER :: iposa_size1          !  size of the first dimension in iposa
INTEGER :: g_blsize1(0:nproc-1)
                                !  X dimension of decomposition on each PE
INTEGER :: g_blsize2(0:nproc-1)
                                !  Y dimension of decomposition on each PE
INTEGER :: igpos1d(0:nproc-1)
INTEGER :: isum                 !  keep partial sums
INTEGER :: imax

CHARACTER(LEN=errormessagelength) :: cmessage              ! error message
CHARACTER (LEN=*), PARAMETER :: RoutineName='SCATTER_FIELD_MPL32'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!-------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

#if defined(RECON_SERIAL)
! We don't have any MPI access and just run on a single PE, so just
! copy data from input to output
local_field(1:local_row_len*local_rows) =                   &
            global_field(1:local_row_len*local_rows)

#else

! 0.0 Error checking.
! Do we recognise the field type?
IF (grid_type  ==  fld_type_unknown) THEN
  WRITE(umMessage,*) 'SCATTER_FIELD_MPL32 : bad field type'
  CALL umPrint(umMessage,src='scatter_field_mpl32')
  WRITE(umMessage,*) 'field will not be scattered.'
  CALL umPrint(umMessage,src='scatter_field_mpl32')
  WRITE (cmessage,*) 'SCATTER_FIELD_MPL32 : bad field type'
  icode = 1
  CALL ereport(routinename, icode, cmessage)
END IF

! We can only do scatters with global processor groups with this
! version currently
IF (proc_group /= gc_all_proc_group) THEN
  WRITE (cmessage, *) 'Can only scatter in an all processor group'
  icode = 2
  CALL ereport(routinename, icode, cmessage)
END IF

! 1.0 Setup MPI datatypes. Can we use the same as last time round?
IF ((global_row_len  /=  old_global_row_len) .OR.                     &
    (global_rows     /=  old_global_rows   ) .OR.                     &
    (proc_group      /=  old_proc_group    ) .OR.                     &
    (scatter_pe      /=  old_scatter_pe    ) .OR.                     &
    (grid_type       /=  old_grid_type     ) .OR.                     &
    (halo_type       /=  old_halo_type     ) .OR.                     &
    (current_decomp_type  /=  old_decomp  )) THEN

  ! Need to create a new receive type
  ! First get rid of the old one if this isn't the first call
  IF (recv_type /= -1234) THEN
    CALL mpl_type_free(recv_type, info)
  END IF

  ! New type is a 2D block with a local row length + 2 * halo stride.
  ! This is used to describe the data size and layout of the local
  ! data without including halos via an MPI datatype to simplify
  ! communications and reduce copies. Will need to specify the real
  ! start of the data in the communications call.
  CALL mpl_type_vector(blsize(2, grid_type),  blsize(1, grid_type),   &
                       local_row_len,  mpl_real4, recv_type, info)
  CALL mpl_type_commit(recv_type, info)

END IF  ! New MPI types

! 2.0 set up the scattering array, counts and displacements
!     Note that global data does not include halos.
!     This copy takes the conventional global field ordering and copies
!     into a per-PE contiguous data layout (all of PE0 data,
!     followed by all of PE1 data, followed by ...).
IF (mype == scatter_pe) THEN
  DO iproc = 0, nproc-1
    g_blsize1(iproc) = g_blsize(1, grid_type, iproc)
    g_blsize2(iproc) = g_blsize(2, grid_type, iproc)
    igpos1d(iproc)   = ( g_datastart_f(2,grid_type,iproc) -2 ) *      &
                         global_row_len - 1 +                         &
                         g_datastart_f(1,grid_type,iproc)
    counts(iproc)    = g_blsize1(iproc) * g_blsize2(iproc)
  END DO ! nprocs
  iposa_size1 = MAXVAL(g_blsize2)

  ALLOCATE( iposa(iposa_size1,0:nproc-1) )

  ipos = 0 ! start of array
  DO iproc = 0, nproc-1
    DO j = 1, g_blsize2(iproc)
      iposa(j,iproc) = ipos
      ipos           = ipos + g_blsize1(iproc)
    END DO
  END DO

  ! Make sure we've got the full array
  IF (ipos /= global_rows * global_row_len) THEN
    WRITE(cmessage, *) 'Addressing gone wrong in field copy'
    icode = 3
    CALL ereport(routinename, icode, cmessage)
  END IF

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( nproc, g_blsize2, igpos1d, global_row_len, g_blsize1,    &
!$OMP         iposa, global_field_copy, global_field )                 &
!$OMP PRIVATE( iproc, j, igpos, i, imax )
  DO iproc = 0, nproc-1
    DO j = 1, g_blsize2(iproc)
      igpos = igpos1d(iproc) + j*global_row_len
      imax = g_blsize1(iproc)
      DO i = 1, imax
        global_field_copy(iposa(j,iproc)+i) = global_field(igpos+i)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  DEALLOCATE( iposa )

END IF

isum = 0
DO iproc = 0, nproc-1
  displacements(iproc) = isum
  isum = isum + counts(iproc)
END DO ! nprocs

! 3.0 Do the exchange of data
CALL gc_get_communicator(my_comm, info)

! ipos is the location in the recv array of the first datapoint
ipos = halosize(2,halo_type) * local_row_len                         &
     + halosize(1,halo_type) + 1
CALL mpl_scatterv(global_field_copy, counts, displacements,          &
                  mpl_real4, local_field(ipos),                       &
                  1, recv_type, scatter_pe, my_comm, info)

#endif
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE scatter_field_mpl32
