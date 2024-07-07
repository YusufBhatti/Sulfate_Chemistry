! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Parallel cv_coarse_grid: Perform 2D data decomposition

MODULE crmstyle_decompose_grid_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_DECOMPOSE_GRID_MOD'

CONTAINS

SUBROUTINE crmstyle_decompose_grid(deco,global_row_len, global_n_rows, &
                           tot_levels,                                 &
                           model_type,                                 &
                           nproc_ew, nproc_ns, maxproc,                &
                           extended_halo_ew,                           &
                           extended_halo_ns,                           &
                           nstart_x, nstart_y,                         &
                           local_row_len, local_n_rows )


USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v, nfld_max
USE decomp_params, ONLY: decomp_standard_atmos
USE UM_ParCore, ONLY: mype, nproc_max
USE UM_Parparams, ONLY: Ndim_max, halo_type_no_halo, halo_type_single,  &
    halo_type_extended, NHalo_max, bc_static
USE crmstyle_cntl_mod, ONLY: l_endgame  

USE decomp_db, ONLY: allocate_decomposition, set_neighbour, &
    Decomp_DB_type, DecompDB

USE gcom_mod
USE UM_ParVars
USE model_domain_mod, ONLY: mt_bi_cyclic_lam

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Description:
! This routine performs a 2D decomposition - taking the input X
! (global_row_len) and input Y (global_n_rows) data sizes and decomposing
! across nproc_ew processors in the X direction and nproc_ns processors
! in the Y direction.
! ( NOTE this routine is based on the UM decompose_full routine.)
!
! Method:
! The local data sizes are calculated and stored in the module block
! DECOMPDB.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! Subroutine Arguments:

INTEGER, INTENT(IN) ::  &
  deco                  & ! decomposition identifier
 ,global_row_len        & ! number of E-W points of entire model
 ,global_n_rows         & ! number of P rows of entire model
 ,tot_levels            & ! total number of levels
 ,model_type            & ! type (Global,LAM etc) of model
 ,nproc_ew              & ! number of processors East-West
 ,nproc_ns              & ! number of processors North-South
 ,maxproc               & ! maximum number of processors/nodes
 ,extended_halo_ew      & ! size of extended EW halo
 ,extended_halo_ns        ! size of extended NS halo

INTEGER, INTENT(IN) ::  &
  nstart_x              & ! start point on old grid in x direction
 ,nstart_y                ! start point on old grid in y direction


INTEGER, INTENT(IN) ::  &
  local_row_len         & ! Number of E-W points on this processor
 ,local_n_rows            ! Number of rows on this processor


! Local variables
INTEGER ::                                                        &
  iproc                                                           &
, iproc_x                                                         &
, iproc_y                                                         &
, ifld                                                            &
, ihalo                                                           &
, iidim                                                           &
, ipt                                                             &
, start                                                           &
, isize                                                           &
, info                                                            &
, icode                                                           &
, smallhalosize                                                   &
, in_atm_decomp                                                   &
, my_comm       ! current communicator


INTEGER  ::                                                       &
  size_x(0:nproc_ew-1)                                            &
, size_y(0:nproc_ns-1)                                            &
, start_x(0:nproc_ew-1)                                           &
, start_y(0:nproc_ns-1)


! Error reporting
INTEGER  ::     errorstatus ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error message
CHARACTER(LEN=*), PARAMETER :: RoutineName='CRMSTYLE_DECOMPOSE_GRID'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype  ==  0) THEN

  ! Checking valid decomposition of grid

  IF (nproc_ew*nproc_ns  >   maxproc) THEN
    errorstatus=1
    WRITE(cmessage,'(a,i3,a,i3,a,i3,a,i3,a)')                     &
    'Cannot run with decomposition ',nproc_ew,'x',nproc_ns,       &
    '(',nproc_ew*nproc_ns,') processors. Maxproc is ',maxproc,    &
    ' processors.'

    CALL ereport(routinename,errorstatus,cmessage)
  END IF

END IF

CALL Allocate_Decomposition(deco)

! ------------------------------------------------------------------
DecompDB(deco)%sb_model_type = model_type

! I only need single point halos
smallHaloSize=1
DecompDB(deco)%halosize(1,halo_type_single)   = smallHaloSize
DecompDB(deco)%halosize(2,halo_type_single)   = smallHaloSize
DecompDB(deco)%halosize(3,halo_type_single)   = 0

! I only need single point halos so setting extended halo to the same
DecompDB(deco)%halosize(1,halo_type_extended) = 1
DecompDB(deco)%halosize(2,halo_type_extended) = 1
DecompDB(deco)%halosize(3,halo_type_extended) = 0

DecompDB(deco)%halosize(1,halo_type_no_halo) = 0
DecompDB(deco)%halosize(2,halo_type_no_halo) = 0
DecompDB(deco)%halosize(3,halo_type_no_halo) = 0

! ------------------------------------------------------------------
! 1.0 Set up global data size - ok for new dynamics may need to
!     check for v at poles
! ------------------------------------------------------------------

DecompDB(deco)%glsize(1,fld_type_p) = global_row_len
DecompDB(deco)%glsize(2,fld_type_p) = global_n_rows
DecompDB(deco)%glsize(3,fld_type_p) = tot_levels

DecompDB(deco)%glsize(1,fld_type_u) = global_row_len
DecompDB(deco)%glsize(2,fld_type_u) = global_n_rows
DecompDB(deco)%glsize(3,fld_type_u) = tot_levels

DecompDB(deco)%glsize(1,fld_type_v) = global_row_len

IF (model_type /= mt_bi_cyclic_lam) THEN
  IF (l_endgame) THEN
    ! ENDGAME : assume that global_n_rows = number of u/p rows
    DecompDB(deco)%glsize(2,fld_type_v) = global_n_rows+1
  ELSE
    DecompDB(deco)%glsize(2,fld_type_v) = global_n_rows-1
  END IF
ELSE

  DecompDB(deco)%glsize(2,fld_type_v) = global_n_rows

END IF
DecompDB(deco)%glsize(3,fld_type_v) = tot_levels

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------

! select processors to use for the data decomposition
DecompDB(deco)%nproc = nproc_ew*nproc_ns
DecompDB(deco)%first_comp_pe = 0
DecompDB(deco)%last_comp_pe = DecompDB(deco)%nproc-1

!     Set the grid size

DecompDB(deco)%gridsize(1) = nproc_ew
DecompDB(deco)%gridsize(2) = nproc_ns
DecompDB(deco)%gridsize(3) = 1

! Work out the decomposition in the East-West direction.
! All processor get same sized area
! Note (not currently adding any halo?)

start= nstart_x
isize= local_row_len

! all processors get the same size grid areas

DO iproc=1,nproc_ew
  start_x(iproc-1)= start
  size_x(iproc-1) = isize

  start=start+size_x(iproc-1)  ! compute start for next processor

END DO

! Work out the decomposition in the North-South direction.
! All processor get same sized area

start= nstart_y
isize = local_n_rows

DO iproc=1,nproc_ns
  start_y(iproc-1) =start
  size_y(iproc-1)=isize
  start=start+size_y(iproc-1)
END DO


! Set the local data shape and offsets of each processor

DO iproc_y=0,nproc_ns-1

  DO iproc_x=0,nproc_ew-1

    iproc=DecompDB(deco)%first_comp_pe+iproc_x+(iproc_y*nproc_ew)

    ! Set the position in the Logical Processor Grid

    DecompDB(deco)%g_gridpos(1,iproc)=iproc_x
    DecompDB(deco)%g_gridpos(2,iproc)=iproc_y
    DecompDB(deco)%g_gridpos(3,iproc)=0

    ! Set the number of local datapoints (blsize) on the processor

    !  Fields on P grid:
    DecompDB(deco)%g_blsize(1,fld_type_p,iproc)= size_x(iproc_x)
    DecompDB(deco)%g_blsize(2,fld_type_p,iproc)= size_y(iproc_y)
    DecompDB(deco)%g_blsize(3,fld_type_p,iproc)= tot_levels

    ! Fields on U grid:
    DecompDB(deco)%g_blsize(1,fld_type_u,iproc)= size_x(iproc_x)
    DecompDB(deco)%g_blsize(2,fld_type_u,iproc)= size_y(iproc_y)
    DecompDB(deco)%g_blsize(3,fld_type_u,iproc)= tot_levels

    ! Fields on V grid:
    DecompDB(deco)%g_blsize(1,fld_type_v,iproc)= size_x(iproc_x)
    DecompDB(deco)%g_blsize(2,fld_type_v,iproc)= size_y(iproc_y)
    DecompDB(deco)%g_blsize(3,fld_type_v,iproc)= tot_levels

    ! Set the number of points including the halos on the processor

    DO ihalo=1,3
      DO ifld=1,nfld_max
        DO iidim=1,ndim_max
          DecompDB(deco)%g_lasize(iidim,ifld,ihalo,iproc)=              &
                  DecompDB(deco)%g_blsize(iidim,ifld,iproc)+            &
                  2*DecompDB(deco)%halosize(iidim,ihalo)
        END DO  ! iidim
      END DO  ! ifld
    END DO  ! ihalo

    ! Set the starting point in the global domain

    DecompDB(deco)%g_datastart(1,iproc)= start_x(iproc_x)
    DecompDB(deco)%g_datastart(2,iproc)= start_y(iproc_y)
    DecompDB(deco)%g_datastart(3,iproc)=1

    DecompDB(deco)%g_datastart_f(1,fld_type_p,iproc)=start_x(iproc_x)
    DecompDB(deco)%g_datastart_f(2,fld_type_p,iproc)=start_y(iproc_y)
    DecompDB(deco)%g_datastart_f(3,fld_type_p,iproc)=1

    DecompDB(deco)%g_datastart_f(1,fld_type_u,iproc)=start_x(iproc_x)
    DecompDB(deco)%g_datastart_f(2,fld_type_u,iproc)=start_y(iproc_y)
    DecompDB(deco)%g_datastart_f(3,fld_type_u,iproc)=1

    DecompDB(deco)%g_datastart_f(1,fld_type_v,iproc)=start_x(iproc_x)
    DecompDB(deco)%g_datastart_f(2,fld_type_v,iproc)=start_y(iproc_y)
    DecompDB(deco)%g_datastart_f(3,fld_type_v,iproc)=1

  END DO ! iproc_x
END DO ! iproc_y

! Set up the pe_index_EW array - for each point along a global row
! it indicates the PE index (along the processor row) which
! contains that point

DO iproc_x=0,nproc_ew-1
  DO ipt=DecompDB(deco)%g_datastart(1,iproc_x),  &
          DecompDB(deco)%g_datastart(1,iproc_x)+size_x(iproc_x)
    DecompDB(deco)%g_pe_index_ew(ipt)=iproc_x
  END DO
END DO

! And fill in the halos at either end

DO ipt=1-extended_halo_ew,0
  DecompDB(deco)%g_pe_index_ew(ipt)=0
END DO

DO ipt=global_row_len+1,global_row_len+extended_halo_ew
  DecompDB(deco)%g_pe_index_ew(ipt)=nproc_ew-1
END DO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

DO iproc_y=0,nproc_ns-1
  DO ipt=DecompDB(deco)%g_datastart(2,iproc_y*nproc_ew),            &
          DecompDB(deco)%g_datastart(2,iproc_y*nproc_ew)+ size_y(iproc_y)
    DecompDB(deco)%g_pe_index_ns(ipt)=iproc_y
  END DO
END DO

! And fill in the halos at either end

DO ipt=1-extended_halo_ns,0
  DecompDB(deco)%g_pe_index_ns(ipt)=0
END DO

DO ipt=global_n_rows+1,global_n_rows+extended_halo_ns
  DecompDB(deco)%g_pe_index_ns(ipt)=nproc_ns-1
END DO

! ------------------------------------------------------------------
! 3.0 Set boundary conditions - none so not important
! ------------------------------------------------------------------
DecompDB(deco)%bound(1) = bc_static
DecompDB(deco)%bound(2) = bc_static
DecompDB(deco)%bound(3) = bc_static
CALL set_neighbour(deco)

! ------------------------------------------------------------------
! 4.0 Return the new data sizes
! ------------------------------------------------------------------

! Set up the GCOM groups (effectively communicators in MPI):
CALL gc_get_communicator(my_comm, info)

! 1) Group of all processors on my row

IF (DecompDB(deco)%gridsize(2)  ==  1 ) THEN         
  DecompDB(deco)%gc_proc_row_group=my_comm
ELSE
  CALL gcg_split(mype,nproc_max,                                    &
      DecompDB(deco)%g_gridpos(2,mype),                             &
      info,                                                         &
      DecompDB(deco)%gc_proc_row_group)
END IF

! 2) Group of all processors on my column

IF ( DecompDB(deco)%gridsize(1)  ==  1) THEN      
  DecompDB(deco)%gc_proc_col_group=my_comm
ELSE
  CALL gcg_split(mype,nproc_max,                                  &
      DecompDB(deco)%g_gridpos(1,mype),                           &
      info,                                                       &
      DecompDB(deco)%gc_proc_col_group)
END IF

! 3) Group of all processors in the atmosphere model
IF (DecompDB(deco)%nproc  ==  nproc_max) THEN    
  DecompDB(deco)%gc_all_proc_group=my_comm
ELSE
  IF ((mype  >=  DecompDB(deco)%first_comp_pe) .AND.       &
      (mype  <=  DecompDB(deco)%last_comp_pe) ) THEN       
    in_atm_decomp=1
  ELSE
    in_atm_decomp=0
  END IF

  CALL gcg_split(mype,nproc_max,in_atm_decomp,info,               &
          DecompDB(deco)%gc_all_proc_group)
END IF


! Set logical indicating this decomposition has been initialised
! and is now ready for use

DecompDB(deco)%set=.TRUE.

! And return the new horizontal dimensions
! Not needed
!local_row_len=DecompDB(deco)%g_blsize(1,fld_type_p,mype)
!local_n_rows=DecompDB(deco)%g_blsize(2,fld_type_p,mype)

! ------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE crmstyle_decompose_grid
END MODULE crmstyle_decompose_grid_mod
