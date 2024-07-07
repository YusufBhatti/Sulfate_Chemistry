! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Decomposition of grid for reconfiguration

MODULE Rcf_Decompose_Mod

!  Subroutine Rcf_Decompose - performs grid decomposition.
!
! Description:
! This routine performs a 2D decomposition - taking the global X and Y
! data sizes and decomposing across nproc_EW processors in the X
! direction and nproc_NS processors  in the Y direction.
!
! Method:
!   Based on UM4.5 code but now uses F90 data-structures.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DECOMPOSE_MOD'

CONTAINS

SUBROUTINE Rcf_Decompose( grid, nproc_EW, nproc_NS, decomp )

USE UM_ParVars      ! Most of this is used
USE Decomp_DB       ! Most of this is used
USE field_types
USE gcom_mod

USE decomp_params, ONLY: &
    decomp_rcf_input,    &
    decomp_rcf_output

USE Rcf_Lsm_Mod, ONLY: &
    Glob_LSM_Out,       &
    Glob_LSM_In,        &
    Local_LSM_Out,      &
    Local_LSM_In

USE Rcf_Set_Neighbour_Mod, ONLY: &
    Rcf_Set_Neighbour

USE Rcf_Grid_Type_Mod, ONLY:     &
    grid_type

USE Rcf_Headaddress_Mod, ONLY: &
    FH_GridStagger_A

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( grid_type ), INTENT(INOUT)  :: grid
INTEGER, INTENT(IN)               :: nproc_EW
INTEGER, INTENT(IN)               :: nproc_NS
INTEGER, INTENT(IN)               :: decomp


! Local variables
INTEGER                     :: iproc
INTEGER                     :: irest
INTEGER                     :: jrest
INTEGER                     :: info
INTEGER                     :: in_atm_decomp
INTEGER                     :: my_comm
INTEGER                     :: ihalo   ! looper
INTEGER                     :: ifld    ! looper
INTEGER                     :: i_dim   ! looper

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DECOMPOSE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL Allocate_Decomposition(decomp)

! ------------------------------------------------------------------
decompDB (decomp) % halosize (1,1) = 0
decompDB (decomp) % halosize (2,1) = 0
decompDB (decomp) % halosize (3,1) = 0

! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

! recon only has one fld type, so second index is 1.
decompDB (decomp) % glsizep (1) = grid % glob_p_row_length
decompDB (decomp) % glsizep (2) = grid % glob_p_rows
decompDB (decomp) % glsizep (3) = grid % model_levels

decompDB (decomp) % glsizeu (1) = grid % glob_u_row_length
decompDB (decomp) % glsizeu (2) = grid % glob_u_rows
decompDB (decomp) % glsizeu (3) = grid % model_levels

decompDB (decomp) % glsizev (1) = grid % glob_v_row_length
decompDB (decomp) % glsizev (2) = grid % glob_v_rows
decompDB (decomp) % glsizev (3) = grid % model_levels

decompDB (decomp) % glsizer (1) = grid % glob_r_row_length
decompDB (decomp) % glsizer (2) = grid % glob_r_rows
decompDB (decomp) % glsizer (3) = 1

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------


! select processors to use for the data decomposition
decompDB (decomp) % nproc          = nproc_EW * nproc_NS
decompDB (decomp) % first_comp_pe  = 0
decompDB (decomp) % last_comp_pe   =  decompDB (decomp) % nproc -1

!     Set the grid size

decompDB (decomp) % gridsize (1) = nproc_EW
decompDB (decomp) % gridsize (2) = nproc_NS
decompDB (decomp) % gridsize (3) = 1

! Calculate the local data shape of each processor.
DO iproc=decompDB (decomp) % first_comp_pe ,                  &
    decompDB (decomp) % last_comp_pe
  !       ! Loop over all processors in this decomposition

  decompDB (decomp) % g_gridpos (3,iproc) = 0

  decompDB (decomp) % g_gridpos (2,iproc) =                   &
      iproc / decompDB (decomp) % gridsize (1)

  decompDB (decomp) % g_gridpos (1,iproc) =                   &
      iproc - decompDB (decomp) % g_gridpos (2,iproc)*        &
      decompDB (decomp) % gridsize (1)

  !  Find the number of grid points (blsizep) in each domain and starting
  !  points in the global domain (datastart) We first try to divide
  !  the total number equally among the processors. The rest is
  !  distributed one by one to first processor in each direction.

  ! The X (East-West) direction:

  decompDB (decomp) % g_blsize (1,fld_type_p,iproc) =         &
      decompDB (decomp) % glsize (1,1) /                      &
      decompDB (decomp) % gridsize (1)
  irest = decompDB (decomp) % glsize (1,1)-                   &
      decompDB (decomp) % g_blsize (1,fld_type_p,iproc)*      &
      decompDB (decomp) % gridsize (1)
  decompDB (decomp) % g_datastart (1,iproc) =                 &
      decompDB (decomp) % g_gridpos (1,iproc)*                &
      decompDB (decomp) % g_blsize (1,fld_type_p,iproc) + 1

  IF (decompDB (decomp) % g_gridpos (1,iproc) <               &
      irest) THEN
    decompDB (decomp) % g_blsize (1,fld_type_p,iproc) =       &
        decompDB (decomp) % g_blsize (1,fld_type_p,iproc)+1
    decompDB (decomp) % g_datastart (1,iproc) =               &
        decompDB (decomp) % g_datastart (1,iproc) +           &
        decompDB (decomp) % g_gridpos (1,iproc)
  ELSE
    decompDB (decomp) % g_datastart (1,iproc) =               &
        decompDB (decomp) % g_datastart (1,iproc) +           &
        irest
  END IF

  ! X direction for river routing

  decompDB (decomp) % g_blsize (1,fld_type_r,iproc) =         &
      decompDB (decomp) % glsizer (1) /                       &
      decompDB (decomp) % gridsize (1)
  irest = decompDB (decomp) % glsizer (1) -                   &
      decompDB (decomp) % g_blsize (1,fld_type_r,iproc)*      &
      decompDB (decomp) % gridsize (1)
  decompDB (decomp) % g_datastartr (1,iproc) =                &
      decompDB (decomp) % g_gridpos (1,iproc) *               &
      decompDB (decomp) % g_blsize (1,fld_type_r,iproc) + 1

  IF (decompDB (decomp) % g_gridpos (1,iproc) <               &
      irest) THEN
    decompDB (decomp) % g_blsize (1,fld_type_r,iproc) =       &
        decompDB (decomp) % g_blsize (1,fld_type_r,iproc) + 1

    decompDB (decomp) % g_datastartr (1,iproc) =              &
        decompDB (decomp) % g_datastartr (1,iproc) +          &
        decompDB (decomp) % g_gridpos (1,iproc)
  ELSE
    decompDB (decomp) % g_datastartr (1,iproc) =              &
        decompDB (decomp) % g_datastartr (1,iproc) +          &
        irest
  END IF

  ! The Y (North-South) direction

  decompDB (decomp) % g_blsize (2,fld_type_p,iproc) =         &
      decompDB (decomp) % glsize (2,1) /                      &
      decompDB (decomp) % gridsize (2)
  jrest = decompDB (decomp) % glsize (2,1)-                   &
      decompDB (decomp) % g_blsize (2,fld_type_p,iproc)*      &
      decompDB (decomp) % gridsize (2)
  decompDB (decomp) % g_datastart (2,iproc) =                 &
      decompDB (decomp) % g_gridpos (2,iproc)*                &
      decompDB (decomp) % g_blsize (2,fld_type_p,iproc) + 1

  IF (decompDB (decomp) % g_gridpos (2,iproc) <               &
      jrest) THEN
    decompDB (decomp) % g_blsize (2,fld_type_p,iproc) =       &
        decompDB (decomp) % g_blsize (2,fld_type_p,iproc)+1
    decompDB (decomp) % g_datastart (2,iproc) =               &
        decompDB (decomp) % g_datastart (2,iproc) +           &
        decompDB (decomp) % g_gridpos (2,iproc)
  ELSE
    decompDB (decomp) % g_datastart (2,iproc) =               &
        decompDB (decomp) % g_datastart (2,iproc) +  jrest
  END IF

  ! Y decomp for rivers

  decompDB (decomp) % g_blsize (2,fld_type_r,iproc) =         &
      decompDB (decomp) % glsizer (2) /                       &
      decompDB (decomp) % gridsize (2)
  jrest = decompDB (decomp) % glsizer (2)-                    &
      decompDB (decomp) % g_blsize (2,fld_type_r,iproc)*      &
      decompDB (decomp) % gridsize (2)
  decompDB (decomp) % g_datastartr (2,iproc) =                &
      decompDB (decomp) % g_gridpos (2,iproc) *               &
      decompDB (decomp) % g_blsize (2,fld_type_r,iproc) + 1

  IF (decompDB (decomp) % g_gridpos (2,iproc) <               &
      jrest) THEN
    decompDB (decomp) % g_blsize (2,fld_type_r,iproc) =       &
        decompDB (decomp) % g_blsize (2,fld_type_r,iproc) + 1

    decompDB (decomp) % g_datastartr (2,iproc) =              &
        decompDB (decomp) % g_datastartr (2,iproc) +          &
        decompDB (decomp) % g_gridpos (2,iproc)
  ELSE
    decompDB (decomp) % g_datastartr (2,iproc) =              &
        decompDB (decomp) % g_datastartr (2,iproc) +  jrest
  END IF

  !  The Z (vertical) direction (no decomposition):

  decompDB (decomp) % g_datastart (3,iproc) = 1
  decompDB (decomp) % g_blsize (3,fld_type_p,iproc)   =  grid % model_levels
  decompDB (decomp) % g_blsize(3,fld_type_r,iproc)    = 1

  ! Setup the g_datastart_f arrays - i.e. field indexed
  ! p fields
  decompDB (decomp) % g_datastart_f(1,fld_type_p,iproc) =     &
               decompDB (decomp) % g_datastart(1,iproc)
  decompDB (decomp) % g_datastart_f(2,fld_type_p,iproc) =     &
               decompDB (decomp) % g_datastart(2,iproc)
  decompDB (decomp) % g_datastart_f(3,fld_type_p,iproc) =     &
               decompDB (decomp) % g_datastart(3,iproc)

  ! u fields
  decompDB (decomp) % g_datastart_f(1,fld_type_u,iproc) =     &
               decompDB (decomp) % g_datastart(1,iproc)
  decompDB (decomp) % g_datastart_f(2,fld_type_u,iproc) =     &
               decompDB (decomp) % g_datastart(2,iproc)
  decompDB (decomp) % g_datastart_f(3,fld_type_u,iproc) =     &
               decompDB (decomp) % g_datastart(3,iproc)

  ! v fields
  decompDB (decomp) % g_datastart_f(1,fld_type_v,iproc) =     &
               decompDB (decomp) % g_datastart(1,iproc)
  decompDB (decomp) % g_datastart_f(2,fld_type_v,iproc) =     &
               decompDB (decomp) % g_datastart(2,iproc)
  decompDB (decomp) % g_datastart_f(3,fld_type_v,iproc) =     &
               decompDB (decomp) % g_datastart(3,iproc)

  ! r fields
  decompDB (decomp) % g_datastart_f(1,fld_type_r,iproc) =     &
               decompDB (decomp) % g_datastartr(1,iproc)
  decompDB (decomp) % g_datastart_f(2,fld_type_r,iproc) =     &
               decompDB (decomp) % g_datastartr(2,iproc)
  decompDB (decomp) % g_datastart_f(3,fld_type_r,iproc) =     &
               decompDB (decomp) % g_datastartr(3,iproc)

  ! Now sort out the relevant u and v grid types.
  ! U works out the box as we choose to keep row length the
  ! same for both ND and EG, even if the data is not used
  ! in the last u (ND) and p/v (EG) col.


  decompDB (decomp) % g_blsize (1,fld_type_u,iproc) =         &
      decompDB (decomp) % g_blsize (1,fld_type_p,iproc)

  decompDB (decomp) % g_blsize (2,fld_type_u,iproc) =         &
      decompDB (decomp) % g_blsize (2,fld_type_p,iproc)

  decompDB (decomp) % g_blsize (3,fld_type_u,iproc) =         &
      decompDB (decomp) % g_blsize (3,fld_type_p,iproc)


  ! One less v row at the north in ND while one more for EG.

  decompDB (decomp) % g_blsize (1,fld_type_v,iproc) =         &
      decompDB (decomp) % g_blsize (1,fld_type_p,iproc)

  IF ( (decompDB (decomp) % g_gridpos (2,iproc) ==            &
      decompDB (decomp) % gridsize (2) - 1) .AND.             &
      grid % glob_v_rows ==  grid % glob_p_rows - 1 ) THEN
    decompDB (decomp) % g_blsize (2,fld_type_v,iproc) =       &
        decompDB (decomp) % g_blsize (2,fld_type_p,iproc) - 1
    ! For ENDGAME there is one more v row
  ELSE IF ( (decompDB (decomp) % g_gridpos (2,iproc) ==       &
      decompDB (decomp) % gridsize (2) - 1) .AND.             &
      grid % glob_v_rows ==  grid % glob_p_rows + 1 ) THEN
    decompDB (decomp) % g_blsize (2,fld_type_v,iproc) =       &
        decompDB (decomp) % g_blsize (2,fld_type_p,iproc) + 1

  ELSE
    decompDB (decomp) % g_blsize (2,fld_type_v,iproc) =       &
        decompDB (decomp) % g_blsize (2,fld_type_p,iproc)
  END IF

  decompDB (decomp) % g_blsize (3,fld_type_v,iproc) =         &
      decompDB (decomp) % g_blsize (3,fld_type_p,iproc)


   ! Set the number of points on the processor including halos
   DO ihalo=1,nhalo_max
     DO ifld=1,nfld_max
       DO i_dim=1,ndim_max
         DecompDB(decomp)%g_lasize(i_dim,ifld,ihalo,iproc)=          &
             DecompDB(decomp)%g_blsize(i_dim,ifld,iproc)+            &
             2*DecompDB(decomp)%halosize(i_dim,ihalo)
       END DO  ! idim
     END DO  ! ifld
   END DO  ! ihalo

END DO ! loop over processors


! ------------------------------------------------------------------
! 3.0 Set boundary conditions
! ------------------------------------------------------------------

! if global or wrapping LAM
IF ( grid % glob_p_row_length ==  grid % glob_u_row_length ) THEN
  decompDB (decomp) % bound (1) = bc_cyclic ! Cyclic East-West bdy
ELSE
  decompDB (decomp) % bound (1) = bc_static ! No East-West wrap around
END IF

IF ( grid % glob_p_rows == grid % glob_v_rows .AND.          &
     grid % grid_stagger /= FH_GridStagger_A  ) THEN
  decompDB (decomp) % bound (2) = bc_cyclic ! cyclic N-S
ELSE
  decompDB (decomp) % bound (2) = bc_static !No North-South wrap around
END IF
decompDB (decomp) % bound (3) = bc_static ! No vertical wrap around

CALL Rcf_Set_Neighbour (decomp)

! ------------------------------------------------------------------
! 4.0 Return the new data sizes and exit subroutine
! ------------------------------------------------------------------

! Set up the GCOM groups (effectively communicators in MPI):
CALL gc_get_communicator(my_comm, info)

! 1) Group of all processors on my row

IF ( decompDB (decomp) % gridsize (2) == 1) THEN
  decompDB (decomp) % gc_proc_row_group = my_comm
ELSE
  CALL gcg_split(mype,nproc_max,                              &
      decompDB (decomp) % g_gridpos (2,mype), info,           &
      decompDB (decomp) % gc_proc_row_group )
END IF

! 2) Group of all processors on my column

IF ( decompDB (decomp) % gridsize (1) == 1) THEN  
  decompDB (decomp) % gc_proc_col_group = my_comm
ELSE
  CALL gcg_split(mype,nproc_max,                              &
      decompDB (decomp) % g_gridpos (1,mype), info,           &
      decompDB (decomp) % gc_proc_col_group )
END IF

! 3) Group of all processors in the atmosphere model
IF (decompDB (decomp) % nproc  == nproc_max) THEN
  decompDB (decomp) % gc_all_proc_group = my_comm
ELSE
  IF ((mype  >=  decompDB (decomp) % first_comp_pe ) .AND.    &
      (mype  <=  decompDB (decomp) % last_comp_pe ) )         &
      THEN
    in_atm_decomp=1
  ELSE
    in_atm_decomp=0
  END IF

  CALL gcg_split(mype,nproc_max,in_atm_decomp,info,           &
      decompDB (decomp) % gc_all_proc_group )
END IF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

decompDB (decomp) % set =.TRUE.

! And return the new horizontal dimensions

grid % loc_p_row_length =                                     &
    decompDB (decomp) % g_blsize (1,fld_type_p,mype)
grid % loc_p_rows       =                                     &
    decompDB (decomp) % g_blsize (2,fld_type_p,mype)
grid % loc_p_field      =                                     &
    grid % loc_p_row_length * grid % loc_p_rows
grid % loc_u_row_length =                                     &
    decompDB (decomp) % g_blsize (1,fld_type_u,mype)
grid % loc_u_rows       =                                     &
    decompDB (decomp) % g_blsize (2,fld_type_u,mype)
grid % loc_u_field      =                                     &
    grid % loc_u_row_length * grid % loc_u_rows
grid % loc_v_row_length =                                     &
    decompDB (decomp) % g_blsize (1,fld_type_v,mype)
grid % loc_v_rows       =                                     &
    decompDB (decomp) % g_blsize (2,fld_type_v,mype)
grid % loc_v_field      =                                     &
    grid % loc_v_row_length * grid % loc_v_rows
grid % loc_r_row_length =                                     &
    decompDB (decomp) % g_blsize (1,fld_type_r,mype)
grid % loc_r_rows       =                                     &
    decompDB (decomp) % g_blsize (2,fld_type_r,mype)
grid % loc_r_field      =                                     &
    grid % loc_r_row_length * grid % loc_r_rows

! Finally allocate space for land-sea-mask
IF (decomp == decomp_rcf_output) THEN
  ALLOCATE( glob_lsm_out                                      &
      (grid % glob_p_row_length * grid % glob_p_rows))
  ALLOCATE( local_lsm_out                                     &
      ( grid % loc_p_row_length * grid % loc_p_rows))

ELSE IF (decomp == decomp_rcf_input) THEN
  ALLOCATE( glob_lsm_in                                       &
      (grid % glob_p_row_length * grid % glob_p_rows))
  ALLOCATE( local_lsm_in                                      &
      ( grid % loc_p_row_length * grid % loc_p_rows))
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Decompose
END MODULE Rcf_Decompose_Mod
