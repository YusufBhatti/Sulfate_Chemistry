! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP.
!
! Provides variables that describe the current decomposition, and a routine to
! switch between decompositions


MODULE UM_ParVars

USE Field_Types, ONLY: &
  nfld_max, fld_type_p, fld_type_u, fld_type_v, fld_type_r

USE UM_ParParams, ONLY: &
  pnorth, psouth, peast, pwest, Ndim_max, Nhalo_max, halo_type_extended, &
  halo_type_single

USE decomp_params, ONLY: decomp_unset, max_decomps

USE UM_ParCore, ONLY: mype, nproc

USE ereport_mod, ONLY: ereport

USE Atmos_Max_Sizes, ONLY: Max_Halo_Size, row_length_max, rows_max

#if !defined(LFRIC)
USE decomp_db, ONLY: decompDB
#endif

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

PRIVATE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! Parameters and globals required by the UM and other execs

INTEGER, PUBLIC :: current_decomp_type=decomp_unset ! current decomposition type
INTEGER, PUBLIC :: Offx                   ! standard halo size in East-West
INTEGER, PUBLIC :: Offy                   ! standard halo size in North-South
INTEGER, PUBLIC :: halo_i                 ! extended halo size in East-West
INTEGER, PUBLIC :: halo_j                 ! extended halo size in North-South

INTEGER, PUBLIC :: gc_proc_row_group      ! GID for procs along a proc row
INTEGER, PUBLIC :: gc_proc_col_group      ! GID for procs along a proc col
INTEGER, PUBLIC :: gc_all_proc_group      ! GID for all procs

! logicals indicating if a processor is at the edge of the LPG
LOGICAL, PUBLIC :: at_extremity(4)

! array with the tids of the four neighbours in the horizontal plane
INTEGER, PUBLIC :: neighbour(4)

! array with the tids of the eight neighbours in the horizontal plane,
! i.e. include the corners.
INTEGER, PUBLIC :: full_neighbour(8)

! position of personal data in global data (in terms of standard
! Fortran array notation
INTEGER, PUBLIC :: datastart(Ndim_max)

INTEGER, PUBLIC :: first_comp_pe       ! top left pe in LPG
INTEGER, PUBLIC :: last_comp_pe        ! bottom right pe in LPG
INTEGER, PUBLIC :: halosize(Ndim_max,NHalo_max)        ! available halo sizes
INTEGER, TARGET, PUBLIC :: glsize(Ndim_max,Nfld_max)   ! global data size
INTEGER, TARGET, PUBLIC :: blsize(Ndim_max,Nfld_max)   ! personal data size
INTEGER, PUBLIC :: lasize(Ndim_max,Nfld_max,NHalo_max) ! local data size

! Generalised version of datastart for *all* fieldtypes
INTEGER, PUBLIC :: datastart_f(Ndim_max,Nfld_max)

! size of the LPG in each dimension
INTEGER, PUBLIC :: gridsize(Ndim_max)

! position of this process in the LPG 0,1,2,...,nproc_x-1 etc.
INTEGER, PUBLIC :: gridpos(Ndim_max)

! domain type
INTEGER, PUBLIC :: sb_model_type

! type of boundary (cyclic or static) in each direction
INTEGER, PUBLIC :: bound(Ndim_max)

INTEGER, PUBLIC :: datastartr(Ndim_max)

! Which processor column a given point is in: 0 -> nproc_x-1
INTEGER, PUBLIC :: g_pe_index_EW(1-Max_Halo_Size:row_length_max+Max_Halo_Size)

! Which processor row a given point is in: 0 -> nproc_y-1
INTEGER, PUBLIC :: g_pe_index_NS(1-Max_Halo_Size:rows_max+Max_Halo_Size)

INTEGER, PUBLIC :: nproc_x=1    ! number of processors in x-direction
INTEGER, PUBLIC :: nproc_y=1    ! number of processors in y-direction


! Short cut names that are used by some rcf code.
INTEGER, POINTER, PUBLIC :: glsizep(:)=>NULL() ! global u data size
INTEGER, POINTER, PUBLIC :: glsizeu(:)=>NULL() ! global u data size
INTEGER, POINTER, PUBLIC :: glsizev(:)=>NULL() ! global v data size
INTEGER, POINTER, PUBLIC :: glsizer(:)=>NULL() ! global river-routing data size
INTEGER, POINTER, PUBLIC :: blsizep(:)=>NULL() ! personal p data area
INTEGER, POINTER, PUBLIC :: blsizeu(:)=>NULL() ! personal u data area
INTEGER, POINTER, PUBLIC :: blsizev(:)=>NULL() ! personal v data area
INTEGER, POINTER, PUBLIC :: blsizer(:)=>NULL() ! personal river-routing data
                                               !   area

LOGICAL, PUBLIC ::       &
    atSouth,             &! process at the bottom of the LPG
    atNorth,             &! process at the top of the LPG
    atWest,              &! process at the left of the LPG
    atEast                ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UM_PARVARS'

#if !defined(LFRIC)

! Public functions and subroutines:
PUBLIC :: g_datastart,              &
          g_lasize,                 &
          g_blsize,                 &
          g_datastartr,             &
          g_datastart_f,            &
          g_gridpos,                &
          g_at_extremity,           &
          g_blsizep,                &
          g_blsizeu,                &
          g_blsizev,                &
          g_blsizer,                &
          change_decomposition

CONTAINS

INTEGER FUNCTION g_datastart(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_datastart=decompDB(current_decomp_type)%g_datastart(dim,proc)
END FUNCTION g_datastart

INTEGER FUNCTION g_lasize(dim,fld,halo,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
INTEGER, INTENT(IN) :: fld
INTEGER, INTENT(IN) :: halo
g_lasize= &
    decompDb(current_decomp_type)%g_lasize(dim,fld,halo,proc)
END FUNCTION g_lasize

INTEGER FUNCTION g_blsize(dim,fld,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
INTEGER, INTENT(IN) :: fld
g_blsize= &
    decompDb(current_decomp_type)%g_blsize(dim,fld,proc)
END FUNCTION g_blsize

INTEGER FUNCTION g_datastartr(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_datastartr= &
    decompDb(current_decomp_type)%g_datastartr(dim,proc)
END FUNCTION g_datastartr

INTEGER FUNCTION g_datastart_f(dim,fld,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: fld
INTEGER, INTENT(IN) :: proc
g_datastart_f= &
    decompDb(current_decomp_type)%g_datastart_f(dim,fld,proc)
END FUNCTION g_datastart_f

INTEGER FUNCTION g_gridpos(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_gridpos= &
    decompDb(current_decomp_type)%g_gridpos(dim,proc)
END FUNCTION g_gridpos

! at_extremity for each processor
LOGICAL FUNCTION g_at_extremity(dir,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dir
INTEGER, INTENT(IN) :: proc
INTEGER             :: errorstatus !
CHARACTER(LEN=errormessagelength)  :: cmessage    ! Error message

SELECT CASE(dir)

CASE (pnorth)
  g_at_extremity=(g_gridpos(2,proc)  ==  (gridsize(2)-1))
CASE (psouth)
  g_at_extremity=(g_gridpos(2,proc)  ==  0)
CASE (peast)
  g_at_extremity=(g_gridpos(1,proc)  ==  (gridsize(1)-1))
CASE (pwest)
  g_at_extremity=(g_gridpos(1,proc)  ==  0)
CASE DEFAULT
  errorstatus=1
  WRITE(cmessage,'(A,I3,A)')                     &
      'Request for extremity in direction ',dir,' is invalid'
  CALL ereport('UM_Parvars:g_at_extrimity',errorstatus,cmessage)
END SELECT
END FUNCTION g_at_extremity

INTEGER FUNCTION g_blsizep(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_blsizep= &
    decompDb(current_decomp_type)%g_blsize(dim,fld_type_p,proc)
END FUNCTION g_blsizep

INTEGER FUNCTION g_blsizeu(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_blsizeu= &
    decompDb(current_decomp_type)%g_blsize(dim,fld_type_u,proc)
END FUNCTION g_blsizeu

INTEGER FUNCTION g_blsizev(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_blsizev= &
    decompDb(current_decomp_type)%g_blsize(dim,fld_type_v,proc)
END FUNCTION g_blsizev

INTEGER FUNCTION g_blsizer(dim,proc)
IMPLICIT NONE
INTEGER, INTENT(IN) :: dim
INTEGER, INTENT(IN) :: proc
g_blsizer= &
    decompDb(current_decomp_type)%g_blsize(dim,fld_type_r,proc)
END FUNCTION g_blsizer

SUBROUTINE change_decomposition(decomp,icode)

! Method:
! If decomp is already the current decomposition, exit and do
! nothing. If decomposition decomp has not been initialised,
! print a message and exit. Otherwise, copy the
! information from the decompDB arrays in the DECOMP_DB module into
! here.

#if defined (RECON)
USE Rcf_Lsm_Mod
USE decomp_params, ONLY: decomp_rcf_input, decomp_rcf_output
#endif
IMPLICIT NONE

INTEGER, INTENT(IN) ::  decomp ! new decomposition to use
INTEGER, OPTIONAL   ::  icode  ! unused, for legacy reasons only

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
INTEGER                       :: errVal
INTEGER                       :: idm
INTEGER                       :: ifld
INTEGER                       :: ihalo
INTEGER                       :: ipt
INTEGER                       :: iside
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHANGE_DECOMPOSITION'
CHARACTER(LEN=errormessagelength)            :: cmessage ! Error message

! --------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
! Initial Checks:
!

IF (((decomp  >   max_decomps) .OR. (decomp  <   1)) .AND. &
    decomp  /=  decomp_unset) THEN
  WRITE(cmessage,'(A,I3)') &
      'Error: Cannot change to out of range decomposition ',decomp
  errVal=10
  CALL ereport("um_parvars:change_decomposition",errVal, &
      cmessage)
END IF

! Check to see if setting decomposition to unset
IF (decomp  ==  decomp_unset) THEN
  current_decomp_type = decomp_unset
END IF

! Check if this decomposition has been initialised

IF ( decomp /= decomp_unset ) THEN
  IF ( .NOT. decompDB(decomp)%set ) THEN
    WRITE(cmessage,'(A,I3)')                                                   &
        'Attempt to select uninitialised decomposition ', decomp
    errVal=10
    CALL ereport("um_parvars:change_decomposition", errVal, cmessage)
  END IF
END IF

! Check if this is already the current decomposition
IF (decomp  /=  current_decomp_type .AND. &
    decomp  /=  decomp_unset) THEN

  ! Set the current decomp early, we implicitly use it later on.
  current_decomp_type=decomp

  ! Set short aliases
  glsizep=>glsize(1:ndim_max,fld_type_p)
  glsizeu=>glsize(1:ndim_max,fld_type_u)
  glsizev=>glsize(1:ndim_max,fld_type_v)
  glsizer=>glsize(1:ndim_max,fld_type_r)
  blsizep=>blsize(1:ndim_max,fld_type_p)
  blsizeu=>blsize(1:ndim_max,fld_type_u)
  blsizev=>blsize(1:ndim_max,fld_type_v)
  blsizer=>blsize(1:ndim_max,fld_type_r)

  !
  ! Copy decomp_db vars into um_parvars
  !

  first_comp_pe   = decompDB(decomp)%first_comp_pe
  last_comp_pe    = decompDB(decomp)%last_comp_pe
  nproc           = decompDB(decomp)%nproc
  nproc_x         = decompDB(decomp)%gridsize(1)
  nproc_y         = decompDB(decomp)%gridsize(2)
  sb_model_type   = decompDB(decomp)%sb_model_type
  offx            = decompDB(decomp)%halosize(1,halo_type_single)
  offy            = decompDB(decomp)%halosize(2,halo_type_single)
  halo_i          = decompDB(decomp)%halosize(1,halo_type_extended)
  halo_j          = decompDB(decomp)%halosize(2,halo_type_extended)

  gc_proc_row_group = decompDB(decomp)%gc_proc_row_group
  gc_proc_col_group = decompDB(decomp)%gc_proc_col_group
  gc_all_proc_group = decompDB(decomp)%gc_all_proc_group

  halosize(:,:)     = decompDB(decomp)%halosize(:,:)
  neighbour(:)      = decompDB(decomp)%neighbour(:)
  full_neighbour(:) = decompDB(decomp)%full_neighbour(:)
  bound(:)          = decompDB(decomp)%bound(:)
  gridsize(:)       = decompDB(decomp)%gridsize(:)
  glsize(:,:)       = decompDB(decomp)%glsize(:,:)

  DO idm=1,ndim_max
    DO ifld=1,nfld_max
      DO ihalo=1,nhalo_max
        lasize(idm,ifld,ihalo) = g_lasize(idm,ifld,ihalo,mype)
      END DO
      blsize(idm,ifld)      = g_blsize(idm,ifld,mype)
      datastart_f(idm,ifld) = g_datastart_f(idm,ifld,mype)
    END DO
    datastart(idm)  = g_datastart(idm,mype)
    datastartr(idm) = g_datastartr(idm,mype)
    gridpos(idm)    = g_gridpos(idm,mype)
  END DO

  DO iside=1,4
    at_extremity(iside)=g_at_extremity(iside,mype)
  END DO

  DO ipt=1-decompDB(decomp)%halosize(1,halo_type_extended),     &
      decompDB(decomp)%glsize(1,fld_type_p)+                    &
      decompDB(decomp)%halosize(1,halo_type_extended)

    g_pe_index_ew(ipt)=decompDB(decomp)%g_pe_index_ew(ipt)

  END DO

  DO ipt=1-decompDB(decomp)%halosize(2,halo_type_extended),     &
      decompDB(decomp)%glsize(2,fld_type_p)+                    &
      decompDB(decomp)%halosize(2,halo_type_extended)

    g_pe_index_ns(ipt)=decompDB(decomp)%g_pe_index_ns(ipt)

  END DO

  atSouth = ( gridpos(2) == 0)
  atNorth = ( gridpos(2) == (gridsize(2)-1))
  atEast  = ( gridpos(1) == (gridsize(1)-1))
  atWest  = ( gridpos(1) == 0)

#if defined(RECON)
  ! Additional info for reconfiguration - for LandSeaMask
  IF (decomp == decomp_rcf_input ) THEN
    glob_atmos_landmask       => glob_lsm_in
    local_atmos_landmask      => local_lsm_in
    local_land_field          => local_land_in
    glob_land_field           => glob_land_in
  ELSE IF (decomp == decomp_rcf_output) THEN
    glob_atmos_landmask       => glob_lsm_out
    local_atmos_landmask      => local_lsm_out
    local_land_field          => local_land_out
    glob_land_field           => glob_land_out
  END IF
#endif
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE change_decomposition
#endif

END MODULE UM_ParVars
