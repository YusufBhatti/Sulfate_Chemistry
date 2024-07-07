! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

MODULE Rcf_Read_Multi_Mod

USE lookup_addresses

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READ_MULTI_MOD'

CONTAINS

! Subroutine Rcf_Read_Multi - parallel interface to buffin
!
! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  reconfiguration. It is used where each process must read in a
!  local section of a global field.
!
! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.

SUBROUTINE Rcf_Read_Multi(  nft, d1, isize, expand_size, len_io,  &
    local_len, iostatus, lookup, fixhd12,                         &
    Stash_Record )

USE io, ONLY:            &
    Buffin,               &
    set_unit_bcast_flag,  &
    clear_unit_bcast_flag

USE UM_ParVars, ONLY: &
    glsizep,          &
    glsize,           &
    g_blsize,         &
    gc_all_proc_group

USE UM_ParCore, ONLY: &
    mype

USE UM_ParParams, ONLY: &
    halo_type_single

USE Rcf_Global_To_Local_Mod, ONLY: &
    Rcf_Get_Fld_Type

USE wgdos_packing_mod, ONLY: wgdos_expand_field
USE Ereport_Mod, ONLY: &
    Ereport

USE umPrintMgr, ONLY:          &
    umPrint,                    &
    umMessage,                  &
    PrintStatus,                &
    PrStatus_Normal

USE Rcf_Average_Polar_Mod, ONLY: &
    Rcf_Average_Polar

USE Rcf_General_Scatter_Field_Mod, ONLY: &
    Rcf_General_Scatter_Field

USE Rcf_Ppx_Info_Mod, ONLY: &
    STM_Record_Type

USE cppxref_mod, ONLY:  &
    ppx_atm_compressed, &
    ppx_atm_tzonal,     &
    ppx_atm_uzonal

USE lookup_addresses

USE missing_data_mod, ONLY: rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE packing_codes_mod, ONLY: PC_WGDOS_Packing, PC_Cray32_Packing

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN)   :: nft        ! Fortain unit number
INTEGER, INTENT(IN)   :: isize      ! no. of words to read (glob)
INTEGER, INTENT(IN)   :: expand_size! size of expanded field
INTEGER, INTENT(IN)   :: Fixhd12    ! 12th element of fixed header
INTEGER, INTENT(IN)   :: Lookup(64) ! lookup table
INTEGER, INTENT(OUT)  :: Len_IO     ! no. of words read in (glob)
INTEGER, INTENT(OUT)  :: Local_Len  ! no. of local field words

REAL,    INTENT(OUT)  :: iostatus   ! Return Code
REAL,    INTENT(OUT)  :: d1(*)      ! Array to read data into

TYPE( STM_Record_Type ), INTENT(IN) :: Stash_Record

! Local variables
INTEGER      :: info            ! GCOM return code
INTEGER      :: dimx            ! x dimension returned from wgdos_expand_field
INTEGER      :: dimy            ! y dimension returned from wgdos_expand_field
INTEGER      :: idum            ! dummy for wgdos_expand_field
INTEGER      :: ErrorStatus     ! Error code returned from wgdos_expand_field
INTEGER      :: local_x         ! local X dimension
INTEGER      :: local_y         ! local Y dimension
INTEGER      :: global_x        ! global X dimension
INTEGER      :: global_y        ! global Y dimension
INTEGER      :: fld_type        ! p/u/v field
LOGICAL      :: averaged        ! Have polar rows been averaged?
LOGICAL      :: expanded        ! Have fields been expanded?
LOGICAL,PARAMETER :: global = .TRUE.
REAL         :: buf( expand_size ) ! Buffer for reading data into
INTEGER      :: int_buf( isize)    ! INTEGER Buffer for reading data into

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READ_MULTI'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------

iostatus    = -1.0
len_io      = isize
local_len   = 0
errorstatus = 0
expanded    = .FALSE.

! First thing to do is to read the field in to PE 0

CALL set_unit_bcast_flag(nft) ! Switch off IO broadcasting

IF (mype == 0) THEN

  ! Check for packed data and call the correct buffin version
  ! accordingly.
  IF (MOD((lookup(lbpack)),10) == PC_Cray32_Packing ) THEN
    ! Data is packed using CRAY 32 bit method

    ! For packed fields we need to read in 2*ISIZE 32 bit words
    ! using BUFFIN32
    ! DEPENDS ON: buffin32_f77
    CALL buffin32_f77(nft,buf,2*isize,len_io,iostatus)
    ! And then halve len_io to satisfy tests against ISIZE
    len_io = len_io/2
    ! We must check to see if it is a 32 bit field on disk, and if
    ! so, expand it before doing anything with it.
    IF (stash_record % grid_type == ppx_atm_compressed .OR.   &
        stash_record % grid_type == ppx_atm_tzonal     .OR.   &
        stash_record % grid_type == ppx_atm_uzonal) THEN
      ! DEPENDS ON: expand32b
      CALL expand32b( lookup(lblrec) , buf, fixhd12 )

      expanded = .TRUE.
    END IF
    ! What about WGDOS packing?
  ELSE IF (MOD((lookup(lbpack)),10) == PC_WGDOS_Packing) THEN
    ! Data is packed using WGDOS packing
    ! Read in as INTEGER, conversion to real done by wgdos_expand_field
    CALL buffin(nft,int_buf,isize,len_io,iostatus)

    CALL wgdos_expand_field( buf, expand_size, int_buf, isize, idum, dimx,    &
                             dimy, idum, rmdi, lookup(item_code), errorstatus)

    IF ( ErrorStatus /= 0 ) THEN
      CALL Ereport( RoutineName, ErrorStatus,  &
           'Problems reported by wgdos_expand_field' )
    END IF
    expanded = .TRUE.
  ELSE
    ! For non-packed data
    CALL buffin(nft,buf,isize,len_io,iostatus)
  END IF

  ! Lets do a polar row check/averaging (on p-grids) for read data only
  IF (lookup(lbhem) == 0 .AND. lookup(lbegin) == 1 .AND.  &
      Stash_Record % grid_type <= 3 .AND. expanded ) THEN
    CALL Rcf_Average_Polar( buf, glsizep(2), glsizep(1), global, averaged)

    !  Print a warning if we've done averaging
    IF (averaged .AND. PrintStatus >= PrStatus_Normal) THEN
      WRITE(umMessage,*) 'Field had had polar rows averaged in Rcf_Read_Multi'
      CALL umPrint(umMessage,src='rcf_read_multi_mod')
    END IF
  END IF

END IF  ! (mype == 0)

CALL clear_unit_bcast_flag(nft)! Restore broadcast behaviour

! Now we can distribute it out to the other processes
IF (MOD((lookup(lbpack)),10) == PC_Cray32_Packing .AND.                    &
    stash_record % grid_type /= ppx_atm_compressed .AND.   &
    stash_record % grid_type /= ppx_atm_tzonal     .AND.   &
    stash_record % grid_type /= ppx_atm_uzonal   ) THEN

  ! Distribute in 32-bit and expand
  fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

  global_x = glsize(1, fld_type)
  global_y = glsize(2, fld_type)
  local_x  = g_blsize(1, fld_type, mype)
  local_y  = g_blsize(2, fld_type, mype)

  local_len = local_x * local_y
  ! DEPENDS ON: scatter_field_mpl32
  CALL Scatter_Field_mpl32(d1,          buf,              &
                           local_x,     local_y,          &
                           global_x,    global_y,         &
                           fld_type,    halo_type_single, &
                           0,  gc_all_proc_group )

  ! DEPENDS ON: expand32b
  CALL expand32b( local_len, d1, fixhd12 )

ELSE
  ! Distribute the expanded data
  CALL Rcf_General_Scatter_Field( d1, buf, &
      local_len, lookup(lblrec), stash_record, 0 )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Read_Multi

END MODULE Rcf_Read_Multi_Mod
