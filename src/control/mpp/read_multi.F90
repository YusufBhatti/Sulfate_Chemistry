! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM interface to BUFFIN - primarily for small execs

MODULE read_multi_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: read_multi

CHARACTER(LEN=*), PARAMETER :: ModuleName='READ_MULTI_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE read_multi(nft,d1,isize,unpack_size,len_io,local_len,  &
                      lookup,fixhd12,                             &
                      grid_code, halo_code, object_type,          &
                      proc_no_code, length, north_code,&
                      south_code, east_code, west_code,           &
                      gridpoint_code, imodl, section, item,       &
                      n_levels,                                   &
                      icode,cmessage,expand)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE UM_ParVars
USE UM_ParCore,  ONLY: nproc
USE Decomp_DB
USE lookup_addresses
USE cppxref_mod, ONLY:                                               &
    ppx_atm_tall, ppx_atm_cuall, ppx_atm_cvall,                      &
    ppx_atm_ozone,ppx_atm_tzonal
USE d1_array_mod, ONLY: prognostic
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: rmdi
USE errormessagelength_mod, ONLY: errormessagelength
USE submodel_mod, ONLY: atmos_im
USE mpp_conf_mod, ONLY: scatter_sync
USE read_serial_mod, ONLY: read_serial
USE Packing_Codes_Mod, ONLY: PC_Cray32_Packing

IMPLICIT NONE


! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  Unified Model. It is used where each process must read in a
!  local section of a global field.

! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields on P/U/V points are distributed in 32-bit (if they
!  were compressed) and then expanded. Other fields are
!  expanded before distribution

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) ::  nft         !  FORTRAN unit number
INTEGER, INTENT(IN) ::  isize       !  no. of words to be read in
                                    !  (global field)
INTEGER, INTENT(IN) ::  unpack_size !  no. of words after any unpacking (global)
INTEGER, INTENT(OUT) :: len_io      !  no. of words read in (global field)
INTEGER, INTENT(OUT) :: local_len   !  no. of words in local field

INTEGER, INTENT(IN) :: expand      ! expansion for small exec

INTEGER, INTENT(IN) :: fixhd12     ! 12th element of fixed length header
INTEGER, INTENT(IN) :: n_levels    ! Number of levels
INTEGER, INTENT(OUT) :: icode      ! Return code

INTEGER, INTENT(IN) :: grid_code      ! grid code
INTEGER, INTENT(IN) :: halo_code      ! halo code
INTEGER, INTENT(IN) :: object_type    ! field type
INTEGER, INTENT(IN) :: proc_no_code   ! Processing Code
INTEGER, INTENT(IN) :: length         ! Record length of field
INTEGER, INTENT(IN) :: north_code     ! northern row of field
INTEGER, INTENT(IN) :: south_code     ! southern row of field
INTEGER, INTENT(IN) :: east_code      ! eastern row of field
INTEGER, INTENT(IN) :: west_code      ! western row of field
INTEGER, INTENT(IN) :: gridpoint_code ! gridpoint info of field
INTEGER, INTENT(IN) :: imodl          ! modle code of field
INTEGER, INTENT(IN) :: section        ! section code of field
INTEGER, INTENT(IN) :: item           ! item code of field

INTEGER, INTENT(INOUT) :: lookup(64)               !   LOOKUP header from dump

CHARACTER(LEN=errormessagelength) ::  cmessage       !  Error message
REAL, INTENT(OUT) ::  d1(:)          !  Array to read data in to


! Local variables
INTEGER :: info          ! GCOM return code
INTEGER :: i             ! loop counter
INTEGER :: map(n_levels) ! Which processor holds which level

INTEGER :: orig_decomp  ! decomposition on entry
INTEGER :: new_decomp   ! decomposition to change to

INTEGER :: grid_type    ! grid type

INTEGER :: local_expand ! do we expand on this PE? 1 = yes, 0 = no

REAL    :: buf_icode    ! return code from BUFFIN

REAL    :: buf(unpack_size)       ! buffer for reading the field into
REAL    :: local_buf(unpack_size) ! and distributed - this is too big

LOGICAL :: scatter32              ! will scatter the field in 32-bit

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_MULTI'

INTEGER :: get_fld_type  ! function
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
len_io=isize
local_len=0

! If it's a compressed REAL field, expand it out if requested
! (except for fields we want to distribute in 32-bit)

scatter32 = .FALSE.

IF ( (grid_code == ppx_atm_tall .OR. grid_code == ppx_atm_cuall .OR.           &
      grid_code == ppx_atm_cvall)                                              &
    .AND. MOD((lookup(lbpack)),10) == PC_Cray32_Packing ) THEN

  scatter32 = .TRUE.

END IF

IF (scatter32) THEN

  local_expand=0

ELSE

  local_expand=expand

END IF

! First thing to do is to read the field in to PE 0

IF (mype  ==  0) THEN

  CALL read_serial(nft, buf, isize, unpack_size, len_io, lookup, fixhd12,      &
                   imodl, section, item, icode, cmessage, local_expand)

END IF ! IF (mype  ==  0)

! For certain 32-bit packed fields we can distribute in 32-bit and
! thus unpack in parallel
IF (scatter32) THEN

  ! Ensure local buffer is initialised so that expand/pack operations
  ! won't hit uninitialised memory.
  local_buf(:) = rmdi

  ! DEPENDS ON: get_fld_type
  grid_type = get_fld_type(grid_code)

  ! Scatter followed by unpack for this case
  ! DEPENDS ON: scatter_field_mpl32
  CALL scatter_field_mpl32( local_buf, buf,                           &
                            lasize(1,grid_type,halo_code),            &
                            lasize(2,grid_type,halo_code),            &
                            glsize(1,grid_type), glsize(2,grid_type), &
                            grid_type, halo_code, 0,                  &
                            gc_all_proc_group)

  local_len=length

  IF (expand /= 0) THEN
    ! DEPENDS ON: expand21
    CALL expand21(local_len, local_buf, d1)
  ELSE
    local_len=local_len/2
    d1(1:local_len)=local_buf(1:local_len)
  END IF

ELSE

  ! Now decompose the field in buf to the local D1 arrays

  ! Create map showing which processor contains which level.
  ! At present, all levels are stored on PE 0
  DO i=1,n_levels
    map(i)=0
  END DO

  ! DEPENDS ON: general_scatter_field
  CALL general_scatter_field( d1, buf, local_len, lookup(lblrec), n_levels,    &
                              grid_code, halo_code, object_type, proc_no_code, &
                              length, 1, north_code, south_code,       &
                              east_code, west_code, gridpoint_code, map,       &
                              icode, cmessage)

  IF (icode  ==  1) THEN
    WRITE(umMessage,'(A)') 'READ_MULTI: Call to GENERAL_SCATTER_FIELD failed'
    CALL umPrint(umMessage,src='read_multi')
    WRITE(umMessage,'(A,I0)') 'Return code was ',icode
    CALL umPrint(umMessage,src='read_multi')
    WRITE(umMessage,'(A,A)') 'Error message was ',cmessage
    CALL umPrint(umMessage,src='read_multi')
    WRITE(umMessage,'(A,I0)') 'Field number ',lookup(item_code)
    CALL umPrint(umMessage,src='read_multi')
    WRITE(umMessage,'(A,I0)') 'Grid type ', grid_code
    CALL umPrint(umMessage,src='read_multi')
    WRITE(umMessage,'(A,I0,A,I0)') 'Dimensioned : ',                           &
                                   lookup(lbnpt),' x ',lookup(lbrow)
    CALL umPrint(umMessage,src='read_multi')
    icode=300
    cmessage='Failure decomposing field'
    GO TO 9999
  END IF

END IF  ! default expand method

! A sync can help avoid MPI congestion here.
IF (scatter_sync) THEN
  CALL gc_gsync(nproc, icode)
END IF

9999  CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_multi

END MODULE read_multi_mod
