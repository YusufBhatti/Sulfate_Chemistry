! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE cray32_packing_mod

USE io,                     ONLY: setpos, buffout, buffin
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE um_types,               ONLY: real32, real64
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


! Description:
!   Routines used for CRAY32 packing.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CRAY32_PACKING_MOD'

CONTAINS

SUBROUTINE lbc_write_packed_field (                                         &
           nftout,lbc_data64,len_lbc_data,disk_address, len_data_disk)


IMPLICIT NONE

INTEGER :: nftout       ! Unit no for boundary file
INTEGER :: len_lbc_data ! Length of data buffer
INTEGER :: disk_address  ! disk address
INTEGER :: len_data_disk ! Number of 64-bit words required to contain the
                         ! 32-bit data

REAL(KIND=real64) :: lbc_data64(len_lbc_data) !  LBC data
REAL(KIND=real32) :: lbc_data32(len_lbc_data) !  LBC data

! Error reporting variables
CHARACTER (LEN=*), PARAMETER :: routinename= 'LBC_WRITE_PACKED_FIELD'
INTEGER   :: len_io        ! length of data actually written
REAL      :: a_io          ! return code from buffout
INTEGER   :: ErrorStatus   ! error flag
CHARACTER (LEN=errormessagelength)  :: CMessage   !  error message
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! DEPENDS ON: pack21
CALL pack21 (len_lbc_data,lbc_data64,lbc_data32)


CALL setpos  (nftout, disk_address, ErrorStatus)

! We double len_data_disk here because it's a number of 64-bit words, and
! buffout32_f77 expects the number of 32-bit words

! DEPENDS ON: buffout32_f77
CALL buffout32_f77 &
    (nftout, lbc_data32, 2*len_data_disk, len_io, a_io)

IF (a_io /= -1.0) THEN
  ErrorStatus = 10
  WRITE (cmessage,'(A)') ' Problem in BUFFOUT for LBCs'

  CALL ereport ( routinename, ErrorStatus, CMessage )

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lbc_write_packed_field




SUBROUTINE ff_write_packed_field(                                            &
           nftout, ff_data64, nx, ny, disk_address, len_data_disk)


IMPLICIT NONE

INTEGER :: nftout       ! Unit no for boundary file
INTEGER :: nx, ny, len_ff_data ! Length of data buffer
INTEGER :: disk_address  ! disk address
INTEGER :: len_data_disk ! Number of 64-bit words required to contain the
                         ! 32-bit data

REAL(KIND=real64) :: ff_data64(nx,ny) !  LBC data
REAL(KIND=real32) :: ff_data32(nx,ny) !  LBC data

! Error reporting variables
CHARACTER (LEN=*), PARAMETER :: routinename= 'FF_WRITE_PACKED_FIELD'
INTEGER   :: len_io        ! length of data actually written
REAL      :: a_io          ! return code from buffout
INTEGER   :: ErrorStatus   ! error flag
CHARACTER (LEN=errormessagelength)  :: CMessage   !  error message
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

len_ff_data = nx * ny

! DEPENDS ON: pack21
CALL pack21 (len_ff_data,ff_data64,ff_data32)

CALL setpos  (nftout, disk_address, ErrorStatus)

! We double len_data_disk here because it's a number of 64-bit words, and
! buffout32_f77 expects the number of 32-bit words

! DEPENDS ON: buffout32_f77
CALL buffout32_f77 &
    (nftout, ff_data32, 2*len_data_disk, len_io, a_io)

IF (a_io /= -1.0) THEN
  ErrorStatus = 10
  WRITE (cmessage,'(A)') ' Problem in BUFFOUT for LBCs'

  CALL ereport ( routinename, ErrorStatus, CMessage )

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ff_write_packed_field


END MODULE cray32_packing_mod

