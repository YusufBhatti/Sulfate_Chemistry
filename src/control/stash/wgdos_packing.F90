! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! PURPOSE:   TO INTERFACE TO SHUMLIB LIBRARY ALLOWING FIELDS
!            TO BE PACKED AND UNPACKED FROM WGDOS FORMAT
!
! -------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

MODULE wgdos_packing_mod

USE f_shum_wgdos_packing_mod,   ONLY: f_shum_read_wgdos_header,   &
                                      f_shum_wgdos_unpack,        &
                                      f_shum_wgdos_pack
USE ereport_mod,                ONLY: ereport
USE umPrintMgr,                 ONLY: newline
USE errormessagelength_mod,     ONLY: errormessagelength
USE UM_types,                   ONLY: integer64, real64, integer32
USE yomhook,                    ONLY: lhook, dr_hook
USE parkind1,                   ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'WGDOS_PACKING_MOD'


CONTAINS

SUBROUTINE wgdos_expand_field(uncomp_r_field,     field_size,       &
                              icomp,              comp_data_size,   &
                              no_comp_vals,                         &
                              row_len,            num_rows,         &
                              accuracy,           rmdi,             &
                              istash,             icode)


IMPLICIT NONE

INTEGER :: field_size        ! Size of an uncompressed field
INTEGER :: comp_data_size    ! Size of compressed data block
INTEGER :: no_comp_vals      ! Total number of compressed values returned.
INTEGER :: row_len, num_rows ! Dimensions of field
INTEGER :: accuracy          ! Accuracy in power of 2
INTEGER :: istash            ! STASH section/item of field
INTEGER :: icode             ! error code returned from external for checking

REAL(KIND=real64) :: uncomp_r_field(field_size) ! uncompressed field
REAL              :: rmdi

INTEGER(KIND=integer32) :: icomp_32(2*comp_data_size) ! 32-bit working array
INTEGER(KIND=integer64) :: icomp(comp_data_size)      ! 64-bit array where each
                                                      ! element holds 2 32-bit
                                                      ! numbers
INTEGER(KIND=integer64), PARAMETER :: mask32 = INT(z'FFFFFFFF', KIND=integer64)

INTEGER :: i,idx  ! loop counter and index value in array
CHARACTER(LEN=errormessagelength)  :: cmessage_temp
CHARACTER(LEN=errormessagelength)  :: cmessage

! Dr. Hook
INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle
CHARACTER(LEN=*), PARAMETER        :: RoutineName='WGDOS_EXPAND_FIELD'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Currently (for historical reasons) the input packed array (icomp) is a
! 64-bit array, so here extract the two 32-bit values that make up each of
! those 64-bit words (since the Shumlib unpacking expects a 32-bit array)
DO i = 1,comp_data_size
  idx = (2*i)-1
  icomp_32(idx)   = INT(ISHFT(icomp(i), -32),   KIND=integer32)
  icomp_32(idx+1) = INT(IAND(icomp(i), mask32), KIND=integer32)
END DO

! Extract scale factor, cols and rows from field header
icode = f_shum_read_wgdos_header(icomp_32, no_comp_vals, accuracy, row_len,  &
                                 num_rows, cmessage_temp)

! Check for errors with reading the header and report
IF (icode /= 0) THEN
  IF (istash /= 0) THEN
    WRITE(cmessage, "(A,I6,A)") "Problem reading packed field header..."     &
     //newline//"STASH:", istash,                                            &
      newline//"Message: "//TRIM(cmessage_temp)
  ELSE
    cmessage = TRIM(cmessage_temp)
  END IF
  CALL ereport("WGDOS Un-Packing (f_shum_read_wgdos_header)", icode, cmessage)
END IF

! Use SHUMLib packing routine to unpack the field (note that this accepts
! icomp_32 - a 32-bit array of integers)
icode = f_shum_wgdos_unpack(icomp_32, rmdi, uncomp_r_field, row_len,         &
                            cmessage_temp)

! Check for errors with the unpacking and report
IF (icode /= 0) THEN
  IF (istash /= 0) THEN
    WRITE(cmessage, "(A,I6,A,I0,A)") "Problem expanding packed field..."     &
      //newline//"STASH:", istash,                                           &
        newline//"Accuracy: ", accuracy,                                     &
        newline//"Message: "//TRIM(cmessage_temp)
  ELSE
    cmessage = TRIM(cmessage_temp)
  END IF
  CALL ereport("WGDOS Un-Packing (f_shum_wgdos_unpack)", icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE wgdos_expand_field

SUBROUTINE wgdos_compress_field(uncomp_r_field,     field_size,       &
                                icomp,              comp_data_size,   &
                                stride,             no_comp_vals,     &
                                accuracy,           rmdi,             &
                                istash,             icode)


IMPLICIT NONE

INTEGER :: field_size        ! Size of an uncompressed field
INTEGER :: stride            ! Stride, or row length, for 2D field held as 1D
                             ! array.
INTEGER :: comp_data_size    ! Size of compressed data block
INTEGER :: no_comp_vals      ! Total number of compressed values returned.
INTEGER :: accuracy          ! Accuracy in power of 2
INTEGER :: istash            ! STASH section/item of field
INTEGER :: icode             ! error code returned from external for checking

REAL(KIND=real64) :: uncomp_r_field(field_size) ! uncompressed field
REAL              :: rmdi

INTEGER(KIND=integer32) :: icomp_32(comp_data_size) ! 32-bit working array
INTEGER(KIND=integer64) :: icomp(comp_data_size)    ! 64-bit array where each
                                                    ! element holds 2 32-bit
                                                    ! numbers

! Local variables
CHARACTER(LEN=errormessagelength)  :: cmessage_temp
CHARACTER(LEN=errormessagelength)  :: cmessage
INTEGER(KIND=integer64), PARAMETER :: mask32 = INT(z'FFFFFFFF', KIND=integer64)

INTEGER :: i,idx  ! loop counter and index value in array

! Dr. Hook
INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle
CHARACTER(LEN=*), PARAMETER        :: RoutineName='WGDOS_COMPRESS_FIELD'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use SHUMLib packing routine to pack the field (note that this returns
! icomp_32 - a 32-bit array of integers)
icode = f_shum_wgdos_pack(uncomp_r_field, stride, accuracy, rmdi, icomp_32, &
                          no_comp_vals, cmessage_temp)

! Check for errors with the packing and report
IF (icode /= 0) THEN
  IF (istash /= 0) THEN
    WRITE(cmessage, "(A,I6,A,I0,A,E18.10,A,E18.10,A)")                 &
     "Problem packing field..."//newline//"STASH:", istash,            &
     newline//"Accuracy: ", accuracy,                                  &
     newline//"Minimum:  ", MINVAL(uncomp_r_field),                    &
     newline//"Maximum:  ", MAXVAL(uncomp_r_field),                    &
     newline//"Message:  " //TRIM(cmessage_temp)
  ELSE
    cmessage = TRIM(cmessage_temp)
  END IF
  CALL ereport("WGDOS Packing (f_shum_wgdos_pack)", icode, cmessage)
END IF

! Currently (for historical reasons) this routine still returns a 64-bit
! array (icomp) so here pack the 32-bit array from above into it, so that
! each 64-bit word contains 2 of the 32-bit words
DO i = 1,(no_comp_vals+1)/2
  idx = (2*i)-1
  icomp(i) = IOR(IAND(INT(icomp_32(idx+1), KIND=integer64), mask32),   &
                 ISHFT(INT(icomp_32(idx),  KIND=integer64), 32))
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE wgdos_compress_field

END MODULE wgdos_packing_mod
