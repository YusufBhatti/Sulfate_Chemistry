! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Unifies two packed fragments of a data field into one
!
! output += input
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

MODULE unite_output_files_mod

IMPLICIT NONE

INTEGER, PARAMETER :: dataOffset=3

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UNITE_OUTPUT_FILES_MOD'

CONTAINS

SUBROUTINE unite_wgdos_files( &
    appendBuffer,  &
    targetBuffer,  &
    num_out,       &
    blockID)

USE fort2c_data_conv_interfaces, ONLY: ibm2ieee, ieee2ibm, integer_type
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod
IMPLICIT NONE

! Description:
!   Combines two blocks of WGDOS packed data into a unified packed
!   data field.

! Method:
!   If the row block is the first one (ID=0), the input buffer is
!   copied (with endianness corrections) into the output buffer.
!   If it is a later one, the input buffer is copied into the output
!   buffer from word 4 (32-bit) onwards. After this, the WGDOS header
!   is updated with the new sizes (see UMDP F3 for WGDOS header/data
!   structure)

REAL(KIND=real64), INTENT(IN)    :: appendBuffer(*)
!                         Buffer holding the current
!                         block of compressed data, which
!                         is to be added to the output  buffer

REAL(KIND=real64), INTENT(INOUT) :: targetBuffer(*)
!                         The output buffer that holds all
!                         the collected, compressed data

INTEGER, INTENT(IN)    :: blockID
!                         Sequence id in the set of WGDOS to accumulate
!                         only relevent if zero in which case the WGDOS
!                         global header is included

INTEGER, INTENT(OUT)   :: num_out
!                         Final Length of the target data buffer
!                         in 32 bit words

! Local variables:
INTEGER :: num_in
INTEGER :: iy_in
INTEGER :: iy_out
INTEGER :: ierr

INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UNITE_WGDOS_FILES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find the length of the additional entry
ierr=ibm2ieee(integer_type,1,appendBuffer(1),0,num_in,1,BIT_SIZE(num_in),32)

! If this is the global_gather_pe, we need the
! control information from the global header
IF (blockID == 0) THEN
  ! The first block has no offset in either the target or append
  ! file. The copy is therefore a straight copy.

  ! Round up to whole underlying 64-bit word boundary, in case there is a
  ! trailing blank upper 32-bit block copied due to little-endianism
  DO i = 1, (num_in+1)/2
    targetBuffer(i) = appendBuffer(i)
  END DO

  num_out = num_in

ELSE
  ! New data to add starts in word 4 (32-bit) of appendBuffer, words 1-3
  ! are the global header. None 64bit aligned offsets may result in
  ! 32 bit word swapped copies.

  ! Find the initial total length in the target buffer.
  ierr=ibm2ieee(integer_type,1,targetBuffer(1),0,num_out,1,BIT_SIZE(num_out),32)

  ! DEPENDS ON: copy_buffer_32
  CALL copy_buffer_32( &
      appendBuffer,    &
      targetBuffer,    &
      num_in,          &
      num_out,         &
      dataOffset)

  num_out = num_out + (num_in - dataOffset)

  ! Update the number of rows in the global header of the target file.
  ierr=ibm2ieee(integer_type,1,targetBuffer(2),16,iy_out,1,BIT_SIZE(iy_out),16)
  ierr=ibm2ieee(integer_type,1,appendBuffer(2),16,iy_in,1,BIT_SIZE(iy_in),16)
  iy_out=iy_out+iy_in

  ierr=ieee2ibm(integer_type,1,targetBuffer(2),16,iy_out,1,BIT_SIZE(iy_out),16)

END IF

! Update the number of 32-bit words in the global header of the target
! file.
ierr=ieee2ibm(integer_type,1,targetBuffer(1),0,num_out,1,BIT_SIZE(num_out),32)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE unite_wgdos_files

END MODULE unite_output_files_mod
