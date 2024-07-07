! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine FIND_MAX_FIELD_SIZE ---------------------------------
!
!    Purpose:  Reads in and searches LOOKUP header for maximum field
!              size
!
!
!    Logical component number: E5
!
!    External Documentation: None
!
!
!    Arguments:--------------------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Small execs

SUBROUTINE find_max_field_size(                                   &
  nftin,len1_lookup,len2_lookup,fixhd,max_field_size              &
  )

USE io
USE io_configuration_mod, ONLY: io_field_padding
USE ereport_mod, ONLY: ereport

USE lookup_addresses
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER ::                                                        &
 nftin                                                            &
                !IN Unit number of file
,len1_lookup                                                      &
                !IN 1st dim of LOOKUP array
,len2_lookup                                                      &
                !IN 2nd dim of LOOKUP array
,fixhd(*)                                                         &
                !IN Fixed length header
,max_field_size !OUT Maximum size of field held on file

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------

INTEGER ::                                                        &
 lookup(len1_lookup,len2_lookup) ! Lookup header

! -------------------------------------------------------------
! Local variables:---------------------------------------------

INTEGER ::                                                        &
 len_io                                                           &
           ! No of words transferred by BUFFIN
,k                                                                &
           ! Loop index

,icode     !Return code from setpos

INTEGER :: err_code
REAL ::                                                           &
 a         ! BUFFIN error code

! -------------------------------------------------------------

!  Internal structure: none

! Move to start of Look Up Table
CALL setpos(nftin,fixhd(150)-1,icode)

! Read in fields from LOOKUP table
CALL buffin(nftin,lookup(:,:),fixhd(151)*fixhd(152),len_io,a)

! Check for I/O errors
IF (a /= -1.0 .OR. len_io /= fixhd(151)*fixhd(152)) THEN
  ! DEPENDS ON: ioerror
  CALL ioerror('buffer in of lookup table',a,len_io,              &
               fixhd(151)*fixhd(152))

  err_code = 1000
  CALL ereport('FIND_MAX_FIELD_SIZE', err_code,                   &
   'buffer in of lookup table wrong size')

END IF

! Find maximum field size
max_field_size=0
DO k=1,len2_lookup
  max_field_size=MAX(max_field_size,lookup(15,k))
  max_field_size=MAX(max_field_size,lookup(30,k))
END DO

RETURN
END SUBROUTINE find_max_field_size
