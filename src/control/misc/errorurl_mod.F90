! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Description
!   This module contains the URL to the Known UM Failure Points
!   documentation. This documentation will provide guidnace on
!   how to deal with the given error.

MODULE errorurl_mod

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Set the URL to Known Error documentation. This could be either a
! webpage, twiki page or UM Trac or UMDP.

! These lines need to be joined by newline from umPrintMgr or ereport_mod.
! This can only be done in the routine where ereport is being called as 
! newline is not set until runtime, not at compile time.

CHARACTER(LEN=errormessagelength) :: errorurl1 = &
  'See the following URL for more information:'
CHARACTER(LEN=errormessagelength) :: errorurl2 = &
  'https://code.metoffice.gov.uk/trac/um/wiki/KnownUMFailurePoints'

END MODULE errorurl_mod

