! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE INIT_FLH --------------------------------------
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    System component: R30
!
!    Purpose:
!             Initialises the fixed length header to IMDI except
!             for all array dimensions which are set to 1.
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!  ------------------------------------------------------------
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O
SUBROUTINE init_flh (fixhd,len_fixhd)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE

INTEGER ::                                                        &
 len_fixhd                                                        &
                  ! IN    Length of fixed length header
,fixhd(len_fixhd) ! INOUT Fixed length header

! Local arrays:------------------------------------------------
! None
! -------------------------------------------------------------
! Local variables:---------------------------------------------
INTEGER :: j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_FLH'
! -------------------------------------------------------------

! 1.0 Initialise to IMDI
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO j = 1,len_fixhd
  fixhd(j) = imdi
END DO

! 2.0 Set all array dimensions to 1
fixhd(101) = 1     !  Integer Constants
fixhd(106) = 1     !  Real Constants
fixhd(111) = 1     !  1st dim - Level dependent constants
fixhd(112) = 1     !  2nd dim - Level dependent constants
fixhd(116) = 1     !  1st dim - Row dependent constants
fixhd(117) = 1     !  2nd dim - Row dependent constants
fixhd(121) = 1     !  1st dim - Column dependent constants
fixhd(122) = 1     !  2nd dim - Column dependent constants
fixhd(126) = 1     !  1st dim - Field of constants
fixhd(127) = 1     !  2nd dim - Field of constants
fixhd(131) = 1     !  Extra constants
fixhd(136) = 1     !  Temp History file
fixhd(141) = 1     !  Compressed field Index 1
fixhd(143) = 1     !  Compressed field Index 2
fixhd(145) = 1     !  Compressed field Index 3
fixhd(151) = 1     !  1st dim - Lookup Table
fixhd(152) = 1     !  2nd dim - Lookup Table
fixhd(161) = 1     !  Data

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_flh
