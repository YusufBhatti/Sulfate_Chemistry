! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE pr_corner_mod

! Description:      Parameters used to define corners of
!                   diagnostic print domain
!
! Method:
!               The required parameters are initialised at start
!               and are used by subroutine  print_dom_corn
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

INTEGER  :: dom_w
INTEGER  :: dom_eo    ! no offset yet
INTEGER  :: dom_s
INTEGER  :: dom_no    ! no offset yet
INTEGER  :: ix
INTEGER  :: jy

END MODULE pr_corner_mod
