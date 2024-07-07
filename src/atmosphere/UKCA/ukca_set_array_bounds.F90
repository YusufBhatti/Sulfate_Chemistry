MODULE ukca_set_array_bounds_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SET_ARRAY_BOUNDS_MOD'

CONTAINS

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To set array dimensions using halo type.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
SUBROUTINE ukca_set_array_bounds(n,i1,i2,j1,j2)
USE ukca_d1_defs
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
IMPLICIT NONE

INTEGER, INTENT(IN)     :: n  ! id of array
INTEGER,  INTENT(OUT)   :: i1,i2,j1,j2 ! array bounds

INTEGER :: errcode
INTEGER :: halox      ! halo EW
INTEGER :: haloy      ! halo EW

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_ARRAY_BOUNDS'



!     set haloes

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

halox=halosize(1,ukcaD1codes(n)%halo_type)
haloy=halosize(2,ukcaD1codes(n)%halo_type)

i1 = 1-halox
i2 = ukcaD1codes(n)%len_dim1+halox
j1 = 1-haloy
j2 = ukcaD1codes(n)%len_dim2+haloy
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_set_array_bounds
END MODULE ukca_set_array_bounds_mod

