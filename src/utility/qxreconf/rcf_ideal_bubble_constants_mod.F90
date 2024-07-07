! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE rcf_ideal_bubble_constants_mod

IMPLICIT NONE
!
! Description: Parameters for bubble initialisation
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3. 

INTEGER, PARAMETER :: omit      = 0           ! No bubble
INTEGER, PARAMETER :: gaussian  = 1
INTEGER, PARAMETER :: block     = 2
INTEGER, PARAMETER :: sine      = 3
INTEGER, PARAMETER :: plume     = 4


INTEGER,PARAMETER  :: idl_max_num_bubbles = 3

END MODULE rcf_ideal_bubble_constants_mod
