! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Information for defining precisions for variables in terms of words

MODULE word_sizes_mod

IMPLICIT NONE

! Description:
!   Allow different precisions to be used with the aim of reducing memory
!   ( From Ben Shipway)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

INTEGER, PARAMETER  ::                  &
  OneByteInt   = SELECTED_INT_KIND(2)   &
 ,TwoByteInt   = SELECTED_INT_KIND(4)   &
 ,FourByteInt  = SELECTED_INT_KIND(9)   &
 ,EightByteInt = SELECTED_INT_KIND(18)

INTEGER, PARAMETER  ::                            &
  FourByteReal  = SELECTED_REAL_KIND(p=6, r=37)   &
 ,EightByteReal = SELECTED_REAL_KIND(p=13, r=307)

INTEGER, PARAMETER  ::    &
  wp = FourByteReal       & ! Real working precision for arrays
 ,iwp = TwoByteInt        & ! Integer working precision (weights values)
 ,crm_wp = FourByteReal     ! Real working precision for CRM sampling arrays

INTEGER, PARAMETER  ::    &
  sp=FourByteReal         & ! single precision
 ,dp=EightByteReal          ! double preciison

!------------------------------------------------------------------------------
END MODULE word_sizes_mod
