! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!             ENDGame semi-Lagrangian parameters constant for run          
!          
!  Code Owner: Please refer to the UM file CodeOwners.txt
!                 This file belongs in section: Top Level
!
! Method:
!   Parameters set from SETCONA/atm_step_consts so only need to be set  
!    at start of run. 
!   ATM_STEP_4A does not need to see or use these parameters.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
MODULE sl_param_mod

IMPLICIT NONE

INTEGER  ::  flat_level

INTEGER, ALLOCATABLE ::  g_row_length(:)
                         ! Table of number of points on a row
INTEGER, ALLOCATABLE ::  g_rows(:)
                         ! Table number of rows in theta field
INTEGER, ALLOCATABLE ::  g_i_pe(:)
              ! proc on my proc-row holding a given value in i direction
INTEGER, ALLOCATABLE ::  g_j_pe(:)
              ! proc on my proc-col holding a given value in j direction

END MODULE sl_param_mod

