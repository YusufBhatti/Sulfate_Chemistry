! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE lookup_table_mod
IMPLICIT NONE

!
! Description: Arrays used to locate departure points
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

INTEGER, POINTER :: i_lkup_u(:), i_lkup_p(:)
INTEGER, POINTER :: j_lkup_v(:), j_lkup_p(:)
REAL             :: delta_xi1_u, delta_xi1_p
REAL             :: delta_xi2_v, delta_xi2_p
REAL             :: base_xi1_u,  base_xi1_p
REAL             :: base_xi2_v,  base_xi2_p

CONTAINS

SUBROUTINE set_look_up_table()
USE horiz_grid_mod
!  USE lookup_table_mod
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows
USE UM_parvars, ONLY: halo_i, halo_j
IMPLICIT NONE
INTEGER  :: IS, ie, js, je

! is and ie are the start and end of the global p array in xi1

IS = 1-halo_i
ie = global_row_length+halo_i

! js and je are the start and end of the global p array in xi2

js = 1-halo_j
je = global_rows+halo_j

CALL calc_look_up(glob_xi1_p, IS,ie, halo_i, i_lkup_p,                     &
                  base_xi1_p, delta_xi1_p)
CALL calc_look_up(glob_xi2_p, js,je, halo_j, j_lkup_p,                     &
                  base_xi2_p, delta_xi2_p)

CALL calc_look_up(glob_xi1_u, IS-1,ie-1, halo_i, i_lkup_u,                 &
                  base_xi1_u, delta_xi1_u)
CALL calc_look_up(glob_xi2_v, js-1,je, halo_j, j_lkup_v,                   &
                  base_xi2_v, delta_xi2_v)

END SUBROUTINE set_look_up_table

SUBROUTINE calc_look_up(xi, istrt, iend, halo, i_look, base, delta)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: istrt, iend, halo
REAL,    INTENT(IN)  :: xi(istrt:iend)
INTEGER, POINTER     :: i_look(:)
REAL,    INTENT(OUT) :: base, delta
INTEGER              :: i, ii, Npnts
INTEGER              :: halo_pnts
REAL                 :: x

delta = 1.0e20
DO i = istrt, iend-1
  delta = MIN( delta, xi(i+1) - xi(i) )
END DO

base  = xi(istrt)
Npnts = ( xi(iend) - xi(istrt) )/delta + 1

halo_pnts = (  xi(istrt+halo) - xi(istrt) )/delta
IF (ASSOCIATED(i_look)) THEN
  DEALLOCATE(i_look)
END IF
ALLOCATE( i_look(0:Npnts) )

i_look = -999999
x = base
ii = istrt - 1
DO i = 0, Npnts !-halo_pnts
  IF ( x >= xi(ii+1) ) THEN
    ii = ii + 1
    IF ( ii == iend ) EXIT   ! Shouldn't need this!
  END IF
  i_look(i) = ii
  x  = x + delta
END DO

END SUBROUTINE calc_look_up
END MODULE lookup_table_mod
