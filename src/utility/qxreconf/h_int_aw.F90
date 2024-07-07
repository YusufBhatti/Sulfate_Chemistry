! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Performs Area weighted horizitontal interpolation
!
! Subroutine Interface:
SUBROUTINE h_int_aw(rows_in,rows_out                              &
,                   row_length_in,row_length_out,global           &
,                   aw_index_targ_lhs,aw_index_targ_top           &
,                   aw_colat_t,aw_long_l,data_in,data_out)

!    System component: S121
!
!    System task: S1
!
!    Purpose:
!
!    Documentation:
!              The interpolation formulae are described in
!              unified model on-line documentation paper S1.
!
IMPLICIT NONE
!
! Description:
!   Carries out bi-linear horizontal interpolation using coefficients
!   and gather indices calculated in subroutine H_INT_CO
!
! Method:
!   See UMDP S1 for full desciption
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::   rows_in          !No of rows on source grid
INTEGER ::   rows_out         !No of rows on target grid
INTEGER ::   row_length_in    !No of pts per row on source grid
INTEGER ::   row_length_out   !No of pts per row on target grid
LOGICAL ::   global           !True if global area required

!   Array  arguments with intent(in):
INTEGER ::   aw_index_targ_lhs(row_length_out+1)
                              !Index of source box overlapping
                              !lhs of target grid-box
INTEGER ::   aw_index_targ_top(rows_out+1)
                              !Index of source box overlapping
                              !top of target grid-box
REAL ::      aw_colat_t(rows_out+1)
                              !Colatitude of top of target grd-box
                              ! (in units of DELTA_LAT_SRCE)
REAL ::      aw_long_l(row_length_out+1)
                              !Left longitude of target grid-box
                              ! (in units of DELTA_LONG_SRCE)
REAL ::      data_in(row_length_in*rows_in)
                              !Data before interpolation

!   Array  arguments with intent(out):
REAL ::      data_out(row_length_out*rows_out)
                              !Data after interpolation

! Local scalars:
INTEGER ::   i

! Local arrays:
REAL ::      boxsum(row_length_out,rows_out)
                              !Sum of data on target grid

!- End of header

!     1. Calculate sum of contribution from gridboxes

! DEPENDS ON: box_sum
CALL box_sum(row_length_in,rows_in,row_length_out,rows_out        &
,            aw_long_l,aw_colat_t                                 &
,            aw_index_targ_lhs,aw_index_targ_top                  &
,            global,data_out,data_in)


RETURN
END SUBROUTINE h_int_aw

