MODULE OASIS_Grad_Calc_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

IMPLICIT NONE

CONTAINS

SUBROUTINE OASIS_Grad_Calc(field_in,grad_ew,grad_ns)
!
! Description:
! Calculate gradients for use in 2nd order conservative regridding
! with the OASIS-MCT coupler. This routine only works for fields
! on T grid points (i.e. non-vector feilds.)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!=======================================================================

USE Field_Types, ONLY: fld_type_p
USE um_types, ONLY: integer32, real64
USE oasis_atm_data_mod, ONLY: GradIndex, oasis_imt, oasis_jmt, fland_loc
USE mpp_conf_mod, ONLY: swap_field_is_scalar

IMPLICIT NONE

REAL (KIND=real64), INTENT(IN) :: field_in(oasis_imt,oasis_jmt) ! The field
                                  ! whose gradient we want to know
REAL (KIND=real64), INTENT(OUT) :: grad_ew(oasis_imt,oasis_jmt) ! Field
                                  ! gradient following lon (W-E) dirn.
REAL (KIND=real64), INTENT(OUT) :: grad_ns(oasis_imt,oasis_jmt) ! Field
                                  ! gradient following lat (S-N) dirn.

REAL (KIND=real64) :: field(0:oasis_imt+1,0:oasis_jmt+1)

REAL (KIND=real64) :: component_row
REAL (KIND=real64) :: component_col

INTEGER (KIND=integer32) :: i, j

! Actually, we  could avoid this copying by judicious addition of halos to 
! field_in (and all the outgoing UM fields!) but here we stick to using a 
! local copybecause there's also a payoff between this and extra memory use 
! for multiple fields (maybe thousands in 3D coupling) some, most or all of 
! which will not require 2nd order terms at all.
DO j = 1, oasis_jmt
  DO i = 1, oasis_imt
    field(i,j) = field_in(i,j)
  END DO
END DO
! Ensure our field has halos fully populated, luckily we only deal with
! fields on the p grid and our halos are fixed at 1 (see above) so we only
! need to know i and j at +/-1. 
CALL Swap_Bounds(field, oasis_imt, oasis_jmt, 1,           &
          1, 1, fld_type_p, swap_field_is_scalar)

! We don't initialise the gradient terms explicitly since all points 
! should acquire values naturally from the calculations below.
! Calculate the gradients along true lines of latitude and longitude because
! grid may be rotated!
DO j = 1, oasis_jmt
  DO i = 1, oasis_imt

    ! Take the difference of values along row (i direction) and divide by
    ! distance.  
    component_row = ((field(GradIndex(i,j)%ip1,j)      &
                     - field(GradIndex(i,j)%im1,j))  &
                     *GradIndex(i,j)%Rdi)

    ! Take the difference of values along col (j direction) and divide by
    ! distance.  
    component_col = ((field(i,GradIndex(i,j)%jp1)      &
                     - field(i,GradIndex(i,j)%jm1))  &
                     * GradIndex(i,j)%Rdj)

    ! Now calculate the W-E contributions from each potential component. 
    grad_ew(i,j) = (component_row * GradIndex(i,j)%dlon_row) + &
                   (component_col * GradIndex(i,j)%dlat_row)

    ! Divide by distance for the gradient
    grad_ew(i,j) = grad_ew(i,j) * GradIndex(i,j)%Rdi
 
    ! Now calculate the S-N contributions from each potential component. 
    ! along lines of longitude.
    grad_ns(i,j) = (component_row * GradIndex(i,j)%dlon_col) + &
                   (component_col * GradIndex(i,j)%dlat_col)

    ! Divide by distance for the gradient
    grad_ns(i,j) = grad_ns(i,j) * GradIndex(i,j)%Rdj

  END DO  ! over i
END DO ! over j

END SUBROUTINE OASIS_Grad_Calc

END MODULE OASIS_Grad_Calc_mod
