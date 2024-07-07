! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Controls horizontal interpolation

MODULE Rcf_H_Int_Nearest_Mod

!  Subroutine Rcf_H_Int_Nearest_Mod - controls horizontal interpolation.
!
! Description:
!   Add neareat method when h_int_method = 3. Changgui Wang (24/9/07)
!
! Method:
!   Chooses method based on h_int_method variable.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_H_INT_NEAREST_MOD'

CONTAINS

SUBROUTINE Rcf_H_Int_Nearest(rows_in,row_length_in,len_field_out,          &
                             bl_index_b_l,bl_index_b_r,data_in,            &
                             weight_b_l,weight_b_r, weight_t_l,weight_t_r, &
                             data_out)
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::   rows_in          !No of rows on source grid
INTEGER ::   row_length_in    !No of pts per row on source grid
INTEGER ::   len_field_out    !No of points on target grid

!   Array  arguments with intent(in):
INTEGER, POINTER ::      bl_index_b_l(:)
                              !Gather index for bottom l.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      bl_index_b_r(:)
                              !Gather index for bottom r.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
REAL          ::         data_in(:)
                              !Data before interpolation
REAL, POINTER ::         weight_b_l(:)  !\ Weights used in
REAL, POINTER ::         weight_b_r(:)  ! \bilinear horizontal
REAL, POINTER ::         weight_t_l(:)  ! /interpolation
REAL, POINTER ::         weight_t_r(:)  !/ 1=P-pts; 2=U-pts
                                        !  3=V-pts; 4=zonal
                                        !             means


!   Array  arguments with intent(out):
REAL ::      data_out(*)
                              !Data after interpolation
INTEGER                      :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_H_INT_NEAREST'

! Local
REAL        ::nest_mask(4)
INTEGER     ::i,NEAREST(1)

!- End of header

!      1. Carry out horizontal interpolation using nearest neighbour
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1,len_field_out
  nest_mask(1)          = weight_b_l(i)
  nest_mask(2)          = weight_b_r(i)
  nest_mask(3)          = weight_t_l(i)
  nest_mask(4)          = weight_t_r(i)
  NEAREST               = MAXLOC( nest_mask )
  nest_mask             = 0.0
  nest_mask(NEAREST(1)) = 1.0

  data_out(i) = nest_mask(1) * data_in(bl_index_b_l(i)) +  &
                nest_mask(2) * data_in(bl_index_b_r(i)) +  &
                nest_mask(3) * data_in(bl_index_b_l(i)  +  &
                row_length_in)                          +  &
                nest_mask(4) * data_in(bl_index_b_r(i)  +  &
                row_length_in)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_H_Int_Nearest
END MODULE Rcf_H_Int_Nearest_Mod
