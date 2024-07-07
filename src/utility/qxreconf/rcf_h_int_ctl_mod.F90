! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Controls horizontal interpolation

MODULE Rcf_H_Int_Ctl_Mod

!  Subroutine rcf_h_int_ctl - controls horizontal interpolation.
!
! Description:
!   Wrapper routine choosing correct horizontal interpolation method
!
! Method:
!   Chooses method based on h_int_method variable.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_H_INT_CTL_MOD'

CONTAINS

SUBROUTINE Rcf_H_INT_CTL(len_field_out, row_length_in, row_length_out, &
                        rows_in, rows_out,global,aw_index_targ_lhs,    &
                        aw_index_targ_top, bl_index_b_l,bl_index_b_r,  &
                        bl_index_t_l,bl_index_t_r,                     &
                        aw_colat_t,aw_long_l,data_in,weight_t_r,       &
                        weight_b_r,weight_t_l,weight_b_l, data_out)

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_method,                  &
    bilinear,                      &
    area_weighted,                 &
    nearest_neighbour

USE Rcf_H_Int_Nearest_Mod, ONLY: &
    Rcf_H_Int_Nearest

USE Ereport_mod, ONLY: &
    Ereport

USE f_shum_horizontal_field_interp_mod, ONLY :                                 &
                                     f_shum_horizontal_field_bi_lin_interp_calc

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER ::   len_field_out    !No of points on target grid
INTEGER ::   rows_in          !No of rows on source grid
INTEGER ::   rows_out         !No of rows on target grid
INTEGER ::   row_length_in    !No of pts per row on source grid
INTEGER ::   row_length_out   !No of pts per row on target grid
LOGICAL ::   global           !True if global area required

!   Array  arguments with intent(in):
INTEGER, POINTER ::      aw_index_targ_lhs(:)
                              !Index of source box overlapping
                              !lhs of target grid-box
INTEGER, POINTER ::      aw_index_targ_top(:)
                              !Index of source box overlapping
                              !top of target grid-box
INTEGER, POINTER ::      bl_index_b_l(:)
                              !Gather index for bottom l.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      bl_index_b_r(:)
                              !Gather index for bottom r.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      bl_index_t_l(:)
                              !Gather index for top l.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      bl_index_t_r(:)
                              !Gather index for top r.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
REAL, POINTER ::         aw_colat_t(:)
                              !Colatitude of top of target grd-box
                              ! (in units of DELTA_LAT_SRCE)
REAL, POINTER ::         aw_long_l(:)
                              !Left longitude of target grid-box
                              ! (in units of DELTA_LONG_SRCE)
REAL          ::         data_in(:)
                              !Data before interpolation
REAL, POINTER ::         weight_t_r(:) !\ Weights used in
REAL, POINTER ::         weight_b_r(:) ! \bilinear horizontal
REAL, POINTER ::         weight_t_l(:) ! /interpolation
REAL, POINTER ::         weight_b_l(:) !/ 1=P-pts; 2=U-pts
                                       !  3=V-pts; 4=zonal
                                       !             means
!   Array  arguments with intent(out):
REAL ::      data_out(*)
                              !Data after interpolation
INTEGER                      :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_H_INT_CTL'
CHARACTER (LEN=errormessagelength)      :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!- End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE ( h_int_method )
CASE ( bilinear )

  ! Bi-linear interpolation requested

  CALL f_shum_horizontal_field_bi_lin_interp_calc                            &
               (rows_in,row_length_in,len_field_out,                         &
                bl_index_b_l,bl_index_b_r,bl_index_t_l,bl_index_t_r,data_in, &
                weight_b_l,weight_b_r, weight_t_l,weight_t_r,                &
                data_out)

CASE ( nearest_neighbour )

  ! Nearest neighbour
  CALL Rcf_H_Int_Nearest(rows_in,row_length_in,len_field_out,          &
                         bl_index_b_l,bl_index_b_r,data_in,            &
                         weight_b_l,weight_b_r, weight_t_l,weight_t_r, &
                         data_out)

CASE ( area_weighted )

  !  Area weighted interpolation

  ! DEPENDS ON: h_int_aw
  CALL h_int_aw(rows_in,rows_out,                    &
                row_length_in,row_length_out,global, &
                aw_index_targ_lhs,aw_index_targ_top, &
                aw_colat_t,aw_long_l,data_in,data_out)




CASE DEFAULT
  Cmessage = 'Unsupported interpolation method'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_H_Int_Ctl



END MODULE Rcf_H_Int_Ctl_Mod
