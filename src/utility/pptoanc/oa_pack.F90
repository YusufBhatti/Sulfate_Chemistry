! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routines: OA_PACK, OA_UNPACK and OA_LEV_CMP
!
!
!    Logical components covered:
!
!
!    Programming standard : FOAM Doc Paper 3/2/1
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Small execs
SUBROUTINE oa_pack(icode, cmessage, ll_ac_tim,                    &
 no_rows_m, no_cols_m, no_levs_m, no_seg, no_cmp_klev,            &
 indx_cmp, indx_exp, indx_to_rows, no_cmp, real_mdi,              &
 klev, i_fld_typ, ll_cyc_m, fld_exp, fld_cmp)
!
!    Purpose:  This subroutine packs one full field of model data
!              at a given level  into compressed form.
!
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
!    ARGUMENT LIST
!
INTEGER :: icode      ! return error code
CHARACTER(LEN=errormessagelength) :: cmessage      ! return error message
LOGICAL :: ll_ac_tim  ! T => time output on input and exit
! Dimensions
INTEGER :: no_rows_m  ! IN number of rows on model grid
INTEGER :: no_cols_m  ! IN number of columns on expanded model grid
INTEGER :: no_levs_m  ! IN number of levels on model grid
INTEGER :: no_seg     ! IN total number of sea segments
INTEGER :: no_cmp_klev! IN number of sea points on this level
! Compression and expansion indices
INTEGER :: indx_cmp(no_seg) ! IN contains position in compressed arra
!                                   of start of each sea segment
INTEGER :: indx_exp(no_seg) ! IN contains position in expanded array
!                                   of start of each sea segment
INTEGER :: indx_to_rows(no_rows_m*no_levs_m) ! IN contains number of
!                          first/next sea segment for each row and level
INTEGER :: no_cmp           ! IN total no of points in a 3D
!                                   compressed array
REAL :: real_mdi            ! IN missing data indicator
! Level and field type
INTEGER :: klev             ! IN model level
INTEGER :: i_fld_typ        ! IN =0 for tracers, =1 for currents
LOGICAL :: ll_cyc_m         ! T => FLD_EXP is cyclic in W-E
! Fields
REAL :: fld_exp(no_cols_m*no_rows_m) ! IN  data on expanded grid
REAL :: fld_cmp(no_cmp_klev)         ! OUT data on compressed grid
!
!    NO PARAMETERS
!
!    NO WORK ARRAYS
!
!    NO EXTERNAL SUBROUTINES CALLED
!
!    LIST OF OTHER VARIABLES
!
INTEGER :: icyc        ! number of columns of expanded field; sec. 1.
INTEGER :: icount      ! index for loop over points in segment
INTEGER :: inc_cyc     ! extra columns on cyclic grid; sec. 1
INTEGER :: ipt_cmp     ! pointer to location in compressed field
INTEGER :: ipt_exp     ! pointer to location in 3D expanded field
INTEGER :: ipt_exp_cyc ! pointer to location in 2D expanded cyclic fi
INTEGER :: ipt_seg     ! segment index number
INTEGER :: iseg        ! index for loop over segments in row
INTEGER :: iseg_st     ! first segment on this level
INTEGER :: ist_cmp_m1  ! index of 1st point in compressed field at th
!                           level, minus one
INTEGER :: ist_exp_m1  ! index of 1st point in expanded field at this
!                           level (for non-cyclic grid), minus one
INTEGER :: j           ! index for loop over rows on level
INTEGER :: jpt         ! point index number
INTEGER :: len_seg     ! number of grid points in current segment
INTEGER :: no_seg_row  ! number of segments in row
!
!-----------------------------------------------------------------------
!
!  1. Prelminaries
!
!  1.1 Set the number of columns of distinct data
IF (ll_cyc_m) THEN
  inc_cyc = 2
  icyc = no_cols_m - 2
ELSE
  inc_cyc = 0
  icyc = no_cols_m
END IF
!
!  1.2 Set offsets for compressed and expanded fields at this level
iseg_st = indx_to_rows(no_rows_m*(klev-1) + 1)
ist_cmp_m1 = indx_cmp(iseg_st) - 1
ist_exp_m1 = (klev-1)*no_rows_m*icyc
!
!  2. Loop over rows (index J)
!
!
!  2.2 Start loop over rows and define the pointer to the row
!
DO j = 1, no_rows_m
  jpt = (klev - 1)*no_rows_m + j
  !
  !  2.3 Calculate the number of sea segments in the row
  !
  IF (jpt  ==  no_levs_m*no_rows_m ) THEN
    no_seg_row = no_seg - indx_to_rows(jpt) + 1
  ELSE
    no_seg_row = indx_to_rows(jpt+1) - indx_to_rows(jpt)
  END IF
  !
  !  2.4 Start loop over sea segments and define pointer to segment
  DO iseg = 1, no_seg_row
    ipt_seg = indx_to_rows(jpt) + iseg - 1
    !
    !  2.5 Calculate the length of the present sea segment
    !
    IF (ipt_seg  <   no_seg) THEN
      len_seg = indx_cmp(ipt_seg+1) - indx_cmp(ipt_seg)
    ELSE
      len_seg = no_cmp - indx_cmp(ipt_seg) + 1
    END IF
    !
    !  2.6 Calculate FLD_CMP for all points in the segment
    !      (the last point may be overwritten in 2.7 if I_FLD_TYP = 1)
    !
    DO icount = 1, len_seg
      ipt_exp = indx_exp(ipt_seg) + icount - 1
      ipt_exp_cyc = ipt_exp - ist_exp_m1 + inc_cyc*(j-1)
      ipt_cmp = indx_cmp(ipt_seg) + icount - 1
      fld_cmp(ipt_cmp - ist_cmp_m1) = fld_exp(ipt_exp_cyc)
    END DO  ! index ICOUNT
    !
    !  2.7 Case of current field:
    !
    IF (i_fld_typ  ==  1) THEN
      !
      !   Last value in segment is only retained if grid is cyclic and
      !   the first point on the row is ICYC-1 points before it.
      IF (ll_cyc_m) THEN
        IF (ipt_exp-indx_exp(ipt_seg+1-iseg)  /=  icyc-1) THEN
          fld_cmp(ipt_cmp - ist_cmp_m1) = real_mdi
        END IF
      ELSE
        fld_cmp(ipt_cmp - ist_cmp_m1) = real_mdi
      END IF
    END IF
    !
  END DO  ! index ISEG
  !
END DO  ! index J
!
!   End loop over rows
!
RETURN
END SUBROUTINE oa_pack
!
