! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE STASH_COMP_GRID----------------------
!
!    Compute grid descriptors
!    out_bdx,out_bdy,out_bzx and out_bzy from inputs
!    1)  samples
!    2)  grid type code
!    3)  ocean (REMOVED)
!    4)  real and integer headers
!    5)  w_mostcol and n_mostrow
!    6) processing code
!
!    Programming standard: U M DOC  Paper NO. 3
!
!    External documentation  C4
!
!

!
!    INTERFACE and ARGUMENTS:------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH
SUBROUTINE stash_comp_grid(                                       &
  out_bzx,out_bzy,out_bdx,out_bdy,                                &
  samples,st_grid,                                                &
  w_mostcol,n_mostrow,                                            &
  realhd,len_realhd,inthd,len_inthd,gr,                           &
  icode,cmessage)
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE lookup_addresses
USE stparam_mod, ONLY: st_riv_grid, st_uv_grid, st_cv_grid, &
                   st_cu_grid, block_size, zonal_mean_base, &
        field_mean_base, global_mean_base, merid_mean_base
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!
INTEGER :: samples                ! IN no of samples in period (times
!
INTEGER ::                                                        &
  icode                                                           &
                    !OUT   Return code from the routine
, len_realhd                                                      &
                    !IN    Length of the Real Constants
, len_inthd                                                       &
                    !IN    Length of integer constants
, inthd(len_inthd)  !IN    Integer constants
!
INTEGER ::                                                        &
  st_grid                                                         &
                    !IN    STASH horizontal grid type
, n_mostrow                                                       &
                    !IN    The most nrthrly row.
, w_mostcol                                                       &
                    !IN    The most westerly column
, gr                !IN    The type of processing done
!
REAL ::                                                           &
  realhd(len_realhd)  !IN  Real header
REAL ::                                                           &
  out_bzx,out_bdx,out_bzy,out_bdy  ! OUT grid descriptors
CHARACTER(LEN=errormessagelength) ::                              &
  cmessage           ! OUT error messages
!
!
!  Local Variables
!
INTEGER :: mean_code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STASH_COMP_GRID'
!
! ---------------------------------------------------------------------

!     Construct PP header
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
  out_bzx=0.0
  out_bdx=0.0
  out_bzy=0.0
  out_bdy=0.0
ELSE
    !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX for
  IF (st_grid == st_riv_grid) THEN
    out_bdy = 180.0/180
    out_bzy = realhd(3) - out_bdy*0.5
    out_bdx = 360.0/360
    out_bzx = realhd(4) - out_bdx*0.5
  ELSE
    IF (st_grid == st_uv_grid .OR. st_grid == st_cv_grid) THEN
      out_bzy=realhd(3)-realhd(2)/2.0 ! UV pts
    ELSE
      out_bzy=realhd(3)-realhd(2) ! Zeroth Lat OUT_BZY
    END IF
    !
    IF (st_grid == st_uv_grid .OR. st_grid == st_cu_grid) THEN
      out_bzx=realhd(4)-realhd(1)/2.0 !UV points
    ELSE
      out_bzx=realhd(4)-realhd(1) ! Zeroth Long OUT_BZX
    END IF
    out_bdx=realhd(1) ! Long intvl OUT_BDX
    out_bdy=realhd(2) ! Lat intvl OUT_BDY

  END IF
  !
  ! Add on offset for fields not starting from the origin
  !
  out_bzy=out_bzy                                                 &
      +(n_mostrow-1)*out_bdy
  out_bzx=out_bzx                                                 &
      +(w_mostcol-1)*out_bdx

  IF (out_bzx >= 360.0)                                           &
     out_bzx=out_bzx-360.0
  !
  ! If horizontal averaging has been applied to the output field,
  ! set OUT_BDX and/or OUT_BDY to the full domain extent
  !
  mean_code=(gr/block_size)*block_size
  IF (mean_code == zonal_mean_base .OR.                           &
      mean_code == field_mean_base .OR.                           &
      mean_code == global_mean_base) THEN
    out_bdx=REAL(inthd(6))*realhd(1)
  END IF
  IF (mean_code == merid_mean_base .OR.                           &
      mean_code == field_mean_base .OR.                           &
      mean_code == global_mean_base) THEN
    out_bdy=REAL(inthd(7))*realhd(2)
  END IF
END IF
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_comp_grid
