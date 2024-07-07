! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE NUM_OBS2 -----------------------------------------------
!
!  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!
!  Project Task : P3
!
!  Purpose : Read in header section of observation file and get
!            number of observations and data values in this obs file.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE num_obs2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NUM_OBS2_MOD'

CONTAINS

SUBROUTINE num_obs2 (unit_no,p_rows,row_length, &
                     ak,bk,realhd1,realhd2,realhd3,realhd4,       &
                     realhd5,realhd6,                             &
                     len_data,nobtyp,tnobs,tndv,                  &
                     icode,cmessage)
! ----------------------------------------------------------------------
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE ac_dump_mod, ONLY: len_inthd, len_realhd, len1_levdepc,       &
                       len2_levdepc, len1_rowdepc, len2_rowdepc,  &
                       len1_coldepc, len2_coldepc, len1_flddepc,  &
                       len2_flddepc, len_extcnst, len_dumphist,   &
                       len_cfi1, len_cfi2, len_cfi3,              &
                       len1_lookup_obs, len2_lookup_obs, fixhd,   &
                       inthd, realhd, levdepc, rowdepc, coldepc,  &
                       flddepc, extcnst, dumphist, cfi1, cfi2,    &
                       cfi3, lookup
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod, ONLY: model_levels, len_fixhd

IMPLICIT NONE

!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
INTEGER :: unit_no       ! IN  : Unit no of Observation file
INTEGER :: p_rows        ! IN  : No of model rows
INTEGER :: row_length    ! IN  : No of points on row
REAL :: ak(model_levels) ! IN  : Vertical grid
REAL :: bk(model_levels)
REAL :: realhd1,realhd2 ! IN  : Horizontal grid
REAL :: realhd3,realhd4
REAL :: realhd5,realhd6
INTEGER :: len_data     ! IN  : Dimension of data section
INTEGER :: nobtyp       ! OUT : No of observation types
INTEGER :: tnobs        ! OUT : Total no of observations
INTEGER :: tndv         ! OUT : Total no of data values
INTEGER :: icode        ! OUT : Return code
CHARACTER(LEN=errormessagelength) :: cmessage!  OUT : Error message if ICODE > 0
!-----------------------------------------------------------------------
!     LEVEL/GRID VARIABLES
!-----------------------------------------------------------------------
INTEGER :: obs_row_length,obs_p_rows,obs_p_levels,obs_q_levels
REAL :: obs_ak(model_levels),obs_bk(model_levels)
REAL :: obs_long_res,obs_lat_res,obs_start_lat
REAL :: obs_start_long,obs_lat_pseudo_pole,obs_long_pseudo_pole
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
INTEGER :: jlev
INTEGER :: start_block

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NUM_OBS2'

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Go to start of obs file

CALL setpos (unit_no,0,icode)

!     Read in headers from obs file
! DEPENDS ON: readhead
CALL readhead (unit_no,                                   &
               fixhd,    len_fixhd,                       &
               inthd,    len_inthd,                       &
               realhd,   len_realhd,                      &
               levdepc,  len1_levdepc,len2_levdepc,       &
               rowdepc,  len1_rowdepc,len2_rowdepc,       &
               coldepc,  len1_coldepc,len2_coldepc,       &
               flddepc,  len1_flddepc,len2_flddepc,       &
               extcnst,  len_extcnst,                     &
               dumphist, len_dumphist,                    &
               cfi1,     len_cfi1,                        &
               cfi2,     len_cfi2,                        &
               cfi3,     len_cfi3,                        &
               lookup,   len1_lookup_obs,len2_lookup_obs, &
               len_data,                                  &
               start_block,icode,cmessage)
IF (icode >  0) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

obs_row_length = inthd(6)         !  No of points on row
obs_p_rows     = inthd(7)         !  No of rows
obs_p_levels   = inthd(8)         !  No of model levels

tnobs  = inthd(28)                !  Total no of observations
tndv   = inthd(29)                !  Total no of data values
nobtyp = inthd(32)                !  No of observation types

DO jlev=1,model_levels
  obs_ak(jlev) = levdepc(jlev+2,nobtyp+1) !  Vertical grid
  obs_bk(jlev) = levdepc(jlev+2,nobtyp+2) !
END DO

obs_long_res         = realhd(1)  !  Horizontal grid
obs_lat_res          = realhd(2)  !
obs_start_lat        = realhd(3)  !
obs_start_long       = realhd(4)  !
obs_lat_pseudo_pole  = realhd(5)  !
obs_long_pseudo_pole = realhd(6)  !

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE num_obs2
!-----------------------------------------------------------------------
END MODULE num_obs2_mod
