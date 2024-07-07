! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE disturb_veg_category_mod

IMPLICIT NONE

! Pointer used to select which type of disturbed veg fraction
! should be passed to JULES routines
REAL, POINTER, SAVE :: disturb_veg_pointer(:)
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DISTURB_VEG_CATEGORY_MOD'
PRIVATE
PUBLIC :: set_disturb_veg_category, disturb_veg_pointer

CONTAINS
 
SUBROUTINE set_disturb_veg_category()
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Vegetation
!
! Code Description:
!  Language: FORTRAN 2003.
!  This code is written to UMDP3 programming standards.
USE atm_fields_real_mod, ONLY: disturb_veg, agr_crop_frac_d1
USE jules_vegetation_mod, ONLY: l_trif_crop
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='SET_DISTURB_VEG_CATEGORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (l_trif_crop) THEN
  ! Get the crop fraction of disturbed veg
  disturb_veg_pointer => agr_crop_frac_d1
ELSE
  ! Get total disturbed fraction of veg
  disturb_veg_pointer => disturb_veg
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_disturb_veg_category
END MODULE disturb_veg_category_mod
