! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Top level routine for setting up common GRIB between original GRIB and
!  fieldsfiles created with GRIB_API.

MODULE rcf_grib2ff_init_mod

! MODULE Rcf_Grib2ff_Init_Mod

! Description: This sets up required information to handle translation from
! fieldsfile created from GRIB_API.

! Method: Initialise module data which is passed down into rcf_generate_heights
! and rcf_grib_interp_tnpstar

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB2FF_INIT_MOD'

CONTAINS

SUBROUTINE rcf_grib2ff_init( hdr_in, grid )

USE rcf_grib_t_n_pstar_h_interp_mod, ONLY:                                     &
  ak, bk

USE rcf_headaddress_mod, ONLY:                                                 &
  ldc_ak_hybrid,                                                               &
  ldc_bk_hybrid

USE rcf_umhead_mod, ONLY:                                                     &
  um_header_type

USE rcf_grid_type_mod, ONLY:                                                  &
  grid_type

USE rcf_generate_heights_mod, ONLY:                                           &
  height_gen_ecmwf_hybrd

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

TYPE (um_header_type), INTENT(IN)    :: hdr_in
TYPE( grid_type ), INTENT(INOUT)   :: grid

INTEGER :: i
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_GRIB2FF_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Ak may have been allocated by old grib conversion section in reconfiguration.
IF (grid % height_gen_method == height_gen_ecmwf_hybrd .AND. &
    .NOT. ALLOCATED(ak)) THEN
  ! Read in ak and bk values from fieldsfile header.
  ! Assuming we have mlindex in levdepc already - in this case just
  ! probably be a sequence from 1 to model_levels.
  ALLOCATE(ak(hdr_in % len1levdepc-1))
  ALLOCATE(bk(hdr_in % len1levdepc-1))

  ! We actually want full levels so we need to find the middle of adjacent
  ! half levels.
  DO i = 1, hdr_in % len1levdepc-1
    ak(i) = 0.5*(hdr_in % levdepc(i,ldc_ak_hybrid)   + &
                 hdr_in % levdepc(i+1,ldc_ak_hybrid))
    bk(i) = 0.5*(hdr_in % levdepc(i,ldc_bk_hybrid)   + &
                 hdr_in % levdepc(i+1,ldc_bk_hybrid))
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_grib2ff_init

END MODULE rcf_grib2ff_init_mod
