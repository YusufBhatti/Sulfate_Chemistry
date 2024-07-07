! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To obtain GLOMAP-MODE input to UKCA_RADAER from d1 and Stashwork
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_radaer_get_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName ='GLOMAP_CLIM_RADAER_GET_MOD'

CONTAINS

SUBROUTINE glomap_clim_radaer_get( ierr, cmessage, ukca_radaer, stashwork54 )

USE atm_fields_real_mod,    ONLY: &
    gc_nd_nuc_sol,                &
    gc_nuc_sol_su,                &
    gc_nuc_sol_oc,                &
    gc_nd_ait_sol,                &
    gc_ait_sol_su,                &
    gc_ait_sol_bc,                &
    gc_ait_sol_oc,                &
    gc_nd_acc_sol,                &
    gc_acc_sol_su,                &
    gc_acc_sol_bc,                &
    gc_acc_sol_oc,                &
    gc_acc_sol_ss,                &
    gc_nd_cor_sol,                &
    gc_cor_sol_su,                &
    gc_cor_sol_bc,                &
    gc_cor_sol_oc,                &
    gc_cor_sol_ss,                &
    gc_nd_ait_ins,                &
    gc_ait_ins_bc,                &
    gc_ait_ins_oc

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE glomap_clim_fields_mod, ONLY: &
    prepare_fields_for_radaer

USE glomap_clim_option_mod, ONLY: &
    i_glomap_clim_setup,          &
    i_gc_sussocbc_5mode

USE nlsizes_namelist_mod,   ONLY: &
    model_levels

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE stash_array_mod,        ONLY: &
    stash_maxlen

USE submodel_mod,           ONLY: &
    atmos_im

USE ukca_radaer_struct_mod, ONLY: &
    ukca_radaer_struct

USE um_stashcode_mod,       ONLY: &
    stashcode_glomap_clim_sec

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook


IMPLICIT NONE

! Arguments

! Error indicator (0 is OK, /=0 error)
INTEGER, INTENT(OUT) :: ierr      ! error code

! Error message if ierr is larger than 0
CHARACTER (LEN=errormessagelength), INTENT(OUT) :: cmessage

! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer

! stashwork54 array
REAL, INTENT(INOUT) :: stashwork54(stash_maxlen(stashcode_glomap_clim_sec,     &
                                                atmos_im))

! Local variables

INTEGER :: i, j ! loop

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GLOMAP_CLIM_RADAER_GET'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Initialise error code
ierr = 0

! Populate ukca_radaer %dry_diam %wet_diam %modal_rho %modal_wtv %comp_vol
! Return stashwork54
CALL prepare_fields_for_radaer( ukca_radaer, stashwork54 )

! Update ukca_radaer %wet_diam %modal_wtv %modal_vol
DO j = 1, ukca_radaer%n_mode
  IF (.NOT. ukca_radaer%l_soluble(j)) THEN
    ! Insoluble modes: copy the dry diameter into the wet
    ! diameter, for consistency, and set the volume of water to zero.
    ukca_radaer%wet_diam(:,:,:,j)  = ukca_radaer%dry_diam(:,:,:,j)
    ukca_radaer%modal_wtv(:,:,:,j) = 0.0e+00
  END IF
  
  ! Initialise modal volume to that of water (the latter has been
  ! set to zero for insoluble modes above).
  ! The other components will be added in the loop below.
  ukca_radaer%modal_vol(:,:,:,j) = ukca_radaer%modal_wtv(:,:,:,j)
END DO

! Update ukca_radaer%modal_vol of each mode by adding the volume
! of each component within that mode.
DO j = 1, ukca_radaer%n_mode
  DO i = 1, ukca_radaer%n_cpnt_in_mode(j)
    ukca_radaer%modal_vol(:,:,:,j) = ukca_radaer%modal_vol(:,:,:,j) +          &
                     ukca_radaer%comp_vol(:,:,:,ukca_radaer%i_cpnt_index(i, j))
  END DO
END DO

! Populate ukca_radaer %modal_nbr %mix_ratio
SELECT CASE(i_glomap_clim_setup)
CASE(i_gc_sussocbc_5mode)
  IF ( (ukca_radaer%n_mode /= 4) .OR. (ukca_radaer%n_cpnt /= 13) ) THEN
    ierr = (ukca_radaer%n_mode * 100) + ukca_radaer%n_cpnt
    cmessage = 'glomap_clim_radaer_get : number components/modes does not match'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,           &
                                                          zhook_handle)
    RETURN
  END IF
  ! There are 4 modes (note RADAER neglects interaction with nucleation mode)
  ukca_radaer%modal_nbr(:,:,:,1)  = gc_nd_ait_sol(:,:,1:model_levels)
  ukca_radaer%modal_nbr(:,:,:,2)  = gc_nd_acc_sol(:,:,1:model_levels)
  ukca_radaer%modal_nbr(:,:,:,3)  = gc_nd_cor_sol(:,:,1:model_levels)
  ukca_radaer%modal_nbr(:,:,:,4)  = gc_nd_ait_ins(:,:,1:model_levels)
  ! There are 13 components
  ukca_radaer%mix_ratio(:,:,:,1)  = gc_ait_sol_su(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,2)  = gc_ait_sol_bc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,3)  = gc_ait_sol_oc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,4)  = gc_acc_sol_su(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,5)  = gc_acc_sol_bc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,6)  = gc_acc_sol_oc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,7)  = gc_acc_sol_ss(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,8)  = gc_cor_sol_su(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,9)  = gc_cor_sol_bc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,10) = gc_cor_sol_oc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,11) = gc_cor_sol_ss(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,12) = gc_ait_ins_bc(:,:,1:model_levels)
  ukca_radaer%mix_ratio(:,:,:,13) = gc_ait_ins_oc(:,:,1:model_levels)
CASE DEFAULT
  ierr = 54
  cmessage = 'glomap_clim_radaer_get : Unrecognised glomap_clim setup integer'
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
  RETURN
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE glomap_clim_radaer_get

END MODULE glomap_clim_radaer_get_mod
