! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_updh_diag ---------------------------------------
!
!   Purpose: Calculates updraught helicity for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No ??
!   Also see: http://journals.ametsoc.org/doi/abs/10.1175/WAF2007106.1 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_updh_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_UPDH_DIAG_MOD'

CONTAINS

SUBROUTINE pws_updh_diag(pws_upd_helicity_5k, icode, cmessage)

USE halo_exchange, ONLY: swap_bounds
USE level_heights_mod, ONLY: r_theta_levels
USE planet_constants_mod, ONLY: planet_radius
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE trignometric_mod, ONLY: cos_theta_latitude

USE nlsizes_namelist_mod, ONLY: model_levels
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s
USE pws_diags_mod, ONLY: pws_pcd_mlev_u, pws_pcd_mlev_v, pws_pcd_mlev_w
USE um_parvars, ONLY: offx, offy
USE Field_Types, ONLY: fld_type_p
USE mpp_conf_mod, ONLY: swap_field_is_scalar
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

! Arguments
REAL, INTENT(INOUT) :: pws_upd_helicity_5k (tdims%i_start:tdims%i_end,  &
                                            tdims%j_start:tdims%j_end)

INTEGER :: icode ! Return code : 0 Normal exit, >0 Error exit
CHARACTER(LEN=errormessagelength) :: cmessage ! Error message

! Local variables
INTEGER :: i, j, k
INTEGER :: levs                  ! number of levels in halo exchange
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_UPDH_DIAG'

! Layer dimensions in metres.
REAL, PARAMETER :: top_of_layer = 5000.0
REAL, PARAMETER :: bottom_of_layer = 2000.0

! thickness of updraught layer
REAL :: thickness (tdims%i_len, tdims%j_len)
! level by level helicity
REAL :: upd_helicity (tdims%i_len, tdims%j_len)
REAL :: DvDx (tdims%i_len, tdims%j_len) ! latitudinal derivative of u
REAL :: DuDy (tdims%i_len, tdims%j_len) ! longitudinal derivative of v
! wind components with halos
REAL :: uhalo(tdims_s%i_start:tdims_s%i_end,               &
              tdims_s%j_start:tdims_s%j_end,2:model_levels-1)
REAL :: vhalo(tdims_s%i_start:tdims_s%i_end,               &
              tdims_s%j_start:tdims_s%j_end,2:model_levels-1)
REAL :: HdelX  ! 1/2*W-E grid spacing
REAL :: HdelY  ! 1/2*S-N grid spacing

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Zero some key arrays
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    thickness           (i,j) = 0.0
    upd_helicity        (i,j) = 0.0
    pws_upd_helicity_5k (i,j) = 0.0
  END DO
END DO

! compute updraught helicity for 2000-5000 m range from wind components 
! using local copies of diagnostic stash code 30,001/2/3 (units m/s)

! Note that r_theta_levels(:,:,k+1), r_theta_levels(:,:,k) 
! are upper and lower layer boundaries for r_rho_levels(:,:,k)

DO k = 2, model_levels-1

! Calculate relative vorticity on model levels

! Copy horizontal winds into arrays with haloes
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      uhalo(i,j,k) = pws_pcd_mlev_u(i,j,k)
      vhalo(i,j,k) = pws_pcd_mlev_v(i,j,k)
    END DO
  END DO
END DO

! Call Swapbounds to fill in haloes ready for taking derivatives
levs = model_levels - 2 ! i.e. from 2 to model_levels-1
CALL swap_bounds(uhalo,tdims%i_len,tdims%j_len,levs,                &
                 offx,offy,fld_type_p,swap_field_is_scalar)

CALL swap_bounds(vhalo,tdims%i_len,tdims%j_len,levs,                &
                 offx,offy,fld_type_p,swap_field_is_scalar)

DO k = 2, model_levels-1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
! Set thickness to zero unless level within required range, 2-5km.
      IF ( (r_theta_levels(i,j,k+1)                               &
          - r_theta_levels(i,j,1)) < bottom_of_layer ) THEN
        thickness (i,j) = 0.0
      ELSE IF ( (r_theta_levels(i,j,k)                            &
          - r_theta_levels(i,j,1)) > top_of_layer) THEN
        thickness (i,j) = 0.0
      ELSE
        thickness (i,j) = r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k)
        HdelX = 0.5/(planet_radius * delta_lambda * cos_theta_latitude(i,j))
        HdelY = 0.5/(planet_radius * delta_phi * cos_theta_latitude(i,j))
        DvDx(i,j) = HdelX * ( vhalo(i+1,j,k) - vhalo(i-1,j,k) )
        DuDy(i,j) = HdelY * (uhalo(i,j+1,k)*cos_theta_latitude(i,j+1) &
                           - uhalo(i,j-1,k)*cos_theta_latitude(i,j-1))
        upd_helicity(i,j) = (DvDx(i,j) - DuDy(i,j)) * pws_pcd_mlev_w(i,j,k)
! Option 1: Zero this level UPdraft helicity if vertical velocity is -ve
        IF (pws_pcd_mlev_w(i,j,k) < 0.0) THEN
          upd_helicity(i,j) = 0.0
        END IF
! Option 2: Zero this level UPdraft helicity if -ve in NH or +ve in SH
! N.B. true_latitude will need to be USEd from trignometric_mod
!        IF (upd_helicity(i,j)*true_latitude(i,j) < 0.0) THEN
!          upd_helicity(i,j) = 0.0
!        END IF
      END IF
    END DO
  END DO

! Sum updraft helicities over the required vertical range
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      pws_upd_helicity_5k(i,j) = pws_upd_helicity_5k(i,j) + &
                                 thickness(i,j) * upd_helicity(i,j)
    END DO
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_updh_diag

END MODULE pws_updh_diag_mod
