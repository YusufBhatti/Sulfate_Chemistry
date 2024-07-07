! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE write_sulpc_oxidants_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WRITE_SULPC_OXIDANTS_MOD'

CONTAINS

SUBROUTINE write_sulpc_oxidants(                                               &
  oh_conc_in, h2o2_mmr_in, ho2_conc_in, o3_mmr_in, hno3_mmr_in,                &
  theta, p_theta_levels, exner_theta_levels,                                   &
  oh_ukca, ho2_ukca, h2o2_ukca, o3_ukca, hno3_ukca)

!---------------------------------------------------------------------
! Purpose: To write to the ukca arrays the concentrations of oxidants
!          after depletion in the sulphur cycle.
!          Called by Atm_Step only if L_SULPC_ONLINE_OXIDANTS is set
!          to .TRUE.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code Description:
!  Language: Fortran
!  This code is written to UMDP3 programming standards
!
!---------------------------------------------------------------------

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s
USE ukca_constants, ONLY: c_oh, c_ho2
USE chemistry_constants_mod, ONLY: avogadro
USE c_rmol, ONLY: rmol

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WRITE_SULPC_OXIDANTS'

! Variables with intent IN:
! No halo on theta points
REAL, INTENT(IN) :: oh_conc_in(tdims%i_start:tdims%i_end,         &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: h2o2_mmr_in(tdims%i_start:tdims%i_end,        &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: ho2_conc_in(tdims%i_start:tdims%i_end,        &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: o3_mmr_in(tdims%i_start:tdims%i_end,          &
                              tdims%j_start:tdims%j_end,          &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) :: hno3_mmr_in(tdims%i_start:tdims%i_end,        &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)

! Small halo on theta points
! Potential temperature
REAL, INTENT(IN) :: theta(tdims_s%i_start:tdims_s%i_end,           &
                          tdims_s%j_start:tdims_s%j_end,           &
                          tdims_s%k_start:tdims_s%k_end)

! Pressure on theta levels
REAL, INTENT(IN) :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,  &
                                   tdims_s%j_start:tdims_s%j_end,  &
                                   tdims_s%k_start:tdims_s%k_end)

! Exner pressure on theta levels
REAL, INTENT(IN) :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end, &
                                       tdims_s%j_start:tdims_s%j_end, &
                                       tdims_s%k_start:tdims_s%k_end)

! Local arrays:

! Arrays with intent INOUT:
! UKCA oxidant fields on theta points
! No halo
REAL, INTENT(INOUT) :: oh_ukca (  tdims%i_start:tdims%i_end,      &
                                  tdims%j_start:tdims%j_end,      &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT(INOUT)  :: ho2_ukca( tdims%i_start:tdims%i_end,      &
                                  tdims%j_start:tdims%j_end,      &
                                  tdims%k_start:tdims%k_end)
! Small halo
REAL, INTENT(INOUT)  :: h2o2_ukca(tdims_s%i_start:tdims_s%i_end,  &
                                  tdims_s%j_start:tdims_s%j_end,  &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT)  :: o3_ukca  (tdims_s%i_start:tdims_s%i_end,  &
                                  tdims_s%j_start:tdims_s%j_end,  &
                                  tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT)  :: hno3_ukca(tdims_s%i_start:tdims_s%i_end,  &
                                  tdims_s%j_start:tdims_s%j_end,  &
                                  tdims_s%k_start:tdims_s%k_end)
! Local scalars:
REAL :: t
REAL :: air_density

INTEGER :: i          ! Loop variable
INTEGER :: j          ! Loop variable
INTEGER :: k          ! Loop variable


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------

! For ENDGAME level 0 is never defined in the input arrays,
! so loop from 1.
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,t,air_density) &
!$OMP SHARED( tdims,theta,exner_theta_levels,p_theta_levels,                  &
!$OMP oh_ukca,ho2_ukca,o3_ukca,hno3_ukca,h2o2_ukca,h2o2_mmr_in,hno3_mmr_in,   &
!$OMP o3_mmr_in,ho2_conc_in,oh_conc_in)
DO k=1,tdims%k_end
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end

      ! Calculate temperature from theta and exner
      t = theta(i,j,k)*exner_theta_levels(i,j,k)

      ! Calculate molecular density from pressure and temperature
      air_density = p_theta_levels(i,j,k) * avogadro                 &
        * 1.0e-06 / (rmol * t)
      
      ! Convert OH and HO2 from molecules/cm3 to mass-mixing ratios
      
      oh_ukca(i,j,k) = oh_conc_in(i,j,k) * c_oh / air_density
      ho2_ukca(i,j,k) = ho2_conc_in(i,j,k) * c_ho2 / air_density

      ! For the other variables just write to output array
      o3_ukca(i,j,k) = o3_mmr_in(i,j,k)
      hno3_ukca(i,j,k) = hno3_mmr_in(i,j,k)
      h2o2_ukca(i,j,k) = h2o2_mmr_in(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE write_sulpc_oxidants
END MODULE write_sulpc_oxidants_mod
