! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE get_sulpc_oxidants_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_SULPC_OXIDANTS_MOD'

CONTAINS

SUBROUTINE get_sulpc_oxidants(                                                 &
  ! Arguments IN
  l_sulpc_online_oxidants, l_ukca, l_ukca_trop, l_ukca_tropisop,               &
  l_ukca_strattrop, l_ukca_raq, l_ukca_raqaero, l_ukca_offline,                &
  theta, p_theta_levels, exner_theta_levels,                                   &
  o3_chem, h2o2_limit, oh, ho2,                                                &
  oh_ukca, ho2_ukca, h2o2_ukca, o3_ukca, hno3_ukca,                            &
  ! Arguments OUT
  o3_mmr_out, hno3_mmr_out, h2o2_mmr_out, oh_conc_out,                         &
  ho2_conc_out)

!---------------------------------------------------------------------
! Purpose: To extract from the D1 array the appropriate concentrations
!          of oxidants for the sulphur cycle, either prescribed
!          oxidants from ancillary file or on-line oxidants from UKCA,
!          depending on the value of L_SULPC_ONLINE_OXIDANTS.
!          Called by Atm_Step.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran
!  This code is written to UMDP3 programming standards

!---------------------------------------------------------------------

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s
USE ukca_constants, ONLY: c_oh, c_ho2
USE ereport_mod, ONLY: ereport
USE chemistry_constants_mod, ONLY: avogadro
USE c_rmol, ONLY: rmol

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Variables with intent IN:

! Logical variables
LOGICAL, INTENT(IN) :: l_sulpc_online_oxidants
LOGICAL, INTENT(IN) :: l_ukca
LOGICAL, INTENT(IN) :: l_ukca_trop
LOGICAL, INTENT(IN) :: l_ukca_tropisop
LOGICAL, INTENT(IN) :: l_ukca_strattrop
LOGICAL, INTENT(IN) :: l_ukca_raq
LOGICAL, INTENT(IN) :: l_ukca_raqaero
LOGICAL, INTENT(IN) :: l_ukca_offline

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

! Oxidant fields from anciliary file  - no halos
REAL, INTENT(IN) :: oh        (tdims%i_start:tdims%i_end,         &
                               tdims%j_start:tdims%j_end,         &
                               tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: ho2       (tdims%i_start:tdims%i_end,         &
                               tdims%j_start:tdims%j_end,         &
                               tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: h2o2_limit(tdims%i_start:tdims%i_end,         &
                               tdims%j_start:tdims%j_end,         &
                               tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: o3_chem   (tdims%i_start:tdims%i_end,         &
                               tdims%j_start:tdims%j_end,         &
                               tdims%k_start:tdims%k_end)

! UKCA oxidant fields on theta points
! No halo
REAL, INTENT(IN) :: oh_ukca (  tdims%i_start:tdims%i_end,      &
                               tdims%j_start:tdims%j_end,      &
                               tdims%k_start:tdims%k_end)
REAL, INTENT(IN)  :: ho2_ukca( tdims%i_start:tdims%i_end,      &
                               tdims%j_start:tdims%j_end,      &
                               tdims%k_start:tdims%k_end)
! Small halo
REAL, INTENT(IN)  :: h2o2_ukca(tdims_s%i_start:tdims_s%i_end,  &
                               tdims_s%j_start:tdims_s%j_end,  &
                               tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: o3_ukca  (tdims_s%i_start:tdims_s%i_end,  &
                               tdims_s%j_start:tdims_s%j_end,  &
                               tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN)  :: hno3_ukca(tdims_s%i_start:tdims_s%i_end,  &
                               tdims_s%j_start:tdims_s%j_end,  &
                               tdims_s%k_start:tdims_s%k_end)

! Variables with intent OUT:
! Oxidant arrays to be passed back to ATM_STEP have no haloes:
REAL, INTENT(OUT) :: o3_mmr_out(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,          &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(OUT) :: hno3_mmr_out(tdims%i_start:tdims%i_end,      &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(OUT) :: h2o2_mmr_out(tdims%i_start:tdims%i_end,      &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(OUT) :: oh_conc_out(tdims%i_start:tdims%i_end,       &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)
REAL, INTENT(OUT) :: ho2_conc_out(tdims%i_start:tdims%i_end,      &
                                tdims%j_start:tdims%j_end,        &
                                tdims%k_start:tdims%k_end)

! Local scalars:
REAL :: t
REAL :: air_density    ! molecules / cm3

INTEGER :: error_code             ! For passing to subroutine EREPORT
INTEGER :: i                      ! Loop variable
INTEGER :: j                      ! Loop variable
INTEGER :: k                      ! Loop variable
INTEGER :: level_off              ! Level offset depending on ND/EG


CHARACTER(LEN=errormessagelength) :: cmessage         ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_SULPC_OXIDANTS'

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that values of logical variables are consistent:
IF ((l_ukca_trop .OR. l_ukca_tropisop .OR. l_ukca_strattrop .OR.               &
  l_ukca_raq .OR. l_ukca_offline .OR. l_ukca_raqaero) .AND.                    &
  (.NOT. l_ukca)) THEN
  cmessage = 'Inconsistency in UKCA logicals.'
  error_code=100
  CALL ereport('GET_SULPC_OXIDANTS',error_code,cmessage)
END IF

IF (l_sulpc_online_oxidants .AND. (.NOT. (l_ukca_trop .OR.                     &
  l_ukca_tropisop .OR. l_ukca_strattrop .OR.                                   &
  l_ukca_raq .OR. l_ukca_offline .OR. l_ukca_raqaero))) THEN
  cmessage = 'UKCA oxidants are selected for sulphur cycle, '                  &
    //'but UKCA is turned off.'
  error_code=101
  CALL ereport('GET_SULPC_OXIDANTS',error_code,cmessage)
END IF

IF (l_sulpc_online_oxidants .AND. l_ukca_offline) THEN
  cmessage = 'UKCA oxidants are selected for sulphur cycle, '                  &
    //'but offline oxidants chemical scheme selected.'
  error_code=102
  CALL ereport('GET_SULPC_OXIDANTS',error_code,cmessage)
END IF

! There are two potential sources of oxidant fields: 
! 1. From UKCA
! 2. From an anciliary file
IF (l_sulpc_online_oxidants) THEN

  ! No level zero data set here as not used in CLASSIC 
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,t,air_density)  &
!$OMP SHARED( tdims,theta,exner_theta_levels,p_theta_levels,                   &
!$OMP oh_conc_out,ho2_conc_out,o3_mmr_out,hno3_mmr_out,h2o2_mmr_out,h2o2_ukca, &
!$OMP hno3_ukca,o3_ukca,ho2_ukca,oh_ukca)
  DO k=1,tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
      
        ! Calculate temperature from theta and exner
        t=theta(i,j,k)*exner_theta_levels(i,j,k)

        ! calculate molecular density from pressure and temperature
        air_density = p_theta_levels(i,j,k) *                        &
          avogadro * 1.0e-06 / (rmol * t)

        ! Convert OH and HO2 from mass-mixing ratios to molecules/cm3
        oh_conc_out(i,j,k) = oh_ukca(i,j,k) * air_density  / c_oh
        ho2_conc_out(i,j,k) = ho2_ukca(i,j,k) * air_density / c_ho2

        ! For the other variables just write to output array
        o3_mmr_out(i,j,k) = o3_ukca(i,j,k) 
        hno3_mmr_out(i,j,k) = hno3_ukca(i,j,k) 
        h2o2_mmr_out(i,j,k) = h2o2_ukca(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE

  ! Prescribed oxidants from ancillary selected. Set prescribed
  ! oxidants from array (OH, HO2, O3 and H2O2), or set to zero (HNO3).
  ! Again, don't set level 0. 
  
  ! Note that the offline oxidant arrays are indexed from 0 when using ENDGAME
  ! but the ancil files do not contian level 0 data.
  ! Calculate a level offset depending on the value of the first model 
  ! level - if EG then this is 0 and offset is -1. If ND then this is
  ! 1 and the offset is zero.
  level_off = tdims%k_start - 1
  
  ! Now loop over grid point from level 1 to top
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)              &
!$OMP SHARED( tdims,level_off,o3_mmr_out,h2o2_mmr_out,oh_conc_out,           &
!$OMP ho2_conc_out,o3_chem,h2o2_limit,oh,ho2,hno3_mmr_out)
  DO k=1,tdims%k_end
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
      
        ! O3
        o3_mmr_out(i,j,k)  = o3_chem(i,j,k + level_off)
        ! H2O2
        h2o2_mmr_out(i,j,k) = h2o2_limit(i,j,k + level_off)
        ! OH
        oh_conc_out(i,j,k)  =  oh(i,j,k + level_off)
        ! HO2
        ho2_conc_out(i,j,k) = ho2(i,j,k + level_off)
        ! HNO3
        hno3_mmr_out(i,j,k) = 0.0

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE get_sulpc_oxidants
END MODULE get_sulpc_oxidants_mod
