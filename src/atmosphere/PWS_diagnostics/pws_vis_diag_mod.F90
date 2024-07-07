! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutines pws_vis_diag, pws_vis2_diag ---------------------------
!
!   Purpose: Calculates visibilities for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No ??
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_vis_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_VIS_DIAG_MOD'

CONTAINS

! --------------------------------------------------------------
! Vis based on 6 bin dust

SUBROUTINE pws_vis_diag(pws_bl_1p5m_vis_tot, pws_bl_1p5m_temp, &
                        pws_1p5m_vis_dust, pws_1p5m_vis_tot,   &
                        flag_1p5m_vis_tot, flag_1p5m_vis_dust)

USE atm_fields_real_mod, ONLY: pstar, dust_div1, dust_div2, dust_div3, &
                                      dust_div4, dust_div5, dust_div6

USE visbty_constants_mod, ONLY: LnLiminalContrast, RecipVisAir
USE planet_constants_mod, ONLY: r  

USE nlsizes_namelist_mod, ONLY: model_levels
USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE dust_parameters_mod, ONLY: pwsdiag_sfc_em
USE pws_diags_mod, ONLY: pws_dustmmr1_em, pws_dustmmr2_em,              &
                         pws_dustmmr3_em, pws_dustmmr4_em,              &
                         pws_dustmmr5_em, pws_dustmmr6_em,              &
                         flag_dustmmr_em

IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: pws_bl_1p5m_vis_tot                                       &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: pws_bl_1p5m_temp                                          &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(INOUT) :: pws_1p5m_vis_tot                                       &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(INOUT) :: pws_1p5m_vis_dust                                      &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
LOGICAL, INTENT(IN) :: flag_1p5m_vis_tot
LOGICAL, INTENT(IN) :: flag_1p5m_vis_dust

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_VIS_DIAG'

REAL :: beta_air      ! extinction in clean air
REAL :: rho           ! air density
REAL :: k_ext(6)      ! specific extinction coeffs at 550nm for each bin
DATA k_ext(1), k_ext(2), k_ext(3), k_ext(4) , k_ext(5), k_ext(6) &
                     / 652.797, 3626.70, 979.290, 260.675, 77.6338, 23.9075 /

! dust mmr to use to calculate extinctions:
REAL :: dust_wrk1(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk3(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk4(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk5(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk6(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
! extinction due to dust
REAL :: beta_dust(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
! running total of the extinction
REAL :: beta_tot (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)

INTEGER :: i,j,k

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! calculate the dust mmr to use as inputs:
k = 1
IF (flag_dustmmr_em) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dust_wrk1(i,j) = (pws_dustmmr1_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div1(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk2(i,j) = (pws_dustmmr2_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div2(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk3(i,j) = (pws_dustmmr3_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div3(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk4(i,j) = (pws_dustmmr4_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div4(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk5(i,j) = (pws_dustmmr5_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div5(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk6(i,j) = (pws_dustmmr6_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div6(i,j,k) * (1.0-pwsdiag_sfc_em))
    END DO
  END DO
ELSE
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dust_wrk1(i,j) = dust_div1(i,j,k)
      dust_wrk2(i,j) = dust_div2(i,j,k)
      dust_wrk3(i,j) = dust_div3(i,j,k)
      dust_wrk4(i,j) = dust_div4(i,j,k)
      dust_wrk5(i,j) = dust_div5(i,j,k)
      dust_wrk6(i,j) = dust_div6(i,j,k)
    END DO
  END DO
END IF

! Calculate the extinction in clean air from RecipVisAir
beta_air = -LnLiminalContrast * RecipVisAir

! Calculate the extinction due to dust
beta_dust(:,:)=0.0
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    rho = pstar(i,j) / (pws_bl_1p5m_temp(i,j) * r)

    beta_dust(i,j)= rho * (             dust_wrk1(i,j)*k_ext(1) + &
              dust_wrk2(i,j)*k_ext(2) + dust_wrk3(i,j)*k_ext(3) + &
              dust_wrk4(i,j)*k_ext(4) + dust_wrk5(i,j)*k_ext(5) + &
                                        dust_wrk6(i,j)*k_ext(6))
  END DO
END DO

! Calculate visibility fields
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    !   Invert visibility to calculate total extinction
    beta_tot(i,j) = -LnLiminalContrast / pws_bl_1p5m_vis_tot(i,j)

    !   Add the extinction from dust to the total
    beta_tot(i,j) = beta_tot(i,j) + beta_dust(i,j)

  END DO
END DO

IF (flag_1p5m_vis_tot) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !   Invert back to get visibilities
      pws_1p5m_vis_tot(i,j)  = -LnLiminalContrast / beta_tot(i,j)
    END DO
  END DO
END IF
IF (flag_1p5m_vis_dust) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !   Include small contribution from air to limit vis model's max value
      pws_1p5m_vis_dust(i,j) = -LnLiminalContrast / (beta_dust(i,j) + beta_air)
    END DO
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_vis_diag

! ---------------------------------------------------------------
! Vis based on 2 bin dust

SUBROUTINE pws_vis2_diag(pws_bl_1p5m_vis_tot, pws_bl_1p5m_temp, &
                         pws_1p5m_vis_dust, pws_1p5m_vis_tot,   &
                         flag_1p5m_vis_tot, flag_1p5m_vis_dust)

USE atm_fields_real_mod, ONLY: pstar, dust_div1, dust_div2

USE visbty_constants_mod, ONLY: LnLiminalContrast, RecipVisAir
USE planet_constants_mod, ONLY: r  

USE nlsizes_namelist_mod, ONLY: model_levels
USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE dust_parameters_mod, ONLY: pwsdiag_sfc_em
USE pws_diags_mod, ONLY: pws_dustmmr1_em, pws_dustmmr2_em,              &
                         pws_dustmmr3_em, pws_dustmmr4_em,              &
                         pws_dustmmr5_em, pws_dustmmr6_em,              &
                         flag_dustmmr_em

IMPLICIT NONE

! Arguments
REAL, INTENT(IN) :: pws_bl_1p5m_vis_tot                                       &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: pws_bl_1p5m_temp                                          &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(INOUT) :: pws_1p5m_vis_tot                                       &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL, INTENT(INOUT) :: pws_1p5m_vis_dust                                      &
                   (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
LOGICAL, INTENT(IN) :: flag_1p5m_vis_tot
LOGICAL, INTENT(IN) :: flag_1p5m_vis_dust
! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_VIS2_DIAG'

REAL :: beta_air      ! extinction in clean air
REAL :: rho           ! air density
REAL :: k_ext(2)      ! specific extinction coeffs at 550nm for each bin
DATA k_ext(1), k_ext(2) / 700.367, 141.453 /
! dust mmr to use to calculate extinctions:
REAL :: dust_wrk1(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL :: dust_wrk2(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
! extinction due to dust
REAL :: beta_dust(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
! running total of the extinction
REAL :: beta_tot (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)

INTEGER :: i,j,k

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the extinction in clean air from RecipVisAir
beta_air = -LnLiminalContrast * RecipVisAir

! calculate the dust mmr to use as inputs:
k = 1
IF (flag_dustmmr_em) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dust_wrk1(i,j) = (pws_dustmmr1_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div1(i,j,k) * (1.0-pwsdiag_sfc_em))
      dust_wrk2(i,j) = (pws_dustmmr2_em(i,j) * pwsdiag_sfc_em) + &
                     (dust_div2(i,j,k) * (1.0-pwsdiag_sfc_em))
    END DO
  END DO
ELSE
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dust_wrk1(i,j) = dust_div1(i,j,k)
      dust_wrk2(i,j) = dust_div2(i,j,k)
    END DO
  END DO
END IF

! Calculate the extinction due to dust
beta_dust(:,:)=0.0
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    rho = pstar(i,j) / (pws_bl_1p5m_temp(i,j) * r)

    beta_dust(i,j)= rho * ( dust_wrk1(i,j)*k_ext(1) + &
                            dust_wrk2(i,j)*k_ext(2) )
  END DO
END DO

! Calculate visibility fields
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end

    !   Invert visibility to calculate total extinction
    beta_tot(i,j) = -LnLiminalContrast / pws_bl_1p5m_vis_tot(i,j)

    !   Add the extinction from dust to the total
    beta_tot(i,j) = beta_tot(i,j) + beta_dust(i,j)

  END DO
END DO

IF (flag_1p5m_vis_tot) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !   Invert back to get visibilities
      pws_1p5m_vis_tot(i,j)  = -LnLiminalContrast / beta_tot(i,j)
    END DO
  END DO
END IF
IF (flag_1p5m_vis_dust) THEN
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !   Include small contribution from air to limit vis model's max value
      pws_1p5m_vis_dust(i,j) = -LnLiminalContrast / (beta_dust(i,j) + beta_air)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_vis2_diag

END MODULE pws_vis_diag_mod
