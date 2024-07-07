! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_cloudturb_pot_mod ---------------------------------------
!
!   Purpose: Calculates PWS cloud turb potential, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_cloudturb_pot_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_CLOUDTURB_POT_MOD'

CONTAINS

SUBROUTINE pws_cloudturb_pot(p_theta_levels,theta,exner_theta_levels,cf_bulk,&
                              press_lev,klev)

! Description
!   1) Loop through the gridpoints and for each find the model level
!      immediately above and below the pressure level of interest.
!   2) Construct a cloud mask array indicating if there is cloud on the
!      model level above or below each grid point.
!   3) Find vertical derivatives of height and equivalent potential temp
!      w.r.t. pressure using cubic spline. Then calculate the vertical
!      derviative of equivalent potential temperature with height.
!   4) Then, for each gridpoint, set the in cloud turbulence predictor
!      If orography makes calculation impossible, set
!      cld_turb_Pred to missing data indicator. If calculation is possible,
!      check that the cloud mask=1. If so, check sign of dthetaedZ. If
!      negative mark predictor as turbulent (set to dthetadz); if positive
!      mark as non-turbulent (zero).

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, tdims_l, pdims
USE level_heights_mod,     ONLY: r_theta_levels
USE planet_constants_mod,  ONLY: planet_radius, g_over_r, lapse_trop, &
                                 cp, rv, r, repsilon
USE conversions_mod,       ONLY: zerodegc
USE water_constants_mod,   ONLY: lc
USE missing_data_mod,      ONLY: imdi, rmdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels

USE pws_diags_mod, ONLY: pws_cloudturb_pot_diag, flag_cloudturb_pot

USE diff_mod, ONLY: diffp

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


REAL, INTENT(IN)  :: theta(tdims_s%i_start:tdims_s%i_end,                    &
                           tdims_s%j_start:tdims_s%j_end,                    &
                           tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,       &
                                        tdims_s%j_start:tdims_s%j_end,       &
                                        tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: cf_bulk(tdims_l%i_start:tdims_l%i_end,                  &
                             tdims_l%j_start:tdims_l%j_end,                  &
                             tdims_l%k_start:tdims_l%k_end)

REAL, INTENT(IN) :: press_lev
INTEGER, INTENT(IN) :: klev

! Local variables
REAL    :: temp(row_length,rows,model_levels) ! temperature
REAL    :: EqPotTemp(row_length,rows,model_levels) ! equiv pot temperature


REAL, PARAMETER :: eso = 611.0    ! Sat. vapour pressure at reference temp

REAL            :: x             ! exponiential part of saturated vapour
                                 !   pressure equation
REAL            :: es            ! saturated vapour pressure in Pa
REAL            :: ws            ! saturated mixing ratio in
REAL            :: theta_e       ! Equivalent potential temperature in K


! Local Variables:
INTEGER              :: i, j, k           ! Loop counters
INTEGER              :: jlwr, jupr, jmid  ! Variables used to find model
                                          !  levels above and below level
                                          !  of interest
INTEGER :: Lower (row_length,rows)       ! Lower model level
INTEGER :: Upper (row_length,rows)       ! Upper model level
INTEGER :: MD_flag (row_length,rows)     ! Missing data flag (used when
                                          !   level is beneath first model
                                          !   level)
INTEGER :: cld_mask (row_length,rows)    ! Cloud mask

REAL :: dZdP(row_length,rows)           ! Vert. Deriv of height wrt p
REAL :: dthetaedP(row_length,rows)      ! Vert. Deriv of equivalent pot.
                                          !   temp wrt p
REAL :: dthetaedZ(row_length,rows)      ! Vert. Deriv of equivalent pot.
                                          !   temp wrt p


CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_CLOUDTURB_POT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -----------------------------------

! calculate the temperature.
DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      temp(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
    END DO 
  END DO 
END DO 


! Calculate saturated vapour pressure,
! saturated mixing ratio then equivalent potential temperature for each
! gridpoint & level where possible. 

! Set equivalent potential temperature
! to missing data indictor if the terrain makes calculation impossible.

DO k = 1, model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end

      !Calculate the saturated vapour pressure for this point
      x = (lc/rv) * ( ( 1/zerodegc ) - ( 1/temp(i,j,k) ) )
      es = eso * EXP(x)

      !Calculate the saturated mixing ratio at this point
      ws = repsilon * ( es / ( p_theta_levels(i,j,k)-es ) )

      !Calculate the equivalent potential temperature at this point
      !and put into the output field element

      EqPotTemp(i,j,k)= theta(i,j,k) * EXP( (lc*ws)/(cp*temp(i,j,k)) )

    END DO
  END DO
END DO

!------------------------------------------------------
!Find the model level immediately above and below each gridpointd 


DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    IF (press_lev >= p_theta_levels(i,j,1)) THEN
      Lower(i,j)=1
      Upper(i,j)=2
    ELSE IF ( press_lev < p_theta_levels(i,j,model_levels)) THEN
      Lower(i,j)=model_levels-1
      Upper(i,j)=model_levels
    ELSE
      Lower(i,j)=-1 ! one level to find
    END IF
  END DO
END DO

DO k=2,model_levels
  DO j = tdims%j_start,tdims%j_end
    DO i = tdims%i_start,tdims%i_end
      IF ((p_theta_levels(i,j,k) <= press_lev) .AND. &
          (Lower(i,j) == -1)) THEN
        Lower(i,j)=k-1
        Upper(i,j)=k
      END IF
    END DO
  END DO
END DO

DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    IF ( p_theta_levels(i,j,Lower(i,j)) < press_lev &
         .OR. Lower(i,j) < 3 ) THEN
      MD_flag(i,j) = 1
    ELSE
      MD_flag(i,j) = 0
    END IF
  END DO
END DO

!------------------------------------------------------
!Construct cloud mask. Do this by checking that, for each gridpoint, there
! either the model level above or below this pressure level; if so mark cld
! if not mark as 0.

DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end

    IF (cf_bulk(i,j,Lower(i,j)) == rmdi .OR. &
        cf_bulk(i,j,Upper(i,j)) == rmdi .OR. &
         MD_flag(i,j) == 1 ) THEN
      cld_mask(i,j) = imdi
    ELSE IF ( cf_bulk(i,j,Lower (i,j)) > 0.0 .OR.         &
              cf_bulk( i,j,Upper (i,j)) > 0.0 ) THEN
      cld_mask(i,j) = 1
    ELSE
      cld_mask(i,j) = 0
    END IF

  END DO
END DO


!----------------------------------------------------------
! Find vertical derivatives of z and equivalent pot temp wrt pressure using
! cubic spline. Note that spline values must be monotonically increasing.
! Pressure level supplied to DiffP must be in Pa.

CALL DiffP(tdims%k_end, press_lev,                    &
           tdims%i_start,tdims%i_end,                 &
           tdims%j_start,tdims%j_end,                 &
           r_theta_levels(tdims%i_start:tdims%i_end,  &
                          tdims%j_start:tdims%j_end,  &
                          1:tdims%k_end),             &
           p_theta_levels(tdims%i_start:tdims%i_end,  &
                          tdims%j_start:tdims%j_end,  &
                          1:tdims%k_end), dZdP)

CALL DiffP(tdims%k_end, press_lev,                    &
           tdims%i_start,tdims%i_end,                 &
           tdims%j_start,tdims%j_end,                 &
           EqPotTemp(:,:,1:tdims%k_end),              &
           p_theta_levels(tdims%i_start:tdims%i_end,  &
                          tdims%j_start:tdims%j_end,  &
                          1:tdims%k_end), dthetaedP)


WHERE (  cld_mask /=  imdi )
  dthetaedZ(:,:) = dthetaedP(:,:)/ dZdP(:,:)
ELSE WHERE
  dthetaedZ(:,:) = rmdi
END WHERE


!--------------------------------------------------------------
!Loop over all gridpoints at this level & set cloud turbulence predictor.
! If orography makes calculation impossible, set cld_turb_Pred
! to missing data indicator. If calculation is possible,
! check cloud mask=1. If so, check sign of dthetaedZ. If negative,
! mark as turbulent (set to dthetadz); if positive mark as
! non-turbulent (zero).

DO j =  tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end

    IF (cld_mask(i,j) == imdi ) THEN
      pws_cloudturb_pot_diag(i,j,klev) = rmdi
    ELSE IF (cld_mask (i,j) == 1 .AND. &
             dthetaedZ(i,j) < 0.0 ) THEN
      pws_cloudturb_pot_diag(i,j,klev) = ABS(dthetaedZ(i,j))
    ELSE
      pws_cloudturb_pot_diag(i,j,klev) = 0.0
    END IF

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_cloudturb_pot
END MODULE pws_cloudturb_pot_mod
