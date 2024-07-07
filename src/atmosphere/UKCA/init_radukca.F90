! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   initialization of radiative feedback
!
! Subroutine Interface:
MODULE init_radukca_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INIT_RADUKCA_MOD'

CONTAINS 

SUBROUTINE init_radukca(                                          &
      ngrgas,grgas_addr)

! use switches for minor gases
USE lw_control_struct
USE ukca_d1_defs,  ONLY: ukca_sect, i_ukca_grg_o3, i_ukca_grg_ch4,&
                   i_ukca_grg_n2o, i_ukca_grg_f11, i_ukca_grg_f12,&
                   i_ukca_grg_f113, i_ukca_grg_h22
USE ukca_option_mod, ONLY: L_ukca_stratcfc, L_ukca_radf113,       &
                     L_ukca_radf12, L_ukca_radf22, L_ukca_radch4, &
                     L_ukca_rado3, L_ukca_radf11, L_ukca_radn2o,  &
                     l_ukca
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims, udims,             &
                                 vdims, wdims, tdims_s,           &
                                 o3dims2, wdims_s
USE atm_d1_indices_mod, ONLY: jtracer, jtr_ukca
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,  ONLY: mype

USE ukca_feedback_mod, ONLY: p_o3, p_ch4, p_n2o, p_f11, &
                             p_f12, p_f113, p_f22
USE submodel_mod, ONLY: submodel_for_sm, atmos_im
USE d1_array_mod, ONLY: d1_object_type, d1_section, d1_item,      &
                        d1_address, prognostic, d1, d1_addr,      &
                        no_obj_d1
USE nlsizes_namelist_mod, ONLY:                                   &
    len_tot, model_levels, n_cca_lev, n_obj_d1_max, ozone_levels, &
    sm_levels, st_levels, theta_off_size, tpps_ozone_levels,      &
    tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca, tr_vars

USE free_tracers_inputs_mod, ONLY: a_tracer_first
USE ukca_tracer_stash, ONLY: a_ukca_first

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
! initialization of radiative feedback
!  Initialise and check address array to feed chemical tracers from
!  UKCA into radiation scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

! Subroutine arguments
!   Scalar arguments with Intent(In):
INTEGER, INTENT(IN) :: ngrgas

!   Array  arguments with Intent(In):
!   Scalar arguments with Intent(InOut):
!   Array  arguments with Intent(InOut):
!   Scalar arguments with Intent(Out):
!   Array  arguments with Intent(Out):
INTEGER, INTENT(INOUT) :: grgas_addr(ngrgas)

! Local parameters:
CHARACTER (LEN=* ), PARAMETER :: RoutineName='INIT_RADUKCA'

! Local scalars:
INTEGER ::  i                               ! Loop counter
INTEGER :: m_atm_modl                       !
INTEGER :: section                          ! model section
INTEGER :: item                             ! item
INTEGER :: address                          ! address in D1
INTEGER :: tr_levs                          ! UKCA tracer levels
INTEGER            :: errcode               ! error code
CHARACTER (LEN=errormessagelength) :: cmessage              ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
tr_levs = tdims%k_len

m_atm_modl = submodel_for_sm(atmos_im)
DO i=1,no_obj_d1(m_atm_modl)
  section = d1_addr(d1_section, i, m_atm_modl)
  IF ((d1_addr(d1_object_type,i,m_atm_modl)==prognostic) .AND.  &
    (section == ukca_sect)) THEN
    item    = d1_addr(d1_item,    i, m_atm_modl)
    address = d1_addr(d1_address, i, m_atm_modl)
    SELECT CASE (item)
    CASE (i_ukca_grg_o3)
      IF (L_ukca_rado3  ) grgas_addr(p_o3  )= address
    CASE (i_ukca_grg_ch4)
      IF (L_ukca_radch4 ) grgas_addr(p_ch4 )= address
    CASE (i_ukca_grg_n2o)
      IF (L_ukca_radn2o ) grgas_addr(p_n2o )= address
    CASE (i_ukca_grg_f11)
      IF (L_ukca_radf11 ) grgas_addr(p_f11 )= address
    CASE (i_ukca_grg_f12)
      IF (L_ukca_radf12 ) grgas_addr(p_f12 )= address
    CASE (i_ukca_grg_f113)
      IF (L_ukca_radf113) grgas_addr(p_f113)= address
    CASE (i_ukca_grg_h22)
      IF (L_ukca_radf22 ) grgas_addr(p_f22 )= address
    END SELECT
  END IF
END DO

! Turn grgas_addr into number of tracer in free_tracers array
IF (ukca_sect == 33) THEN
  grgas_addr = (grgas_addr -                                  &
    jtracer(tdims_s%k_start,a_tracer_first)) /                &
    (theta_off_size * tr_levs) + 1
ELSE IF (ukca_sect == 34) THEN
  grgas_addr = (grgas_addr -                                  &
    jtr_ukca(tdims_s%k_start,a_ukca_first)) /             &
    (theta_off_size * tr_levs) + 1
ELSE
  cmessage=' Tracer section not identified'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  IF ((grgas_addr(p_o3) < 0) .AND. L_ukca_rado3) THEN
    cmessage='WARNING: O3 not found among chemical species.'
    errcode = -1
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_ch4) < 0) .AND. L_ukca_radch4) THEN
    cmessage='WARNING: CH4 not found among chemical species.'
    errcode = -2
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF (ozone_levels /= model_levels) THEN
    cmessage = 'Ozone levels must equal model levels.'
    errcode = 1
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF (ozone_levels /= tr_levels) THEN
    cmessage='Tracer levels must equal model levels.'
    errcode = 2
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_n2o) > 0) .AND.                               &
      (.NOT. lw_control(1)%l_n2o)) THEN
    cmessage='N2O found but absorption by N2O not selected.'
    errcode = -3
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  IF ((grgas_addr(p_f11) > 0) .AND.                               &
      (.NOT. lw_control(1)%l_cfc11)) THEN
    cmessage='CFCl1 found but absorption by CFCl1 not selected.'
    errcode = -4
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_f12) > 0) .AND.                               &
      (.NOT. lw_control(1)%l_cfc12)) THEN
    cmessage='CF2Cl2 found but absorption by CF2Cl2 not selected.'
    errcode = -5
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_f113) > 0) .AND.                              &
      (.NOT. lw_control(1)%l_cfc113)) THEN
    cmessage=                                                     &
    'CFC-113 found but absorption by CFC-113 not selected.'
    errcode = -6
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_f22) > 0) .AND.                               &
      (.NOT. lw_control(1)%l_hcfc22)) THEN
    cmessage=                                                     &
    'CHF2Cl found but absorption by CHF2Cl not selected.'
    errcode = -7
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_f11) > 0) .AND.                               &
      (lw_control(1)%l_cfc11) .AND.                               &
      .NOT. (L_ukca_stratcfc)) THEN
    cmessage=                                                     &
    'L_ukca_stratcfc not selected, CFC-11 not treated properly.'
    errcode = -8
    CALL ereport(RoutineName,errcode,cmessage)
  END IF

  IF ((grgas_addr(p_f12) > 0) .AND.                               &
      (lw_control(1)%l_cfc12) .AND.                               &
      .NOT. (L_ukca_stratcfc)) THEN
    cmessage=                                                     &
    'L_ukca_stratcfc not selected, CFC-12 not treated properly.'
    errcode = -9
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_radukca
END MODULE init_radukca_mod

