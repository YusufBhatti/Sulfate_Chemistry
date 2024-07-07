! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to hold species mass mixing ratios from atmos_physics1 (etc)
!  which will be used to prescribe the surface concentrations
!  in UKCA when the logical L_ukca_prescribech4 is set to true, and when
!  lower boundary conditions are required.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_trace_gas_mixratio

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE

SAVE
PUBLIC

REAL :: um_ch4_for_ukca
REAL :: um_co2_for_ukca
REAL :: um_n2o_for_ukca
REAL :: um_o2_for_ukca
REAL :: um_cfc11_for_ukca
REAL :: um_cfc12_for_ukca
REAL :: um_cfc113_for_ukca
REAL :: um_hcfc22_for_ukca
REAL :: um_hfc125_for_ukca
REAL :: um_hfc134a_for_ukca
REAL :: um_mebr_for_ukca
REAL :: um_mecl_for_ukca
REAL :: um_ch2br2_for_ukca
REAL :: um_h2_for_ukca
REAL :: um_n2_for_ukca
REAL :: um_cfc114_for_ukca
REAL :: um_cfc115_for_ukca
REAL :: um_ccl4_for_ukca
REAL :: um_meccl3_for_ukca
REAL :: um_hcfc141b_for_ukca
REAL :: um_hcfc142b_for_ukca
REAL :: um_h1211_for_ukca
REAL :: um_h1202_for_ukca
REAL :: um_h1301_for_ukca
REAL :: um_h2402_for_ukca
REAL :: um_cos_for_ukca

! These to be initialised from GUI values
REAL :: MeBrMMR
REAL :: MeClMMR
REAL :: CH2Br2MMR
REAL :: h2mmr
REAL :: n2mmr
REAL :: cfc114mmr
REAL :: cfc115mmr
REAL :: CCl4MMR
REAL :: MeCCl3MMR
REAL :: HCFC141bMMR
REAL :: HCFC142bMMR
REAL :: h1211mmr
REAL :: h1202mmr
REAL :: h1301mmr
REAL :: h2402mmr
REAL :: cosmmr

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_TRACE_GAS_MIXRATIO'

CONTAINS

! ######################################################################
! Description:
!  Subroutine to copy species mass mixing ratios from atmos_physics1
!  and UKCA GUI settings to module for use later in UKCA as lower
!  boundary conditions or as species concentrations (CO2, O2, H2 or N2)
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Called from Atmos_Physics1
!
! Method:
!
! ----------------------------------------------------------------------
!
SUBROUTINE ukca_set_trace_gas_mixratio(                           &
     ch4_mix_ratio, co2_mmr, n2o_mix_ratio,                       &
     o2mmr, cfc11_mix_ratio, cfc12_mix_ratio,                     &
     c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr)

USE ukca_option_mod,      ONLY: ukca_MeBrMMR, ukca_MeClMMR,       &
                                ukca_CH2Br2MMR, ukca_H2MMR,       &
                                ukca_N2MMR, ukca_CFC115MMR,       &
                                ukca_CCl4MMR, ukca_MeCCl3MMR,     &
                                ukca_HCFC141bMMR, ukca_H1211MMR,  &
                                ukca_HCFC142bMMR, ukca_H1202MMR,  &
                                ukca_H1301MMR, ukca_H2402MMR,     &
                                ukca_COSMMR, L_ukca_strat,        &
                                L_ukca_strattrop, L_ukca_stratcfc,&
                                L_ukca_achem
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim
USE ereport_mod,          ONLY: ereport
USE missing_data_mod,      ONLY: rmdi, imdi

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!  Mixing ratios from GUI LW radiation settings  - these may be time varying
REAL, INTENT(IN) :: ch4_mix_ratio
REAL, INTENT(IN) :: co2_mmr
REAL, INTENT(IN) :: n2o_mix_ratio
REAL, INTENT(IN) :: o2mmr
REAL, INTENT(IN) :: cfc11_mix_ratio
REAL, INTENT(IN) :: cfc12_mix_ratio
REAL, INTENT(IN) :: c113mmr
REAL, INTENT(IN) :: c114mmr
REAL, INTENT(IN) :: hcfc22mmr
REAL, INTENT(IN) :: hfc125mmr
REAL, INTENT(IN) :: hfc134ammr

LOGICAL, SAVE :: L_first=.TRUE.
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: ierr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_TRACE_GAS_MIXRATIO'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Pick up GUI values for UKCA gases,
!  currently only need to do on first timestep
!  any missing data will be flagged below
IF (L_first) THEN
  MeBrMMR     = ukca_MeBrMMR
  MeClMMR     = ukca_MeClMMR
  CH2Br2MMR   = ukca_CH2Br2MMR
  h2mmr       = ukca_H2MMR
  n2mmr       = ukca_N2MMR
  cfc115mmr   = ukca_CFC115MMR
  CCl4MMR     = ukca_CCl4MMR
  MeCCl3MMR   = ukca_MeCCl3MMR
  HCFC141bMMR = ukca_HCFC141bMMR
  HCFC142bMMR = ukca_HCFC142bMMR
  h1211mmr    = ukca_H1211MMR
  h1202mmr    = ukca_H1202MMR
  h1301mmr    = ukca_H1301MMR
  h2402mmr    = ukca_H2402MMR
  cosmmr      = ukca_COSMMR
END IF

IF (ch4_mix_ratio < 0.0) THEN
  IF (L_first) THEN
    cmessage='Missing value for trace_gas_mix_ratio CH4'
    WRITE(umMessage,'(2A)') 'Set value for CH4 in GUI panel -> ',       &
        'setting to default pre-industrial level'
    CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
    ierr = -1
    CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
  END IF
  um_ch4_for_ukca = 4.360e-07
ELSE
  um_ch4_for_ukca = ch4_mix_ratio
END IF

IF (co2_mmr < 0.0) THEN
  IF (L_first) THEN
    cmessage='Missing value for trace_gas_mix_ratio CO2'
    WRITE(umMessage,'(2A)') 'Set value for CO2 in GUI panel -> ',       &
               'setting to default pre-industrial level'
    CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
    ierr = -2
    CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
  END IF
  um_co2_for_ukca = 4.31980e-04
ELSE
  um_co2_for_ukca = co2_mmr
END IF

IF (h2mmr < 0.0) THEN
  IF (L_first) THEN
    cmessage='Missing value for trace_gas_mix_ratio H2'
    WRITE(umMessage,'(2A)') 'Set value for H2 in GUI panel -> ',        &
               'setting to default pre-industrial level'
    CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
    ierr = -14
    CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
  END IF
  um_h2_for_ukca = 3.452e-8 ! as in cinit
ELSE
  um_h2_for_ukca = h2mmr
END IF

! These species are only set for stratospheric chemistry schemes
IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc) THEN

  IF (n2o_mix_ratio < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio N2O'
      WRITE(umMessage,'(2A)') 'Set value for N2O in GUI panel -> ',     &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -3
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_n2o_for_ukca = 4.180e-07
  ELSE
    um_n2o_for_ukca = n2o_mix_ratio
  END IF

  IF (o2mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio O2'
      WRITE(umMessage,'(2A)') 'Set value for O2 in GUI panel -> ',      &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -4
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_o2_for_ukca = 0.23144 ! from cinit
  ELSE
    um_o2_for_ukca = o2mmr
  END IF


  IF (cfc11_mix_ratio  < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CFC11'
      WRITE(umMessage,'(2A)') 'Set value for CFC11 in GUI panel -> ',   &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -5
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cfc11_for_ukca = 0.0
  ELSE
    um_cfc11_for_ukca = cfc11_mix_ratio
  END IF

  IF (cfc12_mix_ratio < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CFC12'
      WRITE(umMessage,'(2A)') 'Set value for CFC12 in GUI panel -> ',   &
          'setting to pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -6
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cfc12_for_ukca = 0.0
  ELSE
    um_cfc12_for_ukca = cfc12_mix_ratio
  END IF

  IF (c113mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CFC113'
      WRITE(umMessage,'(2A)') 'Set value for CFC113 in GUI panel -> ',  &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -7
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cfc113_for_ukca = 0.0
  ELSE
    um_cfc113_for_ukca = c113mmr
  END IF

  IF (hcfc22mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio HCFC22'
      WRITE(umMessage,'(2A)') 'Set value for HCFC22 in GUI panel -> ',  &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -8
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_hcfc22_for_ukca = 0.0
  ELSE
    um_hcfc22_for_ukca = hcfc22mmr
  END IF

  IF (hfc125mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio HFC125'
      WRITE(umMessage,'(2A)') 'Set value for HFC125 in GUI panel -> ',  &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -9
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_hfc125_for_ukca = 0.0
  ELSE
    um_hfc125_for_ukca = hfc125mmr
  END IF

  IF (hfc134ammr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio HFC134A'
      WRITE(umMessage,'(2A)') 'Set value for HFC134A in GUI panel -> ', &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -10
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_hfc134a_for_ukca = 0.0
  ELSE
    um_hfc134a_for_ukca = hfc134ammr
  END IF

  IF (MeBrMMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio MeBr'
      WRITE(umMessage,'(2A)') 'Set value for MeBr in GUI -> ',          &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -11
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_mebr_for_ukca = 1.64025e-11 ! 5pptv
  ELSE
    um_mebr_for_ukca = MeBrMMR
  END IF

  IF (MeClMMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio MeCl'
      WRITE(umMessage,'(2A)') 'Set value for MeCl in GUI -> ',          &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -12
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_mecl_for_ukca = 8.36736e-10 ! 480pptv
  ELSE
    um_mecl_for_ukca = MeClMMR
  END IF

  IF (CH2Br2MMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CH2Br2'
      WRITE(umMessage,'(2A)') 'Set value for CH2Br2 in GUI -> ',        &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -13
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_ch2br2_for_ukca = 18.0186e-12 ! 3pptv
  ELSE
    um_ch2br2_for_ukca = CH2Br2MMR
  END IF

  IF (n2mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio N2'
      WRITE(umMessage,'(2A)') 'Set value for N2 in GUI -> ',            &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -15
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_n2_for_ukca = 0.754682 ! as in cinit
  ELSE
    um_n2_for_ukca = n2mmr
  END IF

  IF (c114mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CFC114'
      WRITE(umMessage,'(2A)') 'Set value for CFC-114 in GUI -> ',       &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -16
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cfc114_for_ukca = 0.0
  ELSE
    um_cfc114_for_ukca = c114mmr
  END IF

  IF (cfc115mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CFC115'
      WRITE(umMessage,'(2A)') 'Set value for CFC-115 in GUI -> ',       &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -17
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cfc115_for_ukca = 0.0
  ELSE
    um_cfc115_for_ukca = cfc115mmr
  END IF

  IF (CCl4MMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio CCl4'
      WRITE(umMessage,'(2A)') 'Set value for CCl4 in GUI -> ',          &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -18
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_ccl4_for_ukca = 0.0
  ELSE
    um_ccl4_for_ukca = CCl4MMR
  END IF

  IF (MeCCl3MMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio MeCCl3'
      WRITE(umMessage,'(2A)') 'Set value for MeCCl3 in GUI -> ',        &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -18
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_meccl3_for_ukca = 0.0
  ELSE
    um_meccl3_for_ukca = MeCCl3MMR
  END IF

  IF (HCFC141bMMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio HCFC141b'
      WRITE(umMessage,'(2A)') 'Set value for HCFC141b in GUI -> ',      &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -19
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_hcfc141b_for_ukca = 0.0
  ELSE
    um_hcfc141b_for_ukca = HCFC141bMMR
  END IF

  IF (HCFC142bMMR < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio HCFC142b'
      WRITE(umMessage,'(2A)') 'Set value for HCFC142b in GUI -> ',      &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -20
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_hcfc142b_for_ukca = 0.0
  ELSE
    um_hcfc142b_for_ukca = HCFC142bMMR
  END IF

  IF (h1211mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio H1211'
      WRITE(umMessage,'(2A)')'Set value for H1211 in GUI hand-edit -> ',&
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -21
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_h1211_for_ukca = 0.0
  ELSE
    um_h1211_for_ukca = h1211mmr
  END IF

  IF (h1202mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio H1202'
      WRITE(umMessage,'(2A)') 'Set value for H1202 in GUI -> ',         &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -22
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_h1202_for_ukca = 0.0
  ELSE
    um_h1202_for_ukca = h1202mmr
  END IF

  IF (h1301mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio H1301'
      WRITE(umMessage,'(2A)') 'Set value for H1301 in GUI -> ',         &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -23
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_h1301_for_ukca = 0.0
  ELSE
    um_h1301_for_ukca = h1301mmr
  END IF

  IF (h2402mmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio H2402'
      WRITE(umMessage,'(2A)') 'Set value for H2402 in GUI -> ',         &
          'setting to default pre-industrial level'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -24
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_h2402_for_ukca = 0.0
  ELSE
    um_h2402_for_ukca = h2402mmr
  END IF

  IF (cosmmr < 0.0) THEN
    IF (L_first) THEN
      cmessage='Missing value for trace_gas_mix_ratio COS  '
      WRITE(umMessage,'(A,A)') 'Set value for COS in GUI  -> ',         &
          'setting to default value of 520 pptm'
      CALL umPrint(umMessage,src='ukca_trace_gas_mixratio')
      ierr = -25
      CALL ereport('ukca_set_trace_gas_mixratio',ierr,cmessage)
    END IF
    um_cos_for_ukca = 520.0e-12
     ! see Montzka et al, JGR, 109, D22302, (2004)
  ELSE
    um_cos_for_ukca = cosmmr
  END IF

END IF    ! L_ukca_strat etc.


L_first = .FALSE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_trace_gas_mixratio

END MODULE ukca_trace_gas_mixratio
