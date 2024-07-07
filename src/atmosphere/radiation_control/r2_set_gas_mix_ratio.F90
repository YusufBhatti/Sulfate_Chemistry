! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE r2_set_gas_mix_ratio_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'R2_SET_GAS_MIX_RATIO_MOD'
CONTAINS

SUBROUTINE r2_set_gas_mix_ratio(control, spectrum, atm,                        &
! Grid
  n_profile, n_layer, nozone, nd_field, i_gather,                              &
! Gas mass mixing ratios
  h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio,                   &
  c11_mix_ratio, c12_mix_ratio, o2_mix_ratio,                                  &
  c113_mix_ratio, c114_mix_ratio, hcfc22_mix_ratio,                            &
  hfc125_mix_ratio, hfc134a_mix_ratio,                                         &
! 3D CO2
  co2_dim1, co2_dim2, co2_3d, l_co2_3d,                                        &
! Chemical greenhouse gas fields
  ngrgas, grgas_field)

! Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.

USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_spectrum, ONLY: StrSpecData
USE def_atm,      ONLY: StrAtm
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE ereport_mod,  ONLY: ereport
USE gas_list_pcf, ONLY: ip_h2o, ip_co2, ip_o3, ip_n2o,                         &
                        ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2,           &
                        ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12,            &
                        ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a,           &
                        ip_cfc114, ip_tio, ip_vo, ip_h2, ip_he,                &
                        ip_na, ip_k, ip_li, ip_rb, ip_cs
USE ukca_feedback_mod, ONLY: p_ch4, p_n2o, p_f11, p_f12, p_f113, p_f22
USE mmr_BS1999_mod, ONLY: calc_gas_mixing_ratio_BS1999
USE rad_input_mod, ONLY: l_BS1999_abundances, l_extra_top
USE nlsizes_namelist_mod, ONLY: model_levels

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Spectral data:
TYPE (StrSpecData), INTENT(IN)    :: spectrum

! Atmospheric properties:
TYPE(StrAtm),       INTENT(INOUT) :: atm


! Sizes of arrays:
INTEGER, INTENT(IN) ::                                                         &
  n_profile,                                                                   &
!   Number of profiles
  n_layer,                                                                     &
!   Number of radiative layers
  nozone,                                                                      &
!   Number of ozone levels
  nd_field
!   Size of array from UM

! Gathering array:
INTEGER, INTENT(IN) ::                                                         &
  i_gather(nd_field)
!   List of points to be gathered

! 3D CO2
INTEGER, INTENT(IN) :: co2_dim1, co2_dim2
!   Dimensions of CO2_3D field
LOGICAL, INTENT(IN) :: l_co2_3d
!   Controls use of 3D co2 field

! Mixing ratios supplied:
REAL, INTENT(IN) ::                                                            &
  h2o(nd_field, model_levels),                                                 &
!   Mass mixing ratio of water vapour
  co2,                                                                         &
!   Mass mixing ratio of carbon dioxide
  co2_3d(co2_dim1, co2_dim2),                                                  &
!   3D mass mixing ratio of CO2 (full field)
  o3(nd_field, nozone),                                                        &
!   Mass mixing ratio of ozone
  n2o_mix_ratio,                                                               &
!   Mass mixing ratio of nitrous oxide
  ch4_mix_ratio,                                                               &
!   Mass mixing ratio of methane
  so2_mix_ratio,                                                               &
!   Mass mixing ratio of sulphur dioxide
  c11_mix_ratio,                                                               &
!   Mass mixing ratio of CFC11
  c12_mix_ratio,                                                               &
!   Mass mixing ratio of CFC12
  o2_mix_ratio,                                                                &
!   Mass mixing ratio of O2
  c113_mix_ratio,                                                              &
!   Mass mixing ratio of CFC113
  c114_mix_ratio,                                                              &
!   Mass mixing ratio of CFC114
  hcfc22_mix_ratio,                                                            &
!   Mass mixing ratio of HCFC22
  hfc125_mix_ratio,                                                            &
!   Mass mixing ratio of HFC125
  hfc134a_mix_ratio
!   Mass mixing ratio of HFC134a

! Chemical greenhouse gas fields
INTEGER, INTENT(IN) :: ngrgas
REAL, INTENT(IN) :: grgas_field(nd_field, model_levels, ngrgas)


! Local variables.

! Pointers to gases:
INTEGER ::                                                                     &
  iump_h2o,                                                                    &
!   Pointer to Water Vapour
  iump_co2,                                                                    &
!   Pointer to Carbon Dioxide
  iump_o3,                                                                     &
!   Pointer to Ozone
  iump_n2o,                                                                    &
!   Pointer to Nitous Oxide
  iump_co,                                                                     &
!   Pointer to Carbon Monoxide
  iump_ch4,                                                                    &
!   Pointer to Methane
  iump_so2,                                                                    &
!   Pointer to Sulphur Dioxide
  iump_cfc11,                                                                  &
!   Pointer to CFC11
  iump_cfc12,                                                                  &
!   Pointer to CFC12
  iump_o2,                                                                     &
!   Pointer to O2
  iump_nh3,                                                                    &
!   Pointer to Ammonia
  iump_cfc113,                                                                 &
!   Pointer to CFC113
  iump_cfc114,                                                                 &
!   Pointer to CFC114
  iump_hcfc22,                                                                 &
!   Pointer to HCFC22
  iump_hfc125,                                                                 &
!   Pointer to HFC125
  iump_hfc134a,                                                                &
!   Pointer to HFC134a
  iump_tio,                                                                    &
!   Pointer to Titanium Oxide
  iump_vo,                                                                     &
!   Pointer to Vanadium Oxide
  iump_h2,                                                                     &
!   Pointer to Hydrogen
  iump_he,                                                                     &
!   Pointer to Helium
  iump_na,                                                                     &
!   Pointer to Sodium
  iump_k,                                                                      &
!   Pointer to Potassium
  iump_li,                                                                     &
!   Pointer to Lithium
  iump_rb,                                                                     &
!   Pointer to Rubidium
  iump_cs
!   Pointer to Cesium

INTEGER ::                                                                     &
  i, l,                                                                        &
!   Loop variables
  lg,                                                                          &
!   Corresponding ungathered index
  i_top_copy
!   Topmost layer where properties are set by copying the input fields.

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'R2_SET_GAS_MIX_RATIO'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Match the indexing numbers of gaseous species in the spectral
! file with actual types of gases known to the UM.

! Set all pointers to 0 initially to flag missing gases.
iump_h2o=0
iump_co2=0
iump_o3=0
iump_n2o=0
iump_co=0
iump_ch4=0
iump_so2=0
iump_cfc11=0
iump_cfc12=0
iump_o2=0
iump_nh3=0
iump_cfc113=0
iump_cfc114=0
iump_hcfc22=0
iump_hfc125=0
iump_hfc134a=0
iump_tio=0
iump_vo=0
iump_h2=0
iump_he=0
iump_na=0
iump_k=0
iump_li=0
iump_rb=0
iump_cs=0

DO i=1, spectrum%gas%n_absorb
  SELECT CASE(spectrum%gas%type_absorb(i))
  CASE(ip_h2o)
     iump_h2o=i
  CASE(ip_co2)
     iump_co2=i
  CASE(ip_o3)
     iump_o3=i
  CASE(ip_n2o)
     iump_n2o=i
  CASE(ip_co)
     iump_co=i
  CASE(ip_ch4)
     iump_ch4=i
  CASE(ip_so2)
     iump_so2=i
  CASE(ip_cfc11)
     iump_cfc11=i
  CASE(ip_cfc12)
     iump_cfc12=i
  CASE(ip_o2)
     iump_o2=i
  CASE(ip_nh3)
     iump_nh3=i
  CASE(ip_cfc113)
     iump_cfc113=i
  CASE(ip_cfc114)
     iump_cfc114=i
  CASE(ip_hcfc22)
     iump_hcfc22=i
  CASE(ip_hfc125)
     iump_hfc125=i
  CASE(ip_hfc134a)
     iump_hfc134a=i
  CASE(ip_tio)
     iump_tio=i
  CASE(ip_vo)
     iump_vo=i
  CASE(ip_h2)
     iump_h2=i
  CASE(ip_he)
     iump_he=i
  CASE(ip_na)
     iump_na=i
  CASE(ip_k)
     iump_k=i
  CASE(ip_li)
     iump_li=i
  CASE(ip_rb)
     iump_rb=i
  CASE(ip_cs)
     iump_cs=i
  END SELECT
END DO


IF (l_extra_top) THEN
  ! The second radiative layer will be the first to have properties
  ! set by copying input fields.
  i_top_copy=2
ELSE
  ! The first radiative layer will be the first to have properties
  ! set by copying input fields.
  i_top_copy=1
END IF


! Water Vapour
IF (control%l_h2o .AND. iump_h2o > 0) THEN
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_h2o,                    &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_h2o))
  ELSE
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_h2o)=MAX(h2o(lg, model_levels), 0.0e+00)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_h2o)=MAX(h2o(lg, n_layer-i+1), 0.0e+00)
      END DO
    END DO
  END IF
END IF


! Carbon Dioxide
IF (control%l_co2 .AND. iump_co2 > 0) THEN
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_co2,                    &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_co2))
  ELSE IF (l_co2_3d) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_co2)=co2_3d(lg, model_levels)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_co2)=co2_3d(lg, n_layer-i+1)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_co2)=co2
      END DO
    END DO
  END IF
END IF


! Ozone
IF (control%l_o3 .AND. iump_o3 > 0) THEN
  ! The input field of ozone is supplied on NOZONE levels.
  ! These values apply to the upper layers used by the full UM.
  ! If NOZONE is smaller than NLEVS, the mixing ratio on the
  ! bottom level supplied is copied to lower levels. If an
  ! extra top level is used its mixing ratio is set by copying
  ! the value for the top non-radiative level.
  IF (l_extra_top) THEN
    DO l=1, n_profile
      lg=i_gather(l)
      atm%gas_mix_ratio(l, 1, iump_o3)=o3(lg, nozone)
    END DO
  END IF
  DO i=i_top_copy, nozone+i_top_copy-1
    DO l=1, n_profile
      lg=i_gather(l)
      atm%gas_mix_ratio(l, i, iump_o3)=o3(lg, nozone+i_top_copy-i)
    END DO
  END DO
  DO i=nozone+i_top_copy, n_layer
    DO l=1, n_profile
      lg=i_gather(l)
      atm%gas_mix_ratio(l, i, iump_o3)=o3(lg, 1)
    END DO
  END DO
  atm%gas_mix_ratio(1:n_profile, 1:n_layer, iump_o3) =                         &
     MAX(atm%gas_mix_ratio(1:n_profile, 1:n_layer, iump_o3), 0.0)
END IF


! Nitrous Oxide
IF (control%l_n2o .AND. iump_n2o > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (ngrgas >= p_n2o) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_n2o)=grgas_field(lg, model_levels,p_n2o)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_n2o)=grgas_field(lg, n_layer-i+1,p_n2o)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_n2o)=n2o_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_n2o) THEN
  cmessage = 'N2O is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Methane
IF (control%l_ch4 .AND. iump_ch4 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_ch4,                    &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_ch4))
  ELSE IF (ngrgas >= p_ch4) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_ch4)=grgas_field(lg, model_levels,p_ch4)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_ch4)=grgas_field(lg, n_layer-i+1,p_ch4)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_ch4)=ch4_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_ch4) THEN
  cmessage = 'CH4 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! CFC-11
IF (control%l_cfc11 .AND. iump_cfc11 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (ngrgas >= p_f11) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc11)=grgas_field(lg, model_levels,p_f11)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc11)=grgas_field(lg, n_layer-i+1,p_f11)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_cfc11)=c11_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_cfc11) THEN
  cmessage = 'CFC11 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! CFC-12
IF (control%l_cfc12 .AND. iump_cfc12 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (ngrgas >= p_f12) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc12)=grgas_field(lg, model_levels,p_f12)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc12)=grgas_field(lg, n_layer-i+1,p_f12)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_cfc12)=c12_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_cfc12) THEN
  cmessage = 'CFC12 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Oxygen
IF (control%l_o2 .AND. iump_o2 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  DO i=1, n_layer
    DO l=1, n_profile
      atm%gas_mix_ratio(l, i, iump_o2)=o2_mix_ratio
    END DO
  END DO
ELSE IF (control%l_o2) THEN
  cmessage = 'O2 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Sulphur Dioxide
IF (control%l_so2 .AND. iump_so2 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  DO i=1, n_layer
    DO l=1, n_profile
      atm%gas_mix_ratio(l, i, iump_so2)=so2_mix_ratio
    END DO
  END DO
ELSE IF (control%l_so2) THEN
  cmessage = 'SO2 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! CFC-113
IF (control%l_cfc113 .AND. iump_cfc113 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (ngrgas >= p_f113) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc113)=grgas_field(lg,model_levels,p_f113)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_cfc113)=grgas_field(lg,n_layer-i+1,p_f113)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_cfc113)=c113_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_cfc113) THEN
  cmessage = 'CFC113 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! CFC-114
IF (control%l_cfc114 .AND. iump_cfc114 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  DO i=1, n_layer
    DO l=1, n_profile
      atm%gas_mix_ratio(l, i, iump_cfc114)=c114_mix_ratio
    END DO
  END DO
ELSE IF (control%l_cfc114) THEN
  cmessage = 'CFC114 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! HCFC-22
IF (control%l_hcfc22 .AND. iump_hcfc22 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  IF (ngrgas >= p_f22) THEN
    DO i=1, n_layer-model_levels
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_hcfc22)=grgas_field(lg,model_levels,p_f22)
      END DO
    END DO
    DO i=n_layer-model_levels+1, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        atm%gas_mix_ratio(l, i, iump_hcfc22)=grgas_field(lg,n_layer-i+1,p_f22)
      END DO
    END DO
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_hcfc22)=hcfc22_mix_ratio
      END DO
    END DO
  END IF
ELSE IF (control%l_hcfc22) THEN
  cmessage = 'HCFC22 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! HFC-125
IF (control%l_hfc125 .AND. iump_hfc125 > 0) THEN
  ! The gas is in the spectral file and has been selected.
  DO i=1, n_layer
    DO l=1, n_profile
      atm%gas_mix_ratio(l, i, iump_hfc125)=hfc125_mix_ratio
    END DO
  END DO
ELSE IF (control%l_hfc125) THEN
  cmessage = 'HFC125 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! HFC-134a
IF (control%l_hfc134a .AND. iump_hfc134a > 0) THEN
  ! The gas is in the spectral file and has been selected.
  DO i=1, n_layer
    DO l=1, n_profile
      atm%gas_mix_ratio(l, i, iump_hfc134a)=hfc134a_mix_ratio
    END DO
  END DO
ELSE IF (control%l_hfc134a) THEN
  cmessage = 'HFC134a is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Carbon Monoxide
IF (control%l_co .AND. iump_co > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_co,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_co))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_co)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_co) THEN
  cmessage = 'CO is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Ammonia
IF (control%l_nh3 .AND. iump_nh3 > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_nh3,                    &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_nh3))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_nh3)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_nh3) THEN
  cmessage = 'NH3 is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Titanium Oxide
IF (control%l_tio .AND. iump_tio > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_tio,                    &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_tio))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_tio)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_tio) THEN
  cmessage = 'TiO is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Vanadium Oxide
IF (control%l_vo .AND. iump_vo > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_vo,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_vo))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_vo)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_vo) THEN
  cmessage = 'VO is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Hydrogen
IF (control%l_h2 .AND. iump_h2 > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_h2,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_h2))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_h2)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_h2) THEN
  cmessage = 'H2-H2 CIA is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Helium
IF (control%l_he .AND. iump_he > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_he,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_he))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_he)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_he) THEN
  cmessage = 'H2-He CIA is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Sodium
IF (control%l_na .AND. iump_na > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_na,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_na))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_na)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_na) THEN
  cmessage = 'Na is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Potassium
IF (control%l_k .AND. iump_k > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_k,                      &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_k))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_k)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_k) THEN
  cmessage = 'K is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Lithium
IF (control%l_li .AND. iump_li > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_li,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_li))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_li)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_li) THEN
  cmessage = 'Li is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Rubidium
IF (control%l_rb .AND. iump_rb > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_rb,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_rb))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_rb)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_rb) THEN
  cmessage = 'Rb is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


! Cesium
IF (control%l_cs .AND. iump_cs > 0) THEN
  ! The gas is in the spectral file. Currently the only option is to set
  ! abundances suitable for hot jupiter / gas giant atmospheres.
  IF (l_BS1999_abundances) THEN
    CALL calc_gas_mixing_ratio_BS1999(atm%p, atm%t, ip_cs,                     &
      n_profile, n_layer, atm%gas_mix_ratio(:, :, iump_cs))
  ELSE
    DO i=1, n_layer
      DO l=1, n_profile
        atm%gas_mix_ratio(l, i, iump_cs)=0.0e+00
      END DO
    END DO
  END IF
ELSE IF (control%l_cs) THEN
  cmessage = 'Cs is not in the spectral file but has been selected'
  ierr=i_err_fatal
  GO TO 9999
END IF


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE r2_set_gas_mix_ratio
END MODULE r2_set_gas_mix_ratio_mod
