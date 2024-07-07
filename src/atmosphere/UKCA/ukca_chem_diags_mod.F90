! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to contain code to place chemistry time step diagnostics 
!   into the section 50 STASHwork array
!
! Method:
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

MODULE ukca_chem_diags_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEM_DIAGS_MOD'

CONTAINS

SUBROUTINE ukca_chem_diags ( &
! IN model dimensions
  row_length, rows, model_levels,                  &
! IN fields for diagnostics
  nat_psc, trop_ch4_mol, trop_o3_mol, trop_oh_mol, &
  strat_ch4_mol, strat_ch4loss,                    &
  atm_ch4_mol, atm_co_mol, atm_n2o_mol,            &
  atm_cf2cl2_mol, atm_cfcl3_mol, atm_mebr_mol,     &
  atm_h2_mol, so4_sa,                              &
  dj,                                              &
! INOUT stash workspace
  len_stashwork, stashwork)

! Description:
!   To place chemistry time step diagnostics into the section 50 
!   STASHwork array

USE ukca_d1_defs,       ONLY: ukca_diag_sect
USE ukca_option_mod,    ONLY: jppj
USE ukca_chem_defs_mod, ONLY: ratj_defs
USE ukca_all_tracers_copy_mod, ONLY: ukca_q_increment
USE atm_fields_bounds_mod,     ONLY: tdims_s

USE submodel_mod,      ONLY: atmos_im
USE stash_array_mod,   ONLY: len_stlist, stindex, stlist,           &
                             num_stash_levels, stash_levels, si, sf
USE um_stashcode_mod,  ONLY: stashcode_ukca_chem_diag,              &
                             stashcode_ukca_nat,                    &
                             stashcode_ukca_trop_ch4,               &
                             stashcode_ukca_trop_o3,                &
                             stashcode_ukca_trop_oh,                &
                             stashcode_ukca_strat_ch4,              &
                             stashcode_ukca_strt_ch4_lss,           &
                             stashcode_ukca_jo1d,                   &
                             stashcode_ukca_jn2o,                   &
                             stashcode_ukca_atmos_ch4,              &
                             stashcode_ukca_atmos_co,               &
                             stashcode_ukca_atmos_n2o,              &
                             stashcode_ukca_atmos_cfc12,            &
                             stashcode_ukca_atmos_cfc11,            &
                             stashcode_ukca_atmos_ch3br,            &
                             stashcode_ukca_atmos_h2,               &
                             stashcode_ukca_h2o_incr,               &
                             stashcode_ukca_jo2,                    &
                             stashcode_ukca_jo3p,                   &
                             stashcode_ukca_so4_sad
USE um_parvars,   ONLY: at_extremity
USE ereport_mod,  ONLY: ereport
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length        ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! Diagnostic tracers
! Nitric acid trihydrate (kg(nat)/kg(air))
REAL, INTENT(IN) :: nat_psc(row_length,rows,model_levels)
! Trop CH4 burden (moles)  
REAL, INTENT(IN) :: trop_ch4_mol(row_length,rows,model_levels)
! Trop O3 burden (moles)
REAL, INTENT(IN) :: trop_o3_mol(row_length,rows,model_levels)
! Trop OH burden (moles)
REAL, INTENT(IN) :: trop_oh_mol(row_length,rows,model_levels)
! Strat CH4 burden (moles)
REAL, INTENT(IN) :: strat_ch4_mol(row_length,rows,model_levels)
! Strat CH4 loss (Moles/s)
REAL, INTENT(IN) :: strat_ch4loss(row_length,rows,model_levels)
! Photolysis rates (/s)
REAL, INTENT(IN) :: dj (row_length, rows ,model_levels, jppj)
! Atmospheric Burden of CH4 in moles
REAL, INTENT(IN) :: atm_ch4_mol(row_length,rows,model_levels)
! Atmospheric Burden of CO in moles
REAL, INTENT(IN) :: atm_co_mol(row_length,rows,model_levels)
! Atmospheric Burden of Nitrous Oxide (N2O) in moles
REAL, INTENT(IN) :: atm_n2o_mol(row_length,rows,model_levels)
! Atmospheric Burden of CFC-12 in moles
REAL, INTENT(IN) :: atm_cf2cl2_mol(row_length,rows,model_levels)
! Atmospheric Burden of CFC-11 in moles
REAL, INTENT(IN) :: atm_cfcl3_mol(row_length,rows,model_levels)
! Atmospheric Burden of CH3Br in moles
REAL, INTENT(IN) :: atm_mebr_mol(row_length,rows,model_levels)
! Atmospheric Burden of H2 in moles 
REAL, INTENT(IN) :: atm_h2_mol(row_length,rows,model_levels)
! Aerosol surface area used in chemistry
REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels)

! Diagnostics info
INTEGER, INTENT(IN) :: len_stashwork ! Length of diagnostics array
REAL, INTENT(INOUT) :: stashwork(len_stashwork)  ! STASH workspace

! Local variables
INTEGER :: item                                ! STASH item
INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0  ! DrHook tracing entry
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 ! DrHook tracing exit
INTEGER :: icode                               ! error code for EReport
INTEGER :: im_index                            ! internal model index
INTEGER :: ii                                  ! loop variable

INTEGER, SAVE :: i_jo1d = -1   ! index for photolysis rate of O3 to O1D (JO1D)
INTEGER, SAVE :: i_jno2 = -1   ! index for photolysis rate of NO2 (JNO2)
INTEGER, SAVE :: i_jo2  = -1   ! index for photolysis rate of O2 (O2 to O1D+O3P)
INTEGER, SAVE :: i_jo3p = -1   ! index for photolysis rate of O3 to O3P (JO3P)

REAL(KIND=jprb) :: zhook_handle ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEM_DIAGS' ! used for EReport
CHARACTER(LEN=errormessagelength) :: cmessage ! used for EReport

LOGICAL, SAVE :: first = .TRUE.

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise indices of some photolysis rates (e.g. JO1D, JNO2)
! and then find their right indices in the given chemistry scheme.

IF (first) THEN
  i_jo1d = -1
  i_jno2 = -1
  i_jo2  = -1
  i_jo3p = -1

  DO ii = 1, jppj
    ! O3 + hv --> O2 + O(1D) 
    IF (ratj_defs(ii)%fname == 'jo3a      ') THEN
      i_jo1d = ii
    END IF

    ! NO2 + hv --> NO + O(3P)  
    IF (ratj_defs(ii)%fname == 'jno2      ') THEN
      i_jno2 = ii
    END IF

    ! O2 + hv --> O1D + O3P 
    IF (ratj_defs(ii)%fname == 'jo2b      ') THEN
      i_jo2  = ii
    END IF

    ! O3 + hv --> O3P + O2
    IF (ratj_defs(ii)%fname == 'jo3b      ') THEN
      i_jo3p = ii
    END IF

  END DO

  first = .FALSE.
END IF

icode = 0 ! Initialise error status
im_index = 1

! ----------------------------------------------------------------------
!   Copy diagnostic information to STASHwork for STASH processing
! ----------------------------------------------------------------------
! DIAG.50218 Nitric acid trihydrate 
! ----------------------------------------------------------------------
item = stashcode_ukca_nat - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       nat_psc(:,:,:),                                                  &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.218 "
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50220 Trop CH4 burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       trop_ch4_mol (:,:,:),                                            &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.220"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50221 Trop O3 burden in mol 
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_o3 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       trop_o3_mol(:,:,:),                                              &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.221"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50222 Trop OH burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_trop_oh - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       trop_oh_mol(:,:,:),                                              &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.222"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50223 Strat CH4 burden in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_strat_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       strat_ch4_mol(:,:,:),                                            &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.223"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50226 Strat CH4 loss 
! ----------------------------------------------------------------------
item = stashcode_ukca_strt_ch4_lss - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       strat_ch4loss(:,:,:),                                            &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.226"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50228 Photolysis rate JO1D in s-1 
! ----------------------------------------------------------------------
item = stashcode_ukca_jo1d - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

  IF (i_jo1d < 0) THEN
    icode    = 1
    cmessage = "JO1D (STASH item 50.228) unavailable: index not found"
    CALL ereport(RoutineName,icode,cmessage)
  END IF

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       dj(:,:,:,i_jo1d),                                                &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist (1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,   &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Error in copydiag_3d for STASH item 50.228"
    CALL ereport(RoutineName, icode, cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50229 Photolysis rate JNO2 in s-1
! ----------------------------------------------------------------------
item = stashcode_ukca_jn2o - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

  IF (i_jno2 < 0) THEN
    icode    = 2
    cmessage = "JNO2 (STASH item 50.229) unavailable: index not found"
    CALL ereport(RoutineName,icode,cmessage)
  END IF

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       dj(:,:,:,i_jno2),                                                &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist (1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,   &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Error in copydiag_3d for STASH item 50.229"
    CALL ereport (RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50231 Atmospheric Burden of CH4 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_ch4 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_ch4_mol(:,:,:),                                              &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.231"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50232 Atmospheric Burden of CO in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_co - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_co_mol(:,:,:),                                               &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.232"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50233 Atmospheric Burden of N2O in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_n2o - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_n2o_mol(:,:,:),                                              &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.233"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50234 Atmospheric Burden of CFC-12 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_cfc12 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_cf2cl2_mol(:,:,:),                                           &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.233"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50235 Atmospheric Burden of CFC-11 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_cfc11 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_cfcl3_mol(:,:,:),                                            &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.235"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50236 Atmospheric Burden of CH3Br in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_ch3br - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_mebr_mol(:,:,:),                                             &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.236"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50237 Atmospheric Burden of H2 in mol
! ----------------------------------------------------------------------
item = stashcode_ukca_atmos_h2 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       atm_h2_mol(:,:,:),                                               &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.237"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50240 Net specific humidity change from chemistry
! ----------------------------------------------------------------------
item = stashcode_ukca_h2o_incr - 1000*stashcode_ukca_chem_diag

IF (sf(item,stashcode_ukca_chem_diag)) THEN
! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       ukca_q_increment(tdims_s%i_start:tdims_s%i_end,                  &
                        tdims_s%j_start:tdims_s%j_end,                  &
                        1:tdims_s%k_end),                               &
       row_length,rows,model_levels,0,0,                                &
       tdims_s%halo_i,tdims_s%halo_j,at_extremity,                      &
       stlist (1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,   &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,stashcode_ukca_chem_diag,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Error in copydiag_3d for STASH item 50.240"
    CALL ereport (RoutineName,icode,cmessage)
  END IF
END IF   ! sf(item,stashcode_ukca_chem_diag)


! ----------------------------------------------------------------------
! DIAG.50245 Photolysis rate JO2 in s-1 i.e. O2 -> O1D + O3P (jo2b)
! ----------------------------------------------------------------------
item = stashcode_ukca_jo2 - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

  IF (i_jo2 < 0) THEN
    icode    = 2
    cmessage = "JO2 (STASH item 50.245) unavailable: index not found"
    CALL ereport(RoutineName,icode,cmessage)
  END IF

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       dj(:,:,:,i_jo2),                                                 &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist (1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,   &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Error in copydiag_3d for STASH item 50.245"
    CALL ereport (RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50246 Photolysis rate JO3P in s-1 i.e. O3 -> O3P+O2 (jo3b)
! ----------------------------------------------------------------------
item = stashcode_ukca_jo3p - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

  IF (i_jo3p < 0) THEN
    icode    = 2
    cmessage = "JO3P (STASH item 50.246) unavailable: index not found"
    CALL ereport(RoutineName,icode,cmessage)
  END IF

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       dj(:,:,:,i_jo3p),                                                &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist (1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,   &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="Error in copydiag_3d for STASH item 50.246"
    CALL ereport (RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! --------------------------------------------------------------------
! DIAG.50256 Aerosol surface area density as used in heteogenous 
!            chemistry & photolysis
! --------------------------------------------------------------------
item = stashcode_ukca_so4_sad - 1000*stashcode_ukca_chem_diag

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       so4_sa(:,:,:),                                                   &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags : error in copydiag_3d 50.256"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)


! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chem_diags

END MODULE ukca_chem_diags_mod
