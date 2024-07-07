! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE lbc_stashcode_mapping_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! Convert the STASH code for a prognostic field (in section 0 or tracers) into 
! the corresponding STASH code for an LBC field for that quanitity (in section
! 31).

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, PARAMETER :: num_ukca_tracers = 150
INTEGER, PARAMETER :: num_free_tracers = 150

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LBC_STASHCODE_MAPPING_MOD'

CONTAINS

SUBROUTINE lbc_stashcode_mapping(input_prognostic_stashcode,                 &
                                 output_lbc_stashcode)
USE um_stashcode_mod
IMPLICIT NONE
  
INTEGER, INTENT(IN) :: input_prognostic_stashcode
INTEGER, INTENT(OUT) :: output_lbc_stashcode
INTEGER :: mapping(99999)

CHARACTER(LEN=*), PARAMETER :: routinename = 'LBC_STASHCODE_MAPPING'
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: i
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mapping = 0

mapping(stashcode_orog) = stashcode_lbc_orog     
mapping(stashcode_u) = stashcode_lbc_u        
mapping(stashcode_v) = stashcode_lbc_v        
mapping(stashcode_w) = stashcode_lbc_w        
mapping(stashcode_rho) = stashcode_lbc_density  
mapping(stashcode_theta) = stashcode_lbc_theta    
mapping(stashcode_q) = stashcode_lbc_q        
mapping(stashcode_qcl) = stashcode_lbc_qcl      
mapping(stashcode_qcf) = stashcode_lbc_qcf      
mapping(stashcode_exner) = stashcode_lbc_exner    

mapping(stashcode_u_adv) = stashcode_lbc_u_adv    
mapping(stashcode_v_adv) = stashcode_lbc_v_adv    
mapping(stashcode_w_adv) = stashcode_lbc_w_adv    
mapping(stashcode_qcf2) = stashcode_lbc_qcf2     
mapping(stashcode_qrain) = stashcode_lbc_qrain    
mapping(stashcode_qgraup) = stashcode_lbc_qgraup   
mapping(stashcode_bulk_cf) = stashcode_lbc_cf_bulk  
mapping(stashcode_liquid_cf) = stashcode_lbc_cf_liquid
mapping(stashcode_frozen_cf) = stashcode_lbc_cf_frozen
mapping(stashcode_total_aero) = stashcode_lbc_murk     

mapping(stashcode_dust1_mmr) = stashcode_lbc_dust1_mmr
mapping(stashcode_dust2_mmr) = stashcode_lbc_dust2_mmr
mapping(stashcode_dust3_mmr) = stashcode_lbc_dust3_mmr
mapping(stashcode_dust4_mmr) = stashcode_lbc_dust4_mmr
mapping(stashcode_dust5_mmr) = stashcode_lbc_dust5_mmr
mapping(stashcode_dust6_mmr) = stashcode_lbc_dust6_mmr
mapping(stashcode_so2) = stashcode_lbc_so2      
mapping(stashcode_dms) = stashcode_lbc_dms      

mapping(stashcode_mmr_so4_aitken) = stashcode_lbc_so4_aitken
mapping(stashcode_mmr_so4_accum) = stashcode_lbc_so4_accu 
mapping(stashcode_mmr_so4_diss) = stashcode_lbc_so4_diss 
mapping(stashcode_mmr_nh3) = stashcode_lbc_nh3      
mapping(stashcode_mmr_bc_fr) = stashcode_lbc_soot_new 
mapping(stashcode_mmr_bc_ag) = stashcode_lbc_soot_agd 
mapping(stashcode_mmr_bc_cl) = stashcode_lbc_soot_cld 
mapping(stashcode_mmr_smoke_fr) = stashcode_lbc_bmass_new
mapping(stashcode_mmr_smoke_ag) = stashcode_lbc_bmass_agd
mapping(stashcode_mmr_smoke_cl) = stashcode_lbc_bmass_cld
mapping(stashcode_mmr_ocff_fr) = stashcode_lbc_ocff_new 
mapping(stashcode_mmr_ocff_ag) = stashcode_lbc_ocff_agd 
mapping(stashcode_mmr_ocff_cl) = stashcode_lbc_ocff_cld 
mapping(stashcode_mmr_nitr_acc) = stashcode_lbc_nitr_acc 
mapping(stashcode_mmr_nitr_diss) = stashcode_lbc_nitr_diss

DO i = 1, num_free_tracers
  mapping( i + stashcode_tracer_sec * 1000) = stashcode_lbc_free_tracer_1  &
                                               + i - 1
END DO

DO i = 1, num_ukca_tracers
  mapping( i + stashcode_ukca_sec * 1000) = stashcode_lbc_ukca_1 + i - 1
END DO

output_lbc_stashcode = mapping(input_prognostic_stashcode)
IF (output_lbc_stashcode == 0) THEN
  WRITE(cmessage,'(A,I5,A)') 'STASH code ', input_prognostic_stashcode,    &
                   ' does not have a section 31 counterpart'
  icode = 30
  CALL ereport(routinename,icode,cmessage)
  output_lbc_stashcode = input_prognostic_stashcode
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE lbc_stashcode_mapping

END MODULE lbc_stashcode_mapping_mod
