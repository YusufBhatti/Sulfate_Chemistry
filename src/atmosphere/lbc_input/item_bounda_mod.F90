! *****************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*********************************

! Module to set up array containing STASH codes to be read in from LBC
! file, and the corresponding subscripts of the array of lookup headers.
!
MODULE item_bounda_mod

USE ukca_tracer_stash, ONLY: ukca_tr_lbc_stashitem

IMPLICIT NONE

! Description:
!   Copy of the module used in the branch dealing with balanced LBCs:
!   vn8.1_makeBC_reduced_LBCs_balance
!
!   Not necessary for the advected winds case, but using this module
!   will make merging the branches easier.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: lbc_input
!
! Code description:
!   Language: Fortran 95.

! Module variables containing lookup subscript values
! Used by READ_ATMOS_LBCS
INTEGER :: lookup_subscript_u
INTEGER :: lookup_subscript_v
INTEGER :: lookup_subscript_w
INTEGER :: lookup_subscript_rho
INTEGER :: lookup_subscript_theta
INTEGER :: lookup_subscript_q
INTEGER :: lookup_subscript_qcl
INTEGER :: lookup_subscript_qcf
INTEGER :: lookup_subscript_exner
INTEGER :: lookup_subscript_u_adv
INTEGER :: lookup_subscript_v_adv
INTEGER :: lookup_subscript_w_adv

! Stores number of "core" LBCs being used (i.e. before aerosols etc.)
INTEGER, SAVE :: lbc_num_stored

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ITEM_BOUNDA_MOD'

CONTAINS

SUBROUTINE assign_item_bounda(item_bounda)

! Description:
!   See Module description

USE dust_parameters_mod, ONLY:                                    &
     l_dust_div1_lbc,   l_dust_div2_lbc,   l_dust_div3_lbc,       &
     l_dust_div4_lbc,   l_dust_div5_lbc,   l_dust_div6_lbc

USE run_aerosol_mod,  ONLY:                                  &
     l_so2_lbc,        l_dms_lbc,       l_so4_aitken_lbc,    &
     l_so4_accu_lbc,   l_so4_diss_lbc,  l_nh3_lbc,           &
     l_soot_new_lbc,   l_soot_agd_lbc,  l_soot_cld_lbc,      &
     l_bmass_new_lbc,  l_bmass_agd_lbc, l_bmass_cld_lbc,     &
     l_ocff_new_lbc,   l_ocff_agd_lbc,  l_ocff_cld_lbc,      &
     l_nitr_acc_lbc,   l_nitr_diss_lbc

USE mphys_inputs_mod, ONLY:                                       &
     l_mcr_qcf2_lbc,  l_mcr_qgraup_lbc, l_mcr_qrain_lbc

USE murk_inputs_mod,  ONLY: l_murk_lbc

USE cloud_inputs_mod, ONLY: l_pc2_lbc

USE free_tracers_inputs_mod, ONLY: A_TR_lbc_stashitem

USE lbc_read_data_mod, ONLY:  l_int_uvw_lbc, L_old_lbc_file

USE lbc_mod,         ONLY: rim_lookupsa
USE umPrintMgr
USE parkind1,        ONLY: jpim, jprb
USE um_parcore,      ONLY: mype
USE yomhook,         ONLY: lhook, dr_hook

USE nlsizes_namelist_mod, ONLY: &
    tr_lbc_ukca, tr_lbc_vars

IMPLICIT NONE


INTEGER, PARAMETER :: sect36 = 36000
INTEGER, PARAMETER :: sect37 = 37000

! Subroutine arguments
INTEGER, INTENT(INOUT) :: item_bounda(rim_lookupsa) ! Array of STASH codes of LBC
                                                    ! input fields

! Local variables
INTEGER :: i, j, lbc_num

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASSIGN_ITEM_BOUNDA'

!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

item_bounda(1)  = 31001 ! Orography
item_bounda(2)  = 31002 ! U
item_bounda(3)  = 31003 ! V

lookup_subscript_u     = 1
lookup_subscript_v     = 2

IF (L_old_lbc_file) THEN ! Use old LBC file structure
  item_bounda(4)  = 31004 ! W
  item_bounda(5)  = 31005 ! Density
  item_bounda(6)  = 31006 ! Potential temperature
  item_bounda(7)  = 31007 ! Specific humidity
  item_bounda(8)  = 31008 ! QCL
  item_bounda(9)  = 31009 ! QCF
  item_bounda(10) = 31010 ! Exner
  item_bounda(11) = 31011 ! U_Adv
  item_bounda(12) = 31012 ! V_Adv
  item_bounda(13) = 31013 ! W Adv

  lookup_subscript_w     = 3
  lookup_subscript_rho   = 4
  lookup_subscript_theta = 5
  lookup_subscript_q     = 6
  lookup_subscript_qcl   = 7
  lookup_subscript_qcf   = 8
  lookup_subscript_exner = 9
  lookup_subscript_u_adv = 10
  lookup_subscript_v_adv = 11
  lookup_subscript_w_adv = 12

  lbc_num = 13

ELSE ! New LBC file - advected winds optional
  item_bounda(4)  = 31004 ! W
  item_bounda(5)  = 31005 ! Density
  item_bounda(6)  = 31006 ! Potential temperature
  item_bounda(7)  = 31007 ! Specific humidity
  item_bounda(8)  = 31008 ! QCL
  item_bounda(9)  = 31009 ! QCF
  item_bounda(10) = 31010 ! Exner

  lookup_subscript_w     = 3
  lookup_subscript_rho   = 4
  lookup_subscript_theta = 5
  lookup_subscript_q     = 6
  lookup_subscript_qcl   = 7
  lookup_subscript_qcf   = 8
  lookup_subscript_exner = 9

  lbc_num = 10 ! Fewer compulsory fields in new LBC files
               ! Graham's original code had 9 but why?

  IF (.NOT. L_int_uvw_lbc) THEN ! If not interpolating, require advected winds
    item_bounda(lbc_num+1) = 31011 ! U advected
    item_bounda(lbc_num+2) = 31012 ! V advected
    item_bounda(lbc_num+3) = 31013 ! W advected

    lookup_subscript_u_adv = lbc_num
    lookup_subscript_v_adv = lbc_num+1
    lookup_subscript_w_adv = lbc_num+2

    lbc_num = lbc_num+3
  END IF

END IF

lbc_num_stored = lbc_num

! Include additional microphysics lbcs if active

IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31014
END IF
IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31015
END IF
IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31016
END IF

! Setup for additional cloud fraction lbcs IF active
IF (L_pc2_lbc) THEN
  lbc_num = lbc_num+1
  item_bounda(lbc_num) = 31017 ! cf_bulk
  lbc_num = lbc_num+1
  item_bounda(lbc_num) = 31018 ! cf_liquid
  lbc_num = lbc_num+1
  item_bounda(lbc_num) = 31019 ! cf_frozen
END IF

! Include murk aerosol lbcs IF in input lbc file
IF (L_murk_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31020 ! murk aerosol
END IF

! Include dust lbcs IF in input lbc file
IF (L_dust_div1_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31023 ! dust_div1
END IF
! Include dust lbcs IF in input lbc file
IF (L_dust_div2_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31024 ! dust_div2
END IF
! Include dust lbcs IF in input lbc file
IF (L_dust_div3_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31025 ! dust_div3
END IF
! Include dust lbcs IF in input lbc file
IF (L_dust_div4_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31026 ! dust_div4
END IF
! Include dust lbcs IF in input lbc file
IF (L_dust_div5_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31027 ! dust_div5
END IF
! Include dust lbcs IF in input lbc file
IF (L_dust_div6_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31028 ! dust_div6
END IF

! Include so2 lbcs IF in input lbc file
IF (L_so2_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31029 ! so2
END IF

! Include dms lbcs IF in input lbc file
IF (L_dms_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31030 ! dms
END IF

! Include so4_aitken lbcs IF in input lbc file
IF (L_so4_aitken_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31031 ! so4_aitken
END IF
! Include so4_accu lbcs IF in input lbc file
IF (L_so4_accu_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31032 ! so4_accu
END IF
! Include so4_diss lbcs IF in input lbc file
IF (L_so4_diss_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31033 ! so4_diss
END IF

! Include nh3 lbcs IF in input lbc file
IF (L_nh3_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31035 ! nh3
END IF

! Include soot_new lbcs IF in input lbc file
IF (L_soot_new_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31036 ! soot_new
END IF
! Include soot_agd lbcs IF in input lbc file
IF (L_soot_agd_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31037 ! soot_agd
END IF
! Include soot_cld lbcs IF in input lbc file
IF (L_soot_cld_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31038 ! soot_cld
END IF

! Include bmass_new lbcs IF in input lbc file
IF (L_bmass_new_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31039 ! bmass_new
END IF
! Include bmass_agd lbcs IF in input lbc file
IF (L_bmass_agd_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31040 ! bmass_agd
END IF
! Include bmass_cld lbcs IF in input lbc file
IF (L_bmass_cld_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31041 ! bmass_cld
END IF

! Include ocff_new lbcs IF in input lbc file
IF (L_ocff_new_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31042 ! ocff_new
END IF
! Include ocff_agd lbcs IF in input lbc file
IF (L_ocff_agd_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31043 ! ocff_agd
END IF
! Include ocff_cld lbcs IF in input lbc file
IF (L_ocff_cld_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31044 ! ocff_cld
END IF

! Include nitr_acc lbcs IF in input lbc file
IF (L_nitr_acc_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31045 ! nitr_new
END IF
! Include nitr_diss lbcs IF in input lbc file
IF (L_nitr_diss_lbc) THEN
  lbc_num = lbc_num + 1
  item_bounda(lbc_num) = 31046 ! nitr_diss
END IF


! Free tracer LBCs
IF (tr_lbc_vars > 0) THEN

  ! Set item_bounda to STASH item codes read in from the
  ! file so that the user only needs to specify the total
  ! number of tracer lbcs

  DO i = 1, tr_lbc_vars
    item_bounda(lbc_num+i)=sect36 + A_TR_lbc_stashitem(i)
    WRITE(umMessage,'(A,I8,I8)') 'INBOUNDA (free):',lbc_num+i,       &
        item_bounda(lbc_num+i)
    CALL umPrint(umMessage)
  END DO
  lbc_num = lbc_num + tr_lbc_vars

END IF

IF (tr_lbc_ukca > 0) THEN

  ! Set item_bounda to STASH item codes read in from the
  ! file so that the user only needs to specify the total
  ! number of tracer lbcs
  DO i = 1, tr_lbc_ukca
    item_bounda(lbc_num+i)=sect37 + UKCA_TR_lbc_stashitem(i)
    WRITE(umMessage,'(A, I8, I8)') 'INBOUNDA (ukca):',lbc_num+i,     &
        item_bounda(lbc_num+i)
    CALL umPrint(umMessage)
  END DO

END IF

IF (PrintStatus == PrStatus_Diag .AND. mype == 0) THEN
  WRITE(6,'(A,I2)')'lbc_num=',lbc_num
  WRITE(6,*)'ITEM_BOUNDA=',(item_bounda(j), j=1,rim_lookupsa)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE assign_item_bounda
END MODULE item_bounda_mod

