! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the elements required for radiative forcing

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

!- End of header

MODULE coradoca

USE missing_data_mod, ONLY: rmdi
USE umPrintMgr, ONLY: umPrint, umMessage
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Defaults for things to change in diagnostic 2nd call to radiation.
LOGICAL ::  c2c_all = .FALSE.
LOGICAL ::  c2c_wmg = .FALSE.
LOGICAL ::  c2c_o2 = .FALSE.
LOGICAL ::  c2c_o3 = .FALSE.
LOGICAL ::  c2c_co2 = .FALSE.
LOGICAL ::  c2c_n2o = .FALSE.
LOGICAL ::  c2c_ch4 = .FALSE.
LOGICAL ::  c2c_cfc11 = .FALSE.
LOGICAL ::  c2c_cfc12 = .FALSE.
LOGICAL ::  c2c_c113 = .FALSE.
LOGICAL ::  c2c_hcfc22 = .FALSE.
LOGICAL ::  c2c_hfc125 = .FALSE.
LOGICAL ::  c2c_hfc134 = .FALSE.
LOGICAL ::  c2c_aerosol = .FALSE.
LOGICAL ::  c2c_sulpc_d = .FALSE.
LOGICAL ::  c2c_seas_d = .FALSE.
LOGICAL ::  c2c_soot_d = .FALSE.
LOGICAL ::  c2c_bmb_d = .FALSE.
LOGICAL ::  c2c_ocff_d = .FALSE.
LOGICAL ::  c2c_nitr_d = .FALSE.
LOGICAL ::  c2c_dust_d = .FALSE.
LOGICAL ::  c2c_biog_d = .FALSE.
LOGICAL ::  c2c_ukca_d = .FALSE.
LOGICAL ::  c2c_land_s = .FALSE.
LOGICAL ::  c2c_easy_d = .FALSE.

! Gas mass mixing ratio scaling and additive factors for the
! diagnostic call:
REAL :: co2_mmr_scl     = rmdi
REAL :: co2_mmr_add     = rmdi
REAL :: n2o_mmr_scl     = rmdi
REAL :: n2o_mmr_add     = rmdi
REAL :: ch4_mmr_scl     = rmdi
REAL :: ch4_mmr_add     = rmdi
REAL :: o2_mmr_scl      = rmdi
REAL :: o2_mmr_add      = rmdi
REAL :: cfc11_mmr_scl   = rmdi
REAL :: cfc11_mmr_add   = rmdi
REAL :: cfc12_mmr_scl   = rmdi
REAL :: cfc12_mmr_add   = rmdi
REAL :: cfc113_mmr_scl  = rmdi
REAL :: cfc113_mmr_add  = rmdi
REAL :: hcfc22_mmr_scl  = rmdi
REAL :: hcfc22_mmr_add  = rmdi
REAL :: hfc125_mmr_scl  = rmdi
REAL :: hfc125_mmr_add  = rmdi
REAL :: hfc134a_mmr_scl = rmdi
REAL :: hfc134a_mmr_add = rmdi 


NAMELIST / radfcdia / c2c_o2, c2c_o3, c2c_co2, c2c_n2o, c2c_ch4,  &
  c2c_cfc11, c2c_cfc12, c2c_c113, c2c_hcfc22, c2c_hfc125,         &
  c2c_hfc134, c2c_aerosol, c2c_sulpc_d, c2c_seas_d, c2c_soot_d,   &
  c2c_bmb_d, c2c_ocff_d, c2c_land_s, c2c_all,                     &
  c2c_wmg, c2c_nitr_d, c2c_dust_d, c2c_biog_d, c2c_ukca_d,        &
  c2c_easy_d,                                                     &
  co2_mmr_scl, co2_mmr_add, n2o_mmr_scl, n2o_mmr_add,             &
  ch4_mmr_scl, ch4_mmr_add, o2_mmr_scl, o2_mmr_add,               &
  cfc11_mmr_scl, cfc11_mmr_add, cfc12_mmr_scl, cfc12_mmr_add,     &
  cfc113_mmr_scl, cfc113_mmr_add, hcfc22_mmr_scl, hcfc22_mmr_add, &
  hfc125_mmr_scl, hfc125_mmr_add, hfc134a_mmr_scl, hfc134a_mmr_add

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ----------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CORADOCA'

CONTAINS

! Subroutine to set the input values of the control structure.

SUBROUTINE coradoca_defaults

IMPLICIT NONE

REAL(KIND=jprb) :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CORADOCA_DEFAULTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( c2c_all ) THEN
  ! C_O3 = .TRUE.  ! Commented out till we have the ancillary
  ! Best made conditional on the alternate ozone file being specified
  ! if we don't hardwire the latter's name (in which case it'd always
  ! be safe to read from that file, & add nothing much to the cost of
  ! a 2nd call which is going to be made anyway).
  c2c_wmg = .TRUE.
  c2c_aerosol = .TRUE.
  c2c_land_s = .TRUE.
  ! This is used for surface albedo forcing
END IF
IF ( c2c_wmg ) THEN
  c2c_o2 = .TRUE.
  c2c_co2 = .TRUE.
  c2c_n2o = .TRUE.
  c2c_ch4 = .TRUE.
  c2c_cfc11 = .TRUE.
  c2c_cfc12 = .TRUE.
  c2c_c113 = .TRUE.
  c2c_hcfc22 = .TRUE.
  c2c_hfc125 = .TRUE.
  c2c_hfc134 = .TRUE.
END IF
IF ( c2c_aerosol ) THEN
  c2c_sulpc_d = .TRUE.
  c2c_seas_d = .TRUE.
  c2c_soot_d = .TRUE.
  c2c_bmb_d = .TRUE.
  c2c_ocff_d = .TRUE.
  c2c_nitr_d = .TRUE.
  c2c_dust_d = .TRUE.
  c2c_biog_d = .TRUE.
  c2c_ukca_d = .TRUE.
  c2c_easy_d = .TRUE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE coradoca_defaults

SUBROUTINE print_nlist_radfcdia()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RADFCDIA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist radfcdia', &
    src='coradoca')

WRITE(lineBuffer,*)' c2c_o2 = ',c2c_o2
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_o3 = ',c2c_o3
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_co2 = ',c2c_co2
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_n2o = ',c2c_n2o
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_ch4 = ',c2c_ch4
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_cfc11 = ',c2c_cfc11
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_cfc12 = ',c2c_cfc12
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_c113 = ',c2c_c113
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_hcfc22 = ',c2c_hcfc22
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_hfc125 = ',c2c_hfc125
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_hfc134 = ',c2c_hfc134
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_aerosol = ',c2c_aerosol
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_sulpc_d = ',c2c_sulpc_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_seas_d = ',c2c_seas_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_soot_d = ',c2c_soot_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_bmb_d = ',c2c_bmb_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_ocff_d = ',c2c_ocff_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_land_s = ',c2c_land_s
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_all = ',c2c_all
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_wmg = ',c2c_wmg
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_nitr_d = ',c2c_nitr_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_dust_d = ',c2c_dust_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_biog_d = ',c2c_biog_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,*)' c2c_ukca_d = ',c2c_ukca_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,L1)')   ' c2c_easy_d =  ',c2c_easy_d
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' co2_mmr_scl = ',co2_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' co2_mmr_add = ',co2_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' n2o_mmr_scl = ',n2o_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' n2o_mmr_add = ',n2o_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' ch4_mmr_scl = ',ch4_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' ch4_mmr_add = ',ch4_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' o2_mmr_scl = ',o2_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' o2_mmr_add = ',o2_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc11_mmr_scl = ',cfc11_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc11_mmr_add = ',cfc11_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc12_mmr_scl = ',cfc12_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc12_mmr_add = ',cfc12_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc113_mmr_scl = ',cfc113_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' cfc113_mmr_add = ',cfc113_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hcfc22_mmr_scl = ',hcfc22_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hcfc22_mmr_add = ',hcfc22_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hfc125_mmr_scl = ',hfc125_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hfc125_mmr_add = ',hfc125_mmr_add
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hfc134a_mmr_scl = ',hfc134a_mmr_scl
CALL umPrint(lineBuffer,src='coradoca')
WRITE(lineBuffer,'(A,E15.6)')' hfc134a_mmr_add = ',hfc134a_mmr_add
CALL umPrint(lineBuffer,src='coradoca')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='coradoca')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_radfcdia

SUBROUTINE read_nml_radfcdia(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RADFCDIA'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_real = 20
INTEGER, PARAMETER :: n_log = 25

TYPE my_namelist
  SEQUENCE
  REAL :: co2_mmr_scl
  REAL :: co2_mmr_add
  REAL :: n2o_mmr_scl
  REAL :: n2o_mmr_add
  REAL :: ch4_mmr_scl
  REAL :: ch4_mmr_add
  REAL :: o2_mmr_scl
  REAL :: o2_mmr_add
  REAL :: cfc11_mmr_scl
  REAL :: cfc11_mmr_add
  REAL :: cfc12_mmr_scl
  REAL :: cfc12_mmr_add
  REAL :: cfc113_mmr_scl
  REAL :: cfc113_mmr_add
  REAL :: hcfc22_mmr_scl
  REAL :: hcfc22_mmr_add
  REAL :: hfc125_mmr_scl
  REAL :: hfc125_mmr_add
  REAL :: hfc134a_mmr_scl
  REAL :: hfc134a_mmr_add
  LOGICAL :: c2c_o2
  LOGICAL :: c2c_o3
  LOGICAL :: c2c_co2
  LOGICAL :: c2c_n2o
  LOGICAL :: c2c_ch4
  LOGICAL :: c2c_cfc11
  LOGICAL :: c2c_cfc12
  LOGICAL :: c2c_c113
  LOGICAL :: c2c_hcfc22
  LOGICAL :: c2c_hfc125
  LOGICAL :: c2c_hfc134
  LOGICAL :: c2c_aerosol
  LOGICAL :: c2c_sulpc_d
  LOGICAL :: c2c_seas_d
  LOGICAL :: c2c_soot_d
  LOGICAL :: c2c_bmb_d
  LOGICAL :: c2c_ocff_d
  LOGICAL :: c2c_land_s
  LOGICAL :: c2c_all
  LOGICAL :: c2c_wmg
  LOGICAL :: c2c_nitr_d
  LOGICAL :: c2c_dust_d
  LOGICAL :: c2c_biog_d
  LOGICAL :: c2c_ukca_d
  LOGICAL :: c2c_easy_d
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_log_in=n_log, &
                    n_real_in=n_real)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=radfcdia, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RADFCDIA", iomessage)

  my_nml % co2_mmr_scl     = co2_mmr_scl    
  my_nml % co2_mmr_add     = co2_mmr_add    
  my_nml % n2o_mmr_scl     = n2o_mmr_scl    
  my_nml % n2o_mmr_add     = n2o_mmr_add    
  my_nml % ch4_mmr_scl     = ch4_mmr_scl    
  my_nml % ch4_mmr_add     = ch4_mmr_add    
  my_nml % o2_mmr_scl      = o2_mmr_scl     
  my_nml % o2_mmr_add      = o2_mmr_add     
  my_nml % cfc11_mmr_scl   = cfc11_mmr_scl  
  my_nml % cfc11_mmr_add   = cfc11_mmr_add  
  my_nml % cfc12_mmr_scl   = cfc12_mmr_scl  
  my_nml % cfc12_mmr_add   = cfc12_mmr_add  
  my_nml % cfc113_mmr_scl  = cfc113_mmr_scl 
  my_nml % cfc113_mmr_add  = cfc113_mmr_add 
  my_nml % hcfc22_mmr_scl  = hcfc22_mmr_scl 
  my_nml % hcfc22_mmr_add  = hcfc22_mmr_add 
  my_nml % hfc125_mmr_scl  = hfc125_mmr_scl 
  my_nml % hfc125_mmr_add  = hfc125_mmr_add 
  my_nml % hfc134a_mmr_scl = hfc134a_mmr_scl
  my_nml % hfc134a_mmr_add = hfc134a_mmr_add
  ! end of reals
  my_nml % c2c_o2      = c2c_o2
  my_nml % c2c_o3      = c2c_o3
  my_nml % c2c_co2     = c2c_co2
  my_nml % c2c_n2o     = c2c_n2o
  my_nml % c2c_ch4     = c2c_ch4
  my_nml % c2c_cfc11   = c2c_cfc11
  my_nml % c2c_cfc12   = c2c_cfc12
  my_nml % c2c_c113    = c2c_c113
  my_nml % c2c_hcfc22  = c2c_hcfc22
  my_nml % c2c_hfc125  = c2c_hfc125
  my_nml % c2c_hfc134  = c2c_hfc134
  my_nml % c2c_aerosol = c2c_aerosol
  my_nml % c2c_sulpc_d = c2c_sulpc_d
  my_nml % c2c_seas_d  = c2c_seas_d
  my_nml % c2c_soot_d  = c2c_soot_d
  my_nml % c2c_bmb_d   = c2c_bmb_d
  my_nml % c2c_ocff_d  = c2c_ocff_d
  my_nml % c2c_land_s  = c2c_land_s
  my_nml % c2c_all     = c2c_all
  my_nml % c2c_wmg     = c2c_wmg
  my_nml % c2c_nitr_d  = c2c_nitr_d
  my_nml % c2c_dust_d  = c2c_dust_d
  my_nml % c2c_biog_d  = c2c_biog_d
  my_nml % c2c_ukca_d  = c2c_ukca_d
  my_nml % c2c_easy_d  = c2c_easy_d

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  co2_mmr_scl     = my_nml % co2_mmr_scl    
  co2_mmr_add     = my_nml % co2_mmr_add    
  n2o_mmr_scl     = my_nml % n2o_mmr_scl    
  n2o_mmr_add     = my_nml % n2o_mmr_add    
  ch4_mmr_scl     = my_nml % ch4_mmr_scl    
  ch4_mmr_add     = my_nml % ch4_mmr_add    
  o2_mmr_scl      = my_nml % o2_mmr_scl     
  o2_mmr_add      = my_nml % o2_mmr_add     
  cfc11_mmr_scl   = my_nml % cfc11_mmr_scl  
  cfc11_mmr_add   = my_nml % cfc11_mmr_add  
  cfc12_mmr_scl   = my_nml % cfc12_mmr_scl  
  cfc12_mmr_add   = my_nml % cfc12_mmr_add  
  cfc113_mmr_scl  = my_nml % cfc113_mmr_scl 
  cfc113_mmr_add  = my_nml % cfc113_mmr_add 
  hcfc22_mmr_scl  = my_nml % hcfc22_mmr_scl 
  hcfc22_mmr_add  = my_nml % hcfc22_mmr_add 
  hfc125_mmr_scl  = my_nml % hfc125_mmr_scl 
  hfc125_mmr_add  = my_nml % hfc125_mmr_add 
  hfc134a_mmr_scl = my_nml % hfc134a_mmr_scl
  hfc134a_mmr_add = my_nml % hfc134a_mmr_add
  ! end of reals
  c2c_o2      = my_nml % c2c_o2
  c2c_o3      = my_nml % c2c_o3
  c2c_co2     = my_nml % c2c_co2
  c2c_n2o     = my_nml % c2c_n2o
  c2c_ch4     = my_nml % c2c_ch4
  c2c_cfc11   = my_nml % c2c_cfc11
  c2c_cfc12   = my_nml % c2c_cfc12
  c2c_c113    = my_nml % c2c_c113
  c2c_hcfc22  = my_nml % c2c_hcfc22
  c2c_hfc125  = my_nml % c2c_hfc125
  c2c_hfc134  = my_nml % c2c_hfc134
  c2c_aerosol = my_nml % c2c_aerosol
  c2c_sulpc_d = my_nml % c2c_sulpc_d
  c2c_seas_d  = my_nml % c2c_seas_d
  c2c_soot_d  = my_nml % c2c_soot_d
  c2c_bmb_d   = my_nml % c2c_bmb_d
  c2c_ocff_d  = my_nml % c2c_ocff_d
  c2c_land_s  = my_nml % c2c_land_s
  c2c_all     = my_nml % c2c_all
  c2c_wmg     = my_nml % c2c_wmg
  c2c_nitr_d  = my_nml % c2c_nitr_d
  c2c_dust_d  = my_nml % c2c_dust_d
  c2c_biog_d  = my_nml % c2c_biog_d
  c2c_ukca_d  = my_nml % c2c_ukca_d
  c2c_easy_d  = my_nml % c2c_easy_d

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_radfcdia

END MODULE coradoca
