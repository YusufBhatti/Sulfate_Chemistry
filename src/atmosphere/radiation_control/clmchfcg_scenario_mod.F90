! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing tunable parameters used for histories/scenarios
!   of climate change forcings (clm ch fcg).
!
MODULE clmchfcg_scenario_mod

!
! Description:
!   This module contains declarations for tunable parameters
!   used for histories/scenarios of climate change forcings.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Radiation Control
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards UMDP 3 vn8.2
!
USE missing_data_mod, ONLY: imdi, rmdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Number of well-mixed greenhouse gases
INTEGER, PARAMETER :: Nwmghg   = 10

! Number of sulphate loading patterns
INTEGER, PARAMETER :: Nsulpat  = 2

! Maximum length of scenarios
INTEGER, PARAMETER :: LenScen  = 500

! Number of such scenarios, made up of:
INTEGER, PARAMETER :: Nclmfcgs = Nwmghg + Nsulpat

! Indices indicating which scenario corresponds to which forcing:
INTEGER, PARAMETER :: s_co2     = 1
INTEGER, PARAMETER :: s_ch4     = 2
INTEGER, PARAMETER :: s_n2o     = 3
INTEGER, PARAMETER :: s_cfc11   = 4
INTEGER, PARAMETER :: s_cfc12   = 5
INTEGER, PARAMETER :: s_so4     = 6
INTEGER, PARAMETER :: s_cfc113  = 8
INTEGER, PARAMETER :: s_hcfc22  = 9
INTEGER, PARAMETER :: s_hfc125  = 10
INTEGER, PARAMETER :: s_hfc134a = 11
INTEGER, PARAMETER :: s_cfc114  = 12

!  Carbon dioxide (CO2), methane (CH4), nitrous oxide (N2O),
!  trichlorofluoromethane (CCl3F, "CFC-11"),
!  dichlorodifluoromethane (CCl2F2, "CFC-12"), and then the first
!  HadCM2-style anthropogenic sulphate loading pattern - these
!  come at the end as their number in principle may vary.

! Switch for time varying GHG
LOGICAL     :: L_ClmChFcg = .FALSE.

! Switch to make rates apply continuously rather than a step change at
!   the beginning of a year
LOGICAL     :: L_Cts_Fcg_Rates = .FALSE.

! Years at which a rate or level is specified
INTEGER     :: Clim_Fcg_Years(LenScen,Nclmfcgs) = imdi
INTEGER     :: Clim_Fcg_Years_CO2(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_CH4(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_N2O(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_CFC11(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_CFC12(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_SO4(2*LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_CFC113(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_HCFC22(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_HFC125(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_HFC134A(LenScen) = imdi
INTEGER     :: Clim_Fcg_Years_CFC114(LenScen) = imdi

! Number of such years, for each forcing
INTEGER     :: Clim_Fcg_NYears(Nclmfcgs) = imdi
INTEGER     :: Clim_Fcg_NYears_CO2 = imdi
INTEGER     :: Clim_Fcg_NYears_CH4 = imdi
INTEGER     :: Clim_Fcg_NYears_N2O = imdi
INTEGER     :: Clim_Fcg_NYears_CFC11 = imdi
INTEGER     :: Clim_Fcg_NYears_CFC12 = imdi
INTEGER     :: Clim_Fcg_NYears_SO4(2) = imdi
INTEGER     :: Clim_Fcg_NYears_CFC113 = imdi
INTEGER     :: Clim_Fcg_NYears_HCFC22 = imdi
INTEGER     :: Clim_Fcg_NYears_HFC125 = imdi
INTEGER     :: Clim_Fcg_NYears_HFC134A = imdi
INTEGER     :: Clim_Fcg_NYears_CFC114 = imdi

! Values, or rates of increase, for the designated years.
!  See GAS_CALC  in the gui/namelist for details.
REAL        :: Clim_Fcg_Levls(LenScen,Nclmfcgs) = rmdi
REAL        :: Clim_Fcg_Levls_CO2(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_CH4(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_N2O(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_CFC11(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_CFC12(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_SO4(2*LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_CFC113(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_HCFC22(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_HFC125(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_HFC134A(LenScen) = rmdi
REAL        :: Clim_Fcg_Levls_CFC114(LenScen) = rmdi

REAL        :: Clim_Fcg_Rates(LenScen,Nclmfcgs) = rmdi
REAL        :: Clim_Fcg_Rates_CO2(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_CH4(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_N2O(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_CFC11(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_CFC12(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_SO4(2*LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_CFC113(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_HCFC22(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_HFC125(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_HFC134A(LenScen) = rmdi
REAL        :: Clim_Fcg_Rates_CFC114(LenScen) = rmdi

NAMELIST / clmchfcg / Clim_Fcg_NYears_CO2, Clim_Fcg_NYears_CH4, &
                      Clim_Fcg_NYears_N2O, Clim_Fcg_NYears_CFC11, &
                      Clim_Fcg_NYears_CFC12, Clim_Fcg_NYears_SO4, &
                      Clim_Fcg_NYears_CFC113, Clim_Fcg_NYears_HCFC22, &
                      Clim_Fcg_NYears_HFC125, Clim_Fcg_NYears_HFC134A, &
                      Clim_Fcg_NYears_CFC114, &
                      Clim_Fcg_Years_CO2, Clim_Fcg_Years_CH4, &
                      Clim_Fcg_Years_N2O, Clim_Fcg_Years_CFC11, &
                      Clim_Fcg_Years_CFC12, Clim_Fcg_Years_SO4, &
                      Clim_Fcg_Years_CFC113, Clim_Fcg_Years_HCFC22, &
                      Clim_Fcg_Years_HFC125, Clim_Fcg_Years_HFC134A, &
                      Clim_Fcg_Years_CFC114, &
                      Clim_Fcg_Levls_CO2, Clim_Fcg_Levls_CH4, &
                      Clim_Fcg_Levls_N2O, Clim_Fcg_Levls_CFC11, &
                      Clim_Fcg_Levls_CFC12, Clim_Fcg_Levls_SO4, &
                      Clim_Fcg_Levls_CFC113, Clim_Fcg_Levls_HCFC22, &
                      Clim_Fcg_Levls_HFC125, Clim_Fcg_Levls_HFC134A, &
                      Clim_Fcg_Levls_CFC114, &
                      Clim_Fcg_Rates_CO2, Clim_Fcg_Rates_CH4, &
                      Clim_Fcg_Rates_N2O, Clim_Fcg_Rates_CFC11, &
                      Clim_Fcg_Rates_CFC12, Clim_Fcg_Rates_SO4, &
                      Clim_Fcg_Rates_CFC113, Clim_Fcg_Rates_HCFC22, &
                      Clim_Fcg_Rates_HFC125, Clim_Fcg_Rates_HFC134A, &
                      Clim_Fcg_Rates_CFC114, &
                      L_ClmChFcg, L_Cts_Fcg_Rates

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! ----------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CLMCHFCG_SCENARIO_MOD'

CONTAINS

SUBROUTINE read_nml_clmchfcg(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_CLMCHFCG'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = Nclmfcgs + LenScen * Nclmfcgs
INTEGER, PARAMETER :: n_real = 2 * LenScen * Nclmfcgs
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: Clim_Fcg_NYears(Nclmfcgs)
  INTEGER :: Clim_Fcg_Years(LenScen,Nclmfcgs)
  REAL :: Clim_Fcg_Levls(LenScen,Nclmfcgs)
  REAL :: Clim_Fcg_Rates(LenScen,Nclmfcgs)
  LOGICAL :: L_ClmChFcg
  LOGICAL :: L_Cts_Fcg_Rates
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=clmchfcg, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist CLMCHFCG", iomessage)

  ! Populate the grouped-types ready for broadcasting
  Clim_Fcg_Years(:,S_CO2)     = Clim_Fcg_Years_CO2
  Clim_Fcg_Years(:,S_CH4)     = Clim_Fcg_Years_CH4
  Clim_Fcg_Years(:,S_N2O)     = Clim_Fcg_Years_N2O
  Clim_Fcg_Years(:,S_CFC11)   = Clim_Fcg_Years_CFC11
  Clim_Fcg_Years(:,S_CFC12)   = Clim_Fcg_Years_CFC12
  Clim_Fcg_Years(:,S_SO4)     = Clim_Fcg_Years_SO4(1:LenScen)
  Clim_Fcg_Years(:,S_CFC113)  = Clim_Fcg_Years_CFC113
  Clim_Fcg_Years(:,S_HCFC22)  = Clim_Fcg_Years_HCFC22
  Clim_Fcg_Years(:,S_HFC125)  = Clim_Fcg_Years_HFC125
  Clim_Fcg_Years(:,S_HFC134A) = Clim_Fcg_Years_HFC134A
  Clim_Fcg_Years(:,S_CFC114)  = Clim_Fcg_Years_CFC114

  Clim_Fcg_NYears(S_CO2)      = Clim_Fcg_NYears_CO2
  Clim_Fcg_NYears(S_CH4)      = Clim_Fcg_NYears_CH4
  Clim_Fcg_NYears(S_N2O)      = Clim_Fcg_NYears_N2O
  Clim_Fcg_NYears(S_CFC11)    = Clim_Fcg_NYears_CFC11
  Clim_Fcg_NYears(S_CFC12)    = Clim_Fcg_NYears_CFC12
  Clim_Fcg_NYears(S_SO4)      = Clim_Fcg_NYears_SO4(1)
  Clim_Fcg_NYears(S_CFC113)   = Clim_Fcg_NYears_CFC113
  Clim_Fcg_NYears(S_HCFC22)   = Clim_Fcg_NYears_HCFC22
  Clim_Fcg_NYears(S_HFC125)   = Clim_Fcg_NYears_HFC125
  Clim_Fcg_NYears(S_HFC134A)  = Clim_Fcg_NYears_HFC134A
  Clim_Fcg_NYears(S_CFC114)   = Clim_Fcg_NYears_CFC114

  Clim_Fcg_Levls(:,S_CO2)     = Clim_Fcg_Levls_CO2
  Clim_Fcg_Levls(:,S_CH4)     = Clim_Fcg_Levls_CH4
  Clim_Fcg_Levls(:,S_N2O)     = Clim_Fcg_Levls_N2O
  Clim_Fcg_Levls(:,S_CFC11)   = Clim_Fcg_Levls_CFC11
  Clim_Fcg_Levls(:,S_CFC12)   = Clim_Fcg_Levls_CFC12
  Clim_Fcg_Levls(:,S_SO4)     = Clim_Fcg_Levls_SO4(1:LenScen)
  Clim_Fcg_Levls(:,S_CFC113)  = Clim_Fcg_Levls_CFC113
  Clim_Fcg_Levls(:,S_HCFC22)  = Clim_Fcg_Levls_HCFC22
  Clim_Fcg_Levls(:,S_HFC125)  = Clim_Fcg_Levls_HFC125
  Clim_Fcg_Levls(:,S_HFC134A) = Clim_Fcg_Levls_HFC134A
  Clim_Fcg_Levls(:,S_CFC114)  = Clim_Fcg_Levls_CFC114

  Clim_Fcg_Rates(:,S_CO2)     = Clim_Fcg_Rates_CO2
  Clim_Fcg_Rates(:,S_CH4)     = Clim_Fcg_Rates_CH4
  Clim_Fcg_Rates(:,S_N2O)     = Clim_Fcg_Rates_N2O
  Clim_Fcg_Rates(:,S_CFC11)   = Clim_Fcg_Rates_CFC11
  Clim_Fcg_Rates(:,S_CFC12)   = Clim_Fcg_Rates_CFC12
  Clim_Fcg_Rates(:,S_SO4)     = Clim_Fcg_Rates_SO4(1:LenScen)
  Clim_Fcg_Rates(:,S_CFC113)  = Clim_Fcg_Rates_CFC113
  Clim_Fcg_Rates(:,S_HCFC22)  = Clim_Fcg_Rates_HCFC22
  Clim_Fcg_Rates(:,S_HFC125)  = Clim_Fcg_Rates_HFC125
  Clim_Fcg_Rates(:,S_HFC134A) = Clim_Fcg_Rates_HFC134A
  Clim_Fcg_Rates(:,S_CFC114)  = Clim_Fcg_Rates_CFC114

  ! Copy to the broadcast type
  my_nml % Clim_Fcg_Years  = Clim_Fcg_Years
  my_nml % Clim_Fcg_NYears = Clim_Fcg_NYears
  my_nml % Clim_Fcg_Levls  = Clim_Fcg_Levls
  my_nml % Clim_Fcg_Rates  = Clim_Fcg_Rates
  my_nml % L_ClmChFcg      = L_ClmChFcg
  my_nml % L_Cts_Fcg_Rates = L_Cts_Fcg_Rates

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  Clim_Fcg_NYears = my_nml % Clim_Fcg_NYears
  Clim_Fcg_Years  = my_nml % Clim_Fcg_Years
  Clim_Fcg_Levls  = my_nml % Clim_Fcg_Levls
  Clim_Fcg_Rates  = my_nml % Clim_Fcg_Rates
  L_ClmChFcg      = my_nml % L_ClmChFcg
  L_Cts_Fcg_Rates = my_nml % L_Cts_Fcg_Rates

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE read_nml_clmchfcg

!----------------------------------------------------

SUBROUTINE clmchfcg_rates

IMPLICIT NONE

INTEGER :: j, jj

IF ( L_clmchfcg ) THEN
   ! Convert rates from percent to multiplicative factors:
  DO j = 1, nclmfcgs
    ! This is a null loop, as it should be, if clim_fcg_nyears=0
    DO jj = 1, clim_fcg_nyears(j)
      IF ( clim_fcg_rates(jj,j)  >   -100.0 ) THEN
        clim_fcg_rates(jj,j) = 1.0 + 0.01 * clim_fcg_rates(jj,j)
      END IF
    END DO
  END DO
ELSE
  ! If the namelist is not to be read, set number of designated
  !   years to zero for all possible forcings, as this may be
  !   used to test if this system is being used for each forcing.
  DO j = 1, nclmfcgs
    clim_fcg_nyears(j) = 0
  END DO
END IF !  L_clmchfcg


END SUBROUTINE clmchfcg_rates

!----------------------------------------------------

SUBROUTINE print_nlist_clmchfcg()

USE umPrintMgr, ONLY : umPrint

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_CLMCHFCG'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist clmchfcg', &
     src='clmchfcg_scenario_mod')

WRITE(lineBuffer,FMT='(A,L1)') ' L_ClmChFcg = ',L_ClmChFcg
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

WRITE(lineBuffer,FMT='(A,L1)') ' L_Cts_Fcg_Rates = ',L_Cts_Fcg_Rates
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

WRITE(lineBuffer,*) ' Clim_Fcg_Years_CO2 = ',Clim_Fcg_Years_CO2
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_CH4 = ',Clim_Fcg_Years_CH4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_N2O = ',Clim_Fcg_Years_N2O
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_CFC11 = ',Clim_Fcg_Years_CFC11
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_CFC12 = ',Clim_Fcg_Years_CFC12
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_SO4 = ',Clim_Fcg_Years_SO4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_CFC113 = ',Clim_Fcg_Years_CFC113
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_HCFC22 = ',Clim_Fcg_Years_HCFC22
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_HFC125 = ',Clim_Fcg_Years_HFC125
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_HFC134A = ',Clim_Fcg_Years_HFC134A
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Years_CFC114 = ',Clim_Fcg_Years_CFC114
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CO2 = ',Clim_Fcg_NYears_CO2
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CH4 = ',Clim_Fcg_NYears_CH4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_N2O = ',Clim_Fcg_NYears_N2O
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CFC11 = ',Clim_Fcg_NYears_CFC11
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CFC12 = ',Clim_Fcg_NYears_CFC12
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_SO4 = ',Clim_Fcg_NYears_SO4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CFC113 = ',Clim_Fcg_NYears_CFC113
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_HCFC22 = ',Clim_Fcg_NYears_HCFC22
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_HFC125 = ',Clim_Fcg_NYears_HFC125
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_HFC134A = ',Clim_Fcg_NYears_HFC134A
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_NYears_CFC114 = ',Clim_Fcg_NYears_CFC114
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CO2 = ',Clim_Fcg_Levls_CO2
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CH4 = ',Clim_Fcg_Levls_CH4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_N2O = ',Clim_Fcg_Levls_N2O
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CFC11 = ',Clim_Fcg_Levls_CFC11
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CFC12 = ',Clim_Fcg_Levls_CFC12
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_SO4 = ',Clim_Fcg_Levls_SO4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CFC113 = ',Clim_Fcg_Levls_CFC113
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_HCFC22 = ',Clim_Fcg_Levls_HCFC22
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_HFC125 = ',Clim_Fcg_Levls_HFC125
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_HFC134A = ',Clim_Fcg_Levls_HFC134A
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Levls_CFC114 = ',Clim_Fcg_Levls_CFC114
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CO2 = ',Clim_Fcg_Rates_CO2
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CH4 = ',Clim_Fcg_Rates_CH4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_N2O = ',Clim_Fcg_Rates_N2O
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CFC11 = ',Clim_Fcg_Rates_CFC11
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CFC12 = ',Clim_Fcg_Rates_CFC12
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_SO4 = ',Clim_Fcg_Rates_SO4
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CFC113 = ',Clim_Fcg_Rates_CFC113
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_HCFC22 = ',Clim_Fcg_Rates_HCFC22
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_HFC125 = ',Clim_Fcg_Rates_HFC125
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_HFC134A = ',Clim_Fcg_Rates_HFC134A
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')
WRITE(lineBuffer,*) ' Clim_Fcg_Rates_CFC114 = ',Clim_Fcg_Rates_CFC114
CALL umPrint(lineBuffer,src='clmchfcg_scenario_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE print_nlist_clmchfcg

END MODULE clmchfcg_scenario_mod
