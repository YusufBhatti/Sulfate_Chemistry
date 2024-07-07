! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: To initialize the nmspec array
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided by The University of Cambridge,
!  by, University of Leeds, University of Oxford, and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Called from addres
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: Fortran
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_nmspec_mod

USE ukca_tracer_stash,    ONLY: a_max_ukcavars, a_ukca_first, a_ukca_last

IMPLICIT NONE

PRIVATE
PUBLIC ::  nm_spec, ukca_set_nmspec, nm_spec_active, nmspec_len, &
           ukca_name2index

! Names for tracers (aka species)
! Once the  conflicting tracers from the RAQ chemistry are moved
! this should become a parameter array
INTEGER, PARAMETER :: nmspec_len=10
CHARACTER(LEN=nmspec_len), SAVE :: nm_spec(a_max_ukcavars)
! Names for active UKCA tracers - set at run time 
CHARACTER(LEN=nmspec_len), ALLOCATABLE, SAVE :: nm_spec_active(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_NMSPEC_MOD'

CONTAINS

SUBROUTINE ukca_set_nmspec()

USE ukca_option_mod,      ONLY: l_ukca_raq, l_ukca_raqaero, tr_ukca_a, jpctr
USE nlsizes_namelist_mod, ONLY: tr_ukca
USE UM_ParCore,           ONLY: mype
USE umPrintMgr,           ONLY: umPrint, umMessage, PrintStatus, PrStatus_Oper
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i, j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_NMSPEC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  This list must agree with the STASHmaster file.

! This is the generic list of tracers; values for RAQ and RAQ+aerosol chemistry
! are overwritten later.
! Tracers 98,99 & 100 are for lumped Nitrogen, Br and Cl for stratospheric
! chemistry, but can only be renamed in STASHmaster file, not in advt or
! nm_spec.
nm_spec(1:a_max_ukcavars) = (/                                     &
'O3        ','NO        ','NO3       ','NO2       ','N2O5      ',  &
'HO2NO2    ','HONO2     ','H2O2      ','CH4       ','CO        ',  & !10
'HCHO      ','MeOOH     ','HONO      ','C2H6      ','EtOOH     ',  &
'MeCHO     ','PAN       ','C3H8      ','n-PrOOH   ','i-PrOOH   ',  & !20
'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',  &
'O3_S      ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  & !30
'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',  &
'HCOOH     ','MeCO3H    ','MeCO2H    ','H2O       ','ISO2      ',  & !40
'Cl        ','ClO       ','Cl2O2     ','OClO      ','Br        ',  &
'BrO       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  & !50
'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',  &
'CF2Cl2    ','MeBr      ','N         ','O(3P)     ','MACRO2    ',  & !60
'MeCl      ','CF2ClBr   ','CCl4      ','CF2ClCFCl2','CHF2Cl    ',  &
'MeCCl3    ','CF3Br     ','H2OS      ','CH2Br2    ','H2        ',  & !70
'DMS       ','SO2       ','H2SO4     ','MSA       ','DMSO      ',  &
'NH3       ','CS2       ','COS       ','H2S       ','H         ',  & !80
'OH        ','HO2       ','MeOO      ','EtOO      ','MeCO3     ',  &
'n-PrOO    ','i-PrOO    ','EtCO3     ','MeCOCH2OO ','MeOH      ',  & !90
'Monoterp  ','Sec_Org   ','SESQUITERP','SO3       ','AROM      ',  &
'O(3P)_S   ','O(1D)_S   ','NO2       ','BrO       ','HCl       ',  & !100
'Nuc_SOL_ND','Nuc_SOL_SU','Ait_SOL_ND','Ait_SOL_SU','Ait_SOL_BC',  &
'Ait_SOL_OC','Acc_SOL_ND','Acc_SOL_SU','Acc_SOL_BC','Acc_SOL_OC',  & !110
'Acc_SOL_SS','Acc_SOL_DU','Cor_SOL_ND','Cor_SOL_SU','Cor_SOL_BC',  &
'Cor_SOL_OC','Cor_SOL_SS','Cor_SOL_DU','Ait_INS_ND','Ait_INS_BC',  & !120
'Ait_INS_OC','Acc_INS_ND','Acc_INS_DU','Cor_INS_ND','Cor_INS_DU',  &
'Nuc_SOL_OC','Ait_SOL_SS','Nuc_SOL_SO','Ait_SOL_SO','Acc_SOL_SO',  & !130
'Cor_SOL_SO','Nuc_SOL_NH','Ait_SOL_NH','Acc_SOL_NH','Cor_SOL_NH',  &
'Nuc_SOL_NT','Ait_SOL_NT','Acc_SOL_NT','Cor_SOL_NT','XXX       ',  & !140
'Anth_Prec ','Bio_Prec  ','Anth_Cond ','Bio_Cond  ','MSIA      ',  &
'XXX       ','XXX       ','XXX       ','PASSIVE O3','AGE OF AIR',  & !150
'RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ',  & 
'RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ',  & !160
'RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ',  & 
'RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ','RETIRED   ',  & !170
'RETIRED   ','RETIRED   ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !180
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !190
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !200
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !210
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !220
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !230
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !240
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & !250
'XXX       ','XXX       ','XXX       ','XXX       ','XXX       ',  & 
'XXX       '                                                       & !256
 /)

IF (l_ukca_raq .OR. l_ukca_raqaero) THEN
  ! Overwrite some tracers for RAQ and RAQ-AERO chemistry.
  ! If MODE aerosols are used with it but their positions
  ! change in the array then the list needs to be updated.

  ! Version above has underscore O3_S
  nm_spec(26) ='O3S       '

  nm_spec(39) ='MVK       '
  nm_spec(40) ='MVKOOH    '
  nm_spec(60) ='ORGNIT    '
  nm_spec(69) ='CH3OH     '
    
  nm_spec(90) ='RNC2H4    '
  nm_spec(93) ='C3H6      '
  nm_spec(94) ='C4H10     '
  nm_spec(95) ='s-BuOOH   '
  nm_spec(96) ='MEK       '
  nm_spec(97) ='TOLUENE   '
  nm_spec(98) ='MEMALD    '
  nm_spec(99) ='GLY       '
  nm_spec(100)='oXYLENE   '

! Two tracers are inconsistent between the two RAQ schemes because they
! appear in a different order in ukca_chem_raq and ukca_chem_raqaero
  IF (l_ukca_raq) THEN
    nm_spec(91) ='RNC3H6    '
    nm_spec(92) ='C2H4      '
  ELSE ! using l_ukca_raqaero
    nm_spec(81) ='RNC3H6    '
    nm_spec(82) ='C2H4      '
  END IF
END IF ! l_ukca_raq or l_ukca_raqaero

! Mode components: SU: sulphate, BC: black carbon, OC: organic carbon
!                  SS: sea-salt, Du: dust,         SO: organic carbon 2
!                  NH: ammonium, NT: nitrate,      ND: number density

! Make an array for the "real" nm_spec accounting for which 
! tracers are on. This allows us to map from the tracer array to 
! name of species and back
ALLOCATE (nm_spec_active(COUNT(tr_ukca_a)))

! loop over all entries in nm_spec, testing whether they are on or not
! If they are on, add their names into nm_spec_active
j = 0

! If high value of print status, output on PE0 only
IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE(umMessage,'(A)') 'NM spec active: '
    CALL umPrint(umMessage,src=RoutineName)
END IF

DO i = a_ukca_first, a_ukca_last
  IF (tr_ukca_a(i)) THEN
    j = j + 1
    nm_spec_active(j) = nm_spec(i)
    
    ! If high value of print status, output on PE0 only
    IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
        WRITE(umMessage,'(I4, 1X, A)') j, nm_spec_active(j)
        CALL umPrint(umMessage,src=RoutineName)
    END IF
  END IF
END DO 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_set_nmspec


FUNCTION ukca_name2index(tracer_name)

! This function is used in iau.F90 to update
! the Air Quality (AQ) Analysis increments.
! Find the index of a tracer based on the name given.
! if not found, stop with an ereport.
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: ukca_name2index

CHARACTER(LEN=*), INTENT(IN) :: tracer_name

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_NAME2INDEX'

INTEGER :: i
! ErrorStatus
INTEGER                    :: errcode
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check if nm_spec_active has been allocated, 
! if not, call ereport
IF (.NOT. ALLOCATED (nm_spec_active) ) THEN
  errcode = 1
  cmessage = 'Tried to look up contents of nm_spec_active but not allocated'
  CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)
END IF

! Loop over all of nm_spec_active
DO i = 1, SIZE(nm_spec_active)
  ! If names match
  IF (TRIM(ADJUSTL(tracer_name)) == TRIM(ADJUSTL(nm_spec_active(i)))) THEN
    ukca_name2index = i
    IF (lhook) CALL dr_hook(ModuleName//':'// &
                            RoutineName,zhook_out,zhook_handle)
    RETURN ! Exit the function
  END IF
END DO

! Call ereport as we haven't found the name if we get here
errcode = 2
cmessage = 'Failed to find ' // tracer_name // ' in nm_spec_active'
CALL ereport(ModuleName//':'//RoutineName, errcode, cmessage)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION ukca_name2index


END MODULE ukca_nmspec_mod
