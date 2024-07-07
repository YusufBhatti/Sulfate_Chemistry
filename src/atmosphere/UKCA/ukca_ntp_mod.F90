! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to hold derived type for non-transported prognostic variables
!  and subroutines/functions related to their use.
!
!  Non-transported prognostics are variables which are part of the 
!  model state but are not transported. For example short lived radicals
!  such as O(3P) where the value from the previous timestep is used 
!  as an initial value to improve the solution accuracy.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds,
!  University of Oxford and The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_ntp_mod

! Standard UM modules used 
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: umPrint, umMessage

IMPLICIT NONE

! All variables and subroutines are private by default
PRIVATE

! Dr hook variables/parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

! subroutines/functions which are public
PUBLIC ntp_init, name2ntpindex, stash2ntpindex, print_all_ntp, ntp_dealloc

! The size of the all_ntp array is defined here.
! If adding or removing entries remember to change
! the size of dim_ntp 
INTEGER, PARAMETER, PUBLIC :: dim_ntp = 73

! Derived type used to hold all information for each NTP. 
! Section, item, data_3d, l_required, name
TYPE, PUBLIC :: ntp_type
  INTEGER           :: section
  INTEGER           :: item
  REAL,ALLOCATABLE  :: data_3d(:,:,:)   ! only 3D data allowed
  LOGICAL           :: l_required       ! required by UKCA
  CHARACTER(LEN=20) :: varname
END TYPE ntp_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_NTP_MOD'

CONTAINS

! ----------------------------------------------------------------------
SUBROUTINE ntp_init(all_ntp)
! ----------------------------------------------------------------------
! Description:
! initialise values of data in all_ntp structure
! 
! Method:
! Section, item and name have to be set manually here via a call 
! add_ntp_item. 
! Whether a specific entry is required from D1 is calculated at runtime
! and set in add_ntp_item after the other values are set
! Must have exactly the right size - model will stop with an error 
! message if it is too small or too large.
!
! To add an additional non-transported prognostic: 
! 1. increment dim_ntp above
! 2. add a call to add_ntp_item in this subroutine to add it to the array
! 3. unless it is a chemical species in an explicit BE scheme, add
!    the logic to the function ntp_req to control when it is on or off.
!
! ----------------------------------------------------------------------

USE ukca_d1_defs,  ONLY: ukca_sect
USE ukca_cdnc_mod, ONLY: stashc_cdnc, stashc_cdnc3
USE ukca_setup_chem_mod, ONLY: i_ho2_be, i_oh_be
USE um_stashcode_mod, ONLY:      stashcode_dryd_ait_sol,         &
    stashcode_dryd_acc_sol,      stashcode_dryd_cor_sol,         &
    stashcode_dryd_ait_insol,    stashcode_dryd_acc_insol,       &
    stashcode_dryd_cor_insol,    stashcode_wetd_ait_sol,         &
    stashcode_wetd_acc_sol,      stashcode_wetd_cor_sol,         &
    stashcode_rho_ait_sol,       stashcode_rho_acc_sol,          &
    stashcode_rho_cor_sol,       stashcode_rho_ait_insol,        &
    stashcode_rho_acc_insol,     stashcode_rho_cor_insol,        &
    stashcode_pvol_ait_su_sol,   stashcode_pvol_ait_bc_sol,      &
    stashcode_pvol_ait_oc_sol,   stashcode_pvol_ait_so_sol,      &
    stashcode_pvol_ait_no3_sol,  stashcode_pvol_ait_h2o_sol,     &
    stashcode_pvol_acc_su_sol,   stashcode_pvol_acc_bc_sol,      &
    stashcode_pvol_acc_oc_sol,   stashcode_pvol_acc_so_sol,      &
    stashcode_pvol_acc_du_sol,   stashcode_pvol_acc_ss_sol,      &
    stashcode_pvol_acc_no3_sol,  stashcode_pvol_acc_h2o_sol,     &
    stashcode_pvol_cor_su_sol,   stashcode_pvol_cor_bc_sol,      &
    stashcode_pvol_cor_oc_sol,   stashcode_pvol_cor_so_sol,      &
    stashcode_pvol_cor_du_sol,   stashcode_pvol_cor_ss_sol,      &
    stashcode_pvol_cor_no3_sol,  stashcode_pvol_cor_h2o_sol,     &
    stashcode_pvol_ait_bc_insol, stashcode_pvol_ait_oc_insol,    &
    stashcode_pvol_acc_du_insol, stashcode_pvol_cor_du_insol

IMPLICIT NONE

TYPE(ntp_type), INTENT(INOUT) :: all_ntp(dim_ntp)

! Local variables
INTEGER       :: msect                                     ! 1000 * section
INTEGER       :: errcode                                   ! error code
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_INIT'
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

msect = 1000*ukca_sect

! Set all item and section numbers to -999 to do an error check at the end
all_ntp(:)%section    = -999
all_ntp(:)%item       = -999
all_ntp(:)%varname    = 'Not set!!!!!!!!!!!!!'
all_ntp(:)%l_required = .FALSE.

! Initialise metadata for all entries. Where these are chemical compounds
! which are not transported, but are stored in the dump, these names must
! be the same as those in the chch_defs arrays as these are used to check
! if a compound is on and where to put the data obtained from D1

! The order of these is NOT important, but adding them in stashcode order
! makes it easier to read

! Aerosol surface area
CALL add_ntp_item(all_ntp, section=ukca_sect, item=966,          &
  varname='surfarea  ')

! Cube root of cloud droplet number concentration
CALL add_ntp_item(all_ntp, section=ukca_sect, item=stashc_cdnc3, &
  varname='cdnc3     ')

! Cloud droplet number concentration 
CALL add_ntp_item(all_ntp, section=ukca_sect, item=stashc_cdnc,  &
  varname='cdnc      ')

! Stratospheric HO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=969,          &
  varname='HO2S      ')

! Stratospheric OH
CALL add_ntp_item(all_ntp, section=ukca_sect, item=970,          &
  varname='OHS       ')

! Stratospheric O(1D)
CALL add_ntp_item(all_ntp, section=ukca_sect, item=971,          &
  varname='O(1D)S    ')

! Stratospheric O(3P)
CALL add_ntp_item(all_ntp, section=ukca_sect, item=972,          &
  varname='O(3P)S    ')

! Heterogeneous self reaction rate of HO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=973,          &
  varname='het_ho2   ')

! Heterogeneous loss rate of N2O5
CALL add_ntp_item(all_ntp, section=ukca_sect, item=974,          &
  varname='het_n2o5  ')

! TOLP1
CALL add_ntp_item(all_ntp, section=ukca_sect, item=975,          &
  varname='TOLP1     ')

! HOIPO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=976,          &
  varname='HOIPO2    ')

! HOMVKO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=977,          &
  varname='HOMVKO2   ')

! MEMALD1
CALL add_ntp_item(all_ntp, section=ukca_sect, item=978,          &
  varname='MEMALD1   ')

! OXYL1
CALL add_ntp_item(all_ntp, section=ukca_sect, item=979,          &
  varname='OXYL1     ')

! HOC3H6O2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=980,          &
  varname='HOC3H6O2  ')

! HOC2H4O
CALL add_ntp_item(all_ntp, section=ukca_sect, item=981,          &
  varname='HOC2H4O2  ')

! MEKO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=982,          &
  varname='MEKO2     ')

! MeCOCH2OO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=983,          &
  varname='MeCOCH2OO ')

! MeCOCH2OO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=984,          &
  varname='MeCOCH2OO ')

! EtCO3
CALL add_ntp_item(all_ntp, section=ukca_sect, item=985,          &
  varname='EtCO3     ')

! i-PrOO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=986,          &
  varname='i-PrOO    ')

! s-BuOO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=987,          &
  varname='s-BuOO    ')

! n-PrOO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=988,          &
  varname='n-PrOO    ')

! MeCO3
CALL add_ntp_item(all_ntp, section=ukca_sect, item=989,          &
  varname='MeCO3     ')

! EtOO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=990,          &
  varname='EtOO      ')

! MeOO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=991,          &
  varname='MeOO      ')

! HCl
CALL add_ntp_item(all_ntp, section=ukca_sect, item=992,          &
  varname='HCl       ')

! HO2. Note that for coupling to CLASSIC, the item numbers of OH/HO2 must be
! consistent with ukca_setup_chem_mod and so we use parameters from that MODULE
CALL add_ntp_item(all_ntp, section=ukca_sect, item=i_ho2_be, &
  varname='HO2       ')

! BrO
CALL add_ntp_item(all_ntp, section=ukca_sect, item=994,          &
  varname='BrO       ')

! OH. 
CALL add_ntp_item(all_ntp, section=ukca_sect, item=i_oh_be, &
  varname='OH        ')

! NO2
CALL add_ntp_item(all_ntp, section=ukca_sect, item=996,          &
  varname='NO2       ')

! O(1D)
CALL add_ntp_item(all_ntp, section=ukca_sect, item=997,          &
  varname='O(1D)     ')

! O(3P)
CALL add_ntp_item(all_ntp, section=ukca_sect, item=998,          &
  varname='O(3P)     ')

! items for use by RADAER:
! Dry diameter
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_ait_sol - msect, varname='drydiam_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_acc_sol - msect, varname='drydiam_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_cor_sol - msect, varname='drydiam_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_ait_insol - msect, varname='drydiam_ait_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_acc_insol - msect, varname='drydiam_acc_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_dryd_cor_insol - msect, varname='drydiam_cor_insol   ')

! Wet diameter
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_wetd_ait_sol - msect, varname='wetdiam_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_wetd_acc_sol - msect, varname='wetdiam_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_wetd_cor_sol - msect, varname='wetdiam_cor_sol     ')

! Aerosol density
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_ait_sol - msect, varname='aerdens_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_acc_sol - msect, varname='aerdens_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_cor_sol - msect, varname='aerdens_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_ait_insol - msect, varname='aerdens_ait_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_acc_insol - msect, varname='aerdens_acc_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_rho_cor_insol - msect, varname='aerdens_cor_insol   ')

! Partial volume
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_su_sol - msect, varname='pvol_su_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_bc_sol - msect, varname='pvol_bc_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_oc_sol - msect, varname='pvol_oc_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_so_sol - msect, varname='pvol_so_ait_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_h2o_sol - msect, varname='pvol_h2o_ait_sol    ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_su_sol - msect, varname='pvol_su_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_bc_sol - msect, varname='pvol_bc_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_oc_sol - msect, varname='pvol_oc_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_ss_sol - msect, varname='pvol_ss_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_no3_sol - msect, varname='pvol_no3_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_du_sol - msect, varname='pvol_du_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_so_sol - msect, varname='pvol_so_acc_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_h2o_sol - msect, varname='pvol_h2o_acc_sol    ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_su_sol - msect, varname='pvol_su_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_bc_sol - msect, varname='pvol_bc_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_oc_sol - msect, varname='pvol_oc_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_ss_sol - msect, varname='pvol_ss_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_no3_sol - msect, varname='pvol_no3_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_du_sol - msect, varname='pvol_du_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_so_sol - msect, varname='pvol_so_cor_sol     ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_h2o_sol - msect, varname='pvol_h2o_cor_sol    ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_bc_insol - msect, varname='pvol_bc_ait_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_ait_oc_insol - msect, varname='pvol_oc_ait_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_acc_du_insol - msect, varname='pvol_du_acc_insol   ')
CALL add_ntp_item(all_ntp, section=ukca_sect,                                  &
  item=stashcode_pvol_cor_du_insol - msect, varname='pvol_du_cor_insol   ')

! Finally, check metadata for all entries is set
IF (ANY(all_ntp(:)%section == -999) .OR.                &
    ANY(all_ntp(:)%item    == -999) .OR.                &
    ANY(all_ntp(:)%varname == 'Not set!!!!!!!!!!!!!')) THEN 

! If some values are not set, write some useful messages and stop the model
  WRITE(umMessage,'(A)') 'one or more entries in all_ntp not set'
  CALL umPrint(umMessage,src=RoutineName)
  CALL print_all_ntp(all_ntp)
  WRITE(cmessage,'(A)') 'one or more entries in all_ntp not set'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_init

! ----------------------------------------------------------------------
SUBROUTINE add_ntp_item(all_ntp, section, item, varname)
! ----------------------------------------------------------------------
! Description:
! Add data to the all_ntp array. If try and add too many stop the model. 
! 1. Increment the internal counter ntp_index
! 2. Check if the counter is within the the array size
! 3. Add the values passed in 
! 4. Set the logical flag for whether this is required by calling 
!    the function ntp_req
! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(ntp_type), INTENT(INOUT) :: all_ntp(dim_ntp)
INTEGER, INTENT(IN) :: section
INTEGER, INTENT(IN) :: item
CHARACTER(LEN=*), INTENT(IN) :: varname

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ADD_NTP_ITEM'
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

! Counter for position in all_ntp array. Initialised to zero
! and incremented by one on each call to this routine.
INTEGER, SAVE :: ntp_index = 0
INTEGER       :: errcode                                   ! error code
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! increment the location we store information to
ntp_index = ntp_index + 1

! Test that the all_ntp array is big enough - if not, abort.
IF (ntp_index > dim_ntp) THEN
  WRITE(cmessage,'(A,I6)') 'dim_ntp too small in ukca_ntp_mod:', dim_ntp
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF 

! now set all values for this entry

! section, item and variable name are all input arguments to this 
! subroutine
all_ntp(ntp_index)%section = section
all_ntp(ntp_index)%item = item
all_ntp(ntp_index)%varname = varname

! use the logical function ntp_req to test whether this
! is required for the current model run
all_ntp(ntp_index)%l_required = ntp_req(varname)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE add_ntp_item

! ----------------------------------------------------------------------
LOGICAL FUNCTION ntp_req(varname)
! ----------------------------------------------------------------------
! Description:
! Use run time logic to test if varname is on in any specific model run. 
! For chemistry, test nadvt. Need to treat others as special cases. 
!
! ----------------------------------------------------------------------
USE asad_mod, ONLY: nadvt, O1D_in_ss, O3P_in_ss 
USE ukca_option_mod,  ONLY: l_ukca_trophet, l_ukca_arg_act, l_ukca_aie1,      &
  l_ukca_aie2, l_ukca_mode, l_ukca_strat, l_ukca_strattrop, l_ukca_stratcfc,  &
  ukca_int_method, l_ukca_radaer
USE ukca_chem_schemes_mod, ONLY: int_method_be_explicit
USE ukca_mode_setup,       ONLY: mode, component, mode_ait_sol, mode_acc_sol,  &
                                 mode_cor_sol, mode_ait_insol, mode_acc_insol, &
                                 mode_cor_insol, cp_su, cp_bc, cp_oc, cp_cl,   &
                                 cp_no3, cp_du, cp_so

IMPLICIT NONE


CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_REQ'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ntp_req = .FALSE.

! First handle the GLOMAP variables. To allow for the use of allocatable
! arrays for mode and component, handle these together

IF ( varname (1:7) == 'drydiam' .OR. &
     varname (1:7) == 'wetdiam' .OR. &
     varname (1:7) == 'aerdens' .OR. &
     varname (1:4) == 'pvol' ) THEN
     IF (l_ukca_mode) THEN
        SELECT CASE (varname)
        ! RADAER items, check which modes and components are active
        ! Dry diameter
        CASE ('drydiam_ait_sol     ')           ! Aitken-sol dry diameter
          ntp_req = mode(mode_ait_sol) .AND. l_ukca_radaer
        CASE ('drydiam_acc_sol     ')           ! accumulation-sol dry diameter
          ntp_req = mode(mode_acc_sol) .AND. l_ukca_radaer
        CASE ('drydiam_cor_sol     ')           ! coarse-sol dry diameter
          ntp_req = mode(mode_cor_sol) .AND. l_ukca_radaer
        CASE ('drydiam_ait_insol   ')           ! Aitken-sol dry diameter
          ntp_req = mode(mode_ait_insol) .AND. l_ukca_radaer
        CASE ('drydiam_acc_insol   ')           ! accumulation-sol dry diameter
          ntp_req = mode(mode_acc_insol) .AND. l_ukca_radaer
        CASE ('drydiam_cor_insol   ')           ! coarse-sol dry diameter
          ntp_req = mode(mode_cor_insol) .AND. l_ukca_radaer

        ! Wet diameter
        CASE ('wetdiam_ait_sol     ')           ! Aitken-sol wet diameter
          ntp_req = mode(mode_ait_sol) .AND. l_ukca_radaer
        CASE ('wetdiam_acc_sol     ')           ! accumulation-sol wet diameter
          ntp_req = mode(mode_acc_sol) .AND. l_ukca_radaer
        CASE ('wetdiam_cor_sol     ')           ! coarse-sol wet diameter
          ntp_req = mode(mode_cor_sol) .AND. l_ukca_radaer

        ! Aerosol density
        CASE ('aerdens_ait_sol     ')           ! Aitken-sol aerosol density
          ntp_req = mode(mode_ait_sol) .AND. l_ukca_radaer
        CASE ('aerdens_acc_sol     ')           ! accumulation-sol " density
          ntp_req = mode(mode_acc_sol) .AND. l_ukca_radaer
        CASE ('aerdens_cor_sol     ')           ! coarse-sol aerosol density
          ntp_req = mode(mode_cor_sol) .AND. l_ukca_radaer
        CASE ('aerdens_ait_insol   ')           ! Aitken-sol aerosol density
          ntp_req = mode(mode_ait_insol) .AND. l_ukca_radaer
        CASE ('aerdens_acc_insol   ')           ! accumulation-sol " density
          ntp_req = mode(mode_acc_insol) .AND. l_ukca_radaer
        CASE ('aerdens_cor_insol   ')           ! coarse-sol aerosol density
          ntp_req = mode(mode_cor_insol) .AND. l_ukca_radaer

  ! Partial volume
        CASE ('pvol_su_ait_sol     ')           ! Aitken-sol sulphate
          ntp_req = component(mode_ait_sol,cp_su) .AND. l_ukca_radaer
        CASE ('pvol_bc_ait_sol     ')           ! Aitken-sol black carbon
          ntp_req = component(mode_ait_sol,cp_bc) .AND. l_ukca_radaer
        CASE ('pvol_oc_ait_sol     ')           ! Aitken-sol organic carbon
          ntp_req = component(mode_ait_sol,cp_oc) .AND. l_ukca_radaer
        CASE ('pvol_so_ait_sol     ')           ! Aitken-sol secondary organic
          ntp_req = component(mode_ait_sol,cp_so) .AND. l_ukca_radaer
        CASE ('pvol_h2o_ait_sol    ')           ! Aitken-sol H2O
          ntp_req = mode(mode_ait_sol) .AND. l_ukca_radaer
        CASE ('pvol_su_acc_sol     ')           ! accumulation-sol sulphate
          ntp_req = component(mode_acc_sol,cp_su) .AND. l_ukca_radaer
        CASE ('pvol_bc_acc_sol     ')           ! accumulation-sol black carbon
          ntp_req = component(mode_acc_sol,cp_bc) .AND. l_ukca_radaer
        CASE ('pvol_oc_acc_sol     ')           ! accumulatn-sol organic carbon
          ntp_req = component(mode_acc_sol,cp_oc) .AND. l_ukca_radaer
        CASE ('pvol_ss_acc_sol     ')           ! accumulation-sol sea-salt
          ntp_req = component(mode_acc_sol,cp_cl) .AND. l_ukca_radaer
        CASE ('pvol_no3_acc_sol     ')         ! accumulation-sol nitrate
          IF (UBOUND(component,DIM=2) >= cp_no3) &
          ntp_req = component(mode_acc_sol,cp_no3) .AND. l_ukca_radaer
        CASE ('pvol_du_acc_sol     ')           ! accumulation-sol dust
          ntp_req = component(mode_acc_sol,cp_du) .AND. l_ukca_radaer
        CASE ('pvol_so_acc_sol     ')           ! accumulation-sol 2ndy organic
          ntp_req = component(mode_acc_sol,cp_so) .AND. l_ukca_radaer
        CASE ('pvol_h2o_acc_sol    ')           ! accumulation-sol H2O
          ntp_req = mode(mode_acc_sol) .AND. l_ukca_radaer
        CASE ('pvol_su_cor_sol     ')           ! coarse-sol sulphate
          ntp_req = component(mode_cor_sol,cp_su) .AND. l_ukca_radaer
        CASE ('pvol_bc_cor_sol     ')           ! coarse-sol black carbon
          ntp_req = component(mode_cor_sol,cp_bc) .AND. l_ukca_radaer
        CASE ('pvol_oc_cor_sol     ')           ! coarse-sol organic carbon
          ntp_req = component(mode_cor_sol,cp_oc) .AND. l_ukca_radaer
        CASE ('pvol_ss_cor_sol     ')           ! coarse-sol sea-salt
          ntp_req = component(mode_cor_sol,cp_cl) .AND. l_ukca_radaer
        CASE ('pvol_no3_cor_sol     ')           ! coarse-sol nitrate
          IF (UBOUND(component,DIM=2) >= cp_no3) &
          ntp_req = component(mode_cor_sol,cp_no3) .AND. l_ukca_radaer
        CASE ('pvol_du_cor_sol     ')           ! coarse-sol dust
          ntp_req = component(mode_cor_sol,cp_du) .AND. l_ukca_radaer
        CASE ('pvol_so_cor_sol     ')           ! coarse-sol 2ndy organic
          ntp_req = component(mode_cor_sol,cp_so) .AND. l_ukca_radaer
        CASE ('pvol_h2o_cor_sol    ')           ! coarse-sol H2O
          ntp_req = mode(mode_cor_sol) .AND. l_ukca_radaer
        CASE ('pvol_bc_ait_insol   ')           ! coarse-insol black carbon
          ntp_req = component(mode_ait_insol,cp_bc) .AND. l_ukca_radaer
        CASE ('pvol_oc_ait_insol   ')           ! coarse-insol organic carbon
          ntp_req = component(mode_ait_insol,cp_oc) .AND. l_ukca_radaer
        CASE ('pvol_du_acc_insol   ')           ! accumulation-insol dust
          ntp_req = component(mode_acc_insol,cp_du) .AND. l_ukca_radaer
        CASE ('pvol_du_cor_insol   ')           ! coarse-insol dust
          ntp_req = component(mode_cor_insol,cp_du) .AND. l_ukca_radaer
        CASE DEFAULT
          ntp_req = .FALSE.
        END SELECT
     ELSE
       ntp_req = .FALSE.
     END IF

ELSE
    SELECT CASE (varname)

    ! Aerosol surface area stored in dump if using GLOMAP-mode.
    ! Only needed as prognostic under certain circumstances,
    ! however, we want to make it available as a diagnostic whenever
    ! GLOMAP-mode is being used.
    CASE ('surfarea  ')
        ntp_req = l_ukca_mode

    ! CDNC prognostics - on if aerosol indirect effects are on
    ! or ACTIVATE on
    CASE ('cdnc     ')
        ntp_req = (l_ukca_aie1 .OR. l_ukca_aie2 .OR. l_ukca_arg_act)

    ! cube root of CDNC only comes from activate
    CASE ('cdnc3    ')
        ntp_req = l_ukca_arg_act

    ! The heterogeneous loss rates are on if l_ukca_trophet
    ! is true
    CASE ('het_ho2   ','het_n2o5  ')
      ntp_req = l_ukca_trophet

    ! O1D and O3P - special case, used by both NR and BE schemes
    ! and on if they are in SS (some NR schemes have O3P as a tracer
    ! in which case it is stored in the tracer array not in the NTP
    ! structure)
    CASE ('O(1D)     ')
      ntp_req = O1D_in_ss

    CASE ('O(3P)     ')
      ntp_req = O3P_in_ss

    ! Lumped species - on for chemical schemes using lumping
    CASE ('NO2       ','BrO        ','HCl        ')
      ntp_req = l_ukca_strat .OR. l_ukca_strattrop .OR. l_ukca_stratcfc

    ! All others are checked against whether they are in
    ! the nadvt array (only for BE explicit schemes)
    CASE DEFAULT
      ntp_req = ANY(nadvt(:) == varname(1:10)) .AND.                       &
               (ukca_int_method  == int_method_be_explicit)

    END SELECT
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END FUNCTION ntp_req

! ----------------------------------------------------------------------
INTEGER FUNCTION name2ntpindex(all_ntp, varname)
! ----------------------------------------------------------------------
! Description:
! Given a variable name look up where it is in the all_ntp array. 
! Abort here if failed to look up.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Input values are the all_ntp array and the
! variable name to look up
TYPE(ntp_type),   INTENT(IN) :: all_ntp(dim_ntp)
CHARACTER(LEN=*), INTENT(IN) :: varname

! Local variables
INTEGER :: i
INTEGER :: errcode                                   ! error code
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NAME2NTPINDEX'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set a default value. If it is still this at the end, the search has failed.
name2ntpindex = -999

! Search all entries in ntp array to find the one we want
DO i = 1, dim_ntp
  IF (all_ntp(i)%varname == varname) THEN
    name2ntpindex = i  
    EXIT
  END IF
END DO

! If name2ntpindex is -999 then call ereport
IF (name2ntpindex == -999) THEN
  WRITE(cmessage,'(A,A,A)') 'Failed to find: ', varname, ' in NTP structure.'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END FUNCTION name2ntpindex
! ----------------------------------------------------------------------

INTEGER FUNCTION stash2ntpindex(all_ntp, section, item)
! ----------------------------------------------------------------------
! Description:
! Given a stash code look up where it is in the all_ntp array.
! Abort here if failed to look up.
! ----------------------------------------------------------------------

IMPLICIT NONE

! Input values are the all_ntp array and the
! section and item number to look up
TYPE(ntp_type), INTENT(IN) :: all_ntp(dim_ntp)
INTEGER, INTENT(IN)        :: section, item

! Local variables
INTEGER :: i
INTEGER       :: errcode                                   ! error code
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'STASH2NTPINDEX'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! set a default value. If it is still this on exit, the search has failed.
stash2ntpindex = -999

! Search all entries in ntp array to find the one we want
DO i = 1, dim_ntp
  IF (all_ntp(i)%section == section .AND. all_ntp(i)%item == item) THEN
    stash2ntpindex = i  
    EXIT
  END IF
END DO

! If name2ntpindex is -999 then call ereport
IF (stash2ntpindex == -999) THEN
  WRITE(cmessage,'(A,2I6,A)') 'Failed to find: ', section, item,           &
    ' in NTP structure.'
  errcode = 1
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END FUNCTION stash2ntpindex

! ----------------------------------------------------------------------
SUBROUTINE print_all_ntp(all_ntp)
! ----------------------------------------------------------------------
! Description:
! Print all information in the all_ntp array. Used for debugging.
!
! Method:
! Loop through each entry and print values
! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(ntp_type), INTENT(IN) :: all_ntp(dim_ntp)

INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'PRINT_ALL_NTP'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Print header
WRITE(umMessage,'(A)') '  i section   item l_required varname'
CALL umPrint(umMessage,src=RoutineName)

! Loop over all items in ntp, printing metadata and min/max if allocated
DO i =1, dim_ntp
  WRITE(umMessage,'(I4,2X,I6,1X,I6,1X,L1,10X,A)') i , all_ntp(i)%section, &
    all_ntp(i)%item, all_ntp(i)%l_required, all_ntp(i)%varname
  CALL umPrint(umMessage,src=RoutineName)
  IF (ALLOCATED(all_ntp(i)%data_3d)) THEN
    WRITE(umMessage,'(A,2E12.4)') 'Min/max: ',                            &
      MINVAL(all_ntp(i)%data_3d),MAXVAL(all_ntp(i)%data_3d)
    CALL umPrint(umMessage,src=routinename)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_all_ntp

! ----------------------------------------------------------------------
SUBROUTINE ntp_dealloc(all_ntp)
! ----------------------------------------------------------------------
! Description:
! Deallocate all data arrays in all_ntp array.
! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(ntp_type), INTENT(INOUT) :: all_ntp(dim_ntp)

INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'NTP_DEALLOC'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1, SIZE(all_ntp)
  ! If ntp allocated, we need to deallocate
  IF (ALLOCATED(all_ntp(i)%data_3d)) DEALLOCATE(all_ntp(i)%data_3d)
END DO 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ntp_dealloc
! ----------------------------------------------------------------------
END MODULE ukca_ntp_mod

