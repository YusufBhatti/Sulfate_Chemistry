! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Converts dust emissions supplied in bins into modal number and mass
!    emissions for GLOMAP-mode.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and The Met Office.
!  See:  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90 
!
! ######################################################################

MODULE ukca_prim_du_mod

IMPLICIT NONE

 CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_PRIM_DU_MOD'

CONTAINS


 SUBROUTINE ukca_prim_du(row_length, rows, model_levels,           &
                         verbose, dust_flux,                       &
                         aer_num_primdu, aer_mas_primdu,           &
                         aer_mas_du_primdu)
                         
! Subroutine takes dust emission fluxes in kg m-2 s-1 in 6 size-bins and 
! converts them to mass and number fluxes for modes

! References:  Woodward et al JGR (2001), HCTN 87 
! ----------

! Inputs
! ------
! ROW_LENGTH   : columns 
! ROWS         : rows
! MODEL_LEVELS : levels
! VERBOSE      : controls amount of printed output  
! DUST_FLUX    : Dust emissions flux per bin kg(dust) m-2 s-1 

! Outputs
! -------
! AER_NUM_PRIMDU    : number emission flux to modes kg(air) m-2 s-1 
! AER_MAS_PRIMDU    : mass emission flux to modes kg(air) m-2 s-1 
! AER_MAS_DU_PRIMDU : mass emission flux to modes kg(dust) m-2 s-1 


! Input by module UKCA_CONSTANTS
! ---------------------------------
! MM_DA      : Molar mass of dry air
! AVC        : Avogadro's constant (molecules per mole)
! RA         : Dry air gas constant = 287.05 Jkg^-1 K^-1
! ZBOLTZ     : Stefan-Boltzmann constant (kg m2 s-2 K-1 molec-1)
! EMS_EPS    : Value of ems flux below which dont do ems (kg/box/s)
! DN_EPS     : Value of DELN below which do not carry out process
 
! Input by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES     : Number of modes set
! NCP        : Number of components set
! MODE       : Logical variable defining which modes are set.
! COMPONENT  : Logical variable defining which cpt are in which dsts
! RHOCOMP    : Component mass densities (kg/m3)
! MM         : Molar masses of components (kg/mol)
! CP_DU      : Index for dust component

! Input by module dust_parameters_mod
! --------------------------------------
! DMIN      : Minimum diameter of each bin
! DMAX      : Maximum diameter of each bin
! DREP      : Particle diameter for each bin
 
! Local Variables
! ---------------
! MS_PART_TO_KG : converts from molecules/particles to kg 
! DELMDNDCP     : Change in cpt mass conc. for particular emission
 
!--------------------------------------------------------------------------

 USE conversions_mod,        ONLY: pi
 USE ukca_constants,         ONLY: avc, ra, zboltz,                &
                                   ems_eps, dn_eps
 USE ukca_mode_setup,        ONLY: nmodes, ncp, mode, component,   &
                                   cp_du, mm, rhocomp
 USE dust_parameters_mod,    ONLY: dmin, dmax ,drep, ndiv
 USE ereport_mod,            ONLY: ereport
 USE yomhook,                ONLY: lhook, dr_hook
 USE parkind1,               ONLY: jprb, jpim
 USE ereport_mod,            ONLY: ereport
 USE umPrintMgr
 USE errormessagelength_mod, ONLY: errormessagelength

 IMPLICIT NONE

! .. Subroutine interface
 INTEGER, INTENT(IN) :: row_length              ! No of columns
 INTEGER, INTENT(IN) :: rows                    ! No of rows
 INTEGER, INTENT(IN) :: model_levels            ! No of model levels
 INTEGER, INTENT(IN) :: verbose                 ! Sets level of output

 REAL, INTENT(IN)  :: dust_flux(row_length,rows,ndiv) ! kg(dust) m-2 s-1
 
! modal mass and number emission fluxes 
!       
!     aer_num_primdu contains number emission flux in kg(air) m-2 s-1
 REAL, INTENT(INOUT) :: aer_num_primdu(row_length,rows,model_levels,nmodes) 
!     aer_mass_primdu contains mass emission flux in kg(air) m-2 s-1
 REAL, INTENT(INOUT) :: aer_mas_primdu(row_length,rows,model_levels,nmodes) 
!     aer_mass_du_primdu contains mass emission flux in kg(dust) m-2 s-1
 REAL, INTENT(INOUT)::aer_mas_du_primdu(row_length,rows,model_levels,nmodes) 

! Local variables
 INTEGER :: imode                  ! loop counter for modes
 INTEGER :: idiv                   ! loop counter for dust bins

 REAL, PARAMETER :: nm3conv = 1e27 ! converts m^3 to nm^3 
 
 REAL    :: ms_part_to_kg          ! converts molecules(air)/particle to kg
 REAL    :: emcdu_emit(row_length,rows) ! component mass flux (kg-cpt/m2/s)
 REAL    :: delmdndcp(row_length,rows)  ! component mass flux (kg-air/m2/s)
 REAL    :: deln(row_length,rows)       ! component number flux 
 REAL    :: totnum(row_length,rows) ! tot emitted particle num (ptcls/m2/s)
 REAL    :: emcvol(row_length,rows) ! tot emitted particle vol (nm3/m2/s)
 REAL    :: conv_dneps(row_length,rows) ! converts partcls/cc/s to kg/m2/s
 REAL    :: mm_da                       ! Molar mass of dry air (kg/mol)

 REAL    :: fracduem(ndiv,nmodes)   ! fraction dust emitted to each mode

 LOGICAL, SAVE :: first = .TRUE.    ! T for first call

 INTEGER           :: errcode       ! Error code
 CHARACTER(LEN=errormessagelength) :: cmessage

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PRIM_DU'



 IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Molar mass of dry air (kg/mol)
 mm_da = avc*zboltz/ra

! Check that ndiv is 6
 IF (ndiv /= 6) THEN
   cmessage = ' No of dust bins incorrect, should be set to six'
   errcode = 1
   WRITE(umMessage,'(A72,A6,I10)') cmessage,' ndiv=',ndiv
   CALL umPrint(umMessage,src='ukca_prim_du')
   CALL EREPORT('UKCA_PRIM_DU',errcode,cmessage)
 END IF


! Mapping of dust bins to insoluble modes. Currently only bins 2-5 used.
! NB This mapping is not ideal, and will be changed in future versions.

 fracduem(1,:)=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
 fracduem(2,:)=(/0.0,0.0,0.0,0.0,0.0,1.0,0.0/)
 fracduem(3,:)=(/0.0,0.0,0.0,0.0,0.0,0.5,0.5/)
 fracduem(4,:)=(/0.0,0.0,0.0,0.0,0.0,0.0,1.0/)
 fracduem(5,:)=(/0.0,0.0,0.0,0.0,0.0,0.0,1.0/)
 fracduem(6,:)=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
!

! convert from molecules_air/particles to kg(air)
 ms_part_to_kg = mm_da/avc   


! for each bin convert emissions (kg m-2 s-1) to mass and number flux 
! per mode in appropriate units 

 DO idiv = 1,ndiv    ! loop over bins 
             
   IF (verbose > 2) THEN
     WRITE(umMessage,'(A20,I6,3E15.4)') 'dust_flux: idiv: ',idiv,  &
                            MINVAL(dust_flux(:,:,idiv)),           &
                            MAXVAL(dust_flux(:,:,idiv)),           &
                            SUM(dust_flux(:,:,idiv))/              &
                            SIZE(dust_flux(:,:,idiv))
    CALL umPrint(umMessage,src='ukca_prim_du')
   END IF
   
   DO imode = 1,nmodes                    ! loop over UKCA modes
     IF (mode(imode)) THEN
       IF (fracduem(idiv,imode) > 0.0 .AND.                        & 
           ANY(dust_flux(:,:,idiv) > 1.0E-20)) THEN
         emcdu_emit(:,:) = dust_flux(:,:,idiv)*                    &
                           fracduem(idiv,imode)       ! kg-cpt/m2/s

!             dust volume emission flux (m3  per m2 per s)

         emcvol(:,:) = (emcdu_emit(:,:)/rhocomp(cp_du))

!             dust number emission flux in kg(air) m-2 s-1
!             i.e. mass flux of air with same number of molecules as there
!             are particles of dust in the emission flux

         totnum(:,:) = emcvol(:,:) * 6./(pi*(drep(idiv)**3.0))
         deln(:,:) = ms_part_to_kg*totnum(:,:)

         aer_num_primdu(:,:,1,imode) =                             &
            aer_num_primdu(:,:,1,imode) + deln(:,:) 

!             dust mass emission flux in kg(air) m-2 s-1
!             i.e. mass flux of air with same number of mols as there are mols
!             of dust in the emission flux

         delmdndcp(:,:) = emcdu_emit(:,:) * mm_da/mm(cp_du)

         aer_mas_primdu(:,:,1,imode) =                             &
             aer_mas_primdu(:,:,1,imode) + delmdndcp(:,:) 

!             dust mass emission flux in kg(dust) m-2 s-1 
!             for new emissions code
         aer_mas_du_primdu(:,:,1,imode) =                          &
             aer_mas_du_primdu(:,:,1,imode) +emcdu_emit(:,:)
! 

       END IF   ! fracduem(idiv,imode) > 0.0
     END IF   ! mode(imode)
   END DO   ! nmodes
 END DO   ! ndiv
   
 IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
 
 RETURN
 END SUBROUTINE UKCA_PRIM_DU
 
END MODULE ukca_prim_du_mod      
