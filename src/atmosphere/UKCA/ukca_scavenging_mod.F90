! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module for plume scavenging, contains definitions, and contained routines:
!    ukca_plume_scav   (called from convec2)
!    ukca_calc_aqueous (called from ukca_plume_scav)
!    ukca_set_conv_indices (called from ukca_plume_scav)
!
!  Method:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office. See
!  www.ukca.ac.uk
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran
!
! ######################################################################
MODULE ukca_scavenging_mod

USE ukca_d1_defs, ONLY: code
USE ukca_mode_setup,           ONLY: nmodes, ncp
USE cv_run_mod,   ONLY: i_convection_vn, i_convection_vn_6a,            &
                        i_convection_vn_5a
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE
PRIVATE
SAVE

REAL, PARAMETER :: precip_lim = 1e-50   ! No calculation if ppn less than this
                                        ! also used to to test mass flux

! Information about how tracers are scavenged
TYPE tr_scav_struct
  ! Index of first UKCA tracer in convection tracer array
  INTEGER :: i_ukca_first
  ! Index of last UKCA  tracer in convection tracer array
  INTEGER :: i_ukca_last
END TYPE tr_scav_struct

! structure holding information for scavenging of all tracers
TYPE (tr_scav_struct), PUBLIC :: tracer_info

! Index of tracers for mode number, size is the number of modes.
! Holds the tracer number in the tracer_ukca array. 
! Defaults to zero if the mode is not found.
INTEGER :: nmr_index_um(nmodes)

! Index of tracers for mass mixing ratio of component, size is
! no. of modes * no. of components. 
! Holds the appropriate tracer number in the tracer_ukca array. 
! Defaults to -1 if the component is not found.
INTEGER, PUBLIC, ALLOCATABLE :: mmr_index_um(:,:)

REAL    :: rscavn(nmodes)         ! Number scavenging ratios
REAL    :: rscavm(nmodes)         ! Mass scavenging ratios

PUBLIC ukca_plume_scav, ukca_mode_scavcoeff, ukca_set_conv_indices

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SCAVENGING_MOD'

CONTAINS

! ######################################################################
SUBROUTINE ukca_set_conv_indices()
!
! set the entries in tracer_info to allow plume scavenging to 
! correctly locate MODE tracers in the tot_tracer array in convection

! we need to set up mmr_index_um and nmr_index_um - the location
! in the tracer_ukca array of the number and mass mixing ratios
! Call from addres - outside of all OpenMP loops and only called once

USE ereport_mod,               ONLY: ereport
USE errormessagelength_mod,    ONLY: errormessagelength
USE parkind1,                  ONLY: jprb, jpim
USE ukca_all_tracers_copy_mod, ONLY: arrayloc
! Note that all values used from ukca_mode_setup are parameters so will be
! available before call to UKCA
USE ukca_mode_setup,           ONLY: nmodes, ncp, mode_names, cp_su, cp_cl, &
                                     cp_bc, cp_oc, cp_du, cp_so, cp_nh4,    &
                                     cp_no3
USE ukca_nmspec_mod,           ONLY: nm_spec_active
USE um_parcore,                ONLY: mype
USE umprintmgr,                ONLY: umprint, ummessage, printstatus,       &
                                     prstatus_diag
USE yomhook,                   ONLY: lhook, dr_hook

IMPLICIT NONE

! Local variables
INTEGER :: icp                  ! component index
INTEGER :: imode                ! mode index
INTEGER :: icode                ! error code
INTEGER :: i, j                 ! loop counter
CHARACTER(LEN=errormessagelength) :: cmessage   ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SET_CONV_INDICES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(mmr_index_um)) ALLOCATE(mmr_index_um(nmodes,ncp))

! Construct the number and component mass index arrays
! by looping over all the active UKCA tracers
DO i=1,SIZE(nm_spec_active)
    IF ( ANY(mode_names(:) == nm_spec_active(i)(1:7)) ) THEN 
        ! Tracer belongs to UKCA_MODE aerosol scheme
        
        ! first 7 characters of tracer name are the mode name 
        ! imode is the position in the mode_names array
        imode = arrayloc(mode_names, nm_spec_active(i)(1:7))

        ! last 2 characters are the component name
        ! or 'ND' for number density
        SELECT CASE (nm_spec_active(i)(9:10))
        CASE ('ND')
          icp = 0
        CASE ('SU')
          icp = cp_su
        CASE ('SS')
          icp = cp_cl
        CASE ('BC')
          icp = cp_bc
        CASE ('OC')
          icp = cp_oc
        CASE ('DU')
          icp = cp_du
        CASE ('SO')
          icp = cp_so
        CASE ('NH')
          icp = cp_nh4
        CASE ('NT')
          icp = cp_no3
        CASE DEFAULT
          icode = 1
          cmessage = ' Unknown component in CASE statement'
          WRITE(ummessage,'(A40,A10)') cmessage,nm_spec_active(i)
          CALL umPrint(ummessage,src='ukca_mode_setup')
          CALL ereport('UKCA_MODE_SETUP.UKCA_AERO_TRACER_INIT',icode,      &
                        cmessage)
        END SELECT

        ! now use this information to set up the indices in 
        ! all tracers for number and mass
        IF (icp == 0) THEN
          ! This is a number density - set index in nmr_index
          ! We need the index in all_tracers so we have to use
          ! tr_lookup to map from UKCA tracer index to tot_tracer
          ! 
          nmr_index_um(imode) = i
        ELSE IF (icp > 0) THEN
          ! This is a component masss - set index in mmr_index_um
          mmr_index_um(imode,icp) = i
        END IF
  END IF
END DO

! print out nmr and mmr indices
IF (printstatus >= prstatus_diag .AND. mype ==0) THEN
  WRITE(ummessage,'(A)') 'number tracer maps'
  CALL umprint(ummessage,src=routinename)
  DO i=1, nmodes
    WRITE(ummessage,'(I6,1X,I6)') i, nmr_index_um(i)
    CALL umPrint(ummessage,src=routinename)
  END DO
END IF

IF (printstatus >= prstatus_diag .AND. mype ==0) THEN
  WRITE(ummessage,'(A)') 'mass tracer maps'
  CALL umprint(ummessage,src=routinename)
  DO j=1, ncp
      DO i=1,nmodes
        WRITE(ummessage,'(I6,1X,I6,1X,I6)') i, j, mmr_index_um(i, j)
        CALL umprint(ummessage,src=routinename)
      END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_set_conv_indices

! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_plume_scav(k, np_full, ntra, npnts, trapkp1, xpkp1,     &
                           prekp1, flxkp1)

! ----------------------------------------------------------------------
!  Calculate the change in parcel tracer content due to scavenging
!  by precipitation for UKCA tracers.
! ----------------------------------------------------------------------
USE ukca_mode_setup,     ONLY: nmodes, ncp, mode_names
USE ukca_nmspec_mod,     ONLY: nm_spec_active
USE planet_constants_mod, ONLY: g
USE parkind1,            ONLY: jprb, jpim
USE umprintmgr,          ONLY: umprint, ummessage
USE ereport_mod,         ONLY: ereport
USE yomhook,             ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: k                ! Present model layer
INTEGER, INTENT(IN) :: np_full          ! Full vector length
INTEGER, INTENT(IN) :: ntra             ! Number of UKCA tracer variables
INTEGER, INTENT(IN) :: npnts            ! Number of points


REAL, INTENT(IN)    ::  &
  xpkp1(npnts)                ! parcel cloud condensate in layer k+1 (kg/kg)
                              ! liquid water for [6A] convection
                              ! liquid + frozen water for [4A/5A] convection

REAL, INTENT(IN)    ::  &
  prekp1(npnts)               ! precipitation from parcel as it rises from
                              ! layer k to k+1 (kg/m**2/s)
REAL, INTENT(IN)    ::  &
  flxkp1(npnts)               ! parcel massflux in layer k+1 (Pa/s)

REAL, INTENT(INOUT) ::  &
  trapkp1(npnts,ntra)         ! parcel tracer content in layer k+1 (kg/kg)

! Local
INTEGER :: ktra               ! counter
INTEGER :: ktra_fixed         ! loop index 
INTEGER :: imode              ! index for mode 
INTEGER :: icp                ! index for component 

REAL :: traprekp1(npnts,ntra) ! Tracer content in precipitation produced
                              ! during ascent from layer K to K+1
                              ! [TR/kg(H2O), where tracer is TR/kg(air)]
REAL :: removal(npnts,ntra)   ! Removed tracer
REAL :: trapkp1_old(npnts,ntra) ! Temporary array to store original value of
                                ! trapkp1

INTEGER :: ifirst_ukca               ! index of first UKCA tracer in trapkp1
INTEGER :: ilast_ukca                ! index of last UKCA tracer in trapkp1
INTEGER :: iukca                     ! index of tracer in nm_spec_active
INTEGER :: icode                     ! error code
CHARACTER(LEN=errormessagelength) :: cmessage        ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PLUME_SCAV'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

traprekp1(:,:) = 0.0
removal(:,:) = 0.0


! Calculate the aqueous phase mixing ratio for the UKCA tracers
! only if precipitiation is greater than a threshold
IF (ANY(prekp1(:) > precip_lim)) THEN
  IF (i_convection_vn == i_convection_vn_6a) THEN
    CALL ukca_calc_aqueous_6a(k, npnts, ntra,                                &
                              trapkp1(1:npnts, :), xpkp1(:),                 &
                              prekp1(:), flxkp1(:),                          &
                              traprekp1(:,:) )
  ELSE IF (i_convection_vn == i_convection_vn_5a ) THEN
    CALL ukca_calc_aqueous_5a(k, npnts, ntra,                                &
                              trapkp1(1:npnts, :),                           &
                              xpkp1(:) + prekp1(:)*g/flxkp1(:),              &
                              traprekp1(:,:) )
  ELSE
    cmessage = ' Plume scavenging is not supported for this convection version'
    icode = 1
    WRITE(ummessage,'(A65,I4)') cmessage,i_convection_vn
    CALL umPrint(ummessage,src='ukca_scavenging_mod.ukca_plume_scav')
    CALL ereport('ukca_scavenging_mod.ukca_plume_scav',icode,cmessage)
  END IF       ! i_convection_vn
END IF     ! prekp1(:) > precip_lim)


! Remove scavenged tracer along with water

! set tracer indices in to loop over 
ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca  = tracer_info%i_ukca_last

! Check whether there is precipitation - if not do nothing
IF (ANY(prekp1(:) > precip_lim)) THEN

  ! loop over all UKCA tracers
  DO ktra = ifirst_ukca,ilast_ukca
  
    ! for names in nm_spec_active need to index relative
    ! to the ukca_tracer array
    iukca = ktra - ifirst_ukca + 1
    
    ! check whether this is a MODE tracer by matching the 
    ! first 7 characters of the tracer name against 
    ! the list of mode names
    IF ( ANY(mode_names(:) == nm_spec_active(iukca)(1:7)) ) THEN
    
        ! Only scavenge aerosols where there is a minimum
        ! amount of tracer in precip
        WHERE (traprekp1(:,ktra) > 1.0e-30)
        
          ! calculate the loss of tracer
          removal(:,ktra) = traprekp1(:,ktra) * prekp1(:)*g/flxkp1(:)
          ! store the previous value of the tracer 
          trapkp1_old(:,ktra) = trapkp1(1:npnts,ktra)
          ! now reduce the tracer concentration by the removal amount
          trapkp1(1:npnts,ktra) = trapkp1(1:npnts,ktra) - removal(:,ktra)
        END WHERE


        ! If negative tracer results, then limit the removal to give zero.
        WHERE (trapkp1(1:npnts,ktra) < 0.0)
          trapkp1(1:npnts,ktra) = 0.0
          removal(:,ktra) = trapkp1_old(:,ktra)
        END WHERE

    END IF   ! Mode tracer
  END DO     ! ktra
END IF       ! ANY(prekp1 > 0)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_plume_scav

! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_calc_aqueous_5a(k, npnts, ntra, trapkp1, qcl, &
                                aqueous)

!   Calculates the mixing ratio of UKCA tracers in the aqueous phase
!   for scavenging. For [4A/5A] convection scheme.

! UKCA_D1_DEFS etc won't have been initialised on the first call,
!  as convection happens before the first call to ukca_main1, so
!  we must only use constants here. See ukca_wetdp_addr

USE ukca_mode_setup, ONLY: nmodes, ncp, mode, component
USE umPrintMgr,      ONLY: umPrint, ummessage
USE ereport_mod,     ONLY: ereport
USE Control_Max_Sizes
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: k                    ! Level
INTEGER, INTENT(IN) :: npnts                ! Number of grid boxes
INTEGER, INTENT(IN) :: ntra                 ! Number of UKCA MODE tracers
REAL, INTENT(IN)    :: trapkp1(npnts,ntra)  ! parcel tracer content in
                                            ! layer k+1 (kg/kg)
                                            ! n.b. all model tracers
REAL, INTENT(IN)    :: qcl(npnts)           ! Liquid water mixing ratio
                                            ! [kg(liquid water)/kg(air)] 

! Outputs
REAL, INTENT(OUT)   :: aqueous(npnts,ntra)  ! Tracer/water mixing ratio
                                            ! [TR(aq)/kg(water)]

! Local variables
INTEGER :: j                    ! array index
INTEGER :: ktra                 ! counter
INTEGER :: ktra_fixed           ! loop index 
INTEGER :: icp                  ! component index
INTEGER :: imode                ! mode index
! position of first UKCA tracer in all tracers array
INTEGER :: ifirst_ukca
! position of last UKCA tracer in all tracers array
INTEGER :: ilast_ukca
INTEGER :: icode                ! error code
CHARACTER (LEN=* ), PARAMETER :: RoutineName='UKCA_CALC_AQUEOUS_5A'
CHARACTER(LEN=errormessagelength) :: cmessage   ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set this to zero for any other tracers (e.g. chemistry)
aqueous(:,:) = 0.0

ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca = tracer_info%i_ukca_last

DO imode=1,nmodes
  IF (mode(imode)) THEN
    IF (nmr_index_um(imode) > 0) THEN
      ! For number tracers, [tracer] = [ptcl/molec(air)]
      !                              = [ptcl/kg(air)] * Boltzmann / R
      !                                i.e. [TR] = [ptcl] * Boltzmann / R

      ! Only calculate where non-negligible liquid water
      j = nmr_index_um(imode) + ifirst_ukca - 1
      IF (j < ifirst_ukca .OR. j > ilast_ukca) THEN
        cmessage = 'j is out of range for aerosol number'
        icode = ABS(j)
        CALL umPrint(cmessage,src=RoutineName)
        WRITE(ummessage,'(A,I5,1X,A,I5)') 'imode', imode, 'j= ', j
        CALL umPrint(ummessage,src=RoutineName)
        WRITE(ummessage,'(A,I5,1X,I5)') 'should be between ', &
            ifirst_ukca, ilast_ukca
        CALL umPrint(ummessage,src=RoutineName)
        WRITE(ummessage,'(A,I5)') 'icode: ',icode
        CALL umPrint(ummessage,src=RoutineName)
        CALL ereport(RoutineName,icode,cmessage)
      END IF
      WHERE (qcl(1:npnts) > 1e-10)
        aqueous(:,j) = rscavn(imode) *  trapkp1(1:npnts,j) / qcl(:)
      END WHERE
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (mmr_index_um(imode,icp) > 0) THEN
            ! For mass tracers, [tracer] = [kg/kg(air)] i.e. [TR] = [kg]
            j = mmr_index_um(imode,icp) + ifirst_ukca - 1
            IF (j < ifirst_ukca .OR. j > ilast_ukca) THEN
              cmessage = 'j is out of range'
              icode = ABS(j)
              CALL umPrint(cmessage,src=RoutineName)
              WRITE(ummessage,'(A,I6,1X,I6,1x,A,I6)') 'imode, icp',  &
                  imode, icp, 'j= ', j
              CALL umPrint(ummessage,src=RoutineName)
              WRITE(ummessage,'(A,I5,1X,I5)') 'should be between: ', &
                  ifirst_ukca, ilast_ukca
              CALL umPrint(ummessage,src=RoutineName)
              WRITE(ummessage,'(A,I5)') 'icode: ',icode
              CALL umPrint(ummessage,src=RoutineName)
              CALL ereport(RoutineName,icode,cmessage)
            END IF
            WHERE (qcl(:) > 1e-10)
              aqueous(:,j) = rscavm(imode) * trapkp1(1:npnts,j) / qcl(:)
            END WHERE
          END IF
        END IF
      END DO    ! icp
    END IF     ! nmr_index_um > 0
  END IF    ! mode(imode)
END DO     ! imode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_aqueous_5a

! ######################################################################
! Subroutine Interface:
SUBROUTINE ukca_calc_aqueous_6a(k, npnts, ntra, trapkp1,                &
                                xpkp1, prekp1, flxkp1,                  &
                                aqueous)

!   Calculates the mixing ratio of UKCA tracers in the aqueous phase
!   for scavenging. For [6A] convection scheme.

! UKCA_D1_DEFS etc won't have been initialised on the first call,
!  as convection happens before the first call to ukca_main1, so
!  we must only use constants here. See ukca_wetdp_addr

USE ukca_mode_setup, ONLY: nmodes, ncp, mode, component
USE ukca_nmspec_mod, ONLY: nm_spec_active
USE ukca_option_mod, ONLY: mode_aitsol_cvscav, i_mode_nucscav, l_ukca_advh2o
USE planet_constants_mod, ONLY: g
USE umPrintMgr,      ONLY: umPrint, ummessage
USE ereport_mod,     ONLY: ereport
USE Control_Max_Sizes
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: k                    ! Level
INTEGER, INTENT(IN) :: npnts                ! Number of grid boxes
INTEGER, INTENT(IN) :: ntra                 ! Number of UKCA MODE tracers
REAL, INTENT(IN)    :: trapkp1(npnts,ntra)  ! parcel tracer content in
                                            ! layer k+1 (kg/kg)
                                            ! n.b. all model tracers
REAL, INTENT(IN)    ::  &
  xpkp1(npnts)                ! parcel cloud water in layer k+1 (kg/kg)

REAL, INTENT(IN)    ::  &
  prekp1(npnts)               ! precipitation from parcel as it rises from
                              ! layer k to k+1 (kg/m**2/s)
REAL, INTENT(IN)    ::  &
  flxkp1(npnts)               ! parcel massflux in layer k+1 (Pa/s)


! Outputs
REAL, INTENT(OUT)   :: aqueous(npnts,ntra)  ! Tracer/water mixing ratio
                                            ! [TR(aq)/kg(water)]

! Local variables

REAL, PARAMETER :: qcl_min = 1e-15      ! Minimum qcl for calculations
REAL, PARAMETER :: trc_min = 1e-20      ! Minimum tracer for calculations
  
REAL          :: qcl(npnts)             ! Liquid water mixing ratio
                                        ! [TR/kg(air)] 
INTEGER :: j                    ! array index
INTEGER :: ktra                 ! counter
INTEGER :: ktra_fixed           ! loop index 
INTEGER :: icp                  ! component index
INTEGER :: imode                ! mode index
INTEGER :: ifirst_ukca          ! position of first UKCA tracer in
                                ! all tracers array
INTEGER :: ilast_ukca           ! position of last UKCA tracer in
                                ! ukca tracer array
INTEGER :: icode                ! error code
CHARACTER (LEN=* ), PARAMETER :: RoutineName='UKCA_CALC_AQUEOUS_6A'
CHARACTER(LEN=errormessagelength) :: cmessage   ! error message
LOGICAL :: todo(npnts)          ! mask for calculation

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set this to zero for any other tracers (e.g. chemistry)
aqueous(:,:) = 0.0

!  Consistency check to avoid floating point error
WHERE (prekp1 > precip_lim .AND. flxkp1 > precip_lim)
  qcl(:) = xpkp1(:) + prekp1(:)*g/flxkp1(:)
  todo(:) = .TRUE.
ELSE WHERE
  todo(:) = .FALSE.
  qcl(:) = 0.0
END WHERE

ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca = tracer_info%i_ukca_last

DO imode=1,nmodes
  IF (mode(imode)) THEN
    IF (nmr_index_um(imode) > 0) THEN
      ! For number tracers, [tracer] = [ptcl/molec(air)]
      !                              = [ptcl/kg(air)] * Boltzmann / R
      !                                i.e. [TR] = [ptcl] * Boltzmann / R

      j = nmr_index_um(imode) + ifirst_ukca - 1
      IF (j < ifirst_ukca .OR. j > ilast_ukca) THEN
        cmessage = 'j is out of range for aerosol number'
        icode = ABS(j)
        CALL umPrint(cmessage,src=RoutineName)
        WRITE(ummessage,'(A,I5,1X,A,I5)') 'imode', imode, 'j= ', j
        CALL umPrint(ummessage,src=RoutineName)
        WRITE(ummessage,'(A,I5,1X,I5)') 'should be between ', &
            ifirst_ukca, ilast_ukca
        CALL umPrint(ummessage,src=RoutineName)
        WRITE(ummessage,'(A,I5)') 'icode: ',icode
        CALL umPrint(ummessage,src=RoutineName)
        CALL ereport(RoutineName,icode,cmessage)
      END IF
      WHERE (todo .AND. trapkp1(1:npnts,j) > trc_min .AND. qcl(:) > qcl_min)
        aqueous(:,j) = rscavn(imode) * trapkp1(1:npnts,j) / qcl(:)
      END WHERE
      DO icp=1,ncp
        IF (component(imode,icp)) THEN
          IF (mmr_index_um(imode,icp) > 0) THEN
            ! For mass tracers, [tracer] = [kg/kg(air)] i.e. [TR] = [kg]
            j = mmr_index_um(imode,icp) + ifirst_ukca - 1
            IF (j < ifirst_ukca .OR. j > ilast_ukca) THEN
              cmessage = 'j is out of range'
              icode = ABS(j)
              CALL umPrint(cmessage,src=RoutineName)
              WRITE(ummessage,'(A,I6,1X,I6,1x,A,I6)') 'imode, icp',  &
                  imode, icp, 'j= ', j
              CALL umPrint(ummessage,src=RoutineName)
              WRITE(ummessage,'(A,I5,1X,I5)') 'should be between: ', &
                  ifirst_ukca, ilast_ukca
              CALL umPrint(ummessage,src=RoutineName)
              WRITE(ummessage,'(A,I5)') 'icode: ',icode
              CALL umPrint(ummessage,src=RoutineName)
              CALL ereport(RoutineName,icode,cmessage)
            END IF
            WHERE (todo .AND. trapkp1(1:npnts,j) > trc_min .AND.            &
              qcl(:) > qcl_min)
              aqueous(:,j) = rscavm(imode) * trapkp1(1:npnts,j) / qcl(:)
            END WHERE
          END IF
        END IF
      END DO    ! icp
    END IF     ! nmr_index > 0
  END IF    ! mode(imode)
END DO     ! imode

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_aqueous_6a

! ######################################################################
SUBROUTINE ukca_mode_scavcoeff()
! Set the scavenging co-efficients to be used in plume scavenging
! set depending on the switch i_mode_nucscav set in ukca_option_mod
! call after namelists read 
! ######################################################################

USE ereport_mod,               ONLY: ereport
USE umprintmgr,                ONLY: umprint, ummessage
USE errormessagelength_mod, ONLY: errormessagelength
USE ukca_option_mod, ONLY: mode_aitsol_cvscav, i_mode_nucscav
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

IMPLICIT NONE

LOGICAL, PARAMETER :: l_convective = .TRUE. ! True if convective cloud

INTEGER :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=errormessagelength) :: cmessage   ! error message
CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MODE_SCAVCOEFF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (i_mode_nucscav)
  CASE (1,3)
    IF (l_convective) THEN
      ! Original GLOMAP-MODE scavenging ratios, plus scavenging
      ! of part of the soluble Aitken mode due to the stronger
      ! updraughts and higher supersaturations in convective cloud.
      !              1=NS  2=KS   3=AS  4=CS  5=KI  6=AI  7=CI
      rscavn(:) = (/ 0.00, mode_aitsol_cvscav, 1.00, 1.00, 0.00, 0.00, &
                     0.00 /)
      rscavm(:) = (/ 0.00, mode_aitsol_cvscav, 1.00, 1.00, 0.00, 0.00, &
                     0.00 /)
    ELSE
      ! Original GLOMAP-MODE scavenging ratios
      !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI
      rscavn(:) = (/ 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00 /)
      rscavm(:) = rscavn(:)
    END IF
  CASE (2)
    ! ECHAM5-HAM scavenging ratios (Stier et al. 2005)
    IF (l_convective) THEN
      ! for convective mixed cloud
      !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI
      rscavn(:) = (/ 0.20, 0.60, 0.99, 0.99, 0.20, 0.40, 0.40 /)
    ELSE
      ! for stratiform liquid cloud
      !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI
      rscavn(:) = (/ 0.10, 0.25, 0.85, 0.99, 0.20, 0.40, 0.40 /)
    END IF
    rscavm(:) = rscavn(:)
  CASE DEFAULT
    cmessage = 'i_mode_nucscav out of range'
    icode = 1
    WRITE(ummessage,'(A40,I8)') cmessage,i_mode_nucscav
    CALL umPrint(ummessage,src=RoutineName)
    CALL ereport(RoutineName,icode,cmessage)
END SELECT ! i_mode_nucscav

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_mode_scavcoeff

END MODULE ukca_scavenging_mod
