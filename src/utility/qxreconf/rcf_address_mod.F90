! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates addressing of the output dump

MODULE Rcf_Address_Mod

!  Subroutine Rcf_Address  - Loops over items and calls Rcf_Primary
!  Subroutine Rcf_Primary  - Calculated data lengths and addresses
!  Function Rcf_Disct_Lev  - Tests if level is discrete/continuous
!  Subroutine Rcf_PSLevCod - Decodes pseudo-level sizes
!
! Description:
!   Calculates the lengths, levels and addresses for the output
!   dump.
!
! Method:
!   Rcf_Address loops over STASHmaster records, and for those which
!   may be primary calls Rcf_Primary. This determines whether or not
!   the field is requied (via tstmsk) and all the related sizes and
!   addressing.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.

! These are all *current* values from a STASHmaster file.

USE Rcf_Ppx_Info_Mod, ONLY:    &
    STM_Record_Type,            &
    Ppxptr,                     &
    Ppxref_Sections,            &
    Ppxref_Items,               &
    STM_OptCodeLen

! Some variables (which are also in cstash_mod) that are used only
! for addressing purposes (mostly in Address and Primary).

USE Rcf_Address_Vars_Mod, ONLY: &
    ispace,                & ! Space code
    igp,                   & ! Grid of data code
    ilev,                  & ! Level type code
    ibot,                  & ! First level code
    itop,                  & ! Last level code
    iflag,                 & ! Level compression flag
    iopn,                  & ! Sectional option code
    vmsk,                  & ! Integer equiv of bin vers mask
    ipseudo,               & ! Pseudo dimension type
    ipfirst,               & ! First pseudo dim code
    iplast,                & ! Last pseudo dim code
    halo                     ! Halo type code

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE Ereport_Mod, ONLY: &
    Ereport

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ADDRESS_MOD'

CONTAINS

!-------------------------------------------------------------------
! This routine determines the Addressing of the output dump for
! reconfiguration - includes calls to primary etc to determine
! which fields are included.
!--------------------------------------------------------------------

SUBROUTINE Rcf_Address(input_lookup)

USE Rcf_Exppx_Mod, ONLY:      &
    Rcf_Exppx

USE submodel_mod, ONLY:                                                        &
    internal_model_list, n_internal_model_max, n_internal_model,               &
    submodel_for_im, n_submodel_partition_max

USE stash_model_mod, ONLY: len_extra

USE Rcf_NRecon_Mod             ! Use all of this...

USE rcf_nlist_recon_technical_mod, ONLY: &
    Var_Recon,                           &
    select_output_fields,                &
    defined_by_namelist,                 &
    interp_all_fields,                   &
    tstmsk_to_decide,                    &
    max_num_fields,                      &
    output_field_list

USE nlsizes_namelist_mod, ONLY: tr_ukca
USE ukca_option_mod,      ONLY: tr_ukca_a
USE ukca_tracer_stash,    ONLY: a_max_ukcavars
USE um_stashcode_mod,     ONLY: stashcode_prog_sec, stashcode_tracer_sec, &
                                stashcode_ukca_sec, stashcode_glomap_clim_sec

USE ereport_mod, ONLY: ereport
USE lookup_addresses, ONLY: item_code

USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE


! Arguments
INTEGER, INTENT(IN) :: input_lookup(:,:)

! Local scalars:
INTEGER    :: Im_ident  !Internal model identifier (absolute)
INTEGER    :: Im_index  !Internal model index (expt. dependent)
INTEGER    :: Sm_ident  !Submodel identifier (absolute)
INTEGER    :: isec
INTEGER    :: iitm
INTEGER    :: rlevs
INTEGER    :: raddress
INTEGER    :: PIrow
INTEGER    :: i
INTEGER    :: ErrorStatus
INTEGER    :: item_code_val
LOGICAL    :: laddr
LOGICAL    :: lmask

CHARACTER (LEN=errormessagelength) :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = "RCF_ADDRESS"

! Local arrays:
!    Submodel definitions array: stores list of Im_index's
!     for each submodel partition
INTEGER :: SM_def(n_submodel_partition_max,n_internal_model_max)

! Local structures
TYPE (STM_Record_Type) :: STM_Rec
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! 1.  Set STASHIN addresses and input lengths for primary fields
!   The address loop for primary fields is performed for each
!   internal model in turn. Hence, each internal model's primary
!   data occupies a contiguous block in D1. The order of these blocks
!   is the same as the order of the internal models given in the
!   array INTERNAL_MODEL_LIST.
!   User-defined prognostics are included in this primary addressing
!   routine, since they are incorporated into the ppxref lookup
!   arrays PPXI, PPXC in routine GETPPX.

!   Initialisation

PrimDataLen(:)   = 0
SM_def( :, : )   = 0

!   Obtain submodel definitions and store in SMdef array
DO Im_index = 1,n_internal_model
   !   Submodel ident.
  Sm_ident =   submodel_for_im(Im_index)
  !   Internal model index
  SM_def(Sm_ident,Im_index) = Im_index
END DO

!   Primary address loop

!     Loop over submodel partitions
DO  Sm_ident = 1,n_submodel_partition_max
  ! Initialise len_extra
  len_extra(Sm_ident) = 0

  !     Initialise address for reconfiguration
  raddress = 1

  !      Loop over internal models for each SM partition
  DO Im_index = 1,n_internal_model

    !       Test whether current SM contains this IM
    IF (SM_def(Sm_ident,Im_index) > 0) THEN

      !        Obtain internal model identifier
      Im_ident   = internal_model_list(Im_index)

      PIrow  = 0
      ! Now that reconfiguration can take place over any section for VAR loop
      ! over all sections
      DO isec = 0,ppxref_sections
        IF (isec == stashcode_prog_sec .OR. & ! Sect. 0 prognostic variables
            isec == stashcode_tracer_sec .OR. & ! Free tracers
            isec == stashcode_ukca_sec .OR.  & ! UKCA tracers
            isec == stashcode_glomap_clim_sec .OR. & ! GLOMAP_MODE NWP
            Var_Recon  .OR.                   & ! VAR reconfiguration
            select_output_fields /= tstmsk_to_decide) THEN
          !       Loop over wanted section items
          recondat_node => RecondatList( im_index, isec )
          DO iitm   = 1,ppxref_items
            !   Check whether there is a primary field corresponding
            !   to this item number
            IF ( ppxptr( Im_ident, isec, iitm) /= 0) THEN
              STM_Rec = Rcf_Exppx( Im_ident, Isec, Iitm )
              vmsk    = STM_Rec % version_mask
              ispace  = STM_Rec % space_code
              igp     = STM_Rec % grid_type
              ilev    = STM_Rec % lv_code
              ibot    = STM_Rec % lb_code
              itop    = STM_Rec % lt_code
              halo    = STM_Rec % halo_type
              DO i = 1, STM_OptCodeLen / 5
                iopn(i) = STM_Rec % opt_code(i)
              END DO
              iflag   = STM_Rec % lev_flag
              ipseudo = STM_Rec % pt_code
              ipfirst = STM_Rec % pf_code
              iplast  = STM_Rec % pl_code

              IF (select_output_fields == interp_all_fields) THEN
                DO i = LBOUND(input_lookup,2), UBOUND(input_lookup,2)
                  item_code_val = input_lookup(item_code,i)
                  ! Determine if current field is present in input
                  ! dump
                  IF (isec == item_code_val/1000 .AND.  &
                       iitm == MOD(item_code_val,1000)) THEN
                    CALL rcf_primary(isec,iitm,im_index,im_ident, &
                         sm_ident,raddress,pirow)
                    ! Now exit loop.
                    EXIT
                  END IF
                END DO
              ELSE IF (select_output_fields == defined_by_namelist) THEN
                DO i = 1, max_num_fields
                  IF( output_field_list(i) /= imdi ) THEN
                    ! Determine if current field is present in namelist
                    ! array
                    IF (isec == output_field_list(i)/1000 .AND.  &
                        iitm == MOD(output_field_list(i),1000)) THEN
                      CALL rcf_primary(isec,iitm,im_index,im_ident, &
                           sm_ident,raddress,pirow)
                      ! Now exit loop.
                      EXIT
                    END IF
                  END IF
                END DO
              ELSE IF (select_output_fields == tstmsk_to_decide) THEN
                IF ( (ispace == 2) .OR. (ispace == 3) .OR. (ispace == 9) &
                     .OR. (ispace == 5)  .OR. (ispace == 10) &
                     .OR. (ispace == 8) .OR. Var_Recon) THEN !Primary variable
                  ! VAR reconfiguration in 4D-VAR needs to include
                  ! any section in the
                  ! reconfiguation if specified by option code.
                  ! Find out whether the primary is included for this version
                  ! DEPENDS ON: tstmsk
                  CALL tstmsk(Im_ident,isec,lmask,laddr,ErrorStatus,cmessage)

                  ! If this is a UKCA tracer default it to FALSE i.e. off
                  IF ( isec == stashcode_ukca_sec .AND. &
                       iitm <= a_max_ukcavars ) THEN
                    tr_ukca_a(iitm) = .FALSE.
                  END IF
                  ! If tstmsk determines that the current field should be
                  ! present it sets laddr to true
                  IF (laddr) THEN

                    ! If this is a UKCA tracer it is on becuase tstmsk
                    ! returned TRUE for it. Turn the tracer on and
                    ! increment the number of active UKCA tracers by 1
                    IF ( isec == stashcode_ukca_sec .AND. &
                         iitm <= a_max_ukcavars ) THEN
                      tr_ukca_a(iitm) = .TRUE.
                      tr_ukca = tr_ukca + 1
                    END IF

                    CALL Rcf_Primary( isec, iitm, Im_index, Im_ident, &
                         Sm_ident, raddress, PIrow )
                  END IF
                END IF
              ELSE
                 WRITE(cmessage, '(A,I8)') "Invalid setting for select_" // &
                      "output_fields : ", select_output_fields
                 ErrorStatus = 10
                 CALL Ereport( RoutineName, ErrorStatus, Cmessage )
              END IF
            END IF    !  PPXPTR(m,s,i) /= 0
          END DO    !  Loop over items
          recondat_node => NULL()
        END IF ! ISEC or VAR
      END DO    !  Loop over 'primary' sections
    END IF    !  test whether SM contains IM
  END DO     !  Loop over Im_index
END DO      !  Loop over SM partitions

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Address

!===================================================================

! Compute data lengths and addresses for primary fields
! Subroutine Interface:

SUBROUTINE Rcf_Primary( isec, iitm, Im_index, Im_ident, &
                        Sm_ident, raddress, PIrow )

USE Rcf_NRecon_Mod, ONLY:  &
    Recondat_Node,           &
    RecondatList,           &
    DumpProgLevs,       &
    PrimDataLen

USE rcf_nlist_recon_technical_mod, ONLY: &
    tstmsk_to_decide, select_output_fields, var_recon

USE nlsizes_namelist_mod, ONLY: &
    ntiles,             &
    tr_vars

USE Rcf_Level_Code_Mod, ONLY: &
    Rcf_Level_Code

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Address_Length_Mod, ONLY: &
    Rcf_Address_Length

USE jules_sea_seaice_mod, ONLY: &
    nice

USE um_stashcode_mod,     ONLY: stashcode_prog_sec, stashcode_tracer_sec, &
                                stashcode_ukca_sec, stashcode_glomap_clim_sec

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN)    :: isec     ! Current section number
INTEGER, INTENT(IN)    :: iitm     ! Current section 0 item number
INTEGER, INTENT(IN)    :: Im_ident ! Current internal model number
INTEGER, INTENT(IN)    :: Im_index ! pos'n in internal model list
INTEGER, INTENT(IN)    :: Sm_ident ! Submodel identifier (absolute)
INTEGER, INTENT(INOUT) :: raddress   ! Address for reconfiguration
INTEGER, INTENT(INOUT) :: PIrow    ! Counter for ProgItems array

! Local variables
LOGICAL            :: model_lev
INTEGER            :: rlevs      ! No. of levels for reconfiguration
INTEGER            :: rplevs     ! no. of pseudo-levels
INTEGER            :: il1,il2
INTEGER            :: ipl1,ipl2
INTEGER            :: LEN        ! Data length for primary field
INTEGER            :: ErrorStatus
CHARACTER (LEN=errormessagelength) :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_PRIMARY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!- End of Header ---------------------------------------------------

IF (ispace /= 10) THEN

  ! Find address length per level
  CALL Rcf_Address_Length( igp, halo, LEN )

  ! If tracer, multiply length by number of tracer vars for
  ! LBC calculation
  IF ( ( igp  ==  26 ) .AND. ( itop  ==  11 ) ) THEN
    LEN = LEN * tr_vars
  END IF

  ! If level type is 0, is LBC so multiply size by levels
  IF ( ilev  ==  0 ) THEN
    CALL Rcf_Level_Code( ibot, il1, Output_Grid )
    CALL Rcf_Level_Code( itop, il2, Output_Grid )
    LEN = LEN * (il2 - il1 + 1)
  END IF

  model_lev = Rcf_Disct_Lev( ilev )
  IF (model_lev .OR. (ilev == 5 .AND. ipseudo /= 0)) THEN
    ! Field has model levels - decode level codes
    IF (ilev  /=  5) THEN
      CALL Rcf_Level_Code( ibot, il1, Output_Grid )
      CALL Rcf_Level_Code( itop, il2, Output_Grid )
    ELSE
      il1=1
      il2=1
    END IF

    ! No. of model levels (for reconfiguration)
    rlevs=il2-il1+1
    ! Initialise first & last pseudo level indices
    ipl1 =0
    ipl2 =0

    IF (Iflag == 0 .AND. Ipseudo /= 0) THEN
      ! Primary with input on all available pseudo levels -
      !   decode pseudo level codes
      CALL Rcf_PSLevCod( ipfirst, ntiles, nice, ipl1, 'F' )
      CALL Rcf_PSLevCod( iplast , ntiles, nice, ipl2, 'L' )
    END IF

    rplevs=ipl2-ipl1+1
    ! Multiply length per level by no. of levels
    LEN=len*(il2-il1+1)*(ipl2-ipl1+1)

    ! Ignore items with space code 9 unless using var_recon or 
    ! overidding the tstmsk logic that defines which fields should be in
    ! the output dump
    IF (( ispace /= 9 )                                              &
         .OR. var_recon .OR. select_output_fields /= tstmsk_to_decide) THEN
      ! Increment no. of headers
      DumpProgLevs   (Im_ident)=   DumpProgLevs(Im_ident) &
                                    +(il2-il1+1)*(ipl2-ipl1+1)
    END IF

  ELSE  ! Not model levels
    rlevs=1
    rplevs=1
    ! Ignore items with space code 9 unless using var_recon or 
    ! overidding the tstmsk logic that defines which fields should be in
    ! the output dump
    IF (( ispace /= 9 ) .OR.                                          &
       var_recon .OR. select_output_fields /= tstmsk_to_decide) THEN
      DumpProgLevs   (Im_ident)=DumpProgLevs   (Im_ident)+1
    END IF
  END IF

  ! Addresses are set relative to the beginning of the primary data,
  ! since the primary data starts at the beginning of D1.
  ! Ignore items with space code 9 unless using var_recon or 
  ! overriding the tstmsk logic that defines which fields should be in
  ! the output dump
  IF (ispace /= 9 .OR. var_recon .OR.                       &
       select_output_fields /= tstmsk_to_decide) THEN
    ! Start address for this primary field
    ! Increment len_prim by LEN (=data length for this primary field)
    PrimDataLen(Im_ident)      =PrimDataLen(Im_ident)+LEN
  END IF

  ! Store levels, lengths and addresses required for reconfiguration
  !                                                in array Recondat
  ! Ignore items with space code 9 unless using var_recon or 
  ! overidding the tstmsk logic that defines which fields should be in
  ! the output dump
  IF ((ispace /= 9) .OR.                                   &
       var_recon .OR.                                      &
       select_output_fields /= tstmsk_to_decide) THEN
    IF (isec == stashcode_prog_sec .OR. isec == stashcode_tracer_sec .OR.      &
        isec == stashcode_ukca_sec .OR. isec == stashcode_glomap_clim_sec .OR. &
        var_recon .OR.                                                         &
        select_output_fields /= tstmsk_to_decide) THEN
      ALLOCATE(recondat_node % recondat_info)
      recondat_node % recondat_info % sec_item=iitm + 1000*isec
      recondat_node % recondat_info % rlevs=rlevs
      recondat_node % recondat_info % LEN=len
      recondat_node % recondat_info % raddress=raddress
      recondat_node % recondat_info % rplevs=rplevs

      ALLOCATE(recondat_node % next)
      recondat_node => recondat_node % next
    ELSE
      ErrorStatus=1
      Cmessage='Rcf_Primary : Invalid value for ISEC, section no.'
      CALL Ereport( 'Rcf_Primary', ErrorStatus, Cmessage )
    END IF
    raddress = raddress+LEN
  END IF
END IF  ! ISPACE /= 10
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Primary

!====================================================================

! Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
LOGICAL FUNCTION Rcf_Disct_Lev( lev_code )

USE Ereport_Mod, ONLY: &
    Ereport

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Function arguments:
INTEGER, INTENT(IN)          :: lev_code !Level code from STASHmaster

! ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'Disct_Lev'
INTEGER                      :: ErrorStatus
CHARACTER (LEN=errormessagelength)           :: cmessage

!- End of Header ----------------------------------------------

IF (lev_code == 1 .OR. lev_code == 2 .OR. lev_code == 6 .OR. &
                       lev_code == 10) THEN
  Rcf_Disct_Lev=.TRUE.
ELSE IF (lev_code  >=  0 .AND. lev_code  <=  10) THEN
  Rcf_Disct_Lev=.FALSE.
ELSE
  Rcf_Disct_Lev=.FALSE.
  ErrorStatus=1
  cmessage='Rcf_DISCT_LEV : Invalid level type in STASHmaster'

  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

END FUNCTION Rcf_Disct_Lev

!===================================================================

! Decode the STASH pseudo level code
! Subroutine Interface:
SUBROUTINE Rcf_PSLevCod( ilin, ntiles, nice, ilout, swtch )

USE Ereport_Mod, ONLY: &
    Ereport

USE rad_input_mod, ONLY: &
    h_swbands,            &
    h_lwbands
USE aodband_mod, ONLY: aod_bands

! JULES
USE jules_snow_mod, ONLY: nsmax

USE jules_surface_types_mod
USE clmchfcg_scenario_mod, ONLY: nsulpat

IMPLICIT NONE

! Description:
!   Sets ILOUT to an appropriate pseudo level size according
!   to the value of IL
!
! Subroutine arguments:
INTEGER, INTENT(IN)           :: ilin    ! Model pseudo level code
INTEGER, INTENT(IN)           :: ntiles  ! Number of surface tiles
INTEGER, INTENT(IN)           :: nice    ! Number of seaice catagories
CHARACTER (LEN=1), INTENT(IN) :: swtch
INTEGER, INTENT(OUT)          :: ilout   ! An actual pseudo level

! Local variables
INTEGER              :: ErrorStatus
CHARACTER (LEN=80)   :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_PSLEVCOD'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of Header --------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (swtch == 'F') THEN
  IF (ilin == 1) THEN
    ilout=1
    ! Ocean assimilation groups - removed
  ELSE
    WRITE(umMessage,*) &
         'MSG FROM RCF_PSLEVCOD: ', &
         'INAPPROPRIATE FIRST PSEUDO LEVEL CODE FOUND ',ilin
    CALL umPrint(umMessage,src='rcf_address_mod')
    ErrorStatus=2
    Cmessage = 'Error from PSLevCod!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END IF
ELSE IF (swtch == 'L') THEN
  IF (ilin == 1) THEN
    ilout=h_swbands
  ELSE IF (ilin == 2) THEN
    ilout=h_lwbands
  ELSE IF (ilin == 3) THEN
    ! Radiation bands for measuring aerosol optical depth
    ilout=aod_bands
  ELSE IF (ilin == 4) THEN
    ! Last frequency (wave model)
    ! Wave model not used in RCF
    !    ILOUT=NFRE
    ilout=0
  ELSE IF (ilin == 5) THEN
    ! Last wave train (wave model)
    !   ILOUT=NWTRAIN
    ilout=0
  ELSE IF ( ilin  ==  6 ) THEN
    ! Last index for HadCM2 sulphate loading patterns.
    ilout = nsulpat
  ELSE IF ( ilin  ==  7 ) THEN
    ! All surface types
    ilout = ntype
  ELSE IF ( ilin  ==  8 ) THEN
    ! Plant functional types only
    ilout = npft
  ELSE IF ( ilin  ==  9 ) THEN
    ! All tiles
    ilout = ntiles
  ELSE IF ( ilin  ==  10 ) THEN
    ! All seaice catagories
    ilout = nice
  ELSE IF ( ilin  ==  11 ) THEN
    !New snow scheme: all tiles <times> max no. of snow levels
    ilout = ntiles*nsmax
  ELSE
    WRITE(umMessage,*) &
         'MSG FROM RCF_PSLEVCOD: ', &
         'INAPPROPRIATE LAST PSEUDO LEVEL CODE FOUND ',ilin
    CALL umPrint(umMessage,src='rcf_address_mod')
    ErrorStatus=2
    Cmessage = 'Error from PSLevCod'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_PSLevCod

END MODULE Rcf_Address_Mod
