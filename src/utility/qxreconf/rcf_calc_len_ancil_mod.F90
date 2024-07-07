! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate work space required for ancillary data.

MODULE Rcf_calc_len_ancil_Mod

IMPLICIT NONE
!  Subroutine Rcf_Calc_len_ancil  - Calculate work space for anc. data.
!
! Description:
!    Determines the work space required for the ancillary data.
!
! Method:
!    For each ancillary fields to be read in (SOURCE=2 in ITEMS
!    namelist), determine length of data to be read in and accumulate
!    to return overall workspace required.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_LEN_ANCIL_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_len_ancil (P_Field, R_Field, P_Rows, Len_Ancil)

USE Rcf_NRecon_Mod, ONLY: &
    ReconDatList,              &
    Recondat_Node

USE ozone_inputs_mod, ONLY: &
    zon_av_ozone

USE um_stashcode_mod, ONLY: &
    stashcode_ozone,           &
    stashcode_riv_sequence,    &
    stashcode_riv_direction,   &
    stashcode_riv_storage,     &
    stashcode_prog_sec

USE Ancil_Mod, ONLY: &
    ancil_requests,num_ancil_requests

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Diag

USE Ereport_Mod, ONLY:   &
    Ereport

USE UM_ParCore, ONLY: &
    mype

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER  :: P_Field         !  Length of field
INTEGER  :: R_Field         !  Length of River Routing field
INTEGER  :: P_Rows          !  No of rows
INTEGER  :: Len_Ancil       !  Total workspace for anc data

! Local variables
INTEGER  :: i,irec,j1       !  Loop indices
INTEGER  :: N_Levs          !  No of levels
INTEGER  :: N_Pseudo_Levs   !  No of pseudo levels
INTEGER  :: Len_anc_data    !  Length of anc data
INTEGER  :: StashCode       !  Stash Code
INTEGER  :: SectionCode     !  Section Code
INTEGER  :: Sec_Item        !  1000*Section + Stashcode
INTEGER  :: ErrorStatus     !  Error return code

CHARACTER (LEN=errormessagelength) :: CMessage      ! Error return message
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CALC_LEN_ANCIL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0
CMessage = ' '

Len_Ancil = 0

i = 1   ! hard-wired setting of submodel_ident (ATMOS)

DO irec = 1, num_ancil_requests

  len_anc_data = 0

  Recondat_Node => RecondatList(i,ancil_requests(irec)%section)
  ! Recondat is ordered by section and contains a singly link list to each
  ! item number (in order) so only item number less than whats in list needs
  ! to be compared.  This has room for improvement.
  DO WHILE ( ASSOCIATED(Recondat_Node % Recondat_Info) .AND.        &
             Recondat_Node % Recondat_Info % Sec_Item <               &
                               ancil_requests(irec)%stashcode )
    Recondat_Node => Recondat_Node % Next
  END DO

  IF (Recondat_Node % Recondat_Info % Sec_Item ==       &
                                  ancil_requests(irec)%stashcode ) THEN
    N_Levs = Recondat_Node % Recondat_Info % RLevs
    N_Pseudo_Levs = Recondat_Node % Recondat_Info % RPLevs
  ELSE
    ErrorStatus=1
    WRITE (CMessage, '(A, I7, A)') &
    'Attempted to process ancillary with STASHcode: ',   &
    ancil_requests(irec)%stashcode, ' however the prognostic will not be ' // &
    'present in the output dump. Either remove ancil request or check that '//&
    'science selected in namelists is compatible with this ancil field.'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  !       Ancillary data is nearly always read in on full grid.
  !       Special cases are as follows :

  SELECT CASE ( ancil_requests(irec)%section )

  CASE (stashcode_prog_sec)

    SELECT CASE( ancil_requests(irec)%item )

    CASE (stashcode_ozone)
      !If using Zonal Ozone then size is P_Rows not P_Field.
      IF ( zon_av_ozone ) THEN
        Len_anc_data = N_Levs * N_Pseudo_Levs * P_Rows
      ELSE
        Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field
      END IF

    CASE ( stashcode_riv_sequence,      &
      stashcode_riv_direction,     &
      stashcode_riv_storage )
      ! These fields are on an R_Field sized grid.
      Len_anc_data = N_Levs * N_Pseudo_Levs * R_Field

    CASE DEFAULT
      Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field

    END SELECT

  CASE DEFAULT
    Len_anc_data = N_Levs * N_Pseudo_Levs * P_Field

  END SELECT

  IF ( Len_anc_data == 0 ) THEN  !  Prognostic not in output dump
    ErrorStatus=10
    WRITE (CMessage, '(A, I7, A)') &
    ' Ancillary prognostic, Stash Code ',ancil_requests(irec)%stashcode, &
    ' not found in output dump.'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  Len_Ancil = Len_Ancil + Len_anc_data

  IF (PrintStatus > PrStatus_Diag .AND. mype == 0) THEN
    WRITE(umMessage,'(A)')'  '
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' section_code  ',ancil_requests(irec)%section
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' item_code     ',ancil_requests(irec)%item
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' n_levs        ',n_levs
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' n_pseudo_levs ',n_pseudo_levs
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' len_anc_data  ',len_anc_data
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
    WRITE(umMessage,'(A,I0)') ' len_ancil     ',len_ancil
    CALL umPrint(umMessage,src='rcf_calc_len_ancil_mod')
  END IF

END DO    ! irec

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_calc_len_ancil

END MODULE Rcf_calc_len_ancil_Mod
