! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Estimate deep soil moisture fraction for large-scale hydrology

MODULE Rcf_Est_Sthzw_Mod

! Subroutine Rcf_Est_Sthzw
!
! Description:
!   Estimate deep soil moisture fraction for large-scale hydrology when
!   not available in the dump.
!
! Method:
!   Estimates a temporary gridbox mean water table depth by assuming
!   that the soil drainage into the additional deep layer and
!   the baseflow out of this layer balance, i.e:
!   Ksat*theta^(2b+3)=Ksat/f*exp(-f*Z-Ti)
!   This is then used to estimate the soil moisture fraction in the
!   additional deep layer in the large-scale hydrology (LTOP) scheme.
!   When the estimated temporary mean water table is above the deep
!   layer assume that the deep soil moisture fraction is 1.0
!   When it is in the deep layer assume that is linearly dependent
!   on the water table position within the layer.
!   This does not have to be very accurate, its purpose is to reduce
!   spin-up time when the initial fields are not available.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_EST_STHZW_MOD'

CONTAINS

SUBROUTINE Rcf_Est_Sthzw( fields_out, field_count_out, &
                             sthzw, hdr_out )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_clapp_hb,        &
    stashcode_vol_smc_sat,     &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
    stashcode_sthzw,           &
    stashcode_ti_mean,         &
    stashcode_fexp,            &
    stashcode_prog_sec

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE nlstcall_mod, ONLY: &
    LTimer


USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Ereport_mod, ONLY: &
    Ereport

! Replaces c_topog.h
USE jules_hydrology_mod, ONLY: l_top, zw_max

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)  :: sthzw
TYPE( um_header_type ), INTENT(IN) :: hdr_out

INTEGER, INTENT(IN)                :: field_count_out

! Local variables
TYPE( field_type ), POINTER          :: clapp_hb
TYPE( field_type ), POINTER          :: vol_smc_sat
TYPE( field_type ), POINTER          :: unfrozen_soil
TYPE( field_type ), POINTER          :: frozen_soil
TYPE( field_type ), POINTER          :: ti_mean
TYPE( field_type ), POINTER          :: fexp

INTEGER                              :: soil_index  &
                                       (sthzw % level_size)
INTEGER                              :: soil_pts

INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper
INTEGER                              :: n    ! looper
INTEGER                              :: j    ! looper

INTEGER      :: ErrorStatus

CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_EST_STHZW'

REAL                                 :: fnstuf
REAL                                 :: temp_fn
REAL                                 :: dumzw
REAL                                 :: soil_depth

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (LTimer) CALL Timer( RoutineName, 3)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Estimating the initial deep layer water fraction'
  CALL umPrint(umMessage,src='rcf_est_sthzw_mod')
END IF

IF (.NOT. l_top) THEN
  ErrorStatus = 10
  WRITE(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
  CALL Ereport ( RoutineName, ErrorStatus, CMessage)
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------

CALL Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,             &
                 fields_out, field_count_out, pos)
clapp_hb => fields_out(pos)
CALL Rcf_Alloc_Field( clapp_hb )
CALL Rcf_Read_Field( clapp_hb, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,          &
                 fields_out, field_count_out, pos)
vol_smc_sat => fields_out(pos)
CALL Rcf_Alloc_Field( vol_smc_sat )
CALL Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,        &
                 fields_out, field_count_out, pos)
unfrozen_soil => fields_out(pos)
CALL Rcf_Alloc_Field( unfrozen_soil )
CALL Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,          &
                 fields_out, field_count_out, pos)
frozen_soil => fields_out(pos)
CALL Rcf_Alloc_Field( frozen_soil )
CALL Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,              &
                 fields_out, field_count_out, pos)
ti_mean => fields_out(pos)
CALL Rcf_Alloc_Field( ti_mean )
CALL Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_fexp,                 &
                 fields_out, field_count_out, pos)
fexp => fields_out(pos)
CALL Rcf_Alloc_Field( fexp )
CALL Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Set up soil index:
!----------------------------------------------------------------------
soil_pts=0
soil_index(:)=0
DO i=1 , vol_smc_sat % level_size
  IF (vol_smc_sat % DATA(i,1) > 0.0) THEN
    soil_pts = soil_pts + 1
    soil_index(soil_pts) = i
  END IF
END DO

Soil_depth=0.0
DO n = 1,unfrozen_soil % levels
  soil_depth=soil_depth+Output_Grid % soil_depths(n)
END DO

!----------------------------------------------------------------------
! Estimate sthzw by assuming equilibruim:
!----------------------------------------------------------------------
sthzw % DATA(:,:) = 1.0

DO j = 1, soil_pts
  i = soil_index(j)

  IF (fexp % DATA(i,1) <= 0.0) THEN
    ErrorStatus = 20
    WRITE (CMessage, '(A)') 'Exp. decay in Vol_Smc_Sat wrongly/not set'
    CALL Ereport ( RoutineName, ErrorStatus, CMessage)
  END IF

  IF (frozen_soil % DATA(i,unfrozen_soil % levels) < 1.0) THEN

    fnstuf = unfrozen_soil % DATA(i,unfrozen_soil % levels) /  &
      (1.0 - frozen_soil % DATA(i,unfrozen_soil % levels))
    temp_fn = fnstuf**(2.0 * clapp_hb % DATA(i,1) + 3.0) *      &
      fexp % DATA(i,1) * EXP(ti_mean % DATA(i,1)) +            &
      EXP(-fexp % DATA(i,1) * (Zw_Max - Soil_Depth))

    IF (temp_fn > 0.0) THEN
      dumzw =   &
      -1.0/Fexp % DATA(i,1) * LOG(temp_fn) + Soil_Depth
    END IF

    IF (dumzw > Zw_Max) THEN
      sthzw % DATA(i,1)=0.0
    END IF
    IF (dumzw < 0.0) THEN
      Dumzw = 1.0
    END IF

    IF (dumzw > soil_depth) THEN
      sthzw % DATA(i,1)=(Zw_Max-dumzw)/(Zw_Max-soil_depth)
    END IF

  END IF

END DO

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( clapp_hb )
CALL Rcf_Dealloc_Field( vol_smc_sat )
CALL Rcf_Dealloc_Field( unfrozen_soil )
CALL Rcf_Dealloc_Field( frozen_soil )
CALL Rcf_Dealloc_Field( ti_mean )
CALL Rcf_Dealloc_Field( fexp )

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Est_Sthzw

END MODULE Rcf_Est_Sthzw_Mod
