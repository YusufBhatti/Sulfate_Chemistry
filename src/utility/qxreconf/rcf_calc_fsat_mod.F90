! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Set up Fsat/Fwetl for large-scale hydrology

MODULE Rcf_Calc_Fsat_Mod

! Subroutine Rcf_Calc_Fsat
!
! Description:
!   Initialises the surface saturation and/or wetland fraction
!
! Method:
!   The fields are calculated from the initial soil parameters and
!   soil moisture content and water table depth.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_FSAT_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_Fsat( fields_out, field_count_out, &
                          frac_pos, hdr_out )

!Use in relevant JULES routines
USE calc_baseflow_jls_mod, ONLY: calc_baseflow
USE calc_fsat_mod,         ONLY: calc_fsat

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_HeadAddress_Mod, ONLY: &
    IC_SoilMoistLevs

USE um_stashcode_mod, ONLY: &
    stashcode_clapp_hb,        &
    stashcode_ksat,            &
    stashcode_vol_smc_sat,     &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
    stashcode_ti_mean,         &
    stashcode_ti_sig,          &
    stashcode_fexp,            &
    stashcode_gamtot,          &
    stashcode_zw,              &
    stashcode_fsat,            &
    stashcode_fwetl,           &
    stashcode_prog_sec

USE umPrintMgr, ONLY:         &
    umPrint,                   &
    umMessage,                 &
    PrintStatus,               &
    PrStatus_Normal

USE nlstcall_mod, ONLY:   &
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

USE jules_hydrology_mod, ONLY: l_top

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
TYPE( um_header_type ), INTENT(IN) :: hdr_out

INTEGER, INTENT(IN)                :: field_count_out
INTEGER, INTENT(IN)                :: frac_pos

! Local variables
TYPE( field_type ), POINTER          :: frac_field
TYPE( field_type ), POINTER          :: clapp_hb
TYPE( field_type ), POINTER          :: ksat
TYPE( field_type ), POINTER          :: vol_smc_sat
TYPE( field_type ), POINTER          :: unfrozen_soil
TYPE( field_type ), POINTER          :: frozen_soil
TYPE( field_type ), POINTER          :: ti_mean
TYPE( field_type ), POINTER          :: ti_sig
TYPE( field_type ), POINTER          :: fexp
TYPE( field_type ), POINTER          :: gamtot
TYPE( field_type ), POINTER          :: zw

INTEGER                              :: soil_index  &
                                       (fields_out(frac_pos) % level_size)
INTEGER                              :: soil_pts
INTEGER                              :: SIZE
INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper
INTEGER                              :: j    ! looper
INTEGER                              :: n    ! looper

REAL          :: wutot    (fields_out(frac_pos) % level_size)
REAL          :: top_crit (fields_out(frac_pos) % level_size)
REAL          :: qbase    (fields_out(frac_pos) % level_size)

REAL          :: qbase_l (fields_out(frac_pos) % level_size,                   &
                               hdr_out % IntC( IC_SoilMoistLevs ) +1 )
REAL          :: Ksz     (fields_out(frac_pos) % level_size,                   &
                               0: hdr_out % IntC( IC_SoilMoistLevs ) )
REAL          :: zdepth  (0: hdr_out % IntC( IC_SoilMoistLevs) )

! Save the results if we call this routine again.
REAL, ALLOCATABLE, SAVE    :: fsat     (:)
REAL, ALLOCATABLE, SAVE    :: fwetl    (:)

INTEGER       :: ErrorStatus

CHARACTER (LEN=errormessagelength)    :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CALC_FSAT'

LOGICAL :: L_Gamtot=.FALSE.

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
  WRITE(umMessage,*) 'Setting up surface saturation/wetland fraction'
  CALL umPrint(umMessage,src='rcf_calc_fsat_mod')
END IF

IF (.NOT. l_top) THEN
  ErrorStatus = 10
  WRITE(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
  CALL Ereport ( RoutineName, ErrorStatus, CMessage)
END IF

! Setup pointer to output field.
frac_field => fields_out(frac_pos)

! Only calculate on first call when both fields are not allocated.
IF (.NOT. ALLOCATED(fsat) .AND. .NOT. ALLOCATED(fwetl)) THEN
  ! Allocate the SAVE arrays for results.
  ALLOCATE(fsat(frac_field % level_size))
  ALLOCATE(fwetl(frac_field % level_size))
  !----------------------------------------------------------------------
  ! Find required fields in output dump and read them in:
  !----------------------------------------------------------------------
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,           &
                   fields_out, field_count_out, pos)
  clapp_hb => fields_out(pos)
  CALL Rcf_Alloc_Field( clapp_hb )
  CALL Rcf_Read_Field( clapp_hb, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_ksat,               &
                   fields_out, field_count_out, pos)
  ksat => fields_out(pos)
  CALL Rcf_Alloc_Field( ksat )
  CALL Rcf_Read_Field( ksat, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,        &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  CALL Rcf_Alloc_Field( vol_smc_sat )
  CALL Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,      &
                   fields_out, field_count_out, pos)
  unfrozen_soil => fields_out(pos)
  CALL Rcf_Alloc_Field( unfrozen_soil )
  CALL Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,        &
                   fields_out, field_count_out, pos)
  frozen_soil => fields_out(pos)
  CALL Rcf_Alloc_Field( frozen_soil )
  CALL Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,            &
                   fields_out, field_count_out, pos)
  ti_mean => fields_out(pos)
  CALL Rcf_Alloc_Field( ti_mean )
  CALL Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,             &
                   fields_out, field_count_out, pos)
  ti_sig => fields_out(pos)
  CALL Rcf_Alloc_Field( ti_sig )
  CALL Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_fexp,               &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  CALL Rcf_Alloc_Field( fexp )
  CALL Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,             &
                   fields_out, field_count_out, pos)
  gamtot => fields_out(pos)
  CALL Rcf_Alloc_Field( gamtot )
  CALL Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_zw,                 &
                   fields_out, field_count_out, pos)
  zw => fields_out(pos)
  CALL Rcf_Alloc_Field( zw )
  CALL Rcf_Read_Field( zw, hdr_out, decomp_rcf_output )

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

  zdepth(0) = 0.0
  DO n = 1 , unfrozen_soil % levels
    zdepth(n) = zdepth(n-1) + Output_Grid % soil_depths(n)
  END DO

  !----------------------------------------------------------------------
  ! Set up saturated hydraulic conductivity for each soil layer:
  !----------------------------------------------------------------------
  DO j = 1 , soil_pts
    i = soil_index(j)
    DO n = 0 , unfrozen_soil % levels
      Ksz(i,n) = ksat % DATA(i,1)
    END DO
  END DO

  !----------------------------------------------------------------------
  ! Initialise variables to zero:
  !----------------------------------------------------------------------
  fsat(:) = 0
  fwetl(:) = 0
  qbase(:) = 0
  qbase_l(:,:) = 0
  frac_field % DATA (:,:) = 0

  !----------------------------------------------------------------------
  ! Setting wutot to 1.0 as no longer output from calc_baseflow
  !----------------------------------------------------------------------
  wutot(:) = 1.0
  !----------------------------------------------------------------------
  ! Calculate Baseflow:
  !----------------------------------------------------------------------
  CALL Calc_Baseflow(                    &
                  soil_pts,              &
                  soil_index,            &
                  ti_mean % level_size,  &
                  unfrozen_soil % levels,&
                  zdepth,                &
                  Ksz,                   &
                  clapp_hb % DATA,       &
                  fexp % DATA,           &
                  ti_mean % DATA,        &
                  zw % DATA,             &
                  frozen_soil % DATA,    &
                  top_crit,              &
                  qbase,                 &
                  qbase_l)

  !----------------------------------------------------------------------
  ! Calculate surface saturation and wetland fraction:
  !----------------------------------------------------------------------
  CALL Calc_Fsat(                        &
                  L_gamtot,              &
                  soil_pts,              &
                  soil_index,            &
                  ti_mean % level_size,  &
                  ti_mean % DATA,        &
                  ti_sig % DATA,         &
                  wutot,                 &
                  top_crit,              &
                  gamtot % DATA,         &
                  fsat,                  &
                  fwetl)


  !----------------------------------------------------------------------
  ! Tidy Up
  !----------------------------------------------------------------------
  CALL Rcf_Dealloc_Field( clapp_hb )
  CALL Rcf_Dealloc_Field( ksat )
  CALL Rcf_Dealloc_Field( vol_smc_sat )
  CALL Rcf_Dealloc_Field( unfrozen_soil )
  CALL Rcf_Dealloc_Field( frozen_soil )
  CALL Rcf_Dealloc_Field( ti_mean )
  CALL Rcf_Dealloc_Field( ti_sig )
  CALL Rcf_Dealloc_Field( fexp )
  CALL Rcf_Dealloc_Field( gamtot )
  CALL Rcf_Dealloc_Field( zw )

END IF

!----------------------------------------------------------------------
! Set up output to surface saturation or wetland fraction:
!----------------------------------------------------------------------
SELECT CASE(frac_field % stashmaster % item)
CASE (stashcode_fsat)
  frac_field % DATA(:,1)=fsat(:)
  DEALLOCATE(fsat)
CASE (stashcode_fwetl)
  frac_field % DATA(:,1)=fwetl(:)
  DEALLOCATE(fwetl)
END SELECT

IF (LTimer) CALL Timer( RoutineName, 4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_Fsat

END MODULE Rcf_Calc_Fsat_Mod
