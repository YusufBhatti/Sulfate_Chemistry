! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets up the integer constants for output dump

MODULE Rcf_Setup_IntC_Mod

!  Subroutine Rcf_Setup_IntC - initialisation of Integer Constants
!
! Description:
!   Sets up the output dump integer constants in the header.
!
! Method:
!   Uses namelist variables to define.
!   UMDP F3 defines the integer constants.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_INTC_MOD'

CONTAINS

SUBROUTINE Rcf_Setup_IntC( Hdr_In, Hdr_Out )

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Oper

USE Rcf_HeadAddress_Mod, ONLY: &
    IC_XLen,              IC_YLen,            IC_PLevels,       &
    IC_WetLevels,         IC_SoilTlevels,     IC_NoCloudLevels, &
    IC_NoSeaPts,          IC_TracerLevs,      IC_BLevels,       &
    IC_NumLandPoints,     IC_NumOzoneLevs,    IC_ConvectLevs,   &
    ic_mdi,               IC_SoilMoistLevs,   IC_1stConstRho,   &
    IC_TracerVars,        IC_HeightMethod,    IC_RiverRows,     &
    IC_RiverRowLength,    IC_Stochastic_flag, IC_stph_n1,       &
    IC_stph_n2,           IC_stph_seed

USE rcf_headers_Mod, ONLY: &
    inthd

USE river_inputs_mod, ONLY: &
    l_rivers

USE nlsizes_namelist_mod, ONLY: &
    tr_vars

USE missing_data_mod, ONLY: imdi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_header_type), INTENT(IN)   :: Hdr_In
TYPE (Um_header_type), TARGET       :: Hdr_Out

! Local variables
INTEGER, POINTER         :: IntC(:)
INTEGER                  :: i

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_INTC'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Start with clean sheet - IMDI all around
!------------------------------------------------------------------
IntC => Hdr_Out % IntC
IntC(:) = imdi

!------------------------------------------------------------------
! Set values we have numbers for
!------------------------------------------------------------------
IntC( IC_XLen        )   = Output_Grid % glob_p_row_length
IntC( IC_YLen        )   = Output_Grid % glob_p_rows
IntC( IC_PLevels     )   = Output_Grid % model_levels

IntC( IC_WetLevels   )   = Output_Grid % model_levels
IntC( IC_SoilTLevels )   = Output_Grid % st_levels
IntC( IC_NoCloudLevels ) = Output_Grid % cloud_levels
IntC( IC_TracerVars    ) = tr_vars
IntC( IC_BLevels       ) = Output_Grid % bl_levels
IntC( IC_NumOzoneLevs  ) = Output_Grid % ozone_levels
IntC( IC_SoilMoistLevs ) = Output_Grid % sm_levels
IntC( IC_1stConstRho   ) = Output_Grid % first_constant_r_rho_level
IntC( IC_ConvectLevs   ) = Output_Grid % conv_levels
IntC( IC_NumLandPoints ) = Output_Grid % glob_land_field
IntC( IC_HeightMethod  ) = Output_Grid % height_gen_method

IntC( IC_Stochastic_flag  ) = Hdr_In % IntC(IC_Stochastic_flag)
IntC( IC_stph_n1  )         = Hdr_In % IntC(IC_stph_n1)
IntC( IC_stph_n2  )         = Hdr_In % IntC(IC_stph_n2)
IntC( IC_stph_seed  )       = Hdr_In % IntC(IC_stph_seed)

IF ( l_rivers ) THEN
  IntC( IC_RiverRows     )  = Output_Grid % glob_r_rows
  IntC( IC_RiverRowLength)  = Output_Grid % glob_r_row_length
ELSE
  IntC( IC_RiverRows     )  = imdi
  IntC( IC_RiverRowLength)  = imdi
END IF

IntC( IC_TracerLevs    ) = Output_Grid % tr_levels
IntC( ic_mdi           ) = imdi

! Note AMIP/other IntC (35 +) isn't set yet.
! Note also that IntC(15) - number of land tile types - has not
! been set

!------------------------------------------------------------------
! Module overrides
!------------------------------------------------------------------
DO i = 1, Hdr_Out % LenIntC
  IF ( IntHd(i) /= imdi ) THEN
    IF ( PrintStatus >= PrStatus_Oper ) THEN
      WRITE(umMessage,*) 'IntC(',i,') has been reset from ', IntC(i), &
                  ' to ', IntHd(i)
      CALL umPrint(umMessage,src='rcf_setup_intc_mod')
    END IF

    IntC( i ) = IntHd( i )
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_IntC
END MODULE Rcf_Setup_IntC_Mod
