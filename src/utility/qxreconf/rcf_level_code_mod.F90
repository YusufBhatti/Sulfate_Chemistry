! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Decode the STASH level code

MODULE Rcf_Level_Code_Mod

IMPLICIT NONE

! Description:
!   Sets ILOUT to an appropriate level size according to the value of
!   ILIN.
!
! *****************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to levcod.F90 and
! *****************************************************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_LEVEL_CODE_MOD'

CONTAINS

SUBROUTINE Rcf_Level_Code( ilin, ilout, Grid )


USE stash_model_mod, ONLY: stlevgwdrag, botvdifflev, topvdifflev


USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

! The following are retired from the STASHmaster.  Using pseudo levels instead.
USE rad_input_mod, ONLY:    &
    H_SWBands,              &
    H_LWBands

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_HeadAddress_Mod, ONLY: &
    FH_GridStagger_C,           &
    FH_GridStagger_Endgame

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)    :: ilin           ! Model level code
INTEGER, INTENT(OUT)   :: ilout          ! The actual level

TYPE (Grid_Type), INTENT(IN) :: Grid     ! The grid that decoded
                                         ! values should correspond to


! Local variables
INTEGER                      :: ErrorStatus
CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_LEVEL_CODE'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE ( ilin )
CASE ( 1 )                            ! First atmos level
  ilout=1

CASE ( 2 )                            ! Top atmos level
  ilout= Grid % model_levels

CASE ( 3 )                            ! Top wet level
  ! Wet levels has been hardwired to match model levels
  ilout= Grid % model_levels
CASE ( 4 )
  ilout= Grid % model_levels - 1

CASE ( 5 )                            ! First boundary layer level
  ilout=MIN(1,Grid % bl_levels)

CASE ( 6 )                            ! Last boundary layer level
  ilout=Grid % bl_levels

CASE ( 7 )
  ilout= Grid % bl_levels+1

CASE ( 8 )                            ! First soil level
  ilout=MIN(1, Grid % st_levels)

CASE ( 9 )                            ! Last soil level
  ilout= Grid % st_levels

CASE ( 10 )                           ! First tracer level
  ilout= Grid % model_levels - Grid % tr_levels+1

CASE ( 11 )                           ! Last tracer level
  ilout= Grid % model_levels

CASE ( 12 )
  ilout= Grid % model_levels+1

CASE ( 13 )                           ! First gravity wave drag level
  ilout=StLevGWdrag

CASE ( 14 )                           ! Last gravity wave drag level
  ilout= Grid % model_levels

CASE ( 15 )
  ilout=BotVDiffLev

CASE ( 16 )
  ilout=TopVDiffLev-1

CASE ( 17 )
  ilout=TopVDiffLev

CASE ( 18 )
  ilout= Grid % bl_levels-1

CASE ( 19 )
  ilout= Grid % model_levels+1

CASE ( 20 )
  ilout=MIN(2, Grid % st_levels)

  ! Ocean removed at vn7.0 so this is redundant
CASE ( 21 )
  ilout=1

  ! Ocean removed at vn7.0 so this is redundant
CASE ( 22 )
  ilout=0

CASE ( 23 )
  ilout= Grid % ozone_levels

CASE ( 24 )
  ilout= Grid % model_levels*h_swbands

CASE ( 25 )
  ilout=( Grid % model_levels+1)*h_swbands

CASE ( 26 )
  ilout= Grid % model_levels*h_swbands

CASE ( 27 )
  ilout= Grid % model_levels*h_lwbands

CASE ( 28 )
  ilout=( Grid % model_levels+1)*h_lwbands

CASE ( 29 )
  ilout= Grid % model_levels*h_lwbands

CASE ( 30 )
  ilout=2

CASE ( 32 )
  ilout=h_swbands

CASE ( 33 )
  ilout=h_lwbands

CASE ( 34 )
  ilout= Grid % sm_levels

CASE ( 35 )
  ilout= Grid % cloud_levels

CASE ( 36 )                       ! Wave model first level (direction)
  ilout=1

CASE ( 37 )                       ! Wave model last level (direction)
                                  ! No wave model in rcf
  !    ILOUT=NANG
  ilout=0

CASE ( 38 )                       ! Surace theta level
  ilout=0

CASE ( 39 )
  ! Number of ISCCP simulator levels
  ilout=7

CASE ( 40 )
  ! Fields which are different between ENDGAME and New Dynamics regarding the
  ! whether it starts at theta level 0 in ENDGAME but 1 in New Dynamics.  Use
  ! grid stagger value to identify endgame grid - not perfect.
  IF (grid % grid_stagger == FH_GridStagger_Endgame) THEN
    ilout=0
  ELSE
    ilout=1
  END IF

CASE DEFAULT
  WRITE(Cmessage, '(A, I5)') 'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND: ',ilin
  ErrorStatus=1
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Level_Code
END MODULE Rcf_Level_Code_Mod

