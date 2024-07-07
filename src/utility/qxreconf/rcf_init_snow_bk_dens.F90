! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialises the bulk density of the snow pack.

MODULE rcf_init_snow_bk_dens_mod
IMPLICIT NONE

! Description:
!     This subroutine initializes the bulk density of the snow pack.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_SNOW_BK_DENS_MOD'

CONTAINS

SUBROUTINE rcf_init_snow_bk_dens( fields_in, field_count_in, hdr_in,    &
                                  fields_out, field_count_out, hdr_out, &
                                  bulkdensity )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE um_stashcode_mod, ONLY:    &
    stashcode_prog_sec          , &
    stashcode_snowdep_grd_tile  , &
    stashcode_snowpack_bk_dens

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parcore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE jules_snow_mod, ONLY: nsmax, rho_snow_const, rho_snow_fresh

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_in(:), fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_in, hdr_out
INTEGER, INTENT(IN)               :: field_count_in, field_count_out
TYPE( field_type ), INTENT(INOUT) :: bulkdensity

! Internal variables
TYPE( field_type ), POINTER  :: snowdepth

INTEGER      ::  pos         ! position in array

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_INIT_SNOW_BK_DENS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE(umMessage,'(a,i4)') &
      ' Initialising the bulk density of the snow pack.'
  CALL umPrint(umMessage,src='rcf_init_snow_bk_dens')
END IF

!-----------------------------------------------------------------------------
! The depth of snow on the ground will have been set by now and
! can be used to see where we have snow. This is simpler than using
! the tile and ground snow stores where we'd need to know about the
! canopy model.
!-----------------------------------------------------------------------------

CALL rcf_locate( stashcode_prog_sec, stashcode_snowdep_grd_tile, &
           fields_out, field_count_out, pos)
snowdepth => fields_out(pos)
CALL rcf_alloc_field( snowdepth )
CALL rcf_read_field( snowdepth, hdr_out, decomp_rcf_output )

bulkdensity % DATA(:, :) = rho_snow_const

! If using the multilayer scheme, set the default where there is
! no snow to the value for fresh snow.
IF ( nsmax > 0 ) THEN
  WHERE ( snowdepth % DATA < EPSILON(snowdepth % DATA) )
    bulkdensity % DATA(:, :) = rho_snow_fresh
  END WHERE
END IF

!----------------------------------------------------------------------
! Clear up dynamic memory used.
!----------------------------------------------------------------------
CALL rcf_dealloc_field( snowdepth )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_init_snow_bk_dens
END MODULE rcf_init_snow_bk_dens_mod
