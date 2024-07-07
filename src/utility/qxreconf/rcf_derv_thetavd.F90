! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_derv_thetavd_mod

IMPLICIT NONE
!  Subroutine Rcf_Derv_Thetavd_Mod

! Description:
!   Derive virtual dry potential temperature

! Method:
!   As previously implement in VN8.4 atm_step_4A.F90 l. 1101-1107
!   this routine computed theta_vd from theta and
!   mixing ratio moisture (m_v)

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_THETAVD_MOD'

CONTAINS

SUBROUTINE rcf_derv_thetavd( fields_in, field_count_in,    &
                             fields_out, field_count_out,  &
                             hdr_in,hdr_out,               &
                             thetavd )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_mv,  stashcode_theta, stashcode_prog_sec

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parcore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE rcf_grid_type_mod, ONLY:                    &
    grid_type,                                   &
    input_grid, output_grid

USE planet_constants_mod, ONLY: recip_epsilon

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( field_type ), POINTER       :: fields_in(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( um_header_type), INTENT(IN) :: hdr_in
TYPE( field_type ), INTENT(INOUT), TARGET :: thetavd
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_count_in

! Internal variables
TYPE( field_type ), POINTER       :: theta_out
TYPE( field_type ), POINTER       :: mv

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_DERV_THETAVD'

INTEGER                           :: pos   ! position in array
INTEGER                           :: dummy_pos

INTEGER                           :: i

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 'Deriving theta_vd '
  CALL umPrint(umMessage,src='rcf_derv_thetavd')
END IF



CALL rcf_locate(stashcode_prog_sec, stashcode_theta, &
                  fields_out, field_count_out, pos)
theta_out => fields_out(pos)
CALL rcf_alloc_field( theta_out )
CALL rcf_read_field(  theta_out, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mv, &
                  fields_out, field_count_out, pos)
mv => fields_out(pos)
CALL rcf_alloc_field( mv )
CALL rcf_read_field( mv, hdr_out, decomp_rcf_output )

DO i = 1, thetavd % levels
  thetavd % DATA(:,i) = theta_out % DATA(:,i) * &
                    (1.0 + mv % DATA(:,i)*recip_epsilon)
END DO

CALL rcf_dealloc_field( mv )
CALL rcf_dealloc_field( theta_out )

! fix 0th level:

IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE(umMessage,FMT='(a)') 'Thetavd level 0 == level 1'
  CALL umPrint(umMessage,src='rcf_derv_thetavd')
END IF

thetavd % DATA(:,1) =  thetavd % DATA(:,2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_derv_thetavd
END MODULE rcf_derv_thetavd_mod
