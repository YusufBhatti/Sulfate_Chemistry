! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------

! Checks on condensate at the beginning of convection if PC2 being used.

MODULE conv_pc2_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CONV_PC2_INIT_MOD'
CONTAINS

SUBROUTINE  conv_pc2_init(rows, row_length, exner_theta_levels,            &
                          theta_conv, q_conv, qcl_conv, qcf_conv,          &
                          cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

USE atm_fields_bounds_mod, ONLY: tdims_s
USE planet_constants_mod,   ONLY: cp
USE water_constants_mod,   ONLY: lc, lf
USE nlsizes_namelist_mod,  ONLY: model_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Updates convection diagnostics after call to glue each substep of
!   convection.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Convection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN)  ::    &
  row_length               & ! Row length
 ,rows                       ! Number of rows

REAL, INTENT(IN)  ::                                 &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                           &
  theta_conv(row_length, rows, model_levels)     & ! theta (K)
 ,q_conv(row_length, rows, model_levels)         & ! water vapour  (kq/kg)
 ,qcl_conv(row_length, rows, model_levels)       & ! cloud liquid water (kq/kg)
 ,qcf_conv(row_length, rows, model_levels)       & ! cloud ice water (kq/kg)
 ,cf_liquid_conv(row_length, rows, model_levels) & ! liquid cloud fraction
 ,cf_frozen_conv(row_length, rows, model_levels) & ! ice cloud fraction
 ,bulk_cf_conv(row_length, rows, model_levels)     ! total cloud fraction

! ------------------------------------------------------------------------------
! Local declarations:
! ------------------------------------------------------------------------------
INTEGER  ::         &
  i,j,k                ! loop counters

INTEGER  ::         &
  n_qcx_1           &  ! Record  out-of-range input events
  ,n_qcx_2             ! Record  out-of-range input events


! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CONV_PC2_INIT'

!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Prevent negative condensate problems by condensing vapour if needed.
! Input fields are updated without updating increments so change will
! only affect PC2 increments and not the top-level QX_STAR fields.
! ----------------------------------------------------------------------

n_qcx_1 = 0
n_qcx_2 = 0

!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k)                     &
!$OMP& SHARED(rows,row_length,model_levels,exner_theta_levels,        &
!$OMP&        theta_conv,q_conv,qcl_conv,qcf_conv,cf_liquid_conv,     &
!$OMP&        cf_frozen_conv,bulk_cf_conv,cp)                         &
!$OMP& SCHEDULE(STATIC) REDUCTION(+:n_qcx_1,n_qcx_2)
DO k=1, model_levels

  DO j=1, rows
    DO i=1, row_length
      IF (qcf_conv(i,j,k)  <   0.0) THEN
        !             Freeze some liquid to zero ice content
        qcl_conv(i,j,k) = qcl_conv(i,j,k) + qcf_conv(i,j,k)
        theta_conv(i,j,k) = theta_conv(i,j,k) -                       &
                ( (qcf_conv(i,j,k) * lf) / (cp * exner_theta_levels(i,j,k)) )
        qcf_conv(i,j,k) = 0.0
        cf_frozen_conv(i,j,k) = 0.0
        bulk_cf_conv(i,j,k) = cf_liquid_conv(i,j,k)
        n_qcx_1 = n_qcx_1 + 1
      END IF

      IF (qcl_conv(i,j,k)  <   0.0) THEN
        !             Condense some vapour to zero liquid content
        q_conv(i,j,k) = q_conv(i,j,k) + qcl_conv(i,j,k)
        theta_conv(i,j,k) = theta_conv(i,j,k) -                 &
                 ( (qcl_conv(i,j,k) * lc) / (cp * exner_theta_levels(i,j,k)) )
        qcl_conv(i,j,k) = 0.0
        cf_liquid_conv(i,j,k) = 0.0
        bulk_cf_conv(i,j,k) = cf_frozen_conv(i,j,k)
        n_qcx_2 = n_qcx_2 + 1
      END IF

  !     Might also be necessary to place a limit on supersaturation.
  !     Convection copes but cloud scheme response is less predictable.
  !     Would need eg. LS_CLD_C to reduce supersaturation consistently.

      !           Ensure that input values of cloud fields lie within the
      !           bounds of physical possibility (should do nothing).

      bulk_cf_conv(i,j,k)   = MAX(0.0, ( MIN(1.0, bulk_cf_conv(i,j,k)) ))
      cf_liquid_conv(i,j,k) = MAX(0.0, ( MIN(1.0, cf_liquid_conv(i,j,k)) ))
      cf_frozen_conv(i,j,k) = MAX(0.0, ( MIN(1.0, cf_frozen_conv(i,j,k)) ))
    END DO  ! I
  END DO  ! J

END DO  ! K
!$OMP END PARALLEL DO

IF ( PrintStatus >= PrStatus_Normal ) THEN

  IF (n_qcx_1  >   0) THEN
    WRITE(umMessage,'(A,I07,A)')                                            &
      'Qcf < 0 fixed by PC2, (', n_qcx_1, ' occurances)'
    CALL umPrint(umMessage,src='conv_pc2_init')
  END IF
  IF (n_qcx_2  >   0) THEN
    WRITE(umMessage,'(A,I07,A)')                                            &
      'Qcl < 0 fixed by PC2, (', n_qcx_2, ' occurances)'
    CALL umPrint(umMessage,src='conv_pc2_init')
  END IF

END IF

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE
END MODULE conv_pc2_init_mod
