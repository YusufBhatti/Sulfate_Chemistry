! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates RHO for the output dump

MODULE Rcf_Calc_Rho_Mod

!  Subroutine Rcf_Calc_Rho - calculates tho
!
! Description:
!   Calculates RHO for the 5.0/5.1 dump (rho should not be interpolated
!   so as to maintain dynamical balance)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_RHO_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_Rho( theta, q, exner, p, theta_heights,        &
                         rho_heights, rho )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE Rcf_Calc_Exner_Theta_Mod, ONLY: &
    Rcf_Calc_Exner_Theta

USE planet_constants_mod, ONLY: r, repsilon

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(IN)     :: theta
TYPE( field_type ), INTENT(IN)     :: q
TYPE( field_type ), INTENT(IN)     :: p           ! on rho levels
TYPE( field_type ), INTENT(IN)     :: exner       ! on rho levels
TYPE( field_type ), INTENT(INOUT)  :: rho

REAL, INTENT(IN)                   :: rho_heights( :, 0: )
REAL, INTENT(IN)                   :: theta_heights( :, 0: )

! Local variables
TYPE( field_type )                 :: exner_theta ! on theta levels
INTEGER                            :: i
INTEGER                            :: k
INTEGER                            :: k_off
INTEGER                            :: k_rho
REAL                               :: weight1
REAL                               :: weight2
REAL                               :: weight3
REAL                               :: temp
REAL                               :: work_real ( theta % level_size )
REAL                               :: work_real2( theta % level_size )

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CALC_RHO'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Need to calculate exner on theta levels
!-------------------------------------------------------------------
CALL Rcf_field_equals( exner_theta, exner )
exner_theta % levels = theta % levels
CALL Rcf_Alloc_Field( exner_theta )

IF ( theta % bottom_level == 0 ) THEN
  k_off = 1
  ! Have to handle zeroth level theta differently
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels-1, &
                             theta_heights, rho_heights(:,1:),          &
                             exner % DATA, exner_theta % DATA(:,2:) )
  exner_theta % DATA (:,1) = exner_theta % DATA (:,2)
ELSE
  k_off = 0
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                             theta_heights, rho_heights(:,1:),        &
                             exner % DATA, exner_theta % DATA )
END IF

!--------------------------------------------------------------------
! Now do the calculation of rho
!--------------------------------------------------------------------

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*) 'Calculating Rho'
  CALL umPrint(umMessage,src='rcf_calc_rho_mod')
END IF

k = 1+k_off
DO i = 1, theta % level_size
  k_rho = k-k_off
  ! calculate thetav
  work_real(i) = theta % DATA(i,k) * (1.0 +                             &
                                   (1.0/repsilon -1.0) * q % DATA(i,k) )
  ! calculate rho
  rho % DATA(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho) * &
                 p % DATA(i,k_rho) / (r * work_real(i) * exner % DATA(i,k_rho))
END DO

DO k = 2+k_off, theta % levels
  k_rho = k-k_off
  DO i = 1, theta % level_size
    work_real2(i) = work_real(i)
  END DO

  IF (k  <=  q % levels ) THEN
    DO i = 1, theta % level_size
      work_real(i) = theta % DATA(i,k) * (1.0 +                           &
                                     (1.0/repsilon -1.0) * q % DATA(i,k) )
    END DO

  ELSE
    DO i = 1, theta % level_size
      work_real(i) = theta % DATA(i,k)
    END DO
  END IF

  IF (k  /=  theta % levels) THEN
    DO i = 1, theta % level_size
      weight1 = rho_heights(i,k_rho) - theta_heights(i,k_rho-1)
      weight2 = theta_heights(i,k_rho) - rho_heights(i,k_rho)
      weight3 = theta_heights(i,k_rho) - theta_heights(i,k_rho-1)

      !      temp = ( weight2 * work_real(i) +             &
      !               weight1 * work_real2(i) ) /          &
      !               weight3

      temp = ( weight1 * work_real(i) +             &
               weight2 * work_real2(i) ) /          &
               weight3

      rho % DATA(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho)   &
                      * p % DATA(i,k_rho) / (r * temp * exner % DATA(i,k_rho))
    END DO

  ELSE

    DO i= 1, theta % level_size

      temp = work_real2(i) *                                   &
             exner_theta % DATA(i, exner_theta % levels - 1) / &
             exner % DATA(i,k_rho)

      rho % DATA(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho) *      &
                        p % DATA(i,k_rho) / (r * temp * exner % DATA(i,k_rho) )

    END DO

  END IF
END DO

!------------------------------------------------------------------
! Tidy up
!------------------------------------------------------------------
CALL Rcf_DeAlloc_Field( exner_theta )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_Rho
END MODULE Rcf_Calc_Rho_Mod
