! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculations between T and Theta

MODULE Rcf_theta_t_Convs_mod

!  Subroutine Rcf_Conv_Theta_T - converts theta to t
!  Subroutine Rcf_Conv_T_Theta - converts t to theta
!
! Description:
!   Performs calculations between T and Theta
!
! Method:
!   Derived from New Dynamics 2.7 code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE submodel_mod, ONLY: atmos_im

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_THETA_T_CONVS_MOD'

CONTAINS

! *****************************************************************
! Routine to calculate T from Theta
! *****************************************************************
SUBROUTINE Rcf_Conv_theta_t( theta, fields, field_count, hdr, decomp, &
                             rho_heights, theta_heights)

USE Rcf_field_equals_mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Calc_Exner_Theta_Mod, ONLY: &
    Rcf_Calc_Exner_Theta

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_t,               &
    stashcode_theta,           &
    stashcode_prog_sec

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)  ! array of fields
TYPE( field_type ), INTENT(INOUT)  :: theta      ! theta field - convert
TYPE( um_header_type ), INTENT(IN) :: hdr        ! header of dump
INTEGER, INTENT(IN)                :: decomp     ! working decomposition
INTEGER, INTENT(IN)                :: field_count ! # of fields in dump
REAL, INTENT(IN)                   :: theta_heights( :, 0 : )
REAL, INTENT(IN)                   :: rho_heights( :, 0 : )

! Local variables
INTEGER                            :: i
INTEGER                            :: j
INTEGER                            :: pos
TYPE( field_type )                 :: exner_theta
TYPE( field_type )                 :: exner

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_CONV_THETA_T'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Get exner and convert it to P on rho levels
! (for VAR we will only have P - just read it in)
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields, field_count, pos,  zero_ok_arg = .TRUE. )

IF (pos /= 0) THEN
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )

ELSE
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields, field_count, pos )
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )

  CALL Rcf_Conv_P_Exner( exner )
END IF

!-------------------------------------------------------------------
! Convert exner from rho levels to theta levels
!-------------------------------------------------------------------
CALL Rcf_Field_Equals( exner_theta, exner )
exner_theta % levels = theta % levels
CALL Rcf_Alloc_Field( exner_theta )

IF ( theta % bottom_level == 0 ) THEN
  ! Have to handle zeroth level theta differently
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels-1, &
                             theta_heights, rho_heights(:,1:),           &
                             exner % DATA, exner_theta % DATA(:,2:) )
  exner_theta % DATA (:,1) = exner_theta % DATA (:,2)
ELSE
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                             theta_heights, rho_heights(:,1:),         &
                             exner % DATA, exner_theta % DATA )
END IF


!-------------------------------------------------------------------
! Calculate T from theta and exner
!-------------------------------------------------------------------
DO j = 1, theta % levels
  DO i = 1, theta % level_size
    theta % DATA(i,j) = theta % DATA(i,j) * exner_theta % DATA(i,j)
  END DO
END DO

!-------------------------------------------------------------------
! Tidy up
!-------------------------------------------------------------------
CALL Rcf_Dealloc_Field( exner_theta )
CALL Rcf_Dealloc_Field( exner )

! Reset stashmaster to point to a more appropriate one.
theta % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                  MOD(stashcode_t,1000) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Conv_theta_t





! *****************************************************************
! Routine to calculate theta from T
! *****************************************************************
SUBROUTINE Rcf_Conv_t_theta( theta, fields, field_count, hdr, decomp, &
                             rho_heights, theta_heights )

USE Rcf_field_equals_mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx


USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Calc_Exner_Theta_Mod, ONLY: &
    Rcf_Calc_exner_Theta

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_t,               &
    stashcode_theta,           &
    stashcode_prog_sec

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields(:)  ! array of fields
TYPE( field_type ), INTENT(INOUT)  :: theta      ! theta field - convert
TYPE( um_header_type ), INTENT(IN) :: hdr        ! header of dump
INTEGER, INTENT(IN)                :: decomp     ! working decomposition
INTEGER, INTENT(IN)                :: field_count ! # of fields in dump
REAL, INTENT(IN)                   :: theta_heights( : , 0: )
REAL, INTENT(IN)                   :: rho_heights(   : , 0: )

! Local variables
INTEGER                            :: i
INTEGER                            :: j
INTEGER                            :: pos
TYPE( field_type )                 :: exner_theta
TYPE( field_type )                 :: exner
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_CONV_T_THETA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Get exner and convert it to P on rho levels
! (for VAR we will only have P - just read it in)
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_exner, &
                 fields, field_count, pos,  zero_ok_arg = .TRUE. )

IF (pos /= 0) THEN
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )

ELSE
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_p, &
                   fields, field_count, pos )
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )

  CALL Rcf_Conv_P_Exner( exner )
END IF

!-------------------------------------------------------------------
! Convert exner from rho levels to theta levels
!-------------------------------------------------------------------
CALL Rcf_Field_Equals( exner_theta, exner )
exner_theta % levels = theta % levels
CALL Rcf_Alloc_Field( exner_theta )

IF ( theta % bottom_level == 0 ) THEN
  ! Have to handle zeroth level theta differently
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels-1, &
                             theta_heights, rho_heights(:,1:),           &
                             exner % DATA, exner_theta % DATA(:,2:) )
  exner_theta % DATA (:,1) = exner_theta % DATA (:,2)
ELSE
  CALL Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                             theta_heights, rho_heights(:,1:),         &
                             exner % DATA, exner_theta % DATA )
END IF

!-------------------------------------------------------------------
! Calculate theta from T and exner
!-------------------------------------------------------------------
DO j = 1, theta % levels
  DO i = 1, theta % level_size
    theta % DATA(i,j) = theta % DATA(i,j) / exner_theta % DATA(i,j)
  END DO
END DO

!-------------------------------------------------------------------
! Tidy up
!-------------------------------------------------------------------
CALL Rcf_Dealloc_Field( exner_theta )
CALL Rcf_Dealloc_Field( exner )

! Reset stashmaster to point to a more appropriate one.
theta % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                  stashcode_theta )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Conv_t_theta

END MODULE Rcf_theta_t_Convs_mod

