! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates the heights in the SCM.
!
! Subroutine Interface:

SUBROUTINE calc_levels                                                      &
! Input data
  ( orog, height_gen_method                                                 &
  , bl_levels, model_levels, rows, row_length )

USE vertnamelist_mod, ONLY:                                                 &
    z_top_of_model,  first_constant_r_rho_level

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels,                  &
                           eta_theta_levels, eta_rho_levels
USE bl_option_mod, ONLY: z_nl_bl_levels, nl_bl_levels,                      &
                         i_bl_vn, i_bl_vn_9b, i_bl_vn_9c, off, alpha_cd_in, &
                         alpha_cd

USE stochastic_physics_run_mod, ONLY: zmin_pert_theta, minlev_pert_theta,   &
                                      zmax_pert_theta, maxlev_pert_theta,   &
                                      i_pert_theta

USE planet_constants_mod, ONLY: planet_radius

USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY: umPrint, printstatus, prstatus_normal

IMPLICIT NONE

!
! Description:  To set up r_theta_levels and r_rho_levels for the
!               Single Column Model
!
!
! Method:
!   It follows the code used in the main UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!

! Code description:
!      Fortran 90, Written to UM coding standards
!      as specified in UMDP 3, vn8.2

! Inputs
INTEGER ::                                                        &
  bl_levels                                                       &
, model_levels                                                    &
, rows                                                            &
, row_length                                                      &
, height_gen_method

REAL ::                                                           &
  orog(row_length, rows)

CHARACTER(LEN=80) ::                                                 &
       Cmessage              ! Error message if ICODE >0

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'Calc_Levels'

INTEGER :: ErrorStatus

! Local variables
REAL ::                                                           &
  r_ref_theta(model_levels)                                       &
, r_ref_rho(model_levels)

INTEGER ::                                                        &
  i,j,k

INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

!----------------------------------------------------------------------
!     Set up heights
!----------------------------------------------------------------------

! Set reference profile

DO k=1, model_levels
  r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
  r_ref_rho(k)   = eta_rho_levels(k)   * z_top_of_model
END DO

! Set bottom level, ie orography
DO j=1, rows
  DO i=1, row_length
    r_theta_levels(i,j,0) = orog(i,j) + planet_radius
  END DO
END DO
! For constant levels set r to be a constant on the level
DO k=first_constant_r_rho_level, model_levels
  DO j=1, rows
    DO i=1, row_length
      r_theta_levels(i,j,k) = planet_radius + r_ref_theta(k)
      r_rho_levels(i,j,k)   = planet_radius + r_ref_rho(k)
    END DO
  END DO
END DO

SELECT CASE( height_gen_method)
CASE ( height_gen_original)
  ! The original version of height generation used in the SI dynamics
  !
  ! For boundary layer levels set depth to be constant.
  DO k=1, bl_levels
    DO j=1, rows
      DO i=1, row_length
        r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                   r_ref_theta(k)
        r_rho_levels(i,j,k)   = r_theta_levels(i,j,0) +         &
                                   r_ref_rho(k)
      END DO
    END DO
  END DO
  ! For intemediate levels use linear relaxation to constant value.
  ! set orographic heights.
  DO k=bl_levels+1, first_constant_r_rho_level-1
    DO j=1, rows
      DO i=1, row_length

        r_rho_levels(i,j,k) =                                   &
          ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
            r_theta_levels(i,j,bl_levels) ) *                   &
          ( eta_rho_levels(k) -                                 &
            eta_theta_levels(bl_levels) )/                      &
          (eta_rho_levels(first_constant_r_rho_level) -         &
           eta_theta_levels(bl_levels) )                        &
          +  r_theta_levels(i,j,bl_levels)

        r_theta_levels(i,j,k) =                                 &
          ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
            r_theta_levels(i,j,bl_levels) ) *                   &
          ( eta_theta_levels(k) -                               &
            eta_theta_levels(bl_levels) ) /                     &
          ( eta_rho_levels(first_constant_r_rho_level) -        &
            eta_theta_levels(bl_levels) )                       &
          +  r_theta_levels(i,j,bl_levels)

      END DO
    END DO
  END DO

CASE ( height_gen_smooth )
  ! A smooth quadratic height generation
  DO k=1, first_constant_r_rho_level-1
    DO j=1, rows
      DO i=1, row_length
        r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +  &
         planet_radius + orog(i,j) * (1.0 - eta_rho_levels(k)       &
              /eta_rho_levels(first_constant_r_rho_level))**2

        r_theta_levels(i,j,k) = eta_theta_levels(k) *               &
             z_top_of_model + planet_radius + Orog(i,j) *           &
             (1.0 - eta_theta_levels(k) /                           &
              eta_rho_levels(first_constant_r_rho_level))**2
      END DO
    END DO
  END DO

CASE DEFAULT
  ErrorStatus = 10
  WRITE (Cmessage,'(A)') 'Unrecognised height generation method - ',&
                     'Dump needs to be reconfigured'

  CALL ereport( RoutineName, ErrorStatus, Cmessage )
END SELECT

!  ----------------------------------------------------------
!  Set number of levels for non-local PBL scheme
!  ----------------------------------------------------------
IF (i_bl_vn == i_bl_vn_9b .OR. i_bl_vn == i_bl_vn_9c) THEN
  IF (z_nl_bl_levels > 0.0 .AND. z_nl_bl_levels < z_top_of_model) THEN 
  ! nl_bl_levels is the first theta-level at or above z_nl_bl_levels
  ! and no greater than bl_levels
    k=1
    DO WHILE ( k < bl_levels .AND.                                      &
             r_ref_theta(k) < z_nl_bl_levels )
      k=k+1
    END DO
    nl_bl_levels = k 
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(Cmessage,'(A,ES12.4,A,I4)')                                 &
                     ' z_nl_bl_levels=',z_nl_bl_levels,                 &
                     ' which gives nl_bl_levels = ',nl_bl_levels 
      CALL umPrint(Cmessage,src='calc_levels')
    END IF
  ELSE
    ErrorStatus = 1
    WRITE(cmessage,'(A,ES12.4,A)') 'calc_levels: Parameter z_nl_bl_levels=',&
          z_nl_bl_levels,' is outside allowed range'
    CALL ereport(RoutineName,ErrorStatus,cmessage)
  END IF ! test on z_nl_bl_levels
END IF

!  ----------------------------------------------------------
!  Set the values of alpha_cd
!  ----------------------------------------------------------
ALLOCATE(alpha_cd(bl_levels))
alpha_cd(1) = alpha_cd_in(1)
alpha_cd(2:bl_levels) = alpha_cd_in(2)

!  ----------------------------------------------------------
!  Set number of levels for stochastic perturbations
!  ----------------------------------------------------------
IF (i_pert_theta > off) THEN
  IF (zmin_pert_theta >= 0.0 .AND.                                      &
      zmin_pert_theta < z_top_of_model) THEN 
    ! minlev_pert_theta is the highest theta-level at or below 
    ! zmin_pert_theta and no greater than model_levels
    k=2
    DO WHILE ( k < model_levels .AND.                                   &
               r_ref_theta(k) < zmin_pert_theta )
      k=k+1
    END DO
    minlev_pert_theta = k-1
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(Cmessage,'(A,ES12.4,A,I4)')                                 &
             ' zmin_pert_theta=',zmin_pert_theta,                       &
             ' which gives minlev_pert_theta = ',minlev_pert_theta
      CALL umPrint(Cmessage,src='calc_levels')
    END IF
  ELSE
    ErrorStatus = 1
    WRITE(cmessage,'(A,ES12.4,A)')                                      &
            'calc_levels: Parameter zmin_pert_theta=',                  &
            zmin_pert_theta,' is outside allowed range'
    CALL ereport(RoutineName,ErrorStatus,cmessage)
  END IF ! test on zmin_pert_theta

  IF (zmax_pert_theta > 0.0 .AND. zmax_pert_theta < z_top_of_model) THEN 
    ! maxlev_pert_theta is the highest theta-level at or below 
    ! zmax_pert_theta and no greater than model_levels
    k=2
    DO WHILE ( k < model_levels .AND.                                   &
               r_ref_theta(k) < zmax_pert_theta )
      k=k+1
    END DO
    maxlev_pert_theta = k-1
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(Cmessage,'(A,ES12.4,A,I4)')                                 &
                   ' zmax_pert_theta=',zmax_pert_theta,                 &
                   ' which gives maxlev_pert_theta = ',maxlev_pert_theta
      CALL umPrint(Cmessage,src='calc_levels')
    END IF
  ELSE
    ErrorStatus = 1
    WRITE(cmessage,'(A,ES12.4,A)')                                      &
                   'calc_levels: Parameter zmax_pert_theta=',           &
                   zmax_pert_theta,' is outside allowed range'
    CALL ereport(RoutineName,ErrorStatus,cmessage)
  END IF ! test on zmax_pert_theta

  IF (maxlev_pert_theta < minlev_pert_theta) THEN 
    ! will give no perturbations so stop
    ErrorStatus = 1
    WRITE(cmessage,'(A)')                                               &
           'calc_levels: maxlev_pert_theta is less than '//             &
           'minlev_pert_theta so no perturbations will be generated'
    CALL ereport(RoutineName,ErrorStatus,cmessage)
  END IF ! test on maxlev_pert_theta

END IF ! test on i_pert_theta

RETURN

END SUBROUTINE calc_levels

