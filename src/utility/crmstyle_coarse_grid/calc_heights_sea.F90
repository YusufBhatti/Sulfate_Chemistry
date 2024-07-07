! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates the heights in the SCM.
MODULE calc_heights_sea_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_HEIGHTS_SEA_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
!  Calculate model heights as in UM
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------


SUBROUTINE calc_heights_sea                                            &
  ( height_gen_method, model_levels, bl_levels, rows, row_length,      &
                      in_rows, in_cols )


USE  crmwork_arrays_mod, ONLY: orog_full, r_theta_levels_local,        &
  r_rho_levels_local,r_theta_sea, r_rho_sea

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE vertnamelist_mod, ONLY:                                              &
    first_constant_r_rho_level, z_top_of_model                           &
    ,eta_theta, eta_rho

USE planet_constants_mod, ONLY: planet_radius

USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p


USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParCore, ONLY: mype, nproc_max

! For printing out info
USE umPrintMgr

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


! Inputs
INTEGER,INTENT(IN) ::       &
    model_levels            & ! Number of model levels
  , bl_levels               & ! Number of boundary layer levels
  , rows                    & ! Number of rows
  , row_length              & ! Number of columns
  , height_gen_method       & ! Method for generating model heights
  , in_rows                 & ! Full grid number of rows
  , in_cols                   ! Full grid number of columns


! Local variables

CHARACTER(LEN=errormessagelength) ::        &
  Cmessage                    ! Error message if ICODE >0

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CALC_HEIGHTS_SEA'

REAL ::                                                             &
    r_ref_theta(model_levels)                                       &
  , r_ref_rho(model_levels)

INTEGER ::                                                        &
    i,j,k, k1, k2, kloop
INTEGER :: ErrorStatus

INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

INTEGER, PARAMETER ::  &
  offx  = 1            &  ! one point halos
 ,offy  = 1               !


! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
!----------------------------------------------------------------------
!     Set up heights
!----------------------------------------------------------------------
WRITE(umMessage,'(A,F10.2,6I6)') ' Cal Heights ',z_top_of_model,             &
   height_gen_method, bl_levels, model_levels, first_constant_r_rho_level,   &
   rows, row_length
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' eta theta '
CALL umPrint(umMessage,src=RoutineName)

! Problem as number of levels unknown and want all levels. umPrint will
! only print the first 1024 characters of any buffer

kloop=(model_levels+9)/10

DO k=1,kloop
  k1=(k-1)*10+1
  k2=k*10
  IF (k2 > model_levels) THEN
    k2=model_levels
  END IF
  WRITE(umMessage,'(10E25.18)') (eta_theta(j),j=k1,k2)
  CALL umPrint(umMessage,src=RoutineName)
END DO

WRITE(umMessage,'(A)') ' eta rho '
CALL umPrint(umMessage,src=RoutineName)
DO k=1,kloop
  k1=(k-1)*10+1
  k2=k*10
  IF (k2 > model_levels) THEN
    k2=model_levels
  END IF
  WRITE(umMessage,'(10E25.18)') (eta_rho(j),j=k1,k2)
  CALL umPrint(umMessage,src=RoutineName)
END DO

! Set reference profile  - eta on model_levels +1 (first surface = 0.0)
!                          and array NOT (0:levs) but (1:levs+1)

DO k=1, model_levels
  r_ref_theta(k) = eta_theta(k+1) * z_top_of_model
  r_ref_rho(k)   = eta_rho(k)   * z_top_of_model
END DO

! levels above sea
DO k=1, model_levels
  r_theta_sea(k) =  r_ref_theta(k)+ planet_radius
  r_rho_sea(k)   =  r_ref_rho(k)+ planet_radius
END DO

WRITE(umMessage,'(A)')  ' r_theta_sea '
CALL umPrint(umMessage,src=RoutineName)
DO k=1,kloop
  k1=(k-1)*10+1
  k2=k*10
  IF (k2 > model_levels) THEN
    k2=model_levels
  END IF
  WRITE(umMessage,'(10E25.18)') (r_theta_sea(j),j=k1,k2)
  CALL umPrint(umMessage,src=RoutineName)
END DO

IF (mype == 0) THEN
  ! Set bottom level, ie orography
  DO j=1,in_rows
    DO i=1,in_cols
      r_theta_levels(i,j,0) = orog_full(i,j) + planet_radius
    END DO
  END DO

  ! For constant levels set r to be a constant on the level
  DO k=first_constant_r_rho_level, model_levels
    DO j=1,in_rows
      DO i=1,in_cols
        r_theta_levels(i,j,k) =  r_ref_theta(k)+ planet_radius
        r_rho_levels(i,j,k)   =  r_ref_rho(k)+ planet_radius
      END DO
    END DO
  END DO

  SELECT CASE( height_gen_method)
  CASE ( height_gen_original)
    ! The original version of height generation used in the SI dynamics
    !
    ! For boundary layer levels set depth to be constant.
    DO k=1, bl_levels
      DO j=1,in_rows
        DO i=1,in_cols
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
      DO j=1,in_rows
        DO i=1,in_cols
          r_rho_levels(i,j,k) =                                   &
            ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
              r_theta_levels(i,j,bl_levels) ) *                   &
            ( eta_rho(k) - eta_theta(bl_levels+1) )/              &
            (eta_rho(first_constant_r_rho_level) -                &
             eta_theta(bl_levels+1) )                             &
            +  r_theta_levels(i,j,bl_levels)

          r_theta_levels(i,j,k) =                                 &
            ( r_rho_levels(i,j,first_constant_r_rho_level) -      &
              r_theta_levels(i,j,bl_levels) ) *                   &
            ( eta_theta(k+1) - eta_theta(bl_levels+1) ) /         &
            ( eta_rho(first_constant_r_rho_level) -               &
              eta_theta(bl_levels+1) )                            &
            +  r_theta_levels(i,j,bl_levels)

        END DO
      END DO
    END DO

  CASE ( height_gen_smooth )
    ! A smooth quadratic height generation
    DO k=1, first_constant_r_rho_level-1
      DO j=1,in_rows
        DO i=1,in_cols
          r_rho_levels(i,j,k) = eta_rho(k) * z_top_of_model +        &
            planet_radius + orog_full(i,j) * (1.0 - eta_rho(k)       &
                /eta_rho(first_constant_r_rho_level))**2

          r_theta_levels(i,j,k) = eta_theta(k+1) *                      &
               z_top_of_model + planet_radius + Orog_full(i,j) *        &
               (1.0 - eta_theta(k+1) /                                  &
                eta_rho(first_constant_r_rho_level))**2
        END DO
      END DO
    END DO

  CASE DEFAULT
    ErrorStatus = 10
    WRITE (Cmessage,'(A,A)') 'Unrecognised height generation method - ',&
                       'Dump needs to be reconfigured'
  END SELECT

END IF  ! test on mype

! scatter r_theta and r_rho to PEs

DO k=1,model_levels

  ! DEPENDS ON: scatter_field
  CALL scatter_field( r_rho_levels_local(1,1,k), r_rho_levels(1,1,k),  &
                      row_length,rows,                                 &
                      in_cols,in_rows,                                 &
                      fld_type_p,halo_type_no_halo,                    &
                      0,gc_all_proc_group)

END DO

! scatter r_theta and r_rho to PEs
DO k=0,model_levels

  ! DEPENDS ON: scatter_field
  CALL scatter_field( r_theta_levels_local(1,1,k), r_theta_levels(1,1,k), &
                      row_length,rows,                                    &
                      in_cols,in_rows,                                    &
                      fld_type_p,halo_type_no_halo,                       &
                      0,gc_all_proc_group)

END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN

END SUBROUTINE calc_heights_sea
END MODULE calc_heights_sea_mod
