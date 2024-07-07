! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  check the column integral of moist increments
!
MODULE check_dmoist_inc_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CHECK_DMOIST_INC_MOD'

CONTAINS

SUBROUTINE check_dmoist_inc(row_length, rows, model_levels,                 &
                            delta_lambda,delta_phi,timestep,                &
                            p_layer_boundaries,rho_r2,                      &
                            dqtotalbydt,                                    &
                            ep4,cprint,                                     &
                            ep1, ep2, ep3)

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:       &
   tdims, pdims, pdims_s


! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

! Trig information
USE trignometric_mod,  ONLY:                                      &
  FV_cos_theta_latitude

USE timestep_mod, ONLY: timestep_number
USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParCore, ONLY: mype
USE planet_constants_mod, ONLY: g

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr                   ! Everything connected with umMessage

USE Field_Types
USE global_2d_sums_mod, ONLY: global_2d_sums

IMPLICIT NONE
!------------------------------------------------------------------------
! Description:
!  Use 3 different methods to evaluate the column integral of the moist
! increments.
! 1. On rho grid
! 2. On theta grid
! 3. On theta grid but using pressure information rather than density
!
! Either wet or dry density is passed in depending on whether the values are
! in specific humidities or mixing ratios. Note the pressure integral is not
! really valid in the case of mixing ratios as it is effectively using a wet
! density.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine arguments

! Model dimensions
INTEGER, INTENT(IN)  ::        &
    row_length                 & ! Row length
  , rows                       & ! Number of rows
  , model_levels                 ! Number of model levels


REAL, INTENT(IN) ::            &
  delta_lambda                 &
 ,delta_phi                    & ! model grid spacing in radians
 ,timestep                       ! model timestep

REAL, INTENT(IN) ::                            &
  p_layer_boundaries(pdims%i_end,              & ! pressure at layer boundaries.
                      pdims%j_end,             & ! Same as p except at
                               0:pdims%k_end)  & ! bottom level = pstar
 ,rho_r2(pdims_s%i_start:pdims_s%i_end,        & ! density *r*r (kg/m)
        pdims_s%j_start:pdims_s%j_end,         & ! Either wet density (specfic
        pdims_s%k_start:pdims_s%k_end)         & ! case) or dry density (mixing)
 ,dqtotalbydt(tdims%i_start:tdims%i_end,       & ! total moist increment
              tdims%j_start:tdims%j_end,       & ! (kg/kg/s)
                          1:tdims%k_end)


REAL, INTENT(IN) ::               &
  ep4(row_length, rows)             ! diagnostic Evap/Precip (kg/m2/s)

CHARACTER(LEN=60), INTENT(IN) ::  &
  cprint
REAL, INTENT(OUT) ::              &
  ep1(row_length, rows)           & ! Evap/Precip Method 1
 ,ep2(row_length, rows)           & ! Evap/Precip Method 2
 ,ep3(row_length, rows)             ! Evap/Precip Method 3


!----------------------------------------------------------------------
! local variables
!----------------------------------------------------------------------
INTEGER ::            &
 i, j, k                ! Loop counters

INTEGER, PARAMETER :: &
  n_sums =4           &  ! number of global sums
 ,ip_1 = 1            &  ! Method 1
 ,ip_2 = 2            &  ! Method 2
 ,ip_3 = 3            &  ! Method 3
 ,ip_4 = 4               ! Method 4

REAL ::                       &
  weight1, weight2, weight3   & ! weights
 , tempd, tempw, ww2, ww1     &
 , ra2, a2                    &
 , factor                        ! grid resolution factor

REAL ::                                     &
  rho_only(row_length,rows,model_levels)    & ! rho without r*r term
 ,rho_theta(row_length,rows,model_levels)   & ! rho on theta with r*r term
 ,delr_rho(row_length,rows,model_levels)    & ! dr for rho layers
 ,delr_theta(row_length,rows,model_levels)  & ! dr for theta layers
 ,dp_theta(row_length,rows,model_levels)      ! dp for theta layers

REAL ::                                     &
  dqt_rho(row_length,rows,model_levels)       ! dmoist on rho levels

REAL ::                                     &
  vert_int_array(row_length,rows,n_sums)    & ! vertical integrals
 ,sum_results(n_sums)                         ! sum of SUM_ARRAY


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_DMOIST_INC'
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! zero summation arrays
!----------------------------------------------------------------------


DO k=1,n_sums
  sum_results(k)=0.0
END DO

DO k=1,n_sums
  DO j = 1, rows
    DO i=1,row_length
      vert_int_array(i,j,k) = 0.0
    END DO
  END DO
END DO
!----------------------------------------------------------------------
! Calculate layer thickness for rho layers, rho without r2, and dp
!----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,ww1,ww2,weight1,weight2,weight3) &
!$OMP SHARED(model_levels,rows,row_length,delr_rho,rho_only,dp_theta,  &
!$OMP r_theta_levels,rho_r2,r_rho_levels,p_layer_boundaries,delr_theta,  &
!$OMP rho_theta,dqt_rho,dqtotalbydt)

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels
  DO j = 1, rows
    DO i=1,row_length
      delr_rho(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
      rho_only(i,j,k) = rho_r2(i,j,k)/(r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
      dp_theta(i,j,k) = p_layer_boundaries(i,j,k)                         &
                                     - p_layer_boundaries(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END DO

!----------------------------------------------------------------------
! Calculate layer thickness for theta layers
!----------------------------------------------------------------------
k=1
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i=1,row_length
    delr_theta(i,j,k) = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k-1)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k=2,model_levels-1
  DO j = 1, rows
    DO i=1,row_length
      delr_theta(i,j,k) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

k=model_levels
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i=1,row_length
    delr_theta(i,j,k) = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
  END DO
END DO
!$OMP END DO

!----------------------------------------------------------------------
! rho on theta levels * r2
!----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k=1,model_levels-1
  DO j = 1, rows
    DO i=1,row_length
      weight1 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
      weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight3 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      ww1 = weight1/weight3
      ww2 = weight2/weight3
      rho_theta(i,j,k) = (ww2*rho_only(i,j,k+1) + ww1*rho_only(i,j,k)) *   &
                         r_theta_levels(i,j,k)*r_theta_levels(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

k=model_levels
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i=1,row_length
    rho_theta(i,j,k) = rho_only(i,j,k)*                                  &
                       r_theta_levels(i,j,k)*r_theta_levels(i,j,k)
  END DO
END DO
!$OMP END DO


!-------------------------------------------------------------------
! Interpolate increments to rho points
! Note needs to remain consistent with interpolation used in dynamics.
! Using linear interpolation
!
!                      K               for rho
!      K                          K-1  for theta
!      X<--- w1------->X<-- w2--->X
!       <----------w3------------>
!-------------------------------------------------------------------
k=1
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i=1,row_length
    ! assume bottom rho level value equal to bottom theta level value
    dqt_rho(i,j,k)     = dqtotalbydt(i,j,k)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k=2,model_levels
  DO j = 1, rows
    DO i=1,row_length
      weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight2 = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
      weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
      ww1 = weight1/weight3
      ww2 = weight2/weight3
      dqt_rho(i,j,k)     = ww2*dqtotalbydt(i,j,k)  + ww1*dqtotalbydt(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL 
!-------------------------------------------------------------------
! vertical integrals
!-------------------------------------------------------------------

DO k=1,model_levels
  DO j = 1, rows
    DO i=1,row_length

      ! Method 1   - rho levels
      vert_int_array(i,j,ip_1) = vert_int_array(i,j,ip_1) +  &
                             rho_r2(i,j,k)*delr_rho(i,j,k)*dqt_rho(i,j,k)

      ! Method 2   - theta levels
      vert_int_array(i,j,ip_2) = vert_int_array(i,j,ip_2) +  &
                         rho_theta(i,j,k)*delr_theta(i,j,k)*dqtotalbydt(i,j,k)


      ! Method 3   - pressure integral
      vert_int_array(i,j,ip_3) = vert_int_array(i,j,ip_3)- dp_theta(i,j,k)/g * &
                                                          dqtotalbydt(i,j,k)
    END DO
  END DO
END DO

! copy results back to output arrays
DO j = 1, rows
  DO i=1,row_length
    a2  = (r_theta_levels(i,j,0)*r_theta_levels(i,j,0))
    ra2 = 1.0/a2
    ep1(i,j) = vert_int_array(i,j,ip_1)*ra2
    ep2(i,j) = vert_int_array(i,j,ip_2)*ra2
    ep3(i,j) = vert_int_array(i,j,ip_3)
    ! convert to amount (not per unit area)
    vert_int_array(i,j,ip_3) = vert_int_array(i,j,ip_3)*a2
    ! Method 4   - diagnostic copy kg/m2/s so want to convert to kg/s
    vert_int_array(i,j,ip_4) = ep4(i,j)*a2
  END DO
END DO

!----------------------------------------------------------------------
DO k=1,n_sums
  DO j = 1, rows
    DO i=1,row_length
      vert_int_array(i,j,k) = vert_int_array(i,j,k)*                      &
                                     FV_cos_theta_latitude(i,j)
    END DO
  END DO
END DO

CALL global_2d_sums(vert_int_array, row_length, rows, 0, 0,       &
                          n_sums, sum_results, gc_all_proc_group)

factor=delta_lambda*delta_phi   ! already in radians


IF (mype == 0) THEN
  WRITE(umMessage,'(2A)')'Global precip/evap from ',cprint
  CALL umPrint(umMessage,src='check_dmoist_inc')

  WRITE(umMessage,'(A,E25.18)') 'Method 1       ',sum_results(ip_1)*factor
  CALL umPrint(umMessage,src='check_dmoist_inc')
  WRITE(umMessage,'(A,E25.18)') 'Method 2       ',sum_results(ip_2)*factor
  CALL umPrint(umMessage,src='check_dmoist_inc')
  WRITE(umMessage,'(A,E25.18)') 'Method 3       ',sum_results(ip_3)*factor
  CALL umPrint(umMessage,src='check_dmoist_inc')
  WRITE(umMessage,'(A,E25.18)') 'Method 4       ',sum_results(ip_4)*factor
  CALL umPrint(umMessage,src='check_dmoist_inc')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_dmoist_inc

END MODULE check_dmoist_inc_mod
