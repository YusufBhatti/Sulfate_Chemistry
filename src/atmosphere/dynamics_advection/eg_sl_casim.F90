! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_casim_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_CASIM_MOD'

CONTAINS
SUBROUTINE eg_sl_casim(                                                     &
              row_length, rows, model_levels, halo_i, halo_j, datastart,    &
              g_i_pe, high_order_scheme, monotone_scheme, l_high, l_mono,   &
              error_code )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_parcore,               ONLY: mype, nproc
USE nlsizes_namelist_mod,     ONLY: global_row_length
USE um_parvars,               ONLY: nproc_x,nproc_y, at_extremity,          &
                                    gc_proc_row_group, gc_proc_col_group
USE level_heights_mod,        ONLY: eta_theta_levels
USE atm_fields_bounds_mod,    ONLY: tdims
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf


USE departure_pts_mod,  ONLY: depart_xi1_w, depart_xi2_w, depart_xi3_w
USE ereport_mod,        ONLY: ereport
USE Field_Types,        ONLY: fld_type_p, fld_type_w
USE dynamics_input_mod, ONLY: l_sl_bc_correction

USE casim_prognostics, ONLY: cloudnumber, rainnumber, rain3mom,             &
                             icenumber, snownumber, snow3mom,               &
                             graupnumber, graup3mom, activesolliquid,       &
                             activesolrain, activeinsolice,                 &
                             activeinsolliquid, activesolnumber,            &
                             activeinsolnumber, activesolice

USE casim_switches, ONLY: l_mp_cloudnumber, l_mp_rainnumber, l_mp_rain3mom, &
                          l_mp_icenumber,   l_mp_snownumber, l_mp_snow3mom, &
                          l_mp_graupnumber, l_mp_graup3mom,                 &
                          l_mp_activesolliquid, l_mp_activesolrain,         &
                          l_mp_activeinsolice, l_mp_activesolice,           &
                          l_mp_activeinsolliquid, l_mp_activesolnumber,     &
                          l_mp_activeinsolnumber, n_casim_progs

USE mpp_conf_mod,   ONLY: swap_field_is_scalar

IMPLICIT NONE
!
! Description:
!   Find departure point timelevel n for prognostic variables associated
!   with the Cloud Aerosol Interacting Microphysics (CASIM).
!
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_CASIM'

! Model dimensions

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! MPP options

INTEGER, INTENT(IN) :: halo_i ! Size of halo in i.
INTEGER, INTENT(IN) :: halo_j ! Size of halo in j.
INTEGER, INTENT(IN) :: datastart(3) ! First gridpoints held by this processor.
INTEGER, INTENT(IN) :: g_i_pe( 1-halo_i : global_row_length+halo_i )
                       ! processor on my processor-row
                       ! holding a given value in i direction

! Integer parameters for advection
INTEGER, INTENT(IN) :: high_order_scheme  
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present

INTEGER, INTENT(IN) :: monotone_scheme
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.

! Error code
INTEGER, INTENT(INOUT) :: error_code  ! Non-zero on exit if error detected.

LOGICAL, INTENT(IN) :: l_high ! True, if high order interpolation required.
LOGICAL, INTENT(IN) :: l_mono ! True, if interpolation required to be monotone.

! Local variables
INTEGER :: i,j,k, number_of_inputs
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

INTEGER :: mvcnt ! moist variable count in super array
                 ! Avoids use of the Fortran keyword 'COUNT'

! tmp & dummy arrays

REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,    &
                    0:model_levels, n_casim_progs )

REAL :: super_array_out(1:row_length,1:rows,                            &
                        0:model_levels, n_casim_progs)

!------------------------------------------------------------------------------
! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

mvcnt = 0

IF ( l_mp_cloudnumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,cloudnumber)
   DO k = 0, model_levels
     DO j = tdims%j_start, tdims%j_end
       DO i = tdims%i_start, tdims%i_end
         super_array(i,j,k,mvcnt) = cloudnumber(i,j,k)
       END DO
     END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_cloudnumber

IF ( l_mp_rainnumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,rainnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = rainnumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_rainnumber

IF ( l_mp_rain3mom ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,rain3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = rain3mom(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_rain3mom

IF ( l_mp_icenumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,icenumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = icenumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_icenumber

IF ( l_mp_snownumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,snownumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = snownumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_snownumber

IF ( l_mp_snow3mom ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,snow3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = snow3mom(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_snow3mom

IF ( l_mp_graupnumber ) THEN   
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,graupnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = graupnumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_graupnumber

IF ( l_mp_graup3mom ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,graup3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = graup3mom(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_graup3mom

! CASIM microphysics extra activated aerosol prognostics
IF ( l_mp_activesolliquid ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activesolliquid) 
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = activesolliquid(i,j,k) 
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activesolliquid

IF ( l_mp_activesolrain ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activesolrain)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = activesolrain(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activesolrain

IF ( l_mp_activeinsolice ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activeinsolice)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            super_array(i,j,k,mvcnt) = activeinsolice(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activeinsolice

IF ( l_mp_activesolice ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activesolice)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            super_array(i,j,k,mvcnt) = activesolice(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activesolice

IF ( l_mp_activeinsolliquid ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activeinsolliquid)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            super_array(i,j,k,mvcnt) = activeinsolliquid(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activeinsolliquid

IF ( l_mp_activesolnumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activesolnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            super_array(i,j,k,mvcnt) = activesolnumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activesolnumber

IF ( l_mp_activeinsolnumber ) THEN
   mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array,mvcnt,activeinsolnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            super_array(i,j,k,mvcnt) = activeinsolnumber(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF ! l_mp_activeinsolnumber

! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                       &
             super_array, row_length, rows,                             &
             n_casim_progs*(model_levels+1),                            &
             halo_i, halo_j, fld_type_p,swap_field_is_scalar)

number_of_inputs = n_casim_progs

! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

CALL eg_interpolation_eta_pmf(                                          &
                     eta_theta_levels,fld_type_w,                       &
                     number_of_inputs,                                  &
                     row_length, rows, model_levels+1,                  &
                     rows,                                              &
                     row_length, rows, model_levels+1,                  &
                     high_order_scheme, monotone_scheme,                &
                     l_high, l_mono, depart_xi3_w, depart_xi1_w,        &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,       &
                     halo_i, halo_j,                                    &
                     global_row_length, datastart, at_extremity,        &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                     0, 0, error_code,                                  &
                     super_array, super_array_out,                      &
                     k_int_linear_in=k_int_linear)

! Reset moist variable count back to zero ready for unpacking of the 
! super array
mvcnt = 0

! CASIM extra prognostics output from the super array
IF ( l_mp_cloudnumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,cloudnumber)
  DO k = 0, model_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cloudnumber(i,j,k) = super_array_out(i,j,k,mvcnt)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF !l_mp_cloudnumber

! CASIM extra prognostics output from the super array
IF ( l_mp_rainnumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,rainnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            rainnumber(i,j,k) = super_array_out(i,j,k,mvcnt) 
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_rainnumber

IF ( l_mp_rain3mom ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,rain3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            rain3mom(i,j,k) = super_array_out(i,j,k,mvcnt) 
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_rain3mom

IF ( l_mp_icenumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,icenumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            icenumber(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_icenumber

IF ( l_mp_snownumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,snownumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            snownumber(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_snownumber

IF ( l_mp_snow3mom ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,snow3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            snow3mom(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_snow3mom

IF ( l_mp_graupnumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,graupnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
           graupnumber(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_graupnumber

IF ( l_mp_graup3mom ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,graup3mom)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            graup3mom(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_graup3mom
!
! End of the CASIM extra cloud and ice moment prognostics
! Begin adding the activated aerosol to the advection super array
!
IF ( l_mp_activesolliquid ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activesolliquid)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            activesolliquid(i,j,k) = super_array_out(i,j,k,mvcnt) 
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_activesolliquid

IF ( l_mp_activesolrain ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activesolrain)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            activesolrain(i,j,k) =  super_array_out(i,j,k,mvcnt) 
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_activesolrain

IF ( l_mp_activeinsolice ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activeinsolice)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            activeinsolice(i,j,k) = super_array_out(i,j,k,mvcnt)  
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_activeinsolice

IF ( l_mp_activesolice ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activesolice)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            activesolice(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_ActiveSolIce           

IF ( l_mp_activeinsolliquid ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activeinsolliquid)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            activeinsolliquid(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_activeinsolice  

IF ( l_mp_activesolnumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activesolnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            activesolnumber(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_ActiveSolNumber  

IF ( l_mp_activeinsolnumber ) THEN
  mvcnt = mvcnt + 1
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,tdims,super_array_out,mvcnt,activeinsolnumber)
   DO k = 0, model_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end 
            activeinsolnumber(i,j,k) = super_array_out(i,j,k,mvcnt)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF !l_mp_ActiveInSolNumber  

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_casim
END MODULE eg_sl_casim_mod
