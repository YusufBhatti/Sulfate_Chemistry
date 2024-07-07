! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Force idealised vertical sine shaped heating at a location

MODULE idl_local_heat_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='IDL_LOCAL_HEAT_MOD'

CONTAINS

SUBROUTINE idl_local_heat(temp_inc)

! Purpose: To apply a forcing increment to a field at a location for a
!          period of time
!
! Method:  Works out heating for location based on a function of height 
!          and distance from the central location
!

USE idealise_run_mod, ONLY: l_ideal_2d 
USE local_heat_mod, ONLY: local_heat_xoffset,                            &
 local_heat_yoffset, local_heat_amp, local_heat_sigma, local_heat_base,  &
 local_heat_top, local_heat_period

USE conversions_mod,     ONLY: pi, rsec_per_day
USE planet_constants_mod, ONLY: planet_radius

USE timestep_mod, ONLY:                                             &
   timestep, timestep_number

USE level_heights_mod, ONLY:                                        &
   r_theta_levels

USE atm_fields_bounds_mod, ONLY: tdims

! grid info
USE nlsizes_namelist_mod,   ONLY: global_row_length, global_rows
USE horiz_grid_mod,         ONLY: glob_xi1_u, glob_xi2_v, xi1_p, xi2_p

USE umPrintMgr, ONLY: umMessage,umPrint
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE


REAL, INTENT(INOUT)  ::                &
  temp_inc(tdims%i_start:tdims%i_end,  & ! Theta increment in K/s
           tdims%j_start:tdims%j_end,  &
           tdims%k_start:tdims%k_end)

!---------------------------------------------------------------------------
! Local variables

INTEGER ::        &
  i,j,k, ii, jj      ! loop counters

LOGICAL ::        &
  l_inc             ! increments to be applied this timestep

REAL  ::                 &
  factor                 & ! scaling factor to apply to inc
 ,per_sec_factor         & ! Factor converting data in namelist
                           ! to increment per second
 ,temp_height            & ! height of model level above surface (m) 
 ,time_now               & ! current time in (s)
 ,dz                     & ! distance from base of heating
 ,dx                     & ! distance from centre in x direction
 ,dy                     & ! distance from centre in y direction
 ,x0                     & ! Heating centre in x direction (m)
 ,y0                     & ! Heating centre in y direction (m)
 ,dr_over_sig_sq           ! (r/sigma)**2

REAL  ::                                    &
  t_heat_profile(tdims%k_start:tdims%k_end) & ! Sine heating profile on 
                                              ! model levels 
 ,weight(tdims%i_start:tdims%i_end,         & ! weight for heating in horizontal
         tdims%j_start:tdims%j_end)


CHARACTER(LEN=*), PARAMETER :: RoutineName='IDL_LOCAL_HEAT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_inc=.FALSE.

! Applying heat over a period at the begining of the run unitl reached
! local_heat_period

time_now = REAL(timestep_number) * timestep

IF (time_now <= local_heat_period) THEN
  l_inc = .TRUE.
END IF

! Convert from K/day to K/sec
per_sec_factor = 1.0/rsec_per_day

!---------------------------------------------------------------------------
! Heating profile as used in Oliver Halliday's analytic calculations. 
!---------------------------------------------------------------------------

IF (l_inc) THEN

  ! Coordinates of heating centre in m
  x0 = glob_xi1_u(0) + local_heat_xoffset          &
      *(glob_xi1_u(global_row_length)-glob_xi1_u(0))

  y0 = glob_xi2_v(0) + local_heat_yoffset          &
      *(glob_xi2_v(global_rows) - glob_xi2_v(0))

  ! Heating profile as a function of height - assumes a flat surface
  ! In terms of K/s

  factor = pi/(local_heat_top - local_heat_base)
  WRITE(umMessage,'(a,e16.8,2(a,e16.8))')                                   &
  'Applying a sine heating - factor = ',factor, ' centre x = ',x0,' y= ',y0
  CALL umPrint(umMessage,src='idl_force_sine_at')
  
  DO k = tdims%k_start, tdims%k_end
    temp_height = r_theta_levels(1,1,k) - planet_radius  
    IF (temp_height > local_heat_base .AND. temp_height < local_heat_top) THEN
      dz = temp_height - local_heat_base    
      t_heat_profile(k) = local_heat_amp * SIN(factor*dz)*per_sec_factor
    ELSE
      t_heat_profile(k) = 0.0
    END IF
  END DO

  ! Weight   
  ! r direction    F(r) =  exp (-0.5(r/sigma)**2)

  IF (l_ideal_2d) THEN      ! 2d heating same for all y, centred on x0

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dx = (xi1_p(i)-x0)
        dr_over_sig_sq = (dx/local_heat_sigma)**2
        weight(i,j) = EXP(-0.5 * dr_over_sig_sq)
      END DO 
    END DO 

  ELSE                ! 3d heating centred on x0, y0

    DO j = tdims%j_start, tdims%j_end
      dy = xi2_p(j)-y0
      DO i = tdims%i_start, tdims%i_end
        dx = xi1_p(i)-x0
        dr_over_sig_sq = (dx*dx +dy*dy)/(local_heat_sigma*local_heat_sigma)
        weight(i,j) = EXP(-0.5 * dr_over_sig_sq)
      END DO 
    END DO 

  END IF

  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        temp_inc(i,j,k) = t_heat_profile(k) * weight(i,j) * timestep
      END DO
    END DO
  END DO 

ELSE

  ! No heating increment being applied    
  temp_inc(:,:,:) = 0.0

END IF ! test l_inc


! End of routine.
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE idl_local_heat

END MODULE idl_local_heat_mod
