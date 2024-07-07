! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_zlf_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_ZLF_MOD'
!
! Description:
!   This mod and its parameters are for  ZLF scheme (LAM only)
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP programming standards
!
!
! The following parameters and logicals are for controlling the ZLF scheme.
! ZLF is scheme in which mass conservation is enforced at the advection stage.
! It can be applied in both global and LAM configurations. However, this
! was developed primary for use with LAMs. Mass conservation for the global
! model is enforced at the end of the time step. However, there is no 
! theoretical or technical objections why this cannot also be applied 
! to the global model too at the advection stage. 
!
! For LAMs ZLF will zero data in part (zlf_np_zerod) of the Rim zone
! then perform a SL advection and enforce conservation then overwrite
! part (i.e., zlf_np_overwrite) rim zone 
!
! The parameters (integers and logicals) are: 
!
! zlf_np_weights_ones = number of points in RIM with blending weights=1
! zlf_np_zerod        = number of points in RIM to be zeroed for ZLF
! zlf_np_overwrite    = number of points in RIM to be overwritten after 
!                       ZLF
! zlf_np_linear       = number of points in RIM to be SL linear
! zlf_non_zeros_safe_width = number of points non-zeroed in RIM under
!                      which ZLF is safe to operate CFL<= this number
!                            
! zlf_cfl_top_level = Vertical level above which CFL violations are 
!             ignored. This depends on the user specified a height above
!             "zlf_maximum_height"  which defines the height where 
!             there are no significant moisture (practically dry).
!                                       
! zlf_zeros_rim_default  = This logical specifies whether we zeroed fields
!                 in the RIM or not (Global or cyclic domains = F, Lam=T)
!                                                
! zlf_zeros_rim_step  = specifies whether during this step we zeroed the
!                 field or not. This logical permits to switch off ZLF
!            if the CFL doesn't allow the safe use of ZLF for that step
!                          
! zlf_linear_boundaries=T; use linear SL at the boundaries of the LAM
!                          this is the default if ZLF is invoked.
! zlf_zeros_rim_default=T; if ZLF is used and model_type = LAM
! zlf_zeros_rim_default=F; if ZLF is used and model_type = global
!
! zlf_overwrite_rim_default=T (for LAM) F (for global).  
!                          
! This to overwrite the "zlf_np_overwrite" points in the rim with data 
! at start of the time step. This is just to prevent physics from dealing 
! with the artificial zeroing because this region will be overwritten by 
! the LBCs anyway at the end of the time step.
!
! zlf_para_set = This logical is here to allow ZLF setting to happen
!                only once if ZLF is used with several variables  
!

INTEGER :: zlf_np_weights_one
INTEGER :: zlf_np_zerod
INTEGER :: zlf_np_overwrite
INTEGER :: zlf_np_linear
INTEGER :: zlf_non_zeros_safe_width
INTEGER :: zlf_cfl_top_level_moist
INTEGER :: zlf_cfl_top_level_tracers
INTEGER :: zlf_cfl_top_level_theta
LOGICAL :: zlf_zeros_rim_default
LOGICAL :: zlf_zeros_rim_step
LOGICAL :: zlf_overwrite_rim_default
LOGICAL :: zlf_linear_boundaries
LOGICAL :: zlf_para_set = .FALSE.
!
! These routines are related to the ZLF scheme
!
CONTAINS
!
SUBROUTINE eg_zlf_setup(kks,kkf)
!
! kks == start of w (or theta) levels
! kkf == end   of w (or theta) levels
!
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE level_heights_mod,  ONLY: r_theta_levels, eta_theta_levels
USE dynamics_input_mod, ONLY: zlf_maximum_height
USE planet_constants_mod,  ONLY: planet_radius
USE lbc_read_data_mod,  ONLY: rimweightsa
USE model_domain_mod,   ONLY: model_type, mt_lam
USE dynamics_input_mod, ONLY: zlf_conservation_moist_option,   &
                              zlf_conservation_tracers_option, &
                              zlf_conservation_theta_option

IMPLICIT NONE
!
! Description:
!   setup the parameters (integers and logicals) of ZLF 
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='EG_ZLF_SETUP'

INTEGER,INTENT(IN) :: kks,kkf 
INTEGER            :: k
REAL, PARAMETER    :: epsmz=EPSILON(1.0)
REAL               :: h_atmos

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! compute the level above which ZLF ignores Courant number violations
!
! Note that z is approximated as z=eta*H, where H=height of domain, 
! excluding orography. This is so the level "zlf_cfl_top_level" above
! which ZLF will switch off will be unique irrespective of the domain
! decomposition.
! 
! ==> therefore the height corresponding to "zlf_cfl_top_level" level
! will be approximately "zlf_maximum_height" and this is OK for this 
! purpose
!

zlf_cfl_top_level_moist = kks - 1
h_atmos = r_theta_levels(1,1,kkf)-planet_radius
DO k = kks, kkf
  IF ( h_atmos*eta_theta_levels(k) < zlf_maximum_height ) THEN
    zlf_cfl_top_level_moist = k
  END IF
END DO
zlf_cfl_top_level_tracers = kkf
zlf_cfl_top_level_theta   = kkf

zlf_np_weights_one=INT(SUM(rimweightsa,mask=rimweightsa  >= (1.0-epsmz)))
zlf_np_zerod = INT( 0.5*(zlf_np_weights_one - 1) )
zlf_np_overwrite = zlf_np_weights_one
zlf_np_linear = zlf_np_overwrite + 1
zlf_linear_boundaries=.TRUE.

IF ( model_type == mt_lam ) THEN
  zlf_zeros_rim_default     =.TRUE.
  zlf_zeros_rim_step        =.TRUE.
  zlf_overwrite_rim_default =.TRUE.
ELSE
  zlf_zeros_rim_default     =.FALSE.
  zlf_zeros_rim_step        =.FALSE.
  zlf_overwrite_rim_default =.FALSE.     
END IF 
     
zlf_non_zeros_safe_width = zlf_np_weights_one - zlf_np_zerod - 1 + &
                            MIN(zlf_np_linear-zlf_np_weights_one,1 )
zlf_para_set = .TRUE.
                                   
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE eg_zlf_setup
!
!==========================================================================
!

SUBROUTINE eg_zlf_zero_rim(tracers,                                  &
                           iis, iif, jjs, jjf, kks, kkf, n_tracers,  &
                           zlf_cfl_top_level_local,                  &
                           zlf_conserv_option_local                  )
                           

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod
USE horiz_grid_mod
USE departure_pts_mod
USE um_parvars, ONLY: at_extremity
USE um_parparams, ONLY: PSouth,PNorth,PEast,PWest
USE um_parcore, ONLY: nproc
USE umPrintMgr, ONLY: umPrint, umMessage

IMPLICIT NONE

!
! Description:
!   Zero external halos and part of the rim (LAM only)
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='EG_ZLF_ZERO_RIM'

INTEGER, INTENT(IN)     :: iis, iif, jjs, jjf, kks, kkf, n_tracers
INTEGER, INTENT(IN)     :: zlf_cfl_top_level_local
INTEGER, INTENT(INOUT)  :: zlf_conserv_option_local
REAL, INTENT(INOUT)     :: tracers(iis:iif,jjs:jjf,kks:kkf,n_tracers)

REAL    :: rim_cfl, rim_max_cfl ! Non-dimensionalised CFL within the rim
INTEGER :: i,j,k,i_tr           ! Loopers
INTEGER :: ierr                 ! GCOM error code
INTEGER :: np_cfl_check         ! Number of points the CFL check is performed on
REAL    :: h_disp_max, h_disp   ! Horizontal displacements
INTEGER :: ii, jj, im, jm       ! Indices in CFL calculation

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

zlf_zeros_rim_step = zlf_zeros_rim_default

!
! Compute the RIM CFL number to see if its safe to apply ZLF or not
! This is done by finding the horizontal displacement within the rim; 
! this is essentially the number of gridpoints (in either X or Y co-ordinate)
! between the departure point and the arrival point.

rim_max_cfl  = 0.0
np_cfl_check = zlf_np_overwrite + 1
h_disp_max   = 0.0 
 
IF ( at_extremity(PWest) ) THEN 
  i = np_cfl_check
  DO j = tdims%j_start, tdims%j_end
    DO k = tdims%k_start, zlf_cfl_top_level_local                  
      h_disp      = xi1_p(i) - depart_xi1_w(i,j,k)        
      h_disp_max  = MAX(h_disp_max, h_disp)
    END DO
  END DO
  ii = i - 1
  DO WHILE ( h_disp_max > (xi1_p(i)-xi1_p(ii)) )
    ii = ii - 1 
  END DO
  im = ii
  h_disp = xi1_p(i) - xi1_p(im+1)
  rim_cfl = ( h_disp_max - h_disp      ) / &
            ( xi1_p(im+1) - xi1_p(im)  )     
  rim_cfl = REAL(i-im-1) + rim_cfl
  rim_max_cfl = MAX(rim_max_cfl,rim_cfl)           
END IF 
        
IF ( at_extremity(PEast)  ) THEN
  i = tdims%i_end-np_cfl_check+1
  DO j = tdims%j_start, tdims%j_end
    DO k = tdims%k_start, zlf_cfl_top_level_local      
      h_disp     = depart_xi1_w(i,j,k) - xi1_p(i)
      h_disp_max = MAX(h_disp_max, h_disp)             
    END DO
  END DO
  ii = i + 1
  DO WHILE ( h_disp_max > (xi1_p(ii)-xi1_p(i)) )
    ii = ii + 1 
  END DO
  im = ii - 1
  h_disp  = xi1_p(im) - xi1_p(i)
  rim_cfl = ( h_disp_max - h_disp      ) / &
            ( xi1_p(im+1) - xi1_p(im)  )     
  rim_cfl = REAL(im-i) + rim_cfl
  rim_max_cfl = MAX(rim_max_cfl,rim_cfl)    
END IF

IF ( at_extremity(PSouth) ) THEN
  j = np_cfl_check 
  DO i = tdims%i_start, tdims%i_end
    DO k = tdims%k_start, zlf_cfl_top_level_local     
      h_disp      = xi2_p(j)-depart_xi2_w(i,j,k)        
      h_disp_max  = MAX(h_disp_max, h_disp)              
    END DO
  END DO
  jj = j - 1
  DO WHILE ( h_disp_max > (xi2_p(j)-xi2_p(jj)) )
    jj = jj - 1 
  END DO
  jm = jj
  h_disp = xi2_p(j) - xi2_p(jm+1)
  rim_cfl = ( h_disp_max - h_disp      ) / &
            ( xi2_p(jm+1) - xi2_p(jm)  )     
  rim_cfl = REAL(j-jm-1) + rim_cfl
  rim_max_cfl = MAX(rim_max_cfl,rim_cfl)         
END IF
    
IF ( at_extremity(PNorth) ) THEN
  j = tdims%j_end-np_cfl_check+1
  DO i = tdims%i_start, tdims%i_end
    DO k = tdims%k_start, zlf_cfl_top_level_local       
      h_disp      = depart_xi2_w(i,j,k) - xi2_p(j)        
      h_disp_max  = MAX(h_disp_max, h_disp)     
    END DO
  END DO
  jj = j + 1
  DO WHILE ( h_disp_max > (xi2_p(jj)-xi2_p(j)) )       
    jj = jj + 1 
  END DO
  jm = jj - 1
  h_disp  = xi2_p(jm) - xi2_p(j)
  rim_cfl = ( h_disp_max - h_disp      ) / &
            ( xi2_p(jm+1) - xi2_p(jm)  )     
  rim_cfl = REAL(jm-j) + rim_cfl
  rim_max_cfl = MAX(rim_max_cfl,rim_cfl)
END IF

! Ensure all processes are aware of the maximum rim CFL.
! All processes need this to ensure that ZLF is switched off in a 
! consistent manner.
CALL gc_rmax(1,nproc,ierr,rim_max_cfl)

IF ( rim_max_cfl > REAL(zlf_non_zeros_safe_width) ) THEN
  WRITE(ummessage,'(A,F16.4,I0)') &
  'WARNING: RIM CFL > zeroed zone = ',rim_max_cfl,zlf_non_zeros_safe_width
  CALL umPrint(ummessage,src='eg_zlf_zero_rim')
  WRITE(ummessage,'(A)') &
  'WARNING: RIM CFL > zeroed zone => ZLF may not work properly'
  CALL umPrint(ummessage,src='eg_zlf_zero_rim')
  WRITE(ummessage,'(A)') &
  'WARNING: RIM CFL > zeroed zone ==> ZLF is swiched off'
  CALL umPrint(ummessage,src='eg_zlf_zero_rim')
  !
  ! If (cfl > zlf_non_zeros_safe_width) then ZLF may not work properly
  ! In this case set "zlf_conserv_option_local=0" & "zlf_zeros_rim_step=F"
  ! which is equivalent to switching ZLF off
  !  
  zlf_conserv_option_local = 0
  zlf_zeros_rim_step       = .FALSE.

END IF 
 
! If we require it, we now zero the tracers in each of the rims
IF ( zlf_zeros_rim_step ) THEN
!$OMP PARALLEL DEFAULT (NONE) PRIVATE(i,j,k,i_tr)                       &
!$OMP SHARED(at_extremity, n_tracers, kks, kkf, jjs, jjf, iif, iis,     &
!$OMP        zlf_np_zerod, tracers)
  IF ( at_extremity(PWest) ) THEN 
    DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
      DO k = kks, kkf
        DO j = jjs, jjf
          DO i = iis, zlf_np_zerod
            tracers(i,j,k,i_tr) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF 

  IF ( at_extremity(PEast)  ) THEN   
    DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
      DO k = kks, kkf
        DO j = jjs, jjf
          DO i = iif-zlf_np_zerod+iis, iif
            tracers(i,j,k,i_tr) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

  IF ( at_extremity(PSouth) ) THEN
    DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
      DO k = kks, kkf
        DO j = jjs, zlf_np_zerod
          DO i = iis, iif
            tracers(i,j,k,i_tr) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF

  IF ( at_extremity(PNorth) ) THEN   
    DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
      DO k = kks, kkf
        DO j = jjf-zlf_np_zerod+jjs, jjf
          DO i = iis, iif
            tracers(i,j,k,i_tr) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
  END IF
!$OMP END PARALLEL
END IF 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE eg_zlf_zero_rim
!
!==========================================================================
!
SUBROUTINE eg_zlf_overwrite_rim(tracers,  tracers_overwrite,     &
                                its, itf, jts, jtf, kts, ktf,    &
                                ios, iof, jos, jof,              &
                                n_tracers, np_overwrite          )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod
USE um_parvars, ONLY: at_extremity
USE um_parparams, ONLY: PSouth,PNorth,PEast,PWest
USE umPrintMgr, ONLY: umPrint, umMessage

IMPLICIT NONE
!
! Description:
!  Overwrite the rim region  of tracers with tracers_overwrite
!
!

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='EG_ZLF_OVERWRITE_RIM'

INTEGER, INTENT(IN) :: its, itf, jts, jtf, kts, ktf
INTEGER, INTENT(IN) :: ios, iof, jos, jof, np_overwrite, n_tracers
REAL, INTENT(INOUT) ::           tracers(its:itf,jts:jtf,kts:ktf,n_tracers )
REAL, INTENT(IN)    :: tracers_overwrite(ios:iof,jos:jof,kts:ktf,n_tracers )

INTEGER :: i,j,k,i_tr     ! Loopers

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
 
IF ( at_extremity(PWest) ) THEN          
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr,k,j,i)                &
!$OMP SHARED(n_tracers, kts, ktf, tdims, np_overwrite, tracers, &
!$OMP        tracers_overwrite)
  DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
    DO k = kts, ktf
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, np_overwrite       
          tracers(i,j,k,i_tr) = tracers_overwrite(i,j,k,i_tr)
        END DO
      END DO  
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
END IF 
IF ( at_extremity(PEast)  ) THEN        
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr,k,j,i)                &
!$OMP SHARED(n_tracers, kts, ktf, tdims, np_overwrite, tracers, &
!$OMP        tracers_overwrite)
  DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
    DO k = kts, ktf
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_end-np_overwrite+1,tdims%i_end
          tracers(i,j,k,i_tr) = tracers_overwrite(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
END IF
IF ( at_extremity(PSouth) ) THEN        
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr,k,j,i)                &
!$OMP SHARED(n_tracers, kts, ktf, tdims, np_overwrite, tracers, &
!$OMP        tracers_overwrite)
  DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
    DO k = kts, ktf
      DO i = tdims%i_start, tdims%i_end
        DO j = tdims%j_start, np_overwrite    
          tracers(i,j,k,i_tr) = tracers_overwrite(i,j,k,i_tr)
        END DO
      END DO          
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
END IF
IF ( at_extremity(PNorth) ) THEN        
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr,k,j,i)                &
!$OMP SHARED(n_tracers, kts, ktf, tdims, np_overwrite, tracers, &
!$OMP        tracers_overwrite)
  DO i_tr = 1, n_tracers
!$OMP DO SCHEDULE(STATIC)
    DO k = kts, ktf
      DO i = tdims%i_start, tdims%i_end
        DO j = tdims%j_end-np_overwrite+1, tdims%j_end   
          tracers(i,j,k,i_tr) = tracers_overwrite(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE eg_zlf_overwrite_rim
END MODULE eg_zlf_mod
