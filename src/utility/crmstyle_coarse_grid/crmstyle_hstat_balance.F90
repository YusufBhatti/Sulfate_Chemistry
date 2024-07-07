! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine crmstyle_hstat_balance
MODULE crmstyle_hstat_balance_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
           ModuleName = 'CRMSTYLE_HSTAT_BALANCE_MOD'

CONTAINS

SUBROUTINE crmstyle_hstat_balance(mlevs)

! Purpose:
!          Calculates hydrostatically balanced Exner from dry potential 
!          temperature (theta) and moisture variables.
!
!
USE word_sizes_mod, ONLY: iwp, wp   ! Allows use of 4 byte words to reduce
                                    ! memory

USE planet_constants_mod, ONLY: cp, kappa, pref, recip_epsilon, g
USE missing_data_mod, ONLY: rmdi

USE crmstyle_cntl_mod, ONLY:                                                &
   l_qgraup, l_no_orog

USE crmstyle_grid_info_mod, ONLY:                                           &
   local_new_x, local_new_y, local_row_len, local_rows

USE crmwork_arrays_mod, ONLY:                                               &
  h_theta_sea, h_rho_sea, mask, index_col, index_row

USE crmstyle_sample_arrays_mod, ONLY:                                       &
  all_pstar, all_tstar, all_th, p_theta_hydro,                              &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup

USE hires_data_mod, ONLY:                                                   &
  pstar, tstar, orog

USE parkind1,             ONLY: jpim, jprb       !DrHook
USE yomhook,              ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Input
INTEGER, INTENT(IN)    :: mlevs    ! number of input levels

! Local
INTEGER             :: k, i, j, ii, jj, kmax
REAL(wp)            :: exner_lid, scale_ht_lid, w1, qt, theta_vd, mt
REAL(wp)            :: exner_surf, theta_surf, theta_v_surf, tref

INTEGER(iwp)        ::             & ! safe to use 32 bits integers for levels
  k_lowest(local_new_x,local_new_y)  ! Level where gridpoints above surface


LOGICAL              ::                    &
  mean_mask(local_new_x,local_new_y,mlevs)   ! true if above ground

REAL(wp)   ::                              &
  m_v(local_new_x,local_new_y,mlevs)       & ! water vapour mixing ratio
 ,theta_v(local_new_x,local_new_y,0:mlevs) & ! Virtual potential temperature
 ,exner(local_new_x,local_new_y,0:mlevs+1) & ! Exner pressure on rho levels 
 ,exner_th(local_new_x,local_new_y,mlevs)  & ! Exner pressure on theta levels 
 ,mv_surf(local_new_x,local_new_y)         & ! mv at surface
 ,mt_surf(local_new_x,local_new_y)         & ! mt at surface
 ,pstar_lowest(local_new_x,local_new_y)    & ! Lowest surface pstar (will be
                                             ! mean value over the sea).
 ,tstar_lowest(local_new_x,local_new_y)    & ! Lowest surface tstar (will be
                                             ! mean value over the sea).
 ,orog_lowest(local_new_x,local_new_y)     & ! Lowest surface orogaphy (zero
                                             ! over the sea).
 ,n_lowest(local_new_x,local_new_y)          ! number of gridpoints above 
                                             ! surface

! DR HOOK
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'CRMSTYLE_HSTAT_BALANCE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!--------------------------------------------------------------------------
! Work out pstar_lowest and tstar_lowest for coarse grids
! This will be all_pstar and all_tstar if no orography in area, otherwise
! it is the mean of the points below the lowest level where there are
! points not below the land surface.

IF (l_no_orog) THEN  

  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      pstar_lowest(ii,jj) = all_pstar(ii,jj)
      tstar_lowest(ii,jj) = all_tstar(ii,jj)
      orog_lowest(ii,jj)  = 0.0
      k_lowest(ii,jj) = 1   
    END DO
  END DO

!$OMP PARALLEL DO PRIVATE(k, ii, jj) DEFAULT(SHARED)
  DO k=1,mlevs 
    DO jj = 1,local_new_y
      DO ii = 1,local_new_x
        IF( all_q(ii,jj,k) /= rmdi) THEN
          mean_mask(ii,jj,k) = .TRUE.
        ELSE
          mean_mask(ii,jj,k) = .FALSE.
        END IF    
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE           ! Region has orography  so more complex

  tref = 300.0     ! reference temperature

  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      k_lowest(ii,jj) = 0   
      n_lowest(ii,jj) = 0.0
      pstar_lowest(ii,jj) = 0.0
      tstar_lowest(ii,jj) = 0.0 
      orog_lowest(ii,jj)  = 0.0
    END DO
  END DO

  ! Cannot use openmp on K loop as need to do levels in order from bottom
  ! to get correct k lowest
  DO k=1,mlevs 
    DO jj = 1,local_new_y
      DO ii = 1,local_new_x
        IF( all_q(ii,jj,k) /= rmdi) THEN
          mean_mask(ii,jj,k) = .TRUE.
          IF (k_lowest(ii,jj) == 0) THEN
            k_lowest(ii,jj) = k
          END IF
        ELSE
          mean_mask(ii,jj,k) = .FALSE.
        END IF    
      END DO
    END DO
  END DO

  kmax=0
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      IF (k_lowest(ii,jj) > kmax) THEN
        kmax = k_lowest(ii,jj)
      END IF
    END DO
  END DO

!$OMP PARALLEL DO PRIVATE(k, i, j, ii, jj) DEFAULT(SHARED)
  DO k=1,kmax 
    DO j=1,local_rows
      jj = index_row(j)
      DO i=1,local_row_len
        ii = index_col(i)
        IF ( k_lowest(ii,jj) == k) THEN 
          IF (mask(i,j,k)) THEN 
            pstar_lowest(ii,jj) = pstar_lowest(ii,jj) + (pstar(i,j) -pref)
            tstar_lowest(ii,jj) = tstar_lowest(ii,jj) + (tstar(i,j) -tref)
            orog_lowest(ii,jj)  = orog_lowest(ii,jj)  + orog(i,j)
            n_lowest(ii,jj) = n_lowest(ii,jj)  + 1.0         
          END IF
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! work out mean values
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      pstar_lowest(ii,jj) = pstar_lowest(ii,jj)/n_lowest(ii,jj) +pref
      tstar_lowest(ii,jj) = tstar_lowest(ii,jj)/n_lowest(ii,jj) +tref
      orog_lowest(ii,jj)  = orog_lowest(ii,jj)/n_lowest(ii,jj)
    END DO
  END DO

END IF


! Work out water vapour mixing ratio
! Then work out Virtual and virtual dry potential temperature

!$OMP PARALLEL DO PRIVATE(k, ii, jj, qt, mt, theta_vd ) DEFAULT(SHARED)
DO k=1,mlevs
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      exner(ii,jj,k)  = 0.0     ! Initialise exner array to zero
      IF (mean_mask(ii,jj,k) ) THEN
        qt = all_q(ii,jj,k)+all_qcl(ii,jj,k)+all_qcf(ii,jj,k)+     &
                      all_qrain(ii,jj,k)
        IF (l_qgraup) THEN
          qt = qt +all_qgraup(ii,jj,k) 
        END IF
        m_v(ii,jj,k)  = all_q(ii,jj,k)/(1.0-qt)
        ! total mixing ratio sum, mt
        mt = qt/(1.0-qt)
        ! theta_vd consistent with ENDGame definitions
        theta_vd = all_th(ii,jj,k) * (1.0 + recip_epsilon * m_v(ii,jj,k))
        ! Only correct if no condensate - this is used in idealised
        ! model to derive a balance state
        ! theta_v(ii,jj,k)  = theta_vd / (1.0 + m_v(ii,jj,k))
        ! Correct for condensate present
        theta_v(ii,jj,k)  = theta_vd / (1.0 + mt)
        IF (k==1) THEN
          mv_surf(ii,jj) = m_v(ii,jj,1)
          mt_surf(ii,jj) = mt
        ELSE    ! surface above bottom level
          IF(.NOT.mean_mask(ii,jj,k-1)) THEN   
            mv_surf(ii,jj) = m_v(ii,jj,k)
            mt_surf(ii,jj) = mt
          END IF
        END IF
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Exner at the surface
DO jj = 1,local_new_y
  DO ii = 1,local_new_x
    exner_surf   = (pstar_lowest(ii,jj)/pref)**kappa
    theta_surf   = tstar_lowest(ii,jj)/exner_surf
    theta_v_surf = theta_surf * (1.0 + recip_epsilon * mv_surf(ii,jj)) &
                                     /(1.0 + mt_surf(ii,jj)) 
    ! Exner lowest rho level in mean column
    k = k_lowest(ii,jj)
    exner(ii,jj,k) = exner_surf - g*(h_rho_sea(k)-orog_lowest(ii,jj)) / &
                                                   (cp*theta_v_surf)
  END DO
END DO

! Cannot parallelise loop k using openmp as integral so depends on k order
! working up from surface.
DO k = 1, mlevs-1
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      IF (mean_mask(ii,jj,k+1) .AND. mean_mask(ii,jj,k) ) THEN
       
        exner(ii,jj,k+1) = exner(ii,jj,k) - g*(h_rho_sea(k+1)-h_rho_sea(k)) / &
                                               (cp*theta_v(ii,jj,k))
      END IF
    END DO
  END DO
END DO

! Exner at model_levels+1 is such that, when averaged with the value at
! model_levels, it provides Exner at the upper boundary that is
! obtained by linearly extrapolating log(exner)
DO jj = 1,local_new_y
  DO ii = 1,local_new_x
    scale_ht_lid = (LOG(exner(ii,jj,mlevs-1)) - LOG(exner(ii,jj,mlevs))) /  &
                      (h_rho_sea(mlevs)        - h_rho_sea(mlevs-1))

    exner_lid = exner(ii,jj,mlevs)*EXP(-scale_ht_lid*(h_theta_sea(mlevs) -  &
                                                   h_rho_sea(mlevs)))
    exner(ii,jj,mlevs+1) = 2*exner_lid-exner(ii,jj,mlevs)
  END DO
END DO

! Now want Exner hydrostatic at theta levels, use linear interpolation
!$OMP PARALLEL DO PRIVATE(k, ii, jj, w1 ) DEFAULT(SHARED)
DO k = 1, mlevs-1
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      w1 = (h_rho_sea(k+1) - h_theta_sea(k))/(h_rho_sea(k+1)-h_rho_sea(k))
      exner_th(ii,jj,k) = exner(ii,jj,k)* w1 + exner(ii,jj,k+1)*(1.0-w1)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
k=mlevs
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
     exner_th(ii,jj,k) =0.5*( exner(ii,jj,k) + exner(ii,jj,k+1))
    END DO
  END DO

! Now want P hydrostatic at theta levels - using linear interpolation of exner
! This is consistent with what the UM does.
!$OMP PARALLEL DO PRIVATE(k, ii, jj ) DEFAULT(SHARED)
DO k = 1, mlevs
  DO jj = 1,local_new_y
    DO ii = 1,local_new_x
      IF (mean_mask(ii,jj,k)) THEN
        p_theta_hydro(ii,jj,k) = pref *(exner_th(ii,jj,k)**(1.0/kappa)) 
      ELSE
        p_theta_hydro(ii,jj,k) = rmdi
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE crmstyle_hstat_balance

END MODULE crmstyle_hstat_balance_mod
