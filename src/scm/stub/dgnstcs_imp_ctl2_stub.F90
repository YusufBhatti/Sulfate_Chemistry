! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Stub routine of dgnstcs_imp_ctl2 for full model runs

SUBROUTINE dgnstcs_imp_ctl2                                                    &
  ( row_length, rows, rhc_row_length, rhc_rows                                 &
  , bl_levels, cloud_levels                                                    &
  , rhcrit, combined_cloud, nclds, cumulus, ntml                               &
  , plsp, conv_rain, conv_snow, ls_rain, ls_snow, r_u, r_v                     &
  , zh, zht, bl_type_1, bl_type_2                                              &
  , bl_type_3, bl_type_4, bl_type_5, bl_type_6, bl_type_7, bl_alltypes         &
  , fqt, ftl, t, q, cf, cfl, cff, t_earliest, q_earliest, qcl_earliest         &
  , qcf_earliest, cf_earliest, cfl_earliest, cff_earliest, p_theta_levels      &
  , lwp, iwp, z0m, z0h_eff_gb, z0m_eff_gb, sea_ice_htf                         &
  , taux, tauy, area_cloud_fraction, p_star                                    &
  , cca_2d, rho1, qcl, qcf, surf_ht_flux_sice, rhcpt                           &
  , surf_ht_flux_gb, aerosol, nscmdpkgs, l_scmdiags, bl_diag, sf_diag)


USE atm_fields_bounds_mod, ONLY: vdims_s, vdims, udims_s, pdims_s            &
                               , ScmRowLen,ScmRow
USE nlsizes_namelist_mod,  ONLY: model_levels
USE bl_diags_mod,          ONLY: strnewbldiag
USE sf_diags_mod,          ONLY: strnewsfdiag

IMPLICIT NONE


! Description: (STUB Routine)
!
! Much of what is in here is (regrettably) duplication of code
! in diagnostics_lscld and diagnostics_bl. These routines can't
! be called in the SCM because of explicit STASH calls, but we
! still want the diagnostics.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran 95

! ALL PARAMETERS ARE INTENT IN AS IT IS A STUB ROUTINE

INTEGER :: row_length        
INTEGER :: rows              
INTEGER :: rhc_row_length    
INTEGER :: rhc_rows          
INTEGER :: bl_levels         
INTEGER :: cloud_levels      
INTEGER :: nclds

LOGICAL :: cumulus(row_length,rows)     
INTEGER :: ntml(row_length,rows)         ! Height of diagnosed BL top

REAL :: rhcrit(model_levels)               
REAL :: plsp(row_length,rows)            
REAL :: conv_rain(row_length,rows)           
REAL :: conv_snow(row_length,rows)           
REAL :: ls_rain(row_length,rows)             
REAL :: ls_snow(row_length,rows)             

REAL :: r_u (udims_s%i_start:udims_s%i_end                                   &
            ,udims_s%j_start:udims_s%j_end                                   &
            ,model_levels)                   

REAL :: r_v (vdims_s%i_start:vdims_s%i_end                                   &
            ,vdims_s%j_start:vdims_s%j_end                                   &
            ,model_levels)

REAL :: zh  (row_length,rows)
REAL :: zht (row_length,rows)

! Indicators for diagnosed boundary layers types
! Value set to 1.0 for the following types
REAL :: bl_type_1(row_length,rows)   ! Stable b.l.
REAL :: bl_type_2(row_length,rows)   ! Sc over stable surface layer
REAL :: bl_type_3(row_length,rows)   ! Well mixed b.l.
REAL :: bl_type_4(row_length,rows)   ! Decoupled Sc layer (not over cumulus)
REAL :: bl_type_5(row_length,rows)   ! Decoupled Sc layer over cumulus
REAL :: bl_type_6(row_length,rows)   ! Cumulus capped b.l.
REAL :: bl_type_7(row_length,rows)   ! Shear-dominated b.l
REAL :: bl_alltypes(ScmRowLen,ScmRow)

REAL :: fqt (row_length,rows,bl_levels) 
REAL :: ftl (row_length,rows,bl_levels) 
REAL :: z0m (row_length,rows)           ! Roughness length for mom

REAL :: z0m_eff_gb (row_length,rows)    ! Orographic roughness length
REAL :: z0h_eff_gb (row_length,rows)    ! Roughness length for heat

! Variables needed for PC2 diagnostics
REAL :: t   (row_length,rows,model_levels)                           
REAL :: q   (row_length,rows,model_levels)                             
REAL :: cf  (row_length,rows,model_levels)                            
REAL :: cfl (row_length,rows,model_levels)                           
REAL :: cff (row_length,rows,model_levels)                           

! _earliest arrays contain fields at start of imp_ctl
REAL :: T_earliest   (row_length,rows,model_levels)                  
REAL :: q_earliest   (row_length,rows,model_levels)                    
REAL :: qcl_earliest (row_length,rows,model_levels)                  
REAL :: qcf_earliest (row_length,rows,model_levels)                  

! _earliest values contain the values of temperature, water contents
! and cloud fractions before the boundary layer call.
REAL :: cf_earliest    (row_length,rows,model_levels)                   
REAL :: cfl_earliest   (row_length,rows,model_levels)                  
REAL :: cff_earliest   (row_length,rows,model_levels)                  
REAL :: p_theta_levels (row_length,rows,model_levels)

REAL :: lwp (ScmRowLen,ScmRow)                ! Liquid water path  (kg/m2)
REAL :: iwp (ScmRowLen,ScmRow)                ! Ice water path     (kg/m2)

REAL :: sea_ice_htf  (row_length,rows)       

REAL :: taux (row_length,rows,bl_levels)    
REAL :: tauy (row_length,rows,bl_levels)    

REAL :: area_cloud_fraction (row_length,rows,model_levels) 

REAL :: p_star      (row_length,rows)            
REAL :: u10m        (row_length,rows)              
REAL :: v10m        (row_length,rows)              
REAL :: cca_2d      (row_length,rows)            
REAL :: rho1        (row_length,rows)         ! Air density at level 1/kg/m^3
REAL :: qcl         (row_length,rows,model_levels)    
REAL :: qcf         (row_length,rows,model_levels)    

REAL :: surf_ht_flux_sice (row_length,rows) 
REAL :: surf_ht_flux_gb   (row_length,rows)   

REAL :: rhcpt (rhc_row_length,rhc_rows,model_levels)       

REAL :: aerosol (pdims_s%i_start:pdims_s%i_end                               &
                ,pdims_s%j_start:pdims_s%j_end                               &
                ,model_levels)                          
REAL :: combined_cloud(row_length,rows,model_levels)

! Logicals for SCM Diagnostics packages
INTEGER :: nscmdpkgs               ! No of SCM diagnostics packages

LOGICAL :: l_scmdiags(nscmdpkgs)   ! Logicals for SCM diagnostics packages

! Declaration of new BL diagnostics.
TYPE (strnewbldiag) :: bl_diag
TYPE (strnewsfdiag) :: sf_diag

RETURN

END SUBROUTINE dgnstcs_imp_ctl2

