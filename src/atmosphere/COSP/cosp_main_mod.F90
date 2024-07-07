! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_main_mod
  USE planet_constants_mod, ONLY: r, repsilon
  USE cosp_constants_mod, ONLY: i_lscliq, i_lscice, i_lsrain, i_lsciag, &
                                i_cvcliq, i_cvcice, i_cvrain, i_cvsnow, &
                                i_lsgrpl, n_hydro,  parasol_nrefl
  USE cosp_types_mod, cosp_rttov_type => cosp_rttov
  USE cosp_input_mod, ONLY: cosp_csat_vgrid, cosp_nchannels,                   &
                            cosp_nlr, cosp_use_vgrid, cosp_overlap
  USE cosp_mxratio_mod
  USE cosp_diagnostics_mod
  USE mod_cosp
  USE cosp_reff_mod
  USE ereport_mod
  USE mphys_psd_mod,    ONLY: cr, dr, x2i, x3i, x4i, x2ic, x3ic, x4ic,         &
                              x1g, x2g, x4g, ag, bg, cg, dg,                   &
                              ci0, di0, cic0, dic0
  USE mphys_inputs_mod, ONLY: x1r, x2r, ai, bi, l_mcr_qgraup
  USE mphys_constants_mod, ONLY: x1i, x1ic, x4r
  USE atm_fields_bounds_mod, ONLY: tdims, tdims_l

 
USE timer_mod, ONLY: timer
 
  IMPLICIT NONE

! Description:
!   Routine that calls the main COSP routine and passes the outputs
!   to the diagnostics routine/
!
! Method:
!   Control routine that calls the routines that computes effective radii
!   of the requested hydrometeors, allocate the output types,
!   call the COSP routines, and copy the diagnotics to STASH.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: COSP
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
  SUBROUTINE cosp_main(at_extremity,Lcosp_run,row_length,rows,                 &
                  n_rows,model_levels,qrain,qgraup,cosp_frac_agg,              &
                  cosp_cfg,cosp_gbx,cosp_sgx,cosp_sgh,                         &
      STASHwork2)

  USE errormessagelength_mod, ONLY: errormessagelength
  USE model_domain_mod, ONLY: model_type, mt_single_column
  USE nlstcall_mod, ONLY: ltimer

  IMPLICIT NONE


!-Input arguments
  LOGICAL,INTENT(IN) :: at_extremity(4) ! PE at N,S,E or W of the grid
  LOGICAL,INTENT(IN) :: Lcosp_run ! COSP run or only inputs requested?
  INTEGER,INTENT(IN) :: row_length,rows,n_rows,model_levels
  REAL,INTENT (IN) :: qrain(tdims_l%i_start:tdims_l%i_end, &
                            tdims_l%j_start:tdims_l%j_end, &
                            tdims_l%k_start:tdims_l%k_end)
  REAL,INTENT(IN) :: qgraup(tdims_l%i_start:tdims_l%i_end, &
                            tdims_l%j_start:tdims_l%j_end, &
                            tdims_l%k_start:tdims_l%k_end)
  REAL,INTENT(IN) :: cosp_frac_agg(row_length*rows, 1:tdims%k_end)
  TYPE(cosp_config),INTENT(IN)     :: cosp_cfg  ! Configuration options
!-Input/output arguments
  TYPE(cosp_gridbox),INTENT(INOUT)  :: cosp_gbx      ! Gridbox-mean inputs
  TYPE(cosp_subgrid),INTENT(INOUT)  :: cosp_sgx      ! Subgrid variables
  TYPE(cosp_sghydro),INTENT(INOUT)  :: cosp_sgh      ! Subgrid hydrometeors
  REAL,INTENT(INOUT)                :: STASHwork2(*) ! STASH workspace
!-Local variables
  TYPE(cosp_sgradar)    :: cosp_sgrad   ! Output from radar simulator
  TYPE(cosp_sglidar)    :: cosp_sglid   ! Output from lidar simulator
  TYPE(cosp_isccp)      :: cosp_is      ! Output from ISCCP simulator
  TYPE(cosp_misr)       :: cosp_ms      ! Output from MISR simulator
  TYPE(cosp_modis)      :: cosp_mds     ! Output from MODIS simulator
  TYPE(cosp_rttov_type) :: cosp_rttovl  ! Output from RTTOV
  TYPE(cosp_vgrid)      :: cosp_vg      ! Vertical grid of stats
  TYPE(cosp_radarstats) :: cosp_stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats) :: cosp_stlidar ! Summary statistics from lidar
  INTEGER :: icode, i, j, k, i_hydro
  REAL, ALLOCATABLE :: aux(:,:)

! Exponent that controls the temperature dependence of the intercept
! rainfall of the PSD (0.0 -> no dependence)
  REAL, PARAMETER :: X3R = 0.0
  REAL, PARAMETER :: X3G = 0.0
! Exponent of the normalised density in the terminal fall speed (UMDP26)
  REAL, PARAMETER :: GX = 0.4
! Coefficients for the density distribution of rainfall (Homogeneous liquid
! spheres)
  REAL, PARAMETER :: AR = 523.6
  REAL, PARAMETER :: BR = 3.0

  LOGICAL, PARAMETER :: no_precip_flux = .FALSE.

  CHARACTER(LEN=9), PARAMETER :: RoutineName='COSP_MAIN'
  CHARACTER(LEN=errormessagelength) :: cmessage = ' '

  IF (cosp_cfg%Lrttov_sim) THEN
     icode = 1
     cmessage = " COSP RTTOV simulator not functional at this version."// &
                " Contact the code owner for details."
     CALL Ereport(RoutineName,icode,cmessage)
  END IF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Populate input structure
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute effective radius consistent with microphysics
! The mixing ratios are gathered from several places:
! lw_rad for cloud water mixing ratios
! microphys_ctl for LS precip fluxes
! ni_conv_ctl for convective precip
! The LS and CONV liquid effective radius are gathered from lw_rad
! Sunlit is filled in glue_rad
  ALLOCATE(aux(cosp_gbx%npoints,cosp_gbx%nlevels))
! Density
  aux = cosp_gbx%p/(r*cosp_gbx%T*(1.0+((1.0 - repsilon)/repsilon)*cosp_gbx%sh -&
        cosp_gbx%mr_hydro(:,:,I_LSCLIQ) - cosp_gbx%mr_hydro(:,:,I_LSCICE) -    &
        cosp_gbx%mr_hydro(:,:,I_CVCLIQ) - cosp_gbx%mr_hydro(:,:,I_CVCICE)))
! Ice cristals and aggregates Reff. For convection, assume all aggregates.
! All mass comes in lscice from lw_rad, so
! only ice is used. CALIPSO uses aggregates, so arrays for aggregates
! replicate the LSICE.
! Protected by Lcosp_run because the subgrid types are not allocated if COSP is
! not run (diagnostic mode).
  IF (Lcosp_run) THEN
    DO j=1,cosp_gbx%ncolumns
      CALL cosp_reff(no_precip_flux,X1I,X2I,X3I,X4I,AI,BI,CI0,DI0,GX,          &
                cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                  &
                cosp_sgh%mr_hydro(:,j,:,I_LSCICE),cosp_sgh%Reff(:,j,:,I_LSCICE))
      CALL cosp_reff(no_precip_flux,X1I,X2I,X3I,X4I,AI,BI,CI0,DI0,GX,          &
                cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                  &
                cosp_sgh%mr_hydro(:,j,:,I_CVCICE),cosp_sgh%Reff(:,j,:,I_CVCICE))
    END DO
    cosp_sgh%mr_hydro(:,:,:,I_LSCIAG) = cosp_sgh%mr_hydro(:,:,:,I_LSCICE)
    cosp_sgh%Reff(:,:,:,I_LSCIAG)     = cosp_sgh%Reff(:,:,:,I_LSCICE)
  END IF

! Rainfall mixing ratios and effective radii.
  IF (Lcosp_run .AND. (.NOT. cosp_gbx%use_precipitation_fluxes)) THEN
!   LS rain mixing ratios. Snow diagnosed from the procnostic ice, so the
!   snow mass is already contained in the cloud ice variable.
    i_hydro = i_lsrain
    cosp_gbx%mr_hydro(:,:,i_lsrain) =                                          &
             RESHAPE(qrain(tdims%i_start:tdims%i_end,                          &
                     tdims%j_start:tdims%j_end, 1:tdims%k_end),                &
                    (/cosp_gbx%npoints,model_levels/))
    CALL cosp_reff(no_precip_flux,X1R,X2R,X3R,X4R,AR,BR,CR,DR,GX,              &
                 cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,                 &
                 cosp_gbx%mr_hydro(:,:,I_LSRAIN),cosp_gbx%Reff(:,:,I_LSRAIN))
!   Convective precip is always expressed in terms of fluxes, so we need to
!   calculate convective precip mixing ratios and effective radii here
!   when precip fluxes are not used. When use_precipitation_fluxes is true,
!   then the calculations are done within COSP along the large-scale precip.
    i_hydro = i_cvrain
    CALL cosp_mxratio(cosp_gbx%npoints,model_levels,cosp_gbx%p,cosp_gbx%t,     &
                  n_ax(i_hydro),n_bx(i_hydro),alpha_x(i_hydro),c_x(i_hydro),   &
                  d_x(i_hydro),g_x(i_hydro),a_x(i_hydro),b_x(i_hydro),         &
                  gamma_1(i_hydro),gamma_2(i_hydro),gamma_3(i_hydro),          &
                  gamma_4(i_hydro),cosp_gbx%rain_cv,                           &
                  cosp_gbx%mr_hydro(:,:,i_hydro),cosp_gbx%reff(:,:,i_hydro))
    i_hydro = i_cvsnow
    CALL cosp_mxratio(cosp_gbx%npoints,model_levels,cosp_gbx%p,cosp_gbx%t,     &
                  n_ax(i_hydro),n_bx(i_hydro),alpha_x(i_hydro),c_x(i_hydro),   &
                  d_x(i_hydro),g_x(i_hydro),a_x(i_hydro),b_x(i_hydro),         &
                  gamma_1(i_hydro),gamma_2(i_hydro),gamma_3(i_hydro),          &
                  gamma_4(i_hydro),cosp_gbx%snow_cv,                           &
                  cosp_gbx%mr_hydro(:,:,i_hydro),cosp_gbx%reff(:,:,i_hydro))
!   Graupel
    IF (l_mcr_qgraup) THEN
      cosp_gbx%mr_hydro(:,:,i_lsgrpl) =                                        &
             RESHAPE(qgraup(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end, 1:tdims%k_end),                &
                    (/cosp_gbx%npoints,model_levels/))
      CALL cosp_reff(no_precip_flux,X1G,X2G,X3G,X4G,AG,BG,CG,DG,GX,            &
                     cosp_gbx%npoints,model_levels,cosp_gbx%T,aux,             &
                     cosp_gbx%grpl_ls,cosp_gbx%Reff(:,:,I_LSGRPL))
    END IF
  END IF
  DEALLOCATE(aux)

! Frac out (large-scale=1, convective=2)
  IF (Lcosp_run) THEN
    DO k=1,cosp_gbx%nlevels
      DO j=1,cosp_gbx%ncolumns
        DO i=1,cosp_gbx%npoints
          IF ((cosp_sgh%mr_hydro(i,j,k,I_CVCLIQ) > 0.0).OR.                    &
              (cosp_sgh%mr_hydro(i,j,k,I_CVCICE) > 0.0)) THEN
            cosp_sgx%frac_out(i,j,k) = I_CVC
          END IF
          IF ((cosp_sgh%mr_hydro(i,j,k,I_LSCLIQ) > 0.0).OR.                    &
              (cosp_sgh%mr_hydro(i,j,k,I_LSCICE) > 0.0)) THEN
            cosp_sgx%frac_out(i,j,k) = I_LSC
          END IF
        END DO
      END DO
    END DO
  END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Define new vertical grid
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL construct_cosp_vgrid(cosp_gbx,cosp_nlr,cosp_use_vgrid,                  &
                            cosp_csat_vgrid,cosp_vg)

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Allocate memory for other types
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  IF (Lcosp_run) THEN
    CALL construct_cosp_sgradar(cosp_cfg,cosp_gbx%npoints,cosp_gbx%ncolumns,   &
                    model_levels,N_HYDRO,cosp_sgrad)
    CALL construct_cosp_radarstats(cosp_cfg,cosp_gbx%npoints,cosp_gbx%ncolumns,&
                    cosp_vg%Nlvgrid,N_HYDRO,cosp_stradar)
    CALL construct_cosp_sglidar(cosp_cfg,cosp_gbx%npoints,cosp_gbx%ncolumns,   &
                    model_levels,N_HYDRO,PARASOL_NREFL,cosp_sglid)
    CALL construct_cosp_lidarstats(cosp_cfg,cosp_gbx%npoints,cosp_gbx%ncolumns,&
                    cosp_vg%Nlvgrid,N_HYDRO,PARASOL_NREFL,cosp_stlidar)
    CALL construct_cosp_isccp(cosp_cfg,cosp_gbx%npoints,cosp_gbx%ncolumns,     &
                    model_levels,cosp_is)
    CALL construct_cosp_misr(cosp_cfg,cosp_gbx%npoints,cosp_ms)
    CALL construct_cosp_modis(cosp_cfg,cosp_gbx%npoints,cosp_mds)
    CALL construct_cosp_rttov(cosp_cfg,cosp_gbx%npoints,cosp_nchannels,        &
                    cosp_rttovl)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Call simulator
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF (Ltimer) CALL timer ('COSP',5)
    CALL cosp(cosp_overlap,cosp_gbx%ncolumns,cosp_cfg,cosp_vg,cosp_gbx,        &
                  cosp_sgx,cosp_sgh,cosp_sgrad,cosp_sglid,cosp_is,cosp_ms,     &
                  cosp_mds,cosp_rttovl,cosp_stradar,cosp_stlidar)
    IF (Ltimer) CALL timer ('COSP',6)
  ELSE ! Minumum allocation
    CALL construct_cosp_subgrid(1, 1, 1, cosp_sgx)
    CALL construct_cosp_sgradar(cosp_cfg,1, 1, 1,N_HYDRO,cosp_sgrad)
    CALL construct_cosp_radarstats(cosp_cfg,1, 1, 1, 1,cosp_stradar)
    CALL construct_cosp_sglidar(cosp_cfg,1, 1, 1, 1, 1,cosp_sglid)
    CALL construct_cosp_lidarstats(cosp_cfg,1, 1, 1, 1, 1,cosp_stlidar)
    CALL construct_cosp_isccp(cosp_cfg,1, 1, 1,cosp_is)
    CALL construct_cosp_misr(cosp_cfg, 1,cosp_ms)
    CALL construct_cosp_modis(cosp_cfg, 1,cosp_mds)
    CALL construct_cosp_rttov(cosp_cfg, 1, 1,cosp_rttovl)
    CALL construct_cosp_sghydro(1, 1, 1, 1, cosp_sgh)
  END IF ! Lcosp_run

IF (model_type /= mt_single_column) THEN
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Copy COSP diagnostics to stash
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL cosp_diagnostics(at_extremity,row_length,rows,                          &
         n_rows,model_levels,Lcosp_run,cosp_cfg,cosp_vg,cosp_gbx,              &
         cosp_sgx,cosp_sgrad,cosp_sglid,cosp_is,cosp_ms,cosp_mds,cosp_rttovl,  &
         cosp_stradar,cosp_stlidar,STASHwork2)
END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Deallocate memory in derived types.
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL free_cosp_rttov(cosp_rttovl)
  CALL free_cosp_modis(cosp_mds)
  CALL free_cosp_misr(cosp_ms)
  CALL free_cosp_isccp(cosp_is)
  CALL free_cosp_lidarstats(cosp_stlidar)
  CALL free_cosp_sglidar(cosp_sglid)
  CALL free_cosp_radarstats(cosp_stradar)
  CALL free_cosp_sgradar(cosp_sgrad)
  CALL free_cosp_vgrid(cosp_vg)

  RETURN
  END SUBROUTINE cosp_main
END MODULE cosp_main_mod
