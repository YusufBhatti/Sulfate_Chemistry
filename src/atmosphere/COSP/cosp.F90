! (c) British Crown Copyright 2008-2018, the Met Office.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
!
!     * Redistributions of source code must retain the above copyright notice, this list
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used
!       to endorse or promote products derived from this software without specific prior written
!       permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!  Description:  Module that applies some quality control to the inputs,
!                calculates extra inputs if needed, and makes iterative
!                calls to COSP, the routine that calls the individual
!                simulators.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: COSP

MODULE MOD_COSP
  USE cosp_types_mod
  USE mod_cosp_utils, ONLY: cosp_ereport
  USE MOD_COSP_SIMULATOR
  USE MOD_COSP_MODIS_SIMULATOR
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP ---------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sghydro,sgradar,sglidar,    &
                isccp,misr,modis,rttov,stradar,stlidar)

  ! Arguments
  INTEGER,INTENT(IN) :: overlap ! overlap type in SCOPS: 1=max,2=rand,3=max/rand
  INTEGER,INTENT(IN) :: Ncolumns
  TYPE(cosp_config),INTENT(IN) :: cfg   ! Configuration options
  TYPE(cosp_vgrid),INTENT(IN) :: vgrid   ! Information on vertical grid of stats
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx
  TYPE(cosp_subgrid),INTENT(INOUT) :: sgx   ! Subgrid info
  TYPE(cosp_sghydro),INTENT(INOUT) :: sghydro ! Subgrid info for hydrometeors
  TYPE(cosp_sgradar),INTENT(INOUT) :: sgradar ! Output from radar simulator
  TYPE(cosp_sglidar),INTENT(INOUT) :: sglidar ! Output from lidar simulator
  TYPE(cosp_isccp),INTENT(INOUT)   :: isccp   ! Output from ISCCP simulator
  TYPE(cosp_misr),INTENT(INOUT)    :: misr    ! Output from MISR simulator
  TYPE(cosp_modis),INTENT(INOUT)   :: modis   ! Output from MODIS simulator
  TYPE(cosp_rttov),INTENT(INOUT)   :: rttov   ! Output from RTTOV
  TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics from lidar

  ! Local variables
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: Niter     ! Number of calls to cosp_simulator
  INTEGER :: i_first,i_last ! First and last gridbox in each iteration
  INTEGER :: i,Ni
  INTEGER,DIMENSION(2) :: ix,iy
  LOGICAL :: reff_zero
  REAL :: maxp,minp

  ! Types used in one iteration
  TYPE(cosp_gridbox) :: gbx_it
  TYPE(cosp_subgrid) :: sgx_it
  TYPE(cosp_vgrid)   :: vgrid_it
  TYPE(cosp_sgradar) :: sgradar_it
  TYPE(cosp_sglidar) :: sglidar_it
  TYPE(cosp_isccp)   :: isccp_it
  TYPE(cosp_modis)   :: modis_it
  TYPE(cosp_misr)    :: misr_it
  TYPE(cosp_rttov)   :: rttov_it
  TYPE(cosp_radarstats) :: stradar_it
  TYPE(cosp_lidarstats) :: stlidar_it
  ! Error handling
  INTEGER :: icode
  CHARACTER(LEN=200) :: cmessage
  CHARACTER(LEN=*),PARAMETER :: routine_name='COSP'


! +++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

! +++++++++ Depth of model layers ++++++++++++
  DO i=1,Nlevels-1
    gbx%dlev(:,i) = gbx%zlev_half(:,i+1) - gbx%zlev_half(:,i)
  ENDDO
  gbx%dlev(:,Nlevels) = 2.0*(gbx%zlev(:,Nlevels) - gbx%zlev_half(:,Nlevels))

! +++++++++ Apply sanity checks to inputs ++++++++++
  CALL cosp_check_input('longitude',gbx%longitude,min_val=0.0,max_val=360.0)
  CALL cosp_check_input('latitude',gbx%latitude,min_val=-90.0,max_val=90.0)
  CALL cosp_check_input('dlev',gbx%dlev,min_val=0.0)
  CALL cosp_check_input('p',gbx%p,min_val=0.0)
  CALL cosp_check_input('ph',gbx%ph,min_val=0.0)
  CALL cosp_check_input('T',gbx%T,min_val=0.0)
  CALL cosp_check_input('q',gbx%q,min_val=0.0)
  CALL cosp_check_input('sh',gbx%sh,min_val=0.0)
  CALL cosp_check_input('dtau_s',gbx%dtau_s,min_val=0.0)
  CALL cosp_check_input('dtau_c',gbx%dtau_c,min_val=0.0)
  CALL cosp_check_input('dem_s',gbx%dem_s,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('dem_c',gbx%dem_c,min_val=0.0,max_val=1.0)
  ! Point information (Npoints)
  CALL cosp_check_input('land',gbx%land,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('psfc',gbx%psfc,min_val=0.0)
  CALL cosp_check_input('sunlit',gbx%sunlit,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('skt',gbx%skt,min_val=0.0)
  ! TOTAL and CONV cloud fraction for SCOPS
  CALL cosp_check_input('tca',gbx%tca,min_val=0.0,max_val=1.0)
  CALL cosp_check_input('cca',gbx%cca,min_val=0.0,max_val=1.0)
  ! Precipitation fluxes on model levels
  CALL cosp_check_input('rain_ls',gbx%rain_ls,min_val=0.0)
  CALL cosp_check_input('rain_cv',gbx%rain_cv,min_val=0.0)
  CALL cosp_check_input('snow_ls',gbx%snow_ls,min_val=0.0)
  CALL cosp_check_input('snow_cv',gbx%snow_cv,min_val=0.0)
  CALL cosp_check_input('grpl_ls',gbx%grpl_ls,min_val=0.0)
  ! Hydrometeors concentration and distribution parameters
  CALL cosp_check_input('mr_hydro',gbx%mr_hydro,min_val=0.0)
  ! Effective radius [m]. (Npoints,Nlevels,Nhydro)
  CALL cosp_check_input('Reff',gbx%Reff,min_val=0.0)
  reff_zero=.TRUE.
  IF (ANY(gbx%Reff > 1.e-8)) THEN
     reff_zero=.FALSE.
      ! reff_zero == .false.
      !     and gbx%use_reff == .true.  Reff use in radar and lidar
      !     and reff_zero    == .false. Reff use in lidar and set to 0 for radar
  ENDIF
  IF ((.NOT. gbx%use_reff) .AND. (reff_zero)) THEN
        ! No Reff in radar. Default in lidar
        gbx%Reff = DEFAULT_LIDAR_REFF
        icode = -9
        cmessage = 'Using default Reff in lidar simulations'
        CALL cosp_ereport(routine_name,cmessage,icode)
  ENDIF

  ! Aerosols concentration and distribution parameters
  CALL cosp_check_input('conc_aero',gbx%conc_aero,min_val=0.0)
  ! Checks for CRM mode
  IF (Ncolumns == 1) THEN
     IF (gbx%use_precipitation_fluxes) THEN
        icode = 9
        cmessage = 'Use of precipitation fluxes not supported '//              &
                   'in CRM mode (Ncolumns=1)'
        CALL cosp_ereport(routine_name,cmessage,icode)
     ENDIF
     IF ((MAXVAL(gbx%dtau_c) > 0.0).OR.(MAXVAL(gbx%dem_c) > 0.0)) THEN
        icode = 9
        cmessage = 'dtau_c > 0.0 or dem_c > 0.0. In CRM mode (Ncolumns=1), '// &
                   'the optical depth (emmisivity) of all clouds must be '//   &
                   'passed through dtau_s (dem_s)'
        CALL cosp_ereport(routine_name,cmessage,icode)
     ENDIF
  ENDIF

  CALL cosp_iter(overlap,cfg,vgrid,gbx,sgx,sghydro,sgradar,sglidar,            &
                   isccp,misr,modis,rttov,stradar,stlidar)

END SUBROUTINE COSP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------- SUBROUTINE COSP_ITER ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_ITER(overlap,cfg,vgrid,gbx,sgx,sghydro,sgradar,sglidar,   &
                     isccp,misr,modis,rttov,stradar,stlidar)

  ! Arguments
  INTEGER,INTENT(IN) :: overlap ! overlap type in SCOPS: 1=max,2=rand,3=max/rand
  TYPE(cosp_config),INTENT(IN) :: cfg   ! Configuration options
  TYPE(cosp_vgrid),INTENT(IN) :: vgrid   ! Information on vertical grid of stats
  TYPE(cosp_gridbox),INTENT(INOUT) :: gbx
  TYPE(cosp_subgrid),INTENT(INOUT) :: sgx   ! Subgrid info
  TYPE(cosp_sghydro),INTENT(INOUT) :: sghydro ! Subgrid info for hydrometeors
  TYPE(cosp_sgradar),INTENT(INOUT) :: sgradar ! Output from radar simulator
  TYPE(cosp_sglidar),INTENT(INOUT) :: sglidar ! Output from lidar simulator
  TYPE(cosp_isccp),INTENT(INOUT)   :: isccp   ! Output from ISCCP simulator
  TYPE(cosp_misr),INTENT(INOUT)    :: misr    ! Output from MISR simulator
  TYPE(cosp_modis),INTENT(INOUT)   :: modis   ! Output from MODIS simulator
  TYPE(cosp_rttov),INTENT(INOUT)   :: rttov   ! Output from RTTOV
  TYPE(cosp_radarstats),INTENT(INOUT) :: stradar ! Summary statistics from radar
  TYPE(cosp_lidarstats),INTENT(INOUT) :: stlidar ! Summary statistics from lidar

  ! Local variables
  INTEGER :: Npoints   ! Number of gridpoints
  INTEGER :: Ncolumns  ! Number of subcolumns
  INTEGER :: Nlevels   ! Number of levels
  INTEGER :: Nhydro    ! Number of hydrometeors
  INTEGER :: i,j,k
  INTEGER :: I_HYDRO
  REAL,DIMENSION(:,:),POINTER :: column_frac_out ! One column of frac_out
  REAL,DIMENSION(:,:),POINTER :: column_prec_out ! One column of prec_frac
  INTEGER :: scops_debug=0    ! non-zero to print out debug info in SCOPS
  REAL,DIMENSION(:, :),ALLOCATABLE :: ls_p_rate,cv_p_rate
                               ! Large-scale and convective precip rates
  REAL,DIMENSION(:,:),ALLOCATABLE :: frac_ls,prec_ls,frac_cv,prec_cv
                              ! Cloud/Precipitation fraction in each model level
                              ! Levels are from SURFACE to TOA
  REAL,DIMENSION(:,:),ALLOCATABLE :: rho ! (Npoints, Nlevels). Air density


  ! +++++++++ Dimensions ++++++++++++
  Npoints  = gbx%Npoints
  Ncolumns = gbx%Ncolumns
  Nlevels  = gbx%Nlevels
  Nhydro   = gbx%Nhydro

  ! +++++++++ Climate/NWP mode ++++++++++
  IF (Ncolumns > 1) THEN
        ! +++++++++ Subgrid sampling ++++++++++
        ! Allocate arrays before calling SCOPS
        ALLOCATE(frac_ls(Npoints,Nlevels),frac_cv(Npoints,Nlevels),            &
                 prec_ls(Npoints,Nlevels),prec_cv(Npoints,Nlevels))
        ALLOCATE(ls_p_rate(Npoints,Nlevels),cv_p_rate(Npoints,Nlevels))
        ! Initialize to zero
        frac_ls=0.0
        prec_ls=0.0
        frac_cv=0.0
        prec_cv=0.0

        ! Flip frac_out upside down, as expected by prec_scops
        DO j=1,Npoints
          sgx%frac_out(j,:,1:Nlevels)  = sgx%frac_out(j,:,Nlevels:1:-1)
        END DO

        ! temporarily use prec_ls/cv to transfer information about
        ! precipitation flux into prec_scops
        IF (gbx%use_precipitation_fluxes) THEN
            ls_p_rate(:,Nlevels:1:-1)=gbx%rain_ls(:,1:Nlevels) +               &
                                      gbx%snow_ls(:,1:Nlevels) +               &
                                      gbx%grpl_ls(:,1:Nlevels)
            cv_p_rate(:,Nlevels:1:-1)=gbx%rain_cv(:,1:Nlevels) +               &
                                      gbx%snow_cv(:,1:Nlevels)
        ELSE
            ls_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_LSRAIN) + &
                                      gbx%mr_hydro(:,1:Nlevels,I_LSGRPL)
            cv_p_rate(:,Nlevels:1:-1)=gbx%mr_hydro(:,1:Nlevels,I_CVRAIN) + &
                                      gbx%mr_hydro(:,1:Nlevels,I_CVSNOW)
        ENDIF

! DEPENDS ON: prec_scops
        CALL prec_scops(Npoints,Nlevels,Ncolumns,ls_p_rate,cv_p_rate,          &
                        sgx%frac_out,sgx%prec_frac)

        ! Precipitation fraction
        DO j=1,Npoints,1
          DO k=1,Nlevels,1
            DO i=1,Ncolumns,1
                IF (sgx%frac_out (j,i,Nlevels+1-k) == I_LSC)                   &
                    frac_ls(j,k)=frac_ls(j,k)+1.
                IF (sgx%frac_out (j,i,Nlevels+1-k) == I_CVC)                   &
                    frac_cv(j,k)=frac_cv(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 1)                     &
                    prec_ls(j,k)=prec_ls(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 2)                     &
                    prec_cv(j,k)=prec_cv(j,k)+1.
                IF (sgx%prec_frac(j,i,Nlevels+1-k) .eq. 3) THEN
                    prec_cv(j,k)=prec_cv(j,k)+1.
                    prec_ls(j,k)=prec_ls(j,k)+1.
                ENDIF
            ENDDO  !i
            frac_ls(j,k)=frac_ls(j,k)/Ncolumns
            frac_cv(j,k)=frac_cv(j,k)/Ncolumns
            prec_ls(j,k)=prec_ls(j,k)/Ncolumns
            prec_cv(j,k)=prec_cv(j,k)/Ncolumns
          ENDDO  !k
        ENDDO  !j

        ! Levels from SURFACE to TOA.
        DO j=1,Npoints
          sgx%frac_out(j,:,1:Nlevels)  = sgx%frac_out(j,:,Nlevels:1:-1)
          sgx%prec_frac(j,:,1:Nlevels) = sgx%prec_frac(j,:,Nlevels:1:-1)
        ENDDO

        ! Deallocate arrays that will no longer be used
        DEALLOCATE(ls_p_rate,cv_p_rate)

        ! Populate the subgrid arrays
        DO k=1,Ncolumns
            !--------- Mixing ratios for clouds and Reff for Clouds and precip -
            column_frac_out => sgx%frac_out(:,k,:)
            WHERE (column_frac_out == I_LSC)     !+++++++++++ LS clouds ++++++++
                sghydro%Np(:,k,:,I_LSCLIQ)     = gbx%Np(:,:,I_LSCLIQ)
                sghydro%Np(:,k,:,I_LSCICE)     = gbx%Np(:,:,I_LSCICE)
                sghydro%Np(:,k,:,I_LSCIAG)     = gbx%Np(:,:,I_LSCIAG)
            ELSEWHERE (column_frac_out == I_CVC) !+++++++++++ CONV clouds ++++++
                sghydro%Np(:,k,:,I_CVCLIQ)     = gbx%Np(:,:,I_CVCLIQ)
                sghydro%Np(:,k,:,I_CVCICE)     = gbx%Np(:,:,I_CVCICE)
            END WHERE
            column_prec_out => sgx%prec_frac(:,k,:)
            WHERE ((column_prec_out == 1) .OR. (column_prec_out == 3))
                !++++ LS precip ++++
                sghydro%Np(:,k,:,I_LSRAIN)     = gbx%Np(:,:,I_LSRAIN)
                sghydro%Np(:,k,:,I_LSGRPL)     = gbx%Np(:,:,I_LSGRPL)

                sghydro%Reff(:,k,:,I_LSRAIN) = gbx%Reff(:,:,I_LSRAIN)
                sghydro%Reff(:,k,:,I_LSGRPL) = gbx%Reff(:,:,I_LSGRPL)
            ELSEWHERE ((column_prec_out == 2) .OR. (column_prec_out == 3))
                !++++ CONV precip ++++
                sghydro%Np(:,k,:,I_CVRAIN)     = gbx%Np(:,:,I_CVRAIN)
                sghydro%Np(:,k,:,I_CVSNOW)     = gbx%Np(:,:,I_CVSNOW)

                sghydro%Reff(:,k,:,I_CVRAIN) = gbx%Reff(:,:,I_CVRAIN)
                sghydro%Reff(:,k,:,I_CVSNOW) = gbx%Reff(:,:,I_CVSNOW)
            END WHERE
            !--------- Precip -------
            IF (.NOT. gbx%use_precipitation_fluxes) THEN
                WHERE ((column_prec_out == 1) .OR. (column_prec_out == 3))
                    !+++++++++++ LS Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_LSRAIN)= gbx%mr_hydro(:,:,I_LSRAIN)
                    sghydro%mr_hydro(:,k,:,I_LSGRPL)= gbx%mr_hydro(:,:,I_LSGRPL)
                ELSEWHERE ((column_prec_out == 2) .OR. (column_prec_out == 3))
                    !+++++++++++ CONV Precipitation ++++++++
                    sghydro%mr_hydro(:,k,:,I_CVRAIN)= gbx%mr_hydro(:,:,I_CVRAIN)
                    sghydro%mr_hydro(:,k,:,I_CVSNOW)= gbx%mr_hydro(:,:,I_CVSNOW)
                END WHERE
            ENDIF
        ENDDO
        ! convert the mixing ratio and precipitation flux from gridbox mean
        ! to the fraction-based values
        DO k=1,Nlevels
            DO j=1,Npoints
                !--------- Clouds -------
                IF (frac_ls(j,k) .ne. 0.) THEN
                    sghydro%mr_hydro(j,:,k,I_LSCLIQ) =                         &
                           sghydro%mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                    sghydro%mr_hydro(j,:,k,I_LSCICE) =                         &
                           sghydro%mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
                    sghydro%mr_hydro(j,:,k,I_LSCIAG) =                         &
                           sghydro%mr_hydro(j,:,k,I_LSCIAG)/frac_ls(j,k)
                ENDIF
                IF (frac_cv(j,k) .ne. 0.) THEN
                    sghydro%mr_hydro(j,:,k,I_CVCLIQ) =                         &
                           sghydro%mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                    sghydro%mr_hydro(j,:,k,I_CVCICE) =                         &
                           sghydro%mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
                ENDIF
                !--------- Precip -------
                IF (gbx%use_precipitation_fluxes) THEN
                    IF (prec_ls(j,k) .ne. 0.) THEN
                        gbx%rain_ls(j,k) = gbx%rain_ls(j,k)/prec_ls(j,k)
                        gbx%snow_ls(j,k) = gbx%snow_ls(j,k)/prec_ls(j,k)
                        gbx%grpl_ls(j,k) = gbx%grpl_ls(j,k)/prec_ls(j,k)
                    ENDIF
                    IF (prec_cv(j,k) .ne. 0.) THEN
                        gbx%rain_cv(j,k) = gbx%rain_cv(j,k)/prec_cv(j,k)
                        gbx%snow_cv(j,k) = gbx%snow_cv(j,k)/prec_cv(j,k)
                    ENDIF
                ELSE
                    IF (prec_ls(j,k) .ne. 0.) THEN
                        sghydro%mr_hydro(j,:,k,I_LSRAIN) =                     &
                               sghydro%mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                        sghydro%mr_hydro(j,:,k,I_LSGRPL) =                     &
                               sghydro%mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                    ENDIF
                    IF (prec_cv(j,k) .ne. 0.) THEN
                        sghydro%mr_hydro(j,:,k,I_CVRAIN) =                     &
                               sghydro%mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                        sghydro%mr_hydro(j,:,k,I_CVSNOW) =                     &
                               sghydro%mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                    ENDIF
                ENDIF
            ENDDO !k
        ENDDO !j
        DEALLOCATE(frac_ls,prec_ls,frac_cv,prec_cv)

        IF (gbx%use_precipitation_fluxes) THEN
          ! Density
          ALLOCATE(rho(Npoints,Nlevels))
          I_HYDRO = I_LSRAIN
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,1.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%rain_ls,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_CVRAIN
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,2.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%rain_cv,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_CVSNOW
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,2.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%snow_cv,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          I_HYDRO = I_LSGRPL
          CALL cosp_precip_mxratio(Npoints,Nlevels,Ncolumns,gbx%p,gbx%T,       &
                  sgx%prec_frac,1.,n_ax(I_HYDRO),n_bx(I_HYDRO),                &
                  alpha_x(I_HYDRO),c_x(I_HYDRO),d_x(I_HYDRO),g_x(I_HYDRO),     &
                  a_x(I_HYDRO),b_x(I_HYDRO),gamma_1(I_HYDRO),gamma_2(I_HYDRO), &
                  gamma_3(I_HYDRO),gamma_4(I_HYDRO),gbx%grpl_ls,               &
                  sghydro%mr_hydro(:,:,:,I_HYDRO),sghydro%Reff(:,:,:,I_HYDRO))
          IF(ALLOCATED(rho)) DEALLOCATE(rho)
        ENDIF
   ! +++++++++ CRM mode ++++++++++
   ELSE
      sghydro%mr_hydro(:,1,:,:) = gbx%mr_hydro
      sghydro%Reff(:,1,:,:) = gbx%Reff
      sghydro%Np(:,1,:,:) = gbx%Np      ! added by Roj with Quickbeam V3.0

      !--------- Clouds -------
      WHERE ((gbx%dtau_s > 0.0))
        ! Subgrid cloud array. Dimensions (Npoints,Ncolumns,Nlevels)
        sgx%frac_out(:,1,:) = 1
      ENDWHERE
   ENDIF ! Ncolumns > 1

   ! +++++++++ Simulator ++++++++++
    CALL cosp_simulator(gbx,sgx,sghydro,cfg,vgrid,sgradar,sglidar,isccp,       &
                        misr,modis,rttov,stradar,stlidar)

END SUBROUTINE COSP_ITER

END MODULE MOD_COSP
