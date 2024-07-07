! Copyright (c) 2009, Centre National de la Recherche Scientifique
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
!     * Neither the name of the LMD/IPSL/CNRS/UPMC nor the names of its
!       contributors may be used to endorse or promote products derived from this software without
!       specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: COSP

!------------------------------------------------------------------------------------
! Authors: Sandrine Bony and Helene Chepfer (LMD/IPSL, CNRS, UPMC, France).
!------------------------------------------------------------------------------------
MODULE MOD_LMD_IPSL_STATS
  USE mod_cosp_utils, ONLY: cosp_ereport
  USE MOD_LLNL_STATS
  USE cosp_input_mod, ONLY: cosp_sr_cloud
  IMPLICIT NONE

CONTAINS
      SUBROUTINE diag_lidar(npoints,ncol,llm,max_bin,nrefl,tmp,pnorm,          &
                            pnorm_perp,pmol,refl,land,pplay,undef,             &
                            ok_lidar_cfad,cfad2,srbval,ncat,lidarcld,          &
                            lidarcldphase,cldlayer,cldlayerphase,lidarcldtmp,  &
                            parasolrefl)
      IMPLICIT NONE
!
! -----------------------------------------------------------------------------------
! Lidar outputs :
!
! Diagnose cloud fraction (3D cloud fraction + 
! low/middle/high/total cloud fraction)
! and phase cloud fraction (3D, low/mid/high/total and 3D temperature)
! from the lidar signals (ATB, ATBperp and molecular ATB).
!      +
! Compute CFADs of lidar scattering ratio SR and of depolarization index
!
! Authors: Sandrine Bony, Helene Chepfer, and Gregory Cesana
!          (LMD/IPSL, CNRS, UPMC, France).
!

! c inputs :
      integer npoints
      integer ncol
      integer llm
      integer max_bin               ! nb of bins for SR CFADs
      integer ncat                  ! nb of cloud layer types (low,mid,high,total)
      integer nrefl                 ! nb of solar zenith angles for parasol reflectances

      real undef                    ! undefined value
      real pnorm(npoints,ncol,llm)  ! lidar ATB
      real pmol(npoints,llm)        ! molecular ATB
      real land(npoints)            ! Landmask [0 - Ocean, 1 - Land]
      real pplay(npoints,llm)       ! pressure on model levels (Pa)
      LOGICAL ok_lidar_cfad         ! true if lidar CFAD diagnostics
                                    ! need to be computed
      real refl(npoints,ncol,nrefl) ! subgrid parasol reflectance ! parasol
      REAL tmp(npoints,llm)         ! temp at each levels
      REAL pnorm_perp(npoints,ncol,llm)  ! lidar perpendicular ATB

! c outputs :
      real lidarcld(npoints,llm)     ! 3D "lidar" cloud fraction
      REAL sub(npoints,llm)     ! 3D "lidar" indice
      REAL cldlayer(npoints,ncat)    ! "lidar" cloud layer fraction
                                     ! (low, mid, high, total)

      real cfad2(npoints,max_bin,llm) ! CFADs of SR
      real srbval(max_bin)           ! SR bins in CFADs
      real parasolrefl(npoints,nrefl)! grid-averaged parasol reflectance

! c threshold for cloud detection :
      real S_clr
      parameter (S_clr = 1.2)
      real S_cld
      real S_att
      parameter (S_att = 0.01)

! c local variables :
      INTEGER ic,k,i,j
      real x3d(npoints,ncol,llm)
      real x3d_c(npoints,llm),pnorm_c(npoints,llm)
      real xmax

! Output variables
      INTEGER,PARAMETER :: nphase = 6
        ! nb of cloud layer phase types
        ! (ice,liquid,undefined,false ice,false liquid,Percent of ice)
      REAL lidarcldphase(npoints,llm,nphase)
        ! 3D "lidar" phase cloud fraction
      REAL lidarcldtmp(npoints,40,5)
        ! 3D "lidar" phase cloud fraction as a function of temp
      REAL cldlayerphase(npoints,ncat,nphase)
        ! "lidar" phase low mid high cloud fraction

! SR detection threshold
      REAL, PARAMETER  ::  S_cld_att = 30. ! New threshold for undefine cloud
                                           ! phase detection


!
! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
!
      S_cld = cosp_sr_cloud

!  Should be modified in future version
      xmax=undef-1.0

! c -------------------------------------------------------
! c 1- Lidar scattering ratio :
! c -------------------------------------------------------

      do ic = 1, ncol
        pnorm_c = pnorm(:,ic,:)
        where ((pnorm_c.lt.xmax) .and. (pmol.lt.xmax) .and. (pmol.gt. 0.0 ))
            x3d_c = pnorm_c/pmol
        elsewhere
            x3d_c = undef
        end where
        x3d(:,ic,:) = x3d_c
      enddo

! c -------------------------------------------------------
! c 2- Diagnose cloud fractions (3D, low, middle, high, total)
! c from subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------

    CALL COSP_CLDFRAC(npoints,ncol,llm,ncat,nphase,tmp,x3d,pnorm,pnorm_perp,   &
                      pplay, S_att,S_cld,S_cld_att,undef,lidarcld,cldlayer,    &
                      lidarcldphase,sub,cldlayerphase,lidarcldtmp)

! c -------------------------------------------------------
! c 3- CFADs
! c -------------------------------------------------------
      if (ok_lidar_cfad) then
!
! c CFADs of subgrid-scale lidar scattering ratios :
! c -------------------------------------------------------
      CALL COSP_CFAD_SR(npoints,ncol,llm,max_bin,undef, &
                 x3d, &
                 S_att,S_clr,xmax,cfad2,srbval)

      endif   ! ok_lidar_cfad
! c -------------------------------------------------------

! c -------------------------------------------------------
! c 4- Compute grid-box averaged Parasol reflectances
! c -------------------------------------------------------

      parasolrefl(:,:) = 0.0

      do k = 1, nrefl
       do ic = 1, ncol
         parasolrefl(:,k) = parasolrefl(:,k) + refl(:,ic,k)
       enddo
      enddo

      do k = 1, nrefl
        parasolrefl(:,k) = parasolrefl(:,k) / float(ncol)
! if land=1 -> parasolrefl=undef
! if land=0 -> parasolrefl=parasolrefl
        parasolrefl(:,k) = parasolrefl(:,k) * MAX(1.0-land(:),0.0) &
                           + (1.0 - MAX(1.0-land(:),0.0))*undef
      enddo

      RETURN
      END SUBROUTINE diag_lidar


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- FUNCTION COSP_CFAD_SR ------------------------
! Author: Sandrine Bony (LMD/IPSL, CNRS, Paris)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE COSP_CFAD_SR(Npoints,Ncolumns,Nlevels,Nbins,undef, &
                      x,S_att,S_clr,xmax,cfad,srbval)
      IMPLICIT NONE

!--- Input arguments
! Npoints: Number of horizontal points
! Ncolumns: Number of subcolumns
! Nlevels: Number of levels
! Nbins: Number of x axis bins
! xmax: maximum value allowed for x
! S_att: Threshold for full attenuation
! S_clr: Threshold for clear-sky layer
!
!--- Input-Outout arguments
! x: variable to process (Npoints,Ncolumns,Nlevels), mofified where saturation occurs
!
! -- Output arguments
! srbval : values of the histogram bins
! cfad: 2D histogram on each horizontal point

! Input arguments
      integer Npoints,Ncolumns,Nlevels,Nbins
      real xmax,S_att,S_clr,undef
! Input-output arguments
      real x(Npoints,Ncolumns,Nlevels)
! Output :
      real cfad(Npoints,Nbins,Nlevels)
      real srbval(Nbins)
! Local variables
      integer i, j, k, ib
      real srbval_ext(0:Nbins)

! c -------------------------------------------------------
! c 0- Initializations
! c -------------------------------------------------------
      if ( Nbins .lt. 6) return

      srbval(1) =  S_att
      srbval(2) =  S_clr
      srbval(3) =  3.0
      srbval(4) =  5.0
      srbval(5) =  7.0
      srbval(6) = 10.0
      do i = 7, MIN(10,Nbins)
       srbval(i) = srbval(i-1) + 5.0
      enddo
      DO i = 11, MIN(13,Nbins)
       srbval(i) = srbval(i-1) + 10.0
      enddo
      srbval(MIN(14,Nbins)) = 80.0
      srbval(Nbins) = xmax
      cfad(:,:,:) = 0.0

      srbval_ext(1:Nbins) = srbval
      srbval_ext(0) = -1.0
! c -------------------------------------------------------
! c c- Compute CFAD
! c -------------------------------------------------------

      do j = 1, Nlevels
         do ib = 1, Nbins
            do k = 1, Ncolumns
               do i = 1, Npoints
                  if (x(i,k,j) /= undef) then
                     if ((x(i,k,j).gt.srbval_ext(ib-1)).and.(x(i,k,j).le.srbval_ext(ib))) &
                          cfad(i,ib,j) = cfad(i,ib,j) + 1.0
                  else
                     cfad(i,ib,j) = undef
                  endif
               enddo
            enddo
         enddo
      enddo

      where (cfad .ne. undef)  cfad = cfad / float(Ncolumns)

! c -------------------------------------------------------
      RETURN
      END SUBROUTINE COSP_CFAD_SR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------------- SUBROUTINE COSP_CLDFRAC -------------------
! c Purpose: Cloud fraction diagnosed from lidar measurements
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE COSP_CLDFRAC(Npoints,Ncolumns,Nlevels,Ncat,Nphase, &
                  tmp,x,ATB,ATBperp,pplay,S_att,S_cld,S_cld_att,undef,lidarcld,&
                  cldlayer,lidarcldphase,nsub,cldlayerphase,lidarcldtemp)


      IMPLICIT NONE
! Input arguments
      integer Npoints,Ncolumns,Nlevels,Ncat
      real x(Npoints,Ncolumns,Nlevels)


! Local parameters
      INTEGER nphase ! nb of cloud layer phase types
                     ! (ice,liquid,undefined,false ice,
                     ! false liquid,Percent of ice)
      INTEGER,PARAMETER  ::  Ntemp=40 ! indice of the temperature vector
      INTEGER ip, k, iz, ic, ncol, nlev, i, itemp  ! loop indice
      REAL  S_cld_att ! New threshold for undefine cloud phase detection (SR=30)
      INTEGER toplvlsat  ! level of the first cloud with SR>30
      REAL alpha50, beta50, gamma50, delta50, epsilon50, zeta50
        ! Polynomial Coef of the phase discrimination line

! Input variables
      REAL tmp(Npoints,Nlevels) ! temperature
      REAL ATB(Npoints,Ncolumns,Nlevels) ! 3D Attenuated backscatter
      REAL ATBperp(Npoints,Ncolumns,Nlevels) ! 3D perpendicular attenuated
                                             ! backscatter
      real pplay(Npoints,Nlevels)
      real S_att,S_cld
      real undef

! Output variables
      REAL lidarcldtemp(Npoints,Ntemp,5)
        ! 3D Temperature 1=tot,2=ice,3=liq,4=undef,5=ice/ice+liq
      REAL tempmod(Ntemp+1) ! temperature bins
      REAL lidarcldphase(Npoints,Nlevels,Nphase) ! 3D cloud phase fraction
      REAL cldlayerphase(Npoints,Ncat,Nphase) 
        ! low, middle, high, total cloud fractions
        ! for ice liquid and undefine phase
      real lidarcld(Npoints,Nlevels) ! 3D cloud fraction
      real cldlayer(Npoints,Ncat)    ! low, middle, high, total cloud fractions

! Local variables
      REAL tmpi(Npoints,Ncolumns,Nlevels) ! temperature of ice cld
      REAL tmpl(Npoints,Ncolumns,Nlevels) ! temperature of liquid cld
      REAL tmpu(Npoints,Ncolumns,Nlevels) ! temperature of undef cld

      REAL checktemp, ATBperp_tmp ! temporary variable
      REAL checkcldlayerphase, checkcldlayerphase2 ! temporary variable
      REAL sumlidarcldtemp(Npoints,Ntemp) ! temporary variable

      REAL cldlayphase(Npoints,Ncolumns,Ncat,Nphase)
        ! subgrided low mid high phase cloud fraction
      REAL cldlayerphasetmp(Npoints,Ncat) ! temporary variable
      REAL cldlayerphasesum(Npoints,Ncat) ! temporary variable
      REAL lidarcldtempind(Npoints,Ntemp) ! 3D Temperature indice
      REAL lidarcldphasetmp(Npoints,Nlevels)
        ! 3D sum of ice and liquid cloud occurences


! Local variables
      real p1
      real cldy(Npoints,Ncolumns,Nlevels)
      real srok(Npoints,Ncolumns,Nlevels)
      real cldlay(Npoints,Ncolumns,Ncat)
      real nsublay(Npoints,Ncolumns,Ncat), nsublayer(Npoints,Ncat)
      real nsub(Npoints,Nlevels)
! Error handling
      INTEGER :: icode
      CHARACTER(LEN=100) :: cmessage
      CHARACTER(LEN=*),PARAMETER :: routine_name='COSP_CLDFRAC'
#ifdef SYS_SX
      real cldlay1(Npoints,Ncolumns)
      real cldlay2(Npoints,Ncolumns)
      real cldlay3(Npoints,Ncolumns)
      real nsublay1(Npoints,Ncolumns)
      real nsublay2(Npoints,Ncolumns)
      real nsublay3(Npoints,Ncolumns)
#endif

! ---------------------------------------------------------------
! 1- initialization
! ---------------------------------------------------------------

      if ( Ncat .ne. 4 ) then
        icode = 9
        WRITE(cmessage,*) 'Ncat must be 4, not ',Ncat
        CALL cosp_ereport(routine_name,cmessage,icode)
      endif

      lidarcld = 0.0
      nsub = 0.0
      cldlay = 0.0
      nsublay = 0.0

      ATBperp_tmp = 0.
      lidarcldphase(:,:,:) = 0.
      cldlayphase(:,:,:,:) = 0.
      cldlayerphase(:,:,:) = 0.
      tmpi(:,:,:) = 0.
      tmpl(:,:,:) = 0.
      tmpu(:,:,:) = 0.
      cldlayerphasesum(:,:) = 0.
      lidarcldtemp(:,:,:) = 0.
      lidarcldtempind(:,:) = 0.
      sumlidarcldtemp(:,:) = 0.
      toplvlsat=0
      lidarcldphasetmp(:,:) = 0.

! temperature bins
      tempmod=(/-273.15,-90.,-87.,-84.,-81.,-78.,-75.,-72.,-69.,-66.,          &
                -63.,-60.,-57.,-54.,-51.,-48.,-45.,-42.,-39.,-36.,             &
                -33.,-30.,-27.,-24.,-21.,-18.,-15.,-12.,-9.,-6.,               &
                -3.,0.,3.,6.,9.,12.,15.,18.,21.,24.,200. /)

! convert C to K
      tempmod=tempmod+273.15

! Polynomial coefficient of the phase discrimination line used to separate
! liquid from ice (Cesana and Chepfer, JGR, 2013)
! ATBperp = ATB^5*alpha50 + ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50 +
!           + ATB*epsilon50 + zeta50
      alpha50   = 9.0322e+15
      beta50    = -2.1358e+12
      gamma50   = 173.3963e06
      delta50   = -3.9514e03
      epsilon50 = 0.2559
      zeta50    = -9.4776e-07


! ---------------------------------------------------------------
! 2- Cloud detection
! ---------------------------------------------------------------

      do k = 1, Nlevels

! cloud detection at subgrid-scale:
         where ( (x(:,:,k).gt.S_cld) .and. (x(:,:,k).ne. undef) )
           cldy(:,:,k)=1.0
         elsewhere
           cldy(:,:,k)=0.0
         endwhere

! number of usefull sub-columns:
         where ( (x(:,:,k).gt.S_att) .and. (x(:,:,k).ne. undef)  )
           srok(:,:,k)=1.0
         elsewhere
           srok(:,:,k)=0.0
         endwhere

      enddo ! k

! ---------------------------------------------------------------
! 3- grid-box 3D cloud fraction and layered cloud fractions (ISCCP pressure
! categories) :
! ---------------------------------------------------------------
      lidarcld = 0.0
      nsub = 0.0
#ifdef SYS_SX
!! XXX: Use cldlay[1-3] and nsublay[1-3] to avoid bank-conflicts.
      cldlay1 = 0.0
      cldlay2 = 0.0
      cldlay3 = 0.0
      cldlay(:,:,4) = 0.0 !! XXX: Ncat == 4
      nsublay1 = 0.0
      nsublay2 = 0.0
      nsublay3 = 0.0
      nsublay(:,:,4) = 0.0
      do k = Nlevels, 1, -1
       do ic = 1, Ncolumns
        do ip = 1, Npoints

         IF (srok(ip,ic,k).gt.0.) THEN
           ! Computation of the cloud fraction as a function of the temperature
           ! instead of height, for ice,liquid and all clouds
           DO itemp=1,Ntemp
             IF ((tmp(ip,k).ge.tempmod(itemp)).AND.                            &
                 (tmp(ip,k).lt.tempmod(itemp+1)) ) THEN
               lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1.
             END IF
           END DO
         END IF

         IF (cldy(ip,ic,k).eq.1.) THEN
           DO itemp=1,Ntemp
             IF ((tmp(ip,k).ge.tempmod(itemp)).AND.                            &
                 (tmp(ip,k).lt.tempmod(itemp+1)) ) THEN
               lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1.
             END IF
           END DO
         END IF

         p1 = pplay(ip,k)

         if ( p1.gt.0. .and. p1.lt.(440.*100.)) then ! high clouds
            cldlay3(ip,ic) = MAX(cldlay3(ip,ic), cldy(ip,ic,k))
            nsublay3(ip,ic) = MAX(nsublay3(ip,ic), srok(ip,ic,k))
         else if(p1.ge.(440.*100.) .and. p1.lt.(680.*100.)) then  ! mid clouds
            cldlay2(ip,ic) = MAX(cldlay2(ip,ic), cldy(ip,ic,k))
            nsublay2(ip,ic) = MAX(nsublay2(ip,ic), srok(ip,ic,k))
         else
            cldlay1(ip,ic) = MAX(cldlay1(ip,ic), cldy(ip,ic,k))
            nsublay1(ip,ic) = MAX(nsublay1(ip,ic), srok(ip,ic,k))
         endif

         cldlay(ip,ic,4) = MAX(cldlay(ip,ic,4), cldy(ip,ic,k))
         lidarcld(ip,k)=lidarcld(ip,k) + cldy(ip,ic,k)
         nsublay(ip,ic,4) = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
         nsub(ip,k)=nsub(ip,k) + srok(ip,ic,k)
        enddo
       enddo
      enddo
      cldlay(:,:,1) = cldlay1
      cldlay(:,:,2) = cldlay2
      cldlay(:,:,3) = cldlay3
      nsublay(:,:,1) = nsublay1
      nsublay(:,:,2) = nsublay2
      nsublay(:,:,3) = nsublay3
#else
      cldlay = 0.0
      nsublay = 0.0
      do k = Nlevels, 1, -1
       do ic = 1, Ncolumns
        do ip = 1, Npoints

          ! Computation of the cloud fraction as a function of the temperature
          ! instead of height, for ice,liquid and all clouds
          IF (srok(ip,ic,k).gt.0.) THEN
            DO itemp=1,Ntemp
              IF ((tmp(ip,k).ge.tempmod(itemp)).AND.                           &
                (tmp(ip,k).lt.tempmod(itemp+1))) THEN
                lidarcldtempind(ip,itemp)=lidarcldtempind(ip,itemp)+1.
              END IF
            END DO
          END IF

          IF (cldy(ip,ic,k).eq.1.) THEN
            DO itemp=1,Ntemp
              IF((tmp(ip,k).ge.tempmod(itemp)).AND.                            &
                 (tmp(ip,k).lt.tempmod(itemp+1))) THEN
                lidarcldtemp(ip,itemp,1)=lidarcldtemp(ip,itemp,1)+1.
              END IF
            END DO
          END IF
!

          iz=1
          p1 = pplay(ip,k)
          if ( p1.gt.0. .and. p1.lt.(440.*100.)) then ! high clouds
            iz=3
          else if(p1.ge.(440.*100.) .and. p1.lt.(680.*100.)) then  ! mid clouds
            iz=2
         endif

         cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
         cldlay(ip,ic,4) = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
         lidarcld(ip,k)=lidarcld(ip,k) + cldy(ip,ic,k)

         nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
         nsublay(ip,ic,4) = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
         nsub(ip,k)=nsub(ip,k) + srok(ip,ic,k)

        enddo
       enddo
      enddo
#endif

! -- grid-box 3D cloud fraction

      where ( nsub(:,:).gt.0.0 )
         lidarcld(:,:) = lidarcld(:,:)/nsub(:,:)
      elsewhere
         lidarcld(:,:) = undef
      endwhere

! -- layered cloud fractions

      cldlayer = 0.0
      nsublayer = 0.0

      do iz = 1, Ncat
       do ic = 1, Ncolumns

          cldlayer(:,iz)=cldlayer(:,iz) + cldlay(:,ic,iz)
          nsublayer(:,iz)=nsublayer(:,iz) + nsublay(:,ic,iz)

       enddo
      enddo
      where ( nsublayer(:,:).gt.0.0 )
         cldlayer(:,:) = cldlayer(:,:)/nsublayer(:,:)
      elsewhere
         cldlayer(:,:) = undef
      endwhere

  ! ---------------------------------------------------------------
  ! 4- grid-box 3D cloud Phase :
  ! ---------------------------------------------------------------
  ! ---------------------------------------------------------------
  ! 4.1 - For Cloudy pixels with 8.16km < z < 19.2km
  ! ---------------------------------------------------------------
  DO ncol=1,Ncolumns
     DO i=1,Npoints

        DO nlev=Nlevels,18,-1  ! from 19.2km until 8.16km
           p1 = pplay(i,nlev)

           ! Avoid zero values
           IF ((cldy(i,ncol,nlev).eq.1.).AND.(ATBperp(i,ncol,nlev).gt.0.)) THEN
              ! Computation of the ATBperp along the phase discrimination line
              ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 +                    &
                   (ATB(i,ncol,nlev)**4)*beta50 +                     &
                   (ATB(i,ncol,nlev)**3)*gamma50 +                    &
                   (ATB(i,ncol,nlev)**2)*delta50 +                    &
                   ATB(i,ncol,nlev)*epsilon50 + zeta50

              ! ---------------------------------------------------------------
              !4.1.a Ice: ATBperp above the phase discrimination line
              ! ---------------------------------------------------------------
              IF( (ATBperp(i,ncol,nlev)-ATBperp_tmp).ge.0. )THEN   ! Ice clouds
                 ! ICE with temperature above 273.15K = Liquid (false ice)
                 IF(tmp(i,nlev).gt.273.15)THEN ! Temperature above 273,15 K
                    ! Liquid: False ice corrected by the temperature to Liquid
                    ! false ice detection ==> added to Liquid
                    lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1.
                    tmpl(i,ncol,nlev)=tmp(i,nlev)
                    ! keep the information "temperature criterium used"
                    lidarcldphase(i,nlev,5)=lidarcldphase(i,nlev,5)+1.
                    ! to classify the phase cloud
                    cldlayphase(i,ncol,4,2) = 1.  ! tot cloud
                    ! Cloud layers
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,2) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,2) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,2) = 1.
                    END IF
                    cldlayphase(i,ncol,4,5) = 1. ! tot cloud
                    ! Cloud layers
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,5) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,5) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,5) = 1.
                    END IF

                 ELSE
                    ! ICE with temperature below 273.15K
                    lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1.
                    tmpi(i,ncol,nlev)=tmp(i,nlev)
                    cldlayphase(i,ncol,4,1) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,1) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,1) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,1) = 1.
                    END IF

                 END IF

                 !__________________________________________________________
                 !
                 ! 4.1.b Liquid: ATBperp below the phase discrimination line
                 !__________________________________________________________
                 !
              ELSE                                        ! Liquid clouds
                 ! Liquid with temperature above 231.15K
                 IF(tmp(i,nlev).gt.231.15)THEN 
                    lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1.
                    tmpl(i,ncol,nlev)=tmp(i,nlev)
                    cldlayphase(i,ncol,4,2) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,2) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,2) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,2) = 1.
                    END IF

                 ELSE
                    ! Liquid with temperature below 231.15K = Ice (false liquid)
                    tmpi(i,ncol,nlev)=tmp(i,nlev)
                    ! false liquid detection ==> added to ice
                    lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1.
                    ! keep the information "temperature criterium used"
                    lidarcldphase(i,nlev,4)=lidarcldphase(i,nlev,4)+1.
                    ! to classify the phase cloud
                    cldlayphase(i,ncol,4,4) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,4) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,4) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,4) = 1.
                    END IF
                    cldlayphase(i,ncol,4,1) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,1) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,1) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,1) = 1.
                    END IF

                 END IF

              END IF  ! end of discrimination condition 
           END IF  ! end of cloud condition
        END DO ! end of altitude loop



        ! ---------------------------------------------------------------
        ! 4.2 - For Cloudy pixels with 0km < z < 8.16km
        ! ---------------------------------------------------------------

        toplvlsat=0
        DO nlev=17,1,-1  ! from 8.16km until 0km
           p1 = pplay(i,nlev)

           IF( (cldy(i,ncol,nlev).eq.1.).AND.(ATBperp(i,ncol,nlev).gt.0.) )THEN
              ! Phase discrimination line : ATBperp = ATB^5*alpha50 + 
              !   ATB^4*beta50 + ATB^3*gamma50 + ATB^2*delta50
              !   + ATB*epsilon50 + zeta50
              ! Computation of the ATBperp of the phase discrimination line
              ATBperp_tmp = (ATB(i,ncol,nlev)**5)*alpha50 +                    &
                            (ATB(i,ncol,nlev)**4)*beta50 +                     &
                            (ATB(i,ncol,nlev)**3)*gamma50 +                    &
                            (ATB(i,ncol,nlev)**2)*delta50 +                    &
                             ATB(i,ncol,nlev)*epsilon50 + zeta50
              !________________________________________________________
              !
              ! 4.2.a Ice: ATBperp above the phase discrimination line
              !________________________________________________________
              !
              ! ICE with temperature above 273.15K = Liquid (false ice)
              IF( (ATBperp(i,ncol,nlev)-ATBperp_tmp).ge.0. )THEN   ! Ice clouds
                 IF(tmp(i,nlev).gt.273.15)THEN
                    ! false ice ==> liq
                    lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1.
                    tmpl(i,ncol,nlev)=tmp(i,nlev)
                    lidarcldphase(i,nlev,5)=lidarcldphase(i,nlev,5)+1.

                    cldlayphase(i,ncol,4,2) = 1. ! tot cloud
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,2) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,2) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,2) = 1.
                    END IF

                    cldlayphase(i,ncol,4,5) = 1. ! tot cloud
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,5) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,5) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,5) = 1.
                    END IF

                 ELSE
                    ! ICE with temperature below 273.15K
                    lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1.
                    tmpi(i,ncol,nlev)=tmp(i,nlev)

                    cldlayphase(i,ncol,4,1) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,1) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,1) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,1) = 1.
                    END IF

                 END IF

                 !__________________________________________________________
                 !
                 ! 4.2.b Liquid: ATBperp below the phase discrimination line
                 !__________________________________________________________
                 !
              ELSE  
                 ! Liquid with temperature above 231.15K
                 IF(tmp(i,nlev).gt.231.15)THEN 
                    lidarcldphase(i,nlev,2)=lidarcldphase(i,nlev,2)+1.
                    tmpl(i,ncol,nlev)=tmp(i,nlev)

                    cldlayphase(i,ncol,4,2) = 1.
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,2) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,2) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,2) = 1.
                    END IF

                 ELSE
                    ! Liquid with temperature below 231.15K = Ice (false liquid)
                    tmpi(i,ncol,nlev)=tmp(i,nlev)
                    ! false liq ==> ice
                    lidarcldphase(i,nlev,1)=lidarcldphase(i,nlev,1)+1.
                    ! false liq ==> ice
                    lidarcldphase(i,nlev,4)=lidarcldphase(i,nlev,4)+1.

                    cldlayphase(i,ncol,4,4) = 1. ! tot cloud
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,4) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,4) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,4) = 1.
                    END IF

                    cldlayphase(i,ncol,4,1) = 1. ! tot cloud
                    IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                       cldlayphase(i,ncol,3,1) = 1.
                    ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                       cldlayphase(i,ncol,2,1) = 1.
                    ELSE
                       cldlayphase(i,ncol,1,1) = 1.
                    END IF

                 END IF
              END IF  ! end of discrimination condition 

              toplvlsat=0

              ! Find the level of the highest cloud with SR>30
              IF(x(i,ncol,nlev).gt.S_cld_att)THEN  ! SR > 30.
                 toplvlsat=nlev-1
                 GOTO 9999
              END IF

           END IF  ! end of cloud condition
        END DO  ! end of altitude loop

9999    CONTINUE

        !___________________________________________________________________
        !
        ! Undefined phase: For a cloud located below another cloud with SR>30 
        ! see Cesana and Chepfer 2013 Sect.III.2
        !___________________________________________________________________
        !
        IF(toplvlsat.ne.0)THEN
           DO nlev=toplvlsat,1,-1
              p1 = pplay(i,nlev)
              IF(cldy(i,ncol,nlev).eq.1.)THEN
                 lidarcldphase(i,nlev,3)=lidarcldphase(i,nlev,3)+1.
                 tmpu(i,ncol,nlev)=tmp(i,nlev)

                 cldlayphase(i,ncol,4,3) = 1. ! tot cloud
                 IF ( p1.gt.0. .AND. p1.lt.(440.*100.)) THEN
                    cldlayphase(i,ncol,3,3) = 1.
                 ELSE IF(p1.ge.(440.*100.) .AND. p1.lt.(680.*100.)) THEN
                    cldlayphase(i,ncol,2,3) = 1.
                 ELSE
                    cldlayphase(i,ncol,1,3) = 1.
                 END IF

              END IF
           END DO
        END IF

        toplvlsat=0

     END DO
  END DO



  !________________________________________________________
  !
  ! Computation of final cloud phase diagnosis
  !________________________________________________________
  !

  ! Compute the Ice percentage in cloud = ice/(ice+liq) as a function
  ! of the occurrences
  lidarcldphasetmp(:,:)=lidarcldphase(:,:,1)+lidarcldphase(:,:,2);
  WHERE (lidarcldphasetmp(:,:).gt. 0.)
     lidarcldphase(:,:,6)=lidarcldphase(:,:,1)/lidarcldphasetmp(:,:)
  ELSEWHERE
     lidarcldphase(:,:,6) = undef
  END WHERE

  ! Compute Phase 3D Cloud Fraction
  WHERE ( nsub(:,:).gt.0.0 )
     lidarcldphase(:,:,1)=lidarcldphase(:,:,1)/nsub(:,:)
     lidarcldphase(:,:,2)=lidarcldphase(:,:,2)/nsub(:,:)
     lidarcldphase(:,:,3)=lidarcldphase(:,:,3)/nsub(:,:)
     lidarcldphase(:,:,4)=lidarcldphase(:,:,4)/nsub(:,:)
     lidarcldphase(:,:,5)=lidarcldphase(:,:,5)/nsub(:,:)
  ELSEWHERE
     lidarcldphase(:,:,1) = undef
     lidarcldphase(:,:,2) = undef
     lidarcldphase(:,:,3) = undef
     lidarcldphase(:,:,4) = undef
     lidarcldphase(:,:,5) = undef
  END WHERE


  ! Compute Phase low mid high cloud fractions
  DO iz = 1, Ncat
     DO i=1,Nphase-3
        DO ic = 1, Ncolumns
           cldlayerphase(:,iz,i)=cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)
           cldlayerphasesum(:,iz)=cldlayerphasesum(:,iz)+cldlayphase(:,ic,iz,i)
        END DO
     END DO
  END DO

  DO iz = 1, Ncat
     DO i=4,5
        DO ic = 1, Ncolumns
           cldlayerphase(:,iz,i)=cldlayerphase(:,iz,i) + cldlayphase(:,ic,iz,i)
        END DO
     END DO
  END DO

  ! Compute the Ice percentage in cloud = ice/(ice+liq)
  cldlayerphasetmp(:,:)=cldlayerphase(:,:,1)+cldlayerphase(:,:,2)
  WHERE (cldlayerphasetmp(:,:).gt. 0.)
     cldlayerphase(:,:,6)=cldlayerphase(:,:,1)/cldlayerphasetmp(:,:)
  ELSEWHERE
     cldlayerphase(:,:,6) = undef
  END WHERE

  DO i=1,Nphase-1
     WHERE ( cldlayerphasesum(:,:).gt.0.0 )
        cldlayerphase(:,:,i) = (cldlayerphase(:,:,i)/cldlayerphasesum(:,:)) *  &
                                cldlayer(:,:) 
     END WHERE
  END DO


  DO i=1,Npoints
     DO iz=1,Ncat
        checkcldlayerphase=0.
        checkcldlayerphase2=0.

        IF (cldlayerphasesum(i,iz).gt.0.0 )THEN
           DO ic=1,Nphase-3
              checkcldlayerphase=checkcldlayerphase+cldlayerphase(i,iz,ic)
           END DO
           checkcldlayerphase2=cldlayer(i,iz)-checkcldlayerphase

        END IF

     END DO
  END DO

  DO i=1,Nphase-1
     WHERE ( nsublayer(:,:).eq.0.0 )
        cldlayerphase(:,:,i) = undef
     END WHERE
  END DO



  ! Compute Phase 3D as a function of temperature
  DO nlev=1,Nlevels
     DO ncol=1,Ncolumns
        DO i=1,Npoints
           DO itemp=1,Ntemp
              IF(tmpi(i,ncol,nlev).gt.0.)THEN
                 IF( (tmpi(i,ncol,nlev).ge.tempmod(itemp)).AND.                &
                     (tmpi(i,ncol,nlev).lt.tempmod(itemp+1)) )THEN
                    lidarcldtemp(i,itemp,2)=lidarcldtemp(i,itemp,2)+1.
                 END IF
              ELSE IF(tmpl(i,ncol,nlev).gt.0.)THEN
                 IF( (tmpl(i,ncol,nlev).ge.tempmod(itemp)).AND.                &
                     (tmpl(i,ncol,nlev).lt.tempmod(itemp+1)) )THEN
                    lidarcldtemp(i,itemp,3)=lidarcldtemp(i,itemp,3)+1.
                 END IF
              ELSE IF(tmpu(i,ncol,nlev).gt.0.)THEN
                 IF( (tmpu(i,ncol,nlev).ge.tempmod(itemp)).AND.                &
                     (tmpu(i,ncol,nlev).lt.tempmod(itemp+1)) )THEN
                    lidarcldtemp(i,itemp,4)=lidarcldtemp(i,itemp,4)+1.
                 END IF
              END IF
           END DO
        END DO
     END DO
  END DO

  ! Check temperature cloud fraction
  DO i=1,Npoints
     DO itemp=1,Ntemp
        checktemp=lidarcldtemp(i,itemp,2)+lidarcldtemp(i,itemp,3)+             &
                  lidarcldtemp(i,itemp,4)
     END DO
  END DO

  ! Compute the Ice percentage in cloud = ice/(ice+liq)
  !   sumlidarcldtemp=sum(lidarcldtemp(:,:,2:3),3)
  sumlidarcldtemp(:,:)=lidarcldtemp(:,:,2)+lidarcldtemp(:,:,3)

  WHERE(sumlidarcldtemp(:,:)>0.)
     lidarcldtemp(:,:,5)=lidarcldtemp(:,:,2)/sumlidarcldtemp(:,:)
  ELSEWHERE
     lidarcldtemp(:,:,5)=undef
  END WHERE

  DO i=1,4
     WHERE(lidarcldtempind(:,:).gt.0.)
        lidarcldtemp(:,:,i) = lidarcldtemp(:,:,i)/lidarcldtempind(:,:)
     ELSEWHERE
        lidarcldtemp(:,:,i) = undef
     END WHERE
  END DO

  RETURN
  END SUBROUTINE COSP_CLDFRAC
! ---------------------------------------------------------------

END MODULE MOD_LMD_IPSL_STATS
