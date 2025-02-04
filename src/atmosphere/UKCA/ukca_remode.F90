! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Carry out mode-merging algorithm where average mass for
!    mixed nucl,Aitken,accum modes exceeds mid-point mass of next
!    mode up.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
MODULE ukca_remode_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_REMODE_MOD'

CONTAINS

SUBROUTINE ukca_remode(nbox,nd,md,mdt,drydp,wetdp,verbose,        &
  imerge,bud_aer_mas,n_merge_1d,pmid)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Carry out mode-merging algorithm where average mass for
! mixed nucl,Aitken,accum modes exceeds mid-point mass of next
! mode up. Transfer fraction of mode number and mass which is greater
! than this threshold. Calculate fraction using error function to
! evaluate integrals of log-normal functions for number and mass
! to threshold. This then gives the amounts to transfer to next mode.
! Evaluate error function.
! Re-calculates ND,MD,MDT according to amounts to transfer.
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DRYDP    : Median dry diameter for particles in size mode (m)
! WETDP    : Median wet diameter for particles in size mode (m)
! VERBOSE  : Switch for printing out test print statements
! IMERGE   : Switch to use mid-pts (=1), edges (2) or dynamic (=3)
!
! Outputs
! -------
! ND       : Modified number concentration for each mode
! MD       : Modified average cpt particle mass for each mode
! MDT      : Modified average total particle mass for each mode
! BUD_AER_MAS : Aerosol mass budget terms
! N_MERGE_1D  : # of merges (grown out of bounds) in each box, each mode
!
! Local variables
! ---------------
! LNRATN   : Log of ratio of threshold diameter to no. median diameter
! LNRATM   : Log of ratio of threshold diameter to mass median diameter
! ERFNUM   : LNRATN/(sqrt(2)*log(sigmag))
! ERFMAS   : LNRATM/(sqrt(2)*log(sigmag))
! FRAC_N   : Fraction of ptcl number which is within bounds
! FRAC_M   : Fraction of ptcl mass which is within bounds
! DELN     : Number concentration to transfer due to mode-merging
! DM       : Cpt mass concentration to transfer due to mode-merging
! DP       : Number median diameter of mode (m)
! DP2      : Volume median diameter of mode (m)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Logical variable denoting where mode is defined
! COMPONENT: Logical variable denoting where cpt is defined
! DDPMID   : Mid-point of size mode = exp(0.5*(lndp0+lndp1)) (m)
! MMID     : Ptcl mass with dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
! MFRAC_0  : Initial mass fraction to set when no particles.
! SIGMAG   : Geometric standard deviation for each mode
! DDPLIM0  : Lower limit for dry diameter in mode (m)
! DDPLIM1  : Upper limit for dry diameter in mode (m)
! NUM_EPS  : Value of NEWN below which don't recalculate MD
!                                            or carry out process
! CP_SU    : Index of component where SO4    cpt is stored
! CP_BC    : Index of component where BC     cpt is stored
! CP_OC    : Index of component where 1st OC cpt is stored
! CP_CL    : Index of component where NaCl   cpt is stored
! CP_SO    : Index of component where 2nd OC cpt is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
USE ukca_mode_setup, ONLY: nmodes, ncp, mode, component, ddpmid,    &
                           mmid, mfrac_0, sigmag, ddplim0, ddplim1, &
                           num_eps, cp_su, cp_bc, cp_oc, cp_cl,     &
                           cp_so
USE ukca_setup_indices
USE yomhook,   ONLY: lhook, dr_hook
USE parkind1,  ONLY: jprb, jpim
USE umErf_mod, ONLY: umErf

USE um_types, ONLY: integer32
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: nbox
INTEGER, INTENT(IN) :: verbose
INTEGER, INTENT(IN) :: imerge
REAL, INTENT(IN)    :: drydp(nbox,nmodes)
REAL, INTENT(IN)    :: wetdp(nbox,nmodes)
REAL, INTENT(IN)    :: pmid(nbox)

INTEGER (KIND=integer32), INTENT(INOUT) :: n_merge_1d(nbox,nmodes)
REAL, INTENT(INOUT) :: nd(nbox,nmodes)
REAL, INTENT(INOUT) :: md(nbox,nmodes,ncp)
REAL, INTENT(INOUT) :: mdt(nbox,nmodes)
REAL, INTENT(INOUT) :: bud_aer_mas(nbox,0:nbudaer)

! Local variables
INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
INTEGER :: iimode
REAL    :: dp
REAL    :: dp_ip1
REAL    :: dp_thresh1
REAL    :: dp_thresh2
REAL    :: lnratn
REAL    :: erfnum
REAL    :: frac_n
REAL    :: deln
REAL    :: log2sg
REAL    :: lnratm
REAL    :: erfmas
REAL    :: frac_m
REAL    :: dm(ncp)
REAL    :: dmt
REAL    :: dp2
REAL    :: newn
REAL    :: newnp1
!
INTEGER (KIND=integer32) :: nmodemax_merge(nbox)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_REMODE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nmodemax_merge(:)=3
WHERE (pmid(:) < 1.0e4) nmodemax_merge(:)=2

DO jl=1,nbox
  DO imode=1,nmodemax_merge(jl)
    ! do for soluble nucl, Aitken, accum in troposphere
    ! do for soluble nucl, Aitken  only  in stratosphere
    IF (mode(imode)) THEN

      dp=drydp(jl,imode)
      dp_ip1=drydp(jl,imode+1)

      ! IMERGE=1 --> Apply as originally done in M7 (next mode mid-pt)
      IF (imerge == 1) THEN
        dp_thresh1=ddpmid(imode+1)
        dp_thresh2=ddpmid(imode+1)
      END IF

      ! IMERGE=2 --> Apply whenever DP is out of bounds
      IF (imerge == 2) THEN
        dp_thresh1=ddplim0(imode+1)
        dp_thresh2=ddplim0(imode+1)
      END IF

      ! IMERGE=3 --> Keep transferring up top-tail of mode (dynamic)
      IF (imerge == 3) THEN
        dp_thresh1=SQRT(dp*dp_ip1)
        dp_thresh2=SQRT(dp*dp_ip1)
      END IF

      ! .. 1st threshold determines if mode-merging to take place.
      IF ((dp > dp_thresh1) .OR. (imerge == 3)) THEN

        IF (nd(jl,imode) > num_eps(imode)) THEN

          n_merge_1d(jl,imode)=n_merge_1d(jl,imode)+1

          ! .. 2nd threshold determines fraction of number/mass
          lnratn=LOG(dp_thresh2/dp)
          erfnum=lnratn/SQRT(2.0)/LOG(sigmag(imode))

          ! .. FRAC_N is fraction of number remaining in mode
          frac_n=0.5*(1.0+umErf(erfnum))

          ! .. limit DELN to be max = half # of ptcls
          IF (frac_n < 0.5) frac_n=0.5

          deln=nd(jl,imode)*(1.0-frac_n)

          log2sg=LOG(sigmag(imode))*LOG(sigmag(imode))
          dp2=EXP(LOG(dp)+3.0*log2sg) ! volume median diameter

          ! .. 2nd threshold determines fraction of number/mass
          lnratm=LOG(dp_thresh2/dp2)
          erfmas=lnratm/SQRT(2.0)/LOG(sigmag(imode))
          frac_m=0.5*(1.0+umErf(erfmas))

          ! .. limit DELM to be max 99.9%
          IF (frac_m < 0.001) frac_m=0.001

          ! .. calculate new number concs for mode and larger mode
          newn=nd(jl,imode)-deln
          newnp1=nd(jl,imode+1)+deln

          !-----------------------------------------------------------------------
          !
          ! .. This section updates mode from which merging is taking place

          IF (newn > num_eps(imode)) THEN

            DO icp=1,ncp
              IF (component(imode,icp)) THEN

                ! .. calculate cpt mass conc to transfer to next mode (use old no/mass)
                dm(icp)=md(jl,imode,icp)*nd(jl,imode)*(1.0-frac_m)

                IF (imode == 1) THEN
                  IF ((icp == cp_su) .AND. (nmasmergsuintr12 > 0))            &
                         bud_aer_mas(jl,nmasmergsuintr12)=                 &
                         bud_aer_mas(jl,nmasmergsuintr12)+dm(icp)
                  IF ((icp == cp_oc) .AND. (nmasmergocintr12 > 0))            &
                         bud_aer_mas(jl,nmasmergocintr12)=                 &
                         bud_aer_mas(jl,nmasmergocintr12)+dm(icp)
                  IF ((icp == cp_so) .AND. (nmasmergsointr12 > 0))            &
                         bud_aer_mas(jl,nmasmergsointr12)=                 &
                         bud_aer_mas(jl,nmasmergsointr12)+dm(icp)
                END IF

                IF (imode == 2) THEN
                  IF ((icp == cp_su) .AND. (nmasmergsuintr23 > 0))            &
                         bud_aer_mas(jl,nmasmergsuintr23)=                 &
                         bud_aer_mas(jl,nmasmergsuintr23)+dm(icp)
                  IF ((icp == cp_bc) .AND. (nmasmergbcintr23 > 0))            &
                         bud_aer_mas(jl,nmasmergbcintr23)=                 &
                         bud_aer_mas(jl,nmasmergbcintr23)+dm(icp)
                  IF ((icp == cp_oc) .AND. (nmasmergocintr23 > 0))            &
                         bud_aer_mas(jl,nmasmergocintr23)=                 &
                         bud_aer_mas(jl,nmasmergocintr23)+dm(icp)
                  IF ((icp == cp_so) .AND. (nmasmergsointr23 > 0))            &
                         bud_aer_mas(jl,nmasmergsointr23)=                 &
                         bud_aer_mas(jl,nmasmergsointr23)+dm(icp)
                END IF

                IF (imode == 3) THEN
                  IF ((icp == cp_su) .AND. (nmasmergsuintr34 > 0))            &
                         bud_aer_mas(jl,nmasmergsuintr34)=                 &
                         bud_aer_mas(jl,nmasmergsuintr34)+dm(icp)
                  IF ((icp == cp_bc) .AND. (nmasmergbcintr34 > 0))            &
                         bud_aer_mas(jl,nmasmergbcintr34)=                 &
                         bud_aer_mas(jl,nmasmergbcintr34)+dm(icp)
                  IF ((icp == cp_oc) .AND. (nmasmergocintr34 > 0))            &
                         bud_aer_mas(jl,nmasmergocintr34)=                 &
                         bud_aer_mas(jl,nmasmergocintr34)+dm(icp)
                  IF ((icp == cp_cl) .AND. (nmasmergssintr34 > 0))            &
                         bud_aer_mas(jl,nmasmergssintr34)=                 &
                         bud_aer_mas(jl,nmasmergssintr34)+dm(icp)
                  IF ((icp == cp_so) .AND. (nmasmergsointr34 > 0))            &
                         bud_aer_mas(jl,nmasmergsointr34)=                 &
                         bud_aer_mas(jl,nmasmergsointr34)+dm(icp)
                END IF

              ELSE
                dm(icp)=0.0
              END IF
            END DO
            !
            ! .. first remove mass to be transferred from mode IMODE
            mdt(jl,imode)=0.0
            DO icp=1,ncp
              IF (component(imode,icp)) THEN
                md(jl,imode,icp)=                                         &
                    (nd(jl,imode)*md(jl,imode,icp)-dm(icp))/newn
                mdt(jl,imode)=mdt(jl,imode)+md(jl,imode,icp)
              ELSE
                md(jl,imode,icp)=0.0
              END IF ! COMPONENT(IMODE,ICP)
            END DO

            ! .. now set new number to mode IMODE
            nd(jl,imode)=newn ! set particle number to new value
            !
            !-----------------------------------------------------------------------
            !
            ! .. This section updates mode IMODE+1

            ! .. Update MD for this mode
            mdt(jl,imode+1)=0.0
            DO icp=1,ncp
              IF (component(imode+1,icp)) THEN
                md(jl,imode+1,icp)=                                       &
                    (nd(jl,imode+1)*md(jl,imode+1,icp)+dm(icp))/newnp1
                mdt(jl,imode+1)=mdt(jl,imode+1)+md(jl,imode+1,icp)
              ELSE
                md(jl,imode+1,icp)=0.0
              END IF ! COMPONENT(IMODE+1,ICP)
            END DO

            ! .. now set new number to mode IMODE+1
            nd(jl,imode+1)=newnp1

          END IF ! IF NEWNP1>0
          !
          !----------------------------------------------------------------------

        ELSE

          DO icp=1,ncp
            IF (component(imode,icp)) THEN
              md(jl,imode,icp)=mmid(imode)*mfrac_0(imode,icp)
            END IF
          END DO
          mdt(jl,imode)=mmid(imode)

        END IF ! if significant number of particles in lower mode

      END IF ! if DP outside limits (mode-merge criterion)

    END IF ! IF MODE(IMODE)
  END DO ! IMODE=1,3 (only merge for soluble nuc/Ait/acc)
END DO ! JL=1,NBOX

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_remode
END MODULE ukca_remode_mod
