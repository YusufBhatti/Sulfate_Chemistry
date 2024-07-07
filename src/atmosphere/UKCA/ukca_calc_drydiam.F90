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
!   Calculates geometric mean dry diameter for multi-component
!   aerosol population which is lognormally distributed with
!   number concentration ND in mode, component mass concentration
!   MD in mode, component density RHOCOMP and component molecular mass MM
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
MODULE ukca_calc_drydiam_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CALC_DRYDIAM_MOD'

CONTAINS

SUBROUTINE ukca_calc_drydiam(nbox,nd,md,mdt,drydp,dvol)

!-------------------------------------------------------------
!
!     Calculates geometric mean dry diameter for multi-component
!     aerosol population which is lognormally distributed with
!     number concentration ND in mode,
!     cpt mass concentration MD in mode,
!     cpt density RHOCOMP and cpt molecular mass MM
!
!     Calculate dry volume per particle using composition info as:
!
!     DVOL = sum_cp { MD(ICP)*MM(ICP)/(AVC*RHOCOMP(ICP)) }
!
!     Where AVC is Avogadro's constant. Then, from Jacobsen,
!     "Fundamentals of Atmospheric Modeling", pp. 412, have
!
!     dry volume conc. = ND*(PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     i.e. DVOL  = (PI/6)*(Dp^3)*exp{9/2 log^2(sigma_g)}
!
!     where Dp is the number mean dry diameter,
!     and sigma_g is the geometric standard deviation.
!
!     Then calculate Dp as:
!
!     Dp=CBRT( DVOL*(6/PI)/EXP{9.2 log^2(sigma_g)} )
!
!     Inputs
!     ------
!     NBOX      : Number of grid boxes
!     ND        : Aerosol ptcl no. concentration (ptcls per cc)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!
!     Outputs
!     -------
!     DRYDP     : Median particle dry diameter for each mode (m)
!     DVOL      : Median particle dry volume for each mode (m^3)
!     MD        : Component median aerosol mass (molecules per ptcl)
!     MDT       : Total median aerosol mass (molecules per ptcl)
!
!     Local Variables
!     ---------------
!     None
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     AVC       : Avogadro's constant (per mole)
!     MMSUL     : Molar mass of a pure H2SO4 aerosol (kg per mole)
!     RHOSUL    : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES    : Number of aerosol modes
!     NCP       : Number of aerosol components
!     MM        : Molecular mass of each component
!     RHOCOMP   : Densities (dry) of each component (kg/m^3)
!     MODE      : Which modes are being carried
!     COMPONENT : Which components are in each of modes
!     DDPLIM0   : Lower limits for dry diameter for each mode (m)
!     MFRAC_0   : Initial mass fraction to set when no particles.
!     X         : EXP((9/2)*LOG^2(SIGMA_G))
!     NUM_EPS   : Value of NEWN below which do not carry out process
!     MMID      : Mid-point masses for initial radius grid
!     MLO       : Lo-interf masses for initial radius grid
!
!--------------------------------------------------------------------
USE ukca_constants
USE conversions_mod, ONLY: pi
USE ukca_mode_setup, ONLY: nmodes, ncp, mm, rhocomp, mode,          &
                           component, ddplim0, mfrac_0, x, num_eps, &
                           mmid, mlo
USE vectlib_mod,     ONLY: cubrt_v
USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
USE ereport_mod,     ONLY: ereport
USE umPrintMgr,      ONLY: umPrint, umMessage
USE um_types,        ONLY: log_small


USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Interface
INTEGER, INTENT(IN) :: nbox
REAL, INTENT(IN)    :: nd(nbox,nmodes)
REAL, INTENT(INOUT) :: md(nbox,nmodes,ncp)
REAL, INTENT(INOUT) :: mdt(nbox,nmodes)
REAL, INTENT(OUT)   :: drydp(nbox,nmodes)
REAL, INTENT(OUT)   :: dvol(nbox,nmodes)

! Local variables
INTEGER :: jl
INTEGER :: imode
INTEGER :: icp
REAL    :: ddpcub(nbox,nmodes)
REAL    :: sixovrpix(nmodes)
LOGICAL (KIND=log_small) :: mask(nbox)
REAL    :: ratio1(ncp)
REAL    :: ratio2(nmodes)
REAL    :: dp_thresh1
REAL    :: dp
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER           :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_DRYDIAM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Below is over ncp
ratio1(:)=mm(:)/(avc*rhocomp(:))
!
! Below is over NMODES
sixovrpix(:)=6.0/(pi*x(:))
ratio2(:)=mmsul*mmid(:)/(avc*rhosul)
!
DO imode=1,nmodes
  IF (mode(imode)) THEN
    mask(:) =(nd(:,imode) > num_eps(imode))
    WHERE (mask(:))
      dvol(:,imode)=0.0
    ELSEWHERE
      dvol(:,imode)=mmid(imode)*mmsul/(avc*rhosul)
    END WHERE
    !     calculate particle dry volume using composition info
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        WHERE (mask(:))                                                &
         dvol(:,imode)=dvol(:,imode)+ratio1(icp)*md(:,imode,icp)
      END IF
    END DO
    !     DVOL calculates particle dry volume assuming pure H2SO4
  ELSE
    dvol(:,imode)=ratio2(imode)
  END IF
END DO

DO imode=1,nmodes
  ddpcub(:,imode)=sixovrpix(imode)*dvol(:,imode)
END DO
CALL cubrt_v(nmodes*nbox, ddpcub, drydp)

! .. also check whether mean diameter too low for mode
!      DO IMODE=1,NMODES
DO imode=1,3 ! only do check for modes 1, 2 and 3
  !      DO IMODE=1,1 ! only do check for mode 1
  IF (mode(imode)) THEN
    DO jl=1,nbox
      dp    =drydp(jl,imode)
      dp_thresh1=ddplim0(imode)*0.1
      IF (dp < dp_thresh1) THEN
        DO icp=1,ncp
          IF (component(imode,icp)) THEN
            md(jl,imode,icp)=mlo(imode)*mfrac_0(imode,icp)
          END IF
        END DO
        mdt(jl,imode)=mlo(imode)
        dvol(jl,imode)=mlo(imode)*mmsul/(avc*rhosul)
        drydp(jl,imode)=(sixovrpix(imode)*dvol(jl,imode))**(1.0/3.0)
      END IF ! if DP<DPTHRESH
    END DO
  END IF
END DO

DO imode=1,nmodes
  IF ((MINVAL(dvol (:,imode)) <= 0.0) .OR.                          &
     (MINVAL(drydp(:,imode)) <= 0.0)) THEN
    cmessage = ' dvol or drydp <= 0'
    WRITE(umMessage,*) 'In calcdrydiam: drydp (min,max,sum) imode=',imode
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) MINVAL(drydp(:,imode)),MAXVAL(drydp(:,imode)),       &
               SUM(drydp(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) 'Location of min: ',MINLOC(drydp(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) 'Location of max: ',MAXLOC(drydp(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) 'In calcdrydiam: dvol (min,max,sum) imode=',imode
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) MINVAL(dvol(:,imode)),MAXVAL(dvol(:,imode)),         &
               SUM(dvol(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) 'Location of min: ',MINLOC(drydp(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    WRITE(umMessage,*) 'Location of max: ',MAXLOC(drydp(:,imode))
    CALL umPrint(umMessage,src='ukca_calc_drydiam')
    errcode=1

    CALL ereport('ukca_calc_drydiam',errcode,cmessage)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_calc_drydiam
END MODULE ukca_calc_drydiam_mod
