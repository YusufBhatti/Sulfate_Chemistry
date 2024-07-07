! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme Core.
! Subroutine Interface:

MODULE gw_ussp_core_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GW_USSP_CORE_MOD'

CONTAINS

SUBROUTINE gw_ussp_core(rows, row_length, launchlev,                    &
  i_start, i_end, j_start, j_end, k_start, k_end, ussp_launch_factor,   &
  wavelstar, cgw_scale_factor, two_omega, sin_theta_latitude, rho_th,   &
  nbv, udotk, totalppn, fptot, L_add_cgw)
!
! purpose: This subroutine calculates the vertically propagated flux of
!          pseudomomentum due to gravity waves as parametrised by the
!          Ultra Simple Spectral gravity wave Parametrization
!         (originally Warner and McIntyre, since reengineered for use in
!          the UM and extended).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Gravity Wave Drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

USE conversions_mod, ONLY: pi
USE gw_ussp_params_mod, ONLY: idir
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! ----------------------------------------------------------------------+-------
! Subroutines defined from modules
! ----------------------------------------------------------------------+-------
      
! ----------------------------------------------------------------------+-------

IMPLICIT NONE

! ----------------------------------------------------------------------+-------
! Subroutine arguments of gw_ussp_core
! ----------------------------------------------------------------------+-------
! Number of rows for u field
INTEGER, INTENT(IN) :: rows
!
! Number of grid points in row
INTEGER, INTENT(IN) :: row_length
!
! Level for gw launch - in principle this could be launchlev(i,j)
INTEGER, INTENT(IN) :: launchlev
!
! Starting grid point in local row
INTEGER, INTENT(IN) :: i_start
!
! Final grid point in local row
INTEGER, INTENT(IN) :: i_end
!
! Starting row of local grid
INTEGER, INTENT(IN) :: j_start
!
! Final row of local grid 
INTEGER, INTENT(IN) :: j_end
!
! Starting level
INTEGER, INTENT(IN) :: k_start
!
! Top model level
INTEGER, INTENT(IN) :: k_end
!
! Factor enhancement for invariant global wave launch amplitude
REAL, INTENT(IN)    :: ussp_launch_factor 
!
! omega is Angular speed of Earth's rotation, 2 * omega
REAL, INTENT(IN)    :: two_omega
!
! Scale conversion factor from convective rain to GW flux
REAL, INTENT(IN)    :: cgw_scale_factor
!
! Characteristic (spectrum peak) wavelength (m)
REAL, INTENT(IN)    :: wavelstar
!
! Latitudes on grid
REAL, INTENT(IN)    :: sin_theta_latitude(i_start:i_end,j_start:j_end)
!
! Rho (kg m^-3)
REAL, INTENT(IN)    :: rho_th(i_start:i_end,j_start:j_end,k_start:k_end)
!
! Buoyancy [Brunt Vaisala] Frequency (rad s^-1)
REAL, INTENT(IN)    :: nbv(i_start:i_end,j_start:j_end,k_start:k_end)
!
! Component of wind in phi_jdir direction (m s^-1)
REAL, INTENT(IN)    :: udotk(i_start:i_end,j_start:j_end,k_start:k_end, &
                                                                  idir)
!
! Total Precipitation for CGW source
REAL, INTENT(IN)    :: totalppn(i_start:i_end,j_start:j_end)
!
! Pseudomomentum flux integrated over azimuthal sector (kg m^-1 s^-2)
REAL, INTENT(OUT)   :: fptot(i_start:i_end,j_start:j_end,k_start:k_end, &
                                                                   idir)
!
! Flag to indicate CGW flux calculation should be carried out
LOGICAL, INTENT(IN) :: L_add_cgw
!
! ----------------------------------------------------------------------+-------
! Local parameters
! ----------------------------------------------------------------------+-------
!
! Maximum number of iterations of Newton Raphson DO (While) loop
INTEGER, PARAMETER :: maxwhile       = 9
!
! Parameter beta in the launch spectrum total energy equation
REAL, PARAMETER :: beta_e0           = 1.0227987125e-1
!
! Parameter p in B_0(p) for launch spectrum intrinsic frequency
! NOTE: This parameter determines the intrinsic frequency spectrum
!       shape and hence the integral form in 4.1, which is strictly
!       valid only for p > 1. !!IF contemplating changes BE WARNED!!
REAL, PARAMETER :: psat              = 5.0 / 3.0
!
! Psat - 1
REAL, PARAMETER :: psatm1            = psat - 1.0
!
! 2 - Psat
REAL, PARAMETER :: twompsat          = 2.0 - psat
!
! Reciprocal of max wavelength at launch (/m)
REAL, PARAMETER :: lminl             = 1.0/20000
!
! Minimum vertical wavenumber at launch (/m )
REAL, PARAMETER ::  mminl            = 2.0 * pi * lminl
!
! Power s of vertical wavenumber spectrum A_0(s,t) at low m
REAL, PARAMETER :: ss                = 1.0
!
! s + 1
REAL, PARAMETER :: ssp1              = ss + 1.0
!
! Power t=t_sat of vertical wavenumber spectrum at large m due to
! saturation by gravity wave breaking (and shape of chopping fn)
REAL, PARAMETER :: tt                = 3.0
!
! t - 1, t - 2, 1 / (t-2), (t-3) / (t-2), 2 - t
REAL, PARAMETER :: ttm1              = tt - 1.0
REAL, PARAMETER :: ttm2              = tt - 2.0
REAL, PARAMETER :: rttm2             = 1.0 / ttm2
REAL, PARAMETER :: ttrat             = (tt - 3.0) * rttm2
REAL, PARAMETER :: twomtt            = 2.0 - tt
!
! s + t, 1 / (s+t)
REAL, PARAMETER :: ssptt             = ss + tt
REAL, PARAMETER :: rssptt            = 1.0 / ssptt
!
! Weight for (n+1)th guess at mNlX in iteration solution
REAL, PARAMETER :: mweight           = 0.8
!
! Strength coefficient constant for Launch spectrum (CCL / A0)
REAL, PARAMETER :: ccl0              = 3.41910625e-9
!
! Scale conversion factor from convective rain to GW flux (mPa)
REAL, PARAMETER :: cor2flux          = 7.20507e-04
!
! ----------------------------------------------------------------------+-------
! Security parameters
! ----------------------------------------------------------------------+-------
!
! Minimum allowed non-zero value of precipitation (0.1 mm / day)
REAL, PARAMETER ::  rppnmin           = 1.0 / 1.1574e-06
!
! Minimum allowed non-zero value of Curvature Coefficient A
REAL, PARAMETER ::  asecp            =  1.0e-20
REAL, PARAMETER ::  asecn            = -(asecp)
!
! ----------------------------------------------------------------------+-------
! Local Constants (a) Physical
! ----------------------------------------------------------------------+-------
!
! Normalised minimum vertical wavenumber at launch
REAL            ::  mnlmin
!
! Wavenumber at peak in spectrum (/m)
REAL            ::  mstar
!
! Reciprocal of mstar (m) and mstar^2 (m^2)
REAL            ::  rmstar
REAL            ::  rmstarsq
!
! Equatorial planetary vorticity gradient parameter B_eq (/m /s )
REAL            :: beta_eq_rmstar
!
! Normalisation factor for the launch spectrum vertical wavenumber A0/(s+1)(t-1)
REAL            :: a0_r_sp1tm1
!
! Normalisation factor for the launch spectrum vertical wavenumber -A0/(t-1)
REAL            :: a0_r_1mt
!
! Integral segment [(s+1) * mNLmin**(1-t)]
REAL            :: tail_chop2b
!
! Integral segment [(t-1) * mNLmin**(s+1)]
REAL            :: head_chop2a
!
! Isotropic GW Launch Flux = mStar^(-2) * [CCL0 * ussp_launch_factor]
REAL            :: glob_launch_flux
!
! Scale conversion factor from convective rain to GW flux
! [cor2flux * cgw_scale_factor]
REAL            :: cgw_convert_factor
!
! ----------------------------------------------------------------------+-------
! Local variables (scalars) used in GW_USSP
! Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------
!
! Longitude index
INTEGER         :: i
!
! Latitude index
INTEGER         :: j
!
! Level index
INTEGER         :: k
!
! Index values for chop case loops
INTEGER         :: nni, nnj, nnjd
!
! Azimuthal direction index
INTEGER         :: jdir
!
! Counter for while loop
INTEGER         :: jwhile
!
!
! Number of spectra with low m intersect
INTEGER         :: nchop2
!
! Number of spiral A solution points
! INTEGER         :: nspira
!
! Minimum level in launchlev(i,j)
INTEGER         :: minlaunchlev
!
! Maximum level in launchlev(i,j)
INTEGER         :: maxlaunchlev
! ----------------------------------------------------------------------+-------
!
! Azimuthal sector for launch spectrum integral Delta Phi / 2
REAL            :: ddphir2
!
! Inertial frequency at current latitude (rad s^-1)
REAL            :: f_f
!
! High wavenumber intersection point
REAL            :: mnly
!
! Wavenumber where Doppler transformed spectrum is reduced to zero
REAL            :: mkill
!
! omega_min(launch) / N (k)
REAL            :: ominrnbv
!
! Constant component of saturation spectrum
REAL            :: ccs0rmstarsq
!
! Minimum range value of function f
REAL            :: fminus
!
! Intermediate  value of function f
REAL            :: fterm
!
! Maximum range value of function f
REAL            :: fplus
!
! Minimum range value of function g
REAL            :: gminus
!
! Intermediate  value of function g
REAL            :: gterm
!
! Maximum range value of function g
REAL            :: gplus
!
! Indicates spectra with low m intersect
LOGICAL         :: L_chop2
! ----------------------------------------------------------------------+-------
! Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!
! Delta U=udotk(launch)-udotk(k)
REAL            :: ddu_a
!
! Record of total flux
REAL            :: fpfac(i_start:i_end,j_start:j_end,idir)
!
! Coefficient A in intersect point equation
REAL            :: acoeff
!
! Term (A / B) in intersect point equation
REAL            :: curvature(i_start:i_end,j_start:j_end,               &
                                           idir)
!
! Coefficient B in intersect point equation
REAL            :: attenuation(i_start:i_end,j_start:j_end,             &
                                             idir)
!
! Chop function B*[1 + (A/B)]^(t-2)
REAL            :: intercept1(i_start:i_end,j_start:j_end,              &
                                            idir)
!
! Chop function B*[1 + (A/B)*mNlmin]^(t-2)
REAL            :: mintercept(i_start:i_end,j_start:j_end,              &
                                            idir)
!
! Starting value of vertical wavenumber for crossing point search
REAL            :: mguess_a(i_start:i_end,j_start:j_end,                &
                                          idir)
!
! Intersect mNlX estimates
REAL            :: mnlx(idir*rows*row_length,0:maxwhile)
!
! I location of chop type points
REAL            :: indexi(idir*rows*row_length)
!
! J location of chop type points
REAL            :: indexj(idir*rows*row_length)
!
! jdir location of chop type points
REAL            :: indexjd(idir*rows*row_length)
!
! Compressed attenuation array
REAL            :: atte_c(idir*rows*row_length)
!
! Compressed curvature array
REAL            :: curv_c(idir*rows*row_length)
!
! Weighting of n term in iter
REAL            :: wgtn(idir*rows*row_length)
!
! Either f_f or the equatorial minimum frequency, whichever is less (rad s^-1)
REAL            :: omin(i_start:i_end,j_start:j_end)
!
! [CGW Flux]_klaunch
REAL            :: cgwll(i_start:i_end, j_start:j_end)
!
! [Rho . Cl]_klaunch
REAL            :: rhocl(i_start:i_end,j_start:j_end,idir)
!
! [Rho(z) . Csat(z) / m*^2]_k
REAL            :: fsatk_scale(i_start:i_end,j_start:j_end,             &
                                             k_start:k_end)    
!
! Indicator for direction of spiral solution
LOGICAL         :: L_ftheng(idir*rows*row_length)
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
CHARACTER(LEN=*), PARAMETER :: RoutineName='GW_USSP_CORE'
!
!  End of Header
!
! ==Main Block==--------------------------------------------------------+-------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------
! Local Constants (a) Physical
! ----------------------------
mstar    = 2.0 * pi / wavelstar
rmstar   = 1.0 / mstar
rmstarsq = rmstar * rmstar
beta_eq_rmstar = 2.3e-11 * rmstar
!
mnlmin   = wavelstar * lminl
fminus = mnlmin**ssptt
tail_chop2b = ssp1 / (mnlmin**ttm1)
head_chop2a = ttm1 * (mnlmin**ssp1)
a0_r_sp1tm1 = 1.0 / ( ss + tt - head_chop2a )
a0_r_1mt    = (-(ssp1) ) * a0_r_sp1tm1
ddphir2 = pi / idir
ccs0rmstarsq = rmstarsq * beta_e0 * SIN(ddphir2) * psatm1 /(pi*twompsat)
!
minlaunchlev = launchlev
maxlaunchlev = launchlev
!
! ----------------------------------------------------------------------+-------
! Core calculations on Physics Grid now within separate subroutine
! ----------------------------------------------------------------------+-------
!
! L_add_cgw = F : calculate standard USSP isotropic GW launch flux
! L_add_cgw = T : calculate variable CGW launch flux
!
IF (L_add_cgw)  THEN
  glob_launch_flux = 0.0
  cgw_convert_factor = cor2flux * cgw_scale_factor
! Option also to adjust Fsat_k scaling
  ccs0rmstarsq = ccs0rmstarsq * 1.0
ELSE
  glob_launch_flux = rmstarsq * ccl0 * ussp_launch_factor
  cgw_convert_factor = 0.0
END IF
!
! ----------------------------------------------------------------------+-------
! 3.0 Initialize gravity wave spectrum variables
! ----------------------------------------------------------------------+-------
!

! ----------------------------------------------------------------------+-------
! 3.1 Compute minimum intrinsic frequency Omin from inertial frequency
!     squared f_f. See UMDP 34, eqn. (A.7).
!-----------------------------------------------------------------------+-------

Rows_do1c: DO j=j_start,j_end
  Row_length_do1c: DO i=i_start,i_end
    f_f = (two_omega * sin_theta_latitude(i,j))**2
    omin(i,j) = nbv(i,j,launchlev) * beta_eq_rmstar
    omin(i,j) = MAX( omin(i,j), f_f )
    omin(i,j) = SQRT(omin(i,j))
!
! ----------------------------------------------------------------------+-------
!   Convection source of flux at Launch
! ----------------------------------------------------------------------+-------
    cgwll(i,j) = totalppn(i,j) * rppnmin
    cgwll(i,j) = MAX( cgwll(i,j), 1.0 )
    cgwll(i,j) = SQRT(cgwll(i,j)) * cgw_convert_factor
  END DO  Row_length_do1c
END DO  Rows_do1c

! ----------------------------------------------------------------------+-------
! 3.2 Compute [rho(z) . C(z)]_k / m*^2 scaling factor for quasi-saturated
!     spectrum (as per ... UMDP 34: 1.15,1.16)
! ----------------------------------------------------------------------+-------
!
! IF (ABS(psatm1) >= 0.1) THEN
! For current setting of parameter psat this test is always true

Levels_do2c: DO k=minlaunchlev,k_end
  Rows_do2c: DO j=j_start,j_end
    Row_length_do2c: DO i=i_start,i_end
      ominrnbv = omin(i,j) / nbv(i,j,k)
      fsatk_scale(i,j,k) = rho_th(i,j,k) * ccs0rmstarsq *               &
       (nbv(i,j,k))**2 * (ominrnbv**psatm1) *                           &
       (1.0 - (ominrnbv**twompsat)) / (1.0 - (ominrnbv**psatm1))
    END DO  Row_length_do2c
  END DO  Rows_do2c
END DO  Levels_do2c

! ELSE
! Require a different functional form for normalisation factor B0
!   BBS = 1.0 / ALOG(nbv(i,j,k) / omin(i,j))
! END IF

!
! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out for each azimuthal direction
! ----------------------------------------------------------------------+-------
!

IDir_do1c: DO jdir=1,idir
!
  Rows_do3c: DO j=j_start,j_end
    Row_length_do3c: DO i=i_start,i_end
! ----------------------------------------------------------------------+-------
!     Initialize gravity wave spectrum variables for the launch level
! ----------------------------------------------------------------------+-------
!
! ----------------------------------------------------------------------+-------
!     Globally invariant value of Total vertical Flux of horizontal
!     pseudomomentum at Launch level, analytic integral under curve 
!     = rho_l * C_l / m*^2 ... UMDP 34: 1.14
! ----------------------------------------------------------------------+-------
      rhocl(i,j,jdir) = rho_th(i,j,launchlev) *                         &
                              (glob_launch_flux + cgwll(i,j))
    END DO  Row_length_do3c
  END DO  Rows_do3c
END DO  IDir_do1c

IDir_do1a: DO jdir=1,idir
  Levels_do4c: DO k=k_start,(minlaunchlev-1)
    Rows_do4c: DO j=j_start,j_end
      Row_length_do4c: DO i=i_start,i_end
! ----------------------------------------------------------------------+-------
!       Set total flux of horizontal pseudomomentum at bottom
! ----------------------------------------------------------------------+-------
        fptot(i,j,k,jdir) = 0.0
      END DO  Row_length_do4c
    END DO  Rows_do4c
  END DO  Levels_do4c
!
  Levels_do5c: DO k=minlaunchlev,k_end
    Rows_do5c: DO j=j_start,j_end
      Row_length_do5c: DO i=i_start,i_end
! ----------------------------------------------------------------------+-------
!     Note: the total flux at levels above is initialised to the sum of 
!     all sources propagated up through the column.
! ----------------------------------------------------------------------+-------
        IF (k == launchlev)  THEN
          fptot(i,j,k,jdir) = rhocl(i,j,jdir)
        ELSE
          fptot(i,j,k,jdir) = 0.0
        END IF
        fptot(i,j,k,jdir) = fptot(i,j,k-1,jdir) + fptot(i,j,k,jdir)
      END DO  Row_length_do5c
    END DO  Rows_do5c
  END DO  Levels_do5c
END DO  IDir_do1a

! ----------------------------------------------------------------------+-------
!     Loop over directions and levels and calculate horizontal component
!     of the vertical flux of pseudomomentum for each azimuthal
!     direction and for each altitude
! ----------------------------------------------------------------------+-------
  Levels_do92: DO k=minlaunchlev+1,k_end
    L_chop2 = .FALSE.
    IDir_do2a: DO jdir=1, idir
      Rows_do92: DO j=j_start,j_end
        Row_length_do92: DO i=i_start,i_end
! ----------------------------------------------------------------------+-------
!       Initialise MGUESS (start point for iterative searches if needed)
! ----------------------------------------------------------------------+-------
          mguess_a(i,j,jdir) = 0.0
          Fptot_if1: IF (fptot(i,j,k-1,jdir) >  0.0) THEN
! ----------------------------------------------------------------------+-------
! 4.2       Calculate variables that define the Chop Type Cases.
! ----------------------------------------------------------------------+-------
            ddu_a = udotk(i,j,launchlev,jdir) -             &
                    udotk(i,j,k,jdir)
! ----------------------------------------------------------------------+-------
!           UMDP 34: 1.23 coefficient B
!           Using ratio of flux scalings rather than densities is more 
!           robust because total flux imported from other schemes may 
!           have different relationship to launch density than 1.14
!           Note: the total flux is initialised to the sum of all sources
!           propagated up through the column.
! ----------------------------------------------------------------------+-------
            attenuation(i,j,jdir) =                                     &
               ( fsatk_scale(i,j,k) / fptot(i,j,k,jdir) ) *             &
               ( nbv(i,j,launchlev) / nbv(i,j,k) )**ttm1
! ----------------------------------------------------------------------+-------
!           UMDP 34: 1.22 coefficient A = (A/B) * B
! ----------------------------------------------------------------------+-------
            curvature(i,j,jdir) = ddu_a * mstar /                       &
                                                nbv(i,j,launchlev)
!
            acoeff = curvature(i,j,jdir) * attenuation(i,j,jdir)
!
            mintercept(i,j,jdir) = attenuation(i,j,jdir) *              &
                      ( (mnlmin * curvature(i,j,jdir)) + 1.0 )**ttm2
!
            Curv_if1: IF (curvature(i,j,jdir) <  asecn)  THEN
! ----------------------------------------------------------------------+-------
!           Negative Doppler Shift : factor will hit zero (kill point)
! ----------------------------------------------------------------------+-------
              mkill = 1.0 / ABS(curvature(i,j,jdir))
!
              Mkill_if1: IF (mkill <= mnlmin )  THEN
! ----------------------------------------------------------------------+-------
!             Chop Type IV : No flux propagates
! ----------------------------------------------------------------------+-------
                fptot(i,j,k,jdir) = 0.0
              ELSE
                IF (mkill >  1.0)  THEN
                  intercept1(i,j,jdir) = attenuation(i,j,jdir) *        &
                                 ( 1.0 + curvature(i,j,jdir) )**ttm2
                ELSE
! ----------------------------------------------------------------------+-------
!               Doppler factor minimum (kill point) situated below
!               mstar in the low-m part of the launch spectrum
! ----------------------------------------------------------------------+-------
                  intercept1(i,j,jdir) = 0.0
                END IF
!
                Lowend_if1: IF (intercept1(i,j,jdir) >= 1.0) THEN
! ----------------------------------------------------------------------+-------
!               Chop Type I: Intersection in high wavenumber part only
! ----------------------------------------------------------------------+-------
                  mnly = ( attenuation(i,j,jdir)**ttrat -               &
                           attenuation(i,j,jdir)) / acoeff
!
                  fptot(i,j,k,jdir) = fptot(i,j,k,jdir) *               &
               (1.0 - (a0_r_1mt * curvature(i,j,jdir) * mnly**twomtt))
                ELSE
                  IF (mintercept(i,j,jdir) <= fminus)  THEN
! ----------------------------------------------------------------------+-------
!                 Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                    fptot(i,j,k,jdir) = fptot(i,j,k,jdir) *             &
                     a0_r_sp1tm1 * tail_chop2b * mintercept(i,j,jdir)   &
                           * ( (mnlmin * curvature(i,j,jdir)) + 1.0 )
                  ELSE
! ----------------------------------------------------------------------+-------
!                 Chop Type IIa: Low wavenumber intersect only
! ----------------------------------------------------------------------+-------
                    L_chop2 = .TRUE.
                    mguess_a(i,j,jdir) = MIN(mkill, 1.0)
                    fpfac(i,j,jdir) = fptot(i,j,k,jdir) * a0_r_sp1tm1
                    fptot(i,j,k,jdir) = 0.0
                  END IF
                END IF  Lowend_if1
!
              END IF  Mkill_if1
!
            ELSE IF (curvature(i,j,jdir) >  asecp)  THEN
! ----------------------------------------------------------------------+-------
!           Positive Doppler Shift : non-zero factor (no kill point)
! ----------------------------------------------------------------------+-------
              intercept1(i,j,jdir) = attenuation(i,j,jdir) *            &
                                 ( 1.0 + curvature(i,j,jdir) )**ttm2
!
              Chop3_if1: IF (intercept1(i,j,jdir) <  1.0)  THEN
! ----------------------------------------------------------------------+-------
!             Chop Type III: Intersection in both wavenumber parts
! ----------------------------------------------------------------------+-------
                fpfac(i,j,jdir) = fptot(i,j,k,jdir) * a0_r_sp1tm1
!
! ----------------------------------------------------------------------+-------
!               First find intersect in high wavenumber part
!               UMDP 34: 1.25
! ----------------------------------------------------------------------+-------
                mnly = ( attenuation(i,j,jdir)**ttrat -                 &
                         attenuation(i,j,jdir)) / acoeff
!
                fptot(i,j,k,jdir) = fptot(i,j,k,jdir) * a0_r_1mt *      &
                                curvature(i,j,jdir) * mnly**twomtt
!
! ----------------------------------------------------------------------+-------
!               Then find intersect in low wavenumber part to reckon
!               its flux contribution for addition when available
! ----------------------------------------------------------------------+-------
                IF (mintercept(i,j,jdir) <= fminus)  THEN
! ----------------------------------------------------------------------+-------
!               Chop Type IIIb: Low wavenumber intersect below min
! ----------------------------------------------------------------------+-------
                  fptot(i,j,k,jdir) = fptot(i,j,k,jdir) +               &
                                    ( fpfac(i,j,jdir) * tail_chop2b *   &
                                 mintercept(i,j,jdir) *                 &
                      ( (mnlmin * curvature(i,j,jdir)) + 1.0 ) )
                ELSE
! ----------------------------------------------------------------------+-------
!               Chop Type IIIa: Low wavenumber intersect
! ----------------------------------------------------------------------+-------
                  L_chop2 = .TRUE.
                  mguess_a(i,j,jdir) = 1.0
                END IF
!
! ----------------------------------------------------------------------+-------
!             ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
              END IF  Chop3_if1
            ELSE
! ----------------------------------------------------------------------+-------
!           Negligible Doppler shift
! ----------------------------------------------------------------------+-------
!             Strictly this is analytic solution mNLX.  UMDP 34: 1.27
              mnly = attenuation(i,j,jdir)**rssptt
              IF (mnly <= mnlmin)  THEN
! ----------------------------------------------------------------------+-------
!             Chop Type IIb: Low wavenumber intersect only below min
! ----------------------------------------------------------------------+-------
                fptot(i,j,k,jdir) = fptot(i,j,k,jdir) * a0_r_sp1tm1 *   &
                                 tail_chop2b * attenuation(i,j,jdir)
              ELSE
                IF (mnly <  1.0)  fptot(i,j,k,jdir) =                   &
! ----------------------------------------------------------------------+-------
!               Chop Type IIc: Low wavenumber intersect only (analytic)
! ----------------------------------------------------------------------+-------
                  fptot(i,j,k,jdir) * a0_r_sp1tm1 *                     &
                           ( (ssptt * (mnly**ssp1)) - head_chop2a )
! ----------------------------------------------------------------------+-------
!               ELSE Chop Type 0: No intersection (spectrum unaltered)
! ----------------------------------------------------------------------+-------
              END IF
            END IF  Curv_if1
!
          END IF  Fptot_if1
        END DO  Row_length_do92
      END DO  Rows_do92
    END DO  IDir_do2a
!
    Lchop2_if1: IF (L_chop2)  THEN
! ----------------------------------------------------------------------+-------
!   Process low wavenumber contribution: evaluate intersect mNX
! ----------------------------------------------------------------------+-------
      nchop2 = 0
!
      IDir_do2b: DO jdir=1,idir
        Rows_do93: DO j=j_start,j_end
          Row_length_do93: DO i=i_start,i_end
            IF (mguess_a(i,j,jdir) >  0.0)  THEN
              nchop2 = nchop2 + 1
!
              indexjd(nchop2) = jdir
              indexj(nchop2)  = j
              indexi(nchop2)  = i
            END IF
          END DO  Row_length_do93
        END DO  Rows_do93
      END DO  IDir_do2b
!     nspira = 0
!
      Nchop2_do1: DO i=1,nchop2
! ----------------------------------------------------------------------+-------
!     Chop Type IIa : / Full solution required for mNlX
! or  Chop Type IIIa: \
! ----------------------------------------------------------------------+-------
        nnjd = indexjd(i)
        nnj  = indexj(i)
        nni  = indexi(i)
!
        fplus  = mguess_a(nni,nnj,nnjd)**ssptt
        gplus  = intercept1(nni,nnj,nnjd)
!       fminus = mnlmin**ssptt    Defined as a constant
        gminus = mintercept(nni,nnj,nnjd)
        atte_c(i) = attenuation(nni,nnj,nnjd)
        curv_c(i) = curvature(nni,nnj,nnjd)
!
        fterm = ( ((fminus / atte_c(i))**rttm2) - 1.0 ) / curv_c(i)
        gterm = gminus**rssptt
        L_ftheng(i) = .FALSE.
!
        Curv_if2: IF (curvature(nni,nnj,nnjd) >  asecp)  THEN
! ----------------------------------------------------------------------+-------
!       Positive Doppler Shift
! ----------------------------------------------------------------------+-------
          wgtn(i) = 0.0
        ELSE
! ----------------------------------------------------------------------+-------
!       Negative Doppler Shift
! ----------------------------------------------------------------------+-------
          wgtn(i) = 1.0 - mweight
!
          IF (fplus <= gminus  .AND.  gplus >  fminus) THEN
            fterm = (((fplus / atte_c(i))**rttm2) - 1.0)/ curv_c(i)
            gterm = gplus**rssptt
            L_ftheng(i) = (gterm <  fterm)
!
          ELSE IF (fplus >  gminus  .AND.  gplus <= fminus) THEN
            L_ftheng(i) = (gterm >= fterm)
!
          ELSE IF (fplus <= gminus  .AND.  gplus <= fminus) THEN
            L_ftheng(i) = .TRUE.
!
!         ELSE Use default settings
          END IF
        END IF  Curv_if2
!
        IF (L_ftheng(i))  THEN
!         nspira = nspira + 1
          mnlx(i,0) = fterm
        ELSE
          mnlx(i,0) = gterm
        END IF
      END DO  Nchop2_do1
!
      Jwhile_do2: DO jwhile=0,maxwhile-1
        Nchop2_do2: DO i=1,nchop2
!
          IF (L_ftheng(i)) THEN
! ----------------------------------------------------------------------+-------
!         Obtain m_n+1 from g_n+1  = f_n (m_n)
! ----------------------------------------------------------------------+-------
            mnlx(i,jwhile+1) = (                                        &
                 (((mnlx(i,jwhile)**ssptt) / atte_c(i))**rttm2) - 1.0 ) &
                 / curv_c(i)
          ELSE
! ----------------------------------------------------------------------+-------
!         Obtain m_n+1 from f_n+1  = g_n (m_n)
! ----------------------------------------------------------------------+-------
            mnlx(i,jwhile+1) = ( (atte_c(i) *                           &
                ((1.0 + (curv_c(i) * mnlx(i,jwhile)))**ttm2))**rssptt )
          END IF
!
          mnlx(i,jwhile+1) = ((1.0 - wgtn(i)) * mnlx(i,jwhile+1)) +     &
                                    (wgtn(i)  * mnlx(i,jwhile))
!
        END DO  Nchop2_do2
      END DO  Jwhile_do2
!
!CDIR NODEP
      Nchop2_do3: DO i=1,nchop2
        nnjd = indexjd(i)
        nnj  = indexj(i)
        nni  = indexi(i)
!
        fptot(nni,nnj,k,nnjd) = fptot(nni,nnj,k,nnjd) + (               &
        fpfac(nni,nnj,nnjd) * ( ((mnlx(i,maxwhile)**ssp1) * ( ssptt +   &
        (ssp1 * mnlx(i,maxwhile) * curv_c(i)) ) ) - head_chop2a ) )
!
      END DO  Nchop2_do3
    END IF  Lchop2_if1
!
    IDir_do2c: DO jdir=1,idir
      Rows_do10: DO j=j_start,j_end
        Row_length_do10: DO i=i_start,i_end
!-----------------------------------------------------------------------+-------
!       Now correct pseudomomentum flux in the evolved spectrum if the
!       new value is non-physical (pseudomomentum flux cannot increase
!       with altitude)
!-----------------------------------------------------------------------+-------
          fptot(i,j,k,jdir) = MIN(fptot(i,j,k,jdir),fptot(i,j,k-1,jdir))
!
        END DO  Row_length_do10
      END DO  Rows_do10
    END DO  IDir_do2c
!
  END DO  Levels_do92
!
! ----------------------------------------------------------------------+-------
! 4.5   Apply Opaque Upper Boundary to set fluxes to zero at top
! ----------------------------------------------------------------------+-------
IDir_do6c: DO jdir=1,idir
  Rows_do12: DO j=j_start,j_end
    Row_length_do12: DO i=i_start,i_end
      fptot(i,j,k_end,jdir) =  0.0
    END DO  Row_length_do12
  END DO  Rows_do12
END DO  IDir_do6c

!
! ----------------------------------------------------------------------+-------
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
!
END SUBROUTINE gw_ussp_core

END MODULE gw_ussp_core_mod
