! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to calculate NOx lightning emissions for one
!          vertical column based on cloud height above surface,
!          latitude & land/ocean surface.
!
!          Based on light.F from Cambridge TOMCAT model
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:  Cloud height above surface and surface type yield lightning
!          flash frequency. Latitude yields CG/CC ratio (Price & Rind
!          1993), assuming that all dH (zero degrees to cloud top) lie
!          between 5.5-14km. Convert flashes to NOx production for the
!          1/2 hr dynamical period.
!
!          Called from UKCA_LIGHT_CTL.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE ukca_light_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_LIGHT_MOD'

CONTAINS

SUBROUTINE ukca_light(delta_lambda,delta_phi,ppress,niv, &
                      hkmb,hkmt,klc,klt,adlat,           &
                      asfaera,asurf,anox,                &
                      total_flash_rate,                  &
                      cloud2ground_flash_rate,           &
                      cloud2cloud_flash_rate,            &
                      total_N)

USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
USE parcons_mod,     ONLY: rad, deg
USE ukca_constants,  ONLY: avc
USE ukca_option_mod, ONLY: l_ukca_linox_scaling

IMPLICIT NONE

INTEGER, INTENT(IN) :: niv ! No of vertical level
INTEGER, INTENT(IN) :: klc ! Level of cloud base
INTEGER, INTENT(IN) :: klt ! Level of cloud top
INTEGER, INTENT(IN) :: asurf ! Land (1) / sea (0) m

REAL, INTENT(IN) ::  delta_lambda ! gridbox width (radiants)
REAL, INTENT(IN) ::  delta_phi ! gridbox height (radiants)
REAL, INTENT(IN) ::  hkmt ! Height of cloud top
REAL, INTENT(IN) ::  hkmb ! Height of cloud base
REAL, INTENT(IN) ::  adlat ! Latitude
REAL, INTENT(IN) ::  ppress(niv) ! Pressures at model l
REAL, INTENT(IN) ::  asfaera ! Surf area * (radius

!       ...NOx lightning emissions
REAL, INTENT(OUT) :: anox(niv) ! kg(N)/gridcell/s

!       ...number of flashes in a gridcell /min
REAL, INTENT(OUT) :: total_flash_rate

!       ...number of flashes in a gridcell cloud to ground /min
REAL, INTENT(OUT) :: cloud2ground_flash_rate

!       ...number of flashes in a gridcell cloud to cloud /min
REAL, INTENT(OUT) :: cloud2cloud_flash_rate

!       ...lighting N column density in kg(N)/m^2/s
REAL, INTENT(OUT) :: total_N

! Local variables

INTEGER :: jniv ! Loop variable
INTEGER :: k ! Loop variable

!       ...Minimum cloud depth
REAL, PARAMETER ::  Min_clouddepth = 5.0

!       ...distance (km) per degree at equator
REAL, PARAMETER ::  km_per_deg_eq  = 111.11


!       ...NOx production parameters

!       Price et al., J.G.R., 1997 (1 and 2)
!       ...avg energy cld-2-grd flash (J)
!        REAL, PARAMETER ::  E_cg = 6.7E+09
!       ...avg energy cld-2-cld flash (J) (eq 0.1*E_cg)
!        REAL, PARAMETER ::  E_cc = 6.7E+08
!       ...avg NO production rate (molecs(NO) J-1)
!        from Allen and Pickering, JGR, 2002
!        REAL, PARAMETER ::  P_no = 1.0E+17

!       from review of Schumann and Huntrieser, ACP, 2007
!       ...avg energy cld-2-grd flash (J)
REAL, PARAMETER ::  E_cg = 3.0e+09
!       ...avg energy cld-2-cld flash (J)
REAL, PARAMETER ::  E_cc = 0.9e+09
!       ...avg NO production rate (molecs(NO) J-1)
REAL, PARAMETER ::  P_no = 25.0e+16
!      
REAL, PARAMETER ::  N_nitrogen = 1.0e+26

!       ...Molecular mass of N (kg/mol)
REAL, PARAMETER ::  M_n  = 14.01*1.0e-03

REAL :: aflash ! Flash frequency (flashes/min)
REAL :: adh
REAL :: az ! Cloud-cloud/cloud-ground
REAL :: ap ! Cloud-ground flashes / total fl
REAL :: acgfrac ! Cloud-ground flash frequency (f
REAL :: accfrac ! Cloud-cloud flash frequency (fl
REAL :: acgnox
REAL :: accnox
REAL :: dpcg
REAL :: dpcc

REAL :: gb_area_30N ! area of gridbox at 30N (normalisation factor)
REAL :: ew_res_deg  ! gridbox width in degrees
REAL :: ns_res_deg  ! gridbox height in degrees
REAL :: fr_calib_fac ! model resolution calibration factor

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_LIGHT'


!       Initialise variables

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

gb_area_30N  = 0.0
ew_res_deg   = 0.0
ns_res_deg   = 0.0
fr_calib_fac = 0.0

aflash  = 0.0
adh     = 0.0
az      = 0.0
ap      = 0.0
acgfrac = 0.0
accfrac = 0.0
acgnox  = 0.0
accnox  = 0.0
! Initialise to zero and will be value if minimum cloud depth not met.
anox    = 0.0
total_flash_rate = 0.0
cloud2ground_flash_rate = 0.0
cloud2cloud_flash_rate = 0.0
total_N = 0.0

!       ...Calculate resolution-dependent calibration factor
!          spatial calibration factor (N96);
!          cf., Price and Rind, Mon. Weather Rev., 1994 and
!          Allen and Pickering, JGR, 2002, pp. 7-8
ew_res_deg   = delta_lambda*deg
ns_res_deg   = delta_phi*deg
fr_calib_fac = 0.97241*EXP(0.048203*ew_res_deg*ns_res_deg)

!       ...compute gridbox area at 30N (normalisation factor)
gb_area_30N = (ew_res_deg*km_per_deg_eq*0.87)* & ! gridbox width (km)
              (ns_res_deg*km_per_deg_eq)*      & ! gridbox height (km)
              1.0e+06                            ! (km^2 --> m^2)

!       Set minimum cloud depth to 5 km

IF ((hkmt-hkmb) > Min_clouddepth) THEN

  !         ...flashes per minute
  IF (asurf == 0) THEN                  ! Ocean
    aflash = 6.40e-04*(hkmt**1.73)      ! Ocean flash frequency (1/min)
  ELSE                                  ! Land
    aflash = 3.44e-05*(hkmt**4.9)       ! Land flash frequency  (1/min)
  END IF

  !         ...Calculate flash rate in flashes/s/gridbox
  !            note: c.f., Allen and Pickering, J.G.R., 2002

  !         ...calibrate for model resolution
  aflash = aflash*fr_calib_fac*(asfaera/gb_area_30N)

  !         ...convert from flases/gridbox/min to flashes/gridbox/s
  aflash = aflash/60.0

  !         ...work out proportion of flashes that are cloud to ground
  adh = -6.64e-05*(ABS(adlat)**2)-4.73e-03*ABS(adlat)+7.34
  
  az  = 0.021*(adh**4)-0.648*(adh**3)+7.493*(adh**2) &
        - 36.54*adh+63.09
  
  ap  = 1.0/(1.0+az)
  acgfrac = aflash*ap
  accfrac = aflash-acgfrac

  ! convert flashrates back from num/s to num/min as required for diagnostics
  total_flash_rate         = aflash*60.0
  cloud2ground_flash_rate  = acgfrac*60.0
  cloud2cloud_flash_rate   = accfrac*60.0

  !         ...compute NO production in kg(N)/gridbox/s based on flash frequency
  IF (l_ukca_linox_scaling) THEN
     !         ...cloud-to-ground NOx
     acgnox = (acgfrac*N_nitrogen*M_n)/avc
     !         ...cloud-to-cloud NOx
     accnox = (accfrac*N_nitrogen*M_n)/avc
  ELSE 
     !         ...cloud-to-ground NOx
     acgnox = (acgfrac*E_cg*P_no*M_n)/avc
     !         ...cloud-to-cloud NOx
     accnox = (accfrac*E_cc*P_no*M_n)/avc
  END IF

  !         ...total lightning NOx column density in kg(N)/m^2/s
  total_N = (acgnox + accnox)/asfaera


  !         Distribute over the column with each box having same vmr in an
  !         multiply by CONVFAC, conversion factor that gives emissions if
  !         were 100 flashes s**-1

  !         Work out which pressure is closest to 500 hPa

  loop: DO jniv = niv,1,-1
    IF (ppress(jniv) >= 50000.0) EXIT loop
  END DO loop

  !         3 from cloud base to 2 above top
  !         KLT is the level above cloud top

  ! sanity check to prevent dpcg==0
  IF (jniv <= 1) jniv=2
    
  IF (l_ukca_linox_scaling) THEN
    !  DO EVERYTHING LINEARLY IN LOG(PRESSURE)

    dpcg = LOG(ppress(1)) - LOG(ppress(jniv))
    dpcc = LOG(ppress(jniv)) - LOG(ppress(klt))   
    !         ...construct L-NOx profile in kg(N)/gridcell/s
    !         ...first cloud-to-ground L-NOx profiles (kg(N)/gridcell/s)
    IF ((jniv-1) == 1) THEN
       anox(1) = acgnox
    ELSE
      DO k = 1,jniv-1
         anox(k) = acgnox * ((LOG(ppress(k))-LOG(ppress(k+1)))/dpcg) 
      END DO
    END IF
    
    !         ...then cloud-to-cloud L-NOx profiles (kg(N)/gridcell/s)
    IF (LOG(ppress(jniv)) <= LOG(ppress(klt))) THEN
      ! jniv (level of the 500hPa level) is above the
      ! cloud-top-height. In this case, put all C2C N
      ! into the cloud-top level.
      anox(klt-1) = anox(klt-1) + accnox
    ELSE
      ! jniv is greater than the cloud-top-height
      ! Note: anox(k) is also on the RHS of this equation
      DO k = jniv,klt-1
         anox(k) = anox(k) + accnox *                                 &
              ( (LOG(ppress(k)) - LOG(ppress(k+1))) / dpcc ) 
      END DO
    END IF
  ELSE ! .not. l_ukca_linox_scaling
    !  DO EVERYTHING LINEARLY IN PRESSURE
    dpcg = ppress(1) - ppress(jniv)
    dpcc = ppress(jniv) - ppress(klt)
    !         ...construct L-NOx profile in kg(N)/gridcell/s
    !         ...first cloud-to-ground L-NOx profiles (kg(N)/gridcell/s)
    IF ((jniv-1) == 1) THEN
      anox(1) = acgnox
    ELSE
      DO k = 1,jniv-1
         anox(k) = acgnox * ((ppress(k)-ppress(k+1))/dpcg)
      END DO
    END IF
    
    !         ...then cloud-to-cloud L-NOx profiles (kg(N)/gridcell/s)
    IF (ppress(jniv) <= ppress(klt)) THEN
      anox(klt-1) = anox(klt-1) + accnox
    ELSE
      DO k = jniv,klt-1
         anox(k) = accnox * ((ppress(k)-ppress(k+1))/dpcc)
      END DO
    END IF
  END IF ! l_ukca_linox_scaling

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_light
END MODULE ukca_light_mod
