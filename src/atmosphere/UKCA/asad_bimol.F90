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
! Purpose: To calculate bimolecular rate coefficients using
!          data from the bimolecular ratefile ratb.d or
!          ukca_chem1 module
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
!     Method
!     ------
!     Bimolecular rate coefficients are calculated from the Ahrenius
!     expression: k(T) = k * (T/300)^a * exp(-b/T) where T is the
!     temperature (in K). The parameters k, a and b are taken from
!     the bimolecular ratefile ratb.d.
!
!     The reactions CO + OH -> H + CO2 and OH + HONO2 are pressure
!     or density dependent and therefore need to be calculated
!     separately. However, this code should never need changing even
!     if neither of these reaction are included in the chemistry.
!
!     The reactions HO2+MeCO3, MeOO+MeOO and OH+C3H8 have temperature
!     dependent branching ratios. Therefore, they are calculated
!     separately.
!
!     Local Variables
!     ---------------
!     ih2o           Tracer index for H2O if it's an advective tracer
!     iohco          Reaction index for OH + CO.
!     iohhno3        Reaction index for OH + HONO2.
!     iohc3h8a       Reaction index for OH+C3H8 -> n-PrOO + H2O
!     iohc3h8b       Reaction index for OH+C3H8 -> i-PrOO + H2O
!     imeoomeooa     Reaction index for MeOO+MeOO -> HO2+HCHO+HO2+HCHO
!     imeoomeoob     Reaction index for MeOO+MeOO -> MeOH+HCHO
!     imeooho2a      Reaction index for MeOO+HO2 -> MeOOH
!     imeooho2b      Reaction index for MeOO+HO2 -> HCHO
!     idmsoh         Reaction index for DMS+OH  -> 0.6SO2 + 0.4DMSO + MeOO
!     iso3h2o        Reaction index for SO3 + H2O -> H2SO4 + H2O
!     ics2oh         Reaction index for CS2 + OH -> COS + SO2
!     z1,z3,z4       Variables used to calculate pressure
!                    dependent rate coefficients.
!     ratioa2b       Branching ratio of branch A to branch B
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_bimol_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_BIMOL_MOD'

CONTAINS

SUBROUTINE asad_bimol( n_points )

USE asad_mod,        ONLY: t, t300, advt, spb, ab, rk, tnd, p,    &
                           wp, f, peps, nbrkx, jpspb
USE asad_findreaction_mod, ONLY: asad_findreaction
USE ukca_option_mod, ONLY: jpctr, jpbk, l_ukca_offline, l_ukca_offline_be,   &
                           l_ukca_chem_aero, l_ukca_strat, l_ukca_strattrop, &
                           l_ukca_tropisop
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE vectlib_mod, ONLY: exp_v, powr_v, oneover_v
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

! Local variables

INTEGER, SAVE :: iohco = 0
INTEGER, SAVE :: iohhno3 = 0
INTEGER, SAVE :: iho2 = 0
INTEGER, SAVE :: iohc3h8a = 0
INTEGER, SAVE :: iohc3h8b = 0
INTEGER, SAVE :: imeoomeooa = 0
INTEGER, SAVE :: imeoomeoob = 0
INTEGER, SAVE :: imeooho2a = 0
INTEGER, SAVE :: imeooho2b = 0
INTEGER, SAVE :: idmsoh = 0
INTEGER, SAVE :: iso3h2o = 0
INTEGER, SAVE :: ics2oh = 0
INTEGER, SAVE :: ih2o = 0
INTEGER, SAVE :: io2  = 0
INTEGER, SAVE :: iho2no = 0
INTEGER :: jtr
INTEGER :: j
INTEGER :: jr

REAL, ALLOCATABLE :: z1(:)                ! Intermediate result
REAL, ALLOCATABLE :: z3(:)                !          "
REAL, ALLOCATABLE :: z4(:)                !          "
REAL, ALLOCATABLE :: z5(:)                !          "
REAL, ALLOCATABLE :: z6(:)                !          "
REAL, ALLOCATABLE :: alpha(:)             ! Multiplication factor
REAL, ALLOCATABLE :: ratioa2b(:)          ! Branching ratio
REAL, ALLOCATABLE :: ratiob2total(:)      ! Branching ratio

CHARACTER(LEN=10) :: r1
CHARACTER(LEN=10) :: r2
CHARACTER(LEN=10) :: prods(jpspb)

LOGICAL, SAVE :: first = .TRUE.
LOGICAL, SAVE :: first_pass = .TRUE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_BIMOL'

REAL :: tmp(1:n_points)
REAL :: tmp_out(1:n_points)
REAL :: tmp1(1:n_points)
REAL :: inv_t(1:n_points)

! ErrorStatus
INTEGER                             :: errcode=0   ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength)   :: cmessage    ! Error return message

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(z1(1:n_points))
ALLOCATE(z3(1:n_points))
ALLOCATE(z4(1:n_points))
ALLOCATE(z5(1:n_points))
ALLOCATE(z6(1:n_points))
ALLOCATE(ratioa2b(1:n_points))
ALLOCATE(ratiob2total(1:n_points))

z1(:) = 0.0
z3(:) = 0.0
z4(:) = 0.0
z5(:) = 0.0
z6(:) = 0.0
ratioa2b(:) = 0.0
ratiob2total(:) = 0.0

!       1. Calculate bimolecular rate coefficients
!          --------- ----------- ---- ------------

!       Compute intermediate results

t300(1:n_points) = t(1:n_points) / 300.0
CALL oneover_v(n_points, t, inv_t)

!       Check if H2O is an advected tracer

IF (first_pass) THEN
  ! OMP CRITICAL will only allow one thread through this code at a time,
  ! while the other threads are held until completion.
!$OMP CRITICAL (asad_bimol_init)
  IF (first) THEN
    first = .FALSE.
    DO jtr = 1, jpctr
      IF ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
      IF ( advt(jtr)  ==  'O2        ' ) io2 = jtr
    END DO

    !       Look for the reactions which need special treatment.

    r1 = '          '
    r2 = r1
    prods(:) = r1
    r1 = 'OH        '
    r2 = 'C3H8      '
    prods(1) = 'n-PrOO    '
    prods(2) = 'H2O       '
    iohc3h8a = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,    &
                                  jpbk+1, jpspb )
    prods(1) = 'i-PrOO    '
    prods(2) = 'H2O       '
    iohc3h8b = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,    &
                                  jpbk+1, jpspb )
    r1 = 'MeOO      '
    r2 = 'MeOO      '
    prods(1) = 'MeOH      '
    prods(2) = 'HCHO      '
    imeoomeooa = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,  &
                                    jpbk+1, jpspb )

    prods(1) = 'HO2       '
    prods(2) = 'HCHO      '
    prods(3) = 'HO2       '
    prods(4) = 'HCHO      '
    imeoomeoob = asad_findreaction( r1, r2, prods, 4, spb, nbrkx,    &
                                    jpbk+1, jpspb )
    r1 = 'HO2       '
    r2 = 'MeOO      '
    prods(1) = 'MeOOH     '
    imeooho2a = asad_findreaction( r1, r2, prods, 1, spb, nbrkx,     &
                                   jpbk+1, jpspb )
    prods(1) = 'HCHO      '
    imeooho2b = asad_findreaction( r1, r2, prods, 1, spb, nbrkx,     &
                                   jpbk+1, jpspb )

    r1 = 'HO2       '
    r2 = 'NO        '
    prods(1) = 'HONO2     '
    iho2no    = asad_findreaction( r1, r2, prods, 1, spb, nbrkx,  &
                                   jpbk+1, jpspb )

    r1 = 'DMS       '
    r2 = 'OH        '
    prods(1) = 'SO2       '
    ! The products of this reaction vary with chemical scheme
    IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
      !   DMS + OH => SO2 + DMSO
      prods(2) = 'DMSO      '
      idmsoh = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,           &
                                  jpbk+1, jpspb )
!!    ELSE IF (l_ukca_strat .OR. l_ukca_strattrop) THEN !Use original DMS scheme, LER, UC 19.02.19
    ELSE IF (l_ukca_strat) THEN !!LER, UC Feb 2019 (to use Chen scheme)
      !   DMS + OH => SO2 + MSA
      prods(2) = 'MSA       '
      idmsoh = asad_findreaction( r1, r2, prods, 2, spb, nbrkx,           &
                                  jpbk+1, jpspb )
      PRINT*,'TEST1'
!      PRINT('TEST1')
!!    ELSE IF (l_ukca_tropisop) THEN !Use original DMS scheme, LER, UC 19.02.19
    ELSE IF (l_ukca_tropisop .OR. l_ukca_strattrop) THEN !!LER, UC Feb 2019 (to use Chen scheme)
      !   DMS + OH => SO2 + DMSO + MeOO
      prods(2) = 'DMSO      '
      prods(3) = 'MeOO      '
      idmsoh = asad_findreaction( r1, r2, prods, 3, spb, nbrkx,           &
                                  jpbk+1, jpspb )
      PRINT*,'TEST2'
    ELSE
      idmsoh = 0
    END IF
    IF ((l_ukca_chem_aero .AND. (l_ukca_offline .OR. l_ukca_offline_be .OR. &
         l_ukca_strat .OR. l_ukca_strattrop .OR. l_ukca_tropisop))          &
         .AND. idmsoh == 0) THEN
      errcode = 1
      cmessage = ' Rate coefficient for DMS + OH lacks correction'
      CALL ereport('ASAD_BIMOL',errcode,cmessage)
    END IF

    DO j = 1, jpbk
      jr = nbrkx(j)

      IF ( ( spb(j,1) == 'OH     ' .AND. spb(j,2) == 'CO     ' ) .OR.  &
         ( spb(j,1) == 'CO     ' .AND. spb(j,2) == 'OH     ' ) )       &
        iohco = jr
      IF ( ( spb(j,1) == 'OH     ' .AND. spb(j,2) == 'HONO2  ' ) .OR.  &
         (  spb(j,1) == 'HONO2  ' .AND. spb(j,2) == 'OH     ' ) )      &
        iohhno3 = jr
      IF ( spb(j,1) == 'HO2    ' .AND. spb(j,2) == 'HO2    ' )        &
        iho2 = jr

      ! code for stratospheric sulphur scheme
      IF ( spb(j,1) == 'SO3      ' .AND. spb(j,2) == 'H2O   ')    &
        iso3h2o = jr
      IF ( spb(j,1) == 'CS2      ' .AND. spb(j,2) == 'OH    ')    &
        ics2oh = jr

    END DO  ! end of loop (j) over jpbk
  END IF    ! first
  first_pass = .FALSE.
!$OMP END CRITICAL (asad_bimol_init)
END IF      ! first_pass

!       1.2  Compute rates

DO j = 1, jpbk
  jr = nbrkx(j)

  IF ( ABS(ab(j,2)) < peps .AND. ABS(ab(j,3)) < peps ) THEN
    rk(1:n_points,jr) = ab(j,1)
  ELSE IF ( ABS(ab(j,2)) < peps ) THEN
    tmp(1:n_points) = -ab(j,3)*inv_t(1:n_points)
    CALL exp_v(n_points,tmp,tmp_out)
    rk(1:n_points,jr) = ab(j,1) * tmp_out(1:n_points)
  ELSE IF ( ABS(ab(j,3)) < peps ) THEN
    CALL powr_v(n_points,t300,ab(j,2),tmp1)
    rk(1:n_points,jr) = ab(j,1) * tmp1(1:n_points)
  ELSE
    tmp(1:n_points) = -ab(j,3)*inv_t(1:n_points)
    CALL exp_v(n_points,tmp,tmp_out)
    CALL powr_v(n_points,t300,ab(j,2),tmp1)
    rk(1:n_points,jr) = ab(j,1) * tmp1(1:n_points) * tmp_out(1:n_points)
  END IF
END DO  ! end of loop (j) over jpbk

!       2. Dependent reactions.
!          --------- ----------

! OH + CO; updated with IUPAC March 2005 (Paul Young)
! GC: note the original March IUPAC summary had a mistake for this
! reaction. The values given in the datasheet are correct.
! k = k' . (1 + [N2]/4.2E19).. we use TND below instead of [N2] but
! Paul suggests the reaction would probably go with [O2] anyway.

IF ( iohco /= 0 ) THEN
  rk(1:n_points,iohco)=rk(1:n_points,iohco)*                      &
                  (1.0 + tnd(1:n_points)/4.2e19)
END IF

! OH + HONO2; no change with IUPAC Jan 2009 (CJ)
IF ( iohhno3 /= 0 ) THEN
  tmp(1:n_points) = 460.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  z1(:) = 2.4e-14 * tmp_out(1:n_points)
  tmp(1:n_points) = 1335.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  z3(:) = (6.5e-34 * tmp_out(1:n_points)) * tnd(1:n_points)
  tmp(1:n_points) = 2199.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  z4(:) = 2.7e-17 * tmp_out(1:n_points)
  rk(1:n_points,iohhno3) = z1(:) + z3(:)/(1.0+z3(:)/z4(:))
END IF

! HO2 + HO2; no change with IUPAC Nov 2003 (Paul Young)
IF ( ih2o /= 0 .AND. iho2 /= 0) THEN
  !       water is an advected tracer
  tmp(1:n_points) = 2200.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  rk(1:n_points,iho2) = rk(1:n_points,iho2) *                 &
    (1.0+1.4e-21*f(1:n_points,ih2o)*tmp_out(1:n_points))
ELSE IF (ih2o == 0 .AND. iho2 /= 0) THEN
  !       use model water concentration
  tmp(1:n_points) = 2200.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  rk(1:n_points,iho2) = rk(1:n_points,iho2) *                 &
    (1.0+1.4e-21*wp(1:n_points)*tnd(1:n_points)*              &
     tmp_out(1:n_points))
END IF

! HO2 + NO -> HONO2 (with extra temp and pressure dependence)
! Added by Alex 2012
IF ( iho2no /= 0 ) THEN
  rk(1:n_points,iho2no)=rk(1:n_points,iho2no)*                    &
   ((530.0*inv_t(1:n_points)) +                 &
    8.53e-4*(1e-2*p(1:n_points))-1.73)/100.0
END IF

! SO3 + H2O: 2nd H2O molecule dealt with here by multiplying rate by [H2O]
IF ( ih2o /= 0 .AND. iso3h2o /= 0 ) THEN
  !       water is an advected tracer
  rk(1:n_points,iso3h2o) = rk(1:n_points,iso3h2o) *               &
       f(1:n_points,ih2o)
ELSE IF (ih2o == 0 .AND. iso3h2o /= 0) THEN
  !       use model water concentration
  rk(1:n_points,iso3h2o) = rk(1:n_points,iso3h2o) *               &
       wp(1:n_points)*tnd(1:n_points)
END IF


! CS2 + OH rate from JPL evaluation 15.
IF ( ics2oh /= 0 ) THEN
  tmp(1:n_points) = 3400.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  z6(:) = 1.81e-3 * tmp_out(1:n_points)
  rk(1:n_points,ics2oh) =                                         &
        rk(1:n_points,ics2oh)/(t(1:n_points)+z6(:))
END IF

! DMS + OH: Multiply by factor alpha/(1 + alpha) from Pham et al. (1995)
IF ( idmsoh /= 0) THEN
  ALLOCATE(alpha(1:n_points))
  tmp(1:n_points) = 7460.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  alpha(:) = 1.106e-31*tmp_out(1:n_points)* tnd(1:n_points)
  rk(1:n_points,idmsoh) = rk(1:n_points,idmsoh)*                  &
                          alpha(:)/(1.0 + alpha(:))
  ! IUPAC, 2011
  !            rk(1:n_points,idmsoh) = rk(1:n_points,idmsoh)*             &
  !                                    f(1:n_points,io2)/                 &
  !              (1 + 7.5E-29*EXP(5610/t(1:n_points))*f(1:n_points,io2))
END IF

!       3. Temperature-Dependent branching ratios
!          ----------- --------- --------- ------
!       rk above was calculated using the total rate coefficients.
!       Here, rk is reduced according to the branching ratio.


!  OH + C3H8 -> n-PrOO ... Branch A
!  OH + C3H8 -> i-PrOO ... Branch B
IF (iohc3h8a /= 0 .AND. iohc3h8b /=0) THEN
  tmp(1:n_points) = -816.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  CALL powr_v(n_points,t,(-0.64),tmp1)
  ratioa2b(1:n_points) = 226.0 * tmp1(1:n_points)* tmp_out(1:n_points)
  rk(1:n_points,iohc3h8a) = rk(1:n_points,iohc3h8a) *             &
               (ratioa2b(1:n_points)/(ratioa2b(1:n_points)+1.0))
  rk(1:n_points,iohc3h8b) = rk(1:n_points,iohc3h8b)/              &
                              (ratioa2b(1:n_points)+1.0)
END IF

!  MeOO + MeOO -> MeOH + HCHO       ... Branch A
!  MeOO + MeOO -> 2HO2 + 2HCHO      ... Branch B
IF (imeoomeooa /= 0 .AND. imeoomeoob /=0) THEN
  tmp(1:n_points) = 1300.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  ratiob2total(1:n_points) = 1.0/(1.0+ tmp_out(1:n_points)/33.0)
  rk(1:n_points,imeoomeooa) = rk(1:n_points,imeoomeooa)*          &
                               (1.0-ratiob2total(1:n_points))
  rk(1:n_points,imeoomeoob) = rk(1:n_points,imeoomeoob)*          &
                               (ratiob2total(1:n_points))
END IF

!  Added from IUPAC
!  MeOO + HO2 -> MeOOH     ... Branch A
!  MeOO + HO2 -> HCHO      ... Branch B

IF ( imeooho2a /= 0 .AND. imeooho2b /= 0 ) THEN
  tmp(1:n_points) = -1160.0*inv_t(1:n_points)
  CALL exp_v(n_points,tmp,tmp_out)
  ratiob2total(1:n_points) = 1.0 / (1.0 + 498.0* tmp_out(1:n_points))
  rk(1:n_points,imeooho2a) = rk(1:n_points,imeooho2a)*            &
                               (1.0-ratiob2total(1:n_points))
  rk(1:n_points,imeooho2b) = rk(1:n_points,imeooho2b)*            &
                               (ratiob2total(1:n_points))
END IF


IF (ALLOCATED(z1)) DEALLOCATE(z1)
IF (ALLOCATED(z3)) DEALLOCATE(z3)
IF (ALLOCATED(z4)) DEALLOCATE(z4)
IF (ALLOCATED(z5)) DEALLOCATE(z5)
IF (ALLOCATED(z6)) DEALLOCATE(z6)
IF (ALLOCATED(ratioa2b)) DEALLOCATE(ratioa2b)
IF (ALLOCATED(ratiob2total)) DEALLOCATE(ratiob2total)
IF (ALLOCATED(alpha)) DEALLOCATE(alpha)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_bimol
END MODULE asad_bimol_mod
