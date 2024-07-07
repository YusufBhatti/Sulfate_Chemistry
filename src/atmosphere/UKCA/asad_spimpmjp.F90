! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution.
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!   Called from asad_spmjpdriv
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine via the mjpdriv driver.
!
!     Method.
!     -------
!     Solves the non-linear system of equations via a Newton-Raphson
!     technique. The full Jacobian 'fj' is constructed in 'fuljac', and
!     the net change in family concentration is calculated in 'linslv'.
!     The first few iterations (currently 7) are controlled to prevent
!     very rapid initial changes leading to divergence. Convergence is
!     determined by checking that the total concentration error 'errxo'
!     is less than tolerance 'raferr' or that the total rate error 'errpl'
!     is less than tolerance 'rafpml'. If convergence is not achieved
!     after 'nrsteps', or if divergence is encountered, the routine resets
!     the family concentrations to their initial values, and exits with a
!     non-zero value for ndxraf.
!
!     Global variables
!     ----------------
!     rafpml - tolerance - set to  1.0E-10 in input file
!     rafmin - limit for first few iterations, 0.1 in input file
!     rafmax - limit for first few iterations, 1.0E+04 in input file
!     raferr - tolerance (again?) - set to 1.0E-06 in input file
!     rafbig - maximum concentration above which divergence is assumed
!     rafeps - small non-zero concentration
!
!     Local variables
!     ---------------
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     ndxraf        Convergence exit code
!     damp1         Damping factor to apply to the first iteration
!     deltt         Reciprocal time step length  (1./cdt)
!
!  Changes for whole-atmosphere chemistry:
!  a. Increase rafeps to sqrt(peps)
!  b. Disactivate crash if too many negatives occur. They are now
!     fine during the iteration.
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE asad_spimpmjp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_SPIMPMJP_MOD'

CONTAINS

SUBROUTINE asad_spimpmjp(ndxraf, ix, jy, nlev, n_points, location, solver_iter)


USE asad_mod,              ONLY: ptol, peps, cdt,                 &
                                 f, fdot, nrsteps, nitnr, nstst,  &
                                 y, prod, slos, fj, lsvjac,       &
                                 nonzero_map, ltrig
USE asad_sparse_vars,      ONLY: spfj, spfuljac, splinslv2,       &
                                 spresolv2
USE ukca_option_mod,       ONLY: jpctr, l_ukca_quasinewton,       &
                                 i_ukca_quasinewton_start,        &
                                 i_ukca_quasinewton_end
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: &
    printstatus, PrStatus_Oper, PrStatus_Diag, umMessage, umPrint
USE UM_ParCore, ONLY: mype

USE asad_diffun_mod, ONLY: asad_diffun
USE asad_fuljac_mod, ONLY: asad_fuljac
USE asad_steady_mod, ONLY: asad_steady
USE asad_ftoy_mod,   ONLY: asad_ftoy
IMPLICIT NONE


! Subroutine interface
INTEGER, INTENT(IN) :: n_points
INTEGER, INTENT(IN) :: ix                    ! i counter
INTEGER, INTENT(IN) :: jy                    ! j counter
INTEGER, INTENT(IN) :: nlev
INTEGER, INTENT(IN) :: location
INTEGER, INTENT(OUT):: ndxraf
INTEGER, INTENT(OUT):: solver_iter ! No. of iterations

! Local variables
INTEGER, PARAMETER :: maxneg=2000     ! Max No. negatives allowed
INTEGER :: jtr
INTEGER :: jit
INTEGER :: ifi
INTEGER :: i
INTEGER :: itr
INTEGER :: nl
INTEGER :: jl
INTEGER :: kr
INTEGER :: j
INTEGER :: ip

REAL :: ztmp
REAL :: rafpml
REAL :: rafmin
REAL :: rafmax
REAL :: raferr
REAL :: rafbig
REAL :: rafeps                ! low value limit
REAL :: deltt
REAL :: damp1
REAL :: errxo
REAL :: errpl

LOGICAL :: errfl80
LOGICAL :: not_first_call = .FALSE.

REAL :: zf(n_points,jpctr)
REAL :: xoo(n_points,jpctr)
REAL :: fxo(n_points,jpctr)
REAL :: zsum(jpctr)
REAL :: tmprc(n_points,jpctr)
REAL :: bx(jpctr)
!
INTEGER, PARAMETER :: ltrig_jit=51    ! Set to nrsteps if want LTRIG

! variables required for possible quasi-Newton step:
!  array for holding the previous xoo (increment in tracer
!  concentrations) value
REAL :: xold(n_points,jpctr)
!  array which holds the previous value of fxo
REAL :: fxold(n_points,jpctr)
!  temporary array used in computing the value of xoo
REAL :: fxotmp(n_points,jpctr)
!  array containing the difference of fxo-fxold
REAL :: deltaf(n_points,jpctr)
!  array of coefficients that is used to compute xoo
REAL :: coff(n_points)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_SPIMPMJP'

!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
rafpml=1.0e-10
rafmin=1.0e-01
rafmax=1.0e+04
raferr=ptol*10.0   !   care...

rafbig=1.0/SQRT(peps)
rafeps=SQRT(peps)

ndxraf = 0
errxo = 1.0
lsvjac = .FALSE.
deltt = 1.0/cdt
damp1 = 0.5

solver_iter = 0

nl = n_points
!
!  Save values of f at start of step and make linearised first guess
zf = f
WHERE (f<rafeps) f = rafeps

! Call ASAD_STEADY at start of step to initialise deriv properly
IF (nstst /= 0) THEN
  CALL asad_steady( nl )
END IF

CALL spfuljac(n_points)
DO jtr=1,jpctr
  ip = nonzero_map(jtr,jtr)
  DO jl=1,n_points
    IF (spfj(jl,ip) > 0.0) THEN
      f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
    ELSE
      f(jl,jtr) = zf(jl,jtr) + (cdt*fdot(jl,jtr))                 &
                               /(1.0-cdt*spfj(jl,ip))
    END IF
    IF (f(jl,jtr) < rafeps) f(jl,jtr) = rafeps
  END DO
END DO

!  Start Loop - ensure mixing ratios positive and non-zero on entry
DO jit=1,nrsteps
  ifi = 0
  IF (jit == 1) ifi=nitnr
  !
  CALL asad_ftoy(not_first_call,ifi, n_points, ix, jy, nlev)
  IF (nstst /= 0) THEN
    IF (ifi == 0) CALL asad_steady( nl )
  END IF
  !
  IF (ltrig .AND. printstatus >= prstatus_oper) THEN
    DO jl=1,n_points
      WRITE(umMessage,*) 'Point: ',jl
      CALL umPrint(umMessage,src='asad_spimpmjp')
      IF (jit == 1) THEN
        WRITE(umMessage,"(1x,i2,20(1x,1pG12.4))")                          &
                    jit-1,(zf(jl,jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,i2,20(1x,1pG12.4))") jit-1,                   &
                    ((zf(jl,jtr)+cdt*fdot(jl,jtr)),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
      WRITE(umMessage,"(1x,i2,20(1x,1pG12.4))")                           &
                      jit-1,(f(jl,jtr),jtr=1,jpctr),              &
                      (y(jl,i),i=1,2)
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END DO
  END IF
  !
  CALL asad_diffun( nl )
  !
  !  Temporary prod+loss array
  tmprc(1:n_points,:) = prod(1:n_points,1:jpctr) +                &
                        slos(1:n_points,1:jpctr)
  !
  IF (errxo < raferr) THEN
    IF (jit >= ltrig_jit) THEN               ! 51, Set to 50 if want LTRIG  !!
      ndxraf = 4
      IF (printstatus >= prstatus_oper) THEN
        WRITE(umMessage,                                                   &
"('Convergence problems (',i3,1x,'iter) at location=',i3,' pe=',i3)")  &
        jit, location, mype
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
    END IF
   ! Pass back jit-1, since technically converged on previous iteration
   solver_iter = jit - 1
    GO TO 9999
  END IF

  IF (jit == nrsteps) THEN
    IF (printstatus >= prstatus_oper) THEN
      WRITE(umMessage,                    &
    "('Convergence not achieved in spimpmjp (iter',i3,') location=',"//&
    "i3,'  pe=',i3,'; halving step')") jit, location, mype
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END IF
    ndxraf = 2
    f = zf
    IF (jit >= ltrig_jit) THEN               ! 51, Set to 50 if want LTRIG  !!
      ndxraf = 4
      IF (printstatus >= prstatus_oper) THEN  
        WRITE(umMessage,                  &
      "('Convergence problems (',i3,1x,'iter) at location=',i3,'  pe=',i3)") &
      jit, location, mype
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
    END IF
    GO TO 9999
  END IF
  !
  !  Calculate fxo (& save value for possible quasi-Newton step)
  fxo = (f - zf)*deltt - fdot
  fxold = fxo
  !
  !  Test for convergence
  errfl80 = .FALSE.
  errpl = 0.0
  DO jl=1,n_points
    DO jtr=1,jpctr
      IF (ABS(tmprc(jl,jtr)) > rafeps)                             &
             errpl=MAX(errpl,ABS(fxo(jl,jtr)/tmprc(jl,jtr)))
      IF (f(jl,jtr) > rafbig) THEN  
        errfl80 = .TRUE.
      END IF
    END DO
  END DO
  IF (errfl80) THEN
    ndxraf = 3
    f = zf
    IF (jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
      ndxraf = 4
      IF (printstatus >= prstatus_oper) THEN  
        WRITE(umMessage,                  &
"('Convergence problems (',i3,1x,'iter) at location=',i3,'  pe=',i4)") &
        jit, location, mype
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
    END IF
    GO TO 9999
  END IF
  IF (errpl < rafpml) THEN
    IF (jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
      ndxraf = 4
      IF (printstatus >= prstatus_oper) THEN  
        WRITE(umMessage,                  &
"('Convergence problems (',i3,1x,'iter) at location=',i3,'  pe=',i4)")&
        jit, location, mype
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END IF
    END IF
    GO TO 9999
  END IF
  !
  !  Fill in the Jacobian, or just solve if lsvjac = .true.
  IF (lsvjac) THEN
    ! Sparse resolve routine
    CALL spresolv2(fxo,xoo,n_points,rafeps)
  ELSE
    CALL spfuljac(n_points)
    !
    IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
      WRITE(umMessage,*) 'Iteration ',jit
      CALL umPrint(umMessage,src='asad_spimpmjp')
      DO jl=1,n_points
        WRITE(umMessage,*) 'Point: ',jl
        CALL umPrint(umMessage,src='asad_spimpmjp')
        fj(jl,:,:) = 0.0
        DO jtr = 1,jpctr
          DO itr = 1,jpctr
            IF (nonzero_map(jtr,itr) > 0)                                &
              fj(jl,jtr,itr) = spfj(jl,nonzero_map(jtr,itr))
          END DO
        END DO
        DO jtr=1,jpctr
          WRITE(umMessage,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
             (fj(jl,jtr,itr),itr=1,jpctr)
          CALL umPrint(umMessage,src='asad_spimpmjp')
        END DO
      END DO
    END IF

    xoo(:,:)=0.0     ! initialisation required
    CALL splinslv2(fxo,xoo,n_points,rafeps,rafbig)

    IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
      DO jl=1,n_points
        WRITE(umMessage,*) 'Point: ',jl
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))")                        &
            'fxo',(fxo(jl,jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))")                        &
            'fdt',(fdot(jl,jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))")                        &
            'del',((f(jl,jtr)-zf(jl,jtr))*deltt,jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))")                        &
            'f  ',(f(jl,jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))")                        &
            'xoo',(xoo(jl,jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END DO
      CALL asad_fuljac(n_points)
      WRITE(umMessage,*) 'Itertion: ',jit
      CALL umPrint(umMessage,src='asad_spimpmjp')
      DO jl=1,n_points
        DO jtr=1,jpctr
          zsum(jtr) = 0.0
          DO itr=1,jpctr
            zsum(jtr)=zsum(jtr)+fj(jl,jtr,itr)*xoo(jl,itr)
          END DO
          WRITE(umMessage,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
                  (fj(jl,jtr,itr)*xoo(jl,itr),itr=1,jpctr)
          CALL umPrint(umMessage,src='asad_spimpmjp')
        END DO
        WRITE(umMessage,"(1x,a3,20(1x,1pG12.4))") 'sum',                 &
          (zsum(jtr),jtr=1,jpctr)
        CALL umPrint(umMessage,src='asad_spimpmjp')
      END DO
    END IF
    !
  END IF
  !
  errxo=0.0
  ndxraf=0
  DO jtr=1,jpctr
    DO jl=1,n_points
      !  Filter for negatives
      !  Special kick for troublesome convergence
      IF (jit == 1) xoo(jl,jtr) = damp1*xoo(jl,jtr)
      !  Calculate error
      xoo(jl,jtr) = MIN(MAX(xoo(jl,jtr),-rafbig),rafbig)
      IF (ABS(xoo(jl,jtr)) > 1.0e-16)                              &
          errxo = MAX(errxo,ABS(xoo(jl,jtr)/                      &
             MAX(f(jl,jtr),rafeps)))
      !  New mixing ratios
      ztmp = f(jl,jtr) + xoo(jl,jtr)
      !  Put limit on MAXimum correction for first few (6) iterations
      IF (jit < 7)                                                &
        ztmp = MAX(rafmin*f(jl,jtr),MIN(rafmax*f(jl,jtr),ztmp))
      !  Filter negatives and zeros
      IF (ztmp == 0.0) ztmp = rafeps
      IF (ztmp < 0.0) THEN
        ztmp = rafeps
        ndxraf = ndxraf+1
      END IF
      !  Final mixing ratios
      f(jl,jtr) = ztmp
    END DO
  END DO

  ! Perform quasi-Newton (Broyden) Method to reduce number of iterations
  ! This is done on iterations 2 <= jit <= 50, and recommended on steps 2 & 3
  ! This step will not be done if errxo < raferr, i.e. the values are
  ! converged and the routine is about to exit (this is actually tested
  ! at the top of the next loop).
  IF (l_ukca_quasinewton .AND. (errxo >= raferr)) THEN 
    IF((jit >= i_ukca_quasinewton_start) .AND.                        &
       (jit <= i_ukca_quasinewton_end)) THEN
    CALL asad_ftoy(not_first_call,ifi, n_points, ix, jy, nlev)
    CALL asad_steady(n_points)
    CALL asad_diffun(n_points)
    
    xold   = xoo
    fxo    = (f - zf)*deltt - fdot  
    xoo    = 0.0
    deltaf = fxo - fxold
    
    DO jl=1,n_points
       coff(jl) = DOT_PRODUCT(fxo(jl,:),   deltaf(jl,:))             &
            /DOT_PRODUCT(deltaf(jl,:),deltaf(jl,:))
       fxotmp(jl,:) = fxo(jl,:)*(1.0 - coff(jl))
       xold(jl,:)   = coff(jl)*xold(jl,:) 
    END DO
    CALL spresolv2(fxotmp,xoo,n_points,rafeps)

    f = f + xoo
    ! remove negative values. Does not need to be done in
    ! as intelligent a way as above, as we are not exiting
    ! the routine directly after this step.
    f = ABS(f)
    END IF
 END IF ! l_ukca_quasinewton

  ! if 5 or more negatives, drop out and halve step
  IF (ndxraf > maxneg) THEN
    WRITE(umMessage,*) ndxraf, ' Negatives - exceeds maxneg'
    CALL umPrint(umMessage,src='asad_spimpmjp')
    WRITE(umMessage,                                                      &
      "(1x,'Too many negatives (>',i4,') in spimpmjp (iter',i3,"//&
      "')  lon=',i3,'  lat=',i3,'; halving step')")               &
      maxneg,jit, location, mype
    CALL umPrint(umMessage,src='asad_spimpmjp')
    ndxraf=2
    f = zf
    IF (jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
      ndxraf = 4
      WRITE(umMessage,"(1x,'Convergence problems (',i3,1x,"//            &
       "'iter) at lon=',i3,'  lat=',i3)")                         &
        jit, location, mype
      CALL umPrint(umMessage,src='asad_spimpmjp')
    END IF
    GO TO 9999
  END IF

  ndxraf=0
END DO
!
!  Failure to Converge - reset f's and exit
WRITE(umMessage,                                                          &
  "('Convergence not achieved in spimpmjp (iter',i3,')  location=',"// &
  "i3,'  pe=',i3,'; halving step')") jit, location, mype
CALL umPrint(umMessage,src='asad_spimpmjp')
ndxraf = 2
f = zf
IF (jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
  ndxraf = 4
  WRITE(umMessage,                                                    &
"('Convergence problems (',i3,1x,'iter) at location=',i3,'  pe=',i3)")&
       jit, location, mype
  CALL umPrint(umMessage,src='asad_spimpmjp')
END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_spimpmjp
END MODULE asad_spimpmjp_mod
