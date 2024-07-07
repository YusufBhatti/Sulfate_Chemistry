! *****************************COPYRIGHT*******************************
! (c) [University of California] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
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
!       copyright noticed, this list of conditions and the following
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
!   Fast-jx is an updated routine for calculating online photolysis rates
!   This module contains routines that read in the phase factors etc from file
!   Based upon the fast_specs routine, though with large differences caused by
!   differences between fast-j and fast-jx
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE fastjx_specs

USE    fastjx_data
USE umPrintMgr
USE ereport_mod, ONLY: ereport
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE um_parcore,  ONLY: mype
USE setup_namelist, ONLY: setup_nml_type
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! ********************************************
! STRINGS

! String containing description of data set
CHARACTER(LEN=78) :: title0

! String containing cloud/aerosol scattering
CHARACTER(LEN=7)  ::  titlaa(a_)

! String containing species being photolysed
CHARACTER(LEN=7)  ::  titlej(x_)

! Dummy strings containing duplicate species strings
CHARACTER(LEN=7)  ::  titlej2,titlej3

! *********************************************

PUBLIC fastjx_rd_js, fastjx_rd_mie, fastjx_rd_xxx, fastjx_rd_sol

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FASTJX_SPECS'

CONTAINS

! ######################################################################
      ! Reads in quantum yields and labels
      ! Copied from fastj_rd_js routine
SUBROUTINE fastjx_rd_js

USE ukca_chem_defs_mod,   ONLY: ratj_t, ratj_defs
IMPLICIT NONE
!     Read in quantum yield 'jfacta' and fastjx label 'jlabel'

!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use

INTEGER :: i

CHARACTER(LEN=10) :: adjusted_fname
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_JS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
adjusted_fname = ' '

DO i=1,jppj
  jfacta(i)=ratj_defs(i)%jfacta/100.0e0
  adjusted_fname=TRIM(ADJUSTL(ratj_defs(i)%fname))
  jlabel(i)=adjusted_fname(1:7)
  IF (PrintStatus >= PrStatus_Diag) THEN
    WRITE(umMessage,'(I6,E12.3,A12)') i,jfacta(i),jlabel(i)
    CALL umPrint(umMessage,src='fastjx_specs')
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_rd_js


! ######################################################################
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
SUBROUTINE fastjx_rd_mie(nj1,namfil)

IMPLICIT NONE

INTEGER, INTENT(IN) :: nj1  ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN) :: namfil ! Name of scattering data file
                                       ! (e.g., FJX_scat.dat)
CHARACTER (LEN=errormessagelength) :: cmessage        
                                      ! Contains string for error handling
INTEGER :: errcode                    ! error code

INTEGER :: i, j, k                      ! Loop variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_MIE'

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = 2 + 2*a_ + (3+sw_phases)*sw_band_aer*a_
INTEGER, PARAMETER :: n_chars = 78 + 7*a_

TYPE my_namelist
  SEQUENCE
  INTEGER :: naa
  INTEGER :: jtaumx
  REAL :: atau
  REAL :: atau0
  REAL :: raa(a_)
  REAL :: daa(a_)
  REAL :: waa(sw_band_aer,a_)
  REAL :: qaa(sw_band_aer,a_)
  REAL :: saa(sw_band_aer,a_)
  REAL :: paa(sw_phases, sw_band_aer,a_)
  CHARACTER (LEN=78 ) :: title0
  CHARACTER (LEN=7) :: titlaa(a_)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! ***********************************
! End of Header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

!remove optional arguments that correspond to n_type=0
CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_int_in=n_int, n_real_in=n_real,               &
                    n_chars_in=n_chars)

IF (mype == 0) THEN
  ! Open data file containing aerosl/cloud data
  OPEN (nj1,FILE=namfil,STATUS='OLD',FORM='FORMATTED',ACTION='READ')

  ! Read number of data types and title
  READ (nj1,'(i2,a78)') naa,title0

  ! If the number of data types exceeds maximum allowed exit with an error
  IF (naa > a_) THEN
    cmessage = 'Too many scattering data sets'
    errcode = 100
    CALL ereport('FASTJX_RD_MIE',errcode,cmessage)
  END IF

  ! Read Cloud layering variables
  READ (nj1,'(5x,I5,2F10.5)') jtaumx,atau,atau0
  IF (printstatus >= prstatus_oper) THEN  
    WRITE(umMessage,'(a,2F9.5,I5)')         &
        ' atau/atau0/jmx',atau,atau0,jtaumx
    CALL umPrint(umMessage,src='fastjx_specs')
  END IF

  ! Read blank line
  READ (nj1,*)

  ! Loop over aerosol types
  DO j = 1,naa

    ! Read title, effective radius and density
    READ (nj1,'(3x,a20,32x,f5.3,15x,f5.3)')  titlaa(j),raa(j),daa(j)

    ! Loop over 5 wavelength bins
    DO k = 1,5
      ! read wavelength, q, scattering albedo and phases
      READ (nj1,'(f4.0,f7.4,f7.4,7f6.3,1x,f7.3,f8.4)') &
        waa(k,j),qaa(k,j),saa(k,j),(paa(i,k,j),i=2,8)
      ! set first phase to 0
      paa(1,k,j) = 1.0e0

    END DO ! wavelengths
  END DO ! aerosols

  ! Close file
  CLOSE(nj1)

  my_nml % naa    = naa
  my_nml % jtaumx = jtaumx
  my_nml % atau   = atau
  my_nml % atau0  = atau0
  my_nml % raa    = raa
  my_nml % daa    = daa
  my_nml % waa    = waa
  my_nml % qaa    = qaa
  my_nml % saa    = saa
  my_nml % paa    = paa
  my_nml % title0 = title0
  my_nml % titlaa = titlaa

END IF  ! mype == 0

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF ( mype /= 0 ) THEN

  naa    = my_nml % naa
  jtaumx = my_nml % jtaumx
  atau   = my_nml % atau
  atau0  = my_nml % atau0
  raa    = my_nml % raa
  daa    = my_nml % daa
  waa    = my_nml % waa
  qaa    = my_nml % qaa
  saa    = my_nml % saa
  paa    = my_nml % paa
  title0 = my_nml % title0
  titlaa = my_nml % titlaa

END IF

IF (printstatus >= prstatus_oper) THEN
  ! Output some information
  WRITE(umMessage,'(A,5F8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'  &
         ,(waa(k,1),k=1,5)
  CALL umPrint(umMessage,src='fastjx_specs')

  ! Output file title
  WRITE(umMessage,*) title0
  CALL umPrint(umMessage,src='fastjx_specs')

  ! Loop over aerosol types writing radius, density and q
  DO j=1,naa
    WRITE(umMessage,'(i3,1x,a8,7f8.3)') j,titlaa(j),raa(j),daa(j),       &
        (qaa(k,j),k=1,5)
    CALL umPrint(umMessage,src='fastjx_specs')
  END DO
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_rd_mie

!########################################################################
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 (parm_CTM.f) should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (parm_MIE.f).
!           if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!-----------------------------------------------------------------------
SUBROUTINE fastjx_rd_xxx(nj1,namfil)

IMPLICIT NONE

INTEGER, INTENT(IN)       :: nj1 ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN)  :: namfil    ! Name of spectral data file
                                           ! (JX_spec.dat)
CHARACTER (LEN=errormessagelength)        :: cmessage      
                                           ! String for error handling
INTEGER                   :: errcode       ! errror code

INTEGER ::  i, j, jj, k, iw                ! Loop variables
INTEGER ::  nqqq                           ! No of cross sections read in
INTEGER ::  nwww                           ! No of wavelength bins
INTEGER ::  nqrd                           ! number of x-sections
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FASTJX_RD_XXX'

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 6
INTEGER, PARAMETER :: n_real = 3*wx_ +1 + 3*x_ + 9*wx_ +2*wx_*x_
INTEGER, PARAMETER :: n_chars = 78 + 7*x_

TYPE my_namelist
  SEQUENCE
  INTEGER :: njval
  INTEGER :: nqrd
  INTEGER :: nwww
  INTEGER :: nw1
  INTEGER :: nw2
  INTEGER :: nqqq
  REAL :: wl(wx_)
  REAL :: fl(wx_)
  REAL :: qrayl(wx_+1)
  REAL :: tqq(3,x_)
  REAL :: qo2(wx_,3)
  REAL :: qo3(wx_,3)
  REAL :: q1d(wx_,3)
  REAL :: qqq(wx_,2, x_)
  CHARACTER (LEN=78 ) :: title0
  CHARACTER (LEN=7) :: titlej(x_)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! *****************************

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

!remove optional arguments that correspond to n_type=0
CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_int_in=n_int, n_real_in=n_real,               &
                    n_chars_in=n_chars)

! Initialise the temperatures to 0
DO j = 1,x_
  DO k = 1,3
    tqq(k,j) = 0.0e0
  END DO
END DO

!----------spectral data----set for new format data J-ver8.3------------------
!         note that NJVAL = # J-values, but NQQQ (>NJVAL) = # Xsects read in
!         for 2005a data, NJVAL = 62 (including a spare XXXX) and
!              NQQQ = 64 so that 4 wavelength datasets read in for acetone
!         note NQQQ is not used outside this subroutine!
! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

IF (mype == 0) THEN
  ! Open file containing cross sections
  OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted',ACTION='READ')

  ! Read title of file
  READ (nj1,'(a)') title0

  ! Read number of photolysed species, number of x-sections & number of wavelength bins
  READ (nj1,'(10x,5i5)') njval, nqrd, nwww

  ! set maximum and minimum wavelngth bins
  nw1 = 1
  nw2 = nwww

  ! Check that number of photolysed species and number of cross sections
  ! doesn't exceed maximum allowed. If either do then exit
  IF (njval > x_ .OR. nqrd > x_) THEN
    cmessage = 'Number of Cross Sections exceeds Maximum Allowed'
    errcode = 100
    CALL ereport('FASTJX_RD_XXX',errcode,cmessage)
  END IF

  !----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
  READ (nj1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        (wl(iw),iw=1,nwww)
  READ (nj1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        (fl(iw),iw=1,nwww)
  READ (nj1,'(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        (qrayl(iw),iw=1,nwww)

  !---READ O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej(1),tqq(1,1), (qo2(iw,1),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej2,  tqq(2,1), (qo2(iw,2),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej3,  tqq(3,1), (qo2(iw,3),iw=1,nwww)

  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej(2),tqq(1,2), (qo3(iw,1),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej2,  tqq(2,2), (qo3(iw,2),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej3,  tqq(3,2), (qo3(iw,3),iw=1,nwww)

  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej(3),tqq(1,3), (q1d(iw,1),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej2,  tqq(2,3), (q1d(iw,2),iw=1,nwww)
  READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej3,  tqq(3,3), (q1d(iw,3),iw=1,nwww)

  ! WRITE information to stdout
  IF (printstatus >= prstatus_oper) THEN
    DO j = 1,3
      WRITE(umMessage,'(i6,a7,3e10.3)') j,titlej(j),(tqq(i,j),i=1,3)
      CALL umPrint(umMessage,src='fastjx_specs')
    END DO
  END IF

  !---READ remaining species:  X-sections at 2 T_s
  jj = 4
  DO j = 4,nqrd

    READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej(jj),tqq(1,jj),(qqq(iw,1,jj),iw=1,nwww)
    READ (nj1,'(a7,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))')           &
        titlej2,  tqq(2,jj),(qqq(iw,2,jj),iw=1,nwww)

    !---include stratospheric J's (this also includes Cl and Br compounds!)
    IF (w_ == 18 .OR. titlej2(7:7) /= 'x') THEN
      IF (printstatus >= prstatus_oper) THEN
        WRITE(umMessage,'(i6,a7,2e10.3)') jj,titlej(jj), (tqq(i,jj),i=1,2)
        CALL umPrint(umMessage,src='fastjx_specs')
      END IF
      jj = jj+1
    END IF

  END DO
  nqqq = jj-1
  njval = njval + (nqqq - nqrd)

  my_nml % njval  = njval
  my_nml % nqrd   = nqrd
  my_nml % nwww   = nwww
  my_nml % nw1    = nw1
  my_nml % nw2    = nw2
  my_nml % nqqq   = nqqq
  my_nml % wl     = wl
  my_nml % fl     = fl
  my_nml % qrayl  = qrayl
  my_nml % tqq    = tqq
  my_nml % qo2    = qo2
  my_nml % qo3    = qo3
  my_nml % q1d    = q1d
  my_nml % qqq    = qqq
  my_nml % title0 = title0
  my_nml % titlej = titlej

END IF  ! mype == 0

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF ( mype /= 0 ) THEN

  njval  = my_nml % njval
  nqrd   = my_nml % nqrd
  nwww   = my_nml % nwww
  nw1    = my_nml % nw1
  nw2    = my_nml % nw2
  nqqq   = my_nml % nqqq
  wl     = my_nml % wl
  fl     = my_nml % fl
  qrayl  = my_nml % qrayl
  tqq    = my_nml % tqq
  qo2    = my_nml % qo2
  qo3    = my_nml % qo3
  q1d    = my_nml % q1d
  qqq    = my_nml % qqq
  title0 = my_nml % title0
  titlej = my_nml % titlej

END IF

CALL mpl_type_free(mpl_nml_type,icode)

!---truncate number of wavelengths to DO troposphere-only
IF (w_ /= wx_) THEN

  !---TROP-ONLY
  IF (w_ == 12) THEN
    IF (printstatus >= prstatus_oper) THEN  
      CALL umPrint(' >>>TROP-ONLY reduce wavelengths to 12,'//           &
          ' drop strat X-sects',                                         &
          src='fastjx_specs')
    END IF
    nw2 = 12

    ! Remove first four wavelength bins from  total
    DO iw = 1,4
      wl(iw) = wl(iw+4)
      fl(iw) = fl(iw+4)
      qrayl(iw) = qrayl(iw+4)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+4,k)
        qo3(iw,k) = qo3(iw+4,k)
        q1d(iw,k) = q1d(iw+4,k)
      END DO

      DO j = 4,nqqq
        qqq(iw,1,j) = qqq(iw+4,1,j)
        qqq(iw,2,j) = qqq(iw+4,2,j)
      END DO
    END DO

    ! Remove 9/10 wavelength bins from total
    DO iw = 5,12
      wl(iw) = wl(iw+6)
      fl(iw) = fl(iw+6)
      qrayl(iw) = qrayl(iw+6)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+6,k)
        qo3(iw,k) = qo3(iw+6,k)
        q1d(iw,k) = q1d(iw+6,k)
      END DO
      DO j = 4,nqqq
        qqq(iw,1,j) = qqq(iw+6,1,j)
        qqq(iw,2,j) = qqq(iw+6,2,j)
      END DO
    END DO

    !---TROP-QUICK  (must scale solar flux for W=5)
  ELSE IF (w_ == 8) THEN
    IF (printstatus >= prstatus_oper) THEN  
      CALL umPrint(' >>>TROP-QUICK reduce wavelengths to 8, '//          &
          'drop strat X-sects', &
          src='fastjx_specs')
    END IF
    nw2 = 8

    DO iw = 1,1
      wl(iw) = wl(iw+4)
      fl(iw) = fl(iw+4)*2.0e0
      qrayl(iw) = qrayl(iw+4)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+4,k)
        qo3(iw,k) = qo3(iw+4,k)
        q1d(iw,k) = q1d(iw+4,k)
      END DO

      DO j = 4,nqqq
        qqq(iw,1,j) = qqq(iw+4,1,j)
        qqq(iw,2,j) = qqq(iw+4,2,j)
      END DO
    END DO

    DO iw = 2,8
      wl(iw) = wl(iw+10)
      fl(iw) = fl(iw+10)
      qrayl(iw) = qrayl(iw+10)

      DO k = 1,3
        qo2(iw,k) = qo2(iw+10,k)
        qo3(iw,k) = qo3(iw+10,k)
        q1d(iw,k) = q1d(iw+10,k)
      END DO

      DO j = 4,nqqq
        qqq(iw,1,j) = qqq(iw+10,1,j)
        qqq(iw,2,j) = qqq(iw+10,2,j)
      END DO
    END DO

  ELSE
    cmessage = 'Incorrect Number of Wavelength Bins, must be 8, 12 or 18'
    errcode = 100
    CALL ereport('FASTJX_RD_XXX',errcode,cmessage)
  END IF
END IF

! Close file
IF ( mype == 0 ) CLOSE(nj1)

! *************************************************************
! Map local indices to UKCA ones

DO j = 1,njval
  DO k = 1,jppj
    IF ( k == 1 .AND. printstatus >= prstatus_oper) THEN
      WRITE(umMessage,*) 'FASTJX Compare titles ', titlej(j), jlabel(k)
      CALL umPrint(umMessage,src='fastjx_specs')
    END IF
    IF (jlabel(k) == titlej(j)) jind(k) = j
  END DO
END DO

DO k = 1,jppj

  IF (printstatus >= prstatus_oper) THEN
    WRITE(umMessage,*) 'Comparing J-rate for photolysis reaction ',      &
        jlabel(k),' ?'
    CALL umPrint(umMessage,src='fastjx_specs')
    WRITE(umMessage,*) 'Using index ', jind(k), MAXVAL(qqq(:,:,jind(k)))
    CALL umPrint(umMessage,src='fastjx_specs')
  END IF

  IF (jfacta(k) < 1.0e-20 .AND. printstatus >= prstatus_oper) THEN
    WRITE(umMessage,*) 'Not using photolysis reaction ',jlabel(k)
    CALL umPrint(umMessage,src='fastjx_specs')
  END IF

  IF (jind(k) == 0) THEN

    IF (jfacta(k) < 1.0e-20) THEN
      jind(k)=1
    ELSE
      cmessage =                                                  &
       ' Which J-rate for photolysis reaction '//jlabel(k)//' ?'
      errcode = 1
      CALL ereport("FASTJX_RD_XXX",errcode,cmessage)
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fastjx_rd_xxx

!########################################################################
!-----------------------------------------------------------------------
!  Read in solar cycle data
!-----------------------------------------------------------------------
SUBROUTINE fastjx_rd_sol(nj1,namfil)

IMPLICIT NONE

INTEGER, INTENT(IN)       :: nj1 ! Channel number for reading data file
CHARACTER(LEN=*), INTENT(IN)  :: namfil    ! Name of solar cycle data file
CHARACTER (LEN=errormessagelength)        :: cmessage      
                                           ! String for error handling
INTEGER                   :: errcode       ! errror code

INTEGER ::  i, j                           ! Loop variables
INTEGER ::  n                              ! Array length

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_real = 18 + 203 + 1476 + 128

TYPE my_namelist
  SEQUENCE
  REAL :: solcyc_spec(18)
  REAL :: solcyc_quanta(203)
  REAL :: solcyc_ts(1476)
  REAL :: solcyc_av(128)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! *****************************

CALL gc_get_communicator(my_comm, icode)

!remove optional arguments that correspond to n_type=0
CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_real_in=n_real)


IF (mype == 0) THEN
  OPEN (nj1,FILE=namfil,STATUS='old',FORM='formatted',ACTION='READ')

  READ (nj1,'(50x,i2)') n 
  DO i = 1,CEILING(REAL(n)/6.)
    READ (nj1,'(6e10.3)')                             &
            (solcyc_spec(j+(i-1)*6),j=1,6)
  END DO 

  READ (nj1,'(50x,i3)') n
  DO i = 1,CEILING(REAL(n)/7.) 
      READ (nj1,'(7e10.3)')                           &
            (solcyc_quanta(j+(i-1)*7),j=1,7)
  END DO

  READ (nj1,'(50x,i4)') n
  DO i = 1,CEILING(REAL(n)/6.)
      READ (nj1,'(6f9.5)')                             &
            (solcyc_ts(j+(i-1)*6),j=1,6)
  END DO

  READ (nj1,'(50x,i3)') n
  DO i = 1,CEILING(REAL(n)/8.)
    READ (nj1,'(8f9.5)')                             &
            (solcyc_av(j+(i-1)*8),j=1,8)
  END DO

  CLOSE(nj1)
  
  my_nml % solcyc_spec   = solcyc_spec
  my_nml % solcyc_quanta = solcyc_quanta
  my_nml % solcyc_ts     = solcyc_ts
  my_nml % solcyc_av     = solcyc_av

END IF  ! mype == 0

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF ( mype /= 0 ) THEN

  solcyc_spec    = my_nml % solcyc_spec
  solcyc_quanta  = my_nml % solcyc_quanta
  solcyc_ts      = my_nml % solcyc_ts
  solcyc_av      = my_nml % solcyc_av

END IF

CALL mpl_type_free(mpl_nml_type,icode)

END SUBROUTINE fastjx_rd_sol


!#######################################################################

END MODULE fastjx_specs
