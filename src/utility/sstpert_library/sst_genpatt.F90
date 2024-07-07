! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:  Generates random forcing pattern with given power
!               function (gspect)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: SSTpert library

MODULE sst_genpatt_mod

IMPLICIT NONE 

CONTAINS

SUBROUTINE sst_genpatt( global_row_length, global_rows, psif )

USE skeb_forcing_mod, ONLY: skeb_forcing
USE update_dpsidt_mod, ONLY: update_dpsidt
USE stochastic_physics_run_mod, ONLY: firsttimestep_true, Ymn, stph_n2
USE legendre_poly_comp_stph_mod, ONLY: legendre_poly_comp_stph
USE cderived_mod, ONLY: delta_phi

IMPLICIT NONE


INTEGER, INTENT (IN) ::                                                &
    global_row_length                                                  &
                  ! global number of points on a row
,   global_rows
                  ! global number of rows on a theta field

REAL, INTENT (INOUT) ::                                                &
    psif(global_row_length, global_rows)

! ----------------------------------------------------------------
!     NLIM is the spectral truncation. its possible values are
!     constrained by the fft routine which requires it to have
!     no prime factors > 19 and the total number of primes
!     factors (including repetitions) must not exceed 20.

!     NLIM, NLAT are given values based on  model dimensions
! ----------------------------------------------------------------

INTEGER ::                                                             &
     nlim                                                              &
             ! Spectral equivalent global row_length
,    nlat                                                              &
             ! Spectral equivalent number of latitude pairs * 2
,    icode, info
             ! Error return codes
INTEGER, PARAMETER ::                                                  &
     n1 = 1, n2 = 50, zero = 0
             ! Wavenumber range of power-spectrum and PE0 number
             ! Note that further down the routine there is a section 
             ! which sets (hardcodes) the array "gspect" assuming these 
             ! sizes and hence if they are modified that section must
             ! also be updated accordingly

! Allocatable variables (dims depending on row_length and rows)

REAL, ALLOCATABLE ::                                                   &
    dpsidtc(:,:)                                                       &
             ! 2D version of d(psi)/d(t) COS coeffs in Fourier
,   dpsidts(:,:)                                                       &
             ! 2D version of d(psi)/d(t) SIN coeffs in Fourier
,   seqalf  (:,:,:)                                                     
             ! The Legendre Polynomials read by skeb_forcing

REAL, ALLOCATABLE ::                                                   &
    dpsidt(:)                                                          &
             ! d(psi)/d(t) used for Markov process integration
,   gspect(:)
             ! Wave-number dependent noise amplitude

INTEGER, ALLOCATABLE ::                                                &
    iranseed(:)
             ! Random seed size used for positioning read statement
             ! of dpsidtc/s from the seed file (unit=149)

REAL, PARAMETER ::                                                     &
    alpha = 0.5
            ! dummy argument passed to shared routine {0 < alpha < 1}

INTEGER ::                                                             &
    i                                                                  &
            ! loop index over x direction
,   ilat                                                               &
            ! loop index over latitude
,   j                                                                  &
            ! loop index over y direction
,   m                                                                  &
            ! loop index over EW wavespace
,   n                                                                  &
            ! loop index over NS wavespace
,   m1                                                                 &
            ! Loop index used to populate seqalf
,   nn                                                             
            ! Loop index used to populate seqalf

! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------

! Initialize variables from UM data for the sph.harm calculations
! For the UM, not being an spectral model, NLIM:

IF (MOD(global_row_length,2) == 0) THEN
  nlim=global_row_length/2
ELSE
  nlim=(global_row_length+1)/2
END IF

! nlat should be equal to global_rows (and even)
! Uses SCATTER_FIELD for psi => psif
IF (MOD(global_rows,2) == 0) THEN
  nlat = global_rows
ELSE
  nlat = global_rows-1
END IF

! Allocate work variables
ALLOCATE (dpsidt(0:2*nlim+1))
ALLOCATE (gspect(nlim))
ALLOCATE (dpsidtc(0:n2,n1:n2))
ALLOCATE (dpsidts(0:n2,n1:n2))

DO n = n1, n2
  DO m = 0, n2
    dpsidtc(m, n) = 0.0
    dpsidts(m, n) = 0.0
  END DO
END DO

! Power spectrum derived from analysis of SST data - note that n1 and n2 are
! parameters defined at the top of this routine; but the hardcoding of this
! array means they can't be changed without this being changed too
gspect(n1:n2) = (/0.378777,0.704892,1.026114,1.203914,1.255679,        &
                  1.359430,1.478994,1.429783,1.526795,1.500924,        &
                  1.490598,1.465734,1.477826,1.484674,1.448614,        &
                  1.421848,1.410815,1.438613,1.386196,1.348778,        &
                  1.308668,1.316984,1.251467,1.214483,1.196684,        &
                  1.160380,1.132038,1.067881,1.057473,1.035978,        &
                  1.006957,0.979293,0.958451,0.936141,0.914674,        &
                  0.896876,0.869185,0.862772,0.839728,0.821467,        &
                  0.802636,0.768016,0.745082,0.720489,0.697255,        &
                  0.668994,0.647636,0.634320,0.615679,0.599321/)
gspect(n2+1:nlim) = 0.0

CALL update_dpsidt(alpha, dpsidtc, dpsidts, n1, n2, nlim,              &
                   gspect, zero, icode, info, firsttimestep_true)

stph_n2 = n2
CALL legendre_poly_comp_stph(global_rows, delta_phi, .TRUE.)

! Populate the array with the sequence of Legendre Polynomials read by
! the skeb_forcing method.
  
ALLOCATE(seqalf(n1:n2,0:nlim,nlat))
DO ilat = 1, nlat
  DO m1 = 0, 2*nlim, 2
    m  = m1/2
    nn = MAX(n1,m)
    DO n = nn, n2
      seqalf(n,m,ilat)= Ymn(ilat,n,m)
    END DO
  END DO    
END DO

DO ilat = 1, nlat
   
  ! Calculates dpsidt in grid space
  CALL skeb_forcing( dpsidtc, dpsidts, n1, n2, nlim, dpsidt,           &
                       ilat, nlat, seqalf)

  ! Copy dpsidt to 2-D array
  psif(:,ilat)= dpsidt(0:2*nlim-1)
END DO

DEALLOCATE (seqalf)
DEALLOCATE (dpsidt)
DEALLOCATE (dpsidts)
DEALLOCATE (dpsidtc)
DEALLOCATE (gspect)

END SUBROUTINE sst_genpatt

END MODULE sst_genpatt_mod
