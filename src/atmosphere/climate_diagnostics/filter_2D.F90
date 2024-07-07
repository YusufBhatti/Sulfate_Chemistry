! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This routine does the filtering for n<n1 and n>n2
! of the field "vort". The subroutine does:
! - Gathers the field and distributes a number of rows to each processor.
! - A fourier transform is performed using fft_track.
! - Using Spherical harmonic projection (Swarztrauber and Spotz 2003):
!      + fm_T(theta)=Fm x fm(theta) where fm is the FT coefficients over
!        longitude (fm_T is the truncated one)
!      + Fm is the projection matrix where Fm=UxU' where U is the basis of
!        the singular value decomposition of the
!      + Legendre matrix Pm=UST. Fm matrix is computed on
!        the first "track" timestep and then saved
! - Projected field fm_T(theta) is transformed back to
!   gridpoint space using fft_track
! - Gathers the resulting filtered field in
!   each processor and scatters it for STASH.
!
! Arguments:
!----
! vort: Vorticity, input
! vort_T: vorticity filtered
! n1: Lower truncation
! n2: Upper truncation
! nlim: Dimension of the model on N (longitudinal gridpoints / 2)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics

MODULE filter_2D_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FILTER_2D_MOD'

CONTAINS

SUBROUTINE filter_2D(n1, n2, nlim, delta_phi, vort, tc_activate, vort_T)

! Use array dimensions
USE atm_fields_bounds_mod,  ONLY:                                       &
     udims, vdims, vdims_s, array_dims
! Use global dimensions
USE nlsizes_namelist_mod,   ONLY: global_row_length, global_rows

! UM settings for swap-bounds
USE UM_ParVars,             ONLY: gc_all_proc_group
USE UM_ParParams,           ONLY: halo_type_no_halo
USE Field_Types,            ONLY: fld_type_v
! UM settings to print out and error handling
USE umPrintMgr,             ONLY: umPrint, umMessage
USE errormessagelength_mod, ONLY: errormessagelength
! FV-TRACK parameters
USE track_mod,              ONLY: l_hoskins, sm, Ymn, Fm, Fm_TC
! Call routine to compute Associate Legendre Polynomials
USE legendre_poly_comp_mod, ONLY: legendre_poly_comp
! Call routine to compute the projection matrix
USE proj_matrix_comp_mod,   ONLY: proj_matrix_comp
! Call routine FFT
USE fft_track_mod,          ONLY: fft_track

! DrHook parameters
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!Argument variables
INTEGER, INTENT(IN) :: n1          ! Lower truncation wavenumber
INTEGER, INTENT(IN) :: n2          ! Upper truncation wavenumber
INTEGER, INTENT(IN) :: nlim        ! Maximun wavenumber of the grid
REAL,    INTENT(IN) :: delta_phi   ! grid latitude  spacing in radians

REAL,    INTENT(IN) :: vort (udims%i_start:udims%i_end,                 &
                             vdims%j_start:vdims%j_end)
                                   ! Vorticity field

LOGICAL, INTENT(IN) :: tc_activate ! Activates TC-tracking (uses
                                   ! different settings for filtering)

REAL, INTENT(INOUT) :: vort_T(udims%i_start:udims%i_end,                &
                              vdims%j_start:vdims%j_end)
                                   ! Vorticity field truncated

   

! Local variables
INTEGER :: m_in,i,j ! indexing for DO loops

REAL ::  vort_global(global_row_length,global_rows+1)
   ! Vort Global field (all processors)

COMPLEX ::                                                              &
    vort_spectra(global_rows+1,0:nlim)                                  &
  ! Fourier coefficients of latitudinal rows ( global_rows +1 x nlim)

,   vort_spectra_T(global_rows+1,0:nlim)
  ! Truncated Fourier coefficients

LOGICAL, SAVE :: declare = .TRUE.    ! to save Fm array
LOGICAL, SAVE :: declare_TC = .TRUE. ! to save Fm_TC

LOGICAL :: l_direction
! if .TRUE. then do a direct Fourier transformation in fft_track()
! if .FALSE. then do an inverse Fourier transformation

LOGICAL, PARAMETER :: hoskins_off = .FALSE.   ! hoskins is a filter.

CHARACTER(LEN=*), PARAMETER  :: RoutineName='FILTER_2D'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate and compute Fm matrices
IF (.NOT. ALLOCATED(Fm)) ALLOCATE(Fm(global_rows+1,global_rows+1,0:n2))


IF (declare) THEN

  ! Print out the status
  WRITE(umMessage,'(A)')' Setting up Fm matrix for FV-TRACK at vor. 850hPa'
  CALL umPrint(umMessage,src='FILTER_2D')

  ! Compute Associate Legendre Polynomials
  CALL legendre_poly_comp(global_rows,delta_phi)

  ! Compute Fm matrix
  CALL proj_matrix_comp(global_rows,n1,n2,delta_phi,sm,l_hoskins,Fm)
  ! switch off the definition of arrays
  declare=.FALSE.
END IF

! For TC different truncation so another Fm matrix is needed!
IF (tc_activate) THEN

  ! Allocate Fm matrix for TC projection (to ntop_tc)
  IF (.NOT. ALLOCATED(Fm_TC))                                           &
     ALLOCATE(Fm_TC(global_rows+1,global_rows+1,0:n2))

  IF (declare_TC) THEN
    ! ++ Compute Fm matrix
    CALL proj_matrix_comp(global_rows,n1,n2,delta_phi,0.0,hoskins_off,Fm_TC)
    ! Print out the status
    WRITE(umMessage,'(A)')' Setting up Fm matrix for FV-TRACK at TC levs'
    CALL umPrint(umMessage,src='FILTER_2D')

    ! Get rid of the Associate Legendre Polynomials
    IF (ALLOCATED(Ymn)) DEALLOCATE(Ymn)
    ! switch off the definition of arrays
    declare_TC=.FALSE.
  END IF
END IF

! Gather field and split the global field in nproc latitude bands

! Initialize array for global vorticity
DO j=1,global_rows+1
  DO i=1,global_row_length
    vort_global(i,j) = 0.0
  END DO
END DO

! DEPENDS ON: gather_field
CALL gather_field( vort, vort_global,                                    &
                   udims%i_len, vdims%j_len,                             &
                   global_row_length, global_rows+1,                     &
                   fld_type_v, halo_type_no_halo,                        &
                   0, gc_all_proc_group )

! Fourier Transform over the latitudinal rows
l_direction = .TRUE.
CALL fft_track(nlim, l_direction, vort_global, vort_spectra)

! Do Legendre Projection using Fm for each m
! Initialize the truncated
DO m_in=0,nlim
  DO j=1,global_rows+1
    vort_spectra_T(j,m_in)=CMPLX(0.0,0.0)
  END DO
END DO

! For TC
IF (tc_activate) THEN

  DO m_in=0,n2
    DO j=1,global_rows+1
      ! Do projection Fm x am
      ! sum up over all vort_spectra
      DO i=1,global_rows+1
        vort_spectra_T(j,m_in)= vort_spectra_T(j,m_in) +                &
                                Fm_TC(i,j,m_in) * vort_spectra(i,m_in)
      END DO
    END DO
  END DO

ELSE
  ! For 850hPa tracking at ntop_850
  DO m_in=0,n2
    DO j=1,global_rows+1
      ! Do projection Fm x am
      ! sum up over all vort_spectra
      DO i=1,global_rows+1
        vort_spectra_T(j,m_in)= vort_spectra_T(j,m_in) +                &
                                Fm(i,j,m_in) * vort_spectra(i,m_in)
      END DO
    END DO
  END DO
  ! End loop over Legendre projection
END IF

! Inverse Fourier Transform taking vort_spectra back to vort_global
l_direction = .FALSE.
CALL fft_track(nlim, l_direction, vort_global, vort_spectra_T)
! Scatter the global field to all longitude fields
! DEPENDS ON: scatter_field
CALL scatter_field(vort_T,vort_global,                                  &
                   udims%i_len, vdims%j_len,                            &
                   global_row_length, global_rows+1,                    &
                   fld_type_v,halo_type_no_halo,                        &
                   0,gc_all_proc_group)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE filter_2D

END MODULE filter_2D_mod
