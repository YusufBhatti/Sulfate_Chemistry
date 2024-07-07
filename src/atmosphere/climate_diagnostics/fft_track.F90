! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This program computes the FT truncated to wavenumber n_max (as all the other
! n>n_max will be set to zero when truncated)
!
! Input:
! vort_global: Vorticity field (nlat x nlon)
! vort_spectra: Vorticity spectra (nlat x k_max)
! nlim: Truncation wavenumber
! Direction:
!   - If TRUE: It does the direct Fourier Transformation.
!   - If FALSE: It does the inverse transformation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics

MODULE fft_track_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FFT_TRACK_MOD'

CONTAINS

SUBROUTINE fft_track(nlim, L_direction, vort_global, vort_spectra)


USE conversions_mod,      ONLY: pi ! Get pi value
! Get global dimensions
USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows

! Get maximum n (for TC)
USE track_mod,            ONLY:                                           &
    ntop_850, ntop_tc, cos_ki, sin_ki, sin_ki_T, cos_ki_T

! DrHook modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

! Argument variables
INTEGER, INTENT(IN) ::    nlim        ! Truncation wavenumber
LOGICAL, INTENT(IN) ::    L_direction ! Do Direct or inverse FFT
REAL, INTENT(INOUT) ::    vort_global(global_row_length,global_rows+1)
   ! Vort Global field (all processors)
COMPLEX, INTENT(INOUT) :: vort_spectra(global_rows+1,0:nlim)
   ! Fourier coefficients of latitudinal rows ( global_rows x nlim )

!Local variables
INTEGER :: i,j,k ! Indexing for Do loops

INTEGER, SAVE :: n_max
   ! Maximun wavenumber to compute the cos_ki and sin_ki

REAL :: a_k      ! FT coeff for the real part
REAL :: b_k      ! FT coeff for the imaginary part
REAL :: trig_arg ! Argument for sin and cos

COMPLEX :: vort_spectra_T(0:nlim,global_rows+1)
   ! Transpose of Fourier coefficients of latitudinal rows

LOGICAL, SAVE :: declare_ki = .TRUE. ! to save cos and sin arrays

! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FFT_TRACK'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Create cosine and sine arrays for different k and i
IF (declare_ki) THEN

  n_max=MAX(ntop_850,ntop_tc)

  ! Allocate sin for each k and i
  IF (.NOT. ALLOCATED(cos_ki)) THEN
    ALLOCATE(cos_ki(global_row_length,0:n_max))
  END IF

  ! Allocate cos for each k and i
  IF (.NOT. ALLOCATED(sin_ki)) THEN
    ALLOCATE(sin_ki(global_row_length,0:n_max))
  END IF

  ! Compute the cos and sin functions
  DO k=0,n_max
    ! Set the arguments for the cosine and sine here for each k
    trig_arg = 2.0*pi*k/global_row_length
    DO i=1,global_row_length
      cos_ki(i,k) = COS( trig_arg*(i-1) )
      sin_ki(i,k) = SIN( trig_arg*(i-1) )
    END DO
  END DO

  ! Create the transpose for cos_ki and sin_ki (for optimization)
  ! cosine
  IF (.NOT. ALLOCATED(cos_ki_T)) THEN
    ALLOCATE(cos_ki_T(0:n_max,global_row_length))
  END IF

  !sine
  IF (.NOT. ALLOCATED(sin_ki_T)) THEN
    ALLOCATE(sin_ki_T(0:n_max,global_row_length))
  END IF

  ! Do transpose for inverse (for optimization)
  cos_ki_T=TRANSPOSE(cos_ki)
  sin_ki_T=TRANSPOSE(sin_ki)

  ! switch off the definition of arrays so it is not read again
  declare_ki=.FALSE.
END IF ! end if declare_ki

! ============= Direct Transformation  ===================
IF (L_direction) THEN
  ! Loop for latitude arrays
  DO k=0,n_max
    DO j = 1,global_rows+1
      ! Intialize a_k and b_k (Fourier coefficients)
      a_k = 0.0
      b_k = 0.0

      ! Using Cooley-Tukey algorithm (split between odd and even
      ! i integers in the summatory).
      ! Note: longitude definition is shifted +1; thus vort(1) -> 0 deg.
      ! and vort(global_row_length) -> 360 - delta_lambda

      DO i = 1,nlim
        ! Real coefficient (even + odd part)
        a_k = a_k + ( vort_global(2*i-1,j) * cos_ki(2*i-1,k)              &
                    + vort_global(2*i,j) * cos_ki(2*i,k) )
        ! Complex coefficient
        b_k = b_k - ( vort_global(2*i-1,j) * sin_ki(2*i-1,k)              &
                    + vort_global(2*i,j) * sin_ki(2*i,k) )
      END DO

       ! Combine in the complex number and normalize it by 1/N
      vort_spectra(j,k)=CMPLX(a_k,b_k)*1.0/REAL(global_row_length)
    END DO
  END DO

END IF ! end if for direct transformation


! ============= Inverse Transformation  ===================
IF (.NOT. L_direction) THEN

  ! Transpose vort_spectra to optimize the loop
  vort_spectra_T=TRANSPOSE(vort_spectra)

  ! Loop for latitude arrays
  DO j=1,global_rows+1
    DO i=1,global_row_length
      !Initialize the vort_global field with the zonal value
      vort_global(i,j)=REAL(vort_spectra_T(0,j))

      ! Summatory for all wavelengths k
      DO k=1,n_max
        ! Compute the sum of FT coefficients
        ! the 2 factor is for the terms where k> nlim as
        ! Real -> vort_spectra(k,j)= vort_spectra(N-k,j)
        ! Img -> vort_spectra(k,j)= -vort_spectra(N-k,j)
        ! Where N=global_row_length
        vort_global(i,j) = vort_global(i,j)                               &
                         + 2*REAL(vort_spectra_T(k,j))*cos_ki_T(k,i)      &
                         - 2* AIMAG(vort_spectra_T(k,j))*sin_ki_T(k,i)
      END DO
    END DO
  END DO

END IF ! end if for inverse transformation

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE fft_track

END MODULE fft_track_mod
