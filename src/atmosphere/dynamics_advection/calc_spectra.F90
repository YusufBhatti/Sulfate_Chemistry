! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine calc_spectra
SUBROUTINE calc_spectra(local_field, local_norm,                  &
                local_row_length, local_rows,                     &
                global_row_length, global_rows,                   &
                fld_type,halo_type,                               &
                gath_proc)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
! Description:
!   Calculates the square of the norm of the complex Fourier
!   transform.
!
! Method:
!   N.B. Subroutine FFT_2D is currently a dummy routine which returns
!   input field. This will be replaced with a suitable FFT routine.
!
!   Gathers the field and uses routine FFT_2D to obtain the
!   Fourier transforms. The required square of the norm is then
!   scattered over all processors.  Post processing is required
!   to obtain the power spectra
!   i.e If |F(n)|^2 is the square of the norm (where n is the
!   frequency) then the discrete spectral intensity (or energy),
!   E(n), is
!   For odd number of data points
!   E(n)=2.|F(n)|^2  (n=1,Nyquist frequency)
!   For even number of data points
!   E(n)=2.|F(n)|^2  (n=1,Nyquist frequency-1)
!   E(n)=F(n)|^2      n=Nyquist frequency
!
!   Documentation available.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
!  Dates should have leading zeroes in dd/mm/yy format
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables :

INTEGER ::                                                        &
 local_rows                                                       &
,local_row_length                                                 &
,global_rows                                                      &
,global_row_length                                                &
,fld_type                                                         &
,halo_type                                                        &
,gath_proc

REAL ::                                                           &
 local_field(local_row_length*local_rows)                         &
                                          !IN - original 2D field
,local_norm(local_row_length*local_rows)  !OUT - square of norm

! Local variables
INTEGER ::                                                        &
 icount                                                           &
,MINVAL                                                           &
,ifail                                                            &
,k2                                                               &
,klev                                                             &
,i,j,k

PARAMETER(MINVAL=1e-12)

REAL ::                                                           &
 spectra(global_row_length,global_rows)                           &
,spectra_im(global_row_length,global_rows)                        &
,spectra2(global_row_length*global_rows)                          &
,spectra_im2(global_row_length*global_rows)                       &
,trign(2*global_row_length)                                       &
                               !  HOLDS TRIGONOMETRIC TERMS
,trigm(2*global_rows)                                             &
                               !  USED IN FFT'S
,work(2*global_row_length*global_rows)                            &
,spec_energy(global_row_length/2+1,global_rows/2+1)               &
,sq_of_norm(global_row_length*global_rows)                        &
,field(global_row_length,global_rows)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_SPECTRA'


!---------------------------------------------------------------------
! Section 1.  Gather 2D field and prepare field for C06FUE
!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! DEPENDS ON: gather_field
CALL gather_field(local_field,spectra,                            &
                local_row_length,local_rows,                      &
                global_row_length,global_rows,                    &
                fld_type,halo_type,                               &
                gath_proc,gc_all_proc_group )

IF (mype  ==  gath_proc) THEN

  icount=1
  DO k2 = 1,global_rows
    DO k=1,global_row_length
      field(k,k2)=spectra(k,k2)
      spectra_im(k,k2)=0.0
      IF (ABS(spectra(k,k2))  <   MINVAL) THEN
        spectra(k,k2) = 0.0
      END IF
      spectra2(icount)=spectra(k,k2)
      spectra_im2(icount)=spectra_im(k,k2)
      icount=icount+1
    END DO
  END DO

  !---------------------------------------------------------------------
  ! Section 2. Call routine to calculate 2D Fourier transform
  !---------------------------------------------------------------------
  ifail=0
  ! DEPENDS ON: fft_2d
  CALL fft_2d(global_rows,global_row_length,spectra2              &
              ,spectra_im2,'I',trigm,trign,work,ifail)

  DO k =1,global_row_length*global_rows
    !          sq_of_norm(k)=spectra2(k)**2.0+
    !     &               spectra_im2(k)**2.0
    sq_of_norm(k)=spectra2(k)
  END DO

END IF ! on processor gath_pe

!---------------------------------------------------------------------
! Section 3. Scatter the square of the norm
!---------------------------------------------------------------------
! DEPENDS ON: scatter_field
CALL scatter_field(local_norm,sq_of_norm,                         &
                local_row_length,local_rows,                      &
                global_row_length,global_rows,                    &
                fld_type,halo_type,                               &
                gath_proc,gc_all_proc_group)

!!    END OF ROUTINE CALC_SPECTRA
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_spectra
