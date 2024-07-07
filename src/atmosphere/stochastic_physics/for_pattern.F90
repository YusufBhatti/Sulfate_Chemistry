! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  This routine computes the SPT Forcing Pattern with
!               a Gaussian power law, follows similar methodology as
!               "stph_skeb2", also creates a 2nd forcing pattern
!               flipping the field.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics
MODULE for_pattern_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FOR_PATTERN_MOD'

CONTAINS

SUBROUTINE for_pattern(                                                 &
! in
        row_length, rows, model_levels, first_atmstep_call)


USE timestep_mod,      ONLY:  timestep_number,timestep

! Stochastic Physics settings passed in via NameList READ
USE stochastic_physics_run_mod, ONLY:                                   &
    tau_spt, l_spt, stph_n2, spt_top_tap_lev, spt_bot_tap_lev,          &
    stphseed, stph_spt_data_check, stph_spt_data_present

! Type needed for gather field of FFT run on multiple procs
USE mpl, ONLY:                                                          &
    mpl_real

! SPT Forcing pattern field psif and psif2
USE fp_mod  , ONLY:                                                     &
    psif, psif2

! Get planet constants like planet radius or g from a module
USE planet_constants_mod, ONLY: planet_radius
! Get constant Pi from a module
USE conversions_mod, ONLY: pi

! Model level-height modules
USE level_heights_mod,     ONLY:                                        &
    r_theta_levels, eta_theta_levels

! Variables related to MPP
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

USE atm_fields_bounds_mod, ONLY:                                        &
    pdims,tdims

USE stph_closeinput_mod,  ONLY: stph_closeinput
USE stph_closeoutput_mod, ONLY: stph_closeoutput
USE stph_openinput_mod,   ONLY: stph_openinput
USE stph_openoutput_mod,  ONLY: stph_openoutput
USE stph_readentry_mod,   ONLY: stph_readentry
USE stph_writeentry_mod,  ONLY: stph_writeentry

USE pattern_spectrum2_mod, ONLY: pattern_spectrum2
USE update_pattern_mod,    ONLY: update_pattern
USE SH_spect2grid_mod,     ONLY: SH_spect2grid

USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,             ONLY: mype, nproc, nproc_max
USE UM_ParParams,           ONLY: halo_type_no_halo
USE Field_Types,            ONLY: fld_type_p

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE nlstcall_mod, ONLY: ldump
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE dump_headers_mod, ONLY: fdc_spt_coeffc_start,   &
                            fdc_spt_coeffs_start,   &
                            a_flddepc


IMPLICIT NONE

INTEGER ::                                                              &
    row_length                                                          &
                  ! local number of points on a row
,   rows                                                                &
                  ! local number of rows for u
,   model_levels
                  ! model levels

INTEGER,SAVE ::   global_row_length_div_2
                  ! global number of rows on a theta field

! ----------------------------------------------------------------
!     NLIM is the spectral truncation. its possible values are
!     constrained by the fft routine which requires it to have
!     no prime factors > 19 and the total number of primes
!     factors (including repetitions) must not exceed 20.
!     NLIM, NLAT are given values based on  model dimensions
! ----------------------------------------------------------------
INTEGER, SAVE ::                                                        &
    nlim                                                                &
             ! Spectral equivalent global row_length
,   nlat
             ! Spectral equivalent number of latitude pairs * 2

! Allocatable variables (dims depending on row_length, rows and levels)
REAL, ALLOCATABLE, SAVE ::                                              &
    psi(:,:)                                                            &
             ! Global single-level version of the pattern psif
,   psi2(:,:)
             ! Global single_level forcing pattern rotated 180


REAL,  ALLOCATABLE ::                                                   &
    my_coeff(:,:)                                                       &
             ! Local PE pattern values in grid-space
,   my_coeffc(:,:)                                                      &
             ! Local PE Fourier Space COS coeffs
,   my_coeffs(:,:)                                                      &
             ! Local PE Fourier Space SIN coeffs
,   my_coeffr(:,:)                                                      &
             ! Local PE Fourier Space Radius
,   my_phi_spt(:,:)                                                     &
             ! Local PE Fourier Space degree
,   my_phishft_spt(:,:)
             ! Local PE Fourier Space rotation

REAL, ALLOCATABLE ::                                                    &
    g_rand_nums(:,:,:)
             ! Global version of rand_nums, used to generate global
             ! random field which is SCATTERED to each PE to ensure
             ! reproducibility across different PE configurations.


REAL, ALLOCATABLE, SAVE ::                                              &
    coeffc(:,:)                                                         &
             ! spherical harmonic cosine coefficients
,   coeffs(:,:)
             ! spherical harmonic sine coefficients

REAL, ALLOCATABLE, SAVE ::                                              &
    pspect(:)
             ! wavenumber-dependent, noise amplitude
REAL,  ALLOCATABLE ::                                                   &
    coeff(:)
             ! d(psi)/d(t) used for Markov process integration
INTEGER, ALLOCATABLE ::                                                 &
    iranseed(:)
              ! Random seed size used for positioning read statement
              ! of coeffc/s from the seed file (unit=stphseed_unit)
REAL  ::                                                                &
            ! forcing pattern
    kr
            ! Ratio in z direction

REAL, SAVE ::  alpha_stoch_tend
            ! autoregressive process parameter related to tau_spt

INTEGER ::                                                              &
    i                                                                   &
            ! loop index over x direction
,   ilat                                                                &
            ! loop index over latitude
,   j                                                                   &
            ! loop index over y direction
,   k                                                                   &
            ! loop index over z direction
,   m                                                                   &
            ! loop index over EW wavespace
,   n                                                                   &
            ! loop index over NS wavespace
,   ii                                                                  &
            ! loop index ii
,   icode = 0                                                           &
            ! Return code for error reporting
,   info  = 0                                                           &
            ! return code from GC stuff
,   sndcount                                                            &
            ! Size of real array bcast to all proc's using gc_rbcast
,   tam                                                                 &
            ! scalar holding size of the random seed
,   mu
            ! latitude band
            
INTEGER, PARAMETER ::                                                   &
    zero  = 0
            ! used for identifying zero'th PE
LOGICAL ::                                                              &
    first_atmstep_call                                                  &
            ! Is true for first step of: NRUN and each CRUN
,   l_ensure_min_in_cloud_qcf
           ! Reduce CFF when qcf low to ensure minimum qcf/CFF.
CHARACTER(LEN=errormessagelength) ::                                   &
    cmessage     ! OUT error message

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FOR_PATTERN'

!    Local variables for Latitude decomposition in FFT CALL (looping
!     over latitude)
INTEGER :: fft_rows(0:nproc_max), my_rows
INTEGER :: rowcounts(0:nproc_max)
INTEGER :: displs(0:nproc_max), my_displs

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Settings required only at the first time-step
IF (first_atmstep_call) THEN
  ! Initialize variables from UM data for the sph.harm calculations
  ! For the UM, not being an spectral model, NLIM:
  IF (MOD(global_row_length,2) == 0) THEN
    nlim = global_row_length/2
  ELSE
    nlim = (global_row_length+1)/2
  END IF
  
  IF (stph_n2 >= nlim) THEN
    WRITE(umMessage,'("**ERROR**: SPT stph_n2 >= MODEL RES")')
    CALL umPrint(umMessage,src='for_pattern')
    WRITE(umMessage,'("  Section 35: check namelist settings")')
    CALL umPrint(umMessage,src='for_pattern')
    icode = 1
    WRITE (cmessage,'(A,I0,A,I0,A)') 'for_pattern stph_n2 (',stph_n2,   &
                     ') >= N (', nlim, ')'
 
    CALL ereport(routinename, icode, cmessage)
  END IF

  
  ! nlat should be equal to global_rows (and even)
  ! Uses SCATTER_FIELD for psi => psif
  IF (MOD(global_rows,2) == 0) THEN
    nlat = global_rows
    IF (.NOT. ALLOCATED(psi)) THEN
      ALLOCATE (psi(1:2*nlim,nlat))
      DO j = 1, nlat
        DO i = 1, 2*nlim
          psi(i,j) = 0.0
        END DO
      END DO
    END IF
  ELSE
    nlat = global_rows-1
    IF (.NOT. ALLOCATED(psi)) THEN
      ALLOCATE (psi(1:2*nlim,nlat+1))
      DO j = 1, nlat+1
        DO i = 1, 2*nlim
          psi(i,j) = 0.0
        END DO
      END DO
    END IF
  END IF

  ! for psi2 in case SPT is called
  IF (l_spt) THEN
    ! nlat should be equal to global_rows (and even)
    ! Uses SCATTER_FIELD for psi2 => psif2
    IF (MOD(global_rows,2) == 0) THEN
      nlat = global_rows
      IF (.NOT. ALLOCATED(psi2)) THEN
        ALLOCATE (psi2(1:2*nlim,nlat))
        DO j = 1, nlat
          DO i = 1, 2*nlim
            psi2(i,j) = 0.0
          END DO
        END DO
      END IF
    ELSE
      nlat = global_rows-1
      IF (.NOT. ALLOCATED(psi2)) THEN
        ALLOCATE (psi2(1:2*nlim,nlat+1))
        DO j = 1, nlat+1
          DO i = 1, 2*nlim
            psi2(i,j) = 0.0
          END DO
        END DO
      END IF
    END IF
  END IF

  ! The values of these variables are saved so I only need
  ! to calculate/initialise them in the first time-step
  IF (.NOT. ALLOCATED(pspect)) THEN
    ALLOCATE (pspect(nlim))
  END IF
  IF (.NOT. ALLOCATED(coeffc)) THEN
    ALLOCATE (coeffc(0:stph_n2,1:stph_n2))
  END IF
  IF (.NOT. ALLOCATED(coeffs)) THEN
    ALLOCATE (coeffs(0:stph_n2,1:stph_n2))
  END IF
  !  These arrays have sections that may not be properly initialised
  DO n = 1, stph_n2
    DO m = 0, stph_n2
      coeffc(m, n) = 0.0
      coeffs(m, n) = 0.0
    END DO
  END DO

  !Do division of global_row_length for 2nd FP
  IF (l_spt) THEN
    global_row_length_div_2=global_row_length/2
  ELSE
    global_row_length_div_2=0
  END IF

END IF ! first timestep
!------------------------------------------------------------------
!    ! Allocate work variables
IF (.NOT. ALLOCATED(coeff)) THEN
  ALLOCATE (coeff(0:2*nlim+1))
END IF

! Allocate psif SPT pattern for the module
IF (.NOT. ALLOCATED(psif)) THEN
  ALLOCATE (psif(pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,                             &
                 pdims%k_start:pdims%k_end))
END IF

 !Initialize psif, set to 0 everywhere
DO k = pdims%k_start,pdims%k_end
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      psif(i,j,k)   = 0.0
    END DO
  END DO
END DO

 ! IF SPT is called, then a second FP is required
IF (l_spt) THEN
   ! Allocate psif2 SPT pattern for the module
  IF (.NOT. ALLOCATED(psif2)) THEN
    ALLOCATE (psif2(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end))
  END IF
ELSE
  IF (.NOT. ALLOCATED(psif2)) THEN
    ALLOCATE (psif2(1,1))
  END IF
END IF

! -----------------------------------------------------------------

 ! Perform work on PE=0 to avoid extra calls to random_number which
! affects bit reproducibility by changing the random seed according
! to the number of PEs
IF (first_atmstep_call) THEN

  alpha_stoch_tend = 1.0 - EXP(-timestep/tau_spt)

  ! New power Law (Jul2010) Gaussian exp(beta x n x (n+1))
  CALL pattern_spectrum2( pspect, nlim, alpha_stoch_tend)


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Copy SPT values from dump header:
  ! ---------------------------------
  ! This section of code is called the first time the UM executes. Depending on 
  ! the timestep value, we could be starting an NRUN or CRUN. If the timestep 
  ! value is 1, then we are in an NRUN. If the timestep value is greater than 1,
  ! then we are in a CRUN. There are 3 cases where we want to retrieve the 
  ! stochastic physics values from the dump header. The first is if this is an
  ! NRUN and the stphseed value is set to 2. The second is all
  ! CRUNs. The third is if you are running an NRUN as if it were a
  ! CRUN (using the l_nrun_as_crun logical). If any of those conditions is met, 
  ! then we also check that the data is actually in the dump. If these checks
  ! are passed, then we proceed with copying the SPT values from the additional 
  ! parameters section of the dump header to the relevant location here.
  
  IF (((stphseed == 2 .OR. l_nrun_as_crun) .AND.                        &
        timestep_number == 1 .AND.                                      &
        stph_spt_data_check == stph_spt_data_present) .OR.              &
        (timestep_number > 1 .AND.                                      &
        stph_spt_data_check == stph_spt_data_present) ) THEN
    
    IF (l_nrun_as_crun) THEN
        cmessage = 'l_nrun_as_crun is TRUE:' //  &
           'Reading random numbers from the dump (equivalent of stphseed=2)'
        icode = -1
        CALL ereport(RoutineName, icode, cmessage)
    END IF

    ! copy from dump header if available
    IF (mype == 0) THEN
      CALL umPrint('retrieving coeffc/coeffs from dump header array' , &
                   src='for_pattern')

      DO n = 1, stph_n2
        DO m = 0, stph_n2
          i = m + (n-1) * (stph_n2+1)
          coeffc(m, n) = a_flddepc (fdc_spt_coeffc_start + i) 
          coeffs(m, n) = a_flddepc (fdc_spt_coeffs_start + i)
        END DO 
      END DO
    END IF
    ! Scatter restart coeffs to all processors
    sndcount = (stph_n2+1)*(stph_n2-1+1)
    CALL gc_rbcast(3247, sndcount, zero, nproc, icode, coeffc)
    CALL gc_rbcast(3248, sndcount, zero, nproc, icode, coeffs)
    

  END IF       ! reading in coeffc and coeffs
      
END IF         ! first timestep

CALL update_pattern( alpha_stoch_tend, coeffc, coeffs, nlim,            &
                     pspect, zero, icode, info)
                     
! Compute decomposition of latitude loop
! Array "displs" keeps track of the real latitudes offset on each PE
displs(0) = 0
DO i = 0, nproc-1
  fft_rows(i) = nlat/nproc
  IF (i >= nproc - MOD(nlat,nproc)) THEN
    fft_rows(i) = fft_rows(i) + 1
  END IF
  rowcounts(i) = fft_rows(i)
  IF (i > 0) THEN
    displs(i) = SUM(rowcounts(0:i-1))
  END IF
END DO
! Number of rows on this PE
my_rows = fft_rows(mype)
! Number of rows before this on lower PEs
my_displs = displs(mype)
! Allocate local data arrays
! Real space, Fourier Sin, Cos Coeffs
ALLOCATE (my_coeff(1:2*nlim, my_rows))
ALLOCATE (my_coeffc(0:stph_n2,1:stph_n2))
ALLOCATE (my_coeffs(0:stph_n2,1:stph_n2))
ALLOCATE (my_phi_spt(0:stph_n2,1:stph_n2))
ALLOCATE (my_coeffr(0:stph_n2,1:stph_n2))
ALLOCATE (my_phishft_spt(0:stph_n2,1:stph_n2))


!  These arrays have sections that may not be properly initialised
DO n = 1, stph_n2
  DO m = n+1, stph_n2
    my_coeffc(m, n) = 0.0
    my_coeffs(m, n) = 0.0
  END DO
END DO
DO n = 1, stph_n2-1
  DO m = 0, stph_n2
    my_coeffr(m, n) = 0.0
    my_phishft_spt(m, n) = 0.0
    my_phi_spt(m, n) = 0.0
  END DO
END DO

DO n = 1, stph_n2
  DO m = 0, n
    my_coeffr(m, n) = SQRT(coeffc(m, n)**2 + coeffs(m, n)**2)
    !    Determine angle from SIN & COS wave components (single step)
    my_phi_spt(m, n) = ATAN2(coeffs(m, n), coeffc(m, n))
    !    Max shift ranges from 0 <-> pi  for wavenos 1 <-> stph_n2
    my_phishft_spt(m, n) = (stph_n2 - MAX(m, n)) * pi/ (stph_n2-1)
  END DO
END DO

DO k = spt_bot_tap_lev, spt_top_tap_lev
  ! Adjust coefficients in the vertical
  ! Level 1 = no change -> 12km Level = max change (=pi)
  !  cycles around above that level
  kr =  eta_theta_levels(k) * (r_theta_levels(1,1,tdims%k_end) -        &
        planet_radius)/ 12000.0


  DO n = 1, stph_n2
    DO m = 0, n

      my_coeffc(m, n) = my_coeffr(m, n) * COS(my_phi_spt(m, n) +        &
                         kr * my_phishft_spt(m, n))
      my_coeffs(m, n) = my_coeffr(m, n) * SIN(my_phi_spt(m, n) +        &
                         kr * my_phishft_spt(m, n))

    END DO
  END DO

  ii = 1
  DO ilat = 1, my_rows
    mu = ilat + my_displs
    ! Calculates coeff in grid space

    ! Calculate coeff in grid space
    CALL SH_spect2grid( my_coeffc, my_coeffs, stph_n2, nlim, coeff,     &
                       mu, nlat, first_atmstep_call, ii)

    ! Copy coeff to 2-D array
    my_coeff(:,ilat) = coeff(0:2*nlim-1)
  END DO

  ! Gather field split by latitude and then distribute as normal
  CALL gc_get_communicator(icode, info)
  CALL mpl_gatherv( my_coeff(1,1), my_rows*2*nlim, mpl_real,            &
                  psi, rowcounts*2*nlim, displs*2*nlim, mpl_real,       &
                  zero, icode, info)
  ! Scatter global field psi to local array psif on each processor
  cmessage=''


  !Scatter psif to each processor
  ! DEPENDS ON: scatter_field
  CALL scatter_field( psif(pdims%i_start, pdims%j_start,k), psi(1,1),   &
                      row_length, rows, global_row_length,              &
                      global_rows, fld_type_p, halo_type_no_halo,       &
                      zero, gc_all_proc_group)


  ! Create psif2 for SPT slow physics-> FP1 rotated 180 round the pole
  ! FP2 =FP1 (lambda + 180) where lamdba is latitude
  ! NOTE: 2D pattern, uniform in the vertical
  ! Use FP 2D field at spt_bot_tap_lev + 1 (arbitrary choice)
  IF (l_spt .AND. k == spt_bot_tap_lev + 1) THEN
    DO j = 1, global_rows
      ! Put first half into the second one!
      DO i = 1, global_row_length_div_2
        psi2(i,j)=psi(global_row_length_div_2+i,j)
      END DO
      ! Other way around!
      DO i = 1, global_row_length_div_2
        psi2(global_row_length_div_2+i,j)=psi(i,j)
      END DO
    END DO

    ! Decompose psif_global2 into psif2 for each processor
    ! DEPENDS ON: scatter_field
    CALL scatter_field(psif2(pdims%i_start, pdims%j_start),             &
                       psi2(1,1), row_length, rows,                     &
                       global_row_length, global_rows, fld_type_p,      &
                       halo_type_no_halo, zero, gc_all_proc_group)

  END IF

END DO ! End loop over the vertical levels

! Deallocate work arrays
IF (ALLOCATED(my_coeff)) DEALLOCATE(my_coeff)
IF (ALLOCATED(my_coeffc)) DEALLOCATE(my_coeffc)
IF (ALLOCATED(my_coeffs)) DEALLOCATE(my_coeffs)
IF (ALLOCATED(my_phi_spt)) DEALLOCATE(my_phi_spt)
IF (ALLOCATED(my_coeffr)) DEALLOCATE(my_coeffr)
IF (ALLOCATED(my_phishft_spt)) DEALLOCATE(my_phishft_spt)


! Append coeffc/s to seed file (if active) at dump times.
! Generally used for CRUNs.
IF (ldump) THEN
  IF (mype == 0) THEN
    CALL umPrint('Copying SPT parameters to dump headers', &
                 src='for_pattern')
    DO n = 1, stph_n2
      DO m = 0, stph_n2
        i = m + (n-1) * (stph_n2+1)
        a_flddepc (fdc_spt_coeffc_start + i) = coeffc(m, n)
        a_flddepc (fdc_spt_coeffs_start + i) = coeffs(m, n)
      END DO
    END DO
  END IF
  
END IF



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

! --------------------------------------------------------------------
END SUBROUTINE for_pattern

END MODULE for_pattern_mod
