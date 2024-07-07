! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Stochastic Kinetic Energy Backscatter V2 (SKEB2)
!     This code provides a function routine for computing horizontally
!     non-divergent wind increments due to a supposed upscale energy
!     transfer from the mesoscale. The forcing is computed in the
!     spherical domain by expanding a streamfunction forcing function
!     in spherical harmonics. Each harmonic mode has a time variation
!     described by a first-order Markov process with in-built memory
!     corresponding to a typical mesoscale eddy turnover time (tau_skeb
!     ~ 0.5 -> 1 day).
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics
MODULE stph_skeb2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_SKEB2_MOD'

CONTAINS

SUBROUTINE stph_skeb2(                                                  &
! in
           row_length, rows, n_rows, model_levels,                      &
           delta_phi, delta_lambda,                                     &
           rho_r2, u, v, w,                                             &
! in/out
           r_u, r_v,                                                    &
           stashwork35, first_atmstep_call)

USE swap_bounds_mv_mod, ONLY: swap_bounds_mv

! SKEB2 gui settings passed in via NameList READ
USE stochastic_physics_run_mod, ONLY:                                   &
    stph_n1, stph_n2, skeb2_toplev, skeb2_botlev, br                    &
 ,  tot_backscat, tau_skeb, alphac, l_skeb2_psicdisp                    &
 ,  l_skeb2_psisdisp, l_skeb2_skeb1disp                                 &
 ,  skeb2_sdisp, type_smag, type_bihm                                   &
 ,  skeb2_cdisp, type_cape, type_mflx                                   &
 ,  l_skeb2_conv_disp_mod, l_skeb2_velpot                               &
 ,  skeb2_up_flux, skeb2_dwn_flux, skeb2_cape                           &
 ,  sdispfac, cdispfac, kdispfac, nsmooth                               &
 ,  offx_stph, offy_stph, mask_pdamp, mask_smooth                       &
 ,  l_skebsmooth_adv, l_skebprint, stphseed           &
 ,  firsttimestep_true, firsttimestep_false, firsttimestep_crun         &
 ! SPT flag to deactivate the deallocation of mass_flux array
 ,  l_spt, Ymn, stph_skeb2_data_present, stph_skeb2_data_check

! Sructure containing stochastic physics diagnostics
USE stph_diag_mod, ONLY:                                                &
    strstphdiag
! Sructure containing field details for multi-variable SWAP_BOUNDS
USE swapable_field_mod, ONLY:                                           &
    swapable_field_pointer_type
! Type needed for gather field of FFT run on multiple procs
USE mpl, ONLY:                                                          &
    mpl_real

! Model level-height modules
USE level_heights_mod,     ONLY:                                        &
    r_rho_levels, r_theta_levels, eta_theta_levels

! Bounds of arrays
USE atm_fields_bounds_mod, ONLY:                                        &
    pdims, pdims_s,                                                     &
    udims, udims_s,                                                     &
    vdims, vdims_s,                                                     &
    tdims, wdims_s, stphdims_l

! Model grid trigonometry
USE trignometric_mod, ONLY:                                             &
    cos_theta_latitude, sin_theta_latitude,                             &
    cos_v_latitude, sin_v_latitude

! Variables related to MPP
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows

USE timestep_mod,      ONLY:                                            &
    timestep, timestep_number

! Routines for global sums
USE global_2d_sums_mod, ONLY:                                           &
    global_2d_sums

! * tau_skeb:          Decorrelation time (~5.5 hrs in this case)
! * tot_backscat: Global-mean rate of energy backscatter in m**2.s**(-3)
! * br:           Backscatter ratio of dissipation Energy backscattered
! * alphac:       Updraught proportion of gridbox (0.2%)

USE planet_constants_mod, ONLY: g, planet_radius

USE conversions_mod, ONLY: pi

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,           ONLY: mype, nproc, nproc_max
USE UM_ParParams,         ONLY: pnorth, psouth
USE Field_Types,          ONLY: fld_type_p, fld_type_u, fld_type_v
USE backscatter_spectrum_mod, ONLY: backscatter_spectrum
USE diagnostics_stph_mod, ONLY: diagnostics_stph
USE skeb_forcing_mod,     ONLY: skeb_forcing
USE biharm_diss_mod,      ONLY: biharm_diss
USE skeb_smagorinsky_mod, ONLY: skeb_smagorinsky
USE stph_closeinput_mod,  ONLY: stph_closeinput
USE stph_closeoutput_mod, ONLY: stph_closeoutput
USE stph_openinput_mod,   ONLY: stph_openinput
USE stph_openoutput_mod,  ONLY: stph_openoutput
USE stph_readentry_mod,   ONLY: stph_readentry
USE stph_skeb1_mod,       ONLY: stph_skeb1
USE stph_writeentry_mod,  ONLY: stph_writeentry
USE update_dpsidt_mod,    ONLY: update_dpsidt
USE stash_array_mod, ONLY: sf
USE nlstcall_mod, ONLY: ldump
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE dump_headers_mod, ONLY: fdc_skeb2_dpsidtc_start, &
                            fdc_skeb2_dpsidts_start, &
                            a_flddepc

USE mpp_conf_mod,  ONLY: swap_field_is_vector, swap_field_is_scalar

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

! Get tendencies for stph
USE physics_tendencies_mod,  ONLY:                                      &
    du_stph, dv_stph, l_retain_stph_tendencies

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
    row_length                                                          &
             ! local number of points on a row
,   rows                                                                &
             ! local number of rows for u
,   n_rows                                                              &
             ! local number of rows for v
,   model_levels
             ! model levels

REAL, INTENT (IN) ::                                                    &
    delta_lambda                                                        &
             ! EW (x) grid spacing in radians
,   delta_phi                                                           &
             ! NS (y) grid spacing in radians
,   rho_r2(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end,                               &
           pdims_s%k_start:pdims_s%k_end)                               &
             ! density * square of radius of earth
,   u     (udims_s%i_start:udims_s%i_end,                               &
           udims_s%j_start:udims_s%j_end,                               &
           udims_s%k_start:udims_s%k_end)                               &
             ! main prog. u field
,   v     (vdims_s%i_start:vdims_s%i_end,                               &
           vdims_s%j_start:vdims_s%j_end,                               &
           vdims_s%k_start:vdims_s%k_end)                               &
             ! main prog. v field
,   w     (wdims_s%i_start:wdims_s%i_end,                               &
           wdims_s%j_start:wdims_s%j_end,                               &
           wdims_s%k_start:wdims_s%k_end)
          ! main prog. w field
REAL, INTENT (INOUT) ::                                                 &
    r_u   (udims_s%i_start:udims_s%i_end,                               &
           udims_s%j_start:udims_s%j_end,                               &
           udims_s%k_start:udims_s%k_end)                               &
             ! main physics u increment (SKEB2 added to this)
,   r_v   (vdims_s%i_start:vdims_s%i_end,                               &
           vdims_s%j_start:vdims_s%j_end,                               &
           vdims_s%k_start:vdims_s%k_end)                               &
             ! main physics v increment (SKEB2 added to this)
,   stashwork35(*)
             ! Array containing requested sec35 STASH diagnostics

REAL   ::                                                               &
    skeb2_urot(udims%i_start:udims%i_end,                               &
               udims%j_start:udims%j_end,                               &
               udims%k_start:udims%k_end)                               &
             ! rotational u increment from SKEB2
,   skeb2_vrot(vdims%i_start:vdims%i_end,                               &
               vdims%j_start:vdims%j_end,                               &
               vdims%k_start:vdims%k_end)                               &
             ! rotational v increment from SKEB2 on V-points
,   skeb2_udiv(udims%i_start:udims%i_end,                               &
               udims%j_start:udims%j_end,                               &
               udims%k_start:udims%k_end)                               &
             ! divergent u increment from SKEB2
,   skeb2_vdiv(vdims%i_start:vdims%i_end,                               &
               vdims%j_start:vdims%j_end,                               &
               vdims%k_start:vdims%k_end)
             ! divergent v increment from SKEB2 on V-points


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
,   nlat                                                                &
             ! Spectral equivalent number of latitude pairs * 2
,   level2km
             ! model level at ~2km (=10 for L38)

REAL, SAVE ::                                                           &
    max_tot_backscat                                                    &
             ! Maximum allowed backscatter
,   crit_tot_backscat                                                   &
             ! Critical backscatter amount requiring damping
,   logscale                                                            
             ! Multiplication factor for LOG10 arguments (see levfac)

REAL ::                                                                 &
    levfac
             ! LOG10(model_level), used to reduce

! Allocatable work arrays
REAL, ALLOCATABLE ::                                                    &
    my_dpsidtc   (:,:)                                                  &
             ! Local PE Fourier Space COS coeffs
,   my_dpsidts   (:,:)                                                  &
             ! Local PE Fourier Space SIN coeffs
,   my_dpsidtr   (:,:)                                                  &
             ! Local PE Fourier Space Radius
,   my_phi       (:,:)                                                  &
             ! Local PE Fourier Space degree
,   my_phishft   (:,:)                                                  &
             ! Local PE Fourier Space rotation
,   work_z_halo  (:,:,:)                                                &
             ! Dummy variable to swap bounds and hold data on vort-grid
,   vert_int_work(:,:)                                                  &
             ! vertically integrated flux W.m^-2 (working array)
,   dpsidt       (:)                                                    &
             ! d(psi)/d(t) used for Markov process integration
,   sum_temp     (:)                                                    &
             ! Temporary array to hold row sums for rvecsumr
,   psif0        (:,:)                                                  &
             ! streamfunction forcing field with halo
,   sm_tot_disp0 (:,:,:)
             ! holding variable for smoothing of modulating field

! Allocatable GC arrays
REAL, ALLOCATABLE ::                                                    &
    sendbuf   (:,:,:,:)                                                 &
             ! Send the local PE d(psi)/d(t) to other PEs in the same
             ! logical processed grid (LPG) row
,   recvbuf   (:,:,:,:)                                                 
             ! Receive the local PE d(psi)/d(t) from other PEs in the same 
             ! logical processed grid (LPG) row

! Allocatable variables (needing to be saved between timesteps)
REAL, ALLOCATABLE, SAVE ::                                              &
    dpsidtc (:,:)                                                       &
             ! 2D version of d(psi)/d(t) COS coeffs in Fourier
,   dpsidts (:,:)                                                       &
             ! 2D version of d(psi)/d(t) SIN coeffs in Fourier
,   dx_theta(:,:)                                                       &
             ! Delta-X on u points
,   dy_theta(:,:)                                                       &
             ! Delta-Y on u points
,   dx_v    (:,:)                                                       &
             ! Delta-X on v points
,   dy_v    (:,:)                                                       &
             ! Delta-Y on v points
,   darea   (:,:)                                                       &
             ! Horizontal area of gridboxes
,   c_dadt  (:,:)                                                       &
             ! Horizontal area of gridboxes * timestep (inverted)
,   gspect  (:)                                                         &
             ! Wave-number dependent noise amplitude                    
,   seqalf  (:,:,:)
             ! The Legendre Polynomials read by skeb_forcing

INTEGER, ALLOCATABLE ::                                                 &
    iranseed(:)
             ! Random seed size used for positioning read statement
             ! of dpsidtc/s from the seed file (unit=stphseed_unit)

! Local arrays and scalars for dissipation calculations
REAL  ::                                                                &
    psif     (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! streamfunction forcing field (E=tot_backscat)
,   m_psif   (pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)                            &
             ! Modulated stream_function forcing field
             ! Also used to swap bounds and hold data on p-grid
,   rand_nums(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! 3 dimensional random field
,   alpha                                                               &
             ! autoregressive process parameter related to tau_skeb
,   delr_rho (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! thickness between r_rho_levels
,   mass     (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! mass (kg) of 3D grid boxes
,   kr                                                                  &
             ! Ratio in z direction
,   smallp = 1.0e-08                                                    &
             ! Small real number
,   local_crit_tot_backscat                                             &
             ! crit_tot_backscat at local level
,   gltotke                                                             &
             ! Used for global sum of quantity for printing
,   gltotke_tmp(1)
             ! Array in global sum

REAL, ALLOCATABLE  ::                                                   &
    mass_2d  (:,:)
             ! mass (kg) of column

! Loop counters
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
             ! loop index ii inside i-loop
,   jj                                                                  &
             ! loop index jj inside j-loop
,   i1                                                                  &
             ! loop index i1 inside i-loop
,   i2                                                                  &
             ! loop index i2 inside j-loop
,   j1                                                                  &
             ! loop index j1 inside j-loop
,   j2                                                                  &
             ! loop index j2 inside j-loop
,   ip1                                                                 &
             ! i plus one
,   im1                                                                 &
             ! i minus one
,   jp1                                                                 &
             ! j plus one
,   jm1                                                                 &
             ! j minus one
,   ju                                                                  &
             ! j on u-grid (limit to rows if at north pole)
,   jv                                                                  &
             ! j on v-grid (limit to n_rows if at north pole)
,   ismooth                                                             &
             ! loop index over smoothing iterations
,   icode = 0                                                           &
             ! Return code for error reporting
,   info  = 0                                                           &
             ! return code from GC stuff
,   ifield                                                              &
             ! number of fields used in swap_bounds_mv
,   firsttimestep                                                       &
             ! indicates whether first timestep of NRUN or CRUN
,   sndcount                                                            &
             ! Size of real array bcast to all proc's using gc_rbcast
,   tam                                                                 &
             ! scalar holding size of the random seed 
,   mu                                                                  &  
             ! mu is the latitude band
,   m1                                                                  & 
             ! loop index used to populate seqalf
,   nn                                                                  & 
             ! loop index used to populate seqalf
,   pe                                                                  &
             ! loop index between all PEs
,   pe_i                                                                &
             ! loop index between PEs on the same row
,   pe_j                                                              
             ! loop index between PEs on the same column.

REAL ::                                                                 &
    rho      (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! density on each model level
,   sm_tot_disp(pdims%i_start:pdims%i_end,                              &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! smoothed total instantaneous dissipation rate
,   diff_flux(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! difference of the up-dwn mass flux
,   sdisp    (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! dissipitation in SMAGORINSKY bit
,   cdisp    (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! dissipitation from Convection
,   kdisp    (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,                                &
              pdims%k_start:pdims%k_end)                                &
             ! dissipitation in SKEB1-type (KE-based) calculation
,   cape_kk  (pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end)
             ! local modified cape field

! added the resolution dependent factor to mimimize
! convection rate a high resolutions
REAL , SAVE :: conv_res_factor

LOGICAL ::                                                              &
    l_latwgt = .FALSE.                                                  &
             ! Set to TRUE will apply sin(lat)/cos(lat) weights to
             ! SKEB2 rot/div wind increments (default=.FALSE.)
,   first_atmstep_call
             ! Used to set firsttimestep
             ! Is true for first step of: NRUN and each CRUN

CHARACTER(LEN=errormessagelength) :: cmessage
             ! OUT error message

! Parameters
INTEGER, PARAMETER :: zero = 0 ! used for identifying zero'th PE
CHARACTER(LEN=*)       :: RoutineName
PARAMETER( RoutineName='STPH_SKEB2' )

! Local targets
REAL, TARGET ::                                                         &
    skeb2_urot0(udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end)                          &
             ! rotational u increment on U-points
,   skeb2_vrot0(vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end)                          &
             !  rotational v increment on V-points (before smoothing)
,   skeb2_udiv0(udims_s%i_start:udims_s%i_end,                          &
                udims_s%j_start:udims_s%j_end,                          &
                udims_s%k_start:udims_s%k_end)                          &
             ! divergent u increment on U-points
,   skeb2_vdiv0(vdims_s%i_start:vdims_s%i_end,                          &
                vdims_s%j_start:vdims_s%j_end,                          &
                vdims_s%k_start:vdims_s%k_end)
             ! divergent v increment on V-points (before smoothing)

! Local variables for Latitude decomposition in FFT CALL (looping over
! latitude). 
INTEGER ::                                                              &
     my_rows                                                            &
             ! The rows in the local sub-stripe of the FFT decomposition
,    my_displs                                                          &
             ! The real latitudes offset of the local sub-stripe 
,    pe_rows(0:nproc_max)                                               &
             ! The number of rows in the sub-stripe of the FFT decomposition
             ! on each PE
,    pe_displs(0:nproc_max)                                             &
             ! The real latitudes offset of the sub-stripe on each PE
,    pe_offstr(0:nproc_max)                                             &
             ! The offset of the sub-stripe inside a stripe for each PE
,    stripe_rows(0:nproc_max)                                           &
             ! The number of rows of the local stripe on each PE
,    max_row_length                                                     & 
             ! Maximum row_length between all the PEs.
,    max_my_rows                                                        &
             ! Maximum rows between all the sub-stripes (PEs).
,    dstart                                                             &
             ! row_length offset inside the stripe for a PE
,    dsize                                                              &
             ! row_lenght inside the stripe for a PE
,    kdiff    
             ! Different between top level and bottom level in skeb2. 

! GC data communication size
INTEGER ::                                                              &
    sendcount                                                           &
             ! Number of elements sent 
,   recvcount
             ! Number of elements received

! Declaration of Stochastic Physics diagnostics.
TYPE (strstphdiag)                 :: stph_diag
! Declaration of fields used in swap_bounds_mv
TYPE (swapable_field_pointer_type) :: fields_to_swap(8)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! NOTE on Arakawa C-grid
!      Smagorinsky and Convection estimates of Energy Dissipation
!      are stored on P-points. Final streamfunction forcing (S) used
!      to calculate winds needs to be interpolated to the vorticty
!      (Z) points for rot winds to avoid the U|V grids being swapped.
!      Rot winds are from PSI on Z-points: (u;v) = (-dS/dy; dS/dx)
!      Div winds are from PSI on P-points: (u;v) = ( dS/dx; dS/dy)
!
!  LON=  0   1/2dX  dX   3/2dX
!        |     :     |     :
!        |     :     |     :
!      --V-----Z-----V-----Z-- 3/2dY
!        |     :     |     :
!        |     :     |     :
!      ==P=====U=====P=====U== dY
!        |     :     |     :
!        |     :     |     :
!      --V-----Z-----V-----Z-- 1/2dY
!        |     :     |     :
!        |     :     |     :
!      ==P=====U=====P=====U== SOUTH POLE
!
! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Error trap if this is a LAM configuration (SKEB2 should not be called)
IF (model_type /= mt_global) THEN
  WRITE(umMessage,'("**ERROR**: SKEB2 not available in a Limited '      &
       // 'Area Model")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  icode   = 1
  WRITE (cmessage,*) 'STPH_SKEB2: Not configured to run in LAMS. ',     &
                     'Turn scheme off in namelist/gui.'


  CALL ereport(routinename, icode, cmessage)
END IF  ! .NOT. GLOBAL

! Check bounds of SKEB2_TOPLEV and SKEB2_BOTLEV
! SKEB2 cannot extend to top and bottom model levels because of
! calculations involving levels (k-1) and (k+1) in the code
! and physical realism of the thickness of the backscattered layer
IF (skeb2_toplev >= model_levels) THEN
  WRITE(umMessage,'("**ERROR**: SKEB2 TOP LEVEL = OR > MODEL TOP")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  icode   = 1
  WRITE (cmessage,'(A,I0,A,I0,A)') 'STPH_SKEB2: skeb2_toplev (',        &
         skeb2_toplev, ') >= model_levels (', model_levels, ')'

  CALL ereport(routinename, icode, cmessage)
END IF
IF (skeb2_botlev <= 1) THEN
  WRITE(umMessage,'("ERROR: SKEB2 BOTTOM LEVEL = OR < MODEL BASE")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  icode   = 1
  WRITE (cmessage,'(A,I0,A)') 'STPH_SKEB2: skeb2_botlev (',             &
         skeb2_botlev, ') <= 1'

  CALL ereport(routinename, icode, cmessage)
END IF
IF (skeb2_toplev < skeb2_botlev + 3) THEN
  WRITE(umMessage,'("ERROR: SKEB2 VERTICAL EXTENT TOO THIN")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("       SKEB2 TOP LEVEL < BOTTOM LEVEL+4")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("  Section 35: check namelist settings")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  icode   = 1
  WRITE (cmessage,'(A,I0,A,I0,A)')                                      &
        'STPH_SKEB2: Need at least 3 levels between skeb2_toplev (',    &
        skeb2_toplev, ') and skeb2_botlev (',skeb2_botlev,')'

  CALL ereport(routinename, icode, cmessage)
END IF


! Some calculations can be done once and stored
firsttimestep = firsttimestep_false
IF (first_atmstep_call) THEN
  firsttimestep = firsttimestep_true
  IF (timestep_number > 1 )                        &
      firsttimestep = firsttimestep_crun
END IF

IF (firsttimestep == firsttimestep_true .OR.                            &
    firsttimestep == firsttimestep_crun) THEN

  ! Initialize variables from UM data for the sph.harm calculations
  ! For the UM, not being an spectral model, NLIM:

  IF (MOD(global_row_length,2) == 0) THEN
    nlim=global_row_length/2
  ELSE
    nlim=(global_row_length+1)/2
  END IF
  
  IF (stph_n2 >= nlim) THEN
    WRITE(umMessage,'("**ERROR**: SKEB2 stph_n2 >= MODEL RES")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("  Section 35: check namelist settings")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    icode   = 1
    WRITE (cmessage,'(A,I0,A,I0,A)') 'STPH_SKEB2: stph_n2 (',stph_n2,   &
                     ') >= N (', nlim, ')'
  
    CALL ereport(routinename, icode, cmessage)
  END IF

  ! nlat should be equal to global_rows (and even)
  IF (MOD(global_rows,2) == 0) THEN
    nlat = global_rows
  ELSE
    nlat = global_rows-1
  END IF

  ! The values of these variables are saved so I only need
  ! to calculate/initialise them in the first time-step
  IF (.NOT. ALLOCATED(gspect)) THEN
    ALLOCATE (gspect(nlim))
  END IF
  IF (.NOT. ALLOCATED(dpsidtc)) THEN
    ALLOCATE (dpsidtc(0:stph_n2,stph_n1:stph_n2))
  END IF
  IF (.NOT. ALLOCATED(dpsidts)) THEN
    ALLOCATE (dpsidts(0:stph_n2,stph_n1:stph_n2))
  END IF
  
  IF (.NOT. ALLOCATED(dx_theta)) THEN
    ALLOCATE (dx_theta(udims%i_start:udims%i_end,                       &
                       udims%j_start:udims%j_end))
  END IF
  IF (.NOT. ALLOCATED(dy_theta)) THEN
    ALLOCATE (dy_theta(udims%i_start:udims%i_end,                       &
                       udims%j_start:udims%j_end))
  END IF
  IF (.NOT. ALLOCATED(dx_v)) THEN
    ALLOCATE (dx_v    (vdims%i_start:vdims%i_end,                       &
                       vdims%j_start:vdims%j_end))
  END IF
  IF (.NOT. ALLOCATED(dy_v)) THEN
    ALLOCATE (dy_v    (vdims%i_start:vdims%i_end,                       &
                       vdims%j_start:vdims%j_end))
  END IF
  IF (.NOT. ALLOCATED(darea)) THEN
    ALLOCATE (darea   (pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
  END IF
  IF (.NOT. ALLOCATED(c_dadt)) THEN
    ALLOCATE (c_dadt  (pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
  END IF

  !  These arrays have sections that may not be properly initialised
  DO n = stph_n1, stph_n2
    DO m = 0, stph_n2
      dpsidtc(m, n) = 0.0
      dpsidts(m, n) = 0.0
    END DO
  END DO

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(l_latwgt,udims,vdims,pdims,dx_theta,dy_theta,delta_lambda,&
!$OMP&    delta_phi,sin_theta_latitude,planet_radius,dx_v,dy_v,         &
!$OMP&    cos_v_latitude,darea,cos_theta_latitude,c_dadt,timestep,      &
!$OMP&    sin_v_latitude)

  !  Calculate and store grid increment values
  IF (l_latwgt) THEN
    ! Weight wind incr by sin/cos latitude to decrease impact of
    ! rot-wind in tropics and div-wind in high latitudes
    ! Note: dx = dlam * a * cos(lat), so "cos(lat)" term cancels out
!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        dx_theta(i,j) = delta_lambda * planet_radius
        dy_theta(i,j) = (delta_phi * planet_radius)/                    &
                         sin_theta_latitude(1,j)
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dx_v(i,j) = (delta_lambda * planet_radius * cos_v_latitude(i,j))&
                    / sin_v_latitude(i,j)
        dy_v(i,j) = (delta_phi * planet_radius)/ cos_v_latitude(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
  ELSE
    ! No wind incr weighting (this is the default)
!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        dx_theta(i,j) = delta_lambda * planet_radius
        dy_theta(i,j) = delta_phi * planet_radius
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dx_v(i,j) = delta_lambda * planet_radius
        dy_v(i,j) = delta_phi * planet_radius
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF   ! l_latwgt

  !  Limit ratio between dX & dY to 1:3 to eliminate generating large
  !  wind increments near the poles
  !  WHERE(dx_theta < 0.333*dy_theta) dx_theta = 0.333*dy_theta
  !  WHERE(dx_v < 0.333*dy_v) dx_v = 0.333*dy_v

  !  Horizontal area of gridbox (weighted by latitude)
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      darea(i,j) = delta_lambda * planet_radius * delta_phi *           &
                   planet_radius * cos_theta_latitude(i,j)
      !  Set minimum grid-box area (to avoid divide-by-zero at poles)
      IF (darea(i,j) < 1.0e3) darea(i,j) = 1.0e3
      !  Invert array so that it can be used as a multiplication factor
      c_dadt(i,j) = 1.0/(darea(i,j)*timestep)
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  ! Keep a record of the settings for this run in the PE output files
  IF (printstatus  >=  prstatus_normal) THEN
    WRITE(umMessage,*) 'L_SKEB2_PSISDISP = ', l_skeb2_psisdisp
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,*) 'L_SKEB2_PSICDISP = ', l_skeb2_psicdisp
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,*) 'L_SKEB2_SKEB1DISP = ', l_skeb2_skeb1disp
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,I3)') 'SKEB2_CDISP type = ', skeb2_cdisp
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,L3)') 'L_SKEB2_VELPOT = ', l_skeb2_velpot
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    IF (l_skeb2_psicdisp) THEN
      IF (skeb2_cdisp /= type_cape .AND. skeb2_cdisp /= type_mflx) THEN
        WRITE(umMessage,'("**warning**: SKEB2 convective dissipation '  &
             // 'is on but no valid scheme is selected")')
        CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
        WRITE(umMessage,'("  Section 35: check namelist settings")')
        CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
        icode    = -1
        WRITE (cmessage,*) 'STPH_SKEB2: no convective dissipation '     &
                           // 'scheme is selected - will return zeroes'
        CALL ereport(routinename, icode, cmessage)
      END IF
    END IF
  END IF

  ! Find 1st model theta level above 2km
  ! Model levels are constant in time
  DO k = wdims_s%k_start, wdims_s%k_end
    level2km = k
    IF (eta_theta_levels(k) * (r_theta_levels(1,1,wdims_s%k_end) -      &
        planet_radius) > 2000.0) EXIT
  END DO

  ! Convert range [mod_lev=1; mod_lev=level2km] => [1; 10]
  logscale = 10.0/level2km
  IF (printstatus  ==  prstatus_diag) THEN
    WRITE(umMessage,'(" ---------------------------------------- ")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(" Vertical ramp decrease of wind increment ")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(" Level of 2km = ",I4)') level2km
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(" Logscale = ",ES12.4)') logscale
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF

  ! Compute resolution dependent factor:
    ! Let Assume N216 is the right value (as used in MOGREPS)
  IF (l_skeb2_conv_disp_mod) THEN
    conv_res_factor=SQRT(216.0/nlim)
  ELSE
    conv_res_factor= 1.0
  END IF

END IF ! firsttimestep_true .OR. firsttimestep_crun

! Allocate work variables
IF (.NOT. ALLOCATED(vert_int_work)) THEN
  ALLOCATE (vert_int_work(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end))
END IF
IF (.NOT. ALLOCATED(work_z_halo)) THEN
  ALLOCATE (work_z_halo(  vdims_s%i_start:vdims_s%i_end,                &
                          vdims_s%j_start:vdims_s%j_end,                &
                          vdims_s%k_start:vdims_s%k_end))
END IF

! Prepare STASH diagnostics where requested
! Allocate STASH diagnostic arrays using Structure TYPE stph_diag
stph_diag%l_skeb2_u              = sf(1,35)
stph_diag%l_skeb2_v              = sf(2,35)
stph_diag%l_skeb2_u_incr         = sf(3,35)
stph_diag%l_skeb2_v_incr         = sf(4,35)
stph_diag%l_skeb2_u_rot          = sf(5,35)
stph_diag%l_skeb2_v_rot          = sf(6,35)
stph_diag%l_skeb2_u_div          = sf(7,35)
stph_diag%l_skeb2_v_div          = sf(8,35)
stph_diag%l_skeb2_disp_smag      = sf(9,35)
stph_diag%l_skeb2_disp_conv      = sf(10,35)
stph_diag%l_skeb2_disp_skeb1     = sf(11,35)
stph_diag%l_skeb2_smodfield      = sf(12,35)
stph_diag%l_skeb2_streamfunction = sf(13,35)
stph_diag%l_skeb2_random_pattern = sf(14,35)
stph_diag%l_skeb2_ke_psif        = sf(15,35)
stph_diag%l_skeb2_ke_sdisp       = sf(16,35)
stph_diag%l_skeb2_ke_cdisp       = sf(17,35)
stph_diag%l_skeb2_ke_kdisp       = sf(18,35)
stph_diag%l_skeb2_ke_m_psif      = sf(19,35)
stph_diag%l_skeb2_ke_prewindincr = sf(20,35)
stph_diag%l_skeb2_ke_windincr    = sf(21,35)
stph_diag%l_skeb2_ke_postwindincr= sf(22,35)

IF (sf(0,35)) THEN

  ! u after skeb2
  IF (stph_diag%L_skeb2_u) THEN
    ALLOCATE(stph_diag%skeb2_u(udims%i_start:udims%i_end,               &
                               udims%j_start:udims%j_end,               &
                               udims%k_start:udims%k_end))
  ELSE
    ! Code to allocate unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_u(1,1,1))
    stph_diag%skeb2_u(:,:,:) = 0.0
  END IF

  ! v after skeb2
  IF (stph_diag%L_skeb2_v) THEN
    ALLOCATE(stph_diag%skeb2_v(vdims%i_start:vdims%i_end,               &
                               vdims%j_start:vdims%j_end,               &
                               vdims%k_start:vdims%k_end))
  ELSE
    ! Code to allocate unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_v(1,1,1))
    stph_diag%skeb2_v(:,:,:) = 0.0
  END IF

  ! u increment diagnostic (store increment before SKEB2)
  IF (stph_diag%l_skeb2_u_incr) THEN
    ALLOCATE(stph_diag%skeb2_u_incr(udims%i_start:udims%i_end,          &
                                    udims%j_start:udims%j_end,          &
                                    udims%k_start:udims%k_end))
    DO k = udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          stph_diag%skeb2_u_incr(i,j,k) = r_u(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_u_incr(1,1,1))
    stph_diag%skeb2_u_incr(:,:,:) = 0.0
  END IF

  ! v increment diagnostic (store increment before SKEB2)
  IF (stph_diag%L_skeb2_v_incr) THEN
    ALLOCATE(stph_diag%skeb2_v_incr(vdims%i_start:vdims%i_end,          &
                                    vdims%j_start:vdims%j_end,          &
                                    vdims%k_start:vdims%k_end))
    DO k = vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          stph_diag%skeb2_v_incr(i,j,k) = r_v(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_v_incr(1,1,1))
    stph_diag%skeb2_v_incr(:,:,:) = 0.0
  END IF

  ! rotational u increments from SKEB2
  IF (stph_diag%L_skeb2_u_rot) THEN
    ALLOCATE(stph_diag%skeb2_u_rot(udims%i_start:udims%i_end,           &
                                   udims%j_start:udims%j_end,           &
                                   udims%k_start:udims%k_end))
    stph_diag%skeb2_u_rot(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_u_rot(1,1,1))
    stph_diag%skeb2_u_rot(:,:,:) = 0.0
  END IF

  ! rotational v increments from SKEB2
  IF (stph_diag%L_skeb2_v_rot) THEN
    ALLOCATE(stph_diag%skeb2_v_rot(vdims%i_start:vdims%i_end,           &
                                   vdims%j_start:vdims%j_end,           &
                                   vdims%k_start:vdims%k_end))
    stph_diag%skeb2_v_rot(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_v_rot(1,1,1))
    stph_diag%skeb2_v_rot(:,:,:) = 0.0
  END IF

  ! divergent u increments from SKEB2
  IF (stph_diag%L_skeb2_u_div) THEN
    ALLOCATE(stph_diag%skeb2_u_div(udims%i_start:udims%i_end,           &
                                   udims%j_start:udims%j_end,           &
                                   udims%k_start:udims%k_end))
    stph_diag%skeb2_u_div(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_u_div(1,1,1))
    stph_diag%skeb2_u_div(:,:,:) = 0.0
  END IF

  ! divergent v increments from SKEB2
  IF (stph_diag%L_skeb2_v_div) THEN
    ALLOCATE(stph_diag%skeb2_v_div(vdims%i_start:vdims%i_end,           &
                                   vdims%j_start:vdims%j_end,           &
                                   vdims%k_start:vdims%k_end))
    stph_diag%skeb2_v_div(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_v_div(1,1,1))
    stph_diag%skeb2_v_div(:,:,:) = 0.0
  END IF

  ! SKEB2: dissipation field from smagorinsky
  IF (stph_diag%L_skeb2_disp_smag) THEN
    ALLOCATE(stph_diag%skeb2_disp_smag(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end,       &
                                       pdims%k_start:pdims%k_end))
    stph_diag%skeb2_disp_smag(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_disp_smag(1,1,1))
    stph_diag%skeb2_disp_smag(:,:,:) = 0.0
  END IF

  ! SKEB2: dissipation field from convection
  IF (stph_diag%L_skeb2_disp_conv) THEN
    ALLOCATE(stph_diag%skeb2_disp_conv(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end,       &
                                       pdims%k_start:pdims%k_end))
    stph_diag%skeb2_disp_conv(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_disp_conv(1,1,1))
    stph_diag%skeb2_disp_conv(:,:,:) = 0.0
  END IF

  ! SKEB2: dissipation field from SKEB1 KE dissipation
  IF (stph_diag%L_skeb2_disp_skeb1) THEN
    ALLOCATE(stph_diag%skeb2_disp_skeb1(pdims%i_start:pdims%i_end,      &
                                        pdims%j_start:pdims%j_end,      &
                                        pdims%k_start:pdims%k_end))
    stph_diag%skeb2_disp_skeb1(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_disp_skeb1(1,1,1))
    stph_diag%skeb2_disp_skeb1(:,:,:) = 0.0
  END IF

  ! SKEB2: Smoothed Modulating field
  IF (stph_diag%L_skeb2_smodfield) THEN
    ALLOCATE(stph_diag%skeb2_smodfield(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end,       &
                                       pdims%k_start:pdims%k_end))
    stph_diag%skeb2_smodfield(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_smodfield(1,1,1))
    stph_diag%skeb2_smodfield(:,:,:) = 0.0
  END IF

  ! SKEB2: final stream function forcing field
  IF (stph_diag%L_skeb2_streamfunction) THEN
    ALLOCATE(stph_diag%skeb2_streamfunction(pdims%i_start:pdims%i_end,  &
                                            pdims%j_start:pdims%j_end,  &
                                            pdims%k_start:pdims%k_end))
    stph_diag%skeb2_streamfunction(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_streamfunction(1,1,1))
    stph_diag%skeb2_streamfunction(:,:,:) = 0.0
  END IF

  ! SKEB2: initial random pattern
  IF (stph_diag%L_skeb2_random_pattern) THEN
    ALLOCATE(stph_diag%skeb2_random_pattern(pdims%i_start:pdims%i_end,  &
                                            pdims%j_start:pdims%j_end,  &
                                            pdims%k_start:pdims%k_end))
    stph_diag%skeb2_random_pattern(:,:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_random_pattern(1,1,1))
    stph_diag%skeb2_random_pattern(:,:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of initial SF forcing
  IF (stph_diag%L_skeb2_ke_psif) THEN
    ALLOCATE(stph_diag%skeb2_ke_psif(pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end))
    stph_diag%skeb2_KE_psif(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_psif(1,1))
    stph_diag%skeb2_KE_psif(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of numerical diss
  IF (stph_diag%L_skeb2_ke_sdisp) THEN
    ALLOCATE(stph_diag%skeb2_ke_sdisp(pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end))
    stph_diag%skeb2_KE_sdisp(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_sdisp(1,1))
    stph_diag%skeb2_KE_sdisp(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of convection diss
  IF (stph_diag%L_skeb2_ke_cdisp) THEN
    ALLOCATE(stph_diag%skeb2_ke_cdisp(pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end))
    stph_diag%skeb2_KE_cdisp(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_cdisp(1,1))
    stph_diag%skeb2_KE_cdisp(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of 2nd convection diss
  IF (stph_diag%L_skeb2_ke_kdisp) THEN
    ALLOCATE(stph_diag%skeb2_ke_kdisp(pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end))
    stph_diag%skeb2_KE_kdisp(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_kdisp(1,1))
    stph_diag%skeb2_KE_kdisp(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of modulated SF forcing
  IF (stph_diag%L_skeb2_ke_m_psif) THEN
    ALLOCATE(stph_diag%skeb2_ke_m_psif(pdims%i_start:pdims%i_end,       &
                                       pdims%j_start:pdims%j_end))
    stph_diag%skeb2_KE_m_psif(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_m_psif(1,1))
    stph_diag%skeb2_ke_m_psif(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of total wind incr before SKEB2
  IF (stph_diag%L_skeb2_ke_prewindincr) THEN
    ALLOCATE(stph_diag%skeb2_ke_prewindincr(pdims%i_start:pdims%i_end,  &
                                            pdims%j_start:pdims%j_end))
    stph_diag%skeb2_ke_prewindincr(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_prewindincr(1,1))
    stph_diag%skeb2_ke_prewindincr(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of wind incr from SKEB2
  IF (stph_diag%L_skeb2_ke_windincr) THEN
    ALLOCATE(stph_diag%skeb2_ke_windincr(pdims%i_start:pdims%i_end,     &
                                         pdims%j_start:pdims%j_end))
    stph_diag%skeb2_ke_windincr(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_windincr(1,1))
    stph_diag%skeb2_ke_windincr(:,:) = 0.0
  END IF

  ! SKEB2: Vert Integ. KE of total wind incr after SKEB2
  IF (stph_diag%L_skeb2_ke_postwindincr) THEN
    ALLOCATE(stph_diag%skeb2_ke_postwindincr(                           &
                                     pdims%i_start:pdims%i_end,         &
                                     pdims%j_start:pdims%j_end))
    stph_diag%skeb2_ke_postwindincr(:,:) = 0.0
  ELSE
    ! Code to ALLOCATE unity arrays when not
    ! used (for portability)
    ALLOCATE(stph_diag%skeb2_ke_postwindincr(1,1))
    stph_diag%skeb2_ke_postwindincr(:,:) = 0.0
  END IF

END IF ! SF(0,35)

! -----------------------------------------------------------------
!  This bit is to calculate the spherical harm. forcing
!  compute the spectral weighting function g(n) - which contains
!  the power law  "chi= n**(-1.54)" that Glenn saw in the CRM -
!  new recommended value (used in SPBS at ECMWF) is -1.54
!  **************
!  A vertical variation of the pattern is done by changing the wave
!  coefficients with height at an angular speed dependent on wave
!  number. The coefficients evolve in time according to a 1st-order
!  Markov process expression.
! -----------------------------------------------------------------


! Perform work on PE=0 to avoid extra calls to random_number which
! affects bit reproducibility by changing the random seed according
! to the number of PEs
IF (firsttimestep == firsttimestep_true .OR.                            &
    firsttimestep == firsttimestep_crun) THEN

  ! In this version of SKEB2, alpha is not wavenumber dependent
  alpha= 1.0- EXP(-timestep/tau_skeb)

  ! Set maximum backscatter amount
  max_tot_backscat = 1.0e3 * tot_backscat

  ! Set backscatter level requiring damping of forcing pattern
  crit_tot_backscat = 1.0e4 * tot_backscat

  CALL backscatter_spectrum( gspect, nlim, stph_n1, stph_n2,            &
                             timestep, tot_backscat ,alpha)
END IF         ! firsttimestep

IF (((stphseed == 2 .OR. l_nrun_as_crun) .AND.                     & 
    timestep_number == 1 .AND.                                     &
    stph_skeb2_data_check == stph_skeb2_data_present) .OR.         &
    (firsttimestep == firsttimestep_crun .AND.                     &
     stph_skeb2_data_check == stph_skeb2_data_present) ) THEN

  IF (l_nrun_as_crun) THEN
      cmessage = 'l_nrun_as_crun is TRUE:' //  &
         'Reading random numbers from the dump (equivalent of stphseed=2)'
      icode = -1
      CALL ereport(RoutineName, icode, cmessage)
  END IF

  ! copy from dump header if available
  IF (mype == 0) THEN
    CALL umPrint('retrieving dpsidtc/dpsidts from dump header array', &
                 src=RoutineName)
    
    DO n = stph_n1, stph_n2
      DO m = 0, stph_n2
        i = m + (n-stph_n1) * (stph_n2 + 1)
        dpsidtc(m, n) = a_flddepc(fdc_skeb2_dpsidtc_start + i)
        dpsidts(m, n) = a_flddepc(fdc_skeb2_dpsidts_start + i)
      END DO
    END DO

  END IF
  ! Scatter restart coeffs to all processors
  sndcount = (stph_n2+1)*(stph_n2-stph_n1+1)
  CALL gc_rbcast(3247, sndcount, zero, nproc, icode, dpsidtc)
  CALL gc_rbcast(3248, sndcount, zero, nproc, icode, dpsidts)
  
END IF


CALL update_dpsidt( alpha, dpsidtc, dpsidts, stph_n1, stph_n2, nlim,    &
                    gspect, zero, icode, info, firsttimestep)
                    
! Maximum row_length between all the PEs.
max_row_length = 0

! The FFT coefficients are computed using a different decomposition
! than the base logical processed grid (LPG) decomposition, i.e. X/Y
! decomposition. The FFT domain is decomposed latitudinally between PES,
! i.e. only a Y decomposition

! The FFT domain is divided into 'nproc_y' horizontal
! stripes which have the same number of rows of the local LPG Y
! decomposition, i.e. 'rows' variable. This is done to reduce the data
! communication between PEs when moving data from the FFT
! decomposition to the LPG decomposition since the communication is
! happening only between PEs that are on the same LPG row,
! i.e. 'gc_proc_row_group'.

! Example LPG decomposition between 8 PEs in 4 X 2 (nproc_x X nproc_y) grid.
! The global size of the domain is 'global_row_length' X 'global_rows'.
! 
!         X -> global_row_length
! =================
! |   |   |   |   |
! | 4 | 5 | 6 | 7 |   
! |   |   |   |   |
! |   |   |   |   |
! =================   Y -> global_rows
! |   |   |   |   |
! | 0 | 1 | 2 | 3 |
! |   |   |   |   |
! |   |   |   |   |
! =================    

! Example FFT Latitude decomposition between 8 PEs with 2 (nproc_y)
! stripes and 4 (nproc_y) sub-stripes.
! The global size of the domain is 'nlim'*2 X 'nlat'.
!
!         X -> nlim*2
! =================
! |       7       |    
! -----------------
! |       6       |    
! -----------------   Stripe 2
! |       5       |    
! -----------------
! |       4       |
! =================   Y -> nlat
! |       3       |    
! -----------------
! |       2       |
! -----------------   Stripe 1
! |       1       |
! -----------------
! |       0       |
! =================    


! Compute the size and offset of the stripe, i.e. latitude decomposition
DO pe = 0, nproc-1
  ! Number of rows of the horizontal stripe where the specific PE is
  ! located, this is the same number of rows that the PE has in the
  ! LPG decomposition
  stripe_rows(pe) = g_blsize(2, fld_type_p, pe)
  
  ! If 'global_rows' is not even, 'nlat' is equal to 'global_rows' - 1. 
  ! We have to make sure that if we add the number of rows between all
  ! the stripes, this value is not greater than nlat. Thus the
  ! following line reduces the number of rows in the last stride of
  ! 1 if global_rows is not even. 
  ! The value "g_datastart(2, pe) - 1" is the "stripe displacement"
  ! which is used to keep track of the real latitudes offset of each
  ! stripe.
  stripe_rows(pe) = stripe_rows(pe) -                               &
    MAX((g_datastart(2, pe)-1)+stripe_rows(pe)-nlat,0)
  
  ! Maximum row_length between all the PEs.
  max_row_length = MAX(max_row_length, g_blsize(1, fld_type_p, pe))
END DO

! An horizontal stripe is further decomposed latitudinally between the
! PEs on the same LPG row. Thus for each stripe, 'nproc_x' sub-stripes
! are created which are composed of multiple horizontal rows of lenght
! 'nlim*2'. 

! Finally, the full FFT domain is decomposed in 'nproc' sub-stripes,
! one for each PE.

! Maximum number of rows between all the sub-stripes.
max_my_rows = 0

! Compute the size and offset of the sub-stripes
pe_displs(0) = 0
DO pe_j = 0, nproc_y-1
  DO pe_i = 0, nproc_x-1
    ! The rank 
    pe=pe_i+pe_j*nproc_x

    ! Compute the number of rows of the sub-stripe on each PE.
    pe_rows(pe) = stripe_rows(pe)/nproc_x  

    ! If the number of rows in a stripe cannot be evenly divided
    ! between all the sub-stripes, add one extra row in the last few
    ! sub-stripes as needed. Thus the sum of rows between all the
    ! sub-stripes is equal to the number of rows of the stripe.
    IF (pe_i >= nproc_x - MOD(stripe_rows(pe),nproc_x)) THEN
      pe_rows(pe) = pe_rows(pe) + 1
    END IF    
    
    ! Compute the real latitudes (global) offset of the sub-stripe on each PE
    IF (pe > 0) THEN
      pe_displs(pe) = SUM(pe_rows(0:pe-1))
    END IF

    ! Compute the offset of the sub-stripe inside the stripe on each PE.
    IF (pe_i == 0) THEN
      pe_offstr(pe) = 0
    ELSE
      pe_offstr(pe) = SUM(pe_rows(pe_j*nproc_x:pe-1))
    END IF
    
    ! Maximum rows between all the sub-stripes (PEs).
    max_my_rows = MAX(max_my_rows,pe_rows(pe))
  END DO
END DO

! The number of rows and real latitudes offset of the local sub-stripe (PE)
my_rows   = pe_rows(mype)
my_displs = pe_displs(mype)

! Allocate local data arrays
! Real space, Fourier Sin, Cos, Radius Coeffs
ALLOCATE (my_dpsidtr(0:stph_n2,stph_n1:stph_n2))
ALLOCATE (my_phi(0:stph_n2,stph_n1:stph_n2))
ALLOCATE (my_phishft(0:stph_n2,stph_n1:stph_n2))

! Different between top level and bottom level in skeb2.
kdiff = skeb2_toplev - skeb2_botlev + 1

! Allocate GC array
ALLOCATE (sendbuf(max_row_length,max_my_rows,                                 &
  skeb2_botlev:skeb2_toplev,0:nproc_x-1))
ALLOCATE (recvbuf(max_row_length,max_my_rows,                                 &
  skeb2_botlev:skeb2_toplev,0:nproc_x-1))

! The number of elements to send and receive
sendcount=max_row_length*max_my_rows*kdiff
recvcount=max_row_length*max_my_rows*kdiff

! Populate the array with the sequence of Legendre Polynomials read by
! the skeb_forcing method.
IF (first_atmstep_call) THEN    
  IF (.NOT. ALLOCATED(seqalf)) THEN
    ALLOCATE(seqalf(stph_n1:stph_n2,0:nlim,my_rows))
  END IF
  DO ilat = 1, my_rows
    mu= ilat + my_displs
    DO m1 = 0, 2*nlim, 2
      m  = m1/2
      nn = MAX(stph_n1,m)
      DO n = nn, stph_n2
        seqalf(n,m,ilat)= Ymn(mu,n,m)
      END DO
    END DO    
  END DO
END IF

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(n,m,k,dpsidt,kr,ilat,pe_j,pe_i,pe,       &
!$OMP& dstart,dsize,i,j,my_dpsidtc,my_dpsidts)                                 &
!$OMP& SHARED(nlim,nlat,stph_n1, stph_n2,my_dpsidtr,dpsidts,                   &
!$OMP&    dpsidtc,my_phi,my_phishft,skeb2_botlev,skeb2_toplev,tdims,           &
!$OMP&    eta_theta_levels,r_theta_levels,my_rows,seqalf,gridpos,sendbuf,      &
!$OMP&    stripe_rows,mype,rows,row_length,psif,sendcount,recvcount,           &
!$OMP&    gc_proc_row_group,info,nproc_x,pe_rows,pe_offstr,recvbuf,            &
!$OMP&    planet_radius) 

! Allocate private working arrays (one for each OpenMP thread)
IF (.NOT. ALLOCATED(dpsidt)) THEN
  ALLOCATE (dpsidt(0:2*nlim+1))
END IF

IF (.NOT. ALLOCATED(my_dpsidtc)) THEN
  ALLOCATE (my_dpsidtc(0:stph_n2,stph_n1:stph_n2))
END IF

IF (.NOT. ALLOCATED(my_dpsidts)) THEN
  ALLOCATE (my_dpsidts(0:stph_n2,stph_n1:stph_n2))
END IF

!  These private arrays have sections that may not be properly
!  initialised later
DO n = stph_n1, stph_n2
  DO m = n+1, stph_n2
    my_dpsidtc(m, n) = 0.0
    my_dpsidts(m, n) = 0.0
  END DO
END DO

! Initialised shared arrays
!$OMP DO SCHEDULE(STATIC)
DO n = stph_n1, stph_n2-1
  DO m = 0, stph_n2
    my_dpsidtr(m, n) = 0.0
    my_phishft(m, n) = 0.0
    my_phi(m, n) = 0.0
  END DO
END DO
!$OMP END DO

! Initialised shared arrays
!$OMP DO SCHEDULE(STATIC)
DO n = stph_n1, stph_n2
  DO m = 0, n
    my_dpsidtr(m, n) = SQRT(dpsidtc(m, n)**2 + dpsidts(m, n)**2)
    !    Determine angle from SIN & COS wave components (single step)
    my_phi(m, n) = ATAN2(dpsidts(m, n), dpsidtc(m, n))
    !    Max shift ranges from 0 <-> pi  for wavenumbers stph_n1 <-> stph_n2
    my_phishft(m, n) = (stph_n2 - MAX(m, n)) * pi/ (stph_n2 - stph_n1)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  ! Adjust coefficients in the vertical
  ! Level 1 = no change -> 12km Level = max change (=pi)
  !  cycles around above that level
  kr = eta_theta_levels(k) * (r_theta_levels(1,1,tdims%k_end) -         &
      planet_radius)/ 12000.0
  DO n = stph_n1, stph_n2
    DO m = 0, n

      my_dpsidtc(m, n) = my_dpsidtr(m, n) * COS(my_phi(m, n) +          &
                         kr * my_phishft(m, n))
      my_dpsidts(m, n) = my_dpsidtr(m, n) * SIN(my_phi(m, n) +          &
                         kr * my_phishft(m, n))

    END DO
  END DO

  DO ilat = 1, my_rows
    ! Calculates the FFT coefficients in dpsidt using the FFT domain
    CALL skeb_forcing( my_dpsidtc, my_dpsidts, stph_n1, stph_n2, nlim,  &
                       dpsidt, ilat, my_rows, seqalf)
    
    ! Copy the dpsidt data into the send buffer which is used to
    ! transfer the data from from the FFT decomposition into LPG
    ! decomposition.

    ! Attention 'dpsidt' has length 'nlim*2' which could be 1 element
    ! greater than 'global_row_length' if this last one is not
    ! even. However, this extra element is ignored and not copied into
    ! sendbuf.
    pe_j = gridpos(2)
    DO pe_i = 0, nproc_x-1
      ! The rank of the PE in the same LPG row.
      pe = pe_i + pe_j*nproc_x
      ! Compute the row_length (dsize) and X offset (dstart) of the data
      ! inside the stripe for each PE.
      dstart = g_datastart(1, pe)
      dsize  = g_blsize(1, fld_type_p, pe)

      sendbuf(1:dsize,ilat,k,pe_i) = dpsidt(dstart-1:dstart-1+dsize-1)      
    END DO
  END DO
  
  ! If 'global_rows' is not even, 'nlat' is less than 'global_rows',
  ! thus 'stripe_rows' could be less then 'rows'. Set the eventual
  ! untouched area of psif to zero.
  DO j = stripe_rows(mype) + 1, rows
    DO i = 1, row_length
      psif(i, j, k) = 0.0
    END DO
  END DO

END DO ! k
!$OMP END DO 

! Transfer the data at all the level from a FFT decomposition to LPG
! decomposition to between the PEs in the same LPG row.  
!$OMP MASTER
CALL mpl_alltoall(sendbuf, sendcount, mpl_real, &
  recvbuf, recvcount, mpl_real,gc_proc_row_group, info)
!$OMP END MASTER

!$OMP BARRIER

! Copy the received buffer to the local array 'psif' which has the
! normal LPG decomposition.  The variable 'row_length' is always less
! or equal than 'nlim*2'
pe_j = gridpos(2)
DO pe_i = 0, nproc_x-1
  ! The rank of the PE in the same LPG row.
  pe = pe_i + pe_j*nproc_x

!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_botlev, skeb2_toplev 
    DO ilat = 1, pe_rows(pe)
      
      ! Compute the row offset inside the stripe.
      j = ilat + pe_offstr(pe)
      
      DO i = 1, row_length
        psif(i, j, k) = recvbuf(i,ilat,k,pe_i)
      END DO
    END DO
  END DO
!$OMP END DO  NOWAIT
END DO

! Deallocate private working arrays
IF (ALLOCATED(dpsidt))     DEALLOCATE (dpsidt)
IF (ALLOCATED(my_dpsidtc)) DEALLOCATE(my_dpsidtc)
IF (ALLOCATED(my_dpsidts)) DEALLOCATE(my_dpsidts)

!$OMP END PARALLEL

! Deallocate GC arrays
IF (ALLOCATED(sendbuf))    DEALLOCATE(sendbuf)
IF (ALLOCATED(recvbuf))    DEALLOCATE(recvbuf)
 
! Deallocate work arrays
IF (ALLOCATED(my_dpsidtr)) DEALLOCATE(my_dpsidtr)
IF (ALLOCATED(my_phi)) DEALLOCATE(my_phi)
IF (ALLOCATED(my_phishft)) DEALLOCATE(my_phishft)

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(skeb2_botlev,skeb2_toplev,udims,vdims,pdims,    &
!$OMP&    r_u,r_v,psif,m_psif,skeb2_urot,skeb2_udiv,skeb2_vrot,         &
!$OMP&    skeb2_vdiv,delr_rho,r_theta_levels,rho,rho_r2,r_rho_levels,   &
!$OMP&    mass,darea,sdisp,cdisp,kdisp)

! Reset r_u and r_v in EG configuration
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_u(i,j,k)     = 0.0
    END DO
  END DO
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_v(i,j,k)     = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Set outside part of arrays to zero for StdOut and STASH consistency
!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, skeb2_botlev - 1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      psif(i,j,k)       = 0.0
      m_psif(i,j,k)     = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = udims%k_start, skeb2_botlev - 1
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      skeb2_urot(i,j,k) = 0.0
      skeb2_udiv(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = vdims%k_start, skeb2_botlev - 1
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      skeb2_vrot(i,j,k) = 0.0
      skeb2_vdiv(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_toplev + 1, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      psif(i,j,k)       = 0.0
      m_psif(i,j,k)     = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_toplev + 1, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      skeb2_urot(i,j,k) = 0.0
      skeb2_udiv(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_toplev + 1, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      skeb2_vrot(i,j,k) = 0.0
      skeb2_vdiv(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT


! Calculate Density, mass and dZ for each level (k denotes top of level)
! These values are used in vertical integrals of Convective Dissipation
! and Kinetic Energy diagnostics at the end of this subroutine
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      delr_rho(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
      rho(i,j,k) = rho_r2(i,j,k)/(r_rho_levels(i,j,k) *                 &
                   r_rho_levels(i,j,k))
      mass(i,j,k) = darea(i,j) * rho(i,j,k) * delr_rho(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Set outside part of arrays to zero

!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, skeb2_botlev - 1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rho(i,j,k)       = 0.0
      delr_rho(i,j,k)     = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_toplev + 1, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rho(i,j,k)       = 0.0
      delr_rho(i,j,k)     = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Initialise all dissipation fields. They are conditionally filled based
! on logicals l_skeb2_psisdisp, l_skeb2_psicdisp and l_skeb2_skeb1disp
!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, pdims%k_end
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      sdisp(i,j,k) = 0.0
      cdisp(i,j,k) = 0.0
      kdisp(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (printstatus  ==  prstatus_diag) THEN
  WRITE(umMessage,*)' ---------------------------------------- '
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'max(dpsidtc)= ',MAXVAL(dpsidtc)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'min(dpsidtc)= ',MINVAL(dpsidtc)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'max(dpsidts)= ',MAXVAL(dpsidts)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'min(dpsidts)= ',MINVAL(dpsidts)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'max(psif)= ',MAXVAL(psif)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES16.8)') 'min(psif)= ',MINVAL(psif)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)' ---------------------------------------- '
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***** CAPE")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: m^2")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(skeb2_cape)=",      &
       MAXVAL(skeb2_cape),"min(skeb2_cape)=",MINVAL(skeb2_cape)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,2I7,10X,A,2I7)') "maxloc(skeb2_cape)=",           &
       MAXLOC(skeb2_cape),"minloc(skeb2_cape)=",MINLOC(skeb2_cape)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(skeb2_cape)= ",                  &
       SUM(skeb2_cape)/SIZE(skeb2_cape)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("**SKEB2** up_flux")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: m^2")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(skeb2_up_flux)=",   &
       MAXVAL(skeb2_up_flux),"min(skeb2_up_flux)=",MINVAL(skeb2_up_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(skeb2_up_flux)=",        &
       MAXLOC(skeb2_up_flux),"minloc(skeb2_up_flux)=",                  &
       MINLOC(skeb2_up_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(skeb2_up_flux)= ",               &
      SUM(skeb2_up_flux)/SIZE(skeb2_up_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("**SKEB2** dwn_flux")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: m^2")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(skeb2_dwn_flux)=",  &
       MAXVAL(skeb2_dwn_flux),"min(skeb2_dwn_flux)=",                   &
       MINVAL(skeb2_dwn_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(skeb2_dwn_flux)=",       &
       MAXLOC(skeb2_dwn_flux),"minloc(skeb2_dwn_flux)=",                &
       MINLOC(skeb2_dwn_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(skeb2_dwn_flux)= ",              &
      SUM(skeb2_dwn_flux)/SIZE(skeb2_dwn_flux)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** darea")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: m^2")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(darea)=",           &
       MAXVAL(darea),"min(darea)=",MINVAL(darea)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,2I7,10X,A,2I7)') "maxloc(darea)=",                &
       MAXLOC(darea),"minloc(darea)=",MINLOC(darea)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(darea)= ",                       &
       SUM(darea)/SIZE(darea)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** rho")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: kg/m^3")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(rho)=",             &
       MAXVAL(rho),"min(rho)=",MINVAL(rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(rho)=",                  &
       MAXLOC(rho),"minloc(rho)=",MINLOC(rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(rho)= ",                         &
       SUM(rho)/SIZE(rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** dZ")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: m")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(delr_rho)=",        &
       MAXVAL(delr_rho),"min(delr_rho)=",MINVAL(delr_rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(delr_rho)=",             &
       MAXLOC(delr_rho),"minloc(delr_rho)=",MINLOC(delr_rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(delr_rho)= ",                    &
       SUM(delr_rho)/SIZE(delr_rho)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  IF (.NOT. ALLOCATED(mass_2d)) THEN
    ALLOCATE (mass_2d(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end))
  END IF
  ! Compute 2D Mass
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      mass_2d(i,j) = 0.0
    END DO
  END DO
  DO k = skeb2_botlev, skeb2_toplev
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        mass_2d(i,j) = mass_2d(i,j) + mass(i,j,k)
      END DO
    END DO
  END DO
  ! Zero mass at poles reset to unity (to avoid divide-by-zero)
  WHERE (mass_2d < 1.0e3) mass_2d = 1.0

  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** 2D Mass")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** UNITS: kg")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(mass_2d)=",         &
       MAXVAL(mass_2d),"min(mass_2d)=",MINVAL(mass_2d)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,2I7,10X,A,2I7)') "maxloc(mass_2d)=",              &
       MAXLOC(mass_2d),"minloc(mass_2d)=",MINLOC(mass_2d)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') "mean(mass_2d)= ",                     &
       SUM(mass_2d)/SIZE(mass_2d)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***********************************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  IF (ALLOCATED(mass_2d)) DEALLOCATE(mass_2d)
END IF


! --------------------------------------------------------------------
! Calculate energy dissipated from advection (using Smagorinsky
!  diffusion)
! --------------------------------------------------------------------


IF (skeb2_sdisp == type_bihm) THEN

  CALL biharm_diss( rows, row_length, n_rows                          &
,                   model_levels                                      &
,                   skeb2_botlev,skeb2_toplev                         &
,                   u, v, w, rho                                      &
,                   delta_lambda, delta_phi                           &
,                   sdisp                                             &
,                   first_atmstep_call, timestep                      &
                       )

  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,                                                  &
    '("***** 3D Energy Dissipation Rate per volume biharmonic")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF
END IF     ! type_bihm

! Compute Smagorinsky Dnum
IF (skeb2_sdisp == type_smag) THEN
  CALL skeb_smagorinsky(   rows, row_length, n_rows                     &
,                          model_levels                                 &
,                          u, v                                         &
,                          delta_lambda, delta_phi                      &
,                          sdisp                                        &
,                          first_atmstep_call                           &
                       )

  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'(A)') '***** 3D Energy Dissipation Rate per '      &
         // 'volume in skeb_smagorinsky'
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF

END IF     ! type_smag


! WRITE OUT SDISP STATEMENTS
IF (l_skeb2_psisdisp ) THEN
  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'("***** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(sdisp)=",         &
         MAXVAL(sdisp),"min(sdisp)=",MINVAL(sdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(sdisp)=",              &
         MAXLOC(sdisp),"minloc(sdisp)=",MINLOC(sdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15)') "mean(sdisp)= ",                     &
         SUM(sdisp)/SIZE(sdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("***********************************")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF
END IF

! --------------------------------------------------------------------
!  Calculate energy dissipated by convection (CAPE*mass_flux)
!  to re-scale the amount of energy feedback which is
!  10E-4 (see variable: tot_backscat) from the sph.harm. field
! --------------------------------------------------------------------

IF (l_skeb2_psicdisp) THEN
  IF (skeb2_cdisp  ==  type_cape) THEN
    ! ---------------------------------------------------------------
    ! Calculate modulating field (based on CAPE and with a vertical
    ! profile derived from dM/dz)
    ! ---------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(skeb2_botlev,skeb2_toplev,pdims,diff_flux,skeb2_up_flux,  &
!$OMP&    skeb2_dwn_flux,cape_kk,skeb2_cape,cdisp,conv_res_factor,g,    &
!$OMP&    delr_rho,rho)

    ! Calculate diff flux
!$OMP DO SCHEDULE(STATIC)
    DO k = skeb2_botlev - 1, skeb2_toplev + 1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          diff_flux(i,j,k) = skeb2_up_flux(i,j,k) -                     &
                             skeb2_dwn_flux(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! Remove negative/noise values in CAPE array
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        cape_kk(i,j) = skeb2_cape(i,j)
        IF(cape_kk(i,j) < 0.01) cape_kk(i,j) = 0.0
      END DO
    END DO
!$OMP END DO

    ! ---------------------------------------------------------------
    ! Calculate sm_tot_disp =  vertical profile from centred derivative
    !                     of dIF f_flux
    !                     Horizontal (xy) field from cape
    ! --------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
    DO k = skeb2_botlev, skeb2_toplev
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          !   Need to divide by (g*rho*dz) to convert Pa.s^-1 to s^-1
          !   Then: multiply by CAPE (J.kg^-1) => J.kg^-1.s^-1  or  m^2.s^-3
          !   Remember: delr_rho is a model-level thickness and diff_flux
          !             is differenced over 2 levels, => delr_rho@(k)+(k+1)
          cdisp(i,j,k) = conv_res_factor   *                            &
                       (ABS(diff_flux(i,j,k+1)-diff_flux(i,j,k-1))      &
                       /((delr_rho(i,j,k)+delr_rho(i,j,k+1))            &
                       *g*rho(i,j,k)))*cape_kk(i,j)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'(A)') '*** 3D Energy Dissipation Rate per '      &
           // 'volume in Convective M-FLX*CAPE'
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("*** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(cdisp)=",       &
           MAXVAL(cdisp),"min(cdisp)=",MINVAL(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(cdisp)=",            &
           MAXLOC(cdisp),"minloc(cdisp)=",MINLOC(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15)') "mean(cdisp)= ",                   &
           SUM(cdisp)/SIZE(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***********************************")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,*)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF

    !!!   END IF  ! (l_skeb2_cdisp_cape)

  ELSE IF (skeb2_cdisp  ==  type_mflx) THEN
    DO k = skeb2_botlev,skeb2_toplev
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          !   KE= 0.5 * m * v^2
          !     = 1/g*rho * d(upflux)/dz * (upflux/g*rho*alphac)**2
          cdisp(i,j,k)=(MAX(skeb2_up_flux(i,j,k+1) -                    &
                            skeb2_up_flux(i,j,k-1),0.0)/                &
                    (g*rho(i,j,k)*(delr_rho(i,j,k)+delr_rho(i,j,k+1)))) &
                       *(MAX(skeb2_up_flux(i,j,k),0.0)/                 &
                       (g*rho(i,j,k)*alphac))**2
        END DO
      END DO
    END DO

    IF (printstatus  >  prstatus_normal) THEN
      WRITE(umMessage,'(A)') '*** 3D Energy Dissipation Rate per '      &
           // 'volume in Convective M-FLX*CAPE => w'
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("*** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(cdisp)=",       &
           MAXVAL(cdisp),"min(cdisp)=",MINVAL(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(cdisp)=",            &
           MAXLOC(cdisp),"minloc(cdisp)=",MINLOC(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15)') "mean(cdisp)= ",                   &
           SUM(cdisp)/SIZE(cdisp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***********************************")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,*)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF

  END IF  ! (skeb2_cdisp)
END IF  ! (l_skeb2_psicdisp)


! -----------------------------------------------------------
! Calculate assumed energy dissipation (using SKEB1-type KE)
! This is similar to Smagorinsky turbulence, but focuses on
!   the Kinetic Energy field.
! -----------------------------------------------------------

IF (l_skeb2_skeb1disp) THEN
  CALL stph_skeb1( rows, row_length, n_rows                             &
,                          model_levels                                 &
,                          mass, c_dadt                                 &
,                          u, v                                         &
,                          kdisp                                        &
                       )

  IF (printstatus  >  prstatus_normal) THEN
    WRITE(umMessage,'(A)') '*** 3D Energy Dissipation Rate per '        &
         // 'volume in skeb1'
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("*** UNITS: m^2.s^-3 or kg.m^2.s^-2/s/kg")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') "max(kdisp)=",         &
         MAXVAL(kdisp),"min(kdisp)=",MINVAL(kdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,3I7,10X,A,3I7)') "maxloc(kdisp)=",              &
         MAXLOC(kdisp),"minloc(kdisp)=",MINLOC(kdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15)') "mean(kdisp)= ",                     &
         SUM(kdisp)/SIZE(kdisp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("***********************************")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF

END IF  ! (l_skeb2_skeb1disp)


! Calculate total dissipation field
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)        &
!$OMP& SHARED(skeb2_botlev,skeb2_toplev,pdims,sm_tot_disp,cdispfac,     &
!$OMP&    cdisp,sdispfac,sdisp,kdispfac,kdisp)
DO k = skeb2_botlev, skeb2_toplev
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      sm_tot_disp(i,j,k) = cdispfac*cdisp(i,j,k) +                      &
                  sdispfac*sdisp(i,j,k) + kdispfac*kdisp(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (printstatus  ==  prstatus_diag) THEN
  WRITE(umMessage,'("*********** Total Dissipation Field *************")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,F16.3)') 'Convection Diss. Factor (cdispfac) = ', &
       cdispfac
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,F16.3)') 'Numerical Diss. Factor (sdispfac) = ',  &
       sdispfac
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,F16.3)') 'SKEB1-KE Diss. Factor (kdispfac) = ',   &
       kdispfac
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
END IF

! Iterative cycling spreads isolated convective dissipation elements
! more evenly (number of iterations controlled by nsmooth)
IF (.NOT. l_skebsmooth_adv) THEN

  ! Limiter on energy backscattered
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i,j,k)        &
!$OMP& SHARED(pdims,sm_tot_disp,max_tot_backscat)
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (sm_tot_disp(i,j,k) > max_tot_backscat) THEN
          sm_tot_disp(i,j,k) = max_tot_backscat
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  IF (printstatus  ==  prstatus_diag) THEN
    !  Fill perimeter of array with zero for diagnostic print consistency
    DO k = pdims%k_start, skeb2_botlev - 1
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          sm_tot_disp(i,j,k) = 0.0
        END DO
      END DO
    END DO
    DO k = skeb2_toplev + 1, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          sm_tot_disp(i,j,k) = 0.0
        END DO
      END DO
    END DO
    WRITE(umMessage,'(A,ES22.15,12x,A,ES22.15)') 'max(sm_tot_disp)=',   &
         MAXVAL(sm_tot_disp),'min(sm_tot_disp)=',MINVAL(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(sm_tot_disp)=',        &
         MAXLOC(sm_tot_disp),'minloc(sm_tot_disp)=',MINLOC(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15)') 'mean(sm_tot_disp)= ',               &
         SUM(sm_tot_disp)/SIZE(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("*******************************************")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF

  !  Allocate smoothing buffer array
  IF (.NOT. ALLOCATED(sm_tot_disp0)) THEN
    ALLOCATE(sm_tot_disp0(pdims_s%i_start:pdims_s%i_end,                &
                          pdims_s%j_start:pdims_s%j_end,                &
                          pdims_s%k_start:pdims_s%k_end))
  END IF

  DO ismooth = 1, nsmooth

    IF (printstatus  ==  prstatus_diag) THEN
      WRITE(umMessage,*) 'starting ismooth = ', ismooth ,' of ', nsmooth
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF

    !  Fill buffer array with latest data and exchange halos
!$OMP  PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) PRIVATE(i,j,k)       &
!$OMP& SHARED(skeb2_botlev,skeb2_toplev,pdims,sm_tot_disp0,sm_tot_disp, &
!$OMP& at_extremity)
    DO k = skeb2_botlev,skeb2_toplev
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          sm_tot_disp0(i,j,k)=sm_tot_disp(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! DEPENDS ON: swap_bounds
    CALL swap_bounds( sm_tot_disp0, row_length, rows,                   &
                      skeb2_toplev, offx, offy, fld_type_p,             &
                      swap_field_is_scalar )


    !  2D 1-2-1 spatial filter
!$OMP  PARALLEL DO SCHEDULE (STATIC) DEFAULT(NONE) PRIVATE(i,j,k,jm1,  &
!$OMP  jp1,im1,ip1) SHARED(skeb2_botlev,skeb2_toplev,pdims,            &
!$OMP&    sm_tot_disp0,sm_tot_disp,at_extremity)
    DO k = skeb2_botlev, skeb2_toplev
      DO j = pdims%j_start, pdims%j_end
        jm1 = j - 1
        jp1 = j + 1
        IF (at_extremity(psouth)) jm1 = MAX(jm1, pdims%j_start)
        IF (at_extremity(pnorth)) jp1 = MIN(jp1, pdims%j_end)
        DO i = pdims%i_start, pdims%i_end
          im1 = i - 1
          ip1 = i + 1
          sm_tot_disp(i,j,k) = 0.0625*(sm_tot_disp0(im1,jp1,k) +        &
                                       sm_tot_disp0(ip1,jp1,k) +        &
                                       sm_tot_disp0(im1,jm1,k) +        &
                                       sm_tot_disp0(ip1,jm1,k) +        &
                                       2*( sm_tot_disp0(i,jp1,k) +      &
                                           sm_tot_disp0(im1,j,k) +      &
                                           sm_tot_disp0(ip1,j,k) +      &
                                           sm_tot_disp0(i,jm1,k) ) +    &
                                         4*sm_tot_disp0(i,j,k) )
        END DO  ! i
      END DO  ! j
    END DO  ! k
!$OMP END PARALLEL DO

    IF (printstatus  ==  prstatus_diag) THEN
      WRITE(umMessage,                                                  &
           '("*********** Smoothed Dissipation Field ***********")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15,12x,A,ES22.15)') 'max(sm_tot_disp)=', &
           MAXVAL(sm_tot_disp),'min(sm_tot_disp)=',MINVAL(sm_tot_disp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(sm_tot_disp)=',      &
           MAXLOC(sm_tot_disp),'minloc(sm_tot_disp)=',                  &
           MINLOC(sm_tot_disp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'(A,ES22.15)') 'mean(sm_tot_disp)= ',             &
           SUM(sm_tot_disp)/SIZE(sm_tot_disp)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***********************************")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF

  END DO   ! ismooth iteration

  ! DEALLOCATE smoothing buffer array
  DEALLOCATE(sm_tot_disp0)

ELSE

  ! Advanced smoothing option: using pre-defined smoothing array
  ! Removes backscatter in larger region of instability

    !  Allocate smoothing buffer array (using STPH-defined halo)
  IF (.NOT. ALLOCATED(sm_tot_disp0)) THEN
    ALLOCATE(sm_tot_disp0(stphdims_l%i_start:stphdims_l%i_end,          &
                          stphdims_l%j_start:stphdims_l%j_end,          &
                          stphdims_l%k_start:stphdims_l%k_end))
  END IF

  !  Allocate pattern smoothing array mask
  IF (.NOT. ALLOCATED(psif0)) THEN
    ALLOCATE(psif0(stphdims_l%i_start:stphdims_l%i_end,                 &
                   stphdims_l%j_start:stphdims_l%j_end))
  END IF

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,levfac,i1,i2,j1,j2,         &
!$OMP& local_crit_tot_backscat)                                         &
!$OMP& SHARED(pdims,psif0,crit_tot_backscat,sm_tot_disp,l_skebprint,    &
!$OMP&    max_tot_backscat,sm_tot_disp0,skeb2_botlev,skeb2_toplev,      &
!$OMP&    logscale,row_length,rows,offx_stph,offy_stph,stphdims_l,      &
!$OMP&    mask_pdamp,mask_smooth,at_extremity)

  ! Initialise mask to zero
  ! This is a 2-D array and will be set to one if instability is
  ! detected at any level in the column
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      psif0(i,j) = 0.0
    END DO  ! i
  END DO  ! j
!$OMP END DO

  !  Fill buffer array with latest data and exchange halos
  ! When excessive backscatter diagnosed, do not backscatter in this
  ! column by setting forcing pattern to zero in the halo region
  ! [offx_stph:offy_stph]
  ! ... useful diagnostic to find instabilities leading to GPStorms
!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_botlev,skeb2_toplev
    levfac = MAX(0.1, LOG10(logscale*k))   ! Max factor=10
    local_crit_tot_backscat = crit_tot_backscat/levfac
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ! Case of severe instability
        IF (sm_tot_disp(i,j,k) > local_crit_tot_backscat) THEN
          IF (l_skebprint) THEN
            WRITE(umMessage,'(A,ES12.5,A,I0,A,I0,A,I0,A,ES12.5)')       &
                  "WARNING: SKEB Backscatter = ",sm_tot_disp(i,j,k),    &
                  " at [",i,";",j,";",k,"] EXCEEDS limit of ",          &
                  local_crit_tot_backscat
            CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
          END IF
          sm_tot_disp(i,j,k) = 0.0
          psif0(i,j) = 1.0
        END IF
        sm_tot_disp0(i,j,k) = MIN(max_tot_backscat, sm_tot_disp(i,j,k))
        ! Reset to zero for iterative smoothing at next step
        sm_tot_disp(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP MASTER
  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( psif0, row_length, rows, 1,                         &
                    offx_stph, offy_stph, fld_type_p, swap_field_is_scalar  )
!$OMP END MASTER

!$OMP BARRIER

  ! Reduce dissipation in region around instability
  ! Area of influence can cross PE boundaries, so central point of
  ! smoothing can lie in the halo region.
!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_botlev, skeb2_toplev
    DO j = stphdims_l%j_start, stphdims_l%j_end
      DO i = stphdims_l%i_start, stphdims_l%i_end
        IF (psif0(i,j) > 0.5) THEN
          ! Apply at all SKEB levels
          i1 = MAX(pdims%i_start, i-offx_stph)
          i2 = MIN(pdims%i_end,   i+offx_stph)
          j1 = MAX(pdims%j_start, j-offy_stph)
          j2 = MIN(pdims%j_end,   j+offy_stph)
          DO jj = j1, j2
            DO ii = i1, i2
              sm_tot_disp0(ii,jj,k) = mask_pdamp(ii-i,jj-j) *           &
                                          sm_tot_disp0(ii,jj,k)
            END DO      ! ii
          END DO      ! jj
        END IF      ! Unstable point (can also be in halo)
      END DO      ! i
    END DO      ! j
  END DO      ! k
!$OMP END DO

!$OMP MASTER
  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( sm_tot_disp0, row_length, rows,                     &
                    skeb2_toplev, offx_stph, offy_stph, fld_type_p,     &
                    swap_field_is_scalar )
!$OMP END MASTER

!$OMP BARRIER

  !  2D 1-2-1 spatial filter (applied nsmooth times)
!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_botlev, skeb2_toplev
    !  Apply smoother
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        DO jj = -offy_stph, offy_stph
          DO ii = -offx_stph, offx_stph
            sm_tot_disp(i,j,k) = sm_tot_disp(i,j,k) +                   &
               mask_smooth(ii,jj) * sm_tot_disp0(ii+i,jj+j,k)
          END DO  ! ii
        END DO  ! jj
      END DO  ! i
    END DO  ! j
  END DO  ! SKEB levels
!$OMP END DO NOWAIT

  !  Fill perimeter of array with zero for diagnostic print consistency
!$OMP DO SCHEDULE(STATIC)
  DO k = pdims%k_start, skeb2_botlev - 1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        sm_tot_disp(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_toplev + 1, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        sm_tot_disp(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  IF (printstatus  ==  prstatus_diag) THEN
    WRITE(umMessage,'("***** Smoothed Dissipation Field *********")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15,12x,A,ES22.15)') 'max(sm_tot_disp)=', &
         MAXVAL(sm_tot_disp),'min(sm_tot_disp)=',MINVAL(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(sm_tot_disp)=',      &
         MAXLOC(sm_tot_disp),'minloc(sm_tot_disp)=',                  &
         MINLOC(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'(A,ES22.15)') 'mean(sm_tot_disp)= ',             &
         SUM(sm_tot_disp)/SIZE(sm_tot_disp)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("***********************************")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF

  ! DEALLOCATE smoothing buffer and pattern mask arrays
  DEALLOCATE(sm_tot_disp0)
  DEALLOCATE(psif0)

END IF   ! Smoothing option: l_skebsmooth_adv


! --------------------------------------------------------------------
!              Calculate u,v increments
! --------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,jp1,ip1,jv,jm1,im1)         &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,m_psif,psif,br,           &
!$OMP&    sm_tot_disp,tot_backscat,vdims,row_length,rows,offx,offy,     &
!$OMP&    work_z_halo,at_extremity,udims,skeb2_urot0,dy_theta,&
!$OMP&    l_skeb2_velpot,skeb2_udiv0,sin_theta_latitude,dx_theta,       &
!$OMP&    skeb2_vrot0,dx_v,skeb2_vdiv0,sin_v_latitude,dy_v,n_rows)

! Modulate streamfunction forcing field with smoothed E. Diss Rate
!
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      m_psif(i,j,k) = psif(i,j,k) *                                     &
                      SQRT(br*sm_tot_disp(i,j,k)/tot_backscat)
    END DO  ! i
  END DO    ! j
END DO      ! k
!$OMP END DO

! SWAP BOUNDS: used to communicate information from neighbouring
!              processors to make the smoothing correctly

!$OMP MASTER
! DEPENDS ON: swap_bounds
CALL swap_bounds( m_psif, row_length, rows,                             &
                  skeb2_toplev, offx, offy, fld_type_p,                 &
                  swap_field_is_scalar  )
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  DO j = vdims%j_start, vdims%j_end
    jp1 = j + 1
    DO i = vdims%i_start, vdims%i_end
      ip1 = i + 1
      work_z_halo(i,j,k) = 0.25*(                                       &
                       m_psif(i,j,k) + m_psif(ip1,j,k) +                &
                       m_psif(i,jp1,k) + m_psif(ip1,jp1,k))
    END DO  ! i
  END DO    ! j
  ! Set dissipation to zero on Pole row
  IF (at_extremity(psouth)) THEN
    DO i = vdims%i_start, vdims%i_end
      work_z_halo(i,vdims%j_start,k) = 0.0
    END DO  ! i
  END IF
  IF (at_extremity(pnorth)) THEN
    DO i = vdims%i_start, vdims%i_end
      work_z_halo(i,vdims%j_end,k) = 0.0
    END DO  ! i
  END IF
END DO      ! k
!$OMP END DO

!$OMP MASTER
! DEPENDS ON: swap_bounds
CALL swap_bounds( work_z_halo, row_length, n_rows,                      &
                  skeb2_toplev, offx, offy, fld_type_v,                 &
                  swap_field_is_scalar  )
!$OMP END MASTER

!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  ! Wind Incr on U-grid
  DO j = udims%j_start, udims%j_end
    jm1 = j - 1
    jv = j
    DO i = udims%i_start, udims%i_end
      ! -------------------------------------------------------------
      ! These are the ROTATIONAL incr :: from m_psif on vort points
      ! -------------------------------------------------------------
      skeb2_urot0(i,j,k)=(work_z_halo(i,jm1,k)-work_z_halo(i,jv,k))     &
                          /dy_theta(i,j)
    END DO    ! I
    IF (l_skeb2_velpot) THEN
      DO i = udims%i_start, udims%i_end
        ip1 = i + 1
        ! --------------------------------------------------
        ! These are the DIVERGENT increments
        !  Careful: +ve Vel Pot = conv in both hemispheres
        !           +ve StreamF = clockwise in both H
        !           Need (div)convergence in (anti)cyclone
        !           # sin(lat) < 0 in S.H. => (-1)*div-comp
        !           # using single column from sin_theta_latitude
        ! In future releases the Velocity Potential field
        !   will be calculated independently
        ! --------------------------------------------------
        skeb2_udiv0(i,j,k)=((m_psif(i,j,k)-m_psif(ip1,j,k))*  &
                 SIGN(1.0,sin_theta_latitude(1,j)))/ dx_theta(i,j)

      END DO    ! I
    END IF
  END DO     ! J
  ! Wind Incr on V-grid
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      im1 = i - 1
      skeb2_vrot0(i,j,k)=(work_z_halo(i,j,k)-work_z_halo(im1,j,k))      &
                          /dx_v(i,j)

    END DO    ! I
    IF (l_skeb2_velpot) THEN
      jp1 = j + 1
      jv = j
      ! Avoid using m_psif values beyond limits of P-grid at pole
      IF (at_extremity(pnorth)) jp1 = MIN( jp1, pdims%j_end)
      IF (at_extremity(psouth)) jv = MAX( jv, pdims%j_start)
      DO i = vdims%i_start, vdims%i_end
        skeb2_vdiv0(i,j,k)=((m_psif(i,jv,k)-m_psif(i,jp1,k))*           &
                 SIGN(1.0,sin_v_latitude(i,jv)))/ dy_v(i,jv)
      END DO    ! I
    END IF
  END DO     ! J
END DO      ! K
!$OMP END DO NOWAIT

!$OMP END PARALLEL

! Deallocate work arrays
IF (ALLOCATED(work_z_halo)) DEALLOCATE(work_z_halo)

! Use Multi-variable swap-bounds to exploit MPP
ifield = 0
ifield = ifield + 1
fields_to_swap(ifield) % field       => skeb2_urot0(:,:,:)
fields_to_swap(ifield) % field_type  =  fld_type_u
fields_to_swap(ifield) % levels      =  model_levels
fields_to_swap(ifield) % rows        =  rows
fields_to_swap(ifield) % vector      =  .TRUE.

ifield = ifield + 1
fields_to_swap(ifield) % field       => skeb2_vrot0(:,:,:)
fields_to_swap(ifield) % field_type  =  fld_type_v
fields_to_swap(ifield) % levels      =  model_levels
fields_to_swap(ifield) % rows        =  n_rows
fields_to_swap(ifield) % vector      =  .TRUE.

IF (l_skeb2_velpot) THEN
  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_udiv0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_u
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  rows
  fields_to_swap(ifield) % vector      =  .TRUE.

  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_vdiv0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_v
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  n_rows
  fields_to_swap(ifield) % vector      =  .TRUE.
END IF

CALL swap_bounds_mv( fields_to_swap, ifield, row_length, offx, offy )


! --------------------------------------------------------------------
!  Smoothing of the wind fields to reduce noise (1-2-1 filter)
! --------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,jm1,jp1,im1,ip1,levfac)     &
!$OMP& SHARED(udims,skeb2_botlev,skeb2_toplev,at_extremity,skeb2_urot,  &
!$OMP&     skeb2_urot0,timestep,l_skeb2_velpot,skeb2_udiv,skeb2_udiv0,  &
!$OMP&     vdims,skeb2_vrot,skeb2_vrot0,skeb2_vdiv,skeb2_vdiv0,         &
!$OMP&     level2km,printstatus,l_skebsmooth_adv,logscale)

!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  ! u-component
  DO j = udims%j_start, udims%j_end
    jm1= j - 1
    jp1= j + 1
    IF (at_extremity(psouth)) jm1 = MAX(jm1, udims%j_start)
    IF (at_extremity(pnorth)) jp1 = MIN(jp1, udims%j_end)
    DO i = udims%i_start, udims%i_end
      im1= i - 1
      ip1= i + 1
      ! -----------------------
      ! Rotational part
      ! -----------------------
      skeb2_urot(i,j,k)=0.0625*(skeb2_urot0(im1,jp1,k)+                 &
                                  skeb2_urot0(ip1,jp1,k)+               &
                                  skeb2_urot0(im1,jm1,k)+               &
                                  skeb2_urot0(ip1,jm1,k)+               &
                                  2*( skeb2_urot0(i,jp1,k)+             &
                                      skeb2_urot0(im1,j,k)+             &
                                      skeb2_urot0(ip1,j,k)+             &
                                      skeb2_urot0(i,jm1,k) )+           &
                                  4*skeb2_urot0(i,j,k) )*timestep
    END DO !I
    IF (l_skeb2_velpot) THEN
      DO i = udims%i_start, udims%i_end
        im1= i - 1
        ip1= i + 1
        ! -----------------------
        ! Divergent part
        ! -----------------------
        skeb2_udiv(i,j,k)=0.0625*(skeb2_udiv0(im1,jp1,k)+               &
                                      skeb2_udiv0(ip1,jp1,k)+           &
                                      skeb2_udiv0(im1,jm1,k)+           &
                                      skeb2_udiv0(ip1,jm1,k)+           &
                                      2*( skeb2_udiv0(i,jp1,k)+         &
                                          skeb2_udiv0(im1,j,k)+         &
                                          skeb2_udiv0(ip1,j,k)+         &
                                          skeb2_udiv0(i,jm1,k) )+       &
                                      4*skeb2_udiv0(i,j,k) )*timestep
      END DO !I
    ELSE
      ! This array is not initialised otherwise
      DO i = udims%i_start, udims%i_end
        skeb2_udiv(i,j,k)=0.0
      END DO !I
    END IF
  END DO !J
  ! v-component !
  DO j = vdims%j_start, vdims%j_end
    jm1= j - 1
    jp1= j + 1
    IF (at_extremity(psouth)) jm1 = MAX(jm1, vdims%j_start)
    IF (at_extremity(pnorth)) jp1 = MIN(jp1, vdims%j_end)
    DO i = vdims%i_start, vdims%i_end
      im1= i - 1
      ip1= i + 1
      ! -----------------------
      ! Rotational part
      ! -----------------------
      skeb2_vrot(i,j,k)=0.0625*(skeb2_vrot0(im1,jp1,k)+                 &
                                  skeb2_vrot0(ip1,jp1,k)+               &
                                  skeb2_vrot0(im1,jm1,k)+               &
                                  skeb2_vrot0(ip1,jm1,k)+               &
                                  2*( skeb2_vrot0(i,jp1,k)+             &
                                      skeb2_vrot0(im1,j,k)+             &
                                      skeb2_vrot0(ip1,j,k)+             &
                                      skeb2_vrot0(i,jm1,k) )+           &
                                 4*skeb2_vrot0(i,j,k) )*timestep
    END DO !I
    IF (l_skeb2_velpot) THEN
      DO i = vdims%i_start, vdims%i_end
        im1= i - 1
        ip1= i + 1
        ! -----------------------
        ! Divergent part
        ! -----------------------
        skeb2_vdiv(i,j,k)=0.0625*(skeb2_vdiv0(im1,jp1,k)+               &
                                    skeb2_vdiv0(ip1,jp1,k)+             &
                                    skeb2_vdiv0(im1,jm1,k)+             &
                                    skeb2_vdiv0(ip1,jm1,k)+             &
                                    2*( skeb2_vdiv0(i,jp1,k)+           &
                                        skeb2_vdiv0(im1,j,k)+           &
                                        skeb2_vdiv0(ip1,j,k)+           &
                                        skeb2_vdiv0(i,jm1,k) )+         &
                                    4*skeb2_vdiv0(i,j,k) )*timestep
      END DO
    ELSE
      ! This array is not initialised otherwise
      DO i = vdims%i_start, vdims%i_end
        skeb2_vdiv(i,j,k)=0.0
      END DO
    END IF
  END DO
  ! Set value of v-wind incr at pole = zero
  IF (at_extremity(psouth)) THEN
    DO i = vdims%i_start, vdims%i_end
      skeb2_vrot(i,vdims%j_start,k)=0.0
      skeb2_vdiv(i,vdims%j_start,k)=0.0
    END DO
  END IF
  IF (at_extremity(pnorth)) THEN
    DO i = vdims%i_start, vdims%i_end
      skeb2_vrot(i,vdims%j_end,k)=0.0
      skeb2_vdiv(i,vdims%j_end,k)=0.0
    END DO
  END IF
END DO !k
!$OMP END DO

! --------------------------------------------------------------------
!  Reduce wind increments to 0 at level=1 from level~2km using LOG10
!  This avoids creating large values near the top of the boundary lyr
! --------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k = skeb2_botlev, level2km
  levfac = MAX(0.0,LOG10(logscale*k))   ! log decrease to zero
  IF (printstatus  >  prstatus_oper) THEN
    WRITE(umMessage,'(" Factor at level (",I4,") = ",ES12.4)') k, levfac
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      skeb2_urot(i,j,k) = skeb2_urot(i,j,k) * levfac
      skeb2_udiv(i,j,k) = skeb2_udiv(i,j,k) * levfac
    END DO
  END DO
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      skeb2_vrot(i,j,k) = skeb2_vrot(i,j,k) * levfac
      skeb2_vdiv(i,j,k) = skeb2_vdiv(i,j,k) * levfac
    END DO
  END DO
END DO
!$OMP END DO

! Ramp off backscatter (linear) at top 3 levels of SKEB2 range
IF (l_skebsmooth_adv) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = skeb2_toplev-2, skeb2_toplev
    levfac = 0.25 * (skeb2_toplev - k + 1)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        skeb2_urot(i,j,k) = skeb2_urot(i,j,k) * levfac
        skeb2_udiv(i,j,k) = skeb2_udiv(i,j,k) * levfac
      END DO
    END DO
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        skeb2_vrot(i,j,k) = skeb2_vrot(i,j,k) * levfac
        skeb2_vdiv(i,j,k) = skeb2_vdiv(i,j,k) * levfac
      END DO
    END DO
  END DO
!$OMP END DO
END IF   ! Advanced smoothing option

!$OMP END PARALLEL

IF (printstatus  ==  prstatus_diag) THEN
  WRITE(umMessage,*)' ---------------------------------------- '
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)')                          &
       'max(R_u)=',MAXVAL(r_u),'min(R_u)=',MINVAL(r_u)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)')                                  &
       'maxloc(R_u)=',MAXLOC(r_u),'minloc(R_u)=',MINLOC(r_u)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(R_u)= ', SUM(ABS(r_u))/SIZE(r_u)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)')                          &
       'max(R_v)=',MAXVAL(r_v),'min(R_v)=',MINVAL(r_v)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)')                                  &
       'maxloc(R_v)=',MAXLOC(r_v),'minloc(R_v)=',MINLOC(r_v)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(R_v)= ', SUM(ABS(r_v))/SIZE(r_v)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') 'max(urot)=',            &
       MAXVAL(skeb2_urot),'min(urot)=',MINVAL(skeb2_urot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(urot)=',                 &
       MAXLOC(skeb2_urot),'minloc(urot)=',MINLOC(skeb2_urot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(urot)= ',                        &
       SUM(ABS(skeb2_urot))/SIZE(skeb2_urot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') 'max(vrot)=',            &
       MAXVAL(skeb2_vrot),'min(vrot)=',MINVAL(skeb2_vrot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(vrot)=',                 &
       MAXLOC(skeb2_vrot),'minloc(vrot)=',MINLOC(skeb2_vrot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(vrot)= ',                        &
       SUM(ABS(skeb2_vrot))/SIZE(skeb2_vrot)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') 'max(udiv)=',            &
       MAXVAL(skeb2_udiv),'min(udiv)=',MINVAL(skeb2_udiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(udiv)=',                 &
       MAXLOC(skeb2_udiv),'minloc(udiv)=',MINLOC(skeb2_udiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(udiv)= ',                        &
       SUM(ABS(skeb2_udiv))/SIZE(skeb2_udiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)') 'max(vdiv)=',            &
       MAXVAL(skeb2_vdiv),'min(vdiv)=',MINVAL(skeb2_vdiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)') 'maxloc(vdiv)=',                 &
       MAXLOC(skeb2_vdiv),'minloc(vdiv)=',MINLOC(skeb2_vdiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(vdiv)= ',                        &
       SUM(ABS(skeb2_vdiv))/SIZE(skeb2_vdiv)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)')                          &
       'max(sdisp)=',MAXVAL(sdisp),'min(sdisp)=',MINVAL(sdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)')                                  &
       'maxloc(sdisp)=',MAXLOC(sdisp),'minloc(sdisp)=',MINLOC(sdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(sdisp)= ',                       &
       SUM(ABS(sdisp))/SIZE(sdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)')                          &
       'max(cdisp)=',MAXVAL(cdisp),'min(cdisp)=',MINVAL(cdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)')                                  &
       'maxloc(cdisp)=',MAXLOC(cdisp),'minloc(cdisp)=',MINLOC(cdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(cdisp)= ',                       &
       SUM(ABS(cdisp))/SIZE(cdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

  WRITE(umMessage,'(A,ES22.15,12X,A,ES22.15)')                          &
       'max(kdisp)=',MAXVAL(kdisp),'min(kdisp)=',MINVAL(kdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,3I7,10X,A,3I7)')                                  &
       'maxloc(kdisp)=',MAXLOC(kdisp),'minloc(kdisp)=',MINLOC(kdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'(A,ES22.15)') 'mean(kdisp)= ',                       &
       SUM(ABS(kdisp))/SIZE(kdisp)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')

END IF


! SKEB2: Vert Integ. KE of total wind incr before SKEB2

IF (stph_diag%l_skeb2_ke_prewindincr) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j) = 0.0
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( r_u, row_length, rows,                              &
                    skeb2_toplev, offx, offy, fld_type_u,swap_field_is_vector  )

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( r_v, row_length, n_rows,                            &
                    skeb2_toplev, offx, offy, fld_type_v,swap_field_is_vector  )

  ! Interpolate squared velocity increments to pi-points
  DO k = skeb2_botlev, skeb2_toplev
    DO j = pdims%j_start, pdims%j_end
      jm1= j - 1
      ju = j
      jv = j
      DO i = pdims%i_start, pdims%i_end
        im1 = i - 1
        vert_int_work(i,j) = vert_int_work(i,j) + (0.5 * mass(i,j,k) *  &
                   (0.5*(r_u(im1,ju,k)**2 + r_u(i,ju,k)**2) +           &
                    0.5*(r_v( i,jm1,k)**2 + r_v(i,jv,k)**2)))
      END DO !i
    END DO !j
  END DO !k

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      stph_diag%skeb2_ke_prewindincr(i,j) = vert_int_work(i,j) *        &
                                             c_dadt(i,j)
    END DO
  END DO
END IF


! --------------------------------------------------------------------

! Pass increments from SKEB2 to R_u, R_v

! --------------------------------------------------------------------
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k)                         &
!$OMP& SHARED(skeb2_botlev,skeb2_toplev,udims,vdims,r_u,skeb2_urot,     &
!$OMP& skeb2_udiv,r_v,skeb2_vrot,skeb2_vdiv)                            &
!$OMP& SCHEDULE(STATIC)
DO k = skeb2_botlev, skeb2_toplev
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_u(i,j,k) = r_u(i,j,k) + skeb2_urot(i,j,k) + skeb2_udiv(i,j,k)
    END DO ! i
  END DO ! j
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_v(i,j,k) = r_v(i,j,k) + skeb2_vrot(i,j,k) + skeb2_vdiv(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k
!$OMP END PARALLEL DO

! Copy over increments to physics_tendencies_mod
IF (l_retain_stph_tendencies) THEN
  DO k = skeb2_botlev, skeb2_toplev
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du_stph(i,j,k) = skeb2_urot(i,j,k) + skeb2_udiv(i,j,k)
      END DO ! i
    END DO ! j
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv_stph(i,j,k) = skeb2_vrot(i,j,k) +  skeb2_vdiv(i,j,k)
      END DO ! i
    END DO ! j
  END DO  ! k
END IF ! end if over l_retain_stph_tendencies

! --------------------------------------------------------------------

! Output Stash Diagnostics

! --------------------------------------------------------------------


! u after skeb2
IF (stph_diag%L_skeb2_u) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        stph_diag%skeb2_u(i,j,k) = u(i,j,k) + r_u(i,j,k)
      END DO
    END DO
  END DO
END IF


! v after skeb2
IF (stph_diag%L_skeb2_v) THEN
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        stph_diag%skeb2_v(i,j,k) = v(i,j,k) + r_v(i,j,k)
      END DO
    END DO
  END DO
END IF


! u increment diagnostic
IF (stph_diag%L_skeb2_u_incr) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        stph_diag%skeb2_u_incr(i,j,k) = r_u(i,j,k) -                    &
                    stph_diag%skeb2_u_incr(i,j,k)
      END DO
    END DO
  END DO
END IF


! v increment diagnostic
IF (stph_diag%L_skeb2_v_incr) THEN
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        stph_diag%skeb2_v_incr(i,j,k) = R_v(i,j,k) -                    &
                    stph_diag%skeb2_v_incr(i,j,k)
      END DO
    END DO
  END DO
END IF


! rotational u increments from SKEB2
IF (stph_diag%l_skeb2_u_rot) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        stph_diag%skeb2_u_rot(i,j,k) = skeb2_urot(i,j,k)
      END DO
    END DO
  END DO
END IF


! rotational v increments from SKEB2
IF (stph_diag%l_skeb2_v_rot) THEN
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        stph_diag%skeb2_v_rot(i,j,k) = skeb2_vrot(i,j,k)
      END DO
    END DO
  END DO
END IF


! divergent u increments from SKEB2
IF (stph_diag%l_skeb2_u_div) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        stph_diag%skeb2_u_div(i,j,k) = skeb2_udiv(i,j,k)
      END DO
    END DO
  END DO
END IF


! divergent v increments from SKEB2
IF (stph_diag%l_skeb2_v_div) THEN
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        stph_diag%skeb2_v_div(i,j,k) = skeb2_vdiv(i,j,k)
      END DO
    END DO
  END DO
END IF


! dissipation from smagorinsky
IF (stph_diag%l_skeb2_disp_smag) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_disp_smag(i,j,k) = sdisp(i,j,k)
      END DO
    END DO
  END DO
END IF


! dissipation from convection
IF (stph_diag%l_skeb2_disp_conv) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_disp_conv(i,j,k) = cdisp(i,j,k)
      END DO
    END DO
  END DO
END IF


! Dissipation from SKEB1
IF (stph_diag%l_skeb2_disp_skeb1) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_disp_skeb1(i,j,k) = kdisp(i,j,k)
      END DO
    END DO
  END DO
END IF


! Smoothed dissipation field
IF (stph_diag%l_skeb2_smodfield) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_smodfield(i,j,k) = sm_tot_disp(i,j,k)
      END DO
    END DO
  END DO
END IF


! Streamfunction (modulated stream function forcing field)
IF (stph_diag%l_skeb2_streamfunction) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_streamfunction(i,j,k) = m_psif(i,j,k)
      END DO
    END DO
  END DO
END IF


! Streamfunction Forcing Field (initial from FFT)
IF (stph_diag%l_skeb2_random_pattern) THEN
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_random_pattern(i,j,k) = psif(i,j,k)
      END DO
    END DO
  END DO
END IF


! SKEB2: Mass-weighted Vert Integ. of modulated SF forcing field
IF (stph_diag%l_skeb2_ke_psif) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j) = 0.0
    END DO
  END DO
  DO k = skeb2_botlev, skeb2_toplev
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j)=vert_int_work(i,j) +                         &
                (m_psif(i,j,k) * mass(i,j,k))
      END DO
    END DO
  END DO
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      stph_diag%skeb2_ke_psif(i,j) = vert_int_work(i,j) * c_dadt(i,j)
    END DO
  END DO
END IF


! SKEB2: Mass-weighted Vert Integ. of numerical diss
IF (stph_diag%l_skeb2_ke_sdisp) THEN
  IF (.NOT. l_skeb2_psisdisp) THEN
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("**WARNING**: SKEB2 Numerical Dissipation off")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("  Requested STASH diagnostic (0,35,16) not available")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  ELSE
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,vert_int_work,sdisp,mass, &
!$OMP&     stph_diag,c_dadt)
   
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
! The main index of the loop is "j" so that it can be parallelised
! using OpenMP.
    DO j = pdims%j_start, pdims%j_end
      DO k = skeb2_botlev, skeb2_toplev
        DO i = pdims%i_start, pdims%i_end
          vert_int_work(i,j)=vert_int_work(i,j) +                       &
                 (sdisp(i,j,k) * mass(i,j,k))
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
        stph_diag%skeb2_ke_sdisp(i,j) = vert_int_work(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ! Print global mean value of Numerical Energy Dissipation
    IF (l_skebprint) THEN

      CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,     &
                          gltotke_tmp)

      gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
      WRITE(umMessage,*)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***** SKEB2 Global-total sdisp *****")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("tot(W/m^2)= ",ES22.15)') gltotke
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF    ! l_skebprint
  END IF
END IF


! SKEB2: Mass-weighted Vert Integ. of convection diss
IF (stph_diag%l_skeb2_ke_cdisp) THEN
  IF (.NOT. l_skeb2_psicdisp) THEN
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("**WARNING**: SKEB2 Convective Dissipation off")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("  Requested STASH diagnostic (0,35,17) not available")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  ELSE
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,vert_int_work,cdisp,mass, &
!$OMP&     stph_diag,c_dadt)
   
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
! The main index of the loop is "j" so that it can be parallelised
! using OpenMP.
    DO j = pdims%j_start, pdims%j_end
      DO k = skeb2_botlev, skeb2_toplev
        DO i = pdims%i_start, pdims%i_end
          vert_int_work(i,j)=vert_int_work(i,j) +                       &
                   (cdisp(i,j,k) * mass(i,j,k))
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
        stph_diag%skeb2_ke_cdisp(i,j) = vert_int_work(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ! Print global mean value of Numerical Energy Dissipation
    IF (l_skebprint) THEN
      CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,     &
                          gltotke_tmp)

      gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
      WRITE(umMessage,*)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***** SKEB2 Global-total cdisp *****")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("tot(W/m^2)= ",ES22.15)') gltotke
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF    ! l_skebprint
  END IF
END IF


! SKEB2: Mass-weighted Vert Integ. of SKEB1 KE dissipation
IF (stph_diag%l_skeb2_ke_kdisp) THEN
  IF (.NOT. l_skeb2_skeb1disp) THEN
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("**WARNING**: SKEB1 Dissipation l_skeb2_skeb1disp=F")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("  Requested STASH diagnostic (0,35,18) not available")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        stph_diag%skeb2_ke_kdisp(i,j) = 0.0
      END DO
    END DO
  ELSE
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,vert_int_work,kdisp,mass, &
!$OMP&     stph_diag,c_dadt)
   
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
! The main index of the loop is "j" so that it can be parallelised
! using OpenMP.
    DO j = pdims%j_start, pdims%j_end
      DO k = skeb2_botlev, skeb2_toplev
        DO i = pdims%i_start, pdims%i_end
          vert_int_work(i,j)=vert_int_work(i,j) +                       &
                 (kdisp(i,j,k) * mass(i,j,k))
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
        stph_diag%skeb2_ke_kdisp(i,j) = vert_int_work(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Print global mean value of Numerical Energy Dissipation
    IF (l_skebprint) THEN
      CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,     &
                          gltotke_tmp)

      gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
      WRITE(umMessage,*)
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("***** SKEB2 Global-total kdisp *****")')
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
      WRITE(umMessage,'("tot(W/m^2)= ",ES22.15)') gltotke
      CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    END IF    ! l_skebprint
  END IF
END IF


! SKEB2: Mass-weighted Vert Integ. of smoothed dissipation field
! Values for this diagnostic are printed to PE output for monitoring
! purposes
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)                            &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,vert_int_work,mass,       &
!$OMP&     c_dadt,sm_tot_disp)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    vert_int_work(i,j) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
! The main index of the loop is "j" so that it can be parallelised
! using OpenMP.
DO j = pdims%j_start, pdims%j_end
  DO k = skeb2_botlev, skeb2_toplev
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j)=vert_int_work(i,j) +                           &
               (sm_tot_disp(i,j,k) * mass(i,j,k))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

! Print global mean value of Numerical Energy Dissipation
IF (l_skebprint) THEN
  CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,         &
                      gltotke_tmp)

  gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
  WRITE(umMessage,*)
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("***** SKEB2 Global-total energy dissipation *****")')
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  WRITE(umMessage,'("tot(W/m^2)= ",ES22.15)') gltotke
  CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
END IF    ! l_skebprint
IF (stph_diag%l_skeb2_ke_m_psif) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      stph_diag%skeb2_ke_m_psif(i,j) = vert_int_work(i,j)
    END DO
  END DO
END IF


! SKEB2: Mass-weighted Vert Integ. KE of total wind incr after SKEB2
IF (stph_diag%l_skeb2_ke_postwindincr) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j) = 0.0
    END DO
  END DO

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( r_u, row_length, rows,                              &
                    skeb2_toplev, offx, offy, fld_type_u,swap_field_is_vector  )

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds( r_v, row_length, n_rows,                            &
                    skeb2_toplev, offx, offy, fld_type_v,swap_field_is_vector  )

  ! Interpolate squared windspeed (i.e. scalar field) to P-gridpoints
  !  (to avoid over-smoothing of a vector field)
  DO k = skeb2_botlev, skeb2_toplev
    DO j = pdims%j_start, pdims%j_end
      jm1 = j - 1
      ju = j
      jv = j
      DO i = pdims%i_start, pdims%i_end
        im1 = i - 1
        vert_int_work(i,j) = vert_int_work(i,j) +                       &
                   (0.5 * mass(i,j,k) *                                 &
                   (0.5*(r_u(im1,ju,k)**2 + r_u(i,ju,k)**2) +           &
                    0.5*(r_v( i,jm1,k)**2 + r_v(i,jv,k)**2)))
      END DO !I
    END DO !J
  END DO !K

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      stph_diag%skeb2_ke_postwindincr(i,j) = vert_int_work(i,j) *       &
                                             c_dadt(i,j)
    END DO
  END DO
END IF


! SKEB2: Vert Integ. KE of wind incr from SKEB2
IF (stph_diag%l_skeb2_ke_windincr) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k)                         &
!$OMP& SHARED(udims,vdims,skeb2_urot0,skeb2_udiv0,skeb2_vrot0,          &
!$OMP&     skeb2_vdiv0,skeb2_vdiv,skeb2_urot,skeb2_udiv,skeb2_vrot)     &
!$OMP  SCHEDULE(STATIC)
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        skeb2_urot0(i,j,k) = skeb2_urot(i,j,k)
        skeb2_udiv0(i,j,k) = skeb2_udiv(i,j,k)
      END DO
    END DO
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        skeb2_vrot0(i,j,k) = skeb2_vrot(i,j,k)
        skeb2_vdiv0(i,j,k) = skeb2_vdiv(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! Use Multi-variable swap-bounds to exploit MPP
  ifield = 0
  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_urot0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_u
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  rows
  fields_to_swap(ifield) % vector      =  .TRUE.

  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_udiv0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_u
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  rows
  fields_to_swap(ifield) % vector      =  .TRUE.

  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_vrot0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_v
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  n_rows
  fields_to_swap(ifield) % vector      =  .TRUE.

  ifield = ifield + 1
  fields_to_swap(ifield) % field       => skeb2_vdiv0(:,:,:)
  fields_to_swap(ifield) % field_type  =  fld_type_v
  fields_to_swap(ifield) % levels      =  model_levels
  fields_to_swap(ifield) % rows        =  n_rows
  fields_to_swap(ifield) % vector      =  .TRUE.

  CALL swap_bounds_mv( fields_to_swap, ifield, row_length,              &
                          offx, offy )

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i,j,k,jm1,ju,jv,im1)              &
!$OMP& SHARED(pdims,skeb2_botlev,skeb2_toplev,vert_int_work,mass,       &
!$OMP&     stph_diag,c_dadt,skeb2_urot0,skeb2_udiv0,skeb2_vrot0,        &
!$OMP&     skeb2_vdiv0,at_extremity)

  ! KE = 0.5*M*V**2, but because the u and v components are on different
  !      grids, we convert it all to the P-grid (see diagram at top)
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
! The main index of the loop is "j" so that it can be parallelised
! using OpenMP.
  DO j = pdims%j_start, pdims%j_end
    jm1 = j - 1
    ju = j
    jv = j
    DO k = skeb2_botlev, skeb2_toplev
      DO i = pdims%i_start, pdims%i_end
        im1= i - 1
        vert_int_work(i,j) = vert_int_work(i,j) +                       &
             (0.5 * mass(i,j,k) *                                       &
             (0.5*((skeb2_urot0(im1,ju,k)+ skeb2_udiv0(im1,ju,k))**2 +  &
                   (skeb2_urot0(i,ju,k)  + skeb2_udiv0(i,ju,k))**2)  +  &
              0.5*((skeb2_vrot0(i,jm1,k) + skeb2_vdiv0(i,jm1,k))**2  +  &
                   (skeb2_vrot0(i,jv,k)  + skeb2_vdiv0(i,jv,k))**2)))
      END DO !I
    END DO !K
  END DO !J
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      vert_int_work(i,j) = vert_int_work(i,j) * c_dadt(i,j)
      stph_diag%skeb2_ke_windincr(i,j) = vert_int_work(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  ! Print global mean value of Numerical Energy Dissipation
  IF (l_skebprint) THEN
    CALL global_2d_sums(vert_int_work, row_length, rows, 0, 0, 1,       &
                        gltotke_tmp)

    gltotke = gltotke_tmp(1)/(global_row_length*global_rows)
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("***** SKEB2 Global-total incr KE *****")')
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
    WRITE(umMessage,'("tot(W/m^2)= ",ES22.15)') gltotke
    CALL umPrint(umMessage,src='stph_skeb2-stph_skeb2')
  END IF
END IF

! Deallocate work arrays
IF (ALLOCATED(vert_int_work)) DEALLOCATE(vert_int_work)

IF (sf(0,35)) THEN
  CALL  diagnostics_stph( row_length, rows, model_levels,               &
                          n_rows, at_extremity, stph_diag,              &
                          stashwork35)
  ! ------------------------
  ! Tidy allocatable arrays
  ! ------------------------
  DEALLOCATE(stph_diag%skeb2_u)
  DEALLOCATE(stph_diag%skeb2_v)
  DEALLOCATE(stph_diag%skeb2_u_incr)
  DEALLOCATE(stph_diag%skeb2_v_incr)
  DEALLOCATE(stph_diag%skeb2_u_rot)
  DEALLOCATE(stph_diag%skeb2_v_rot)
  DEALLOCATE(stph_diag%skeb2_u_div)
  DEALLOCATE(stph_diag%skeb2_v_div)
  DEALLOCATE(stph_diag%skeb2_disp_smag)
  DEALLOCATE(stph_diag%skeb2_disp_conv)
  DEALLOCATE(stph_diag%skeb2_disp_skeb1)
  DEALLOCATE(stph_diag%skeb2_smodfield)
  DEALLOCATE(stph_diag%skeb2_streamfunction)
  DEALLOCATE(stph_diag%skeb2_random_pattern)
  DEALLOCATE(stph_diag%skeb2_ke_psif)
  DEALLOCATE(stph_diag%skeb2_ke_sdisp)
  DEALLOCATE(stph_diag%skeb2_ke_cdisp)
  DEALLOCATE(stph_diag%skeb2_ke_kdisp)
  DEALLOCATE(stph_diag%skeb2_ke_m_psif)
  DEALLOCATE(stph_diag%skeb2_ke_prewindincr)
  DEALLOCATE(stph_diag%skeb2_ke_windincr)
  DEALLOCATE(stph_diag%skeb2_ke_postwindincr)
END IF
! Deallcate skeb2_up_flux only if SPT is not called.
IF (.NOT. l_spt) THEN
  IF (ALLOCATED(skeb2_up_flux)) DEALLOCATE(skeb2_up_flux)
END IF

IF (ALLOCATED(skeb2_dwn_flux)) DEALLOCATE(skeb2_dwn_flux)
IF (ALLOCATED(skeb2_cape)) DEALLOCATE(skeb2_cape)

! Append dpsidtc/s to seed file (if active) at dump times.
! Generally used for CRUNs.
IF (ldump) THEN
  IF (mype == 0) THEN
    !copy data to dump header array
    CALL umPrint('stph_skeb2: Copying SKEB2 data to dump headers', &
                 src='stph_skeb2')
    DO n = stph_n1, stph_n2
      DO m = 0, stph_n2
        i = m + (n-stph_n1) * (stph_n2 + 1)
        a_flddepc(fdc_skeb2_dpsidtc_start + i) = dpsidtc(m, n) 
        a_flddepc(fdc_skeb2_dpsidts_start + i) = dpsidts(m, n) 
      END DO
    END DO
  END IF
  
END IF



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

! --------------------------------------------------------------------

END SUBROUTINE stph_skeb2
END MODULE stph_skeb2_mod
