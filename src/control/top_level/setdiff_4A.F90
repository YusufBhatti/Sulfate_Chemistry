! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Setdiff
!
! Purpose : to set up polar filtering and diffusion
!
! Language: FORTRAN 77 + common extensions also in Fortran 90.
! Programming standard; Unified Model Documentation Paper No. 3
! version 7.2, dated 5/2/98
!
! Documentation : Unified Model Documentation Paper No P0
!
!----------------------------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level
!
SUBROUTINE eg_setdiff(                                            &
! Input size and control variables
            global_rows, row_length, rows, n_rows, model_levels,        &
            at_extremity, datastart, offx, offy, mype, nproc_y,         &
            max_model_levels, max_121_rows, max_sweeps,                 &
! other constants
            delta_lambda, delta_phi,                                    &
            polar_cap, scale_ratio,                                     &
            ref_lat_deg, diff_coeff_ref,                                &
            cos_theta_latitude, sin_theta_latitude,                     &
            cos_v_latitude, sin_v_latitude,                             &
! Output data
             global_u_filter, global_v_filter,                          &
             u_sweeps, v_sweeps,                                        &
             u_begin, u_end, v_begin, v_end,                            &
             diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
             diffusion_coefficient_thermo, diffusion_order_thermo,      &
             diffusion_coefficient_wind, diffusion_order_wind,          &
             diff_order_thermo, diff_order_wind,                        &
             diff_timescale_thermo, diff_timescale_wind,                &
             L_sponge, sponge_power, sponge_ew, sponge_ns,              &
             sponge_wts_ew, sponge_wts_ns, L_subfilter_horiz,           &
             L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
             L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
             L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs,         &
             L_diff_wind_ew_only, halo_i, halo_j,xi2_p, xi2_v)


USE turb_diff_ctl_mod, ONLY: turb_ustart, turb_uend,              &
                             turb_vstart, turb_vend,              &
                             turb_phi_st, turb_phi_end
USE conversions_mod, ONLY: pi, pi_over_180
USE ereport_mod, ONLY: ereport
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength


USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


INTEGER, INTENT(IN) ::                                            &
       global_rows,                                               &
                           ! IN total number of rows in model
       row_length,                                                &
       rows,                                                      &
       n_rows,                                                    &
       model_levels,                                              &
       offx,                                                      &
       offy,                                                      &
       mype,                                                      &
       nproc_y                                                    &
, max_model_levels                                                &
                    ! Max no. model levels
, max_121_rows                                                    &
                ! Max no. of rows 1-2-1 filtered in a hemisphere
, max_sweeps                                                      &
                ! Max no. of sweeps of 1-2-1 filter
, diff_timescale_thermo                                           &
, diff_timescale_wind                                             &
, diff_order_thermo                                               &
, diff_order_wind                                                 &
, datastart(3)                                                    &
                     ! First gridpoints held by this processor
, diffusion_order_thermo(max_model_levels)                        &
, diffusion_order_wind(max_model_levels)                          &
, sponge_power                                                    &
, sponge_ew                                                       &
, sponge_ns

REAL, INTENT(IN) ::                                               &
  delta_lambda                                                    &
               ! EW (x) grid spacing
, delta_phi                                                       &
               ! NS (y) grid spacing
, scale_ratio                                                     &
, diffusion_coefficient_thermo(max_model_levels)                  &
, diffusion_coefficient_wind(max_model_levels)

LOGICAL, INTENT(IN) ::                                            &
  at_extremity(4)                                                 &
                   ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid
, L_subfilter_horiz                                               &
              ! switch for horiz turb diffusion
, L_diffusion                                                     &
              ! switch for old diffusion
, L_cdiffusion                                                    &
               ! switch for conservative diffusion
, L_pfincs                                                        &
           ! switch for polar filter
, L_pftheta                                                       &
                 ! switch for polar filter for theta
, L_pfw                                                           &
                 ! switch for polar filter for w
, L_pfuv                                                          &
                 ! switch for polar filter for horizontal winds
, L_pfexner                                                       &
                 ! switch for polar filter for Exner pressure
, L_diff_w                                                        &
                 ! switch for horiz. diffusion of w
, L_diff_auto                                                     &
                 ! T if automatic calculation of diffusion parms
, L_sponge       ! T if sponge zone active at lateral boundaries

REAL, INTENT(IN) ::                                               &
  sin_theta_latitude (row_length, rows)                           &
, sin_v_latitude (row_length, n_rows)                             &
, cos_theta_latitude(1-offx:row_length+offx,                      &
                     1-offy:rows+offy )                           &
, cos_v_latitude(1-offx:row_length+offx, 1-offy:n_rows+offy )

INTEGER, INTENT(IN) :: halo_i, halo_j

REAL, INTENT(IN) ::                                               &
  xi2_p(1-halo_j:rows+halo_j),                                    &
  xi2_v(-halo_j:n_rows-1+halo_j)


! Output arguments for diffusion/filtering control

LOGICAL, INTENT(INOUT) ::                                         &
  L_filter                                                        &
           ! switch for polar filter
, L_pfcomb                                                        &
                 ! switch for combined polar filter/diffusion
, L_diff_thermo                                                   &
                 ! switch for horiz. diffusion of theta
, L_diff_wind                                                     &
                 ! switch for horiz. diffusion of u,v 
, L_diff_wind_ew_only                                             &
                 ! switch to omit horiz. diffusion in NS dirn.
, L_diff_incs
                 ! switch for horiz. diffusion of incs

INTEGER, INTENT(INOUT) ::                                         &
  global_u_filter                                                 &
                   ! number of u rows filtered in a hemisphere
, global_v_filter                                                 &
                   ! number of v rows filtered in a hemisphere
, u_sweeps(max_121_rows)                                          &
                          ! sweeps for 1-2-1 filter
, v_sweeps(max_121_rows)                                          &
                          ! sweeps for 1-2-1 filter
, u_begin(0:max_121_rows)                                         &
                           ! row pointers for 1-2-1 filter
, u_end(0:max_121_rows)                                           &
                           ! row pointers for 1-2-1 filter
, v_begin(0:max_121_rows)                                         &
                           ! row pointers for 1-2-1 filter
, v_end(0:max_121_rows)    ! row pointers for 1-2-1 filter

REAL, INTENT(INOUT) ::                                            &
  diff_coeff_ref                                                  &
, polar_cap                                                       &
, ref_lat_deg                                                     &
, diff_coeff_phi                                                  &
                   ! NS diffusion coeff
!    diffusion coefficient for u rows
      , diff_coeff_u(1-offx:row_length+offx, 1-offy:rows+offy)          &
!    diffusion coefficient for v rows
      , diff_coeff_v(1-offx:row_length+offx, 1-offy:n_rows+offy)        &
      , sponge_wts_ew(sponge_ew)                                        &
      , sponge_wts_ns(sponge_ns)


! Local variables
REAL ::                                                           &
  angle                                                           &
, c1, c2                                                          &
, cos_lat1                                                        &
, cos_lat2                                                        &
, cos_ref_lat                                                     &
, denom                                                           &
, diff_new                                                        &
, diff_time_equ                                                   &
, diff_time_print                                                 &
, diff_coeff_thprint                                              &
, diff_coeff_uvprint                                              &
, diff_coeff_temp                                                 &
, dphi_half                                                       &
, exp1                                                            &
, filter_lat                                                      &
, first_lat                                                       &
, lat_new                                                         &
, Pi_over_2                                                       &
, pioverr                                                         &
, ref_lat_rad                                                     &
, rpower                                                          &
, scale_test                                                      &
, small                                                           &
, TINY                                                            &
, weight

! Local arrays
REAL ::                                                           &
 diff_scale(0:global_rows)

INTEGER ::                                                        &
  diff_order                                                      &
, diff_start_row                                                  &
, diff_timescale                                                  &
, filt_row                                                        &
, filt_NP                                                         &
, filt_SP                                                         &
, first_row                                                       &
, last_row                                                        &
, global_filter                                                   &
, half_grows                                                      &
, n_sweeps                                                        &
, power

INTEGER ::                                                        &
  sweeps(2 * max_121_rows)

LOGICAL ::                                                        &
  L_rows_even                                                     &
, L_diff_message                                                  &
, L_test

INTEGER ::                                                        &
        ! Mostly loop counters, but j also used for interp points
  i, j, gj, k                                                     &
, j_start_u, j_stop_u                                             &
, j_start_v, j_stop_v                                             &
, j_begin, j_end, j_stop                                          &
, js, je, icode

CHARACTER(LEN=*),PARAMETER :: RoutineName = 'EG_SETDIFF'
CHARACTER(LEN=errormessagelength)          :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! 1.0 Set filtering/diffusion parameters
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
L_filter = .FALSE.

IF ( L_pfincs .OR. L_pftheta .OR. L_pfuv .OR. L_pfw .OR. L_pfexner)&
    L_pfcomb = .TRUE.

IF (model_type /= mt_global) THEN
  IF ( L_pfcomb) THEN
    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,'(A)')                &
            ' POLAR filter is not allowed. Set false in namelist '
    CALL umPrint(umMessage,src='setdiff_4A')
  END IF ! L_pfcomb
  filter_lat = 0.0
END IF ! model_type /= mt_global

turb_ustart = 1
turb_uend = rows
turb_vstart = 1
turb_vend = n_rows
turb_phi_st = 1
turb_phi_end = rows

DO j = 0, max_121_rows
  u_begin(j) = -1
  v_begin(j) = -1
  u_end(j) = -3
  v_end(j) = -3
END DO ! j = 0, max_121_rows
DO j = 1, max_121_rows
  u_sweeps(j) = 0
  v_sweeps(j) = 0
END DO ! j = 0, max_121_rows

!  Main block for setting filtering and diffusion for global domain
IF (model_type == mt_global) THEN
  TINY = 0.000001
  small = 0.001
  Pi_over_2 = 0.5 * pi
  pioverr = pi / scale_ratio
  ref_lat_rad = ref_lat_deg * pi_over_180
  cos_ref_lat = COS(ref_lat_rad)
  dphi_half = 0.5 * delta_phi
  filter_lat = pi_over_180 * polar_cap
  global_filter = 0

  !  Initialise sweep boundaries for filtering and diffusion
  j_start_u = 1
  j_stop_u = rows
  j_start_v = 1
  j_stop_v = n_rows
  IF (at_extremity(PSouth)) THEN
    j_start_v = 2
  END IF ! at_extremity(PSouth))
  IF (at_extremity(PNorth)) THEN
    j_stop_v = n_rows - 1
  END IF ! at_extremity(PNorth)

  !  Diffusion coefficients calculated for 1 hemisphere only, so ....
  !  .... set pointers needed for switching hemispheres
  half_grows = global_rows / 2
  IF ( 2 * half_grows < global_rows ) THEN
    L_rows_even = .FALSE.
  ELSE !  2 * half_grows = global_rows
    L_rows_even = .TRUE.
  END IF !  2 * half_grows < global_rows

  ! Set diffusion coefficient = 0.25 on all rows . Loop limits
  !   in filtering routines will determine whether they are used or not.
  !   Values for diffusion will be calculated and overwritten later
  DO j = 1, rows + offy ! extend to halo - TomM
    DO i = 1-offx, row_length
      diff_coeff_u(i,j) = 0.25
    END DO
  END DO !    j = 1, rows
  DO j = 1, n_rows + offy ! extend to halo - TomM
    DO i = 1-offx, row_length
      diff_coeff_v(i,j) = 0.25
    END DO
  END DO !    j = 1, n_rows

  ! ----------------------------------------------------------------------
  ! 1.2 Initialise scanning structure for 1-2-1 filter
  !  Do the calculations on all processors to avoid communication
  !  NB ... thus cannot use pre-calculated trigs since polar values
  !  might not be on this processor
  ! ----------------------------------------------------------------------
  IF ( L_pfcomb ) THEN
    L_filter = .TRUE.
    diff_coeff_temp = 0.25
  ELSE ! no polar filter wanted
    ! calculate diffusion coefficients for all rows if no polar filter
    diff_coeff_temp = diff_coeff_ref
  END IF  !   L_pfcomb

  IF (L_filter) THEN
    !  Use positive latitudes (i.e. N Pole, S Pole is mirror image)
    first_lat = Pi_over_2 - dphi_half
    filter_lat = pi_over_180 * polar_cap
    ! If polar filter has been switched on make sure the end latitude
    !   after the row next to the pole
    IF ( filter_lat > first_lat ) THEN
      cos_lat1  = COS(first_lat)
      c1 = COS(pioverr * cos_lat1 / cos_ref_lat)
      lat_new = first_lat - dphi_half
      cos_lat2  = COS(lat_new)
      c2 = COS(pioverr * cos_lat2 / cos_ref_lat)
      diff_new = c2 - c1*c1
      sweeps(1) = 1

      IF ( diff_new > 0 ) THEN
        j = 2
        filter_lat  = first_lat
      ELSE ! diff_new < 0
        sweeps(2) = 2
        j = 3
        DO    !  cycle loop
          j = j + 1
          c1 = c2
          lat_new = lat_new - dphi_half
          cos_lat2  = COS(lat_new)
          c2 = COS(pioverr * cos_lat2 / cos_ref_lat)
          diff_new = c2 - c1 * c1
          IF ( diff_new > 0 ) EXIT
        END DO  !  cycle loop
      END IF ! diff_new > 0

      filter_lat  =  lat_new
      angle = filter_lat / pi_over_180
      global_filter = j - 1

      WRITE(umMessage,*) 'global_filter = ', global_filter
      CALL umPrint(umMessage,src='setdiff_4A')
      IF ( global_filter > 2 * max_121_rows) THEN
        cmessage = 'max_121_rows must be increased in CMAXSIZE '
        icode = 1
        CALL ereport(RoutineName, icode, cmessage)
      END IF ! global_filter > 2 * max_121_rows

      polar_cap = filter_lat / pi_over_180
      WRITE(umMessage,912) polar_cap
      CALL umPrint(umMessage,src='setdiff_4A')

      ! Determine the number of additional sweeps for each row
      power = 2
      ! First 2 rows done
      lat_new = filter_lat + dphi_half
      DO i = 3, global_filter
        c2 = c1
        lat_new = lat_new + dphi_half
        cos_lat1  = COS(lat_new)
        c1 = COS(pioverr * cos_lat1 / cos_ref_lat)
        DO   !  cycle loop
          power = power + 1
          IF ( power >= max_sweeps ) EXIT
          scale_test = c1 ** power
          IF ( c2 > scale_test ) EXIT
        END DO !   cycle loop
        sweeps(i) = power
        c1 = scale_test
      END DO  ! i = 3, global_filter

    ELSE  ! filter_lat set by user

      lat_new = first_lat
      global_filter  = 0
      ! Set number of sweeps for each row
      WRITE(umMessage,*) ' Setting polar filtering and diffusion settings'
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*) '     ****  ASSUMING V-AT-THE-POLES   ****'
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,920) polar_cap, first_lat / pi_over_180
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,921) polar_cap
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,'(A)')            &
' Set polar_cap = 90.0 in namelist if filtering is to be set automatically.'
      CALL umPrint(umMessage,src='setdiff_4A')
      power = 0
      DO
        IF ( lat_new < filter_lat ) EXIT
        global_filter = global_filter + 1
        IF ( global_filter > 2 * max_121_rows) THEN
          cmessage = 'max_121_rows must be increased in CMAXSIZE '
          icode = 1
          CALL ereport(RoutineName, icode, cmessage)
        ELSE
          power = power + 1
          sweeps(global_filter) = power
          lat_new  = lat_new - dphi_half
        END IF ! global_filter > 2 * max_121_rows
      END DO !   while ( diff_new < 0.0 )
      WRITE(umMessage,*) 'global_filter = ', global_filter
      CALL umPrint(umMessage,src='setdiff_4A')
      filter_lat  = first_lat - (global_filter - 1.0) * dphi_half
      polar_cap = filter_lat / pi_over_180
      WRITE(umMessage,912) polar_cap
      CALL umPrint(umMessage,src='setdiff_4A')

    END IF  ! filter_lat > first_lat

    ! Now set required values for u and v rows
    global_v_filter = global_filter / 2
    global_u_filter = global_filter - global_v_filter

    ! sweeps are cumulative to set as additional sweeps
    IF ( global_filter == 2 * global_v_filter ) THEN
      v_sweeps(1) = sweeps(1)
      u_sweeps(1) = sweeps(2)
      DO  i = 2, global_v_filter
        v_sweeps(i) = sweeps(2*i - 1) - v_sweeps(i-1)
      END DO ! i = 2, global_v_filter
      DO  i = 2, global_u_filter
        u_sweeps(i) = sweeps(2*i) - u_sweeps(i-1)
      END DO ! i = 2, global_u_filter
    ELSE IF ( global_v_filter > 0 ) THEN
      u_sweeps(1) = sweeps(1 )
      v_sweeps(1) = sweeps(2)
      DO  i = 2, global_u_filter
        u_sweeps(i) = sweeps(2*i - 1) - u_sweeps(i-1)
      END DO ! i = 1, global_u_filter
      DO  i = 2, global_v_filter
        v_sweeps(i) = sweeps(2*i) - v_sweeps(i-1)
      END DO ! i = 1, global_v_filter
    ELSE  !
      u_sweeps(1) = sweeps(1)
    END IF ! global_filter == 2 * global_v_filter

    !   Set loop limits on active processors.

    first_row = datastart(2)
    last_row = first_row + rows - 1
    filt_NP = global_rows - global_u_filter + 1

    IF ( first_row <= global_u_filter ) THEN
      ! S Hem rows need filtering

      ! Are we filtering across multiple processors ?
      IF ( last_row >= global_u_filter ) THEN
        ! northern-most processor row being filtered, set loop limits
        je = global_u_filter - first_row + 2
        j = 1
      ELSE ! last_row < global_u_filter
        ! other processor rows being filtered
        je = rows
        DO  j = 1, global_u_filter - last_row + 1
          u_begin(j) = 1
          u_end(j) = je
        END DO ! j = 1, global_u_filter - last_row + 1
        j = global_u_filter - last_row + 2
      END IF ! last_row > global_u_filter

      DO  ! cycle over rest of filter rows
        je = je - 1
        IF (je < 1 ) EXIT
        u_begin(j) = 1
        u_end(j) = je
        j = j + 1
      END DO  ! cycle over filter rows

    ELSE IF ( last_row >= filt_NP) THEN
      ! N Hem rows need filtering

      ! Are we filtering across multiple processors ?
      IF ( first_row <= filt_NP ) THEN
        ! southern-most processor row being filtered, set loop limits
        js = filt_NP - first_row
        j = 1
      ELSE ! first_row > filt_NP
        ! other processor rows being filtered
        js = 1
        DO  j = 1, first_row - filt_NP + 1
          u_begin(j) = js
          u_end(j) = rows
        END DO ! j = 1, first_row - filt_NP + 1
        j = first_row - filt_NP + 2

      END IF ! last_row > global_u_filter

      DO  ! cycle over rest of filter rows
        js = js + 1
        IF (js > rows ) EXIT
        u_begin(j) = js
        u_end(j) = rows
        j = j + 1
      END DO  ! cycle over filter rows

    END IF ! first_row <= global_u_filter

    !  Repeat for v rows

    last_row = first_row + rows - 1
    filt_SP = global_v_filter + 1
    filt_NP = global_rows - global_v_filter + 1

    IF ( first_row <= filt_SP ) THEN
      ! S Hem rows need filtering (but not at pole)
      js = 1
      IF ( first_row == 1 ) js = 2

      ! Are we filtering across multiple processors ?
      IF ( last_row >= filt_SP ) THEN
        ! northern-most processor row being filtered, set loop limits

        je = filt_SP - first_row + 2
        j = 1

      ELSE ! last_row < filt_SP
        ! other processor rows being filtered

        je = n_rows
        DO  j = 1, filt_SP - last_row + 1
          v_begin(j) = js
          v_end(j) = je
        END DO ! j = 1, filt_SP - last_row + 1
        j = filt_SP - last_row + 2

      END IF ! last_row > filt_SP

      DO  ! cycle over rest of filter rows
        je = je - 1
        IF (je < js ) EXIT
        v_begin(j) = js
        v_end(j) = je
        j = j + 1
      END DO  ! cycle over filter rows

    ELSE IF ( last_row >= filt_NP) THEN
      ! N Hem rows need filtering
      je = rows

      ! Are we filtering across multiple processors ?
      IF ( first_row <= filt_NP ) THEN
        ! southern-most processor row being filtered, set loop limits

        js = filt_NP - first_row
        j = 1

      ELSE ! first_row > filt_NP
        ! other processor rows being filtered

        js = 1
        DO  j = 1, first_row - filt_NP + 1
          v_begin(j) = js
          v_end(j) = je
        END DO ! j = 1, first_row - filt_NP + 1
        j = first_row - filt_NP + 2

      END IF ! last_row > global_u_filter

      DO  ! cycle over rest of filter rows
        js = js + 1
        IF (js > je ) EXIT
        v_begin(j) = js
        v_end(j) = je
        j = j + 1
      END DO  ! cycle over filter rows

    END IF ! first_row <= global_u_filter

    WRITE(umMessage,*) '  '
    CALL umPrint(umMessage,src='setdiff_4A')
    IF ( global_u_filter > max_121_rows ) THEN
      WRITE(umMessage,*) '  ******    WARNING      ******* '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*) ' You need to increase the PARAMETER '           &
             ,'max_121_rows in include file CMAXSIZE for array '  &
             , global_u_filter
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*)' OR set polar_cap closer to 90 degrees.'
      CALL umPrint(umMessage,src='setdiff_4A')
      icode = 1
      cmessage = "max_121_rows needs to be increased or "//       &
                   "polar_cap set closer to 90 degrees"
      CALL ereport(RoutineName, icode, cmessage)

    ELSE !   global_u_filter =< max_121_rows

      WRITE(umMessage,*) ' Polar filter details for processor ', mype
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*) ' global_u_filter = ', global_u_filter,          &
               ' global_v_filter = ', global_v_filter
      CALL umPrint(umMessage,src='setdiff_4A')
      DO i = 1, global_u_filter
        WRITE(umMessage,*) ' Pass ', i, ' sweeps = ', u_sweeps(i),        &
               ' u_begin = ', u_begin(i),' u_end = ', u_end(i)
        CALL umPrint(umMessage,src='setdiff_4A')
      END DO
      DO i = 1, global_v_filter
        WRITE(umMessage,*) ' Pass ', i, ' sweeps = ', v_sweeps(i),        &
               ' v_begin = ', v_begin(i),' v_end = ', v_end(i)
        CALL umPrint(umMessage,src='setdiff_4A')
      END DO
      WRITE(umMessage,910) polar_cap
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,911) global_v_filter
      CALL umPrint(umMessage,src='setdiff_4A')

    END IF !  global_u_filter > max_121_rows

    IF ( L_subfilter_horiz ) THEN
      !  Set turb_diff start and end rows if there is polar filtering
      IF ( u_end(1) > 0 .OR. v_end(1) > 0 ) THEN
        !  polar filtering on this processor

        IF ( v_begin(1) < 3 ) THEN
          !  S. Hem processor
          turb_ustart = v_end(1) + 1
          turb_uend = n_rows
        ELSE   ! v_begin(1) >= 3
          !  N. Hem processor
          turb_vstart = 1
          turb_vend = v_begin(1) - 1
        END IF ! v_begin(1) < 3
        IF ( u_begin(1) < 2 ) THEN
          !  S. Hem processor
          turb_ustart = u_end(1) + 1
          turb_uend = rows
          turb_phi_st = 2
        ELSE   ! u_begin(1) >= 2
          !  N. Hem processor
          turb_ustart = 1
          turb_uend = u_begin(1) - 1
          turb_phi_end = n_rows - 1
        END IF ! u_begin(1) < 2

      END IF ! u_end(1) > 0 .OR. v_end(1) > 0

    END IF ! L_subfilter_horiz

  END IF  !   L_filter

END IF !  model_type == mt_global

! ----------------------------------------------------------------------
! 1.4 Calculate diffusion coefficients given e-folding timesteps, order
! ----------------------------------------------------------------------
L_diff_message = .TRUE.
L_filter = .FALSE.
diff_order = MAX( diff_order_thermo, diff_order_wind )
IF ( diff_order < 1 ) THEN
  diff_timescale = 0
  L_diff_thermo = .FALSE.
  L_diff_wind = .FALSE.
  L_diff_incs = .FALSE.
END IF ! diff_order < 1

IF ( L_diff_thermo .OR. L_diff_wind .OR. L_diff_incs  ) THEN
  L_filter = .TRUE.
  !  Only 1st order diffusion allowed for combi diffusion
  diff_order = 1
  diff_timescale = MAX(diff_timescale_thermo,diff_timescale_wind)
  IF ( diff_timescale < 1 ) THEN
    diff_timescale = 1
    WRITE(umMessage,*)' Diffusion timescale was set to < 1 timestep'      &
    ,' Did you mean to do this? '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,'(2A)')                &
     ' Diffusion timescale has been reset to 1 timestep' &
    ,' but you might want to change this in the namelist'
    CALL umPrint(umMessage,src='setdiff_4A')
  END IF ! diff_timescale < 1
  IF ( diff_coeff_ref < TINY ) THEN
    WRITE(umMessage,'(2A)')                &
     ' Reference diffusion coefficient is too small or 0'&
    ,' Reset this in the namelist or switch off diffusion '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*)' Reference diffusion coefficient has been reset to'&
    ,' 0.000001 to avoid error message in deriving damping time.'
    CALL umPrint(umMessage,src='setdiff_4A')
    diff_coeff_ref = TINY
  END IF ! diff_coeff_ref < tiny
END IF ! L_diff_thermo .or. L_diff_wind

! Determine diffusion coefficients
IF (L_filter .AND. model_type == mt_global) THEN

  ! Initialise values for horizontal diffusion
  DO j = 1, global_filter
    diff_scale(j) = 0.25
  END DO ! j = 1, global_filter
  ! Set starting row for diffusion (1st after polar filter)
  diff_start_row = global_filter + 1

  IF ( L_diff_auto ) THEN

    IF ( global_filter > 0 ) THEN
      diff_start_row = global_filter
      WRITE(umMessage,*)'Horizontal diffusion will be applied equatorwards'&
             ,' of ',filter_lat / pi_over_180,' degrees'
      CALL umPrint(umMessage,src='setdiff_4A')
    ELSE ! no polar filter
      diff_start_row = 1
      lat_new = Pi_over_2 - dphi_half
      DO WHILE ( lat_new > filter_lat )
        diff_scale(diff_start_row) = 0.0
        lat_new = lat_new - dphi_half
        diff_start_row = diff_start_row + 1
      END DO !  while  lat_new > filter_lat
      filter_lat = lat_new
      WRITE(umMessage,*)'Horizontal diffusion will be applied equatorwards'&
             ,' of ',filter_lat / pi_over_180,' degrees'
      CALL umPrint(umMessage,src='setdiff_4A')
      diff_scale(diff_start_row) = diff_coeff_temp
    END IF ! global_filter > 0

    cos_lat1 = COS(filter_lat)
    c1 = SIN( pioverr * cos_lat1 / cos_ref_lat )
    lat_new = filter_lat - dphi_half
    j = diff_start_row
    ! This loop scales relative to ref_lat
    IF ( lat_new < ref_lat_rad ) THEN
      diff_scale(j) = diff_coeff_temp
      cos_lat2 = cos_lat1
      WRITE(umMessage,*)'  row    lat    ref_lat_rad    diffusion coefficient'
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*)j, lat_new, ref_lat_rad, diff_scale(j)
      CALL umPrint(umMessage,src='setdiff_4A')
    ELSE
      WRITE(umMessage,*)'row      lat             diffusion coefficient '
      CALL umPrint(umMessage,src='setdiff_4A')
    END IF ! lat_new < ref_lat_rad
    DO WHILE ( lat_new > ref_lat_rad - small )
      j = j + 1
      IF ( lat_new < 0.0 ) lat_new = 0.0
      cos_lat2 = COS(lat_new)
      c2 = SIN( pioverr * cos_lat2 / cos_ref_lat )
      diff_scale(j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
      WRITE(umMessage,*)j, lat_new, diff_scale(j)
      CALL umPrint(umMessage,src='setdiff_4A')
      lat_new = lat_new - dphi_half
    END DO !  while ( lat_new > ref_lat_rad - small )

    c2 = SIN(pioverr)
    denom = 1.0 / ( c2 * c2 )
    ! This loop scales r-grid wave to the same physical scale
    ! as next row polewards (if there are any rows left)
    IF ( lat_new > 0.0 - small ) THEN
      WRITE(umMessage,*)'  Now match successive rows '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*)'  row    lat      diff_scale '
      CALL umPrint(umMessage,src='setdiff_4A')
      DO WHILE ( lat_new > 0.0 - small )
        j = j + 1
        cos_lat1 = cos_lat2
        IF ( lat_new < 0.0 ) lat_new = 0.0
        cos_lat2 = COS(lat_new)
        c1 = SIN(pioverr * cos_lat1 /cos_lat2)
        diff_scale(j) = diff_scale(j - 1) * c1 * c1 * denom
        WRITE(umMessage,*)j, lat_new, diff_scale(j)
        CALL umPrint(umMessage,src='setdiff_4A')
        lat_new = lat_new - dphi_half
      END DO ! while ( lat_new > 0.0 - small )
    END IF  ! lat_new > 0.0 - small

    !  NS diffusion coefficient is the scaled relative to
    !EW diffusion at equator
    diff_coeff_phi = diff_scale(j) * delta_phi * delta_phi /      &
                                     (delta_lambda * delta_lambda)
    !  Now copy appropriate values for this processor
    v_begin(0) = -1
    v_end(0) = -2
    IF ( sin_v_latitude(1,1) + TINY < 0.0 ) THEN
      !  First row is in Southern  hemisphere
      IF ( v_begin(1) > 0 ) THEN
        j_begin = v_begin(1)
      ELSE
        j_begin = j_start_v
      END IF ! v_begin(1) > 0
      j_end = j_stop_v
      j_stop = j_end

      IF ( nproc_y == 1 ) j_stop = half_grows + 1

    ELSE ! sin_v_latitude(1,1) + tiny > 0.0
      !  First row is in Northern  hemisphere
      IF ( v_end(1) > 0 ) THEN
        j_end = v_end(1)
      ELSE
        j_end = j_stop_v
      END IF ! v_begin(1) > 0
      j_begin = j_start_v
      j_stop = j_end
    END IF ! sin_v_latitude(1,1) + tiny < 0.0

    L_test = .FALSE.
    j = j_begin - 1
    DO !  j = j_begin, j_stop
      j = j + 1
      IF ( j > j_stop ) EXIT
      IF ( sin_v_latitude(1,j) > 0.0 ) EXIT
      gj = 2 * ( datastart(2) + j - 2)
      IF ( v_begin(0) < 1 .AND. diff_scale(gj) > TINY ) THEN
        v_begin(0) = j
        v_end(0) = j_stop
        L_test = .TRUE.
      END IF ! v_begin(0) < 1 .and. diff_scale(gj) > tiny
      IF ( L_test .AND. diff_scale(gj) < TINY ) THEN
        v_end(0) = j - 1
        L_test = .FALSE.
      END IF ! L_test .and. diff_scale(gj) < tiny
      DO i = 1-offx, row_length

        ! Overwrite using ENDGame xi_2 array
        IF ( ABS(xi2_v(j-1)) <  filter_lat ) THEN
          c2 = SIN( pioverr * COS(xi2_v(j-1)) / cos_ref_lat )
          diff_coeff_v(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
        ELSE ! Polar filter
          diff_coeff_v(i,j) = 0.25
        END IF

      END DO
    END DO  !  j cycle
    !  now do N. Hem values if needed (loop will exit if not)
    j = j - 1  ! reset row counter
    DO !  j = j_stop + 1, j_end
      j = j + 1
      IF ( j > j_end ) EXIT
      gj = 2 * (global_rows - datastart(2) - j + 1)
      IF ( v_begin(0) < 1 .AND. diff_scale(gj) > TINY ) THEN
        v_begin(0) = j
        v_end(0) = j_stop_v
        L_test = .TRUE.
      END IF ! v_begin(0) < 1 .and. diff_scale(gj) > tiny
      IF ( L_test .AND. diff_scale(gj) < TINY ) THEN
        v_end(0) = j - 1
        L_test = .FALSE.
      END IF ! L_test .and. diff_scale(gj) < tiny
      DO i = 1-offx, row_length

        IF ( ABS(xi2_v(j-1)) <  filter_lat ) THEN
          ! Overwrite using ENDGame xi_2 array
          c2 = SIN( pioverr * COS(xi2_v(j-1)) / cos_ref_lat )
          diff_coeff_v(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
        ELSE ! Polar filter
          diff_coeff_v(i,j) = 0.25
        END IF
      END DO
    END DO  !  j cycle

    ! Reset v_begin(0) and v_end(0) to determine rows for NS diffusion only
    v_begin(0) = j_start_v
    v_end(0) = j_stop_v
    WRITE(umMessage,*)'On processor ', mype
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*)'v_begin(0) = ', v_begin(0),' v_end(0) = ', v_end(0)
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*)'v_begin(1) = ', v_begin(1),' v_end(1) = ', v_end(1)
    CALL umPrint(umMessage,src='setdiff_4A')

    ! EW diffusion rows are contained within NS diffusion rows and row
    ! limits are combined with polar filtering rows

    !  when  nproc_y == 1  S hem filtered/diffused first then N.Hem
    ! but   diff_coeff_v  contains coeff for all rows
    IF ( nproc_y == 1) THEN
      v_end(0) = j_stop
    END IF ! nproc_y == 1

    u_begin(0) = -1
    u_end(0) = -2
    IF ( sin_theta_latitude(1,1) + TINY < 0.0 ) THEN
      !  First row is in Southern  hemisphere
      IF ( u_begin(1) > 0 ) THEN
        j_begin = u_begin(1)
      ELSE
        j_begin = 1
      END IF ! u_begin(1) > 0
      j_end = rows
      j_stop = j_end
      IF ( nproc_y == 1 ) THEN
        IF ( L_rows_even) THEN
          j_stop = half_grows
        ELSE ! global_rows is odd
          j_stop = half_grows + 1
        END IF !  L_rows_even
      END IF ! nproc_y == 1

    ELSE   !( sin_theta_latitude(1,1) + tiny > 0.0
      !  First row is in Northern  hemisphere
      IF ( u_end(1) > 0 ) THEN
        j_end = u_end(1)
      ELSE
        j_end = rows
      END IF ! u_end(1) > 0
      j_begin = 1
      j_stop = j_end
    END IF !  sin_theta_latitude(1,1) + tiny < 0.0

    L_test = .FALSE.
    j = j_begin - 1
    DO !  j = j_begin, j_stop
      j = j + 1
      IF ( j > j_stop ) EXIT
      IF ( sin_theta_latitude(1,j) > 0.0 ) EXIT
      gj = 2 * ( datastart(2) + j - 1) - 1
      IF ( u_begin(0) < 1 .AND. diff_scale(gj) > TINY ) THEN
        u_begin(0) = j
        u_end(0) = j_stop
        L_test = .TRUE.
      END IF ! u_begin(0) < 1 .and. diff_scale(gj) > tiny
      IF ( L_test .AND. diff_scale(gj) < TINY ) THEN
        u_end(0) = j - 1
        L_test = .FALSE.
      END IF ! L_test .and. diff_scale(gj) < tiny
      DO i = 1-offx, row_length

        IF ( ABS(xi2_p(j)) <  filter_lat ) THEN
          ! Overwrite using ENDGame xi_2 array
          c2 = SIN( pioverr * COS(xi2_p(j)) / cos_ref_lat )
          diff_coeff_u(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
        ELSE ! Polar filter
          diff_coeff_u(i,j) = 0.25
        END IF

      END DO
    END DO !  j cycle
    !  now do N. Hem values if needed (loop will exit if not)
    j = j - 1  ! reset row counter
    DO !  j = j_stop + 1, j_end
      j = j + 1
      IF ( j > j_end ) EXIT
      gj = 2 * (global_rows - datastart(2) - j) + 1
      IF ( u_begin(0) < 1 .AND. diff_scale(gj) > TINY ) THEN
        u_begin(0) = j
        u_end(0) = rows
        L_test = .TRUE.
      END IF ! u_begin(0) < 1 .and. diff_scale(gj) > tiny
      IF ( L_test .AND. diff_scale(gj) < TINY ) THEN
        u_end(0) = j - 1
        L_test = .FALSE.
      END IF ! L_test .and. diff_scale(gj) < tiny
      DO i = 1-offx, row_length

        IF ( ABS(xi2_p(j)) <  filter_lat ) THEN
          ! Overwrite using ENDGame xi_2 array
          c2 = SIN( pioverr * COS(xi2_p(j)) / cos_ref_lat )
          diff_coeff_u(i,j) = diff_coeff_temp * c1 * c1 / ( c2 * c2 )
        ELSE ! Polar filter
          diff_coeff_u(i,j) = 0.25
        END IF

      END DO
    END DO  !  j cycle

    !           WRITE(6,*) ' '
    !           Do j = 1-offy,rows+offy
    !             WRITE(6,*) 'j, diff_u = ', j,diff_coeff_u(11-offx,j)
    !           End Do
    !           WRITE(6,*) ' '
    !           Do j = 1-offy,n_rows+offy
    !             WRITE(6,*) 'j, diff_v = ', j,diff_coeff_v(11-offx,j)
    !           End Do
    !           WRITE(6,*) ' '
    !           Do j = 1,global_rows
    !             WRITE(6,*) 'j, diff_scale = ', j,diff_scale(j)
    !           End Do
    !           WRITE(6,*) ' '
    !           gj = 1
    !           Do j = 1,rows
    !             WRITE(6,*) ' v ',j,diff_coeff_v(11-offx,j),diff_scale(gj)
    !             gj=gj+1
    !             WRITE(6,*) ' u ',j,diff_coeff_u(11-offx,j),diff_scale(gj)
    !             gj=gj+1
    !           End Do
    !           WRITE(6,*) ' v ',j,diff_coeff_v(11-offx,n_rows),diff_scale(gj)
    !           WRITE(6,*) ' '
    !           flush(6)

    ! Reset u_begin(0) and u_end(0) to determine rows for NS diffusion only
    u_begin(0) = j_start_u
    u_end(0) = j_stop_u
    WRITE(umMessage,*)'u_begin(0) = ', u_begin(0),' u_end(0) = ', u_end(0)
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*)'u_begin(1) = ', u_begin(1),' u_end(1) = ', u_end(1)
    CALL umPrint(umMessage,src='setdiff_4A')

    ! EW diffusion rows are contained within NS diffusion rows and row
    ! limits are combined with polar filtering rows

    !  when  nproc_y == 1  S hem filtered/diffused first then N.Hem
    ! but   diff_coeff_u  contains coeff for all rows
    IF ( nproc_y == 1) THEN
      u_end(0) = j_stop
    END IF ! nproc_y == 1

    diff_time_equ = -1.0 / LOG ( 1.0 -                            &
                          4.0 * diff_scale(global_rows - 1))
    diff_time_print = -1.0 / LOG ( 1.0 -                          &
                        4.0 * diff_scale(diff_start_row + 1))
    IF ( mype == 0 ) THEN
      WRITE(umMessage,*) ' '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*) '  ****   Horizontal diffusion is ACTIVE ****'
      CALL umPrint(umMessage,src='setdiff_4A')
      IF ( 2 * (diff_start_row/2) == diff_start_row ) THEN
        WRITE(umMessage,918) diff_start_row/2 + 1, diff_time_print
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,919) diff_time_equ
        CALL umPrint(umMessage,src='setdiff_4A')
      ELSE !  diffusion starts on v row
        WRITE(umMessage,917) diff_start_row/2 + 1, diff_time_print
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,919) diff_time_equ
        CALL umPrint(umMessage,src='setdiff_4A')
      END IF ! 2 * (diff_start_row/2) == diff_start_row
    END IF ! mype == 0

  ELSE ! L_diff_auto false, using user defined parameters


    L_test = .FALSE.
    exp1 = EXP(1.0)
    rpower = -1.0 / diff_timescale
    diff_coeff_temp = 0.25 * (1.0 - exp1**rpower)
    v_begin(0) = j_start_v
    v_end(0) = j_stop_v
    DO j = j_start_v, j_stop_v
      DO i = 1-offx, row_length
        diff_coeff_v(i,j) = diff_coeff_temp
      END DO
    END DO

    u_begin(0) = 1
    u_end(0) = rows
    DO j = 1, rows
      DO i = 1-offx, row_length
        diff_coeff_u(i,j) = diff_coeff_temp
      END DO
    END DO

    ! Set the N-S diffusion coefficent to that
    ! at the equator. Should be spatially constant
    IF (L_diff_wind_ew_only) THEN
       diff_coeff_phi=0.0
    ELSE
       diff_coeff_phi=diff_coeff_u(1,rows/2)
    END IF

    !  NS Diffusion coeff scaled from equator EW coefficient
    diff_coeff_thprint = diff_coeff_ref
    diff_coeff_uvprint = diff_coeff_ref

    IF ( L_diff_thermo) THEN
      diff_time_equ =  -1.0 / LOG( 1.0 -                          &
                                   (8.0 * diff_coeff_ref)         &
                                       ** diff_order)
      WRITE(umMessage,*) ' '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*) '  ****  Horizontal diffusion of theta ',        &
               '    is ACTIVE **** '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,903) diff_order, diff_coeff_thprint
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,904) diff_timescale
      CALL umPrint(umMessage,src='setdiff_4A')
      IF (model_type == mt_global) THEN
        WRITE(umMessage,905) polar_cap
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,906) diff_coeff_ref
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,907) diff_time_equ
        CALL umPrint(umMessage,src='setdiff_4A')
      ELSE!  model_type  /=  mt_global
        WRITE(umMessage,*)' This is the time taken to damp',              &
                 ' the 2 gridlength wave by 1/e '
        CALL umPrint(umMessage,src='setdiff_4A')
      END IF !  model_type == mt_global
    ELSE ! L_diff_thermo = .false.
      L_diff_message = .TRUE.
      IF ( L_diffusion .OR. L_cdiffusion ) THEN
        k = 1
        DO ! k = 1, model_levels but end on EXIT
          IF ( diffusion_coefficient_thermo(k) > TINY .AND.       &
                diffusion_order_thermo(k) > 0 ) THEN
            L_diff_message = .FALSE.
            EXIT
          END IF
          k = k + 1
          IF ( k > model_levels) EXIT
        END DO
      END IF ! L_diffusion .or. L_cdiffusion
      IF ( L_diff_message ) THEN
        WRITE(umMessage,*) '  '
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,*) ' *** Horizontal diffusion is  NOT ACTIVE ***'
        CALL umPrint(umMessage,src='setdiff_4A')
      END IF ! L_diff_message
    END IF ! L_diff_thermo

    IF (L_diff_wind ) THEN
      diff_time_equ =  -1.0 / LOG( 1.0 -                          &
                                      (8.0 * diff_coeff_ref)      &
                                          ** diff_order)
      WRITE(umMessage,*) ' '
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,*)'  ****   Horizontal diffusion of horizontal wind',&
                  ' is ACTIVE ****'
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,903) diff_order, diff_coeff_uvprint
      CALL umPrint(umMessage,src='setdiff_4A')
      WRITE(umMessage,901) diff_timescale
      CALL umPrint(umMessage,src='setdiff_4A')
      IF (model_type == mt_global) THEN
        WRITE(umMessage,905) polar_cap
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,906) diff_coeff_ref
        CALL umPrint(umMessage,src='setdiff_4A')
        WRITE(umMessage,907) diff_time_equ
        CALL umPrint(umMessage,src='setdiff_4A')
      ELSE !  model_type  /=  mt_global
        WRITE(umMessage,*) ' This is the time taken to damp',             &
                 ' the 2 gridlength wave by 1/e '
        CALL umPrint(umMessage,src='setdiff_4A')
      END IF ! model_type == mt_global
    ELSE ! L_diff_wind = .false.
      L_diff_message = .TRUE.
      IF ( L_diffusion .OR. L_cdiffusion ) THEN
        k = 1
        DO ! k = 1, model_levels but end on EXIT
          IF ( diffusion_coefficient_wind(k) > TINY .AND.         &
                diffusion_order_wind(k) > 0 ) THEN
            L_diff_message = .FALSE.
            EXIT
          END IF
          k = k + 1
          IF ( k > model_levels) EXIT
        END DO
      END IF ! L_diffusion .or. L_cdiffusion
      IF ( L_diff_message ) THEN
        WRITE(umMessage,*) ' *** Horizontal diffusion for winds is ',     &
                     ' NOT ACTIVE ***'
        CALL umPrint(umMessage,src='setdiff_4A')
      END IF ! L_diff_message
    END IF ! L_diff__wind

    IF ( L_diff_incs ) THEN
      WRITE(umMessage,*) 'L_diff_incs = ',L_diff_incs,                    &
              ' so increments will be horizontally diffused.'
      CALL umPrint(umMessage,src='setdiff_4A')
    END IF ! L_diff_incs
    IF ( L_diff_w ) THEN
      WRITE(umMessage,*)'L_diff_w = ',L_diff_w,                           &
             ' so w will be horizontally diffused like theta'
      CALL umPrint(umMessage,src='setdiff_4A')
    END IF ! L_diff_w

  END IF  !  L_diff_auto

  !  LAM domains can have diffusion so ...
ELSE IF (L_filter) THEN
  j_start_u = 1
  j_stop_u = rows
  j_start_v = 1
  j_stop_v = n_rows
  IF (at_extremity(PSouth)) THEN
    j_start_u = 2
    j_start_v = 2
  END IF ! at_extremity(PSouth)
  IF (at_extremity(PNorth)) THEN
    j_stop_u = rows - 1
    j_stop_v = n_rows - 1
  END IF ! at_extremity(PNorth)

  ! Now set loop bounds for one sweep of diffusion
  u_begin(0) = j_start_u
  u_end(0) = j_stop_u
  v_begin(0) = j_start_v
  v_end(0) = j_stop_v

  exp1 = EXP(1.0)
  rpower = -1.0 / diff_timescale
  diff_coeff_temp = 0.25 * (1.0 - exp1**rpower)
  diff_coeff_phi = diff_coeff_temp

  DO j = j_start_u, j_stop_u
    DO i = 1-offx, row_length
      diff_coeff_u(i,j) = diff_coeff_temp
    END DO
  END DO

  DO j = 1, n_rows
    DO i = 1-offx, row_length
      diff_coeff_v(i,j) = diff_coeff_temp
    END DO
  END DO

  WRITE(umMessage,*) ' Check filter sweeping set for non-global domains'
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*)' u_begin(0) = ', u_begin(0), ' u_end(0) = ', u_end(0)
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*)' v_begin(0) = ', v_begin(0), ' v_end(0) = ', v_end(0)
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,'(2A)')                   &
           ' Only 1st order diffusion allowed - '               &
         , ' (overwriting namelist value supplied). '
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*) ' 1st order diffusion e-folding time scale = '       &
         , diff_timescale, ' time steps'
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*) ' diffusion coefficient = ', diff_coeff_temp
  CALL umPrint(umMessage,src='setdiff_4A')

END IF ! L_filter .and. model_type  /=  mt_global

IF ( L_sponge ) THEN
  WRITE(umMessage,*) '  Sponge zone is active '
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*) '  Sides sponge width  = ', sponge_ew
  CALL umPrint(umMessage,src='setdiff_4A')
  WRITE(umMessage,*) '  Northern/Southern sponge widths = ' &
                             , sponge_ns
  CALL umPrint(umMessage,src='setdiff_4A')
  IF ( sponge_power == 1 ) THEN
    weight = 1.0 / (sponge_ew + 1)
    sponge_wts_ew(1) = 1.0 - weight
    DO i = 2, sponge_ew
      sponge_wts_ew(i) = sponge_wts_ew(i-1) - weight
    END DO
    weight = 1.0 / (sponge_ns + 1)
    sponge_wts_ns(1) = 1.0 - weight
    DO i = 2, sponge_ns
      sponge_wts_ns(i) = sponge_wts_ns(i-1) - weight
    END DO
    WRITE(umMessage,*) ' sponge_power = ', sponge_power,    &
      '  therefore sponge blending is linear'
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) ' East/West sponge weights '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) sponge_wts_ew
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) ' North/South sponge weights '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) sponge_wts_ns
    CALL umPrint(umMessage,src='setdiff_4A')
  ELSE IF ( sponge_power == 2 ) THEN
    weight = 1.0 / (sponge_ew + 1)
    sponge_wts_ew(1) = 1.0 - weight
    DO i = 2, sponge_ew
      sponge_wts_ew(i) = sponge_wts_ew(i-1) - weight
      sponge_wts_ew(i-1) = sponge_wts_ew(i-1) * sponge_wts_ew(i-1)
    END DO
    sponge_wts_ew(sponge_ew) = sponge_wts_ew(sponge_ew) *        &
                               sponge_wts_ew(sponge_ew)
    weight = 1.0 / (sponge_ns + 1)
    sponge_wts_ns(1) = 1.0 - weight
    DO i = 2, sponge_ns
      sponge_wts_ns(i) = sponge_wts_ns(i-1) - weight
      sponge_wts_ns(i-1) = sponge_wts_ns(i-1) * sponge_wts_ns(i-1)
    END DO
    sponge_wts_ns(sponge_ns) = sponge_wts_ns(sponge_ns) *        &
                               sponge_wts_ns(sponge_ns)
    WRITE(umMessage,*) ' sponge_power = ', sponge_power,    &
      '  therefore sponge blending is quadratic'
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) ' East/West sponge weights '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) sponge_wts_ew
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) ' North/South sponge weights '
    CALL umPrint(umMessage,src='setdiff_4A')
    WRITE(umMessage,*) sponge_wts_ns
    CALL umPrint(umMessage,src='setdiff_4A')
  ELSE
    WRITE(umMessage,*) ' sponge_power = ', sponge_power,    &
      '  is NOT PERMITTED'
    CALL umPrint(umMessage,src='setdiff_4A')
  END IF ! sponge_power == 1
END IF ! L_sponge

!   quick fix to fix bit-reprod problem with filtering increments
! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                               &
                 diff_coeff_u, row_length, rows, 1,             &
                 offx, offy, fld_type_u, swap_field_is_scalar)
! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                               &
                 diff_coeff_v, row_length, n_rows, 1,           &
                 offx, offy, fld_type_u, swap_field_is_scalar)

WRITE(umMessage,*) ' ENDGame diffusion settings:'
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diffusion', L_diffusion
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_cdiffusion', L_cdiffusion
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_filter', L_filter
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_auto', L_diff_auto
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pfcomb', L_pfcomb
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pfexner', L_pfexner
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pftheta', L_pftheta
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pfuv', L_pfuv
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pfw', L_pfw
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_pfincs',L_pfincs
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_thermo', L_diff_thermo
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_wind', L_diff_wind
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_wind_ew_only', L_diff_wind_ew_only
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_w', L_diff_w
CALL umPrint(umMessage,src='setdiff_4A')
WRITE(umMessage,*) '        L_diff_incs',L_diff_incs
CALL umPrint(umMessage,src='setdiff_4A')

! ----------------------------------------------------------------------
! Section 2. Print FORMATTING
! ----------------------------------------------------------------------

901  FORMAT(' Diffusion timescale for wind is ',i4,' timesteps ')
903  FORMAT(' Diffusion order is ',i2,                                 &
            ' with diffusion coefficient = ',F7.4)
904  FORMAT(' Diffusion timescale for theta is ',i4,' timesteps ')
905  FORMAT( ' This is the time taken to damp the 2 gridlength wave',  &
             ' at ',f6.2,' degrees by 1/e ')
906  FORMAT( ' The diffusion coefficient at the equator = ',ES9.2)
907  FORMAT( ' This is a damping timescale of ',f7.1,' timesteps')
910  FORMAT( '  Polar filter activates polewards of ',f6.2,' degrees')
911  FORMAT( '  and needs a maximum of ',i3,                           &
             ' pass(es) of the 1-2-1 filter')
912  FORMAT( ' polar_cap reset to ',f6.2,' degrees')
917  FORMAT('Diffusion timescale for v-row ',i4                        &
                                        ,'  is ',f6.1,' timesteps')
918  FORMAT('Diffusion timescale for u-row ',i4                        &
                                        ,'  is ',f6.1,' timesteps')
919  FORMAT('Diffusion timescale at Equator is ',f6.1,' timesteps')
920  FORMAT(' By setting polar_cap = ',f6.2,' <',f6.2,                 &
            ' (row next to pole), ')
921  FORMAT(' polar filtering will stop at row before ',f6.2)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_Setdiff
