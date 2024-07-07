! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE r2_cloud_level_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'R2_CLOUD_LEVEL_DIAG_MOD'
CONTAINS

SUBROUTINE r2_cloud_level_diag(control, dimen, atm, cld, cdnc                  &
   , nclds, i_gather, i_cloud, l_all_temps, l_radius                           &
   , weighted_diag, sum_weight_diag                                            &
   , col_list, row_list, row_length, rows, nd_field)

! Subroutine to calculate an observed effective radius.
!
! Purpose:
!   Calculates a 2-D cloud-top diagnostic from the input 3-D version
!   of the diagnostic.
!
! Method:
!   For each type of cloud containing water in any layer the input diagnostic
!   is weighted with the product of the area of the cloud and the probability
!   that light emitted from the cloud reaches the observing instrument.
!

USE rad_pcf
USE def_dimen,   ONLY: StrDim
USE def_control, ONLY: StrCtrl
USE def_atm,     ONLY: StrAtm
USE def_cld,     ONLY: StrCld
USE water_constants_mod, ONLY: tm
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength


IMPLICIT NONE


! Controlling options:
TYPE (StrCtrl), INTENT(IN) :: control

! Dimensions:
TYPE (StrDim), INTENT(IN)  :: dimen

! Atmospheric properties:
TYPE(StrAtm), INTENT(IN)   :: atm

! Cloud properties:
TYPE(StrCld), INTENT(IN)   :: cld
REAL, INTENT(IN) :: cdnc(dimen%nd_profile, dimen%id_cloud_top : dimen%nd_layer,&
                         dimen%nd_cloud_component)
! Cloud-droplet number concentration (m-3)

! Dimensions of arrays:
INTEGER, INTENT(IN) :: row_length
!   Number of grid-points in EW-direction in the local domain
INTEGER, INTENT(IN) :: rows
!   Number of grid-points in NS-direction in the local domain
INTEGER, INTENT(IN) :: nd_field
!   Size of arrays passed from main code

! Actual sizes used:
INTEGER, INTENT(IN) :: nclds
!   Number of cloudy levels
INTEGER, INTENT(IN) :: i_gather(nd_field)
!   List of gathered points
INTEGER, INTENT(IN) :: col_list(nd_field)
!   EW indices of gathered points in the 2-D domain
INTEGER, INTENT(IN) :: row_list(nd_field)
!   NS indices of gathered points in the 2-D domain

! Logical flags for diagnostics
LOGICAL, INTENT(IN) :: l_all_temps
!   If TRUE, the routine has been called to obtain diagnostics
!   for clouds consisting of liquid water at any temperature
!   (as done in MODIS retrievals). If FALSE, only clouds with
!   temperatures above freezing are to be diagnosed (as done
!   in AVHRR retrievals).
LOGICAL, INTENT(IN) :: l_radius
!   If TRUE, the routine has been called for calculating cloud
!   droplet effective radius at cloud top. If FALSE then it's
!   been called for calculating cloud droplet number concentration
!   at cloud top.

! Representation of clouds
INTEGER, INTENT(IN) :: i_cloud
!   Treatment of overlaps

REAL, INTENT(INOUT) :: weighted_diag(row_length, rows)
!   Weighted sum of cloud-top diagnostic and weighting function
REAL, INTENT(INOUT) :: sum_weight_diag(row_length, rows)
!   Sum of weights for cloud-top diagnostic



! Local variables:
INTEGER :: ierr
!   Error flag

INTEGER :: i, l
!   Loop variables
INTEGER :: i_inv
!   Inverted loop index

REAL ::                                                                        &
  trans_overlying_space(dimen%nd_profile),                                     &
!   Probability of a photon in clear air in the level above
!   the current one reaching space
  area_exposed(dimen%nd_profile),                                              &
!   Total area of cloud in the current layer exposed to
!   clear air in the layer above
  area_exposed_st(dimen%nd_profile),                                           &
!   Total area of stratiform cloud in the current layer
!   exposed to clear air in the layer above
  area_exposed_cnv(dimen%nd_profile),                                          &
!   Total area of convective cloud in the current layer
!   exposed to clear air in the layer above
  area_clear_above(dimen%nd_profile),                                          &
!   Area of the clear sky region in the layer above
  area_strat(dimen%nd_profile),                                                &
!   Area of stratiform cloud in the current layer
  area_strat_above(dimen%nd_profile),                                          &
!   Area of stratiform cloud in the layer above
  area_conv(dimen%nd_profile),                                                 &
!   Area of convective cloud in the current layer
  area_conv_above(dimen%nd_profile),                                           &
!   Area of convective cloud in the layer above
  area_clear_clear(dimen%nd_profile),                                          &
!   Area of boundary where clear sky overlies clear sky
  area_clear(dimen%nd_profile),                                                &
!   Area of clear sky in the current layer down to a level
  area_uncorrelated(dimen%nd_profile),                                         &
!   Uncorrelated region on the interface
  weighted_diag_g(dimen%nd_profile),                                           &
!   Weighted sum of cloud-top diagnostic and weighting function
  sum_weight_diag_g(dimen%nd_profile)
!   Sum of weights for cloud-top diagnostic

! Variables for gathering
INTEGER :: l_list(dimen%nd_profile)
!   Indices of points in list

! Indicator function
REAL :: chi_cnv(dimen%nd_profile)
!   Convective indicator function
REAL :: chi_st(dimen%nd_profile)
!   Stratiform indicator function

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'R2_CLOUD_LEVEL_DIAG'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr=i_normal

! Initialization of local fields.
DO l=1, atm%n_profile
  weighted_diag_g(l)=0.0e+00
  sum_weight_diag_g(l)=0.0e+00
END DO

! Initialize the transmision above clouds.
DO l=1, atm%n_profile
  trans_overlying_space(l)=1.0e+00
  area_clear_above(l)=1.0e+00
END DO
IF (i_cloud == ip_cloud_triple) THEN
  DO l=1, atm%n_profile
    area_strat_above(l)=0.0e+00
    area_conv_above(l)=0.0e+00
  END DO
END IF

! Step down through the atmosphere calculating contributions to
! the diagnostics and subsequently allowing for transmission
! through the current layer.

DO i=nclds, 1, -1
  i_inv=atm%n_layer+1-i

  DO l=1, atm%n_profile
    area_clear(l)=1.0e+00-cld%w_cloud(l, i_inv)
  END DO

  ! Calculate the local area of cloud radiating into clear air.
  IF (i_cloud == ip_cloud_mix_random) THEN
    DO l=1, atm%n_profile
      area_exposed(l)=cld%w_cloud(l, i_inv)                                    &
         *area_clear_above(l)
    END DO
  ELSE IF ( (i_cloud == ip_cloud_mix_max) .OR.                                 &
            (i_cloud == ip_cloud_triple) ) THEN
    DO l=1, atm%n_profile
      area_exposed(l)=MAX(0.0e+00, (cld%w_cloud(l, i_inv)                      &
         +area_clear_above(l)-1.0e+00))
    END DO
  END IF


  IF (control%i_cloud_representation == ip_cloud_conv_strat) THEN

    IF ( (i_cloud == ip_cloud_mix_max) .OR.                                    &
         (i_cloud == ip_cloud_mix_random) ) THEN

      ! If the overlap of convective cloud is not assumed
      ! to be coherent the overall exposed area may be
      ! partitioned according to the fractional
      ! contributions of cloud in the current layer.

      DO l=1, atm%n_profile
        area_exposed_st(l)=area_exposed(l)                                     &
           *cld%frac_cloud(l, i_inv, ip_cloud_type_strat)
        area_exposed_cnv(l)=area_exposed(l)                                    &
           *cld%frac_cloud(l, i_inv, ip_cloud_type_conv)
      END DO

    ELSE IF (i_cloud == ip_cloud_triple) THEN

      ! Here, the different types of clouds overlap
      ! coherently so stratiform cloud will be exposed
      ! only if there is less stratiform cloud in the
      ! layer above and more clear air in the layer above:
      ! under these conditions the non-correlated areas
      ! overlap randomly.

      DO l=1, atm%n_profile
        area_strat(l)=cld%w_cloud(l, i_inv)                                    &
           *cld%frac_cloud(l, i_inv, ip_cloud_type_strat)
        area_conv(l)=cld%w_cloud(l, i_inv)                                     &
          *cld%frac_cloud(l, i_inv, ip_cloud_type_conv)
        area_uncorrelated(l)=1.0e+00                                           &
          -MIN(area_clear(l), area_clear_above(l))                             &
          -MIN(area_strat(l), area_strat_above(l))                             &
          -MIN(area_conv(l), area_conv_above(l))
        ! First find the area of uncorrelated stratiform cloud.
        area_exposed_st(l)=MAX(0.0e+00                                         &
           , (area_strat(l)-area_strat_above(l)))
        area_exposed_st(l)=MAX(0.0e+00, area_exposed_st(l)                     &
           *(area_clear_above(l)-area_clear(l)))
        ! Now normalize within the uncorrelated region.
        ! If the uncorrelated area is 0 the exposed area
        ! must be 0, so no second branch of the IF-test
        ! is required.
        IF (area_uncorrelated(l) >  0.0e+00)                                   &
          area_exposed_st(l)                                                   &
            =area_exposed_st(l)/area_uncorrelated(l)
        area_exposed_cnv(l)                                                    &
           =area_exposed(l)-area_exposed_st(l)
      END DO
    ELSE
      cmessage =                                                               &
         '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '                   &
         //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
      ierr=i_err_fatal
      GO TO 9999
    END IF

    ! The indicator functions for liquid water in
    ! convective or straiform clouds are set to 1
    ! if there is any liquid water and to 0 otherwise.
    DO l=1, atm%n_profile
      IF (cld%condensed_mix_ratio(l, i_inv, ip_clcmp_cnv_water)                &
          >  0.0e+00) THEN
        chi_cnv(l)=1.0e+00
      ELSE
        chi_cnv(l)=0.0e+00
      END IF
      IF (cld%condensed_mix_ratio(l, i_inv, ip_clcmp_st_water)                 &
          >  0.0e+00) THEN
        chi_st(l)=1.0e+00
      ELSE
        chi_st(l)=0.0e+00
      END IF
    END DO

    ! Include contributions from convective and stratiform water clouds.
    DO l=1, atm%n_profile
     IF (l_radius) THEN
      weighted_diag_g(l)=weighted_diag_g(l)                                    &
         +trans_overlying_space(l)                                             &
         *(area_exposed_cnv(l)*chi_cnv(l)                                      &
         *cld%condensed_dim_char(l, i_inv, ip_clcmp_cnv_water)                 &
         +area_exposed_st(l)*chi_st(l)                                         &
         *cld%condensed_dim_char(l, i_inv, ip_clcmp_st_water))
     ELSE
      weighted_diag_g(l)=weighted_diag_g(l)                                    &
         +trans_overlying_space(l)                                             &
         *(area_exposed_cnv(l)*chi_cnv(l)                                      &
         *cdnc(l, i_inv, ip_clcmp_cnv_water)                                   &
         +area_exposed_st(l)*chi_st(l)                                         &
         *cdnc(l, i_inv, ip_clcmp_st_water))
     END IF
      sum_weight_diag_g(l)=sum_weight_diag_g(l)                                &
         +trans_overlying_space(l)                                             &
         *(area_exposed_cnv(l)*chi_cnv(l)                                      &
         +area_exposed_st(l)*chi_st(l))
    END DO

  ELSE IF ((control%i_cloud_representation == ip_cloud_csiw) .OR.              &
          (control%i_cloud_representation == ip_cloud_ice_water)) THEN

    IF ( (i_cloud == ip_cloud_mix_max) .OR.                                    &
         (i_cloud == ip_cloud_mix_random) ) THEN

      ! If the overlap of convective cloud is not assumed
      ! to be coherent the overall exposed area may be
      ! partitioned according to the fractional
      ! contributions of cloud in the current layer.
      ! the exposed areas include only the parts of the
      ! clouds containing water droplets.

      DO l=1, atm%n_profile
        area_exposed_st(l)=area_exposed(l)                                     &
           *cld%frac_cloud(l, i_inv, ip_cloud_type_sw)
        area_exposed_cnv(l)=area_exposed(l)                                    &
           *cld%frac_cloud(l, i_inv, ip_cloud_type_cw)
      END DO

    ELSE IF (i_cloud == ip_cloud_triple) THEN

      ! Here, the different types of clouds overlap
      ! coherently so stratiform cloud will be exposed
      ! only if there is less stratiform cloud in the
      ! layer above and more clear air in the layer above:
      ! under these conditions the non-correlated areas
      ! overlap randomly.
      ! The actual exposed areas of convective or
      ! stratiform cloud must then be weighted by factors
      ! representing the liquid portion of each cloud, since
      ! nothing is retrieved over ice. (The horizontal
      ! arrangement of ice and water within either type of
      ! cloud is random).

      DO l=1, atm%n_profile

        area_strat(l)=cld%w_cloud(l, i_inv)                                    &
           *(cld%frac_cloud(l, i_inv, ip_cloud_type_sw)                        &
           +cld%frac_cloud(l, i_inv, ip_cloud_type_si))
        area_conv(l)=cld%w_cloud(l, i_inv)                                     &
           *(cld%frac_cloud(l, i_inv, ip_cloud_type_cw)                        &
           +cld%frac_cloud(l, i_inv, ip_cloud_type_ci))
        area_uncorrelated(l)=1.0e+00                                           &
           -MIN(area_clear(l), area_clear_above(l))                            &
           -MIN(area_strat(l), area_strat_above(l))                            &
           -MIN(area_conv(l), area_conv_above(l))
        area_exposed_st(l)=MAX(0.0e+00                                         &
           , (area_strat(l)-area_strat_above(l)))
        IF (area_uncorrelated(l) >  0.0e+00) THEN
          area_exposed_st(l)                                                   &
             =MAX(0.0e+00, area_exposed_st(l)                                  &
             *(area_clear_above(l)-area_clear(l)))                             &
             /area_uncorrelated(l)
        ELSE
          area_exposed_st(l)=0.0e+00
        END IF
        area_exposed_cnv(l)                                                    &
           =area_exposed(l)-area_exposed_st(l)

        IF (cld%frac_cloud(l, i_inv, ip_cloud_type_cw)                         &
            >  0.0e+00) THEN
          area_exposed_cnv(l)=area_exposed_cnv(l)                              &
             /(1.0e+00                                                         &
             +cld%frac_cloud(l, i_inv, ip_cloud_type_ci)                       &
             /cld%frac_cloud(l, i_inv, ip_cloud_type_cw))
        ELSE
          area_exposed_cnv(l)=0.0e+00
        END IF

        IF (cld%frac_cloud(l, i_inv, ip_cloud_type_sw)                         &
            >  0.0e+00) THEN
          area_exposed_st(l)=area_exposed_st(l)                                &
             /(1.0e+00                                                         &
             +cld%frac_cloud(l, i_inv, ip_cloud_type_si)                       &
             /cld%frac_cloud(l, i_inv, ip_cloud_type_sw))
        ELSE
          area_exposed_st(l)=0.0e+00
        END IF

      END DO
    ELSE
      cmessage =                                                               &
         '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '                   &
         //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
      ierr=i_err_fatal
      GO TO 9999
    END IF


    DO l=1, atm%n_profile

      IF ((atm%t(l, i_inv)  >   tm) .OR. l_all_temps) THEN
       IF (l_radius) THEN
        weighted_diag_g(l)=weighted_diag_g(l)                                  &
           +trans_overlying_space(l)                                           &
           *(area_exposed_cnv(l)                                               &
           *cld%condensed_dim_char(l, i_inv, ip_clcmp_cnv_water)               &
           +area_exposed_st(l)                                                 &
           *cld%condensed_dim_char(l, i_inv, ip_clcmp_st_water))
       ELSE
        weighted_diag_g(l)=weighted_diag_g(l)                                  &
           +trans_overlying_space(l)                                           &
           *(area_exposed_cnv(l)                                               &
           *cdnc(l, i_inv, ip_clcmp_cnv_water)                                 &
           +area_exposed_st(l)                                                 &
           *cdnc(l, i_inv, ip_clcmp_st_water))
       END IF
        sum_weight_diag_g(l)=sum_weight_diag_g(l)                              &
           +trans_overlying_space(l)                                           &
           *(area_exposed_cnv(l)+area_exposed_st(l))
      END IF
    END DO

  END IF


  ! Advance the stored quantities refferring to overlying layers.

  ! The transmission to space currently holds the probability that
  ! a photon travelling upwards in the clear air in the layer above
  ! will escape to space without encountering a cloud. To advance
  ! this to the current layer it must be multiplied by a factor
  ! representing the overlap assumption at the top of the present
  ! layer.

  IF (i_cloud == ip_cloud_mix_random) THEN

    DO l=1, atm%n_profile
      trans_overlying_space(l)=trans_overlying_space(l)                        &
         *area_clear_above(l)
    END DO

  ELSE IF ( (i_cloud == ip_cloud_mix_max) .OR.                                 &
            (i_cloud == ip_cloud_triple) ) THEN

    DO l=1, atm%n_profile
      area_clear_clear(l)=MIN(area_clear(l)                                    &
         , area_clear_above(l))
      IF (area_clear(l) >  0.0e+00) THEN
        trans_overlying_space(l)=trans_overlying_space(l)                      &
           *area_clear_clear(l)/area_clear(l)
      ELSE
        trans_overlying_space(l)=0.0e+00
      END IF
    END DO

  END IF

  ! Advance the areas of cloud.
  DO l=1, atm%n_profile
    area_clear_above(l)=area_clear(l)
  END DO
  IF (control%i_cloud_representation == ip_cloud_conv_strat) THEN
    DO l=1, atm%n_profile
      area_strat_above(l)=cld%w_cloud(l, i_inv)                                &
         *cld%frac_cloud(l, i_inv, ip_cloud_type_strat)
    END DO
  ELSE IF ((control%i_cloud_representation == ip_cloud_csiw) .OR.              &
           (control%i_cloud_representation == ip_cloud_ice_water)) THEN
    DO l=1, atm%n_profile
      area_strat_above(l)=area_strat(l)
      area_conv_above(l)=area_conv(l)
    END DO
  END IF

END DO



! Scatter the diagnostics back to the output arrays and convert
! units as appropriate (effective radius to microns, CDNC to cm-3)
! to avoid fields being corrupted by packing.
DO l=1, atm%n_profile
 IF (l_radius) THEN
  weighted_diag(col_list(l), row_list(l))                                      &
    =1.0e06*weighted_diag_g(l)
 ELSE
  weighted_diag(col_list(l), row_list(l))                                      &
    =1.0e-06*weighted_diag_g(l)
 END IF
  sum_weight_diag(col_list(l), row_list(l))                                    &
    =sum_weight_diag_g(l)
END DO


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE r2_cloud_level_diag
END MODULE r2_cloud_level_diag_mod
