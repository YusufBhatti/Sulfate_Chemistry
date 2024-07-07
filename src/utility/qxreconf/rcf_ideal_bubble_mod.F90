! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised bubbles

MODULE rcf_ideal_bubble_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_BUBBLE_MOD'

CONTAINS

! Subroutine rcf_ideal_bubble
!
! Description:
!   Adds a 3D temperature anomaly (bubble) to the potential temperature
!   field (warm or cold). There are a number of options for different
!   shaped bubbles. There is also the option to saturate the bubble.
!
! Method:
!   Sets up x,y,z functions containing the distance from the centre
!   of the bubble, then applies a function to generate the bubble
!   amplitude at all points, and adds this to the model field.
!   The bubble is calculated on the global grid one level at a time,
!   before being scattered to all processors to perturb local fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_bubble ( hdr_out, xi3_at_theta,                    &
                              fields_out, field_count_out )

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type
                                 
USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Locate_Mod, ONLY: &
    rcf_locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    rcf_alloc_field,           &
    rcf_dealloc_field

USE Rcf_Read_Field_Mod, ONLY: &
    rcf_read_field

USE Rcf_Write_Field_Mod, ONLY: &
    rcf_write_field

USE rcf_calc_coords_mod, ONLY: &
    rcf_calc_coords

USE rcf_gather_field_mod, ONLY: &
    rcf_gather_field_real

USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE rcf_interp_weights_mod, ONLY: &
    intw_w2rho, intw_rho2w

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_parvars, ONLY: &
    gc_all_proc_group

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_dry_rho,      &
    stashcode_thetavd,      &
    stashcode_mv,           &
    stashcode_exner_surf

USE rcf_nlist_recon_idealised_mod, ONLY: &
    idl_bubble_depth,                    &
    idl_bubble_height,                   &
    idl_bubble_max,                      &
    idl_bubble_option,                   &
    idl_bubble_width,                    &
    idl_bubble_xoffset,                  &
    idl_bubble_yoffset,                  &
    l_saturate_bubble

USE rcf_ideal_bubble_constants_mod, ONLY: &
    block,                                &
    gaussian,                             &
    idl_max_num_bubbles,                  &
    omit,                                 &
    sine,                                 &
    plume

USE idealise_run_mod,              ONLY: &
    l_ideal_2d

USE planet_constants_mod, ONLY: &
    kappa,                      &
    planet_radius,              &
    pref,                       &
    r,                          &
    recip_kappa,                &
    recip_epsilon

USE lam_config_inputs_mod, ONLY: &
    delta_lat,                   &
    delta_lon

USE conversions_mod, ONLY: &
    pi

USE missing_data_mod, ONLY: &
    imdi

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE ereport_mod, ONLY: &
    ereport

USE um_parcore, ONLY: &
    nproc

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE(um_header_type), INTENT(IN) :: hdr_out      ! Output dump header
REAL, INTENT(IN)          :: xi3_at_theta(output_grid % loc_p_field,         &
                                        0:output_grid % model_levels+1)
                                          ! heights at theta levels
TYPE(field_type), POINTER :: fields_out(:)       ! output fields
INTEGER, INTENT(IN)       :: field_count_out     ! no. of output fields

! Local variables
INTEGER :: i
INTEGER :: j
INTEGER :: k
INTEGER :: bubble
INTEGER :: icode
INTEGER :: nlen

REAL    :: xi3_at_theta_global(output_grid % glob_p_row_length,   &
                               output_grid % glob_p_rows) ! A single level
REAL    :: xi3_at_theta_global_ml(output_grid % glob_p_row_length,   &
                                  output_grid % glob_p_rows) ! At model_levels

REAL    :: xi1_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi2_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi1_u_global(output_grid % glob_u_row_length, &
                        output_grid % glob_u_rows)
REAL    :: xi2_v_global(output_grid % glob_v_row_length, &
                        output_grid % glob_v_rows)

! Bubble calculations:
REAL    :: glob_bubble_perturbation(output_grid % glob_p_row_length, &
                                    output_grid % glob_p_rows)
                                   ! 2d single level
REAL    :: loc_bubble_perturbation(output_grid % loc_p_field, &
                                   0:output_grid % model_levels)
                                   ! 1d per level, multiple levels

REAL    :: x0   ! centre coordinate
REAL    :: y0   ! centre coordinate
REAL    :: z0   ! centre coordinate
REAL    :: dist ! distance from bubble centre
REAL    :: xi   ! temp variable
REAL    :: bubble_fn ! function for calculating bubble perturbation
REAL    :: t
REAL    :: rh   ! relative humidity
REAL    :: hght ! height above surface

REAL, ALLOCATABLE  :: exner_theta(:,:)  ! Exner on theta levels
REAL, ALLOCATABLE  :: p_on_theta(:,:)   ! pressure on theta levels  (Pa)
REAL, ALLOCATABLE  :: t_on_theta(:,:)   ! temperature on theta levels (K)
REAL, ALLOCATABLE  :: mv_sat(:,:)       ! saturated mixing ratio (kg/kg)

CHARACTER(LEN=errormessagelength) :: cmessage

! Indices for locating fields
INTEGER                            :: pos_thetavd
INTEGER                            :: pos_dry_rho
INTEGER                            :: pos_exner
INTEGER                            :: pos_mv
INTEGER                            :: pos_exner_surf
INTEGER                            :: pe
INTEGER                            :: iend
INTEGER                            :: div
INTEGER                            :: blk
INTEGER                            :: rlen

LOGICAL                            :: l_plume_bubble  ! different saturation
                                                      ! option

! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: dry_rho_out
TYPE( field_type ), POINTER        :: mv_out
TYPE( field_type ), POINTER        :: exner_surf_out

CHARACTER (LEN=*),  PARAMETER :: routinename='RCF_IDEAL_BUBBLE'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Extract some fields
CALL rcf_locate( stashcode_prog_sec, stashcode_thetavd,         &
               fields_out, field_count_out, pos_thetavd)
thetavd_out => fields_out(pos_thetavd)
CALL rcf_alloc_field( thetavd_out )
CALL rcf_read_field( thetavd_out, hdr_out, decomp_rcf_output )

CALL rcf_locate( stashcode_prog_sec, stashcode_dry_rho,             &
               fields_out, field_count_out, pos_dry_rho)
dry_rho_out => fields_out(pos_dry_rho)
CALL rcf_alloc_field( dry_rho_out )
! No read of current state required; field is completely overwritten.

CALL rcf_locate( stashcode_prog_sec, stashcode_exner,           &
               fields_out, field_count_out, pos_exner)
exner_out => fields_out(pos_exner)
CALL rcf_alloc_field( exner_out )
CALL rcf_read_field( exner_out, hdr_out, decomp_rcf_output )

! Calculate global grid info and scatter to generate local data
CALL rcf_calc_coords(hdr_out, output_grid, xi1_p_global, xi2_p_global, &
                     xi1_u_global, xi2_v_global)

IF (l_saturate_bubble) THEN
  l_plume_bubble = .FALSE.   ! only used for saturated bubbles

  ! Then applying perturbation to theta and not theta_vd

  ! Need water vapour field and exner_surface
  CALL rcf_locate( stashcode_prog_sec, stashcode_mv,                 &
               fields_out, field_count_out, pos_mv)
  mv_out => fields_out(pos_mv)
  CALL rcf_alloc_field( mv_out )
  CALL rcf_read_field( mv_out, hdr_out, decomp_rcf_output )
  CALL rcf_locate( stashcode_prog_sec, stashcode_exner_surf,         &
               fields_out, field_count_out, pos_exner_surf)
  exner_surf_out => fields_out(pos_exner_surf)
  CALL rcf_alloc_field( exner_surf_out )
  CALL rcf_read_field( exner_surf_out, hdr_out, decomp_rcf_output )

  ! Need theta instead of theta_vd
  DO k = 1, output_grid % model_levels+1
    DO i = 1, output_grid % loc_p_field
      thetavd_out % data(i,k) = thetavd_out % data(i,k)        &
                       /(1.0 + recip_epsilon * mv_out % data(i,k))
    END DO
  END DO 
END IF


pe = 0
! Sinusoidal bubbles need to know xi3_theta_global at model_levels, so:
IF (ANY(idl_bubble_option == sine)) THEN
  ! Gather on PE0, broadcast to everything else.
  CALL rcf_gather_field_real(xi3_at_theta(:,output_grid % model_levels), &
                             xi3_at_theta_global_ml,                     &
                             output_grid % loc_p_row_length,             &
                             output_grid % loc_p_rows,                   &
                             output_grid % glob_p_row_length,            &
                             output_grid % glob_p_rows, pe,              &
                             gc_all_proc_group )

  rlen = output_grid % glob_p_row_length * output_grid % glob_p_rows
  CALL gc_rbcast(1, rlen, pe, nproc, icode,  xi3_at_theta_global_ml)
END IF


! We are only concerned with 0:model_levels of xi3_at_theta:
div = (output_grid % model_levels + 1) / nproc
IF ( nproc * div  <   output_grid % model_levels + 1) div = div + 1

DO blk = 1, div

  IF (blk == div ) THEN
    iend = output_grid % model_levels
  ELSE
    iend =  blk * nproc - 1
  END IF


  DO k = ((blk-1) * nproc ), iend


    CALL rcf_gather_field_real(xi3_at_theta(:,k),                   &
                               xi3_at_theta_global,                 &
                               output_grid % loc_p_row_length,      &
                               output_grid % loc_p_rows,            &
                               output_grid % glob_p_row_length,     &
                               output_grid % glob_p_rows, pe,       &
                               gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO


    ! Calculate the bubble perturbation for the global grid on this level:

    ! Begin bubble initialisation
    DO j = 1, output_grid % glob_p_rows
      DO i = 1, output_grid % glob_p_row_length
        glob_bubble_perturbation(i,j) = 0.0
      END DO
    END DO

    ! Loop over bubbles
    DO bubble = 1, idl_max_num_bubbles

      IF (idl_bubble_option(bubble) == imdi .OR. &
          idl_bubble_option(bubble) == omit) CYCLE

      ! Set coordinates of the bubble's centre.
      ! u- and v-grids for ENDGame are offset by 1 in the reconfiguration, so
      ! that e.g. xi1_u(0,1) in the model is xi1_u(1,1) in the reconfiguration.
      ! This means that what the reconfiguration calls xi1_u(M,N) is actually
      ! xi1_u(M-1,N), and so delta_lon must be added to the final point to
      ! calculate the centre correctly. Similarly for xi2_v.
      x0 = xi1_u_global(1,1) + idl_bubble_xoffset(bubble) *  &
       (xi1_u_global(output_grid % glob_u_row_length,output_grid % glob_u_rows)&
        + delta_lon - xi1_u_global(1,1))

      y0 = xi2_v_global(1,1) + idl_bubble_yoffset(bubble) *  &
       (xi2_v_global(output_grid % glob_v_row_length,output_grid % glob_v_rows)&
        + delta_lat - xi2_v_global(1,1))

      z0 = idl_bubble_height(bubble)

      DO j = 1, output_grid % glob_p_rows
        DO i = 1, output_grid % glob_p_row_length

          dist = SQRT( (xi1_p_global(i,j) - x0)**2  &
                     + (xi2_p_global(i,j) - y0)**2  &
                     + (xi3_at_theta_global(i,j) - planet_radius - z0)**2 )

          SELECT CASE(idl_bubble_option(bubble))

          CASE(gaussian)

            xi = (dist - idl_bubble_depth(bubble)) / idl_bubble_width(bubble)

            IF (xi <= 1.0) THEN
              bubble_fn = 1.0
            ELSE
              bubble_fn = EXP(-xi**2)
            END IF

          CASE(plume)   ! No height dependence

            l_plume_bubble = .TRUE.

            hght = xi3_at_theta_global(i,j)-planet_radius 
            IF (hght <= idl_bubble_depth(bubble)) THEN
              IF (l_ideal_2d) THEN
                dist = SQRT( (xi1_p_global(i,j) - x0)**2 )
                xi = dist / idl_bubble_width(bubble)
                bubble_fn=EXP(-xi**2)
              ELSE
                dist = SQRT( (xi1_p_global(i,j) - x0)**2  &
                           + (xi2_p_global(i,j) - y0)**2 )

                xi = dist / idl_bubble_width(bubble)
                bubble_fn=EXP(-xi**2)
              END IF
            ELSE
              bubble_fn = 0.0
            END IF

          CASE(block)

            IF (ABS(xi1_p_global(i,j) - x0) <= 0.5 * idl_bubble_width(bubble)  &
            .AND. ABS(xi2_p_global(i,j) - y0) <= 0.5 * idl_bubble_width(bubble)&
            .AND. ABS(xi3_at_theta_global(i,j) - planet_radius - z0)           &
                                     <= 0.5 * idl_bubble_depth(bubble)) THEN
              bubble_fn = 1.0
            ELSE
              bubble_fn = 0.0
            END IF

          CASE(sine)

            dist = 1.0                                                        &
                 + ((xi1_p_global(i,j) - x0) / idl_bubble_width(bubble))**2   &
                 + ((xi2_p_global(i,j) - y0) / idl_bubble_depth(bubble))**2
            z0   = (xi3_at_theta_global(i,j) - planet_radius) /               &
                   (xi3_at_theta_global_ml(i,j) - planet_radius)

            bubble_fn = SIN(pi*z0)/dist

          CASE DEFAULT

            icode = 10
            WRITE(cmessage, '(A, I0)') &
              'Invalid idl_bubble_option:',idl_bubble_option
            CALL ereport('rcf_ideal_bubble', icode, cmessage)

          END SELECT

          glob_bubble_perturbation(i,j) = glob_bubble_perturbation(i,j)       &
                                          + idl_bubble_max(bubble)*bubble_fn

        END DO ! i
      END DO  ! j

    END DO  ! bubble no.

! Scatter the perturbation back to all PEs:
  pe = 0
  DO k = ((blk-1) * nproc ), iend
    CALL rcf_scatter_field_real(loc_bubble_perturbation(:,k),        &
                                glob_bubble_perturbation(:,:),       &
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO ! k
END DO ! blk


! Perturb fields.
! Perturbation array defined over 0:model_levels, data array 1:model_levels+1
DO k = 1, output_grid % model_levels+1
  DO i = 1, output_grid % loc_p_field
    thetavd_out % data(i,k) = thetavd_out % data(i,k)                        &
                            + loc_bubble_perturbation(i,k-1)
  END DO
END DO


! saturated bubbles where theta perturbation >0.1K
IF (l_saturate_bubble) THEN

  ! Need pressure and temperature.
  ! Calculations ignoring top model level - assumption that a bubble
  ! will not include that level. Don't have exner for top level.
  ALLOCATE(p_on_theta(output_grid % loc_p_field,output_grid % model_levels))
  ALLOCATE(exner_theta(output_grid % loc_p_field,output_grid % model_levels))
  ALLOCATE(t_on_theta(output_grid % loc_p_field, output_grid % model_levels))
  ALLOCATE(mv_sat(output_grid % loc_p_field, output_grid % model_levels))
  ! Surface
  DO i = 1, output_grid % loc_p_field
    p_on_theta(i,1) = pref*exner_surf_out % data(i,1)**recip_kappa
    exner_theta(i,1) = exner_surf_out % data(i,1)
    t_on_theta(i,1) = thetavd_out % data(i,1)*exner_surf_out % data(i,1)
  END DO
  ! Model theta levels but not top 
  DO k = 1, output_grid % model_levels-1
    DO i = 1, output_grid % loc_p_field
      exner_theta(i,k+1) = intw_rho2w(k,1) * exner_out % data(i,k+1)         &
                           + intw_rho2w(k,2) * exner_out % data(i,k)
      p_on_theta(i,k+1) = pref *exner_theta(i,k+1) **recip_kappa
      t_on_theta(i,k+1) = thetavd_out % data(i,k+1) * exner_theta(i,k+1) 
    END DO
  END DO

  nlen = (output_grid % model_levels)*output_grid % loc_p_field
  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(mv_sat,t_on_theta, p_on_theta, nlen ,.TRUE. )
 
  IF (l_plume_bubble) THEN
    ! Specail case to reproduce the type of plume initialisation done by a
    ! branch to new dynamics as a completely saturated bubble can generate
    ! a ring of convection around the edge of the bubble in the horizontal.
    ! Assumption that all bubbles are type plume and have the same magnitude
    ! For bubble region go from RH =80% to RH=100% in the centre
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_p_field
        IF (loc_bubble_perturbation(i,k-1) > 0.01) THEN
          rh = 0.8 + 0.2*loc_bubble_perturbation(i,k-1)/idl_bubble_max(1)
          mv_out % data(i,k) = rh* mv_sat(i,k)
        END IF
      END DO
    END DO
  ELSE   ! standard saturation of all bubbles to 100% relative humidity
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_p_field
        IF (loc_bubble_perturbation(i,k-1) > 0.1) THEN
          mv_out % data(i,k) = mv_sat(i,k)
        END IF
      END DO
    END DO
  END IF
  DO k = 1, output_grid % model_levels+1
    DO i = 1, output_grid % loc_p_field
      ! reset virtual dry potential temperature using modifed mv field
      thetavd_out % data(i,k) = thetavd_out % data(i,k) *            &
                                (1.0+recip_epsilon * mv_out % data(i,k))
    END DO
  END DO

  ! write out new mv field
  CALL rcf_write_field(fields_out(pos_mv), hdr_out, decomp_rcf_output)
  ! Deallocate space 
  CALL rcf_dealloc_field(fields_out(pos_mv))
  CALL rcf_dealloc_field(fields_out(pos_exner_surf))

  DEALLOCATE(mv_sat)
  DEALLOCATE(t_on_theta)
  DEALLOCATE(exner_theta)
  DEALLOCATE(p_on_theta)

END IF

! Correct dry rho:
! (thetavd_out declared over 1:model_levels+1, not 0:model_levels, so +1)
DO k = 1, output_grid % model_levels
  DO i = 1, output_grid % loc_p_field
    t = (intw_w2rho(k,1) * thetavd_out % data(i,k+1) +                       &
         intw_w2rho(k,2) * thetavd_out % data(i,k)) * exner_out % data(i,k)
    dry_rho_out % data(i,k) = pref * exner_out % data(i,k)**(1.0/kappa)/(r*t)
  END DO
END DO

! Write out and deallocate the modified fields
CALL rcf_write_field   (fields_out(pos_thetavd),    hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_thetavd))
CALL rcf_write_field   (fields_out(pos_dry_rho),    hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_dry_rho))
CALL rcf_dealloc_field (fields_out(pos_exner))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_bubble
END MODULE rcf_ideal_bubble_mod
