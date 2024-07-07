! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised vapour mixing ratio calculations (3d)

MODULE rcf_ideal_vapour_calc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_VAPOUR_CALC_MOD'

CONTAINS

! Subroutine rcf_ideal_vapour_calc
!
! Description:
!   Calculates vapour mixing ratios for selected qprofile numbers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_vapour_calc ( hdr_out, xi3_at_theta, fields_out,        &
                                   field_count_out)

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

USE rcf_interp_weights_mod, ONLY: &
    intw_rho2w

USE decomp_params, ONLY: &
    decomp_rcf_output

USE rcf_nlist_recon_idealised_mod, ONLY: &
    qprofile_number

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_mv,           &
    stashcode_thetavd

USE planet_constants_mod, ONLY: &
    p_zero,                     &
    planet_radius,              &
    recip_kappa

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    PrintStatus,      &
    PrStatus_Diag,    &
    PrStatus_Normal

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE ereport_mod, ONLY: &
    ereport

USE rcf_ideal_qprofile_mod, ONLY: &
    qp_BS1999,                    &
    qp_qsat

USE mmr_BS1999_mod, ONLY: &
    calc_gas_mixing_ratio_BS1999

USE gas_list_pcf, ONLY: &
    ip_h2o

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
INTEGER :: k
INTEGER :: ErrorStatus

! Temporary variables for idealised moisture calculations
REAL :: sat_z0
REAL :: sat_z1
REAL :: sat_pc
REAL :: sat_min
REAL :: sat_grad
REAL :: t_tmp      ! T for a single grid point
REAL :: p_tmp      ! p for a single grid point

REAL, ALLOCATABLE :: p_on_theta(:,:)    ! p on theta levels
REAL, ALLOCATABLE :: t_on_theta(:,:)    ! t on theta levels
REAL, ALLOCATABLE :: mv_BS1999(:,:)     ! Mixing ratio for gas giants

CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_IDEAL_VAPOUR_CALC'

! Indices for locating fields
INTEGER                            :: pos_exner
INTEGER                            :: pos_exner_surf
INTEGER                            :: pos_thetavd
INTEGER                            :: pos_mv

! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_surf_out
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: mv_out

INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Extract (and read) some fields
CALL rcf_locate( stashcode_prog_sec, stashcode_exner,           &
               fields_out, field_count_out, pos_exner)
exner_out => fields_out(pos_exner)
CALL rcf_alloc_field( exner_out )
CALL rcf_read_field(exner_out, hdr_out, decomp_rcf_output)

CALL rcf_locate( stashcode_prog_sec, stashcode_thetavd,         &
               fields_out, field_count_out, pos_thetavd)
thetavd_out => fields_out(pos_thetavd)
CALL rcf_alloc_field( thetavd_out )
CALL rcf_read_field(thetavd_out, hdr_out, decomp_rcf_output)

CALL rcf_locate( stashcode_prog_sec, stashcode_mv,              &
               fields_out, field_count_out, pos_mv)
mv_out => fields_out(pos_mv)
CALL rcf_alloc_field( mv_out )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SELECT CASE (qprofile_number)

CASE (qp_BS1999)

  CALL umprint('Setting moisture profile for gas giant atmosphere', &
               level=PrStatus_Normal, src=RoutineName)

  ALLOCATE(p_on_theta(output_grid % loc_p_field, output_grid % model_levels))
  ALLOCATE(t_on_theta(output_grid % loc_p_field, output_grid % model_levels))
  ALLOCATE(mv_BS1999(output_grid % loc_p_field, output_grid % model_levels))

  DO k = 1, output_grid % model_levels
    DO i = 1, exner_out % level_size
      ! Compute p on theta levels
      p_on_theta(i,k) = p_zero * (intw_rho2w(k,1) * exner_out % data(i,k+1)   &
                                + intw_rho2w(k,2) * exner_out % data(i,k))    &
                                  **recip_kappa
      ! Compute T on theta levels
      t_on_theta(i,k) = thetavd_out % data(i,k) *                             &
                       (intw_rho2w(k,1) * exner_out % data(i,k+1)             &
                      + intw_rho2w(k,2) * exner_out % data(i,k))
    END DO
  END DO

  CALL calc_gas_mixing_ratio_BS1999(p_on_theta, t_on_theta, ip_h2o,           &
                                    output_grid % loc_p_field,                &
                                    output_grid % model_levels, mv_BS1999)

  ! Data array runs from 1:model_levels+1, not 0:model_levels.
  DO k = 2, output_grid % model_levels + 1
    DO i = 1, output_grid % loc_p_field
      mv_out % data(i,k) = mv_BS1999(i,k-1)
    END DO
  END DO

  DO i = 1, mv_out % level_size
    mv_out % data(i,1) = mv_out % data(i,2)
  END DO

  IF (PrintStatus >= PrStatus_Diag) THEN
    CALL umprint('Moisture profile of first point:', src=RoutineName)
    WRITE(ummessage,'(A)') 'k   xi3_at_theta(1,k)  mv(1,k)'
    CALL umprint(ummessage, src=RoutineName)
    DO k = 0, output_grid % model_levels
      WRITE(ummessage,'(I0,2(1X,E16.8))') &
            k, xi3_at_theta(1,k), mv_out % data(1,k+1)
    END DO
  END IF

  DEALLOCATE(mv_BS1999)
  DEALLOCATE(t_on_theta)
  DEALLOCATE(p_on_theta)

  ! Write out new mv field
  CALL rcf_write_field (fields_out(pos_mv), hdr_out, decomp_rcf_output)

CASE (qp_qsat)

  ! The following code is faulty because:
  ! 1. It's a hard-wired profile; mv_init_data should be used instead.
  ! 2. It uses the qsat function which returns specific humidity,
  !    but a mixing ratio is required here.
  ! Therefore abort to stop this code being used.

  ErrorStatus = 100
  WRITE(cmessage,'(A,I0,A)') 'Running with qprofile_number = ',qp_qsat, &
      ' is currently disabled: this code requires an update.'
  CALL ereport(RoutineName, ErrorStatus, cmessage)

  CALL umprint('Setting idealised moisture profile', level=PrStatus_Normal,  &
               src=RoutineName)

  ! Exner_surf also required:
  CALL rcf_locate( stashcode_prog_sec, stashcode_exner_surf,      &
                 fields_out, field_count_out, pos_exner_surf)
  exner_surf_out => fields_out(pos_exner_surf)
  CALL rcf_alloc_field( exner_surf_out )
  CALL rcf_read_field(exner_surf_out, hdr_out, decomp_rcf_output)

  sat_z0 = 1000.0 + planet_radius
  sat_z1 = 20000.0 + planet_radius
  sat_pc = 0.95
  sat_min = 0.01
  sat_grad = (sat_pc - sat_min) / (sat_z1-sat_z0)

  DO k = 1, output_grid % model_levels + 1
    DO i = 1, output_grid % loc_p_field

      ! Compute T & p on theta levels
      IF (k == 1) THEN
        t_tmp = thetavd_out % data(i,k) * exner_surf_out % data(i,1)
        p_tmp = p_zero * exner_surf_out % data(i,1)**recip_kappa
      ELSE
        ! theta runs from 0:model_levels -> data(1:model_levels+1), so use k+1:
        t_tmp = thetavd_out % data(i,k+1) *                                  &
                    (intw_rho2w(k,1) * exner_out % data(i,k+1) +             &
                     intw_rho2w(k,2) * exner_out % data(i,k))
        p_tmp = p_zero * (intw_rho2w(k,1) * exner_out % data(i,k+1) +        &
                          intw_rho2w(k,2) * exner_out % data(i,k))**recip_kappa
      END IF
      ! DEPENDS ON: qsat
      CALL qsat(mv_out % data(i,k),t_tmp, p_tmp, 1)

      ! Again: xi3_at_theta runs from 0:ML+1, so xi3_at_theta(i,0) <--> mv(i,1).
      IF (xi3_at_theta(i,k-1) < sat_z0) THEN
        mv_out % data(i,k)  = sat_pc * mv_out % data(i,k)
      ELSE IF (xi3_at_theta(i,k-1) < sat_z1) THEN
        mv_out % data(i,k) = (sat_pc                                         &
                            - sat_grad * (xi3_at_theta(i,k-1) - sat_z0)      &
                            - sat_z0) * mv_out % data(i,k)
      ELSE
        mv_out % data(i,k) = 0.0
      END IF
    END DO
  END DO

  ! Write out new mv field
  CALL rcf_write_field (fields_out(pos_mv), hdr_out, decomp_rcf_output)

  IF (PrintStatus >= PrStatus_Diag) THEN
    CALL umprint('Moisture profile of first point:', src=RoutineName)
    WRITE(ummessage,'(A)') 'k   xi3_at_theta(1,k)  mv(1,k)'
    CALL umprint(ummessage, src=RoutineName)
    DO k = 0, output_grid % model_levels
      WRITE(ummessage,'(I0,2(1X,E16.8))') &
            k, xi3_at_theta(1,k), mv_out % data(1,k+1)
    END DO
  END IF

END SELECT

! Deallocate the relevant fields
CALL rcf_dealloc_field (fields_out(pos_mv))
CALL rcf_dealloc_field (fields_out(pos_exner))
CALL rcf_dealloc_field (fields_out(pos_thetavd))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_vapour_calc
END MODULE rcf_ideal_vapour_calc_mod
