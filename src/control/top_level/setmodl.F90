! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Defines submodel and section/version configuration
!
! Subroutine Interface:
SUBROUTINE setmodl(ErrorStatus,cmessage)

#if defined(RECON)
USE Rcf_Grid_Type_Mod, ONLY:                                     &
    Output_Grid

USE Rcf_Lsm_Mod, ONLY:                                           &
    Glob_Land_out

USE model_domain_mod, ONLY: output_grid_stagger

USE Rcf_HeadAddress_Mod, ONLY: &
    FH_GridStagger_C,           &
    FH_GridStagger_Endgame

#endif

USE stash_model_mod, ONLY:                                                     &
    h_atmos, mean_number, h_strat, h_global, h_a_ewspace, h_a_nsspace,         &
    h_a_firstlat, h_a_firstlong, h_a_polelat, h_a_polelong, len_extra,         &
    stlevgwdrag, botvdifflev, topvdifflev, h_orog_rough


USE rad_input_mod, ONLY:                                          &
    a_sw_radstep_diag,                                            &
    a_sw_radstep_prog,                                            &
    a_lw_radstep_diag,                                            &
    a_lw_radstep_prog,                                            &
    H_SWBands,                                                    &
    H_LWBands

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE decomp_params,     ONLY: decomp_standard_atmos
USE Decomp_DB
USE ukca_option_mod,   ONLY: tr_ukca_a
USE ukca_tracer_stash, ONLY: a_max_ukcavars

USE nlstgen_mod, ONLY: meanfreqim
USE submodel_mod, ONLY:                                                        &
    internal_model_list, n_internal_model_max, atmos_sm, atmos_im
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE h_vers_mod, ONLY: h_vers_init, h_vers
USE lam_config_inputs_mod, ONLY: polelata, polelona, frstlata,           &
    frstlona, delta_lat, delta_lon
USE free_tracers_inputs_mod, ONLY: a_max_trvars, tracer_a, i_free_tracer
USE nlsizes_namelist_mod, ONLY:                                          &
    cloud_levels, land_field, model_levels, row_length, rows,            &
    tr_levels, tr_vars, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_cyclic_lam       &
                          , mt_bi_cyclic_lam, mt_smexe

IMPLICIT NONE

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!

! Array arguments with intent(out):
CHARACTER(LEN=errormessagelength) :: cmessage    ! Error return message

! Error status:
INTEGER ::     ErrorStatus ! +ve = fatal error


! Local scalars
REAL :: ASteps !Atmos timesteps per hour
INTEGER :: i,j
INTEGER :: Im

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETMODL'

!- End of Header ---------------------------------------------------

!  Define submodel configuration

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
h_atmos = 'N'

DO i = 1,n_internal_model_max
  IF (internal_model_list(i) == atmos_im) THEN
    h_atmos ='Y'
  END IF
END DO

DO i = 1,n_internal_model_max
  h_vers(i,0)= 1
END DO

IF (h_atmos  ==  'Y') THEN
  ! Initialise the submodel configurations as required for version mask
  CALL h_vers_init()
END IF

! Submodel independent
#if !defined(RECON)
DO i=1,n_internal_model_max
  mean_number(i)=0
  DO j=1,4
    IF (meanfreqim(j,i) >  0) THEN
      mean_number(i)=mean_number(i)+1
    END IF
  END DO
END DO
#endif

Im = atmos_im
IF (h_atmos  ==  'Y') THEN
  ! Atmos model included
  h_strat='N'
  DO i=1,a_max_trvars   ! Up to A_MAX_TRVARS=150 free tracers
    IF (i_free_tracer(i) == 0) THEN
      tracer_a(i)=.FALSE.
    ELSE
      tracer_a(i)=.TRUE.
    END IF
  END DO

  IF (model_type == mt_global) THEN
    ! Atmos global model
    h_global(atmos_im)='Y'

#if defined(RECON)
    h_a_ewspace=360.0/ Output_Grid % glob_p_row_length
    IF (output_grid_stagger == FH_GridStagger_C) THEN
      h_a_nsspace=180.0/( Output_Grid % glob_p_rows-1)
    ELSE
      h_a_nsspace=180.0/( Output_Grid % glob_p_rows)
    END IF
#else
    h_a_ewspace=360.0/                                            &
      decompDB(decomp_standard_atmos)%glsize(1,fld_type_p)
      h_a_nsspace=180.0/                                          &
        (decompDB(decomp_standard_atmos)%glsize(2,fld_type_p))
#endif

    h_a_firstlat=-90.0       ! S to N
    h_a_firstlong=0.0
    h_a_polelat=90.0
    h_a_polelong=0.0
  ELSE IF (model_type == mt_lam .OR.        &
           model_type == mt_cyclic_lam .OR. &
           model_type == mt_bi_cyclic_lam) THEN
    ! Atmos LAM
    h_global(atmos_im)='N'
    h_a_ewspace=delta_lon
    h_a_nsspace=delta_lat
    h_a_firstlat=frstlata
    h_a_firstlong=frstlona
    IF (h_a_firstlong <  0) h_a_firstlong=h_a_firstlong+360.0
    h_a_polelat=polelata
    h_a_polelong=polelona
  ELSE IF (model_type /= mt_smexe) THEN
    WRITE(umMessage,'(A,I3)')                                    &
   'Setmodl: UNEXPECTED ATMOSPHERIC AREA CODE model_type',model_type
    CALL umPrint(umMessage,src='setmodl')
  END IF

#if defined(RECON)
  len_extra(atmos_sm) = (1+ 3 * Output_Grid % model_levels ) *             &
                         Output_Grid % glob_p_rows *                       &
                         Output_Grid % glob_p_row_length
#else
  len_extra(atmos_sm) = (1+3*model_levels)*rows*row_length
#endif

ELSE   ! Atmosphere model not included

  h_strat       ='N'
  h_global(atmos_im)='N'

#if defined(RECON)
  glob_land_out =0
  Output_Grid % cloud_levels       =0
  Output_Grid % tr_levels          =0
  tr_vars                          =0
#else
  land_field    =0
  rows          =0
  model_levels  =0
  cloud_levels  =0
  tr_levels     =0
  tr_vars       =0
#endif
  DO i=1,a_max_trvars   ! Up to A_MAX_TRVARS=150 free tracers
    tracer_a(i) =.FALSE.
  END DO
  DO i=1,a_max_ukcavars ! Up to A_MAX_UKCAVARS=150 UKCA tracers
    tr_ukca_a(i) =.FALSE.
  END DO
  StLevGWdrag   =0
  BotVDiffLev   =0
  TopVDifflev   =0
  a_sw_radstep_diag  = 0
  a_sw_radstep_prog  = 0
  a_lw_radstep_diag  = 0
  a_lw_radstep_prog  = 0
  h_swbands     =0
  h_lwbands     =0
  h_a_ewspace   =0.0
  h_a_nsspace   =0.0
  h_a_firstlat  =0.0
  h_a_firstlong =0.0
  h_a_polelat   =0.0
  h_a_polelong  =0.0
END IF


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setmodl
