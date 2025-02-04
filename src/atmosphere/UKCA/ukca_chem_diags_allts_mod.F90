! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module to contain code to place chemistry-related diagnostics 
!   that are updated at all time steps into the section 50 
!   STASHwork array. 
!
! Method:
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds, University of Oxford, and
!  The Met Office. See www.ukca.ac.uk.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

MODULE ukca_chem_diags_allts_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CHEM_DIAGS_ALLTS_MOD'

CONTAINS

SUBROUTINE ukca_chem_diags_allts ( &
! IN model dimensions
  row_length, rows, model_levels, n_use_tracers, &
! IN fields for diagnostics
  p_tropopause, all_tracers, p_theta_levels,     &
  p_layer_boundaries, z_top_of_model,            &
! INOUT stash workspace
  len_stashwork, stashwork)

! Description:
!   To place chemistry-related diagnostics that are updated at 
!   all time steps into the section 50 STASHwork array

USE ukca_d1_defs,       ONLY: ukca_diag_sect
USE ukca_strat_update,  ONLY: ukca_calc_ozonecol
USE ukca_cspecies,      ONLY: n_o3
USE ukca_constants,     ONLY: c_o3, dobson

USE submodel_mod,      ONLY: atmos_im
USE stash_array_mod,   ONLY: len_stlist, stindex, stlist,           &
                             num_stash_levels, stash_levels, si, sf
USE um_parvars,   ONLY: at_extremity
USE ereport_mod,  ONLY: ereport
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN) :: row_length        ! Model dimensions
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: n_use_tracers

! Diagnostic tracers
! pv-theta Trop surface (Pa)
REAL, INTENT(IN) :: p_tropopause(row_length,rows)
! Tracer array
REAL, INTENT(IN) :: all_tracers(row_length,rows,model_levels,n_use_tracers)
! pressure on theta levels
REAL, INTENT(IN) :: p_theta_levels(row_length,rows,model_levels)
! interface pressures
REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,0:model_levels)
! height at top of model (m)
REAL, INTENT(IN) :: z_top_of_model

! Diagnostics info
INTEGER, INTENT(IN) :: len_stashwork ! Length of diagnostics array
REAL, INTENT(INOUT) :: stashwork(len_stashwork)  ! STASH workspace

! Local variables
INTEGER :: item                                ! STASH item
INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0  ! DrHook tracing entry
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 ! DrHook tracing exit
INTEGER :: icode                               ! error code for EReport
INTEGER :: im_index                            ! internal model index
INTEGER :: ii                                  ! loop variable

! Ozone col (from model level to top of atmosphere) in Dobson Unit
REAL            :: o3col_du(row_length,rows,model_levels)

REAL(KIND=jprb) :: zhook_handle ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEM_DIAGS_ALLTS'
                                              ! used for EReport
CHARACTER(LEN=errormessagelength) :: cmessage ! used for EReport

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0 ! Initialise error status
im_index = 1

! ----------------------------------------------------------------------
!   Copy diagnostic information to STASHwork for STASH processing
! ----------------------------------------------------------------------
! DIAG.50219 Ozone column in DU 
! ----------------------------------------------------------------------
item = 219

IF (sf(item,ukca_diag_sect)) THEN

  ! Calculate ozone column from actual post-chemistry
  ! ozone, for diagnostic purposes. As ozone is tracer,
  ! this makes the ozone column available on all timesteps
  IF (n_o3 > 0) THEN
     CALL ukca_calc_ozonecol(model_levels, rows, row_length,     &
           z_top_of_model,p_layer_boundaries, p_theta_levels,    &
           all_tracers(:,:,:,n_o3)/c_o3,                         &
           o3col_du(:,:,:)) 
     ! convert to Dobson Units from molecules/cm2
     o3col_du(:,:,:) = o3col_du(:,:,:)/(dobson * 1.0e-4) 
                                    
  END IF   ! n_o3

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,UKCA_diag_sect,im_index)),        &
       o3col_du(:,:,:),                                                 &
       row_length,rows,model_levels,0,0,0,0,at_extremity,               &
       stlist(1,stindex(1,item,UKCA_diag_sect,im_index)),len_stlist,    &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags_allts : error in copydiag_3d 50.219"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------
! DIAG.50224 pv-theta Trop surface
! ----------------------------------------------------------------------
item = 224

IF (sf(item,ukca_diag_sect)) THEN

! DEPENDS ON: copydiag
  CALL copydiag (stashwork(si(item,UKCA_diag_sect,im_index)),           &
       p_tropopause(:,:),                                               &
       row_length,rows,0,0,0,0,at_extremity,                            &
       atmos_im,ukca_diag_sect,item,icode,cmessage)

  IF (icode >  0) THEN
    cmessage="ukca_chem_diags_allts : error in copydiag 50.224"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF  ! sf(item,UKCA_diag_sect)

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chem_diags_allts

END MODULE ukca_chem_diags_allts_mod
