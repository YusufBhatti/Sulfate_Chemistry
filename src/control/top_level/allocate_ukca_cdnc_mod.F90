! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Purpose: To ALLOCATE ukca_cdnc%cdnc & ukca_cdnc%cdnc3
!            This is here to keep atm_step_4A tidy

MODULE allocate_ukca_cdnc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_UKCA_CDNC_MOD'

CONTAINS

SUBROUTINE allocate_ukca_cdnc( ukca_cdnc,                                      &
                               cdnc_dim1,                                      &
                               cdnc_dim2,                                      &
                               cdnc_dim3,                                      &
                               stashwork54 )

USE atm_step_local,         ONLY: &
    first_atmstep_call

USE ereport_mod,            ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE glomap_clim_cndc_mod,   ONLY: &
    glomap_clim_arg_act_get_cdnc

USE glomap_clim_option_mod, ONLY: &
    l_glomap_clim_arg_act,        &
    l_glomap_clim_aie1,           &
    l_glomap_clim_aie2,           &
    l_glomap_clim_radaer,         &
    l_glomap_mode_clim,           &
    i_glomap_clim_setup,          &
    i_gc_sussocbc_5mode

USE nlsizes_namelist_mod,   ONLY: &
    row_length,                   &
    rows,                         &
    model_levels

USE parkind1,               ONLY: &
    jprb,                         &
    jpim

USE stash_array_mod,        ONLY: &
    stash_maxlen

USE submodel_mod,           ONLY: &
    atmos_im

USE ukca_cdnc_mod,          ONLY: &
    ukca_cdnc_get,                &
    ukca_cdnc_struct

USE ukca_mode_setup,        ONLY: &
    ukca_mode_sussbcoc_5mode

USE ukca_option_mod,        ONLY: &
    l_ukca

USE um_stashcode_mod,       ONLY: &
    stashcode_glomap_clim_sec

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook

IMPLICIT NONE

! Arguments

TYPE (ukca_cdnc_struct), INTENT(INOUT) :: ukca_cdnc
INTEGER,                 INTENT(INOUT) :: cdnc_dim1
INTEGER,                 INTENT(INOUT) :: cdnc_dim2
INTEGER,                 INTENT(INOUT) :: cdnc_dim3
REAL, OPTIONAL,          INTENT(INOUT) :: stashwork54(stash_maxlen(            &
                                          stashcode_glomap_clim_sec, atmos_im))

! Local variables

INTEGER :: errorstatus ! errorstatus

CHARACTER(LEN=errormessagelength) :: Cmessage

CHARACTER (LEN=*),  PARAMETER :: RoutineName = 'ALLOCATE_UKCA_CDNC'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays to receive UKCA output.
IF (.NOT. ALLOCATED(ukca_cdnc%cdnc))                                           &
           ALLOCATE(ukca_cdnc%cdnc  ( row_length, rows, model_levels ) )
IF (.NOT. ALLOCATED(ukca_cdnc%cdnc3))                                          &
           ALLOCATE(ukca_cdnc%cdnc3 ( row_length, rows, model_levels ) )

cdnc_dim1 = row_length
cdnc_dim2 = rows
cdnc_dim3 = model_levels

IF ( l_ukca ) THEN
  
  ! Get current UKCA results
  CALL ukca_cdnc_get( first_atmstep_call, errorstatus, Cmessage )
  
  IF (errorstatus /= 0) THEN
    IF (ALLOCATED(ukca_cdnc%cdnc))  DEALLOCATE(ukca_cdnc%cdnc)
    IF (ALLOCATED(ukca_cdnc%cdnc3)) DEALLOCATE(ukca_cdnc%cdnc3)
    CALL ereport( RoutineName, errorstatus, Cmessage )
  END IF
  
  ! Set minimum CDNC value (equivalent to 5 cm-3)
  WHERE (ukca_cdnc%cdnc < 5.0e06) ukca_cdnc%cdnc = 5.0e06
  
ELSE IF ( l_glomap_mode_clim ) THEN
  
  IF ( first_atmstep_call .AND. .NOT. l_glomap_clim_radaer ) THEN
    ! CALL mode setup routine to set modes, components...
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      CALL ukca_mode_sussbcoc_5mode
    CASE DEFAULT
      
      IF ( i_glomap_clim_setup == 0 ) THEN
        errorstatus = 54999
      ELSE
        errorstatus = ABS( i_glomap_clim_setup )
      END IF
      
      Cmessage    = 'Unknown glomap_clim setup integer'
      CALL ereport(RoutineName,errorstatus,Cmessage)
      
    END SELECT
  END IF
  
  ! Currently, only Abdul-Razzak & Ghan used to calculate cdnc
  IF ( l_glomap_clim_arg_act ) THEN
    CALL glomap_clim_arg_act_get_cdnc ( ukca_cdnc, stashwork54 )
  ELSE
    ! Later add method described in A.Jones et al, 1994, Nature 370 450-453
    errorstatus = 54999
    Cmessage    = 'Only Abdul-Razzak & Ghan method calculates cdnc'
    CALL ereport(RoutineName,errorstatus,Cmessage)
  END IF
  
  ! Set minimum CDNC value (equivalent to 5 cm-3)
  WHERE (ukca_cdnc%cdnc < 5.0e06) ukca_cdnc%cdnc = 5.0e06
  
  IF ( l_glomap_clim_arg_act .AND. .NOT.                                       &
       ( l_glomap_clim_aie1 .OR. l_glomap_clim_aie2 ) ) THEN
    
    ! Allocate dummy ukca_cdnc structure to pass to atmos_physics1
    IF (ALLOCATED(ukca_cdnc%cdnc))  DEALLOCATE(ukca_cdnc%cdnc)
    ALLOCATE(ukca_cdnc%cdnc(1,1,1))
    ukca_cdnc%cdnc(:,:,:)  = 0.0
    IF (ALLOCATED(ukca_cdnc%cdnc3)) DEALLOCATE(ukca_cdnc%cdnc3)
    ALLOCATE(ukca_cdnc%cdnc3(1,1,1))
    ukca_cdnc%cdnc3(:,:,:) = 0.0
    
    cdnc_dim1 = 1
    cdnc_dim2 = 1
    cdnc_dim3 = 1
  END IF ! glomap_clim ARG not aie1/aie2
  
ELSE ! neither l_ukca nor l_glomap_mode_clim
  
  ! It should not be possible to get this error
  errorstatus = 54999
  Cmessage    = 'Neither l_ukca nor l_glomap_mode_clim.'
  CALL ereport(RoutineName,errorstatus,Cmessage)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE allocate_ukca_cdnc

END MODULE allocate_ukca_cdnc_mod
