! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Update LBCs
SUBROUTINE lbc_update(                                                    &
     submodel,icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE lbc_read_data_mod, ONLY: albc_num, albc_swapstep
USE submodel_mod, ONLY: atmos_im
USE nlstcall_mod, ONLY: Num_ALBCs
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup,      &
    len_dumphist, len_fixhd, len_tot,                                  &
    model_levels, mpp_len1_lookup, n_cca_lev, n_obj_d1_max, sm_levels, &
    st_levels, tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca,   &
    tr_vars

USE model_time_mod, ONLY: &
    stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

! Description
!  Does an LBC update

! Method
!  This file used to control the updating of LBCs when using coupling - where
! a driving model generated LBCs which a LAM then received in parallel. This
! functionality has now been retired, and this routine now only deals with
! updating the LBCs from a file without the overhead of the coupling with a
! second model.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC input

! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v8.2 programming standards.


INTEGER,           INTENT(IN) :: submodel! Submodel id
INTEGER,           INTENT(OUT) :: icode   ! =0 normal exit; >0 error exit
CHARACTER(LEN=errormessagelength),INTENT(OUT) :: cmessage! Error message


CHARACTER(LEN=8)  :: ch_date2    ! Date returned from date_and_time
CHARACTER(LEN=10) :: ch_time2    ! Time returned from date_and_time
INTEGER           :: lbc_ntimes  ! No of BCs in communication file.
INTEGER           :: ms_ntimes   ! No of BCs required in mesoscale.
INTEGER           :: len_wait_tot!  Total wait for availability of BCs
INTEGER           :: um_lbc_wait_usec ! Total number of microseconds to sleep.

LOGICAL           :: l_active    !  T : Output stream active for LBCs.

! Error reporting
CHARACTER(LEN=*) ::     RoutineName

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
PARAMETER (RoutineName='LBC_UPDATE')

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (model_type /= mt_global) THEN

  ! Swap to second atmos boundary file?
  IF (num_albcs == 2 .AND. stepim(atmos_im) == albc_swapstep) THEN
    albc_num = 2
    IF (printstatus >= prstatus_normal) THEN
      WRITE(umMessage,'(A)') ''
      CALL umPrint(umMessage,src='lbc_update')
      WRITE(umMessage,'(A)') 'U_MODEL: Swapping to 2nd atmos boundary '//    &
           'file'
      CALL umPrint(umMessage,src='lbc_update')
      WRITE(umMessage,'(A,I8)') '         Step = ', stepim(atmos_im)
      CALL umPrint(umMessage,src='lbc_update')
      WRITE(umMessage,'(A)') ''
      CALL umPrint(umMessage,src='lbc_update')
    END IF
  END IF

  IF (num_albcs == 2) THEN
    IF (stepim(atmos_im) == albc_swapstep) THEN ! Swap to 2nd bndy file

       ! DEPENDS ON: inbounda
      CALL inbounda(                                                         &
           a_len1_levdepc,a_len2_levdepc,                                    &
           a_len1_rowdepc,a_len2_rowdepc,                                    &
           a_len1_coldepc,a_len2_coldepc )

    END IF
  END IF
END IF   ! if .not. GLOBAL

! DEPENDS ON: up_bound
CALL up_bound(submodel,                                                      &
     icode,cmessage)

IF (icode  /=  0) CALL ereport(routinename,icode,cmessage)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lbc_update
