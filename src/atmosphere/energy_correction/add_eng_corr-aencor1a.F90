! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE ADD_ENG_CORR
!
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              - TO ADD IN TEMPERATURE CORRECTION TO
!                GLOBAL TEMPERATURE FIELD SO TO
!                CONSERVE TOTAL ENERGY GLOBALLY
!
!    NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!
!    DOCUMENTATION :
!
!----------------------------------------------------------------------
!
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Energy Correction
MODULE add_eng_corr_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ADD_ENG_CORR_MOD'

CONTAINS

SUBROUTINE add_eng_corr (energy_corr,t,                           &
                         STASHwork14)

USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE ereport_mod,            ONLY: ereport
USE um_parvars,             ONLY: at_extremity
USE submodel_mod,           ONLY: atmos_im
USE stash_array_mod,        ONLY: len_stlist, stindex, stlist,    & 
                                  num_stash_levels, stash_levels, &
                                  si, sf
USE atm_fields_bounds_mod,  ONLY: tdims
USE timestep_mod,           ONLY: tstep => timestep
USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod,       ONLY: model_type, mt_single_column

IMPLICIT NONE
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
REAL :: energy_corr         ! IN ENERGY CORRECTION
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!
! INOUT sum of temperature increments
REAL :: t(tdims%i_start:tdims%i_end,                              &
          tdims%j_start:tdims%j_end,                              &
                      1:tdims%k_end)
!  diagnostics  out
REAL :: STASHwork14(*)   ! STASH workspace
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
REAL ::                                                           &
  work(tdims%i_start:tdims%i_end,                                 &
       tdims%j_start:tdims%j_end,                                 &
                   1:tdims%k_end)

INTEGER :: i,j,k,                                                 &
                             ! LOOP COUNTERs
  icode,                                                          &
                   ! return code
  im_index         ! model index

CHARACTER(LEN=errormessagelength) ::                              &
  cmessage        ! return message
CHARACTER(LEN=*), PARAMETER :: RoutineName ='ADD_ENG_CORR'


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
!----------------------------------------------------------------------
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode=0      ! initialise return code to zero

!----------------------------------------------------------------------
! CORRECT TEMPERATURE FOR ERROR IN ENERGY BUDGET OF THE
! PREVIOUS DAY
!----------------------------------------------------------------------
!
DO k=1,tdims%k_end
  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      t(i,j,k) = t(i,j,k) + energy_corr*tstep
    END DO
  END DO
END DO

IF (model_type /= mt_single_column) THEN
  !-----------------------------------------------------------------------
  ! Stash diagnostics
  !-----------------------------------------------------------------------
  ! 14 181  T increment on model levels

  im_index = 1

  IF (sf(181,14)) THEN
    DO k=1,tdims%k_end
      DO j=tdims%j_start,tdims%j_end
        DO i=tdims%i_start,tdims%i_end

          work(i,j,k) = energy_corr*tstep

        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork14(si(181,14,im_index)),             &
         work,tdims%i_len,                                          &
         tdims%j_len,                                               &
         tdims%k_end,                                               &
         0,0,0,0, at_extremity,                                     &
         stlist(1,stindex(1,181,14,im_index)),len_stlist,           &
         stash_levels,num_stash_levels+1,atmos_im,14,181,           &
         icode,cmessage)

    IF (icode >  0) THEN
      cmessage=":error in copydiag_3d(item 181)"//cmessage
    END IF

  END IF

  ! Note old 14201 would be expensive to calculate here as it is
  ! now defined as cv*dt*(column integral of dry mass) and is therefore
  ! calculated and output from section 30.

  IF (icode /= 0) THEN
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF ! model_type
!----------------------------------------------------------------------
! ADD ENERGY CORRECTION INTO SUM OF DIABATIC FLUXES
!----------------------------------------------------------------------
! From UM 5.1 Moved elsewhere in code

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE add_eng_corr

END MODULE
