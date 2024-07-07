! *****************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*********************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input

MODULE lbc_calc_size_mod

! Proves the SUBROUTINE lbc_calc_size, which can be used to calculate the sizes
! of LBC rows and lengths, as well as their start and end points on this PE.

IMPLICIT NONE
PRIVATE

PUBLIC :: lbc_calc_size

!------------------------------------------------------------------------------!
CONTAINS
!------------------------------------------------------------------------------!

SUBROUTINE lbc_calc_size(iside, full_lbc_row_len, full_lbc_nrows,              &
                         decomp_lbc_row_len, decomp_lbc_nrows,                 &
                         full_lbc_start_pt, decomp_lbc_start_pt, fld_type,     &
                         halo_type, rim_type, proc)

USE UM_ParVars, ONLY:                                                          &
    glsize, halosize, g_lasize, g_at_extremity, g_datastart

USE lbc_mod, ONLY:                                                             &
    g_lbc_starta, rimwidtha, global_lbc_starta

USE um_parparams, ONLY:                                                        &
    pnorth, peast, psouth, pwest

USE ereport_mod, ONLY:                                                         &
    ereport

USE errormessagelength_mod, ONLY:                                              &
    errormessagelength

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) ::                                                         &
  iside,                                                                       &
  proc,                                                                        &
  fld_type,                                                                    &
  halo_type,                                                                   &
  rim_type

INTEGER, INTENT(OUT) ::                                                        &
  full_lbc_row_len,                                                            &
  full_lbc_nrows,                                                              &
  decomp_lbc_row_len,                                                          &
  decomp_lbc_nrows,                                                            &
  full_lbc_start_pt,                                                           &
  decomp_lbc_start_pt

! Local Variables

INTEGER ::                                                                     &
  first_lbc_pt,                                                                &
  first_lbc_row,                                                               &
  icode

CHARACTER(LEN=errormessagelength) ::                                           &
  cmessage

CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='LBC_CALC_SIZE'

! The subroute calculates the following values for the size of the LBC:
!
! * full_lbc_row_len    : East-West dimension of the full lbc side
!
! * full_lbc_nrows      : North-South dimension of the full lbc side
!
! * decomp_lbc_row_len  : East-West dimension of the decomposed lbc side on
!                         this PE.
!
! * decomp_lbc_nrows    : North-South dimension of the decomposed lbc side on
!                         this PE.
!
! * full_lbc_start_pt   : First point of the decomposed lbc side inside the 1d
!                         full_lbc array
!
! * decomp_lbc_start_pt : First point of the decomposed lbc side inside the 1d
!                         decomposed lbc array

SELECT CASE (iside)

  CASE (pnorth, psouth)
    ! Calculate size of FULL_LBC

    full_lbc_row_len = glsize(1,fld_type) + 2*halosize(1,halo_type)
    full_lbc_nrows   = halosize(2,halo_type) + rimwidtha(rim_type)

    ! Calculate size of DECOMP_LBC

    decomp_lbc_row_len = g_lasize(1,fld_type,halo_type,proc)
    decomp_lbc_nrows   = halosize(2,halo_type) + rimwidtha(rim_type)

    ! Calculate first point of DECOMP_LBC in FULL_LBC

    first_lbc_pt  = g_datastart(1,proc)
    first_lbc_row = 1

  CASE (peast, pwest)
    ! Calculate size of FULL_LBC

    full_lbc_row_len = halosize(1,halo_type) + rimwidtha(rim_type)
    full_lbc_nrows   = glsize(2,fld_type) - 2*rimwidtha(rim_type)

    ! Calculate size of DECOMP_LBC

    decomp_lbc_row_len = halosize(1,halo_type) + rimwidtha(rim_type)
    decomp_lbc_nrows = g_lasize(2,fld_type,halo_type,proc)

    IF (g_at_extremity(pnorth,proc)) THEN
      decomp_lbc_nrows = decomp_lbc_nrows - halosize(2,halo_type)              &
                         - rimwidtha(rim_type)
    END IF

    IF (g_at_extremity(psouth,proc)) THEN
      decomp_lbc_nrows = decomp_lbc_nrows - halosize(2,halo_type)              &
                         - rimwidtha(rim_type)
    END IF

    ! Calculate first point of DECOMP_LBC in FULL_LBC

    first_lbc_pt  = 1
    first_lbc_row = g_datastart(2,proc)

    IF (.NOT. g_at_extremity(psouth,proc)) THEN
      first_lbc_row = first_lbc_row - rimwidtha(rim_type)                      &
                      - halosize(2,halo_type)
    END IF

  CASE DEFAULT
    ! Error - invalid side.
    WRITE(cmessage,'(A,I0,A)') RoutineName // ': iside = ',                    &
                               iside, ' is not a valid LBC side code'
    icode = 1
    CALL ereport(RoutineName,icode,cmessage)

END SELECT

full_lbc_start_pt = global_lbc_starta(iside,fld_type,halo_type,rim_type)       &
                    + (first_lbc_row-1)*full_lbc_row_len + first_lbc_pt - 1

decomp_lbc_start_pt = g_lbc_starta(iside,fld_type,halo_type,rim_type,proc)

END SUBROUTINE lbc_calc_size

!------------------------------------------------------------------------------!

END MODULE lbc_calc_size_mod