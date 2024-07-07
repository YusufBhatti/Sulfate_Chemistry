! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_lam_domain_kind_mod
IMPLICIT NONE
INTEGER, ALLOCATABLE, SAVE :: lam_domain_kind(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_LAM_DOMAIN_KIND_MOD'

CONTAINS
SUBROUTINE eg_init_lam_domain_kind(halo_i, halo_j, row_length, rows, &
                                   datastart, global_rows, global_row_length)
! Description:
!              Set flag for each grid to distinguish interior and rim region.
!
! Method: ENDGame formulation version 4.xx
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics advection
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.


USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s
USE lbc_mod,               ONLY: rimwidtha
USE rimtypes,              ONLY: rima_type_norm
USE um_parparams,          ONLY: pwest, peast, pnorth, psouth

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE
INTEGER, INTENT(IN)  :: halo_i
INTEGER, INTENT(IN)  :: halo_j
INTEGER, INTENT(IN)  :: row_length
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: datastart(3)
INTEGER, INTENT(IN)  :: global_rows
INTEGER, INTENT(IN)  :: global_row_length

! local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INIT_LAM_DOMAIN_KIND'
INTEGER              :: nrim
INTEGER              :: i, j
INTEGER              :: gi, gj
LOGICAL, PARAMETER   :: l_only_adjacent = .FALSE.
REAL                 :: i_distance

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (model_type == mt_lam) THEN

  ALLOCATE( lam_domain_kind(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j) )

  nrim = rimwidtha(rima_type_norm)
!$OMP  PARALLEL DO DEFAULT(NONE)                                     &
!$OMP& PRIVATE(j, i, gi, gj, i_distance)                             &
!$OMP& SHARED(halo_j, rows, halo_i ,row_length, datastart,           &
!$OMP&        global_row_length, global_rows, nrim,                  &
!$OMP&        lam_domain_kind)
  DO j = 1-halo_j, rows+halo_j
    DO i = 1-halo_i ,row_length+halo_i
      gi=i+datastart(1)-1
      gj=j+datastart(2)-1
      i_distance = MIN( gi,                         &
                        global_row_length - gi + 1, &
                        gj,                         &
                        global_rows       - gj + 1  )
      IF (i_distance > 2*nrim) THEN
        IF (l_only_adjacent) THEN
          lam_domain_kind(i,j) = 2
        ELSE
          lam_domain_kind(i,j) = 1
        END IF
      ELSE IF (i_distance > nrim) THEN
        lam_domain_kind(i,j) = 1
      ELSE
        lam_domain_kind(i,j) = 0
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_init_lam_domain_kind
END MODULE eg_lam_domain_kind_mod
