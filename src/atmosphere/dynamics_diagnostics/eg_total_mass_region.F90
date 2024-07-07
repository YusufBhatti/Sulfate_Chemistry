! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_total_mass_region_mod
IMPLICIT NONE

! Description:
!             This routine gives indices of the region to compute `total' mass
!
! Method: ENDGame formulation version 4.xx
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dynamics diagnostics
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.

CONTAINS
SUBROUTINE eg_total_mass_region(IS, ie, js, je, l_exclude_rim)

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_s
USE um_parvars,            ONLY: at_extremity
USE lbc_mod,               ONLY: rimwidtha
USE rimtypes,              ONLY: rima_type_norm
USE um_parparams,          ONLY: pwest, peast, pnorth, psouth

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE
LOGICAL, INTENT(IN)  :: l_exclude_rim
INTEGER, INTENT(OUT) :: IS, ie, js, je

! local variables
INTEGER              :: nrim

nrim = rimwidtha(rima_type_norm)
IS   = tdims%i_start
ie   = tdims%i_end
js   = tdims%j_start
je   = tdims%j_end

IF (model_type == mt_lam .AND. l_exclude_rim) THEN
  IF (at_extremity(pwest)) THEN
    IS = IS + nrim
  END IF
  IF (at_extremity(peast)) THEN
    ie = ie - nrim
  END IF
  IF (at_extremity(psouth)) THEN
    js = js + nrim
  END IF
  IF (at_extremity(pnorth)) THEN
    je = je - nrim
  END IF
END IF

RETURN
END SUBROUTINE eg_total_mass_region
END MODULE eg_total_mass_region_mod
