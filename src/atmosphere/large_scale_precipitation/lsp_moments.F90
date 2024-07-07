! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Conversion between moments of PSD

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

!Note regarding variable precision:
!This module is used in various places beyond the lsp scheme which may require
!either 64 or 32 bit calculation. Therefore, we can't use the real_lsprec 
!parameter as seen generally in the lsp_* modules. Instead, we create an 
!INTERFACE with both 32 and 64 bit versions available.

MODULE lsp_moments_mod

USE um_types, ONLY: real64, real32

IMPLICIT NONE

PRIVATE
PUBLIC ::  lsp_moments

INTERFACE lsp_moments
  MODULE PROCEDURE lsp_moments_64b, lsp_moments_32b
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_MOMENTS_MOD'

CONTAINS

SUBROUTINE lsp_moments_64b( points, rho, t, qcf, cficei, n_out, moment_out )

  ! microphysics modules
USE mphys_inputs_mod, ONLY: l_psd_global, ai_mod => ai, bi_mod => bi

! General atmosphere modules
USE conversions_mod,  ONLY: zerodegc

! Dr Hook modules
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE vectlib_mod,      ONLY: powr_v    => powr_v_interface,                    &
                            oneover_v => oneover_v_interface,                 &
                            exp_v     => exp_v_interface

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_MOMENTS_64B'
INTEGER, PARAMETER :: prec = real64
#include "lsp_moments.h"
RETURN
END SUBROUTINE lsp_moments_64b

!==============================================================================

SUBROUTINE lsp_moments_32b( points, rho, t, qcf, cficei, n_out, moment_out )

  ! microphysics modules
USE mphys_inputs_mod, ONLY: l_psd_global, ai_mod => ai, bi_mod => bi

! General atmosphere modules
USE conversions_mod,  ONLY: zerodegc => zerodegc_32

! Dr Hook modules
USE yomhook,          ONLY: lhook, dr_hook
USE parkind1,         ONLY: jprb, jpim
USE vectlib_mod,      ONLY: powr_v    => powr_v_interface,                    &
                            oneover_v => oneover_v_interface,                 &
                            exp_v     => exp_v_interface

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_MOMENTS_32B'
INTEGER, PARAMETER :: prec = real32
#include "lsp_moments.h"
RETURN
END SUBROUTINE lsp_moments_32b

END MODULE lsp_moments_mod
