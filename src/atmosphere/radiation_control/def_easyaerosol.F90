! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module defining and managing structures to store EasyAerosol 
! distributions.
!
! Method:
!
!  Provide structure and memory allocation/deallocation routines.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Contained subroutines in this module:
!   allocate_easyaerosol_rad           (Public)
!   allocate_easyaerosol_cdnc          (Public)
!   deallocate_easyaerosol_rad         (Public)
!   deallocate_easyaerosol_cdnc        (Public)
!
! --------------------------------------------------------------------------
MODULE def_easyaerosol

IMPLICIT NONE
  
!
! Optical properties: 3D distributions with also a dependence on
! spectral waveband.
!
TYPE t_easyaerosol_rad
  INTEGER :: dim1
  INTEGER :: dim2
  INTEGER :: dim3
  INTEGER :: dim4
  REAL, ALLOCATABLE :: extinction(:,:,:,:)
  REAL, ALLOCATABLE :: absorption(:,:,:,:)
  REAL, ALLOCATABLE :: asymmetry(:,:,:,:)
END TYPE t_easyaerosol_rad
  
!
! Cloud droplet number concentrations: 3D distributions
!
TYPE t_easyaerosol_cdnc
  INTEGER dim1
  INTEGER dim2
  INTEGER dim3
  REAL, ALLOCATABLE :: cdnc(:,:,:)
END TYPE t_easyaerosol_cdnc
 
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEF_EASYAEROSOL'
 
CONTAINS
  
!
! Memory allocation routines
!
SUBROUTINE allocate_easyaerosol_rad(rad, row_length, rows, model_levels, &
                                    n_wavebands)
 
  USE yomhook,         ONLY: lhook, dr_hook
  USE parkind1,        ONLY: jprb, jpim
 
  IMPLICIT NONE
    
  TYPE (t_easyaerosol_rad), INTENT(INOUT) :: rad
  INTEGER, INTENT(IN) :: row_length
  INTEGER, INTENT(IN) :: rows
  INTEGER, INTENT(IN) :: model_levels
  INTEGER, INTENT(IN) :: n_wavebands
   
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_EASYAEROSOL_RAD'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_in, zhook_handle)
 
  rad%dim1 = row_length
  rad%dim2 = rows
  rad%dim3 = model_levels
  rad%dim4 = n_wavebands
    
  IF (.NOT. ALLOCATED(rad%extinction)) THEN
    ALLOCATE(rad%extinction(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
  END IF
  rad%extinction(:,:,:,:) = 0.0

  IF (.NOT. ALLOCATED(rad%absorption)) THEN
    ALLOCATE(rad%absorption(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
  END IF
  rad%absorption(:,:,:,:) = 0.0

  IF (.NOT. ALLOCATED(rad%asymmetry)) THEN
    ALLOCATE(rad%asymmetry(rad%dim1, rad%dim2, rad%dim3, rad%dim4))
  END IF
  rad%asymmetry(:,:,:,:) = 0.0

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_out, zhook_handle)
  RETURN

END SUBROUTINE allocate_easyaerosol_rad

SUBROUTINE allocate_easyaerosol_cdnc(easy_cdnc, row_length, rows, model_levels)
    
  USE yomhook,         ONLY: lhook, dr_hook
  USE parkind1,        ONLY: jprb, jpim 
    
  IMPLICIT NONE
    
  TYPE (t_easyaerosol_cdnc), INTENT(INOUT) :: easy_cdnc
  INTEGER, INTENT(IN) :: row_length
  INTEGER, INTENT(IN) :: rows
  INTEGER, INTENT(IN) :: model_levels

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ALLOCATE_EASYAEROSOL_CDNC'

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_in, zhook_handle)
    
  easy_cdnc%dim1 = row_length
  easy_cdnc%dim2 = rows
  easy_cdnc%dim3 = model_levels
    
  IF (.NOT. ALLOCATED(easy_cdnc%cdnc)) THEN
    ALLOCATE(easy_cdnc%cdnc(easy_cdnc%dim1, easy_cdnc%dim2, easy_cdnc%dim3))
  END IF
  easy_cdnc%cdnc(:,:,:) = 0.0

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, &
                          zhook_out, zhook_handle)
  RETURN

END SUBROUTINE allocate_easyaerosol_cdnc

!
! Memory deallocation routines
!
SUBROUTINE deallocate_easyaerosol_rad(rad)

  IMPLICIT NONE
    
  TYPE (t_easyaerosol_rad), INTENT(INOUT) :: rad
   
  IF (ALLOCATED(rad%extinction)) DEALLOCATE(rad%extinction)
  IF (ALLOCATED(rad%absorption)) DEALLOCATE(rad%absorption)
  IF (ALLOCATED(rad%asymmetry))  DEALLOCATE(rad%asymmetry)
 
  RETURN
 
END SUBROUTINE deallocate_easyaerosol_rad

SUBROUTINE deallocate_easyaerosol_cdnc(easy_cdnc)
 
  IMPLICIT NONE

  TYPE (t_easyaerosol_cdnc), INTENT(INOUT) :: easy_cdnc

  IF (ALLOCATED(easy_cdnc%cdnc)) DEALLOCATE(easy_cdnc%cdnc)
 
  RETURN
 
END SUBROUTINE deallocate_easyaerosol_cdnc

END MODULE def_easyaerosol
