! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:  Global data module to store physical tendencies
!               from MetUM parametrizations. It is called by 
!               schemes such as the SPT which perturbs the physical
!               tendencies, or the PV-tracer scheme.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Top-level
MODULE physics_tendencies_mod

IMPLICIT NONE

REAL, ALLOCATABLE ::                                                    &
! Rad (SW and LW)
 dt_sw(:,:,:)                                                           &
,dq_sw(:,:,:)                                                           &
,dt_lw(:,:,:)                                                           &
,dq_lw(:,:,:)                                                           &
! Microphysics
,dt_mic(:,:,:)                                                          &
,dq_mic(:,:,:)                                                          &
! GWD
,dt_gwd(:,:,:)                                                          &
,du_gwd(:,:,:)                                                          &
,dv_gwd(:,:,:)                                                          &
! All slow physics combined
,dtheta_ph1(:,:,:)                                                      &
! Conv
,dt_conv(:,:,:)                                                         &
,dq_conv(:,:,:)                                                         &
,du_conv(:,:,:)                                                         &
,dv_conv(:,:,:)                                                         &
! Boundary layer
,dt_bl(:,:,:)                                                           &
,du_bl(:,:,:)                                                           &
,dv_bl(:,:,:)                                                           &
,dq_cl_bl(:,:,:)                                                        &
! stochastic physics
,dtheta_stph(:,:,:)                                                     &
,du_stph(:,:,:)                                                         &
,dv_stph(:,:,:)                                                         &
! cloud rebalancing
,dt_cld(:,:,:)                                                          &
! nudging
,dtheta_nud(:,:,:)                                                      &
,du_nud(:,:,:)                                                          &
,dv_nud(:,:,:)

! Set flags to retain tendencies
LOGICAL, SAVE :: l_retain_rad_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_mic_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_gwd_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_slow_tendencies     = .FALSE.
LOGICAL, SAVE :: l_retain_ph1_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_conv_all_tendencies = .FALSE.
LOGICAL, SAVE :: l_retain_conv_tendencies     = .FALSE.
LOGICAL, SAVE :: l_retain_conv_mom_tendencies = .FALSE.
LOGICAL, SAVE :: l_retain_bl_tendencies       = .FALSE.
LOGICAL, SAVE :: l_retain_stph_tendencies     = .FALSE.
LOGICAL, SAVE :: l_retain_cld_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_nud_tendencies      = .FALSE.
LOGICAL, SAVE :: l_retain_q_cl_bl_tendencies  = .FALSE.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHYSICS_TENDENCIES_MOD'

CONTAINS

SUBROUTINE  init_slowphys_tendencies()
! Allocate arrays to store slow physics tendencies of theta, q, u and v

USE atm_fields_bounds_mod, ONLY: udims, vdims, tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_SLOWPHYS_TENDENCIES'

!------------------------------------------------------------------------------


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! +++ Allocate arrays
! SW ten
IF (.NOT. ALLOCATED(dt_sw))                                           &
  ALLOCATE ( dt_sw(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   1:tdims%k_end))

IF (.NOT. ALLOCATED(dq_sw))                                           &
  ALLOCATE ( dq_sw(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   1:tdims%k_end))

! LW tend
IF (.NOT. ALLOCATED(dt_lw))                                           &
  ALLOCATE ( dt_lw(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   1:tdims%k_end))

IF (.NOT. ALLOCATED(dq_lw))                                           &
  ALLOCATE ( dq_lw(tdims%i_start:tdims%i_end,                         &
                   tdims%j_start:tdims%j_end,                         &
                   1:tdims%k_end))

! microphysics (LS rain)
IF (.NOT. ALLOCATED(dt_mic))                                          &
  ALLOCATE ( dt_mic(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    1:tdims%k_end))

IF (.NOT. ALLOCATED(dq_mic))                                          &
  ALLOCATE ( dq_mic(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    1:tdims%k_end))
! Gravity wave drag
IF (.NOT. ALLOCATED(dt_gwd))                                          &
  ALLOCATE ( dt_gwd(tdims%i_start:tdims%i_end,                        &
                    tdims%j_start:tdims%j_end,                        &
                    1:tdims%k_end))

IF (.NOT. ALLOCATED(du_gwd))                                          &
  ALLOCATE ( du_gwd(udims%i_start:udims%i_end,                        &
                    udims%j_start:udims%j_end,                        &
                    udims%k_start:udims%k_end) )

IF (.NOT. ALLOCATED(dv_gwd))                                          &
  ALLOCATE ( dv_gwd(vdims%i_start:vdims%i_end,                        &
                    vdims%j_start:vdims%j_end,                        &
                    vdims%k_start:vdims%k_end) )

! Initialize all to zero
dt_sw = 0.0
dq_sw = 0.0

dt_lw = 0.0
dq_lw = 0.0

dq_mic= 0.0

dt_gwd = 0.0
du_gwd = 0.0
dv_gwd = 0.0

! Allocate all Theta for slow physics only if requested by the scheme
IF (l_retain_ph1_tendencies) THEN
  IF (.NOT. ALLOCATED(dtheta_ph1))                                    &
    ALLOCATE(dtheta_ph1(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end))
  ! Initialize
  dtheta_ph1 = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE init_slowphys_tendencies
!-----------------------------------------------------------------------

SUBROUTINE init_convection_tendencies
! Allocate arrays to store convective tendencies of theta, q, u and v

USE atm_fields_bounds_mod, ONLY:                                        &
udims, vdims, tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i,j,k !indexes for looping

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_CONVECTION_TENDENCIES'

!------------------------------------------------------------------------------


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! +++++ ALLOCATE T +++++++++++
IF (.NOT. ALLOCATED(dt_conv)) THEN
  ALLOCATE ( dt_conv(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,                         &
                     1:tdims%k_end))
END IF


! Initialize
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dt_conv(i,j,k) = 0.0
    END DO ! i
  END DO ! j
END DO ! k


! +++++ ALLOCATE q +++++++++++
IF (.NOT. ALLOCATED(dq_conv)) THEN
  ALLOCATE ( dq_conv(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,                         &
                     1:tdims%k_end))
END IF
! Initialize
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dq_conv(i,j,k) = 0.0
    END DO ! i
  END DO ! j
END DO ! k

! +++++ ALLOCATE u +++++++++++
IF (.NOT. ALLOCATED(du_conv)) THEN
  ALLOCATE ( du_conv(udims%i_start:udims%i_end,                         &
                     udims%j_start:udims%j_end,                         &
                     udims%k_start:udims%k_end) )

END IF

DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      du_conv(i,j,k) = 0.0
    END DO ! i
  END DO ! j
END DO ! k

! +++++ ALLOCATE v +++++++++++
IF (.NOT. ALLOCATED(dv_conv)) THEN
  ALLOCATE ( dv_conv(vdims%i_start:vdims%i_end,                         &
                     vdims%j_start:vdims%j_end,                         &
                     vdims%k_start:vdims%k_end) )
END IF

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      dv_conv(i,j,k) = 0.0
    END DO ! i
  END DO ! j
END DO ! k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE init_convection_tendencies

SUBROUTINE init_bl_tendencies()
! Allocate arrays to store BL tendencies of theta, u, v

USE atm_fields_bounds_mod, ONLY:                                        &
udims, vdims, tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_BL_TENDENCIES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! U wind
IF (.NOT. ALLOCATED(du_bl)) THEN
  ALLOCATE(du_bl( udims%i_start:udims%i_end,                            &
                  udims%j_start:udims%j_end,                            &
                  udims%k_start:udims%k_end) )
END IF

! V wind
IF (.NOT. ALLOCATED(dv_bl)) THEN
  ALLOCATE(dv_bl( vdims%i_start:vdims%i_end,                            &
                  vdims%j_start:vdims%j_end,                            &
                  vdims%k_start:vdims%k_end) )
END IF

! temperature
IF (.NOT. ALLOCATED(dt_bl)) THEN
  ALLOCATE(dt_bl( tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  1:tdims%k_end) )
END IF

! Initialize
du_bl = 0.0
dv_bl = 0.0
dt_bl = 0.0

IF (l_retain_q_cl_bl_tendencies) THEN
  IF (.NOT. ALLOCATED(dq_cl_bl))                                      &
    ALLOCATE ( dq_cl_bl(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        1:tdims%k_end))

  dq_cl_bl = 0.0
END IF 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE init_bl_tendencies

SUBROUTINE init_stph_tendencies()
! Allocate arrays to store Stochastic Physics tendencies

USE atm_fields_bounds_mod, ONLY:                                        &
udims, vdims, tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_STPH_TENDENCIES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! temperature
IF (.NOT. ALLOCATED(dtheta_stph))                                       &
  ALLOCATE(dtheta_stph( tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,                      &
                        1:tdims%k_end) )

IF (.NOT. ALLOCATED(du_stph))                                           &
  ALLOCATE(du_stph( udims%i_start:udims%i_end,                          &
                    udims%j_start:udims%j_end,                          &
                    udims%k_start:udims%k_end) )


! V wind
IF (.NOT. ALLOCATED(dv_stph))                                           &
  ALLOCATE(dv_stph( vdims%i_start:vdims%i_end,                          &
                    vdims%j_start:vdims%j_end,                          &
                    vdims%k_start:vdims%k_end) )

dtheta_stph  = 0.0
du_stph      = 0.0
dv_stph      = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE init_stph_tendencies

SUBROUTINE init_cld_tendencies()
! Allocate arrays to store cloud tendencies of theta

USE atm_fields_bounds_mod, ONLY: tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_CLD_TENDENCIES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! temperature
IF (.NOT. ALLOCATED(dt_cld)) THEN
  ALLOCATE(dt_cld( tdims%i_start:tdims%i_end,                            &
                   tdims%j_start:tdims%j_end,                            &
                   1:tdims%k_end) )
END IF

! Initialize
dt_cld = 0.0



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE init_cld_tendencies

SUBROUTINE init_nud_tendencies()
! Allocate arrays to store cloud tendencies of theta
! and water species

USE atm_fields_bounds_mod, ONLY: tdims, udims, vdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_NUD_TENDENCIES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(dtheta_nud)) THEN
  ALLOCATE(dtheta_nud( tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,                       &
                       1:tdims%k_end) )
END IF

IF (.NOT. ALLOCATED(du_nud)) THEN
  ALLOCATE(du_nud( udims%i_start:udims%i_end,                           &
                   udims%j_start:udims%j_end,                           &
                   udims%k_start:udims%k_end) )
END IF

IF (.NOT. ALLOCATED(dv_nud)) THEN
  ALLOCATE(dv_nud( vdims%i_start:vdims%i_end,                           &
                   vdims%j_start:vdims%j_end,                           &
                   vdims%k_start:vdims%k_end) )
END IF

! Initialise
dtheta_nud = 0.0
du_nud = 0.0
dv_nud = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------
RETURN
END SUBROUTINE init_nud_tendencies

SUBROUTINE destroy_phy_tend()
! Deallocate physical tendencies at the end of the model run.

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DESTROY_PHY_TEND'
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_retain_q_cl_bl_tendencies) THEN
  IF (ALLOCATED(dq_cl_bl))   DEALLOCATE(dq_cl_bl)
END IF

IF (l_retain_nud_tendencies) THEN
  IF (ALLOCATED(dv_nud))      DEALLOCATE(dv_nud)
  IF (ALLOCATED(du_nud))      DEALLOCATE(du_nud)
  IF (ALLOCATED(dtheta_nud))  DEALLOCATE(dtheta_nud)
END IF

IF (l_retain_cld_tendencies) THEN
  IF (ALLOCATED(dt_cld))       DEALLOCATE(dt_cld)
END IF

IF (l_retain_stph_tendencies) THEN
  IF (ALLOCATED(dv_stph))      DEALLOCATE(dv_stph)
  IF (ALLOCATED(du_stph))      DEALLOCATE(du_stph)  
  IF (ALLOCATED(dtheta_stph))  DEALLOCATE(dtheta_stph)
END IF

IF (l_retain_bl_tendencies) THEN
  IF (ALLOCATED(dv_bl))        DEALLOCATE(du_bl)
  IF (ALLOCATED(du_bl))        DEALLOCATE(dv_bl)
  IF (ALLOCATED(dt_bl))        DEALLOCATE(dt_bl)
END IF

IF (l_retain_conv_all_tendencies) THEN
  IF (ALLOCATED(dv_conv))  DEALLOCATE(dv_conv)
  IF (ALLOCATED(du_conv))  DEALLOCATE(du_conv)
  IF (ALLOCATED(dq_conv))  DEALLOCATE(dq_conv)
  IF (ALLOCATED(dt_conv))  DEALLOCATE(dt_conv)
END IF

IF (l_retain_slow_tendencies) THEN
  IF (l_retain_ph1_tendencies) THEN
    IF (ALLOCATED(dtheta_ph1))   DEALLOCATE(dtheta_ph1)
  END IF

  IF (ALLOCATED(dv_gwd))   DEALLOCATE(dv_gwd)
  IF (ALLOCATED(du_gwd))   DEALLOCATE(du_gwd)
  IF (ALLOCATED(dt_gwd))   DEALLOCATE(dt_gwd)

  IF (ALLOCATED(dv_gwd))   DEALLOCATE(dv_gwd)
  IF (ALLOCATED(du_gwd))   DEALLOCATE(du_gwd)
  IF (ALLOCATED(dt_gwd))   DEALLOCATE(dt_gwd)

  IF (ALLOCATED(dq_mic))   DEALLOCATE(dq_mic)
  IF (ALLOCATED(dt_mic))   DEALLOCATE(dt_mic)

  IF (ALLOCATED(dq_lw))    DEALLOCATE(dq_lw)
  IF (ALLOCATED(dt_lw))    DEALLOCATE(dt_lw)

  IF (ALLOCATED(dq_sw))    DEALLOCATE(dq_sw)
  IF (ALLOCATED(dt_sw))    DEALLOCATE(dt_sw)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE destroy_phy_tend

END MODULE physics_tendencies_mod
