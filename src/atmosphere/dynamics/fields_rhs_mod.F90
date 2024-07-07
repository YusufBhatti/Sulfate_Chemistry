MODULE fields_rhs_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame departure points
!
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

IMPLICIT NONE

! Horizontal coordinates


REAL,SAVE, ALLOCATABLE ::                                             &
   r_u      (:,:,:),&
   r_v      (:,:,:),&
   r_w      (:,:,:),&
   r_u_p2   (:,:,:),&
   r_v_p2   (:,:,:),&
   r_w_p2   (:,:,:),&
   r_u_p2_n (:,:,:),&
   r_v_p2_n (:,:,:),&
   r_w_p2_n (:,:,:),&
   r_theta  (:,:,:),&
   r_rho    (:,:,:),&
   r_m_v    (:,:,:),&
   r_m_cl   (:,:,:),&
   r_m_cf   (:,:,:),&
   r_m_r    (:,:,:),&
   r_m_gr   (:,:,:),&
   r_m_cf2  (:,:,:),&
   r_u_d    (:,:,:),&
   r_v_d    (:,:,:),&
   r_w_d    (:,:,:),&
   r_theta_d(:,:,:),&
   r_rho_d  (:,:,:),&
   r_m_v_d  (:,:,:),&
   r_m_cl_d (:,:,:),&
   r_m_cf_d (:,:,:),&
   r_m_r_d  (:,:,:),&
   r_m_gr_d (:,:,:),&
   r_m_cf2_d(:,:,:),&
   s_u      (:,:,:),&
   s_v      (:,:,:),&
   s_w      (:,:,:),&
   s_thetav (:,:,:),&
   s_m_v    (:,:,:),&
   s_m_cl   (:,:,:),&
   s_m_cf   (:,:,:),&
   s_m_r    (:,:,:),&
   s_m_gr   (:,:,:),&
   s_m_cf2  (:,:,:),&
   r_u_skeb (:,:,:),&
   r_v_skeb (:,:,:),&
 ! Add SPT increments
   theta_spt(:,:,:),&
   q_spt   (:,:,:),&
   r_u_spt (:,:,:),&
   r_v_spt (:,:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FIELDS_RHS_MOD'

CONTAINS

SUBROUTINE init_r_skeb()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod, ONLY : tdims_s, vdims_s, udims_s, wdims, pdims_s

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_R_SKEB'

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (    r_u_skeb(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_r_skeb",ierr, "Unable to allocate1.")

ALLOCATE (    r_v_skeb(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_r_skeb",ierr, "Unable to allocate2.")

r_u_skeb = 0.0
r_v_skeb = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_r_skeb

SUBROUTINE destroy_r_skeb()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DESTROY_R_SKEB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE ( r_v_skeb )
DEALLOCATE ( r_u_skeb )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE destroy_r_skeb

SUBROUTINE init_r_spt()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_R_SPT'

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (    theta_spt(tdims_s%i_start:tdims_s%i_end,                     &
                        tdims_s%j_start:tdims_s%j_end,                     &
                        tdims_s%k_start:tdims_s%k_end),stat=ierr)


IF (ierr/=0) CALL Ereport("init_r_spt",ierr, "Unable to allocate theta_spt.")

ALLOCATE (    q_spt (tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:udims_s%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_r_spt",ierr, "Unable to allocate q_spt.")

ALLOCATE (    r_u_spt (udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_r_spt",ierr, "Unable to allocate r_u_spt.")

ALLOCATE (    r_v_spt (vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_r_spt",ierr, "Unable to allocate r_v_spt.")

theta_spt = 0.0
q_spt     = 0.0
r_u_spt   = 0.0
r_v_spt   = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_r_spt

SUBROUTINE destroy_r_spt()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod, ONLY : tdims

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DESTROY_R_SPT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE ( theta_spt )
DEALLOCATE ( q_spt )
DEALLOCATE ( r_v_spt )
DEALLOCATE ( r_u_spt )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE destroy_r_spt

SUBROUTINE init_fields_rhs(l_skeb2,l_spt)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_FIELDS_RHS'

LOGICAL, INTENT (IN) :: l_skeb2,l_spt

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (         r_u(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                 r_u_d(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                r_u_p2(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
              r_u_p2_n(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                   s_u(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                   r_v(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                 r_v_d(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                r_v_p2(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
              r_v_p2_n(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                   s_v(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                   r_w(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                 r_w_d(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                r_w_p2(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
              r_w_p2_n(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                   s_w(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
               r_theta(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
             r_theta_d(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
              s_thetav(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
                 r_rho(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
               r_rho_d(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
                 r_m_v(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                r_m_cl(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                r_m_cf(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
               r_m_cf2(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                r_m_gr(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                 r_m_r(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
               r_m_v_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
              r_m_cl_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
              r_m_cf_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
             r_m_cf2_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
              r_m_gr_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
               r_m_r_d(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                 s_m_v(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                s_m_cl(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                s_m_cf(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
               s_m_cf2(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                s_m_gr(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),                        &
                 s_m_r(tdims%i_start:tdims%i_end,                         &
                       tdims%j_start:tdims%j_end,                         &
                       tdims%k_start:tdims%k_end),stat=ierr)

IF (ierr/=0) CALL Ereport("init_fields_rhs",ierr, "Unable to allocate.")

IF (l_skeb2) CALL init_r_skeb()

IF (l_spt) CALL init_r_spt()

! These fields may not be initialised with certain physics switches
! but are then used
! in eg_helm_rhs_star (for example) without the flags being tested!
! Therefore we do initialise them here to zero:

s_m_r  (:,:,:) = 0.0
s_m_gr (:,:,:) = 0.0
s_m_cf2(:,:,:) = 0.0
s_m_cf (:,:,:) = 0.0
s_m_cl (:,:,:) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_fields_rhs
END MODULE fields_rhs_mod
