! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_conserv_moist_mod
IMPLICIT NONE

! Description:
!            This routine enforces the mass conservation for
!            the 6 moisture variables. It should work with
!            either mixing ratios and dry density or with
!            specific humdities and wet density.
!            If L_conserve_mass=.FALSE., then the routine make
!            sure that the moisture variables are above the
!            minimum values required only (doesn't impose
!            any mass conservation constraint).
!
! Method: ENDGame formulation version 3.02
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code description:
! Language: Fortran 90.
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CONSERV_MOIST_MOD'

CONTAINS
SUBROUTINE eg_correct_moisture_fix(rho_n, rho_np1,   &
               q1_n, q2_n, q3_n, q4_n, q5_n, q6_n,   &
   q1_np1, q2_np1, q3_np1, q4_np1, q5_np1, q6_np1,   &
              q1_s, q2_s, q3_s, q4_s, q5_s, q6_s,    &
    q1min, L_q4, L_q5, L_q6, L_conserve_mass, mype,  &
    pseudo_lbflux, number_qs_in,                     &
                      L_conserv_smooth_lap           )

USE eg_check_conserv_mod
USE eg_mass_conserv_mod
USE umPrintMgr
USE Field_Types
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims,        &
                                 tdims_s, array_dims
USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_MOISTURE_FIX'

INTEGER, INTENT(IN) ::  mype
LOGICAL, INTENT(IN) ::  L_conserv_smooth_lap

REAL, INTENT(IN)   ::   rho_n(pdims_s%i_start:pdims_s%i_end,   &
                              pdims_s%j_start:pdims_s%j_end,   &
                              pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN)   :: rho_np1(pdims_s%i_start:pdims_s%i_end,   &
                              pdims_s%j_start:pdims_s%j_end,   &
                              pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN)   ::    q1_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q2_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q3_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q4_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN)   ::    q5_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q6_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q1_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q2_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q3_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q4_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q5_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::     q6_s (tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(INOUT) :: q1_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q2_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q3_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q4_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q5_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q6_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(IN)     :: q1min
LOGICAL, INTENT(IN) :: L_q4, L_q5, L_q6, L_conserve_mass
INTEGER, INTENT(IN) :: number_qs_in
REAL, INTENT(IN)    :: pseudo_lbflux(number_qs_in)


! locals

INTEGER :: number_qs, i, j, k, ie, je, pt_q4, pt_q5, pt_q6

REAL,    ALLOCATABLE ::   qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE ::   qs_s(:,:,:,:)
REAL,    ALLOCATABLE ::  qsmin(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( L_conserve_mass ) THEN

  number_qs = 3
  IF ( L_q4  ) THEN
    number_qs = number_qs + 1
    pt_q4     = number_qs
  END IF
  IF ( L_q5 ) THEN
    number_qs = number_qs + 1
    pt_q5     = number_qs
  END IF
  IF ( L_q6   ) THEN
    number_qs = number_qs + 1
    pt_q6     = number_qs
  END IF

  ALLOCATE(    qs_n(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(  qs_np1(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(    qs_s(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(     qsmin(number_qs) )

  ! copy the moisture variables and sources into 4D superarrays

  ! Copying in halos removes need to halo exchange for qs_np1
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        qs_np1(i,j,k,1) = q1_np1(i,j,k)
        qs_np1(i,j,k,2) = q2_np1(i,j,k)
        qs_np1(i,j,k,3) = q3_np1(i,j,k)


        IF( l_q4 ) THEN
           qs_np1(i,j,k,pt_q4) = q4_np1(i,j,k)
        END IF

        IF( l_q5 ) THEN
          qs_np1(i,j,k,pt_q5) = q5_np1(i,j,k)
        END IF

        IF( l_q6 ) THEN
          qs_np1(i,j,k,pt_q6) = q6_np1(i,j,k)
        END IF
      END DO
    END DO
  END DO

  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,1) = q1_n(i,j,k)
        qs_n(i,j,k,2) = q2_n(i,j,k)
        qs_n(i,j,k,3) = q3_n(i,j,k)

        qs_s(i,j,k,1) = q1_s(i,j,k)
        qs_s(i,j,k,2) = q2_s(i,j,k)
        qs_s(i,j,k,3) = q3_s(i,j,k)


        IF( l_q4 ) THEN
           qs_n(i,j,k,pt_q4) = q4_n(i,j,k)
           qs_s(i,j,k,pt_q4) = q4_s(i,j,k)
        END IF

        IF( l_q5 ) THEN
          qs_n(i,j,k,pt_q5) = q5_n(i,j,k)
          qs_s(i,j,k,pt_q5) = q5_s(i,j,k)
        END IF

        IF( l_q6 ) THEN
          qs_n(i,j,k,pt_q6) = q6_n(i,j,k)
          qs_s(i,j,k,pt_q6) = q6_s(i,j,k)
        END IF
      END DO
    END DO
  END DO
  ! set min-values for each tracer

  qsmin(1)           = q1min
  qsmin(2:number_qs) = 0.0

  ! print mass conservation error before correction if needed

  IF (printstatus > prstatus_normal) THEN
    IF ( mype == 0) THEN
      CALL umPrint( 'Error in mass conservation for' //        &
           ' moisture before correction',src='eg_correct_moisture')
    END IF
    CALL eg_check_mass_conservation_fix(rho_n, rho_np1, qs_n,  &
                                    qs_np1, qs_s,              &
                      number_qs, mype,'eg_correct_moisture_fix1')
  END IF

  ! enfore mass conservation

  CALL eg_mass_conservation_fix(rho_n, rho_np1, qs_n, qs_np1, qs_s, &
                            qsmin, number_qs,                       &
                            L_conserv_smooth_lap, pseudo_lbflux     )

  ! may be there is no need for this swap-bounds if the moisture
  ! fileds are swap-bounded later on -- but it's safer to leave it

  i  = tdims%i_len
  j  = tdims%j_len
  k  = tdims_s%k_len*number_qs
  ie = tdims%i_start - tdims_s%i_start
  je = tdims%j_start - tdims_s%j_start

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds(qs_np1,i,j,k,ie,je,fld_type_p,swap_field_is_scalar)

  ! put the corrected moisture variables back into their original arrays
  ! unpack them from the superarray

  DO k = tdims_s%k_start, tdims_s%k_end
    q1_np1(:,:,k) = qs_np1(:,:,k,1)
    q2_np1(:,:,k) = qs_np1(:,:,k,2)
    q3_np1(:,:,k) = qs_np1(:,:,k,3)
  END DO
  IF ( L_q4 ) THEN
    DO k = tdims_s%k_start, tdims_s%k_end
      q4_np1(:,:,k) = qs_np1(:,:,k,pt_q4)
    END DO
  END IF
  IF ( L_q5 ) THEN
    DO k = tdims_s%k_start, tdims_s%k_end
      q5_np1(:,:,k) = qs_np1(:,:,k,pt_q5)
    END DO
  END IF
  IF ( L_q6 ) THEN
    DO k = tdims_s%k_start, tdims_s%k_end
      q6_np1(:,:,k) = qs_np1(:,:,k,pt_q6)
    END DO
  END IF

  ! print mass conservation error after correction if needed

  IF (printstatus > prstatus_normal) THEN
    IF ( mype == 0) THEN
      CALL umPrint('Error in mass conservation for' //         &
           ' moisture after correction',src='eg_correct_moisture')
    END IF
    CALL eg_check_mass_conservation_fix(rho_n, rho_np1, qs_n,      &
                                    qs_np1, qs_s,                  &
                          number_qs, mype,'eg_correct_moisture_fix2')
  END IF


  DEALLOCATE (qsmin)
  DEALLOCATE (qs_s)
  DEALLOCATE (qs_np1)
  DEALLOCATE (qs_n)

ELSE

  ! if L_conserve_mass==.false. then simply make sure that the moisture
  ! variables are above the minimum values required

  q1_np1 = MAX(q1_np1, q1min)
  q2_np1 = MAX(q2_np1, 0.0  )
  q3_np1 = MAX(q3_np1, 0.0  )
  IF( L_q4 )  q4_np1 = MAX(q4_np1, 0.0  )
  IF( L_q5 )  q5_np1 = MAX(q5_np1, 0.0  )
  IF( L_q6 )  q6_np1 = MAX(q6_np1, 0.0  )

END IF  ! L_conserve_mass


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_correct_moisture_fix

!=====================================================================
!
! Version to be used if l_fix_conserv = .FALSE.
!
SUBROUTINE eg_correct_moisture(rho_n, rho_np1,       &
               q1_n, q2_n, q3_n, q4_n, q5_n, q6_n,   &
   q1_np1, q2_np1, q3_np1, q4_np1, q5_np1, q6_np1,   &
              q1_s, q2_s, q3_s, q4_s, q5_s, q6_s,    &
    q1min, L_q4, L_q5, L_q6, L_conserve_mass, mype,  &
   pseudo_lbflux, number_qs_in,                      &
                      L_conserv_smooth_lap           )

USE eg_check_conserv_mod
USE eg_mass_conserv_mod
USE umPrintMgr
USE Field_Types
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims,        &
                                 tdims_s, array_dims

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
IMPLICIT NONE
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CORRECT_MOISTURE'

INTEGER, INTENT(IN) ::  mype
LOGICAL, INTENT(IN) :: L_conserv_smooth_lap

REAL, INTENT(IN)   ::   rho_n(pdims_s%i_start:pdims_s%i_end,   &
                              pdims_s%j_start:pdims_s%j_end,   &
                              pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN)   :: rho_np1(pdims_s%i_start:pdims_s%i_end,   &
                              pdims_s%j_start:pdims_s%j_end,   &
                              pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(IN)   ::    q1_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q2_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q3_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q4_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN)   ::    q5_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q6_n(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::      q1_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q2_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q3_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q4_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::      q5_s(tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)
REAL, INTENT(IN) ::     q6_s (tdims%i_start:tdims%i_end,       &
                              tdims%j_start:tdims%j_end,       &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(INOUT) :: q1_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q2_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q3_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q4_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q5_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(INOUT) :: q6_np1(tdims_s%i_start:tdims_s%i_end,   &
                              tdims_s%j_start:tdims_s%j_end,   &
                              tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(IN)     :: q1min
LOGICAL, INTENT(IN) ::  L_q4, L_q5, L_q6, L_conserve_mass
INTEGER, INTENT(IN) :: number_qs_in
REAL, INTENT(IN)    :: pseudo_lbflux(number_qs_in)

! locals

INTEGER :: number_qs, k,i,j
INTEGER :: pt_q4, pt_q5, pt_q6

REAL,    ALLOCATABLE ::   qs_n(:,:,:,:)
REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
REAL,    ALLOCATABLE ::   qs_s(:,:,:,:)

REAL,    ALLOCATABLE ::     qsmin(:)
INTEGER, ALLOCATABLE :: qswitches(:)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

number_qs = 3
IF ( L_q4  ) THEN
  number_qs = number_qs + 1
  pt_q4  = number_qs
END IF

IF ( L_q5 ) THEN
  number_qs = number_qs + 1
  pt_q5  = number_qs
END IF

IF ( L_q6   ) THEN
  number_qs = number_qs + 1
  pt_q6  = number_qs
END IF

IF ( L_conserve_mass ) THEN

  ALLOCATE(    qs_n(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(  qs_np1(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(    qs_s(tdims_s%i_start:tdims_s%i_end,             &
                    tdims_s%j_start:tdims_s%j_end,             &
                    tdims_s%k_start:tdims_s%k_end,number_qs))

  ALLOCATE(     qsmin(number_qs) )
  ALLOCATE( qswitches(number_qs) )

  ! copy the 6 moisture variable fields and sources into
  ! superarrays of dimension 6

  ! Copying in halos removes need to halo exchange for qs_np1
  DO k = tdims_s%k_start, tdims_s%k_end
    DO j = tdims_s%j_start, tdims_s%j_end
      DO i = tdims_s%i_start, tdims_s%i_end
        qs_np1(i,j,k,1) = q1_np1(i,j,k)
        qs_np1(i,j,k,2) = q2_np1(i,j,k)
        qs_np1(i,j,k,3) = q3_np1(i,j,k)


        IF( l_q4 ) THEN
           qs_np1(i,j,k,pt_q4) = q4_np1(i,j,k)
        END IF

        IF( l_q5 ) THEN
          qs_np1(i,j,k,pt_q5) = q5_np1(i,j,k)
        END IF

        IF( l_q6 ) THEN
          qs_np1(i,j,k,pt_q6) = q6_np1(i,j,k)
        END IF
      END DO
    END DO
  END DO

  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qs_n(i,j,k,1) = q1_n(i,j,k)
        qs_n(i,j,k,2) = q2_n(i,j,k)
        qs_n(i,j,k,3) = q3_n(i,j,k)

        qs_s(i,j,k,1) = q1_s(i,j,k)
        qs_s(i,j,k,2) = q2_s(i,j,k)
        qs_s(i,j,k,3) = q3_s(i,j,k)


        IF( l_q4 ) THEN
           qs_n(i,j,k,pt_q4) = q4_n(i,j,k)
           qs_s(i,j,k,pt_q4) = q4_s(i,j,k)
        END IF

        IF( l_q5 ) THEN
          qs_n(i,j,k,pt_q5) = q5_n(i,j,k)
          qs_s(i,j,k,pt_q5) = q5_s(i,j,k)
        END IF

        IF( l_q6 ) THEN
          qs_n(i,j,k,pt_q6) = q6_n(i,j,k)
          qs_s(i,j,k,pt_q6) = q6_s(i,j,k)
        END IF
      END DO
    END DO
  END DO

  ! Make sure the boundary condition filed(surface)=field(level 1) is imposed.
  ! This may be not necessary here because it may already be done higher up
  ! but for safety we keep it.

  qs_n(:,:,tdims_s%k_start,:) =   qs_n(:,:,tdims_s%k_start+1,:)
  qs_s(:,:,tdims_s%k_start,:) =   qs_s(:,:,tdims_s%k_start+1,:)
  qs_np1(:,:,tdims_s%k_start,:) = qs_np1(:,:,tdims_s%k_start+1,:)

  ! Set switches to 1 (means impose mass conservation). Usually the first
  ! 3 moisture variables are always there but the last 3 are optional.

  qswitches(:)   = 0
  qswitches(1:3) = 1

  IF ( L_q4 )  qswitches(pt_q4) = 1
  IF ( L_q5 )  qswitches(pt_q5) = 1
  IF ( L_q6 )  qswitches(pt_q6) = 1

  qsmin(1)           = q1min
  qsmin(2:number_qs) = 0.0

  ! print mass conservation error before correction if needed

  IF (printstatus > prstatus_normal) THEN
    IF ( mype == 0) THEN
      CALL umPrint( 'Error in mass conservation for' //        &
          ' moisture before correction',src='eg_correct_moisture')
    END IF
    CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,      &
                                    qs_np1, qs_s,              &
               qswitches, number_qs, mype,'eg_correct_moisture')
  END IF

  ! enfore mass conservation

  CALL eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1, qs_s, &
                            qsmin, qswitches, number_qs,        &
                            L_conserv_smooth_lap, pseudo_lbflux )

  ! put the corrected moisture variables back into their original arrays
  ! unpack them from the superarray

  DO k = tdims_s%k_start, tdims_s%k_end
    q1_np1(:,:,k) = qs_np1(:,:,k,1)
    q2_np1(:,:,k) = qs_np1(:,:,k,2)
    q3_np1(:,:,k) = qs_np1(:,:,k,3)
    IF( l_q4 ) q4_np1(:,:,k) = qs_np1(:,:,k,pt_q4)
    IF( l_q5 ) q5_np1(:,:,k) = qs_np1(:,:,k,pt_q5)
    IF( l_q6 ) q6_np1(:,:,k) = qs_np1(:,:,k,pt_q6)
  END DO

  ! print mass conservation error after correction if needed

  IF (printstatus > prstatus_normal) THEN
    IF ( mype == 0) THEN
      CALL umPrint('Error in mass conservation for' //         &
          ' moisture after correction',src='eg_correct_moisture')
    END IF
    CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,      &
                                    qs_np1, qs_s,              &
               qswitches, number_qs, mype,'eg_correct_moisture')
  END IF


  DEALLOCATE (qswitches)
  DEALLOCATE (qsmin)
  DEALLOCATE (qs_s)
  DEALLOCATE (qs_np1)
  DEALLOCATE (qs_n)

ELSE

  ! if L_conserve_mass==.false. then simply make sure that the moisture
  ! variables are above the minimum values required

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tdims_s, q1_np1,q2_np1,q3_np1,q4_np1,q5_np1,q6_np1,q1min,   &
!$OMP&        l_q4, l_q5, l_q6)
  DO k = tdims_s%k_start,tdims_s%k_end
    DO j = tdims_s%j_start,tdims_s%j_end
      DO i = tdims_s%i_start,tdims_s%i_end
        q1_np1(i,j,k) = MAX(q1_np1(i,j,k), q1min)
        q2_np1(i,j,k) = MAX(q2_np1(i,j,k), 0.0  )
        q3_np1(i,j,k) = MAX(q3_np1(i,j,k), 0.0  )
      END DO
    END DO

    IF( l_q4 ) THEN
      DO j = tdims_s%j_start,tdims_s%j_end
        DO i = tdims_s%i_start,tdims_s%i_end
          q4_np1(i,j,k) = MAX(q4_np1(i,j,k), 0.0  )
        END DO
      END DO
    END IF

    IF( l_q5 ) THEN
      DO j = tdims_s%j_start,tdims_s%j_end
        DO i = tdims_s%i_start,tdims_s%i_end
          q5_np1(i,j,k) = MAX(q5_np1(i,j,k), 0.0  )
        END DO
      END DO
    END IF

    IF( l_q6 ) THEN
      DO j = tdims_s%j_start,tdims_s%j_end
        DO i = tdims_s%i_start,tdims_s%i_end
          q6_np1(i,j,k) = MAX(q6_np1(i,j,k), 0.0  )
        END DO
      END DO
    END IF

  END DO
!$OMP END PARALLEL DO

END IF  ! L_conserve_mass

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_correct_moisture

END MODULE eg_conserv_moist_mod
