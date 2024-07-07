! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculates parcel vertical velocity diagnostic (w)

MODULE calc_w_eqn_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_W_EQN_MOD'

CONTAINS

SUBROUTINE calc_w_eqn ( k, npnts, np_full, nlev, bterm, blowst    &
                      , ekp14, ekp34, zkm1, zk, zkp12, zkp1       &
                      , thek,  thekp1,  thpk,  thpkp1             &
                      , qek,   qekp1,   qpk,   qpkp1              &
                      , qclek, qclekp1, qclpk, qclpkp1            &
                      , qcfek, qcfekp1, qcfpk, qcfpkp1            &
                      , w2p_km1, w2p_k, w2p_kp1 )

USE planet_constants_mod, ONLY: c_virtual, g

USE cv_run_mod, ONLY: cnv_wat_load_opt
USE cv_param_mod,                                               &
    ONLY: gamma_in_w_eqn, watload_opt, w2pi                     &
        , cumulus_r, drag_coeff, k2_const, fix_alpha, alpha_opt &
        , wCalcMethod, NegBuoyOpt, gamma_b, NegBuoyMinW

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Description:
!   Calculates parcel vertical velocity (w) based on
!   [Simpson & Wiggert 1969]. Equation is discretised so that
!   a forward difference is used for the first layer in the parcel
!   profile. The remaining layers used centred difference.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.0 programming standards.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

! Vector lengths and loop counters

INTEGER, INTENT(IN) :: k            ! Current model layer
INTEGER, INTENT(IN) :: npnts        ! Vector length (nconv)
INTEGER, INTENT(IN) :: np_full      ! Full field length (n_dp/n_md/n_sh)
INTEGER, INTENT(IN) :: nlev         ! Number of model levels

REAL, INTENT(IN)    :: ekp14(npnts) ! Entrainment rate on k+1/4 (-)
REAL, INTENT(IN)    :: ekp34(npnts) ! Entrainment rate on k+3/4 (-)


! Masks for points where ...
LOGICAL, INTENT(IN) :: bterm(npnts) ! ...parcels terminate in layer k+1
LOGICAL, INTENT(IN) :: blowst(npnts)! ...stability is low enough for
                                    !    convection to occur


! Layer variables
REAL, INTENT(IN) :: zkm1  (npnts) ! Height of theta level k-1
REAL, INTENT(IN) :: zk    (npnts) ! Height of theta level k
REAL, INTENT(IN) :: zkp12 (npnts) ! Height of theta level k+1/2
REAL, INTENT(IN) :: zkp1  (npnts) ! Height of theta level k+1

! Potential Temperature (Theta) on ...
REAL, INTENT(IN) :: thek  (npnts) ! ...layer centre k   (env)
REAL, INTENT(IN) :: thekp1(npnts) ! ...layer centre k+1 (env)
REAL, INTENT(IN) :: thpk  (npnts) ! ...layer centre k   (par)
REAL, INTENT(IN) :: thpkp1(npnts) ! ...layer centre k+1 (par)

! Specific humidity on ...
REAL, INTENT(IN) :: qek   (npnts) ! ...layer centre k   (env)
REAL, INTENT(IN) :: qekp1 (npnts) ! ...layer centre k+1 (env)
REAL, INTENT(IN) :: qpk   (npnts) ! ...layer centre k   (par)
REAL, INTENT(IN) :: qpkp1 (npnts) ! ...layer centre k+1 (par)

! Liquid water content on ...
REAL, INTENT(IN) :: qclek   (npnts) ! ...layer centre k   (env)
REAL, INTENT(IN) :: qclekp1 (npnts) ! ...layer centre k+1 (env)
REAL, INTENT(IN) :: qclpk   (npnts) ! ...layer centre k   (par)
REAL, INTENT(IN) :: qclpkp1 (npnts) ! ...layer centre k+1 (par)

! Ice water content on ...
REAL, INTENT(IN) :: qcfek   (npnts) ! ...layer centre k   (env)
REAL, INTENT(IN) :: qcfekp1 (npnts) ! ...layer centre k+1 (env)
REAL, INTENT(IN) :: qcfpk   (npnts) ! ...layer centre k   (par)
REAL, INTENT(IN) :: qcfpkp1 (npnts) ! ...layer centre k+1 (par)



! Output variables
! (Parcel vertical velocity)^2 on ...
REAL, INTENT(IN)    :: w2p_km1(npnts) ! ...layer centre k-1 [(m/s)^2]
REAL, INTENT(INOUT) :: w2p_k  (npnts) ! ...layer centre k   [(m/s)^2]
REAL, INTENT(OUT)   :: w2p_kp1(npnts) ! ...layer centre k+1 [(m/s)^2]

!=====================================================================

  ! Local variables
  !===================
REAL :: alpha  (npnts)
REAL :: coeff1
REAL :: dz
REAL :: thvpk
REAL :: thvek
REAL :: thvpkp1
REAL :: thvekp1
REAL :: term1
REAL :: term2
REAL :: r_min
REAL :: BuoyExcess

LOGICAL :: parcel_neg_buoy
INTEGER :: i


! Local parameters so as not to use magic numbers
! (development)
!===================================================================
!
! For method use to calculate w or w^2
INTEGER, PARAMETER :: method_simpsonwiggert = 1
INTEGER, PARAMETER :: method_buoyscale      = 2

! For setting calcuation of alpha
INTEGER, PARAMETER :: alpha_simpsonwiggert  = 1
INTEGER, PARAMETER :: alpha_um              = 2
INTEGER, PARAMETER :: alpha_const           = 3

! For treatment of negative buoyancy occurance
INTEGER, PARAMETER :: NegBuoyLimit = 1
INTEGER, PARAMETER :: NegBuoyZero  = 2

! For treatment of water loading
INTEGER, PARAMETER :: watload_simpsonwiggert = 1
INTEGER, PARAMETER :: watload_um             = 2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_W_EQN'

!---------------------------------------------------------------------
! End of declarations
!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! Current development notes (Subject to further testing)
  ! =======================================================
  ! The following notes apply to some sensitivity tests
  ! on using a fixing value of alpha in the SCM. These are
  ! NOT intended as general advice or recommendations
  ! =======================================================
  ! alpha_opt = um
  ! fix_alpha = 0.01   ! This alpha dies
  ! fix_alpha = 0.005  ! This alpha seems noise creeping in, effectively dead
  ! fix_alpha = 0.002  ! Works okay
  ! fix_alpha = 0.001  ! Nice and smooth some buoyancy violations
  ! fix_alpha = 0.0001 ! Smooth, similar to none, peaks higher up, possibly
                       ! too low

w2p_kp1(:) = 0.0

SELECT CASE (wCalcMethod)

CASE (method_simpsonwiggert)

  SELECT CASE (alpha_opt)
  CASE (alpha_SimpsonWiggert)
    alpha(:) = (0.75*k2_const+drag_coeff)*(3.0/(8.0*cumulus_r))

  CASE (alpha_um)
    ! Use convection scheme entrainment coefficents
    alpha(:) =  0.5*(  ekp14(:)/(zkp12(:) - zk(:)   ) &
                     + ekp34(:)/(zkp1(:)  - zkp12(:)) )
  CASE (alpha_const)
    alpha(:) = fix_alpha

  CASE DEFAULT
    alpha(:) = 0.0
  END SELECT

  coeff1 = 2.0*g/(1+gamma_in_w_eqn)

  ! Note: Could use module to get the layer thicknesses,
  !       though the height in the modules have halos etc, and these
  !       pointare compressed... not a problem if the indexes are available
  !       though how to relate the haloes.?
  !
  ! For now pass down height levels like the rest of the convection scheme.

  ! Give parcel at surface an initial kick to start it off?

  ! Calculate limit of Radius, to small and the velocity squared goes -ve
  ! which is not possible. i.e. the drag term (for the given seetings)
  ! is too large.

  ! Calculate min. radius required to prevent -ve w^2 for a given level

  DO i=1,npnts

    ! Calculate w2p if this point is marked as convecting
    IF ( .NOT. bterm(i) ) THEN

      ! Calculate Parcel/Environment Virtual Theta for level k
      !-------------------------------------------------------
      IF (cnv_wat_load_opt == 1) THEN

        ! How to apply water loading?
        IF (watload_opt == watload_simpsonwiggert) THEN
          thvpk = (thpk(i)*(1.0 + c_virtual*qpk(i)))*(1-qclpk(i)-qcfpk(i))
          thvek = (thek(i)*(1.0 + c_virtual*qek(i)))*(1-qclek(i)-qcfek(i))
        ELSE IF (watload_opt == watload_um) THEN
          thvpk =  thpk(i)*(1.0 + c_virtual*qpk(i)-qclpk(i)-qcfpk(i))
          thvek =  thek(i)*(1.0 + c_virtual*qek(i)-qclek(i)-qcfek(i))
        END IF

      ELSE
        ! Do not apply water loading
        thvpk = thpk(i) * (1.0 + c_virtual*qpk(i))
        thvek = thek(i) * (1.0 + c_virtual*qek(i))
      END IF


      ! Take into account of water loading in parcel
      ! Assume no water loading for now
      term2 = (thvpk/thvek) - 1.0


      ! This is the first level of a convective layer?
      ! If so, assume initial w=w2pi for this convecting layer
      IF (blowst(i)) w2p_k(i) = w2pi


      ! Calculate w2p on level k+1 using forward difference
      dz = zkp1(i) - zk(i)
      w2p_kp1(i) = dz*(coeff1*term2 - alpha(i)*w2p_k(i)) + w2p_k(i)

    END IF ! bterm test
  END DO ! npnts

CASE (method_buoyscale)

  ! Calculate w via scaled buoyancy on level k+1
  ! This can only be done here because the call to
  ! calc_w_eqn is at the end of convec2 when the input
  ! parcel properties from convec2 at level k+1 are now
  ! known

  !=======================================================
  ! NOTE: May need to add a check for the initiation of
  !       convection so as to apply scaling to initial
  !       buoyancy excess
  !=======================================================
  !
  DO i=1,npnts

    ! Calculate w2p if this point is marked as convecting
    IF ( .NOT. bterm(i) ) THEN

      ! Calculate Parcel/Environment Virtual Theta for level k+1
      !---------------------------------------------------------
      IF (cnv_wat_load_opt == 1) THEN
        thvpkp1 = thpkp1(i) * (1.0+c_virtual*qpkp1(i)-qclpkp1(i)-qcfpkp1(i))
        thvekp1 = thekp1(i) * (1.0+c_virtual*qekp1(i)-qclekp1(i)-qcfekp1(i))
      ELSE
        thvpkp1 = thpkp1(i) * (1.0+c_virtual*qpkp1(i))
        thvekp1 = thekp1(i) * (1.0+c_virtual*qekp1(i))
      END IF

      BuoyExcess = g * (thvpkp1 - thvekp1) / thvekp1

      w2p_kp1(i) = (gamma_b*BuoyExcess) * (gamma_b*BuoyExcess)

      IF (BuoyExcess < 0.0) THEN
        SELECT CASE(NegBuoyOpt)
        CASE (NegBuoyLimit)
          w2p_kp1(i) = (gamma_b*NegBuoyMinW) * (gamma_b*NegBuoyMinW)
        CASE (NegBuoyZero)
          w2p_kp1(i) = 0.0
        END SELECT
      END IF

    END IF
  END DO

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_w_eqn

END MODULE calc_w_eqn_mod
