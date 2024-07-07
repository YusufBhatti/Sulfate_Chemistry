! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates theta and q of detraining air in layer k
!
MODULE thetar_6a_mod

IMPLICIT NONE

!
! Description:
!   Calculates the potential temperature and the humidity of the parcel
!   undergoing forced detrainment in layer k
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='THETAR_6A_MOD'

CONTAINS

SUBROUTINE thetar_6a(npnts, exk, exkp1, pk, pkp1, thek, thekp1, qek, qekp1, &
                     qpk, watldek, watldekp1, watldpk, watldpkp1,           &
                     bwk, bgmk, thrk, qrk)

USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls
USE planet_constants_mod, ONLY: c_virtual, rv, recip_kappa, pref
USE cv_run_mod, ONLY: fdet_opt

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    l_new_qsat_conv !Currently defaults to FALSE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: npnts        ! Number of points

REAL,INTENT(IN) :: exk(npnts)      ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: exkp1(npnts)    ! Exner ratio at mid-point of layer k+1
REAL,INTENT(IN) :: pk(npnts)       ! pressure at mid-point of layer k (Pa)
REAL,INTENT(IN) :: pkp1(npnts)     ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(IN) :: thek(npnts)     ! Env. p. temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)   ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qek(npnts)      ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)    ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qpk(npnts)      ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: watldek(npnts)  ! Env. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldekp1(npnts)! Env. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldpk(npnts)  ! Par. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldpkp1(npnts)! Par. water loading in layer k+1 (kg/kg)

LOGICAL,INTENT(IN) :: bwk(npnts)   ! Mask for whether condensate is
                                   ! liquid in layer k+1
LOGICAL,INTENT(IN) :: bgmk(npnts)  ! Mask for parcels which are
                                   ! saturated in layer k

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: thrk(npnts)    ! Pot. temperature of forced detrained
                                   ! parcel in layer k (K)
REAL,INTENT(OUT) :: qrk(npnts)     ! Specific humidity of forced detrained
                                   ! parcel in layer k (kg/kg)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i, j ! loop counters

REAL :: trk(npnts)      ! estimate of the detrained parcel's temp (K)
REAL :: thve(npnts)     ! env. thetav that sets the detrained parcel's thv (K)
REAL :: qsrk(npnts)     ! qsat at trk
REAL :: watldr(npnts)   ! Par. water loading of detrained parcel (kg/kg)
REAL :: exr(npnts)      ! Exner ratio for the detrained parcel
REAL :: pr(npnts)       ! Pressure of the detrained parcel (Pa)
REAL :: dqsdth          ! Rate of change of qsat with potential temperature
REAL :: el              ! Latent heat of gas-to-whatever-condenses
                        ! PC2 defn (J/kg)

REAL, PARAMETER :: a_smth = 0.5 ! The weighting between at k and k-1 for 
                                ! fdet_opt=2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='THETAR_6A'


!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (fdet_opt == 3) THEN
  DO i=1,npnts
    ! Set the env. properties to the weighted average of k and k+1
    thve(i)   = a_smth *(thek(i)  *(1.0+c_virtual*qek(i)  - watldek(i)))  &
            +(1-a_smth)*(thekp1(i)*(1.0+c_virtual*qekp1(i)- watldekp1(i)))
    exr(i)    = a_smth*exk(i)     +(1.0-a_smth)*exkp1(i)
    ! Set the water loading of the detrained parcel to the weighted
    ! average of k and k+1
    watldr(i) = a_smth*watldpk(i) +(1.0-a_smth)*watldpkp1(i)

    IF ( bgmk(i) ) THEN
      ! If saturated then an iterative calculation is required
      ! Set the first estimate of the detrained p.temp to that of the
      ! environment
      thrk(i)   = a_smth*thek(i)  +(1.0-a_smth)*thekp1(i)
    ELSE
      ! If not saturated the detrained p.temp and humidity can be directly
      ! calculated
      qrk(i)  = 0.5*(qpk(i) + a_smth * qek(i) +(1.0-a_smth)*qekp1(i))
      thrk(i) = thve(i) / (1.0 + c_virtual*qrk(i) - watldr(i))      
    END IF
    !update the parcel's temperature
    trk(i)  = thrk(i) * exk(i)
  END DO

  DO i=1,npnts
    ! calculate pressure at forced detrainment to be consistent with exner
    pr(i)     = pref*exr(i)**recip_kappa
  END DO

  IF ( l_new_qsat_conv ) THEN
    CALL qsat_new(qsrk, trk, pr, npnts)
  ELSE
    ! DEPENDS ON: qsat
    CALL qsat (qsrk, trk, pr, npnts)
  END IF

  DO j=1,3  !Three iterations should be sufficient
    DO i=1,npnts
      IF ( bgmk(i) ) THEN !if saturated
        IF ( bwk(i) ) THEN
          el=lc
        ELSE
          el=ls
        END IF

        !Calculate the gradient of qsat. nb qsrk is saturated humidity at trk
        dqsdth = el * qsrk(i) / ( rv * exr(i) * thrk(i) * thrk(i) )

        !Calculate the next estimate of the parcel's p.temp
        thrk(i) = ( thve(i)                                                 &
                  - thrk(i)*c_virtual*(qsrk(i) - thrk(i)*dqsdth) )          &
                  / ( 1.0 + c_virtual*thrk(i)*dqsdth - watldr(i) )

        !update the parcel's temperature
        trk(i)  = thrk(i) * exr(i)
      END IF    !if saturated
    END DO      !i

    ! Assume that the detrained parcel is saturated and update its humidity.
    IF ( l_new_qsat_conv ) THEN
      CALL qsat_new(qsrk, trk, pr, npnts)
    ELSE
      ! DEPENDS ON: qsat
      CALL qsat (qsrk, trk, pr, npnts)
    END IF

  END DO  !j
ELSE ! fdet_opt != 3
  DO i=1,npnts
    IF ( bgmk(i) ) THEN
      ! If saturated then an iterative calculation is required
      ! Set the first estimate of the detrained p.temp to that of the
      ! environment
      thrk(i) = thek(i)
    ELSE
      ! If not saturated the detrained p.temp and humidity can be directly
      ! calculated
      IF (fdet_opt == 2) THEN
        thrk(i) = thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))             &
                  /(1.0 + c_virtual*0.5*(qpk(i) + qek(i)) - watldpk(i))
        qrk(i)  = 0.5*(qpk(i) + qek(i))
      ELSE
        thrk(i) = thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))             &
                  /(1.0 + c_virtual*qpk(i) - watldpk(i))
        qrk(i)  = qpk(i)
      END IF
    END IF
    !update the parcel's temperature
    trk(i)  = thrk(i) * exk(i)
  END DO

  IF ( l_new_qsat_conv ) THEN
    CALL qsat_new(qsrk, trk, pk, npnts)
  ELSE
    ! DEPENDS ON: qsat
    CALL qsat (qsrk, trk, pk, npnts)
  END IF


  DO j=1,3  !Three iterations should be sufficient
    DO i=1,npnts
      IF ( bgmk(i) ) THEN !if saturated
        IF ( bwk(i) ) THEN
          el=lc
        ELSE
          el=ls
        END IF

        !Calculate the gradient of qsat. nb qsrk is saturated humidity at trk
        dqsdth = el * qsrk(i) / ( rv * exk(i) * thrk(i) * thrk(i) )

        !Calculate the next estimate of the parcel's p.temp
        thrk(i) = ( thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))           &
                  - thrk(i)*c_virtual*(qsrk(i) - thrk(i)*dqsdth) )          &
                  / ( 1.0 + c_virtual*thrk(i)*dqsdth - watldpk(i) )

        !update the parcel's temperature
        trk(i)  = thrk(i) * exk(i)
      END IF    !if saturated
    END DO      !i

    ! Assume that the detrained parcel is saturated and update its humidity.
    IF ( l_new_qsat_conv ) THEN
      CALL qsat_new(qsrk, trk, pk, npnts)
    ELSE
      ! DEPENDS ON: qsat
      CALL qsat (qsrk, trk, pk, npnts)
    END IF

  END DO  !j
END IF

DO i=1,npnts
  IF ( bgmk(i) ) THEN
    ! If saturated update qrk to saturated humidity at trk
    qrk(i) = qsrk(i)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE thetar_6a
END MODULE thetar_6a_mod
