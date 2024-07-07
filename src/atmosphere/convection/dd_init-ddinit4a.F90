! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to initialise the downdraught
!
MODULE dd_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DD_INIT_MOD'
CONTAINS

SUBROUTINE dd_init(npnts, np_full                                           &
                  ,th_ud_k, q_ud_k, the_k, qe_k, qse_k, pk, exk             &
                  ,thdd_k, qdd_k, deltd, delqd                              &
                  ,bdd_start, k, bddi, bdd_on                               &
                  ,l_tracer                                                 &
                  ,ntra, tra_ud_k, trae_k, tradd_k, deltrad)


USE planet_constants_mod, ONLY: c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE satcal_mod, ONLY: satcal
IMPLICIT NONE

!
! Description: Routine to initialise the downdraught
!
! Method: UM documentataion paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,ntra                 &  ! Number of tracers
 ,k                       ! Present model layer

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers

LOGICAL, INTENT(IN) :: &
  bddi(npnts)          & ! Mask for points where downdraught may initiate
 ,bdd_on(npnts)          ! mask for those points where downdraught is on

REAL, INTENT(IN) ::  &
  th_ud_k(npnts)     & ! Parcel potential temperature of updraught, layer k (K)
 ,q_ud_k(npnts)      & ! Parcel mixing ratio of updraught, layer k (kg/kg)
 ,the_k(npnts)       & ! Potential temperature of environment in layer k (K)
 ,qe_k(npnts)        & ! Mixing ratio of environment in layer k (kg/kg)
 ,qse_k(npnts)       & ! Mixing ratio of environment qsat in layer k (kg/kg)
 ,pk(npnts)          & ! Pressure of layer k (Pa)
 ,exk(npnts)           ! Exner ratio layer k

REAL, INTENT(IN) ::     &
  trae_k(np_full,ntra)  & ! Tracer content of environment in layer k  (kg/kg)
 ,tra_ud_k(np_full,ntra)  ! Parcel tracer content of updraught, layer k (kg/kg)

LOGICAL, INTENT(INOUT) ::  &
  bdd_start(npnts)           ! input mask for those points where DD may start ?
                             ! output set to true if DD started on level k

REAL, INTENT(OUT) ::     &
  thdd_k(npnts)          & ! Downdraught potential temperature of layer k (K)
 ,qdd_k(npnts)           & ! Downdraught mixing ratio of layer k (kg/kg)
 ,tradd_k(np_full,ntra)  & ! Downdraught tracer content of layer k (kg/kg)
 ,deltd(npnts)           & ! Cooling necessary to achieve saturation
 ,delqd(npnts)           & ! Moistening necessary to achieve saturation
 ,deltrad(np_full,ntra)    ! Depletion of environment tracer due to formation
                           ! of Downdraught
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::       &
  i, ktra          ! Loop counters

INTEGER ::       &
  errorstatus

CHARACTER (LEN=errormessagelength) ::  cmessage

REAL ::               &
  th_mean(npnts)      & ! Mean potential temperature used in calculation of
                        ! saturated downdraught potential temperature  in
                        ! layer k (K)
 ,q_mean(npnts)       & ! Mean mixing ratio used in calculation of
                        ! saturated downdraught mixing ratio in layer k (kg/kg)
 ,t_mean(npnts)       & ! Mean temperature used in calculation of
                        ! saturated downdraught potential temperature in
                        ! layer k (kg/kg)
 ,tra_mean(npnts,ntra)& ! Mean tracer used as initial tracer content of
                        ! downdraught in layer k (kg/kg)
 ,thdds(npnts)        & ! Saturated downdraught potential temperature in
                        ! layer k (K)
 ,qdds(npnts)         & ! Saturated downdraught mixing ratio in
                        ! layer k (kg/kg)
 ,buoy(npnts)           ! Buoyancy of parcel in layer k

REAL ::             &
  thdd_v            &  ! Virtual potential temperature of parcel in layer k (K)
 ,the_v                ! Virtual potential temperature of environment in
                       ! layer k (K)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DD_INIT'


!-----------------------------------------------------------------------
! Calculate mean temperature, mixing ratio, U, V  and tracer.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1,npnts
  th_mean(i) = (the_k(i)+th_ud_k(i))*0.5
  q_mean(i)  = (qe_k(i)+q_ud_k(i))*0.5
  t_mean(i)  = th_mean(i)*exk(i)
END DO


! For tracers

IF (l_tracer) THEN

  DO ktra=1,ntra
    DO i=1,npnts
      tra_mean(i,ktra) = (trae_k(i,ktra)+tra_ud_k(i,ktra))*0.5
    END DO
  END DO

END IF

!-----------------------------------------------------------------------
! Calculate saturated downdraught potential temperature for layer k
!-----------------------------------------------------------------------

CALL satcal(npnts,k,t_mean,th_mean,pk,exk,q_mean,the_k,qse_k,qdds,thdds)

!-----------------------------------------------------------------------
! Is saturated parcel negatively buoyant compared to environment?
!-----------------------------------------------------------------------

DO i=1,npnts
  IF (.NOT. bdd_on(i) .AND. bddi(i) ) THEN
    thdd_v = thdds(i)*(1.0+c_virtual*qdds(i))
    the_v  = the_k(i)*(1.0+c_virtual*qe_k(i))
    buoy(i) = thdd_v - the_v

    IF (buoy(i)  <   0.5 ) THEN

      !-----------------------------------------------------------------------
      ! Initiate downdraught
      !-----------------------------------------------------------------------

      thdd_k(i) = thdds(i)
      qdd_k(i) = qdds(i)
      bdd_start(i) = .TRUE.

      !-----------------------------------------------------------------------
      ! Calculate cooling and moistening to achieve saturation
      !-----------------------------------------------------------------------

      deltd(i) = thdds(i)-the_k(i)
      delqd(i) = qdds(i)-qe_k(i)
    END IF
  END IF
END DO


IF (l_tracer) THEN

  DO ktra=1,ntra
    DO i=1,npnts
      IF (.NOT. bdd_on(i) .AND. bddi(i) .AND. k >= 4) THEN
        IF (buoy(i) <  0.5) THEN
          tradd_k(i,ktra) = tra_mean(i,ktra)
          deltrad(i,ktra) = tradd_k(i,ktra)-trae_k(i,ktra)
        END IF
      END IF
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE dd_init

END MODULE dd_init_mod
