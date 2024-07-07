! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE mix_inc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MIX_INC_MOD'
CONTAINS

SUBROUTINE mix_inc                                                      &
  ( np_field, npnts, ncpts, nlev, nbl, nwml,                            &
    dthbydt, dqbydt, dubydt, dvbydt, l_tracer, ntra,                    &
    dtrabydt, p_layer_boundaries, p_layer_centres, index4 )

USE cv_run_mod, ONLY: l_mom
USE yomhook,    ONLY: lhook, dr_hook
USE parkind1,   ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------
! Description:
!   To well-mix convective increments in the boundary layer.
!   Suitable for single column model use.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-----------------------------------------------------------------------
! Subroutine Arguments


INTEGER, INTENT(IN) :: np_field ! Length of data
                                ! (also used to specify starting point of
                                !  data passed in)
INTEGER, INTENT(IN) :: npnts    ! Full vector length
INTEGER, INTENT(IN) :: ncpts    ! Vector length
INTEGER, INTENT(IN) :: nlev     ! Number of model levels
INTEGER, INTENT(IN) :: nbl      ! Number of model layers
                                ! potentially in the boundary layer

INTEGER, INTENT(IN) :: nwml(npnts)  ! Number of model layers for which
                                    ! increments are to be well-mixed.

REAL, INTENT(INOUT) :: dthbydt(np_field,nlev ) ! IN  Increment to potential
                                               !     temperature due to
                                               !     convection
                                               ! OUT Modified increment
REAL, INTENT(INOUT) :: dqbydt(np_field,nlev  ) ! IN  Increment to mixing ratio
                                               !     due to convection
                                               ! OUT Modified increment
REAL, INTENT(INOUT) :: dubydt(np_field,nlev  ) ! IN  Increment to x-component
                                               !     of wind due to convection
                                               ! OUT Modified increment
REAL, INTENT(INOUT) :: dvbydt(np_field,nlev  ) ! IN  Increment to y-component
                                               !     of wind due to convection
                                               ! OUT Modified increment

LOGICAL, INTENT(IN) :: l_tracer     ! Logical switch for inclusion
                                    ! of convective mixing of tracers
INTEGER, INTENT(IN) :: ntra         ! Number of tracer fields

REAL, INTENT(INOUT) :: dtrabydt(npnts,nlev,ntra)   ! IN  Increments to tracers
                                                   !     due to convection
                                                   ! OUT Modified increment

REAL, INTENT(IN) ::  p_layer_boundaries (np_field,0:nlev)
REAL, INTENT(IN) ::  p_layer_centres    (np_field,0:nlev)

INTEGER, INTENT(IN) :: index4(npnts)

!----------------------------------------------------------------------
! Local Variables
!----------------------------------------------------------------------
INTEGER :: i,k,itra       ! Loop counters

REAL :: thsum  (ncpts)    ! Summation of increments to potential temperature
                          ! due to convection in the vertical, weighted
                          ! according to mass.

REAL :: qsum   (ncpts)    ! Summation of increments to mixing ratio
                          ! due to convection in the vertical, weighted
                          ! according to mass.

REAL :: usum   (ncpts)    ! Summation of increments to x_component of wind
                          ! due to convection in the vertical,
                          ! weighted according to mass.

REAL :: vsum   (ncpts)    ! Summation of increments to y_component of wind
                          ! due to convection in the vertical,
                          ! weighted according to mass.

REAL :: trasum (ncpts)    ! Summation of increments to tracers
                          ! due to convection in the vertical,
                          ! weighted according to mass.

REAL :: delpsum(ncpts)      ! Sum of model layer thicknesses with height. (Pa)
REAL :: delpsum_uv(ncpts)   ! Sum of model layer thicknesses with height. (Pa)
REAL :: delpk(ncpts,nbl)    ! Difference in pressure across a layer (Pa)
REAL :: delpk_uv(ncpts,nbl) ! Difference in PRESSURE ACROSS a layer (Pa)

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MIX_INC'


!
!----------------------------------------------------------------------
!  Sum up mass weighted rates of change of variables
!  and layer thicknesses for model layers 1 to nwml.
!----------------------------------------------------------------------
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i=1, ncpts
  delpk(i,1) = p_layer_boundaries(index4(i),0) -                        &
               p_layer_boundaries(index4(i),1)

  thsum(i)   = dthbydt(index4(i),1) * delpk(i,1)
  qsum(i)    = dqbydt(index4(i),1)  * delpk(i,1)
  delpsum(i) = delpk(i,1)
END DO

DO k=2, nbl
  DO i=1, ncpts
    IF ( k  <=  nwml(index4(i))) THEN
      delpk(i,k) = p_layer_boundaries(index4(i),k-1) -                  &
                   p_layer_boundaries(index4(i),k)
      thsum(i)   = thsum(i) + dthbydt(index4(i),k) * delpk(i,k)
      qsum(i)    = qsum(i)  + dqbydt(index4(i),k)  * delpk(i,k)
      delpsum(i) = delpsum(i) + delpk(i,k)
    END IF
  END DO
END DO

!
!----------------------------------------------------------------------
! Reset potential temperature and humidity increments in layers 1 to
! nwml to mean values over these layers.
!----------------------------------------------------------------------
!

DO k=1, nbl
  DO i=1, ncpts
    IF ( k <= nwml(index4(i)) ) THEN
      dthbydt(index4(i),k) = thsum(i) / delpsum(i)
      dqbydt(index4(i),k)  = qsum(i)  / delpsum(i)
    END IF
  END DO
END DO

IF (l_mom) THEN

  !
  !----------------------------------------------------------------------
  !  Sum up mass weighted rates of change of variables
  !  for model layers 1 to nwml.
  !----------------------------------------------------------------------
  !

  DO i=1, ncpts
    delpk_uv(i,1) = p_layer_centres(index4(i),0) -                      &
                    p_layer_centres(index4(i),1)
    delpsum_uv(i) = delpk_uv(i,1)
    usum(i) = dubydt(index4(i),1) * delpk_uv(i,1)
    vsum(i) = dvbydt(index4(i),1) * delpk_uv(i,1)
  END DO

  DO k=2, nbl
    DO i=1, ncpts
      IF ( k <= nwml(index4(i)) ) THEN
        delpk_uv(i,k) = p_layer_centres(index4(i),k-1) -                &
                        p_layer_centres(index4(i),k)
        delpsum_uv(i) = delpsum_uv(i) + delpk_uv(i,k)
        usum(i) = usum(i) + dubydt(index4(i),k) * delpk_uv(i,k)
        vsum(i) = vsum(i) + dvbydt(index4(i),k) * delpk_uv(i,k)
      END IF
    END DO
  END DO

  !
  !----------------------------------------------------------------------
  ! Reset wind component increments in layers 1 to NWML to mean values
  ! over these layers.
  !----------------------------------------------------------------------
  !

  DO k=1, nbl
    DO i=1, ncpts
      IF ( k <= nwml(index4(i)) ) THEN
        dubydt(index4(i),k) = usum(i) / delpsum_uv(i)
        dvbydt(index4(i),k) = vsum(i) / delpsum_uv(i)
      END IF
    END DO
  END DO

END IF  ! l_mom

IF (l_tracer) THEN

  !
  !----------------------------------------------------------------------
  !  Sum up mass weighted rates of change of each tracer in turn
  !  for model layers 1 to nwml.
  !----------------------------------------------------------------------
  !

  DO itra=1, ntra
    DO i=1, ncpts
      trasum(i) = dtrabydt(index4(i),1,itra) * delpk(i,1)
    END DO

    DO k=2, nbl
      DO i=1, ncpts
        IF ( k <= nwml(index4(i)) ) THEN
          trasum(i) = trasum(i) + dtrabydt(index4(i),k,itra) * delpk(i,k)
        END IF
      END DO
    END DO

    !
    !----------------------------------------------------------------------
    ! Reset tracer increments in layers 1 to nwml to mean values
    ! over these layers.
    !----------------------------------------------------------------------
    !

    DO k=1, nbl
      DO i=1, ncpts
        IF ( k <= nwml(index4(i)) ) THEN
          dtrabydt(index4(i),k,itra) = trasum(i) / delpsum(i)
        END IF
      END DO
    END DO

  END DO  ! Loop over tracers
END IF  ! l_tracer

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mix_inc
END MODULE mix_inc_mod
