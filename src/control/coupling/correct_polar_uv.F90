! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE correct_polar_uv(               &
   delta_lambda, global_row_length,        &
   ocn_uv_in, ocn_uv_out,                  &
   oasis_jmt_in, oasis_jmt_out)

USE UM_ParVars
USE UM_ParCore,   ONLY: mype
USE UM_ParParams, ONLY: PNorth
USE Control_Max_Sizes
USE oasis_atm_data_mod, ONLY: oasis_imt

USE trignometric_mod, ONLY: sin_theta_longitude,  &
                            cos_theta_longitude,  &
                            sin_u_longitude,      &
                            cos_u_longitude

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!
! Description:
! Calculate consistent values for polar row velocities based on velocity pair 
! values half a grid point below the polar row. We only expect to be dealing 
! with ENDGAME grids so in effect we calculate V at pole based on the 
! northeren-most row of U data (which is offset half a grid point to the 
! south in ENDGAME or V at Poles.)
! This process is necessary to ensure meaningful polar values in coupled
! fields which otherwise tend to contain junk as a result of SCRIP/OASIS
! regridding failing to deal effectively with values on the polar
! singularity.
! This routine only needs to deal with the North Pole since the South pole
! plays no part of ocean-atmos coupling. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!=======================================================================


!     Subroutine arguments
REAL, INTENT(IN) :: delta_lambda ! longitude per grid box

INTEGER, INTENT(IN) :: global_row_length

INTEGER, INTENT(IN) :: oasis_jmt_in   ! Size of incoming velocity array
INTEGER, INTENT(IN) :: oasis_jmt_out  ! Aize of outgoing velocity array,
                          ! typically oasis_jmt_in+1 because it covers the
                          ! N pole.

REAL, INTENT(IN)    :: ocn_uv_in(oasis_imt,oasis_jmt_in)
REAL, INTENT(INOUT) :: ocn_uv_out(oasis_imt,oasis_jmt_out)


!     Local variables
INTEGER :: i, gi, levs   ! Loop counters/indices
INTEGER :: halo_none     ! Size of mpp halo (always zero)

REAL :: mag_vector_np(1),                                &
        dir_vector_np(1),                                &
        mag_vector_sp(1),                                &
        dir_vector_sp(1)

REAL :: posneg           ! +/-1 depending on whether we are calculating
                         ! U based on V values or V based on U values.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CORRECT_POLAR_UV'
REAL(KIND=jprb)               :: aicenmin

! Initialise potential incoming fields
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! First compute total U component (EG).
! Note the fields we are dealing with here are those returned from
! the coupler - i.e. they do not contain halos, so we must indicate
! this accordingly in the following subroutine calls.
! We also need to ensure that the dimensions supplied for the
! input array are similarly devoid of halos.

! We only deal with a single vertical level in surface coupling
! fields, hence the following.
levs=1

halo_none=0

! We have V points at the N pole
! i.e. we're using EG and consequently need to calculate V values based
! on U (half a row below the pole.) Note, therefore, the use of
! longitudes at U points in the following.

posneg = -1.0

! DEPENDS ON: polar_vector_wind_n
CALL Polar_Vector_Wind_n(                                            &
                       ocn_uv_in,                                    &
                       sin_u_longitude,                              &
                       cos_u_longitude, oasis_imt,                   &
                       oasis_jmt_in, levs, mag_vector_np,            &
                       dir_vector_np, mag_vector_sp,                 &
                       dir_vector_sp,                                &
                       halo_none, halo_none, global_row_length ,     &
                       gc_proc_row_group, at_extremity)


IF (at_extremity(PNorth) ) THEN

   ! We now have a single value of magnitude and direction
   ! for cross polar U from which we can work out the full set of
   ! V vectors along the polar row.
   ! Note we only do this for the N pole.
   ! Direction at each point is given by:
   !       longitude(i)-direction of maximum magnitude
   !             (offset horizontally by half a d-lambda due to
   !              U/V point staggering)
   ! and using a sine (not a cosine), because our reference direction
   ! by definition, must be a point where V = 0 (ignoring the half
   ! grid point offset for the moment) because its's where v=maximum.
   ! i.e. V and U maxima are out of phase by 90 degrees.

  DO i = 1, oasis_imt
     ! Other areas in the UM which use this approach
     ! employ something called "l_datastart" to figure out the
     ! local data start position in terms of global indices
     ! which seems to add an unnecessary level of obfuscation.
     ! Here we use g_datastart directly!

    gi = g_datastart(1,mype) + i - 1

    ! As our V points are positioned half a grid point
    ! to the West of our U points we subtract half a d-lambda
    ! when calculating the sine of the longitude at the
    ! equivalent (i,j) V point.
    ocn_uv_out(i,oasis_jmt_out) = posneg * (mag_vector_np(1) *         &
             SIN(((gi-.5)*delta_lambda) - dir_vector_np(1)))
  END DO
END IF  ! at_extremity(PNorth)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE correct_polar_uv
