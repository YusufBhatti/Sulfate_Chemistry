! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine up_ancil in Module up_ancil_mod
!
!   Programing standard : UMDP3
!
!   Purpose: The routine is entered when any of the ancillary
!           fields have to be updated. The list of fields is
!           searched for update requirements. The position of
!           the data required is found from the header information
!           read in from subroutine IN_ANCIL. The data is read in
!           and updates the existing information.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: ancillaries

MODULE up_ancil_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UP_ANCIL_MOD'

CONTAINS

SUBROUTINE up_ancil(i_ao, icode, cmessage)

!Use in relevant JULES routines
USE surf_couple_ancil_update_mod, ONLY: surf_couple_ancil_update

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE atm_fields_mod, ONLY: tstar, land, tstar_anom, frac_land,  &
                          tstar_land, tstar_sea, tstar_sice,   &
                          ice_fraction
USE d1_array_mod, ONLY: d1
USE dump_headers_mod, ONLY: a_realhd
USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE Control_Max_Sizes
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
USE submodel_mod, ONLY: atmos_im
USE nlstcall_mod, ONLY: ancil_reftime, lcal360
USE umPrintMgr, ONLY: umPrint, umMessage

USE ancilcta_namelist_mod, ONLY: nancil_lookupsa

USE nlsizes_namelist_mod, ONLY:                                        &
    land_field, len_dumphist, len_fixhd, len_tot,                      &
    pp_len_inthd, pp_len_realhd, len1_lookup, rows, sm_levels,         &
    theta_field_size, u_field_size, v_field_size

USE model_time_mod, ONLY:                                                    &
    basis_time_days, basis_time_secs, i_day, i_day_number, i_hour, i_minute, &
    i_month, i_second, i_year, stepim

USE replanca_mod, ONLY: replanca

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!   Arguments
!

INTEGER :: i_ao      ! Sub-model indicator = 1 Atmosphere
INTEGER :: icode     ! Return code

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message
CHARACTER(LEN=*),   PARAMETER :: RoutineName = 'UP_ANCIL'

! Local Storage
INTEGER :: ancil_ref_days,ancil_ref_secs
INTEGER :: ancil_offset_steps ! offset of ref. from basis time
LOGICAL :: smc_updated  ! T if soil moisture updated
INTEGER :: a_steps_per_hr
REAL :: dz_soil(sm_levels)  ! soil level thicknesses
LOGICAL :: lland(theta_field_size)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=' '

!   Convert ancillary reference time to days & secs
! DEPENDS ON: time2sec
CALL time2sec(ancil_reftime(1),ancil_reftime(2),                &
              ancil_reftime(3),ancil_reftime(4),                &
              ancil_reftime(5),ancil_reftime(6),                &
              0,0,ancil_ref_days,ancil_ref_secs,lcal360)

!   Compute offset in timesteps of basis time from ancillary ref.time
! DEPENDS ON: tim2step
CALL tim2step(basis_time_days-ancil_ref_days,                   &
              basis_time_secs-ancil_ref_secs,                   &
              steps_per_periodim(i_ao),secs_per_periodim(i_ao), &
              ancil_offset_steps)

IF (i_ao == 1) THEN
  ! Compute A_STEPS_PER_HR for use in REPLANCA
  a_steps_per_hr = 3600*steps_per_periodim(atmos_im)/               &
                         secs_per_periodim(atmos_im)

  lland = TRANSFER(land,lland,theta_field_size)

  CALL replanca(i_year,i_month,i_day,i_hour,i_minute,i_second,    &
                i_day_number,ancil_reftime,ancil_offset_steps,    &
                theta_field_size,rows,u_field_size,               &
                v_field_size,                                     &
                d1,lland,                                          &
                stepim(i_ao),land_field,a_steps_per_hr,           &
                frac_land,                                        &
                tstar_land,tstar_sea,                             &
                tstar_sice,                                       &
                ice_fraction,tstar,                               &
                tstar_anom,sm_levels,                             &
                dz_soil,smc_updated,                              &
                a_realhd(2),a_realhd(3),                          &
                len1_lookup,len_fixhd,                            &
                pp_len_inthd,pp_len_realhd,len_tot,               &
                nancil_lookupsa,                                  &
                icode,cmessage,lcal360)

  IF (icode  /=  0) THEN
    WRITE(umMessage,*) ' UP_ANCIL : Error in REPLANCA.'
    CALL umPrint(umMessage,src='up_ancil')

    CALL Ereport (RoutineName,icode,Cmessage)
  END IF

  !Update JULES
  CALL surf_couple_ancil_update(smc_updated, dz_soil)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE up_ancil

END MODULE up_ancil_mod
