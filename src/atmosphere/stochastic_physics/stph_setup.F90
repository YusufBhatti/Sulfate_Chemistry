! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_setup_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_SETUP_MOD'

CONTAINS


SUBROUTINE stph_read_ensmem(length)
! Reads ensemble member from ENVIRONMENT variable

USE stochastic_physics_run_mod, ONLY: stph_nens, ens_member
USE get_env_var_mod, ONLY: get_env_var
USE umPrintMgr

IMPLICIT NONE

INTEGER :: length         ! length of env var contents

ens_member = '00000000'
! allow_missing=.TRUE. & allow_empty=.TRUE. required for calls from initial_4A
CALL get_env_var('ENS_MEMBER', ens_member, allow_missing=.TRUE.,        &
                 allow_empty=.TRUE., length=length)
ens_member = TRIM(ens_member)
IF (length > 0) THEN
  WRITE(umMessage,'(A,A8)') 'Successfully retrieved ens_member = ',     &
        ens_member
  CALL umPrint(umMessage,src='stph_setup')

  ! Convert STRING to INTEGER value
  READ(ens_member, '(I8)') stph_nens
END IF

END SUBROUTINE stph_read_ensmem


SUBROUTINE stph_setup( )

! Generate and store or read Random Seed for STOCHASTIC PHYSICS
!  Driven by logicals: l_rp2, l_skeb2, l_spt
!            integer : stphseed
!
! --------------------------------------------------------------------
! Sets up a seed that should vary from run to run, i.e. to ensure that
!  each ensemble member gets a different set of random numbers
!  Every RANDOM call after this changes the seed in a predictable way.
!
! Notes on the use of the random seed:
! stphseed: ! 0 => Use ensemble member and date/time of dump
!             1 => Use date/time from computer clock
!             2 => Use seed value stored in dump file
!
! In all instances, the SEED value is written to standard output.
!
! Creates masks used for smoothing and damping SKEB2 fields if set
!  by logicals. These are saved in memory to reduce runtime.
! 

! Main SKEB2 switch
USE stochastic_physics_run_mod,  ONLY:                                  &
    ! Seed variables
    stphseed,                                                           &
    ! SKEB variables
    l_skeb2, nsmooth, offx_stph, offy_stph, mask_smooth,                &
    mask_pdamp, l_skebsmooth_adv, l_skeb2_psisdisp,                     &
    ! RP variables
    l_rp2,                                                              &
    ! Add SPT variables
    l_spt, nsmooth_spt,offx_spt, offy_spt, mask_smooth_spt

! Value of pi
USE conversions_mod, ONLY: pi

! ENDGame compatible array-bounds pointers
USE atm_fields_bounds_mod, ONLY:                                        &
         pdims, stphdims_l

USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,           ONLY: mype, nproc
USE Control_Max_Sizes
USE stph_closeinput_mod,  ONLY: stph_closeinput
USE stph_closeoutput_mod, ONLY: stph_closeoutput
USE stph_openinput_mod,   ONLY: stph_openinput
USE stph_openoutput_mod,  ONLY: stph_openoutput
USE stochastic_physics_run_mod, ONLY: stph_nens, ens_member
USE stph_seed_mod,        ONLY: &
    stph_seed_copy_from_dump, stph_seed_gen_from_date
USE timestep_mod,         ONLY: timestep_number

USE model_time_mod, ONLY: &
    i_day, i_hour, i_minute, i_month, i_year
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! This include contains the current model date

! parameter for finding PE0
INTEGER, PARAMETER ::  stph_zero = 0   ! Zero PE

INTEGER :: length   ! length of string setting ens. member no.

! local Arrays for Stochastic Physics Random Number generation
INTEGER :: icode    ! error code
INTEGER :: errcode  ! Error code for ereport
INTEGER :: dt(8)    ! date/time info for the random seed
INTEGER :: i, j, n, ip1, im1, jp1, jm1

! local arrays for smoothing mask
REAL              :: r_dist   ! Radial distance for mask_pdamp
REAL, ALLOCATABLE :: mask_smooth_tmp(:,:)
REAL, ALLOCATABLE :: mask_smooth_tmp_spt(:,:)

! local temporary arrays
CHARACTER(LEN=errormessagelength)       :: cmessage      ! out error message
CHARACTER(LEN=10000)     :: lineBuffer    ! message for umPrint
CHARACTER(LEN=*), PARAMETER  :: RoutineName='STPH_SETUP'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Keep a record of the settings for this run in the PE output files

IF (printstatus  >=  prstatus_normal) THEN
  CALL umPrint('',src='stph_setup')
  CALL umPrint('*** STPH SETUP ***',src='stph_setup')
  WRITE(umMessage,'(A,L1)') 'L_RP2 =', l_rp2
  CALL umPrint(umMessage,src='stph_setup')
  WRITE(umMessage,'(A,L1)') 'L_SKEB2 = ', l_skeb2
  CALL umPrint(umMessage,src='stph_setup')
  WRITE(umMessage,'(A,L1)') 'L_SPT = ', l_spt
  CALL umPrint(umMessage,src='stph_setup')
  WRITE(umMessage,'(A,I4)') 'STPHSEED = ', stphseed
  CALL umPrint(umMessage,src='stph_setup')
END IF


IF (l_skebsmooth_adv) THEN
  ! Set Halo size for spatial smoothing to number of smoothing iterations
  offx_stph = nsmooth
  offy_stph = nsmooth
  stphdims_l%i_start = pdims%i_start - offx_stph
  stphdims_l%i_end   = pdims%i_end   + offx_stph
  stphdims_l%i_len   = stphdims_l%i_end - stphdims_l%i_start + 1
  stphdims_l%j_start = pdims%j_start - offy_stph
  stphdims_l%j_end   = pdims%j_end   + offy_stph
  stphdims_l%j_len   = stphdims_l%j_end - stphdims_l%j_start + 1
  stphdims_l%k_start = pdims%k_start
  stphdims_l%k_end   = pdims%k_end
  stphdims_l%k_len   = stphdims_l%k_end - stphdims_l%k_start + 1
  stphdims_l%halo_i  = offx_stph
  stphdims_l%halo_j  = offy_stph

  ! Setup pattern damping mask ( 0 in middle -> 1 at [offy_stph+1] )
  IF (.NOT. ALLOCATED(mask_pdamp))                                      &
    ALLOCATE(mask_pdamp(-offx_stph:offx_stph, -offy_stph:offy_stph))
  DO j = -offy_stph, offy_stph
    DO i = -offx_stph, offx_stph
      r_dist = pi * SQRT(REAL(i**2 + j**2))/REAL(offx_stph + 1)
      r_dist = MIN(r_dist, pi)   ! Maximum radius
      mask_pdamp(i,j) = 0.5 - 0.5 * COS(r_dist)
    END DO
  END DO
  IF (printstatus  ==  prstatus_diag) THEN
    WRITE(lineBuffer,'(441ES11.4)') mask_pdamp
    CALL umPrint(lineBuffer,src='stph_setup')
  END IF

  ! Setup spatial smoothing mask
  IF (.NOT. ALLOCATED(mask_smooth))                                     &
    ALLOCATE(mask_smooth(-offx_stph-1:offx_stph+1,                      &
                         -offy_stph-1:offy_stph+1))
  IF (.NOT. ALLOCATED(mask_smooth_tmp))                                 &
    ALLOCATE(mask_smooth_tmp(-offx_stph-1:offx_stph+1,                  &
                           -offy_stph-1:offy_stph+1))
  DO j = -offy_stph-1, offy_stph+1
    DO i = -offx_stph-1, offx_stph+1
      mask_smooth(i,j) = 0.0
      mask_smooth_tmp(i,j) = 0.0
    END DO
  END DO
  mask_smooth(0,0) = 1.0
  DO n = 1, nsmooth
    IF (printstatus  ==  prstatus_diag) THEN
      WRITE(umMessage,'(A,ES22.15)') 'SUM(mask_smooth)= ',              &
            SUM(mask_smooth)
      CALL umPrint(umMessage,src='stph_setup')
      WRITE(lineBuffer,'(529ES11.4)') mask_smooth
      CALL umPrint(lineBuffer,src='stph_setup')
    END IF
    ! Keep smoothed values in temporary array
    DO j = -offy_stph, offy_stph
      jm1 = j - 1
      jp1 = j + 1
      DO i = -offx_stph, offx_stph
        im1 = i - 1
        ip1 = i + 1
        mask_smooth_tmp(i,j) = 0.0625*(mask_smooth(im1,jp1) +           &
                                       mask_smooth(ip1,jp1) +           &
                                       mask_smooth(im1,jm1) +           &
                                       mask_smooth(ip1,jm1) +           &
                                       2*( mask_smooth(i,jp1) +         &
                                           mask_smooth(im1,j) +         &
                                           mask_smooth(ip1,j) +         &
                                           mask_smooth(i,jm1) ) +       &
                                         4*mask_smooth(i,j) )

      END DO
    END DO
    ! Update main array
    DO j = -offy_stph-1, offy_stph+1
      DO i = -offx_stph-1, offx_stph+1
        mask_smooth(i,j) = mask_smooth_tmp(i,j)
      END DO
    END DO
  END DO

  IF (printstatus  ==  prstatus_diag) THEN
    WRITE(umMessage,'(A,ES22.15)') 'SUM(mask_smooth)= ',                &
          SUM(mask_smooth)
    CALL umPrint(umMessage,src='stph_setup')
    WRITE(lineBuffer,'(529ES11.4)') mask_smooth
    CALL umPrint(lineBuffer,src='stph_setup')
  END IF

  DEALLOCATE(mask_smooth_tmp)

END IF ! l_skebsmooth_adv

! +++++++++ Build SPT smoothing field if requested

IF (l_spt) THEN
  ! Set Halo size for spatial smoothing to number of smoothing iterations
  offx_spt = nsmooth_spt
  offy_spt = nsmooth_spt

  ! Setup spatial smoothing mask
  IF (.NOT. ALLOCATED(mask_smooth_spt))                                 &
    ALLOCATE(mask_smooth_spt(-offx_spt-1:offx_spt+1,                    &
                         -offy_spt-1:offy_spt+1))
  IF (.NOT. ALLOCATED(mask_smooth_tmp_spt))                             &
    ALLOCATE(mask_smooth_tmp_spt(-offx_spt-1:offx_spt+1,                &
                            -offy_spt-1:offy_spt+1))
  DO j = -offy_spt-1, offy_spt+1
    DO i = -offx_spt-1, offx_spt+1
      mask_smooth_spt(i,j) = 0.0
      mask_smooth_tmp_spt(i,j) = 0.0
    END DO
  END DO

  mask_smooth_spt(0,0) = 1.0
  DO n = 1, nsmooth_spt

    ! Keep smoothed values in temporary array
    DO j = -offy_spt, offy_spt
      jm1 = j - 1
      jp1 = j + 1
      DO i = -offx_spt, offx_spt
        im1 = i - 1
        ip1 = i + 1
        mask_smooth_tmp_spt(i,j) = 0.0625*(mask_smooth_spt(im1,jp1) +   &
                                        mask_smooth_spt(ip1,jp1) +      &
                                        mask_smooth_spt(im1,jm1) +      &
                                        mask_smooth_spt(ip1,jm1) +      &
                                        2*( mask_smooth_spt(i,jp1) +    &
                                            mask_smooth_spt(im1,j) +    &
                                            mask_smooth_spt(ip1,j) +    &
                                            mask_smooth_spt(i,jm1) ) +  &
                                          4*mask_smooth_spt(i,j) )

      END DO
    END DO
    ! Update main array
    DO j = -offy_spt-1, offy_spt+1
      DO i = -offx_spt-1, offx_spt+1
        mask_smooth_spt(i,j) = mask_smooth_tmp_spt(i,j)
      END DO
    END DO
  END DO

  IF (printstatus  ==  prstatus_diag) THEN
    WRITE(umMessage,'(A,ES22.15)') 'SUM(mask_smooth_spt)= ',            &
                                     SUM(mask_smooth_spt)
    CALL umPrint(umMessage,src='stph_setup')
    WRITE(lineBuffer,'(529ES11.4)') mask_smooth_spt
    CALL umPrint(lineBuffer,src='stph_setup')
  END IF

  DEALLOCATE(mask_smooth_tmp_spt)

END IF ! l_spt

! Random seed is only set on PE=0, as all calls to random are done on
! this PE and then broadcast to other PEs
IF (mype == stph_zero) THEN

  IF (stphseed == 2 .OR. l_nrun_as_crun) THEN
    CALL stph_seed_copy_from_dump

    IF (l_nrun_as_crun) THEN
        cmessage = 'l_nrun_as_crun is TRUE:' //  &
           'Reading random number seed from the dump (equivalent of stphseed=2)'
        icode = -1
        CALL ereport(RoutineName, icode, cmessage)
    END IF
  END IF

  IF (stphseed == 1 .AND. .NOT. l_nrun_as_crun) THEN
    CALL stph_seed_gen_from_date
  END IF

! Read ENS_MEMBER from environment
  CALL stph_read_ensmem(length)

  IF (stphseed == 0 .AND. .NOT. l_nrun_as_crun) THEN
    ! -----------------------------------------------------------
    ! Use model date and ENVIRONMENT variable ENS_MEMBER
    ! -----------------------------------------------------------
    IF (length < 1) THEN
      ! Force abort if ENS_MEMBER cannot be read for option stphseed=0
      CALL umPrint( '**ERROR**: Stochastic Physics Initial Random Seed')
      CALL umPrint( '  Problem retrieving ENS_MEMBER from ENVIRONMENT')
      CALL umPrint( '  Section 35: ENS_MEMBER should be set when using')
      CALL umPrint( '    seed option stphseed = 0                     ')
      WRITE (cmessage,'(A)') 'STPH_SETUP: Cannot retrieve environment'  &
                             //' variable ENS_MEMBER.'
      icode = 10

      CALL ereport(routinename, icode, cmessage)
    END IF
    IF (stph_nens < 0) THEN
      ! Force abort if ENS_MEMBER is negative
      CALL umPrint( '**ERROR**: SKEB2 Initial Random Seed')
      CALL umPrint( '  User provided ENS_MEMBER cannot be negative')
      CALL umPrint( '  Section 35: check namelist settings')
      WRITE (cmessage,'(A)') 'STPH_SETUP: Environment variable'         &
                             //' ENS_MEMBER cannot be negative.'

      icode = 10

      CALL ereport(routinename, icode, cmessage)

    ELSE
      ! Stochastic physics in non-ensemble model requires ENS_MEMBER=0
      !  to be set as an environment variable,
      !  in the case of no seed file being used. This is set later
      !  on all PEs, so need to bcast ens_member (nens) here from PE0

      dt(1) = i_year
      dt(2) = i_month
      dt(3) = i_day
      dt(4) = 0                ! Shift from UTC
      dt(5) = i_hour + 1       ! Set range 1 - 24
      dt(6) = i_minute
      dt(7) = stph_nens        ! Ens mem into second dimension
      dt(8) = stph_nens + 100  ! Ens mem +100 into millisecond dim

      CALL stph_seed_gen_from_date(dt)

    END IF    ! stph_nens >= 0
  END IF
END IF ! mype == stph_zero

! BCAST ens_member to all PEs
CALL gc_ibcast(1,1,stph_zero,nproc,icode,stph_nens)
IF (icode /= 0) THEN
  WRITE (cmessage,'(A)') 'STPH_SETUP: ens_member not bcast correctly'
  errcode = ABS(icode)
  CALL ereport(routinename, errcode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE stph_setup
END MODULE stph_setup_mod
