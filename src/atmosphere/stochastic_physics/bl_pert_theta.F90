! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Add random perturbations to the near surface theta field
!
MODULE bl_pert_theta_mod

IMPLICIT NONE
SAVE

! Description:
!  Randomly perturb the near surface potential temperature field to
!  help initiate convection in high-resolution models. Perturbations
!  are added up to (2/3)*NTML on points diagnosed as CUMULUS only.
!  The magnitude of the perturbations is either specified via the
!  namelist or based on the surface buoyancy flux.  An additional option to
!  use time correlated pertubations is specified through the namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic_physics
!
! Code description:
!  Language: Fortran 95.
!  This code is written to UMDP3 standards.

! random numbers held here from 1st to last EG cycle
REAL, ALLOCATABLE :: RandomNumbers(:,:)
REAL, ALLOCATABLE :: rand_numb(:,:) ! Array of random numbers
REAL, ALLOCATABLE :: pert_flag(:,:) ! Array of flags: 0.0 or 1.0

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BL_PERT_THETA_MOD'

CONTAINS

! subroutine interface
SUBROUTINE bl_pert_theta(cycleno, theta, theta_star_surf,               &
                         qv, qv_star_surf, ntml, cumulus,               &
                         z_theta, stashwork3)

USE atm_fields_bounds_mod, ONLY: tdims, tdims_s
USE dynamics_input_mod, ONLY: numcycles
USE model_domain_mod, ONLY: model_type, mt_global
USE nlsizes_namelist_mod,  ONLY: global_row_length, global_rows
USE stochastic_physics_run_mod, ONLY: mag_pert_theta,                   &
     minlev_pert_theta, maxlev_pert_theta,                              &
     npts_pert_theta, l_pert_all_points, l_pert_shape, i_pert_theta,    &
     pert_theta_mag, pert_theta_star, pert_theta_and_moist,             &
     decorr_ts_pert_theta, i_pert_theta_type, pert_theta_random_seq,    &
     pert_theta_correl_seq, stphseed
USE update_pert_mod, ONLY: update_pert
USE UM_ParParams, ONLY: PNorth, PSouth
USE um_parvars, ONLY: at_extremity, datastart
USE um_parcore, ONLY: mype, nproc
USE umPrintMgr, ONLY: umPrint, umMessage
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE stash_array_mod, ONLY: si, sf
USE submodel_mod, ONLY: atmos_im
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE timestep_mod, ONLY: timestep, timestep_number
USE atm_fields_real_mod, ONLY: bl_pert_rand_fld, bl_pert_flag

IMPLICIT NONE

REAL, INTENT(INOUT) ::                                                  &
 stashwork3(*)
                   ! STASH workspace for section 3 (b layer)

! Arguments with INTENT IN. ie: Input variables.
INTEGER, INTENT(IN) :: cycleno
INTEGER, INTENT(IN) :: ntml(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end)
LOGICAL, INTENT(IN) :: cumulus(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: theta_star_surf(tdims%i_start:tdims%i_end,          &
                                    tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: qv_star_surf(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: z_theta(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            1:tdims%k_end)

! Arguments with INTENT INOUT
REAL, INTENT (INOUT) ::                                                 &
    theta(tdims%i_start:tdims%i_end,                                    &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end),                                   &
    qv(tdims%i_start:tdims%i_end,                                       &
          tdims%j_start:tdims%j_end,                                    &
          tdims%k_start:tdims%k_end)

! Local Variables
LOGICAL :: l_add_pert(tdims%i_start:tdims%i_end,                        &  
                      tdims%j_start:tdims%j_end)     
                               ! Array of logical switches set to TRUE
                               ! when perturbations are applied,
                               ! FALSE otherwise
LOGICAL, SAVE :: l_firstcall = .TRUE.                       

REAL, SAVE:: auto_corr_coeff   ! Auto-correlation coefficient for time
                               ! correlation of theta perturbations
REAL :: rand_no                ! Random number for use in update of perturbation
REAL :: rand_min, rand_max     ! Minimum/maximum random number
REAL :: shock_amp              ! Shock amplitude for perturbation update
REAL :: theta_pert(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                               ! perturbation to theta (K)
REAL :: qv_pert(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                               ! perturbation to qv (kg/kg)
REAL :: mag_pert_qv ! Maximum perturbation to qv (10% of qv)
REAL :: zh                     ! Height of PBL top
REAL :: zfac                   ! Height factor for perturbation
REAL :: tfac                   ! Timestep dependence

INTEGER :: i,j,k    ! loop variables
INTEGER :: seed_len ! Length of random number seed vector
INTEGER :: j_start,j_end ! variables to avoid poles
INTEGER :: gi,gj    ! global loop variables
INTEGER :: istat=0  ! communication status

INTEGER, ALLOCATABLE :: seed(:) ! Vector random number seed

INTEGER :: item, icode, im_index
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BL_PERT_THETA'

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( cycleno == 1 ) THEN
  ALLOCATE (rand_numb(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end))
  ALLOCATE (pert_flag(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end))
END IF

IF ( i_pert_theta_type == pert_theta_correl_seq ) THEN
  ! For time correlated perturbations, initialise bl_pert_flag and 
  ! auto_corr_coeff the first time the module is called.
  ! For the case of CRUNS or NRUNS-with-stphseed=2, then rand_numb
  ! and pert_flag are read in from the previous dump

  IF (l_firstcall) THEN
    auto_corr_coeff = EXP(-timestep/decorr_ts_pert_theta)

    IF (l_nrun_as_crun) THEN
        rand_numb = bl_pert_rand_fld
        pert_flag = bl_pert_flag
        cmessage = 'l_nrun_as_crun is TRUE:' //  &
           'Reading random numbers from the dump (equivalent of stphseed=2)'
        icode = -1
        CALL ereport(RoutineName, icode, cmessage)
    ELSE
      IF ( timestep_number == 1 .AND. stphseed < 2 ) THEN
        ! Initialise pert_flag to be 0.0 everywhere
        pert_flag(:,:) = 0.0
        rand_numb(:,:) = 0.0
      ELSE
        ! CRUN or NRUN-with-stphseed=2
        ! Set rand_numb from values read in from the dump
        rand_numb = bl_pert_rand_fld
        pert_flag = bl_pert_flag
      END IF
    END IF

    l_firstcall = .FALSE.
  ELSE
    ! Assign rand_numb and pert_flag values stored
    ! in prognostics bl_pert_rand_fld and bl_pert_flag
    rand_numb = bl_pert_rand_fld
    pert_flag = bl_pert_flag
  END IF

END IF ! pert_theta_correl_seq

IF (model_type==mt_global .AND. at_extremity(PSouth)) THEN
  j_start=2
ELSE
  j_start=1
END IF
IF (model_type==mt_global .AND. at_extremity(PNorth)) THEN
  j_end=tdims%j_end-1
ELSE
  j_end=tdims%j_end
END IF

! Set appropriate values for local parameters used in time-correlating
! random number sequence
rand_min = -1.0
rand_max =  1.0
shock_amp = 1.0

! Very crude attempt to introduce a timestep dependence, because we are 
! adding an increment every timestep (60s is the UKV timestep so this
! shouldn't change the answers there)
IF ( i_pert_theta_type == pert_theta_correl_seq ) THEN
  tfac=timestep/60.0
ELSE
  tfac=1.0
END IF

! Initialise perturbations and logicals
l_add_pert(:,:)=.FALSE.
theta_pert(:,:)=0.0
qv_pert(:,:)=0.0

! initialise the random numbers at first EG cycle only
IF (cycleno == 1) THEN

  ! Initialise random number seed.
  CALL RANDOM_SEED (SIZE = seed_len)
  ALLOCATE (seed(seed_len))

  ! Get the current random seed on PE0 and scatter it to the rest
  IF (mype == 0) THEN
    CALL RANDOM_SEED (GET = seed)
  END IF
  CALL gc_rbcast(1234,seed_len,0,nproc,istat,seed)
  IF (mype /= 0) THEN
    CALL RANDOM_SEED (PUT = seed)
  END IF

  WRITE(umMessage,*)'Perturbing theta field with seed: ',seed
  CALL umPrint(umMessage,src='bl_pert_theta')

  DEALLOCATE (seed)

  ! Get random numbers in the range 0-1:
  ALLOCATE(RandomNumbers(global_row_length/npts_pert_theta+1,           &
                         global_rows/npts_pert_theta+1))

  CALL RANDOM_NUMBER (RandomNumbers)

END IF !1st cycle

!---------------------
! Perturb theta field
!---------------------

IF (i_pert_theta == pert_theta_star) THEN
  ! perturbations proportional to theta_star
  DO j=j_start, j_end
    gj=j+datastart(2)-1
    DO i=tdims%i_start, tdims%i_end
      gi=i+datastart(1)-1
      IF ( theta_star_surf(i,j) > 0.0 .AND.                             &
          (cumulus(i,j) .OR. l_pert_all_points) ) THEN
        ! Calulate perturbation
        l_add_pert(i,j) = .TRUE.

        ! First get random number and set pert_amp
        rand_no =                                                       & 
           RandomNumbers((gi+npts_pert_theta-1)/npts_pert_theta,        &
                         (gj+npts_pert_theta-1)/npts_pert_theta)
        IF ( i_pert_theta_type == pert_theta_random_seq                 &
              .OR. ( pert_flag(i,j) == 0.0 ) ) THEN
          rand_numb(i,j) = 2.0*(rand_no - 0.5)
          pert_flag(i,j) = 1.0
        ELSE IF (cycleno == 1) THEN
          ! Time correlate random numbers
          CALL update_pert( rand_numb(i,j), auto_corr_coeff,            &
                            rand_no, shock_amp, rand_min, rand_max )
        END IF

        theta_pert(i,j) = rand_numb(i,j) *                              &
               MIN( tfac*theta_star_surf(i,j),mag_pert_theta )

      END IF   ! theta_star gt 0 & (cumulus or all_points) 
    END DO
  END DO


  DO j=j_start, j_end
    DO i=tdims%i_start, tdims%i_end
      IF ( l_add_pert(i,j) ) THEN
        DO k=minlev_pert_theta, MIN( maxlev_pert_theta , 2*ntml(i,j)/3 )

          ! Perturb theta
          theta(i,j,k) = theta(i,j,k) + theta_pert(i,j)

        END DO 
      END IF   ! l_add_pert
    END DO 
  END DO

ELSE IF (i_pert_theta == pert_theta_and_moist) THEN
  ! perturbations to theta and moisture
  ! proportional to theta_star and qv_star

  DO j=j_start, j_end
    gj=j+datastart(2)-1
    DO i=tdims%i_start, tdims%i_end
      gi=i+datastart(1)-1
      IF ( theta_star_surf(i,j) > 0.0 .AND.                             &
           (cumulus(i,j) .OR. l_pert_all_points) ) THEN
        ! Calulate perturbation
        l_add_pert(i,j) = .TRUE.

        ! First get random number
        rand_no =                                                       & 
            RandomNumbers((gi+npts_pert_theta-1)/npts_pert_theta,       &
                          (gj+npts_pert_theta-1)/npts_pert_theta)
        IF ( i_pert_theta_type == pert_theta_random_seq                 &
             .OR. ( pert_flag(i,j) == 0.0 ) ) THEN
          rand_numb(i,j) = 2.0*(rand_no - 0.5)
          pert_flag(i,j) = 1.0
        ELSE IF (cycleno == 1) THEN
          ! Time correlate random numbers
          CALL update_pert( rand_numb(i,j), auto_corr_coeff,            &
                  rand_no, shock_amp, rand_min, rand_max )
        END IF

        ! Calculate theta and moisture perturbationa
        theta_pert(i,j) = rand_numb(i,j) *                              &
                   MIN( tfac*theta_star_surf(i,j),mag_pert_theta )
        mag_pert_qv = 0.1*qv(i,j,1)
        qv_pert(i,j) = rand_numb(i,j) *                                 &
                   MIN( tfac*qv_star_surf(i,j), mag_pert_qv)

      END IF   ! theta_star gt 0 & (cumulus or all_points)
    END DO
  END DO



  DO j=j_start, j_end
    DO i=tdims%i_start, tdims%i_end
      IF ( l_add_pert(i,j) ) THEN
        IF (l_pert_shape) THEN

          DO k=minlev_pert_theta, MIN( maxlev_pert_theta , ntml(i,j) )
            zh = z_theta(i,j,ntml(i,j)+1)
            zfac = MIN( 2.0*z_theta(i,j,k)/zh,                          &
                        2.0*(1.0-z_theta(i,j,k)/zh) )
            ! Perturb theta
            theta(i,j,k) = theta(i,j,k) + theta_pert(i,j)*zfac
            ! Perturb qv
            qv(i,j,k) = qv(i,j,k) + qv_pert(i,j)*zfac
          END DO

        ELSE

          DO k=minlev_pert_theta, MIN( maxlev_pert_theta , 2*ntml(i,j)/3 )
            ! Perturb theta
            theta(i,j,k) = theta(i,j,k) + theta_pert(i,j)
            ! Perturb qv
            qv(i,j,k) = qv(i,j,k) + qv_pert(i,j)
          END DO

        END IF   ! l_pert_shape
      END IF   ! l_add_pert
    END DO
  END DO

ELSE IF (i_pert_theta == pert_theta_mag) THEN
  ! perturbations used fixed input magnitude, mag_pert_theta

  DO j=j_start, j_end
    gj=j+datastart(2)-1
    DO i=tdims%i_start, tdims%i_end
      gi=i+datastart(1)-1
      IF ( cumulus(i,j) .OR. l_pert_all_points ) THEN
        ! Calulate perturbation
        l_add_pert(i,j) = .TRUE.

        ! First get random number
        rand_no =                                                       & 
           RandomNumbers((gi+npts_pert_theta-1)/npts_pert_theta,        &
                         (gj+npts_pert_theta-1)/npts_pert_theta)   
        IF ( i_pert_theta_type == pert_theta_random_seq                 &
                   .OR. ( pert_flag(i,j) == 0.0 ) ) THEN
          rand_numb(i,j) = 2.0*(rand_no - 0.5)
          pert_flag(i,j) = 1.0
        ELSE IF (cycleno == 1) THEN
          ! Time correlate perturbations
          CALL update_pert( rand_numb(i,j), auto_corr_coeff,            &
                       rand_no, shock_amp, rand_min, rand_max )
        END IF
        
        theta_pert (i,j) = tfac * rand_numb(i,j) * mag_pert_theta

      END IF   ! cumulus or all_points
    END DO
  END DO


  DO j=j_start, j_end
    DO i=tdims%i_start, tdims%i_end
      IF ( l_add_pert(i,j) ) THEN
        DO k=minlev_pert_theta, MIN( maxlev_pert_theta , 2*ntml(i,j)/3 )

          ! Perturb theta
          theta(i,j,k) = theta(i,j,k) + theta_pert(i,j)

        END DO
      END IF   ! l_add_pert
    END DO
  END DO

END IF

IF ( i_pert_theta_type == pert_theta_correl_seq ) THEN
  ! Pass rand_numb to bl_pert_rand_fld to be written out to
  ! the dump for use in future cruns or nruns-with-stphseed=2

  bl_pert_rand_fld = rand_numb
  bl_pert_flag     = pert_flag

END IF

! random perturbations of theta
item=520
icode=0
im_index = 1
IF (icode <= 0 .AND. sf(item,3) .AND. cycleno==numcycles) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork3(si(item,3,im_index)),theta_pert,             &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                   &
       atmos_im,3,item,                                                 &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag(item 520)"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF
! random perturbations of qv
item=521
icode=0
im_index = 1
IF (i_pert_theta == pert_theta_and_moist .AND. icode <= 0               &
    .AND. sf(item,3) .AND. cycleno==numcycles) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(stashwork3(si(item,3,im_index)),qv_pert,                &
       tdims%i_end,tdims%j_end,0,0,0,0, at_extremity,                   &
       atmos_im,3,item,                                                 &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage=": error in copydiag(item 521)"
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF

! deallocate random numbers on last cycle
IF (cycleno == numcycles) THEN
  DEALLOCATE(RandomNumbers)
  DEALLOCATE(rand_numb)
  DEALLOCATE(pert_flag)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bl_pert_theta

END MODULE bl_pert_theta_mod
