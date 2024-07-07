! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  Module to make the all_tracers array from tracer_ukca and q_um
!  and copy back after doing chemistry and aerosols
!
!  Method:
!  We make a lookup array on the first call and use this to do the copying
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran
!
! ----------------------------------------------------------------------
MODULE ukca_all_tracers_copy_mod

! length of character strings in nmspec
USE ukca_nmspec_mod,    ONLY: nmspec_len
USE ukca_option_mod,    ONLY: l_ukca_h2o_feedback

IMPLICIT NONE

PRIVATE

PUBLIC :: ukca_all_tracers_copy_in, ukca_all_tracers_copy_out, &
          ukca_make_tr_map, arrayloc

! Increment of q before/after chemistry
REAL, ALLOCATABLE, PUBLIC :: ukca_q_increment(:,:,:)

! Array to be populated from tracer_ukca_um and q_um
! all_tracers has the chemistry items first,
! all in the same order as advt, and then the MODE
! items in STASH code order.
! Other tracers e.g. age of air come last
REAL, ALLOCATABLE, PUBLIC :: all_tracers(:,:,:,:)

INTEGER, ALLOCATABLE, SAVE, PUBLIC :: tr_lookup(:)

CHARACTER(LEN=nmspec_len), ALLOCATABLE, SAVE, PUBLIC :: tr_names(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: modulename='UKCA_ALL_TRACERS_COPY_MOD'

CONTAINS 

SUBROUTINE ukca_all_tracers_copy_in(tracer_ukca_um,  q_um)

USE atm_fields_bounds_mod, ONLY: tdims_s, tdims_l
USE ukca_d1_defs, ONLY: n_use_tracers
USE nlsizes_namelist_mod, ONLY: tr_ukca
USE um_parcore, ONLY: mype
USE umPrintMgr
USE yomhook,                   ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

! Populate all_tracers array from tracer_ukca_um

! Input data 
! UKCA Tracers from the UM code
REAL, INTENT(IN) :: tracer_ukca_um                              &
              (tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,1:tr_ukca)
! Specific humidity from the UM code
REAL, INTENT(IN) :: q_um                                        &
               (tdims_l%i_start:tdims_l%i_end,                  &
                tdims_l%j_start:tdims_l%j_end,                  &
                tdims_l%k_start:tdims_l%k_end)

LOGICAL, SAVE :: l_first = .TRUE.
CHARACTER(LEN=*), PARAMETER :: routinename='UKCA_ALL_TRACERS_COPY_IN'

INTEGER :: i                                    ! loop counter

! Variables for use with Dr Hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! on first call to this routine set up the mapping from tracer_ukca_um 
! to all_tracers including an array of names
IF (l_first) THEN 
   CALL ukca_make_tr_map()
   l_first = .FALSE.
END IF 

! Allocate arrays here. Need to deallocate later
IF (.NOT. ALLOCATED(all_tracers))                                             &
     ALLOCATE(all_tracers(tdims_s%i_start:tdims_s%i_end,                      &
                          tdims_s%j_start:tdims_s%j_end,                      &
                          tdims_s%k_start:tdims_s%k_end,1:n_use_tracers))


! Loop over the all_tracers array and copy in data from q_um and tracer_ukca_um
DO i=1,n_use_tracers
   IF (tr_lookup(i) == -1) THEN
      ! Note inconsistent halo sizes here - just copy small halo fields.
      all_tracers(tdims_s%i_start:tdims_s%i_end,                &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end, i) =              &
               q_um(tdims_s%i_start:tdims_s%i_end,              &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end) 

     IF (l_ukca_h2o_feedback) THEN
       ! Store q_um in increment field
       IF (.NOT. ALLOCATED(ukca_q_increment))                          &
           ALLOCATE(ukca_q_increment(tdims_l%i_start:tdims_l%i_end,    &
                                     tdims_l%j_start:tdims_l%j_end,    &
                                     tdims_l%k_start:tdims_l%k_end))
       ukca_q_increment(:,:,:) = q_um(:,:,:)
     END IF
   ELSE
      all_tracers(:,:,:,i) =  tracer_ukca_um(:,:,:,tr_lookup(i))
   END IF 
END DO 

! print the data in all_tracers for debugging purposes
IF (PrintStatus >= PrStatus_Diag) THEN
  DO i=1,n_use_tracers
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'troptracer: ',mype,i,               &
      MAXVAL(all_tracers(:,:,:,i)),MINVAL(all_tracers(:,:,:,i))
    CALL umPrint(umMessage,src=routinename)
  END DO
END IF

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

RETURN

END SUBROUTINE ukca_all_tracers_copy_in

SUBROUTINE ukca_all_tracers_copy_out(tracer_ukca_um,  q_um)

USE atm_fields_bounds_mod
USE cv_run_mod,   ONLY: qmin_conv
USE ukca_d1_defs, ONLY: n_use_tracers
USE nlsizes_namelist_mod, ONLY: tr_ukca
USE ukca_option_mod, ONLY: l_ukca_h2o_feedback
USE level_heights_mod,    ONLY: r_theta_levels
USE planet_constants_mod, ONLY: planet_radius
USE yomhook,                   ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

! Populate tracer_ukca_um array from all_tracers


! UKCA tracer from UM to be updated with data from
! all_tracers
REAL, INTENT(INOUT) :: tracer_ukca_um                           &
              (tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,1:tr_ukca)
! Specific humidity array to be updated with data from
! all_tracers if required
REAL, INTENT(INOUT) :: q_um                                     &
               (tdims_l%i_start:tdims_l%i_end,                  &
                tdims_l%j_start:tdims_l%j_end,                  &
                tdims_l%k_start:tdims_l%k_end)
INTEGER :: i
INTEGER :: j, l, k            ! loop counters
INTEGER :: i1, j1, k1

REAL, PARAMETER :: qmin_buffer = 0.01
! arbitrary buffer, to ensure that returned q is slightly above qmin_conv
REAL, PARAMETER :: qmin_zmin = 7.5E4
! minimum height for correcting very low values of q

CHARACTER(LEN=*), PARAMETER :: routinename='UKCA_ALL_TRACERS_COPY_OUT'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! Loop over the all_tracers array and copy back data to q and tracer_ukca
DO i=1, n_use_tracers
   IF (tr_lookup(i) == -1) THEN
     IF (l_ukca_h2o_feedback) THEN
        ! Note inconsistent halo sizes here - just copy small halo fields.
        q_um(tdims_s%i_start:tdims_s%i_end,                       &
             tdims_s%j_start:tdims_s%j_end,                       &
             tdims_s%k_start:tdims_s%k_end) =                     &
            all_tracers(tdims_s%i_start:tdims_s%i_end,            &
                        tdims_s%j_start:tdims_s%j_end,            &
                        tdims_s%k_start:tdims_s%k_end, i) 

         ! Ensure level 0 is set to values of level 1
           q_um(tdims_l%i_start:tdims_l%i_end,                     &
                tdims_l%j_start:tdims_l%j_end,                     &
                tdims_l%k_start) =                                 &
           q_um(tdims_l%i_start:tdims_l%i_end,                     &
                tdims_l%j_start:tdims_l%j_end, 1)

         ! Reset very low values of q as they provoke copious output
         ! from the convection scheme when q < qmin_conv. Only used
         ! if height exceeds qmin_zmin
         DO k=tdims_l%k_start,tdims_l%k_end
           DO j= tdims_l%j_start,tdims_l%j_end
             DO l= tdims_l%i_start,tdims_l%i_end
               IF ((r_theta_levels(l,j,k) - planet_radius > qmin_zmin) .AND. &
                   q_um(l,j,k) < qmin_conv ) THEN
                 q_um(l,j,k) = qmin_conv*(1.0 + qmin_buffer)
               END IF
             END DO
           END DO
         END DO

         ! Calculate incremental change in q
          ukca_q_increment(:,:,:) = q_um(:,:,:) - ukca_q_increment(:,:,:)
     END IF 
   ELSE
      tracer_ukca_um(:,:,:,tr_lookup(i)) = all_tracers(:,:,:,i) 
   END IF 
END DO

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END SUBROUTINE ukca_all_tracers_copy_out

SUBROUTINE ukca_make_tr_map()

! Create a mapping array from UKCA tracer order to the correct
! order for all_tracers. Called once and values saved
!
! all_tracers has the chemistry items first, 
! all in the same order as advt, and then the MODE
! items in STASH code order. 
! Other tracers e.g. age of air come last
! For BE schemes, the chemistry items are in stashcode order and q isn't
! included.
! For NR schemes the chemistry tracers are not always in stashcode
! order and q may or may not be included
! Note that advt(i) is only dimensioned to jpctr - i.e. the chemistry
! 'tracer' items. 

USE ukca_mode_setup, ONLY: mode_names
USE ukca_nmspec_mod,    ONLY: nm_spec, nm_spec_active
! the length of nm_spec is a_max_ukcavars
USE ukca_tracer_stash,     ONLY: a_max_ukcavars, a_ukca_first, a_ukca_last
! tr_ukca_a is a logical array defining which of the 
! tracers in nm_spec are on
USE ukca_option_mod, ONLY: tr_ukca_a, jpctr
USE ukca_d1_defs, ONLY: n_use_tracers, n_mode_tracers
USE ukca_cspecies, ONLY: n_age

! also information needed on the chemistry scheme
USE asad_mod,                ONLY: advt

USE ereport_mod,     ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr
USE yomhook,                   ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE


CHARACTER (LEN=*), PARAMETER   :: RoutineName = 'UKCA_MAKE_TR_MAP'
CHARACTER (LEN=errormessagelength)            :: cmessage   ! Error message
INTEGER                        :: errcode    ! Variable passed to ereport
INTEGER :: i, j, tr_index, icount_mode, icount_other

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! create arrays to hold indices in tracer_ukca for each entry 
! of all_tracers and names for each entry in all_tracers
ALLOCATE (tr_lookup(n_use_tracers))
ALLOCATE (tr_names(n_use_tracers))

! Set tr_lookup to -999 to check for missing entries 
tr_lookup(:) = -999

! initialise counters for MODE tracers and other
icount_mode=0
icount_other=0

! Loop over advected species. 
DO i=1, SIZE(nm_spec_active)

    ! If in chemistry array, set to be in same place as in advt
    IF (ANY (advt(:) == nm_spec_active(i))) THEN
       tr_index = arrayloc(advt, nm_spec_active(i))
       ! check that not too big (<= jpctr)
       IF (tr_index < 1 .OR. tr_index > jpctr) THEN
          WRITE(umMessage,'(A)') 'Problem with chemistry index'
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I4)') 'tracer index : ', tr_index
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I4, A, I4)') 'Should be between:', &
             1, ' and ', jpctr
          CALL umPrint(umMessage,src=routinename)
          cmessage='Problem with chemistry index'
          errcode = 1
          CALL ereport(routinename,errcode,cmessage)
       END IF 

    ! If the first 7 characters of the name match one of the UKCA
    ! modes this is a MODE tracer -
    ! slot these incrementally into MODE section
    ELSE IF (ANY(mode_names(:) == nm_spec_active(i)(1:7))) THEN
       icount_mode = icount_mode  + 1
       tr_index = jpctr + icount_mode
       ! check that not too big or small (>jpctr <= jpctr+n_mode_Tracers)
       IF (tr_index <= jpctr .OR. tr_index > jpctr + n_mode_tracers) THEN
          WRITE(umMessage,'(A)') 'Problem with MODE index'
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I4)') 'tracer index : ', tr_index
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I4, A, I4)') 'Should be between:', &
             jpctr + 1, ' and ', jpctr + n_mode_tracers
          CALL umPrint(umMessage,src=routinename)
          cmessage='Problem with MODE index'
          errcode = 2
          CALL ereport(routinename,errcode,cmessage)
       END IF 
    ELSE
       icount_other = icount_other  + 1
       tr_index= jpctr + n_mode_tracers + icount_other
       ! check that not too big or small (>jpctr+n_mode_Tracers <=n_use_tracers)
       IF (tr_index <= jpctr + n_mode_tracers .OR. &
               tr_index > n_use_tracers) THEN
          WRITE(umMessage,'(A)') 'Problem with other tracer index'
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I6, 1X, I6)') 'i, tracer index : ', i, tr_index
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A,A10)') 'tracer name : ', nm_spec_active(i)
          CALL umPrint(umMessage,src=routinename)
          WRITE(umMessage,'(A, I6, A, I6)') 'Should be between:', &
             jpctr + n_mode_tracers + 1, ' and ', n_use_tracers
          CALL umPrint(umMessage,src=routinename)
          cmessage='Problem with other tracer index'
          errcode = 3
          CALL ereport(routinename,errcode,cmessage)
       END IF 
       ! Set age of air index
       IF(nm_spec_active(i) == 'AGE OF AIR') n_age = tr_index
    END IF 

    ! Set the index and name using tr_index from the above
    ! block
    tr_lookup(tr_index) = i
    tr_names(tr_index) = nm_spec_active(i)

END DO

! If H2O in advt set the right slot for that to -1
IF (ANY(advt(:) == 'H2O')) THEN
     tr_index=arrayloc(advt(:), 'H2O')
     tr_lookup(tr_index) = -1
     tr_names(tr_index) = 'H2O'
END IF

! Finally check no entries in tr_lookup are -999
IF (ANY(tr_lookup(:) == -999)) THEN
    cmessage='Failed to set tr_lookup'
    errcode=1
    CALL ereport(routinename,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END SUBROUTINE ukca_make_tr_map

INTEGER FUNCTION arrayloc(ch_arr,ch_test)

! takes as input a 1D CHARACTER arrray and 
! a CHARACTER variable
! If the value is in the array, return the index,
! if it isn't return -99
! Note that the match is case sensitive but ignores
! all leading and trailing white space
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i
CHARACTER(LEN=*),INTENT(IN) :: ch_arr(:)
CHARACTER(LEN=*),INTENT(IN) :: ch_test

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: routinename='ARRAYLOC'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

! start with default value
arrayloc = -99 

DO i = 1, SIZE (ch_arr)
   IF (TRIM(ADJUSTL(ch_arr(i))) == TRIM(ADJUSTL(ch_test))) THEN
     arrayloc = i
     EXIT
   END IF 
END DO 

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

END FUNCTION arrayloc

END MODULE ukca_all_tracers_copy_mod
