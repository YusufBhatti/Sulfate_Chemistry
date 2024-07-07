! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash.

MODULE close_unneeded_stash_files_mod

! Description:
! Closes files which no longer need to be open, because all the required stash
! entries have already been written to them.

! Method:
! loop over the stashlist to find all the stash entries written to a file. If
! every entry has been written for every model step to a given file, close and
! reset it.

IMPLICIT NONE
PRIVATE

PUBLIC ::                                                                      &
  close_unneeded_stash_files,                                                  &
  close_unneeded_stash_files_freq,                                             &
  check_close_unneeded_stash_files

INTEGER, PARAMETER :: close_unneeded_stash_files_freq = 1
LOGICAL, PARAMETER :: check_close_unneeded_stash_files = .TRUE.

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName='CLOSE_UNNEEDED_STASH_FILES_MOD'

INTEGER, ALLOCATABLE :: stlist_il(:)

!------------------------------------------------------------------------------!
CONTAINS
!------------------------------------------------------------------------------!

SUBROUTINE close_unneeded_stash_files

USE umPrintMgr, ONLY:                                                          &
    ummessage, umPrint, printstatus, prstatus_oper

USE file_manager, ONLY:                                                        &
    init_file_loop, um_file_type

USE stash_array_mod, ONLY:                                                     &
    totitems, stlist, nsttims, sttabl

USE model_time_mod, ONLY:                                                      &
    stepim

USE submodel_mod, ONLY:                                                        &
    atmos_im

USE stparam_mod, ONLY:                                                         &
    st_freq_code, st_end_of_list, st_end_time_code, st_infinite_time,          &
    st_proc_no_code, st_output_code, st_output_type, st_fieldsfile, st_netcdf

USE yomhook, ONLY:                                                             &
    lhook, dr_hook

USE parkind1, ONLY:                                                            &
    jprb, jpim

USE model_file, ONLY:                                                          &
    model_file_close

USE um_parcore, ONLY:                                                          &
    mype

USE umnetcdf_mod, ONLY:                                                        &
    nc_file_close

IMPLICIT NONE

TYPE(um_file_type), POINTER ::                                                 &
  pp_file

LOGICAL ::                                                                     &
  closable_unit, found_nonfile, found_closable

LOGICAL, ALLOCATABLE ::                                                        &
  closable_units_array(:)

INTEGER ::                                                                     &
  output_code, output_type, il, it, i, ntab, il_size

INTEGER, ALLOCATABLE ::                                                        &
  used_units(:), stlist_il_new(:), output_type_array(:)

INTEGER, PARAMETER :: not_file = -1

INTEGER, PARAMETER :: num_output_file_types = 2

! Dr-hook
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CLOSE_UNNEEDED_STASH_FILES"
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//":"//RoutineName,zhook_in,zhook_handle)

IF (printstatus>=prstatus_oper) THEN
  WRITE(ummessage,'(A)') "Checking for stash output files which" //            &
                         " are no longer required"
  CALL umprint(ummessage)
END IF

! Create a local stashlist, which can be trimmed as the corresponding files are
! closed

IF (.NOT.ALLOCATED(stlist_il)) THEN

  ALLOCATE(stlist_il(totitems))

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(il)                                                             &
!$OMP& SHARED(stlist_il, totitems)
  DO il =1,totitems

    stlist_il(il)=il

  END DO
!$OMP END PARALLEL DO

END IF

il_size = SIZE(stlist_il)

ALLOCATE(closable_units_array(il_size))
ALLOCATE(used_units(il_size))
ALLOCATE(output_type_array(il_size))

! loop over the stashlist:
!
!  * generate a list of units output to by each stash entry. (used_units)
!
!  * generate a list of logicals indicating which stash entries have been
!    written for every required step (closable_units_array)
!
! if a unit is in used_units, but every corresponding entry in
! closable_units_array is .TRUE., it can be safely closed

found_closable = .FALSE.
found_nonfile = .FALSE.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(il, it, closable_unit, output_code, output_type, ntab)          &
!$OMP& SHARED(used_units, closable_units_array, il_size, stlist, stepim,       &
!$OMP&        nsttims,sttabl,stlist_il,output_type_array) &
!$OMP& REDUCTION(.OR.:found_closable,found_nonfile)
DO il =1,il_size

  output_code = stlist(st_output_code,stlist_il(il))
  output_type = stlist(st_output_type,stlist_il(il))

  IF (output_code >=  0) THEN
    ! Not output to file

    used_units(il) = not_file
    output_type_array(il) = not_file
    closable_units_array(il) = .FALSE.
    found_nonfile = .TRUE.

  ELSE
    ! Output to file

    used_units(il) = -output_code
    output_type_array(il) = output_type

    closable_unit = .FALSE.

    ! Check diagnostics which are not active
    IF (stlist(st_proc_no_code,stlist_il(il)) == 0) THEN

      ! Stash diagnostic is not active
      closable_unit = .TRUE.

    ELSE
      IF (stlist(st_freq_code,stlist_il(il)) > 0) THEN
        ! Stash diagnostic is active; timesteps are from a sequence

        ! Check if we have passed the last output step
        IF (stepim(atmos_im) > stlist(st_end_time_code,stlist_il(il))          &
            .AND. .NOT. stlist(st_end_time_code,stlist_il(il))                 &
                                                          == st_infinite_time) &
        THEN
          ! we are past the last output step
          closable_unit = .TRUE.
        END IF

      ELSE IF (stlist(st_freq_code,stlist_il(il)) <  0) THEN
        ! Stash diagnostic is active; timesteps are from a list

        ! Find the list of output steps
        ntab=-stlist(st_freq_code,stlist_il(il))

        ! assume we are past the last output step,
        ! until we find an output step in the future
        closable_unit = .TRUE.
        DO it=1,nsttims
          ! did we get to the end of the output step list already?
          IF (sttabl(it,ntab) == st_end_of_list) EXIT

          ! check if this output step is in the future
          IF (stepim(atmos_im) <= sttabl(it,ntab)) THEN
            ! output step is in the future
            closable_unit = .FALSE.
            EXIT
          END IF
        END DO
      END IF
    END IF

    IF (closable_unit) THEN
      closable_units_array(il) = .TRUE.
      found_closable = .TRUE.
    ELSE
      closable_units_array(il) = .FALSE.
    END IF

  END IF

END DO
!$OMP END PARALLEL DO

! If we found any closable units, now close them.

IF (found_closable) THEN

  DO i = 1, num_output_file_types
    NULLIFY(pp_file)
    SELECT CASE(i)
      CASE(st_fieldsfile)
        pp_file => init_file_loop(handler="portio")
        output_type = st_fieldsfile
      CASE(st_netcdf)
        pp_file => init_file_loop(handler="netcdf")
        output_type = st_netcdf
    END SELECT

    DO WHILE (ASSOCIATED(pp_file))

      IF (printstatus>=prstatus_oper) THEN
        WRITE(ummessage,'(A)') "Found file: " // pp_file % filename
        CALL umprint(ummessage)
      END IF

      ! Are there any units in used_units corresponding to this one?
      IF ( ANY(used_units==pp_file % unit .AND. &
               output_type_array==output_type) ) THEN

        ! Select out all the corresponding entries in closable_units_array
        ! If they are all .TRUE. we can safely close the file
        IF ( ALL(PACK(closable_units_array,used_units==pp_file % unit .AND.    &
                      output_type_array==output_type)) ) THEN

          IF (printstatus>=prstatus_oper) THEN
            WRITE(ummessage,'(A)') 'file "'//pp_file%filename//'" is closable'
            CALL umprint(ummessage)
          END IF

          IF (output_type == st_netcdf) THEN
            IF (mype == 0) THEN
              CALL nc_file_close(pp_file)
            END IF
          ELSE
            CALL model_file_close(pp_file % unit, pp_file % filename)
          END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(il)                                                             &
!$OMP& SHARED(pp_file, stlist,stlist_il, il_size, output_type)
          DO il =1,il_size
            IF (-stlist(st_output_code,stlist_il(il)) == pp_file % unit .AND.  &
                 stlist(st_output_type,stlist_il(il)) == output_type) THEN
              ! Now reset this stashlist entry's output_code to zero, 
              ! so we don't try to close an already closed unit
              stlist(st_output_code,stlist_il(il)) = 0
            END IF
          END DO
!$OMP END PARALLEL DO

        END IF

      END IF

      ! Increment file loop pointer for next iteration
      pp_file => pp_file % next

    END DO

  END DO

END IF

! If we found any entries in the stashlist which don't correspond to files,
! we can now trim them out.

IF (found_nonfile) THEN
  il_size = COUNT(used_units/=not_file)
  ALLOCATE(stlist_il_new(il_size))
  stlist_il_new = PACK(stlist_il,used_units/=not_file)
  CALL MOVE_ALLOC(stlist_il_new,stlist_il)
END IF

DEALLOCATE(used_units)
DEALLOCATE(output_type_array)
DEALLOCATE(closable_units_array)

IF (lhook) CALL dr_hook(ModuleName//":"//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE close_unneeded_stash_files

!------------------------------------------------------------------------------!

END MODULE close_unneeded_stash_files_mod

