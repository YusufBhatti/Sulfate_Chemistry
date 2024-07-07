#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_grid(l_fullset)
!
! Description: This routine sets up the grid and transient
!              field definitions used when coupling using
!              OASIS3-MCT. All interface integers are
!              expected to be 32 bits.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!==================================================================

! OASIS3-MCT PSMILe subroutines followed by
! essential control variables
USE mod_prism, ONLY: prism_def_partition_proto, &
                     prism_def_var_proto, &
                     prism_enddef_proto, &
                     prism_out,       &
                     prism_real,      &
                     prism_in

USE oasis_grad_index_mod, ONLY:   oasis_grad_index

USE oasis_atm_data_mod, ONLY:     icpl,            &
                                  icpl_atm,        &
                                  partitionbox,    &
                                  partitionserial, &
                                  fld_type_u,      &
                                  fld_type_p,      &
                                  fld_type_v,      &
                                  prism_nsec,      &
                                  transient_a2o, transient_o2a,  &
                                  transient_a2c, transient_c2a,  &
                                  transient_a2o_count, transient_o2a_count,  &
                                  transient_a2c_count, transient_c2a_count

USE um_types
USE um_parvars
USE um_parparams, ONLY: halo_type_no_halo
USE Field_Types,  ONLY: Nfld_Max
USE ereport_mod, ONLY: ereport

USE coupling_control_mod,  ONLY: l_senior, l_junior
USE oasis_timers, ONLY: l_oasis_timers, gridstarts, gridends
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE



LOGICAL, INTENT(IN) :: l_fullset ! Indicates whether to perform
                        ! a null set of calls for processes which are
                        ! ultimately not involved in coupling, or the
                        ! full range of oasis set-up for processes
                        ! which are involved in coupling.
                        ! TRUE gives you the full set.

! Local variables...
INTEGER (KIND=integer32) :: ierror ! Return code for OASIS3-MCT calls

! This module expects to define an Arakawa C grid:
!   - this is an instance of PRISM_reglonlatvrt
!   - 3D grid
!   - points are created for the first time here
!   - each grid cell is a cube w/ 8 vertices
!   - C grid has t, u, v points

INTEGER (KIND=integer64), PARAMETER :: nGridDims = 3

INTEGER (KIND=integer32) :: valid_shape(2,nGridDims,Nfld_Max)

INTEGER (KIND=integer32) :: var_nodims(2)
INTEGER (KIND=integer32) :: var_shape(4)

INTEGER (KIND=integer32) :: k, tc
INTEGER (KIND=integer32) :: coupling_intent

INTEGER (KIND=integer32), ALLOCATABLE :: il_paral(:)
! Decomposition for each proc
INTEGER :: ig_parsize    ! Size of array decomposition

INTEGER (KIND=integer32):: partition_id_tu ! Local T/U partition ID
INTEGER (KIND=integer32):: partition_id_v  ! Local V partition ID
INTEGER (KIND=integer32):: partition_id_null  ! Local null partition ID

INTEGER (KIND=integer32):: partition

CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=*) :: RoutineName

INTEGER (KIND=integer32):: type_index  ! Index to grid point type

LOGICAL (KIND=logical32):: l_2nd_order_reqd

PARAMETER (RoutineName='OASIS3_GRID')

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (l_oasis_timers) CALL gridstarts

l_2nd_order_reqd = .FALSE.

IF (l_fullset) THEN

  ! Initialise all valid shape components to 1 for all
  ! dimensions and grid point types.
  valid_shape(:,:,:) = 1



  DO type_index = 1, Nfld_Max

    ! Valid shapes are the local dimsions for OASIS3-MCT
    valid_shape(2,1,type_index) = lasize(1,type_index,halo_type_no_halo)
    valid_shape(2,2,type_index) = lasize(2,type_index,halo_type_no_halo)
    valid_shape(2,3,type_index) = 1

  END DO

  ierror = 0  ! Initialise error flag to zero because some OASIS
              ! routines do not do this explicitily which
              ! can cause us problems.

  ! Each PE couples in parallel, hence they all need to call this.
  ig_parsize = 5
  ALLOCATE(il_paral(ig_parsize))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NOTE: In the following we only define two different
  ! OASIS partitions, despite the fact that the rest of our code
  ! caters for the possibility of there being three
  ! separate grid point types (T, U and V.)
  ! The underlying assumption is that T and U points always use
  ! the same partition regardless of dynamical core choice.
  ! That's true for ND and Endgame but may not be always thus.
  ! However, this is assumed at this point because of the extra
  ! complexity a third separate grid would bring to the OASIS
  ! control files (e.g. namcouple contents, grids and rmp files etc.)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set up partitions in TWO dimensions NOT 1 as in the
  ! standard OASIS3-MCT toy models, to allow visualisation.

  ! See OASIS3-MCT user guide for the meanings of the following
  ! array components and possible alternatives.
  ! 1 = partition type (we use box)
  ! 2 = data offset
  ! 3 = local x dimension
  ! 4 = local y dimension
  ! 5 = global x dimension

  il_paral(1) = PartitionBox
  il_paral(2) = (datastart(2)-1)*glsize(1,fld_type_p)+datastart(1)-1
  il_paral(3) = valid_shape(2,1,fld_type_p)
  il_paral(4) = valid_shape(2,2,fld_type_p)
  il_paral(5) = glsize(1,fld_type_p)

  ! We have to define two separate partitions -
  ! The first one for the T and U points, the second one for
  ! the V points which has one row fewer.
  CALL PRISM_def_partition_proto (partition_id_tu, il_paral, ierror)

  IF (ierror /= 0 ) THEN 
    Cmessage = "PRISM_def_partition failure TU (slave)"
    CALL Ereport("OASIS3_GRID", ierror, Cmessage)
  END IF

  il_paral(4) = valid_shape(2,2,fld_type_v)

  CALL PRISM_def_partition_proto (partition_id_v, il_paral, ierror)

  IF (ierror /= 0 ) THEN 
    Cmessage = "PRISM_def_partition failure V (slave)"
    CALL Ereport("OASIS3_GRID", ierror, Cmessage)
  END IF


  DEALLOCATE(il_paral)

  var_nodims(1) = 2
  var_nodims(2) = 1

  var_shape(1) = 1
  var_shape(3) = 1

  !==================================================================
  ! Now we define outgoing fields for ocean/ice components:
  !==================================================================
  DO tc = 1, transient_a2o_count

    IF (transient_a2o(tc)%indx > 0) THEN

      IF (transient_a2o(tc)%grid == "U") THEN
        var_shape(2) = valid_shape(2,1,fld_type_u)
        var_shape(4) = valid_shape(2,2,fld_type_u)
        partition = partition_id_tu
      ELSE
        IF (transient_a2o(tc)%grid == "V") THEN
          var_shape(2) = valid_shape(2,1,fld_type_v)
          var_shape(4) = valid_shape(2,2,fld_type_v)
          partition = partition_id_v
        ELSE
          var_shape(2) = valid_shape(2,1,fld_type_p)
          var_shape(4) = valid_shape(2,2,fld_type_p)
          partition = partition_id_tu
        END IF
      END IF

      ! Define output field.
      WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define ',tc,         &
            transient_a2o(tc)%NAME,' for ',PRISM_Out
      CALL umPrint(umMessage,src='oasis3_grid')

      ! If this field requires special 2nd order processing, 
      ! set the flag to tell us to do some calculations later on.
      IF (transient_a2o(tc)%order == 2) l_2nd_order_reqd = .TRUE.

      !-------------------------------------------------
      ! Loop over vertical levels (slice), if any.
      !-------------------------------------------------
      DO k=1,transient_a2o(tc)%n_slices
        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define ',tc,  &
             transient_a2o(tc)%name(k),' for ',PRISM_Out
        CALL umPrint(umMessage,src='oasis3_grid')
        !-------------------------------------------------
        ! Determine the oasis_id
        !-------------------------------------------------
        CALL PRISM_def_var_proto(transient_a2o(tc)%oasis_id(k),    &
                      transient_a2o(tc)%name(k), partition,        &
                      var_nodims, PRISM_Out,  var_shape,           &
                      PRISM_Real, ierror )
        !-------------------------------------------------
        ! Check for error
        !-------------------------------------------------
        IF (ierror /= 0 ) THEN
          Cmessage = "PRISM_def_var failure for OUTPUT " // &
                       transient_a2o(tc)%name(k)
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF
      END DO
    END IF
  END DO  ! tc over a2o fields


  !==================================================================
  ! Now we can define incoming ocean/seaice fields:
  !==================================================================
  DO tc = 1, transient_o2a_count

    IF (transient_o2a(tc)%indx > 0) THEN
      IF (transient_o2a(tc)%grid == "U") THEN
        var_shape(2) = valid_shape(2,1,fld_type_u)
        var_shape(4) = valid_shape(2,2,fld_type_u)
        partition = partition_id_tu
      ELSE
        IF (transient_o2a(tc)%grid == "V") THEN
          var_shape(2) = valid_shape(2,1,fld_type_v)
          var_shape(4) = valid_shape(2,2,fld_type_v)
          partition = partition_id_v
        ELSE
          var_shape(2) = valid_shape(2,1,fld_type_p)
          var_shape(4) = valid_shape(2,2,fld_type_p)
          partition = partition_id_tu
        END IF
      END IF

      !-------------------------------------------------
      ! Loop over vertical levels (slice), if any.
      !-------------------------------------------------
      DO k=1,transient_o2a(tc)%n_slices
        !-------------------------------------------------
        ! Write slice name to output
        !-------------------------------------------------
        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define ',TC,  &
            transient_o2a(tc)%name(k),' for ',prism_in
        CALL umPrint(umMessage,src='oasis3_grid')
        !-------------------------------------------------
        ! Determine the oasis_id
        !-------------------------------------------------
        CALL PRISM_def_var_proto (transient_o2a(tc)%oasis_id(k),    &
                  transient_o2a(tc)%name(k), partition,             &
                  var_nodims, PRISM_In, var_shape,                  &
                  PRISM_Real, ierror )
        !-------------------------------------------------
        ! Check for error
        !-------------------------------------------------
        IF (ierror /= 0 ) THEN
          Cmessage = "PRISM_def_var failure for INPUT " // &
                      transient_o2a(tc)%name(k)
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF
      END DO
    END IF
  END DO  ! Over tc


  !==================================================================
  ! Now we deal with defining any chemistry coupling fields
  ! Remember, that if this is the senior model A2C fields are output
  ! and C2A are input. If this is the Junior model A2C are input and
  ! C2A are output.
  !==================================================================

  IF (l_junior) THEN
     coupling_intent = PRISM_In
  ELSE IF (l_senior) THEN
     coupling_intent = PRISM_Out
  END IF

  DO tc = 1, transient_a2c_count

    IF (transient_a2c(tc)%indx > 0) THEN

      IF (transient_a2c(tc)%grid == "U") THEN
        var_shape(2) = valid_shape(2,1,fld_type_u)
        var_shape(4) = valid_shape(2,2,fld_type_u)
        partition = partition_id_tu
      ELSE
        IF (transient_a2c(tc)%grid == "V") THEN
          var_shape(2) = valid_shape(2,1,fld_type_v)
          var_shape(4) = valid_shape(2,2,fld_type_v)
          partition = partition_id_v
        ELSE
          var_shape(2) = valid_shape(2,1,fld_type_p)
          var_shape(4) = valid_shape(2,2,fld_type_p)
          partition = partition_id_tu
        END IF
      END IF

      !-------------------------------------------------
      ! If this field requires special 2nd order processing, 
      ! set the flag to tell us to do some calculations later on.
      !-------------------------------------------------
      IF (transient_a2c(tc)%order == 2) l_2nd_order_reqd = .TRUE.

      !-------------------------------------------------
      ! Loop over vertical levels (slice), if any
      !-------------------------------------------------
      DO k=1,transient_a2c(tc)%n_slices
        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define ',tc,  &
             transient_a2c(tc)%name(k),' for ',coupling_intent
        CALL umPrint(umMessage,src='oasis3_grid')
        !-------------------------------------------------
        ! Determine the oasis_id
        !-------------------------------------------------
        CALL PRISM_def_var_proto(transient_a2c(tc)%oasis_id(k),  &
                      transient_a2c(tc)%name(k), partition,        &
                      var_nodims, coupling_intent,  var_shape,     &
                      PRISM_Real, ierror )
        !-------------------------------------------------
        ! Check for error
        !-------------------------------------------------
        IF (ierror /= 0 ) THEN
          Cmessage = "PRISM_def_var failure for OUTPUT " // &
                       transient_a2c(tc)%name(k) 
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF
      END DO
    END IF
  END DO  ! tc over a2c fields

  IF (l_junior) THEN
     coupling_intent = PRISM_Out
  ELSE IF (l_senior) THEN
     coupling_intent = PRISM_In
  END IF

  DO tc = 1, transient_c2a_count

    IF (transient_c2a(tc)%indx > 0) THEN

      IF (transient_c2a(tc)%grid == "U") THEN
        var_shape(2) = valid_shape(2,1,fld_type_u)
        var_shape(4) = valid_shape(2,2,fld_type_u)
        partition = partition_id_tu
      ELSE
        IF (transient_c2a(tc)%grid == "V") THEN
          var_shape(2) = valid_shape(2,1,fld_type_v)
          var_shape(4) = valid_shape(2,2,fld_type_v)
          partition = partition_id_v
        ELSE
          var_shape(2) = valid_shape(2,1,fld_type_p)
          var_shape(4) = valid_shape(2,2,fld_type_p)
          partition = partition_id_tu
        END IF
      END IF

      !-------------------------------------------------
      ! If this field requires special 2nd order processing, 
      ! set the flag to tell us to do some calculations later on.
      !-------------------------------------------------
      IF (transient_c2a(tc)%order == 2) l_2nd_order_reqd = .TRUE.

      !-------------------------------------------------
      ! Loop over verical levels (slice), if any.
      !-------------------------------------------------
      DO k=1,transient_c2a(tc)%n_slices
        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define ',tc,  &
             transient_c2a(tc)%name(k),' for ',coupling_intent
        CALL umPrint(umMessage,src='oasis3_grid')
        !-------------------------------------------------
        ! Determine the oasis_id
        !-------------------------------------------------
        CALL PRISM_def_var_proto(transient_c2a(tc)%oasis_id(k),  &
                      transient_c2a(tc)%name(k), partition,      &
                      var_nodims, coupling_intent,  var_shape,   &
                      PRISM_Real, ierror )
        !-------------------------------------------------
        ! Check for error
        !-------------------------------------------------
        IF (ierror /= 0 ) THEN
          Cmessage = "PRISM_def_var failure for OUTPUT " // &
                       transient_c2a(tc)%name(k) 
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF
      END DO
    END IF
  END DO  ! tc over c2a fields


  !=================================================================
  ! NOTE: 
  ! If your name here doesn't correspond with the name in namcouple
  ! it will not necessarily tell you, or crash. However this is 
  ! mitigated against by having the UM extract names direct from the
  ! namcouple in normal model configurations.  
  !=================================================================

  ! Finish the PRISM definition phase and perform inter component
  ! integrity checking.

  CALL PRISM_enddef_proto(ierror)

  IF (ierror /= 0 ) THEN 
    Cmessage = "PRISM_enddef failure" 
    CALL Ereport("OASIS3_GRID", ierror, Cmessage)
  END IF

  ! Set the start time. For OASIS3-MCT this is simply based on the
  ! number of seconds this run, so at the start of the
  ! run it will be zero regardless of the true model date
  ! It may be that we need to do some cross checking with other
  ! components in a coupled model system, in which case this
  ! is one suitable point to employ such a test.
  PRISM_nsec = 0


  ! Set up indexing and distances for use in gradient calculation
  ! if 2nd order conservative coupling is employed with OASIS3-MCT
  IF (l_2nd_order_reqd) THEN
    CALL OASIS_Grad_Index()
  END IF

ELSE   ! l_fullset = FALSE

  ! We need all the non-atmos processes to perform a limited set of
  ! prism routines for initialisation purposes. The positioning of these
  ! calls is critical to avoid deadlock.
  IF ( icpl /= icpl_atm ) THEN

    ierror = 0  ! Initialise error flag to zero because some OASIS
                ! routines do not do this explicitly which
                ! can cause us problems.

    ig_parsize = 3
    ALLOCATE(il_paral(ig_parsize))

    ! OASIS3-MCT documentation asserts that we should set
    ! il_paral(:) = 0, which is all well and good, but it fails
    ! to tell us what size this array needs to be for the null case!
    ! Here we assume serial is OK, despite the fact that for
    ! processes actually involved in coupling (see above)
    ! we use a Box partition.

    il_paral(1) = PartitionSerial

    ! All the partition dimensions are zero
    il_paral(2:3) = 0

    ! Although we have to define two separate partitions in the main code
    ! for our null subset case we can use the same partition for all
    ! def_var calls.
    CALL PRISM_def_partition_proto (partition_id_null, il_paral, ierror)

    IF (ierror /= 0 ) THEN 
      Cmessage = "PRISM_def_partition failure (IOS)" 
      CALL Ereport("OASIS3_GRID", ierror, Cmessage)
    END IF

    DEALLOCATE(il_paral)

    ! Now we define outgoing fields:
    var_nodims(1) = 2
    var_nodims(2) = 1

    ! All field dimensions are zero.
    var_shape(1:4) = 0

    ! Null def_var calls actually don't need the names to correspond
    ! to the main calls, but since we have them avaialable to us, we
    ! do things "properly" anyway in that respect. We don't need to
    ! distinguish between fields on U, T and V grids so we can just
    ! use our null partition for everything
    partition = partition_id_null

    ! Now we define outgoing fields:
    DO tc = 1, transient_a2o_count

      DO k=1,transient_a2o(tc)%n_slices

        ! Define output field.
        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define Null ',tc,   &
            transient_a2o(tc)%NAME(k), ' for ',PRISM_Out
        CALL umPrint(umMessage,src='oasis3_grid')

        CALL PRISM_def_var_proto(transient_a2o(tc)%oasis_id(k),   &
                  transient_a2o(tc)%NAME(k), partition,          &
                  var_nodims, PRISM_Out,  var_shape,          &
                  PRISM_Real, ierror )

        IF (ierror /= 0 ) THEN 
          Cmessage = "PRISM_def_var failure (IOS) for OUTPUT " &
                            // transient_a2o(tc)%NAME(k)
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF

      END DO  ! Over K

    END DO  ! Over tc

    ! Now we define incoming fields:
    DO tc = 1, transient_o2a_count

      DO k=1,transient_o2a(tc)%n_slices

        WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define Null ',tc,  &
            transient_o2a(tc)%NAME(k),' for ',prism_in
        CALL umPrint(umMessage,src='oasis3_grid')

        CALL PRISM_def_var_proto (transient_o2a(tc)%oasis_id(k),    &
                transient_o2a(tc)%NAME(k), partition,              &
                var_nodims, PRISM_In, var_shape,               &
                PRISM_Real, ierror )

        IF (ierror /= 0 ) THEN 
          Cmessage = "PRISM_def_var failure (IOS) for INPUT " &
                            // transient_o2a(tc)%NAME(k)
          CALL Ereport("OASIS3_GRID", ierror, Cmessage)
        END IF

      END DO  ! Over k

    END DO  ! Over tc


    !==================================================================
    ! Now we deal with any chemistry coupling fields if required.
    ! Remember, that if this is the senior model A2C fields are output
    ! and C2A are input. If this is the Junior model A2C are input and
    ! C2A are output.
    !==================================================================

    IF (l_junior) THEN
       coupling_intent = PRISM_In
    ELSE
       coupling_intent = PRISM_Out
    END IF

    DO tc = 1, transient_a2c_count

      IF (transient_a2c(tc)%indx > 0) THEN

        !-------------------------------------------------
        ! Loop over vertical levels (slice)
        !-------------------------------------------------
        DO k=1,transient_a2c(tc)%n_slices
          !-------------------------------------------------
          ! Name needs changing for greater than first level
          !-------------------------------------------------
          IF (k > 1) THEN
            WRITE(transient_a2c(tc)%name(k)(6:8),'(I3.3)') k
          END IF
          !-------------------------------------------------
          ! Define output field.
          !-------------------------------------------------
          WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define Null ', &
            tc, transient_a2c(tc)%name(k), ' for ',coupling_intent
          CALL umPrint(umMessage,src='oasis3_grid')
          !-------------------------------------------------
          ! Determine the oasis_id
          !-------------------------------------------------
          CALL PRISM_def_var_proto(transient_a2c(tc)%oasis_id(k),  &
                 transient_a2c(tc)%name(k), partition,             &
                 var_nodims, coupling_intent,  var_shape,                &
                 PRISM_Real, ierror )
          !-------------------------------------------------
          ! Check for error
          !-------------------------------------------------
          IF (ierror /= 0 ) THEN
            Cmessage = "PRISM_def_var failure (IOS) for OUTPUT " &
                       // transient_a2c(tc)%name(k)
            CALL Ereport("OASIS3_GRID", ierror, Cmessage)
          END IF
        END DO
      END IF
    END DO  ! Over tc

    IF (l_junior) THEN
       coupling_intent = PRISM_Out
    ELSE
       coupling_intent = PRISM_In
    END IF

    DO tc = 1, transient_c2a_count

      IF (transient_c2a(tc)%indx > 0) THEN

        !-------------------------------------------------
        ! Loop over vertical levels (slice)
        !-------------------------------------------------
        DO k=1,transient_c2a(tc)%n_slices
          !-------------------------------------------------
          ! Name needs changing for greater than first level
          !-------------------------------------------------
          IF (k > 1) THEN
            WRITE(transient_c2a(tc)%name(k)(6:8),'(I3.3)') k
          END IF
          !-------------------------------------------------
          ! Define output field.
          !-------------------------------------------------
          WRITE(umMessage,'(A,1X,I4,1X,A,1X,A,1X,I6)') 'Define Null ', &
            tc, transient_c2a(tc)%name(k), ' for ',coupling_intent
          CALL umPrint(umMessage,src='oasis3_grid')
          !-------------------------------------------------
          ! Determine the oasis_id
          !-------------------------------------------------
          CALL PRISM_def_var_proto(transient_c2a(tc)%oasis_id(k),  &
                 transient_c2a(tc)%name(k), partition,             &
                 var_nodims, coupling_intent,  var_shape,                &
                 PRISM_Real, ierror )
          !-------------------------------------------------
          ! Check for error
          !-------------------------------------------------
          IF (ierror /= 0 ) THEN
            Cmessage = "PRISM_def_var failure (IOS) for OUTPUT " &
                       // transient_c2a(tc)%name(k)
            CALL Ereport("OASIS3_GRID", ierror, Cmessage)
          END IF
        END DO
      END IF
    END DO  ! Over tc


    ! Finish the PRISM definition phase and perform inter component
    ! integrity checking.

    CALL PRISM_enddef_proto(ierror)

    IF (ierror /= 0 ) THEN 
      Cmessage = "PRISM_enddef failure (IOS)" 
      CALL Ereport("OASIS3_GRID", ierror, Cmessage)
    END IF

  END IF ! IOS procs only.

END IF ! Null case for non-coupling procs.

IF (l_oasis_timers) CALL gridends
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE OASIS3_grid
#endif
