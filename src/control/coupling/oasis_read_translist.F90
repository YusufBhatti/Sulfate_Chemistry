! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_read_translist()

!
! Description:
! Read translist to help us set up a full list of
! fields involved in coupling.
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!==================================================================

USE oasis_atm_data_mod, ONLY: transcount,    &
                              c_in,          &
                              n_in,          &
                              c_out,         &
                              n_out,         &
                              i_number,      &
                              f_id,          &
                              map_type,      &
                              transfld,      &
                              max_lev,       &
                              max_transients, &
                              transient_a2o, transient_o2a, &
                              transient_a2c, transient_c2a, &
                              transient_c2a_count, transient_a2c_count, &
                              transient_o2a_count, transient_a2o_count, &
                              ocn_component, atm_component, jnr_component

USE oasis_operations_mod, ONLY: read_nml_transcount, &
                                read_nml_transfld

USE umPrintMgr, ONLY:  umPrint, umMessage
USE ereport_mod, ONLY: ereport
USE um_parcore, ONLY: mype
USE coupling_control_mod,  ONLY: l_junior, l_senior
USE Field_Types, ONLY: fld_type_u, fld_type_v

IMPLICIT NONE

! Note: since this routine gets called before the file manager is 
!       initialised, we have to use a non-reservable unit
INTEGER, PARAMETER :: lb=8
INTEGER :: iostatus
INTEGER :: indx_a2c, indx_a2o, indx_c2a, indx_o2a
INTEGER :: k

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_READ_TRANSLIST'
CHARACTER (LEN=3) :: my_component

! This routine may potentially be called in various places to cater for
! optional components. We don't want to read the input file multiple times
! etc so to ensure the code here is only performed once, we check to see if
! things have already been set up. 

IF (mype == 0) THEN
  OPEN(lb,FILE="translist",STATUS="OLD",ACTION='READ',IOSTAT=iostatus)
END IF

! Obtain the maximum index used by our coupling transients
CALL read_nml_transcount(lb,iostatus) 

! Set the name of this component
IF (l_junior) THEN
  my_component=jnr_component

  ! For the time being, Junior does no ocean coupling so we reset 
  ! ocean array sizes to zero in the junior component to avoid 
  ! unnecessarily allocating space and wasting time on irrelevant 
  ! loop and array processing. 
  transient_a2o_count=0
  transient_o2a_count=0
ELSE
  my_component=atm_component
END IF

ALLOCATE (transient_a2c(transient_a2c_count))
ALLOCATE (transient_c2a(transient_c2a_count))

ALLOCATE (transient_a2o(transient_a2o_count))
ALLOCATE (transient_o2a(transient_o2a_count))

! Initialise all our transient details with
! default conditions.
transient_a2c(:)%indx = -999
transient_c2a(:)%indx = -999
transient_a2o(:)%indx = -999
transient_o2a(:)%indx = -999

transient_a2c(:)%grid = ""
transient_c2a(:)%grid = ""
transient_a2o(:)%grid = ""
transient_o2a(:)%grid = ""

transient_a2c(:)%field_id = -999
transient_c2a(:)%field_id = -999
transient_a2o(:)%field_id = -999
transient_o2a(:)%field_id = -999

transient_a2c(:)%order = -999
transient_c2a(:)%order = -999
transient_a2o(:)%order = -999
transient_o2a(:)%order = -999

transient_a2c(:)%n_slices = 0
transient_c2a(:)%n_slices = 0
transient_a2o(:)%n_slices = 0
transient_o2a(:)%n_slices = 0

! Set default processing conditions:
transient_a2c(:)%ice_min = .FALSE.
transient_c2a(:)%ice_min = .FALSE.
transient_a2o(:)%ice_min = .FALSE.
transient_o2a(:)%ice_min = .FALSE.

indx_a2o = 0
indx_a2c = 0
indx_o2a = 0
indx_c2a = 0

! Read all the transient fields from our namelist and save the
! types and names in the appropriate in or out arrays.
DO WHILE ( iostatus == 0 )

  CALL read_nml_transfld(lb,iostatus)

  IF (iostatus == 0) THEN

    IF (my_component == atm_component) THEN

      ! Process ocean fields involved with atmosphere
      IF ((c_in(1:3) == atm_component).AND.   & 
          (c_out(1:3) == ocn_component)) THEN

        ! Compile list of arrays coming from the ocean. 
        indx_o2a = indx_o2a + 1

        ALLOCATE(transient_o2a(indx_o2a)%name(1:max_lev))
        ALLOCATE(transient_o2a(indx_o2a)%oasis_id(1:max_lev))

        ! Source component is ocean, target component atmos
        transient_o2a(indx_o2a)%grid     = c_in(4:4)
        transient_o2a(indx_o2a)%name(:)  = n_in
        transient_o2a(indx_o2a)%indx     = i_number
        transient_o2a(indx_o2a)%field_id = f_id
        transient_o2a(indx_o2a)%order    = map_type
        transient_o2a(indx_o2a)%c_from   = c_out(1:3)
        transient_o2a(indx_o2a)%c_to     = c_in(1:3)
        transient_o2a(indx_o2a)%n_slices = max_lev

        IF (max_lev > 1) THEN
          ! Set up full names for any multiple (3D) fields 
          DO k = 1, max_lev
            WRITE(transient_o2a(indx_o2a)%name(k)(6:8),'(I3.3)') k
          END DO
        END IF

      ELSE IF ((c_out(1:3) == atm_component).AND.  &
               (c_in(1:3) == ocn_component)) THEN

        ! Compile list of arrays going to the ocean. 
        indx_a2o = indx_a2o + 1 

        ALLOCATE(transient_a2o(indx_a2o)%name(1:max_lev))
        ALLOCATE(transient_a2o(indx_a2o)%oasis_id(1:max_lev))

        ! Source component is atmos, target component ocean
        transient_a2o(indx_a2o)%grid     = c_out(4:4)
        transient_a2o(indx_a2o)%name(:)  = n_out
        transient_a2o(indx_a2o)%indx     = i_number
        transient_a2o(indx_a2o)%field_id = f_id
        transient_a2o(indx_a2o)%order    = map_type
        transient_a2o(indx_a2o)%c_from   = c_out(1:3)
        transient_a2o(indx_a2o)%c_to     = c_in(1:3)
        transient_a2o(indx_a2o)%n_slices = max_lev

        IF (max_lev > 1) THEN
          ! Set up full names for any multiple (3D) fields 
          DO k = 1, max_lev
            WRITE(transient_a2o(indx_a2o)%name(k)(6:8),'(I3.3)') k
          END DO
        END IF
  
      END IF

    END IF ! Atmosphere processing only

    ! Process lists of arrays involved in chemistry coupling
    ! This code is common to both senior (UM atmos) and junior (chem)  

    IF ((c_out(1:3) == jnr_component).AND.  &
        (c_in(1:3) == atm_component)) THEN

      ! Source component is junior (chem), target component atmos  
      indx_c2a = indx_c2a + 1

      ALLOCATE(transient_c2a(indx_c2a)%name(1:max_lev))
      ALLOCATE(transient_c2a(indx_c2a)%oasis_id(1:max_lev))

      IF (l_junior) THEN
        ! This is an output field for the junior component
        transient_c2a(indx_c2a)%grid     = c_out(4:4)
        transient_c2a(indx_c2a)%name(:)  = n_out
      ELSE
        ! This is an input field for the senior component
        transient_c2a(indx_c2a)%grid     = c_in(4:4)
        transient_c2a(indx_c2a)%name(:)  = n_in
      END IF

      transient_c2a(indx_c2a)%indx     = i_number
      transient_c2a(indx_c2a)%field_id = f_id
      transient_c2a(indx_c2a)%order    = map_type
      transient_c2a(indx_c2a)%c_from   = c_out(1:3)
      transient_c2a(indx_c2a)%c_to     = c_in(1:3)
      transient_c2a(indx_c2a)%n_slices = max_lev

      IF (max_lev > 1) THEN
        ! Set up full names for multiple (3D) fields 
        DO k = 1, max_lev
          WRITE(transient_c2a(indx_c2a)%name(k)(6:8),'(I3.3)') k
        END DO
      END IF

    ELSE IF ((c_out(1:3) == atm_component).AND.  &
             (c_in(1:3) == jnr_component)) THEN

      ! Source component is atmos, target is junior (chem) 
      indx_a2c = indx_a2c + 1

      ALLOCATE(transient_a2c(indx_a2c)%name(1:max_lev))
      ALLOCATE(transient_a2c(indx_a2c)%oasis_id(1:max_lev))

      ! Source component is atmos, target component junior (chem) 
      IF (l_junior) THEN
        ! This is an input field for the junior component
        transient_a2c(indx_a2c)%grid     = c_in(4:4)
        transient_a2c(indx_a2c)%name(:)  = n_in
      ELSE
        ! This is an output field for the senior component
        transient_a2c(indx_a2c)%grid     = c_out(4:4)
        transient_a2c(indx_a2c)%name(:)  = n_out
      END IF

      transient_a2c(indx_a2c)%indx     = i_number
      transient_a2c(indx_a2c)%field_id = f_id
      transient_a2c(indx_a2c)%order    = map_type
      transient_a2c(indx_a2c)%c_from   = c_out(1:3)
      transient_a2c(indx_a2c)%c_to     = c_in(1:3)
      transient_a2c(indx_a2c)%n_slices = max_lev

      IF (max_lev > 1) THEN
        ! Set up full names for multiple (3D) fields 
        DO k = 1, max_lev
          WRITE(transient_a2c(indx_a2c)%name(k)(6:8),'(I3.3)') k
        END DO
      END IF

    END IF

  END IF ! iostatus OK.

END DO

IF (mype == 0) THEN
  CLOSE(lb)
END IF

END SUBROUTINE oasis_read_translist
