! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Obtain UKCA-MODE input to UKCA_RADAER from D1.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
MODULE ukca_radaer_get_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_RADAER_GET_MOD'

CONTAINS

SUBROUTINE ukca_radaer_get(            &
                           ierr        &
  ,                        cmessage    &
  ,                        first_call, &
                           ukca_radaer &
  )

USE submodel_mod, ONLY: submodel_for_sm, atmos_im
USE d1_array_mod, ONLY: d1_object_type, d1_section, d1_item, d1_address,   &
                        d1_length, d1_no_levels, d1_halo_type, prognostic, &
                        d1, d1_addr, no_obj_d1

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE atm_fields_bounds_mod, ONLY: tdims_s, tdims

USE ukca_d1_defs, ONLY: ukca_sect
USE ukca_radaer_struct_mod, ONLY: ukca_radaer_struct

USE UM_ParVars, ONLY: halosize
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, model_levels, n_obj_d1_max, row_length, rows

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!
! Arguments with intent(in)
!
!
! Error indicator (0 is OK, >0 error)
!
INTEGER :: ierr      ! error code

! Error message if ierr is larger than 0
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER :: klev1     ! first model level

! Logical indicating first call: complete setup has to be done
LOGICAL :: first_call

!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer

!
! Local variables
!
INTEGER :: i, &
        j
INTEGER :: i_obj

!
! Atmosphere submodel index
!
INTEGER :: m_atm_modl

!
! STASH section where UKCA diagnostics are expected to reside.
!  Mass-mixing ratios/numbers and misc diagnostics (diameters,
!  densities, and volumes) should both be in section 34.

INTEGER, PARAMETER :: ukca_section_mmr = ukca_sect

!
! Logical for checking whether all diags have been found
!
LOGICAL :: l_diag_missing

!
! Local buffer read from D1 and its size.
!
REAL, ALLOCATABLE :: buffer(:,:,:)
INTEGER :: buffer_size

!
! Halo sizes
!
INTEGER :: halo_x
INTEGER :: prev_halo_x
INTEGER :: halo_y
INTEGER :: prev_halo_y

! Tracer Levels
INTEGER :: tr_levs

!
! Variables for looking for tagged diagnostics in D1.
!
INTEGER :: tag, ptd1, section, item, levs, dlen, addr, halotyp

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_RADAER_GET'

!
!
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

m_atm_modl = submodel_for_sm(atmos_im)

! Set tracer levels
tr_levs = tdims_s%k_len
!
! For first call, check that all diagnostics required are available
! at the expected dimensions. Retain their D1 address for later
! calls.
!

IF (first_call) THEN

  !
  ! Initialise all d1 addresses in mode and component
  ! information structures to "not found" (-1)
  !
  DO j = 1, ukca_radaer%n_mode
    ukca_radaer%d1_address_dry(j) = -1
    ukca_radaer%d1_address_wet(j) = -1
    ukca_radaer%d1_address_rho(j) = -1
  END DO ! j

  DO i = 1, ukca_radaer%n_cpnt
    ukca_radaer%d1_address_mmr(i) = -1
    ukca_radaer%d1_address_cvl(i) = -1
  END DO ! i

  !
  ! Go through D1 and retain some information about the
  ! prognostics required.
  DO i_obj=1,no_obj_D1(m_atm_modl)

    IF ( d1_addr(D1_object_type,i_obj,m_atm_modl) == prognostic ) THEN
      section = d1_addr(D1_section,i_obj,m_atm_modl)

      IF (section == ukca_section_mmr) THEN

        !
        ! Get addresses for the component mass-mixing
        ! ratios and modal number mixing ratios
        !

        item = d1_addr(D1_item,i_obj,m_atm_modl)
        levs = d1_addr(d1_no_levels,i_obj,m_atm_modl)
        dlen = d1_addr(d1_length,i_obj,m_atm_modl)
        addr    = d1_addr(d1_address,i_obj,m_atm_modl)
        halotyp = d1_addr(d1_halo_type,i_obj,m_atm_modl)

        DO i = 1, ukca_radaer%n_cpnt

          IF (item == ukca_radaer%stashcode_mmr(i)) THEN

            ukca_radaer%d1_address_mmr(i) = addr
            ukca_radaer%d1_nlevs_mmr(i)   = levs
            ukca_radaer%d1_length_mmr(i)  = dlen
            ! Mass-mixing ratios have halos
            ukca_radaer%d1_halo_type_mmr(i) = halotyp

          END IF

        END DO ! i (components)

        DO j = 1, ukca_radaer%n_mode

          IF (item == ukca_radaer%stashcode_nbr(j)) THEN

            ukca_radaer%d1_address_nbr(j) = addr
            ukca_radaer%d1_nlevs_nbr(j)   = levs
            ukca_radaer%d1_length_nbr(j)  = dlen
            ! Modal numbers have halos
            ukca_radaer%d1_halo_type_nbr(j) = halotyp

          END IF


        END DO

        ! Next find the addresses for dry and wet diameters,
        ! modal densities, water volume, and component volumes.

        DO j = 1, ukca_radaer%n_mode

          IF (item == ukca_radaer%stashcode_dry(j)) THEN

            ukca_radaer%d1_address_dry(j) = addr
            ukca_radaer%d1_nlevs_dry(j)   = levs
            ukca_radaer%d1_length_dry(j)  = dlen

          END IF

          IF (ukca_radaer%l_soluble(j)) THEN

            IF (item == ukca_radaer%stashcode_wet(j)) THEN

              ukca_radaer%d1_address_wet(j) = addr
              ukca_radaer%d1_nlevs_wet(j)   = levs
              ukca_radaer%d1_length_wet(j)  = dlen

            ELSE IF (item == ukca_radaer%stashcode_wtv(j)) THEN


              ukca_radaer%d1_address_wtv(j) = addr
              ukca_radaer%d1_nlevs_wtv(j)   = levs
              ukca_radaer%d1_length_wtv(j)  = dlen

            END IF

          END IF

          IF (item == ukca_radaer%stashcode_rho(j)) THEN

            ukca_radaer%d1_address_rho(j) = addr
            ukca_radaer%d1_nlevs_rho(j)   = levs
            ukca_radaer%d1_length_rho(j)  = dlen

          END IF

        END DO ! j (modes)

        DO i = 1, ukca_radaer%n_cpnt

          IF (item == ukca_radaer%stashcode_cvl(i)) THEN

            ukca_radaer%d1_address_cvl(i) = addr
            ukca_radaer%d1_nlevs_cvl(i)   = levs
            ukca_radaer%d1_length_cvl(i)  = dlen

          END IF

        END DO ! i (components)

      END IF  ! section = ukca_section_mmr

    END IF ! prognostic

  END DO ! i_obj (D1 items)

  !
  ! Check that all diagnostics have been found.
  !
  l_diag_missing = .FALSE.

  DO j = 1, ukca_radaer%n_mode

    IF (ukca_radaer%d1_address_dry(j) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                        &
            ukca_radaer%stashcode_dry(j), ' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    ELSE IF (ukca_radaer%l_soluble(j) .AND.                       &
        ukca_radaer%d1_address_wet(j) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                            &
             ukca_radaer%stashcode_wet(j), ' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    ELSE IF (ukca_radaer%l_soluble(j) .AND.                       &
        ukca_radaer%d1_address_wtv(j) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                          &
              ukca_radaer%stashcode_wtv(j),' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    ELSE IF (ukca_radaer%d1_address_rho(j) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                          &
              ukca_radaer%stashcode_rho(j),' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    ELSE IF (ukca_radaer%d1_address_nbr(j) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                          &
               ukca_radaer%stashcode_nbr(j),' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    END IF

  END DO ! j

  IF (l_diag_missing) THEN
    ierr = 702
    cmessage =                                                    &
      'ukca_radaer_get: Diag(s) needed for UKCA are missing from D1.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  DO i = 1, ukca_radaer%n_cpnt

    IF (ukca_radaer%d1_address_mmr(i) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                          &
           ukca_radaer%stashcode_mmr(i),' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    ELSE IF (ukca_radaer%d1_address_cvl(i) == -1) THEN
      WRITE(umMessage,'(A,I6,A)') 'Diagnostic ',                          &
           ukca_radaer%stashcode_cvl(i), ' not found in D1.'
      CALL umPrint(umMessage,src='ukca_radaer_get')
      l_diag_missing = .TRUE.

    END IF

  END DO ! i

  IF (l_diag_missing) THEN
    ierr = 702
    cmessage =                                                    &
      'ukca_radaer_get: Diag(s) needed for UKCA are missing from D1.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

END IF ! first_call

!
! At this point, we know where to look for the requested D1 items.
! Retrieve. The retrieval takes place into a temporary buffer that
! is reshaped.
!

buffer_size = row_length * rows * (model_levels+1)

klev1 = tdims%k_start
ALLOCATE(buffer(row_length, rows, klev1:model_levels))

DO j = 1, ukca_radaer%n_mode

  !
  ! Modal dry diameter.
  ! Check number of levels and total length.
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_dry(j) /= SIZE(buffer,dim=3)) THEN
    ierr = 705
    WRITE(umMessage,'(3(A,I5))')'Expecting ',SIZE(buffer,dim=3),        &
          ' levels, got ', ukca_radaer%d1_nlevs_dry(j),' for index: ',j
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage = 'ukca_radaer_get: Unexpected number of levels in D1'//   &
               ' diag. for dry diam.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                            zhook_handle)
    RETURN
  END IF

  IF (ukca_radaer%d1_length_dry(j) /= buffer_size) THEN
    ierr = 706
    WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,      &
             ' elements, got ', ukca_radaer%d1_length_dry(j)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, &
                            zhook_handle)
    RETURN
  END IF

  buffer = RESHAPE(d1(ukca_radaer%d1_address_dry(j) :             &
                      ukca_radaer%d1_address_dry(j) +             &
                      ukca_radaer%d1_length_dry(j) - 1),          &
                   (/row_length, rows, SIZE(buffer,dim=3)/))
  ukca_radaer%dry_diam(:,:,1:model_levels,j) = buffer(:,:,1:model_levels)

  IF (ukca_radaer%l_soluble(j)) THEN

    !
    ! Modal wet diameter.
    ! Check number of levels and total length.
    ! If ok, read into buffer and reshape into target array.
    !
    IF (ukca_radaer%d1_nlevs_wet(j) /= SIZE(buffer,dim=3)) THEN
      ierr = 705
      WRITE(umMessage,'(3(A,I5))') 'Expecting ', SIZE(buffer,dim=3),      &
            ' levels, got ', ukca_radaer%d1_nlevs_wet(j),' for index ',j
      CALL umPrint(umMessage,src='ukca_radaer_get')
      cmessage = 'ukca_radaer_get: Unexpected number of levels in D1'//   &
                 ' diag. for wet diam.'
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                              zhook_handle)
      RETURN
    END IF

    IF (ukca_radaer%d1_length_wet(j) /= buffer_size) THEN
      ierr = 706
      WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,            &
               ' elements, got ',ukca_radaer%d1_length_wet(j)
      CALL umPrint(umMessage,src='ukca_radaer_get')
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected total size of D1 diag.'
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                              zhook_handle)
      RETURN
    END IF

    buffer = RESHAPE(d1(ukca_radaer%d1_address_wet(j) :           &
                        ukca_radaer%d1_address_wet(j) +           &
                        ukca_radaer%d1_length_wet(j) - 1),        &
                   (/row_length, rows, SIZE(buffer,dim=3)/))
    ukca_radaer%wet_diam(:,:,1:model_levels,j) = buffer(:,:,1:model_levels)

    !
    ! Volume of water.
    ! Check number of levels and total length.
    ! If ok, read into buffer and reshape into target array.
    !
    IF (ukca_radaer%d1_nlevs_wtv(j) /= SIZE(buffer,dim=3)) THEN
      ierr = 705
      WRITE(umMessage,'(3(A,I5))') 'Expecting ', SIZE(buffer,dim=3),      &
            ' levels, got ',ukca_radaer%d1_nlevs_wtv(j),' for index: ',j
      CALL umPrint(umMessage,src='ukca_radaer_get')
      cmessage = 'ukca_radaer_get: Unexpected number of levels in D1'//   &
                 ' diag. for water volume'
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                              zhook_handle)
      RETURN
    END IF

    IF (ukca_radaer%d1_length_wtv(j) /= buffer_size) THEN
      ierr = 706
      WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,            &
               ' elements, got ',ukca_radaer%d1_length_wtv(j)
      CALL umPrint(umMessage,src='ukca_radaer_get')
      cmessage =                                                  &
        'ukca_radaer_get: Unexpected total size of D1 diag.'
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                              zhook_handle)
      RETURN
    END IF

    buffer = RESHAPE(d1(ukca_radaer%d1_address_wtv(j) :           &
                        ukca_radaer%d1_address_wtv(j) +           &
                        ukca_radaer%d1_length_wtv(j) - 1),        &
                   (/row_length, rows, SIZE(buffer,dim=3)/))
    ukca_radaer%modal_wtv(:,:,1:model_levels,j) = buffer(:,:,1:model_levels)

  ELSE

    !
    ! Insoluble modes: copy the dry diameter into the wet
    ! diameter, for consistency, and set the volume of water
    ! to zero.
    !
    ukca_radaer%wet_diam(:,:,:,j) = ukca_radaer%dry_diam(:,:,:,j)
    ukca_radaer%modal_wtv(:,:,:,j) = 0.0e+00

  END IF ! l_soluble

  !
  ! Modal density.
  ! Check number of levels and total length.
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_rho(j) /= SIZE(buffer,dim=3)) THEN
    ierr = 705
    WRITE(umMessage,'(3(A,I5))') 'Expecting ', SIZE(buffer,dim=3),        &
            ' levels, got ',ukca_radaer%d1_nlevs_rho(j),' for index: ',j
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                            &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  IF (ukca_radaer%d1_length_rho(j) /= buffer_size) THEN
    ierr = 706
    WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,            &
            ' elements, got ', ukca_radaer%d1_length_rho(j)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                            zhook_handle)
    RETURN
  END IF

  buffer = RESHAPE(d1(ukca_radaer%d1_address_rho(j) :             &
                      ukca_radaer%d1_address_rho(j) +             &
                      ukca_radaer%d1_length_rho(j) - 1),          &
                   (/row_length, rows, SIZE(buffer,dim=3)/))
  ukca_radaer%modal_rho(:,:,1:model_levels,j) = buffer(:,:,1:model_levels)

  !
  ! Initialise modal volume to that of water (the latter has been
  ! set to zero for insoluble modes above).
  ! The other components will be added in the loop below.
  !
  ukca_radaer%modal_vol(:,:,:,j) = ukca_radaer%modal_wtv(:,:,:,j)

END DO ! j

!
! Variables depending on components now.
!
DO i = 1, ukca_radaer%n_cpnt

  !
  ! Component fractional volumes.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_cvl(i) /= SIZE(buffer,dim=3)) THEN
    ierr = 705
    WRITE(umMessage,'(3(A,I5))') 'Expecting ', SIZE(buffer,dim=3),       &
           ' levels, got ', ukca_radaer%d1_nlevs_cvl(i),' for index: ',i
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage = 'ukca_radaer_get: Unexpected number of levels in D1'//    &
               ' diag. for fractional volume'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,     &
                            zhook_handle)
    RETURN
  END IF

  IF (ukca_radaer%d1_length_cvl(i) /= buffer_size) THEN
    ierr = 706
    WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,            &
         ' elements, got ', ukca_radaer%d1_length_cvl(i)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,    &
                            zhook_handle)
    RETURN
  END IF

  buffer = RESHAPE(d1(ukca_radaer%d1_address_cvl(i) :             &
                      ukca_radaer%d1_address_cvl(i) +             &
                      ukca_radaer%d1_length_cvl(i) - 1),          &
                   (/row_length, rows, SIZE(buffer,dim=3)/))
  ukca_radaer%comp_vol(:,:,1:model_levels,i) = buffer(:,:,1:model_levels)

END DO ! i

!
! Update the volume of each mode by adding the volume
! of each component within that mode.
!
DO j = 1, ukca_radaer%n_mode

  DO i = 1, ukca_radaer%n_cpnt_in_mode(j)

    ukca_radaer%modal_vol(:,:,:,j) = &
      ukca_radaer%modal_vol(:,:,:,j) + &
      ukca_radaer%comp_vol(:,:,:,ukca_radaer%i_cpnt_index(i, j))

  END DO ! i

END DO ! j

!
! Deallocate the buffer. It will be reallocated to account
! for halos in mass mixing ratios and modal number.
!
DEALLOCATE(buffer)
buffer_size = 0

prev_halo_x = -1
prev_halo_y = -1

DO i = 1, ukca_radaer%n_cpnt

  !
  ! Mass mixing ratio.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_mmr(i) /= tr_levs ) THEN
    ierr = 707
    WRITE(umMessage,'(2(A,I5))') 'Expecting ', tr_levs, ' levels, got ',  &
                ukca_radaer%d1_nlevs_mmr(i)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  halo_x = halosize(1, ukca_radaer%d1_halo_type_mmr(i))
  halo_y = halosize(2, ukca_radaer%d1_halo_type_mmr(i))

  IF (halo_x /= prev_halo_x .AND.                                 &
      halo_y /= prev_halo_y) THEN
    !
    ! Buffer needs to be allocated/reallocated.
    !
    IF (buffer_size /= 0) THEN
      DEALLOCATE(buffer)
    END IF

    buffer_size =  tdims_s%i_len  *                               &
                   tdims_s%j_len  *                               &
                  tr_levs
    ALLOCATE(buffer(tdims_s%i_start:tdims_s%i_end,                &
                    tdims_s%j_start:tdims_s%j_end,                &
                    tdims_s%k_start:tdims_s%k_end))

    prev_halo_x = halo_x
    prev_halo_y = halo_y

  END IF

  IF (ukca_radaer%d1_length_mmr(i) /= buffer_size) THEN
    ierr = 708
    WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,              &
            ' elements, got ',ukca_radaer%d1_length_mmr(i)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  buffer = RESHAPE(d1(ukca_radaer%d1_address_mmr(i) :             &
                      ukca_radaer%d1_address_mmr(i) +             &
                      ukca_radaer%d1_length_mmr(i) - 1),          &
                   (/  tdims_s%i_len ,                            &
                       tdims_s%j_len ,                            &
                      tr_levs /))

  ukca_radaer%mix_ratio(:, :, :, i) =                             &
                buffer(1:row_length, 1:rows, 1:model_levels)

END DO ! i

DO j = 1, ukca_radaer%n_mode

  !
  ! Modal number concentrations.
  ! Check number of levels, get halo size and check
  ! total length
  ! If ok, read into buffer and reshape into target array.
  !
  IF (ukca_radaer%d1_nlevs_nbr(j) /= tr_levs) THEN
    ierr = 707
    WRITE(umMessage,'(2(A,I5))') 'Expecting ', tr_levs, ' levels, got ',  &
                ukca_radaer%d1_nlevs_nbr(j)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected number of levels in D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  halo_x = halosize(1, ukca_radaer%d1_halo_type_nbr(j))
  halo_y = halosize(2, ukca_radaer%d1_halo_type_nbr(j))

  IF (halo_x /= prev_halo_x .AND.                                 &
      halo_y /= prev_halo_y) THEN
    !
    ! Buffer needs to be allocated/reallocated.
    !
    IF (buffer_size /= 0) THEN
      DEALLOCATE(buffer)
    END IF

    buffer_size =  tdims_s%i_len  *                               &
                   tdims_s%j_len  *                               &
                  tr_levs
    ALLOCATE(buffer(tdims_s%i_start:tdims_s%i_end,                &
                    tdims_s%j_start:tdims_s%j_end,                &
                    tdims_s%k_start:tdims_s%k_end))

    prev_halo_x = halo_x
    prev_halo_y = halo_y

  END IF

  IF (ukca_radaer%d1_length_nbr(j) /= buffer_size) THEN
    ierr = 708
    WRITE(umMessage,'(2(A,I10))') 'Expecting ', buffer_size,              &
         ' elements, got ', ukca_radaer%d1_length_nbr(j)
    CALL umPrint(umMessage,src='ukca_radaer_get')
    cmessage =                                                    &
      'ukca_radaer_get: Unexpected total size of D1 diag.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out,      &
                            zhook_handle)
    RETURN
  END IF

  buffer = RESHAPE(d1(ukca_radaer%d1_address_nbr(j) :             &
                      ukca_radaer%d1_address_nbr(j) +             &
                      ukca_radaer%d1_length_nbr(j) - 1),          &
                   (/  tdims_s%i_len ,                            &
                       tdims_s%j_len ,                            &
                      tr_levs /))

  ukca_radaer%modal_nbr(:, :, :, j) =                             &
                buffer(1:row_length, 1:rows, 1:model_levels)

END DO ! j

DEALLOCATE(buffer)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_get
END MODULE ukca_radaer_get_mod
