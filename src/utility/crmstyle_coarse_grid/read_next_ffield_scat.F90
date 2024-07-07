! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Read requested fields from a UM ancillary fieldsfile, regrid and scatter

MODULE read_next_ffield_scat_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READ_NEXT_FFIELD_SCAT_MOD'

CONTAINS

SUBROUTINE read_next_ffield_scat (levels, istart,iend,iproc_typ,istashcode, &
                                  ff_hdr,                                   &  
                                  field_out,                                &  
                     ErrorStatus) 


! UM constants
USE planet_constants_mod, ONLY:                                      &
  kappa, pref, planet_radius
USE conversions_mod,     ONLY:                                       &
  pi_over_180

USE missing_data_mod, ONLY: rmdi

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words

USE crmwork_arrays_mod, ONLY:                                            &
  xcoslat_full, th_km_level_full, th_weight_full, prec_full
USE crmstyle_pp_data_mod, ONLY:                                          &
  bdx,bdy
USE hires_data_mod , ONLY:                                               &
  dpdx, dpdy

USE crmstyle_grid_info_mod, ONLY:                                 &
  nprocs,local_row_len, local_rows

USE crmstyle_cntl_mod, ONLY:                                      &
  in_cols, in_rows, l_whole_grid, l_all_sea, l_ENDGame, l_class_col   

! UM parallel info and grids
USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p
USE UM_ParCore, ONLY: mype

USE io

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength


! subroutines
USE put_on_fixed_heights_mod, ONLY: put_on_fixed_heights

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE packing_codes_mod, ONLY: PC_BitMask_CompressType

IMPLICIT NONE

! Description:
!   This subroutine reads a given number of fields, matching given
!   criteria from an open UM fieldsfile.  The fields are returned in field_out
!   Fields are read on PE 0, regridded if required and scattered
!   back to all the processors.
!
! Method:
!   See online documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   This code is written to UM programming standards version 8.3

! Subroutine arguments:

INTEGER, INTENT(IN) ::  &
 levels                 & ! Number of model levels required
,istart                 & ! first field required
,iend                   & ! last field required
,iproc_typ              & ! How to process the input field
,istashcode               ! stashcode of field

TYPE(UM_Header_type), INTENT(IN) :: ff_hdr

REAL(wp), INTENT(OUT)  ::                    &
  field_out(local_row_len,local_rows,levels)  ! Field for just region required
                                              ! by this PE on levels required

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "READ_NEXT_FFIELD_SCAT"

! Local Variables:
INTEGER :: i,j,k                      ! local counter

INTEGER ::           &
  NumCols            &
 ,Numrows            &
 ,NumStore           &
 ,icode                ! error code

LOGICAL :: LDecode = .TRUE.
LOGICAL :: PPHdrMod = .FALSE.
LOGICAL ::   &
  l_inter

REAL, POINTER :: TempArray(:,:)

REAL ::                                      &
  field_out_full(in_cols,in_rows,levels)     & ! work array
 ,field_out_local(local_row_len,local_rows)    ! work array

REAL(wp), ALLOCATABLE ::  &
  work1(:,:,:)            & ! work array
 ,work2(:,:,:)              ! work array

REAL, ALLOCATABLE ::                                      &
  dpdx_full(:,:)          & ! work array
 ,dpdy_full(:,:)            ! work array

REAL ::                                 &
  dy                                      ! used in dp/dy cal

TYPE(PP_Field_type) :: fields(levels)

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
IF (l_all_sea) THEN
  l_inter = .FALSE.
ELSE
  l_inter = .TRUE.
END IF



IF ( .NOT. ASSOCIATED(ff_hdr % Lookup) ) THEN
  WRITE(umMessage,'(A)') "Cannot get field since file uninitialised - skipping"
  CALL umPrint(umMessage,src=RoutineName)
  errorstatus = statuswarning
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Stops I/O from one PE being broadcast to all PEs (at present it is setup to
! do this for all the files).

CALL set_unit_bcast_flag(ff_hdr % UnitNum)
IF (mype == 0) THEN

  !-----------------------------------------------------------------------
  ! Loop through required fields
  !-----------------------------------------------------------------------

  DO k = 1,levels

    ! Clear space for field
    IF ( ASSOCIATED( Fields(k) % RData ) ) THEN
      DEALLOCATE( Fields(k) % RData )
      NULLIFY( Fields(k) % RData )
    END IF

    !----------------------------------------
    !  Search lookup for the required FIELD
    !----------------------------------------
    ! Search on LBTYP, LBLEV, LBFT - similar check in copyflds

    Fields(k) % LookupPos = istart + k -1

    Fields(k) % Hdr = ff_hdr % Lookup( Fields(k) % LookupPos )

    ! Allocate space for unpacked field
    ! Need to check if land packed.
    IF ( MOD(Fields(k) % Hdr % LBPack / 10, 10) == PC_BitMask_CompressType )   &
    THEN
      IF (ASSOCIATED(LSMField % RData)) THEN
        Fields(k) % Hdr % NumCols = LSMField % Hdr % NumCols
        Fields(k) % Hdr % Numrows = LSMField % Hdr % Numrows
      ELSE
        ErrorStatus = StatusWarning

        CALL EReport( RoutineName, ErrorStatus,                             &
                  "Cannot store land/sea packed field without LSM - skipping.")
        ! Reset Errorstatus but return warning back to calling routine at end.
        ErrorStatus = StatusOk
      END IF
    END IF

    NumCols = Fields(k) % Hdr % NumCols
    Numrows = Fields(k) % Hdr % Numrows

    ALLOCATE ( Fields(k) % RData( NumCols*Numrows,1 ) )

    !----------------------------------------
    ! Retrieving data
    !----------------------------------------

    ! DEPENDS ON: readfld
    CALL ReadFld( ff_hdr, LDecode, PPHdrMod, Fields(k), ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      WRITE(umMessage,'(A,I10)') ' Error in read field ',errorStatus
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(I10)') ff_hdr%lookup(k)
      CALL umPrint(umMessage,src=RoutineName)
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                              zhook_handle)
      RETURN
    END IF

    ! Now we need to reshape the read in data to be the correct shape of the
    ! unpacked field.
    TempArray => Fields(k) % RData
    NULLIFY(Fields(k) % RData)
    ALLOCATE( Fields(k) % RData( NumCols,Numrows ) )
    Fields(k) % RData = RESHAPE(source = TempArray,                   &
                                 SHAPE  = (/NumCols,Numrows/))
    DEALLOCATE(TempArray)
    NULLIFY(TempArray)

    ! store full surface precipitation on PE0
    IF (l_class_col .AND. (istashcode == 4203 .OR. istashcode == 4204) ) THEN
      ! Add precip fields 
      DO j=1,in_rows
        DO i=1,in_cols
          prec_full(i,j) = prec_full(i,j) + Fields(k) % RData(i,j)
        END DO
      END DO
    END IF
    !-----------------------------------------------------------------------
    ! Nothing in the headers indicates whether grid ENDGame or not.
    !-----------------------------------------------------------------------
    ! Decide how to process input field both U and V for LAM have same rows
    !-----------------------------------------------------------------------
    !   New dynamics LAM grid                           ENDGame  LAM grid
    !
    !   1   1   2   2   in_cols                     1   1   2   2    in_cols
    !   X   U   X   U    X   U   1                      V       V       V
    !
    !   V       V        V       1                  U   T   U   T   U   T
    !
    !   X   U   X   U    X   U   2                      V       V       V
    !
    !   V       V        V       2                  U   T   U   T   U   T
    !
    !   X   U   X   U    X   U   3                      V       V       V
    !
    !   V       V        V       3                  U   T   U   T   U   T
    !
    !   X   U   X   U    X   U in_rows                  V       V       V
    !
    !   V       V        V     in_rows              U   T   U   T   U   T
    !-----------------------------------------------------------------------

    SELECT CASE (iproc_typ)
    CASE (0)  ! pgrid
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
      DO j = 1,in_rows
        DO i = 1,in_cols
          field_out_full(i,j,k) = Fields(k) % RData(i,j)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (1)    ! U wind
      ! Put on p required grid
      IF (l_ENDGame) THEN   ! Input ENDGame grid
        IF (l_whole_grid) THEN ! assuming bicyclic
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                          &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
          DO j=1,in_rows
            field_out_full(in_cols,j,k) = 0.5* (Fields(k) % RData(1,j) +    &
                                        Fields(k) % RData(in_cols,j))
            DO i=1,in_cols-1
              field_out_full(i,j,k) =0.5* (Fields(k) % RData(i,j) +    &
                                               Fields(k) % RData(i+1,j))
            END DO
          END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
          DO j=1,in_rows
            field_out_full(in_cols,j,k) = 0.0   ! unset but will not be wanted
            DO i=1,in_cols-1
              field_out_full(i,j,k) =0.5* (Fields(k) % RData(i,j) +    &
                                               Fields(k) % RData(i+1,j))
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF

      ELSE                  ! Input New dynamics grid

        IF (l_whole_grid) THEN ! assuming bicyclic
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                          &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
          DO j=1,in_rows
            field_out_full(1,j,k) = 0.5* (Fields(k) % RData(1,j) +       &
                                        Fields(k) % RData(in_cols,j))
            DO i=2,in_cols
              field_out_full(i,j,k) =0.5* (Fields(k) % RData(i-1,j) +    &
                                               Fields(k) % RData(i,j))
            END DO
          END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
          DO j=1,in_rows
            field_out_full(1,j,k) = 0.0   ! unset but will not be wanted
            DO i=2,in_cols
              field_out_full(i,j,k) =0.5* (Fields(k) % RData(i-1,j) +    &
                                               Fields(k) % RData(i,j))
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
      END IF       ! test input grid

    CASE (2)    ! V wind
      IF (l_ENDGame) THEN   ! Input ENDGame grid
        IF (l_whole_grid) THEN ! assuming bicyclic
               ! Put on p required grid
          DO i=1,in_cols
            field_out_full(i,1,k) = 0.5* (Fields(k)%RData(i,1) +    &
                                               Fields(k)%RData(i,2))
            field_out_full(i,in_rows,k)= 0.5*(Fields(k)%RData(i,in_rows) + &
                                               Fields(k)%RData(i,1))
          END DO

        ELSE
          DO i=1,in_cols
            field_out_full(i,1,k) = 0.0         ! unset but will not be wanted
            field_out_full(i,in_rows,k) = 0.0   ! unset but will not be wanted
          END DO
        END IF
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
        DO j=2,in_rows-1
          DO i=1,in_cols
            field_out_full(i,j,k) =0.5* (Fields(k) % RData(i,j) +    &
                                               Fields(k) % RData(i,j+1))
          END DO
        END DO
!$OMP END PARALLEL DO

      ELSE                  ! Input New dynamics grid
        IF (l_whole_grid) THEN ! assuming bicyclic
               ! Put on p required grid
          DO i=1,in_cols
            field_out_full(i,1,k) = 0.5* (Fields(k)%RData(i,in_rows) +    &
                                               Fields(k)%RData(i,1))
            field_out_full(i,in_rows,k)= 0.5*(Fields(k)%RData(i,in_rows-1) + &
                                               Fields(k)%RData(i,in_rows))
          END DO

        ELSE
          DO i=1,in_cols
            field_out_full(i,1,k) = 0.0         ! unset but will not be wanted
            field_out_full(i,in_rows,k) = 0.0   ! unset but will not be wanted
          END DO
        END IF
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields)
        DO j=2,in_rows-1
          DO i=1,in_cols
            field_out_full(i,j,k) =0.5* (Fields(k) % RData(i,j-1) +    &
                                               Fields(k) % RData(i,j))
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

    CASE (3)    ! pressure - gradients required as well
               ! Convert to exner on full grid

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(k, in_rows, in_cols, field_out_full, Fields, pref, kappa)
      DO j = 1,in_rows
        DO i = 1,in_cols
          field_out_full(i,j,k) = (Fields(k) % RData(i,j) /pref)**kappa
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE DEFAULT

      WRITE(umMessage,'(A)') ' Unallowed option processing code '
      CALL umPrint(umMessage,src=RoutineName)

    END SELECT ! test on stashcode

    ! Release space
    DEALLOCATE( Fields(k) % RData )

  END DO   ! Loop over levels required


  IF (iproc_typ == 3) THEN
    ALLOCATE (work2(in_cols,in_rows,levels))
    ALLOCATE (work1(in_cols,in_rows,levels))

    ! At this point have a full 3d exner field
    ! Need 32 bit input field for put_on_fixed_heights
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, levels, work1, field_out_full)
    DO k=1,levels
      DO j = 1,in_rows
        DO i = 1,in_cols
          work1(i,j,k) = field_out_full(i,j,k)
        END DO
      END DO
    END DO    ! k loop
!$OMP END PARALLEL DO
      ! Put on fixed heights - full grid work2 = exner on fixed heights
    CALL put_on_fixed_heights(in_cols,in_rows,levels,th_km_level_full,     &
                              work1, th_weight_full,l_inter, work2)

    DEALLOCATE(work1)

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                            &
!$OMP& SHARED(levels, in_rows, in_cols, th_km_level_full, field_out_full, &
!$OMP&         work2, pref, kappa)
    DO k=1,levels
      ! Convert from exner back to P - full grid
      DO j=1,in_rows
        DO i=1,in_cols
          IF (th_km_level_full(i,j,k) >= 0) THEN  ! above surface
            field_out_full(i,j,k) = pref*(work2(i,j,k)**(1.0/kappa))
          ELSE
            field_out_full(i,j,k) = rmdi
          END IF
        END DO
      END DO

    END DO     ! k loop
!$OMP END PARALLEL DO

    DEALLOCATE(work2)
  END IF

END IF  ! mype == 0

! Allow all processors/nodes to do IO again
CALL clear_unit_bcast_flag(ff_hdr % UnitNum)

icode = 0
! Force synchronisation before trying to scatter back fields
CALL  gc_gsync(nprocs,icode)

! scatter field back into output array
DO k=1,levels

  ! DEPENDS ON: scatter_field
  CALL scatter_field( field_out_local, field_out_full(1,1,k),         &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                   0,gc_all_proc_group)

  ! copy to 32 bit field
  DO j = 1,local_rows
    DO i = 1,local_row_len
      field_out(i,j,k) = field_out_local(i,j)
    END DO
  END DO

END DO   ! k

! scatter back dp/dx & dpdy to processors
IF (iproc_typ == 3) THEN
  ! Need dp/dx and dp/dy

  dy=planet_radius*bdy*pi_over_180
  ALLOCATE(dpdx_full(in_cols,in_rows))
  ALLOCATE(dpdy_full(in_cols,in_rows))

  DO k=1,levels
    ! Calculate dp/dx  - first and last values on a row set to zero
    IF (mype == 0) THEN
      IF (l_whole_grid) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k,in_rows, in_cols, th_km_level_full, field_out_full,        &
!$OMP& dpdx_full, dy, xcoslat_full)
        DO j=1,in_rows
          ! Initialise to zero
          i=1
          IF ((th_km_level_full(i+1,j,k) >= 0) .AND.                       &
              (th_km_level_full(i,j,k) >= 0)  .AND.                        &
              (th_km_level_full(in_cols,j,k) >= 0) ) THEN
            dpdx_full(i,j) = 0.5*(field_out_full(i+1,j,k)                 &
                                      - field_out_full(in_cols,j,k))/     &
                                                (dy*xcoslat_full(i,j))
          END IF
          i=in_cols
          IF ((th_km_level_full(1,j,k) >= 0) .AND.                         &
              (th_km_level_full(i,j,k) >= 0)  .AND.                        &
              (th_km_level_full(in_cols,j,k) >= 0) ) THEN
            dpdx_full(i,j) = 0.5*(field_out_full(1,j,k)                   &
                                      - field_out_full(in_cols,j,k))/     &
                                                (dy*xcoslat_full(i,j))
          END IF

          DO i=2,in_cols-1
            IF ((th_km_level_full(i+1,j,k) >= 0) .AND.                       &
                (th_km_level_full(i,j,k) >= 0)  .AND.                        &
                (th_km_level_full(i-1,j,k) >= 0) ) THEN
              dpdx_full(i,j) = 0.5*(field_out_full(i+1,j,k)                 &
                                        - field_out_full(i-1,j,k))/       &
                                                  (dy*xcoslat_full(i,j))
            END IF
          END DO
        END DO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k,in_rows, in_cols, th_km_level_full, field_out_full,        &
!$OMP& dpdx_full, dy, xcoslat_full)
        DO j=1,in_rows
          ! Initialise to zero incase below surface
          DO i=1,in_cols
            dpdx_full(i,j) = 0.0
          END DO

          DO i=2,in_cols-1
            IF ((th_km_level_full(i+1,j,k) >= 0) .AND.                       &
                (th_km_level_full(i,j,k) >= 0)  .AND.                        &
                (th_km_level_full(i-1,j,k) >= 0) ) THEN
              dpdx_full(i,j) = 0.5*(field_out_full(i+1,j,k)                 &
                                        - field_out_full(i-1,j,k))/       &
                                                  (dy*xcoslat_full(i,j))
            END IF
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF
    END IF
    icode = 0
    CALL  gc_gsync(nprocs,icode)


    ! DEPENDS ON: scatter_field
    CALL scatter_field( field_out_local, dpdx_full,                   &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                   0,gc_all_proc_group )

    ! copy to 32 bit field
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k, local_rows, local_row_len, dpdx, field_out_local)
    DO j = 1,local_rows
      DO i = 1,local_row_len
        dpdx(i,j,k) = field_out_local(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

        ! Calculate  dp/dy - first and last local_rows set to zero

    IF (mype == 0) THEN

      IF (l_whole_grid) THEN
        j=1
        DO i=1,in_cols
          IF ((th_km_level_full(i,j+1,k) >= 0) .AND.                       &
              (th_km_level_full(i,j,k) >= 0)   .AND.                       &
              (th_km_level_full(i,in_rows,k) >= 0) ) THEN
            dpdy_full(i,j) = 0.5*( field_out_full(i,j+1,k)                 &
                                      - field_out_full(i,in_rows,k) )/dy
          END IF
        END DO
        j=in_rows
        DO i=1,in_cols
          IF ((th_km_level_full(i,1,k) >= 0) .AND.                         &
              (th_km_level_full(i,j,k) >= 0)   .AND.                       &
              (th_km_level_full(i,j-1,k) >= 0) ) THEN
            dpdy_full(i,j) = 0.5*( field_out_full(i,1,k)                   &
                                      - field_out_full(i,j-1,k) )/dy
          END IF
        END DO
      ELSE
        ! Initialise edges to zero
        DO i=1,in_cols
          dpdy_full(i,1) = 0.0
          dpdy_full(i,in_rows) = 0.0
        END DO
      END IF

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k,in_rows, in_cols, th_km_level_full, field_out_full,        &
!$OMP& dpdy_full, dy)
      DO j=2,in_rows-1
        DO i=1,in_cols
          IF ((th_km_level_full(i,j+1,k) >= 0) .AND.                       &
              (th_km_level_full(i,j,k) >= 0)   .AND.                       &
              (th_km_level_full(i,j-1,k) >= 0) ) THEN
            dpdy_full(i,j) = 0.5*( field_out_full(i,j+1,k)                 &
                                      - field_out_full(i,j-1,k) )/dy
          ELSE
            ! Cannot evaluate as at least some of points below surface.
            dpdy_full(i,j) = 0.0
          END IF
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF
    icode = 0
    CALL  gc_gsync(nprocs,icode)

    ! DEPENDS ON: scatter_field
    CALL scatter_field( field_out_local, dpdy_full,                   &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                   0,gc_all_proc_group )

    ! copy to 32 bit field

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k, local_rows, local_row_len, dpdy, field_out_local)
    DO j = 1,local_rows
      DO i = 1,local_row_len
        dpdy(i,j,k) = field_out_local(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

  END DO    ! k loop

  DEALLOCATE(dpdy_full)
  DEALLOCATE(dpdx_full)

END IF    ! iproc_typ

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-----------------------------------------------------------------------------

RETURN
END SUBROUTINE read_next_ffield_scat

END MODULE read_next_ffield_scat_mod
