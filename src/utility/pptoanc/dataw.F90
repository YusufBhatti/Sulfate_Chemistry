! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine interface:
SUBROUTINE dataw(rows,columns,fieldsize,nlevels,levn,len_extra,   &
 fieldn,len1_lookup_all,lookup_all,fixhd,                         &
 len_cfi, cfi1, cfi2, cfi3, fldsizelev,ftin1,ftout,               &
 tracer_grid,add_wrap_pts,ibm_to_cray,compress,rmdi_input,wave,   &
 lsmask, l_bit_32,                                                &
 icode)

USE um_types, ONLY: real32
USE mask_compression, ONLY: compress_to_mask
USE Decomp_DB, ONLY: decompose
USE near_equal_real_mod, ONLY: near_equal_real
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: rmdi 
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE


!
! Description: This writes the data out using WRITFLD.
!               If compress oa_pack is used.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

! Subroutine arguments
!   Scalar arguments with intent(in):

INTEGER :: rows       ! number of rows in input pp field
INTEGER :: columns    ! number of columns in input pp field
INTEGER :: fieldsize  ! number of points in output anc. field
INTEGER :: nlevels    ! number of levels
INTEGER :: levn       ! current level number
INTEGER :: len_extra
INTEGER :: ftin1
INTEGER :: ftout
INTEGER :: fieldn
INTEGER :: len1_lookup_all
INTEGER :: icode          ! error status

REAL :: rmdi_input

LOGICAL :: tracer_grid
LOGICAL :: add_wrap_pts
LOGICAL :: ibm_to_cray
LOGICAL :: compress
LOGICAL :: l_bit_32

!   Array  arguments with intent(in):

INTEGER :: lookup_all(len1_lookup_all,*)
INTEGER :: fixhd(*)

INTEGER :: len_cfi(3)         ! dimensions of arrays
INTEGER :: cfi1(len_cfi(1))   !   compressed
INTEGER :: cfi2(len_cfi(2))   !   field index
INTEGER :: cfi3(len_cfi(3))   !   arrays
INTEGER :: fldsizelev(nlevels)   ! size of output field on each level

LOGICAL :: lsmask(rows*columns)

! local arrays
REAL(KIND=real32) :: datain(rows*columns)
REAL :: data_field(rows*columns)
REAL :: field_wrap(columns+2,rows)
REAL :: field_to_write(fieldsize)
REAL :: extra_data(len_extra+1)  ! space for extra data

! Local Scalars
INTEGER :: i
INTEGER :: j,istart,iend,ii
INTEGER :: field_type          ! 0 for tracers; 1 for velocities
INTEGER :: no_cmp   ! # of pts in full compressed field (all levels)
INTEGER :: no_rows_m   ! number of rows east-west on model grid
INTEGER :: n_sea_points

LOGICAL :: LTimer      ! timer switch (set to false)
LOGICAL :: cyclic_grid ! T => input field to OA_PACK has
                    !      overlap points
LOGICAL :: wave ! creating wave dump

LOGICAL :: l_skip ! If T, the data is read, but nothing is passed back.

CHARACTER(LEN=errormessagelength) :: cmessage      ! error message


! Input arguments for decompose_smexe
INTEGER ::                                                        &

  tot_levels        ! IN  :total number of levels


!- End of header

l_skip = .FALSE.
!   1. Read data and do number format conversion if needed
! DEPENDS ON: readdata
CALL readdata( rows, columns, ftin1, ibm_to_cray, len_extra,      &
               l_bit_32, l_skip, data_field, extra_data )


!  1.1 Convert real missing data indicators
IF ( rmdi_input  /=  rmdi) THEN
  i=0
  DO j = 1,rows*columns
    IF ( near_equal_real(data_field(j), rmdi_input) ) THEN
      !         if ( rmdi_input  >   0.0 ) then
      data_field(j) = rmdi
      i=i+1
    END IF
  END DO
  IF (i >  0) THEN
    WRITE(umMessage,*) i,' RMDI converted from ',rmdi_input,' to ',rmdi
    CALL umPrint(umMessage,src='dataw')
  END IF
END IF

!  2. Add in wrap points when add_wrap_pts=t

IF (add_wrap_pts) THEN

  DO j=1,rows
    DO i=1,columns
      field_wrap(i,j)=data_field(i+(j-1)*columns)
    END DO
  END DO
  DO j=1,rows
    field_wrap(columns+1,j)=field_wrap(1,j)
    field_wrap(columns+2,j)=field_wrap(2,j)
  END DO

  !  3. Pack data using compression indices when compress=t

  IF (compress) THEN

    IF (.NOT. wave) THEN

      IF (tracer_grid) THEN
        field_type = 0
        no_rows_m = rows
      ELSE
        field_type = 1
        no_rows_m = rows + 1
      END IF

      no_cmp = 0
      DO i = 1, nlevels   ! do not use levn in this loop
        no_cmp = no_cmp + fldsizelev(i)
      END DO

      cyclic_grid = .TRUE.   ! input pp fields do not
                              ! have wrap-points

      LTimer = .FALSE.
      icode  = 0

      ! DEPENDS ON: oa_pack
      CALL oa_pack(icode, cmessage, LTimer,                        &
        no_rows_m, columns+2, nlevels, len_cfi(1), fieldsize,      &
        cfi1, cfi2, cfi3, no_cmp, rmdi,                            &
        levn, field_type, cyclic_grid, field_wrap,                 &
        field_to_write)


      IF (icode  >   0) THEN
        WRITE(umMessage,*) 'error from OA_PACK:', cmessage
        CALL umPrint(umMessage,src='dataw')
        GO TO 9999
      END IF

    ELSE       ! add_wrap .and. compress .and. wave

      ! compress using SLMASK for wave model - use SEA POINTS set to TRUE
      ! a value for n-SEA-points is returned from this subroutine

      !!!!!!!!! This needs attention

      CALL compress_to_mask(data_field,field_to_write,lsmask,       &
          rows*columns,n_SEA_points)

      WRITE(umMessage,*)'after to land points no_cmp is ',n_sea_points
      CALL umPrint(umMessage,src='dataw')
      no_cmp=n_sea_points

    END IF

  ELSE      ! add_wrap .and. .not. compress

    DO j=1,rows
      DO i=1,columns+2
        field_to_write(i+(j-1)*(columns+2) ) =  field_wrap(i,j)
      END DO
    END DO

  END IF

ELSE         ! .not. add_wrap

  !  3.1 Pack data using compression indices when compress=t

  IF (compress) THEN

    IF (.NOT. wave) THEN

      IF (tracer_grid) THEN
        field_type = 0
        no_rows_m = rows
      ELSE
        field_type = 1
        no_rows_m = rows + 1
      END IF

      no_cmp = 0
      DO i = 1, nlevels   ! do not use levn in this loop
        no_cmp = no_cmp + fldsizelev(i)
      END DO

      cyclic_grid = .FALSE.   ! input pp fields do not
                              ! have wrap-points

      LTimer = .FALSE.
      icode  = 0

      ! DEPENDS ON: oa_pack
      CALL oa_pack(icode, cmessage, LTimer,                        &
        no_rows_m, columns, nlevels, len_cfi(1), fieldsize,        &
        cfi1, cfi2, cfi3, no_cmp, rmdi,                            &
        levn, field_type, cyclic_grid, data_field, field_to_write)


      IF (icode  >   0) THEN
        WRITE(umMessage,*) 'error from OA_PACK:', cmessage
        CALL umPrint(umMessage,src='dataw')
        GO TO 9999
      END IF

    ELSE         ! .not. add_wrap .and. compress .and. wave

      ! compress using SLMASK for wave model - use SEA POINTS set to TRUE
      ! a value for n-SEA-points is returned from this subroutine

      CALL compress_to_mask(data_field,field_to_write,lsmask,       &
          rows*columns,n_SEA_points)

      WRITE(umMessage,*)'after to land points no_cmp is ',n_sea_points
      CALL umPrint(umMessage,src='dataw')
      no_cmp=n_sea_points

    END IF

  ELSE         ! .not. add_wrap .and. .not. compress

    DO j = 1, fieldsize
      field_to_write(j) =  data_field(j)
    END DO

  END IF

END IF

!  5. Output data using WRITFLDS

!C TEMP print out data LSMASK for wave dump

IF (lookup_all(23,fieldn) == 38) THEN
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='dataw')
  WRITE(umMessage,*)'before writing data array'
  CALL umPrint(umMessage,src='dataw')
  istart=1
  iend=istart+columns-1
  DO i=rows,1,-1
    WRITE(umMessage,*) (field_to_write(ii),ii=istart,iend)
    CALL umPrint(umMessage,src='dataw')
    istart=istart+columns
    iend=iend+columns
  END DO
END IF

CALL decompose(columns, rows,0,0,1)

! DEPENDS ON: writflds
CALL writflds(ftout,1,fieldn,lookup_all,                          &
              len1_lookup_all,field_to_write,fieldsize,           &
              fixhd,                                              &
              icode,cmessage )

IF (icode  >   0) THEN
  WRITE(umMessage,*) 'error from WRITFLDS:', cmessage
  CALL umPrint(umMessage,src='dataw')
  GO TO 9999
END IF

9999  CONTINUE
RETURN
END SUBROUTINE dataw
! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:

!  Skip namelists in f90 compiled UM code removing need for
!  assign -f 77 g:sf  in script
!
! Subroutine Interface:
