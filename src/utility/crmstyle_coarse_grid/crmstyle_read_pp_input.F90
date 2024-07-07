! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  read from pp input file

MODULE crmstyle_read_pp_input_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_READ_PP_INPUT_MOD'

CONTAINS

SUBROUTINE crmstyle_read_pp_input(pp_env,error_code)

USE hires_data_mod , ONLY:                                                   &
  landsea, precip, rain, snow, zh, lh, sh, pstar, tstar,                     &
  theta, thetav, t, rh, q, qcl, qcf,qrain, qgraup,                           &
  p_theta_lev,  u, v, w, dpdx, dpdy,   density,                              &
  dt1, dt2, dt4, dt9, dt12, dq4, dq9, dq12,                                  &
  dqcl4, dqcl9, dqcl12,  dqcf4, dqcf3, dqcf12

USE crmstyle_cntl_mod, ONLY:                                                 &
  mlevs, in_rows, in_cols, l_all_sea, l_ENDGame, l_class_col

USE crmstyle_grid_info_mod, ONLY:                                           &
  nprocs, local_row_len, local_rows

USE crmstyle_pp_data_mod, ONLY:                                              &
  iyear,imon,iday,ihour,imin,isec,isyear,ismon,isday,ishour,ismin,issec,     &
  bdx,bdy

USE crmwork_arrays_mod, ONLY:                                                &
  xcoslat_full, th_km_level_full, th_weight_full, prec_full,                 &
  h_theta_sea, mask, uv_km_level, th_km_level, uv_weight, th_weight

USE field_flags_mod, ONLY:                                                   &
   l_u, l_v, l_w, l_theta, l_ptheta, l_q, l_qcl, l_qcf, l_qrain, l_qgraup,   &
   l_dt1, l_dt2, l_dt4, l_dt9, l_dt12,  l_dq4, l_dq9, l_dq12,                &
   l_dqcl4, l_dqcl9, l_dqcl12,                                               &
   l_dqcf4, l_dqcf3, l_dqcf12,  l_got_fields,                                &
   l_rain, l_snow, l_precip, l_sh, l_lh, l_zh, l_tstar, l_pstar

USE planet_constants_mod, ONLY:                                              &
  kappa, cp, pref, r,  c_virtual, planet_radius
USE conversions_mod,     ONLY:                                               &
  pi_over_180

USE missing_data_mod, ONLY: rmdi

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words

USE filenamelength_mod,     ONLY: filenamelength
USE file_manager,           ONLY: assign_file_unit, release_file_unit
USE errormessagelength_mod, ONLY: errormessagelength

USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p
USE UM_ParCore, ONLY: mype

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE ereport_mod, ONLY: ereport, ereport_finalise

! subroutines
USE put_on_fixed_heights_mod, ONLY: put_on_fixed_heights
USE read_next_pp_field_mod,   ONLY: read_next_pp_field
USE reset_field_flags_mod,    ONLY: reset_field_flags
USE get_env_var_mod,          ONLY: get_env_var

USE qsat_mod, ONLY: qsat

USE umPrintMgr, ONLY: umprint, ummessage, newline


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Read in pp file of model output
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

CHARACTER(LEN=filenamelength), INTENT(IN)  :: &
  pp_env                                     ! pp file environment variable

INTEGER, INTENT(OUT) :: error_code   ! error code (non-zero if fails)

!-------------------------------------------------------------------------
! Local variables
CHARACTER(LEN=filenamelength) :: filename ! pp file filename

INTEGER ::               &
  i, j, k, l               ! loop counters

INTEGER ::               &
  unit_in                & ! unit number for pp_file
 ,errorstatus            & ! error code
 ,read_code              & ! read code
 ,stashcode              & ! stash code
 ,dmin                   & ! interval for single fields
 ,levels                 & ! Number of levels
 ,mesg                   & ! Identifier for GCOM message
 ,icode                  & ! Error code
 ,num_fields_found       & ! number of different fields found so far
 ,num_fields_wanted        ! number of different fields wanted

INTEGER ::               &
  iihead(45)             & ! integer pp header
 ,isend(15)                ! Information to send from PE 0 and received by all

LOGICAL ::               &
  lsend                    ! Sould pp field read in be scattered to other PEs
LOGICAL ::               &
  l_inter_th             & ! interpolate height fields theta input
 ,l_inter_uv               ! interpolate height fields uv input

REAL(wp) ::      &
  dy             & !
 ,exner            ! exner

REAL(wp)  ::                       &
  field(in_cols,in_rows,mlevs)     & ! Full input field but just required levels
 ,data_grid(local_row_len,local_rows,mlevs)   ! Grid being processed

REAL(wp), ALLOCATABLE ::  &
  work1(:,:,:)            & ! work array
 ,work2(:,:,:)              ! work array

REAL ::                   &
  rrhead(19)                ! real pp header

REAL ::                                      &
  field_out_full(in_cols,in_rows,mlevs)      & ! work array
 ,field_out_local(local_row_len,local_rows)    ! work array

REAL, ALLOCATABLE ::                  &
  t64(:,:,:), p64(:,:,:),rh64(:,:,:)  & ! Full precision
 ,dpdx_full(:,:)                      & ! work array
 ,dpdy_full(:,:)                        ! work array

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CRMSTYLE_READ_PP_INPUT'
CHARACTER(LEN=errormessagelength) :: iomessage

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Initialise variables

num_fields_found  = 0   ! fields found so far

num_fields_wanted = 31  ! added single level fields

read_code  = 0          ! return code
error_code = 0

! reset all the field flags to false

CALL reset_field_flags

! Only reading pp file on PE/Node 0 and scattering info to other nodes.

IF (mype == 0) THEN
  CALL get_env_var(pp_env, filename)

  CALL assign_file_unit(filename, unit_in, handler="fortran")

  OPEN(UNIT= unit_in, FILE=TRIM(filename), FORM= "UNFORMATTED",            &
      ACCESS = "SEQUENTIAL", ACTION = "READ", IOSTAT=  ErrorStatus,        &
      STATUS = "OLD", IOMSG=iomessage )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusFatal
    CALL EReport( RoutineName, ErrorStatus,"Cannot read pp file:"// newline//&
                  'IoMsg: '//TRIM(iomessage) )
  END IF
END IF

l_inter_uv = .TRUE.

IF (l_all_sea) THEN
  l_inter_th = .FALSE.
ELSE
  l_inter_th = .TRUE.
END IF

!----------------------------------------------------------------------------
! Do forever loop - until hit end of file or read in all fields for date/time
!----------------------------------------------------------------------------


major: DO
  ! read in next 2d or 3d field

  IF (mype == 0) THEN

    CALL read_next_pp_field(unit_in,stashcode,read_code,dmin,iihead,rrhead, &
                             field)

    IF (read_code /= 0 ) THEN
      WRITE(umMessage,'(A,I6)') 'Failed to read next field from PP file ', &
                                 read_code
      CALL umPrint(umMessage,src=RoutineName)
      ErrorStatus =  StatusFatal
      CALL EReport( RoutineName, ErrorStatus,"Cannot read next field" )
    END IF

    ! At this point the field is on the full model grid in 32 bit
    ! Want to scatter required part back to each PE/node after processing
    !-----------------------------------------------------------------------
    ! Decide how to process input field
    !-----------------------------------------------------------------------

    lsend  = .FALSE.
    levels = 1

    SELECT CASE (stashcode)

    CASE (2)    ! U wind
      IF (.NOT. l_u) THEN
        IF (l_ENDGame) THEN   ! Input ENDGame grid
          ! Put on p required grid
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, Field)
          DO k=1,mlevs
            DO j=1,in_rows
              field_out_full(in_cols,j,k) = 0.0   ! unset but will not be wanted
              DO i=1,in_cols-1
                field_out_full(i,j,k) =0.5* (Field(i,j,k) + Field(i+1,j,k))
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        ELSE                  ! Input New dynamics grid

          ! Put on p required grid
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, Field)
          DO k=1,mlevs
            DO j=1,in_rows
              field_out_full(1,j,k) = 0.0   ! unset but will not be wanted
              DO i=2,in_cols
                field_out_full(i,j,k) =0.5* (Field(i-1,j,k) + Field(i,j,k))
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF  ! test on grid
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (3)    ! v wind
      IF (.NOT. l_v) THEN
        IF (l_ENDGame) THEN   ! Input ENDGame grid
                 ! Put on p required grid
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, Field)
          DO k=1,mlevs
            DO i=1,in_cols
              field_out_full(i,1,k) = 0.0       ! unset but will not be wanted
              field_out_full(i,in_rows,k) = 0.0 ! unset but will not be wanted
            END DO
            DO j=2,in_rows-1
              DO i=1,in_cols
                field_out_full(i,j,k) =0.5* (Field(i,j,k) + Field(i,j+1,k))
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        ELSE                  ! Input New dynamics grid
                 ! Put on p required grid
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, Field)
          DO k=1,mlevs
            DO i=1,in_cols
              field_out_full(i,1,k) = 0.0       ! unset but will not be wanted
              field_out_full(i,in_rows,k) = 0.0 ! unset but will not be wanted
            END DO
            DO j=2,in_rows-1
              DO i=1,in_cols
                field_out_full(i,j,k) =0.5* (Field(i,j-1,k) + Field(i,j,k))
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF    ! test on grid
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4)    ! Theta
      IF (.NOT. l_theta) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (10)    ! q
      IF (.NOT. l_q) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12)    ! qcf
      IF (.NOT. l_qcf) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (24)   ! tstar
      IF (.NOT. l_tstar) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        ELSE
          WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END IF

    CASE (25)   ! boundary layer depth
      IF (.NOT. l_zh) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        ELSE
          WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END IF

    CASE (150)   ! w
      IF (.NOT. l_w) THEN
        ! extract required area
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs,in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (254)   ! qcl
      IF (.NOT. l_qcl) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (272)   ! qrain
      IF (.NOT. l_qrain) THEN
        ! extract required area
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (273)   ! qgraup
      IF (.NOT. l_qgraup) THEN
        ! extract required area
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (408)   ! p_theta_lev

      IF (.NOT. l_ptheta) THEN
        ! In the UM exner is interpolated not P so convert P to exner before
        ! interpolation to fixed heights.
        ! full grid - calculate exner pressure on theta levels work1
        ALLOCATE (work1(in_cols,in_rows,mlevs))
        ALLOCATE (work2(in_cols,in_rows,mlevs))

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED( mlevs, in_rows, in_cols, work1, Field, pref, kappa)
        DO k = 1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              work1(i,j,k) = (Field(i,j,k) /pref)**kappa
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
              ! Put on fixed heights - full grid work2 = exner on fixed heights
        CALL put_on_fixed_heights(in_cols,in_rows,mlevs,th_km_level_full,   &
                                work1, th_weight_full,l_inter_th,work2)

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                            &
!$OMP& SHARED(mlevs, in_rows, in_cols, th_km_level_full, field_out_full,  &
!$OMP&         work2, pref, kappa)
        DO k=1,mlevs
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
        DEALLOCATE(work1)

        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (409)   ! pstar
      IF (.NOT. l_pstar) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field)
        DO j = 1,in_rows
          DO i = 1,in_cols
            field_out_full(i,j,1) = field(i,j,1)
          END DO
        END DO
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = 1
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (1181)  ! dt1
      IF (.NOT. l_dt1) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (2181)  ! dt2
      IF (.NOT. l_dt2) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (3184)  ! dqcf3
      IF (.NOT. l_dqcf3) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(levels, mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,levels
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (3217)  ! sensible heat

      IF (.NOT. l_sh) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        ELSE
          WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END IF

    CASE (3234)  ! latent heat

      IF (.NOT. l_lh) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        ELSE
          WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
          CALL umPrint(umMessage,src=RoutineName)
        END IF
      END IF

    CASE (4203)  ! rain   Assuming 5 min output

      IF (.NOT. l_rain) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field, l_class_col, prec_full)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
              IF (l_class_col) THEN
                prec_full(i,j) = prec_full(i,j) +  field(i,j,1)
              END IF
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        END IF
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4204)  ! snow   Assuming 5 min output
      IF (.NOT. l_snow) THEN
        IF (dmin == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                         &
!$OMP& SHARED(in_rows, in_cols, field_out_full, field, l_class_col, prec_full)
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,1) = field(i,j,1)
              IF (l_class_col) THEN
                prec_full(i,j) = prec_full(i,j) +  field(i,j,1)
              END IF
            END DO
          END DO
!$OMP END PARALLEL DO
          lsend  = .TRUE.
          levels = 1
        END IF
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4181)  ! dt4
      IF (.NOT. l_dt4) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4182)  ! dq4
      IF (.NOT. l_dq4) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(levels, mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,levels
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4183)  ! dqcl4
      IF (.NOT. l_dqcl4) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (4184)  ! dqcf4
      IF (.NOT. l_dqcf4) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (9181)  ! dt9
      IF (.NOT. l_dt9) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (9182)  ! dq9
      IF (.NOT. l_dq9) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (9183)  ! dqcl9
      IF (.NOT. l_dqcl9) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12181)  ! dt12
      IF (.NOT. l_dt12) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12182)  ! dq12
      IF (.NOT. l_dq12) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12183)  ! dqcl12
      IF (.NOT. l_dqcl12) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12184)  ! dqcf12
      IF (.NOT. l_dqcf12) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                         &
!$OMP& SHARED(mlevs, in_rows, in_cols, field_out_full, field)
        DO k=1,mlevs
          DO j = 1,in_rows
            DO i = 1,in_cols
              field_out_full(i,j,k) = field(i,j,k)
            END DO
          END DO
        END DO    ! k loop
!$OMP END PARALLEL DO
        lsend  = .TRUE.
        levels = mlevs
      ELSE
        WRITE(umMessage,'(A,I6)') ' Field read in more ', stashcode
        CALL umPrint(umMessage,src=RoutineName)
      END IF


    CASE DEFAULT

      WRITE(umMessage,'(A,I6)') ' Field not wanted ', stashcode
      CALL umPrint(umMessage,src=RoutineName)

    END SELECT ! test on stashcode

    ! send back integer array with lsend, stashcode, and date info for field
    IF (lsend) THEN
      isend(1) = 1
      isend(2) = stashcode
      isend(3) = levels
      DO i=1,12
        isend(i+3) = iihead(i)
      END DO
    ELSE
      isend(:) = 0
    END IF

  END IF ! mype = 0


  icode = 0

  ! Force synchronisation before trying to scatter back fields and send
  ! integer information

  CALL  gc_gsync(nprocs,icode)

  ! Sends array isend to all PE/nodes to array isend
  mesg = 9001
  CALL gc_ibcast(mesg,15,0,nprocs,icode,isend)

  stashcode = isend(2)    ! set stashcode on all PEs/nodes

  IF (isend(1) == 1 .AND. icode == 0) THEN  ! complete processing of field

    levels = isend(3)

    ! scatter field back into output array
    DO k=1,isend(3)

      ! DEPENDS ON: scatter_field
      CALL scatter_field( field_out_local, field_out_full(1,1,k),    &
                 local_row_len,local_rows,                           &
                 in_cols,in_rows,                                    &
                 fld_type_p,halo_type_no_halo,                       &
                  0,gc_all_proc_group)

      ! copy to 32 field
      DO j=1,local_rows
        DO i=1,local_row_len
          data_grid(i,j,k) = field_out_local(i,j)
        END DO
      END DO

    END DO    ! k loop

    ! also need to calculate dpdx & dpdy on full grid and scatter
    ! back to PE/nodes
    IF (isend(2) == 408) THEN

      ! Need dp/dx and dp/dy
      dy=planet_radius*bdy*pi_over_180

      ALLOCATE(dpdx_full(in_cols,in_rows))
      ALLOCATE(dpdy_full(in_cols,in_rows))

      DO k=1,levels

        ! Calculate dp/dx  - first and last values on a row set to zero
        ! Full grid calculation can only be done on PE0 where all the
        ! information is held

        IF (mype == 0) THEN
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k,in_rows, in_cols, th_km_level_full, field_out_full,        &
!$OMP& dpdx_full, dy, xcoslat_full)
          DO j=1,in_rows
            ! Initialise to zero
            DO i=1,in_cols
              dpdx_full(i,j) = 0.0
            END DO
            DO i=2,in_cols-1
              IF ((th_km_level_full(i+1,j,k) >= 0) .AND.                   &
                  (th_km_level_full(i,j,k) >= 0)  .AND.                    &
                  (th_km_level_full(i-1,j,k) >= 0) ) THEN
                dpdx_full(i,j) = 0.5*(field_out_full(i+1,j,k)              &
                                    - field_out_full(i-1,j,k))/            &
                                           (dy*xcoslat_full(i,j))
              END IF
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
        icode = 0
        CALL  gc_gsync(nprocs,icode)

        ! DEPENDS ON: scatter_field
        CALL scatter_field( field_out_local, dpdx_full,               &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                   0,gc_all_proc_group)

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
          ! Initialise to zero
          DO j=1,in_rows
            DO i=1,in_cols
              dpdy_full(i,j) = 0.0
            END DO
          END DO
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                               &
!$OMP& SHARED(k,in_rows, in_cols, th_km_level_full, field_out_full,        &
!$OMP& dpdy_full, dy)
          DO j=2,in_rows-1
            DO i=1,in_cols
              IF ((th_km_level_full(i,j+1,k) >= 0) .AND.                     &
                  (th_km_level_full(i,j,k) >= 0)   .AND.                     &
                  (th_km_level_full(i,j-1,k) >= 0) ) THEN
                dpdy_full(i,j) = 0.5*( field_out_full(i,j+1,k)             &
                                      - field_out_full(i,j-1,k) )/dy
              END IF
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
        icode = 0
        CALL  gc_gsync(nprocs,icode)

        ! DEPENDS ON: scatter_field
        CALL scatter_field( field_out_local, dpdy_full,               &
                  local_row_len,local_rows,                           &
                  in_cols,in_rows,                                    &
                  fld_type_p,halo_type_no_halo,                       &
                   0,gc_all_proc_group)

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

    END IF   ! specail case of pressure


     ! complete processing of field on each PE/node based on stashcode
    SELECT CASE (isend(2))

    CASE (2)    ! U wind
      IF (.NOT. l_u) THEN
        ! Put on fixed heights
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,uv_km_level, &
                                  data_grid, uv_weight,l_inter_uv,u)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,             &
                                  num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_u = .TRUE.
      END IF

    CASE (3)    ! v wind
      IF (.NOT. l_v) THEN
        ! Put on fixed heights
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,uv_km_level, &
                                  data_grid, uv_weight,l_inter_uv,v)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_v = .TRUE.
      END IF

    CASE (4)    ! Theta
      IF (.NOT. l_theta) THEN
        ! For this field copy header info as correct P grid
        iyear = isend(3+1)
        imon  = isend(3+2)
        iday  = isend(3+3)
        ihour = isend(3+4)
        imin  = isend(3+5)
        isec  = isend(3+6)
        isyear = isend(3+7)
        ismon  = isend(3+8)
        isday  = isend(3+9)
        ishour = isend(3+10)
        ismin  = isend(3+11)
        issec  = isend(3+12)
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,theta)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_theta = .TRUE.
      END IF

    CASE (10)    ! q
      IF (.NOT. l_q) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,q)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_q = .TRUE.
      END IF

    CASE (12)    ! qcf
      IF (.NOT. l_qcf) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,qcf)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_qcf = .TRUE.
      END IF

    CASE (24)   ! tstar
      IF (.NOT. l_tstar) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            tstar(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_tstar = .TRUE.
      END IF

    CASE (25)   ! boundary layer depth
      IF (.NOT. l_zh) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            zh(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_zh = .TRUE.
      END IF

    CASE (150)   ! w
      IF (.NOT. l_w) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,w)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_w = .TRUE.
      END IF

    CASE (254)   ! qcl
      IF (.NOT. l_qcl) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,qcl)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_qcl = .TRUE.
      END IF

    CASE (272)   ! qrain
      IF (.NOT. l_qrain) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,qrain)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_qrain = .TRUE.
      END IF

    CASE (273)   ! qgraup
      IF (.NOT. l_qgraup) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,qgraup)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_qgraup = .TRUE.
      END IF

    CASE (408)   ! p_theta_lev
      IF (.NOT. l_ptheta) THEN
        ! Copy scattered field back (already on correct levels).
!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                               &
!$OMP& SHARED(mlevs, local_rows, local_row_len, p_theta_lev, data_grid)
        DO k=1,mlevs
          DO j=1,local_rows
            DO i=1,local_row_len
              p_theta_lev(i,j,k) = data_grid(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_ptheta = .TRUE.
      END IF

    CASE (409)   ! pstar
      IF (.NOT. l_pstar) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            pstar(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_pstar = .TRUE.
      END IF

    CASE (1181)  ! dt1
      IF (.NOT. l_dt1) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dt1)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dt1 = .TRUE.
      END IF

    CASE (2181)  ! dt2
      IF (.NOT. l_dt2) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dt2)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dt2 = .TRUE.
      END IF

    CASE (3184)  ! dqcf3
      IF (.NOT. l_dqcf3) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dqcf3)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcf3 = .TRUE.
      END IF

    CASE (3217)  ! sensible heat
      IF (.NOT. l_sh) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            sh(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_sh = .TRUE.
      END IF

    CASE (3234)  ! latent heat
      IF (.NOT. l_lh) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            lh(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_lh = .TRUE.
      END IF

    CASE (4203)  ! rain   Assuming 5 min output
      IF (.NOT. l_rain) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            rain(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_rain = .TRUE.
      END IF

    CASE (4204)  ! snow   Assuming 5 min output
      IF (.NOT. l_snow) THEN
        DO j=1,local_rows
          DO i=1,local_row_len
            snow(i,j) = data_grid(i,j,1)
          END DO
        END DO
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_snow = .TRUE.
      END IF

    CASE (4181)  ! dt4
      IF (.NOT. l_dt4) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dt4)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dt4 = .TRUE.
      END IF

    CASE (4182)  ! dq4
      IF (.NOT. l_dq4) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dq4)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dq4 = .TRUE.
      END IF

    CASE (4183)  ! dqcl4
      IF (.NOT. l_dqcl4) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dqcl4)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcl4 = .TRUE.
      END IF

    CASE (4184)  ! dqcf4
      IF (.NOT. l_dqcf4) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dqcf4)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcf4 = .TRUE.
      END IF

    CASE (9181)  ! dt9
      IF (.NOT. l_dt9) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dt9)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dt9 = .TRUE.
      END IF

    CASE (9182)  ! dq9
      IF (.NOT. l_dq9) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dq9)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dq9 = .TRUE.
      END IF

    CASE (9183)  ! dqcl9
      IF (.NOT. l_dqcl9) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dqcl9)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcl9 = .TRUE.
      END IF

    CASE (12181)  ! dt12
      IF (.NOT. l_dt12) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level, &
                                  data_grid,th_weight,l_inter_th,dt12)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dt12 = .TRUE.
      END IF

    CASE (12182)  ! dq12
      IF (.NOT. l_dq12) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,  &
                                  data_grid, th_weight,l_inter_th, dq12)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dq12 = .TRUE.
      END IF

    CASE (12183)  ! dqcl12
      IF (.NOT. l_dqcl12) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,  &
                                  data_grid, th_weight,l_inter_th, dqcl12)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcl12 = .TRUE.
      ELSE
        WRITE(umMessage,'(A,I6,I4)') ' Field read in more ', stashcode,  &
                                      num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
      END IF

    CASE (12184)  ! dqcf12
      IF (.NOT. l_dqcf12) THEN
        CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,  &
                                  data_grid, th_weight,l_inter_th, dqcf12)
        WRITE(umMessage,'(A,I6,I4)') ' Field read in ',stashcode,num_fields_found
        CALL umPrint(umMessage,src=RoutineName)
        num_fields_found = num_fields_found +1
        l_dqcf12 = .TRUE.
      END IF

    CASE DEFAULT

      WRITE(umMessage,'(A)') ' Field not wanted '
      CALL umPrint(umMessage,src=RoutineName)

    END SELECT ! test on stashcode

  END IF   ! test on whether field wanted so scattered back

  IF (num_fields_found == num_fields_wanted) THEN
    l_got_fields = .TRUE.
  END IF
  IF (read_code == 1 .OR. l_got_fields) EXIT major

END DO major

IF (mype == 0) THEN
  CLOSE(unit_in)
  CALL release_file_unit(unit_in, handler="fortran")
END IF

! calculate precip Not reading in now

DO j=1,local_rows
  DO i=1,local_row_len
    precip(i,j) = rain(i,j) + snow(i,j)
  END DO
END DO

!  thetav, density and t

!$OMP PARALLEL DO PRIVATE(i,j,k, exner) DEFAULT(NONE)                          &
!$OMP& SHARED(mlevs, local_rows, local_row_len, mask, p_theta_lev, t, density, &
!$OMP&     thetav, theta, q, qcl, qcf, qrain, qgraup, l_qgraup, pref, kappa,   &
!$OMP&     r, c_virtual)
DO k=1,mlevs
  DO j=1,local_rows
    DO i=1,local_row_len
      IF (mask(i,j,k)) THEN
        exner          = (p_theta_lev(i,j,k)/pref)**kappa
        t(i,j,k)       = theta(i,j,k)*exner
        density(i,j,k) = p_theta_lev(i,j,k)/(r*t(i,j,k))
        IF (l_qgraup) THEN
          thetav(i,j,k) = theta(i,j,k)*( 1.0 + c_virtual*q(i,j,k)        &
                              -qcl(i,j,k) - qcf(i,j,k) -qrain(i,j,k)     &
                              -qgraup(i,j,k) )
        ELSE
          thetav(i,j,k) = theta(i,j,k)*( 1.0 + c_virtual*q(i,j,k)        &
                              -qcl(i,j,k) - qcf(i,j,k) -qrain(i,j,k) )
        END IF
      ELSE
        t(i,j,k) = 300.0         ! set for call to qsat
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! relative humidity
! Problem call to qsat mix expects 64 bit numbers for rh, t, p

ALLOCATE (p64(local_row_len,local_rows,mlevs))
ALLOCATE (t64(local_row_len,local_rows,mlevs))
ALLOCATE (rh64(local_row_len,local_rows,mlevs))

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                 &
!$OMP& SHARED(mlevs, local_rows, local_row_len, t64, p64, t, p_theta_lev)

DO k=1,mlevs
  DO j=1,local_rows
    DO i=1,local_row_len
      t64(i,j,k) = t(i,j,k)
      p64(i,j,k) = p_theta_lev(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL qsat(rh64,t64,p64,local_row_len,local_rows,mlevs)

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                 &
!$OMP& SHARED(mlevs, local_rows, local_row_len, q, rh, rh64, mask)

DO k=1,mlevs
  DO j=1,local_rows
    DO i=1,local_row_len
      IF (mask(i,j,k)) THEN
        rh(i,j,k) =100.0*q(i,j,k)/rh64(i,j,k)  ! q/qsat
      ELSE
        rh(i,j,k) =0.0
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(rh64)
DEALLOCATE(t64)
DEALLOCATE(p64)

WRITE(umMessage,'(A,i4)') 'read in all fields required ', num_fields_found
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(3A50)')                                              &
                  ' theta     w         u        v         q        ', &
                  ' qcl       qcf       qrain    qgraup    t        ', &
                  ' p         rh        dpdx     dpdy      density  '

CALL umPrint(umMessage,src=RoutineName)
DO k=1,mlevs
  WRITE(umMessage,'(15E10.3)')                                               &
        theta(50,50,k),w(50,50,k),u(50,50,k),v(50,50,k),                     &
        q(50,50,k),qcl(50,50,k),qcf(50,50,k),qrain(50,50,k),qgraup(50,50,k), &
        t(50,50,k),p_theta_lev(50,50,k),rh(50,50,k),dpdx(50,50,k),           &
        dpdy(50,50,k),density(50,50,k)
  CALL umPrint(umMessage,src=RoutineName)
END DO
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_read_pp_input

END MODULE crmstyle_read_pp_input_mod
