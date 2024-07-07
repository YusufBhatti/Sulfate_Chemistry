! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  read from pp input file

MODULE crmstyle_read_ff_input_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_READ_FF_INPUT_MOD'

CONTAINS

SUBROUTINE crmstyle_read_ff_input(date_required,date_requiredp1,num_ff,ff_hdr)

USE hires_data_mod , ONLY:                                             &
  landsea, precip, rain, snow, zh, lh, sh, pstar, tstar,               &
  theta, thetav, t, rh, q, qcl, qcf,qrain, qgraup,                     &
  p_theta_lev,  u, v, w, dpdx, dpdy,   density,                        &
  dt1, dt2, dt4, dt9, dt12, dt30, dq4, dq9, dq12, dq30,                &
  dqcl4, dqcl9, dqcl12, dqcl30, dqcf4, dqcf3, dqcf12, dqcf30,          &
  drho, dqrain30, dqgr30

USE crmstyle_grid_info_mod, ONLY:                                          &
  local_row_len, local_rows

USE crmstyle_cntl_mod, ONLY:                                                &
  nx_start, ny_start, mlevs, model_levels, iprint, new_res,                 &
  num_want, stash_list, lev_list, date_typ, proc_list, l_qgraup,            &
  l_all_sea, l_sect30, num_want_no30

USE crmstyle_pp_data_mod, ONLY:                                             &
  iyear,imon,iday,ihour,imin,isec,isyear,ismon,isday,ishour,ismin,issec,    &
  itim,bdx,bdy

USE crmwork_arrays_mod, ONLY:                                            &
  h_theta_sea, mask, uv_km_level, th_km_level,                           &
  uv_weight, th_weight

USE planet_constants_mod, ONLY:                                      &
  kappa, cp, pref, r,  c_virtual

USE missing_data_mod, ONLY: rmdi

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE ereport_mod, ONLY: ereport, ereport_finalise

USE IO_Mod, ONLY: UM_Header_type,PP_Header_type

! subroutines
USE put_on_fixed_heights_mod, ONLY: put_on_fixed_heights
USE read_next_ffield_scat_mod, ONLY: read_next_ffield_scat
USE locate_fields_mod,        ONLY: locate_fields

USE qsat_mod, ONLY: qsat

USE umPrintMgr                              ! for writing output

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ---------------------------------------------------------------------------
! Description:
!   Read in pp file of model output
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3

! ---------------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::  &
  date_required(6)     &  ! date time for prognostics (time T)
 ,date_requiredp1(6)   &  ! date time for prognostics (time T+1 step)
 ,num_ff                  ! Number of fieldsfiles

TYPE(UM_Header_type), INTENT(INOUT) :: ff_hdr(num_ff) ! UM Headers: fieldsfile

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, ii,jj, ij, k,l   &  ! loop counters
 ,num_get                   ! Fields to get

INTEGER ::               &
  unit_in                & ! unit number for pp_file
 ,errorstatus            & ! error code
 ,read_code                ! read code

INTEGER ::               &
  ifile                  & ! file header number
 ,ilevs                    ! number of levels for field

LOGICAL ::               &
  l_inter_th             & ! interpolate height fields theta input
 ,l_inter_uv               ! interpolate height fields uv input


REAL(wp)  ::                            &
  field(local_row_len,local_rows,mlevs)  ! Full input field but just required
                                         ! levels


REAL(wp) ::      &
  dy             & !
 ,exner            ! exner

INTEGER ::              &
  fstart_num(num_want)  & ! start position of required field
 ,fend_num(num_want)    & ! end position of required fields
 ,file_num(num_want)      ! fieldsfile required fild is in

REAL, ALLOCATABLE :: &
  t64(:,:,:), p64(:,:,:),rh64(:,:,:)       ! Full precision

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CRMSTYLE_READ_FF_INPUT'

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

l_inter_uv = .TRUE.

IF (l_all_sea) THEN
  l_inter_th = .FALSE.
ELSE
  l_inter_th = .TRUE.
END IF

IF (l_sect30) THEN
  num_get = num_want
ELSE
  num_get = num_want_no30
END IF


! Already read in all fieldsfiles headers
! Need to find out where the required fields for this timestep are.

CALL locate_fields(num_ff,num_want,stash_list,lev_list, date_typ,         &
                   date_required, date_requiredp1, ff_hdr,                &
                   fstart_num,fend_num,file_num)

! Formats ok provided num_want <= 40 values
WRITE(umMessage,'(A)') ' locations of fields : start '
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (fstart_num(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' locations of fields : end '
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (fend_num(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' locations of fields : file '
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (file_num(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' locations of fields : proc '
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (proc_list(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)




! Loop to read in required fields and process

DO l = 1,num_get

  ifile = file_num(l)
  ilevs = lev_list(l)
  ErrorStatus = 0      ! initilise error status to none
  CALL read_next_ffield_scat(ilevs,fstart_num(l),fend_num(l),proc_list(l),  &
                             stash_list(l), ff_hdr(ifile),field,            &
                             ErrorStatus)

  ! Process  etc based on stashcode
  IF (ErrorStatus == 0 ) THEN

    SELECT CASE (stash_list(l))

    CASE (2)    ! U wind
      ! Put on fixed heights
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,uv_km_level,   &
                                field,uv_weight,l_inter_uv,u)

    CASE (3)    ! v wind

      ! Put on fixed heights
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,uv_km_level,   &
                                field,uv_weight,l_inter_uv,v)

    CASE (4)    ! Theta
      ! For this field copy header info as correct date time
      iyear = ff_hdr(ifile)%lookup(fstart_num(l))%ValidYear
      imon  = ff_hdr(ifile)%lookup(fstart_num(l))%ValidMonth
      iday  = ff_hdr(ifile)%lookup(fstart_num(l))%ValidDate
      ihour = ff_hdr(ifile)%lookup(fstart_num(l))%ValidHour
      imin  = ff_hdr(ifile)%lookup(fstart_num(l))%ValidMin
      isec  = ff_hdr(ifile)%lookup(fstart_num(l))%ValidSec
      isyear = ff_hdr(ifile)%lookup(fstart_num(l))%DataYear
      ismon  = ff_hdr(ifile)%lookup(fstart_num(l))%DataMonth
      isday  = ff_hdr(ifile)%lookup(fstart_num(l))%DataDate
      ishour = ff_hdr(ifile)%lookup(fstart_num(l))%DataHour
      ismin  = ff_hdr(ifile)%lookup(fstart_num(l))%DataMin
      issec  = ff_hdr(ifile)%lookup(fstart_num(l))%DataSec
      itim   = ff_hdr(ifile)%lookup(fstart_num(l))%LBTim

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,theta)

    CASE (10)    ! q

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,q)

    CASE (12)    ! qcf

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,qcf)

    CASE (24)   ! tstar

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, tstar, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          tstar(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (25)   ! boundary layer depth

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, zh, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          zh(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (150)   ! w

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,w)

    CASE (254)   ! qcl

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,qcl)

    CASE (272)   ! qrain

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,qrain)

    CASE (273)   ! qgraup

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,qgraup)


    CASE (407)   ! p_rho_lev
      WRITE(umMessage,'(A,i6)') ' Unwanted Field read in ', stash_list(l)
      CALL umPrint(umMessage,src=RoutineName)

    CASE (408)   ! p_theta_lev

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE) &
!$OMP& SHARED(mlevs, local_rows, local_row_len,  &
!$OMP&          mask, field, p_theta_lev, pref)
           ! Copy p into array for just required grid
      DO k=1,mlevs
        DO j=1,local_rows
          DO i=1,local_row_len
            IF (mask(i,j,k)) THEN
              p_theta_lev(i,j,k) = field(i,j,k)
            ELSE
              p_theta_lev(i,j,k) = pref    ! need values for qsat cal
            END IF
          END DO
        END DO
      END DO

        ! Calculate dp/dx and dp/dy  - done in reading ff
!$OMP END PARALLEL DO

    CASE (409)   ! pstar

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, pstar, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          pstar(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (1181)  ! dt1

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt1)

    CASE (2181)  ! dt2

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt2)

    CASE (3184)  ! dqcf3

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcf3)

    CASE (3217)  ! sensible heat

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, sh, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          sh(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (3234)  ! latent heat

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, lh, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          lh(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (4203)  ! rain

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, rain, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          rain(i,j) = field(i,j,1)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (4204)  ! snow
!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, snow, field)
      DO j=1,local_rows
        DO i=1,local_row_len
          snow(i,j) = field(i,j,12)
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (4181)  ! dt4

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt4)

    CASE (4182)  ! dq4

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dq4)

    CASE (4183)  ! dqcl4

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcl4)

    CASE (4184)  ! dqcf4
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcf4)

    CASE (9181)  ! dt9
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt9)

    CASE (9182)  ! dq9

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dq9)

    CASE (9183)  ! dqcl9
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcl9)

    CASE (12181)  ! dt12

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt12)

    CASE (12182)  ! dq12

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dq12)

    CASE (12183)  ! dqcl12
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcl12)

    CASE (12184)  ! dqcf12

      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcf12)

    CASE (30181)  ! dt30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dt30)

    CASE (30182)  ! dq30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dq30)

    CASE (30183)  ! dqcl30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                              field,th_weight,l_inter_th,dqcl30)

    CASE (30184)  ! dqcf30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                                field,th_weight,l_inter_th,dqcf30)

    CASE (30188)  ! drho
      ! Put on fixed heights from rho levels rather than theta
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,uv_km_level,   &
                                field,uv_weight,l_inter_uv,drho)

    CASE (30189)  ! dqrain30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                                field,th_weight,l_inter_th,dqrain30)

    CASE (30190)  ! dqgr30
      CALL put_on_fixed_heights(local_row_len,local_rows,mlevs,th_km_level,   &
                                field,th_weight,l_inter_th,dqgr30)

    CASE DEFAULT

      WRITE(umMessage,'(A,I6)') ' Field not wanted ', stash_list(l)
      CALL umPrint(umMessage,src=RoutineName)

    END SELECT ! test on stashcode

  END IF  ! test on read return code

END DO ! l fields required


! calculate total precip

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(local_rows, local_row_len, precip, rain, snow)
DO j=1,local_rows
  DO i=1,local_row_len
    precip(i,j) = rain(i,j) + snow(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

ALLOCATE (p64(local_row_len,local_rows,mlevs))
ALLOCATE (t64(local_row_len,local_rows,mlevs))
ALLOCATE (rh64(local_row_len,local_rows,mlevs))

!  thetav, density and t
! relative humidity
! Problem call to qsat mix expects 64 bit numbers for rh, t, p

!$OMP PARALLEL DO PRIVATE(i,j,k, exner) DEFAULT(NONE)                          &
!$OMP& SHARED(mlevs, local_rows, local_row_len, mask, p_theta_lev, t, density, &
!$OMP&     thetav, theta, q, qcl, qcf, qrain, qgraup, t64, p64, l_qgraup,      &
!$OMP&     pref, kappa, r, c_virtual)
DO k=1,mlevs
  DO j=1,local_rows
    DO i=1,local_row_len
      IF (mask(i,j,k)) THEN
        exner = (p_theta_lev(i,j,k)/pref)**kappa
        t(i,j,k) = theta(i,j,k)*exner
        density(i,j,k) = p_theta_lev(i,j,k)/(r*t(i,j,k))
        IF (l_qgraup) THEN
          thetav(i,j,k) = theta(i,j,k)*( 1.0 + c_virtual*q(i,j,k)    &
                          -qcl(i,j,k) - qcf(i,j,k) -qrain(i,j,k)     &
                          -qgraup(i,j,k) )
        ELSE
          thetav(i,j,k) = theta(i,j,k)*( 1.0 + c_virtual*q(i,j,k)    &
                             -qcl(i,j,k) - qcf(i,j,k) -qrain(i,j,k) )
        END IF
      ELSE
        t(i,j,k) = 300.0         ! set for call to qsat
        density(i,j,k) = 1.0     ! initialise to 1.0, should never be used
      END IF
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

IF (iprint == 1) THEN
  ! Prints info from a point in the centre of the find grid for the first
  ! coarse area as a check that data is being read in as sensible values.
  i=new_res(1)/2

  WRITE(umMessage,'(4A50)')                                              &
                    ' theta     w         u        v         q        ', &
                    ' qcl       qcf       qrain    qgraup    t        ', &
                    ' p         rh        dpdx     dpdy      density  ', &
                    ' thetav                                          '
  CALL umPrint(umMessage,src=RoutineName)

  DO k=1,mlevs
    WRITE(umMessage,'(16E10.3)') theta(i,i,k),w(i,i,k),u(i,i,k),v(i,i,k), &
        q(i,i,k),qcl(i,i,k),qcf(i,i,k),qrain(i,i,k),qgraup(i,i,k),        &
        t(i,i,k),p_theta_lev(i,i,k),rh(i,i,k),dpdx(i,i,k),                &
        dpdy(i,i,k),density(i,i,k),thetav(i,i,k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
  WRITE(umMessage,'(3A50)')                                              &
                    ' dt1       dt2       dt4      dt9       dt12     ', &
                    ' dq4       dq9       dq12     dqcl4     dqcl9    ', &
                    ' dqcl12    dqcf4     dqcf3    dqcf12             '

  CALL umPrint(umMessage,src=RoutineName)
  DO k=1,mlevs
    WRITE(umMessage,'(14E10.3)')                                         &
        dt1(i,i,k),dt2(i,i,k),dt4(i,i,k),dt9(i,i,k),                     &
        dt12(i,i,k),dq4(i,i,k),dq9(i,i,k),dq12(i,i,k),dqcl4(i,i,k),      &
        dqcl9(i,i,k),dqcl12(i,i,k),dqcf4(i,i,k),dqcf3(i,i,k),            &
        dqcf12(i,i,k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  WRITE(umMessage,'(2A50)')                                              &
                    ' rain      snow      total    pstar     tstar    ', &
                    ' lh        sh        zh       landsea            '

  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(9E10.3)') rain(i,i),snow(i,i),precip(i,i),pstar(i,i),    &
        tstar(i,i),lh(i,i),sh(i,i),zh(i,i),landsea(i,i)
  CALL umPrint(umMessage,src=RoutineName)

  IF (l_sect30) THEN
    WRITE(umMessage,'(3A50)')                                              &
                    ' dt30      dq30      dqcl30   dqcf30    drho     ',   &
                    ' dqrain30  dqgr30   '
    CALL umPrint(umMessage,src=RoutineName)
    DO k=1,mlevs
      WRITE(umMessage,'(7E10.3)')                                        &
        dt30(i,i,k),dq30(i,i,k),dqcl30(i,i,k),dqcf30(i,i,k),             &
        drho(i,i,k),dqrain30(i,i,k),dqgr30(i,i,k)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF
END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_read_ff_input

END MODULE crmstyle_read_ff_input_mod
