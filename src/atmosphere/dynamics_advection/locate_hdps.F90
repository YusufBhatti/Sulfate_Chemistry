! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE locate_hdps_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LOCATE_HDPS_MOD'

CONTAINS
SUBROUTINE locate_hdps(i_out, j_out, wght_xi1, wght_xi2,                    &
                       xi1_out, xi2_out, dsm1, i_off, h_f,                  &
                       pnt_type, dim_i, dim_j, dim_k, halo_i, halo_j)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_types, ONLY: integer32

USE horiz_grid_mod
USE lookup_table_mod
USE UM_ParParams
USE um_parcore,             ONLY: mype, nproc
USE model_domain_mod,       ONLY: l_regular, model_type, mt_global
USE umPrintMgr
USE ereport_mod,            ONLY: ereport, newline
USE errormessagelength_mod, ONLY: errormessagelength
USE errorurl_mod,           ONLY: errorurl1, errorurl2

USE dynamics_testing_mod,   ONLY: L_trap_uv

IMPLICIT NONE

!
! Description:
!
!          Locate the position of the departure points
!          and return the integer location of the pivot point
!          and offset/weight for use in interpolation.
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER, INTENT(IN)  :: dim_i, dim_j, dim_k, pnt_type, dsm1

INTEGER(KIND=integer32), INTENT(OUT) ::                                     &
                        i_out(dim_i, dim_j, dim_k),                         &
                        j_out(dim_i, dim_j, dim_k)

REAL,    INTENT(IN)  :: xi1_out(dim_i, dim_j, dim_k),                       &
                        xi2_out(dim_i, dim_j, dim_k)

REAL,    INTENT(OUT) :: wght_xi1(dim_i, dim_j, dim_k),                      &
                        wght_xi2(dim_i, dim_j, dim_k)

INTEGER, INTENT(IN)  :: halo_i, halo_j, i_off, h_f

INTEGER              :: i, j, k, id
REAL                 :: x_off, y_off, temp
REAL                 :: rdxi1, rdxi2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER        :: RoutineName='LOCATE_HDPS'
CHARACTER(LEN=errormessagelength)  :: cmessage

INTEGER :: no_lim_pnts(4), ierr, icode

REAL    :: real_temp_var
INTEGER :: my_floor
my_floor(real_temp_var) = INT(real_temp_var - 0.5 + SIGN(0.5,real_temp_var))


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF ( l_regular ) THEN

  rdxi1 = 1.0/delta_xi1
  rdxi2 = 1.0/delta_xi2

  SELECT CASE(pnt_type)
  CASE (fld_type_u)
    x_off = 1.0
    y_off = 0.5
  CASE (fld_type_v)
    x_off = 0.5
    y_off = 1.0
  CASE (fld_type_w, fld_type_p)
    x_off = 0.5
    y_off = 0.5
  END SELECT

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(static) PRIVATE(k,j,i, temp)    &
!$OMP              SHARED(dim_k, dim_j, dim_i, dsm1,                      &
!$OMP                     xi1_out,base_xi1,rdxi1,i_out,wght_xi1,x_off,    &
!$OMP                     xi2_out,base_xi2,rdxi2,j_out,wght_xi2,y_off)
  !cdir collapse
  DO k = 1, dim_k
    DO j = 1, dim_j
      DO i = 1, dim_i
        temp            = (xi1_out(i,j,k)-base_xi1)*rdxi1
        i_out(i,j,k)    = my_floor(x_off + temp)
        wght_xi1(i,j,k) = temp - ( i_out(i,j,k) - x_off )
      END DO
    END DO

    DO j = 1, dim_j
      DO i = 1, dim_i
        temp            = (xi2_out(i,j,k) - base_xi2)*rdxi2
        j_out(i,j,k)    = my_floor(y_off + temp)
        wght_xi2(i,j,k) = temp - ( j_out(i,j,k) - y_off)
        j_out(i,j,k)    = j_out(i,j,k) - dsm1
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE  ! variable resolution

  IF ( pnt_type == fld_type_u ) THEN
    rdxi1 = 1.0/delta_xi1_u
  ELSE
    rdxi1 = 1.0/delta_xi1_p
  END IF

  IF ( pnt_type == fld_type_v ) THEN
    rdxi2 = 1.0/delta_xi2_v
  ELSE
    rdxi2 = 1.0/delta_xi2_p
  END IF

  ! DO xi1 direction first

!$OMP  PARALLEL DO  DEFAULT(NONE) SCHEDULE(static) PRIVATE(k,j,i,temp,id)   &
!$OMP               SHARED(dim_k, dim_j, dim_i, dsm1,                       &
!$OMP                      xi1_out,base_xi1_u,base_xi1_p, rdxi1,            &
!$OMP                      xi2_out,base_xi2_v,base_xi2_p, rdxi2,            &
!$OMP                      i_out,wght_xi1, j_out,wght_xi2,                  &
!$OMP                      glob_xi1_u, glob_xi1_p, glob_xi2_v, glob_xi2_p,  &
!$OMP                      glob_rdxi1_u, glob_rdxi1_p, glob_rdxi2_v,        &
!$OMP                      glob_rdxi2_p, i_lkup_u, i_lkup_p,                &
!$OMP                      j_lkup_v, j_lkup_p, pnt_type)
  DO k = 1, dim_k
    IF (pnt_type == fld_type_u) THEN
      DO j = 1, dim_j
        DO i = 1, dim_i
          temp            = (xi1_out(i,j,k)-base_xi1_u)*rdxi1
          id              = i_lkup_u( INT(temp) )
          wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_u(id))               &
                            *glob_rdxi1_u(id)
          IF ( wght_xi1(i,j,k) > 1.0 ) THEN
            id              = id + 1
            wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_u(id))             &
                                *glob_rdxi1_u(id)
          END IF
          i_out(i,j,k) = id + 1 ! offset for u -> p
        END DO
      END DO
    ELSE
      DO j = 1, dim_j
        DO i = 1, dim_i
          temp            = (xi1_out(i,j,k)-base_xi1_p)*rdxi1
          id              = i_lkup_p( INT(temp) )
          wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_p(id))               &
                            *glob_rdxi1_p(id)
          IF ( wght_xi1(i,j,k) > 1.0 ) THEN
            id              = id + 1
            wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_p(id))             &
                                *glob_rdxi1_p(id)
          END IF
          i_out(i,j,k) = id
        END DO
      END DO
    END IF

    ! Now xi2

    IF (pnt_type == fld_type_v) THEN
      DO j = 1, dim_j
        DO i = 1, dim_i
          temp            = (xi2_out(i,j,k)-base_xi2_v)*rdxi2
          id              = j_lkup_v( INT(temp) )
          wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_v(id))               &
                            *glob_rdxi2_v(id)
          IF ( wght_xi2(i,j,k) > 1.0 ) THEN
            id              = id + 1
            wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_v(id))             &
                              *glob_rdxi2_v(id)
          END IF
          j_out(i,j,k) = id + 1 - dsm1 ! offset for v -> p
        END DO
      END DO
    ELSE
      DO j = 1, dim_j
        DO i = 1, dim_i
          temp            = (xi2_out(i,j,k)-base_xi2_p)*rdxi2
          id              = j_lkup_p( INT(temp) )
          wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_p(id))               &
                            *glob_rdxi2_p(id)
          IF ( wght_xi2(i,j,k) > 1.0 ) THEN
            id              = id + 1
            wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_p(id))             &
                              *glob_rdxi2_p(id)
          END IF
          j_out(i,j,k) = id - dsm1
        END DO
      END DO
    END IF
  END DO
!$OMP END PARALLEL DO

END IF ! l_regular

IF ( L_trap_uv ) THEN
  no_lim_pnts(:) = 0

  ! Now limit departure points
!$OMP  PARALLEL DO  DEFAULT(NONE) SCHEDULE(static) PRIVATE(k,j,i,id)   &
!$OMP               SHARED(dim_k, dim_j, dim_i, i_off, h_f,            &
!$OMP                      halo_i, halo_j, i_out, j_out,               &
!$OMP                      wght_xi1, wght_xi2, model_type)             &
!$OMP               REDUCTION(+:no_lim_pnts)
  DO k = 1, dim_k
    DO j = 1, dim_j
      DO i = 1, dim_i

! find global value of i (j is already local due to no coms-on-demand
!                         in the N/S directions)

        id = i + i_off - 1

! Note when clipping h_f is the stencil length (depending on interpolation
! type) in the positive direction and h_f-1 is the length in the minus
! direction.
! The interpolation weights are altered for clipped points to be 0.5
! which should reduce biasing but this is a nicety and we could have
! left them unaltered (setting them to 0.0 or 1.0 would obvious bias
! the interpolations)

! If global don't limit E-W

        IF( model_type /= mt_global ) THEN

! Clip East if required
          IF ( i_out(i,j,k) - id > halo_i-h_f ) THEN
            i_out(i,j,k) = id + halo_i - h_f
            wght_xi1(i,j,k) = 0.5
            no_lim_pnts(1) = no_lim_pnts(1) + 1
          END IF

! Clip West if required
          IF ( id - i_out(i,j,k) > halo_i-h_f+1 ) THEN
            i_out(i,j,k) = id - halo_i + h_f-1
            wght_xi1(i,j,k) = 0.5
            no_lim_pnts(2) = no_lim_pnts(2) + 1
          END IF

        END IF

! Clip North if required
        IF ( j_out(i,j,k) - j > halo_j-h_f ) THEN
          j_out(i,j,k) = j + halo_j - h_f
          wght_xi2(i,j,k) = 0.5
          no_lim_pnts(3) = no_lim_pnts(3) + 1
        END IF

! Clip South if required
        IF ( j - j_out(i,j,k) > halo_j-h_f+1 ) THEN
          j_out(i,j,k) = j - halo_j + h_f-1
          wght_xi2(i,j,k) = 0.5
          no_lim_pnts(4) = no_lim_pnts(4) + 1
        END IF

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

! Output slipping diagnostics if required; i.e. the accumulated number
! of clipped points in each direction. Note the printstatus level here
! needs to match that in eg_adjust_vert_bounds to avoid deadlock

  IF (PrintStatus == PrStatus_Diag) THEN
    CALL gc_isum(4,nproc,ierr,no_lim_pnts)
    IF ( mype == 0 .AND. SUM(no_lim_pnts) > 0 ) THEN
      SELECT CASE(pnt_type)
      CASE (fld_type_u)
        WRITE(umMessage,                                                    &
              '("Limiting u-pnts  E:", I8," W:",I8," N:",I8," S",I8)')      &
                          no_lim_pnts(:)
      CASE (fld_type_v)
        WRITE(umMessage,                                                    &
              '("Limiting v-pnts  E:", I8," W:",I8," N:",I8," S",I8)')      &
                          no_lim_pnts(:)
      CASE (fld_type_w)
        WRITE(umMessage,                                                    &
              '("Limiting w-pnts  E:", I8," W:",I8," N:",I8," S",I8)')      &
                          no_lim_pnts(:)
      CASE (fld_type_p)
        WRITE(umMessage,                                                    &
              '("Limiting p-pnts  E:", I8," W:",I8," N:",I8," S",I8)')      &
                          no_lim_pnts(:)
      END SELECT
      CALL umPrint(umMessage,src=RoutineName)

    END IF
  END IF

ELSE ! If global and not trapping the wind, we will  fail if North/South 
     ! point is outside of the halo

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                           &
!$OMP             PRIVATE(i,j,k,cmessage,icode)                            &
!$OMP             SHARED(dim_k, dim_j, dim_i, j_out, halo_j, h_f,          &
!$OMP                    model_type, newline, errorurl1, errorurl2)
  DO k = 1, dim_k
    DO j = 1, dim_j
      DO i = 1, dim_i
        IF ( j_out(i,j,k) - j > halo_j-h_f  .OR.                           &
             j - j_out(i,j,k) > halo_j-h_f+1 ) THEN
             icode = 15
             IF( model_type == mt_global ) THEN
                cmessage='North/South halos too small for advection.'
             ELSE
                cmessage='Halos too small for advection: set l_trap_uv'    &
                        //' or increase halo size.'
             END IF
             ! Include link to documentation for known UM failure points.
             cmessage = TRIM(cmessage)//newline                            &
                        //TRIM(errorurl1)//newline//TRIM(errorurl2)
             CALL ereport('LOCATE_HDPS',icode,cmessage)
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE locate_hdps
END MODULE locate_hdps_mod
