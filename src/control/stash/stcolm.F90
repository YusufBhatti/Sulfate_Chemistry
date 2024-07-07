! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: STCOLM ---------------------------------------------------
!
!  Purpose: Calculate weighted column mean within a region specified
!           by a lower left hand and upper right hand corner.
!           (STASH service routine).
!
!  Programming standard: UM Doc Paper 3
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!  ---------------------------------------------------------------
!
!  Interface and arguments: ------------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH

SUBROUTINE stcolm(fieldin,vx,vy,vz,fld_type,halo_type,            &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  fieldout,index_lev,zsize,                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE sterr_mod, ONLY: st_illegal_weight, unknown_weight
USE stparam_mod, ONLY: stash_weight_null_code, &
                 stash_weight_area_code, stash_weight_mass_code

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
INTEGER ::                                                        &
    vx,vy,vz,                                                     &
                                          ! IN  input field size
    fld_type,                                                     &
                                          ! IN  field type(u/v/p)
    halo_type,                                                    &
                                          ! IN  halo type
    xstart,ystart,                                                &
                                          ! IN  lower LH corner
    xend,yend,                                                    &
                                          ! IN  upper RH corner
    zsize,                                                        &
                                ! IN no of horiz levels to process
    index_lev(zsize),                                             &
                                ! IN offset for each horiz level
    level_code,                                                   &
                                          ! IN  input level code
    mask_code,                                                    &
                                          ! IN  masking code
    weight_code,                                                  &
                                          ! IN  weighting code
    icode                                 ! OUT error return code
CHARACTER(LEN=errormessagelength) ::                              &
    cmessage                              ! OUT error return msg
LOGICAL ::                                                        &
    lwrap,                                                        &
                                          ! IN  TRUE if wraparound
    lmasswt,                                                      &
                                          ! IN  TRUE if masswts OK
    mask(vx+1,vy)                         ! IN  mask array
REAL ::                                                           &
    fieldin(vx,vy,vz),                                            &
                                          ! IN  input field
    fieldout(xstart:xend,ystart:yend),                            &
                                          ! OUT output field
    pstar_weight(vx+1,vy,zsize),                                  &
                                          ! IN  mass weight factor
    area_weight(vx+1,vy),                                         &
                                          ! IN  area weight factor
    rmdi                                  ! IN  missing data indic
! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER :: i,j,k,ii,kk       ! ARRAY INDICES FOR VARIABLE

REAL :: sumctop(xstart:xend,ystart:yend)
REAL :: sumcbot(xstart:xend,ystart:yend)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STCOLM'

! ----------------------------------------------------------------------
!  0. Initialise sums
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO i=xstart,xend
  DO j=ystart,yend
    sumctop(i,j)=0.0
    sumcbot(i,j)=0.0
  END DO
END DO
! ----------------------------------------------------------------------
!  1. Form column sums
!
!  1.1 NULL weighting or area weighting
!
IF (weight_code == stash_weight_null_code .OR.                     &
    weight_code == stash_weight_area_code) THEN
  DO kk=1,zsize
    k=index_lev(kk)
    DO i=xstart,xend

      IF ( lwrap .AND.                                            &
          (i  >   (lasize(1,fld_type,halo_type)-                  &
                   halosize(1,halo_type)))) THEN
        ! miss halos on wrap around
        ii=i-blsize(1,fld_type)
      ELSE
        ii=i
      END IF

      DO j=ystart,yend
        sumcbot(i,j)=sumcbot(i,j)+1.0
        sumctop(i,j)=sumctop(i,j)+fieldin(ii,j,k)
      END DO
    END DO
  END DO
  !
  !  1.2 mass weighting
  !
ELSE IF (weight_code == stash_weight_mass_code) THEN
  IF (.NOT. lmasswt) THEN
    ! Mass-weighting on level types with no mass-weight defined is not
    ! supported - should be prevented by UI
    cmessage='STCOLM  : mass-weights not defined for this diag'
    icode=st_illegal_weight
    GO TO 999
  ELSE
    DO kk=1,zsize
      k=index_lev(kk)
      DO i=xstart,xend
        IF ( lwrap .AND.                                            &
            (i  >   (lasize(1,fld_type,halo_type)-                  &
                     halosize(1,halo_type)))) THEN
          ! miss halos on wrap around
          ii=i-blsize(1,fld_type)
        ELSE
          ii=i
        END IF
        DO j=ystart,yend
          sumcbot(i,j)=sumcbot(i,j) +  pstar_weight(ii,j,kk)
          sumctop(i,j)=sumctop(i,j) +  fieldin(ii,j,k)*           &
                                       pstar_weight(ii,j,kk)
        END DO
      END DO
    END DO
  END IF
ELSE
  cmessage='STCOLM  : Invalid weighting code detected'
  icode=unknown_weight
  GO TO 999
END IF
! ----------------------------------------------------------------------
!  2. Perform masking (set missing data at masked points) - compute mean
!
DO i=xstart,xend
  IF ( lwrap .AND.                                            &
      (i  >   (lasize(1,fld_type,halo_type)-                  &
               halosize(1,halo_type)))) THEN
    ! miss halos on wrap around
    ii=i-blsize(1,fld_type)
  ELSE
    ii=i
  END IF
  DO j=ystart,yend
    IF (mask(ii,j)) THEN
      fieldout(i,j)=sumctop(i,j)/sumcbot(i,j)
    ELSE
      fieldout(i,j)=rmdi
    END IF
  END DO
END DO
!
999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stcolm
