! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: STEXTC ---------------------------------------------------
!
! Purpose: Extracts a weighted subfield within a region specified
!          by a lower left hand and upper right hand corner.
!          Single level at a time. (STASH service routine).
!
! External documentation:
!   Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                system (STASH)
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH
SUBROUTINE stextc(fieldin,vx,vy,fld_type,halo_type,               &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  fieldout,                                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
                  icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
INTEGER ::                                                        &
    vx,vy,                                                        &
                                          ! IN  input field size
    fld_type,                                                     &
                                          ! IN  field type(u/v/p)
    halo_type,                                                    &
                                          ! IN  halo type
    xstart,ystart,                                                &
                                          ! IN  lower LH corner
    xend,yend,                                                    &
                                          ! IN  upper RH corner
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
    fieldin(vx,vy),                                               &
                                          ! IN  input field
    in_aux(vx+1,vy),                                              &
                                          ! IN  input field
    fieldout(xstart:xend,ystart:yend),                            &
                                          ! OUT output field
    out_aux(xstart:xend+1,ystart:yend),                           &
    pstar_weight(vx+1,vy),                                        &
                                          ! IN  pstar mass weight
! (already interpolated to the correct grid and
!  set to 1.0 where no mass weighting is required)
    area_weight(vx+1,vy),                                         &
    rmdi                                  ! IN  missing data indic
! ----------------------------------------------------------------------
!
! Local variables
!
INTEGER :: i,j,ii   ! ARRAY INDICES FOR VARIABLE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STEXTC'


! ----------------------------------------------------------------------

! Calculate the output field, by multiplying the input field by
! pstar_weight and area_weight. These arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
in_aux(1:vx,:) = fieldin(1:vx,:)

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
      out_aux(i,j) =                                            &
        in_aux(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
    ELSE
      out_aux(i,j)=rmdi
    END IF
  END DO

END DO

fieldout(:,:) = out_aux(xstart:xend,ystart:yend)


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE stextc
