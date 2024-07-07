! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine to expand ozone from the ancillary to a full field.
!     Subroutine Interface:
SUBROUTINE O3_to_3D(lexpand_ozone, i_ozone_int,                   &
  rows, row_length, model_levels, ozone_levels,                   &
  halo_i, halo_j, off_x, off_y, at_extremity,                     &
  z_top_of_model,                                                 &
  theta, exner_theta_levels,                                      &
  rho,   exner_rho_levels,                                        &
  nd_o3, ozone_in,                                                &
  min_trop_level, max_trop_level,                                 &
  L_O3_trop_level,L_O3_trop_height,L_T_trop_level,L_T_trop_height,&
  O3_trop_level,O3_trop_height,T_trop_level,T_trop_height,        &
  proc_row_group,                                                 &
  global_row_length,                                              &
  ozone3D,                                                        &
  ErrorStatus, cmessage                                           &
  )
!
USE level_heights_mod, ONLY:                                      &
                r_theta_levels, r_rho_levels,                     &
                eta_theta_levels, eta_rho_levels
USE atm_fields_bounds_mod, ONLY: o3dims2
USE global_2d_sums_mod, ONLY: global_2d_sums
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE Field_Types
USE o3crits_mod, ONLY: o3_grad_crit, o3_conc_crit, o3_strat_crit
USE o3intp_mod, ONLY: io3_3dspec, io3_2dspec, io3_2dmasscon,      &
                      io3_trop_map, io3_trop_map_masscon

USE errormessagelength_mod, ONLY: errormessagelength

USE tropin_mod, ONLY: tropin
IMPLICIT NONE
!
! Description:
!   This routine takes the ozone field supplied in the ancillary
!   and converts it to a full 2-D (zonal mean) or 3-D field
!   as directed by the option selected.
!
! Method:
!   Essentially, the routine must carry out vertical interpolation
!   of the ozone field supplied, typically from an ancillary file,
!   on to the vertical levels at a grid-point. Various options are
!   permitted, the newer ones allowing interpolation in height to
!   match the new vertical structure of the model.
!      This code is largely developmental and will subsequently be
!   combined with changes and extensions to the ancillary system
!   for ozone.
!      The ozone concentrations are found on theta levels,
!   and the ozone tropopause is found on a rho level.
!   Rho level 2 is found between theta levels 1 and 2.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code description:
!   FORTRAN 90
!   This code is written to the programming standards of version 6
!   of UMDP3.
!
!     ----------------------------------------------------------------
!
!     Input arguments:
!
INTEGER, INTENT(IN) :: proc_row_group, global_row_length
! Specify the total number of gridpoints in the east-west direction
! and which processor takes care of which row.
LOGICAL, INTENT(IN) :: lexpand_ozone
!       Flag for expansion of ozone from the ancillary
INTEGER, INTENT(IN) :: i_ozone_int
!         Method of expanding the ozone in the dump to a full 3-D field
INTEGER, INTENT(IN) :: rows
!         Number of EW rows
INTEGER, INTENT(IN) :: row_length
!         Number of points on each row
INTEGER, INTENT(IN) :: halo_i
!         Size of large EW Halo
INTEGER, INTENT(IN) :: halo_j
!         Size of large NS Halo
INTEGER, INTENT(IN) :: off_x
!         Size of EW Halo
INTEGER, INTENT(IN) :: off_y
!         Size of NS Halo
LOGICAL, INTENT(IN) :: at_extremity(4)
!         Flags to indicate whether the processor is at the boundary
!         of the domain
!
LOGICAL, INTENT(IN) :: L_O3_trop_level, L_O3_trop_height
LOGICAL, INTENT(IN) :: L_T_trop_level, L_T_trop_height
! Flags indicating which diagnostics are on or off.
INTEGER, INTENT(IN) :: model_levels
!         Number of vertical levels in the atmosphere
INTEGER, INTENT(IN) :: ozone_levels
!         Number of levels on which ozone is specified: these are
!         contiguous levels reaching to the top of the model. If
!         used with the options IO3_3DSPEC or IO3_2DSPEC the value
!         on the lowest level of this array is used to set the mixing
!         ratio on lower layers of the model. Its interpretation under
!         other options is not yet finalized (but note that it is
!         used in radiation, so any changes here will need to be
!         accounted for there).
INTEGER, INTENT(IN) :: nd_o3
!         Size of the array of ozone supplied from D1
!
INTEGER, INTENT(IN) :: min_trop_level
!         Minimum permitted level of the tropopause
INTEGER, INTENT(IN) :: max_trop_level
!         Maximum permitted level of the tropopause
!
REAL, INTENT(IN)    :: z_top_of_model
!         Height of the top of the model
!
REAL, INTENT(IN)    :: exner_theta_levels(                        &
                         1-off_x:row_length+off_x,                &
                         1-off_y:rows+off_y,                      &
                         model_levels)
!         Exner function on theta levels
REAL, INTENT(IN)    :: theta(1-off_x:row_length+off_x,            &
                             1-off_y:rows+off_y,                  &
                             model_levels)
!         Potential temperatures on theta-levels
!
REAL, INTENT(IN)    :: exner_rho_levels(                          &
                         1-off_x:row_length+off_x,                &
                         1-off_y:rows+off_y,                      &
                         model_levels)
!         Exner function on rho levels
!
REAL, INTENT(IN)    :: rho(1-off_x:row_length+off_x,              &
                           1-off_y:rows+off_y,                    &
                           model_levels)
!         Atmospheric densities on rho-levels
!
REAL, INTENT(IN)    :: ozone_in(nd_o3)
!         Input ozone field in D1
!
! SCM Dummy variables to keep call to tropin consistent.
REAL ::                                                          &
scm_dummy_1d(1,1)                                                &
,scm_dummy_2d(1,1,0:model_levels)

!     Output arguments:

REAL, INTENT(OUT)   :: ozone3D(row_length, rows, ozone_levels)
!         Expanded ozone field
!
REAL, INTENT(OUT) :: T_trop_level(row_length,rows)
REAL, INTENT(OUT) :: O3_trop_level(row_length,rows)
! Points to the lower boundary of the layer containing the thermal and
! ozone tropopause.
REAL, INTENT(OUT) :: O3_trop_height(row_length,rows)
REAL, INTENT(OUT) :: T_trop_height(row_length,rows)
! Height of the ozone and thermal tropopause.
!
!     Error Status:
INTEGER, INTENT(INOUT):: ErrorStatus
!         Error code
CHARACTER(LEN=errormessagelength), INTENT(INOUT) :: cmessage
!         Short error message
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!     Local variables:
INTEGER :: i
!         Loop variable for Latitude
INTEGER :: j
!         Loop variable for Longitude
INTEGER :: k
!         Loop variable for vertical
INTEGER :: jk
!         Loop variable
INTEGER :: jkp1
!         Loop variable
INTEGER :: ijk
!         Loop variable
INTEGER :: k_o3_trop(row_length)
!         Index of theta level just below ozone tropopause
LOGICAL :: k_o3_trop_done(row_length)
INTEGER :: k_anc(row_length)
!         Index of immediate theta level below the interpolated height
INTEGER :: k_off
!         Offset of the ozone levels from the model levels
INTEGER :: trindx(row_length, rows)
!         Points to the lower boundary of the layer containing the
!         thermal tropopause.
REAL :: height_above_surf
!         Height above the surface of the current grid-point
!         (on theta level).
REAL :: z_trop_O3(row_length)
!         Height of the ozone tropopause above the surface in the
!         ancillary profile.
REAL :: z_trop_pt(row_length)
!         Height of the thermal tropopause above the surface at the
!         current grid-point.
REAL :: z_int(row_length)
!         Interpolated height in the ancillary profile (on no level).
REAL :: T_n(row_length, rows, model_levels)
!         Temperatures on theta levels
REAL :: rho_anc(model_levels)
REAL :: layer_mass(row_length,rows)
! Mass on rho_levels
REAL :: avg_layer_mass(rows,model_levels)
! Mass on rho_levels averaged over longitude.
REAL :: ww1(rows)
! Working Array
!         Atmospheric densities on rho-levels in the ancillary profile
REAL :: o3_mass_anc(0: model_levels)
!         Column integrated mass of ozone in the ancillary profile
!         on rho levels, the unit is kg per kg per m^2
REAL :: o3_mass_cumul(row_length)
!         Column mass of ozone at the grid-point integrated from
!         the surface to the top of the current layer
REAL :: o3_mass_cumul_below(row_length)
!         Column mass of ozone at the grid-point integrated from
!         the surface to the top of the layer below the current one
REAL :: o3_gradient_anc(row_length, rows, model_levels)
!         Ozone concentration gradient on rho-levels in ancillary file.
!
!
!- End of header
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!
! Additional loop variable for vectorization

INTEGER :: kl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='O3_TO_3D'
!
!     Preliminary calculations:-
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF ( (i_ozone_int  ==  io3_trop_map) .OR.                          &
     (i_ozone_int  ==  io3_trop_map_masscon) ) THEN
  !
  !       The (thermal) tropopause is required: this code is
  !       copied from RAD_CTL, T_n being replaced by
  !       its definition.
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        T_n(i,j,k)                                            &
          =theta(i,j,k)*exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  !
  CALL tropin (T_n, exner_rho_levels, exner_theta_levels,         &
                 row_length, rows, model_levels, off_x, off_y,    &
                 at_extremity,scm_dummy_1d,scm_dummy_2d,          &
                 min_trop_level, max_trop_level, trindx )
  !

  ! Calculate the mass in each layer between rho-levels.

  DO k = 1,model_levels

    DO j = 1, rows
      DO i = 1,row_length
        IF (k == 1) THEN
          layer_mass(i,j)=rho(i,j,k+1)                          &
             *(r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1))
        END IF
        IF ((k >  1) .AND. (k <  model_levels)) THEN
          layer_mass(i,j)=0.5*(rho(i,j,k+1)+rho(i,j,k))         &
             *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
        END IF
        IF (k == model_levels) THEN
          layer_mass(i,j)=rho(i,j,k)                            &
             *(r_theta_levels(i,j,k)-r_rho_levels(i,j,k))
        END IF
      END DO
    END DO

    ! Average the mass per layer over longitude.

    CALL global_2d_sums(layer_mass, row_length, 1, 0, 0, rows,     &
                        ww1, proc_row_group)

    DO j=1, rows
      avg_layer_mass(j,k)=ww1(j)/global_row_length
    END DO
  END DO
  !
END IF
!
!======================================================
! Set SCM dummy values to zero
scm_dummy_1d(:,:)   = 0.0
scm_dummy_2d(:,:,:) = 0.0
!======================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!
SELECT CASE(i_ozone_int)
  !
CASE (io3_3dspec)
  !
  !       In this case we should have a 3D ancillary file so
  !       we return an error message if not.
  IF (lexpand_ozone) THEN
    ErrorStatus=123

    CALL Ereport("O3_to_3D", ErrorStatus,                         &
       "A 2D ozone ancillary has been specified with a" //        &
       "3D ozone option.")
  END IF
  !
  !       The ozone field supplied is the intended 3-D field and is
  !       simply copied across. (N.B. This is inefficient: in the future
  !       we should be able to avoid this with a pointer and an
  !       allocatable array.)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,ijk)                 &
!$OMP SHARED(ozone_levels,rows,row_length,o3dims2,ozone3D,ozone_in)
  DO k = 1, ozone_levels
    DO j = 1, rows
      DO i = 1, row_length
        ijk = i + (j-1)*row_length +                              &
              (k-o3dims2%k_start)*row_length*rows
        ozone3D(i,j,k)=ozone_in(ijk)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
  !
CASE (io3_2dspec)
  !
  !       In this case we should have a 2D ancillary file so
  !       we return an error message if not.
  IF (.NOT. lexpand_ozone) THEN
    ErrorStatus=132

    CALL Ereport("O3_to_3D", ErrorStatus,                         &
       "A 3D ozone ancillary has been specified with a" //        &
       "2D ozone option.")
  END IF
  !
  !       The ozone is copied round a latitude circle. This does not
  !       conserve the total amount of ozone in the atmosphere because
  !       it does not allow for orography.
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,jk)                 &
!$OMP SHARED(ozone_levels,rows,row_length,o3dims2,ozone3D,ozone_in)
  DO k = 1, ozone_levels
    DO j = 1, rows
      DO i = 1, row_length
        jk = j + (k-o3dims2%k_start)*rows
        ozone3D(i,j,k)=ozone_in(jk)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
CASE (io3_2dmasscon)
  !
  ! ! !   Since z_top_of_model * eta levels doesn't work in new dynamics
  ! ! !   This option is not possible and this option will be identical
  ! ! !   with the option above.
  !
  !       The ozone field is expanded to three dimensions in such a way
  !       that mass-loading of each atmospheric layer is the same at
  !       each grid-point on a latitude circle, regardless of the
  !       orography.
  !
  !       The ozone levels do not start at the bottom of the model.

  k_off=model_levels-ozone_levels
  !
  !       In this provisional code the difference in density between
  !       the ancillary profile and the profile at the grid-point at
  !       the same level is ignored.
  !
  DO k = 1, ozone_levels
    DO j = 1, rows
      DO i = 1, row_length
        jk = j + (k-o3dims2%k_start)*rows
        ozone3D(i,j,k)=ozone_in(jk)
        !     &          * z_top_of_model
        !     &          * (eta_theta_levels(k+k_off)
        !     &           - eta_theta_levels(k-1+k_off))
        !     &          / (r_theta_levels(i, j, k+k_off)
        !     &           - r_theta_levels(i, j, k-1+k_off))

      END DO
    END DO
  END DO
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
CASE (io3_trop_map)
  !
  !       Still unoptimizied code: this option has been
  !       developed, but might not be the definite solution.
  !
  DO j=1, rows
    !

    !          ! Initialize the pointer to theta level just below the ozone
    !          ! tropopause level.
    k_o3_trop(:) = 0
    k_o3_trop_done(:) = .FALSE.
    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
    DO k=min_trop_level, max_trop_level
      !           ! Compute the ozone gradient in the ancillary file.
      !           ! In order to access data on row 4, you need to index
      !           ! j + 3 i.e. j + (k-1)*rows, hence the strange indexing
      jk   = j + (k-o3dims2%k_start)*rows
      jkp1 = j + (k-o3dims2%k_start+1)  *rows

      DO i=1, row_length
        !           ! min_trop_level should be a rho-level, and rho levels
        !           ! with the same index as theta-levels are higher up.

        !           ! Set the pointer to the current rho-level. When the
        !           ! tropopause criteria are met, this will be saved and will
        !           ! indicate the rho-level, between two ozone levels that
        !           ! contains the tropopause.
        IF (.NOT. k_o3_trop_done(i)) THEN
          k_o3_trop(i) = k
          !
          o3_gradient_anc(i,j,k) =                                &
              ( ozone_in(jkp1) - ozone_in(jk) ) /                 &
              ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
          !
          !           ! test whether we have found the ozone tropopause !
          IF ( o3_gradient_anc(i,j,k) > O3_grad_crit .AND.        &
             ozone_in(jk)           > O3_conc_crit .AND.          &
             ozone_in(jkp1)         > O3_strat_crit     )         &
             k_o3_trop_done(i) = .TRUE.
        END IF
      END DO
    END DO


    !          ! Calculate the height of the ozone tropopause above the
    !          ! surface. This is assumed to lie on a rho-level, and
    !          ! is used for the scaling of the ozone concentration
    !          ! in the model.
    !
    DO i=1, row_length
      z_trop_O3(i) = r_rho_levels(i, j, k_o3_trop(i))             &
                - r_theta_levels(i, j, 0)

      !
      IF (L_O3_trop_height) THEN
        O3_trop_height(i,j)=z_trop_O3(i)
      END IF

      IF (L_O3_trop_level) THEN
        O3_trop_level(i,j)=REAL(k_o3_trop(i))
      END IF
      !
      !           Set the thermal tropopause at the grid-point. Note that as
      !           TRINDX is counted upward from the surface, while the first
      !           rho level is omitted from the physics, no offset to the
      !           final index of r_rho_levels is required.
      z_trop_pt(i) = r_rho_levels(i, j, trindx(i, j))             &
                - r_theta_levels(i, j, 0)
      !

      IF (L_T_trop_height) THEN
        T_trop_height(i,j)=z_trop_pt(i)
      END IF

      IF (L_T_trop_level) THEN
        T_trop_level(i,j)=REAL(trindx(i,j))
      END IF
    END DO
    !
    !           Initialize the pointer that points to the lower boundary of
    !           the layer containing the height we are currently inter-
    !           polating ozone concentrations for.
    !          k_anc(:) =0
    !
    DO k=1, model_levels
      DO i=1, row_length
        !
        !           ! Find the height above the surface to the center of the
        !           ! current model layer. I.E. find the height to the theta
        !           ! level we want ozone on.
        height_above_surf                                         &
          = r_theta_levels(i, j, k) - r_theta_levels(i, j, 0)
        !
        IF (height_above_surf  <   z_trop_pt(i)) THEN
          !               We are in the thermal troposphere.
          !               Calculate the corresponding height in the ancillary
          !               file using similarity.
          z_int(i) = height_above_surf                            &
            * z_trop_O3(i) / z_trop_pt(i)
        ELSE
          !               We are in the thermal stratosphere.
          z_int(i) = z_trop_O3(i)                                 &
            + (height_above_surf - z_trop_pt(i))                  &
            * (z_top_of_model - z_trop_O3(i))                     &
            / (z_top_of_model - z_trop_pt(i))
        END IF
      END DO
      !
      !           ! Increment k_anc until it points to the theta level in
      !           ! the ancillary profile just below the interpolated height.
      !
      !          Do i=1, row_length
      !              Do While ( (r_theta_levels(i,j,k_anc(i))
      !     &                  - r_theta_levels(i, j, 0)  )   <
      !     &                   z_int(i) .AND. k_anc(i)  <=  model_levels )
      !                k_anc(i) = k_anc(i) + 1
      !              Enddo
      !!             In order to point to theta level below int. height.
      !              k_anc(i) = k_anc(i) - 1
      !            EndDo

      k_anc(:)=-1
      DO kl=0, model_levels
        DO i=1, row_length
          IF ( (                                                  &
              ( r_theta_levels(i,j,kl)-r_theta_levels(i, j, 0))   &
               >=   z_int(i) )                                    &
              .AND. ( k_anc(i)  ==  -1 ) ) THEN
            k_anc(i)=kl-1
          END IF
        END DO
      END DO

      WHERE (k_anc==-1) k_anc=model_levels

      DO i=1, row_length
        !
        !             Interpolate within this layer.
        jk=j+(k_anc(i)-o3dims2%k_start)*rows

        IF ( k_anc(i) == 0 ) THEN
          ozone3D(i, j, k) = ozone_in( j)
        ELSE IF ( k_anc(i) == model_levels ) THEN
          ozone3D(i, j, k) = ozone_in(jk)
        ELSE
          ozone3D(i, j, k) = ozone_in(jk)                           &
            + ( ozone_in(jk+rows) - ozone_in(jk) )                  &
            *  (z_int(i) - ( r_theta_levels(i,j,k_anc(i))           &
                        - r_theta_levels(i,j,0)     ) )             &
            / (r_theta_levels(i,j,k_anc(i)+1)                       &
                 -  r_theta_levels(i,j,k_anc(i)) )
        END IF
      END DO
      !
    END DO
    !
  END DO
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
CASE (io3_trop_map_masscon)
  !
  !       In this case we should have a 3D ancillary file so
  !       we return an error message if not.
  IF (.NOT. lexpand_ozone) THEN
    ErrorStatus=123

    CALL Ereport("O3_to_3D", ErrorStatus,                         &
       "A 3D ozone ancillary has been specified with a" //        &
       "2D ozone option.")
  END IF
  !
  !       Provisional unoptimizied code: these options are still being
  !       developed.
  !
  !       Interpolation using cumulative ozone amounts is preferred to
  !       conserved column masses.
  !
  DO j=1, rows
    !
    !         Calculate cumulative amounts of ozone to the top of each
    !         layer in the field supplied (i.e. rho level).
    o3_mass_anc(0)=0.0
    !
    !         When ozone_levels is less than model_levels,
    !         the lowest levels all have the same ozone mass mixing ratio.
    !

    !
    DO k=1, model_levels-ozone_levels
      o3_mass_anc(k) = o3_mass_anc(k-1)                           &
        + ozone_in(j) * avg_layer_mass(j,k)
    END DO
    DO k=model_levels-ozone_levels+1, model_levels
      jk=j+(k-o3dims2%k_start)*rows
      o3_mass_anc(k) = o3_mass_anc(k-1)                           &
        + ozone_in(jk) * avg_layer_mass(j,k)

    END DO
    !
    !
    !
    !          ! Initialize the pointer to theta level just below the ozone
    !          ! tropopause level.
    k_o3_trop(:) = 0
    k_o3_trop_done(:) = .FALSE.
    !
    DO k=min_trop_level+1, max_trop_level

      !           ! Compute the ozone gradient in the ancillary file.
      !           ! In order to access data on row 4, you need to index
      !           ! j + 3 i.e. j + (k-1)*rows, hence the strange indexing

      jk   = j + (k-o3dims2%k_start)*rows
      jkp1 = j + (k-o3dims2%k_start+1)  *rows

      !cdir novector
      DO i=1, row_length

        !           ! min_trop_level should be a rho-level, and rho levels
        !           ! with the same index as theta-levels are higher up.
        !           !
        !           ! Set the pointer to the current rho_level. When the
        !           ! tropopause criteria are met, this will be saved and will
        !           ! indicate the rho-level, between two ozone levels that
        !           ! contains the tropopause.

        IF (.NOT. k_o3_trop_done(i)) THEN
          k_o3_trop(i) = k
          !
          o3_gradient_anc(i,j,k) =                                &
              ( ozone_in(jkp1) - ozone_in(jk) ) /                 &
              ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
          !
          !           ! test whether we have found the ozone tropopause !
          IF (o3_gradient_anc(i,j,k) > O3_grad_crit .AND.         &
            ozone_in(jk)           > O3_conc_crit .AND.           &
            ozone_in(jkp1)         > O3_strat_crit     )          &
            k_o3_trop_done(i) = .TRUE.
        END IF
      END DO
    END DO

    DO i=1, row_length

      !          ! Calculate the height of the ozone tropopause above the
      !          ! surface. This is assumed to lie on a rho-level, and
      !          ! is used for the scaling of the ozone concentration
      !          ! in the model.
      !
      z_trop_O3(i) = r_rho_levels(i, j, k_o3_trop(i))             &
                - r_theta_levels(i, j, 0)

      IF (L_O3_trop_height) THEN
        O3_trop_height(i,j)=z_trop_O3(i)
      END IF

      IF (L_O3_trop_level) THEN
        O3_trop_level(i,j)=REAL(k_o3_trop(i))
      END IF

      !           Set the thermal tropopause at the grid-point. Note that as
      !           TRINDX is counted upward from the surface, while the first
      !           rho level is omitted from the physics, no offset to the
      !           final index of r_rho_levels is required.
      z_trop_pt(i) = r_rho_levels(i, j, trindx(i, j))             &
                - r_theta_levels(i, j, 0)
      !
      !
      IF (L_T_trop_height) THEN
        T_trop_height(i,j)=z_trop_pt(i)
      END IF

      IF (L_T_trop_level) THEN
        T_trop_level(i,j)=REAL(trindx(i,j))
      END IF
    END DO

    !           Initialize the pointer, that points to the lower boundary
    !           of the layer containing the height we are currently inter-
    !           polating ozone concentrations for.
    !
    !          k_anc(:)=1
    !
    !           At the surface there is no ozone below the current level.
    o3_mass_cumul_below(:)=0.0
    !
    DO k=1, model_levels

      IF ((k >  1) .AND. (k <  model_levels)) THEN

        DO i=1, row_length
          !
          !           ! Find the height above the surface to the center of the
          !           ! current model layer. I.E. find the height to the theta
          !           ! level we want ozone on.

          height_above_surf                                         &
            = r_rho_levels(i, j, k+1) - r_theta_levels(i, j, 0)
          IF (height_above_surf  <=  z_trop_pt(i)) THEN
            !           !   We are in the thermal troposphere.
            !           !   Calculate the corresponding height in the ancillary
            !           !   file using similarity.
            z_int(i) = height_above_surf                            &
              * z_trop_O3(i) / z_trop_pt(i)
          ELSE
            !           !   We are in the thermal stratosphere.
            z_int(i) = z_trop_O3(i)                                 &
              + (height_above_surf - z_trop_pt(i))                  &
              * (z_top_of_model - z_trop_O3(i))                     &
              / (z_top_of_model - z_trop_pt(i))
          END IF
        END DO
        !
        !           ! Increment k_anc until it points to the theta level in
        !           ! the ancillary profile just below the interpolated height.
        !
        !            Do i=1, row_length
        !
        !              Do While ( ( r_theta_levels(i,j,k_anc(i))
        !     &                   - r_theta_levels(i,j,0)   )  <
        !     &                   z_int(i) )
        !                if(k_anc(i) <  model_levels) then
        !                  k_anc(i) = k_anc(i) + 1
        !                else
        !                  exit
        !                endif
        !              Enddo
        !!!             In order to point to theta level below int. height.
        !              k_anc(i) = k_anc(i) - 1
        !            EndDo


        k_anc(:)=-1
        DO kl=0, model_levels
          DO i=1, row_length
            IF ( (                                                  &
                ( r_theta_levels(i,j,kl)-r_theta_levels(i, j, 0))   &
                 >=   z_int(i) )                                    &
                .AND. ( k_anc(i)  ==  -1 ) ) THEN
              k_anc(i)=kl-1
            END IF
          END DO
        END DO

        WHERE (k_anc==-1) k_anc=model_levels




      END IF
      !
      !             Interpolate within this layer.

      DO i=1, row_length

        IF ((k >  1) .AND. (k <  model_levels)) THEN

          !             Calculate the o3_mass_cumul at z_int using the o3_mass_anc
          !             at the theta levels below and above z_int
          !
          o3_mass_cumul(i) = o3_mass_anc(k_anc(i)-1)              &
          + ( o3_mass_anc(k_anc(i)) - o3_mass_anc(k_anc(i)-1) )   &
          *   (z_int(i) - ( r_rho_levels(i,j,k_anc(i))            &
                       - r_theta_levels(i,j,0)   ) )              &
            / ( r_rho_levels(i,j,k_anc(i)+1)                      &
              - r_rho_levels(i,j,k_anc(i))   )
          !
        END IF
        IF (k == 1) THEN
          o3_mass_cumul(i)=o3_mass_anc(1)
        END IF
        IF (k == model_levels) THEN
          o3_mass_cumul(i)=o3_mass_anc(model_levels)
        END IF
        !
        !           ! Calculate the ozone contentration on the theta level
        !           ! at the model gridpoint, using
        !           ! qO3 = delta_QO3 / (delta_rho * delta_z)
        !
        IF (k == 1) THEN
          ozone3D(i, j, k) =                                      &
                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
                  1.0 /( rho(i, j, k+1)                            &
                  * (r_rho_levels(i, j, k+1)                      &
                  - r_theta_levels(i, j, k-1) ) )
        END IF
        IF ((k >  1) .AND. (k <  model_levels)) THEN
          ozone3D(i, j, k) =                                      &
                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
                  2.0 / ( (rho(i, j, k) + rho(i, j, k+1) )         &
                  * (r_rho_levels(i, j, k+1)                      &
                  - r_rho_levels(i, j, k) ) )
          !
        END IF
        IF (k == model_levels) THEN
          ozone3D(i, j, k) =                                      &
                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
                  1.0/( rho(i, j, k)                              &
                  * (r_theta_levels(i ,j ,k)                      &
                  - r_rho_levels(i, j, k) ) )
        END IF

        o3_mass_cumul_below(i)=o3_mass_cumul(i)
        !
      END DO
      !
    END DO
    !
  END DO
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
CASE DEFAULT
  !
  cmessage='*** Error: Unrecognized expansion of ozone.'
  ErrorStatus=123
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
  !
END SELECT
!

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE O3_to_3D
