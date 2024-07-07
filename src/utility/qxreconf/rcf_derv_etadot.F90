! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_derv_etadot_mod

IMPLICIT NONE
!  Subroutine Rcf_Derv_Etadot_Mod

! Description:
!
!   Derive Etadot as previously done at runtime.
!
! Method:
!

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_ETADOT_MOD'

CONTAINS

SUBROUTINE rcf_derv_etadot( fields_in, field_count_in,    &
                              fields_out, field_count_out,  &
                              hdr_in,hdr_out,               &
                              etadot )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_mv,  stashcode_theta, stashcode_prog_sec, &
    stashcode_orog, stashcode_thetavd, stashcode_exner, &
    stashcode_u,stashcode_v,stashcode_w,stashcode_dry_rho

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parvars, ONLY: &
    gc_all_proc_group,g_datastart

USE um_parcore, ONLY: &
    mype,nproc

USE rcf_gather_field_mod, ONLY: &
    rcf_gather_field_real
USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE decomp_params, ONLY: &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE rcf_interp_weights_mod, ONLY:               &
  h_int_active

USE rcf_v_int_ctl_mod, ONLY:                    &
  v_int_active

USE rcf_set_interp_flags_mod, ONLY:             &
  interp_h_only,                               &
  interp_v_only,                               &
  interp_copy,                                 &
  interp_all,                                  &
  interp_no_op

USE rcf_grid_type_mod, ONLY:                    &
    grid_type,                                   &
    input_grid, output_grid

USE rcf_interpolate_mod, ONLY:                  &
    rcf_interpolate

USE rcf_field_equals_mod, ONLY:                 &
    rcf_field_equals

USE rcf_headaddress_mod, ONLY: &
    fh_gridstagger_endgame

USE planet_constants_mod, ONLY: planet_radius

USE rcf_generate_heights_mod, ONLY: rcf_generate_heights
USE latlon_eq_rotation_mod, ONLY: rotate_eq_to_latlon
USE set_metric_terms_4a_mod, ONLY: eg_phi, tny

USE rcf_headaddress_mod, ONLY:                    &
    rc_polelong,                  rc_longspacing,  &
    rc_latspacing,                rc_polelat,      &
    rc_firstlong,                 rc_firstlat,     &
    rc_modeltop,                  rc_swldeg,       &
    rc_wedgedeg

USE ereport_mod, ONLY: &
    ereport

USE conversions_mod, ONLY: pi

USE cppxref_mod, ONLY: ppx_atm_tall,    &
                       ppx_atm_cvall,   &
                       ppx_atm_cuall

USE stparam_mod, ONLY: &
    st_levels_model_theta, st_levels_model_rho

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( field_type ), POINTER       :: fields_in(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( um_header_type), INTENT(IN) :: hdr_in
TYPE( field_type ), INTENT(INOUT), TARGET :: etadot
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_count_in

! Internal variables
TYPE( field_type ), POINTER       :: u
TYPE( field_type ), POINTER       :: v
TYPE( field_type ), POINTER       :: w
TYPE( field_type ), POINTER       :: orog
TYPE( field_type )                :: urho
TYPE( field_type )                :: vrho
TYPE( field_type )                :: u_at_w
TYPE( field_type )                :: v_at_w

CHARACTER (LEN=errormessagelength):: cmessage       ! used for EReport
CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_DERV_ETADOT'
INTEGER                           :: errorstatus

INTEGER                           :: interp_option ! interpolation to perform


INTEGER                           :: pos      ! position in array
INTEGER                           :: i,j,k,ii,jj ! loop index

REAL :: intw_w2rho

REAL                               :: xi2_p                            &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)
REAL                               :: dxi2_v                           &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)
REAL                               :: xi1_p                            &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)
REAL                               :: dxi1_u                           &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)

REAL, ALLOCATABLE :: u_level(:,:)
REAL, ALLOCATABLE :: v_level(:,:)
REAL, ALLOCATABLE :: p_level(:,:)
REAL, ALLOCATABLE :: du_level(:,:)
REAL, ALLOCATABLE :: dv_level(:,:)
REAL, ALLOCATABLE :: up_level(:,:)
REAL, ALLOCATABLE :: vp_level(:,:)
REAL, ALLOCATABLE :: dxi3_at_u_w(:,:)
REAL, ALLOCATABLE :: dxi3_at_v_w(:,:)
REAL, ALLOCATABLE :: xi1_p_global(:,:)
REAL, ALLOCATABLE :: xi2_p_global(:,:)
REAL, ALLOCATABLE :: xi1_u_global(:,:)
REAL, ALLOCATABLE :: xi2_v_global(:,:)
REAL, ALLOCATABLE :: dxi1_u_global(:,:)
REAL, ALLOCATABLE :: dxi2_v_global(:,:)
REAL, ALLOCATABLE :: xi3_at_rho(:,:,:)
REAL, ALLOCATABLE :: xi3_at_theta(:,:,:)
REAL, ALLOCATABLE :: xi3_at_u(:,:,:)
REAL, ALLOCATABLE :: xi3_at_v(:,:,:)
REAL, ALLOCATABLE :: xi3_at_u_w(:,:,:)
REAL, ALLOCATABLE :: xi3_at_v_w(:,:,:)


INTEGER :: div,jend,pe,ij,blk

REAL :: phi,rdxi2,rdxi1,rdeta,h1_p_eta,h2_p_eta,&
     h3_p_eta,deta_xi3_theta,dxi1_xi3,dxi2_xi3
REAL :: intw_rho2w

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 'Deriving etadot '
  CALL umPrint(umMessage,src='rcf_derv_etadot')
  WRITE (umMessage,FMT='(2a)') 'NOTE: assuming global, spherical, deep,',&
                       ' eccentricity==0!'
  CALL umPrint(umMessage,src='rcf_derv_etadot')
END IF

CALL rcf_locate(stashcode_prog_sec, stashcode_u, &
                  fields_out, field_count_out, pos)
u => fields_out(pos)
CALL rcf_alloc_field( u )
CALL rcf_read_field(  u, hdr_out, decomp_rcf_output )

CALL rcf_locate(stashcode_prog_sec, stashcode_v, &
                  fields_out, field_count_out, pos)
v => fields_out(pos)
CALL rcf_alloc_field( v )
CALL rcf_read_field( v, hdr_out, decomp_rcf_output )

CALL rcf_locate(stashcode_prog_sec, stashcode_w, &
                  fields_out, field_count_out, pos)
w => fields_out(pos)
CALL rcf_alloc_field( w )
CALL rcf_read_field( w, hdr_out, decomp_rcf_output )

CALL rcf_field_equals( u_at_w, fields_out(pos) )
CALL rcf_alloc_field( u_at_w )
CALL rcf_field_equals( v_at_w, fields_out(pos) )
CALL rcf_alloc_field( v_at_w )


CALL rcf_locate(stashcode_prog_sec, stashcode_dry_rho, &
                  fields_out, field_count_out, pos)
CALL rcf_field_equals( urho, fields_out(pos) )
CALL rcf_alloc_field( urho )
CALL rcf_field_equals( vrho, fields_out(pos) )
CALL rcf_alloc_field( vrho )



SELECT CASE(output_grid % grid_stagger)
CASE (fh_gridstagger_endgame)

CASE DEFAULT
  ! Print error if grid staggering is not catered for in this routine.
  cmessage = 'Grid staggering method is not supported.'
  errorstatus= 12
  CALL ereport( routinename, errorstatus, cmessage )
END SELECT


ALLOCATE( u_level(  output_grid % glob_u_row_length, &
                    output_grid % glob_u_rows ) )
ALLOCATE( v_level(  output_grid % glob_v_row_length, &
                    output_grid % glob_v_rows ) )
ALLOCATE( p_level(  output_grid % glob_p_row_length, &
                    output_grid % glob_p_rows ) )
ALLOCATE( up_level( output_grid % glob_p_row_length,&
                    output_grid % glob_p_rows ) )
ALLOCATE( vp_level( output_grid % glob_p_row_length,&
                    output_grid % glob_p_rows ) )


! So do Gather 1 level per PE and then Scatter again.
div = output_grid % model_levels / nproc

IF ( nproc * div  <   output_grid % model_levels ) div = div + 1

pe = 0

ALLOCATE(  xi1_p_global( output_grid % glob_p_row_length, &
                         output_grid % glob_p_rows ) )
ALLOCATE(  xi2_p_global( output_grid % glob_p_row_length, &
                         output_grid % glob_p_rows ) )
ALLOCATE(  xi1_u_global( output_grid % glob_u_row_length, &
                         output_grid % glob_u_rows ) )
ALLOCATE(  xi2_v_global( output_grid % glob_v_row_length, &
                         output_grid % glob_v_rows ) )
ALLOCATE( dxi1_u_global( output_grid % glob_p_row_length, &
                         output_grid % glob_p_rows ) )
ALLOCATE( dxi2_v_global( output_grid % glob_p_row_length, &
                         output_grid % glob_p_rows ) )

divdo : DO blk = 1, div

  IF ( blk == div ) THEN
    jend = output_grid % model_levels
  ELSE
    jend =  blk * nproc
  END IF

  DO j = ((blk-1) * nproc) + 1, jend
    ! Will gather level j on PE pe
    CALL rcf_gather_field_real( u % DATA( :, j), u_level,             &
                                u % row_len,                          &
                                u % rows,                             &
                                u % glob_row_len,                     &
                                u % glob_rows, pe,                    &
                                gc_all_proc_group )

    CALL rcf_gather_field_real( v % DATA( :, j), v_level,             &
                                v % row_len,                          &
                                v % rows,                             &
                                v % glob_row_len,                     &
                                v % glob_rows, pe,                    &
                                gc_all_proc_group )

    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO


  !--------------------------------------------------------------------
  ! Set up a latitude/xi1_p field for the grid.
  !--------------------------------------------------------------------

  DO j = 1, output_grid % loc_p_rows
    DO i = 1, output_grid % loc_p_row_length

      xi2_p (i,j) = hdr_out % realc(rc_firstlat) +  &
                  (g_datastart(2,mype)+j-1.5) *  &
                   hdr_out % realc(rc_latspacing)

      xi1_p(i,j) = hdr_out % realc(rc_firstlong) +  &
                  (g_datastart(1,mype)+i-1.5)  *  &
                   hdr_out % realc(rc_longspacing)

    END DO
  END DO

  !--------------------------------------------------------------------
  ! For rotated grids, get true lats and longs.
  !--------------------------------------------------------------------

  IF (output_grid % rotated) THEN

    ! Use the same arrays to store true lats & longs

    CALL rotate_eq_to_latlon(xi2_p, xi1_p, xi2_p, xi1_p,  &
                hdr_out % realc(rc_polelat), hdr_out % realc(rc_polelong), &
                output_grid % loc_p_field)

  END IF


  xi1_p(:,:)=xi1_p(:,:)/180.0*pi
  xi2_p(:,:)=xi2_p(:,:)/180.0*pi

  ! now compute dxi1_u and dxi2_v

  pe=0
  DO j = ((blk-1) * nproc) + 1, jend
    ! we need this on every cpu!
    CALL rcf_gather_field_real( xi1_p, xi1_p_global, &
                              output_grid % loc_p_row_length,      &
                              output_grid % loc_p_rows,            &
                              output_grid % glob_p_row_length,     &
                              output_grid % glob_p_rows, pe,       &
                             gc_all_proc_group )

    CALL rcf_gather_field_real( xi2_p, xi2_p_global, &
                              output_grid % loc_p_row_length,      &
                              output_grid % loc_p_rows,            &
                              output_grid % glob_p_row_length,     &
                              output_grid % glob_p_rows, pe,       &
                             gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO

  ! assuming regular global grid for now ....


  DO i=2,output_grid % glob_u_row_length
    DO j=1,output_grid % glob_u_rows

      xi1_u_global(i,j) = 0.5*(xi1_p_global(i-1,j)+xi1_p_global(i,j))

    END DO
  END DO

  IF (output_grid % global) THEN
    i=1
    DO j=1,output_grid % glob_u_rows
      xi1_u_global(i,j) = 0.5*(xi1_p_global(i,j)+&
                               xi1_p_global(output_grid%glob_p_row_length,j) &
                               -2.0*pi)
    END DO
  ELSE
    i=1
    DO j=1,output_grid % glob_u_rows
      xi1_u_global(i,j) = xi1_p_global(i,j)
    END DO
  END IF

  DO i=1,output_grid % glob_v_row_length
    DO j=2,output_grid % glob_v_rows-1

      xi2_v_global(i,j) = 0.5*(xi2_p_global(i,j-1)+xi2_p_global(i,j))

    END DO
  END DO


  IF (output_grid % global) THEN
    ! assume that this means global and LAM otherwise
    DO i=1,output_grid % glob_v_row_length
      j=1
      xi2_v_global(i,j) = -pi/2.0

      j=output_grid % glob_v_rows
      xi2_v_global(i,j) = pi/2.0
    END DO
  ELSE
    DO i=1,output_grid % glob_v_row_length
      j=1
      xi2_v_global(i,j) = xi2_v_global(i,2)

      j=output_grid % glob_v_rows
      xi2_v_global(i,j) = xi2_v_global(i,j-1)
    END DO
  END IF

  ! Now compute the difference:

  DO i=1,output_grid % glob_p_row_length-1
    DO j=1,output_grid % glob_p_rows
      dxi1_u_global(i,j) = xi1_u_global(i+1,j)-xi1_u_global(i,j)
    END DO
  END DO


  IF (output_grid % global) THEN
    ! assume that this means global and LAM otherwise
    DO j=1,output_grid % glob_p_rows
      dxi1_u_global(output_grid % glob_p_row_length,j) = &
             (xi1_u_global(1,j)+2.0*pi)-&
              xi1_u_global(output_grid % glob_u_row_length,j)
    END DO
  ELSE
    DO j=1,output_grid % glob_p_rows
      dxi1_u_global(output_grid % glob_p_row_length,j) = &
      dxi1_u_global(output_grid % glob_p_row_length-1,j)
    END DO
  END IF

  DO i=1,output_grid % glob_p_row_length
    DO j=1,output_grid % glob_p_rows
      dxi2_v_global(i,j) = xi2_v_global(i,j+1)-xi2_v_global(i,j)
    END DO
  END DO

  pe=0
  DO j = ((blk-1) * nproc) + 1, jend
    ! we need this on every cpu - this is not the nicest
    ! way of doing this, but works for now ...
    CALL rcf_scatter_field_real( dxi1_u, dxi1_u_global, &
                              output_grid % loc_p_row_length,      &
                              output_grid % loc_p_rows,            &
                              output_grid % glob_p_row_length,     &
                              output_grid % glob_p_rows, pe,       &
                             gc_all_proc_group )

    CALL rcf_scatter_field_real( dxi2_v, dxi2_v_global, &
                              output_grid % loc_p_row_length,      &
                              output_grid % loc_p_rows,            &
                              output_grid % glob_p_row_length,     &
                              output_grid % glob_p_rows, pe,       &
                             gc_all_proc_group )
    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO

  ! assuming regular global grid for now ....

  DO i=1,output_grid % glob_p_row_length-1
    DO j=1,output_grid % glob_p_rows

      up_level(i,j) = 0.5*(u_level(i,j)+u_level(i+1,j))
      vp_level(i,j) = 0.5*(v_level(i,j)+v_level(i,j+1))

    END DO
  END DO

  i=output_grid % glob_p_row_length
  DO j=1,output_grid % glob_p_rows
    up_level(i,j) = 0.5*(u_level(i,j)+u_level(1,j))
    vp_level(i,j) = 0.5*(v_level(i,j)+v_level(i,j+1))
  END DO


  pe = 0

  DO j = ((blk-1) * nproc) + 1, jend

    CALL rcf_scatter_field_real( urho%DATA( :, j), up_level, &
                              urho % row_len,                &
                              urho % rows,                   &
                              urho % glob_row_len,           &
                              urho % glob_rows, pe,          &
                              gc_all_proc_group )

    CALL rcf_scatter_field_real( vrho%DATA( :, j), vp_level, &
                              vrho % row_len,                &
                              vrho % rows,                   &
                              vrho % glob_row_len,           &
                              vrho % glob_rows, pe,          &
                             gc_all_proc_group )

    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO
END DO divdo

DEALLOCATE( dxi2_v_global )
DEALLOCATE( dxi1_u_global )


DEALLOCATE( xi2_v_global )
DEALLOCATE( xi1_u_global )
DEALLOCATE( xi2_p_global )
DEALLOCATE( xi1_p_global )


ALLOCATE( dxi3_at_u_w (  output_grid % glob_u_row_length * &
                         output_grid % glob_u_rows,        &
                         output_grid % model_levels) )
ALLOCATE( dxi3_at_v_w (  output_grid % glob_v_row_length * &
                         output_grid % glob_v_rows,        &
                         output_grid % model_levels) )

DEALLOCATE( vp_level )
DEALLOCATE( up_level )


! orogoraphy in output dump, must exist...
CALL rcf_locate( stashcode_prog_sec, stashcode_orog,       &
               fields_out,field_count_out,pos)
orog => fields_out(pos)
CALL rcf_alloc_field( orog )
CALL rcf_read_field( orog, hdr_out, decomp_rcf_output )


ALLOCATE (xi3_at_rho (output_grid % loc_p_row_length ,  &
                      output_grid % loc_p_rows,         &
                      0:output_grid % model_levels+1))
ALLOCATE (xi3_at_theta (output_grid % loc_p_row_length ,&
                        output_grid % loc_p_rows,       &
                        0:output_grid % model_levels+1))
ALLOCATE (xi3_at_u (output_grid % loc_u_row_length ,    &
                    output_grid % loc_u_rows,           &
                    0:output_grid % model_levels+1))
ALLOCATE (xi3_at_v (output_grid % loc_v_row_length ,    &
                    output_grid % loc_v_rows,           &
                    0:output_grid % model_levels+1))
ALLOCATE (xi3_at_u_w (output_grid % loc_u_row_length ,  &
                      output_grid % loc_u_rows,         &
                      0:output_grid % model_levels+1))
ALLOCATE (xi3_at_v_w (output_grid % loc_v_row_length ,  &
                      output_grid % loc_v_rows,         &
                      0:output_grid % model_levels+1))


! from rcf_post_process_mod.F90
CALL rcf_generate_heights( output_grid,orog,                       &
                           ppx_atm_tall, st_levels_model_theta,    &
                           xi3_at_theta , etadot % level_size)

CALL rcf_generate_heights( output_grid,orog,                       &
                           ppx_atm_tall, st_levels_model_rho,      &
                           xi3_at_rho , etadot % level_size )

CALL rcf_generate_heights( output_grid,orog,                       &
                           ppx_atm_cvall, st_levels_model_rho,     &
                           xi3_at_v ,  v % level_size )

CALL rcf_generate_heights( output_grid,orog,                       &
                           ppx_atm_cuall, st_levels_model_rho,     &
                           xi3_at_u ,  u % level_size )

CALL rcf_dealloc_field( orog )

u_at_w%DATA(:,:) = 0.0
v_at_w%DATA(:,:) = 0.0


DO k=2,output_grid % model_levels

  intw_rho2w = ( output_grid % eta_theta_levels(k-1)      &
                -output_grid % eta_rho_levels(k-1) ) /    &
               ( output_grid % eta_rho_levels(k)          &
                -output_grid % eta_rho_levels(k-1) )

  u_at_w%DATA(:, k )  = (1.0-intw_rho2w)*urho%DATA(:,k-1)+ &
                            intw_rho2w *urho%DATA(:,k)

  v_at_w%DATA(:, k )  = (1.0-intw_rho2w)*vrho%DATA(:,k-1)+ &
                            intw_rho2w *vrho%DATA(:,k)

END DO

ALLOCATE( dv_level( output_grid % glob_p_row_length, &
                    output_grid % glob_p_rows ) )

ALLOCATE( du_level( output_grid % glob_p_row_length, &
                    output_grid % glob_p_rows ) )


divdo3 : DO blk = 1, div

  IF (blk == div ) THEN
    jend = output_grid % model_levels
  ELSE
    jend =  blk * nproc
  END IF

  pe = 0

  DO j = ((blk-1) * nproc) + 1, jend
    ! Will gather level j on PE pe
    CALL rcf_gather_field_real( xi3_at_theta( :,:, j), p_level,      &
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )


    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO

  ! compute a r_at_v_w as in ./control/top_level/setcona_4A.F90

  ! assuming regular global grid for now ....

  DO j=1,output_grid % glob_u_rows
    DO i=2,output_grid % glob_u_row_length

      u_level(i,j) = 0.5*(p_level(i-1,j)+p_level(i,j))

    END DO
    IF (output_grid % global) THEN
      i = output_grid % glob_p_row_length
      u_level(1,j) = 0.5*(p_level(1,j)+p_level(i,j))
    ELSE
      u_level(1,j) = p_level(1,j)
    END IF
  END DO

  DO i=1,output_grid % glob_v_row_length
    DO j=2,output_grid % glob_v_rows-1
      v_level(i,j) = 0.5*(p_level(i,j-1)+p_level(i,j))
    END DO
    IF (output_grid % global) THEN
      ii = i + output_grid % glob_p_row_length/2
      IF (ii  >   output_grid % glob_p_row_length) &
         ii = ii - output_grid % glob_p_row_length
      j=1
      v_level(i,j) = 0.5*(p_level(i,j)+p_level(ii,j))

      j=output_grid % glob_v_rows
      v_level(i,j) = 0.5*(p_level(ii,output_grid % glob_p_rows)&
                         +p_level(i,output_grid % glob_p_rows))
    ELSE
      j=1
      v_level(i,j) = p_level(i,j)

      j=output_grid % glob_v_rows
      v_level(i,j) = p_level(i,output_grid % glob_p_rows)
    END IF
  END DO

  ! now compute the difference .... if we just had halos ....

  ! assuming regular global grid for now ....

  DO j=1,output_grid % glob_p_rows
    DO i=1,output_grid % glob_p_row_length-1
      du_level(i,j) = (u_level(i+1,j)-u_level(i,j))
    END DO
    IF (output_grid % global) THEN
      i=output_grid % glob_p_row_length
      du_level(i,j) = (u_level(1,j)&
                      -u_level(output_grid % glob_u_row_length,j))
    ELSE
      i=output_grid % glob_p_row_length
      du_level(i,j) = u_level(i-1,j)
    END IF
  END DO

  DO i=1,output_grid % glob_p_row_length
    DO j=1,output_grid % glob_p_rows
      dv_level(i,j) = (v_level(i,j+1)-v_level(i,j))
    END DO
  END DO

  pe = 0

  DO j = ((blk-1) * nproc) + 1, jend

    CALL rcf_scatter_field_real(dxi3_at_u_w( :, j), du_level,        &
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )

    CALL rcf_scatter_field_real(dxi3_at_v_w( :, j), dv_level,        &
                                output_grid % loc_p_row_length,      &
                                output_grid % loc_p_rows,            &
                                output_grid % glob_p_row_length,     &
                                output_grid % glob_p_rows, pe,       &
                                gc_all_proc_group )

    pe = pe + 1
    IF (pe == nproc) pe = 0
  END DO
END DO divdo3

DEALLOCATE( dv_level )
DEALLOCATE( du_level )
DEALLOCATE(  v_level )
DEALLOCATE(  u_level )


! h1_p(i,j,k)     = 1.0
! h2_p(i,j,k)     = 1.0
! h3_p_eta = 1. if cartesian,

! h1_p_eta(i,j,k) = planet_radius * cos(xi2_p(j))
! h2_p_eta(i,j,k) = planet_radius
! h3_p_eta = 1. if shallow,

! eps = eccentricity   == 0 for all runs      ! Hardwired for spherical coords
! a   = planet_radius
! phi =  eg_phi(xi2_p(j),xi3_at_theta(i,j,k),a,eps)
! tmp = SQRT( 1.0 - (eps*SIN(phi))**2 )

! h1_p_eta(i,j,k) = xi3_at_theta(i,j,k)*COS(phi)/tmp
! h2_p_eta(i,j,k) = xi3_at_theta(i,j,k)*SIN(2.0*phi+tny)          &
!                          /(tmp*SIN(2.0*xi2_p(j)+tny))
! h3_p_eta(i,j,k) = tmp


DO k=2,output_grid % model_levels
  DO i=1, output_grid % loc_p_row_length
    DO j=1, output_grid % loc_p_rows

      ij = (j-1)*output_grid % loc_p_row_length + i

      phi =  eg_phi(xi2_p(i,j),xi3_at_theta(i,j,k-1),planet_radius,0.0)

      h1_p_eta = xi3_at_theta(i,j,k-1)*COS(phi)
      h2_p_eta = xi3_at_theta(i,j,k-1)*SIN(2.0*phi+tny)&
                 /(SIN(2.0*xi2_p(i,j)+tny))
      h3_p_eta = 1.0

      rdeta = 1.0/( output_grid % eta_rho_levels(k) - &
                    output_grid % eta_rho_levels(k-1) )
      deta_xi3_theta = ( xi3_at_rho(i,j,k) -  xi3_at_rho(i,j,k-1) )*rdeta

      rdxi2    = 1.0/(dxi2_v(i,j))
      rdxi1    = 1.0/(dxi1_u(i,j))

      ! computed in set_metric_terms_4A
      dxi1_xi3 = dxi3_at_u_w(ij,k-1) * rdxi1
      dxi2_xi3 = dxi3_at_v_w(ij,k-1) * rdxi2

      etadot%DATA(ij,k) = (    w%DATA(ij,k)/h3_p_eta  &
                           - u_at_w%DATA(ij,k)* dxi1_xi3/ h1_p_eta &
                           - v_at_w%DATA(ij,k)* dxi2_xi3/ h2_p_eta)&
                           / deta_xi3_theta

    END DO
  END DO
END DO

DEALLOCATE( xi3_at_v_w )
DEALLOCATE( xi3_at_u_w )
DEALLOCATE( xi3_at_v )
DEALLOCATE( xi3_at_u )
DEALLOCATE( xi3_at_theta )
DEALLOCATE( xi3_at_rho )

DEALLOCATE( dxi3_at_u_w )
DEALLOCATE( dxi3_at_v_w )


etadot%DATA(:,1) = 0.0
etadot%DATA(:,output_grid % model_levels+1) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_derv_etadot
END MODULE rcf_derv_etadot_mod
