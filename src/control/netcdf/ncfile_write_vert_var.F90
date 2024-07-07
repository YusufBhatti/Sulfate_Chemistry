! *****************************COPYRIGHT*******************************
! 
! Copyright 2017-2018 University of Reading
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation 
! and/or other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors 
! may be used to endorse or promote products derived from this software without 
! specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: NetCDF output
!
!  Purpose: Writes data for coordinate variables, auxiliary coordinate variables
!           and bounds variables created in ncfile_write_vert_dim


MODULE ncfile_write_vert_var_mod

USE netcdf,                 ONLY: nf90_max_name ! External netCDF library module
USE yomhook,                ONLY: lhook,dr_hook
USE parkind1,               ONLY: jprb,jpim
USE umprintmgr,             ONLY: umprint,ummessage,printstatus,prdiag
USE file_manager,           ONLY: um_file_type
USE cstash_mod,             ONLY: idom_b,dompro
USE ppxlook_mod,            ONLY: exppxi
USE submodel_mod,           ONLY: atmos_im
USE profilename_length_mod, ONLY: profilename_length
USE stextend_mod,           ONLY: levlst_s,rlevlst_s,llistty
USE version_mod,            ONLY: nlevp_s
USE stash_array_mod,        ONLY: stlist,stash_levels
USE nlsizes_namelist_mod,   ONLY: model_levels
USE dump_headers_mod,       ONLY: a_levdepc
USE cppxref_mod,            ONLY: &
    ppx_lv_code,ppx_pt_code
USE stparam_mod,            ONLY: &
    st_diag_address,st_output_bottom,st_output_top,st_gridpoint_code, &
    global_mean_base,vert_mean_base,st_sect_no_code,st_item_code, &
    block_size,st_input_bottom,st_input_top, st_levels_model_rho, &
    st_levels_model_theta, st_levels_single, st_levels_deep_soil
USE atm_d1_indices_mod,     ONLY: &
    jetatheta,jetarho,jzseak_theta,jzseak_rho,jck_theta,jck_rho,jsoil_thickness
USE nc_dimension_id_mod,    ONLY: &
    nc_vert_id_model_level, &
    nc_vert_id_eta_rho,nc_vert_id_zsea_rho,nc_vert_id_C_rho, &
    nc_vert_id_eta_theta,nc_vert_id_zsea_theta,nc_vert_id_C_theta
USE umnetcdf_mod,           ONLY: &
    nc_inq_ndims,nc_inq_dim,nc_get_varid,nc_get_att,nc_put_var
USE cf_metadata_mod,        ONLY: cf_extra_info

IMPLICIT NONE

! Dr hook
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NCFILE_WRITE_VERT_VAR_MOD'

CONTAINS

SUBROUTINE ncfile_write_vert_var(um_file)

IMPLICIT NONE

TYPE(um_file_type), INTENT(IN), POINTER :: um_file

!
!  Local variables
!
INTEGER :: ndims             ! Number of dimensions in NetCDF file
INTEGER :: dimlen            ! Dimensions size
INTEGER :: varid             ! Variable id of dimension
INTEGER :: ix                ! STASH item index
INTEGER :: i                 ! NetCDF dimension index
INTEGER :: jj                ! Level dimension index
INTEGER :: j                 ! Model level number for index jj
INTEGER :: jb                ! Bottom model level number for vertical mean
INTEGER :: jt                ! Top model level number for vertical mean
INTEGER :: idiag             ! STASH diagnostic address
INTEGER :: section           ! STASH section number
INTEGER :: item              ! STASH item number
INTEGER :: levbot            ! < 0 points to levels list -levbot
                             ! > 0 start level for model level output
INTEGER :: levtop            ! levbot > 0 last level for model level output
INTEGER :: lev_list(nlevp_s) ! Array of output model level numbers
INTEGER :: nlevs             ! Number of input vertical levels for vertical mean
INTEGER :: start(1)          ! Start position of dimension output data
INTEGER :: count1(1)         ! Number of items of dimension output data
INTEGER :: start_bnds(2)     ! Start position of dimension bounds output data
INTEGER :: count_bnds(2)     ! Number of items of dimension bounds output data
INTEGER :: mean_opt          ! Meaning option
INTEGER :: level_code        ! STASH level type code
INTEGER :: cf_extra_info_len ! Length of string containing CF extra information
INTEGER :: ibounds(2,1)      ! Integer bounds values for vertical means

REAL    :: zval(1)               ! Real value of czval
REAL    :: rbounds(2,1)          ! Real bounds values for vertical means
REAL    :: zseak(nlevp_s)        ! Array of output zsea values
REAL    :: zseak_bnds(2,nlevp_s) ! Array of output zsea bounds values
REAL    :: ck(nlevp_s)           ! Array of output C values
REAL    :: ck_bnds(2,nlevp_s)    ! Array of output C bounds values
REAL    :: etak(nlevp_s)         ! Array of output eta values
REAL    :: etak_bnds(2,nlevp_s)  ! Array of output eta bounds values
REAL    :: soil(nlevp_s)         ! Array of output soil level values
REAL    :: soil_bnds(2,nlevp_s)  ! Array of output soil level bounds values
REAL    :: soiltop               ! Top of soil layer
REAL    :: soilbot               ! Bottom of soil layer

LOGICAL :: lexists    ! TRUE if selected attribute exists
LOGICAL :: lvert_mean ! TRUE if field has a vertical mean

CHARACTER(LEN=*), PARAMETER       :: RoutineName = 'NCFILE_WRITE_VERT_VAR'
CHARACTER(LEN=profilename_length) :: domname ! Domain profile name
CHARACTER(LEN=nf90_max_name)      :: dimname ! NetCDF dimension name
CHARACTER(LEN=nf90_max_name)      :: varname ! NetCDF variable name
CHARACTER(LEN=1)                  :: axis    ! Dimension axis type
CHARACTER(LEN=5)                  :: czval   ! String containing height values
                                             ! read from CF extra information

! Dr-Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get numnber of dimensions in NetCDF file
CALL nc_inq_ndims(um_file,ndims)

! Loop over each dimension, i is netCDF dimid
DO i=1,ndims

  ! Get dimension name and length for this dimid
  CALL nc_inq_dim(um_file,i,dimname,dimlen)

  ! Not interested in bounds dimension
  IF (dimname(1:6) == 'bounds') CYCLE

  ! Get id of variable with same name as dimension
  CALL nc_get_varid(um_file,dimname,varid)

  ! Get axis attribute
  CALL nc_get_att(um_file,varid,'axis',axis,lexists)
  IF (.NOT. lexists) CYCLE

  ! This routine is for processing vertical dimensions, so only select those
  IF (axis == 'Z') THEN

    ! Get stash index value
    ix = um_file % nc_meta % stash_index_dim(i)

    IF (printstatus >= prdiag) THEN
      WRITE(umMessage,'(A,A,A,I6,A,I4,A)') &
            'Vertical dimension name ',TRIM(dimname), &
            ' of size ',dimlen,' on unit ',um_file % UNIT,' found'
      CALL umPrint(umMessage,src=RoutineName)
      WRITE(umMessage,'(A,A,A,I5)') &
                'STASH index for vertical dimension ',TRIM(dimname),' is ',ix
      CALL umPrint(umMessage,src=RoutineName)
    END IF

    ! Intialise start and count variables needed to write netCDF data
    start(1) = 1
    count1(1) = dimlen
    start_bnds(1) = 1
    start_bnds(2) = 1
    count_bnds(1) = 2
    count_bnds(2) = dimlen

    ! Get STASH and levels information and meaning options
    section = stlist(st_sect_no_code,ix)
    item = stlist(st_item_code,ix)
    levbot = stlist(st_output_bottom,ix)
    levtop = stlist(st_output_top,ix)
    mean_opt = (stlist(st_gridpoint_code,ix)/block_size)*block_size
    lvert_mean = mean_opt == vert_mean_base .OR. mean_opt == global_mean_base
    level_code = exppxi(atmos_im,section,item,ppx_lv_code)

    ! STASH stores vertical levels as a list of values
    ! Integer values for model levels (including soil model levels)
    ! Real values for all other level types
    SELECT CASE (level_code)
    CASE (st_levels_model_rho,   &
          st_levels_model_theta, &
          st_levels_deep_soil) 
      ! Need to use model level number auxilliary coordinate variable
      ! to match the integer STASH list values
      idiag = stlist(st_diag_address,ix)
      domname = dompro(idom_b(idiag))
      varname = TRIM(domname)//'_'//nc_vert_id_model_level
    CASE DEFAULT
      ! The coordinate variable matches the real STASH list values
      ! for other level types
      varname = dimname
    END SELECT

    ! Get variable id and write level values
    CALL nc_get_varid(um_file,varname,varid)

    IF (level_code == st_levels_single) THEN
      ! Single height level, get height value from STASH to CF metadata
      cf_extra_info_len = LEN_TRIM(cf_extra_info(section,item))
      czval = cf_extra_info(section,item)(8:cf_extra_info_len-1)
      READ(czval,'(F5.0)') zval(1)
      CALL nc_put_var(um_file,varid,zval,start,count1)
    ELSE IF (levbot >= 0) THEN ! model level range
      DO j=levbot,levbot+dimlen-1
        lev_list(j-levbot+1) = j
      END DO
      CALL nc_put_var(um_file,varid,lev_list,start,count1)
    ELSE IF (levbot < 0) THEN ! levels defined using STASH list
      j = -levbot
      IF (llistty(j) == 'R') THEN
        CALL nc_put_var(um_file,varid,rlevlst_s(2:dimlen+1,j),start,count1)
      ELSE IF (llistty(j) == 'I') THEN
        CALL nc_put_var(um_file,varid,levlst_s(2:dimlen+1,j),start,count1)
      END IF
    END IF

    ! If vertical mean specified need bounds values
    IF (lvert_mean) THEN

      ! Get variable id
      varname = TRIM(varname)//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)

      ! Write bounds values
      IF (levbot >= 0) THEN
        ibounds(1,1) = levbot
        ibounds(2,1) = levtop
        CALL nc_put_var(um_file,varid,ibounds,start_bnds,count_bnds)
      ELSE IF (levbot < 0) THEN
        nlevs = stash_levels(1,-levbot)
        IF (llistty(j) == 'R') THEN
          rbounds(1,1) = rlevlst_s(2,-levbot)
          rbounds(2,1) = rlevlst_s(nlevs+1,-levbot)
          CALL nc_put_var(um_file,varid,rbounds,start_bnds,count_bnds)
        ELSE IF (llistty(j) == 'I') THEN
          ibounds(1,1) = levlst_s(2,-levbot)
          ibounds(2,1) = levlst_s(nlevs+1,-levbot)
          CALL nc_put_var(um_file,varid,ibounds,start_bnds,count_bnds)
        END IF
      END IF
    END IF

    ! Write etak,zsea,C values for hybrid height coordinates
    ! and their associated bounds values
    SELECT CASE (level_code)
    CASE (st_levels_model_rho)

      IF (lvert_mean) THEN

        ! jb = bottom model level, jt = top model level of vertical mean
        IF (levbot >= 0) THEN
          jb = levbot
          jt = levtop
        ELSE
          nlevs = stash_levels(1,-levbot)
          jb = levlst_s(2,-levbot)
          jt = levlst_s(nlevs+1,-levbot)
        END IF

        ! Lower bounds values
        IF (jb <= 0) THEN
          ! Level below 1st rho level is surface level
          etak_bnds(1,1) = 0.0
          zseak_bnds(1,1) = 0.0
          ck_bnds(1,1) = 1.0
        ELSE
          etak_bnds(1,1) = a_levdepc(jetatheta+jb-1)
          zseak_bnds(1,1) = a_levdepc(jzseak_theta+jb-1)
          ck_bnds(1,1) = a_levdepc(jck_theta+jb-1)
        END IF

        ! Upper bounds values
        IF (jt > model_levels) THEN ! above top level
          etak_bnds(2,1) = a_levdepc(jetatheta+model_levels)
          zseak_bnds(2,1) = 2.0*a_levdepc(jzseak_theta+model_levels) - &
                                a_levdepc(jzseak_rho+model_levels-1)
          ck_bnds(2,1) = 2.0*a_levdepc(jck_theta+model_levels) - &
                             a_levdepc(jck_rho+model_levels-1)
        ELSE
          etak_bnds(2,1) = a_levdepc(jetatheta+jt)
          zseak_bnds(2,1) = a_levdepc(jzseak_theta+jt)
          ck_bnds(2,1) = a_levdepc(jck_theta+jt)
        END IF

        ! Use bottom level for coordinate values for vertical mean
        IF (jb <= 0) THEN
          ! Level below 1st rho level is surface level
          etak(1) = 0.0
          zseak(1) = 0.0
          ck(1) = 1.0
        ELSE
          etak(1) = a_levdepc(jetarho+jb-1)
          zseak(1) = a_levdepc(jzseak_rho+jb-1)
          ck(1) = a_levdepc(jck_rho+jb-1)
        END IF
      ELSE
        DO jj=1,dimlen

          ! Get model level number j for index jj
          IF (levbot >= 0) THEN
            j = levbot+jj-1
          ELSE
            j = levlst_s(jj+1,-levbot)
          END IF

          ! Lower bounds values
          IF (j <= 0) THEN ! Surface level
            etak_bnds(1,jj) = 0.0
            zseak_bnds(1,jj) = 0.0
            ck_bnds(1,jj) = 1.0
          ELSE
            etak_bnds(1,jj) = a_levdepc(jetatheta+j-1)
            zseak_bnds(1,jj) = a_levdepc(jzseak_theta+j-1)
            ck_bnds(1,jj) = a_levdepc(jck_theta+j-1)
          END IF

          ! Upper bounds values
          IF (j > model_levels) THEN ! above top level
            etak_bnds(2,jj) = a_levdepc(jetatheta+model_levels)
            zseak_bnds(2,jj) = 2.0*a_levdepc(jzseak_theta+model_levels) - &
                                   a_levdepc(jzseak_rho+model_levels-1)
            ck_bnds(2,jj) = 2.0*a_levdepc(jck_theta+model_levels) - &
                                a_levdepc(jck_rho+model_levels-1)
          ELSE
            etak_bnds(2,jj) = a_levdepc(jetatheta+j)
            zseak_bnds(2,jj) = a_levdepc(jzseak_theta+j)
            ck_bnds(2,jj) = a_levdepc(jck_theta+j)
          END IF

          ! Coordinate values
          IF (j <= 0) THEN ! Surface level
            etak(jj) = 0.0
            zseak(jj) = 0.0
            ck(jj) = 1.0
          ELSE IF (j > model_levels) THEN ! above top level
            etak(jj) = etak_bnds(2,jj)
            zseak(jj) = zseak_bnds(2,jj)
            ck(jj) = ck_bnds(2,jj)
          ELSE
            etak(jj) = a_levdepc(jetarho+j-1)
            zseak(jj) = a_levdepc(jzseak_rho+j-1)
            ck(jj) = a_levdepc(jck_rho+j-1)
          END IF
        END DO
      END IF

      ! Get variable id and write eta rho values
      varname = TRIM(domname)//'_'//nc_vert_id_eta_rho
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,etak,start,count1)

      ! Get variable id and write zsea rho values
      varname = TRIM(domname)//'_'//nc_vert_id_zsea_rho
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,zseak,start,count1)

      ! Get variable id and write C rho values
      varname = TRIM(domname)//'_'//nc_vert_id_C_rho
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,ck,start,count1)

      ! Get variable id and write eta rho bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_eta_rho//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,etak_bnds,start_bnds,count_bnds)

      ! Get variable id and write zsea rho bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_zsea_rho//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,zseak_bnds,start_bnds,count_bnds)

      ! Get variable id and write C rho bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_C_rho//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,ck_bnds,start_bnds,count_bnds)

    CASE (st_levels_model_theta)

      IF (lvert_mean) THEN

        ! jb = bottom model level, jt = top model level of vertical mean
        IF (levbot >= 0) THEN
          jb = levbot
          jt = levtop
        ELSE
          nlevs = stash_levels(1,-levbot)
          jb = levlst_s(2,-levbot)
          jt = levlst_s(nlevs+1,-levbot)
        END IF ! above top level
  
        ! Lower bounds values
        IF (jb <= 1) THEN
          etak_bnds(1,1) = 0.0
          zseak_bnds(1,1) = 0.0
          ck_bnds(1,1) = 1.0
        ELSE
          etak_bnds(1,1) = a_levdepc(jetarho+jb-1)
          zseak_bnds(1,1) = a_levdepc(jzseak_rho+jb-1)
          ck_bnds(1,1) = a_levdepc(jck_rho+jb-1)
        END IF

        ! Upper bounds values
        IF (jt >= model_levels) THEN ! above top level
          etak_bnds(2,1) = a_levdepc(jetatheta+model_levels)
          zseak_bnds(2,1) = 2.0*a_levdepc(jzseak_theta+model_levels) - &
                                a_levdepc(jzseak_rho+model_levels-1)
          ck_bnds(2,1) = 2.0*a_levdepc(jck_theta+model_levels) - &
                             a_levdepc(jck_rho+model_levels-1)
        ELSE
          etak_bnds(2,1) = a_levdepc(jetarho+jt)
          zseak_bnds(2,1) = a_levdepc(jzseak_rho+jt)
          ck_bnds(2,1) = a_levdepc(jck_rho+jt)
        END IF

        ! Use bottom level for coordinate values for vertical mean
        IF (jb > model_levels) THEN
          etak(1) = etak_bnds(2,1)
          zseak(1) = zseak_bnds(2,1)
          ck(1) = ck_bnds(2,1)
        ELSE
          etak(1) = a_levdepc(jetatheta+jb)
          zseak(1) = a_levdepc(jzseak_theta+jb)
          ck(1) = a_levdepc(jck_theta+jb)
        END IF
      ELSE
        DO jj=1,dimlen

          ! Get model level number j for index jj
          IF (levbot >= 0) THEN
            j = levbot+jj-1
          ELSE
            j = levlst_s(jj+1,-levbot)
          END IF

          ! Lower bounds values
          ! Note that the lowest theta level (=1) has a lower layer boundary at
          ! the surface, as set explicitly in the model interface to physics.
          IF (j <= 1) THEN
            etak_bnds(1,jj) = 0.0
            zseak_bnds(1,jj) = 0.0
            ck_bnds(1,jj) = 1.0
          ELSE
            etak_bnds(1,jj) = a_levdepc(jetarho+j-1)
            zseak_bnds(1,jj) = a_levdepc(jzseak_rho+j-1)
            ck_bnds(1,jj) = a_levdepc(jck_rho+j-1)
          END IF

          ! Upper bounds values
          IF (j >= model_levels) THEN ! above top level
            etak_bnds(2,jj) = a_levdepc(jetatheta+model_levels)
            zseak_bnds(2,jj) = 2.0*a_levdepc(jzseak_theta+model_levels) - &
                                   a_levdepc(jzseak_rho+model_levels-1)
            ck_bnds(2,jj) = 2.0*a_levdepc(jck_theta+model_levels) - &
                                a_levdepc(jck_rho+model_levels-1)
          ELSE
            etak_bnds(2,jj) = a_levdepc(jetarho+j)
            zseak_bnds(2,jj) = a_levdepc(jzseak_rho+j)
            ck_bnds(2,jj) = a_levdepc(jck_rho+j)
          END IF

          ! Coordinate values
          IF (j > model_levels) THEN
            etak(jj) = etak_bnds(2,jj)
            zseak(jj) = zseak_bnds(2,jj)
            ck(jj) = ck_bnds(2,jj)
          ELSE
            etak(jj) = a_levdepc(jetatheta+j)
            zseak(jj) = a_levdepc(jzseak_theta+j)
            ck(jj) = a_levdepc(jck_theta+j)
          END IF
        END DO
      END IF

      ! Get variable id and write eta theta values
      varname = TRIM(domname)//'_'//nc_vert_id_eta_theta
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,etak,start,count1)

      ! Get variable id and write zsea theta values
      varname = TRIM(domname)//'_'//nc_vert_id_zsea_theta
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,zseak,start,count1)

      ! Get variable id and write C theta values
      varname = TRIM(domname)//'_'//nc_vert_id_C_theta
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,ck,start,count1)

      ! Get variable id and write eta theta bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_eta_theta//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,etak_bnds,start_bnds,count_bnds)

      ! Get variable id and write zsea theta bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_zsea_theta//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,zseak_bnds,start_bnds,count_bnds)

      ! Get variable id and write C theta bounds values
      varname = TRIM(domname)//'_'//nc_vert_id_C_theta//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,ck_bnds,start_bnds,count_bnds)

    CASE (st_levels_deep_soil)

      ! Work out layer boundaries from soil thickness
      ! Coordinate values are layer mid-points
      IF (levbot > 0) THEN
        soiltop = 0.0
        DO j=1,levtop
          soilbot = soiltop + a_levdepc(jsoil_thickness+j-1)
          IF (j >= levbot) THEN
            soil(j-levbot+1) = 0.5*(soiltop+soilbot)
            soil_bnds(1,j-levbot+1) = soiltop
            soil_bnds(2,j-levbot+1) = soilbot
          END IF
          soiltop = soilbot
        END DO
      ELSE
        jj = 1
        soiltop = 0.0
        DO j=1,nlevp_s
          soilbot = soiltop + a_levdepc(jsoil_thickness+j-1)
          IF (j == levlst_s(jj+1,-levbot)) THEN
            soil(jj) = 0.5*(soiltop+soilbot)
            soil_bnds(1,jj) = soiltop
            soil_bnds(2,jj) = soilbot
            IF (jj == dimlen) EXIT
            jj = jj+1
          END IF
          soiltop = soilbot
        END DO
      END IF

      ! Get variable id and write soil level values
      CALL nc_get_varid(um_file,dimname,varid)
      CALL nc_put_var(um_file,varid,soil,start,count1)

      ! Get variable id and write soil level bounds values
      varname = TRIM(dimname)//'_bounds'
      CALL nc_get_varid(um_file,varname,varid)
      CALL nc_put_var(um_file,varid,soil_bnds,start_bnds,count_bnds)

    END SELECT

  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ncfile_write_vert_var

END MODULE ncfile_write_vert_var_mod
