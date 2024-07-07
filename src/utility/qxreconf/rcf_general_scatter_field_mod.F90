! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Scatters any type of field from one processor to many processors

MODULE Rcf_General_Scatter_Field_Mod
IMPLICIT NONE

!  Subroutine Rcf_General_Scatter_Field
!
! Description:
!   Chooses scattering method for a field based on grid_code.
!
! Method:
!   Land compressed fields are uncompressed, scattered and recompressed
!   Zonal fields have specialist subroutines
!   Other fields are "normal" and can be dealt with by Scatter_Field
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_GENERAL_SCATTER_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_General_Scatter_Field ( local_field,  global_field , &
                                       local_size,   global_size ,  &
                                       stash_record, scatter_pe )

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Lsm_Mod, ONLY:   &
    Local_Land_Field,     &
    Local_Atmos_Landmask, &
    Glob_Atmos_Landmask

USE UM_ParVars, ONLY:                          &
    blsizep,                blsizeu,           &
    blsizev,                glsizeu,           &
    glsizev,                glsizep,           &
    lasize,                 gc_all_proc_group, &
    glsize,                 g_blsize

USE UM_ParCore, ONLY: &
    mype

USE Field_Types, ONLY: &
    fld_type_p

USE Rcf_Scatter_Zonal_Field_Mod, ONLY: &
    Rcf_Scatter_Zonal_Field

USE Rcf_Global_To_Local_Mod, ONLY: &
    Rcf_Get_Fld_Type

USE Rcf_Ppx_Info_Mod, ONLY: &
    STM_Record_Type

USE mask_compression, ONLY: compress_to_mask, expand_from_mask

USE UM_ParParams

USE cppxref_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)   :: Global_Size  ! size of global_field
INTEGER, INTENT(IN)   :: Scatter_PE   ! PE on which global field lives
INTEGER, INTENT(OUT)  :: Local_Size   ! size of local_field

REAL, INTENT(IN)      :: Global_Field( Global_Size ) ! field to scatter
REAL, INTENT(OUT)     :: Local_Field( * )

TYPE( STM_Record_Type ), INTENT(IN) :: Stash_Record

! Local variables
INTEGER               :: dummy       ! ignored argument
INTEGER               :: fld_type    ! P or U or V
INTEGER               :: local_x     ! local x size
INTEGER               :: local_y     ! local y size
INTEGER               :: global_x    ! global x size
INTEGER               :: global_y    ! global y size
INTEGER               :: ErrorStatus ! error code
CHARACTER (LEN=*), PARAMETER   :: RoutineName='RCF_GENERAL_SCATTER_FIELD'//&
        'Field'
CHARACTER (LEN=errormessagelength)             :: Cmessage

REAL                  :: buf_expand( &
    glsizep(1) * &
    glsizep(2) )
REAL                  :: buf_expand_local     &
    ( lasize(1,fld_type_p,halo_type_single) * &
      lasize(2,fld_type_p,halo_type_single) )
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!===================================================================

! Choose action depending on grid code
SELECT CASE( Stash_Record % grid_type )

  !-------------------------------------------------------------------
  ! Land compressed fields
  !-------------------------------------------------------------------
CASE ( ppx_atm_compressed )

  IF (mype  ==  scatter_pe) THEN

    CALL expand_from_mask(buf_expand, global_field,    &
                          glob_atmos_landmask,         &
                          glsizep(1)*                  &
                          glsizep(2),dummy)

    ! SCATTER_PE now contains the expanded version of the full field

  END IF

  ! Now scatter this to all the other processors, putting the local
  ! part of the field into the array buf_expand_local

!DEPENDS ON: scatter_field
  CALL Scatter_Field(buf_expand_local , buf_expand,         &
                     lasize(1,fld_type_p,halo_type_single), &
                     lasize(2,fld_type_p,halo_type_single), &
                     glsizep(1),     glsizep(2),            &
                     fld_type_p,     halo_type_single,      &
                     scatter_pe,     gc_all_proc_group )

  ! Pack the local field down to local land points and put
  ! the packed field into LOCAL_FIELD

  CALL compress_to_mask(buf_expand_local,local_field,           &
                      local_atmos_landmask,                     &
                      lasize(1,fld_type_p,halo_type_single)*    &
                      lasize(2,fld_type_p,halo_type_single),    &
                      dummy)

  local_size = local_land_field

  !-------------------------------------------------------------------
  ! Zonal fields
  !-------------------------------------------------------------------
CASE ( ppx_atm_tzonal, ppx_atm_uzonal )

  IF ( Stash_Record % grid_type == ppx_atm_tzonal ) THEN
    global_y = glsizep(2)
    local_y  = blsizep(2)

  ELSE                ! U grid
    global_y = glsizeu(2)
    local_y  = blsizeu(2)
  END IF

  local_size = local_y

  CALL Rcf_Scatter_Zonal_Field( Local_Field, Global_Field,           &
                                local_y,     global_y,               &
                                1,          Stash_Record % grid_type,&
                                Scatter_PE )

  !-------------------------------------------------------------------
  ! Normal fields
  !-------------------------------------------------------------------
CASE &
  !     atmosphere grids
         ( ppx_atm_tall,   &! Atmos T points
           ppx_atm_tland,  &! Atmos T land points
           ppx_atm_tsea,   &! Atmos T sea points
           ppx_atm_uall,   &! Atmos U points
           ppx_atm_uland,  &! Atmos U land points
           ppx_atm_usea,   &! Atmos U sea points
           ppx_atm_cuall,  &! Atmos C grid U pts
           ppx_atm_cvall,  &! Atmos C grid V pts
           ppx_atm_ozone,  &! Atmos ozone field
           ppx_atm_river)  ! Atmos river routing field


  fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

  global_x = glsize(1, fld_type)
  global_y = glsize(2, fld_type)
  local_x  = g_blsize(1, fld_type, mype)
  local_y  = g_blsize(2, fld_type, mype)

  local_size = local_x * local_y
!DEPENDS ON: scatter_field
  CALL Scatter_Field(Local_Field, Global_Field,     &
                     local_x,     local_y,          &
                     global_x,    global_y,         &
                     fld_type,    halo_type_single, &
                     Scatter_PE,  gc_all_proc_group )

  !-------------------------------------------------------------------
  ! Any other type of field
  !-------------------------------------------------------------------
CASE DEFAULT

  ErrorStatus = 10
  Cmessage='Field type not recognized for Scatter'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_General_Scatter_Field
END MODULE Rcf_General_Scatter_Field_Mod
