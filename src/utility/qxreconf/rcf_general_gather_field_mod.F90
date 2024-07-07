! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Gathers any type of field from many processors to one processor

MODULE Rcf_General_Gather_Field_Mod

IMPLICIT NONE
!  Subroutine Rcf_General_Gather_Field - generalised gather field
!
! Description:
!   GAthers any type of field from many processors to one processor by
!   consideration of type - only copes with full fields
!
! Method:
!   Land compressed fields are uncompressed, gathered and recompressed
!   Zonal fields have their own routine
!   All other fields are "normal" and use rcf_gather_field
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_GENERAL_GATHER_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_General_Gather_Field( local_field,  global_field , &
                                     local_size,   global_size ,  &
                                     stash_record, gather_pe )

USE Rcf_Lsm_Mod, ONLY:     &
    Local_Land_Field,        &
    Local_Atmos_Landmask,    &
    Glob_Atmos_Landmask

USE UM_ParVars, ONLY:                          &
    blsizep,                blsizeu,           &
    glsizeu,                glsizep,           &
    lasize,                 gc_all_proc_group, &
    glsize,                 g_blsize

USE UM_ParCore, ONLY: &
    mype

USE Field_Types, ONLY: &
    fld_type_p

USE Ereport_Mod, ONLY: &
          Ereport

USE Rcf_Gather_Zonal_Field_Mod, ONLY: &
          Rcf_Gather_Zonal_Field

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
INTEGER, INTENT(IN)   :: Global_Size   ! size of global_field
INTEGER, INTENT(IN)   :: Gather_PE     ! PE to collect global data on
INTEGER, INTENT(OUT)  :: Local_Size    ! size of local field

REAL, INTENT(IN)      :: Local_Field(*)      ! local part of field
REAL, INTENT(OUT)     :: Global_Field( Global_Size )

TYPE( STM_Record_Type ), INTENT(IN) :: Stash_Record

! Local variables
INTEGER       :: dummy            ! ignored argument
INTEGER       :: fld_type         ! P or U or V
INTEGER       :: local_x          ! local x size
INTEGER       :: local_y          ! local y size
INTEGER       :: global_x         ! global x size
INTEGER       :: global_y         ! global y size
INTEGER       :: ErrorStatus      ! error code

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_GENERAL_GATHER_FIELD'
CHARACTER (LEN=errormessagelength)   :: Cmessage

REAL          :: buf_expand( glsizep(1) * glsizep(2) )
REAL          :: buf_expand_local(          &
    lasize(1,fld_type_p,halo_type_single) * &
    lasize(2,fld_type_p,halo_type_single) )
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!===================================================================

! Choose action depending on grid code
SELECT CASE( Stash_Record % grid_type )

  !------------------------------------------------------------------
  ! Land compressed fields
  !------------------------------------------------------------------

CASE ( ppx_atm_compressed )

  ! Unpack the local field out to full (local) field size and
  ! put this into the array buf_expand_local

  CALL expand_from_mask(buf_expand_local, local_field, &
                        local_atmos_landmask,          &
                        lasize(1,fld_type_p,halo_type_single)* &
                        lasize(2,fld_type_p,halo_type_single), &
                        dummy)

  ! Now gather in all the processors local fields into the global
  ! field (array buf_expand)

  CALL Gather_Field(buf_expand_local, buf_expand,          &
                    lasize(1,fld_type_p,halo_type_single), &
                    lasize(2,fld_type_p,halo_type_single), &
                    glsizep(1),        glsizep(2),         &
                    fld_type_p,        halo_type_single,   &
                    gather_pe,         gc_all_proc_group )

  ! And now pack the global field (buf_expand) back to land points
  ! and put into the array GLOBAL_FIELD.

  IF (mype  ==  0) THEN
    CALL compress_to_mask(buf_expand, global_field, &
                        glob_atmos_landmask,      &
                        glsizep(1)*glsizep(2),dummy)
  END IF

  Local_Size = local_land_field

  !-------------------------------------------------------------------
  ! Zonal fields
  !-------------------------------------------------------------------
CASE ( ppx_atm_tzonal, ppx_atm_uzonal )

  IF ( Stash_Record % grid_type == ppx_atm_tzonal ) THEN   ! P grid
    global_y = glsizep(2)
    local_y  = blsizep(2)

  ELSE
    global_y = glsizeu(2)
    local_y  = blsizeu(2)
  END IF

  CALL Rcf_Gather_Zonal_Field( Local_Field, Global_Field,            &
                               local_y,     global_y,                &
                               1,           Stash_Record % grid_type,&
                               Gather_PE )

  Local_Size = local_y
  !-------------------------------------------------------------------
  ! "Normal" fields
  !-------------------------------------------------------------------
CASE &
  !     atmosphere grids
         ( ppx_atm_tall,    &! Atmos T points
           ppx_atm_tland,   &! Atmos T land points
           ppx_atm_tsea,    &! Atmos T sea points
           ppx_atm_uall,    &! Atmos U points
           ppx_atm_uland,   &! Atmos U land points
           ppx_atm_usea,    &! Atmos U sea points
           ppx_atm_cuall,   &! Atmos C grid U pts
           ppx_atm_cvall,   &! Atmos C grid V pts
           ppx_atm_ozone,   &! Atmos ozone field
           ppx_atm_river)    ! Atmos river-routing field

  fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

  global_x = glsize(1, fld_type)
  global_y = glsize(2, fld_type)
  local_x  = g_blsize(1, fld_type, mype)
  local_y  = g_blsize(2, fld_type, mype)

! DEPENDS ON: gather_field
  CALL Gather_Field( Local_Field, Global_Field,      &
                     local_x,     local_y,           &
                     global_x,    global_y,          &
                     fld_type,    halo_type_single,  &
                     Gather_PE,   gc_all_proc_group )

  Local_Size = local_x * local_y

  !-------------------------------------------------------------------
  ! Any other type of field
  !-------------------------------------------------------------------
CASE DEFAULT

  ErrorStatus = 10
  Cmessage = 'Field type not recognised for Gather'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_General_Gather_Field
END MODULE Rcf_General_Gather_Field_Mod



