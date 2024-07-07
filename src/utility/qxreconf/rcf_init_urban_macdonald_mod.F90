! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_Init_Urban_MacDonald_Mod

IMPLICIT NONE

! Subroutine Rcf_Init_Urban_MacDonald
!
! Description:
!   Initialises urban roughness length and displacement height prognostics
!
! Method:
!   Uses MacDonald 1998 and replaces the equivalent UM_JULES init_urban section.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_INIT_URBAN_MACDONALD_MOD'

CONTAINS

SUBROUTINE Rcf_Init_Urban_MacDonald( fields_out, field_count_out, hdr_out, &
                                     data_source )

USE Rcf_Locate_Mod, ONLY:       &
    Rcf_Locate

USE decomp_params, ONLY:        &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY:   &
    Field_Type

USE Rcf_UMhead_Mod, ONLY:       &
    um_header_type

USE um_stashcode_mod, ONLY:                                 &
    stashcode_prog_sec, stashcode_urbhgt, stashcode_urbhwr, &
    stashcode_urbwrr, stashcode_urbdisp, stashcode_urbztm

USE umPrintMgr, ONLY:           &
    umPrint,                    &
    umMessage,                  &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY:   &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY:  &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Write_Field_Mod, ONLY:  &
    Rcf_Write_Field

USE Rcf_Data_Source_Mod, ONLY:  &
    Already_Processed, data_source_type

USE urban_param, ONLY:          &
    a, cdz, kappa2, z0m_mat

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_out
TYPE( data_source_type ), POINTER  :: data_source( : )

INTEGER, INTENT(IN)                :: field_count_out

! Local variables
INTEGER                            :: pos
REAL, ALLOCATABLE ::                                                     &
   sc_hwr(:),       & ! working variable
   d_h(:)             ! working variable

TYPE( field_type ), POINTER :: hgt, hwr, wrr, ztm, disp

INTEGER                              :: l    ! looper
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_INIT_URBAN_MACDONALD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,'(a,2(2x,i6))')                              &
     'Initialising urban roughness length and ' //             &
     'displacement height using MacDonald 1998: Stashcodes',   &
     stashcode_urbztm, stashcode_urbdisp
  CALL umPrint(umMessage,src='Rcf_Init_Urban_MacDonald')
END IF

! The following fields should be present as l_urban2t
! Read fields used to calculate ztm and disp
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_urbhgt,        &
                  fields_out, field_count_out, pos )
hgt => fields_out( pos )
CALL Rcf_Alloc_Field( hgt )
CALL Rcf_Read_Field( hgt, hdr_out, decomp_rcf_output )

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_urbhwr,        &
                  fields_out, field_count_out, pos )
hwr => fields_out( pos )
CALL Rcf_Alloc_Field( hwr )
CALL Rcf_Read_Field( hwr, hdr_out, decomp_rcf_output )

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_urbwrr,        &
                  fields_out, field_count_out, pos )
wrr => fields_out( pos )
CALL Rcf_Alloc_Field( wrr )
CALL Rcf_Read_Field( wrr, hdr_out, decomp_rcf_output )

! Allocate fields for output
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_urbztm,        &
                  fields_out, field_count_out, pos )
ztm => fields_out( pos )
CALL Rcf_Alloc_Field( ztm )
! Set to already processed to avoid being called twice as routine calculates
! both disp and ztm at the same time.
data_source( pos ) % source = already_processed

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_urbdisp,        &
                  fields_out, field_count_out, pos )
disp => fields_out( pos )
CALL Rcf_Alloc_Field( disp )
! Set to already processed to avoid being called twice as routine calculates
! both disp and ztm at the same time.
data_source( pos ) % source = already_processed

ALLOCATE( sc_hwr( hgt % level_size ) )
ALLOCATE( d_h( hgt % level_size ) )

sc_hwr(:) = 0.5 * ( hwr % DATA(:,1) / (2.0 * ATAN(1.0)) )
d_h(:)    = 1.0 - wrr % DATA(:,1) * ( a**(wrr % DATA(:,1) - 1.0) )
DO l = 1, hgt % level_size
  ! Urban present (set to 1.0 otherwise) > 0 also as a check
  IF ( wrr % DATA(l,1) > 0.0 .AND. wrr % DATA(l,1) < 1.0 ) THEN
    disp % DATA(l,1)   = d_h(l) * hgt % DATA(l,1)
    ztm % DATA(l,1)    = (cdz * (1.0 - d_h(l)) *                            &
       sc_hwr(l) * wrr % DATA(l,1) / kappa2)**(-0.5)
    ztm % DATA(l,1)    = (1.0 - d_h(l))*EXP(-ztm % DATA(l,1))
    ztm % DATA(l,1)    = ztm % DATA(l,1) * hgt % DATA(l,1)
    ztm % DATA(l,1)    = MAX(ztm % DATA(l,1),z0m_mat)
  ELSE
    disp % DATA(l,1) = 0.0
    ztm % DATA(l,1)  = 0.0
  END IF
END DO

CALL Rcf_Write_Field( ztm,  hdr_out, decomp_rcf_output )
CALL Rcf_Write_Field( disp, hdr_out, decomp_rcf_output )

DEALLOCATE( sc_hwr )
DEALLOCATE( d_h )
CALL Rcf_Dealloc_Field( hgt )
CALL Rcf_Dealloc_Field( hwr )
CALL Rcf_Dealloc_Field( wrr )
CALL Rcf_Dealloc_Field( ztm )
CALL Rcf_Dealloc_Field( wrr )
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_Urban_MacDonald

END MODULE Rcf_Init_Urban_MacDonald_Mod
