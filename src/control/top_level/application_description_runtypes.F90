! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: Application_Description_Runtypes ---------------------------
!
!  Purpose: The application names and run-time types.
!           Used to enable the UM code to determine
!           the executable target (program main),
!           build type, sizes, and parallel nature. This facilitates
!           reduction in preprocessor based compilation in favour of
!           software switching
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE Application_Description_Runtypes

IMPLICIT NONE

INTEGER, PARAMETER :: exe_unknown           = 1
INTEGER, PARAMETER :: exe_UM                = 2
INTEGER, PARAMETER :: exe_RCF               = 3
INTEGER, PARAMETER :: exe_scm               = 4
! Pumf (formerly having value 5) is now retired 
! ConvIEEE (formerly having value 6) is now retired 
INTEGER, PARAMETER :: exe_combine           = 7
INTEGER, PARAMETER :: exe_merge             = 8
! Cumf (formerly having value 9) is now retired 
INTEGER, PARAMETER :: exe_hreset            = 10
! Fieldcos (formerly having value 11) is now retired
! FieldOp (formerly having value 12) is now retired 
! FieldCalc (formerly having value 13) is now retired
INTEGER, PARAMETER :: exe_convpp            = 14
INTEGER, PARAMETER :: exe_setup             = 15
INTEGER, PARAMETER :: exe_fldmod            = 16
INTEGER, PARAMETER :: exe_hprint            = 17
INTEGER, PARAMETER :: exe_pptoanc           = 18
INTEGER, PARAMETER :: exe_pickup            = 19
! MakeBC (formerly having value 20) is now retired
! Frames (formerly having value 21) is now retired
INTEGER, PARAMETER :: exe_vomext            = 22
INTEGER, PARAMETER :: exe_crmstyle_coarse_grid = 23
INTEGER, PARAMETER :: exe_createbc          = 24
INTEGER, PARAMETER :: num_exe_types         = 24

CHARACTER (LEN=*),PARAMETER::&
    UM_name=&
    "Application Unknown     "//&
    "Unified Model           "//&
    "Reconfiguration         "//&
    "Single Column Model     "//&
    "                        "//&
    "                        "//&
    "Combine                 "//&
    "Merge                   "//&
    "                        "//&
    "HReset                  "//&
    "                        "//&
    "                        "//&
    "                        "//&
    "Convert PP              "//&
    "Setup                   "//&
    "Field Mod               "//&
    "HPrint                  "//&
    "PP to Ancilary          "//&
    "Pickup                  "//&
    "                        "//&
    "                        "//&
    "Vom Extract             "//&
    "CRM style regridding    "//&
    "CreateBC                "
!    123456789012345678901234 <- 24

! This is the length of each string subsection above
INTEGER, PARAMETER  :: UM_Name_Len = 24

INTEGER :: exe_type         =exe_unknown
LOGICAL :: exe_is_parallel  =.FALSE.
LOGICAL :: exe_addr_64      =.FALSE.
LOGICAL :: exe_data_64      =.FALSE.


END MODULE Application_Description_Runtypes
