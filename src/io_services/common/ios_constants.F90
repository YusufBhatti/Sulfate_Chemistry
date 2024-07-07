! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Module defining interface constants and variables for clients and servers

MODULE IOS_Constants

USE UM_Types, ONLY: integer32
USE mpl, ONLY: mpl_integer4
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!------------------------------------------------------------------
!  IO Server - parameters
!------------------------------------------------------------------

! Global server limit
INTEGER, PARAMETER  :: IOS_maxServers                       = 512

! Actions
INTEGER, PARAMETER  :: IOS_Action_Unset                     = 0
INTEGER, PARAMETER  :: IOS_Action_Open                      = 1
INTEGER, PARAMETER  :: IOS_Action_Close                     = 2
INTEGER, PARAMETER  :: IOS_Action_Read64                    = 3
INTEGER, PARAMETER  :: IOS_Action_Read32                    = 4
INTEGER, PARAMETER  :: IOS_Action_Write64                   = 5
INTEGER, PARAMETER  :: IOS_Action_Write32                   = 6
INTEGER, PARAMETER  :: IOS_Action_Setpos                    = 7
INTEGER, PARAMETER  :: IOS_Action_Flush                     = 8
INTEGER, PARAMETER  :: IOS_Action_Finish                    = 9
INTEGER, PARAMETER  :: IOS_Action_Process                   = 10
INTEGER, PARAMETER  :: IOS_Action_StashInitPPLookup         = 11
INTEGER, PARAMETER  :: IOS_Action_Write_PP_Prepacked        = 12
INTEGER, PARAMETER  :: IOS_Action_StashWritePPData          = 13
INTEGER, PARAMETER  :: IOS_Action_StashWritePPLookup        = 14
INTEGER, PARAMETER  :: IOS_Action_StashSetPos               = 15
INTEGER, PARAMETER  :: IOS_Action_StashInitModel            = 16
INTEGER, PARAMETER  :: IOS_Action_StashInitHeader           = 17
INTEGER, PARAMETER  :: IOS_Action_StashSetHeader            = 18
INTEGER, PARAMETER  :: IOS_Action_Config                    = 19
INTEGER, PARAMETER  :: IOS_Action_Read64_Integer            = 20
INTEGER, PARAMETER  :: IOS_Action_Read32_Integer            = 21
INTEGER, PARAMETER  :: IOS_Action_Enquire                   = 22
INTEGER, PARAMETER  :: IOS_Action_Sync                      = 23
INTEGER, PARAMETER  :: IOS_Action_Sync_Barrier              = 24
INTEGER, PARAMETER  :: IOS_Action_Fence                     = 25
INTEGER, PARAMETER  :: IOS_Action_StashWriteDumpData        = 26
INTEGER, PARAMETER  :: IOS_Action_DumpInitModel             = 27
INTEGER, PARAMETER  :: IOS_Action_Assign_Unit               = 28
INTEGER, PARAMETER  :: IOS_Action_LoadStatus                = 29
INTEGER, PARAMETER  :: IOS_Action_MergePPLookup             = 30
INTEGER, PARAMETER  :: IOS_Action_FileOp                    = 31
INTEGER, PARAMETER  :: IOS_Num_Actions                      = 31


CHARACTER (LEN=*),PARAMETER::                                                &
    IOS_Action_strings=                                                      &
    "IOS_Action_Open                "//                                      &
    "IOS_Action_Close               "//                                      &
    "IOS_Action_Read64              "//                                      &
    "IOS_Action_Read32              "//                                      &
    "IOS_Action_Write64             "//                                      &
    "IOS_Action_Write32             "//                                      &
    "IOS_Action_Setpos              "//                                      &
    "IOS_Action_Flush               "//                                      &
    "IOS_Action_Finish              "//                                      &
    "IOS_Action_Process             "//                                      &
    "IOS_Action_StashInitPPLookup   "//                                      &
    "IOS_Action_Write_PP_Prepacked  "//                                      &
    "IOS_Action_StashWritePPData    "//                                      &
    "IOS_Action_StashWritePPLookup  "//                                      &
    "IOS_Action_StashSetPos         "//                                      &
    "IOS_Action_StashInitModel      "//                                      &
    "IOS_Action_StashInitHeader     "//                                      &
    "IOS_Action_StashSetHeder       "//                                      &
    "IOS_Action_Config              "//                                      &
    "IOS_Action_Read64_Integer      "//                                      &
    "IOS_Action_Read32_Integer      "//                                      &
    "IOS_Action_Enquire             "//                                      &
    "IOS_Action_Sync                "//                                      &
    "IOS_Action_Sync_Barrier        "//                                      &
    "IOS_Action_Fence               "//                                      &
    "IOS_Action_StashWriteDumpData  "//                                      &
    "IOS_Action_DumpInitModel       "//                                      &
    "IOS_Action_Assign_Unit         "//                                      &
    "IOS_Action_LoadStatus          "//                                      &
    "IOS_Action_MergePPLookup       "//                                      &
    "IOS_Action_FileOp              "//                                      &
    "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"//                                      &
    "YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY"//                                      &
    "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
!    1234567890123456789012345678901 <- 31

! This is the length of each string subsection above
! Since the strings have a similar name to the parameter we should
! limit the size of the parameter name to 31 since Fortran standard has
! a maximum length of 31 for names.
INTEGER, PARAMETER  :: IOS_Action_Strlen         = 31

! Units
INTEGER, PARAMETER  :: IOS_Unit_Closed      = 0
INTEGER, PARAMETER  :: IOS_Unit_Open_RO     = 1
INTEGER, PARAMETER  :: IOS_Unit_Open_RW     = 2

! Server
INTEGER, PARAMETER  :: IOS_No_Server        =-1
INTEGER, PARAMETER  :: IOS_All_Units        =-1
INTEGER, PARAMETER  :: IOS_No_Location      =-1
INTEGER, PARAMETER  :: IOS_No_Unit          =-1

! MPI Tags
INTEGER, PARAMETER  :: IOS_Request_Tag_Gap       = 2000
INTEGER, PARAMETER  :: IOS_Request_Tag_Base      = 2000
INTEGER, PARAMETER  :: IOS_Request_Tag_Str       =                           &
    IOS_Request_Tag_Base    +IOS_Request_Tag_Gap
INTEGER, PARAMETER  :: IOS_Request_Tag_Payload   =                           &
    IOS_Request_Tag_Str     +IOS_Request_Tag_Gap
INTEGER, PARAMETER  :: IOS_Request_Tag_Express   =                           &
    IOS_Request_Tag_Payload +IOS_Request_Tag_Gap
INTEGER, PARAMETER  :: IOS_ControlSync_Tag  = 74

! Metadata sizes
INTEGER, PARAMETER  :: IOS_string_max       = 256
INTEGER, PARAMETER  :: IOS_md_len           = 10

! Action policy specification
INTEGER, PARAMETER  :: IOS_Policy_Strict      = 1
INTEGER, PARAMETER  :: IOS_Policy_Consistent  = 2
INTEGER, PARAMETER  :: IOS_Policy_Tolerant    = 3
INTEGER, PARAMETER  :: IOS_Policy_Lax         = 4
! Policy specifies the degree of sloppyness that will be allowed
! before the IO server aborts with an error. This ranges from
! Strict (everything must be consistent and correctly sequenced)
! to Lax. The specifics will vary
! between actions, so inspect IOS_writer() for details. This parameter
! is a hint, and so may be freely ignored by any particular action, setting
! policy_strict is not a guarantee of correctness.



  ! Code for timer activation
INTEGER,PARAMETER  :: IOS_Timer_Activate  = 6

! Type information, the client queue will use 4-byte integers
! as the unit (TU) for transactions, so define some helpful quantities
INTEGER, PARAMETER :: IOS_BytesPerTU             = 4
INTEGER, PARAMETER :: IOS_tutype                 = mpl_integer4
INTEGER, PARAMETER :: IOS_tukind                 = integer32

INTEGER, PARAMETER :: IOS_BytesPerWord32         = IOS_BytesPerTU
INTEGER, PARAMETER :: IOS_BytesPerWord64         = IOS_BytesPerWord32*2
INTEGER, PARAMETER :: IOS_tuperword32            =                           &
    IOS_BytesPerWord32/IOS_BytesPerTU
INTEGER, PARAMETER :: IOS_tuperword64            =                           &
    IOS_BytesPerWord64/IOS_BytesPerTU

!!States that a send queue slot can be in
INTEGER, PARAMETER :: IOS_queue_slot_unused      = 0
INTEGER, PARAMETER :: IOS_queue_slot_initialized = 1
INTEGER, PARAMETER :: IOS_queue_slot_partfilled  = 2
INTEGER, PARAMETER :: IOS_queue_slot_dispatched  = 3
INTEGER, PARAMETER :: IOS_queue_slot_ready       = 4

!!Optional values for IO Server to unit allocation
INTEGER, PARAMETER :: IOS_Unit_Alloc_Static         = 1
INTEGER, PARAMETER :: IOS_Unit_Alloc_AtFirstUse     = 2
INTEGER, PARAMETER :: IOS_Unit_Alloc_Dynamic_Rotate = 3
INTEGER, PARAMETER :: IOS_Unit_Alloc_Dynamic_LB     = 4

!!Length of character strings for messages
INTEGER, PARAMETER :: IOS_err_strlen             = errormessagelength

!!Length of character strings for MPI-IO hints
INTEGER, PARAMETER :: IOS_info_strlen            = 20

!!Parameters for describing server load.
INTEGER, PARAMETER :: loadBalanceDataSize      = 2
INTEGER, PARAMETER :: loadBalanceQueueLen      = 1
INTEGER, PARAMETER :: loadBalanceQueueData     = 2

!!What things can we do with reads?
INTEGER, PARAMETER :: IOS_Read_NoBroadCast          = 0
INTEGER, PARAMETER :: IOS_Read_Broadcast            = 1

END MODULE IOS_constants
