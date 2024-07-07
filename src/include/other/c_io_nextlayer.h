#if !defined(C_IO_NEXTLAYER_H)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#define C_IO_NEXTLAYER_H

#include "c_io_layers.h"
#include "c_io.h"
#include "um_compile_diag_suspend.h"

/* Sanity Check the layer macros                                              */

#if (CIO_THIS_LAYER < 0) || (CIO_THIS_LAYER >= CIO_NUMBERLAYERS)               \
                         || !defined(CIO_THIS_LAYER)                           \
                         || ((CIO_THIS_LAYER == CIO_NO_LAYER)                  \
                             && !defined(CIO_ROOT_LAYER))

#error This layer is invalid CIO_THIS_LAYER

#endif

#if defined(CIO_NEXT_LAYER)
#undef CIO_NEXT_LAYER
#endif

/* For now, BYTESWAP_IO should always be defined */
#if !defined(BYTESWAP_IO)
#define BYTESWAP_IO
#endif

/*----------------------------------------------------------------------------*/
/*  Meta-layers                                                               */
/*----------------------------------------------------------------------------*/
/* Meta-layers are defined from bottom upwards                                */

/*
 * The following meta-layer cascade scheme is used:
 * CIO_TOP -> CIO_INFO -> CIO_BUFFER -> CIO_IOAPI -> CIO_BOTTOM
 *
 * This is always the cascade used for the meta-layers, but the actual (real)
 * implemented layers used by each meta-layer will depend on the compile-time
 * options used. It may be that some meta-layers do not consist of any real
 * layers, or in fact consist of several real layers.
 *
 */

/*
 * CIO_BOTTOM meta-layer
 *
 * do we use stadard libc file access, or POSIX extensions?
 * POSIX is the default unless LIBC_IO is set
 *
 */

#if defined(LIBC_IO)
#define CIO_BOTTOM CIO_LIBC
#else
#define CIO_BOTTOM CIO_UNIX
#endif

/*
 * CIO_IOAPI meta-layer
 *
 * do we use a file system API or filesystem specific tuning layer?
 * No API is the default, but the following are available: -
 *   - Lustre API (llapi), if LUSTRE_IO is set
 *   - Write Black Hole (used for debugging and timings excluding disk
 *     performance), if BLACKHOLE_IO is set
 *   - Rate limiting throttle, if THROTTLE_IO is set
 *
 */

#if defined(BLACKHOLE_IO)
#define CIO_IOAPI CIO_BLACKHOLE
#elif defined(LUSTRE_IO)
#define CIO_IOAPI CIO_LUSTRE_API
#elif defined(THROTTLE_IO)
#define CIO_IOAPI CIO_THROTTLE
#endif

#if !defined(CIO_IOAPI)
#define CIO_IOAPI CIO_BOTTOM
#endif

/*
 * CIO_BUFFER meta-layer
 *
 * do we use buffering?
 * the default is no, unless BUFRD_IO is set
 *
 */
#if defined(BUFRD_IO)
#define CIO_BUFFER CIO_WBUFFERING
#else
#define CIO_BUFFER CIO_IOAPI
#endif

/*
 * CIO_INFO meta-layer
 *
 * do we want tracing info?
 * the default is no, unless TRACE_IO is set
 *
 */
#if defined(TRACE_IO)
#define CIO_INFO CIO_TRACE
#else
#define CIO_INFO CIO_BUFFER
#endif

/*
 * CIO_TOP meta-layer
 *
 * what is the top IO layer?
 * the default is CIO_TIMING
 *
 */
#define CIO_TOP CIO_TIMING

/*----------------------------------------------------------------------------*/
/* Next Layers                                                                */
/*----------------------------------------------------------------------------*/

/*  Next we have macros to allow plugin layers
 *  You must have defined CIO_THIS_LAYER in your code to one
 *  of the prescribed values (in c_io_layers.h) before you including this header
 *
 *  The purpose here is to redirect the IO from the
 *  file including this header, into the next selected
 *  layer. Polymorphism by force.
 */

/* deal with special cases */

/* CIO_RBUFFERING always follows CIO_WBUFFERING */
#if (CIO_THIS_LAYER == CIO_WBUFFERING)
/* CIO_NEXT_LAYER => CIO_RBUFFERING */
#define CIO_NEXT_LAYER CIO_RBUFFERING
#endif

/* CIO_BLACKHOLE is followed by CIO_LUSTRE_API or CIO_THROTTLE if enabled */
#if (CIO_THIS_LAYER == CIO_BLACKHOLE)
#if defined(LUSTRE_IO)
/* CIO_NEXT_LAYER => CIO_LUSTRE_API */
#define CIO_NEXT_LAYER CIO_LUSTRE_API
#elif defined(THROTTLE_IO)
/* CIO_NEXT_LAYER => CIO_THROTTLE */
#define CIO_NEXT_LAYER CIO_THROTTLE
#endif
#endif

/* CIO_LUSTRE_API is followed by CIO_THROTTLE if enabled */
#if (CIO_THIS_LAYER == CIO_LUSTRE_API)
#if defined(THROTTLE_IO)
/* CIO_NEXT_LAYER => CIO_THROTTLE */
#define CIO_NEXT_LAYER CIO_THROTTLE
#endif
#endif

/* CIO_TIMING is followed by CIO_BYTESWAP if enabled */
#if (CIO_THIS_LAYER == CIO_TIMING)
#if defined(BYTESWAP_IO)
/* CIO_NEXT_LAYER => CIO_BYTESWAP */
#define CIO_NEXT_LAYER CIO_BYTESWAP
#endif
#endif

/* If not a special case, move on to the next meta-layer */

#if !defined(CIO_NEXT_LAYER)
/* define as next meta-layer */

/* The expansion of the META_LAYER_LOOKUP() macro has some unreachable code,
 * corresponding to the non-selected layers. This is by design, and is part of
 * the compile-time selection method.
 *
 * However, some compilers throw warnings for this behaviour, as it could be an
 * indicator of erroneous behaviour and/or a programming mistake. We will use
 * the um_compile_diag_suspend mechanism to suppress the warnings in this case,
 * as the code is correct.
 */

#if defined(__clang__)

#define LAYER_DIAG_OFF um_compile_diag_suspend(-Wunreachable-code)
#define LAYER_DIAG_ON  um_compile_diag_resume

#else

#define LAYER_DIAG_OFF
#define LAYER_DIAG_ON

#endif

#define CIO_NEXT_LAYER LAYER_DIAG_OFF                                          \
                       META_LAYER_LOOKUP(GET_META_LAYER(CIO_THIS_LAYER)+1)     \
                       LAYER_DIAG_ON

#endif

/*----------------------------------------------------------------------------*/
/* Layers Proceedure Macros                                                   */
/*----------------------------------------------------------------------------*/

#if (CIO_NEXT_LAYER == CIO_RBUFFERING)

/* Use CIO_RBUFFERING as the next layer */
#include "c_io_rbuffering.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_rbuffering_init(params)

#define C_IO_FINI                                                              \
        c_io_rbuffering_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_rbuffering_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_rbuffering_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_rbuffering_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_rbuffering_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_rbuffering_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_rbuffering_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_rbuffering_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_rbuffering_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_rbuffering_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_WBUFFERING)

/* Use CIO_WBUFFERING as the next layer */
#include "c_io_wbuffering.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_wbuffering_init(params)

#define C_IO_FINI                                                              \
        c_io_wbuffering_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_wbuffering_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_wbuffering_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_wbuffering_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_wbuffering_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_wbuffering_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_wbuffering_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_wbuffering_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_wbuffering_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_wbuffering_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_TIMING)

/* Use CIO_TIMING as the next layer */
#include "c_io_timing.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_timing_init(params)

#define C_IO_FINI                                                              \
        c_io_timing_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_timing_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_timing_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_timing_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_timing_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_timing_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_timing_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_timing_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_timing_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_timing_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_TRACE)

/* Use CIO_TRACE as the next layer */
#include "c_io_trace.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_trace_init(params)

#define C_IO_FINI                                                              \
        c_io_trace_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_trace_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_trace_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_trace_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_trace_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_trace_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_trace_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_trace_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_trace_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_trace_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_BYTESWAP)

/* Use CIO_BYTESWAP as the next layer */
#include "c_io_byteswap.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_byteswap_init(params)

#define C_IO_FINI                                                              \
        c_io_byteswap_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_byteswap_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_byteswap_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_byteswap_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_byteswap_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_byteswap_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_byteswap_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_byteswap_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_byteswap_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_byteswap_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_LIBC)

/* Use CIO_LIBC as the next layer */
#include "c_io_libc.h"

#include <stdio.h>

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_dummyInit(params);

#define C_IO_FINI()                                                            \
        c_io_isSuccessful();

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_libc_in((unit),(array),(len),(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_libc_out((unit),(array),(len),(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_libc_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_libc_close((unit))

#define C_IO_GETPOS(unit)                                                      \
        ftello(c_io_attributes[ (unit) ].fileHandle)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_libc_setpos((unit), (pos) ,(word_len))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_libc_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_libc_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_libc_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if ( CIO_NEXT_LAYER == CIO_UNIX )

/* Use CIO_UNIX as the next layer */
#include "c_io_unix.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_dummyInit((void *)(params))

#define C_IO_FINI()                                                            \
        c_io_isSuccessful()

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_unix_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_unix_close(unit)

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_unix_in((unit),(array),(len),(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_unix_out((unit),(array),(len),(word_len))

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_unix_setpos((unit), (pos) ,(word_len))

#define C_IO_GETPOS(unit)                                                      \
        lseek(*(int *)c_io_attributes[(unit)].fileHandle, (off_t)(0), SEEK_CUR)

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_unix_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_unix_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_unix_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_LUSTRE_API)

/* Use CIO_LUSTRE_API as the next layer */
#include "c_io_lustreapi.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_lustreapi_init(params)

#define C_IO_FINI                                                              \
        c_io_lustreapi_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_lustreapi_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_lustreapi_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_lustreapi_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_lustreapi_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_lustreapi_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_lustreapi_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_lustreapi_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_lustreapi_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_lustreapi_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_BLACKHOLE)

/* Use CIO_BLACKHOLE as the next layer */
#include "c_io_blackhole.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_blackhole_init(params)

#define C_IO_FINI                                                              \
        c_io_blackhole_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_blackhole_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_blackhole_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_blackhole_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_blackhole_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_blackhole_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_blackhole_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_blackhole_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_blackhole_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_blackhole_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#if (CIO_NEXT_LAYER == CIO_THROTTLE)

/* Use CIO_THROTTLE as the next layer */
#include "c_io_throttle.h"

/* layer macros */
#define C_IO_INIT(params)                                                      \
        c_io_throttle_init(params)

#define C_IO_FINI                                                              \
        c_io_throttle_fini

#define C_IO_IN(unit,array,len,word_len)                                       \
        c_io_throttle_in((unit), (array) ,(len) ,(word_len))

#define C_IO_OUT(unit,array,len,word_len)                                      \
        c_io_throttle_out((unit), (array) ,(len) ,(word_len))

#define C_IO_OPEN(unit,filename,filestatus,filemode)                           \
        c_io_throttle_open((unit), (filename) ,(filestatus) ,(filemode))

#define C_IO_CLOSE(unit)                                                       \
        c_io_throttle_close(unit)

#define C_IO_SETPOS(unit,pos,word_len)                                         \
        c_io_throttle_setpos((unit),(pos),(word_len))

#define C_IO_GETPOS(unit)                                                      \
        c_io_throttle_getpos((unit))

#define C_IO_CHANGE_MODE(unit,newMode)                                         \
        c_io_throttle_change_mode((unit),(newMode))

#define C_IO_SYNC(unit)                                                        \
        c_io_throttle_sync((unit))

#define C_IO_QUERY(unit, query, in_pack, out_pack)                             \
        c_io_throttle_query((unit), (query), (in_pack), (out_pack))

#endif

/* -------------------------------------------------------------------------- */

#endif
