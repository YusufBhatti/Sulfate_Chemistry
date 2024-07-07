#if !defined(PIO_BYTESWAP_OPT)
#define PIO_BYTESWAP_OPT

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This header contains all the macros and pre-processing keys to select the  */
/* correct optimised version of the byteswapping code for your machine.       */

/* -------------------------------------------------------------------------- */
/* bswap() family of macros                                                   */
/* -------------------------------------------------------------------------- */

/* Scalar Macros */

/*
 * There are several possible implementations of scalar bswap() we may use.
 *
 * i) If we are using a compiler which implements GNU GCC extensions
 *    (defines __GNUC__ >= 2), but is eariler than version 4.3, we can use the
 *    non-standard <byteswap.h> implementations, which may expand the macro to
 *    optimised inline assembly [dependant on architecture and version].
 *
 * ii) If we are using a compiler which implements GNU GCC extensions at version
 *     4.3 or later, but before version 4.8, we can use the the non-standard
 *     intrinsic __builtin_bswap*() except for bswap_16, which we implement
 *     directly as a bit-shift based macro.
 *
 * iii) If we are using a compiler which implements GNU GCC extensions at
 *      version 4.8 or later, we can use the the non-standard intrinsic
 *      __builtin_bswap*() functions for all sizes of bswap.
 *
 * iv) If we are using the Cray compiler (defines _CRAYC) with GNU extensions
 *     enabled (-hgnu): we can use GCC extensions, but not inline assembly.
 *     Therefore, fall-back to the non-standard intrinsic __builtin_bswap*()
 *     functions, except for __builtin_bswap16, which isn't implemented, and
 *     so therefore use a directly implemented macro for bswap_16
 *
 * v) If none of the above are possible, implement the macros directly as
 *    bit-shift operations. This is much less efficient, but uses only
 *    standard C99, so will always work.
 *
 * vi) If the user defines one of the macros PIO_USE_BSWAP_BITSHIFT,
 *     PIO_USE_BSWAP_BUILTINS, or PIO_USE_BSWAP_BYTESWAP_H, override the above
 *     cases and use the direct bit-shift based macros, the __builtin_bswap*()
 *     non-standard intrinsic functions, or the <byteswap.h> headers
 *     repectively, as per that choice.
 *
 */

/* ensure test macros are unset */
#undef PIO_USE_BSWAP_USERCHOICE
#undef PIO_USE_BSWAP_64_BITSHIFT
#undef PIO_USE_BSWAP_32_BITSHIFT
#undef PIO_USE_BSWAP_16_BITSHIFT
#undef PIO_USE_BSWAP_BUILTINS_64
#undef PIO_USE_BSWAP_BUILTINS_32
#undef PIO_USE_BSWAP_BUILTINS_16

/* test cases */

/* test for user choice */

#if defined PIO_USE_BSWAP_BYTESWAP_H

/* user choice is <byteswap.h> */
#define PIO_USE_BSWAP_USERCHOICE

#elif defined PIO_USE_BSWAP_BUILTINS

/* user choice is __builtin_bswap*() */
#define PIO_USE_BSWAP_USERCHOICE
#define PIO_USE_BSWAP_BUILTINS_64
#define PIO_USE_BSWAP_BUILTINS_32
#define PIO_USE_BSWAP_BUILTINS_16

#elif defined PIO_USE_BSWAP_BITSHIFT

/* user choice is bit-shift macros */
#define PIO_USE_BSWAP_USERCHOICE
#define PIO_USE_BSWAP_64_BITSHIFT
#define PIO_USE_BSWAP_32_BITSHIFT
#define PIO_USE_BSWAP_16_BITSHIFT

#endif

/* test for used case */

#if !defined(PIO_USE_BSWAP_USERCHOICE)

#if defined __GNUC__ && __GNUC__ >= 2

#if defined __GNUC__ && __GNUC__       >= 4 \
                     && __GNUC_MINOR__ >= 3 \
                     && defined _CRAYC

/* case iv) */
#define PIO_USE_BSWAP_BUILTINS_64
#define PIO_USE_BSWAP_BUILTINS_32
#define PIO_USE_BSWAP_16_BITSHIFT

#elif defined __GNUC__ && __GNUC__       >= 4 \
                       && __GNUC_MINOR__ >= 8

/* case iii) */
#define PIO_USE_BSWAP_BUILTINS_64
#define PIO_USE_BSWAP_BUILTINS_32
#define PIO_USE_BSWAP_BUILTINS_16

#elif defined __GNUC__ && __GNUC__       >= 4 \
                       && __GNUC_MINOR__ >= 3

/* case ii) */
#define PIO_USE_BSWAP_BUILTINS_64
#define PIO_USE_BSWAP_BUILTINS_32
#define PIO_USE_BSWAP_16_BITSHIFT

#else

/* case i) */
#define PIO_USE_BSWAP_BYTESWAP_H

#endif

#else

/* case v) */
#define PIO_USE_BSWAP_64_BITSHIFT
#define PIO_USE_BSWAP_32_BITSHIFT
#define PIO_USE_BSWAP_16_BITSHIFT

#endif

#else

/* case vi) */

#endif

/* implement macros */

#if defined(PIO_USE_BSWAP_BUILTINS_64)
#define bswap_64(x) __builtin_bswap64(x)
#endif

#if defined(PIO_USE_BSWAP_BUILTINS_32)
#define bswap_32(x) __builtin_bswap32(x)
#endif

#if defined(PIO_USE_BSWAP_BUILTINS_16)
#define bswap_16(x) __builtin_bswap16(x)
#endif

#if defined(PIO_USE_BSWAP_BYTESWAP_H)
#include <byteswap.h>
#endif

#if defined(PIO_USE_BSWAP_64_BITSHIFT)

#define bswap_64(x) \
               ((((x) & UINT64_C(0xff00000000000000)) >> 56) \
              | (((x) & UINT64_C(0x00ff000000000000)) >> 40) \
              | (((x) & UINT64_C(0x0000ff0000000000)) >> 24) \
              | (((x) & UINT64_C(0x000000ff00000000)) >> 8)  \
              | (((x) & UINT64_C(0x00000000ff000000)) << 8)  \
              | (((x) & UINT64_C(0x0000000000ff0000)) << 24) \
              | (((x) & UINT64_C(0x000000000000ff00)) << 40) \
              | (((x) & UINT64_C(0x00000000000000ff)) << 56))

#endif

#if defined(PIO_USE_BSWAP_32_BITSHIFT)

#define bswap_32(x) \
               ((((x) & UINT32_C(0xff000000)) >> 24) \
              | (((x) & UINT32_C(0x00ff0000)) >>  8) \
              | (((x) & UINT32_C(0x0000ff00)) <<  8) \
              | (((x) & UINT32_C(0x000000ff)) << 24))

#endif

#if defined(PIO_USE_BSWAP_16_BITSHIFT)

#define bswap_16(x) \
               ((((x) & UINT16_C(0xff00)) >>  8) \
              | (((x) & UINT16_C(0x00ff)) <<  8))

#endif

#endif
