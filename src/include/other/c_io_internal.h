#if !defined(C_IO_INTERNAL_H)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#define C_IO_INTERNAL_H

#include "c_io.h"
#include <inttypes.h>

/* -------------------------------------------------------------------------- */
/*  Print function decorators                                                 */
/* -------------------------------------------------------------------------- */

/* The following decorators are used to specify additional attributes about
 * the function.
 *
 * They use a compiler extension supported by GCC and Clang. Other compilers
 * will use an empty decorator.
 *
 * They specify to the compiler that the function expects formatting arguments.
 * This allows additional checking to be performed on the correctness of those
 * arguments. For example the following would fail to compile with appropriate
 * warning messages with GCC and clang once the decorators are used: -
 *
 *    c_io_printf(1,"this expects an integer (%d)\n","but I am a string!");
 *
 */

/* the "common" decorator is used for c_io_printf_common() */

#if defined(__GNUC__) || defined(__clang__)
/* the C_IO_PRINTF_COMMON_DECORATOR specifies that the decorated function is a
 * printf() variant which expects a variadic format specifier.
 */
#define C_IO_PRINTF_COMMON_DECORATOR __attribute__                             \
                                     ((__format__ (__printf__, 2, 0)))
#endif

#if !defined(C_IO_PRINTF_COMMON_DECORATOR)
#define C_IO_PRINTF_COMMON_DECORATOR
#endif

/* the "base" decorator is used for c_io_printf() and c_io_printf_stderr() */

#if defined(__GNUC__) || defined(__clang__)
/* the C_IO_PRINTF_BASE_DECORATOR specifies that the decorated function is a
 * printf() variant which expects explicit format specifiers.
 */
#define C_IO_PRINTF_BASE_DECORATOR __attribute__                               \
                                   ((__format__ (__printf__, 2, 3)))
#endif

#if !defined(C_IO_PRINTF_BASE_DECORATOR)
#define C_IO_PRINTF_BASE_DECORATOR
#endif

/* -------------------------------------------------------------------------- */
/*  c_io_internal.h Exported Prototypes                                       */
/* -------------------------------------------------------------------------- */

extern void    c_io_printf(int64_t,
                           const char *fmt,
                           ...)
                           C_IO_PRINTF_BASE_DECORATOR;

extern void    c_io_printf_stderr(int64_t unit,
                                  const char *fmt,
                                  ...)
                                  C_IO_PRINTF_BASE_DECORATOR;

extern int64_t c_io_init(int64_t, int64_t, int64_t, int64_t,
                         int64_t, bool, int64_t, bool);
extern int64_t c_io_fini(void);
extern void    c_io_attachHelper(int64_t role);
extern void    c_io_detachAllHelpers(void);

#endif

