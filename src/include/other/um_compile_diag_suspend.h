#if !defined(UM_COMPILE_DIAG_SUSPEND_H)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This header defines macros to turn off compile-time checking in a compiler
 * agnostic way.
 *
 * The macros are:
 *    um_compile_diag_global_suspend(flag)
 *    um_compile_diag_suspend(flag)
 *    um_compile_diag_resume
 *
 * um_compile_diag_global_suspend() is used to suspend a checking flag globally
 * if the compiler does not support scoped suspension. This must be called
 * outside of a function body. "flag" is the compiler flag to suspend.
 *
 * um_compile_diag_suspend() is used to suspend a checking flag in a scoped
 * region. "flag" is the compiler flag to suspend. um_compile_diag_resume is
 * to end the scoped region and resume compiler checking. These may be uesed
 * anywhere in the code, but pairs must be strictly nested.
 */

#define UM_COMPILE_DIAG_SUSPEND_H

#define UM_EXEC_PRAGMA(pragma_string) _Pragma(#pragma_string)
#define UM_EXEC_PRAGMA_STRINGIFY(string) #string

/* Deal with GCC */

#if defined(__GNUC__) && !defined(__clang__)

/* we are using GCC */

#if (__GNUC__ < 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ <8))

/* older versions of GCC (<4.8) cannot suspend checking, so we must turn it off
 * completely. Call um_compile_diag_global_suspend from the top of the unit,
 * after including this header.
 */
#define um_compile_diag_global_suspend(flag)                                   \
_Pragma("message \"GCC version too low to modify diagnostic checking\"")       \
UM_EXEC_PRAGMA(message UM_EXEC_PRAGMA_STRINGIFY(flag permanently turned off!)) \
UM_EXEC_PRAGMA(GCC diagnostic ignored #flag)

#else

/* newer versions of GCC can suspend checking, so we do that as it is preferable
 * to continue to check other code in the same file.
 */
#define um_compile_diag_suspend(flag) _Pragma("GCC diagnostic push")           \
                                    UM_EXEC_PRAGMA(GCC diagnostic ignored #flag)

#define um_compile_diag_resume _Pragma("GCC diagnostic pop")

#endif

#endif

/* Deal with Clang */

#if defined(__clang__)

/* we are using Clang - which supports suspending checking */

#define um_compile_diag_suspend(flag) _Pragma("clang diagnostic push")         \
                                  UM_EXEC_PRAGMA(clang diagnostic ignored #flag)

#define um_compile_diag_resume _Pragma("clang diagnostic pop")

#endif

/* Deal with other compilers */

/* if we have not specifically set the macros yet, set them blank */

#if !defined(um_compile_diag_suspend)
#define um_compile_diag_suspend(flag)
#endif

#if !defined(um_compile_diag_resume)
#define um_compile_diag_resume
#endif

#if !defined(um_compile_diag_global_suspend)
#define um_compile_diag_global_suspend(flag)
#endif

#endif
