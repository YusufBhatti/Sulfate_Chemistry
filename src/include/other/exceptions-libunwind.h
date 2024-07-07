#if !defined(EXCEPTIONS_LIBUNWIND_H)
#define EXCEPTIONS_LIBUNWIND_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: exceptions-libunwind.o */

#include <ucontext.h>

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

void signal_do_backtrace_libunwind(int, siginfo_t *, void *);

#endif
