#if !defined(EXCEPTIONS_IBM_H)
#define EXCEPTIONS_IBM_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: exceptions-ibm.o */

#include <ucontext.h>
#include "exceptions.h"

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

extern void xl__trce                 (int, siginfo_t *, void *);
extern void xl__trcedump             (int, siginfo_t *, void *);
extern void signal_do_backtrace_ibm  (int, siginfo_t *, void *);
extern void signal_do_core_ibm       (int, siginfo_t *, void *);
extern void signal_trap_extra_fp_ibm (fpetraps_t *);

#endif
