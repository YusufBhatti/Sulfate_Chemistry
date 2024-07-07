#if !defined(EXCEPTIONS_LINUX_H)
#define EXCEPTIONS_LINUX_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: exceptions-linux.o */

#include <ucontext.h>
#include "exceptions.h"

/*----------------------------------------------------------------------------*/
/* Prepocessor Defines                                                        */
/*----------------------------------------------------------------------------*/

#define xstr(s) #s
#define str(s) xstr(s)

#if   defined(__i386__)
#define INSTR_PTR eip
#elif defined(__x86_64__)
#define INSTR_PTR rip
#elif defined(__arm__)
#define INSTR_PTR arm_pc
#elif defined(__aarch64__)
#define INSTR_PTR pc
#else
#error - no specific handling of this architecture in exceptions-linux.h
#endif

#if defined(INSTR_PTR)
#define INSTR_PTR_STR str(INSTR_PTR)
#endif

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

extern void signal_do_backtrace_linux  (int, siginfo_t *, void *);
extern void signal_trap_extra_fp_linux (fpetraps_t *);

#endif
