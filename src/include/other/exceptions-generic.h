#if !defined(EXCEPTIONS_GENERIC_H)
#define EXCEPTIONS_GENERIC_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: exceptions-generic.o */

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

extern void signal_do_backtrace_generic (int, siginfo_t *, void *);
extern void signal_do_core_generic      (int, siginfo_t *, void *);

#endif
