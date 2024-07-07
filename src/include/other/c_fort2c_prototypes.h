#if !defined(C_FORT2C_PROTOTYPES_H)
#define C_FORT2C_PROTOTYPES_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>

/* Description:                                                               */
/*                                                                            */
/* Definitions needed for Fortran to C Interface                              */
/*                                                                            */
/* Information:                                                               */
/*                                                                            */
/* Provides definitions needed for Fortran to C interface.                    */

/*----------------------------------------------------------------------------*/
/* Define the functions that output the text string for messages.             */
/*----------------------------------------------------------------------------*/

/*
 *  Note that the trailing newline character is to be supplied by the print
 *  routine, and is no longer in the string. Similarly, leading new lines are
 *  now handled by the print routine - typically a newline is inserted for each
 *  change of unit.
 */

/* DEPENDS ON: fort2c_umprint_interfaces.o */
/* DEPENDS ON: ereport_mod.o */

extern int c_maxLineLen(void);

/* NB: MAX_OUTSTR is c_maxLineLen+1 to allow for the training NULL character. */
#define MAX_OUTSTR (size_t)(c_maxLineLen()+1)

extern void f_ereport(const char *, int64_t, const char *);
#define ereport(a,b,c) f_ereport(a,b,c)

extern void f_umPrint(const char * const,
                      const char * const,
                      const int  * const,
                      const int  * const,
                      const char * const,
                      const char * const,
                      const char * const,
                      const bool * const,
                      const int  * const);

/*----------------------------------------------------------------------------*/
/* Define other prototypes                                                    */
/*----------------------------------------------------------------------------*/

/* f_get_file_striping() is defined in lustre_control_mod.F90 */
void f_get_file_striping(int64_t *,
                         int64_t *,
                         int64_t *,
                         int64_t *,
                         const char *);

#endif
