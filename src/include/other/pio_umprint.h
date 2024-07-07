#if !defined(PIO_UMPRINT_H)
#define PIO_UMPRINT_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "um_compile_diag_suspend.h"
#include "c_fort2c_prototypes.h"
#include <stdbool.h>

/* We must suspend the -Woverride-init checking on older versions of GCC
 * [See later for more info]
 */
#if defined(__GNUC__)
um_compile_diag_global_suspend(-Woverride-init)
#endif

/* DEPENDS ON: pio_umprint.o */

typedef struct UMPRINT_OPTIONS
{
  const char *line;
  const char *src;
  const int  level;
  const int  pe;
  const char *model;
  const char *UsrPrefix;
  const char *HangIndent;
  const bool stdErrorToo;
  const char dummy;
}
umprint_options;

extern int level_set(void);
extern int pe_set(void);
extern int model_set(void);
extern int UsrPrefix_set(void);
extern int HangIndent_set(void);
extern void c_umPrint(umprint_options);

/* The umPrint implementation provides defaults for all possible arguments -
 * but then overrides them if they are provided by the user - through the use
 * of a compound literal
 *
 * However, this can cause warning messages to be omitted by compilers, as
 * re-initialisation is often an indicator of a programming error.
 *
 * Therefore, the um_compile_diag_suspend mechanism is used for suppressing
 * these messages in this case, as it is correct usage.
 */

#if defined(__GNUC__) && !defined(__clang__)

#define printdiag_suspend um_compile_diag_suspend(-Woverride-init)
#define printdiag_resume  um_compile_diag_resume

#elif defined(__clang__)

#define printdiag_suspend um_compile_diag_suspend(-Winitializer-overrides)
#define printdiag_resume  um_compile_diag_resume

#else

#define printdiag_suspend
#define printdiag_resume

#endif

/* define umPrint()
 * As the C standard requires at least one argument to be passed into variadic
 * macros, we will indirect this by one level of expansion.
 * (only line and src are required arguments, all the rest are optional)
 */

#define umPrint(a,...) XumPrint(a,__VA_ARGS__,.dummy='x')

/* define XumPrint()
 * note: this definition contains a single iteration do{}while(0) loop.
 * This is done to protect the macro expansion in all cases (commonly known as
 * "swallowing the semicolon").
 * e.g. if() umPrint(); else {} is a valid calling method
 */
#define XumPrint(a,b,...) do{                                                  \
                          printdiag_suspend                                    \
                          c_umPrint((umprint_options){.line=(a),               \
                                                      .src=(b),                \
                                                      .level=0,                \
                                                      .pe=-1,                  \
                                                      .model=NULL,             \
                                                      .UsrPrefix=NULL,         \
                                                      .HangIndent=NULL,        \
                                                      .stdErrorToo=false,      \
                                                      __VA_ARGS__});           \
                          printdiag_resume                                     \
                          }while(0)

#endif
