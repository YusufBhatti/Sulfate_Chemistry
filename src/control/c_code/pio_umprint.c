/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_fort2c_prototypes.h"
#include "pio_umprint.h"
#include <stdlib.h>

/* C interface for umPrint */
void c_umPrint(umprint_options options)
{
  int setsum=0;

  /* first, encode the options which have been changed from their defaults */

  if(options.level!=0)
  {
    setsum|=level_set();
  }

  if(options.pe!=-1)
  {
    setsum|=pe_set();
  }

  if(options.model!=NULL)
  {
    setsum|=model_set();
  }

  if(options.UsrPrefix!=NULL)
  {
    setsum|=UsrPrefix_set();
  }

  if(options.HangIndent!=NULL)
  {
    setsum|=HangIndent_set();
  }

  /* call the Fortran interface for umPrint */

  f_umPrint(options.line,
            options.src,
            &options.level,
            &options.pe,
            options.model,
            options.UsrPrefix,
            options.HangIndent,
            &options.stdErrorToo,
            &setsum);
}
