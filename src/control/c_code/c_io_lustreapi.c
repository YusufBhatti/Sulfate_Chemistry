#if defined(C95_2B) && defined(LUSTRE_IO)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_LUSTRE_API

#include "c_io_nextlayer.h"

#include "c_io_internal.h"
#include "c_fort2c_prototypes.h"
#include "c_io_lustreapi.h"

#include <lustre/lustreapi.h>
#include <fcntl.h>

/* Functions in this file:

   Layer interfaces:
     init
     fini
     open
     close
     in
     out
     getpos
     setpos
     change_mode
     sync
     query

   Implementation routines:
     get_stripe
     apply_stripe

   To make a new implementation of an IO layer, copy this
   file, replacing lustreapi with your identifier.
*/

/*----------------------------------------------------------------------------*/
/* Typedefs                                                                   */
/*----------------------------------------------------------------------------*/

typedef struct {
  int64_t count;
  int64_t size;
  int64_t offset;
  int64_t pattern;
} lustrestripe;

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

static lustrestripe get_stripe   (char *filename);

static bool         apply_stripe (int64_t unit,
                                  lustrestripe *stripe_info,
                                  char *filename);

/*----------------------------------------------------------------------------*/
/* Procedures                                                                 */
/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_init(c_io_init_t *init_params)
{
  c_io_printf(-1,"init: lustreapi: initialised");
  c_io_init_next_layer(CIO_NEXT_LAYER);
  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_fini(void)
{
  return C_IO_FINI();
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_in(int64_t unit, char *array, int64_t len,
                          int64_t word_len)
{
    return C_IO_IN(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_out(int64_t unit, char *array, int64_t len,
                           int64_t word_len)
{
  return C_IO_OUT(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_setpos (int64_t unit, int64_t pos, int64_t word_len)
{
    return C_IO_SETPOS(unit,pos,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_getpos (int64_t unit)
{
  return C_IO_GETPOS(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_open(int64_t unit, char *filename, int64_t fileStatus,
                            int64_t fileMode)
{
  /*
   * if we are opening a new file, check for striping information,
   * else open normaly
   */

  int64_t updatedStatus = fileStatus;

  if (fileStatus==newFile)
  {
    lustrestripe stripe_info;

    stripe_info=get_stripe(filename);

    if(apply_stripe(unit,&stripe_info,filename))
    {
      /* after striping is applied, the file exists, so is now "old" */
      updatedStatus = oldFile;
    }

  }

  return C_IO_OPEN(unit,filename,updatedStatus,fileMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_close(int64_t unit)
{
  return C_IO_CLOSE(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_change_mode(int64_t unit, int64_t newMode)
{
  /* newMode is unused in some bottom layers, so avoid unused variable */
  (void)newMode;

  return C_IO_CHANGE_MODE(unit,newMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_sync(int64_t unit)
{
  return C_IO_SYNC(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_lustreapi_query(int64_t unit,
                             query_type_t query,
                             void *in_pack,
                             void *out_pack)
{
  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*----------------------------------------------------------------------------*/

static lustrestripe get_stripe(char *filename)
{
  lustrestripe result;

  /* query the loaded lustre settings in fortran */
  f_get_file_striping(&result.count,
                      &result.size,
                      &result.offset,
                      &result.pattern,
                      filename);

  return result;
}

/*----------------------------------------------------------------------------*/

static bool apply_stripe(int64_t unit,
                         lustrestripe *stripe_info,
                         char *filename)
{
  /*
   * We cannot use llapi_file_create(), as this would use the wrong file
   * permissions. The following method is functionally equivalent to a call
   * to llapi_file_create()
   */

  bool success=0;
  int open_return;

  open_return = llapi_file_open(filename,
                                O_CREAT|O_WRONLY,
                                S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH,
                                (unsigned long long)stripe_info->size,
                                (int)stripe_info->offset,
                                (int)stripe_info->count,
                                (int)stripe_info->pattern);

  if (open_return < 0 )
  {
    /* Something went wrong. */
    c_io_printf(unit, "open: lustreapi: ERROR: %s", strerror(-open_return));
  }
  else
  {
    /* Striping information was set. */
    close(open_return);
    success=1;
  }

  return success;
}

#endif
