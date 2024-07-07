#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_BLACKHOLE

#include "c_io_nextlayer.h"

#include "c_io_internal.h"
#include "c_io_blackhole.h"

/* This unit allows the "blackholing" of I/O writes. When this is enabled, any
 * I/O operation which would modify the contents of the disk is discarded. This
 * allows the emulation of the complete I/O stack, but eliminates the time lost
 * waiting for the disk.
 * This may be useful for debugging or benchmarking code on systems with poor
 * I/O performance.
 */

/*  Functions in this file:
 *
 *  Layer interfaces:
 *    out
 *    open
 *    sync
 *    init
 *    fini
 *    change_mode
 *    getpos
 *    close
 *    setpos
 *    in
 *    query
 *
 */

#define BLACKHOLE_HANDLE -1

/*----------------------------------------------------------------------------*/
/* Implentation routines                                                      */
/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_init(c_io_init_t *init_params)
{
  c_io_init_next_layer(CIO_NEXT_LAYER);

  c_io_printf(-1,
              "init: blackhole: File creation, writes, and syncs "
              "will have no effect.");

  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_open(int64_t unit,
                            char    *filename,
                            int64_t fileStatus,
                            int64_t fileMode)
{
  /*
   * If we are creating a new file, blackhole here to prevent its creation on
   * disk.
   */
  if (fileStatus==newFile)
  {
    c_io_attributes[unit].fileHandle = malloc(sizeof(int));
    *(int *)c_io_attributes[unit].fileHandle = BLACKHOLE_HANDLE;
    return c_io_success;
  }

  return C_IO_OPEN(unit,filename,fileStatus,fileMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_out (int64_t unit,
                            char    *array,
                            int64_t len,
                            int64_t word_len)
{
  /* out calls are blackholed - return without continuing down the IO stack */

  (void) unit;
  (void) array;
  (void) word_len;

  return len;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_sync(int64_t unit)
{
  /* sync calls are blackholed - return without continuing down the IO stack */

  (void) unit;

  return c_io_success;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_in (int64_t unit,
                           char    *array,
                           int64_t len,
                           int64_t word_len)
{
  if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
  {

    return c_io_success;
  }

  return C_IO_IN(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_setpos (int64_t unit,
                               int64_t pos,
                               int64_t word_len)
{

  if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
  {
    return c_io_success;
  }

  return C_IO_SETPOS(unit,pos,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_getpos (int64_t unit)
{
  if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
  {
    return 0;
  }

  return C_IO_GETPOS(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_fini(void)
{
  return C_IO_FINI();
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_change_mode (int64_t unit,
                                    int64_t newMode)
{
  if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
  {
    return c_io_success;
  }

  return C_IO_CHANGE_MODE(unit,newMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_query (int64_t      unit,
                              query_type_t query,
                              void         *in_pack,
                              void         *out_pack)
{
  /* Queries can happen on unopened units, so we must check if we have issued
   * a handle before checking if it is the BLACKHOLE_HANDLE pseudo-handle.
   */
  if (c_io_attributes[unit].fileHandle != NULL)
  {
    if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
    {
      return c_io_success;
    }
  }

  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_blackhole_close (int64_t unit)
{
  if (*(int *)c_io_attributes[unit].fileHandle == BLACKHOLE_HANDLE)
  {
    free(c_io_attributes[unit].fileHandle);
    return c_io_success;
  }

  return C_IO_CLOSE(unit);
}

/*----------------------------------------------------------------------------*/

#endif
