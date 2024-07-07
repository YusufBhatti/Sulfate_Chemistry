#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_BYTESWAP

#include "c_io_nextlayer.h"

#include "c_io_internal.h"

#include "c_fort2c_prototypes.h"
#include "c_io_byteswap.h"
#include "pio_byteswap.h"
#include "c_fort2c_types.h"
#include <inttypes.h>

/*  Functions in this file:
 *
 *  Layer interfaces:
 *    init          : Determine endianness and pass through
 *    fini          : pass through
 *    open          : pass through
 *    close         : pass through
 *    in            : pass through, then byteswap
 *    out           : byteswap, then pass through
 *    getpos        : pass through
 *    setpos        : pass through
 *    change_mode   : pass through
 *    sync          : pass through
 *    query         : pass through
 *
 *  Implentation auxiliary routines:
 *    c_io_shouldByteSwap     : report whether a byteswap is needed
 *    c_io_byteswap           : byteswap a supplied array with specified
 *                              word length and element count
 *
 */

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

bool    c_io_shouldByteSwap(void);
void    c_io_byteswap(char *, int64_t, int64_t);

static endianness machineEndian;
static const endianness fileEndian=bigEndian;

/*----------------------------------------------------------------------------*/

int64_t c_io_byteswap_init(c_io_init_t *init_params)
{
  machineEndian=get_machine_endianism();

  c_io_init_next_layer(CIO_NEXT_LAYER);
  if (machineEndian==bigEndian){
    c_io_printf(-1,"init: byteswap: Machine is big endian");
  } else
    c_io_printf(-1,"init: byteswap: Machine is little endian");

  if (c_io_shouldByteSwap())
    c_io_printf(-1,"init: byteswap: IO will byteswap on this machine");
  return C_IO_INIT(init_params);
}

int64_t c_io_byteswap_fini(void)
{
  return C_IO_FINI();
}

int64_t c_io_byteswap_open   (int64_t unit,
                              char *filename,
                              int64_t fileStatus,
                              int64_t fileMode)
{
  return C_IO_OPEN(unit,filename,fileStatus,fileMode);
}

int64_t c_io_byteswap_close  (int64_t unit)
{
  return C_IO_CLOSE(unit);
}


int64_t c_io_byteswap_in     (int64_t unit,
                              char *array,
                              int64_t len,
                              int64_t word_len)
{
  int64_t ammountRead;
  ammountRead=C_IO_IN(unit,array,len,word_len);
  c_io_byteswap(array,ammountRead,word_len);
  return ammountRead;
}

int64_t c_io_byteswap_out    (int64_t unit,
                              char *array,
                              int64_t len,
                              int64_t word_len)
{
  int64_t ammountWritten;
  c_io_byteswap(array,len,word_len);
  ammountWritten=C_IO_OUT(unit,array,len,word_len);
  c_io_byteswap(array,len,word_len);
  return ammountWritten;
}

int64_t c_io_byteswap_setpos (int64_t unit,
                              int64_t pos,
                              int64_t word_len)
{
  return C_IO_SETPOS(unit,pos,word_len);
}

int64_t c_io_byteswap_getpos (int64_t unit)
{
  return C_IO_GETPOS(unit);
}

int64_t c_io_byteswap_change_mode(int64_t unit,int64_t newMode)
{
  //Macro expansion of C_IO_CHANGE_MODE may remove newMode
  //Therefore, avoid unused variable
  (void)(newMode);

  return C_IO_CHANGE_MODE(unit,newMode);
}

int64_t c_io_byteswap_sync(int64_t unit)
{
  return C_IO_SYNC(unit);
}

int64_t c_io_byteswap_query(int64_t unit,
                            query_type_t query,
                            void *in_pack,
                            void *out_pack)
{
  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*
 *
 * local routines below, should only be accessed by this file
 *
 */

bool c_io_shouldByteSwap(void)
{
  return (fileEndian!=machineEndian);
}

void c_io_byteswap(char *array, int64_t len, int64_t word_len)
{

  /* There is no need to byteswap single byte data */
  if (word_len == WORD8BYTES)
  {
    return;
  }

  if (!c_io_shouldByteSwap()) return;

  if (pio_byteswap(array,len,word_len) != 0)
  {
    c_io_printf(-1,
                "Byteswap: Don't know how to swap a word size of %" PRId64,
                word_len);
    ereport("c_io_byteswap",1,"Unsupported word length");
  }

  return;
}

#endif
