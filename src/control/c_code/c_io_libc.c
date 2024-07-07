#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_LIBC

#include "c_io_nextlayer.h"

#include "c_io_internal.h"
#include "c_io_libc.h"
#include "c_io.h"

#include <stdio.h>
#include <libgen.h>
#include <fcntl.h>

/*  Functions in this file:
 *
 *  Layer interfaces:
 *    open
 *    close
 *    setpos
 *    in
 *    out
 *    change_mode
 *    sync
 *    query
 *
 *  Layer interface implemented as macros (not in this file):
 *    init
 *    fini
 *    getpos
 *
 */

/*---------------------------------------------------------------------------*/
/* Prototypes                                                                */
/*---------------------------------------------------------------------------*/

static int64_t call_fsync(int);

/*---------------------------------------------------------------------------*/

/*
 *
 *  Implentation routines
 *
 */

int64_t c_io_libc_open(int64_t unit,
                       char *filename,
                       int64_t fileStatus,
                       int64_t fileMode)
{
  int dir_fd=-1;

  if (fileStatus==newFile)
  {
    char *dir_name = NULL;
    char *filename2 = NULL;

    /* calls to dirname() mangle the original string, so we must create a
     * copy first so we can retain the original filename.
     */
    filename2 = malloc(strlen(filename)+1);
    snprintf(filename2,strlen(filename)+1,"%s",filename);
    dir_name = dirname(filename2);
    dir_fd = open(dir_name, O_RDONLY);

    if (dir_fd==-1)
    {
      perror("open(2) [dir] :");
      free(filename2);
      c_io_printf(unit, "c_io_libc: open: open call failed on directory "
                  "for unit %" PRId64, unit);
      return c_io_fail;
    }

    free(filename2);
  }

  c_io_attributes[(unit)].fileHandle=
    fopen(filename,modeText[fileStatus][fileMode]);

  if (c_io_attributes[(unit)].fileHandle==NULL)
  {
    perror("open(2) [file] :");
    c_io_printf(unit, "c_io_libc: open: open call failed on file "
                "for unit %" PRId64, unit);
    return c_io_fail;
  }

  if (fileStatus==newFile)
  {
    /* fsync directory:
     *   Whilst the file may be committed to disk, its existance in the
     *   directory may not be. Force a directory inode update with fsync()
     */
    int64_t fsync_ret = call_fsync(dir_fd);

    close(dir_fd);

    if (fsync_ret == c_io_fail)
    {
       c_io_printf(unit, "c_io_libc: open: fsync call failed on directory "
                   "for unit %" PRId64, unit);
    }

    return fsync_ret;
  }

  return c_io_success;
}


int64_t c_io_libc_setpos  (int64_t unit,int64_t pos,int64_t word_len)
{
  int64_t retVal;

  if (pos==c_io_file_end)
  {
    retVal = fseeko(c_io_attributes[unit].fileHandle, (off_t)(0), SEEK_END);
    c_io_attributes[unit].end_of_file=true;
  }
  else
  {
    retVal = fseeko(c_io_attributes[unit].fileHandle,
                    (off_t)(pos*word_len), SEEK_SET);
    c_io_attributes[unit].end_of_file=false;
  }

  if (retVal==0)
  {
    return c_io_success;
  }
  else
  {
    c_io_printf(unit, "c_io_libc: setpos: fseeko failed with return code %"
                PRId64 ", when seeking to position %" PRId64 " of unit %"
                PRId64, retVal, pos, unit);
    return c_io_fail;
  }
}

int64_t c_io_libc_close  (int64_t unit)
{
  /* If we have written to file, ensure the contents are synced to disk */
  if (c_io_attributes[unit].modified_on_disk)
  {
    c_io_printf(unit,"c_io_libc: close: syncing file pointer fp=%" PRIdPTR,
                (uintptr_t)c_io_attributes[ (unit) ].fileHandle);
    if (c_io_libc_sync(unit))
    {
      /* fsync failed */
      return c_io_fail;
    }

  }

  fclose(c_io_attributes[(unit)].fileHandle);
  return c_io_success;
}

/*
 * Calls to fread() return the number elements transfered. This is less than
 * that requested only if a read error or end-of-file (EOF) occurs. We
 * therefore need to test on the length of the return to trap errors and EOF.
 */
int64_t c_io_libc_in (int64_t unit,
                      char    *array,
                      int64_t len,
                      int64_t word_len)
{
  const c_io_attributes_t *attrib=&c_io_attributes[unit];
  int64_t elements_transfered;

  c_io_attributes[unit].end_of_file=false;

  elements_transfered = (int64_t)fread(array, (size_t)word_len,
                                       (size_t)len, attrib->fileHandle);

  if(elements_transfered<len)
  {
    /* Either EOF or Error; so check end-of-file marker */
    if(feof(attrib->fileHandle))
    {
      c_io_printf(unit, "c_io_libc: in: %" PRId64
                  " bytes were requested, but only %" PRId64 " were transfered"
                  " - we have reached end-of-file : fp=%" PRIdPTR,
                  len, elements_transfered, (uintptr_t)attrib->fileHandle);
      c_io_attributes[unit].end_of_file=true;
    }
    else
    {
      c_io_printf(unit, "c_io_libc: ERROR: "
                  "Error in reading file : unit = %" PRId64 ", fp=%" PRIdPTR,
                  unit, (uintptr_t)attrib->fileHandle);
    }

  }

  return elements_transfered;
}

/*
 * Calls to fwrite() return the number elements transfered. This is less than
 * that requested only if a write error occurs. We therefore need to test on
 * the length of the return to trap errors.
 */
int64_t c_io_libc_out (int64_t unit,
                       char    *array,
                       int64_t len,
                       int64_t word_len)
{
  const c_io_attributes_t *attrib=&c_io_attributes[unit];
  int64_t elements_transfered;

  elements_transfered = (int64_t)fwrite(array, (size_t)word_len,
                                        (size_t)len, attrib->fileHandle);

  if(elements_transfered<len)
  {
    /* There was an error on write */
    c_io_printf(unit,"c_io_libc: ERROR: "
                "Error in writing file : unit = %" PRId64 ", fp=%" PRIuPTR,
                unit, (uintptr_t)attrib->fileHandle);
  }

  return elements_transfered;
}

int64_t c_io_libc_change_mode(int64_t unit,int64_t newMode)
{
  int returnval;

  /* avoid unused variable */
  (void)(newMode);

  returnval=fflush(c_io_attributes[unit].fileHandle);

  if(returnval==0)
  {
    return c_io_success;
  }
  else
  {
    c_io_printf(unit, "c_io_libc: change_mode: fflush failed with return code "
                "%d, when flushing unit %" PRId64, returnval, unit);
    return c_io_fail;
  }

}

int64_t c_io_libc_sync(int64_t unit)
{
  int returnval;

  /* we must fflush the unit fileHandle before we fsync, or some data will be
   * held in the application stream buffer (and therefore miss the fsync).
   */
  returnval = fflush(c_io_attributes[unit].fileHandle);

  if (returnval != 0)
  {
    c_io_printf(unit, "c_io_libc: c_io_libc_sync: fflush failed with return code "
                "%d, when flushing unit %" PRId64, returnval, unit);
    return c_io_fail;
  }

  return call_fsync(fileno(c_io_attributes[unit].fileHandle));
}

int64_t c_io_libc_query(int64_t unit,
                        query_type_t query,
                        void *in_pack,
                        void *out_pack)
{
  switch(query)
  {
    case qry_file_exist:
      if (access(((in_pack_fname_t *)in_pack)->filename,F_OK)==0)
      {
        ((out_pack_bool_t *)out_pack)->res = true;
      }
      else
      {
        ((out_pack_bool_t *)out_pack)->res = false;
      }
      break;
  }

  /* avoid unused variable */
  (void)(unit);

  return c_io_success;
}

/*---------------------------------------------------------------------------*/
/* Local Routines                                                            */
/*---------------------------------------------------------------------------*/

static int64_t call_fsync(int fp)
{
  int fsync_return;

  fsync_return=fsync(fp);

  if(fsync_return!=0)
  {
    perror("c_io_libc_sync: ERROR: ");
    return c_io_fail;
  }

  return c_io_success;
}

/*---------------------------------------------------------------------------*/

#endif
