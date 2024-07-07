#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Compliance to Posix specification IEEE Std 1003.1-2001 or later is requred
 * by this unit in order to access some functionality (e.g. access() )
 */
#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_UNIX

#include "c_io_nextlayer.h"

#include "c_io_internal.h"
#include "c_io_unix.h"

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <libgen.h>

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
/* Implentation routines                                                     */
/*---------------------------------------------------------------------------*/

int64_t c_io_unix_open(int64_t unit,
                       char *filename,
                       int64_t fileStatus,
                       int64_t fileMode)
{
  int flags=0;
  int mode;
  int dir_fd = -1;

  if (fileStatus==newFile)
  {
    char *dir_name = NULL;
    char *filename2 = NULL;

    if (fileMode==readonly)
    {
      c_io_printf(unit,"c_io_unix: open:"
                  " ERROR: cannot create a new file in read-only mode.");
      return c_io_fail;
    }
    else
    {
      flags=flags | O_CREAT | O_TRUNC;
    }

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
      c_io_printf(unit, "c_io_unix: open: open call failed on directory "
                  "for unit %" PRId64, unit);
      free(filename2);
      return c_io_fail;
    }

    free(filename2);
  }

  mode=S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH;

  if (fileMode==readonly)
    flags=flags | O_RDONLY;
  else if (fileMode==readwrite)
    flags=flags | O_RDWR;
  else if (fileMode==writeonly)
    flags=flags | O_WRONLY;

  c_io_attributes[ (unit) ].fileHandle=malloc(sizeof(int));

  *(int *)c_io_attributes[unit].fileHandle = open(filename,flags,mode);

  if (*(int *)c_io_attributes[unit].fileHandle==-1)
  {
    perror("open(2) [file] :");
    c_io_printf(unit, "c_io_unix: open: open call failed on file "
                "for unit %" PRId64, unit);
    return c_io_fail;
  }

  if (fileStatus==newFile)
  {
    /* fsync directory:
     *   Whilst the file may be committed to disk, its existance in the
     *   directory may not be. Force a directory inode update with fsync()
     */
    if (fsync(dir_fd)==0)
    {
      close(dir_fd);
    }
    else
    {
      perror("fsync(2) :");
      c_io_printf(unit, "c_io_unix: open: fsync call failed on directory "
                  "for unit %" PRId64, unit);
      close(dir_fd);
      return c_io_fail;
    }
  }

  return c_io_success;
}

/*---------------------------------------------------------------------------*/

int64_t c_io_unix_close  (int64_t unit)
{
  int retVal;

  /* If we have written to file, ensure the contents are synced to disk */
  if (c_io_attributes[unit].modified_on_disk)
  {
    c_io_printf(unit,"c_io_unix: close: syncing file descriptor fd=%d",
                *(int *)c_io_attributes[unit].fileHandle);

    if (c_io_unix_sync(unit))
    {
      /* fsync failed */
      return c_io_fail;
    }
  }

  retVal=close(*(int *)c_io_attributes[unit].fileHandle);

  if (retVal==-1)
  {
    perror("close(2) :");
    c_io_printf(unit,
                "c_io_unix: close: ERROR in close(2): fd=%d ret=%d",
                *(int *)c_io_attributes[unit].fileHandle,
                retVal);
    return c_io_fail;
  }

  free(c_io_attributes[unit].fileHandle);

  return c_io_success;
}

/*---------------------------------------------------------------------------*/

int64_t c_io_unix_setpos  (int64_t unit,int64_t pos,int64_t word_len)
{
  int64_t retVal;

  if (pos==c_io_file_end)
  {
    retVal=lseek(*(int *)c_io_attributes[unit].fileHandle,
                 (off_t)(0), SEEK_END);

    c_io_attributes[unit].end_of_file=true;
  }
  else
  {
    retVal=lseek(*(int *)c_io_attributes[unit].fileHandle,
                 (off_t)((pos)*(word_len)), SEEK_SET);
    c_io_attributes[unit].end_of_file=false;
  }

  if (retVal!=-1)
  {
    retVal=c_io_success;
  }
  else
  {
    retVal=c_io_fail;
  }

  return retVal;
}

/*---------------------------------------------------------------------------*/

/*
 * The POSIX standard allows calls to write() to return "successfully" as long
 * as at least some of the data is transfered. A failure only occurs if *NO*
 * data is transfered (indicated by a return value of -1). We must therefore
 * loop until all the data is transfered, as it may take several calls to
 * write() to complete the transfer.
 */
int64_t c_io_unix_out (int64_t unit,
                       char    *array,
                       int64_t len,
                       int64_t word_len)
{
  int64_t words_written=0;
  ssize_t bytes_written=0;
  size_t  bytes_transfered=0;
  size_t  bytes_requested=(size_t)(word_len*len);

  const int unit_fd = *(int *)c_io_attributes[unit].fileHandle;

  while(
          /* we have not transfered all the data... */
          (bytes_written >= 0 && bytes_transfered < (size_t)(word_len*len))

          /* ... or there is a temporary failure, and we can retry */
          || (bytes_written == -1 && errno == EINTR)
       )
  {
    bytes_written = write(unit_fd, array, bytes_requested);

    /* bytes_written will be: -
     * +ve for a successfull (maybe partial) write
     * -ve for a failure.
     */
    if (bytes_written > 0)
    {
      bytes_transfered += (size_t)bytes_written;
      array += bytes_written;
      bytes_requested -= (size_t)bytes_written;
    }
    else
    {
      if (errno != EINTR)
      {
        perror("c_io_unix_out");
      }
      else
      {
        perror("EINTR: c_io_unix_out");
      }
    }
  }

  words_written = (int64_t)(bytes_transfered)/word_len;

  return words_written;
}

/*---------------------------------------------------------------------------*/

/*
 * The POSIX standard allows calls to read() to return "successfully" as long
 * as at least some of the data is transfered. A failure only occurs if *NO*
 * data is transfered (indicated by a return value of -1). We must therefore
 * loop until all the data is transfered, as it may take several calls to
 * read() to complete the transfer. A read length of 0 is returned when we
 * get to end-of-file (EOF).
 */
int64_t c_io_unix_in (int64_t unit,
                      char    *array,
                      int64_t len,
                      int64_t word_len)
{
  int64_t words_read=0;
  ssize_t bytes_read=0;
  size_t bytes_transfered=0;
  size_t bytes_requested=(size_t)(word_len*len);

  const int unit_fd = *(int *)c_io_attributes[unit].fileHandle;

  c_io_attributes[unit].end_of_file=false;

  while(
          /* we have not transfered all the data... */
          (bytes_read >= 0 && bytes_transfered < (size_t)len)

          /* ... or there is a temporary failure, and we can retry */
          || (bytes_read == -1 && errno == EINTR)
       )
  {
    bytes_read = read(unit_fd, array, bytes_requested);

    /* bytes_read will be: -
     * +ve for a successfull (maybe partial) read
     * -ve for a failure.
     */
    if (bytes_read >= 0)
    {
      /* Have we reached EOF? */
      if(bytes_read == 0)
      {
        /* Break from the loop, as there is nothing else there */
        c_io_printf(unit, "c_io_unix: in: %" PRId64 " bytes were requested,"
                    " but only %zu were transfered"
                    " - we have reached end-of-file : fd=%d",
                    (word_len*len), bytes_transfered, unit_fd);
        c_io_attributes[unit].end_of_file=true;
        break;
      }
      else
      {
        /* Update the counters so we can continue reading */
        bytes_transfered += (size_t)bytes_read;
        array += bytes_read;
        bytes_requested -= (size_t)bytes_read;
      }

    }
    else
    {
      if (errno != EINTR)
      {
        perror("c_io_unix_in");
      }
      else
      {
        perror("EINTR: c_io_unix_in");
      }
    }
  }

  words_read = (int64_t)(bytes_transfered)/word_len;

  return words_read;
}

/*---------------------------------------------------------------------------*/

int64_t c_io_unix_change_mode(int64_t unit,int64_t newMode)
{
  off_t returnval;

  /* avoid unused variable */
  (void)(newMode);

  returnval = lseek(*(int *)c_io_attributes[unit].fileHandle,
                    (off_t)(0),
                    SEEK_CUR);

  if(returnval!=-1)
  {
    return c_io_success;
  }
  else
  {
    return c_io_fail;
  }
}

/*---------------------------------------------------------------------------*/

int64_t c_io_unix_sync(int64_t unit)
{
  int fsync_return;

  fsync_return=fsync(*(int *)c_io_attributes[unit].fileHandle);

  if(fsync_return!=0)
  {
    perror("c_io_unix_sync: ERROR: ");
      return c_io_fail;
  }

  return c_io_success;
}

/*---------------------------------------------------------------------------*/

int64_t c_io_unix_query(int64_t unit,
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

#endif
