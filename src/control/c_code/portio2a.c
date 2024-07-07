#if defined(C95_2A)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* C language routines for portable version of UM */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

/* Standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <malloc.h>
#include <stdbool.h>
#include <unistd.h>
#include <libgen.h>
#include <fcntl.h>

/* Header files for Fortran to C interface */
#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"

/* Header files for I/O */
#include "c_portio.h"
#include "pio_umprint.h"

/* Header files for I/O timing */
#include "c_pio_timer.h"

/* byteswap headers/optimised macro */
#include "pio_byteswap.h"

/* Prototypes for internally called service functions */

static void change_endian(void *,
                          size_t,
                          size_t);

static void output_buffer(integer *,
                          um_data_t [],
                          integer *,
                          integer *,
                          integer *,
                          integer *);

static void sync_to_cache(integer *);

static inline void getpos_sized(int64_t *,
                                int64_t *,
                                const size_t * const);

static inline void setpos_sized(int64_t *,
                                int64_t *,
                                int64_t *,
                                const size_t * const);

#if defined(BUFRD_IO)
/* Prototypes for buffered I/O routines. */

static void flush_unit_buffer(integer *);

static void use_buffer(integer *,
                       void *,
                       integer *,
                       integer *,
                       integer *,
                       real *);

#endif

static void auto_portioinit(void);

/* Global variables used for I/O  */

/* portio_initialised is used to check if a portioinit() call has been made */
static bool portio_initialised=0;

static FILE *pf[MAX_UNITS]={NULL};

/* Record of the previous file operation for flushing purposes */
static enum fileop { UM_READ, UM_WRITE, UM_UNSET } prev_fileop [MAX_UNITS] ;

/* The file position in bytes */
#if defined(LFS)
  static off_t io_position[MAX_UNITS];
#else
  static long io_position[MAX_UNITS];
#endif


#if defined(BUFRD_IO)

/* Buffers to hold data before it is written */

static um_data_t *unit_buffer[MAX_UNITS]={NULL};

/* Pointer to the current address and offset (size in fortran default
 * integer/reals) within the buffer */
static integer unit_offset[MAX_UNITS];

#endif

static integer unit_intent[MAX_UNITS];

/* Define the buffer size - not in def as we want the routine
   to change it to be open */
static size_t buffer_size=512*1024;

static integer readonly = 0;

static bool open_flag[MAX_UNITS] = {false};

/* Unit properies table - one word per unit, one bit per
property at present */

/* integer file_properties;
moved to mppio.F90 as attributes array top support remote IO
*/

/* Initialise 'printstatus' control variable.  This will
be set from its Fortran equivalent. */
static integer printstatus = 0;

/* -------------------------------------------------------------------------- */

/*
 * auto_portioinit() is a safety-net to ensure code correctly initialises.
 *
 *        *** IT SHOULD NOT BE CALLED UNDER NORMAL CIRCUMSTANCES ***
 *
 * Code which leads to calls to auto_portioinit() should be fixed by adding a
 * call to portioinit() before any other portio procedures are called, and a
 * call to portioshutdown() after the last portio procedure call.
 *
 */
void auto_portioinit(void)
{
  int64_t auto_buffer_size=(int64_t)buffer_size;

  portioinit(&auto_buffer_size,NULL,NULL,NULL,NULL,NULL);

  ereport("portio2a:auto_portioinit", -1,
          "portio2a has called auto_portioinit().\n"
          "This will be because a call to a portio2a procedure was made "
          "before portioinit() has been called.\nPlease insert calls to "
          "portioinit() and portioshutdown() into your code!");
}

void set_printstatus(integer *printstatus_f /* IN: fortran equivalent */)
{
printstatus = *printstatus_f;
}


void buffin64_single
(
integer *unit,    /* Fortran unit                                          */
void *array,      /* Array into which data is read (void * to um_data_t[]) */
integer *maxlen,  /* Number of real numbers to be read                     */
integer *length,  /* Number of real numbers actually read                  */
real *status      /* Return code                                           */
)
{
struct timeval first;
struct timeval second;

int k;

/* detect if we have been called without a prior call to portioinit()         */
if (!portio_initialised) auto_portioinit();

if(open_flag[*unit])
{

/* Only need to flush if previous file operation was not a READ */

    if(prev_fileop[*unit] == UM_WRITE)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_READ;

    if (io_timer_active){
      gettimeofday( &first, NULL );
    }

    *length = (integer)(fread(array,WORD64BYTES,(size_t)(*maxlen),
                              pf[(size_t)(*unit)]));

    if (io_timer_active){
      gettimeofday( &second, NULL );
      io_update_timer(first, second, (u_integer)(*length)*WORD64BYTES, IO_B_IN);
    }

change_endian(array, WORD64BYTES, (size_t)(*maxlen));

*status=-1.0;
k=feof(pf[*unit]);

    if(k != 0)
    {
      char message[MAX_OUTSTR];

      perror("\nBUFFIN: Read Failed");
      snprintf(message, MAX_OUTSTR ,
       "BUFFIN: C I/O Error - Return code = %d", k);
      umPrint(message,"portio2a");
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      char message[MAX_OUTSTR];

      perror("\nBUFFIN: Read Failed");
      snprintf(message, MAX_OUTSTR ,
       "BUFFIN: C I/O Error - Return code = %d", k);
      umPrint(message,"portio2a");
      *status=1.0;
    }
   }
   else
        *status=3.0;

    io_position[(size_t)(*unit)]=
           io_position[(size_t)(*unit)]+(*length * WORD64BYTES);
}

void buffout64_single
(integer *unit,    /* Fortran unit                            */
void *array, /* Array from which data is written        */
integer *maxlen,   /* Number of real numbers to be written    */
integer *length,   /* Number of real numbers actually written */
real *status)      /* Return code                             */
{
  integer elsize;

  /* Default size in the size of buffer element */
  elsize = WORD64BYTES;

  #if !defined(BUFRD_IO)
    integer tmpstatus;
  #endif

  if(open_flag[*unit])
  {

/* Only need to flush if previous file operation was not a WRITE */
    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_WRITE;


    change_endian(array, WORD64BYTES, (size_t)(*maxlen));

#if defined(BUFRD_IO)
    use_buffer(unit, array, maxlen, length, &elsize, status);
#else
    output_buffer(unit, array, maxlen, length, &elsize, &tmpstatus);
    *status=(real)tmpstatus;
#endif

    change_endian(array, WORD64BYTES, (size_t)(*maxlen));
  }
  else
    *status=3.0;
}

#if defined(BUFRD_IO)
void use_buffer
(integer *unit,    /* Fortran unit                            */
void *array,       /* Array from which data is written        */
integer *maxlen,   /* Number of real numbers to be written    */
integer *length,   /* Number of real numbers actually written */
integer *elsize,   /* Size of elements to write               */
real *status)      /* Return code                             */
{
  size_t record_length;
  size_t record_offset;
  size_t i;

  um_data_t *buffer_address;

/* Set up a few values from the input */
  record_length=(size_t)(*maxlen)*(((size_t)(*elsize)<sizeof(um_data_t))
                                    +(size_t)(*elsize)/sizeof(um_data_t));
  record_offset=0;
  buffer_address=unit_buffer[*unit];

/* Buffered I/O - check if the new record will fit in the buffer */

  while((size_t)unit_offset[*unit]+record_length > buffer_size) {

/* Buffer is too small - store what we can and output that */
#if defined(NEC)
#pragma vdir nodep
#endif
    for(i=(size_t)unit_offset[*unit]; i<buffer_size; i++) {
      buffer_address[i]=((um_data_t *)array)[record_offset];
      record_length=record_length-1;
      record_offset=record_offset+1;
    }

/* Update the unit_offset */

    unit_offset[*unit]=(integer)buffer_size;
    flush_unit_buffer(unit);
  }

/* Now store the remainder of the record */

#if defined(NEC)
#pragma vdir nodep
#endif
  for(i=0; i<record_length; i++) {
    buffer_address[i+(size_t)unit_offset[*unit]]
            =((um_data_t *)array)[record_offset];
    record_offset=record_offset+1;
  }


/* Update the pointers */

  unit_offset[*unit]=unit_offset[*unit]+(integer)record_length;

  *length=*maxlen;
  *status=-1.0;

}
#endif

void output_buffer
(integer *unit,    /* Fortran unit                            */
um_data_t array[], /* Array from which data is written        */
integer *maxlen,   /* Number of real numbers to be written    */
integer *length,   /* Number of real numbers actually written */
integer *elsize,   /* Size of elements to write               */
integer *status)      /* Return code                             */
{
int k;
struct timeval first;
struct timeval second;

/* detect if we have been called without a prior call to portioinit()         */
if (!portio_initialised) auto_portioinit();

    if (io_timer_active){
      gettimeofday( &first, NULL );
    }

    *length = (integer)fwrite(array,(size_t)(*elsize),
                              (size_t)(*maxlen),pf[(size_t)(*unit)]);

    if (io_timer_active){
      gettimeofday( &second, NULL);
      io_update_timer( first, second, (u_integer)((*elsize) * (*length)),
                       IO_B_OUT);
    }

    *status=-1;
    k=feof(pf[*unit]);

    if(k != 0)
    {
      char message[MAX_OUTSTR];

      perror("\nBUFFOUT: Write Failed");
      snprintf(message, MAX_OUTSTR ,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      umPrint(message,"portio2a");
      *status=0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      char message[MAX_OUTSTR];

      perror("\nBUFFOUT: Write Failed");
      snprintf(message, MAX_OUTSTR ,
       "BUFFOUT: C I/O Error - Return code = %d", k);
      umPrint(message,"portio2a");
      *status=1;
    }

#if defined(BUFRD_IO)
    unit_offset[*unit]=0;
#endif

    io_position[*unit]=io_position[*unit]+((*length) * (*elsize));
}


void setpos_single(
int64_t *unit,        /* Fortran unit                                     */
int64_t *word_address, /* Number of default fortran int/real into file     */
int64_t *err)          /* Error checking, err = 0 no errors,err = 1 errors */
{
   const size_t word_size=sizeof(um_data_t);

  setpos_sized(unit, word_address, err, &word_size);

}

#if defined(BUFRD_IO)
void flush_unit_buffer(integer *unit)
{
  integer k;

  integer len;
  integer elsize;

  integer stat;

/* Set size to size of the elements in the buffer. */
  elsize = sizeof(um_data_t);
/* Check if there is anything to flush */

  if(unit_offset[*unit] > 0) {


    k=unit_offset[*unit];

    output_buffer(unit, unit_buffer[*unit],
                  &(unit_offset[*unit]), &len, &elsize, &stat);

/* Check for an error */

    if( (stat != -1) || (len != k) ) {

      fprintf(stderr,
       "\nFLUSH_UNIT_BUFFER: Error Flushing Buffered Data on PE %d\n",
       0);
      fprintf(stderr, "FLUSH_UNIT_BUFFER: Status is %7.1f\n", (real)stat);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length Requested was %zu\n", (size_t)k);
      fprintf(stderr,
       "FLUSH_UNIT_BUFFER: Length written   was %zu\n\n", (size_t)len);
      ereport("portio2a:flush_unit_buffer",1,"Failed in output_buffer()");
    }
    else {
      unit_offset[*unit]=0;
    }
  }
}
#endif

void open_single
(
  int64_t *unit,             /* Fortran unit                                  */
  void    *file_name,        /* File name or environ. var. (void * to char[]) */
  int64_t *char_len,         /* No of chars in file name                      */
  int64_t *intent,           /* =0 read only,!=0 read and write               */
  int64_t *environ_var_flag, /* =0 file name in environment var,              */
                             /*!=0 explicit file name                         */
  int64_t *err               /* =0 file opened,                               */
                             /* =1 file failed to open                        */
)
{
   char *fname;
   struct timeval first;
   struct timeval second;

   char *gname;
   char message[MAX_OUTSTR];

   #if defined(BUFRD_IO)
     size_t k;
     int malloc_error;
   #endif


   enum   filestat { old, new };
   enum   filestat filestatus;

   /* detect if we have been called without a prior call to portioinit()      */
   if (!portio_initialised) auto_portioinit();

   fname = calloc((size_t)(*char_len) + 1,1);

   pf[*unit] = NULL;
/* convert file name to C format */

   strncpy( fname, file_name, (size_t)(*char_len) );
   fname[ *char_len ] = '\0';
   sscanf( fname, "%s", fname );

   if (io_timer_active){
     gettimeofday( &first, NULL );
   }

   if ( *environ_var_flag == 0 )  /* File name held in environ var */
   {  gname = getenv( fname );
      if ( gname == NULL ) {
        snprintf(message, MAX_OUTSTR ,
         "OPEN:  WARNING: Environment variable %s not set",
         fname);
        umPrint(message,"portio2a");
        open_flag[*unit]=false;
        *err=1;
        free (fname);
        return;
      }
   }
   else                           /* get file name from argmt fname */
      gname = fname;


   /* Check if file exists */

   if ( access( gname, 0 ) == 0 )  {   /* file exists */

      if (printstatus >= PRINTSTATUS_NORMAL)
      {
        snprintf(message, MAX_OUTSTR ,
       "OPEN:  File %s to be Opened on Unit %d Exists",
       gname, (int) *unit );
      umPrint(message,"portio2a");
      }      /* test on printstatus */
      filestatus = old;
   }
   else  {   /* non-existent file */

      if (printstatus >= PRINTSTATUS_NORMAL)
      {
      snprintf(message, MAX_OUTSTR ,
       "OPEN:  File %s to be Opened on Unit %d does not Exist",
       gname, (int) *unit);
      umPrint(message,"portio2a");
      }      /* test on printstatus */
      filestatus = new;
   }


   if ( filestatus == old )  {

      if ( *intent == readonly )  {

         if ( ( pf[*unit] = fopen( gname, "rb" ) ) == NULL )
         {
            perror("OPEN:  File Open Failed");
            snprintf(message, MAX_OUTSTR ,
              "OPEN:  Unable to Open File %s for Reading", gname );
            umPrint(message,"portio2a");
         }
      }
      else  {   /*  *intent == read_and_write )  */

         if ( ( pf[*unit] = fopen( gname, "r+b" ) ) == NULL )
         {
            perror("OPEN:  File Open Failed");
            snprintf(message, MAX_OUTSTR ,
              "OPEN:  Unable to Open File %s for Read/Write", gname );
            umPrint(message,"portio2a");
         }
      }
   }


/* New file - check for write */
   if ( filestatus == new )  {

/* Initialise the file control word to NULL */
      pf[*unit] = NULL;

      if ( *intent == readonly )
      {
         snprintf(message, MAX_OUTSTR , "OPEN:  **WARNING: FILE NOT FOUND" );
         umPrint(message,"portio2a");
         snprintf(message, MAX_OUTSTR ,
          "OPEN:  Ignored Request to Open File %s for Reading",
           gname );
         umPrint(message,"portio2a");
      }
      else
      {
        /*  *intent == read_and_write   */

        char *dir_name = NULL;
        char *gname2 = NULL;
        int dir_fd = -1;

        /* get file handle on directory, to fsync later.
         * Calls to dirname() mangle the original string, so we must create a
         * copy first so we can retain the original filename.
         */
        gname2 = malloc(strlen(gname)+1);
        snprintf(gname2,strlen(gname)+1,"%s",gname);
        dir_name = dirname(gname2);
        dir_fd = open(dir_name, O_RDONLY);

        /* File size not given - just open the file normally */

        if (dir_fd != -1)
        {
          pf[*unit] = fopen( gname, "w+b" );

          if ( pf[*unit] == NULL )
          {
            perror("OPEN:  File Creation Failed");
            snprintf(message, MAX_OUTSTR ,
             "OPEN:  Unable to Open File %s for Read/Write", gname );
            umPrint(message,"portio2a");
          }
          else
          {
            snprintf(message, MAX_OUTSTR , "OPEN:  File %s Created on Unit %d",
                     gname, (int) *unit );
            umPrint(message,"portio2a");
          }

          /* fsync directory:
           *   Whilst the file may be committed to disk, its existance in the
           *   directory may not be. Force a directory inode update with fsync()
           */

          if (fsync(dir_fd)!=0)
          {
            perror("OPEN: directory fsync failed");
            snprintf(message, MAX_OUTSTR , "OPEN: WARNING: Failed to fsync "
                     "directory file descriptor for file %s",
                     gname);
            umPrint(message,"portio2a",.stdErrorToo=true);
            snprintf(message, MAX_OUTSTR , "OPEN: WARNING: There may be a delay"
                     " before the file appears on disk.");
            umPrint(message,"portio2a",.stdErrorToo=true);
            snprintf(message, MAX_OUTSTR , "OPEN: WARNING: In the event of "
                     "failure data may be lost. [Even after successful run "
                     "completion.]");
            umPrint(message,"portio2a",.stdErrorToo=true);
          }
          close(dir_fd);
        }
        else
        {
          perror("OPEN: File directory open failed");
          snprintf(message, MAX_OUTSTR , "OPEN: Error: Failed to open "
                   "directory file descriptor for file %s",
                   gname);
          umPrint(message,"portio2a",.stdErrorToo=true);
        }

        free(gname2);

      }
   }


   /* Set error code and open flag used by buffin and buffout */

   if( pf[*unit] == NULL )  {
      *err = 1;
      open_flag[*unit]=false;
   }
   else  {
      if ( setvbuf( pf[*unit], NULL, _IOFBF, BUFSIZ ) != 0 )  {
         perror("\n**Warning: setvbuf failed");
         *err=1;
         open_flag[*unit]=false;
       }
       else
       {
         *err = 0;
         open_flag[*unit]=true;

/*    set buffer to default size to force buffer alloc on heap */
  /*  setvbuf(pf[*unit],NULL,_IOFBF,BUFSIZ);  See above */
        }
   }
   io_position[*unit]=0;
   prev_fileop[*unit]=UM_UNSET;



/* Check if we are just reading this unit */

  if(*intent != readonly) {

#if defined(BUFRD_IO)
/* Check if this unit has a buffer already */


    if(unit_buffer[*unit] == NULL) {

/* No buffer for this file yet - create one */

/* Compute the number of bytes */
      k=buffer_size*sizeof(um_data_t);

/* Claim memory */

      unit_buffer[*unit]=malloc(k);

/* Check the error response */

      if( unit_buffer[*unit] == NULL) {
        malloc_error=errno;
        fprintf(stderr,
         "OPEN: PE %d had C I/O Error: failed in MALLOC\n",
         0);
        fprintf(stderr,
         "OPEN: Return code = %d, while claiming %zu bytes\n",
         malloc_error, k);
        ereport("portio2a:open",1,"Failure in malloc()");
      }
      else
      {
        snprintf(message, MAX_OUTSTR ,
         "OPEN:  Claimed %zu Bytes (%zu Words) for Buffering",
         k, k/sizeof(um_data_t));
        umPrint(message,"portio2a");
        snprintf(message, MAX_OUTSTR ,
         "OPEN:  Buffer Address is %p", unit_buffer[*unit]);
        umPrint(message,"portio2a");
      }

/* Set the flags for this buffer */

      unit_offset[*unit]=0;
    }
#endif

    unit_intent[*unit]=1;
  }
  else
  {
    unit_intent[*unit]=0;
  }

   if (io_timer_active){
     gettimeofday( &second, NULL );
     io_update_timer( first, second, 0, IO_F_OP );
   }

free (fname);
}

void close_single
(
int64_t *unit,             /* Fortran unit                                  */
void    *file_name,        /* File name or environ. var. (void * to char[]) */
int64_t *char_len,         /* No of chars in file name                      */
int64_t *environ_var_flag, /* =0 file name in environment var,              */
                           /*!=0 explicit file name                         */
int64_t *delete,           /* =0 do not delete file,!=0 delete file         */
int64_t *err               /* ERROR CHECKING err = 0 no errors,             */
                           /*                err = 1 Errors                 */
)
{
char *fname;
char *gname;
struct timeval first;
struct timeval second;

int i;
integer k;
char message[MAX_OUTSTR];

/* detect if we have been called without a prior call to portioinit()         */
if (!portio_initialised) auto_portioinit();

*err = 1;

fname = calloc((size_t)(*char_len) + 1,1);
/* first check to see if unit has been closed already (or not opened)*/
if(open_flag[*unit]) /* unit currently open  */
{

/* close file */
      if (io_timer_active){
        gettimeofday( &first, NULL );
      }

/*
 * Flush the file contents to disk.
 * If we are buffering, this will first include a flushing of the buffer.
 */
if(unit_intent[*unit]!=readonly)
{
  integer integer_unit, icode_tmp;

  integer_unit=(integer)*unit;

  sync_single(&integer_unit, &icode_tmp);

  if(icode_tmp!=0)
  {
    fprintf(stderr,
            "CLOSE: Cannot synchronise unit %" PRId64
            " with disk before close. Data will be lost!",*unit);
    ereport("portio2a:close", 1,
            "Failure in sync_single() call [via close_single()]");
  }

}

      k=fclose(pf[*unit]);

/* convert file name to C format */
        strncpy(fname,file_name,(size_t)(*char_len));
        fname[*char_len] = '\0';
        for (i=0; i<*char_len; i++){

            if (fname[i] == ' '){
               fname[i] = '\0';
               break;
            }
         }

        if(*environ_var_flag == 0)
          { gname = getenv( fname );
            if ( gname == NULL )
            {
              snprintf(message, MAX_OUTSTR ,
               "CLOSE: WARNING: Environment variable %s not set",
               fname);
              umPrint(message,"portio2a");
            open_flag[*unit]=false;
            free( fname );
            return;
            }
          }
        else
          gname=fname;

      if(k==0){

/* delete file */
        if(*delete != 0){
          k=remove(gname);
          if( k != 0){
            fprintf(stderr,
             "CLOSE: Cannot Delete File %s",gname);
            ereport("portio2a:open",1,"Failure in remove() (aka delete)");
          }
          else{  /*normal end to delete so:*/
            open_flag[*unit]=false; /* set unit flag to closed */
            *err = 0;
            if (printstatus >= PRINTSTATUS_NORMAL)
            {
              snprintf(message, MAX_OUTSTR ,
               "CLOSE: File %s Deleted", gname);
              umPrint(message,"portio2a");
            }
          }

        }
        else{
/* file closed */
           open_flag[*unit]=false; /* set unit flag to closed */
                                 *err = 0;
          if (printstatus >= PRINTSTATUS_NORMAL)
          {
            snprintf(message, MAX_OUTSTR ,
             "CLOSE: File %s Closed on Unit %d",
             gname, (int) *unit);
            umPrint(message,"portio2a");
          }
        }
      }
/* file not closed */
    else
    {
      snprintf(message, MAX_OUTSTR ,
           "CLOSE: Cannot Close File %s on Unit %d",
           gname, (int) *unit);
          umPrint(message,"portio2a");
    }

    if (io_timer_active){
      gettimeofday( &second, NULL );
      io_update_timer( first, second, 0, IO_F_CL );
    }

}   /* end of test for unit already closed */

else {
/* unit either closed already or not open yet */
      if (printstatus >= PRINTSTATUS_NORMAL)
      {
snprintf(message, MAX_OUTSTR ,
   "CLOSE: WARNING: Unit %d Not Opened", (int) *unit);
  umPrint(message,"portio2a");
      }      /* test on printstatus */
}

free( fname );

}



void buffout32_single
(integer *unit, /* Fortran unit                            */
 void *array, /* Array from which data is written        */
 integer *maxlen, /* Number of real numbers to be written    */
integer * length,/* Number of real numbers actually written */
 real *status) /* Return code                             */
{
  integer nsize;                          /* size of array */

  integer elsize;  /* Size of elements in array */
  integer istatus;

  /* Default size of element is size of array elements */
  elsize = WORD32BYTES; /* Size of 32 bit word */

  if (open_flag[*unit])
  {
/* Only need to flush if previous file operation was not a WRITE */
    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }
    prev_fileop[*unit] = UM_WRITE;

#if defined (FRL8)
    nsize = (*maxlen+1)/2;
#else
    nsize = *maxlen;
#endif

    change_endian(array,  WORD32BYTES, (size_t)nsize);

    /* Output array (or possibly buffer it) */
#if defined(BUFRD_IO)
    if (*maxlen % (integer)(sizeof(um_data_t)/WORD32BYTES) == 0)
    {
      /* We can happily store 32 bit words exactly into buffer.  We can use
       * nsize here since this is the number of Fortran default reals*/
      use_buffer(unit, array, &nsize, length, &elsize, status);
    }
    else
    {
      /* This is not an equal multiple of buffer size. */
      flush_unit_buffer(unit);
      output_buffer(unit, array, maxlen, length, &elsize, &istatus);
      *status=(real)(istatus);
    }
#else
    output_buffer(unit, array, maxlen, length, &elsize, &istatus);
    *status=(real)(istatus);
#endif

    change_endian(array,  WORD32BYTES, (size_t)nsize);

#if defined (BUFRD_IO)
    /* We would have used buffer if possible to fit into it */
    if (*maxlen % (integer)(sizeof(um_data_t)/WORD32BYTES) == 0)
    {
      *length = *length * (integer)(sizeof(um_data_t)/WORD32BYTES);
    }
#endif
  }
}

void buffin32_single
(
integer *unit,   /* Fortran unit                                           */
void *array,     /* Array from which data is read (void * to um_data_t[])  */
integer *maxlen, /* Number of real numbers to be read                      */
integer *length, /* Number of real numbers actually read                   */
real *status     /* Return code                                            */
)
{
  struct timeval first;
  struct timeval second;

  int nsize;

  int k;

  /* detect if we have been called without a prior call to portioinit()       */
  if (!portio_initialised) auto_portioinit();

    if (open_flag[*unit])
    {
/* Only need to flush if previous file operation was not a READ */
      if(prev_fileop[*unit] == UM_WRITE)
      {
        sync_to_cache(unit);
      }
      prev_fileop[*unit] = UM_READ;

      if (io_timer_active){
         gettimeofday( &first, NULL );
      }

      *length = (integer)fread(array,WORD32BYTES,(size_t)(*maxlen),pf [*unit]);

      if (io_timer_active){
        gettimeofday( &second, NULL );
        io_update_timer( first, second, (u_integer)(WORD32BYTES * (*maxlen)),
                         IO_B_IN );
      }


#if defined (FRL8)
      nsize = (int)(*maxlen+1)/2;
#else
      nsize = *maxlen;
#endif
      change_endian(array, WORD32BYTES, (size_t)nsize);


      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        char message[MAX_OUTSTR];

        snprintf(message, MAX_OUTSTR ,"C I/O Error: failed in BUFFIN32");
        umPrint(message,"portio2a");
        snprintf(message, MAX_OUTSTR ,"Return code = %d",k);
        umPrint(message,"portio2a");
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        char message[MAX_OUTSTR];

        snprintf(message, MAX_OUTSTR , "C I/O Error: failed in BUFFIN32");
        umPrint(message,"portio2a");
        snprintf(message, MAX_OUTSTR , "Return code = %d",k);
        umPrint(message,"portio2a");
        *status=1.0;
      }
    }
    else
      *status=3.0;
    /* Position of file must be updated */
    io_position[*unit]=io_position[*unit]+(*length*WORD32BYTES);
}

void buffin16_single(integer *unit,   /* Fortran unit                         */
                     void    *array,  /* Array from which data is read        */
                     integer *maxlen, /* Number of real numbers to be read    */
                     integer *length, /* Number of real numbers actually read */
                     real    *status) /* Return code                          */
{
  struct timeval first;
  struct timeval second;

  int k;

  /* detect if we have been called without a prior call to portioinit()       */
  if (!portio_initialised) auto_portioinit();

    if (open_flag[*unit])
    {
/* Only need to flush if previous file operation was not a READ */
      if(prev_fileop[*unit] == UM_WRITE)
      {
        sync_to_cache(unit);
      }
      prev_fileop[*unit] = UM_READ;

      if (io_timer_active){
         gettimeofday( &first, NULL );
      }

      *length = (integer)fread(array,WORD16BYTES,(size_t)(*maxlen),pf [*unit]);

      if (io_timer_active){
        gettimeofday( &second, NULL );
        io_update_timer( first, second, (u_integer)(WORD16BYTES * (*maxlen)),
                         IO_B_IN );
      }


      change_endian(array, WORD16BYTES, (size_t)(*maxlen));


      *status=-1.0;
      k=feof(pf[*unit]);
      if(k != 0)
      {
        char message[MAX_OUTSTR];

        snprintf(message, MAX_OUTSTR , "C I/O Error: failed in BUFFIN16");
        umPrint(message,"portio2a");
        snprintf(message, MAX_OUTSTR , "Return code = %d", k);
        umPrint(message,"portio2a");
        *status=0.0;
      }
      k=ferror(pf[*unit]);
      if(k != 0)
      {
        char message[MAX_OUTSTR];

        snprintf(message, MAX_OUTSTR , "C I/O Error: failed in BUFFIN16");
        umPrint(message,"portio2a");
        snprintf(message, MAX_OUTSTR , "Return code = %d", k);
        umPrint(message,"portio2a");
        *status=1.0;
      }
    }
    else
      *status=3.0;
    /* Position of file must be updated */
    io_position[*unit]=io_position[*unit]+(*length*WORD16BYTES);
}

void buffin8_single
(
integer *unit,   /* Fortran unit                                      */
void *array,     /* Array into which data is read (void * to char[] ) */
integer *maxlen, /* Number of bytes to be read                        */
integer *length, /* Number of bytes actually read                     */
real *status     /* Return code                                       */
)
{
  int k;

  /* detect if we have been called without a prior call to portioinit()       */
  if (!portio_initialised) auto_portioinit();

  if(open_flag[*unit])
  {
    *length = (integer)fread(array,1,(size_t)(*maxlen),pf[*unit]);

        *status=-1.0;
        k=feof(pf[*unit]);
    if(k != 0)
    {
      char message[MAX_OUTSTR];

      snprintf(message, MAX_OUTSTR , "C I/O Error: failed in BUFFIN8");
      umPrint(message,"portio2a");
      snprintf(message, MAX_OUTSTR , "Return code = %d",k);
      umPrint(message,"portio2a");
      *status=0.0;
    }
        k=ferror(pf[*unit]);
    if(k != 0)
    {
      char message[MAX_OUTSTR];

      snprintf(message, MAX_OUTSTR ,"C I/O Error: failed in BUFFIN8");
      umPrint(message,"portio2a");
      snprintf(message, MAX_OUTSTR ,"Return code = %d",k);
      umPrint(message,"portio2a");
      *status=1.0;
    }
   }
   else
   {
        *status=3.0;
   }

       /* Position of file must be updated */
       io_position[*unit]=io_position[*unit]+(*length*WORD8BYTES);

}

void buffout8_single(integer *unit,   /* Fortran unit                         */
                     void    *array,  /* Array from which data is written     */
                     integer *maxlen, /* Number of bytes to be written        */
                     integer *length, /* Number of bytes actually written     */
                     real    *status) /* Return code                          */
{
  integer elsize;  /* Size of elements in array */
  integer istatus;

  /* Default size of element is size of array elements */
  elsize = WORD8BYTES; /* Size of 8 bit word */

  if(open_flag[*unit])
  {

    /* Only need to flush if previous file operation was not a WRITE */
    if(prev_fileop[*unit] == UM_READ)
    {
      sync_to_cache(unit);
    }

    prev_fileop[*unit] = UM_WRITE;

    #if defined(BUFRD_IO)
    flush_unit_buffer(unit);
    #endif

    output_buffer(unit, array, maxlen, length, &elsize, &istatus);
    *status=(real)(istatus);

  }
  else
  {
    *status=3.0;
  }

}

void setpos8_single
(
int64_t *unit,         /* Fortran unit                         */
int64_t *byte_address, /* Number of bytes into file            */
int64_t *err           /* 0: no error                          */
                       /* 1: error occured                     */
)
{
  const size_t word_size=WORD8BYTES;

  setpos_sized(unit, byte_address, err, &word_size);
}

void getpos8
(int64_t *unit,        /* Fortran unit                         */
int64_t *byte_address) /* Number of bytes into file            */
{
  const size_t word_size=WORD8BYTES;

  getpos_sized(unit, byte_address, &word_size);
}

void setpos32_single(
 int64_t *unit,           /* Fortran unit                         */
 int64_t *word32_address, /* Number of 32bit ints into the file   */
 int64_t *err)            /* 0: no error                          */
                          /* 1: error occured                     */
{
  const size_t word_size=WORD32BYTES;

  setpos_sized(unit, word32_address, err, &word_size);
}

/* -------------------------------------------------------------------------- */

static inline void setpos_sized(int64_t *unit,         /* Fortran unit        */
                                int64_t *word_address, /* Number of words     */
                                                       /* into file           */
                                int64_t *err,          /* Error checking:     */
                                                       /* err = 0 no errors,  */
                                                       /* err = 1 errors      */
                                const size_t * const word_size)
{
  struct timeval first;
  struct timeval second;
  int k;

#if defined(LFS)
  off_t byte_address;
#else
  long byte_address;
#endif

  *err = 0;       /* default a successful return */

  /* detect if we have been called without a prior call to portioinit()       */
  if (!portio_initialised) auto_portioinit();

  if (io_timer_active)
  {
    gettimeofday(&first, NULL);
  }

  if (open_flag[*unit])
  {
#if defined(BUFRD_IO)
    /* check if this file is for writing */

    if(unit_intent[*unit] != readonly)
    {

      /* Does the new address match our current block address */
      if (io_position[*unit] + unit_offset[*unit] * (integer)(sizeof(um_data_t))
           == *word_address * (int64_t)(*word_size))
      {
        /* We already have a block position matching - no action required */

        if (io_timer_active)
        {
          gettimeofday(&second, NULL);
          io_update_timer(first, second, 0, IO_F_SK);
        }

        return;
      }

      /* Else we need to flush the buffer */
      integer loc_unit=(integer)*unit;
      flush_unit_buffer(&loc_unit);

    }
#endif

    /* Only do the seek if we need to */
    if (io_position[*unit] != *word_address * (int64_t)(*word_size))
    {

#if defined(LFS)
      byte_address=(off_t)(*word_address)*(off_t)(*word_size);
      k = fseeko(pf[*unit],byte_address,SEEK_SET);
#else /* LFS */
      byte_address=(long)(*word_address)*(long)(*word_size);
      k = fseek(pf[*unit],byte_address,SEEK_SET);
#endif /* LFS */

      if (k!=0)
      {
        char message[MAX_OUTSTR];

        snprintf(message, MAX_OUTSTR, "SETPOS: Seek Failed");
        umPrint(message, "portio2a");
        snprintf(message, MAX_OUTSTR ,
                 "SETPOS: Unit %d word_address = %lld",
                 (int) *unit, (long long) *word_address);
        umPrint(message, "portio2a");
        snprintf(message, MAX_OUTSTR, "SETPOS: Return code from fseek = %d", k);
        umPrint(message, "portio2a");
        *err = 1;
        ereport("portio2a:setpos", 1, "Failed in fseek[o]");
      }

      io_position[*unit]=byte_address;
    }
  }

  if (io_timer_active)
  {
    gettimeofday(&second, NULL);
    io_update_timer(first, second, 0, IO_F_SK);
  }

}

/* -------------------------------------------------------------------------- */

void setpos64_single(
              int64_t *unit,           /* Fortran unit                        */
              int64_t *word64_address, /* Number of 64bit ints into the file  */
              int64_t *err             /* 0: no error                         */
                                       /* 1: error occured                    */
             )
{
  const size_t word_size=WORD64BYTES;

  setpos_sized(unit, word64_address, err, &word_size);

}

void getpos32
(int64_t *unit,          /* Fortran unit                         */
int64_t *word32_address) /* Number of 32bit ints into the file   */
{
  const size_t word_size=WORD32BYTES;

  getpos_sized(unit, word32_address, &word_size);
}

void getpos64(
              int64_t *unit,          /* Fortran unit                         */
              int64_t *word64_address /* Number of 64bit ints into the file   */
             )
{
  const size_t word_size=WORD64BYTES;

  getpos_sized(unit, word64_address, &word_size);
}

void getpos
(int64_t *unit, /* Fortran unit                                     */
int64_t *word_address)/* Number of default fortran int/reals into file    */
{
  const size_t word_size=sizeof(um_data_t);

  getpos_sized(unit, word_address, &word_size);
}

static inline void getpos_sized(
int64_t *unit,          /* Fortran unit */
int64_t *word_address,  /* Number of default fortran int/reals into file */
const size_t * const word_size /* Fortran word size */
)
{
  if (open_flag[*unit])
  {

#if defined(LFS)
    off_t byte_address;
#else
    long byte_address;
#endif

    /* detect if we have been called without a prior call to portioinit()     */
    if (!portio_initialised) auto_portioinit();

#if defined(LFS)
    byte_address = ftello(pf[*unit]);
#else /* LFS */
    byte_address = ftell(pf[*unit]);
#endif /* LFS */
    *word_address = byte_address/(int64_t)(*word_size);

    /*
       Check that the file pointer alignment is compatible with
       the word indexing requested
    */

    if (byte_address  != *word_address * (int64_t)(*word_size))
    {
      char message[MAX_OUTSTR];

      snprintf(message, MAX_OUTSTR,
              "GETPOS: byte address (%16lX) is not aligned"
              " on a %lu byte boundary!",
              (unsigned long) byte_address, (unsigned long) (*word_size));
      umPrint(message,"portio2a");
      ereport("portio2a:getpos",1,"File pointer not on appropriate boundary");
    }

    /*
       If data is buffered then we need to return where the file pointer
       would be if the data were commited
    */

#if defined(BUFRD_IO)
    /*
       If data is buffered then we need to return where the file pointer
       would be if the data were commited
    */
    if (unit_offset[*unit]!=0)
    {
      *word_address += (unit_offset[*unit] * (integer)sizeof(um_data_t)
                                           / (integer)(*word_size) );
    }
#endif

    if(byte_address != io_position[*unit])
    {
      fprintf(stderr,
              "GETPOS: IO_POSITION is %lld, but FTELL gives %lld",
              (long long) io_position[*unit], (long long) byte_address);
      ereport("portio2a:getpos",1,"Mismatch between internal state and ftell");
    }
  }
  else
  {
    *word_address = 0;
  }

}

/******************************************************************************/
/* is_unit_open_local() determines the status (open or closed) of a unit.     */
/* Returned values are boolean "ret_code":                                    */
/*  - True  - unit is open                                                    */
/*  - False - unit is closed                                                  */
/******************************************************************************/
void is_unit_open_local(int64_t *unit, bool *ret_code)
{
  *ret_code = open_flag[*unit];
  return;
}

void
sync_to_cache
(integer * unit)
{
  if(open_flag[*unit])
  {
#if defined(BUFRD_IO)
    /* Flush our private UM buffering */
    flush_unit_buffer(unit);
#endif
    /* Now flush to disk - we could have used fseek but this is tidier */
    fflush(pf[*unit]);
  }
}


void sync_single (integer *unit, integer *icode)
{
  int fsync_return=0;
  int fileno_unit = fileno(pf[*unit]);

  *icode=0;
  sync_to_cache(unit);
  fsync_return=fsync(fileno_unit);
  if (fsync_return!=0)
  {
    /* Something has gone wrong - we failed to fsync the file */
    char message[MAX_OUTSTR];

    perror("sync_single: ERROR: ");

    snprintf(message, MAX_OUTSTR,
             "fsync() on unit %" PRId64 " (fd=%i) has failed.",
             (int64_t)*unit, fileno_unit);
    umPrint(message,"portio2a");

    *icode=fsync_return;
  }

  if (printstatus > PRINTSTATUS_NORMAL)
  {
    char message[MAX_OUTSTR];

    snprintf(message, MAX_OUTSTR ,
             "sync_single:  unit %" PRId64 " (fd=%i) has fsynced.",
             (int64_t)*unit, fileno_unit);
    umPrint(message,"portio2a");
  }
}


/****************************************************************/
/* The following function determines the size in bytes of the   */
/* file attached to a given unit.                               */
/****************************************************************/
void getextent(integer *unit,
               integer *extent,
               integer *icode)
{
  struct timeval first;
  struct timeval second;

#if defined(LFS)
  off_t c_pos;
  off_t e_pos;
#else
  long int c_pos;
  long int e_pos;
#endif

  *icode=0;

#if defined (BUFRD_IO)
  flush_unit_buffer(unit);
#endif

  if (io_timer_active)
  {
    gettimeofday(&first, NULL);
  }

#if defined (LFS)
  c_pos=ftello(pf[*unit]);
  fseeko(pf[*unit],0,SEEK_END);
  e_pos=ftello(pf[*unit]);
  fseeko(pf[*unit],c_pos,SEEK_SET);
#else
  c_pos=ftell(pf[*unit]);
  fseek(pf[*unit],0,SEEK_END);
  e_pos=ftell(pf[*unit]);
  fseek(pf[*unit],c_pos,SEEK_SET);
#endif

  if (io_timer_active)
  {
    gettimeofday(&second, NULL);
    io_update_timer(first, second, 0, IO_F_SK);
  }

  *extent=(integer)e_pos;
}

/******************************************************************************/
/* portioshutdown() will perform any finialising of portio.                   */
/* It should be called once per PE, after the final call to any portio        */
/* procedure. It is not thread-safe, so should only be called by one thread   */
/* on each PE.                                                                */
/******************************************************************************/
void portioshutdown(void)
{
  portio_initialised=0;
}

/******************************************************************************/
/* portioinit() will perform any initialising of portio.                      */
/* It allows the user to override the default size of the I/O buffer          */
/*                                                                            */
/* It is normally to be called from Fortran. It should be called once per PE, */
/* after the first call to any portio procedure. It is not thread-safe, so    */
/* should only be called by one thread on each PE.                            */
/******************************************************************************/
void portioinit
(int64_t *b_size,
 int64_t *dummyA,
 int64_t *dummyB,
 bool    *dummyC,
 int64_t *dummyD,
 bool    *io_timing)
{
    buffer_size = (size_t)*b_size;

#if defined(BUFRD_IO)
    /* For buffered I/O we will report the buffer size */
    char message[MAX_OUTSTR];

    snprintf(message, MAX_OUTSTR,
             "Buffered I/O active. Buffer size set to %" PRId64
             ", %lu byte words ",
             *b_size, (unsigned long) sizeof(um_data_t));
    umPrint(message,"portio2a");
#endif

    //Avoid unused variables
    (void)(dummyA);
    (void)(dummyB);
    (void)(dummyC);
    (void)(dummyD);
    (void)(io_timing);

    /* record that portio has been initialised */
    portio_initialised=1;
}

void portioaddhelper(integer *role)
{
  /*
   *  Stub in 2A
   */

  //Avoid Unused Arguments
  (void)(role);
}

void portiodetachallhelpers(void)
{
  /*
   *  Stub in 2A
   */
}

/* In order to be consistent with present systems all data files      */
/* will be saved in big endian format.                                */
/*                                                                    */
/* The following functions are used to convert byte order (endian)    */
/* of data on little-endian machines                                  */
/*                                                                    */

void change_endian (void *ptr_array, size_t nbytes, size_t nsize)
{

  if (get_machine_endianism()!=littleEndian) return;

  {
    int64_t error=0;

    /* Swap byte order of *ptr_array */
    if (nbytes == WORD32BYTES)
    {
      error = pio_byteswap((char *)ptr_array,
                           (int64_t)(nsize*sizeof(integer)/WORD32BYTES),
                           (int64_t)nbytes);
    }
    else
    {
      error = pio_byteswap((char *)ptr_array, (int64_t)nsize, (int64_t)nbytes);
    }

    if(error!=0)
    {
      ereport("portio2a:change_endian", error,
              "A fatal error has occured during a call to pio_byteswap.\n"
              "Data may have become corrupted.");
    }
  }
}

#else
extern int unused_unit;
#endif
