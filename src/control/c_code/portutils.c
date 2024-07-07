/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#if !defined(_BSD_SOURCE)
#define _BSD_SOURCE
#endif

#include <time.h>
#include <memory.h>
#include <errno.h>
#include <utime.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>
#include <stdbool.h>
#include <inttypes.h>

#include "c_fort2c_prototypes.h"
#include "pio_umprint.h"
#include "c_fort2c_types.h"
#include "portutils.h"

#define FNAME_SIZE 16

/******************************************************************************/
/* Miscellaneous wrappers to c and other utility                              */
/******************************************************************************/

void date_time (integer *year,   /* year                  */
                integer *month,  /* month   Current date  */
                integer *day,    /* day      and time     */
                integer *hour,   /* hour                  */
                integer *minute, /* minute                */
                integer *second) /* second                */
{
  char s[5];
  time_t a;

  a = time(NULL);

  strftime(s,5,"%Y",localtime(&a));
  *year=atoi(s);
  strftime(s,5,"%m",localtime(&a));
  *month=atoi(s);
  strftime(s,5,"%d",localtime(&a));
  *day=atoi(s);
  strftime(s,5,"%H",localtime(&a));
  *hour=atoi(s);
  strftime(s,5,"%M",localtime(&a));
  *minute=atoi(s);
  strftime(s,5,"%S",localtime(&a));
  *second=atoi(s);
}

/******************************************************************************/

void word_length (integer *length)
{
  *length=sizeof(real);
}

/******************************************************************************/

void get_file (integer  *unit,      /* Fortran unit number         */
               char     filename[], /* File name                   */
               integer  *file_len,  /* Dimension of filename       */
               integer  *err)       /* Error checking:             */
                                   /*     err = 0 no errors,      */
                                   /*     err = 1 env var not set */
{
  char fname[FNAME_SIZE];
  char *gname;
  char message[MAX_OUTSTR];

  size_t i;
  size_t k;

  *err=0;
  /* construct environment variable name         */
  /* in form UNITnn, where nn is Fortran unit no */

  if ( *unit < 100)
  {
    snprintf(fname,FNAME_SIZE,"UNIT%02" PRIdinteger, *unit);
  }
  else
  {
    snprintf(fname,FNAME_SIZE,"UNIT%03" PRIdinteger, *unit);
  }

  /* get file name stored in environment variable UNITnn */
  gname=getenv(fname);
  if ( gname == NULL) {
    snprintf(message, MAX_OUTSTR ,
             "GET_FILE: Environment Variable %s not Set", fname);
    umPrint(message,"portutils");
    filename[0] = '\0';
    for (i=1; i<(size_t)(*file_len); i++){
      filename[i] = ' ';
    }
    *err=1;
    return;
  }
  k=strlen(gname);
  if(k >  (size_t)(*file_len)){
    snprintf(message, MAX_OUTSTR ,
             "GET_FILE: File Name too long for Allocated Storage");
    umPrint(message,"portutils");
    snprintf(message, MAX_OUTSTR ,
             "GET_FILE: Environment Variable %s", fname);
    umPrint(message,"portutils");
    snprintf(message, MAX_OUTSTR ,
             "GET_FILE: File Name %s", gname);
    umPrint(message,"portutils");
    *err = 1;
    ereport("get_file",1,"File Name too long for Allocated Storage");
  }

  /* convert file name to Fortran format */
  *err = 0;
  strncpy(filename,gname,k);
  for (i=k; i<(size_t)(*file_len); i++){
    filename[i] = ' ';
  }

}

/******************************************************************************/

/* The following function is to be called from Fortran.
 * It is a simple wrapper around usleep to allow a process to relinquish
 * control.
 */

void um_sleep (integer *sleep_time)
{
//usleep is "obsolete" according to POSIX.1-2001
//use nanosleep() instead

struct timespec sleeptime;

//sleep_time is usec; tv_nsec is nano sec [1000 * usec]
sleeptime.tv_nsec = (*sleep_time % 1000000) * 1000;
sleeptime.tv_sec = ((*sleep_time * 1000) - sleeptime.tv_nsec) / 1000000000;
nanosleep(&sleeptime , NULL);
}

/******************************************************************************/

/* 4 flavours of data copy for 32, 64 and byte and fortran default lengths */

/* i) byte */
void um_memcpy (char *dest, char *src, integer *n)
{
  memcpy((void*)dest,(void*)src,(size_t)(*n));
}

/******************************************************************************/

/* ii) fortran default */
void um_memcpy_f (char *dest, char *src, integer *n)
{
  memcpy((void*)dest,(void*)src,(size_t)(*n) * sizeof(integer));
}

/******************************************************************************/

/* iii) 64-bit */
void um_memcpy64 (void *dest, void *src, int64_t n)
{
memcpy(dest,src,(size_t)(n*WORD64BYTES));
}

/******************************************************************************/

/* iv) 32-bit */
void um_memcpy32 (void *dest,void *src, int64_t n)
{
memcpy(dest,src,(size_t)(n*WORD32BYTES));
}

/******************************************************************************/

/* Split str into separate strings by delimiter splitter */
int um_str_split (char *str, char splitter, char ***strlist)
{
  unsigned int numargs=1;
  size_t start;
  size_t end;
  size_t arg;
  size_t i;

  for (i=0;i<strlen(str);i++)
    if (str[i]==splitter)numargs++;

  *strlist = malloc ((size_t)numargs*sizeof(char *));

  start=0;
  arg=0;
   for (i=0;i<=strlen(str);i++){
    if (str[i]==splitter || i==strlen(str)){
      end=i;
      (*strlist)[arg]=malloc(end-start+1);
      memcpy((*strlist)[arg],&str[start],end-start);
      (*strlist)[arg][end-start]='\0';
      start=i+1;
      arg++;
    }
  }
  return (int)numargs;
}

/******************************************************************************/

integer file_op_cp (char *src, char *dst)
{
  FILE *f1;
  FILE *f2;
  integer ret_code=0;
  int cha;
  int count;
  char message[MAX_OUTSTR];

  snprintf(message, MAX_OUTSTR ,
           "file_op: Copy: \"%s\" to \"%s\"",
           src,dst);
  umPrint(message,"portutils");

  f1 = fopen(src,"r");
  f2 = fopen(dst,"w");
  if(f1==NULL){
    snprintf(message, MAX_OUTSTR ,
             "file_op: Copy: could not open source file: %s ",
             src);
    umPrint(message,"portutils");
    ret_code=1;
  }
  else if(f2==NULL){
    snprintf(message, MAX_OUTSTR ,
             "file_op: Copy: could not open destination file: %s ",
             dst);
    umPrint(message,"portutils");
    fclose(f1);
    ret_code=1;
  }
  else{
    count=0;
    cha=getc(f1);
    while( cha != EOF ){
      putc(cha,f2);
      count++;
      cha=getc(f1);
    }
    fclose(f1);
    fclose(f2);
    snprintf(message, MAX_OUTSTR ,
             "file_op: Copy: Complete, %d bytes",count);
    umPrint(message,"portutils");
  }
  return ret_code;
}

/******************************************************************************/

integer file_op_touch (char *touch)
{
  int ret_code=0;
  struct utimbuf ubuf;
  time_t t;

  char message[MAX_OUTSTR];

  ret_code=open(touch,O_CREAT | O_EXCL |O_NOCTTY  , S_IRUSR |S_IWUSR  );
  if (ret_code==-1){
    /* the file already exists, jsut update its timestamps */
    t=time(NULL);
    ubuf.actime = ubuf.modtime = t;
    ret_code=utime(touch,&ubuf);
    if (ret_code!=0){
      perror("file_op_touch: utime:");
      snprintf(message, MAX_OUTSTR ,
               "file_op: Updating timestamp on file: \"%s\" failed",
               touch);
      umPrint(message,"portutils");
    }else{
      snprintf(message, MAX_OUTSTR ,
               "file_op: Existing file: \"%s\" timestamp set to %ld",
               touch,t);
      umPrint(message,"portutils");
    }

  } else {
    close(ret_code);
    ret_code=0;
    snprintf(message, MAX_OUTSTR ,
             "file_op: Created: \"%s\"",
             touch);
    umPrint(message,"portutils");
  }
  return (integer)ret_code;
}

/******************************************************************************/

integer file_op_rm (char *opt, char *in)
{

  integer ret_code=0;
  char message[MAX_OUTSTR];

  bool silent=false;

  if (strcmp(opt,"-f")==0) silent=true;

  snprintf(message, MAX_OUTSTR ,
           "file_op: Delete: \"%s\"",
           in);
  umPrint(message,"portutils");

  if (unlink(in)==0){
    snprintf(message, MAX_OUTSTR ,"file_op: Delete: File deleted");
    umPrint(message,"portutils");
  }
  else{ /* report an error unless the file didn't exist and silent was set */
    if ( ! ((access(in,0)!=0) && silent)  ){
      perror("File Delete:");
      snprintf(message, MAX_OUTSTR ,"file_op: Delete: Failed");
      umPrint(message,"portutils");
      ret_code=1;
    }
  }
  return ret_code;
}

/******************************************************************************/

integer file_op_mkdir (char *in)
{
  integer ret_code=0;

  char message[MAX_OUTSTR];

  snprintf(message, MAX_OUTSTR ,
           "file_op: Create Directory: \"%s\"",
           in);
  umPrint(message,"portutils");
  if (mkdir(in,0755)==0){
    snprintf(message, MAX_OUTSTR ,
             "file_op: Create Directory: Directory Created");
    umPrint(message,"portutils");
  }
  else
  {
    /* Need to check errno before we use perror. */
    if (errno == EEXIST)
    {
      /* check that it is a directory not a file */
      struct stat s;
      if( stat(in,&s) == 0 )
      {
        if( s.st_mode & S_IFDIR )
        {
          ret_code=0;
        }
        else if( s.st_mode & S_IFREG )
        {
          /* it's a file */
          ret_code=-1;
          snprintf(message, MAX_OUTSTR ,
                   "file_op: Create Directory: A file of that name exists.");
          umPrint(message,"portutils");
        }
        else
        {
          /* it's .... ???  */
          ret_code=-1;
          snprintf(message, MAX_OUTSTR ,
                   "file_op: Create Directory: Something of that name exists.");
          umPrint(message,"portutils");
        }
      }
      else
      {
        perror("stat:");
        ret_code=-1;
        snprintf(message, MAX_OUTSTR ,
                 "file_op: Create Directory failed: stat failed on EEXIST");
        umPrint(message,"portutils");
      }
    }
    else{
      ret_code=1;
      perror("mkdir:");
      snprintf(message, MAX_OUTSTR ,"file_op: Create Directory: Failed.");
      umPrint(message,"portutils");
    }
  }
  return ret_code;
}

/******************************************************************************/

/* generic file operation and manipulation engine */

integer file_op (char *operation, integer *len)
{
/*
  operation is a colon separated command resembling a posix shell command
  that we will translate to a c99 operation
*/

  char *command;
  char **strlist;
  int num_args;
  int i;
  integer ret_code;

  char message[MAX_OUTSTR];

/* Set return code to zero */
  ret_code=0;

  /* convert file name to C format */
  command = malloc((size_t)(*len) + 1);
  strncpy( command, operation, (size_t)(*len) );
  command[ *len ] = '\0';

  /* Break the single string encoding into fields */
  num_args=um_str_split(command,':',&strlist);

  if      (strlen(strlist[0])==2 &&
           strcmp(strlist[0],"cp")==0){
    ret_code=file_op_cp(strlist[2],strlist[3]);
  }
  else if (strlen(strlist[0])==2 &&
           strcmp(strlist[0],"rm")==0){
    ret_code=file_op_rm(strlist[1],strlist[2]);
  }
  else if (strlen(strlist[0])==5 &&
           strcmp(strlist[0],"mkdir")==0){
    ret_code=file_op_mkdir(strlist[2]);
  }
  else if (strlen(strlist[0])==5 &&
           strcmp(strlist[0],"touch")==0){
    ret_code=file_op_touch(strlist[2]);
  }
  else{
    snprintf(message, MAX_OUTSTR ,
             "file_op: unrecognised operation: \"%s\"",strlist[0]);
    umPrint(message,"portutils");
    ret_code=2;
  }

  /* free up buffers */
#if defined(__INTEL_COMPILER)
  /* This loop cannot be vectorised due to the presence of free() calls.
   * Suppress any vectorisation attempt with a novector #pragma, to ensure no
   * compiler warnings are generated about failed vectorisation attempts.
   */
  #pragma novector
#endif
  for (i=0;i<num_args;i++)
  {
    free(strlist[i]);
  }
  free(strlist);
  free(command);

  return ret_code;
}

/******************************************************************************/

void shell (char *command, integer *command_len)
{
  char *fname;
  int i;
  char message[MAX_OUTSTR];

  /* convert file name to C format */

  fname = calloc( (size_t)(*command_len) + 1,1);
  strncpy(fname,command,(size_t)(*command_len));
  fname[*command_len]='\0';

  /* execute command */
  i=system(fname);
  if ( i == -1 ){
    /* command failed */
    snprintf(message, MAX_OUTSTR , "C Error: failed in SHELL");
    umPrint(message,"portutils");
  }
  free( fname );
}
