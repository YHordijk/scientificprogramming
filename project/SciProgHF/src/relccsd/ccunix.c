/*
 *     Copyright (c) 2019 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://www.diracprogram.org
 */

 /*
 *
 * FILE    : ccunix.c
 *
        Access to C I/O routrines from FORTAN.
        The original idea comes from a MOLPRO routine that was incorporated in a skeleton (test)
        coupled cluster code from Lee and Jayatilaka that formed the basis for RELCCSD.
        The original routines were since then (1995) stripped and revised by Luuk Visscher.

I/O routines:   openc(unit,fname,size,status)
                character*(*) fname
                integer unit,size,status
                        open file fname with unit number unit
                        integer status codes as defined below

                closec(unit)

                rdabsf(unit,a,l,p)
                wrabsf(unit,a,l,p)
                integer unit,l,p
                double precision a
                        read,write respectively l words on unit with buffer a at
                        offset p words relative to beginning of file
                        all counting done in double words

        */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#if !defined (SYS_WINDOWS)
#include <sys/times.h>
#endif
/*
#include <sys/param.h>
*/

extern char  *getenv(); 

#define MAXLENGTH       15
#define MAXUNIT         99
#define INTEGER_WORD     4       

#define PATHMAX         512 
#define SCRATCH         0
#define UNKNOWN         1
#define NEW             2
#define OLD             3

/* Check that the large file option is enabled, if not use normal lseek */

#if defined (SMALL_CFILES) || (SYS_DEC) || (SYS_LINUX) || (SYS_DARWIN)
typedef off_t off64_t;
#define lseek64 lseek
#elif defined (SYS_HPUX)
#define lseek64 lseek
#endif

typedef struct
{ int  fd ;
    off64_t size, addr;
    char fname[MAXLENGTH];
} FILE_DEFINITION;

FILE_DEFINITION files[MAXUNIT];

#define MACRO_QUIT quitc_

void MACRO_QUIT();

#define PAM_GLOBAL_SYMBOLS_UNDERSCORED

#ifdef PAM_GLOBAL_SYMBOLS_UNDERSCORED
void openc_(int *unit, char *fname, off64_t* size, int *status)
#else
void openc(unit,fcd,size,status) off64_t *size; int *unit, *status; _fcd fcd;
#endif

{ char *cp, *env ;

  char  name[PATHMAX]; 

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"openc: Unit number out of range (%d)\n",*unit);
          return;
        }

   cp=fname+ strlen( fname )-1;

 /*  if (*cp==' ') { while (*cp--==' '); *(cp+2)=NULL; }  */
  if (*cp==' ') { while (*cp--==' '); *(cp+2)='\0'; }
  if (strlen(fname)-1>MAXLENGTH) fprintf(stderr,"openc: filename too int >15\n");
       /*  Get alternate file name from environment, 
           if it exists, use it, if not, use program supplied name 
        */  
  strcpy(name,fname); 

  if((env=getenv(fname)) !=NULL) 
  strcpy(name,env); 

  switch(*status)
{ case SCRATCH: sprintf(files[*unit].fname,"Tmp%d",getpid());
                files[*unit].fd=open(name,O_RDWR|O_CREAT,0666);
                unlink(files[*unit].fname);
                break;
  case UNKNOWN: if ((files[*unit].fd=open(name,O_RDWR|O_CREAT,0666))==-1)
                { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                                }
                break;
  case NEW:     if ((files[*unit].fd=open(name,O_RDWR|O_CREAT|O_TRUNC,0666))==-1)
                        { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                        }
                break;
  case OLD:     if ((files[*unit].fd=open(name,O_RDWR))==-1)
                        { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                        }
                break;
  default:      fprintf(stderr,"openc: Unknown status\n"); 
                          MACRO_QUIT();
                          exit(1);
}
  files[*unit].size= *size;
  /* strncpy(files[*unit].fname,fname,14); *(files[*unit].fname+14)=NULL;  */
  strncpy(files[*unit].fname,fname,14); *(files[*unit].fname+14)='\0';

  *size=(lseek64(files[*unit].fd,0L,2)+511)/512;
  files[*unit].addr= -1;

}

#ifdef PAM_GLOBAL_SYMBOLS_UNDERSCORED
void wrabsf_(int* unit, char* a, int* l, off64_t* p)
#else
void wrabsf(unit,fcd,l,p) off64_t *p; int *unit, *l; _fcd fcd;
#endif

{ off64_t addr, m, n, temp;

 int pi;

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"wrabs: Unit number out of range (%d)\n",*unit);
          MACRO_QUIT();
          return;
        }
  if (!*files[*unit].fname)
        { fprintf(stderr,"wrabs: write without open file unit=%d\n",*unit);
          MACRO_QUIT();
          return;
        }
  temp = *p;
  addr= temp * INTEGER_WORD;
  temp = *l;
  m= temp * INTEGER_WORD;

  pi = (long int) p;

  if (addr!=files[*unit].addr)
        if (lseek64(files[*unit].fd,addr,0)==-1)
        { fprintf(stderr,"wrabs: Error in lseek64 of (%d:%s) -> pointer=%XH\n",*unit,files[*unit].fname,pi);
          files[*unit].addr= -1;
          MACRO_QUIT();
          return;
        }

  if ((n=write (files[*unit].fd,a,m))!=m)
      { fprintf(stderr,"wrabs: Error in writing %d words in file %s with unit %d\n",*l,files[*unit].fname,*unit);
        MACRO_QUIT();
      }
  files[*unit].addr=addr+n;
}

#ifdef PAM_GLOBAL_SYMBOLS_UNDERSCORED
void rdabsf_(int* unit, char* a, int* l, off64_t* p)
#else
void rdabsf(unit,fcd,l,p) off64_t *p; int *unit, *l; _fcd fcd;
#endif

{ off64_t addr, m, n, temp;
  int pi;

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"rdabs: Unit number out of range (%d)\n",*unit);
          MACRO_QUIT();
          return;
        }
  if (!*files[*unit].fname)
        { fprintf(stderr,"rdabs: read without open file unit=%d\n",*unit);
          MACRO_QUIT();
          return;
        }
  temp = *p;
  addr= temp * INTEGER_WORD;
  temp = *l;
  m= temp * INTEGER_WORD;
 
  pi = (long int)p; /* Convert for displaying the address p is pointing to */

  if (addr!=files[*unit].addr)
        if (lseek64(files[*unit].fd,addr,0)==-1)
        { fprintf(stderr,"rdabs: Error in lseek64 of (%d:%s) -> pointer=%XH\n",*unit,files[*unit].fname,pi);
          files[*unit].addr= -1;
          MACRO_QUIT();
          return;
        }

  if ((n=read (files[*unit].fd,a,m))!=m)
      { fprintf(stderr,"rdabs: Error in reading %d words in file %s with unit %d\n",*l,files[*unit].fname,*unit);
        MACRO_QUIT();
      }
  files[*unit].addr=addr+n;
}

#ifdef PAM_GLOBAL_SYMBOLS_UNDERSCORED
void closec_(int* unit)
#else
void closec(unit) int *unit;
#endif

{ if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"closec: unit out of range\n");
           MACRO_QUIT();
           return;
        }
  if (files[*unit].fd)
        { close(files[*unit].fd);
          files[*unit].fd=files[*unit].addr=files[*unit].size=0;
          /* *files[*unit].fname=NULL; */
          *files[*unit].fname='\0';
        }
  return;
}




#if defined (SYS_AIX)
void
unlink_(char* fname)
{
  unlink(fname); /* ... and pray it does not crash.                  */
                 /* alternatively, check out ts_fmm version of gpc.c */
}
#endif /* SYS_AIX */
