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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if !defined (SYS_WINDOWS)
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <math.h>
typedef int64_t fortran_int8_t;
typedef int32_t fortran_int4_t;
#ifdef INT_STAR8
  typedef fortran_int8_t fortran_int_t;
#else
  typedef fortran_int4_t fortran_int_t;
#endif

#if defined (VAR_MPI)
#include <mpi.h>
#endif

void printaddr_(void *p)
{
   fprintf(stderr,"Address is %p\n",p);
}


/*
 * Take string of tokens and return new ordering
 */

void creadreorder(char *string, int *ireord, int imoord[], int *mxreord)

{
   char *token;
   int itemp;

   *ireord = 0;
   
   token = strtok(string,",");

   while ( token != NULL )
   {
      itemp = atoi(token);
      if ( itemp <= 0 ) return;
      if ( (*ireord = *ireord + 1) > (*mxreord) ) 
      {
         *ireord = -1;
         return;
      }
      imoord[(*ireord)-1] = atoi(token);
      token = strtok(NULL,",");
   }
}

static void print_time(time_t secs, FILE *dst)
{
    if (secs > 60*60*24)
    {
	int days = secs / (60*60*24);
	fprintf(dst,"%idays",days);
	secs -= days*60*60*24;
    }
    if (secs > 60*60)
    {
	int hrs = secs / (60*60);
	fprintf(dst,"%ih",hrs);
	secs -= hrs*60*60;
    }
    if (secs > 60)
    {
	int mins = secs / (60);
	fprintf(dst,"%imin",mins);
	secs -= mins*60;
    }
    fprintf(dst,"%is",(int)secs);
}

/* Start timers, or stop and print cpu and wall time*/
void stopwatch_(fortran_int_t *start, fortran_int_t *node_nr)
{
#if !defined (SYS_WINDOWS)
    static struct timeval start_wall;
    if (*start)
    {
	gettimeofday(&start_wall, NULL);
    }
    else
    {
	struct rusage rusage;
	if (getrusage(RUSAGE_SELF,&rusage) == 0)
	{
	  /* These are the fields supported by linux, we may add others on other systems. See the man page of getrusage*/
         /* miro: deactivated, because often it's beeing printed among master data, hurting thus output

	    printf(">>>> Node %i, utime: %lld, stime: %lld, minflt: %ld, majflt: %ld, nvcsw: %ld, nivcsw: %ld\n",
		   *node_nr, (long long int)(rusage.ru_utime.tv_sec), (long long int)(rusage.ru_stime.tv_sec), 
		   rusage.ru_minflt, rusage.ru_majflt, rusage.ru_nvcsw, rusage.ru_nivcsw, rusage.ru_maxrss);
          */
	}

	if (*node_nr == 0) /* print wall time only on the master */
	{
	    struct timeval end_wall;
	    gettimeofday(&end_wall, NULL);

	    printf(">>>> Node %i, utime: %lld, stime: %lld, minflt: %ld, majflt: %ld, nvcsw: %ld, nivcsw: %ld, maxrss: %ld\n",
		   (int)*node_nr, (long long int)(rusage.ru_utime.tv_sec), (long long int)(rusage.ru_stime.tv_sec), 
		   rusage.ru_minflt, rusage.ru_majflt, rusage.ru_nvcsw, rusage.ru_nivcsw, rusage.ru_maxrss);

	    printf(">>>> Total WALL time used in DIRAC: ");
	    print_time(end_wall.tv_sec - start_wall.tv_sec,stdout);
	    printf("\n");
	}
    }
#endif
}

#define quit        quit_
    /* NOTE: if this does not link with the FORTRAN subroutine QUIT,
     * then try quit without the underscore */
extern void quit(const char* str, fortran_int_t len);

void compoff_(double *wrk, double *work, fortran_int8_t *k_offset) {
  *k_offset = wrk-work+1;
/* Aug 2008: compoff_ shall always return k_offset as a FORTRAN INTEGER*8
 * so that we don't get nonsense if address difference is bigger than 31
 * bit. /HJAaJ
 *
*/
  if ( sizeof(*k_offset) != 8 )
  {
    printf( "ERROR in COMPOFF in gpc.c: long is not I*8");
    printf( "compoff: wrk=%p work=%p diff=%li\n", wrk, work, (long)(wrk-work+1));
    printf( "k_offset=%li\n", (long)*k_offset );
    printf( "int %i long %i short int %i fortran_int8_t %i\n", (int)sizeof(int), (int)sizeof(long), (int)sizeof(short int),(int)sizeof(fortran_int8_t));
    quit( "ERROR in COMPOFF in gpc.c: long is not I*8",42);
    exit(1);
  }
/*
 *printf( "compoff: wrk=%p work=%p diff=%p\n", wrk, work, wrk-work+1);
 *printf( "k_offset=%d\n", *k_offset );
 *printf( "int %d long %d short int %d\n", sizeof(int), sizeof(long),
 *sizeof(short int));
*/
}

void distribute_env(const char* name, int mytid)
{

 /* Routine to set the value of environment variables on the slaves
    equal to that of the master if it was not defined on the slave
    Written by Luuk Visscher, october 2002                          */

#if !defined (VAR_MPI)
  return;
#else
  char *env_local, *env_master, *newname;
  int len_local, len_master, len_name;

  len_local = 0;
  if( (env_local=getenv(name)) != NULL) len_local=strlen(env_local)+1;

  if ( mytid == 0 ) len_master = len_local; 
  MPI_Bcast(&len_master, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ( mytid == 0 ) 
     env_master = env_local; 
  else
     env_master = (char *) malloc(len_master);

  MPI_Bcast(env_master, len_master, MPI_CHAR, 0, MPI_COMM_WORLD);

  if ( mytid != 0 && len_local == 0 && len_master != 0) 
  { 
    len_name = strlen(name) + len_master + 1;
    newname = (char *) malloc(len_name);

    /*
    printf("Environment variable %s not set on node %i taking value %s from master\n",
            name,mytid,env_master);
    */
    sprintf(newname,"%s=%s",name,env_master);
    putenv(newname);
    free(env_master);
  }
#endif
}

int icommsize_(const fortran_int_t *istart, const fortran_int_t *iend)
{
  return iend-istart+1;
}


/* Try to cause a segmentation fault, to enable debuggers to
   trap this and inspect the program.
*/
void crash_(void)
{
  *(int *)0 = 0;
}

/*
 Often  we need to detect "ill" real numbers - NaN, +Inf, -Inf.
 For that we use C <math.h> functions isnan(x), isinf(x).
 since not all Fortran compilers have these, IMHO(MI) functions.
*/
int is_ill_number_(double x)
{
if (isnan(x) != 0 || isinf(x) ==-1 || isinf(x) ==+1)
{
 return 1;
}
else{return 0;}
}


static double *mem_start = 0;
static double *mem_high = 0;

void set_mem_start_(double *here)
{
  mem_start = here;
  mem_high = here;
}

void push_water_mark_(double *here)
{
  if (mem_high < here)
    mem_high = here;
}

int high_water_mark_(void)
{
  return mem_high - mem_start;
}
