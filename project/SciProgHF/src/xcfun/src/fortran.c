#include "xcfun.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
  Functions used by the Fortran interface. Use single underscore 
  convention. Do _not_ put underscores in the names except for
  the final one.

  note about pgi compilation:
  pgcc 10.2-0 does not like FSYM(foo)FCSYM(FOO)
  instead please use FSYM(foo) FCSYM(FOO)

  If fortran is using 64 bit integers ("-i8") you must
  define INT_STAR8 when compiling this file.
 */

#ifdef FTN_UPPERCASE
#define FSYM(name)
#define FCSYM(name) name
#elif defined FTN_UPPERCASE_UNDERSCORE
#define FSYM(name)
#define FCSYM(name) name##_
#else
#define FSYM(name) name##_
#define FCSYM(name)
#endif

#ifdef INT_STAR8
typedef long long int fortranint_t;
#else
typedef int fortranint_t;
#endif

#define MAX_FORTRAN_FUNCTIONALS 7

static xc_functional fortran_functionals[MAX_FORTRAN_FUNCTIONALS] = {0};

double FSYM(xcfuve) FCSYM(XCFUVE)(void)
{
  return xcfun_version();
}

fortranint_t FSYM(xcnewf) FCSYM(XCNEWF)(void)
{
  int i;
  for (i=0;i<MAX_FORTRAN_FUNCTIONALS;i++)
    {
      if (fortran_functionals[i] == 0)
	{
	  fortran_functionals[i] = xc_new_functional();
	  return i;
	}
    }
  return -1;
}

void FSYM(xcfree) FCSYM(XCFREE)(fortranint_t *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_free_functional(fortran_functionals[*fun]);
  fortran_functionals[*fun] = 0;
}

/*
  Figure out the memory distance between points
  from the first and second pointers.
 */
void FSYM(xceval) FCSYM(XCEVAL)(fortranint_t *fun,  fortranint_t *order, 
		  fortranint_t *nr_points, 
		  double *first_density,
		  double *second_density,
		  double *first_result,
		  double *second_result)
{
  fortranint_t dpitch, rpitch;
  int ord = *order, np = *nr_points;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  dpitch = second_density - first_density;
  rpitch = second_result - first_result;
  xc_eval_vec(fortran_functionals[*fun], ord, np,
	      first_density, dpitch,
	      first_result, rpitch);
}


void FSYM(xcpotential) FCSYM(XCPOTENTIAL)(fortranint_t *fun, double *density, double *energy, double *potential)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_potential(fortran_functionals[*fun], density, energy, potential);
}

void FSYM(xcsmod) FCSYM(XCSMOD)(fortranint_t *fun, fortranint_t *mode)
{
  int m = *mode;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_set_mode(fortran_functionals[*fun], m);
}

fortranint_t FSYM(xcgett) FCSYM(XCGETT)(fortranint_t *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_get_type(fortran_functionals[*fun]);
}

fortranint_t FSYM(xcmord) FCSYM(XCMORD)(fortranint_t *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_max_order(fortran_functionals[*fun]);
}

fortranint_t FSYM(xcinle) FCSYM(XCINLE)(fortranint_t *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_input_length(fortran_functionals[*fun]);
}

fortranint_t FSYM(xcoule) FCSYM(XCOULE)(fortranint_t *fun, fortranint_t *order)
{
  int ord = *order;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_output_length(fortran_functionals[*fun],ord);
}

fortranint_t FSYM(xcdind) FCSYM(XCDIND)(fortranint_t *fun, const fortranint_t *derivative)
{
  int *d;
  int i, n = xc_input_length(fortran_functionals[*fun]);
  d = malloc(sizeof(int)*n);
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);  
  for (i=0;i<n;i++)
    d[i] = derivative[i];
  i = xc_derivative_index(fortran_functionals[*fun],d);
  free(d);
  return i;
}

static void str2ints(fortranint_t ints[], int len, const char *s)
{
  int i = 0;
  while (*s && i < len)
    ints[i++] = *s++;
  if (*s)
    ints[len-1] = 0;
  else
    ints[i] = 0;
}

#if 0
static char *ints2str(int ints[])
{
  int len = 0;
  while (ints[len])
    len++;
  char *s = new char[len+1];
  for (int i=0;i<=len;i++)
    s[i] = ints[i];
  return s;
}
#endif

void FSYM(xcspla) FCSYM(XCSPLA)(fortranint_t *text, fortranint_t *len)
{
  int l = *len;
  str2ints(text,l,xcfun_splash());
}

void FSYM(xcsnam) FCSYM(XCSNAM)(fortranint_t *dst, fortranint_t *dstlen, fortranint_t *n)
{
  const char *s;
  int nn = *n, dl = *dstlen;
  s= xc_name(nn-1);
  if (s)
    str2ints(dst,dl,s);
  else
    dst[0] = 0;
}

void FSYM(xcssho) FCSYM(XCSSHO)(fortranint_t *dst, fortranint_t *dstlen, fortranint_t *n)
{
  const char *s;
  int nn = *n, dl = *dstlen;
  s = xc_short_description(nn-1);
  str2ints(dst,dl,s);
}

void FSYM(xcslon) FCSYM(XCSLON)(fortranint_t *dst, fortranint_t *dstlen, fortranint_t *n)
{
  const char *s;   
  int nn = *n, dl = *dstlen;
  s = xc_long_description(nn-1);
  str2ints(dst,dl,s);
}

fortranint_t FSYM(xcisfu) FCSYM(XCISFU)(fortranint_t *n)
{
  int nn = *n;
  return xc_is_functional(nn-1);
}

void FSYM(xcsets) FCSYM(XCSETS)(fortranint_t *fun, fortranint_t *n, double *value)
{
  int nn = *n;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_set_param(fortran_functionals[*fun],nn-1,*value);
}

double FSYM(xcgets) FCSYM(XCGETS)(fortranint_t *fun, fortranint_t *n)
{
  int nn = *n;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_get_param(fortran_functionals[*fun],nn-1);
}
