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
 * GLLBhole
 *
 * The hole part of GLLB-potential for the SAOP asymptotic correction
 * this should not be used without the response potential in SAOP and .ALDA
 * (both are activated by the keyword .SAOP!)
 *
 * O. Gritsenko, R. van Leeuwen, E. van Lenthe, E. J. Baerends, Phys. Rev. A 51, 1944 (1995).
 *
 * also have a look at equations (2.2) - (2.6) in J. Chem. Phys. 112, 1344 (2000).
 *
 * implementation: Radovan Bast - last modification 21/04/2005
 *
 */

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define __CVERSION__

#include "functionals.h"

/*
 *
 * interface
 *
 */

static int  gllbhole_isgga(void) { return 1; }
static int  gllbhole_read(const char* conf_line);
static real gllbhole_energy(const FunDensProp* dens_prop);
static void gllbhole_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dens_prop);
static void gllbhole_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);

static void gllbhole_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp);

Functional GLLBholeFunctional = {
  "GLLBhole",       /* name */
  gllbhole_isgga,   /* gga-corrected */
  gllbhole_read,    /* set common blocks */
  NULL,             /* reporter */
  gllbhole_energy,
  gllbhole_first,
  gllbhole_second,
  gllbhole_third
};

/*
 *
 * implementation
 *
 */

static int
gllbhole_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real GLLBhole_THRESHOLD = 1e-14;

static real
gllbhole_energy(const FunDensProp* dp)
{
  return SlaterFunctional.func(dp)+BeckeFunctional.func(dp)+PW91cFunctional.func(dp);
}

static void
gllbhole_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real vx;

  vx = 2.0*(SlaterFunctional.func(dp)+BeckeFunctional.func(dp)+PW91cFunctional.func(dp));

  ds->df1000 += vx*factor;
  ds->df0100 += vx*factor;
}

void
gllbholepot_(real *ds, const real* weight, const real* rho, const real* grad)
{
  FunFirstFuncDrv drvs;
  FunDensProp dp;

  dp.rhoa = dp.rhob = *rho *0.5;
  dp.grada = dp.gradb = *grad *0.5;
  dp.gradab = dp.grada * dp.gradb;
  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;
  drv1_clear(&drvs);

  gllbhole_first(&drvs, *weight, &dp);

  ds[0] = drvs.df1000/ *rho;
  ds[1] = 0.0;
}

static void
gllbhole_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  fun_printf("gllbhole_second not implemented"); exit(1);
}

static void
gllbhole_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  fun_printf("gllbhole_third not implemented"); exit(1);
}
