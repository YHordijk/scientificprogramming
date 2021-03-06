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

/* fun-lb94.c:
   implementation of Exchange-correlation potential with correct 
   asymptotic behavior by R. van Leeuwen and E. J. Baerends:

   [ van Leeuwen and EJ Baerends, Phys Rev A 49, 2421 (1994)]
   See also comments in Gisbergen et al, JCP 105(8) 3142.

   (c) P. Salek, oct 2003 - the working implementation.
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  lb94_isgga(void) { return 1; }
static int  lb94_read(const char* conf_line);
static real lb94_energy(const FunDensProp* dens_prop);
static void lb94_first(FunFirstFuncDrv *ds, real factor, 
                       const FunDensProp* dens_prop);
static void lb94_second(FunSecondFuncDrv *ds, real factor,
                        const FunDensProp* dens_prop);

static void lb94_third(FunThirdFuncDrv *ds, real factor,
                       const FunDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94_fourth(FunFourthFuncDrv *ds, real factor,
                        const const FunDensProp* dens_prop);
#endif

Functional LB94Functional = {"LB94",      /* name */
                             lb94_isgga,  /* gga-corrected */
                             lb94_read,   /* set common blocks */
                             NULL,         /* reporter */
                             lb94_energy, 
                             lb94_first,
                             lb94_second,
                             lb94_third
#ifdef FOURTH_ORDER_DERIVATIVES
                             ,lb94_fourth,
#endif
};

/* IMPLEMENTATION PART */

static int
lb94_read(const char* conf_line)
{
    fun_set_hf_weight(0.0);
    return 1;
}

/* lb94_energy:
   lb94 threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real LB94_THRESHOLD = 1e-14;
static const real BETA = 0.05;

static real
lb94_energy(const FunDensProp* dp)
{
    return SlaterFunctional.func(dp)+VWNFunctional.func(dp);
}

static void
lb94_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real vx;

  /*real rho    = dp->rhoa + dp->rhob;*/
    real rho    = dp->rhoa;
    real rho13 = pow(rho, 1.0/3.0);
  /*real grad = dp->grada + dp->gradb;*/
    real grad = dp->grada;
    real rho43=rho*rho13;
    real scaled_grad, sg2;
    scaled_grad = grad/(rho43>1e-13 ? rho43 : 1e-13);
    sg2   = scaled_grad*scaled_grad;


    vx = -BETA*rho13*sg2/
        (1+3*BETA*scaled_grad*asinh(scaled_grad));

    ds->df1000 += vx*factor;
    ds->df0100 += vx*factor;

    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional.first(ds,   factor, dp);
}

void
lb94pot_(real *ds, const real* weight, const real* rho, const real* grad)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    /* HP does not grok C99's { *rhoa, *rhob, *grada, *gradb }; */
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    drv1_clear(&drvs);
    lb94_first(&drvs, *weight, &dp);
    /* we should really use a separate structure here
     * but this time, we cheat */
    ds[0] = drvs.df1000;
    ds[1] = drvs.df0010 + 0.5*drvs.df00001* *grad;
}

static void
lb94_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
/* according to the authors, it is equivalent to ALDA for LR and higher.
 * See comments in Gisbergen et al, JCP 105(8) 3142. */
#if 0
    real
        x1 = pow(dp->gradb,2.0)+2.0*dp->gradab+pow(dp->grada,2.0),
        x2 = dp->rhob+dp->rhoa,
        x3 = 1/pow(x2,2.333333333333334),
        x4 = sqrt(x1),
        x5 = 1/pow(x2,1.333333333333333),
        x6 = asinh(x4*x5),
        x7 = 0.15*x4*x5*x6+1.0,
        x8 = 1/x7,
        x9 = -0.05*x1*x3*x8,
        x10 = 1/pow(x2,2.666666666666667),
        x11 = 1/sqrt(x1*x10+1.0),
        x12 = 1/pow(x7,2.0),
        x13 = 0.11666666666667*x1*x8/pow(x2,3.333333333333334)
        +0.05*x1*x12*x3*(-0.2*x4*x3*x6-0.2*x1*x11/pow(x2,3.666666666666667)),
        x14 = 1/x4,
        x15 = 0.05*x1*x3*(0.15*dp->grada*x14*x5*x6+0.15*dp->grada*x10*x11)*x12
        -0.1*dp->grada*x3*x8,
        x16 = 0.05*x1*x3*(0.15*dp->gradb*x14*x5*x6+0.15*dp->gradb*x10*x11)*x12
        -0.1*dp->gradb*x3*x8;
    
    ds->df1000 += x9*factor;
    ds->df0100 += x9*factor;

    ds->df2000 += x13*factor;
    ds->df0200 += x13*factor;
    ds->df1100 += x13*factor;

    ds->df1010 += x15*factor;
    ds->df1001 += x16*factor;
    ds->df10001 += (0.05*x1*x3*(0.15*x14*x5*x6+0.15*x10*x11)*x12-0.1*x3*x8)
        *factor;
    ds->df0110 += x15*factor;
    ds->df0101 += x16*factor;
#endif    
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional.second(ds,   factor, dp);
}

 
static void
lb94_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional.third(ds,   factor, dp);
}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lb94_fourth(FourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    SlaterFunctional.fourth(ds, factor, dp);
    VWNFunctional.fourth(ds,   factor, dp);
}
#endif
