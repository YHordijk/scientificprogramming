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

/* Automatically generated functional code: kin_tf
   Maxima input:
    >> PI: 3.14159265358979312;
    >> ctf: (3/10)*(3*PI^2)^(2/3);
    >> 
    >> rho(rhoa,rhob):= rhoa+rhob;
    >> 
    >> // Thomas-Fermi kinetic energy functional 
    >> tf(rhoa,rhob):= ctf*rho(rhoa,rhob)^(5/3);
    >> 
    >> K(rhoa,grada,rhob,gradb,gradab):= tf(rhoa,rhob);
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int kin_tf_read(const char* conf_line);
static real kin_tf_energy(const FunDensProp* dp);
static void kin_tf_first(FunFirstFuncDrv *ds, real factor, 
                         const FunDensProp* dp);
static void kin_tf_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void kin_tf_third(FunThirdFuncDrv *ds, real factor,
                         const FunDensProp* dp);

Functional TF_KinFunctional = {
  "kin_tf",
  fun_false,
  kin_tf_read,
  NULL,
  kin_tf_energy,
  kin_tf_first,
  kin_tf_second,
  kin_tf_third
};

/* IMPLEMENTATION PART */
static int
kin_tf_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
kin_tf_energy(const FunDensProp* dp)
{
    return 2.871234000188191 * pow(dp->rhob+dp->rhoa, 5.0/3.0);
}

static void
kin_tf_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t;
    
    t = 4.785390000313651 * pow(dp->rhob+dp->rhoa, 2.0/3.0);
    
    ds->df1000 += factor * t;
    ds->df0100 += factor * t;

}

static void
kin_tf_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t1, t2, t3;

    t1 = dp->rhob + dp->rhoa;
    t2 = 4.785390000313651 * pow(t1, 2.0/3.0);
    t3 = 3.190260000209101 / pow(t1, 1.0/3.0);

    ds->df1000 += factor*t2;
    ds->df0100 += factor*t2;
    ds->df2000 += factor*t3;
    ds->df1100 += factor*t3;
    ds->df0200 += factor*t3;
}

static void
kin_tf_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t1, t2, t3, t4;
 
    t1 = dp->rhob + dp->rhoa;
    t2 = 4.785390000313651 * pow(t1, 2.0/3.0);
    t3 = 3.190260000209101 / pow(t1, 1.0/3.0);
    t4 = -1.0634200000697  / pow(t1, 4.0/3.0);

    ds->df1000 += factor * t2;
    ds->df0100 += factor * t2;
    ds->df2000 += factor * t3;
    ds->df1100 += factor * t3;
    ds->df0200 += factor * t3;
    ds->df3000 += factor * t4; 
    ds->df2100 += factor * t4; 
    ds->df1200 += factor * t4; 
    ds->df0300 += factor * t4;

}
