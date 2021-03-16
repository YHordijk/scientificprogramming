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

/* Automatically generated functional code: kin_vw
   Maxima input:
    >> PI: 3.14159265358979312;
    >> ctf: (3/10)*(3*PI^2)^(2/3);
    >> 
    >> rho(rhoa, rhob) := rhoa + rhob;
    >> grad_square(grada, gradb, gradab) := grada^2 + gradb^2 + 2*gradab;
    >> 
    >> // Thomas-Fermi part 
    >> tf(rhoa,rhob):= ctf*rho(rhoa,rhob)^(5/3);
    >> 
    >> // Von Weizsaecker part 
    >> vw(rhoa,grada,rhob,gradb,gradab) := (1/8) * grad_square(grada, gradb, gradab) / rho(rhoa, rhob) ;
    >> 
    >> // Thomas-Fermi + Von Weizsaecker kinetic energy functional 
    >> K(rhoa,grada,rhob,gradb,gradab) := tf(rhoa,rhob) + (1/9) * vw(rhoa,grada,rhob,gradb,gradab) ;
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int kin_vw_read(const char* conf_line);
static real kin_vw_energy(const FunDensProp* dp);
static void kin_vw_first(FunFirstFuncDrv *ds, real factor, 
                         const FunDensProp* dp);
static void kin_vw_second(FunSecondFuncDrv *ds, real factor,
                        const FunDensProp* dp);
static void kin_vw_third(FunThirdFuncDrv *ds, real factor,
                       const FunDensProp* dp);

Functional VW_KinFunctional = {
  "kin_vw",
  fun_true,
  kin_vw_read,
  NULL,
  kin_vw_energy,
  kin_vw_first,
  kin_vw_second,
  kin_vw_third
};

/* IMPLEMENTATION PART */
static int
kin_vw_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
kin_vw_energy(const FunDensProp* dp)
{
    real t1;

    t1 = dp->rhob + dp->rhoa;
    return 2.871234000188191 * pow(t1, 5.0/3.0) + 
           (1.0/72.0)*(  dp->grada*dp->grada + dp->gradb*dp->gradb + 2.0*dp->gradab ) / t1;
   
}

static void
kin_vw_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t1, t2, t3;

    t1 = dp->rhob + dp->rhoa;
    t2 = 4.785390000313651*pow(t1, 2.0/3.0) 
         - (1.0/72.0)*(dp->grada*dp->grada + dp->gradb*dp->gradb + 2.0*dp->gradab) / (t1*t1);
    t3 = (1.0/36.0) * (1/t1);

    ds->df1000 += factor * t2;
    ds->df0100 += factor * t2;
    ds->df0010 += factor * t3 * dp->grada;
    ds->df0001 += factor * t3 * dp->gradb;
    ds->df00001 += factor * t3;
}

static void
kin_vw_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[11];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;
    real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;
    real d2fdraab, d2fdrbab;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(gradb,2.0)+2.0*gradab+pow(grada,2.0);
    t[2] = rhob+rhoa;
    t[3] = 1/pow(t[2],2.0);
    t[4] = 4.785390000313651*pow(t[2], 2.0/3.0) - (1.0/72.0)*t[1]*t[3];
    t[5] = 1/t[2];
    t[6] = (1.0/36.0)*t[5];
    t[7] = (1.0/36.0)*t[1]/pow(t[2], 3.0) + 3.190260000209101/pow(t[2], 1.0/3.0);
    t[8] = -(1.0/36.0)*t[3]*grada;
    t[9] = -(1.0/36.0)*t[3]*gradb;
    t[10] = -(1.0/36.0)*t[3];
    dfdra =t[4];
    dfdrb = t[4];
    dfdga = (1.0/36.0)*t[5]*grada;
    dfdgb = (1.0/36.0)*t[5]*gradb;
    dfdab = t[6];
    d2fdrara = t[7];
    d2fdrarb = t[7];
    d2fdraga = t[8];
    d2fdragb = t[9];
    d2fdrbrb = t[7];
    d2fdraab = t[10];
    d2fdrbab = t[10];
    d2fdgaga = t[6];
    d2fdgbgb = t[6];
    d2fdrbga = t[8];
    d2fdrbgb = t[9];
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001 += factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001 += factor*d2fdrbab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0002 += factor*d2fdgbgb;
}

static void
kin_vw_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t[16];
    real dfdra, dfdrb, dfdga, dfdgb, dfdab;
    real d2fdraga, d2fdrara, d2fdrarb, d2fdragb, d2fdrbrb;
    real d2fdrbgb, d2fdgaga, d2fdgbgb, d2fdrbga;
    real d2fdraab, d2fdrbab;
    real d3fdraraga, d3fdraragb, d3fdraraab, d3fdrbrbab;
    real d3fdrarara, d3fdrararb, d3fdragaga, d3fdrarbrb;
    real d3fdragbgb, d3fdrarbgb, d3fdrarbab, d3fdgagaga;
    real d3fdrbrbrb, d3fdrbrbga, d3fdrbrbgb, d3fdrbgbgb;
    real d3fdrbgbga, d3fdrarbga, d3fdrbgaga;
    real rhoa = dp->rhoa;
    real rhob = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;
    real gradab = dp->gradab;

    t[1] = pow(gradb,2.0)+2.0*gradab+pow(grada,2.0);
    t[2] = rhob+rhoa;
    t[3] = 1/pow(t[2],2.0);
    t[4] = 4.785390000313651*pow(t[2], 2.0/3.0) - (1.0/72.0)*t[1]*t[3];
    t[5] = 1/t[2];
    t[6] = (1.0/36.0) * t[5];
    t[7] = 1/pow(t[2],3.0);
    t[8] = (1.0/36.0)*t[1]*t[7]+3.190260000209101/pow(t[2], 1.0/3.0);
    t[9] = -(1.0/36.0)*t[3]*grada;
    t[10] = -(1.0/36.0)*t[3]*gradb;
    t[11] = -(1.0/36.0)*t[3];
    t[12] = -(1.0/12.0)*t[1]/pow(t[2],4.0)-1.0634200000697/pow(t[2], 4.0/3.0);
    t[13] = (1.0/18.0) * t[7]*grada;
    t[14] = (1.0/18.0) *t[7]*gradb;
    t[15] = (1.0/18.0) *t[7];
    dfdra = t[4];
    dfdrb = t[4];
    dfdga = (1.0/36.0)*t[5]*grada;
    dfdgb = (1.0/36.0)*t[5]*gradb;
    dfdab = t[6];
    d2fdrara = t[8];
    d2fdrarb = t[8];
    d2fdraga = t[9];
    d2fdragb = t[10];
    d2fdrbrb = t[8];
    d2fdraab = t[11];
    d2fdrbab = t[11];
    d2fdgaga = t[6];
    d2fdgbgb = t[6];
    d2fdrbga = t[9];
    d2fdrbgb = t[10];
    d3fdrararb = t[12];
    d3fdraraga = t[13];
    d3fdraragb = t[14];
    d3fdrbrbab = t[15];
    d3fdraraab = t[15];
    d3fdrarbrb = t[12];
    d3fdrarbga = t[13];
    d3fdrarbgb = t[14];
    d3fdrarbab = t[15];
    d3fdragaga = t[11];
    d3fdragbgb= t[11];
    d3fdrarara = t[12];
    d3fdrbrbrb = t[12];
    d3fdrbrbga = t[13];
    d3fdrbrbgb = t[14];
    d3fdrbgaga = t[11];
    d3fdrbgbga = t[11];
    d3fdrbgbgb = t[11];
    d3fdgagaga = 0.0;
    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001 += factor*dfdab;
    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001 += factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001 += factor*d2fdrbab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0002 += factor*d2fdgbgb;
    ds->df2010 += factor*d3fdraraga;
    ds->df2001 += factor*d3fdraragb;
    ds->df1101 += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020 += factor*d3fdragaga;
    ds->df1002 += factor*d3fdragbgb;
    ds->df3000 += factor*d3fdrarara;
    ds->df2100 += factor*d3fdrararb;
    ds->df20001 += factor*d3fdraraab;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df1200 += factor*d3fdrarbrb;
    ds->df1110 += factor*d3fdrarbga;
    ds->df0300 += factor*d3fdrbrbrb;
    ds->df0210 += factor*d3fdrbrbga;
    ds->df0201 += factor*d3fdrbrbgb;
    ds->df0120 += factor*d3fdrbgaga;
    ds->df0102 += factor*d3fdrbgbgb;
    ds->df0030 += factor*d3fdgagaga;
}
