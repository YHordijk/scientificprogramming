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

/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-revpbex.c:

   Automatically generated code implementing REVPBEX functional and
   its derivatives. It is generated by func-codegen.pl being a part of
   a "Automatic code generation framework for analytical functional
   derivative evaluation", Pawel Salek, 2005

    This functional is connected by making following changes:
    1. add "extern Functional revpbexFunctional;" to 'functionals.h'
    2. add "&revpbexFunctional," to 'functionals.c'
    3. add "fun-revpbex.c" to 'Makefile.am', 'Makefile.in' or 'Makefile'.

    This functional has been generated from following input:
    ------ cut here -------
 The revised revPBE functional is one improved version of PBE for details see
    Y. Zhang and W. Yang, Phys. Rev. Lett. 80, 890 ~1998.i                   

pi:3.14159265358979312;

xa:sqrt(grada*grada)/rhoa^(4/3);
xb:sqrt(gradb*gradb)/rhob^(4/3);

 parameters for pbex 
R:1.245;                 
d:0.066725;              
mu:d*pi^2/3;             
Sa:xa/(2*(6*pi^2)^(1/3));
Sb:xb/(2*(6*pi^2)^(1/3));

 functions for pbex 
F(S):=1+R-R/(1+mu*S^2/R);
Ea(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sa);
Eb(n):=-3/(4*pi)*(3*pi^2)^(1/3)*n^(4/3)*F(Sb);

 kernel 
K(rhoa,grada,rhob,gradb,gradab):=0.5*(Ea(2*rhoa)+Eb(2*rhob));

    ------ cut here -------
*/

 
/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stddef.h>
 
#define __CVERSION__
 
#include "functionals.h"
 
/* INTERFACE PART */
static int revpbex_isgga(void) { return 1; } /* FIXME: detect! */
static int revpbex_read(const char *conf_line);
static real revpbex_energy(const FunDensProp* dp);
static void revpbex_first(FunFirstFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void revpbex_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
 
Functional RevpbexFunctional = {
  "REVPBEX",       /* name */
  revpbex_isgga,   /* gga-corrected */
  revpbex_read,
  NULL,
  revpbex_energy,
  revpbex_first,
  revpbex_second
};
 
/* IMPLEMENTATION PART */
static int
revpbex_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}



static real
revpbex_energy(const FunDensProp *dp)
{
    real res;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1;

    t1 = pow(2.0,0.33333333333333);

   /* code */
    res = 0.5*(-1.477117532764045*t1*pow(rhob,1.333333333333333)*
        (2.245-1.245/(0.0029013741221733*pow(gradb,2.0)/pow(rhob,2.666666666666667)+
        1.0))-1.477117532764045*t1*pow(rhoa,1.333333333333333)*(2.245-
        1.245/(0.0029013741221733*pow(grada,2.0)/pow(rhoa,2.666666666666667)+
        1.0)));

    return res;
}

static void
revpbex_first_helper(real rhoa, real grada, real *res)
{    real t1, t2, t3, t4;

    t1 = pow(2.0,0.33333333333333);
    t2 = pow(grada,2.0);
    t3 = 0.0029013741221733*t2/pow(rhoa,2.666666666666667)+
        1.0;
    t4 = 1/pow(t3,2.0);

   /* code */
    res[0] = 0.5*(0.014228426342101*t1*t2*t4/pow(rhoa,2.333333333333334)-
        1.969490043685393*t1*(2.245-1.245/t3)*pow(rhoa,0.33333333333333));
    res[1] = -0.0053356598782878*t1*t4*grada/pow(rhoa,1.333333333333333);
}

static void
revpbex_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    real res[2];

    revpbex_first_helper(dp->rhoa, dp->grada, res);
   /* Final assignment */
    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        revpbex_first_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

}

static void
revpbex_second_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8;

    t1 = pow(2.0,0.33333333333333);
    t2 = pow(grada,2.0);
    t3 = 0.0029013741221733*t2/pow(rhoa,2.666666666666667)+
        1.0;
    t4 = 1/pow(t3,2.0);
    t5 = 1/pow(rhoa,2.333333333333334);
    t6 = 2.245-1.245/t3;
    t7 = 1/pow(rhoa,1.333333333333333);
    t8 = 1/pow(t3,3.0);

   /* code */
    res[0] = 0.5*(0.014228426342101*t1*t2*t4*t5-1.969490043685393*
        t1*t6*pow(rhoa,0.33333333333333));
    res[1] = -0.0053356598782878*t1*grada*t4*t7;
    res[2] = 0.5*(2.2017060260384294E-4*t1*t8*pow(grada,4.0)/
        pow(rhoa,6.0)-0.014228426342101*t1*t2*t4/pow(rhoa,3.333333333333334)-
        0.65649668122846*t1*t6/pow(rhoa,0.66666666666667));
    res[3] = 0.5*(0.014228426342101*t1*grada*t4*t5-1.6512795195288221E-4*
        t1*t8*pow(grada,3.0)/pow(rhoa,5.0));
    res[4] = 6.1922981982330839E-5*t1*t2*t8/pow(rhoa,4.0)-
        0.0053356598782878*t1*t4*t7;

}

static void
revpbex_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[5];
 
    revpbex_second_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        revpbex_second_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

}
