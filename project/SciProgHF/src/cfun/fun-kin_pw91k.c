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

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#define LOG log
#define EXP exp
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static int pw91k_read(const char* conf_line);
static real pw91k_energy(const FunDensProp* dp);

Functional PW91k_KinFunctional = {
  "PW91k",
  fun_true,
  pw91k_read,
  NULL,
  pw91k_energy,
  NULL,
  NULL,
  NULL
};

/* IMPLEMENTATION PART */
static int
pw91k_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
pw91k_energy(const FunDensProp* dp)
{
    real zk;
    real rhoa  = dp->rhoa;
    real rhob  = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    real A1    =   0.093907;
    real A2    =   0.26608;
    real A3    =   0.0809615;
    real A4    = 100.00;
    real A     =  76.320;
    real B1    =   0.57767e-4;

    real Cf    = (3.0/10.0)*pow((3.0/(M_PI*M_PI)),2.0/3.0);

    real kfa = pow(6*M_PI*M_PI*rhoa,1.0/3.0); 
    real kfb = pow(6*M_PI*M_PI*rhob,1.0/3.0);

    real ya  = ABS(grada)/(2*kfa*rhoa);
    real yb  = ABS(grada)/(2*kfb*rhob);

    real Fya = ( 1 + A1*ya*ASINH(A*ya) + ya*ya*(A2 - A3*EXP(-A4*ya*ya)) ) / ( 1 + A1*ya*ASINH(A*ya) + B1*ya*ya*ya*ya );
    real Fyb = ( 1 + A1*yb*ASINH(A*yb) + yb*yb*(A2 - A3*EXP(-A4*yb*yb)) ) / ( 1 + A1*yb*ASINH(A*yb) + B1*yb*yb*yb*yb );

    zk = pow(2,2.0/3.0) * Cf * ( pow(rhoa,5.0/3.0)*Fya + pow(rhob,5.0/3.0)*Fyb );

    return zk;
}

