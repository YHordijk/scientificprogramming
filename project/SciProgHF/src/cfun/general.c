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
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* Use BSD's strncasecmp(); if there is a platform that has no strncasecmp()
 * ask pawsa@theochem.kth.se for replacement */
#define _DEFAULT_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

#include "general.h"
#include "functionals.h"
#ifdef PRG_DIRAC
#include "inforb.h"
#include "dcbham.h"
#endif

#ifdef VAR_MPI
#include <mpi.h>
#include "infpar.h"
#endif

/* C-wide constants */
const int  ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
const real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;

/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
static char* EmbXcConfString = NULL;
static char* EmbKinConfString = NULL;

#ifdef PRG_DIRAC
extern void dftsetcam_(real *w, real *b);

static void
dft_set_hf_weight(real hfxfac)
{
  dcrham_.hfxfac = hfxfac;
}

static real
dft_get_hf_weight(void)
{
    return dcrham_.hfxfac;
}

static void
dft_set_cam_param(real w, real be) { dftsetcam_(&w, &be);}
#endif

static int
dft_func_input (const char* line, int* inperr, int len, char** ConfString,
                Functional* avail[], Functional** func);

int
FSYM(dftinput)(const char* line, int * inperr, int len)
{
    fun_printf        = fort_print;
#ifdef PRG_DIRAC
    fun_set_hf_weight = dft_set_hf_weight;
    fun_get_hf_weight = dft_get_hf_weight;
    fun_set_cam_param = dft_set_cam_param;
#endif
    return dft_func_input (line, inperr, len, &DftConfString,
                           available_functionals, &selected_func); 
}

int
FSYM(dft_emb_set_xcfunc)(const char* line, int * inperr, int len)
{
    return dft_func_input (line, inperr, len, &EmbXcConfString,
                           available_functionals, &emb_xc_func); 
}

int
FSYM(dft_emb_set_kinfunc)(const char* line, int * inperr, int len)
{
    return dft_func_input (line, inperr, len, &EmbKinConfString,
                           available_kin_functionals, &emb_kin_func); 
}

static int
dft_func_input (const char* line, int* inperr, int len, char** ConfString,
                Functional* avail[], Functional** func) 
{
    char func_name[20];
    int i, off;

    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(*ConfString) free(*ConfString);
    i++;
    for(off=0; line[off] && isspace((int)line[off]); off++)
        ;
    *ConfString = malloc(i+1-off);
    strncpy(*ConfString, line+off, i-off); 
    (*ConfString)[i-off] = '\0';
    sscanf(*ConfString,"%20s", func_name);
    for(i=0; avail[i]; i++)
        if(strcasecmp(avail[i]->name, func_name)==0) {
            int ok;
            *func = avail[i];
            /* running read function of functional if available */
            ok = (*func)->read ?
                 (*func)->read((*ConfString)+strlen(func_name)) : 1; 
	    if(!ok) { /* read functional failed */
               fort_print("Read function for '%s' failed. Aborting.\n", func_name);
	       (*inperr)++;
            }
            return ok;
        }

    fort_print("Unknown functional '%s'. Aborting.\n", func_name);
    funlistfuncs(avail);
    (*inperr)++;

    return 0; /* failed to find */
}

static int do_xalda   = 0;
static int do_alda_hs = 0; /* ALDA for      density part (     hermitian part), off by default */
static int do_alda_ha = 0; /* ALDA for spin density part (anti-hermitian part), off by default */

void
setxalda_(const int *xalda)
{
  do_xalda = *xalda;
}

void
setaldahs_(const int *aldahs)
{
  do_alda_hs = *aldahs;
}

void
setaldaha_(const int *aldaha)
{
  do_alda_ha = *aldaha;
}

void
dftpot0_(FirstDrv *ds,
         const real* weight,
         const real* rho,
         const real* grad)
{
  FunFirstFuncDrv drvs;
  FunDensProp dp;

  real grada;
  real bygrada;
  real by2, by4;

  drv1_clear(&drvs);

  dp.rhoa   = dp.rhob  = *rho *0.5;
  dp.grada  = dp.gradb = *grad*0.5;
  dp.gradab = dp.grada * dp.gradb;

  grada  = *grad*0.5;

  bygrada  = 1.0/grada;

  by2   = 0.50000000;
  by4   = 0.25000000;

  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;

  selected_func->first(&drvs, *weight, &dp);

  ds->fR   = drvs.df1000;
  ds->fZ   = drvs.df0010*by4*bygrada +drvs.df00001*by4;
}

void
dftpot1_(SecondDrv *ds,
         const real* weight,
         const real* rho,
         const real* grad)
{
  FunSecondFuncDrv drvs;
  FunDensProp dp;

  real grada, grada2, grada3;
  real bygrada, bygrada2, bygrada3;
  real by2, by4, by8, by16, by32;
  real x;

  drv2_clear(&drvs);

  dp.rhoa   = dp.rhob  = *rho *0.5;
  dp.grada  = dp.gradb = *grad*0.5;
  dp.gradab = dp.grada * dp.gradb;

  grada  = *grad*0.5;
  grada2 = grada*grada;
  grada3 = grada*grada2;

  bygrada  = 1.0/grada;
  bygrada2 = 1.0/grada2;
  bygrada3 = 1.0/grada3;

  by2   = 0.50000000;
  by4   = 0.25000000;
  by8   = 0.12500000;
  by16  = 0.06250000;
  by32  = 0.03125000;

  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;

  if (do_alda_hs) {
    if (do_xalda)
      x = 1.0 - fun_get_hf_weight();
    else
      x = 1.0;
    SlaterFunctional.second(&drvs, *weight*x, &dp);
    VWNFunctional.second(&drvs, *weight, &dp); }
  else
    selected_func->second(&drvs, *weight, &dp);

  ds->fR   = drvs.df1000;
  ds->fRR  = drvs.df2000*by2 +drvs.df1100*by2;

  if (!do_alda_hs) {
     ds->fZ   = drvs.df0010*by4*bygrada +drvs.df00001*by4;
     ds->fRZ  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada +drvs.df10001*by4;
     ds->fZZ  = drvs.df00101*by8*bygrada +drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 +drvs.df00002*by16;
  }
}

void
dftpot2_(ThirdDrv *ds,
         const real* weight,
         const real* rho, 
         const real* grad)
{
/***********************************************************************/
/*                                                                     */
/*    Function that will return the functional derivatives needed for  */
/*    quadratic response calculations.                                 */
/*                                                                     */
/*      E_{xc} = \int e_{xc} d\tau                                     */
/*                                                                     */
/*    This routine will return the derivatives of e_{xc}  w.r.t.       */
/*    the electron density \rho and the entity                         */
/*      \zeta = \nabla\rho \cdot \nabla\rho,                           */
/*    which is the quadratic response equations are expressed in.      */
/*                                                                     */
/*    Below:                                                           */
/*    R corresponds to derivative w.r.t. \rho                          */
/*    Z corresponds to derivative w.r.t. \zeta                         */
/*                                                                     */
/*    Input:                                                           */
/*    w       --- weight                                               */
/*    rho     --- \rho in current grid point                           */
/*    grad    --- |\nabla\rho| in current grid point                   */
/*    triplet --- triplet flag                                         */
/*                                                                     */
/*    Output:                                                          */
/*    *ds --- contains the required derivatives up to third order.     */
/*                                                                     */
/*     Written by johhe Jan 2007                                       */
/*                                                                     */
/***********************************************************************/

    /* Calculated different powers of the norm of the gradient needed */
    real grad2 = *grad* *grad;
    real grad3 = grad2* *grad;
    real grad4 = grad3* *grad;
    real grad5 = grad4* *grad;
  real x;

    FunThirdFuncDrv drvs;
    FunDensProp dp;
    drv3_clear(&drvs);

    /* HP does not grok C99's { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5; 
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;

    /* Get the functional derivatives for the appropriate functional */
    if (do_alda_hs) {
      if (do_xalda)
        x = 1.0 - fun_get_hf_weight();
      else
        x = 1.0;
      SlaterFunctional.third(&drvs, *weight*x, &dp);
      VWNFunctional.third(&drvs, *weight, &dp); }
    else
      selected_func->third(&drvs, *weight, &dp);
    
    /* Evaluate the different derivatives */
    ds->fR   = drvs.df1000;
    ds->fZ   = 0.5* (drvs.df0010/ *grad + 0.5*drvs.df00001);
    ds->fRR  = 0.5* (drvs.df2000 + drvs.df1100);
    ds->fRZ  = 0.25* ((drvs.df1010 + drvs.df1001)/ *grad + drvs.df10001);
    ds->fZZ  = 0.125* (drvs.df0020 + drvs.df0011)/grad2 - 0.25*drvs.df0010/grad3;
    ds->fRRR = 0.25*(drvs.df3000 + 3*drvs.df2100);
    ds->fRRZ = 0.125* ((drvs.df2010 + drvs.df2001 + 2.0*drvs.df1110)/ *grad +
                       drvs.df20001 + drvs.df11001);
    ds->fRZZ = 0.0625* (drvs.df1020 + drvs.df1002 + 2*drvs.df1011)/ grad2 -
               0.125* (drvs.df1010 + drvs.df1001)/ grad3;
    ds->fZZZ = 0.03125* (drvs.df0030 + 3* drvs.df0021)/ grad3 -
               0.0625*3* (drvs.df0020 + drvs.df0011)/ grad4 +
               0.125*3* drvs.df0010/ grad5;
}

void 
sdftpot1_(SecondDrv *ds,
          const real* weight,
          const real* rho, 
          const real* grad)
{
  FunSecondFuncDrv drvs;
  FunDensProp dp;

  real grada, grada2, grada3;
  real bygrada, bygrada2, bygrada3;
  real by2, by4, by8, by16, by32;
  real x;

  drv2_clear(&drvs);

  dp.rhoa   = dp.rhob  = *rho *0.5;
  dp.grada  = dp.gradb = *grad*0.5;
  dp.gradab = dp.grada * dp.gradb;

  grada  = *grad*0.5;
  grada2 = grada*grada;
  grada3 = grada*grada2;

  bygrada  = 1.0/grada;
  bygrada2 = 1.0/grada2;
  bygrada3 = 1.0/grada3;

  by2   = 0.50000000;
  by4   = 0.25000000;
  by8   = 0.12500000;
  by16  = 0.06250000;
  by32  = 0.03125000;

  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;

  if (do_alda_ha) {
    if (do_xalda)
      x = 1.0 - fun_get_hf_weight();
    else
      x = 1.0;
    SlaterFunctional.second(&drvs, *weight*x, &dp);
    VWNFunctional.second(&drvs, *weight, &dp); }
  else
    selected_func->second(&drvs, *weight, &dp);
       
//ds->fR   = drvs.df1000;
//ds->fRR  = drvs.df2000*by2 +drvs.df1100*by2;
  ds->fSS  = drvs.df2000*by2 -drvs.df1100*by2;

  if (!do_alda_ha) {
//   ds->fZ   = drvs.df0010*by4*bygrada +drvs.df00001*by4;
//   ds->fZZ  = drvs.df00101*by8*bygrada +drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 +drvs.df00002*by16;
//   ds->fRZ  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada +drvs.df10001*by4;
     ds->fX   = drvs.df0010*by4*bygrada -drvs.df00001*by4;
     ds->fYY  = drvs.df0020*by8*bygrada2 -drvs.df0011*by8*bygrada2 -drvs.df0010*by8*bygrada3;
     ds->fSY  = drvs.df1010*by4*bygrada -drvs.df1001*by4*bygrada;
  }

/*ds->fX   = drvs.df0010*by4*bygrada -drvs.df00001*by4;

  ds->fSS  = drvs.df2000*by2 -drvs.df1100*by2;
  ds->fYY  = drvs.df0020*by8*bygrada2 -drvs.df0011*by8*bygrada2 -drvs.df0010*by8*bygrada3;
  ds->fXX  = drvs.df0020*by32*bygrada2 -drvs.df00101*by8*bygrada +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 +drvs.df00002*by16;
  ds->fRX  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada -drvs.df10001*by4;
  ds->fSY  = drvs.df1010*by4*bygrada -drvs.df1001*by4*bygrada;
  ds->fZX  = drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 -drvs.df00002*by16;*/
}

void 
sdftpot2_(ThirdDrv *ds,
          const real* weight,
          const real* rho, 
          const real* grad)
{
  FunThirdFuncDrv drvs;
  FunDensProp dp;

  real grada, grada2, grada3, grada4, grada5;
  real bygrada, bygrada2, bygrada3, bygrada4, bygrada5;
  real by2, by4, by8, by16, by32, by64, by128, by256;
  real x;

  drv3_clear(&drvs);

//#define GET_REF_DATA
#ifdef GET_REF_DATA
//  this is useful for obtaining reference data
//  of course you have to make sure that the *.c
//  implementation is correct

    dp.rhoa   = 0.39E+02;
    dp.rhob   = 0.38E+02;
    dp.grada  = sqrt(0.81E+06);
    dp.gradb  = sqrt(0.82E+06);
    dp.gradab = 0.82E+06;

//  the derivatives below can be obtained with ->second under dftpot1
    selected_func->third(&drvs, 1.0, &dp);

    fun_printf("drvs.df00000 =%20.12e", selected_func->func(&dp));
    fun_printf("drvs.df10000 =%20.12e", drvs.df1000);
    fun_printf("drvs.df01000 =%20.12e", drvs.df0100);
    fun_printf("drvs.df00100 =%20.12e", drvs.df0010*0.5/dp.grada);
    fun_printf("drvs.df00010 =%20.12e", drvs.df0001*0.5/dp.gradb);
    fun_printf("drvs.df00001 =%20.12e", drvs.df00001);
    fun_printf("drvs.df20000 =%20.12e", drvs.df2000);
    fun_printf("drvs.df11000 =%20.12e", drvs.df1100);
    fun_printf("drvs.df10100 =%20.12e", drvs.df1010*0.5/dp.grada);
    fun_printf("drvs.df10010 =%20.12e", drvs.df1001*0.5/dp.gradb);
    fun_printf("drvs.df10001 =%20.12e", drvs.df10001);
    fun_printf("drvs.df02000 =%20.12e", drvs.df0200);
    fun_printf("drvs.df01100 =%20.12e", drvs.df0110*0.5/dp.grada);
    fun_printf("drvs.df01010 =%20.12e", drvs.df0101*0.5/dp.gradb);
    fun_printf("drvs.df01001 =%20.12e", drvs.df01001);
#endif

  dp.rhoa   = dp.rhob  = *rho *0.5;
  dp.grada  = dp.gradb = *grad*0.5;
  dp.gradab = dp.grada * dp.gradb;

  grada  = *grad*0.5;
  grada2 = grada*grada;
  grada3 = grada*grada2;
  grada4 = grada*grada3;
  grada5 = grada*grada4;

  bygrada  = 1.0/grada;
  bygrada2 = 1.0/grada2;
  bygrada3 = 1.0/grada3;
  bygrada4 = 1.0/grada4;
  bygrada5 = 1.0/grada5;

  by2   = 0.50000000;
  by4   = 0.25000000;
  by8   = 0.12500000;
  by16  = 0.06250000;
  by32  = 0.03125000;
  by64  = 0.01562500;
  by128 = 0.00781250;
  by256 = 0.00390625;

  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;

  if (do_alda_ha) {
    if (do_xalda)
      x = 1.0 - fun_get_hf_weight();
    else
      x = 1.0;
    SlaterFunctional.third(&drvs, *weight*x, &dp);
    VWNFunctional.third(&drvs, *weight, &dp); }
  else
    selected_func->third(&drvs, *weight, &dp);

  ds->fR   = drvs.df1000;
  ds->fZ   = drvs.df0010*by4*bygrada +drvs.df00001*by4;

  ds->fRR  = drvs.df2000*by2 +drvs.df1100*by2;
  ds->fZZ  =                           drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3;
  ds->fRZ  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada +drvs.df10001*by4;

  ds->fRRR = drvs.df3000*by4 +3.0*drvs.df2100*by4;
  ds->fZZZ = drvs.df0030*by256*bygrada3 +3.0*drvs.df0021*by256*bygrada3                                     
            -3.0*drvs.df0020*by256*bygrada4 -3.0*drvs.df0011*by256*bygrada4 +3.0*drvs.df0010*by256*bygrada5;
  ds->fRRZ = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada +drvs.df1110*by8*bygrada +drvs.df20001*by8 +drvs.df11001*by8;
  ds->fRZZ =                                                       drvs.df1020*by64*bygrada2 +drvs.df1011*by32*bygrada2
            +drvs.df1002*by64*bygrada2 -drvs.df1010*by64*bygrada3 -drvs.df1001*by64*bygrada3;
/*ds->fR   = drvs.df1000;
  ds->fZ   = drvs.df0010*by4*bygrada +drvs.df00001*by4;

  ds->fRR  = drvs.df2000*by2 +drvs.df1100*by2;
  ds->fZZ  = drvs.df00101*by8*bygrada +drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 +drvs.df00002*by16;
  ds->fRZ  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada +drvs.df10001*by4;

  ds->fRRR = drvs.df3000*by4 +3.0*drvs.df2100*by4;
  ds->fZZZ = 3.0*drvs.df00102*by64*bygrada +3.0*drvs.df00201*by128*bygrada2 +3.0*drvs.df00111*by128*bygrada2
            +drvs.df0030*by256*bygrada3 +3.0*drvs.df0021*by256*bygrada3 -3.0*drvs.df00101*by128*bygrada3
            -3.0*drvs.df0020*by256*bygrada4 -3.0*drvs.df0011*by256*bygrada4 +3.0*drvs.df0010*by256*bygrada5 +drvs.df00003*by64;
  ds->fRRZ = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada +drvs.df1110*by8*bygrada +drvs.df20001*by8 +drvs.df11001*by8;
  ds->fRZZ = drvs.df10101*by16*bygrada +drvs.df10011*by16*bygrada +drvs.df1020*by64*bygrada2 +drvs.df1011*by32*bygrada2
            +drvs.df1002*by64*bygrada2 -drvs.df1010*by64*bygrada3 -drvs.df1001*by64*bygrada3 +drvs.df10002*by16;*/

  ds->fX   = drvs.df0010*by4*bygrada -drvs.df00001*by4;

  ds->fSS  = drvs.df2000*by2 -drvs.df1100*by2;
  ds->fYY  = drvs.df0020*by8*bygrada2 -drvs.df0011*by8*bygrada2 -drvs.df0010*by8*bygrada3;
  ds->fRX  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada -drvs.df10001*by4;
  ds->fSY  = drvs.df1010*by4*bygrada -drvs.df1001*by4*bygrada;
  ds->fZX  = drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3;

  ds->fRSS = drvs.df3000*by4 -drvs.df2100*by4;
  ds->fRYY = drvs.df1020*by16*bygrada2 -drvs.df1011*by8*bygrada2 +drvs.df1002*by16*bygrada2 -drvs.df1010*by16*bygrada3 -drvs.df1001*by16*bygrada3;
  ds->fSSZ = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada -drvs.df1110*by8*bygrada +drvs.df20001*by8 -drvs.df11001*by8;
  ds->fZYY =                                                         drvs.df0030*by64*bygrada3 -drvs.df0021*by64*bygrada3
                                        -3.0*drvs.df0020*by64*bygrada4 +drvs.df0011*by64*bygrada4 +3.0*drvs.df0010*by64*bygrada5;
  ds->fRSY = drvs.df2010*by8*bygrada -drvs.df2001*by8*bygrada;
  ds->fSZY =                                                       drvs.df1020*by32*bygrada2 -drvs.df1002*by32*bygrada2
            -drvs.df1010*by32*bygrada3 +drvs.df1001*by32*bygrada3;

/*ds->fX   = drvs.df0010*by4*bygrada -drvs.df00001*by4;

  ds->fSS  = drvs.df2000*by2 -drvs.df1100*by2;
  ds->fYY  = drvs.df0020*by8*bygrada2 -drvs.df0011*by8*bygrada2 -drvs.df0010*by8*bygrada3;
  ds->fXX  = drvs.df0020*by32*bygrada2 -drvs.df00101*by8*bygrada +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 +drvs.df00002*by16;
  ds->fRX  = drvs.df1010*by8*bygrada +drvs.df1001*by8*bygrada -drvs.df10001*by4;
  ds->fSY  = drvs.df1010*by4*bygrada -drvs.df1001*by4*bygrada;
  ds->fZX  = drvs.df0020*by32*bygrada2 +drvs.df0011*by32*bygrada2 -drvs.df0010*by32*bygrada3 -drvs.df00002*by16;

  ds->fXXX = 3.0*drvs.df00102*by64*bygrada -3.0*drvs.df00201*by128*bygrada2 -3.0*drvs.df00111*by128*bygrada2
            +drvs.df0030*by256*bygrada3 +3.0*drvs.df0021*by256*bygrada3 +3.0*drvs.df00101*by128*bygrada3
            -3.0*drvs.df0020*by256*bygrada4 -3.0*drvs.df0011*by256*bygrada4 +3.0*drvs.df0010*by256*bygrada5 -drvs.df00003*by64;
  ds->fRRX = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada +drvs.df1110*by8*bygrada -drvs.df20001*by8 -drvs.df11001*by8;
  ds->fRSS = drvs.df3000*by4 -drvs.df2100*by4;
  ds->fRYY = drvs.df1020*by16*bygrada2 -drvs.df1011*by8*bygrada2 +drvs.df1002*by16*bygrada2 -drvs.df1010*by16*bygrada3 -drvs.df1001*by16*bygrada3;
  ds->fRXX = drvs.df1020*by64*bygrada2 -drvs.df10101*by16*bygrada -drvs.df10011*by16*bygrada +drvs.df1011*by32*bygrada2
            +drvs.df1002*by64*bygrada2 -drvs.df1010*by64*bygrada3 -drvs.df1001*by64*bygrada3 +drvs.df10002*by16;
  ds->fSSZ = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada -drvs.df1110*by8*bygrada +drvs.df20001*by8 -drvs.df11001*by8;
  ds->fSSX = drvs.df2010*by16*bygrada +drvs.df2001*by16*bygrada -drvs.df1110*by8*bygrada -drvs.df20001*by8 +drvs.df11001*by8;
  ds->fZZX = drvs.df00201*by128*bygrada2 -drvs.df00102*by64*bygrada +drvs.df00111*by128*bygrada2 +drvs.df0030*by256*bygrada3
            +3.0*drvs.df0021*by256*bygrada3 -drvs.df00101*by128*bygrada3 -3.0*drvs.df0020*by256*bygrada4
            -3.0*drvs.df0011*by256*bygrada4 +3.0*drvs.df0010*by256*bygrada5 -drvs.df00003*by64;
  ds->fZYY = drvs.df00201*by32*bygrada2 -drvs.df00111*by32*bygrada2 +drvs.df0030*by64*bygrada3 -drvs.df0021*by64*bygrada3
            -drvs.df00101*by32*bygrada3 -3.0*drvs.df0020*by64*bygrada4 +drvs.df0011*by64*bygrada4 +3.0*drvs.df0010*by64*bygrada5;
  ds->fZXX = drvs.df0030*by256*bygrada3 -drvs.df00102*by64*bygrada -drvs.df00201*by128*bygrada2 -drvs.df00111*by128*bygrada2
            +3.0*drvs.df0021*by256*bygrada3 +drvs.df00101*by128*bygrada3 -3.0*drvs.df0020*by256*bygrada4
            -3.0*drvs.df0011*by256*bygrada4 +3.0*drvs.df0010*by256*bygrada5 +drvs.df00003*by64;
  ds->fYYX = drvs.df00111*by32*bygrada2 -drvs.df00201*by32*bygrada2 +drvs.df0030*by64*bygrada3 -drvs.df0021*by64*bygrada3
            +drvs.df00101*by32*bygrada3 -3.0*drvs.df0020*by64*bygrada4 +drvs.df0011*by64*bygrada4 +3.0*drvs.df0010*by64*bygrada5;
  ds->fRSY = drvs.df2010*by8*bygrada -drvs.df2001*by8*bygrada;
  ds->fRZX = drvs.df1020*by64*bygrada2 +drvs.df1011*by32*bygrada2 +drvs.df1002*by64*bygrada2 -drvs.df1010*by64*bygrada3
            -drvs.df1001*by64*bygrada3 -drvs.df10002*by16;
  ds->fSZY = drvs.df10101*by16*bygrada -drvs.df10011*by16*bygrada +drvs.df1020*by32*bygrada2 -drvs.df1002*by32*bygrada2
            -drvs.df1010*by32*bygrada3 +drvs.df1001*by32*bygrada3;
  ds->fSYX = drvs.df10011*by16*bygrada -drvs.df10101*by16*bygrada +drvs.df1020*by32*bygrada2 -drvs.df1002*by32*bygrada2
            -drvs.df1010*by32*bygrada3 +drvs.df1001*by32*bygrada3;*/
}

extern void quit_(const char* str, int len);
void
dalton_quit(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;
 
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
#ifdef PRG_DIRAC
    quit_(line, len);
#endif
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
#ifdef INT_STAR8
extern void fortwrt_(const char *line, const long *linelen);
#else
extern void fortwrt_(const char *line, const int *linelen);
#endif /* ifdef INT_STAR8 */
int
fort_print(const char* format, ...)
{
    char line[128];
#if defined(INT_STAR8)
    long len; 
#else 
    int len;
#endif
    va_list a;

/* DEC does not recognize the vsnprintf command */
/* unless when linked with -ldb */
#ifndef SYS_DEC
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    fortwrt_(line, &len);
    return len;
#endif
}

/* ################################################################### */
/* MPI related functions                                               */
/* ################################################################### */
/* dft_sync_func:
   synchronizes functional data between nodes.
   Needs to be called only once (after functional initialization/change).
   NO-OP for serial code.
   In any case, do it only once, otherwise the functional becomes
   a sum of functionals.
*/
void
dftsyncfunc_(int* amImaster)
{
#ifdef VAR_MPI
    static int Initialized = 0;
    int success; /* jvs: this is not used !!! */
    int dft_conf_len =  DftConfString ? strlen(DftConfString)+1 : 0;
    
    MPI_Bcast(&dcrham_.hfxmu,1,MPI_DOUBLE,diracinfpar_.mparid,MPI_COMM_WORLD);
    if(diracinfpar_.nodes==0||Initialized) return; 
    MPI_Bcast(&do_alda_hs,   1, MPI_INT, diracinfpar_.mparid, MPI_COMM_WORLD);
    MPI_Bcast(&do_alda_ha,   1, MPI_INT, diracinfpar_.mparid, MPI_COMM_WORLD);
    MPI_Bcast(&do_xalda,     1, MPI_INT, diracinfpar_.mparid, MPI_COMM_WORLD);
    MPI_Bcast(&dft_conf_len, 1, MPI_INT, diracinfpar_.mparid, MPI_COMM_WORLD);
    if(dft_conf_len>0) {
        char* line = malloc(dft_conf_len);
        if(*amImaster)
            strcpy(line, DftConfString);
       
        MPI_Bcast(line, dft_conf_len, MPI_CHAR, diracinfpar_.mparid, MPI_COMM_WORLD);

        if(!*amImaster)
            dftinput_(line, &success, dft_conf_len);

        free(line);
        Initialized = 1;
    }
#endif /* VAR_MPI */
}
