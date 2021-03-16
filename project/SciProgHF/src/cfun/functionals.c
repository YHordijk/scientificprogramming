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
/* functionals.c:
   Program for computing and testing functional routines in the DFT module.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02

 */

#define _DEFAULT_SOURCE 1

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define __CVERSION__

#include "functionals.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &KTFunctional,
    &LB94Functional,
    &LBalphaFunctional,
    &GLLBholeFunctional,
    &LYPFunctional,
    &OPTXFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &PW91xFunctional,
    &PW91cFunctional,
    &PZ81Functional,
    &PBEcFunctional,
    &PBExFunctional,
    &RpbexFunctional,
    &RevpbexFunctional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNIFunctional,
    &VWNFunctional,
    &XAlphaFunctional,
    /* mixed functionals */
    &B3LYPFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86GFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &PP86Functional,
    &BPW91Functional,
    &Camb3lypFunctional,
    &GGAKeyFunctional,
    &CamFunctional,
    &HseFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &OLYPFunctional,
    &PBE0Functional,
    &PBE38Functional,
    &PBEFunctional,
    &RPBEFunctional,
    &SVWN3Functional,
    &SVWN5Functional,
    /* range-separated fnuctionals */
    &SRVWN5Functional,
    &HjsxFunctional,
    &SrpbecFunctional,
    NULL
};

/* kinetic energy functionals */
Functional* available_kin_functionals[] = {
    &TF_KinFunctional,
    &VW_KinFunctional,
    NULL
};

static int my_printf(const char *fmt, ...)
{
    int i;va_list ap; va_start(ap, fmt); i= vprintf(fmt, ap); va_end(ap);
    puts("");
    return i;
}
 
static void set_hf_weight(real w)         {}
static real get_hf_weight(void)           {return 0;}
static void set_cam_param(real w, real b) {}

Functional* selected_func = &LDAFunctional;
Functional* emb_xc_func   = &LDAFunctional;
Functional* emb_kin_func  = NULL;
int (*fun_printf)(const char *fmt, ...) = my_printf;
void (*fun_set_hf_weight)(real w)         = set_hf_weight;
real (*fun_get_hf_weight)(void)           = get_hf_weight;
void (*fun_set_cam_param)(real w, real b) = set_cam_param;

/* =================================================================== */
void
drv1_clear(FunFirstFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv2_clear(FunSecondFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv3_clear(FunThirdFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

enum FunError
fun_select_by_name(const char *conf_string)
{
    int ok, i;
    char func_name[20];

    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name)) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

int fun_true(void)  { return 1; }
int fun_false(void) { return 0; }

/* Fortran interface. We specify different names suffixes so that
 * library can be linked with code compiled with different compilers
 * or different compilation options. */

void
funset(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}
void
funset_(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}

int
funisgga(void)
{ return selected_func->is_gga(); }
int
funisgga_(void)
{ return selected_func->is_gga(); }

/* end of backward compatiblity section */

int
dftisemb_(void)
{
    return !strcmp(selected_func->name, "embedding");
    /* or return selected_func == EmbFunctional */
}

real
funenergy(const FunDensProp *dp)
{ return selected_func->func(dp); }
real
funenergy_(const FunDensProp *dp)
{ return selected_func->func(dp); }

void
funfirst(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }
void
funfirst_(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }

void
funlistfuncs(Functional *l[])
{
    int i;
    
    fun_printf("\nAvailable functionals:");
    for(i=0; l[i]; i++)
        fun_printf(l[i]->name);
}

/* =================================================================== */
/*           fortran (and not only) functional stub routines           */
/* =================================================================== */

/* dftreport:
   report the selected functional and its configuration.
*/
void
dftreport_(void)
{
    fun_printf("* Kohn-Sham calculation using the xc functional: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();
}

void
dftlistfuncs_(void)
{
    funlistfuncs(available_functionals);
}

/* declare both known fortran name-mangled variants */
int
dftisgga_(void)
{ return selected_func->is_gga(); }

int
dftisgga__(void)
{ return selected_func->is_gga(); }

real
dftenergy_(const real* rho, const real* grad)
{
    FunDensProp dp;
    /* HP does not grok { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    if (*rho > 1e-20)	
        dp.rhoa  = dp.rhob  = *rho *0.5;
    else dp.rhoa = dp.rhob = 1e-20;
    if (*grad > 1e-20) 
        dp.grada = dp.gradb = *grad*0.5;
    else dp.grada = dp.gradb = 1e-20; 
    dp.gradab= dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    return selected_func->func(&dp);
}
