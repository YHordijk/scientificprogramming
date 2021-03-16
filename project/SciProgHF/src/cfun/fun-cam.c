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
/** @file fun-cam.c General CAM functional. 

    Often called a range-separated exchange method.

    Pawel Salek, 2004.06, Himmelbjerg-initial implementation.
*/

 
/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
/* Use BSD's strncasecmp() */
#define _DEFAULT_SOURCE 1

#include <assert.h>
#include <math.h>

#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
 
#define __CVERSION__
 
#include "functionals.h"

#define ELEMENTS(arr) (sizeof(arr)/sizeof(arr[0]))

/* RGFirstDrv, RGSecondDrv and RGThirdDrv hold derivatives of a function
 * with respect to density R (R) and density gradient (G).
 * Field dfAB contains d^A/dR^a d^B/dG^B f.
 */
typedef struct {
    real df10, df01;
} RGFirstDrv;
typedef struct {
    real df10, df01, df20, df11, df02;
} RGSecondDrv;
typedef struct {
    real df10, df01, df20, df11, df02;
    real df30, df21, df12, df03;
} RGThirdDrv;

 
/* INTERFACE PART */
#define THR 1e-40

static int camb3lyp_read(const char *conf_line);
static void camb3lyp_report(void);

static int cam_isgga(void) { return 1; }
static int cam_read(const char *conf_line);
static void cam_report(void);
static real cam_energy(const FunDensProp* dp);
static void cam_first(FunFirstFuncDrv *ds,   real factor,
                      const FunDensProp* dp);
static void cam_second(FunSecondFuncDrv *ds, real factor,
                       const FunDensProp* dp);
static void cam_third(FunThirdFuncDrv *ds, real factor,
                      const FunDensProp* dp);

static int hse_read(const char *conf_line);
static void hse_report(void);

Functional Camb3lypFunctional = {
  "Camb3lyp",       /* name */
  cam_isgga,   /* gga-corrected */
  camb3lyp_read,
  camb3lyp_report,
  cam_energy,
  cam_first,
  cam_second,
  cam_third
};

Functional CamFunctional = {
  "Cam",       /* name */
  cam_isgga,   /* gga-corrected */
  cam_read,
  cam_report,
  cam_energy,
  cam_first,
  cam_second,
  cam_third
};

Functional HseFunctional = {
  "HSE",       /* name */
  cam_isgga,   /* gga-corrected */
  hse_read,
  hse_report,
  cam_energy,
  cam_first,
  cam_second,
  cam_third
};
 
/* IMPLEMENTATION PART */
struct FunctionalList {
    Functional *func;
    real weight;
    struct FunctionalList *next;
};

/** The module uses program-wide configuration. It uses following
    range separation of the HF exchange:
    HF_RS_Exch = (alpha + beta*erf(mu*r))*HF_exchange

    This means that the DFT exchange becomes:
    1 - HF_RS_Exch
*/ 

static struct FunctionalList *exchangeFunctionals = NULL;
static struct FunctionalList *correlationFunctionals = NULL;
static real camAlpha = 0.0; 
static real camBeta  = 0.0;
static real camMu    = 0.0;

static struct FunctionalList*
newFunc(Functional *f, real weight,
        struct FunctionalList *next)
{
    struct FunctionalList *t = (struct FunctionalList*)
        malloc(sizeof(struct FunctionalList));
    t->func = f;
    t->weight = weight;
    t->next = next;
    return t;
}

static void
cam_free_config(void)
{
    struct FunctionalList *fi, *t;

    for(fi=exchangeFunctionals; fi; fi = t) {
        t = fi->next;
        free(fi);
    }
    for(fi=correlationFunctionals; fi; fi = t) {
        t = fi->next;
        free(fi);
    }
    exchangeFunctionals = NULL;
    correlationFunctionals = NULL;
}

/* CAM-B3LYP part. */
static int
parse_table(const char *func, const char *str,
            int cnt, const char *keywords[], double *weights)
{
  int res=1;
  while(*str) {
    int i;
    while(*str && isspace((int)*str)) str++; /* skip whitespace */
    if(*str =='\0') break; /* line ended by whitespace */
    for(i=0; i<cnt; i++) {
      int len = strlen(keywords[i]);
      if(strncasecmp(keywords[i], str, len)==0 &&
         str[len] == '=') {
        if(sscanf(str+len+1,"%lg", &weights[i]) != 1) {
          fun_printf("%s: %s not followed by the weight: ",
                     func, keywords[i]);
          res = 0;
        }
        break;
      }
    }
    if(i==cnt) {
      fun_printf("%s: unknown string: '%s'", func, str);
      res = 0;
    }
    while(*str && !isspace((int)*str)) str++; /* skip nonws */
  }
  return res;
}

#if 1
/* yanai-consistent */
#define BECKE88_CORR_WEIGHT 1
#else
/* B3-lyp consistent */
#define BECKE88_CORR_WEIGHT 0.9
#endif
#define LYP_WEIGHT 0.81
#define VWN_WEIGHT 0.19

static const char *cam_keywords[] = { "alpha", "beta", "mu" };
static int
camb3lyp_read(const char *conf_line)
{
    return cam_read("p:alpha=0.19 p:beta=0.46 p:mu=0.33 x:slater=1 x:becke=1 c:lyp=0.81 c:vwn5=0.19");
}


static void
camb3lyp_report(void)
{
    fun_printf("  T. Yanai, D. Tew and N. C. Handy, Chem. Phys. Lett., 393 (2004) 51");
    cam_report();
    return;
}

/* General CAM part */

/** Read the configuration. The configuration consists of three types
    of terms that follow general pattern:

    (p|x|c):(FUNCTIONAL)=(weight)

    p prefix is followed by either 'alpha', 'beta' or 'mu' parameters
    and corresponding weights. x prefix defines an exchange functional
    - no actual check is performed! c allows to add an additive
    functional, usually a correlation one.

    Example configuration for CAM-B3LYP is:

    CAM p:alpha=0.19 p:beta=0.46 p:mu=0.33 x:slater=1 x:becke=1 
        c:lyp=0.81 c:vwn5=0.19

    We obviously need to carefully exclude the recursive case of cam
    functional built from another cam functional....

    @returns 0 on failure, 1 on success.
*/
    
static int
cam_read(const char *conf_line)
{
    static const struct {
        const char *s;
        real *p;
    } Params[] = {
        { "alpha=", &camAlpha },
        { "beta=",  &camBeta },
        { "mu=",    &camMu }
    };
    const char * str;
    struct FunctionalList *t, **funcList;

    cam_free_config();

    for(str = conf_line; *str; ) {
        int i;
        char prefix;
        while(*str && isspace((int)*str))
            str++; /* skip whitespace */

        if(*str =='\0') break; /* line ended by whitespace */

        prefix = *str;
        if(*++str != ':') {
            fun_printf("':' is expected to follow the one-character long "
                       "type prefix '%c'.", prefix);
            return 0;
        }

        if(*++str == '\0') {
            fun_printf("FUNC=WEIGHT is expected to follow '%c:' prefix.",
                       prefix);
            return 0;
        }

        if(prefix == 'p') {
            unsigned p;
            int l;
            for(p=0; p<ELEMENTS(Params); p++) {
                l = strlen(Params[p].s);
                if (strncmp(str, Params[p].s, l) == 0)
                    break;
            }
            if (p==ELEMENTS(Params)) {
                fun_printf("unknown parameter %s following p: prefix.",
                           str);
                return 0;
            }

            if (sscanf(str + l, "%lf", Params[p].p) != 1) {
                fun_printf("cannot extract weight from %s.", str);
                return 0;
            }
        } else {
            if (prefix == 'x') 
                funcList = &exchangeFunctionals;
            else if (prefix == 'c')
                funcList = &correlationFunctionals;
            else {
                fun_printf("unknown prefix '%c' in %s.", prefix, str);
                return 0;
            }
            for(i=0; available_functionals[i]; i++) {
                real weight;
                int len;
                if(available_functionals[i] == &CamFunctional)
                    continue; /* Prevent functional definition recursion. */
                
                len = strlen(available_functionals[i]->name);
                if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
                   str[len] == '=') {
                    if(sscanf(str+len+1,"%lf", &weight) != 1) {
                        fun_printf("CAM: %s not followed by a weight.",
                                   available_functionals[i]->name);
                        return 0;
                    }

                    *funcList = newFunc(available_functionals[i], weight,
                                        *funcList);
                    break;
                }
            }

            if(!available_functionals[i]) {
                fun_printf("CAM: unknown functional: '%s'", str);
                return 0;
            }
        }

        while(*str && !isspace((int)*str)) str++; /* skip nonws */
    }
    fun_set_hf_weight(camAlpha);
    fun_set_cam_param(camMu, camBeta);
    return 1;
}


static void
cam_report(void)
{
    struct FunctionalList *l;
    fun_printf("  CAM functional with alpha=%5.3g beta=%5.3g mu=%5.3g",
               camAlpha, camBeta, camMu);
    fun_printf("    and exchange functionals and weights:");
    for(l=exchangeFunctionals; l; l = l->next)
        fun_printf("%20s : %lg", l->func->name, l->weight);
    fun_printf("    and correlation functionals and weights:");
    for(l=correlationFunctionals; l; l = l->next)
        fun_printf("%20s : %lg", l->func->name, l->weight);
    fun_printf("  Short-range exchange functionals are obtained according to:");
    fun_printf("    H. Iikura, T. Tsuneda, T. Yanai, K. Hirao, J. Chem. Phys. 115(2001) 3540");
}

static int
hse_read(const char *conf_line)
{
    /* HSE functional is defined as in J. Heyd, G.E. Scuseria and
       M. Ernzerhof, J. Chem. Phys. 118 (2003), p. 8207.
       The implementation is not fully tested and there is no regression test.
       H atom with 6-311++G(3df,3pd) basis set gives          -0.50368291542.
       With HSE's erratum [J. Chem. Phys. 124 219906 (2006)]: -0.51377005236
       With HSE's erratum and mu_HF:= mu_PBE/sqrt(2):         -0.51020009241
    */
    return cam_read("p:alpha=0.25 p:beta=-0.25 p:mu=0.1890 c:pbec=1 x:pbex=1");
}
static void
hse_report(void)
{
    fun_printf("  J. Heyd, G.E. Scuseria and M. Ernzerhof, J. Chem. Phys. 118 (2003) 8207");
    cam_report();
    return;
}

/* ===================================================================
 * compute a and its derivatives.
 * =================================================================== */
/* consider usign M_2_SQRTPI defined as 2/sqrt(pi) */
#define SQRT_PI 1.77245385090552
static real
fun_a(real rho, real ex)
{
    real a = camMu*sqrt(-2.0*ex)/(6.0*SQRT_PI*rho);
    return a;
}

static void
fun_a_first(real rho, real a, real ex, RGFirstDrv *fun1, 
            RGFirstDrv *res)
{
    memset(res, 0, sizeof(RGFirstDrv));
    if(fabs(a)<THR || fabs(ex)<THR) return;
    res->df10 = (fun1->df10/(2*ex)-1.0/rho)*a;
    res->df01 = fun1->df01/(2*ex)*a;
}

static void
fun_a_second(real rho, real a, real ex, RGSecondDrv *f2,
             RGSecondDrv *res)
{
    real f10, f20, f01;

    memset(res, 0, sizeof(RGSecondDrv));
    if(fabs(a)<THR || fabs(ex)<THR) return;
    f10 = f2->df10/(2*ex)-1.0/rho;
    f20 = f2->df20/(2*ex) - f2->df10*f2->df10/(2*ex*ex) +1/(rho*rho);
    f01 = f2->df01/(2*ex);
    res->df10 = a*f10;
    res->df01 = a*f01;
    res->df20 = a*(f10*f10 + f20);
    res->df11 = a*(f10*f01 + f2->df11/(2*ex) -
                   f2->df10*f2->df01/(2*ex*ex));
    res->df02 = a*(f2->df02/(2*ex)-f2->df01*f2->df01/(4*ex*ex)); 
}

static void
fun_a_third(real rho, real a, real ex, RGThirdDrv *f3, RGThirdDrv *res)
{
    real f10, f01, f20, f11, f02, f30, f21, f12, f03;

    memset(res, 0, sizeof(RGThirdDrv));
    if(fabs(a)<1e-15 || fabs(ex)<1e-15) return;
    f10 = f3->df10/(2*ex)-1.0/rho;
    f01 = f3->df01/(2*ex);
    f20 = f3->df20/(2*ex) - f3->df10*f3->df10/(2*ex*ex) +1/(rho*rho);
    f11 = f3->df11/(2*ex) - f3->df10*f3->df01/(2*ex*ex);
    f02 = f3->df02/(2*ex) - f3->df01*f3->df01/(2*ex*ex);
    f30 = f3->df30/(2*ex) - 1.5*f3->df10*f3->df20/(ex*ex)
          + f3->df10*f3->df10*f3->df10/(ex*ex*ex) - 2/(rho*rho*rho);
    f21 = f3->df21/(2*ex)
        - f3->df20*f3->df01/(2*ex*ex) - 2*f3->df10*f3->df11/(2*ex*ex)
        + f3->df10*f3->df10*f3->df01/(ex*ex*ex);
    f12 = f3->df12/(2*ex)
        - 2*f3->df11*f3->df01/(2*ex*ex) - f3->df10*f3->df02/(2*ex*ex) 
        + f3->df10*f3->df01*f3->df01/(ex*ex*ex);
    f03 = f3->df03/(2*ex) - 1.5*f3->df02*f3->df01/(ex*ex)
        + f3->df01*f3->df01*f3->df01/(ex*ex*ex);
    res->df10 = a*f10;
    res->df01 = a*f01;
    res->df20 = a*(f10*f10 + f20);
    res->df11 = a*(f10*f01 + f11);
    res->df02 = a*(f01*f01 + f02);

    /* terms will partially cancel out - think how to simplify them */
    res->df30 = a*(f10*f10*f10 + 3*f10*f20 + f30);
    res->df21 = a*(f10*f10*f01 + f20*f01 + 2*f10*f11 + f21);
    res->df12 = a*(f10*f01*f01 + 2*f11*f01 + f10*f02 + f12);
    res->df03 = a*(f01*f01*f01 + 3*f01*f02 + f03);
}

/* optimized with help of maxima */

/* ===================================================================
 * The expansion of the B-factor and its first derivative for small
 * values of a.
 * =================================================================== */
static real
cam_b_energy_small(real a)
{
    real res;
    a = 2*a; /* the expension derived for different a; correct for this. */
    res = 1-4.0/3.0*SQRT_PI*a + 2 * a*a - 2.0/3.0*a*a*a*a;
    return 1-res;
}


static real
cam_b_first_small(real a)
{
    real res;
    a = 2*a; /* the expension derived for different a; correct for this. */
    res = 4.0/3.0*(-SQRT_PI + 3 * a +(2*exp(-1/(a*a)) - 2.0)*a*a*a);
    return 2*res;
}

/* ===================================================================
 * The expansion of the B factor and its first derivative for large
 * values of a.
 * =================================================================== */
#define MAX_LARGE_COEFS 5
static const real large_coefs[] = { 9, 60, 420, 3240, 27720 };
static real
cam_b_energy_large(real a)
{
    real res, ac, a2;
    int i;
    
    a = 2*a; /* the expension derived for different a; correct for this. */
    a2 = a*a;
    res = 0;
    for(i=0, ac = a2; i<MAX_LARGE_COEFS; i++, ac *= -a2)
        res += 1.0/(large_coefs[i]*ac);
    return 1-res;
}


static const real large_coefs1[] = { 4.5, 15, 70, 405 };
static real
cam_b_first_large(real a)
{
    real tmp;
    real ac, a2;
    int i;

    a = 2*a; /* the expension derived for different a; correct for this. */
    a2  = a*a;
    tmp = 0;
    for(i=0, ac = -a2*a; i<MAX_LARGE_COEFS-1; i++) {
        tmp += 1.0/(large_coefs1[i]*ac);
        ac *= -a2;
    }
    return 2*tmp;
}

/* ===================================================================
 * The expansion of the B factor and its first, second and third
 * derivatives for medium values of a. This part uses the full
 * expression.
 * 8/3*a*(sqrt(%PI)*erf(1/(2*a))+2*a*(b-c))
 * =================================================================== */
static real
cam_b_energy_medium(real a)
{
    real b = exp(-1/(4*a*a))-1;
    real c = 2*a*a*b + 0.5;
    real res= 8.0/3.0*a*(SQRT_PI*erf(1/(2*a))+2*a*(b-c));
    return res;
}

/* tested version */
static real
cam_b_first_medium(real a)
{
    real t1 = 1/a;
    real t2 = a*a;
    real t3 = 1/t2;
    real t4 = exp(-0.25*t3);
    real t5 = t4-1;
    real t6 = t4-2*t2*t5-1.5;
    real res = -2.666666666666667*a
        *(2*a*(t4/(2*pow(a,3.0))-4*a*t5-t1*t4)+2*t6-t3*t4)
        -2.666666666666667*(2*a*t6+SQRT_PI*erf(0.5*t1));
    return res;
}

static real
evaluate_series(int n, const real*coefs, real lambda)
{
    real res = 0, ac =1.0;
    int i;
    for(i=0; i<n; i++) {
        res += 1.0/(coefs[i]*ac);
        ac *= lambda;
    }
    return res;
}
static real
cam_b_second_medium(real a)
{
    real t1, t2, t3, t4, t5, t6, t7;

    if(fabs(a)<THR) return 16.0;
    t1 = a*a;
    if(a>=5)  {
        static const real large_coefs[] = {  6, -48, 640,  -11520 };
        return evaluate_series(ELEMENTS(large_coefs), large_coefs, t1)/
            (t1*t1);
    }
    t2 = 1/t1;
    t3 = exp(-0.25*t2);
    t4 = pow(a,-3.0);
    t5 = t3-1.0;
    t6 = -t2*t3;
    t7 = -t3/a+0.5*t4*t3-4*a*t5;
    return -(8*a*(2*a*(t3/(4*pow(a,6.0))-2*t3*pow(a,-4.0)+t6-4*t5)
                  -t3/(2*pow(a,5.0))+4*t7+2*t4*t3)/3.0
             +16*(2.0*a*t7+2.0*(t3-2*t1*t5-1.5)+t6)/3.0);
}

static real
cam_b_third_medium(real a)
{
    real t1, t2, t3, t4, t5, t6, t7, t8;

    if(fabs(a)<THR) return 0;
    if(a>=5)  {
        static const real large_coefs[] = {  -1.5, 8, -80,  1152 };
        real a2 = a*a;
        return evaluate_series(ELEMENTS(large_coefs), large_coefs, a2)/
            (a2*a2*a);
    }
    t1 = pow(a,-2.0);
    t2 = exp(-0.25*t1);
    t3 = pow(a,-6.0);
    t4 = pow(a,-4.0);
    t5 = pow(a,-5.0);
    t6 = t2-1;
    t7 = -t1*t2-2*t4*t2+t3*t2/4-4*t6;
    t8 = pow(a,-3.0);
    return -8*
        (a*(2*a*(t2/(8*pow(a,9.0))-5*t2/(2*pow(a,7.0))+15*t5*t2/2)
                 -t2/(4*pow(a,8.0))+6*t7-6*t4*t2+7*t3*t2/2)/3
         +(4*(-t2/a+t8*t2/2-4*a*t6)+2*a*t7+2*t8*t2-t5*t2/2));
}


#define FAC M_SQRT2
#define EVALUATOR(a,type) \
  ((a<0.14) ? cam_b_ ## type ## _small(a)  : \
  ((a<4.25) ? cam_b_ ## type ## _medium(a) : \
   cam_b_ ## type ## _large(a)))

static real
cam_energy_sigma(real rho, real ex)
{
    real a       = fun_a(rho, ex);
    real bfactor = EVALUATOR(a, energy);
    return ex*bfactor;
}

/* Exchange energy should be always negative. Some functionals like
   PBE are exteremely numerically sensitive and contributions that
   should given exactly 0 do no. In general we should implement such
   functionals more carefully. For now, we just stop if we detect a
   suspicious case. */

static const real ENERGY_THR = 1e-13;

static real
cam_energy(const FunDensProp *dp)
{
    real res, ea, eb, ex;
    FunDensProp dsigma;
    struct FunctionalList *f;

    dsigma.rhoa  = dsigma.rhob  = dp->rhoa;
    dsigma.grada = dsigma.gradb = dp->grada;
    dsigma.gradab = dsigma.grada*dsigma.gradb;

    
    for(ex = 0.0, f=exchangeFunctionals; f; f = f->next)
        ex += 0.5 * f->func->func(&dsigma) * f->weight;
    if(ex>0) {
        assert(ex < ENERGY_THR);
        return 0.0;
    }

    ea = ex*(1-camAlpha) - camBeta * cam_energy_sigma(dp->rhoa, ex);

    if(fabs(dp->rhoa-dp->rhob)>THR || fabs(dp->grada-dp->gradb)>THR) {
        dsigma.rhoa  = dsigma.rhob  = dp->rhob;
        dsigma.grada = dsigma.gradb = dp->gradb;
	dsigma.gradab = dsigma.grada*dsigma.gradb;
        
        for(ex = 0.0, f=exchangeFunctionals; f; f = f->next)
            ex += 0.5 * f->func->func(&dsigma) * f->weight;

        if(ex>0) {
            assert(ex < ENERGY_THR);
            return 0.0;
        }
        eb = ex*(1-camAlpha) - camBeta * cam_energy_sigma(dp->rhob, ex);
    } else eb = ea;
    res = ea + eb;

    
    for(f=correlationFunctionals; f; f = f->next)
        res += f->func->func(dp) * f->weight;

    return res;
}

static void
cam_first_sigma(real rho, real ex,
                     RGFirstDrv *ds, RGFirstDrv *res)
{
    real             a = fun_a(rho, ex);
    real       bfactor = -camBeta*EVALUATOR(a, energy);
    real bfactor_first =  camBeta*EVALUATOR(a, first);
    RGFirstDrv ader;

    fun_a_first(rho, a, ex, ds, &ader);

    res->df10 = ds->df10*bfactor + ex*bfactor_first*ader.df10;
    res->df01 = ds->df01*bfactor + ex*bfactor_first*ader.df01;
}

static void
cam_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    RGFirstDrv res, dfun;
    FunDensProp dsigma;
    real ex, weight_lr = 1-camAlpha;
    FunFirstFuncDrv fun1;
    struct FunctionalList *f;

    dsigma.rhoa  = dsigma.rhob  = dp->rhoa;
    dsigma.grada = dsigma.gradb = dp->grada;
    dsigma.gradab = dsigma.grada*dsigma.gradb;

    memset(&fun1, 0, sizeof(fun1));

    for(ex = 0.0, f=exchangeFunctionals; f; f = f->next) {
        ex += 0.5 * f->func->func(&dsigma) * f->weight;
        f->func->first(&fun1, f->weight, dp);
    }
    if(ex>0) {
        assert(ex < ENERGY_THR);
        return;
    }

    dfun.df10 = fun1.df1000; dfun.df01 = fun1.df0010;
    cam_first_sigma(dp->rhoa, ex, &dfun, &res);
    ds->df1000 += factor*(weight_lr*dfun.df10 + res.df10);
    ds->df0010 += factor*(weight_lr*dfun.df01 + res.df01);

    if(dp->rhob>10e-13) {
        if(fabs(dp->rhoa-dp->rhob)>THR || fabs(dp->grada-dp->gradb)>THR) {
            dsigma.rhoa  = dsigma.rhob  = dp->rhob;
            dsigma.grada = dsigma.gradb = dp->gradb;
	    dsigma.gradab = dsigma.grada*dsigma.gradb;

            for(ex = 0.0, f=exchangeFunctionals; f; f = f->next)
                ex += 0.5 * f->func->func(&dsigma) * f->weight;

            dfun.df10 = fun1.df0100; dfun.df01 = fun1.df0001;
            if(ex>0) {
                assert(ex < ENERGY_THR);
                return;
            }

            cam_first_sigma(dp->rhob, ex, &dfun, &res);
        }
        ds->df0100 += factor*(weight_lr*dfun.df10 + res.df10);
        ds->df0001 += factor*(weight_lr*dfun.df01 + res.df01);
    }

    for(f=correlationFunctionals; f; f = f->next)
        f->func->first(ds, f->weight * factor, dp);

}

static void
cam_second_sigma(real rho, real ex, RGSecondDrv *f2,
                      RGSecondDrv *res)
{
    real bfactor, b_first, b_second;
    RGSecondDrv ader;
    real a;
    
    if(rho<1e-13) {
        res->df10 = res->df01 = 0;
        res->df20 = res->df11 =  res->df02 = 0;
        return;
    }

    a = fun_a(rho, ex);
    bfactor  = -camBeta*EVALUATOR(a, energy);
    b_first  =  camBeta*EVALUATOR(a, first);
    b_second =  camBeta*cam_b_second_medium(a);

    fun_a_second(rho, a, ex, f2, &ader);
    
    res->df10 = f2->df10*bfactor + ex*b_first*ader.df10;
    res->df01 = f2->df01*bfactor + ex*b_first*ader.df01;

    res->df20 = f2->df20*bfactor + 2*f2->df10*b_first*ader.df10 +
        ex*b_second*ader.df10*ader.df10 + ex*b_first*ader.df20;
    res->df11 = f2->df11*bfactor + f2->df10*b_first*ader.df01 + 
        f2->df01*b_first*ader.df10 + 
        ex*b_second*ader.df10*ader.df01 +
        ex*b_first*ader.df11;
    res->df02 = f2->df02*bfactor + 2*f2->df01*b_first*ader.df01 +
        ex*b_second*ader.df01*ader.df01 + ex*b_first*ader.df02;
}

static void
cam_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FunSecondFuncDrv f2;
    RGSecondDrv res, dfun;
    FunDensProp dsigma;
    real ex, weight_lr = 1-camAlpha;
    struct FunctionalList *f;

    dsigma.rhoa  = dsigma.rhob  = dp->rhoa;
    dsigma.grada = dsigma.gradb = dp->grada;
    dsigma.gradab = dsigma.grada*dsigma.gradb;

    memset(&f2, 0, sizeof(f2));
    for(ex = 0.0, f=exchangeFunctionals; f; f = f->next) {
        ex += 0.5 * f->func->func(&dsigma) * f->weight;
        f->func->second(&f2, f->weight, dp);
    }
    if(ex>0) {
        assert(ex < ENERGY_THR);
        return;
    }

    dfun.df10 = f2.df1000; dfun.df20 = f2.df2000; 
    dfun.df01 = f2.df0010; dfun.df02 = f2.df0020;
    dfun.df11 = f2.df1010; 
    cam_second_sigma(dp->rhoa, ex, &dfun, &res);

    ds->df1000 += factor*(weight_lr*dfun.df10 + res.df10);
    ds->df0010 += factor*(weight_lr*dfun.df01 + res.df01);
    ds->df2000 += factor*(weight_lr*dfun.df20 + res.df20);
    ds->df1010 += factor*(weight_lr*dfun.df11 + res.df11);
    ds->df0020 += factor*(weight_lr*dfun.df02 + res.df02);

    if(fabs(dp->rhoa-dp->rhob)>THR || fabs(dp->grada-dp->gradb)>THR) {
        dsigma.rhoa  = dsigma.rhob  = dp->rhob;
        dsigma.grada = dsigma.gradb = dp->gradb;
	dsigma.gradab = dsigma.grada*dsigma.gradb;

        for(ex = 0.0, f=exchangeFunctionals; f; f = f->next)
            ex += 0.5 * f->func->func(&dsigma) * f->weight;

        if(ex>0) {
            assert(ex < ENERGY_THR);
            return;
        }

        dfun.df10 = f2.df0100; dfun.df20 = f2.df0200; 
        dfun.df01 = f2.df0001; dfun.df02 = f2.df0002;
        dfun.df11 = f2.df0101; 
        cam_second_sigma(dp->rhob, ex, &dfun, &res);
    }

    ds->df0100 += factor*(weight_lr*dfun.df10 + res.df10);
    ds->df0001 += factor*(weight_lr*dfun.df01 + res.df01);
    ds->df0200 += factor*(weight_lr*dfun.df20 + res.df20);
    ds->df0101 += factor*(weight_lr*dfun.df11 + res.df11);
    ds->df0002 += factor*(weight_lr*dfun.df02 + res.df02);

    for(f=correlationFunctionals; f; f = f->next)
        f->func->second(ds, f->weight * factor, dp);
}

/* ===================================================================
   Third order derivatives specific code.
   =================================================================== */

static void
cam_third_sigma(real rho, real ex, RGThirdDrv *f3,
                     RGThirdDrv *res)
{
    real bfactor, b_first, b_second, b_third, a;
    RGThirdDrv ader;

    a = fun_a(rho, ex);
    bfactor  = -camBeta*EVALUATOR(a, energy);
    b_first  =  camBeta*EVALUATOR(a, first);
    b_second =  camBeta*cam_b_second_medium(a);
    b_third  =  camBeta*cam_b_third_medium(a);
    fun_a_third(rho, a, ex, f3, &ader);

    res->df10 = f3->df10*bfactor + ex*b_first*ader.df10;
    res->df01 = f3->df01*bfactor + ex*b_first*ader.df01;

    res->df20 = f3->df20*bfactor + 2*f3->df10*b_first*ader.df10 +
        ex*b_second*ader.df10*ader.df10 + ex*b_first*ader.df20;
    res->df11 = f3->df11*bfactor + f3->df10*b_first*ader.df01 + 
        f3->df01*b_first*ader.df10 + 
        ex*b_second*ader.df10*ader.df01 +
        ex*b_first*ader.df11;
    res->df02 = f3->df02*bfactor + 2*f3->df01*b_first*ader.df01 +
        ex*b_second*ader.df01*ader.df01 + ex*b_first*ader.df02;

    res->df30 = f3->df30*bfactor + 3*f3->df20*b_first*ader.df10
        + 3*f3->df10*b_second*ader.df10*ader.df10
        + 3*f3->df10*b_first *ader.df20 + 3*ex*b_second*ader.df10*ader.df20
        + ex*b_third*ader.df10*ader.df10*ader.df10
        + ex*b_first*ader.df30;

    res->df21 = f3->df21*bfactor + ex*b_first*ader.df21 +
        ex*b_third*ader.df10*ader.df10*ader.df01 +
        b_first*(f3->df20*ader.df01 + f3->df01*ader.df20 + 
                 2*f3->df11*ader.df10 + 2*f3->df10*ader.df11) + 
        b_second*(ex*ader.df20*ader.df01 + f3->df01*ader.df10*ader.df10 +
                  2*ex*ader.df10*ader.df11 + 2*f3->df10*ader.df10*ader.df01);

    res->df12 = f3->df12*bfactor + ex*b_first*ader.df12 +
        ex*b_third*ader.df10*ader.df01*ader.df01 +
        b_first*(f3->df02*ader.df10 + f3->df10*ader.df02 + 
                 2*f3->df11*ader.df01 + 2*f3->df01*ader.df11) +
        b_second*(ex*ader.df10*ader.df02 + f3->df10*ader.df01*ader.df01 +
                  2*ex*ader.df01*ader.df11 + 2*f3->df01*ader.df10*ader.df01);

    res->df03 = f3->df03*bfactor + 3*f3->df02*b_first*ader.df01
        + 3*f3->df01*b_second*ader.df01*ader.df01
        + 3*f3->df01*b_first *ader.df02 + 3*ex*b_second*ader.df01*ader.df02
        + ex*b_third*ader.df01*ader.df01*ader.df01
        + ex*b_first*ader.df03;
}


static void
cam_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    FunThirdFuncDrv f3;
    RGThirdDrv res, dfun;
    FunDensProp dsigma;
    real ex, weight_lr = 1-camAlpha;
    struct FunctionalList *f;

    dsigma.rhoa  = dsigma.rhob  = dp->rhoa;
    dsigma.grada = dsigma.gradb = dp->grada;
    dsigma.gradab = dsigma.grada*dsigma.gradb;

    memset(&f3, 0, sizeof(f3));
    for(ex = 0.0, f=exchangeFunctionals; f; f = f->next) {
        ex += 0.5 * f->func->func(&dsigma) * f->weight;
        f->func->third(&f3, f->weight, dp);
    }

    if(ex>0) {
        assert(ex < ENERGY_THR);
        return;
    }

    dfun.df10 = f3.df1000; dfun.df20 = f3.df2000; dfun.df30 = f3.df3000; 
    dfun.df01 = f3.df0010; dfun.df02 = f3.df0020; dfun.df03 = f3.df0030;
    dfun.df11 = f3.df1010; dfun.df21 = f3.df2010; dfun.df12 = f3.df1020;
    cam_third_sigma(dp->rhoa, ex, &dfun, &res);

    ds->df1000 += factor*(weight_lr*dfun.df10 + res.df10);
    ds->df0010 += factor*(weight_lr*dfun.df01 + res.df01);
    ds->df2000 += factor*(weight_lr*dfun.df20 + res.df20);
    ds->df1010 += factor*(weight_lr*dfun.df11 + res.df11);
    ds->df0020 += factor*(weight_lr*dfun.df02 + res.df02);

    ds->df3000 += factor*(weight_lr*dfun.df30 + res.df30);
    ds->df2010 += factor*(weight_lr*dfun.df21 + res.df21);
    ds->df1020 += factor*(weight_lr*dfun.df12 + res.df12);
    ds->df0030 += factor*(weight_lr*dfun.df03 + res.df03);

    if(fabs(dp->rhoa-dp->rhob)>THR || fabs(dp->grada-dp->gradb)>THR) {
        dsigma.rhoa  = dsigma.rhob  = dp->rhob;
        dsigma.grada = dsigma.gradb = dp->gradb;
	dsigma.gradab = dsigma.grada*dsigma.gradb;


        for(ex = 0.0, f=exchangeFunctionals; f; f = f->next)
            ex += 0.5 * f->func->func(&dsigma) * f->weight;

        if(ex>0) {
            assert(ex < ENERGY_THR);
            return;
        }

        dfun.df10 = f3.df0100; dfun.df20 = f3.df0200; dfun.df30 = f3.df0300; 
        dfun.df01 = f3.df0001; dfun.df02 = f3.df0002; dfun.df03 = f3.df0003;
        dfun.df11 = f3.df0101; dfun.df21 = f3.df0201; dfun.df12 = f3.df0102;
        cam_third_sigma(dp->rhob, ex, &dfun, &res);
    }

    ds->df0100 += factor*(weight_lr*dfun.df10 + res.df10);
    ds->df0001 += factor*(weight_lr*dfun.df01 + res.df01);
    ds->df0200 += factor*(weight_lr*dfun.df20 + res.df20);
    ds->df0101 += factor*(weight_lr*dfun.df11 + res.df11);
    ds->df0002 += factor*(weight_lr*dfun.df02 + res.df02);

    ds->df0300 += factor*(weight_lr*dfun.df30 + res.df30);
    ds->df0201 += factor*(weight_lr*dfun.df21 + res.df21);
    ds->df0102 += factor*(weight_lr*dfun.df12 + res.df12);
    ds->df0003 += factor*(weight_lr*dfun.df03 + res.df03);

    for(f=correlationFunctionals; f; f = f->next)
        f->func->third(ds, f->weight * factor, dp);
}



int
fun_get_cam_param(real *alpha, real *beta, real *mu)
{
    if(selected_func == &CamFunctional ||
       selected_func == &Camb3lypFunctional ||
       selected_func == &HseFunctional) {
        *alpha = camAlpha;
        *beta  = camBeta;
        /* Yes, HSE functional is inconsistently defined! */
        *mu    = selected_func == &HseFunctional
            ? camMu/M_SQRT2 : camMu;
        return 1;
    } else {
        *alpha = 0;
        *beta  = 0;
        *mu    = 0;
        return 0;
    }
}
/* ===================================================================
   TEST CODE
   =================================================================== */
#ifdef TEST
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char *argv[])
{
    FunDensProp dp;
    RGFirstDrv fa;
    real a, bfactor, bfactor_first, bfactor_second;
    if(argc<4) { printf("rhoa grada mu\n"); return 1;}
    dp.rhoa  = dp.rhob  = atof(argv[1]);
    dp.grada = dp.gradb = atof(argv[2]);
    CamMuFactor = atof(argv[3]);
    CamAlpha = 0.19;
    CamBeta  = 0.46;
    printf("mu=%f: energy(%g,%g) = %g\n", CamMuFactor, dp.rhoa*2, dp.grada*2,
           cam_energy(&dp));
    fun_a_first(1,1,&fa);
    printf("fun a derivatives: %g %g\n", fa.df10, fa.df01);
    for(a=0.1; a<10; a += 0.1) {
        bfactor_first  = cam_b_second_medium(a*FAC);/* EVALUATOR(a, second); */
        real b1 = EVALUATOR((a+1e-6), first);
        real b2 = EVALUATOR((a-1e-6), first);
        printf("%20.15g %20.15g %20.15g %20.15g\n", a,
               bfactor, bfactor_first*FAC, (b1-b2)/2e-6);
    }
    return 1;
}

void dftsethf_(real* s) {}
void dftsetcam_(real* s, real *b) {}
void fun_printf(const char *fmt, ...){}
#endif



