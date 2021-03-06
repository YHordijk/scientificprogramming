/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.h: general definitions needed by C code:
   function prototypes, often used constants, etc.

   (c) Pawel Salek, pawsa@theochem.kth.se, feb 2002
*/
#ifndef _GENERAL_H_
#define _GENERAL_H_

#include <stdlib.h>

#if !defined(RESTRICT)
#define RESTRICT
#endif

/* define the basic floating-point variable type used by the Fortran code */

#if !defined(__CVERSION)
#define __CVERSION__
#endif

#include "functionals.h"

/* Match Fortran name mangling. If the Fortran compiler does not
 * mangle names, define NO_UNDERSCORE in CFLAGS.  g77 and compaq fort
 * (cryptically referred to with HAVE_GCPP below) for linux-alpha both
 * insert a second underscore if routine name contains at least one
 * underscore /hjaaj Oct04 */
#ifdef NO_UNDERSCORE
#define FSYM(a) a
#define FSYM2(a) a
#else
#define FSYM(a) a ## _
#if defined(VAR_G77) || defined(HAVE_GCPP)
#define FSYM2(a) a ## __
#else
#define FSYM2(a) a ## _
#endif
#endif

#if defined(VAR_PGF77)
#define __FUNCTION__ "PGI_does_not_define__FUNCTION__"
#endif
#if defined(SYS_SUN)
#define __FUNCTION__ "SUNs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_IRIX)
#define __FUNCTION__ "SGIs CC compiler_does_not_define__FUNCTION__"
#endif
#if defined(SYS_DEC)
#define __FUNCTION__ "DEC CC compiler does not define __FUNCTION__"
#endif

#define ELEMENTS(arr) (sizeof(arr)/sizeof(arr[0]))

typedef struct FirstDrv_  FirstDrv;
typedef struct SecondDrv_ SecondDrv;
typedef struct ThirdDrv_  ThirdDrv;

/* Density evaluators */
typedef struct DftDensity_  DftDensity;
typedef struct DftGrid_     DftGrid;

typedef void (*DftDensEvaluator)(DftDensity* dens, FunDensProp* dp,
                                 DftGrid* grid, real* tmp_vec);
struct DftDensity_ {
    DftDensEvaluator evaluate;
    real *dmata, *dmatb;
};


/*
 *     R = \rho_\alpha + \rho_\beta
 *     S = \rho_\alpha - \rho_\beta
 *    
 *     Z = \nabla R \cdot \nabla R
 *     Y = \nabla R \cdot \nabla S
 *     X = \nabla S \cdot \nabla S
 */

struct FirstDrv_ {
       real fR; 
       real fZ; 
};

struct SecondDrv_ {
       real fR; 
       real fZ; 
       real fX;
 
       real fRR;
       real fSS;
       real fZZ;
       real fYY;
       real fRZ;
       real fSY;
};

struct ThirdDrv_ {
       real fR; 
       real fZ; 
       real fX;
 
       real fRR;
       real fSS;
       real fZZ;
       real fYY;
       real fRZ;
       real fSY;
       real fRX;
       real fZX;
 
       real fRRR;
       real fZZZ;
       real fRSS;
       real fRZZ;
       real fRYY;
       real fRRZ;
       real fSSZ;
       real fRSY;
       real fSZY;
       real fZYY;
};

void dftpot0_(FirstDrv *ds, 
              const real* weight,
              const real* rho, 
              const real* grad);

void dftpot1_(SecondDrv *ds, 
              const real* weight,
              const real* rho, 
              const real* grad);

void dftpot2_(ThirdDrv *ds, 
              const real* weight,
              const real* rho, 
              const real* grad);

void sdftpot1_(SecondDrv *ds, 
               const real* weight,
               const real* rho, 
               const real* grad);

void sdftpot2_(ThirdDrv *ds, 
               const real* weight,
               const real* rho, 
               const real* grad);

/* Property evaluators */
void dft_kohn_sham_(real* dmat, real* ksm, real *edfty,
                    real* work, int *lwork, int* ipr);
void dft_lin_resp_(real* fmat, real *cmo, real *zymat, int *trplet, 
		   int *ksymop, real* work,int* lwork);
void FSYM2(dft_lin_respf)(int *nosim, real* fmat, real *cmo, real *zymat,
                          int *trplet, int *ksymop, real* work,int* lwork);
void dft_mol_grad_(real* dmat, real* work, int* lwork, int* iprint);
void dftqrcf_(real* fi, real* cmo, real* kappaY, int* symY, int* spinY,
              real* kappaZ, int* symZ, int* spinZ, int* addfock,
              real* work, int* lwork);


void dft_kohn_shamab_(real* dmat, real* ksm, real *edfty, 
                      real* work, int *lwork, int* ipr);
void dft_lin_respab_(real* fmatc, real* fmato,  real *cmo, real *zymat, 
                     int *trplet, int *ksymop, real* work,int* lwork);
void dftmolgradab_(real* work, int* lwork, int* iprint);
typedef void (*DFTPropEvalMaster)(void);
typedef void (*DFTPropEvalSlave)(real* work, int* lwork, const int* iprint);
#if defined(VAR_MPI)
#include <mpi.h>
void dft_kohn_sham_slave(real* work, int* lwork, const int* iprint);
void dft_lin_resp_slave (real* work, int* lwork, const int* iprint);
void dft_lin_respf_slave (real* work, int* lwork, const int* iprint);
void dft_kohn_shamab_slave(real* work, int* lwork, const int* iprint);
void dft_lin_respab_slave (real* work, int* lwork, const int* iprint);
void dft_mol_grad_slave (real* work, int* lwork, const int* iprint);
void dft_qr_resp_slave  (real* work, int* lwork, const int* iprint);
void dft_wake_slaves(DFTPropEvalMaster);
typedef struct {
    void*        data;
    int          count;
    MPI_Datatype type;
} SyncData;
void mpi_sync_data(const SyncData* data, int count);
#else
#define dft_wake_slaves(a)
#endif

void* dal_malloc_(size_t sz, const char *func, int line);
#define dal_malloc(sz) dal_malloc_((sz),__FUNCTION__, __LINE__)

int fort_print(const char* format, ...);
/* FORTRAN FUNCTION PROTOTYPES */
void dzero_(real* arr, const int* len);
void dunit_(real* arr, const int* len);
void outmat_(const real* mat, const int* rowlow, const int* rowhi,
	     const int* collow, const int* colhi,
             const int* rowdim, const int* coldim);
void getrho_(const real*dmat, const real* atv, real* rho, real* dmagao, 
	     const real* densthr);
void dftgrd_(real* work, int* lwork, const int* d1, const int* log1);
void dftdns_(real* dmat, real* work,int *lwork,int* iprint);
void gtdmso_(real* udv, real* cmo, real* di, real* dv, real* work);
void dftdnsab_(real* dmata,real* dmatb, real* work, int* lwork, int* iprint);
void udftmolgrdab_(real* gao, real* damta, real* dmatb, real* rha, real* rhb, 
                   real* vra, real* vrb, real* vza, real* vzb, real* vzg); 
int FSYM2(ishell_cnt)(void);
void dalton_quit(const char* format, ...);

/* BLAS and other linear algebra routines */
real dsum_(const int* cnt, const real* v, const int* stride);
void dscal_(const int* cnt, const real* fac, real* v, const int* stride);
real dnorm2_(const int* cnt, const real* v, const int* stride);
real ddot_(const int* cnt, const real* v1, const int* stride1,
	   const real* v2, const int* stride2);
real dcopy_(const int* cnt, const real* v1, const int* stride1,
	   const real* v2, const int* stride2);
void daxpy_(const int* cnt, const real* alpha, const real* v1, 
	    const int* stride1, const real* v2, const int* stride2);

void dger_(const int* m, const int* n, const real* alpha, 
	   const real* x, const int* incx, const real* y, const int* incy,
	   real* a, const int* lda);

void dgemv_(const char* tr, const int* nrows, const int* ncols, 
	    const real* alpha, const real* a, const int* lda, const real* vb,
	    const int* strideb, const real* beta, 
	    real* c, const int* stridec);

void dgemm_(const char* transa, const char* transb, const int* nrowc, 
	    const int* ncolc, const int* ncolopa, const real* alpha,
	    const real* a, const int* lda, const real* b, const int* ldb,
	    const real* beta, real* c, const int* ldc);

real* alloc_mat_MO(int cnt);

void FSYM2(dft_get_ao_dens_mat)(const real* cmo, real* dmat,
                                real* work, int* lwork);
void FSYM2(dft_get_ao_dens_matab)(real* cmo, real* dmata, real* dmatb,
                                  real* work, int* lwork);

/* useful  constants for fortran interfacing */
extern const int  ZEROI, ONEI, THREEI, FOURI;
extern const real ZEROR, ONER, TWOR, FOURR;

#if !defined __inline__
/* inline some stuff whenever possible */
#define __inline__
#endif

#endif
