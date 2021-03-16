/*-----------------------------------------------------------------------------------------------------
  ------                                                                                         ------
  ------     This C module is a 'behind-the-curtain' stack system that only works                ------
  ------     with **real** numbers which is sufficient in most cases because complex arrays      ------
  ------     are simply stored as real streams of twice the length.                              ------
  ------                                                                                         ------
  ------     Extending the module beyond one simple data type would take away much of the        ------
  ------     code coherence.                                                                     ------
  ------                                                                                         ------
  ------     There is one special case in this module, namely the routine 'qstack_directx_c'     ------
  ------     where two real number streams of length n are interpreted as two complex            ------
  ------     streams of length n/2. These two streams are internally multiplied in the form      ------
  ------     sum_i dconjg(a(i)) * b(i) and the resulting complex number is returned.             ------
  ------     This was done for minimizing transactions with large vectors.                       ------
  ------                                                                                         ------
  -----------------------------------------------------------------------------------------------------*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include<complex.h>

#define MAXNUMSTACK 100
/*
 * the following type definitions are introduced to provide
 * functionality with integer*4 and integer*8 compilations.
 *
 * The same took place in the fortran module interface
 * polprp_qstackm.F90
 *
 */

#if defined (INT_STAR8)
#define CINT48  int64_t
#else
#define CINT48  int
#endif

/*
 *     Local variables managing the heap space for the vector stack.
 *     Visible only in this source file.
 *     Values are kept during program execution.
 *     All data are on the heap.
 */

static CINT48 initialized = 0;
static CINT48 numofstacks = 0;
static CINT48 verbose = 0;
static double *** qsbase = NULL;    /* stack distributor */
static CINT48 ** qslen = NULL;  /* length of each vector on each stack */
static CINT48 * qsdep = NULL; /* depth of stack #i */

/*
 * local variable for string output
 */

const char *taginf = "*** Qstack info:";
const char *tagerr = "*** Qstack error:";
const char *tagnit = "*** Qstack runtime error: Please initialize before use.";


/*
 *               function qstack_init
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_init(CINT48 *nost, CINT48 *vbflag){

  CINT48 i;
  CINT48 nos;

  if(initialized){
    printf("%s Stack already initialized. No action taken.\n",taginf);
    return 0;
  }
 
  nos = *nost;

  if(nos<1 || nos > MAXNUMSTACK){
    printf("%s max. %d stacks allowed!\n",tagerr,MAXNUMSTACK);
    return -1;
  }

/* 
 * now generate requested number of stacks;
 * initialize depths, length array and main distributor
 * array
 */

  qsdep = (CINT48 *) malloc(nos*sizeof(CINT48));
  memset(qsdep, 0, nos*sizeof(CINT48));

  qslen = (CINT48 **)malloc(nos*sizeof(CINT48 *));
  if(qslen == NULL){
#if defined (INT_STAR8)
    printf("%s Could not generate %lld stacks.\n",tagerr,nos);
#else
    printf("%s Could not generate %d stacks.\n",tagerr,nos);
#endif
    return -1;
  }
  for(i=0;i<nos;i++){
    qslen[i] = NULL;
  }

  qsbase = (double ***)malloc(nos*sizeof(double **));
  if(qsbase==NULL){
#if defined (INT_STAR8)
    printf("%s Could not generate %lld stacks.\n",tagerr,nos);
#else
    printf("%s Could not generate %d stacks.\n",tagerr,nos);
#endif
    return -1;
  }
  for(i=0;i<nos;i++){
    qsbase[i] = NULL;
  }
  initialized = 1;
  numofstacks = (CINT48)nos;
  verbose = (*vbflag)?1:0;

#if defined (INT_STAR8)
  if(verbose) printf("%s Created %lld stacks.\n",taginf,nos);
#else
  if(verbose) printf("%s Created %d stacks.\n",taginf,nos);
#endif
  return nos;
}



/*
 *               function qstack_push
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_push(CINT48 *stack, CINT48 *len, double *arr){ 

  CINT48 i,stackn;
  CINT48 curdep;
  CINT48 * newlen;
  double * curvec;
  double ** newvp;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }
  stackn = *stack - 1;  // C indexing starts at zero!

  if(*len < 1){
    printf("*** Qstack runtime problem: Illegal data length!\n");
    printf("*** No action taken.\n");
    return -1;
  }

/* first store vector item on the heap */

  curvec = (double *)malloc(*len * sizeof(double));
  memcpy(curvec,arr,*len * sizeof(double));

/* insert this address now to the corresponding stack and
 * adjust stack depths and update corresponding length
 * array.
 * create intermediate structure for copying existing pointers
 * and link it into main array at the end.
 */

  curdep = qsdep[stackn];
  curdep++;
  newvp = (double **)malloc(curdep*sizeof(double *));
  newlen = (CINT48 *)malloc(curdep*sizeof(CINT48));

  for(i=0;i<curdep-1;i++) {
    newvp[i]  = qsbase[stackn][i];   
    newlen[i] = qslen[stackn][i];
  }
  newvp[curdep-1] = curvec;
  newlen[curdep-1] = *len;
  free(qsbase[stackn]);   /* NULL pointer save! */
  free(qslen[stackn]);   /* NULL pointer save! */
  qsbase[stackn] = newvp;
  qslen[stackn]  = newlen;

  qsdep[stackn] = curdep;
#if defined (INT_STAR8)
  if(verbose) printf("%s PUSH: Stack #%lld has now %lld vectors stored.\n",taginf,*stack,curdep);
#else
  if(verbose) printf("%s PUSH: Stack #%d has now %d vectors stored.\n",taginf,*stack,curdep);
#endif

  return *len;
}



/*
 *               function qstack_popf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_popf(CINT48 *stack, double *retvec){
  CINT48 i;
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;
  double ** newvp;
  CINT48 * newlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* fetch >>first<< vector from stack (FIFO) and copy to retvec */

  curlen = qslen[stackn][0];
  memcpy(retvec,*qsbase[stackn],curlen*sizeof(double));

/* shrink stack from below and check special case for 1 vector
   on stack */

  curdep = qsdep[stackn];
  curdep--;
  free(qsbase[stackn][0]);  /* free lowest vector on selected stack */
  if(curdep == 0){    /*there was only one element, no loop required */
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn]  = NULL;
    qslen[stackn] = NULL;
  } else {
    newvp = (double **)malloc(curdep*sizeof(double *));
    newlen = (CINT48 *)malloc(curdep*sizeof(CINT48));
    for(i=0;i<curdep;i++) {
      newvp[i]  = qsbase[stackn][i+1];
      newlen[i] = qslen[stackn][i+1];
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn] = newvp;
    qslen[stackn] = newlen;
  }
  qsdep[stackn] = curdep;
#if defined (INT_STAR8)
  if(verbose) printf("%s POPF: Stack #%lld has now %lld vectors stored.\n",taginf,*stack,curdep);
#else
  if(verbose) printf("%s POPF: Stack #%d has now %d vectors stored.\n",taginf,*stack,curdep);
#endif
  return curlen;
}



/*
 *               function qstack_popl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_popl(CINT48 *stack, double *retvec){
  CINT48 i;
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;
  double ** newvp;
  CINT48 * newlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* fetch >>last<< vector from stack (LIFO) and copy to retvec */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep-1];
  memcpy(retvec,qsbase[stackn][curdep-1],curlen*sizeof(double));
/*
  for(i=0;i<curlen;i++) {
    retvec[i] = *(qsbase[stackn][curdep-1] + i);
  }
*/

/* shrink stack from above and check special case for 1 vector
   on stack */

  curdep--;
  free(qsbase[stackn][curdep]);  /* free topmost vector on selected stack */
  if(curdep == 0){    /*there was only one element, no loop required */
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn]  = NULL;
    qslen[stackn] = NULL;
  } else {
    newvp = (double **)malloc(curdep*sizeof(double *));
    newlen = (CINT48 *)malloc(curdep*sizeof(CINT48));
    for(i=0;i<curdep;i++) {
      newvp[i]  = qsbase[stackn][i];
      newlen[i] = qslen[stackn][i];
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn] = newvp;
    qslen[stackn] = newlen;
  }
  qsdep[stackn] = curdep;
#if defined (INT_STAR8)
  if(verbose) printf("%s POPL: Stack # %lld has now %lld vectors stored.\n",taginf,*stack,curdep);
#else
  if(verbose) printf("%s POPL: Stack # %d has now %d vectors stored.\n",taginf,*stack,curdep);
#endif
  return curlen;
}




/*
 *               function qstack_peekf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_peekf(CINT48 *stack, double *retvec){
  CINT48 stackn;
  CINT48 curlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and do not touch retvec array */

  if(qsdep[stackn] == 0) return 0;

/* copy >>first<< vector from stack (FIFO) to retvec, do not adjust stack */

  curlen = qslen[stackn][0];
  memcpy(retvec,qsbase[stackn][0],curlen*sizeof(double));

  return curlen;
}




/*
 *               function qstack_peekl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_peekl(CINT48 *stack, double *retvec){
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and do not touch retvec array */

  if(qsdep[stackn] == 0) return 0;

/* copy >>last<< vector from stack (LIFO) to retvec, do not adjust stack */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep-1];
  memcpy(retvec,qsbase[stackn][curdep-1],curlen*sizeof(double));

  return curlen;
}




/*
 *               function qstack_peekn
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_peekn(CINT48 *stack, double *retvec, CINT48 *shard){
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and do not touch retvec array */

  if(qsdep[stackn] == 0) return 0;

/* copy >>last<< vector from stack (LIFO) to retvec, do not adjust stack */

  curdep = qsdep[stackn];      /* number of shards on stack stackn */
  if(*shard < 1  ||  *shard > curdep){
    printf("*** Qstack runtime problem: Shard number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }
  curdep = *shard;
  curlen = qslen[stackn][curdep-1];
  memcpy(retvec,qsbase[stackn][curdep-1],curlen*sizeof(double));

  return curlen;
}






/*
 *               function qstack_compl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_compl(CINT48 *stack, double *retvec){
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return -1 and do not touch retvec array */

  if(qsdep[stackn] == 0) return -1;

/* compare >>last<< vector from stack (LIFO) to retvec, do not adjust stack */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep - 1];
  return memcmp(retvec,qsbase[stackn][curdep - 1],curlen*sizeof(double));
}






/*
 *               function qstack_getinf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_getinf(void){
   return (CINT48) numofstacks;
}





/*
 *               function qstack_meminf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
double qstack_meminf(void){

  double usedmem=0.0;
  uint64_t umem=0;
  CINT48 i,stackn;

  if(!initialized){
    printf("*** Qstack runtime problem: Stack not initialized!\n");
    printf("*** No action taken.\n");
    return 0.0;
  }

// contribution from qsbase

  for(stackn=0;stackn < numofstacks;stackn++){
    for(i=0;i<qsdep[stackn];i++){
      umem+=qslen[stackn][i]*sizeof(double);
    }
    umem+=qsdep[stackn]*sizeof(double *);
  }
  umem+=numofstacks*sizeof(double **);

// contribution from qslen

  for(stackn=0;stackn < numofstacks;stackn++){
    umem+=qsdep[stackn]*sizeof(CINT48);
  }
  umem+=numofstacks*sizeof(CINT48 *);

// contribution from qsdep

  umem+=numofstacks*sizeof(CINT48);

  usedmem=(double)(umem)/(1024.0 * 1024.0);  // mem in MB as a double
  if(verbose) printf("%s total allocated memory: %llu bytes = (%.2f Megabytes).\n",taginf,umem,usedmem);

  return usedmem;
}



/*
 *               function qstack_drop
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_drop(CINT48 * stack){

  CINT48 stackn,i,curdep;
  double ** curstack;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stackline, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* clear stack line # stackn */

  curstack = qsbase[stackn];
  for(i=0;i<qsdep[stackn];i++){
    free(curstack[i]);
  }
  free(qsbase[stackn]);
  free(qslen[stackn]);
  curdep = qsdep[stackn];
  qsdep[stackn] = 0;

  qsbase[stackn] = NULL;
  qslen[stackn] = NULL;

#if defined (INT_STAR8)
  if(verbose) printf("%s dropped %lld elements from stack %lld\n",taginf,curdep,*stack);
#else
  if(verbose) printf("%s dropped %d elements from stack %d\n",taginf,curdep,*stack);
#endif

  return curdep;
}






/*
 *               function qstack_directx_r
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_directx_r(CINT48 *stack, CINT48 *lenv, double *bravec, CINT48 *shard, double *result){
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;
  CINT48 i;
  double * ketadd;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and set result to zero */

  if(qsdep[stackn] == 0){
    printf("*** Qstack runtime problem: No ket vectors available!\n");
    *result = 0.0;
    return -1;
  }

/* do consistency checks */

  curdep = qsdep[stackn];      /* number of shards on stack stackn */
  if(*shard < 1  ||  *shard > curdep){
    printf("*** Qstack runtime problem: Shard number out of range!\n");
    *result = 0.0;
    return -1;
  }
  curdep = *shard;
  curlen = qslen[stackn][curdep-1];
  if(*lenv != curlen){
    printf("*** Qstack runtime problem: bra/ket lengths do not match!\n");
    *result = 0.0;
    return -1;
  }

/* do bra/ket multiplication */

  ketadd = qsbase[stackn][curdep-1];
  *result = 0.0;
  for(i=0;i<*lenv;i++){
    *result += *bravec++ * *ketadd++;
  }
  
  return *lenv;
}




/*
 *               function qstack_directx_c
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
CINT48 qstack_directx_c(CINT48 *stack, CINT48 *lenv, double *bravec, CINT48 *shard, double complex *result){
  CINT48 stackn;
  CINT48 curlen;
  CINT48 curdep;
  CINT48 i;
  double complex zsum;
  double complex * zbra, * zket;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and set result to zero */

  if(qsdep[stackn] == 0){
    printf("*** Qstack runtime problem: No ket vectors available!\n");
    *result = 0.0;
    return -1;
  }

/* do consistency checks */

  curdep = qsdep[stackn];      /* number of shards on stack stackn */
  if(*shard < 1  ||  *shard > curdep){
    printf("*** Qstack runtime problem: Shard number out of range!\n");
    *result = 0.0;
    return -1;
  }

/* The transferred number *lenv is the number of 
   --> complex numbers <-- on the stack stored as a real stream of twice the size.
   directx_c therefore only applies if the transferred length is exactly half of the
   storage length in the addressed shard  */

  curdep = *shard;
  curlen = qslen[stackn][curdep-1];
  if(  (curlen%2)  || (*lenv != curlen/2)){
    printf("*** Qstack runtime problem: vector lengths inappropriate for complex multiplication!\n");
    *result = 0.0;
    return -1;
  }

/* do complex bra/ket multiplication by pointer recast. Tested for
   GNU and INTEL compilers. Storage order is real/imaginary part in the data stream. 
   zket and zbra point to double complex numbers, offsetting by unity 
   therefore advances in steps of double complex numbers!  */

  zket = (double complex *) qsbase[stackn][curdep-1];
  zbra = (double complex *) bravec;

  zsum = 0.0;
  for(i=0;i<*lenv;i++){
    zsum += conj(*zbra++) * *zket++;
  }
  
  *result = zsum;

  return *lenv;
}






/*
 *               function qstack_shutdown
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
void qstack_shutdown(void){

  CINT48 i,stackn;
  double ** curstack;

  if(!initialized){
    printf("*** Qstack runtime problem: Stack not initialized!\n");
    printf("*** Qstack: No action taken.\n");
    return;
  }
 
  if(verbose) printf("%s Shutting down stack system. Clearing heap space.\n",taginf);

  for(stackn=0;stackn < numofstacks;stackn++){
    curstack = qsbase[stackn];       //clear now stack # stackn
    for(i=0;i<qsdep[stackn];i++){
      free(curstack[i]);
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
  }
  free(qsbase);
  free(qslen);
  free(qsdep);

  qsbase = NULL;
  qslen = NULL;
  qsdep = NULL;

  numofstacks = 0;
  initialized = 0;
  verbose = 0;

  return;
}
